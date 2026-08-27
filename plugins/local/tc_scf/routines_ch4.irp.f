

  BEGIN_PROVIDER [ double precision, integral_charge_ch4, (ao_num,ao_num)]
  implicit none
 double precision :: R(3,2),charge(2)
  integer          :: num_A, num_B, power_A(3), power_B(3)
  integer          :: i, j, k, l, n_pt_in, m
  double precision :: alpha, beta
  double precision :: A_center(3),B_center(3),C_center(3)
  double precision :: overlap_x,overlap_y,overlap_z,overlap,dx,NAI_pol_mult
 R(1:3,1) = nucl_coord_transp(1:3,2)
 R(1,1) += 0.28
 R(1:3,2) = nucl_coord_transp(1:3,4)
 R(3,2) -= 1.3
 charge(1) = 0.32
 charge(2) = -1.4
 n_pt_in = n_pt_max_integrals


 integral_charge_ch4 = 0.d0
 do j = 1, ao_num
   num_A = ao_nucl(j)
   power_A(1:3)= ao_power(j,1:3)
   A_center(1:3) = nucl_coord(num_A,1:3)

   do i = 1, ao_num

     num_B = ao_nucl(i)
     power_B(1:3)= ao_power(i,1:3)
     B_center(1:3) = nucl_coord(num_B,1:3)

     do l=1,ao_prim_num(j)
       alpha = ao_expo_ordered_transp(l,j)

       do m=1,ao_prim_num(i)
         beta = ao_expo_ordered_transp(m,i)

         double precision               :: c, c1
         c = 0.d0

         do  k = 1,2 
           double precision               :: Z
           Z = charge(k)

           C_center(1:3) = R(1:3,k)

           !print *, ' '
           !print *, A_center, B_center, C_center, power_A, power_B
           !print *, alpha, beta

           c1 = NAI_pol_mult( A_center, B_center, power_A, power_B &
                            , alpha, beta, C_center, n_pt_in )

           !print *, ' c1 = ', c1

           c = c - Z * c1

         enddo
         integral_charge_ch4(i,j) = integral_charge_ch4(i,j)  &
             + ao_coef_normalized_ordered_transp(l,j)             &
             * ao_coef_normalized_ordered_transp(m,i) * c
       enddo
     enddo
   enddo
 enddo
!do i = 1, ao_num
! write(*,'(100(F16.10,X))')integral_charge_ch4(i,:)
!enddo
END_PROVIDER

BEGIN_PROVIDER [ integer, n_orb_extract]
 implicit none
  n_orb_extract = 3
END_PROVIDER
BEGIN_PROVIDER [ integer, list_orb_extract, (n_orb_extract)]
 implicit none
  list_orb_extract(1) = 3
  list_orb_extract(2) = 4
  list_orb_extract(3) = 5
END_PROVIDER

subroutine ao_to_mo_general(A_ao,LDA_ao,A_mo,n_mo,mo_coef_local)
  implicit none
  BEGIN_DOC
  ! Transform A from the |AO| basis to the |MO| basis
  !
  ! $C^\dagger.A_{ao}.C$
  END_DOC
  integer, intent(in)            :: LDA_ao,n_mo
  double precision, intent(in)   :: A_ao(LDA_ao,ao_num),mo_coef_local(ao_num, n_mo) 
  double precision, intent(out)  :: A_mo(n_mo,n_mo)
  double precision, allocatable  :: T(:,:)
  integer :: LDA_mo
  LDA_mo = n_mo 
  allocate ( T(ao_num,n_mo) )
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: T

  call dgemm('N','N', ao_num, n_mo, ao_num,                    &
      1.d0, A_ao,LDA_ao,                                             &
      mo_coef_local, size(mo_coef_local,1),                                      &
      0.d0, T, size(T,1))

  call dgemm('T','N', n_mo, n_mo, ao_num,                &
      1.d0, mo_coef_local,size(mo_coef_local,1),                                 &
      T, ao_num,                                                     &
      0.d0, A_mo, size(A_mo,1))

  deallocate(T)
  end
