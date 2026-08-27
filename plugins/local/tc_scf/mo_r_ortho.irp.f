
BEGIN_PROVIDER [ double precision, overlap_orb_extract, (n_orb_extract,n_orb_extract)]
 implicit none
 integer :: i,ii,j,jj
 do i = 1, n_orb_extract
  ii = list_orb_extract(i)
  do j = 1, n_orb_extract
   jj = list_orb_extract(j)
   overlap_orb_extract(i,j) = overlap_mo_r(ii,jj)
  enddo
 enddo
!print*,'overlap_orb_extract'
!do i = 1, n_orb_extract
! write(*,'(100(F16.10,X))')overlap_orb_extract(i,:)
!enddo
END_PROVIDER 

BEGIN_PROVIDER [ double precision, eigenvec_overlap_orb_extract, (n_orb_extract,n_orb_extract)]
 implicit none
 integer :: i
 double precision, allocatable :: eigvalues(:)
 allocate(eigvalues(n_orb_extract))
 call lapack_diagd(eigvalues,eigenvec_overlap_orb_extract,overlap_orb_extract,n_orb_extract,n_orb_extract)
!print*,'eigval overlap_orb_extract'
 do i = 1, n_orb_extract
! print*,eigvalues(i)
  eigenvec_overlap_orb_extract(1:n_orb_extract,i) *= 1.d0/dsqrt(eigvalues(i))
 enddo
!print*,'eigvec overlap_orb_extract'
!do i = 1, n_orb_extract
! write(*,'(100(F16.10,X))')eigenvec_overlap_orb_extract(i,:)
!enddo
END_PROVIDER 

BEGIN_PROVIDER [ double precision, ao_eigenvec_overlap_orb_extract, (ao_num,n_orb_extract)] 
 implicit none
 integer :: i,j,k
 ao_eigenvec_overlap_orb_extract = 0.d0
 print*,'ao_eigenvec_overlap_orb_extract'
 do i = 1, n_orb_extract
  do j = 1, ao_num
   do k = 1, n_orb_extract
     ao_eigenvec_overlap_orb_extract(j,i) += eigenvec_overlap_orb_extract(k,i) * mo_r_coef(j,list_orb_extract(k))
   enddo
  enddo
 enddo
!do i = 1, n_orb_extract
!   write(*,'(100(F16.10,X))') ao_eigenvec_overlap_orb_extract(:,i)
!enddo
END_PROVIDER 

BEGIN_PROVIDER [ double precision, integral_charge_ch4_mo, (n_orb_extract,n_orb_extract)]
 implicit none
 call ao_to_mo_general(integral_charge_ch4,ao_num,integral_charge_ch4_mo,n_orb_extract,ao_eigenvec_overlap_orb_extract) 
END_PROVIDER 

BEGIN_PROVIDER [ double precision, overlap_eigenvec_overlap_orb_extract, (n_orb_extract,n_orb_extract)]
 implicit none
 call ao_to_mo_general(ao_overlap,ao_num,  & 
                      overlap_eigenvec_overlap_orb_extract,n_orb_extract,ao_eigenvec_overlap_orb_extract) 
 integer :: i
!print*,'overlap_eigenvec_overlap_orb_extract'
!do i = 1,n_orb_extract
! write(*,'(100(F16.10,X))')overlap_eigenvec_overlap_orb_extract(i,:)
!enddo
END_PROVIDER 


 BEGIN_PROVIDER [ double precision, eigevec_integral_charge_ch4_mo, (n_orb_extract,n_orb_extract)]
&BEGIN_PROVIDER [ double precision, eigeval_integral_charge_ch4_mo, (n_orb_extract)]
 implicit none
 integer :: i
 call lapack_diagd(eigeval_integral_charge_ch4_mo,eigevec_integral_charge_ch4_mo,integral_charge_ch4_mo,n_orb_extract,n_orb_extract)
 do i = 1, n_orb_extract
  print*,'eigeval_integral_charge_ch4_mo',eigeval_integral_charge_ch4_mo(i)
 enddo
 
END_PROVIDER

BEGIN_PROVIDER [ double precision, ao_coef_eigevec_integral, (ao_num,n_orb_extract)]
 implicit none
 integer :: i,j,k
 ao_coef_eigevec_integral = 0.d0
 do i = 1, n_orb_extract
  do j = 1, ao_num
   do k = 1, n_orb_extract
     ao_coef_eigevec_integral(j,i) += eigevec_integral_charge_ch4_mo(k,i) * ao_eigenvec_overlap_orb_extract(j,k)
!   print*,eigevec_integral_charge_ch4_mo(k,i) , ao_eigenvec_overlap_orb_extract(j,k)
   enddo
  enddo
 enddo
END_PROVIDER


BEGIN_PROVIDER [ double precision, overlap_ao_coef_eigevec_integral, (n_orb_extract,n_orb_extract)]
 implicit none
 integer :: i,j,k,l
 print*,'overlap between the eigenvectors of strange operator'
 call ao_to_mo_general(ao_overlap,ao_num,  & 
                      overlap_ao_coef_eigevec_integral,n_orb_extract,ao_coef_eigevec_integral) 
 print*,'should be a diagonal matrix'
 do i = 1, n_orb_extract
  write(*,'(100(F16.10,X))')overlap_ao_coef_eigevec_integral(i,:)
 enddo
 
END_PROVIDER

