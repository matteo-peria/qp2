
BEGIN_PROVIDER [ double precision, overlap_orb_read_extract, (n_orb_extract,n_orb_extract)]
 implicit none
 integer :: i,ii,j,jj
 do i = 1, n_orb_extract
  ii = list_orb_extract(i)
  do j = 1, n_orb_extract
   jj = list_orb_extract(j)
   overlap_orb_read_extract(i,j) = overlap_mo_r_read_read(ii,jj)
  enddo
 enddo
 print*,'overlap_orb_read_extract'
 do i = 1, n_orb_extract
  write(*,'(100(F16.10,X))')overlap_orb_read_extract(i,:)
 enddo
END_PROVIDER 

BEGIN_PROVIDER [ double precision, eigenvec_overlap_orb_read_extract, (n_orb_extract,n_orb_extract)]
 implicit none
 integer :: i
 double precision, allocatable :: eigvalues(:)
 allocate(eigvalues(n_orb_extract))
 call lapack_diagd(eigvalues,eigenvec_overlap_orb_read_extract,overlap_orb_read_extract,n_orb_extract,n_orb_extract)
 print*,'eigval overlap_orb_read_extract'
 do i = 1, n_orb_extract
  print*,eigvalues(i)
  eigenvec_overlap_orb_read_extract(1:n_orb_extract,i) *= 1.d0/dsqrt(eigvalues(i))
 enddo
 print*,'eigvec overlap_orb_read_extract'
 do i = 1, n_orb_extract
  write(*,'(100(F16.10,X))')eigenvec_overlap_orb_read_extract(i,:)
 enddo
END_PROVIDER 

BEGIN_PROVIDER [ double precision, ao_eigenvec_overlap_orb_read_extract, (ao_num,n_orb_extract)] 
 implicit none
 integer :: i,j,k
 ao_eigenvec_overlap_orb_read_extract = 0.d0
 do i = 1, n_orb_extract
  do j = 1, ao_num
   do k = 1, n_orb_extract
     ao_eigenvec_overlap_orb_read_extract(j,i) += eigenvec_overlap_orb_read_extract(k,i) * mo_r_coef_read(j,list_orb_extract(k))
   enddo
  enddo
 enddo
!do i = 1, n_orb_extract
!   write(*,'(100(F16.10,X))') ao_eigenvec_overlap_orb_read_extract(:,i)
!enddo
END_PROVIDER 

BEGIN_PROVIDER [ double precision, integral_charge_ch4_mo_read, (n_orb_extract,n_orb_extract)]
 implicit none
 call ao_to_mo_general(integral_charge_ch4,ao_num,integral_charge_ch4_mo_read,n_orb_extract,ao_eigenvec_overlap_orb_read_extract) 
END_PROVIDER 

BEGIN_PROVIDER [ double precision, overlap_eigenvec_overlap_orb_read_extract, (n_orb_extract,n_orb_extract)]
 implicit none
 call ao_to_mo_general(ao_overlap,ao_num,  & 
                      overlap_eigenvec_overlap_orb_read_extract,n_orb_extract,ao_eigenvec_overlap_orb_read_extract) 
 integer :: i
print*,'overlap_eigenvec_overlap_orb_read_extract'
do i = 1,n_orb_extract
 write(*,'(100(F16.10,X))')overlap_eigenvec_overlap_orb_read_extract(i,:)
enddo
END_PROVIDER 


 BEGIN_PROVIDER [ double precision, eigevec_integral_charge_ch4_mo_read, (n_orb_extract,n_orb_extract)]
&BEGIN_PROVIDER [ double precision, eigeval_integral_charge_ch4_mo_read, (n_orb_extract)]
 implicit none
 integer :: i
 call lapack_diagd(eigeval_integral_charge_ch4_mo_read,eigevec_integral_charge_ch4_mo_read,integral_charge_ch4_mo_read,n_orb_extract,n_orb_extract)
 do i = 1, n_orb_extract
  print*,'eigeval_integral_charge_ch4_mo_read',eigeval_integral_charge_ch4_mo_read(i)
 enddo
 
END_PROVIDER

BEGIN_PROVIDER [ double precision, ao_coef_eigevec_read_integral, (ao_num,n_orb_extract)]
 implicit none
 integer :: i,j,k
 ao_coef_eigevec_read_integral = 0.d0
 do i = 1, n_orb_extract
  do j = 1, ao_num
   do k = 1, n_orb_extract
     ao_coef_eigevec_read_integral(j,i) += eigevec_integral_charge_ch4_mo_read(k,i) * ao_eigenvec_overlap_orb_read_extract(j,k)
!   print*,eigevec_integral_charge_ch4_mo_read(k,i) , ao_eigenvec_overlap_orb_read_extract(j,k)
   enddo
  enddo
 enddo
END_PROVIDER


BEGIN_PROVIDER [ double precision, overlap_ao_coef_eigevec_read_integral, (n_orb_extract,n_orb_extract)]
 implicit none
 integer :: i,j,k,l
 print*,'overlap between the eigenvectors of strange operator'
 call ao_to_mo_general(ao_overlap,ao_num,  & 
                      overlap_ao_coef_eigevec_read_integral,n_orb_extract,ao_coef_eigevec_read_integral) 
 print*,'should be a diagonal matrix'
 do i = 1, n_orb_extract
  write(*,'(100(F16.5,X))')overlap_ao_coef_eigevec_read_integral(i,:)
 enddo
 
END_PROVIDER

BEGIN_PROVIDER [ double precision, overlap_charge_op_read, (n_orb_extract,n_orb_extract)]
 implicit none
 integer :: i,j,k,l
 print*,'overlap between the eigenvectors of strange operator in the read and no read basis'
 overlap_charge_op_read = 0.d0
 do i = 1,n_orb_extract 
  do j = 1,n_orb_extract 
   do k = 1, ao_num
    do l = 1, ao_num
      overlap_charge_op_read(j,i) += ao_coef_eigevec_read_integral(l,j) * ao_overlap(l,k) * ao_coef_eigevec_integral(k,i)
    enddo
   enddo
  enddo
 enddo
 do i = 1, n_orb_extract
  write(*,'(100(F16.10,x))')overlap_charge_op_read(i,:)
 enddo
END_PROVIDER

