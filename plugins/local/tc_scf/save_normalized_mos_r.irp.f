program save_normalized_mos
 implicit none
 integer :: i,j 
 double precision, allocatable :: mo_coef_new(:,:)
 allocate(mo_coef_new(ao_num, mo_num))
 do i = 1, mo_num
  do j = 1, ao_num
   mo_coef_new(j,i) = mo_r_coef(j,i) * 1.d0/dsqrt(overlap_mo_r(i,i))
  enddo
 enddo
 mo_coef = mo_coef_new
 SOFT_TOUCH mo_coef
 call save_mos

end
