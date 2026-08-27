subroutine save_tc_mos_right_normalized
 implicit none
 integer :: i
BEGIN_DOC
 ! routine that normalizes the right mos to 1 
 ! 
 ! and correspondingly the left mos such that bi orthonormality is conserved
END_DOC
 do i = 1, mo_num
   mo_r_coef(1:ao_num,i) *= 1.d0/dsqrt(overlap_mo_r(i,i))
   mo_l_coef(1:ao_num,i) *= dsqrt(overlap_mo_r(i,i))
 enddo
 touch mo_r_coef mo_l_coef 

 call ezfio_set_bi_ortho_mos_mo_l_coef(mo_l_coef)
 call ezfio_set_bi_ortho_mos_mo_r_coef(mo_r_coef)

end
