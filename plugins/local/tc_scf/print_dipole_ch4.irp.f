program print_dipole_ch4
 implicit none
 integer :: i,j
 double precision :: norm
 do i =2, 5
  print*,'orbital ',i
  print*,mo_dipole_x(i,i),mo_dipole_y(i,i),mo_dipole_z(i,i)
  norm = mo_dipole_x(i,i)**2+ mo_dipole_y(i,i)**2 + mo_dipole_z(i,i)**2
  norm = dsqrt(norm)
  print*,'norm = ',norm
 enddo
  
end
