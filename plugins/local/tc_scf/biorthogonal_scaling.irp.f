BEGIN_PROVIDER [ double precision, fock_tc_reigvec_ao_balanced, (ao_num, mo_num) ]
&BEGIN_PROVIDER [ double precision, fock_tc_leigvec_ao_balanced, (ao_num, mo_num) ]

  use biorthogonal_scaling_module

  implicit none

  integer :: info
  double precision :: left_scale(mo_num)

  PROVIDE fock_tc_reigvec_ao fock_tc_leigvec_ao ao_overlap

  call balance_biorthogonal_pairs( fock_tc_leigvec_ao, &
                                   fock_tc_reigvec_ao, ao_overlap, &
                                   fock_tc_leigvec_ao_balanced, &
                                   fock_tc_reigvec_ao_balanced, &
                                   left_scale, info )

  if(info == scaling_success) then
    if(maxval(abs(log(left_scale))) > log(10.d0)) then
      print *, ' TC-SCF left/right orbital scaling range = ', &
               minval(left_scale), maxval(left_scale)
    endif
  else
    print *, ' TC-SCF left/right orbital balancing skipped, status = ', info
  endif

END_PROVIDER
