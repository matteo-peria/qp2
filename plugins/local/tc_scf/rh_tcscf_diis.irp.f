! ---

subroutine rh_tcscf_diis()

  use diis_history_module

  implicit none

  type(diis_history)          :: history
  integer                        :: info, inserted_slot, it
  integer                        :: history_count
  integer, parameter             :: max_level_shift_tries = 5
  integer                        :: n_retry
  logical                        :: converged
  double precision               :: etc_tot, etc_1e, etc_2e, etc_3e
  double precision               :: e_save, e_delta
  double precision               :: level_shift_save, rate_th
  double precision               :: t0, t1
  double precision               :: er_DIIS, er_delta, er_save
  double precision               :: trial_residual
  double precision               :: rcond
  double precision, allocatable  :: coefficients(:)
  double precision, allocatable  :: fock_before_diis(:,:)
  double precision, allocatable  :: fock_extrapolated(:,:)
  double precision, allocatable  :: mo_r_coef_save(:,:)
  double precision, allocatable  :: mo_l_coef_save(:,:)

  double precision, allocatable  :: x_left(:,:), x_right(:,:)
  double precision               :: ao_residual
  logical, external              :: qp_stop

  PROVIDE level_shift_TCSCF
  PROVIDE mo_l_coef mo_r_coef

  it      = 0
  e_save  = 0.d0
  rate_th = 0.1d0

  allocate(mo_r_coef_save(ao_num,mo_num))
  allocate(mo_l_coef_save(ao_num,mo_num))
  allocate(fock_before_diis(ao_num,ao_num))
  allocate(fock_extrapolated(ao_num,ao_num))
  mo_l_coef_save = 0.d0
  mo_r_coef_save = 0.d0

  allocate(x_left(ao_num,mo_num), x_right(ao_num,mo_num))
  if (mo_num == ao_ortho_canonical_num) then
    x_left = ao_ortho_canonical_coef(:,1:mo_num)
    x_right = x_left
  else
    ! A truncated MO calculation tests only the supplied left/right spaces.
    call tc_orthogonalize_span(mo_l_coef, x_left)
    call tc_orthogonalize_span(mo_r_coef, x_right)
  endif

  call diis_initialize(history, ao_num, max_dim_DIIS_TCSCF, info)
  if(info /= diis_success) then
        stop 'TC-SCF DIIS history initialization failed'
  endif

  call write_time(6)

  write(6, '(A4,1X, A16,1X, A16,1X, A16,1X, A16,1X, A16,1X, A16,1X, A16,1X, A4, 1X, A8)')                      &
    '====', '================', '================', '================', '================', '================' &
          , '================', '================', '====', '========'
  write(6, '(A4,1X, A16,1X, A16,1X, A16,1X, A16,1X, A16,1X, A16,1X, A16,1X, A4, 1X, A8)')                      &
    ' it ', '  SCF TC Energy ', '      E(1e)     ', '      E(2e)     ', '      E(3e)     ', '   energy diff  ' &
          , '    DIIS error  ', '  level shift   ', 'DIIS', '  WT (m)'
  write(6, '(A4,1X, A16,1X, A16,1X, A16,1X, A16,1X, A16,1X, A16,1X, A16,1X, A4, 1X, A8)')                      &
    '====', '================', '================', '================', '================', '================' &
          , '================', '================', '====', '========'

  call wall_time(t0)

  etc_tot = TC_HF_energy
  etc_1e  = TC_HF_one_e_energy
  etc_2e  = TC_HF_two_e_energy
  etc_3e  = TC_HF_three_e_energy
  er_DIIS = maxval(abs(FQS_SQF_mo))
  e_delta = dabs(etc_tot - e_save)

  e_save  = etc_tot
  er_save = er_DIIS

  history_count = diis_count(history)
  call wall_time(t1)
  write(6, '(I4,1X, F16.10,1X, F16.10,1X, F16.10,1X, F16.10,1X, F16.10,1X, F16.10,1X, F16.10,1X, I4,1X, F8.2)')  &
    it, etc_tot, etc_1e, etc_2e, etc_3e, e_delta, er_DIIS, level_shift_tcscf, history_count, (t1-t0)/60.d0

  call tc_ao_residual_norm(x_left, x_right, ao_residual)
  write(6,'(A,I5,1X,ES22.14)') ' TC-SCF AO residual ', it, ao_residual

  PROVIDE FQS_SQF_ao Fock_matrix_tc_ao_tot

  converged = .false.
  call ezfio_set_tc_scf_converged_tcscf(converged)

  do while(.not. converged)
    call wall_time(t0)

    it += 1
    if(it > n_it_TCSCF_max) then
      print *, ' max of TCSCF iterations is reached ', n_it_TCSCF_max
      exit
    endif

    ! This snapshot is independent of whether the current pair is inserted.
    fock_before_diis(:,:) = Fock_matrix_tc_ao_tot(:,:)
    ! A level-shift retry must restart from the orbitals entering this step.
    mo_l_coef_save(1:ao_num,1:mo_num) = mo_l_coef(1:ao_num,1:mo_num)
    mo_r_coef_save(1:ao_num,1:mo_num) = mo_r_coef(1:ao_num,1:mo_num)
    level_shift_save = level_shift_TCSCF
    inserted_slot = 0
    rcond = 0.d0

    if(dabs(e_delta) > 1.d-12) then
      call diis_push( history, Fock_matrix_tc_ao_tot, FQS_SQF_ao, &
                         inserted_slot, info )
      if(info /= diis_success) then
        stop 'TC-SCF DIIS history insertion failed'
      endif

      call diis_extrapolate( history, fock_extrapolated, &
                                coefficients, rcond, info )
      if(info == diis_success) then
        Fock_matrix_tc_ao_tot(:,:) = fock_extrapolated(:,:)
        call ao_to_mo_bi_ortho( Fock_matrix_tc_ao_tot, &
                                size(Fock_matrix_tc_ao_tot, 1), &
                                Fock_matrix_tc_mo_tot, &
                                size(Fock_matrix_tc_mo_tot, 1) )
        TOUCH Fock_matrix_tc_mo_tot fock_matrix_tc_diag_mo_tot
      else
        ! A rejected solve keeps the unextrapolated Fock and starts fresh.
        call diis_reset(history)
        print *, ' TC-DIIS history reset after rejected solve, info = ', info
        print *, ' TC-DIIS reciprocal condition estimate = ', rcond
      endif
    endif

    mo_l_coef(1:ao_num,1:mo_num) = &
      fock_tc_leigvec_ao_balanced(1:ao_num,1:mo_num)
    mo_r_coef(1:ao_num,1:mo_num) = &
      fock_tc_reigvec_ao_balanced(1:ao_num,1:mo_num)
    TOUCH mo_l_coef mo_r_coef

    er_delta = maxval(abs(FQS_SQF_mo)) - er_save

    if((er_delta > rate_th * er_save) .and. (it > 1)) then
      Fock_matrix_tc_ao_tot(1:ao_num,1:ao_num) = &
        fock_before_diis(1:ao_num,1:ao_num)
      call ao_to_mo_bi_ortho( Fock_matrix_tc_ao_tot, &
                              size(Fock_matrix_tc_ao_tot, 1), &
                              Fock_matrix_tc_mo_tot, &
                              size(Fock_matrix_tc_mo_tot, 1) )
      TOUCH Fock_matrix_tc_mo_tot fock_matrix_tc_diag_mo_tot

      mo_l_coef(1:ao_num,1:mo_num) = &
        fock_tc_leigvec_ao_balanced(1:ao_num,1:mo_num)
      mo_r_coef(1:ao_num,1:mo_num) = &
        fock_tc_reigvec_ao_balanced(1:ao_num,1:mo_num)
      TOUCH mo_l_coef mo_r_coef
    endif

    er_delta = maxval(abs(FQS_SQF_mo)) - er_save
    n_retry = 0
    if((er_delta > rate_th * er_save) .and. (it > 1)) then
      print *, ' big or bad step '
      print *, er_delta, rate_th * er_save

      ! Discard entries belonging to the step that just failed, once.
      call diis_reset(history)

      ! Retry from the incoming orbitals. Keep a trial only if it improves
      ! their residual; the unextrapolated fallback already failed this test.
      do while(n_retry < max_level_shift_tries)
        n_retry += 1

        mo_l_coef(1:ao_num,1:mo_num) = &
          mo_l_coef_save(1:ao_num,1:mo_num)
        mo_r_coef(1:ao_num,1:mo_num) = &
          mo_r_coef_save(1:ao_num,1:mo_num)
        if(level_shift_TCSCF <= .1d0) then
          level_shift_TCSCF = 1.d0
        else
          level_shift_TCSCF = level_shift_TCSCF * 3.d0
        endif
        TOUCH mo_l_coef mo_r_coef level_shift_TCSCF

        mo_l_coef(1:ao_num,1:mo_num) = &
          fock_tc_leigvec_ao_balanced(1:ao_num,1:mo_num)
        mo_r_coef(1:ao_num,1:mo_num) = &
          fock_tc_reigvec_ao_balanced(1:ao_num,1:mo_num)
        TOUCH mo_l_coef mo_r_coef

        trial_residual = maxval(abs(FQS_SQF_mo))
        if(trial_residual < er_save) exit
      enddo

      if(.not. (trial_residual < er_save)) then
        mo_l_coef(:,:) = mo_l_coef_save(:,:)
        mo_r_coef(:,:) = mo_r_coef_save(:,:)
        TOUCH mo_l_coef mo_r_coef
        print *, ' TC-SCF retries did not improve the residual; restored incoming orbitals'
      endif
    endif

    ! Retries use temporary shifts. Keep the existing per-iteration shift
    ! decay, but do not let a trial shift leak into the next SCF step.
    if(n_retry > 0) then
      level_shift_TCSCF = level_shift_save
      SOFT_TOUCH level_shift_TCSCF
    endif
    level_shift_TCSCF = level_shift_TCSCF * 0.5d0
    SOFT_TOUCH level_shift_TCSCF

    etc_tot = TC_HF_energy
    etc_1e  = TC_HF_one_e_energy
    etc_2e  = TC_HF_two_e_energy
    etc_3e  = TC_HF_three_e_energy
    er_DIIS = maxval(abs(FQS_SQF_mo))
    e_delta = dabs(etc_tot - e_save)
    er_delta = er_DIIS - er_save

    e_save           = etc_tot
    er_save          = er_DIIS

    call tc_ao_residual_norm(x_left, x_right, ao_residual)
    write(6,'(A,I5,1X,ES22.14)') ' TC-SCF AO residual ', it, ao_residual
    converged = ao_residual .lt. dsqrt(thresh_tcscf)
    history_count = diis_count(history)

    call wall_time(t1)
    write(6, '(I4,1X, F16.10,1X, F16.10,1X, F16.10,1X, F16.10,1X, F16.10,1X, F16.10,1X, F16.10,1X, I4,1X, F8.2)')  &
      it, etc_tot, etc_1e, etc_2e, etc_3e, e_delta, er_DIIS, level_shift_tcscf, history_count, (t1-t0)/60.d0

    call lock_io
    if(it == 1) then
      write(json_unit, json_dict_uopen_fmt)
    else
      write(json_unit, json_dict_close_uopen_fmt)
    endif
    write(json_unit, json_int_fmt)  ' iteration      ', it
    write(json_unit, json_real_fmt) ' SCF TC Energy  ', etc_tot
    write(json_unit, json_real_fmt) ' E(1e)          ', etc_1e
    write(json_unit, json_real_fmt) ' E(2e)          ', etc_2e
    write(json_unit, json_real_fmt) ' E(3e)          ', etc_3e
    write(json_unit, json_real_fmt) ' AO residual    ', ao_residual
    write(json_unit, json_real_fmt) ' delta Energy   ', e_delta
    write(json_unit, json_real_fmt) ' DIIS error     ', er_DIIS
    write(json_unit, json_real_fmt) ' level_shift    ', level_shift_tcscf
    write(json_unit, json_int_fmtx) ' DIIS           ', history_count
    write(json_unit, json_real_fmt) ' Wall time (min)', (t1-t0)/60.d0
    call unlock_io

    if(er_delta .lt. 0.d0) then
      call ezfio_set_tc_scf_tcscf_energy(etc_tot)
      call ezfio_set_bi_ortho_mos_mo_l_coef(mo_l_coef)
      call ezfio_set_bi_ortho_mos_mo_r_coef(mo_r_coef)
      write(json_unit, json_true_fmt) 'saved'
    else
      write(json_unit, json_false_fmt) 'saved'
    endif
    call lock_io

    if(converged) then
      write(json_unit, json_true_fmtx) 'converged'
      call ezfio_set_tc_scf_converged_tcscf(converged)
    else
      write(json_unit, json_false_fmtx) 'converged'
    endif
    call unlock_io
    if(qp_stop()) exit
  enddo

  write(json_unit, json_dict_close_fmtx)

  if(converged) then
    print *, ' TCSCF DIIS converged !'
  else
    print *, 'Warning, TCSCF DIIS did not converge ...'
  endif
  call write_time(6)

  deallocate(mo_r_coef_save, mo_l_coef_save)
  deallocate(fock_before_diis, fock_extrapolated)

  call ezfio_set_tc_scf_tcscf_energy(TC_HF_energy)
  call ezfio_set_bi_ortho_mos_mo_l_coef(mo_l_coef)
  call ezfio_set_bi_ortho_mos_mo_r_coef(mo_r_coef)

end

! ---

