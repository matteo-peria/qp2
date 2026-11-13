! ---

program diis_history_qp_test

  use diis_history_module

  implicit none

  type(diis_history) :: history
  double precision      :: fock(2,2), error(2,2)
  integer               :: expected, info, inserted_slot, iteration

  call diis_initialize(history, 2, 3, info)
  if(info /= diis_success) then
    stop 'QP DIIS test: initialization failed'
  endif

  do iteration = 1, 5
    fock(:,:) = dble(iteration)
    error(:,:) = -dble(iteration)
    call diis_push(history, fock, error, inserted_slot, info)
    if(info /= diis_success) then
      stop 'QP DIIS test: insertion failed'
    endif
    expected = modulo(iteration - 1, 3) + 1
    if(inserted_slot /= expected) then
      stop 'QP DIIS test: wrong physical insertion slot'
    endif
  enddo

  if(diis_count(history) /= 3) then
    stop 'QP DIIS test: wrong count after wraparound'
  endif
  if(diis_next_slot(history) /= 3) then
    stop 'QP DIIS test: wrong next slot after wraparound'
  endif

  do iteration = 1, 3
    call diis_get_entry(history, iteration, fock, error, info)
    if(info /= diis_success) then
      stop 'QP DIIS test: logical retrieval failed'
    endif
    expected = iteration + 2
    if(fock(1,1) /= dble(expected)) then
      stop 'QP DIIS test: wrong logical Fock chronology'
    endif
    if(error(1,1) /= -dble(expected)) then
      stop 'QP DIIS test: Fock/error pair separated'
    endif
  enddo

  call diis_reset(history)
  if(diis_count(history) /= 0) then
    stop 'QP DIIS test: reset did not clear logical history'
  endif
  if(diis_next_slot(history) /= 1) then
    stop 'QP DIIS test: reset did not restore first slot'
  endif

  print *, 'PASS: QP/IRPF90 TC-SCF DIIS history integration test'

end

! ---
