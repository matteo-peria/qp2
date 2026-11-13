program test_diis_history

  use diis_history_module

  implicit none

  call test_wraparound_and_pairing()
  call test_reset_and_skipped_insertion()
  call test_gram_and_extrapolation()
  call test_ill_conditioned_rejection()

  print *, 'PASS: diis_history regression tests'

contains

  subroutine assert_true(condition, message)

    implicit none

    logical, intent(in) :: condition
    character(*), intent(in) :: message

    if(.not. condition) then
      print *, 'FAIL: ', trim(message)
      error stop 1
    endif

  end subroutine assert_true


  subroutine assert_close(actual, expected, tolerance, message)

    implicit none

    double precision, intent(in) :: actual, expected, tolerance
    character(*), intent(in) :: message

    call assert_true( abs(actual - expected) <= tolerance, &
                      trim(message) )

  end subroutine assert_close


  subroutine push_label(history, label)

    implicit none

    type(diis_history), intent(inout) :: history
    integer, intent(in) :: label
    double precision :: fock(2,2), error(2,2)
    integer :: info, inserted_slot

    fock = dble(label)
    error = -dble(label)
    call diis_push(history, fock, error, inserted_slot, info)
    call assert_true(info == diis_success, 'push_label failed')

  end subroutine push_label


  subroutine assert_history_labels(history, expected, message)

    implicit none

    type(diis_history), intent(in) :: history
    integer, intent(in) :: expected(:)
    character(*), intent(in) :: message
    double precision :: fock(2,2), error(2,2)
    integer :: i, info

    call assert_true( diis_count(history) == size(expected), &
                      trim(message)//': wrong count' )
    do i = 1, size(expected)
      call diis_get_entry(history, i, fock, error, info)
      call assert_true(info == diis_success, &
                       trim(message)//': get_entry failed')
      call assert_close(fock(1,1), dble(expected(i)), 0.d0, &
                        trim(message)//': wrong Fock chronology')
      call assert_close(error(1,1), -dble(expected(i)), 0.d0, &
                        trim(message)//': Fock/error pair separated')
    enddo

  end subroutine assert_history_labels


  subroutine test_wraparound_and_pairing()

    implicit none

    type(diis_history) :: history
    integer :: info, iteration

    call diis_initialize(history, 2, 4, info)
    call assert_true(info == diis_success, 'initialization failed')

    do iteration = 1, 4
      call push_label(history, iteration)
    enddo
    call assert_history_labels(history, [1, 2, 3, 4], &
                               'history before wraparound')
    call assert_true(diis_next_slot(history) == 1, &
                     'next slot before wraparound')

    call push_label(history, 5)
    call assert_history_labels(history, [2, 3, 4, 5], &
                               'iteration 5 must evict iteration 1')
    call assert_true(diis_next_slot(history) == 2, &
                     'next slot after first wraparound')

    do iteration = 6, 8
      call push_label(history, iteration)
    enddo
    call assert_history_labels(history, [5, 6, 7, 8], &
                               'history after full wraparound')

  end subroutine test_wraparound_and_pairing


  subroutine test_reset_and_skipped_insertion()

    implicit none

    type(diis_history) :: history
    integer :: info

    call diis_initialize(history, 2, 4, info)
    call push_label(history, 17)
    call push_label(history, 18)
    call diis_reset(history)
    call assert_true(diis_count(history) == 0, 'reset count')
    call assert_true(diis_next_slot(history) == 1, 'reset next slot')

    call push_label(history, 19)
    call push_label(history, 20)
    call assert_history_labels(history, [19, 20], &
                               'fresh chronology after reset')

    ! A skipped insertion is represented by not calling diis_push.
    call assert_history_labels(history, [19, 20], &
                               'skipped insertion must not change history')

  end subroutine test_reset_and_skipped_insertion


  subroutine test_gram_and_extrapolation()

    implicit none

    type(diis_history) :: history
    double precision :: fock(2,2), error(2,2), extrapolated(2,2)
    double precision, allocatable :: gram(:,:), coefficients(:)
    double precision :: rcond
    integer :: info, inserted_slot

    call diis_initialize(history, 2, 2, info)

    ! Fill slot 1 with a stale pair. The two useful entries below then occupy
    ! physical slots 2 and 1, so extrapolation must follow logical chronology.
    fock = -999.d0
    error = -999.d0
    call diis_push(history, fock, error, inserted_slot, info)

    fock = 0.d0
    fock(1,1) = 10.d0
    fock(2,2) = 10.d0
    error = 0.d0
    error(1,1) = 1.d0
    call diis_push(history, fock, error, inserted_slot, info)

    fock = 0.d0
    fock(1,1) = 20.d0
    fock(2,2) = 20.d0
    error = 0.d0
    error(2,2) = 2.d0
    call diis_push(history, fock, error, inserted_slot, info)

    call diis_build_gram(history, gram, info)
    call assert_true(info == diis_success, 'Gram construction failed')
    call assert_close(gram(1,1), 1.d0, 1.d-14, 'Gram(1,1)')
    call assert_close(gram(1,2), 0.d0, 1.d-14, 'Gram(1,2)')
    call assert_close(gram(2,1), 0.d0, 1.d-14, 'Gram(2,1)')
    call assert_close(gram(2,2), 4.d0, 1.d-14, 'Gram(2,2)')

    call diis_extrapolate( history, extrapolated, coefficients, &
                              rcond, info )
    call assert_true(info == diis_success, 'extrapolation failed')
    call assert_close(coefficients(1), 0.8d0, 1.d-13, 'coefficient 1')
    call assert_close(coefficients(2), 0.2d0, 1.d-13, 'coefficient 2')
    call assert_close(sum(coefficients), 1.d0, 1.d-13, &
                      'coefficient sum')
    call assert_close(extrapolated(1,1), 12.d0, 1.d-12, &
                      'extrapolated Fock(1,1)')
    call assert_close(extrapolated(2,2), 12.d0, 1.d-12, &
                      'extrapolated Fock(2,2)')

  end subroutine test_gram_and_extrapolation


  subroutine test_ill_conditioned_rejection()

    implicit none

    double precision :: gram(2,2), rcond
    double precision, allocatable :: coefficients(:)
    integer :: info

    gram = 1.d0
    call diis_solve_coefficients(gram, coefficients, rcond, info)
    call assert_true( info == diis_lapack_failure .or. &
                      info == diis_ill_conditioned, &
                      'singular DIIS system must be rejected' )

  end subroutine test_ill_conditioned_rejection

end program test_diis_history
