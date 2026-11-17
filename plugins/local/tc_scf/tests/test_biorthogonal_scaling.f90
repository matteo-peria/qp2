program test_biorthogonal_scaling

  use biorthogonal_scaling_module

  implicit none

  call test_reciprocal_balance_and_invariants()
  call test_reject_overlap_growth()
  call test_reject_zero_norm_without_partial_output()
  call test_reject_nonsymmetric_metric()

  print *, 'PASS: biorthogonal scaling tests'

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

    call assert_true(abs(actual - expected) <= tolerance, trim(message))

  end subroutine assert_close


  subroutine test_reciprocal_balance_and_invariants()

    implicit none

    double precision :: left(2,2), right(2,2), metric(2,2)
    double precision :: left_balanced(2,2), right_balanced(2,2)
    double precision :: left_scale(2), overlap_before(2,2), overlap_after(2,2)
    double precision :: projector_before(2,2), projector_after(2,2)
    double precision :: left_norm, right_norm
    integer :: info, i

    metric = 0.d0
    metric(1,1) = 4.d0
    metric(2,2) = 1.d0
    left = 0.d0
    right = 0.d0
    left(1,1) = 1.d8
    left(2,2) = 1.d-8
    right(1,1) = 2.5d-9
    right(2,2) = 1.d8

    overlap_before = matmul(transpose(left), matmul(metric, right))
    projector_before = matmul(right, transpose(left))

    call balance_biorthogonal_pairs(left, right, metric, &
                                    left_balanced, right_balanced, &
                                    left_scale, info)
    call assert_true(info == scaling_success, 'valid reciprocal scaling rejected')

    overlap_after = matmul(transpose(left_balanced), &
                           matmul(metric, right_balanced))
    projector_after = matmul(right_balanced, transpose(left_balanced))
    call assert_true(maxval(abs(overlap_after-overlap_before)) < 1.d-13, &
                     'cross overlap changed')
    call assert_true(maxval(abs(projector_after-projector_before)) < 1.d-13, &
                     'projector changed')

    do i = 1, 2
      left_norm = sqrt(dot_product(left_balanced(:,i), &
                                   matmul(metric,left_balanced(:,i))))
      right_norm = sqrt(dot_product(right_balanced(:,i), &
                                    matmul(metric,right_balanced(:,i))))
      call assert_close(left_norm, right_norm, 1.d-13, &
                        'left/right metric norms not balanced')
      call assert_true(left_scale(i) > 0.d0, 'invalid scale factor')
    enddo

  end subroutine test_reciprocal_balance_and_invariants


  subroutine test_reject_overlap_growth()

    implicit none

    double precision :: left(2,2), right(2,2), metric(2,2)
    double precision :: left_balanced(2,2), right_balanced(2,2)
    double precision :: left_scale(2), overlap_before(2,2), overlap_after(2,2)
    integer :: info

    metric = 0.d0
    metric(1,1) = 1.d0
    metric(2,2) = 1.d0
    left = 0.d0
    left(1,1) = 100.d0
    left(2,2) = 1.d-2
    right = 0.d0
    right(1,1) = 1.d-2
    right(2,1) = 1.d-4
    right(1,2) = 1.d-6
    right(2,2) = 100.d0
    overlap_before = matmul(transpose(left), matmul(metric, right))

    call balance_biorthogonal_pairs(left, right, metric, &
                                    left_balanced, right_balanced, &
                                    left_scale, info)
    call assert_true(info == scaling_overlap_changed, &
                     'balancing that grows off-diagonal overlap was accepted')
    call assert_true(maxval(abs(left_balanced-left)) <= tiny(1.d0), &
                     'rejected left matrix was partially changed')
    call assert_true(maxval(abs(right_balanced-right)) <= tiny(1.d0), &
                     'rejected right matrix was partially changed')
    call assert_true(maxval(abs(left_scale-1.d0)) <= epsilon(1.d0), &
                     'rejected scale factors were not reset')

    overlap_after = matmul(transpose(left_balanced), &
                           matmul(metric, right_balanced))
    call assert_true(maxval(abs(overlap_after-overlap_before)) <= tiny(1.d0), &
                     'rejected overlap changed')

  end subroutine test_reject_overlap_growth


  subroutine test_reject_zero_norm_without_partial_output()

    implicit none

    double precision :: left(2,2), right(2,2), metric(2,2)
    double precision :: left_balanced(2,2), right_balanced(2,2)
    double precision :: left_scale(2)
    integer :: info

    metric = 0.d0
    metric(1,1) = 1.d0
    metric(2,2) = 1.d0
    left = 0.d0
    left(1,1) = 1.d0
    right = 0.d0
    right(1,1) = 1.d0
    right(2,2) = 1.d0

    call balance_biorthogonal_pairs(left, right, metric, &
                                    left_balanced, right_balanced, &
                                    left_scale, info)
    call assert_true(info == scaling_bad_overlap, &
                     'zero pair overlap was not rejected')
    call assert_true(maxval(abs(left_balanced-left)) <= tiny(1.d0), &
                     'invalid input changed left output')
    call assert_true(maxval(abs(right_balanced-right)) <= tiny(1.d0), &
                     'invalid input changed right output')
    call assert_true(maxval(abs(left_scale-1.d0)) <= epsilon(1.d0), &
                     'invalid input changed scales')

  end subroutine test_reject_zero_norm_without_partial_output


  subroutine test_reject_nonsymmetric_metric()

    implicit none

    double precision :: left(2,1), right(2,1), metric(2,2)
    double precision :: left_balanced(2,1), right_balanced(2,1)
    double precision :: left_scale(1)
    integer :: info

    left = 0.d0
    right = 0.d0
    left(1,1) = 1.d0
    right(1,1) = 1.d0
    metric = 0.d0
    metric(1,1) = 1.d0
    metric(2,2) = 1.d0
    metric(1,2) = 0.1d0

    call balance_biorthogonal_pairs(left, right, metric, &
                                    left_balanced, right_balanced, &
                                    left_scale, info)
    call assert_true(info == scaling_bad_metric, &
                     'nonsymmetric metric was accepted')
    call assert_true(maxval(abs(left_balanced-left)) <= tiny(1.d0), &
                     'bad metric changed left output')
    call assert_true(maxval(abs(right_balanced-right)) <= tiny(1.d0), &
                     'bad metric changed right output')

  end subroutine test_reject_nonsymmetric_metric

end program test_biorthogonal_scaling
