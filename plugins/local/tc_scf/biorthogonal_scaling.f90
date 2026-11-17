module biorthogonal_scaling_module

  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

  implicit none
  private

  integer, parameter, public :: scaling_success = 0
  integer, parameter, public :: scaling_bad_shape = 1
  integer, parameter, public :: scaling_bad_metric = 2
  integer, parameter, public :: scaling_bad_norm = 3
  integer, parameter, public :: scaling_bad_overlap = 4
  integer, parameter, public :: scaling_nonfinite = 5
  integer, parameter, public :: scaling_overlap_changed = 6

  public :: balance_biorthogonal_pairs

contains

  pure subroutine balance_biorthogonal_pairs(left, right, metric, &
                                              left_balanced, right_balanced, &
                                              left_scale, info)

    implicit none

    double precision, intent(in) :: left(:,:), right(:,:), metric(:,:)
    double precision, intent(out) :: left_balanced(:,:), right_balanced(:,:)
    double precision, intent(out) :: left_scale(:)
    integer, intent(out) :: info

    double precision, allocatable :: overlap_before(:,:), overlap_after(:,:)
    double precision, allocatable :: left_trial(:,:), right_trial(:,:)
    double precision, allocatable :: scale_trial(:)
    double precision :: left_norm_sq, right_norm_sq
    double precision :: left_norm, right_norm, alpha
    double precision :: metric_scale, metric_asymmetry
    double precision :: overlap_scale, overlap_tolerance
    double precision :: offdiag_before, offdiag_after
    double precision :: diagonal_change
    integer :: i, j, n, m

    info = scaling_bad_shape
    n = size(left, 1)
    m = size(left, 2)
    if((n < 1) .or. (m < 1)) return
    if((size(right,1) /= n) .or. (size(right,2) /= m)) return
    if((size(metric,1) /= n) .or. (size(metric,2) /= n)) return
    if((size(left_balanced,1) /= n) .or. &
       (size(left_balanced,2) /= m)) return
    if((size(right_balanced,1) /= n) .or. &
       (size(right_balanced,2) /= m)) return
    if(size(left_scale) /= m) return

    left_balanced = left
    right_balanced = right
    left_scale = 1.d0

    info = scaling_nonfinite
    if(.not. all(ieee_is_finite(left))) return
    if(.not. all(ieee_is_finite(right))) return
    if(.not. all(ieee_is_finite(metric))) return

    metric_scale = max(1.d0, maxval(abs(metric)))
    metric_asymmetry = maxval(abs(metric - transpose(metric)))
    info = scaling_bad_metric
    if(metric_asymmetry > 100.d0 * epsilon(1.d0) * metric_scale) return

    allocate(overlap_before(m,m), overlap_after(m,m))
    allocate(left_trial(n,m), right_trial(n,m))
    allocate(scale_trial(m))
    left_trial = left
    right_trial = right
    scale_trial = 1.d0
    overlap_before = matmul(transpose(left), matmul(metric, right))
    info = scaling_nonfinite
    if(.not. all(ieee_is_finite(overlap_before))) return

    info = scaling_bad_overlap
    do i = 1, m
      if(abs(overlap_before(i,i)) <= tiny(1.d0)) return
    enddo

    do i = 1, m
      left_norm_sq = dot_product(left(:,i), matmul(metric, left(:,i)))
      right_norm_sq = dot_product(right(:,i), matmul(metric, right(:,i)))
      info = scaling_nonfinite
      if(.not. ieee_is_finite(left_norm_sq)) return
      if(.not. ieee_is_finite(right_norm_sq)) return

      info = scaling_bad_norm
      if((left_norm_sq <= 0.d0) .or. (right_norm_sq <= 0.d0)) return

      left_norm = sqrt(left_norm_sq)
      right_norm = sqrt(right_norm_sq)
      alpha = sqrt(right_norm) / sqrt(left_norm)
      info = scaling_nonfinite
      if((.not. ieee_is_finite(alpha)) .or. (alpha <= 0.d0)) return

      left_trial(:,i) = alpha * left(:,i)
      right_trial(:,i) = right(:,i) / alpha
      scale_trial(i) = alpha
    enddo

    if(.not. all(ieee_is_finite(left_trial))) return
    if(.not. all(ieee_is_finite(right_trial))) return

    overlap_after = matmul(transpose(left_trial), &
                           matmul(metric, right_trial))
    if(.not. all(ieee_is_finite(overlap_after))) return

    overlap_scale = max(1.d0, maxval(abs(overlap_before)))
    overlap_tolerance = 100.d0 * epsilon(1.d0) * overlap_scale
    diagonal_change = 0.d0
    offdiag_before = 0.d0
    offdiag_after = 0.d0
    do j = 1, m
      diagonal_change = max(diagonal_change, &
                            abs(overlap_after(j,j) - overlap_before(j,j)))
      do i = 1, m
        if(i == j) cycle
        offdiag_before = max(offdiag_before, abs(overlap_before(i,j)))
        offdiag_after = max(offdiag_after, abs(overlap_after(i,j)))
      enddo
    enddo

    info = scaling_overlap_changed
    if(diagonal_change > overlap_tolerance) then
      left_balanced = left
      right_balanced = right
      left_scale = 1.d0
      return
    endif
    if(offdiag_after > offdiag_before + overlap_tolerance) then
      left_balanced = left
      right_balanced = right
      left_scale = 1.d0
      return
    endif

    info = scaling_success
    left_balanced = left_trial
    right_balanced = right_trial
    left_scale = scale_trial

  end subroutine balance_biorthogonal_pairs

end module biorthogonal_scaling_module
