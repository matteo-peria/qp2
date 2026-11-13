module diis_history_module

  implicit none
  private

  integer, parameter, public :: diis_success = 0
  integer, parameter, public :: diis_bad_state = 1
  integer, parameter, public :: diis_bad_shape = 2
  integer, parameter, public :: diis_ill_conditioned = 3
  integer, parameter, public :: diis_lapack_failure = 4

  type, public :: diis_history
    private
    integer :: matrix_size = 0
    integer :: capacity = 0
    integer :: count = 0
    integer :: next_slot = 1
    double precision, allocatable :: matrix(:,:,:)
    double precision, allocatable :: error(:,:,:)
  end type diis_history

  public :: diis_initialize
  public :: diis_reset
  public :: diis_push
  public :: diis_count
  public :: diis_next_slot
  public :: diis_logical_slots
  public :: diis_get_entry
  public :: diis_build_gram
  public :: diis_solve_coefficients
  public :: diis_combine_matrix
  public :: diis_extrapolate

  interface
    subroutine dgetrf(m, n, a, lda, ipiv, info)
      integer, intent(in) :: m, n, lda
      double precision, intent(inout) :: a(lda,*)
      integer, intent(out) :: ipiv(*), info
    end subroutine dgetrf

    subroutine dgecon(norm, n, a, lda, anorm, rcond, work, iwork, info)
      character, intent(in) :: norm
      integer, intent(in) :: n, lda
      double precision, intent(in) :: a(lda,*), anorm
      double precision, intent(out) :: rcond, work(*)
      integer, intent(out) :: iwork(*), info
    end subroutine dgecon

    subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
      integer, intent(in) :: n, nrhs, lda, ldb
      double precision, intent(inout) :: a(lda,*), b(ldb,*)
      integer, intent(out) :: ipiv(*), info
    end subroutine dgesv
  end interface

contains

  subroutine diis_initialize(history, matrix_size, capacity, info)

    implicit none

    type(diis_history), intent(out) :: history
    integer, intent(in) :: matrix_size, capacity
    integer, intent(out) :: info

    if((matrix_size < 1) .or. (capacity < 1)) then
      info = diis_bad_shape
      return
    endif

    allocate(history%matrix(matrix_size,matrix_size,capacity))
    allocate(history%error(matrix_size,matrix_size,capacity))

    history%matrix_size = matrix_size
    history%capacity = capacity
    history%count = 0
    history%next_slot = 1
    history%matrix = 0.d0
    history%error = 0.d0
    info = diis_success

  end subroutine diis_initialize


  pure subroutine diis_reset(history)

    implicit none

    type(diis_history), intent(inout) :: history

    history%count = 0
    history%next_slot = 1

  end subroutine diis_reset


  pure subroutine diis_push(history, matrix, error, inserted_slot, info)

    implicit none

    type(diis_history), intent(inout) :: history
    double precision, intent(in) :: matrix(:,:), error(:,:)
    integer, intent(out) :: inserted_slot, info

    inserted_slot = 0
    if(.not. allocated(history%matrix)) then
      info = diis_bad_state
      return
    endif
    if( (size(matrix,1) /= history%matrix_size) .or. &
        (size(matrix,2) /= history%matrix_size) .or. &
        (size(error,1) /= history%matrix_size) .or. &
        (size(error,2) /= history%matrix_size) ) then
      info = diis_bad_shape
      return
    endif

    inserted_slot = history%next_slot
    history%matrix(:,:,inserted_slot) = matrix(:,:)
    history%error(:,:,inserted_slot) = error(:,:)

    history%count = min(history%count + 1, history%capacity)
    history%next_slot = modulo(inserted_slot, history%capacity) + 1
    info = diis_success

  end subroutine diis_push


  pure integer function diis_count(history)

    implicit none

    type(diis_history), intent(in) :: history

    diis_count = history%count

  end function diis_count


  pure integer function diis_next_slot(history)

    implicit none

    type(diis_history), intent(in) :: history

    diis_next_slot = history%next_slot

  end function diis_next_slot


  pure function diis_logical_slots(history) result(slots)

    implicit none

    type(diis_history), intent(in) :: history
    integer, allocatable :: slots(:)
    integer :: i, physical_slot

    allocate(slots(history%count))
    if(history%count == 0) return

    ! next_slot points to the oldest entry once the circular buffer is full.
    physical_slot = modulo( history%next_slot - history%count - 1, &
                            history%capacity ) + 1
    do i = 1, history%count
      slots(i) = physical_slot
      physical_slot = modulo(physical_slot, history%capacity) + 1
    enddo

  end function diis_logical_slots


  pure subroutine diis_get_entry(history, logical_index, matrix, error, info)

    implicit none

    type(diis_history), intent(in) :: history
    integer, intent(in) :: logical_index
    double precision, intent(out) :: matrix(:,:), error(:,:)
    integer, intent(out) :: info
    integer, allocatable :: slots(:)
    integer :: physical_slot

    if( (logical_index < 1) .or. &
        (logical_index > history%count) ) then
      info = diis_bad_state
      return
    endif
    if( (size(matrix,1) /= history%matrix_size) .or. &
        (size(matrix,2) /= history%matrix_size) .or. &
        (size(error,1) /= history%matrix_size) .or. &
        (size(error,2) /= history%matrix_size) ) then
      info = diis_bad_shape
      return
    endif

    slots = diis_logical_slots(history)
    physical_slot = slots(logical_index)
    matrix(:,:) = history%matrix(:,:,physical_slot)
    error(:,:) = history%error(:,:,physical_slot)
    info = diis_success

  end subroutine diis_get_entry


  pure subroutine diis_build_gram(history, gram, info)

    implicit none

    type(diis_history), intent(in) :: history
    double precision, allocatable, intent(out) :: gram(:,:)
    integer, intent(out) :: info
    integer, allocatable :: slots(:)
    integer :: i, j

    if(history%count < 1) then
      allocate(gram(0,0))
      info = diis_bad_state
      return
    endif

    slots = diis_logical_slots(history)
    allocate(gram(history%count,history%count))
    do j = 1, history%count
      do i = 1, history%count
        gram(i,j) = sum( history%error(:,:,slots(i)) * &
                         history%error(:,:,slots(j)) )
      enddo
    enddo
    info = diis_success

  end subroutine diis_build_gram


  subroutine diis_solve_coefficients(gram, coefficients, rcond, info, rcond_min)

    implicit none

    double precision, intent(in) :: gram(:,:)
    double precision, allocatable, intent(out) :: coefficients(:)
    double precision, intent(out) :: rcond
    integer, intent(out) :: info
    double precision, intent(in), optional :: rcond_min
    double precision, allocatable :: augmented(:,:), factorized(:,:)
    double precision, allocatable :: rhs(:), work(:)
    integer, allocatable :: ipiv(:), iwork(:)
    integer :: i, lapack_info, n, system_size
    double precision :: anorm, minimum_rcond

    n = size(gram,1)
    allocate(coefficients(n))
    coefficients = 0.d0
    rcond = 0.d0

    if((n < 1) .or. (size(gram,2) /= n)) then
      info = diis_bad_shape
      return
    endif

    minimum_rcond = 1.d-14
    if(present(rcond_min)) minimum_rcond = rcond_min

    system_size = n + 1
    allocate(augmented(system_size,system_size))
    allocate(factorized(system_size,system_size))
    allocate(rhs(system_size), ipiv(system_size))
    allocate(work(4*system_size), iwork(system_size))

    augmented = 0.d0
    augmented(1:n,1:n) = gram(:,:)
    augmented(1:n,system_size) = -1.d0
    augmented(system_size,1:n) = -1.d0
    rhs = 0.d0
    rhs(system_size) = -1.d0

    anorm = 0.d0
    do i = 1, system_size
      anorm = max(anorm, sum(abs(augmented(:,i))))
    enddo

    factorized(:,:) = augmented(:,:)
    call dgetrf( system_size, system_size, factorized, system_size, &
                 ipiv, lapack_info )
    if(lapack_info /= 0) then
      info = diis_lapack_failure
      return
    endif

    call dgecon( '1', system_size, factorized, system_size, anorm, &
                 rcond, work, iwork, lapack_info )
    if(lapack_info /= 0) then
      info = diis_lapack_failure
      return
    endif
    if(rcond < minimum_rcond) then
      info = diis_ill_conditioned
      return
    endif

    call dgesv( system_size, 1, augmented, system_size, ipiv, rhs, &
                system_size, lapack_info )
    if(lapack_info /= 0) then
      info = diis_lapack_failure
      return
    endif

    coefficients(:) = rhs(1:n)
    info = diis_success

  end subroutine diis_solve_coefficients


  pure subroutine diis_combine_matrix(history, coefficients, matrix, info)

    implicit none

    type(diis_history), intent(in) :: history
    double precision, intent(in) :: coefficients(:)
    double precision, intent(out) :: matrix(:,:)
    integer, intent(out) :: info
    integer, allocatable :: slots(:)
    integer :: i

    if( (size(coefficients) /= history%count) .or. &
        (size(matrix,1) /= history%matrix_size) .or. &
        (size(matrix,2) /= history%matrix_size) ) then
      info = diis_bad_shape
      return
    endif

    slots = diis_logical_slots(history)
    matrix = 0.d0
    do i = 1, history%count
      matrix(:,:) = matrix(:,:) + coefficients(i) * &
                  history%matrix(:,:,slots(i))
    enddo
    info = diis_success

  end subroutine diis_combine_matrix


  subroutine diis_extrapolate(history, matrix, coefficients, rcond, info, rcond_min)

    implicit none

    type(diis_history), intent(in) :: history
    double precision, intent(out) :: matrix(:,:)
    double precision, allocatable, intent(out) :: coefficients(:)
    double precision, intent(out) :: rcond
    integer, intent(out) :: info
    double precision, intent(in), optional :: rcond_min
    double precision, allocatable :: gram(:,:)

    call diis_build_gram(history, gram, info)
    if(info /= diis_success) return

    if(present(rcond_min)) then
      call diis_solve_coefficients( gram, coefficients, rcond, &
                                       info, rcond_min )
    else
      call diis_solve_coefficients( gram, coefficients, rcond, info )
    endif
    if(info /= diis_success) return

    call diis_combine_matrix(history, coefficients, matrix, info)

  end subroutine diis_extrapolate

end module diis_history_module
