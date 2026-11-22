
! ---

BEGIN_PROVIDER [double precision, Q_alpha, (ao_num, ao_num) ]

  BEGIN_DOC
  !
  ! Q_alpha = mo_r_coef x eta_occ_alpha x mo_l_coef.T
  !
  ! [Q_alpha]_ij = \sum_{k=1}^{elec_alpha_num} [mo_r_coef]_ik [mo_l_coef]_jk
  !
  END_DOC

  implicit none

  Q_alpha = 0.d0
  call dgemm( 'N', 'T', ao_num, ao_num, elec_alpha_num, 1.d0               &
            , mo_r_coef, size(mo_r_coef, 1), mo_l_coef, size(mo_l_coef, 1) &
            , 0.d0, Q_alpha, size(Q_alpha, 1) )

END_PROVIDER
  
! ---
    
BEGIN_PROVIDER [ double precision, Q_beta, (ao_num, ao_num) ]

  BEGIN_DOC
  !
  ! Q_beta = mo_r_coef x eta_occ_beta x mo_l_coef.T
  !
  ! [Q_beta]_ij = \sum_{k=1}^{elec_beta_num} [mo_r_coef]_ik [mo_l_coef]_jk
  !
  END_DOC

  implicit none

  Q_beta = 0.d0
  call dgemm( 'N', 'T', ao_num, ao_num, elec_beta_num, 1.d0                &
            , mo_r_coef, size(mo_r_coef, 1), mo_l_coef, size(mo_l_coef, 1) &
            , 0.d0, Q_beta, size(Q_beta, 1) )

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, Q_matrix, (ao_num, ao_num) ]

  BEGIN_DOC
  !
  ! Q_matrix = 2 mo_r_coef x eta_occ x mo_l_coef.T
  ! 
  ! with: 
  !                        | 1   if i = j = 1, ..., nb of occ orbitals
  !        [eta_occ]_ij =  |     
  !                        | 0   otherwise
  !
  ! the diis error is defines as:
  !                         e = F_ao x Q x ao_overlap - ao_overlap x Q x F_ao
  ! with: 
  !       mo_l_coef.T x ao_overlap x mo_r_coef = I
  !       F_mo = mo_l_coef.T x F_ao x mo_r_coef
  !       F_ao = (ao_overlap x mo_r_coef) x F_mo x (ao_overlap x mo_l_coef).T
  !
  ! ==> e = 2 ao_overlap x mo_r_coef x [ F_mo x eta_occ - eta_occ x F_mo ] x (ao_overlap x mo_l_coef).T
  !
  !      at convergence:
  !                                      F_mo x eta_occ - eta_occ x F_mo = 0
  !                                  ==> [F_mo]_ij ([eta_occ]_ii - [eta_occ]_jj) = 0  
  !                                  ==> [F_mo]_ia = [F_mo]_ai = 0 where: i = occ and a = vir
  !                                  ==> Brillouin conditions
  !
  END_DOC

  implicit none

  if(elec_alpha_num == elec_beta_num) then
    Q_matrix = Q_alpha + Q_alpha
  else
    Q_matrix = Q_alpha + Q_beta
  endif

END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, FQS_SQF_ao, (ao_num, ao_num)]

  implicit none
  integer                       :: i, j
  double precision              :: t0, t1
  double precision, allocatable :: tmp(:,:)
  double precision, allocatable :: F(:,:)

  PROVIDE Fock_matrix_tc_ao_tot

  allocate(F(ao_num,ao_num))
  do i = 1, ao_num
    do j = 1, ao_num
      F(j,i) = Fock_matrix_tc_ao_tot(j,i)
    enddo
  enddo

  allocate(tmp(ao_num,ao_num))

  ! F x Q
  call dgemm( 'N', 'N', ao_num, ao_num, ao_num, 1.d0     &
            , F, size(F, 1), Q_matrix, size(Q_matrix, 1) &
            , 0.d0, tmp, size(tmp, 1) )

  ! F x Q x S
  call dgemm( 'N', 'N', ao_num, ao_num, ao_num, 1.d0             &
            , tmp, size(tmp, 1), ao_overlap, size(ao_overlap, 1) &
            , 0.d0, FQS_SQF_ao, size(FQS_SQF_ao, 1) )

  ! S x Q
  tmp = 0.d0
  call dgemm( 'N', 'N', ao_num, ao_num, ao_num, 1.d0                       &
            , ao_overlap, size(ao_overlap, 1), Q_matrix, size(Q_matrix, 1) &
            , 0.d0, tmp, size(tmp, 1) )

  ! F x Q x S - S x Q x F
  call dgemm( 'N', 'N', ao_num, ao_num, ao_num, -1.d0 &
            , tmp, size(tmp, 1), F, size(F, 1)        &
            , 1.d0, FQS_SQF_ao, size(FQS_SQF_ao, 1) )

  deallocate(tmp)
  deallocate(F)

END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, FQS_SQF_mo, (mo_num, mo_num)]

  implicit none
  double precision :: t0, t1

  PROVIDE mo_r_coef mo_l_coef
  PROVIDE FQS_SQF_ao

  call ao_to_mo_bi_ortho( FQS_SQF_ao, size(FQS_SQF_ao, 1) &
                        , FQS_SQF_mo, size(FQS_SQF_mo, 1) )

END_PROVIDER

! ---


subroutine tc_orthogonalize_span(coef, x)
  implicit none
  double precision, intent(in) :: coef(ao_num,mo_num)
  double precision, intent(out) :: x(ao_num,mo_num)
  double precision :: scale
  double precision, allocatable :: gram(:,:)
  integer :: i, info

  ! Strip arbitrary column scales before forming the overlap matrix.
  do i=1,mo_num
    scale = maxval(abs(coef(:,i)))
    if (scale == 0.d0) error stop 'Empty TC orbital column'
    x(:,i) = coef(:,i)/scale
  enddo
  allocate(gram(mo_num,mo_num))
  gram = matmul(transpose(x),matmul(ao_overlap,x))
  call dpotrf('U',mo_num,gram,mo_num,info)
  if (info /= 0) error stop 'Cannot orthonormalize truncated TC orbital space'
  call dtrsm('R','U','N','N',ao_num,mo_num,1.d0,gram,mo_num,x,ao_num)
end

subroutine tc_ao_residual_norm(x_left, x_right, residual)
  use ieee_arithmetic, only: ieee_is_finite
  implicit none
  double precision, intent(in) :: x_left(ao_num,mo_num), x_right(ao_num,mo_num)
  double precision, intent(out) :: residual
  double precision, allocatable :: projected(:,:), singular(:), work(:)
  double precision :: unused_u(1), unused_v(1)
  integer :: info

  allocate(projected(mo_num,mo_num), singular(mo_num), work(5*mo_num))
  projected = matmul(transpose(x_left),matmul(FQS_SQF_ao,x_right))
  if (.not.all(ieee_is_finite(projected))) error stop 'Nonfinite AO residual'
  ! Only the singular values are needed; LAPACK returns the largest first.
  call dgesvd('N','N',mo_num,mo_num,projected,mo_num,singular, &
              unused_u,1,unused_v,1,work,size(work),info)
  if (info /= 0) error stop 'AO residual singular value calculation failed'
  residual = singular(1)
end
