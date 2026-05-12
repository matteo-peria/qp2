
subroutine non_hrmt_bieig_degen(n, A, thr_d, thr_nd, leigvec, reigvec, n_real_eigv, eigval)

  BEGIN_DOC
  ! 
  ! routine which returns the EIGENVALUES and corresponding LEFT/RIGHT eigenvetors 
  ! of a non hermitian matrix A(n,n)
  !
  ! n_real_eigv is the number of real eigenvalues, which might be smaller than the dimension "n" 
  !
  ! In the case where some eigenvalues have a complex component, 
  ! 
  ! H Phi_1 = (epsilon_r + i epsilon_i) Phi_1, H Phi_2 = (epsilon_r - i epsilon_i) Phi_2
  ! 
  ! Phi_1 = Phi_r + i Phi_i, Phi_2 = Phi_r - i Phi_i, 
  ! 
  ! Xhi_1 = Xhi_r + i Xhi_i, Xhi_2 = Xhi_r - i Xhi_i, 
  ! 
  ! it returns ONLY the real-part of the eigenvalues (i.e. epsilon_r) 
  !
  ! and the eigenvectors are written as 
  !
  ! Phi_1 = Phi_r + Phi_i, Phi_2 = Phi_r - Phi_i, 
  ! 
  ! Xhi_1 = Xhi_r + Xhi_i, Xhi_2 = Xhi_r - Xhi_i, 
  !
  ! This choice of vectors corresponds to imposing 
  ! 
  ! (Xhi_i|H|Phi_j) = epsilon_r \delta_ij, and 
  ! 
  ! (Xhi_i|Phi_j) = \delta_ij.
  END_DOC

  implicit none
  integer,          intent(in)  :: n
  double precision, intent(in)  :: A(n,n)
  double precision, intent(in)  :: thr_d, thr_nd
  integer,          intent(out) :: n_real_eigv
  double precision, intent(out) :: reigvec(n,n), leigvec(n,n), eigval(n)

  integer                       :: i, j,k 
  integer                       :: n_real, n_im
  double precision              :: thr, thr_cut, thr_diag, thr_norm
  double precision              :: accu_d, accu_nd

  integer,          allocatable :: list_real(:), iorder_real(:), deg_num(:), iorder(:)
  integer,          allocatable :: list_im(:), iorder_im(:)
  double precision, allocatable :: WR(:), WI(:), VL(:,:), VR(:,:), eigv_real(:), eigv_im(:), eigv_im_im(:)
  double precision, allocatable :: S(:,:), eigv_srtd_real(:), eigv_srtd_im(:)
  double precision, allocatable :: reigvec_srtd_real(:,:), reigvec_srtd_im(:,:)
  double precision, allocatable :: leigvec_srtd_real(:,:), leigvec_srtd_im(:,:)
  double precision, allocatable  :: phi_1_tilde(:),phi_2_tilde(:),chi_1_tilde(:),chi_2_tilde(:)
  logical, allocatable :: is_real(:)

  allocate(phi_1_tilde(n),phi_2_tilde(n),chi_1_tilde(n),chi_2_tilde(n))

  allocate(WR(n), WI(n), VL(n,n), VR(n,n), is_real(n),iorder(n) ) 

  call lapack_diag_non_sym(n, A, WR, WI, VL, VR)
!  

  thr_diag = 1d-06
  thr_norm = 1d+10

  ! ---

  ! track & sort the real eigenvalues 

  allocate(list_real(n), list_im(n), eigv_real(n), eigv_im(n),eigv_im_im(n))
  thr    = 0.d0
  n_real = 0
  n_im   = 0
  do i = 1, n
    if(dabs(WI(i)) .le. thr) then
      is_real(i) = .True.
      n_real +=1 
      list_real(n_real) = i
      eigv_real(n_real) = WR(i)
    else
      print*, 'Found an imaginary component to eigenvalue on i = ', i
      print*, 'Re(i) + Im(i)', WR(i), WI(i)
      is_real(i) = .False.
      n_im +=1 
      eigv_im(n_im) = WR(i)
      eigv_im_im(n_im) = WI(i)
      list_im(n_im) = i
    endif
  enddo
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! REAL PART ONLY

 ! sort the real eigenvalues and copy the eigenvalues and eigenvectors 
 ! in eigv_srtd_real and (leigvec_srtd_real , reigvec_srtd_real)
  allocate(leigvec_srtd_real(n,n_real),reigvec_srtd_real(n,n_real),eigv_srtd_real(n_real)) 
  do i = 1, n_real
   iorder(i) = list_real(i)
  enddo
  call dsort(eigv_real, iorder, n_real)
 !print*,'In real part'
  do i = 1, n_real
   eigv_srtd_real(i) = eigv_real(i)
   do j = 1, n
    reigvec_srtd_real(j,i) = VR(j,iorder(i))
    leigvec_srtd_real(j,i) = VL(j,iorder(i))
   enddo
  !write(*,'(100(F16.10,X))')reigvec_srtd_real(1:n,i)
  !write(*,'(100(F16.10,X))')leigvec_srtd_real(1:n,i)
  enddo

  if(n_im>0)then
   print*,'Number of imaginary eigenvalues found = ',n_im
    n_im_even = n_im/2
    integer :: n_im_even, ii
    complex*16, allocatable :: right_vec(:,:), left_vec(:,:), smatrix(:,:)
    complex*16, allocatable :: inv_right_vec(:,:)
    double precision, allocatable :: eigval_tmp(:), Phi_r(:), Phi_i(:), Xhi_r(:), Xhi_i(:),sphichi(:,:)
    double precision :: ddot
    allocate(Phi_r(n), Phi_i(n), Xhi_r(n), Xhi_i(n),eigval_tmp(n_im_even))
    allocate(right_vec(n,2), left_vec(n,2), smatrix(n,2),inv_right_vec(2,n),sphichi(2,2))
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! IMAGINARY COMPONENTS 
   ! sort the complex eigenvalues by their real part and copy the eigenvalues and eigenvectors 
   ! in eigv_srtd_im and (leigvec_srtd_im , reigvec_srtd_im)
    if (iand(n_im,1)==1)then
     call qp_bug(irp_here,n_im,'n_im is odd ')
    endif
    allocate(leigvec_srtd_im(n,n_im),reigvec_srtd_im(n,n_im),eigv_srtd_im(n_im)) 
    do i = 1, n_im, 2
     iorder(i/2+1) = list_im(i)
     eigval_tmp(i/2+1) = eigv_im(i)
    enddo
    call dsort(eigval_tmp, iorder, n_im_even)
    ! first couple 
    do i = 1, n_im
     eigv_srtd_im(i) = eigval_tmp((i-1)/2+1)
    enddo
    do ii = 1,n_im/2
     print*,'eigenvalue number ',ii,eigval_tmp(ii)
     right_vec(1:n,1)   = cmplx(VR(1:n,(iorder(ii))) , VR(1:n,(iorder(ii))+1))
     right_vec(1:n,2)   = cmplx(VR(1:n,(iorder(ii))) ,-1.d0* VR(1:n,(iorder(ii))+1))
!      call get_inverse_complex(right_vec,n,n,inv_right_vec,n)
     call get_pseudo_inverse_complex(right_vec,size(right_vec,1),n,2,inv_right_vec,size(inv_right_vec,1),0.d0)
     ! need to take complex conjug because of convention of L^2 inner in get_pseudo_inverse_complex
     inv_right_vec = conjg(inv_right_vec)
     do i = 1, 2
      do j = 1, n
       left_vec(j,i) = inv_right_vec(i,j)
      enddo
     enddo
       
     smatrix = cmplx(0.d0,0.d0)
     do i= 1, 2
      do j = 1, 2
       do k = 1, n
        smatrix(j,i) += (left_vec(k,j)) * right_vec(k,i)
       enddo
      enddo
     enddo
     double precision :: accu
     accu = abs(smatrix(1,2))+abs(smatrix(2,1))
     if(dabs(accu).gt.1.d-10)then
      call qp_bug(irp_here,int(dlog(accu)), & 
             'the 2x2 smatrix of complex-valued eigenvectors is not diagonal')
     endif
     accu = real(smatrix(1,1)+smatrix(2,2),8)
     if(dabs(2.d0-accu).gt.1.d-10)then
      call qp_bug(irp_here,int(dlog(dabs(2.d0-accu))), & 
             'the trace of the 2x2 smatrix of complex-valued eigenvectors is not 2')
     endif
     accu = aimag(smatrix(1,1)+smatrix(2,2))
     if(dabs(accu).gt.1.d-10)then
      call qp_bug(irp_here,int(dlog(dabs(accu))), & 
             'the trace of the 2x2 smatrix of complex-valued eigenvectors is not 2')
     endif
   
!     print*,'smatrix = '
!     do i = 1, 2
!      write(*,'(100(F6.3,SP,F6.3,X))')smatrix(i,1:2)
!     enddo
     Phi_r(1:n) = real(right_vec(1:n,1),8)
     Phi_i(1:n) = aimag(right_vec(1:n,1))
     Xhi_r(1:n) = real(left_vec(1:n,1),8)
     Xhi_i(1:n) = aimag(left_vec(1:n,1))
     sphichi = 0.d0 
     sphichi(1,1) = ddot(n,Xhi_r,1,Phi_r,1)
     sphichi(2,1) = ddot(n,Xhi_i,1,Phi_r,1)
     sphichi(1,2) = ddot(n,Xhi_r,1,Phi_i,1)
     sphichi(2,2) = ddot(n,Xhi_i,1,Phi_i,1)

     accu = dabs(sphichi(1,2))+dabs(sphichi(2,1))
     if(dabs(accu).gt.1.d-10)then
      call qp_bug(irp_here,int(dlog(accu)), & 
             'the 2x2 sphichi of complex-valued eigenvectors is not diagonal')
     endif
     accu = sphichi(1,1)+sphichi(2,2)
     if(dabs(accu).gt.1.d-10)then
      call qp_bug(irp_here,int(dlog(dabs(2.d0-accu))), & 
             'the trace of the 2x2 sphichi of complex-valued eigenvectors is not 0')
     endif
     if(dabs(sphichi(1,1)-0.5d0).gt.1.d-10)then
      call qp_bug(irp_here,int(dlog(dabs(sphichi(1,1)-0.5d0))), & 
             'the diagonal element of the 2x2 sphichi of complex-valued eigenvectors is not 1/2')
     endif

!     print*,'Smatrix sphichi'
!     do i = 1, 2
!      write(*,'(100(F16.10,X))')sphichi(i,:)
!     enddo
      ! Phi_1 = Phi_r + Phi_i
      reigvec_srtd_im(1:n,2*ii-1)   = Phi_r(1:n)
      reigvec_srtd_im(1:n,2*ii-1)  += Phi_i(1:n)
      ! Phi_2 = Phi_r - Phi_i
      reigvec_srtd_im(1:n,2*ii) = Phi_r(1:n)
      reigvec_srtd_im(1:n,2*ii)-= Phi_i(1:n)
      
      ! Xhi_1 = Xhi_r - Xhi_i
      leigvec_srtd_im(1:n,2*ii-1)   = Xhi_r(1:n)
      leigvec_srtd_im(1:n,2*ii-1)  -= Xhi_i(1:n)
      ! Xhi_2 = Xhi_r + Xhi_i
      leigvec_srtd_im(1:n,2*ii) = Xhi_r(1:n)
      leigvec_srtd_im(1:n,2*ii)+= Xhi_i(1:n)
   enddo
  endif
  
  !! Sorting the union of the real and imaginary eigenvalues 
  double precision, allocatable :: real_eigv_total(:), reigvec_total(:,:), leigvec_total(:,:)
  allocate(real_eigv_total(n), reigvec_total(n,n), leigvec_total(n,n))
  do i = 1, n_real
   real_eigv_total(i) = eigv_real(i)
   reigvec_total(1:n,i) = reigvec_srtd_real(1:n,i)
   leigvec_total(1:n,i) = leigvec_srtd_real(1:n,i)
  enddo
  do ii = 1, n_im
   i=n_real + ii
   real_eigv_total(i) = eigv_srtd_im(ii) 
   reigvec_total(1:n,i) = reigvec_srtd_im(1:n,ii)
   leigvec_total(1:n,i) = leigvec_srtd_im(1:n,ii)
  enddo
  do i = 1, n
   iorder(i) = i
  enddo
  call dsort(real_eigv_total, iorder, n)
  do i = 1, n
   eigval(i) = real_eigv_total(i)
   print*,i,eigval(i)
   reigvec(1:n,i) = reigvec_total(1:n,iorder(i))
   leigvec(1:n,i) = leigvec_total(1:n,iorder(i))
  enddo
  double precision, allocatable :: leigvec_tmp(:,:)
  allocate(leigvec_tmp(n,n))
!  call get_inverse(reigvec,size(reigvec,1),n,leigvec_tmp,size(leigvec_tmp,1))
  call get_pseudo_inverse(reigvec,size(reigvec,1),n,n,leigvec_tmp,size(leigvec_tmp,1),0.d0)
  do i = 1, n
   do j =1, n
    leigvec(i,j) = leigvec_tmp(j,i)
   enddo
  enddo
 

!  ! check bi-orthogonality
!
  thr_diag = 10.d0
  thr_norm = 1d+10
!
  n_real_eigv = n
  allocate( S(n_real_eigv,n_real_eigv) )
  call check_biorthog(n, n_real_eigv, leigvec, reigvec, accu_d, accu_nd, S, thr_d, thr_nd, .false.)
!
  if( (accu_nd .lt. thr_nd) .and. (dabs(accu_d-dble(n_real_eigv))/dble(n_real_eigv) .lt. thr_d) ) then

    print *, ' lapack vectors are normalized and bi-orthogonalized'
    deallocate(S)
    return

  ! accu_nd is modified after adding the normalization
  elseif( (accu_nd .lt. thr_nd) .and. (dabs(accu_d-dble(n_real_eigv))/dble(n_real_eigv) .gt. thr_d) ) then

    print *, ' lapack vectors are not normalized but bi-orthogonalized'
    call check_biorthog_binormalize(n, n_real_eigv, leigvec, reigvec, thr_d, thr_nd, .true.)

    call check_biorthog(n, n_real_eigv, leigvec, reigvec, accu_d, accu_nd, S, thr_d, thr_nd, .true.)
    call check_EIGVEC(n, n, A, eigval, leigvec, reigvec, thr_diag, thr_norm, .true.)

    deallocate(S)
    return

  else

    print *, ' lapack vectors are not normalized neither bi-orthogonalized'

    allocate(deg_num(n))
    call reorder_degen_eigvec(n, deg_num, eigval, leigvec, reigvec)
    call impose_biorthog_degen_eigvec(n, deg_num, eigval, leigvec, reigvec)
    deallocate(deg_num)

    call check_biorthog(n, n_real_eigv, leigvec, reigvec, accu_d, accu_nd, S, thr_d, thr_nd, .false.)
    if( (accu_nd .lt. thr_nd) .and. (dabs(accu_d-dble(n_real_eigv)) .gt. thr_d) ) then
      call check_biorthog_binormalize(n, n_real_eigv, leigvec, reigvec, thr_d, thr_nd, .true.)
    endif
    call check_biorthog(n, n_real_eigv, leigvec, reigvec, accu_d, accu_nd, S, thr_d, thr_nd, .true.)

    deallocate(S)

  endif

  return

end
! ---

subroutine non_hrmt_bieig(n, A, thr_d, thr_nd, leigvec, reigvec, n_real_eigv, eigval)

  BEGIN_DOC
  ! 
  ! routine which returns the sorted REAL EIGENVALUES ONLY and corresponding LEFT/RIGHT eigenvetors 
  ! of a non hermitian matrix A(n,n)
  !
  ! n_real_eigv is the number of real eigenvalues, which might be smaller than the dimension "n" 
  !
  END_DOC

  implicit none
  integer,          intent(in)  :: n
  double precision, intent(in)  :: A(n,n)
  double precision, intent(in)  :: thr_d, thr_nd
  integer,          intent(out) :: n_real_eigv
  double precision, intent(out) :: reigvec(n,n), leigvec(n,n), eigval(n)

  integer                       :: i, j,k 
  integer                       :: n_good
  double precision              :: thr, thr_cut, thr_diag, thr_norm
  double precision              :: accu_d, accu_nd

  integer,          allocatable :: list_good(:), iorder(:), deg_num(:)
  double precision, allocatable :: WR(:), WI(:), VL(:,:), VR(:,:)
  double precision, allocatable :: S(:,:)
  double precision, allocatable  :: phi_1_tilde(:),phi_2_tilde(:),chi_1_tilde(:),chi_2_tilde(:)

  allocate(phi_1_tilde(n),phi_2_tilde(n),chi_1_tilde(n),chi_2_tilde(n))

  allocate(WR(n), WI(n), VL(n,n), VR(n,n)) 

  call lapack_diag_non_sym(n, A, WR, WI, VL, VR)

  thr_diag = 1d-06
  thr_norm = 1d+10

  ! ---

  ! track & sort the real eigenvalues 

  n_good = 0
  thr    = Im_thresh_tc
  do i = 1, n
    if(dabs(WI(i)) .lt. thr) then
      n_good += 1
    else
      print*, 'Found an imaginary component to eigenvalue on i = ', i
      print*, 'Re(i) + Im(i)', WR(i), WI(i)
    endif
  enddo

  if(n_good.ne.n) then
    print*,'there are some imaginary eigenvalues '
    thr_diag = 1d-03
    n_good = n
  endif

  allocate(list_good(n_good), iorder(n_good))

  n_good = 0
  do i = 1, n
    n_good += 1
    list_good(n_good) = i
    eigval(n_good) = WR(i)
  enddo

  deallocate( WR, WI )

  n_real_eigv = n_good 
  do i = 1, n_good
    iorder(i) = i
  enddo
  call dsort(eigval, iorder, n_good)
      
  reigvec(:,:) = 0.d0 
  leigvec(:,:) = 0.d0 
  do i = 1, n_real_eigv
    do j = 1, n
      reigvec(j,i) = VR(j,list_good(iorder(i)))
      leigvec(j,i) = VL(j,list_good(iorder(i)))
    enddo
  enddo

  deallocate( list_good, iorder )
  deallocate( VL, VR )

  ASSERT(n==n_real_eigv)

  ! ---

  ! check bi-orthogonality

  thr_diag = 10.d0
  thr_norm = 1d+10

  allocate( S(n_real_eigv,n_real_eigv) )
  call check_biorthog(n, n_real_eigv, leigvec, reigvec, accu_d, accu_nd, S, thr_d, thr_nd, .false.)

  if( (accu_nd .lt. thr_nd) .and. (dabs(accu_d-dble(n_real_eigv))/dble(n_real_eigv) .lt. thr_d) ) then

    print *, ' lapack vectors are normalized and bi-orthogonalized'
    deallocate(S)
    return

  ! accu_nd is modified after adding the normalization
  elseif( (accu_nd .lt. thr_nd) .and. (dabs(accu_d-dble(n_real_eigv))/dble(n_real_eigv) .gt. thr_d) ) then

    print *, ' lapack vectors are not normalized but bi-orthogonalized'
    call check_biorthog_binormalize(n, n_real_eigv, leigvec, reigvec, thr_d, thr_nd, .true.)

    call check_biorthog(n, n_real_eigv, leigvec, reigvec, accu_d, accu_nd, S, thr_d, thr_nd, .true.)
    call check_EIGVEC(n, n, A, eigval, leigvec, reigvec, thr_diag, thr_norm, .true.)

    deallocate(S)
    return

  else

    print *, ' lapack vectors are not normalized neither bi-orthogonalized'

    allocate(deg_num(n))
    call reorder_degen_eigvec(n, deg_num, eigval, leigvec, reigvec)
    call impose_biorthog_degen_eigvec(n, deg_num, eigval, leigvec, reigvec)
    deallocate(deg_num)

    call check_biorthog(n, n_real_eigv, leigvec, reigvec, accu_d, accu_nd, S, thr_d, thr_nd, .false.)
    if( (accu_nd .lt. thr_nd) .and. (dabs(accu_d-dble(n_real_eigv)) .gt. thr_d) ) then
      call check_biorthog_binormalize(n, n_real_eigv, leigvec, reigvec, thr_d, thr_nd, .true.)
    endif
    call check_biorthog(n, n_real_eigv, leigvec, reigvec, accu_d, accu_nd, S, thr_d, thr_nd, .true.)

    deallocate(S)

  endif

  return

end

! ---

subroutine check_bi_ortho(reigvec, leigvec, n, S, accu_nd)

  BEGIN_DOC
  ! retunrs the overlap matrix S = Leigvec^T Reigvec 
  !
  ! and the square root of the sum of the squared off-diagonal elements of S
  END_DOC

  implicit none
  integer,          intent(in)  :: n
  double precision, intent(in)  :: reigvec(n,n), leigvec(n,n)
  double precision, intent(out) :: S(n,n), accu_nd

  integer :: i,j

  ! S = VL x VR
  call dgemm( 'T', 'N', n, n, n, 1.d0 &
            , leigvec, size(leigvec, 1), reigvec, size(reigvec, 1)  &
            , 0.d0, S, size(S, 1) )

  accu_nd = 0.d0
  do i = 1, n
    do j = 1, n
      if(i.ne.j) then
        accu_nd = accu_nd + S(j,i) * S(j,i)
      endif
    enddo
  enddo
  accu_nd = dsqrt(accu_nd)

end


