! +--INTDIFF----------------------------------------------------------+
! | task: integration and differentiation                             |
! |       build-up the Jacobian for every single experiment           |
! +-------------------------------------------------------------------+

!> \routine intdiff
!> \brief integration and differentiation routine
!! bulid up of the Jacobian, right hand side and seed matrices
!> \author Robert Kircheis

subroutine intdiff(iexp, T, MSTIMES, STATES, PARAM, SR, ndiff, nvar, dimr, &
                   dimrloc, multidimr, dims, nmess, nbed, nt, nms, nsp, &
                   JACOBIAN, RHS, MS, MP, MR, cont, least, CONTR, nglob, IBEL, &
                   maxnb, maxsr, maxnv, maxms, infor)

  implicit none

  integer :: iexp, &
             ndiff, &
             nvar, &
             dimr, &
             dimrloc, &
             multidimr, &
             dims, &
             nmess, &
             nbed, &
             nt, &
             nms, &
             nsp, &
             nglob, &
             IBEL(maxnv,maxms), &
             maxnb, &
             maxsr, &
             maxnv, &
             maxms, & 
             infor
  real*8 :: T(nt), &
            MSTIMES(nms), &
            STATES(nvar,nms), &
            PARAM(nglob), &
            SR(dimr), &
            JACOBIAN(maxnb,maxsr+nglob), &
            RHS(maxnb), &
            MS(maxnv,maxsr,maxms), &
            MP(maxnv,nglob,maxms), &
            MR(maxnv,maxms), &
            cont, &
            least, &
            contr(maxsr)

  integer :: i,j,k,l,ii, & !< numerator
             idum, & !< auxiliary variable
             iflag, & !< error flag
             ld, & !< flag if shooting node
             n, & !< nvar + nglob
             ndir, & !< number of directions; dimr + nglob + 1
             nalg, & !< number of algebraic variables
             nfix, & !< number of fixed state variables
             ndx !< auxiliary variable
  real*8 :: rdum, & !< ???
            s, & !< auxiliary variable
            hmax, h, hopt, & !< stepsizes controls for integration
            dummy !< auxiliary variable
  real*8, allocatable :: GA0(:,:), & !< block of the Jacobian for the computation
                                     !! of the initial seed matrices; dim(nvar,nalg+nglob+dimr)
                         D0(:), & !< block of rhs for the computation
                                  !! of the initial seed matrices; dim(nvar)
                         XDIR(:,:), & !< unit seed matrix for x-direction;
                                      !! dim(nvar,nvar)
                         PDIR(:,:), & !< unit seed matrix for p-direction;
                                      !! dim(nglob,nglob)
                         SDIR(:,:), & !< unit seed matrix for sr-direction;
                                      !! dim(dims,dims+nalg)
                         SZR(:), & !< complete s-vector (s_z,s_r)
                         TAU(:), & !< auxiliary array for QR-decomposition; dim(nvar)
                         X(:), & !< trajectory at current time step; dim(nvar)
                         U(:), & !< auxiliary array; dim(nvar)
                         MST(:,:), & !< MS transposed; dim(dims,nvar)
                         MPT(:,:), & !< MP transposed; dim(nglob,nvar)
                         GA(:,:), & !< directional derivatives of the Wronskian matrix;
                                    !! dim(nvar,ndir)
                         GADIR(:,:), & !< directions of the Wronskian matrix GA;
                                       !! dim(n,ndir)
                         ZI(:), & !< residium for the measurement; dim(nmess)
                         DMDX(:,:), & !< derivatives of the mfcns wrt x;
                                      !! dim(nmess,ndx)
                         DMDP(:,:), & !< derivatives of the mfcns wrt p;
                                      !! dim(nmess,nglob)
                         GR(:,:), & !< derivatives of the rfcns wrt p, s and x;
                                    !! dim(dimr,max(dims,nglob)
                         R(:) !< return value of rfcn; dim(dimr)
  integer, allocatable :: NODE(:) !< flag for node type; dim(3)

!   call dinit2(nbed*(dims+nglob), 0.D0, JACOBIAN,1)
!   call dinit2(nbed, 0.D0, RHS, 1)
!   call dinit2(nvar*dims*nms, 0.D0, MS, 1)
!   call dinit2(nvar*nglob*nms,0.D0, MP, 1)
!   call dinit2(nvar*nms, 0.D0, MR, 1)

  n    = nvar + nglob
  ndir = nglob + dims + 1
  nalg = nvar-ndiff
  hmax = DABS( MSTIMES(nms) - MSTIMES(1) )
  hopt = 0.D0
  ndx  = max0(dims,nglob,1)
  nfix = 0

  call getnvarfix(iexp-1, nfix)

  allocate( GA0(nvar,n+dims) )
  allocate( D0(nvar) )
  allocate( XDIR(nvar,nvar) )
  allocate( PDIR(nglob,nglob) )
  allocate( SDIR(dims,nalg+dims) )
  allocate( SZR(nalg+dims) )
  allocate( TAU(nvar) )
  allocate( X(nvar) )
  allocate( U(nvar) )
  allocate( MST(dims,nvar) )
  allocate( MPT(nglob,nvar) )
  allocate( GA(nvar,ndir) )
  allocate( GADIR(n,ndir) )
  allocate( ZI(nmess) )
  allocate( DMDX(nmess,ndx) )
  allocate( DMDP(nmess,nglob) )
  allocate( GR(dimr,ndx) )
  allocate( R(dimr) )
  allocate( NODE(3) )

  call dinit2(nvar*(n+dims), 0.D0, GA0, 1)
  call dinit2(nvar, 0.D0, D0, 1)
  call dinit2(nvar*nvar, 0.D0, XDIR, 1)
  call dinit2(nglob*nglob, 0.D0, PDIR, 1)
  call dinit2(dims*(nalg+dims), 0.D0, SDIR, 1)
  call dinit2(nalg+dims, 0.D0, SZR, 1)
  call dinit2(nvar*ndir, 0.D0, GA, 1)
  call dinit2(n*ndir, 0.D0, GADIR, 1)
  call dinit2(nmess, 0.D0, ZI, 1)
  call dinit2(nmess*ndx, 0.D0, DMDX, 1)
  call dinit2(nmess*nglob, 0.D0, DMDP, 1)
  call dinit2(dimr*ndx, 0.D0, GR, 1)
  call dinit2(dimr, 0.D0, R, 1)

  cont  = 0.d0
  least = 0.d0

  do i = 1,dimr
    contr(i) = 0.d0
  enddo

  do i = 1,nalg
    SZR(i) = STATES(ndiff+k,1)
  enddo

  do i = 1,dims
    SZR(nalg+i) = SR(i)
  enddo

999 format(90(1x,1pe15.5))

! initialization; Computation of the block of the Jacobian t=t_0

  do i = 1,ndiff
    GA0(i,i) = 1.D0
  enddo

  call initg(iexp-1, STATES, PARAM, SZR, GA0(1,ndiff+1), nvar, D0)

  do i = 1,ndiff
    cont = cont + D0(i)*D0(i)
  enddo

  do i = 1,nvar
    XDIR(i,i) = 1.D0
  enddo

  do i = 1,nglob
    PDIR(i,i) = 1.D0
  enddo

  do i = 1,dims
    SDIR(i,nalg+i) = 1.D0
  enddo

  if ( nalg > 0 ) then

    call vplgfcn(iexp, 0, T(1), STATES, D0(ndiff+1), PARAM, rdum, idum, iflag)

    call x_vplgfcn(iexp, 0, T(1), STATES, D0(ndiff+1), PARAM, nalg, nvar, &
                   XDIR, nvar, GA0(ndiff+1,1), nvar, rdum, idum, iflag)

    call p_vplgfcn(iexp, 0, T(1), STATES, D0(ndiff+1), PARAM, nalg, nglob, &
                   PDIR, nglob, GA0(ndiff+1,nvar+1), nvar, rdum, idum, &
                   iflag)

! QR-decopmposition of GA0(1:nvar,1:nvar) to buildup the initial seed matrices

    call qrdec(GA0, nvar, TAU, iflag)

    if ( iflag < 0 ) then
      write(*,111) -iflag
      infor = -100
      return
    endif
111 format('Error, startup matrix is not definit! This might be an index problem.'/ &
           'Entry ', i3, ' has an illegal value.')

! loop to compute MR(-1)

    do j = 1,nvar

      s = 0.D0

      do k = j,nvar
        s = s + GA0(k,j)*D0(k)
      enddo

      do k = j,nvar
        D0(k) = D0(k) - GA0(k,j)*s
      enddo

    enddo

    call dcopy(nvar,D0,1,MR,1)

    do i = nvar,1,-1

      do j = i+1,nvar
        MR(i,1) = MR(i,1) - GA0(i,j)*MR(j,1)
      enddo

      MR(i,1) = MR(i,1)/TAU(i)
    enddo

    do j = 1,nvar
      MR(j,1) = - MR(j,1)
    enddo

! loop to compute MS(-1)

    do i = 1,dims

      do j = 1,nvar

        s = 0.D0

        do k = j,nvar
          s = s + GA0(k,j)*GA0(k,n+i)
        enddo

        do k = j,nvar
          GA0(k,n+i) = GA0(k,n+i) - GA0(k,j)*s
        enddo

      enddo

    enddo

    if ( dims > 0 ) call dcopy(nvar*dims,GA0(1,n+1),1,MS,1)

    do l = 1,dims

      do i = nvar,1,-1

        do j = i+1,nvar
          MS(i,l,1) = MS(i,l,1) - GA0(i,j)*MS(j,l,1)
        enddo

        MS(i,l,1) = MS(i,l,1)/TAU(i)

      enddo

    enddo

    do i = 1,nvar
      do j = 1,dims
        MS(i,j,1) = - MS(i,j,1)
      enddo
    enddo

! loop to compute MP(-1)

    do i = 1,nglob

      do j = 1,nvar

        s = 0.D0

        do k = j,nvar
          s = s + GA0(k,j)*GA0(k,nvar+i)
        enddo

        do k = j,nvar
          GA0(k,nvar+i) = GA0(k,nvar+i) - GA0(k,j)*s
        enddo

      enddo

    enddo

    call dcopy(nvar*nglob,GA0(1,nvar+1),1,MP,1)

    do l = 1,nglob

      do i = nvar,1,-1

        do j = i+1,nvar
          MP(i,l,1) = MP(i,l,1) - GA0(i,j)*MP(j,l,1)
        enddo

        MP(i,l,1) = MP(i,l,1)/TAU(i)

      enddo

    enddo

    do i = 1,nvar
      do j = 1,nglob
        MP(i,j,1) = - MP(i,j,1)
      enddo
    enddo

  else

    do i = 1,nvar

      MR(i,1) = - D0(i)

      do j = 1,dims
        MS(i,j,1) = - GA0(i,nvar+nglob+j)
      enddo

      do j = 1,nglob
        MP(i,j,1) = - GA0(i,nvar+j)
      enddo

    enddo

  endif

  call dinit2(nvar*(n+dims), 0.D0, GA0, 1)
  call dinit2(nvar, 0.D0, D0, 1)


! loop over all stopping points (shooting and measuring times)
! integration and computation of derivatives

  i  = 1 ! which shooting interval
  j  = 1 ! number of current interval
  ii = 1

  call getnodetype(iexp-1, 0, NODE)

  if ( nms == 1 ) then

    if ( i == 1 ) call dcopy(nvar,STATES(1,i),1,X,1)

    ld = -2

! update of the righthandside

    if ( NODE(1) == 1 ) then ! MSNODE?

      call x_vplgfcn(iexp, j-1, T(j), X, D0, PARAM, nalg, 1, MR(1,ii), 1, &
                     GA0, nvar, rdum, idum, iflag)

      do k = 1,nalg

        cont = cont + D0(k)*D0(k)

        RHS(dimr+k) = RHS(dimr+k) + GA0(k,1) + D0(k)

      enddo

      call dinit2(nvar*(n+dims), 0.D0, GA0, 1)
      call dinit2(nvar, 0.D0, D0, 1)

    endif

! update of part (S) of the Jacobian

    if ( NODE(1) == 1 ) then ! MSNODE?

      do k = 1,nvar
        do l = 1,dims
            MST(l,k) = MS(k,l,ii)
        enddo
      enddo

      call x_vplgfcn(iexp, j-1, T(j), X, D0, PARAM, nalg, dims, MST, dims, &
                     GA0, nvar, rdum, idum, iflag)

      do k = 1,nalg
        do l = 1,dims
          JACOBIAN(dimr+k,l) = JACOBIAN(dimr+k,l) + GA0(k,l)
        enddo
      enddo

      call dinit2(nvar*(n+dims), 0.D0, GA0, 1)
      call dinit2(nvar, 0.D0, D0,1)

    endif

! update of the part (P) of the Jacobian

    if ( NODE(1) == 1 ) then ! MSNODE?

      do k = 1,nvar
        do l = 1,nglob
          MPT(l,k) = MP(k,l,ii)
        enddo
      enddo

      call x_vplgfcn(iexp, j-1, T(j), X, D0, PARAM, nalg, nglob, MPT, &
                     nglob, GA0, nvar, rdum, idum, iflag)

      do k = 1,nalg
        do l = 1,nglob
          JACOBIAN(dimr+k,dims+l) = JACOBIAN(dimr+k,dims+l) + GA0(k,l)
        enddo
      enddo

      call dinit2(nvar*(n+dims), 0.D0, GA0, 1)
      call dinit2(nvar, 0.D0, D0,1)

      call p_vplgfcn(iexp, j-1, T(j), X, D0, PARAM, nalg, nglob, PDIR, &
                     nglob, GA0, nvar, rdum, idum, iflag)

      do k = 1,nalg
        do l = 1,nglob
          JACOBIAN(dimr+k,dims+l) = JACOBIAN(dimr+k,dims+l) + GA0(k,l)
        enddo
      enddo

      call dinit2(nvar*(n+dims), 0.D0, GA0, 1)
      call dinit2(nvar, 0.D0, D0,1)

    endif

  endif

  do while ( i < nms )

    if (  NODE(1) == 1 ) then

      if ( i == 1 ) call dcopy(nvar, STATES(1,i), 1, X, 1)

!       do k = 1,nalg
!         SZR(k) = STATES(ndiff+k,i)
!       enddo

      ld = -2
      ii = i

    else

      ld = -1
      ii = i + 1

    endif

! update of the righthandside

    if ( NODE(3) == 1 ) then ! RNODE?

!       call vplrfcn(iexp, j-1, ld, T(j), X, PARAM, SZR, R, dimr, rdum, idum, &
!                    iflag)

      call x_vplrfcn(iexp, j-1, ld, T(j), X, PARAM, SZR, R, dimr, 1, MR(1,ii), &
                     1, GR, dimr, rdum, idum, iflag)

      do k = 1,dimr

        contr(k) = contr(k) + R(k)

        RHS(k) = RHS(k) + GR(k,1) + R(k)

      enddo

      call dinit2(dimr*(ndx), 0.D0, GR, 1)
      call dinit2(dimr, 0.D0, R, 1)

    endif

    if ( NODE(1) == 1 ) then ! MSNODE?

      call x_vplgfcn(iexp, j-1, T(j), X, D0, PARAM, nalg, 1, MR(1,ii), 1, &
                     GA0, nvar, rdum, idum, iflag)

      do k = 1,nalg

        cont = cont + D0(k)*D0(k)

        RHS(dimr+k) = RHS(dimr+k) + GA0(k,1) + D0(k)

      enddo

      call dinit2(nvar*(n+dims), 0.D0, GA0, 1)
      call dinit2(nvar, 0.D0, D0, 1)

    endif

    if ( NODE(2) == 1 ) then ! MEASNODE?

      call vplmfcn(iexp, j-1, T(j), X, PARAM, ZI, nmess, rdum, idum, iflag)

      call z_vplmfcn(iexp, j-1, ld,  nmess, T(j), X, PARAM, MR(1,ii), 1, &
                     dummy, 0, DMDX, nmess, dummy, 1, rdum, idum, iflag)

      do k = 1,nmess

        least = least + ZI(k)*ZI(k)

        RHS(dimr+nalg+k) = RHS(dimr+nalg+k) + DMDX(k,1) + ZI(k)

      enddo

      call dinit2(nmess*ndx, 0.D0, DMDX, 1)
      call dinit2(nmess, 0.D0, ZI, 1)

    endif

! update of part (S) of the Jacobian

    if ( NODE(3) == 1 ) then ! RNODE?

      do k = 1,nvar
        do l = 1,dims
          MST(l,k) = MS(k,l,ii)
        enddo
      enddo

      call x_vplrfcn(iexp, j-1, ld, T(j), X, PARAM, SZR, R, dimr, dims, MST, &
                     dims, GR, dimr, rdum, idum, iflag)

      do k = 1,dimr
        do l = 1,dims
          JACOBIAN(k,l) = JACOBIAN(k,l) + GR(k,l)
        enddo
      enddo

      call dinit2(dimr, 0.D0, R, 1)
      call dinit2(dimr*ndx, 0.D0, GR, 1)

      call s_vplrfcn(iexp, j-1, ld, T(j), X, PARAM, SZR, R, dimr, dims, SDIR, &
                     dims, GR, dimr, rdum, idum, iflag)

      do k = 1,dimr
        do l = 1,dims
          JACOBIAN(k,l) = JACOBIAN(k,l) + GR(k,l)
        enddo
      enddo

      call dinit2(dimr, 0.D0, R, 1)
      call dinit2(dimr*ndx, 0.D0, GR, 1)

    endif

    if ( NODE(1) == 1 ) then ! MSNODE?

      do k = 1,nvar
        do l = 1,dims
            MST(l,k) = MS(k,l,ii)
        enddo
      enddo

      call x_vplgfcn(iexp, j-1, T(j), X, D0, PARAM, nalg, dims, MST, dims, &
                     GA0, nvar, rdum, idum, iflag)

      do k = 1,nalg
        do l = 1,dims
          JACOBIAN(dimr+k,l) = JACOBIAN(dimr+k,l) + GA0(k,l)
        enddo
      enddo

      call dinit2(nvar*(n+dims), 0.D0, GA0, 1)
      call dinit2(nvar, 0.D0, D0,1)

    endif

    if ( NODE(2) == 1 ) then ! MEASNODE?

      call z_vplmfcn(iexp, j-1, ld, nmess, T(j), X, PARAM, MS(1,1,ii), dims, &
                     dummy, 0, DMDX, nmess, dummy, 1, rdum, idum, iflag)

      do k = 1,nmess
        do l = 1,dims
          JACOBIAN(dimr+nalg+k,l) = JACOBIAN(dimr+nalg+k,l) + DMDX(k,l)
        enddo
      enddo

      call dinit2(nmess*ndx, 0.D0, DMDX, 1)

    endif

! update of the part (P) of the Jacobian

    if ( NODE(3) == 1 ) then ! RNODE?

      do k = 1,nvar
        do l = 1,nglob
          MPT(l,k) = MP(k,l,ii)
        enddo
      enddo

      call x_vplrfcn(iexp, j-1, ld, T(j), X, PARAM, SZR, R, dimr, nglob, MPT, &
                     nglob, GR, dimr, rdum, idum, iflag)

      do k = 1,dimr
        do l = 1,nglob
          JACOBIAN(k,dims+l) = JACOBIAN(k,dims+l) + GR(k,l)
        enddo
      enddo

      call dinit2(dimr, 0.D0, R, 1)
      call dinit2(dimr*ndx, 0.D0, GR, 1)

      call p_vplrfcn(iexp, j-1, ld, T(j), X, PARAM, SZR, R, dimr, nglob, PDIR, &
                     nglob, GR, dimr, rdum, idum, iflag)

      do k = 1,dimr
        do l = 1,nglob
          JACOBIAN(k,dims+l) = JACOBIAN(k,dims+l) + GR(k,l)
        enddo
      enddo

      call dinit2(dimr, 0.D0, R, 1)
      call dinit2(dimr*ndx, 0.D0, GR, 1)

    endif

    if ( NODE(1) == 1 ) then ! MSNODE?

      do k = 1,nvar
        do l = 1,nglob
          MPT(l,k) = MP(k,l,ii)
        enddo
      enddo

      call x_vplgfcn(iexp, j-1, T(j), X, D0, PARAM, nalg, nglob, MPT, &
                     nglob, GA0, nvar, rdum, idum, iflag)

      do k = 1,nalg
        do l = 1,nglob
          JACOBIAN(dimr+k,dims+l) = JACOBIAN(dimr+k,dims+l) + GA0(k,l)
        enddo
      enddo

      call dinit2(nvar*(n+dims), 0.D0, GA0, 1)
      call dinit2(nvar, 0.D0, D0,1)

      call p_vplgfcn(iexp, j-1, T(j), X, D0, PARAM, nalg, nglob, PDIR, &
                     nglob, GA0, nvar, rdum, idum, iflag)

      do k = 1,nalg
        do l = 1,nglob
          JACOBIAN(dimr+k,dims+l) = JACOBIAN(dimr+k,dims+l) + GA0(k,l)
        enddo
      enddo

      call dinit2(nvar*(n+dims), 0.D0, GA0, 1)
      call dinit2(nvar, 0.D0, D0,1)

    endif

    if ( NODE(2) == 1 ) then ! MEASNODE?

      call z_vplmfcn(iexp, j-1, ld, nmess, T(j), X, PARAM, MP(1,1,ii), nglob, &
                     PDIR, nglob, DMDX, nmess, DMDP, nmess, rdum, idum, iflag)

      do k = 1,nmess
        do l = 1,nglob
          JACOBIAN(dimr+nalg+k,dims+l) = JACOBIAN(dimr+nalg+k,dims+l) &
                                        + DMDX(k,l) + DMDP(k,l)
        enddo
      enddo

      call dinit2(nmess*ndx, 0.D0, DMDX, 1)
      call dinit2(nmess*nglob, 0.D0, DMDP, 1)

    endif

! integration until the next stopping point
! and computation of the new seed matrices

    h = dabs( MSTIMES(i+1) - MSTIMES(i) )

    do k = 1,dims

      do l = 1,nvar

        GADIR(l,k) = MS(l,k,ii)
        GA(l,k)    = MS(l,k,ii)

      enddo

    enddo

    do k = 1,nglob

      do l = 1,nvar

        GADIR(l,dims+k) = MP(l,k,ii)
        GA(l,dims+k)    = MP(l,k,ii)

      enddo

      do l = 1,nglob
        GADIR(nvar+l,dims+k) = PDIR(l,k)
      enddo

    enddo

    do l = 1,nvar

      GADIR(l,dims+nglob+1) = MR(l,ii)
      GA(l,dims+nglob+1)    = MR(l,ii)

    enddo

    call vpldaepar(iexp, j, nvar, ndiff, n, GADIR, n, ndir, T(j), T(j+1), X, &
                   PARAM, 1.D-6, hmax, h, hopt, 1, GA, nvar, rdum, idum, &
                   iflag)

    if ( iflag < 0 ) then
      write(*,*) 'No integretor-convergence.'
      infor = -iexp
      return
    endif

! copy the result to MS

    do k = 1,dims
      do l = 1,nvar
        MS(l,k,i+1) = GA(l,k)
      enddo
    enddo

! copy the result to MP

    do k = 1,nglob
      do l = 1,nvar
        MP(l,k,i+1) = GA(l,dims+k)
      enddo
    enddo

! copy the result to MR

    do l = 1,nvar
      MR(l,i+1) = GA(l,dims+nglob+1)
    enddo

    call dinit2(n*ndir, 0.D0, GADIR, 1)
    call dinit2(nvar*ndir, 0.D0, GA, 1)

    j = j + 1

    call getnodetype(iexp-1, j-1, NODE)

    if ( NODE(1) == 1 ) then

      i = i + 1

      if ( i == nms ) then

        call dcopy(nvar, X, 1, STATES(1,nms), 1)

        do k = 1,nvar
          IBEL(k,nms) = 1
        enddo

      endif

      call dcopy(nvar, X, 1, U, 1)

      do k = 1,nvar

        if ( IBEL(k,i) == 1 ) then

!           call dcopy(nvar, STATES(1,i), 1, X, 1)
          X(k) = STATES(k,i)

        else

!           call dcopy(nvar, X, 1, STATES(1,i), 1)
          STATES(k,i) = X(k)

          IBEL(k,i) = 1

        endif

      enddo

! evaluation of the continuity conditions: X - STATES(i)

      do k = nfix+1,ndiff

        cont = cont + ( U(k) - STATES(k,i) )*( U(k) - STATES(k,i) )

        MR(k,i) = MR(k,i) + ( U(k) - STATES(k,i) )

      enddo

      do k = 1,nfix

        MR(k,i) = 0.d0

        do l = 1,nglob
          MP(k,l,i) = 0.d0
        enddo

        do l = 1,dims
          MS(k,l,i) = 0.d0
        enddo

      enddo

!       call vplgfcn(iexp, nt, T(nt), STATES(1,i), D0, PARAM, rdum, idum, iflag)
! 
!       do k = 1,nalg
!         MR(ndiff+k,i) = MR(ndiff+k,i) + D0(k)
!       enddo
! 
!       call dinit2(nvar, 0.D0, D0, 1)

    endif

  enddo

! The shooting node at the endpoint of the time interval is removed from
! the system!
! But the constraints at the last time point still have to be evaluated

!   do k = 1,nalg
!     SZR(k) = STATES(ndiff+k,nms)
!   enddo

  if ( NODE(3) == 1 ) then

    call x_vplrfcn(iexp, nt-1, -2, T(nt), STATES(1,nms), PARAM, SZR, R, dimr, &
                   1, MR(1,nms), 1, GR, dimr, rdum, idum, iflag)

    do k = 1,dimr

      contr(k) = contr(k) + R(k)

      RHS(k) = RHS(k) + R(k) + GR(k,1)

    enddo

    call dinit2(dimr, 0.D0, R, 1)
    call dinit2(dimr*ndx, 0.D0, GR, 1)

    do k = 1,nvar
      do l = 1,dims
        MST(l,k) = MS(k,l,nms)
      enddo
    enddo

    call x_vplrfcn(iexp, nt-1, -2, T(nt), STATES(1,nms), PARAM, SZR, R, dimr, &
                   dims, MST, dims, GR, dimr, rdum, idum, iflag)

    do k = 1,dimr
      do l = 1,dims
        JACOBIAN(k,l) = JACOBIAN(k,l) + GR(k,l)
      enddo
    enddo

    call dinit2(dimr, 0.D0, R, 1)
    call dinit2(dimr*ndx, 0.D0, GR, 1)

    call s_vplrfcn(iexp, nt-1, -2, T(nt), STATES(1,nms), PARAM, SZR, R, dimr, &
                   dims, SDIR, dims, GR, dimr, rdum, idum, iflag)

    do k = 1,dimr
      do l = 1,dims
        JACOBIAN(k,l) = JACOBIAN(k,l) + GR(k,l)
      enddo
    enddo

    call dinit2(dimr, 0.D0, R, 1)
    call dinit2(dimr*ndx, 0.D0, GR, 1)

    do k = 1,nvar
      do l = 1,nglob
        MPT(l,k) = MP(k,l,nms)
      enddo
    enddo

    call x_vplrfcn(iexp, nt-1, -2, T(nt), STATES(1,nms), PARAM, SZR, R, dimr, &
                   nglob, MPT, nglob, GR, dimr, rdum, idum, iflag)

    do k = 1,dimr
      do l = 1,nglob
        JACOBIAN(k,dims+l) = JACOBIAN(k,dims+l) + GR(k,l)
      enddo
    enddo

    call dinit2(dimr, 0.D0, R, 1)
    call dinit2(dimr*ndx, 0.D0, GR, 1)

    call p_vplrfcn(iexp, nt-1, -2, T(nt), STATES(1,nms), PARAM, SZR, R, dimr, &
                   nglob, PDIR, nglob, GR, dimr, rdum, idum, iflag)

    do k = 1,dimr
      do l = 1,nglob
        JACOBIAN(k,dims+l) = JACOBIAN(k,dims+l) + GR(k,l)
      enddo
    enddo

    call dinit2(dimr, 0.D0, R, 1)
    call dinit2(dimr*ndx, 0.D0, GR, 1)

  endif

  if ( NODE(2) == 1 ) then

    call vplmfcn(iexp, nt-1, T(nt), STATES(1,nms), PARAM, ZI, nmess, rdum, &
                 idum, iflag)

    call z_vplmfcn(iexp, nt-1, -2, nmess, T(nt), STATES(1,nms), PARAM, &
                   MR(1,nms), 1, dummy, 0, DMDX, nmess, dummy, 1, rdum, idum, &
                   iflag)

    do k = 1,nmess

      least = least + ZI(k)*ZI(k)

      RHS(dimr+nalg+k) = RHS(dimr+nalg+k) + ZI(k) + DMDX(k,1)

    enddo

    call dinit2(nmess, 0.D0, ZI, 1)
    call dinit2(nmess*(ndx), 0.D0, DMDX, 1)

    if ( dims > 0 ) then

      call z_vplmfcn(iexp, nt-1, -2, nmess, T(nt), STATES(1,nms), PARAM, &
                     MS(1,1,nms), dims, dummy, 0, DMDX, nmess, dummy, 1, &
                     rdum, idum, iflag)

      do k = 1,nmess
         do l = 1,dims
           JACOBIAN(dimr+nalg+k,l) = JACOBIAN(dimr+nalg+k,l) + DMDX(k,l)
         enddo
      enddo

      call dinit2(nmess*(ndx), 0.D0, DMDX, 1)

    endif

    if ( nglob > 0 ) then

      call z_vplmfcn(iexp, nt-1, -2, nmess, T(nt), STATES(1,nms), PARAM, &
                     MP(1,1,nms), nglob, PDIR, nglob, DMDX, nmess, DMDP, &
                     nmess, rdum, idum, iflag)

      do k = 1,nmess
        do l = 1,nglob
          JACOBIAN(dimr+nalg+k,dims+l) = JACOBIAN(dimr+nalg+k,dims+l) &
                                               + DMDX(k,l) + DMDP(k,l)
        enddo
      enddo

      call dinit2(nmess*ndx, 0.D0, DMDX, 1)
      call dinit2(nmess*nglob, 0.D0, DMDP, 1)

    endif

  endif

  deallocate(NODE)
  deallocate(R)
  deallocate(GR)
  deallocate(DMDP)
  deallocate(DMDX)
  deallocate(ZI)
  deallocate(GADIR)
  deallocate(GA)
  deallocate(MPT)
  deallocate(MST)
  deallocate(U)
  deallocate(X)
  deallocate(TAU)
  deallocate(SZR)
  deallocate(SDIR)
  deallocate(PDIR)
  deallocate(XDIR)
  deallocate(D0)
  deallocate(GA0)

end ! subroutine intdiff

!--------------------------------------------------------------------------------!

! QR-Decomposition

!> \routine qrdec
!> \brief QR-decomposition of the matrix A
!> \author Robert Kircheis

subroutine qrdec(A,nvar,TAU,iflag)

  implicit none

  integer :: nvar, &
             iflag
  real*8 :: A(nvar,nvar), &
            TAU(nvar)

  integer :: i,j,k
  real*8 :: s, &
            fak

  iflag = 0

  do j = 1,nvar

    s = 0.D0

    do i = j,nvar
      s = s + A(i,j)*A(i,j)
    enddo

    s = dsqrt(s)

    if ( A(j,j) > 0.0 ) then
      TAU(j) = -s
    else
      TAU(j) = s
    endif

    fak = dsqrt( s*( s + DABS( A(j,j) ) ) )

    if ( fak < 1.D-8 ) then
      iflag = -j
      write(*,123) 
123 format(''/'Columns of the initial seed matrix &
           seem to be linearly dependent!'/'')
      exit
    endif

    A(j,j) = A(j,j) - TAU(j)

    do k = j,nvar
      A(k,j) = A(k,j)/fak
    enddo

    do i = j+1,nvar

      s = 0.D0

      do k = j,nvar
        s = s + A(k,j)*A(k,i)
      enddo

      do k = j,nvar
        A(k,i) = A(k,i) - A(k,j)*s
      enddo

    enddo

  enddo

end ! subroutine qrdec
