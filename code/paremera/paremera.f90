!-------------------------------------------------------------------------------+
! Parameter Estimation Tool for Multiple Experiments using a Reduced Approach   |
! Implemented by R. Kircheis 10.2011                                            |
!-------------------------------------------------------------------------------+
!> @subroutine paremera
!! @file paremera.f95
!! @brief parameter estimation tool using a reduced approach
!> @author Robert Kircheis

subroutine paremera(X1, nexp, NMS, NSP, DIMR, DIMRLOC, multidimr, DIMS, NDIFF, &
                    NVAR, nglob, NMESS, NBED, dimx, eps, itmax, kout, cond0, &
                    index1, fashort, farel, famax, fastart, PUPBND, PLOBND, bcheck, &
                    sflag, lmflag, pr)

  implicit none

  integer :: nexp, &  !< number of experiments
             NMS(nexp), & !< number of shooting nodes for each experiment
             NSP(nexp), & !< number of 'stoping' points for each experiment
             DIMR(nexp), & !< number of rfcns ( local + multi ) for each experiment
             DIMRLOC(nexp), & !< number of local rfcn functions for each experiment
             multidimr, & !< number of multiexperiment rfcn functions
             DIMS(nexp), & !< length of the s-vector (-nalg)  for each experiment
             NDIFF(nexp), & !< number of differential states for each experiment
             NVAR(nexp), & !< number of states for each experiment
             nglob, & !< number of parameters
             NMESS(nexp), & !< number of measurements for each experiment
             NBED(nexp), & !< dimr + nmess + nglob
             dimx, & !< number of all variables
             itmax, & !< maximum number of iterations
             kout, & !< error index1
             bcheck, & !< flag if the parameter bounds check is active or not
             sflag, & !< startflag for parameter estimation
             lmflag, & !< flag for levenberg marquardt algorithm
             pr !< printlevel
  real*8 :: X1(dimx), & !< array of all variables (in and output variable)
            eps, & !< input for the termination accuracy
            index1, & !< accuracy for the index-one condition
            cond0, & !< accuracy measure for the condition of the Jacobian
            fashort, & !< scalar for the step length
            farel, & !<  minimum step length
            famax, & !<  maximum step length
            fastart, & !< first step length (will be divided by the norm of the increment)
            PUPBND(nglob), & !< upper bound for the parameters
            PLOBND(nglob), & !< lower bound for the parameters
            lvdelta, & !< Step bound in Levenberg Marquardt
            lvsigma, & !< for Levenberg Marquardt, error in norm(D*p)
            lvlambda !< The Levenberg Marquardt Parameter

  logical :: convergence, & !< convergence yes or no
             vagt30 !< variable a greater than 10-30
  integer :: i, j, k, l, & !< numerators
             iter, & !< iteration counter
             ifct, & !< counter for function calls
             np, & !< for the call of getnp
             iposp, & !< position of the parameters in the array X1 (from getnp)
             nvtot, & !< length of the array X1
             dimrtot, & !< sum over all experiments of dimrloc
             maxnd, & !< maximum number of differential states per experiment
             maxnv, & !< maximum number of states per experiment
             maxms, & !< maximum number of shooting nodes per experiment
             maxsp, & !< maximum number of 'stopping' points per experiment
             maxnt, & !< maximum number of time steps per experiment
             maxnb, & !< maximum number of conditions; max(nbed)
             maxsr, & !< maximum number of s-variables per experiment
             maxsloc, & !< maximum number of local rfcns per experiment
             rank, & !< rank of the Jacobian
             nbedtotal, & !< sum(nbed)
             intflag, & !< flag, if integration was succesfull
             expsum ! summe über NMESS
  real*8 :: fa, & !< relaxation factor for current iteration
            fa0, & !< relaxation factor for previous iteration
            fpro, & !< 1.D0/( omega*dsqrt( dxrnorm ) ) for RMT
            dxnorm, & !< norm of the increment
            dxrnorm, & !< norm of the reduced increment
            dpnorm, & !< norm of the increment for the parameters
            dxnormold, & !< dxnorm of the previous iteration
            dxhilfnorm, & !< norm of the intermediate step for the RMT
            omega, & !< omega(fa); see Bock2000
            kappa, & !< kappa estimation
            sigma, & !< fa0*fpro
            conv, & !< convergence rate
            epcon, & !< breakup accuracy (eps*eps*nvtot)
            equality, & !< continuity gap + consistence violation
            leastsquares, & !< residual
            eqold, & !< equality of the previous iteration
            lsold, & !< leastsquares of the previous iteration
            A,& !< auxillary variable
            rho,& !< Measures agreement between linear model and nonlinear function
            normf, & !< auxiliary variable
            normfplus,& !< auxiliary variable
            normJp,&
            normDp,&
            lvgamma ,&
            lvmu, &
            normDx ,&
            XTOL, &
            FTOL

  real*8, allocatable :: MS(:,:,:,:), & !< seed matrix for s-directions;
                                        !! dim(maxnv,maxsr,maxms,nexp)
                         MP(:,:,:,:), & !< seed matrix for p-directions;
                                        !! dim(maxnv,nglob,maxms,nexp)
                         MR(:,:,:), & !< seed matrix for r-direction;
                                      !! dim(maxnv,maxms,nexp)
                         MShilf(:,:,:,:), & !< see MS
                                            !! dim(maxnv,maxsr,maxms,nexp)
                         MPhilf(:,:,:,:), & !< see MP;
                                            !! dim(maxnv,nglob,maxms,nexp)
                         MRhilf(:,:,:), & !< see MR;
                                          !! dim(maxnv,maxms,nexp)
                         JACOBIAN(:,:,:), & !< Jacobian;
                                            !! dim(maxnb,maxsr+nglob,nexp)
                         JACOBIANmatrix(:,:),& !< store Jacobian in 2-d Matrix dim(expsum, maxsr+nglob)
                         RHS(:,:), & !< right hand side; dim(maxnb,nexp)
                         RHSvector(:), & !< right hand side as one vector; dim(expsum)
                         GLOBALMATRIX(:,:), & !< auxiliary matrix for the computation
                                              !! of delta_p; dim(nbedtotal,multidimr+nglob)
                         JACOBIANhilf(:,:,:), & !< auxiliary Jacobian;
                                                !! dim(maxnb,maxsr+nglob,nexp)
                         RHShilf(:,:), & !< auxiliary right hand side;
                                         !! dim(maxnb,nexp)
                         Xhilf(:), & !< auxiliary X-array; dim(nvtot
                         SCAL(:), & !< scaling vector; dim(nglob)
                         DX(:), & !< increments; dim(nvtot)
                         DXR(:), & !< reduced increments; dim(dimrtot+nglob)
                         DXhilf(:), & !< auxiliary vector of increments; dim(nvtot)
                         DXRhilf(:), & !< auxiliary vector of reducedincrements;
                                       !! dim(dimrtot+nglob)
                         X0(:), & !< vector of variables of the previous iteration;
                                  !! dim(nvtot)
                         PARAM(:), & !< vector of parameters; dim(nglob)
                         EXPWISEEQ(:), & !< equality per experiment; dim(exp)
                         EXPWISELS(:), & !< leastsquares per experiment; dim(exp)
                         TAU(:), & !< vector for QR-decomposition with Householder;
                                   !! dim(nglob)
                         PSEUDOINVERSE(:,:), & !< Pseudoinverse of Deuflhard's algorithm;
                                               !! dim(nglob,nglob)
                         LMDIAGONAL(:)   ! Scaling Matrix in Levenberg Marquardt algorithm dim(nglob)
  integer, allocatable :: IBEL(:,:,:), & !< flag if shooting nodes are initialized;
                                      !! dim(maxnv,maxms,nexp)
                          IPL(:,:), & !< pivot matrix (total); dim(maxsr,nexp)
                          IPR(:), & !< pivot matrix (row); dim(multidimr)
                          IPC(:), & !<  pivot matrix (columnwise); dim(nglob)
                          iposx(:), & !< pointer for the states of each experiment
                          ipossr(:) !< pointer for the rfcn of each experiment


! 1. initialization and allocation of the variables

  vagt30(A) = dabs(A) >= 1.d-30

  convergence  = .false.

  fa         = fastart
  fa0        = 0.d0
  fpro       = 0.d0
  ifct       = 0
  dxnorm     = 0.d0
  dxrnorm    = 0.d0
  dxnormold  = 0.d0
  dxhilfnorm = 0.d0
  omega      = 0.d0
  kappa      = 0.d0
  eqold      = 0.d0
  lsold      = 0.d0

  maxnd   = 1
  maxnv   = 1
  maxms   = 1
  maxsp   = 1
  maxnt   = 1
  maxnb   = 1
  maxsr   = 1
  maxsloc = 1

  lvdelta = 10000 !initialize delta_0
  lvsigma = 0.1
  lvlambda = 1 !lambda_0
  XTOL = 1.d-8  
  FTOl = 1.d-8

  expsum = sum(NMESS)

!values used in first iteration
  rho = 1
  normf=1
  normfplus = 0.5
  normJp = 1
  normDp = 1
  lvgamma = 1
!until here


  do i = 1,nexp

    maxnd   = max0( maxnd, NDIFF(i) )
    maxnv   = max0( maxnv, NVAR(i) )
    maxms   = max0( maxms, NMS(i) )
    maxsp   = max0( maxsp, NSP(i) )
    maxnt   = max0( maxnt, NMS(i)+NSP(i) )
    maxnb   = max0( maxnb, NBED(i) )
    maxsr   = max0( maxsr, DIMS(i) )
    maxsloc = max0( maxsloc, DIMRLOC(i) )

  enddo

  nbedtotal = 0
  nvtot     = nglob
  dimrtot   = multidimr

  do i = 1,nexp

    nbedtotal = nbedtotal + NBED(i) - DIMR(i)
    nvtot     = nvtot + NMS(i)*NVAR(i) + DIMS(i)
    dimrtot   = dimrtot + DIMRLOC(i)

  enddo

  nbedtotal = nbedtotal + nglob + multidimr
  epcon = eps*eps*nvtot

  allocate( MS(maxnv,maxsr,maxms,nexp) )
  allocate( MP(maxnv,nglob,maxms,nexp) )
  allocate( MR(maxnv,maxms,nexp) )
  allocate( MShilf(maxnv,maxsr,maxms,nexp) )
  allocate( MPhilf(maxnv,nglob,maxms,nexp) )
  allocate( MRhilf(maxnv,maxms,nexp) )
  allocate( JACOBIAN(maxnb,maxsr+nglob,nexp) )
  allocate( RHS(maxnb,nexp) )
  allocate( GLOBALMATRIX(nbedtotal,multidimr+nglob) )
  allocate( JACOBIANhilf(maxnb,maxsr+nglob,nexp) )
  allocate( RHShilf(maxnb,nexp) )
  allocate( Xhilf(nvtot) )
  allocate( SCAL(nglob) )
  allocate( DX(nvtot) )
  allocate( DXR(dimrtot+nglob) )
  allocate( DXhilf(nvtot) )
  allocate( DXRhilf(dimrtot+nglob) )
  allocate( X0(nvtot) )
  allocate( PARAM(nglob) )
  allocate( EXPWISEEQ(nexp) )
  allocate( EXPWISELS(nexp) )
  allocate( IPL(maxsr,nexp) )
  allocate( IPR(multidimr) )
  allocate( IPC(nglob) )
  allocate( TAU(nglob) )
  allocate( PSEUDOINVERSE(nglob,nglob) )
  allocate( IBEL(maxnv,maxms,nexp) )
  allocate( iposx(nexp) )
  allocate( ipossr(nexp) )

  if ( lmflag /= 0 ) then
    allocate( LMDIAGONAL(nglob))
    allocate( RHSvector(expsum))
    allocate( JACOBIANmatrix(expsum, nglob))

    !make LMDIAGONAL the Identity
    do j=1,nglob 
      LMDIAGONAL(j) = 1
    end do
  endif

  call dinit2(nvtot, 0.d0, DX, 1)
  call dinit2(dimrtot+nglob, 0.d0, DXR, 1)
  call dinit2(nvtot, 0.d0, X0, 1)
  call dinit2(nvtot, 0.d0, Xhilf, 1)
  call dinit2(nglob, 0.d0, PARAM, 1)
  call dinit2(nvtot, 0.d0, DXhilf, 1)
  call dinit2(dimrtot+nglob, 0.d0, DXRhilf, 1)

  do k = 1,nexp
    do i = 1,maxnv
      do j = 1,maxms
        IBEL(i,j,k) = 0
      enddo
    enddo
  enddo

  iposx(1) = 1
  do i = 1,nexp-1
    iposx(i+1) = iposx(i) + NMS(i)*NVAR(i)
  enddo

! 2. Starting the generalized Gauss-Newton

  open(unit=4,file="parameter.txt", action="write", status='replace')
  close(4)

  write(*,123)
123 format(5x,47('*') / 5x,'*',5x, 'Parameter estimation with PAREMERA'&
           ,6x,'*' / 5x,47('*') / )

  call getnp( np, iposp )

  ipossr(1) = iposp + nglob
  do i = 1,nexp-1
    ipossr(i+1) = ipossr(i) + DIMS(i)
  enddo

  do iter = 1, itmax

    do i = 1,nglob
      PARAM(i) = X1(iposp+i-1)
    enddo

    open(unit=4,file="parameter.txt", action="write", position='append')

    write(4,124) iter, (PARAM(i),i=1,nglob)
124 format(i3, 20(1pe15.5))

    close(4)

    equality        = 0.0d0
    leastsquares    = 0.0d0

    if ( lmflag == 0 .or. iter == 1 ) then

      call metse2(iter, nexp, X1, iposx, ipossr, JACOBIAN, RHS, MS, MP, MR, &
                EXPWISEEQ, EXPWISELS, NBED, nglob, multidimr, maxnd, maxnv, &
                maxnt, maxms, maxsp, maxnb, maxsr, nvtot, 0, IBEL, sflag, &
                pr, kout)

      sflag = 0

      ifct = ifct + 1

    endif

    if ( ( kout < 0 ) .and. ( iter == 1 ) ) then
      write(*,125)
125 format(/ 'Integration failed at first iteration.' / &
             'This is definitely a problem in the model.' / '')
      return
    endif

    if ( itmax == 1 ) return

    if ( lmflag == 0 ) then
      call scaljacobian(JACOBIAN, SCAL, NBED, DIMS, nglob, nexp, maxnb, maxsr, &
                        cond0)
    else
      do i = 1,nglob
        SCAL(i) = 1.d0
      enddo
    endif

    do i = 1,nexp

      equality     = equality + EXPWISEEQ(i)
      leastsquares = leastsquares + EXPWISELS(i)

    enddo

    if ( lmflag == 1 ) then
      ! 2.5 update Delta und D for Levenberg Marquardt

      ! write RHS in one vector
      j=0
      do i=1,nexp
        do k = 1,NMESS(i)
          RHSvector(j+k) = RHS(DIMR(i) + NVAR(i) - NDIFF(i) + k, i)
        enddo
        j = j + NMESS(i)
      enddo

      !write Jacobian in one Matrix
      l = 0
      do i = 1,nexp
        do j = 1,NMESS(i)
          do k = 1,(nglob)
              JACOBIANmatrix(l + j,k) = JACOBIAN(j,DIMS(i)+k,i)
          enddo
        enddo
        l = l + NMESS(i)
      enddo
     
      !update Delta part
      if (rho <= 0.25) then
          if (normfplus < normf) then
              lvmu = 0.5
          elseif (normfplus > 10 * normf ) then
              lvmu = 0.1
          else
              lvgamma = - ( (normJp/normf)**2 + (sqrt(lvlambda) * normDp /normf)**2)
              lvmu = 0.5* lvgamma / ( lvgamma + 0.5 * (1 - (normfplus/normf)**2))
          endif
          lvdelta = lvdelta * lvmu
      elseif ( 0.25 < rho .and. rho < 0.75 .and. lvlambda == 0 .or. rho >=0.75 .and. iter>1 ) then 
          lvdelta = 2* normDp
      endif  
      
      normf = sqrt(leastsquares) 

      !update D part
      do i=1,nglob    
        LMDIAGONAL(i) = max(LMDIAGONAL(i), sqrt(sum(JACOBIANmatrix(1:,i)**2)))
      enddo
    endif


! 3. Decompose the Jacobian

    rank = nglob
    intflag = -1

    do while ( intflag < 0 )

!       DJACOBIAN(1:maxnb,1:(maxsr+nglob),1:nexp) = &
!                JACOBIAN(1:maxnb,1:(maxsr+nglob),1:nexp)
!       DRHS(1:maxnb,1:nexp) = RHS(1:maxnb,1:nexp)

      call decjacobian(JACOBIAN, maxnb, GLOBALMATRIX, nbedtotal, PARAM, SCAL, &
                       PSEUDOINVERSE, TAU, rank, DIMR, DIMRLOC, DIMS, &
                       multidimr, maxsr, maxsloc, nglob, NBED, nexp, cond0, &
                       IPL, IPR, IPC, pr, kout)

      if ( kout < 0 ) then
        write(*,*) 'Termination of PAREMERA'
        return
      endif


! 4. Compute the increments for s_r and p

      call compredinc(JACOBIAN, maxnb, GLOBALMATRIX, nbedtotal, PARAM, maxsr, &
                      maxsloc, PSEUDOINVERSE, rank, TAU, SCAL, RHS, DXR, &
                      dimrtot, DIMR, DIMRLOC, DIMS, multidimr, nglob, NBED, &
                      nexp, IPL, IPR, IPC, LMDIAGONAL, lvdelta, lvsigma, &
                      lvlambda, RHSvector, JACOBIANmatrix, expsum, lmflag, kout)

      dxrnorm = 0.d0
      dpnorm  = 0.d0

      do j = 1,nglob+dimrtot
        if( vagt30( DXR(j) ) ) dxrnorm = dxrnorm + DXR(j)*DXR(j)
      enddo

      do j = 1,nglob
        if( vagt30( DXR(j) ) ) dpnorm = dpnorm + DXR(j)*DXR(j)
      enddo


! 5 Compute the increments for s_y and s_z

      call compinc(DX, DXR, MS, MP, MR, NMS, NVAR, NDIFF, DIMS, multidimr, &
                   nglob, NBED, nexp, nvtot, dimrtot, maxnv, maxsr, maxms, &
                   kout)

      dxnormold = dxnorm
      dxnorm = 0.0d0

      do J=1,nvtot
        if( vagt30(DX(j)) ) dxnorm = dxnorm + DX(j)*DX(j)
      enddo

      if ( iter == 1 ) then

        fa = fastart/dsqrt( dxrnorm )
        fa = dmin1( fa, fastart )

      else if ( lmflag == 0 ) then ! in Levenberg-Marquardt skip RMT

! 6. RMT; just one iteration
!         evaluation of omega to costly

        fa0 = fa

        do j = 1,nvtot
          Xhilf(j) = X1(j) + fa*DX(j)
        enddo

        if ( bcheck > 0 ) then
          call boundcheck( Xhilf, X1, DX, nvtot, fa, fashort, farel, PUPBND, &
                          PLOBND, nglob, kout )
        endif

        call metse2(iter, nexp, Xhilf, iposx, ipossr, JACOBIANhilf, RHShilf,&
                    MShilf, MPhilf, MRhilf, EXPWISEEQ, EXPWISELS, NBED, &
                    nglob, multidimr, maxnd, maxnv, maxnt, maxms, maxsp, &
                    maxnb, maxsr, nvtot, 1, IBEL, sflag, pr, kout)

        if ( kout < 0 ) then

          write(*,126)
126 format(/ 'Integration for stepsize control failed.' &
           / 'We try to go on with minimal stepsize.' )

          fa = farel

          kout = 0

        else

          ifct = ifct + 1

          call compredinc(JACOBIAN, maxnb, GLOBALMATRIX, nbedtotal, PARAM, &
                          maxsr, maxsloc, PSEUDOINVERSE, rank, TAU, SCAL, &
                          RHShilf, DXRhilf, dimrtot, DIMR, DIMRLOC, DIMS, &
                          multidimr, nglob, NBED, nexp, IPL, IPR, IPC, &
                          LMDIAGONAL, lvdelta, lvsigma, lvlambda, RHSvector, &
                          JACOBIANmatrix, expsum, lmflag, kout)

!           call COMPINC(DXhilf, DXRhilf, MShilf, MPhilf, MRhilf, NMS, NVAR, &
!                        NDIFF, DIMS, multidimr, nglob, NBED, nexp, nvtot, &
!                        dimrtot, maxnv, maxsr, maxms, kout)

          dxhilfnorm = 0.d0

          do j = 1,nglob+dimrtot

            DXRhilf(j) = DXRhilf(j) - (1.d0-fa)*DXR(j)

            if( vagt30( DXRhilf(j) ) ) dxhilfnorm = dxhilfnorm &
                                                   + DXRhilf(j)*DXRhilf(j)

          enddo

! Use only the reduced increment to compute the curvature information
!           do j=1,nvtot
!             DXhilf(j) = DXhilf(j) - (1.d0-fa)*DX(j)
!             if( vagt30(DX(j)) ) dxhilfnorm = dxhilfnorm + DXhilf(j)*DXhilf(j)
!           enddo

          omega = 2*dsqrt( dxhilfnorm )/( fa0*fa0*dxrnorm )
          fpro  = 1.d0/( omega*dsqrt( dxrnorm ) + 1.d-10)

          sigma = fa0/fpro
          conv  = dsqrt( dxnorm/dxnormold )

          fa = dmax1( farel, dmin1( famax, fpro ) )
!           fa = dmax1(farel,fa)

        endif

      endif

      ! Enforce fixed steplength

      if ( farel == famax ) then
         fa = famax
      endif

! 7. Test if problem is integrable at the new iterate

      do j = 1,nvtot
        Xhilf(j) = X1(j) + fashort*fa*DX(j)
      enddo

      if ( bcheck > 0 ) then
        call boundcheck( Xhilf, X1, DX, nvtot, fa, fashort, farel, PUPBND, &
                         PLOBND, nglob, kout )
      endif

      if ( kout < 0 ) then

        write(*,135)

        return
      endif
135   format( 5x, 'Start iteration from another starting point.' / )

!      call inttest( nexp, Xhilf, iposx, nvtot, pr, intflag )
     intflag = 0


! 8. If integration fails, try a shortened steplength

      do while ( intflag < 0 )

        fa = 0.5d0*fa

        if ( fa < farel ) exit

        do j=1,nvtot
          Xhilf(j) = X1(j) + fa*DX(j)
        enddo

        call inttest( nexp, Xhilf, iposx, nvtot, pr, intflag )

      enddo


! 9. If this fails to, reduce rank and compute a new increment.

      if ( intflag < 0 ) then

        write(*,132)
132 format(/ 'Cant find feasible iterate with the current increment.' / &
           / 'The rank will be reduced to compute a new one.')

        rank = rank - 1
        fa = fa0

      endif

      if ( rank == 0 ) then
        write(*,127)
127 format(/ 'Cant find feasible iterate for the integration.' &
             / 'Abbortion of PAREMERA' / '')
        return
      endif

    enddo

! Update the iterate

!     if ( rank < nglob ) fa = farel

    if ( lmflag == 0 ) then  !skip RMT in LM
        do j=1,nvtot

          X0(j) = X1(j)
          X1(j) = X1(j) + fashort*fa*DX(j)

        enddo
    else
    ! Levenbergmarquardt
        do j=1,nvtot

          X0(j) = X1(j)
          X1(j) = X1(j) +DX(j)

        enddo

        !compute rho

        call metse2(iter+1, nexp, X1, iposx, ipossr, JACOBIANhilf, RHShilf, MS, MP, MR,&
                    EXPWISEEQ, EXPWISELS, NBED, nglob, multidimr, maxnd, maxnv,&
                    maxnt, maxms, maxsp, maxnb, maxsr, nvtot, 0, IBEL, sflag, &
                    pr, kout)
        ifct = ifct + 1

        !write RHS in one vector
        j=0
        do i=1,nexp
          do k = 1,NMESS(i)
            RHSvector(j+k) = RHShilf(DIMR(i) + NVAR(i) - NDIFF(i) + k, i)
          enddo
          j = j + NMESS(i)
        enddo

        !write Jacobian in one Matrix
        l = 0
        do i = 1,nexp
          do j = 1,NMESS(i)
            do k = 1,nglob
                JACOBIANmatrix(l + j,k) = JACOBIANhilf(j,DIMS(i)+k,i)
            enddo
          enddo
          l = l + NMESS(i)
        enddo

        !compute normfplus, normJp, normDp
 
        normfplus = sqrt(sum(RHSvector**2))
        normDp = sqrt(sum((LMDIAGONAL * DXR(dimrtot + 1:))**2))  
        normJp = sqrt(sum(matmul(JACOBIANmatrix , DXR(dimrtot +1:))**2))

        if (normfplus > normf) then
            rho = 0
        else
            rho = (1 - (normfplus/normf)**2) / ( (normJp/normf)**2 + 2 * (sqrt(lvlambda) * normDp /normf)**2)
        endif

        if (rho < 0.0001) then
          X1 = X0 !reject step 
        else 
          JACOBIAN = JACOBIANhilf
          RHS = RHShilf
        endif

    endif

    call vpl_save_x( nvtot, X1 )

    if ( pr == 0 ) then

      if ( iter == 1 ) then
        write(*,*)
        write(*,*) 'ITER |  EQUALITY   | LEAST SQUARES |   &
                    INCREM    |  RELAX | RANK '
      endif

      write(*,133) iter, equality, leastsquares, dxnorm, fa, rank
133 format(i4,3x,1pe12.5,3x,1pe12.5,3x,1pe12.5,1pe10.2,2x,i3)

    else

      write(*,128) iter,ifct,fa,equality,leastsquares,dxnorm, dpnorm
128 format(/5x,62('-')/5x,'|',2x,'ITERATION',i5,' FCALLS', &
            i5,' RELAXATIONS-FACTOR ',1pe10.3,2x,'|' / &
            5x,'|',60x,'|' / &
            5x,'|',2x,'LEVEL FUNCTIONS',43x,'|' / &
            5x,'|',4x,'EQUALITY',2x,'|LEAST SQUARES|',4x,'INCREM',3x'|', &
            3x,'INCREMP',7x,'|' / &
            5x,'| ',3(1pE12.5,1X,'|'),1pE12.5,5X,'|' / 5x,62('-') / '')

    endif

    if ( ( iter > 1 ).and.( pr > 0 ) ) then
      write(*,129) omega, conv, fpro, sigma, fa, fa0
    endif
129 format(2x,'APRIORI ESTIMATION FROM OLD STEP DATA' / &
           1x,' OMEGA',7x,'CONV',8x,'FPRO',8x,'SIGMA',6x,' NEW FA',6x,'OLD FA' / &
           6(1x, 1pe11.3) / 2x,74('-') / )


! 11. Convergence test

    if (lmflag == 0) then

      if ( dxnorm <= epcon ) then

        convergence = .true.
        exit

      endif
    else 
      normDx = sqrt(sum((LMDiagonal *PARAM)**2)) 
      if ((lvdelta <= XTOL* normDx .and. (normJp/normf)**2 + 2 * (sqrt(lvlambda) * normDp /normf)**2 <= FTOL)&
          .or. normf < FTOL) then !or is needed in case f has zero as  residual
        convergence = .true.
        exit
      endif
    endif 
    eqold = equality
    lsold = leastsquares

  enddo

  if ( convergence ) then

    write(*,130) iter, ifct, leastsquares, equality
130 format('', / 2x,71('*') / &
           2x,'*',3x,'PAREMERA converges after ',i3,' iterations and ',i3, &
           ' function calls.',3x,'*'/ &
           2x,'*',3x,'Least-Squares: ', 1pe13.7,5x,'Equality: ',1pe13.7,10x,&
           '*'/ 2x,71('*') /)

  else

    write(*,131) iter-1
131 format(/ 'PAREMERA diverges after ',i3,' iterations.'/ &
           'Try a different set of initial parameters or increase the number &
            of iterations. ' / )

999 format(90(1x,1pe15.5))

  endif

  if ( lmflag /= 0 ) then
    deallocate(JACOBIANmatrix)
    deallocate(RHSvector)
    deallocate(LMDIAGONAL)
  endif

  deallocate(ipossr)
  deallocate(iposx)
  deallocate(IBEL)
  deallocate(TAU)
  deallocate(IPC)
  deallocate(IPR)
  deallocate(IPL)
  deallocate(EXPWISELS)
  deallocate(EXPWISEEQ)
  deallocate(PARAM)
  deallocate(X0)
  deallocate(DXRhilf)
  deallocate(DXhilf)
  deallocate(DXR)
  deallocate(DX)
  deallocate(SCAL)
  deallocate(Xhilf)
  deallocate(RHShilf)
  deallocate(JACOBIANhilf)
  deallocate(GLOBALMATRIX)
  deallocate(RHS)
  deallocate(JACOBIAN)
  deallocate(MRhilf)
  deallocate(MPhilf)
  deallocate(MShilf)
  deallocate(MR)
  deallocate(MP)
  deallocate(MS)

end ! subroutine paremera

!----------------------------------------------------------------------------------!

subroutine scaljacobian(JACOBIAN, SCAL, NBED, DIMS, nglob, nexp, maxnb, maxsr, &
                        cond0)

  implicit none

  integer :: maxsr, &
             maxnb, &
             nexp, &
             nglob, &
             DIMS(nexp), &
             NBED(nexp)

  real*8 :: SCAL(nglob), &
            JACOBIAN(maxnb,maxsr+nglob,nexp), &
            cond0

  integer :: i,j,k
  real*8 :: minnorm

  call dinit2(nglob, 0.d0, SCAL, 1)

  minnorm = 1.d10

  do i = 1,nexp
    do j = 1,nglob
      do k = 1,NBED(i)
        SCAL(j) = SCAL(j) + JACOBIAN(k,DIMS(i)+j,i)*JACOBIAN(k,DIMS(i)+j,i)
      enddo
    enddo
  enddo

  do j = 1,nglob
    if ( SCAL(j) < minnorm ) minnorm = SCAL(j)
  enddo

  minnorm = dsqrt( minnorm )

  do i = 1,nglob

    if ( minnorm < 1.d0/cond0 ) then!     if ( SCAL(i) < 1.D0/cond0) then
      SCAL(i) = 1.d0
    else
      SCAL(i) = 1.d0/dsqrt( SCAL(i) )
    endif

  enddo

  do i = 1,nexp
    do j = 1,nglob
      do k = 1,NBED(i)
        JACOBIAN(k,DIMS(i)+j,i) = JACOBIAN(k,DIMS(i)+j,i)*SCAL(j)
      enddo
    enddo
  enddo

999 format(20(1x,E20.10))


end ! subroutine scaljacobian

!----------------------------------------------------------------------------------!

!> @routine dinit2
!> @brief initializing the variable X with value C
!> @author Robert Kircheis

subroutine dinit2(n, c, X, ix)

  implicit none

  integer :: n, &
             ix
  real*8  :: c, &
             X(n*ix - ix + 1)

  integer :: j, &
             jx

  if( ix == 1 ) then

    do j = 1,n
      X(j) = c
    enddo

  else

    jx = 1

    do j = 1,n

      X(jx) = c
      jx    = jx + ix

    enddo

  endif

end ! subroutine dinit2
