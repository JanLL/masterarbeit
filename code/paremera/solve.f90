!  +-COMPREDINC-----------------------------------------------------------+
!  | tasks: Compute the increments of the reduced system                  |
!  |                  min |u_1 + E_1*s + P_1*p|                           |
!  |                s.t.   u_2 + E_2*s + P_2*p = 0                        |
!  |        make use of the decomposition computed in DECOMPOSE           |
!  +----------------------------------------------------------------------+

!> \routine compredinc
!> \brief compute the increments of the reduced system by using the decomposed
!! Jacobian
!> \return DXR
!> \author Robert Kircheis

subroutine compredinc(JACOBIAN, maxnb, GLOBALMATRIX, nbedtotal, PARAM, maxsr, &
                      maxsloc, PSEUDOINVERSE, rank, TAU, SCAL, RHS, DXR, &
                      dimrtot, DIMR, DIMRLOC, DIMS, multidimr, nglob, NBED, &
                      nexp, IPL, IPR, IPC, LMDIAGONAL, lvdelta, lvsigma,&
                      lvlambda, RHSvector, JACOBIANmatrix, expsum, lmflag, kout)

  implicit none

  integer :: maxnb, &
             nbedtotal, &
             maxsr, &
             maxsloc, &
             rank, &
             DIMR(nexp), &
             DIMRLOC(nexp), &
             DIMS(nexp), &
             multidimr, &
             dimrtot, &
             nglob, &
             NBED(nexp), &
             nexp, &
             IPL(maxsr,nexp), &
             IPR(multidimr), &
             IPC(nglob), &
             kout, &
             expsum, & 
             lmflag
  real*8 :: JACOBIAN(maxnb,maxsr+nglob,nexp), &
            GLOBALMATRIX(nbedtotal,multidimr+nglob), &
            PARAM(nglob), &
            PSEUDOINVERSE(nglob,nglob), &
            TAU(nglob), &
            SCAL(nglob), &
            RHS(maxnb,nexp), &
            DXR(dimrtot+nglob), &
            LMDIAGONAL(nglob), &
            lvdelta, &
            lvsigma, &
            lvlambda, &
            RHSvector(expsum), &
            JACOBIANmatrix(expsum,nglob)

  integer :: i,j,k,l, & !< numerators
             ipos, ipos2, nhilf !< position pointer
  real*8 :: eta !< auxiliary variable
  real*8, allocatable :: GRHS(:), & !< corresponding righthand side to GLOBALMATRIX
                         HILF(:), & !< auxiliary array for intermediate results
                         DELTAP(:), & !< increment for the parameters p
                         HELP(:,:), & !< auxiliary variable
                         DXM(:), & !< increment for s^m
                         DXS(:),& !< increment for s
                         GLOBALHELP(:,:), & !< auxiliary variable
                         HELP2(:), &  !< auxiliary variable
                         COPYLMD (:) !< copy of LMDIAGONAL

  call dinit2(dimrtot+nglob, 0.d0, DXR, 1)

  nhilf = max0(nbedtotal,nglob)


  allocate( GRHS(nhilf) )
  allocate( HILF(nhilf) )
  allocate( DELTAP(nglob) )
  allocate( HELP(nglob,nglob) )
  allocate( DXM(multidimr) )
  allocate( DXS(dimrtot) )
  allocate( GLOBALHELP(nbedtotal,nglob))
  allocate( HELP2(nglob))
  allocate( COPYLMD(nglob))

  ipos = multidimr

  call dinit2(nhilf,0.d0,GRHS,1)
  call dinit2(nglob,0.d0,DELTAP,1)
  call dinit2(multidimr,0.d0,DXM,1)
  call dinit2(dimrtot,0.d0,DXS,1)

  ! adjust right-hand side to LU - decomposed system
  ! and build-up the 'global' right-hand side

  do l= 1,nexp

    do k = 1,DIMRLOC(l)

      do j = k+1,NBED(l)

        eta = JACOBIAN(j,k,l)/JACOBIAN(k,k,l)

        RHS(j,l) = RHS(j,l) - eta*RHS(k,l)

      enddo

    enddo

    do j = 1,multidimr
      GRHS(j) = GRHS(j) + RHS(DIMRLOC(l)+j,l)
    enddo

    do j = 1,NBED(l)-DIMR(l)
      GRHS(ipos+j) = RHS(DIMR(l)+j,l)
    enddo

    ipos = ipos + NBED(l) - DIMR(l)

  enddo
   
  call initpr(PARAM, HELP, nglob, GRHS(ipos+1))

  ! adjust 'global' right-hand side to LU - decomposed system

  do k = 1,multidimr

    if ( IPR(k) /= k ) then       !PIVOTING

      eta          = GRHS(k)
      GRHS(k)      = GRHS(IPR(k))
      GRHS(IPR(k)) = eta

    endif

    do j = k+1,nbedtotal

      eta = GLOBALMATRIX(j,k)/GLOBALMATRIX(k,k)

      GRHS(j) = GRHS(j) - eta*GRHS(k)

    enddo

  enddo

  ! adjust GRHS to QR - decomposition

  do i = 1,rank   ! F' = Q^T * F

    call dcopy(nbedtotal,GRHS,1,HILF,1)

    do j = multidimr+i,nbedtotal
      do k = multidimr+i,nbedtotal
        HILF(j) = HILF(j) &
                 - GLOBALMATRIX(j,multidimr+i) &
                  * GLOBALMATRIX(k,multidimr+i) &
                  * GRHS(k)
      enddo
    enddo

    call dcopy(nbedtotal,HILF,1,GRHS,1)

  enddo

  call dinit2(nbedtotal,0.D0,HILF,1)
  
  ! Compute delta_p

  if ( rank == nglob ) then ! use QR-decomposition

    do i = nglob,1,-1    ! delta_p = R^(-1) * F'

      DELTAP(i) = GRHS(multidimr+i)

      do j = i+1,nglob
        DELTAP(i) = DELTAP(i) - GLOBALMATRIX(multidimr+i,multidimr+j) &
                               * DELTAP(j)
      enddo

      DELTAP(i) = DELTAP(i)/TAU(i)

    enddo

  else ! use Deuflhard's algorithm for rank-deficient systems

    do i =  rank,1,-1  ! b' = U_q^(-1) * b

      do j = i+1,rank
        GRHS(multidimr+i) = GRHS(multidimr+i) &
                           - GLOBALMATRIX(multidimr+i,multidimr+j) &
                            * GRHS(multidimr+j)
      enddo

      GRHS(multidimr+i) = GRHS(multidimr+i)/TAU(i)

    enddo

    do i = rank+1,nglob
      GRHS(multidimr+i) = 0.d0
    enddo

    do i = 1,nglob-rank  ! b'' = (j'^T -I) * b'
      do j = 1,rank
        HILF(i) = HILF(i) + GLOBALMATRIX(multidimr+j,multidimr+rank+i) &
                           * GRHS(j)
      enddo
    enddo

    ! b''' = (L*L^T)^(-1) * b''

    do i = 1,nglob-rank

      do k = 1,i-1
        HILF(i) = HILF(i) - PSEUDOINVERSE(i,k)*HILF(k)
      enddo

      HILF(i) = HILF(i)/PSEUDOINVERSE(i,i)

    enddo

    do i = nglob-rank,1,-1

      do k = nglob-rank,i+1,-1
        HILF(i) = HILF(i) - PSEUDOINVERSE(k,i)*HILF(k)
      enddo

      HILF(i) = HILF(i)/PSEUDOINVERSE(i,i)

    enddo

    ! delta_p = (j'^T -I)^T * b'''

    do i = 1,rank
      do j = 1,nglob-rank
        DELTAP(i) = DELTAP(i) + GLOBALMATRIX(multidimr+i,multidimr+rank+j) &
                               * HILF(j)
      enddo
    enddo

    do i = 1,nglob-rank
      DELTAP(rank+i) = -HILF(i)
    enddo

    do i = 1,nglob

      if ( i <= rank ) then
        DELTAP(i) = GRHS(multidimr+i) - DELTAP(i)
      else
        DELTAP(i) = -DELTAP(i)
      endif

    enddo

  endif

  do i = 1,nglob
    DXR(IPC(i)) = DELTAP(i)*SCAL(IPC(i))
  enddo

  ! Compute delta_s^r for each experiment

  do i = multidimr,1,-1

    DXM(i) = GRHS(i)

    do j = 1,nglob
      DXM(i) = DXM(i) - GLOBALMATRIX(i,multidimr+j)*DELTAP(j)
    enddo

    do j = i+1,multidimr
      DXM(i) = DXM(i) - GLOBALMATRIX(i,j)*DXM(j)
    enddo

    DXM(i) = DXM(i)/GLOBALMATRIX(i,i)

  enddo

  ipos  = 0
  ipos2 = 0

  do l = 1,nexp

    do i = DIMRLOC(l),1,-1

      DXS(ipos+i) = RHS(i,l)

      do j = 1,nglob
        DXS(ipos+i) = DXS(ipos+i) - JACOBIAN(i,DIMS(l)+IPC(j),l)*DELTAP(j)
      enddo

      do j = 1,DIMS(l)-DIMRLOC(l)
        DXS(ipos+i) = DXS(ipos+i) - JACOBIAN(i,DIMRLOC(l)+j,l)*DXM(ipos2+j)
      enddo

      do j = i+1,DIMRLOC(l)
        DXS(ipos+i) = DXS(ipos+i) - JACOBIAN(i,j,l)*DXS(ipos+j)
      enddo

      DXS(ipos+i) = DXS(ipos+i)/JACOBIAN(i,i,l)

    enddo

    do i = 1, DIMS(l) - DIMRLOC(l)
      DXS(ipos+DIMRLOC(l)+i) = DXM(ipos2+i)
    enddo

    ipos  = ipos + DIMS(l)
    ipos2 = ipos2 + DIMS(l) - DIMRLOC(l)

  enddo

  ! store everything in DXR???

  ipos  = nglob
  ipos2 = 0

  do l = 1,nexp

    do i = 1,DIMS(l)
      DXR(ipos+IPL(i,l)) = DXS(ipos2+i)
    enddo

    ipos  = ipos + DIMS(l)
    ipos2 = ipos2 + DIMS(l)

  enddo

  ! LM part
  if(lmflag == 1) then

    do i =1,nglob
      HELP2(i) = (LMDIAGONAL(i)*DELTAP(IPC(i)))**2
    enddo


    if (sqrt(sum(HELP2)) < (1 + lvsigma)*lvdelta  ) then
      lvlambda = 0
      do i = 1,nglob 
        DXR(IPC(i)) = DELTAP(i)
      enddo

    else
    
    !Calculate the Levenberg Marquartd increment


    !GLOBALHELP becomes (R 0 0)' in lambdasolve  -> (R 0 D_alpha)

      GLOBALHELP = 0.d0
      do i = 1,nglob
        do j = i, nglob
          GLOBALHELP(i, j) = GLOBALMATRIX(i, j)
        enddo
      enddo

      do i = 1,nglob
        GLOBALHELP(i,i) = tau(i)
      enddo

      call lambdasolve(JACOBIANmatrix, RHSvector, GRHS, GLOBALHELP, LMDIAGONAL,DELTAP,&
                       IPC, lvdelta, lvsigma, nbedtotal, nglob, expsum, lvlambda) 
      do i = 1,nglob 
        DXR(i) = DELTAP(i)
      enddo

    endif



  endif
  !End LM part

  do i = 1,dimrtot+nglob
    DXR(i) = - DXR(i)
  enddo

999 format(200(1x,1pe20.10))


  deallocate( HELP2)
  deallocate( GLOBALHELP)
  deallocate( DXS )
  deallocate( DXM )
  deallocate( HELP )
  deallocate( DELTAP )
  deallocate( HILF )
  deallocate( GRHS )

END ! subroutine compredinc

!----------------------------------------------------------------------------------!
!  +-COMPINC----------------------------------------------------------------+
!  | tasks: Compute the remaining increments                                |
!  |        delta_s = MS*delta_sr + MP*delta_p + MR                         |
!  +------------------------------------------------------------------------+

!> \routine compinc
!> \brief computes the remaining increments
!> \return DX
!> \author Robert Kircheis

subroutine compinc(DX, DXR, MS, MP, MR, NMS, NVAR, NDIFF, DIMS, multidimr, &
                   nglob, NBED, nexp, nvtot, dimrtot, maxnv, maxsr, maxms, &
                   kout)

  implicit none

  integer :: NMS(nexp), &
             NVAR(nexp), &
             NDIFF(nexp), &
             DIMS(nexp), &
             multidimr, &
             nglob, &
             NBED(nexp), &
             nexp, &
             nvtot, &
             dimrtot, &
             maxnv, &
             maxsr, &
             maxms, &
             kout
  real*8 :: DX(nvtot), &
            DXR(dimrtot+multidimr+nglob)

  real*8 :: MS(maxnv, maxsr, maxms, nexp), &
            MP(maxnv, nglob, maxms, nexp), &
            MR(maxnv, maxms, nexp)
  integer :: ipos, ipos2, & !< position pointer
             nms1, & ! help variable for algebraic problems
             i,j,k,l !< numerators


! delta_s^j = MS^(j-1) * delta_s^r + MP^(j-1) * delta_p + MR^(j-1)
! for j = 0, nms

  ipos    = 0
  ipos2   = nglob

  call dinit2(nvtot,0.d0,DX,1)

  do l = 1,nexp

    nms1 = max0( NMS(l)-1, 1 )

    do k = 1,nms1

      do j = 1,NVAR(l)

        DX(ipos+j) = MR(j,k,l)

        do i = 1,nglob
          DX(ipos+j) = DX(ipos+j) + MP(j,i,k,l)*DXR(i)
        enddo

        do i = 1,DIMS(l)
          DX(ipos+j) = DX(ipos+j) + MS(j,i,k,l)*DXR(ipos2+i)
        enddo

      enddo

      if ( NMS(l) > 1 ) ipos = ipos + NVAR(l)

    enddo

    ipos  = ipos + NVAR(l)
    ipos2 = ipos2 + DIMS(l)

  enddo

  do i = 1,dimrtot+nglob
    DX(ipos+i) = DXR(i)
  enddo

999 format(20(1x,1pe20.10))

END ! subroutine compinc
