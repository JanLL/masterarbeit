!  +-Decompose-----------------------------------------------------------+
!  | tasks: Decompose the Jacobian                                       |
!  |        Use LU-Decomposition for the sr-block and a QR-Decomposition |
!  |        for the p-block                                              |
!  +---------------------------------------------------------------------+

!> \routine decjacobian
!> \brief use a suited combination of LU and QR decomposition to decompose
!! the Jacobian
!> \author Robert Kircheis

subroutine decjacobian(JACOBIAN, maxnb, GLOBALMATRIX, nbedtotal, PARAM, SCAL, &
                       PSEUDOINVERSE, TAU, rank, DIMR, DIMRLOC, DIMS, &
                       multidimr, maxsr, maxsloc, nglob, NBED, nexp, cond0, &
                       IPL, IPR, IPC, pr, kout)

  implicit none

  integer :: rank, &
             DIMR(nexp), &
             DIMRLOC(nexp), &
             DIMS(nexp), &
             multidimr, &
             nglob, &
             NBED(nexp), &
             nexp, &
             maxnb, &
             maxsr, &
             maxsloc, &
             nbedtotal, &
             pr, &
             kout, &
             IPL(maxsr,nexp), &
             IPR(multidimr), &
             IPC(nglob)
  real*8 :: JACOBIAN(maxnb,maxsr+nglob,nexp), &
            GLOBALMATRIX(nbedtotal,multidimr+nglob), &
            PARAM(nglob), &
            PSEUDOINVERSE(nglob,nglob), &
            TAU(nglob), &
            SCAL(nglob), &
            cond0

  integer :: i, j, k, l, & !< numerators
             ipos, iposs, jj, kk, & !< position pointer, numerators
             hilfrank !< rank of the jacobian, determined during decomposition
  real*8 :: eta, & !< auxiliary variable for the Householder transformation
            pivot, &
            s, & !< auxiliary variable for the Householder transformation
            norm !< auxiliary variable for the Householder transformation
  real*8, allocatable :: HELP(:) !< auxiliary array
  
  allocate( HELP(nglob) )

  call dinit2(nglob*nglob,0.D0,PSEUDOINVERSE,1)
  call dinit2(nbedtotal*(multidimr+nglob),0.D0,GLOBALMATRIX,1)
  call dinit2(nglob,0.D0,TAU,1)

  do i = 1,nglob
    IPC(i) = i
  enddo

  hilfrank = rank

  !LU-Decomposition on local rfcn-block

  ipos  = multidimr
  iposs = 0

  do l = 1,nexp

    do k = 1,DIMS(l)
      IPL(k,l) = k
    enddo

    do k = 1,DIMRLOC(l)

      jj = k

      pivot = dabs( JACOBIAN(k,k,l) )

      do i = k+1, DIMS(l)

        if ( pivot < dabs( JACOBIAN(k,i,l) ) ) then

          pivot = dabs( JACOBIAN(k,i,l) )

          jj = i

        endif

      enddo

      if ( pivot < 1.e-08 ) then

        write(*,111) l, DIMRLOC(l)

        kout = -l
        return
      endif
111 format('For experiment',i2,' pivot element equal to zero.'/ &
           'The local rfcns do not have rank',i3,'.')

      if ( IPL(k,l) /=  jj ) then

        kk        = IPL(jj,l)
        IPL(jj,l) = IPL(k,l)
        IPL(k,l)  = kk

        call DSWAP(maxnb, JACOBIAN(1,k,l), 1, JACOBIAN(1,jj,l), 1)

      endif

      do j = k+1,NBED(l)

        eta = JACOBIAN(j,k,l)/JACOBIAN(k,k,l)

        do i = k+1,DIMS(l)+nglob
          JACOBIAN(j,i,l) = JACOBIAN(j,i,l) - eta*JACOBIAN(k,i,l)
        enddo

      enddo

    enddo

    do i = 1,DIMS(l)-DIMRLOC(l)
  
      do j = 1,multidimr ! copy dr^m/ds
        GLOBALMATRIX(j,iposs+i) = JACOBIAN(DIMRLOC(l)+j,DIMRLOC(l)+i,l)
      enddo

      do j = 1,NBED(l)-DIMR(l) ! copy dg/ds, dh/ds

        if ( ipos + j > nbedtotal) then

          write(*,*) 'Error: Globalmatrix too small)'

          kout = -l
          return
        endif

        GLOBALMATRIX(ipos+j,iposs+i) = JACOBIAN(DIMR(l)+j,DIMRLOC(l)+i,l)

      enddo

    enddo

    iposs = iposs + DIMS(l) - DIMRLOC(l)

    do i = 1,nglob ! copy dg/dp, dh/dp

      do j = 1,NBED(l)-DIMR(l)
        GLOBALMATRIX(ipos+j,multidimr+i) = JACOBIAN(DIMR(l)+j,DIMS(l)+i,l)
      enddo

    enddo

    ipos = ipos + NBED(l) - DIMR(l)

    do i = 1,nglob ! copy dr^m/dp

      do j = 1,multidimr
        GLOBALMATRIX(j,multidimr+i) = GLOBALMATRIX(j,multidimr+i)         &
                                     + JACOBIAN(DIMRLOC(l)+j,DIMS(l)+i,l)
      enddo

    enddo

  enddo

  call initpr(PARAM, GLOBALMATRIX(ipos+1,multidimr+1), nbedtotal, HELP)

  ! LU-Decomposition on multirfcn-block

  do k = 1,multidimr

    IPR(k) = k

    pivot = dabs( GLOBALMATRIX(k,k) )

    do j = k+1,multidimr

      if ( pivot < dabs( GLOBALMATRIX(j,k) ) ) then

        pivot = dabs( GLOBALMATRIX(j,k) )
        IPR(k)    = j

      endif

    enddo

    if ( pivot < 1.e-08 ) then

      write(*,112)

      kout = -k
      return
    endif
112 format('Globalmatrix does not have full rank for multirfcns.' / &
           'Please check your model!')

    if ( IPR(k) /=  k ) then
      call DSWAP(multidimr+nglob, GLOBALMATRIX(k,1), nbedtotal, &
                 GLOBALMATRIX(IPR(k),1), nbedtotal)
    endif

    do j = k+1,nbedtotal

      eta = GLOBALMATRIX(j,k)/GLOBALMATRIX(k,k)

      do i = k+1,multidimr+nglob
        GLOBALMATRIX(j,i) = GLOBALMATRIX(j,i) - eta*GLOBALMATRIX(k,i)
      enddo

    enddo

  enddo

  ! QR-Decomposition on p-block

  do i = 1,rank ! do 1 = 1,nglob

    norm   = 0.D0
    jj     = i

    do j = multidimr+i,nbedtotal
      norm   = norm + GLOBALMATRIX(j,multidimr+i) * GLOBALMATRIX(j,multidimr+i)
    enddo

    do j = i+1,nglob

      s = 0.D0

      do k = multidimr+i,nbedtotal
        s = s + GLOBALMATRIX(k,multidimr+j) * GLOBALMATRIX(k,multidimr+j)
      enddo

      if ( s > norm ) then
        norm = s
        jj   = j
      endif

    enddo

    norm = dsqrt(norm)

! If Globalmatrix is ill-conditioned switch to Deuflhard's algorithm

    if ( norm < (1.D0/cond0) ) then

      hilfrank = i-1

      if ( hilfrank == 0 ) then
        write(*,333) 
333 format( / 'Error: Jacobimatrix has rank zero.' / '')
        kout = -1000
        return
      endif

      if ( pr > 0 ) write(*,222) IPC(i)
222 format(/ 2x,49('*') / &
             2x,'* Warning: Matrix rank deficent.',16x,'*'/ &
             2x,'* Parameter ',i2,' et seqq. poorly estimable.',7x,'*'/ &
             2x,'* Continue by using a rank reduction algorithm. *'/ &
             2x,49('*'))

! The scaling has to be revoked because of the rank deficiency
! No, it is working fine with scaling. It just might be the case,
! that the rank-deficiency is related to another parameter!

      do j = i,nglob
        do k = 1,rank
          GLOBALMATRIX(multidimr+k,multidimr+j) &
                    = GLOBALMATRIX(multidimr+k,multidimr+j)/SCAL(IPC(j))
        enddo
      enddo

      do j = 1,i-1

        do k = 1,j-1
          GLOBALMATRIX(multidimr+k,multidimr+j) &
                    = GLOBALMATRIX(multidimr+k,multidimr+j)/SCAL(IPC(j))
        enddo

        TAU(j) = TAU(j)/SCAL(IPC(j))

      enddo

      do j = 1,nexp
        do k = 1,nglob
          do l = 1,dims(j)
            JACOBIAN(l,dims(j)+k,j) = JACOBIAN(l,dims(j)+k,j)/SCAL(k)
          enddo
        enddo
      enddo

      do j = 1,nglob
        SCAL(j) = 1.D0
      enddo

      exit
    endif

    if ( IPC(i) /= jj) then

      kk      = IPC(jj)
      IPC(jj) = IPC(i)
      IPC(i)  = kk

      call DSWAP(nbedtotal, GLOBALMATRIX(1,multidimr+i), 1, &
                 GLOBALMATRIX(1,multidimr+jj), 1)

    endif

    if ( GLOBALMATRIX(multidimr+i,multidimr+i) >= 0.D0 ) then
      TAU(i) = -norm
    ELSE
      TAU(i) = norm
    endif

    eta = dsqrt( norm*( norm + dabs(GLOBALMATRIX(multidimr+i,multidimr+i)) ) )

    GLOBALMATRIX(multidimr+i,multidimr+i) &
              = GLOBALMATRIX(multidimr+i,multidimr+i) - TAU(i)

    do j = i,nbedtotal-multidimr
      GLOBALMATRIX(multidimr+j,multidimr+i) &
                = GLOBALMATRIX(multidimr+j,multidimr+i)/eta
    enddo

    do j = i+1,nglob

      s = 0

      do k = i,nbedtotal-multidimr
        s = s + GLOBALMATRIX(multidimr+k,multidimr+i) &
               * GLOBALMATRIX(multidimr+k,multidimr+j)
      enddo

      do k = multidimr+i,nbedtotal
        GLOBALMATRIX(k,multidimr+j) = GLOBALMATRIX(k,multidimr+j) &
                                     - GLOBALMATRIX(k,multidimr+i)*s
      enddo

    enddo

  enddo

  ! for details see Deuflhard

  rank = hilfrank

  if ( rank < nglob ) then

    !GAUSS

    do i = rank+1,nglob

      do j = rank,1,-1

        do k = j+1,rank
          GLOBALMATRIX(multidimr+j,multidimr+i) &
                    = GLOBALMATRIX(multidimr+j,multidimr+i) &
                     - GLOBALMATRIX(multidimr+j,multidimr+k) &
                      * GLOBALMATRIX(multidimr+k,multidimr+i)
        enddo

        GLOBALMATRIX(multidimr+j,multidimr+i) &
                  = GLOBALMATRIX(multidimr+j,multidimr+i)/TAU(j)

      enddo

    enddo

    !PSEUDO-INVERSE

    do i = 1,nglob-rank

      PSEUDOINVERSE(i,i) = PSEUDOINVERSE(i,i) + 1.D0

      do j = 1,nglob-rank
        do k = 1,rank
          PSEUDOINVERSE(i,j) = PSEUDOINVERSE(i,j) &
                               + GLOBALMATRIX(multidimr+k,multidimr+rank+i) &
                                * GLOBALMATRIX(multidimr+k,multidimr+rank+j)
        enddo
      enddo

    enddo

    !CHOLESKY

    do i = 1,nglob-rank

      do j = 1,i-1

        s = PSEUDOINVERSE(i,j)

        do k = 1,j-1
          s = s - PSEUDOINVERSE(i,k)*PSEUDOINVERSE(j,k)
        enddo

        PSEUDOINVERSE(i,j) = s/PSEUDOINVERSE(j,j)

      enddo

      s = PSEUDOINVERSE(i,i)

      do k = 1,i-1
        s = s - PSEUDOINVERSE(i,k)*PSEUDOINVERSE(i,k)
      enddo

      if ( s <= 0 ) then
        write(*,444)

        kout = -1
        return
      ELSE
        PSEUDOINVERSE(i,i) = dsqrt(s)
      endif
444 format(/ 'Error: PSEUDOINVERSE is not positiv definit.' / '')

    enddo

  endif

  if ( pr > 0 ) write(*,1234) ( IPC(i), i=1,nglob )
1234 format( / 'Order of Parameters:' / 20(1x,i2), '', / '')

999 format(50(1x,1pe20.10))

  deallocate(HELP)

end ! subroutine decjacobian
