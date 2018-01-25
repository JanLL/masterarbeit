!  +--INTTEST ---------------------------------------------------------+
!  |  task:  Test if integration works for the new iterative           |
!  +-------------------------------------------------------------------+

subroutine inttest( nexp, X, iposx, nvtot, pr, intflag )


  use omp_lib

  implicit none


  integer, intent(in) :: nexp, &
                         iposx(nexp), &
                         nvtot, &
                         pr
  integer :: intflag
  real*8, intent(in) :: X(nvtot)

  integer :: i,j,k,& !< numerators
             iexp, & !< numerator
             iposp, & !< pointer to the paramteres in array X
             iposx2, & !< pointer to shooting nodes in array X
             nms, & !< number of shooting nodes for the current experiment
             nsp, & !< number of 'stoping' points for the current experiment
             nmeas, & !< number of measurements for the current experiment
             nvar, & !< number of states for the current experiment
             ndiff, & !< number of differential states for the current experiment
             np, & !< number of parameters
             nt, & !< number of time points for the current experiment
             idum, & !< auxiliary variable
             iflag(nexp), & !< integration error flag
             start, ende, rate, cmax !< variables for time measurements
  integer, allocatable :: NODE(:) !< array of flags that shows whether the current
                                  !! time point is a MSNode, RNode os MEASNode
  real*8, allocatable :: T(:), & !< array of time points; dim(maxnt)
                         PARAM(:), & !< parameters; dim(nglob)
                         STATES(:,:), & !< shooting variables for the current
                                        !! experiment; dim(nvar,nms)
                         Y(:) !< auxiliary variable for the integration
  real*8 :: ga, gadir, & !< sensitivity and seed matrix; here empty
            h, rdum, & !< auxiliary variables
            time !< run time
  logical :: return_from_OMP = .false. !< exit OMP Loop on error

  ga = 0.d0
  gadir = 0.d0

  rate = 10
  cmax = 1000000

  intflag = 0

  do i = 1,nexp
    iflag(i) = 0
  enddo

  call getnp(np, iposp)

  allocate( PARAM(np) )
  call dinit2(np, 0.d0, PARAM, 1)

  call dcopy(np, X(iposp+1), 1, PARAM, 1)

  call system_clock(start, rate, cmax)

  !$OMP PARALLEL DEFAULT(SHARED) &
  !$OMP&         PRIVATE( iexp, nms, nsp, nmeas, nt, &
  !$OMP&                  T, STATES, Y, NODE, iposx2, i, j, k )

  !$OMP DO SCHEDULE(DYNAMIC)
  do iexp = 1,nexp
    if (.not.return_from_OMP ) then

      if ( pr > 1 ) write(*,1233) iexp
1233 format('Test integration for Experiment ',i3)

      call getdims(iexp-1, nms, nsp, nmeas)
      call getndyn(iexp-1, nvar, ndiff)

      nt = nms + nsp

      allocate( T(nt) )
      allocate( STATES(nvar,nms) )
      allocate( Y(nvar) )
      allocate( NODE(3) )

      call dinit2(nt, 0.d0, T, 1)
      call dinit2(nvar*nms, 0.d0, STATES, 1)
      call dinit2(nvar, 0.d0, Y, 1)

      call gettimes(iexp-1, T, nt)

      nt = nt - 1

      iposx2 = iposx(iexp)

      do i = 1,nms

        call dcopy(nvar, X(iposx2), 1, STATES(1,i),  1)

        iposx2 = iposx2 + nvar

      enddo

      j = 1
      k = 1

      call getnodetype( iexp-1, j-1, NODE )

      do while ( k < nms )

        if ( NODE(1) == 1 ) call dcopy(nvar, STATES(1,k), 1, Y, 1)

        h = 0.0d0

        call vpldaepar( iexp, j, nvar, ndiff, nvar+np, gadir, 1, 0, T(j), &
                        T(j+1), Y, PARAM, 1.d-06, 0.d0, h, 0.d0, 0, &
                        ga, 1, rdum, idum, iflag(iexp) )

        if ( iflag(iexp) < 0 ) then

          intflag = -iexp

          write(*,444) iexp

          !$OMP CRITICAL
          return_from_OMP = .true.
          !$OMP END CRITICAL
        endif
444 format('*** Error in vpldaepar for Experiment ***',i3)

        j = j + 1

        call getnodetype( iexp-1, j-1, NODE)

        if ( NODE(1) == 1 ) k = k+1

    enddo

    deallocate(NODE)
    deallocate(Y)
    deallocate(STATES)
    deallocate(T)

    endif

  enddo
  !$OMP END DO NOWAIT

  !$OMP END PARALLEL

  call system_clock(ende)

  time = float(ende - start)/float(rate)   ! evtl. cmax beachten!

!   if ( pr > 1 ) write(*,333) time
! 333 format(/ 'Time in seconds for inttest: ',1pe10.3, '', / '' )

  deallocate(PARAM)

999 format(1pe20.10)

end ! subroutine inttest