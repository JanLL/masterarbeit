!  +--METSE -----------------------------------------------------------+
!  |  task: connection from Multiple Experiment To Single Experiment   |
!  |        call of Intdiff for every single Experiment                |
!  +-------------------------------------------------------------------+

!> \routine metse
!> \brief connection from multiple experiments to a single experiment
!! call integration and differentiation routine for every single experiment
!> \author Robert Kircheis

subroutine metse2(iter, nexp, X, iposx, ipossr, JACOBIAN, RHS, MS, MP, MR, &
                  EQ, LS, NBED, nglob, multidimr, maxnd, maxnv, maxnt, maxms, &
                  maxsp, maxnb, maxsr, nvtot, rmtflag, IBEL, sflag, pr, infor)


  use omp_lib

  implicit none

  integer :: iter, &
             nexp, &
             iposx(nexp), &
             ipossr(nexp), &
             nglob, &
             multidimr, &
             maxnd, &
             maxnv, &
             maxnt, &
             maxms, &
             maxsp, &
             maxnb, &
             maxsr, &
             nvtot, &
             rmtflag, &
             sflag, &
             pr, &
             infor
  integer :: NBED(nexp), &
             IBEL(maxnv,maxms,nexp)
  real*8 :: X(nvtot), & !< current iterate
            JACOBIAN(maxnb,maxsr+nglob,nexp), &
            RHS(maxnb,nexp), &
            MS(maxnv,maxsr,maxms,nexp), &
            MP(maxnv,nglob,maxms,nexp), &
            MR(maxnv,maxms,nexp), &
            EQ(nexp), &
            LS(nexp)

  integer :: i,j,k, & !< numerators
             iposp, & !< pointer to the paramteres in array X
             iposx2, & !< pointer to shooting nodes in array X
                       !! auxiliary variables
             nms, & !< number of shooting nodes for the current experiment
             nsp, & !< number of 'stoping' points for the current experiment
             nmeas, & !< number of measurements for the current experiment
             nvar, & !< number of states for the current experiment
             ndiff, & !< number of differential states for the current experiment
             np, & !< number of parameters
             dimr, & !< number of rfcns ( local + global ) for the current experiment
             dimrloc, & !< number of rfcn functions for the current experiment
             dims, & !< length of the s-vector (-nalg) for the current experiment
             nt, & !< number of time points for the current experiment
             nhelp, & !< auxiliary variables for allocation
             iflag(nexp) !< integration error flag
  real*8, allocatable :: T(:), & !< array of time points; dim(maxnt)
                         MSTIMES(:), & !< array of shooting time points; dim(maxms)
                         STATES(:,:), & !< shooting variables for the current
                                        !! experiment; dim(maxnv,maxms)
                         PARAM(:), & !< parameters; dim(nglob)
                         SR(:), & !< Sr variables of the current experiment; dim(maxsr)
                         EQR(:), & !< violation of the rfcns
                         multi_cont(:) !< violation of the multirfcns
  real*8 :: cont, &!< violation of the rfcns
            finish, start !< for time measurement
  logical :: return_from_OMP = .false. !< exit OMP Loop on error


  allocate( multi_cont(multidimr) )
  allocate( PARAM(max0(1,nglob)) )

  call dinit2(maxnb*(maxsr+nglob)*nexp, 0.D0, JACOBIAN,1)
  call dinit2(maxnb*nexp, 0.D0, RHS, 1)
  call dinit2(maxnv*maxsr*maxms*nexp, 0.D0, MS, 1)
  call dinit2(maxnv*nglob*maxms*nexp,0.D0, MP, 1)
  call dinit2(maxnv*maxms*nexp, 0.D0, MR, 1)

  infor = 0

  cont = 0.d0

  do i=1,multidimr
    multi_cont(i) = 0.d0
  enddo

  do i=1,nexp
    iflag(i) = 0
  enddo

  if ( sflag == 1 ) then

    write(*,*) 'Reading variables from file.'
    call vpl_read_x_from_file(nvtot, X)

  endif

  call getnp(np, iposp)

  if ( nglob > 0 ) call dcopy(nglob, X(iposp), 1, PARAM, 1)

  if ( ( rmtflag == 0 ).and.( pr > 0 ) ) then

    write(*,*) 'Parameter:'
    do i=1,nglob
      write(*,111) I,':',PARAM(i)
    enddo
    write(*,*)
111  format(1x,i4,a1,2x,1pe15.5)

  endif

! 1. Run integration and differentiation for each experiment

  call cpu_time(start)

  !$OMP PARALLEL DEFAULT(SHARED) &
  !$OMP&         PRIVATE( i, j, k, iposx2, STATES, T, MSTIMES, SR, &
  !$OMP&         ndiff, nvar, dimr, dimrloc, dims, nmeas, EQR, nt, &
  !$OMP&         nms, nsp, np, multi_cont )

  !$OMP DO SCHEDULE(DYNAMIC)
  do i = 1,nexp

    if ( .not.return_from_OMP ) then

      if ( ( rmtflag == 0 ).and.( pr > 1 ) ) write(*,222) i
222  format('Integrating and Differentiating Experiment ',I3)

      call getdims(i-1, nms, nsp, nmeas)
      call getndyn(i-1, nvar, ndiff)

      nt = nms + nsp
      np = nglob

      allocate( T(nt) )
      allocate( MSTIMES(nms) )
      allocate( STATES(nvar,nms) )

      call dinit2(nms, 0.D0, MSTIMES, 1)
      call dinit2(nt, 0.D0, T, 1)
      call dinit2(nvar*nms, 0.d0, STATES, 1)

      call getdim_rfcn(i-1, dimrloc)
      call getdim_s(i-1, dims)

      dimr = dimrloc + multidimr
      dims = dims - ( nvar - ndiff )

      nhelp = max0(1,dims)

      allocate( SR(nhelp) )
      allocate( EQR(dimr) )

      do j = 1,dimr
        EQR(j) = 0.d0
      enddo

      if ( dims > 0 ) call dinit2(dims, 0.d0, SR, 1)

      do j = 1,nvar
        IBEL(j,1,i) = 1
      enddo

      call getmstimes(i-1, MSTIMES, nms)

      call gettimes(i-1,T,nt)

      nt = nt - 1

      iposx2 = iposx(i)

      do j = 1,nms

        call dcopy(nvar, X(iposx2), 1, STATES(1,j), 1)

        iposx2 = iposx2 + nvar

      enddo

      if ( dimr > 0 ) call dcopy(dims, X(ipossr(i)), 1, SR(1), 1)

      EQ(i) = 0.D0
      LS(i) = 0.D0

      if ( (iter == 1) .and. (sflag == 0) ) then
        call initmsvals(i-1, STATES, nvar, IBEL(1,1,i), maxnv)
      endif

      call intdiff( i, T, MSTIMES, STATES, PARAM, SR, ndiff, nvar, dimr, &
                    dimrloc, multidimr, dims, nmeas, NBED(i), nt, nms, nsp, &
                    JACOBIAN(1,1,i), RHS(1,i), MS(1,1,1,i), MP(1,1,1,i), &
                    MR(1,1,i), EQ(i), LS(i), EQR, np, IBEL(1,1,i), &
                    maxnb, maxsr, maxnv, maxms, iflag(i) )

      if ( iflag(i) < 0 ) then

        infor = -i

        write(*,444) i

        !$OMP CRITICAL
        return_from_OMP = .true.
        !$OMP END CRITICAL
      endif
444 format('*** Error in intdiff for Experiment ***',i3)

      iposx2 = iposx(i)

      do j = 1,nms

        call dcopy(nvar, STATES(1,j), 1, X(iposx2), 1)

        iposx2 = iposx2 + nvar

      enddo

      do j = 1,dimrloc
        EQ(i) = EQ(i) + EQR(j)*EQR(j)
      enddo

      !$OMP CRITICAL
      do j = 1,multidimr
        multi_cont(j) =  multi_cont(j) + EQR(dimrloc+j)
      enddo
      !$OMP END CRITICAL

      deallocate(EQR)
      deallocate(SR)
      deallocate(STATES)
      deallocate(MSTIMES)
      deallocate(T)

    endif

  enddo
  !$OMP END DO NOWAIT

  !$OMP END PARALLEL

  do i=1,multidimr
    cont = cont + multi_cont(i)*multi_cont(i)
  enddo

  EQ(nexp) = EQ(nexp) + cont

  call cpu_time(finish)
  if ( pr > 0 ) print '("Time = ",1pe9.3," seconds.")',finish - start

  deallocate(PARAM)
  deallocate(multi_cont)

999 format(40(1x,1pe20.10))

end ! subroutine         metse2
