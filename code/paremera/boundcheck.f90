! +--BOUNDCHECK-------------------------------------------------------+
! | task: checks if parameter constraints are violated                |
! |       if so the step size is reduced                              |
! +-------------------------------------------------------------------+

!> \routine intdiff
!> \brief checks if parameter constraints are violated
!! if so the step size is reduced so that the new iterate lies within
!! the constraints
!> \author Robert Kircheis

subroutine boundcheck( X1, X0, DX, nvtot, fa, fashort, farel, &
                       PUPBND, PLOBND, nglob, kout )

  implicit none

  integer :: nvtot, &
             nglob, &
             kout
  real*8 :: X1(nvtot), &
            X0(nvtot), &
            DX(nvtot), &
            fa, &
            fashort, &
            farel, &
            PUPBND(nglob), &
            PLOBND(nglob)
  character(len=5) :: bound

  integer :: i, ident, & !< numerator
             iposp !< position of the parameters in the vector X1
  real*8 :: fa1, fa2, & !< auxiliary variables to compute the new step size
            dummy

  call getnp(dummy, iposp)

  fa1 = fa

  do i = 0,nglob-1

    if ( X1(iposp+i) < PLOBND(i+1) ) then

      fa2 = dabs( ( PLOBND(i+1) - X0(iposp+i) )/( DX(iposp+i) ) )

      fa1 = dmin1( fa1, fa2 )

      ident = i+1
      bound = "lower"

    elseif ( X1(iposp+i) > PUPBND(i+1) ) then

      fa2 = dabs( ( X0(iposp+i) - PUPBND(i+1) )/( DX(iposp+i) ) )

      fa1 = dmin1( fa1, fa2 )

      ident = i+1
      bound = "upper"

    endif    

  enddo

  if ( fa1 < fa ) then

    write(*,123) ident, bound
123 format( / 5x, 'Warning: Boundcheck active for parameter', I3,' at ', &
           A5,' bound' / )

    fa = 0.5d0*fa1

    do i = 1,nvtot
      X1(i) = X0(i) + fashort*fa*DX(i)
    enddo

  endif

  if ( fa < 1.e-5 ) then

    write(*,124)

    kout = -111

    return
  endif
124 format( / 5x, 'Abbortion, step size to small' )

end ! subroutine boundcheck