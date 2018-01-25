!  +-LAMBDASOLVE----------------------------------------------------------+
!  | tasks: Compute the increments of the system                          |
!  |                  min |f + J*p|                                       |
!  |                s.t.  |D p  |  <= delta                               |
!  |        make use of the decomposition computed in DECOMPOSE           |
!  +----------------------------------------------------------------------+

!> \routine lambdasolve
!> \brief compute the increments of the linear constraint subproblem
!!  by using the decomposed Jacobian and Givens rotations.
!> \return DELTAP, lvlambda
!> \author Elias Roeger

subroutine lambdasolve(Ja, f, GRHS, GLOBALHELP, LMDIAGONAL,DELTAP, IPC, &
                         lvdelta, lvsigma, nbedtotal, nglob, dimm, lvlambda)

  implicit none
  integer :: dimm, & !< dimm = summe (nmess)
            nbedtotal, &
            nglob, &
            IPC(nglob), &
            i, j, k, counter
  real*8 :: Ja(dimm,nglob) ,&  
            f(dimm), &
            GRHS(nbedtotal), &
            GLOBALHELP(nbedtotal, nglob), &
            LMDIAGONAL(nglob), &
            DELTAP(nglob),&
            lvdelta, &
            lvsigma, &
            lowerbound, &
            upperbound, &
            phi0, &
            phi_dev, &
            norm, & !auxiliary variable
            eta, &
            alpha, &
            lvlambda


  real*8, allocatable :: JACtilde (:,:), & !< =Ja*D^-1 dim(dimm ,nglob)
                         JACtilde2 (:,:), &  !< Jactilde'
                         qalpha(:), & !hilfsvariable die D*p speichert
                         y(:), &
                         COPYGRHS(:), &                                    !raus
                         COPYGLOBALHELP(:,:), &                            !raus
                         COPY(:), & ! auxiliary variables for Permutations
                         COPYhelp(:), & !""
                         COPYLMD(:) ! ""

  allocate (JACtilde (dimm, nglob))
  allocate (JACtilde2(nglob, dimm))
  allocate (qalpha(nglob))
  allocate (y(nglob))
  allocate (COPYGRHS(nbedtotal))
  allocate (COPYGLOBALHELP(nbedtotal, nglob))
  allocate (COPY(nglob))
  allocate (COPYhelp(nglob))
  allocate (COPYLMD(nglob))

  alpha=1 ! alpha=alpha0


  COPYhelp = DELTAP
  do i = 1, nglob  
    if ( IPC(i) /= i ) then       !PIVOTING
      COPYhelp(IPC(i)) = DELTAP(i)
    endif
  enddo  
  DELTAP = COPYhelp
  



  do i = 1, nglob
    JACtilde(1:dimm, i) = Ja(1:dimm, i)* (1/LMDIAGONAL(i))
  enddo

  JACtilde2 = transpose(JACtilde)
  upperbound = sqrt(sum(matmul(JACtilde2, f)**2)) / lvdelta


  do i=1,nglob
    qalpha(i) = LMDIAGONAL(i) * DELTAP(i)
  enddo


 !calculate phi(0) and phi_dev(0)


  phi0 = sqrt(sum(qalpha**2)) - lvdelta


  norm = sqrt(sum(qalpha**2))

  !Copy <-E' * D * q /  norm(qalpha)
  COPY = LMDIAGONAL * qalpha

  COPYhelp = COPY
  do i = 1, nglob  
    if ( IPC(i) /= i ) then       !PIVOTING
      COPYhelp(i) = COPY(IPC(i))
    endif
  enddo  
  COPY = COPYhelp

  COPY = COPY/norm
 

  do i = 1, nglob
    y(i) = COPY(i)
    do j = 1,i-1
      y(i) = y(i) - GLOBALHELP(j,i)*y(j)
    enddo
    y(i)= y(i) / GLOBALHELP(i,i)
  enddo


  phi_dev= - norm * (sum(y**2))

! calc lower bound

  if (abs(phi_dev) >= 1.d+10)  then
    lowerbound = 0
  else 
    lowerbound = - phi0 / phi_dev
  endif

  counter = 0
  if (phi0 > 0) then
    do k =1,20 
  !Algo 5.5
  !a)     

      if (alpha < lowerbound .or. alpha > upperbound ) then
        alpha=max( 0.001*upperbound, sqrt(lowerbound*upperbound))
      endif 


  !b)
    !now calculate p_alpha

      COPYhelp=LMDIAGONAL

      do i = 1,nglob  
        if ( IPC(i) /= i ) then       !PIVOTING
          COPYhelp(i) = LMDIAGONAL(IPC(i))
        endif
      enddo  
      COPYLMD = COPYhelp


      GLOBALHELP (nbedtotal - nglob:nglob ,1:nglob) = 0.d0
      do i = 1,nglob
        GLOBALHELP (nbedtotal - nglob +i ,i) = sqrt(alpha) * COPYLMD(i)   !COPYLMD = E'*D*E, eigentlich nur E*D, oder doch richtig weil D nur vektor
      enddo


      COPYGLOBALHELP = GLOBALHELP

      COPYGRHS = GRHS

      call givens(COPYGLOBALHELP, COPYGRHS, nbedtotal, nglob)
    !weiter gehts mit löse COPYGLOBALHELP *y = COPYGRHS und p = - E * y

      !back substitution wieder nach Goĺub, (in kleinem Programm getestet, läuft)

      do j = nglob ,2 , -1
        COPYGRHS(j) = (COPYGRHS(j) / COPYGLOBALHELP(j,j))
        COPYGRHS(1 : j-1) = COPYGRHS(1 : j-1) - COPYGRHS(j)* COPYGLOBALHELP(1 : j-1, j)
      enddo 
      COPYGRHS(1) = COPYGRHS(1) / COPYGLOBALHELP(1,1)
      
      COPYhelp = COPYGRHS(1:nglob)

      do i = 1, nglob
        COPYHELP(IPC(i))=COPYGRHS(i)
      enddo 
  !! Also E'*vector realisiert mit hilfvektor(i) = vector(IPC(i)) und E*vector mit hilfvektor(IPC(I)) = vector(i)

      COPYGRHS(1:nglob)=COPYhelp

      do i = 1,nglob
        qalpha(i) = LMDIAGONAL(i) *(- COPYGRHS(i)) 
      enddo 

      norm = sqrt(sum(qalpha**2))
      phi0 = norm - lvdelta



      !Copy <-E' * D * q /  norm(qalpha)
      COPY = LMDIAGONAL * qalpha   
      COPYhelp = COPY
      do i = 1,nglob     
        if ( IPC(i) /= i ) then       !PIVOTING
          COPYhelp(i)      = COPY(IPC(i))
        endif
      enddo  
      COPY = COPYhelp /norm

      do i = 1, nglob
        y(i) = COPY(i)
        do j = 1,i-1
          y(i) = y(i) - COPYGLOBALHELP(j,i)*y(j)
        enddo
        y(i)= y(i) / COPYGLOBALHELP(i,i)
      enddo



      phi_dev = - norm * (sum(y**2))


      if (phi0 <= 0) then    !(<=> if phi(alpha) < 0)
        upperbound = alpha
      endif

      lowerbound = max(lowerbound , alpha - phi0 / phi_dev) 

  !Teste Abbruchbedingung

      if (abs(phi0) <= lvsigma * lvdelta) then
        exit
      endif

  !c)

      alpha = alpha - (phi0 + lvdelta)/lvdelta * (phi0/phi_dev)
      counter = counter+1

    enddo 

    DELTAP = COPYGRHS(1 : nglob)

    lvlambda = alpha 
  else
    lvlambda=0
  endif 

  deallocate (COPYLMD)
  deallocate (COPYhelp)
  deallocate (COPY)
  deallocate (COPYGLOBALHELP)
  deallocate (COPYGRHS)
  deallocate (y)
  deallocate (qalpha)
  deallocate (JACtilde2)
  deallocate (JACtilde)


end subroutine lambdasolve
