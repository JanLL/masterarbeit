!  +-Decompose-----------------------------------------------------------+
!  | tasks: Decompose the Matrix (R 0 D)'                                |
!  |        Use Givens Rotations to get rid of the D Block               |
!  +---------------------------------------------------------------------+

!> \routine givens
!> \brief use Givens QR decomposition to decompose
!! the Matrix (R 0 D)'
!> \author Elias Roeger


subroutine givens(Matrix, vectorb, m, n) 

    implicit none

    Integer i, j    !Counters
    Integer m, n    ! Matrix is m x n 
    Real*8 ,Dimension(m,n) :: Matrix
    Real*8 ,Dimension(m) :: vectorb
    Real*8 ,Dimension(2,2) :: Rotationmatrix      
    Real*8 :: c ,s  
    Real*8 :: help1, help2

    do j = 1,n
        do i = m , j+1 , -1
            help1 = Matrix(i-1, j)
            help2 = Matrix(i,j)
            call givensvektor( help1 , help2, c, s)
            Rotationmatrix = reshape( (/c, s, -s, c/), (/2, 2/) )
            vectorb(i-1 :i) = matmul(Rotationmatrix, vectorb(i-1:i))
            Matrix(i-1 :i , j :n)= matmul(Rotationmatrix, Matrix(i-1:i, j:n))
        enddo
    enddo

end subroutine givens

subroutine givensvektor ( a, b, c, s) 
    implicit none
    real*8 :: a, b
    real*8 ,intent(out) :: c ,s      !we need c and s in givens
    real*8 :: tau               !helping variable
    if (b == 0 ) then
        c=1
        s=0
    else
        if (abs(b) > abs(a)) then
            tau = - a/b
            s = 1/(sqrt(1+ tau**2))
            c = s*tau
        else
            tau = -b/a
            c = 1/(sqrt(1+ tau**2))
            s = c*tau
        endif
    endif

end subroutine givensvektor 


