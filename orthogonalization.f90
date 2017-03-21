module orthogonalization
      use, intrinsic :: iso_fortran_env, only: rk => real64
      use prep, only: equals0
      implicit none

      contains

    
      subroutine gs(A, Q, R, E, n,p, ncsf)
          !A: n x p
          !Q: n x p
          !R: p x p
          !A = QR

          !p could change, only output non-zero unit vectors in E
          integer, intent(in) :: n, p
          real(rk), intent(in) :: A(:, :)
          real(rk), allocatable :: Q(:, :), R(:, :), E(:, :)
          integer :: k, i, j, ncsf
          if(.not.allocated(Q)) then
              allocate(Q(n, p))
          endif
          if(.not.allocated(R)) then 
              allocate(R(p, p))
          endif

          if(.not.allocated(E)) then
              allocate(E(n, p))
          endif

          write(*, *) '========================= Orthogonalization ================================'
          j = 1
          do k = 1, p
              Q(:, k) = A(:, k)
              if(k.ne.1) then
                  do i = 1, k-1
                      R(i, k) = dot_product(Q(:, i), Q(:, k))
                  enddo
                  !R(1:k-1, k) = matmul(transpose(Q(:, k-1:k-1)),Q(:, k))
                  Q(:, k) = Q(:, k) - matmul(Q(:, 1:k-1), R(1:k-1, k))
              endif
              R(k, k) = norm2(Q(:, k))
              if(equals0(R(k, k))) then
                  do i = 1, n
                      Q(i, k) = 0
                  enddo
              else
                  Q(:, k) = Q(:, k)/R(k, k)
                  E(:, j) = Q(:, k)
                  j = j + 1
              endif
          enddo
          ncsf = j - 1
          write(*, *) 'Input matrix A'
          write(*, *) A
          write(*, *) 'Matrix Q'
          write(*, *) Q
          write(*, *) 'Matrix R'
          write(*, *) R
          write(*, *) 'Matrix E'
          write(*, *) E

        
      end subroutine gs
end module orthogonalization
