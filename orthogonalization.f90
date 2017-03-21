module orthogonalization
      use, intrinsic :: iso_fortran_env, only: rk => real64
      use prep, only: equals0
      implicit none

      contains
      subroutine gs_modify(A, n, p)
          !A: n x p  n is nbasis p is ncsf
          !Q: n x p
          !R: p x p
          !A = QR

          !p could change, only output non-zero unit vectors in E
          !directly modifies table and p
          integer, intent(in) :: n
          real(rk), allocatable :: A(:, :),  Q(:, :), R(:, :), E(:, :)
          integer :: k, i, j,  p
          if(.not.allocated(Q)) then
              allocate(Q(n, p))
          endif
          if(.not.allocated(R)) then 
              allocate(R(p, p))
          endif


          write(*, *) '========================= Orthogonalization ================================'
          write(*, *) 'Input matrix A'
          do i = 1, n
              write(*, *) A(i, 1:p)
          enddo

          j = 1
          do k = 1, p
              Q(1:n, k) = A(1:n, k)
              if(k.ne.1) then
                  do i = 1, k-1
                      R(i, k) = dot_product(Q(1:n, i), Q(1:n, k))
                  enddo
                  !R(1:k-1, k) = matmul(transpose(Q(:, k-1:k-1)),Q(:, k))
                  Q(1:n, k) = Q(1:n, k) - matmul(Q(1:n, 1:k-1), R(1:k-1, k))
              endif
              R(k, k) = norm2(Q(1:n, k))
              if(equals0(R(k, k))) then
                  do i = 1, n
                      Q(i, k) = 0
                  enddo
              else
                  Q(1:n, k) = Q(1:n, k)/R(k, k)
                  A(1:n, j) = Q(1:n, k)!j<=k
                  j = j + 1
              endif
          enddo

          write(*, *) 'Matrix Q'
          do i = 1, n
              write(*, *) Q(i, 1:p)
          enddo
          write(*, *) 'Matrix R'
          do i = 1, p
              write(*, *) R(i, 1:p)
          enddo

          p = j - 1
          write(*, *) 'Output matrix A'
          do i = 1, n
              write(*, *) A(i, 1:p)
          enddo
          write(*, *) 'Output p =', p
          
          deallocate(Q)
          deallocate(R)

          
        
      end subroutine gs_modify

    
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
              Q(1:n, k) = A(1:n, k)
              if(k.ne.1) then
                  do i = 1, k-1
                      R(i, k) = dot_product(Q(1:n, i), Q(1:n, k))
                  enddo
                  !R(1:k-1, k) = matmul(transpose(Q(:, k-1:k-1)),Q(:, k))
                  Q(1:n, k) = Q(1:n, k) - matmul(Q(1:n, 1:k-1), R(1:k-1, k))
              endif
              R(k, k) = norm2(Q(1:n, k))
              if(equals0(R(k, k))) then
                  do i = 1, n
                      Q(i, k) = 0
                  enddo
              else
                  Q(1:n, k) = Q(1:n, k)/R(k, k)
                  E(1:n, j) = Q(1:n, k)
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
