module gramschimidt
!Modified based on
!     Numerical Analysis:
!     Mathematics of Scientific Computing
!     Third Edition
!     D.R. Kincaid & E.W. Cheney
!     Brooks/Cole Publ., 2002
!     Copyright (c) 1996
!     Section 5.3
!     Example of modified Gram-Schmidt algorithm
!     file: qrshif.f

      use prep, only : i16b
      use, intrinsic :: iso_fortran_env, only: rk => real64
      use projection, only : REAL_MIN

      implicit none
      contains

      subroutine mgs(a, q, t, m, n)
          !a(m, n), q(m, n), t(n, n)
          real(rk) :: a(:, :), q(:, :), t(:, :), z
          integer, intent(in) :: m, n
          integer :: i, j, k

          if(size(a, 1).ne.m.OR.size(a, 2).ne.n.OR.size(q, 1).ne.m.OR.&
              &size(q, 2).ne.n.OR.size(t, 1).ne.n.OR.size(t, 2).ne.n) then
            write(*, *) 'matrix size not right in subroutine mgs'
            return
          endif

          do j = 1, n
              do i = 1, n
                  t(i, j) = 0.
              enddo
              do i = 1, m
                  q(i, j) = a(i, j)
              enddo
          enddo

          do k = 1, n
              z = 0.
              do i = 1, m
                  z = z + q(i, k)**2
              end do
              t(k, k) = sqrt(z)
              do i = 1, m
                  q(i, k) = q(i, k) / t(k, k)
              end do
              do j = k+1, n
                  z = 0.
                  do i = 1, m
                      z = z + q(i, k)
                  end do
                  t(k, j) = z
                  do i = 1, m
                      q(i, j) = q(i, j) - t(k, j) * q(i, k)
                  enddo
              enddo
          enddo
      subroutine mgs


      subroutine mult(A, B, n1, m, n2, C)
      ! Matrix product C(n1, n2) = A(n1, m)*B(m, n2)
      integer, intent(in) :: n1, m, n2
      real(rk), intent(in) :: A(n1, m), B(m, n2)
      real(rk) :: C(n1, n2), tp
      integer :: i, j, k

        do i = 1, n1
            do j = 1, n2
                tp = 0.0
                do k = 1, m
                    tp = tp + A(i, k) * B(k, j)
                enddo
                C(i, j) = tp
            end do
        end do
      end subroutine mult



