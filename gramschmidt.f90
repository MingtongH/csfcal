module gramschmidt
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
!Original program only accepts square matrices

      use prep, only : i16b
      use, intrinsic :: iso_fortran_env, only: rk => real64
      use projection, only : REAL_MIN

      implicit none
      contains

      subroutine orth(a, q, r, m, n)
          real(rk) :: a(:, :), q(:, :), r(:, :), z
          integer, intent(in) :: m, n
          integer :: i, k
          do k = 1, m
              print *,' Matrix A'
              call prtmtx(n,a)
              z = a(n,n)
              do i=1,n
                a(i,i)= a(i,i) - z
              end do
              call mgs(a, q, r, n, n)
              print *,' Matrix Q'
              call prtmtx(n,q)
              print *,' Matrix R'
              call prtmtx(n,r)
              call mult(r,q, n, n, n, a)
              do i=1,n
                 a(i,i)=a(i,i) + z
              end do
          enddo
      end subroutine orth
      
      subroutine prtmtx(n,a)
          real(rk), intent(in):: a(:,:)
          integer, intent(in) :: n

          integer :: i, j
          do i=1,n
              write(*, '(4f10.6)') (a(i,j),j=1,n)
          enddo
      end subroutine

   

      subroutine mgs(a, q, t, m, n)
          !a(m, n), q(m, n), t(n, n)
          real(rk) :: a(:, :), q(:, :), t(:, :), z
          integer, intent(in) :: m, n
          integer :: i, j, k

          !if(size(a, 1).ne.m.OR.size(a, 2).ne.n.OR.size(q, 1).ne.m.OR.&
          !    &size(q, 2).ne.n.OR.size(t, 1).ne.n.OR.size(t, 2).ne.n) then
          !  write(*, *) 'matrix size not right in subroutine mgs'
          !  return
          !endif

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
                      z = z + q(i, j) * q(i, k)
                  end do
                  t(k, j) = z
                  do i = 1, m
                      q(i, j) = q(i, j) - t(k, j) * q(i, k)
                  enddo
              enddo
          enddo
      end subroutine mgs


      subroutine mult(A, B, n1, m, n2, C)
      ! Matrix product C(n1, n2) = A(n1, m)*B(m, n2)
      integer, intent(in) :: n1, m, n2
      real(rk), intent(in) :: A(:, :), B(:, :)
      real(rk) :: C(:, :), tp
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

end module gramschmidt
