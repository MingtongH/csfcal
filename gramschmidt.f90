module gramschimidt
      use prep, only : i16b
      use, intrinsic :: iso_fortran_env, only: rk => real64
      use projection, only : REAL_MIN

      implicit none
      contains


      subroutine multip(A, B, n1, m, n2, C)
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


      end subroutine multip



