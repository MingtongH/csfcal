program functiontests
      use prep
      use detgen
      use projection
      use gramschmidt, only: prtmtx, mgs, mult, orth

      implicit none


      integer(i16b) :: det(2)
      integer:: i, j, n, m, k
      real(rk), allocatable ::z,  A(:, :), Q(:, :), R(:, :)
      !Lz_det Szt2_det tested, correct.
      !det(1) = 0
      !det(2) = 0
      !do i = 1, 10
      !    det(1) = det(1) + 1
      !    do j = 1, 10
      !        det(2) = det(2) + 1
      !        write(*, '(2b16, 2I8)') det(1:2), Lz_det(det), Szt2_det(det)
      !    enddo
      !enddo
      n = 4
      m = 10
      allocate(A(n, n), Q(n, n), R(n, n))
      A(1, :) = (/1., 2., 3., 4./)
      A(2, :) = (/5.,6.,7.,8./)
      A(3, :) = (/0.,9.,10.,11./)
      A(4, :) = (/0.,0.,12.,13./)
      n = 4
      m = 10
      call orth(A, Q, R, m, n)





end program functiontests

              
