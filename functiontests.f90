program functiontests
      use prep
      use detgen
      use projection
      implicit none

      integer(i16b) :: det(2)
      integer:: i, j
      !Lz_det Szt2_det tested, correct.
      det(1) = 0
      det(2) = 0
      do i = 1, 10
          det(1) = det(1) + 1
          do j = 1, 10
              det(2) = det(2) + 1
              write(*, '(2b16, 2I8)') det(1:2), Lz_det(det), Szt2_det(det)
          enddo
      enddo
end program functiontests

              
