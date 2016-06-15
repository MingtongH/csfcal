program test
      use prep
      use detgen

      implicit none
      !integer :: Lz_test
      !integer(i16b) :: det(1, 2)
      !det(1, 1) = 5
      !det(1, 2) = 3
      !Lz_test = Lz(det)
      integer, allocatable :: configs (:, :)
      integer(i16b) :: det
      
      allocate(configs(3, 3))
      configs(1, 1) = 2
      configs(1, 2) = 0
      configs(1, 3) = 1
      configs(2, 1) = 2
      configs(2, 2) = 1
      configs(2, 3) = 1
      configs(3, 1) = 3
      configs(3, 2) = 2
      configs(3, 3) = 1
      !write(*, *) isHalfFull(configs(1, 1:3))
      det = 0
      call assign_shell(det, 0, 0, 0, 0, 0, 0, configs, 0, 3)


end program test
