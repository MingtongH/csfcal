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
      integer(i16b) :: det, subdet
      integer :: nocc
      integer, allocatable :: occ_up(:), occ_dn(:)
      allocate(configs(3, 3))
      configs(1, 1:3) = (/2, 0, 1/)
      configs(2, 1:3) = (/2, 1, 5/)
      configs(3, 1:3) = (/3, 2, 1/)
      !write(*, *) Lzmax(configs)
     
      call assign_shell(det, det, 0, 0, 0, 0, 0, 1, configs, 0)
      !write(*, *) isHalfFull(configs(1, 1:3))
      !allocate(occ_up(3)) 
      !allocate(occ_dn(3))
      !do det = 0, 8
      !  occ_up(1:3) = (/0,0,0/)
      !  occ_dn(1:3) = (/0,0,0/)
      !  nocc = count_occ_orbs(det)*2
      !  
      !  write(*, '("det = ", 1I5)') det
      !  write(*, '("nocc = ", 1I5)') nocc
      !  call get_occ_orbs(det, det, occ_up, occ_dn, nocc)
      !  write(*, '("nocc = ", 1I5)') nocc
      !  write(*, *) occ_up
      !  write(*, *) occ_dn
      !enddo
      
      !bit manipulation POS starts from idex 0
      !det = 2
      !write(*, *) BTEST(det, 0)
      !write(*, *) BTEST(det, 1)
      !write(*, *) BTEST(det, 2)
      
     ! det = 2
     ! write(*, *) trailz(det) 
     ! det = 26
     ! write(*, *) count_occ_orbs(det)

      !write(*, *) Lz_updn(0, 3, det)
      !det = 4
      !write(*, *) count_occ_orbs(det)

      !write(*, *) Lz_updn(0, 2, det)
      !det = 8
      !write(*, *) count_occ_orbs(det)
 
     !write(*, *) Lz_updn(0, 2, det)
     ! Lz_unit correct
     ! Lz_updn correct
     !write(*, *) locate_det(3, 2, -1)
     !det = 8
     !subdet = 7
     !call detmod_shell(3, 1, -1, subdet, det)
     !write(*, *) det
end program test
