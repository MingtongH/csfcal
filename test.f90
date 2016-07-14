program test
      use prep
      use detgen
      use projection

      implicit none
      !integer :: Lz_test
      !integer(i16b) :: det(1, 2)
      !det(1, 1) = 5
      !det(1, 2) = 3
      !Lz_test = Lz(det)
      integer, allocatable :: configs (:, :)
      integer(i16b) :: det, subdet
      integer :: nocc, count_tot, i, n, l, m
      integer, allocatable :: occ_up(:), occ_dn(:), eposup(:), eposdn(:), tmpepos
      integer(i16b), allocatable :: detlist(:, :)
      real(rk) :: coef
      logical :: iszero
      logical, allocatable :: iszerolist(:)
      integer(i16b), allocatable :: basislist(:, :)
      real(rk), allocatable :: coeflist(:)
      integer, allocatable :: eposlist(:, :, :)
      integer :: num

      allocate(configs(3, 3))
      configs(1, 1:3) = (/2, 0, 1/)
      configs(2, 1:3) = (/2, 1, 5/)
      configs(3, 1:3) = (/3, 2, 1/)
      !write(*, *) Lzmax(configs)
     count_tot = 0
     call assign_shell(det, det, 0, 0, 0, 0, 0, 1, configs, 0, count_tot, detlist)
      write(*, *) count_tot
      do i = 1, count_tot
        write(*, '(2B16)') detlist(i, 1:2)
      enddo
      
      call eposinit(eposup, eposdn, detlist(6, 1), detlist(6, 2))
      write(*, '(15I4)') eposup(1:15)
      write(*, '(15I4)') eposdn(1:15)
      !do i = 0, 15
      !  det = detlist(1, 2)
      !  write(*, '("i = ", 1I5)') i
      !  write(*, '(1B16)') det
      !  coef = 1
      !  iszero = .false.
      !  !call lpls(i, det, eposdn, coef)
      !  call lms(i, det, eposdn, coef, iszero)
      !  write(*, *) coef, iszero
      !  write(*, '(1B16)') det
      !end do
      coef = 1
      iszero = .false.
      num = 0

      call Lplus_single(detlist(6, 1), detlist(6, 2), eposup, eposdn, coef, iszero, &
          & basislist, coeflist, eposlist, iszerolist, num)
     ! do i = 0, DET_MAX_LENGTH
     !   call delocate(i, n, l, m)
     !   write(*, '(4I5)') i, n, l, m
     !   call pos2nlm(i, n, l, m)
     !   write(*, '(4I5)') i, n, l, m
     ! enddo


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
      
      !bit manipulation POS starts from idex 0, btest
      !det = 0
      !write(*, *) ibset(det, 0)
      !write(*, *) ibset(det, 1)
      !write(*, *) ibset(det, 2)
     ! det = 14
     ! det = ishft(det, -1)
     ! write(*, *) det
     ! det = ishft(det, -2)
     ! write(*, *) det
     ! det = 2
    !  write(*, *) trailz(det) 
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
     !det = 12
     !subdet = 4
     !call detmod_shell(3, 1, -1, subdet, det)
     !write(*, '(3B16)') 12, subdet, det
     !det = 20
     !write(*,*) detdisplay(det)
end program test
