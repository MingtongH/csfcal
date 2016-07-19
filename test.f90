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
      integer(i16b) :: det, subdet, detup, detdn
      integer :: nocc, count_tot, i, n, l, m
      integer, allocatable :: occ_up(:), occ_dn(:), eposup(:), eposdn(:), tmpepos
      integer(i16b), allocatable :: detlist(:, :)
      real(rk) :: coef
      logical :: iszero
      logical, allocatable :: iszerolist(:), outiszerolist(:), iszerolist1(:),&
          & iszerolist2(:), iszerolist3(:)
      integer(i16b), allocatable :: basislist(:, :), outbasislist(:, :), basislist1(:, :),&
          & basislist2(:, :), basislist3(:, :), allbasis(:, :)
      real(rk), allocatable :: coeflist(:), outcoeflist(:), coeflist1(:), coeflist2(:), coeflist3(:), coeftable(:, :)
      integer, allocatable :: eposlist(:, :, :), outeposlist(:, :, :), eposlist1(:, :, :), eposlist2(:, :, :), eposlist3(:, :, :)
      integer :: num, outnum, num1, num2, num3, ndets, ncsf, j

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
      
      !call eposinit(eposup, eposdn, detlist(9, 1), detlist(9, 2))
      !write(*, '(15I4)') eposup(1:15)
      !write(*, '(15I4)') eposdn(1:15)
      
      write(*, '("************series of operators *********************")')
      detup = 9
      detdn = 4
      call eposinit(eposup, eposdn, detup, detdn)
      write(*, '(2B16)') detup, detdn
      write(*, '(5I4)') eposup(1:5)
      write(*, '(5I4)') eposdn(1:5)
      call getsign(detup, detdn, eposup, eposdn, coef)

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
      !TODO Attention! Initialize these
      coef = 1
      iszero = .false.
      num = 0

      call Lminus_single(detup, detdn, eposup, eposdn, coef, iszero, &
          & basislist, coeflist, eposlist, iszerolist, num)
     !call Lplus_single(detlist(9, 1), detlist(9, 2), eposup, eposdn, coef, iszero, &
     !     & basislist, coeflist, eposlist, iszerolist, num)
     !call Sminus_single(detlist(9, 1), detlist(9, 2), eposup, eposdn, coef, iszero, &
     !    & basislist, coeflist, eposlist, iszerolist, num)
     !call Splus_single(detlist(9, 1), detlist(9, 2), eposup, eposdn, coef, iszero, &
     !    & basislist, coeflist, eposlist, iszerolist, num)
     !TODO Attention! Initialize this
      outnum = 0
      num1 = 0
      num2 = 0
      num3 = 0
      call Lplus_multiple(basislist, coeflist, eposlist, iszerolist, num, &
           & basislist1, coeflist1, eposlist1, iszerolist1, num1)
      call Sminus_multiple(basislist1, coeflist1, eposlist1, iszerolist1, num1, &
           & basislist2, coeflist2, eposlist2, iszerolist2, num2)
      call Splus_multiple(basislist2, coeflist2, eposlist2, iszerolist2, num2, &
           & basislist3, coeflist3, eposlist3, iszerolist3, num3)
      call getallsigns(basislist3, coeflist3, eposlist3, iszerolist3, num3)
      do i = 1, num3
        write(*, '(1I3, 2b16, 1f18.15)') i,  basislist3(i, :), coeflist3(i)
        write(*,'(10I3)') eposlist3(i, 1, 1:10)
        write(*, '(10I3)') eposlist3(i, 2, 1:10)
      enddo
      ndets = 0
      ncsf = 0
      call collect_csf(basislist3, coeflist3, iszerolist3, num3, &
          & allbasis, coeftable, ndets, ncsf)
      call collect_csf(basislist2, coeflist2, iszerolist2, num2, &
          & allbasis, coeftable, ndets, ncsf)

      do i = 1, ndets
          write(*, '(2b16, 2F15.9)') allbasis(i, 1:2), coeftable(i, 1:2)
      enddo
      call normalizetable(coeftable, ndets, ncsf)
     do i = 1, ndets
          write(*, '(2b16, 2F15.9)') allbasis(i, 1:2), coeftable(i, 1:2)
      enddo

       !call getsign(basislist3(2, 1), basislist3(2, 2), eposlist3(2, 1, :), eposlist3(2, 2, :), coef)
      !eposup(1:5) = (/4, 2, 1, 3, 5/)  
      !outnum = countswps(eposup, 5) 


      !call Sminus_multiple(basislist, coeflist, eposlist, iszerolist, num, &
      !    & outbasislist, outcoeflist, outeposlist, outiszerolist, outnum)
      !call Splus_multiple(basislist, coeflist, eposlist, iszerolist, num, &
      !    & outbasislist, outcoeflist, outeposlist, outiszerolist, outnum)
      
      !call Lplus_multiple(basislist, coeflist, eposlist, iszerolist, num, &
      !    & outbasislist, outcoeflist, outeposlist, outiszerolist, outnum)
      !call Lminus_multiple(basislist, coeflist, eposlist, iszerolist, num, &
      !    & outbasislist, outcoeflist, outeposlist, outiszerolist, outnum)

      
      
      
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
