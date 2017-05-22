program test0320_Y2Z
      use checkcsfs, only: sortBasisCoefTable_removeDups
      use convertharmonics
      implicit none

      integer(i16b) ::detplus, detminus
      integer(i16b), allocatable :: zbasislist(:, :), indet(:, :)
      real(rk), allocatable :: zcoefs(:, :), zcoeftable(:, :), realcoefs(:, :)
      integer :: pos, nzbasis, i, ncsf, nybasis
      !allocate(indet(2))
      !allocate(zbasislist(10, 2))
      !allocate(zcoefs(10, 2))

      !pos =13
      !write(*, *) sign_m(pos, detplus, detminus)
      !sign_m tested all correct
     
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Need to change neact in convertharmonics
      !1e tests
      !indet(1:2) = (/21_i16b, 0_i16b/)!Correct
      !indet(1:2) = (/0_i16b, 21_i16b/)!Correct
      !indet(1:2) = (/0_i16b, 512_i16b/)!Correct for 1 e from 1 - 512
      !indet(1:2) = (/512_i16b, 0_i16b/)!Correct

      !2e tests
      !indet(1:2) = (/3_i16b, 0_i16b/)!Correct all m=0  
      !indet(1:2) = (/0_i16b, 3_i16b/)! Correct
      !indet(1:2) = (/43_i16b, 7_i16b/)
      !indet(1:2) = (/0_i16b, 20_i16b/)
      !call Y2Z_det(indet, 1._rk,  zbasislist, zcoefs, 0)
      !write(*, *) sign_ordets(indet(1), indet(2))
      ! indet(1:2) = (/43_i16b, 7_i16b/)
      

      !test0309
      ! indet(1:2) = (/43_i16b, 7_i16b/)
      !test0316 Y2Z_singledet
      ! indet(1:2) = (/31_i16b, 0_i16b/) !good
      !indet(1:2) = (/52_i16b, 2_i16b/) !good
      !indet(1) =2_i16b**13 + 2_i16b**11+2_i16b**10 + 2_i16b**9
      !indet(2) = 4_i16b !good
      ncsf = 2
      nybasis = 2
      allocate(indet(nybasis, 2))
      indet(1, 1:2) = (/5122_i16b, 2_i16b/)
      indet(2, 1:2) = (/8706_i16b, 2_i16b/)
      nzbasis = 0 
      allocate(realcoefs(nybasis, ncsf))
      realcoefs(1:2, 1) = (/ 0.89442719099991896, -1.7888543819998302/)
      realcoefs(1:2, 2) = (/ -0.89442719099991896, 1.7888543819998302/)
      
      call Y2Z_appendtable(indet, realcoefs, nybasis, ncsf, zbasislist, zcoeftable, nzbasis)
      !call Y2Zcsf_append1row(indet, realcoefs, zbasislist, zcoeftable, nzbasis, ncsf)
     ! call Y2Z_singledet(indet, 3._rk, zbasislist, zcoefs, nzbasis) 
     ! call sortBasisCoefTable_removeDups(zbasislist, zcoefs, nzbasis, 2)

      !write(*, *) popcnt(indet(1)), popcnt(indet(2)), popcnt(13_i16b)
      !write(*, *) poppar(indet(1)), poppar(indet(2)), poppar(13_i16b)
      open(UNIT=1, FILE='Zcsfs_test0320')
      write(1, *) nzbasis, 'nzbasis'
      write(1, *) ncsf, 'ncsf'
      do i = 1, nzbasis
          write(1, '(2B16)') zbasislist(i, 1:2) 
          write(1, *) zcoeftable(i, 1:2*ncsf)
      enddo

      deallocate(zbasislist)
      deallocate(zcoeftable)

end program test0320_Y2Z
