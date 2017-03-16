program test0202
      use convertharmonics
      implicit none

      integer(i16b) ::detplus, detminus
      integer(i16b), allocatable :: zbasislist(:, :), indet(:)
      real(rk), allocatable :: zcoefs(:, :)
      real(rk) :: zoefs
      integer :: pos, nzbasis
      allocate(indet(2))
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
      indet(1) =2_i16b**13 + 2_i16b**11+2_i16b**10 + 2_i16b**9
      indet(2) = 4_i16b !good
      nzbasis = 0 
      call Y2Z_singledet(indet, 1._rk, zbasislist, zcoefs, nzbasis) 
      !write(*, *) popcnt(indet(1)), popcnt(indet(2)), popcnt(13_i16b)
      !write(*, *) poppar(indet(1)), poppar(indet(2)), poppar(13_i16b)

      deallocate(zbasislist)
      deallocate(zcoefs)

end program test0202
