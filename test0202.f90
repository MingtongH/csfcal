program test0202
      use convertharmonics
      implicit none

      integer(i16b) :: indet(2), detplus, detminus
      integer(i16b), allocatable :: zbasislist(:, :)
      real(rk), allocatable :: zcoefs(:, :)
      real(rk) :: zoefs
      integer :: pos
      allocate(zbasislist(10, 2))
      allocate(zcoefs(10, 2))

      !pos =13
      !write(*, *) sign_m(pos, detplus, detminus)
      !sign_m tested all correct
      
      indet(1:2) = (/513_i16b, 0_i16b/)
      call Y2Z_det(indet, zbasislist, zcoefs)

      deallocate(zbasislist)
      deallocate(zcoefs)

end program test0202
