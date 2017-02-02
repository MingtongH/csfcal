module convertharmonics
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Step 1: create new basis list, with 2 coef columns for each csf
      !Step 2: calculate coefs for each csf
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 use, intrinsic :: iso_fortran_env, only: rk => real64
 use prep, only: i16b, neact
 use detgen, only: locate_det, delocate
    implicit none
    contains

     subroutine init_Ztable(ynbasis, ncsf, zbasis, zcoeftable, znbasis)
          integer, intent(in) :: ynbasis, ncsf
          integer(i16b), allocatable :: zbasis(:, :)
          real(rk), allocatable :: zcoeftable(:, :)
          integer :: znbasis, nrows
          znbasis = 0
          nrows = 2**neact * ynbasis
          if(.not.allocated(zbasis)) then
              allocate(zbasis(nrows, 2))
          endif
          if(.not.allocated(zcoeftable)) then
              allocate(zcoeftable(nrows, ncsf*2))
          endif
      end subroutine init_Ztable

      subroutine Y2Z_det(det, zbasislist, zcoefs)
          integer(i16b), intent(in) :: det(:)
          integer(i16b) :: zbasislist(:, :), tpup, tpdn, tpplus, tpminus
          real(rk) :: zcoefs(:, :)
          integer(i16b) :: singledets(:)
          real(rk) :: singlecoefs(:)
          integer :: pos, i, j,msign
          !zbasislist and zcoefs have to be allocated already
         
    
          tpup = det(1)
          tpdn = det(2)
          allocate(singledets(2*neact, 2))
          allocate(singlecoefs(2*neact, 2))
          i = 0
          j = neact

          do while(tpup.ne.0)
            i = i + 1
            j = j + 1
            pos = trailz(tpup)

            msign = sign_m(pos, tpplus, tpminus)
            singledets(i, 1) = tpplus
            singledets(i, 2) = 0_i16b
            singledets(j, 1) = tpminus
            singledets(j, 2) = 0_i16b
            singlecoefts(i, 1:2) =(/1.,0./)
            singlecoefts(j, 1:2) = (/0., msign/)

            write(*, *) i, singledets(i, 1:2), singlecoefts(i, 1:2)
            write(*, *) j, singledets(j, 1:2), singlecoefts(j, 1:2)
            tpup = ibclr(tpup, pos)
          enddo
           
          do while(tpdn.ne.0)
            i = i + 1
            j = j + 1
            pos = trailz(tpdn)
 
             msign = sign_m(pos, tpplus, tpminus)
             singledets(i, 2) = tpplus
             singledets(i, 1) = 0_i16b
             singledets(j, 2) = tpminus
             singledets(j, 1) = 0_i16b
             singlecoefs(i, 1:2) = (/1.,0./)
             singlecoefs(j, 1:2) = (/0., msign/)
             write(*, *) i, singledets(i, 1:2), singlecoefts(i, 1:2)
             write(*, *) j, singledets(j, 1:2), singlecoefts(j, 1:2)

             tpdn = ibclr(tpdn, pos)
          enddo

          !Now fill zbasislist and zcoefs
          nbasis = 2**neact
          do i = 1, nbasis/2
              zbasislist(i, 1:2) = singledets(1, 1:2)
              zcoefs(i, 1:2) = singlecoefs(1, 1:2)
              j = i + nbasis/2
              zbasislist(j, 1:2) = singledets(neact + 1, 1:2) 
              zcoefs(j, 1:2) = singlecoefs(neact + 1, 1:2)
          enddo
          do k = 2, neact
              i = 1
              do while(i.le.nbasis)
                  !fill in det plus
                  do j = 1, 2**(neact - k)
                       zbasislist(i, 1) = ior(zbasislist(i, 1), singledets(k, 1)) 
                       zbasislist(i, 2) = ior(zbasislist(i, 2), singledets(k, 2))
                       call product_cmplx(zoefs(i, 1:2), singlecoefs(k, 1:2), zoefs(i, 1:2))
                       i = i + 1
                  end do
                  !fill in det minus
                  do j = 1, 2**(neact - k)
                       zbasislist(i, 1) = ior(zbasislist(i, 1), singledets(k+neact, 1))
                       zbasislist(i, 2) = ior(zbasislist(i, 2), singledets(k+neact, 2))
                       call product_cmplx(zoefs(i, 1:2), singlecoefs(k+neact, 1:2), zoefs(i, 1:2))
                       i = i + 1
                  enddo
               enddo!end while
          enddo
          do i = 1, nbasis
              write(*,*) zbasislist(i, 1:2), zcoefs(i, 1:2)
          deallocate(singlecoefs)
          deallocate(zoefs)
      end subroutine Y2Z_det


      subroutine product_cmplx(c1, c2, nout)
          real(rk) :: c1(:), c2(:)
          real(rk) :: n1(2), n2(2), nout(:)
          n1(1:2) = c1(1:2)
          n2(1:2) = c2(1:2)
          nout(1) = n1(1)*n2(1) - n1(2)*n2(2)
          nout(2) = n1(1)*n2(2)) + n1(2)*n2(1)
      end subroutine product_cmplx
          
      integer function sign_m(pos, detplus, detminus)
          ! for a single det not a pair
          integer, intent(in) :: pos
          integer(i16b) :: detplus, detminus
          integer :: pos2
          call delocate(pos, n, l, m)
          if (m.eq.0)
              detplus = ibset(0_i16b, pos)
              detminus = ibset(0_i16b, pos)
          elseif(mout.gt.0)
              sign_m = 1
              detplus = ibset(0_i16b, pos)
              pos2 =  locate_det(n, l, -m)
              detminus = ibset(0_i16b, pos2)
          else
              sign_m = -1
              detminus = ibset(0_i16b, pos)
              pos2 = locate_det(n, l, -m)
              detplus = ibset(0_i16b, pos2)
          endif
          write(*, '(1I4, 2b16)') pos, detplus, detminus
      end function sign_m



end module convertharmonics
