module convertharmonics
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Step 1: create new basis list, with 2 coef columns for each csf
      !Step 2: calculate coefs for each csf
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 use, intrinsic :: iso_fortran_env, only: rk => real64
 use prep, only: i16b, ARRAY_SHORT_LENGTH
 !, neact TODO for testing purposed this is not imported, will need to afterwards
 use detgen, only: locate_det, delocate
 use checkcsfs, only: sortBasisCoefTable_removeDups
    implicit none
    contains

     subroutine init_Ztable(ynbasis, ncsf, zbasis, zcoeftable, znbasis)
          integer, intent(in) :: ynbasis, ncsf
          integer(i16b), allocatable :: zbasis(:, :)
          real(rk), allocatable:: zcoeftable(:, :)
          integer :: znbasis, nrows
          integer :: neact 

          znbasis = 0
          nrows = 2**neact * ynbasis
          if(.not.allocated(zbasis)) then
              allocate(zbasis(nrows, 2))
          endif
          if(.not.allocated(zcoeftable)) then
              allocate(zcoeftable(nrows, ncsf*2))
          endif
      end subroutine init_Ztable

      !subroutine Y2Z_det(det, nshell_1config, econfigs, coef, allzbasislist, allzcoefs, nzbasis)
      !    integer(i16b), intent(in) :: det(:)
      !    integer, intent(in) :: nshell_1config, econfigs(:, :)
      !    real(rk), intent(in):: coef
      !    integer(i16b)

      integer function sign_ordets(ldet, rdet)
          integer(i16b), intent(in) :: ldet, rdet !single det same spin, not a pair
          !ldet is the newcomer that should be bigger, to be inserted to rdet from the left side
          integer :: ltail, rhead
          !In later use, as a convention, should put det already process at rdet, 
          !but new single det to be incorporated at ldet
          integer(i16b) ::tpl, subbit

          write(*, '("Entering sign_ordets, inputs", 2B16)') ldet, rdet
          if(iand(ldet, rdet).ne.0) then 
              sign_ordets = 0
          else
              !position of the leading 1 in rdet
              rhead = leadz(0_i16b) - leadz(rdet) - 1 !position in bit starts from 0
              write(*, *) "rhead pos = ", rhead
              sign_ordets = 1
              tpl = ldet

              ltail = trailz(tpl)
              do while(ltail.lt.rhead.AND.tpl.ne.0)
                  subbit = ibits(rdet, ltail, rhead - ltail + 1)
                  write(*, '(1B16)') subbit
                  if(poppar(subbit).eq.1) then 
                      sign_ordets = - sign_ordets
                  endif
                  tpl = ibclr(tpl, ltail)
                  ltail = trailz(tpl)
              enddo
          endif

          write(*, '("sign_ordets = ", 1I6)') sign_ordets
      end function sign_ordets


      subroutine Y2Z_appendtable(ybasislist, ycsftable, ynbasis, ncsf, &
              &zbasislist, zcoeftable, nzbasis)
          integer(i16b), intent(in) :: ybasislist(:, :)
          real(rk), intent(in) :: ycsftable(:, :)
          integer, intent(in) :: ynbasis, ncsf
          integer(i16b), allocatable :: zbasislist(:, :)
          real(rk), allocatable :: zcoeftable(:, :)
          integer :: nzbasis, i
          write(*, *) '********************* Converting csftable Y2Z_appendtable *******************'
          do i = 1, ynbasis
              call Y2Zcsf_append1row(ybasislist(i, 1:2), ycsftable(i, 1:ncsf), &
                  &zbasislist, zcoeftable, nzbasis, ncsf)
          enddo

          call sortBasisCoefTable_removeDups(ybasislist, zcoeftable, nzbasis, 2*ncsf)

          write(*, *) '********************* Final number of rows in the table *********************'
          write(*, *) ncsf

      end subroutine Y2Z_appendtable

      subroutine Y2Zcsf_append1row(det, realcoefs, zbasislist, zcoeftable, nzbasis, ncsf)
          !nzbasis will be incremented here
          integer(i16b), intent(in) :: det(:)
          real(rk), intent(in) :: realcoefs(:)
          integer(i16b), allocatable :: tpzbasislist(:, :), &
              &zbasislist(:, :)!two columns, append to this list

          real(rk), allocatable :: tpzcoefs(:, :),&
              &zcoeftable(:, :)!2*ncsf columns, real imag for each csf
          integer :: nzbasis, i, j, tpnz
          integer, intent(in) :: ncsf

          write(*, *) '*******************  Converting 1 row of Ycsf to Z ***********************'
          write(*, *) 'Input:'
          write(*, '(2B16)') det(1:2)
          write(*, *) 'Number of ncsf = ', ncsf
          write(*, *) 'Number of rows in zbasislist before =', nzbasis
          call Y2Z_singledet(det, 1._rk, tpzbasislist, tpzcoefs, tpnz)
          call sortBasisCoefTable_removeDups(tpzbasislist, tpzcoefs, tpnz, 2)

          if(.not.allocated(zbasislist)) then
              allocate(zbasislist(ARRAY_SHORT_LENGTH, 2))
          endif
          if(.not.allocated(zcoeftable)) then
              allocate(zcoeftable(ARRAY_SHORT_LENGTH, 2*ncsf))
          endif

          write(*, *) 'Appending the following rows to the list and table'
          do i = 1, tpnz
              zbasislist(nzbasis + i, 1:2) = tpzbasislist(i, 1:2)
              write(*, *) '-------------------------------'
              write(*, '(2B16)') zbasislist(nzbasis + i, 1:2)
              write(*, *) '--------------------------------'

              do j = 1, ncsf
                  zcoeftable(nzbasis + i, 2*j - 1) = realcoefs(j) * tpzcoefs(i, 1) 
                  zcoeftable(nzbasis + i, 2*j) = realcoefs(j) * tpzcoefs(i, 2) 
                  write(*, *) zcoeftable(nzbasis + i, (2*j-1):(2*j))
              enddo
          enddo
          nzbasis = nzbasis + tpnz
          write(*, *) 'Now nzbasis =', nzbasis
          deallocate(tpzbasislist)
          deallocate(tpzcoefs)

      end subroutine Y2Zcsf_append1row


      subroutine Y2Z_singledet(det, realcoef, zbasislist, zcoefs, nzbasis)
          integer(i16b), intent(in) :: det(:)
          real(rk), intent(in) :: realcoef
          integer(i16b), allocatable:: zbasislist(:, :) ! Append to this list
          real(rk), allocatable:: zcoefs(:, :) !two columnsn real, imag
          integer :: nzbasis ! number of rows already in zbasislist
          integer(i16b), allocatable :: singledets(:, :) !four columns up dn for |m|, -|m| 
                                                         !             det = 0, m = 0
          real(rk), allocatable :: singlecoefs(:, :) ! four columns real, imag for each det
          integer(i16b) :: tpup, tpdn, tpplus, tpminus
          integer :: nelec, pos, pos2, msign, i, istart, tpsign, j, k, iend, tpspin, nsingleup

           write(*, *) ' '
           write(*, *) '******************* Entering Y2Z_singledet *****************'
           write(*, *) 'input det pair is'
           write(*, '(2B16)') det(1:2)
           tpup = det(1)
           tpdn = det(2)
           nelec = popcnt(tpup) + popcnt(tpdn)
           write(*, *) 'total number of electrons = ', nelec


           allocate(singledets(nelec, 4))
           allocate(singlecoefs(nelec, 4))
           write(*, *) '=================Generating single electron det list================='

           i = 1
           write(*, *) '------------------Generate single dets from det up-------------------'
           do while(tpup.ne.0)
               pos = trailz(tpup)
               call sign_m(pos, pos2, tpplus, tpminus, msign)
               !singledets(i, 1:4) = (/tpplus, 0_i16b, tpminus, 0_i16b/)
               !singlecoefs(i, 1:4) = (/1._rk, 0._rk, 0._rk, msign*1._rk/)
               singledets(i, 1) = tpplus
               singledets(i, 2) = 0_i16b
               singledets(i, 3) = tpminus
               singledets(i, 4) = 0_i16b
               singlecoefs(i, 1:3) = (/1._rk, 0._rk, 0._rk/)
               singlecoefs(i, 4) = msign*1._rk
               write(*, '(4B16)') singledets(i, 1:4)
               write(*, *) singlecoefs(i, 1:4)
               i = i + 1
               tpup = ibclr(tpup, pos) !not labeling m, -m coexistance, may need to later
           end do
           nsingleup = i - 1
           write(*, *) '------------------generate single dets for det dn---------------------'
           do while(tpdn.ne.0)
               pos = trailz(tpdn)
               call sign_m(pos, pos2, tpplus, tpminus, msign)
               !singledets(i, 1:4) = (/ 0_i16b,tpplus, 0_i16b, tpminus/)
               !singlecoefs(i, 1:4) = (/1._rk, 0._rk, 0._rk, msign*1._rk/)
               singledets(i, 1) = 0_i16b
               singledets(i, 2) = tpplus
               singledets(i, 3) = 0_i16b
               singledets(i, 4) = tpminus
               singlecoefs(i, 1:3) = (/1._rk, 0._rk, 0._rk/)
               singlecoefs(i, 4) = msign*1._rk
               write(*, '(4B16)') singledets(i, 1:4)
               write(*, *) singlecoefs(i, 1:4)
               i = i + 1
               tpdn = ibclr(tpdn, pos) !not labeling m, -m coexistance, may need to later
           end do
           !now i is index for next row
           nelec = i - 1
           write(*, *) 'Confirm total number of electrons = ', nelec

           !code above this line tested, all correct 
           if(.not.allocated(zbasislist)) then
               allocate(zbasislist(ARRAY_SHORT_LENGTH, 2))
               write(*, *) 'allocate zbasislist'
           endif

           if(.not.allocated(zcoefs)) then
               allocate(zcoefs(ARRAY_SHORT_LENGTH, 2))
               nzbasis = 0
               write(*, *) 'allocate zcoefs'
           endif
           istart = nzbasis + 1
           !zbasislist(istart, 1:2) = (/0_i16b, 0_i16b/)
           !zcoefs(istart, 1:2) = (/1._rk, 0._rk/)
           iend = istart !range to do ior
           i = istart !

           write(*, *) '======================= Generate Z basislist ============================'
           do j = 1, nelec ! index for processing singledets
            
               !write(*, *) 'Now in the  do loop, j = ', j
               !write(*, *) 'nzbasis =', nzbasis, 'istart=', istart, 'iend=', iend
               !if(j.eq.1) then 
                !   zbasislist(istart, 1) = singledets(j, 1)
                !   zbasislist(istart, 2) = singledets(j, 2)
                !   zcoefs(istart, 1) = singlecoefs(j, 1)
                !   zcoefs(istart, 2) = singlecoefs(j, 2)
                !   nzbasis = nzbasis + 1
                !   write(*, *) 'j = 1 directly copy det and coef -------------------------------'
                !   write(*, '(2B16)') zbasislist(i, 1:2)
                !   write(*, *) zcoefs(i, 1:2)
                 !  cycle
              ! endif
              

               write(*, *) 'j = ', j, '------------------------------------------------------'
               if(j.le.nsingleup) then 
                   tpspin = 1
               else
                   tpspin = 2
               endif
               write(*, *) 'tpspin =', tpspin
               !case m == 0  
              !=====================================================================================
                            
               if(singledets(j, 3).eq.0.AND.singledets(j, 4).eq.0) then

                      write(*, *) '>>>> det m = 0'
                      if(j.eq.1) then 
                           zbasislist(istart, 1) = singledets(j, 1)
                           zbasislist(istart, 2) = singledets(j, 2)
                           zcoefs(istart, 1) = singlecoefs(j, 1)
                           zcoefs(istart, 2) = singlecoefs(j, 2)
                           nzbasis = nzbasis + 1
                           !iend doesn't change still 1, nzbasis 0 to 1
                           write(*, *) 'j = 1 directly copy det &
                               & and coef '
                           write(*, '(2B16)') zbasislist(istart, 1:2)
                           write(*, *) zcoefs(istart, 1:2)
   
                           cycle
                      endif
           

                      write(*, *) 'Directly ior existing dets, coefs x 1 not to be changed'
                      do i = istart, iend
                         write(*, *) 'processing det No', i, 'in the appended list'
                         tpsign = sign_ordets(singledets(j, tpspin), zbasislist(i, tpspin))
                         zbasislist(i, tpspin) = ior(singledets(j, tpspin), zbasislist(i, tpspin))
                         !call product_cmplx(singlecoefs(j, 1:2), zcoefs(i, 1:2), &
                         !    &zcoefs(i, 1:2))
                         zcoefs(i, 1) = tpsign * zcoefs(i, 1)
                         zcoefs(i, 2) = tpsign * zcoefs(i, 2)
                         write(*, '(2B16)') zbasislist(i, 1:2)
                         write(*, *) zcoefs(i, 1:2)

                      enddo
                      !iend doesn't change, nzbasis doesn't change
                      write(*, *) ' '
                      write(*, *) 'Processed dets number', istart, iend
                      do i = istart, iend
                            write(*, '(2B16)') zbasislist(i, 1:2)
                            write(*, *) zcoefs(i, 1:2)
                      enddo

                !Note: for case det j already set, because the order of process, &
                !the plus det is always on top, so first det in zbasislist is always m>=0

            !cases m!=0
            !====================================================================================
                else 
                    write(*, *) '>>>> det m!=0'
                    if(j.eq.1) then 
                           zbasislist(istart, 1:2) = singledets(j, 1:2)
                           zbasislist(istart+1, 1:2) = singledets(j, 3:4)
                           zcoefs(istart, 1:2) = singlecoefs(j, 1:2)
                           zcoefs(istart+1, 1:2) = singlecoefs(j, 3:4)
                           iend = istart + 1
                           nzbasis = iend
                           write(*, *) 'j = 1 directly copy  2 dets &
                               & and coefs '
                           write(*, *) 'Added dets number', istart, iend
                            do i = istart, iend
                              write(*, '(2B16)') zbasislist(i, 1:2)
                              write(*, *) zcoefs(i, 1:2)
                            enddo
                           cycle
                      endif
           

          
                    tpsign = sign_ordets(singledets(j, tpspin), zbasislist(istart, tpspin))

                    !----------------------------------------------------------------------
                    !bit already set, m, -m coexist
                         
                    if(tpsign.eq.0) then
                          write(*, *) 'm, -m coexist has been previously set-----------------'
                          do i = istart, iend
                              !tpsign =  sign_ordets(singledets(j, 1), zbasislist(i, 1)) * &
                              !  &sign_ordets(singledets(j, 2), zbasislist(i, 2))
                              tpsign = sign_ordets(singledets(j, tpspin), zbasislist(i, tpspin))

                              if(tpsign.eq.0) then !or -|m|
                                  tpsign = sign_ordets(singledets(j, 2+tpspin), &
                                      &zbasislist(i, tpspin))
                                  zbasislist(i, 1) = ior(singledets(j, 3), zbasislist(i, 1))
                                  zbasislist(i, 2) = ior(singledets(j, 4), zbasislist(i, 2))
                                  call product_cmplx(singlecoefs(j, 3:4), zcoefs(i, 1:2), &
                                    &zcoefs(i, 1:2))
                                  zcoefs(i, 1) = tpsign * zcoefs(i, 1)
                                  zcoefs(i, 2) = tpsign * zcoefs(i, 2)

                              else !or |m|
                                  zbasislist(i, 1) = ior(singledets(j, 1), zbasislist(i, 1))
                                  zbasislist(i, 2) = ior(singledets(j, 2), zbasislist(i, 2))
                                  call product_cmplx(singlecoefs(j, 1:2), zcoefs(i, 1:2), &
                                      & zcoefs(i, 1:2))
                                  zcoefs(i, 1) = tpsign * zcoefs(i, 1)
                                  zcoefs(i, 2) = tpsign * zcoefs(i, 2)
                              endif

                              write(*, *) 'processed det No', i, 'in the appended list'
                              write(*, '(2B16)') zbasislist(i, 1:2)
                              write(*, *) zcoefs(i, 1:2)

                          enddo
                          !idend, nzbasis unchanged
                          write(*, *) ' '
                          write(*, *) 'Processed dets number', istart, iend
                          do i = istart, iend
                            write(*, '(2B16)') zbasislist(i, 1:2)
                            write(*, *) zcoefs(i, 1:2)
                          enddo




                    !--------------------------------------------------------------------------
                    else !bit not set, copy all previous, or seperately
                        write(*, *) 'bit not set, need to copy previous'
                        k = iend + 1
                        do i = istart, iend
                            zbasislist(k, 1:2) = zbasislist(i, 1:2)
                            zcoefs(k, 1:2) = zcoefs(i, 1:2)

                            tpsign = sign_ordets(singledets(j, tpspin), zbasislist(i, tpspin))
                            zbasislist(i, 1) = ior(singledets(j, 1), zbasislist(i, 1))
                            zbasislist(i, 2) = ior(singledets(j, 2) , zbasislist(i, 2))
                            call product_cmplx(singlecoefs(j, 1:2), zcoefs(i, 1:2), &
                                &zcoefs(i, 1:2))
                            write(*, *) zcoefs(i, 1:2)
                            zcoefs(i, 1) = tpsign * zcoefs(i, 1)
                            zcoefs(i, 2) = tpsign * zcoefs(i, 2)
                            write(*, *) zcoefs(i, 1:2)

                            tpsign = sign_ordets(singledets(j, 2 + tpspin), zbasislist(k, tpspin))
                            zbasislist(k, 1) = ior(singledets(j, 3), zbasislist(k, 1))
                            zbasislist(k, 2) = ior(singledets(j, 4), zbasislist(k, 2))
                            call product_cmplx(singlecoefs(j, 3:4), zcoefs(k, 1:2), &
                                &zcoefs(k, 1:2))
                            write(*, *) zcoefs(k, 1:2)
                            zcoefs(k, 1) = tpsign * zcoefs(k, 1)
                            zcoefs(k, 2) = tpsign * zcoefs(k, 2)
                            write(*, *) zcoefs(k, 1:2)
                            write(*, *) 'Processed dets number', i, k
                            write(*, '(2B16)') zbasislist(i, 1:2)
                            write(*, *) zcoefs(i, 1:2)
                            write(*, '(2B16)') zbasislist(k, 1:2)
                            write(*, *) zcoefs(k, 1:2)

                            k = k + 1
                        enddo
                        iend = k - 1
                        nzbasis = iend
                        write(*, *) ' '
                          write(*, *) 'Processed dets number', istart, iend
                          do i = istart, iend
                            write(*, '(2B16)') zbasislist(i, 1:2)
                            write(*, *) zcoefs(i, 1:2)
                          enddo




                     endif
                     !endif for bit set/notset in cases m!=0
                     !--------------------------------------------------------------------------
                 endif!endif for m ==0/ m!=0
              !===================================================================================

           enddo
           !end do j 

           write(*, *) 'nzbasis =', nzbasis, 'istart=', istart, 'iend=', iend
           do i = istart, iend
               zcoefs(i, 1) = zcoefs(i, 1) * realcoef
               zcoefs(i, 2) = zcoefs(i, 2) * realcoef
           enddo
           write(*, *) '========== Appended the following dets and coefs into the zlist ============'
           do i = istart, iend
               write(*, '(2B16)') zbasislist(i, 1:2)
               write(*, *) zcoefs(i, 1:2)
           enddo
           deallocate(singlecoefs)
           deallocate(singledets)

       end subroutine Y2Z_singledet




      subroutine Y2Z_det_incomplete(det, coef, zbasislist, zcoefs, nzbasis)
          integer(i16b), intent(in) :: det(:)
          real(rk), intent(in):: coef
          integer(i16b) ::  tpup, tpdn, tpplus, tpminus
          integer(i16b), allocatable :: singledets(:, :), zbasislist(:, :)
          real(rk), allocatable :: singlecoefs(:, :), zcoefs(:, :)
          integer :: pos, pos2, i, j,msign, k, nbasis0, npairnot0, nzbasis, nznot0
          !number of dets with m == 0, or m!= 0
          !zbasislist and zcoefs have to be allocated already
         
          integer :: neact !TODO for testing convenience, delete later
          neact = 2
    
          write(*, *) 'Entering Y2Z_det, input det pair is '
          write(*, '(2B16)') det(1:2)
          tpup = det(1)
          tpdn = det(2)
          allocate(singledets(3*neact, 2))
          allocate(singlecoefs(3*neact, 2))
          i = 0
          j = neact
          k = 2*neact

          write(*, *) 'Generating single electron det list'
          write(*, *) 'neact = ', neact
          do while(tpup.ne.0)
            pos = trailz(tpup)
            call sign_m(pos, pos2, tpplus, tpminus, msign)
            
            if(msign.eq.0) then
                k = k + 1
                singledets(k, 1) = tpplus !spin up
                singledets(k, 2) = 0_i16b
                singlecoefs(k, 1:2) = (/sqrt(2.), 0./)
                write(*, *)  'm = 0   This index ( 2*neact + ) = ', k
                write(*, '(2B16)') singledets(k, 1:2)
                write(*, *) singlecoefs(k, 1:2)
                tpup = ibclr(tpup, pos)


            elseif(BTEST(tpup, pos2)) then !msign not 0 and -m is set
                k = k + 1
                singledets(k, 1) = ior(tpplus, tpminus) 
                singledets(k, 2) = 0_i16b
                singlecoefs(k, 1:2) = (/0., -2./) ! -2i
                write(*, *) 'm, -m together This index ( 2*neact + ) = ', k
                write(*, '(2B16)') singledets(k, 1:2)
                write(*, *) singlecoefs(k, 1:2)
                tpup = ibclr(tpup, pos)
                tpup = ibclr(tpup, pos2)
            else

                i = i + 1
                j = j + 1
                singledets(i, 1) = tpplus
                singledets(i, 2) = 0_i16b
                singledets(j, 1) = tpminus
                singledets(j, 2) = 0_i16b
                singlecoefs(i, 1:2) =(/1.,0./)
                singlecoefs(j, 1:2) = (/0., msign*1. /)
                write(*, *) 'm != 0', '|m| index = ', i, '-|m| index = ', j
                write(*, '(4B16)') singledets(i, 1:2), singledets(j, 1:2)
                write(*, *) singlecoefs(i, 1:2), singlecoefs(j, 1:2)
                tpup = ibclr(tpup, pos)

            endif

          enddo
           
          do while(tpdn.ne.0)
             pos = trailz(tpdn)
             call sign_m(pos, pos2, tpplus, tpminus, msign)
 
             if(msign.eq.0) then
                 k = k + 1
                 singledets(k, 2) = tpplus !spin dn
                 singledets(k, 1) = 0_i16b
                 singlecoefs(k, 1:2) = (/sqrt(2.), 0./)
                 write(*, *)  'm = 0   This index ( 2*neact + ) = ', k
                 write(*, '(2B16)') singledets(k, 1:2) 
                 write(*, *) singlecoefs(k, 1:2)
                 tpdn = ibclr(tpdn, pos)
               

            elseif(BTEST(tpdn, pos2)) then !msign not 0 and -m is set
                k = k + 1
                singledets(k, 2) = ior(tpplus, tpminus) 
                singledets(k, 1) = 0_i16b
                singlecoefs(k, 1:2) = (/0., -2./) ! -2i
                write(*, *) 'm, -m together This index ( 2*neact + ) = ', k
                write(*, '(2B16)') singledets(k, 1:2)
                write(*, *) singlecoefs(k, 1:2)
                tpdn = ibclr(tpdn, pos)
                tpdn = ibclr(tpdn, pos2)
 

             else
 
                 i = i + 1
                 j = j + 1
                 singledets(i, 2) = tpplus
                 singledets(i, 1) = 0_i16b
                 singledets(j, 2) = tpminus
                 singledets(j, 1) = 0_i16b
                 singlecoefs(i, 1:2) =(/1.,0./)
                 singlecoefs(j, 1:2) = (/0., msign*1. /)
 
                 write(*, *) 'm != 0', '|m| index = ', i, '-|m| index = ', j
                 write(*, '(4B16)') singledets(i, 1:2), singledets(j, 1:2)
                 write(*, *) singlecoefs(i, 1:2), singlecoefs(j, 1:2)
                 tpdn = ibclr(tpdn, pos)
     
             endif
          enddo
          nbasis0 = k - 2*neact !number of dets that are real
          npairnot0 = i ! number of |m|, -|m| pairs
          nznot0 = 2**npairnot0
          !nadd = max(1, nznot0)   ! number of zdet pairs to generate

          allocate(zbasislist(nznot0, 2))
          allocate(zcoefs(nznot0, 2))
          !Now fill zbasislist and zcoefs
          if(npairnot0.gt.0) then
              do i = 1, nznot0/2
                  zbasislist(i, 1:2) = singledets(1, 1:2)
                  zcoefs(i, 1:2) = singlecoefs(1, 1:2)
                  j = i + nznot0/2
                  zbasislist(j, 1:2) = singledets(neact + 1, 1:2) 
                  zcoefs(j, 1:2) = singlecoefs(neact + 1, 1:2)
              enddo

              do k = 2, npairnot0
                  i = 0
                  do while(i.lt.(npairnot0-1))
                      !fill in det plus
                      do j = 1, 2**(npairnot0 - k)
                            i = i + 1
                            call product_cmplx(zcoefs(i, 1:2), singlecoefs(k, 1:2), zcoefs(i, 1:2))
                            !if the det to be combined is smaller, there will be a minus sign 
                            if(zbasislist(i, 1).gt.singledets(k, 1)) then
                                zcoefs(i, 1) = - zcoefs(i, 1)
                                zcoefs(i, 2) = - zcoefs(i, 2)
                            endif
                            if(zbasislist(i, 2).gt.singledets(k, 2)) then
                                zcoefs(i, 1) = -zcoefs(i, 1) 
                                zcoefs(i, 2) = -zcoefs(i, 1) 
                            endif
                            zbasislist(i, 1) = ior(zbasislist(i, 1), singledets(k, 1)) 
                            zbasislist(i, 2) = ior(zbasislist(i, 2), singledets(k, 2))
 
                      end do
                      !fill in det minus
                      do j = 1, 2**(npairnot0 - k)
                            i = i + 1
                            call product_cmplx(zcoefs(i, 1:2), singlecoefs(k+neact, 1:2), &
                                & zcoefs(i, 1:2))
                            if(zbasislist(i, 1).gt.singledets(k + neact, 1)) then
                                zcoefs(i, 1) = - zcoefs(i, 1)
                                zcoefs(i, 2) = - zcoefs(i, 2)
                            endif
                            if(zbasislist(i, 2).gt.singledets(k, 2)) then
                                zcoefs(i, 1) = -zcoefs(i, 1) 
                                zcoefs(i, 2) = -zcoefs(i, 1) 
                            endif
 
                            zbasislist(i, 1) = ior(zbasislist(i, 1), singledets(k+neact, 1))
                            zbasislist(i, 2) = ior(zbasislist(i, 2), singledets(k+neact, 2))
 
                       enddo
                   enddo!end while
                   write(*, *) 'For k =', k, 'last i =', i
              enddo

              if(nbasis0.ne.0) then
                  do i = 1, nznot0
                      do j = 1, nbasis0
                          zbasislist(i, 1) = ior(zbasislist(i, 1), singledets(j + neact*2, 1))
                          zbasislist(i, 2) = ior(zbasislist(i, 2), singledets(j + neact*2, 2))
                          call product_cmplx(zcoefs(i, 1:2), singlecoefs(j+neact*2, 1:2), &
                                & zcoefs(i, 1:2))
 
                      enddo
                  enddo
              endif

              ! end of if npairnot0 > 0
          else if(nbasis0.ne.0) then
              ! npairnot0 = 0 all m = 0, or m, -m together 
              ! put all original ydet into zbasis
              zbasislist(1, 1:2) = det(1:2)
              zcoefs(1, 1:2) = (/1., 0./)
              do i = 1, nbasis0
                  call product_cmplx(zcoefs(1, 1:2), singlecoefs(i + neact*2, 1:2), zcoefs(1, 1:2))
              enddo
          endif
    


          write(*, *) 'Final Z zbasislist, zcoefs for this csf'
          do i = 1, nznot0
              write(*, '(2B16)') zbasislist(i, 1:2)
              write(*, *) zcoefs(i, 1), ' + i', zcoefs(i, 2)
          enddo
          deallocate(singlecoefs)
          deallocate(singledets)
      end subroutine Y2Z_det_incomplete


      subroutine product_cmplx(c1, c2, nout)
          real(rk) :: c1(:), c2(:)
          real(rk) :: n1(2), n2(2), nout(:)
          n1(1:2) = c1(1:2)
          n2(1:2) = c2(1:2)
          nout(1) = n1(1)*n2(1) - n1(2)*n2(2)
          nout(2) = n1(1)*n2(2)+ n1(2)*n2(1)
      end subroutine product_cmplx
          
      subroutine sign_m( pos, pos2, detplus, detminus, msign)
          ! for a single det not a pair
          integer, intent(in) :: pos
          integer(i16b) :: detplus, detminus
          integer :: pos2, n, l, m, msign
          call delocate(pos, n, l, m)
          if (m.eq.0) then
              msign = 0
              detplus = ibset(0_i16b, pos)
              detminus = 0_i16b
              pos2 = pos
          elseif(m.gt.0) then
              msign = 1
              detplus = ibset(0_i16b, pos)
              pos2 =  locate_det(n, l, -m)
              detminus = ibset(0_i16b, pos2)
          else
              msign = -1
              detminus = ibset(0_i16b, pos)
              pos2 = locate_det(n, l, -m)
              detplus = ibset(0_i16b, pos2)
          endif
          write(*, *) 'subroutine sign_m : pos, pos2, detplus, detminus'
          write(*, '(2I5, 2b16)') pos, pos2, detplus, detminus, msign
      end subroutine sign_m



end module convertharmonics
