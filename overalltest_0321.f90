program overalltest
       use prep, only : i16b, rk, Ldes, Lzdes, Sdes_t2, econfigs, nconfig, nshell, neact, &
            !parameters and variables 
           & read_input, checkSmin, checkSmax, checkILmax,checkLmin, Lmax, Lmin, &
               & Smin_t2, Smax_t2, ARRAY_SHORT_LENGTH
            !subroutines and functions
       use detgen, only : assign_shell
       use projection, only : initlists, Proj_L, Proj_S, getallsigns, collect_csf, normalizetable, initlists_fromlists
       use checkcsfs, only: sortBasisCoefTable, Lsq_multiple, Ssq_multiple
       !use gramschmidt, only : orth
       use orthogonalization, only: gs_modify

       implicit none
       call main()
 
 
       contains
 
       subroutine main
           integer :: totdets, minLp, maxLp, innum, num1, num, i, nbasis, ncsf, iconf,&
               & iend, istart, min2Sp, max2Sp
           integer(i16b), allocatable :: detlist(:, :), inbasis(:, :), &
               &basislist1(:, :),basislist(:, :), allbasis(:, :), allbasis1(:, :)
           real(rk), allocatable :: incoefs(:), coeflist1(:), coeflist(:), coeftable(:, :), Q(:, :), R(:, :)
           integer, allocatable :: ineposes(:, :, :), eposlist1(:, :, :), eposlist(:, :, :)
           logical, allocatable :: iniszeros(:), iszerolist1(:), iszerolist(:)
           character(10) :: line 

           call read_input()
           call checkSmin()
           call checkSmax()
           call checkILmax()
           call checkLmin()
           totdets = 0
           nbasis = 0
           ncsf = 0
           iend = 6
           do iconf = 4, 4
               istart = iend + 1
               iend = istart + nshell(iconf) - 1
               write(*, *) '%%%%%%%%%%%%%%%%%%%%%%  Assign_shell for Config', iconf, '  %%%%%%%%%%%%%%%%%%%'


               totdets = 0
               call assign_shell(0_i16b, 0_i16b, 0, 0, 0, 0,&
                  & Lzdes, Sdes_t2, econfigs(istart:iend, 1:3), 0, totdets, detlist)
               write(*, *) totdets
               do i = 1, totdets
                   write(*, '(2B16)') detlist(i, 1:2)
               enddo
               minLp = Lmin(econfigs(istart:iend, 1:3))
               maxLp = Lmax(econfigs(istart:iend, 1:3))
               min2Sp = Smin_t2(neact)
               max2Sp = Smax_t2(econfigs(istart:iend, 1:3))
               ncsf = 0 !Each configs has a csftable,n 
                        !so set size to 0 before filling for current config
               do i = 1, totdets

                   write(*, *) 'xxxxxxxxxxxxxxxxxxxxx  Projection of det #', i,'  xxxxxxxxxxxxxxxxxxxx' 
                   call initlists(detlist(i, 1:2), inbasis, incoefs, ineposes, iniszeros, innum)
                   !inbasis... always initialized to arrays with only one det
                   call Proj_L(minLp, maxLp, Ldes, &
                   & inbasis, incoefs, ineposes, iniszeros, innum, &
                   & basislist1, coeflist1, eposlist1, iszerolist1, num1)

                   call Proj_S(min2Sp, max2Sp, Sdes_t2, &
                       &basislist1, coeflist1, eposlist1, iszerolist1, num1, &
                       &basislist, coeflist, eposlist, iszerolist, num)


                   !PAUSE 4
                   ! write(*,'(''paused, type [enter] to continue'')')
                   !       read (*,*) line


                   call getallsigns(basislist, coeflist, eposlist, iszerolist, num)
                   call collect_csf(basislist, coeflist, iszerolist, num, &
                     & allbasis, coeftable, nbasis, ncsf)
               enddo !for each det

               call sortBasisCoefTable(allbasis, coeftable, nbasis, ncsf)
 

               call normalizetable(coeftable, nbasis, ncsf)
               
               write(*, *) nbasis, ncsf
               do i = 1, nbasis
                 write(*, *) coeftable(i, 1:ncsf)
               enddo
               if(nbasis.ne.ncsf) then
                   write(*, *) 'coeftable not square'
                   exit
               endif

               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               !All correct up to this point
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               

               call gs_modify(coeftable, nbasis, ncsf)

           enddo ! for each config




       end subroutine main
 end program overalltest
