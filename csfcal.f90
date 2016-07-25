program csfcal
       use prep, only : i16b, rk, Ldes, Lzdes, Sdes_t2, econfigs, nshell, &
            !parameters and variables 
           & read_input, checkSmin, checkSmax, checkILmax, Lmax
            !subroutines and functions
       use detgen, only : assign_shell
       use projection, only : initlists, Proj_L, getallsigns, collect_csf, normalizetable

       implicit none
       call main()
 
 
       contains
 
       subroutine main
           integer :: totdets, minLp, maxLp, innum, num, i, nbasis, ncsf, nconf, iend, istart
           integer(i16b), allocatable :: detlist(:, :), inbasis(:, :), basislist(:, :), allbasis(:, :)
           real(rk), allocatable :: incoefs(:), coeflist(:), coeftable(:, :)
           integer, allocatable :: ineposes(:, :, :), eposlist(:, :, :)
           logical, allocatable :: iniszeros(:), iszerolist(:)
           

           call read_input()
           call checkSmin()
           call checkSmax()
           call checkILmax()
          !call checkLmin() !Lmin not that simple
           totdets = 0
           nbasis = 0
           ncsf = 0
           iend = 6
           nconf = 4!TODO loop
           istart = iend + 1
           iend = istart + nshell(nconf) - 1
           write(*, *) '%%%%%%%%%%%%%%%%%%%%%%  Assign_shell for Config', nconf, '  %%%%%%%%%%%%%%%%%%%'
           call assign_shell(0_i16b, 0_i16b, 0, 0, 0, 0,&
           & Lzdes, Sdes_t2, econfigs(istart:iend, 1:3), 0, totdets, detlist)
           write(*, *) totdets
           do i = 1, totdets
               write(*, '(2B16)') detlist(i, 1:2)
           enddo
           do i = 1, totdets
             write(*, *) 'xxxxxxxxxxxxxxxxxxxxx  Projection of det #', i,'  xxxxxxxxxxxxxxxxxxxx' 

             call initlists(detlist(i, 1:2), inbasis, incoefs, ineposes, iniszeros, innum)
             minLp = 0 
             maxLp = Lmax(econfigs(istart:iend, 1:3))
             call Proj_L(minLp, maxLp, Ldes, &
               & inbasis, incoefs, ineposes, iniszeros, innum, &
               & basislist, coeflist, eposlist, iszerolist, num)

             !TODO call Proj_L

             call getallsigns(basislist, coeflist, eposlist, iszerolist, num)
             call collect_csf(basislist, coeflist, iszerolist, num, &
               & allbasis, coeftable, nbasis, ncsf)
           enddo

           call normalizetable(coeftable, nbasis, ncsf)
           write(*, *) coeftable(1:nbasis, 1:ncsf)




       end subroutine main
 end program csfcal

