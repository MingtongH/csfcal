program csfcal
       use prep, only : i16b, Lzmax, Lzdes, Sdes_t2, econfigs, nshell, &
           & read_input, checkSmin, checkSmax, checkILmax
       use detgen, only : assign_shell

       implicit none
       call main()
 
 
       contains
 
       subroutine main
           integer :: totdets 
           integer(i16b), allocatable :: detlist(:, :) 

           call read_input()
           call checkSmin()
           call checkSmax()
           call checkILmax()
          !call checkLmin() !Lmin not that simple
           totdets = 0
           write(*, *) '===============Assign_shell for Config 1===================='
           call assign_shell(0_i16b, 0_i16b, 0, 0, 0, 0,&
           & Lzdes, Sdes_t2, econfigs(1:nshell(1), 1:3), 0, totdets, detlist)
       end subroutine main
 end program csfcal

