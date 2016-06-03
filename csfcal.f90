program csfcal
       use prep
       implicit none
       call main()
 
 
       contains
 
       subroutine main
           call read_input()
           call checkSmin()
           call checkSmax()
           call checkILmax()
           call checkLmin()
       end subroutine main
 end program csfcal

