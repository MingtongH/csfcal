module prep
  
  use, intrinsic :: iso_fortran_env, only: rk => real64

  implicit none
  
  ! Constants
  integer, parameter :: STRING_MAX_LENGTH = 16
  integer, parameter :: i16b = SELECTED_INT_KIND(38)

  ! Input variables
  integer :: nconfig, Ldes, Sdes_t2, Ides
  !Lzdes, Szdes_t2 not defined here
  !Ides is read in as 1/-1, but then converted to 0/1 (even/odd)
  integer :: neinact ! not used
  integer :: neact
  integer, allocatable :: econfigs (:, :) !n, l, ne
  integer, allocatable :: nshell(:) ! number of shells occupied in each config
  
  !call read_input()
  !call checkSmin() !tested
  !call checkSmax() !tested
  !write(*, '(">>>testing ", 1I5)') Lmin(4, 6)
  !call checkILmax()
    
!  testarray(1:5) = 1
!  testarray(1) = 8
!  testarray(2) = 11
!  testarray(3) = 6
!  testarray(4) = 5
!  write(*, *) minsum(testarray)
  !call checkLmin()
  contains

  subroutine read_input()
      
    implicit none 
    !integer, allocatable :: ne(:)
    integer :: i, j, k, tmpne, tmpn, tmpl
    character(1) :: tmporb
    character(STRING_MAX_LENGTH) :: line

    write(*, '("==== Start Reading Input ====")')
    read(*, *) line
    write(*, '("Atom: ", A)') line
    read(*, *) nconfig
    write(*, '("Number of configs: ", 1I5)') nconfig
    read(*, *) neact
    write(*, '("Number of active electrons: ", 1I5)') neact
    read(*, *) neinact
    write(*, '("Number of inactive electrons: ", 1I5)') neinact
    allocate(econfigs(nconfig * neact, 3 )) !not full
    allocate(nshell(nconfig))
    k = 1
    do i = 1, nconfig
         read(*, *)
         write(*, '(">>> Reading Config ", 1I5)') i
         read(*, *) nshell(i)
         write(*, '("Number of shells occupied: ", 1I5)') nshell(i)
    !     allocate(ne(nshell(i))
         do j = 1, nshell(i)
             read(*, *) tmpn, tmporb, tmpne !ne(j)
             tmpl = orbTol(tmporb)
             econfigs(k, 1) = tmpn
             econfigs(k, 2) = tmpl
             econfigs(k, 3) = tmpne !ne(j)Smin*2 = ", 1I5)') Smin_t2(neact://www.phy.ornl.gov/csep/pl/node17.html

             write(*, '("shell # : n, l, ne ", 4I5)') j, tmpn, tmpl, tmpne

             k = k + 1
             !do k = 1, ne(j)
             !    tmpne = (i-1)*neact + sum(ne(1:j-1)) + k
             !    write(*, '("Electron # : n, l ", 3I5)') tmpne, tmpn, tmpl
             !    econfigs(tmpne, 1) = tmpn
            !    econfigs(tmpne, 2) = tmpl
            ! enddo !this shell
         enddo!all shells
         !deallocate(ne)
     enddo

     read(*, *) 
     write(*, '(">>> Reading other parameters")')
     read(*, *) Sdes_t2
     read(*, *) Ldes
     !read(*, *) Szdes_t2
     !read(*, *) Lzdes
     read(*, *) Ides
     !write(*, '("desired S*2 L Sz*2 Lz I :", 5I5)') Sdes_t2, Ldes, Szdes_t2,  Lzdes, Ides
     write(*, '("desired S*2 L Sz*2 Lz I :", 5I5)') Sdes_t2, Ldes, Ides
     if(Ides.eq.1) then 
         Ides = 0
     elseif(Ides.eq.-1) then
         Ides = 1
     else
         write(*, '("Error: Input for desired I not correct. 1 for even and -1 for odd ")') 
         stop 1
     endif
     write(*, '("Ides converted")')
     write(*, '("==== End Reading Input ====")')

  end subroutine read_input
  
  
  function orbTol(orb)
      character(1), intent(in) :: orb
      integer :: l, orbTol
      select case (orb)
        case ('s')
            l = 0
        case ('p')
            l = 1 
        case ('d')
            l = 2
        case ('f')
            l = 3
       end select
       orbTol = l
  end function orbTol
  

  function Smin_t2(ne)
      integer, intent(in) :: ne
      integer :: Smin_t2
      Smin_t2 = mod(ne, 2)
  end function Smin_t2
  
  
  subroutine checkSmin()
      if(Sdes_t2.lt.Smin_t2(neact)) then
          write(*, '("Error: S_des < S_min, current Smin*2 = ", 1I5)') Smin_t2(neact)
          stop 1
      endif
      write(*, '("Smin checked. ")')
  end subroutine checkSmin
  !function Smax_t2()
  
  function Smax_t2_oneshell(l, N)
      integer, intent(in) :: l, N
      integer :: num_orbs, Smax_t2_oneshell
      num_orbs = 2*l + 1
      if(N.le.num_orbs) then
          Smax_t2_oneshell = N
      else if (N.gt.(2*num_orbs)) then
          write(*, '("Error: number of electrons more than full on orbital of l = : ", 1I5)') l
          stop 1
      else
          Smax_t2_oneshell = num_orbs - (N - num_orbs)
      endif
  end function Smax_t2_oneshell



  function Smax_t2(istart, iend)
      integer, intent(in) :: istart, iend
      integer :: i, Smax_t2
      Smax_t2 = 0
      do i = istart, iend
        Smax_t2 = Smax_t2 + Smax_t2_oneshell(econfigs(i, 2), econfigs(i, 3))
      end do 
  end function Smax_t2

  subroutine checkSmax()
      integer :: istart, iend, tpSmax_t2, i
      iend = 0
      do i = 1, nconfig
        istart = iend + 1
        iend = istart + nshell(i) - 1
        tpSmax_t2 = Smax_t2(istart, iend)
        if(Sdes_t2.gt.tpSmax_t2) then
            write(*, '("Error: S_des > S_max, in config #", 1I5)') i
            write(*, '("Current Smax * 2 = ", 1I5)') tpSmax_t2
            stop 1
        endif
      enddo
      write(*, '("Smax checked")')
  end subroutine checkSmax

  function Lmax(istart, iend)
      integer, intent(in) :: istart, iend
      integer :: Lmax, i
      Lmax = 0
      do i = istart, iend
        Lmax = Lmax + econfigs(i, 2)*econfigs(i, 3)
      end do
  end function Lmax
  
  subroutine checkILmax()
      integer :: istart, iend, tpLmax, i, tpI
      iend = 0
      do i = 1, nconfig
        istart = iend + 1
        iend = istart + nshell(i) - 1
        tpLmax = Lmax(istart, iend)
        tpI = mod(tpLmax, 2) !0 for even, 1 for odd

        if(Ldes.gt.tpLmax) then 
            write(*, '("Error: L_des > L_max, in config #", 1I5)') i
            write(*, '("Current Lmax  = ", 1I5)') tpLmax
            stop 1
        endif
        if(Ides.ne.tpI) then
            write(*, '("Error: Wrong parity in config #", 1I5)') i
            stop 3
        endif
      enddo
      write(*, '("I checked. ")')
      write(*, '("Lmax checked. ")')
  end subroutine checkILmax

   function minsum(array)
       implicit none
       integer :: n, array(:), i, j, minsum, tmpsum
       integer(i16b) :: b
       b = 0
       n = size(array)
       minsum = 100*n
       do i = 1, 2**n
           tmpsum = 0
           do j = 1, n
               if (btest(b,j)) then
                   tmpsum = tmpsum - array(j)
               else
                   tmpsum = tmpsum + array(j)
               endif
           enddo
           tmpsum = abs(tmpsum)
           if (tmpsum.lt.minsum) then
               minsum =tmpsum
           end if
           b = b + 1
       enddo
 
   end function minsum


  
  function Lmin(istart, iend)
      integer, intent(in) :: istart, iend
      integer :: i, j, k, Lmin
      integer, dimension(neact) :: larray
      k = 1
      do i = istart, iend
        do j = 1, econfigs(i, 3) 
          larray(k) = econfigs(i, 2)
         ! write(*, '("electron #: , l=", 2I5)') k, econfigs(i, 2)
          k = k + 1
        enddo
      enddo
      Lmin = minsum(larray)
  end function Lmin




  subroutine checkLmin()
      integer :: istart, iend, tpLmin, i
      iend = 0
      do i = 1, nconfig
        istart = iend + 1
        iend = istart + nshell(i) - 1
        !write(*, *) istart, iend, Lmin(istart, iend)
        tpLmin = Lmin(istart, iend)
        if(Ldes.lt.tpLmin) then 
            write(*, '("Error: L_des < L_min, in config #", 1I5)') i
            write(*, '("Current Lmin = ", 1I5)') tpLmin
            stop 1
        endif
      enddo
      write(*, '("Lmin checked. ")')
  end subroutine checkLmin


end module prep 

    
 

