program csfcal
  
  use, intrinsic :: iso_fortran_env, only: rk => real64

  implicit none
  
  ! Constants
  integer, parameter :: STRING_MAX_LENGTH = 16
  ! Input variables
  integer :: nconfig, Ldes, Sdes_t2, Lzdes, Szdes_t2
  integer :: neinact ! not used
  integer :: neact
  integer, allocatable :: econfigs (:, :) !
  integer, allocatable :: nshell(:)
  
  integer, dimension(5) :: testarray
  call read_input()
  !call checkSmin() !tested
  !call checkSmax() !tested
  !write(*, '(">>>testing ", 1I5)') Lmin(4, 6)
  !call checkLmax()
  testarray(1:5) = 1
  write(*, *) minsum(testarray, 5)
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
             econfigs(k, 3) = tmpne !ne(j)

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
     read(*, *) Szdes_t2
     read(*, *) Lzdes
     write(*, '("desired S*2 L Sz*2 Lz:", 4I5)') Sdes_t2, Ldes, Szdes_t2,  Lzdes
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
          write(*, '("Error: S_des < S_min, change neact, Smin*2 = ", 1I5)') Smin_t2(neact)
          stop 1
      endif
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
            write(*, '("S_des * 2 <= ", 1I5)') tpSmax_t2
            stop 1
        endif
      enddo
  end subroutine checkSmax

  function Lmax(istart, iend)
      integer, intent(in) :: istart, iend
      integer :: Lmax, i
      Lmax = 0
      do i = istart, iend
        Lmax = Lmax + econfigs(i, 2)*econfigs(i, 3)
      end do
  end function Lmax
  
  subroutine checkLmax()
      integer :: istart, iend, tpLmax, i
      iend = 0
      do i = 1, nconfig
        istart = iend + 1
        iend = istart + nshell(i) - 1
        tpLmax = Lmax(istart, iend)
        if(Ldes.gt.tpLmax) then 
            write(*, '("Error: L_des > L_max, in config #", 1I5)') i
            write(*, '("L_des  <= ", 1I5)') tpLmax
            stop 1
        endif
      enddo
  end subroutine checkLmax

  function minsum(array, n)
      integer, intent(in) :: n, array(n)
      integer :: i, minsum
      minsum = sum(array(:))
      !! TO DO
  end function minsum


  
  function Lmin(istart, iend)
      integer, intent(in) :: istart, iend
      integer :: i, j, k, Lmin
      integer, dimension(neact) :: larray
      k = 1
      do i = istart, iend
        do j = 1, econfigs(i, 3) 
          larray(k) = econfigs(i, 2)
          write(*, '("electron #: , l=", 2I5)') k, econfigs(i, 2)
          k = k + 1
        enddo
      enddo
      Lmin = minsum(larray, neact)
  end function Lmin




  subroutine checkLmin()
      !!TO DO
  end subroutine checkLmin


end program csfcal 

    
 


