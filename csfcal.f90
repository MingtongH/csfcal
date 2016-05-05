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
  call read_input()
  call checkSmin()
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
          write(*, '("Error: S_des < S_min, Smin*2 = ", 1I5)') Smin_t2(neact)
          stop 1
      endif
  end subroutine checkSmin
  !function Smax_t2()
  function Lmax(iconfig)
      integer, intent(in) :: iconfig
      integer Lmax, i,  istart, iend
      Lmax = 0
      if(iconfig.eq.1) then
          istart = 1
      else
          istart = sum(nshell(1:iconfig-1))+1
      end if
      iend = istart + nshell(iconfig)-1

      do i = istart, iend
        Lmax = Lmax + econfigs(i, 2)*econfigs(i, 3)
      end do
  end function Lmax
  end program csfcal 

    
 


