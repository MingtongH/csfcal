program csfcal

  use, intrinsic :: iso_fortran_env, only: rk => real64

  implicit none
  
  ! Configuration
  integer, parameter :: TITLE_MAX_LENGTH = 16
  integer, parameter :: BUF_LENGTH = 64
  ! Input variables
  character(TITLE_MAX_LENGTH) :: title
  integer :: nconfig
  integer :: neinact ! not used
  integer :: neact
  integer, allocatable :: econfigs (:, :) !
  call read_input()
  contains

  subroutine read_input()
    
    implicit none 
    integer :: nshell
    integer, allocatable :: ne(:)
    integer :: i, j, k, tmpne, tmpn, tmpl
    write(*, '("==== Start Reading Input ====")')
    read(*, *) title
    write(*, '("Title: ", A)') title
    read(*, *) nconfig
    write(*, '("Number of configurations: ", A)') trim(str(nconfig))
    read(*, *) neinact
    write(*, '("Number of inactive electrons: ", A)') trim(str(neinact))
    read(*, *) neact
    write(*, '("Number of active electrons: ", A)') trim(str(neact))
    allocate(econfigs(nconfig * neact, 2 ))
   
    do i = 1, nconfig 
        read(*, *)
        write(*, '(">>> Reading Config ", A)') trim(str(i))
        read(*, *) nshell
        allocate(ne(nshell))
        write(*, '("Number of shells occupied: ", A)') trim(str(nshell))
        do j = 1, nshell
            read(*, *) ne(j), tmpn, tmpl
            do k = 1, ne(j)
                !tmpne = (i-1)*neact + sum(ne(1:j)) - ne(j) + k
                tmpne = (i-1)*neact + sum(ne(1:j-1)) + k
                write(*, '("Electron : n, l ", 3I5)') tmpne, tmpn, tmpl
                econfigs(tmpne, 1) = tmpn
                econfigs(tmpne, 2) = tmpl
            enddo !this shell
        enddo!all shells
        deallocate(ne)
    enddo
    write(*, '("==== End Reading Input ====")')
  end subroutine read_input

  function str(num_in)
     integer, intent(in) :: num_in
     character(BUF_LENGTH) :: str
     write (str, *) num_in
     str = adjustl(str)
 end function str

end program csfcal 

    
 


