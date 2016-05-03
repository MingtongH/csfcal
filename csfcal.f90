program csfcal
  
  use, intrinsic :: iso_fortran_env, only: rk => real64

  implicit none
  
  ! Constants
  integer, parameter :: STRING_MAX_LENGTH = 16
  ! Input variables
  character(STRING_MAX_LENGTH) :: line
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
    read(*, *) line
    write(*, '("Atom: ", A)') line
    read(*, *) nconfig
    write(line, *) nconfig
    write(*, '("Number of configs: ", A)') adjustl(line)
    read(*, *) neact
    write(line, *) neact
    write(*, '("Number of active electrons: ", A)') adjustl(line)
    read(*, *) neinact
    write(line, *) neinact
    write(*, '("Number of inactive electrons: ", A)') adjustl(line)








  end subroutine read_input


end program csfcal 

    
 


