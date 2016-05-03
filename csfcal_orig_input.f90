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
    write(*, '("Title: ", A)') line
    read(*, *) line
    line = trim(line)
    write(*, '("Line 2: ", A)') line
    






  end subroutine read_input


end program csfcal 

    
 


