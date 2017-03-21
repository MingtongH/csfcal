program matMulProduct
   use orthogonalization
   real(rk), allocatable :: A(:, :), Q(:, :), R(:, :), E(:, :)
   integer :: i, j, n, p, ncsf
   n = 2
   p = 3
   allocate(A(n, p))
   A(:, 1) = (/2,-1/)
   A(:, 2) = (/-2, 1/)
   A(:, 3) = (/3, 0/)
  ! call gs(A, Q, R, E, n, p, ncsf)
   call gs_modify(A, n, p)



   !write(*, *) 'Q = '
   !do i = 1, n
   !    write(*, *) Q(i, :)
   !enddo

  ! write(*, *) 'R = '
  ! do i = 1, p
  !     write(*, *) R(i, :)
  ! enddo

!   write(*, *) 'E = '
!   do i = 1, n
!       write(*, *) E(i, 1:ncsf)
!   enddo

       
end program matMulProduct
