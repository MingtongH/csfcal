module detgen
  use prep, only: i16b, Lmax, Smax_t2
  implicit none
  contains

      function isHalfFull(shellconfig)
          implicit none
          integer, dimension(1,3), intent(in) :: shellconfig
          integer :: l, ne
          logical :: isHalfFull
          l = shellconfig(1, 2)
          ne = shellconfig(1, 3)
          if (ne.ge.(2*l+1)) then 
              isHalfFull = .true.
          else
              isHalfFull = .false.
          endif
      end function isHalfFull

      recursive subroutine assign_shell(prdet, prLz, pr2Sz, maxremLz, maxrem2Sz, desLz, des2Sz, configs, ncurr, nshell_tot)
          implicit none
          integer, intent(in) :: prLz, pr2Sz, maxremLz, maxrem2Sz, ncurr
          integer :: nshell_tot, desLz, des2Sz
          integer, intent(in) :: configs(:, :)
          integer(i16b), intent(in) :: prdet
          integer(i16b) :: curdet
          integer :: curlz, cur2sz, curlmax,cursmax_t2, tpLz, tp2Sz, tpremLz, tpremSz
        
          write(*, '(">>>>>>>>>>>>>>>>Now at shell :",1I5)') ncurr
          write(*, '("prdet, prLz, pr2Sz, maxremLz, maxrem2Sz", 5I5)') prdet, prLz, pr2Sz, maxremLz, maxrem2Sz

          !last shell
          if(ncurr.eq.nshell_tot) then
            curlz = desLz - prLz
            cur2sz = des2Sz - pr2Sz
            if(curlz.gt.maxremLz.OR.curlz.lt.(-maxremLz)) then
                return
            else if (cur2sz.ne.maxrem2Sz.OR.cur2sz.ne.(-maxrem2Sz)) then
                return
            else 
                write(*, '("isHalfFull = ")') isHalfFull(configs(ncurr:ncurr, 1:3))
                return
            endif
          !first shell
          else if(ncurr.eq.0) then
              curdet = 0!
              call assign_shell(curdet, 0, 0, Lmax(configs), Smax_t2(configs), desLz, des2Sz, configs, ncurr+1, nshell_tot)
          !recur
          else 
              curlmax = Lmax(configs(ncurr:ncurr, 1:3))
              cursmax_t2 = Smax_t2(configs(ncurr:ncurr, 1:3))
              do curlz = -curlmax, curlmax!
                curdet = prdet + 1 !
                cur2sz = Smax_t2(configs(ncurr:ncurr, 1:3))!
                tpLz = prLz + curLz
                tp2Sz = pr2Sz + cur2sz
                tpremLz = maxremLz - curlmax
                tpremSz = maxrem2Sz - cursmax_t2
                write(*, '("*******At shell, Entering branch curlz = ", 2I5)') ncurr, curlz
                call assign_shell(curdet, tpLz, tp2Sz, tpremLz, tpremSz, desLz, des2Sz, configs, ncurr+1, nshell_tot)
              enddo
          endif
      end subroutine assign_shell

                

            

          

      !function Lz(det) 
      !  implicit none
      !  integer(i16b), intent(in) :: det(1, 2)
      !  integer :: Lz
      !  Lz = 0
      !  write(*, *) det(1, 1)
      !  write(*, *) det(1, 2)
      !end function Lz


      !subroutine gen_alldets(config, ndets, alldets)
       ! implicit none
       ! integer, intent(in) :: config(:, 3)
       ! integer :: ndets
       ! integer, allocatable :: alldets(:, :)

    
end module detgen
