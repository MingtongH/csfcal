module detgen
  use prep, only: i16b, Lzmax, Smax_t2
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
          !input maxrem includes all shells not assigned
          !including the one to be assigned in this recursion
          !only calculated by Lzmax(), no consideration of desLz
          implicit none
          integer, intent(in) :: prLz, pr2Sz, maxremLz, maxrem2Sz, ncurr, nshell_tot
          ! ncurr = current shell number
          integer, intent(in) :: desLz, des2Sz
          integer, intent(in) :: configs(:, :)
          integer(i16b), intent(in) :: prdet
          integer(i16b) :: curdet, tpdet, tpdetlen, ecount
          integer :: curlz, cur2sz, i !iterator
          integer :: curLzmax, curLzmin, m_max, m_min, cur2Szmax, cur2Szmin
          integer :: tpLz, tp2Sz, tpmaxremLz, tpmaxrem2Sz
          integer :: necur ! number of electrons in current shell
          
          write(*, '("nshell_tot = ", 1I5)') nshell_tot
          write(*, '(">>>>>>>>>>> Now at shell :",1I5)') ncurr
          write(*, '("prdet, prLz, pr2Sz, maxremLz, maxrem2Sz", 5I5)') prdet, prLz, pr2Sz, maxremLz, maxrem2Sz

          !*************************last shell
          if(ncurr.eq.nshell_tot) then
            necur = configs(ncurr, 3) 
            curlz = desLz - prLz
            cur2sz = des2Sz - pr2Sz
            if(curlz.gt.maxremLz.OR.curlz.lt.(-maxremLz)) then
                write(*, '("Cannot assign for des with curlz = ", 1I5)') curlz
                return
            !else if (cur2sz.gt.maxrem2Sz.OR.cur2sz.lt.(-maxrem2Sz)) then
            !    return
            else 
                write(*, '("isHalfFull = ")') isHalfFull(configs(ncurr:ncurr, 1:3))
                !!TODO if there is valid assignment, update curdet
                return
            endif
          !**************************first shell(only need inputs as 0)
          else if(ncurr.eq.0) then
              curdet = 0!
              call assign_shell(curdet, 0, 0, Lzmax(configs), Smax_t2(configs), desLz, des2Sz, configs, ncurr+1, nshell_tot)
          !*****************************recur
          else
              necur = configs(ncurr, 3)
              tpLz = Lzmax(configs(ncurr:ncurr, 1:3))
              tpmaxremLz = maxremLz - tpLz ! tpmaxremLz no longer include current shell
              curLzmax = min(tpLz, desLz - prLz + tpmaxremLz)
              curLzmin = max(-tpLz, desLz - prLz - tpmaxremLz)
              m_max = configs(ncurr, 2) 
              m_min = - configs(ncurr, 2) 
              !!TODO May need to change m_max, m_min for pruning
              
              !!TODO check Sz beforehand to save time
              !cursmax_t2 = Smax_t2(configs(ncurr:ncurr, 1:3))
              tpmaxrem2Sz = 1
              
              !! TODO For each assignment, if within lz range and sz range, do the following:
              tpdetlen = m_max - m_min + 1
              tpdet = 0
              do i = 2*tpdetlen - necur + 1, 2*tpdetlen
                tpdet = ibset(tpdet, i)
              enddo
              ! use just one combined det dn_det-up_det
              ! tpdet is the last curdet with largest value

              do curdet = 1, tpdet
              !TODO if occ_count = necur
              !      if Lz within range
                 ! call next shell recursive function


                !call get_occ_orbs(ibits(curdet, tpdetlen,tpdetlen), ibits(curdet, 2*tpdetlen, tpdetlen), occ_up, occ_dn, ecount)
               ! if (ecount.eq.necur) then
                    !Calculate curLz, check range
               !     curLz = 
                    

              !do curlz = m_min, m_max!
              !  cur2sz = Smax_t2(configs(ncurr:ncurr, 1:3))!
              !  tpLz = prLz + curLz
              !  tp2Sz = pr2Sz + cur2sz
              !  curdet = prdet + 1 !
                !Update curdet
              !  write(*, '("*** Branch shell, curlz = ", 2I5)') ncurr, curlz
              !  call assign_shell(curdet, tpLz, tp2Sz, tpmaxremLz, tpmaxrem2Sz, desLz, des2Sz, configs, ncurr+1, nshell_tot)
              enddo
            !******************************end of recur
          endif
      end subroutine assign_shell

                

            

          
  integer function count_occ_orbs(det_up, det_dn)
      integer(i16b), intent(in) :: det_up, det_dn
      integer(i16b) :: tmp_det
      integer :: i, nocc
    tmp_det = det_up
    !i_elec = 0
    nocc = 0
    do while (tmp_det .ne. 0)
      i = trailz(tmp_det)+1
   !   i_elec = i_elec + 1
   !   occ_up(i_elec) = i
      if(i.gt.0) then 
          nocc = nocc + 1
      endif
      tmp_det = ibclr(tmp_det,i-1)
    enddo

    tmp_det = det_dn
    !i_elec = 0
    do while (tmp_det .ne. 0)
      i = trailz(tmp_det)+1
    !  i_elec = i_elec + 1
    !  occ_dn(i_elec) = i
      if(i.gt.0) then 
          nocc = nocc + 1
      endif   
      tmp_det = ibclr(tmp_det,i-1)
    enddo
    count_occ_orbs = nocc
 end function count_occ_orbs


  subroutine get_occ_orbs(det_up,det_dn,occ_up,occ_dn, nocc)
  ! Get lists of which orbitals are occupied in given configuration
  ! A Holmes, 21 Mar 2015

!#ifdef NUM_ORBITALS_GT_127
!  type(ik_vec),intent(in) :: det_up,det_dn
!  type(ik_vec) :: tmp_det
!#else
  !integer(ik),intent(in) :: det_up,det_dn
  !integer(ik) :: tmp_det
  integer(i16b), intent(in) :: det_up, det_dn
  integer(i16b) :: tmp_det
!#endif
  integer,intent(out) :: occ_up(:),occ_dn(:), nocc
  integer :: i_elec,i

    tmp_det = det_up
    i_elec = 0
    nocc = 0
    do while (tmp_det .ne. 0)
      i = trailz(tmp_det)+1
      i_elec = i_elec + 1
      occ_up(i_elec) = i
      if(i.gt.0) then 
          nocc = nocc + 1
      endif
      tmp_det = ibclr(tmp_det,i-1)
    enddo

    tmp_det = det_dn
    i_elec = 0
    do while (tmp_det .ne. 0)
      i = trailz(tmp_det)+1
      i_elec = i_elec + 1
      occ_dn(i_elec) = i
      if(i.gt.0) then 
          nocc = nocc + 1
      endif   
      tmp_det = ibclr(tmp_det,i-1)
    enddo

  end subroutine get_occ_orbs

  integer function Lz_unit(m_min, det)
      !det has one shell one spin
      integer, intent(in) :: m_min
      integer(i16b), intent(in) :: det
      integer(i16b) :: tmp_det
      integer :: i_elec, i
      Lz_unit = 0
      tmp_det = det
      i_elec = 0
      do while (tmp_det.ne.0)
        i = trailz(tmp_det) ! number of trailing 0s = index of the first 1(starting from 0) 
        Lz_unit = Lz_unit + m_min + i
        !i_elec = i_elec + 1
        tmp_det = ibclr(tmp_det, i)
      enddo
  end function Lz_unit

  integer function Lz_updn(m_min, m_max, det)
      !det_dn(m_max...m_min)--det_up(m_max...m_min)
      integer, intent(in) :: m_min, m_max
      integer :: halflen
      integer(i16b), intent(in) :: det
      integer(i16b) :: det_up, det_dn
      halflen = m_max - m_min + 1
      det_up = ibits(det, 0, halflen)
      det_dn = ibits(det, halflen, halflen)
      Lz_updn = Lz_unit(m_min, det_up) + Lz_unit(m_min, det_dn)
  end function Lz_updn

end module detgen
