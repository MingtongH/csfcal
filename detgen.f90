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

      recursive subroutine assign_shell(prdet_up,prdet_dn,prLz,pr2Sz,maxremLz,maxrem2Sz,desLz,des2Sz,configs,ncurr)
          !input maxrem includes all shells not assigned
          !including the one to be assigned in this recursion
          !only calculated by Lzmax(), no consideration of desLz
          implicit none
          integer, intent(in) :: prLz, pr2Sz, maxremLz, maxrem2Sz, ncurr
          ! ncurr = current shell number
          integer, intent(in) :: desLz, des2Sz
          integer, intent(in) :: configs(:, :)
          integer :: nshell_tot
          integer(i16b), intent(in) :: prdet_up, prdet_dn
          integer(i16b) :: curdet, tpdet, det_range!det_dn--det_up
          integer(i16b) :: detup, detdn
          integer :: curLz, cur2Sz, i, halflen
          integer :: curLzmax, curLzmin, m_max, m_min, cur2Szmax, cur2Szmin
          integer :: tpLz, tp2Sz, tpmaxremLz, tpmaxrem2Sz
          integer :: necur ! number of electrons in current shell
          nshell_tot = size(configs, 1)
          !write(*, '("nshell_tot = ", 1I5)') nshell_tot
          write(*, '(">>>>>>>>>>> Now at shell :",1I5)') ncurr
          write(*, '("prdet_up", B16)') prdet_up
          write(*, '("prdet_dn", B16)') prdet_dn
          write(*, '("prLz, pr2Sz, maxremLz, maxrem2Sz", 4I5)') prLz, pr2Sz, maxremLz, maxrem2Sz
          !**************************first shell(only need inputs as 0)
          if(ncurr.eq.0) then !1
              detup = 0
              detdn = 0
              call assign_shell(detup, detdn, 0, 0, Lzmax(configs), Smax_t2(configs), desLz, des2Sz, configs, ncurr+1)

          !*************************last shell
          else if(ncurr.eq.nshell_tot) then !if1
            necur = configs(ncurr, 3) 
            tpLz = desLz - prLz
            tp2Sz = des2Sz - pr2Sz
            tpmaxrem2Sz = Smax_t2(configs(ncurr:ncurr, 1:3)) 
            !No need of this line if this is calculated for all previous levels
            if(tpLz.gt.maxremLz.OR.tpLz.lt.(-maxremLz)) then!if2
                write(*, '("Cannot assign for des with Lz = ", 1I5)') tpLz
                return
            else if (tp2Sz.gt.tpmaxrem2Sz.OR.tp2Sz.lt.(-tpmaxrem2Sz)) then
                return
            else !if there is valid assignment
                !!construct and store complete det_up, prdet_dn:
                ! find all curdet that has curLz and cur2Sz
              write(*, '(" At last shell ")')
              necur = configs(ncurr, 3)
              m_max = configs(ncurr, 2) !lz 
              m_min = - configs(ncurr, 2) 
              !!>>>>  May need to change m_max, m_min for pruning
              halflen = m_max - m_min + 1
              det_range = 0
              do i = 2*halflen - necur, 2*halflen-1
                det_range = ibset(det_range, i)
              enddo
              !------------------------------------------------
              ! use just one combined det dn_det-up_det
              ! det_range is the last curdet with largest value
              !------------------------------------------------
             ! write(*, '(" Halflen = ", 1I5)') halflen
              write(*, '(" Searching det range 1 to ",B16)') det_range
              do curdet = 1, det_range !3
              !-------------------------------------------
              ! !det_dn(m_max...m_min)--det_up(m_max...m_min)
              ! !if occ_count = necur
              !    if Lz = curLz
              !        call next shell recursive functioni
              !-------------------------------------------
                 if (count_occ_orbs(curdet) .eq. necur) then !if 4
                     curLz = Lz_oneshell(m_min, halflen, curdet)
                     if(tpLz.eq.curLz) then
                         cur2Sz = Szt2_oneshell(halflen, curdet)
                         if(tp2Sz .eq. cur2Sz) then

                             write(*, '("*** Branch last shell, curLz =, cur2Sz = ", 3I5)') ncurr, curLz, cur2Sz
                             detup = prdet_up 
                             call detmod_shell(configs(ncurr, 1), configs(ncurr, 2), m_min, ibits(curdet, 0, halflen), detup)
                             detdn = prdet_dn
                             write(*, '("curdet = ", B16)') curdet
                             call  detmod_shell(configs(ncurr, 1), configs(ncurr, 2), m_min, ibits(curdet, halflen, halflen), detdn)
                             write(*, '("!!!Success!!! Generated detup:", B16)') detup 
                             write(*, '("!!!Success!!! Generated detdn:", B16)') detdn
                         endif
                     endif
                 endif !4
              enddo !end search 3

              return
              !break out of the recursion
            endif !2
          !*****************************recur
          else !1
              necur = configs(ncurr, 3)
              tpLz = Lzmax(configs(ncurr:ncurr, 1:3))
             ! write(*, '("Lzmax of this shell = ", 1I5)') tpLz
              tpmaxremLz = maxremLz - tpLz ! tpmaxremLz no longer include current shell
             ! write(*, '("maxremLz of next recur = ", 1I5)') tpmaxremLz
              curLzmax = min(tpLz, desLz - prLz + tpmaxremLz)
              curLzmin = max(-tpLz, desLz - prLz - tpmaxremLz)
              m_max = configs(ncurr, 2) !lz 
              m_min = - configs(ncurr, 2) 
              
              !!>>>>  May need to change m_max, m_min for pruning
              !!>>>>  check Sz beforehand to save time
              !cursmax_t2 = Smax_t2(configs(ncurr:ncurr, 1:3))
              !tpmaxrem2Sz = 1
              
              halflen = m_max - m_min + 1
              det_range = 0
              do i = 2*halflen - necur, 2*halflen-1
                det_range = ibset(det_range, i)
              enddo
              !------------------------------------------------
              ! use just one combined det dn_det-up_det
              ! det_range is the last curdet with largest value
              !------------------------------------------------
              write(*, '(" Searching det range 1 to ", B16)') det_range
              do curdet = 1, det_range
              !-------------------------------------------
              ! !det_dn(m_max...m_min)--det_up(m_max...m_min)
              ! !if occ_count = necur
              !    if Lz within range
              !        call next shell recursive functioni
              !-------------------------------------------
                 if (count_occ_orbs(curdet) .eq. necur) then
                     curLz = Lz_oneshell(m_min, halflen, curdet)
                     if(curLz.le.curLzmax .AND. curLz.ge.curLzmin) then
                         tpLz = prLz + curLz
                         cur2Sz = Szt2_oneshell(halflen, curdet)
                         tp2Sz = pr2Sz + cur2Sz
                         tpmaxrem2Sz = 0!Change later 
                         detup = prdet_up 
                         call detmod_shell(configs(ncurr, 1), configs(ncurr, 2), m_min, ibits(curdet, 0, halflen), detup)
                         detdn = prdet_dn
                         call  detmod_shell(configs(ncurr, 1), configs(ncurr, 2), m_min, ibits(curdet, halflen, halflen), detdn)
                         write(*, '("*** Branch shell#, curLz ", 2I5)') ncurr, curLz
                         call assign_shell(detup, detdn, tpLz, tp2Sz, tpmaxremLz, tpmaxrem2Sz, desLz, des2Sz, configs, ncurr+1)
                     endif
                 endif
              enddo

            !******************************end of recur
          endif!1
      end subroutine assign_shell


  subroutine detmod_shell(n, l, m_min, subdet, det)
      integer, intent(in) :: n, l, m_min
      integer(i16b), intent(in) :: subdet
      integer(i16b) :: det ! single spin det, inout
      integer(i16b) :: tpdet
      integer :: pos, i
      tpdet = subdet
      do while(tpdet.ne.0)
        i = trailz(tpdet)
        pos = locate_det(n, l, i + m_min)
        det = ibset(det, pos)
        tpdet = ibclr(tpdet, i)
      enddo
  end subroutine detmod_shell
                

  integer function locate_det(n, l, m)
      !Starting from 0
      integer, intent(in) :: n, l, m
      locate_det = (n-1)*n*(2*n-1)/6 + l**2 + m + l
  end function locate_det

  
  integer function Szt2_oneshell(halflen, det)
  !det_dn(m_max...m_min)--det_up(m_max...m_min)

      integer(i16b), intent(in) :: det
      integer(i16b) :: det_up, det_dn
      integer, intent(in) :: halflen
      det_up = ibits(det, 0, halflen)
      det_dn = ibits(det, halflen, halflen)
      Szt2_oneshell = count_occ_orbs(det_up) - count_occ_orbs(det_dn)
  end function Szt2_oneshell
          
  
  integer function count_occ_orbs(det)
      integer(i16b), intent(in) :: det
      integer(i16b) :: tmp_det
      integer :: i
    tmp_det = det
    !i_elec = 0
    count_occ_orbs = 0
    do while (tmp_det .ne. 0)
      i = trailz(tmp_det)
      count_occ_orbs = count_occ_orbs + 1
      tmp_det = ibclr(tmp_det,i)
    enddo
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
      integer :: i
      Lz_unit = 0
      tmp_det = det
      do while (tmp_det.ne.0)
        i = trailz(tmp_det) ! number of trailing 0s = index of the first 1(starting from 0) 
        Lz_unit = Lz_unit + m_min + i
        tmp_det = ibclr(tmp_det, i)
      enddo
  end function Lz_unit

  
  integer function Lz_oneshell(m_min, halflen, det)
      !det_dn(m_max...m_min)--det_up(m_max...m_min)
      integer, intent(in) :: m_min, halflen
      integer(i16b), intent(in) :: det
      integer(i16b) :: det_up, det_dn
      det_up = ibits(det, 0, halflen)
      det_dn = ibits(det, halflen, halflen)
      Lz_oneshell = Lz_unit(m_min, det_up) + Lz_unit(m_min, det_dn)
  end function Lz_oneshell
  
  
 ! character(LEN=16) function detdisplay(det)
 !     integer(i16b), intent(in) :: det
 !     integer(i16b) :: tmpdet
 !     integer :: i, t
!
 !     detdisplay= ''
 !     tmpdet = det
 !     do while (tmpdet.ne.0)
 !       i = trailz(tmpdet)
 !       tmpdet = ishft(tmpdet, -(i+1))
 !       do t = 1, i
 !         detdisplay = trim(adjustl(detdisplay))//char(0)
 !       enddo
 !       detdisplay = trim(adjustl(detdisplay))//char(1)
 !     enddo
 ! end function detdisplay 

end module detgen
