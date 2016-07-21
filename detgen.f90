module detgen
!Generate all determinants of Lz = 0 and Sz = S(from readinput) given one configuration
!such as 2s2, 2p5, 3d1
!Mingtong Han, June 27 2016
use prep, only: i16b, Lzmax, Smax_t2, ARRAY_START_LENGTH, DET_MAX_LENGTH
use, intrinsic :: iso_fortran_env, only: rk => real64
  implicit none
  contains
!>>>>>>> functions and subroutines
!isHalfFull(shellconfig), 
!assign_shell(prdet_up,prdet_dn,prLz,pr2Sz,maxremLz,maxrem2Sz, desLz,des2Sz,configs,ncurr, tot)
!detmod_shell(n, l, m_min, subdet, det), locate_det(n, l, m)
!Szt2_oneshell(halflen, det), count_occ_orbs(det)
!Lz_unit(m_min, det), Lz_oneshell(m_min, halflen, det)

!operator subroutines
!lpls, lms, eposinit, spls, sms
      subroutine spls(pos, detup, detdn, eposup, eposdn, coef, iszero)
          !s can only be 1/2, sz can only be 1/2 or -1/2
          !pos refers to the position of the electron in operation in detdn 
          !pos starts from 0
          integer, intent(in) :: pos
          integer(i16b) :: detup, detdn
          integer, allocatable :: eposup(:), eposdn(:)
          real(rk) :: coef, csq !Need initializion (1) at first run
          logical :: iszero !Need initializion(.false.) at first run
          if(.not.btest(detdn, pos)) then 
              return
          endif
          if(iszero) then
              return
          endif
          if(btest(detup, pos)) then
              !pos up is occupied
              coef = 0.0
              iszero = .true.
              return
          else
              write(*, '(15I4)') eposup(1:15)
              write(*, '(15I4)') eposdn(1:15)
              csq = 1 !1/2(1/2+1)-(-1/2)(-1/2 + 1)
              coef = coef*sqrt(real(csq))
              detup = ibset(detup, pos)
              detdn = ibclr(detdn, pos)
              eposup(pos+1) = eposdn(pos+1)
              eposdn(pos+1) = 0
          endif
      end subroutine spls

      subroutine sms(pos, detup, detdn, eposup, eposdn, coef, iszero)
          !s can only be 1/2, sz can only be 1/2 or -1/2
          !pos refers to the position of the electron in operation in detdn 
          !pos starts from 0
          integer, intent(in) :: pos
          integer(i16b) :: detup, detdn
          integer, allocatable :: eposup(:), eposdn(:)
          real(rk) :: coef, csq !Need initializion (1) at first run
          logical :: iszero !Need initializion(.false.) at first run
          if(.not.btest(detup, pos)) then 
              return
          endif
          if(iszero) then
              return
          endif
          if(btest(detdn, pos)) then
              !pos up is occupied
              coef = 0.0
              iszero = .true.
              return
          else
              write(*, '(15I4)') eposup(1:15)
              write(*, '(15I4)') eposdn(1:15)
              csq = 1 !1/2(1/2+1)-(-1/2)(-1/2 + 1)
              coef = coef*sqrt(real(csq))
              detdn = ibset(detdn, pos)
              detup = ibclr(detup, pos)
              eposdn(pos+1) = eposup(pos+1)
              eposup(pos+1) = 0
          endif
      end subroutine sms






      subroutine lpls(pos, det, epos, coef, iszero)
          !pos is the position of the electron that l+ is appled on, starting from 0
          integer, intent(in) :: pos
          integer(i16b) :: det
          integer, allocatable :: epos(:)
          real(rk) :: coef ! Should be initialized as 1 if used for the first time
          logical :: iszero ! Should be initialized as .false. if used for the first time
          integer :: n, l, lz, csq

          if(.not.btest(det, pos)) then
              return
          endif
          if(iszero) then 
              return
          endif
          call pos2nlm(pos, n, l, lz)
          csq = l*(l+1) - lz*(lz+1)
          if(csq.eq.0) then
              coef = 0.0
              iszero = .true.
              return
          elseif(btest(det, pos+1)) then
              coef = 0.0
              iszero = .true.
              return
          else
              !write(*, '(15I4)') epos(1:15)
              coef = sqrt(real(csq))*coef
              det = ibset(det, pos + 1)
              det = ibclr(det, pos)
              epos(pos+2) = epos(pos+1)
              epos(pos+1) = 0 !index of epos array starts from 1
             ! write(*, '(15I3)') epos(1:15)
             ! write(*, *) coef
          endif


      end subroutine lpls

      subroutine lms(pos, det, epos, coef, iszero)
          integer(i16b) :: det
          integer, allocatable :: epos(:) !should have been initialized
          real(rk) :: coef
          integer, intent(in) :: pos
          integer :: n, l, lz, csq
          logical :: iszero
          if(.not.btest(det, pos)) then
              return
          endif
          if(iszero) then 
              return
          endif
          call pos2nlm(pos, n, l, lz)
          csq = l*(l+1) - lz*(lz-1)
          if(csq.eq.0) then
              coef = 0.0
              iszero = .true.
              return
          elseif(btest(det, pos-1)) then
              coef = 0.0
              iszero = .true.
              return
          else
              !write(*, '(15I4)') epos(1:15)
              coef = sqrt(real(csq)) * coef
              det = ibset(det, pos - 1)
              det = ibclr(det, pos)
              epos(pos) = epos(pos+1)
              epos(pos+1) = 0 !index of epos array starts from 1
              !write(*, '(15I4)') epos(1:15)
          endif
      end subroutine lms


          

      subroutine eposinit(eposup, eposdn, detup, detdn)
          integer, allocatable :: eposup(:), eposdn(:)
          integer(i16b), intent(in) :: detup, detdn
          integer(i16b) :: tpdet
          integer :: i, l
          allocate(eposup(DET_MAX_LENGTH))
          allocate(eposdn(DET_MAX_LENGTH))
          do i = 1, DET_MAX_LENGTH
            eposup(i) = 0
            eposdn(i) = 0
          enddo
          
          l = 1
          tpdet = detup
          do while(tpdet.ne.0)
            i = trailz(tpdet) !Starting from 0
            eposup(i+1) = l
            l = l+1
            tpdet = ibclr(tpdet, i)
          enddo

          tpdet = detdn
          do while(tpdet.ne.0)
            i = trailz(tpdet)
            eposdn(i+1) = l
            l = l+1
            tpdet = ibclr(tpdet, i)
          enddo
      end subroutine eposinit


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

      recursive subroutine assign_shell(prdet_up,prdet_dn,prLz,pr2Sz,maxremLz, &
              & maxrem2Sz,desLz,des2Sz,configs,ncurr,tot, detlist)
          !input maxrem includes all shells not assigned
          !including the one
          !to be assigned in this recursion
          !only calculated by Lzmax(), no consideration of desLz
          implicit none
          integer, intent(in) :: prLz, pr2Sz, maxremLz, maxrem2Sz, ncurr
          ! ncurr = current shell number
          integer, intent(in) :: desLz, des2Sz
          integer, intent(in) :: configs(:, :)
          integer :: nshell_tot, tot
          integer(i16b), intent(in) :: prdet_up, prdet_dn
          integer(i16b) :: curdet, tpdet, det_range!det_dn--det_up
          integer(i16b) :: detup, detdn
          integer :: curLz, cur2Sz, i, halflen
          integer :: curLzmax, curLzmin, m_max, m_min, cur2Szmax, cur2Szmin
          integer :: tpLz, tp2Sz, tpmaxremLz, tpmaxrem2Sz
          integer :: necur ! number of electrons in current shell
          integer(i16b), allocatable :: detlist(:, :)
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
              if(.not.allocated(detlist)) then
                  allocate(detlist(ARRAY_START_LENGTH, 2))
              endif
              call assign_shell(detup, detdn, 0, 0, Lzmax(configs), Smax_t2(configs),&
                  & desLz, des2Sz, configs, ncurr+1, tot, detlist)

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

                             write(*, '("*** Branch last shell, curLz =, cur2Sz = ",&
                                 & 3I5)') ncurr, curLz, cur2Sz
                             detup = prdet_up 
                             call detmod_shell(configs(ncurr, 1), configs(ncurr, 2),&
                                 & m_min, ibits(curdet, 0, halflen), detup)
                             detdn = prdet_dn
                             write(*, '("curdet = ", B16)') curdet
                             call  detmod_shell(configs(ncurr, 1), configs(ncurr, 2),&
                                 & m_min, ibits(curdet, halflen, halflen), detdn)
                             write(*, '("!!!Success!!! Generated detup:", B16)') detup 
                             write(*, '("!!!Success!!! Generated detdn:", B16)') detdn
                             tot = tot + 1
                             detlist(tot, 1:2) = (/detup, detdn/)
                             write(*, *) tot
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
                         call detmod_shell(configs(ncurr, 1), configs(ncurr, 2),&
                             & m_min, ibits(curdet, 0, halflen), detup)
                         detdn = prdet_dn
                         call  detmod_shell(configs(ncurr, 1), configs(ncurr, 2),&
                             & m_min, ibits(curdet, halflen, halflen), detdn)
                         write(*, '("*** Branch shell#, curLz ", 2I5)') ncurr, curLz
                         call assign_shell(detup, detdn, tpLz, tp2Sz, tpmaxremLz,&
                             & tpmaxrem2Sz, desLz, des2Sz, configs, ncurr+1, tot, detlist)
                     endif
                 endif
              enddo

            !******************************end of recur
          endif!1
      end subroutine assign_shell


  subroutine detmod_shell(n, l, m_min, subdet, det)
      !modify the a full single spin determinant with a segment determinant of !
      !n, l starting with m_min on the very right position
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

  subroutine delocate(pos, nout, lout, mout)
      integer, intent(in) :: pos
      integer :: n, l, tp, nout, lout, mout
      do n = 1, 20
        if (pos.eq.(n-1)*n*(2*n-1)/6) then
            nout = n
            lout = 0
            mout = 0
            return
        elseif(pos.lt.(n-1)*n*(2*n-1)/6) then
            nout = n - 1
            tp = pos - (nout-1)*nout*(2*nout-1)/6
            do l = 0, 10
              if(tp.eq.l**2) then
                  lout = l
                  mout = -lout 
                  return
              elseif(tp.lt.l**2) then
                  lout = l - 1
                  mout = tp - lout**2 - lout
                  return
              endif
            enddo
        endif
      enddo
      
  end subroutine

  subroutine pos2nlm(pos, n, l, m)
      integer, intent(in) :: pos
      integer :: n, l, m
      select case(pos)
        case(0)
            n = 1
            l = 0
            m = 0
        case(1)
            n = 2
            l = 0
            m = 0
        case default
            call delocate(pos, n, l, m)
      end select
  end subroutine pos2nlm

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

!This subroutine is not used 
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
      !det contains only  one subshell one spin, with the first(rightmost) digit for m_min
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
      !det = det_dn(m_max...m_min)--det_up(m_max...m_min)
      integer, intent(in) :: m_min, halflen
      integer(i16b), intent(in) :: det
      integer(i16b) :: det_up, det_dn
      det_up = ibits(det, 0, halflen)
      det_dn = ibits(det, halflen, halflen)
      Lz_oneshell = Lz_unit(m_min, det_up) + Lz_unit(m_min, det_dn)
  end function Lz_oneshell
  
 integer function Lz_det(det)
     integer(i16b), intent(in) :: det
     integer(i16b) :: tpdet 
     integer :: i, n, l, lz
     tpdet = det
     Lz_det = 0
     do while (tpdet.ne.0)
       i = trailz(tpdet)
       call delocate(i, n, l, lz)
       Lz_det = Lz_det + lz
       tpdet = ibclr(tpdet, i)
     end do
 end function Lz_det



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
