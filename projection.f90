module projection
      !Generate csf by projection
      use prep, only: i16b, ARRAY_START_LENGTH, ARRAY_SHORT_LENGTH, DET_MAX_LENGTH
      use, intrinsic :: iso_fortran_env, only: rk => real64
      use detgen, only: lpls, lms, eposinit, spls, sms

      implicit none
      
      real(rk), parameter :: REAL_MIN = 1.0e-31

      contains
      logical function equals0(coef)
          real(rk), intent(in) :: coef
          if(abs(coef).le.REAL_MIN) then
              equals0 = .true.
          else
              equals0 = .false.
          endif
      end function equals0

      subroutine initlists(det, basis, coefs, eposes, iszeros, num)
          integer(i16b), intent(in) :: det(:)
          integer(i16b), allocatable :: basis(:, :)
          real(rk), allocatable :: coefs(:)
          integer, allocatable :: eposes(:, :, :)
          logical, allocatable :: iszeros(:)
          integer :: i, num
          if(.not.allocated(basis)) then
              allocate(basis(ARRAY_SHORT_LENGTH, 2))
          endif
          if(.not.allocated(coefs)) then
              allocate(coefs(ARRAY_SHORT_LENGTH))
          endif
          if(.not.allocated(eposes)) then
              allocate(eposes(ARRAY_SHORT_LENGTH, 2, DET_MAX_LENGTH))
          endif
          if(.not.allocated(iszeros)) then
              allocate(iszeros(ARRAY_SHORT_LENGTH))
          endif
          basis(1, 1:2) = det(1:2)
          call eposinit(eposes(1, 1, :), eposes(1, 2, :), det(1), det(2))
          coefs(1) = 1.0
          iszeros(1) = .false.
          num = 1
      end subroutine initlists


      subroutine Proj_L(Lmin, Lmax, Ldes, &
              & inbasis, incoefs, inepos, iniszeros, innum,  &
              & basislist, coeflist, eposlist, iszerolist, num)
          integer(i16b), intent(in) :: inbasis(:, :)
          real(rk), intent(in) :: incoefs(:)
          integer, intent(in) :: inepos(:, :, :)
          logical, intent(in) :: iniszeros(:)
          integer, intent(in) :: Lmax, Lmin, Ldes, innum
          integer(i16b), allocatable :: basislist(:, :), tpbasis(:, :), basislist1(:, :)
          real(rk), allocatable :: coeflist(:), tpcoefs(:), coeflist1(:) 
          integer, allocatable :: eposlist(:, :, :), tpepos(:, :, :), eposlist1(:, :, :)
          logical, allocatable :: iszerolist(:), tpiszeros(:), iszerolist1(:)
          integer :: num, tpnum, i, Lp, num1
          real(rk) :: tail
          write(*, *) '>>>>>>>>>>>>>>>>>>>>>>> Projection P_L '
          allocate(tpbasis(ARRAY_START_LENGTH, 2))
          allocate(tpcoefs(ARRAY_START_LENGTH))
          allocate(tpepos(ARRAY_START_LENGTH, 2, DET_MAX_LENGTH))
          allocate(tpiszeros(ARRAY_START_LENGTH))
          !----------Initialize tp.. as in..
          do i = 1, innum
              tpbasis(i, 1:2) = inbasis(i, 1:2)
              tpcoefs(i) = incoefs(i)
              tpepos(i, 1, :) = inepos(i, 1, :)
              tpepos(i, 2, :) = inepos(i, 2, :)
              tpiszeros(i) = iniszeros(i)
          enddo
          tpnum = innum
          !---------- Apply operators for each Lp -----------
          do Lp = Lmin, Lmax
              if(Lp.eq.Ldes) then
                  continue
              endif

              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              !TODO doesn't include Lz^2 + Lz yet
              !implement Lz first then change PL
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              
              !L^2 = L-L+ + Lz^2 + Lz
              write(*, *) ' '
              write(*,  '("====================== Now at Lp =", 1I3)') Lp
              !-------set intermediates (basislist1...) and output ( basislist ... ) empty
              num1 = 0
              num = 0
              !----------L^2 = L-L+ + Lz^2 + Lz
              call Lplus_multiple(tpbasis, tpcoefs, tpepos, tpiszeros, tpnum, &
                  & basislist1, coeflist1, eposlist1, iszerolist1, num1)
              call Lminus_multiple(basislist1, coeflist1, eposlist1, iszerolist1, num1, &
                  & basislist, coeflist, eposlist, iszerolist, num)
              !TODO When adding Lz^2 + Lz
              !Lminus_multiple output basislist2 etc
              !Lz^2+Lz takes in basislist2 etc, output basislist etc

              !-------------------(- Lp(Lp+1))---------------------
              write(*, *) ' ------------------- Append tail to basislist '
              tail = real(-Lp * (Lp + 1))
              write(*, '("tail = ", 1f15.10)') tail            
              call scalar_append(tail, tpbasis, tpcoefs, tpepos, tpiszeros, tpnum, &
                  & basislist, coeflist, eposlist, iszerolist, num)
              write(*, '("Total number of dets in basislist", 1I4)') num

              !----------------------------tp.. <- output ( basislist, ...)
              tpnum = 0
              do i = 1, num
                  tpbasis(i, 1:2) = basislist(i, 1:2)
                  tpcoefs(i) = coeflist(i)
                  tpepos(i, 1, :) = eposlist(i, 1, :)
                  tpepos(i, 2, :) = eposlist(i, 2, :)
                  tpiszeros(i) = iszerolist(i)
              enddo
              tpnum = num
          enddo
          deallocate(tpbasis)
          deallocate(tpcoefs)
          deallocate(tpepos)
          deallocate(tpiszeros)

      end subroutine Proj_L


      subroutine normalizetable(coeftable, ndets, ncsf)
          !call this after coeftable is complete
          real(rk) :: coeftable(:, :), sumcol
          integer, intent(in) :: ndets, ncsf
          integer :: i, j
          write(*, *) '-------- Normalize coef table ------------'
          do i = 1, ncsf
            sumcol = 0
            do j = 1, ndets
              sumcol = sumcol + coeftable(j, i)**2
            enddo
            sumcol = sqrt(sumcol)
            do j = 1, ndets
              coeftable(j, i) = coeftable(j, i)/sumcol
            enddo
          enddo
      end subroutine normalizetable

      subroutine getallsigns(basislist, coeflist, eposlist, iszerolist, num)
          integer(i16b), intent(in) :: basislist(:, :)
          real(rk) :: coeflist(:)
          integer, intent(in) :: eposlist(:, :, :), num
          logical, intent(in) :: iszerolist(:)
          integer :: i

          write(*, *) '>>>>>>>>>>>>>>>>>>>> Get all signs of basislist <<<<<<<<<<<<<<<<<<<'
          do i = 1, num
              if(.not.iszerolist(i)) then
                  call getsign(basislist(i, 1), basislist(i, 2), eposlist(i, 1, :),&
                      & eposlist(i, 2, :), coeflist(i))
              endif
          enddo
      end subroutine getallsigns


      subroutine getsign(detup, detdn, eposup, eposdn, coef)
          integer(i16b), intent(in) :: detup, detdn
          integer, intent(in) :: eposup(:), eposdn(:)
          integer, allocatable :: tparray(:)
          real(rk) :: coef
          integer(i16b) :: tpdet
          integer :: num, i, swaps
          write(*, *) 'Getsign()'
          write(*, '(2b16)') detup, detdn
          write(*, '(15I3)') eposup(1:15)
          write(*, '(15I3)') eposdn(1:15)

          allocate(tparray(2*size(eposup)))
          tpdet = detup
          num = 0
          do while (tpdet.ne.0)
              i = trailz(tpdet)
              num = num + 1
              tparray(num) = eposup(i+1)
              tpdet = ibclr(tpdet, i)
          enddo
          tpdet = detdn
          do while (tpdet.ne.0)
              i = trailz(tpdet)
              num = num + 1
              tparray(num) = eposdn(i+1)
              tpdet = ibclr(tpdet, i)
          enddo
          swaps = countswps(tparray, num)
          i = mod(swaps, 2)
          if(i.ne.0) then
              coef = - coef
          endif
          write(*, '("Sign:", 1I4)') i
          deallocate(tparray)
      end subroutine getsign
      
      
      integer function  countswps(array, num)
            integer, allocatable:: array(:)
            integer:: num, i, tp
            countswps = 0
            write(*, '("Start array:")')
            write(*, *) array(1:num)
            do i = 1, num
                do while (i.ne.array(i))
                    tp = array(i)
                    array(i) = array(tp)
                    array(tp) = tp
                    countswps = countswps + 1
                enddo
            enddo
            write(*, '("End array:, number of swps")')
            write(*, *) array(1:num)
            write(*, '("Count swps = ", 1I4)') countswps
      end function countswps


      subroutine collect_csf(basislist, coeflist, iszerolist, num,&
              & allbasis, coeftable, ndets, ncsf)
          !add a csf to coeftable
          !num = dimension of the basislist
          !ndets, ncsf are coeftable sizes
          !coeflist after calling getsign
          integer(i16b), intent(in):: basislist(:, :)
          real(rk), intent(in) :: coeflist(:)
          logical, intent(in) :: iszerolist(:)
          integer, intent(in) :: num
          integer(i16b), allocatable :: allbasis(:, :)
          real(rk), allocatable :: coeftable(:, :)!ndets, ncsf
          integer :: ndets, ncsf, i,j, tp
          write(*, *) '>>>>>>>>>>>>>>>>>>>>>>>>> Collect one csf <<<<<<<<<<<<<<<<<<<<<<<<'

          !Firsty entry initialize
          if(ndets.eq.0 .OR. ncsf.eq.0) then
              if(.not.allocated(allbasis)) then
                  allocate(allbasis(ARRAY_START_LENGTH, 2))
                  write(*, '("Allocated allbasis")')
              endif
              if(.not.allocated(coeftable)) then
                  allocate(coeftable(ARRAY_START_LENGTH, ARRAY_SHORT_LENGTH))
                  write(*, '("Allocated coeftable")')
              endif
          endif
          !Initialize new column to be 0, ndets 0s
          if(ndets.ne. 0) then
              do i = 1, ndets
                coeftable(i, ncsf+1) = .0
              enddo
              write(*, '("Initialized new column to 0")')
          endif

          do i = 1, num
              write(*, '("Processing det# from list:", 1I3)') i
              if(iszerolist(i)) then
                  cycle
              endif

              tp = findindetlist(basislist(i, :), allbasis(1:ndets, :), ndets)
              write(*, '("findindetlist()=", 1I4)') tp
              !If not in the list, add new det
              !! now this is only appending, may need to order dets later TODO
              if(tp.eq.0) then
                  allbasis(ndets + 1, :) = basislist(i, :)
                  do j = 1, ncsf
                    coeftable(ndets + 1, j) = .0
                  enddo
                  write(*, *) 'Initialized new row'
                  coeftable(ndets +1, ncsf+1) = coeflist(i)
                  ndets = ndets + 1
                  write(*, '("Added new det row", 1I5)') ndets
              else!If in the list, add coef
                  coeftable(tp, ncsf+1) = coeftable(tp, ncsf+1) + coeflist(i)
                  write(*, '("Combined coef with existing det")')
              endif
          enddo
          ncsf = ncsf + 1
          write(*, '("--------------- Collected csf #:", 1I3)') ncsf
          write(*, '("Basis:      detup         detdn     coefs of new csf")')
          do i = 1, ndets
             write(*, '(2b16, 1f18.15)') allbasis(ncsf, 1:2), coeftable(i, ncsf)
          enddo
          !--------------------------------------
          !if all coefs in the coef list are zero
          !ncsf will not be incremented, so even if the new column is initialized 
          !it will still be overwritten by the next loop, or accessible if no next loop 
      end subroutine collect_csf


      integer function findindetlist(det, allbasis, ndets)
          integer(i16b), intent(in) :: det(:), allbasis(:, :)
          integer, intent(in) :: ndets
          integer :: i
          do i = 1, ndets
            if(det(1).eq.allbasis(i, 1)) then
                if(det(2).eq.allbasis(i, 2)) then
                    findindetlist = i
                    return
                endif
            endif
          enddo
          findindetlist = 0
      end function findindetlist


      subroutine scalar_append(const, inbasis, incoefs, inepos, iniszeros, innum, &
              & basislist, coeflist, eposlist, iszerolist, num)
          !append to output list
          integer(i16b), intent(in) :: inbasis(:, :)
          real(rk),  intent(in) :: incoefs(:), const
          integer, intent(in) :: inepos(:, :, :)
          logical, intent(in) :: iniszeros(:)
          integer, intent(in) :: innum
          integer(i16b), allocatable :: basislist(:, :)
          real(rk), allocatable :: coeflist(:)
          integer, allocatable :: eposlist(:, :, :)
          logical, allocatable :: iszerolist(:)
          integer :: num, i, tp
          if(equals0(const)) then
              write(*, *) '0.0, exit append'
              return
          endif
          if(num.eq.0) then
              if(.not.allocated(basislist)) then 
                  allocate(basislist(ARRAY_START_LENGTH, 2))
              endif 
              if(.not.allocated(coeflist)) then
                  allocate(coeflist(ARRAY_START_LENGTH))
              endif
              if(.not.allocated(eposlist)) then
                  allocate(eposlist(ARRAY_START_LENGTH, 2, DET_MAX_LENGTH))
              endif
              if(.not.allocated(iszerolist)) then
                  allocate(iszerolist(ARRAY_START_LENGTH))
              endif
          endif

          do i = 1, innum
            write(*, '("At input basis det #", 1I4)') i
            if(iniszeros(i)) then
                write(*, '("iszero, onto next det")')
                cycle
            endif
            !combine if same basis
            !tp = findindetlist(inbasis(i, 1:2), basislist(1:num, :), num)
            !write(*, '("find in inbasis = ", 1I4)') tp
            !if(tp.ne.0) then
            !    !TODO check epos equal? 
            !    coeflist(tp) = incoefs(tp) * const
            !    write(*, *) 'det in basislist, changed coef'
            !    cycle
            !endif
            num = num + 1
            basislist(num, 1:2) = inbasis(i, 1:2)
            coeflist(num) = incoefs(i)*const
            eposlist(num, 1, :) = inepos(i, 1, :)
            eposlist(num, 2, :) = inepos(i, 2, :)
            iszerolist(num) = iszerolist(i)
            write(*, '("Appended det #", 1I4)') num
          end do
      end subroutine scalar_append


      subroutine Sminus_multiple(inbasis, incoefs, inepos, iniszeros, innum, &
              & basislist, coeflist, eposlist, iszerolist, num)
          !append to output list
          integer(i16b), intent(in) :: inbasis(:, :)
          real(rk),  intent(in) :: incoefs(:)
          integer, intent(in) :: inepos(:, :, :)
          logical, intent(in) :: iniszeros(:)
          integer, intent(in) :: innum
          integer(i16b), allocatable :: basislist(:, :)
          real(rk), allocatable :: coeflist(:)
          integer, allocatable :: eposlist(:, :, :)
          logical, allocatable :: iszerolist(:)
          integer :: num, i

          if(num.eq.0) then
              if(.not.allocated(basislist)) then 
                  allocate(basislist(ARRAY_START_LENGTH, 2))
              endif 
              if(.not.allocated(coeflist)) then
                  allocate(coeflist(ARRAY_START_LENGTH))
              endif
              if(.not.allocated(eposlist)) then
                  allocate(eposlist(ARRAY_START_LENGTH, 2, DET_MAX_LENGTH))
              endif
              if(.not.allocated(iszerolist)) then
                  allocate(iszerolist(ARRAY_START_LENGTH))
              endif
          endif
          do i = 1, innum
            write(*, *)
            write(*, '("=================== Applying S- ta multiple dets, now at No.", 1I5)') i 
            call Sminus_single(inbasis(i,1), inbasis(i, 2), inepos(i, 1, :), &
                &inepos(i, 2, :), incoefs(i), iniszeros(i), &
                & basislist, coeflist, eposlist, iszerolist, num)
          enddo
      end subroutine Sminus_multiple
      
      subroutine Splus_multiple(inbasis, incoefs, inepos, iniszeros, innum, &
              & basislist, coeflist, eposlist, iszerolist, num)
          !append to output list
          integer(i16b), intent(in) :: inbasis(:, :)
          real(rk), intent(in) :: incoefs(:)
          integer, intent(in) :: inepos(:, :, :)
          logical, intent(in) :: iniszeros(:)
          integer, intent(in) :: innum
          integer(i16b), allocatable :: basislist(:, :)
          real(rk), allocatable :: coeflist(:)
          integer, allocatable :: eposlist(:, :, :)
          logical, allocatable :: iszerolist(:)
          integer :: num, i

          if(num.eq.0) then
              if(.not.allocated(basislist)) then 
                  allocate(basislist(ARRAY_START_LENGTH, 2))
              endif 
              if(.not.allocated(coeflist)) then
                  allocate(coeflist(ARRAY_START_LENGTH))
              endif
              if(.not.allocated(eposlist)) then
                  allocate(eposlist(ARRAY_START_LENGTH, 2, DET_MAX_LENGTH))
              endif
              if(.not.allocated(iszerolist)) then
                  allocate(iszerolist(ARRAY_START_LENGTH))
              endif
          endif
          do i = 1, innum
            write(*, *)
            write(*, '("--------------- Applying S+ to multiple dets, now at No.", 1I5)') i 
            call Splus_single(inbasis(i,1), inbasis(i, 2), inepos(i, 1, :), &
                &inepos(i, 2, :), incoefs(i), iniszeros(i), &
                & basislist, coeflist, eposlist, iszerolist, num)
          enddo
      end subroutine Splus_multiple



      subroutine Lminus_multiple(inbasis, incoefs, inepos, iniszeros, innum, &
              & basislist, coeflist, eposlist, iszerolist, num)
          !append to output list
          integer(i16b), intent(in) :: inbasis(:, :)
          real(rk), intent(in) :: incoefs(:)
          integer, intent(in) :: inepos(:, :, :)
          logical, intent(in) :: iniszeros(:)
          integer, intent(in) :: innum
          integer(i16b), allocatable :: basislist(:, :)
          real(rk), allocatable :: coeflist(:)
          integer, allocatable :: eposlist(:, :, :)
          logical, allocatable :: iszerolist(:)
          integer :: num, i

          if(num.eq.0) then
              if(.not.allocated(basislist)) then 
                  allocate(basislist(ARRAY_START_LENGTH, 2))
              endif 
              if(.not.allocated(coeflist)) then
                  allocate(coeflist(ARRAY_START_LENGTH))
              endif
              if(.not.allocated(eposlist)) then
                  allocate(eposlist(ARRAY_START_LENGTH, 2, DET_MAX_LENGTH))
              endif
              if(.not.allocated(iszerolist)) then
                  allocate(iszerolist(ARRAY_START_LENGTH))
              endif
          endif
          do i = 1, innum
            write(*, *)
            write(*, '("---------------------- Applying L- to multiple dets, now at No.", 1I5)') i 
            call Lminus_single(inbasis(i,1), inbasis(i, 2), inepos(i, 1, :), &
                &inepos(i, 2, :), incoefs(i), iniszeros(i), &
                & basislist, coeflist, eposlist, iszerolist, num)
          enddo
      end subroutine Lminus_multiple

      subroutine Lplus_multiple(inbasis, incoefs, inepos, iniszeros, innum, &
              & basislist, coeflist, eposlist, iszerolist, num)
          !append to outaut list
          integer(i16b), intent(in) :: inbasis(:, :)
          real(rk),  intent(in) :: incoefs(:)
          integer, intent(in) :: inepos(:, :, :)
          logical, intent(in) :: iniszeros(:)
          integer, intent(in) :: innum
          integer(i16b), allocatable :: basislist(:, :)
          real(rk), allocatable :: coeflist(:)
          integer, allocatable :: eposlist(:, :, :)
          logical, allocatable :: iszerolist(:)
          integer :: num, i

          if(num.eq.0) then
              if(.not.allocated(basislist)) then 
                  allocate(basislist(ARRAY_START_LENGTH, 2))
              endif 
              if(.not.allocated(coeflist)) then
                  allocate(coeflist(ARRAY_START_LENGTH))
              endif
              if(.not.allocated(eposlist)) then
                  allocate(eposlist(ARRAY_START_LENGTH, 2, DET_MAX_LENGTH))
              endif
              if(.not.allocated(iszerolist)) then
                  allocate(iszerolist(ARRAY_START_LENGTH))
              endif
          endif
          do i = 1, innum
            write(*, *)
            write(*, '("-------------------- Applying L+ to multiple dets, now at No.", 1I5)') i 
            call Lplus_single(inbasis(i,1), inbasis(i, 2), inepos(i, 1, :), &
                &inepos(i, 2, :), incoefs(i), iniszeros(i), &
                & basislist, coeflist, eposlist, iszerolist, num)
          enddo
      end subroutine Lplus_multiple
      
      subroutine Sminus_single(detup, detdn, eposup, eposdn, coef, iszero, &
              & basislist, coeflist, eposlist, iszerolist, num)
          integer(i16b), intent(in) :: detup, detdn
          integer, intent(in) :: eposup(:), eposdn(:)
          real(rk), intent(in) :: coef
          logical, intent(in) :: iszero
          logical, allocatable :: iszerolist(:)
          integer(i16b), allocatable :: basislist(:, :)
          real(rk), allocatable :: coeflist(:)
          integer, allocatable :: eposlist(:, :, :)!#, up/down, electron position
          integer :: num !number of basis determinants
          integer(i16b) :: tpdet, outdetup, outdetdn
          integer :: i
          integer, allocatable :: outeposup(:), outeposdn(:)
          logical :: outis0
          real(rk) :: outcoef
          write(*, '("----------------------- Apply S- on detup, detdn: ", 2B16)') detup, detdn
          write(*, '(15I3)') eposup(1:15)
          write(*, '(15I3)') eposdn(1:15)
          write(*, *) coef, iszero, num

          if(num.eq.0) then
              if(.not.allocated(basislist)) then 
                  allocate(basislist(ARRAY_START_LENGTH, 2))
              endif 
              if(.not.allocated(coeflist)) then
                  allocate(coeflist(ARRAY_START_LENGTH))
              endif
              if(.not.allocated(eposlist)) then
                  allocate(eposlist(ARRAY_START_LENGTH, 2, DET_MAX_LENGTH))
              endif
              if(.not.allocated(iszerolist)) then
                  allocate(iszerolist(ARRAY_START_LENGTH))
              endif
          endif
          allocate(outeposup(DET_MAX_LENGTH))
          allocate(outeposdn(DET_MAX_LENGTH))
          !apply s- on detup
          tpdet = detup
          do while(tpdet.ne.0)
              i = trailz(tpdet) ! starting from 0
              outdetup = detup
              outdetdn = detdn
              outeposup(:) = eposup(:)
              outeposdn(:) = eposdn(:)
              outcoef = coef
              outis0 = iszero.or.equals0(outcoef)
              call sms(i, outdetup, outdetdn, outeposup, outeposdn, outcoef, outis0)
              !outdet outepos outcoef outis0 are all updated here
              !Original det is not added to the list here when outis0 = true
              write(*, '("s- on e at up pos:", 1I4)') i
              write(*, *) outis0
              if(.not.outis0) then
                  num = num + 1
                  basislist(num, 1:2) = (/outdetup, outdetdn/)
                  coeflist(num) = outcoef
                  eposlist(num, 1, :) = outeposup(:)
                  eposlist(num, 2, :) = outeposdn(:)
                  iszerolist(num) = outis0
                  write(*, '("Output detup, detdn: ", 2B16)') basislist(num, 1:2)
                  write(*, '(15I3)') eposlist(num, 1, 1:15)
                  write(*, '(15I3)') eposlist(num, 2, 1:15)
                  write(*, *) coeflist(num), iszerolist(num), num
              endif

              tpdet = ibclr(tpdet, i) 
          enddo
          deallocate(outeposup)
          deallocate(outeposdn)
      end subroutine Sminus_single

      
      subroutine Splus_single(detup, detdn, eposup, eposdn, coef, iszero, &
              & basislist, coeflist, eposlist, iszerolist, num)
          integer(i16b), intent(in) :: detup, detdn
          integer, intent(in) :: eposup(:), eposdn(:)
          real(rk), intent(in) :: coef
          logical, intent(in) :: iszero
          logical, allocatable :: iszerolist(:)
          integer(i16b), allocatable :: basislist(:, :)
          real(rk), allocatable :: coeflist(:)
          integer, allocatable :: eposlist(:, :, :)!#, up/down, electron position
          integer :: num !number of basis determinants
          integer(i16b) :: tpdet, outdetup, outdetdn
          integer :: i
          integer, allocatable :: outeposup(:), outeposdn(:)
          logical :: outis0
          real(rk) :: outcoef
          write(*, *) ' '
          write(*, '("--------------------- Apply S+ on detup, detdn: ", 2B16)') detup, detdn
          write(*, '(15I3)') eposup(1:15)
          write(*, '(15I3)') eposdn(1:15)
          write(*, *) coef, iszero, num

          if(num.eq.0) then
              if(.not.allocated(basislist)) then 
                  allocate(basislist(ARRAY_START_LENGTH, 2))
              endif 
              if(.not.allocated(coeflist)) then
                  allocate(coeflist(ARRAY_START_LENGTH))
              endif
              if(.not.allocated(eposlist)) then
                  allocate(eposlist(ARRAY_START_LENGTH, 2, DET_MAX_LENGTH))
              endif
              if(.not.allocated(iszerolist)) then
                  allocate(iszerolist(ARRAY_START_LENGTH))
              endif
          endif
          allocate(outeposup(DET_MAX_LENGTH))
          allocate(outeposdn(DET_MAX_LENGTH))
          !apply s+ on detdn
          tpdet = detdn
          do while(tpdet.ne.0)
              i = trailz(tpdet) ! starting from 0
              outdetup = detup
              outdetdn = detdn
              outeposup(:) = eposup(:)
              outeposdn(:) = eposdn(:)
              outcoef = coef
              outis0 = iszero.or.equals0(outcoef)
              call spls(i, outdetup, outdetdn, outeposup, outeposdn, outcoef, outis0)
              !outdet outepos outcoef outis0 are all updated here
              !Original det is not added to the list here when outis0 = true
              write(*, '("s+ on e at up pos:", 1I4)') i
              write(*, *) outis0
              if(.not.outis0) then
                  num = num + 1
                  basislist(num, 1:2) = (/outdetup, outdetdn/)
                  coeflist(num) = outcoef
                  eposlist(num, 1, :) = outeposup(:)
                  eposlist(num, 2, :) = outeposdn(:)
                  iszerolist(num) = outis0
                  write(*, '("Output detup, detdn: ", 2B16)') basislist(num, 1:2)
                  write(*, '(15I3)') eposlist(num, 1, 1:15)
                  write(*, '(15I3)') eposlist(num, 2, 1:15)
                  write(*, *) coeflist(num), iszerolist(num), num
              endif

              tpdet = ibclr(tpdet, i) 
          enddo
          deallocate(outeposup)
          deallocate(outeposdn)
      end subroutine Splus_single


 


          !TODO add original det in L+single, L-single? may not need to
      subroutine Lminus_single(detup, detdn, eposup, eposdn, coef, iszero, &
              & basislist, coeflist, eposlist, iszerolist, num )
          integer(i16b), intent(in) :: detup, detdn
          integer, intent(in) :: eposup(:), eposdn(:)
          real(rk), intent(in) :: coef
          logical, intent(in) :: iszero
          logical, allocatable :: iszerolist(:)
          integer(i16b), allocatable :: basislist(:, :)
          real(rk), allocatable :: coeflist(:)
          integer, allocatable :: eposlist(:, :, :)!#, up/down, electron position
          integer :: num !number of basis determinants
          integer(i16b) :: tpdet, outdet
          integer :: i
          integer, allocatable :: outepos(:)
          logical :: outis0
          real(rk) :: outcoef
          write(*, '("---------------------- Apply L- on detup, detdn: ", 2B16)') detup, detdn
          write(*, '(15I3)') eposup(1:15)
          write(*, '(15I3)') eposdn(1:15)
          write(*, *) coef, iszero, num

          if(num.eq.0) then
              if(.not.allocated(basislist)) then 
                  allocate(basislist(ARRAY_START_LENGTH, 2))
              endif 
              if(.not.allocated(coeflist)) then
                  allocate(coeflist(ARRAY_START_LENGTH))
              endif
              if(.not.allocated(eposlist)) then
                  allocate(eposlist(ARRAY_START_LENGTH, 2, DET_MAX_LENGTH))
              endif
              if(.not.allocated(iszerolist)) then
                  allocate(iszerolist(ARRAY_START_LENGTH))
              endif
          endif
          allocate(outepos(DET_MAX_LENGTH))
          !apply l- on detup
          tpdet = detup
          do while(tpdet.ne.0)
              i = trailz(tpdet) ! starting from 0
              outdet = detup
              outepos(:) = eposup(:)
              outcoef = coef
              outis0 = iszero.or.equals0(outcoef)
              call lms(i, outdet, outepos, outcoef, outis0)
              !outdet outepos outcoef outis0 are all updated here
              !Original det is not added to the list here when outis0 = true
              write(*, '("l- on e at up pos:", 1I4)') i
              write(*, *) outis0
              if(.not.outis0) then
                  num = num + 1
                  basislist(num, 1:2) = (/outdet, detdn/)
                  coeflist(num) = outcoef
                  eposlist(num, 1, :) = outepos(:)
                  eposlist(num, 2, :) = eposdn(:)
                  iszerolist(num) = outis0
                  write(*, '("Output detup, detdn: ", 2B16)') basislist(num, 1:2)
                  write(*, '(15I3)') eposlist(num, 1, 1:15)
                  write(*, '(15I3)') eposlist(num, 2, 1:15)
                  write(*, *) coeflist(num), iszerolist(num), num
              endif

              tpdet = ibclr(tpdet, i) 
          enddo
           !apply l+ on detdn
          tpdet = detdn
          do while(tpdet.ne.0)
              i = trailz(tpdet) ! starting from 0
              outdet = detdn
              outepos(:) = eposdn(:)
              outcoef = coef
              outis0 = iszero.or.equals0(coef)
              call lms(i, outdet, outepos, outcoef, outis0)
              !Original det is not added to the list here when outis0 = true
              write(*, '("l- on e at dn pos:", 1I4)') i
              write(*, *) outis0
              if(.not.outis0) then
                  num = num + 1
                  basislist(num, 1:2) = (/detup, outdet/)
                  coeflist(num) = outcoef
                  eposlist(num, 1, :) = eposup(:)
                  eposlist(num, 2, :) = outepos(:)
                  iszerolist(num) = outis0
                  write(*, '("Output detup, detdn: ", 2B16)') basislist(num, 1:2)
                  write(*, '(15I3)') eposlist(num, 1, 1:15)
                  write(*, '(15I3)') eposlist(num, 2, 1:15)
                  write(*, *) coeflist(num), iszerolist(num), num
              endif

              tpdet = ibclr(tpdet, i) 
          enddo
         deallocate(outepos)

     end subroutine Lminus_single

      subroutine Lplus_single(detup, detdn, eposup, eposdn, coef, iszero, &
              & basislist, coeflist, eposlist, iszerolist, num )
          integer(i16b), intent(in) :: detup, detdn
          integer, intent(in) :: eposup(:), eposdn(:)
          real(rk), intent(in) :: coef
          logical, intent(in) :: iszero
          logical, allocatable :: iszerolist(:)
          integer(i16b), allocatable :: basislist(:, :)
          real(rk), allocatable :: coeflist(:)
          integer, allocatable :: eposlist(:, :, :)!#, up/down, electron position
          integer :: num !number of basis determinants
          integer(i16b) :: tpdet, outdet
          integer :: i
          integer, allocatable :: outepos(:)
          logical :: outis0
          real(rk) :: outcoef
          write(*, '("---------------------- Apply L+ on detup, detdn: ", 2B16)') detup, detdn
          write(*, '(15I3)') eposup(1:15)
          write(*, '(15I3)') eposdn(1:15)
          write(*, *) coef, iszero, num

          if(num.eq.0) then
              if(.not.allocated(basislist)) then 
                  allocate(basislist(ARRAY_START_LENGTH, 2))
              endif 
              if(.not.allocated(coeflist)) then
                  allocate(coeflist(ARRAY_START_LENGTH))
              endif
              if(.not.allocated(eposlist)) then
                  allocate(eposlist(ARRAY_START_LENGTH, 2, DET_MAX_LENGTH))
              endif
              if(.not.allocated(iszerolist)) then
                  allocate(iszerolist(ARRAY_START_LENGTH))
              endif
          endif
          allocate(outepos(DET_MAX_LENGTH))
          !apply l+ on detup
          tpdet = detup
          do while(tpdet.ne.0)
              i = trailz(tpdet) ! starting from 0
              outdet = detup
              outepos(:) = eposup(:)
              outcoef = coef
              outis0 = iszero.or.equals0(coef)
              call lpls(i, outdet, outepos, outcoef, outis0)
              !Original det is not added to the list here when outis0 = true
              write(*, '("l+ on e at up pos:", 1I4)') i
              write(*, *) outis0
              if(.not.outis0) then
                  num = num + 1
                  basislist(num, 1:2) = (/outdet, detdn/)
                  coeflist(num) = outcoef
                  eposlist(num, 1, :) = outepos(:)
                  eposlist(num, 2, :) = eposdn(:)
                  iszerolist(num) = outis0
                  write(*, '("Output detup, detdn: ", 2B16)') basislist(num, 1:2)
                  write(*, '(15I3)') eposlist(num, 1, 1:15)
                  write(*, '(15I3)') eposlist(num, 2, 1:15)
                  write(*, *) coeflist(num), iszerolist(num), num
              endif

              tpdet = ibclr(tpdet, i) 
          enddo
           !apply l+ on detdn
          tpdet = detdn
          do while(tpdet.ne.0)
              i = trailz(tpdet) ! starting from 0
              outdet = detdn
              outepos(:) = eposdn(:)
              outcoef = coef
              outis0 = iszero.or.equals0(outcoef)
              call lpls(i, outdet, outepos, outcoef, outis0)
              !Original det is not added to the list here when outis0 = true
              write(*, '("l+ on e at dn pos:", 1I4)') i
              write(*, *) outis0
              if(.not.outis0) then
                  num = num + 1
                  basislist(num, 1:2) = (/detup, outdet/)
                  coeflist(num) = outcoef
                  eposlist(num, 1, :) = eposup(:)
                  eposlist(num, 2, :) = outepos(:)
                  iszerolist(num) = outis0
                  write(*, '("Output detup, detdn: ", 2B16)') basislist(num, 1:2)
                  write(*, '(15I3)') eposlist(num, 1, 1:15)
                  write(*, '(15I3)') eposlist(num, 2, 1:15)
                  write(*, *) coeflist(num), iszerolist(num), num
              endif

              tpdet = ibclr(tpdet, i) 
          enddo
         deallocate(outepos)

     end subroutine Lplus_single
          

end module projection
      

