module projection
      !Generate csf by projection
      use prep, only: i16b, ARRAY_START_LENGTH, DET_MAX_LENGTH
      use, intrinsic :: iso_fortran_env, only: rk => real64
      use detgen, only: lpls, lms, eposinit, spls, sms

      implicit none
      contains

!TODO collect all  determinants and coefs

      subroutine getallsigns(basislist, coeflist, eposlist, iszerolist, num)
          integer(i16b), intent(in) :: basislist(:, :)
          real(rk) :: coeflist(:)
          integer, intent(in) :: eposlist(:, :, :), num
          logical, intent(in) :: iszerolist(:)
          integer :: i
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
            write(*, *) countswps
            end function countswps

          





      !TODO collect determinants


      subroutine Sminus_multiple(inbasis, incoefs, inepos, iniszeros, innum, &
              & basislist, coeflist, eposlist, iszerolist, num)
          !append to output list
          integer(i16b), allocatable, intent(in) :: inbasis(:, :)
          real(rk), allocatable, intent(in) :: incoefs(:)
          integer, allocatable, intent(in) :: inepos(:, :, :)
          logical, allocatable, intent(in) :: iniszeros(:)
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
            write(*, '(">>>>>>>>>>>>> Applying S- to multiple dets, now at No.", 1I3)') i 
            call Sminus_single(inbasis(i,1), inbasis(i, 2), inepos(i, 1, :), &
                &inepos(i, 2, :), incoefs(i), iniszeros(i), &
                & basislist, coeflist, eposlist, iszerolist, num)
          enddo
      end subroutine Sminus_multiple
      
      subroutine Splus_multiple(inbasis, incoefs, inepos, iniszeros, innum, &
              & basislist, coeflist, eposlist, iszerolist, num)
          !append to output list
          integer(i16b), allocatable, intent(in) :: inbasis(:, :)
          real(rk), allocatable, intent(in) :: incoefs(:)
          integer, allocatable, intent(in) :: inepos(:, :, :)
          logical, allocatable, intent(in) :: iniszeros(:)
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
            write(*, '(">>>>>>>>>>>>> Applying S+ to multiple dets, now at No.", 1I3)') i 
            call Splus_single(inbasis(i,1), inbasis(i, 2), inepos(i, 1, :), &
                &inepos(i, 2, :), incoefs(i), iniszeros(i), &
                & basislist, coeflist, eposlist, iszerolist, num)
          enddo
      end subroutine Splus_multiple



      subroutine Lminus_multiple(inbasis, incoefs, inepos, iniszeros, innum, &
              & basislist, coeflist, eposlist, iszerolist, num)
          !append to output list
          integer(i16b), allocatable, intent(in) :: inbasis(:, :)
          real(rk), allocatable, intent(in) :: incoefs(:)
          integer, allocatable, intent(in) :: inepos(:, :, :)
          logical, allocatable, intent(in) :: iniszeros(:)
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
            write(*, '(">>>>>>>>>>>>> Applying L- to multiple dets, now at No.", 1I3)') i 
            call Lminus_single(inbasis(i,1), inbasis(i, 2), inepos(i, 1, :), &
                &inepos(i, 2, :), incoefs(i), iniszeros(i), &
                & basislist, coeflist, eposlist, iszerolist, num)
          enddo
      end subroutine Lminus_multiple

      subroutine Lplus_multiple(inbasis, incoefs, inepos, iniszeros, innum, &
              & basislist, coeflist, eposlist, iszerolist, num)
          !append to outaut list
          integer(i16b), allocatable, intent(in) :: inbasis(:, :)
          real(rk), allocatable, intent(in) :: incoefs(:)
          integer, allocatable, intent(in) :: inepos(:, :, :)
          logical, allocatable, intent(in) :: iniszeros(:)
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
            write(*, '(">>>>>>>>>>>>> Applying L+ to multiple dets, now at No.", 1I3)') i 
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
          write(*, '("-----------Apply S- on detup, detdn: ", 2B16)') detup, detdn
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
              outis0 = iszero
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
          write(*, '("-----------Apply S+ on detup, detdn: ", 2B16)') detup, detdn
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
              outis0 = iszero
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
          write(*, '("-----------Apply L- on detup, detdn: ", 2B16)') detup, detdn
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
              outis0 = iszero
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
              outis0 = iszero
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
          write(*, '("-----------Apply L+ on detup, detdn: ", 2B16)') detup, detdn
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
              outis0 = iszero
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
              outis0 = iszero
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
            

          
      !TODO Lz
      !TODO

end module projection
      

