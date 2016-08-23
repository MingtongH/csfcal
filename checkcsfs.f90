module checkcsfs
       use prep, only: i16b, ARRAY_START_LENGTH, ARRAY_SHORT_LENGTH, DET_MAX_LENGTH
       use, intrinsic :: iso_fortran_env, only: rk => real64

 
       implicit none
    
       contains

       subroutine checkcsfLsq(inbasislist, incoeflist, L, nbasis, ncsf)
           integer(i16b), intent(in) :: inbasislist(:, :)
           real(rk), intent(in) :: incoeflist(:)
           integer, intent(in) :: nbasis, ncsf, L
           integer :: i, j, k
 !          do i = 1, ncsf
 !              do j = 1, nbasis
 !                 call initlists(inbasislist(j, 1:2), baisslist1, coeflist1, eposes1,       iszeros1,   
       endsubroutine checkcsfLsq
 
       subroutine sortBasisCoefTable(inbasislist, incoeftable, nbasis, ncsf)
           integer(i16b) :: inbasislist(:, :)
           integer(i16b), allocatable:: tpbasis(:, :)
           real(rk) :: incoeftable(:, :)
           real(rk), allocatable :: tpcoefs(:, :)
           integer, intent(in) :: nbasis, ncsf
           integer :: i, j, nordered, k
           integer, allocatable :: order(:) !indices of dets, from small to large
           allocate(tpbasis(nbasis, 2))
           allocate(tpcoefs(nbasis, ncsf))
           allocate(order(nbasis))
           do i = 1, nbasis
               tpbasis(i, 1:2) = inbasislist(i, 1:2)
               tpcoefs(i, 1:ncsf) = incoeftable(i, 1:ncsf)
           enddo
           order(1) = 1 !
           nordered = 1 ! number of ordered dets
           do i = 2, nbasis
                do j = 1, nordered

                    if(tpbasis(i, 1).lt.tpbasis(order(j), 1)) then  !1<
                        do k = nordered, j, -1
                            order(k + 1) = order(k)
                        enddo
                        order(j) = i
                        nordered = nordered + 1
                        exit

                    else if(tpbasis(i, 1).eq.tpbasis(order(j), 1)) then !1=
                        if(tpbasis(i, 2).lt.tpbasis(order(j), 2)) then  !2<
                            do k = nordered, j, -1
                                order(k + 1) = order(k)
                            enddo
                            order(j) = i
                            nordered = nordered + 1
                            exit
                        elseif(tpbasis(i, 2).gt.tpbasis(order(j), 2)) then !2>
                            !i has same updet as order(j) but bigger dndet
                            !should be placed behind order(j)
                            if(nordered.eq.j) then   !3
                                order(j+1) = i
                                nordered = nordered + 1
                                exit
                            else
                                do k = nordered, j+1, -1
                                    order(k + 1) = order(k)
                                enddo
                                order(j + 1) = i
                                nordered = nordered + 1
                                exit
                            endif        !3
                        endif        !2

                    else!1 basis i, 1 > basis order(j)
                        if(nordered.eq.j) then
                            order(j+1) = i
                            nordered = nordered + 1
                            exit
                        endif
                    endif !1
                    
                enddo !j
            enddo !i

            write(*, *) '----------------------- Sorting csf table by dets----------------'

            do i = 1, nbasis
                write(*, *) order(i)
            enddo

            do i = 1, nbasis
                inbasislist(i, 1:2) = tpbasis(order(i), 1:2)
                incoeftable(i, 1:ncsf) = tpcoefs(order(i), 1:ncsf)
            enddo
            deallocate(tpbasis)
            deallocate(tpcoefs)
            deallocate(order)

       endsubroutine sortBasisCoefTable
end module checkcsfs
