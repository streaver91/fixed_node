module tools
  use types, only: ik, ik_vec, i16b, rk
  implicit none
  save
  private
  public :: merge_sort_real

contains

recursive subroutine merge_sort_real(key_real, iorder, nelts, temp_real, temp_i_2)
! ===========================================================================================================
! Description   : Merge sort in ascending order.
!               : It sorts nwalk items in key which is real.  Index is then used outside this routine to sort auxilliary items.
!               : temp_real is an auxilliary array of length (nwalk+1)/2 of same type as key.
!               : temp_i_2 is an auxilliary array of length (nwalk+1)/2 of same type as iorder.
!               : In the present version key is real(rk) and iorder is integer.
!               : temp_i16 and temp_i_2 are passed in to avoid overhead of automatic arrays or allocations.
! Adapted by    : Cyrus Umrigar
! -----------------------------------------------------------------------------------------------------------
  use types, only: ik
  implicit none
  real(rk), intent(inout) :: key_real(nelts)
  integer, intent(inout) :: iorder(nelts)
  integer, intent(in) :: nelts
  real(rk), dimension((nelts+1)/2), intent (out) :: temp_real
  integer,  dimension((nelts+1)/2), intent (out) :: temp_i_2     ! temporary array actually neither in nor out
! Local variables
  real(rk)  t0
  integer :: na,nb,i

  if (nelts < 2) return
  if (nelts == 2) then
    if (key_real(1) > key_real(2)) then ! If key_real's in wrong order, sort on key_real
      t0 = key_real(1)
      key_real(1) = key_real(2)
      key_real(2) = t0
      i = iorder(1)
      iorder(1) = iorder(2)
      iorder(2) = i
    endif
    return
  endif

  na=(nelts+1)/2
  nb=nelts-na
  call merge_sort_real(key_real, iorder, na, temp_real, temp_i_2)
  call merge_sort_real(key_real(na+1), iorder(na+1), nb, temp_real, temp_i_2)

  if (key_real(na) > key_real(na+1)) then
    temp_real(1:na)=key_real(1:na)
    temp_i_2(1:na)=iorder(1:na)
    call merge_real(temp_real, temp_i_2, na, key_real(na+1), iorder(na+1), nb, key_real, iorder, nelts)
  endif

  return
end subroutine merge_sort_real
!-----------------------------------------------------------------------------------------------

subroutine merge_real(a,a2,na,b,b2,nb,c,c2,nc)
! ===============================================================================================
! Description   : Called by merge_sort2_up_dn
! Adapted by    : Cyrus Umrigar
! Comments      : Pulled out from semistoch_mod.f90 and chemistry.f90 and put in tools.f90 for general purpose use
! -----------------------------------------------------------------------------------------------
  use types, only: ik
  implicit none
  integer, intent(in) :: na,nb,nc                                      ! na+nb = nc
  real(rk), intent(inout) :: a(na), c(nc) ! b_up overlays c_up(na+1:nc), and similarly for b_dn, c_dn
  integer    , intent(inout) :: a2(na), c2(nc)                         ! b2   overlays c2(na+1:nc)
  real(rk), intent(in)    :: b(nb)
  integer    , intent(in)    :: b2(nb)
  integer :: i,j,k

  i = 1; j = 1; k = 1;
  do while(i <= na .and. j <= nb)
    if (a(i) < b(j)) then
      c(k) = a(i)
      c2(k) = a2(i)
      i = i+1
    else
      c(k) = b(j)
      c2(k) = b2(j)
      j = j+1
    endif
    k = k + 1
  enddo

  do while (i <= na)
    c(k) = a(i)
    c2(k) = a2(i)
    i = i + 1
    k = k + 1
  enddo

  return
end subroutine merge_real
!----------------------------------------------------------------------------------------------------------------------
end module tools
