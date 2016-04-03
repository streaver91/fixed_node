module types
  !---------------------------------------------------------------------------
  ! Description : Define various kinds of integers of various bit sizes
  ! Created     : Frank Petruzielo and Cyrus Umrigar
  !---------------------------------------------------------------------------

  implicit none ; save ; private
  public :: i1b, i4b, i8b, i16b, ik, rk , ik_vec, num_words, bits_per_word

  integer, parameter :: i1b =  SELECTED_INT_KIND(2)  ! 1 byte = 8 bits
  integer, parameter :: i2b =  SELECTED_INT_KIND(4)  ! 2 bytes = 16 bits
  integer, parameter :: i4b =  SELECTED_INT_KIND(9)  ! 4 bytes = 32 bits (default)
  integer, parameter :: i8b =  SELECTED_INT_KIND(18) ! 8 bytes = 64 bits
  integer, parameter :: i16b = SELECTED_INT_KIND(38) ! 16 bytes = 128 bits, approx 10^38 (unfortunately ifort does not have this, but gfortran does)

  integer, parameter :: r4b =  SELECTED_REAL_KIND(6)  ! 4 bytes = 32 bits (default)
  integer, parameter :: r8b =  SELECTED_REAL_KIND(15) ! 8 bytes = 64 bits (double precision what we usually use)
  integer, parameter :: r16b = SELECTED_REAL_KIND(30) ! 16 bytes = 128 bits (unfortunately gfortran does not have this)

! When we change precision, the only lines that need changing are the following two:
  integer, parameter :: ik = i16b
  integer, parameter :: rk = r8b

  type ik_vec
    integer(ik),dimension(2) :: v
  end type ik_vec

#ifdef NUM_ORBITALS_GT_127
  integer, parameter :: num_words=2
#else
  integer,parameter  :: num_words=1
#endif
  integer,parameter  :: bits_per_word=127

! Examples of using it:
! integer(i1b)  :: i1
! integer(i2b)  :: i2
! integer(i4b)  :: i4
! integer(i8b)  :: i8
! integer(i16b) :: i16
! real(r4b)     :: r4
! real(r8b)     :: r8
!!real(r16b)    :: r16

! write(6,*) huge(i1), huge(i2), huge(i4), huge(i8), huge(i16), huge(r4), huge(r8)
!!gives:
!!127  32767  2147483647  9223372036854775807 170141183460469231731687303715884105727  3.40282347E+38  1.79769313486231571E+308
! write(6,*) tiny(r4), tiny(r8)
! 1.17549435E-38  2.22507385850720138E-308

! Examples of using it:
! #ifdef NUM_ORBITALS_GT_127
! type(ik_vec)          :: hubbard_config
! #else
! integer(ik)           :: hubbard_config
! #endif

! i = 2_ik**20 * 3_ik**10
! real(rk), allocatable :: hamiltonian(:,:)
! complex(rk)           :: phase


! Outdated:
! integer, parameter :: rk = kind(1.0_rk)
! integer, parameter :: rkc = kind((1.0_rk,1.0_rk))

end module types
