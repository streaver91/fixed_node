module more_tools
  use types, only: ik, ik_vec, rk
  use tools, only: merge_sort_real
  implicit none
  save
  private
  public :: real_symmetric_diagonalize, real_general_diagonalize_ow_ham

contains

!=====================================================================================================================
 subroutine real_symmetric_diagonalize(n,M,ham,eigenvectors,eigenvalues)
 ! Created by : Hitesh Changlani, March 12,2012
!=====================================================================================================================

 implicit none

 integer,intent(in)   :: n, M
 real(rk),intent(in)  :: ham(M,n)
 real(rk),intent(out) :: eigenvectors(M,n),eigenvalues(n)

 integer              :: len_work,info
 real(rk),allocatable :: work(:)

 len_work = 3 * n -1
 allocate(work(len_work))

 eigenvectors=ham

! On input eigenvectors has the matrix, and on output the eigenvectors
 call dsyev('V', 'U', n, eigenvectors, M, eigenvalues, work, len_work, info)
 deallocate(work)

 if (info /= 0) then
        write(6,*) "info = ", info
        stop 'Diagonalization Failed!'
 endif

 end subroutine real_symmetric_diagonalize

!=====================================================================================================================
 subroutine real_general_diagonalize_ow_ham(n,M,ham,re_eigenvalues,im_eigenvalues)
 ! Created by : Hitesh Changlani, October 17,2012
!=====================================================================================================================

 implicit none

 integer,intent(in)     :: n, M
 real(rk),intent(inout) :: ham(M,n)
 real(rk),intent(out)   :: re_eigenvalues(n)
 real(rk)               :: im_eigenvalues(n)
 real(rk)               :: vl(M,n),vr(M,n)
 integer                :: i,len_work,info
 real(rk)               :: work(8*n)
 integer                :: iorder(n)
 real(rk), dimension((n+1)/2) :: temp_real
 integer,  dimension((n+1)/2) :: temp_i_2     ! temporary array actually neither in nor out

 len_work = 8* n

 call dgeev('N', 'V', n, ham, M, re_eigenvalues, im_eigenvalues,vl,M,vr,M, work, len_work, info)
 ! Only the right eigenvectors are computed
 ! However unlike the symmetric case, the eigenvalues are not sorted

 if (info /= 0) then
        write(6,*) "info = ", info
        stop 'Diagonalization Failed!'
 endif

 do i=1,n
    iorder(i)=i
 enddo

 !write (6,*) "Before sort"
 !call flush(6) 
 !do i=1,n
 !   write(6,*) re_eigenvalues(i),im_eigenvalues(i)
 !   call flush(6)
 !enddo

 call merge_sort_real(re_eigenvalues, iorder, n, temp_real, temp_i_2)

 im_eigenvalues(1:n)=im_eigenvalues(iorder(1:n))
 do i=1,n
    ham(:,i)=vr(:,iorder(i))
 enddo
 
 !write (6,*) "After sort"
 !call flush(6) 
 !do i=1,n
 !   write(6,*) re_eigenvalues(i),im_eigenvalues(i)
 !   call flush(6)
 !enddo

 end subroutine real_general_diagonalize_ow_ham
!=====================================================================================================================
end module more_tools
