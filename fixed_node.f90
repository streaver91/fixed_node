module fn

! This program tests/demonstrates a few different things:

! 1) 
! Calculates the fraction of the energy missing in VMC that is recovered by 3 different fixed-node (FN) approximations.
! The 2 new ones have the virtue that the projector elements never become negative when tau is increased.
! FN1 is the standard fixed-node that moves off-diagonal negative elements to the diagonal.
! FN3 does a heat-bath using absolute values of the column of the importance-sampled projector setting node-crossing elements to 0,
!     and an accept-reject, and, reweights by the local projector.
!     It is exact when Psi_G is exact because VMC is exact and the local projector is 1 when Psi_G=Psi_0.
! FN4 does a heat-bath using absolute values of the column of the importance-sampled projector setting node-crossing elements to 0,
!     and, reweights by the local projector.
!     It is exact when there is no sign problem since it is then sampling from the true projector.
! FN5 Removes the off-diagonal negative elements and rescales the off-diagonal positive elements so that the column sum is unchanged.
! FN1 has no time-step error, but has a maximum time-step (which can be very small) beyond which there is a sign-problem because the diagonal element becomes negative.
! FN3 and FN4 never have a sign problem, but have a time-step error.  In the limit of small time steps they reduce to FN1.  The value of the time-step error depends on E_T
! FN5 has no time-step error, but it can give energies that are lower than the true energy as well energies that are
! considerably higher than FN1, and so it is of no interest.

! 2)
! Calculates the parameter derivative of the FN1 energy, both exactly and using 2 different approximations.
! One of the approximations is to replace the 2nd term by the integral over the chosen path of the reweighting term.
! This is the derivative version of the approximation of FN-DMC forces Claudia and I introduced in 2000.
! The other approximation is to replace the 2nd term by the 3rd term.
! Find empirically that
! a) When psi_G is close to exact, the 2nd and 3rd terms are equal to linear order in psi_G-psi_0.
!    Since the 2nd term is unknown, one can replace the sum of the 2nd and 3rd terms by 2 times the 3rd term.
!    So, the 1st term plus 2 times the 3rd term give the approximate energy derivative.
!    This has the advantage that the 2nd order divergence in the 1st term is exactly cancelled as pointed out by Sorella.
! b) When there is no sign problem the derivative of the energy of course is zero.
!    In that case, one can show that the 2nd and 3rd terms are equal to each other and equal to -1/2 times the first term.
! c) When there is no sign problem, but psi_G has the wrong sign on some components, then the entire space
!    breaks up into 2 disjoint spaces as expected, and, there are 2 degenerate dominant eigenvalues, one in each space.

! FN = fixed-node
! FD = finite difference

 
use types, only: rk
use more_tools, only: real_symmetric_diagonalize, real_general_diagonalize_ow_ham
implicit none

integer, parameter :: MDIM = 100
real(rk), parameter :: eps = 2.e-15, eps2 = 1.e-7, eps_FD = 1.e-7, eps_FD_2 = 1.e-4
integer, dimension(4) :: irand_seed
integer ndim, i, j, iter, itimes, info
real(rk) :: H(MDIM,MDIM), H_imp(MDIM,MDIM), H_imp_i_off(MDIM), &
P(MDIM,MDIM), P_imp(MDIM,MDIM), P_FN1(MDIM,MDIM), P_FN3(MDIM,MDIM), P_FN4(MDIM,MDIM), P_FN5(MDIM,MDIM), P_propos(MDIM,MDIM), P_markov(MDIM,MDIM), &
psio(MDIM), psin(MDIM), psi0(MDIM), psiG(MDIM), psi0_psiG(MDIM), psiG_psiFN(MDIM), &
psiFN1(MDIM), psiFN3(MDIM), psiFN4(MDIM), psiFN5(MDIM), &
P_L(MDIM), P_abs_L(MDIM), &
tau_in, tau, diagonal_dump, tau_av_expon_proj(MDIM), tau_av_expon_proj_new(MDIM), fraction_positive, off_diag_mag, psiG_noise, &
E_T, E_G, E_0, E_FN1, E_FN3, E_FN4, E_FN5, norm, &
eigenval_p, eigenval_h, eigenval_sav, &
overlap_G, overlap_FN1, overlap_FN3, overlap_FN4, overlap_FN5, overlap_G_FN1, &
sign_flip, rejected, sum_P_eig, sum_abs_P_eig, sign_condition_number, &
eigenvalues(MDIM), eigenvalues_imag(MDIM), eigenvectors(MDIM,MDIM), eigenvectors_left(1,1), eigenvectors_right(MDIM,MDIM), work(4*MDIM), &
sum_pos, sum_neg, ratio
integer iparm
real(rk) :: E_G_increm, E_FN1_increm, psiG_increm(MDIM), psiG_psiFN1_increm(MDIM), psiFN1_increm(MDIM), psiFN1_der(MDIM), E_L(MDIM),  E_L_increm(MDIM), E_L_der_increm(MDIM), P_L_increm(MDIM), PP_FN1(MDIM,MDIM), q_FN1(MDIM,MDIM), W(MDIM), W_der(MDIM), W_der_by_W_sav, psiFN1_by_psiG(MDIM), psiFN1_by_psiG_increm(MDIM), psiFN1_by_psiG_der(MDIM), &
P_FN1_P_L_der_by_P_L(MDIM,MDIM), &
tmp_vec1(MDIM), tmp_vec2(MDIM), tmp_vec3(MDIM), tmp_vec1_i(MDIM), tmp_vec1_j(MDIM),tmp_vec4(MDIM),&
H_minus_E_psiFN1(MDIM), H_minus_E_psiFN1_der(MDIM), H_minus_E_psiFN1_by_psiG(MDIM), H_minus_E_psiG_der(MDIM), H_minus_E_psi0_min_psiFN1(MDIM), &
test5a,test5b,test5c,test5d,test5e,test5f,test5g,test7
real(rk) :: a_2,psiG_2der(MDIM,MDIM,MDIM),gradient_FN1_exact(MDIM),hessian_FN1(MDIM,MDIM),&
     &psiG_der(MDIM,MDIM),hessian_VMC_FD(MDIM,MDIM),gradient_FN1_approx(MDIM),&
     &hessian_FN_FD(MDIM,MDIM),gradient_FN_FD(MDIM),gradient_VMC(MDIM),E_L_der(MDIM,MDIM),&
     &hessian_vmc_in_dmc(MDIM,MDIM),gradient_VMC_FD(MDIM),gradient_vmc_in_dmc(MDIM),&
     &hessian_vmc(MDIM,MDIM),E_L_2der(MDIM,MDIM,MDIM),hessian_vmc_err(MDIM,MDIM),&
     &P_L_der(MDIM,MDIM),P_L_der_by_P_L(MDIM,MDIM),P_L_2der(MDIM,MDIM,MDIM),&
     &P_L_2der_by_P_L(MDIM,MDIM,MDIM),W_der_by_W(MDIM),W_2der_by_W(MDIM,MDIM),gradient_FU(MDIM),&
     &W_der_by_W_full(MDIM,MDIM),W_2der_by_W_full(MDIM,MDIM,MDIM),&
     &PP_FN1_minus_I(MDIM,MDIM)

real(rk) :: rannyu
character(len=80) :: fmt, hamiltonian
logical psiG_correct_sign, plot_deriv

integer npsiG, ipsiG,k
real(rk) :: psiG_sav(MDIM)

external rannyu

contains

subroutine read_input

open(1,file='summary1')
write(1,'(''Calculate fraction of energy missing in VMC that is recovered by FN-PMC and 1 minus the overlap of Psi_FN*Psi_G with Psi_0*Psi_G'')')
write(1,'(''--------------------------------------------------------------------------------------------------------------------------'')')
write(1,'(''                                                                                   1 -      1 -      1 -                  '')')
write(1,'('' E_G        E_FN1     E_FN3     E_FN4     E_0       rec_FN1   rec_FN3   rec_FN4   ovlp_FN1 ovlp_FN3 ovlp_FN4 sign_cond_num'')')
write(1,'(''--------------------------------------------------------------------------------------------------------------------------'')')

open(2,file='summary2')
write(2,'(''psiG_noise, 1-overlap_FN1 ipsig grad_FN_exact(1)  grad_FU(1)  grad_FN_approx(1) grad_VMC(1) grad_VMC_in_DMC(1) grad_FN_exact(iparm)  grad_FU(iparm)  grad_FN_approx(iparm) grad_VMC(iparm) grad_VMC_in_DMC(iparm) grad_FN_exact(iparm+1)  grad_FU(iparm+1)  grad_FN_approx(iparm+1) grad_vmc(iparm+1) grad_vmc_in_dmc(iparm+1)    '')')

open(3,file='summary3')
write(3,'(''psiG_noise, 1-overlap_FN1 ipsig h_fd(1,1) h_approx(1,1) h_VMC(1,1) h_VMC_in_DMC(1,1) h_fd(iparm,iparm)  h_approx(iparm,iparm)  h_VMC(iparm,iparm) h_VMC_in_DMC(iparm,iparm) h_fd(iparm,iparm+1)  h_approx(iparm,iparm+1)  h_vmc(iparm,iparm+1) h_vmc_in_dmc(iparm,iparm+1)    '')')

! We are presently reading E_T, but not using it.  We set it to E_0 in code for P_imp, and to E_FN1 for the other projectors.
write(6,'(''Input hamiltonian'')')
read(5,*) hamiltonian
write(6,'(''hamiltonian='',a)') hamiltonian
write(6,'(''Input ndim, E_T, tau_in'')')
read(5,*) ndim, E_T, tau_in
write(6,'(''ndim, E_T, tau_in='',i5,f10.6,9f9.6)') ndim, E_T, tau_in
tau=tau_in
write(6,'(''Input fraction_positive, off_diag_mag, psiG_noise'')')
read(5,*) fraction_positive, off_diag_mag, psiG_noise
write(6,'(''fraction_positive, off_diag_mag,psiG_noise='',9f9.6)') fraction_positive, off_diag_mag, psiG_noise
write(6,'(''psiG_correct_sign'')')
read(5,*) psiG_correct_sign
write(6,'(''psiG_correct_sign='',l2)') psiG_correct_sign
read(5,*) plot_deriv
write(6,'(''plot_deriv='',l2)') plot_deriv

if(fraction_positive.lt..5_rk .or. fraction_positive.gt.1._rk) stop '(fraction_positive shoud be between .5 and 1)'
if(ndim.gt.MDIM) stop 'ndim > MDIM'

write(6,'(''Input irand_seed'')')
read(5,'(4i4)') irand_seed
write(6,'(''irand_seed='',4i4)') irand_seed
call setrn(irand_seed)

do itimes=1,1
write(6,'(/,''Data set'',i5)') itimes

! Define hamiltonian, H.  The projector, P, will be defined once we calculate the ground state energy, E_0.

do j=1,ndim
  do i=j,ndim
    if(i.eq.j) then
      if(hamiltonian.eq.'random') then
        H(j,i)=real(-ndim+2*i-1,rk)/(ndim-1)
      elseif(hamiltonian.eq.'hubbard_r_like') then
        if(rannyu(0).lt.0.7_rk) then
          H(j,i)=0
        else
          H(j,i)=4
        endif
      elseif(hamiltonian.eq.'hubbard_k_like') then
        H(j,i)=(i-1)**2
      endif
    else
      if(hamiltonian.eq.'random') then
        H(j,i)=off_diag_mag*(rannyu(0)-fraction_positive) ! If fraction_positive=1, there is no sign problem
      elseif(hamiltonian.eq.'hubbard_r_like') then
        if(rannyu(0).lt.fraction_positive) then
          H(j,i)=-1
        else
          H(j,i)=1
        endif
      elseif(hamiltonian.eq.'hubbard_k_like') then
        if(rannyu(0).lt.fraction_positive) then
          H(j,i)=-4
        else
          H(j,i)=4
        endif
      endif
    endif
    H(i,j)=H(j,i)
  enddo
enddo

write (fmt,'(a,i2,a)') '(/,''H='',/,(', ndim, 'f9.4))'
write(6,fmt) ((H(j,i),i=1,ndim),j=1,ndim)

call real_symmetric_diagonalize(ndim,MDIM,H,eigenvectors,eigenvalues)
eigenval_h=eigenvalues(1)
E_0=eigenval_h
E_T=E_0
eigenval_p=1+tau*(E_0-eigenval_h)
psi0(1:ndim)=eigenvectors(1:ndim,1)
if(psi0(1).lt.0) psi0(1:ndim)=-psi0(1:ndim)
norm=sqrt(sum(psi0(1:ndim)**2))
psi0=psi0/norm
write(6,'(/,''eigenval_p,eigenval_h='',9f10.6)') eigenval_p, eigenval_h
write (6,'(''psi for P    ='',100f9.4)') (psi0(i),i=1,ndim)
write(6,'(/,''eigenvalues of P='',100f10.6)') (1+tau*(E_0-eigenvalues(i)),i=ndim,1,-1)
write(6,'(''eigenvalues of H='',100f10.6)') (eigenvalues(i),i=ndim,1,-1)

write(6,'(''tau, tau_max='',9f9.6)') tau, 1/(eigenvalues(ndim)-eigenvalues(1))
!tau=1/(eigenvalues(ndim)-eigenvalues(1))

! Check that tau is smaller than 1/2 the maximum allowed by the spectral range.
! We put in the factor of 1/2 so that none of the diagonal elements become negative, though this is not needed for deterministic projection.
if(tau.gt.1/(eigenvalues(ndim)-eigenvalues(1))) stop 'tau is too large for this spectral range'

! Define projector P
do j=1,ndim
  do i=j,ndim
    if(i.eq.j) then
      P(j,i)=1+tau*(E_0-H(j,i))
    else
      P(j,i)=-tau*H(j,i)
    endif
    P(i,j)=P(j,i)
  enddo
enddo

write (fmt,'(a,i2,a)') '(''P='',/,(', ndim, 'f9.4))'
write(6,fmt) ((P(j,i),i=1,ndim),j=1,ndim)

! Calculate the sign-condition number
! Later, I construct the sign-condition number for the importance-sampled projector and its eigenvector
sum_P_eig=0 ; sum_abs_P_eig=0
do j=1,ndim
  do i=1,ndim
    sum_abs_P_eig=sum_abs_P_eig+abs(P(j,i)*psi0(i))
  enddo
  sum_P_eig=sum_P_eig+abs(psi0(j))
enddo
sign_condition_number=(sum_P_eig/sum_abs_P_eig)**(1/tau)
write(6,'(/,''sum_P_eig,sum_abs_P_eig,sign_condition_number='',9f10.6)') sum_P_eig,sum_abs_P_eig,sign_condition_number

! To see sign-flip terms more clearly, reconstruct P and H with sign of basis fns. changed to make eigenstate have nonnegative elements
do j=1,ndim
  if (psi0(j).lt.0) then
    psi0(j)=-psi0(j)
    do i=1,ndim
      if(i.ne.j) then
        H(j,i)=-H(j,i)
        H(i,j)=-H(i,j)
        P(j,i)=-P(j,i)
        P(i,j)=-P(i,j)
      endif
    enddo
  endif
enddo
write (fmt,'(a,i2,a)') '(/,''P after making eigenvec +ve ='',/,(', ndim, 'f9.4))'
write(6,fmt) ((P(j,i),i=1,ndim),j=1,ndim)
write(6,'(/,''1+(E_0*1-H) (i.e. P with tau=1) after making eigenvec +ve ='')')
do j=1,ndim
  do i=1,ndim
    if(i.eq.j) then
      write(6,'(f9.4)',advance='no') 1+E_0-H(j,i)
    else
      write(6,'(f9.4)',advance='no') H(j,i)
    endif
  enddo
  write(6,*)
enddo

! Define approximate guiding fn
do i=1,ndim
  if(psiG_correct_sign) then
    psiG(i)=psi0(i)+psiG_noise*sign(abs(rannyu(0)-0.5_rk),psi0(i)) ! Add psiG_noise maintaining correct signs
  else
    psiG(i)=psi0(i)+psiG_noise*(rannyu(0)-0.5_rk) ! Add noise with either sign
  endif
enddo

iparm=2

!Assume that the WF is:
!
! a_1*e_1 + a_3*exp(a_2)*e_2 + a_3 * e_3 + a_4 * e_4 +...
!
! where a_i are the parameters and e_1 is a unit vector, such that
! the derivative of the linear parameters is just 1 on that site

a_2 = log(psiG(iparm)/psiG(iparm+1))


! This is for plotting the parameter derivative of E_VMC and E_FN versus the parameter value
if(plot_deriv) then
  npsiG=11
  psiG_sav=psiG
!  psiG_sav=psi0
else
  npsiG=1
endif

do ipsiG=1,npsiG

if(plot_deriv) then
  psiG=psiG_sav
  psiG(1) = psi0(1)+.01*(ipsiG-6)
  psiG(iparm+1) = psi0(iparm+1)+.01*(ipsiG-6)
  a_2 = log(psi0(iparm)/psi0(iparm+1))
  psiG(iparm) = psiG(iparm+1)*exp(a_2+.00000001*(ipsiG-6))
endif

norm=sqrt(sum(psiG(1:ndim)**2))
psiG=psiG/norm
overlap_G=overlap(ndim,psi0,psiG)

write (6,'(''psi0         ='',100f9.4)') (psi0(i),i=1,ndim)
write (6,'(''psiG         ='',100f9.4)') (psiG(i),i=1,ndim)
write (6,'(''1-overlap_G  ='',100f9.4)') 1-overlap_G

! Calculate the variational energy for psiG
E_G=0
do j=1,ndim
  do i=1,ndim
    E_G=E_G+psiG(j)*H(j,i)*psiG(i)
  enddo
enddo
write(6,'(/,''E_G='',f10.6)') E_G

! Define importance-sampled projector
do j=1,ndim
  do i=1,ndim
    H_imp(j,i)=psiG(j)*H(j,i)/psiG(i)
    P_imp(j,i)=psiG(j)*P(j,i)/psiG(i)
  enddo
enddo
write (fmt,'(a,i2,a)') '(/,''H_imp='',/,(', ndim, 'f9.4))'
write(6,fmt) ((H_imp(j,i),i=1,ndim),j=1,ndim)
write (fmt,'(a,i2,a)') '(/,''P_imp='',/,(', ndim, 'f9.4))'
write(6,fmt) ((P_imp(j,i),i=1,ndim),j=1,ndim)

! Calculate the average time for an exponential projector
do i=1,ndim
  H_imp_i_off(i)=0
  do j=1,ndim
    if(j.ne.i) H_imp_i_off(i)=H_imp_i_off(i)-abs(H_imp(j,i))
  enddo
  tau_av_expon_proj(i)=-1/H_imp_i_off(i)
  tau_av_expon_proj_new(i)=-1/(E_T-H_imp(i,i))
enddo
write(6,'(''tau_av_expon_proj    ='',100f9.6)') tau_av_expon_proj(1:ndim)
write(6,'(''tau_av_expon_proj_new='',100f9.6)') tau_av_expon_proj_new(1:ndim)
write(6,'(''factor=               '',100f9.6)') (exp((E_T-H_imp(i,i)-H_imp_i_off(i))/abs(H_imp_i_off(i))), i=1,ndim)
!write(6,'(''factor2='',100(3f9.4,2x))') (E_T,H_imp(i,i),H_imp_i_off(i), i=1,ndim)
write(6,'(''factor2=              '',100f9.4)') (H_imp_i_off(i)/(E_T-H_imp(i,i)), i=1,ndim)

! Make sure that P_imp has the same eigenvalue as P

! I was using power method to find dominant state of P but now use LAPACK instead
!eigenval_sav=9.e99_rk
!do
!  call mat_vec_mul(ndim, MDIM, P_imp, psio, psin)
!  eigenval_p=sqrt(sum(psin(1:ndim)**2))*sign(1._rk,psio(1)*psin(1))
!  psio=psin/eigenval_p
!  if(abs(eigenval_p-eigenval_sav).lt.eps) exit
!  eigenval_sav=eigenval_p
!enddo

! Eigenvectors is the matrix on input and the eigenvectors on output
eigenvectors=P_imp
call real_general_diagonalize_ow_ham(ndim,MDIM,eigenvectors,eigenvalues,eigenvalues_imag)
eigenval_p=eigenvalues(ndim)
eigenval_h=(1-eigenval_p+tau*E_T)/tau
psio(1:ndim)=eigenvectors(1:ndim,ndim)
if(psio(1).lt.0) psio(1:ndim)=-psio(1:ndim)
psio=psio/sum(abs(psio(1:ndim))) ! Since psio is already the product of 2 wavefns., normalize by the sum rather than the sqrt of the sum of squares
write(6,'(/,''eigenvalues of P_imp='',100f10.6)') (eigenvalues(i),i=1,ndim)
write(6,'(''eigenvalues of H_imp='',100f10.6)') ((1-eigenvalues(i)+tau*E_0)/tau,i=1,ndim)
write(6,'(/,''eigenval_p,eigenval_h='',9f10.6)') eigenval_p, eigenval_h
write (6,'(''psi0_psiG    ='',100f9.4)') (psio(i),i=1,ndim)

psi0_psiG=psio

! Sign-condition number for importance-sampled projector
sum_P_eig=0 ; sum_abs_P_eig=0
do j=1,ndim
  do i=1,ndim
    sum_abs_P_eig=sum_abs_P_eig+abs(P_imp(j,i)*psio(i))
  enddo
  sum_P_eig=sum_P_eig+abs(psio(j))
enddo
sign_condition_number=(sum_P_eig/sum_abs_P_eig)**(1/tau)
write(6,'(/,''imp-sampled: sum_P_eig,sum_abs_P_eig,sign_condition_number='',9f10.6)') sum_P_eig,sum_abs_P_eig,sign_condition_number

! Make sure tau is small enough that the elements of P_FN1 are all nonnegative
tau=tau_in
do j=1,ndim
  diagonal_dump=H_imp(j,j)
  do i=1,ndim
    if(i.ne.j .and. H_imp(i,j).gt.0) diagonal_dump=diagonal_dump+H_imp(i,j)
  enddo
  !write(6,'(''1/(diagonal_dump-E_0)'',9f10.6)') 1/(diagonal_dump-E_0), 1/(H_imp(j,j)-E_0)
  tau=min(tau,1/(diagonal_dump-E_0))
enddo
if(tau.lt.tau_in) then
  write(6,'(''Warning: tau_in='',f10.6,'' Reduce tau to'',f10.6)') tau_in, tau
  stop 'tau too large for P_FN1 to have all nonnegative elements'
endif

! Make FN1 projector
do i=1,ndim
  sign_flip=0
  do j=1,ndim
    if(j.ne.i) then
      if(P_imp(j,i).ge.0) then
        P_FN1(j,i)=P_imp(j,i)
      else
        sign_flip=sign_flip+P_imp(j,i)
        P_FN1(j,i)=0
      endif
    endif
  enddo
  P_FN1(i,i)=P_imp(i,i)+sign_flip
  if(P_FN1(i,i).lt.0) write(6,'(''Warning: for itimes, i='',2i4,'', P_FN1(i,i)='',es12.4)') itimes, i, P_FN1(i,i)
enddo
write (fmt,'(a,i2,a)') '(/,''P_FN1='',/,(', ndim, 'f9.4))'
write(6,fmt) ((P_FN1(j,i),i=1,ndim),j=1,ndim)

! Diagonalize to calculate eigevector and eigenvalues for P_FN1 and deduce eigenvalue for H_FN1.
! Eigenvectors is the matrix on input and the eigenvectors on output
eigenvectors=P_FN1
call real_general_diagonalize_ow_ham(ndim,MDIM,eigenvectors,eigenvalues,eigenvalues_imag)
eigenval_p=eigenvalues(ndim)
eigenval_h=(1-eigenval_p+tau*E_0)/tau
E_FN1=eigenval_h
E_T=E_FN1
psio(1:ndim)=eigenvectors(1:ndim,ndim)
if(psio(1).lt.0) psio(1:ndim)=-psio(1:ndim)
psio=psio/sum(abs(psio(1:ndim))) ! Since psio is already the product of 2 wavefns., normalize by the sum rather than the sqrt of the sum of squares
psiG_psiFN(1:ndim)=psio(1:ndim)
psiFN1(1:ndim)=psiG_psiFN(1:ndim)/psiG(1:ndim)
psiFN1(1:ndim)=psiFN1(1:ndim)/sqrt(sum(psiFN1(1:ndim)**2))
overlap_G_FN1=dot_product(psiG(1:ndim),psiFN1(1:ndim))
write(6,'(/,''eigenvalues of P_FN1='',100f10.6)') (eigenvalues(i),i=1,ndim)
write(6,'(''eigenvalues of H_FN1='',100f10.6)') ((1-eigenvalues(i)+tau*E_0)/tau,i=1,ndim)
write(6,'(/,''eigenval_p,eigenval_h='',9f10.6)') eigenval_p, eigenval_h
write (6,'(''psiFN1_psiG = '',100f9.4)') (psio(i),i=1,ndim)
write (6,'(''psiFN1      = '',100f9.4)') (psiFN1(i),i=1,ndim)
write (6,'(''overlap_G_FN1= '',100f9.4)') overlap_G_FN1
!write (6,'(''norm_psiFN1      = '',100f19.14)') sum(psiFN1(1:ndim)**2)

!overlap_FN1=overlap(ndim,psi0_psiG,psio)
overlap_FN1=overlap(ndim,psi0,psiFN1)

! Redefine diagonal of P_FN1 so that the largest eigenvalue of P_FN1 is 1
do i=1,ndim
  P_FN1(i,i)=P_FN1(i,i)+tau*(E_FN1-E_0)
enddo


! Calculate local projector, local abs projector, and, proposal matrix for Metropolis allowing node crossing
!do i=1,ndim
!  P_L(i)=0
!  P_abs_L(i)=0
!  do j=1,ndim
!    P_L(i)=P_L(i)+P_imp(j,i)
!    P_abs_L(i)=P_abs_L(i)+abs(P_imp(j,i))
!  enddo
!  do j=1,ndim
!    P_propos(j,i)=abs(P_imp(j,i))/P_abs_L(i)
!  enddo
!enddo
!write (6,'(/,''P_L     ='',100f9.4)') (P_L(i),i=1,ndim)
!write (6,'(''P_abs_L ='',100f9.4)') (P_abs_L(i),i=1,ndim)
!write (fmt,'(a,i2,a)') '(/,''P_propos='',/,(', ndim, 'f9.4))'
!write(6,fmt) ((P_propos(j,i),i=1,ndim),j=1,ndim)
!
!! Calculate Markov matrix for Metropolis allowing node crossing
!do i=1,ndim
!  rejected=0
!  do j=1,ndim
!    if(P_propos(i,j)*psiG(j)**2/(P_propos(j,i)*psiG(i)**2).ge.1) then ! If acceptance prob >= 1
!      P_markov(j,i)=P_propos(j,i)
!    else
!      rejected=rejected+P_propos(j,i)-P_propos(i,j)*(psiG(j)/psiG(i))**2
!      P_markov(j,i)=P_propos(i,j)*(psiG(j)/psiG(i))**2
!    endif
!  enddo
!  P_markov(i,i)=P_markov(i,i)+rejected
!  if(P_markov(i,i).lt.0) write(6,'(''Warning: for i='',i4,'', P_markov(i,i)='',es12.4)') i, P_markov(i,i)
!enddo
!write (fmt,'(a,i2,a)') '(/,''P_markov='',/,(', ndim, 'f9.4))'
!write(6,fmt) ((P_markov(j,i),i=1,ndim),j=1,ndim)
!
!! Calculate eigenvalue for P_markov to make sure we get the variational energy for psiG
!! Eigenvectors is the matrix on input and the eigenvectors on output
!eigenvectors=P_markov
!call real_general_diagonalize_ow_ham(ndim,MDIM,eigenvectors,eigenvalues,eigenvalues_imag)
!eigenval_p=eigenvalues(ndim)
!eigenval_h=(1-eigenval_p+tau*E_T)/tau
!psio(1:ndim)=eigenvectors(1:ndim,ndim)
!if(psio(1).lt.0) psio(1:ndim)=-psio(1:ndim)
!write(6,'(/,''eigenvalues of P_markov='',100f10.6)') (eigenvalues(i),i=1,ndim)
!write(6,'(''eigenvalues of H_markov='',100f10.6)') ((1-eigenvalues(i)+tau*E_0)/tau,i=1,ndim)
!write(6,'(/,''eigenval_p,eigenval_h='',9f10.6)') eigenval_p, eigenval_h
!write (6,'(''psi for P_mar='',100f9.4)') (psio(i),i=1,ndim)
!write (6,'(''psiG**2      ='',100f9.4)') (psiG(i)**2/sqrt(sum(psiG(1:ndim)**4)),i=1,ndim)
!write (6,'(''Above 2 lines should be equal since Markov matrix has dominant state psiG**2'')')
!------------------------------------------

! Calculate local energy, local projector, local abs projector, and, proposal matrix for Metropolis not allowing node crossing
! Note we P_L would be the same using P_FN1 or P_imp if E_T is the same.
! However, we use E_T=E_0 for P_imp and E_T=E_FN1 for P_FN1, so P_L is different by a constant for these 2 options.
! I use P_imp for defining P_abs_L because P_abs_L is used to construct P_propos from P_imp, and,
! I use P_FN1 for defining P_L because P_L is used to calculate psiFN1_by_psiG.
do i=1,ndim
  E_L(i)=0
  P_L(i)=0
  P_abs_L(i)=0
  do j=1,ndim
    E_L(i)=E_L(i)+H_imp(j,i)
    !P_L(i)=P_L(i)+P_imp(j,i)
    P_L(i)=P_L(i)+P_FN1(j,i)
    if(P_imp(j,i).gt.0) P_abs_L(i)=P_abs_L(i)+P_imp(j,i)
  enddo
  do j=1,ndim
    if(P_imp(j,i).gt.0) then
      P_propos(j,i)=abs(P_imp(j,i))/P_abs_L(i)
    else
      P_propos(j,i)=0
    endif
  enddo
enddo
write(6,'(/,''E_0, E_G, E_FN1, E_G-E_0, E_FN1-E_0, rms fluctuation of E_L'',9f10.5)') E_0, E_G, E_FN1, E_G-E_0, E_FN1-E_0, rms_deviation(ndim,E_L,psiG_psiFN)
write (6,'(''E_L     ='',100f9.4)') (E_L(i),i=1,ndim)
write (6,'(''P_L     ='',100f9.4)') (P_L(i),i=1,ndim)
write (6,'(''P_abs_L ='',100f9.4)') (P_abs_L(i),i=1,ndim)
write (fmt,'(a,i2,a)') '(/,''P_propos='',/,(', ndim, 'f9.4))'
write(6,fmt) ((P_propos(j,i),i=1,ndim),j=1,ndim)
write(6,'(''sum P_propos'',9f10.6)') (sum(P_propos(:,i)),i=1,ndim)

write(6,'(''E_FN1='',9f10.6)') E_FN1, sum(psio(1:ndim)*E_L(1:ndim))/sum(psio(1:ndim))
if(abs(E_FN1-sum(psio(1:ndim)*E_L(1:ndim))/sum(psio(1:ndim))).gt.eps2) then
  write(6,'(''E_FN1, sum(psio(1:ndim)*E_L(1:ndim))/sum(psio(1:ndim))'',9f10.6)') E_FN1, sum(psio(1:ndim)*E_L(1:ndim))/sum(psio(1:ndim))
  stop 'growth and mixed estimators for FN1 are not consistent'
endif

! Calculate script psiFN1_by_psiG = psiG_psiFN/(psiG^2)
! Note: the psiFN1_by_psiG here is the script W = P_L * W in Eqs. 102-105 of my notes.
!norm=sum(psiG(1:ndim)**2*P_L(1:ndim))
!psiFN1_by_psiG(1:ndim) = psiG_psiFN(1:ndim)*norm/(psiG(1:ndim)**2*P_L(1:ndim))
norm=sum(psiG(1:ndim)**2)
psiFN1_by_psiG(1:ndim) = psiG_psiFN(1:ndim)*norm/(psiG(1:ndim)**2)

write (6,'(''psiG**2=      '',100f9.4)') (psiG(i)**2,i=1,ndim)
write (6,'(''psiG**2 * P_L='',100f9.4)') (psiG(i)**2 * P_L(i),i=1,ndim)
write (6,'(''psiG_psiFN=   '',100f9.4)') (psiG_psiFN(i),i=1,ndim)
write (6,'(''P_L=          '',100f9.4)') (P_L(i),i=1,ndim)
write (6,'(''psiFN1_by_psiG=       '',100f9.4)') (psiFN1_by_psiG(i),i=1,ndim)

!write(6,'(2f9.4)') (P_L(i),psiFN1_by_psiG(i),i=1,ndim)

! Calculate Markov matrix for Metropolis not allowing node crossing
do i=1,ndim
  rejected=0
  do j=1,ndim
    if(P_propos(j,i)*psiG(i)**2.ne.0) then
      if(P_propos(i,j)*psiG(j)**2/(P_propos(j,i)*psiG(i)**2).ge.1) then ! If acceptance prob >= 1
        P_markov(j,i)=P_propos(j,i)
      else
        rejected=rejected+P_propos(j,i)-P_propos(i,j)*(psiG(j)/psiG(i))**2
        P_markov(j,i)=P_propos(i,j)*(psiG(j)/psiG(i))**2
      endif
    else
      P_markov(j,i)=0
    endif
  enddo
  P_markov(i,i)=P_markov(i,i)+rejected
  if(P_markov(i,i).lt.0) write(6,'(''Warning: for itimes, i='',2i4,'', P_markov(i,i)='',es12.4)') itimes, i, P_markov(i,i)
enddo
write (fmt,'(a,i2,a)') '(/,''P_markov='',/,(', ndim, 'f9.4))'
write(6,fmt) ((P_markov(j,i),i=1,ndim),j=1,ndim)

! Calculate eigenvalue/vector for P_markov to make sure we get the variational energy for psiG and the distribution psiG^2
! Eigenvectors is the matrix on input and the eigenvectors on output
eigenvectors=P_markov
call real_general_diagonalize_ow_ham(ndim,MDIM,eigenvectors,eigenvalues,eigenvalues_imag)
eigenval_p=eigenvalues(ndim)
eigenval_h=(1-eigenval_p+tau*E_T)/tau
psio(1:ndim)=eigenvectors(1:ndim,ndim)
if(psio(1).lt.0) psio(1:ndim)=-psio(1:ndim)
psio=psio/sum(abs(psio(1:ndim))) ! Since psio is already the product of 2 wavefns., normalize by the sum rather than the sqrt of the sum of squares
write(6,'(/,''eigenvalues of P_markov='',100f10.6)') (eigenvalues(i),i=1,ndim)
write(6,'(''eigenvalues of H_markov='',100f10.6)') ((1-eigenvalues(i)+tau*E_0)/tau,i=1,ndim)
write(6,'(/,''eigenval_p,eigenval_h='',9f10.6)') eigenval_p, eigenval_h
!write (6,'(''psi for P_mar='',100f9.4)') (psio(i),i=1,ndim)
!write (6,'(''psiG**2      ='',100f9.4)') (psiG(i)**2/sqrt(sum(psiG(1:ndim)**4)),i=1,ndim)
psio(1:ndim)=psio(1:ndim)/sum(psio(1:ndim)) ! Since psio is already the product of 2 wavefns., normalize by the sum rather than the sqrt of the sum of squares
write (6,'(''psi for P_mar='',100f9.4)') (psio(i),i=1,ndim)
write (6,'(''psiG**2      ='',100f9.4)') (psiG(i)**2,i=1,ndim)
write (6,'(''Above 2 lines should be equal since Markov matrix has dominant state psiG**2'')')
write (6,'(''The exception to this occurs when there is no sign problem, but some components of psiG have the wrong sign, in which case one has 2 decoupled spaces with degenerate dominant eigenvalues'')')
if(dot_product(psio(1:ndim)-psiG(1:ndim)**2,psio(1:ndim)-psiG(1:ndim)**2).gt.eps2) then
  write (6,*) dot_product(psio(1:ndim)-psiG(1:ndim)**2,psio(1:ndim)-psiG(1:ndim)**2)
  stop 'eigenvec of Markov Matrix should be psiG**2, but is not'
endif

! Define P_FN3
do j=1,ndim
  do i=1,ndim
    P_FN3(j,i)=P_markov(j,i)*P_L(i)
    !P_FN3(j,i)=P_markov(j,i)*sqrt(P_L(i)*P_L(j))
    !P_FN3(j,i)=P_markov(j,i)*P_L(j)
  enddo
enddo
write (fmt,'(a,i2,a)') '(/,''P_FN3='',/,(', ndim, 'f9.4))'
write(6,fmt) ((P_FN3(j,i),i=1,ndim),j=1,ndim)

! Diagonalize to calculate eigevector and eigenvalues for P_FN3 and deduce eigenvalue for H_FN3.
! Eigenvectors is the matrix on input and the eigenvectors on output
eigenvectors=P_FN3
call real_general_diagonalize_ow_ham(ndim,MDIM,eigenvectors,eigenvalues,eigenvalues_imag)
eigenval_p=eigenvalues(ndim)
eigenval_h=(1-eigenval_p+tau*E_FN1)/tau
E_FN3=eigenval_h
psio(1:ndim)=eigenvectors(1:ndim,ndim)
if(psio(1).lt.0) psio(1:ndim)=-psio(1:ndim)
psio=psio/sum(abs(psio(1:ndim))) ! Since psio is already the product of 2 wavefns., normalize by the sum rather than the sqrt of the sum of squares
write(6,'(''E_FN3='',9f10.6)') E_FN3, sum(psio(1:ndim)*E_L(1:ndim))/sum(psio(1:ndim))
psiFN3(1:ndim)=psio(1:ndim)/psiG(1:ndim)
psiFN3(1:ndim)=psiFN3(1:ndim)/sqrt(sum(psiFN3(1:ndim)**2))
! Warning: tmp
!if(abs(E_FN3-sum(psio(1:ndim)*E_L(1:ndim))/sum(psio(1:ndim))).gt.eps2) stop 'growth and mixed estimators for FN3 are not consistent'
write(6,'(/,''eigenvalues of P_FN3='',100f10.6)') (eigenvalues(i),i=1,ndim)
write(6,'(''eigenvalues of H_FN3='',100f10.6)') ((1-eigenvalues(i)+tau*E_FN1)/tau,i=1,ndim)
write(6,'(/,''eigenval_p,eigenval_h='',9f10.6)') eigenval_p, eigenval_h
write (6,'(''psiFN3_psiG  ='',100f9.4)') (psio(i),i=1,ndim)
write (6,'(''psiFN3      = '',100f9.4)') (psiFN3(i),i=1,ndim)

!overlap_FN3=overlap(ndim,psi0_psiG,psio)
overlap_FN3=overlap(ndim,psi0,psiFN3)

! Define P_FN4
do j=1,ndim
  do i=1,ndim
    P_FN4(j,i)=P_propos(j,i)*P_L(i)
    !P_FN4(j,i)=P_propos(j,i)*sqrt(P_L(i)*P_L(j))
    !P_FN4(j,i)=P_propos(j,i)*P_L(j)
  enddo
enddo
write (fmt,'(a,i2,a)') '(/,''P_FN4='',/,(', ndim, 'f9.4))'
write(6,fmt) ((P_FN4(j,i),i=1,ndim),j=1,ndim)

! Diagonalize to calculate eigevector and eigenvalues for P_FN4 and deduce eigenvalue for H_FN4.
! Eigenvectors is the matrix on input and the eigenvectors on output
eigenvectors=P_FN4
call real_general_diagonalize_ow_ham(ndim,MDIM,eigenvectors,eigenvalues,eigenvalues_imag)
eigenval_p=eigenvalues(ndim)
eigenval_h=(1-eigenval_p+tau*E_FN1)/tau
E_FN4=eigenval_h
psio(1:ndim)=eigenvectors(1:ndim,ndim)
if(psio(1).lt.0) psio(1:ndim)=-psio(1:ndim)
psio=psio/sum(abs(psio(1:ndim))) ! Since psio is already the product of 2 wavefns., normalize by the sum rather than the sqrt of the sum of squares
write(6,'(''E_FN4='',9f10.6)') E_FN4, sum(psio(1:ndim)*E_L(1:ndim))/sum(psio(1:ndim))
psiFN4(1:ndim)=psio(1:ndim)/psiG(1:ndim)
psiFN4(1:ndim)=psiFN4(1:ndim)/sqrt(sum(psiFN4(1:ndim)**2))
if(abs(E_FN4-sum(psio(1:ndim)*E_L(1:ndim))/sum(psio(1:ndim))).gt.eps2) stop 'growth and mixed estimators for FN4 are not consistent'
write(6,'(/,''eigenvalues of P_FN4='',100f10.6)') (eigenvalues(i),i=1,ndim)
write(6,'(''eigenvalues of H_FN4='',100f10.6)') ((1-eigenvalues(i)+tau*E_FN1)/tau,i=1,ndim)
write(6,'(/,''eigenval_p,eigenval_h='',9f10.6)') eigenval_p, eigenval_h
write (6,'(''psiFN4_psiG  ='',100f9.4)') (psio(i),i=1,ndim)
write (6,'(''psiFN4      = '',100f9.4)') (psiFN4(i),i=1,ndim)

!overlap_FN4=overlap(ndim,psi0_psiG,psio)
overlap_FN4=overlap(ndim,psi0,psiFN4)

!! Define P_FN5
!do i=1,ndim
!  sum_pos=0
!  sum_neg=0
!! do j=1,ndim
!!   if(P_imp(j,i).gt.0) sum_pos=sum_pos+P_imp(j,i)
!!   if(P_imp(j,i).lt.0) sum_neg=sum_neg+P_imp(j,i)
!! enddo
!  do j=1,ndim
!    if(j.ne.i) then
!      if(-H_imp(j,i).gt.0) sum_pos=sum_pos-H_imp(j,i)
!      if(-H_imp(j,i).lt.0) sum_neg=sum_neg-H_imp(j,i)
!    endif
!  enddo
!  ratio=sum_neg/sum_pos
!  if(abs(sum_neg).gt.sum_pos) write(6,'(''Warning: for FN5, sum_neg,sum_pos,ratio'',9f10.6)') sum_neg,sum_pos,ratio
!  do j=1,ndim
!    if(j.eq.i) then
!      P_FN5(j,i)=P_imp(j,i)
!    else
!      if(P_imp(j,i).gt.0) then
!        P_FN5(j,i)=P_imp(j,i) - tau*H_imp(j,i)*ratio
!      else
!        P_FN5(j,i)=0
!      endif
!    endif
!  enddo
!enddo
!write (fmt,'(a,i2,a)') '(/,''P_FN5='',/,(', ndim, 'f9.4))'
!write(6,fmt) ((P_FN5(j,i),i=1,ndim),j=1,ndim)
!
!! Diagonalize to calculate eigevector and eigenvalues for P_FN5 and deduce eigenvalue for H_FN5.
!! Eigenvectors is the matrix on input and the eigenvectors on output
!eigenvectors=P_FN5
!call real_general_diagonalize_ow_ham(ndim,MDIM,eigenvectors,eigenvalues,eigenvalues_imag)
!eigenval_p=eigenvalues(ndim)
!eigenval_h=(1-eigenval_p+tau*E_0)/tau
!E_FN5=eigenval_h
!psio(1:ndim)=eigenvectors(1:ndim,ndim)
!if(psio(1).lt.0) psio(1:ndim)=-psio(1:ndim)
!psio=psio/sum(abs(psio(1:ndim))) ! Since psio is already the product of 2 wavefns., normalize by the sum rather than the sqrt of the sum of squares
!write(6,'(''E_FN5='',9f10.6)') E_FN5, sum(psio(1:ndim)*E_L(1:ndim))/sum(psio(1:ndim))
!if(abs(E_FN5-sum(psio(1:ndim)*E_L(1:ndim))/sum(psio(1:ndim))).gt.eps2) stop 'growth and mixed estimators for FN5 are not consistent'
!write(6,'(/,''eigenvalues of P_FN5='',100f10.6)') (eigenvalues(i),i=1,ndim)
!write(6,'(''eigenvalues of H_FN5='',100f10.6)') ((1-eigenvalues(i)+tau*E_0)/tau,i=1,ndim)
!write(6,'(/,''eigenval_p,eigenval_h='',9f10.6)') eigenval_p, eigenval_h
!write (6,'(''psiFN5_psiG  ='',100f9.4)') (psio(i),i=1,ndim)
!
!overlap_FN5=overlap(ndim,psi0_psiG,psio)

! Define P_FN5
do i=1,ndim
  sum_pos=0
  do j=1,ndim
    if(j.ne.i) then
      if(H_imp(j,i).gt.0) sum_pos=sum_pos+H_imp(j,i)
    endif
  enddo
  write(6,'(''FN5: sum_pos'',9f10.6)') sum_pos
  do j=1,ndim
    if(j.eq.i) then
      !P_FN5(j,i)=P_imp(j,i)/(1+tau*sum_pos+(tau**2)*(sum_pos**2-(E_0-H_imp(j,i))*sum_pos))
      !P_FN5(j,i)=1/(1-tau*(E_0-H_imp(j,i)-sum_pos)+(tau*(E_0-H_imp(j,i)-sum_pos))**2)
      !P_FN5(j,i)=1/(1-tau*(E_0-H_imp(j,i)-sum_pos)*(1-tau*(E_0-H_imp(j,i)-sum_pos)))
      P_FN5(j,i)=1/(1-tau*(E_0-H_imp(j,i)-sum_pos)*(1-tau*(E_0-H_imp(j,i)-sum_pos)*(1-tau*(E_0-H_imp(j,i)-sum_pos))))
    else
      if(P_imp(j,i).gt.0) then
        P_FN5(j,i)=P_imp(j,i)
      else
        P_FN5(j,i)=0
      endif
    endif
  enddo
enddo
write (fmt,'(a,i2,a)') '(/,''P_FN5='',/,(', ndim, 'f9.4))'
write(6,fmt) ((P_FN5(j,i),i=1,ndim),j=1,ndim)

! Diagonalize to calculate eigevector and eigenvalues for P_FN5 and deduce eigenvalue for H_FN5.
! Eigenvectors is the matrix on input and the eigenvectors on output
eigenvectors=P_FN5
call real_general_diagonalize_ow_ham(ndim,MDIM,eigenvectors,eigenvalues,eigenvalues_imag)
eigenval_p=eigenvalues(ndim)
eigenval_h=(1-eigenval_p+tau*E_0)/tau
E_FN5=eigenval_h
psio(1:ndim)=eigenvectors(1:ndim,ndim)
if(psio(1).lt.0) psio(1:ndim)=-psio(1:ndim)
psio=psio/sum(abs(psio(1:ndim))) ! Since psio is already the product of 2 wavefns., normalize by the sum rather than the sqrt of the sum of squares
write(6,'(''E_FN5='',9f10.6)') E_FN5, sum(psio(1:ndim)*E_L(1:ndim))/sum(psio(1:ndim))
psiFN5(1:ndim)=psio(1:ndim)/psiG(1:ndim)
psiFN5(1:ndim)=psiFN5(1:ndim)/sqrt(sum(psiFN5(1:ndim)**2))
! Warning: tmp 100*
!if(abs(E_FN5-sum(psio(1:ndim)*E_L(1:ndim))/sum(psio(1:ndim))).gt.100*eps2) stop 'growth and mixed estimators for FN5 are not consistent'
write(6,'(/,''eigenvalues of P_FN5='',100f10.6)') (eigenvalues(i),i=1,ndim)
write(6,'(''eigenvalues of H_FN5='',100f10.6)') ((1-eigenvalues(i)+tau*E_0)/tau,i=1,ndim)
write(6,'(/,''eigenval_p,eigenval_h='',9f10.6)') eigenval_p, eigenval_h
write (6,'(''psiFN5_psiG  ='',100f9.4)') (psio(i),i=1,ndim)
write (6,'(''psiFN5      = '',100f9.4)') (psiFN5(i),i=1,ndim)

!overlap_FN5=overlap(ndim,psi0_psiG,psio)
overlap_FN5=overlap(ndim,psi0,psiFN5)

! Print the fraction of the energy missing in VMC recovered and 1-overlap for various FN algorithms
write(6,'(/,''E_G, E_FN1, E_FN3, E_FN4, E_0, recover_FN1, recover_FN3, recover_FN4, 1-ovlp_FN1, 1-ovlp_FN3, 1-ovlp_FN4, sign_cond_num='',5f10.6,3f6.3,3es7.0,f6.3)') E_G, E_FN1, E_FN3, E_FN4, E_0, (E_G-E_FN1)/(E_G-E_0), (E_G-E_FN3)/(E_G-E_0), (E_G-E_FN4)/(E_G-E_0), 1-overlap_FN1, 1-overlap_FN3, 1-overlap_FN4, sign_condition_number
write(1,'(5f10.6,3f10.3,x,3es9.1,f7.3)') E_G, E_FN1, E_FN3, E_FN4, E_0, (E_G-E_FN1)/(E_G-E_0), (E_G-E_FN3)/(E_G-E_0), (E_G-E_FN4)/(E_G-E_0), 1-overlap_FN1, 1-overlap_FN3, 1-overlap_FN4, sign_condition_number
!write(1,'(5f10.6,4f10.3,x,3es9.1,f7.3)') E_G, E_FN1, E_FN3, E_FN4, E_0, (E_G-E_FN1)/(E_G-E_0), (E_G-E_FN3)/(E_G-E_0), (E_G-E_FN4)/(E_G-E_0), (E_G-E_FN5)/(E_G-E_0), 1-overlap_FN1, 1-overlap_FN3, 1-overlap_FN4, sign_condition_number

! Verify that the proposal part of the FN1 projector has eigenstate psiG^2 * P_L
! Note this proposal matrix is different from the one before, because now the negative elements in a column have been added to the diagonal
do i=1,ndim
  do j=1,ndim
    P_propos(j,i)=P_FN1(j,i)/P_L(i)
  enddo
enddo
write (fmt,'(a,i2,a)') '(/,''P_propos2='',/,(', ndim, 'f9.4))'
write(6,fmt) ((P_propos(j,i),i=1,ndim),j=1,ndim)

! Calculate eigenvalue/vector for P_FN1_propos to make sure largest eigenvalue is 1 and the distribution psiG^2 * P_L
! Eigenvectors is the matrix on input and the eigenvectors on output
eigenvectors=P_propos
call real_general_diagonalize_ow_ham(ndim,MDIM,eigenvectors,eigenvalues,eigenvalues_imag)
eigenval_p=eigenvalues(ndim)
eigenval_h=(1-eigenval_p+tau*E_T)/tau
psio(1:ndim)=eigenvectors(1:ndim,ndim)
if(psio(1).lt.0) psio(1:ndim)=-psio(1:ndim)
psio=psio/sum(abs(psio(1:ndim))) ! Since psio is already the product of 2 wavefns., normalize by the sum rather than the sqrt of the sum of squares
write(6,'(/,''eigenvalues of P_FN1_propos='',100f10.6)') (eigenvalues(i),i=1,ndim)
psio(1:ndim)=psio(1:ndim)/sum(psio(1:ndim)) ! Since psio is already the product of 2 wavefns., normalize by the sum rather than the sqrt of the sum of squares
write (6,'(''psi for P_FN1_propos='',100f9.4)') (psio(i),i=1,ndim)
!write (6,'(''psiG**2 * P_L       ='',100f9.4)') (psiG(i)**2 * P_L(i)/sqrt(sum((psiG(1:ndim)**2 * P_L(1:ndim))**2)),i=1,ndim)
write (6,'(''psiG**2 * P_L       ='',100f9.4)') (psiG(i)**2 * P_L(i)/sum(psiG(1:ndim)**2 * P_L(1:ndim)),i=1,ndim)
write (6,'(''Above 2 lines should be equal since proposal part of P_FN1 has dominant state psiG**2 * P_L'')')
write (6,'(''The exception to this occurs when there is no sign problem, but some components of psiG have the wrong sign, in which case one has 2 decoupled spaces with degenerate dominant eigenvalues'')')
if(dot_product(psio(1:ndim)-psiG(1:ndim)**2 * P_L(1:ndim)/sum(psiG(1:ndim)**2 * P_L(1:ndim)), psio(1:ndim)-psiG(1:ndim)**2 * P_L(1:ndim)/sum(psiG(1:ndim)**2 * P_L(1:ndim))).gt.eps2) then
  stop 'eigenvec of proposal part of P_FN1 should be psiG**2 * P_L, but is not'
endif



! Now we calculate gradients of VMC and FN-DMC energies

! Old: Take the wavefn parameters to be just the coefficients of each basis state
! Old: Then the derivative wrt a parameter is 1 on that state and zero on other states

! New: All parameters are as above, but 



psiG_der = 0
psiG_2der = 0

!a_3 = psiG(iparm+1)

! New: Treat the parameters as above, but take the iparm(in this case, iparm=2) parameter to be:
!    a_3 exp(a_2)
! So that the derivative is
!    a_3 exp(a_2)

!Calculate a_2
!a_3 = psiG(iparm+1)

!Calculate derivative wrt every parameter
!All the linear parameters will have a 1 when their parameter and state line up
do i=1,ndim
   if(i.ne.iparm) then
      psiG_der(i,i) = 1
   endif
enddo

!The first argument is the state, the second argument is the derivative
psiG_der(iparm,iparm)   = psiG(iparm+1)*exp(a_2)
psiG_der(iparm,iparm+1) = exp(a_2)

!Calculate 2nd derivatives on every site
!The second derivatives are zero for every parameter
!except a special few:
! d^2 Y/ da_2^2 = a_3 exp(a_2) [0 1 0]
! d^2 Y/ da_2 da_3 = exp(a_2) [0 1 0]
!
!And, these are only defined on site iaprm
!
psiG_2der(iparm,iparm,iparm) = psiG(iparm+1)*exp(a_2)
psiG_2der(iparm,iparm,iparm+1) = exp(a_2)
psiG_2der(iparm,iparm+1,iparm) = exp(a_2)

! Warning: tmp2
!psiG_der(1:ndim) = psiG_der(1:ndim) - dot_product(psiG(1:ndim),psiG_der(1:ndim))*psiG(1:ndim)

!Calculate FD gradients
!*_v means vector (as opposed to the old version, where the gradient is a scalar)
gradient_VMC_FD   = 0
gradient_FN_FD = 0
gradient_FN1_exact   = 0

! Calculate E_{L,iparm}(k), i.e., derivative of the local energy on state i wrt the iparm parameter
! E_{L,iparm} = H psiG_der(iparm)/psiG - E_L psiG_der(iparm)/psiG, i.e.,
! E_{L,i}(k) = \sum_j H(i,j) psiG_der(j.k)/psiG(i) - E_L(i) psiG_der(i,k)/psiG(i)
! i is derivative, k is state
do k=1,ndim
   do i=1,ndim
      E_L_der(k,i)=(dot_product(h(k,1:ndim),psiG_der(1:ndim,i))-E_L(k)*psiG_der(k,i)) / psiG(k)
      P_L_der(k,i)=-tau*E_L_der(k,i)
      P_L_der_by_P_L(k,i) = P_L_der(k,i)/P_L(k)
   enddo
enddo

! Calculate the 2nd derivative wrt i,j of E_L
do k=1,ndim
   do i=1,ndim
      do j=1,ndim
         E_L_2der(k,i,j) = (dot_product(h(k,1:ndim),psiG_2der(1:ndim,i,j))-E_L(k)*psiG_2der(k,i,j)) &
              &/psiG(k) - (psiG_der(k,i)*E_L_der(k,j)/psiG(k) + &
              &psiG_der(k,j)*E_L_der(k,i)/psiG(k))

         P_L_2der(k,i,j)=-tau*E_L_2der(k,i,j)
         P_L_2der_by_P_L(k,i,j)=P_L_2der(k,i,j)/P_L(k)
      enddo
   enddo
enddo

!Calculate W_der_by_W

!Set the initial projector to I
PP_FN1 = 0
do k=1,ndim
   PP_FN1(k,k) = 1
enddo

!Project many times
do iter=1,20000
   call mat_mul2(ndim,MDIM,P_FN1,PP_FN1)
enddo

do i=1,ndim!derivative i
   !Calculate tmp_vec1
   call vec_vec_elem_mul(ndim,P_L_der_by_P_L(1:ndim,i),psiG_psiFN,tmp_vec1)
   !Project tmp_vec1
   call mat_vec_mul(ndim,MDIM,PP_FN1,tmp_vec1,tmp_vec2)
   W_der_by_W(i) = dot_product(E_L(1:ndim)-E_FN1,tmp_vec2(1:ndim))
   call mat_vec_mul(ndim,MDIM,PP_FN1,tmp_vec1,tmp_vec2)
   !Full means not dotted with anything; we need the 'full' version for W_ij/W
   W_der_by_W_full(1:ndim,i) = tmp_vec2(1:ndim)
enddo

!get P+P^2+P^3+...
do i=1,ndim
   do j=1,ndim
      PP_FN1_minus_I(i,j) = PP_FN1(i,j)
   enddo
enddo

do k=1,ndim
   PP_FN1_minus_I(k,k) = PP_FN1(k,k) - 1
enddo


!calculate W_2der_by_W
do i=1,ndim!derivative i
   do j=1,ndim!derivative j
      tmp_vec1   = 0
      tmp_vec1_i = 0
      tmp_vec1_j = 0
      !Calculate tmp_vec1
      call vec_vec_elem_mul(ndim,P_L_2der_by_P_L(1:ndim,i,j),psiG_psiFN,tmp_vec1)
      call vec_vec_elem_mul(ndim,P_L_der_by_P_L(1:ndim,i),psiG_psiFN,tmp_vec1_i)
      call vec_vec_elem_mul(ndim,P_L_der_by_P_L(1:ndim,j),psiG_psiFN,tmp_vec1_j)

      call mat_vec_mul(ndim,MDIM,PP_FN1,tmp_vec1,tmp_vec2)

      W_2der_by_W_full(1:ndim,i,j) = tmp_vec2(1:ndim)

      !Do mixed w_i w_j term
      !Project P_L_i term
      call mat_vec_mul(ndim,MDIM,PP_FN1_minus_I,tmp_vec1_i,tmp_vec2)

      !Elementwise multiplication with P_L_j
      call vec_vec_elem_mul(ndim,P_L_der_by_P_L(1:ndim,j),tmp_vec2,tmp_vec3)

      !Project this term
      call mat_vec_mul(ndim,MDIM,PP_FN1,tmp_vec3,tmp_vec4)

      !Add to W_2der_by_W_full
      W_2der_by_W_full(1:ndim,i,j) = W_2der_by_W_full(1:ndim,i,j) + tmp_vec4(1:ndim)

      !Do mixed w_j w_i term
      !Project P_L_j term
      call mat_vec_mul(ndim,MDIM,PP_FN1_minus_I,tmp_vec1_j,tmp_vec2)

      !Elementwise multiplication with P_L_i
      call vec_vec_elem_mul(ndim,P_L_der_by_P_L(1:ndim,i),tmp_vec2,tmp_vec3)

      !Project this term
      call mat_vec_mul(ndim,MDIM,PP_FN1,tmp_vec3,tmp_vec4)

      !Add to W_2der_by_W_full
      W_2der_by_W_full(1:ndim,i,j) = W_2der_by_W_full(1:ndim,i,j) + tmp_vec4(1:ndim)

      W_2der_by_W(i,j) = dot_product(E_L(1:ndim)-E_FN1,W_2der_by_W_full(1:ndim,i,j))
   enddo
enddo

write(6,*)'W_2der_by_W'
do j=1,ndim
   write(6,'(100g14.6)') (W_2der_by_W(i,j),i=1,ndim)
enddo

psig_increm = psiG
!Increment the linear parameters
do i=1,ndim
   if((i.ne.iparm).and.(i.ne.(iparm+1))) then
      !These are the strictly linear parameters
      !I know that they will be 1, but this is a check.
      psiG_increm(i) = psiG(i) + eps_FD
      !Get the incremented energies
      call calc_energy(psiG_increm,E_G_increm,E_FN1_increm)
      gradient_VMC_FD(i)   = (E_G_increm-E_G)/eps_FD
      gradient_FN_FD(i) = (E_FN1_increm-E_FN1)/eps_FD

      psiG_psiFN1_increm(1:ndim)=eigenvectors(1:ndim,ndim)
      if(psiG_psiFN1_increm(1).lt.0) psiG_psiFN1_increm(1:ndim)=-psiG_psiFN1_increm(1:ndim)
      psiG_psiFN1_increm=psiG_psiFN1_increm/sum(abs(psiG_psiFN1_increm(1:ndim))) ! Since psiG_psiFN1_increm is already the product of 2 wavefns., normalize by the sum rather than the sqrt of the sum of squares
      norm=sum(psiG_increm(1:ndim)**2)
      psiFN1_by_psiG_increm(1:ndim) = psiG_psiFN1_increm(1:ndim)*norm/(psiG_increm(1:ndim)**2)

      psiFN1_by_psiG_der(1:ndim)=(psiFN1_by_psiG_increm(1:ndim)-psiFN1_by_psiG(1:ndim))/eps_FD

      !Also calculate gradient_FN1_exact
      do j=1,ndim
         gradient_FN1_exact(i)=gradient_FN1_exact(i)+psiG_psiFN(j)*((2*(psiG_der(j,i)/psiG(j))+psiFN1_by_psiG_der(j)/psiFN1_by_psiG(j))*(E_L(j)-E_FN1)+E_L_der(j,i))
      enddo

      !Reset psiG_increm

      psiG_increm(i) = psiG(i)

   endif
enddo

psiG_increm = psiG

!This is the nonlinear parameter, exp(a_2)
psiG_increm(iparm) = psiG(iparm+1)*exp(a_2+eps_FD)

call calc_energy(psiG_increm,E_G_increm,E_FN1_increm)
gradient_VMC_FD(iparm)   = (E_G_increm-E_G)/eps_FD
gradient_FN_FD(iparm) = (E_FN1_increm-E_FN1)/eps_FD

psiG_psiFN1_increm(1:ndim)=eigenvectors(1:ndim,ndim)
if(psiG_psiFN1_increm(1).lt.0) psiG_psiFN1_increm(1:ndim)=-psiG_psiFN1_increm(1:ndim)
psiG_psiFN1_increm=psiG_psiFN1_increm/sum(abs(psiG_psiFN1_increm(1:ndim))) ! Since psiG_psiFN1_increm is already the product of 2 wavefns., normalize by the sum rather than the sqrt of the sum of squares
norm=sum(psiG_increm(1:ndim)**2)
psiFN1_by_psiG_increm(1:ndim) = psiG_psiFN1_increm(1:ndim)*norm/(psiG_increm(1:ndim)**2)

psiFN1_by_psiG_der(1:ndim)=(psiFN1_by_psiG_increm(1:ndim)-psiFN1_by_psiG(1:ndim))/eps_FD

do j=1,ndim
   gradient_FN1_exact(iparm)=gradient_FN1_exact(iparm)+psiG_psiFN(j)*((2*(psiG_der(j,iparm)/psiG(j))+psiFN1_by_psiG_der(j)/psiFN1_by_psiG(j))*(E_L(j)-E_FN1)+E_L_der(j,iparm))
enddo

!Reset psiG_increm
psiG_increm(iparm) = psiG(iparm)

!The shared linear parameter a_3
psiG_increm(iparm) = (psiG(iparm+1)+eps_FD)*exp(a_2)
psiG_increm(iparm+1) = psiG(iparm+1) + eps_FD

call calc_energy(psiG_increm,E_G_increm,E_FN1_increm)
gradient_VMC_FD(iparm+1)   = (E_G_increm-E_G)/eps_FD
gradient_FN_FD(iparm+1) = (E_FN1_increm-E_FN1)/eps_FD

psiG_psiFN1_increm(1:ndim)=eigenvectors(1:ndim,ndim)
if(psiG_psiFN1_increm(1).lt.0) psiG_psiFN1_increm(1:ndim)=-psiG_psiFN1_increm(1:ndim)
psiG_psiFN1_increm=psiG_psiFN1_increm/sum(abs(psiG_psiFN1_increm(1:ndim))) ! Since psiG_psiFN1_increm is already the product of 2 wavefns., normalize by the sum rather than the sqrt of the sum of squares
norm=sum(psiG_increm(1:ndim)**2)
psiFN1_by_psiG_increm(1:ndim) = psiG_psiFN1_increm(1:ndim)*norm/(psiG_increm(1:ndim)**2)

psiFN1_by_psiG_der(1:ndim)=(psiFN1_by_psiG_increm(1:ndim)-psiFN1_by_psiG(1:ndim))/eps_FD

do j=1,ndim
   gradient_FN1_exact(iparm+1)=gradient_FN1_exact(iparm+1)+psiG_psiFN(j)*((2*(psiG_der(j,iparm+1)/psiG(j))+psiFN1_by_psiG_der(j)/psiFN1_by_psiG(j))*(E_L(j)-E_FN1)+E_L_der(j,iparm+1))
enddo

!Calculate the finite-difference Hessian
psiG_increm    = psiG
hessian_FN_FD  = 0
hessian_VMC_FD = 0

!Do the linear terms
do i=1,ndim
   do j=1,ndim
      psiG_increm = psiG
      if((i.ne.iparm).and.(j.ne.iparm).and.(i.ne.iparm+1).and.(j.ne.iparm+1)) then
!Calculate the both increm part
         psiG_increm(i) = psiG(i) + eps_FD_2
         psiG_increm(j) = psiG(j) + eps_FD_2
         call calc_energy(psiG_increm,E_G_increm,E_FN1_increm)
         hessian_VMC_FD(i,j) = E_G_increm 
         hessian_FN_FD(i,j)  = E_FN1_increm 

!Calculate the i increm j decrem part
         psiG_increm(i) = psiG(i) + eps_FD_2
         psiG_increm(j) = psiG(j) - eps_FD_2
         call calc_energy(psiG_increm,E_G_increm,E_FN1_increm)
         hessian_VMC_FD(i,j) = hessian_VMC_FD(i,j) - E_G_increm 
         hessian_FN_FD(i,j)  = hessian_FN_FD(i,j) - E_FN1_increm 

!Calculate the i decrem j decrem part
         psiG_increm(i) = psiG(i) - eps_FD_2
         psiG_increm(j) = psiG(j) + eps_FD_2
         call calc_energy(psiG_increm,E_G_increm,E_FN1_increm)
         hessian_VMC_FD(i,j) = hessian_VMC_FD(i,j) - E_G_increm 
         hessian_FN_FD(i,j)  = hessian_FN_FD(i,j) - E_FN1_increm 
         
!Calculate the both decrem part
         psiG_increm(i) = psiG(i) - eps_FD_2
         psiG_increm(j) = psiG(j) - eps_FD_2
         call calc_energy(psiG_increm,E_G_increm,E_FN1_increm)
         hessian_VMC_FD(i,j) = (hessian_VMC_FD(i,j) + E_G_increm)/(4*eps_FD_2**2)
         hessian_FN_FD(i,j)  = (hessian_FN_FD(i,j) + E_FN1_increm)/(4*eps_FD_2**2)
      endif
   enddo
enddo

!Linear terms with same derivative and mixed wrt iparm+1
do i=1,ndim
   psiG_increm = psiG
   if((i.ne.iparm).and.(j.ne.iparm+1)) then
      !calculate derivatie wrt i,i
      !Calculate the forward increm part
      psiG_increm(i) = psiG(i)+eps_FD_2
      call calc_energy(psiG_increm,E_G_increm,E_FN1_increm)

      hessian_VMC_FD(i,i) = E_G_increm - 2*E_G
      hessian_FN_FD(i,i)  = E_FN1_increm - 2*E_FN1

      !Calculate the backward increm part
      psiG_increm(i) = psiG(i)-eps_FD_2
      call calc_energy(psiG_increm,E_G_increm,E_FN1_increm)

      hessian_VMC_FD(i,i) = (hessian_VMC_FD(i,i) + E_G_increm)/(eps_FD_2**2)
      hessian_FN_FD(i,i)  = (hessian_FN_FD(i,i) + E_FN1_increm)/(eps_FD_2**2)

      !Mixed terms wrt i,iparm
      !Calculate the both increm part
      psiG_increm    = psiG
      psiG_increm(i) = psiG(i) + eps_FD_2
      psiG_increm(iparm) = psiG(iparm+1)*exp(a_2+eps_FD_2)
      call calc_energy(psiG_increm,E_G_increm,E_FN1_increm)
      hessian_VMC_FD(i,iparm) = E_G_increm 
      hessian_FN_FD(i,iparm)  = E_FN1_increm 

      !Calculate the i increm j decrem part
      psiG_increm(i) = psiG(i) + eps_FD_2
      psiG_increm(iparm) = psiG(iparm+1)*exp(a_2-eps_FD_2)
      call calc_energy(psiG_increm,E_G_increm,E_FN1_increm)
      hessian_VMC_FD(i,iparm) = hessian_VMC_FD(i,iparm) - E_G_increm 
      hessian_FN_FD(i,iparm)  = hessian_FN_FD(i,iparm) - E_FN1_increm 

      !Calculate the i decrem j decrem part
      psiG_increm(i) = psiG(i) - eps_FD_2
      psiG_increm(iparm) = psiG(iparm+1)*exp(a_2+eps_FD_2)
      call calc_energy(psiG_increm,E_G_increm,E_FN1_increm)
      hessian_VMC_FD(i,iparm) = hessian_VMC_FD(i,iparm) - E_G_increm 
      hessian_FN_FD(i,iparm)  = hessian_FN_FD(i,iparm) - E_FN1_increm 

      !Calculate the both decrem part

      psiG_increm(i) = psiG(i) - eps_FD_2
      psiG_increm(iparm) = psiG(iparm+1)*exp(a_2-eps_FD_2)
      call calc_energy(psiG_increm,E_G_increm,E_FN1_increm)
      hessian_VMC_FD(i,iparm) = (hessian_VMC_FD(i,iparm) + E_G_increm)/(4*eps_FD_2**2)
      hessian_FN_FD(i,iparm)  = (hessian_FN_FD(i,iparm) + E_FN1_increm)/(4*eps_FD_2**2)
      hessian_VMC_FD(iparm,i) = hessian_VMC_FD(i,iparm)
      hessian_FN_FD(iparm,i)  = hessian_FN_FD(i,iparm)

      !Mixed wrt iparm+1
      
      !Calculate the both increm part
      psiG_increm    = psiG
      psiG_increm(i) = psiG(i) + eps_FD_2
      psiG_increm(iparm+1) = psiG(iparm+1)+eps_FD_2
      psiG_increm(iparm)   = (psiG(iparm+1)+eps_FD_2)*exp(a_2)
      call calc_energy(psiG_increm,E_G_increm,E_FN1_increm)
      hessian_VMC_FD(i,iparm+1) = E_G_increm 
      hessian_FN_FD(i,iparm+1)  = E_FN1_increm 

      !Calculate the i increm j decrem part
      psiG_increm(i) = psiG(i) + eps_FD_2
      psiG_increm(iparm+1) = psiG(iparm+1)-eps_FD_2
      psiG_increm(iparm)   = (psiG(iparm+1)-eps_FD_2)*exp(a_2)
      call calc_energy(psiG_increm,E_G_increm,E_FN1_increm)
      hessian_VMC_FD(i,iparm+1) = hessian_VMC_FD(i,iparm+1) - E_G_increm 
      hessian_FN_FD(i,iparm+1)  = hessian_FN_FD(i,iparm+1) - E_FN1_increm 

      !Calculate the i decrem j increm part
      psiG_increm(i) = psiG(i) - eps_FD_2
      psiG_increm(iparm+1) = psiG(iparm+1)+eps_FD_2
      psiG_increm(iparm)   = (psiG(iparm+1)+eps_FD_2)*exp(a_2)
      call calc_energy(psiG_increm,E_G_increm,E_FN1_increm)
      hessian_VMC_FD(i,iparm+1) = hessian_VMC_FD(i,iparm+1) - E_G_increm 
      hessian_FN_FD(i,iparm+1)  = hessian_FN_FD(i,iparm+1) - E_FN1_increm 

      !Calculate the both decerm part
      psiG_increm(i) = psiG(i) - eps_FD_2
      psiG_increm(iparm+1) = psiG(iparm+1)-eps_FD_2
      psiG_increm(iparm)   = (psiG(iparm+1)-eps_FD_2)*exp(a_2)
      call calc_energy(psiG_increm,E_G_increm,E_FN1_increm)
      hessian_VMC_FD(i,iparm+1) = (hessian_VMC_FD(i,iparm+1) + E_G_increm)/(4*eps_FD_2**2)
      hessian_FN_FD(i,iparm+1)  = (hessian_FN_FD(i,iparm+1) + E_FN1_increm)/(4*eps_FD_2**2)
      hessian_VMC_FD(iparm+1,i) = hessian_VMC_FD(i,iparm+1)
      hessian_FN_FD(iparm+1,i)  = hessian_FN_FD(i,iparm+1)
   endif
enddo


!2nd wrt iparm+1,iparm+1
!Calculate the forward increm part
psiG_increm = psiG
psiG_increm(iparm+1) = psiG(iparm+1)+eps_FD_2
psiG_increm(iparm)   = (psiG(iparm+1)+eps_FD_2)*exp(a_2)

call calc_energy(psiG_increm,E_G_increm,E_FN1_increm)

hessian_VMC_FD(iparm+1,iparm+1) = E_G_increm - 2*E_G
hessian_FN_FD(iparm+1,iparm+1)  = E_FN1_increm - 2*E_FN1

!Calculate the backward increm part
psiG_increm(iparm+1) = psiG(iparm+1)-eps_FD_2
psiG_increm(iparm)   = (psiG(iparm+1)-eps_FD_2)*exp(a_2)
call calc_energy(psiG_increm,E_G_increm,E_FN1_increm)

hessian_VMC_FD(iparm+1,iparm+1) = (hessian_VMC_FD(iparm+1,iparm+1) + E_G_increm)/(eps_FD_2**2)
hessian_FN_FD(iparm+1,iparm+1)  = (hessian_FN_FD(iparm+1,iparm+1) + E_FN1_increm)/(eps_FD_2**2)


!Calculate the mixed 2nd derivative wrt a_2 and a_3

psiG_increm = psiG
!Calculate both incremented
psiG_increm(iparm)   = (psiG(iparm+1)+eps_FD_2)*exp(a_2+eps_FD_2)
psiG_increm(iparm+1) = psiG(iparm+1)+eps_FD_2

call calc_energy(psiG_increm,E_G_increm,E_FN1_increm)

hessian_VMC_FD(iparm,iparm+1) = E_G_increm 
hessian_FN_FD(iparm,iparm+1)  = E_FN1_increm 

!Calculate a_2 increm, a_3 decrem
psiG_increm(iparm)   = (psiG(iparm+1)-eps_FD_2)*exp(a_2+eps_FD_2)
psiG_increm(iparm+1) = psiG(iparm+1)-eps_FD_2

call calc_energy(psiG_increm,E_G_increm,E_FN1_increm)

hessian_VMC_FD(iparm,iparm+1) = hessian_VMC_FD(iparm,iparm+1) - E_G_increm
hessian_FN_FD(iparm,iparm+1)  = hessian_FN_FD(iparm,iparm+1) - E_FN1_increm

!Calculate a_3 increm, a_2 decrem
psiG_increm(iparm)   = (psiG(iparm+1)+eps_FD_2)*exp(a_2-eps_FD_2)
psiG_increm(iparm+1) = psiG(iparm+1)+eps_FD_2

call calc_energy(psiG_increm,E_G_increm,E_FN1_increm)

hessian_VMC_FD(iparm,iparm+1) = hessian_VMC_FD(iparm,iparm+1) - E_G_increm
hessian_FN_FD(iparm,iparm+1)  = hessian_FN_FD(iparm,iparm+1) - E_FN1_increm

!Calculate both decrem
psiG_increm(iparm)   = (psiG(iparm+1)-eps_FD_2)*exp(a_2-eps_FD_2)
psiG_increm(iparm+1) = psiG(iparm+1)-eps_FD_2

call calc_energy(psiG_increm,E_G_increm,E_FN1_increm)

hessian_VMC_FD(iparm,iparm+1) = (hessian_VMC_FD(iparm,iparm+1) + E_G_increm)/(4*eps_FD_2**2)
hessian_FN_FD(iparm,iparm+1)  = (hessian_FN_FD(iparm,iparm+1) + E_FN1_increm)/(4*eps_FD_2**2)

hessian_FN_FD(iparm+1,iparm) = hessian_FN_FD(iparm,iparm+1)
hessian_VMC_FD(iparm+1,iparm) = hessian_VMC_FD(iparm,iparm+1)


!Calculate the second derivative wrt a_2^2

!Calculate the forward increm part
psiG_increm = psiG
psiG_increm(iparm) = psiG(iparm+1)*exp(a_2+eps_FD_2)
call calc_energy(psiG_increm,E_G_increm,E_FN1_increm)

hessian_VMC_FD(iparm,iparm) = E_G_increm - 2*E_G
hessian_FN_FD(iparm,iparm)  = E_FN1_increm - 2*E_FN1

!Calculate the backward increm part
psiG_increm(iparm) = psiG(iparm+1)*exp(a_2-eps_FD_2)

call calc_energy(psiG_increm,E_G_increm,E_FN1_increm)

hessian_VMC_FD(iparm,iparm) = (hessian_VMC_FD(iparm,iparm) + E_G_increm)/(eps_FD_2**2)
hessian_FN_FD(iparm,iparm)  = (hessian_FN_FD(iparm,iparm) + E_FN1_increm)/(eps_FD_2**2)


gradient_FN1_approx = 0
gradient_vmc        = 0
gradient_vmc_in_dmc = 0
gradient_FU         = 0
write(6,'(''W_der_by_W      ='',100f11.6)') W_der_by_W(1:ndim)
!Calculate the approximate gradient
do i=1,ndim
   do j=1,ndim
      !FN1_1 part
      gradient_FN1_approx(i) = gradient_FN1_approx(i) + psiG_psiFN(j)*((2*(psiG_der(j,i)/psiG(j)))*(E_L(j)-E_FN1))
      !FN1_3 part
      gradient_FN1_approx(i) = gradient_FN1_approx(i) + 2*psiG_psiFN(j)*E_L_der(j,i)
      gradient_FU(i) = gradient_FU(i) + psiG_psiFN(j)*((2*(psiG_der(j,i)/psiG(j)))*(E_L(j)-E_FN1) +&
           &E_L_der(j,i)) + W_der_by_W_full(j,i)*(E_L(j)-E_FN1)
      gradient_vmc(i) = gradient_vmc(i) + 2*psiG(j)**2*(psiG_der(j,i)/psiG(j)*(E_L(j)-E_G))
      gradient_vmc_in_dmc(i) = gradient_vmc_in_dmc(i) + 2*psiG_psiFN(j)*&
           &(psiG_der(j,i)/psiG(j)*(E_L(j)-E_G))
    enddo

enddo

hessian_vmc_in_dmc = 0
hessian_vmc = 0
hessian_FN1 = 0

!Calculate the hessian
do i=1,ndim
   do j=1,ndim
      do k=1,ndim

         hessian_vmc_in_dmc(i,j) = hessian_vmc_in_dmc(i,j)&
              &+psiG_psiFN(k)*(2*((psiG_2der(k,i,j)/psiG(k)+psiG_der(k,i)*psiG_der(k,j)&
              &/(psiG(k)**2))*(E_L(k)-E_G) - psiG_der(k,i)/psiG(k)*gradient_vmc_in_dmc(j)&
              &-psiG_der(k,j)/psiG(k)*gradient_vmc_in_dmc(i))&
              &+psiG_der(k,i)/psiG(k)*E_L_der(k,j)&
              &+psiG_der(k,j)/psiG(k)*E_L_der(k,i))


         hessian_vmc(i,j) = hessian_vmc(i,j)&
              &+psiG(k)**2*(2*((psiG_2der(k,i,j)/psiG(k)+psiG_der(k,i)*psiG_der(k,j)&
              &/(psiG(k)**2))*(E_L(k)-E_G) - psiG_der(k,i)/psiG(k)*gradient_vmc(j)&
              &-psiG_der(k,j)/psiG(k)*gradient_vmc(i))&
              &+psiG_der(k,i)/psiG(k)*E_L_der(k,j)&
              &+psiG_der(k,j)/psiG(k)*E_L_der(k,i))

         hessian_FN1(i,j) = hessian_FN1(i,j) + psiG_psiFN(k)*(&
              &2*((psiG_2der(k,i,j)/psiG(k)+psiG_der(k,i)*psiG_der(k,j)/psiG(k)**2)&
              &*(E_L(k)-E_FN1) - psiG_der(k,i)/psiG(k)*gradient_FN1_approx(j)&
              &- psiG_der(k,j)/psiG(k)*gradient_FN1_approx(i)) + &
              &2*(psiG_der(k,j)/psiG(k)*E_L_der(k,i)+psiG_der(k,i)/psiG(k)*E_L_der(k,j))&
              &+E_L_2der(k,i,j))&
              &+2*(psiG_der(k,i)/psiG(k)*W_der_by_W_full(k,j) + &
              &psiG_der(k,j)/psiG(k)*W_der_by_W_full(k,i))*(E_L(k)-E_FN1) + &
              &W_der_by_W_full(k,i)*E_L_der(k,j)+W_der_by_W_full(k,j)*E_L_der(k,i)&
              &- W_der_by_W_full(k,i)*gradient_FN1_approx(j) - &
              &W_der_by_W_full(k,j)*gradient_FN1_approx(i)
         
      enddo
      hessian_FN1(i,j) = hessian_FN1(i,j) + W_2der_by_W(i,j)

      !Calculate the error
      hessian_vmc_err(i,j) = (hessian_FN1(i,j)-hessian_fn_fd(i,j))!/hessian_fn_fd(i,j)

   enddo
enddo

write(6,'(/,''gradient_VMC_FD='',100f14.6)') (gradient_VMC_FD(i),i=1,ndim)
write(6,'(/,''gradient_VMC='',100f14.6)') (gradient_VMC(i),i=1,ndim)
write(6,'(/,''gradient_VMC_in_dmc='',100f14.6)') (gradient_VMC_in_dmc(i),i=1,ndim)

write(6,'(/,''gradient_FN1_FD='',100f14.6)') (gradient_FN_FD(i),i=1,ndim)
write(6,'(/,''gradient_FN1_exact='',100f14.6)') (gradient_FN1_exact(i),i=1,ndim)
write(6,'(/,''gradient_FN1_approx='',100f14.6)') (gradient_FN1_approx(i),i=1,ndim)
write(6,'(/,''gradient_FU='',100f14.6)') (gradient_FU(i),i=1,ndim)


write(6,*)'hessian_vmc_in_dmc'
do j=1,ndim
   write(6,'(100g14.6)') (hessian_vmc_in_dmc(i,j),i=1,ndim)
enddo

write(6,*)'hessian_VMC_FD'
do j=1,ndim
   write(6,'(100g14.6)') (hessian_VMC_FD(i,j),i=1,ndim)
enddo

write(6,*)'hessian_vmc'
do j=1,ndim
   write(6,'(100g14.6)') (hessian_vmc(i,j),i=1,ndim)
enddo

write(6,*)'hessian_FN_FD'
do j=1,ndim
   write(6,'(100g14.6)') (hessian_FN_FD(i,j),i=1,ndim)
enddo

write(6,*)'hessian_FN1'
do j=1,ndim
   write(6,'(100g14.6)') (hessian_FN1(i,j),i=1,ndim)
enddo

write(6,*)'hessian_FN_err'
do j=1,ndim
   write(6,'(100g14.6)') (hessian_vmc_err(i,j),i=1,ndim)
enddo



! call mat_vec_mul(ndim,MDIM,H,psiFN1,H_minus_E_psiFN1)
! H_minus_E_psiFN1(1:ndim)=H_minus_E_psiFN1(1:ndim)-E_FN1*psiFN1(1:ndim)
! call mat_vec_mul(ndim,MDIM,H,psiFN1_der,H_minus_E_psiFN1_der)
! H_minus_E_psiFN1_der(1:ndim)=H_minus_E_psiFN1_der(1:ndim)-E_FN1*psiFN1_der(1:ndim)
! call mat_vec_mul(ndim,MDIM,H,psiFN1_by_psiG,H_minus_E_psiFN1_by_psiG)
! H_minus_E_psiFN1_by_psiG(1:ndim)=H_minus_E_psiFN1_by_psiG(1:ndim)-E_FN1*psiFN1_by_psiG(1:ndim)
! call mat_vec_mul(ndim,MDIM,H,psiG_der,H_minus_E_psiG_der)
! H_minus_E_psiG_der(1:ndim)=H_minus_E_psiG_der(1:ndim)-E_FN1*psiG_der(1:ndim)
! call mat_vec_mul(ndim,MDIM,H,psi0-psiFN1,H_minus_E_psi0_min_psiFN1)
! H_minus_E_psi0_min_psiFN1(1:ndim)=H_minus_E_psi0_min_psiFN1(1:ndim)-E_FN1*(psi0(1:ndim)-psiFN1(1:ndim))
!test5a=dot_product(H_minus_E_psiFN1_der(1:ndim),psiG_der(1:ndim))/overlap_G_FN1
!test5b=-gradient_FN1*dot_product(psiFN1_der(1:ndim),psiG(1:ndim))/overlap_G_FN1
!test5c=-gradient_FN1*dot_product(psiFN1(1:ndim),psiG_der(1:ndim))/overlap_G_FN1
!test5d=dot_product(H_minus_E_psiG_der(1:ndim),psiG_der(1:ndim))
! test5e=dot_product(H_minus_E_psiFN1_by_psiG(1:ndim),psiG_der(1:ndim)**2)/overlap_G_FN1
! test5f=dot_product(H_minus_E_psiG_der(1:ndim),psiG_der(1:ndim)*psiFN1(1:ndim)/psiG(1:ndim))/overlap_G_FN1
! test5g=dot_product(H_minus_E_psiG_der(1:ndim),psiFN1_der(1:ndim))/overlap_G_FN1
! test7=dot_product(H_minus_E_psiG_der(1:ndim),psiFN1(1:ndim)-psi0(1:ndim))/overlap_G_FN1/overlap_G_FN1
! write(6,'(/,''gradient_FN1_1='',100f10.6)') (psiG_psiFN(i)*((2*(psiG_der(i)/psiG(i)))*(E_L(i)-E_FN1)),i=1,ndim)
! write(6,'(/,''gradient_FN1_2='',100f10.6)') (psiG_psiFN(i)*((psiFN1_by_psiG_der(i)/psiFN1_by_psiG(i))*(E_L(i)-E_FN1)),i=1,ndim)
! write(6,'(/,''gradient_FN1_3='',100f10.6)') (psiG_psiFN(i)*(E_L_der(i)),i=1,ndim)
! write(6,'(/,''E_FN1, E_FN1_increm, gradient_FN1_FD, gradient_FN1='',9f10.6)') E_FN1, E_FN1_increm, gradient_FN1_FD, gradient_FN1
! !write(6,'(/,''gradient_FN1_FD, gradient_FN1, 1st+2*2nd terms, gradient_FN1_1, gradient_FN1_2, gradient_FN1_3='',9f10.6)') gradient_FN1_FD, gradient_FN1, gradient_FN1_1+2*gradient_FN1_3, gradient_FN1_1, gradient_FN1_2, gradient_FN1_3
! write(6,'(/,''gradient_FN1_FD, gradient_FN1, 1st+2*2nd terms, gradient_FU2000, gradient_FN1_1, gradient_FN1_2, gradient_FN1_3, W_der_by_W='',9f10.6)') gradient_FN1_FD, gradient_FN1, gradient_FN1_1+2*gradient_FN1_3, gradient_FN1_1+W_der_by_W+gradient_FN1_3, gradient_FN1_1, gradient_FN1_2, gradient_FN1_3, W_der_by_W
! !write(6,'(''W_der_by_W='',t136,9f10.6)') W_der_by_W
! !write(6,'(/,i2,'' gradient_vmc, gradient_vmc_eval_in_dmc, gradient_FN1, gradient_FN1_1+2*gradient_FN1_3, gradient_FN1_1+W_der_by_W+gradient_FN1_3='',9f10.6)') itimes, gradient_vmc, gradient_vmc_eval_in_dmc, gradient_FN1, gradient_FN1_1+2*gradient_FN1_3, gradient_FN1_1+W_der_by_W+gradient_FN1_3
! !write(6,'(/,2i3,'' gradient_vmc, gradient_vmc_eval_in_dmc, gradient_FN1, gradient_FN1_1+2*gradient_FN1_3, gradient_FN1_1+W_der_by_W+gradient_FN1_3='',9f10.6)') itimes, ipsiG, gradient_vmc, gradient_vmc_eval_in_dmc, gradient_FN1, gradient_FN1_1+2*gradient_FN1_3, gradient_FN1_1+W_der_by_W+gradient_FN1_3
! write(6,'(/,2i3,'' gradient_vmc, gradient_vmc_eval_in_dmc, gradient_FN1, gradient_FN1_1+2*gradient_FN1_3, gradient_FN1_1+W_der_by_W+gradient_FN1_3='',5f10.6,'' <psi_FNi(H-e)psi>, <psi_i(H-e)psiFN>='', 10f10.6)') itimes, ipsiG, gradient_vmc, gradient_vmc_eval_in_dmc, gradient_FN1, gradient_FN1_1+2*gradient_FN1_3, gradient_FN1_1+W_der_by_W+gradient_FN1_3, dot_product(H_minus_E_psiFN1_der(1:ndim),psiG(1:ndim))/dot_product(psiG(1:ndim),psiFN1(1:ndim)), dot_product(H_minus_E_psiG_der(1:ndim),psiFN1(1:ndim))/dot_product(psiG(1:ndim),psiFN1(1:ndim)), test7, dot_product(psi0(1:ndim)-psiFN1(1:ndim),psiG_der(1:ndim))
! write(6,'(/,i2,'' Ratio: Error replace 2nd by 3rd term  to  Error replace 2nd by W_der_by_W='',9f10.6)') itimes, (gradient_FN1_3-gradient_FN1_2)/(W_der_by_W-gradient_FN1_2)
! if(abs(gradient_FN1_FD-gradient_FN1).gt.1.e-5_rk) then
!   write(6,'(''gradient_FN1_FD, gradient_FN1, gradient_FN1_FD-gradient_FN1='',2f10.6,es9.1)') gradient_FN1_FD, gradient_FN1, gradient_FN1_FD-gradient_FN1
!   stop 'gradient_FN1_FD .ne. gradient_FN1'
! endif
write(6,'(''FN1*G='',99f14.10)') (psiFN1(i)*psiG(i),i=1,ndim),sum(psiFN1(1:ndim)*psiG(1:ndim))
write(6,'(''FN1_G='',99f14.10)') (psiG_psiFN(i),i=1,ndim),sum(psiG_psiFN(1:ndim))
! write(6,'(''test='',9f10.6)') sum(psiFN1(1:ndim)*psiG_der(1:ndim)*E_L_der(1:ndim)), &
!                              sum(psiG_psiFN(1:ndim)*psiG_der(1:ndim)*E_L_der(1:ndim)/psiG(1:ndim))
! write(6,'(''test2='',6f10.6,9f16.12)') &
!    vec_mat_vec_mul(ndim,MDIM,psiFN1_der,H,psiG), vec_mat_vec_mul(ndim,MDIM,psiG_der,H,psiFN1), &
!    vec_mat_vec_mul(ndim,MDIM,psiFN1,H,psiG_der), vec_mat_vec_mul(ndim,MDIM,psiG,H,psiFN1_der), &
!    E_FN1*dot_product(psiFN1(1:ndim),psiG_der(1:ndim)), E_FN1*dot_product(psiFN1_der(1:ndim),psiG(1:ndim)), &
!    vec_mat_vec_mul(ndim,MDIM,psiFN1_der,H,psiG) - vec_mat_vec_mul(ndim,MDIM,psiG_der,H,psiFN1) + &
!    E_FN1*(dot_product(psiFN1(1:ndim),psiG_der(1:ndim))-dot_product(psiFN1_der(1:ndim),psiG(1:ndim))), &
!    vec_mat_vec_mul(ndim,MDIM,psiFN1_der,H,psiG) - vec_mat_vec_mul(ndim,MDIM,psiG_der,H,psiFN1) + &
!    E_0*(dot_product(psiFN1(1:ndim),psiG_der(1:ndim))-dot_product(psiFN1_der(1:ndim),psiG(1:ndim))), E_FN1-E_0
! write(6,'(''test3='',6f10.6)') dot_product(psiFN1(1:ndim),psiFN1(1:ndim)), dot_product(psiFN1_der(1:ndim),psiFN1_der(1:ndim)), &
!                                dot_product(psiG(1:ndim),psiG(1:ndim)), dot_product(psiG_der(1:ndim),psiG_der(1:ndim)), &
!                                dot_product(psiFN1(1:ndim),psiG_der(1:ndim)), dot_product(psiFN1_der(1:ndim),psiG(1:ndim))
! write(6,'(''test4='',6f10.6)') dot_product(H_minus_E_psiFN1(1:ndim),H_minus_E_psiFN1(1:ndim)), dot_product(H_minus_E_psiFN1(1:ndim),psiG_der(1:ndim))
! write(6,'(''test5='',10f10.6)') 2*(test5a+test5b+test5c),test5a,test5b,test5c,test5d,test5e,test5f,test5g
! write(6,'(''test6a='',10f10.6)') psiG_der(1:ndim)
! write(6,'(''test6b='',10f10.6)') psiG_der(1:ndim)*psiFN1(1:ndim)/psiG(1:ndim)
! write(6,'(''test6c='',10f10.6)') psiFN1_der(1:ndim)
! write(6,'(/,2i3,'' <psi_FNi(H-e)psi>, <psi_i(H-e)psiFN>, gradient_FN1 ='',10f10.6)') itimes, ipsiG, dot_product(H_minus_E_psiFN1_der(1:ndim),psiG(1:ndim))/dot_product(psiG(1:ndim),psiFN1(1:ndim)), dot_product(H_minus_E_psiG_der(1:ndim),psiFN1(1:ndim))/dot_product(psiG(1:ndim),psiFN1(1:ndim)), gradient_FN1, test7, dot_product(psi0(1:ndim)-psiFN1(1:ndim),psiG_der(1:ndim)), dot_product(H_minus_E_psi0_min_psiFN1(1:ndim),H_minus_E_psi0_min_psiFN1(1:ndim)), dot_product(H_minus_E_psiG_der(1:ndim),H_minus_E_psiG_der(1:ndim))
! write(6,'(''test7='',10f10.6)') test7,dot_product(H_minus_E_psiG_der(1:ndim),H_minus_E_psiG_der(1:ndim))
! write(6,'(''test8='',10f10.6)') gradient_FN1, 2*dot_product(H_minus_E_psiG_der(1:ndim),psiFN1(1:ndim))/overlap_G_FN1, dot_product(H_minus_E_psiG_der(1:ndim),psiFN1(1:ndim)**2/psiG(1:ndim)), &
!  dot_product(H_minus_E_psiG_der(1:ndim),psiG(1:ndim)), &
!  dot_product(H_minus_E_psiG_der(1:ndim),psiFN1(1:ndim)**2/psiG(1:ndim))+dot_product(H_minus_E_psiG_der(1:ndim),psiG(1:ndim)), &
!  2*dot_product(H_minus_E_psiG_der(1:ndim),psiFN1(1:ndim))/overlap_G_FN1 - dot_product(H_minus_E_psiG_der(1:ndim),psiFN1(1:ndim)**2/psiG(1:ndim)) - dot_product(H_minus_E_psiG_der(1:ndim),psiG(1:ndim)), &
!  2*dot_product(H_minus_E_psiG_der(1:ndim),psiFN1(1:ndim))/overlap_G_FN1-gradient_FN1, &
!  dot_product(H_minus_E_psiG_der(1:ndim),psiFN1(1:ndim)**2/psiG(1:ndim))+dot_product(H_minus_E_psiG_der(1:ndim),psiG(1:ndim))-gradient_FN1
! !write(6,'(''test999='',10f10.6)') dot_product(H_minus_E_psiG_der(1:ndim),psiG(1:ndim)), dot_product(H_minus_E_psiG(1:ndim),psiG_der(1:ndim))
! ! Test again that <E_Li> is 0 for VMC distribution
! write(6,'(''test99='',10f10.6)') dot_product(psiG_psiFN(1:ndim),overlap_G_FN1*E_L_der(1:ndim)/psiFN1_by_psiG(1:ndim))
! ! The error in the mean of our approximate DMC 
! write(6,'(''test9='',10f10.6)') 2*error_of_mean(ndim,psiG_psiFN(1:ndim),H_minus_E_psiG_der(1:ndim)/psiG(1:ndim)), &
!   2*error_of_mean(ndim,psiG_psiFN(1:ndim),H_minus_E_psiG_der(1:ndim)/psiG(1:ndim)-overlap_G_FN1*E_L_der(1:ndim)/psiFN1_by_psiG(1:ndim)), &
!   error_of_mean(ndim,psiG_psiFN(1:ndim),overlap_G_FN1*psiFN1_by_psiG(1:ndim)*H_minus_E_psiG_der(1:ndim)/psiG(1:ndim)+(overlap_G_FN1/psiFN1_by_psiG(1:ndim))*(E_L(1:ndim)-E_FN1)*psiG_der(1:ndim)/psiG(1:ndim))


write(6,'(''------------------------------------------------------------'')')

write(2,'(f8.4,es12.2,i4,15f16.6)') psiG_noise, 1-overlap_FN1, ipsiG, gradient_FN1_exact(1), gradient_FU(1), gradient_FN1_approx(1), gradient_vmc(1), gradient_vmc_in_dmc(1), gradient_FN1_exact(iparm), gradient_FU(iparm), gradient_FN1_approx(iparm), gradient_vmc(iparm), gradient_vmc_in_dmc(iparm), gradient_FN1_exact(iparm+1), gradient_FU(iparm+1), gradient_FN1_approx(iparm+1), gradient_vmc(iparm+1), gradient_vmc_in_dmc(iparm+1)

write(3,'(f8.4,es12.2,i4,15f16.6)') psiG_noise, 1-overlap_FN1, ipsiG, hessian_fn_fd(1,1), hessian_fn1(1,1), hessian_vmc(1,1), hessian_vmc_in_dmc(1,1), hessian_fn_fd(iparm,iparm), hessian_fn1(iparm,iparm), hessian_vmc(iparm,iparm), hessian_vmc_in_dmc(iparm,iparm), hessian_fn_fd(iparm,iparm+1), hessian_fn1(iparm,iparm+1), hessian_vmc(iparm,iparm+1), hessian_vmc_in_dmc(iparm,iparm+1)


enddo ! ipsiG

enddo ! itimes

end subroutine read_input
!-----------------------------------------------------------------------

subroutine calc_energy(psiG_input,E_G_out,E_FN_out)
  real(rk), intent(in)  :: psiG_input(MDIM)
  real(rk), intent(out) :: E_G_out,E_FN_out
  integer               :: i_loc,j_loc

  ! Calculate the variational energy for input psiG
  E_G_out=0
  do j_loc=1,ndim
     do i_loc=1,ndim
        E_G_out=E_G_out+psiG_input(j_loc)*H(j_loc,i_loc)*psiG_input(i_loc)
     enddo
  enddo
  E_G_out=E_G_out/dot_product(psiG_input(1:ndim),psiG_input(1:ndim))

  ! Calculate finite-difference E_{L,iparm}(i), i.e., derivative of the local energy on state i wrt the iparm parameter (just for check, can be removed)
  ! Define incremented importance-sampled projector
  do j_loc=1,ndim
     do i_loc=1,ndim
        H_imp(j_loc,i_loc)=psiG_input(j_loc)*H(j_loc,i_loc)/psiG_input(i_loc)
        P_imp(j_loc,i_loc)=psiG_input(j_loc)*P(j_loc,i_loc)/psiG_input(i_loc)
     enddo
  enddo

! write (fmt,'(a,i2,a)') '(/,''P_imp='',/,(', ndim, 'f9.4))'
! write(6,fmt) ((P_imp(j_loc,i_loc),i_loc=1,ndim),j_loc=1,ndim)

  ! Make incremented FN1 projector
  do i_loc=1,ndim
     sign_flip=0
     do j_loc=1,ndim
        if(j_loc.ne.i_loc) then
           if(P_imp(j_loc,i_loc).ge.0) then
              P_FN1(j_loc,i_loc)=P_imp(j_loc,i_loc)
           else
              sign_flip=sign_flip+P_imp(j_loc,i_loc)
              P_FN1(j_loc,i_loc)=0
           endif
        endif
     enddo
     !P_FN1(i,i)=P_imp(i,i)+sign_flip+tau(E_FN1-E_0)
     P_FN1(i_loc,i_loc)=P_imp(i_loc,i_loc)+sign_flip                 ! This is correct because in calculating eigenval_h below we use E_T=E_0 rather than E_FN1
     if(P_FN1(i_loc,i_loc).lt.0) write(6,'(''Warning: for itimes, i='',2i4,'', P_FN1(i,i)='',es12.4)') itimes, i_loc, P_FN1(i_loc,i_loc)
  enddo
  ! write (fmt,'(a,i2,a)') '(/,''P_FN1='',/,(', ndim, 'f9.4))'
  ! write(6,fmt) ((P_FN1(j_loc,i_loc),i_loc=1,ndim),j_loc=1,ndim)

! Diagonalize to calculate eigevector and eigenvalues for P_FN1 and deduce eigenvalue for H_FN1.
! Eigenvectors is the matrix on input and the eigenvectors on output
  eigenvectors=P_FN1
  call real_general_diagonalize_ow_ham(ndim,MDIM,eigenvectors,eigenvalues,eigenvalues_imag)
  eigenval_p=eigenvalues(ndim)
  eigenval_h=(1-eigenval_p+tau*E_0)/tau
  E_FN_out=eigenval_h

end subroutine calc_energy

!-----------------------------------------------------------------------
subroutine vec_vec_elem_mul(ndim,vec1_in,vec2_in,vec_out)
! Elements of output vector are products of the corresponding elements of the 2 input vectors
integer, intent(in) :: ndim
real(rk), intent(in) :: vec1_in(MDIM), vec2_in(MDIM)
real(rk), intent(out) :: vec_out(MDIM)

  vec_out(1:ndim)=vec1_in(1:ndim)*vec2_in(1:ndim)

end subroutine vec_vec_elem_mul
!-----------------------------------------------------------------------
subroutine mat_vec_mul(ndim,MDIM,mat,vec_in,vec_out)
! Matrix-vector multiply
integer, intent(in) :: ndim, MDIM
real(rk), intent(in) :: mat(MDIM,MDIM), vec_in(MDIM)
real(rk), intent(out) :: vec_out(MDIM)
integer i,j

vec_out(1:ndim)=0
do i=1,ndim
  do j=1,ndim
    vec_out(i)=vec_out(i)+mat(i,j)*vec_in(j)
  enddo
enddo
end subroutine mat_vec_mul
!-----------------------------------------------------------------------
function vec_mat_vec_mul(ndim,MDIM,vec1,mat,vec2)
! vector-matrix-vector multiply
integer, intent(in) :: ndim, MDIM
real(rk), intent(in) :: vec1(MDIM), mat(MDIM,MDIM), vec2(MDIM)
real(rk) vec_tmp(ndim)
real(rk) vec_mat_vec_mul

call mat_vec_mul(ndim,MDIM,mat,vec2,vec_tmp)
vec_mat_vec_mul=dot_product(vec1(1:ndim),vec_tmp(1:ndim))

end function vec_mat_vec_mul
!-----------------------------------------------------------------------
subroutine mat_mul_rt_diag(ndim,MDIM,mat_in,vec,mat_out)
! Matrix-diagonal_matrix multiply
! The matrix to the right is diagonal and is therefore inputted as a vector
integer, intent(in) :: ndim, MDIM
real(rk), intent(in) :: mat_in(MDIM,MDIM), vec(MDIM)
real(rk), intent(out) :: mat_out(MDIM,MDIM)
integer i,j

do i=1,ndim
  do j=1,ndim
    mat_out(i,j)=mat_in(i,j)*vec(j)
  enddo
enddo
end subroutine mat_mul_rt_diag
!-----------------------------------------------------------------------
subroutine mat_mul_lt_diag(ndim,MDIM,vec,mat_in,mat_out)
! Diagonal_matrix-Matrix multiply
! The matrix to the left is diagonal and is therefore inputted as a vector
integer, intent(in) :: ndim, MDIM
real(rk), intent(in) :: mat_in(MDIM,MDIM), vec(MDIM)
real(rk), intent(out) :: mat_out(MDIM,MDIM)
integer i,j

do i=1,ndim
  do j=1,ndim
    mat_out(i,j)=vec(i)*mat_in(i,j)
  enddo
enddo
end subroutine mat_mul_lt_diag
!-----------------------------------------------------------------------

subroutine mat_mul(ndim,MDIM,mat1,mat2)
! Calculates mat2 <- mat1 * mat2
integer, intent(in) :: ndim, MDIM
real(rk), intent(in) :: mat1(MDIM,MDIM)
real(rk), intent(inout) :: mat2(MDIM,MDIM)
real(rk) :: mattmp(ndim,ndim)
integer i,j,k

mattmp(1:ndim,1:ndim)=0
do i=1,ndim
  do j=1,ndim
    do k=1,ndim
      mattmp(i,j)=mattmp(i,j)+mat1(i,k)*mat2(k,j)
    enddo
  enddo
enddo
mat2(1:ndim,1:ndim)=mattmp(1:ndim,1:ndim)
end subroutine mat_mul
!-----------------------------------------------------------------------

subroutine mat_mul2(ndim,MDIM,mat1,mat2)
! Calculates mat2 <- 1 + mat1*mat2
integer, intent(in) :: ndim, MDIM
real(rk), intent(in) :: mat1(MDIM,MDIM)
real(rk), intent(inout) :: mat2(MDIM,MDIM)
real(rk) :: mattmp(ndim,ndim)
integer i,j,k

mattmp(1:ndim,1:ndim)=0
do i=1,ndim
  do j=1,ndim
    do k=1,ndim
      mattmp(i,j)=mattmp(i,j)+mat1(i,k)*mat2(k,j)
    enddo
  enddo
  mattmp(i,i)=1+mattmp(i,i)
enddo
mat2(1:ndim,1:ndim)=mattmp(1:ndim,1:ndim)
end subroutine mat_mul2
!-----------------------------------------------------------------------

function rms_deviation(ndim,E_L,rho)
! Calculates rms_deviation of E_L wrt distribution rho
integer, intent(in) :: ndim
real(rk), intent(in) :: E_L(ndim), rho(MDIM)
real(rk) rho_av, E_L_av, E_L2_av, rms_deviation
integer i

rho_av=0; E_L_av=0; E_L2_av=0
do i=1,ndim
  rho_av=rho_av+rho(i)
  E_L_av=E_L_av+rho(i)*E_L(i)
  E_L2_av=E_L2_av+rho(i)*E_L(i)**2
enddo
E_L_av=E_L_av/rho_av
E_L2_av=E_L2_av/rho_av
rms_deviation=sqrt(E_L2_av-E_L_av**2)

end function rms_deviation
!-----------------------------------------------------------------------

function error_of_mean(ndim,w,x)
! Calculates the error of the weighted population mean
! n_eff = sum_w**2/sum_w2
! 1) The standard deviation of the population is
! sqrt(sum(w(1:ndim)*(x(1:ndim)-x_av)**2)/sum_w)
! 2) The unbiased estimate of the standard deviation from a sample is
! sqrt(sum(w(1:ndim)*(x(1:ndim)-x_av)**2)/(sum_w*(1-sum_w2/sum_w**2)))
! To calculate the error of the weighted population mean, divide the standard
! deviation by n_eff, so,
! 3) The error of the weighted population mean is
! sqrt(sum(w(1:ndim)*(x(1:ndim)-x_av)**2)/(sum_w**3/sum_w2))
! 4) The unbiased estimate of the error of the weighted mean from a sample is
! sqrt(sum(w(1:ndim)*(x(1:ndim)-x_av)**2)/(sum_w*((sum_w**2/sum_w2)-1)))
! Here we are computing 3)
integer, intent(in) :: ndim
real(rk), intent(in) :: w(ndim), x(ndim)
real(rk) error_of_mean, sum_w, sum_w2, x_av

sum_w=sum(w(1:ndim))
sum_w2=sum(w(1:ndim)**2)
x_av=sum(w(1:ndim)*x(1:ndim))/sum_w

!error_of_mean=sqrt(sum(w(1:ndim)*(x(1:ndim)-x_av)**2)/(sum_w**2-sum_w2))
!error_of_mean=sqrt(sum(w(1:ndim)*(x(1:ndim)-x_av)**2)/sum_w)
error_of_mean=sqrt(sum(w(1:ndim)*(x(1:ndim)-x_av)**2)/(sum_w**3/sum_w2))

end function error_of_mean
!-----------------------------------------------------------------------

function overlap(ndim,a,b)
integer, intent(in) :: ndim
real(rk), intent(in) :: a(ndim), b(ndim)
real(rk) overlap

overlap=sum(a(:)*b(:))/sqrt(sum(a(:)*a(:))*sum(b(:)*b(:)))

end function overlap
!-----------------------------------------------------------------------

end module fn
!-----------------------------------------------------------------------

program fixed_node
use fn, only: read_input

call read_input

end program fixed_node
