module imp_sam

use types, only: rk
!use more_tools, only: real_symmetric_diagonalize, real_general_diagonalize_ow_ham
implicit none

integer, parameter :: MPTS = 10**7 ! 10000000
!integer, dimension(4) :: irand_seed
integer npts, i, ifunc
real(rk) :: a, b, xmax, eps, e, e_num, dpsi_da_h_psi_av_by_dpsi_da_psi_av, dpsi_db_h_psi_av_by_dpsi_db_psi_av, dpsi_da_h_psi_av_by_psisq_av, dpsi_db_h_psi_av_by_psisq_av, dpsi_db_h_minus_e_psi_av_by_psisq_av, rms, &
 x(MPTS), psi(MPTS), grad_psi(MPTS), h_psi(MPTS), h_minus_e_psi(MPTS), E_L(MPTS), dpsi_da(MPTS), dpsi_db(MPTS), one_vec(MPTS), imp_fn(MPTS)

real(rk) :: rannyu

!external rannyu

contains

subroutine read_input

real(rk) tmp
! npts is # of integration pts.  The more the points the closer they get to the singularity
! xmax can be 1 or 2.
! If xmax=2, npts should be even and a=0 for points to get closer to node at 1 as npts is increased
! If xmax=1 and b=1, the singularity is at the edge of the interval at 1
! If xmax=2 and b=1, and a=1/8, the singularity is at 1/
read(5,*) npts, a, b, eps
write(6,'(''npts, a, b, eps='',i8,3f10.6)') npts, a, b, eps
if(npts.gt.MPTS) stop 'npts>MPTS'
read(5,*) ifunc, xmax
write(6,'(''ifunc, xmax='',i2,f8.4)') ifunc, xmax
if(ifunc.lt.1 .or. ifunc.gt.3) stop 'ifunc should be 1 or 2 or 3'
if((ifunc.eq.1 .or. ifunc.eq.2) .and. abs(xmax-1._rk).gt.1.e-10) stop 'xmax should be 1 for ifunc = 1 or 2'
if(ifunc.eq.3 .and. abs(xmax-2._rk).gt.1.e-10) stop 'xmax should be 2 for ifunc = 3'
if(abs(b-1._rk).gt.1.e-10) stop 'b should be 1 so that psi goes to 0 at 1 for ifunc=1,2 and at 2 for ifunc=3'
if(ifunc.eq.3 .and. abs(a-0._rk).gt.1.e-4) stop 'a should be 0 so that psi goes to 0 at 1 (midpoint of interval)'

! Check that derivatives are correctly programmed
a=a+1.e-6_rk
tmp=psi_func(.3_rk)
a=a-1.e-6_rk
write(6,'(''dpsi_da'',9f10.6)') (tmp-psi_func(.3_rk))/1.e-6_rk, dpsi_da_func(.3_rk)
b=b+1.e-6_rk
tmp=psi_func(.3_rk)
b=b-1.e-6_rk
write(6,'(''dpsi_db'',9f10.6)') (tmp-psi_func(.3_rk))/1.e-6_rk, dpsi_db_func(.3_rk)

! Place points to avoid end pts.
do i=1,npts
  x(i)=(i-.5_rk)/npts
  psi(i)=psi_func(x(i))
  grad_psi(i)=grad_psi_func(x(i))
  h_psi(i)=h_psi_func(x(i))
  E_L(i)=h_psi_func(x(i))/psi_func(x(i))
  dpsi_da(i)=dpsi_da_func(x(i))
  dpsi_db(i)=dpsi_db_func(x(i))
  one_vec(i)=1
enddo
e=dot_product(psi(1:npts),h_psi(1:npts))/dot_product(psi(1:npts),psi(1:npts))
e_num=dot_product(psi(1:npts),h_psi(1:npts))
write(6,'(''e_num, one_vec_int='',9f10.6)') e_num, real(npts)
h_minus_e_psi(1:npts)=h_psi(1:npts)-e*psi(1:npts)
dpsi_da_h_psi_av_by_dpsi_da_psi_av=dot_product(dpsi_da(1:npts),h_psi(1:npts))/dot_product(dpsi_da(1:npts),psi(1:npts))
dpsi_db_h_psi_av_by_dpsi_db_psi_av=dot_product(dpsi_db(1:npts),h_psi(1:npts))/dot_product(dpsi_db(1:npts),psi(1:npts))
dpsi_da_h_psi_av_by_psisq_av=dot_product(dpsi_da(1:npts),h_psi(1:npts))/dot_product(psi(1:npts),psi(1:npts))
dpsi_db_h_psi_av_by_psisq_av=dot_product(dpsi_db(1:npts),h_psi(1:npts))/dot_product(psi(1:npts),psi(1:npts))
dpsi_db_h_minus_e_psi_av_by_psisq_av=dot_product(dpsi_db(1:npts),h_minus_e_psi(1:npts))/dot_product(psi(1:npts),psi(1:npts))
write(6,'(''Variational energy E, parameter deriv E_b='',9f10.6)') e, dpsi_db_h_minus_e_psi_av_by_psisq_av

!write(6,'(''x      ='',1000f10.5)') x(1:npts)
!write(6,'(''psi    ='',1000f10.5)') psi(1:npts)
!write(6,'(''h_psi  ='',1000f10.5)') h_psi(1:npts)
if(npts.le.100) write(6,'(''psi    ='',1000f14.5)') psi(1:npts)
if(npts.le.100) write(6,'(''h_psi  ='',1000f14.5)') h_psi(1:npts)
if(npts.le.100) write(6,'(''E_L    ='',1000f14.5)') E_L(1:npts)
if(npts.le.100) write(6,'(''dpsi_da='',1000f14.5)') dpsi_da(1:npts)
if(npts.le.100) write(6,'(''dpsi_db='',1000f14.5)') dpsi_db(1:npts)

!a=a+.0001
!do i=1,npts
!  dpsi_da(i)=psi_func(x(i))
!  dpsi_da(i)=(dpsi_da(i)-psi(i))/.0001
!enddo
!write(6,'(''dpsi_da='',1000f10.5)') dpsi_da(1:npts)

write(6,'(''1: <psi H psi> / <psi^2>'',f10.6)') e
imp_fn(1:npts)=1
rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)**2/imp_fn(1:npts),psi(1:npts)**2/imp_fn(1:npts))
write(6,'(''imp_fn=1,                rms='',f18.5)') rms
imp_fn(1:npts)=psi(1:npts)**.5
rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)**2/imp_fn(1:npts),psi(1:npts)**2/imp_fn(1:npts))
write(6,'(''imp_fn=psi**.5,          rms='',f18.5)') rms
imp_fn(1:npts)=psi(1:npts)
rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)**2/imp_fn(1:npts),psi(1:npts)**2/imp_fn(1:npts))
write(6,'(''imp_fn=psi,              rms='',f18.5)') rms
imp_fn(1:npts)=psi(1:npts)**2
rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)**2/imp_fn(1:npts),psi(1:npts)**2/imp_fn(1:npts))
write(6,'(''imp_fn=psi**2,           rms='',f18.5)') rms
!imp_fn(1:npts)=psi(1:npts)*(1+psi(1:npts))
!rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)**2/imp_fn(1:npts),psi(1:npts)**2/imp_fn(1:npts))
!write(6,'(''imp_fn=psi+psi^2,        rms='',f18.5)') rms
imp_fn(1:npts)=1+psi(1:npts)
rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)**2/imp_fn(1:npts),psi(1:npts)**2/imp_fn(1:npts))
write(6,'(''imp_fn=1+psi,            rms='',f18.5)') rms
imp_fn(1:npts)=1+psi(1:npts)**2
rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)**2/imp_fn(1:npts),psi(1:npts)**2/imp_fn(1:npts))
write(6,'(''imp_fn=1+psi^2           rms='',f18.5)') rms
!imp_fn(1:npts)=.01+psi(1:npts)**.5
!rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)**2/imp_fn(1:npts),psi(1:npts)**2/imp_fn(1:npts))
!write(6,'(''imp_fn=.01+psi^.5        rms='',f18.5)') rms
!imp_fn(1:npts)=psi(1:npts)**.5 + .01*psi(1:npts)
!rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)**2/imp_fn(1:npts),psi(1:npts)**2/imp_fn(1:npts))
!write(6,'(''imp_fn=psi^.5+.01psi,    rms='',f18.5)') rms
!imp_fn(1:npts)=abs(psi(1:npts)*h_psi(1:npts)-e*psi(1:npts)**2)+.01*psi(1:npts)**.5
!rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)**2/imp_fn(1:npts),psi(1:npts)**2/imp_fn(1:npts))
!write(6,'(''imp_fn=Optimal+.01psi^.5 rms='',f18.5)') rms
call imp_fn_sub1(psi,grad_psi,imp_fn)
rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)**2/imp_fn(1:npts),psi(1:npts)**2/imp_fn(1:npts))
write(6,'(''imp_fn=imp_fn1           rms='',f18.5)') rms
call imp_fn_sub2(psi,grad_psi,imp_fn)
rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)**2/imp_fn(1:npts),psi(1:npts)**2/imp_fn(1:npts))
write(6,'(''imp_fn=imp_fn2           rms='',f18.5)') rms
imp_fn(1:npts)=abs(psi(1:npts)*h_psi(1:npts)-e*psi(1:npts)**2)
rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)**2/imp_fn(1:npts),psi(1:npts)**2/imp_fn(1:npts))
write(6,'(''imp_fn=Optimal -----     rms='',f18.5)') rms
!imp_fn(1:npts)=abs(psi(1:npts)**1.1*h_psi(1:npts)-e*psi(1:npts)**2.1)
!rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)**2/imp_fn(1:npts),psi(1:npts)**2/imp_fn(1:npts))
!write(6,'(''imp_fn=?,                rms='',f18.5)') rms
write(6,*)

write(6,'(''2: <psi_a H psi> / <psi_a psi>'',f10.6)') dpsi_da_h_psi_av_by_dpsi_da_psi_av
imp_fn(1:npts)=1
rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)*dpsi_da(1:npts)/imp_fn(1:npts),psi(1:npts)*dpsi_da(1:npts)/imp_fn(1:npts))
write(6,'(''imp_fn=1,                rms='',f18.5)') rms
imp_fn(1:npts)=psi(1:npts)**.5
rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)*dpsi_da(1:npts)/imp_fn(1:npts),psi(1:npts)*dpsi_da(1:npts)/imp_fn(1:npts))
write(6,'(''imp_fn=psi**.5,          rms='',f18.5)') rms
imp_fn(1:npts)=psi(1:npts)
rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)*dpsi_da(1:npts)/imp_fn(1:npts),psi(1:npts)*dpsi_da(1:npts)/imp_fn(1:npts))
write(6,'(''imp_fn=psi,              rms='',f18.5)') rms
imp_fn(1:npts)=psi(1:npts)**2
rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)*dpsi_da(1:npts)/imp_fn(1:npts),psi(1:npts)*dpsi_da(1:npts)/imp_fn(1:npts))
write(6,'(''imp_fn=psi**2,           rms='',f18.5)') rms
!imp_fn(1:npts)=psi(1:npts)*(1+psi(1:npts))
!rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)*dpsi_da(1:npts)/imp_fn(1:npts),psi(1:npts)*dpsi_da(1:npts)/imp_fn(1:npts))
!write(6,'(''imp_fn=psi+psi^2,        rms='',f18.5)') rms
imp_fn(1:npts)=1+psi(1:npts)
rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)*dpsi_da(1:npts)/imp_fn(1:npts),psi(1:npts)*dpsi_da(1:npts)/imp_fn(1:npts))
write(6,'(''imp_fn=1+psi,            rms='',f18.5)') rms
imp_fn(1:npts)=1+psi(1:npts)**2
rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)*dpsi_da(1:npts)/imp_fn(1:npts),psi(1:npts)*dpsi_da(1:npts)/imp_fn(1:npts))
write(6,'(''imp_fn=1+psi^2           rms='',f18.5)') rms
!imp_fn(1:npts)=.01+psi(1:npts)**.5
!rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)*dpsi_da(1:npts)/imp_fn(1:npts),psi(1:npts)*dpsi_da(1:npts)/imp_fn(1:npts))
!write(6,'(''imp_fn=.01+psi^.5        rms='',f18.5)') rms
!imp_fn(1:npts)=psi(1:npts)**.5 + .01*psi(1:npts)
!rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)*dpsi_da(1:npts)/imp_fn(1:npts),psi(1:npts)*dpsi_da(1:npts)/imp_fn(1:npts))
!write(6,'(''imp_fn=psi^.5+.01psi,    rms='',f18.5)') rms
!imp_fn(1:npts)=abs(dpsi_da(1:npts)*h_psi(1:npts)-dpsi_da_h_psi_av_by_dpsi_da_psi_av*dpsi_da(1:npts)*psi(1:npts))+.01*psi(1:npts)**.5
!rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)*dpsi_da(1:npts)/imp_fn(1:npts),psi(1:npts)*dpsi_da(1:npts)/imp_fn(1:npts))
!write(6,'(''imp_fn=Optimal+.01psi^.5 rms='',f18.5)') rms
call imp_fn_sub1(psi,grad_psi,imp_fn)
rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)*dpsi_da(1:npts)/imp_fn(1:npts),psi(1:npts)*dpsi_da(1:npts)/imp_fn(1:npts))
write(6,'(''imp_fn=imp_fn1           rms='',f18.5)') rms
call imp_fn_sub2(psi,grad_psi,imp_fn)
rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)*dpsi_da(1:npts)/imp_fn(1:npts),psi(1:npts)*dpsi_da(1:npts)/imp_fn(1:npts))
write(6,'(''imp_fn=imp_fn2           rms='',f18.5)') rms
imp_fn(1:npts)=abs(dpsi_da(1:npts)*h_psi(1:npts)-dpsi_da_h_psi_av_by_dpsi_da_psi_av*dpsi_da(1:npts)*psi(1:npts))
rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)*dpsi_da(1:npts)/imp_fn(1:npts),psi(1:npts)*dpsi_da(1:npts)/imp_fn(1:npts))
write(6,'(''imp_fn=Optimal -----     rms='',f18.5)') rms
!imp_fn(1:npts)=abs(dpsi_da(1:npts)**1.1*h_psi(1:npts)-dpsi_da_h_psi_av_by_dpsi_da_psi_av*dpsi_da(1:npts)**1.1*psi(1:npts))
!rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)*dpsi_da(1:npts)/imp_fn(1:npts),psi(1:npts)*dpsi_da(1:npts)/imp_fn(1:npts))
!write(6,'(''imp_fn=?,                rms='',f18.5)') rms
write(6,*)

write(6,'(''3: <psi_a H psi> / <psi^2>'',f10.6)') dpsi_da_h_psi_av_by_psisq_av
imp_fn(1:npts)=1
rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)*dpsi_da(1:npts)/imp_fn(1:npts),psi(1:npts)**2/imp_fn(1:npts))
write(6,'(''imp_fn=1,                rms='',f18.5)') rms
imp_fn(1:npts)=psi(1:npts)**.5
rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)*dpsi_da(1:npts)/imp_fn(1:npts),psi(1:npts)**2/imp_fn(1:npts))
write(6,'(''imp_fn=psi**.5,          rms='',f18.5)') rms
imp_fn(1:npts)=psi(1:npts)
rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)*dpsi_da(1:npts)/imp_fn(1:npts),psi(1:npts)**2/imp_fn(1:npts))
write(6,'(''imp_fn=psi,              rms='',f18.5)') rms
imp_fn(1:npts)=psi(1:npts)**2
rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)*dpsi_da(1:npts)/imp_fn(1:npts),psi(1:npts)**2/imp_fn(1:npts))
write(6,'(''imp_fn=psi**2,           rms='',f18.5)') rms
!imp_fn(1:npts)=psi(1:npts)*(1+psi(1:npts))
!rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)*dpsi_da(1:npts)/imp_fn(1:npts),psi(1:npts)**2/imp_fn(1:npts))
!write(6,'(''imp_fn=psi+psi^2,        rms='',f18.5)') rms
imp_fn(1:npts)=1+psi(1:npts)
rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)*dpsi_da(1:npts)/imp_fn(1:npts),psi(1:npts)**2/imp_fn(1:npts))
write(6,'(''imp_fn=1+psi,            rms='',f18.5)') rms
imp_fn(1:npts)=1+psi(1:npts)**2
rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)*dpsi_da(1:npts)/imp_fn(1:npts),psi(1:npts)**2/imp_fn(1:npts))
write(6,'(''imp_fn=1+psi^2           rms='',f18.5)') rms
!imp_fn(1:npts)=.01+psi(1:npts)**.5
!rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)*dpsi_da(1:npts)/imp_fn(1:npts),psi(1:npts)**2/imp_fn(1:npts))
!write(6,'(''imp_fn=.01+psi^.5        rms='',f18.5)') rms
!imp_fn(1:npts)=psi(1:npts)**.5 + .01*psi(1:npts)
!rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)*dpsi_da(1:npts)/imp_fn(1:npts),psi(1:npts)**2/imp_fn(1:npts))
!write(6,'(''imp_fn=psi^.5+.01psi,    rms='',f18.5)') rms
!imp_fn(1:npts)=abs(dpsi_da(1:npts)*h_psi(1:npts)-dpsi_da_h_psi_av_by_psisq_av*psi(1:npts)**2)+.01*psi(1:npts)**.5
!rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)*dpsi_da(1:npts)/imp_fn(1:npts),psi(1:npts)**2/imp_fn(1:npts))
!write(6,'(''imp_fn=Optimal+.01psi^.5 rms='',f18.5)') rms
call imp_fn_sub1(psi,grad_psi,imp_fn)
rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)*dpsi_da(1:npts)/imp_fn(1:npts),psi(1:npts)**2/imp_fn(1:npts))
write(6,'(''imp_fn=imp_fn1           rms='',f18.5)') rms
call imp_fn_sub2(psi,grad_psi,imp_fn)
rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)*dpsi_da(1:npts)/imp_fn(1:npts),psi(1:npts)**2/imp_fn(1:npts))
write(6,'(''imp_fn=imp_fn2           rms='',f18.5)') rms
imp_fn(1:npts)=abs(dpsi_da(1:npts)*h_psi(1:npts)-dpsi_da_h_psi_av_by_psisq_av*psi(1:npts)**2)
rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)*dpsi_da(1:npts)/imp_fn(1:npts),psi(1:npts)**2/imp_fn(1:npts))
write(6,'(''imp_fn=Optimal -----     rms='',f18.5)') rms
write(6,*)

write(6,'(''4: <psi_b H psi> / <psi_b psi>'',f10.6)') dpsi_db_h_psi_av_by_dpsi_db_psi_av
imp_fn(1:npts)=1
rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)*dpsi_db(1:npts)/imp_fn(1:npts),psi(1:npts)*dpsi_db(1:npts)/imp_fn(1:npts))
write(6,'(''imp_fn=1,                rms='',f18.5)') rms
imp_fn(1:npts)=psi(1:npts)**.5
rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)*dpsi_db(1:npts)/imp_fn(1:npts),psi(1:npts)*dpsi_db(1:npts)/imp_fn(1:npts))
write(6,'(''imp_fn=psi**.5,          rms='',f18.5)') rms
imp_fn(1:npts)=psi(1:npts)
rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)*dpsi_db(1:npts)/imp_fn(1:npts),psi(1:npts)*dpsi_db(1:npts)/imp_fn(1:npts))
write(6,'(''imp_fn=psi,              rms='',f18.5)') rms
imp_fn(1:npts)=psi(1:npts)**2
rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)*dpsi_db(1:npts)/imp_fn(1:npts),psi(1:npts)*dpsi_db(1:npts)/imp_fn(1:npts))
write(6,'(''imp_fn=psi**2,           rms='',f18.5)') rms
!imp_fn(1:npts)=psi(1:npts)*(1+psi(1:npts))
!rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)*dpsi_db(1:npts)/imp_fn(1:npts),psi(1:npts)*dpsi_db(1:npts)/imp_fn(1:npts))
!write(6,'(''imp_fn=psi+psi^2,        rms='',f18.5)') rms
imp_fn(1:npts)=1+psi(1:npts)
rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)*dpsi_db(1:npts)/imp_fn(1:npts),psi(1:npts)*dpsi_db(1:npts)/imp_fn(1:npts))
write(6,'(''imp_fn=1+psi,            rms='',f18.5)') rms
imp_fn(1:npts)=1+psi(1:npts)**2
rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)*dpsi_db(1:npts)/imp_fn(1:npts),psi(1:npts)*dpsi_db(1:npts)/imp_fn(1:npts))
write(6,'(''imp_fn=1+psi^2           rms='',f18.5)') rms
!imp_fn(1:npts)=.01+psi(1:npts)**.5
!rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)*dpsi_db(1:npts)/imp_fn(1:npts),psi(1:npts)*dpsi_db(1:npts)/imp_fn(1:npts))
!write(6,'(''imp_fn=.01+psi^.5        rms='',f18.5)') rms
!imp_fn(1:npts)=psi(1:npts)**.5 + .01*psi(1:npts)
!rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)*dpsi_db(1:npts)/imp_fn(1:npts),psi(1:npts)*dpsi_db(1:npts)/imp_fn(1:npts))
!write(6,'(''imp_fn=psi^.5+.01psi,    rms='',f18.5)') rms
!imp_fn(1:npts)=abs(dpsi_db(1:npts)*h_psi(1:npts)-dpsi_db_h_psi_av_by_dpsi_db_psi_av*dpsi_db(1:npts)*psi(1:npts))+.01*psi(1:npts)**.5
!rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)*dpsi_db(1:npts)/imp_fn(1:npts),psi(1:npts)*dpsi_db(1:npts)/imp_fn(1:npts))
!write(6,'(''imp_fn=Optimal+.01psi^.5 rms='',f18.5)') rms
call imp_fn_sub1(psi,grad_psi,imp_fn)
rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)*dpsi_db(1:npts)/imp_fn(1:npts),psi(1:npts)*dpsi_db(1:npts)/imp_fn(1:npts))
write(6,'(''imp_fn=imp_fn1           rms='',f18.5)') rms
call imp_fn_sub2(psi,grad_psi,imp_fn)
rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)*dpsi_db(1:npts)/imp_fn(1:npts),psi(1:npts)*dpsi_db(1:npts)/imp_fn(1:npts))
write(6,'(''imp_fn=imp_fn2           rms='',f18.5)') rms
imp_fn(1:npts)=abs(dpsi_db(1:npts)*h_psi(1:npts)-dpsi_db_h_psi_av_by_dpsi_db_psi_av*dpsi_db(1:npts)*psi(1:npts))
rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)*dpsi_db(1:npts)/imp_fn(1:npts),psi(1:npts)*dpsi_db(1:npts)/imp_fn(1:npts))
write(6,'(''imp_fn=Optimal -----     rms='',f18.5)') rms
!imp_fn(1:npts)=abs(dpsi_db(1:npts)**1.1*h_psi(1:npts)-dpsi_db_h_psi_av_by_dpsi_db_psi_av*dpsi_db(1:npts)**1.1*psi(1:npts))
!rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)*dpsi_db(1:npts)/imp_fn(1:npts),psi(1:npts)*dpsi_db(1:npts)/imp_fn(1:npts))
!write(6,'(''imp_fn=?,                rms='',f18.5)') rms
write(6,*)

write(6,'(''5: <psi_b H psi> / <psi^2>'',f10.6)') dpsi_db_h_psi_av_by_psisq_av
imp_fn(1:npts)=1
rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)*dpsi_db(1:npts)/imp_fn(1:npts),psi(1:npts)**2/imp_fn(1:npts))
write(6,'(''imp_fn=1,                rms='',f18.5)') rms
imp_fn(1:npts)=psi(1:npts)**.5
rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)*dpsi_db(1:npts)/imp_fn(1:npts),psi(1:npts)**2/imp_fn(1:npts))
write(6,'(''imp_fn=psi**.5,          rms='',f18.5)') rms
imp_fn(1:npts)=psi(1:npts)
rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)*dpsi_db(1:npts)/imp_fn(1:npts),psi(1:npts)**2/imp_fn(1:npts))
write(6,'(''imp_fn=psi,              rms='',f18.5)') rms
imp_fn(1:npts)=psi(1:npts)**2
rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)*dpsi_db(1:npts)/imp_fn(1:npts),psi(1:npts)**2/imp_fn(1:npts))
write(6,'(''imp_fn=psi**2,           rms='',f18.5)') rms
!imp_fn(1:npts)=psi(1:npts)*(1+psi(1:npts))
!rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)*dpsi_db(1:npts)/imp_fn(1:npts),psi(1:npts)**2/imp_fn(1:npts))
!write(6,'(''imp_fn=psi+psi^2,        rms='',f18.5)') rms
imp_fn(1:npts)=1+psi(1:npts)
rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)*dpsi_db(1:npts)/imp_fn(1:npts),psi(1:npts)**2/imp_fn(1:npts))
write(6,'(''imp_fn=1+psi,            rms='',f18.5)') rms
imp_fn(1:npts)=1+psi(1:npts)**2
rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)*dpsi_db(1:npts)/imp_fn(1:npts),psi(1:npts)**2/imp_fn(1:npts))
write(6,'(''imp_fn=1+psi^2           rms='',f18.5)') rms
!imp_fn(1:npts)=.01+psi(1:npts)**.5
!rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)*dpsi_db(1:npts)/imp_fn(1:npts),psi(1:npts)**2/imp_fn(1:npts))
!write(6,'(''imp_fn=.01+psi^.5        rms='',f18.5)') rms
!imp_fn(1:npts)=psi(1:npts)**.5 + .01*psi(1:npts)
!rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)*dpsi_db(1:npts)/imp_fn(1:npts),psi(1:npts)**2/imp_fn(1:npts))
!write(6,'(''imp_fn=psi^.5+.01psi,    rms='',f18.5)') rms
!imp_fn(1:npts)=abs(dpsi_db(1:npts)*h_psi(1:npts)-dpsi_db_h_psi_av_by_psisq_av*psi(1:npts)**2)+.01*psi(1:npts)**.5
!rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)*dpsi_db(1:npts)/imp_fn(1:npts),psi(1:npts)**2/imp_fn(1:npts))
!write(6,'(''imp_fn=Optimal+.01psi^.5 rms='',f18.5)') rms
call imp_fn_sub1(psi,grad_psi,imp_fn)
rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)*dpsi_db(1:npts)/imp_fn(1:npts),psi(1:npts)**2/imp_fn(1:npts))
write(6,'(''imp_fn=imp_fn1           rms='',f18.5)') rms
call imp_fn_sub2(psi,grad_psi,imp_fn)
rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)*dpsi_db(1:npts)/imp_fn(1:npts),psi(1:npts)**2/imp_fn(1:npts))
write(6,'(''imp_fn=imp_fn2           rms='',f18.5)') rms
imp_fn(1:npts)=abs(dpsi_db(1:npts)*h_psi(1:npts)-dpsi_db_h_psi_av_by_psisq_av*psi(1:npts)**2)
rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)*dpsi_db(1:npts)/imp_fn(1:npts),psi(1:npts)**2/imp_fn(1:npts))
write(6,'(''imp_fn=Optimal -----     rms='',f18.5)') rms
write(6,*)

write(6,'(''6: <psi_b (H-E) psi> / <psi^2>'',f10.6)') dpsi_db_h_psi_av_by_psisq_av
imp_fn(1:npts)=1
rms=rms_deviation_ratio(npts,imp_fn(1:npts),(E_L(1:npts)-e)*psi(1:npts)*dpsi_db(1:npts)/imp_fn(1:npts),psi(1:npts)**2/imp_fn(1:npts))
write(6,'(''imp_fn=1,                rms='',f18.5)') rms
imp_fn(1:npts)=psi(1:npts)**.5
rms=rms_deviation_ratio(npts,imp_fn(1:npts),(E_L(1:npts)-e)*psi(1:npts)*dpsi_db(1:npts)/imp_fn(1:npts),psi(1:npts)**2/imp_fn(1:npts))
write(6,'(''imp_fn=psi**.5,          rms='',f18.5)') rms
imp_fn(1:npts)=psi(1:npts)
rms=rms_deviation_ratio(npts,imp_fn(1:npts),(E_L(1:npts)-e)*psi(1:npts)*dpsi_db(1:npts)/imp_fn(1:npts),psi(1:npts)**2/imp_fn(1:npts))
write(6,'(''imp_fn=psi,              rms='',f18.5)') rms
imp_fn(1:npts)=psi(1:npts)**2
rms=rms_deviation_ratio(npts,imp_fn(1:npts),(E_L(1:npts)-e)*psi(1:npts)*dpsi_db(1:npts)/imp_fn(1:npts),psi(1:npts)**2/imp_fn(1:npts))
write(6,'(''imp_fn=psi**2,           rms='',f18.5)') rms
!imp_fn(1:npts)=psi(1:npts)*(1+psi(1:npts))
!rms=rms_deviation_ratio(npts,imp_fn(1:npts),(E_L(1:npts)-e)*psi(1:npts)*dpsi_db(1:npts)/imp_fn(1:npts),psi(1:npts)**2/imp_fn(1:npts))
!write(6,'(''imp_fn=psi+psi^2,        rms='',f18.5)') rms
imp_fn(1:npts)=1+psi(1:npts)
rms=rms_deviation_ratio(npts,imp_fn(1:npts),(E_L(1:npts)-e)*psi(1:npts)*dpsi_db(1:npts)/imp_fn(1:npts),psi(1:npts)**2/imp_fn(1:npts))
write(6,'(''imp_fn=1+psi,            rms='',f18.5)') rms
imp_fn(1:npts)=1+psi(1:npts)**2
rms=rms_deviation_ratio(npts,imp_fn(1:npts),(E_L(1:npts)-e)*psi(1:npts)*dpsi_db(1:npts)/imp_fn(1:npts),psi(1:npts)**2/imp_fn(1:npts))
write(6,'(''imp_fn=1+psi^2           rms='',f18.5)') rms
!imp_fn(1:npts)=.01+psi(1:npts)**.5
!rms=rms_deviation_ratio(npts,imp_fn(1:npts),(E_L(1:npts)-e)*psi(1:npts)*dpsi_db(1:npts)/imp_fn(1:npts),psi(1:npts)**2/imp_fn(1:npts))
!write(6,'(''imp_fn=.01+psi^.5        rms='',f18.5)') rms
!imp_fn(1:npts)=psi(1:npts)**.5 + .01*psi(1:npts)
!rms=rms_deviation_ratio(npts,imp_fn(1:npts),(E_L(1:npts)-e)*psi(1:npts)*dpsi_db(1:npts)/imp_fn(1:npts),psi(1:npts)**2/imp_fn(1:npts))
!write(6,'(''imp_fn=psi^.5+.01psi,    rms='',f18.5)') rms
!imp_fn(1:npts)=abs(dpsi_db(1:npts)*h_psi(1:npts)-dpsi_db_h_psi_av_by_psisq_av*psi(1:npts)**2)+.01*psi(1:npts)**.5
!rms=rms_deviation_ratio(npts,imp_fn(1:npts),(E_L(1:npts)-e)*psi(1:npts)*dpsi_db(1:npts)/imp_fn(1:npts),psi(1:npts)**2/imp_fn(1:npts))
!write(6,'(''imp_fn=Optimal+.01psi^.5 rms='',f18.5)') rms
call imp_fn_sub1(psi,grad_psi,imp_fn)
rms=rms_deviation_ratio(npts,imp_fn(1:npts),(E_L(1:npts)-e)*psi(1:npts)*dpsi_db(1:npts)/imp_fn(1:npts),psi(1:npts)**2/imp_fn(1:npts))
write(6,'(''imp_fn=imp_fn1           rms='',f18.5)') rms
call imp_fn_sub2(psi,grad_psi,imp_fn)
rms=rms_deviation_ratio(npts,imp_fn(1:npts),(E_L(1:npts)-e)*psi(1:npts)*dpsi_db(1:npts)/imp_fn(1:npts),psi(1:npts)**2/imp_fn(1:npts))
write(6,'(''imp_fn=imp_fn2           rms='',f18.5)') rms
imp_fn(1:npts)=abs(dpsi_db(1:npts)*h_minus_e_psi(1:npts)-dpsi_db_h_minus_e_psi_av_by_psisq_av*psi(1:npts)**2)
rms=rms_deviation_ratio(npts,imp_fn(1:npts),(E_L(1:npts)-e)*psi(1:npts)*dpsi_db(1:npts)/imp_fn(1:npts),psi(1:npts)**2/imp_fn(1:npts))
write(6,'(''imp_fn=Optimal -----     rms='',f18.5)') rms
write(6,*)

write(6,'(''8: <psi H psi> / <1>'',f10.6)') e
imp_fn(1:npts)=1
rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)**2/imp_fn(1:npts),one_vec(1:npts)**2/imp_fn(1:npts))
write(6,'(''imp_fn=1,                rms='',f18.5)') rms
imp_fn(1:npts)=psi(1:npts)**.5
rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)**2/imp_fn(1:npts),one_vec(1:npts)**2/imp_fn(1:npts))
write(6,'(''imp_fn=psi**.5,          rms='',f18.5)') rms
imp_fn(1:npts)=psi(1:npts)
rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)**2/imp_fn(1:npts),one_vec(1:npts)**2/imp_fn(1:npts))
write(6,'(''imp_fn=psi,              rms='',f18.5)') rms
imp_fn(1:npts)=psi(1:npts)**2
rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)**2/imp_fn(1:npts),one_vec(1:npts)**2/imp_fn(1:npts))
write(6,'(''imp_fn=psi**2,           rms='',f18.5)') rms
!imp_fn(1:npts)=psi(1:npts)*(1+psi(1:npts))
!rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)**2/imp_fn(1:npts),one_vec(1:npts)**2/imp_fn(1:npts))
!write(6,'(''imp_fn=psi+psi^2,        rms='',f18.5)') rms
imp_fn(1:npts)=1+psi(1:npts)
rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)**2/imp_fn(1:npts),one_vec(1:npts)**2/imp_fn(1:npts))
write(6,'(''imp_fn=1+psi,            rms='',f18.5)') rms
imp_fn(1:npts)=1+psi(1:npts)**2
rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)**2/imp_fn(1:npts),one_vec(1:npts)**2/imp_fn(1:npts))
write(6,'(''imp_fn=1+psi^2           rms='',f18.5)') rms
!imp_fn(1:npts)=.01+psi(1:npts)**.5
!rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)**2/imp_fn(1:npts),one_vec(1:npts)**2/imp_fn(1:npts))
!write(6,'(''imp_fn=.01+psi^.5        rms='',f18.5)') rms
!imp_fn(1:npts)=psi(1:npts)**.5 + .01*psi(1:npts)
!rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)**2/imp_fn(1:npts),one_vec(1:npts)**2/imp_fn(1:npts))
!write(6,'(''imp_fn=psi^.5+.01psi,    rms='',f18.5)') rms
!imp_fn(1:npts)=abs(psi(1:npts)*h_psi(1:npts)-e*psi(1:npts)**2)+.01*psi(1:npts)**.5
!rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)**2/imp_fn(1:npts),one_vec(1:npts)**2/imp_fn(1:npts))
!write(6,'(''imp_fn=Optimal+.01psi^.5 rms='',f18.5)') rms
call imp_fn_sub1(psi,grad_psi,imp_fn)
rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)**2/imp_fn(1:npts),one_vec(1:npts)**2/imp_fn(1:npts))
write(6,'(''imp_fn=imp_fn1           rms='',f18.5)') rms
call imp_fn_sub2(psi,grad_psi,imp_fn)
rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)**2/imp_fn(1:npts),one_vec(1:npts)**2/imp_fn(1:npts))
write(6,'(''imp_fn=imp_fn2           rms='',f18.5)') rms
imp_fn(1:npts)=abs(psi(1:npts)*h_psi(1:npts)-e_num/npts)
rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)**2/imp_fn(1:npts),one_vec(1:npts)**2/imp_fn(1:npts))
write(6,'(''imp_fn=Optimal1-----     rms='',f18.5)') rms
imp_fn(1:npts)=abs(psi(1:npts)*h_psi(1:npts))
rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)**2/imp_fn(1:npts),one_vec(1:npts)**2/imp_fn(1:npts))
write(6,'(''imp_fn=Optimal2-----     rms='',f18.5)') rms
!imp_fn(1:npts)=abs(psi(1:npts)**1.1*h_psi(1:npts)-e*psi(1:npts)**2.1)
!rms=rms_deviation_ratio(npts,imp_fn(1:npts),E_L(1:npts)*psi(1:npts)**2/imp_fn(1:npts),psi(1:npts)**2/imp_fn(1:npts))
!write(6,'(''imp_fn=?,                rms='',f18.5)') rms
write(6,*)

end subroutine read_input
!-----------------------------------------------------------------------

subroutine vec_vec_elem_mul(ndim,vec1_in,vec2_in,vec_out)
! Elements of output vector are products of the corresponding elements of the 2 input vectors
integer, intent(in) :: ndim
real(rk), intent(in) :: vec1_in(ndim), vec2_in(ndim)
real(rk), intent(out) :: vec_out(ndim)

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
real(rk), intent(in) :: E_L(ndim), rho(ndim)
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

write(6,'(''E_L_av='',9f10.6)') E_L_av

end function rms_deviation
!-----------------------------------------------------------------------

subroutine covariance(ndim,rho,f,g, rho_av, f_av, f2_var, g_av, g2_var, fg_var)
! Given an unnormalized probability density rho, calculate the average, variance and covariance of 2 random variables, f and g.
integer, intent(in) :: ndim
real(rk), intent(in) :: rho(ndim), f(ndim), g(ndim)
real(rk), intent(out) :: rho_av, f_av, f2_var, g_av, g2_var, fg_var
real(rk) f2_av, g2_av, fg_av
integer i

rho_av=0; f_av=0; f2_av=0; g_av=0; g2_av=0; fg_av=0
do i=1,ndim
  rho_av=rho_av+rho(i)
  f_av=f_av+rho(i)*f(i)
  f2_av=f2_av+rho(i)*f(i)**2
  g_av=g_av+rho(i)*g(i)
  g2_av=g2_av+rho(i)*g(i)**2
  fg_av=fg_av+rho(i)*f(i)*g(i)
enddo
f_av=f_av/rho_av
f2_av=f2_av/rho_av
g_av=g_av/rho_av
g2_av=g2_av/rho_av
fg_av=fg_av/rho_av
f2_var=f2_av-f_av**2
g2_var=g2_av-g_av**2
fg_var=fg_av-f_av*g_av

! Check that the ratio of expectation values is indep of the importance function
!write(6,'(''f_av/g_av='',9f10.6)') f_av/g_av

end subroutine covariance
!-----------------------------------------------------------------------

function rms_deviation_ratio(ndim,rho,f,g)
! Given an unnormalized probability density rho, calculate the rms deviation of the ratio of 2 random variables, f and g.
integer, intent(in) :: ndim
real(rk), intent(in) :: rho(ndim), f(ndim), g(ndim)
real(rk) rho_av, f_av, f2_var, g_av, g2_var, fg_var, rms_deviation_ratio

call covariance(ndim,rho,f,g, rho_av, f_av, f2_var, g_av, g2_var, fg_var)
rms_deviation_ratio=abs(f_av/g_av)*sqrt(f2_var/f_av**2+g2_var/g_av**2-2*fg_var/(f_av*g_av))

end function rms_deviation_ratio
!-----------------------------------------------------------------------

function error_of_mean(ndim,w,x)
! Warning this routine is probably wrong, but is not used so it does not matter.
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

function psi_func(x)
! Trial wavefunction
! Goes to zero linearly at 0 and 1, but if a=b=1 then it goes to zero quadratically at 1
real(rk), intent(in) :: x
real(rk) psi_func, pi

pi=4*datan(1._rk)

if(ifunc.eq.1) then
  psi_func=sin(pi*(x+a*(x**2-b*x**3)))
elseif(ifunc.eq.2) then
  psi_func=sin(pi*(x+a*(x**2-2*x**3+b*x**4)))
else
  psi_func=sin(pi*(x+a*(8*x-b*x**4)))
endif

end function psi_func
!-----------------------------------------------------------------------
function grad_psi_func(x)
! Gradient of wavefunction
real(rk), intent(in) :: x
real(rk) grad_psi_func, pi

pi=4*datan(1._rk)

if(ifunc.eq.1) then
  grad_psi_func=pi*((1+a*(2*x-3*b*x**2))*cos(pi*(x+a*(x**2-b*x**3))))
elseif(ifunc.eq.2) then
  grad_psi_func=pi*((1+a*(2*x-6*x**2+4*b*x**3))*cos(pi*(x+a*(x**2-2*x**3+b*x**4))))
else
  grad_psi_func=pi*(1+a*(8-4*b*x**3))*cos(pi*(x+a*(8*x-b*x**4)))
endif

end function grad_psi_func
!-----------------------------------------------------------------------
function h_psi_func(x)
! Hamiltonian acting on trial wavefunction
real(rk), intent(in) :: x
real(rk) h_psi_func, pi

pi=4*datan(1._rk)

if(ifunc.eq.1) then
  h_psi_func=.5_rk*pi*(pi*(1+a*(2*x-3*b*x**2))**2 *sin(pi*(x+a*(x**2-b*x**3))) - a*(2-6*x)*cos(pi*(x+a*(x**2-b*x**3))))
elseif(ifunc.eq.2) then
  h_psi_func=.5_rk*pi*(pi*(1+a*(2*x-6*x**2+4*b*x**3))**2 *sin(pi*(x+a*(x**2-2*x**3+b*x**4))) - a*(2-12*x+12*x**2)*cos(pi*(x+a*(x**2-2*x**3+b*x**4))))
else
  h_psi_func=.5_rk*pi*(pi*(1+a*(8-4*b*x**3))**2* sin(pi*(x+a*(8*x-b*x**4))) + 12*a*b*x**2*cos(pi*(x+a*(8*x-b*x**4))))
endif

end function h_psi_func
!-----------------------------------------------------------------------
function dpsi_da_func(x)
! Hamiltonian acting on trial wavefunction
real(rk), intent(in) :: x
real(rk) dpsi_da_func, pi

pi=4*datan(1._rk)

if(ifunc.eq.1) then
  dpsi_da_func=pi*(x**2-b*x**3)*cos(pi*(x+a*(x**2-b*x**3)))
elseif(ifunc.eq.2) then
  dpsi_da_func=pi*(x**2-2*x**3+b*x**4)*cos(pi*(x+a*(x**2-2*x**3+b*x**4)))
else
  dpsi_da_func=pi*(8*x-b*x**4)*cos(pi*(x+a*(8*x-b*x**4)))
endif

end function dpsi_da_func
!-----------------------------------------------------------------------
function dpsi_db_func(x)
! Hamiltonian acting on trial wavefunction
real(rk), intent(in) :: x
real(rk) dpsi_db_func, pi

pi=4*datan(1._rk)

if(ifunc.eq.1) then
  dpsi_db_func=-pi*a*x**3*cos(pi*(x+a*(x**2-b*x**3)))
elseif(ifunc.eq.2) then
  dpsi_db_func=pi*a*x**4*cos(pi*(x+a*(x**2-2*x**3+b*x**4)))
else
  dpsi_db_func=-pi*a*x**4*cos(pi*(x+a*(8*x-b*x**4)))
endif

end function dpsi_db_func
!-----------------------------------------------------------------------

subroutine imp_fn_sub1(psi,grad_psi,imp_fn)
! Importance function that is nonzero at nodes of psi_T
real(rk), intent(in) :: psi(npts),grad_psi(npts)
real(rk), intent(out) :: imp_fn(npts)
real(rk) d,f
integer i
f(d)=.5_rk*eps*(1/d+d/eps**2)

do i=1,npts
  d=abs(psi(i)/grad_psi(i))
!  write(6,*) d,f(d)
  if(d.ge.eps) then
    imp_fn(i)=psi(i)
  else
    imp_fn(i)=f(d)*psi(i)
  endif
enddo
!write(6,'(''psi=   '',1000f10.6)') psi(1:npts)**2
!write(6,'(''imp_fn='',1000f10.6)') imp_fn(1:npts)

end subroutine imp_fn_sub1
!-----------------------------------------------------------------------

subroutine imp_fn_sub2(psi,grad_psi,imp_fn)
! Importance function that is nonzero at nodes of psi_T
real(rk), intent(in) :: psi(npts),grad_psi(npts)
real(rk), intent(out) :: imp_fn(npts)
real(rk) d,f
integer i
f(d)=.5_rk*eps*(1/d+d/eps**2)

do i=1,npts
  d=abs(psi(i)/grad_psi(i))
!  write(6,*) d,f(d)
  if(d.ge.eps) then
    imp_fn(i)=psi(i)**2
  else
    imp_fn(i)=(f(d)*psi(i))**2
  endif
enddo
!write(6,'(''psi=   '',1000f10.6)') psi(1:npts)**2
!write(6,'(''imp_fn='',1000f10.6)') imp_fn(1:npts)

end subroutine imp_fn_sub2
!-----------------------------------------------------------------------

end module imp_sam
!-----------------------------------------------------------------------

program imp_sampl
use imp_sam, only: read_input

call read_input

end program imp_sampl
