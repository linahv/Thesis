module alpha_modules

use bmad                 ! Define the structures we need to know about.

implicit none


contains

subroutine FNdisp(ring, N, n_slices, nux, F0, F1)

use bmad

implicit none

type (lat_struct) ring   ! This structure holds the lattice info
integer N
integer n_slices
REAL(rp) nux
REAL(rp):: F0(0:N-1)
REAL(rp) :: F1(0:N-1)
!REAL(rp):: phi(0:nmax*n_sicles)

type (ele_struct), pointer :: ele
type (ele_struct) :: ele2
integer nmax
integer i,j
REAL(rp) delta_s ,length_s, s
REAL(rp) phi, beta, rho, k, nu
REAL(rp) k2l, k2

REAL(rp) :: F1dip(0:N-1)
REAL(rp) :: F1quad(0:N-1)
REAL(rp) :: F1sext(0:N-1)
REAL(rp) :: F1multi(0:N-1)
logical :: err, use_last
err = .false.
use_last = .true.

! Propagate Twiss parameters
if (ring%param%geometry == closed$) call twiss_at_start (ring)
call twiss_propagate_all (ring)      

nmax = ring%n_ele_max
!print *, nux

!Initialisation
do n=0,N-1
			F0(n) = 0
		    	F1(n) = 0
			F1dip(n) = 0
			F1quad(n) = 0
			F1sext(n) = 0
			F1multi(n) = 0
enddo

do j=1,nmax

	if (ring%ele(j)%sub_key == sbend$) then
	s = ring%ele(j)%s
	length_s = ring%ele(j)%value(l$)
	delta_s = length_s/n_slices
		do i=1,n_slices
		s = ring%ele(j)%s -(i-1)*delta_s
			call twiss_and_track_at_s(ring,s,ele2)
			rho = ring%ele(j)%value(l$)/ring%ele(j)%value(angle$)
			k = ring%ele(j)%value(k1$)
			phi = ele2%a%phi/nux
			beta= ele2%a%beta
			do n=0,N-1
			F0(n) = F0(n) + 1/rho*sqrt(beta)*cos(n*phi)*delta_s*nux/pi
		    	F1dip(n) = F1dip(n) + sqrt(beta)*(ele2%a%eta)*((1/rho)**2+k-2*k/rho*(ele2%a%eta))*cos(n*phi)*delta_s*nux/pi 
			enddo
		enddo 
	else if (ring%ele(j)%key == quadrupole$) then
	s = ring%ele(j)%s
	length_s = ring%ele(j)%value(l$)
	delta_s = length_s/n_slices
		do i=1,n_slices
		s = ring%ele(j)%s -(i-1)*delta_s
		CALL twiss_and_track_at_s(ring,s,ele2)
			if(err) exit
			phi = ele2%a%phi/nux
			beta= ele2%a%beta
			k = ele2%value(k1$)
			do n=0,N-1
		    F1quad(n) = F1quad(n) + sqrt(beta)*(ele2%a%eta)*k*cos((n)*phi)*delta_s*nux/pi
			enddo
		enddo
	else if (ring%ele(j)%key == sextupole$) then
	s = ring%ele(j)%s
	length_s = ring%ele(j)%value(l$)
	delta_s = length_s/n_slices
		do i=1,n_slices
		s = ring%ele(j)%s -(i-1)*delta_s
			CALL twiss_and_track_at_s(ring,s,ele2)
			phi = ele2%a%phi/nux
			beta= ele2%a%beta
			k2 = ele2%value(k2$)
			!print *, k2
			do n=0,N-1
		    F1sext(n) = F1sext(n) - sqrt(beta)*(ele2%a%eta)**2*k2*cos(n*phi)/2*nux/pi*delta_s
			enddo
		enddo
	else if (ring%ele(j)%key == multipole$) then
	s = ring%ele(j)%s
	length_s = ring%ele(j)%value(l$)
	!print *, ring%ele(j)%name
		if (abs(ring%ele(j)%a_pole(2))>0.001) then
			!print *, ring%ele(j)%a_pole(:)
			CALL twiss_and_track_at_s(ring,s,ele2)
			phi = ele2%a%phi/nux
			beta= ele2%a%beta
			k2l = ring%ele(j)%a_pole(2)
			do n=0,N-1
		    F1multi(n) = F1multi(n) - sqrt(beta)*(ele2%a%eta)**2*k2l*cos(n*phi)/2*nux/pi
			enddo
		endif
	else
	s = ring%ele(j)%s
	endif
enddo

F1(:)  = F1dip(:)  + F1multi(:)  + F1quad(:)  + F1sext(:) 
F0(0) = F0(0)/2
F1(0) = F1(0)/2
endsubroutine

subroutine F2N(ring, N, n_slices, nux, F0, F1, F2)

use bmad

implicit none

type (lat_struct) ring   ! This structure holds the lattice info
integer N
integer n_slices
REAL(rp) nux
REAL(rp):: F0(0:N-1)
REAL(rp) :: F1(0:N-1)
REAL(rp) :: F2(0:N-1)
!REAL(rp):: phi(0:nmax*n_sicles)

type (ele_struct), pointer :: ele
type (ele_struct) :: ele2
integer nmax
integer i,j,l
REAL(rp) delta_s ,length_s, s, d1
REAL(rp) phi, beta, rho, k, nu
REAL(rp) k2l, k2

logical :: err, use_last
err = .false.
use_last = .true.

! Propagate Twiss parameters
if (ring%param%geometry == closed$) call twiss_at_start (ring)
call twiss_propagate_all (ring)      

nmax = ring%n_ele_max

!Initialisation
do l=0,N-1
			F2(n) = 0
enddo

do j=1,nmax
	if (ring%ele(j)%sub_key == sbend$) then
	s = ring%ele(j)%s
	length_s = ring%ele(j)%value(l$)
	delta_s = length_s/n_slices
		do i=1,n_slices
		d1 = 0
		s = ring%ele(j)%s -(i-1)*delta_s
			call twiss_and_track_at_s(ring,s,ele2)
			rho = ring%ele(j)%value(l$)/ring%ele(j)%value(angle$)
			k = ring%ele(j)%value(k1$)
			phi = ele2%a%phi/nux
			beta= ele2%a%beta
			do l=0,N-1
			d1 = d1 + sqrt(beta)*(F1(l)-F0(l))*cos(l*phi)/(nux**2-l**2)
			enddo
			do n=0,N-1
			F2(n) = F2(n) + d1*sqrt(beta)*(k-4*k/rho*ele2%a%eta)*cos(n*phi)*delta_s*nux/pi	
			enddo
		enddo 
	else if (ring%ele(j)%key == quadrupole$) then
	s = ring%ele(j)%s
	length_s = ring%ele(j)%value(l$)
	delta_s = length_s/n_slices
		do i=1,n_slices
		d1 = 0
		s = ring%ele(j)%s -(i-1)*delta_s
		CALL twiss_and_track_at_s(ring,s,ele2)
			if(err) exit
			phi = ele2%a%phi/nux
			beta= ele2%a%beta
			k = ele2%value(k1$)
			do n=0,N-1
		    F2(n) = F2(n) + sqrt(beta)*d1*k*cos(n*phi)*delta_s*nux/pi
			enddo
		enddo
	else if (ring%ele(j)%key == sextupole$) then
	s = ring%ele(j)%s
	length_s = ring%ele(j)%value(l$)
	delta_s = length_s/n_slices
		do i=1,n_slices
		d1 = 0
		s = ring%ele(j)%s -(i-1)*delta_s
			CALL twiss_and_track_at_s(ring,s,ele2)
			phi = ele2%a%phi/nux
			beta= ele2%a%beta
			k2 = ele2%value(k2$)
			!print *, k2
			do l=0,N-1
			d1 = d1 + sqrt(beta)*(F1(l)-F0(l))*cos(l*phi)/(nux**2-l**2)
			enddo
			do n=0,N-1
		    F2(n) = F2(n) -sqrt(beta)*d1*k2*ele2%a%eta*cos(n*phi)*delta_s*nux/pi
			enddo
		enddo
	else if (ring%ele(j)%key == multipole$) then
	s = ring%ele(j)%s
	length_s = ring%ele(j)%value(l$)
	!print *, ring%ele(j)%name
		if (abs(ring%ele(j)%a_pole(2))>0.001) then
		d1 = 0
			!print *, ring%ele(j)%a_pole(:)
			CALL twiss_and_track_at_s(ring,s,ele2)
			phi = ele2%a%phi/nux
			beta= ele2%a%beta
			k2l = ring%ele(j)%a_pole(2)
			do l=0,N-1
			d1 = d1 + sqrt(beta)*(F1(l)-F0(l))*cos(l*phi)/(nux**2-l**2)
			enddo
			
			do n=0,N-1
		    	F2(n) = F2(n) -sqrt(beta)*d1*k2l*ele2%a%eta*cos(n*phi)*nux/pi
			enddo
		endif
	else
	s = ring%ele(j)%s
	endif
enddo
F2(0) = F2(0)/2
endsubroutine


subroutine alpha0calc(ring, F0, F1, F2, N, n_slices, nux, alpha0, alpha1, alpha2)

type (lat_struct) ring
integer N
integer n_slices
REAL(rp):: F0(0:N-1)
REAL(rp):: F1(0:N-1)
REAL(rp):: F2(0:N-1)
REAL(rp) alpha0, alpha1, alpha2, nux

type (ele_struct), pointer :: ele
type (ele_struct) :: ele2
integer nmax, j, i
integer l
REAL(rp) delta_s ,length_s, s
REAL(rp) phi, rho, beta, alpha
REAL(rp) a, a1, a2, a3, a4, C0, d, dp, d1, dp1, dp2, d2
REAL(rp) disp0, disp1, disp2, disp1p, disp0p2, disp0pdisp1p
REAL(rp):: ones(n_slices:1)

if (ring%param%geometry == closed$) call twiss_at_start (ring)
call twiss_propagate_all (ring)      

nmax = ring%n_ele_max
C0 = ring%param%total_length
a = 0
a1 = 0
a2 =0
a3 = 0
a4 = 0

do j=1,nmax
	disp0 = 0
	disp1 = 0
	disp0p2 = 0
	disp2 = 0
	disp0pdisp1p = 0
	s = ring%ele(j)%s
	length_s = ring%ele(j)%value(l$)
	delta_s = length_s/n_slices
	do i=1,n_slices
		dp = 0
		d=0
		dp = 0
		d1 = 0
		dp1 = 0
		dp2 = 0
		d2 = 0
		dp = 0
		s = ring%ele(j)%s - (i-1)*delta_s
			call twiss_and_track_at_s(ring,s,ele2)
			phi = ele2%a%phi/nux
			beta= ele2%a%beta
			alpha= ele2%a%alpha
			do l=0,N-1
			d= d + sqrt(beta)*F0(l)*cos(l*phi)/(nux**2-l**2)
			dp = dp - l*sqrt(beta)*F0(l)*sin(l*phi)/(nux**2-l**2)/nux/beta
			d1 = d1 + sqrt(beta)*(F1(l)-F0(l))*cos(l*phi)/(nux**2-l**2)
			dp1 = dp1 - l*sqrt(beta)*(F1(l))*sin(l*phi)/(nux**2-l**2)/nux/beta
			d2 = d2 + sqrt(beta)*(F0(l)-F1(l)+F2(l))*cos(l*phi)/(nux**2-l**2)
			enddo
			dp = dp -alpha/beta*d
			dp1 = dp1-dp -alpha/beta*(d + d1)
			disp0 = disp0 + d
			disp1 = disp1 + d1
			disp0p2 = disp0p2 + (dp)**2
			disp2 = disp2 + d2
			disp0pdisp1p = disp0pdisp1p + dp*dp1
	enddo		
	a1 = a1 + (disp0p2)/2*length_s/n_slices
	a3 = a3 + disp0pdisp1p*length_s/n_slices
	if (ring%ele(j)%sub_key == sbend$) then
	rho = ring%ele(j)%value(l$)/ring%ele(j)%value(angle$)
	a = a + disp0*length_s/rho/n_slices
	a2 = a2 + disp1/rho*length_s/n_slices
	a4 = a4 + (disp2- disp0*disp0p2/2)/rho*length_s/n_slices
	endif
enddo

alpha0 = a/C0
alpha1 = (a1+a2)/C0
alpha2 = (a3+a4)/C0
endsubroutine

subroutine getalpha(ring, N, n_slices, nux, alpha0, alpha1, alpha2)
type (lat_struct) ring
integer N
integer n_slices
REAL(rp) alpha0, alpha1, alpha2, nux

REAL(rp):: F0(0:N-1)
REAL(rp):: F1(0:N-1)
REAL(rp):: F2(0:N-1)

call FNdisp(ring, N,n_slices, nux, F0,F1)
call F2N(ring, N, n_slices, nux, F0, F1, F2)
call alpha0calc(ring, F0, F1, F2, N, n_slices, nux, alpha0, alpha1, alpha2)
endsubroutine

endmodule