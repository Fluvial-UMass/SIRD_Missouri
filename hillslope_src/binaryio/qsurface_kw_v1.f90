module surfacerunoff
use variables
!use mpi
contains
subroutine qsurface(j,alp,bet)
	real alp,bet
	real alat,asum,qe,qlatp,qave,qup
	real bem,ben,albet,albem,dtxa,dtx
	real error,daq,aest,adev,fder,sder,bb,sc,stem,etem,x1,x2,ad1,ad2
	integer i,j,k1,iter
!
!	Hillslope Routing for Surface Runoff
!	Type of routing to be performed: Plane bordering a channel segment
!	j	subunit identifier
!	k1 	space step index
!	plane parameters
	dtx=dtis*ndx/length_p(j)
	qup = 0.
    shbar = 0.
	
!   begin the kinematic wave routing calculation for the current time step.
    	do k1=1,ndx	! loop over downstream distance steps
!   lateral inflow
        	qlatp=0.5*(ex_s(j)+qlat_s_old(j)) !added 11/19/09
        	alat=qlatp*dtis

!   asum is omega in b-10.
        	asum=alat+y_pl_s(j,k1)+dtx*qup

!   if the rhs of b-9 is essentially zero, that is the lateral inflows, and
!   upstream inflows for this space step are essentially zero, then q(n+1,j+1) is assumed
!   to be zero. The upstream input for the current space step, qup, is also zero.  
!   the water routing routine wrout will be skiped for this dx.
        	if(asum.le.1.0e-06) then
     	    		qe=0.
     	    		y_pl_s(j,k1)=0.
     	    		qup=qe
     	    		goto 350
     	    endif
     	    
!   wrout is the solver for area and discharge for the current dx.
!   call wrout (j,i,k1,dtx,alp,bet,qe,asum,alat) !nonlinear k-wave
!
!   subroutine wrout (j,i,k,dtx,alp,bet,qe,asum,alat)
!   based on li, simons and stevens, 1975, wrr.    
!   this the routine that impliments li and steven's solution for the
!   channel kinematic wave routing of water. it is called to calculate 
!   q(x,t)=q(j+1,n+1) see eq b-2 (eggert, 1980).  nonlinear scheme is 
!   provides an explicit solution to the kinematic wave, and uses a 
!   linear first approximation to speed convergence of the solution 
!   of nonlinear equation for q(j+1,n+1), qe in the code.
!
!   call parameter list:
!   j		= j+1 in eq b-2 (eggert,1980 (for convenience))
!   k		= kth stream reach in the model - so that the correct
!       	  parameters, inflows and lateral inflows are used.
!   dtx	= dtx = dt/dx for this reach
!   alp	= alpha' in eq b-3
!   bet	= beta' in eq b-3
!   qup	= q(j,n+1)
!   qe	= the desired discharge, q(j+1,n+1)
!
!   set up the linear scheme and a few convenient parameters
  	    	qave=0.5*(qup+q_pl_s(j,k1))
  	    	bem=bet-1.
  	    	ben=bem-1.
  	    	albet=alp*bet
  	    	albem=alp*bet*bem
  	    	dtxa=dtx+alp
  	    	error=eps*asum
			
!   linear scheme to find the first approximation initialize the iteration counter
  	    	iter=0
!   first guess... 
!   if the existing and entering flows in the reach are small
!   use a simpler first guess, other wise compute rhs of eq b-25,
!   for the first guess, and start iterating for the current q, qe
  	    	if (qave.gt.0.000001) then
!           	the next two lines are eq b-25
     	    	daq=albet*qave**bem
     	    	qe=(alat+dtx*qup+daq*q_pl_s(j,k1))/(dtx+daq)
	    	else
  		        qe=asum/dtxa
	    	endif
			

!    nonlinear scheme to refine the solution
!    iteration loop
115     	iter=iter+1

!   calculate the new rhs of eq b-17, as the next estimate of qe
  	    	aest=dtx*qe+alp*qe**bet
!   adev is the difference between omega and qe
  	    	adev=asum-aest
!write(*,*) imax, 'imax'	,iter,'iter',j,iyear,imonth,iday		
!   return if converged, with qe as the desired discharge...
  	    	if (abs(adev).le.error) go to 120
  	    	if (iter.ge.imax) then
  		        write(*,*)j,ex_s(j),iyear,imonth,iday,jday,tt
  		        stop 'wrout: max iterations'  
  	    	endif
!   Stop if the maximum iterations is exceeded.
!   otherwise, continue the solution, by forming the rhs of eq b-18,
!   first derivative of f(qe)
116	    	fder=dtx+albet*qe**bem
!   and rhs of eq b-19, 2nd derivative of f(qe),
  	    	sder=albem*qe**ben
!   begin to form the terms of eq b-16 for the next guess of qe
  	    	bb=fder/sder
  	    	sc=2.*adev/sder
  	    	stem=bb*bb+sc
  	    	if (stem.ge.0.) go to 117
!   if root will be imaginary, make a first order truncation guess
  	    	qe=qe+adev/fder
  	    	go to 115
!   otherwise continue with the second order scheme...
117	    	stem=sqrt(stem)
  	    	if (adev.gt.0.) go to 119
  	    	etem=bb+stem
  	    	qe=qe-etem
!   it is possible that the next discharge estimate might be < 0., 
!   in that case, enter the loop below, add fractions of etem back
!   onto qe!!!!!!!if qe is negative, but very near 0,
!   one could spend a lot of time in the 118 loop..., there must be 
!   a better way.  there is no iteration control on this loop.
  	    	if (qe.gt.0.) go to 115
118	    	etem=0.5*etem
  	    	qe=qe+etem
  	    	if (qe.gt.0.) go to 115
  	    	go to 118
!   this section picks the root that produces the smallest error as the next guess
119	    	x1=qe-bb-stem
  	    	x2=qe-bb+stem
  	    	ad1=abs(asum-dtx*x1-alp*x1**bet)
  	    	ad2=abs(asum-dtx*x2-alp*x2**bet)
  	    	qe=x1
  	    	if (ad1.gt.ad2) qe=x2
  	    	go to 115
120	    	continue
 	    	if(qe.lt.0.) qe=0. 
!   End wrout()

!   assign the area and discharge calculated to the arrays storing
!   a and q for the spatial increment.
350	    	y_pl_s(j,k1)=alp*qe**bet
!	q_pl_p will be the previous discharge per unit width or q(j+1,n) from the last 
!	timestep, when the solution proceeds to the next time step.  If k1 = ndx,
!	then q_pl_p(j,i,ndx) is the outflow discharge per unit width of the plane for 
! 	the current time step.
        	q_pl_s(j,k1)=qe
!	Set qup for the next space step
        	qup=qe
            shbar = shbar+y_pl_s(j,k1)
    	enddo ! loop over downstream distance steps

        !if(storageR.eq.1) sv(j) = sv(j)+((shbar/float(ndx))*(304.8)*(1/24.)) !GRACE ft to mm add 1-hr storage for monthly average
        
        qlat_s_old(j)=ex_s(j) !added 11/19/09
        !if(q_pl_s(j,ndx).gt.0.) write(*,*) q_pl_s(j,ndx), y_pl_s(j,ndx)
    	
    	return
end subroutine qsurface

end module surfacerunoff
