module channel
use variables
contains

subroutine route_mc_ch(j)
! 	this subroutine routes the main channel for the current reach for a set time
!	and distance step. save the initial and final inflow discharge values for the main
! 	channel sections.  they are needed to compute the initial and final channel
! 	storage values and to perform the mass balance calculations.
!
! 	Flow routing in the main channel.
	integer j, k
			
	do k=1,ndx	!k is the space index, and indicates the position being calculated.
        
		q_in_old(j)=old_q_ch_in(j,k)   !old discharges defined in the previous time step
		q_out_old(j)=old_q_ch_out(j,k)	!old discharges defined in the previous time step
       	
       	call muskin_cp(j,k)   !using constant parameter run
    	
    	old_q_ch_in(j,k)=q_in(j)	!to be used in the next time step
    	old_q_ch_out(j,k)=q_out(j)	!to be used in the next time step
    	
    	q_in(j)=q_out(j)
	enddo

    if (updateMode .eq. 'normal')then
	    old_q(j,ndx)=q_out(j)
    else if (updateMode .eq. 'assimUpdate')then
        write(*,*) j,old_q(j,ndx)
        !old_q(j,ndx)=old_q_res(j,ndx)
    else
        write(*,*) 'invalid argument',updateMode
    end if

	return
end subroutine route_mc_ch
!
!***************************************************************************************************************************
!
subroutine muskin_cp(j,k)
	integer j, k
! this subroutine performs channel routing using the muskingum-cunge method
! with constant-parameter. need to setup the following parameters 
!	topwdt, depth, cross_area, qave, using wide channel assumption
	!mc_topwdt(j)=width_ch(j)	
	!y=(q_ave_ch(j)/(C1_ch(j)*width_ch(j)))**(0.6)
	!mc_cross_area(j)=mc_topwdt(j)*y
	!celert=(5./3.)*q_ave_ch(j)/mc_cross_area(j)
	!tv=dtis/sreach
	!sreach=slen3(j)/real(ndx)
	!c= celert * tv
	!d= q_ave_ch(j) / mc_topwdt(j) / (slope_ch(j) * celert * sreach)
    !cdenom= 1 + c + d
    !cc1= (-1 + c + d) / cdenom
    !cc2= (1 + c - d) / cdenom
    !cc3= (1 - c + d) / cdenom
    !cc4= 2.0 * c/cdenom (Ponce book, P323)
    
    !For constant Parameter, can pre-calculate cc1-cc4 (read in as input)
    !write(*,*) "cc1(j),cc2(j),cc3(j),cc4(j)"
    !write(*,*) cc1(j),cc2(j),cc3(j),cc4(j)
    !write(*,*) "j,k,q_in(j),q_in_old(j),q_out_old(j)"
    !write(*,*) j,k,q_in(j),q_in_old(j),q_out_old(j)
    q_out(j)=cc1(j)*q_in(j)+cc2(j)*q_in_old(j)+cc3(j)*q_out_old(j)+cc4(j)*qlat_ch_ave*(length_ch(j)/ndx)
    
 			
	!Never allow negative discharge
	if (q_out(j).lt.0.0)then
		!write(*,*) j, q_out(j), t
		q_out(j) = 0.
		!stop
	endif
	return
end subroutine muskin_cp

end module channel
