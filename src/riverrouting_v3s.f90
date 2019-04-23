module riverrouting
use variables
use channel
contains

subroutine rivers(rstart,rend)
integer rstart, rend
!Loop over all channels (planes must be done first, see above)
!This loop is dependent on upstream inflows, use domains to parallelize
do j=rstart,rend    

	qlat_ch=ex_s(j)*length_p(j)*2.0 !(fps*Lp)/ft/plane

	qlat_ch_ave=0.5*(qlat_ch+qlat_ch_old(j)) !added 11/19/09
    
	!inflow = will be sum of upstream outflows;
	if (nup(j).eq.0)then
		q_in(j)= 0. !old_q(j,1) !first DX q
	elseif(nup(j).eq.1)then
		q_in(j)=old_q(uppfaf(j,1),ndx)
	elseif(nup(j).eq.2)then
		q_in(j)=old_q(uppfaf(j,1),ndx)+old_q(uppfaf(j,2),ndx)
	elseif(nup(j).eq.3)then
		q_in(j)=old_q(uppfaf(j,1),ndx)+old_q(uppfaf(j,2),ndx)+old_q(uppfaf(j,3),ndx)
	elseif(nup(j).eq.4)then
		q_in(j)=old_q(uppfaf(j,1),ndx)+old_q(uppfaf(j,2),ndx)+old_q(uppfaf(j,3),ndx)+old_q(uppfaf(j,4),ndx)
	endif
    
	call route_mc_ch(j)
    
    qlat_ch_old(j)=qlat_ch 
    
enddo

return
end subroutine rivers

end module riverrouting
