module hillsloperouting
use variables
use surfacerunoff

!use mpi
contains

subroutine hillslopes(hstart,hend)
integer hstart, hend
!Loop over all planes
!determine surface runoff and subsurface runoff to channels
!This loop can be sent to all process without waiting for results
 do j=hstart,hend
    !surface routing
	!check is surface routing is needed for plane 1		
	if(ex_s(j).eq.0.AND.y_pl_s(j,ndx).lt.0.001)then
		q_pl_s(j,ndx)=0
	else 
	   call qsurface(j,alp_pl_s(j),bet_pl_s(j))
	endif
 enddo

return
end subroutine hillslopes

end module hillsloperouting
