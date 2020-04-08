module excesscalc
use variables
use outputdata
use inputdata

contains

subroutine excess1() !daily WBM output, 1-hr Routing
integer j
real getpval, qsoff
character*100 fname180
character*4  cYYYY
character  sval

    !Open Qlat file data
    if (tt.eq.1.AND.t.eq.1) then
        open(180,file=roffFile)
        read(180,*) sval
    endif
        
    do j=1,pfafunits
            !get qlat data from file
            !getpval units m3/s/m
            read(180,*) dval, dval, getpval
            if(getpval.eq.-9999) then
                !qsoff = 0.3/(1000**2)/0.3048 !m3/s/km2 to ft/s
                getpval = 0.3 !mm/day
            else
                !qsoff=(getpval*(length_ch(j)*0.3048))/(A(j)*1000**2)/0.3048 !m3/s/km2 to ft/s
                !qsoff=getpval/(1000**2)/0.3048 !m3/s/km2 to ft/s
                getpval=getpval !mm/day
            endif
            
            !getpval = 1.0 !mm/day
            
            qsoff=(getpval/10/2.54/12/24) !runoff mm/day to ft/hr
            qsoff=qsoff/dtis !ft/DT to ft/s
            
            !write(*,*) qsoff
            ex_s(j)=qsoff
    enddo 
  
         
return
end subroutine excess1

end module excesscalc
