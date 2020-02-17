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
        open(180,file=trim(roffFile))
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
                getpval=getpval*qlat_all !mm/day
            endif
            
            !getpval = 1.0 !mm/day
            
            qsoff=(getpval/10/2.54/12/24) !runoff mm/day to ft/hr
            qsoff=qsoff/dtis !ft/DT to ft/s
            
            !write(*,*) qsoff
            ex_s(j)=qsoff
    enddo 
  
    !update date/time for next DT only used after spinup period
    if(tt.eq.1) then 
            iday=iday+1
            jday=jday+1
            if(iday.eq.29.AND.imonth.eq.2) then
                leapyrint=int(iyear/4.) !integer
                leapyrft=float(iyear)/4. !float
                leapyrdiff=leapyrft-leapyrint !on leap year diff = 0.0
                if(leapyrdiff.eq.0) then
                    iday=29 !do leap year, feb 29
                else
                    iday=1
                    imonth=imonth+1
                endif
            elseif(iday.eq.30.AND.imonth.eq.2) then !just finished Feb 29
                    iday=1
                    imonth=imonth+1
            elseif(iday.eq.31) then
                    Select Case (imonth)
                            case (1,3,5,7,8,10,12)
                                    iday=iday
                            case (4,6,9,11)
                                    iday=1
                                    imonth=imonth+1
                    End Select
            elseif(iday.eq.32) then
                    Select Case (imonth)
                            case (1,3,5,7,8,10)
                                    iday=1
                                    imonth=imonth+1
                            case (12)
                                    iday=1
                                    imonth=1
                                    iyear=iyear+1
                                    jday=1
                    End Select
            endif
    endif
     
return
end subroutine excess1

end module excesscalc
