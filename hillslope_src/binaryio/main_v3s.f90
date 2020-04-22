!****************************************************************************
!  PROGRAM:  Hillslope River Routing (HRR) Model Verison 2b Parallel
!  PURPOSE:  Entry point for the console application.
!****************************************************************************
program PFAF_Model
    use variables
    use inputdata
    use outputdata
    use excesscalc
    use hillsloperouting
    use riverrouting
    use channel

    character(128)    ::  mode

    !get mode type from command line
    call getarg(1,mode)
	!mode = 'restart'

    !get main model parmaters and build arrays
    call inputdata1()
    call allocatenow()
    call setzero()

    !get planes and channel data
    call inputdata2(mode)

    !Prepare output files for results
    call outputdata1()

    !write(*,*) 'Starting HRR'
    !##############################################################################
    !Start of Main time loop
    daynum=1 !start run on day 1
    tt=1 !hour of day at start of the run
    jday=1

    do t=1,ndt
        !Get P and ET; determine surface runoff and drainage to GW
        if (tt.eq.1) then
            call excess1()
        endif

        !Loop over all planes: surface and subsurface runoff to channels
        call hillslopes(1,pfafunits)
        call rivers(1,pfafunits)

        !Write out hydrographs at select nodes: cfs at each DT
        !do i=1, numout
            !for daily output
        !    old_q_day(idout(i)) = old_q_day(idout(i)) + old_q(idout(i),ndx)
        !enddo

        do i=1, pfafunits
            !for daily output
            old_q_day(i) = old_q_day(i) + old_q(i,ndx)
        enddo

        if(tt.eq.24)then
            !for daily output
            !do i=1, numout
            !        old_q_day(idout(i)) = (old_q_day(idout(i))/24.)*0.3048**3
            !enddo

            do i=1, pfafunits
                    old_q_day(i) = (old_q_day(i)/24.)*0.3048**3
                    if(old_q_day(i).ge.1000000000000.0) write(*,*) i, &
                    iyear, imonth, iday, old_q_day(i)
            enddo

            tt = 0 !will be 1 after endif
            daynum = daynum+1

            ! if(iday.eq.1.AND.imonth.eq.1)then
                ! write(*,*) 'yr mo day ', iyear, imonth, iday, 'outlet Q(cms) ', old_q_day(28931)
            ! endif
            call restartout()
            call resultsout()
            ! do i=1, numout
            !     old_q_day(idout(i)) = 0.
            ! enddo
            do i=1, pfafunits
                old_q_day(i) = 0.
            enddo
        endif
        tt = tt + 1 !hr counter

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


enddo
!##############################################################################
!End Main time loop

end program PFAF_Model
