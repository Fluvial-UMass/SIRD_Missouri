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

    !get main model parmaters and build arrays
    call inputdata1()
    call allocatenow()
    call setzero()

    !get planes and channel data
    call inputdata2(mode)

    !Prepare output files for results
    call outputdata1()

    !write(*,*) 'Starting HRR'
    write(*, *) "read in:"
    write(*, *) trim(roffFile)
    write(*, *) "restart from:"
    write(*, *) trim(outDir)//'/'//trim(restartFile)
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
        do i=1, numout
            !for daily output
            old_q_day(idout(i)) = old_q_day(idout(i)) + old_q(idout(i),ndx)
        enddo

        if(tt.eq.24)then
            !for daily output
            do i=1, numout
                    old_q_day(idout(i)) = (old_q_day(idout(i))/24.)*0.3048**3
            enddo

            tt = 0 !will be 1 after endif
            daynum = daynum+1

            !if(iday.eq.1.AND.imonth.eq.1)then
                write(*,*) 'yr mo day ', iyear, imonth, iday, 'outlet Q(cms) ', old_q_day(idout(pfafunits))
            !endif
            call restartout()
            call resultsout()
            do i=1, numout
                old_q_day(idout(i)) = 0.
            enddo
        endif
        tt = tt + 1 !hr counter
enddo
!##############################################################################
!End Main time loop
write(*, *) ""
end program PFAF_Model
