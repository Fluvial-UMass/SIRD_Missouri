module outputdata
use variables
contains


subroutine outputdata1()
    CHARACTER(LEN=80) :: FMT

    !Get output locations
    fname11 = trim(rootDir)//'/'//'output_calibration.txt'
    open(11,file=fname11,status='OLD')
    read(11,*) numout
    do k=1,numout
            read(11,*) idout(k), dvalint, idgauge(k)
    enddo
    close(11)

    !fname12 = 'discharge_cms.txt'
    fname12 = trim(outDir)//'/discharge_cms.txt'
    open(12,file=fname12,status='REPLACE')
    FMT = "(00000i15)"
    write(FMT(2:6),'(i5.5)') numout
    write(12,FMT) (idgauge(i),i=1,numout)
    close(12)

    return
end subroutine outputdata1

subroutine resultsout()
    CHARACTER(LEN=80) :: FMT
    !Open results file for specified model units
    !fname12 = 'discharge_cms.txt'
    fname12 = trim(outDir)//'/discharge_cms.txt'
    open(12,file=fname12,status='OLD',POSITION='APPEND')
    FMT = "(00000F15.3)"
    write(FMT(2:6),'(i5.5)') numout
    write(12,FMT) (old_q_day(idout(i)),i=1,numout)
    close(12)    
    
    return
end subroutine resultsout

subroutine restartout()
    character(256) :: buf, resname
    real(kind=JPRM),allocatable :: qlat_ch_old_2d(:,:)
    allocate (qlat_ch_old_2d(1:pfafunits,1:ndx))
    !write restart file
    do i=1,pfafunits
        qlat_ch_old_2d(i,:) = qlat_ch_old(i)
    enddo
    resname = trim(outDir)//'/'//trim(restartFile)
    open(12,file=resname,form='unformatted',access='direct',recl=4*pfafunits*ndx)
    write(12,rec=1) old_q
    write(12,rec=2) old_q_ch_in
    write(12,rec=3) old_q_ch_out
    write(12,rec=4) qlat_ch_old_2d
    close(12)

    return
end subroutine restartout

end module outputdata
