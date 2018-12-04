module outputdata
use variables
contains


subroutine outputdata1()
    CHARACTER(LEN=80) :: FMT

    !Get output locations
    fname11 = 'output_calibration.txt'
    open(11,file=fname11,status='OLD')
    read(11,*) numout
    do k=1,numout
            read(11,*) idout(k), dvalint, idgauge(k)
    enddo
    close(11)

    fname12 = 'discharge_cms.txt'
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
    fname12 = 'discharge_cms.txt'
    open(12,file=fname12,status='OLD',POSITION='APPEND')
    FMT = "(00000F15.3)"
    write(FMT(2:6),'(i5.5)') numout
    write(12,FMT) (old_q_day(idout(i)),i=1,numout)
    close(12)    
    
    return
end subroutine resultsout

subroutine restartout()
    character(256) :: buf, resname
    !write restart file
    resname = restartFile
    open(128,file=resname,status='REPLACE')
    write(128,'(a)') "i,k,old_q,old_q_ch_in,old_q_ch_out,qlat_ch_old"
    do i=1, numout
        do k=1, ndx
            write(buf,*) i,",",k,",",old_q(i,ndx),",",old_q_ch_in(i,k),",",old_q_ch_out(i,k),",",qlat_ch_old(i)
            call del_spaces(buf)
            write(128,'(a)') trim(buf)
        enddo
    enddo
    close(12)

    return
end subroutine restartout

subroutine del_spaces(s)
    character (*), intent (inout) :: s
    character (len=len(s)) tmp
    integer i, j
    j = 1
    do i = 1, len(s)
      if (s(i:i)==' ') cycle
      tmp(j:j) = s(i:i)
      j = j + 1
    end do
    s = tmp(1:j-1)
end subroutine del_spaces

end module outputdata
