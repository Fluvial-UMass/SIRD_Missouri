module inputdata
use variables

contains
!***Read the input.txt for those parameters
subroutine inputdata1()
        character(120):: fname2, fname3
        call getarg(2, fname2)
        !fname2 = 'input.txt'
        open(2,file=fname2)
        read(2,'(a)') updateMode
        read(2,'(a)') rootDir
        read(2,'(a)') roffFile
        read(2,'(a)') channelsFile
        read(2,'(a)') planesFile
        read(2,'(a)') restartFile
        read(2,'(a)') outDir
        !write(*,*) "read in:",roffFile
        read(2,*) pfafunits
        read(2,*) ndx
        read(2,*) ndt
        read(2,*) dtis
        read(2,*) iyear
        read(2,*) imonth
        read(2,*) iday
        read(2,*) jday
        read(2,*) ks_all
        read(2,*) setfsub_rate !used to set baseflow (x mm/day) 1.3
        read(2,*) n_ch_all !channel roughness
        read(2,*) qlat_all
        close(2)

        fname3 = trim(rootDir)//'/'//'output_calibration.txt'
        open(3,file=fname3)
        read(3,*) numout
        close(3)
        return
end subroutine inputdata1

!***read planes and channels data
subroutine inputdata2(mode)
        character(120):: fname3,fname4,mode,restart
        real    Lch_min, setfsub, qperdx, qtop, cbar
        integer :: ji,ki
        real(kind=JPRM),allocatable :: qlat_ch_old_2d(:,:)
        real(kind=JPRM),allocatable :: qlat_s_old_2d(:,:)
        real(kind=JPRM),allocatable :: dummy(:,:)
        allocate (qlat_ch_old_2d(1:pfafunits,1:ndx))
        allocate (qlat_s_old_2d(1:pfafunits,1:ndx))
        allocate (dummy(1:pfafunits,1:ndx))


        if (mode.eq.'restart') then
            restart = trim(outDir)//'/'//trim(restartFile)
            !write(*,*) "restart from:",restart
            open(150,file=restart,form='unformatted',access='direct',recl=4*pfafunits*ndx)
            read(150,rec=1) old_q_res(:,:)
            read(150,rec=2) old_q_ch_in_res(:,:)
            read(150,rec=3) old_q_ch_out_res(:,:)
            read(150,rec=4) qlat_ch_old_2d(:, :)
            read(150,rec=5) y_pl_s(:,:)
            read(150,rec=6) q_pl_s(:,:)
            read(150,rec=7) qlat_s_old_2d(:,:)
            read(150,rec=8) dummy(:,:)
            close(150)
            qlat_ch_old(:) = qlat_ch_old_2d(:, 1)
            qlat_s_old(:) = qlat_ch_old_2d(:,1)
        endif

        fname3 = trim(planesFile)
        fname4 = trim(channelsFile)

        open(3,file=fname3)
        open(4,file=fname4)

        !load the plane data
        cnt = 1
        do j=1,pfafunits

            read(3,*) id(j), slope_p(j),ksurf(j)
            !Added the channel read here to get areas
            read(4,*) id(j),downpfaf(j),nup(j),uppfaf(j,1),uppfaf(j,2),uppfaf(j,3),uppfaf(j,4), &
                A(j),Aup(j),length_ch(j),slope_ch(j),n_ch(j),width_ch(j),Qr_ch(j),dvalint
            !Qr_ch(j)=(1./10./100.)*(Aup(j)*1000**2.0)/24./3600.  !1mm/day to cms
            if(ks_all.ne.1.0) ksurf(j)=ksurf(j)*ks_all
            length_p(j)= A(j)/2.0/length_ch(j) !km
            length_p(j)=length_p(j)*1000.0/0.3048 !km to ft

            if (mode.eq.'initial') then
                !setup initial channel flow
                setfsub=(setfsub_rate/10./2.54/12./3600./24.) !mm/day to ft/s
                qperdx=(A(j)*1000**2/0.3048**2)*setfsub/ndx
                qtop=(Aup(j)-A(j))*(1000**2/0.3048**2)*setfsub
                do k = 1,ndx
                    !setup initial discharge in all channels
                    old_q(j,k)=qtop+qperdx*k
                    old_q_ch_in(j,k)=qtop+qperdx*k !ft3/s
                    old_q_ch_out(j,k)=qtop+qperdx*k !ft3/s
                enddo
            else
                do k = 1,ndx
                    !setup discharge in all channels from restart file
                    old_q(j,k)=old_q_res(j,k)
                    old_q_ch_in(j,k)=old_q_ch_in_res(j,k)
                    old_q_ch_out(j,k)=old_q_ch_out_res(j,k)
                enddo
            end if

            !channel setting
            if(n_ch_all.ne.1.0) n_ch(j)=n_ch_all*n_ch(j)

            length_ch(j) = length_ch(j)*1000./0.3048 !km to ft
            width_ch(j) = width_ch(j)/0.3048 !m to ft
            Qr_ch(j) = Qr_ch(j)/(0.3048**3) !cms to cfs

            !determine MC constant parameter method coefficients
            C1 = (1.486*slope_ch(j)**0.5)/n_ch(j)
            y = (Qr_ch(j)/(C1*width_ch(j)))**0.6
            Ax = width_ch(j)*y
            celert = (5./3.)*Qr_ch(j)/Ax
            sreach = (length_ch(j)/ndx)
            chv1(j) = width_ch(j)*C1 !used to solve for depth in MC model
            tv = dtis/sreach
            c = celert*tv
            d = Qr_ch(j)/width_ch(j)/(slope_ch(j)*celert*sreach)
            cdenom = 1 + c + d
            cc1(j) = (-1 + c + d)/cdenom
            cc2(j) = (1 + c - d)/cdenom
            cc3(j) = (1 - c + d)/cdenom
            cc4(j) = 2.*c/cdenom
    enddo
    close(4)

        return
end subroutine inputdata2

end module inputdata
