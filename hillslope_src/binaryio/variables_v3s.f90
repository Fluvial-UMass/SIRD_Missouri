module variables
        !misc file names
        character(120)::    fname10,fname11,fname12,fname13,fname15
        character(120)::    fname100
        character(120)::    fname20,fname21,fname22,fname23,fname24,fname25
        character(256)::    roffFile, restartFile, updateMode, rootDir, outDir
        character(256)::    planesFile, channelsFile

        integer,parameter :: JPRM = SELECTED_REAL_KIND(6,37)

        !output parameters
        integer,allocatable    :: idout(:)      !model unit ID for output (1 to pfafunits)
        integer,allocatable    :: idgauge(:)    !Gauge ID for output PFAF ID

        !Model unit information
        integer,allocatable    :: id(:)         !model unit ID (1 to pfafunits)
        integer,allocatable    :: pfaf(:)       !model unit pfaf ID (9 to 1)
        integer,allocatable    :: downpfaf(:)   !model unit downstream
        integer,allocatable    :: nup(:)        !number of upstream model units
        integer,allocatable    :: uppfaf(:,:)   !model unit ID for 1 or 2 unpstream units
        real(kind=JPRM),allocatable       :: A(:)          !model unit area (km2)
        real(kind=JPRM),allocatable       :: Aup(:)        !total upstream pfaf area (km2)
        real(kind=JPRM),allocatable       :: chv1(:)       !part of mannings equ, for getting depth from Q

        !WBM parameters
	    real(kind=JPRM),allocatable    :: ex_s(:)          !surface excess from WBM (ft/s)

        !channel parameters
	    real(kind=JPRM),allocatable    :: Qr_ch(:)         !reference discharge (cfs)
	    real(kind=JPRM),allocatable    :: n_ch(:)          !Manning's n of the channel
	    real(kind=JPRM),allocatable    :: length_ch(:)     !length of the channel, ft
	    !real,allocatable    :: length_p(:)     !length of the plane, ft
	    real(kind=JPRM),allocatable    :: slope_ch(:)      !slope of the channel, ft/ft
	    real(kind=JPRM),allocatable    :: width_ch(:)      !width of the channel, ft

	    real(kind=JPRM),allocatable    :: cc1(:)           !MC - Constant Parameter, ch: 1
	    real(kind=JPRM),allocatable    :: cc2(:)           !MC - Constant Parameter, ch: 2
	    real(kind=JPRM),allocatable    :: cc3(:)           !MC - Constant Parameter, ch: 3
	    real(kind=JPRM),allocatable    :: cc4(:)           !MC - Constant Parameter, ch: 4

	    real(kind=JPRM),allocatable    :: q_out(:)         !model unit ch Q (cfs) out
	    real(kind=JPRM),allocatable    :: q_in(:)          !model unit ch Q (cfs) in
	    real(kind=JPRM),allocatable    :: q_in_old(:)      !model unit ch Q (cfs) old in
	    real(kind=JPRM),allocatable    :: q_out_old(:)     !model unit ch Q (cfs) old out
	    real(kind=JPRM),allocatable    :: old_q(:,:)       !model unit ch Q (cfs) old for each dx
	    real(kind=JPRM),allocatable    :: old_q_day(:)     !model unit ch Q (cfs) for daily average
      ! real(kind=JPRM),allocatable    :: old_q_day2(:)     !model unit ch Q (cfs) for daily average
	    real(kind=JPRM),allocatable    :: old_q_ch_in(:,:)    !model unit ch Q (cfs) old for each dx
	    real(kind=JPRM),allocatable    :: old_q_ch_out(:,:)    !model unit ch Q (cfs) old for each dx
	    real(kind=JPRM),allocatable    :: old_q_res(:,:)   !model unit ch Q (cfs) old for each dx at last time step
	    real(kind=JPRM),allocatable    :: old_q_ch_in_res(:,:) !model unit ch Q (cfs) old for each dx at last time step
	    real(kind=JPRM),allocatable    :: old_q_ch_out_res(:,:) !model unit ch Q (cfs) old for each dx at last time step
	    real(kind=JPRM),allocatable    :: qlat_ch_old(:)   !old lateral flow into channel (cfs/ft)

      !Plane parameters
  	  real(kind=JPRM),allocatable    :: length_p(:)      !down slope plane length (ft)
  	  real(kind=JPRM),allocatable    :: slope_p(:)       !down slope plane  ft/ft)
  	  real(kind=JPRM),allocatable    :: ksurf(:)         !surface roughness

  	  real(kind=JPRM),allocatable    :: alp_pl_s(:)      !surface runoff alpha; y=alha*q^beta
  	  real(kind=JPRM),allocatable    :: bet_pl_s(:)      !surface runoff beat; y=alha*q^beta

  	  real(kind=JPRM),allocatable    :: y_pl_s(:,:)      !surface runoff depth (ft)
  	  real(kind=JPRM),allocatable    :: q_pl_s(:,:)      !surface runoff (cfs/ft of channel)

  	  real(kind=JPRM),allocatable    :: qlat_s_old(:)    !old lateral flow into surface runoff (ft/s)

        !Misc. parameters
        integer j, t, k, tt, ii, krec     !counters
        integer numout          !total number of model units to output data
        integer pfafunits       !total number of model units
        integer ndx             !number of dx steps
        integer ndt             !number of dt steps
        integer imax            !max number of iterations in the surface_kwave routing
        integer dtis            !seconds in a dt step
        integer dvalint         !blank value used to read data tables

        real dval               !blank value used to read data tables
        integer iyear, imonth, iday, daynum, jday !counters
        real    leapyrft,leapyrdiff

        !Channel Paramters
        real(kind=JPRM) Qr_ref                 !Q reference split for q_bank between ch and fp
        real(kind=JPRM) C1, y, Ax, celert, sreach, tv, c, d, cdenom   !MC channel parameters
        real(kind=JPRM) qlat_ch, qlat_ch_ave   !combined lateral inflow (qs + qss from both planes) to a channel (cfs/ft)

        !plane parameters
        real(kind=JPRM) eps                !min error for accepting surface_kwave routing solution

        !Uniform adjustment factors
        real(kind=JPRM) n_ch_all          ! Manning's n, uniform value for all channels
        real(kind=JPRM) ks_all            ! surface runoughness, uniform value for all planes
        real(kind=JPRM) qlat_all          ! parameter  adjusting qlat
        real(kind=JPRM) n_ch_min, n_ch_max
        real(kind=JPRM) Lch_min_slope, Lch_max_slope
        real(kind=JPRM) setfsub_rate

contains

subroutine allocatenow
        allocate        (idout(1:pfafunits))
        allocate        (idgauge(1:pfafunits))

        allocate        (id(1:pfafunits))
        allocate        (pfaf(1:pfafunits))
        allocate        (downpfaf(1:pfafunits))
        allocate        (nup(1:pfafunits))
        allocate        (uppfaf(1:pfafunits,1:4))

        allocate        (A(1:pfafunits))
        allocate        (Aup(1:pfafunits))
        allocate        (chv1(1:pfafunits))

        allocate        (ex_s(1:pfafunits))

        allocate        (slope_p(1:pfafunits))
    	allocate        (ksurf(1:pfafunits))
        allocate        (bet_pl_s(1:pfafunits))
        allocate        (alp_pl_s(1:pfafunits))
        allocate        (y_pl_s(1:pfafunits,1:ndx))
        allocate        (q_pl_s(1:pfafunits,1:ndx))
        allocate        (qlat_s_old(1:pfafunits))

        allocate        (Qr_ch(1:pfafunits))
        allocate        (length_ch(1:pfafunits))
        allocate        (length_p(1:pfafunits))
        allocate        (slope_ch(1:pfafunits))
        allocate        (width_ch(1:pfafunits))
        allocate        (n_ch(1:pfafunits))
        allocate        (cc1(1:pfafunits))
        allocate        (cc2(1:pfafunits))
        allocate        (cc3(1:pfafunits))
        allocate        (cc4(1:pfafunits))

        allocate        (q_out(1:pfafunits))
        allocate        (q_in(1:pfafunits))
        allocate        (q_in_old(1:pfafunits))
        allocate        (q_out_old(1:pfafunits))
        allocate        (old_q(1:pfafunits,1:ndx))
        allocate        (old_q_res(1:pfafunits,1:ndx))
        allocate        (old_q_day(1:pfafunits))
        ! allocate        (old_q_day2(1:pfafunits))
        allocate        (old_q_ch_in(1:pfafunits,1:ndx))
        allocate        (old_q_ch_in_res(1:pfafunits,1:ndx))
        allocate        (old_q_ch_out(1:pfafunits,1:ndx))
        allocate        (old_q_ch_out_res(1:pfafunits,1:ndx))
        allocate        (qlat_ch_old(1:pfafunits))

        return
end subroutine allocatenow

subroutine setzero()
        tt=1.
        qlat_ch = 0.
        qlat_ch_ave = 0. !added 11/19/09
        krec =1
        imax = 100
        eps = 0.001

        do j = 1,pfafunits

                idout(j) = 0.
                idgauge(j) = 0.

                id(j) = 0.
                pfaf(j) = 0.
                downpfaf(j) = 0.
                nup(j) = 0.
                uppfaf(j,1) = 0.
                uppfaf(j,2) = 0.
                uppfaf(j,3) = 0.
                uppfaf(j,4) = 0.

                slope_p(j) = 0.
                alp_pl_s(j) = 0.
                bet_pl_s(j) = 0.
        		qlat_s_old(j) = 0.
        		ksurf(j) = 0.

                Qr_ch(j) = 0.
                length_ch(j) = 0.
                length_p(j) = 0.
                slope_ch(j) = 0.
                width_ch(j) = 0.
                n_ch(j) = 0.
                chv1(j) = 0.

                cc1(j) = 0.
                cc2(j) = 0.
                cc3(j) = 0.
                cc4(j) = 0.

                ex_s(j) = 0.

                q_in(j) = 0.
                q_out(j) = 0.
                old_q_day(j)=0.
                ! old_q_day2(j)=0.

                do k = 1,ndx
                    y_pl_s(j,k)=0.
                    q_pl_s(j,k)=0.
                    old_q(j,k)=0.
                    old_q_ch_in(j,k)=0.
                    old_q_ch_out(j,k)=0.
                enddo
        enddo

        return
end subroutine setzero

subroutine deallocatenow()
        deallocate      (idout)
        deallocate      (idgauge)

        deallocate      (id)
        deallocate      (pfaf)
        deallocate      (downpfaf)
        deallocate      (nup)
        deallocate      (uppfaf)

        deallocate      (A)
        deallocate      (Aup)
        deallocate      (chv1)

        deallocate      (ex_s)

        deallocate      (ksurf)
        deallocate      (slope_p)
        deallocate      (bet_pl_s)
        deallocate      (alp_pl_s)

    	deallocate      (y_pl_s)
        deallocate      (q_pl_s)

        deallocate      (Qr_ch)
        deallocate      (length_ch)
        deallocate      (length_p)
        deallocate      (slope_ch)
        deallocate      (width_ch)
        deallocate      (n_ch)
        deallocate      (cc1)
        deallocate      (cc2)
        deallocate      (cc3)
        deallocate      (cc4)

        deallocate      (q_out)
        deallocate      (q_in)
        deallocate      (q_in_old)
        deallocate      (q_out_old)
        deallocate      (old_q)
        deallocate      (old_q_day)
        ! deallocate      (old_q_day2)
        deallocate      (old_q_ch_in)
        deallocate      (old_q_ch_out)

    return
end subroutine deallocatenow

end module variables
