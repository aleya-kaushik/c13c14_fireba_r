!---------------------------------------------------------------------
subroutine output_init( )
!---------------------------------------------------------------------

!Initializes output information

use module_pparams, only:  &
    secs_per_day, days_per_month
use module_sibconst, only: &
     subcount, nlu, npoolpft, nsoil, spinup, nsib, subset
use module_io
use module_sib, only: sib
use module_time

use kinds
use netcdf
implicit none

!...local variables
integer(i4) :: xmon
integer(i4) :: daysref, daystot, monref, yearref
integer(i4) :: ncid, status, vid, schema, dtsib_saved, count_saved
integer(i4) :: ng_saved, nlu_saved, nsib_saved, di, j, global_ref, l
integer(i4) :: held_schema, held_len
integer(i4), parameter :: nheld_diag=10
real(r8), allocatable :: held_scalars(:,:,:), held_gain(:,:,:,:), &
     held_loss(:,:,:,:)
real(r8), allocatable :: global_alloc(:,:,:), global_alloc_moist(:,:,:), &
     global_alloc_temp(:,:,:)
real(r8), allocatable :: global_grid(:,:), global_lu(:,:,:)
integer(i4), allocatable :: global_sif(:,:)
integer(i4), allocatable :: refcheck(:)
logical :: exists_g, exists_lu

    ! Last active allocation fractions are retained across exact restarts.
    ! Restore after sibtype_setup(), which otherwise overwrites this held
    ! state with current phenological allocation.  The saved fractions only
    ! affect output while dormant; active allocation is recalculated normally.
    if (ic_restart) then
       status = nf90_open(trim(ic_file),nf90_nowrite,ncid)
       if (status /= nf90_noerr) &
          error stop 'Allocation continuation: cannot open restart file'
       status = nf90_inq_dimid(ncid,'nsib',di)
       if (status /= nf90_noerr) &
          error stop 'Allocation continuation: restart has no spatial dimension'
       status = nf90_inquire_dimension(ncid,di,len=nsib_saved)
       if (status /= nf90_noerr) &
          error stop 'Allocation continuation: cannot read spatial dimension'
       if (nsib_saved /= subcount .and. nsib_saved /= nsib) &
          error stop 'Allocation continuation: incompatible spatial dimension'
       status = nf90_inq_varid(ncid,'restart_alloc_current',vid)
       if (status /= nf90_noerr) &
          error stop 'Allocation continuation: checkpoint lacks allocation state'
       allocate(global_alloc(nsib_saved,nlu,npoolpft))
       allocate(global_alloc_moist(nsib_saved,nlu,npoolpft))
       allocate(global_alloc_temp(nsib_saved,nlu,npoolpft))
       status = nf90_get_var(ncid,vid,global_alloc)
       if (status /= nf90_noerr) &
          error stop 'Allocation continuation: cannot read allocation state'
       status = nf90_inq_varid(ncid,'restart_alloc_moist',vid)
       if (status /= nf90_noerr) &
          error stop 'Allocation continuation: checkpoint lacks moist allocation'
       status = nf90_get_var(ncid,vid,global_alloc_moist)
       if (status /= nf90_noerr) &
          error stop 'Allocation continuation: cannot read moist allocation'
       status = nf90_inq_varid(ncid,'restart_alloc_temp',vid)
       if (status /= nf90_noerr) &
          error stop 'Allocation continuation: checkpoint lacks thermal allocation'
       status = nf90_get_var(ncid,vid,global_alloc_temp)
       if (status /= nf90_noerr) &
          error stop 'Allocation continuation: cannot read thermal allocation'
       do j=1,subcount
          global_ref=j
          if (nsib_saved == nsib) global_ref=subset(j)
          do l=1,sib%g(j)%g_nlu
             sib%g(j)%l(l)%poollt%alloc(:) = global_alloc(global_ref,l,:)
             sib%g(j)%l(l)%poollt%alloc_moist(:) = global_alloc_moist(global_ref,l,:)
             sib%g(j)%l(l)%poollt%alloc_temp(:) = global_alloc_temp(global_ref,l,:)
          enddo
       enddo
       deallocate(global_alloc,global_alloc_moist,global_alloc_temp)
       status = nf90_close(ncid)
       if (status /= nf90_noerr) &
          error stop 'Allocation continuation: cannot close checkpoint'
    endif

    !------------------------------------------------------------
    ! Rehydrate held diagnostic values after sibtype_setup(), which
    ! otherwise initializes them to zero.  Exact restart only: do not
    ! change any scientific update conditions or fresh requib starts.
    ! Processor-local checkpoint dimensions may be subcount, or nsib
    ! in a merged checkpoint.  The same spatial mapping as allocation
    ! continuation is applied to every held field.
    !------------------------------------------------------------
    if (ic_restart) then
       status=nf90_open(trim(ic_file),nf90_nowrite,ncid)
       if (status /= nf90_noerr) &
          error stop 'Held diagnostic continuation: cannot open checkpoint'
       held_schema=-1
       status=nf90_get_att(ncid,nf90_global,'held_diag_schema',held_schema)
       if (status /= nf90_noerr .or. held_schema /= 1) &
          error stop 'Held diagnostic continuation: unsupported checkpoint schema'
       status=nf90_inq_dimid(ncid,'nsib',di)
       if (status /= nf90_noerr) &
          error stop 'Held diagnostic continuation: missing nsib'
       status=nf90_inquire_dimension(ncid,di,len=nsib_saved)
       if (status /= nf90_noerr) &
          error stop 'Held diagnostic continuation: invalid nsib'
       if (nsib_saved /= subcount .and. nsib_saved /= nsib) &
          error stop 'Held diagnostic continuation: incompatible site count'
       status=nf90_inq_dimid(ncid,'held_diag_count',di)
       if (status /= nf90_noerr) &
          error stop 'Held diagnostic continuation: no scalar dimension'
       status=nf90_inquire_dimension(ncid,di,len=held_len)
       if (status /= nf90_noerr .or. held_len /= nheld_diag) &
          error stop 'Held diagnostic continuation: scalar count mismatch'
       status=nf90_inq_dimid(ncid,'nlu',di)
       if (status /= nf90_noerr) &
          error stop 'Held diagnostic continuation: missing nlu dimension'
       status=nf90_inquire_dimension(ncid,di,len=held_len)
       if (status /= nf90_noerr .or. held_len /= nlu) &
          error stop 'Held diagnostic continuation: incompatible nlu'
       status=nf90_inq_dimid(ncid,'npoolpft',di)
       if (status /= nf90_noerr) &
          error stop 'Held diagnostic continuation: missing npoolpft dimension'
       status=nf90_inquire_dimension(ncid,di,len=held_len)
       if (status /= nf90_noerr .or. held_len /= npoolpft) &
          error stop 'Held diagnostic continuation: incompatible npoolpft'
       status=nf90_inq_dimid(ncid,'nsoil',di)
       if (status /= nf90_noerr) &
          error stop 'Held diagnostic continuation: missing nsoil dimension'
       status=nf90_inquire_dimension(ncid,di,len=held_len)
       if (status /= nf90_noerr .or. held_len /= nsoil) &
          error stop 'Held diagnostic continuation: incompatible nsoil'
       allocate(held_scalars(nsib_saved,nlu,nheld_diag))
       allocate(held_gain(nsib_saved,nlu,npoolpft,nsoil))
       allocate(held_loss(nsib_saved,nlu,npoolpft,nsoil))
       status=nf90_inq_varid(ncid,'held_diag_scalars',vid)
       if (status /= nf90_noerr) &
          error stop 'Held diagnostic continuation: missing scalar values'
       status=nf90_get_var(ncid,vid,held_scalars)
       if (status /= nf90_noerr) &
          error stop 'Held diagnostic continuation: scalar read failed'
       status=nf90_inq_varid(ncid,'held_assim_dgain',vid)
       if (status /= nf90_noerr) &
          error stop 'Held diagnostic continuation: missing gain values'
       status=nf90_get_var(ncid,vid,held_gain)
       if (status /= nf90_noerr) &
          error stop 'Held diagnostic continuation: gain read failed'
       status=nf90_inq_varid(ncid,'held_assim_dloss',vid)
       if (status /= nf90_noerr) &
          error stop 'Held diagnostic continuation: missing loss values'
       status=nf90_get_var(ncid,vid,held_loss)
       if (status /= nf90_noerr) &
          error stop 'Held diagnostic continuation: loss read failed'
       do j=1,subcount
          global_ref=j
          if (nsib_saved == nsib) global_ref=subset(j)
          do l=1,sib%g(j)%g_nlu
             sib%g(j)%l(l)%fract%d13cassim=held_scalars(global_ref,l,1)
             sib%g(j)%l(l)%fract%d13cassim_nog=held_scalars(global_ref,l,2)
             sib%g(j)%l(l)%fract%d13cresp_growc13=held_scalars(global_ref,l,3)
             sib%g(j)%l(l)%fract%kiecps_k1=held_scalars(global_ref,l,4)
             sib%g(j)%l(l)%fract%kiecps_k2=held_scalars(global_ref,l,5)
             sib%g(j)%l(l)%fract%kiecps_k3=held_scalars(global_ref,l,6)
             sib%g(j)%l(l)%fract%kiecps_k4=held_scalars(global_ref,l,7)
             sib%g(j)%l(l)%fract%kiecps_k5=held_scalars(global_ref,l,8)
             sib%g(j)%l(l)%fract%kiecps_nog=held_scalars(global_ref,l,9)
             sib%g(j)%l(l)%sift%sif_x=held_scalars(global_ref,l,10)
             sib%g(j)%l(l)%poollt%assim_dgain(:,:)=held_gain(global_ref,l,:,:)
             sib%g(j)%l(l)%poollt%assim_dloss(:,:)=held_loss(global_ref,l,:,:)
          enddo
       enddo
       deallocate(held_scalars,held_gain,held_loss)
       status=nf90_close(ncid)
       if (status /= nf90_noerr) &
          error stop 'Held diagnostic continuation: checkpoint close failed'
    endif

    !-----------------------------------------
    !...Initialize output variables
    IF (restart_savef .OR. (write_final_restart .AND. (.NOT. spinup)) .OR. &
         hr_savegf .OR. hr_saveluf .OR. &
         pbp_savegf .OR. pbp_saveluf .OR. &
         qp_savegf .OR. qp_saveluf) THEN
        print *, 'Setting Output'
    ENDIF

    !-----------------------------------------
    !...Restart Files
    !-----------------------------------------
    if (restart_savef) then
       if (restart_dtsib > 0) then
          if (mod(restart_dtsib, 86400) .ne. 0.) then
             !...Hourly (not allowed since SiB4 starts at 0 GMT)
             print('(a)'),'Restart Output Must Be Days or Months.  Stopping.'
             STOP
           else
              !..Daily
              daystot = restart_dtsib/86400.
              print('(a,i6)'), '    SiB Restart Files Written (days) = ', daystot
              restart_ntpersave = daystot * steps_per_day
              restart_countt = 1

              ! Positive-day restart filenames are assigned at the actual write
              ! time in output_control(), using the post-time_manager year/doy.
              ! This keeps the filename synchronized with the checkpoint state.
              restart_filename = ''
           endif

       else !monthly
           print('(a,i4)'),'    SiB restart files written (month) = ', -restart_dtsib

           daystot = 0
           monref = month-1
           yearref = year
           do xmon=1,abs(restart_dtsib)
              monref = monref + 1
              if (monref .GT. 12) then
                 monref = 1
                 yearref = yearref + 1
              endif

              daysref = days_per_month(monref)
              if ((monref .eq. 2) .and. &
                  (((mod(endyear,4) .eq. 0) .and. &
                   (mod(endyear,100) .ne. 0)) .or. &
                   (mod(endyear,400) .eq. 0))) then
                  daysref = daysref + 1
              endif
              daystot = daystot + daysref 
           enddo

           restart_ntpersave = daystot * steps_per_day
           restart_countt = 1

           !...update mon/year ref for output file name
           monref = monref+1
           if (monref .gt. 12) then
               monref=1
               yearref=yearref+1
           endif

           if (printrank) then
              write( restart_filename, '(a,i4.4,i2.2,a,i5.5,a)') trim(out_path)//'sib_r', &
                     yearref,monref,'p',rank,'.nc'
           else
              write( restart_filename, '(a,i4.4,i2.2,a)') trim(out_path)//'sib_r', &
                     yearref,monref,'.nc'
           endif
       endif !restart_dtsib
    endif !restart_savef

    if (write_final_restart .and. (.not. spinup)) then
       print('(a)'),'    SiB final transient restart enabled'
    endif

    !-----------------------------------------
    !...HR Diagnostic Output (Typically Hourly)
    !-----------------------------------------
    if ((hr_savegf) .OR. (hr_saveluf)) then
       if (hr_dtsib > secs_per_day) then
           !...Multi-Day HR Output
           print('(a)'),'Not Expecting Multi-Day HR Output.'
           print('(a)'),'Stopping.'
           stop
       elseif (hr_dtsib > izero) then
           !...Daily and Hourly HR Output
           print('(a,i6)'), '    SiB HR files written (s) = ', hr_dtsib
         
           hr_start = doy-1
           hr_step = hr_dtsib / dble(secs_per_day)
           hr_ntpersave = hr_dtsib / dtsib
           hr_nsaveperout = dtoutwrite / hr_dtsib
           if (hr_nsaveperout <= 0 .or. mod(secs_per_day,hr_dtsib) /= 0) &
              error stop 'HR continuation: unsupported output interval'
           monref = month
           yearref = year
           call diagnostic_filename('hr',yearref,monref, &
                hr_filenameg,hr_filenamelu)
           exists_g = .false.
           exists_lu = .false.
           if (ic_restart .and. day > 1) then
              if (hr_savegf) inquire(file=trim(hr_filenameg),exist=exists_g)
              if (hr_saveluf) inquire(file=trim(hr_filenamelu),exist=exists_lu)
              if ((hr_savegf .and. hr_saveluf) .and. &
                  (exists_g .neqv. exists_lu)) &
                  error stop 'HR continuation: partial-month grid/LU files disagree'
           endif
           if (exists_g .or. exists_lu) then
              hr_start = doy-day
              hr_nsaveperfile = curday_per_mon*hr_nsaveperout
           else
              daystot = curday_per_mon-day+1
              if ((year .eq. endyear) .and. &
                  ((doy+daystot-1) .gt. endtime)) &
                  daystot = endtime-doy+1
              hr_nsaveperfile = daystot*hr_nsaveperout
           endif

       elseif (hr_dtsib > -13) then
           !...Monthly HR Output
           print('(a,i4)'), '    SiB HR files written (mon) = ', abs(hr_dtsib)

            daystot = izero
            monref = month-1
            yearref = year
            do xmon=1,abs(hr_dtsib)
               monref = monref + 1
               if (monref .GT. 12) then
                  monref = 1
                  yearref = yearref + 1
               endif

               daystot = daystot + curday_per_mon
           enddo

           hr_start = (doy-1) + daystot/dble(2.)
           if (hr_start .gt. curday_per_year) &
               hr_start = hr_start - curday_per_year

           hr_step = 0
           hr_ntpersave = daystot*steps_per_day
           hr_nsaveperout = 1
           hr_nsaveperfile = 1

       else
           !...Yearly HR Output
           print('(a)'),'Not Expecting Yearly HR Output.'
           print('(a)'),'Stopping.'
           stop
       endif

       !...General HR counters/holders
       hr_countt = 1
       hr_countsave = 1
       hr_wtpersave = 1. / real(hr_ntpersave)

       call diagnostic_filename('hr', yearref, monref, &
            hr_filenameg, hr_filenamelu)

        IF (hr_savegf) THEN
           allocate( hr_g(subcount,hr_nvarg,hr_nsaveperout))
           allocate( hr_idvarg(hr_nvarg))
           hr_g = dzero
           hr_idvarg = izero
        ENDIF
        IF (hr_saveluf) THEN
           allocate( hr_lu(subcount,nlu,hr_nvarlu,hr_nsaveperout))
           allocate( hr_idvarlu(hr_nvarlu))
           hr_lu = dzero
           allocate(hr_lu_pref_time(subcount,nlu,hr_nsaveperout))
           allocate(hr_lu_area_time(subcount,nlu,hr_nsaveperout))
           hr_lu_pref_time = izero
           hr_lu_area_time = rzero
           hr_idvarlu = izero
        ENDIF
        
        !...General HR file info
        hr_openf = .false.
        hr_idg = izero
        hr_idlu = izero

        hr_writef = .false.
        hr_countwrite = 1
        if (ic_restart .and. day > 1 .and. hr_dtsib > 0 .and. &
            (exists_g .or. exists_lu)) then
           hr_countwrite = (day-1)*hr_nsaveperout + 1
           hr_openf = .true.
        endif

        !...SIF HR file info
        allocate(hr_sif_satcount(subcount,hr_nsaveperout,2))
        hr_sif_satcount = izero
    
    endif  !hrintg or hrintlu

    !-----------------------------------------
    !...QP Diagnostic Output (Typically Monthly)
    !-----------------------------------------
    if ((qp_savegf) .OR. (qp_saveluf)) then
       if (qp_dtsib > secs_per_day) then
           !...Multi-Day QP Output
           print('(a)'),'Not Expecting Multi-Day QP Output.'
           print('(a)'),'Stopping.'
           stop
       elseif (qp_dtsib > izero) then
           !...Daily and Hourly QP Output
           print('(a,i6)'), '    SiB QP files written (s) = ', qp_dtsib
         
           qp_start = doy-1
           qp_step = qp_dtsib / dble(secs_per_day)
           qp_ntpersave = qp_dtsib / dtsib
           qp_nsaveperout = dtoutwrite / qp_dtsib
           qp_nsaveperfile = curday_per_mon*qp_nsaveperout
           if ((year .eq. endyear) .and. &
               ((doy+curday_per_mon) .gt. endtime)) then
               qp_nsaveperfile = (endtime - doy + 1)*qp_nsaveperout
           endif
           monref = month
           yearref = year

       elseif (qp_dtsib > -13) then
           !...Monthly QP Output
           print('(a,i4)'), '    SiB QP files written (mon) = ', abs(qp_dtsib)

            daystot = izero
            monref = month-1
            yearref = year
            do xmon=1,abs(qp_dtsib)
               monref = monref + 1
               if (monref .GT. 12) then
                  monref = 1
                  yearref = yearref + 1
               endif

               daystot = daystot + curday_per_mon
           enddo

           ! If continuing in the middle of a month, the first QP
           ! average still represents the FULL calendar month.  Its
           ! unfinished raw accumulation is recovered below from ic_file.
           if (qp_dtsib == -1 .and. day > 1) then
              daystot = curday_per_mon
              qp_start = dble(doy-day) + dble(curday_per_mon)/2.d0
           else
              qp_start = (doy-1) + daystot/dble(2.)
           endif
           if (qp_start .gt. curday_per_year) &
               qp_start = qp_start - curday_per_year

           qp_step = 0
           qp_ntpersave = daystot*steps_per_day
           qp_nsaveperout = 1
           qp_nsaveperfile = 1

       else
           !...Yearly QP Output
           print('(a)'),'Not Expecting Yearly QP Output.'
           print('(a)'),'Stopping.'
           stop
       endif

       !...General QP counters/holders
       qp_countt = 1
       qp_countsave = 1
       qp_wtpersave = 1. / real(qp_ntpersave)

       call diagnostic_filename('qp', yearref, monref, &
            qp_filenameg, qp_filenamelu)

        IF (qp_savegf) THEN
           allocate( qp_g(subcount,qp_nvarg,qp_nsaveperout))
           allocate( qp_idvarg(qp_nvarg))
           qp_g = dzero
           qp_idvarg = izero
        ENDIF
        IF (qp_saveluf) THEN
           allocate( qp_lu(subcount,nlu,qp_nvarlu,qp_nsaveperout))
           allocate( qp_idvarlu(qp_nvarlu))
           qp_lu = dzero
           allocate(qp_lu_pref_time(subcount,nlu,qp_nsaveperout))
           allocate(qp_lu_area_time(subcount,nlu,qp_nsaveperout))
           qp_lu_pref_time = izero
           qp_lu_area_time = rzero
           qp_idvarlu = izero
        ENDIF

        !...General QP file info
        qp_openf = .false.
        qp_idg = izero
        qp_idlu = izero

        qp_writef = .false.
        qp_countwrite = 1

        !...SIF QP file info
        allocate(qp_sif_satcount(subcount,qp_nsaveperout,2))
        qp_sif_satcount = izero

        ! Exact restart, midmonth: recover raw QP sums/counters.  Do not
        ! reconstruct month totals from rounded daily diagnostics.  The
        ! checkpoint records the POST-timestep QP counter (output_control
        ! writes the checkpoint after advancing diagnostic counters).
        if (ic_restart .and. day > 1 .and. qp_dtsib == -1) then
           status = nf90_open(trim(ic_file),nf90_nowrite,ncid)
           if (status /= nf90_noerr) &
              error stop 'QP continuation: cannot open exact restart file'
           status = nf90_get_att(ncid,nf90_global,'qp_continuation_schema',schema)
           if (status /= nf90_noerr) &
              error stop 'QP continuation: checkpoint lacks partial-month QP state; regenerate checkpoint'
           if (schema /= 1) error stop 'QP continuation: unsupported schema'
           status = nf90_get_att(ncid,nf90_global,'qp_continuation_dtsib',dtsib_saved)
           if (status /= nf90_noerr .or. dtsib_saved /= qp_dtsib) &
              error stop 'QP continuation: output cadence differs from checkpoint'
           status = nf90_get_att(ncid,nf90_global,'qp_continuation_countt',count_saved)
           if (status /= nf90_noerr) error stop 'QP continuation: missing QP counter'
           if (count_saved /= (day-1)*steps_per_day+1) &
              error stop 'QP continuation: checkpoint date and QP counter disagree'
           status = nf90_get_att(ncid,nf90_global,'qp_continuation_nvarg',ng_saved)
           if (status /= nf90_noerr .or. ng_saved /= qp_nvarg) &
              error stop 'QP continuation: grid variable count changed'
           status = nf90_get_att(ncid,nf90_global,'qp_continuation_nvarlu',nlu_saved)
           if (status /= nf90_noerr .or. nlu_saved /= qp_nvarlu) &
              error stop 'QP continuation: LU variable count changed'
           status = nf90_inq_dimid(ncid,'nsib',di)
           if (status /= nf90_noerr) error stop 'QP continuation: missing nsib'
           status = nf90_inquire_dimension(ncid,di,len=nsib_saved)
           if (status /= nf90_noerr) error stop 'QP continuation: cannot read nsib'
           if (nsib_saved /= subcount .and. nsib_saved /= nsib) &
              error stop 'QP continuation: incompatible spatial dimension'
           allocate(global_sif(nsib_saved,2))
           status = nf90_inq_varid(ncid,'qp_continuation_sif_counts',vid)
           if (status /= nf90_noerr) error stop 'QP continuation: no SIF counts'
           status = nf90_get_var(ncid,vid,global_sif)
           if (status /= nf90_noerr) error stop 'QP continuation: cannot read SIF counts'
           do j=1,subcount
              global_ref=j
              if (nsib_saved == nsib) global_ref=subset(j)
              qp_sif_satcount(j,1,:)=global_sif(global_ref,:)
           enddo
           deallocate(global_sif)
           if (qp_savegf) then
              status = nf90_inq_varid(ncid,'qp_continuation_g_ids',vid)
              if (status /= nf90_noerr) error stop 'QP continuation: no grid ID inventory'
              allocate(refcheck(qp_nvarg))
              status = nf90_get_var(ncid,vid,refcheck)
              if (status /= nf90_noerr) error stop 'QP continuation: cannot read grid IDs'
              if (any(refcheck /= qp_vrefg)) &
                 error stop 'QP continuation: grid outopts changed across restart'
              deallocate(refcheck)
              status = nf90_inq_varid(ncid,'qp_continuation_g',vid)
              if (status /= nf90_noerr) error stop 'QP continuation: no grid sums'
              allocate(global_grid(nsib_saved,qp_nvarg))
              status = nf90_get_var(ncid,vid,global_grid)
              if (status /= nf90_noerr) error stop 'QP continuation: cannot read grid sums'
              do j=1,subcount
                 global_ref=j
                 if (nsib_saved == nsib) global_ref=subset(j)
                 qp_g(j,:,1)=global_grid(global_ref,:)
              enddo
              deallocate(global_grid)
           endif
           if (qp_saveluf) then
              status = nf90_inq_varid(ncid,'qp_continuation_lu_ids',vid)
              if (status /= nf90_noerr) error stop 'QP continuation: no LU ID inventory'
              allocate(refcheck(qp_nvarlu))
              status = nf90_get_var(ncid,vid,refcheck)
              if (status /= nf90_noerr) error stop 'QP continuation: cannot read LU IDs'
              if (any(refcheck /= qp_vreflu)) &
                 error stop 'QP continuation: LU outopts changed across restart'
              deallocate(refcheck)
              status = nf90_inq_varid(ncid,'qp_continuation_lu',vid)
              if (status /= nf90_noerr) error stop 'QP continuation: no LU sums'
              allocate(global_lu(nsib_saved,nlu,qp_nvarlu))
              status = nf90_get_var(ncid,vid,global_lu)
              if (status /= nf90_noerr) error stop 'QP continuation: cannot read LU sums'
              do j=1,subcount
                 global_ref=j
                 if (nsib_saved == nsib) global_ref=subset(j)
                 qp_lu(j,:,:,1)=global_lu(global_ref,:,:)
              enddo
              deallocate(global_lu)
           endif
           status = nf90_close(ncid)
           if (status /= nf90_noerr) error stop 'QP continuation: closing restart failed'
           qp_countt = count_saved
        endif

    endif  !qp_savegf or qp_saveluf

    !-----------------------------------------
    !...PBP Diagnostic Output (Typically Daily)
    !-----------------------------------------
    if ((pbp_savegf) .or. (pbp_saveluf)) then

       if (pbp_dtsib > secs_per_day) then
           !...Multi-Day PBP Output
           print('(a)'),'Not Expecting Multi-Day PBP Output.'
           print('(a)'),'Stopping.'
           stop
       elseif (pbp_dtsib > izero) then
           !...Daily and Hourly PBP Output
           print('(a,i6)'), '    SiB PBP files written (s) = ', pbp_dtsib
         
           pbp_start = doy-1
           pbp_step = pbp_dtsib / dble(secs_per_day)
           pbp_ntpersave = pbp_dtsib / dtsib
           pbp_nsaveperout = dtoutwrite / pbp_dtsib
           if (pbp_nsaveperout <= 0 .or. mod(secs_per_day,pbp_dtsib) /= 0) &
              error stop 'PBP continuation: unsupported output interval'

           ! A restarted segment may share an output directory with
           ! its predecessor.  Only append when BOTH selected products
           ! have the expected pre-restart month file; the writer checks
           ! existing time records before touching either file.
           exists_g = .false.
           exists_lu = .false.
           monref = month
           yearref = year
           call diagnostic_filename('pbp',yearref,monref, &
                pbp_filenameg,pbp_filenamelu)
           if (ic_restart .and. day > 1 .and. pbp_dtsib > 0) then
              if (pbp_savegf) inquire(file=trim(pbp_filenameg),exist=exists_g)
              if (pbp_saveluf) inquire(file=trim(pbp_filenamelu),exist=exists_lu)
              if ((pbp_savegf .and. pbp_saveluf) .and. &
                  (exists_g .neqv. exists_lu)) &
                  error stop 'PBP continuation: one of grid/LU partial-month files is absent'
           endif
           if (exists_g .or. exists_lu) then
              pbp_start = doy-day
              pbp_nsaveperfile = curday_per_mon*pbp_nsaveperout
           else
              pbp_start = doy-1
              daystot = curday_per_mon-day+1
              if ((year .eq. endyear) .and. &
                  ((doy+daystot-1) .gt. endtime)) then
                 daystot = endtime-doy+1
              endif
              pbp_nsaveperfile = daystot*pbp_nsaveperout
           endif
           monref = month
           yearref = year

       elseif (pbp_dtsib > -13) then
           !...Monthly PBP Output
           print('(a,i4)'), '    SiB PBP files written (mon) = ', abs(pbp_dtsib)

            daystot = izero
            monref = month-1
            yearref = year
            do xmon=1,abs(pbp_dtsib)
               monref = monref + 1
               if (monref .GT. 12) then
                  monref = 1
                  yearref = yearref + 1
               endif

               daystot = daystot + curday_per_mon
           enddo

           pbp_start = (doy-1) + daystot/dble(2.)
           if (pbp_start .gt. curday_per_year) &
               pbp_start = pbp_start - curday_per_year

           pbp_step = 0
           pbp_ntpersave = daystot*steps_per_day
           pbp_nsaveperout = 1
           pbp_nsaveperfile = 1

       else
           !...Yearly PBP Output
           print('(a)'),'Not Expecting Yearly PBP Output.'
           print('(a)'),'Stopping.'
           stop
       endif

       !...General PBP counters/holders
       pbp_countt = 1
       pbp_countsave = 1
       pbp_wtpersave = 1. / real(pbp_ntpersave)

       call diagnostic_filename('pbp', yearref, monref, &
            pbp_filenameg, pbp_filenamelu)

        IF (pbp_savegf) THEN
           allocate( pbp_g(npbp,pbp_nvarg,pbp_nsaveperout))
           allocate(pbp_idvarg(pbp_nvarg))
           pbp_g = dzero
           pbp_idvarg = izero
        ENDIF
        IF (pbp_saveluf) THEN
           allocate( pbp_lu(npbp,nlu,pbp_nvarlu,pbp_nsaveperout))
           allocate(pbp_idvarlu(pbp_nvarlu))
            pbp_lu = dzero
           allocate(pbp_lu_pref_time(npbp,nlu,pbp_nsaveperout))
           allocate(pbp_lu_area_time(npbp,nlu,pbp_nsaveperout))
           pbp_lu_pref_time = izero
           pbp_lu_area_time = rzero
            pbp_idvarlu = izero
         ENDIF
         
        !...General PBP file info
        pbp_openf = .false.
        pbp_idg = izero
        pbp_idlu = izero

        pbp_writef = .false.
        pbp_countwrite = 1
        if (ic_restart .and. day > 1 .and. pbp_dtsib > 0 .and. &
            (exists_g .or. exists_lu)) then
           pbp_countwrite = (day-1)*pbp_nsaveperout + 1
           pbp_openf = .true.
        endif

        !...SIF PBP file info
        allocate(pbp_sif_satcount(subcount,pbp_nsaveperout,2))
        pbp_sif_satcount = izero
        
    endif  !pbpintg or pbpintlu

end subroutine output_init

!---------------------------------------------------------------------
subroutine diagnostic_filename( &
           cout_type, yr, mon, filenameg, filenamelu) 
!---------------------------------------------------------------------

use kinds
use module_io

implicit none

!...input variables
character(len=*), intent(in) :: cout_type
integer(i4), intent(in) :: yr, mon
character(len=iofnlen), intent(inout) :: &
    filenameg, filenamelu

!...PBP Files
IF (cout_type .eq. 'pbp') THEN

   if (printrank) then
       write( filenameg, '(a,i4.4,i2.2,a,i5.5,a)') &
              trim(out_path)//pbp_prefix, &
              yr,mon,'p',rank,trim(pbp_suffixg)
       write( filenamelu, '(a,i4.4,i2.2,a,i5.5,a)') &
              trim(out_path)//pbp_prefix, &
              yr,mon,'p',rank,trim(pbp_suffixlu)
    else
       write( filenameg, '(a,i4.4,i2.2,a)') &
              trim(out_path)//pbp_prefix, &
              yr,mon,trim(pbp_suffixg)
       write( filenamelu, '(a,i4.4,i2.2,a)') &
               trim(out_path)//pbp_prefix, &
               yr,mon,trim(pbp_suffixlu)
     endif

ELSEIF (cout_type .eq. 'qp') THEN

   if (printrank) then
       write( filenameg, '(a,i4.4,i2.2,a,i5.5,a)') &
              trim(out_path)//qp_prefix, &
              yr,mon,'p',rank,trim(qp_suffixg)
       write( filenamelu, '(a,i4.4,i2.2,a,i5.5,a)') &
              trim(out_path)//qp_prefix, &
              yr,mon,'p',rank,trim(qp_suffixlu)
    else
       write( filenameg, '(a,i4.4,i2.2,a)') &
              trim(out_path)//qp_prefix, &
              yr,mon,trim(qp_suffixg)
       write( filenamelu, '(a,i4.4,i2.2,a)') &
               trim(out_path)//qp_prefix, &
               yr,mon,trim(qp_suffixlu)
     endif

ELSEIF (cout_type .eq. 'hr') THEN

   if (printrank) then
       write( filenameg, '(a,i4.4,i2.2,a,i5.5,a)') &
              trim(out_path)//hr_prefix, &
              yr,mon,'p',rank,trim(hr_suffixg)
       write( filenamelu, '(a,i4.4,i2.2,a,i5.5,a)') &
              trim(out_path)//hr_prefix, &
              yr,mon,'p',rank,trim(hr_suffixlu)
    else
       write( filenameg, '(a,i4.4,i2.2,a)') &
              trim(out_path)//hr_prefix, &
              yr,mon,trim(hr_suffixg)
       write( filenamelu, '(a,i4.4,i2.2,a)') &
               trim(out_path)//hr_prefix, &
               yr,mon,trim(hr_suffixlu)
     endif

ELSE
    print*,'Unknown Output File Type:',cout_type
    print*,'Stopping.'
    stop
ENDIF


end subroutine diagnostic_filename
