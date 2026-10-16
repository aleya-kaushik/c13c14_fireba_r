
! Reads in sibdrv control variables.
subroutine read_namel( namel_name )

use module_sibconst
use module_io
use module_time

implicit none

!...input variables
character(len=350) :: namel_name

!...local variables
character(len=256) :: dr_path, fr_path, frba_path
!character(len=512) :: icfile

!...NAMELISTS
namelist /CONTROL_LIST_SIBDRV/ & ! USER SELECTED TIMES
    nsib, starttime, startyear, endtime, endyear, &
    dtsib, restart_dtsib, write_final_restart, qp_dtsib, pbp_dtsib, hr_dtsib
namelist /IO_LIST_SIBDRV/ & ! USER SELETED I/O OPTIONS
    pft_info, pool_info, &
    aero_file, pgdd_file, pstg_file, &
    phys_file, pool_file, isodata_file, firep_file, &
    vs_file, ic_file, ic_mode, minveg_diag, minveg_clim_policy, dr_path, tm5_path, fr_path, frba_path, &
    out_path, out_info, out_rinfo
namelist /SPINUP_LIST_SIBDRV/ & ! USER SELECTED SPINUP OPTIONS
    spinup, spinup_default, spinup_continue, spinup_numyrs, spinup_maxiter, &
    spinup_threshold, spinup_writediag, spinup_writetxtf
namelist /SUBGRID_SIBDRV/ & ! SUBGRID INFO FOR NON-POINT RUNS
    minlon, maxlon, minlat, maxlat
namelist /PBP_LIST_SIBDRV/ & ! USER DEFINED PBP DIAGNOSTIC LOCATIONS
    npbp
namelist /BALAN_LIST_SIBDRV/ & ! USER SELECTED BALANCE CHECKING OPTIONS
    badtc_print, badtc_stop, &
    canb_print, canb_stop, canb_thresh, &
    carbonb_print, carbonb_stop, carbonb_thresh, &
    carbonb_threshc13, carbonb_threshc14, &
    fireb_print, fireb_stop, fireb_thresh, &
    snocmbn_print, snocmbn_stop, snocmbn_thresh, &
    bnum_allow, energyb_print, energyb_stop, energyb_thresh, &
    waterb_print, waterb_stop, waterb_thresh
namelist /PRINT_LIST_SIBDRV/ & ! USER SELECTED PRINTING OPTIONS
    print_avec, print_driver, print_fire, print_harvest, &
    print_pftinfo, print_pooll, print_soil, print_sscol, &
    print_veg, print_stop
namelist /SWITCH_LIST_SIBDRV/ & ! USER SELECTED SWITCHES
    cornsoy_switch, fire_switch, firegfed_switch, grazing_switch, green_switch, &
    eqclear_switch, leapyr_switch, updatelst_switch, tm5mr_switch, &
    soilogee_switch, varco2_switch, varciso_switch, varcisom_switch, &
    pftswap_enable, pftswap_open_year, pftswap_open_doy

!print*,'icfile:',icfile
    !-----------------------------------------------------------------------
    ! read in namel_sibdrv
    !-----------------------------------------------------------------------
    print*,''
    print*,'Reading in namelist file: ',trim(namel_name)
    open(unit=2,file=trim(namel_name),form='formatted')
    read (2,CONTROL_LIST_SIBDRV)
    read (2,IO_LIST_SIBDRV)
    read (2,SPINUP_LIST_SIBDRV)
    read (2,SUBGRID_SIBDRV)
    read (2,PBP_LIST_SIBDRV)
    IF (npbp .GT. 0) THEN
       allocate (pbp_lonlat(2,npbp))
       pbp_lonlat = 0.0
       read(2,*,err=919)pbp_lonlat
       919  continue
    ELSE
       allocate(pbp_lonlat(2,nsib))
       
    ENDIF

    read (2,BALAN_LIST_SIBDRV)
    ! The C14 tolerance applies to the decay-corrected *residual*.
    ! Default scale is approximate natural-abundance 14C/13C (~1e-10).
    ! Prefer an explicit, tested carbonb_threshc14 in production namelists.
    if (carbonb_threshc14 .lt. 0.0_r4) then
       carbonb_threshc14 = max(1.0e-20_r4, carbonb_threshc13 * 1.0e-10_r4)
       print '(a,es13.5)', 'Default C14 balance threshold (mol 14C/m2): ', &
             carbonb_threshc14
    endif
    if (carbonb_threshc14 .le. 0.0_r4) then
       error stop 'carbonb_threshc14 must be positive'
    endif
    read (2,PRINT_LIST_SIBDRV)
    read (2,SWITCH_LIST_SIBDRV)

    close(2)

    ! Interpret the initial-condition mode once here so downstream startup
    ! routines do not have to infer semantics from the filename or date.
    ic_mode = adjustl(ic_mode)
    select case (trim(ic_mode))
    case ('requib')
       ic_restart = .false.
    case ('restart')
       ic_restart = .true.
    case default
       print*,'ERROR: unknown ic_mode in IO_LIST_SIBDRV: ',trim(ic_mode)
       print*,'       Valid values are: requib, restart'
       stop
    end select

    if (minveg_diag) then
       if (spinup .or. ic_restart) then
          print*,'ERROR: minveg_diag requires transient ic_mode=requib (not spinup/restart)'
          stop
       endif
       if (trim(minveg_clim_policy) /= 'preserve' .and. &
           trim(minveg_clim_policy) /= 'zero') then
          print*,'ERROR: minveg_clim_policy must be preserve or zero'
          stop
       endif
       print*,'MINVEG DIAGNOSTIC enabled; climatology policy: ',trim(minveg_clim_policy)
    endif

    ! Step 1C is strictly a one-site, post-fire, restart-based experiment.
    if (pftswap_enable) then
       if (.not. ic_restart .or. spinup .or. minveg_diag) then
          error stop 'PFTSWAP opening requires transient ic_mode=restart and minveg_diag=.false.'
       endif
       if (.not. fire_switch .or. firegfed_switch) then
          error stop 'PFTSWAP opening requires standard-fire only for this experiment'
       endif
       if (pftswap_open_year .lt. startyear .or. &
           pftswap_open_year .gt. endyear .or. &
           pftswap_open_doy .lt. 1 .or. pftswap_open_doy .gt. 366) then
          error stop 'PFTSWAP requires an opening year and valid DOY within segment'
       endif
       print '(a,2i7)', 'PFTSWAP armed for open year/DOY: ', &
             pftswap_open_year, pftswap_open_doy
    endif

    if (spinup .and. ic_restart) then
       print*,'ERROR: ic_mode=restart is reserved for transient continuation.'
       print*,'       Use the existing spinup_continue pathway for spinup continuation.'
       stop
    endif

    print ('(a,2i6)'), '   Starting Time (DOY/Year): ',starttime, startyear
    print ('(a,2i6)'), '   Ending Time (DOY/Year): ',endtime, endyear
    print ('(a,i6)'),  '   SiB time step (s) = ',dtsib
    print ('(2a)'),    '   IC mode = ',trim(ic_mode)

    driver_path = dr_path
    fire_path = fr_path
    firegfed_path = frba_path
!    ic_file = icfile       

    if ((spinup) .and. (.not. spinup_writediag)) then
       restart_savef = .false.
       hr_savegf   = .false.
       hr_saveluf  = .false.
       qp_savegf   = .false.
       qp_saveluf  = .false.
       pbp_savegf  = .false.
       pbp_saveluf = .false.
    else
       restart_savef = restart_dtsib /= 0
       hr_savegf  = (hr_dtsib /= 0)
       hr_saveluf = hr_savegf
       qp_savegf  = (qp_dtsib /= 0)
       qp_saveluf = qp_savegf
       pbp_savegf = (pbp_dtsib /= 0)
       pbp_saveluf = pbp_savegf
    endif

    single_pt = .false.
    if (nsib == 1) then
        single_pt=.true. 
        printrank = .false.
    endif

end subroutine read_namel
