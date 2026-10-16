subroutine restart_read()

    !============================================================
    !  This subroutine reads in the initial conditions file and
    !    stores the variables into the sib structure.
    !=============================================================

    use netcdf 
    use kinds
    use module_io
    use module_pparams, only:  &
       bcosm, p0_sfc, secs_per_day
    use module_sibconst, only: &
       eqclear_switch,   &
       nsib, npft, nlu, &
       npoolpft, npoollu, npoolcan, &
       nsoil, nsnow, ntot, npstgmax, &
       subcount, subset, subpref, sublarea, &
       spinup, spinup_default, &
       npoolcanc13, npoolcanc14, &
       spinup_continue, fire_switch, firegfed_switch
    use module_sib, only: &
       sib, gprog_type
    use module_pftinfo, only: &
       pft_num
    use module_phosib, only: &
        pressure
    use module_time, only: &
        startyear 

    implicit none

    !...netcdf variables
    integer(i4) :: ncid, status
    integer(i4) :: dimid, varid
    integer(i4) :: numdims
 
    !...local variables
    integer(i4) :: nsibf, ntest
    character(len=16) :: restart_fire_type
    integer(i4) :: iatt
    real(r8) :: ratt
    logical :: nsibr_flag
    logical :: file_exist, var_exist
    
    integer(i4) :: i,l,p,v
    integer(i4) :: gref
    integer(i4) :: swschema, swphaseid, swsourceid, swliveid, swdeadid
    integer(i4) :: swibuf(1,1)
    real(r8), allocatable :: swlivebuf(:,:,:,:), swdeadbuf(:,:,:,:)

    integer(byte) :: btvar
    integer(i4), dimension(:,:), allocatable :: ivar2d
    real(r4), dimension(:,:), allocatable :: tvar2d
    integer(i4), dimension(:,:,:), allocatable :: ivar3d

    !real(r4), dimension(:), allocatable :: tvar
    !real(r4), dimension(:,:,:), allocatable :: tvar3d
    !real(r4), dimension(:,:,:,:), allocatable :: tvar4d

    real(r8), dimension(:), allocatable :: dvar
    real(r8), dimension(:,:), allocatable :: dvar2d
    real(r8), dimension(:,:,:), allocatable :: dvar3d
    real(r8), dimension(:,:,:,:), allocatable :: dvar4d

!    type(gprog_type), intent(in) :: gprogt
!    pressure = gprogt%ps * 100.0    
!--------------------------------------------
print*,'Setting SiB4 Restart Variables'

! Enquire if initial conditions file exists
file_exist = .true.
inquire(file=trim(ic_file),exist=file_exist)
if (spinup .and. (spinup_default) .and. &
    (.not. (spinup_continue))) file_exist=.false.
if (spinup .and. (spinup_default) .and. &
    (spinup_continue)) file_exist=.true.

! Open and read initial conditions
if (.not. file_exist) then
   print*,'   Using Default Initial Conditions'
   ncid = -1
else
   print*,'   Reading IC File: '
   print*,'     ',trim(ic_file)
   status = nf90_open( trim(ic_file), nf90_nowrite, ncid )
   if (status /= nf90_noerr) then
      print*,'Error Opening IC File.'
      stop
   endif

   if (ic_restart) then
      ! Exact-continuation driver control state.  driver_init/fire_init
      ! establish open files first; these attributes then restore the exact
      ! record counters, bracket times, and read/switch flags of the checkpoint.
      status = nf90_get_att(ncid,nf90_global,'driver_step',driver_step)
      if (status .ne. nf90_noerr) stop 'Missing continuation driver_step attribute'
      status = nf90_get_att(ncid,nf90_global,'driver_perday',driver_perday)
      if (status .ne. nf90_noerr) stop 'Missing continuation driver_perday attribute'
      status = nf90_get_att(ncid,nf90_global,'driver_permon',driver_permon)
      if (status .ne. nf90_noerr) stop 'Missing continuation driver_permon attribute'
      status = nf90_get_att(ncid,nf90_global,'driver_recnum',driver_recnum)
      if (status .ne. nf90_noerr) stop 'Missing continuation driver_recnum attribute'
      status = nf90_get_att(ncid,nf90_global,'driver_year',driver_year)
      if (status .ne. nf90_noerr) stop 'Missing continuation driver_year attribute'
      status = nf90_get_att(ncid,nf90_global,'driver_month',driver_month)
      if (status .ne. nf90_noerr) stop 'Missing continuation driver_month attribute'
      status = nf90_get_att(ncid,nf90_global,'driver_day',driver_day)
      if (status .ne. nf90_noerr) stop 'Missing continuation driver_day attribute'
      status = nf90_get_att(ncid,nf90_global,'driver_hour',driver_hour)
      if (status .ne. nf90_noerr) stop 'Missing continuation driver_hour attribute'
      status = nf90_get_att(ncid,nf90_global,'driver_seccur',ratt)
      if (status .ne. nf90_noerr) stop 'Missing continuation driver_seccur attribute'
      driver_seccur = int(ratt,kind=long)
      status = nf90_get_att(ncid,nf90_global,'driver_secnext',ratt)
      if (status .ne. nf90_noerr) stop 'Missing continuation driver_secnext attribute'
      driver_secnext = int(ratt,kind=long)
      status = nf90_get_att(ncid,nf90_global,'driver_readf',iatt)
      if (status .ne. nf90_noerr) stop 'Missing continuation driver_readf attribute'
      driver_readf = (iatt.ne.0)
      status = nf90_get_att(ncid,nf90_global,'driver_switchf',iatt)
      if (status .ne. nf90_noerr) stop 'Missing continuation driver_switchf attribute'
      driver_switchf = (iatt.ne.0)
      status = nf90_get_att(ncid,nf90_global,'driver_updatef',iatt)
      if (status .ne. nf90_noerr) stop 'Missing continuation driver_updatef attribute'
      driver_updatef = (iatt.ne.0)

      restart_fire_type = ''
      status = nf90_get_att(ncid,nf90_global,'fire_driver_type',restart_fire_type)
      if (status.ne.nf90_noerr) stop 'Missing continuation fire_driver_type attribute'
      if (fire_switch .and. firegfed_switch) then
         if (trim(adjustl(restart_fire_type)).ne.'both') stop 'Restart fire-driver mode does not match namelist'
      elseif (fire_switch) then
         if (trim(adjustl(restart_fire_type)).ne.'standard') stop 'Restart fire-driver mode does not match namelist'
      elseif (firegfed_switch) then
         if (trim(adjustl(restart_fire_type)).ne.'ba') stop 'Restart fire-driver mode does not match namelist'
      else
         if (trim(adjustl(restart_fire_type)).ne.'none') stop 'Restart fire-driver mode does not match namelist'
      endif

      if (fire_switch .or. firegfed_switch) then
         status = nf90_get_att(ncid,nf90_global,'fire_step',fire_step)
         if (status .ne. nf90_noerr) stop 'Missing continuation fire_step attribute'
         status = nf90_get_att(ncid,nf90_global,'fire_perday',fire_perday)
         if (status .ne. nf90_noerr) stop 'Missing continuation fire_perday attribute'
         status = nf90_get_att(ncid,nf90_global,'fire_permon',fire_permon)
         if (status .ne. nf90_noerr) stop 'Missing continuation fire_permon attribute'
         status = nf90_get_att(ncid,nf90_global,'fire_recnum',fire_recnum)
         if (status .ne. nf90_noerr) stop 'Missing continuation fire_recnum attribute'
         status = nf90_get_att(ncid,nf90_global,'fire_year',fire_year)
         if (status .ne. nf90_noerr) stop 'Missing continuation fire_year attribute'
         status = nf90_get_att(ncid,nf90_global,'fire_month',fire_month)
         if (status .ne. nf90_noerr) stop 'Missing continuation fire_month attribute'
         status = nf90_get_att(ncid,nf90_global,'fire_day',fire_day)
         if (status .ne. nf90_noerr) stop 'Missing continuation fire_day attribute'
         status = nf90_get_att(ncid,nf90_global,'fire_hour',fire_hour)
         if (status .ne. nf90_noerr) stop 'Missing continuation fire_hour attribute'
         status = nf90_get_att(ncid,nf90_global,'fire_seccur',ratt)
         if (status .ne. nf90_noerr) stop 'Missing continuation fire_seccur attribute'
         fire_seccur = int(ratt,kind=long)
         status = nf90_get_att(ncid,nf90_global,'fire_secnext',ratt)
         if (status .ne. nf90_noerr) stop 'Missing continuation fire_secnext attribute'
         fire_secnext = int(ratt,kind=long)
      endif
      if (fire_switch) then
         status = nf90_get_att(ncid,nf90_global,'fire_readf',iatt)
         if (status .ne. nf90_noerr) stop 'Missing continuation fire_readf attribute'
         fire_readf = (iatt.ne.0)
         status = nf90_get_att(ncid,nf90_global,'fire_switchf',iatt)
         if (status .ne. nf90_noerr) stop 'Missing continuation fire_switchf attribute'
         fire_switchf = (iatt.ne.0)
         status = nf90_get_att(ncid,nf90_global,'fire_updatef',iatt)
         if (status .ne. nf90_noerr) stop 'Missing continuation fire_updatef attribute'
         fire_updatef = (iatt.ne.0)
      endif
      if (firegfed_switch) then
         status = nf90_get_att(ncid,nf90_global,'firegfed_readf',iatt)
         if (status .ne. nf90_noerr) stop 'Missing continuation firegfed_readf attribute'
         firegfed_readf = (iatt.ne.0)
         status = nf90_get_att(ncid,nf90_global,'firegfed_switchf',iatt)
         if (status .ne. nf90_noerr) stop 'Missing continuation firegfed_switchf attribute'
         firegfed_switchf = (iatt.ne.0)
         status = nf90_get_att(ncid,nf90_global,'firegfed_updatef',iatt)
         if (status .ne. nf90_noerr) stop 'Missing continuation firegfed_updatef attribute'
         firegfed_updatef = (iatt.ne.0)
      endif
   endif
   
   ! Read and check dimensions
   status = nf90_inq_dimid( ncid, 'nsib', dimid )
   if (status == nf90_noerr) then 
       status = nf90_inquire_dimension( ncid, dimid, len=nsibf )
   endif
   !....currently require restart file to have nsib or subcount elements...
   if (nsibf .eq. nsib) then
      nsibr_flag=.true.
   elseif (nsibf .eq. subcount) then
      nsibr_flag=.false.
   else
       print*, 'NSIB INCORRECT IN RESTART FILE (IC_PATH)'
       print*,' Sim nsib=',nsib, ' File nsib=',ntest, 'Subcount=',subcount
       stop
   endif

   status = nf90_inq_dimid( ncid, 'npft', dimid )
   if (status == nf90_noerr) then
      status = nf90_inquire_dimension( ncid, dimid, len=ntest )
   endif
   if (npft /= ntest) then
       print*, 'NPFT IN RESTART FILE (IC_PATH) INCORRECT'
       print*,'  Sim npft=',npft, ' File npft=',ntest
       !stop
   endif

   status = nf90_inq_dimid( ncid, 'nlu', dimid )
   if (status == nf90_noerr) then
       status = nf90_inquire_dimension( ncid, dimid, len=ntest )
   endif
   if (nlu /= ntest) then
       print*, 'NLU IN RESTART FILE (IC_PATH) INCORRECT'
       print*,' Sim nlu=',nlu, ' File nlu=',ntest
       stop
   endif

   status = nf90_inq_dimid( ncid, 'npoolpft', dimid )
   if (status == nf90_noerr) then
       status = nf90_inquire_dimension( ncid, dimid, len=ntest )
   endif
   if (npoolpft /= ntest) then
       print*, 'NPOOLPFT IN RESTART FILE (IC_PATH) INCORRECT'
       print*,' Sim npoolpft=',npoolpft, ' File npoolpft=',ntest
       stop
   endif

   status = nf90_inq_dimid( ncid, 'npoollu', dimid )
   if (status == nf90_noerr) then
      status = nf90_inquire_dimension( ncid, dimid, len=ntest )
   endif
   if (npoollu /= ntest) then
       print*, 'NPOOLLU IN RESTART FILE (IC_PATH) INCORRECT'
       print*,' Sim npoollu=',npoollu, ' File npoollu=',ntest
       stop
   endif

   status = nf90_inq_dimid( ncid, 'nsoil', dimid )
   if (status == nf90_noerr) then
      status = nf90_inquire_dimension( ncid, dimid, len=ntest )
   endif
   if (nsoil /= ntest) then
       print*, 'NSOIL IN RESTART FILE (IC_PATH) INCORRECT'
       print*,' Sim nsoil=',nsoil, ' File nsoil=',ntest
       stop
   endif

   status = nf90_inq_dimid( ncid, 'nsnow', dimid )
   if (status == nf90_noerr) then
       status = nf90_inquire_dimension( ncid, dimid, len=ntest )
   endif
   if (nsnow /= ntest) then
       print*, 'NSNOW IN RESTART FILE (IC_PATH) INCORRECT'
       print*,' Sim nsnow=',nsnow, ' File nsnow=',ntest
       stop
   endif
endif !ic file exists test

!---------------------------------------------------------------
! Restore the current runtime vegetation configuration for exact
! transient continuation.  sib_vs remains the authoritative initial
! vegetation specification for fresh/requib starts; for ic_restart,
! the checkpoint PFT identity and LU area override the sib_vs values.
!
! Step-1 dynamic-vegetation scope: the number of allocated LUs is still
! established by sib_vs in sibtype_init().  This supports an in-place
! PFT transition of an existing LU (the Figure 8A end member), but not
! yet creation of additional active LUs for mixed-PFT transitions.
if (file_exist .and. ic_restart) then
   status = nf90_inq_varid(ncid, 'lu_pftref', varid)
   if (status .ne. nf90_noerr) then
      stop 'Missing lu_pftref in continuation restart'
   endif
   allocate(ivar2d(nsibf,nlu))
   status = nf90_get_var(ncid, varid, ivar2d)
   if (status .ne. nf90_noerr) then
      stop 'Error reading lu_pftref from continuation restart'
   endif

   status = nf90_inq_varid(ncid, 'lu_area', varid)
   if (status .ne. nf90_noerr) then
      stop 'Missing lu_area in continuation restart'
   endif
   allocate(tvar2d(nsibf,nlu))
   status = nf90_get_var(ncid, varid, tvar2d)
   if (status .ne. nf90_noerr) then
      stop 'Error reading lu_area from continuation restart'
   endif

   do i=1,subcount
      if (nsibr_flag) then
         gref = subset(i)
      else
         gref = i
      endif

      do l=1,sib%g(i)%g_nlu
         if ((ivar2d(gref,l) .le. 0) .or. &
             (pft_num(ivar2d(gref,l)) .le. 0)) then
            stop 'Invalid active lu_pftref in continuation restart'
         endif
         if ((tvar2d(gref,l) .lt. rzero) .or. &
             (tvar2d(gref,l) .gt. rone)) then
            stop 'Invalid active lu_area in continuation restart'
         endif

         sib%g(i)%l(l)%ipft  = ivar2d(gref,l)
         sib%g(i)%l(l)%larea = tvar2d(gref,l)

         ! Keep legacy/static output metadata synchronized at the start
         ! of a restarted segment.  Time-varying within-file vegetation
         ! diagnostics are a later Figure 8A output change.
         subpref(i,l)  = ivar2d(gref,l)
         sublarea(i,l) = tvar2d(gref,l)
      enddo
   enddo

   deallocate(ivar2d,tvar2d)
endif

!---------------------------------------------------------------
! Loop over SiB4 variables, setting where appropriate
do v = 1, sibr_nvar
   status = nf90_inq_varid(ncid, trim(sibr_vname(v)), varid)
   var_exist = (status .eq. nf90_noerr)

   select case (sibr_vref(v))

   case (1)
       if (var_exist) then
           allocate(dvar(nsibf))
           status = nf90_get_var(ncid, varid, dvar)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              sib%g(i)%gprogt%tmd = dvar(gref)
           enddo
           deallocate(dvar)
       else
           if (sibr_vd(v) .lt. rzero) then
              do i=1, subcount
                 sib%g(i)%gprogt%tmd =  sib%g(i)%gprogt%tm
              enddo
           else
              do i=1, subcount
                 sib%g(i)%gprogt%tmd = sibr_vd(v)
              enddo
           endif
       endif

       
   case (2)
       if (var_exist) then
           allocate(dvar(nsibf))
           status = nf90_get_var(ncid, varid, dvar)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              sib%g(i)%gprogt%seas_precip = dvar(gref)
           enddo
           deallocate(dvar)
       else
           if (sibr_vd(v) .lt. rzero) then
              do i=1, subcount
                 sib%g(i)%gprogt%seas_precip = MAX(0.01, &
                    (sib%g(i)%gprogt%lspr + sib%g(i)%gprogt%cupr)*secs_per_day)
              enddo
           else
              do i=1, subcount
                 sib%g(i)%gprogt%seas_precip = MAX(0.01, sibr_vd(v))
              enddo
           endif
       endif

   case (3)
       if (var_exist) then
           allocate(dvar(nsibf))
           status = nf90_get_var(ncid, varid, dvar)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              sib%g(i)%gprogt%seas_tm = dvar(gref)
           enddo
           deallocate(dvar)
       else
           do i=1, subcount
              if (sibr_vd(v) .lt. rzero) then
                  sib%g(i)%gprogt%seas_tm = sibr_vd(v)
              else
                  sib%g(i)%gprogt%seas_tm = sib%g(i)%gprogt%tm
              endif
           enddo
       endif
       
   case (4)
       if (var_exist) then
           allocate(dvar(nsibf))
           status = nf90_get_var(ncid, varid, dvar)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              print*,'  npts/varid: ',nsibf,varid
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              sib%g(i)%gprogt%clim_cupr = dvar(gref)
           enddo
           deallocate(dvar)
       else
           if (sibr_vd(v) .lt. rzero) then
              do i=1, subcount
                 sib%g(i)%gprogt%clim_cupr = MAX(0.01, &
                    sib%g(i)%gprogt%cupr*secs_per_day)
              enddo
           else
              do i=1, subcount
                 sib%g(i)%gprogt%clim_cupr = MAX(0.01, sibr_vd(v))
              enddo
           endif
       endif

   case (5)
       if (var_exist) then
           allocate(dvar(nsibf))
           status = nf90_get_var(ncid, varid, dvar)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              sib%g(i)%gprogt%clim_precip = dvar(gref)
           enddo
           deallocate(dvar)
       else
           if (sibr_vd(v) .lt. rzero) then
              do i=1, subcount
                 sib%g(i)%gprogt%clim_precip = MAX(0.01, &
                    (sib%g(i)%gprogt%lspr + sib%g(i)%gprogt%cupr)*secs_per_day)
              enddo
           else
              do i=1, subcount
                 sib%g(i)%gprogt%clim_precip = MAX(0.01, sibr_vd(v))
              enddo
           endif
       endif
       
   case (6)
       if (var_exist) then
           allocate(dvar(nsibf))
           status = nf90_get_var(ncid, varid, dvar)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              sib%g(i)%gprogt%clim_tm = dvar(gref)
           enddo
           deallocate(dvar)
       else
           do i=1, subcount
              if (sibr_vd(v) .lt. rzero) then
                  sib%g(i)%gprogt%clim_tm = sibr_vd(v)
              else
                  sib%g(i)%gprogt%clim_tm = sib%g(i)%gprogt%tm
              endif
           enddo
       endif

   case (10)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%cast%eacas = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%cast%eacas = sibr_vd(v)
              enddo
           enddo
       endif

   case (11)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%cast%shcas = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           if (sibr_vd(v) .lt. rzero) then
              do i=1, subcount
                 do l=1,sib%g(i)%g_nlu
                    sib%g(i)%l(l)%cast%shcas = sib%g(i)%gprogt%sh
                 enddo
              enddo
           else
              do i=1, subcount
                 do l=1,sib%g(i)%g_nlu
                     sib%g(i)%l(l)%cast%eacas = sibr_vd(v)
                 enddo
              enddo
           endif
       endif

   case (12)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%cast%tcas = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           if (sibr_vd(v) .lt. rzero) then
              do i=1, subcount
                 do l=1,sib%g(i)%g_nlu
                    sib%g(i)%l(l)%cast%tcas = sib%g(i)%gprogt%tm
                 enddo
              enddo
           else
              do i=1, subcount
                 do l=1,sib%g(i)%g_nlu
                     sib%g(i)%l(l)%cast%tcas = sibr_vd(v)
                 enddo
              enddo
           endif
       endif

   case (13)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%cast%tkecas = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%cast%tkecas = sibr_vd(v)
              enddo
           enddo
       endif

   case (14)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%cast%tc = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           if (sibr_vd(v) .lt. rzero) then
              do i=1, subcount
                 do l=1,sib%g(i)%g_nlu
                    sib%g(i)%l(l)%cast%tc = sib%g(i)%gprogt%tm
                 enddo
              enddo
           else
              do i=1, subcount
                 do l=1,sib%g(i)%g_nlu
                     sib%g(i)%l(l)%cast%tc = sibr_vd(v)
                 enddo
              enddo
           endif
       endif
       
   case (15)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%cast%tcmin = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           if (sibr_vd(v) .lt. rzero) then
              do i=1, subcount
                 do l=1,sib%g(i)%g_nlu
                    sib%g(i)%l(l)%cast%tcmin = sib%g(i)%gprogt%tm
                 enddo
              enddo
           else
              do i=1, subcount
                 do l=1,sib%g(i)%g_nlu
                     sib%g(i)%l(l)%cast%tcmin = sibr_vd(v)
                 enddo
              enddo
           endif
       endif

   case (20)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%co2t%co2m = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           if (sibr_vd(v) .lt. rzero) then
              do i=1, subcount
                 do l=1,sib%g(i)%g_nlu
                    sib%g(i)%l(l)%co2t%co2m = sib%g(i)%gprogt%co2m
                 enddo
              enddo
           else
              do i=1, subcount
                 do l=1,sib%g(i)%g_nlu
                     sib%g(i)%l(l)%co2t%co2m = sibr_vd(v)
                 enddo
              enddo
           endif
       endif

   !case (20)
   !    if (var_exist) then
   !        allocate(dvar2d(nsibf,nlu))
   !        status = nf90_get_var(ncid, varid, dvar2d)
   !        if (status .ne. nf90_noerr) then
   !           print*,'Error reading restart variable: ', &
   !               trim(sibr_vname(v))
   !           stop
   !        endif

   !        do i=1, subcount
   !           if (nsibr_flag) then
   !               gref = subset(i)
   !           else
   !               gref = i
   !           endif

   !           do l=1,sib%g(i)%g_nlu
   !               sib%g(i)%l(l)%co2t%co2m = dvar2d(gref,l)
   !           enddo
   !        enddo
   !        deallocate(dvar2d)
   !    else
   !        if (sibr_vd(v) .lt. rzero) then
   !           do i=1, subcount
   !              do l=1,sib%g(i)%g_nlu
   !                 sib%g(i)%l(l)%co2t%co2m = dble(bco2m)
   !              enddo
   !           enddo
   !        else
   !           do i=1, subcount
   !              do l=1,sib%g(i)%g_nlu
   !                  sib%g(i)%l(l)%co2t%co2m = sibr_vd(v)
   !              enddo
   !           enddo
   !        endif
   !    endif
       
   !case (20)
   !    if (var_exist) then
   !        allocate(dvar(nsibf))
   !        status = nf90_get_var(ncid, varid, dvar)
   !        if (status .ne. nf90_noerr) then
   !           print*,'Error reading restart variable: ', &
   !               trim(sibr_vname(v))
   !           stop
   !        endif

   !        do i=1, subcount
   !           if (nsibr_flag) then
   !               gref = subset(i)
   !           else
   !               gref = i
   !           endif

   !           sib%g(i)%gprogt%co2m = dvar(gref)
   !        enddo
   !        deallocate(dvar)
   !    else
   !        if (sibr_vd(v) .lt. rzero) then
   !           do i=1,subcount
   !              !sib%g(i)%gprogt%pco2m = &
   !              !    (bco2m * p0_sfc/1.E6)
   !              !sib%g(i)%gprogt%pco2m = &
   !              !    dble(bco2m * pressure/1.E6)
   !              sib%g(i)%gprogt%co2m = dble(bco2m) 
   !              !print*,'sibr_vd(v) :',sibr_vd(v)
   !              !print*,'bco2m :',bco2m
   !              !print*,'p0_sfc :',p0_sfc
   !              !print*,'pressure :',pressure
   !              !print*,'pco2m rr 1 :',sib%g(i)%gprogt%pco2m
   !           enddo
   !        else
   !           do i=1, subcount
   !              sib%g(i)%gprogt%co2m = sibr_vd(v)
   !              !print*,'sibr_vd(v) :',sibr_vd(v)
   !              !print*,'pco2m rr 2 :',sib%g(i)%gprogt%pco2m
   !           enddo
   !        endif
   !    endif
       
   case (21)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%co2t%pco2cas = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%co2t%pco2cas = sibr_vd(v)
              enddo
            enddo
       endif

   case (22)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%co2t%pco2c = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%co2t%pco2c = sibr_vd(v)
              enddo
            enddo
         endif

    case (23)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%co2t%pco2i = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%co2t%pco2i = sibr_vd(v)
              enddo
            enddo
         endif

   case (24)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%co2t%pco2s = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%co2t%pco2s = sibr_vd(v)
              enddo
            enddo
       endif
       
   case (25)
       if (var_exist) then
           allocate(dvar(nsibf))
           status = nf90_get_var(ncid, varid, dvar)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              sib%g(i)%gprogt%pcosm = dvar(gref)
           enddo
           deallocate(dvar)
       else
           if (sibr_vd(v) .lt. rzero) then
              do i=1, subcount
                 sib%g(i)%gprogt%pcosm = &
                      dble(bcosm * p0_sfc/1.E12)
              enddo
           else
              do i=1, subcount
                 sib%g(i)%gprogt%pcosm = sibr_vd(v)
              enddo
           endif
       endif
       
   case (26)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%cost%coscasp = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%cost%coscasp = sibr_vd(v)
              enddo
            enddo
       endif

   case (27)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%cost%coscas = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%cost%coscas = sibr_vd(v)
              enddo
            enddo
       endif

   case (28)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%co2t%co2cas = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%co2t%co2cas = sibr_vd(v)
              enddo
            enddo
       endif
       
   case (29)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif  
              
           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif
                  
              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%co2t%co2c = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%co2t%co2c = sibr_vd(v)
              enddo
            enddo
       endif

   case (30)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif  
              
           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif
                  
              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%co2t%co2i = dvar2d(gref,l)
              enddo
           enddo  
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%co2t%co2i = sibr_vd(v)
              enddo
            enddo
       endif

   case (31)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif  
              
           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif
                  
              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%co2t%co2s = dvar2d(gref,l)
              enddo
           enddo  
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%co2t%co2s = sibr_vd(v)
              enddo
            enddo
       endif

   case (35)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%fluxt%ra = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%fluxt%ra = sibr_vd(v)
              enddo
           enddo
       endif

   case (36)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%fluxt%rb = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%fluxt%rb = sibr_vd(v)
              enddo
           enddo
       endif
       
   case (37)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%co2t%rst = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%co2t%rst = sibr_vd(v)
              enddo
           enddo
       endif

   case (40)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%co2t%assim = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%co2t%assim = sibr_vd(v)
              enddo
           enddo
       endif

       
   case (41)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%co2t%assimd = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%co2t%assimd = sibr_vd(v)
              enddo
           enddo
       endif

   case (42)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%resp_nveg = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%poollt%resp_nveg = sibr_vd(v)
              enddo
           enddo
       endif
       
   case (43)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%resp_leaf = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%poollt%resp_leaf = sibr_vd(v)
              enddo
           enddo
       endif

   case (44)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%resp_root = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%poollt%resp_root = sibr_vd(v)
              enddo
           enddo
       endif
       
   case (45)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%resp_grow = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%poollt%resp_grow = sibr_vd(v)
              enddo
           enddo
       endif
       
   case (46)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%resp_mntn = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%poollt%resp_mntn = sibr_vd(v)
              enddo
           enddo
       endif
       
   case (47)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%resp_auto = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%poollt%resp_auto = sibr_vd(v)
              enddo
           enddo
       endif

   case (48)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%pooldt%resp_soil = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%pooldt%resp_soil = sibr_vd(v)
              enddo
           enddo
       endif

   case (49)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%pooldt%resp_het = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%pooldt%resp_het = sibr_vd(v)
              enddo
           enddo
       endif
      
   case (50)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%resp_grz = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%poollt%resp_grz = sibr_vd(v)
              enddo
           enddo
       endif
 
   case (51)
       ! Historically this variable was written but not restored.  Preserve
       ! that requib behavior; exact continuation restores it explicitly.
       if (ic_restart) then
          if (.not. var_exist) then
             print*,'ERROR: continuation restart variable is missing: ', &
                    trim(sibr_vname(v))
             stop
          endif
          allocate(dvar2d(nsibf,nlu))
          status = nf90_get_var(ncid, varid, dvar2d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1, subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%poollt%resp_hrvst = dvar2d(gref,l)
             enddo
          enddo
          deallocate(dvar2d)
       endif

   case (60)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%hydrost%capacc_liq = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%hydrost%capacc_liq = sibr_vd(v)
              enddo
           enddo
       endif

   case (61)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%hydrost%capacc_snow = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%hydrost%capacc_snow = sibr_vd(v)
              enddo
           enddo
       endif

   case (62)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%hydrost%capacg = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%hydrost%capacg = sibr_vd(v)
              enddo
           enddo
       endif

   case (63)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%hydrovt%clim_pawfrw = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%hydrovt%clim_pawfrw = sibr_vd(v)
              enddo
           enddo
       endif

   case (64)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%hydrovt%clim_tawfrw = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%hydrovt%clim_tawfrw = sibr_vd(v)
              enddo
           enddo
       endif

   case (70)
       if (var_exist) then
           allocate(ivar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, ivar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  btvar=bnegone*abs(int(ivar2d(gref,l),kind=byte))
                  sib%g(i)%l(l)%sscolt%nsl = btvar
              enddo
           enddo
           deallocate(ivar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%sscolt%nsl = bnegone*abs(int(sibr_vd(v),kind=byte))
              enddo
           enddo
       endif

   case (71)
       if (var_exist) then
           allocate(dvar3d(nsibf,nlu,nsnow))
           status = nf90_get_var(ncid, varid, dvar3d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%sscolt%dz(-nsnow+1:0) = dvar3d(gref,l,1:nsnow)
              enddo
           enddo
           deallocate(dvar3d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%sscolt%dz(-nsnow+1:0) = sibr_vd(v)
              enddo
           enddo
       endif

   case (72)
       if (var_exist) then
           allocate(dvar3d(nsibf,nlu,ntot))
           status = nf90_get_var(ncid, varid, dvar3d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%sscolt%td(-nsnow+1:nsoil) = dvar3d(gref,l,1:ntot)
              enddo
           enddo
           deallocate(dvar3d)
       else
           if (sibr_vd(v) .lt. rzero) then
              do i=1, subcount
                 do l=1,sib%g(i)%g_nlu
                    btvar = sib%g(i)%l(l)%sscolt%nsl
                    sib%g(i)%l(l)%sscolt%td(btvar+1:nsoil) = &
                         sib%g(i)%gprogt%tm
                 enddo
              enddo
           else
              do i=1, subcount
                 do l=1,sib%g(i)%g_nlu
                    btvar = sib%g(i)%l(l)%sscolt%nsl
                    sib%g(i)%l(l)%sscolt%td(btvar+1:nsoil) = sibr_vd(v)
                 enddo
              enddo
           endif
       endif

   case (73)
       if (var_exist) then
           allocate(dvar3d(nsibf,nlu,ntot))
           status = nf90_get_var(ncid, varid, dvar3d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%sscolt%www_liq(-nsnow+1:nsoil) = dvar3d(gref,l,1:ntot)
              enddo
           enddo
           deallocate(dvar3d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%sscolt%www_liq(-nsnow+1:nsoil) = sibr_vd(v)
              enddo
           enddo
       endif

   case (74)
       if (var_exist) then
           allocate(dvar3d(nsibf,nlu,ntot))
           status = nf90_get_var(ncid, varid, dvar3d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%sscolt%www_ice(-nsnow+1:nsoil) = dvar3d(gref,l,1:ntot)
              enddo
           enddo
           deallocate(dvar3d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%sscolt%www_ice(-nsnow+1:nsoil) = sibr_vd(v)
              enddo
           enddo
       endif

   case (80)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%co2t%clim_assim = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%co2t%clim_assim = sibr_vd(v)
              enddo
           enddo
       endif

   case (81)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%vegt%clim_lai = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%vegt%clim_lai = sibr_vd(v)
              enddo
           enddo
       endif

   case (82)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%phent%phenave_assim = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%phent%phenave_assim = sibr_vd(v)
              enddo
           enddo
       endif

   case (83)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%phent%phenave_assimsm = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%phent%phenave_assimsm = sibr_vd(v)
              enddo
           enddo
       endif

   case (84)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%phent%phenave_pr = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%phent%phenave_pr = sibr_vd(v)
              enddo
           enddo
       endif

   case (85)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%phent%phenave_prsm = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%phent%phenave_prsm = sibr_vd(v)
              enddo
           enddo
       endif

   case (86)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%phent%phenave_prsdoy = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%phent%phenave_prsdoy = sibr_vd(v)
              enddo
           enddo
       endif

   case (87)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%phent%phenave_prcdoy = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%phent%phenave_prcdoy = sibr_vd(v)
              enddo
           enddo
       endif
       
   case (88)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%phent%phenave_tawftop = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%phent%phenave_tawftop = sibr_vd(v)
              enddo
           enddo
       endif

   case (89)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%phent%phenave_tm = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           if (sibr_vd(v) .lt. rzero) then
              do i=1, subcount
                 do l=1,sib%g(i)%g_nlu
                    sib%g(i)%l(l)%phent%phenave_tm = &
                        sib%g(i)%gprogt%tm
                 enddo
              enddo
           else
              do i=1, subcount
                 do l=1,sib%g(i)%g_nlu
                    sib%g(i)%l(l)%phent%phenave_tm = sibr_vd(v)
                 enddo
              enddo
           endif
       endif

   case (90)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%phent%phenc_climp = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%phent%phenc_climp = sibr_vd(v)
              enddo
           enddo
       endif

   case (91)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%phent%phenc_laimax = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%phent%phenc_laimax = sibr_vd(v)
              enddo
           enddo
       endif

   case (92)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%phent%phenc_laimin = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%phent%phenc_laimin = sibr_vd(v)
              enddo
           enddo
       endif

   case (93)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%phent%phenave_env = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%phent%phenave_env = sibr_vd(v)
              enddo
           enddo
       endif

   case (94)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%phent%phenave_wa = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%phent%phenave_wa = sibr_vd(v)
              enddo
           enddo
       endif

   case (95)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%phent%phenave_wacsm = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%phent%phenave_wacsm = sibr_vd(v)
              enddo
           enddo
       endif

   case (96)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%phent%phen_pi = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%phent%phen_pi = sibr_vd(v)
              enddo
           enddo
       endif
       
   case (97)
       if (var_exist) then
           allocate(ivar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, ivar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%phent%phen_istage = ivar2d(gref,l)
              enddo
           enddo
           deallocate(ivar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%phent%phen_istage = int(sibr_vd(v))
              enddo
           enddo
       endif
       
   case (100)
       if (var_exist) then
           allocate(ivar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, ivar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%phent%dapd = ivar2d(gref,l)
              enddo
           enddo
           deallocate(ivar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%phent%dapd = int(sibr_vd(v))
              enddo
           enddo
       endif

   case (101)
       if (var_exist) then
           allocate(ivar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, ivar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%phent%dapdaf = ivar2d(gref,l)
              enddo
           enddo
           deallocate(ivar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%phent%dapdaf = int(sibr_vd(v))
              enddo
           enddo
       endif
       
   case (102)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%phent%gdd = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%phent%gdd = sibr_vd(v)
              enddo
           enddo
       endif

   case (103)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%phent%seed_pool = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%phent%seed_pool = sibr_vd(v)
              enddo
           enddo
       endif

    case (110)
       if (.not. var_exist) then
          status = nf90_inq_varid(ncid, 'poolpft_lay', varid)
          var_exist = (status .eq. nf90_noerr)
       endif
       
       if (var_exist) then
          status = nf90_inquire_variable(ncid, varid, ndims=numdims)

          if (numdims .eq. 3) then
             allocate(dvar3d(nsibf,nlu,npoolpft))
             status = nf90_get_var(ncid, varid, dvar3d)
             if (status .ne. nf90_noerr) then
                print*,'Error reading restart variable: ', &
                    trim(sibr_vname(v))
                stop
             endif

             do i=1, subcount
                if (nsibr_flag) then
                    gref = subset(i)
                else
                    gref = i
                endif

                do l=1,sib%g(i)%g_nlu
                    sib%g(i)%l(l)%poollt%poolpft(:) = dvar3d(gref,l,:)
                enddo
             enddo
             deallocate(dvar3d)
          endif

          if (numdims .eq. 4) then
              allocate(dvar4d(nsibf,nlu,npoolpft,nsoil))
              status = nf90_get_var(ncid, varid, dvar4d)
              if (status .ne. nf90_noerr) then
                   print*,'Error reading poolpft_lay.'
                   stop
              endif

              do i=1, subcount
                 if (nsibr_flag) then
                     gref = subset(i)
                 else
                     gref = i
                 endif

                  do l=1,sib%g(i)%g_nlu
                     do p=1,npoolpft
                        sib%g(i)%l(l)%poollt%poolpft_lay(p,:) = &
                              dvar4d(gref,l,p,:)
                         sib%g(i)%l(l)%poollt%poolpft(p) = &
                              sum(dvar4d(gref,l,p,:))
                       enddo
                    enddo
                enddo
                deallocate(dvar4d)
             endif
          else  !var does not exist
           do i=1, subcount
               do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%poolpft(:) = sibr_vd(v)
               enddo
           enddo
       endif

    case (111)
       if (var_exist) then
          allocate(dvar4d(nsibf,nlu,npoolpft,nsoil))
          status = nf90_get_var(ncid, varid, dvar4d)
          if (status .ne. nf90_noerr) then
               print*,'Error reading poolpft_dgain.'
               stop
          endif

          do i=1, subcount
             if (nsibr_flag) then
                 gref = subset(i)
             else
                 gref = i
             endif

              do l=1,sib%g(i)%g_nlu
                 do p=1,npoolpft
                    sib%g(i)%l(l)%poollt%poolpft_dgain(p,:) = &
                          dvar4d(gref,l,p,:)
                  enddo
               enddo
           enddo
           deallocate(dvar4d)
             
       else  !var does not exist
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%poolpft_dgain(:,:) = sibr_vd(v)
               enddo
           enddo
       endif

    case (112)
       if (var_exist) then
          allocate(dvar4d(nsibf,nlu,npoolpft,nsoil))
          status = nf90_get_var(ncid, varid, dvar4d)
          if (status .ne. nf90_noerr) then
               print*,'Error reading poolpft_dloss.'
               stop
          endif

          do i=1, subcount
             if (nsibr_flag) then
                 gref = subset(i)
             else
                 gref = i
             endif

              do l=1,sib%g(i)%g_nlu
                 do p=1,npoolpft
                    sib%g(i)%l(l)%poollt%poolpft_dloss(p,:) = &
                          dvar4d(gref,l,p,:)
                  enddo
               enddo
           enddo
           deallocate(dvar4d)
             
       else  !var does not exist
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%poolpft_dloss(:,:) = sibr_vd(v)
               enddo
           enddo
       endif

    case (113)
       if (.not. var_exist) then
          status = nf90_inq_varid(ncid, 'gain_assim', varid)
          var_exist = (status .eq. nf90_noerr)
       endif
       
       if (var_exist) then
          allocate(dvar3d(nsibf,nlu,npoolpft))
          status = nf90_get_var(ncid, varid, dvar3d)
          if (status .ne. nf90_noerr) then
               print*,'Error reading gain_assim.'
               stop
          endif

          do i=1, subcount
             if (nsibr_flag) then
                 gref = subset(i)
             else
                 gref = i
             endif

              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%poollt%gain_assim = dvar3d(gref,l,:)
              enddo
           enddo
           deallocate(dvar3d)
             
       else  !var does not exist
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%gain_assim(:) = sibr_vd(v)
               enddo
           enddo
       endif

    case (114)
       if (.not. var_exist) then
          status = nf90_inq_varid(ncid, 'gain_seed', varid)
          var_exist = (status .eq. nf90_noerr)
       endif
       
       if (var_exist) then
          allocate(dvar3d(nsibf,nlu,npoolpft))
          status = nf90_get_var(ncid, varid, dvar3d)
          if (status .ne. nf90_noerr) then
               print*,'Error reading gain_seed.'
               stop
          endif

          do i=1, subcount
             if (nsibr_flag) then
                 gref = subset(i)
             else
                 gref = i
             endif

              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%poollt%gain_seed = dvar3d(gref,l,:)
              enddo
           enddo
           deallocate(dvar3d)
             
       else  !var does not exist
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%gain_seed(:) = sibr_vd(v)
               enddo
           enddo
       endif

    case (115)
       if (.not. var_exist) then
          status = nf90_inq_varid(ncid, 'loss_gresp', varid)
          var_exist = (status .eq. nf90_noerr)
       endif
       
       if (var_exist) then
          allocate(dvar3d(nsibf,nlu,npoolpft))
          status = nf90_get_var(ncid, varid, dvar3d)
          if (status .ne. nf90_noerr) then
               print*,'Error reading loss_gresp.'
               stop
          endif

          do i=1, subcount
             if (nsibr_flag) then
                 gref = subset(i)
             else
                 gref = i
             endif

              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%poollt%loss_gresp = dvar3d(gref,l,:)
              enddo
           enddo
           deallocate(dvar3d)
             
       else  !var does not exist
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%loss_gresp(:) = sibr_vd(v)
               enddo
           enddo
       endif

    case (116)
       if (.not. var_exist) then
          status = nf90_inq_varid(ncid, 'loss_mresp_lay', varid)
          var_exist = (status .eq. nf90_noerr)
       endif
       
       if (var_exist) then
          allocate(dvar4d(nsibf,nlu,npoolpft,nsoil))
          status = nf90_get_var(ncid, varid, dvar4d)
          if (status .ne. nf90_noerr) then
               print*,'Error reading loss_mresp_lay.'
               stop
          endif

          do i=1, subcount
             if (nsibr_flag) then
                 gref = subset(i)
             else
                 gref = i
             endif

             do l=1,sib%g(i)%g_nlu
                do p=1,npoolpft
                   sib%g(i)%l(l)%poollt%loss_mresp_lay(p,:) = dvar4d(gref,l,p,:)
                enddo
             enddo
           enddo
           deallocate(dvar4d)
             
       else  !var does not exist
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%loss_mresp_lay(:,:) = sibr_vd(v)
               enddo
           enddo
       endif

    case (117)
       if (.not. var_exist) then
          status = nf90_inq_varid(ncid, 'loss_trans_lay', varid)
          var_exist = (status .eq. nf90_noerr)
       endif
       
       if (var_exist) then
          allocate(dvar4d(nsibf,nlu,npoolpft,nsoil))
          status = nf90_get_var(ncid, varid, dvar4d)
          if (status .ne. nf90_noerr) then
               print*,'Error reading loss_trans_lay.'
               stop
          endif

          do i=1, subcount
             if (nsibr_flag) then
                 gref = subset(i)
             else
                 gref = i
             endif

             do l=1,sib%g(i)%g_nlu
                do p=1,npoolpft
                   sib%g(i)%l(l)%poollt%loss_trans_lay(p,:) = dvar4d(gref,l,p,:)
                enddo
             enddo
           enddo
           deallocate(dvar4d)
             
       else  !var does not exist
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%loss_trans_lay(:,:) = sibr_vd(v)
               enddo
           enddo
       endif

    case (118)
       if (.not. var_exist) then
          status = nf90_inq_varid(ncid, 'loss_fire_lay', varid)
          var_exist = (status .eq. nf90_noerr)
       endif
       
       if (var_exist) then
          allocate(dvar4d(nsibf,nlu,npoolpft,nsoil))
          status = nf90_get_var(ncid, varid, dvar4d)
          if (status .ne. nf90_noerr) then
               print*,'Error reading loss_fire_lay.'
               stop
          endif

          do i=1, subcount
             if (nsibr_flag) then
                 gref = subset(i)
             else
                 gref = i
             endif

             do l=1,sib%g(i)%g_nlu
                do p=1,npoolpft
                   sib%g(i)%l(l)%poollt%loss_fire_lay(p,:) = dvar4d(gref,l,p,:)
                enddo
             enddo
           enddo
           deallocate(dvar4d)
             
       else  !var does not exist
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%loss_fire_lay(:,:) = sibr_vd(v)
               enddo
           enddo
       endif

    case (119)
       if (.not. var_exist) then
          status = nf90_inq_varid(ncid, 'loss_grz', varid)
          var_exist = (status .eq. nf90_noerr)
       endif
       
       if (var_exist) then
          allocate(dvar3d(nsibf,nlu,npoolcan))
          status = nf90_get_var(ncid, varid, dvar3d)
          if (status .ne. nf90_noerr) then
               print*,'Error reading loss_grz.'
               stop
          endif

          do i=1, subcount
             if (nsibr_flag) then
                 gref = subset(i)
             else
                 gref = i
             endif

             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%poollt%loss_grz(:) = dvar3d(gref,l,:)
             enddo
           enddo
           deallocate(dvar3d)
             
       else  !var does not exist
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%loss_grz(:) = sibr_vd(v)
               enddo
           enddo
       endif

    case (120)
       if (.not. var_exist) then
          status = nf90_inq_varid(ncid, 'loss_hrvst_lay', varid)
          var_exist = (status .eq. nf90_noerr)
       endif
       
       if (var_exist) then
          allocate(dvar4d(nsibf,nlu,npoolpft,nsoil))
          status = nf90_get_var(ncid, varid, dvar4d)
          if (status .ne. nf90_noerr) then
               print*,'Error reading loss_hrvst_lay.'
               stop
          endif

          do i=1, subcount
             if (nsibr_flag) then
                 gref = subset(i)
             else
                 gref = i
             endif

             do l=1,sib%g(i)%g_nlu
                do p=1,npoolpft
                   sib%g(i)%l(l)%poollt%loss_hrvst_lay(p,:) = dvar4d(gref,l,p,:)
                enddo
             enddo
           enddo
           deallocate(dvar4d)
             
       else  !var does not exist
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%loss_hrvst_lay(:,:) = sibr_vd(v)
               enddo
           enddo
       endif
       
   case (121)
       !if (.not. eqclear_switch) then
          if (var_exist) then
              allocate(dvar3d(nsibf,nlu,npoolpft))
              status = nf90_get_var(ncid, varid, dvar3d)
              if (status .ne. nf90_noerr) then
                 print*,'Error reading restart variable: ', &
                     trim(sibr_vname(v))
                 stop
              endif

              do i=1, subcount
                 if (nsibr_flag) then
                     gref = subset(i)
                 else
                     gref = i
                 endif

                 do l=1,sib%g(i)%g_nlu
                     sib%g(i)%l(l)%equiblt%poolpft_totgain(:) = dvar3d(gref,l,:)
                 enddo
              enddo
              deallocate(dvar3d)
          else
              do i=1, subcount
                 do l=1,sib%g(i)%g_nlu
                    sib%g(i)%l(l)%equiblt%poolpft_totgain(:) = sibr_vd(v)
                 enddo
              enddo
          endif
      !endif !not eqclear_switch

   case (122)
       !if (.not. eqclear_switch) then
          if (var_exist) then
              allocate(dvar3d(nsibf,nlu,npoolpft))
              status = nf90_get_var(ncid, varid, dvar3d)
              if (status .ne. nf90_noerr) then
                 print*,'Error reading restart variable: ', &
                     trim(sibr_vname(v))
                 stop
              endif

              do i=1, subcount
                 if (nsibr_flag) then
                     gref = subset(i)
                 else
                     gref = i
                 endif

                 do l=1,sib%g(i)%g_nlu
                     sib%g(i)%l(l)%equiblt%poolpft_totloss(:) = dvar3d(gref,l,:)
                 enddo
              enddo
              deallocate(dvar3d)
          else
              do i=1, subcount
                 do l=1,sib%g(i)%g_nlu
                    sib%g(i)%l(l)%equiblt%poolpft_totloss(:) = sibr_vd(v)
                 enddo
              enddo
          endif
      !endif !not eqclear_switch

    case (123)
       if (var_exist) then
          allocate(dvar4d(nsibf,nlu,npoolpft,nsoil))
          status = nf90_get_var(ncid, varid, dvar4d)
          if (status .ne. nf90_noerr) then
               print*,'Error reading poolpft_lay.'
               stop
          endif

          do i=1, subcount
             if (nsibr_flag) then
                 gref = subset(i)
             else
                 gref = i
             endif

              do l=1,sib%g(i)%g_nlu
                 do p=1,npoolpft
                    sib%g(i)%l(l)%poollt%poolpft_lay(p,:) = &
                          dvar4d(gref,l,p,:)
                  enddo
               enddo
           enddo
           deallocate(dvar4d)

       else  !var does not exist
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%poolpft_lay(:,:) = sibr_vd(v)
               enddo
           enddo
       endif

    case (124)
       if (.not. var_exist) then
          status = nf90_inq_varid(ncid, 'loss_grzc13', varid)
          var_exist = (status .eq. nf90_noerr)
       endif

       if (var_exist) then
          allocate(dvar3d(nsibf,nlu,npoolcanc13))
          status = nf90_get_var(ncid, varid, dvar3d)
          if (status .ne. nf90_noerr) then
               print*,'Error reading loss_grzc13.'
               stop
          endif

          do i=1, subcount
             if (nsibr_flag) then
                 gref = subset(i)
             else
                 gref = i
             endif

             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%poollt%loss_grzc13(:) = dvar3d(gref,l,:)
             enddo
           enddo
           deallocate(dvar3d)

       else  !var does not exist
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%loss_grzc13(:) = sibr_vd(v)
               enddo
           enddo
       endif

   case (125)
          if (var_exist) then
              allocate(dvar3d(nsibf,nlu,npoolpft))
              status = nf90_get_var(ncid, varid, dvar3d)
              if (status .ne. nf90_noerr) then
                 print*,'Error reading restart variable: ', &
                     trim(sibr_vname(v))
                 stop
              endif

              do i=1, subcount
                 if (nsibr_flag) then
                     gref = subset(i)
                 else
                     gref = i
                 endif

                 do l=1,sib%g(i)%g_nlu
                     sib%g(i)%l(l)%equiblt%poolpft_equib(:) = dvar3d(gref,l,:)
                 enddo
              enddo
              deallocate(dvar3d)
          else
              do i=1, subcount
                 do l=1,sib%g(i)%g_nlu
                    sib%g(i)%l(l)%equiblt%poolpft_equib(:) = sibr_vd(v)
                 enddo
              enddo
          endif

   case (126)
          if (var_exist) then
              allocate(dvar3d(nsibf,nlu,npoolpft))
              status = nf90_get_var(ncid, varid, dvar3d)
              if (status .ne. nf90_noerr) then
                 print*,'Error reading restart variable: ', &
                     trim(sibr_vname(v))
                 stop
              endif

              do i=1, subcount
                 if (nsibr_flag) then
                     gref = subset(i)
                 else
                     gref = i
                 endif

                 do l=1,sib%g(i)%g_nlu
                     sib%g(i)%l(l)%poollt%poolpftmin(:) = dvar3d(gref,l,:)
                 enddo
              enddo
              deallocate(dvar3d)
          else
              do i=1, subcount
                 do l=1,sib%g(i)%g_nlu
                    sib%g(i)%l(l)%poollt%poolpftmin(:) = sibr_vd(v)
                 enddo
              enddo
          endif

   case (127)
          if (var_exist) then
              allocate(dvar3d(nsibf,nlu,npoolpft))
              status = nf90_get_var(ncid, varid, dvar3d)
              if (status .ne. nf90_noerr) then
                 print*,'Error reading restart variable: ', &
                     trim(sibr_vname(v))
                 stop
              endif

              do i=1, subcount
                 if (nsibr_flag) then
                     gref = subset(i)
                 else
                     gref = i
                 endif

                 do l=1,sib%g(i)%g_nlu
                     sib%g(i)%l(l)%equiblt%poolpft_init(:) = dvar3d(gref,l,:)
                 enddo
              enddo
              deallocate(dvar3d)
          else
              do i=1, subcount
                 do l=1,sib%g(i)%g_nlu
                    sib%g(i)%l(l)%equiblt%poolpft_init(:) = sibr_vd(v)
                 enddo
              enddo
          endif

   case (128)
          if (var_exist) then
              allocate(dvar3d(nsibf,nlu,npoolpft))
              status = nf90_get_var(ncid, varid, dvar3d)
              if (status .ne. nf90_noerr) then
                 print*,'Error reading restart variable: ', &
                     trim(sibr_vname(v))
                 stop
              endif

              do i=1, subcount
                 if (nsibr_flag) then
                     gref = subset(i)
                 else
                     gref = i
                 endif

                 do l=1,sib%g(i)%g_nlu
                     sib%g(i)%l(l)%equiblt%poolpft_end(:) = dvar3d(gref,l,:)
                 enddo
              enddo
              deallocate(dvar3d)
          else
              do i=1, subcount
                 do l=1,sib%g(i)%g_nlu
                    sib%g(i)%l(l)%equiblt%poolpft_end(:) = sibr_vd(v)
                 enddo
              enddo
          endif

    case (180)
       if (.not. var_exist) then
          status = nf90_inq_varid(ncid, 'loss_transf_pft_lay', varid)
          var_exist = (status .eq. nf90_noerr)
       endif

       if (var_exist) then
          allocate(dvar4d(nsibf,nlu,npoolpft,nsoil))
          status = nf90_get_var(ncid, varid, dvar4d)
          if (status .ne. nf90_noerr) then
               print*,'Error reading loss_transf_pft_lay.'
               stop
          endif

          do i=1, subcount
             if (nsibr_flag) then
                 gref = subset(i)
             else
                 gref = i
             endif

             do l=1,sib%g(i)%g_nlu
                do p=1,npoolpft
                   sib%g(i)%l(l)%poollt%loss_transf_lay(p,:) = dvar4d(gref,l,p,:)
                enddo
             enddo
           enddo
           deallocate(dvar4d)

       else  !var does not exist
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%loss_transf_lay(:,:) = sibr_vd(v)
               enddo
           enddo
       endif

    case (130)
       if (.not. var_exist) then
          status = nf90_inq_varid(ncid, 'poollu_lay', varid)
          var_exist = (status .eq. nf90_noerr)
       endif
       
       if (var_exist) then
          status = nf90_inquire_variable(ncid, varid, ndims=numdims)

         if (numdims .eq. 3) then
             allocate(dvar3d(nsibf,nlu,npoollu))
             status = nf90_get_var(ncid, varid, dvar3d)
             if (status .ne. nf90_noerr) then
                print*,'Error reading restart variable: ', &
                    trim(sibr_vname(v))
                stop
             endif

             do i=1, subcount
                if (nsibr_flag) then
                    gref = subset(i)
                else
                    gref = i
                endif

                do l=1,sib%g(i)%g_nlu
                    sib%g(i)%l(l)%pooldt%poollu(:) = dvar3d(gref,l,:)
                enddo
             enddo
             deallocate(dvar3d)
          endif

          if (numdims .eq. 4) then
              allocate(dvar4d(nsibf,nlu,npoollu,nsoil))
              status = nf90_get_var(ncid, varid, dvar4d)
              if (status .ne. nf90_noerr) then
                   print*,'Error reading poollu_lay.'
                   stop
              endif

              do i=1, subcount
                 if (nsibr_flag) then
                     gref = subset(i)
                 else
                     gref = i
                 endif

                  do l=1,sib%g(i)%g_nlu
                     do p=1,npoollu
                        sib%g(i)%l(l)%pooldt%poollu_lay(p,:) = &
                              dvar4d(gref,l,p,:)
                         sib%g(i)%l(l)%pooldt%poollu(p) = &
                              sum(dvar4d(gref,l,p,:))
                       enddo
                    enddo
                enddo
                deallocate(dvar4d)
             endif
            
       else  !var does not exist
           do i=1, subcount
               do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%pooldt%poollu(:) = sibr_vd(v)
               enddo
           enddo
       endif

    case (131)
       if (var_exist) then
           allocate(dvar4d(nsibf,nlu,npoollu,nsoil))
           status = nf90_get_var(ncid, varid, dvar4d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading poollu_dgain.'
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

               do l=1,sib%g(i)%g_nlu
                  do p=1,npoollu
                     sib%g(i)%l(l)%pooldt%poollu_dgain(p,:) = &
                              dvar4d(gref,l,p,:)
                   enddo
               enddo
            enddo
            deallocate(dvar4d)
            
       else  !var does not exist
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%pooldt%poollu_dgain(:,:) = sibr_vd(v)
              enddo
           enddo
       endif

    case (132)
       if (var_exist) then
           allocate(dvar4d(nsibf,nlu,npoollu,nsoil))
           status = nf90_get_var(ncid, varid, dvar4d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading poollu_dloss.'
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

               do l=1,sib%g(i)%g_nlu
                  do p=1,npoollu
                     sib%g(i)%l(l)%pooldt%poollu_dloss(p,:) = &
                              dvar4d(gref,l,p,:)
                   enddo
               enddo
            enddo
            deallocate(dvar4d)
            
       else  !var does not exist
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%pooldt%poollu_dloss(:,:) = sibr_vd(v)
              enddo
           enddo
       endif

    case (133)
       if (.not. var_exist) then
          status = nf90_inq_varid(ncid, 'gain_grz_lay', varid)
          var_exist = (status .eq. nf90_noerr)
       endif
       
       if (var_exist) then
           allocate(dvar4d(nsibf,nlu,npoollu,nsoil))
           status = nf90_get_var(ncid, varid, dvar4d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading gain_grz_lay.'
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

               do l=1,sib%g(i)%g_nlu
                  do p=1,npoollu
                     sib%g(i)%l(l)%pooldt%gain_grz_lay(p,:) = &
                              dvar4d(gref,l,p,:)
                   enddo
               enddo
            enddo
            deallocate(dvar4d)
            
       else  !var does not exist
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%pooldt%gain_grz_lay(:,:) = sibr_vd(v)
              enddo
           enddo
       endif

    case (134)
       if (.not. var_exist) then
          status = nf90_inq_varid(ncid, 'gain_hrvst_lay', varid)
          var_exist = (status .eq. nf90_noerr)
       endif
       
       if (var_exist) then
           allocate(dvar4d(nsibf,nlu,npoollu,nsoil))
           status = nf90_get_var(ncid, varid, dvar4d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading gain_hrvst_lay.'
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

               do l=1,sib%g(i)%g_nlu
                  do p=1,npoollu
                     sib%g(i)%l(l)%pooldt%gain_hrvst_lay(p,:) = &
                              dvar4d(gref,l,p,:)
                   enddo
               enddo
            enddo
            deallocate(dvar4d)
            
       else  !var does not exist
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%pooldt%gain_hrvst_lay(:,:) = sibr_vd(v)
              enddo
           enddo
       endif

    case (135)
       if (.not. var_exist) then
          status = nf90_inq_varid(ncid, 'gain_transl_lay', varid)
          var_exist = (status .eq. nf90_noerr)
       endif
       
       if (var_exist) then
           allocate(dvar4d(nsibf,nlu,npoollu,nsoil))
           status = nf90_get_var(ncid, varid, dvar4d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading gain_transl_lay.'
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

               do l=1,sib%g(i)%g_nlu
                  do p=1,npoollu
                     sib%g(i)%l(l)%pooldt%gain_transl_lay(p,:) = &
                              dvar4d(gref,l,p,:)
                   enddo
               enddo
            enddo
            deallocate(dvar4d)
            
       else  !var does not exist
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%pooldt%gain_transl_lay(:,:) = sibr_vd(v)
              enddo
           enddo
       endif

    case (136)
       if (.not. var_exist) then
          status = nf90_inq_varid(ncid, 'gain_transd_lay', varid)
          var_exist = (status .eq. nf90_noerr)
       endif
       
       if (var_exist) then
           allocate(dvar4d(nsibf,nlu,npoollu,nsoil))
           status = nf90_get_var(ncid, varid, dvar4d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading gain_transd_lay.'
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

               do l=1,sib%g(i)%g_nlu
                  do p=1,npoollu
                     sib%g(i)%l(l)%pooldt%gain_transd_lay(p,:) = &
                              dvar4d(gref,l,p,:)
                   enddo
               enddo
            enddo
            deallocate(dvar4d)
            
       else  !var does not exist
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%pooldt%gain_transd_lay(:,:) = sibr_vd(v)
              enddo
           enddo
       endif

    case (137)
       if (.not. var_exist) then
          status = nf90_inq_varid(ncid, 'loss_fire_lay', varid)
          var_exist = (status .eq. nf90_noerr)
       endif
       
       if (var_exist) then
           allocate(dvar4d(nsibf,nlu,npoollu,nsoil))
           status = nf90_get_var(ncid, varid, dvar4d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading loss_fire_lay.'
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

               do l=1,sib%g(i)%g_nlu
                  do p=1,npoollu
                     sib%g(i)%l(l)%pooldt%loss_fire_lay(p,:) = &
                              dvar4d(gref,l,p,:)
                   enddo
               enddo
            enddo
            deallocate(dvar4d)
            
       else  !var does not exist
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%pooldt%loss_fire_lay(:,:) = sibr_vd(v)
              enddo
           enddo
       endif

    case (138)
       if (.not. var_exist) then
          status = nf90_inq_varid(ncid, 'loss_resp_lay', varid)
          var_exist = (status .eq. nf90_noerr)
       endif
       
       if (var_exist) then
           allocate(dvar4d(nsibf,nlu,npoollu,nsoil))
           status = nf90_get_var(ncid, varid, dvar4d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading loss_resp_lay.'
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

               do l=1,sib%g(i)%g_nlu
                  do p=1,npoollu
                     sib%g(i)%l(l)%pooldt%loss_resp_lay(p,:) = &
                              dvar4d(gref,l,p,:)
                   enddo
               enddo
            enddo
            deallocate(dvar4d)
            
       else  !var does not exist
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%pooldt%loss_resp_lay(:,:) = sibr_vd(v)
              enddo
           enddo
       endif

    case (139)
       if (.not. var_exist) then
          status = nf90_inq_varid(ncid, 'loss_trans_lay', varid)
          var_exist = (status .eq. nf90_noerr)
       endif
       
       if (var_exist) then
           allocate(dvar4d(nsibf,nlu,npoollu,nsoil))
           status = nf90_get_var(ncid, varid, dvar4d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading loss_trans_lay.'
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

               do l=1,sib%g(i)%g_nlu
                  do p=1,npoollu
                     sib%g(i)%l(l)%pooldt%loss_trans_lay(p,:) = &
                              dvar4d(gref,l,p,:)
                   enddo
               enddo
            enddo
            deallocate(dvar4d)
            
       else  !var does not exist
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%pooldt%loss_trans_lay(:,:) = sibr_vd(v)
              enddo
           enddo
       endif
       
   case (140)
       !if (.not. eqclear_switch) then
          if (var_exist) then
              allocate(dvar3d(nsibf,nlu,npoollu))
              status = nf90_get_var(ncid, varid, dvar3d)
              if (status .ne. nf90_noerr) then
                 print*,'Error reading restart variable: ', &
                     trim(sibr_vname(v))
                 stop
              endif

              do i=1, subcount
                 if (nsibr_flag) then
                     gref = subset(i)
                 else
                     gref = i
                 endif

                 do l=1,sib%g(i)%g_nlu
                     sib%g(i)%l(l)%equibdt%poollu_totgain(:) = dvar3d(gref,l,:)
                 enddo
              enddo
              deallocate(dvar3d)
          else
              do i=1, subcount
                 do l=1,sib%g(i)%g_nlu
                    sib%g(i)%l(l)%equibdt%poollu_totgain(:) = sibr_vd(v)
                 enddo
              enddo
          endif
      !endif !not eqclear_switch

   case (141)
       !if (.not. eqclear_switch) then
          if (var_exist) then
              allocate(dvar3d(nsibf,nlu,npoollu))
              status = nf90_get_var(ncid, varid, dvar3d)
              if (status .ne. nf90_noerr) then
                 print*,'Error reading restart variable: ', &
                     trim(sibr_vname(v))
                 stop
              endif

              do i=1, subcount
                 if (nsibr_flag) then
                     gref = subset(i)
                 else
                     gref = i
                 endif

                 do l=1,sib%g(i)%g_nlu
                     sib%g(i)%l(l)%equibdt%poollu_totloss(:) = dvar3d(gref,l,:)
                 enddo
              enddo
              deallocate(dvar3d)
          else
              do i=1, subcount
                 do l=1,sib%g(i)%g_nlu
                    sib%g(i)%l(l)%equibdt%poollu_totloss(:) = sibr_vd(v)
                 enddo
              enddo
          endif
      !endif !not eqclear_switch

    case (142)
       if (var_exist) then
           allocate(dvar4d(nsibf,nlu,npoollu,nsoil))
           status = nf90_get_var(ncid, varid, dvar4d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading poollu_lay.'
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

               do l=1,sib%g(i)%g_nlu
                  do p=1,npoollu
                     sib%g(i)%l(l)%pooldt%poollu_lay(p,:) = &
                              dvar4d(gref,l,p,:)
                   enddo
               enddo
            enddo
            deallocate(dvar4d)

       else  !var does not exist
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%pooldt%poollu_lay(:,:) = sibr_vd(v)
              enddo
           enddo
       endif

   case (143)
          if (var_exist) then
              allocate(dvar3d(nsibf,nlu,npoollu))
              status = nf90_get_var(ncid, varid, dvar3d)
              if (status .ne. nf90_noerr) then
                 print*,'Error reading restart variable: ', &
                     trim(sibr_vname(v))
                 stop
              endif

              do i=1, subcount
                 if (nsibr_flag) then
                     gref = subset(i)
                 else
                     gref = i
                 endif

                 do l=1,sib%g(i)%g_nlu
                     sib%g(i)%l(l)%equibdt%poollu_equib(:) = dvar3d(gref,l,:)
                 enddo
              enddo
              deallocate(dvar3d)
          else
              do i=1, subcount
                 do l=1,sib%g(i)%g_nlu
                    sib%g(i)%l(l)%equibdt%poollu_equib(:) = sibr_vd(v)
                 enddo
              enddo
          endif

   case (144)
          if (var_exist) then
              allocate(dvar3d(nsibf,nlu,npoollu))
              status = nf90_get_var(ncid, varid, dvar3d)
              if (status .ne. nf90_noerr) then
                 print*,'Error reading restart variable: ', &
                     trim(sibr_vname(v))
                 stop
              endif

              do i=1, subcount
                 if (nsibr_flag) then
                     gref = subset(i)
                 else
                     gref = i
                 endif

                 do l=1,sib%g(i)%g_nlu
                     sib%g(i)%l(l)%equibdt%poollu_init(:) = dvar3d(gref,l,:)
                 enddo
              enddo
              deallocate(dvar3d)
          else
              do i=1, subcount
                 do l=1,sib%g(i)%g_nlu
                    sib%g(i)%l(l)%equibdt%poollu_init(:) = sibr_vd(v)
                 enddo
              enddo
          endif

   case (145)
          if (var_exist) then
              allocate(dvar3d(nsibf,nlu,npoollu))
              status = nf90_get_var(ncid, varid, dvar3d)
              if (status .ne. nf90_noerr) then
                 print*,'Error reading restart variable: ', &
                     trim(sibr_vname(v))
                 stop
              endif

              do i=1, subcount
                 if (nsibr_flag) then
                     gref = subset(i)
                 else
                     gref = i
                 endif

                 do l=1,sib%g(i)%g_nlu
                     sib%g(i)%l(l)%equibdt%poollu_end(:) = dvar3d(gref,l,:)
                 enddo
              enddo
              deallocate(dvar3d)
          else
              do i=1, subcount
                 do l=1,sib%g(i)%g_nlu
                    sib%g(i)%l(l)%equibdt%poollu_end(:) = sibr_vd(v)
                 enddo
              enddo
          endif

    case (181)
       if (.not. var_exist) then
          status = nf90_inq_varid(ncid, 'gain_transf_lu_lay', varid)
          var_exist = (status .eq. nf90_noerr)
       endif

       if (var_exist) then
           allocate(dvar4d(nsibf,nlu,npoollu,nsoil))
           status = nf90_get_var(ncid, varid, dvar4d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading gain_transf_lu_lay.'
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

               do l=1,sib%g(i)%g_nlu
                  do p=1,npoollu
                     sib%g(i)%l(l)%pooldt%gain_transf_lay(p,:) = &
                              dvar4d(gref,l,p,:)
                   enddo
               enddo
            enddo
            deallocate(dvar4d)

       else  !var does not exist
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%pooldt%gain_transf_lay(:,:) = sibr_vd(v)
              enddo
           enddo
       endif

   case (147)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%fract%kiecps = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%fract%kiecps = sibr_vd(v)
              enddo
           enddo
       endif

   case (148)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%fract%d13cca = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%fract%d13cca = sibr_vd(v)
              enddo
           enddo
       endif

   case (149)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%fract%rcassim = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%fract%rcassim = sibr_vd(v)
              enddo
           enddo
       endif

   case (150)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%fract%c13assim = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%fract%c13assim = sibr_vd(v)
              enddo
           enddo
       endif

   case (151)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%fract%c13assimd = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%fract%c13assimd = sibr_vd(v)
              enddo
           enddo
       endif

   case (152)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%resp_nvegc13 = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%poollt%resp_nvegc13 = sibr_vd(v)
              enddo
           enddo
       endif

   case (153)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%resp_leafc13 = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%poollt%resp_leafc13 = sibr_vd(v)
              enddo
           enddo
       endif

   case (154)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%resp_rootc13 = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%poollt%resp_rootc13 = sibr_vd(v)
              enddo
           enddo
       endif

   case (155)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%resp_growc13 = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%poollt%resp_growc13 = sibr_vd(v)
              enddo
           enddo
       endif

   case (156)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%resp_mntnc13 = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%poollt%resp_mntnc13 = sibr_vd(v)
              enddo
           enddo
       endif

   case (157)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%resp_autoc13 = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%poollt%resp_autoc13 = sibr_vd(v)
              enddo
           enddo
       endif

   case (158)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%pooldt%resp_soilc13 = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%pooldt%resp_soilc13 = sibr_vd(v)
              enddo
           enddo
       endif

   case (159)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%pooldt%resp_hetc13 = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%pooldt%resp_hetc13 = sibr_vd(v)
              enddo
           enddo
       endif

   case (160)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%resp_grzc13 = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%poollt%resp_grzc13 = sibr_vd(v)
              enddo
           enddo
       endif

   case (161)
          if (var_exist) then
              allocate(dvar3d(nsibf,nlu,npoolpft))
              status = nf90_get_var(ncid, varid, dvar3d)
              if (status .ne. nf90_noerr) then
                 print*,'Error reading restart variable: ', &
                     trim(sibr_vname(v))
                 stop
              endif

              do i=1, subcount
                 if (nsibr_flag) then
                     gref = subset(i)
                 else
                     gref = i
                 endif

                 do l=1,sib%g(i)%g_nlu
                     sib%g(i)%l(l)%poollt%rcpoolpft(:) = dvar3d(gref,l,:)
                 enddo
              enddo
              deallocate(dvar3d)
          else
              do i=1, subcount
                 do l=1,sib%g(i)%g_nlu
                    sib%g(i)%l(l)%poollt%rcpoolpft(:) = sibr_vd(v)
                 enddo
              enddo
          endif

    case (162)
       if (var_exist) then
          allocate(dvar4d(nsibf,nlu,npoolpft,nsoil))
          status = nf90_get_var(ncid, varid, dvar4d)
          if (status .ne. nf90_noerr) then
               print*,'Error reading poolpft_dgain.'
               stop
          endif

          do i=1, subcount
             if (nsibr_flag) then
                 gref = subset(i)
             else
                 gref = i
             endif

              do l=1,sib%g(i)%g_nlu
                 do p=1,npoolpft
                    sib%g(i)%l(l)%poollt%rcpoolpft_lay(p,:) = &
                          dvar4d(gref,l,p,:)
                  enddo
               enddo
           enddo
           deallocate(dvar4d)

       else  !var does not exist
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%rcpoolpft_lay(:,:) = sibr_vd(v)
               enddo
           enddo
       endif

   case (163)
          if (var_exist) then
              allocate(dvar3d(nsibf,nlu,npoollu))
              status = nf90_get_var(ncid, varid, dvar3d)
              if (status .ne. nf90_noerr) then
                 print*,'Error reading restart variable: ', &
                     trim(sibr_vname(v))
                 stop
              endif

              do i=1, subcount
                 if (nsibr_flag) then
                     gref = subset(i)
                 else
                     gref = i
                 endif

                 do l=1,sib%g(i)%g_nlu
                     sib%g(i)%l(l)%pooldt%rcpoollu(:) = dvar3d(gref,l,:)
                 enddo
              enddo
              deallocate(dvar3d)
          else
              do i=1, subcount
                 do l=1,sib%g(i)%g_nlu
                    sib%g(i)%l(l)%pooldt%rcpoollu(:) = sibr_vd(v)
                 enddo
              enddo
          endif

    case (164)
       if (var_exist) then
           allocate(dvar4d(nsibf,nlu,npoollu,nsoil))
           status = nf90_get_var(ncid, varid, dvar4d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading poollu_lay.'
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

               do l=1,sib%g(i)%g_nlu
                  do p=1,npoollu
                     sib%g(i)%l(l)%pooldt%rcpoollu_lay(p,:) = &
                              dvar4d(gref,l,p,:)
                   enddo
               enddo
            enddo
            deallocate(dvar4d)

       else  !var does not exist
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%pooldt%rcpoollu_lay(:,:) = sibr_vd(v)
              enddo
           enddo
       endif

   case (165)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%fract%rcpoolfirec13 = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%fract%rcpoolfirec13 = sibr_vd(v)
              enddo
           enddo
       endif

   case (166)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%fract%rcassimfac = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%fract%rcassimfac = sibr_vd(v)
              enddo
           enddo
       endif

   case (167)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%fract%poolemistotC = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%fract%poolemistotC = sibr_vd(v)
              enddo
           enddo
       endif

   case (168)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%fract%poolemisc13 = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%fract%poolemisc13 = sibr_vd(v)
              enddo
           enddo
       endif

   case (170)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%fract%c12ca = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%fract%c12ca = sibr_vd(v)
              enddo
           enddo
       endif

   case (171)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%fract%c13ca = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%fract%c13ca = sibr_vd(v)
              enddo
           enddo
       endif

   case (172)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%fract%c12assim = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%fract%c12assim = sibr_vd(v)
              enddo
           enddo
       endif

   case (173)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%fract%c13resptot = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%fract%c13resptot = sibr_vd(v)
              enddo
           enddo
       endif

   case (174)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%fract%c12resptot = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%fract%c12resptot = sibr_vd(v)
              enddo
           enddo
       endif

   case (175)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%resp_hrvstc13 = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%poollt%resp_hrvstc13 = sibr_vd(v)
              enddo
           enddo
       endif

    case (182)
       if (.not. var_exist) then
          status = nf90_inq_varid(ncid, 'loss_transf_lu_lay', varid)
          var_exist = (status .eq. nf90_noerr)
       endif

       if (var_exist) then
           allocate(dvar4d(nsibf,nlu,npoollu,nsoil))
           status = nf90_get_var(ncid, varid, dvar4d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading loss_transf_lu_lay.'
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

               do l=1,sib%g(i)%g_nlu
                  do p=1,npoollu
                     sib%g(i)%l(l)%pooldt%loss_transf_lay(p,:) = &
                              dvar4d(gref,l,p,:)
                   enddo
               enddo
            enddo
            deallocate(dvar4d)

       else  !var does not exist
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%pooldt%loss_transf_lay(:,:) = sibr_vd(v)
              enddo
           enddo
       endif

    case (188)
       if (.not. var_exist) then
          status = nf90_inq_varid(ncid, 'loss_raddecay_lay', varid)
          var_exist = (status .eq. nf90_noerr)
       endif

       if (var_exist) then
          allocate(dvar4d(nsibf,nlu,npoolpft,nsoil))
          status = nf90_get_var(ncid, varid, dvar4d)
          if (status .ne. nf90_noerr) then
               print*,'Error reading loss_raddecay_lay.'
               stop
          endif

          do i=1, subcount
             if (nsibr_flag) then
                 gref = subset(i)
             else
                 gref = i
             endif

             do l=1,sib%g(i)%g_nlu
                do p=1,npoolpft
                   sib%g(i)%l(l)%poollt%loss_raddecay_lay(p,:) = dvar4d(gref,l,p,:)
                enddo
             enddo
           enddo
           deallocate(dvar4d)

       else  !var does not exist
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%loss_raddecay_lay(:,:) = sibr_vd(v)
               enddo
           enddo
       endif

    case (189)
       if (.not. var_exist) then
          status = nf90_inq_varid(ncid, 'loss_raddecay_lay', varid)
          var_exist = (status .eq. nf90_noerr)
       endif

       if (var_exist) then
           allocate(dvar4d(nsibf,nlu,npoollu,nsoil))
           status = nf90_get_var(ncid, varid, dvar4d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading loss_raddecay_lay.'
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

               do l=1,sib%g(i)%g_nlu
                  do p=1,npoollu
                     sib%g(i)%l(l)%pooldt%loss_raddecay_lay(p,:) = &
                              dvar4d(gref,l,p,:)
                   enddo
               enddo
            enddo
            deallocate(dvar4d)

       else  !var does not exist
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%pooldt%loss_raddecay_lay(:,:) = sibr_vd(v)
              enddo
           enddo
       endif

    case (190)
       if (.not. var_exist) then
          status = nf90_inq_varid(ncid, 'loss_grzc14', varid)
          var_exist = (status .eq. nf90_noerr)
       endif

       if (var_exist) then
          allocate(dvar3d(nsibf,nlu,npoolcanc14))
          status = nf90_get_var(ncid, varid, dvar3d)
          if (status .ne. nf90_noerr) then
               print*,'Error reading loss_grzc14.'
               stop
          endif

          do i=1, subcount
             if (nsibr_flag) then
                 gref = subset(i)
             else
                 gref = i
             endif

             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%poollt%loss_grzc14(:) = dvar3d(gref,l,:)
             enddo
           enddo
           deallocate(dvar3d)

       else  !var does not exist
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%loss_grzc14(:) = sibr_vd(v)
               enddo
           enddo
       endif

   case (191)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%fract%rcpoolfirec14 = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%fract%rcpoolfirec14 = sibr_vd(v)
              enddo
           enddo
       endif

   case (192)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%fract%rcassimfacc14 = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%fract%rcassimfacc14 = sibr_vd(v)
              enddo
           enddo
       endif

   case (193)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%fract%poolemisc14 = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%fract%poolemisc14 = sibr_vd(v)
              enddo
           enddo
       endif

   case (194)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%fract%c14alpha = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%fract%c14alpha = sibr_vd(v)
              enddo
           enddo
       endif

   case (195)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%fract%d14cca = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%fract%d14cca = sibr_vd(v)
              enddo
           enddo
       endif

   case (196)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%fract%rcassimc14 = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%fract%rcassimc14 = sibr_vd(v)
              enddo
           enddo
       endif

   case (197)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%fract%c14assim = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%fract%c14assim = sibr_vd(v)
              enddo
           enddo
       endif

   case (198)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%fract%c14assimd = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%fract%c14assimd = sibr_vd(v)
              enddo
           enddo
       endif

   case (199)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%resp_nvegc14 = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%poollt%resp_nvegc14 = sibr_vd(v)
              enddo
           enddo
       endif

   case (200)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%resp_leafc14 = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%poollt%resp_leafc14 = sibr_vd(v)
              enddo
           enddo
       endif

   case (201)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%resp_rootc14 = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%poollt%resp_rootc14 = sibr_vd(v)
              enddo
           enddo
       endif

   case (202)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%resp_growc14 = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%poollt%resp_growc14 = sibr_vd(v)
              enddo
           enddo
       endif

   case (203)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%resp_mntnc14 = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%poollt%resp_mntnc14 = sibr_vd(v)
              enddo
           enddo
       endif

   case (204)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%resp_autoc14 = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%poollt%resp_autoc14 = sibr_vd(v)
              enddo
           enddo
       endif

   case (205)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%pooldt%resp_soilc14 = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%pooldt%resp_soilc14 = sibr_vd(v)
              enddo
           enddo
       endif

   case (206)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%pooldt%resp_hetc14 = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%pooldt%resp_hetc14 = sibr_vd(v)
              enddo
           enddo
       endif

   case (207)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%resp_grzc14 = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%poollt%resp_grzc14 = sibr_vd(v)
              enddo
           enddo
       endif

   case (208)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%resp_hrvstc14 = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%poollt%resp_hrvstc14 = sibr_vd(v)
              enddo
           enddo
       endif

   case (209)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%fract%c14ca = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%fract%c14ca = sibr_vd(v)
              enddo
           enddo
       endif

   case (210)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%fract%c14resptot = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%fract%c14resptot = sibr_vd(v)
              enddo
           enddo
       endif

   case (211)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(dvar2d(nsibf,nlu))
          status = nf90_get_var(ncid, varid, dvar2d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%phent%phenave_assimpot = dvar2d(gref,l)
             enddo
          enddo
          deallocate(dvar2d)
       endif

   case (213)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(dvar2d(nsibf,nlu))
          status = nf90_get_var(ncid, varid, dvar2d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%phent%phenave_prpot = dvar2d(gref,l)
             enddo
          enddo
          deallocate(dvar2d)
       endif

   case (222)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(dvar2d(nsibf,nlu))
          status = nf90_get_var(ncid, varid, dvar2d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%phent%phens_dayl = dvar2d(gref,l)
             enddo
          enddo
          deallocate(dvar2d)
       endif

   case (223)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(dvar2d(nsibf,nlu))
          status = nf90_get_var(ncid, varid, dvar2d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%phent%phens_grw = dvar2d(gref,l)
             enddo
          enddo
          deallocate(dvar2d)
       endif

   case (224)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(dvar2d(nsibf,nlu))
          status = nf90_get_var(ncid, varid, dvar2d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%phent%phenave_wac = dvar2d(gref,l)
             enddo
          enddo
          deallocate(dvar2d)
       endif

   case (225)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(dvar2d(nsibf,nlu))
          status = nf90_get_var(ncid, varid, dvar2d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%phent%phens_wx = dvar2d(gref,l)
             enddo
          enddo
          deallocate(dvar2d)
       endif

   case (212)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(ivar2d(nsibf,nlu))
          status = nf90_get_var(ncid, varid, ivar2d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%phent%phenflag_assimlow = (ivar2d(gref,l) .ne. izero)
             enddo
          enddo
          deallocate(ivar2d)
       endif

   case (214)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(ivar2d(nsibf,nlu))
          status = nf90_get_var(ncid, varid, ivar2d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%phent%phenflag_precip = (ivar2d(gref,l) .ne. izero)
             enddo
          enddo
          deallocate(ivar2d)
       endif

   case (215)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(ivar2d(nsibf,nlu))
          status = nf90_get_var(ncid, varid, ivar2d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%phent%phenflag_moist = (ivar2d(gref,l) .ne. izero)
             enddo
          enddo
          deallocate(ivar2d)
       endif

   case (216)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(ivar2d(nsibf,nlu))
          status = nf90_get_var(ncid, varid, ivar2d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%phent%phenflag_temp = (ivar2d(gref,l) .ne. izero)
             enddo
          enddo
          deallocate(ivar2d)
       endif

   case (217)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(ivar2d(nsibf,nlu))
          status = nf90_get_var(ncid, varid, ivar2d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%phent%phenflag_daylen = (ivar2d(gref,l) .ne. izero)
             enddo
          enddo
          deallocate(ivar2d)
       endif

   case (218)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(ivar2d(nsibf,nlu))
          status = nf90_get_var(ncid, varid, ivar2d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%phent%phenflag_gsspass = (ivar2d(gref,l) .ne. izero)
             enddo
          enddo
          deallocate(ivar2d)
       endif

   case (219)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(ivar2d(nsibf,nlu))
          status = nf90_get_var(ncid, varid, ivar2d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%phent%nd_dormant = ivar2d(gref,l)
             enddo
          enddo
          deallocate(ivar2d)
       endif

   case (220)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(ivar2d(nsibf,nlu))
          status = nf90_get_var(ncid, varid, ivar2d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%phent%nd_gs = ivar2d(gref,l)
             enddo
          enddo
          deallocate(ivar2d)
       endif

   case (226)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(ivar2d(nsibf,nlu))
          status = nf90_get_var(ncid, varid, ivar2d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%phent%ipd = ivar2d(gref,l)
             enddo
          enddo
          deallocate(ivar2d)
       endif

   case (221)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(ivar3d(nsibf,nlu,npstgmax))
          status = nf90_get_var(ncid, varid, ivar3d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%phent%nd_stg(:) = &
                     ivar3d(gref,l,1:size(sib%g(i)%l(l)%phent%nd_stg))
             enddo
          enddo
          deallocate(ivar3d)
       endif

   case (230)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(dvar3d(nsibf,nlu,npoolpft))
          status = nf90_get_var(ncid, varid, dvar3d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%poollt%poolpftp(:) = dvar3d(gref,l,:)
             enddo
          enddo
          deallocate(dvar3d)
       endif

   case (231)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(dvar3d(nsibf,nlu,npoollu))
          status = nf90_get_var(ncid, varid, dvar3d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%pooldt%poollup(:) = dvar3d(gref,l,:)
             enddo
          enddo
          deallocate(dvar3d)
       endif

   case (232)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(dvar3d(nsibf,nlu,npoolpft))
          status = nf90_get_var(ncid, varid, dvar3d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%poollt%alloc_phen(:) = dvar3d(gref,l,:)
             enddo
          enddo
          deallocate(dvar3d)
       endif

   case (233)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(ivar2d(nsibf,nlu))
          status = nf90_get_var(ncid, varid, ivar2d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%poollt%aadj_moist = (ivar2d(gref,l) .ne. izero)
             enddo
          enddo
          deallocate(ivar2d)
       endif

   case (234)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(ivar2d(nsibf,nlu))
          status = nf90_get_var(ncid, varid, ivar2d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%poollt%aadj_temp = (ivar2d(gref,l) .ne. izero)
             enddo
          enddo
          deallocate(ivar2d)
       endif

   case (235)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(dvar2d(nsibf,nlu))
          status = nf90_get_var(ncid, varid, dvar2d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%poollt%tfl_pstage = dvar2d(gref,l)
             enddo
          enddo
          deallocate(dvar2d)
       endif

   case (236)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(dvar2d(nsibf,nlu))
          status = nf90_get_var(ncid, varid, dvar2d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%poollt%nd_fire = dvar2d(gref,l)
             enddo
          enddo
          deallocate(dvar2d)
       endif

   case (237)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(dvar2d(nsibf,nlu))
          status = nf90_get_var(ncid, varid, dvar2d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%poollt%nd_grz = dvar2d(gref,l)
             enddo
          enddo
          deallocate(dvar2d)
       endif

   case (238)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(dvar2d(nsibf,nlu))
          status = nf90_get_var(ncid, varid, dvar2d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%pooldt%mhrt_soil_hot = dvar2d(gref,l)
             enddo
          enddo
          deallocate(dvar2d)
       endif

   case (239)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(dvar3d(nsibf,nlu,nsoil))
          status = nf90_get_var(ncid, varid, dvar3d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%pooldt%resp_soilnr_lay(:) = dvar3d(gref,l,:)
             enddo
          enddo
          deallocate(dvar3d)
       endif

   case (240)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(dvar2d(nsibf,nlu))
          status = nf90_get_var(ncid, varid, dvar2d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%hydrovt%pawfrw = dvar2d(gref,l)
             enddo
          enddo
          deallocate(dvar2d)
       endif

   case (241)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(dvar2d(nsibf,nlu))
          status = nf90_get_var(ncid, varid, dvar2d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%hydrovt%pawftop = dvar2d(gref,l)
             enddo
          enddo
          deallocate(dvar2d)
       endif

   case (242)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(dvar2d(nsibf,nlu))
          status = nf90_get_var(ncid, varid, dvar2d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%hydrovt%pawfzw = dvar2d(gref,l)
             enddo
          enddo
          deallocate(dvar2d)
       endif

   case (243)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(dvar2d(nsibf,nlu))
          status = nf90_get_var(ncid, varid, dvar2d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%hydrovt%tawfrw = dvar2d(gref,l)
             enddo
          enddo
          deallocate(dvar2d)
       endif

   case (244)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(dvar2d(nsibf,nlu))
          status = nf90_get_var(ncid, varid, dvar2d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%hydrovt%tawftop = dvar2d(gref,l)
             enddo
          enddo
          deallocate(dvar2d)
       endif

   case (245)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(dvar2d(nsibf,nlu))
          status = nf90_get_var(ncid, varid, dvar2d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%hydrovt%tawfzw = dvar2d(gref,l)
             enddo
          enddo
          deallocate(dvar2d)
       endif

   case (246)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(dvar2d(nsibf,nlu))
          status = nf90_get_var(ncid, varid, dvar2d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%hydrost%snow_gmass = dvar2d(gref,l)
             enddo
          enddo
          deallocate(dvar2d)
       endif

   case (247)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(dvar2d(nsibf,nlu))
          status = nf90_get_var(ncid, varid, dvar2d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%hydrost%snow_gdepth = dvar2d(gref,l)
             enddo
          enddo
          deallocate(dvar2d)
       endif

   case (248)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(dvar2d(nsibf,nlu))
          status = nf90_get_var(ncid, varid, dvar2d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%hydrost%snow_cvfc = dvar2d(gref,l)
             enddo
          enddo
          deallocate(dvar2d)
       endif

   case (249)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(dvar2d(nsibf,nlu))
          status = nf90_get_var(ncid, varid, dvar2d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%hydrost%wetfracc = dvar2d(gref,l)
             enddo
          enddo
          deallocate(dvar2d)
       endif

   case (250)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(dvar2d(nsibf,nlu))
          status = nf90_get_var(ncid, varid, dvar2d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%hydrost%wetfracg = dvar2d(gref,l)
             enddo
          enddo
          deallocate(dvar2d)
       endif

   case (255)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(dvar2d(nsibf,nlu))
          status = nf90_get_var(ncid, varid, dvar2d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%vegt%lai = dvar2d(gref,l)
             enddo
          enddo
          deallocate(dvar2d)
       endif

   case (256)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(dvar2d(nsibf,nlu))
          status = nf90_get_var(ncid, varid, dvar2d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%vegt%lait = dvar2d(gref,l)
             enddo
          enddo
          deallocate(dvar2d)
       endif

   case (257)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(dvar2d(nsibf,nlu))
          status = nf90_get_var(ncid, varid, dvar2d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%vegt%fpar = dvar2d(gref,l)
             enddo
          enddo
          deallocate(dvar2d)
       endif

   case (258)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(dvar2d(nsibf,nlu))
          status = nf90_get_var(ncid, varid, dvar2d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%vegt%vcover = dvar2d(gref,l)
             enddo
          enddo
          deallocate(dvar2d)
       endif

   case (259)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(dvar2d(nsibf,nlu))
          status = nf90_get_var(ncid, varid, dvar2d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%vegt%green = dvar2d(gref,l)
             enddo
          enddo
          deallocate(dvar2d)
       endif

   case (260)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(dvar2d(nsibf,nlu))
          status = nf90_get_var(ncid, varid, dvar2d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%vegt%z0 = dvar2d(gref,l)
             enddo
          enddo
          deallocate(dvar2d)
       endif

   case (261)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(dvar2d(nsibf,nlu))
          status = nf90_get_var(ncid, varid, dvar2d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%vegt%z0d = dvar2d(gref,l)
             enddo
          enddo
          deallocate(dvar2d)
       endif

   case (262)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(dvar2d(nsibf,nlu))
          status = nf90_get_var(ncid, varid, dvar2d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%vegt%zp_dispd = dvar2d(gref,l)
             enddo
          enddo
          deallocate(dvar2d)
       endif

   case (263)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(dvar2d(nsibf,nlu))
          status = nf90_get_var(ncid, varid, dvar2d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%vegt%zpd_adj = dvar2d(gref,l)
             enddo
          enddo
          deallocate(dvar2d)
       endif

   case (264)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(dvar2d(nsibf,nlu))
          status = nf90_get_var(ncid, varid, dvar2d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%vegt%zztemp = dvar2d(gref,l)
             enddo
          enddo
          deallocate(dvar2d)
       endif

   case (265)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(dvar2d(nsibf,nlu))
          status = nf90_get_var(ncid, varid, dvar2d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%vegt%zzwind = dvar2d(gref,l)
             enddo
          enddo
          deallocate(dvar2d)
       endif

   case (266)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(dvar2d(nsibf,nlu))
          status = nf90_get_var(ncid, varid, dvar2d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%vegt%cc1 = dvar2d(gref,l)
             enddo
          enddo
          deallocate(dvar2d)
       endif

   case (267)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(dvar2d(nsibf,nlu))
          status = nf90_get_var(ncid, varid, dvar2d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%vegt%cc2 = dvar2d(gref,l)
             enddo
          enddo
          deallocate(dvar2d)
       endif

   case (268)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(dvar2d(nsibf,nlu))
          status = nf90_get_var(ncid, varid, dvar2d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%vegt%gmudmu = dvar2d(gref,l)
             enddo
          enddo
          deallocate(dvar2d)
       endif

   case (269)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(dvar2d(nsibf,nlu))
          status = nf90_get_var(ncid, varid, dvar2d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%vegt%vmax = dvar2d(gref,l)
             enddo
          enddo
          deallocate(dvar2d)
       endif

   case (275)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(dvar2d(nsibf,nlu))
          status = nf90_get_var(ncid, varid, dvar2d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%co2t%rstfac(1) = dvar2d(gref,l)
             enddo
          enddo
          deallocate(dvar2d)
       endif

   case (276)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(dvar2d(nsibf,nlu))
          status = nf90_get_var(ncid, varid, dvar2d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%co2t%rstfac(2) = dvar2d(gref,l)
             enddo
          enddo
          deallocate(dvar2d)
       endif

   case (277)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(dvar2d(nsibf,nlu))
          status = nf90_get_var(ncid, varid, dvar2d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%co2t%rstfac(3) = dvar2d(gref,l)
             enddo
          enddo
          deallocate(dvar2d)
       endif

   case (278)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(dvar2d(nsibf,nlu))
          status = nf90_get_var(ncid, varid, dvar2d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%co2t%rstfac(4) = dvar2d(gref,l)
             enddo
          enddo
          deallocate(dvar2d)
       endif

   case (279)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(dvar3d(nsibf,nlu,nsoil))
          status = nf90_get_var(ncid, varid, dvar3d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%cost%cos_s(:) = dvar3d(gref,l,:)
             enddo
          enddo
          deallocate(dvar3d)
       endif

   case (280:337)
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: V3 forcing-bracket restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with the V3 executable.'
          stop
       endif
       if (var_exist) then
          allocate(dvar(nsibf))
          status = nf90_get_var(ncid,varid,dvar)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             select case (sibr_vref(v))
             case (280)
                sib%g(i)%gprogt%tm1 = dvar(gref)
             case (281)
                sib%g(i)%gprogt%tm2 = dvar(gref)
             case (282)
                sib%g(i)%gprogt%sh1 = dvar(gref)
             case (283)
                sib%g(i)%gprogt%sh2 = dvar(gref)
             case (284)
                sib%g(i)%gprogt%ps1 = dvar(gref)
             case (285)
                sib%g(i)%gprogt%ps2 = dvar(gref)
             case (286)
                sib%g(i)%gprogt%spdm1 = dvar(gref)
             case (287)
                sib%g(i)%gprogt%spdm2 = dvar(gref)
             case (288)
                sib%g(i)%gprogt%lspr1 = dvar(gref)
             case (289)
                sib%g(i)%gprogt%lspr2 = dvar(gref)
             case (290)
                sib%g(i)%gprogt%cupr1 = dvar(gref)
             case (291)
                sib%g(i)%gprogt%cupr2 = dvar(gref)
             case (292)
                sib%g(i)%gprogt%dlwbot1 = dvar(gref)
             case (293)
                sib%g(i)%gprogt%dlwbot2 = dvar(gref)
             case (294)
                sib%g(i)%gprogt%sw_dwn1 = dvar(gref)
             case (295)
                sib%g(i)%gprogt%sw_dwn2 = dvar(gref)
             case (296)
                sib%g(i)%gprogt%firec1 = dvar(gref)
             case (297)
                sib%g(i)%gprogt%firec2 = dvar(gref)
             case (298)
                sib%g(i)%gprogt%fireco21 = dvar(gref)
             case (299)
                sib%g(i)%gprogt%fireco22 = dvar(gref)
             case (300)
                sib%g(i)%gprogt%ba_dbg1 = dvar(gref)
             case (301)
                sib%g(i)%gprogt%ba_dbg2 = dvar(gref)
             case (302)
                sib%g(i)%gprogt%ba_en21 = dvar(gref)
             case (303)
                sib%g(i)%gprogt%ba_en22 = dvar(gref)
             case (304)
                sib%g(i)%gprogt%ba_en31 = dvar(gref)
             case (305)
                sib%g(i)%gprogt%ba_en32 = dvar(gref)
             case (306)
                sib%g(i)%gprogt%ba_dnf1 = dvar(gref)
             case (307)
                sib%g(i)%gprogt%ba_dnf2 = dvar(gref)
             case (308)
                sib%g(i)%gprogt%ba_eb11 = dvar(gref)
             case (309)
                sib%g(i)%gprogt%ba_eb12 = dvar(gref)
             case (310)
                sib%g(i)%gprogt%ba_eb21 = dvar(gref)
             case (311)
                sib%g(i)%gprogt%ba_eb22 = dvar(gref)
             case (312)
                sib%g(i)%gprogt%ba_db11 = dvar(gref)
             case (313)
                sib%g(i)%gprogt%ba_db12 = dvar(gref)
             case (314)
                sib%g(i)%gprogt%ba_db21 = dvar(gref)
             case (315)
                sib%g(i)%gprogt%ba_db22 = dvar(gref)
             case (316)
                sib%g(i)%gprogt%ba_db31 = dvar(gref)
             case (317)
                sib%g(i)%gprogt%ba_db32 = dvar(gref)
             case (318)
                sib%g(i)%gprogt%ba_shb1 = dvar(gref)
             case (319)
                sib%g(i)%gprogt%ba_shb2 = dvar(gref)
             case (320)
                sib%g(i)%gprogt%ba_sha1 = dvar(gref)
             case (321)
                sib%g(i)%gprogt%ba_sha2 = dvar(gref)
             case (322)
                sib%g(i)%gprogt%ba_c3a1 = dvar(gref)
             case (323)
                sib%g(i)%gprogt%ba_c3a2 = dvar(gref)
             case (324)
                sib%g(i)%gprogt%ba_c3g1 = dvar(gref)
             case (325)
                sib%g(i)%gprogt%ba_c3g2 = dvar(gref)
             case (326)
                sib%g(i)%gprogt%ba_c4g1 = dvar(gref)
             case (327)
                sib%g(i)%gprogt%ba_c4g2 = dvar(gref)
             case (328)
                sib%g(i)%gprogt%ba_c3c1 = dvar(gref)
             case (329)
                sib%g(i)%gprogt%ba_c3c2 = dvar(gref)
             case (330)
                sib%g(i)%gprogt%ba_c4c1 = dvar(gref)
             case (331)
                sib%g(i)%gprogt%ba_c4c2 = dvar(gref)
             case (332)
                sib%g(i)%gprogt%ba_mze1 = dvar(gref)
             case (333)
                sib%g(i)%gprogt%ba_mze2 = dvar(gref)
             case (334)
                sib%g(i)%gprogt%ba_soy1 = dvar(gref)
             case (335)
                sib%g(i)%gprogt%ba_soy2 = dvar(gref)
             case (336)
                sib%g(i)%gprogt%ba_wwt1 = dvar(gref)
             case (337)
                sib%g(i)%gprogt%ba_wwt2 = dvar(gref)
             end select
          enddo
          deallocate(dvar)
       endif

   case (338:339)
       ! Exact-continuation diagnostic inputs are required whenever
       ! ic_mode='restart' is selected.
       if (ic_restart .and. (.not. var_exist)) then
          print*,'ERROR: diagnostic-continuation restart variable is missing: ', &
                 trim(sibr_vname(v))
          print*,'       Regenerate the checkpoint with a continuation-capable executable.'
          stop
       endif
       if (var_exist) then
          allocate(dvar2d(nsibf,nlu))
          status = nf90_get_var(ncid,varid,dvar2d)
          if (status .ne. nf90_noerr) then
             print*,'Error reading restart variable: ',trim(sibr_vname(v))
             stop
          endif
          do i=1,subcount
             if (nsibr_flag) then
                gref = subset(i)
             else
                gref = i
             endif
             do l=1,sib%g(i)%g_nlu
                select case (sibr_vref(v))
                case (338)
                   sib%g(i)%l(l)%cast%vpd = dvar2d(gref,l)
                case (339)
                   sib%g(i)%l(l)%fluxt%ect = dvar2d(gref,l)
                end select
             enddo
          enddo
          deallocate(dvar2d)
       endif

   end select
enddo  !v=1,sibr_nvar

! Source-cohort provenance is used only by exact transient continuation.
! Fresh requib ICs initialize this storage to zero via sibtype_init.
if (file_exist .and. ic_restart) then
   status = nf90_get_att(ncid,nf90_global,'pftswap_schema_version',swschema)
   if (status .eq. nf90_noerr) then
      if (swschema .ne. 1) stop 'Unsupported pftswap checkpoint schema'
      status = nf90_inq_varid(ncid,'pftswap_phase',swphaseid)
      if (status .ne. nf90_noerr) stop 'Missing pftswap_phase in checkpoint'
      status = nf90_inq_varid(ncid,'pftswap_source_pftref',swsourceid)
      if (status .ne. nf90_noerr) stop 'Missing pftswap_source_pftref'
      status = nf90_inq_varid(ncid,'pftswap_legacy_live_lay',swliveid)
      if (status .ne. nf90_noerr) stop 'Missing pftswap_legacy_live_lay'
      status = nf90_inq_varid(ncid,'pftswap_legacy_dead_lay',swdeadid)
      if (status .ne. nf90_noerr) stop 'Missing pftswap_legacy_dead_lay'
      allocate(swlivebuf(1,1,npoolpft,nsoil), swdeadbuf(1,1,npoollu,nsoil))
      do i=1,subcount
         if (nsibr_flag) then
            gref = subset(i)
         else
            gref = i
         endif
         do l=1,sib%g(i)%g_nlu
            status = nf90_get_var(ncid,swphaseid,swibuf, &
                       start=(/gref,l/),count=(/1,1/))
            if (status .ne. nf90_noerr) stop 'Reading pftswap_phase failed'
            sib%g(i)%l(l)%pftswap_legacy%phase = swibuf(1,1)
            status = nf90_get_var(ncid,swsourceid,swibuf, &
                       start=(/gref,l/),count=(/1,1/))
            if (status .ne. nf90_noerr) stop 'Reading pftswap_source failed'
            sib%g(i)%l(l)%pftswap_legacy%source_pftref = swibuf(1,1)
            status = nf90_get_var(ncid,swliveid,swlivebuf, &
                 start=(/gref,l,1,1/),count=(/1,1,npoolpft,nsoil/))
            if (status .ne. nf90_noerr) stop 'Reading pftswap_live failed'
            status = nf90_get_var(ncid,swdeadid,swdeadbuf, &
                 start=(/gref,l,1,1/),count=(/1,1,npoollu,nsoil/))
            if (status .ne. nf90_noerr) stop 'Reading pftswap_dead failed'
            sib%g(i)%l(l)%pftswap_legacy%live_lay(:,:) = swlivebuf(1,1,:,:)
            sib%g(i)%l(l)%pftswap_legacy%dead_lay(:,:) = swdeadbuf(1,1,:,:)
            select case (sib%g(i)%l(l)%pftswap_legacy%phase)
            case (0)
               ! Fully vegetated pre-transition checkpoint.
            case (1)
               if (sib%g(i)%l(l)%ipft .ne. 1) &
                  error stop 'PFTSWAP open checkpoint must use bare-ground PFT 1'
               if (sib%g(i)%l(l)%pftswap_legacy%source_pftref .ne. 11) &
                  error stop 'PFTSWAP open checkpoint source must be shrub PFT 11'
               if (maxval(abs(sib%g(i)%l(l)%poollt%poolpft_lay)) .gt. 1.d-12) &
                  error stop 'PFTSWAP open checkpoint has active live C'
            case default
               error stop 'Unsupported PFTSWAP transition phase'
            end select
         enddo
      enddo
      deallocate(swlivebuf,swdeadbuf)
   else
      print*,'PFTSWAP: earlier checkpoint has no legacy schema; initial state retained.'
   endif
endif

status = nf90_close(ncid)

end subroutine restart_read


