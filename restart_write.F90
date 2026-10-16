#include "nc_util.h"

subroutine restart_write( equib )
    !---------------------------------------------
    ! This subroutine creates a netcdf restart file
    ! for all prognostic variables.
    !---------------------------------------------

    use kinds
    use netcdf

    use module_io
    use module_pftinfo, only: &
        clen, pft_name, pft_ref
    use module_sibconst, only: &
        slen, npft, nlu, &
        npoolpft, npoollu, npoolcan, &
        nsoil, nsnow, ntot, npstgmax, &
        subcount, sublatsib, sublonsib, &
        subsitename, &
        npoolcanc13, npoolcanc14, &
        fire_switch, firegfed_switch
    use module_sib, only: sib

    implicit none
    
    !---input variables
    logical :: equib  !equilibrium pool IC file

    !---netcdf variables
    integer(i4) :: ncid, status
    integer(i4) :: varid
    
    integer(i4) :: nsibdid, npftdid, nludid, &
                   nppdid, nlpdid, ncpdid, &
                   nsdid, nsnowdid, ntotdid, nstgdid, &
                   nclendid, nslendid, &
                   ncpdidc13, ncpdidc14
    integer(i4) :: lonid, latid, &
                   pnameid, prefid, snameid, &
                   areaid, pftrefid

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
    
    !---local variables
    integer(i4) :: i,l,p,v,iatt
    integer(i4) :: swphaseid, swsourceid, swliveid, swdeadid
    integer(i4) :: qp_sifdim, qp_ngdim, qp_nludim
    integer(i4) :: qp_sifvar, qp_gvar, qp_luvar, qp_gidvar, qp_luidvar
    integer(i4) :: alloc_continue_var, alloc_moist_var, alloc_temp_var
    integer(i4), parameter :: nheld_diag=10
    integer(i4) :: held_dim, held_var, held_gain_var, held_loss_var
    real(r8), allocatable :: held_scalars(:,:,:), held_gain(:,:,:,:), &
          held_loss(:,:,:,:)
    real(r8), allocatable :: alloc_continue(:,:,:), alloc_moist_continue(:,:,:), &
         alloc_temp_continue(:,:,:)
    integer(i4) :: swibuf(1,1)
    real(r8), allocatable :: swlivebuf(:,:,:,:), swdeadbuf(:,:,:,:)
    character(len=16) :: fire_driver_type
    logical :: savesite
    logical, dimension(:), allocatable :: dovar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Determine if saving site names
if (subsitename(1) .ne. '') then
    savesite = .true.
else
    savesite = .false.
endif
    
! Write out the data to a file
allocate(dovar(sibr_nvar))
if (equib) then
    status = nf90_create( requib_filename, ncid=ncid, &
             cmode=or(nf90_clobber,nf90_64bit_offset))
    if (status .ne. nf90_noerr) then
        print('(a)'),'Error creating requib file: '
        print('(2a)'),'  ',requib_filename
        stop
    else
        print('(a)'),'Writing Equilibrium Restart File: '
        print('(2a)'),'  ', trim(requib_filename)
    endif
    dovar = sibreq_doref
else
    status = nf90_create( restart_filename, ncid=ncid, &
             cmode=or(nf90_clobber,nf90_64bit_offset))
    if (status .ne. nf90_noerr) then
        print('(a)'),'Error creating file: '
        print('(2a)'),'  ',restart_filename
        stop
    !else
    !    print('(a)'),'Writing Restart File: '
    !    print('(2a)'),'  ', trim(restart_filename)
    endif
    dovar = sibr_doref
endif

! Transient checkpoints carry exact-continuation forcing/control state below.
! requib files retain their historical fresh-initialization semantics.

! Store external-driver control state for exact continuation.
! These attributes are only consumed when ic_mode='restart'.  Integer(long)
! second counters are written as doubles because the classic 64-bit-offset
! restart format does not support netCDF int64, while IEEE double represents
! the relevant integer-second range exactly.
if (.not. equib) then
   CHECK( nf90_put_att(ncid, nf90_global, 'pftswap_schema_version', 1) )
   CHECK( nf90_put_att(ncid, nf90_global, 'driver_step', driver_step) )
   CHECK( nf90_put_att(ncid, nf90_global, 'driver_perday', driver_perday) )
   CHECK( nf90_put_att(ncid, nf90_global, 'driver_permon', driver_permon) )
   CHECK( nf90_put_att(ncid, nf90_global, 'driver_recnum', driver_recnum) )
   CHECK( nf90_put_att(ncid, nf90_global, 'driver_year', driver_year) )
   CHECK( nf90_put_att(ncid, nf90_global, 'driver_month', driver_month) )
   CHECK( nf90_put_att(ncid, nf90_global, 'driver_day', driver_day) )
   CHECK( nf90_put_att(ncid, nf90_global, 'driver_hour', driver_hour) )
   CHECK( nf90_put_att(ncid, nf90_global, 'driver_seccur', real(driver_seccur,r8)) )
   CHECK( nf90_put_att(ncid, nf90_global, 'driver_secnext', real(driver_secnext,r8)) )
   iatt = merge(1,0,driver_readf)
   CHECK( nf90_put_att(ncid, nf90_global, 'driver_readf', iatt) )
   iatt = merge(1,0,driver_switchf)
   CHECK( nf90_put_att(ncid, nf90_global, 'driver_switchf', iatt) )
   iatt = merge(1,0,driver_updatef)
   CHECK( nf90_put_att(ncid, nf90_global, 'driver_updatef', iatt) )

   if (fire_switch .and. firegfed_switch) then
      fire_driver_type = 'both'
   elseif (fire_switch) then
      fire_driver_type = 'standard'
   elseif (firegfed_switch) then
      fire_driver_type = 'ba'
   else
      fire_driver_type = 'none'
   endif
   CHECK( nf90_put_att(ncid, nf90_global, 'fire_driver_type', trim(fire_driver_type)) )

   if (fire_switch .or. firegfed_switch) then
      CHECK( nf90_put_att(ncid, nf90_global, 'fire_step', fire_step) )
      CHECK( nf90_put_att(ncid, nf90_global, 'fire_perday', fire_perday) )
      CHECK( nf90_put_att(ncid, nf90_global, 'fire_permon', fire_permon) )
      CHECK( nf90_put_att(ncid, nf90_global, 'fire_recnum', fire_recnum) )
      CHECK( nf90_put_att(ncid, nf90_global, 'fire_year', fire_year) )
      CHECK( nf90_put_att(ncid, nf90_global, 'fire_month', fire_month) )
      CHECK( nf90_put_att(ncid, nf90_global, 'fire_day', fire_day) )
      CHECK( nf90_put_att(ncid, nf90_global, 'fire_hour', fire_hour) )
      CHECK( nf90_put_att(ncid, nf90_global, 'fire_seccur', real(fire_seccur,r8)) )
      CHECK( nf90_put_att(ncid, nf90_global, 'fire_secnext', real(fire_secnext,r8)) )
   endif

   if (fire_switch) then
      iatt = merge(1,0,fire_readf)
      CHECK( nf90_put_att(ncid, nf90_global, 'fire_readf', iatt) )
      iatt = merge(1,0,fire_switchf)
      CHECK( nf90_put_att(ncid, nf90_global, 'fire_switchf', iatt) )
      iatt = merge(1,0,fire_updatef)
      CHECK( nf90_put_att(ncid, nf90_global, 'fire_updatef', iatt) )
   endif
   if (firegfed_switch) then
      iatt = merge(1,0,firegfed_readf)
      CHECK( nf90_put_att(ncid, nf90_global, 'firegfed_readf', iatt) )
      iatt = merge(1,0,firegfed_switchf)
      CHECK( nf90_put_att(ncid, nf90_global, 'firegfed_switchf', iatt) )
      iatt = merge(1,0,firegfed_updatef)
      CHECK( nf90_put_att(ncid, nf90_global, 'firegfed_updatef', iatt) )
   endif
endif


!...dimensions...
CHECK( nf90_def_dim( ncid, 'nsib', subcount, nsibdid ) )
CHECK( nf90_def_dim( ncid, 'nlu', nlu, nludid ) )
CHECK( nf90_def_dim( ncid, 'npft', npft, npftdid ) )
CHECK( nf90_def_dim( ncid, 'npoolcan', npoolcan, ncpdid ) )
CHECK( nf90_def_dim( ncid, 'npoolcanc13', npoolcanc13, ncpdidc13 ) )
CHECK( nf90_def_dim( ncid, 'npoolcanc14', npoolcanc14, ncpdidc14 ) )
CHECK( nf90_def_dim( ncid, 'npoolpft', npoolpft, nppdid ) )
CHECK( nf90_def_dim( ncid, 'npoollu', npoollu, nlpdid ) )
CHECK( nf90_def_dim( ncid, 'nsnow', nsnow, nsnowdid ) )
CHECK( nf90_def_dim( ncid, 'nsoil', nsoil, nsdid ) )
CHECK( nf90_def_dim( ncid, 'ntot', ntot, ntotdid ) )
CHECK( nf90_def_dim( ncid, 'npstgmax', npstgmax, nstgdid ) )
CHECK( nf90_def_dim( ncid, 'clen', clen, nclendid ) )
if (savesite) CHECK( nf90_def_dim( ncid, 'slen', slen, nslendid ) )

!...define and write info variables...
CHECK( nf90_def_var( ncid, 'lonsib', nf90_float, nsibdid, lonid ) )
CHECK( nf90_def_var( ncid, 'latsib', nf90_float, nsibdid, latid ) )
CHECK( nf90_def_var( ncid, 'pftnames', nf90_char, (/nclendid,npftdid/), pnameid))
CHECK( nf90_def_var( ncid, 'pftrefs', nf90_int, (/npftdid/), prefid ) )
if (savesite) then
    CHECK( nf90_def_var( ncid, 'site_names', nf90_char, (/nslendid,nsibdid/), snameid))
endif
CHECK( nf90_def_var( ncid, 'lu_area', nf90_float, (/nsibdid,nludid/), areaid ) )
CHECK( nf90_def_var( ncid, 'lu_pftref', nf90_int, (/nsibdid,nludid/), pftrefid ) )

CHECK( nf90_enddef(ncid) )
CHECK( nf90_put_var( ncid, lonid, sublonsib ) )
CHECK( nf90_put_var( ncid, latid, sublatsib ) )
CHECK( nf90_put_var( ncid, pnameid, pft_name ) ) 
CHECK( nf90_put_var( ncid, prefid, pft_ref ) )
if (savesite) CHECK( nf90_put_var( ncid, snameid, subsitename ) )
! Write the current runtime vegetation configuration, not the original
! sib_vs-derived metadata.  This makes a transient vegetation/PFT change
! restart-authoritative while leaving sib_vs as the fresh-start specification.
allocate(tvar2d(subcount,nlu))
allocate(ivar2d(subcount,nlu))
tvar2d = rzero
ivar2d = izero
do i=1,subcount
   do l=1,sib%g(i)%g_nlu
      tvar2d(i,l) = sib%g(i)%l(l)%larea
      ivar2d(i,l) = sib%g(i)%l(l)%ipft
   enddo
enddo
CHECK( nf90_put_var( ncid, areaid, tvar2d ) )
CHECK( nf90_put_var( ncid, pftrefid, ivar2d ) )
deallocate(tvar2d,ivar2d)

! Goal 1 Step 1B: explicit source-provenance restart fields.
! Write these even when disabled so a 2023 checkpoint is ready for
! future controlled testing.  No active carbon pools are modified.
if (.not. equib) then
   CHECK( nf90_redef(ncid) )
   CHECK( nf90_def_var(ncid, 'pftswap_phase', nf90_int, &
         (/nsibdid,nludid/), swphaseid) )
   CHECK( nf90_def_var(ncid, 'pftswap_source_pftref', nf90_int, &
         (/nsibdid,nludid/), swsourceid) )
   CHECK( nf90_def_var(ncid, 'pftswap_legacy_live_lay', nf90_double, &
         (/nsibdid,nludid,nppdid,nsdid/), swliveid) )
   CHECK( nf90_def_var(ncid, 'pftswap_legacy_dead_lay', nf90_double, &
         (/nsibdid,nludid,nlpdid,nsdid/), swdeadid) )
   CHECK( nf90_enddef(ncid) )
   allocate(swlivebuf(1,1,npoolpft,nsoil), swdeadbuf(1,1,npoollu,nsoil))
   do i=1,subcount
      do l=1,nlu
         swibuf = 0
         swlivebuf = dzero
         swdeadbuf = dzero
         if (l .le. sib%g(i)%g_nlu) then
            swibuf(1,1) = sib%g(i)%l(l)%pftswap_legacy%phase
         endif
         CHECK( nf90_put_var(ncid,swphaseid,swibuf, &
                    start=(/i,l/),count=(/1,1/)) )
         swibuf = 0
         if (l .le. sib%g(i)%g_nlu) then
            swibuf(1,1) = sib%g(i)%l(l)%pftswap_legacy%source_pftref
            swlivebuf(1,1,:,:) = sib%g(i)%l(l)%pftswap_legacy%live_lay(:,:)
            swdeadbuf(1,1,:,:) = sib%g(i)%l(l)%pftswap_legacy%dead_lay(:,:)
         endif
         CHECK( nf90_put_var(ncid,swsourceid,swibuf, &
                    start=(/i,l/),count=(/1,1/)) )
         CHECK( nf90_put_var(ncid,swliveid,swlivebuf, &
                    start=(/i,l,1,1/),count=(/1,1,npoolpft,nsoil/)) )
         CHECK( nf90_put_var(ncid,swdeadid,swdeadbuf, &
                    start=(/i,l,1,1/),count=(/1,1,npoollu,nsoil/)) )
      enddo
   enddo
   deallocate(swlivebuf,swdeadbuf)
endif

!...define variables...
do v=1, sibr_nvar

   if (dovar(v)) then
      select case (sibr_vref(v))

      case (1)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar(subcount))
          do i=1,subcount
             dvar(i) = sib%g(i)%gprogt%tmd
          enddo
          status = nf90_put_var(ncid,varid,dvar)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar)
          
      case (2)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar(subcount))
          do i=1,subcount
             dvar(i) = sib%g(i)%gprogt%seas_precip
          enddo
          status = nf90_put_var(ncid,varid,dvar)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar)

      case (3)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar(subcount))
          do i=1,subcount
             dvar(i) = sib%g(i)%gprogt%seas_tm
          enddo
          status = nf90_put_var(ncid,varid,dvar)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar)
          
      case (4)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar(subcount))
          do i=1,subcount
             dvar(i) = sib%g(i)%gprogt%clim_cupr
          enddo
          status = nf90_put_var(ncid,varid,dvar)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar)

      case (5)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar(subcount))
          do i=1,subcount
             dvar(i) = sib%g(i)%gprogt%clim_precip
          enddo
          status = nf90_put_var(ncid,varid,dvar)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar)
          
      case (6)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar(subcount))
          do i=1,subcount
             dvar(i) = sib%g(i)%gprogt%clim_tm
          enddo
          status = nf90_put_var(ncid,varid,dvar)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar)

     case (10)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%cast%eacas
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (11)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%cast%shcas
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (12)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%cast%tcas
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (13)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%cast%tkecas
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (14)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%cast%tc
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (15)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%cast%tcmin
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)
          
      case (20)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%co2t%co2m
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)
         
      case (21)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%co2t%pco2cas
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (22)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%co2t%pco2c
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (23)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%co2t%pco2i
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (24)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%co2t%pco2s
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)
          
    case (25)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar(subcount))
          do i=1,subcount
             dvar(i) = sib%g(i)%gprogt%pcosm
          enddo
          status = nf90_put_var(ncid,varid,dvar)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar)
          
      case (26)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%cost%coscasp
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

       case (27)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%cost%coscas
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)             

      case (28)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%co2t%co2cas
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (29)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%co2t%co2c
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (30)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
              
          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%co2t%co2i
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (31)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
              
          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%co2t%co2s
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)
          
      case (35)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%fluxt%ra
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)
          
      case (36)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%fluxt%rb
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)
          
      case (37)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%co2t%rst
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (40)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%co2t%assim
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)
          
      case (41)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%co2t%assimd
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (42)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%poollt%resp_nveg
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)
          
      case (43)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%poollt%resp_leaf
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (44)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%poollt%resp_root
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)
          
      case (45)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%poollt%resp_grow
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (46)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%poollt%resp_mntn
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)
          
      case (47)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%poollt%resp_auto
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (48)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%pooldt%resp_soil
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)
          
      case (49)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%pooldt%resp_het
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (50)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%poollt%resp_grz
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (51)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%poollt%resp_hrvst
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (60)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%hydrost%capacc_liq
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (61)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%hydrost%capacc_snow
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (62)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%hydrost%capacg
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (63)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%hydrovt%clim_pawfrw
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (64)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%hydrovt%clim_tawfrw
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (70)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_int, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(ivar2d(subcount,nlu))
          ivar2d(:,:) = izero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                ivar2d(i,l) = int(sib%g(i)%l(l)%sscolt%nsl)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,ivar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(ivar2d)

      case (71)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nsnowdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar3d(subcount,nlu,nsnow))
          dvar3d(:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar3d(i,l,:) = sib%g(i)%l(l)%sscolt%dz(-nsnow+1:0)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar3d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar3d)

      case (72)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,ntotdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar3d(subcount,nlu,ntot))
          dvar3d(:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar3d(i,l,:) = sib%g(i)%l(l)%sscolt%td(:)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar3d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar3d)

      case (73)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,ntotdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar3d(subcount,nlu,ntot))
          dvar3d(:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar3d(i,l,:) = sib%g(i)%l(l)%sscolt%www_liq(:)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar3d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar3d)

      case (74)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,ntotdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar3d(subcount,nlu,ntot))
          dvar3d(:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar3d(i,l,:) = sib%g(i)%l(l)%sscolt%www_ice(:)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar3d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar3d)

      case (80)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%co2t%clim_assim
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (81)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%vegt%clim_lai
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)


      case (82)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%phent%phenave_assim
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (83)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%phent%phenave_assimsm
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (84)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%phent%phenave_pr
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (85)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%phent%phenave_prsm
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (86)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%phent%phenave_prsdoy
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (87)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%phent%phenave_prcdoy
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)
          
      case (88)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%phent%phenave_tawftop
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (89)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%phent%phenave_tm
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (90)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%phent%phenc_climp
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (91)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%phent%phenc_laimax
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (92)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%phent%phenc_laimin
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (93)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%phent%phenave_env
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (94)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%phent%phenave_wa
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (95)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%phent%phenave_wacsm
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (96)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%phent%phen_pi
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)          

      case (97)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_int, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(ivar2d(subcount,nlu))
          ivar2d(:,:) = izero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                ivar2d(i,l) = sib%g(i)%l(l)%phent%phen_istage
             enddo
          enddo
          status = nf90_put_var(ncid,varid,ivar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(ivar2d)

      case (100)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_int, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(ivar2d(subcount,nlu))
          ivar2d(:,:) = izero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                ivar2d(i,l) = sib%g(i)%l(l)%phent%dapd
             enddo
          enddo
          status = nf90_put_var(ncid,varid,ivar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(ivar2d)

      case (101)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_int, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(ivar2d(subcount,nlu))
          ivar2d(:,:) = izero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                ivar2d(i,l) = sib%g(i)%l(l)%phent%dapdaf
             enddo
          enddo
          status = nf90_put_var(ncid,varid,ivar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(ivar2d)
          
      case (102)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%phent%gdd
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (103)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%phent%seed_pool
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (110)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nppdid,nsdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar4d(subcount,nlu,npoolpft,nsoil))
          dvar4d(:,:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                do p=1,npoolpft
                   dvar4d(i,l,p,:) = sib%g(i)%l(l)%poollt%poolpft_lay(p,:)
                enddo
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar4d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar4d)
          
      case (111)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nppdid,nsdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar4d(subcount,nlu,npoolpft,nsoil))
          dvar4d(:,:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                do p=1,npoolpft
                   dvar4d(i,l,p,:) = sib%g(i)%l(l)%poollt%poolpft_dgain(p,:)
                enddo
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar4d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar4d)

      case (112)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nppdid,nsdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar4d(subcount,nlu,npoolpft,nsoil))
          dvar4d(:,:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                do p=1,npoolpft
                   dvar4d(i,l,p,:) = sib%g(i)%l(l)%poollt%poolpft_dloss(p,:)
                enddo
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar4d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar4d)

      case (113)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nppdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar3d(subcount,nlu,npoolpft))
          dvar3d(:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar3d(i,l,:) = sib%g(i)%l(l)%poollt%gain_assim(:)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar3d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar3d)

      case (114)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nppdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar3d(subcount,nlu,npoolpft))
          dvar3d(:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar3d(i,l,:) = sib%g(i)%l(l)%poollt%gain_seed(:)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar3d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar3d)
          
      case (115)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nppdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar3d(subcount,nlu,npoolpft))
          dvar3d(:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar3d(i,l,:) = sib%g(i)%l(l)%poollt%loss_gresp(:)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar3d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar3d)

      case (116)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nppdid,nsdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar4d(subcount,nlu,npoolpft,nsoil))
          dvar4d(:,:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                do p=1,npoolpft
                   dvar4d(i,l,p,:) = sib%g(i)%l(l)%poollt%loss_mresp_lay(p,:)
                enddo
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar4d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar4d)

      case (117)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nppdid,nsdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar4d(subcount,nlu,npoolpft,nsoil))
          dvar4d(:,:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                do p=1,npoolpft
                   dvar4d(i,l,p,:) = sib%g(i)%l(l)%poollt%loss_trans_lay(p,:)
                enddo
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar4d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar4d)

      case (118)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nppdid,nsdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar4d(subcount,nlu,npoolpft,nsoil))
          dvar4d(:,:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                do p=1,npoolpft
                   dvar4d(i,l,p,:) = sib%g(i)%l(l)%poollt%loss_fire_lay(p,:)
                enddo
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar4d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar4d)

      case (119)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,ncpdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar3d(subcount,nlu,npoolcan))
          dvar3d(:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar3d(i,l,:) = sib%g(i)%l(l)%poollt%loss_grz(:)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar3d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar3d)
          
      case (120)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nppdid,nsdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar4d(subcount,nlu,npoolpft,nsoil))
          dvar4d(:,:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                do p=1,npoolpft
                   dvar4d(i,l,p,:) = sib%g(i)%l(l)%poollt%loss_hrvst_lay(p,:)
                enddo
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar4d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar4d)

      case (121)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nppdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar3d(subcount,nlu,npoolpft))
          dvar3d(:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar3d(i,l,:) = sib%g(i)%l(l)%equiblt%poolpft_totgain(:)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar3d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar3d)

      case (122)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nppdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar3d(subcount,nlu,npoolpft))
          dvar3d(:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar3d(i,l,:) = sib%g(i)%l(l)%equiblt%poolpft_totloss(:)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar3d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar3d)

      case (123)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nppdid,nsdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar4d(subcount,nlu,npoolpft,nsoil))
          dvar4d(:,:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                do p=1,npoolpft
                   dvar4d(i,l,p,:) = sib%g(i)%l(l)%poollt%poolpft_lay(p,:)
                enddo
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar4d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar4d)

      case (124)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,ncpdidc13/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar3d(subcount,nlu,npoolcanc13))
          dvar3d(:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar3d(i,l,:) = sib%g(i)%l(l)%poollt%loss_grzc13(:)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar3d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar3d)

      case (125)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nppdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar3d(subcount,nlu,npoolpft))
          dvar3d(:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar3d(i,l,:) = sib%g(i)%l(l)%equiblt%poolpft_equib(:)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar3d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar3d)

      case (126)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nppdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar3d(subcount,nlu,npoolpft))
          dvar3d(:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar3d(i,l,:) = sib%g(i)%l(l)%poollt%poolpftmin(:)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar3d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar3d)

      case (127)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nppdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar3d(subcount,nlu,npoolpft))
          dvar3d(:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar3d(i,l,:) = sib%g(i)%l(l)%equiblt%poolpft_init(:)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar3d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar3d)

      case (128)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nppdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar3d(subcount,nlu,npoolpft))
          dvar3d(:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar3d(i,l,:) = sib%g(i)%l(l)%equiblt%poolpft_end(:)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar3d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar3d)

      case (180)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nppdid,nsdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar4d(subcount,nlu,npoolpft,nsoil))
          dvar4d(:,:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                do p=1,npoolpft
                   dvar4d(i,l,p,:) = sib%g(i)%l(l)%poollt%loss_transf_lay(p,:)
                enddo
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar4d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar4d)

      case (130)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nlpdid,nsdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar4d(subcount,nlu,npoollu,nsoil))
          dvar4d(:,:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                do p=1,npoollu
                    dvar4d(i,l,p,:) = sib%g(i)%l(l)%pooldt%poollu_lay(p,:)
                 enddo
              enddo   
          enddo
          status = nf90_put_var(ncid,varid,dvar4d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar4d)

      case (131)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nlpdid,nsdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar4d(subcount,nlu,npoollu,nsoil))
          dvar4d(:,:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                do p=1,npoollu
                   dvar4d(i,l,p,:) = sib%g(i)%l(l)%pooldt%poollu_dgain(p,:)
                enddo
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar4d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar4d)

      case (132)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nlpdid,nsdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar4d(subcount,nlu,npoollu,nsoil))
          dvar4d(:,:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                do p=1,npoollu
                   dvar4d(i,l,p,:) = sib%g(i)%l(l)%pooldt%poollu_dloss(p,:)
                enddo
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar4d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar4d)

      case (133)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nlpdid,nsdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar4d(subcount,nlu,npoollu,nsoil))
          dvar4d(:,:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                do p=1,npoollu
                   dvar4d(i,l,p,:) = sib%g(i)%l(l)%pooldt%gain_grz_lay(p,:)
                enddo
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar4d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar4d)

      case (134)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nlpdid,nsdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar4d(subcount,nlu,npoollu,nsoil))
          dvar4d(:,:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                do p=1,npoollu
                   dvar4d(i,l,p,:) = sib%g(i)%l(l)%pooldt%gain_hrvst_lay(p,:)
                enddo
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar4d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar4d)

      case (135)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nlpdid,nsdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar4d(subcount,nlu,npoollu,nsoil))
          dvar4d(:,:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                do p=1,npoollu
                   dvar4d(i,l,p,:) = sib%g(i)%l(l)%pooldt%gain_transl_lay(p,:)
                enddo
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar4d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar4d)

      case (136)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nlpdid,nsdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar4d(subcount,nlu,npoollu,nsoil))
          dvar4d(:,:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                do p=1,npoollu
                   dvar4d(i,l,p,:) = sib%g(i)%l(l)%pooldt%gain_transd_lay(p,:)
                enddo
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar4d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar4d)

      case (137)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nlpdid,nsdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar4d(subcount,nlu,npoollu,nsoil))
          dvar4d(:,:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                do p=1,npoollu
                   dvar4d(i,l,p,:) = sib%g(i)%l(l)%pooldt%loss_fire_lay(p,:)
                enddo
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar4d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar4d)

      case (138)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nlpdid,nsdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar4d(subcount,nlu,npoollu,nsoil))
          dvar4d(:,:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                do p=1,npoollu
                   dvar4d(i,l,p,:) = sib%g(i)%l(l)%pooldt%loss_resp_lay(p,:)
                enddo
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar4d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar4d)

      case (139)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nlpdid,nsdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar4d(subcount,nlu,npoollu,nsoil))
          dvar4d(:,:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                do p=1,npoollu
                   dvar4d(i,l,p,:) = sib%g(i)%l(l)%pooldt%loss_trans_lay(p,:)
                enddo
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar4d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar4d)
          
      case (140)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nlpdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar3d(subcount,nlu,npoollu))
          dvar3d(:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar3d(i,l,:) = sib%g(i)%l(l)%equibdt%poollu_totgain(:)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar3d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar3d)

      case (141)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nlpdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar3d(subcount,nlu,npoollu))
          dvar3d(:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar3d(i,l,:) = sib%g(i)%l(l)%equibdt%poollu_totloss(:)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar3d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar3d)

      case (142)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nlpdid,nsdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar4d(subcount,nlu,npoollu,nsoil))
          dvar4d(:,:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                do p=1,npoollu
                   dvar4d(i,l,p,:) = sib%g(i)%l(l)%pooldt%poollu_lay(p,:)
                enddo
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar4d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar4d)

      case (143)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nlpdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar3d(subcount,nlu,npoollu))
          dvar3d(:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar3d(i,l,:) = sib%g(i)%l(l)%equibdt%poollu_equib(:)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar3d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar3d)

      case (144)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nlpdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar3d(subcount,nlu,npoollu))
          dvar3d(:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar3d(i,l,:) = sib%g(i)%l(l)%equibdt%poollu_init(:)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar3d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar3d)

      case (145)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nlpdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar3d(subcount,nlu,npoollu))
          dvar3d(:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar3d(i,l,:) = sib%g(i)%l(l)%equibdt%poollu_end(:)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar3d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar3d)

      case (181)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nlpdid,nsdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar4d(subcount,nlu,npoollu,nsoil))
          dvar4d(:,:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                do p=1,npoollu
                   dvar4d(i,l,p,:) = sib%g(i)%l(l)%pooldt%gain_transf_lay(p,:)
                enddo
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar4d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar4d)


      case (147)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%fract%kiecps
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (148)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%fract%d13cca
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (149)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%fract%rcassim
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (150)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%fract%c13assim
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (151)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%fract%c13assimd
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (152)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%poollt%resp_nvegc13
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (153)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%poollt%resp_leafc13
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (154)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%poollt%resp_rootc13
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (155)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%poollt%resp_growc13
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (156)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%poollt%resp_mntnc13
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (157)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%poollt%resp_autoc13
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (158)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%pooldt%resp_soilc13
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (159)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%pooldt%resp_hetc13
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (160)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%poollt%resp_grzc13
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (161) 
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nppdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar3d(subcount,nlu,npoolpft))
          dvar3d(:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar3d(i,l,:) = sib%g(i)%l(l)%poollt%rcpoolpft(:)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar3d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar3d)

      case (162) 
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nppdid,nsdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar4d(subcount,nlu,npoolpft,nsoil))
          dvar4d(:,:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                do p=1,npoolpft
                   dvar4d(i,l,p,:) = sib%g(i)%l(l)%poollt%rcpoolpft_lay(p,:)
                enddo
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar4d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar4d)

      case (163)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nlpdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar3d(subcount,nlu,npoollu))
          dvar3d(:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar3d(i,l,:) = sib%g(i)%l(l)%pooldt%rcpoollu(:)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar3d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar3d)

      case (164) 
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nlpdid,nsdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar4d(subcount,nlu,npoollu,nsoil))
          dvar4d(:,:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                do p=1,npoollu
                   dvar4d(i,l,p,:) = sib%g(i)%l(l)%pooldt%rcpoollu_lay(p,:)
                enddo
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar4d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar4d)

      case (165)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%fract%rcpoolfirec13
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (166)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%fract%rcassimfac
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (167)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%fract%poolemistotC
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (168)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount 
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%fract%poolemisc13
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (170)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%fract%c12ca
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (171)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount 
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%fract%c13ca
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (172)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount 
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%fract%c12assim
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (173)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%fract%c13resptot
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (174)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount 
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%fract%c12resptot
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (175)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%poollt%resp_hrvstc13
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)





      case (182)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nlpdid,nsdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar4d(subcount,nlu,npoollu,nsoil))
          dvar4d(:,:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                do p=1,npoollu
                   dvar4d(i,l,p,:) = sib%g(i)%l(l)%pooldt%loss_transf_lay(p,:)
                enddo
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar4d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar4d)

      case (188)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nppdid,nsdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar4d(subcount,nlu,npoolpft,nsoil))
          dvar4d(:,:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                do p=1,npoolpft
                   dvar4d(i,l,p,:) = sib%g(i)%l(l)%poollt%loss_raddecay_lay(p,:)
                enddo
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar4d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar4d)

      case (189)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nlpdid,nsdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar4d(subcount,nlu,npoollu,nsoil))
          dvar4d(:,:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                do p=1,npoollu
                   dvar4d(i,l,p,:) = sib%g(i)%l(l)%pooldt%loss_raddecay_lay(p,:)
                enddo
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar4d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar4d)

      case (190)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,ncpdidc14/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar3d(subcount,nlu,npoolcanc14))
          dvar3d(:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar3d(i,l,:) = sib%g(i)%l(l)%poollt%loss_grzc14(:)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar3d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar3d)

      case (191)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%fract%rcpoolfirec14
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (192)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%fract%rcassimfacc14
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (193)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%fract%poolemisc14
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (194)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%fract%c14alpha
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (195)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%fract%d14cca
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (196)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%fract%rcassimc14
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (197)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%fract%c14assim
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (198)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%fract%c14assimd
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (199)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%poollt%resp_nvegc14
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (200)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%poollt%resp_leafc14
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (201)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%poollt%resp_rootc14
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (202)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%poollt%resp_growc14
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (203)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%poollt%resp_mntnc14
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (204)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%poollt%resp_autoc14
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (205)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%pooldt%resp_soilc14
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (206)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%pooldt%resp_hetc14
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (207)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%poollt%resp_grzc14
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (208)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%poollt%resp_hrvstc14
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (209)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%fract%c14ca
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (210)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%fract%c14resptot
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (211)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%phent%phenave_assimpot
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (213)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%phent%phenave_prpot
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (222)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%phent%phens_dayl
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (223)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%phent%phens_grw
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (224)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%phent%phenave_wac
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (225)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%phent%phens_wx
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (212)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_int, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(ivar2d(subcount,nlu))
          ivar2d(:,:) = izero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                ivar2d(i,l) = merge(ione,izero,sib%g(i)%l(l)%phent%phenflag_assimlow)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,ivar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(ivar2d)

      case (214)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_int, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(ivar2d(subcount,nlu))
          ivar2d(:,:) = izero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                ivar2d(i,l) = merge(ione,izero,sib%g(i)%l(l)%phent%phenflag_precip)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,ivar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(ivar2d)

      case (215)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_int, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(ivar2d(subcount,nlu))
          ivar2d(:,:) = izero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                ivar2d(i,l) = merge(ione,izero,sib%g(i)%l(l)%phent%phenflag_moist)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,ivar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(ivar2d)

      case (216)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_int, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(ivar2d(subcount,nlu))
          ivar2d(:,:) = izero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                ivar2d(i,l) = merge(ione,izero,sib%g(i)%l(l)%phent%phenflag_temp)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,ivar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(ivar2d)

      case (217)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_int, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(ivar2d(subcount,nlu))
          ivar2d(:,:) = izero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                ivar2d(i,l) = merge(ione,izero,sib%g(i)%l(l)%phent%phenflag_daylen)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,ivar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(ivar2d)

      case (218)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_int, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(ivar2d(subcount,nlu))
          ivar2d(:,:) = izero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                ivar2d(i,l) = merge(ione,izero,sib%g(i)%l(l)%phent%phenflag_gsspass)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,ivar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(ivar2d)

      case (219)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_int, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(ivar2d(subcount,nlu))
          ivar2d(:,:) = izero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                ivar2d(i,l) = sib%g(i)%l(l)%phent%nd_dormant
             enddo
          enddo
          status = nf90_put_var(ncid,varid,ivar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(ivar2d)

      case (220)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_int, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(ivar2d(subcount,nlu))
          ivar2d(:,:) = izero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                ivar2d(i,l) = sib%g(i)%l(l)%phent%nd_gs
             enddo
          enddo
          status = nf90_put_var(ncid,varid,ivar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(ivar2d)

      case (226)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_int, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(ivar2d(subcount,nlu))
          ivar2d(:,:) = izero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                ivar2d(i,l) = sib%g(i)%l(l)%phent%ipd
             enddo
          enddo
          status = nf90_put_var(ncid,varid,ivar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(ivar2d)

      case (221)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_int, &
                        (/nsibdid,nludid,nstgdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(ivar3d(subcount,nlu,npstgmax))
          ivar3d(:,:,:) = izero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                ivar3d(i,l,1:size(sib%g(i)%l(l)%phent%nd_stg)) = &
                     sib%g(i)%l(l)%phent%nd_stg(:)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,ivar3d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(ivar3d)

      case (230)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nppdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(dvar3d(subcount,nlu,npoolpft))
          dvar3d(:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar3d(i,l,:) = sib%g(i)%l(l)%poollt%poolpftp(:)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar3d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar3d)

      case (231)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nlpdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(dvar3d(subcount,nlu,npoollu))
          dvar3d(:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar3d(i,l,:) = sib%g(i)%l(l)%pooldt%poollup(:)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar3d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar3d)

      case (232)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nppdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(dvar3d(subcount,nlu,npoolpft))
          dvar3d(:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar3d(i,l,:) = sib%g(i)%l(l)%poollt%alloc_phen(:)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar3d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar3d)

      case (233)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_int, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(ivar2d(subcount,nlu))
          ivar2d(:,:) = izero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                ivar2d(i,l) = merge(ione,izero,sib%g(i)%l(l)%poollt%aadj_moist)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,ivar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(ivar2d)

      case (234)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_int, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(ivar2d(subcount,nlu))
          ivar2d(:,:) = izero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                ivar2d(i,l) = merge(ione,izero,sib%g(i)%l(l)%poollt%aadj_temp)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,ivar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(ivar2d)

      case (235)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%poollt%tfl_pstage
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (236)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%poollt%nd_fire
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (237)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%poollt%nd_grz
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (238)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%pooldt%mhrt_soil_hot
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (239)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nsdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(dvar3d(subcount,nlu,nsoil))
          dvar3d(:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar3d(i,l,:) = sib%g(i)%l(l)%pooldt%resp_soilnr_lay(:)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar3d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar3d)

      case (240)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%hydrovt%pawfrw
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (241)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%hydrovt%pawftop
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (242)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%hydrovt%pawfzw
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (243)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%hydrovt%tawfrw
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (244)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%hydrovt%tawftop
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (245)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%hydrovt%tawfzw
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (246)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%hydrost%snow_gmass
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (247)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%hydrost%snow_gdepth
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (248)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%hydrost%snow_cvfc
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (249)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%hydrost%wetfracc
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (250)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%hydrost%wetfracg
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (255)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%vegt%lai
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (256)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%vegt%lait
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (257)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%vegt%fpar
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (258)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%vegt%vcover
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (259)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%vegt%green
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (260)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%vegt%z0
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (261)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%vegt%z0d
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (262)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%vegt%zp_dispd
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (263)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%vegt%zpd_adj
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (264)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%vegt%zztemp
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (265)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%vegt%zzwind
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (266)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%vegt%cc1
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (267)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%vegt%cc2
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (268)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%vegt%gmudmu
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (269)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%vegt%vmax
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (275)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%co2t%rstfac(1)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (276)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%co2t%rstfac(2)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (277)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%co2t%rstfac(3)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (278)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%co2t%rstfac(4)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (279)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nsdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(dvar3d(subcount,nlu,nsoil))
          dvar3d(:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar3d(i,l,:) = sib%g(i)%l(l)%cost%cos_s(:)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar3d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar3d)

      case (280:337)
          ! V3 exact-continuation external forcing brackets.
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(dvar(subcount))
          dvar(:) = dzero
          select case (sibr_vref(v))
          case (280)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%tm1
             enddo
          case (281)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%tm2
             enddo
          case (282)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%sh1
             enddo
          case (283)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%sh2
             enddo
          case (284)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%ps1
             enddo
          case (285)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%ps2
             enddo
          case (286)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%spdm1
             enddo
          case (287)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%spdm2
             enddo
          case (288)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%lspr1
             enddo
          case (289)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%lspr2
             enddo
          case (290)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%cupr1
             enddo
          case (291)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%cupr2
             enddo
          case (292)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%dlwbot1
             enddo
          case (293)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%dlwbot2
             enddo
          case (294)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%sw_dwn1
             enddo
          case (295)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%sw_dwn2
             enddo
          case (296)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%firec1
             enddo
          case (297)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%firec2
             enddo
          case (298)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%fireco21
             enddo
          case (299)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%fireco22
             enddo
          case (300)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%ba_dbg1
             enddo
          case (301)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%ba_dbg2
             enddo
          case (302)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%ba_en21
             enddo
          case (303)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%ba_en22
             enddo
          case (304)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%ba_en31
             enddo
          case (305)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%ba_en32
             enddo
          case (306)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%ba_dnf1
             enddo
          case (307)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%ba_dnf2
             enddo
          case (308)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%ba_eb11
             enddo
          case (309)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%ba_eb12
             enddo
          case (310)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%ba_eb21
             enddo
          case (311)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%ba_eb22
             enddo
          case (312)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%ba_db11
             enddo
          case (313)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%ba_db12
             enddo
          case (314)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%ba_db21
             enddo
          case (315)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%ba_db22
             enddo
          case (316)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%ba_db31
             enddo
          case (317)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%ba_db32
             enddo
          case (318)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%ba_shb1
             enddo
          case (319)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%ba_shb2
             enddo
          case (320)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%ba_sha1
             enddo
          case (321)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%ba_sha2
             enddo
          case (322)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%ba_c3a1
             enddo
          case (323)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%ba_c3a2
             enddo
          case (324)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%ba_c3g1
             enddo
          case (325)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%ba_c3g2
             enddo
          case (326)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%ba_c4g1
             enddo
          case (327)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%ba_c4g2
             enddo
          case (328)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%ba_c3c1
             enddo
          case (329)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%ba_c3c2
             enddo
          case (330)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%ba_c4c1
             enddo
          case (331)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%ba_c4c2
             enddo
          case (332)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%ba_mze1
             enddo
          case (333)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%ba_mze2
             enddo
          case (334)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%ba_soy1
             enddo
          case (335)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%ba_soy2
             enddo
          case (336)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%ba_wwt1
             enddo
          case (337)
             do i=1,subcount
                dvar(i) = sib%g(i)%gprogt%ba_wwt2
             enddo
          end select
          status = nf90_put_var(ncid,varid,dvar)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar)

      case (338:339)
          ! V4 exact-continuation previous-timestep diagnostic inputs.
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                select case (sibr_vref(v))
                case (338)
                   dvar2d(i,l) = sib%g(i)%l(l)%cast%vpd
                case (339)
                   dvar2d(i,l) = sib%g(i)%l(l)%fluxt%ect
                end select
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      !case default
      !   print*,'Invalid Restart Variable: ', sibr_vname(v)
      !   stop

      end select
   endif !dovar
enddo


!---------------------------------------------------------------
! Preserve unfinished monthly QP accumulation in the exact transient
! checkpoint, independently of sibr_outopts (which describes physical
! state, not the output schema).  Values here are RAW timestep sums,
! not the normalized monthly averages.  Requib retains old semantics.
!---------------------------------------------------------------
if ((.not. equib) .and. qp_dtsib == -1 .and. &
    (qp_savegf .or. qp_saveluf)) then
   CHECK(nf90_redef(ncid))
   CHECK(nf90_put_att(ncid,nf90_global,'qp_continuation_schema',1))
   CHECK(nf90_put_att(ncid,nf90_global,'qp_continuation_dtsib',qp_dtsib))
   CHECK(nf90_put_att(ncid,nf90_global,'qp_continuation_countt',qp_countt))
   CHECK(nf90_put_att(ncid,nf90_global,'qp_continuation_nvarg',qp_nvarg))
   CHECK(nf90_put_att(ncid,nf90_global,'qp_continuation_nvarlu',qp_nvarlu))
   CHECK(nf90_def_dim(ncid,'qp_continuation_sif2',2,qp_sifdim))
   CHECK(nf90_def_var(ncid,'qp_continuation_sif_counts',nf90_int, &
         (/nsibdid,qp_sifdim/),qp_sifvar))
   if (qp_savegf) then
      CHECK(nf90_def_dim(ncid,'qp_continuation_ng',qp_nvarg,qp_ngdim))
      CHECK(nf90_def_var(ncid,'qp_continuation_g',nf90_double, &
            (/nsibdid,qp_ngdim/),qp_gvar))
      CHECK(nf90_def_var(ncid,'qp_continuation_g_ids',nf90_int, &
            (/qp_ngdim/),qp_gidvar))
   endif
   if (qp_saveluf) then
      CHECK(nf90_def_dim(ncid,'qp_continuation_nluvars',qp_nvarlu,qp_nludim))
      CHECK(nf90_def_var(ncid,'qp_continuation_lu',nf90_double, &
            (/nsibdid,nludid,qp_nludim/),qp_luvar))
      CHECK(nf90_def_var(ncid,'qp_continuation_lu_ids',nf90_int, &
            (/qp_nludim/),qp_luidvar))
   endif
   CHECK(nf90_enddef(ncid))
   CHECK(nf90_put_var(ncid,qp_sifvar,qp_sif_satcount(:,1,:)))
   if (qp_savegf) then
      CHECK(nf90_put_var(ncid,qp_gvar,qp_g(:,:,1)))
      CHECK(nf90_put_var(ncid,qp_gidvar,qp_vrefg))
   endif
   if (qp_saveluf) then
      CHECK(nf90_put_var(ncid,qp_luvar,qp_lu(:,:,:,1)))
      CHECK(nf90_put_var(ncid,qp_luidvar,qp_vreflu))
   endif
endif

! Retain the last successfully computed allocation fractions across exact
! continuation, including dormant intervals when pool_assim returns early.
! This is independent of the experiment-specific sib_outopts selection and
! is NOT included in equilibrium/requib initialization.
if (.not. equib) then
   allocate(alloc_continue(subcount,nlu,npoolpft))
   allocate(alloc_moist_continue(subcount,nlu,npoolpft))
   allocate(alloc_temp_continue(subcount,nlu,npoolpft))
   alloc_continue = dzero
   alloc_moist_continue = dzero
   alloc_temp_continue = dzero
   do i=1,subcount
      do l=1,sib%g(i)%g_nlu
         alloc_continue(i,l,:) = sib%g(i)%l(l)%poollt%alloc(:)
         alloc_moist_continue(i,l,:) = sib%g(i)%l(l)%poollt%alloc_moist(:)
         alloc_temp_continue(i,l,:) = sib%g(i)%l(l)%poollt%alloc_temp(:)
      enddo
   enddo
   CHECK(nf90_redef(ncid))
   CHECK(nf90_def_var(ncid,'restart_alloc_current',nf90_double, &
         (/nsibdid,nludid,nppdid/),alloc_continue_var))
   CHECK(nf90_def_var(ncid,'restart_alloc_moist',nf90_double, &
         (/nsibdid,nludid,nppdid/),alloc_moist_var))
   CHECK(nf90_def_var(ncid,'restart_alloc_temp',nf90_double, &
         (/nsibdid,nludid,nppdid/),alloc_temp_var))
   CHECK(nf90_enddef(ncid))
   CHECK(nf90_put_var(ncid,alloc_continue_var,alloc_continue))
   CHECK(nf90_put_var(ncid,alloc_moist_var,alloc_moist_continue))
   CHECK(nf90_put_var(ncid,alloc_temp_var,alloc_temp_continue))
   deallocate(alloc_continue,alloc_moist_continue,alloc_temp_continue)
endif

! Held, non-prognostic diagnostics needed for bitwise continuation.
! These fields retain the previous timestep's value when their respective
! calculations are inactive.  We preserve native values rather than reset
! diagnostics or change their scientific update conditions.
!
! Scalar channel order (held_diag_count = 10):
!  1 d13cassim; 2 d13cassim_nog; 3 d13cresp_growc13;
!  4..8 kiecps_k1..k5; 9 kiecps_nog; 10 sif_x.
! Assimilation gain/loss scratch arrays are stored in full, retaining both
! C and isotope pool/layer bookkeeping.  Requib files are unaffected.
if (.not. equib) then
   allocate(held_scalars(subcount,nlu,nheld_diag))
   allocate(held_gain(subcount,nlu,npoolpft,nsoil))
   allocate(held_loss(subcount,nlu,npoolpft,nsoil))
   held_scalars=dzero
   held_gain=dzero
   held_loss=dzero
   do i=1,subcount
      do l=1,sib%g(i)%g_nlu
         held_scalars(i,l,1)=sib%g(i)%l(l)%fract%d13cassim
         held_scalars(i,l,2)=sib%g(i)%l(l)%fract%d13cassim_nog
         held_scalars(i,l,3)=sib%g(i)%l(l)%fract%d13cresp_growc13
         held_scalars(i,l,4)=sib%g(i)%l(l)%fract%kiecps_k1
         held_scalars(i,l,5)=sib%g(i)%l(l)%fract%kiecps_k2
         held_scalars(i,l,6)=sib%g(i)%l(l)%fract%kiecps_k3
         held_scalars(i,l,7)=sib%g(i)%l(l)%fract%kiecps_k4
         held_scalars(i,l,8)=sib%g(i)%l(l)%fract%kiecps_k5
         held_scalars(i,l,9)=sib%g(i)%l(l)%fract%kiecps_nog
         held_scalars(i,l,10)=sib%g(i)%l(l)%sift%sif_x
         held_gain(i,l,:,:)=sib%g(i)%l(l)%poollt%assim_dgain(:,:)
         held_loss(i,l,:,:)=sib%g(i)%l(l)%poollt%assim_dloss(:,:)
      enddo
   enddo
   CHECK(nf90_redef(ncid))
   CHECK(nf90_put_att(ncid,nf90_global,'held_diag_schema',1))
   CHECK(nf90_def_dim(ncid,'held_diag_count',nheld_diag,held_dim))
   CHECK(nf90_def_var(ncid,'held_diag_scalars',nf90_double, &
         (/nsibdid,nludid,held_dim/),held_var))
   CHECK(nf90_def_var(ncid,'held_assim_dgain',nf90_double, &
         (/nsibdid,nludid,nppdid,nsdid/),held_gain_var))
   CHECK(nf90_def_var(ncid,'held_assim_dloss',nf90_double, &
         (/nsibdid,nludid,nppdid,nsdid/),held_loss_var))
   CHECK(nf90_enddef(ncid))
   CHECK(nf90_put_var(ncid,held_var,held_scalars))
   CHECK(nf90_put_var(ncid,held_gain_var,held_gain))
   CHECK(nf90_put_var(ncid,held_loss_var,held_loss))
   deallocate(held_scalars,held_gain,held_loss)
endif

!...close the file
CHECK( nf90_close( ncid ) )


end subroutine restart_write
