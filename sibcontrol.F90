!----------------------------------------------------------------------
subroutine sib_control ( gnum, gref,  &
       lonsib, latsib, daynew, daylmax, &
       doy, sibg )

    use kinds
    use module_param, only: &
       phencon, physcon, poolcon
    use module_pftinfo
    use module_poolinfo, only: pool_indx_lay
    use module_oparams, only: canlai_min
    use module_sib, only: &
        gridcell_type, sib
    use module_sibconst, only: nsoil, npoolpft, npoollu, ntpool, subcount
    use module_time, only: year, startyear, starttime, &
         sec_day, sec_tot, new_day, dayfrac_lst, new_day_lst
    use module_io, only: pftswap_enable, pftswap_open_year, pftswap_open_doy, &
         minveg_diag, ic_restart, rank, out_path, driver_recnum, driver_year, &
         driver_month, driver_day, driver_hour, driver_seccur, &
         driver_secnext, driver_readf, driver_switchf, driver_updatef, &
         hr_saveluf, hr_countt, hr_countsave, hr_countwrite, &
         hr_nvarlu, hr_vreflu, hr_lu

    implicit none

    !---INPUT VARIABLES
    integer(i4), intent(in) :: gnum, gref, doy
    real(r4), intent(in) :: lonsib, latsib, daylmax
    logical, intent(in) :: daynew
    type (gridcell_type), intent(inout) :: sibg

    !---LOCAL VARIABLES
    integer(i4) :: l,pref,pnum, hnum
    integer(i4) :: vegtype
    logical :: isbare, iscrop, isgrass, isveg
    logical :: hold_minveg_first_day
    real(r8) :: legacy_rootf(nsoil)
    integer(i4), save :: trace_unit = -1
    logical, save :: trace_opened = .false.
    
!-------------------------------------------------------------------------

! Loop over land units
do l=1,sibg%g_nlu

   if ((.not. sibg%l(l)%equibdt%lupft_spunup)) then

       pref = sibg%l(l)%ipft
       pnum = pft_num(pref)
       vegtype = pft_type(pnum)

       isbare = (vegtype .eq. type_bare)
       iscrop = (vegtype .eq. type_crop)
       isgrass = (vegtype .eq. type_grass)
       isveg = (.not. isbare)
       if (pftswap_enable .and. (subcount .ne. 1 .or. sibg%g_nlu .ne. 1)) &
            error stop 'PFTSWAP Step1C supports only one single-PFT site'

       call trace_restart('ENTRY_AFTER_DRIVER', l)

       if (isveg) then
          ! Match the C2/C3 phen_update interface.  Never hold stage 1
          ! for true restarts or runs with the minveg diagnostic disabled.
          hold_minveg_first_day = minveg_diag .and. (.not. ic_restart) &
               .and. (year .eq. startyear) .and. (doy .eq. starttime)
          !...Update phenology
          !...(potentials every time-step, stage once daily)
          call phen_update( &
               gnum, l, pnum, sibg%l(l)%ipft, &
               phencon(pnum), poolcon(pnum), &
               daynew, daylmax, sibg%gdiagt%daylen, &
               sibg%gdiagt%daylendt, &
               sibg%gprogt%clim_cupr, sibg%gprogt%clim_precip, &
               sibg%gprogt%cupr, sibg%gprogt%lspr, sibg%gdiagt%tmdf, &
               sibg%gprogt%tm, sibg%l(l)%co2t%assim, &
               sibg%l(l)%co2t%rstfac(2), sibg%l(l)%co2t%rstfac(4), &
               sibg%l(l)%vegt%lai, sibg%l(l)%hydrovt, sibg%l(l)%phent, &
               sibg%l(l)%pooldt, sibg%l(l)%poollt, sibg%l(l)%vegt%vmax, &
               sibg%l(l)%fract%c13assim, physcon(pnum), sibg%l(l)%fract, &
               hold_minveg_first_day, &
               minveg_diag .and. (.not. ic_restart) .and. &
               (year .eq. startyear) .and. (doy .ge. starttime) .and. &
               (doy .le. starttime + 7), doy)

          call trace_restart('AFTER_PHEN_UPDATE', l)
       endif ! Existing vegetation phenology (open ground has none).

       if (daynew .and. (isveg .or. &
            sibg%l(l)%pftswap_legacy%phase .eq. 1)) then

             !...Update pools (once daily)
             hnum = pnum
             if (sibg%l(l)%pftswap_legacy%phase .eq. 1) &
                hnum = pft_num(sibg%l(l)%pftswap_legacy%source_pftref)
             call pool_update(gref, lonsib, latsib, pref, &
                  sibg%l(l)%co2t%assim, sibg%l(l)%fract%c13assim, &
                  sibg%l(l)%fract%c14assim, &
                  sibg%l(l)%vegt%lai, sibg%l(l)%vegt%fpar, &
                  sibg%l(l)%equibdt, sibg%l(l)%pooldt, &
                  sibg%l(l)%equiblt, sibg%l(l)%poollt, &
                  poolcon(hnum)%graze_trans, &
                  sibg%l(l)%co2t, &
                  sibg%l(l)%fract)

             call trace_restart('AFTER_POOL_UPDATE', l)

             if (pftswap_enable .and. isveg .and. &
                  sibg%l(l)%pftswap_legacy%phase .eq. 0 .and. &
                  year .eq. pftswap_open_year .and. doy .eq. pftswap_open_doy) then
                if (pref .ne. 11) error stop 'PFTSWAP opening source is not shrub (11)'
                if (abs(sibg%gprogt%fireco2) .gt. 1.d-14) &
                    error stop 'PFTSWAP refuses opening while fire forcing nonzero'
                call pftswap_open_transaction(gnum, l, sibg%l(l))
                pref = sibg%l(l)%ipft
                pnum = pft_num(pref)
                isbare = .true.
                isveg = .false.
                isgrass = .false.
                iscrop = .false.
             endif

             !...Update vegetation (once daily)
             call veg_update(doy, &
                  gref, lonsib, latsib, &
                  pnum, pref, iscrop, isgrass, &
                  physcon(pnum), sibg%l(l)%pooldt%poollu, &
                  sibg%l(l)%hydrost%snow_cvfc, &
                  sibg%l(l)%sscolt%node_z(1:nsoil), &
                  sibg%l(l)%poollt, sibg%l(l)%vegt)

             if (sibg%l(l)%pftswap_legacy%phase .eq. 1) then
                if (sibg%l(l)%vegt%lai .ge. canlai_min) &
                    error stop 'PFTSWAP bare LAI above photosynthetic canopy threshold'
             endif
             call trace_restart('AFTER_VEG_UPDATE', l)
       endif

       !...Call SiB Main
       call sib_main(gref, lonsib, latsib, &
            pref, pnum, physcon(pnum), &
            sibg%gprogt, sibg%gdiagt, sibg%l(l))

       call trace_restart('AFTER_SIB_MAIN', l)

       !...For phase 1, phonveg() has set assimilation to zero.
       ! Retain these native routines to RESET autotrophic flux diagnostics
       ! and prevent stale previous-shrub respiration from reaching balance.
       !...Calculate Pool Assimilation Gains
       call pool_assim( &
            gref, lonsib, latsib, pref, &
            physcon(pnum)%hhti, physcon(pnum)%hlti, &
            physcon(pnum)%shti, physcon(pnum)%slti, &
            poolcon(pnum)%gr_frac, sibg%l(l)%co2t%assim, &
            sibg%l(l)%fract%c13assim, sibg%l(l)%fract%c14assim,&
            sibg%l(l)%co2t%rstfac(2), &
            sibg%gprogt%tm, sibg%l(l)%poollt, sibg%l(l)%co2t, &
            sibg%l(l)%fract)

       call trace_restart('AFTER_POOL_ASSIM', l)

       !...Autotrophic Respiration
       call pool_auto_resp(poolcon(pnum), &
            sibg%l(l)%co2t%clim_assim, sibg%l(l)%vegt%clim_lai, &
            sibg%l(l)%co2t%assimd, sibg%l(l)%fract%c13assimd, & 
            sibg%l(l)%fract%c14assimd, sibg%l(l)%vegt%lai, &
            sibg%l(l)%cast%tc, sibg%l(l)%sscolt%td(1:nsoil), &
            sibg%l(l)%vegt%rootf, sibg%l(l)%poollt, &
            sibg%l(l)%pooldt%resp_soil, sibg%l(l)%pooldt%resp_soil_lay, &
            sibg%l(l)%pooldt%resp_soilc13, sibg%l(l)%pooldt%resp_soilc13_lay, &
            sibg%l(l)%pooldt%resp_soilc14, sibg%l(l)%pooldt%resp_soilc14_lay, &
            sibg%l(l)%fract)

       call trace_restart('AFTER_AUTO_RESP', l)

       !...Autotrophic Transfer
       call pool_auto_tran(poolcon(pnum), &
            sibg%gdiagt%daylen, sibg%gdiagt%daylendt, &
            daylmax, sibg%l(l)%cast%tc, sibg%l(l)%hydrovt%pawfrw, &
            sibg%l(l)%vegt%rootf, sibg%l(l)%poollt, &
            sibg%l(l)%pooldt%gain_transl_lay, &
            sibg%l(l)%pooldt%poollu_dgain, sibg%l(l)%fract) 

       call trace_restart('AFTER_AUTO_TRAN', l)

       !...Grazing
       IF (physcon(pnum)%pftgraze .and. &
           sibg%l(l)%pftswap_legacy%phase .ne. 1) THEN
          !call pool_graze(poolcon(pnum)%poolpft_min, &
          !     poolcon(pnum)%graze_trans, &
          !     sibg%l(l)%vegt%clim_lai, sibg%l(l)%vegt%lai, &
          !     sibg%l(l)%poollt, sibg%l(l)%pooldt)
          call pool_graze(poolcon(pnum), &
               sibg%l(l)%vegt%clim_lai, sibg%l(l)%vegt%lai, &
               sibg%l(l)%poollt, sibg%l(l)%pooldt, sibg%l(l)%fract)
       ENDIF

       ! Existing shrub-origin dead pools remain in native LU soil storage,
       ! but retain SHRUB kinetics during the open phase. No double-counting
       ! in the separate legacy arrays (zero until distinct cohorts exist).
       hnum = pnum
       if (sibg%l(l)%pftswap_legacy%phase .eq. 1) then
          hnum = pft_num(sibg%l(l)%pftswap_legacy%source_pftref)
          if (hnum .le. 0) error stop 'PFTSWAP invalid source PFT tag'
          sibg%l(l)%pooldt%resp_soil = dzero
          sibg%l(l)%pooldt%resp_soilc13 = dzero
          sibg%l(l)%pooldt%resp_soilc14 = dzero
          sibg%l(l)%pooldt%resp_soil_lay = dzero
          sibg%l(l)%pooldt%resp_soilc13_lay = dzero
          sibg%l(l)%pooldt%resp_soilc14_lay = dzero
       endif
       !...Heterotrophic Respiration/Transfer
       ! In open phase the native dead pools still contain shrub-origin
       ! carbon.  Do not supply the new bare-PFT rooting profile to the
       ! source-shrub decay and surface-to-soil transfer operator.
       legacy_rootf = sibg%l(l)%vegt%rootf
       if (sibg%l(l)%pftswap_legacy%phase .eq. 1) then
          call setup_rootf(physcon(hnum)%kroot, physcon(hnum)%rootd, &
               sibg%l(l)%sscolt%dz(1:nsoil), &
               sibg%l(l)%sscolt%layer_z(1:nsoil), &
               sibg%l(l)%sscolt%node_z(1:nsoil), legacy_rootf)
          if (minval(legacy_rootf) .lt. dzero .or. &
              abs(sum(legacy_rootf)-done) .gt. 1.d-5) &
               error stop 'PFTSWAP invalid reconstructed source rooting profile'
          call pftswap_set_dead_layer_metadata(sibg%l(l)%pooldt, legacy_rootf)
       endif
       call pool_het_retr(poolcon(hnum), &
            sibg%l(l)%soilt%zm, sibg%l(l)%soilt%woptzm, &
            sibg%l(l)%soilt%wsat, &
            sibg%gprogt%seas_precip, sibg%gprogt%clim_precip, &
            sibg%l(l)%co2t%clim_assim, sibg%l(l)%co2t%assimd, &
            legacy_rootf, &
            sibg%l(l)%hydrovt%pawfrac_lay(1:nsoil), &
            sibg%l(l)%sscolt%td(1:nsoil), &
            sibg%l(l)%sscolt%satfrac_lay, sibg%l(l)%pooldt, &
            sibg%l(l)%fract)   

       call trace_restart('AFTER_HET_RETR', l)

       !...C13 isotope signature calculations
       call c13_iso_calc(poolcon(hnum),sibg%l(l)%poollt, &
            sibg%l(l)%pooldt, sibg%l(l)%fract)

       call trace_restart('AFTER_C13', l)

       !...C14 isotope signature calculations
       call c14_iso_calc(poolcon(hnum),sibg%l(l)%poollt, &
            sibg%l(l)%pooldt, sibg%l(l)%fract)

       call trace_restart('AFTER_C14', l)

       !...Save diagnostics
       call diagnostic_save_lall(gnum, l, pnum, &
            sibg%l(l)%soilt, &
            sibg%l(l)%cast, sibg%l(l)%co2t,  &
            sibg%l(l)%cost, sibg%l(l)%fract, sibg%l(l)%fluxt, &
            sibg%l(l)%hydrost, sibg%l(l)%hydrovt, &
            sibg%l(l)%phent, sibg%l(l)%pooldt, &
            sibg%l(l)%poollt, sibg%l(l)%radt, &
            sibg%l(l)%sift, sibg%l(l)%sscolt, &
            sibg%l(l)%vegt)

       call trace_restart('AFTER_DIAG_SAVE_LU', l)

   endif !.not. lu_spunup

enddo !l=1,g_nlu

contains

! Prescribed complete shrub removal AFTER the once-daily pool_update has
! settled the fire and carbon gain/loss accumulators. This is additional
! prescribed non-combustion mortality; it does not alter prior fire emissions.
subroutine pftswap_open_transaction(g, lu, land)
    use module_sib, only: lu_type
    implicit none
    integer(i4), intent(in) :: g, lu
    type(lu_type), intent(inout) :: land
    integer(i4) :: sender, receiver, block, si, ri, ls, lr, k, nsrc
    real(r8) :: weight, source_weight_sum, root_total, tol
    real(r8) :: before(3), after(3), incoming(3), leaving(3)
    real(r8) :: adddead(npoollu,nsoil)
    real(r8) :: original(npoolpft,nsoil)
    real(r8) :: source_root(nsoil)

    if (land%ipft .ne. 11 .or. land%pftswap_legacy%phase .ne. 0) &
          error stop 'PFTSWAP attempted opening from incorrect vegetation state'
    if (maxval(abs(land%poollt%poolpft_dgain)) .gt. 1.d-12 .or. &
        maxval(abs(land%poollt%poolpft_dloss)) .gt. 1.d-12 .or. &
        maxval(abs(land%pooldt%poollu_dgain)) .gt. 1.d-12 .or. &
        maxval(abs(land%pooldt%poollu_dloss)) .gt. 1.d-12) &
          error stop 'PFTSWAP pending pool transfers: open boundary is unsafe'
    if (npoolpft .ne. 15 .or. npoollu .ne. 18 .or. ntpool .ne. 33) &
          error stop 'PFTSWAP Step1C expects 5+6 pools, each in C/C13/C14 blocks'
    if (maxval(abs(land%pftswap_legacy%live_lay)) .gt. 0.d0 .or. &
        maxval(abs(land%pftswap_legacy%dead_lay)) .gt. 0.d0) &
          error stop 'PFTSWAP legacy inventories must be empty before opening'

    nsrc = pft_num(11)
    if (nsrc .le. 0) error stop 'PFTSWAP shrub PFT parameter not defined'
    original = land%poollt%poolpft_lay
    if (minval(original) .lt. dzero) &
         error stop 'PFTSWAP negative live inventory before opening'
    source_root = land%vegt%rootf
    root_total = sum(source_root)
    if (root_total .le. 0.d0 .or. minval(source_root) .lt. 0.d0) &
         error stop 'PFTSWAP source rooting distribution invalid'
    adddead = dzero
    before = dzero
    after = dzero
    incoming = dzero
    leaving = dzero

    do block=0,2
      before(block+1) = sum(original(block*5+1:block*5+5,:)) &
          + sum(land%pooldt%poollu_lay(block*6+1:block*6+6,:))
      do sender=1,5
        si = sender + block*5
        ! Source fractions are stored as real(r4) parameters: an r4 SUM
        ! may round to 1 even when the r8 sum of r4 entries differs by
        ! a few 1.e-8.  Check the r8 sum against the intended unity,
        ! then normalize ONLY that harmless parameter roundoff so the
        ! one-time C/C13/C14 inventory transfer closes in real(r8).
        ! Reject inconsistent parameter sets rather than conceal them.
        source_weight_sum = sum(real(poolcon(nsrc)%pool_trans_frac( &
               sender+block*11,6+block*11:11+block*11),r8))
        if (abs(source_weight_sum-done) .gt. 1.d-6) then
           print *, 'PFTSWAP invalid source transfer-fraction sum: ', &
                    sender,source_weight_sum
           error stop 'PFTSWAP cannot conserve residual shrub carbon with these fractions'
        endif
        do receiver=6,11
          ri = receiver-5+block*6
          weight = real(poolcon(nsrc)%pool_trans_frac( &
                    sender+block*11,receiver+block*11),r8) / source_weight_sum
          if (weight .eq. dzero) cycle
          if (weight .lt. dzero .or. weight .gt. done) &
                error stop 'PFTSWAP invalid negative/excess transfer fraction'
          ls = int(pool_indx_lay(sender+block*11))
          lr = int(pool_indx_lay(receiver+block*11))
          if (ls .eq. 1 .and. lr .eq. 1) then
             adddead(ri,1) = adddead(ri,1) + original(si,1)*weight
          elseif (ls .eq. 1 .and. lr .eq. nsoil) then
             do k=1,nsoil
                adddead(ri,k) = adddead(ri,k) &
                   + original(si,1)*weight*source_root(k)/root_total
             enddo
          elseif (ls .eq. nsoil .and. lr .eq. nsoil) then
             adddead(ri,:) = adddead(ri,:) + original(si,:)*weight
          else
             error stop 'PFTSWAP unsupported shrub live-to-dead layer routing'
          endif
        enddo
      enddo
      incoming(block+1) = sum(original(block*5+1:block*5+5,:))
      leaving(block+1) = sum(adddead(block*6+1:block*6+6,:))
      tol = max(1.d-24,1.d-11*abs(incoming(block+1)))
      if (abs(incoming(block+1)-leaving(block+1)) .gt. tol) then
          print *, 'PFTSWAP ROUTING DIAGNOSTIC: year/DOY:',year,doy
          print *, 'PFTSWAP tracer block (1=C,2=C13,3=C14):',block+1
          print '(a,4es24.15)', 'PFTSWAP input/output/residual/tolerance: ', &
               incoming(block+1),leaving(block+1), &
               incoming(block+1)-leaving(block+1),tol
          print '(a,3es24.15)', 'PFTSWAP source-root min/max/sum: ', &
               minval(source_root),maxval(source_root),root_total
          do sender=1,5
              si=sender+block*5
              ls=int(pool_indx_lay(sender+block*11))
              weight=sum(real(poolcon(nsrc)%pool_trans_frac( &
                   sender+block*11,6+block*11:11+block*11),r8))
              print '(a,2i4,4es24.15)', 'PFTSWAP sender/layers, layer1/other/all/weight:', &
                   sender,ls,original(si,1),sum(original(si,2:nsoil)), &
                   sum(original(si,:)),weight
          enddo
          do receiver=1,6
              ri=receiver+block*6
              print '(a,i4,es24.15)', 'PFTSWAP receiving dead pool/added: ', &
                   receiver,sum(adddead(ri,:))
          enddo
          error stop 'PFTSWAP C/C13/C14 shrub-to-dead routing does not close'
      endif
    enddo

    land%pooldt%poollu_lay = land%pooldt%poollu_lay + adddead
    ! Heterotrophic isotope losses at this SAME timestep use the cached
    ! per-layer C13/C14 fractions. Refresh after this one-time inventory
    ! transfer; otherwise they would describe the pre-opening composition.
    call pftswap_refresh_dead_ratios(land%pooldt)
    do k=1,npoollu
       land%pooldt%poollu(k) = sum(land%pooldt%poollu_lay(k,:))
       land%pooldt%poollup(k) = land%pooldt%poollu(k)
    enddo
    land%poollt%poolpft_lay = dzero
    land%poollt%poolpft = dzero
    land%poollt%poolpftp = dzero
    land%poollt%poolpft_dgain = dzero
    land%poollt%poolpft_dloss = dzero
    land%pftswap_legacy%source_pftref = 11
    land%pftswap_legacy%phase = 1

    do block=0,2
       after(block+1) = sum(land%poollt%poolpft_lay(block*5+1:block*5+5,:)) &
          + sum(land%pooldt%poollu_lay(block*6+1:block*6+6,:))
       tol = max(1.d-24,1.d-11*abs(before(block+1)))
       if (abs(before(block+1)-after(block+1)) .gt. tol) &
            error stop 'PFTSWAP total C/C13/C14 budget failed at opening'
    enddo
    print '(a,2i7,a,3es20.11)', 'PFTSWAP_OPEN: year/DOY=',year,doy, &
             ' C/C13/C14 transferred to dead (mol/m2):',incoming
    print '(a,3es20.11)', 'PFTSWAP_OPEN: closure residuals:',after-before

    land%ipft = 1
    call rebuild_pft_derived_state(g,lu)
    ! Preserve soil hydrology/snow temperatures; bare-ground parameters
    ! affect the next normal sib_main call. Complete C3 recruitment is NOT active.
end subroutine pftswap_open_transaction

! Pool provenance in phase 1 is attached to the ONE native dead-pool
! inventory via source_pftref=11. The separate cohort dead_lay is zero.
! Reconstruct shrub-origin layer metadata after bare-PFT rebuilding and
! after restart; do not move any carbon or alter the live bare-PFT profile.
subroutine pftswap_set_dead_layer_metadata(pooldt, rootf_src)
    use module_sib, only: poold_type
    implicit none
    type(poold_type), intent(inout) :: pooldt
    real(r8), intent(in) :: rootf_src(nsoil)
    integer(i4) :: j, ref, block

    do block=0,2
       do j=1,npoollu/3
          ref = j+5+11*block
          pooldt%poollu_flay(j+block*6,:) = dzero
          if (pool_indx_lay(ref) .eq. 1) then
             pooldt%poollu_flay(j+block*6,1) = done
          elseif (pool_indx_lay(ref) .eq. nsoil) then
             pooldt%poollu_flay(j+block*6,:) = rootf_src
          else
             error stop 'PFTSWAP unsupported legacy dead-pool layer count'
          endif
       enddo
    enddo
end subroutine pftswap_set_dead_layer_metadata

! Refresh mixed isotope fractions at the opening instant. Only the
! existing native dead-pool inventory is affected; no carbon is copied
! into the reserved pftswap_legacy arrays. The native isotope routines
! continue to evolve these ratios on subsequent timesteps.
subroutine pftswap_refresh_dead_ratios(pooldt)
    use module_sib, only: poold_type
    implicit none
    type(poold_type), intent(inout) :: pooldt
    integer(i4) :: j, tr, s, ix
    real(r8) :: tc, iso, totc, totiso, ratio
    do tr=1,2
       do j=1,npoollu/3
          ix = j+tr*(npoollu/3)
          totc = dzero
          totiso = dzero
          do s=1,nsoil
             tc = pooldt%poollu_lay(j,s)
             iso = pooldt%poollu_lay(ix,s)
             if (tc .lt. dzero .or. iso .lt. dzero .or. &
                 iso .gt. tc*(1.d0+1.d-10)) &
                  error stop 'PFTSWAP invalid dead carbon isotope inventory'
             ratio = dzero
             if (tc .gt. dzero) ratio = iso/tc
             pooldt%rcpoollu_lay(ix,s) = ratio
             totc = totc+tc
             totiso = totiso+iso
          enddo
          pooldt%rcpoollu(ix) = dzero
          if (totc .gt. dzero) pooldt%rcpoollu(ix) = totiso/totc
       enddo
    enddo
end subroutine pftswap_refresh_dead_ratios

subroutine trace_restart(stage, lnum)

    character(len=*), intent(in) :: stage
    integer(i4), intent(in) :: lnum
    integer(i4) :: j, k, idx_arm
    character(len=600) :: trace_name
    logical :: in_window

    ! Diagnostic restart trace disabled for production runs.
    ! Re-enable for BA exact-continuation testing.
    return

    ! Diagnostic-only V3.2 trace for the Step-0 restart test.
    ! Trace UTC DOY 201 at 00:00-02:00 and 07:00-10:00.
    in_window = (year == 1850 .and. doy == 201 .and. &
         ((sec_day >= 0 .and. sec_day <= 7200) .or. &
          (sec_day >= 25200 .and. sec_day <= 36000)))
    if (.not. in_window) return

    if (.not. trace_opened) then
       write(trace_name,'(a,"restart_trace_rank",i4.4,".dat")') &
            trim(out_path), rank
       open(newunit=trace_unit, file=trim(trace_name), status='replace', &
            action='write')
       write(trace_unit,'(a)') &
            'year doy sec_day sec_tot gnum lnum stage variable i1 i2 value'
       trace_opened = .true.
    endif

    call trace_val(stage, lnum, 'daynew_local', 0, 0, &
         real(merge(1,0,daynew),r8))
    call trace_val(stage, lnum, 'new_day_utc', 0, 0, &
         real(merge(1,0,new_day),r8))
    call trace_val(stage, lnum, 'dayfrac_lst', 0, 0, dayfrac_lst(gnum))
    call trace_val(stage, lnum, 'new_day_lst', 0, 0, &
         real(merge(1,0,new_day_lst(gnum)),r8))

    call trace_val(stage, lnum, 'driver_recnum', 0, 0, real(driver_recnum,r8))
    call trace_val(stage, lnum, 'driver_year', 0, 0, real(driver_year,r8))
    call trace_val(stage, lnum, 'driver_month', 0, 0, real(driver_month,r8))
    call trace_val(stage, lnum, 'driver_day', 0, 0, real(driver_day,r8))
    call trace_val(stage, lnum, 'driver_hour', 0, 0, real(driver_hour,r8))
    call trace_val(stage, lnum, 'driver_seccur', 0, 0, real(driver_seccur,r8))
    call trace_val(stage, lnum, 'driver_secnext', 0, 0, real(driver_secnext,r8))
    call trace_val(stage, lnum, 'driver_readf', 0, 0, &
         real(merge(1,0,driver_readf),r8))
    call trace_val(stage, lnum, 'driver_switchf', 0, 0, &
         real(merge(1,0,driver_switchf),r8))
    call trace_val(stage, lnum, 'driver_updatef', 0, 0, &
         real(merge(1,0,driver_updatef),r8))

    call trace_val(stage, lnum, 'hr_countt', 0, 0, real(hr_countt,r8))
    call trace_val(stage, lnum, 'hr_countsave', 0, 0, real(hr_countsave,r8))
    call trace_val(stage, lnum, 'hr_countwrite', 0, 0, real(hr_countwrite,r8))

    call trace_val(stage, lnum, 'g_tm', 0, 0, sibg%gprogt%tm)
    call trace_val(stage, lnum, 'g_tm1', 0, 0, sibg%gprogt%tm1)
    call trace_val(stage, lnum, 'g_tm2', 0, 0, sibg%gprogt%tm2)
    call trace_val(stage, lnum, 'g_sh', 0, 0, sibg%gprogt%sh)
    call trace_val(stage, lnum, 'g_sh1', 0, 0, sibg%gprogt%sh1)
    call trace_val(stage, lnum, 'g_sh2', 0, 0, sibg%gprogt%sh2)
    call trace_val(stage, lnum, 'g_ps', 0, 0, sibg%gprogt%ps)
    call trace_val(stage, lnum, 'g_ps1', 0, 0, sibg%gprogt%ps1)
    call trace_val(stage, lnum, 'g_ps2', 0, 0, sibg%gprogt%ps2)
    call trace_val(stage, lnum, 'g_spdm', 0, 0, sibg%gprogt%spdm)
    call trace_val(stage, lnum, 'g_sw_dwn', 0, 0, sibg%gprogt%sw_dwn)
    call trace_val(stage, lnum, 'g_dlwbot', 0, 0, sibg%gprogt%dlwbot)
    call trace_val(stage, lnum, 'g_lspr', 0, 0, sibg%gprogt%lspr)
    call trace_val(stage, lnum, 'g_cupr', 0, 0, sibg%gprogt%cupr)
    call trace_val(stage, lnum, 'g_tmd', 0, 0, sibg%gprogt%tmd)
    call trace_val(stage, lnum, 'g_seas_precip', 0, 0, sibg%gprogt%seas_precip)
    call trace_val(stage, lnum, 'g_seas_tm', 0, 0, sibg%gprogt%seas_tm)
    call trace_val(stage, lnum, 'g_clim_cupr', 0, 0, sibg%gprogt%clim_cupr)
    call trace_val(stage, lnum, 'g_clim_precip', 0, 0, sibg%gprogt%clim_precip)
    call trace_val(stage, lnum, 'g_clim_tm', 0, 0, sibg%gprogt%clim_tm)
    call trace_val(stage, lnum, 'g_pco2m', 0, 0, sibg%gprogt%pco2m)
    call trace_val(stage, lnum, 'g_pcosm', 0, 0, sibg%gprogt%pcosm)
    call trace_val(stage, lnum, 'g_cosz', 0, 0, sibg%gdiagt%cosz)
    call trace_val(stage, lnum, 'g_daylen', 0, 0, sibg%gdiagt%daylen)
    call trace_val(stage, lnum, 'g_daylendt', 0, 0, sibg%gdiagt%daylendt)
    call trace_val(stage, lnum, 'g_tmdf', 0, 0, sibg%gdiagt%tmdf)

    call trace_val(stage, lnum, 'ipft', 0, 0, real(sibg%l(lnum)%ipft,r8))
    call trace_val(stage, lnum, 'larea', 0, 0, real(sibg%l(lnum)%larea,r8))
    call trace_val(stage, lnum, 'tc', 0, 0, sibg%l(lnum)%cast%tc)
    call trace_val(stage, lnum, 'eacas', 0, 0, sibg%l(lnum)%cast%eacas)
    call trace_val(stage, lnum, 'shcas', 0, 0, sibg%l(lnum)%cast%shcas)
    call trace_val(stage, lnum, 'tcas', 0, 0, sibg%l(lnum)%cast%tcas)
    call trace_val(stage, lnum, 'tkecas', 0, 0, sibg%l(lnum)%cast%tkecas)

    call trace_val(stage, lnum, 'assim', 0, 0, sibg%l(lnum)%co2t%assim)
    call trace_val(stage, lnum, 'assimd', 0, 0, sibg%l(lnum)%co2t%assimd)
    call trace_val(stage, lnum, 'clim_assim', 0, 0, &
         sibg%l(lnum)%co2t%clim_assim)
    call trace_val(stage, lnum, 'assimpot', 0, 0, sibg%l(lnum)%co2t%assimpot)
    call trace_val(stage, lnum, 'co2cas', 0, 0, sibg%l(lnum)%co2t%co2cas)
    call trace_val(stage, lnum, 'pco2cas', 0, 0, sibg%l(lnum)%co2t%pco2cas)
    call trace_val(stage, lnum, 'rst', 0, 0, sibg%l(lnum)%co2t%rst)
    do j=1,4
       call trace_val(stage, lnum, 'rstfac', j, 0, sibg%l(lnum)%co2t%rstfac(j))
    enddo

    call trace_val(stage, lnum, 'ra', 0, 0, sibg%l(lnum)%fluxt%ra)
    call trace_val(stage, lnum, 'rb', 0, 0, sibg%l(lnum)%fluxt%rb)
    call trace_val(stage, lnum, 'cos_assim', 0, 0, sibg%l(lnum)%cost%cos_assim)
    call trace_val(stage, lnum, 'cos_soil', 0, 0, sibg%l(lnum)%cost%cos_soil)

    ! Explicit isotope continuity diagnostics.
    call trace_val(stage, lnum, 'd13cca', 0, 0, sibg%l(lnum)%fract%d13cca)
    call trace_val(stage, lnum, 'd14cca', 0, 0, sibg%l(lnum)%fract%d14cca)
    call trace_val(stage, lnum, 'c13assim', 0, 0, sibg%l(lnum)%fract%c13assim)
    call trace_val(stage, lnum, 'c14assim', 0, 0, sibg%l(lnum)%fract%c14assim)
    call trace_val(stage, lnum, 'c13assimd', 0, 0, sibg%l(lnum)%fract%c13assimd)
    call trace_val(stage, lnum, 'c14assimd', 0, 0, sibg%l(lnum)%fract%c14assimd)

    call trace_val(stage, lnum, 'capacc_liq', 0, 0, &
         sibg%l(lnum)%hydrost%capacc_liq)
    call trace_val(stage, lnum, 'capacc_snow', 0, 0, &
         sibg%l(lnum)%hydrost%capacc_snow)
    call trace_val(stage, lnum, 'capacg', 0, 0, sibg%l(lnum)%hydrost%capacg)
    call trace_val(stage, lnum, 'snow_gmass', 0, 0, &
         sibg%l(lnum)%hydrost%snow_gmass)
    call trace_val(stage, lnum, 'snow_gdepth', 0, 0, &
         sibg%l(lnum)%hydrost%snow_gdepth)
    call trace_val(stage, lnum, 'snow_cvfc', 0, 0, &
         sibg%l(lnum)%hydrost%snow_cvfc)
    call trace_val(stage, lnum, 'wetfracc', 0, 0, sibg%l(lnum)%hydrost%wetfracc)
    call trace_val(stage, lnum, 'wetfracg', 0, 0, sibg%l(lnum)%hydrost%wetfracg)

    call trace_val(stage, lnum, 'pawfrw', 0, 0, sibg%l(lnum)%hydrovt%pawfrw)
    call trace_val(stage, lnum, 'pawftop', 0, 0, sibg%l(lnum)%hydrovt%pawftop)
    call trace_val(stage, lnum, 'pawfzw', 0, 0, sibg%l(lnum)%hydrovt%pawfzw)
    call trace_val(stage, lnum, 'tawfrw', 0, 0, sibg%l(lnum)%hydrovt%tawfrw)
    call trace_val(stage, lnum, 'tawftop', 0, 0, sibg%l(lnum)%hydrovt%tawftop)
    call trace_val(stage, lnum, 'tawfzw', 0, 0, sibg%l(lnum)%hydrovt%tawfzw)
    call trace_val(stage, lnum, 'clim_pawfrw', 0, 0, &
         sibg%l(lnum)%hydrovt%clim_pawfrw)
    call trace_val(stage, lnum, 'clim_tawfrw', 0, 0, &
         sibg%l(lnum)%hydrovt%clim_tawfrw)

    call trace_val(stage, lnum, 'phenave_assim', 0, 0, &
         sibg%l(lnum)%phent%phenave_assim)
    call trace_val(stage, lnum, 'phenave_assimpot', 0, 0, &
         sibg%l(lnum)%phent%phenave_assimpot)
    call trace_val(stage, lnum, 'phenflag_assimlow', 0, 0, &
         real(merge(1,0,sibg%l(lnum)%phent%phenflag_assimlow),r8))
    call trace_val(stage, lnum, 'phenave_pr', 0, 0, &
         sibg%l(lnum)%phent%phenave_pr)
    call trace_val(stage, lnum, 'phenave_prpot', 0, 0, &
         sibg%l(lnum)%phent%phenave_prpot)
    call trace_val(stage, lnum, 'phenave_tm', 0, 0, &
         sibg%l(lnum)%phent%phenave_tm)
    call trace_val(stage, lnum, 'phenave_tawftop', 0, 0, &
         sibg%l(lnum)%phent%phenave_tawftop)
    call trace_val(stage, lnum, 'phenave_env', 0, 0, &
         sibg%l(lnum)%phent%phenave_env)
    call trace_val(stage, lnum, 'phenave_wa', 0, 0, &
         sibg%l(lnum)%phent%phenave_wa)
    call trace_val(stage, lnum, 'phenave_wac', 0, 0, &
         sibg%l(lnum)%phent%phenave_wac)
    call trace_val(stage, lnum, 'phens_grw', 0, 0, sibg%l(lnum)%phent%phens_grw)
    call trace_val(stage, lnum, 'phens_wx', 0, 0, sibg%l(lnum)%phent%phens_wx)
    call trace_val(stage, lnum, 'phen_pi', 0, 0, sibg%l(lnum)%phent%phen_pi)
    call trace_val(stage, lnum, 'phen_istage', 0, 0, &
         real(sibg%l(lnum)%phent%phen_istage,r8))
    call trace_val(stage, lnum, 'nd_gs', 0, 0, real(sibg%l(lnum)%phent%nd_gs,r8))
    call trace_val(stage, lnum, 'nd_dormant', 0, 0, &
         real(sibg%l(lnum)%phent%nd_dormant,r8))

    call trace_val(stage, lnum, 'lai', 0, 0, sibg%l(lnum)%vegt%lai)
    call trace_val(stage, lnum, 'lait', 0, 0, sibg%l(lnum)%vegt%lait)
    call trace_val(stage, lnum, 'fpar', 0, 0, sibg%l(lnum)%vegt%fpar)
    call trace_val(stage, lnum, 'vcover', 0, 0, sibg%l(lnum)%vegt%vcover)
    call trace_val(stage, lnum, 'green', 0, 0, sibg%l(lnum)%vegt%green)
    call trace_val(stage, lnum, 'clim_lai', 0, 0, sibg%l(lnum)%vegt%clim_lai)
    call trace_val(stage, lnum, 'vmax', 0, 0, sibg%l(lnum)%vegt%vmax)
    call trace_val(stage, lnum, 'z0', 0, 0, sibg%l(lnum)%vegt%z0)
    call trace_val(stage, lnum, 'zp_dispd', 0, 0, sibg%l(lnum)%vegt%zp_dispd)

    call trace_val(stage, lnum, 'resp_auto', 0, 0, &
         sibg%l(lnum)%poollt%resp_auto)
    call trace_val(stage, lnum, 'resp_leaf', 0, 0, &
         sibg%l(lnum)%poollt%resp_leaf)
    call trace_val(stage, lnum, 'resp_root', 0, 0, &
         sibg%l(lnum)%poollt%resp_root)
    call trace_val(stage, lnum, 'resp_grow', 0, 0, &
         sibg%l(lnum)%poollt%resp_grow)
    call trace_val(stage, lnum, 'resp_mntn', 0, 0, &
         sibg%l(lnum)%poollt%resp_mntn)
    call trace_val(stage, lnum, 'resp_het', 0, 0, sibg%l(lnum)%pooldt%resp_het)
    call trace_val(stage, lnum, 'resp_soil', 0, 0, sibg%l(lnum)%pooldt%resp_soil)
    call trace_val(stage, lnum, 'mhrt_soil_hot', 0, 0, &
         sibg%l(lnum)%pooldt%mhrt_soil_hot)

    call trace_val(stage, lnum, 'nsl', 0, 0, real(sibg%l(lnum)%sscolt%nsl,r8))
    do j=lbound(sibg%l(lnum)%sscolt%td,1),ubound(sibg%l(lnum)%sscolt%td,1)
       call trace_val(stage, lnum, 'td', j, 0, sibg%l(lnum)%sscolt%td(j))
       call trace_val(stage, lnum, 'www_liq', j, 0, sibg%l(lnum)%sscolt%www_liq(j))
       call trace_val(stage, lnum, 'www_ice', j, 0, sibg%l(lnum)%sscolt%www_ice(j))
    enddo

    do j=1,size(sibg%l(lnum)%poollt%poolpft)
       call trace_val(stage, lnum, 'poolpft', j, 0, &
            sibg%l(lnum)%poollt%poolpft(j))
       call trace_val(stage, lnum, 'poolpftp', j, 0, &
            sibg%l(lnum)%poollt%poolpftp(j))
    enddo
    do j=1,size(sibg%l(lnum)%pooldt%poollu)
       call trace_val(stage, lnum, 'poollu', j, 0, &
            sibg%l(lnum)%pooldt%poollu(j))
       call trace_val(stage, lnum, 'poollup', j, 0, &
            sibg%l(lnum)%pooldt%poollup(j))
    enddo

    if (hr_saveluf .and. allocated(hr_lu)) then
       idx_arm = 0
       do k=1,hr_nvarlu
          if (hr_vreflu(k) == 487) idx_arm = k
       enddo
       if (idx_arm > 0 .and. hr_countsave >= 1 .and. &
            hr_countsave <= size(hr_lu,4)) then
          call trace_val(stage, lnum, 'hr_arm_gain_accum', 0, 0, &
               hr_lu(gnum,lnum,idx_arm,hr_countsave))
       endif
    endif

end subroutine trace_restart

subroutine trace_val(stage, lnum, vname, i1, i2, value)

    character(len=*), intent(in) :: stage, vname
    integer(i4), intent(in) :: lnum, i1, i2
    real(r8), intent(in) :: value

    write(trace_unit,'(i4,1x,i3,1x,i6,1x,i12,1x,i4,1x,i3,1x,a,1x,a,1x,i3,1x,i3,1x,es25.17e3)') &
         year, doy, sec_day, sec_tot, gnum, lnum, trim(stage), trim(vname), i1, i2, value

end subroutine trace_val

end subroutine sib_control
