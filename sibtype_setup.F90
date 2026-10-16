! Sets up the SiB4 diagnostic variables for a simulation
subroutine sibtype_setup()

use kinds

use module_io, only: &
     ic_restart, minveg_diag, minveg_clim_policy
use module_oparams, only: &
     cos_casd_min
use module_pparams, only: &
     pi180, tref
use module_param
use module_pftinfo, only: &
     pft_num, pft_type, &
     type_bare, type_crop, type_grass
use module_phosib, only: co2_casd_min
use module_sib, only: sib
use module_sibconst, only: &
     nsoil, subcount, subset, &
     sublonsib, sublatsib, &
     tan_dec
use module_sibvs, only: sibvs
use module_time, only: &
     daylenmax, starttime

implicit none

!local variables
integer(i4) :: gref
real(r4) :: glon, glat
integer(i4) :: g,l

integer(byte) :: ptyperef
integer(i4) :: dnum, pref, pnum
logical :: isbare, iscrop, isgrass


!------------------------------------------
print*,'Setting SiB4 Variables'

! Set the diagnostic variables
do g=1, subcount
   gref = subset(g)
   glon = sublonsib(g)
   glat = sublatsib(g)

   !...Daylength
   call day_lengthpt(pi180, tan_dec, glat, &
        sib%g(g)%gdiagt%daylen, sib%g(g)%gdiagt%daylendt)

   IF (abs(sib%g(g)%gdiagt%daylendt) .lt. 1.E-6) THEN
       call day_lengthpt(pi180, tan_dec, glat, &
            sib%g(g)%gdiagt%daylen, sib%g(g)%gdiagt%daylendt)
   ENDIF

   do l=1, sib%g(g)%g_nlu

      pref = sib%g(g)%l(l)%ipft
      pnum = pft_num(pref)
      ptyperef = pft_type(pnum)
      isbare = (ptyperef .eq. type_bare)
      iscrop = (ptyperef .eq. type_crop)
      isgrass = (ptyperef .eq. type_grass)

      !...Time Invariant
      !Set soil properties
      call setup_soilt( &
          sibvs(gref)%clayfrac(l),  &
          sibvs(gref)%sandfrac(l),  &
          sibvs(gref)%soref_vis(l), &
          sibvs(gref)%soref_nir(l), &
          sib%g(g)%l(l)%soilt)

      !Set soil column variables.  This is geometry/physical-state setup and
      !is not part of a runtime PFT transition rebuild.
      call setup_sscolt( &
          sib%g(g)%l(l)%soilt%poros, &
          sib%g(g)%l(l)%sscolt)

      !Set the quantities that depend on the active PFT but can be rebuilt
      !without resetting prognostic soil/snow state, carbon pools, or
      !phenology/history state.  This same helper is the entry point for
      !future runtime PFT activation.
      call rebuild_pft_derived_state(g, l)

      !...Initial-condition dependent setup
      call setup_poolt( poolcon(pnum), &
           sib%g(g)%l(l)%equibdt, sib%g(g)%l(l)%pooldt, &
           sib%g(g)%l(l)%equiblt, sib%g(g)%l(l)%poollt, &
           sib%g(g)%l(l)%fract, ic_restart)

      ! Diagnostic branch C: remove equilibrated live biomass only, leaving
      ! soil/dead pools and prognostic physical state untouched.  This is
      ! deliberately NOT a carbon-conserving ecological transition: printed
      ! removed C inventories define the externally imposed IC adjustment.
      if (minveg_diag) then
          call reset_minveg_diagnostic(g, l, pnum)
      endif

      if (.not. ic_restart) then
         ! Fresh spinup/requib initialization retains the historical full
         ! setup behavior.
         sib%g(g)%l(l)%cost%cos_s = 2d-8 ! ~500 ppt

         call setup_hydrost( physcon(pnum)%z1, &
              physcon(pnum)%z2, &
              sib%g(g)%l(l)%sscolt, &
              sib%g(g)%l(l)%hydrost)
         call hydro_sets(sib%g(g)%l(l)%soilt, &
              sib%g(g)%l(l)%hydrost, sib%g(g)%l(l)%sscolt)

         call hydro_setv(  &
              sib%g(g)%l(l)%soilt%fc_eff, &
              sib%g(g)%l(l)%soilt%wp_eff, &
              sib%g(g)%l(l)%vegt%rootf, &
              sib%g(g)%l(l)%sscolt, sib%g(g)%l(l)%hydrovt)

         dnum = MAX(starttime-1, 0)
         call veg_update( dnum, &
             gref, glon, glat, &
             pnum, pref, iscrop, isgrass, &
             physcon(pnum), &
             sib%g(g)%l(l)%pooldt%poollu, &
             sib%g(g)%l(l)%hydrost%snow_cvfc, &
             sib%g(g)%l(l)%sscolt%node_z(1:nsoil), &
             sib%g(g)%l(l)%poollt, sib%g(g)%l(l)%vegt)

         if (.not. isbare) then
             call setup_phent( &
                  g, l, pnum, pref, &
                  phencon(pnum), poolcon(pnum),  &
                  sib%g(g)%gdiagt%daylen, sib%g(g)%gdiagt%daylendt, &
                  daylenmax(g), &
                  sib%g(g)%l(l)%vegt%lai, sib%g(g)%gprogt%tm, &
                  sib%g(g)%l(l)%vegt%vmax, sib%g(g)%l(l)%pooldt, &
                  sib%g(g)%l(l)%poollt, sib%g(g)%l(l)%phent, &
                  physcon(pnum), .false.)
             if (minveg_diag) then
                ! C3 check: show whether setup_phent/phen_gss altered the
                ! same-PFT equilibrated weather-reference histories.
                associate(ph => sib%g(g)%l(l)%phent)
                   write(*,'(a,3(1x,i4),4(1x,es16.8))') &
                     'MINVEG_C3_AFTER_SETUP g,l,pft,env,wa,wac,wacsm', &
                     g,l,pnum,ph%phenave_env,ph%phenave_wa, &
                     ph%phenave_wac,ph%phenave_wacsm
                end associate
                if (sib%g(g)%l(l)%phent%phen_istage .ne. 1) then
                   print *, 'MINVEG ERROR: startup phenology is not stage 1',g,l,pnum, &
                        sib%g(g)%l(l)%phent%phen_istage
                   stop 1
                endif
                print *, 'MINVEG startup stage=1, g/l/PFT=',g,l,pnum, &
                     ' vmax=',sib%g(g)%l(l)%vegt%vmax, &
                     ' alloc=',sib%g(g)%l(l)%poollt%alloc_phen(:)
             endif
          endif
      endif
   enddo !l=1, g_nlu

   !Set grid cell diagnostic variables
   sib%g(g)%gdiagt%tmdf = 32.0 &
       + (sib%g(g)%gprogt%tmd-tref * 1.8)

enddo !g=1, subcount


end subroutine sibtype_setup


! One-time minimum-live-vegetation diagnostic for a new transient launched
! from equilibrated initial conditions. Not used by exact transient restart.
! This helper never changes source PFT, LU fraction, soil/dead C or forcing.
subroutine reset_minveg_diagnostic(g, l, pnum)
use kinds
use module_io, only: minveg_clim_policy
use module_param, only: poolcon
use module_sib, only: sib
use module_sibconst, only: npoolpft
implicit none
integer(i4), intent(in) :: g, l, pnum
real(r8) :: removed(3)
integer(i4) :: n, j, k
real(r8) :: layer_total

n = npoolpft/3
removed(1) = sum(sib%g(g)%l(l)%poollt%poolpft(1:n)) &
           - sum(poolcon(pnum)%poolpft_min(1:n))
removed(2) = sum(sib%g(g)%l(l)%poollt%poolpft(n+1:2*n)) &
           - sum(poolcon(pnum)%poolpft_min(n+1:2*n))
removed(3) = sum(sib%g(g)%l(l)%poollt%poolpft(2*n+1:3*n)) &
           - sum(poolcon(pnum)%poolpft_min(2*n+1:3*n))

! pool_update evolves poolpft_lay and then replaces poolpft with its sum.
! Resetting only poolpft/poolpftp therefore reintroduces equilibrated biomass
! on the first day.  Initialize all layered inventories consistently, for
! total C and both isotopes, using the already-built PFT layer fractions.
associate(pl => sib%g(g)%l(l)%poollt, eq => sib%g(g)%l(l)%equiblt)
   pl%poolpft(:) = poolcon(pnum)%poolpft_min(:)
   pl%poolpftp(:) = poolcon(pnum)%poolpft_min(:)
   do j=1,npoolpft
      pl%poolpft_lay(j,:) = pl%poolpft(j) * pl%poolpft_flay(j,:)
      layer_total = sum(pl%poolpft_lay(j,:))
      if (abs(layer_total-pl%poolpft(j)) > &
          1.d-11*max(1.d0,abs(pl%poolpft(j)))) then
         print *, 'ERROR MINVEG layer distribution mismatch',g,l,pnum,j, &
              layer_total,pl%poolpft(j)
         stop 1
      endif
      eq%poolpft_init(j) = pl%poolpft(j)
      eq%poolpft_min(j) = pl%poolpft(j)
      eq%poolpft_max(j) = pl%poolpft(j)
   enddo
   ! Previously equilibrated isotope ratios and working pool snapshots
   ! must not describe the removed biomass.
   pl%rcpoolpft(:) = 0.d0
   pl%rcpoolpft_lay(:,:) = 0.d0
   pl%curpoolpft(:) = pl%poolpft(:)
   do k=1,2
      do j=1,n
         if (pl%poolpft(j) > 0.d0) then
            pl%rcpoolpft(j+k*n) = pl%poolpft(j+k*n)/pl%poolpft(j)
         endif
         where (pl%poolpft_lay(j,:) > 0.d0)
            pl%rcpoolpft_lay(j+k*n,:) = &
                 pl%poolpft_lay(j+k*n,:)/max(pl%poolpft_lay(j,:),tiny(1.d0))
         end where
      enddo
   enddo
end associate
! C3 same-PFT sensitivity: retain that PFT's equilibrated
! environmental/available-water histories and seasonal combined maximum.
! This tests severe biomass reduction of an existing PFT; it is NOT an
! initialization rule for establishment of a different PFT. Keep the
! stage-1 diagnostic treatment and other C2 reset behavior unchanged.
associate(ph => sib%g(g)%l(l)%phent)
   write(*,'(a,3(1x,i4),4(1x,es16.8))') &
     'MINVEG_C3_BEFORE_RESET g,l,pft,env,wa,wac,wacsm', &
     g,l,pnum,ph%phenave_env,ph%phenave_wa, &
     ph%phenave_wac,ph%phenave_wacsm
   ph%phen_istage = 1
   ph%phen_pi = 1.d0
   ph%phens_dayl = 1.d0
   ph%phenave_assim = 0.d0
   ph%phenave_assimsm = 0.d0
   ph%phenave_assimpot = 0.d0
   ph%phenave_pr = 0.d0
   ph%phenave_prsm = 0.d0
   ph%phenave_prsdoy = 0.d0
   ph%phenave_prcdoy = 0.d0
   ph%phenave_prpot = 0.d0
   ph%phenave_tawftop = 0.d0
   ph%phenave_tm = 0.d0
   ! Keep ph%phenave_env, ph%phenave_wa and ph%phenave_wacsm
   ! from the same-PFT requib IC; rebuild derived wac with the
   ! verified corrected wa+env expression, never wac+env.
   ph%phenave_wac = 0.5d0*(ph%phenave_wa + ph%phenave_env)
   write(*,'(a,3(1x,i4),4(1x,es16.8))') &
     'MINVEG_C3_AFTER_RESET g,l,pft,env,wa,wac,wacsm', &
     g,l,pnum,ph%phenave_env,ph%phenave_wa, &
     ph%phenave_wac,ph%phenave_wacsm
   ph%phenc_climp = 0.d0
   ph%phenc_laimax = 0.d0
   ph%phenc_laimin = 0.d0
   ph%phens_grw = 0.d0
   ph%phens_wx = 0.d0
   ph%phenflag_assimlow = .false.
   ph%phenflag_precip = .false.
   ph%phenflag_moist = .false.
   ph%phenflag_temp = .false.
   ph%phenflag_daylen = .false.
   ph%phenflag_gsspass = .false.
   ph%nd_dormant = 0
   ph%nd_gs = 0
   if (allocated(ph%nd_stg)) ph%nd_stg(:) = 0
   ph%ipd = 0
   ph%dapd = 0
   ph%dapdaf = 0
   ph%gdd = 0.d0
   ph%seed_pool = 0.d0
end associate

if (trim(minveg_clim_policy) == 'zero') then
   sib%g(g)%l(l)%vegt%clim_lai = 0.d0
   sib%g(g)%l(l)%co2t%clim_assim = 0.d0
   sib%g(g)%l(l)%co2t%assimd = 0.d0
endif
print '(a,3i5,a,3es23.14)', &
 'MINVEG IC adjustment (g,l,pft): ',g,l,pnum, &
 ' removed total-C/C13/C14 (mol m-2 PFT): ', removed(:)
end subroutine reset_minveg_diagnostic


!------------------------------------------
! SET-UP TIME-INVARIANT VARIABLES
!------------------------------------------

!----------------------------------------
!Routine to set the soil variables
!
!...These are set using %sand and %clay,
!...using relationships from:
!   Clapp and Hornberger, 1978; 
!   Cosby et al., 1984; and
!   Lawrence and Slater, 2008
subroutine setup_soilt( &
     clayfrac, sandfrac, &
     soref_vis, soref_nir, &
     soilt)

use kinds
use module_pparams, only: &
    wpotfc, wpotwp
use module_sib, only: &
    soil_type
use module_sibconst, only: nsoil, &
    print_soil, print_stop

implicit none

!...input variables
real(r8), intent(in) :: clayfrac, sandfrac
real(r8), intent(in) :: soref_vis, soref_nir
type(soil_type), intent(inout) :: soilt

!...local variables
integer(i4) :: k
real(r8) :: cfrac, sfrac
real(r8) :: tkmsoil  ! soil mineral conductivity (W/m K)
real(r8) :: bd       ! bulk density of dry soil material (kg/m3)

     !...basic soil properties
     soilt%clayfrac = clayfrac
     soilt%sandfrac = sandfrac
     soilt%soref_vis = soref_vis
     soilt%soref_nir = soref_nir

     cfrac = clayfrac*100.
     sfrac = sandfrac*100.

     soilt%poros = 0.489-0.00126*sfrac
     soilt%satco = 0.0070556*10**(-0.884+0.0153*sfrac)/1000.

     !...misc soil variables with depth
     soilt%csolid = (2.128*sandfrac + 2.385*clayfrac) / &
                      (sandfrac+clayfrac)*1.0E6

     bd = (1.0 - soilt%poros) * 2.7E3
     soilt%tkdry = (0.135*bd + 64.7) / (2.7E3 - 0.947*bd)

     tkmsoil = (8.8*sandfrac + 2.92*clayfrac ) / &
               (sandfrac + clayfrac)
     soilt%tkmg = tkmsoil**(1.0-soilt%poros) 
     soilt%tksat = soilt%tkmg*0.57**soilt%poros

     soilt%bee = 2.91+0.159*cfrac
     soilt%phsat = -10.*10.**(1.88-0.0131*sfrac)/1000.
     soilt%fieldcap = soilt%poros * &
            ((wpotfc/9.8)/soilt%phsat) ** (-1.0/soilt%bee)
     soilt%vwcmin = soilt%poros * &
            ((wpotwp/9.8)/soilt%phsat) ** (-1.0/soilt%bee)

     soilt%wopt = (-0.08*clayfrac**2 + &
                    0.22*clayfrac+0.59)*100.
     soilt%wsat = 0.25*clayfrac+0.5
     soilt%zm = -2*clayfrac**3 - 0.4491*clayfrac**2 + &
                    0.2101*clayfrac+0.3478
     soilt%woptzm = (soilt%wopt/100.)**soilt%zm

     !...print if requested
     if (print_soil) then 
         print*,''
         print('(a,2F6.3)'),' Clay/Sand Fracs: ', clayfrac, sandfrac
         print('(a,2F6.3)'),' Soil Reflectivity (vis/nir): ', &
              soref_vis, soref_nir
         print('(a,F6.3)'),' Soil Tension at Sat: ', soilt%phsat
         print('(a,F6.3)'),' Porosity: ',soilt%poros
         print('(a,2F6.3)'),' Field Capacity/Wilting Point: ', &
              soilt%fieldcap, soilt%vwcmin
         print*,''

         if (print_stop) stop
     endif

end subroutine setup_soilt


!----------------------------------------
!Routine to set the soil column variables
subroutine setup_sscolt(poros, sscolt)

use kinds
use module_oparams, only: wsat_default
use module_pparams, only: denh2o
use module_sib, only: sscol_type
use module_sibconst, only: &
    nsnow, nsoil, &
    print_sscol, print_stop

implicit none

!...input variables
real(r8), intent(in) :: poros
type(sscol_type), intent(inout) :: sscolt

!...local variables
integer(i4) :: j,k
real(r8) :: scalez
real(r8) :: www_tot


      !...chose soil layer vertical scaling factor
      !scalez = 0.025 !...used for SiBCrop
      scalez = 0.073

      !...soil layers
      do k=1, nsoil
          sscolt%node_z(k) = (scalez*(exp(0.5*(k-0.5))-1.0))
      enddo

      sscolt%dz(1) = 0.5*(sscolt%node_z(1) + sscolt%node_z(2))
      do k=2,nsoil-1
         sscolt%dz(k) = 0.5*(sscolt%node_z(k+1) - sscolt%node_z(k-1))
      enddo
      sscolt%dz(nsoil) = sscolt%node_z(nsoil) - sscolt%node_z(nsoil-1)

      do k=1,nsoil-1
          sscolt%layer_z(k) = 0.5*(sscolt%node_z(k) + sscolt%node_z(k+1))
      enddo
      sscolt%layer_z(nsoil) = sscolt%node_z(nsoil) + 0.5*sscolt%dz(nsoil)

      !...snow layers
      do j=0,sscolt%nsl+1,-1
         sscolt%node_z(j) = sscolt%layer_z(j) - 0.5*sscolt%dz(j)
         sscolt%layer_z(j-1) = sscolt%node_z(j) - 0.5*sscolt%dz(j)
      enddo

      !...snow/soil column moisture
      www_tot = sum(sscolt%www_liq(1:nsoil)) + sum(sscolt%www_ice(1:nsoil))
      if (www_tot .le. 1.E-12) then
          !print*,'   Saturated Soil Moisture'
          do k=sscolt%nsl+1,nsoil
              sscolt%www_liq(k) = &
                  sscolt%dz(k) * denh2o * poros * wsat_default
          enddo
      endif

      !...print soil/snow layer info
      IF (print_sscol) THEN
          print*, '     Number of Snow Layers: ', abs(sscolt%nsl)
          print'(2a)','        lev   layer (m)  node (m)  dz (m)      ', &
                      'td (K)      liq (kg/m2)    ice (kg/m2)'
          do k=-nsnow+1,nsoil
              print'(a,i6,3f10.5,3f14.5)','     ',k, &
                  sscolt%layer_z(k), sscolt%node_z(k), &
                  sscolt%dz(k), sscolt%td(k), &
                  sscolt%www_liq(k), sscolt%www_ice(k)
          enddo
          if (print_stop) stop
      ENDIF

end subroutine setup_sscolt


!----------------------------------------
!Routine to set the vegetation rooting fraction
subroutine setup_rootf( &
      kroot, rootd, &
      dz, layer_z, node_z, rootf)

use kinds
use module_sibconst, only: nsoil

implicit none

!...input variables
real(r4), intent(in) :: kroot, rootd
real(r8), dimension(nsoil), intent(in) :: &
      dz, layer_z, node_z
real(r8), dimension(nsoil), intent(inout) :: &
      rootf

!...local variables
integer(i4) :: k
real(r8) :: rtot
real(r8) :: zbot, ztop
real(r8) :: totalroot

      !...calculate root distribution
      ztop = dzero
      if (kroot > 0.) then
          totalroot = (1.0 - exp(-kroot*layer_z(nsoil)))/kroot

          do k=1,nsoil
              zbot = ztop + dz(k)
              rootf(k) = (exp(-kroot*ztop) - exp(-kroot*zbot)) &
                          / (kroot*totalroot)
              ztop = zbot
          enddo
      endif

      !...modify root distribution
      if (rootd < layer_z(nsoil)) then
          rtot = dzero

          do k=1,nsoil
              if (node_z(k) > rootd) then
                  rootf(k) = dzero
              else
                  rtot = rtot + rootf(k)
              endif
          enddo

          if (rtot > 0.) then
              do k=1,nsoil
                  rootf(k) = rootf(k)/rtot
              enddo
          endif
      endif

end subroutine setup_rootf


!----------------------------------------
! Set PFT-dependent soil-water thresholds only.
!
! The underlying soil hydraulic properties are site/soil properties and
! remain unchanged during a vegetation transition.  Only the effective
! field-capacity and wilting-point thresholds depend on the active PFT.
subroutine setup_soil_pft(fc_min, wp_min, soilt)

use kinds
use module_sib, only: soil_type
use module_sibconst, only: nsoil

implicit none

real(r4), dimension(nsoil/2), intent(in) :: fc_min, wp_min
type(soil_type), intent(inout) :: soilt

integer(i4) :: count, k

count=1
do k=1,nsoil
   soilt%fc_eff(k) = MIN(soilt%fieldcap, fc_min(count))
   soilt%wp_eff(k) = MIN(soilt%vwcmin, wp_min(count))
   if (mod(k,2) .eq. 0) count=count+1
enddo

end subroutine setup_soil_pft


!----------------------------------------
! Rebuild carbon-pool vertical-fraction metadata for the active PFT.
!
! This routine changes only fraction metadata.  It does not move carbon,
! change pool masses, or alter previous/equilibrium pool state.  That
! separation is important for the future legacy-carbon bookkeeping.
subroutine setup_pool_flay(rootf, pooldt, poollt)

use kinds
use module_poolinfo, only: pool_indx_lay
use module_sib, only: poold_type, pooll_type
use module_sibconst, only: npoollu, npoolpft, nsoil

implicit none

real(r8), dimension(nsoil), intent(in) :: rootf
type(poold_type), intent(inout) :: pooldt
type(pooll_type), intent(inout) :: poollt

integer(i4) :: p, dlay

!...total-C dead pools
do p=1,npoollu/3
   dlay = pool_indx_lay(p+npoolpft/3)
   pooldt%poollu_flay(p,:) = dzero
   if (dlay .eq. 1) then
      pooldt%poollu_flay(p,1) = done
   else
      pooldt%poollu_flay(p,:) = rootf(:)
   endif
enddo

!...total-C live pools
do p=1,npoolpft/3
   dlay = pool_indx_lay(p)
   poollt%poolpft_flay(p,:) = dzero
   if (dlay .eq. 1) then
      poollt%poolpft_flay(p,1) = done
   else
      poollt%poolpft_flay(p,:) = rootf(:)
   endif
enddo

!...C13 dead pools
do p=npoollu/3+1,2*npoollu/3
   dlay = pool_indx_lay(p+2*npoolpft/3)
   pooldt%poollu_flay(p,:) = dzero
   if (dlay .eq. 1) then
      pooldt%poollu_flay(p,1) = done
   else
      pooldt%poollu_flay(p,:) = rootf(:)
   endif
enddo

!...C13 live pools
do p=npoolpft/3+1,2*npoolpft/3
   dlay = pool_indx_lay(p+npoolpft/3+1)
   poollt%poolpft_flay(p,:) = dzero
   if (dlay .eq. 1) then
      poollt%poolpft_flay(p,1) = done
   else
      poollt%poolpft_flay(p,:) = rootf(:)
   endif
enddo

!...C14 dead pools
do p=2*npoollu/3+1,npoollu
   dlay = pool_indx_lay(p+npoolpft)
   pooldt%poollu_flay(p,:) = dzero
   if (dlay .eq. 1) then
      pooldt%poollu_flay(p,1) = done
   else
      pooldt%poollu_flay(p,:) = rootf(:)
   endif
enddo

!...C14 live pools
do p=2*npoolpft/3+1,npoolpft
   dlay = pool_indx_lay(p+2*npoolpft/3+2)
   poollt%poolpft_flay(p,:) = dzero
   if (dlay .eq. 1) then
      poollt%poolpft_flay(p,1) = done
   else
      poollt%poolpft_flay(p,:) = rootf(:)
   endif
enddo

end subroutine setup_pool_flay


!----------------------------------------
! Rebuild quantities that are derived from the currently active PFT and are
! safe to change without resetting prognostic physical or carbon state.
!
! This intentionally does NOT modify:
!   * soil/snow temperature or water
!   * carbon-pool masses or previous-pool state
!   * phenology/history variables
!   * LAI/FPAR/aerodynamic canopy state
!   * climatological LAI or assimilation
!
! The canopy/phenology activation step will be handled separately when the
! destination PFT is actually established.
subroutine rebuild_pft_derived_state(g, l)

use kinds
use module_oparams, only: cos_casd_min
use module_param, only: physcon
use module_pftinfo, only: pft_num
use module_phosib, only: co2_casd_min
use module_sib, only: sib
use module_sibconst, only: nsoil

implicit none

integer(i4), intent(in) :: g, l
integer(i4) :: pref, pnum

pref = sib%g(g)%l(l)%ipft
pnum = pft_num(pref)

if (pnum .le. 0) then
   print*,'ERROR: rebuild_pft_derived_state invalid PFT ref at g/l: ', &
          g, l, pref
   stop
endif

!...PFT-dependent water-availability thresholds.
call setup_soil_pft( &
     physcon(pnum)%fc_min, physcon(pnum)%wp_min, &
     sib%g(g)%l(l)%soilt)

!...PFT-dependent rooting profile using the existing soil-column geometry.
call setup_rootf( &
     physcon(pnum)%kroot, physcon(pnum)%rootd, &
     sib%g(g)%l(l)%sscolt%dz(1:nsoil), &
     sib%g(g)%l(l)%sscolt%layer_z(1:nsoil), &
     sib%g(g)%l(l)%sscolt%node_z(1:nsoil), &
     sib%g(g)%l(l)%vegt%rootf)

!...PFT-dependent canopy-air-space depths.
sib%g(g)%l(l)%co2t%casd = MAX(dble(co2_casd_min), &
                               dble(physcon(pnum)%z2))
sib%g(g)%l(l)%cost%cos_casd = MAX(dble(cos_casd_min), &
                                   dble(physcon(pnum)%z2))

!...Vertical-fraction metadata only; no carbon is moved here.
call setup_pool_flay( &
     sib%g(g)%l(l)%vegt%rootf, &
     sib%g(g)%l(l)%pooldt, sib%g(g)%l(l)%poollt)

end subroutine rebuild_pft_derived_state


!------------------------------------------
! SET-UP TIME-STEP VARYING VARIABLES
!------------------------------------------

!-----------------------------------------------
!Routine to set hydrological snow/soil variables
subroutine setup_hydrost( z1, z2, &
      sscolt, hydrost)
 
use kinds
use module_pparams, only:  &
     cwlim, gwlim, gwctog, &
     denice
use module_sib, only: &
     soil_type, sscol_type, &
     hydros_type
use module_sibconst, only: nsnow

implicit none

!...input variables
real(r4), intent(in) :: z1, z2
type(sscol_type), intent(in) :: sscolt
type(hydros_type), intent(inout) :: hydrost

!...local variables
integer(i4) :: k

      hydrost%satcapc = cwlim
      hydrost%satcapg = gwlim
      hydrost%wetfracc = &
          MAX(dzero, MIN(done, &
              (hydrost%capacc_liq + &
               hydrost%capacc_snow) / &
               hydrost%satcapc))
      hydrost%wetfracg = &
          MAX(dzero, MIN(done, gwctog * &
              (hydrost%capacg / &
               hydrost%satcapg)))

      do k=-nsnow+1,0
          hydrost%snow_gmass = hydrost%snow_gmass + &
              sscolt%www_liq(k) + sscolt%www_ice(k)
          hydrost%snow_gdepth = hydrost%snow_gdepth + &
              sscolt%dz(k)
      enddo

      hydrost%www_tot = sum(sscolt%www_liq) + &
          sum(sscolt%www_ice)

      hydrost%snow_cvfc = &
          ((hydrost%capacc_snow/denice)*5. - z1) &
            / (z2-z1)
      hydrost%snow_cvfc = MIN(done, &
          MAX(dzero, hydrost%snow_cvfc))

end subroutine setup_hydrost


!-------------------------------------
!Routine to set the pool variables
subroutine setup_poolt( poolcont, &
      equibdt, pooldt, &
      equiblt, poollt, fract, continuation_restart)
 
use kinds
use module_poolinfo, only: &
    pool_indx_lay
use module_sib, only: &
    equibd_type, equibl_type, &
    poold_type, pooll_type, &
    fract_type
use module_sibconst, only: &
    npoollu, npoolpft, nsoil, &
    spinup, spinup_default, &
    spinup_continue
use module_param, only: &
    pool_param

implicit none

!...input variables
!real(r8), dimension(npoolpft), intent(in) :: poolpft_min
type(pool_param), intent(in) :: poolcont
type(equibd_type), intent(inout) :: equibdt
type(poold_type), intent(inout) :: pooldt
type(equibl_type), intent(inout) :: equiblt
type(pooll_type), intent(inout) :: poollt
type(fract_type), intent(in) :: fract
logical, intent(in) :: continuation_restart

!...local variables
integer(i4) :: p

!----------------------------------------
!...reset pools if spinning up
if ((spinup) .and. (spinup_default) .and. &
    (.not. spinup_continue)) then
    pooldt%poollu(:) = dzero
    poollt%poolpft(:) = poolcont%poolpft_min
    poollt%poolpftp(:) = poolcont%poolpft_min
endif

!...PFT-dependent vertical-fraction metadata are set by
!...rebuild_pft_derived_state(), which is safe to call at runtime without
!...changing carbon mass or prognostic physical state.


!...set pool vertical distribution and
!...equilibrium information
do p=1,npoollu/3 !1,6 dead pools
   if (.not. continuation_restart) then
      pooldt%poollup(p) = pooldt%poollu(p)
      if (sum(pooldt%poollu_lay(p,:)) .ne. pooldt%poollu(p)) then
          pooldt%poollu_lay(p,:) = pooldt%poollu(p) &
              * pooldt%poollu_flay(p,:)
      endif

      equibdt%poollu_init(p) = pooldt%poollu(p)
      equibdt%poollu_min(p) = pooldt%poollu(p)
      equibdt%poollu_max(p) = pooldt%poollu(p)
   endif
enddo

!.. same as above but for C13
do p=npoollu/3+1,2*npoollu/3 !7,12 C13 dead pools
   if (.not. continuation_restart) then
      pooldt%poollup(p) = pooldt%poollu(p)
      if (sum(pooldt%poollu_lay(p,:)) .ne. pooldt%poollu(p)) then
          pooldt%poollu_lay(p,:) = pooldt%poollu(p) &
              * pooldt%poollu_flay(p,:)
      endif

      equibdt%poollu_init(p) = pooldt%poollu(p)
      equibdt%poollu_min(p) = pooldt%poollu(p)
      equibdt%poollu_max(p) = pooldt%poollu(p)
   endif
enddo

!.. same as above but for C14
do p=2*npoollu/3+1,npoollu !13,18 C14 dead pools
   if (.not. continuation_restart) then
      pooldt%poollup(p) = pooldt%poollu(p)
      if (sum(pooldt%poollu_lay(p,:)) .ne. pooldt%poollu(p)) then
          pooldt%poollu_lay(p,:) = pooldt%poollu(p) &
              * pooldt%poollu_flay(p,:)
      endif

      equibdt%poollu_init(p) = pooldt%poollu(p)
      equibdt%poollu_min(p) = pooldt%poollu(p)
      equibdt%poollu_max(p) = pooldt%poollu(p)
   endif
enddo


!...total C
do p=1,npoolpft/3 !1,5 live pools
   if (.not. continuation_restart) then
      poollt%poolpftp(p) = poollt%poolpft(p)
      if (sum(poollt%poolpft_lay(p,:)) .ne. poollt%poolpft(p)) then
          poollt%poolpft_lay(p,:) = poollt%poolpft(p) &
              * poollt%poolpft_flay(p,:)
      endif

      equiblt%poolpft_init(p) = poollt%poolpft(p)
      equiblt%poolpft_min(p) = poollt%poolpft(p)
      equiblt%poolpft_max(p) = poollt%poolpft(p)
   endif
enddo


!... same as above but for C13
do p=npoolpft/3+1,2*npoolpft/3 !6,10 live pools
   if (.not. continuation_restart) then
      poollt%poolpftp(p) = poollt%poolpft(p)
      if (sum(poollt%poolpft_lay(p,:)) .ne. poollt%poolpft(p)) then
          poollt%poolpft_lay(p,:) = poollt%poolpft(p) &
              * poollt%poolpft_flay(p,:)
      endif

      equiblt%poolpft_init(p) = poollt%poolpft(p)
      equiblt%poolpft_min(p) = poollt%poolpft(p)
      equiblt%poolpft_max(p) = poollt%poolpft(p)
   endif
enddo


!... same as above but for C14
do p=2*npoolpft/3+1,npoolpft !11,15 live pools
   if (.not. continuation_restart) then
      poollt%poolpftp(p) = poollt%poolpft(p)
      if (sum(poollt%poolpft_lay(p,:)) .ne. poollt%poolpft(p)) then
          poollt%poolpft_lay(p,:) = poollt%poolpft(p) &
              * poollt%poolpft_flay(p,:)
      endif

      equiblt%poolpft_init(p) = poollt%poolpft(p)
      equiblt%poolpft_min(p) = poollt%poolpft(p)
      equiblt%poolpft_max(p) = poollt%poolpft(p)
   endif
enddo


end subroutine setup_poolt


!-------------------------------------
!Routine to set the phenology variables
subroutine setup_phent( &
      gnum, lnum, pnum, pref, &
      phencont, poolcont, &
      daylen, daylendt, daylenmax, &
      lai, tm, vmax, &
      pooldt, poollt, phent, physcont, continuation_restart)

use kinds
use module_param, only: &
   phen_param, pool_param, phys_param
use module_pftinfo, only: &
   pft_pmeth, pmeth_gdd, pmeth_stg, pmeth_nvg
use module_pparams, only: tref
use module_sib, only: &
   poold_type, pooll_type, phen_type

implicit none

!...input variables
integer(i4), intent(in) :: gnum, lnum, pnum, pref
type(phen_param), intent(in) :: phencont
type(pool_param), intent(in) :: poolcont
type(phys_param), intent(in) :: physcont
logical, intent(in) :: continuation_restart
real(r4), intent(in) :: daylenmax
real(r8), intent(in) :: daylen, daylendt
real(r8), intent(in) :: lai, tm

real(r8), intent(inout) :: vmax
type(poold_type), intent(inout) :: pooldt
type(pooll_type), intent(inout) :: poollt
type(phen_type), intent(inout) :: phent

!...local variables
integer(i4) :: ipft, ips
real(r8) :: tmdf

!...setup local variables
ipft = pref
tmdf = ((tm - tref)*1.8) + 32.0
! phenave_wac is a derived combination of the restored weather and
! assimilation-memory terms.  Use the same definition as phen_update.F90.
phent%phenave_wac = 0.5 &
    * (phent%phenave_wa + phent%phenave_env)

if (.not. continuation_restart) then
   ! A new requib/default-IC simulation performs the historical phenology
   ! setup.  A true restart must not execute an extra phenology update before
   ! the first actual model timestep.
   call phen_gss(phencont, daylenmax, daylen, daylendt, phent)

   !...setup phenology stage index and phenology stage
   if (phent%phen_istage .le. 0) then
      if (pft_pmeth(pnum) .eq. pmeth_stg) then
         !print('(a)'),'   Setting Up Dynamic Phenology'
         call phen_dynamic( phencont, &
              daylenmax, daylen, daylendt, lai, &
              phent)
      elseif (pft_pmeth(pnum) .eq. pmeth_gdd) then
         !print('(a)'),'   Setting Up Defined Phenology'
         call phen_defined( &
              gnum, lnum, ipft, pnum, &
              phencont, poolcont, &
              lai, tmdf, phent, pooldt, poollt, physcont)
      elseif (pft_pmeth(pnum) .ne. pmeth_nvg) then
             print('(a)'),'Unknown phenology method.'
             print('(a,i3)'),'Method Selected: ',pft_pmeth(pnum)
             print('(a,i3)'),'Dynamic Phenology Method == ',pmeth_stg
             print('(a,i3)'),'Defined Phenology Method == ',pmeth_gdd
             print('(a)'),'Stopping.'
             stop
      endif
   endif
else
   if (phent%phen_istage .le. 0) then
      print*,'ERROR: continuation restart has invalid phen_istage at g/l: ', &
             gnum, lnum
      stop
   endif
endif

!...setup associated variables
ips = phent%phen_istage

!.....Photosynthetic Allocation
poollt%alloc_phen(:) = phencont%allocp(:,ips)
poollt%alloc(:) = poollt%alloc_phen(:)

!.....Dynamical Allocation
!.....not allowed during leaf-out or senescence
poollt%aadj_moist = phencont%adj_moist
poollt%aadj_temp = phencont%adj_temp
IF ((ips .le. 1) .or. (ips .eq. phencont%npstg)) THEN
    poollt%aadj_moist = .false.
    poollt%aadj_temp = .false.
ENDIF

!.....Canopy pool transfer rate
poollt%tfl_pstage = phencont%lptransfer(ips)

!.....Vmax
vmax = phencont%vmax0(ips)


end subroutine setup_phent
