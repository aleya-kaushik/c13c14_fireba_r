!--------------------------------------------------------------
subroutine firegfed_interp(indx, lon, lat, &
           sibg)

!--------------------------------------------------------------
!
! This subroutine interpolates the sibdrv fire emissions
!     between their read times
!
!--------------------------------------------------------------

use kinds
use module_io, only: &
   fire_step, fire_seccur, fire_secnext
use module_oparams, only: &
    fire_leaff, fire_stwdf, &
    fire_cdbf, fire_metlf, fire_strlf
use module_param, only: &
    poolcon, firecon
use module_pparams, only: &
     mwc, mol_to_umol, &
     month_names, pdb, &
     drytoc
!use module_pftinfo, only: pft_num
use module_pftinfo
use module_poolinfo
use module_sib, only: &
    gridcell_type
use module_sibconst, only: &
   npoolpft, npoollu, nsoil, ntpool, &
   fireb_print, fireb_stop, fireb_thresh
use module_time, only: &
   month, day, year, &
   dtisib, dtsib, sec_tot, wt_daily

implicit none

!...input variables
integer(i4), intent(in) :: indx
real(r4), intent(in) :: lon, lat
type (gridcell_type), intent(inout) :: sibg
!type (sscol_type), intent(in) :: sscolt
!type (soil_type), intent(in) :: soilt

!...parameters
integer(i4), parameter :: isave=5
real(r8) :: dnzero=1.E-10

!...interpolation variables
real(r8) :: facsibdrv  ! scaling factor between driver data points
real(r8) :: totemis, curemis, pcemis
real(r8) :: curemisc13, curemisc14, curemis_dm, curemisc13_dm, curemisc14_dm

!...distribution variables
integer(i4) :: ntpft
integer(i4), dimension(:), allocatable :: tpref, tpnum
real(r4), dimension(:), allocatable :: tparea
!real(r8), dimension(:), allocatable :: tpagb, tpagbtemp
real(r8), dimension(:), allocatable :: tpbiomass
character(len=clen), dimension(:), allocatable :: tpname

real(r8) :: tfarea
integer(i4), dimension(:), allocatable :: tsortref
!real(r8), dimension(:), allocatable :: flossb
real(r8), dimension(:,:), allocatable :: flosspft, flosslu
real(r8), dimension(:,:), allocatable :: flosspft_dm, flosslu_dm
real(r8), dimension(:,:), allocatable :: fgainlu
!real(r8), dimension(:), allocatable :: flossbc13

!...net balance variables
logical :: fb_err

!...misc variables
integer :: ido, idoc, l, p, s, myl, mylp
integer :: loc_c3g, loc_c4g
real(r8) :: myemis, myba, myfl, mymort, mytrans, tempemis, tempemisc13, tempemisc14
real(r8) :: mort_fract, c_compl
real(r8) :: tmppooltot, tmppoolc13, tempemispoolc13, rcemispoolc13
real(r8) :: tmppoolc14, tempemispoolc14, rcemispoolc14
real(r8) :: st_scalar, fl_scalar
real(r8) :: sm_sib, sm_scalar, sib_sm_sf
real(r8) :: sib_vpd_sf, sib_sf, vpd_scale_strength

integer(i4), dimension(1) :: tempstore
integer(i4) :: lp,frp,crp,wp,pp
integer(i4) :: lpc13,frpc13,crpc13,wpc13,ppc13
integer(i4) :: lpc14,frpc14,crpc14,wpc14,ppc14
integer(i4) :: cdbp, metlp, strlp, slitp, slowp, armp
integer(i4) :: cdbpc13, metlpc13, strlpc13, slitpc13, slowpc13, armpc13
integer(i4) :: cdbpc14, metlpc14, strlpc14, slitpc14, slowpc14, armpc14
integer(i4) :: j,jref,k,kref,n,nref,m,mref
integer(i4) :: flp, ffrp, fcrp, fwp, fpp
integer(i4) :: fcdbp, fmetlp, fstrlp, fslitp, fslowp, farmp
!real(r8) :: rcpooltest

!... scalars, more detail to be added later
!st_scalar = 0.5 ! moisture stress scaler for mortality and combustion completeness
!fl_scalar = 0 ! forest loss scaler for combustion completeness

!-------Calculatae quasi-equlibrium pools------------
lp =  pool_indx_leaf
frp = pool_indx_froot
crp = pool_indx_croot
wp =  pool_indx_stwd
pp =  pool_indx_prod

lpc13 =  pool_indx_leaf_c13-npoollu/3
frpc13 = pool_indx_froot_c13-npoollu/3
crpc13 = pool_indx_croot_c13-npoollu/3
wpc13 =  pool_indx_stwd_c13-npoollu/3
ppc13 =  pool_indx_prod_c13-npoollu/3

lpc14 =  pool_indx_leaf_c14-(2*(npoollu/3))
frpc14 = pool_indx_froot_c14-(2*(npoollu/3))
crpc14 = pool_indx_croot_c14-(2*(npoollu/3))
wpc14 =  pool_indx_stwd_c14-(2*(npoollu/3))
ppc14 =  pool_indx_prod_c14-(2*(npoollu/3))

cdbp  = pool_indx_cdb-npoolpft/3
metlp = pool_indx_metl-npoolpft/3
strlp = pool_indx_strl-npoolpft/3
slitp = pool_indx_slit-npoolpft/3
slowp = pool_indx_slow-npoolpft/3
armp  = pool_indx_arm-npoolpft/3

cdbpc13  = pool_indx_cdb_c13 - (2*(npoolpft/3))
metlpc13  = pool_indx_metl_c13 - (2*(npoolpft/3))
strlpc13  = pool_indx_strl_c13 - (2*(npoolpft/3))
slitpc13  = pool_indx_slit_c13 - (2*(npoolpft/3))
slowpc13 = pool_indx_slow_c13 - (2*(npoolpft/3))
armpc13  = pool_indx_arm_c13 - (2*(npoolpft/3))

cdbpc14  = pool_indx_cdb_c14 - npoolpft
metlpc14  = pool_indx_metl_c14 - npoolpft
strlpc14  = pool_indx_strl_c14 - npoolpft
slitpc14  = pool_indx_slit_c14 - npoolpft
slowpc14 = pool_indx_slow_c14 - npoolpft
armpc14  = pool_indx_arm_c14 - npoolpft

!indices for firecon variables
flp=1
ffrp=2
fcrp=3
fwp=4
fpp=5
fcdbp=6
fmetlp=7
fstrlp=8
fslitp=9
fslowp=10
farmp=11

!! print statements to check poolpft_loss and poolpft_lay
!print*,' '
!print*,'code: fire_interp'
!print*,'poolpft_dloss(1,1/2/3):',&
!    sibg%l(l)%poollt%poolpft_dloss(1,1),sibg%l(l)%poollt%poolpft_dloss(1,2),sibg%l(l)%poollt%poolpft_dloss(1,3)
!print*,'poolpft_dloss(2,1/2/3):',&
!    sibg%l(l)%poollt%poolpft_dloss(2,1),sibg%l(l)%poollt%poolpft_dloss(2,2),sibg%l(l)%poollt%poolpft_dloss(2,3)
!print*,'poolpft(1) :', sibg%l(l)%poollt%poolpft(1)
!print*,'poolpft(2) :', sibg%l(l)%poollt%poolpft(2)
!print*,'poolpft(3) :', sibg%l(l)%poollt%poolpft(3)
!print*,' '
 


!-----------------------------------------------
! reset values
do l=1, sibg%g_nlu
   sibg%l(l)%poollt%emis_fire = dzero
   sibg%l(l)%poollt%emis_firec13 = dzero
   sibg%l(l)%poollt%emis_firec14 = dzero
   sibg%l(l)%fract%d13cemisfire_pool = dzero
   !sibg%l(l)%poollt%loss_fire_lay(1:5,:) = dzero
   !sibg%l(l)%pooldt%loss_fire_lay(1:6,:) = dzero
   sibg%l(l)%poollt%loss_fire_lay(:,:) = dzero
   sibg%l(l)%pooldt%loss_fire_lay(:,:) = dzero
   sibg%l(l)%poollt%loss_transf_lay(:,:) = dzero
   sibg%l(l)%pooldt%loss_transf_lay(:,:) = dzero
   sibg%l(l)%pooldt%gain_transf_lay(:,:) = dzero
enddo
   
! only continue if fire emissions are being used
IF (fire_step .le. izero) RETURN

! get scaling factors
facsibdrv = dble(fire_seccur-sec_tot) / dble(fire_step)

! only continue if fire emissions are valid at this time
IF ((facsibdrv .GT. 1) .OR. (facsibdrv .LT. 0)) THEN
   sibg%gprogt%ba_dbg = dzero
   sibg%gprogt%ba_en2 = dzero
   sibg%gprogt%ba_en3 = dzero
   sibg%gprogt%ba_dnf = dzero
   sibg%gprogt%ba_eb1 = dzero
   sibg%gprogt%ba_eb2 = dzero
   sibg%gprogt%ba_db1 = dzero
   sibg%gprogt%ba_db2 = dzero
   sibg%gprogt%ba_db3 = dzero
   sibg%gprogt%ba_shb = dzero
   sibg%gprogt%ba_sha = dzero
   sibg%gprogt%ba_c3a = dzero
   sibg%gprogt%ba_c3g = dzero
   sibg%gprogt%ba_c4g = dzero
   sibg%gprogt%ba_c3c = dzero
   sibg%gprogt%ba_c4c = dzero
   sibg%gprogt%ba_mze = dzero
   sibg%gprogt%ba_soy = dzero
   sibg%gprogt%ba_wwt = dzero
   RETURN
ENDIF

! interpolate BA
sibg%gprogt%ba_dbg = facsibdrv*sibg%gprogt%ba_dbg1 &
         + (1.-facsibdrv) * sibg%gprogt%ba_dbg2
sibg%gprogt%ba_en2 = facsibdrv*sibg%gprogt%ba_en21 &
         + (1.-facsibdrv) * sibg%gprogt%ba_en22
sibg%gprogt%ba_en3 = facsibdrv*sibg%gprogt%ba_en31 &
         + (1.-facsibdrv) * sibg%gprogt%ba_en32
sibg%gprogt%ba_dnf = facsibdrv*sibg%gprogt%ba_dnf1 &
         + (1.-facsibdrv) * sibg%gprogt%ba_dnf2
sibg%gprogt%ba_eb1 = facsibdrv*sibg%gprogt%ba_eb11 &
         + (1.-facsibdrv) * sibg%gprogt%ba_eb12
sibg%gprogt%ba_eb2 = facsibdrv*sibg%gprogt%ba_eb21 &
         + (1.-facsibdrv) * sibg%gprogt%ba_eb22
sibg%gprogt%ba_db1 = facsibdrv*sibg%gprogt%ba_db11 &
         + (1.-facsibdrv) * sibg%gprogt%ba_db12
sibg%gprogt%ba_db2 = facsibdrv*sibg%gprogt%ba_db21 &
         + (1.-facsibdrv) * sibg%gprogt%ba_db22
sibg%gprogt%ba_db3 = facsibdrv*sibg%gprogt%ba_db31 &
         + (1.-facsibdrv) * sibg%gprogt%ba_db32
sibg%gprogt%ba_shb = facsibdrv*sibg%gprogt%ba_shb1 &
         + (1.-facsibdrv) * sibg%gprogt%ba_shb2
sibg%gprogt%ba_sha = facsibdrv*sibg%gprogt%ba_sha1 &
         + (1.-facsibdrv) * sibg%gprogt%ba_sha2
sibg%gprogt%ba_c3a = facsibdrv*sibg%gprogt%ba_c3a1 &
         + (1.-facsibdrv) * sibg%gprogt%ba_c3a2
sibg%gprogt%ba_c3g = facsibdrv*sibg%gprogt%ba_c3g1 &
         + (1.-facsibdrv) * sibg%gprogt%ba_c3g2
sibg%gprogt%ba_c4g = facsibdrv*sibg%gprogt%ba_c4g1 &
         + (1.-facsibdrv) * sibg%gprogt%ba_c4g2
sibg%gprogt%ba_c3c = facsibdrv*sibg%gprogt%ba_c3c1 &
         + (1.-facsibdrv) * sibg%gprogt%ba_c3c2
sibg%gprogt%ba_c4c = facsibdrv*sibg%gprogt%ba_c4c1 &
         + (1.-facsibdrv) * sibg%gprogt%ba_c4c2
sibg%gprogt%ba_mze = facsibdrv*sibg%gprogt%ba_mze1 &
         + (1.-facsibdrv) * sibg%gprogt%ba_mze2
sibg%gprogt%ba_soy = facsibdrv*sibg%gprogt%ba_soy1 &
         + (1.-facsibdrv) * sibg%gprogt%ba_soy2
sibg%gprogt%ba_wwt = facsibdrv*sibg%gprogt%ba_wwt1 &
         + (1.-facsibdrv) * sibg%gprogt%ba_wwt2

! new SiB-GFED fire calculation
! Emis = pool     * (MF*CC) * BA  * EF
! g    = kg DM m-2 *    %   *   * g(kgDM)^-1
! note: BA units are BA % PFT⁻¹ s⁻¹ (converted in fire_read_ code)

! current code calculation:
! mol C m-2 = mol C m-2  *   fraction  * fraction
! BA units: BA % PFT⁻¹ s⁻¹: converted from /3h to /s in fire_read_ code
 
   !totemis = sibg%gprogt%firec * dtsib
   ntpft = sibg%g_nlu
   ! number of land units/PFTs per cell -> =1 for Hyy when run as 1.0ENF
   allocate(tpref(ntpft),tpnum(ntpft),tparea(ntpft))
   !allocate(tpagb(ntpft),tpbiomass(ntpft))
   allocate(tpname(ntpft))
   tpref(:) = sibg%l(1:ntpft)%ipft
   tpnum(:) = pft_num(tpref)
   do l=1,ntpft
     tpname(l) = pft_name(tpnum(l))
!     print*,'l, tpname: ',l, tpname(l)
   enddo
   tparea(:) = sibg%l(1:ntpft)%larea

   !...remove carbon from all PFTs
   allocate(flosspft(ntpft,npoolpft))
   allocate(flosslu(ntpft,npoollu))
   allocate(fgainlu(ntpft,npoollu))
   allocate(flosspft_dm(ntpft,npoolpft))
   allocate(flosslu_dm(ntpft,npoollu))
   flosspft(:,:) = dzero
   flosslu(:,:) = dzero
   fgainlu(:,:) = dzero
   flosspft_dm(:,:) = dzero
   flosslu_dm(:,:) = dzero

!print*,'tpref: ',tpref
!print*,'tpnum: ',tpnum
!print*,'pft_name: ',pft_name(tpnum)

   do l=1,ntpft

   curemis = dzero
   curemis_dm = dzero
   curemisc13 = dzero
   curemisc13_dm = dzero
   curemisc14 = dzero
   curemisc14_dm = dzero

      myl=l
      mylp=tpnum(l)
      !... Iniital pool totals come from c13_iso_calc 
      !... for sum of (lp+wp+cbdp+metlp+strlp)
      !tmppooltot = sibg%l(myl)%fract%poolemistotC
      !tmppoolc13 = sibg%l(myl)%fract%poolemisc13

      !... read in BA based on pftname
      !.. note: multiply BA s-1 by dtsib to remove seconds unit

      !.. below is for v14, where global BA totals included pft_area
      !if (tpname(l) .eq. "dbg") myba=sibg%gprogt%ba_dbg*dtsib*tparea(l)
      !if (tpname(l) .eq. "en2") myba=sibg%gprogt%ba_en2*dtsib*tparea(l)
      !if (tpname(l) .eq. "en3") myba=sibg%gprogt%ba_en3*dtsib*tparea(l)     
      !if (tpname(l) .eq. "dnf") myba=sibg%gprogt%ba_dnf*dtsib*tparea(l)
      !if (tpname(l) .eq. "eb1") myba=sibg%gprogt%ba_eb1*dtsib*tparea(l)  
      !if (tpname(l) .eq. "eb2") myba=sibg%gprogt%ba_eb2*dtsib*tparea(l)  
      !if (tpname(l) .eq. "db1") myba=sibg%gprogt%ba_db1*dtsib*tparea(l)  
      !if (tpname(l) .eq. "db2") myba=sibg%gprogt%ba_db2*dtsib*tparea(l)  
      !if (tpname(l) .eq. "db3") myba=sibg%gprogt%ba_db3*dtsib*tparea(l)  
      !if (tpname(l) .eq. "shb") myba=sibg%gprogt%ba_shb*dtsib*tparea(l)  
      !if (tpname(l) .eq. "sha") myba=sibg%gprogt%ba_sha*dtsib*tparea(l)  
      !if (tpname(l) .eq. "c3a") myba=sibg%gprogt%ba_c3a*dtsib*tparea(l)  
      !if (tpname(l) .eq. "c3g") myba=sibg%gprogt%ba_c3g*dtsib*tparea(l)  
      !if (tpname(l) .eq. "c4g") myba=sibg%gprogt%ba_c4g*dtsib*tparea(l)  
      !if (tpname(l) .eq. "c3c") myba=sibg%gprogt%ba_c3c*dtsib*tparea(l)  
      !if (tpname(l) .eq. "c4c") myba=sibg%gprogt%ba_c4c*dtsib*tparea(l)  
      !if (tpname(l) .eq. "mze") myba=sibg%gprogt%ba_mze*dtsib*tparea(l)  
      !if (tpname(l) .eq. "soy") myba=sibg%gprogt%ba_soy*dtsib*tparea(l)  
      !if (tpname(l) .eq. "wwt") myba=sibg%gprogt%ba_wwt*dtsib*tparea(l)  
      
      !v15.. without pft_area
      if (tpname(l) .eq. "dbg") myba=sibg%gprogt%ba_dbg*dtsib
      if (tpname(l) .eq. "en2") myba=sibg%gprogt%ba_en2*dtsib
      if (tpname(l) .eq. "en3") myba=sibg%gprogt%ba_en3*dtsib     
      if (tpname(l) .eq. "dnf") myba=sibg%gprogt%ba_dnf*dtsib
      if (tpname(l) .eq. "eb1") myba=sibg%gprogt%ba_eb1*dtsib  
      if (tpname(l) .eq. "eb2") myba=sibg%gprogt%ba_eb2*dtsib  
      if (tpname(l) .eq. "db1") myba=sibg%gprogt%ba_db1*dtsib  
      if (tpname(l) .eq. "db2") myba=sibg%gprogt%ba_db2*dtsib  
      if (tpname(l) .eq. "db3") myba=sibg%gprogt%ba_db3*dtsib  
      if (tpname(l) .eq. "shb") myba=sibg%gprogt%ba_shb*dtsib  
      if (tpname(l) .eq. "sha") myba=sibg%gprogt%ba_sha*dtsib  
      if (tpname(l) .eq. "c3a") myba=sibg%gprogt%ba_c3a*dtsib  
      if (tpname(l) .eq. "c3g") myba=sibg%gprogt%ba_c3g*dtsib  
      if (tpname(l) .eq. "c4g") myba=sibg%gprogt%ba_c4g*dtsib  
      if (tpname(l) .eq. "c3c") myba=sibg%gprogt%ba_c3c*dtsib  
      if (tpname(l) .eq. "c4c") myba=sibg%gprogt%ba_c4c*dtsib  
      if (tpname(l) .eq. "mze") myba=sibg%gprogt%ba_mze*dtsib  
      if (tpname(l) .eq. "soy") myba=sibg%gprogt%ba_soy*dtsib  
      if (tpname(l) .eq. "wwt") myba=sibg%gprogt%ba_wwt*dtsib                                                           
   !print*,'myl: ',myl
   !print*,'mylp: ',mylp
   !print*,'pftname(l): ',tpname(l)
   !print*,'pftname(l), MFL:',tpname(l),firecon(mylp)%mfl(flp),firecon(mylp)%mfl(frp),firecon(mylp)%mfl(crp),firecon(mylp)%mfl(wp),firecon(mylp)%mfl(pp)
   !print*,'pftname(l), MFH:',tpname(l),firecon(mylp)%mfh(flp),firecon(mylp)%mfh(frp),firecon(mylp)%mfh(crp),firecon(mylp)%mfh(wp),firecon(mylp)%mfh(pp)

!print*,'indexes :',lp,frp,crp,wp,pp,cdbp,metlp,strlp,slitp,slowp,armp

   !print*,'pftname(l), CCL:',tpname(l),firecon(mylp)%ccl(flp),firecon(mylp)%ccl(frp),firecon(mylp)%ccl(crp),firecon(mylp)%ccl(wp),firecon(mylp)%ccl(pp),firecon(mylp)%ccl(cdbp),firecon(mylp)%ccl(metlp),firecon(mylp)%ccl(strlp),firecon(mylp)%ccl(slitp),firecon(mylp)%ccl(slowp),firecon(mylp)%ccl(armp)
   !print*,'pftname(l), CCH:',tpname(l),firecon(mylp)%cch(flp),firecon(mylp)%cch(frp),firecon(mylp)%cch(crp),firecon(mylp)%cch(wp),firecon(mylp)%cch(pp),firecon(mylp)%cch(cdbp),firecon(mylp)%cch(metlp),firecon(mylp)%cch(strlp),firecon(mylp)%cch(slitp),firecon(mylp)%cch(slowp),firecon(mylp)%cch(armp)


    ! Soil moisture scalar
    sm_sib    = sibg%l(l)%sscolt%www_ice(1) + sibg%l(l)%sscolt%www_liq(1)
    sm_scalar = 1.0 - (sm_sib / sibg%l(l)%soilt%poros)
    sm_scalar = MAX(0.1, MIN(1.0, sm_scalar))

    ! SiB soil moisture stress factor
    sib_sm_sf = sm_scalar
    st_scalar = sib_sm_sf

    ! next step.. 
    ! sib_vpd_sf...
    !...
    ! sib_sf = sib_sm_sf + vpd_scale_strength*sib_vpd_sf
    ! sib_sf = MAX(0.0,MIN(1.0,sib_sf))
    ! st_scalar = sib_sf

      !-----Conservative burned-area fire removal-----
      ! Apply the BA/mortality/combustion formulation to carbon that is
      ! actually still available after gains and losses accumulated in the
      ! current update period.  Multi-layer pools retain the existing layer
      ! fractions, but each layer is capped independently so fire cannot
      ! overdraw a depleted layer.  C13 and C14 losses are derived from the
      ! actual capped total-C removal in each layer.
      CALL burn_live_single(lp, lpc13, lpc14, flp)
      CALL burn_live_single(wp, wpc13, wpc14, fwp)
      CALL burn_dead_single(metlp, metlpc13, metlpc14, fmetlp)
      CALL burn_dead_single(strlp, strlpc13, strlpc14, fstrlp)
      CALL burn_dead_single(cdbp, cdbpc13, cdbpc14, fcdbp)
      CALL burn_live_single(pp, ppc13, ppc14, fpp)
      CALL burn_live_multi(frp, frpc13, frpc14, ffrp)
      CALL burn_live_multi(crp, crpc13, crpc14, fcrp)
      CALL burn_dead_multi(slitp, slitpc13, slitpc14, fslitp)
      CALL burn_dead_multi(slowp, slowpc13, slowpc14, fslowp)
      CALL burn_dead_multi(armp, armpc13, armpc14, farmp)

      !-----Carbon Transfers from live pools-----
      !...C pools
      do n=1,npoolpft/3
         do m=npoolpft/3+1,ntpool/3
            mref=m-npoolpft/3
            if (poolcon(myl)%pool_trans_frac(n,m) > dzero) then
               if ((pool_indx_lay(n) .eq. 1) .and. (pool_indx_lay(m) .eq. 1)) then
                  sibg%l(myl)%pooldt%gain_transf_lay(mref,1) = &
                       sibg%l(myl)%pooldt%gain_transf_lay(mref,1) + &
                       sibg%l(myl)%poollt%loss_transf_lay(n,1) * &
                       poolcon(myl)%pool_trans_frac(n,m)
               elseif ((pool_indx_lay(n) .eq. 1) .and. (pool_indx_lay(m) .eq. nsoil)) then
                  do s=1,nsoil
                     sibg%l(myl)%pooldt%gain_transf_lay(mref,s) = &
                          sibg%l(myl)%pooldt%gain_transf_lay(mref,s) + &
                          sibg%l(myl)%poollt%loss_transf_lay(n,1) * &
                          sibg%l(myl)%vegt%rootf(s) * poolcon(myl)%pool_trans_frac(n,m)
                  enddo
               elseif ((pool_indx_lay(n) .eq. nsoil) .and. (pool_indx_lay(m) .eq. nsoil)) then
                  do s=1,nsoil
                     sibg%l(myl)%pooldt%gain_transf_lay(mref,s) = &
                          sibg%l(myl)%pooldt%gain_transf_lay(mref,s) + &
                          sibg%l(myl)%poollt%loss_transf_lay(n,s) * &
                          poolcon(myl)%pool_trans_frac(n,m)
                  enddo
               else
                  print*,'Mismatching levels between pool transfers.'
                  print*,'Stopping in firegfed_interp live transf.'
                  stop
               endif
            endif
         enddo
      enddo

      !...C13 pools
      do n=npoolpft/3+1,2*npoolpft/3
         nref=n+npoollu/3
         do m=ntpool/3+npoolpft/3+1,2*ntpool/3
            mref=m-2*npoolpft/3
            if (poolcon(myl)%pool_trans_frac(nref,m) > dzero) then
               if ((pool_indx_lay(nref) .eq. 1) .and. (pool_indx_lay(m) .eq. 1)) then
                  sibg%l(myl)%pooldt%gain_transf_lay(mref,1) = &
                       sibg%l(myl)%pooldt%gain_transf_lay(mref,1) + &
                       sibg%l(myl)%poollt%loss_transf_lay(n,1) * &
                       poolcon(myl)%pool_trans_frac(nref,m)
               elseif ((pool_indx_lay(nref) .eq. 1) .and. (pool_indx_lay(m) .eq. nsoil)) then
                  do s=1,nsoil
                     sibg%l(myl)%pooldt%gain_transf_lay(mref,s) = &
                          sibg%l(myl)%pooldt%gain_transf_lay(mref,s) + &
                          sibg%l(myl)%poollt%loss_transf_lay(n,1) * &
                          sibg%l(myl)%vegt%rootf(s) * poolcon(myl)%pool_trans_frac(nref,m)
                  enddo
               elseif ((pool_indx_lay(nref) .eq. nsoil) .and. (pool_indx_lay(m) .eq. nsoil)) then
                  do s=1,nsoil
                     sibg%l(myl)%pooldt%gain_transf_lay(mref,s) = &
                          sibg%l(myl)%pooldt%gain_transf_lay(mref,s) + &
                          sibg%l(myl)%poollt%loss_transf_lay(n,s) * &
                          poolcon(myl)%pool_trans_frac(nref,m)
                  enddo
               else
                  print*,'Mismatching levels between pool transfers.'
                  print*,'Stopping in firegfed_interp live transf c13.'
                  stop
               endif
            endif
         enddo
      enddo

      !...C14 pools
      do n=2*npoolpft/3+1,npoolpft
         nref=n+2*npoollu/3
         do m=2*ntpool/3+npoolpft/3+1,ntpool
            mref=m-npoolpft
            if (poolcon(myl)%pool_trans_frac(nref,m) > dzero) then
               if ((pool_indx_lay(nref) .eq. 1) .and. (pool_indx_lay(m) .eq. 1)) then
                  sibg%l(myl)%pooldt%gain_transf_lay(mref,1) = &
                       sibg%l(myl)%pooldt%gain_transf_lay(mref,1) + &
                       sibg%l(myl)%poollt%loss_transf_lay(n,1) * &
                       poolcon(myl)%pool_trans_frac(nref,m)
               elseif ((pool_indx_lay(nref) .eq. 1) .and. (pool_indx_lay(m) .eq. nsoil)) then
                  do s=1,nsoil
                     sibg%l(myl)%pooldt%gain_transf_lay(mref,s) = &
                          sibg%l(myl)%pooldt%gain_transf_lay(mref,s) + &
                          sibg%l(myl)%poollt%loss_transf_lay(n,1) * &
                          sibg%l(myl)%vegt%rootf(s) * poolcon(myl)%pool_trans_frac(nref,m)
                  enddo
               elseif ((pool_indx_lay(nref) .eq. nsoil) .and. (pool_indx_lay(m) .eq. nsoil)) then
                  do s=1,nsoil
                     sibg%l(myl)%pooldt%gain_transf_lay(mref,s) = &
                          sibg%l(myl)%pooldt%gain_transf_lay(mref,s) + &
                          sibg%l(myl)%poollt%loss_transf_lay(n,s) * &
                          poolcon(myl)%pool_trans_frac(nref,m)
                  enddo
               else
                  print*,'Mismatching levels between pool transfers.'
                  print*,'Stopping in firegfed_interp live transf c14.'
                  stop
               endif
            endif
         enddo
      enddo

      !-----Carbon Transfers from dead pools-----
      !...C pools
      do j=1,npoollu/3
         jref=j+npoolpft/3
         do k=1,npoollu/3
            kref=k+npoolpft/3
            if (poolcon(myl)%pool_trans_frac(jref,kref) > dzero) then
               if ((pool_indx_lay(jref) .eq. 1) .and. (pool_indx_lay(kref) .eq. 1)) then
                  sibg%l(myl)%pooldt%gain_transf_lay(k,1) = &
                       sibg%l(myl)%pooldt%gain_transf_lay(k,1) + &
                       sibg%l(myl)%pooldt%loss_transf_lay(j,1) * &
                       poolcon(myl)%pool_trans_frac(jref,kref)
               elseif ((pool_indx_lay(jref) .eq. 1) .and. (pool_indx_lay(kref) .eq. nsoil)) then
                  do s=1,nsoil
                     sibg%l(myl)%pooldt%gain_transf_lay(k,s) = &
                          sibg%l(myl)%pooldt%gain_transf_lay(k,s) + &
                          sibg%l(myl)%pooldt%loss_transf_lay(j,1) * &
                          sibg%l(myl)%vegt%rootf(s) * poolcon(myl)%pool_trans_frac(jref,kref)
                  enddo
               elseif ((pool_indx_lay(jref) .eq. nsoil) .and. (pool_indx_lay(kref) .eq. nsoil)) then
                  sibg%l(myl)%pooldt%gain_transf_lay(k,:) = &
                       sibg%l(myl)%pooldt%gain_transf_lay(k,:) + &
                       sibg%l(myl)%pooldt%loss_transf_lay(j,:) * &
                       poolcon(myl)%pool_trans_frac(jref,kref)
               else
                  print*,'Mismatching levels between pool transfers.'
                  print*,'Stopping in firegfed_interp dead transf.'
                  stop
               endif
            endif
         enddo
      enddo

      !...C13 pools
      do j=npoollu/3+1,2*npoollu/3
         jref=j+2*npoolpft/3
         do k=npoollu/3+1,2*npoollu/3
            kref=k+2*npoolpft/3
            if (poolcon(myl)%pool_trans_frac(jref,kref) > dzero) then
               if ((pool_indx_lay(jref) .eq. 1) .and. (pool_indx_lay(kref) .eq. 1)) then
                  sibg%l(myl)%pooldt%gain_transf_lay(k,1) = &
                       sibg%l(myl)%pooldt%gain_transf_lay(k,1) + &
                       sibg%l(myl)%pooldt%loss_transf_lay(j,1) * &
                       poolcon(myl)%pool_trans_frac(jref,kref)
               elseif ((pool_indx_lay(jref) .eq. 1) .and. (pool_indx_lay(kref) .eq. nsoil)) then
                  do s=1,nsoil
                     sibg%l(myl)%pooldt%gain_transf_lay(k,s) = &
                          sibg%l(myl)%pooldt%gain_transf_lay(k,s) + &
                          sibg%l(myl)%pooldt%loss_transf_lay(j,1) * &
                          sibg%l(myl)%vegt%rootf(s) * poolcon(myl)%pool_trans_frac(jref,kref)
                  enddo
               elseif ((pool_indx_lay(jref) .eq. nsoil) .and. (pool_indx_lay(kref) .eq. nsoil)) then
                  sibg%l(myl)%pooldt%gain_transf_lay(k,:) = &
                       sibg%l(myl)%pooldt%gain_transf_lay(k,:) + &
                       sibg%l(myl)%pooldt%loss_transf_lay(j,:) * &
                       poolcon(myl)%pool_trans_frac(jref,kref)
               else
                  print*,'Mismatching levels between pool transfers.'
                  print*,'Stopping in firegfed_interp dead transf c13.'
                  stop
               endif
            endif
         enddo
      enddo

      !...C14 pools
      do j=2*npoollu/3+1,npoollu
         jref=j+npoolpft
         do k=2*npoollu/3+1,npoollu
            kref=k+npoolpft
            if (poolcon(myl)%pool_trans_frac(jref,kref) > dzero) then
               if ((pool_indx_lay(jref) .eq. 1) .and. (pool_indx_lay(kref) .eq. 1)) then
                  sibg%l(myl)%pooldt%gain_transf_lay(k,1) = &
                       sibg%l(myl)%pooldt%gain_transf_lay(k,1) + &
                       sibg%l(myl)%pooldt%loss_transf_lay(j,1) * &
                       poolcon(myl)%pool_trans_frac(jref,kref)
               elseif ((pool_indx_lay(jref) .eq. 1) .and. (pool_indx_lay(kref) .eq. nsoil)) then
                  do s=1,nsoil
                     sibg%l(myl)%pooldt%gain_transf_lay(k,s) = &
                          sibg%l(myl)%pooldt%gain_transf_lay(k,s) + &
                          sibg%l(myl)%pooldt%loss_transf_lay(j,1) * &
                          sibg%l(myl)%vegt%rootf(s) * poolcon(myl)%pool_trans_frac(jref,kref)
                  enddo
               elseif ((pool_indx_lay(jref) .eq. nsoil) .and. (pool_indx_lay(kref) .eq. nsoil)) then
                  sibg%l(myl)%pooldt%gain_transf_lay(k,:) = &
                       sibg%l(myl)%pooldt%gain_transf_lay(k,:) + &
                       sibg%l(myl)%pooldt%loss_transf_lay(j,:) * &
                       poolcon(myl)%pool_trans_frac(jref,kref)
               else
                  print*,'Mismatching levels between pool transfers.'
                  print*,'Stopping in firegfed_interp dead transf c14.'
                  stop
               endif
            endif
         enddo
      enddo

      !...Save actual fire emissions.  curemis/curemisc13/curemisc14 are
      !...the sums of the capped pool-level emissions calculated above.
      sibg%l(myl)%poollt%emis_fire = curemis*dtisib
      sibg%l(myl)%poollt%emis_firec13 = curemisc13*dtisib
      sibg%l(myl)%poollt%emis_firec14 = curemisc14*dtisib

      !...Pool-fire delta diagnostic from the carbon actually emitted.
      sibg%l(myl)%fract%d13cemisfire_pool = dzero
      if ((curemis .GT. dnzero) .AND. &
          (curemisc13 .GT. dzero) .AND. &
          (curemisc13 .LT. curemis)) then
         rcemispoolc13 = curemisc13/curemis
         sibg%l(myl)%fract%d13cemisfire_pool = &
          (((rcemispoolc13/(1.0D0-rcemispoolc13))/pdb)-1.0D0)*1000.0D0
      endif

     !...Accumulate gains for pool updates
     sibg%l(myl)%pooldt%poollu_dgain = sibg%l(myl)%pooldt%poollu_dgain &!(npoollu,nsoil)
        + sibg%l(myl)%pooldt%gain_transf_lay*dtsib

   enddo !loop thrugh PFTs

!section below is done for individual pools above
!pooldt%poollu_dloss = pooldt%poollu_dloss & !(npoollu,nsoil)
!   + pooldt%loss_resp_lay*dtsib + pooldt%loss_trans_lay*dtsib


      !...No separate isotope reconstruction is needed: emis_fire* already
      !...contains the actual C/C13/C14 emitted from the capped pool losses.

   !...Print Result s
   !IF ((fb_err) .O R. (fireb_print)) THEN
   !   print*,''    
   !   print*,'---F IRE CARBON---'
   !   IF (fb_err)  THEN
   !      print('(a )'),'!!Fire Carbon Imbalance!!'
   !      print*,'F ire Emissions Mismatch (mol C/m2): ', curemis
   !   ENDIF        
                    
   !   print('(a,a, i3,a,i4)'), &
   !       '      D ate: ', trim(month_names(month)), day, ', ', year
   !   print('(a,i6 ,2f8.2)'),  '      Point/Lon/Lat: ', indx, lon, lat
   !   print('(a,i14)'),       '      Current Model Second: ', sec_tot
   !   print('(a,2i12)'),      '      Current/Next Fire Second: ', &
   !           fire_seccur, fire_secnext

   !   print*,''
   !   print('(a,i4)'),     '                                ntpft: ', &
   !        ntpft
   !   print*,''
   !   print('(a,f18.8)'),     '                            tfarea: ', &
   !        tfarea

   !   print*,''
   !   print('(a,f18.8)'),     '             sumtotal biomass(m-2): ', &
   !        sum(tpbiomass)
!  !    print*,''
!  !    print('(a,f6.2,a,f6.2,a,f6.2,a,f6.2,a,f6.2)'),&
!  !         'biomass(m-2): ', &
!  !         tpbiomass(1),' ',tpbiomass(2),' ',tpbiomass(3),' ',&
!  !         tpbiomass(4),' ',tpbiomass(5)
   !   print*,''
   !   print('(a,f18.8)'),     'sumtotal above-ground biomass(m-2): ', &
   !        sum(tpagb)
!  !    print*,''
!  !    print('(a,f6.2,a,f6.2,a,f6.2,a,f6.2,a,f6.2)'),&
!  !         'above-ground biomass(m-2): ', &
!  !         tpagb(1),' ',tpagb(2),' ',tpagb(3),' ',&
!  !         tpagb(4),' ',tpagb(5)

   !   print*,''
   !   print('(a,f18.8)'),     '      Fire C Emissions (umol/m2/s): ', &
   !        sibg%gprogt%firec*mol_to_umol
   !   print('(a,f18.8)'),     '      Time-Step C Losses (mol/m2):', totemis

   !   tempemis = sum(flosspft(:,1:5)) + sum(flosslu(:,1:6))
   !   print('(a,f18.8)'),     '      SiB4 C Removal (mol/m2):    ', tempemis
   !   print('(a)'),           '         PFT    Loss          %-BioBurned   Fraction'
   !      do l=1,ido
   !         myl = tsortref(l)
   !         curemis = sum(flosslu(myl,1:6)) + sum(flosspft(myl,1:5))
   !         if (tpbiomass(myl) .gt. dzero) then
   !            pcemis = curemis/tpbiomass(myl)*100.
   !         else
   !            pcemis = dzero
   !         endif
   !         print('(a,i2,2f14.8,a,f6.2)'), '          ', tpref(myl),  &
   !            curemis, pcemis, '  ',tparea(myl)/tfarea
   !      enddo
   !   print('(a,f12.8)'),     '      Non-Matched C Respired: ', sum(flossb)

   !   IF (fb_err .AND. (fireb_stop)) STOP
   !ENDIF  !print

   deallocate(tpref,tpnum,tparea)
   !deallocate(tpagb,tpagbtemp)
   deallocate(flosspft,flosslu,fgainlu)
   deallocate(flosspft_dm,flosslu_dm)

contains

   subroutine burn_live_single(pidx, p13idx, p14idx, fidx)
      integer(i4), intent(in) :: pidx, p13idx, p14idx, fidx
      real(r8) :: avail, mortfrac, ccomp, mort, emis, trans
      real(r8) :: r13, r14, emis13, emis14

      avail = MAX(dzero, &
           sibg%l(myl)%poollt%poolpft_lay(pidx,1) &
           + sibg%l(myl)%poollt%poolpft_dgain(pidx,1) &
           - sibg%l(myl)%poollt%poolpft_dloss(pidx,1) &
           - poolcon(tpnum(myl))%poolpft_min(pidx))

      mortfrac = firecon(mylp)%mfl(fidx) + &
           st_scalar*(firecon(mylp)%mfh(fidx) - firecon(mylp)%mfl(fidx))
      ccomp = firecon(mylp)%ccl(fidx) + &
           st_scalar*(firecon(mylp)%cch(fidx) - firecon(mylp)%ccl(fidx))

      mort = MIN(avail, MAX(dzero, avail*myba*mortfrac))
      emis = mort*ccomp
      trans = mort*(1.0_r8-ccomp)

      r13 = sibg%l(myl)%poollt%rcpoolpft_lay(p13idx,1)
      r14 = sibg%l(myl)%poollt%rcpoolpft_lay(p14idx,1)
      emis13 = r13*emis
      emis14 = r14*emis

      sibg%l(myl)%poollt%loss_fire_lay(pidx,1) = emis*dtisib
      sibg%l(myl)%poollt%loss_transf_lay(pidx,1) = trans*dtisib
      sibg%l(myl)%poollt%poolpft_dloss(pidx,1) = &
           sibg%l(myl)%poollt%poolpft_dloss(pidx,1) + mort

      sibg%l(myl)%poollt%loss_fire_lay(p13idx,1) = emis13*dtisib
      sibg%l(myl)%poollt%loss_transf_lay(p13idx,1) = r13*trans*dtisib
      sibg%l(myl)%poollt%poolpft_dloss(p13idx,1) = &
           sibg%l(myl)%poollt%poolpft_dloss(p13idx,1) + r13*mort

      sibg%l(myl)%poollt%loss_fire_lay(p14idx,1) = emis14*dtisib
      sibg%l(myl)%poollt%loss_transf_lay(p14idx,1) = r14*trans*dtisib
      sibg%l(myl)%poollt%poolpft_dloss(p14idx,1) = &
           sibg%l(myl)%poollt%poolpft_dloss(p14idx,1) + r14*mort

      flosspft(myl,pidx) = emis
      flosspft(myl,p13idx) = emis13
      flosspft(myl,p14idx) = emis14
      flosspft_dm(myl,pidx) = emis*drytoc
      flosspft_dm(myl,p13idx) = emis13*drytoc
      flosspft_dm(myl,p14idx) = emis14*drytoc

      curemis = curemis + emis
      curemis_dm = curemis_dm + emis*drytoc
      curemisc13 = curemisc13 + emis13
      curemisc13_dm = curemisc13_dm + emis13*drytoc
      curemisc14 = curemisc14 + emis14
      curemisc14_dm = curemisc14_dm + emis14*drytoc
   end subroutine burn_live_single


   subroutine burn_live_multi(pidx, p13idx, p14idx, fidx)
      integer(i4), intent(in) :: pidx, p13idx, p14idx, fidx
      integer(i4) :: ss
      real(r8) :: avail_sum, mortfrac, ccomp, target_mort
      real(r8) :: layer_mort, layer_emis, layer_trans
      real(r8) :: emis, emis13, emis14, r13, r14
      real(r8), dimension(nsoil) :: avail_lay

      avail_lay(:) = dzero
      do ss=1,nsoil
         avail_lay(ss) = MAX(dzero, &
              sibg%l(myl)%poollt%poolpft_lay(pidx,ss) &
              + sibg%l(myl)%poollt%poolpft_dgain(pidx,ss) &
              - sibg%l(myl)%poollt%poolpft_dloss(pidx,ss) &
              - poolcon(tpnum(myl))%poolpft_min(pidx) &
                * sibg%l(myl)%poollt%poolpft_flay(pidx,ss))
      enddo
      avail_sum = sum(avail_lay)

      mortfrac = firecon(mylp)%mfl(fidx) + &
           st_scalar*(firecon(mylp)%mfh(fidx) - firecon(mylp)%mfl(fidx))
      ccomp = firecon(mylp)%ccl(fidx) + &
           st_scalar*(firecon(mylp)%cch(fidx) - firecon(mylp)%ccl(fidx))
      target_mort = MIN(avail_sum, MAX(dzero, avail_sum*myba*mortfrac))

      emis = dzero
      emis13 = dzero
      emis14 = dzero
      do ss=1,nsoil
         layer_mort = MIN(avail_lay(ss), MAX(dzero, &
              target_mort*sibg%l(myl)%poollt%poolpft_flay(pidx,ss)))
         layer_emis = layer_mort*ccomp
         layer_trans = layer_mort*(1.0_r8-ccomp)
         r13 = sibg%l(myl)%poollt%rcpoolpft_lay(p13idx,ss)
         r14 = sibg%l(myl)%poollt%rcpoolpft_lay(p14idx,ss)

         sibg%l(myl)%poollt%loss_fire_lay(pidx,ss) = layer_emis*dtisib
         sibg%l(myl)%poollt%loss_transf_lay(pidx,ss) = layer_trans*dtisib
         sibg%l(myl)%poollt%poolpft_dloss(pidx,ss) = &
              sibg%l(myl)%poollt%poolpft_dloss(pidx,ss) + layer_mort

         sibg%l(myl)%poollt%loss_fire_lay(p13idx,ss) = &
              r13*layer_emis*dtisib
         sibg%l(myl)%poollt%loss_transf_lay(p13idx,ss) = &
              r13*layer_trans*dtisib
         sibg%l(myl)%poollt%poolpft_dloss(p13idx,ss) = &
              sibg%l(myl)%poollt%poolpft_dloss(p13idx,ss) + r13*layer_mort

         sibg%l(myl)%poollt%loss_fire_lay(p14idx,ss) = &
              r14*layer_emis*dtisib
         sibg%l(myl)%poollt%loss_transf_lay(p14idx,ss) = &
              r14*layer_trans*dtisib
         sibg%l(myl)%poollt%poolpft_dloss(p14idx,ss) = &
              sibg%l(myl)%poollt%poolpft_dloss(p14idx,ss) + r14*layer_mort

         emis = emis + layer_emis
         emis13 = emis13 + r13*layer_emis
         emis14 = emis14 + r14*layer_emis
      enddo

      flosspft(myl,pidx) = emis
      flosspft(myl,p13idx) = emis13
      flosspft(myl,p14idx) = emis14
      flosspft_dm(myl,pidx) = emis*drytoc
      flosspft_dm(myl,p13idx) = emis13*drytoc
      flosspft_dm(myl,p14idx) = emis14*drytoc

      curemis = curemis + emis
      curemis_dm = curemis_dm + emis*drytoc
      curemisc13 = curemisc13 + emis13
      curemisc13_dm = curemisc13_dm + emis13*drytoc
      curemisc14 = curemisc14 + emis14
      curemisc14_dm = curemisc14_dm + emis14*drytoc
   end subroutine burn_live_multi


   subroutine burn_dead_single(pidx, p13idx, p14idx, fidx)
      integer(i4), intent(in) :: pidx, p13idx, p14idx, fidx
      real(r8) :: avail, ccomp, emis, r13, r14, emis13, emis14

      avail = MAX(dzero, &
           sibg%l(myl)%pooldt%poollu_lay(pidx,1) &
           + sibg%l(myl)%pooldt%poollu_dgain(pidx,1) &
           - sibg%l(myl)%pooldt%poollu_dloss(pidx,1))
      ccomp = firecon(mylp)%ccl(fidx) + &
           st_scalar*(firecon(mylp)%cch(fidx) - firecon(mylp)%ccl(fidx))
      emis = MIN(avail, MAX(dzero, avail*myba*ccomp))

      r13 = sibg%l(myl)%pooldt%rcpoollu_lay(p13idx,1)
      r14 = sibg%l(myl)%pooldt%rcpoollu_lay(p14idx,1)
      emis13 = r13*emis
      emis14 = r14*emis

      sibg%l(myl)%pooldt%loss_fire_lay(pidx,1) = emis*dtisib
      sibg%l(myl)%pooldt%loss_transf_lay(pidx,1) = dzero
      sibg%l(myl)%pooldt%poollu_dloss(pidx,1) = &
           sibg%l(myl)%pooldt%poollu_dloss(pidx,1) + emis

      sibg%l(myl)%pooldt%loss_fire_lay(p13idx,1) = emis13*dtisib
      sibg%l(myl)%pooldt%loss_transf_lay(p13idx,1) = dzero
      sibg%l(myl)%pooldt%poollu_dloss(p13idx,1) = &
           sibg%l(myl)%pooldt%poollu_dloss(p13idx,1) + emis13

      sibg%l(myl)%pooldt%loss_fire_lay(p14idx,1) = emis14*dtisib
      sibg%l(myl)%pooldt%loss_transf_lay(p14idx,1) = dzero
      sibg%l(myl)%pooldt%poollu_dloss(p14idx,1) = &
           sibg%l(myl)%pooldt%poollu_dloss(p14idx,1) + emis14

      flosslu(myl,pidx) = emis
      flosslu(myl,p13idx) = emis13
      flosslu(myl,p14idx) = emis14
      flosslu_dm(myl,pidx) = emis*drytoc
      flosslu_dm(myl,p13idx) = emis13*drytoc
      flosslu_dm(myl,p14idx) = emis14*drytoc

      curemis = curemis + emis
      curemis_dm = curemis_dm + emis*drytoc
      curemisc13 = curemisc13 + emis13
      curemisc13_dm = curemisc13_dm + emis13*drytoc
      curemisc14 = curemisc14 + emis14
      curemisc14_dm = curemisc14_dm + emis14*drytoc
   end subroutine burn_dead_single


   subroutine burn_dead_multi(pidx, p13idx, p14idx, fidx)
      integer(i4), intent(in) :: pidx, p13idx, p14idx, fidx
      integer(i4) :: ss
      real(r8) :: avail_sum, ccomp, target_emis, layer_emis
      real(r8) :: emis, emis13, emis14, r13, r14
      real(r8), dimension(nsoil) :: avail_lay

      avail_lay(:) = dzero
      do ss=1,nsoil
         avail_lay(ss) = MAX(dzero, &
              sibg%l(myl)%pooldt%poollu_lay(pidx,ss) &
              + sibg%l(myl)%pooldt%poollu_dgain(pidx,ss) &
              - sibg%l(myl)%pooldt%poollu_dloss(pidx,ss))
      enddo
      avail_sum = sum(avail_lay)

      ccomp = firecon(mylp)%ccl(fidx) + &
           st_scalar*(firecon(mylp)%cch(fidx) - firecon(mylp)%ccl(fidx))
      target_emis = MIN(avail_sum, MAX(dzero, avail_sum*myba*ccomp))

      emis = dzero
      emis13 = dzero
      emis14 = dzero
      do ss=1,nsoil
         layer_emis = MIN(avail_lay(ss), MAX(dzero, &
              target_emis*sibg%l(myl)%pooldt%poollu_flay(pidx,ss)))
         r13 = sibg%l(myl)%pooldt%rcpoollu_lay(p13idx,ss)
         r14 = sibg%l(myl)%pooldt%rcpoollu_lay(p14idx,ss)

         sibg%l(myl)%pooldt%loss_fire_lay(pidx,ss) = layer_emis*dtisib
         sibg%l(myl)%pooldt%loss_transf_lay(pidx,ss) = dzero
         sibg%l(myl)%pooldt%poollu_dloss(pidx,ss) = &
              sibg%l(myl)%pooldt%poollu_dloss(pidx,ss) + layer_emis

         sibg%l(myl)%pooldt%loss_fire_lay(p13idx,ss) = &
              r13*layer_emis*dtisib
         sibg%l(myl)%pooldt%loss_transf_lay(p13idx,ss) = dzero
         sibg%l(myl)%pooldt%poollu_dloss(p13idx,ss) = &
              sibg%l(myl)%pooldt%poollu_dloss(p13idx,ss) + r13*layer_emis

         sibg%l(myl)%pooldt%loss_fire_lay(p14idx,ss) = &
              r14*layer_emis*dtisib
         sibg%l(myl)%pooldt%loss_transf_lay(p14idx,ss) = dzero
         sibg%l(myl)%pooldt%poollu_dloss(p14idx,ss) = &
              sibg%l(myl)%pooldt%poollu_dloss(p14idx,ss) + r14*layer_emis

         emis = emis + layer_emis
         emis13 = emis13 + r13*layer_emis
         emis14 = emis14 + r14*layer_emis
      enddo

      flosslu(myl,pidx) = emis
      flosslu(myl,p13idx) = emis13
      flosslu(myl,p14idx) = emis14
      flosslu_dm(myl,pidx) = emis*drytoc
      flosslu_dm(myl,p13idx) = emis13*drytoc
      flosslu_dm(myl,p14idx) = emis14*drytoc

      curemis = curemis + emis
      curemis_dm = curemis_dm + emis*drytoc
      curemisc13 = curemisc13 + emis13
      curemisc13_dm = curemisc13_dm + emis13*drytoc
      curemisc14 = curemisc14 + emis14
      curemisc14_dm = curemisc14_dm + emis14*drytoc
   end subroutine burn_dead_multi

end subroutine firegfed_interp

