!=======================================================
subroutine phen_update( &
           gnum, lnum, pnum, ipft, &
           phencont, poolcont, &
           daynew, daylmax, daylen, daylendt, &
           clim_cupr, clim_precip, &
           cupr, lspr, tmdf, tm, &
           assim, rstfac2, rstfac4, &
           lai, hydrovt, phent, pooldt, poollt, vmax, &
           c13assim, physcont, fract, hold_stage1, trace_minveg_stage, trace_doy)
!=======================================================

! Description
! -----------
! - Updates all phenology variables
! - Determines the phenology stage
! - Set stage-dependent variables
!-----------

use kinds
use module_param, only: &
    phen_param, pool_param, &
    phys_param
use module_pftinfo, only: &
    pft_pmeth, pmeth_stg
use module_sib, only: &
    hydrov_type, phen_type, &
    poold_type, pooll_type,&
    fract_type
use module_time, only: &
     steps_per_day, dtsib, wt_clim

implicit none

!...input variables

integer(i4), intent(in) :: gnum, lnum, pnum
integer(i4), intent(inout) :: ipft
type(phen_param), intent(in) :: phencont
type(pool_param), intent(in) :: poolcont
logical, intent(in) :: daynew
real(r4), intent(in) :: daylmax
real(r8), intent(in) :: daylen, daylendt
real(r8), intent(in) :: clim_cupr, clim_precip
real(r8), intent(in) :: cupr, lspr, tmdf, tm
real(r8), intent(in) :: assim, rstfac2, rstfac4, lai, c13assim
real(r8), intent(inout) :: vmax

type(hydrov_type), intent(in) :: hydrovt
type(phen_type), intent(inout) :: phent
type(poold_type), intent(inout) :: pooldt
type(pooll_type), intent(inout) :: poollt
type(phys_param), intent(in) :: physcont
type(fract_type), intent(in) :: fract
logical, intent(in) :: hold_stage1, trace_minveg_stage
integer(i4), intent(in) :: trace_doy

!...local variables
integer(i4) :: ips, stage_before
real(r8) :: lai_min_before, lai_max_before, phen_pi_before
real(r8) :: climp, climw, wapot


!----------------------------------------------------------
!...Set local variables
select case (phencont%cwa_type)
   case (1)
       climw = clim_cupr
   case (2)
       climw = clim_precip
   case (3)
       climw = hydrovt%clim_pawfrw
   case (4)
       climw = hydrovt%clim_tawfrw
   case default
       climw = dzero
   end select
climp = phencont%climp_a*(phencont%climp_b**climw) &
         + phencont%climp_c*(climw - phencont%climp_d)

select case (phencont%pswx_type)
   case (1)
       wapot = hydrovt%pawftop
   case (2)
       wapot = hydrovt%pawfzw
   case (3)
       wapot = MAX(0.2,hydrovt%pawfzw)
       wapot = MIN(1.0, 2.*wapot)
   case (4)
       if (hydrovt%pawfzw .gt. 0.) then
           wapot = 1.0
       else
           wapot = 0.0
       endif
    case (5)
       wapot = hydrovt%tawftop
    case (6)
       wapot = hydrovt%tawfzw
    case (7)
       wapot = rstfac2
    case (8)
       wapot = rstfac4
    case (9)
       wapot = hydrovt%pawfrw
    case default
       wapot = dzero
end select

!------------------------------------------
!...print check

!print*,' '
!print*,'assim from phen_update: ',assim
!print*,'c13assim from phen_update: ',c13assim 
!print*,' '

!...Update phenology running-mean variables
phent%phenave_assim = assim*phencont%wt_assim &
        + phent%phenave_assim*(1. - phencont%wt_assim)
phent%phenave_assimsm = &
    MAX(phent%phenave_assim, phent%phenave_assimsm)

phent%phenave_pr = &
     (cupr + lspr)*phencont%wt_precip*dtsib*steps_per_day &
     + phent%phenave_pr*(1. - phencont%wt_precip)
IF (clim_precip .GT. dzero) THEN
    phent%phenave_prpot = MIN(done, phent%phenave_pr / clim_precip)
ENDIF
 
phent%phenave_tawftop = hydrovt%tawftop*phencont%wt_tawftop &
    + phent%phenave_tawftop*(1. - phencont%wt_tawftop)
phent%phenave_tm = tm*phencont%wt_tm &
    + phent%phenave_tm*(1. - phencont%wt_tm)

! C2 DIAGNOSTIC: seed phenological LAI reference from the current
! environmentally defined climatological suitability, using exactly the
! relationships used by the normal running-mean update below.
! This applies only on the first day of a new minimum-vegetation experiment.
! No stage-selection equation, live-carbon pool, or soil state is modified.
if (hold_stage1 .and. daynew .and. &
    pft_pmeth(pnum) .eq. pmeth_stg) then
   phent%phenc_climp = MAX(MIN(climp,phencont%climp_max), &
                           phencont%climp_min)
   phent%phenc_laimax = phent%phenc_climp*phencont%clai_coef &
                        + phencont%clai_offg
   phent%phenc_laimin = phent%phenc_climp*phencont%clai_coef &
                        + phencont%clai_offl
   if (trace_minveg_stage .and. gnum .eq. 1 .and. lnum .eq. 1) then
      write(*,'(a,1x,4(es16.8,1x))') &
       'MINVEG_C2_SEED climp,refclimp,laimin,laimax',climp, &
        phent%phenc_climp,phent%phenc_laimin,phent%phenc_laimax
   endif
endif

phent%phenc_climp = (1.-wt_clim)*phent%phenc_climp &
    + wt_clim*MAX(MIN(climp,phencont%climp_max), phencont%climp_min)
phent%phenc_laimax = (1.-wt_clim)*phent%phenc_laimax &
    + wt_clim*(phent%phenc_climp*phencont%clai_coef + phencont%clai_offg)
phent%phenc_laimin = (1.-wt_clim)*phent%phenc_laimin &
    + wt_clim*(phent%phenc_climp*phencont%clai_coef + phencont%clai_offl)

phent%phenave_env = rstfac4*phencont%wt_pswx &
    + phent%phenave_env*(1.-phencont%wt_pswx)
phent%phenave_wa = wapot*phencont%wt_pswx &
    + phent%phenave_wa*(1.-phencont%wt_pswx)
phent%phenave_wac = 0.5 &
    * (phent%phenave_wa + phent%phenave_env)
phent%phenave_wacsm = &
    MAX(phent%phenave_wac, phent%phenave_wacsm)


!---------------------------------
!...Update phenology once daily
if (daynew) then
   !...Calculate growing season start
   call phen_gss(phencont, daylmax, daylen, daylendt, phent)

   !...Calculate the phenology stage index
   !.....and set the phenology stage
   if (pft_pmeth(pnum) .eq. pmeth_stg) then
       stage_before = phent%phen_istage
       lai_min_before = phent%phenc_laimin
       lai_max_before = phent%phenc_laimax
       phen_pi_before = phent%phen_pi
       if (.not. hold_stage1) then
          call phen_dynamic( phencont, &
               daylmax, daylen, daylendt, lai, phent)
       else
          ! Intentional only on diagnostic initialization day.  Keep
          ! stage/PI from setup_phent while history potentials update.
          phent%phen_istage = 1
          phent%phen_pi = 1.d0
       endif
       if (trace_minveg_stage .and. gnum .eq. 1 .and. lnum .eq. 1) then
          write(*,'(a,1x,i4,1x,i4,1x,i3,1x,i3,1x,l1)') &
               'MINVEG_STAGE day,pft,before,after,held', &
               trace_doy,pnum,stage_before,phent%phen_istage,hold_stage1
          write(*,'(a,1x,9(es16.8,1x))') &
               'MINVEG_STAGE pi_before,pi,dayl,grw,wx,lai,laimin,laimax,psgmin', &
               phen_pi_before,phent%phen_pi,phent%phens_dayl, &
               phent%phens_grw,phent%phens_wx,lai, &
               lai_min_before,lai_max_before,phencont%psg_min
          write(*,'(a,1x,4(es16.8,1x))') &
               'MINVEG_STAGE threshp_1_4',phencont%threshp(1:4)
          write(*,'(a,1x,7(es16.8,1x))') &
               'MINVEG_STAGE phenave_wac,wacsm,env,wa,climp,daylen,daylendt', &
               phent%phenave_wac,phent%phenave_wacsm, &
               phent%phenave_env,phent%phenave_wa, &
               phent%phenc_climp,daylen,daylendt
          write(*,'(a,1x,3(es16.8,1x))') &
               'MINVEG_STAGE laimin_after,laimax_after,assimhist', &
               phent%phenc_laimin,phent%phenc_laimax, &
               phent%phenave_assim
       endif
   else
       call phen_defined( &
             gnum, lnum, ipft, pnum, &
             phencont, poolcont, &
             lai, tmdf, phent, pooldt, &
             poollt, physcont)
   endif !phenology method choices

   !...Set phenologically-varying variables
   ips = phent%phen_istage

   !.....Growing Season Counters
   IF ((ips .ge. ione) .and. (ips .le. phencont%npstg-1)) THEN
      phent%nd_gs = phent%nd_gs + ione
      phent%nd_stg(ips) = phent%nd_stg(ips) + ione
   ELSE
      phent%nd_dormant = phent%nd_dormant + ione
   ENDIF

   !.....Photosynthetic Allocation
   poollt%alloc_phen(:) = phencont%allocp(:,ips)

   !.....Dynamical Allocation
   !.....not allowed during leaf-out or senescence
   poollt%aadj_moist = phencont%adj_moist
   poollt%aadj_temp = phencont%adj_temp
   IF ((ips .le. 1) .or. (ips .eq. phencont%npstg)) THEN
       poollt%aadj_moist = .false.
       poollt%aadj_temp = .false.
   ENDIF

   !.....Leaf pool transfer
   poollt%tfl_pstage = phencont%lptransfer(ips)
   
   !.....Vmax
   vmax = phencont%vmax0(ips)

endif !daynew

end subroutine phen_update
