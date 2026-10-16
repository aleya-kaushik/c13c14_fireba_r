!-----------------------------------------------------------------------
! Standard GFED emissions-driven fire, S0-W within-PFT sensitivity.
! Initial PFT demand remains area-proportional. Original fire_*f factors
! specify PREFERRED pool shares, not hard maxima on shares of emissions.
! Any unfulfilled demand is reallocated among accessible SURFACE pools
! within the SAME PFT before any cross-PFT redistribution. Soil and
! root combustion are disabled. Availability-limited upper-bound test:
! no burned-area or pool-specific completeness cap is inferred.
!
! Emission requests and pool losses are mol C m-2 of GRID CELL per step,
! retaining the established SiB4 LU fire accounting.  No additional
! larea multiplier is applied to emis_fire or the losses below.
!
! Initial requests follow the areas of all non-bare LUs; priority affects
! traversal, not the independently assigned shares.  Shortfalls are
! redistributed among LUs with accessible SURFACE fuel before any soil
! fallback.  Root combustion is disabled in both configurations.
!-----------------------------------------------------------------------
subroutine fire_interp(indx, lon, lat, sibg)
   use kinds
   use module_io, only: fire_step, fire_seccur, fire_secnext
   use module_oparams, only: fire_leaff, fire_stwdf, &
                             fire_cdbf, fire_metlf, fire_strlf
   use module_param, only: poolcon
   use module_pparams, only: month_names, pdb
   use module_pftinfo
   use module_poolinfo
   use module_sib, only: gridcell_type
   use module_sibconst, only: npoolpft, npoollu, nsoil, &
                              fireb_print, fireb_stop, fireb_thresh
   use module_time, only: month, day, year, dtisib, dtsib, &
                          sec_tot, wt_daily
   implicit none

   integer(i4), intent(in) :: indx
   real(r4), intent(in) :: lon, lat
   type(gridcell_type), intent(inout) :: sibg

   ! S0-W: never combust soil or roots. S1 fallback remains disabled.
   logical, parameter :: fire_soil_fallback = .false.
   integer(i4), parameter :: bare_pft_id = 1
   real(r8), parameter :: tiny_fire = 1.0d-12
   integer(i4), parameter :: max_redistribution_passes = 32

   integer(i4) :: ntpft, neligible, l, k, myl, best, pass
   integer(i4) :: lp, wp, pp, cdbp, metlp, strlp, slitp, slowp, armp
   integer(i4) :: lpc13, wpc13, ppc13, cdbpc13, metlpc13, strlpc13
   integer(i4) :: slitpc13, slowpc13, armpc13
   integer(i4) :: lpc14, wpc14, ppc14, cdbpc14, metlpc14, strlpc14
   integer(i4) :: slitpc14, slowpc14, armpc14
   integer(i4), allocatable :: tpref(:), tpnum(:), tsortref(:)
   real(r4), allocatable :: tparea(:)
   character(len=clen), allocatable :: tpname(:)
   logical, allocatable :: chosen(:), available(:), can_redistribute(:)
   real(r8), allocatable :: request(:), share(:)
   real(r8), allocatable :: flosspft(:,:), flosslu(:,:)
   real(r8) :: facsibdrv, totemis, tfarea, remainder, before
   real(r8) :: allocated, removed, area_sum
   real(r8) :: actualemis, actualemisc13, actualemisc14
   real(r8) :: total_actual, rcemispoolc13, residual, close_c
   logical :: fb_err

   ! Keep the original C/C13/C14 indexing conventions.
   lp = pool_indx_leaf
   wp = pool_indx_stwd
   pp = pool_indx_prod
   lpc13 = pool_indx_leaf_c13 - npoollu/3
   wpc13 = pool_indx_stwd_c13 - npoollu/3
   ppc13 = pool_indx_prod_c13 - npoollu/3
   lpc14 = pool_indx_leaf_c14 - 2*(npoollu/3)
   wpc14 = pool_indx_stwd_c14 - 2*(npoollu/3)
   ppc14 = pool_indx_prod_c14 - 2*(npoollu/3)

   cdbp = pool_indx_cdb - npoolpft/3
   metlp = pool_indx_metl - npoolpft/3
   strlp = pool_indx_strl - npoolpft/3
   slitp = pool_indx_slit - npoolpft/3
   slowp = pool_indx_slow - npoolpft/3
   armp = pool_indx_arm - npoolpft/3
   cdbpc13 = pool_indx_cdb_c13 - 2*(npoolpft/3)
   metlpc13 = pool_indx_metl_c13 - 2*(npoolpft/3)
   strlpc13 = pool_indx_strl_c13 - 2*(npoolpft/3)
   slitpc13 = pool_indx_slit_c13 - 2*(npoolpft/3)
   slowpc13 = pool_indx_slow_c13 - 2*(npoolpft/3)
   armpc13 = pool_indx_arm_c13 - 2*(npoolpft/3)
   cdbpc14 = pool_indx_cdb_c14 - npoolpft
   metlpc14 = pool_indx_metl_c14 - npoolpft
   strlpc14 = pool_indx_strl_c14 - npoolpft
   slitpc14 = pool_indx_slit_c14 - npoolpft
   slowpc14 = pool_indx_slow_c14 - npoolpft
   armpc14 = pool_indx_arm_c14 - npoolpft

   do l=1,sibg%g_nlu
      sibg%l(l)%poollt%emis_fire = dzero
      sibg%l(l)%poollt%fire_initial_request = dzero
      sibg%l(l)%poollt%fire_surface_initial = dzero
      sibg%l(l)%poollt%fire_surface_redist = dzero
      sibg%l(l)%poollt%fire_soil_fallback = dzero
      sibg%l(l)%poollt%rmmd_fire = dzero
      sibg%l(l)%poollt%emis_firec13 = dzero
      sibg%l(l)%poollt%rmmd_firec13 = dzero
      sibg%l(l)%poollt%emis_firec14 = dzero
      sibg%l(l)%poollt%rmmd_firec14 = dzero
      sibg%l(l)%fract%d13cemisfire_pool = dzero
      sibg%l(l)%poollt%loss_fire_lay(:,:) = dzero
      sibg%l(l)%pooldt%loss_fire_lay(:,:) = dzero
   enddo

   if (fire_step <= izero) return
   facsibdrv = dble(fire_seccur-sec_tot)/dble(fire_step)
   if ((facsibdrv > 1.0d0) .or. (facsibdrv < 0.0d0)) then
      sibg%gprogt%firec = dzero
      sibg%gprogt%fireco2 = dzero
      return
   endif
   sibg%gprogt%firec = facsibdrv*sibg%gprogt%firec1 + &
        (1.0d0-facsibdrv)*sibg%gprogt%firec2
   sibg%gprogt%fireco2 = facsibdrv*sibg%gprogt%fireco21 + &
        (1.0d0-facsibdrv)*sibg%gprogt%fireco22
   if (sibg%gprogt%firec <= dzero) return

   totemis = sibg%gprogt%firec*dtsib
   ntpft = sibg%g_nlu
   if (ntpft <= 0) return
   allocate(tpref(ntpft),tpnum(ntpft),tparea(ntpft),tpname(ntpft))
   allocate(tsortref(ntpft),chosen(ntpft),available(ntpft), &
            can_redistribute(ntpft))
   allocate(request(ntpft),share(ntpft))
   allocate(flosspft(ntpft,npoolpft),flosslu(ntpft,npoollu))
   flosspft = dzero
   flosslu = dzero
   tpref = sibg%l(1:ntpft)%ipft
   tpnum = pft_num(tpref)
   tparea = sibg%l(1:ntpft)%larea
   do l=1,ntpft
      tpname(l) = pft_name(tpnum(l))
   enddo
   chosen = .false.
   available = .false.
   can_redistribute = .false.
   request = dzero
   share = dzero

   ! Eligibility means vegetated land area; fuel constraints are applied
   ! later.  Bare LU area is excluded, and all non-bare LUs are retained.
   tfarea = dzero
   do l=1,ntpft
      if (tpref(l) == bare_pft_id .or. tparea(l) <= 0.0_r4) cycle
      available(l) = .true.
      tfarea = tfarea + dble(tparea(l))
   enddo

   if (tfarea <= tiny_fire) then
      ! There is no vegetation to burn.  Record only an unmatched target;
      ! no bare LU carbon is removed or assigned as realized emissions.
      sibg%l(1)%poollt%rmmd_fire = totemis*dtisib
      sibg%l(1)%poollt%rmmd_firec13 = &
           sibg%l(1)%fract%rcpoolfirec13*totemis*dtisib
      sibg%l(1)%poollt%rmmd_firec14 = &
           sibg%l(1)%fract%rcpoolfirec14*totemis*dtisib
      return
   endif

   ! Rank: both grasses form the same priority class. C4 before C3 is
   ! solely a stable traversal convention, not preferential allocation.
   ! All other PFTs follow in descending areal coverage (stable tie order).
   neligible = 0
   do l=1,ntpft
      if (.not.available(l)) cycle
      if (trim(tpname(l)) /= 'c4g') cycle
      neligible = neligible+1
      tsortref(neligible) = l
      chosen(l) = .true.
   enddo
   do l=1,ntpft
      if (.not.available(l) .or. chosen(l)) cycle
      if (trim(tpname(l)) /= 'c3g') cycle
      neligible = neligible+1
      tsortref(neligible) = l
      chosen(l) = .true.
   enddo
   do
      best=0
      do l=1,ntpft
         if (.not.available(l) .or. chosen(l)) cycle
         if (best == 0) then
            best=l
         else if (tparea(l) > tparea(best)) then
            best=l
         endif
      enddo
      if (best == 0) exit
      neligible=neligible+1
      tsortref(neligible)=best
      chosen(best)=.true.
   enddo

   ! Pass 1: ALL vegetated PFTs receive an area-proportional request.
   ! Their processing order cannot change another PFT's assigned share.
   do k=1,neligible
      myl=tsortref(k)
      sibg%l(myl)%poollt%nd_fire = &
           sibg%l(myl)%poollt%nd_fire + wt_daily
      request(myl) = totemis*dble(tparea(myl))/tfarea
      sibg%l(myl)%poollt%fire_initial_request = request(myl)*dtisib
      call burn_surface(myl,request(myl),removed)
      sibg%l(myl)%poollt%fire_surface_initial = removed*dtisib
      ! A PFT that cannot fulfill its area-based share AFTER within-PFT
      ! redistribution donates only the remaining unfulfilled shortfall.
      can_redistribute(myl) = (request(myl)-removed <= tiny_fire) .and. &
                               (surface_capacity(myl) > tiny_fire)
   enddo

   ! Pass 2: distribute unmet GRID-CELL emissions among remaining surface
   ! fuels.  Each pass assigns independent shares before burning; loop
   ! order does not advantage C4 relative to C3 or any other PFT.
   ! Only PFTs that fulfilled their assigned request and still have
   ! accessible surface fuel can receive more. S0-W already redistributed
   ! inside each PFT before any cross-PFT share is calculated.
   do pass=1,max_redistribution_passes
      total_actual = sum(flosspft(:,1:npoolpft/3)) + &
                     sum(flosslu(:,1:npoollu/3))
      remainder = max(dzero,totemis-total_actual)
      if (remainder <= tiny_fire) exit
      area_sum=dzero
      do k=1,neligible
         myl=tsortref(k)
         if (can_redistribute(myl) .and. &
             surface_capacity(myl) > tiny_fire) &
            area_sum=area_sum+dble(tparea(myl))
      enddo
      if (area_sum <= tiny_fire) exit
      share=dzero
      do k=1,neligible
         myl=tsortref(k)
         if (.not.can_redistribute(myl)) cycle
         if (surface_capacity(myl) <= tiny_fire) cycle
         share(myl)=remainder*dble(tparea(myl))/area_sum
      enddo
      before=total_actual
      do k=1,neligible
         myl=tsortref(k)
         if (share(myl) <= tiny_fire) cycle
         call burn_surface(myl,share(myl),removed)
         sibg%l(myl)%poollt%fire_surface_redist = &
              sibg%l(myl)%poollt%fire_surface_redist + removed*dtisib
         if (share(myl)-removed > tiny_fire .or. &
             surface_capacity(myl) <= tiny_fire) &
             can_redistribute(myl)=.false.
      enddo
      total_actual = sum(flosspft(:,1:npoolpft/3)) + &
                     sum(flosslu(:,1:npoollu/3))
      if (total_actual-before <= tiny_fire) exit
   enddo

   ! S1 only: soil litter/slow/arm carbon may provide residual emissions,
   ! AFTER every surface allocation and redistribution.  No roots burn.
   if (fire_soil_fallback) then
      do pass=1,max_redistribution_passes
         total_actual = sum(flosspft(:,1:npoolpft/3)) + &
                        sum(flosslu(:,1:npoollu/3))
         remainder = max(dzero,totemis-total_actual)
         if (remainder <= tiny_fire) exit
         area_sum=dzero
         do k=1,neligible
            myl=tsortref(k)
            if (soil_capacity(myl) > tiny_fire) &
               area_sum=area_sum+dble(tparea(myl))
         enddo
         if (area_sum <= tiny_fire) exit
         share=dzero
         do k=1,neligible
            myl=tsortref(k)
            if (soil_capacity(myl) <= tiny_fire) cycle
            share(myl)=remainder*dble(tparea(myl))/area_sum
         enddo
         before=total_actual
         do k=1,neligible
            myl=tsortref(k)
            if (share(myl) <= tiny_fire) cycle
            call burn_soil(myl,share(myl),removed)
            sibg%l(myl)%poollt%fire_soil_fallback = &
                 sibg%l(myl)%poollt%fire_soil_fallback + removed*dtisib
         enddo
         total_actual = sum(flosspft(:,1:npoolpft/3)) + &
                        sum(flosslu(:,1:npoollu/3))
         if (total_actual-before <= tiny_fire) exit
      enddo
   endif

   ! Final per-LU output: emissions equal actual C/C13/C14 pool losses.
   ! The remaining GRID-CELL target is now redistributed only as a
   ! diagnostic mismatch (by eligible LU area), never combusted twice.
   total_actual = dzero
   do k=1,neligible
      myl=tsortref(k)
      actualemis=sum(flosspft(myl,1:npoolpft/3)) + &
                 sum(flosslu(myl,1:npoollu/3))
      actualemisc13=sum(flosspft(myl,npoolpft/3+1:2*npoolpft/3)) + &
                    sum(flosslu(myl,npoollu/3+1:2*npoollu/3))
      actualemisc14=sum(flosspft(myl,2*npoolpft/3+1:npoolpft)) + &
                    sum(flosslu(myl,2*npoollu/3+1:npoollu))
      total_actual=total_actual+actualemis
      sibg%l(myl)%poollt%emis_fire = actualemis*dtisib
      sibg%l(myl)%poollt%emis_firec13 = actualemisc13*dtisib
      sibg%l(myl)%poollt%emis_firec14 = actualemisc14*dtisib
      sibg%l(myl)%fract%d13cemisfire_pool=dzero
      if (actualemis > tiny_fire .and. actualemisc13 > dzero .and. &
          actualemisc13 < actualemis) then
         rcemispoolc13=actualemisc13/actualemis
         sibg%l(myl)%fract%d13cemisfire_pool = &
              (((rcemispoolc13/(1.0d0-rcemispoolc13))/pdb)-1.0d0)*1000.0d0
      endif
   enddo

   ! Confirm the staged diagnostic totals are consistent with actual fire losses.
   residual=max(dzero,totemis-total_actual)
   do k=1,neligible
      myl=tsortref(k)
      allocated=residual*dble(tparea(myl))/tfarea
      sibg%l(myl)%poollt%rmmd_fire = allocated*dtisib
      sibg%l(myl)%poollt%rmmd_firec13 = &
           sibg%l(myl)%fract%rcpoolfirec13*allocated*dtisib
      sibg%l(myl)%poollt%rmmd_firec14 = &
           sibg%l(myl)%fract%rcpoolfirec14*allocated*dtisib
   enddo

   close_c=totemis-total_actual-residual
   fb_err=abs(close_c) > fireb_thresh
   if (fb_err .or. fireb_print) then
      print*,' '
      print*,'---FIRE CARBON (staged allocation)---'
      print('(a,a,i3,a,i4)'), ' Date: ',trim(month_names(month)),day,', ',year
      print('(a,i6,2f8.2)'), ' Point/Lon/Lat: ',indx,lon,lat
      print('(a,i14)'), ' Model second: ',sec_tot
      print('(a,2i12)'), ' Current/Next Fire Second: ',fire_seccur,fire_secnext
      print('(a,l2)'), ' Soil sensitivity active: ',fire_soil_fallback
      print('(a,i4)'), ' Selected vegetated LUs: ',neligible
      print('(a,f18.8)'), ' Eligible area: ',tfarea
      print('(a,f18.8)'), ' Prescribed fire C (mol/m2): ',totemis
      print('(a,f18.8)'), ' Realized fire C (mol/m2): ',total_actual
      print('(a,f18.8)'), ' Unmet fire C (mol/m2): ',residual
      print('(a,f18.8)'), ' Carbon bookkeeping residual: ',close_c
      print('(a)'), ' PFT / area / realized C (mol/m2):'
      do k=1,neligible
         myl=tsortref(k)
         print('(i5,2f18.8)'),tpref(myl),dble(tparea(myl)), &
              sibg%l(myl)%poollt%emis_fire*dtsib
      enddo
      if (fb_err) print*,'!! Fire carbon allocation imbalance !!'
      if (fb_err .and. fireb_stop) stop
   endif

contains

   ! Pool availability is evaluated from CURRENT prognostic state minus
   ! prior losses in the same timestep.  Only positive-combustion-factor
   ! aboveground/surface pools count toward redistribution capacity.
   real(r8) function surface_capacity(i)
      integer(i4), intent(in) :: i
      surface_capacity=dzero
      if (fire_leaff > dzero) surface_capacity=surface_capacity+ &
           pft_avail(i,lp)
      if (fire_stwdf > dzero) surface_capacity=surface_capacity+ &
           pft_avail(i,wp)
      if (fire_metlf > dzero) surface_capacity=surface_capacity+ &
           lu_avail(i,metlp)
      if (fire_strlf > dzero) surface_capacity=surface_capacity+ &
           lu_avail(i,strlp)
      if (fire_cdbf > dzero) surface_capacity=surface_capacity+ &
           lu_avail(i,cdbp)
      surface_capacity=surface_capacity+pft_avail(i,pp)
   end function surface_capacity

   real(r8) function pft_avail(i,p)
      integer(i4), intent(in) :: i,p
      pft_avail=max(dzero,sibg%l(i)%poollt%poolpft(p) &
           -sibg%l(i)%poollt%poolpft_dloss(p,1) &
           -poolcon(tpnum(i))%poolpft_min(p))
   end function pft_avail

   real(r8) function lu_avail(i,p)
      integer(i4), intent(in) :: i,p
      lu_avail=max(dzero,sibg%l(i)%pooldt%poollu(p) &
           -sibg%l(i)%pooldt%poollu_dloss(p,1))
   end function lu_avail

   real(r8) function soil_capacity(i)
      integer(i4), intent(in) :: i
      soil_capacity=layer_available(i,slitp) &
                  +layer_available(i,slowp) &
                  +layer_available(i,armp)
   end function soil_capacity

   real(r8) function layer_available(i,p)
      integer(i4), intent(in) :: i,p
      integer(i4) :: j
      layer_available=dzero
      do j=1,nsoil
         layer_available=layer_available+max(dzero, &
              sibg%l(i)%pooldt%poollu_lay(p,j) &
             +sibg%l(i)%pooldt%poollu_dgain(p,j) &
             -sibg%l(i)%pooldt%poollu_dloss(p,j))
      enddo
   end function layer_available

   subroutine burn_surface(i,amount,taken)
      ! S0-W sensitivity: first honor the original pool preferences,
      ! then redistribute any deficit WITHIN this PFT to remaining
      ! eligible aboveground/dead/litter pools. The GFED request (amount)
      ! remains the hard overall emissions cap. poolpft_min still protects
      ! the existing minimum live C pools; roots and soil are untouched.
      ! This is NOT a combustion-completeness or burned-area model.
      integer(i4), intent(in) :: i
      real(r8), intent(in) :: amount
      real(r8), intent(out) :: taken
      real(r8) :: left, use_c, cap(5), wt(5), burn(5)
      real(r8) :: extra, weight_sum, used
      integer(i4) :: iter

      left=max(dzero,amount)
      if (left <= tiny_fire) then
         taken=dzero
         return
      endif

      ! Preferred first allocation: identical five-pool fractions to S0.
      use_c=min(left,max(dzero,dble(fire_leaff)*amount),pft_avail(i,lp))
      call take_pft(i,lp,lpc13,lpc14,use_c)
      left=left-use_c
      use_c=min(left,max(dzero,dble(fire_stwdf)*amount),pft_avail(i,wp))
      call take_pft(i,wp,wpc13,wpc14,use_c)
      left=left-use_c
      use_c=min(left,max(dzero,dble(fire_metlf)*amount),lu_avail(i,metlp))
      call take_lu_surface(i,metlp,metlpc13,metlpc14,use_c)
      left=left-use_c
      use_c=min(left,max(dzero,dble(fire_strlf)*amount),lu_avail(i,strlp))
      call take_lu_surface(i,strlp,strlpc13,strlpc14,use_c)
      left=left-use_c
      use_c=min(left,max(dzero,dble(fire_cdbf)*amount),lu_avail(i,cdbp))
      call take_lu_surface(i,cdbp,cdbpc13,cdbpc14,use_c)
      left=left-use_c

      ! Reassign unmet demand to ACTUALLY AVAILABLE surface material
      ! within the same PFT, before considering another PFT. Weights
      ! preserve the original relative preferences among usable pools;
      ! iterative clipping frees shortfalls from exhausted pools.
      ! A zero fire_*f factor continues to disable combustion of its pool.
      wt=(/max(dzero,dble(fire_leaff)), &
           max(dzero,dble(fire_stwdf)), &
           max(dzero,dble(fire_metlf)), &
           max(dzero,dble(fire_strlf)), &
           max(dzero,dble(fire_cdbf))/)
      do iter=1,6
         if (left <= tiny_fire) exit
         cap=(/pft_avail(i,lp),pft_avail(i,wp), &
               lu_avail(i,metlp),lu_avail(i,strlp), &
               lu_avail(i,cdbp)/)
         where (cap <= tiny_fire .or. wt <= dzero)
            cap=dzero
         end where
         weight_sum=sum(wt,mask=cap>dzero)
         if (weight_sum <= dzero) exit
         extra=min(left,sum(cap))
         if (extra <= tiny_fire) exit
         burn=dzero
         where (cap > dzero)
            burn=min(cap,extra*wt/weight_sum)
         end where
         used=sum(burn)
         if (used <= tiny_fire) exit
         call take_pft(i,lp,lpc13,lpc14,burn(1))
         call take_pft(i,wp,wpc13,wpc14,burn(2))
         call take_lu_surface(i,metlp,metlpc13,metlpc14,burn(3))
         call take_lu_surface(i,strlp,strlpc13,strlpc14,burn(4))
         call take_lu_surface(i,cdbp,cdbpc13,cdbpc14,burn(5))
         left=max(dzero,left-used)
      enddo

      ! Product remains an explicitly tracked last-resort surface pool.
      ! Do not access fine/coarse roots, soil litter, slow, or armored C.
      if (left > tiny_fire) then
         use_c=min(left,pft_avail(i,pp))
         call take_pft(i,pp,ppc13,ppc14,use_c)
         left=left-use_c
      endif
      taken=amount-max(dzero,left)
   end subroutine burn_surface

   subroutine take_pft(i,p,p13,p14,amount)
      integer(i4), intent(in) :: i,p,p13,p14
      real(r8), intent(in) :: amount
      real(r8) :: a13,a14
      if (amount <= dzero) return
      a13=amount*sibg%l(i)%poollt%rcpoolpft_lay(p13,1)
      a14=amount*sibg%l(i)%poollt%rcpoolpft_lay(p14,1)
      sibg%l(i)%poollt%loss_fire_lay(p,1) = &
           sibg%l(i)%poollt%loss_fire_lay(p,1)+amount*dtisib
      sibg%l(i)%poollt%poolpft_dloss(p,1) = &
           sibg%l(i)%poollt%poolpft_dloss(p,1)+amount
      sibg%l(i)%poollt%loss_fire_lay(p13,1) = &
           sibg%l(i)%poollt%loss_fire_lay(p13,1)+a13*dtisib
      sibg%l(i)%poollt%poolpft_dloss(p13,1) = &
           sibg%l(i)%poollt%poolpft_dloss(p13,1)+a13
      sibg%l(i)%poollt%loss_fire_lay(p14,1) = &
           sibg%l(i)%poollt%loss_fire_lay(p14,1)+a14*dtisib
      sibg%l(i)%poollt%poolpft_dloss(p14,1) = &
           sibg%l(i)%poollt%poolpft_dloss(p14,1)+a14
      flosspft(i,p)=flosspft(i,p)+amount
      flosspft(i,p13)=flosspft(i,p13)+a13
      flosspft(i,p14)=flosspft(i,p14)+a14
   end subroutine take_pft

   subroutine take_lu_surface(i,p,p13,p14,amount)
      integer(i4), intent(in) :: i,p,p13,p14
      real(r8), intent(in) :: amount
      real(r8) :: a13,a14
      if (amount <= dzero) return
      a13=amount*sibg%l(i)%pooldt%rcpoollu_lay(p13,1)
      a14=amount*sibg%l(i)%pooldt%rcpoollu_lay(p14,1)
      sibg%l(i)%pooldt%loss_fire_lay(p,1) = &
           sibg%l(i)%pooldt%loss_fire_lay(p,1)+amount*dtisib
      sibg%l(i)%pooldt%poollu_dloss(p,1) = &
           sibg%l(i)%pooldt%poollu_dloss(p,1)+amount
      sibg%l(i)%pooldt%loss_fire_lay(p13,1) = &
           sibg%l(i)%pooldt%loss_fire_lay(p13,1)+a13*dtisib
      sibg%l(i)%pooldt%poollu_dloss(p13,1) = &
           sibg%l(i)%pooldt%poollu_dloss(p13,1)+a13
      sibg%l(i)%pooldt%loss_fire_lay(p14,1) = &
           sibg%l(i)%pooldt%loss_fire_lay(p14,1)+a14*dtisib
      sibg%l(i)%pooldt%poollu_dloss(p14,1) = &
           sibg%l(i)%pooldt%poollu_dloss(p14,1)+a14
      flosslu(i,p)=flosslu(i,p)+amount
      flosslu(i,p13)=flosslu(i,p13)+a13
      flosslu(i,p14)=flosslu(i,p14)+a14
   end subroutine take_lu_surface

   subroutine burn_soil(i,amount,taken)
      integer(i4), intent(in) :: i
      real(r8), intent(in) :: amount
      real(r8), intent(out) :: taken
      real(r8) :: left, use_c
      left=max(dzero,amount)
      use_c=min(left,layer_available(i,slitp))
      call take_lu_layered(i,slitp,slitpc13,slitpc14,use_c)
      left=left-use_c
      use_c=min(left,layer_available(i,slowp))
      call take_lu_layered(i,slowp,slowpc13,slowpc14,use_c)
      left=left-use_c
      use_c=min(left,layer_available(i,armp))
      call take_lu_layered(i,armp,armpc13,armpc14,use_c)
      left=left-use_c
      taken=amount-max(dzero,left)
   end subroutine burn_soil

   subroutine take_lu_layered(i,p,p13,p14,amount)
      integer(i4), intent(in) :: i,p,p13,p14
      real(r8), intent(in) :: amount
      real(r8) :: avail(nsoil),burn(nsoil),whole,a13,a14
      integer(i4) :: j
      if (amount <= dzero) return
      avail=dzero
      do j=1,nsoil
         avail(j)=max(dzero,sibg%l(i)%pooldt%poollu_lay(p,j) &
              +sibg%l(i)%pooldt%poollu_dgain(p,j) &
              -sibg%l(i)%pooldt%poollu_dloss(p,j))
      enddo
      whole=sum(avail)
      if (whole <= dzero) return
      burn=amount*avail/whole
      do j=1,nsoil
         if (burn(j) <= dzero) cycle
         a13=burn(j)*sibg%l(i)%pooldt%rcpoollu_lay(p13,j)
         a14=burn(j)*sibg%l(i)%pooldt%rcpoollu_lay(p14,j)
         sibg%l(i)%pooldt%loss_fire_lay(p,j) = &
              sibg%l(i)%pooldt%loss_fire_lay(p,j)+burn(j)*dtisib
         sibg%l(i)%pooldt%poollu_dloss(p,j) = &
              sibg%l(i)%pooldt%poollu_dloss(p,j)+burn(j)
         sibg%l(i)%pooldt%loss_fire_lay(p13,j) = &
              sibg%l(i)%pooldt%loss_fire_lay(p13,j)+a13*dtisib
         sibg%l(i)%pooldt%poollu_dloss(p13,j) = &
              sibg%l(i)%pooldt%poollu_dloss(p13,j)+a13
         sibg%l(i)%pooldt%loss_fire_lay(p14,j) = &
              sibg%l(i)%pooldt%loss_fire_lay(p14,j)+a14*dtisib
         sibg%l(i)%pooldt%poollu_dloss(p14,j) = &
              sibg%l(i)%pooldt%poollu_dloss(p14,j)+a14
         flosslu(i,p)=flosslu(i,p)+burn(j)
         flosslu(i,p13)=flosslu(i,p13)+a13
         flosslu(i,p14)=flosslu(i,p14)+a14
      enddo
   end subroutine take_lu_layered

end subroutine fire_interp
