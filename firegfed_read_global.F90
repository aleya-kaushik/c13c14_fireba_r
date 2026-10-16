#include "nc_util.h"

subroutine firegfed_read_global()

!****--------------------------------------------------------------------
!    This subroutines reads the fire burned area data for the current
!    time step for global (group/multi-point) simulations.
!    If required, it closes the current month's data file and opens the 
!    next month's data file.
!****--------------------------------------------------------------------

use kinds
use netcdf

use module_io
use module_sib, only: sib
use module_sibconst, only: &
   nsib, subset, subcount, &
   print_fire, print_stop
use module_time, only: &
   month, year, day, hour
use module_pparams, only:  &
    secs_per_day

implicit none

!...data variables
real(r8), dimension(nsib) :: ba_dbgin, ba_en2in, ba_en3in, ba_dnfin, ba_eb1in
real(r8), dimension(nsib) :: ba_eb2in, ba_db1in, ba_db2in, ba_db3in, ba_shbin
real(r8), dimension(nsib) :: ba_shain, ba_c3ain, ba_c3gin, ba_c4gin, ba_c3cin
real(r8), dimension(nsib) :: ba_c4cin, ba_mzein, ba_soyin, ba_wwtin

!...netcdf variables
integer :: status
integer :: dimid, dimlen
integer :: ncyid, ncmid, ncdid
integer :: fba_dbgid, fba_en2id, fba_en3id, fba_dnfid, fba_eb1id, fba_eb2id
integer :: fba_db1id, fba_db2id, fba_db3id, fba_shbid, fba_shaid, fba_c3aid
integer :: fba_c3gid, fba_c4gid, fba_c3cid, fba_c4cid, fba_mzeid, fba_soyid, fba_wwtid
character(len=20) :: dim_name
integer, dimension(2) :: mstart, mcount

!...local variables
integer :: i
integer :: xyear, xmonth, xday
!real(r8)::xyear, xmonth, xday
logical :: fire_existf

!****--------------------------------------------------------------------
! storing previous time steps data
do i=1,subcount
    sib%g(i)%gprogt%ba_dbg1 = sib%g(i)%gprogt%ba_dbg2
    sib%g(i)%gprogt%ba_en21 = sib%g(i)%gprogt%ba_en22
    sib%g(i)%gprogt%ba_en31 = sib%g(i)%gprogt%ba_en32
    sib%g(i)%gprogt%ba_dnf1 = sib%g(i)%gprogt%ba_dnf2
    sib%g(i)%gprogt%ba_eb11 = sib%g(i)%gprogt%ba_eb12
    sib%g(i)%gprogt%ba_eb21 = sib%g(i)%gprogt%ba_eb22
    sib%g(i)%gprogt%ba_db11 = sib%g(i)%gprogt%ba_db12
    sib%g(i)%gprogt%ba_db21 = sib%g(i)%gprogt%ba_db22
    sib%g(i)%gprogt%ba_db31 = sib%g(i)%gprogt%ba_db32
    sib%g(i)%gprogt%ba_shb1 = sib%g(i)%gprogt%ba_shb2
    sib%g(i)%gprogt%ba_sha1 = sib%g(i)%gprogt%ba_sha2
    sib%g(i)%gprogt%ba_c3a1 = sib%g(i)%gprogt%ba_c3a2
    sib%g(i)%gprogt%ba_c3g1 = sib%g(i)%gprogt%ba_c3g2
    sib%g(i)%gprogt%ba_c4g1 = sib%g(i)%gprogt%ba_c4g2
    sib%g(i)%gprogt%ba_c3c1 = sib%g(i)%gprogt%ba_c3c2
    sib%g(i)%gprogt%ba_c4c1 = sib%g(i)%gprogt%ba_c4c2
    sib%g(i)%gprogt%ba_mze1 = sib%g(i)%gprogt%ba_mze2
    sib%g(i)%gprogt%ba_soy1 = sib%g(i)%gprogt%ba_soy2
    sib%g(i)%gprogt%ba_wwt1 = sib%g(i)%gprogt%ba_wwt2
enddo

!print*,'firegfed_read_global top fire_year:',fire_year
!print*,'firegfed_read_global top fire_month:',fire_month
!print*,'firegfed_read_global top fire_day:',fire_day

! switch files if needed
fire_existf = .false.
if (firegfed_switchf) then
   status = nf90_close(fireid)

   write(firegfed_filename,fmt='(a,i4.4,i2.2,a3)') trim(firegfed_path), &
         fire_year, fire_month, '.nc'
   inquire(file=trim(firegfed_filename), exist=fire_existf)
   if (fire_existf) then
      STATUS = nf90_open(trim(firegfed_filename), nf90_nowrite, fireid)
   else
      if (firefile_stop) then
          print*,''
          print*,'Missing/Non-Existant Fire Emissions File!'
          print*,'File: ',trim(firegfed_filename)
          print*,'Please check fr_path in namel_sibdrv.'
          print*,''
          STOP
        else
          fire_step = 0
          fireid = 0
          sib%g(:)%gprogt%ba_dbg2 = dzero
          sib%g(:)%gprogt%ba_en22 = dzero
          sib%g(:)%gprogt%ba_en32 = dzero
          sib%g(:)%gprogt%ba_dnf2 = dzero
          sib%g(:)%gprogt%ba_eb12 = dzero
          sib%g(:)%gprogt%ba_eb22 = dzero
          sib%g(:)%gprogt%ba_db12 = dzero
          sib%g(:)%gprogt%ba_db22 = dzero
          sib%g(:)%gprogt%ba_db32 = dzero
          sib%g(:)%gprogt%ba_shb2 = dzero
          sib%g(:)%gprogt%ba_sha2 = dzero
          sib%g(:)%gprogt%ba_c3a2 = dzero
          sib%g(:)%gprogt%ba_c3g2 = dzero
          sib%g(:)%gprogt%ba_c4g2 = dzero
          sib%g(:)%gprogt%ba_c3c2 = dzero
          sib%g(:)%gprogt%ba_c4c2 = dzero
          sib%g(:)%gprogt%ba_mze2 = dzero
          sib%g(:)%gprogt%ba_soy2 = dzero
          sib%g(:)%gprogt%ba_wwt2 = dzero

          RETURN
        endif
    endif
endif !fire_switchf

! read data 
!...check nsib
CHECK(nf90_inq_dimid(fireid, trim(fnsibname), dimid))
CHECK(nf90_inquire_dimension(fireid, dimid, dim_name, dimlen))
if (dimlen /= nsib) then
   print*,'Mismatching Fire File!'
   print*,' File nsib: ',dimlen,' Sim nsib: ',nsib
   STOP
endif

!...check time values
ENSURE_VAR(fireid,'year',ncyid)
ENSURE_VAR(fireid,'month',ncmid)
ENSURE_VAR(fireid,'day',ncdid)
mstart(1) = fire_recnum

!CHECK(nf90_get_var(fireid, ncyid, xyear, mstart(1:1)))
!CHECK(nf90_get_var(fireid, ncmid, xmonth, mstart(1:1)))
!CHECK(nf90_get_var(fireid, ncdid, xday, mstart(1:1)))
!print*,'firegfed_read_global mid fire_year:',fire_year
!print*,'firegfed_read_global mid file xyear:',xyear
!print*,'firegfed_read_global mid fire_month:',fire_month
!print*,'firegfed_read_global mid file xmonth:',xmonth
!print*,'firegfed_read_global mid fire_day:',fire_day
!print*,'firegfed_read_global mid file xday:',xday

CHECK(nf90_get_var(fireid, ncyid, xyear, mstart(1:1)))
if (xyear .ne. fire_year) then
   print*,'Fire year in file does not match simulation.  Stopping.'
   print*,'   File: ',xyear,' Sim: ',year
   stop
endif

CHECK(nf90_get_var(fireid, ncmid, xmonth, mstart(1:1)))
!CHECK(nf90_get_var(fireid, ncdid, xday, mstart(1:1)))
if (xmonth .ne. fire_month) then
   print*,'Fire month in file does not match simulation.  Stopping.'
   print*,'   File: ',xmonth,' Sim: ',month, '   fire_month:',fire_month
!   print*,'   File: ',xyear,' Sim: ',year, 'fire_year:',fire_year
!   print*,'fire_recnum:',fire_recnum
!   print*,'fireid: ',fireid
!   print*,'mstart(1:1) :',mstart(1:1)
!   print*,'  File: ',xday,' Sim: ', day
   stop
endif

CHECK(nf90_get_var(fireid, ncdid, xday, mstart(1:1)))
if (xday .ne. fire_day) then
   print*,'Fire day in file does not match simulation. Stopping.'
   print*,'  File: ',xday,' Sim: ', day
!   print*,'   File: ',xyear,' Sim: ',year
   stop
endif

!CHECK(nf90_get_var(fireid, ncdid, xday, mstart(1:1)))
!if (xday .ne. fire_day) then
!   print*,'Fire day in file does not match simulation. Stopping.'
!   print*,'  File: ',xday,' Sim: ', day
!   stop
!endif

!...get variable id's
! corresponds to burned area for each PFT
ENSURE_VAR( fireid, trim(ba_dbgname), fba_dbgid) 
ENSURE_VAR( fireid, trim(ba_en2name), fba_en2id) 
ENSURE_VAR( fireid, trim(ba_en3name), fba_en3id)
ENSURE_VAR( fireid, trim(ba_dnfname), fba_dnfid)
ENSURE_VAR( fireid, trim(ba_eb1name), fba_eb1id)
ENSURE_VAR( fireid, trim(ba_eb2name), fba_eb2id)
ENSURE_VAR( fireid, trim(ba_db1name), fba_db1id)
ENSURE_VAR( fireid, trim(ba_db2name), fba_db2id)
ENSURE_VAR( fireid, trim(ba_db3name), fba_db3id)
ENSURE_VAR( fireid, trim(ba_shbname), fba_shbid)
ENSURE_VAR( fireid, trim(ba_shaname), fba_shaid)
ENSURE_VAR( fireid, trim(ba_c3aname), fba_c3aid)
ENSURE_VAR( fireid, trim(ba_c3gname), fba_c3gid)
ENSURE_VAR( fireid, trim(ba_c4gname), fba_c4gid)
ENSURE_VAR( fireid, trim(ba_c3cname), fba_c3cid)
ENSURE_VAR( fireid, trim(ba_c4cname), fba_c4cid)
ENSURE_VAR( fireid, trim(ba_mzename), fba_mzeid)
ENSURE_VAR( fireid, trim(ba_soyname), fba_soyid)
ENSURE_VAR( fireid, trim(ba_wwtname), fba_wwtid)

!...get data
mstart=(/1,fire_recnum/); mcount=(/nsib,1/)
STATUS = nf90_get_var(fireid, fba_dbgid, ba_dbgin, mstart, mcount)
IF (status .ne. nf90_noerr) THEN
   print*,'Error Getting DBG Burned Area!'
   stop
ENDIF

status = nf90_get_var(fireid, fba_en2id, ba_en2in, mstart, mcount)
IF (status .ne. nf90_noerr) THEN
   print*,'Error Getting EN2 Burned Area!'
   stop
ENDIF

status = nf90_get_var(fireid, fba_en3id, ba_en3in, mstart, mcount)
IF (status .ne. nf90_noerr) THEN
   print*,'Error Getting EN3 Burned Area!'
   stop
ENDIF

status = nf90_get_var(fireid, fba_dnfid, ba_dnfin, mstart, mcount)
IF (status .ne. nf90_noerr) THEN
   print*,'Error Getting DNF Burned Area!'
   stop
ENDIF

status = nf90_get_var(fireid, fba_eb1id, ba_eb1in, mstart, mcount)
IF (status .ne. nf90_noerr) THEN
   print*,'Error Getting EB1 Burned Area!'
   stop
ENDIF

status = nf90_get_var(fireid, fba_eb2id, ba_eb2in, mstart, mcount)
IF (status .ne. nf90_noerr) THEN
   print*,'Error Getting EB2 Burned Area!'
   stop
ENDIF

status = nf90_get_var(fireid, fba_db1id, ba_db1in, mstart, mcount)
IF (status .ne. nf90_noerr) THEN
   print*,'Error Getting DB1 Burned Area!'
   stop
ENDIF

status = nf90_get_var(fireid, fba_db2id, ba_db2in, mstart, mcount)
IF (status .ne. nf90_noerr) THEN
   print*,'Error Getting DB2 Burned Area!'
   stop
ENDIF

status = nf90_get_var(fireid, fba_db3id, ba_db3in, mstart, mcount)
IF (status .ne. nf90_noerr) THEN
   print*,'Error Getting DB3 Burned Area!'
   stop
ENDIF

status = nf90_get_var(fireid, fba_shbid, ba_shbin, mstart, mcount)
IF (status .ne. nf90_noerr) THEN
   print*,'Error Getting SHB Burned Area!'
   stop
ENDIF

status = nf90_get_var(fireid, fba_shaid, ba_shain, mstart, mcount)
IF (status .ne. nf90_noerr) THEN
   print*,'Error Getting SHA Burned Area!'
   stop
ENDIF

status = nf90_get_var(fireid, fba_c3aid, ba_c3ain, mstart, mcount)
IF (status .ne. nf90_noerr) THEN
   print*,'Error Getting C3A Burned Area!'
   stop
ENDIF

status = nf90_get_var(fireid, fba_c3gid, ba_c3gin, mstart, mcount)
IF (status .ne. nf90_noerr) THEN
   print*,'Error Getting C3G Burned Area!'
   stop
ENDIF

status = nf90_get_var(fireid, fba_c4gid, ba_c4gin, mstart, mcount)
IF (status .ne. nf90_noerr) THEN
   print*,'Error Getting C4G Burned Area!'
   stop
ENDIF

status = nf90_get_var(fireid, fba_c3cid, ba_c3cin, mstart, mcount)
IF (status .ne. nf90_noerr) THEN
   print*,'Error Getting C3C Burned Area!'
   stop
ENDIF

status = nf90_get_var(fireid, fba_c4cid, ba_c4cin, mstart, mcount)
IF (status .ne. nf90_noerr) THEN
   print*,'Error Getting C4C Burned Area!'
   stop
ENDIF

status = nf90_get_var(fireid, fba_mzeid, ba_mzein, mstart, mcount)
IF (status .ne. nf90_noerr) THEN
   print*,'Error Getting MZE Burned Area!'
   stop
ENDIF

status = nf90_get_var(fireid, fba_soyid, ba_soyin, mstart, mcount)
IF (status .ne. nf90_noerr) THEN
   print*,'Error Getting SOY Burned Area!'
   stop
ENDIF

status = nf90_get_var(fireid, fba_wwtid, ba_wwtin, mstart, mcount)
IF (status .ne. nf90_noerr) THEN
   print*,'Error Getting WWT Burned Area!'
   stop
ENDIF

!...pull out points in subdomain
do i=1, subcount
    sib%g(i)%gprogt%ba_dbg2 = ba_dbgin(subset(i))*(1./(3.*60.*60.))
    sib%g(i)%gprogt%ba_en22 = ba_en2in(subset(i))*(1./(3.*60.*60.))
    sib%g(i)%gprogt%ba_en32 = ba_en3in(subset(i))*(1./(3.*60.*60.))
    sib%g(i)%gprogt%ba_dnf2 = ba_dnfin(subset(i))*(1./(3.*60.*60.))
    sib%g(i)%gprogt%ba_eb12 = ba_eb1in(subset(i))*(1./(3.*60.*60.))
    sib%g(i)%gprogt%ba_eb22 = ba_eb2in(subset(i))*(1./(3.*60.*60.))
    sib%g(i)%gprogt%ba_db12 = ba_db1in(subset(i))*(1./(3.*60.*60.))
    sib%g(i)%gprogt%ba_db22 = ba_db2in(subset(i))*(1./(3.*60.*60.))
    sib%g(i)%gprogt%ba_db32 = ba_db3in(subset(i))*(1./(3.*60.*60.))
    sib%g(i)%gprogt%ba_shb2 = ba_shbin(subset(i))*(1./(3.*60.*60.))
    sib%g(i)%gprogt%ba_sha2 = ba_shain(subset(i))*(1./(3.*60.*60.))
    sib%g(i)%gprogt%ba_c3a2 = ba_c3ain(subset(i))*(1./(3.*60.*60.))
    sib%g(i)%gprogt%ba_c3g2 = ba_c3gin(subset(i))*(1./(3.*60.*60.))
    sib%g(i)%gprogt%ba_c4g2 = ba_c4gin(subset(i))*(1./(3.*60.*60.))
    sib%g(i)%gprogt%ba_c3c2 = ba_c3cin(subset(i))*(1./(3.*60.*60.))
    sib%g(i)%gprogt%ba_c4c2 = ba_c4cin(subset(i))*(1./(3.*60.*60.))
    sib%g(i)%gprogt%ba_mze2 = ba_mzein(subset(i))*(1./(3.*60.*60.))
    sib%g(i)%gprogt%ba_soy2 = ba_soyin(subset(i))*(1./(3.*60.*60.))
    sib%g(i)%gprogt%ba_wwt2 = ba_wwtin(subset(i))*(1./(3.*60.*60.))

    !sib%g(i)%gprogt%ba_dbg2 = ba_dbgin(subset(i))
    !sib%g(i)%gprogt%ba_en22 = ba_en2in(subset(i))
    !sib%g(i)%gprogt%ba_en32 = ba_en3in(subset(i))
    !sib%g(i)%gprogt%ba_dnf2 = ba_dnfin(subset(i))
    !sib%g(i)%gprogt%ba_eb12 = ba_eb1in(subset(i))
    !sib%g(i)%gprogt%ba_eb22 = ba_eb2in(subset(i))
    !sib%g(i)%gprogt%ba_db12 = ba_db1in(subset(i))
    !sib%g(i)%gprogt%ba_db22 = ba_db2in(subset(i))
    !sib%g(i)%gprogt%ba_db32 = ba_db3in(subset(i))
    !sib%g(i)%gprogt%ba_shb2 = ba_shbin(subset(i))
    !sib%g(i)%gprogt%ba_sha2 = ba_shain(subset(i))
    !sib%g(i)%gprogt%ba_c3a2 = ba_c3ain(subset(i))
    !sib%g(i)%gprogt%ba_c3g2 = ba_c3gin(subset(i))
    !sib%g(i)%gprogt%ba_c4g2 = ba_c4gin(subset(i))
    !sib%g(i)%gprogt%ba_c3c2 = ba_c3cin(subset(i))
    !sib%g(i)%gprogt%ba_c4c2 = ba_c4cin(subset(i))
    !sib%g(i)%gprogt%ba_mze2 = ba_mzein(subset(i))
    !sib%g(i)%gprogt%ba_soy2 = ba_soyin(subset(i))
    !sib%g(i)%gprogt%ba_wwt2 = ba_wwtin(subset(i))
enddo

!...print out the new data if requested
if (print_fire) then
    print*,'     ------------------------------------'
    print*,'     !!!READING FIRE EMISSIONS!!!'

    if (firegfed_switchf) then
       print*,'     Opening fire file: '
       print*,'       ',trim(firegfed_filename)
    endif
    print*,'   REAL TIMES:' 
    print*,'   Year   Mon   Day   Hour'
    print('(a,i5,a,i4,a,i4,a,i4)'), '   ', &
         year,'  ',month,'  ',day,'  ',hour
    print*,''
    print*,'   FIRE TIMES:'
    print*,'   Year    Mon   Day   Hour'
    print('(a,i5,a,i4,a,i4,a,f6.2)'), '   ', &
          fire_year,'  ',fire_month,'  ',fire_day, &
          '  ',fire_hour
    print*,'   FIRE RECNUM/TOTNUM: ', fire_recnum, fire_permon
    !print('(a,2e16.6)'),'     Min/Max Fire C Loss (umol C/m2/s)    : ', &
    !      minval(firecin), maxval(firecin)
    !print('(a,2e16.6)'),'     Min/Max Fire CO2 Resp (umol CO2/m2/s): ', &
    !      minval(fireco2in), maxval(fireco2in)
    print('(a,2e16.6)'),'   Min/Max  BA DBG BA (m2 PFT⁻¹ day⁻¹)    : ', &
              minval(ba_dbgin), maxval(ba_dbgin)
    print('(a,2e16.6)'),'   Min/Max  BA EN2 BA (m2 PFT⁻¹ day⁻¹)    : ', &
              minval(ba_en2in), maxval(ba_en2in)
    print('(a,2e16.6)'),'   Min/Max  BA EN3 BA (m2 PFT⁻¹ day⁻¹)    : ', &
              minval(ba_en3in), maxval(ba_en3in)
    print('(a,2e16.6)'),'   Min/Max  BA DNF BA (m2 PFT⁻¹ day⁻¹)    : ', &
              minval(ba_dnfin), maxval(ba_dnfin)
    print('(a,2e16.6)'),'   Min/Max  BA EB1 BA (m2 PFT⁻¹ day⁻¹)    : ', &
              minval(ba_eb1in), maxval(ba_eb1in)
    print('(a,2e16.6)'),'   Min/Max  BA EB2 BA (m2 PFT⁻¹ day⁻¹)    : ', &
              minval(ba_eb2in), maxval(ba_eb2in)
    print('(a,2e16.6)'),'   Min/Max  BA DB1 BA (m2 PFT⁻¹ day⁻¹)    : ', &
              minval(ba_db1in), maxval(ba_db1in)
    print('(a,2e16.6)'),'   Min/Max  BA DB2 BA (m2 PFT⁻¹ day⁻¹)    : ', &
              minval(ba_db2in), maxval(ba_db2in)
    print('(a,2e16.6)'),'   Min/Max  BA DB3 BA (m2 PFT⁻¹ day⁻¹)    : ', &
              minval(ba_db3in), maxval(ba_db3in)
    print('(a,2e16.6)'),'   Min/Max  BA SHB BA (m2 PFT⁻¹ day⁻¹)    : ', &
              minval(ba_shbin), maxval(ba_shbin)
    print('(a,2e16.6)'),'   Min/Max  BA SHA BA (m2 PFT⁻¹ day⁻¹)    : ', &
              minval(ba_shain), maxval(ba_shain)
    print('(a,2e16.6)'),'   Min/Max  BA C3A BA (m2 PFT⁻¹ day⁻¹)    : ', &
              minval(ba_c3ain), maxval(ba_c3ain)
    print('(a,2e16.6)'),'   Min/Max  BA C3G BA (m2 PFT⁻¹ day⁻¹)    : ', &
              minval(ba_c3gin), maxval(ba_c3gin)
    print('(a,2e16.6)'),'   Min/Max  BA C4G BA (m2 PFT⁻¹ day⁻¹)    : ', &
              minval(ba_c4gin), maxval(ba_c4gin)
    print('(a,2e16.6)'),'   Min/Max  BA C3C BA (m2 PFT⁻¹ day⁻¹)    : ', &
              minval(ba_c3cin), maxval(ba_c3cin)
    print('(a,2e16.6)'),'   Min/Max  BA C4C BA (m2 PFT⁻¹ day⁻¹)    : ', &
              minval(ba_c4cin), maxval(ba_c4cin)
    print('(a,2e16.6)'),'   Min/Max  BA MZE BA (m2 PFT⁻¹ day⁻¹)    : ', &
              minval(ba_mzein), maxval(ba_mzein)
    print('(a,2e16.6)'),'   Min/Max  BA SOY BA (m2 PFT⁻¹ day⁻¹)    : ', &
              minval(ba_soyin), maxval(ba_soyin)
    print('(a,2e16.6)'),'   Min/Max  BA WWT BA (m2 PFT⁻¹ day⁻¹)    : ', &
              minval(ba_wwtin), maxval(ba_wwtin)

    print*,'     ------------------------------------'
 
    if (print_stop) stop
endif

end subroutine firegfed_read_global
