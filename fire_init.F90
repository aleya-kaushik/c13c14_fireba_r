!---------------------------------------------------------------------
subroutine fire_init( )
!---------------------------------------------------------------------

!---------------------------------------------------------------------
! Sets initial fire information
!
! Expects standard driver output:
!   Year, DOY, Hour
!   Emis_C (micromoles C/m2/s)
!   Emis_CO2 (micromoles C/m2/s)
!
! For time specification in the file:
!   Mid-point of valid time range 
!
! For format:
!   Global: Netcdf File
!   Site: Text File
!
!---------------------------------------------------------------------

use kinds
use module_io
use module_pparams, only: &
    days_per_month, secs_per_day, mol_to_umol
use module_sib, only: sib
use module_sibconst
use module_time

use netcdf
implicit none

!parameters
logical, parameter :: fire_required = .true.

!fire data step setting variables
integer(i4) :: status, varid
integer(i4) :: fyr, fdoy, fmon, fmonref, dpermon
integer(i4) :: fire_recnumstart
integer(i4) :: restart_ncid, restart_recnum, restart_year, &
    restart_month, restart_switch
integer(long) :: fsec_tot
real(r4), dimension(2) :: fhr
character(len=256) :: record
logical :: exist

!local variables
integer :: frnum
character(len=iofnlen) :: fire_tempname

!----------------------------------------------------------------------

if ((spinup) .or. (.not. fire_switch)) then
   !Only use fires for production runs
   fire_step = izero
   RETURN
else
   print*,'Setting Fire Data'
endif

!...find first fire file
fire_year = izero
fire_month = izero
do fyr=init_year, end_year
   if (fyr .eq. init_year) then
      fmonref = init_month
      if (fmonref .gt. 1) then
          dpermon = days_per_month(fmonref-1)
          if ( (leapyr_switch) .and. &
               (fmonref .eq. 2) .and. &
               (((mod(fyr,4) .eq. 0) .and. &
               (mod(fyr,100) .ne. 0)) .or. &
               (mod(fyr,400) .eq. 0)) ) then
              dpermon = dpermon + 1
           endif
      else
          dpermon = izero
      endif
      fsec_tot = dpermon*secs_per_day
   else
      fmonref = 1
   endif
         
   do fmon=fmonref,12
      if (single_pt) then
          write(fire_tempname,fmt='(a,i4.4,i2.2,a4)') &
               trim(fire_path), fyr, fmon, '.dat'
      else
         write(fire_tempname,fmt='(a,i4.4,i2.2,a3)') &
              trim(fire_path), fyr, fmon, '.nc'
      endif
      inquire(file=trim(fire_tempname), exist=exist)

      dpermon = days_per_month(fmon)
      if ( (leapyr_switch) .and. &
           (fmon .eq. 2) .and. &
           (((mod(fyr,4) .eq. 0) .and. &
           (mod(fyr,100) .ne. 0)) .or. &
           (mod(fyr,400) .eq. 0)) ) then
         dpermon = dpermon + 1
      endif
      
      if ((exist) .AND. (fire_year .eq. izero)) then
         fire_year = fyr
         fire_month = fmon
         fire_seccur = fsec_tot
         fire_filename = fire_tempname
      else
         fsec_tot = fsec_tot + dpermon*secs_per_day
      endif
   enddo
enddo


IF (fire_year .eq. izero) THEN
    if (fire_required) THEN
        print*,''
        print*,'Missing/Non-Existant Fire Emissions File!'
        print*,'Directory: ',trim(fire_path)
        print*,'Please check fr_path in namel_sibdrv.'
        print*,''
        STOP
     else
        print*,'  Missing Fire Data, Using Zero Emissions.'
        fire_step = 0
        RETURN
     endif
ENDIF


!...set fire timestep
if (single_pt) then !site simulation
    fireid = 43
    open(unit=fireid, file=trim(fire_filename), &
         status='old', form='formatted', iostat=status)
    DO   ! read until not a comment
         read(fireid,'(a)',iostat=status) record
         IF (status .ne. 0) THEN
            print*,'Stopping due to fire data error!'
            stop
         ENDIF
         IF (record(1:1) .ne. '#') THEN
            exit
         ENDIF
     ENDDO

     read(unit=record,fmt=*) fyr, fdoy, fhr(1)
     read(fireid,'(a)',iostat=status) record
     read(unit=record,fmt=*) fyr, fdoy, fhr(2)
     close(fireid)

else !global simulation
     status = nf90_open(trim(fire_filename), nf90_nowrite, &
                  fireid)
     status = nf90_inq_varid(fireid,'hour',varid)
     status = nf90_get_var(fireid,varid,fhr) 
     IF (status .ne. nf90_noerr) THEN
         print*,'Error Reading Hour in Global Fire Data!'
         STOP
     ENDIF
     status = nf90_close(fireid)
endif


!...set and check fire step
fire_step = int((fhr(2)-fhr(1))*3600)
if (fire_step .le. 0) then
   print*,'Missing Fire Time Step'
   STOP
elseif (mod(secs_per_day,fire_step) /= 0) then
   print*,'Fire Time Step Does Not Divide Evenly Into Day.'
   STOP
endif


!...set fire variables
fire_perday = secs_per_day/fire_step
fire_permon = curday_per_mon * fire_perday
fire_hour = fhr(1)

if ((fire_year .eq. init_year) .and. &
    (fire_month .eq. init_month)) then
     fire_recnumstart = (day-1)*fire_perday + 1
     fire_day = day
else
     fire_recnumstart = 1
     fire_day = 1
endif
! For an exact transient restart, the saved fire_recnum names the NEXT
! unconsumed record.  Position the sequential site file immediately
! before that record.  In particular, do not consume it during fire_init:
! restart_read restores the checkpoint bracket and asks fire_update to
! read that same record when appropriate.
if (ic_restart .and. single_pt) then
    status = nf90_open(trim(ic_file), nf90_nowrite, restart_ncid)
    if (status .ne. nf90_noerr) then
        print*, 'Cannot open restart checkpoint for fire positioning: ', &
             trim(nf90_strerror(status))
        stop
    endif
    status = nf90_get_att(restart_ncid, nf90_global, &
                         'fire_recnum', restart_recnum)
    if (status .ne. nf90_noerr) stop 'Missing restart fire_recnum'
    status = nf90_get_att(restart_ncid, nf90_global, &
                         'fire_year', restart_year)
    if (status .ne. nf90_noerr) stop 'Missing restart fire_year'
    status = nf90_get_att(restart_ncid, nf90_global, &
                         'fire_month', restart_month)
    if (status .ne. nf90_noerr) stop 'Missing restart fire_month'
    status = nf90_get_att(restart_ncid, nf90_global, &
                         'fire_switchf', restart_switch)
    if (status .ne. nf90_noerr) stop 'Missing restart fire_switchf'
    status = nf90_close(restart_ncid)
    if (status .ne. nf90_noerr) stop 'Error closing fire restart checkpoint'

    if (restart_switch .eq. 0) then
        if ((restart_year .ne. fire_year) .or. &
            (restart_month .ne. fire_month)) then
            print*, 'Fire restart month differs from initial fire file: ', &
                 restart_year, restart_month, fire_year, fire_month
            stop
        endif
        if ((restart_recnum .lt. 1) .or. &
            (restart_recnum .gt. fire_permon)) then
            print*, 'Invalid next fire record in restart: ', &
                 restart_recnum, fire_permon
            stop
        endif
        fire_recnumstart = restart_recnum - 1
        print*, 'Positioning site fire file before next restart record: ', &
             restart_recnum
    endif
endif

fire_secnext = fire_seccur + fire_step
fire_updatef = .true.
fire_readf = .false.


!...read in initial fire data
fire_switchf = .true.
fire_recnum = 1
if ((single_pt) .and. (fire_recnumstart .eq. 0)) then
    ! A checkpoint may be positioned just before the first monthly record.
    ! Open the file without consuming its first data record.
    open(unit=fireid, file=trim(fire_filename), status='old', &
         form='formatted', iostat=status)
    if (status .ne. 0) stop 'Cannot open initial fire file for restart'
endif

do frnum=1, fire_recnumstart
    if (single_pt) then
       call fire_read_single(sib%g(1)%gprogt)
    else
       call fire_read_global()
    endif
    ! The first read opens the monthly file; later reads must advance within it.
    fire_switchf = .false.
    fire_hour = fire_hour + real(fire_step/3600.)
    fire_recnum = fire_recnum + 1
    fire_seccur = fire_seccur + fire_step
enddo
fire_switchf = .false.

sib%g(:)%gprogt%firec1 = sib%g(:)%gprogt%firec2
sib%g(:)%gprogt%firec  = sib%g(:)%gprogt%firec2
sib%g(:)%gprogt%fireco21 = sib%g(:)%gprogt%fireco22
sib%g(:)%gprogt%fireco2  = sib%g(:)%gprogt%fireco22

end subroutine fire_init


