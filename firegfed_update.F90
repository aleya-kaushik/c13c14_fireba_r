!-------------------------------------------------------------------------------
subroutine firegfed_update()
!-------------------------------------------------------------------------------

!Updates the time/control variables for the fire data.

!-------------------------------------------------------------------------------

use kinds
use module_pparams, only:  &
    secs_per_day, secs_per_hr, &
    hrs_per_day
use module_io
use module_sib, only: sib
use module_sibconst, only: single_pt, subcount
use module_time

implicit none

! local
integer(long) :: local_secahead

!-----------------------------------------------------

!print*,'firegfed_update top fire_year:',fire_year
!print*,'firegfed_update top fire_month:',fire_month
!print*,'firegfed_update top fire_day:',fire_day

if (firegfed_readf) then
    !...read in fire data needed
    if (single_pt) then
       call firegfed_read_single(sib%g(1)%gprogt)
    else
       call firegfed_read_global()
    endif

    fire_recnum = fire_recnum + 1
    fire_seccur = fire_seccur + fire_step
    fire_secnext = fire_secnext + fire_step
    local_secahead = fire_secnext + fire_step
    fire_hour = fire_hour + real(fire_step/3600.)

    firegfed_readf = .false.
    firegfed_switchf = .false.
    if (fire_recnum .gt. fire_permon) then
       !...switch fire data for new month
       if (local_secahead .lt. end_second) then
          fire_recnum = 1
          fire_permon = curday_per_mon_day_ahead * fire_perday
          firegfed_switchf = .true.
       else
          firegfed_updatef = .false.
       endif
    endif
    
    if (fire_hour >= hrs_per_day) then
        fire_day = fire_day + 1
        fire_hour = fire_hour - hrs_per_day
     endif
    if (fire_day > curday_per_mon) then
       fire_day = 1
       fire_month = fire_month + 1
    endif
    if (fire_month > 12) then
       fire_month = 1
       fire_year = fire_year+1
    endif
    if (fire_year .gt. endyear) fire_year=startyear
    
endif !firegfed_readf

!print*,'firegfed_update end fire_year:',fire_year
!print*,'firegfed_update end fire_month:',fire_month
!print*,'firegfed_update end fire_day:',fire_day
 
end subroutine firegfed_update
