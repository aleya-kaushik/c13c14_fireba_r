#include "nc_util.h"

!--------------------------------------------
!Includes routines to create new output files
!   for both global and land unit output.
!---------------------------------------------

subroutine diagnostic_createg( &
    year, curday_per_year,     &
    filename, numvars, subcount, &
    lonsib, latsib, sitesib, &
    period_start, period_step, &
    nsaveperfile, &
    namevars, listvars, &
    fileid, varids, timeid, append_count)

    use module_io, only: iofnlen, ic_restart
    use module_sibconst, only: slen
    use netcdf
    use kinds

    implicit none

    ! input variables
    integer(i4), intent(in) :: year            !current year of simulation
    integer(i4), intent(in) :: curday_per_year !number of days in current year
    character(len=iofnlen), intent(in) :: filename   ! filename for output
    integer(i4), intent(in) :: numvars           ! max # of variables written to file
    integer(i4), intent(in) :: subcount          ! number of landpoints on subset
    real(r4), dimension(subcount), intent(in) :: lonsib   ! array of longitude coordinates
    real(r4), dimension(subcount), intent(in) :: latsib   ! array of latitude coordinates
    character(len=slen), dimension(subcount), intent(in) :: sitesib  ! site names
    real(r8), intent(in) :: period_start    ! start of averaged period
    real(r8), intent(in) :: period_step     ! timestep for averaged period
    integer(i4), intent(in) :: nsaveperfile ! number of times saved per file
    character(len=21), dimension(numvars), intent(in) :: namevars   ! names of vars
    character(len=100), dimension(numvars), intent(in) :: listvars  ! variable descriptions
    integer(i4), intent(out) :: fileid           ! file id#
    integer(i4), dimension(numvars), intent(out):: varids  ! variable id#s 
    integer(i4), intent(out) :: timeid           ! time variable id#
    integer(i4), intent(in) :: append_count ! verified prior PBP records

    ! netcdf variables
    integer(i4) :: nscid   ! landpoints dimension id #
    integer(i4) :: slenid  ! site name string length dimension id #
    integer(i4) :: timedid ! time dimension id #

    integer(i4) :: latsibid ! latitude variable id #
    integer(i4) :: lonsibid ! longitude variable id #
    integer(i4) :: snameid  ! site name variable id #

    ! local variables
    logical :: savesite
    integer(i4) :: t, nv, syeari
    integer(i4) :: status          ! netcdf status
    integer(i4) :: old_count, nrecord, dimid, old_vid
    real(r8) :: old_time, expected_time
    logical :: already_exists
    character(len=4)  :: syear     !year string
    character(len=40) :: units          ! variable units
    character(len=80) :: longname       ! variable description
    character(len=60) :: titlename      ! variable title

    integer(i4) :: step
    real(r8) :: ddoy, dtime

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! determine if saving site names
    if (sitesib(1) .ne. '') then
       savesite = .true.
    else
       savesite = .false.
    endif

    ! Append ONLY to an exact, validated partial-month PBP file.
    old_count = 0
    old_count = append_count
    if (old_count > 0) then
       call diagnostic_open_existing(filename,old_count,nsaveperfile, &
            period_start,period_step,numvars,namevars,fileid,varids,timeid, &
            subcount, 0)
       return
    endif
    ! make sure fileid is not tied to any open file
    status = nf90_close( fileid )
    if (ic_restart) then
       inquire(file=trim(filename),exist=already_exists)
       if (already_exists) &
          error stop 'Diagnostic continuation: refusing to clobber existing output file'
    endif
    CHECK( nf90_create( trim(filename), nf90_clobber, fileid) )

    ! define dimensions
    CHECK( nf90_def_dim( fileid, 'landpoints', subcount, nscid ) )
    CHECK( nf90_def_dim( fileid, 'time', nf90_unlimited, timedid ) )
    IF (savesite) CHECK( nf90_def_dim( fileid, 'slen', slen, slenid ) )    

    ! define control variables
    CHECK( nf90_def_var( fileid, 'lonsib', nf90_float, (/nscid/), lonsibid ) )
    CHECK( nf90_put_att( fileid, lonsibid, 'long_name', 'SiB Point Longitude' ) )
    CHECK( nf90_put_att( fileid, lonsibid, 'units', 'degrees_east' ) )
    CHECK( nf90_put_att( fileid, lonsibid, 'quantity', 'longitude' ) )

    CHECK( nf90_def_var( fileid, 'latsib', nf90_float, (/nscid/), latsibid ) )
    CHECK( nf90_put_att( fileid, latsibid, 'long_name', 'SiB Point Latitude' ) )
    CHECK( nf90_put_att( fileid, latsibid, 'units', 'degrees_north' ) )
    CHECK( nf90_put_att( fileid, latsibid, 'quantity', 'latitude' ) )

    IF (sitesib(1) .NE. '') THEN
       CHECK( nf90_def_var( fileid, 'site_names', nf90_char, (/slenid,nscid/), snameid ) )
    ENDIF

    CHECK( nf90_def_var( fileid, 'time', nf90_double, (/timedid/), timeid ) )
    CHECK( nf90_put_att( fileid, timeid, 'quantity', 'time' ) )
    syeari = index(trim(filename),'_',back=.true.)
    syeari = syeari + 1
    write(syear,'(a4)') filename(syeari:syeari+3)
    CHECK( nf90_put_att( fileid, timeid, 'units', 'days since '//syear//'-01-01'  ) )
    
    ! define data variables
    do nv = 1, numvars
       status = nf90_def_var(fileid,trim(namevars(nv)),nf90_double,(/nscid,timedid/),varids(nv))
       if (status .ne. nf90_noerr) then
          print*,'Error defining variables for diagnostic_createg.F90. Stopping.'
          stop
       endif
       call get_units( listvars(nv), longname, units )
       call get_title( longname, titlename )

       CHECK( nf90_put_att( fileid, varids(nv), 'long_name', trim(longname) ) )
       CHECK( nf90_put_att( fileid, varids(nv), 'title', trim(titlename) ) )
       CHECK( nf90_put_att( fileid, varids(nv), 'units', trim(units) ) )
    enddo

    ! switch from definition mode to data mode
    CHECK( nf90_enddef( fileid ) )

    ! assign values to variables not variant with time
    CHECK( nf90_put_var( fileid, lonsibid, lonsib ) )
    CHECK( nf90_put_var( fileid, latsibid, latsib ) )
    IF (savesite) CHECK( nf90_put_var( fileid, snameid, sitesib) )

    ! write time variables
    ddoy = dble(period_start) + (period_step * 0.5)
    step=0
    do t=1,nsaveperfile
       dtime = year + ddoy / dble(curday_per_year)
       step = step + 1
       CHECK( nf90_put_var( fileid, timeid, ddoy, (/step/) ) )
       ddoy = ddoy + period_step
   enddo
    
end subroutine diagnostic_createg

!------------------------------------------------------------

subroutine diagnostic_createlu( &
    year, cdaypyr, &
    filename, numvars, subcount, &
    lonsib, latsib, sitesib, &
    lgarea, lgref, &
    period_start, period_step, &
    nsaveperfile,  &
    namevars, listvars, & 
    fileid, varids, timeid, append_count)

    use netcdf
    use kinds
    use module_io, only: iofnlen, ic_restart
    use module_pftinfo, only: &
         clen, pft_name, pft_ref
    use module_sibconst, only: &
         nlu, npft, slen

    implicit none

    ! input variables
    integer(i4), intent(in) :: year     !current year of simulation
    integer(i4), intent(in) :: cdaypyr  !number of days in current year
    character(len=iofnlen), intent(in) :: filename
    integer(i4), intent(in) :: numvars     ! # of variables to write out
    integer(i4), intent(in) :: subcount    ! # of landpoints simulated
    real(r4), dimension(subcount), intent(in) :: lonsib   ! array of longitude coordinates
    real(r4), dimension(subcount), intent(in) :: latsib   ! array of latitude coordinates
    character(len=slen), dimension(subcount), intent(in) :: sitesib ! site names of landpoints
    real(r4), dimension(subcount,nlu), intent(in) :: lgarea    ! % areal coverage of lu/gridcell
    integer(i4), dimension(subcount,nlu), intent(in) :: lgref  !reference of lu/gridcell 
    real(r8), intent(in) :: period_start    ! start of averaged period
    real(r8), intent(in) :: period_step     ! timestep for averaged period
    integer(i4), intent(in) :: nsaveperfile ! number of times saved per file
    character(len=21), dimension(numvars), intent(in) ::  namevars  ! variable names
    character(len=100), dimension(numvars), intent(in) :: listvars ! variable descriptions
    integer(i4), intent(out) :: fileid     ! file id #
    integer(i4), dimension(numvars), intent(out) :: varids  ! variable id #s
    integer(i4), intent(out) :: timeid       ! time variable id #
    integer(i4), intent(in) :: append_count ! verified prior PBP/HR records

    ! netcdf variables
    integer(i4) :: nscid  ! landpoints dimension id #
    integer(i4) :: nluid  ! nlu dimension id #
    integer(i4) :: npftid ! npft dimension id #
    integer(i4) :: clenid ! clen dimension id #
    integer(i4) :: slenid ! slen dimension id #
    integer(i4) :: tdid   ! time dimension id #
    integer(i4) :: pnameid, prefid    ! pft name and reference id #
    integer(i4) :: lonsibid, latsibid ! landpoints lat/lon variable id #
    integer(i4) :: snameid ! sitename variable id #
    integer(i4) :: areaid  ! percent coverage of gridcell
    integer(i4) :: rupid   ! Static PFT reference number
    integer(i4) :: dynpftid, dynareaid ! Time-varying LU PFT and area

    ! local variables
    logical :: savesite   ! flag for writing out site names
    integer(i4) :: n, t, syeari      ! index variable
    integer(i4) :: status            ! status of netcdf
    integer(i4) :: old_count, nrecord, dimid, old_vid
    real(r8) :: old_time, expected_time
    logical :: already_exists 
    character(len=4)  :: syear       ! year string
    character(len=40) :: units       ! variable units
    character(len=80) :: longname    ! variable description
    character(len=60) :: titlename   ! variable title

    integer(i4) :: step
    real(r8) :: ddoy, dtime

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! determine if saving site names
    if (sitesib(1) .ne. '') then
       savesite = .true.
    else
       savesite = .false.
    endif

    ! Append ONLY to an exact, validated partial-month PBP file.
    old_count = 0
    old_count = append_count
    if (old_count > 0) then
       call diagnostic_open_existing(filename,old_count,nsaveperfile, &
            period_start,period_step,numvars,namevars,fileid,varids,timeid, &
            subcount,nlu)
       return
    endif
    ! make sure fileid is not tied to any open file
    status = nf90_close( fileid )
    if (ic_restart) then
       inquire(file=trim(filename),exist=already_exists)
       if (already_exists) &
          error stop 'Diagnostic continuation: refusing to clobber existing output file'
    endif
    CHECK( nf90_create( trim(filename), nf90_clobber, fileid) )

    ! define dimensions
    CHECK( nf90_def_dim( fileid, 'landpoints', subcount, nscid ) )
    CHECK( nf90_def_dim( fileid, 'nlu', nlu, nluid ) )
    CHECK( nf90_def_dim( fileid, 'npft', npft, npftid ) )    
    CHECK( nf90_def_dim( fileid, 'time', nf90_unlimited, tdid ) )
    CHECK( nf90_def_dim( fileid, 'clen', clen, clenid ) )
    IF (savesite) CHECK( nf90_def_dim( fileid, 'slen', slen, slenid ) )  
    
    ! define control variables
    CHECK( nf90_def_var( fileid, 'lonsib', nf90_float, (/nscid/), lonsibid ) )
    CHECK( nf90_put_att( fileid, lonsibid, 'long_name', 'SiB Point Longitude' ) )
    CHECK( nf90_put_att( fileid, lonsibid, 'units', 'degrees_east' ) )

    CHECK( nf90_def_var( fileid, 'latsib', nf90_float, (/nscid/), latsibid ) )
    CHECK( nf90_put_att( fileid, latsibid, 'long_name', 'SiB Point Latitude' ) )
    CHECK( nf90_put_att( fileid, latsibid, 'units', 'degrees_north' ) )

    IF (sitesib(1) .NE. '') THEN
       CHECK( nf90_def_var( fileid, 'site_names', nf90_char, (/slenid,nscid/), snameid ) )
    ENDIF

    CHECK( nf90_def_var( fileid, 'pft_names', nf90_char, (/clenid,npftid/), pnameid ) )
    CHECK( nf90_put_att( fileid, pnameid, 'long_name', 'PFT Name'))

    CHECK( nf90_def_var( fileid, 'pft_refnums', nf90_byte, (/npftid/), prefid ) )
    CHECK( nf90_put_att( fileid, prefid, 'long_name', 'PFT Reference Numbers' ) )
    
    CHECK( nf90_def_var( fileid, 'lu_area', nf90_float, (/nscid, nluid/), areaid ) )
    CHECK( nf90_put_att( fileid, areaid, 'long_name', 'Land Unit Areal Coverage' ) )
    CHECK ( nf90_put_att( fileid, areaid, 'units', 'Fraction (0-1)' ) )

    CHECK( nf90_def_var( fileid, 'lu_pref', nf90_int, (/nscid, nluid/), rupid ) )
    CHECK( nf90_put_att( fileid, rupid, 'long_name', 'Land Unit PFT Reference' ) )
    CHECK( nf90_put_att( fileid, rupid, 'units', 'Number') )
    CHECK( nf90_put_att( fileid, rupid, 'min_val', 1 ) )
    CHECK( nf90_put_att( fileid, rupid, 'max_val', maxval(pft_ref) ) )
    CHECK( nf90_put_att( fileid, rupid, 'note', &
         'Static initial/reference PFT at output initialization; use lu_pref_time for transitions' ) )
    CHECK( nf90_put_att( fileid, areaid, 'note', &
         'Static initial/reference area at output initialization; use lu_area_time for transitions' ) )

    ! These are SNAPSHOTS at the END of each output interval. Unlike
    ! flux diagnostics, PFT IDs and fractional areas must never be averaged.
    ! 0 identifies an unoccupied LU slot in lu_pref_time.
    CHECK( nf90_def_var( fileid, 'lu_pref_time', nf90_int, &
         (/nscid,nluid,tdid/), dynpftid ) )
    CHECK( nf90_put_att( fileid, dynpftid, 'long_name', &
         'Active land-unit PFT reference at end of output interval' ) )
    CHECK( nf90_put_att( fileid, dynpftid, 'units', 'Number' ) )
    CHECK( nf90_put_att( fileid, dynpftid, 'inactive_value', 0 ) )
    CHECK( nf90_def_var( fileid, 'lu_area_time', nf90_float, &
         (/nscid,nluid,tdid/), dynareaid ) )
    CHECK( nf90_put_att( fileid, dynareaid, 'long_name', &
         'Active land-unit area fraction at end of output interval' ) )
    CHECK( nf90_put_att( fileid, dynareaid, 'units', 'Fraction (0-1)' ) )
    
    CHECK( nf90_def_var( fileid, 'time', nf90_double, (/tdid/), timeid ) )
    CHECK( nf90_put_att( fileid, timeid, 'quantity', 'time' ) )
    syeari = index(trim(filename),'_',back=.true.)
    syeari = syeari + 1
    write(syear,'(a4)') filename(syeari:syeari+3)
    CHECK( nf90_put_att( fileid, timeid, 'units', 'days since '//syear//'-01-01'  ) )


    ! define data variables
    do n = 1, numvars
        status = nf90_def_var(fileid,trim(namevars(n)),nf90_double,(/nscid,nluid,tdid/),varids(n))
        if (status .ne. nf90_noerr) then
           print*,'Error defining variables for diagnostic_createlu.F90. Stopping.'
           stop
        endif
        call get_units( listvars(n), longname, units )
        call get_title( longname, titlename )

        CHECK( nf90_put_att( fileid, varids(n), 'long_name', trim(longname) ) )
        CHECK( nf90_put_att( fileid, varids(n), 'title', trim(titlename) ) )
        CHECK( nf90_put_att( fileid, varids(n), 'units', trim(units) ) )
    enddo

    ! switch from definition mode to data mode
    CHECK( nf90_enddef( fileid ) )

    ! assign values to variables not variant with time
    CHECK( nf90_put_var( fileid, latsibid, latsib ) )
    CHECK( nf90_put_var( fileid, lonsibid, lonsib ) )
    IF (savesite) CHECK( nf90_put_var( fileid, snameid, sitesib) )
    CHECK( nf90_put_var( fileid, pnameid, pft_name) )
    CHECK( nf90_put_var( fileid, prefid, pft_ref) )
    CHECK( nf90_put_var( fileid, areaid, lgarea ) )
    CHECK( nf90_put_var( fileid, rupid, lgref ) )

    ! write time variables
    ddoy = dble(period_start) + (period_step * 0.5)
    step=0
    do t=1,nsaveperfile
       dtime = year + ddoy / dble(cdaypyr)
       step = step + 1
       CHECK( nf90_put_var( fileid, timeid, ddoy, (/step/) ) )
       ddoy = ddoy + period_step
    enddo
    
end subroutine diagnostic_createlu

!-----------------------------------------------------------------------
subroutine diagnostic_close(ncid)

    use netcdf
    use kinds

    implicit none

    ! input variables
    integer(i4), intent(inout) :: ncid

    ! local variables
    integer(i4) :: status

    status = nf90_noerr
    IF (ncid > 0) status = nf90_close(ncid)

    IF (status /= nf90_noerr) print*,'Error Closing Diagnostic Output'
    ncid = 0

end subroutine diagnostic_close

!-----------------------------------------------------------------------
subroutine get_units(description, longname, units)

    !-----------------------------------------------------------------------
    !  Purpose:
    !   Extracts a string enclosed within parentheses. 
    !   Used for units which are contained in a general description string.
    !
    !   Returns the units string (units) and
    !   the description string with the units removed (longname).
    !-----------------------------------------------------------------------

    character(len=*), intent(in) :: description
    character(len=*), intent(out) :: units
    character(len=*), intent(out) :: longname

    integer :: n, start_paren, end_paren, paren_count

    paren_count = 0
    start_paren = len_trim(description)
    end_paren = len_trim(description)

    do n = len(description), 1, -1
        if (description(n:n)==")") then
            if (paren_count == 0) then
                end_paren = n
            endif
            paren_count = paren_count + 1
        else if (description(n:n) == "(") then
            paren_count = paren_count - 1
            if (paren_count == 0) then
                start_paren = n
                exit
            endif
        end if
    end do

    !   in case of confusion, clear units and return unaltered description
    !   note: start_paren > end_paren should not be possible, but just in case...
    !         start_paren = end_paren occurs when there are no units
    !         start_paren = end_paren-1 occurs when units are "()"
    !   FIXME: n==1 is too limiting - what if I wanted only units?
    if (n == 1 .or. start_paren >= end_paren) then   ! no units
        units = " "
        longname = trim(description)
    else if (start_paren == (end_paren-1)) then      ! "()" case
        units = " "
        longname = trim(description(:start_paren-1))// &
            description(end_paren+1:)
    else                                             ! normal units
        units = description(start_paren+1:end_paren-1)
        longname = trim(description(:start_paren-1))// &
            description(end_paren+1:)
    end if

end subroutine get_units

!-----------------------------------------------------------------------
subroutine get_title(description, titlename)

    !-----------------------------------------------------------------------
    !  Purpose:
    !   Extracts the string following a colon.
    !
    !   Returns the title string (titlename).
    !-----------------------------------------------------------------------

    character(len=*), intent(in) :: description
    character(len=*), intent(out) :: titlename

    integer :: n, colon_pos, end_pos

    colon_pos = 0
    end_pos = len_trim(description)

    do n = len(description), 1, -1
        if (description(n:n)==":") then
           colon_pos = n
        end if
    end do

    if (colon_pos == 0) then
        titlename = trim(description)
    else
        titlename = trim(description(colon_pos+2:end_pos))
    end if

end subroutine get_title

!-----------------------------------------------------------------------
subroutine diagnostic_open_existing(filename, old_count, total_count, &
                                   period_start, period_step, nvars, names, &
                                   fileid, varids, timeid, expected_npts, expected_nlu)
    use kinds
    use netcdf
    use module_io, only: iofnlen, ic_restart
    implicit none
    character(len=iofnlen), intent(in) :: filename
    integer(i4), intent(in) :: old_count, total_count, nvars
    integer(i4), intent(in) :: expected_npts, expected_nlu
    real(r8), intent(in) :: period_start, period_step
    character(len=21), intent(in) :: names(nvars)
    integer(i4), intent(out) :: fileid, timeid, varids(nvars)
    integer(i4) :: status, dimid, ntime, v, t, ndim
    integer(i4) :: land_dimid, nlu_dimid, npts, file_nlu
    integer(i4) :: all_vars, time_vars, vdims(NF90_MAX_VAR_DIMS)
    real(r8) :: prior, expected, tm
    if (old_count < 1 .or. old_count >= total_count) &
       error stop 'PBP append: invalid previous record count'
    status = nf90_open(trim(filename),nf90_write,fileid)
    if (status /= nf90_noerr) error stop 'PBP append: existing file cannot be opened'
    status = nf90_inq_dimid(fileid,'time',dimid)
    if (status /= nf90_noerr) error stop 'Diagnostic append: existing file has no time dimension'
    status = nf90_inquire_dimension(fileid,dimid,len=ntime)
    if (status /= nf90_noerr .or. ntime /= old_count) &
       error stop 'Diagnostic append: previous record count differs from restart day'
    status = nf90_inq_dimid(fileid,'landpoints',land_dimid)
    if (status /= nf90_noerr) error stop 'Diagnostic append: no landpoints dimension'
    status = nf90_inquire_dimension(fileid,land_dimid,len=npts)
    if (status /= nf90_noerr .or. npts /= expected_npts) &
       error stop 'Diagnostic append: landpoint dimensions disagree'
    if (expected_nlu > 0) then
       status = nf90_inq_dimid(fileid,'nlu',nlu_dimid)
       if (status /= nf90_noerr) error stop 'Diagnostic append: no nlu dimension'
       status = nf90_inquire_dimension(fileid,nlu_dimid,len=file_nlu)
       if (status /= nf90_noerr .or. file_nlu /= expected_nlu) &
          error stop 'Diagnostic append: LU dimensions disagree'
    endif
    status = nf90_inquire(fileid,nVariables=all_vars)
    if (status /= nf90_noerr) error stop 'Diagnostic append: cannot inspect variables'
    time_vars=0
    do v=1,all_vars
       status = nf90_inquire_variable(fileid,v,ndims=ndim,dimids=vdims)
       if (status /= nf90_noerr) error stop 'Diagnostic append: cannot inspect variable dimensions'
       if (any(vdims(1:ndim)==dimid)) time_vars=time_vars+1
    enddo
    if (expected_nlu > 0) then
       if (time_vars /= nvars+3) &
          error stop 'Diagnostic append: previous LU output selection differs'
    else
       if (time_vars /= nvars+1) &
          error stop 'Diagnostic append: previous grid output selection differs'
    endif
    status = nf90_inq_varid(fileid,'time',timeid)
    if (status /= nf90_noerr) error stop 'PBP append: no time coordinate'
    expected = period_start + 0.5d0*period_step
    do t=1,old_count
       status = nf90_get_var(fileid,timeid,prior,start=(/t/))
       if (status /= nf90_noerr) error stop 'Diagnostic append: cannot read existing time'
       if (prior /= expected) &
          error stop 'Diagnostic append: prior timestamps do not match'
       expected = expected + period_step
    enddo
    do v=1,nvars
       status = nf90_inq_varid(fileid,trim(names(v)),varids(v))
       if (status /= nf90_noerr) &
          error stop 'PBP append: expected diagnostic variable missing'
    enddo
    ! Do not pre-create future time records: a second restart within this
    ! same month must see only records that were actually completed.  The
    ! output driver writes the timestamp alongside each finished interval.
end subroutine diagnostic_open_existing

!-----------------------------------------------------------------------
subroutine diagnostic_write_time(ncid, timeid, period_start, period_step, &
                                 first_record, nrecords)
    use netcdf
    use kinds
    implicit none
    integer(i4), intent(in) :: ncid, timeid, first_record, nrecords
    real(r8), intent(in) :: period_start, period_step
    integer(i4) :: i, status
    real(r8) :: timestamp
    if (nrecords < 1 .or. first_record < 1) &
       error stop 'Diagnostic write time: invalid record range'
    ! Reproduce precisely the floating-point accumulation in the original
    ! diagnostic_create writer (not a direct multiplication of record index).
    timestamp = period_start + 0.5d0*period_step
    do i=1,first_record+nrecords-1
       if (i >= first_record) then
          status = nf90_put_var(ncid,timeid,timestamp,start=(/i/))
          if (status /= nf90_noerr) &
             error stop 'Diagnostic write time: cannot append time coordinate'
       endif
       timestamp = timestamp + period_step
    enddo
end subroutine diagnostic_write_time
