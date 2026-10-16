subroutine firegfed_read_single(gprogt)

!****--------------------------------------------------------------------
!    This subroutines reads the fire burned area file 
!    for the current time step for single site (point) simulations.
!    If required, it closes the current month's data file and opens the 
!    next month's data file.
!****--------------------------------------------------------------------

use kinds
use module_io
use module_sib, only: gprog_type
use module_sibconst, only: &
    print_fire, print_stop
use module_time, only: &
     year, month, day, hour
use module_pparams, only:  &
    secs_per_day

implicit none

!...input variables
type(gprog_type), intent(inout) :: gprogt

!...local variables
integer(i4) :: yr, doy
real(r4) :: hr
real(r8) :: ba_dbgin, ba_en2in, ba_en3in, ba_dnfin, ba_eb1in, ba_eb2in
real(r8) :: ba_db1in, ba_db2in, ba_db3in, ba_shbin, ba_shain, ba_c3ain
real(r8) :: ba_c3gin, ba_c4gin, ba_c3cin, ba_c4cin, ba_mzein, ba_soyin, ba_wwtin

!...misc variables
integer(i4) :: status
logical :: opened, exist
character(len=1024) :: record

!...storing previous data
gprogt%ba_dbg1 = gprogt%ba_dbg2
gprogt%ba_en21 = gprogt%ba_en22
gprogt%ba_en31 = gprogt%ba_en32
gprogt%ba_dnf1 = gprogt%ba_dnf2
gprogt%ba_eb11 = gprogt%ba_eb12
gprogt%ba_eb21 = gprogt%ba_eb22
gprogt%ba_db11 = gprogt%ba_db12
gprogt%ba_db21 = gprogt%ba_db22
gprogt%ba_db31 = gprogt%ba_db32
gprogt%ba_shb1 = gprogt%ba_shb2
gprogt%ba_sha1 = gprogt%ba_sha2
gprogt%ba_c3a1 = gprogt%ba_c3a2
gprogt%ba_c3g1 = gprogt%ba_c3g2
gprogt%ba_c4g1 = gprogt%ba_c4g2
gprogt%ba_c3c1 = gprogt%ba_c3c2
gprogt%ba_c4c1 = gprogt%ba_c4c2
gprogt%ba_mze1 = gprogt%ba_mze2
gprogt%ba_soy1 = gprogt%ba_soy2
gprogt%ba_wwt1 = gprogt%ba_wwt2

!...switch files if needed
if (firegfed_switchf) then !switch and read
    inquire(unit=fireid, exist=exist,opened=opened)
    if (opened) close(fireid, iostat=status)

    write(unit=firegfed_filename,fmt='(a,i4.4,i2.2,a4)') &
        trim(firegfed_path), fire_year, fire_month, '.dat'
    inquire(file=trim(firegfed_filename), exist=exist)

    if (exist) then
       open(unit=fireid,file=trim(firegfed_filename), status='old', &
            form='formatted', iostat=status)
       if (status > 0) then
          print*,'!!!Error opening fire file!!!'
          stop
       endif
     else !fire file does not exist
         if (firefile_stop) then
             print*,''
             print*,'Stopping due to non-existent fire file: '
             print*,' ',trim(firegfed_filename)
             stop
         else
             fire_step = 0
             fireid = 0
             gprogt%ba_dbg2 = dzero
             gprogt%ba_en22 = dzero
             gprogt%ba_en32 = dzero
             gprogt%ba_dnf2 = dzero
             gprogt%ba_eb12 = dzero
             gprogt%ba_eb22 = dzero
             gprogt%ba_db12 = dzero
             gprogt%ba_db22 = dzero
             gprogt%ba_db32 = dzero
             gprogt%ba_shb2 = dzero
             gprogt%ba_sha2 = dzero
             gprogt%ba_c3a2 = dzero
             gprogt%ba_c3g2 = dzero
             gprogt%ba_c4g2 = dzero
             gprogt%ba_c3c2 = dzero
             gprogt%ba_c4c2 = dzero
             gprogt%ba_mze2 = dzero
             gprogt%ba_soy2 = dzero
             gprogt%ba_wwt2 = dzero

             RETURN
         endif
           
     endif 
endif !switch files

if (fire_recnum .gt. 0) then
   DO  !Read until not a comment
       read(fireid,'(a)',iostat=status) record
       IF (status /= 0) THEN
           print*,'Stopping due to error in single fire data.'
           stop
       ENDIF
       IF (record(1:1) .ne. '#') THEN
          exit
       ENDIF
   ENDDO

   read(unit=record,fmt=*) yr, doy, hr, &
        ba_dbgin, ba_en2in, ba_en3in, ba_dnfin, ba_eb1in, ba_eb2in, ba_db1in, &
        ba_db2in, ba_db3in, ba_shbin, ba_shain, ba_c3ain, ba_c3gin, ba_c4gin, &
        ba_c3cin, ba_c4cin, ba_mzein, ba_soyin, ba_wwtin 
else
   ba_dbgin=dzero
   ba_en2in=dzero
   ba_en3in=dzero
   ba_dnfin=dzero
   ba_eb1in=dzero
   ba_eb2in=dzero
   ba_db1in=dzero
   ba_db2in=dzero
   ba_db3in=dzero
   ba_shbin=dzero
   ba_shain=dzero
   ba_c3ain=dzero
   ba_c3gin=dzero
   ba_c4gin=dzero
   ba_c3cin=dzero
   ba_c4cin=dzero
   ba_mzein=dzero
   ba_soyin=dzero
   ba_wwtin=dzero
endif

!...set the new variables
!...note: for BA, driver file units are BA % PFT⁻¹ 3h⁻¹, convert 3h⁻¹ to sec⁻¹ 
!... test removing this in case it's casuing a problem
gprogt%ba_dbg2 = ba_dbgin*(1./(3.*60.*60.))
gprogt%ba_en22 = ba_en2in*(1./(3.*60.*60.))
gprogt%ba_en32 = ba_en3in*(1./(3.*60.*60.))
gprogt%ba_dnf2 = ba_dnfin*(1./(3.*60.*60.))
gprogt%ba_eb12 = ba_eb1in*(1./(3.*60.*60.))
gprogt%ba_eb22 = ba_eb2in*(1./(3.*60.*60.))
gprogt%ba_db12 = ba_db1in*(1./(3.*60.*60.))
gprogt%ba_db22 = ba_db2in*(1./(3.*60.*60.))
gprogt%ba_db32 = ba_db3in*(1./(3.*60.*60.))
gprogt%ba_shb2 = ba_shbin*(1./(3.*60.*60.))
gprogt%ba_sha2 = ba_shain*(1./(3.*60.*60.))
gprogt%ba_c3a2 = ba_c3ain*(1./(3.*60.*60.))
gprogt%ba_c3g2 = ba_c3gin*(1./(3.*60.*60.))
gprogt%ba_c4g2 = ba_c4gin*(1./(3.*60.*60.))
gprogt%ba_c3c2 = ba_c3cin*(1./(3.*60.*60.))
gprogt%ba_c4c2 = ba_c4cin*(1./(3.*60.*60.))
gprogt%ba_mze2 = ba_mzein*(1./(3.*60.*60.))
gprogt%ba_soy2 = ba_soyin*(1./(3.*60.*60.))
gprogt%ba_wwt2 = ba_wwtin*(1./(3.*60.*60.))

!gprogt%ba_dbg2 = ba_dbgin
!gprogt%ba_en22 = ba_en2in
!gprogt%ba_en32 = ba_en3in
!gprogt%ba_dnf2 = ba_dnfin
!gprogt%ba_eb12 = ba_eb1in
!gprogt%ba_eb22 = ba_eb2in
!gprogt%ba_db12 = ba_db1in
!gprogt%ba_db22 = ba_db2in
!gprogt%ba_db32 = ba_db3in
!gprogt%ba_shb2 = ba_shbin
!gprogt%ba_sha2 = ba_shain
!gprogt%ba_c3a2 = ba_c3ain
!gprogt%ba_c3g2 = ba_c3gin
!gprogt%ba_c4g2 = ba_c4gin
!gprogt%ba_c3c2 = ba_c3cin
!gprogt%ba_c4c2 = ba_c4cin
!gprogt%ba_mze2 = ba_mzein
!gprogt%ba_soy2 = ba_soyin
!gprogt%ba_wwt2 = ba_wwtin

!...print out the new data if requested
if (print_fire) then
    print*,'     ------------------------------------'
    print*,'     !!!READING BURNED AREA!!!'

    if (firegfed_switchf) then
        print*,'   Opening fire file: '
        print*,'       ',trim(firegfed_filename)
    endif

    print*,'   REAL TIMES:'
    print*,'   Year    Mon   Day   Hour'
    print('(a,i5,a,i4,a,i4,a,i4)'), '   ', &
         year,'  ',month,'  ',day,'  ',hour
    print*,''
    print*,'   FIRE TIMES:'
    print*,'   Year    Mon   Day   Hour'
    print('(a,i5,a,i4,a,i4,a,f6.2)'), '   ', &
          fire_year,'  ',fire_month,'  ',fire_day, &
          '  ',fire_hour
    print*,'   FIRE RECNUM/TOTNUM: ', fire_recnum, fire_permon
    !print('(a,e16.6)'),'     Fire C Loss (umol C/m2/s)    : ', &
    !          firecin
    !print('(a,e16.6)'),'     Fire CO2 Resp (umol CO2/m2/s): ', &
    !          fireco2in
    print('(a,e16.6)'),'     BA DBG BA (m2 PFT⁻¹ day⁻¹)    : ', &
              ba_dbgin
    print('(a,e16.6)'),'     BA EN2 BA (m2 PFT⁻¹ day⁻¹)    : ', &
              ba_en2in
    print('(a,e16.6)'),'     BA EN3 BA (m2 PFT⁻¹ day⁻¹)    : ', &
              ba_en3in
    print('(a,e16.6)'),'     BA DNF BA (m2 PFT⁻¹ day⁻¹)    : ', &
              ba_dnfin
    print('(a,e16.6)'),'     BA EB1 BA (m2 PFT⁻¹ day⁻¹)    : ', &
              ba_eb1in
    print('(a,e16.6)'),'     BA EB2 BA (m2 PFT⁻¹ day⁻¹)    : ', &
              ba_eb2in
    print('(a,e16.6)'),'     BA DB1 BA (m2 PFT⁻¹ day⁻¹)    : ', &
              ba_db1in
    print('(a,e16.6)'),'     BA DB2 BA (m2 PFT⁻¹ day⁻¹)    : ', &
              ba_db2in
    print('(a,e16.6)'),'     BA DB3 BA (m2 PFT⁻¹ day⁻¹)    : ', &
              ba_db3in
    print('(a,e16.6)'),'     BA SHB BA (m2 PFT⁻¹ day⁻¹)    : ', &
              ba_shbin
    print('(a,e16.6)'),'     BA SHA BA (m2 PFT⁻¹ day⁻¹)    : ', &
              ba_shain
    print('(a,e16.6)'),'     BA C3A BA (m2 PFT⁻¹ day⁻¹)    : ', &
              ba_c3ain
    print('(a,e16.6)'),'     BA C3G BA (m2 PFT⁻¹ day⁻¹)    : ', &
              ba_c3gin
    print('(a,e16.6)'),'     BA C4G BA (m2 PFT⁻¹ day⁻¹)    : ', &
              ba_c4gin
    print('(a,e16.6)'),'     BA C3C BA (m2 PFT⁻¹ day⁻¹)    : ', &
              ba_c3cin
    print('(a,e16.6)'),'     BA C4C BA (m2 PFT⁻¹ day⁻¹)    : ', &
              ba_c4cin
    print('(a,e16.6)'),'     BA MZE BA (m2 PFT⁻¹ day⁻¹)    : ', &
              ba_mzein
    print('(a,e16.6)'),'     BA SOY BA (m2 PFT⁻¹ day⁻¹)    : ', &
              ba_soyin
    print('(a,e16.6)'),'     BA WWT BA (m2 PFT⁻¹ day⁻¹)    : ', &
              ba_wwtin
    print*,'     ------------------------------------'
 
    if (print_stop) stop
endif


end subroutine firegfed_read_single
