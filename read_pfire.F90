
! Opens and reads in SiB pool parameters.
subroutine read_pfire()

use kinds
use module_sibconst, only: &
    cornsoy_switch, &
    npft, ngroup, &
    npoolpft, npoollu, ntpool
use module_io, only: &
    firep_file, firepid
use module_pparams, only: &
    secs_per_day, days_per_year, &
    drytoc, mwc
use module_pftinfo, only: &
    clen,npft_gdd, &
    pft_mze, pft_soy, &
    pft_group, group_grass, group_crop
use module_poolinfo, only: &
    pool_indx_leaf, pool_indx_leaf_c13
use module_param, only: &
    physcon, firecon
!use module_phosib, only: c4

implicit none

!...file variables
integer(i4) :: finpft, fingddpft
integer(i4) :: finpool, fingroup
integer(i4) :: num
character(len=clen), dimension(npft) :: pftname
character(len=10), dimension(ntpool) :: poolname
character(len=100) :: trash
logical :: iscomment1, iscomment2

real(r4), dimension(npoollu+2) :: graze_trans, harvest_trans

!...misc variables
integer(i4) :: i,j
integer(byte) :: groupref
real(r8) :: poolval

!...C, C13 and C14 are stored explicitly as three equal pool blocks.
if ((mod(npoolpft,3) /= 0) .or. (mod(ntpool,3) /= 0)) then
    print*,'Fire parameter reader requires C/C13/C14 pool blocks of equal size.'
    stop
endif

!-------------------
!...Initialize the pool variables
allocate(firecon(npft))
call init_firecon(npft, npoolpft, npoollu, firecon)

!...Open file
print*,'Reading Fire Parameter File: '
print*,'  ',trim(firep_file)
open(unit=firepid,file=trim(firep_file),form='formatted')

read(firepid,*) finpft
if (finpft /= npft) then
    print*,''
    print('(a)'),'!!!Fire param File Mistmatch!!!'
    print('(a,i4,a,i4)'),'  Firep file npft: ', &
          finpft,' Sim npft: ',npft
    print*,''
    stop
endif

read(firepid,*) fingddpft
if (fingddpft /= npft_gdd) then
    print*,''
    print('(a)'),'!!!Fire param File Mistmatch!!!'
    print('(a,i4,a,i4)'),'  Firep file npft_gdd: ', &
          fingddpft,' Sim npft_gdd: ',npft_gdd
    print*,''
    stop
endif
read(firepid,*) finpool
read(firepid,*) fingroup

if (finpool /= ntpool) then
    print*,''
    print('(a)'),'!!!Fire param file does not match C/C13/C14 simulation!!!'
    print('(a,i4,a,i4)'),'  Firep file npool: ',finpool, &
         ' Sim npool: ',ntpool
    print*,''
    stop
endif

if (fingroup /= ngroup) then
    print*,''
    print('(a)'),'!!!Fire param file does not match simulation!!!'
    print('(a,i4,a,i4)'),'  Firep file ngroup: ',fingroup,' Sim ngroup: ',ngroup
    print*,''
    stop
endif

iscomment1=.true.
iscomment2=.true.
do while ((iscomment1) .or. (iscomment2))
    read(firepid,*) trash
    if (index(trash,'****') .gt. 0) then
        if (iscomment1) then
           iscomment1=.false.
        else 
           iscomment2=.false.
        endif
    endif
enddo


read(firepid,*) trash
read(firepid,*) trash
read(firepid,*) trash
do i=1, npft
   read(firepid,*) num, pftname(i), firecon(i)%mfl(1:npoolpft)
!   print*,'pftname, MFL: ',pftname(i),firecon(i)%mfl(1),firecon(i)%mfl(2),firecon(i)%mfl(3),firecon(i)%mfl(4),firecon(i)%mfl(5)
enddo

read(firepid,*) trash
read(firepid,*) trash
read(firepid,*) trash
do i=1, npft
   read(firepid,*) num, pftname(i), firecon(i)%mfh(1:npoolpft)
!   print*,'pftname, MFH: ',pftname(i),firecon(i)%mfh(1),firecon(i)%mfh(2),firecon(i)%mfh(3),firecon(i)%mfh(4),firecon(i)%mfh(5)
enddo

read(firepid,*) trash
read(firepid,*) trash
read(firepid,*) trash
do i=1, npft
   read(firepid,*) num, pftname(i), firecon(i)%ccl(1:ntpool)
!   print*,'pftname, CCL:',pftname(i),firecon(i)%ccl(1),firecon(i)%ccl(2),firecon(i)%ccl(3),firecon(i)%ccl(4),firecon(i)%ccl(5),firecon(i)%ccl(6),firecon(i)%ccl(7),firecon(i)%ccl(8),firecon(i)%ccl(9),firecon(i)%ccl(10),firecon(i)%ccl(11)
enddo

read(firepid,*) trash
read(firepid,*) trash
read(firepid,*) trash
do i=1, npft
   read(firepid,*) num, pftname(i), firecon(i)%cch(1:ntpool)
!   print*,'pftname, CCH:',pftname(i),firecon(i)%cch(1),firecon(i)%cch(2),firecon(i)%cch(3),firecon(i)%cch(4),firecon(i)%cch(5),firecon(i)%cch(6),firecon(i)%cch(7),firecon(i)%cch(8),firecon(i)%cch(9),firecon(i)%cch(10),firecon(i)%cch(11)

enddo

end subroutine read_pfire

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Initializes the pool parameters.
subroutine init_firecon(npft, npoolpft, npoollu, &
     firecon)

  use kinds
  use module_param, only: fire_param

  !...input variables
  integer(i4), intent(in) :: npft, npoolpft, npoollu
  type(fire_param), dimension(npft), intent(inout) :: firecon

  !...local variables
  integer(i4) :: i, ntpool

  !...set local variables
  ntpool = npoolpft + npoollu

  !...initialize pool parameters
  do i=1, npft

     allocate(firecon(i)%mfl(npoolpft))
     firecon(i)%mfl(:) = dzero
     allocate(firecon(i)%mfh(npoolpft))
     firecon(i)%mfh(:) = dzero

     allocate(firecon(i)%ccl(ntpool))
     firecon(i)%ccl(:) = dzero
     allocate(firecon(i)%cch(ntpool))
     firecon(i)%cch(:) = dzero

  enddo !i=1,npft

end subroutine init_firecon
