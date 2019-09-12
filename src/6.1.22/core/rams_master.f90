!
! Copyright (C) 1991-2004  ; All Rights Reserved ; Colorado State University
! Colorado State University Research Foundation ; ATMET, LLC
! 
! This file is free software; you can redistribute it and/or modify it under the
! terms of the GNU General Public License as published by the Free Software 
! Foundation; either version 2 of the License, or (at your option) any later version.
! 
! This software is distributed in the hope that it will be useful, but WITHOUT ANY 
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
! PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with this 
! code; if not, write to the Free Software Foundation, Inc., 
! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
!======================================================================================

Subroutine rams_master (iparall,nproc,taskids,master_num,name_name)

use rpara
use mem_grid
use mem_oda, only:if_oda
use io_params, only:hfilin

implicit none

integer taskids(*)
character(len=*) :: name_name
integer :: iparall,nproc,master_num
integer :: i,ifm,ns,nndtflg,namelength
real :: w1,w2,w3,t1,t2,wtime_start
real, external :: walltime

wtime_start = walltime(0.)
w1 = walltime(wtime_start)

iparallel=iparall

!          Read grid point and options information
!          ------------------------------------------------

!Read the RAMSIN namelist file
print*,'READ RAMSIN'
CALL read_nl (name_name)

!Set flags regarding history restart or initialization and set some
!initial condition flags with not a history restart.
hrestart=0
if(trim(runtype)=='HISTORY') then
 CALL read_hheader (hfilin)
 initial=3
 if(ngrids <= ngridsh) hrestart=1
 if(ngrids  > ngridsh) hrestart=2
else
 hrestart=0
 if(initial==3) hrestart=2
 initorig=0
 time=0.
 ngbegun(1:ngrids) = 0
endif

do ns=1,nproc
   if(iprntstmt>=1 .and. print_msg)print*,'CPU:',taskids(ns),master_num,iparall
enddo
   if(iprntstmt>=1 .and. print_msg)print*,'FILES:',name_name(1:len_trim(name_name))

!          Reset timmax for an ERROR run
if (trim(runtype) == 'ERROR') timmax = 0.

CALL eng_params ()  ! Various option settings that should normally 
                 !    not be changed

!          Reset parallel flag if necessary
!          --------------------------------------
if (trim(runtype)=='MAKESFC'   .or. trim(runtype)=='MAKEVFILE' .or. &
    trim(runtype)=='MAKEHFILE' .or. trim(runtype)=='ERROR') then
 if(iparallel /= 0) then
   if(iprntstmt>=1 .and. print_msg) &
      print*,'Resetting IPARALLEL=0 for a ',trim(runtype),' run'
   iparallel=0
 endif
endif

!          Print initial banner
!          ------------------------------------------------
!
namelength=len_trim(expnme)
write(6,'(a1,78a1)') ' ',('*',i=1,78)
write(6,'(2a1,a42)') ' ','*','    RAMS - Version 6.0'
write(6,'(2a1)') ' ','*'
write(6,'(3a1,A)') ' ','*',' ',trim(expnme)
write(6,'(2a1)') ' ','*'
write(6,'(2a1,A,A)') ' ','*',' RUNTYPE = ',trim(runtype)
write(6,'(a1,78a1)') ' ',('*',i=1,78)

! First check of options, mainly for numbers of grid points

CALL opspec1 ()

! Basic grid coordinate setup

CALL grid_setup (1)

! Additional checks, mainly for nesting 

CALL opspec2 ()      

! Check sfc,sst,ndvi files; remake if needed, even if RUNTYPE=ERROR

CALL make_sfcfiles ()

! Exit if a MAKESFC run

if (trim(runtype) == 'MAKESFC') then
   print*, 'MAKESFC run complete'
   go to 1000
endif

! If we are doing a "MAKEVFILE" run, call ISAN;
! If we are doing a "MAKEHFILE" run, call HISTORY NUDGING;
! then exit.
if (trim(runtype) == 'MAKEVFILE') then
   CALL isan_driver (name_name)
   print*,' ISAN complete '
   go to 1000
endif   
if (trim(runtype) == 'MAKEHFILE') then
   CALL nudh_driver (name_name)
   print*,' History Nudging Varfiles Complete '
   go to 1000
endif


!-----------------------------------------------------------
! If we got here, we are doing an actual simulation
!-----------------------------------------------------------

! Initialize micro arrays. May need to change some settings which affect memory.
CALL jnmbinit ()

! Allocate main memory

print*, '---------------------------------------------------'
if (iparallel == 0) then
   CALL rams_mem_alloc (0) !     Allocate new data types
else
   CALL rams_mem_alloc (1) !     Allocate new data types
endif
print*, '---------------------------------------------------'

!          call the main initialization driver
!          -----------------------------------

CALL initlz (name_name)

! Compute Courant numbers cflxy and cflz.

do ifm = 1,ngrids
   CALL newgrid (ifm)
   CALL cfl (nzp,nxp,nyp,0,0,0)
enddo

! Initialize dtlongn, nndtrat, and nnacoust, and compute the timestep
! schedule for all grid operations.

CALL dtset (nndtflg)
CALL modsched (isched,maxsched,ngrids,nxtnest,nndtrat,nsubs)

CALL timing (1,t1)
w2 = walltime(wtime_start)
if(iprntstmt>=1  .and. print_msg) &
   print'(a,2f12.3)','CPU - wall time: master init: ',t1,w2-w1

if(iparallel == 1) then

!   ---------------------------------------------------------------
!     Initialize parallel processes with all relevant information
!   --------------------------------------------------------------
      if(iprntstmt>=1 .and. print_msg) &
         print*, '++++rams_master start init',nzg,nzs,npatch,nnzp(1)

   CALL par_ready (nmachs,machnum,777)
      if(iprntstmt>=1 .and. print_msg) &
         print*, '++++rams_master-ready1',nzg,nzs,npatch,nnzp(1)

   CALL masterput_processid (nproc,taskids,master_num)
      if(iprntstmt>=1 .and. print_msg) &
         print*, '++++rams_master1-processid',nzg,nzs,npatch,nnzp(1)

   CALL masterput_nl ()
      if(iprntstmt>=1 .and. print_msg) &
         print*, '++++rams_master-nl',nzg,nzs,npatch,nnzp(1)

   CALL masterput_gridinit ()
      if(iprntstmt>=1 .and. print_msg) &
         print*, '++++rams_master-gridinit',nzg,nzs,npatch,nnzp(1)

   CALL node_decomp ()
      if(iprntstmt>=1 .and. print_msg) &
         print*, '++++rams_master-decomp',nzg,nzs,npatch,nnzp(1)

   CALL masterput_grid_dimens ()
      if(iprntstmt>=1 .and. print_msg) &
         print*, '++++rams_master-grid_dimens'

   CALL masterput_gridset ()
      if(iprntstmt>=1 .and. print_msg) &
         print*, '++++rams_master-gridset',nzg,nzs,npatch,nnzp(1)

   CALL masterput_cofnest ()
      if(iprntstmt>=1 .and. print_msg) &
         print*, '++++rams_master-cofnest',nzg,nzs,npatch,nnzp(1)

   CALL masterput_micphys ()
      if(iprntstmt>=1 .and. print_msg) &
         print*, '++++rams_master-micphys',nzg,nzs,npatch,nnzp(1)

   if (if_oda == 1) CALL masterput_oda ()
      if(iprntstmt>=1 .and. print_msg) &
         print*, '++++rams_master-oda',nzg,nzs,npatch,nnzp(1)

   CALL masterput_misc ()
      if(iprntstmt>=1 .and. print_msg) &
         print*, '++++rams_master-misc',nzg,nzs,npatch,nnzp(1)

   CALL master_sendinit ()
      if(iprntstmt>=1 .and. print_msg) &
         print*, '++++rams_master-sendinit',nzg,nzs,npatch,nnzp(1)

   CALL par_ready (nmachs,machnum,777)
      if(iprntstmt>=1 .and. print_msg) &
         print*, '++++rams_master-ready2',nzg,nzs,npatch,nnzp(1)

endif

CALL timing (1,t2)
w3 = walltime(wtime_start)
if(iprntstmt>=1 .and. print_msg) &
   print '(a,2f12.3)','CPU - wall time: node init: ',t2-t1,w3-w2
if(iprntstmt>=1 .and. print_msg) print*, 'rams_master-done init'

! Exit if doing a zero time run
if (time >= timmax .or. trim(runtype) == 'ERROR') go to 1000

!  call the model time integration driver
!  --------------------------------------

if(iparallel == 1) then
   CALL par_model ()
else
   CALL model ()
endif

!  RAMS finished, clean up some last things...
!  -----------------------------------------------------------

1000 continue

if (trim(runtype) == 'ERROR') then
   print*
   print*,'---------------------------------------------------------'
   print*,'|  ERROR run completed successfully. No fatal errors.'
   print*,'---------------------------------------------------------'
   print*
endif


if (iparallel == 1) CALL par_exit ()

return
END SUBROUTINE rams_master

!##############################################################################
Subroutine comm_time (isendflg,isendlite,isendmean,isendboth)

use mem_varinit
use mem_cuparm
use io_params
use mem_grid

implicit none

integer :: isendflg,isendlite,isendmean,isendboth

real :: timemf
integer :: ifm

!         ISENDFLG designates whether nodes should send back
!            stuff things it normally doesn't have to
!            at the end of timestep for history/analysis write.

!         isendflg  = the usual RAMS stuff
!         isendlite = the "lite" variables
!         isendmean = the "mean" variasbles
!         isendboth = Both the "mean" and "lite" variables

!            Determines whether nodes send stuff back at the END of the
!            timestep!!!!!

timemf = time + dtlongn(1)
isendflg = 0
isendlite = 0
isendmean = 0
isendboth = 0

if(frqlite > 0.) then
   if (mod(timemf,frqlite) < dtlongn(1)) isendlite=1
endif

if (frqmean > 0.) then
   if(avgtim > 0.)then
      if(mod(timemf-avgtim/2.,frqmean) < dtlongn(1) .and.  &
         timemf >= avgtim) isendmean=1
   elseif(avgtim < 0.)then
      if(mod(timemf,frqmean) < dtlongn(1)) isendmean=1
   endif
endif

if (frqboth > 0.) then
   if(avgtim > 0.)then
      if(mod(timemf-avgtim/2.,frqboth) < dtlongn(1) .and.  &
         timemf >= avgtim) isendboth=1
   elseif(avgtim < 0.)then
      if(mod(timemf,frqboth) < dtlongn(1)) isendboth=1
   endif
endif

if (ioutput  /=  0) then
   do ifm = 1,ngrids
      if ( mod(timemf,frqstate(ifm))  <  dtlongn(1) ) then
            isendflg = 1
            return
      endif
   enddo
endif

if( timemf  >=  timmax - .01*dtlongn(1) ) then
   isendflg = 1
   return
endif

if  ( nud_type == 1 .and. timemf  >=  vtime2  &
    .and. timemf  <  timmax) then
   isendflg = 1
   return
endif

if  ( nud_cond == 1 .and. timemf  >=  vtime2  &
    .and. timemf  <  timmax) then
   isendflg = 1
   return
endif

if (iupdsst  ==  1 ) then
   do ifm = 1,ngrids
      if (isstflg(ifm)  ==  1) then
         if (timemf  >=  ssttime2(ifm) .and.  &
            timemf  <  timmax) then
            isendflg = 1
            return
         endif
      endif
   enddo
endif

if (iupdndvi  ==  1 ) then
   do ifm = 1,ngrids
      if (ndviflg(ifm)  ==  1) then
         if (timemf  >=  ndvitime2(ifm) .and.  &
            timemf  <  timmax) then
            isendflg = 1
            return
         endif
      endif
   enddo
endif

return
END SUBROUTINE comm_time

!##############################################################################
Subroutine rams_output ()

use mem_leaf
use mem_varinit
use mem_cuparm
use io_params
use mem_grid

implicit none

integer :: ierr,ifm
logical :: analwrite
real :: timmaxr


do ifm = 1,ngrids
   CALL newgrid (ifm)
   timmaxr = timmax - .01*dtlongn(ifm)

! Implement call to NON_SCALAR_BC if calling ANLWRT.
! This sets all BC's for non-scalar (non-advected) quantities such
! as THETA, RV, and variables in LEAF, RADIATION, TURBULENCE.

!        For testing the mean fields, skip this rediagnosis
!        GOTO 10

   if (ioutput  ==  0) go to 10

   analwrite=.false.

   !FOR STATE FILE OUTPUT
   if(mod(time,frqstate(ifm)) <  dtlongn(ifm)               .or.  &
                        time  >= timmaxr                    .or.  &
                        iflag == 1                                &
      ) analwrite=.true.
   !FOR LITE FILE OUTPUT
   if(frqlite > 0.) then
      if( mod(time,frqlite) < dtlongn(1).or.  &
                         time  >=  timmaxr ) then
        analwrite=.true.
      endif
   endif
   !FOR MEAN FILE OUTPUT
   if (frqmean > 0.) then
      if(avgtim > 0.0.and.mod(time-avgtim/2.,frqmean) < dtlongn(1)  &
         .and.time >= avgtim) then
        analwrite=.true.
      endif
      if(avgtim < 0.0.and.mod(time,frqmean) < dtlongn(1)) then
        analwrite=.true.
      endif
   endif
   !FOR BOTH LITE AND MEAN FILE OUTPUT
   if (frqboth>0.) then
      if(avgtim > 0.0.and.mod(time-avgtim/2.,frqboth) < dtlongn(1)  &
         .and.time >= avgtim) analwrite=.true.
      if(avgtim < 0.0.and.mod(time,frqboth) < dtlongn(1))  &
         analwrite=.true.
   endif

   !If we are writing any output files, then do the non-scalar BC's.
   !We have to do this call within the ifm grid loop after call to 
   !"newgrid" so we have the grid appropriate nxp,nyp,nzp.
   if (analwrite) then
     CALL non_scalar_bc (nzp,nxp,nyp,nzg,nzs)
   endif

10      continue

enddo

timmaxr = timmax - .01*dtlongn(1)

analwrite=.false.
do ifm = 1,ngrids
   if(mod(time,frqstate(ifm)) < dtlongn(1).or.  &
      time  >=  timmaxr .or.  &
      iflag == 1) analwrite=.true.
enddo

if (analwrite) CALL anal_write ('INST')

!     call the analysis writing routine again for the other var types
if(frqlite > 0.) then
   if( mod(time,frqlite) < dtlongn(1).or.  &
                      time  >=  timmaxr ) then
    CALL anal_write ('LITE')
   endif
endif

if (frqmean > 0.) then
   if(iprntstmt>=1 .and. print_msg) &
      print*, 'check avg',time,frqmean,avgtim,dtlongn(1)
   if(avgtim > 0.0.and.mod(time-avgtim/2.,frqmean) < dtlongn(1)  &
      .and.time >= avgtim) then
      CALL anal_write ('MEAN')
   endif
      
   if(avgtim < 0.0.and.mod(time,frqmean) < dtlongn(1)) then
      CALL anal_write ('MEAN')
   endif
endif

if (frqboth>0.) then
   if(avgtim > 0.0.and.mod(time-avgtim/2.,frqboth) < dtlongn(1)  &
      .and.time >= avgtim) CALL anal_write ('BOTH')
   if(avgtim < 0.0.and.mod(time,frqboth) < dtlongn(1))  &
      CALL anal_write ('BOTH')
endif

timmaxr = time+.00001

if (iupdsst  ==  1 .and. timmaxr < timmax) then
   do ifm = 1,ngrids
      CALL sst_read (3,ifm,ierr)
   enddo
endif

if (iupdndvi  ==  1 .and. timmaxr < timmax) then
   do ifm = 1,ngrids
      CALL ndvi_read (3,ifm,ierr)
   enddo   
endif

if(nud_type == 1 .and. time >= vtime2 .and. timmaxr < timmax) then
   CALL varf_read (2)
endif

if (iflag  ==  1) stop 'IFLAG'

return
END SUBROUTINE rams_output
