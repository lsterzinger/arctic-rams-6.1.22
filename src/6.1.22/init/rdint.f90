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

Subroutine initlz ()

use mem_leaf
use mem_sib
use mem_grid
use mem_tracer
use mem_scratch
use mem_basic
use mem_micro
use var_tables
use mem_varinit
use mem_cuparm
use mem_oda
use io_params
use micphys

implicit none

!---------------------------------------------------------------------
!     *** This routine is the driver of all initialization packages
!---------------------------------------------------------------------

integer :: ifm,icm,ngr,nv,ierr,nsc

! Initialize aerosol density and vanthoff factors if they are used
if(iaerosol>0 .or. idust>0 .or. isalt>0) &
 CALL aerosol_init ()

if (trim(runtype) == 'INITIAL' .or. &
    trim(runtype) == 'ERROR'   .or. &
    trim(runtype) == 'HISTORY') then

!----------------------------------------------------------------------
!                 Initial startup
!----------------------------------------------------------------------
   
   print*,'  Initial start:', runtype, initial


   ! Read surface, topo, sst, and ndvi files for all grids. All the files
   !   were checked earlier, so they must be correct.,

   do ifm = 1,ngrids
      CALL top_read (ifm)
   enddo

   do ifm = 1,ngrids
      CALL sfc_read (ifm)
   enddo

   !     Define grid topography, transform, latitude-longitude,
   !        and map factor arrays.

   CALL grid_setup (2)
   
   ! read SST files

   CALL sst_read (1,ifm,ierr)
   if (ierr /= 0) then
      print*,'rdint: Error in  sst surface files'
      stop 'rdint: sst surface file error'
   endif

   ! read NDVI files

   CALL ndvi_read (1,ifm,ierr)
   if (ierr /= 0) then
      print*,'rdint: Error in  ndvi surface files'
      stop 'rdint: ndvi surface file error'
   endif
   
   ! The following things will be done for INITIAL = 1 or 3...
   
   if(initial==1 .or. initial==3) then

      ! If horizontally homogeneous initialization, 
      !    routine INITHH loops through all grids and initializes 
      !    those for which nxtnest = 0.

      if(initial == 1) then
         initorig=1
         print*,'----------------------------------------------------'
         print*,'Horizontally-homogeneous-INITIAL start of grid- 1' 
         print*,'----------------------------------------------------'
         CALL inithh ()
      endif
   
      ! If "history" initialization or restart, call INITHIS.  
      !      This will define initial fields and reference state on grid 1 from
      !      history file. Other grids will be interpolated as in a INITIAL=1 start.

      if(initial==3) then
         if(print_msg) then
          print*,'----------------------------------------------------'
         endif
         if(hrestart == 1) then
          CALL history_start ()
          print*,'History-RESTART'
         endif
         if(hrestart == 2) then
          CALL inithis (1)
          print*,'History-INITIALIZATION'
         endif
         if(print_msg) then
          print*,'----------------------------------------------------'
          print*, ''
         endif
      endif

      !  On all fine grids, initialize the surface layer characteristics,
      !  the 1-D reference state arrays, the 3-D reference state arrays,
      !  and the prognostic atmospheric fields by interpolation.

      CALL fmrefs1d (2,ngrids)

      do ifm = 2,ngrids
         icm = nxtnest(ifm)
         if (icm  >=  1) then
            CALL fmrefs3d (ifm)

            !Do not interpolate from grid-1 if doing history initialization
            if(initial==1) CALL prgintrp (nnzp(icm),nnxp(icm),nnyp(icm),ifm,1)

            CALL fmdn0 (ifm)
            print*,'----------------------------------------------------'            
            print*,'Model Start Initial interpolation of grid-',ifm
            print*,'----------------------------------------------------'
         endif
      enddo
 

   elseif(initial == 2) then
   
      ! If "variable initialization", do it all here
      initorig = 2
      CALL varf_read (0)

   endif

!     Initialize past time level velocity and perturbation Exner func
!     on all grids.

   do ifm=1,ngrids
      CALL newgrid (ifm)
      if(initial==1 .or. initial==2) then
        CALL fldinit (1)
        CALL negadj1 (nzp,nxp,nyp)
        CALL thermo (nzp,nxp,nyp,1,nxp,1,nyp)
      endif
      !Initialize aerosols
      if(initial==1 .or. initial==2 .or. (initial==3.and.iaerohist==1)) then
       if(iaerosol > 0) CALL init_ccn (nzp,nxp,nyp    &
          ,micro_g(ifm)%cccnp (1,1,1)  &
          ,micro_g(ifm)%cccmp (1,1,1)  &
          ,basic_g(ifm)%dn0   (1,1,1),ifm)
       if(iaerosol > 0) CALL init_gccn (nzp,nxp,nyp   &
          ,micro_g(ifm)%gccnp (1,1,1)  &
          ,micro_g(ifm)%gccmp (1,1,1)  &
          ,basic_g(ifm)%dn0   (1,1,1),ifm)
       if(idust  > 0)  CALL init_dust (nzp,nxp,nyp   &
          ,micro_g(ifm)%md1np (1,1,1)  &
          ,micro_g(ifm)%md2np (1,1,1)  &
          ,micro_g(ifm)%md1mp (1,1,1)  &
          ,micro_g(ifm)%md2mp (1,1,1),ifm)
       if(isalt  > 0)  CALL init_salt (nzp,nxp,nyp   &
          ,micro_g(ifm)%salt_film_np (1,1,1)  &
          ,micro_g(ifm)%salt_jet_np  (1,1,1)  &
          ,micro_g(ifm)%salt_spum_np (1,1,1)  &
          ,micro_g(ifm)%salt_film_mp (1,1,1)  &
          ,micro_g(ifm)%salt_jet_mp  (1,1,1)  &
          ,micro_g(ifm)%salt_spum_mp (1,1,1),ifm)
       if(level == 3) then
         if(ipris >= 5 .and. (iifn==1.or.iifn==2)) then
           CALL init_ifn (nzp,nxp,nyp    &
            ,micro_g(ifm)%cifnp (1,1,1)  &
            ,basic_g(ifm)%dn0   (1,1,1),ifm)
         endif
       endif
      endif

      !Initialize tracers
      if(initial==1 .or. initial==2 .or. (initial==3.and.itrachist==1)) then
       if(itracer > 0) then
        do nsc=1,itracer
           CALL init_tracer (nzp,nxp,nyp,tracer_g(nsc,ifm)%tracerp(1,1,1) &
             ,basic_g(ifm)%dn0(1,1,1),ifm,nsc)
        enddo
       endif
      endif

      !Initialize micro hydrometeor internal energy
      if (level == 3 .and. (initial==1.or.initial==2)) then
         CALL initqin (nzp,nxp,nyp        &
            ,micro_g(ifm)%q2    (1,1,1)  &
            ,micro_g(ifm)%q6    (1,1,1)  &
            ,micro_g(ifm)%q7    (1,1,1)  &
            ,basic_g(ifm)%pi0   (1,1,1)  &
            ,basic_g(ifm)%pp    (1,1,1)  &
            ,basic_g(ifm)%theta (1,1,1))
      endif

      !Initialize some SiB fields (CO2, etc)
      if (isfcl == 2 .and. (initial==1.or.initial==2)) then
         CALL init_sib1 (nzp,nxp,nyp        &
            ,sib_g(ifm)%rco2p   (1,1,1),ifm)
         CALL init_sib2 (nzp,nxp,nyp,npatch &
            ,sib_g(ifm)%rco2p   (1,1,1)     &
            ,sib_g(ifm)%pco2ap  (1,1,1)     &
            ,basic_g(ifm)%pi0   (1,1,1)     &
            ,basic_g(ifm)%pp    (1,1,1))
      endif

      !Initialize bubble perturbation via RAMSIN flags for HH run only
      if(initial==1) then
       if((ibubble==1.or.ibubble==2.or.ibubble==3) .and. ifm==ibubgrd) then
         CALL bubble (nzp,nxp,nyp &
           ,basic_g(ifm)%thp(1,1,1),basic_g(ifm)%rtp(1,1,1))
       endif
      endif

   enddo !end do for all grid initialization

! If initializing LEAF surface fields from history file.

   if((initial==1.or.initial==2).and.ipast_sfc == 1) then
     CALL inithis (0)
   endif

! Fill land surface data for all grids that have no standard input files.

   CALL sfcdata ()

! Initialize various surface variables.

   if((initial==1.or.initial==2).and.ipast_sfc == 0) then
     CALL geonest_nofile (1,ngrids)
   endif

! Reinitialize certain surface variables to prevent inconsistencies
! that arise from horizontal interpolation related to surface patches.
! This is not invoked if current and history grids match.

   if(initial==3 .or. ((initial==1.or.initial==2).and.ipast_sfc==1)) then
     if(hrestart == 2) CALL sfcinit_hstart ()
   endif

else

   stop 'WRONG RUNTYPE IN INITLZ'

endif

CALL micro_master ()

!       Fill latitude-longitude, map factor, and Coriolis arrays.

do ifm = 1,ngrids
   CALL newgrid (ifm)
   CALL fcorio (nxp,nyp           &
      ,basic_g(ifm)%fcoru (1,1)  &
      ,basic_g(ifm)%fcorv (1,1)  &
      ,grid_g(ifm)%glat   (1,1)  )
enddo


!  If we are doing one-way nesting or varfile nudging, inventory, 
!     prepare history/varfile files
!     and fill past/future nudging arrays for start of simulation

if(nud_type == 1) then
   CALL varf_read (1)
endif

! Process and read observations for ODA - observational data assimilation

if (if_oda == 1) CALL oda_read (nnzp(1),nnxp(1),nnyp(1) &
                              ,basic_g(1)%pi0(1,1,1)  &
                              ,scratch%scr1(1))

! Print locations of all grids

CALL gridloc_prt ()

!                  Save initial fields on history and analysis files
!                  -------------------------------------------------
CALL anal_write ('INST')
if(frqlite > 0.) CALL anal_write ('LITE')

!                  Save initial fields into the averaged arrays
!                  -------------------------------------------------
if(avgtim /= 0.)then

   do ngr=1,ngrids

      do nv=1,num_var(ngr)
         if(vtab_r(nv,ngr)%imean == 1) then
            CALL atob (vtab_r(nv,ngr)%npts,vtab_r(nv,ngr)%var_p  &
                     ,vtab_r(nv,ngr)%var_m)
         endif 
      enddo
   enddo
endif

!                  Print the header information and initial fields
!                  -----------------------------------------------

CALL prtopt ()

CALL nameout ()

CALL opspec3 ()

return
END SUBROUTINE initlz

!##############################################################################
Subroutine read_nl (file)

use mem_grid
use io_params

implicit none

character(len=*) :: file
real :: tfact
integer :: nf,n

! Initialize some variables to known values so we can reset them
!     after the read if they are not set in the namelist

frqstate(1:maxgrds)=-1.

! Read grid point and options information

open(1,status='OLD',file=file)

CALL namein (1,'$MODEL_GRIDS')
CALL namein (1,'$MODEL_FILE_INFO')

CALL namein (1,'$MODEL_OPTIONS')
CALL namein (1,'$MODEL_SOUND')
close(1)

!         Change some input time specifications into seconds

if(timeunit == 'd'.or.timeunit == 'D') tfact=86400.
if(timeunit == 'h'.or.timeunit == 'H') tfact=3600.
if(timeunit == 'm'.or.timeunit == 'M') tfact=60.
if(timeunit == 's'.or.timeunit == 'S') tfact=1.

timmax=timmax*tfact

! Set some variables if they weren't specified in namelist

! frqstate:
!----------------------------------
nf=0
do n = 1,ngrids
   if(iprntstmt>=1)print*,'Analysis file freqency on grid(',n,')',frqstate(n)
   if (frqstate(n) > 0.) nf = n
enddo

if(ioutput > 0) then
   if(nf == 0) then
      print*,'No good values for frqstate'
      stop 'read_nl: frqstate bad'
   elseif (nf /= ngrids) then
      print*,'Defaulting frqstate for grids ',nf+1,' to ',ngrids  &
            ,' to coarser grid value:',frqstate(nf)
      do n = nf+1 , ngrids
         frqstate(n) = frqstate(nf)
      enddo
   endif   
endif

return
END SUBROUTINE read_nl
