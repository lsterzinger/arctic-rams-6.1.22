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

Subroutine model ()

use mem_grid

implicit none

!   +------------------------------------------------------------------
!   ! This routine drives the entire time integration process
!   !   for a non-parallel run.
!   +------------------------------------------------------------------

integer :: npass,nndtflg,icm,ifm,nfeed
real :: wtime_start,begtime,t1,wtime1,wtime2,t2,wtime_tot,timeh
real, external :: walltime

print*, 'starting routine MODEL'

wtime_start = walltime(0.)
istp = 0

!         Start the timesteps

do while (time < timmax)

   istp = istp + 1
   begtime=time

!            CPU timing information
   CALL timing (1,t1)
   wtime1 = walltime(wtime_start)

! Examine Courant numbers in case model needs to be stopped
! or (if ideltat < 0), to update dtlongn, nndtrat,
! nnacoust, sspct and isched.

   CALL dtset (nndtflg)
   if (nndtflg .gt. 0) then
      CALL modsched (isched,maxsched,ngrids,nxtnest,nndtrat,nsubs)
   endif

!----------------------------------------------------------------------------
!                  Loop through all grids and advance a 'DTLONG' timestep.
!----------------------------------------------------------------------------
!                  Start the timestep schedule

   do npass=1,nsubs

      isstp=isched(npass,3)
      ngrid=isched(npass,1)
      CALL newgrid (ngrid)

!---------------------------------------------------------------------
!                  Advance this grid forward by the appropriate timestep.

      time=begtime + (isched(npass,5)-1) * dtlt

      !Call main timestep driver
      CALL onenode ()
      CALL timestep ()

      ngbegun(ngrid)=1

!----------------------------------------------------------------------------
!                  Is it time to send the coarse domain points to the nested
!                  nodes to interpolate a nested grid's boundaries?

      if(isched(npass,2).ne.0) then
         CALL newgrid (isched(npass,2))
         isstp=isched(npass,3)
         icm=nxtnest(ngrid)
         CALL prgintrp (nnzp(icm),nnxp(icm),nnyp(icm),ngrid,0)
      endif

!----------------------------------------------------------------------------
!                  Is it time to feedback fields from fine grids to
!                     coarser grids?

     !Saleeby(2012): added flag for choosing 1 or 2 way nesting
     !Nesting feedback call
     if(inesting.eq.1) then 
      if(isched(npass,4).ne.0) then
         ifm = isched(npass,1)
         icm = nxtnest(ifm)
         do nfeed = 1,isched(npass,4)
            CALL newgrid (ifm)
            CALL nstfeed (ifm,icm)
            ifm = icm
            icm = nxtnest(ifm)
         enddo
      endif
     endif

   enddo
!                  End the timestep schedule

! All nesting operations (interpolation and feedback) have been carried out.

!----------------------------------------------------------------------------
! Compute Courant numbers cflxy and cflz and do averaging.

   do ifm = 1,ngrids
      CALL newgrid (ifm)
      CALL cfl (nzp,nxp,nyp,0,0,1)
      !THETAM and RVM have not been updated after nesting feedback
      !This means that these variables are really a timestep
      !behind the instantaneous variables.
      !Calculate the means
      CALL anlavg (nzp,nxp,nyp)
   enddo

!----------------------------------------------------------------------------
! Send timing info.

   wtime2 = walltime(wtime_start)
   CALL timing (2,T2)

!---------------------------------------------------------------------------
! Update main time variable by a long timestep.

   time=begtime+dtlongn(1)

!---------------------------------------------------------------------------
! Print timestep info and see if we write output now.

   timeh = TIME / 3600.
   PRINT 201, ISTP,TIME,timeh,T2-T1,wtime2-wtime1
201     FORMAT(' Timestep-',I5,'   Sim time(sec)=',F9.1  &
        ,'  (hr)=',F6.2,'  CPU(sec)=',F7.2,'  Wall(sec)=',F7.2)

   CALL rams_output ()

enddo

wtime_tot = walltime(wtime_start)
print '(//,a,f10.0)'  &
     ,' -----Total elapsed time: ',wtime_tot

return
END SUBROUTINE model

!##############################################################################
Subroutine par_model ()

use mem_grid
use rpara

implicit none
!   +------------------------------------------------------------------
!   ! This routine drives the entire time integration process
!   !   for a parallel run.
!   +------------------------------------------------------------------

integer :: isendflg,isendlite,isendmean,isendboth,nndtflg,ntsend,n
real :: wtime_start,begtime,t1,wtime1,wtime2,t2,wtime_tot,pcpu,pwall,timeh
real, external :: walltime

print*, 'starting routine par_MODEL', time,timmax

wtime_start = walltime(0.)
istp = 0

isendflg=0
isendlite = 0
isendmean = 0
isendboth = 0

!         Start the main timestep loop

do while (time < timmax)

   istp = istp + 1
   begtime=time

   if(isendflg.eq.1) then
      CALL master_sendinit ()
   endif

! CPU timing information
   CALL timing (1,t1)
   wtime1 = walltime(wtime_start)

! ISENDFLG designates whether nodes should send back
! stuff things it normally does not have to
! at the end of timestep for history/analysis write.

! Determines whether nodes send stuff back at the END of the timestep!!

   CALL comm_time (isendflg,isendlite,isendmean,isendboth)

! Examine Courant numbers in case model needs to be stopped or
! (if ideltat < 0), to update dtlongn, nndtrat,  nnacoust, sspct and isched.

   CALL dtset (nndtflg)
   if (iflag > 0) then
      isendflg = 1
      isendlite = 1
      isendmean = 1
      isendboth = 1
   endif
   if (nndtflg > 0) then
      CALL modsched (isched,maxsched,ngrids,nxtnest,nndtrat,nsubs)
   endif
   
   ! Send timestep schedule and timesteps to nodes only if they have changed.
   !   Need to do it in first timestep regardless...
   ntsend=0
   if(istp == 1 .or. nndtflg > 0) ntsend=1
   CALL master_putdtsched (isendflg,isendlite,isendmean,isendboth,ntsend)

!---------------------------------------------------------------------
! Bypass the timestep schedule, then update the main time variable.
!         do npass=1,nsubs
!         enddo

! Wait for cpu time, wallclock time, and Courant numbers
! cflxy and cflz from nodes.

   CALL master_getcflcpu ()

! Do the following calls to get the domain max W for the convergence
! forcing code used to initiate idealized convection.
! Wait for max vertical velocity from nodes, compute the domain max
! and then send this back to the nodes. Sub-domain max is computed
! within routine cfll in core/modsched.f90.
   CALL master_getvertvel ()
   CALL master_putvertvel ()

! Wait for whole subdomains from nodes
   if(isendflg.eq.1) then
      CALL master_getall ()
      if(iprntstmt>=1)print*,'calling par_ready',99999
      CALL par_ready (nmachs,machnum,99999)
   endif

   if(isendlite.eq.1) then
      CALL master_getanl ('LITE')
      if(iprntstmt>=1)print*,'calling par_ready',99998
      CALL par_ready (nmachs,machnum,99998)
   endif

   if(isendmean.eq.1) then
      CALL master_getanl ('MEAN')
      if(iprntstmt>=1)print*,'calling par_ready',99997
      CALL par_ready (nmachs,machnum,99997)
   endif

   if(isendboth.eq.1) then
      CALL master_getanl ('BOTH')
      if(iprntstmt>=1)print*,'calling par_ready',99996
      CALL par_ready (nmachs,machnum,99996)
   endif

!----------------------------------------------------------------------------
! Get timing info.

   wtime2 = walltime(wtime_start)
   CALL timing (2,T2)

!---------------------------------------------------------------------------
! Update main time variable by a long timestep.

   time=begtime+dtlongn(1)

!---------------------------------------------------------------------------
! Print timestep info and see if we write output now.

   timeh = TIME / 3600.
      PRINT 200, ISTP,TIME,timeh,T2-T1,wtime2-wtime1
200        FORMAT(' Master step-',I5,'  Sim time(sec)=',F9.1  &
           ,'  (hr)=',F6.2,'  CPU(sec)=',F6.2,'  Wall(sec)=',F6.2)

      pcpu=0
      pwall=0
      do n=1,nmachs
         pcpu=pcpu+ptimes(n,1)
         pwall=pwall+ptimes(n,2)
      enddo

    ! If printing is enabled, add in the individual node times
    if(iprntstmt>=1 .and. print_msg) then
      print 201,(ptimes(n,1),n=1,nmachs)
      print 202,(ptimes(n,2),n=1,nmachs)
201       format(' Node---CPU(sec)=',2000f8.3)
202       format(' Node--wall(sec)=',2000f8.3)
      print*,'----------------------------------------------'
    endif

      CALL rams_output ()
enddo

wtime_tot = walltime(wtime_start)
print '(//,a,f10.0)'  &
     ,' -----Total elapsed time: ',wtime_tot

return
END SUBROUTINE par_model
