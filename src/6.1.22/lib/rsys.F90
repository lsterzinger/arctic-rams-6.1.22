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

!-------------------------------------------------------------------
!                 MACHINE/OS DEPENDENT ROUTINES
!  This contains short utility routines that are not
!  of the FORTRAN 77 standard and may differ from machine to machine
!  or OS to OS. These include bit manipulation, I/O, JCL calls, and 
!  vector funcs.
!-------------------------------------------------------------------

Subroutine usystem (arg)

!MAY NEED MACHINE DEPENDENT ARGUMENTS OR USE STATEMENTS FOR THIS ROUTINE

implicit none

      character(len=*) :: arg
      
      !Routine to get command line argument

#if defined (PC_LINUX1)
      CALL system (trim(arg))
#else
   print*,"You specified machine/OS other than PC_LINUX1"
   print*,"You need to modify filelist.F90 to add your machine/OS"
   stop
#endif

return
END SUBROUTINE usystem

!##############################################################################
Subroutine ugetarg (i,arg)

!MAY NEED MACHINE DEPENDENT ARGUMENTS OR USE STATEMENTS FOR THIS ROUTINE

implicit none

      integer :: i
      character(len=*) :: arg
      
      !Routine to get command line argument

#if defined (PC_LINUX1)
      CALL getarg (i,arg)
#else
   print*,"You specified machine/OS other than PC_LINUX1"
   print*,"You need to modify filelist.F90 to add your machine/OS"
   stop
#endif

return
END SUBROUTINE ugetarg

!##############################################################################
Subroutine timing (icall,t1)

!MAY NEED MACHINE DEPENDENT ARGUMENTS OR USE STATEMENTS FOR THIS ROUTINE

implicit none

      integer :: icall
      real :: t1

      !Routine returns CPU time.  Called with icall=1 at beginning
      !of timestep, icall=2 at end of timestep.

      real :: et(2)
      real :: etime
      
      !If statement is here in case a machine/OS has a beginning and end 
      !timing command structure. Linux does NOT.

      if(icall.eq.1) then
      
#if defined (PC_LINUX1) 
        T1=ETIME(et)
#else
   print*,"You specified machine/OS other than PC_LINUX1"
   print*,"You need to modify filelist.F90 to add your machine/OS"
   stop
#endif

      ELSEIF(icall.EQ.2) THEN

#if defined (PC_LINUX1)
        T1=ETIME(et)
#else
   print*,"You specified machine/OS other than PC_LINUX1"
   print*,"You need to modify filelist.F90 to add your machine/OS"
   stop
#endif

      endif
      
return
END SUBROUTINE timing

!##############################################################################
Subroutine rams_filelist (fnames,file_prefix,nfile)

! Newer version that just uses ls and C to form unique filenames

use grid_dims

implicit none

integer :: nfile
character(len=*) :: fnames(*),file_prefix
character(len=strl1) :: file_pref,chome
character(len=1000000) :: fstring
integer :: iflag,iprelen

! this version uses nfile as flag for whether to stop if no files exist
! if nfile.ge.0, then stop

iflag=nfile

nfile = 0
!print*,'RAMS_filelist: Checking prefix: ',trim(file_prefix)

iprelen=len_trim(file_prefix)
if(iprelen == 0) iprelen=len(file_prefix)

! Process leading shell variables
if (file_prefix(1:5) == '$HOME') then
   CALL getenv ('HOME',chome)
   file_pref=trim(chome)//trim(file_prefix(6:))
else
   file_pref=trim(file_prefix)
endif
 
      
#if defined (PC_LINUX1)

   ! Call glob'bing routine
   fstring = ' '
   CALL c_listfile (trim(file_pref)//char(0),fstring)
   ! print*,'ffffs:',nfile,trim(fstring)
   
   ! Parse file name string 
   CALL tokenize2 (fstring,fnames,nfile,':')
   
   !do nc=1,nfile
   !   print*,'flist:',nc,trim(fnames(nc))
   !enddo

#else

   print*,"You specified machine/OS other than PC_LINUX1"
   print*,"You need to modify filelist.F90 to add your machine/OS"
   stop  

#endif

if (nfile == 0) then
   print *, 'No RAMS files for prefix:',trim(file_pref)
   if(iflag >= 0) stop 'RAMS_filelist-no_files'
endif

return
END SUBROUTINE rams_filelist
