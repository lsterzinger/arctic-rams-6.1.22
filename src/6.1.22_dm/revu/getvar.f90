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

integer Function rams_getvar (string,itype,ngrd,a,flnm)

use an_header
use hdf5_utils
use rcommons

implicit none

real :: a(*)
integer :: itype,ngrd,ni,npts,iword
character(len=*) :: flnm,string
character(len=1) :: cgrid
character(len=strl1) :: flng,errmsg
logical :: there

! First see if data file for this grid/time exists...

write(cgrid,'(i1)') ngrd
flng=trim(flnm)//'-g'//cgrid//'.h5'
inquire(file=flng,exist=there)

if(.not.there) then
   errmsg='File not found - '//flng
   CALL error_mess (errmsg)
   rams_getvar=2
   return
endif

! Now search table for actual variable

do ni=1,nvbtab
   if(string == anal_table(ni)%string.and.ngrd == anal_table(ni)%ngrid) then
   
      npts=anal_table(ni)%nvalues
      itype=anal_table(ni)%idim_type
      iword=anal_table(ni)%npointer

         CALL shdf5_open (trim(flng),'R')
         !call shdf5_info (string,ndims,idims)
         !npts_chk=product(idims(1:ndims)) 
         !if (npts /= npts_chk) then
         !   print*,'No. of points in anal table and in hdf5 file do not match.'
         !   print*,'   anal field:',string
         !   print*,'   anal table:',npts
         !   print*,'   hdf5 file :',idims(1:ndims)
         !   print*,'   hdf5 file :',npts_chk
         !   stop
         !endif
         CALL shdf5_irec (trim(string),rvara=a)
         CALL shdf5_close ()

      rams_getvar=0
      ifound=ifound+1
      return

   endif
enddo

errmsg='Variable not available in this run - '//string
CALL error_mess (errmsg)
rams_getvar=1
ierr_getvar=1

return
END FUNCTION rams_getvar
