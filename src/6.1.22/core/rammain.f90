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

Program main

use grid_dims
use mem_grid, only:print_msg

implicit none

integer :: machsize=0,machnum=0
integer :: taskids(maxmach),nn,icall,nproc,ipara

integer :: i,numarg,iargc,bad=0
character(len=strl1) :: name_name,arg,cargs(0:maxargs)

! argument defaults
name_name='RAMSIN'

! read arguments
do i=1,maxargs
 cargs(i)=''
enddo

numarg=iargc()

do i=0,numarg
   CALL ugetarg (i,arg)
   cargs(i)=trim(arg)//char(0) !Null terminate the string to pass to C
enddo

if(numarg > 0)then
 if(numarg == 2)then
   if(cargs(1)(1:1)=='-' .and. cargs(1)(1:2)=='-f') then
    if(len_trim(cargs(2)).gt.80) then
       print*,'Max filename length = 80 characters'
       bad=1
    endif
    ! We just null terminated the file name (cargs(2)) so it can be
    ! passed to C (par_init_fortran) as a proper C string.
    ! Trim off the trailing null for the FORTRAN string nl_fname
    name_name=cargs(2)(1:len_trim(cargs(2))-1)
   else
    bad=1
   endif
 else
   bad=1
 endif
endif

if (bad > 0) then
   print*,'RAMS usage: ''exec name'' '
   print*,'  [-f ''Namelist file''] '
   stop 'bad command line arguments'
endif

numarg=numarg+1 !Total number of arguments on command line

CALL par_init_fortran (numarg,cargs,len(cargs),machnum,machsize)

icall=0
if (machnum .ne. 0) icall=1
nproc=machsize-1
do nn=1,nproc
   taskids(nn)=nn
enddo

ipara=0
if (machsize > 0 ) then
   ipara = 1
   nproc=machsize-1
   if(nproc == 0) ipara=0
endif

!print*,'+-----------------------------------------------------'
!print*,'! RAMS input namelist file: ',trim(name_name)
!print*,'! RAMS call (master=0,node=1):', icall
!print*,'! Parallel info,machnum,machsize,ipara:', machnum,machsize,ipara
!print*,'+-----------------------------------------------------'

if (icall == 0) then

   print_msg = .true.
   CALL rams_master (ipara,nproc,taskids,machnum,name_name)

else if (icall == 1) then

   CALL rams_node ()

endif

END PROGRAM main
