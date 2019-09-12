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

use node_mod

implicit none

integer :: i,numarg,iargc,bad
character(len=strl1) :: nl_fname,arg,cargs(0:maxargs)

print*,'STARTING RAMS'
! argument defaults
nl_fname='RAMSIN'

! read arguments
do i=1,maxargs
 cargs(i)=''
enddo

numarg=iargc()

do i=0,numarg
   CALL ugetarg (i,arg)
   cargs(i)=trim(arg)//char(0) !Null terminate the string to pass to C
enddo

bad = 0
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
    nl_fname=cargs(2)(1:len_trim(cargs(2))-1)
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

! par_init_fortran() will return
!
!  my_mpi_num -> 0
!  nmachs -> 1
!
! for both the case where rams is called directly from the command
! line and the case where rams is called using mpirun -np 1 (one
! process). For these cases, rams will run in "sequential" mode.
! Otherwise, rams will run in "parallel" mode.
CALL par_init_fortran (numarg,cargs,len(cargs),my_mpi_num,nmachs)

! set global vars in node_mod
!
! MPI always numbers nodes from 0 to n-1
! RAMS numbers nodes from 1 to n
! Therefore, make the RAMS node (my_rams_num) equal to the MPI node (my_mpi_num) plus one
!
! Arbitrarily choose node 0 to be the "main" node
! which will take on the following extra duties:
!     1) Read and distribute the namelist groups
!     2) Decompose the domain, and distribute the results
!     3) Printing messages
!
mainnum = 1        ! rams node number 1, which is MPI node number 0
my_rams_num = my_mpi_num + 1

not_a_node = -1
do i = 1, nmachs
  machnum(i) = i - 1  ! MPI process numbers go from 0 to n-1
enddo

! run the model
CALL rams_model (nl_fname)

! MPI bookend call to par_init_fortran
CALL par_exit ()

END PROGRAM main
