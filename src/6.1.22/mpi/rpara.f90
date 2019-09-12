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

Module rpara

use grid_dims

implicit none

!---------------------------------------------------------------------------
integer                                         :: mainnum,nmachs,iparallel
integer, dimension(maxmach)                     :: machnum,nbuff_nest1  &
                                                  ,newbuff_nest1
integer, dimension(maxmach,maxgrds)             :: nxbeg,nxend,nybeg  &
                                                  ,nyend,nxbegc,nxendc  &
                                                  ,nybegc,nyendc,ixoff  &
                                                  ,iyoff,npxy,ibcflg  &
                                                  ,nestflg,ifeednode  &
                                                  ,ixb,ixe,iyb,iye
integer, dimension(maxmach,maxgrds,4)           :: nextnode
integer, dimension(5,7,maxgrds,maxmach,maxmach) :: inode_paths_master
integer, dimension(6,maxgrds,maxmach,maxmach)   :: iget_paths_master
integer, dimension(2,maxmach,maxmach)           :: lbc_buffs
real, dimension(maxmach)                        :: hperf
real, dimension(maxmach,maxgrds)                :: perf
real, dimension(maxmach,2)                      :: ptimes
!---------------------------------------------------------------------------

END MODULE rpara
