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

Module node_mod

use grid_dims

implicit none

!---------------------------------------------------------------------------
integer :: mxp,myp,mzp,ia,iz,ja,jz,i0,j0,master_num,mchnum,ibcon,ipara  &
          ,ia_1,ia_2,ia_3,ia1,ia2,ia3,ja_1,ja_2,ja_3,ja1,ja2,ja3  &
          ,iz_1,iz_2,iz_3,iz1,iz2,iz3,jz_1,jz_2,jz_3,jz1,jz2,jz3  &
          ,izu,jzv,mynum
!---------------------------------------------------------------------------
integer, target, dimension(maxgrds) :: mmxp,mmyp,mmzp
integer, dimension(maxgrds) :: mia,miz,mja,mjz  &
                              ,mi0,mj0,mibcon,mnestflg,mfeednode
!---------------------------------------------------------------------------
integer                                 :: nmachs
integer, dimension(maxmach)             :: machs
integer, dimension(maxmach,maxgrds)     :: nodemxp,nodemyp,nodemzp  &
                                          ,nodeia,nodeiz,nodeja,nodejz  &
                                          ,nodei0,nodej0,nodeibcon  &
                                          ,nodenestflg,nodefeednode
integer, dimension(maxmach,maxgrds,4)   :: nodeconn
integer, dimension(5,7,maxgrds,maxmach) :: ipaths
integer, dimension(6,maxgrds,maxmach)   :: iget_paths
!---------------------------------------------------------------------------
integer                       :: newbuff_feed,nbuff_feed,newbuff_nest  &
                                ,nbuff_nest
!---------------------------------------------------------------------------
integer, dimension(maxmach) :: irecv_req,isend_req
!---------------------------------------------------------------------------

type lbc_buffs
   real, allocatable :: lbc_send_buff(:),lbc_recv_buff(:)
   integer :: nsend,nrecv
end type

type (lbc_buffs) :: node_buffs(maxmach)

END MODULE node_mod
