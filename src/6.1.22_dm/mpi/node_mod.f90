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
integer :: mainnum,nmachs,my_mpi_num,not_a_node,my_rams_num
integer, dimension(maxmach)            :: machnum,nbuff_nest1,newbuff_nest1
integer, dimension(maxmach,maxgrds)    :: ixb,ixe,iyb,iye
integer, dimension(maxmach,maxgrds,4)  :: nextnode
real, dimension(maxmach)               :: hperf
real, dimension(maxmach,maxgrds)       :: perf
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
integer :: mxp,myp,mzp,ia,iz,ja,jz,i0,j0,ibcon  &
          ,ia_1,ia_2,ia_3,ia1,ia2,ia3,ja_1,ja_2,ja_3,ja1,ja2,ja3  &
          ,iz_1,iz_2,iz_3,iz1,iz2,iz3,jz_1,jz_2,jz_3,jz1,jz2,jz3  &
          ,izu,jzv,mxyzp,mxysp,mxyp
!---------------------------------------------------------------------------
integer, dimension(maxgrds) :: mmxp,mmyp,mmzp,mmxyzp,mmxysp,mmxyp,mmxyzbp & 
      ,mia,miz,mja,mjz,mi0,mj0,mibcon,mxbeg,mxend,mybeg,myend &
      ,fd_sw_num,fd_nw_num,fd_se_num,fd_ne_num
!---------------------------------------------------------------------------
type io_descrip
  integer :: xblock
  integer :: yblock
  integer :: xoff
  integer :: yoff
end type io_descrip

integer, dimension(maxgrds) :: file_xchunk, file_ychunk
type (io_descrip), dimension(maxgrds) :: mem_read,mem_write,file_read,file_write
!---------------------------------------------------------------------------
integer, dimension(5,7,maxgrds,maxmach) :: ipaths
integer, dimension(6,maxgrds,maxmach)   :: iget_paths
!---------------------------------------------------------------------------
integer                       :: newbuff_feed,nbuff_feed,newbuff_nest,nbuff_nest
!---------------------------------------------------------------------------
integer, dimension(maxmach) :: irecv_req,isend_req
!---------------------------------------------------------------------------

type lbc_buff_type
   real, allocatable :: lbc_send_buff(:),lbc_recv_buff(:)
   integer :: nsend,nrecv
end type

type (lbc_buff_type) :: node_buffs(maxmach)

END MODULE node_mod
