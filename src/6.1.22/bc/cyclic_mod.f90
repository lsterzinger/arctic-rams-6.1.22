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

Module cyclic_mod

implicit none

integer, save ::   &
    npts_cyc=0     & ! number of ij columns in icycpaths array
   ,nbuffsend_cyc  & ! size of send buffer
   ,nbuffrecv_cyc    ! size of receive buffer (for each node received from)

integer, save, dimension(2) ::         &
    ndns_cyc=0     & ! number of destination nodes this source node sends to
   ,nsns_cyc=0       ! number of source nodes that send to this destination node

integer, save, allocatable, dimension(:,:) ::  &
    ndn_cyc        & ! counter over destination nodes actually sent to
   ,nsn_cyc        & ! counter over source nodes actually received from
   ,msn_cyc        & ! list of source nodes that send to this destination node
   ,mdn_cyc        & ! list of destination node numbers sent to
   ,ijrecv_cyc     & ! list of ij columns received in buffer
   ,nijsendt_cyc   & ! number of t ij columns sent to each destination node
   ,nijsendu_cyc   & ! number of u ij columns sent to each destination node 
   ,nijsendv_cyc   & ! number of v ij columns sent to each destination node 
   ,nijrecvt_cyc   & ! number of t ij columns received from each source node
   ,nijrecvu_cyc   & ! number of u ij columns received from each source node
   ,nijrecvv_cyc   & ! number of v ij columns received from each source node
   ,buffsend_cyc   & ! send buffer
   ,buffrecv_cyc     ! receive buffer

integer, save, allocatable, dimension(:,:,:) ::  &
    ipathst_cyc    & ! cyclic paths array
   ,ipathsu_cyc    & ! cyclic paths array
   ,ipathsv_cyc    & ! cyclic paths array
   ,isend_req_cyc  & ! send request array
   ,irecv_req_cyc    ! receive request array

integer, save, allocatable, dimension(:,:,:,:) ::  &
    ijsendt_cyc    & ! list of ij columns sent to each destination node
   ,ijsendu_cyc    & ! list of ij columns sent to each destination node
   ,ijsendv_cyc      ! list of ij columns sent to each destination node

END MODULE cyclic_mod

