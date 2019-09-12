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

Module io_params

use grid_dims

implicit none

character(len=32) :: lite_vars(maxlite)
character(len=strl1) :: hfilin,afilepref

integer :: ipast_sfc
!-------------------------------------------------------------------------------
integer :: ioutput,iclobber,nlite_vars
real    :: frqstate(maxgrds),avgtim,frqlite,frqmean,frqboth,frqst_keep  
!-------------------------------------------------------------------------------
integer, dimension(maxgrds) :: itoptflg,isstflg,ivegtflg,isoilflg  &
                              ,ndviflg,itopsflg,iz0flg
real                        :: z0fact
real, dimension(maxgrds)    :: z0max,toptenh,toptwvl
!-------------------------------------------------------------------------------
character(len=strl1), dimension(maxgrds) :: itoptfn,isstfn,ivegtfn,isoilfn  &
                                        ,ndvifn
!-------------------------------------------------------------------------------
character(len=strl1) :: sstfpfx,sfcfiles,topfiles
character(len=strl1), dimension(maxsstfiles,maxgrds) :: fnames_sst
character(len=14), dimension(maxsstfiles,maxgrds) :: itotdate_sst
!-------------------------------------------------------------------------------
integer                             :: iupdsst,isstcyclic,isstcycdata
integer,dimension(maxgrds)          :: nsstfiles,isstflp,isstflf
real,dimension(maxgrds)             :: ssttime1,ssttime2
!-------------------------------------------------------------------------------
character(len=strl1) :: ndvifpfx
character(len=strl1), dimension(maxndvifiles,maxgrds) :: fnames_ndvi
character(len=14), dimension(maxndvifiles,maxgrds) :: itotdate_ndvi
!-------------------------------------------------------------------------------
integer                             :: iupdndvi,indvicyclic,indvicycdata
integer,dimension(maxgrds)          :: nndvifiles,indviflp,indviflf
real,dimension(maxgrds)             :: ndvitime1,ndvitime2

END MODULE io_params
