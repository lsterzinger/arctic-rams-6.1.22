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

Module rcommons

use grid_dims

implicit none

!**************************************************
! For rainit.f90 routines
!**************************************************
character(len=strl1), save :: fnames(maxfiles)
integer, save :: nfgrids(maxfiles),ifdates(maxfiles),iftimes(maxfiles)
real, save :: ftimes(maxfiles),startutc
!**************************************************
! FROM REVUIN
!**************************************************
character(len=strl1) :: anpref,revpref
character(len=8) :: anatype
character(len=20) :: xvar,yvar,zvar,tvar
character(len=strl1) :: revuvar(maxrevu)
integer :: igrid,iztran
!**************************************************

character(len=1) :: ftran
integer,save :: maxmem
integer :: ierr_getvar,ifound,ivar_type,iany

integer :: nii,nib,nie,niinc,njj,njb,nje,njinc,nnb,nne,nninc &
          ,itbeg,itstep,itend,ixbeg,ixstep,ixend  &
          ,iybeg,iystep,iyend,izbeg,izstep,izend

integer, dimension(nxpmax,maxgrds) :: ipm,jpm,nrz,kpm
integer, parameter :: nplevs=37
integer, save :: iplevs(nplevs)

data iplevs/1000,  975,  950,  925,  900, &
             875,  850,  825,  800,  775, &
             750,  725,  700,  675,  650, &
             625,  600,  575,  550,  525, &
             500,  475,  450,  425,  400, &
             375,  350,  325,  300,  275, &
             250,  225,  200,  175,  150, &
             125,  100/

END MODULE rcommons
