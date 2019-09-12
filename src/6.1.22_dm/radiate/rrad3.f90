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

Module rrad3

implicit none

!---------------------------------------------------------------------------
integer                         :: nrad,narad,nsolb,nb,ng,jday
integer, parameter              :: mb=8,mg=3,mk=7,ncog=5,ncb=2,npartob=13  &
                                  ,npartg=7,namax=10
integer, dimension(mg,mb)       :: npsb
real                            :: solfac
real, dimension(mb)             :: nuum,ralcs,solar1,solar0,a0,a1,a2,a3  &
                                  ,wlenlo,wlenhi
real, dimension(150)            :: exptabc
real, dimension(mg,mb)          :: prf,trf,ulim
real, dimension(mg,mk,mb)       :: wght,xp,alpha,beta
real, dimension(ncog,mb,npartob) :: ocoef
real, dimension(ncb,mb,npartob) :: bcoef
real, dimension(ncog,mb,npartg) :: gcoef
!---------------------------------------------------------------------------

END MODULE rrad3
