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

Module ref_sounding

use grid_dims

implicit none

   integer :: iref,jref
   real :: topref
   real, dimension(nzpmax,maxgrds) :: u01dn,v01dn,pi01dn,th01dn,dn01dn,rt01dn
   
   integer                    :: ipsflg,itsflg,irtsflg,iusflg,nsndg
   real, dimension(maxsndg)   :: us,vs,ts,thds,ps,hs,rts

END MODULE ref_sounding
