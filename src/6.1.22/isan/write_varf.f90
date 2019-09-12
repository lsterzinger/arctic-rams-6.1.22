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

Subroutine write_varf ()

use isan_coms
use mem_grid
use io_params
use hdf5_utils

implicit none

integer :: ng,ndims,idims(4),isan_file_ver
character(len=strl1) :: locfn
character(len=3) :: csuff


do ng=1,nigrids
   nxyzp=nnxp(ng)*nnyp(ng)*nnzp(ng)
   nxyp =nnxp(ng)*nnyp(ng)
   write(csuff,'(a1,i1)') 'g',ng

   CALL makefnam (locfn,varpfx,0,iyear,imonth,idate  &
                 ,ihour*100,'V',csuff,'h5 ')
 
   CALL shdf5_open (locfn,'W',iclobber)

   isan_file_ver=3
   ndims=1 ; idims(1)=1
   CALL shdf5_orec (ndims,idims,'version',ivars=isan_file_ver)
   CALL shdf5_orec (ndims,idims,'year'   ,ivars=iyear)     
   CALL shdf5_orec (ndims,idims,'month'  ,ivars=imonth)    
   CALL shdf5_orec (ndims,idims,'day'    ,ivars=idate)     
   CALL shdf5_orec (ndims,idims,'hour'   ,ivars=ihour)     
   CALL shdf5_orec (ndims,idims,'nx'     ,ivars=nnxp(ng))  
   CALL shdf5_orec (ndims,idims,'ny'     ,ivars=nnyp(ng))  
   CALL shdf5_orec (ndims,idims,'nz'     ,ivars=nnzp(ng))  
   CALL shdf5_orec (ndims,idims,'polelat',rvars=polelat)        
   CALL shdf5_orec (ndims,idims,'polelon',rvars=polelon)       
   CALL shdf5_orec (ndims,idims,'dx'     ,rvars=deltaxn(ng))         
   CALL shdf5_orec (ndims,idims,'dz'     ,rvars=deltazn(ng))     
   CALL shdf5_orec (ndims,idims,'dzrat'  ,rvars=dzrat)     
   CALL shdf5_orec (ndims,idims,'dzmax'  ,rvars=dzmax)

   ndims=3 ; idims(1)=nnxp(ng); idims(2)=nnyp(ng); idims(3)=nnzp(ng)
   CALL rearrange (nnzp(ng),nnxp(ng),nnyp(ng),is_grids(ng)%rr_u,rr_scr3)
   CALL shdf5_orec (ndims,idims,'UP',rvara=rr_scr3)
   CALL rearrange (nnzp(ng),nnxp(ng),nnyp(ng),is_grids(ng)%rr_v,rr_scr3)
   CALL shdf5_orec (ndims,idims,'VP',rvara=rr_scr3)
   CALL rearrange (nnzp(ng),nnxp(ng),nnyp(ng),is_grids(ng)%rr_p,rr_scr3)
   CALL shdf5_orec (ndims,idims,'PI',rvara=rr_scr3)
   CALL rearrange (nnzp(ng),nnxp(ng),nnyp(ng),is_grids(ng)%rr_t,rr_scr3)
   CALL shdf5_orec (ndims,idims,'THETA',rvara=rr_scr3)
   CALL rearrange (nnzp(ng),nnxp(ng),nnyp(ng),is_grids(ng)%rr_r,rr_scr3)
   CALL shdf5_orec (ndims,idims,'RV',rvara=rr_scr3)

   CALL rearrange (nnzp(ng),nnxp(ng),nnyp(ng),is_grids(ng)%rr_cond,rr_scr3)
   CALL shdf5_orec (ndims,idims,'COND',rvara=rr_scr3)
   
   ndims=2 ; idims(1)=nnxp(ng); idims(2)=nnyp(ng)
   CALL vmissw (is_grids(ng)%rr_soilmoist1(1,1),nxyp,rr_vt2da(1),1E30,-9999.)
   CALL shdf5_orec (ndims,idims,'SOILMOIST1',rvara=rr_vt2da)
   CALL vmissw (is_grids(ng)%rr_soilmoist2(1,1),nxyp,rr_vt2da(1),1E30,-9999.)
   CALL shdf5_orec (ndims,idims,'SOILMOIST2',rvara=rr_vt2da)
   CALL vmissw (is_grids(ng)%rr_soiltemp1(1,1),nxyp,rr_vt2da(1),1E30,-9999.)
   CALL shdf5_orec (ndims,idims,'SOILTEMP1',rvara=rr_vt2da)
   CALL vmissw (is_grids(ng)%rr_soiltemp2(1,1),nxyp,rr_vt2da(1),1E30,-9999.)
   CALL shdf5_orec (ndims,idims,'SOILTEMP2',rvara=rr_vt2da)
   CALL vmissw (is_grids(ng)%rr_snowmass(1,1),nxyp,rr_vt2da(1),1E30,-9999.)
   CALL shdf5_orec (ndims,idims,'SNOWMASS',rvara=rr_vt2da)
   CALL vmissw (is_grids(ng)%rr_snowdepth(1,1),nxyp,rr_vt2da(1),1E30,-9999.)
   CALL shdf5_orec (ndims,idims,'SNOWDEPTH',rvara=rr_vt2da)

   CALL shdf5_close ()

enddo
CALL makefnam (locfn,varpfx,0,iyear,imonth,idate  &
              ,ihour*100,'V','$','tag')
CALL rams_f_open (2,locfn,'FORMATTED','REPLACE','WRITE',iclobber)
write(2,*) nigrids
close(2)

return
END SUBROUTINE write_varf
