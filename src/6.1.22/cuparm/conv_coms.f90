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

Module conv_coms

implicit none

integer, parameter :: nkp=100

integer ::      kcon,klcl,klfc,ketl,kct,igo,kmt  &
               ,icprtfl,icpltfl
               
real    ::      zmid,cdzmin,dzlow,dzhigh,plcl,tlcl,dzlcl,zlcl,garea  &
               ,wconmin,contim,preff,envshr,supply,cptime,cprecip

real, dimension(nkp) :: ucon,vcon,wcon,thtcon ,rvcon,prcon,picon,tmpcon  &
               ,dncon,zcon,zzcon  &
               ,upe,vpe,wpe,ze,te,the,pe,rte,pke,rhoe,thve,zc  &
               ,rve,thee,qvct1,qvct2,qvct3,qvct4  &
               ,vheat,vmois,vmdry,frcon,ftcon,tcon,rcon &
               ,theu,rsu,thu,tu,thd,wtd,thcon,rtcon

END MODULE conv_coms
