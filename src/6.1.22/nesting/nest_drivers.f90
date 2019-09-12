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

Subroutine nstbdriv ()

use mem_tend
use var_tables
use mem_basic
use mem_nestb
use node_mod
use mem_grid

implicit none

integer :: n
real :: tymeinvv,tymeinvs

tymeinvv = 1.0 / (dtlv * 0.5 * float(nndtrat(ngrid)+2-isstp))
tymeinvs = 1.0 / (dtlt * float(nndtrat(ngrid)+1-isstp))

!         print*, 'calling nstbtnd for u'
!         print*, 'jdim,ibcon,ia,iz,ja,jz'
!         print*,  jdim,ibcon,ia,iz,ja,jz


CALL nstbtnd (mzp,mxp,myp,ia,iz,ja,jz,ibcon  &
   ,basic_g(ngrid)%up(1,1,1),tend%ut(1)  &
   ,nbounds(ngrid)%bux(1,1,1),nbounds(ngrid)%buy(1,1,1)  &
   ,nbounds(ngrid)%buz(1,1,1)  &
   ,'u',tymeinvv,nstbot,nsttop,jdim)

CALL nstbtnd (mzp,mxp,myp,ia,iz,ja,jz,ibcon  &
   ,basic_g(ngrid)%vp(1,1,1),tend%vt(1)  &
   ,nbounds(ngrid)%bvx(1,1,1),nbounds(ngrid)%bvy(1,1,1)  &
   ,nbounds(ngrid)%bvz(1,1,1)  &
   ,'v',tymeinvv,nstbot,nsttop,jdim)

CALL nstbtnd (mzp,mxp,myp,ia,iz,ja,jz,ibcon  &
   ,basic_g(ngrid)%wp(1,1,1),tend%wt(1)  &
   ,nbounds(ngrid)%bwx(1,1,1),nbounds(ngrid)%bwy(1,1,1)  &
   ,nbounds(ngrid)%bwz(1,1,1)  &
   ,'w',tymeinvv,nstbot,nsttop,jdim)

CALL nstbtnd (mzp,mxp,myp,ia,iz,ja,jz,ibcon  &
   ,basic_g(ngrid)%pp(1,1,1),tend%pt(1)  &
   ,nbounds(ngrid)%bpx(1,1,1),nbounds(ngrid)%bpy(1,1,1)  &
   ,nbounds(ngrid)%bpz(1,1,1)  &
   ,'p',tymeinvv,nstbot,nsttop,jdim)

do n = 1,num_scalar(ngrid)
   CALL nstbtnd (mzp,mxp,myp,ia,iz,ja,jz,ibcon  &
      ,scalar_tab(n,ngrid)%var_p,scalar_tab(n,ngrid)%var_t  &
      ,nbounds(ngrid)%bsx(1,1,1,n),nbounds(ngrid)%bsy(1,1,1,n)  &
      ,nbounds(ngrid)%bsz(1,1,1,n)  &
      ,'t',tymeinvs,nstbot,nsttop,jdim)
enddo

return
END SUBROUTINE nstbdriv

!##############################################################################
Subroutine prgintrp (m1c,m2c,m3c,ifm,initflg)

use var_tables
use mem_scratch
use mem_basic
use mem_grid
use mem_nestb

implicit none

integer :: m1c,m2c,m3c,ifm,initflg,icm,nf,nc


!     Temporarily fill VT2DA with interpolated topography from coarser grid

icm = nxtnest(ifm)
if (icm .eq. 0) return
if (initflg .eq. 1) then
   CALL fillscr (1,maxnxp,maxnyp,1,m2c,m3c,1,1  &
      ,scratch%scr1(1),grid_g(icm)%topt(1,1))
   CALL eintp (scratch%scr1(1),scratch%scr2(1)  &
      ,1,maxnxp,maxnyp,1,nnxp(ifm),nnyp(ifm),ifm,2,'t',0,0)
   CALL fillvar (1,maxnxp,maxnyp,1,nnxp(ifm),nnyp(ifm),1,1  &
      ,scratch%scr2(1),scratch%vt2da(1))
endif

!     Interpolate atmospheric variables

CALL fmint3 (m1c,m2c,m3c,nnzp(ifm),nnxp(ifm),nnyp(ifm),maxnzp,maxnxp,maxnyp  &
   ,ifm,icm,nnstbot(ifm),nnsttop(ifm),jdim,initflg,1,1,'u'  &
   ,basic_g(icm)%uc(1,1,1),basic_g(ifm)%uc(1,1,1)  &
   ,basic_g(icm)%dn0u(1,1,1),basic_g(ifm)%dn0u(1,1,1)  &
   ,scratch%scr1(1),scratch%scr2(1)  &
   ,grid_g(ifm)%topt(1,1),scratch%vt2da(1)  &
   ,nbounds(ifm)%bux(1,1,1),nbounds(ifm)%buy(1,1,1)  &
   ,nbounds(ifm)%buz(1,1,1))

if (jdim .eq. 1 .or. icorflg .eq. 1)  &
CALL fmint3 (m1c,m2c,m3c,nnzp(ifm),nnxp(ifm),nnyp(ifm),maxnzp,maxnxp,maxnyp  &
   ,ifm,icm,nnstbot(ifm),nnsttop(ifm),jdim,initflg,1,1,'v'  &
   ,basic_g(icm)%vc(1,1,1),basic_g(ifm)%vc(1,1,1)  &
   ,basic_g(icm)%dn0v(1,1,1),basic_g(ifm)%dn0v(1,1,1)  &
   ,scratch%scr1(1),scratch%scr2(1)  &
   ,grid_g(ifm)%topt(1,1),scratch%vt2da(1)  &
   ,nbounds(ifm)%bvx(1,1,1),nbounds(ifm)%bvy(1,1,1)  &
   ,nbounds(ifm)%bvz(1,1,1))
CALL fmint3 (m1c,m2c,m3c,nnzp(ifm),nnxp(ifm),nnyp(ifm),maxnzp,maxnxp,maxnyp  &
   ,ifm,icm,nnstbot(ifm),nnsttop(ifm),jdim,initflg,1,1,'w'  &
   ,basic_g(icm)%wc(1,1,1),basic_g(ifm)%wc(1,1,1)  &
   ,basic_g(icm)%dn0(1,1,1),basic_g(ifm)%dn0(1,1,1)  &
   ,scratch%scr1(1),scratch%scr2(1)  &
   ,grid_g(ifm)%topt(1,1),scratch%vt2da(1)  &
   ,nbounds(ifm)%bwx(1,1,1),nbounds(ifm)%bwy(1,1,1)  &
   ,nbounds(ifm)%bwz(1,1,1))
CALL fmint3 (m1c,m2c,m3c,nnzp(ifm),nnxp(ifm),nnyp(ifm),maxnzp,maxnxp,maxnyp  &
   ,ifm,icm,nnstbot(ifm),nnsttop(ifm),jdim,initflg,0,1,'t'  &
   ,basic_g(icm)%pc(1,1,1),basic_g(ifm)%pc(1,1,1)  &
   ,basic_g(icm)%dn0(1,1,1),basic_g(ifm)%dn0(1,1,1)  &
   ,scratch%scr1(1),scratch%scr2(1)  &
   ,grid_g(ifm)%topt(1,1),scratch%vt2da(1)  &
   ,nbounds(ifm)%bpx(1,1,1),nbounds(ifm)%bpy(1,1,1)  &
   ,nbounds(ifm)%bpz(1,1,1))

do nf = 1,num_scalar(ifm)
   do nc = 1,num_scalar(icm)
      if (scalar_tab(nf,ifm)%name == scalar_tab(nc,icm)%name) then
         CALL fmint3 (m1c,m2c,m3c,nnzp(ifm),nnxp(ifm),nnyp(ifm)  &
            ,maxnzp,maxnxp,maxnyp,ifm,icm  &
            ,nnstbot(ifm),nnsttop(ifm),jdim,initflg,1,1,'t'  &
            ,scalar_tab(nc,icm)%var_p &
            ,scalar_tab(nf,ifm)%var_p  &
            ,basic_g(icm)%dn0(1,1,1),basic_g(ifm)%dn0(1,1,1)  &
            ,scratch%scr1(1),scratch%scr2(1)  &
            ,grid_g(ifm)%topt(1,1),scratch%vt2da(1)  &
            ,nbounds(ifm)%bsx(1,1,1,nf),nbounds(ifm)%bsy(1,1,1,nf)  &
            ,nbounds(ifm)%bsz(1,1,1,nf))
      endif
   enddo
enddo

return
END SUBROUTINE prgintrp

!##############################################################################
Subroutine nstfeed (ifm,icm)

use var_tables
use mem_basic
use mem_scratch
use node_mod
use mem_grid

implicit none

integer :: ifm,icm

integer :: nf,nc
real, pointer :: scalarf, scalarc

!     feed back the finer mesh to the coarser mesh.

if (icm .eq. 0) return

CALL fdback (basic_g(icm)%uc   (1,1,1)  ,basic_g(ifm)%uc   (1,1,1)  &
           ,basic_g(icm)%dn0u (1,1,1)  ,basic_g(ifm)%dn0u (1,1,1)  &
   ,nnzp(icm),nnxp(icm),nnyp(icm),nnzp(ifm),nnxp(ifm),nnyp(ifm)    &
   ,ifm,'u',scratch%scr1(1))

if (nnstbot(icm) .eq. 1) CALL botset (nnzp(icm),nnxp(icm),nnyp(icm)  &
   ,basic_g(icm)%uc(1,1,1),'U')

if (nnsttop(icm) .eq. 1) CALL topset (nnzp(icm),nnxp(icm),nnyp(icm)  &
   ,basic_g(icm)%uc(1,1,1),'U')

CALL fdback (basic_g(icm)%vc   (1,1,1)  ,basic_g(ifm)%vc   (1,1,1)  &
           ,basic_g(icm)%dn0v (1,1,1)  ,basic_g(ifm)%dn0v (1,1,1)  &
   ,nnzp(icm),nnxp(icm),nnyp(icm),nnzp(ifm),nnxp(ifm),nnyp(ifm)    &
   ,ifm,'v',scratch%scr1(1))

if (nnstbot(icm) .eq. 1) CALL botset (nnzp(icm),nnxp(icm),nnyp(icm)  &
   ,basic_g(icm)%vc(1,1,1),'V')

if (nnsttop(icm) .eq. 1) CALL topset (nnzp(icm),nnxp(icm),nnyp(icm)  &
   ,basic_g(icm)%vc(1,1,1),'V')

CALL fdback (basic_g(icm)%wc  (1,1,1)  ,basic_g(ifm)%wc  (1,1,1)  &
           ,basic_g(icm)%dn0 (1,1,1)  ,basic_g(ifm)%dn0 (1,1,1)  &
   ,nnzp(icm),nnxp(icm),nnyp(icm),nnzp(ifm),nnxp(ifm),nnyp(ifm)  &
   ,ifm,'w',scratch%scr1(1))

CALL fdback (basic_g(icm)%pc  (1,1,1)  ,basic_g(ifm)%pc   (1,1,1)  &
           ,basic_g(icm)%dn0 (1,1,1)  ,basic_g(ifm)%dn0  (1,1,1)  &
   ,nnzp(icm),nnxp(icm),nnyp(icm),nnzp(ifm),nnxp(ifm),nnyp(ifm)   &
   ,ifm,'p',scratch%scr1(1))

if (nnstbot(icm) .eq. 1) CALL botset (nnzp(icm),nnxp(icm),nnyp(icm)  &
   ,basic_g(icm)%pc(1,1,1),'P')

if (nnsttop(icm) .eq. 1) CALL topset (nnzp(icm),nnxp(icm),nnyp(icm)  &
   ,basic_g(icm)%pc(1,1,1),'P')

do nf = 1,num_scalar(ifm)
   do nc = 1,num_scalar(icm)
      if (scalar_tab(nf,ifm)%name == scalar_tab(nc,icm)%name) then
         scalarc => scalar_tab(nc,icm)%var_p
         scalarf => scalar_tab(nf,ifm)%var_p

         CALL fdback (scalarc,scalarf  &
            ,basic_g(icm)%dn0(1,1,1),basic_g(ifm)%dn0(1,1,1)  &
            ,nnzp(icm),nnxp(icm),nnyp(icm),nnzp(ifm),nnxp(ifm)  &
            ,nnyp(ifm),ifm,'t',scratch%scr1(1))

         if (nnstbot(icm) .eq. 1) CALL botset (nnzp(icm),nnxp(icm),nnyp(icm)  &
            ,scalarc,'T')

         if (nnsttop(icm) .eq. 1) CALL topset (nnzp(icm),nnxp(icm),nnyp(icm)  &
            ,scalarc,'T')

      endif
   enddo
enddo

return
END SUBROUTINE nstfeed

