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

Subroutine diffuse ()

! +-----------------------------------------------------------------+
! \     this routine is the subdriver to compute tendencies due to  \
! \       subgrid-scale turbulence.                                 \
! +-----------------------------------------------------------------+

use mem_tend
use mem_basic
use var_tables
use mem_turb
use mem_grid
use mem_leaf
use mem_micro
use mem_scratch
use node_mod
use ke_coms
use micphys

implicit none

integer :: mxyzp,ind,n
real :: s1,s2,s3
real, pointer :: scalarp, scalart
integer :: i,j,k,ksf

mxyzp = mxp * myp * mzp

CALL strain (mzp,mxp,myp,ia,iz,ja,jz                          &
      ,ia_1,ja_1,iz1,jz1,jdim                                &
      ,basic_g(ngrid)%up (1,1,1) ,basic_g(ngrid)%vp (1,1,1)  &
      ,basic_g(ngrid)%wp (1,1,1) ,scratch%vt3da     (1)      &
      ,scratch%vt3db     (1)     ,scratch%vt3dc     (1)      &
      ,scratch%vt3dd     (1)     ,scratch%vt3de     (1)      &
      ,scratch%vt3df     (1)     ,scratch%vt3dg     (1)      &
      ,scratch%vt3dh     (1)     ,scratch%vt3di     (1)      &
      ,scratch%vt3dn     (1)     ,scratch%scr2      (1)      &
      ,idiffk(ngrid))

if (level <= 1) CALL azero (mxyzp,scratch%vt3dp(1))
if (level >= 2) CALL ae1 (mxyzp,scratch%vt3dp(1),micro_g(ngrid)%rcp(1,1,1))

CALL bruvais (mzp,mxp,myp,ia,iz,ja,jz                          &
   ,basic_g(ngrid)%theta (1,1,1) ,basic_g(ngrid)%rtp (1,1,1)  &
   ,basic_g(ngrid)%rv    (1,1,1) ,scratch%vt3dp(1)            &
   ,basic_g(ngrid)%pp    (1,1,1) ,basic_g(ngrid)%pi0 (1,1,1)  &
   ,scratch%vt3dj        (1)     ,grid_g(ngrid)%rtgt (1,1))

if (idiffk(ngrid) <= 3) then
   CALL mxdefm (mzp,mxp,myp,ia,iz,ja,jz                      &
      ,scratch%vt3dh      (1)     ,scratch%vt3di      (1)    &
      ,scratch%vt3dj      (1)     ,scratch%vt3dk      (1)    &
      ,scratch%scr1       (1)     ,scratch%scr2       (1)    &
      ,basic_g(ngrid)%dn0 (1,1,1) ,grid_g(ngrid)%rtgt (1,1)  &
      ,grid_g(ngrid)%dxt  (1,1))
endif

if (idiffk(ngrid) == 1) then
   CALL tkemy (mzp,mxp,myp,ia,iz,ja,jz,jdim  &
      ,turb_g(ngrid)%tkep   (1,1,1) ,tend%tket            (1)      &
      ,scratch%vt3dh        (1)     ,scratch%vt3di        (1)      &
      ,scratch%vt3dj        (1)     ,scratch%scr1         (1)      &
      ,grid_g(ngrid)%rtgt   (1,1)   ,basic_g(ngrid)%theta (1,1,1)  &
      ,basic_g(ngrid)%dn0   (1,1,1) ,basic_g(ngrid)%up    (1,1,1)  &
      ,basic_g(ngrid)%vp    (1,1,1) ,basic_g(ngrid)%wp    (1,1,1)  &
      ,turb_g(ngrid)%sflux_u(1,1)   ,turb_g(ngrid)%sflux_v(1,1)    &
      ,turb_g(ngrid)%sflux_w(1,1)   ,turb_g(ngrid)%sflux_t(1,1),vctr34)
endif

if (idiffk(ngrid) == 4) then
   CALL mxtked (mzp,mxp,myp,ia,iz,ja,jz  &
      ,turb_g(ngrid)%tkep   (1,1,1) ,tend%tket            (1)      &
      ,scratch%vt3da        (1)     ,scratch%vt3dc        (1)      &
      ,scratch%vt3dh        (1)     ,scratch%vt3dj        (1)      &
      ,scratch%scr1         (1)     ,scratch%scr2         (1)      &
      ,grid_g(ngrid)%dxt    (1,1)   ,grid_g(ngrid)%rtgt   (1,1))
endif

!_STC............................................................
!_STC Call to routine tkescl for E-l closure
!_STC (S. Trini Castelli)
!_STC............................................................
if (idiffk(ngrid) == 5) then
   CALL tkescl (mzp,mxp,myp,npatch,ia,iz,ja,jz  &
         ,turb_g(ngrid)%tkep(1,1,1),tend%tket(1)  &
         ,turb_g(ngrid)%epsp(1,1,1)  &
         ,scratch%vt3da(1),scratch%vt3dc(1)  &
         ,scratch%vt3dh(1),scratch%vt3di(1)  &
         ,scratch%vt3dj(1),scratch%scr1(1)  &
         ,scratch%scr2(1) ,grid_g(ngrid)%rtgt(1,1)  &
         ,scratch%vt3dd(1),scratch%vt3de(1),grid_g(ngrid)%dxt(1,1)  &
         ,leaf_g(ngrid)%ustar(1,1,1),leaf_g(ngrid)%patch_area(1,1,1) &
         ,basic_g(ngrid)%dn0(1,1,1)  )
endif
!_STC............................................................
!_STC Call to routine tkeeps for E-eps closure
!_STC (S. Trini Castelli)
!_STC............................................................
if (idiffk(ngrid) == 6) then
   CALL tkeeps (mzp,mxp,myp,npatch,ia,iz,ja,jz  &
         ,turb_g(ngrid)%tkep(1,1,1),tend%tket(1)  &
         ,turb_g(ngrid)%epsp(1,1,1),tend%epst(1)  &
         ,scratch%vt3da(1),scratch%vt3dc(1)  &
         ,scratch%vt3dh(1),scratch%vt3di(1)  &
         ,scratch%vt3dj(1),scratch%scr1(1)  &
         ,scratch%scr2(1) ,grid_g(ngrid)%rtgt(1,1)  &
         ,leaf_g(ngrid)%ustar(1,1,1),leaf_g(ngrid)%patch_area(1,1,1) &
         ,basic_g(ngrid)%dn0(1,1,1)  )
endif
!_STC..................................................
!_STC    Note: from routines TKESCL, TKEEPS :
!_STC           VT3DI=Ke
!_STC           SCR1=Km
!_STC           VT3DH = Kh
!_STC           SCR2 = SCR1 = Km
!_STC..................................................
!_STC............................................................

 CALL klbnd (mzp,mxp,myp,ibcon,jdim  &
    ,scratch%scr1 (1),basic_g(ngrid)%dn0(1,1,1))
 CALL klbnd (mzp,mxp,myp,ibcon,jdim  &
    ,scratch%scr2 (1),basic_g(ngrid)%dn0(1,1,1))
 CALL klbnd (mzp,mxp,myp,ibcon,jdim  &
    ,scratch%vt3dh(1),basic_g(ngrid)%dn0(1,1,1))
!_STC ....... boundary conditions even on Ke diffusion coefficient
if(idiffk(ngrid) ==  5 .or. idiffk(ngrid) == 6) &
 CALL klbnd (mzp,mxp,myp,ibcon,jdim  &
    ,scratch%vt3di(1),basic_g(ngrid)%dn0(1,1,1))

!bob  swap new hkm, vkm, and vkh with past time level:  lagged K's have
!bob  internal lateral boundary values from neighboring nodes

ind = 0
do j = 1,mmyp(ngrid)
   do i = 1,mmxp(ngrid)
      do k = 1,mmzp(ngrid)
         ind = ind + 1
         s1 = scratch%scr2(ind)
         s2 = scratch%scr1(ind)
         s3 = scratch%vt3dh(ind)
         scratch%scr2(ind) = turb_g(ngrid)%hkm(k,i,j)
         scratch%scr1(ind) = turb_g(ngrid)%vkm(k,i,j)
         scratch%vt3dh(ind) = turb_g(ngrid)%vkh(k,i,j)
!! also for vt3di = K(tke) ?????    22 March 02
!!         scratch%vt3di(ind) = turb_g(ngrid)%vke(k,i,j)
         turb_g(ngrid)%hkm(k,i,j) = s1
         turb_g(ngrid)%vkm(k,i,j) = s2
         turb_g(ngrid)%vkh(k,i,j) = s3
      enddo
   enddo
enddo

CALL diffvel (mzp,mxp,myp,ia,iz,ja,jz,jdim                         &
      ,iz1,jz1,izu,jzv,idiffk(ngrid)                               &
      ,basic_g(ngrid)%up    (1,1,1) ,basic_g(ngrid)%vp    (1,1,1)  &
      ,basic_g(ngrid)%wp    (1,1,1) ,tend%ut              (1)      &
      ,tend%vt              (1)     ,tend%wt              (1)      &
      ,scratch%vt3da        (1)     ,scratch%vt3db        (1)      &
      ,scratch%vt3dc        (1)     ,scratch%vt3dd        (1)      &
      ,scratch%vt3de        (1)     ,scratch%vt3df        (1)      &
      ,scratch%vt3dg        (1)     ,scratch%vt3dj        (1)      &
      ,scratch%vt3dk        (1)     ,scratch%vt3dl        (1)      &
      ,scratch%vt3dm        (1)     ,scratch%vt3dn        (1)      &
      ,scratch%vt3do        (1)     ,grid_g(ngrid)%rtgu   (1,1)    &
      ,grid_g(ngrid)%rtgv   (1,1)   ,grid_g(ngrid)%rtgt   (1,1)    &
      ,turb_g(ngrid)%sflux_u(1,1)   ,turb_g(ngrid)%sflux_v(1,1)    &
      ,turb_g(ngrid)%sflux_w(1,1)   ,basic_g(ngrid)%dn0   (1,1,1)  &
      ,basic_g(ngrid)%dn0u  (1,1,1) ,basic_g(ngrid)%dn0v  (1,1,1)  &
      ,scratch%scr1         (1)     ,scratch%scr2         (1)      )

! Convert momentum K's to scalar K's, if necessary

if (idiffk(ngrid) <= 3) then
   do ind = 1,mxyzp
      scratch%scr2(ind) = scratch%scr2(ind) * xkhkm(ngrid)
   enddo
elseif (idiffk(ngrid) == 4) then
   do ind = 1,mxyzp
      scratch%vt3di(ind) = 2. * scratch%scr1(ind)
   enddo
endif

do n = 1,num_scalar(ngrid)

   scalarp => scalar_tab(n,ngrid)%var_p
   scalart => scalar_tab(n,ngrid)%var_t

   CALL azero (mxp*myp,scratch%vt2da(1))
   if (nstbot == 1) then
      if (scalar_tab(n,ngrid)%name == 'THP') then
         CALL atob (mxp*myp,turb_g(ngrid)%sflux_t(1,1),scratch%vt2da(1))
      elseif (scalar_tab(n,ngrid)%name == 'RTP') then
         CALL atob (mxp*myp,turb_g(ngrid)%sflux_r(1,1),scratch%vt2da(1))
      endif
   endif

! 3/10/01 - Define ksf below, the "K scalar flag", to let routine diffsclr
! know which vertical K is being passed to it.  If diffsclr sees that it's
! a different K from the previous one, diffsclr will re-compute the tridiff
! matrix coefficients.  In order to use vertical scalar K's other than
! vt3dh and vt3di, use ksf = 3, ksf = 4, etc. for each different K.

!_STC..................................................
!_STC Corrections to account for the new idiffk options
!_STC for E-l and E-eps closure. Isotropy hypothesis.
!_STC (S. Trini Castelli)
!_STC..................................................

   ! Need to adjust the vkkh and hkkh arguments of diffsclr() according
   ! to which scheme we are using. These are the last two arguments:
   !
   !    diffsclr(..., vkkh, hkkh)
   !
   if (scalar_tab(n,ngrid)%name == 'TKEP') then
      ksf = 1   
      if (idiffk(ngrid) >= 4) then
        ! vkkh <- scratch%vt3di(1)
        ! hkkh <- scratch%vt3di(1)
        CALL diffsclr (mzp,mxp,myp,ia,iz,ja,jz,jdim,n,ksf         &
         ,scalarp,scalart            ,scratch%vt3da(1)            &
         ,scratch%vt3db      (1)     ,scratch%vt3df      (1)      &
         ,scratch%vt3dg      (1)     ,scratch%vt3dj      (1)      &
         ,scratch%vt3dk      (1)     ,scratch%vt3do      (1)      &
         ,scratch%vt3dc      (1)     ,scratch%vt3dd      (1)      &
         ,scratch%vt3dl      (1)     ,scratch%vt3dm      (1)      &
         ,scratch%vt2db      (1)     ,grid_g(ngrid)%rtgt (1,1)    &
         ,scratch%vt2da      (1)     ,basic_g(ngrid)%dn0 (1,1,1)  &
         ,scratch%vt3di      (1)     ,scratch%vt3di      (1)      )
      else
        ! vkkh <- scratch%vt3di(1)
        ! hkkh <- scratch%scr2(1)
        CALL diffsclr (mzp,mxp,myp,ia,iz,ja,jz,jdim,n,ksf         &
         ,scalarp,scalart            ,scratch%vt3da(1)            &
         ,scratch%vt3db      (1)     ,scratch%vt3df      (1)      &
         ,scratch%vt3dg      (1)     ,scratch%vt3dj      (1)      &
         ,scratch%vt3dk      (1)     ,scratch%vt3do      (1)      &
         ,scratch%vt3dc      (1)     ,scratch%vt3dd      (1)      &
         ,scratch%vt3dl      (1)     ,scratch%vt3dm      (1)      &
         ,scratch%vt2db      (1)     ,grid_g(ngrid)%rtgt (1,1)    &
         ,scratch%vt2da      (1)     ,basic_g(ngrid)%dn0 (1,1,1)  &
         ,scratch%vt3di      (1)     ,scratch%scr2       (1)      )
      endif
   elseif (scalar_tab(n,ngrid)%name == 'EPSP') then
      ksf = 3
      ! Convert Ktke to Keps; it will be converted back after use below
      if (idiffk(ngrid) >= 4) then
        ! vkkh <- scratch%vt3di(1)
        ! hkkh <- scratch%vt3di(1)
        CALL ae1t0 (mxyzp,scratch%vt3di,scratch%vt3di,ALF_EPS/ALF_TKE)
        CALL ae1t0 (mxyzp,scratch%vt3di,scratch%vt3di,ALF_EPS/ALF_TKE)

        CALL diffsclr (mzp,mxp,myp,ia,iz,ja,jz,jdim,n,ksf         &
         ,scalarp,scalart            ,scratch%vt3da(1)            &
         ,scratch%vt3db      (1)     ,scratch%vt3df      (1)      &
         ,scratch%vt3dg      (1)     ,scratch%vt3dj      (1)      &
         ,scratch%vt3dk      (1)     ,scratch%vt3do      (1)      &
         ,scratch%vt3dc      (1)     ,scratch%vt3dd      (1)      &
         ,scratch%vt3dl      (1)     ,scratch%vt3dm      (1)      &
         ,scratch%vt2db      (1)     ,grid_g(ngrid)%rtgt (1,1)    &
         ,scratch%vt2da      (1)     ,basic_g(ngrid)%dn0 (1,1,1)  &
         ,scratch%vt3di      (1)     ,scratch%vt3di      (1)      )

        CALL ae1t0 (mxyzp,scratch%vt3di,scratch%vt3di,ALF_TKE/ALF_EPS)
        CALL ae1t0 (mxyzp,scratch%vt3di,scratch%vt3di,ALF_TKE/ALF_EPS)
      else
        ! vkkh <- scratch%vt3di(1)
        ! hkkh <- scratch%scr2(1)
        CALL ae1t0 (mxyzp,scratch%vt3di,scratch%vt3di,ALF_EPS/ALF_TKE)
        CALL ae1t0 (mxyzp,scratch%scr2,scratch%scr2,ALF_EPS/ALF_TKE)

        CALL diffsclr (mzp,mxp,myp,ia,iz,ja,jz,jdim,n,ksf         &
         ,scalarp,scalart            ,scratch%vt3da(1)            &
         ,scratch%vt3db      (1)     ,scratch%vt3df      (1)      &
         ,scratch%vt3dg      (1)     ,scratch%vt3dj      (1)      &
         ,scratch%vt3dk      (1)     ,scratch%vt3do      (1)      &
         ,scratch%vt3dc      (1)     ,scratch%vt3dd      (1)      &
         ,scratch%vt3dl      (1)     ,scratch%vt3dm      (1)      &
         ,scratch%vt2db      (1)     ,grid_g(ngrid)%rtgt (1,1)    &
         ,scratch%vt2da      (1)     ,basic_g(ngrid)%dn0 (1,1,1)  &
         ,scratch%vt3di      (1)     ,scratch%scr2       (1)      )

        CALL ae1t0 (mxyzp,scratch%vt3di,scratch%vt3di,ALF_TKE/ALF_EPS)
        CALL ae1t0 (mxyzp,scratch%scr2,scratch%scr2,ALF_TKE/ALF_EPS)
      endif
   else
      ksf = 2
      if (idiffk(ngrid) >= 4) then
        ! vkkh <-  scratch%vt3dh(1)
        ! hkkh <-  scratch%vt3dh(1)
        CALL diffsclr (mzp,mxp,myp,ia,iz,ja,jz,jdim,n,ksf         &
         ,scalarp,scalart            ,scratch%vt3da(1)            &
         ,scratch%vt3db      (1)     ,scratch%vt3df      (1)      &
         ,scratch%vt3dg      (1)     ,scratch%vt3dj      (1)      &
         ,scratch%vt3dk      (1)     ,scratch%vt3do      (1)      &
         ,scratch%vt3dc      (1)     ,scratch%vt3dd      (1)      &
         ,scratch%vt3dl      (1)     ,scratch%vt3dm      (1)      &
         ,scratch%vt2db      (1)     ,grid_g(ngrid)%rtgt (1,1)    &
         ,scratch%vt2da      (1)     ,basic_g(ngrid)%dn0 (1,1,1)  &
         ,scratch%vt3dh      (1)     ,scratch%vt3dh      (1)      )
      else
        ! vkkh <-  scratch%vt3dh(1)
        ! hkkh <-  scratch%scr2(1)
        CALL diffsclr (mzp,mxp,myp,ia,iz,ja,jz,jdim,n,ksf         &
         ,scalarp,scalart            ,scratch%vt3da(1)            &
         ,scratch%vt3db      (1)     ,scratch%vt3df      (1)      &
         ,scratch%vt3dg      (1)     ,scratch%vt3dj      (1)      &
         ,scratch%vt3dk      (1)     ,scratch%vt3do      (1)      &
         ,scratch%vt3dc      (1)     ,scratch%vt3dd      (1)      &
         ,scratch%vt3dl      (1)     ,scratch%vt3dm      (1)      &
         ,scratch%vt2db      (1)     ,grid_g(ngrid)%rtgt (1,1)    &
         ,scratch%vt2da      (1)     ,basic_g(ngrid)%dn0 (1,1,1)  &
         ,scratch%vt3dh      (1)     ,scratch%scr2       (1)      )
      endif
   endif

enddo

return
END SUBROUTINE diffuse

!##############################################################################
Subroutine strain (m1,m2,m3,ia,iz,ja,jz,ia_1,ja_1,iz1,jz1  &
   ,jd,up,vp,wp,vt3da,vt3db,vt3dc,vt3dd,vt3de  &
   ,vt3df,vt3dg,vt3dh,vt3di,vt3dn,scr2,idiffk)

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,ia_1,ja_1,iz1,jz1,jd,idiffk,i,j,k
real, dimension(m1,m2,m3) :: up,vp,wp,vt3da,vt3db,vt3dc,vt3dd,vt3de,vt3df &
   ,vt3dg,vt3dh,vt3di,vt3dn,scr2

CALL grad (m1,m2,m3,ia,iz1,ja,jz,up,vt3da,'XDIR','UPNT')
CALL grad (m1,m2,m3,ia_1,iz,ja_1,jz,vp,vt3db,'XDIR','VPNT')
CALL grad (m1,m2,m3,ia_1,iz,ja,jz,wp,vt3df,'XDIR','WPNT')

CALL grad (m1,m2,m3,ia_1,iz,ja_1,jz,up,vt3dn,'YDIR','UPNT')
CALL grad (m1,m2,m3,ia,iz,ja,jz1,vp,vt3dc,'YDIR','VPNT')
CALL grad (m1,m2,m3,ia,iz,ja_1,jz,wp,vt3dg,'YDIR','WPNT')

CALL grad (m1,m2,m3,ia_1,iz,ja,jz,up,vt3dd,'ZDIR','UPNT')
CALL grad (m1,m2,m3,ia,iz,ja_1,jz,vp,vt3de,'ZDIR','VPNT')
if(idiffk.ge.3)then
   CALL grad (m1,m2,m3,ia,iz,ja,jz,wp,scr2,'ZDIR','WPNT')
endif

if (idiffk .le. 2) then
   do j = ja,jz
      do i = ia,iz
         do k = 2,m1-1
            vt3dh(k,i,j) =2. * (vt3da(k,i,j) * vt3da(k,i,j)  &
               + vt3dc(k,i,j) * vt3dc(k,i,j))  &
               + .0625 * (vt3db(k,i,j) + vt3db(k,i-1,j)  &
               + vt3db(k,i,j-jd) + vt3db(k,i-1,j-jd)  &
               + vt3dn(k,i,j) + vt3dn(k,i-1,j)  &
               + vt3dn(k,i,j-jd) + vt3dn(k,i-1,j-jd)) ** 2
            vt3di(k,i,j) = .0625 * ((vt3dd(k,i,j) + vt3dd(k-1,i,j)  &
               + vt3dd(k,i-1,j) + vt3dd(k-1,i-1,j)) ** 2  &
               + (vt3de(k,i,j) + vt3de(k-1,i,j)  &
               + vt3de(k,i,j-jd) + vt3de(k-1,i,j-jd)) ** 2)
         enddo
      enddo
   enddo
else
   do j = ja,jz
      do i = ia,iz
         do k = 2,m1-1
            vt3da(k,i,j) = 2. * vt3da(k,i,j)
            vt3dc(k,i,j) = 2. * vt3dc(k,i,j)
            scr2(k,i,j) = 2. * scr2(k,i,j)
            vt3db(k,i,j) = vt3db(k,i,j) + vt3dn(k,i,j)
            vt3dn(k,i,j) = vt3db(k,i,j)
            vt3dd(k,i,j) = vt3dd(k,i,j) + vt3df(k,i,j)
            vt3de(k,i,j) = vt3de(k,i,j) + vt3dg(k,i,j)
            vt3di(k,i,j) = 0.333333  &
               * (vt3da(k,i,j) + vt3dc(k,i,j) + scr2(k,i,j))
         enddo
      enddo

      do k = 2,m1-1
         vt3da(k,iz1,j) = 2. * vt3da(k,iz1,j)
         vt3db(k,ia_1,j) = vt3db(k,ia_1,j) + vt3dn(k,ia_1,j)
         vt3dn(k,ia_1,j) = vt3db(k,ia_1,j)
         vt3dd(k,ia_1,j) = vt3dd(k,ia_1,j) + vt3df(k,ia_1,j)
      enddo
   enddo

   do i = ia_1,iz
      do k = 2,m1-1
         vt3dc(k,i,jz1) = 2. * vt3dc(k,i,jz1)
         vt3db(k,i,ja_1) = vt3db(k,i,ja_1) + vt3dn(k,i,ja_1)
         vt3dn(k,i,ja_1) = vt3db(k,i,ja_1)
         vt3de(k,i,ja_1) = vt3de(k,i,ja_1) + vt3dg(k,i,ja_1)
      enddo
   enddo

   do j = ja,jz
      do i = ia,iz
         do k = 2,m1-1
            vt3dh(k,i,j) = .5 * (  &
               (vt3da(k,i,j) - vt3di(k,i,j)) ** 2  &
               + (vt3dc(k,i,j) - vt3di(k,i,j)) ** 2  &
               + ( scr2(k,i,j) - vt3di(k,i,j)) ** 2)  &
               + .0625 * ((vt3db(k,i,j) + vt3db(k,i-1,j)  &
               + vt3db(k,i,j-jd) + vt3db(k,i-1,j-jd)) ** 2  &
               + (vt3dd(k,i,j) + vt3dd(k,i-1,j)  &
               + vt3dd(k-1,i,j) + vt3dd(k-1,i-1,j)) ** 2  &
               + (vt3de(k,i,j) + vt3de(k-1,i,j)  &
               + vt3de(k,i,j-jd) + vt3de(k-1,i,j-jd)) ** 2)
            vt3di(k,i,j) = vt3dh(k,i,j)
         enddo
      enddo
   enddo
endif

return
END SUBROUTINE strain

!##############################################################################
Subroutine bruvais (m1,m2,m3,ia,iz,ja,jz,theta,rtp,rv,rcp,pp,pi0,en2,rtgt)

use mem_scratch
use micphys
use mem_grid
use rconstants

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz
real, dimension(m1,m2,m3) :: theta,rtp,rcp,rv,pp,pi0,en2
real, dimension(m2,m3) :: rtgt
integer :: i,j,k,iweten,iwdiffk,ki
real :: c1,c2,c3,ci1,ci2,ci3,rvlsi,rvii
real, dimension(nzpmax) :: pi,temp,prt,rvls,rc

!     calculate brunt-vaisalla frequency squared (en2)

iweten = 1

iwdiffk = 0
c1 = alvl / rgas
c2 = ep * alvl ** 2 / (cp * rgas)
c3 = alvl / cp
ci1 = alvi / rgas
ci2 = ep * alvi ** 2 / (cp * rgas)
ci3 = alvi / cp
ki=m1+1
!     calculate potential temperature profile

do j = ja,jz
   do i = ia,iz

      do k = 1,m1
         vctr11(k) = theta(k,i,j)
         vctr12(k) = theta(k,i,j)
         vctr32(k) = 0.
      enddo
      if (level .ge. 1) then
         do k = 1,m1
            vctr12(k) = vctr11(k) * (1. + .61 * rv(k,i,j))
            vctr32(k) = (rtp(k,i,j) - rv(k,i,j))
         enddo
      endif

!     check for saturation if level is 2 or greater.

      if (level .ge. 2 .and. iweten .eq. 1) then
         do k = 1,m1
            pi(k) = (pp(k,i,j) + pi0(k,i,j)) / cp
            temp(k) = theta(k,i,j) * pi(k)
            prt(k) = p00 * pi(k) ** cpor
         enddo
         CALL mrsl (m1,prt(1),temp(1),rvls(1))
         do k = 1,m1
            vctr2(k) = c1
            vctr3(k) = c2
            vctr4(k) = c3
         enddo
         ki = m1 + 1

!     if any ice phase microphysics are activated ....

         if ((ipris .ge. 1 .or. isnow .ge. 1 .or.  &
             igraup .ge. 1 .or. iaggr .ge. 1 .or.  &
             ihail .ge. 1) .and. level .eq. 3) then

!            find level of -20 c.  assume ice saturation above this
!              level.

            do k = 1,m1
               if (temp(k) .le. 253.16) then
                  ki = k
                  go to 10
               endif
            enddo
            ki = m1 + 1
10               continue

            CALL mrsi (m1-ki+1,prt(ki),temp(ki),rvls(ki))
            if(ki > 1) CALL mrsi (1,prt(ki-1),temp(ki-1),rvlsi)
            do k = ki,m1
               vctr2(k) = ci1
               vctr3(k) = ci2
               vctr4(k) = ci3
            enddo
         endif

         if (level .eq. 3) then
            do k = 1,m1
               rc(k) = rcp(k,i,j)
            enddo
         else
            do k = 1,m1
               rc(k) = max(rv(k,i,j) / rvls(k) - .999,0.)
            enddo
         endif

      endif

      do k = 2,m1-1
         vctr1(k) = g / ((zt(k+1) - zt(k-1)) * rtgt(i,j))
      enddo
      if (level .ge. 2 .and. iweten .eq. 1) then
         do k = 2,m1-1
            if (rc(k) .gt. 0.) then
               rvii = rvls(k-1)
               if (k .eq. ki) rvii = rvlsi
               en2(k,i,j) = vctr1(k) * (  &
                  (1. + vctr2(k) * rvls(k) / temp(k))  &
                  / (1. + vctr3(k) * rvls(k) / temp(k) ** 2)  &
                  * ((vctr11(k+1) - vctr11(k-1)) / vctr11(k)  &
                  + vctr4(k) / temp(k) * (rvls(k+1) - rvii))  &
                  - (rtp(k+1,i,j) - rtp(k-1,i,j)))
            else
               en2(k,i,j) = vctr1(k)*((vctr12(k+1)-vctr12(k-1))  &
                  / vctr12(k) - (vctr32(k+1) - vctr32(k-1)))
            endif

         enddo
      else
         do k = 2,m1-1
            en2(k,i,j) = vctr1(k) * ((vctr12(k+1)-vctr12(k-1))  &
               / vctr12(k) - (vctr32(k+1) - vctr32(k-1)))
         enddo
      endif
      en2(1,i,j) = en2(2,i,j)
      en2(nzp,i,j)=en2(nz,i,j)

   enddo
enddo

return
END SUBROUTINE bruvais

!##############################################################################
Subroutine mxdefm (m1,m2,m3,ia,iz,ja,jz  &
   ,vt3dh,vt3di,vt3dj,vt3dk,scr1,scr2,dn0,rtgt,dxt)

!     +-------------------------------------------------------------+
!     \   this routine calculates the mixing coefficients with a    \
!     \     smagorinsky-type deformational based k with an optional \
!     \     unstable brunt-vaisala enhancement and an optional      \
!     \     richardson number modification.                         \
!     +-------------------------------------------------------------+

use mem_scratch
use mem_grid
use mem_turb
use rconstants

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz
real, dimension(m1,m2,m3) :: vt3dh,vt3di,vt3dj,vt3dk,scr1,scr2,dn0
real, dimension(m2,m3) :: rtgt,dxt
integer :: i,j,k,irich,ienfl
real :: csx2,sq300,enfl,rchmax,c1,c2,c3,c4,akm,rmin,rmax

irich = 1
ienfl = 1

csx2 = csx(ngrid) * csx(ngrid)
sq300 = 90000.
if (idiffk(ngrid) .eq. 2 .or. idiffk(ngrid) .eq. 3) then
   rmin = -100.
   rmax = 1. / zkhkm(ngrid)
   do j = ja,jz
      do i = ia,iz
         do k = 2,m1-1
            vt3dk(k,i,j) = max(min(vt3dj(k,i,j)  &
               / max(vt3di(k,i,j),1.e-15),rmax),rmin)
         enddo
      enddo
   enddo
   enfl = float(ienfl)
   rchmax = 1.0 + 9.0 * float(irich)
   do k = 2, m1
      vctr1(k) = csz(ngrid) * (zm(k) - zm(k-1))
      vctr2(k) = vctr1(k) * vctr1(k)
   enddo
endif

if (idiffk(ngrid) .eq. 1) then
   do j = ja,jz
      do i = ia,iz
         c2 = 1.0 / (dxt(i,j) * dxt(i,j))
         c3 = csx2 * c2
         akm = akmin(ngrid) * 0.075 * c2 ** (0.666667)
         do k = 2,m1-1
            scr2(k,i,j) = dn0(k,i,j)  &
               * max(akm,c3*sqrt(vt3dh(k,i,j)))
         enddo
      enddo
   enddo
elseif (idiffk(ngrid) .eq. 2) then
   do j = ja,jz
      do i = ia,iz
         c1 = rtgt(i,j) * rtgt(i,j)
         c2 = 1.0 / (dxt(i,j) * dxt(i,j))
         c3 = csx2 * c2
         akm = akmin(ngrid) * 0.075 * c2 ** (0.666667)
         c4 = vonk * vonk * c1
         do k = 2,m1-1

! old csz*dz len  scr1(k,i,j) = dn0(k,i,j) * c1 * vctr2(k)
! asymptotic vertical scale length from bjorn with modifications:
! c3 is (csx * dx)^2, c1*vctr2(k) is (csz * dz)^2, sq300 is the square
! of 300 meters (used as a limit for horizontal grid spacing influence
! on vertical scale length), ambda is (asymptotic_vertical_length_scale)^2,
! and vkz2 is (vonk * height_above_surface)^2.

!Saleeby: THESE FOLLOWING 3 LINES OF CODE NEED SOME WORK.
!            ambda = max(c1 * vctr2(k),min(sq300,c3))
!            vkz2 = c4 * zt(k) * zt(k)
!            scr1(k,i,j) = dn0(k,i,j) * vkz2 / (vkz2 / ambda + 1)  &

             scr1(k,i,j) = dn0(k,i,j) * c1 * vctr2(k) & !(Saleeby: added this back)
               * (sqrt(vt3di(k,i,j))  &
               + enfl * sqrt(max(0.,-vt3dj(k,i,j))))*min(rchmax  &
               ,sqrt(max(0.,(1.-zkhkm(ngrid)*vt3dk(k,i,j)))))

            scr2(k,i,j) = dn0(k,i,j)  &
               * max(akm,c3*sqrt(vt3dh(k,i,j)))
            vt3dh(k,i,j) = scr1(k,i,j) * zkhkm(ngrid)

         enddo
      enddo
   enddo

elseif (idiffk(ngrid) .eq. 3) then
   do j = ja,jz
      do i = ia,iz
         c1 = rtgt(i,j) * rtgt(i,j)
         do k = 2,m1-1
            scr1(k,i,j) = dn0(k,i,j) * c1 * vctr2(k)  &
               * (sqrt(vt3dh(k,i,j))  &
               + enfl * sqrt(max(0.,-vt3dj(k,i,j))))*min(rchmax  &
               ,sqrt(max(0.,(1.-zkhkm(ngrid)*vt3dk(k,i,j)))))
            scr2(k,i,j) = scr1(k,i,j)
            vt3dh(k,i,j) = scr1(k,i,j) * zkhkm(ngrid)
         enddo
      enddo
   enddo
endif

return
END SUBROUTINE mxdefm

!##############################################################################
Subroutine klbnd (m1,m2,m3,ibcon,jd,akay,dn0)

implicit none

integer :: m1,m2,m3,ibcon,jd
real, dimension(m1,m2,m3) :: akay,dn0
integer ::i,j,k

!     boundary conditions on a mixing coefficient

do j = 1,m3
   do i = 1,m2
      akay(1,i,j) = akay(2,i,j) * dn0(1,i,j) / dn0(2,i,j)
      akay(m1,i,j) = akay(m1-1,i,j) * dn0(m1,i,j)  &
         / dn0(m1-1,i,j)
   enddo
enddo

if (iand(ibcon,1) .ne. 0) then
   do j = 1,m3
      do k = 1,m1
         akay(k,1,j) = akay(k,2,j)
      enddo
   enddo
endif

if (iand(ibcon,2) .ne. 0) then
   do j = 1,m3
      do k = 1,m1
         akay(k,m2,j) = akay(k,m2-1,j)
      enddo
   enddo
endif

if (jd .eq. 1) then
   if (iand(ibcon,4) .ne. 0) then
      do i = 1,m2
         do k = 1,m1
            akay(k,i,1) = akay(k,i,2)
         enddo
      enddo
   endif

   if (iand(ibcon,8) .ne. 0) then
      do i = 1,m2
         do k = 1,m1
            akay(k,i,m3) = akay(k,i,m3-1)
         enddo
      enddo
   endif
endif

return
END SUBROUTINE klbnd
