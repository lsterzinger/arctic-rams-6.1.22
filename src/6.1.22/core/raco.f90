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

Subroutine acoustic ()
!--------------------------------------------------------------------
!                  Acoustic terms small time-step driver
!
!     This routine calls all the necessary routines to march the model
!     through the small timesteps.
!-----------------------------------------------------------------------
use mem_tend
use mem_grid
use mem_basic
use mem_scratch
use node_mod

implicit none

CALL acoust (mzp,mxp,myp,scratch%scr1(1),scratch%scr2(1),scratch%vt3da(1)&
           ,scratch%vt3db(1),scratch%vt3dc(1),scratch%vt3dd(1)  &
           ,scratch%vt3de(1),scratch%vt3df(1),scratch%vt3dg(1)  &
           ,scratch%vt3dh(1),scratch%vt2da(1)  &
           ,basic_g(ngrid)%dn0(1,1,1),basic_g(ngrid)%pi0(1,1,1)  &
           ,basic_g(ngrid)%th0(1,1,1),basic_g(ngrid)%up(1,1,1)  &
           ,basic_g(ngrid)%vp(1,1,1),basic_g(ngrid)%wp(1,1,1)  &
           ,basic_g(ngrid)%pp(1,1,1)  &
           ,tend%ut(1),tend%vt(1),tend%wt(1),tend%pt(1)  &
           ,grid_g(ngrid)%topt(1,1),grid_g(ngrid)%topu(1,1)  &
           ,grid_g(ngrid)%topv(1,1),grid_g(ngrid)%rtgt(1,1)  &
           ,grid_g(ngrid)%rtgu(1,1),grid_g(ngrid)%f13u(1,1)  &
           ,grid_g(ngrid)%dxu(1,1),grid_g(ngrid)%rtgv(1,1)  &
           ,grid_g(ngrid)%dyv(1,1)  &
           ,grid_g(ngrid)%f23v(1,1),grid_g(ngrid)%f13t(1,1)  &
           ,grid_g(ngrid)%f23t(1,1),grid_g(ngrid)%fmapui(1,1)  &
           ,grid_g(ngrid)%fmapvi(1,1),grid_g(ngrid)%dxt(1,1)  &
           ,grid_g(ngrid)%dyt(1,1),grid_g(ngrid)%fmapt(1,1))

return
END SUBROUTINE acoustic

!##############################################################################
Subroutine acoust (m1,m2,m3,scr1,scr2,vt3da,vt3db,vt3dc,vt3dd  &
                 ,vt3de,vt3df,vt3dg,vt3dh,vt2da  &
                 ,dn0,pi0,th0,up,vp,wp,pp,ut,vt,wt,pt  &
                 ,topt,topu,topv,rtgt,rtgu,f13u,dxu,rtgv,dyv  &
                 ,f23v,f13t,f23t,fmapui,fmapvi,dxt,dyt,fmapt)
!--------------------------------------------------------------------
!                  Acoustic terms small time-step driver
!
!     This routine calls all the necessary routines to march the model
!     through the small timesteps.
!-----------------------------------------------------------------------
use mem_grid
use mem_scratch
use node_mod

implicit none

integer :: m1,m2,m3
real, dimension(m1,m2,m3) ::dn0,pi0,th0,up,vp,wp,pp
real, dimension(m2,m3) ::   topt,topu,topv,rtgt,rtgu,f13u,dxu,rtgv,dyv  &
                           ,f23v,f13t,f23t,fmapui,fmapvi,dxt,dyt,fmapt
real, dimension(*) ::       scr1,scr2,vt3da,vt3db,vt3dc,vt3dd  &
                           ,vt3de,vt3df,vt3dg,vt3dh,vt2da,ut,vt,wt,pt

real :: a1da2

integer :: iter


do iter=1,nnacoust(ngrid)

!     Get coefficients for computations

dts = 2. * dtlt / nnacoust(ngrid)

if (iter .eq. 1)  &
   CALL coefz (mzp,mxp,myp,ia,iz,ja,jz  &
      ,vt3dc(1),vt3dd(1)  &
      ,vt3de(1),dn0(1,1,1)  &
      ,pi0(1,1,1),th0(1,1,1)  &
      ,rtgt(1,1),a1da2,vt3df(1)  &
      ,vt3dg(1),scr2(1)  &
      ,vctr1,vctr2)

CALL prdctu (mzp,mxp,myp,ia,izu,ja,jz  &
   ,up(1,1,1),ut(1)  &
   ,pp(1,1,1),vt3da(1)  &
   ,th0(1,1,1),vt3db(1)  &
   ,f13u(1,1),rtgu(1,1)  &
   ,rtgt(1,1),dxu(1,1)  &
   ,vt3dh(1),topu(1,1))

! Update u overlap regions
if (ipara .eq. 1) then
  CALL node_sendlbc (2)
  CALL node_getlbc (2)
endif

! Update u cyclic boundaries (u -> 2)
if (ngrid .eq. 1) CALL update_cyclic (2)

CALL prdctv (mzp,mxp,myp,ia,iz,ja,jzv  &
   ,vp(1,1,1),vt(1)  &
   ,pp(1,1,1),vt3da(1)  &
   ,th0(1,1,1),vt3db(1)  &
   ,f23v(1,1),rtgv(1,1)  &
   ,rtgt(1,1),dyv(1,1)  &
   ,vt3dh(1),topv(1,1))

! Update v overlap regions
if (ipara .eq. 1) then
  CALL node_sendlbc (3)
  CALL node_getlbc (3)
endif

! Update v cyclic boundaries (v -> 3)
if (ngrid .eq. 1) CALL update_cyclic (3)

CALL prdctw1 (mzp,mxp,myp,ia,iz,ja,jz  &
   ,wp(1,1,1),wt(1)  &
   ,pp(1,1,1),vt3dc(1)  &
   ,a1da2,vt3dh(1),rtgt(1,1)  &
   ,topt(1,1))

CALL prdctp1 (mzp,mxp,myp,ia,iz,ja,jz  &
   ,pp(1,1,1),up(1,1,1)  &
   ,vp(1,1,1),pi0(1,1,1)  &
   ,dn0(1,1,1),th0(1,1,1)  &
   ,pt(1),vt3da(1)  &
   ,vt3db(1),f13t(1,1)  &
   ,f23t(1,1),rtgt(1,1)  &
   ,rtgu(1,1),rtgv(1,1)  &
   ,vt2da(1),fmapui(1,1)  &
   ,fmapvi(1,1),dxt(1,1)  &
   ,dyt(1,1),fmapt(1,1))

CALL prdctw2 (mzp,mxp,myp,ia,iz,ja,jz  &
   ,wp(1,1,1),pp(1,1,1)  &
   ,vt3dc(1),vt3dg(1)  &
   ,scr1(1),scr2(1)  &
   ,rtgt(1,1),vt2da(1))

CALL prdctw3 (mzp,mxp,myp,ia,iz,ja,jz  &
   ,wp(1,1,1),scr1(1)  &
   ,vt3df(1),impl)

! Update w overlap regions
if (ipara .eq. 1) then
  CALL node_sendlbc (5)
  CALL node_getlbc (5)
endif

! Update v cyclic boundaries (w -> 5)
if (ngrid .eq. 1) CALL update_cyclic (5)

CALL prdctp2 (mzp,mxp,myp,ia,iz,ja,jz  &
   ,pp(1,1,1),wp(1,1,1)  &
   ,vt3dd(1),vt3de(1))

! Update p overlap regions
if (ipara .eq. 1) then
  CALL node_sendlbc (4)
  CALL node_getlbc (4)
endif

! Update p cyclic boundaries (p -> 4)
if (ngrid .eq. 1) CALL update_cyclic (4)

enddo

return
END SUBROUTINE acoust

!##############################################################################
Subroutine prdctu (m1,m2,m3,ia,iz,ja,jz  &
   ,up,ut,pp,vt3da,th0,dpdx,f13u,rtgu,rtgt,dxu,vt3dh,topu)

use mem_grid

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,i,j,k
real :: dxl
real, dimension(m1,m2,m3) :: up,ut,th0,dpdx,pp,vt3da,vt3dh
real, dimension(m2,m3) ::  f13u,rtgu,rtgt,dxu,topu

!     U prediction

CALL azero (m1*m2*m3,dpdx)

!     Calculate acoustic tendency (horizontal pressure gradient)

do j = ja,jz
   do i = ia,iz
      do k = 1,m1-1
         vt3da(k,i,j) = (pp(k,i,j) + pp(k+1,i,j)  &
            + pp(k,i+1,j) + pp(k+1,i+1,j)) * hw4(k)
      enddo
   enddo
enddo

do j = ja,jz
   do i = ia,iz
      dxl = dxu(i,j) / rtgu(i,j)
      do k = 2,m1-1
         dpdx(k,i,j) = -(th0(k,i,j) + th0(k,i+1,j)) * .5  &
            * ((pp(k,i+1,j) * rtgt(i+1,j) - pp(k,i,j) * rtgt(i,j)) * dxl  &
            + (vt3da(k,i,j) - vt3da(k-1,i,j)) * dzt(k) * f13u(i,j))
      enddo
   enddo
enddo

if (distim .ne. 0.) then
   CALL rayf (1,m1,m2,m3,ia,iz,ja,jz,up,th0,vt3dh,rtgu,topu)
endif

do j = 1,m3
   do i = 1,m2
      do k = 1,m1
         up(k,i,j) = up(k,i,j) + dts * (dpdx(k,i,j) + ut(k,i,j))
      enddo
   enddo
enddo

if (nstbot .eq. 1 .and. itopo .eq. 1)  &
     CALL botset (m1,m2,m3,up,'U')

return
END SUBROUTINE prdctu

!##############################################################################
Subroutine prdctv (m1,m2,m3,ia,iz,ja,jz  &
   ,vp,vt,pp,vt3da,th0,dpdy,f23v,rtgv,rtgt,dyv,vt3dh,topv)
   
use mem_grid

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,i,j,k

real :: dyl
real, dimension (m1,m2,m3) :: vp,vt,th0,dpdy,pp,vt3da,vt3dh
real, dimension(m2,m3) :: f23v,rtgv,rtgt,dyv,topv

!     V prediction

CALL azero (m1*m2*m3,dpdy)

if (jdim .eq. 1) then

!       calculate acoustic tendency (horizontal pressure gradient)

   do j = ja,jz
      do i = ia,iz
         do k = 1,m1-1
            vt3da(k,i,j) =(pp(k,i,j) + pp(k+1,i,j)  &
               + pp(k,i,j+1) + pp(k+1,i,j+1)) * hw4(k)
         enddo
      enddo
   enddo

   do j = ja,jz
      do i = ia,iz
         dyl = dyv(i,j) / rtgv(i,j)
         do k = 2,m1-1
            dpdy(k,i,j) = -(th0(k,i,j) + th0(k,i,j+1)) * .5  &
               * ((pp(k,i,j+1) * rtgt(i,j+1) - pp(k,i,j) * rtgt(i,j)) * dyl  &
               + (vt3da(k,i,j) - vt3da(k-1,i,j)) * dzt(k) * f23v(i,j))
         enddo
      enddo
   enddo

   if (distim .ne. 0.) then
      CALL rayf (2,m1,m2,m3,ia,iz,ja,jz,vp,th0,vt3dh,rtgv,topv)
   endif

   do j = 1,m3
      do i = 1,m2
         do k = 1,m1
            vp(k,i,j) = vp(k,i,j) + dts * (dpdy(k,i,j) + vt(k,i,j))
         enddo
      enddo
   enddo

   if (nstbot .eq. 1 .and. itopo .eq. 1)  &
        CALL botset (m1,m2,m3,vp,'V')

endif

return
END SUBROUTINE prdctv

!##############################################################################
Subroutine prdctw1 (m1,m2,m3,ia,iz,ja,jz  &
   ,wp,wt,pp,acoc,a1da2,vt3dh,rtgt,topt)

use mem_grid

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,i,j,k
real :: a1da2
real, dimension(m1,m2,m3) :: wp,wt,pp,acoc,vt3dh
real, dimension(m2,m3) :: rtgt,topt

!     First part of prediction at I,J point

!     Compute forward part of Crank-Nickelson scheme. This will be total
!     W prediction for explicit case.

if (distim .ne. 0.) then
   CALL rayf (3,m1,m2,m3,ia,iz,ja,jz,wp,wp,vt3dh,rtgt,topt)
endif

!      do j=ja,jz
!         do i=ia,iz

do j = 1,m3
   do i = 1,m2
      do k = 1,m1-2
         !This is Wk* in RAMS tech manual
         wp(k,i,j) = wp(k,i,j) + dts * wt(k,i,j)
      enddo
   enddo
enddo

do j = ja,jz
   do i = ia,iz
      do k = 1,m1-2
         !This is dk(PIk') in RAMS tech manual (dk = dt * theta0 / dz)
         wp(k,i,j) = wp(k,i,j) + a1da2 * acoc(k,i,j) * (pp(k,i,j)-pp(k+1,i,j))
      enddo
   enddo
enddo

return
END SUBROUTINE prdctw1

!##############################################################################
Subroutine prdctw2 (m1,m2,m3,ia,iz,ja,jz  &
   ,wp,pp,acoc,amof,amog,acoaa,rtgt,heatfx1)

use mem_grid

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,i,j,k
real, dimension(m1,m2,m3) :: wp,pp,acoc,amof,amog,acoaa
real, dimension(m2,m3) :: heatfx1,rtgt(m2,m3)

if (nstbot .eq. 1) then
   do j = ja,jz
      do i = ia,iz
         wp(1,i,j) = -heatfx1(i,j) * rtgt(i,j)
      enddo
   enddo
endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!    new trial fix 1/11/96

if (nsttop .eq. 1) then
   do j = ja,jz
      do i = ia,iz
         wp(nz,i,j) = 0.
      enddo
   enddo
endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

if (impl .eq. 1) then

!         First implicit part of the w prediction

   do j = ja,jz
      do i = ia,iz
         do k = 2,m1-2
            wp(k,i,j) = wp(k,i,j) - (pp(k+1,i,j) - pp(k,i,j)) * acoc(k,i,j)
         enddo
      enddo
   enddo

   do j = ja,jz
      do i = ia,iz
         amog(1,i,j) = -wp(1,i,j) / amof(1,i,j)
      enddo
      do k = 2,m1-2
         do i = ia,iz
            amog(k,i,j) = (-wp(k,i,j) - acoaa(k,i,j) * amog(k-1,i,j))  &
               / amof(k,i,j)
         enddo
      enddo
   enddo

endif

return
END SUBROUTINE prdctw2

!##############################################################################
Subroutine prdctw3 (m1,m2,m3,ia,iz,ja,jz,wp,amog,amoe,impl)

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,impl,i,j,k

real, dimension(m1,m2,m3) :: wp,amog,amoe

!     Conclusion of implicit w prediction

if (impl .eq. 1) then
   do k = m1-2,2,-1
      do j = ja,jz
         do i = ia,iz
            wp(k,i,j) = amog(k,i,j) - amoe(k,i,j) * wp(k+1,i,j)
         enddo
      enddo
   enddo
endif

return
END SUBROUTINE prdctw3

!##############################################################################
Subroutine prdctp1 (m1,m2,m3,ia,iz,ja,jz  &
   ,pp,up,vp,pi0,dn0,th0,pt,heatdv,heatfx,f13t,f23t,rtgt,rtgu,rtgv  &
   ,heatfx1,fmapui,fmapvi,dxt,dyt,fmapt)

use mem_grid
use rconstants

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,i,j,k
real :: rocvpct
real, dimension(m1,m2,m3) :: pp,up,vp,pi0,heatdv,pt,heatfx,dn0,th0
real, dimension(m2,m3) :: f13t,f23t,rtgt,heatfx1,rtgu,rtgv,fmapui  &
   ,fmapvi,dxt,dyt,fmapt

CALL azero (m1*m2*m3,heatdv)
rocvpct =rocv *sspct ** 2

!     Divergence calculations for topographical transformation

!     First calculate vertically transformed heat flux

do j = ja,jz
   do i = ia,iz
      do k = 1,m1
         heatfx(k,i,j) = ((up(k,i,j) + up(k,i-1,j)) * f13t(i,j)  &
            + (vp(k,i,j) + vp(k,i,j-jdim)) * f23t(i,j)  &
            ) * dn0(k,i,j) * th0(k,i,j)
      enddo
   enddo
enddo
do j = ja,jz
   do i = ia,iz
      do k = 1,m1-1
         heatfx(k,i,j) = (heatfx(k,i,j) + heatfx(k+1,i,j)) * hw4(k)
      enddo
      heatfx1(i,j) = heatfx(1,i,j) / (.5 * (dn0(1,i,j) * th0(1,i,j)  &
         + dn0(2,i,j) * th0(2,i,j)))
   enddo
enddo

do j = ja,jz
   do i = ia,iz
      do k = 2,m1-1
         heatdv(k,i,j) = (heatfx(k,i,j) - heatfx(k-1,i,j)) * dzt(k)
      enddo
   enddo
enddo
do j = 1,m3
   do i = 1,m2
      do k = 1,m1
         heatfx(k,i,j) = dn0(k,i,j) * th0(k,i,j)
      enddo
   enddo
enddo
do j = ja,jz
   do i = ia,iz
      do k = 2,m1-1

         heatdv(k,i,j) = -rocvpct * pi0(k,i,j) / heatfx(k,i,j)  &
            * (heatdv(k,i,j) + fmapt(i,j)  &
            * ((up(k,i,j) * rtgu(i,j) * fmapui(i,j)  &
            * (heatfx(k,i,j) + heatfx(k,i+1,j))  &
            - up(k,i-1,j) * rtgu(i-1,j) * fmapui(i-1,j)  &
            * (heatfx(k,i,j) + heatfx(k,i-1,j))) * dxt(i,j) * .5  &
            + (vp(k,i,j) * rtgv(i,j) * fmapvi(i,j)  &
            * (heatfx(k,i,j) + heatfx(k,i,j+jdim))  &
            - vp(k,i,j-jdim) * rtgv(i,j-jdim)  &
            * fmapvi(i,j-jdim)  &
            * (heatfx(k,i,j) + heatfx(k,i,j-jdim)))  &
            * dyt(i,j) * .5) / rtgt(i,j))

      enddo
   enddo
enddo


do j = ja,jz
   do i = ia,iz
      do k = 1,m1
         pp(k,i,j) = pp(k,i,j) + (pt(k,i,j) + heatdv(k,i,j)) * dts
      enddo
   enddo
enddo

return
END SUBROUTINE prdctp1

!##############################################################################
Subroutine prdctp2 (m1,m2,m3,ia,iz,ja,jz,pp,wp,acof,acog)

use mem_grid

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,i,j,k
real, dimension(m1,m2,m3) :: pp,wp,acof,acog

!           Finish pressure prediction

do j = ja,jz
   do i = ia,iz
      do k = 2,m1-1
         pp(k,i,j) = pp(k,i,j)  &
            + (wp(k,i,j) * acof(k,i,j) + wp(k-1,i,j) * acog(k,i,j))
      enddo
   enddo
enddo

if (nstbot .eq. 1) CALL botset (m1,m2,m3,pp,'P')

return
END SUBROUTINE prdctp2

!##############################################################################
Subroutine coefz (m1,m2,m3,ia,iz,ja,jz  &
   ,acoc,acof,acog,dn0,pi0,th0,rtgt,a1da2,amoe,amof,acoaa,acobb,acocc)

use mem_grid
use mem_scratch
use rconstants

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,i,j,k

real :: dt2al2,a1da2,rdto2cv,dt2al2r,rdtr
real, dimension(m1,m2,m3) :: th0,pi0,dn0,acoc,acof,acog,amoe,amof,acoaa
real, dimension(m2,m3) :: rtgt
real, dimension(*) :: acobb,acocc

! +--------------------------------------------------------------------+
! \   Calculate coefficients for the vertical pressure gradient        \
! \     and divergence terms.  These will be combined later for the    \
! \     implicit computations.                                         \
! +--------------------------------------------------------------------+

if (impl .eq. 1) then
   dt2al2 = dts * .75
   a1da2 = 1. / 3.
else
   dt2al2 = dts
   a1da2 = 1.
endif
rdto2cv = sspct ** 2 * rgas * dts / (2.0 * cv)

do j = ja,jz
   do i = ia,iz

!         Coefficient for the vertical pressure gradient term

      dt2al2r = .5 * dt2al2 / rtgt(i,j)
      do k = 1,m1-1
         acoc(k,i,j) = dt2al2r * dzm(k) * (th0(k,i,j) + th0(k+1,i,j))
      enddo

!         Coefficients for the vertical divergence term

      rdtr = rdto2cv / rtgt(i,j)
      do k = 2,m1
         vctr12(k) = dn0(k,i,j) * th0(k,i,j)
         vctr11(k) = rdtr * pi0(k,i,j) * dzt(k) / vctr12(k)
      enddo
      vctr12(1) = dn0(1,i,j) * th0(1,i,j)
      do k = 2,m1-1
         acof(k,i,j) = -vctr11(k) * (vctr12(k) + vctr12(k+1))
         acog(k,i,j) = vctr11(k) * (vctr12(k) + vctr12(k-1))
      enddo
      acog(m1,i,j) = vctr11(nzp) * (vctr12(nzp) + vctr12(nz))

      do k = 2,m1-1
         acoaa(k,i,j) = acoc(k,i,j) * acog(k,i,j)
         acobb(k) = acoc(k,i,j) * (acof(k,i,j) - acog(k+1,i,j)) - 1.
         acocc(k) = -acoc(k,i,j) * acof(k+1,i,j)
      enddo
      acobb(1) = -1.
      acocc(1) = 0.
      acoaa(m1,i,j) = 0.
      acobb(m1) = -1.

      amof(1,i,j) = acobb(1)
      amoe(1,i,j) = acocc(1) / amof(1,i,j)
      do k = 2,m1
         amof(k,i,j) = acobb(k) - acoaa(k,i,j) * amoe(k-1,i,j)
         amoe(k,i,j) = acocc(k) / amof(k,i,j)
      enddo

   enddo
enddo

return
END SUBROUTINE coefz

!##############################################################################
Subroutine buoyancy ()

use mem_tend
use mem_basic
use mem_scratch
use mem_grid
use node_mod
use micphys

implicit none

CALL boyanc (mzp,mxp,myp,ia,iz,ja,jz,tend%wt  &
   ,basic_g(ngrid)%theta,basic_g(ngrid)%rtp   &
   ,basic_g(ngrid)%rv,basic_g(ngrid)%th0,dtlt)

return
END SUBROUTINE buoyancy

!##############################################################################
Subroutine boyanc (m1,m2,m3,ia,iz,ja,jz,wt,theta,rtc,rv,th0,dtlt)

use rconstants
use micphys
use mem_basic
use mem_grid, only:ngrid

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,i,j,k
real, dimension(m1,m2,m3) :: wt,theta,rtc,rv,th0
real, dimension(:,:,:), allocatable :: vtemp,wpbuoytheta,wpbuoycond,wpadvdif
real :: dtlt

allocate(vtemp(m1,m2,m3))

if(imbudget>=1) then
 allocate(wpbuoytheta(m1,m2,m3))
 allocate(wpbuoycond(m1,m2,m3))
 allocate(wpadvdif(m1,m2,m3))
 CALL azero (m1*m2*m3,wpbuoytheta)
 CALL azero (m1*m2*m3,wpbuoycond)
 CALL azero (m1*m2*m3,wpadvdif)
endif

if (level .ge. 1) then
   do j = ja,jz
      do i = ia,iz
         do k = 2,m1-1
            vtemp(k,i,j) = gg * ((theta(k,i,j) * (1. + .61 * rv(k,i,j))  &
               - th0(k,i,j)) / th0(k,i,j) - (rtc(k,i,j) - rv(k,i,j)) )
            !calculate partial W buoyancy budgets
            if(imbudget>=1) then
              wpbuoytheta(k,i,j) = gg * ((theta(k,i,j)*(1.+.61*rv(k,i,j)) &
                - th0(k,i,j)) / th0(k,i,j))
              wpbuoycond(k,i,j)  = gg * (-1.0*(rtc(k,i,j) - rv(k,i,j)))
            endif
         enddo
      enddo
   enddo
else
   do j = ja,jz
      do i = ia,iz
         do k = 2,m1-1
            vtemp(k,i,j) = gg * (theta(k,i,j) / th0(k,i,j) - 1.)
            !calculate partial W buoyancy budgets
            if(imbudget>=1) wpbuoytheta(k,i,j) = vtemp(k,i,j)
         enddo
      enddo
   enddo
endif

do j = ja,jz
   do i = ia,iz
      do k = 2,m1-2
         !Calculate W buoyancy budgets (m/s)
         !Calculate W due to advection and diffustion (current wt)
         !Multiply by 2*dt for leapfrog timestep t-dt to t+dt
         if(imbudget>=1) then         
           wpadvdif(k,i,j)    = 2.0 * dtlt * wt(k,i,j)
           wpbuoytheta(k,i,j) = 2.0 * dtlt * (wpbuoytheta(k,i,j) &
                                            + wpbuoytheta(k+1,i,j))
           wpbuoycond(k,i,j)  = 2.0 * dtlt * (wpbuoycond(k,i,j) &
                                            + wpbuoycond(k+1,i,j))
         endif
         wt(k,i,j) = wt(k,i,j) + vtemp(k,i,j) + vtemp(k+1,i,j)
      enddo
   enddo
enddo

deallocate(vtemp)

!Copy local W buoyancy budgets (m/s) to global budget variables
if(imbudget>=1) then
 basic_g(ngrid)%wp_buoy_theta = wpbuoytheta
 basic_g(ngrid)%wp_buoy_cond  = wpbuoycond
 basic_g(ngrid)%wp_advdif     = wpadvdif
 deallocate(wpbuoytheta)
 deallocate(wpbuoycond)
 deallocate(wpadvdif)
endif

return
END SUBROUTINE boyanc
