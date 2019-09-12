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

Subroutine non_scalar_bc (mzp,mxp,myp,mzg,mzs)

use mem_grid
use mem_basic
use mem_radiate
use mem_leaf
use mem_sib
use mem_turb
use mem_micro
use micphys, only:level

implicit none

integer :: mzp,mxp,myp,mzg,mzs

if (ilwrtyp+iswrtyp > 0) then
   if (mod(time + .001,radfrq) .lt. dtlt .or. time .lt. .001) then
     CALL rad_bcond (mzp,mxp,myp,radiate_g(ngrid))
   endif
endif

if (ilwrtyp+iswrtyp > 0) CALL sfcrad_bcond (mxp,myp,radiate_g(ngrid))

CALL leaf_bcond (mxp,myp,mzg,mzs,npatch,leaf_g(ngrid))

CALL turb_bcond (mxp,myp,turb_g(ngrid))

CALL therm_bcond (mzp,mxp,myp,basic_g(ngrid))

if (level==3) CALL micro_bcond (mzp,mxp,myp,micro_g(ngrid))

if (isfcl==2) CALL sib_bcond (mxp,myp,npatch,sib_g(ngrid))

return
END SUBROUTINE non_scalar_bc

!##############################################################################
Subroutine micro_bcond (m1,m2,m3,micro)

use mem_grid
use mem_micro
use micphys

implicit none

type (micro_vars) :: micro

integer :: m1,m2,m3,i,j,k

!Lateral boundary set
 do j = 1,m3
  do k = 1,m1
   if(idriz  >= 1) micro%pcpvd(k,1,j) = 0.
   if(irain  >= 1) micro%pcpvr(k,1,j) = 0.
   if(ipris  >= 1) micro%pcpvp(k,1,j) = 0.
   if(isnow  >= 1) micro%pcpvs(k,1,j) = 0.
   if(iaggr  >= 1) micro%pcpva(k,1,j) = 0.
   if(igraup >= 1) micro%pcpvg(k,1,j) = 0.
   if(ihail  >= 1) micro%pcpvh(k,1,j) = 0.
  enddo
  if(idriz  >= 1) micro%pcprd(1,j) = 0.
  if(irain  >= 1) micro%pcprr(1,j) = 0.
  if(ipris  >= 1) micro%pcprp(1,j) = 0.
  if(isnow  >= 1) micro%pcprs(1,j) = 0.
  if(iaggr  >= 1) micro%pcpra(1,j) = 0.
  if(igraup >= 1) micro%pcprg(1,j) = 0.
  if(ihail  >= 1) micro%pcprh(1,j) = 0.
 enddo

 do j = 1,m3
  do k = 1,m1
   if(idriz  >= 1) micro%pcpvd(k,m2,j) = 0.
   if(irain  >= 1) micro%pcpvr(k,m2,j) = 0.
   if(ipris  >= 1) micro%pcpvp(k,m2,j) = 0.
   if(isnow  >= 1) micro%pcpvs(k,m2,j) = 0.
   if(iaggr  >= 1) micro%pcpva(k,m2,j) = 0.
   if(igraup >= 1) micro%pcpvg(k,m2,j) = 0.
   if(ihail  >= 1) micro%pcpvh(k,m2,j) = 0.
  enddo
   if(idriz  >= 1) micro%pcprd(m2,j) = 0.
   if(irain  >= 1) micro%pcprr(m2,j) = 0.
   if(ipris  >= 1) micro%pcprp(m2,j) = 0.
   if(isnow  >= 1) micro%pcprs(m2,j) = 0.
   if(iaggr  >= 1) micro%pcpra(m2,j) = 0.
   if(igraup >= 1) micro%pcprg(m2,j) = 0.
   if(ihail  >= 1) micro%pcprh(m2,j) = 0.
 enddo

if(jdim .eq. 1) then
  do i = 1,m2
    do k = 1,m1
     if(idriz  >= 1) micro%pcpvd(k,i,1) = 0.
     if(irain  >= 1) micro%pcpvr(k,i,1) = 0.
     if(ipris  >= 1) micro%pcpvp(k,i,1) = 0.
     if(isnow  >= 1) micro%pcpvs(k,i,1) = 0.
     if(iaggr  >= 1) micro%pcpva(k,i,1) = 0.
     if(igraup >= 1) micro%pcpvg(k,i,1) = 0.
     if(ihail  >= 1) micro%pcpvh(k,i,1) = 0.
    enddo
    if(idriz  >= 1) micro%pcprd(i,1) = 0.
    if(irain  >= 1) micro%pcprr(i,1) = 0.
    if(ipris  >= 1) micro%pcprp(i,1) = 0.
    if(isnow  >= 1) micro%pcprs(i,1) = 0.
    if(iaggr  >= 1) micro%pcpra(i,1) = 0.
    if(igraup >= 1) micro%pcprg(i,1) = 0.
    if(ihail  >= 1) micro%pcprh(i,1) = 0.
  enddo
endif

if(jdim .eq. 1) then
  do i = 1,m2
    do k = 1,m1
     if(idriz  >= 1) micro%pcpvd(k,i,m3) = 0.
     if(irain  >= 1) micro%pcpvr(k,i,m3) = 0.
     if(ipris  >= 1) micro%pcpvp(k,i,m3) = 0.
     if(isnow  >= 1) micro%pcpvs(k,i,m3) = 0.
     if(iaggr  >= 1) micro%pcpva(k,i,m3) = 0.
     if(igraup >= 1) micro%pcpvg(k,i,m3) = 0.
     if(ihail  >= 1) micro%pcpvh(k,i,m3) = 0.
    enddo
    if(idriz  >= 1) micro%pcprd(i,m3) = 0.
    if(irain  >= 1) micro%pcprr(i,m3) = 0.
    if(ipris  >= 1) micro%pcprp(i,m3) = 0.
    if(isnow  >= 1) micro%pcprs(i,m3) = 0.
    if(iaggr  >= 1) micro%pcpra(i,m3) = 0.
    if(igraup >= 1) micro%pcprg(i,m3) = 0.
    if(ihail  >= 1) micro%pcprh(i,m3) = 0.
  enddo
endif

return
END SUBROUTINE micro_bcond

!##############################################################################
Subroutine therm_bcond (m1,m2,m3,basic)

use mem_grid
use mem_basic

implicit none

type (basic_vars) :: basic

integer :: m1,m2,m3,i,j,k

!Lateral boundary set
do j = 1,m3
 do k = 1,m1
  basic%theta(k,1,j)  = basic%theta(k,2,j)
  basic%theta(k,m2,j) = basic%theta(k,m2-1,j)
  basic%rv(k,1,j)  = basic%rv(k,2,j)
  basic%rv(k,m2,j) = basic%rv(k,m2-1,j)
 enddo
enddo

if(jdim == 1)then
  do i = 1,m2
   do k = 1,m1
     basic%theta(k,i,1)   = basic%theta(k,i,2)
     basic%theta(k,i,m3)  = basic%theta(k,i,m3-1)
     basic%rv(k,i,1)   = basic%rv(k,i,2)
     basic%rv(k,i,m3)  = basic%rv(k,i,m3-1)
   enddo
  enddo
endif

return
END SUBROUTINE therm_bcond

!##############################################################################
Subroutine leaf_bcond (m2,m3,mzg,mzs,npat,leaf)

use mem_grid
use mem_leaf

implicit none

type (leaf_vars) :: leaf

integer :: m2,m3,mzg,mzs,npat,ipat,i,j,k

!Set BCs for Leaf variables. Do not include veg_ndvif, leaf_class,
! patch_area, and soil_text since these are read from the surface
! files and are boundary set when files are read in by the model.

do ipat = 1,npat

   do j = 1,m3

      leaf%ustar          (1,j,ipat) = leaf%ustar            (2,j,ipat)
      leaf%tstar          (1,j,ipat) = leaf%tstar            (2,j,ipat)
      leaf%rstar          (1,j,ipat) = leaf%rstar            (2,j,ipat)
      leaf%veg_albedo     (1,j,ipat) = leaf%veg_albedo       (2,j,ipat)
      leaf%veg_fracarea   (1,j,ipat) = leaf%veg_fracarea     (2,j,ipat)
      leaf%veg_lai        (1,j,ipat) = leaf%veg_lai          (2,j,ipat)
      leaf%veg_tai        (1,j,ipat) = leaf%veg_tai          (2,j,ipat)
      leaf%veg_rough      (1,j,ipat) = leaf%veg_rough        (2,j,ipat)
      leaf%veg_height     (1,j,ipat) = leaf%veg_height       (2,j,ipat)
      leaf%patch_rough    (1,j,ipat) = leaf%patch_rough      (2,j,ipat)
      leaf%soil_rough     (1,j,ipat) = leaf%soil_rough       (2,j,ipat)
      leaf%sfcwater_nlev  (1,j,ipat) = leaf%sfcwater_nlev    (2,j,ipat)
      leaf%stom_resist    (1,j,ipat) = leaf%stom_resist      (2,j,ipat)
      leaf%ground_rsat    (1,j,ipat) = leaf%ground_rsat      (2,j,ipat)
      leaf%ground_rvap    (1,j,ipat) = leaf%ground_rvap      (2,j,ipat)
      leaf%veg_water      (1,j,ipat) = leaf%veg_water        (2,j,ipat)
      leaf%veg_temp       (1,j,ipat) = leaf%veg_temp         (2,j,ipat)
      leaf%can_rvap       (1,j,ipat) = leaf%can_rvap         (2,j,ipat)
      leaf%can_temp       (1,j,ipat) = leaf%can_temp         (2,j,ipat)
      leaf%veg_ndvip      (1,j,ipat) = leaf%veg_ndvip        (2,j,ipat)
      leaf%veg_ndvic      (1,j,ipat) = leaf%veg_ndvic        (2,j,ipat)
   
      leaf%ustar         (m2,j,ipat) = leaf%ustar         (m2-1,j,ipat)
      leaf%tstar         (m2,j,ipat) = leaf%tstar         (m2-1,j,ipat)
      leaf%rstar         (m2,j,ipat) = leaf%rstar         (m2-1,j,ipat)
      leaf%veg_albedo    (m2,j,ipat) = leaf%veg_albedo    (m2-1,j,ipat)
      leaf%veg_fracarea  (m2,j,ipat) = leaf%veg_fracarea  (m2-1,j,ipat)
      leaf%veg_lai       (m2,j,ipat) = leaf%veg_lai       (m2-1,j,ipat)
      leaf%veg_tai       (m2,j,ipat) = leaf%veg_tai       (m2-1,j,ipat)
      leaf%veg_rough     (m2,j,ipat) = leaf%veg_rough     (m2-1,j,ipat)
      leaf%veg_height    (m2,j,ipat) = leaf%veg_height    (m2-1,j,ipat)
      leaf%patch_rough   (m2,j,ipat) = leaf%patch_rough   (m2-1,j,ipat)
      leaf%soil_rough    (m2,j,ipat) = leaf%soil_rough    (m2-1,j,ipat)
      leaf%sfcwater_nlev (m2,j,ipat) = leaf%sfcwater_nlev (m2-1,j,ipat)
      leaf%stom_resist   (m2,j,ipat) = leaf%stom_resist   (m2-1,j,ipat)
      leaf%ground_rsat   (m2,j,ipat) = leaf%ground_rsat   (m2-1,j,ipat)
      leaf%ground_rvap   (m2,j,ipat) = leaf%ground_rvap   (m2-1,j,ipat)
      leaf%veg_water     (m2,j,ipat) = leaf%veg_water     (m2-1,j,ipat)
      leaf%veg_temp      (m2,j,ipat) = leaf%veg_temp      (m2-1,j,ipat)
      leaf%can_rvap      (m2,j,ipat) = leaf%can_rvap      (m2-1,j,ipat)
      leaf%can_temp      (m2,j,ipat) = leaf%can_temp      (m2-1,j,ipat)
      leaf%veg_ndvip     (m2,j,ipat) = leaf%veg_ndvip     (m2-1,j,ipat)
      leaf%veg_ndvic     (m2,j,ipat) = leaf%veg_ndvic     (m2-1,j,ipat)
   
      do k = 1,mzg
         leaf%soil_water       (k,1,j,ipat) = leaf%soil_water         (k,2,j,ipat)
         leaf%soil_energy      (k,1,j,ipat) = leaf%soil_energy        (k,2,j,ipat)
         leaf%soil_water      (k,m2,j,ipat) = leaf%soil_water      (k,m2-1,j,ipat)
         leaf%soil_energy     (k,m2,j,ipat) = leaf%soil_energy     (k,m2-1,j,ipat)
      enddo

      do k = 1,mzs
         leaf%sfcwater_mass    (k,1,j,ipat) = leaf%sfcwater_mass      (k,2,j,ipat)
         leaf%sfcwater_energy  (k,1,j,ipat) = leaf%sfcwater_energy    (k,2,j,ipat)
         leaf%sfcwater_depth   (k,1,j,ipat) = leaf%sfcwater_depth     (k,2,j,ipat)
         leaf%sfcwater_mass   (k,m2,j,ipat) = leaf%sfcwater_mass   (k,m2-1,j,ipat)
         leaf%sfcwater_energy (k,m2,j,ipat) = leaf%sfcwater_energy (k,m2-1,j,ipat)
         leaf%sfcwater_depth  (k,m2,j,ipat) = leaf%sfcwater_depth  (k,m2-1,j,ipat)
      enddo

   enddo   

   if (jdim == 1) then

      do i = 1,m2
         leaf%ustar          (i,1,ipat) = leaf%ustar            (i,2,ipat)
         leaf%tstar          (i,1,ipat) = leaf%tstar            (i,2,ipat)
         leaf%rstar          (i,1,ipat) = leaf%rstar            (i,2,ipat)
         leaf%veg_albedo     (i,1,ipat) = leaf%veg_albedo       (i,2,ipat)
         leaf%veg_fracarea   (i,1,ipat) = leaf%veg_fracarea     (i,2,ipat)
         leaf%veg_lai        (i,1,ipat) = leaf%veg_lai          (i,2,ipat)
         leaf%veg_tai        (i,1,ipat) = leaf%veg_tai          (i,2,ipat)
         leaf%veg_rough      (i,1,ipat) = leaf%veg_rough        (i,2,ipat)
         leaf%veg_height     (i,1,ipat) = leaf%veg_height       (i,2,ipat)
         leaf%patch_rough    (i,1,ipat) = leaf%patch_rough      (i,2,ipat)
         leaf%soil_rough     (i,1,ipat) = leaf%soil_rough       (i,2,ipat)
         leaf%sfcwater_nlev  (i,1,ipat) = leaf%sfcwater_nlev    (i,2,ipat)
         leaf%stom_resist    (i,1,ipat) = leaf%stom_resist      (i,2,ipat)
         leaf%ground_rsat    (i,1,ipat) = leaf%ground_rsat      (i,2,ipat)
         leaf%ground_rvap    (i,1,ipat) = leaf%ground_rvap      (i,2,ipat)
         leaf%veg_water      (i,1,ipat) = leaf%veg_water        (i,2,ipat)
         leaf%veg_temp       (i,1,ipat) = leaf%veg_temp         (i,2,ipat)
         leaf%can_rvap       (i,1,ipat) = leaf%can_rvap         (i,2,ipat)
         leaf%can_temp       (i,1,ipat) = leaf%can_temp         (i,2,ipat)
         leaf%veg_ndvip      (i,1,ipat) = leaf%veg_ndvip        (i,2,ipat)
         leaf%veg_ndvic      (i,1,ipat) = leaf%veg_ndvic        (i,2,ipat)
   
         leaf%ustar         (i,m3,ipat) = leaf%ustar         (i,m3-1,ipat)
         leaf%tstar         (i,m3,ipat) = leaf%tstar         (i,m3-1,ipat)
         leaf%rstar         (i,m3,ipat) = leaf%rstar         (i,m3-1,ipat)
         leaf%veg_albedo    (i,m3,ipat) = leaf%veg_albedo    (i,m3-1,ipat)
         leaf%veg_fracarea  (i,m3,ipat) = leaf%veg_fracarea  (i,m3-1,ipat)
         leaf%veg_lai       (i,m3,ipat) = leaf%veg_lai       (i,m3-1,ipat)
         leaf%veg_tai       (i,m3,ipat) = leaf%veg_tai       (i,m3-1,ipat)
         leaf%veg_rough     (i,m3,ipat) = leaf%veg_rough     (i,m3-1,ipat)
         leaf%veg_height    (i,m3,ipat) = leaf%veg_height    (i,m3-1,ipat)
         leaf%patch_rough   (i,m3,ipat) = leaf%patch_rough   (i,m3-1,ipat)
         leaf%soil_rough    (i,m3,ipat) = leaf%soil_rough    (i,m3-1,ipat)
         leaf%sfcwater_nlev (i,m3,ipat) = leaf%sfcwater_nlev (i,m3-1,ipat)
         leaf%stom_resist   (i,m3,ipat) = leaf%stom_resist   (i,m3-1,ipat)
         leaf%ground_rsat   (i,m3,ipat) = leaf%ground_rsat   (i,m3-1,ipat)
         leaf%ground_rvap   (i,m3,ipat) = leaf%ground_rvap   (i,m3-1,ipat)
         leaf%veg_water     (i,m3,ipat) = leaf%veg_water     (i,m3-1,ipat)
         leaf%veg_temp      (i,m3,ipat) = leaf%veg_temp      (i,m3-1,ipat)
         leaf%can_rvap      (i,m3,ipat) = leaf%can_rvap      (i,m3-1,ipat)
         leaf%can_temp      (i,m3,ipat) = leaf%can_temp      (i,m3-1,ipat)
         leaf%veg_ndvip     (i,m3,ipat) = leaf%veg_ndvip     (i,m3-1,ipat)
         leaf%veg_ndvic     (i,m3,ipat) = leaf%veg_ndvic     (i,m3-1,ipat)
   
         do k = 1,mzg
            leaf%soil_water       (k,i,1,ipat) = leaf%soil_water         (k,i,2,ipat)
            leaf%soil_energy      (k,i,1,ipat) = leaf%soil_energy        (k,i,2,ipat)
            leaf%soil_water      (k,i,m3,ipat) = leaf%soil_water      (k,i,m3-1,ipat)
            leaf%soil_energy     (k,i,m3,ipat) = leaf%soil_energy     (k,i,m3-1,ipat)
         enddo

         do k = 1,mzs
            leaf%sfcwater_mass    (k,i,1,ipat) = leaf%sfcwater_mass      (k,i,2,ipat)
            leaf%sfcwater_energy  (k,i,1,ipat) = leaf%sfcwater_energy    (k,i,2,ipat)
            leaf%sfcwater_depth   (k,i,1,ipat) = leaf%sfcwater_depth     (k,i,2,ipat)
            leaf%sfcwater_mass   (k,i,m3,ipat) = leaf%sfcwater_mass   (k,i,m3-1,ipat)
            leaf%sfcwater_energy (k,i,m3,ipat) = leaf%sfcwater_energy (k,i,m3-1,ipat)
            leaf%sfcwater_depth  (k,i,m3,ipat) = leaf%sfcwater_depth  (k,i,m3-1,ipat)
         enddo

      enddo   

   endif

enddo

return
END SUBROUTINE leaf_bcond

!##############################################################################
Subroutine rad_bcond (m1,m2,m3,radiate)

use mem_grid
use mem_radiate

implicit none

type (radiate_vars) :: radiate

integer :: m1,m2,m3,i,j,k
real :: dzmr,dztr

dzmr = dzm(m1-2) / dzm(m1-1)
dztr = dzt(m1-2) / dzt(m1-1)

!Top and bottom boundary set
!Swup, swdn, lwup, lwdn top and bottom boundary conditions are done
!in radcalc3 (Harrington radiation) - these are on "m" levels and 
!are based on the fluxes.
do j = 1,m3
 do i = 1,m2
    radiate%fthrd(m1,i,j) = max(0.,radiate%fthrd(m1-1,i,j) &
           +dzmr*(radiate%fthrd(m1-1,i,j)-radiate%fthrd(m1-2,i,j)))
    radiate%fthrd(1,i,j)  = radiate%fthrd(2,i,j)
    if(ilwrtyp == 3 .or. iswrtyp == 3) then
      radiate%bext(m1,i,j) = max(0.,radiate%bext(m1-1,i,j) &
           +dzmr*(radiate%bext(m1-1,i,j)-radiate%bext(m1-2,i,j)))
      radiate%bext(1,i,j)  = radiate%bext(2,i,j)
    endif
 enddo
enddo

!Lateral boundary set
do j = 1,m3
 do k = 1,m1
   radiate%fthrd(k,1,j)  = radiate%fthrd(k,2,j)
   if(ilwrtyp == 3 .or. iswrtyp == 3) radiate%bext(k,1,j)   = radiate%bext(k,2,j)
   if(iswrtyp == 3) then
     radiate%swup(k,1,j) = radiate%swup(k,2,j)
     radiate%swdn(k,1,j) = radiate%swdn(k,2,j)
   endif
   if(ilwrtyp == 3) then
     radiate%lwup(k,1,j) = radiate%lwup(k,2,j)
     radiate%lwdn(k,1,j) = radiate%lwdn(k,2,j)
   endif
   radiate%fthrd(k,m2,j) = radiate%fthrd(k,m2-1,j)
   if(ilwrtyp == 3 .or. iswrtyp == 3) radiate%bext(k,m2,j)  = radiate%bext(k,m2-1,j)
   if(iswrtyp == 3) then
     radiate%swup(k,m2,j) = radiate%swup(k,m2-1,j)
     radiate%swdn(k,m2,j) = radiate%swdn(k,m2-1,j)
   endif
   if(ilwrtyp == 3) then
     radiate%lwup(k,m2,j) = radiate%lwup(k,m2-1,j)
     radiate%lwdn(k,m2,j) = radiate%lwdn(k,m2-1,j)
   endif
 enddo
enddo

if(jdim == 1)then
  do i = 1,m2
   do k = 1,m1
     radiate%fthrd(k,i,1)   = radiate%fthrd(k,i,2)
     if(ilwrtyp == 3 .or. iswrtyp == 3) radiate%bext(k,i,1)    = radiate%bext(k,i,2)
     if(iswrtyp == 3) then
       radiate%swup(k,i,1) = radiate%swup(k,i,2)
       radiate%swdn(k,i,1) = radiate%swdn(k,i,2)
     endif
     if(ilwrtyp == 3) then
       radiate%lwup(k,i,1) = radiate%lwup(k,i,2)
       radiate%lwdn(k,i,1) = radiate%lwdn(k,i,2)
     endif
     radiate%fthrd(k,i,m3)  = radiate%fthrd(k,i,m3-1)
     if(ilwrtyp == 3 .or. iswrtyp == 3) radiate%bext(k,i,m3)   = radiate%bext(k,i,m3-1)
     if(iswrtyp == 3) then
       radiate%swup(k,i,m3) = radiate%swup(k,i,m3-1)
       radiate%swdn(k,i,m3) = radiate%swdn(k,i,m3-1)
     endif
     if(ilwrtyp == 3) then
       radiate%lwup(k,i,m3) = radiate%lwup(k,i,m3-1)
       radiate%lwdn(k,i,m3) = radiate%lwdn(k,i,m3-1)
     endif
   enddo
  enddo
endif

return
END SUBROUTINE rad_bcond

!##############################################################################
Subroutine turb_bcond (m2,m3,turb)

use mem_grid
use mem_turb

implicit none

type (turb_vars) :: turb

integer :: m2,m3,i,j

do j = 1,m3
  turb%sflux_u(1,j) = turb%sflux_u(2,j)
  turb%sflux_v(1,j) = turb%sflux_v(2,j)
  turb%sflux_w(1,j) = turb%sflux_w(2,j)
  turb%sflux_t(1,j) = turb%sflux_t(2,j)
  turb%sflux_r(1,j) = turb%sflux_r(2,j)
  turb%sflux_u(m2,j) = turb%sflux_u(m2-1,j)
  turb%sflux_v(m2,j) = turb%sflux_v(m2-1,j)
  turb%sflux_w(m2,j) = turb%sflux_w(m2-1,j)
  turb%sflux_t(m2,j) = turb%sflux_t(m2-1,j)
  turb%sflux_r(m2,j) = turb%sflux_r(m2-1,j)
enddo

if(jdim == 1)then
  do i = 1,m2
     turb%sflux_u(i,1) = turb%sflux_u(i,2)
     turb%sflux_v(i,1) = turb%sflux_v(i,2)
     turb%sflux_w(i,1) = turb%sflux_w(i,2)
     turb%sflux_t(i,1) = turb%sflux_t(i,2)
     turb%sflux_r(i,1) = turb%sflux_r(i,2)
     turb%sflux_u(i,m3) = turb%sflux_u(i,m3-1)
     turb%sflux_v(i,m3) = turb%sflux_v(i,m3-1)
     turb%sflux_w(i,m3) = turb%sflux_w(i,m3-1)
     turb%sflux_t(i,m3) = turb%sflux_t(i,m3-1)
     turb%sflux_r(i,m3) = turb%sflux_r(i,m3-1)
  enddo
endif

return
END SUBROUTINE turb_bcond

!##############################################################################
Subroutine sfcrad_bcond (m2,m3,radiate)

use mem_grid
use mem_radiate

implicit none

type (radiate_vars) :: radiate

integer :: m2,m3,i,j

do j = 1,m3
  radiate%rshort(1,j)  = radiate%rshort(2,j)
  radiate%rlong(1,j)   = radiate%rlong(2,j)
  radiate%rlongup(1,j) = radiate%rlongup(2,j)
  radiate%albedt(1,j)  = radiate%albedt(2,j)
  radiate%cosz(1,j)    = radiate%cosz(2,j)

  radiate%rshort(m2,j)  = radiate%rshort(m2-1,j)
  radiate%rlong(m2,j)   = radiate%rlong(m2-1,j)
  radiate%rlongup(m2,j) = radiate%rlongup(m2-1,j)
  radiate%albedt(m2,j)  = radiate%albedt(m2-1,j)
  radiate%cosz(m2,j)    = radiate%cosz(m2-1,j)
enddo

if(jdim == 1)then
  do i = 1,m2
     radiate%rshort(i,1)  = radiate%rshort(i,2)
     radiate%rlong(i,1)   = radiate%rlong(i,2)
     radiate%rlongup(i,1) = radiate%rlongup(i,2)
     radiate%albedt(i,1)  = radiate%albedt(i,2)
     radiate%cosz(i,1)    = radiate%cosz(i,2)

     radiate%rshort(i,m3)  = radiate%rshort(i,m3-1)
     radiate%rlong(i,m3)   = radiate%rlong(i,m3-1)
     radiate%rlongup(i,m3) = radiate%rlongup(i,m3-1)
     radiate%albedt(i,m3)  = radiate%albedt(i,m3-1)
     radiate%cosz(i,m3)    = radiate%cosz(i,m3-1)
   enddo
endif

return
END SUBROUTINE sfcrad_bcond

!##############################################################################
Subroutine sib_bcond (m2,m3,npat,sib)

use mem_grid
use mem_sib

implicit none

type (sib_vars) :: sib

integer :: m2,m3,npat,ipat,i,j

!Set BCs for SiB variables.

do ipat = 1,npat

   do j = 1,m3

      sib%snow1(1,j,ipat)   = sib%snow1(2,j,ipat)
      sib%snow2(1,j,ipat)   = sib%snow2(2,j,ipat)
      sib%capac1(1,j,ipat)  = sib%capac1(2,j,ipat)
      sib%capac2(1,j,ipat)  = sib%capac2(2,j,ipat)
      sib%pco2ap(1,j,ipat)  = sib%pco2ap(2,j,ipat)
      sib%co2flx(1,j,ipat)  = sib%co2flx(2,j,ipat)
      sib%sfcswa(1,j,ipat)  = sib%sfcswa(2,j,ipat)
      sib%uplwrf(1,j,ipat)  = sib%uplwrf(2,j,ipat)
      sib%assimn(1,j,ipat)  = sib%assimn(2,j,ipat)
      sib%respg(1,j,ipat)   = sib%respg(2,j,ipat)
      sib%rstfac1(1,j,ipat) = sib%rstfac1(2,j,ipat)
      sib%rstfac2(1,j,ipat) = sib%rstfac2(2,j,ipat)
      sib%rstfac3(1,j,ipat) = sib%rstfac3(2,j,ipat)
      sib%ect(1,j,ipat)     = sib%ect(2,j,ipat)
      sib%eci(1,j,ipat)     = sib%eci(2,j,ipat)
      sib%egi(1,j,ipat)     = sib%egi(2,j,ipat)
      sib%egs(1,j,ipat)     = sib%egs(2,j,ipat)
      sib%hc(1,j,ipat)      = sib%hc(2,j,ipat)
      sib%hg(1,j,ipat)      = sib%hg(2,j,ipat)
      sib%ra(1,j,ipat)      = sib%ra(2,j,ipat)
      sib%rb(1,j,ipat)      = sib%rb(2,j,ipat)
      sib%rc(1,j,ipat)      = sib%rc(2,j,ipat)
      sib%rd(1,j,ipat)      = sib%rd(2,j,ipat)
      sib%roff(1,j,ipat)    = sib%roff(2,j,ipat)
      sib%green(1,j,ipat)   = sib%green(2,j,ipat)
      sib%apar(1,j,ipat)    = sib%apar(2,j,ipat)
      sib%ventmf(1,j,ipat)  = sib%ventmf(2,j,ipat)
      sib%pco2c(1,j,ipat)   = sib%pco2c(2,j,ipat)
      sib%pco2i(1,j,ipat)   = sib%pco2i(2,j,ipat)
      sib%pco2s(1,j,ipat)   = sib%pco2s(2,j,ipat)
      sib%pco2m(1,j,ipat)   = sib%pco2m(2,j,ipat)
      sib%ea(1,j,ipat)      = sib%ea(2,j,ipat)
      sib%em(1,j,ipat)      = sib%em(2,j,ipat)
      sib%rha(1,j,ipat)     = sib%rha(2,j,ipat)
      sib%radvbc(1,j,ipat)  = sib%radvbc(2,j,ipat)
      sib%radvdc(1,j,ipat)  = sib%radvdc(2,j,ipat)
      sib%radnbc(1,j,ipat)  = sib%radnbc(2,j,ipat)
      sib%radndc(1,j,ipat)  = sib%radndc(2,j,ipat)
      sib%psy(1,j,ipat)     = sib%psy(2,j,ipat)

      sib%snow1(m2,j,ipat)   = sib%snow1(m2-1,j,ipat)
      sib%snow2(m2,j,ipat)   = sib%snow2(m2-1,j,ipat)
      sib%capac1(m2,j,ipat)  = sib%capac1(m2-1,j,ipat)
      sib%capac2(m2,j,ipat)  = sib%capac2(m2-1,j,ipat)
      sib%pco2ap(m2,j,ipat)  = sib%pco2ap(m2-1,j,ipat)
      sib%co2flx(m2,j,ipat)  = sib%co2flx(m2-1,j,ipat)
      sib%sfcswa(m2,j,ipat)  = sib%sfcswa(m2-1,j,ipat)
      sib%uplwrf(m2,j,ipat)  = sib%uplwrf(m2-1,j,ipat)
      sib%assimn(m2,j,ipat)  = sib%assimn(m2-1,j,ipat)
      sib%respg(m2,j,ipat)   = sib%respg(m2-1,j,ipat)
      sib%rstfac1(m2,j,ipat) = sib%rstfac1(m2-1,j,ipat)
      sib%rstfac2(m2,j,ipat) = sib%rstfac2(m2-1,j,ipat)
      sib%rstfac3(m2,j,ipat) = sib%rstfac3(m2-1,j,ipat)
      sib%ect(m2,j,ipat)     = sib%ect(m2-1,j,ipat)
      sib%eci(m2,j,ipat)     = sib%eci(m2-1,j,ipat)
      sib%egi(m2,j,ipat)     = sib%egi(m2-1,j,ipat)
      sib%egs(m2,j,ipat)     = sib%egs(m2-1,j,ipat)
      sib%hc(m2,j,ipat)      = sib%hc(m2-1,j,ipat)
      sib%hg(m2,j,ipat)      = sib%hg(m2-1,j,ipat)
      sib%ra(m2,j,ipat)      = sib%ra(m2-1,j,ipat)
      sib%rb(m2,j,ipat)      = sib%rb(m2-1,j,ipat)
      sib%rc(m2,j,ipat)      = sib%rc(m2-1,j,ipat)
      sib%rd(m2,j,ipat)      = sib%rd(m2-1,j,ipat)
      sib%roff(m2,j,ipat)    = sib%roff(m2-1,j,ipat)
      sib%green(m2,j,ipat)   = sib%green(m2-1,j,ipat)
      sib%apar(m2,j,ipat)    = sib%apar(m2-1,j,ipat)
      sib%ventmf(m2,j,ipat)  = sib%ventmf(m2-1,j,ipat)
      sib%pco2c(m2,j,ipat)   = sib%pco2c(m2-1,j,ipat)
      sib%pco2i(m2,j,ipat)   = sib%pco2i(m2-1,j,ipat)
      sib%pco2s(m2,j,ipat)   = sib%pco2s(m2-1,j,ipat)
      sib%pco2m(m2,j,ipat)   = sib%pco2m(m2-1,j,ipat)
      sib%ea(m2,j,ipat)      = sib%ea(m2-1,j,ipat)
      sib%em(m2,j,ipat)      = sib%em(m2-1,j,ipat)
      sib%rha(m2,j,ipat)     = sib%rha(m2-1,j,ipat)
      sib%radvbc(m2,j,ipat)  = sib%radvbc(m2-1,j,ipat)
      sib%radvdc(m2,j,ipat)  = sib%radvdc(m2-1,j,ipat)
      sib%radnbc(m2,j,ipat)  = sib%radnbc(m2-1,j,ipat)
      sib%radndc(m2,j,ipat)  = sib%radndc(m2-1,j,ipat)
      sib%psy(m2,j,ipat)     = sib%psy(m2-1,j,ipat)

   enddo

   if (jdim == 1) then

    do i = 1,m2
      sib%snow1(i,1,ipat)   = sib%snow1(i,2,ipat)
      sib%snow2(i,1,ipat)   = sib%snow2(i,2,ipat)
      sib%capac1(i,1,ipat)  = sib%capac1(i,2,ipat)
      sib%capac2(i,1,ipat)  = sib%capac2(i,2,ipat)
      sib%pco2ap(i,1,ipat)  = sib%pco2ap(i,2,ipat)
      sib%co2flx(i,1,ipat)  = sib%co2flx(i,2,ipat)
      sib%sfcswa(i,1,ipat)  = sib%sfcswa(i,2,ipat)
      sib%uplwrf(i,1,ipat)  = sib%uplwrf(i,2,ipat)
      sib%assimn(i,1,ipat)  = sib%assimn(i,2,ipat)
      sib%respg(i,1,ipat)   = sib%respg(i,2,ipat)
      sib%rstfac1(i,1,ipat) = sib%rstfac1(i,2,ipat)
      sib%rstfac2(i,1,ipat) = sib%rstfac2(i,2,ipat)
      sib%rstfac3(i,1,ipat) = sib%rstfac3(i,2,ipat)
      sib%ect(i,1,ipat)     = sib%ect(i,2,ipat)
      sib%eci(i,1,ipat)     = sib%eci(i,2,ipat)
      sib%egi(i,1,ipat)     = sib%egi(i,2,ipat)
      sib%egs(i,1,ipat)     = sib%egs(i,2,ipat)
      sib%hc(i,1,ipat)      = sib%hc(i,2,ipat)
      sib%hg(i,1,ipat)      = sib%hg(i,2,ipat)
      sib%ra(i,1,ipat)      = sib%ra(i,2,ipat)
      sib%rb(i,1,ipat)      = sib%rb(i,2,ipat)
      sib%rc(i,1,ipat)      = sib%rc(i,2,ipat)
      sib%rd(i,1,ipat)      = sib%rd(i,2,ipat)
      sib%roff(i,1,ipat)    = sib%roff(i,2,ipat)
      sib%green(i,1,ipat)   = sib%green(i,2,ipat)
      sib%apar(i,1,ipat)    = sib%apar(i,2,ipat)
      sib%ventmf(i,1,ipat)  = sib%ventmf(i,2,ipat)
      sib%pco2c(i,1,ipat)   = sib%pco2c(i,2,ipat)
      sib%pco2i(i,1,ipat)   = sib%pco2i(i,2,ipat)
      sib%pco2s(i,1,ipat)   = sib%pco2s(i,2,ipat)
      sib%pco2m(i,1,ipat)   = sib%pco2m(i,2,ipat)
      sib%ea(i,1,ipat)      = sib%ea(i,2,ipat)
      sib%em(i,1,ipat)      = sib%em(i,2,ipat)
      sib%rha(i,1,ipat)     = sib%rha(i,2,ipat)
      sib%radvbc(i,1,ipat)  = sib%radvbc(i,2,ipat)
      sib%radvdc(i,1,ipat)  = sib%radvdc(i,2,ipat)
      sib%radnbc(i,1,ipat)  = sib%radnbc(i,2,ipat)
      sib%radndc(i,1,ipat)  = sib%radndc(i,2,ipat)
      sib%psy(i,1,ipat)     = sib%psy(i,2,ipat)

      sib%snow1(i,m3,ipat)   = sib%snow1(i,m3-1,ipat)
      sib%snow2(i,m3,ipat)   = sib%snow2(i,m3-1,ipat)
      sib%capac1(i,m3,ipat)  = sib%capac1(i,m3-1,ipat)
      sib%capac2(i,m3,ipat)  = sib%capac2(i,m3-1,ipat)
      sib%pco2ap(i,m3,ipat)  = sib%pco2ap(i,m3-1,ipat)
      sib%co2flx(i,m3,ipat)  = sib%co2flx(i,m3-1,ipat)
      sib%sfcswa(i,m3,ipat)  = sib%sfcswa(i,m3-1,ipat)
      sib%uplwrf(i,m3,ipat)  = sib%uplwrf(i,m3-1,ipat)
      sib%assimn(i,m3,ipat)  = sib%assimn(i,m3-1,ipat)
      sib%respg(i,m3,ipat)   = sib%respg(i,m3-1,ipat)
      sib%rstfac1(i,m3,ipat) = sib%rstfac1(i,m3-1,ipat)
      sib%rstfac2(i,m3,ipat) = sib%rstfac2(i,m3-1,ipat)
      sib%rstfac3(i,m3,ipat) = sib%rstfac3(i,m3-1,ipat)
      sib%ect(i,m3,ipat)     = sib%ect(i,m3-1,ipat)
      sib%eci(i,m3,ipat)     = sib%eci(i,m3-1,ipat)
      sib%egi(i,m3,ipat)     = sib%egi(i,m3-1,ipat)
      sib%egs(i,m3,ipat)     = sib%egs(i,m3-1,ipat)
      sib%hc(i,m3,ipat)      = sib%hc(i,m3-1,ipat)
      sib%hg(i,m3,ipat)      = sib%hg(i,m3-1,ipat)
      sib%ra(i,m3,ipat)      = sib%ra(i,m3-1,ipat)
      sib%rb(i,m3,ipat)      = sib%rb(i,m3-1,ipat)
      sib%rc(i,m3,ipat)      = sib%rc(i,m3-1,ipat)
      sib%rd(i,m3,ipat)      = sib%rd(i,m3-1,ipat)
      sib%roff(i,m3,ipat)    = sib%roff(i,m3-1,ipat)
      sib%green(i,m3,ipat)   = sib%green(i,m3-1,ipat)
      sib%apar(i,m3,ipat)    = sib%apar(i,m3-1,ipat)
      sib%ventmf(i,m3,ipat)  = sib%ventmf(i,m3-1,ipat)
      sib%pco2c(i,m3,ipat)   = sib%pco2c(i,m3-1,ipat)
      sib%pco2i(i,m3,ipat)   = sib%pco2i(i,m3-1,ipat)
      sib%pco2s(i,m3,ipat)   = sib%pco2s(i,m3-1,ipat)
      sib%pco2m(i,m3,ipat)   = sib%pco2m(i,m3-1,ipat)
      sib%ea(i,m3,ipat)      = sib%ea(i,m3-1,ipat)
      sib%em(i,m3,ipat)      = sib%em(i,m3-1,ipat)
      sib%rha(i,m3,ipat)     = sib%rha(i,m3-1,ipat)
      sib%radvbc(i,m3,ipat)  = sib%radvbc(i,m3-1,ipat)
      sib%radvdc(i,m3,ipat)  = sib%radvdc(i,m3-1,ipat)
      sib%radnbc(i,m3,ipat)  = sib%radnbc(i,m3-1,ipat)
      sib%radndc(i,m3,ipat)  = sib%radndc(i,m3-1,ipat)
      sib%psy(i,m3,ipat)     = sib%psy(i,m3-1,ipat)

    enddo   

   endif

enddo

return
END SUBROUTINE sib_bcond

