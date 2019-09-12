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

Subroutine radiate (mzp,mxp,myp,ia,iz,ja,jz)

use mem_tend
use mem_grid
use mem_leaf
use mem_sib
use mem_radiate
use mem_basic
use mem_scratch
use mem_micro
use rconstants
use rrad3
use micphys
use ref_sounding

implicit none

integer :: mzp,mxp,myp,ia,iz,ja,jz

real, save :: prsnz,prsnzp

integer, save :: ncall=0

if (ilwrtyp + iswrtyp .eq. 0) return

CALL tend_accum (mzp,mxp,myp,ia,iz,ja,jz,tend%tht(1)  &
   ,radiate_g(ngrid)%fthrd(1,1,1))

if (mod(time + .001,radfrq) .lt. dtlt .or. time .lt. 0.001) then

   if(iprntstmt>=1 .and. print_msg) &
      print 90,time,time/3600.+(itime1/100+mod(itime1,100)/60.)
90      format(' Radiation Tendencies Updated Time =',F10.1,  &
        '  UTC TIME (HRS) =',F6.1)

! Compute solar zenith angle, multiplier for solar constant, sfc albedo,
! and surface upward longwave radiation.

   CALL radprep (mxp,myp,nzg,nzs,npatch,ia,iz,ja,jz,jday   &

      ,leaf_g(ngrid)%soil_water      (1,1,1,1)  &
      ,leaf_g(ngrid)%soil_energy     (1,1,1,1)  &
      ,leaf_g(ngrid)%soil_text       (1,1,1,1)  &
      ,leaf_g(ngrid)%sfcwater_energy (1,1,1,1)  &
      ,leaf_g(ngrid)%sfcwater_depth  (1,1,1,1)  &
      ,leaf_g(ngrid)%leaf_class      (1,1,1)    &
      ,leaf_g(ngrid)%veg_fracarea    (1,1,1)    &
      ,leaf_g(ngrid)%veg_height      (1,1,1)    &
      ,leaf_g(ngrid)%veg_albedo      (1,1,1)    &
      ,leaf_g(ngrid)%patch_area      (1,1,1)    &
      ,leaf_g(ngrid)%sfcwater_nlev   (1,1,1)    &
      ,leaf_g(ngrid)%veg_temp        (1,1,1)    &
      ,leaf_g(ngrid)%can_temp        (1,1,1)    &
      ,solfac                                   &
      ,grid_g(ngrid)%glat            (1,1)      &
      ,grid_g(ngrid)%glon            (1,1)      &
      ,radiate_g(ngrid)%rshort       (1,1)      &
      ,radiate_g(ngrid)%rlong        (1,1)      &
      ,radiate_g(ngrid)%rlongup      (1,1)      &
      ,radiate_g(ngrid)%albedt       (1,1)      &
      ,radiate_g(ngrid)%cosz         (1,1)      &
      ,sib_g(ngrid)%sfcswa           (1,1,1)    &
      ,sib_g(ngrid)%uplwrf           (1,1,1)    )

   CALL azero (mzp*mxp*myp,radiate_g(ngrid)%fthrd(1,1,1))

   if (ilwrtyp .le. 2 .or. iswrtyp .le. 2) then
   
! If using Mahrer-Pielke and/or Chen-Cotton radiation, call radcomp.

      CALL radcomp (mzp,mxp,myp,ia,iz,ja,jz,solfac  &
         ,basic_g(ngrid)%theta     (1,1,1)  &
         ,basic_g(ngrid)%pi0       (1,1,1)  &
         ,basic_g(ngrid)%pp        (1,1,1)  &
         ,basic_g(ngrid)%rv        (1,1,1)  &
         ,basic_g(ngrid)%dn0       (1,1,1)  &
         ,basic_g(ngrid)%rtp       (1,1,1)  &
         ,radiate_g(ngrid)%fthrd   (1,1,1)  &
         ,grid_g(ngrid)%rtgt       (1,1)    &
         ,grid_g(ngrid)%f13t       (1,1)    &
         ,grid_g(ngrid)%f23t       (1,1)    &
         ,grid_g(ngrid)%glon       (1,1)    &
         ,radiate_g(ngrid)%rshort  (1,1)    &
         ,radiate_g(ngrid)%rlong   (1,1)    &
         ,radiate_g(ngrid)%albedt  (1,1)    &
         ,radiate_g(ngrid)%cosz    (1,1)    &
         ,radiate_g(ngrid)%rlongup (1,1))

   endif

   if (iswrtyp .eq. 3 .or. ilwrtyp .eq. 3) then

! Using Harrington radiation

      if (ncall .eq. 0) then

! If first call for this node, initialize several quantities & Mclatchy
! sounding data.

         CALL radinit (ng,nb,nsolb,npsb,nuum,prf,alpha,trf,beta  &
            ,xp,wght,wlenlo,wlenhi,solar0,ralcs,a0,a1,a2,a3  &
            ,exptabc,ulim,npartob,npartg,ncog,ncb  &
            ,ocoef,bcoef,gcoef,gnu)

         prsnz  = (pi01dn(nnzp(1)-1,1) / cp) ** cpor * p00
         prsnzp = (pi01dn(nnzp(1)  ,1) / cp) ** cpor * p00

        CALL mclatchy (1,mzp  &
            ,prsnz,prsnzp  &
            ,grid_g(ngrid)%glat       (1,1)  &
            ,grid_g(ngrid)%rtgt       (1,1)  &
            ,grid_g(ngrid)%topt       (1,1)  &
            ,radiate_g(ngrid)%rlongup (1,1)  &
            ,zm,zt,vctr1,vctr2,vctr3,vctr4,vctr5,vctr6,vctr7  &
            ,vctr8,vctr9,vctr10,vctr11,vctr12)

         ncall = ncall + 1
      endif

! For any call, interpolate the mclatchy sounding data by latitude and
! season.

     CALL mclatchy (2,mzp  &
         ,prsnz,prsnzp  &
         ,grid_g(ngrid)%glat       (1,1)  &
         ,grid_g(ngrid)%rtgt       (1,1)  &
         ,grid_g(ngrid)%topt       (1,1)  &
         ,radiate_g(ngrid)%rlongup (1,1)  &
         ,zm,zt,vctr1,vctr2,vctr3,vctr4,vctr5,vctr6,vctr7  &
         ,vctr8,vctr9,vctr10,vctr11,vctr12)

! If using Harrington radiation with moisture complexity LEVEL < 3,
! call radcomp3 which is a substitute driving structure to the bulk
! microphysics.

      if (level <= 2) then
         ! The only difference between the calls to radcomp3 for level == 2
         ! and level != 2 is the variable given for argument #22 ("rcp" in 
         ! radcomp3)

         if (level == 2) then
           ! Use actual rcp values for argument #22
           CALL radcomp3 (mzp,mxp,myp,ia,iz,ja,jz  &
            ,grid_g(ngrid)%glat       (1,1)    &
            ,grid_g(ngrid)%rtgt       (1,1)    &
            ,grid_g(ngrid)%topt       (1,1)    &
            ,radiate_g(ngrid)%albedt  (1,1)    &
            ,radiate_g(ngrid)%cosz    (1,1)    &
            ,radiate_g(ngrid)%rlongup (1,1)    &
            ,radiate_g(ngrid)%rshort  (1,1)    &
            ,radiate_g(ngrid)%rlong   (1,1)    &
            ,basic_g(ngrid)%rv        (1,1,1)  &
            ,basic_g(ngrid)%dn0       (1,1,1)  &
            ,radiate_g(ngrid)%fthrd   (1,1,1)  &
            ,basic_g(ngrid)%pi0       (1,1,1)  &
            ,basic_g(ngrid)%pp        (1,1,1)  &
            ,basic_g(ngrid)%theta     (1,1,1)  &
            ,micro_g(ngrid)%rcp       (1,1,1)  &
            ,radiate_g(ngrid)%bext    (1,1,1)  &
            ,radiate_g(ngrid)%swup    (1,1,1)  &
            ,radiate_g(ngrid)%swdn    (1,1,1)  &
            ,radiate_g(ngrid)%lwup    (1,1,1)  &
            ,radiate_g(ngrid)%lwdn    (1,1,1)  &
            ,micro_g(ngrid)%cccnp     (1,1,1)  &
            ,micro_g(ngrid)%cccmp     (1,1,1)  &
            ,micro_g(ngrid)%gccnp     (1,1,1)  &
            ,micro_g(ngrid)%gccmp     (1,1,1)  &
            ,micro_g(ngrid)%md1np     (1,1,1)  &
            ,micro_g(ngrid)%md1mp     (1,1,1)  &
            ,micro_g(ngrid)%md2np     (1,1,1)  &
            ,micro_g(ngrid)%md2mp     (1,1,1)  &
            ,micro_g(ngrid)%salt_film_np (1,1,1)  &
            ,micro_g(ngrid)%salt_film_mp (1,1,1)  &
            ,micro_g(ngrid)%salt_jet_np  (1,1,1)  &
            ,micro_g(ngrid)%salt_jet_mp  (1,1,1)  &
            ,micro_g(ngrid)%salt_spum_np (1,1,1)  &
            ,micro_g(ngrid)%salt_spum_mp (1,1,1) )
         else
            ! Use zeros for argument #22 ("rcp" in radcomp3)
            CALL azero (mzp*mxp*myp,scratch%vt3dp(1))
            CALL radcomp3 (mzp,mxp,myp,ia,iz,ja,jz  &
            ,grid_g(ngrid)%glat       (1,1)    &
            ,grid_g(ngrid)%rtgt       (1,1)    &
            ,grid_g(ngrid)%topt       (1,1)    &
            ,radiate_g(ngrid)%albedt  (1,1)    &
            ,radiate_g(ngrid)%cosz    (1,1)    &
            ,radiate_g(ngrid)%rlongup (1,1)    &
            ,radiate_g(ngrid)%rshort  (1,1)    &
            ,radiate_g(ngrid)%rlong   (1,1)    &
            ,basic_g(ngrid)%rv        (1,1,1)  &
            ,basic_g(ngrid)%dn0       (1,1,1)  &
            ,radiate_g(ngrid)%fthrd   (1,1,1)  &
            ,basic_g(ngrid)%pi0       (1,1,1)  &
            ,basic_g(ngrid)%pp        (1,1,1)  &
            ,basic_g(ngrid)%theta     (1,1,1)  &
            ,scratch%vt3dp            (1)      &
            ,radiate_g(ngrid)%bext    (1,1,1)  &
            ,radiate_g(ngrid)%swup    (1,1,1)  &
            ,radiate_g(ngrid)%swdn    (1,1,1)  &
            ,radiate_g(ngrid)%lwup    (1,1,1)  &
            ,radiate_g(ngrid)%lwdn    (1,1,1)  &
            ,micro_g(ngrid)%cccnp     (1,1,1)  &
            ,micro_g(ngrid)%cccmp     (1,1,1)  &
            ,micro_g(ngrid)%gccnp     (1,1,1)  &
            ,micro_g(ngrid)%gccmp     (1,1,1)  &
            ,micro_g(ngrid)%md1np     (1,1,1)  &
            ,micro_g(ngrid)%md1mp     (1,1,1)  &
            ,micro_g(ngrid)%md2np     (1,1,1)  &
            ,micro_g(ngrid)%md2mp     (1,1,1)  &
            ,micro_g(ngrid)%salt_film_np (1,1,1)  &
            ,micro_g(ngrid)%salt_film_mp (1,1,1)  &
            ,micro_g(ngrid)%salt_jet_np  (1,1,1)  &
            ,micro_g(ngrid)%salt_jet_mp  (1,1,1)  &
            ,micro_g(ngrid)%salt_spum_np (1,1,1)  &
            ,micro_g(ngrid)%salt_spum_mp (1,1,1) )
         endif

      endif

   endif

endif

return
END SUBROUTINE radiate

!##############################################################################
Subroutine tend_accum (m1,m2,m3,ia,iz,ja,jz,at,at2)

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,i,j,k
real, dimension(m1,m2,m3) :: at,at2

do j = ja,jz
   do i = ia,iz
      do k = 1,m1
         at(k,i,j) = at(k,i,j) + at2(k,i,j)
      enddo
   enddo
enddo

return
END SUBROUTINE tend_accum

!##############################################################################
Subroutine radprep (m2,m3,mzg,mzs,np,ia,iz,ja,jz,jday   &
   ,soil_water       ,soil_energy      ,soil_text      &
   ,sfcwater_energy  ,sfcwater_depth   ,leaf_class     &
   ,veg_fracarea     ,veg_height       ,veg_albedo     &
   ,patch_area                                         &
   ,sfcwater_nlev    ,veg_temp         ,can_temp       & 
   ,solfac,glat,glon,rshort,rlong,rlongup,albedt,cosz  &
   ,sfcswa,uplwrf)

use rconstants
use mem_leaf, only: isfcl
use mem_grid, only: initial,time

implicit none

integer :: m2,m3,mzg,mzs,np,ia,iz,ja,jz,jday,runsfcrad
real :: solfac
real, dimension(m2,m3) :: glat,glon,rshort,rlong,rlongup,albedt,cosz
real, dimension(mzg,m2,m3,np) :: soil_water,soil_energy,soil_text
real, dimension(mzs,m2,m3,np) :: sfcwater_energy,sfcwater_depth
real, dimension(m2,m3,np) :: leaf_class,veg_fracarea,veg_height,veg_albedo  &
   ,patch_area,sfcwater_nlev,veg_temp,can_temp,sfcswa,uplwrf

integer :: ip,i,j

! Compute solar zenith angle [cosz(i,j)] & solar constant factr [solfac].

CALL zen (m2,m3,ia,iz,ja,jz,jday,glat,glon,cosz,solfac)

! Compute patch-averaged surface albedo [albedt(i,j)] and up longwave
! radiative flux [rlongup(i,j)].

CALL azero2 (m2*m3,albedt,rlongup)

do ip = 1,np

   !Only run sfcrad if using LEAF3 or it is the first model timestep
   !or if its a water patch.
   !If running SIB land surface, use SiB surface albedo and surface
   !upward longwave after first timestep.
   runsfcrad=0
   if    (isfcl<=1) then
      runsfcrad=1
   elseif(isfcl==2) then
      if(ip==1 .or. (initial<=2 .and. time < 0.001)) runsfcrad=1
   endif

   do j = ja,jz
      do i = ia,iz

       if(runsfcrad==1) &
         CALL sfcrad (mzg,mzs,ip                                  &
          ,soil_energy    (1,i,j,ip) ,soil_water      (1,i,j,ip)  &
          ,soil_text      (1,i,j,ip) ,sfcwater_energy (1,i,j,ip)  &
          ,sfcwater_depth (1,i,j,ip) ,patch_area      (i,j,ip)    &
          ,can_temp       (i,j,ip)   ,veg_temp        (i,j,ip)    &
          ,leaf_class     (i,j,ip)   ,veg_height      (i,j,ip)    &
          ,veg_fracarea   (i,j,ip)   ,veg_albedo      (i,j,ip)    &
          ,sfcwater_nlev  (i,j,ip)   ,rshort          (i,j)       &
          ,rlong          (i,j)      ,albedt          (i,j)       &
          ,rlongup        (i,j)      ,cosz            (i,j)       )

       if(isfcl==2 .and. ip>=2) then
         albedt(i,j)  = albedt(i,j)  + patch_area(i,j,ip)*sfcswa(i,j,ip)
         rlongup(i,j) = rlongup(i,j) + patch_area(i,j,ip)*uplwrf(i,j,ip)
       endif

      enddo
   enddo
enddo

return
END SUBROUTINE radprep

!##############################################################################
Subroutine radcomp (m1,m2,m3,ia,iz,ja,jz,solfac  &
   ,theta,pi0,pp,rv,dn0,rtp,fthrd  &
   ,rtgt,f13t,f23t,glon,rshort,rlong,albedt,cosz,rlongup)

use mem_grid
use mem_scratch
use mem_radiate
use rconstants

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,jday,oyr,omn,ody,otm
integer, external :: julday

real :: solfac,cdec,declin,dzsdx,dzsdy,dlon,a1,a2,dayhr,gglon,otmf  &
   ,dayhrr,hrangl,sinz,sazmut,slazim,slangl,cosi,pisolar,eqt_julian_nudge
real, dimension(nzpmax) :: rvr,rtr,dn0r,pird,prd,fthrl,dzmr,dztr,fthrs
real, dimension(nzpmax+1) :: temprd
real, dimension(m1,m2,m3) :: theta,pi0,pp,rv,dn0,rtp,fthrd
real, dimension(m2,m3) :: rtgt,f13t,f23t,glon,rshort,rlong,cosz  &
   ,albedt,rlongup
integer :: i,j,k

! Find the hour angle and julian day, then get cosine of zenith angle.
CALL date_add_to (iyear1,imonth1,idate1,itime1*100,time,'s',oyr,omn,ody,otm)
otmf=float(otm)
dayhr = (otmf-mod(otmf,10000.))/10000. + (mod(otmf,10000.)/6000.)
jday = julday(omn,ody,oyr)

! Find day of year equation of time adjustment (Saleeby2008: improve accuracy)
pisolar=3.1415927
eqt_julian_nudge=0.0
if(jday>=1 .and. jday<=106) then
 eqt_julian_nudge = -14.2 * sin(pisolar * (jday +   7) / 111.)
elseif(jday>=107 .and. jday<=166) then
 eqt_julian_nudge =   4.0 * sin(pisolar * (jday - 106) /  59.)
elseif(jday>=167 .and. jday<=246) then
 eqt_julian_nudge =  -6.5 * sin(pisolar * (jday - 166) /  80.)
elseif(jday>=247 .and. jday<=366) then
 eqt_julian_nudge =  16.4 * sin(pisolar * (jday - 247) / 113.)
endif
eqt_julian_nudge = eqt_julian_nudge / 60.0
! cdec - cosine of declination
declin = -23.5 * cos(6.283 / 365. * (jday + 9)) * pi180
cdec = cos(declin)

do j = ja,jz
   do i = ia,iz
      do k = 1,m1
         pird(k) = (pp(k,i,j) + pi0(k,i,j)) / cp
         temprd(k) = theta(k,i,j) * pird(k)
         rvr(k) = max(0.,rv(k,i,j))
         rtr(k) = max(rvr(k),rtp(k,i,j))
! Convert the next 7 variables to cgs for now.
         prd(k) = pird(k) ** cpor * p00 * 10.
         dn0r(k) = dn0(k,i,j) * 1.e-3
         dzmr(k) = dzm(k) / rtgt(i,j) * 1.e-2
         dztr(k) = dzt(k) / rtgt(i,j) * 1.e-2

         if(rvr(k)    <   0. .or.  &
            rtr(k)    <   0. .or.  &
            prd(k)    <   0. .or.  &
            dn0r(k)   <   0. .or.  &
            temprd(k) < 173.) then
            print*, 'Temperature too low or negative value of'
            print*, 'density, vapor, pressure, or ozone'
            print*, 'when calling Chen-Cotton/Mahrer-Pielke radiation'
            print*, 'at (k,i,j),ngrid = ',k,i,j,ngrid
            print*, 'rvr(k)  rtr(k)  prd(k)  dn0r(k)  temprd(k)'
            print*, rvr(k), rtr(k), prd(k), dn0r(k), temprd(k)
            print*, 'stopping model'
            stop 'radiation call'
         endif

      enddo
      temprd(1) = (rlongup(i,j) / stefan) ** 0.25
      temprd(m1+1) = temprd(m1)

! Call the longwave parameterizations.

      if (ilwrtyp .eq. 2) then
         CALL lwradp (m1,temprd,rvr,dn0r,dztr,pird,vctr1,fthrl,rlong(i,j))
      elseif (ilwrtyp .eq. 1) then
         CALL lwradc (m1+1,rvr,rtr,dn0r,temprd,prd,pird,dztr,fthrl,rlong(i,j)  &
            ,vctr1,vctr2,vctr3,vctr4,vctr5,vctr6,vctr7,vctr8,vctr9,vctr10  &
            ,vctr11,vctr12,vctr13,vctr14,vctr15)
      endif

! The shortwave parameterizations are only valid if the cosine
!    of the zenith angle is greater than .03 .

      if (cosz(i,j) .gt. .03) then

         if (iswrtyp .eq. 2) then
            CALL shradp (m1,rvr,dn0r,dzmr,vctr1,pird,cosz(i,j)  &
               ,albedt(i,j),solar*1e3*solfac,fthrs,rshort(i,j))
         elseif (iswrtyp .eq. 1) then
            CALL shradc (m1+1,rvr,rtr,dn0r,dztr,prd,pird,vctr1  &
              ,albedt(i,j),solar*1.e3*solfac,cosz(i,j),fthrs,rshort(i,j))
         endif

! Modify the downward surface shortwave flux by considering
!    the slope of the topography.

         if (itopo .eq. 1) then
            dzsdx = f13t(i,j) * rtgt(i,j)
            dzsdy = f23t(i,j) * rtgt(i,j)

! The y- and x-directions must be true north and east for
! this correction. the following rotates the model y/x
! to the true north/east.

!Saleeby(2013):Check into this statement below
! The following rotation seems to be incorrect, so call this instead:
! routine uvtoueve (u,v,ue,ve,qlat,qlon,polelat,polelon)

            dlon = (polelon - glon(i,j)) * pi180
            a1 = dzsdx*cos(dlon) + dzsdy * sin(dlon)
            a2 = -dzsdx*sin(dlon) + dzsdy * cos(dlon)
            dzsdx = a1
            dzsdy = a2

            gglon = glon(i,j)
            if (lonrad .eq. 0) gglon = centlon(1)
!Saleeby(2009): Add eqt_julian_nudge to improve hour angle accuracy
            dayhrr = mod(dayhr+gglon/15.+24.,24.) + eqt_julian_nudge
            hrangl = 15. * (dayhrr - 12.) * pi180
            sinz = sqrt(1. - cosz(i,j) ** 2)
            sazmut = asin(max(-1.,min(1.,cdec*sin(hrangl)/sinz)))
            if (abs(dzsdx) .lt. 1e-15) dzsdx = 1.e-15
            if (abs(dzsdy) .lt. 1e-15) dzsdy = 1.e-15
            slazim = 1.571 - atan2(dzsdy,dzsdx)
            slangl = atan(sqrt(dzsdx*dzsdx+dzsdy*dzsdy))
            cosi = cos(slangl) * cosz(i,j) + sin(slangl) * sinz  &
               * cos(sazmut-slazim)
!Saleeby(08-10-2008): Check for cosi greater than zero
            if (cosi > 0.) then
               rshort(i,j) = rshort(i,j) * cosi / cosz(i,j)
            else
               rshort(i,j) =  0.
            endif
         endif

      else
         do k = 1,nzp
            fthrs(k) = 0.
         enddo
         rshort(i,j) = 0.
      endif
      
      ! Add fluxes
      do k = 2,m1-1
         fthrd(k,i,j) = fthrl(k) + fthrs(k)
      enddo

! Convert the downward flux at the ground to SI.

      rshort(i,j) = rshort(i,j) * 1.e-3 / (1. - albedt(i,j))
      rlong(i,j) = rlong(i,j) * 1.e-3
      fthrd(1,i,j) = fthrd(2,i,j)

   enddo
enddo

return
END SUBROUTINE radcomp

!##############################################################################
Subroutine radcomp3 (m1,m2,m3,ia,iz,ja,jz  &
   ,glat,rtgt,topt,albedt,cosz,rlongup,rshort,rlong  &
   ,rv,dn0,fthrd,pi0,pp,theta,rcp &
   ,bext,swup,swdn,lwup,lwdn &
   ,cccnp,cccmp,gccnp,gccmp,md1np,md1mp,md2np,md2mp &
   ,salt_film_np,salt_film_mp,salt_jet_np,salt_jet_mp &
   ,salt_spum_np,salt_spum_mp)

use mem_grid
use mem_micro
use mem_radiate
use rconstants
use rrad3
use micphys

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,mcat,i,j,k

real :: cfmasi,cparmi,glg,glgm,picpi
real, dimension(m2,m3) :: glat,rtgt,topt,cosz,albedt,rlongup,rshort,rlong
real, dimension(m1,m2,m3) :: dn0,rv,fthrd,pi0,pp,theta,rcp
real, dimension(m1,m2,m3) :: bext,swup,swdn,lwup,lwdn
real, dimension(m1,m2,m3) :: cccnp,cccmp,gccnp,gccmp,md1np,md1mp,md2np,md2mp &
  ,salt_film_np,salt_film_mp,salt_jet_np,salt_jet_mp,salt_spum_np,salt_spum_mp
real, external :: gammln

! Fill cloud parameters if not running microphysics
cparmi = 1. / cparm
if (level <= 1) then
   mcat = 0
elseif (level == 2) then
   mcat = 1
   pwmas(1) = 3.
   pwmasi(1) = 1. / pwmas(1)
   cfmasi = 1. / 524.
   emb0(1) = 5.24e-16
   emb1(1) = 3.35e-11
   glg  = gammln(gnu(1))
   glgm = gammln(gnu(1) + pwmas(1))
   dnfac(1) = (cfmasi * exp(glg - glgm)) ** pwmasi(1)
   do k = 2,m1-1
      jhcat(k,1) = 1
   enddo
endif

! Loop over columns

do j = ja,jz
   do i = ia,iz

      ! To be used in mclatchy call
      do k = 2,m1-1
         picpi = (pi0(k,i,j) + pp(k,i,j)) * cpi
         press(k) = p00 * picpi ** cpor
         tair(k) = theta(k,i,j) * picpi
      enddo

      ! Finish cloud parameter computation
      if (level == 2) then
         do k = 2,m1-1
            emb(k,1) = max(emb0(1),min(emb1(1),rcp(k,i,j) * cparmi))
            cx(k,1) = rcp(k,i,j) / emb(k,1)
         enddo
      elseif (level == 3) then
         CALL cloud_prep (m1,i,j,ngrid,dn0(1,i,j))
         mcat = 7 !Saleeby(2008): Increase to 8 for adding drizzle
      endif

      ! Call the sub-driver
      if(iaerorad==1)then
        CALL aero_copy (1,m1 &
           ,cccnp(1,i,j),cccmp(1,i,j) &
           ,gccnp(1,i,j),gccmp(1,i,j) &
           ,md1np(1,i,j),md1mp(1,i,j) &
           ,md2np(1,i,j),md2mp(1,i,j) &
           ,salt_film_np(1,i,j),salt_film_mp(1,i,j) &
           ,salt_jet_np(1,i,j) ,salt_jet_mp(1,i,j)  &
           ,salt_spum_np(1,i,j),salt_spum_mp(1,i,j))
      endif

      CALL radcalc3 (m1,maxnzp,mcat,iswrtyp,ilwrtyp  &
         ,glat(i,j),rtgt(i,j),topt(i,j),albedt(i,j),cosz(i,j)  &
         ,rlongup(i,j),rshort(i,j),rlong(i,j)  &
         ,zm,zt,rv(1,i,j),dn0(1,i,j),pi0(1,i,j),pp(1,i,j) &
         ,fthrd(1,i,j),i,j,ngrid &
         ,bext(1,i,j),swup(1,i,j),swdn(1,i,j),lwup(1,i,j),lwdn(1,i,j))

   enddo
enddo

return
END SUBROUTINE radcomp3

!##############################################################################
Subroutine zen (m2,m3,ia,iz,ja,jz,jday,glat,glon,cosz,solfac)

use mem_grid
use mem_radiate
use rconstants

implicit none

integer :: m2,m3,ia,iz,ja,jz,jday,i,j,oyr,omn,ody,otm
integer, external :: julday

real :: solfac,sdec,cdec,declin,d0,d02,dayhr,radlat,cslcsd,snlsnd  &
   ,gglon,dayhrr,hrangl,pisolar,eqt_julian_nudge,otmf
real, dimension(m2,m3) :: glat,glon,cosz

! Find the hour angle and julian day, then get cosine of zenith angle.
CALL date_add_to (iyear1,imonth1,idate1,itime1*100,time,'s',oyr,omn,ody,otm)
otmf=float(otm)
dayhr = (otmf-mod(otmf,10000.))/10000. + (mod(otmf,10000.)/6000.)
jday = julday(omn,ody,oyr)

!      sdec - sine of declination, cdec - cosine of declination
declin = -23.5 * cos(6.283 / 365. * (jday + 9)) * pi180
sdec = sin(declin)
cdec = cos(declin)

! Find the factor, solfac, to multiply the solar constant to correct
! for Earth's varying distance to the sun.

d0 = 6.2831853 * float(jday-1) / 365.
d02 = d0 * 2.
solfac = 1.000110 + 0.034221 * cos (d0) + 0.001280 * sin(d0)  &
   + 0.000719 * cos(d02) + 0.000077 * sin(d02)

! Find day of year equation of time adjustment (Saleeby2008: improve accuracy)
pisolar=3.1415927
eqt_julian_nudge=0.0
if(jday>=1 .and. jday<=106) then
 eqt_julian_nudge = -14.2 * sin(pisolar * (jday +   7) / 111.)
elseif(jday>=107 .and. jday<=166) then
 eqt_julian_nudge =   4.0 * sin(pisolar * (jday - 106) /  59.)
elseif(jday>=167 .and. jday<=246) then
 eqt_julian_nudge =  -6.5 * sin(pisolar * (jday - 166) /  80.)
elseif(jday>=247 .and. jday<=366) then
 eqt_julian_nudge =  16.4 * sin(pisolar * (jday - 247) / 113.)
endif
eqt_julian_nudge = eqt_julian_nudge / 60.0

do j = ja,jz
   do i = ia,iz
      radlat = glat(i,j) * pi180
      if (lonrad .eq. 0) radlat = centlat(1) * pi180
      if (radlat .eq. declin) radlat = radlat + 1.e-5
      cslcsd = cos(radlat) * cdec
      snlsnd = sin(radlat) * sdec
      gglon = glon(i,j)
      if (lonrad .eq. 0) gglon = centlon(1)
!Saleeby(2009): Add eqt_julian_nudge to improve hour angle accuracy
      dayhrr = mod(dayhr+gglon/15.+24.,24.) + eqt_julian_nudge
      hrangl = 15. * (dayhrr - 12.) * pi180
      cosz(i,j) = snlsnd + cslcsd * cos(hrangl)
   enddo
enddo

return
END SUBROUTINE zen

!##############################################################################
Subroutine radcalc3 (m1,maxnzp,mcat,iswrtyp,ilwrtyp  &
   ,glat,rtgt,topt,albedt,cosz,rlongup,rshort,rlong  &
   ,zm,zt,rv,dn0,pi0,pp,fthrd,i,j,ngrid &
   ,bext,swup,swdn,lwup,lwdn)

!-----------------------------------------------------------------------------
! radcalc3: column driver for twostream radiation code
! variables used within routine radcalc3:
! ==================================================================
! Variables in rrad3 parameter statement
!  mb               : maximum allowed number of bands [=8]
!  mg                  : maximum allowed number of gases [=3]
!  mk               : maximum number of pseudobands allowed for any gas [=7]
!  ncog             : number of fit coefficients (omega and asym) [=5]
!  ncb              : number of fit coefficients (extinction) [=2]
!  npartob          : number of hydrometeor categories (including different habits)
!  npartg           : number of hydrometeor categories used for gc coefficients [=7]
!  nrad                  : total number of radiation levels used (m1 - 1 + narad)
!  narad            : number of radiation levels added above model
!  nsolb            : active number of solar bands
!  nb               : active number of bands
!  ng                : active number of gases
!  jday             : julian day
!  solfac           : solar constant multiplier for variable E-S distance
!  ralcs (mb)       : rayleigh scattering integration constants
!  solar1 (mb)      : solar fluxes at top of atmosphere - corrected for ES distance
!  solar0 (mb)      : solar fluxes at top of atmosphere - uncorrected for ES distance
!  nuum (mb)        :    continuum flags
!  a0,a1,a2,a3 (mb) : Planck func fit coefficients
!  npsb (mg,mb)     : number of pseudo bands
!  trf (mg,mb)      : reference temperature for xp and wght coefficients
!  prf (mg,mb)      : reference pressure for xp and wght coefficients
!  ulim (mg,mb)     : upper bound on pathlength for gases
!  xp (mg,mk,mb)    : coefficient used in computing gaseous optical depth
!  alpha (mg,mk,mb) : pressure correction factor exponent
!  beta (mg,mk,mb)  : temperature correction factor exponent
!  wght (mg,mk,mb)  : pseudo band weight
!  exptabc (150)    : table of exponential func values
!  ocoef(ncog,mb,npartob)  : fit coefficients for hyd. single scatter.
!  bcoef(ncb,mb ,npartob)  : fit coefficients for hyd. extinction coefficient.
!  gcoef(ncog,mb,npartg)   : fit coefficients for hyd. asymmetry parameter.
!
! Input variables from model
!
!  m1               : number of vertical levels in model grid
!  ncat             : max number of hydrometeor categories [=7]
!  mcat             : actual number of hydrometeor categories [= 0, 1, or 7]
!  nhcat            : number of hydrometeor categories including ice habits [=15]
!  iswrtyp          : shortwave radiation parameterization selection flag
!  ilwrtyp          : longwave radiation parameterization selection flag
!  glat             : latitude
!  rtgt             : terrain-following coordinate metric factor
!  topt             : topography height
!  albedt          : surface albedo
!  cosz             : solar zenith angle
!  rlongup          : upward longwave radiation at surface (W/m^2)
!  rshort           : downward shortwave radiation at surface (W/m^2)
!  rlong            : downward longwave radiation at surface (W/m^2)
!  jnmb (ncat)      : microphysics category flag
!  dnfac (nhcat)    : factor for computing dn from emb
!  pwmasi (nhcat)   : inverse of mass power law exponent for hydrometeors
!  zm (m1)          : model physical heights of W points (m)
!  zt (m1)          : model physical heights of T points (m)
!  press (nzpmax)   : model pressure (Pa)
!  tair (nzpmax)    : model temperature (K)
!  rv (m1)          : model vapor mixing ratio (kg/kg)
!  dn0 (m1)         : model air density (kg/m^3)
!  fthrd (m1)       : theta_il tendency from radiation
!  jhcat (nzpmax,ncat)  : microphysics category array
!  cx (nzpmax,ncat) : hydrometeor number concentration (#/kg)
!  emb (nzpmax,ncat): hydrometeor mean mass (kg)
!
! Variables input from model scratch space (redefined internally on each call)
!
!  zml (nrad)       : physical heights of W points of all radiation levels (m)
!  ztl (nrad)       : physical heights of T points of all radiation levels (m)
!  dzl (nrad)       : delta-z (m) of all radiation levels
!  pl (nrad)        : pressure (Pa)
!  tl (nrad)        : temperature (K)
!  dl (nrad)        : air density of all radiation levels (kg/m^3)
!  rl (nrad)        : vapor density of all radiation levels (kg/m^3)
!  vp (nrad)        : vapor pressure (Pa)
!  o3l (nrad)       : stores the calculated ozone profile (g/m^3)
!  flxu (nrad)      : Total upwelling flux (W/m^2)
!  flxd (nrad)      : Total downwelling flux (W/m^2)
!  t (nrad)         : layer transmission func
!  r (nrad)         : layer reflection func
!  tc (nrad)        : cumulative optical depth
!  sigu (nrad)      : upwelling layer source func
!  sigd (nrad)      : downwelling layer source func
!  re (nrad)        : cumulative reflection func
!  vd (nrad)        : multi-scat. diffuse downwelling contributions
!                         from source func
!  td (nrad)        : inverse of cumulative transmission fnct
!  vu (nrad)        : multi-scat. diffuse upwelling contributions
!                         from source func
!  tg (nrad)        : gaseous optical depth
!  tcr (nrad)       : continuum/Rayleigh optical depth
!  src (nrad)       : Planck func source for each band
!  fu (nrad*6)      : upwelling fluxes for pseudo-bands (W/m^2)
!  fd (nrad*6)      : downwelling fluxes for pseudo-bands (W/m^2)
!  u (nrad*mg)      : path-length for gases (H_2O, CO_2, O_3)  (Pa)
!  tp (nrad*mb)     : optical depth of hydrometeors (m^-1)
!  omgp (nrad*mb)   : Single scatter albedo of hydrometeors
!  gp (nrad*mb)     : Asymmetry factor of hydrometeors
!
! Locally-defined variables
!
!  ngass (mg)       : flags indicating if H20, CO2, O3 are active for solar wavelengths
!  ngast (mg)       : flags indicating if H20, CO2, O3 are active for long wavelengths
!  prsnz,prsnzp     : pressure in top two model reference state levels
!
! Additional variables used only within routine mclatchy:
! ==================================================================
! namax            : maximum allowed number of added rad levels above model top[=10]
!                       used for oc and bc coefficients [=13]
! mcdat (33,9,6)    : Mclatchy sounding data (33 levels, 9 soundings, 6 vars)
! mclat (33,9,6)    : mcdat interpolated by season to latitude bands
! mcol (33,6)       : mclat interpolated to lat-lon of grid column
!
! Additional variables used only within routine cloud_opt:
! ==================================================================
!  ib .......... band number
!  dn .......... characteristic diameter (m)
!  oc .......... scattering albedo fit coefficients
!  bc .......... extinction fit coefficients
!  gc .......... asymmetery fit coefficients
!  kradcat ..... cross-reference table giving Jerry's 13 hydrometeor category
!                   numbers as a func of 15 microphysics category numbers
!
! Particle Numbers describe the following particles:
!
!     Harrington radiation code             RAMS microphysics
! ----------------------------------------------------------------
!  1:   cloud drops                 1.  cloud drops
!  2:   rain                        2.  rain
!  3:   pristine ice columns        3.  pristine ice columns
!  4:   pristine ice rosettes       4.  snow columns
!  5:   pristine ice plates         5.  aggregates
!  6:   snow columns                6.  graupel
!  7:   snow rosettes               7.  hail
!  8:   snow plates                 8.  pristine ice hexagonal plates
!  9:   aggregates columns          9.  pristine ice dendrites
!  10:  aggregates rosettes        10.  pristine ice needles
!  11:  aggregates plates          11.  pristine ice rosettes
!  12:  graupel                    12.  snow hexagonal plates
!  13:  hail                       13.  snow dendrites
!                                  14.  snow needles
!                                  15.  snow rosettes
!
! for the asymmetery parameter, since we only have spherical
! particles, there are only 7 particle types...
!  1:   cloud drops
!  2:   rain
!  3:   pristine ice
!  4:   snow
!  5:   aggregates
!  6:   graupel
!  7:   hail
!---------------------------------------------------------------------------

use rconstants
use rrad3
use micphys

implicit none

integer m1,maxnzp,mcat,ngrid
integer :: iswrtyp,ilwrtyp
integer i,j,k
integer, save :: ncall = 0,nradmax
integer, save :: ngass(mg)=(/1, 1, 1/),ngast(mg)=(/1, 1, 1/)
!     one can choose the gases of importance here,
!       ngas = 1    gas active
!            = 0    gas not active
!
!       ngas(1) = H2O
!       ngas(2) = CO2
!       ngas(3) =  O3

real, save :: eps=1.e-15
real :: prsnz,prsnzp
real :: glat,rtgt,topt,cosz,albedt,rlongup,rshort,rlong
real :: zm(m1),zt(m1),dn0(m1),rv(m1),pi0(m1),pp(m1),fthrd(m1)
real :: bext(m1),swup(m1),swdn(m1),lwup(m1),lwdn(m1)

real, allocatable, save, dimension(:)     :: zml,ztl,dzl,pl,tl,dl,rl,o3l  &
                                      ,vp,flxus,flxds,tg,tcr,src,t,r,tc  &
                                      ,flxul,flxdl  &
                                      ,sigu,sigd,re,vd,td,vu  &
                                      ,u,fu,fd,tp,omgp,gp

real :: exner(m1) ! non-dimensional pressure

!Saleeby(2011):Variables for radiatively active aerosols
real :: relh(m1)
real, external :: rslf

if (ncall == 0) then
   ncall = 1
   nradmax = maxnzp + namax
   allocate(zml  (nradmax) ,ztl  (nradmax) ,dzl  (nradmax) ,pl (nradmax)  &
           ,tl   (nradmax) ,dl   (nradmax) ,rl   (nradmax) ,o3l(nradmax)  &
           ,vp   (nradmax) ,flxus(nradmax) ,flxds(nradmax) ,tg (nradmax)  &
           ,flxul(nradmax) ,flxdl(nradmax)                                &
           ,tcr  (nradmax) ,src  (nradmax) ,t    (nradmax) ,r  (nradmax)  &
           ,tc   (nradmax) ,sigu (nradmax) ,sigd (nradmax) ,re (nradmax)  &
           ,vd   (nradmax) ,td   (nradmax) ,vu   (nradmax))
   allocate(u(nradmax*mg),fu(nradmax*6),fd(nradmax*6))
   allocate(tp(nradmax*mb),omgp(nradmax*mb),gp(nradmax*mb))
   tg=0.
endif

nrad = m1 - 1 + narad

! rlongup used to set tl(1): stephan*tl^4=rlongup
 CALL mclatchy (3,m1  &
   ,prsnz,prsnzp,glat,rtgt,topt,rlongup  &
   ,zm,zt,press,tair,dn0,rv,zml,ztl,pl,tl,dl,rl,o3l,dzl)

! calculate non-dimensional pressure
do k=1,m1
  exner(k) = (pi0(k)+pp(k))/cp
enddo

! zero out scratch arrays
 CALL azero (nrad*mg,u)
 CALL azero (nrad*6,fu)
 CALL azero (nrad*6,fd)
 CALL azero (nrad*mb,tp)
 CALL azero (nrad*mb,omgp)
 CALL azero (nrad*mb,gp)
 CALL azero (nrad   ,tg)

! Saleeby(2011): Aerosol radiative impacts section
! This must be run before cloud_opt
! Only run this for level=3 microphysics
if(iaerorad==1) then
 do k=1,m1
   relh(k) = rv(k)/rslf(pl(k),tl(k))
 enddo
 CALL aerorad (mb,nb,nrad,m1,dzl,relh,tp,omgp,gp,dn0)
endif

! Compute hydrometeor optical properties
 CALL cloud_opt (mb,nb,nrad,m1,mcat,dzl  &
   ,dn0,tp,omgp,gp &
   ,ocoef,bcoef,gcoef,ncog,ncb,npartob,npartg)

! Sum up attenutation by aerosols and hydrometeors
 CALL sum_opt (mb,nrad,nb,m1,tp,omgp,gp,bext,dzl)

! Get the path lengths for the various gases...
 CALL path_lengths (nrad,u,rl,dzl,dl,o3l,vp,pl,eps)

do k = 1,nrad
   if (rl(k) <   0. .or.  &
       dl(k) <   0. .or.  &
       pl(k) <   0. .or.  &
      o3l(k) <   0. .or.  &
       tl(k) < 160.) then

      print*, 'Temperature too low or negative value of'
      print*, 'density, vapor, pressure, or ozone'
      print*, 'when calling Harrington radiation'
      print*, 'at k,i,j = ',k,i,j,'   ngrid=',ngrid
      print*, 'stopping model'
      print*, 'rad: k, rl(k), dl(k), pl(k), o3l(k), tl(k)'
      print'(i4,5g15.6)', k, rv(k), dl(k), pl(k), o3l(k), tl(k)
      stop 'stop: radiation call'
   endif
enddo

! call shortwave and longwave schemes...

if (iswrtyp == 3 .and. cosz > 0.03) then
   CALL azero2 (nrad,flxus,flxds)

   CALL swrad (nrad,nb,nsolb,npsb,   &         !  counters
      u,pl,tl,dzl,                  &      !  model variables
      xp,alpha,beta,wght,prf,trf,ralcs,  &   !  band specifics
      solar1,ngass,                      &    !        "
      albedt,cosz,               & !  boundaries
      tp,omgp,gp,fu,fd,         & !       "
      flxus,flxds,ulim)                     !  sw fluxes

   rshort = flxds(1)

   do k = 2,m1-1
      !divide by exner to get potential temp heating rate
      fthrd(k) = fthrd(k)  &
         + (flxds(k) - flxds(k-1) + flxus(k-1) - flxus(k)) &
            / (dl(k) * dzl(k) * cp * exner(k))
      swup(k) = flxus(k)
      swdn(k) = flxds(k)
    enddo
    !lower and upper boundary conditions on swup and swdn
    swup(1) = flxus(1)
    swup(m1) = flxus(nrad) ! use the top radiation value rather than m1 value here
    swdn(1) = flxds(1)
    swdn(m1) = flxds(nrad) ! use the top radiation value rather than m1 value here

else

   swup   = 0.
   swdn   = 0.
   rshort = 0.

endif

if (ilwrtyp == 3) then
   CALL azero2 (nrad,flxul,flxdl)

   CALL lwrad (nrad,nb,nsolb,npsb,nuum,   &    !  counters
      u,pl,tl,vp,                       & !  model variables
      xp,alpha,beta,wght,prf,trf,ralcs,     &!  band specifics
      a0,a1,a2,a3,                          &!        "
      exptabc,ngast,                        &!  boundaries
      tp,omgp,gp,fu,fd,flxul,flxdl,ulim)  !  fluxes, output

   rlong = flxdl(1)

   do k = 2,m1-1
      !divide by exner to get potential temp heating rate
      fthrd(k) = fthrd(k)  &
         + (flxdl(k) - flxdl(k-1) + flxul(k-1) - flxul(k)) &
            / (dl(k) * dzl(k) * cp * exner(k))
      lwup(k) = flxul(k)
      lwdn(k) = flxdl(k) 
   enddo
   !lower and upper boundary conditions on lwup and lwdn
   lwup(1) = flxul(1)
   lwup(m1) = flxul(nrad) ! use the top radiation value rather than m1 value here
   lwdn(1) = flxdl(1)
   lwdn(m1) = flxdl(nrad) ! use the top radiation value rather than m1 value here

endif

return
END SUBROUTINE radcalc3

!##############################################################################
Subroutine cloud_opt (mb,nb,nrad,m1,mcat,dzl  &
   ,dn0,tp,omgp,gp,oc,bc,gc,ncog,ncb,npartob,npartg)

! computing properties of spherical liquid water and irregular ice
! using fits to adt theory
!
! ib .......... band number
! mb .......... maximum number of bands
! nb .......... total number of bands
! m1 .......... number of vertical levels
! dzl .......... delta z in each level (m)
! dn .......... characteristic diameter (m)
! emb ......... mean hydrometeor mass (kg)
! cx .......... hydrometeor concentration (#/kg)
! tp .......... optical depth
! omgp ........ scattering albedo
! gp .......... asymmetry parameter
! oc .......... scattering albedo fit coefficients
! bc .......... extinction fit coefficients
! gc .......... asymmetry fit coefficients
! ncog ........ number of fit coefficients (omega and asym)
! ncb ......... number of fit coefficients (extinction)
! kradcat ..... cross-reference table giving Jerry's 13 hydrometeor category
!                 numbers as a func of 15 microphysics category numbers
! dn0 ......... model air density (kg/m^3)
! dnfac ....... factor for computing dn from emb
! pwmasi ...... inverse of power used in mass power law
! npartob ..... number of hydrometeor categories (including different habits)
!                 used for oc and bc coefficients
! npartg ...... number of hydrometeor categories used for gc coefficients
!
! Saleeby(2008): would need to modify dnmin,dnmax,kradcat for drizzle

use micphys

implicit none

integer mb,nb,ib,nrad,m1,ncog,ncb,krc,npartob,npartg
integer icat,mcat,k,ihcat

integer kradcat(15)
real dzl(nrad),tp(nrad,mb),omgp(nrad,mb),gp(nrad,mb),dn0(m1)
real oc(ncog,mb,npartob),bc(ncb,mb,npartob),gc(ncog,mb,npartg)
real ext,om,gg,dn

real dnmin(7),dnmax(7)
data dnmin /   1.,   10.,   1.,  125.,   10.,   10.,   10./
data dnmax /1000.,10000., 125.,10000.,10000.,10000.,10000./

data kradcat/1,2,3,6,10,12,13,5,5,3,4,8,8,6,7/

do icat = 1,mcat
   if (jnmb(icat) .gt. 0) then

      do k = 2,m1-1

         if (cx(k,icat) .gt. 1.e-9) then
            ihcat = jhcat(k,icat)
            krc = kradcat(ihcat)
            dn = dnfac(ihcat) * emb(k,icat) ** pwmasi(ihcat) * 1.e6
            dn = max(dnmin(icat),min(dnmax(icat),dn))

            do ib = 1,nb

               ext = cx(k,icat) * dn0(k) * dzl(k)  &
                  * bc(1,ib,krc) * dn ** bc(2,ib,krc)
               om = oc(1,ib,krc)  &
                  + oc(2,ib,krc) * exp(oc(3,ib,krc) * dn)  &
                  + oc(4,ib,krc) * exp(oc(5,ib,krc) * dn)
               gg = gc(1,ib,icat)  &
                  + gc(2,ib,icat) * exp(gc(3,ib,icat) * dn)  &
                  + gc(4,ib,icat) * exp(gc(5,ib,icat) * dn)

               tp(k,ib)   = tp(k,ib)   + ext
               omgp(k,ib) = omgp(k,ib) + om * ext
               gp(k,ib)   = gp(k,ib)   + gg * om * ext

            enddo

         endif
      enddo

   endif
enddo

return
END SUBROUTINE cloud_opt

!##############################################################################
Subroutine sum_opt (mb,nrad,nb,m1,tp,omgp,gp,bext,dzl)

implicit none

integer nb,ib,m1,k,mb,nrad
real tp(nrad,mb),omgp(nrad,mb),gp(nrad,mb),bext(m1),dzl(m1)

! Combine the optical properties....

do ib = 1,nb
   do k = 2,m1-1

      if (tp(k,ib) .gt. 0.0) then
         gp(k,ib) = gp(k,ib) / omgp(k,ib)
         !Saleeby(2010): Need this min/max func to prevent 'gp'
         ! from being unphysical. Needed for RCE simulations.
         gp(k,ib) = MAX(MIN(gp(k,ib),1.),0.)
         omgp(k,ib) = omgp(k,ib) / tp(k,ib)
         !Saleeby(2010): Need this min/max func to prevent 'omgp'
         ! from being unphysical. Needed for RCE simulations.
         omgp(k,ib) = MAX(MIN(omgp(k,ib),1.),0.)
      else
         omgp(k,ib) = 0.0
         gp(k,ib) = 0.0
      endif

      !Check for validity of opt values before calling radiation
      if (tp(k,ib) .lt. 0) then
         print*, 'tp(k,ib) less than zero for k,ib = ',k,ib
         print*, 'tp(k,ib) = ',tp(k,ib)
         stop 'opt1'
      endif
      if (omgp(k,ib) .lt. 0. .or. omgp(k,ib) .gt. 1.) then
         print*, 'omgp(k,ib) out of range [0,1] for k,ib = ',k,ib
         print*, 'omgp(k,ib) = ',omgp(k,ib)
         stop 'opt2'
      endif
      if (gp(k,ib) .lt. 0. .or. gp(k,ib) .gt. 1.) then
         print*, 'gp(k,ib) out of range [0,1] for k,ib = ',k,ib
         print*, 'gp(k,ib) = ',gp(k,ib)
         stop 'opt3'
      endif

   enddo
enddo

! Calculating visual range (km) using the Koschmeider equation
! Units: tp(k,ib) : optical thickness in level k, band ib (dimensionless)
!        dzl(k) : [m]
!        final bext(k) : [km]
! Consider only over band #3 (245-700 nm)
 do k = 2,m1-1
  bext(k) = 0.
  bext(k) = tp(k,3) / dzl(k) !Compute extinction coefficient
  !Prevent infinite visibility: constrain max vis to 1000km
  if(bext(k) .lt. 3.912e-6) bext(k) = 3.912e-6
  bext(k) = 3.912/bext(k)/1000.
 enddo

return
END SUBROUTINE sum_opt

!##############################################################################
Subroutine path_lengths (nrad,u,rl,dzl,dl,o3l,vp,pl,eps)

! Get the path lengths for the various gases...

implicit none

integer :: nrad
real, dimension(nrad) :: rl,dzl,dl,o3l,vp,pl
real, dimension(nrad,3) :: u
real :: rvk0,rvk1,dzl9,rmix,eps
integer :: k

u(1,1) = .5 * (rl(2) + rl(1)) * 9.81 * dzl(1)
u(1,2) = .5 * (dl(2) + dl(1)) * .45575e-3 * 9.81 * dzl(1)
u(1,3) = o3l(1) * 9.81 * dzl(1)

rvk0 = rl(1)
do k = 2,nrad
   rvk1 = (rl(k) + 1.e-6)
   dzl9 = 9.81 * dzl(k)
   rmix = rvk1 / dl(k)
   vp(k) = pl(k) * rmix / (.622 + rmix)
   u(k,1) = (rvk1 - rvk0) / (log(rvk1 / rvk0) + eps) * dzl9
   u(k,2) = (dl(k) - dl(k-1)) / (log(dl(k) / dl(k-1)) + eps)  &
       * dzl9 * 0.45575e-3
   u(k,3) = 0.5 * dzl9 * (o3l(k) + o3l(k-1))
   rvk0 = rvk1
enddo

return
END SUBROUTINE path_lengths

!##############################################################################
Subroutine cloud_prep (m1,i,j,ng,dn0)

use micphys
use mem_micro

implicit none

integer :: m1,i,j,ng
real :: dn0(m1)

logical, save :: first_call=.true.
integer :: k,lcat

! Prepare cloud arrays for radiation. Requires LEVEL == 3

! Call micro initialization

if (first_call) then

   first_call = .false.
   CALL micinit ()
   
endif

!=========================================================
! Code from each_call.  Define constant emb for species with specified diameter.

do lcat = 1,7  !Saleeby(2008): Increase to 8 for adding drizzle
   if (jnmb(lcat) == 2) then
      do k = 2,m1-1
         emb(k,lcat) = cfmas(lcat) * parm(lcat) ** pwmas(lcat)
      enddo
   endif
   do k = 2,m1-1
      jhcat(k,lcat) = lcat
   enddo
enddo

!==========================================================
! Code from range_check: Fill mixing ratio (rx) and number conc (cx) arrays

do lcat = 1,7 !Saleeby(2008): Increase to 8 for adding drizzle
   do k = 2,m1-1
      rx(k,lcat) = 0.
      cx(k,lcat) = 0.
   enddo

   if (jnmb(lcat) >= 3) then
      do k = 2,m1-1
         emb(k,lcat) = 0.
      enddo
   endif

enddo


! fill scratch arrays for cloud water

if (jnmb(1) >= 1) then
   do k = 2,m1-1
         rx(k,1) = micro_g(ng)%rcp(k,i,j)
         if (jnmb(1) >= 5) cx(k,1) = micro_g(ng)%ccp(k,i,j)
   enddo
endif

! fill scratch arrays for rain

if (jnmb(2) >= 1) then
   do k = 2,m1-1
         rx(k,2) = micro_g(ng)%rrp(k,i,j)
         if (jnmb(2) >= 5) cx(k,2) = micro_g(ng)%crp(k,i,j)
   enddo
endif

! fill scratch arrays for pristine ice

if (jnmb(3) >= 1) then
   do k = 2,m1-1
         rx(k,3) = micro_g(ng)%rpp(k,i,j)
         if (jnmb(3) >= 5) cx(k,3) = micro_g(ng)%cpp(k,i,j)
   enddo
endif

! fill scratch arrays for snow

if (jnmb(4) >= 1) then
   do k = 2,m1-1
         rx(k,4) = micro_g(ng)%rsp(k,i,j)
         if (jnmb(4) >= 5) cx(k,4) = micro_g(ng)%csp(k,i,j)
   enddo
endif

! fill scratch arrays for aggregates

if (jnmb(5) >= 1) then
   do k = 2,m1-1
         rx(k,5) = micro_g(ng)%rap(k,i,j)
         if (jnmb(5) >= 5) cx(k,5) = micro_g(ng)%cap(k,i,j)
   enddo
endif

! fill scratch arrays for graupel

if (jnmb(6) >= 1) then
   do k = 2,m1-1
         rx(k,6) = micro_g(ng)%rgp(k,i,j)
         if (jnmb(6) >= 5) cx(k,6) = micro_g(ng)%cgp(k,i,j)
   enddo
endif

! fill scratch arrays for hail

if (jnmb(7) >= 1) then
   do k = 2,m1-1
         rx(k,7) = micro_g(ng)%rhp(k,i,j)
         if (jnmb(7) >= 5) cx(k,7) = micro_g(ng)%chp(k,i,j)
   enddo
endif

!Saleeby(2008): Add IF block here for drizzle (jnmb(8)>=1)

! Code from micphys: complete diagnosis of mean particle mass and number conc.

do lcat = 1,7 !Saleeby(2008): Increase to 8 for adding drizzle
   if (jnmb(lcat) >= 1) then
      CALL rad_enemb (m1,lcat,dn0(1))
   endif
enddo

return
END SUBROUTINE cloud_prep

!##############################################################################
Subroutine rad_enemb (m1,lcat,dn0)

use micphys

implicit none

integer :: m1,lcat,k,lhcat
real :: embi,parmi
real, dimension(m1) :: dn0

if (jnmb(lcat) == 2) then
   embi = 1. / emb(2,lcat)
   do k = 2,m1-1
      cx(k,lcat) = rx(k,lcat) * embi
   enddo
elseif (jnmb(lcat) == 3) then
   do k = 2,m1-1
      lhcat = jhcat(k,lcat)
      emb(k,lcat) = cfemb0(lhcat) * (dn0(k) * rx(k,lcat)) ** pwemb0(lhcat)
      cx(k,lcat) = cfen0(lhcat) / dn0(k)  &
         * (dn0(k) * rx(k,lcat)) ** pwen0(lhcat)
   enddo
elseif (jnmb(lcat) == 4) then
   parmi = 1. / parm(lcat)
   do k = 2,m1-1
      emb(k,lcat) = max(emb0(lcat),min(emb1(lcat),rx(k,lcat) * parmi))
      cx(k,lcat) = rx(k,lcat) / emb(k,lcat)
   enddo
elseif (jnmb(lcat) >= 5) then
   do k = 2,m1-1
      emb(k,lcat) = max(emb0(lcat),min(emb1(lcat),rx(k,lcat)  &
         / max(1.e-9,cx(k,lcat))))
      cx(k,lcat) = rx(k,lcat) / emb(k,lcat)
   enddo
endif

return
END SUBROUTINE rad_enemb

