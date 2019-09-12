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

Subroutine dust_sources (m1,m2,m3,ia,iz,ja,jz &
   ,lclass,rtgm,parea,vegrough,u,v  &
   ,md1np,md2np,md1mp,md2mp,swater,stext,glat,glon,dn0)

use micphys
use aero_include
use mem_grid
use mem_micro

implicit none

integer, parameter :: ncls=7
integer :: m1,m2,m3,ia,iz,ja,jz,i,j,z0limit,veglimit,dustloftsource
real, dimension(m2,m3) :: u10,v10,rtgm,num_dust_s,num_dust_l,glat,glon
real, dimension(m1,m2,m3) :: u,v,md1np,md2np,md1mp,md2mp,dn0
real, dimension(nzg,m2,m3,npatch) :: swater,stext
real, dimension(m2,m3,npatch) :: lclass,parea,vegrough
real :: dust1_size,dust2_size,aux1,aux2,aux3

 !******************************************************************
 !Dust lofting option flags
 !******************************************************************
 !Type of dust lofting
 dustloftsource=0   !0 = idealized(Fecan,Pierre), 1 = Ginoux(2001)
 !Limit dust flux due to vegetative roughness
 z0limit=1          !0 = off 1 = on (only for dustloftsource=0)
 !Limit dust lofting to specific bare or low grass surface types
 !3-Desert,bare soil, 8-Short-grass, 9-Tall-grass, 10-Semi-desert,
 !15-Crops, 18-Wooded-grassland
 veglimit=0         !0 = off 1 = on (only for dustloftsource=0)
 !******************************************************************

 !Checks:
 if(dustloftsource==1 .and. (z0limit==1 .or. veglimit==1)) then
   PRINT*,'SHOULD NOT RUN GINOUX(2001) DUST LOFTING WITH'
   PRINT*,'PIERRE(2012) VEGETATION ROUGHNESS LIMIT OR'
   PRINT*,'VEGETATION TYPE LIMITATION TO SHORT GRASS AND DESERT'
   PRINT*,'SEE FLAGS AT TOP OF ROUTINE DUST_SOURCES'
   STOP
 endif

 !Reads Ginoux's source only during the first call
 if(dustloftsource==1 .and. time.eq.0.0) then
   !Retrieve the dust source data from Ginoux's 1deg lat/lon data file
   open(150,file = '../etc/DustEmission.dat', status = 'old')
   do j = 1,ny_source
   do i = 1,nx_source
     read(150,*) aux1, aux2, aux3
     !Latitude and longitude locations of source data
     lon_source(i) = aux1
     lat_source(j) = aux2
     source(i,j)   = aux3
   enddo
   enddo
   close(150)
   !For visualization in Grads format
   !open(unit=21,file='source.gra',status='unknown', &
   !  form='unformatted',access='direct',recl=4*nx_source*ny_source)
   !write(21,rec=1) ((source(i,j),i=1,nx_source),j=1,ny_source)
   !close(21)
 endif

 !Get the 10m wind components
 CALL get_u10_v10 (m1,m2,m3,ia,iz,ja,jz,u,v,u10,v10,zm,rtgm)

 !Compute dust sources one gridcell at a time
 CALL dustflux (u10,v10,swater,stext,nzg,npatch,m1,m2,m3,ia,iz,ja,jz &
       ,lclass,parea,vegrough,zm,rtgm,dtlt,num_dust_s,num_dust_l,ncls &
       ,z0limit,veglimit,dustloftsource,glat,glon)

 !Set mean mass weighted radii of dust modes that corressponds to dust func.
 !Bin radii(um): r1=0.15,r2=0.265,r3=0.471,r4=0.838,r5=1.5,r6=2.65,r7=4.71
 !Apportion flux into size bins with sp, where sp are bin mass fractions.
 !For the 4 small clay size bins, each bin contains the following fractions:
 ! 1=0.9%, 2=8.1%, 3=23.4%, 4=67.6%
 !However, the 4 small clay bins represent only 1/10th the total dust silt,
 ! so the parameter "sp(1:4)" represents decimal-percentage/10
 !The three large classes are weight equally at 33.333% each.
 !sp(1:4)     = (/0.0009,0.0081,0.0234,0.0676/)
 !sp(5:ncls)  = 0.30
 !dust1_size = mass weighted mean fine dust radius (r1 - r4)
 !dust2_size = mass weighted mean course dust radius (r5 - r7)
 !fine: 0.15*0.009 + 0.265*0.081 + 0.471*0.234 + 0.838*0.676 = 0.699 um
 !coarse: 1.5*0.333 + 2.65*0.333 + 4.71*0.333 = 2.95038 microns
 dust1_size=0.699e-6 !small mode (meters)
 dust2_size=2.95e-6  !large mode (meters)

 !Update dust particles due to dust sources. Dust source model outputs in 
 !#/cm3 so we convert to #/kg and compute dust mass mixing ratio in kg/kg
 do j = ja,jz
  do i = ia,iz

    !Convert dust number (#/cm3) to (#/kg)
    num_dust_s(i,j) = num_dust_s(i,j) / dn0(2,i,j) * 1.e6
    num_dust_l(i,j) = num_dust_l(i,j) / dn0(2,i,j) * 1.e6

    !Update number of dust species
    md1np(2,i,j) = md1np(2,i,j) + num_dust_s(i,j)
    md1np(1,i,j) = md1np(2,i,j)

    md2np(2,i,j) = md2np(2,i,j) + num_dust_l(i,j)
    md2np(1,i,j) = md2np(2,i,j)

    !Update mass of dust species using mean mass radii rather than median radii
    md1mp(2,i,j) = md1mp(2,i,j) + &
        (dust1_size**3.)*num_dust_s(i,j)/(0.23873/aero_rhosol(3))
    md1mp(1,i,j) = md1mp(2,i,j)

    md2mp(2,i,j) = md2mp(2,i,j) + &
        (dust2_size**3.)*num_dust_l(i,j)/(0.23873/aero_rhosol(4))
    md2mp(1,i,j) = md2mp(2,i,j)

  enddo
 enddo

return
END SUBROUTINE dust_sources

!##############################################################################
Subroutine salt_sources (m1,m2,m3,ia,iz,ja,jz &
   ,lclass,rtgm,parea,u,v &
   ,salt_film_np,salt_jet_np,salt_spum_np &
   ,salt_film_mp,salt_jet_mp,salt_spum_mp,dn0)

use micphys
use aero_include
use mem_grid
use mem_micro

implicit none

!Assumed salt masses based on radii below and density of salt 2.165 g/cm3
!Assumes all aerosols in each category are all the same size
!Need to alter these if we assume a certain distribution
!Values are the mass of a single particle with radii "film_med" & "jet_med"
!For assumed masses and sizes of large and small salt and dust,
!see values set in aero_include.f90

!Salt mass and number arrays
real, dimension(m1,m2,m3) :: salt_film_np,salt_jet_np,salt_spum_np
real, dimension(m1,m2,m3) :: salt_film_mp,salt_jet_mp,salt_spum_mp

integer :: m1,m2,m3,ia,iz,ja,jz,i,j
real, dimension(m2,m3)    :: u10,v10,rtgm
real, dimension(m1,m2,m3) :: u,v,dn0
real, dimension(m2,m3)    :: film_source,jet_source,spm_source
real, dimension(m2,m3,npatch) :: lclass,parea
real :: sf_size,sj_size,ss_size

!Get the 10m wind components
 CALL get_u10_v10 (m1,m2,m3,ia,iz,ja,jz,u,v,u10,v10,zm,rtgm)

!Get Salt emission
 CALL salt_flux (u10,v10,m1,m2,m3,ia,iz,ja,jz,npatch,lclass,parea,dtlt &
               ,salt_film_np,salt_jet_np,salt_spum_np &
               ,film_source,jet_source,spm_source,dn0)

!Set median radii of salt modes that corresponds to salt func
sf_size=0.10e-6 !film mode radius (m3)
sj_size=1.00e-6 !jet mode radius (m3)
ss_size=6.00e-6 !spume mode radius (m3)

!Nudge diagnostic values over sea
do j = ja,jz
 do i = ia,iz

    !Convert dust number (#/m3) to (#/kg)
    film_source(i,j) = film_source(i,j) / dn0(2,i,j)
    jet_source(i,j)  = jet_source(i,j)  / dn0(2,i,j)
    spm_source(i,j)  = spm_source(i,j)  / dn0(2,i,j)

    !Update number of salt species
    salt_film_np(2,i,j) = salt_film_np(2,i,j) + film_source(i,j)
    salt_film_np(1,i,j) = salt_film_np(2,i,j)

    salt_jet_np(2,i,j)  = salt_jet_np(2,i,j)  + jet_source(i,j)
    salt_jet_np(1,i,j)  = salt_jet_np(2,i,j)

    salt_spum_np(2,i,j) = salt_spum_np(2,i,j) + spm_source(i,j)
    salt_spum_np(1,i,j) = salt_spum_np(2,i,j)
    
    !Update mass of salt species
    salt_film_mp(2,i,j) = salt_film_mp(2,i,j) + &
        ((sf_size*aero_rg2rm(5))**3.)*film_source(i,j)/(0.23873/aero_rhosol(5))
    salt_film_mp(1,i,j) = salt_film_mp(2,i,j)

    salt_jet_mp(2,i,j) = salt_jet_mp(2,i,j) + &
        ((sj_size*aero_rg2rm(6))**3.)*jet_source(i,j) /(0.23873/aero_rhosol(6))
    salt_jet_mp(1,i,j) = salt_jet_mp(2,i,j)

    salt_spum_mp(2,i,j) = salt_spum_mp(2,i,j) + &
        ((ss_size*aero_rg2rm(7))**3.)*spm_source(i,j) /(0.23873/aero_rhosol(7))
    salt_spum_mp(1,i,j) = salt_spum_mp(2,i,j)

 enddo
enddo

return
END SUBROUTINE salt_sources

!##############################################################################
Subroutine get_u10_v10 (m1,m2,m3,ia,iz,ja,jz,uu,vv,u10,v10,zm,rtgm)

implicit none

integer :: m1,m2,m3,i,j,ia,iz,ja,jz
real, parameter :: z0 = 0.05
real :: dz,aux1,aux2
real, dimension(m1) :: zm
real, dimension(m2,m3) :: rtgm
real, dimension(m2,m3) :: u10,v10
real, dimension(m1,m2,m3) :: uu,vv

aux2 = log(10./z0)

!Note that zm(1) is always 0.0m above ground and zm(2) is height
! the first atmospheric level on the "m" grid.
!This is definitely NOT valid for anywhere but bare soil surfaces
! Fine for now as the source func is only valid for bare sfcs
! but if the source func is changed in the future, the roughness
! length from LEAF2 should be used instead of an assumed constant.
do j=ja,jz
do i=ia,iz
  dz = zm(2) * rtgm(i,j)
  aux1 = 1./log(dz/z0)
  u10(i,j) = aux1 * uu(1,i,j) * aux2
  v10(i,j) = aux1 * vv(1,i,j) * aux2
enddo
enddo

return
END SUBROUTINE get_u10_v10

!##############################################################################
Subroutine salt_flux (u10,v10,m1,m2,m3,ia,iz,ja,jz,npatch,lclass,parea &
 ,dtlt,salt_film_np,salt_jet_np,salt_spum_np,film_source,jet_source   &
 ,spm_source,dn0)

!Code to calculate seasalt number flux on global scale based on the O'dowd
!1997, 1999 where u10 is the 10-meter wind speed m/s, salt_film is the number
!of particles /m3 in the submicron mode (ccn), and salt_jet is the number of
!particles /m3 in the supermicron mode (gccn), and salt_spume is the number
!particles /m3 in the ultra-giant mode (ultra-gccn), and film and jet medians
!are user defined in centimeters.
!The total number concentrations for each mode were observed to be:
! log Nfilm = 0.095 U10 + 6.2830, 0.1 microns mode radius ----> CCN
! log Njet = 0.0422 U10 + 5.7122, 1 micron mode radius -------> GCCN
! log Nspume= 0.069 U10 + 0.1900, 6 micron mode radius -------> Super giant
! O’Dowd C. D., Smith M. H. and Jennings S. G. (1993) Submicron aerosol, radon and
! soot carbon characteristics over the North East Atlantic. J. geophys. Res. 98,
! 1132-1136. Updated for O'Dowd 1997,1999.

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,i,j,npatch,np
real :: dtlt,time_salt
real, dimension(m1,m2,m3) :: salt_film_np,salt_jet_np,salt_spum_np,dn0
real, dimension(m2,m3) :: salt_film_dia,salt_jet_dia,salt_spm_dia
real, dimension(m2,m3) :: wind,u10,v10,film_source,jet_source,spm_source
double precision, dimension(m2,m3)::salt_film_dia2,salt_jet_dia2 &
  ,salt_spm_dia2,wind_2
real, dimension(m2,m3,npatch) :: lclass,parea

time_salt = dtlt !relaxation time in seconds

do j = ja,jz
do i = ia,iz

 !Zero out the final source func prior to updated timestep computation
 film_source(i,j) = 0.0
 jet_source(i,j)  = 0.0
 spm_source(i,j)  = 0.0

 !Compute salt source term
 !Limit wind speed to 30m/s (Fan and Toon 2011)
 do np = 1,npatch
  if(nint(lclass(i,j,np)) .eq. 0 .and. parea(i,j,np) .gt. 0.009) then
    !Compute wind speed from 10m u,v winds in [m/s]
    wind(i,j) = min(30.,sqrt(u10(i,j)*u10(i,j) + v10(i,j)*v10(i,j)))
    wind_2(i,j) = dble(wind(i,j))
    !Compute predicted salt concentration based on wind speed (/m3)
    salt_film_dia2(i,j) = 10.**(0.0950*wind_2(i,j) + 6.2830)
    salt_film_dia(i,j)  = sngl(salt_film_dia2(i,j))
    salt_jet_dia2(i,j)  = 10.**(0.0422*wind_2(i,j) + 5.7122)
    salt_jet_dia(i,j)   = sngl(salt_jet_dia2(i,j))
    salt_spm_dia2(i,j)  = 10.**(0.0690*wind_2(i,j) + 0.1900)
    salt_spm_dia(i,j)   = sngl(salt_spm_dia2(i,j))
    !Note that sea salt predicted (/m3) while we carry (/kg) in the model.
    !So do correct density conversion here for adding particles.
    if(salt_film_dia(i,j) > salt_film_np(2,i,j)*dn0(2,i,j)) &
     film_source(i,j)=film_source(i,j) + (salt_film_dia(i,j) &
        - salt_film_np(2,i,j)*dn0(2,i,j)) * dtlt / time_salt * parea(i,j,np)
    if(salt_jet_dia(i,j)  > salt_jet_np(2,i,j)*dn0(2,i,j)) &
     jet_source(i,j) =jet_source(i,j)  + (salt_jet_dia(i,j)  &
        - salt_jet_np(2,i,j)*dn0(2,i,j))  * dtlt / time_salt * parea(i,j,np)
    if(salt_spm_dia(i,j)  > salt_spum_np(2,i,j)*dn0(2,i,j)) &
     spm_source(i,j) =spm_source(i,j)  + (salt_spm_dia(i,j)  &
        - salt_spum_np(2,i,j)*dn0(2,i,j)) * dtlt / time_salt * parea(i,j,np)
  endif
 enddo

enddo
enddo

return
END SUBROUTINE salt_flux

!##############################################################################
Subroutine dustflux (u10,v10,swater,stext,nzg,npatch,m1,m2,m3,ia,iz,ja,jz &
         ,xlclass,xparea,xvegrough,zm,rtgm,dtlt,num_dust_s,num_dust_l,ncls &
         ,z0limit,veglimit,dustloftsource,glat,glon)

!Code to calculate dust flux on regional-global scales based on dust code
!from Peter Colarco at Goddard. It calculates fluxes using a 1x1 degree
!horizontal resolution source func developed by Paul Ginoux. Source
!func is equivalent to the fraction of the grid cell emitting dust.
!Fluxes are calculated based on formulas given in Marticorena & Bermagetti,
!1995 using 10-m winds and soil moisture. Updated for Fecan et al.(1999)
!and Pierre et al.(2012).

use aero_include

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,ncls,nzg,npatch,icls,i,j,np,itype,sclass
integer :: z0limit,veglimit,vegstop,dustloftsource,isource,jsource

!Set number of bins, air density (g/cm3), gravity (cm/s2)
!and "fact" as parameters. "fact" scales mass emissions to 1851 tg
real :: rhoa,grav,fact,pi,S
parameter(rhoa=1.25d-3, grav=981., fact=1.153d-16, pi=3.141593)

real :: dz,dtlt,rmrat,rbmin,rho_p,uth75,rmassmin,Ez0,gwet
real, dimension(ncls) :: sp,den,uth,diam,rmass,r,utwet &
                        ,mass_flux,mass_conc,num_conc
real, dimension(m1) :: zm
real, dimension(m2,m3) :: num_dust_s,num_dust_l,wind,u10,v10,rtgm,glat,glon
real, dimension(nzg,m2,m3,npatch) :: swater,stext
real, dimension(m2,m3,npatch) :: xlclass,xparea,xvegrough
real, dimension(12) :: clayfrac,slmsts,wprime

!Define clay percentage for soil textural classes
!  1 sand            5 loam                 9 sandy clay
!  2 loamy sand      6 sandy clay loam      10 silty clay
!  3 sandy loam      7 silty clay loam      11 clay
!  4 silt loam       8 clay loam            12 peat
data clayfrac/0.0,5.0,10.0,12.0,18.0,28.0,33.0,33.0,42.0,48.0,70.0,0.0/
data slmsts/.395,.410,.435,.485,.451,.420,.477,.476,.426,.492,.482,.863/
!Using Fecan et al.(1999) paramerization based on clay percentage
!w'=0.0014(clayfrac)^2 + 0.17(clayfrac)
data wprime/0.00000,0.885000,1.84000,2.24160,3.51360,5.85760,7.13460 &
           ,7.13460,9.60960,11.3856,18.7600,0.00000/

rmrat = (100**3)**(1./8.) !rmrat=volume ratio between bins
rbmin  = 1.e-5*((1.+rmrat)/2.)**(1./3.) !radius of smallest bin [cm]
rho_p = 2.65 !particle density [g/cm3]

!Set up the bins/apportionment (Ginoux et al.2001)
!r1=0.15, r2=0.265, r3=0.471, r4=0.838, r5=1.5, r6=2.65, r7=4.71 microns
rmassmin = 4./3.*pi*rho_p*rbmin**3.
do icls = 1,ncls
  rmass(icls) = rmassmin*rmrat**(icls-1)
  r(icls) = (rmass(icls)/rho_p/(4./3.*pi))**(1./3.) !units of cm
enddo

!Apportion flux into size bins with sp, where sp are bin mass fractions
!For the 4 small clay size bins, each bin contains the following fractions:
! 1=0.9%, 2=8.1%, 3=23.4%, 4=67.6%
!However, the 4 small clay bins represent only 1/10th the total dust silt,
! so the parameter "sp(1:4)" represents decimal-percentage/10
!The three large classes are weight equally at 30%.
sp(1:4)     = (/0.0009,0.0081,0.0234,0.0676/)
sp(5:ncls)  = 0.30

diam        = 2*r
den(1:4)    = 2.50 !mass density of clay
den(5:ncls) = 2.65 !mass density of silt

!Threshold velocity calcs [cm/s] from
!Marticorena & Bermagetti(1995) Equations 3-6
do icls = 5,ncls
  uth(icls) = 0.13*sqrt(den(icls)*grav*diam(icls)/rhoa) &
              *sqrt(1. + 0.006/(den(icls)*grav*diam(icls)**2.5)) &
              /sqrt(1.928*(1331.*diam(icls)**1.56 + 0.38)**0.092 - 1.)
enddo

!For sub-micron particles assume uth is that of an reff=0.75 micron particle
!Gives higher uth for smaller particle sizes
uth75 = 0.13*sqrt(2.5*grav*(1.5e-4)/rhoa) &
        *sqrt(1. + 0.006/(2.5*grav*(1.5e-4)**2.5)) &
        /sqrt(1.928*(1331.*(1.5e-4)**1.56 + 0.38)**0.092 - 1.)
uth(1:4) = uth75

!Calculate the flux, mass and # concentration in the domain
!Compute fluxes and output as number distributions with median radii (cm)
do j = ja,jz
do i = ia,iz

 !Zero out the final source func prior to updated timestep computation
 num_dust_s(i,j) = 0.0
 num_dust_l(i,j) = 0.0

 !Look for the source to be used in gridcell (i,j)
 !Source data is in 1 degree increment
 !I-points are incremented westward from -179.5 longitude
 !J-points are incremented northward from -89.5 latitude
 if(dustloftsource==1) then
   isource = nint(glon(i,j)-(-179.50))+1
   jsource = nint(glat(i,j)-(-089.50))+1
   S = source(isource,jsource)
 else
   S = 1.0
 endif

 do np = 2,npatch

  if(xparea(i,j,np) .gt. 0.009) then
   itype = nint(xlclass(i,j,np))
   sclass = floor(stext(nzg,i,j,np))

   !If using Pierre et al.(2012) surface zoughness limit on dust lofting
   !if patch roughness > 3.10e-3 cm
   if(z0limit==1.and.xvegrough(i,j,np)*100. > 3.10e-3) then
      Ez0 = 0.7304 - (0.0804 * log10(xvegrough(i,j,np)*100.))
   else
      Ez0 = 1.0
   endif

   !If limiting dust lofting to only low or limited vegetation patches
   if(veglimit==1 .and. itype.ne.3  .and. itype.ne.8  .and. itype.ne.9 .and. & 
       itype.ne.10 .and. itype.ne.15 .and. itype.ne.18) then
     vegstop=1
   else
     vegstop=0
   endif

   !No dust lofting at all for water, ice, marsh, or urban
   if(itype==0.or.itype==1.or.itype==2.or.itype==17.or. &
      itype==19.or.itype==21) vegstop=1

   !Proceed to dust lofting if vegetation type allowed
   if(vegstop==0) then
    !Seigel (10-8-2010): Threshold velocity adjustment, accounting for soil 
    !moisture AND soil textural class. Note that swater is volumetric soil 
    !moisture m3/m3 and not total saturation. Ginoux et al.(2001) uses 
    !total saturation with values from 0.001 to 1.0. But Fecan et al.(1999)
    !uses volumetric soil moisture (m3/m3). So be careful switching between 
    !the two. For RAMS, the maximum volumetric soil moisture is slmsts(sclass).

    !Set soil moisture
    if(dustloftsource==1) then
      !If using Ginoux et al.(2001)
      gwet = max(0.0,min(1.0,swater(nzg,i,j,np)/slmsts(sclass)))
    else
      !If using Fecan et al. (1999)
      gwet = swater(nzg,i,j,np)
    endif

    do icls = 1,ncls
     if(dustloftsource==1) then
       !If using Ginoux et al.(2001) 
       if(gwet .lt. 0.5) then
          utwet(icls) = uth(icls)*(1.2 + 0.2 * alog10(gwet))
       else
          utwet(icls) = 1.e6 !essentially infinity for no lofting
       endif
     else
       !If using Fecan et al. (1999)
       if(gwet*100. .lt. wprime(sclass)) then
          utwet(icls) = uth(icls)
       else
          utwet(icls) = uth(icls)*sqrt(1 + 1.21* &
                           (gwet*100.-wprime(sclass))**0.68)
       endif
     endif
    enddo

    !Compute wind speed from 10m u,v winds in [cm/s]
    wind(i,j) = sqrt(u10(i,j)*u10(i,j) + v10(i,j)*v10(i,j))*100.

    !Note that zm(2) is the top of the first grid cell an zm(1) is surface
    dz = zm(2) * rtgm(i,j)

    !Flux dust if soil is dry
    if(gwet.lt.0.5) then
      do icls = 1,ncls
        !Calculate the mass flux in the domain [g cm-2 s-1]
        !Flux equation from Ginoux et al. 2001 equation 2.
        mass_flux(icls)=Ez0*S*fact*sp(icls)*(wind(i,j)**2.)*(wind(i,j)-utwet(icls))
        if(mass_flux(icls) <= 0.) then
          mass_flux(icls) = 0.
          mass_conc(icls) = 0.
          num_conc(icls)  = 0.
        endif
        if(mass_flux(icls) > 0.) then
          !Convert mass flux to mass concentration [g cm-3]
          !For mass concentration divide by dz in units of [cm]
          mass_conc(icls) = mass_flux(icls) * dtlt / (dz*100.)
          !Then convert to # concentration (particles/cm3) in each bin
          num_conc(icls)  = mass_conc(icls) / &
                               (den(icls)*(4./3.)*pi*r(icls)**3)
        endif
      enddo
      !Sum the sub and super micron bins for num_dust_s & num_dust_l
      do icls = 1,4
         num_dust_s(i,j) = num_dust_s(i,j) + num_conc(icls) * xparea(i,j,np)
      enddo
      do icls = 5,ncls
         num_dust_l(i,j) = num_dust_l(i,j) + num_conc(icls) * xparea(i,j,np)
      enddo
    endif !IF SOURCE EXISTS AND SOIL IS DRY 

   endif !IF patch is capable of releasing dust
  endif !IF patch area is big enough
 enddo !LOOP over surface patches

enddo !LOOP OVER J POINTS
enddo !LOOP OVER I POINTS

return
END SUBROUTINE dustflux
