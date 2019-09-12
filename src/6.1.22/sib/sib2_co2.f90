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

!SiB VERSION 2.5/3, MODIFIED FOR USE WITH RAMS

!Note: "prog" below indicates a SiB prognostic variable
!      "diag" below indicated a SiB diagnostic variable
!      "output" below indicates an output variable passed back to RAMS
!         for use in the main model carbon and surface fluxes.
!      "atmos" indicates fields given by LEAF and/or RAMS for use
!         in driving SiB.
!      "BCs" below indicates soil and vegetation conditions set by the
!         user or MAKESFC surface file conditions from surface characteristics
!         data files.

Subroutine sib_2pt5 ( &
  !THIS SECTION: INPUT GRID AND ATMOSPHERIC CONDITIONS FROM RAMS TO SIB
   igrp        & !I grid point
  ,jgrp        & !J grid point
  ,mzg         & !number of soil levels
  ,mzs         & !number of snow levels
  ,dt          & !timestep
  ,doy         & !time in julian date (Day Of Year)
  ,latitude    & !latitude
  ,thm         & !atmos: -lowest layer theta (potential temperature K)
  ,sh_x        & !atmos: -lowest layer rv (water vapor mixing ratio, g/kg)
  ,ps_x        & !atmos: -lowest layer air pressure (hPa=mb) RAMS uses (Pa)
  ,ros         & !atmos: -lowest layer base state density (kg/m3)
  ,ts_x        & !atmos: -lowest layer temperature (K)
  ,spdm        & !atmos: -lowest layer wind speed (m/s)
  ,zwind       & !atmos: -lowest model level height t-grid (m)
  ,cupr        & !atmos: -cuparm precip rate (conprr) (mm/sec)
  ,lspr        & !atmos: -microphysics precip rate (pcpg / dt) (mm/sec)
  ,dswbot      & !atmos: -surface incident shortwave radiation (W/m^2)
  ,dlwbot      & !atmos: -surface incident longwave radiation (W/m^2)
  ,cosz        & !atmos: -cosine solar zenith angle
  !THIS SECTION: SURFACE FIELD INFO FROM RAMS TO SIB
  ,patcharea   & !BCs: -fractional patch area
  ,biome_f     & !BCs: -leaf vegetation class
  ,cndvi       & !BCs: -current month ndvi
  ,pndvi       & !BCs: -previous/past month ndvi
  ,soiltype_f  & !BSs: -soil textural class
  !THIS SECTION: PROGNOSTIC VALUES FROM SIB TO FEEDBACK TO RAMS
  ,rst         & !prog: -CANOPY/STOMATAL RESISTANCE (S M-1)
  ,ta          & !prog: -canopy (CAS) temperature (K)
  ,sha         & !prog: -canopy (CAS) water vapor mixing ratio (kg/kg)
  ,tc          & !prog: -vegetation temperature (K)
  ,tempk       & !prog: -soil and snow temperature (K)
  ,soil_water  & !prog: -volumetric soil moisture (m3/m3)
  ,sfcswa      & !prog: -surface abledo (fraction)
  ,uplwrf      & !prog: -upwelling longwave radiation (W/m2)
  ,snow1       & !prog: -vegetation snow (kg/m2)
  ,snow2       & !prog: -ground surface snow (kg/m2)
  ,capac1      & !prog: -vegetation liquid store (kg/m^2)
  ,capac2      & !prog: -ground surface liquid interception store (kg/m^2)
  ,pco2ap      & !prog: -canopy air space pCO2 concentration (Pa)
  ,co2flx      & !prog: -sfc CO2 flux between CAS and ref lev(mol/m2/sec)
  ,ustar       & !prog: -u-star
  ,rstar       & !prog: -r-star
  ,tstar       & !prog: -t-star
  !THIS SECTION: INPUT CURRENT CO2 CONENTRATION TO SIB
  ,atm_co2_ppm & !lowest atm level CO2, in ppm
  !THIS SECTION: DIAGNOSTIC QUANTITIES FROM SIB TO FILL STANDARD LEAF3 VARS
  ,sfcwaterlev & !diag: number of sfc water levels (0 or 1 from SiB)
  ,sfcwatermas & !diag: sfc water mass (snow+liquid) (kg/m2)
  ,sfcwaterdep & !diag: sfc water depth (snow+liquid) (m)
  ,vegwatermas & !diag: vegetation water mass (snow+liquid) (kg/m2)
  !THIS SECTION: DIAGNOSTIC QUANTITIES UNIQUE TO SIB
  ,vegalb_out  & !diag: vegetation albedo (fraction)
  ,vcover_out  & !diag: vegetation fractional area
  ,vlt_out     & !diag: vegetation LAI
  ,zlt_out     & !diag: total LAI
  ,veght_out   & !diag: canopy top height (m)
  ,vrough_out  & !diag: vegetation roughness (m)
  ,assimn_out  & !diag: net co2 assim by plants (umol/m^2/sec)
  ,respg_out   & !diag: ground respiration flux (umol/m^2/sec)
  ,rstfac1_out & !diag: CANOPY RESISTANCE STRESS1 : leaf sfc humidity
  ,rstfac2_out & !diag: CANOPY RESISTANCE STRESS2 : soil moisture
  ,rstfac3_out & !diag: CANOPY RESISTANCE STRESS3 : temperature
  ,ect_out     & !diag: transpiration flux (W/m^2)
  ,eci_out     & !diag: canopy interception flux (W/m^2)
  ,egi_out     & !diag: ground interception flux (W/m^2)
  ,egs_out     & !diag: ground surface layer evap (W/m^2)
  ,hc_out      & !diag: canopy (veg) sensible heat flux (W/m^2)
  ,hg_out      & !diag: ground surface sensible heat flux (W/m^2)
  ,ra_out      & !diag: CAS-RAMS aerodynamic resistance (sec/m)
  ,rb_out      & !diag: leaf sfc-CAS resistance (sec/m)
  ,rc_out      & !diag: total canopy resistance (sec/m)
  ,rd_out      & !diag: ground-CAS resistance (sec/m)
  ,roff_out    & !diag: runoff (surface and subsurface) (mm)
  ,green_out   & !diag: greenness fraction (-)
  ,apar_out    & !diag: absorbed fraction of PAR
  ,ventmf_out  & !diag: ventilation mass flux (kg/m^2/sec)
  ,pco2c_out   & !diag: leaf chloroplast CO2 concentration (Pa)
  ,pco2i_out   & !diag: leaf internal CO2 concentration (Pa)
  ,pco2s_out   & !diag: leaf surface CO2 concentration (Pa)
  ,pco2m_out   & !diag: lowest atm level CO2 concentration (Pa)
  ,ea_out      & !diag: CAS water vapor pressure (hPa)
  ,em_out      & !diag: ref level vapor pressure (hPa)
  ,rha_out     & !diag: CAS relative humidity
  ,radvbc_out  & !diag: radiation: visible beam (W/m^2)
  ,radvdc_out  & !diag: radiation: visible diffuse (W/m^2)
  ,radnbc_out  & !diag: radiation: nir beam (W/m^2)
  ,radndc_out  & !diag: radiation: nir diffuse (W/m^2)
  ,psy_out     ) !diag: psychrometric constant (hPa deg^-1)

use mem_grid, only:ngrid
use mem_sib
use node_mod, only:mi0,mj0

implicit none

  !Variables with the '_x' subscript were modified to avoid
  !conflicts with RAMS common block variables.

  !-------------------------------------------------------------------
  ! REFERENCES: Sato, N., P. J. Sellers, D. A. Randall, E. K. Schneider,
  !     J. Shukla, J. L Kinter III, Y-T, Hou, and Albertazzi (1989)
  !     "Effects of implementing the simple biosphere model in a general
  !     circulation model. J. Atmos. Sci., 46, 2767-2782.
  !            Sellers, P. J., D. A. Randall, C. J. Collatz, J. A. Berry,
  !     C. B. Field, D. A. Dazlich, C. Zhang, G. Collelo (1996) A revise
  !     land-surface parameterization (SiB2) for atmospheric GCMs. Part 1:
  !     Model formulation. (accepted by JCL)
  !  MODIFICATIONS:
  !   - changed VQSAT call to VNQSAT.  kwitt 10/23
  !   - added in the prognostic stomatal conductance in addinc. changan
  !   - moved sib diagnostics accumulation from dcontrol to bldif
  !     dd 950202
  !  ROUTINES called:  VNQSAT, SNOW1, balan, VNTLAT
  !       DELHF, DELEF, NETRAD, SIBSLV, endtem, updat2, addinc
  !       inter2, balan, soilprop, soiltherm, begtem, rnload
  !  FUNCS called:
  !       none
  !-------------------------------------------------------------------

  !FIXING LEN and NSIB AT 1-WILL RUN AS A SINGLE POINT.
  INTEGER, PARAMETER :: len = 1  ! Run as single point with RAMS interface
  INTEGER, PARAMETER :: nsib = 1 ! Run as single point with RAMS interface
  INTEGER, PARAMETER :: nsoil = nzg_sib - 1  ! number of soil layers
  INTEGER, PARAMETER :: ioffset = 0 ! subdomain offset
  INTEGER :: sibprint,i,j,k,n,ksoil,l

  !PI
  REAL, PARAMETER :: num_pi = 3.1415926
  !gravity (m/s2)
  REAL, PARAMETER :: grav = 9.81
  !specific heat of air at const pres (J kg-1 deg-1)
  REAL, PARAMETER :: cp = 1004.
  !specific heat of air at const volume (J kg-1 deg-1)
  REAL, PARAMETER :: cv = 1952.
  !universal gas constant
  REAL, PARAMETER :: rgas = 287.
  REAL, PARAMETER :: hltm = 2.52E6
  REAL, PARAMETER :: delta = 0.608
  !conversion for kg water to snow depth (16.7)
  REAL, PARAMETER :: asnow = 16.7
  !R/cp
  REAL, PARAMETER :: kapa = 0.2861328125
  !latent heat of fusion for ice (J m^-3)
  REAL, PARAMETER :: snomel = 3.705185e8
  REAL, PARAMETER :: clai = 4.186*1000.0*0.2
  REAL, PARAMETER :: cww = 4.186*1000.0*1000.
  !von karman constant
  REAL, PARAMETER :: vkrmn = 0.35
  REAL, PARAMETER :: po2m = 20900.
  !stefan-boltzmann constant
  REAL, PARAMETER :: stefan = 5.67e-8
  REAL, PARAMETER :: grav2 = grav *0.01
  REAL, PARAMETER :: tice = 273.16
  REAL, PARAMETER :: snofac = hltm / ( hltm + SNOMEL * 1.E-3 )

  !VARIABLES INPUT/OUTPUT FROM/TO RAMS MODEL
  INTEGER :: mzs,mzg,igrp,jgrp,doy
  REAL, DIMENSION(mzg+mzs) :: tempk
  REAL, DIMENSION(mzg) :: soil_water
  REAL ::             &
        dt            &! time step (s)
       ,latitude      &! latitude of current grid cell
       ,thm(len)      &! mixed layer potential temperature (K)
       ,sh_x(len)     &! mixed layer water vapor mixing ratio (kg/kg)
       ,ps_x(len)     &! surface pressure (hPa=mb) RAMS uses (Pa)
       ,ros (len)     &! surface air density (kg/m^3)
       ,ts_x(len)     &! surface mixed layer air temperature (K)
       ,spdm(len)     &! boundary layer wind speed (m/s)
       ,zwind         &! lowest model level height t-grid (m)
       ,cupr(len)     &! convective-parm precipitation rate (mm/s)
       ,lspr(len)     &! microphysics precipitation rate (mm/s)
       ,dlwbot(len)   &! surface incident longwave radiation (W/m^2)
       ,dswbot(len)   &! surface incident shortwave radiation (W/m^2)
       ,cosz(len)     &! cosine of solar zenith angle
       ,patcharea     &! land patch fractional area
       ,biome_f       &! leaf vegetation class
       ,pndvi         &! past    value of ndvi for the gridcell
       ,cndvi         &! current value of ndvi for the gridcell
       ,soiltype_f    &! soil textural class
       ,rst(len)      &! canopy/stomatal resistance (S M-1)
       ,ta(len)       &! CAS temperature (K)
       ,sha(len)      &! CAS water vapor mixing ratio (kg/kg)
       ,tc(len)       &! canopy (vegetation) temperature (K)
       ,sfcswa        &! surface abledo (fraction)
       ,uplwrf        &! upwelling longwave radiation (W/m2)
       ,snow1(len)    &! vegetation snow (kg/m2)
       ,snow2(len)    &! ground surface snow (kg/m2)
       ,capac1(len)   &! vegetation liquid store (kg/m^2)
       ,capac2(len)   &! ground surface liquid interception store (kg/m^2)
       ,pco2ap(len)   &! canopy air space pCO2 (Pa)
       ,co2flx(len)   &! sfc CO2 flux between CAS and ref lev(mol/m^2/sec)
       ,ustar(len)    &! friction velocity (m/s)
       ,rstar(len)    &! moisture exchange parameter based on U* (kg/kg)/(m/s)
       ,tstar(len)    &! heat exchange parameter based on U* (K)/(m/s)
       ,atm_co2_ppm   &! lowest atm level CO2, in ppm
       ,sfcwaterlev   &! number of surface water layers (0 or 1)
       ,sfcwatermas   &! sfc water mass (snow+liquid) (kg/m2)
       ,sfcwaterdep   &! sfc water depth (snow+liquid) (m)
       ,vegwatermas   &! vegetation water mass (snow+liquid) (kg/m2)
       ,vegalb_out    &!diag: vegetation albedo (fraction)
       ,vcover_out    &!diag: vegetation fractional area
       ,vlt_out       &!diag: vegetation LAI
       ,zlt_out       &!diag: total LAI
       ,veght_out     &!diag: canopy top height (m)
       ,vrough_out    &!diag: vegetation roughness (m)
       ,assimn_out    &!diag: net co2 assim by plants (umol/m^2/sec)
       ,respg_out     &!diag: ground respiration flux (umol/m^2/sec)
       ,rstfac1_out   &!diag: CANOPY RESISTANCE STRESS1 :leaf sfc humidity
       ,rstfac2_out   &!diag: CANOPY RESISTANCE STRESS2 :soil moisture
       ,rstfac3_out   &!diag: CANOPY RESISTANCE STRESS3 :temperature
       ,ect_out       &!diag: transpiration flux (W/m^2)
       ,eci_out       &!diag: canopy interception flux (W/m^2)
       ,egi_out       &!diag: ground interception flux (W/m^2)
       ,egs_out       &!diag: ground surface layer evap (W/m^2)
       ,hc_out        &!diag: canopy (veg) sensible heat flux (W/m^2)
       ,hg_out        &!diag: ground surface sensible heat flux (W/m^2)
       ,ra_out        &!diag: CAS-RAMS resistance (sec/m)
       ,rb_out        &!diag: leaf sfc-CAS resistance (sec/m)
       ,rc_out        &!diag: total canopy resistance (sec/m)
       ,rd_out        &!diag: ground-CAS resistance (sec/m)
       ,roff_out      &!diag: runoff (surface and subsurface) (mm)
       ,green_out     &!diag: greenness fraction (-)
       ,apar_out      &!diag: absorbed fraction of PAR
       ,ventmf_out    &!diag: ventilation mass flux (kg/m^2/sec)
       ,pco2c_out     &!diag: leaf chloroplast CO2 concentration (Pa)
       ,pco2i_out     &!diag: leaf internal CO2 concentration (Pa)
       ,pco2s_out     &!diag: leaf surface CO2 concentration (Pa)
       ,pco2m_out     &!diag: lowest atm level CO2 concentration (Pa)
       ,ea_out        &!diag: CAS water vapor pressure (hPa)
       ,em_out        &!diag: ref level vapor pressure (hPa)
       ,rha_out       &!diag: CAS relative humidity
       ,radvbc_out    &!diag: radiation: visible beam (W/m^2)
       ,radvdc_out    &!diag: radiation: visible diffuse (W/m^2)
       ,radnbc_out    &!diag: radiation: nir beam (W/m^2)
       ,radndc_out    &!diag: radiation: nir diffuse (W/m^2)
       ,psy_out        !diag: psychrometric constant (hPa deg^-1)

  !SIB STATIC SURFACE PARAMETERS
  REAL ::             &
        z0d(len)      &! surface roughness length (m)
       ,z0(len)       &! sfc rough length corrected for canopy snow (m)
       ,zlt(len)      &! leaf area index
       ,z1(len)       &! canopy bottom height (meters)
       ,z2(len)       &! canopy top height (meters)
       ,cc1(len)      &! RB Coefficient (c1) = rbc
       ,cc2(len)      &! RC Coefficient (c2) = rdc
       ,dd(len)       &! Zero plane displacement
       ,poros(len)    &! soil porosity
       ,zdepth(len,3) &! porosity * soil hydrology model layer depths (m)
       ,phsat(len)    &! Soil tension at saturation (units)
       ,bee(len)      &! Clapp & Hornberge 'B' exponent
       ,respcp(len)   &! respiration fraction of Vmax
       ,vmax0(len)    &! rubisco velocity of sun leaf (mol/m2/s)
       ,green(len)    &! Canopy greeness fraction of LAI
       ,tran(len,2,2) &! Leaf transmittance
       ,ref(len,2,2)  &! Leaf reflectance
       ,gmudmu(len)   &! Time-mean leaf projection (leaf orientation to par flux)
       ,trop(len)     &! temperature coefficient in GS-A model (K)
       ,phc(len)      &! one-half critical leaf-water potential limit(m)
       ,trda(len)     &! slope of high temp inhibition (leaf resp,1/K)
       ,trdm(len)     &! half point of high temp inhibition (leaf resp,K)
       ,slti(len)     &! slope of low temperature inhibition (1/K)
       ,shti(len)     &! slope of high temperature inhibition (1/K)
       ,hlti(len)     &! half piont of low temp inhibition (K)
       ,hhti(len)     &! half point of high temp inhibition (K)
       ,effcon(len)   &! quantum efficiency (mol/mol)
       ,binter(len)   &! conductance-photosynthesis intercept (mol/m2/s)
       ,gradm(len)    &! conductance-photosynthesis slope parm (mol/m2/s)
       ,atheta(len)   &! wc,we coupling parameter
       ,btheta(len)   &! wp,ws coupling parameter
       ,aparc(len)    &! Canopy absorbed fraction of PAR
       ,wopt(len)     &! Factor coeff for moisture effect on soil respiration
       ,zmx(len)      &! Power coeff for moisture effect on soil respiration
       ,wsat(len)     &! respiration at soil water saturation?
       ,vcover(len)   &! vegetation cover fraction
       ,sodep(len)    &! total soil depth (meters)
       ,rootd(len)    &! rooting depth (meters)
       ,soref(len,2)  &! Soil reflectance
       ,thermk(len)   &! canopy gap fraction for TIR radiation
       ,satco(len)    &! soil tension at 1/2 assimilation value (true?) units?
       ,slope(len)    &! slope
       ,chil(len)     &! leaf angle distribution factor
       ,ztdep(len,nsoil) ! soil thermal model layer depths (m)

  !    Intent: in/out and some variable local copies
  REAL ::                 & 
        cas_cap_heat(len) &! CAS heat capacity (J/K m^2)
       ,cas_cap_vap(len)  &! CAS vapor capacity
       ,cas_cap_co2(len)  &! CAS CO2 capacity (m/m^2) (moles air / m^2 in phosib)
       ,tg(len)           &! surface boundary temperature (K)
       ,td(len,nsoil)     &! deep soil temperature (K)
       ,www(len,3)        &! soil wetness
       ,wwwtem(len,3)     &! soil wetness copy
       ,snow(len,2)       &! snow cover (kg/m2)
       ,capac(len,2)      &! liquid interception store (kg/m^2)
       ,cuprt(len)        &! copy of cupr
       ,lsprt(len)        &! copy of lspr
       ,thmtem(len)       &! copy of thm
       ,shtem(len)        &! copy of sh
       ,zzwind(len)       &! intermediate surface wind for roughness length
       ,zztemp(len)       &! intermediate surface wind for roughness length
       ,ztemp              ! Used for ratio of reference height (zwind/ztemp)

  !respFactor is the annual total accumulation of carbon in the
  !previous year at each grid cell (annual total ASSIMN).
  !divided by the annual total of soilScale at the same grid pt.
  !respFactor*soilScale is the rate of release of CO2 by the soil.
  REAL :: respfactor(len,nzg_sib)
  !soilScale is a diagnostic of the instantaneous rate of
  !soil respiration (derived by Jim Collatz, similar to TEM)
  REAL :: soilscale(len,nzg_sib)
  REAL :: soilq10(len,nzg_sib)

  !MOSTLY VARIABLES FOR SURFACE FLUXES
  REAL ::                 & 
        fss(len)          &! surface sensible heat flux (W/m^2)
       ,fws(len)          &! surface evaporation (kg/m^2/s)
       ,cflux(len)        &! new formulation of CO2 flux (phosib)(mol/m^2/sec)
       ,cu(len)           &! momentum transfer coefficient (-)
       ,ct(len)           &! thermal transfer coefficient (-)
       ,ventmf(len)       &! ventilation mass flux (kg/m^2/sec)
       ,thvgm(len)        &! delta theta-v between atm and canopy
       ,bps(len)          &! (ps/1000)**kapa
       ,psy(len)          &! psychrometric 'constant'
       ,tha(len)          &! canopy airspace potential temperature (K)
       ,ea(len)           &! canopy airspace water vapor pressure (hPa)
       ,em(len)           &! mixed layer water vapor pressure (hPa)
       ,d(len)            &! dd corrected for snow covered canopy
       ,rbc(len)          &! cc1 corrected for snow covered canopy
       ,rdc(len)          &! cc2 corrected for snow covered canopy
       ,etmass(len)       &! evapotranspiration
       ,egmass(len)       &! ground evaporation (mm)
       ,ecmass(len)       &! canopy evaporation (mm)
       ,totwb(len)        &! total surface and soil water at begiN of timestep
       ,chf(len)          &! canopy heat flux (W/m^2)
       ,ahf(len)          &! CAS heat flux (W/m^2)
       ,shf(len)          &! soil heat flux (W/m^2)
       ,ect(len)          &! transpiration flux (J m^-2 for the timestep)
       ,eci(len)          &! canopy interception evap flux (veg-CAS) (J m^-2)
       ,egs(len)          &! soil/ground evaporation flux (J m^-2)
       ,egi(len)          &! ground interception evaporation flux (J m^-2)
       ,hc(len)           &! canopy (veg) sensible heat flux (W/m^2)
       ,hg(len)           &! ground surface sensible heat flux (W/m^2)
       ,hs(len)           &! snow surface sensible heat flux (W/m^2)
       ,heaten(len)       &! energy to heat snow to ground temp (J m^-2)
       ,hflux(len)        &! sensible heat flux (W/m2)
       ,tgs(len)          &! bare ground and snow surface mean temperature (K)
       ,tsnow(len)        &! Snow temp is lesser of ice or ground temp (K)
       ,czc(len)          &! canopy heat capacity
       ,etc(len)          &! vapor pressure (e*) of the canopy at (Tc) (Pa)
       ,etg(len)          &! vapor pressure (e*) of the ground sfc at (Tg) (Pa)
       ,etgs(len)         &! function result = E(TGS(i))
       ,btc(len)          &! derivatives of ETC
       ,btg(len)          &! derivatives of ETG
       ,rstfac(len,4)     &! canopy resistance stress factors
       !GROUND AND CANOPY STRESS AND PHOSIB VARIABLES
       ,rsoil(len)        &! SOIL SURFACE RESISTANCE (S M-1)
       ,hr(len)           &! SOIL SURFACE RELATIVE HUMIDITY
       ,wc(len)           &! CANOPY WETNESS FRACTION
       ,wg(len)           &! GROUND WETNESS FRACTION
       ,areas(len)        &! fractional snow coverage (0 to 1)
       ,csoil(len)        &! soil heat capacity (J m^-2 deg^-1)
       ,slamda(len,nsoil) &! soil thermal conductivities
       ,shcap(len,nsoil)  &! soil heat capacities
       ,czh(len)          &! surface layer heat capacity
       ,gect(len)         &! dry fraction of canopy/(Rst + 2Rb)
       ,geci(len)         &! wetted fraction of canopy/2Rb
       ,gegs(len)         &! dry fraction of ground/(fg*rsoil + Rd)
       ,gegi(len)         &! wet fraction of ground/Rd
       ,rb(len)           &! leaf sfc to CAS aerodynmaic resistance (sec/m)
       ,rd(len)           &! ground to CAS aerodynamic resistance (sec/m)
       ,rds(len)          &! rsoil + rd
       ,ra(len)           &! CAS to RAMS atmos aerodynamic resistance (sec/m)
       ,rib(len)          &! bulk richardson number
       ,rc(len)           &! total canopy resistance (sec/m)
       ,ggl(len)          &! overall leaf conductance
       ,hrr(len)          &! SOIL SURFACE LAYER RELATIVE HUMIDITY
       ,bintc(len)        &! (B*ZLT)  : EQUATION (35) , SE-92A
       ,aparkk(len)       &! (PI)     : EQUATION (31) , SE-92A
       ,wsfws(len)        &! Water stress
       ,wsfht(len)        &! High temperature stress
       ,wsflt(len)        &! Low temperature stress
       ,wci(len)          &! Intermediate assimilation weighted Ci
       ,whs(len)          &! Intermediate assimilation weighted RH stress factor
       ,wags(len)         &! Intermediate assimilation weighted stomatal conductance
       ,wegs(len)         &! Intermediate evaporation weighted stomatal conductance
       ,omepot(len)       &! Potential light limitation
       ,assim(len)        &! gross primary productivity (mol/m^2/s)
       ,assimpot(len)     &! Final Potential top leaf photosynthesis
       ,assimci(len)      &! Final stress limited top leaf photosynthesis
       ,assimnp(len)      &! Make assimn a top leaf, not the canopy
       ,antemp(len)       &! Bottom stopped assimilation
       ,ansqr(len)        &! Bottom stopped assimilation (squared)
       ,pfd(len)          &! 4.6E-6 * GMUDMU * (RADN(i,1,1)+RADN(i,1,2))
       ,zmstscale(len,2)  &! soil scaling parameter for shallow and root zone
       ,drst(len)         &! stomatal resistance increment
       !SENSIBLE HEAT FLUX DERIVATIVES AND SUCH
       ,dtg(len,2)        &! surface ground and snow temperature increments (K)
       ,dtc(len)          &! canopy temperature increment (K)
       ,dta(len)          &! CAS temperature increment (K)
       ,hgdtg(len)        &! dHG/dTG
       ,hgdta(len)        &! dHg/dTa
       ,hsdts(len)        &! dHS/dTS
       ,hsdta(len)        &! dHS/dTA
       ,hcdtc(len)        &! dHc/dTc
       ,hcdta(len)        &! dHc/dTa
       ,hadta(len)        &! dHA/dTA
       ,fc(len)           &! canopy range function? (not used here)
       ,fg(len)           &! ground humidity range function (0 or 1)
       !LONGWAVE RADIATIVE HEAT FLUX DERIVATIVES AND SUCH
       ,lcdtc(len)        &! dLC/dTC
       ,lcdtg(len)        &! dLC/dTG
       ,lcdts(len)        &! dLC/dTS
       ,lgdtg(len)        &! dLG/dTG
       ,lgdtc(len)        &! dLG/dTC
       ,lsdts(len)        &! dLS/dTS
       ,lsdtc(len)        &! dLS/dTC
       !LATENT HEAT FLUX GROUND AND CANOPY PARTIAL DERIVATIVE AND SUCH
       ,eg(len)           &! EGS + EGI
       ,ec(len)           &! ECT + ECI
       ,es(len)           &
       ,egdtg(len)        &! dEG/dTGS
       ,ecdtg(len)        &! dEC/dTGS
       ,ecdtc(len)        &! dEC/dTC
       ,egdtc(len)        &! dEG/dTC
       ,ecdea(len)        &! for the canopy leaves vapor pressure: W/ (m2* K)
       ,egdea(len)        &! for ground latent heat fluxes: W/ (m2* K)
       ,esdts(len)        &! for snow latent heat fluxes: W/ (m2* K)
       ,esdea(len)        &! for snow latent heat fluxes: W/ (m2 * Pa)
       ,eadea(len)        &! for CAS latent heat fluxes: W/ (m2* Pa)
       ,radt(len,3)       &! canopy, ground, and snow net radiation (W/m2)
       !INCREMENTS FOR INTEGRATION
       ,dea(len)          &! CAS moisture increment (Pa)
       ,dtd(len,nsoil)    &! deep soil temperature increments (K)
       ,q3l(len)          &! 'Liston' drainage from bottom of soillayer 3 (mm)
       ,q3o(len)          &! gravitational drainage out of soillayer 3 (mm)
       ,qqq(len,3)        &! soil layer drainage (mm m^-2 timestep)
       ,evt(len)          &
       ,eastar(len)       &! canopy saturation vapor pressure (hPa)
       ,rha(len)          &! canopy airspace relative humidity (%)
       !SURFACE WATER VARS
       ,zmelt(len)        &! total depth of melted water (m)
       ,zmelt1(len)       &! depth of melted water, main calculation (updat2)(m)
       ,zmelt2(len)       &! depth of melted water, from excess energy (inter2)(m)
       ,satcap(len,2)     &! saturation capacity of ground and vegetation (kg/m^2)
       ,exo(len)          &! total soil water excess of saturation (m)
       ,roffo(len)        &! runoff overland flow contribution (m)
       ,roff(len)         &! total runoff (surface and subsurface) (mm)
       !RADIATION VARIABLES
       ,salb(len,2,2)     &! surface albedos
       ,valb(len,2,2)     &! vegetation albedos
       ,nalb(len,2,2)     &! non-veg albedos
       ,canex(len)        &! [1.-( SNOWw(i,1)*5.-Z1(i))/(Z2(i)-Z1(i))]
       ,fac1(len)         &! effective ground cover for thermal radiation
       ,thgeff(len)       &! Tgeff(I) / BPS(I)
       ,tgeff(len)        &! effective (combined) skin temp from sfc thermal rad(K)
       ,shgeff(len)       &! saturation mixing ratio w.r.t tgeff
       ,tgeff4(len)       &! effective surface radiative temperature (K)
       ,radfac(len,2,2,2) &! radiation absorption factors
       ,radvbc(len)       &! surface incident visible direct beam (W/m^2)
       ,radnbc(len)       &! surface incident near IR direct beam (W/m^2)
       ,radvdc(len)       &! surface incident visible diffuse beam (W/m^2)
       ,radndc(len)       &! surface incident near IR diffuse beam (W/m^2)
       ,radc3(len,2)      &! SUM OF ABSORBED RADIATIVE FLUXES (W M-2)
       ,radn(len,2,2)     &! INCIDENT RADIATION FLUXES (W M-2)
       ,closs(len)        &! vegetation IR loss
       ,gloss(len)        &! ground IR loss
       ,sloss(len)        &! snow IR loss
       ,dtc4(len)         &! 1st derivative of vegetation T^4
       ,dtg4(len)         &! 1st derivative of ground T^4
       ,dts4(len)         &! 1st derivative of snow T^4
       !RESPIRATION VARIABLES
       ,assimn(len)       &! net co2 assimilation by plants (mol/m^2/sec)
       ,respg(len)        &! ground respiration flux (mol/m^2/sec)
       ,pco2i(len)        &! leaf internal pCO2 (Pa)
       ,pco2c(len)        &! chloroplast pCO2 (Pa)
       ,pco2s(len)        &! leaf surface pCO2 (Pa)
       ,pco2m             &! lowest atm level CO2 (Pa)
       ,co2cap(len)        ! moles of air in the canopy (moles/CAS)

  TYPE biome_morph_var
     REAL zc        ! Canopy inflection height (m)
     REAL lwidth    ! Leaf width
     REAL llength   ! Leaf length
     REAL laimax    ! Maximum LAI
     REAL stems     ! Stem area index
     REAL ndvimax   ! Maximum NDVI
     REAL ndvimin   ! Minimum NDVI
     REAL srmax     ! Maximum simple ratio
     REAL srmin     ! Minimum simple ratio
  END TYPE biome_morph_var
  TYPE(biome_morph_var) morphtab

  TYPE aero_var
     REAL zo       ! Canopy roughness coeff
     REAL zp_disp  ! Zero plane displacement
     REAL rbc      ! RB Coefficient
     REAL rdc      ! RC Coefficient
  END TYPE aero_var

  ! aerodynamic interpolation tables
  TYPE(aero_var),DIMENSION(50,50) :: aerovar

  TYPE time_dep_var
     REAL fpar    ! Canopy absorbed fraction of PAR
     REAL lai     ! Leaf-area index
     REAL green   ! Canopy greeness fraction of LAI
     REAL zo      ! Canopy roughness coeff
     REAL zp_disp ! Zero plane displacement
     REAL rbc     ! RB Coefficient (c1)
     REAL rdc     ! RC Coefficient (c2)
     REAL gmudmu  ! Time-mean leaf projection
  END TYPE time_dep_var
  TYPE(time_dep_var) timevar

  !Send in vegetation biosphere and soil class info and coordinate
  !with SiB equivalent classes of vegetation and soil. Note that
  !RAMS clases of vegetation go from 0:20 or 21 total classes, but
  !RAMS LEAF class=0 is ocean. Do not need this for SiB.
  !Zero is water. Should not be running SiB for water patches.
  !Alert model to report error.
  INTEGER :: biome   &  ! biome type, sent in as a single value
            ,soiltype   ! soil type, sent in as a single value
  INTEGER, DIMENSION(0:20) :: leaf_biome_map
  DATA leaf_biome_map /0,0,13,11,4,5,2,1,6,6,9,10,9,9,3,12,12,12,7,11,1/

  !LEAF soil classes do not match the order of classes in SiB so we need
  !to translate these correctly. LEAF uses the classes from Clapp & Hornberger
  !and SiB uses modified USDA classes. LEAF use of FAO data only assigns 
  !classes 2-8 (Clapp & Hornberger) via routine "datp_datsoil". Others
  !default to sandy-clay-loam.
  INTEGER, DIMENSION(1:12) :: soil_type_map
  DATA soil_type_map /1,2,3,4,6,7,10,9,8,11,12,7/

print*,'***************************************************************'
print*,'  SIB MODEL NOT YET AVAILABLE IN THIS RELEASE' 
print*,'***************************************************************'
stop

return
END SUBROUTINE sib_2pt5
