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

Subroutine aerorad (mb,nb,nrad,m1,dzl,relh,tp,omgp,gp,dn0)

! ib .......... band number
! mb .......... maximum number of bands (8)
! nb .......... total number of bands
! m1 .......... number of vertical levels
! dzl .......... delta z in each level (m)
! nrad - number of vertical levels in radiation calculations
! naerbins - number of aerosol bins for radiative calculations (17)
! nrh - number of relative humidity levels for radiative calculations (20)
! relh(m1) - relative humidity at each vertical level
! iaerorad - flag to turn on aerosol radiation
! aerocon(m1,aerocat) - aerosol category number concentration (#/kg)
! aeromas(m1,aerocat) - aerosol category mass (kg/kg)
! irh     - the relative humidity level (1=80%,20=99%)
! aerocat  - number of aerosol species in the model
! aerotype - aerosol type to be used through the radiation code
! rnaer(m1,naerbins) - aerosol number concentration in a given bin (#/m3)
!***************************************************************************
use micphys

implicit none

! Variable declaration:
INTEGER mb,nb,ib,nrad,m1,ibns,naerbins,irh,nrh,k
REAL dzl(nrad),tp(nrad,mb),omgp(nrad,mb),gp(nrad,mb),relh(m1)
REAL rma(17),ext_sum,om_sum,asym_sum,con_sum,ext,om,gg
REAL gfact,qext,qscat,gasym,dm,pi

PARAMETER(naerbins=17,nrh=20)

INTEGER aerotype(aerocat)
REAL rnaer(m1,naerbins),aeroradius(m1,aerocat),dn0(m1)

PARAMETER (pi=3.141593)

! Drop radius at bin center in centimeters
DATA rma /2.0358125E-08,3.2316509E-08,5.1299262E-08                &
         ,8.1432499E-08,1.2926605E-07,2.0519707E-07,3.2573003E-07  &
         ,5.1706418E-07,8.2078826E-07,1.3029201E-06,2.0682569E-06  &
         ,3.2831531E-06,5.2116811E-06,8.2730281E-06,1.3132613E-05  &
         ,2.0846725E-05,3.3092112E-05/

! Initializing aerotype:
! Although there are multiple aerosol variables, there are currently only
! three types that are coded. Aerotype indicies are as follows:
! 1 - Ammonium sulfate
! 2 - Sea Salt
! 3 - Mineral Dust
do acat=1,aerocat
  aerotype(acat) = 0
enddo

if(iaerosol>0) then
   aerotype(1) = 1   ! CCN
   aerotype(2) = 1   ! GCCN
endif
if(idust>0) then
   aerotype(3) = 3   ! Small dust mode
   aerotype(4) = 3   ! Large dust mode
endif
if(isalt>0) then
   aerotype(5) = 2   ! Salt film mode
   aerotype(6) = 2   ! Salt jet mode
   aerotype(7) = 2   ! Salt spume mode
endif
if(iccnlev>=2) then
   aerotype(8) = 1   ! Small regenerated aerosol
   aerotype(9) = 1   ! Large regenerated aerosol
endif

!Loop thru all aerosol species and set number and size
do acat=1,aerocat
 if(aerotype(acat)>0) then
  do k = 2,m1-1
   if(aerocon(k,acat)>mincon .and. aeromas(k,acat)>=minmas) then
    aeroradius(k,acat)=((0.23873/aero_rhosol(acat) &
        *aeromas(k,acat)/aerocon(k,acat))**(1./3.))/aero_rg2rm(acat)
   else
    aeroradius(k,acat)=0.005e-6
   endif
  enddo
 endif
enddo

! Looping over each aerosol species. If the aerosol is turned off, the
! code checks for the next active aerosol:
DO acat=1,aerocat
 IF(aerotype(acat).gt.0) then
  ! For each aerosol, the code requires a total number concentration of 
  ! the aerosol at each grid point, along with either a total mass at the
  ! grid point or a median radius.
  CALL aerobin (m1,naerbins,aerocon(1,acat),aeroradius(1,acat),rnaer,dn0)
  ! This is the routine that calculates the optical properties of the 
  ! aerosols within the radiation routine:
  ! Doing over all vertical levels (not including radiation levels)
  DO k = 2,m1-1
    ! Locate RH as a percentage for table lookup
    irh=INT(100.*(relh(k)-0.80)) + 1
    irh=MAX(1,MIN(irh,20) )

    DO ib=1,nb

      ! Resetting temporary storage variabes to zero for next set of calcs:
      ext_sum = 0.0
      om_sum = 0.0
      asym_sum = 0.0
      con_sum = 0.0
      ext = 0.0
      om  = 0.0
      gg  = 0.0

      ! Calculating optical properties at over each aerosol bin:
      DO ibns = 1,naerbins
        ! Retrieving values for the growth factor, the extinction coefficient,
        ! the scattering coefficient and the asymmetry parameter, from their
        ! appropriate look up tables:
        CALL aerogrowth (aerotype(acat),ibns,irh,gfact)
        CALL aeroqext   (aerotype(acat),ib,ibns,irh,qext)
        CALL aeroqscat  (aerotype(acat),ib,ibns,irh,qscat)
        CALL aerogasym  (aerotype(acat),ib,ibns,irh,gasym)
        ! Determining mean diameter of bin (meters), with deliquescence growth factor
        dm = 2. * rma(ibns) * gfact
        ! Updating temporary storage variables:
        if(qext>0.0 .and. rnaer(k,ibns)>0.0)then
          ext_sum  = ext_sum  + (pi/4. * dm**2 * rnaer(k,ibns) * qext)  !Bext
          om_sum   = om_sum   + (pi/4. * dm**2 * rnaer(k,ibns) * qscat) !Bscat
          asym_sum = asym_sum + (rnaer(k,ibns) * gasym) !Asymmetry parameter
          con_sum  = con_sum  + rnaer(k,ibns) !Total aerosol number
        endif
      ENDDO ! Aerosol bins (ibns)
      !Set up ext, om, gg same as Harrington's cloud_top for consistency
      ext = ext_sum * dzl(k) !Tau = Bext * dz
      if(ext>0.0) then
        om   = om_sum / ext_sum  !Omega = Bscat / Bext
        gg   = asym_sum/(con_sum+1.E-30) !normalize asym_sum by total number
      else
        ext = 0.0
      endif
      ! Obtain values of optical depth (tp), single scatter albedo (omgp),
      ! asymmetry parameter (gp) in a form to interact with cloud_opt.
      ! Some values of dzl may be <0; in this case set tp to 0.
      ! Factor in layer depth for optical depth extinction
      tp(k,ib)   = tp(k,ib)    + ext
      omgp(k,ib) = omgp(k,ib)  + om * ext
      gp(k,ib)   = gp(k,ib)    + gg * om * ext
    ENDDO ! Radiation band loop
  ENDDO ! Vertical level loop
 ENDIF !Aerosol type if
ENDDO ! Aerotype loop

return
END SUBROUTINE aerorad

!##############################################################################
Subroutine aerobin (m1,naerbins,numconc,medianrad,rnaer,dn0)

! The purpose of this routine was to take a number concentration
! and median radius for any given aerosol and turn it into a binned
! distribution
! Input Variables: m1 - Integer, number of vertical levels
!    naerbins  - Integer number of bins used in aerosol radiation
!                calculations.  Set to 17 in radcalc3.
!    numconc   - Number concentration of aerosol (#/kg)
!    medianrad - Median radius of the given aerosol type (m)
!    interp    - Determines beginning point for interpolation calculations
!    radlim    - Radii at which aerosol distribution was explicitly
!                  binned (in cm)
!    radparam  - Percentage of aerosol mass in each of the 17 bins,
!                as a func of median radius.
!    rg        - Median radius of aerosol (m)
!    L(4)      - Values used in the 3rd order Lagrangian polynomial
!                interpolation scheme to determine percentage of aerosol
!                within a given bin
!    m         - Used to define which four median radii should
!                 be used in interpolation
! Output variable: rnaer(m1,naerbins) - Number concentrations for
!                     all 17 bins at every vertical level (#/m3)
!**************************************************************************

implicit none

! Incoming variable declaration:
        INTEGER m1,naerbins
        REAL numconc(m1),medianrad(m1),dn0(m1)

! Internal variable declaraion:
        INTEGER bin,i,k,n,m,interp
        REAL radlim(26),radparam(26,17),rg,L(4)

! Outgoing variable declaration:
        REAL rnaer(m1,naerbins)

! Data statement containing the radii for which the binning routine was
! run off-line...i.e., explicit binning done at these radii:
        DATA radlim/0.00945E-06, 0.01191E-06, 0.01500E-06, 0.01890E-06  &
                   ,0.02382E-06, 0.03000E-06, 0.03780E-06, 0.04760E-06  &
                   ,0.06000E-06, 0.07560E-06, 0.09520E-06, 0.12000E-06  &
                   ,0.15120E-06, 0.19050E-06, 0.24000E-06, 0.30240E-06  &
                   ,0.38100E-06, 0.48000E-06, 0.60500E-06, 0.76200E-06  &
                   ,0.96000E-06, 1.20900E-06, 1.52400E-06, 1.92000E-06  &
                   ,2.41900E-06, 3.04800E-06/ ! in meters

! Percentage of mass within each bin at the above given radii:
        DATA (radparam(1,bin),bin=1,17)  & ! at 9.45nm
        /0.1580, 0.0489, 0.0084, 0.0008, 0.0000, 0.0000  &
        ,0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000  &
        ,0.0000, 0.0000, 0.0000, 0.0000, 0.0000/

        DATA (radparam(2,bin),bin=1,17)  & ! at 11.91nm
        /0.2280, 0.0947, 0.0218, 0.0028, 0.0002, 0.0000  &
        ,0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000  &
        ,0.0000, 0.0000, 0.0000, 0.0000, 0.0000/

        DATA (radparam(3,bin),bin=1,17)  & ! at 15nm
        /0.2839, 0.1580, 0.0488, 0.0083, 0.0008, 0.0000  &
        ,0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000  &
        ,0.0000, 0.0000, 0.0000, 0.0000, 0.0000/

        DATA (radparam(4,bin),bin=1,17)  & ! at 18.9nm
        /0.3055, 0.2280, 0.0945, 0.0217, 0.0028, 0.0002  &
        ,0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000  &
        ,0.0000, 0.0000, 0.0000, 0.0000, 0.0000/

        DATA (radparam(5,bin),bin=1,17)  & ! at 23.82nm
        /0.2839, 0.2841, 0.1581, 0.0487, 0.0084, 0.0008  &
        ,0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000  &
        ,0.0000, 0.0000, 0.0000, 0.0000, 0.0000/

        DATA (radparam(6,bin),bin=1,17)  & ! at 30nm
        /0.2279, 0.3056, 0.2280, 0.0943, 0.0217, 0.0028  &
        ,0.0002, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000  &
        ,0.0000, 0.0000, 0.0000, 0.0000, 0.0000/

        DATA (radparam(7,bin),bin=1,17)  & ! at 37.8nm
        /0.1579, 0.2840, 0.2841, 0.1577, 0.0488, 0.0084  &
        ,0.0008, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000  &
        ,0.0000, 0.0000, 0.0000, 0.0000, 0.0000/

        DATA (radparam(8,bin),bin=1,17)  & ! at 47.6nm
        /0.0946, 0.2280, 0.3057, 0.2276, 0.0943, 0.0217  &
        ,0.0028, 0.0002, 0.0000, 0.0000, 0.0000, 0.0000  &
        ,0.0000, 0.0000, 0.0000, 0.0000, 0.0000/

        DATA (radparam(9,bin),bin=1,17)  & ! at 60nm
        /0.0488, 0.1579, 0.2481, 0.2838, 0.1579, 0.0488  &
        ,0.0084, 0.0008, 0.0000, 0.0000, 0.0000, 0.0000  &
        ,0.0000, 0.0000, 0.0000, 0.0000, 0.0000/

        DATA (radparam(10,bin),bin=1,17)  & ! at 75.6nm
        /0.0218, 0.0945, 0.2280, 0.3055, 0.2278, 0.0946  &
        ,0.0218, 0.0028, 0.0002, 0.0000, 0.0000, 0.0000  &
        ,0.0000, 0.0000, 0.0000, 0.0000, 0.0000/

        DATA (radparam(11,bin),bin=1,17)  & ! at 95.2nm
        /0.0084, 0.0489, 0.1581, 0.2841, 0.2838, 0.1578  &
        ,0.0488, 0.0083, 0.0008, 0.0000, 0.0000, 0.0000  &
        ,0.0000, 0.0000, 0.0000, 0.0000, 0.0000/

        DATA (radparam(12,bin),bin=1,17)  & ! at 120nm
        /0.0028, 0.0217, 0.0945, 0.2280, 0.3055, 0.2279  &
        ,0.0946, 0.0217, 0.0028, 0.0002, 0.0000, 0.0000  &
        ,0.0000, 0.0000, 0.0000, 0.0000, 0.0000/

        DATA (radparam(13,bin),bin=1,17)  & ! at 151.2nm
        /0.0008, 0.0084, 0.0488, 0.1580, 0.2839, 0.2840  &
        ,0.1581, 0.0487, 0.0083, 0.0008, 0.0000, 0.0000  &
        ,0.0000, 0.0000, 0.0000, 0.0000, 0.0000/

        DATA (radparam(14,bin),bin=1,17)  & ! at 190.5nm
        /0.0002, 0.0028, 0.0217, 0.0945, 0.2279, 0.3056  &
        ,0.2281, 0.0944, 0.0217, 0.0028, 0.0002, 0.0000  &
        ,0.0000, 0.0000, 0.0000, 0.0000, 0.0000/

        DATA (radparam(15,bin),bin=1,17)  & ! at 240nm
        /0.0000, 0.0008, 0.0083, 0.0488, 0.1580, 0.2840  &
        ,0.2841, 0.1578, 0.0487, 0.0084, 0.0008, 0.0000  &
        ,0.0000, 0.0000, 0.0000, 0.0000, 0.0000/

        DATA (radparam(16,bin),bin=1,17)  & ! at 302.4 nm
        /0.0000, 0.0002, 0.0028, 0.0218, 0.0945, 0.2279  &
        ,0.3057, 0.2278, 0.0944, 0.0218, 0.0028, 0.0002  &
        ,0.0000, 0.0000, 0.0000, 0.0000, 0.0000/

        DATA (radparam(17,bin),bin=1,17)  & ! at 381nm
        /0.0000, 0.0000, 0.0008, 0.0084, 0.0488, 0.1579  &
        ,0.2841, 0.2840, 0.1578, 0.0488, 0.0084, 0.0008  &
        ,0.0000, 0.0000, 0.0000, 0.0000, 0.0000/

        DATA (radparam(18,bin),bin=1,17)  & ! at 480nm
        /0.0000, 0.0000, 0.0002, 0.0028, 0.0215, 0.0945  &
        ,0.2279, 0.3056, 0.2277, 0.0945, 0.0218, 0.0028  &
        ,0.0002, 0.0000, 0.0000, 0.0000, 0.0000/

        DATA (radparam(19,bin),bin=1,17)  & ! at 605nm
        /0.0000, 0.0000, 0.0000, 0.0008, 0.0084, 0.0487  &
        ,0.1587, 0.2840, 0.2839, 0.1580, 0.0489, 0.0084  &
        ,0.0008, 0.0000, 0.0000, 0.0000, 0.0000/

        DATA (radparam(20,bin),bin=1,17)  & ! at 762nm
        /0.0000, 0.0000, 0.0000, 0.0002, 0.0028, 0.0217  &
        ,0.0944, 0.2280, 0.3055, 0.2279, 0.0946, 0.0218  &
        ,0.0028, 0.0002, 0.0000, 0.0000, 0.0000/

        DATA (radparam(21,bin),bin=1,17)  & ! at 960nm
        /0.0000, 0.0000, 0.0000, 0.0000, 0.0008, 0.0084  &
        ,0.0488, 0.1580, 0.2840, 0.2839, 0.1580, 0.0488  &
        ,0.0083, 0.0008, 0.0000, 0.0000, 0.0000/

        DATA (radparam(22,bin),bin=1,17)  & ! at 1.209 microns
        /0.0000, 0.0000, 0.0000, 0.0000, 0.0002, 0.0028  &
        ,0.0218, 0.0946, 0.2281, 0.3055, 0.2278, 0.0945  &
        ,0.0216, 0.0028, 0.0002, 0.0000, 0.0000/

        DATA (radparam(23,bin),bin=1,17)  & ! at 1.524 microns
        /0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0008  &
        ,0.0083, 0.0488, 0.1580, 0.2839, 0.2840, 0.1581  &
        ,0.0486, 0.0084, 0.0008, 0.0000, 0.0000/

        DATA (radparam(24,bin),bin=1,17)  & ! at 1.92 microns
        /0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0002  &
        ,0.0028, 0.0217, 0.0945, 0.2279, 0.3056, 0.2281  &
        ,0.0943, 0.0217, 0.0028, 0.0002, 0.0000/

        DATA (radparam(25,bin),bin=1,17)  & ! at 2.419 microns
        /0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000  &
        ,0.0008, 0.0084, 0.0488, 0.1580, 0.2840, 0.2842  &
        ,0.1577, 0.0487, 0.0084, 0.0008, 0.0000/

        DATA (radparam(26,bin),bin=1,17)  & ! at 3.048 microns
        /0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000  &
        ,0.0002, 0.0028, 0.0218, 0.0945, 0.2279, 0.3058  &
        ,0.2277, 0.0944, 0.0218, 0.0028, 0.0002/

! Vertical level loop:
  DO k = 1,m1

          rg = medianrad(k)

          !Determining which four points to use in the interpolation scheme
          interp = 0
          DO n=1,25
            IF(rg.GE.radlim(n).AND.rg.LE.radlim(n+1)) THEN
              interp = n
            ENDIF
          ENDDO
          if(interp == 0) then
            IF(rg.LT.radlim(1))   interp = 1
            IF(rg.GT.radlim(25)) interp = 25
          endif
          m = interp - 1
          IF(interp.EQ.1) m = 1
          IF(interp.EQ.25) m = 23

          ! Determining the coefficients to be used in the 3rd order 
          !Lagrangian polynomial interpolation scheme:
          L(1) = (rg-radlim(m+3))*(rg-radlim(m+2))*  &
                 (rg-radlim(m+1))/(radlim(m)-radlim(m+3))/  &
                 (radlim(m)-radlim(m+2))/(radlim(m)-radlim(m+1))
          L(2) = (rg-radlim(m+3))*(rg-radlim(m+2))*  &
                 (rg-radlim(m))/(radlim(m+1)-radlim(m+3))/  &
                 (radlim(m+1)-radlim(m+2))/(radlim(m+1)-radlim(m))
          L(3) = (rg-radlim(m+3))*(rg-radlim(m+1))*  &
                 (rg-radlim(m))/(radlim(m+2)-radlim(m+3))/  &
                 (radlim(m+2)-radlim(m+1))/(radlim(m+2)-radlim(m))
          L(4) = (rg-radlim(m+2))*(rg-radlim(m+1))*  &
                 (rg-radlim(m))/(radlim(m+3)-radlim(m+2))/  &
                 (radlim(m+3)-radlim(m+1))/(radlim(m+3)-radlim(m))

          !At each bin, determine the number concentration of aerosol
          DO i = 1,naerbins
            rnaer(k,i) = L(1)*radparam(m,i) + L(2)*radparam(m+1,i) +  &
                         L(3)*radparam(m+2,i) + L(4)*radparam(m+3,i)
            IF(rnaer(k,i).LT.0.) rnaer(k,i) = 0.
            ! Conversion from #/kg to #/m3:
            rnaer(k,i) = rnaer(k,i) * numconc(k) * dn0(k)
          ENDDO ! Bin loop (i)

  ENDDO   ! Vertical level loop (k)

return
END SUBROUTINE aerobin

!##############################################################################
Subroutine aerogrowth (aerotype,bin,rh,value)

implicit none

        INTEGER aerotype,bin,rh,ib
        REAL growth(3,17,20),value

! Ammonium sulfate with 90% insoluble inclusion (dust)
        DATA(growth(1,ib,1),ib=1,17)  & ! Amm. Sulf. at  80% RH
        /1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000  &
        ,1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000  &
        ,1.0000, 1.0000, 1.0000, 1.0000, 1.0000/

        DATA(growth(1,ib,2),ib=1,17)  & ! Amm. Sulf. at  81% RH
        /1.0024, 1.0027, 1.0030, 1.0031, 1.0033, 1.0033  &
        ,1.0034, 1.0034, 1.0034, 1.0034, 1.0035, 1.0035  &
        ,1.0035, 1.0035, 1.0035, 1.0035, 1.0035/

        DATA(growth(1,ib,3),ib=1,17)  & ! Amm. Sulf. at  82% RH
        /1.0050, 1.0057, 1.0062, 1.0066, 1.0068, 1.0070  &
        ,1.0071, 1.0072, 1.0072, 1.0073, 1.0073, 1.0073  &
        ,1.0073, 1.0073, 1.0073, 1.0073, 1.0073/

        DATA(growth(1,ib,4),ib=1,17)  & ! Amm. Sulf. at  83% RH
        /1.0078, 1.0089, 1.0098, 1.0104, 1.0108, 1.0111  &
        ,1.0112, 1.0114, 1.0114, 1.0115, 1.0115, 1.0115  &
        ,1.0115, 1.0115, 1.0115, 1.0115, 1.0115/

        DATA(growth(1,ib,5),ib=1,17)  & ! Amm. Sulf. at  84% RH
        /1.0108, 1.0125, 1.0138, 1.0146, 1.0152, 1.0156  &
        ,1.0159, 1.0160, 1.0161, 1.0162, 1.0162, 1.0162  &
        ,1.0163, 1.0163, 1.0163, 1.0163, 1.0163/

        DATA(growth(1,ib,6),ib=1,17)  & ! Amm. Sulf. at  85% RH
        /1.0142, 1.0165, 1.0182, 1.0193, 1.0201, 1.0207  &
        ,1.0210, 1.0212, 1.0214, 1.0215, 1.0215, 1.0216  &
        ,1.0216, 1.0216, 1.0216, 1.0216, 1.0216/

        DATA(growth(1,ib,7),ib=1,17)  & ! Amm. Sulf. at  86% RH
        /1.0179, 1.0209, 1.0231, 1.0246, 1.0257, 1.0264  &
        ,1.0268, 1.0271, 1.0273, 1.0274, 1.0275, 1.0276  &
        ,1.0276, 1.0276, 1.0276, 1.0276, 1.0276/

        DATA(growth(1,ib,8),ib=1,17)  & ! Amm. Sulf. at  87% RH
        /1.0220, 1.0258, 1.0286, 1.0306, 1.0320, 1.0329  &
        ,1.0335, 1.0339, 1.0341, 1.0343, 1.0344, 1.0344  &
        ,1.0345, 1.0345, 1.0345, 1.0345, 1.0345/

        DATA(growth(1,ib,9),ib=1,17)  & ! Amm. Sulf. at  88% RH
        /1.0266, 1.0313, 1.0349, 1.0374, 1.0392, 1.0404  &
        ,1.0411, 1.0416, 1.0419, 1.0421, 1.0422, 1.0423  &
        ,1.0424, 1.0424, 1.0424, 1.0424, 1.0425/

        DATA(growth(1,ib,10),ib=1,17)  & ! Amm. Sulf. at  89% RH
        /1.0317, 1.0376, 1.0421, 1.0453, 1.0475, 1.0490  &
        ,1.0500, 1.0506, 1.0510, 1.0513, 1.0514, 1.0515  &
        ,1.0516, 1.0516, 1.0517, 1.0517, 1.0517/

        DATA(growth(1,ib,11),ib=1,17)  & ! Amm. Sulf. at  90% RH
        /1.0375, 1.0448, 1.0504, 1.0545, 1.0573, 1.0592  &
        ,1.0604, 1.0612, 1.0617, 1.0621, 1.0623, 1.0624  &
        ,1.0625, 1.0625, 1.0626, 1.0626, 1.0626/

        DATA(growth(1,ib,12),ib=1,17)  & ! Amm. Sulf. at  91% RH
        /1.0441, 1.0531, 1.0602, 1.0653, 1.0689, 1.0713  &
        ,1.0729, 1.0739, 1.0746, 1.0750, 1.0752, 1.0754  &
        ,1.0755, 1.0756, 1.0756, 1.0756, 1.0757/

        DATA(growth(1,ib,13),ib=1,17)  & ! Amm. Sulf. at  92% RH
        /1.0518, 1.0629, 1.0718, 1.0783, 1.0829, 1.0860  &
        ,1.0880, 1.0894, 1.0902, 1.0907, 1.0911, 1.0913  &
        ,1.0914, 1.0915, 1.0916, 1.0916, 1.0916/

        DATA(growth(1,ib,14),ib=1,17)  & ! Amm. Sulf. at  93% RH
        /1.0606, 1.0745, 1.0859, 1.0943, 1.1002, 1.1042  &
        ,1.1069, 1.1086, 1.1097, 1.1104, 1.1108, 1.1111  &
        ,1.1113, 1.1114, 1.1115, 1.1115, 1.1116/

        DATA(growth(1,ib,15),ib=1,17)  & ! Amm. Sulf. at  94% RH
        /1.0711, 1.0887, 1.1033, 1.1143, 1.1221, 1.1275  &
        ,1.1310, 1.1333, 1.1347, 1.1357, 1.1363, 1.1367  &
        ,1.1369, 1.1371, 1.1371, 1.1372, 1.1372/

        DATA(growth(1,ib,16),ib=1,17)  & ! Amm. Sulf. at  95% RH
        /1.0837, 1.1063, 1.1255, 1.1402, 1.1508, 1.1581  &
        ,1.1630, 1.1661, 1.1682, 1.1695, 1.1703, 1.1708  &
        ,1.1711, 1.1714, 1.1715, 1.1716, 1.1716/

        DATA(growth(1,ib,17),ib=1,17)  & ! Amm. Sulf. at  96% RH
        /1.0993, 1.1289, 1.1548, 1.1753, 1.1903, 1.2007  &
        ,1.2077, 1.2123, 1.2152, 1.2171, 1.2183, 1.2190  &
        ,1.2195, 1.2198, 1.2200, 1.2201, 1.2202/

        DATA(growth(1,ib,18),ib=1,17)  & ! Amm. Sulf. at  97% RH
        /1.1190, 1.1591, 1.1959, 1.2260, 1.2486, 1.2644  &
        ,1.2752, 1.2822, 1.2868, 1.2898, 1.2916, 1.2928  &
        ,1.2936, 1.2940, 1.2943, 1.2945, 1.2947/

        DATA(growth(1,ib,19),ib=1,17)  & ! Amm. Sulf. at  98% RH
        /1.1451, 1.2022, 1.2585, 1.3071, 1.3448, 1.3720  &
        ,1.3907, 1.4032, 1.4113, 1.4165, 1.4199, 1.4220  &
        ,1.4233, 1.4242, 1.4247, 1.4250, 1.4253/

        DATA(growth(1,ib,20),ib=1,17)  & ! Amm. Sulf. at  99% RH
        /1.1816, 1.2705, 1.3696, 1.4645, 1.5442, 1.6047  &
        ,1.6476, 1.6766, 1.6958, 1.7082, 1.7162, 1.7213  &
        ,1.7245, 1.7265, 1.7278, 1.7286, 1.7291/

! Sea salt with no insoluble inclusion
        DATA(growth(2,ib,1),ib=1,17)  & ! Sea Salt at  80% RH
        /1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000  &
        ,1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000  &
        ,1.0000, 1.0000, 1.0000, 1.0000, 1.0000/

        DATA(growth(2,ib,2),ib=1,17)  & ! Sea Salt at  81% RH
        /1.0151, 1.0142, 1.0134, 1.0136, 1.0137, 1.0137  &
        ,1.0138, 1.0138, 1.0138, 1.0138, 1.0138, 1.0138  &
        ,1.0138, 1.0138, 1.0138, 1.0138, 1.0138/

        DATA(growth(2,ib,3),ib=1,17)  & ! Sea Salt at  82% RH
        /1.0295, 1.0281, 1.0276, 1.0279, 1.0281, 1.0282  &
        ,1.0283, 1.0283, 1.0283, 1.0284, 1.0284, 1.0284  &
        ,1.0284, 1.0284, 1.0284, 1.0284, 1.0284/

        DATA(growth(2,ib,4),ib=1,17)  & ! Sea Salt at  83% RH
        /1.0437, 1.0428, 1.0425, 1.0430, 1.0433, 1.0435  &
        ,1.0437, 1.0438, 1.0438, 1.0439, 1.0439, 1.0439  &
        ,1.0439, 1.0439, 1.0439, 1.0439, 1.0439/

        DATA(growth(2,ib,5),ib=1,17)  & ! Sea Salt at  84% RH
        /1.0588, 1.0584, 1.0584, 1.0592, 1.0597, 1.0600  &
        ,1.0602, 1.0603, 1.0604, 1.0605, 1.0605, 1.0605  &
        ,1.0606, 1.0606, 1.0606, 1.0606, 1.0606/

        DATA(growth(2,ib,6),ib=1,17)  & ! Sea Salt at  85% RH
        /1.0748, 1.0750, 1.0756, 1.0767, 1.0774, 1.0778  &
        ,1.0781, 1.0783, 1.0784, 1.0785, 1.0785, 1.0786  &
        ,1.0786, 1.0786, 1.0786, 1.0786, 1.0786/

        DATA(growth(2,ib,7),ib=1,17)  & ! Sea Salt at  86% RH
        /1.0920, 1.0931, 1.0942, 1.0957, 1.0966, 1.0973  &
        ,1.0977, 1.0979, 1.0981, 1.0982, 1.0982, 1.0983  &
        ,1.0983, 1.0983, 1.0983, 1.0983, 1.0983/

        DATA(growth(2,ib,8),ib=1,17)  & ! Sea Salt at  87% RH
        /1.1107, 1.1127, 1.1146, 1.1166, 1.1179, 1.1187  &
        ,1.1192, 1.1196, 1.1198, 1.1199, 1.1200, 1.1201  &
        ,1.1201, 1.1201, 1.1201, 1.1201, 1.1201/

        DATA(growth(2,ib,9),ib=1,17)  & ! Sea Salt at  88% RH
        /1.1312, 1.1345, 1.1372, 1.1398, 1.1415, 1.1426  &
        ,1.1433, 1.1437, 1.1440, 1.1442, 1.1443, 1.1444  &
        ,1.1444, 1.1444, 1.1444, 1.1445, 1.1445/

        DATA(growth(2,ib,10),ib=1,17)  & ! Sea Salt at  89% RH
        /1.1539, 1.1588, 1.1626, 1.1659, 1.1681, 1.1695  &
        ,1.1704, 1.1710, 1.1713, 1.1715, 1.1717, 1.1718  &
        ,1.1718, 1.1719, 1.1719, 1.1719, 1.1719/

        DATA(growth(2,ib,11),ib=1,17)  & ! Sea Salt at  90% RH
        /1.1794, 1.1862, 1.1914, 1.1956, 1.1984, 1.2002  &
        ,1.2013, 1.2020, 1.2025, 1.2028, 1.2029, 1.2030  &
        ,1.2031, 1.2032, 1.2032, 1.2032, 1.2032/

        DATA(growth(2,ib,12),ib=1,17)  & ! Sea Salt at  91% RH
        /1.2083, 1.2177, 1.2245, 1.2298, 1.2333, 1.2356  &
        ,1.2370, 1.2379, 1.2385, 1.2388, 1.2390, 1.2392  &
        ,1.2393, 1.2393, 1.2394, 1.2394, 1.2394/

        DATA(growth(2,ib,13),ib=1,17)  & ! Sea Salt at  92% RH
        /1.2417, 1.2541, 1.2631, 1.2698, 1.2742, 1.2770  &
        ,1.2788, 1.2799, 1.2807, 1.2811, 1.2814, 1.2816  &
        ,1.2817, 1.2818, 1.2818, 1.2819, 1.2819/

        DATA(growth(2,ib,14),ib=1,17)  & ! Sea Salt at  93% RH
        /1.2808, 1.2971, 1.3088, 1.3173, 1.3229, 1.3265  &
        ,1.3288, 1.3302, 1.3311, 1.3317, 1.3321, 1.3323  &
        ,1.3324, 1.3325, 1.3326, 1.3326, 1.3327/

        DATA(growth(2,ib,15),ib=1,17)  & ! Sea Salt at  94% RH
        /1.3275, 1.3489, 1.3642, 1.3751, 1.3823, 1.3869  &
        ,1.3899, 1.3917, 1.3929, 1.3937, 1.3941, 1.3944  &
        ,1.3946, 1.3947, 1.3948, 1.3948, 1.3949/

        DATA(growth(2,ib,16),ib=1,17)  & ! Sea Salt at  95% RH
        /1.3846, 1.4131, 1.4334, 1.4478, 1.4572, 1.4632  &
        ,1.4671, 1.4696, 1.4711, 1.4721, 1.4727, 1.4731  &
        ,1.4734, 1.4735, 1.4736, 1.4737, 1.4737/

        DATA(growth(2,ib,17),ib=1,17)  & ! Sea Salt at  96% RH
        /1.4570, 1.4959, 1.5235, 1.5431, 1.5559, 1.5642  &
        ,1.5696, 1.5729, 1.5751, 1.5764, 1.5773, 1.5778  &
        ,1.5782, 1.5784, 1.5785, 1.5786, 1.5786/

        DATA(growth(2,ib,18),ib=1,17)  & ! Sea Salt at  97% RH
        /1.5534, 1.6089, 1.6489, 1.6770, 1.6956, 1.7077  &
        ,1.7154, 1.7204, 1.7235, 1.7255, 1.7267, 1.7275  &
        ,1.7280, 1.7283, 1.7285, 1.7287, 1.7287/

        DATA(growth(2,ib,19),ib=1,17)  & ! Sea Salt at  98% RH
        /1.6927, 1.7790, 1.8425, 1.8875, 1.9175, 1.9371  &
        ,1.9497, 1.9578, 1.9629, 1.9662, 1.9682, 1.9695  &
        ,1.9704, 1.9709, 1.9712, 1.9714, 1.9715/

        DATA(growth(2,ib,20),ib=1,17)  & ! Sea Salt at  99% RH
        /1.9269, 2.0871, 2.2119, 2.3031, 2.3656, 2.4071  &
        ,2.4341, 2.4515, 2.4626, 2.4696, 2.4741, 2.4769  &
        ,2.4787, 2.4798, 2.4805, 2.4809, 2.4812/

        value = growth(aerotype,bin,rh)

! No growth factor for dust at RH under 100%:
        IF(AEROTYPE.EQ.3) value = 1.0

return
END SUBROUTINE aerogrowth

!##############################################################################
Subroutine aeroqext (aerotype,radband,bin,rh,value)

implicit none

        INTEGER aerotype,radband,bin,rh,ib
        REAL qext(3,8,17,20),value

! Setting for aerosol with no growth (completely insoluble):
        IF(AEROTYPE.EQ.3) rh = 1

! Qext for Ammonium Sulfate (with a 90% insoluble dust inclusion):
        DATA(qext(1,1,ib,1),ib=1,17)  & ! Band 1, RH = 80%
        /0.00149, 0.00242, 0.00397, 0.00702, 0.01392, 0.04336  &
        ,0.16334, 0.54998, 1.39850, 2.52940, 2.77530, 2.71620  &
        ,2.41690, 2.20950, 2.17190, 2.12740, 2.09270/

        DATA(qext(1,1,ib,2),ib=1,17)  & ! Band 1, RH = 81%
        /0.00150, 0.00243, 0.00401, 0.00708, 0.01406, 0.04380  &
        ,0.16471, 0.55391, 1.40580, 2.53920, 2.78010, 2.71520  &
        ,2.41230, 2.21190, 2.17210, 2.12720, 2.09250/

        DATA(qext(1,1,ib,3),ib=1,17)  & ! Band 1, RH = 82%
        /0.00151, 0.00245, 0.00404, 0.00714, 0.01421, 0.04428  &
        ,0.16621, 0.55824, 1.41380, 2.55000, 2.78520, 2.71450  &
        ,2.40880, 2.21310, 2.17130, 2.12680, 2.09230/

        DATA(qext(1,1,ib,4),ib=1,17)  & ! Band 1, RH = 83%
        /0.00153, 0.00247, 0.00408, 0.00722, 0.01438, 0.04482  &
        ,0.16789, 0.56306, 1.42280, 2.56170, 2.79030, 2.71480  &
        ,2.40350, 2.21140, 2.16910, 2.12630, 2.09200/

        DATA(qext(1,1,ib,5),ib=1,17)  & ! Band 1, RH = 84%
        /0.00154, 0.00250, 0.00412, 0.00730, 0.01457, 0.04543  &
        ,0.16976, 0.56845, 1.43290, 2.57470, 2.79550, 2.71600  &
        ,2.39450, 2.20870, 2.16830, 2.12590, 2.09170/

        DATA(qext(1,1,ib,6),ib=1,17)  & ! Band 1, RH = 85%
        /0.00155, 0.00252, 0.00417, 0.00739, 0.01478, 0.04611  &
        ,0.17188, 0.57450, 1.44430, 2.58890, 2.80060, 2.71690  &
        ,2.38340, 2.21340, 2.16940, 2.12460, 2.09120/

        DATA(qext(1,1,ib,7),ib=1,17)  & ! Band 1, RH = 86%
        /0.00157, 0.00255, 0.00422, 0.00749, 0.01501, 0.04688  &
        ,0.17428, 0.58136, 1.45740, 2.60470, 2.80550, 2.71630  &
        ,2.37520, 2.22400, 2.16840, 2.12420, 2.09060/

        DATA(qext(1,1,ib,8),ib=1,17)  & ! Band 1, RH = 87%
        /0.00158, 0.00258, 0.00427, 0.00760, 0.01528, 0.04776  &
        ,0.17702, 0.58920, 1.47230, 2.62230, 2.81030, 2.71180  &
        ,2.36970, 2.22790, 2.16710, 2.12250, 2.08990/

        DATA(qext(1,1,ib,9),ib=1,17)  & ! Band 1, RH = 88%
        /0.00160, 0.00261, 0.00434, 0.00772, 0.01559, 0.04878  &
        ,0.18021, 0.59823, 1.48980, 2.64210, 2.81470, 2.70850  &
        ,2.35840, 2.22120, 2.16830, 2.12150, 2.08920/

        DATA(qext(1,1,ib,10),ib=1,17)  & ! Band 1, RH = 89%
        /0.00162, 0.00265, 0.00441, 0.00787, 0.01594, 0.04998  &
        ,0.18393, 0.60878, 1.51030, 2.66460, 2.81890, 2.70760  &
        ,2.34990, 2.22930, 2.16630, 2.11990, 2.08840/

        DATA(qext(1,1,ib,11),ib=1,17)  & ! Band 1, RH = 90%
        /0.00164, 0.00270, 0.00450, 0.00804, 0.01636, 0.05140  &
        ,0.18836, 0.62123, 1.53480, 2.69040, 2.82340, 2.69900  &
        ,2.34020, 2.23700, 2.16640, 2.11830, 2.08770/

        DATA(qext(1,1,ib,12),ib=1,17)  & ! Band 1, RH = 91%
        /0.00167, 0.00275, 0.00459, 0.00824, 0.01686, 0.05312  &
        ,0.19370, 0.63617, 1.56440, 2.72060, 2.82970, 2.68780  &
        ,2.32630, 2.22930, 2.16360, 2.11730, 2.08720/

        DATA(qext(1,1,ib,13),ib=1,17)  & ! Band 1, RH = 92%
        /0.00170, 0.00281, 0.00471, 0.00847, 0.01746, 0.05524  &
        ,0.20029, 0.65445, 1.60100, 2.75700, 2.83770, 2.68130  &
        ,2.31590, 2.23730, 2.16100, 2.11630, 2.08670/

        DATA(qext(1,1,ib,14),ib=1,17)  & ! Band 1, RH = 93%
        /0.00173, 0.00287, 0.00484, 0.00876, 0.01821, 0.05793  &
        ,0.20863, 0.67733, 1.64690, 2.80240, 2.84230, 2.66420  &
        ,2.30400, 2.24000, 2.15680, 2.11670, 2.08580/

        DATA(qext(1,1,ib,15),ib=1,17)  & ! Band 1, RH = 94%
        /0.00177, 0.00296, 0.00501, 0.00912, 0.01917, 0.06144  &
        ,0.21952, 0.70689, 1.70590, 2.86090, 2.83750, 2.63840  &
        ,2.29430, 2.23600, 2.15490, 2.11590, 2.08410/

        DATA(qext(1,1,ib,16),ib=1,17)  & ! Band 1, RH = 95%
        /0.00181, 0.00305, 0.00521, 0.00958, 0.02044, 0.06625  &
        ,0.23443, 0.74677, 1.78360, 2.93450, 2.81970, 2.60950  &
        ,2.28910, 2.22880, 2.15040, 2.11410, 2.08210/

        DATA(qext(1,1,ib,17),ib=1,17)  & ! Band 1, RH = 96%
        /0.00187, 0.00318, 0.00548, 0.01020, 0.02224, 0.07325  &
        ,0.25619, 0.80417, 1.88990, 3.01840, 2.79740, 2.56840  &
        ,2.27980, 2.21950, 2.15160, 2.10840, 2.08070/

        DATA(qext(1,1,ib,18),ib=1,17)  & ! Band 1, RH = 97%
        /0.00194, 0.00334, 0.00584, 0.01109, 0.02502, 0.08445  &
        ,0.29133, 0.89605, 2.04820, 3.11820, 2.73070, 2.49010  &
        ,2.26510, 2.19690, 2.14400, 2.10610, 2.07700/

        DATA(qext(1,1,ib,19),ib=1,17)  & ! Band 1, RH = 98%
        /0.00202, 0.00356, 0.00638, 0.01253, 0.02998, 0.10540  &
        ,0.35938, 1.07060, 2.32480, 3.22630, 2.61070, 2.36850  &
        ,2.24080, 2.18350, 2.13650, 2.09820, 2.07280/

        DATA(qext(1,1,ib,20),ib=1,17)  & ! Band 1, RH = 99%
        /0.00214, 0.00390, 0.00730, 0.01549, 0.04218, 0.15907  &
        ,0.54778, 1.48410, 2.85730, 3.16580, 2.46260, 2.24560  &
        ,2.24020, 2.15030, 2.11870, 2.08630, 2.06390/

        DATA(qext(1,2,ib,1),ib=1,17)  & ! Band 2, RH = 80%
        /0.00188, 0.00355, 0.00926, 0.03712, 0.17250, 0.70887  &
        ,2.02340, 3.42880, 3.06140, 2.60440, 2.33320, 2.26490  &
        ,2.18920, 2.14060, 2.10170, 2.07490, 2.05500/

        DATA(qext(1,2,ib,2),ib=1,17)  & ! Band 2, RH = 81%
        /0.00187, 0.00355, 0.00927, 0.03727, 0.17317, 0.71054  &
        ,2.02350, 3.42780, 3.06090, 2.60320, 2.33340, 2.26490  &
        ,2.18890, 2.14040, 2.10140, 2.07480, 2.05490/

        DATA(qext(1,2,ib,3),ib=1,17)  & ! Band 2, RH = 82%
        /0.00187, 0.00354, 0.00928, 0.03745, 0.17391, 0.71241  &
        ,2.02380, 3.42680, 3.06020, 2.60170, 2.33270, 2.26420  &
        ,2.18820, 2.13990, 2.10130, 2.07470, 2.05480/

        DATA(qext(1,2,ib,4),ib=1,17)  & ! Band 2, RH = 83%
        /0.00186, 0.00353, 0.00929, 0.03764, 0.17474, 0.71450  &
        ,2.02410, 3.42580, 3.05930, 2.59990, 2.33180, 2.26250  &
        ,2.18720, 2.13950, 2.10110, 2.07430, 2.05460/

        DATA(qext(1,2,ib,5),ib=1,17)  & ! Band 2, RH = 84%
        /0.00185, 0.00352, 0.00931, 0.03785, 0.17567, 0.71686  &
        ,2.02460, 3.42490, 3.05800, 2.59810, 2.33210, 2.26150  &
        ,2.18680, 2.13900, 2.10070, 2.07430, 2.05450/

        DATA(qext(1,2,ib,6),ib=1,17)  & ! Band 2, RH = 85%
        /0.00184, 0.00351, 0.00932, 0.03810, 0.17672, 0.71955  &
        ,2.02530, 3.42410, 3.05620, 2.59560, 2.33120, 2.26080  &
        ,2.18620, 2.13880, 2.10050, 2.07390, 2.05430/

        DATA(qext(1,2,ib,7),ib=1,17)  & ! Band 2, RH = 86%
        /0.00183, 0.00349, 0.00934, 0.03838, 0.17791, 0.72262  &
        ,2.02620, 3.42360, 3.05390, 2.59220, 2.33080, 2.25870  &
        ,2.18510, 2.13780, 2.10020, 2.07370, 2.05410/

        DATA(qext(1,2,ib,8),ib=1,17)  & ! Band 2, RH = 87%
        /0.00182, 0.00348, 0.00936, 0.03870, 0.17928, 0.72617  &
        ,2.02740, 3.42350, 3.05080, 2.58880, 2.32980, 2.25760  &
        ,2.18450, 2.13770, 2.09980, 2.07350, 2.05380/

        DATA(qext(1,2,ib,9),ib=1,17)  & ! Band 2, RH = 88%
        /0.00181, 0.00347, 0.00939, 0.03908, 0.18086, 0.73033  &
        ,2.02910, 3.42380, 3.04690, 2.58440, 2.33040, 2.25510  &
        ,2.18280, 2.13650, 2.09940, 2.07300, 2.05360/

        DATA(qext(1,2,ib,10),ib=1,17)  & ! Band 2, RH = 89%
        /0.00180, 0.00345, 0.00942, 0.03952, 0.18272, 0.73524  &
        ,2.03140, 3.42500, 3.04210, 2.58090, 2.33060, 2.25290  &
        ,2.18150, 2.13590, 2.09920, 2.07260, 2.05330/

        DATA(qext(1,2,ib,11),ib=1,17)  & ! Band 2, RH = 90%
        /0.00178, 0.00343, 0.00946, 0.04005, 0.18493, 0.74113  &
        ,2.03470, 3.42720, 3.03680, 2.57490, 2.33030, 2.25160  &
        ,2.18050, 2.13440, 2.09850, 2.07210, 2.05300/

        DATA(qext(1,2,ib,12),ib=1,17)  & ! Band 2, RH = 91%
        /0.00177, 0.00341, 0.00952, 0.04070, 0.18761, 0.74832  &
        ,2.03930, 3.43070, 3.03190, 2.56800, 2.32880, 2.24900  &
        ,2.17930, 2.13300, 2.09790, 2.07150, 2.05250/

        DATA(qext(1,2,ib,13),ib=1,17)  & ! Band 2, RH = 92%
        /0.00175, 0.00339, 0.00958, 0.04150, 0.19092, 0.75725  &
        ,2.04600, 3.43560, 3.02700, 2.56100, 2.32750, 2.24590  &
        ,2.17700, 2.13130, 2.09730, 2.07070, 2.05200/

        DATA(qext(1,2,ib,14),ib=1,17)  & ! Band 2, RH = 93%
        /0.00173, 0.00337, 0.00967, 0.04253, 0.19510, 0.76865  &
        ,2.05610, 3.44170, 3.01910, 2.55060, 2.32680, 2.24070  &
        ,2.17520, 2.12960, 2.09560, 2.06980, 2.05140/

        DATA(qext(1,2,ib,15),ib=1,17)  & ! Band 2, RH = 94%
        /0.00171, 0.00334, 0.00979, 0.04388, 0.20056, 0.78361  &
        ,2.07180, 3.44860, 3.00600, 2.53720, 2.32560, 2.23450  &
        ,2.17290, 2.12670, 2.09400, 2.06870, 2.05060/

        DATA(qext(1,2,ib,16),ib=1,17)  & ! Band 2, RH = 95%
        /0.00168, 0.00331, 0.00997, 0.04574, 0.20801, 0.80407  &
        ,2.09750, 3.45520, 2.98110, 2.51670, 2.32710, 2.22770  &
        ,2.17010, 2.12290, 2.09200, 2.06740, 2.04970/

        DATA(qext(1,2,ib,17),ib=1,17)  & ! Band 2, RH = 96%
        /0.00165, 0.00328, 0.01024, 0.04847, 0.21878, 0.83363  &
        ,2.14150, 3.46280, 2.94540, 2.49200, 2.32280, 2.21880  &
        ,2.16770, 2.12000, 2.08910, 2.06590, 2.04830/

        DATA(qext(1,2,ib,18),ib=1,17)  & ! Band 2, RH = 97%
        /0.00161, 0.00325, 0.01068, 0.05282, 0.23577, 0.88021  &
        ,2.21860, 3.48830, 2.88970, 2.45290, 2.31840, 2.21050  &
        ,2.16160, 2.11590, 2.08590, 2.06310, 2.04650/

        DATA(qext(1,2,ib,19),ib=1,17)  & ! Band 2, RH = 98%
        /0.00156, 0.00322, 0.01150, 0.06079, 0.26694, 0.96671  &
        ,2.35780, 3.52060, 2.77450, 2.37700, 2.31220, 2.19790  &
        ,2.15500, 2.11260, 2.08080, 2.05910, 2.04360/

        DATA(qext(1,2,ib,20),ib=1,17)  & ! Band 2, RH = 99%
        /0.00150, 0.00323, 0.01343, 0.07997, 0.34651, 1.20200  &
        ,2.70630, 3.51110, 2.55470, 2.24690, 2.28820, 2.19750  &
        ,2.12670, 2.09840, 2.07060, 2.05250, 2.03840/

        DATA(qext(1,3,ib,1),ib=1,17)  & ! Band 3, RH = 80%
        /0.00551, 0.01885, 0.08808, 0.39024, 1.29860, 2.89670  &
        ,3.70940, 2.63090, 2.41080, 2.25560, 2.20910, 2.15790  &
        ,2.11870, 2.08870, 2.06480, 2.04770, 2.03510/

        DATA(qext(1,3,ib,2),ib=1,17)  & ! Band 3, RH = 81%
        /0.00551, 0.01891, 0.08848, 0.39154, 1.30010, 2.89550  &
        ,3.70770, 2.63050, 2.41050, 2.25380, 2.20850, 2.15760  &
        ,2.11860, 2.08850, 2.06470, 2.04760, 2.03500/

        DATA(qext(1,3,ib,3),ib=1,17)  & ! Band 3, RH = 82%
        /0.00551, 0.01898, 0.08893, 0.39298, 1.30180, 2.89420  &
        ,3.70590, 2.62920, 2.41030, 2.25390, 2.20810, 2.15700  &
        ,2.11840, 2.08830, 2.06460, 2.04750, 2.03500/

        DATA(qext(1,3,ib,4),ib=1,17)  & ! Band 3, RH = 83%
        /0.00551, 0.01905, 0.08943, 0.39457, 1.30380, 2.89300  &
        ,3.70390, 2.62760, 2.40910, 2.25360, 2.20740, 2.15680  &
        ,2.11800, 2.08810, 2.06440, 2.04740, 2.03490/

        DATA(qext(1,3,ib,5),ib=1,17)  & ! Band 3, RH = 84%
        /0.00551, 0.01914, 0.08998, 0.39636, 1.30600, 2.89170  &
        ,3.70190, 2.62610, 2.40740, 2.25150, 2.20660, 2.15630  &
        ,2.11750, 2.08800, 2.06420, 2.04730, 2.03480/

        DATA(qext(1,3,ib,6),ib=1,17)  & ! Band 3, RH = 85%
        /0.00551, 0.01923, 0.09060, 0.39837, 1.30850, 2.89040  &
        ,3.70000, 2.62440, 2.40800, 2.25050, 2.20620, 2.15550  &
        ,2.11720, 2.08760, 2.06400, 2.04710, 2.03460/

        DATA(qext(1,3,ib,7),ib=1,17)  & ! Band 3, RH = 86%
        /0.00551, 0.01933, 0.09131, 0.40064, 1.31150, 2.88920  &
        ,3.69760, 2.62270, 2.40740, 2.24890, 2.20580, 2.15510  &
        ,2.11670, 2.08730, 2.06370, 2.04690, 2.03450/

        DATA(qext(1,3,ib,8),ib=1,17)  & ! Band 3, RH = 87%
        /0.00551, 0.01945, 0.09211, 0.40325, 1.31500, 2.88800  &
        ,3.69400, 2.61990, 2.40630, 2.24880, 2.20490, 2.15420  &
        ,2.11610, 2.08680, 2.06340, 2.04670, 2.03430/

        DATA(qext(1,3,ib,9),ib=1,17)  & ! Band 3, RH = 88%
        /0.00551, 0.01959, 0.09302, 0.40626, 1.31910, 2.88700  &
        ,3.68970, 2.61740, 2.40600, 2.24650, 2.20380, 2.15380  &
        ,2.11560, 2.08660, 2.06310, 2.04650, 2.03420/

        DATA(qext(1,3,ib,10),ib=1,17)  & ! Band 3, RH = 89%
        /0.00551, 0.01975, 0.09409, 0.40979, 1.32400, 2.88630  &
        ,3.68580, 2.61230, 2.40590, 2.24480, 2.20240, 2.15230  &
        ,2.11490, 2.08610, 2.06270, 2.04620, 2.03400/

        DATA(qext(1,3,ib,11),ib=1,17)  & ! Band 3, RH = 90%
        /0.00551, 0.01993, 0.09535, 0.41396, 1.32990, 2.88600  &
        ,3.68140, 2.60580, 2.40550, 2.24300, 2.20110, 2.15100  &
        ,2.11390, 2.08550, 2.06240, 2.04590, 2.03380/

        DATA(qext(1,3,ib,12),ib=1,17)  & ! Band 3, RH = 91%
        /0.00551, 0.02015, 0.09686, 0.41900, 1.33730, 2.88670  &
        ,3.67640, 2.59860, 2.40280, 2.24230, 2.20020, 2.15010  &
        ,2.11310, 2.08490, 2.06180, 2.04550, 2.03350/

        DATA(qext(1,3,ib,13),ib=1,17)  & ! Band 3, RH = 92%
        /0.00552, 0.02042, 0.09869, 0.42519, 1.34650, 2.88870  &
        ,3.67070, 2.58810, 2.39860, 2.23810, 2.19920, 2.14890  &
        ,2.11150, 2.08400, 2.06120, 2.04510, 2.03320/

        DATA(qext(1,3,ib,14),ib=1,17)  & ! Band 3, RH = 93%
        /0.00552, 0.02076, 0.10097, 0.43299, 1.35840, 2.89280  &
        ,3.66370, 2.57630, 2.39500, 2.23540, 2.19610, 2.14750  &
        ,2.11050, 2.08340, 2.06040, 2.04460, 2.03280/

        DATA(qext(1,3,ib,15),ib=1,17)  & ! Band 3, RH = 94%
        /0.00554, 0.02118, 0.10388, 0.44312, 1.37410, 2.90010  &
        ,3.65460, 2.55910, 2.39170, 2.23010, 2.19550, 2.14500  &
        ,2.10840, 2.08210, 2.05950, 2.04390, 2.03230/

        DATA(qext(1,3,ib,16),ib=1,17)  & ! Band 3, RH = 95%
        /0.00555, 0.02174, 0.10775, 0.45682, 1.39590, 2.91340  &
        ,3.64110, 2.53880, 2.39300, 2.22590, 2.19220, 2.14330  &
        ,2.10700, 2.08030, 2.05830, 2.04300, 2.03160/

        DATA(qext(1,3,ib,17),ib=1,17)  & ! Band 3, RH = 96%
        /0.00558, 0.02250, 0.11313, 0.47637, 1.42780, 2.93820  &
        ,3.61530, 2.51250, 2.39590, 2.21830, 2.19110, 2.14040  &
        ,2.10350, 2.07690, 2.05710, 2.04190, 2.03080/

        DATA(qext(1,3,ib,18),ib=1,17)  & ! Band 3, RH = 97%
        /0.00562, 0.02360, 0.12118, 0.50645, 1.48020, 2.98580  &
        ,3.57030, 2.45720, 2.39760, 2.21380, 2.19370, 2.13760  &
        ,2.09990, 2.07450, 2.05490, 2.04030, 2.02960/

        DATA(qext(1,3,ib,19),ib=1,17)  & ! Band 3, RH = 98%
        /0.00569, 0.02534, 0.13461, 0.55830, 1.58210, 3.08140  &
        ,3.49730, 2.36320, 2.40250, 2.21830, 2.19310, 2.13570  &
        ,2.09610, 2.06970, 2.05130, 2.03780, 2.02780/

        DATA(qext(1,3,ib,20),ib=1,17)  & ! Band 3, RH = 99%
        /0.00582, 0.02851, 0.16184, 0.66856, 1.82450, 3.29450  &
        ,3.25540, 2.21460, 2.38670, 2.24740, 2.15090, 2.11170  &
        ,2.08330, 2.06190, 2.04540, 2.03320, 2.02430/

        DATA(qext(1,4,ib,1),ib=1,17)  & ! Band 4, RH = 80%
        /0.00254, 0.00407, 0.00650, 0.01036, 0.01535, 0.02445  &
        ,0.03903, 0.06275, 0.10298, 0.17894, 0.35874, 0.90169  &
        ,2.24330, 3.18390, 2.89810, 2.58850, 2.44120/

        DATA(qext(1,4,ib,2),ib=1,17)  & ! Band 4, RH = 81%
        /0.00256, 0.00410, 0.00654, 0.01044, 0.01547, 0.02465  &
        ,0.03935, 0.06328, 0.10386, 0.18051, 0.36192, 0.90814  &
        ,2.24520, 3.17780, 2.89510, 2.58700, 2.44020/

        DATA(qext(1,4,ib,3),ib=1,17)  & ! Band 4, RH = 82%
        /0.00257, 0.00412, 0.00659, 0.01052, 0.01561, 0.02487  &
        ,0.03971, 0.06386, 0.10483, 0.18224, 0.36544, 0.91524  &
        ,2.24730, 3.17110, 2.89170, 2.58540, 2.43900/

        DATA(qext(1,4,ib,4),ib=1,17)  & ! Band 4, RH = 83%
        /0.00259, 0.00416, 0.00665, 0.01061, 0.01575, 0.02511  &
        ,0.04010, 0.06450, 0.10590, 0.18416, 0.36934, 0.92309  &
        ,2.24960, 3.16390, 2.88800, 2.58360, 2.43780/

        DATA(qext(1,4,ib,5),ib=1,17)  & ! Band 4, RH = 84%
        /0.00261, 0.00419, 0.00671, 0.01071, 0.01592, 0.02538  &
        ,0.04054, 0.06521, 0.10710, 0.18630, 0.37370, 0.93181  &
        ,2.25200, 3.15600, 2.88390, 2.58160, 2.43630/

        DATA(qext(1,4,ib,6),ib=1,17)  & ! Band 4, RH = 85%
        /0.00263, 0.00423, 0.00677, 0.01082, 0.01610, 0.02568  &
        ,0.04103, 0.06601, 0.10843, 0.18870, 0.37859, 0.94156  &
        ,2.25480, 3.14750, 2.87940, 2.57940, 2.43480/

        DATA(qext(1,4,ib,7),ib=1,17)  & ! Band 4, RH = 86%
        /0.00265, 0.00427, 0.00685, 0.01095, 0.01631, 0.02602  &
        ,0.04158, 0.06691, 0.10994, 0.19140, 0.38412, 0.95255  &
        ,2.25780, 3.13810, 2.87430, 2.57700, 2.43300/

        DATA(qext(1,4,ib,8),ib=1,17)  & ! Band 4, RH = 87%
        /0.00267, 0.00431, 0.00693, 0.01109, 0.01654, 0.02640  &
        ,0.04220, 0.06794, 0.11166, 0.19449, 0.39043, 0.96501  &
        ,2.26130, 3.12790, 2.86860, 2.57420, 2.43100/

        DATA(qext(1,4,ib,9),ib=1,17)  & ! Band 4, RH = 88%
        /0.00270, 0.00436, 0.00702, 0.01125, 0.01681, 0.02684  &
        ,0.04292, 0.06911, 0.11363, 0.19804, 0.39770, 0.97927  &
        ,2.26520, 3.11660, 2.86210, 2.57110, 2.42870/

        DATA(qext(1,4,ib,10),ib=1,17)  & ! Band 4, RH = 89%
        /0.00273, 0.00442, 0.00713, 0.01143, 0.01711, 0.02734  &
        ,0.04375, 0.07046, 0.11591, 0.20216, 0.40616, 0.99575  &
        ,2.26960, 3.10410, 2.85470, 2.56750, 2.42600/

        DATA(qext(1,4,ib,11),ib=1,17)  & ! Band 4, RH = 90%
        /0.00276, 0.00449, 0.00725, 0.01164, 0.01747, 0.02793  &
        ,0.04471, 0.07205, 0.11859, 0.20701, 0.41614, 1.01500  &
        ,2.27490, 3.09030, 2.84620, 2.56340, 2.42290/

        DATA(qext(1,4,ib,12),ib=1,17)  & ! Band 4, RH = 91%
        /0.00280, 0.00457, 0.00739, 0.01189, 0.01789, 0.02863  &
        ,0.04586, 0.07394, 0.12178, 0.21280, 0.42809, 1.03790  &
        ,2.28110, 3.07500, 2.83630, 2.55850, 2.41930/

        DATA(qext(1,4,ib,13),ib=1,17)  & ! Band 4, RH = 92%
        /0.00285, 0.00465, 0.00756, 0.01219, 0.01839, 0.02947  &
        ,0.04724, 0.07622, 0.12565, 0.21985, 0.44269, 1.06540  &
        ,2.28870, 3.05790, 2.82450, 2.55260, 2.41490/

        DATA(qext(1,4,ib,14),ib=1,17)  & ! Band 4, RH = 93%
        /0.00290, 0.00476, 0.00776, 0.01255, 0.01901, 0.03050  &
        ,0.04894, 0.07904, 0.13044, 0.22862, 0.46091, 1.09920  &
        ,2.29840, 3.03870, 2.81020, 2.54540, 2.40960/

        DATA(qext(1,4,ib,15),ib=1,17)  & ! Band 4, RH = 94%
        /0.00296, 0.00489, 0.00801, 0.01300, 0.01978, 0.03179  &
        ,0.05109, 0.08261, 0.13653, 0.23985, 0.48434, 1.14180  &
        ,2.31120, 3.01700, 2.79240, 2.53650, 2.40290/

        DATA(qext(1,4,ib,16),ib=1,17)  & ! Band 4, RH = 95%
        /0.00303, 0.00504, 0.00832, 0.01357, 0.02077, 0.03348  &
        ,0.05390, 0.08730, 0.14457, 0.25478, 0.51567, 1.19720  &
        ,2.32910, 2.99200, 2.76940, 2.52510, 2.39430/

        DATA(qext(1,4,ib,17),ib=1,17)  & ! Band 4, RH = 96%
        /0.00312, 0.00524, 0.00872, 0.01434, 0.02211, 0.03578  &
        ,0.05775, 0.09375, 0.15572, 0.27574, 0.55989, 1.27230  &
        ,2.35630, 2.96150, 2.73870, 2.50990, 2.38290/

        DATA(qext(1,4,ib,18),ib=1,17)  & ! Band 4, RH = 97%
        /0.00323, 0.00551, 0.00928, 0.01541, 0.02404, 0.03911  &
        ,0.06340, 0.10331, 0.17240, 0.30764, 0.62751, 1.38080  &
        ,2.40100, 2.92120, 2.69690, 2.48850, 2.36690/

        DATA(qext(1,4,ib,19),ib=1,17)  & ! Band 4, RH = 98%
        /0.00337, 0.00587, 0.01010, 0.01708, 0.02710, 0.04453  &
        ,0.07271, 0.11925, 0.20078, 0.36342, 0.74568, 1.55280  &
        ,2.47950, 2.86420, 2.63600, 2.45590, 2.34270/

        DATA(qext(1,4,ib,20),ib=1,17)  & ! Band 4, RH = 99%
        /0.00357, 0.00644, 0.01151, 0.02016, 0.03309, 0.05554  &
        ,0.09218, 0.15362, 0.26458, 0.49546, 1.01500, 1.87470  &
        ,2.61310, 2.76300, 2.53870, 2.39880, 2.30010/

        DATA(qext(1,5,ib,1),ib=1,17)  & ! Band 5, RH = 80%
        /0.00242, 0.00394, 0.00635, 0.01020, 0.01324, 0.02120  &
        ,0.03410, 0.05576, 0.09577, 0.18722, 0.45630, 1.20090  &
        ,2.49710, 3.34910, 2.66470, 2.46510, 2.32040/

        DATA(qext(1,5,ib,2),ib=1,17)  & ! Band 5, RH = 81%
        /0.00246, 0.00401, 0.00647, 0.01040, 0.01357, 0.02174  &
        ,0.03498, 0.05718, 0.09817, 0.19146, 0.46379, 1.20960  &
        ,2.49280, 3.33110, 2.66380, 2.46090, 2.31910/

        DATA(qext(1,5,ib,3),ib=1,17)  & ! Band 5, RH = 82%
        /0.00250, 0.00408, 0.00660, 0.01062, 0.01393, 0.02233  &
        ,0.03594, 0.05875, 0.10080, 0.19614, 0.47202, 1.21900  &
        ,2.48820, 3.31180, 2.66260, 2.45650, 2.31770/

        DATA(qext(1,5,ib,4),ib=1,17)  & ! Band 5, RH = 83%
        /0.00255, 0.00416, 0.00675, 0.01086, 0.01433, 0.02298  &
        ,0.03700, 0.06049, 0.10372, 0.20131, 0.48109, 1.22930  &
        ,2.48340, 3.29090, 2.66090, 2.45200, 2.31630/

        DATA(qext(1,5,ib,5),ib=1,17)  & ! Band 5, RH = 84%
        /0.00260, 0.00425, 0.00691, 0.01113, 0.01477, 0.02370  &
        ,0.03817, 0.06242, 0.10696, 0.20705, 0.49114, 1.24070  &
        ,2.47820, 3.26840, 2.65870, 2.44730, 2.31480/

        DATA(qext(1,5,ib,6),ib=1,17)  & ! Band 5, RH = 85%
        /0.00265, 0.00435, 0.00708, 0.01143, 0.01526, 0.02451  &
        ,0.03948, 0.06457, 0.11059, 0.21348, 0.50236, 1.25320  &
        ,2.47260, 3.24390, 2.65590, 2.44240, 2.31310/

        DATA(qext(1,5,ib,7),ib=1,17)  & ! Band 5, RH = 86%
        /0.00271, 0.00446, 0.00727, 0.01176, 0.01581, 0.02541  &
        ,0.04096, 0.06700, 0.11468, 0.22073, 0.51494, 1.26710  &
        ,2.46670, 3.21730, 2.65230, 2.43720, 2.31130/

        DATA(qext(1,5,ib,8),ib=1,17)  & ! Band 5, RH = 87%
        /0.00277, 0.00458, 0.00749, 0.01213, 0.01643, 0.02644  &
        ,0.04264, 0.06975, 0.11932, 0.22895, 0.52915, 1.28270  &
        ,2.46030, 3.18810, 2.64780, 2.43180, 2.30920/

        DATA(qext(1,5,ib,9),ib=1,17)  & ! Band 5, RH = 88%
        /0.00284, 0.00472, 0.00774, 0.01256, 0.01714, 0.02761  &
        ,0.04455, 0.07289, 0.12463, 0.23837, 0.54535, 1.30010  &
        ,2.45340, 3.15620, 2.64210, 2.42610, 2.30690/

        DATA(qext(1,5,ib,10),ib=1,17)  & ! Band 5, RH = 89%
        /0.00292, 0.00488, 0.00802, 0.01305, 0.01796, 0.02895  &
        ,0.04676, 0.07653, 0.13077, 0.24926, 0.56397, 1.31990  &
        ,2.44600, 3.12090, 2.63500, 2.42000, 2.30430/

        DATA(qext(1,5,ib,11),ib=1,17)  & ! Band 5, RH = 90%
        /0.00301, 0.00505, 0.00834, 0.01361, 0.01891, 0.03052  &
        ,0.04933, 0.08078, 0.13796, 0.26200, 0.58561, 1.34250  &
        ,2.43800, 3.08180, 2.62620, 2.41350, 2.30130/

        DATA(qext(1,5,ib,12),ib=1,17)  & ! Band 5, RH = 91%
        /0.00311, 0.00526, 0.00872, 0.01427, 0.02002, 0.03238  &
        ,0.05238, 0.08581, 0.14649, 0.27713, 0.61107, 1.36850  &
        ,2.42950, 3.03820, 2.61510, 2.40630, 2.29760/

        DATA(qext(1,5,ib,13),ib=1,17)  & ! Band 5, RH = 92%
        /0.00323, 0.00549, 0.00916, 0.01505, 0.02135, 0.03460  &
        ,0.05604, 0.09188, 0.15678, 0.29537, 0.64147, 1.39880  &
        ,2.42040, 2.98930, 2.60130, 2.39830, 2.29320/

        DATA(qext(1,5,ib,14),ib=1,17)  & ! Band 5, RH = 93%
        /0.00337, 0.00577, 0.00969, 0.01600, 0.02297, 0.03731  &
        ,0.06053, 0.09933, 0.16944, 0.31784, 0.67841, 1.43460  &
        ,2.41080, 2.93430, 2.58380, 2.38920, 2.28780/

        DATA(qext(1,5,ib,15),ib=1,17)  & ! Band 5, RH = 94%
        /0.00352, 0.00610, 0.01033, 0.01716, 0.02499, 0.04071  &
        ,0.06618, 0.10872, 0.18545, 0.34620, 0.72430, 1.47770  &
        ,2.40100, 2.87210, 2.56160, 2.37850, 2.28100/

        DATA(qext(1,5,ib,16),ib=1,17)  & ! Band 5, RH = 95%
        /0.00371, 0.00651, 0.01114, 0.01865, 0.02759, 0.04510  &
        ,0.07350, 0.12094, 0.20636, 0.38324, 0.78291, 1.53050  &
        ,2.39160, 2.80170, 2.53300, 2.36540, 2.27230/

        DATA(qext(1,5,ib,17),ib=1,17)  & ! Band 5, RH = 96%
        /0.00394, 0.00703, 0.01219, 0.02061, 0.03105, 0.05102  &
        ,0.08342, 0.13760, 0.23498, 0.43381, 0.86060, 1.59740  &
        ,2.38360, 2.72230, 2.49650, 2.34870, 2.26100/

        DATA(qext(1,5,ib,18),ib=1,17)  & ! Band 5, RH = 97%
        /0.00423, 0.00771, 0.01361, 0.02334, 0.03596, 0.05950  &
        ,0.09777, 0.16184, 0.27693, 0.50762, 0.96915, 1.68600  &
        ,2.37880, 2.63290, 2.45060, 2.32640, 2.24580/

        DATA(qext(1,5,ib,19),ib=1,17)  & ! Band 5, RH = 98%
        /0.00460, 0.00864, 0.01569, 0.02749, 0.04362, 0.07297  &
        ,0.12084, 0.20129, 0.34589, 0.62753, 1.13410, 1.81250  &
        ,2.37980, 2.53210, 2.39340, 2.29560, 2.22470/

        DATA(qext(1,5,ib,20),ib=1,17)  & ! Band 5, RH = 99%
        /0.00511, 0.01006, 0.01915, 0.03493, 0.05805, 0.09927  &
        ,0.16713, 0.28246, 0.49044, 0.86851, 1.42940, 2.02170  &
        ,2.39120, 2.41690, 2.32250, 2.25050, 2.19270/

        DATA(qext(1,6,ib,1),ib=1,17)  & ! Band 6, RH = 80%
        /0.00499, 0.00790, 0.01251, 0.01984, 0.03252, 0.05194  &
        ,0.08390, 0.13958, 0.24905, 0.49482, 0.98441, 1.69740  &
        ,2.35350, 2.59790, 2.44710, 2.28060, 2.21050/

        DATA(qext(1,6,ib,2),ib=1,17)  & ! Band 6, RH = 81%
        /0.00498, 0.00787, 0.01247, 0.01978, 0.03242, 0.05176  &
        ,0.08362, 0.13916, 0.24846, 0.49408, 0.98344, 1.69790  &
        ,2.35730, 2.60170, 2.44690, 2.28000, 2.21020/

        DATA(qext(1,6,ib,3),ib=1,17)  & ! Band 6, RH = 82%
        /0.00497, 0.00785, 0.01243, 0.01971, 0.03230, 0.05158  &
        ,0.08332, 0.13870, 0.24782, 0.49327, 0.98243, 1.69860  &
        ,2.36140, 2.60580, 2.44650, 2.27940, 2.20980/

        DATA(qext(1,6,ib,4),ib=1,17)  & ! Band 6, RH = 83%
        /0.00495, 0.00783, 0.01239, 0.01964, 0.03217, 0.05137  &
        ,0.08299, 0.13819, 0.24713, 0.49240, 0.98136, 1.69950  &
        ,2.36610, 2.61030, 2.44600, 2.27870, 2.20940/

        DATA(qext(1,6,ib,5),ib=1,17)  & ! Band 6, RH = 84%
        /0.00494, 0.00780, 0.01234, 0.01956, 0.03203, 0.05114  &
        ,0.08263, 0.13764, 0.24637, 0.49147, 0.98025, 1.70050  &
        ,2.37140, 2.61540, 2.44530, 2.27800, 2.20900/

        DATA(qext(1,6,ib,6),ib=1,17)  & ! Band 6, RH = 85%
        /0.00492, 0.00777, 0.01229, 0.01947, 0.03188, 0.05089  &
        ,0.08223, 0.13703, 0.24554, 0.49046, 0.97909, 1.70170  &
        ,2.37740, 2.62090, 2.44440, 2.27720, 2.20850/

        DATA(qext(1,6,ib,7),ib=1,17)  & ! Band 6, RH = 86%
        /0.00490, 0.00774, 0.01223, 0.01937, 0.03170, 0.05061  &
        ,0.08179, 0.13635, 0.24462, 0.48937, 0.97790, 1.70320  &
        ,2.38420, 2.62710, 2.44330, 2.27630, 2.20790/

        DATA(qext(1,6,ib,8),ib=1,17)  & ! Band 6, RH = 87%
        /0.00488, 0.00770, 0.01217, 0.01926, 0.03151, 0.05030  &
        ,0.08129, 0.13560, 0.24361, 0.48818, 0.97670, 1.70520  &
        ,2.39210, 2.63410, 2.44180, 2.27530, 2.20730/

        DATA(qext(1,6,ib,9),ib=1,17)  & ! Band 6, RH = 88%
        /0.00486, 0.00766, 0.01209, 0.01914, 0.03130, 0.04995  &
        ,0.08073, 0.13475, 0.24248, 0.48691, 0.97551, 1.70760  &
        ,2.40130, 2.64190, 2.43990, 2.27430, 2.20650/

        DATA(qext(1,6,ib,10),ib=1,17)  & ! Band 6, RH = 89%
        /0.00484, 0.00762, 0.01201, 0.01901, 0.03105, 0.04955  &
        ,0.08010, 0.13380, 0.24123, 0.48553, 0.97440, 1.71080  &
        ,2.41220, 2.65080, 2.43740, 2.27320, 2.20560/

        DATA(qext(1,6,ib,11),ib=1,17)  & ! Band 6, RH = 90%
        /0.00481, 0.00757, 0.01192, 0.01885, 0.03078, 0.04910  &
        ,0.07938, 0.13271, 0.23982, 0.48406, 0.97344, 1.71490  &
        ,2.42520, 2.66090, 2.43410, 2.27210, 2.20450/

        DATA(qext(1,6,ib,12),ib=1,17)  & ! Band 6, RH = 91%
        /0.00478, 0.00751, 0.01182, 0.01867, 0.03046, 0.04858  &
        ,0.07855, 0.13147, 0.23825, 0.48251, 0.97279, 1.72040  &
        ,2.44110, 2.67260, 2.42970, 2.27110, 2.20320/

        DATA(qext(1,6,ib,13),ib=1,17)  & ! Band 6, RH = 92%
        /0.00475, 0.00744, 0.01170, 0.01847, 0.03009, 0.04797  &
        ,0.07759, 0.13003, 0.23647, 0.48093, 0.97270, 1.72790  &
        ,2.46090, 2.68610, 2.42330, 2.26990, 2.20150/

        DATA(qext(1,6,ib,14),ib=1,17)  & ! Band 6, RH = 93%
        /0.00471, 0.00737, 0.01157, 0.01823, 0.02965, 0.04726  &
        ,0.07645, 0.12835, 0.23449, 0.47939, 0.97367, 1.73840  &
        ,2.48600, 2.70200, 2.41390, 2.26910, 2.19960/

        DATA(qext(1,6,ib,15),ib=1,17)  & ! Band 6, RH = 94%
        /0.00467, 0.00728, 0.01140, 0.01794, 0.02913, 0.04640  &
        ,0.07510, 0.12638, 0.23230, 0.47812, 0.97656, 1.75390  &
        ,2.51900, 2.72070, 2.40030, 2.26880, 2.19720/

        DATA(qext(1,6,ib,16),ib=1,17)  & ! Band 6, RH = 95%
        /0.00462, 0.00718, 0.01121, 0.01760, 0.02850, 0.04537  &
        ,0.07347, 0.12405, 0.22996, 0.47757, 0.98317, 1.77750  &
        ,2.56350, 2.74130, 2.38190, 2.26880, 2.19440/

        DATA(qext(1,6,ib,17),ib=1,17)  & ! Band 6, RH = 96%
        /0.00456, 0.00705, 0.01097, 0.01718, 0.02772, 0.04409  &
        ,0.07147, 0.12131, 0.22775, 0.47887, 0.99733, 1.81550  &
        ,2.62500, 2.75890, 2.35230, 2.26910, 2.19110/

        DATA(qext(1,6,ib,18),ib=1,17)  & ! Band 6, RH = 97%
        /0.00449, 0.00690, 0.01068, 0.01665, 0.02674, 0.04248  &
        ,0.06902, 0.11822, 0.22660, 0.48500, 1.02820, 1.88150  &
        ,2.71410, 2.77030, 2.31190, 2.26620, 2.18620/

        DATA(qext(1,6,ib,19),ib=1,17)  & ! Band 6, RH = 98%
        /0.00440, 0.00671, 0.01030, 0.01598, 0.02548, 0.04049  &
        ,0.06614, 0.11546, 0.23018, 0.50585, 1.10060, 2.01490  &
        ,2.86900, 2.73110, 2.26860, 2.24470, 2.17470/

        DATA(qext(1,6,ib,20),ib=1,17)  & ! Band 6, RH = 99%
        /0.00429, 0.00646, 0.00981, 0.01515, 0.02399, 0.03839  &
        ,0.06412, 0.11859, 0.25934, 0.59425, 1.31360, 2.37590  &
        ,3.08790, 2.48080, 2.31080, 2.20770, 2.15760/

        DATA(qext(1,7,ib,1),ib=1,17)  & ! Band 7, RH = 80%
        /0.00480, 0.00761, 0.01208, 0.01921, 0.03078, 0.04945  &
        ,0.08113, 0.14118, 0.28491, 0.73426, 1.97470, 3.23210  &
        ,3.10530, 2.50640, 2.43460, 2.32000, 2.23810/

        DATA(qext(1,7,ib,2),ib=1,17)  & ! Band 7, RH = 81%
        /0.00479, 0.00761, 0.01208, 0.01919, 0.03077, 0.04942  &
        ,0.08109, 0.14112, 0.28480, 0.73326, 1.96640, 3.22440  &
        ,3.10860, 2.50420, 2.43450, 2.31940, 2.23750/

        DATA(qext(1,7,ib,3),ib=1,17)  & ! Band 7, RH = 82%
        /0.00479, 0.00760, 0.01207, 0.01918, 0.03075, 0.04940  &
        ,0.08105, 0.14105, 0.28468, 0.73217, 1.95740, 3.21590  &
        ,3.11220, 2.50190, 2.43440, 2.31870, 2.23690/

        DATA(qext(1,7,ib,4),ib=1,17)  & ! Band 7, RH = 83%
        /0.00479, 0.00760, 0.01206, 0.01916, 0.03073, 0.04936  &
        ,0.08100, 0.14098, 0.28454, 0.73094, 1.94760, 3.20650  &
        ,3.11610, 2.49950, 2.43430, 2.31790, 2.23630/

        DATA(qext(1,7,ib,5),ib=1,17)  & ! Band 7, RH = 84%
        /0.00479, 0.00759, 0.01205, 0.01914, 0.03071, 0.04933  &
        ,0.08093, 0.14088, 0.28438, 0.72958, 1.93680, 3.19600  &
        ,3.12030, 2.49690, 2.43420, 2.31710, 2.23550/

        DATA(qext(1,7,ib,6),ib=1,17)  & ! Band 7, RH = 85%
        /0.00478, 0.00758, 0.01203, 0.01912, 0.03068, 0.04928  &
        ,0.08086, 0.14077, 0.28419, 0.72803, 1.92490, 3.18420  &
        ,3.12500, 2.49410, 2.43410, 2.31610, 2.23470/

        DATA(qext(1,7,ib,7),ib=1,17)  & ! Band 7, RH = 86%
        /0.00478, 0.00758, 0.01202, 0.01910, 0.03065, 0.04923  &
        ,0.08078, 0.14064, 0.28396, 0.72629, 1.91170, 3.17090  &
        ,3.13000, 2.49110, 2.43390, 2.31510, 2.23380/

        DATA(qext(1,7,ib,8),ib=1,17)  & ! Band 7, RH = 87%
        /0.00477, 0.00757, 0.01200, 0.01907, 0.03061, 0.04916  &
        ,0.08067, 0.14048, 0.28369, 0.72429, 1.89710, 3.15580  &
        ,3.13550, 2.48800, 2.43380, 2.31400, 2.23270/

        DATA(qext(1,7,ib,9),ib=1,17)  & ! Band 7, RH = 88%
        /0.00477, 0.00756, 0.01198, 0.01903, 0.03056, 0.04908  &
        ,0.08055, 0.14029, 0.28336, 0.72199, 1.88060, 3.13850  &
        ,3.14160, 2.48470, 2.43360, 2.31270, 2.23160/

        DATA(qext(1,7,ib,10),ib=1,17)  & ! Band 7, RH = 89%
        /0.00476, 0.00754, 0.01196, 0.01899, 0.03050, 0.04899  &
        ,0.08039, 0.14006, 0.28296, 0.71932, 1.86210, 3.11850  &
        ,3.14820, 2.48130, 2.43340, 2.31130, 2.23020/

        DATA(qext(1,7,ib,11),ib=1,17)  & ! Band 7, RH = 90%
        /0.00476, 0.00753, 0.01193, 0.01894, 0.03043, 0.04887  &
        ,0.08020, 0.13976, 0.28247, 0.71619, 1.84100, 3.09520  &
        ,3.15540, 2.47790, 2.43300, 2.30960, 2.22860/

        DATA(qext(1,7,ib,12),ib=1,17)  & ! Band 7, RH = 91%
        /0.00475, 0.00751, 0.01190, 0.01888, 0.03034, 0.04872  &
        ,0.07996, 0.13939, 0.28187, 0.71247, 1.81690, 3.06760  &
        ,3.16320, 2.47450, 2.43250, 2.30760, 2.22670/

        DATA(qext(1,7,ib,13),ib=1,17)  & ! Band 7, RH = 92%
        /0.00474, 0.00749, 0.01186, 0.01881, 0.03023, 0.04853  &
        ,0.07965, 0.13892, 0.28111, 0.70800, 1.78910, 3.03460  &
        ,3.17130, 2.47150, 2.43170, 2.30530, 2.22450/

        DATA(qext(1,7,ib,14),ib=1,17)  & ! Band 7, RH = 93%
        /0.00473, 0.00747, 0.01181, 0.01872, 0.03008, 0.04828  &
        ,0.07925, 0.13832, 0.28015, 0.70256, 1.75670, 2.99450  &
        ,3.17930, 2.46920, 2.43040, 2.30230, 2.22180/

        DATA(qext(1,7,ib,15),ib=1,17)  & ! Band 7, RH = 94%
        /0.00471, 0.00743, 0.01174, 0.01860, 0.02989, 0.04795  &
        ,0.07872, 0.13752, 0.27893, 0.69585, 1.71850, 2.94500  &
        ,3.18610, 2.46820, 2.42810, 2.29850, 2.21840/

        DATA(qext(1,7,ib,16),ib=1,17)  & ! Band 7, RH = 95%
        /0.00469, 0.00739, 0.01166, 0.01844, 0.02963, 0.04751  &
        ,0.07800, 0.13645, 0.27737, 0.68746, 1.67330, 2.88290  &
        ,3.18940, 2.46930, 2.42370, 2.29340, 2.21410/

        DATA(qext(1,7,ib,17),ib=1,17)  & ! Band 7, RH = 96%
        /0.00467, 0.00734, 0.01155, 0.01823, 0.02926, 0.04689  &
        ,0.07701, 0.13501, 0.27544, 0.67690, 1.61940, 2.80360  &
        ,3.18460, 2.47420, 2.41530, 2.28590, 2.20840/

        DATA(qext(1,7,ib,18),ib=1,17)  & ! Band 7, RH = 97%
        /0.00464, 0.00727, 0.01140, 0.01794, 0.02874, 0.04602  &
        ,0.07562, 0.13310, 0.27332, 0.66393, 1.55590, 2.70240  &
        ,3.16540, 2.48320, 2.40070, 2.27470, 2.20060/

        DATA(qext(1,7,ib,19),ib=1,17)  & ! Band 7, RH = 98%
        /0.00460, 0.00717, 0.01118, 0.01753, 0.02798, 0.04476  &
        ,0.07372, 0.13090, 0.27257, 0.65053, 1.48570, 2.58440  &
        ,3.13330, 2.48800, 2.37430, 2.25700, 2.18850/

        DATA(qext(1,7,ib,20),ib=1,17)  & ! Band 7, RH = 99%
        /0.00455, 0.00702, 0.01086, 0.01696, 0.02694, 0.04324  &
        ,0.07208, 0.13183, 0.28505, 0.66024, 1.44830, 2.52130  &
        ,3.07560, 2.44530, 2.33400, 2.22740, 2.16620/

        DATA(qext(1,8,ib,1),ib=1,17)  & ! Band 8, RH = 80%
        /0.00163, 0.00259, 0.00411, 0.00654, 0.01041, 0.01685  &
        ,0.02841, 0.05378, 0.12576, 0.32198, 0.82147, 1.75220  &
        ,2.80160, 2.74470, 2.40400, 2.22710, 2.17840/

        DATA(qext(1,8,ib,2),ib=1,17)  & ! Band 8, RH = 81%
        /0.00163, 0.00259, 0.00411, 0.00654, 0.01042, 0.01686  &
        ,0.02845, 0.05398, 0.12665, 0.32480, 0.82901, 1.76620  &
        ,2.81420, 2.74120, 2.40490, 2.22820, 2.17680/

        DATA(qext(1,8,ib,3),ib=1,17)  & ! Band 8, RH = 82%
        /0.00163, 0.00259, 0.00411, 0.00655, 0.01042, 0.01687  &
        ,0.02850, 0.05420, 0.12765, 0.32794, 0.83738, 1.78160  &
        ,2.82800, 2.73720, 2.40590, 2.22950, 2.17510/

        DATA(qext(1,8,ib,4),ib=1,17)  & ! Band 8, RH = 83%
        /0.00163, 0.00259, 0.00412, 0.00655, 0.01043, 0.01689  &
        ,0.02855, 0.05446, 0.12876, 0.33145, 0.84670, 1.79860  &
        ,2.84300, 2.73280, 2.40700, 2.23120, 2.17340/

        DATA(qext(1,8,ib,5),ib=1,17)  & ! Band 8, RH = 84%
        /0.00163, 0.00259, 0.00412, 0.00655, 0.01043, 0.01690  &
        ,0.02861, 0.05474, 0.13002, 0.33540, 0.85717, 1.81760  &
        ,2.85950, 2.72770, 2.40810, 2.23330, 2.17160/

        DATA(qext(1,8,ib,6),ib=1,17)  & ! Band 8, RH = 85%
        /0.00163, 0.00259, 0.00412, 0.00656, 0.01044, 0.01692  &
        ,0.02868, 0.05507, 0.13144, 0.33989, 0.86900, 1.83900  &
        ,2.87760, 2.72180, 2.40900, 2.23590, 2.16970/

        DATA(qext(1,8,ib,7),ib=1,17)  & ! Band 8, RH = 86%
        /0.00163, 0.00259, 0.00412, 0.00657, 0.01045, 0.01694  &
        ,0.02876, 0.05544, 0.13308, 0.34503, 0.88249, 1.86310  &
        ,2.89760, 2.71470, 2.40920, 2.23890, 2.16790/

        DATA(qext(1,8,ib,8),ib=1,17)  & ! Band 8, RH = 87%
        /0.00163, 0.00260, 0.00413, 0.00657, 0.01046, 0.01697  &
        ,0.02885, 0.05588, 0.13498, 0.35098, 0.89800, 1.89050  &
        ,2.91970, 2.70590, 2.40830, 2.24230, 2.16640/

        DATA(qext(1,8,ib,9),ib=1,17)  & ! Band 8, RH = 88%
        /0.00163, 0.00260, 0.00413, 0.00658, 0.01047, 0.01700  &
        ,0.02897, 0.05640, 0.13720, 0.35793, 0.91602, 1.92210  &
        ,2.94410, 2.69470, 2.40600, 2.24610, 2.16530/

        DATA(qext(1,8,ib,10),ib=1,17)  & ! Band 8, RH = 89%
        /0.00164, 0.00260, 0.00414, 0.00659, 0.01048, 0.01704  &
        ,0.02911, 0.05702, 0.13983, 0.36618, 0.93722, 1.95870  &
        ,2.97090, 2.68010, 2.40290, 2.25010, 2.16490/

        DATA(qext(1,8,ib,11),ib=1,17)  & ! Band 8, RH = 90%
        /0.00164, 0.00260, 0.00414, 0.00660, 0.01050, 0.01709  &
        ,0.02927, 0.05778, 0.14300, 0.37610, 0.96251, 2.00160  &
        ,3.00020, 2.66140, 2.39900, 2.25370, 2.16550/

        DATA(qext(1,8,ib,12),ib=1,17)  & ! Band 8, RH = 91%
        /0.00164, 0.00261, 0.00415, 0.00662, 0.01053, 0.01715  &
        ,0.02948, 0.05871, 0.14690, 0.38828, 0.99316, 2.05240  &
        ,3.03170, 2.63890, 2.39260, 2.25630, 2.16720/

        DATA(qext(1,8,ib,13),ib=1,17)  & ! Band 8, RH = 92%
        /0.00164, 0.00261, 0.00416, 0.00664, 0.01056, 0.01723  &
        ,0.02976, 0.05991, 0.15178, 0.40358, 1.03100, 2.11360  &
        ,3.06510, 2.61260, 2.38150, 2.25620, 2.16960/

        DATA(qext(1,8,ib,14),ib=1,17)  & ! Band 8, RH = 93%
        /0.00164, 0.00262, 0.00418, 0.00667, 0.01061, 0.01733  &
        ,0.03012, 0.06147, 0.15808, 0.42334, 1.07880, 2.18790  &
        ,3.10060, 2.57890, 2.36740, 2.25160, 2.17130/

        DATA(qext(1,8,ib,15),ib=1,17)  & ! Band 8, RH = 94%
        /0.00165, 0.00263, 0.00419, 0.00670, 0.01067, 0.01748  &
        ,0.03061, 0.06360, 0.16650, 0.44979, 1.14070, 2.27980  &
        ,3.13870, 2.53230, 2.34560, 2.23960, 2.16910/

        DATA(qext(1,8,ib,16),ib=1,17)  & ! Band 8, RH = 95%
        /0.00165, 0.00264, 0.00422, 0.00676, 0.01076, 0.01770  &
        ,0.03133, 0.06666, 0.17829, 0.48689, 1.22340, 2.39660  &
        ,3.17410, 2.47980, 2.31710, 2.21940, 2.15990/

        DATA(qext(1,8,ib,17),ib=1,17)  & ! Band 8, RH = 96%
        /0.00166, 0.00266, 0.00426, 0.00683, 0.01090, 0.01804  &
        ,0.03244, 0.07138, 0.19586, 0.54214, 1.33890, 2.55320  &
        ,3.18710, 2.41090, 2.28240, 2.19910, 2.14960/

        DATA(qext(1,8,ib,18),ib=1,17)  & ! Band 8, RH = 97%
        /0.00167, 0.00268, 0.00432, 0.00696, 0.01114, 0.01863  &
        ,0.03436, 0.07949, 0.22470, 0.63126, 1.51360, 2.77130  &
        ,3.15860, 2.34700, 2.25630, 2.20040, 2.15260/

        DATA(qext(1,8,ib,19),ib=1,17)  & ! Band 8, RH = 98%
        /0.00168, 0.00272, 0.00441, 0.00719, 0.01160, 0.01979  &
        ,0.03834, 0.09623, 0.28030, 0.79283, 1.82160, 3.05890  &
        ,2.99010, 2.32150, 2.27670, 2.20150, 2.13750/

        DATA(qext(1,8,ib,20),ib=1,17)  & ! Band 8, RH = 99%
        /0.00170, 0.00279, 0.00462, 0.00771, 0.01277, 0.02302  &
        ,0.05034, 0.14638, 0.43704, 1.19370, 2.46260, 3.34130  &
        ,2.49980, 2.33230, 2.21850, 2.16760, 2.12290/

! Qext for Sea Salt:
        DATA(qext(2,1,ib,1),ib=1,17)  & ! Band 1, RH = 80%
        /0.00349, 0.00585, 0.01063, 0.02441, 0.07376, 0.26299  &
        ,0.84952, 2.05200, 3.35470, 2.85810, 2.47760, 2.33190  &
        ,2.20670, 2.15340, 2.10590, 2.08120, 2.05990/

        DATA(qext(2,1,ib,2),ib=1,17)  & ! Band 1, RH = 81%
        /0.00357, 0.00599, 0.01089, 0.02518, 0.07646, 0.27096  &
        ,0.86962, 2.08480, 3.36910, 2.82990, 2.48150, 2.33470  &
        ,2.21820, 2.14500, 2.11100, 2.07930, 2.06000/

        DATA(qext(2,1,ib,3),ib=1,17)  & ! Band 1, RH = 82%
        /0.00365, 0.00611, 0.01117, 0.02600, 0.07936, 0.27951  &
        ,0.89098, 2.11900, 3.38140, 2.80200, 2.47430, 2.32240  &
        ,2.21670, 2.14250, 2.10880, 2.08250, 2.05990/

        DATA(qext(2,1,ib,4),ib=1,17)  & ! Band 1, RH = 83%
        /0.00372, 0.00625, 0.01146, 0.02689, 0.08251, 0.28878  &
        ,0.91393, 2.15500, 3.39170, 2.76950, 2.48360, 2.31850  &
        ,2.21360, 2.13720, 2.11040, 2.08130, 2.05790/

        DATA(qext(2,1,ib,5),ib=1,17)  & ! Band 1, RH = 84%
        /0.00380, 0.00639, 0.01177, 0.02786, 0.08596, 0.29893  &
        ,0.93888, 2.19320, 3.40050, 2.73050, 2.48410, 2.31380  &
        ,2.23340, 2.13990, 2.11120, 2.07980, 2.05830/

        DATA(qext(2,1,ib,6),ib=1,17)  & ! Band 1, RH = 85%
        /0.00389, 0.00655, 0.01211, 0.02894, 0.08979, 0.31019  &
        ,0.96636, 2.23450, 3.40910, 2.68780, 2.47120, 2.30470  &
        ,2.22760, 2.14050, 2.10520, 2.07710, 2.05760/

        DATA(qext(2,1,ib,7),ib=1,17)  & ! Band 1, RH = 86%
        /0.00397, 0.00671, 0.01248, 0.03014, 0.09408, 0.32282  &
        ,0.99705, 2.27990, 3.41940, 2.65090, 2.46890, 2.30380  &
        ,2.22180, 2.14070, 2.10180, 2.07730, 2.05520/

        DATA(qext(2,1,ib,8),ib=1,17)  & ! Band 1, RH = 87%
        /0.00407, 0.00689, 0.01290, 0.03150, 0.09895, 0.33721  &
        ,1.03180, 2.33070, 3.43100, 2.61480, 2.46280, 2.28200  &
        ,2.21590, 2.14180, 2.10560, 2.07630, 2.05630/

        DATA(qext(2,1,ib,9),ib=1,17)  & ! Band 1, RH = 88%
        /0.00418, 0.00709, 0.01336, 0.03307, 0.10453, 0.35383  &
        ,1.07160, 2.38810, 3.43980, 2.56950, 2.44850, 2.27090  &
        ,2.20540, 2.13770, 2.10780, 2.07520, 2.05480/

        DATA(qext(2,1,ib,10),ib=1,17)  & ! Band 1, RH = 89%
        /0.00429, 0.00731, 0.01389, 0.03490, 0.11102, 0.37334  &
        ,1.11760, 2.45290, 3.44060, 2.52130, 2.42610, 2.26750  &
        ,2.19970, 2.13440, 2.10170, 2.07460, 2.05400/

        DATA(qext(2,1,ib,11),ib=1,17)  & ! Band 1, RH = 90%
        /0.00442, 0.00757, 0.01451, 0.03707, 0.11868, 0.39663  &
        ,1.17140, 2.52490, 3.43240, 2.48520, 2.39120, 2.23590  &
        ,2.16620, 2.12900, 2.09940, 2.07080, 2.05470/

        DATA(qext(2,1,ib,12),ib=1,17)  & ! Band 1, RH = 91%
        /0.00457, 0.00786, 0.01524, 0.03968, 0.12784, 0.42494  &
        ,1.23440, 2.60350, 3.42150, 2.43800, 2.35710, 2.22370  &
        ,2.17250, 2.13740, 2.09600, 2.07170, 2.05300/

        DATA(qext(2,1,ib,13),ib=1,17)  & ! Band 1, RH = 92%
        /0.00473, 0.00819, 0.01611, 0.04291, 0.13902, 0.45998  &
        ,1.30850, 2.68960, 3.40300, 2.40630, 2.31340, 2.21730  &
        ,2.16670, 2.13510, 2.09370, 2.06860, 2.05170/

        DATA(qext(2,1,ib,14),ib=1,17)  & ! Band 1, RH = 93%
        /0.00493, 0.00859, 0.01718, 0.04701, 0.15294, 0.50417  &
        ,1.39610, 2.78890, 3.35590, 2.37480, 2.28100, 2.21960  &
        ,2.17270, 2.12740, 2.09420, 2.06590, 2.05000/

        DATA(qext(2,1,ib,15),ib=1,17)  & ! Band 1, RH = 94%
        /0.00516, 0.00908, 0.01855, 0.05238, 0.17077, 0.56086  &
        ,1.50230, 2.90860, 3.29510, 2.36550, 2.23540, 2.23320  &
        ,2.16420, 2.11760, 2.08870, 2.06420, 2.04830/

        DATA(qext(2,1,ib,16),ib=1,17)  & ! Band 1, RH = 95%
        /0.00544, 0.00970, 0.02037, 0.05975, 0.19450, 0.63487  &
        ,1.63810, 3.03920, 3.18750, 2.38200, 2.23630, 2.23270  &
        ,2.15860, 2.11720, 2.08660, 2.06270, 2.04640/

        DATA(qext(2,1,ib,17),ib=1,17)  & ! Band 1, RH = 96%
        /0.00579, 0.01051, 0.02294, 0.07055, 0.22793, 0.73404  &
        ,1.82270, 3.18160, 3.02990, 2.41900, 2.28590, 2.20170  &
        ,2.15960, 2.11370, 2.08370, 2.05950, 2.04270/

        DATA(qext(2,1,ib,18),ib=1,17)  & ! Band 1, RH = 97%
        /0.00626, 0.01167, 0.02696, 0.08791, 0.27982, 0.87708  &
        ,2.07510, 3.32470, 2.79050, 2.45920, 2.32460, 2.21040  &
        ,2.14170, 2.10430, 2.07820, 2.05680, 2.04480/

        DATA(qext(2,1,ib,19),ib=1,17)  & ! Band 1, RH = 98%
        /0.00695, 0.01356, 0.03432, 0.12030, 0.37724, 1.12820  &
        ,2.46260, 3.40740, 2.49130, 2.38910, 2.24180, 2.17780  &
        ,2.12330, 2.09340, 2.07070, 2.05270, 2.04050/

        DATA(qext(2,1,ib,20),ib=1,17)  & ! Band 1, RH = 99%
        /0.00814, 0.01752, 0.05311, 0.20048, 0.63541, 1.67050  &
        ,3.07500, 3.11800, 2.39820, 2.25790, 2.21080, 2.15350  &
        ,2.11390, 2.08310, 2.06090, 2.04480, 2.03340/

        DATA(qext(2,2,ib,1),ib=1,17)  & ! Band 2, RH = 80%
        /0.00073, 0.00478, 0.02966, 0.15050, 0.53722, 1.54670  &
        ,3.05970, 3.42670, 2.44540, 2.26250, 2.22890, 2.15640  &
        ,2.12460, 2.08730, 2.06720, 2.04990, 2.03650/

        DATA(qext(2,2,ib,2),ib=1,17)  & ! Band 2, RH = 81%
        /0.00077, 0.00503, 0.03101, 0.15593, 0.55377, 1.57760  &
        ,3.09020, 3.40500, 2.43890, 2.27260, 2.22590, 2.15460  &
        ,2.12160, 2.09360, 2.06470, 2.04870, 2.03610/

        DATA(qext(2,2,ib,3),ib=1,17)  & ! Band 2, RH = 82%
        /0.00081, 0.00527, 0.03248, 0.16174, 0.57138, 1.61060  &
        ,3.12080, 3.37910, 2.43450, 2.27830, 2.21480, 2.16280  &
        ,2.12660, 2.09620, 2.06490, 2.04830, 2.03510/

        DATA(qext(2,2,ib,4),ib=1,17)  & ! Band 2, RH = 83%
        /0.00085, 0.00555, 0.03409, 0.16803, 0.59022, 1.64620  &
        ,3.15150, 3.34960, 2.43090, 2.29100, 2.19980, 2.15830  &
        ,2.11930, 2.08830, 2.06620, 2.04760, 2.03490/

        DATA(qext(2,2,ib,5),ib=1,17)  & ! Band 2, RH = 84%
        /0.00090, 0.00585, 0.03587, 0.17489, 0.61054, 1.68490  &
        ,3.18270, 3.31850, 2.42300, 2.30130, 2.19480, 2.16140  &
        ,2.12130, 2.08860, 2.06390, 2.04680, 2.03520/

        DATA(qext(2,2,ib,6),ib=1,17)  & ! Band 2, RH = 85%
        /0.00095, 0.00618, 0.03785, 0.18244, 0.63264, 1.72720  &
        ,3.21510, 3.28760, 2.42060, 2.30780, 2.19060, 2.16150  &
        ,2.11790, 2.08930, 2.06310, 2.04700, 2.03460/

        DATA(qext(2,2,ib,7),ib=1,17)  & ! Band 2, RH = 86%
        /0.00100, 0.00656, 0.04009, 0.19087, 0.65691, 1.77390  &
        ,3.25010, 3.25270, 2.41610, 2.31980, 2.18760, 2.16470  &
        ,2.11480, 2.08460, 2.06290, 2.04650, 2.03430/

        DATA(qext(2,2,ib,8),ib=1,17)  & ! Band 2, RH = 87%
        /0.00107, 0.00700, 0.04265, 0.20038, 0.68382, 1.82560  &
        ,3.28870, 3.20910, 2.40750, 2.33030, 2.18670, 2.16050  &
        ,2.11590, 2.08410, 2.06090, 2.04530, 2.03320/

        DATA(qext(2,2,ib,9),ib=1,17)  & ! Band 2, RH = 88%
        /0.00114, 0.00751, 0.04562, 0.21127, 0.71403, 1.88290  &
        ,3.33030, 3.15620, 2.40130, 2.33820, 2.18590, 2.15390  &
        ,2.11220, 2.08320, 2.06180, 2.04550, 2.03290/

        DATA(qext(2,2,ib,10),ib=1,17)  & ! Band 2, RH = 89%
        /0.00123, 0.00810, 0.04911, 0.22392, 0.74842, 1.94670  &
        ,3.37220, 3.10120, 2.39480, 2.35070, 2.18800, 2.15340  &
        ,2.11070, 2.08280, 2.06080, 2.04420, 2.03280/

        DATA(qext(2,2,ib,11),ib=1,17)  & ! Band 2, RH = 90%
        /0.00133, 0.00882, 0.05327, 0.23887, 0.78821, 2.01810  &
        ,3.41150, 3.03950, 2.38500, 2.34780, 2.19880, 2.15390  &
        ,2.11260, 2.08450, 2.06070, 2.04320, 2.03230/

        DATA(qext(2,2,ib,12),ib=1,17)  & ! Band 2, RH = 91%
        /0.00145, 0.00971, 0.05833, 0.25688, 0.83522, 2.09940  &
        ,3.44850, 2.96190, 2.37200, 2.35100, 2.20660, 2.15880  &
        ,2.11170, 2.08000, 2.05730, 2.04270, 2.03170/

        DATA(qext(2,2,ib,13),ib=1,17)  & ! Band 2, RH = 92%
        /0.00161, 0.01081, 0.06459, 0.27911, 0.89214, 2.19430  &
        ,3.48840, 2.88390, 2.35320, 2.33630, 2.21150, 2.15930  &
        ,2.10720, 2.07590, 2.05600, 2.04210, 2.03100/

        DATA(qext(2,2,ib,14),ib=1,17)  & ! Band 2, RH = 93%
        /0.00180, 0.01223, 0.07253, 0.30737, 0.96298, 2.30810  &
        ,3.52580, 2.78820, 2.32890, 2.31120, 2.20410, 2.14200  &
        ,2.10330, 2.07500, 2.05530, 2.04220, 2.03110/

        DATA(qext(2,2,ib,15),ib=1,17)  & ! Band 2, RH = 94%
        /0.00206, 0.01413, 0.08288, 0.34459, 1.05340, 2.44450  &
        ,3.54420, 2.68990, 2.29220, 2.28850, 2.19630, 2.13620  &
        ,2.10310, 2.07400, 2.05410, 2.03940, 2.02820/

        DATA(qext(2,2,ib,16),ib=1,17)  & ! Band 2, RH = 95%
        /0.00242, 0.01676, 0.09691, 0.39573, 1.17090, 2.60460  &
        ,3.54950, 2.58520, 2.25350, 2.27360, 2.19760, 2.13130  &
        ,2.09650, 2.06900, 2.05110, 2.03880, 2.02800/

        DATA(qext(2,2,ib,17),ib=1,17)  & ! Band 2, RH = 96%
        /0.00293, 0.02066, 0.11694, 0.46932, 1.32710, 2.80500  &
        ,3.50020, 2.49110, 2.22740, 2.26050, 2.17740, 2.12660  &
        ,2.09100, 2.06630, 2.04970, 2.03680, 2.02690/

        DATA(qext(2,2,ib,18),ib=1,17)  & ! Band 2, RH = 97%
        /0.00374, 0.02698, 0.14786, 0.58054, 1.55110, 3.05690  &
        ,3.36840, 2.42170, 2.26450, 2.20950, 2.14900, 2.11740  &
        ,2.08720, 2.06380, 2.04620, 2.03400, 2.02590/

        DATA(qext(2,2,ib,19),ib=1,17)  & ! Band 2, RH = 98%
        /0.00520, 0.03881, 0.20254, 0.76450, 1.91460, 3.35310  &
        ,3.06610, 2.38510, 2.34530, 2.19040, 2.15060, 2.10930  &
        ,2.08170, 2.05850, 2.04290, 2.03050, 2.02300/

        DATA(qext(2,2,ib,20),ib=1,17)  & ! Band 2, RH = 99%
        /0.00854, 0.06789, 0.33643, 1.18010, 2.58590, 3.53320  &
        ,2.54050, 2.23250, 2.26170, 2.19240, 2.12710, 2.09350  &
        ,2.06850, 2.05000, 2.03640, 2.02810, 2.02020/

        DATA(qext(2,3,ib,1),ib=1,17)  & ! Band 3, RH = 80%
        /0.01390, 0.07864, 0.33227, 1.04460, 2.28450, 3.49860  &
        ,2.98590, 2.19520, 2.34920, 2.22500, 2.14170, 2.11100  &
        ,2.08200, 2.06070, 2.04250, 2.03250, 2.02220/

        DATA(qext(2,3,ib,2),ib=1,17)  & ! Band 3, RH = 81%
        /0.01464, 0.08201, 0.34275, 1.06840, 2.31620, 3.50960  &
        ,2.95390, 2.20050, 2.34220, 2.22620, 2.13860, 2.10850  &
        ,2.08410, 2.05850, 2.04410, 2.03120, 2.02240/

        DATA(qext(2,3,ib,3),ib=1,17)  & ! Band 3, RH = 82%
        /0.01537, 0.08541, 0.35400, 1.09360, 2.34930, 3.52010  &
        ,2.92280, 2.20660, 2.32180, 2.21840, 2.14490, 2.10860  &
        ,2.07680, 2.05670, 2.04130, 2.03020, 2.02270/

        DATA(qext(2,3,ib,4),ib=1,17)  & ! Band 3, RH = 83%
        /0.01612, 0.08908, 0.36614, 1.12070, 2.38480, 3.53090  &
        ,2.89170, 2.21720, 2.31400, 2.20760, 2.14590, 2.10120  &
        ,2.07690, 2.05470, 2.04240, 2.03040, 2.02310/

        DATA(qext(2,3,ib,5),ib=1,17)  & ! Band 3, RH = 84%
        /0.01694, 0.09310, 0.37934, 1.14990, 2.42310, 3.53880  &
        ,2.85880, 2.22360, 2.29950, 2.20100, 2.14770, 2.09830  &
        ,2.07260, 2.05660, 2.04210, 2.02990, 2.02290/

        DATA(qext(2,3,ib,6),ib=1,17)  & ! Band 3, RH = 85%
        /0.01784, 0.09754, 0.39382, 1.18200, 2.46470, 3.54480  &
        ,2.82240, 2.23140, 2.28320, 2.19710, 2.14640, 2.10400  &
        ,2.07280, 2.05520, 2.04110, 2.02970, 2.02150/

        DATA(qext(2,3,ib,7),ib=1,17)  & ! Band 3, RH = 86%
        /0.01885, 0.10249, 0.40984, 1.21740, 2.50980, 3.54840  &
        ,2.78170, 2.24250, 2.26770, 2.18710, 2.15200, 2.09650  &
        ,2.07650, 2.05450, 2.03980, 2.02960, 2.02250/

        DATA(qext(2,3,ib,8),ib=1,17)  & ! Band 3, RH = 87%
        /0.01999, 0.10807, 0.42771, 1.25700, 2.55930, 3.54920  &
        ,2.73530, 2.25490, 2.25220, 2.18260, 2.14880, 2.10060  &
        ,2.07430, 2.05380, 2.03970, 2.03000, 2.02120/

        DATA(qext(2,3,ib,9),ib=1,17)  & ! Band 3, RH = 88%
        /0.02130, 0.11445, 0.44784, 1.30140, 2.61360, 3.55050  &
        ,2.68520, 2.27430, 2.23320, 2.16970, 2.14060, 2.10390  &
        ,2.07340, 2.05340, 2.03930, 2.02810, 2.02020/

        DATA(qext(2,3,ib,10),ib=1,17)  & ! Band 3, RH = 89%
        /0.02282, 0.12184, 0.47072, 1.35140, 2.67260, 3.54930  &
        ,2.63090, 2.28640, 2.21800, 2.16130, 2.12830, 2.09900  &
        ,2.07080, 2.05160, 2.03850, 2.02900, 2.02030/

        DATA(qext(2,3,ib,11),ib=1,17)  & ! Band 3, RH = 90%
        /0.02461, 0.13049, 0.49706, 1.40850, 2.73990, 3.54170  &
        ,2.57320, 2.30470, 2.20080, 2.15640, 2.13060, 2.09570  &
        ,2.06830, 2.05220, 2.03740, 2.02830, 2.02120/

        DATA(qext(2,3,ib,12),ib=1,17)  & ! Band 3, RH = 91%
        /0.02676, 0.14077, 0.52789, 1.47490, 2.81590, 3.52420  &
        ,2.51200, 2.32930, 2.18840, 2.16050, 2.12120, 2.09040  &
        ,2.06970, 2.05110, 2.03810, 2.02760, 2.01930/

        DATA(qext(2,3,ib,13),ib=1,17)  & ! Band 3, RH = 92%
        /0.02940, 0.15318, 0.56477, 1.55310, 2.90210, 3.49850  &
        ,2.44450, 2.34220, 2.18330, 2.16790, 2.11720, 2.08720  &
        ,2.06750, 2.05040, 2.03660, 2.02660, 2.02010/

        DATA(qext(2,3,ib,14),ib=1,17)  & ! Band 3, RH = 93%
        /0.03270, 0.16840, 0.61014, 1.64650, 2.99490, 3.45990  &
        ,2.36890, 2.36290, 2.18670, 2.17930, 2.12060, 2.09160  &
        ,2.06580, 2.04780, 2.03540, 2.02620, 2.01870/

        DATA(qext(2,3,ib,15),ib=1,17)  & ! Band 3, RH = 94%
        /0.03695, 0.18755, 0.66800, 1.75840, 3.09880, 3.39210  &
        ,2.30260, 2.37550, 2.20480, 2.17800, 2.12780, 2.08800  &
        ,2.06190, 2.04560, 2.03420, 2.02420, 2.01850/

        DATA(qext(2,3,ib,16),ib=1,17)  & ! Band 3, RH = 95%
        /0.04262, 0.21252, 0.74478, 1.89610, 3.21790, 3.30210  &
        ,2.23270, 2.37720, 2.22460, 2.16150, 2.11120, 2.08180  &
        ,2.06340, 2.04430, 2.03370, 2.02530, 2.01850/

        DATA(qext(2,3,ib,17),ib=1,17)  & ! Band 3, RH = 96%
        /0.05057, 0.24710, 0.84932, 2.07510, 3.34410, 3.16020  &
        ,2.18690, 2.36780, 2.23490, 2.14050, 2.10870, 2.07700  &
        ,2.05890, 2.04380, 2.03110, 2.02410, 2.01650/

        DATA(qext(2,3,ib,18),ib=1,17)  & ! Band 3, RH = 97%
        /0.06256, 0.29951, 0.99765, 2.31660, 3.47200, 2.94570  &
        ,2.19220, 2.32250, 2.21080, 2.13720, 2.10380, 2.07590  &
        ,2.05420, 2.04030, 2.03010, 2.02140, 2.01610/

        DATA(qext(2,3,ib,19),ib=1,17)  & ! Band 3, RH = 98%
        /0.08290, 0.38942, 1.23930, 2.67350, 3.53380, 2.61790  &
        ,2.28460, 2.20780, 2.15350, 2.12660, 2.09820, 2.06710  &
        ,2.04930, 2.03730, 2.02680, 2.01970, 2.01480/

        DATA(qext(2,3,ib,20),ib=1,17)  & ! Band 3, RH = 99%
        /0.12530, 0.57315, 1.71420, 3.21310, 3.30290, 2.21790  &
        ,2.38380, 2.23300, 2.15730, 2.10010, 2.08190, 2.05830  &
        ,2.04370, 2.03130, 2.02260, 2.01610, 2.01320/

        DATA(qext(2,4,ib,1),ib=1,17)  & ! Band 4, RH = 80%
        /0.00639, 0.01033, 0.01658, 0.02650, 0.04061, 0.06489  &
        ,0.10424, 0.17010, 0.28871, 0.52953, 1.02380, 1.76560  &
        ,2.43410, 2.64750, 2.49550, 2.36480, 2.27530/

        DATA(qext(2,4,ib,2),ib=1,17)  & ! Band 4, RH = 81%
        /0.00652, 0.01052, 0.01687, 0.02697, 0.04136, 0.06610  &
        ,0.10622, 0.17347, 0.29498, 0.54260, 1.04820, 1.79500  &
        ,2.45320, 2.64710, 2.49060, 2.36200, 2.27320/

        DATA(qext(2,4,ib,3),ib=1,17)  & ! Band 4, RH = 82%
        /0.00664, 0.01071, 0.01718, 0.02747, 0.04214, 0.06736  &
        ,0.10829, 0.17701, 0.30162, 0.55646, 1.07370, 1.82530  &
        ,2.47220, 2.64620, 2.48570, 2.35910, 2.27100/

        DATA(qext(2,4,ib,4),ib=1,17)  & ! Band 4, RH = 83%
        /0.00676, 0.01091, 0.01750, 0.02798, 0.04297, 0.06869  &
        ,0.11047, 0.18076, 0.30869, 0.57126, 1.10070, 1.85660  &
        ,2.49110, 2.64460, 2.48050, 2.35600, 2.26870/

        DATA(qext(2,4,ib,5),ib=1,17)  & ! Band 4, RH = 84%
        /0.00689, 0.01112, 0.01784, 0.02854, 0.04384, 0.07011  &
        ,0.11281, 0.18477, 0.31628, 0.58721, 1.12950, 1.88930  &
        ,2.51000, 2.64230, 2.47520, 2.35290, 2.26630/

        DATA(qext(2,4,ib,6),ib=1,17)  & ! Band 4, RH = 85%
        /0.00702, 0.01134, 0.01821, 0.02913, 0.04478, 0.07163  &
        ,0.11532, 0.18911, 0.32454, 0.60458, 1.16030, 1.92360  &
        ,2.52880, 2.63920, 2.46960, 2.34950, 2.26370/

        DATA(qext(2,4,ib,7),ib=1,17)  & ! Band 4, RH = 86%
        /0.00716, 0.01158, 0.01860, 0.02977, 0.04580, 0.07328  &
        ,0.11804, 0.19384, 0.33359, 0.62366, 1.19380, 1.95990  &
        ,2.54750, 2.63530, 2.46380, 2.34590, 2.26100/

        DATA(qext(2,4,ib,8),ib=1,17)  & ! Band 4, RH = 87%
        /0.00731, 0.01183, 0.01902, 0.03046, 0.04691, 0.07509  &
        ,0.12104, 0.19905, 0.34363, 0.64486, 1.23020, 1.99850  &
        ,2.56620, 2.63030, 2.45760, 2.34200, 2.25800/

        DATA(qext(2,4,ib,9),ib=1,17)  & ! Band 4, RH = 88%
        /0.00748, 0.01212, 0.01949, 0.03123, 0.04814, 0.07709  &
        ,0.12436, 0.20485, 0.35490, 0.66865, 1.27040, 2.03990  &
        ,2.58470, 2.62410, 2.45100, 2.33780, 2.25480/

        DATA(qext(2,4,ib,10),ib=1,17)  & ! Band 4, RH = 89%
        /0.00766, 0.01243, 0.02002, 0.03209, 0.04951, 0.07933  &
        ,0.12808, 0.21139, 0.36770, 0.69565, 1.31500, 2.08430  &
        ,2.60300, 2.61660, 2.44400, 2.33320, 2.25130/

        DATA(qext(2,4,ib,11),ib=1,17)  & ! Band 4, RH = 90%
        /0.00787, 0.01278, 0.02061, 0.03306, 0.05106, 0.08186  &
        ,0.13231, 0.21886, 0.38243, 0.72666, 1.36490, 2.13220  &
        ,2.62080, 2.60740, 2.43640, 2.32810, 2.24740/

        DATA(qext(2,4,ib,12),ib=1,17)  & ! Band 4, RH = 91%
        /0.00810, 0.01318, 0.02128, 0.03417, 0.05283, 0.08476  &
        ,0.13716, 0.22749, 0.39963, 0.76273, 1.42120, 2.18400  &
        ,2.63770, 2.59620, 2.42810, 2.32250, 2.24320/

        DATA(qext(2,4,ib,13),ib=1,17)  & ! Band 4, RH = 92%
        /0.00836, 0.01364, 0.02205, 0.03545, 0.05488, 0.08813  &
        ,0.14283, 0.23766, 0.42008, 0.80532, 1.48540, 2.23990  &
        ,2.65320, 2.58300, 2.41920, 2.31620, 2.23840/

        DATA(qext(2,4,ib,14),ib=1,17)  & ! Band 4, RH = 93%
        /0.00867, 0.01418, 0.02296, 0.03696, 0.05731, 0.09213  &
        ,0.14958, 0.24986, 0.44492, 0.85650, 1.55930, 2.30050  &
        ,2.66600, 2.56720, 2.40920, 2.30910, 2.23290/

        DATA(qext(2,4,ib,15),ib=1,17)  & ! Band 4, RH = 94%
        /0.00903, 0.01482, 0.02406, 0.03878, 0.06024, 0.09697  &
        ,0.15782, 0.26494, 0.47599, 0.91943, 1.64560, 2.36600  &
        ,2.67460, 2.54850, 2.39810, 2.30080, 2.22660/

        DATA(qext(2,4,ib,16),ib=1,17)  & ! Band 4, RH = 95%
        /0.00947, 0.01561, 0.02541, 0.04106, 0.06390, 0.10306  &
        ,0.16825, 0.28427, 0.51639, 0.99916, 1.74810, 2.43730  &
        ,2.67710, 2.52640, 2.38520, 2.29120, 2.21930/

        DATA(qext(2,4,ib,17),ib=1,17)  & ! Band 4, RH = 96%
        /0.01001, 0.01661, 0.02716, 0.04401, 0.06870, 0.11107  &
        ,0.18212, 0.31048, 0.57190, 1.10440, 1.87310, 2.51460  &
        ,2.67010, 2.49990, 2.36980, 2.27940, 2.21030/

        DATA(qext(2,4,ib,18),ib=1,17)  & ! Band 4, RH = 97%
        /0.01074, 0.01797, 0.02957, 0.04812, 0.07542, 0.12242  &
        ,0.20207, 0.34912, 0.65467, 1.25180, 2.03180, 2.59460  &
        ,2.64740, 2.46830, 2.35040, 2.26450, 2.19900/

        DATA(qext(2,4,ib,19),ib=1,17)  & ! Band 4, RH = 98%
        /0.01178, 0.01999, 0.03324, 0.05452, 0.08603, 0.14061  &
        ,0.23485, 0.41494, 0.79538, 1.47760, 2.24250, 2.66530  &
        ,2.59810, 2.42900, 2.32400, 2.24430, 2.18350/

        DATA(qext(2,4,ib,20),ib=1,17)  & ! Band 4, RH = 99%
        /0.01349, 0.02361, 0.04019, 0.06707, 0.10742, 0.17855  &
        ,0.30676, 0.56758, 1.10240, 1.87870, 2.52380, 2.67720  &
        ,2.50420, 2.37300, 2.28190, 2.21220, 2.15910/

        DATA(qext(2,5,ib,1),ib=1,17)  & ! Band 5, RH = 80%
        /0.01209, 0.01960, 0.03151, 0.05042, 0.07698, 0.12332  &
        ,0.19914, 0.32792, 0.55770, 0.95268, 1.49610, 2.03990  &
        ,2.35980, 2.37690, 2.29660, 2.23220, 2.17940/

        DATA(qext(2,5,ib,2),ib=1,17)  & ! Band 5, RH = 81%
        /0.01237, 0.02001, 0.03214, 0.05143, 0.07859, 0.12591  &
        ,0.20340, 0.33514, 0.57015, 0.97128, 1.51680, 2.05430  &
        ,2.36200, 2.37290, 2.29380, 2.23010, 2.17790/

        DATA(qext(2,5,ib,3),ib=1,17)  & ! Band 5, RH = 82%
        /0.01263, 0.02042, 0.03279, 0.05248, 0.08026, 0.12863  &
        ,0.20787, 0.34272, 0.58320, 0.99062, 1.53800, 2.06890  &
        ,2.36410, 2.36880, 2.29090, 2.22800, 2.17620/

        DATA(qext(2,5,ib,4),ib=1,17)  & ! Band 5, RH = 83%
        /0.01288, 0.02084, 0.03348, 0.05359, 0.08203, 0.13149  &
        ,0.21258, 0.35074, 0.59698, 1.01090, 1.56000, 2.08370  &
        ,2.36610, 2.36470, 2.28800, 2.22580, 2.17460/

        DATA(qext(2,5,ib,5),ib=1,17)  & ! Band 5, RH = 84%
        /0.01315, 0.02128, 0.03420, 0.05477, 0.08390, 0.13453  &
        ,0.21760, 0.35929, 0.61167, 1.03220, 1.58290, 2.09880  &
        ,2.36780, 2.36050, 2.28510, 2.22360, 2.17280/

        DATA(qext(2,5,ib,6),ib=1,17)  & ! Band 5, RH = 85%
        /0.01344, 0.02175, 0.03498, 0.05602, 0.08591, 0.13779  &
        ,0.22299, 0.36850, 0.62746, 1.05490, 1.60700, 2.11440  &
        ,2.36940, 2.35620, 2.28200, 2.22120, 2.17100/

        DATA(qext(2,5,ib,7),ib=1,17)  & ! Band 5, RH = 86%
        /0.01374, 0.02226, 0.03581, 0.05738, 0.08807, 0.14132  &
        ,0.22884, 0.37851, 0.64458, 1.07920, 1.63230, 2.13040  &
        ,2.37070, 2.35170, 2.27880, 2.21880, 2.16910/

        DATA(qext(2,5,ib,8),ib=1,17)  & ! Band 5, RH = 87%
        /0.01407, 0.02281, 0.03671, 0.05886, 0.09044, 0.14518  &
        ,0.23525, 0.38950, 0.66334, 1.10560, 1.65930, 2.14710  &
        ,2.37180, 2.34700, 2.27550, 2.21610, 2.16710/

        DATA(qext(2,5,ib,9),ib=1,17)  & ! Band 5, RH = 88%
        /0.01442, 0.02341, 0.03771, 0.06049, 0.09305, 0.14944  &
        ,0.24234, 0.40169, 0.68408, 1.13430, 1.68820, 2.16430  &
        ,2.37260, 2.34210, 2.27190, 2.21330, 2.16490/

        DATA(qext(2,5,ib,10),ib=1,17)  & ! Band 5, RH = 89%
        /0.01481, 0.02407, 0.03881, 0.06230, 0.09596, 0.15420  &
        ,0.25028, 0.41536, 0.70723, 1.16580, 1.71940, 2.18230  &
        ,2.37300, 2.33690, 2.26810, 2.21030, 2.16250/

        DATA(qext(2,5,ib,11),ib=1,17)  & ! Band 5, RH = 90%
        /0.01524, 0.02481, 0.04005, 0.06434, 0.09923, 0.15957  &
        ,0.25926, 0.43087, 0.73334, 1.20080, 1.75320, 2.20110  &
        ,2.37280, 2.33140, 2.26400, 2.20710, 2.15990/

        DATA(qext(2,5,ib,12),ib=1,17)  & ! Band 5, RH = 91%
        /0.01572, 0.02565, 0.04146, 0.06667, 0.10297, 0.16570  &
        ,0.26955, 0.44870, 0.76313, 1.24000, 1.79020, 2.22070  &
        ,2.37200, 2.32560, 2.25960, 2.20350, 2.15710/

        DATA(qext(2,5,ib,13),ib=1,17)  & ! Band 5, RH = 92%
        /0.01628, 0.02661, 0.04309, 0.06935, 0.10729, 0.17282  &
        ,0.28153, 0.46949, 0.79755, 1.28440, 1.83080, 2.24120  &
        ,2.37030, 2.31920, 2.25470, 2.19950, 2.15400/

        DATA(qext(2,5,ib,14),ib=1,17)  & ! Band 5, RH = 93%
        /0.01692, 0.02774, 0.04499, 0.07251, 0.11238, 0.18123  &
        ,0.29574, 0.49422, 0.83796, 1.33520, 1.87590, 2.26250  &
        ,2.36750, 2.31230, 2.24920, 2.19510, 2.15050/

        DATA(qext(2,5,ib,15),ib=1,17)  & ! Band 5, RH = 94%
        /0.01767, 0.02907, 0.04727, 0.07631, 0.11851, 0.19140  &
        ,0.31302, 0.52435, 0.88637, 1.39450, 1.92640, 2.28440  &
        ,2.36330, 2.30460, 2.24300, 2.19000, 2.14640/

        DATA(qext(2,5,ib,16),ib=1,17)  & ! Band 5, RH = 95%
        /0.01858, 0.03071, 0.05008, 0.08102, 0.12614, 0.20413  &
        ,0.33475, 0.56232, 0.94596, 1.46520, 1.98390, 2.30640  &
        ,2.35720, 2.29590, 2.23580, 2.18420, 2.14180/

        DATA(qext(2,5,ib,17),ib=1,17)  & ! Band 5, RH = 96%
        /0.01972, 0.03279, 0.05369, 0.08712, 0.13609, 0.22081  &
        ,0.36348, 0.61250, 1.02220, 1.55190, 2.05010, 2.32730  &
        ,2.34830, 2.28570, 2.22720, 2.17710, 2.13610/

        DATA(qext(2,5,ib,18),ib=1,17)  & ! Band 5, RH = 97%
        /0.02122, 0.03558, 0.05864, 0.09558, 0.14998, 0.24434  &
        ,0.40440, 0.68366, 1.12520, 1.66280, 2.12760, 2.34560  &
        ,2.33560, 2.27300, 2.21630, 2.16820, 2.12900/

        DATA(qext(2,5,ib,19),ib=1,17)  & ! Band 5, RH = 98%
        /0.02334, 0.03972, 0.06616, 0.10869, 0.17179, 0.28185  &
        ,0.47054, 0.79670, 1.27790, 1.81370, 2.21820, 2.35610  &
        ,2.31650, 2.25560, 2.20140, 2.15600, 2.11930/

        DATA(qext(2,5,ib,20),ib=1,17)  & ! Band 5, RH = 99%
        /0.02685, 0.04709, 0.08028, 0.13427, 0.21559, 0.35945  &
        ,0.60968, 1.02060, 1.54900, 2.04310, 2.31850, 2.34380  &
        ,2.28470, 2.22730, 2.17750, 2.13660, 2.10400/

        DATA(qext(2,6,ib,1),ib=1,17)  & ! Band 6, RH = 80%
        /0.00213, 0.00345, 0.00556, 0.00893, 0.01376, 0.02276  &
        ,0.04070, 0.08782, 0.23711, 0.63424, 1.54460, 2.85800  &
        ,3.34830, 2.16410, 2.41000, 2.23520, 2.15250/

        DATA(qext(2,6,ib,2),ib=1,17)  & ! Band 6, RH = 81%
        /0.00218, 0.00353, 0.00567, 0.00911, 0.01405, 0.02327  &
        ,0.04172, 0.09039, 0.24392, 0.65093, 1.57260, 2.88340  &
        ,3.32840, 2.15330, 2.39860, 2.22840, 2.15190/

        DATA(qext(2,6,ib,3),ib=1,17)  & ! Band 6, RH = 82%
        /0.00222, 0.00360, 0.00579, 0.00930, 0.01436, 0.02381  &
        ,0.04279, 0.09314, 0.25116, 0.66887, 1.60200, 2.90880  &
        ,3.30480, 2.14640, 2.38540, 2.22100, 2.15110/

        DATA(qext(2,6,ib,4),ib=1,17)  & ! Band 6, RH = 83%
        /0.00227, 0.00367, 0.00591, 0.00950, 0.01468, 0.02438  &
        ,0.04394, 0.09611, 0.25891, 0.68830, 1.63310, 2.93380  &
        ,3.27840, 2.13970, 2.37090, 2.21320, 2.14980/

        DATA(qext(2,6,ib,5),ib=1,17)  & ! Band 6, RH = 84%
        /0.00232, 0.00375, 0.00604, 0.00971, 0.01503, 0.02498  &
        ,0.04518, 0.09934, 0.26730, 0.70949, 1.66610, 2.95830  &
        ,3.25100, 2.13040, 2.35460, 2.20510, 2.14820/

        DATA(qext(2,6,ib,6),ib=1,17)  & ! Band 6, RH = 85%
        /0.00237, 0.00384, 0.00618, 0.00993, 0.01540, 0.02563  &
        ,0.04652, 0.10289, 0.27645, 0.73278, 1.70100, 2.98230  &
        ,3.22320, 2.12250, 2.33660, 2.19700, 2.14620/

        DATA(qext(2,6,ib,7),ib=1,17)  & ! Band 6, RH = 86%
        /0.00242, 0.00393, 0.00633, 0.01018, 0.01580, 0.02634  &
        ,0.04800, 0.10684, 0.28653, 0.75853, 1.73830, 3.00700  &
        ,3.19320, 2.12200, 2.31720, 2.18910, 2.14370/

        DATA(qext(2,6,ib,8),ib=1,17)  & ! Band 6, RH = 87%
        /0.00248, 0.00402, 0.00649, 0.01044, 0.01624, 0.02712  &
        ,0.04965, 0.11129, 0.29776, 0.78720, 1.77840, 3.03400  &
        ,3.15690, 2.12510, 2.29540, 2.18190, 2.14100/

        DATA(qext(2,6,ib,9),ib=1,17)  & ! Band 6, RH = 88%
        /0.00254, 0.00413, 0.00666, 0.01074, 0.01672, 0.02799  &
        ,0.05151, 0.11635, 0.31041, 0.81926, 1.82210, 3.06570  &
        ,3.11110, 2.12550, 2.27320, 2.17590, 2.13810/

        DATA(qext(2,6,ib,10),ib=1,17)  & ! Band 6, RH = 89%
        /0.00261, 0.00425, 0.00686, 0.01106, 0.01726, 0.02897  &
        ,0.05363, 0.12220, 0.32482, 0.85527, 1.87090, 3.10360  &
        ,3.05750, 2.13550, 2.24940, 2.17170, 2.13520/

        DATA(qext(2,6,ib,11),ib=1,17)  & ! Band 6, RH = 90%
        /0.00269, 0.00438, 0.00708, 0.01143, 0.01787, 0.03008  &
        ,0.05609, 0.12905, 0.34147, 0.89586, 1.92680, 3.14610  &
        ,2.99950, 2.15290, 2.22670, 2.16980, 2.13260/

        DATA(qext(2,6,ib,12),ib=1,17)  & ! Band 6, RH = 91%
        /0.00277, 0.00453, 0.00733, 0.01185, 0.01857, 0.03137  &
        ,0.05899, 0.13720, 0.36101, 0.94179, 1.99290, 3.18670  &
        ,2.93040, 2.17320, 2.20670, 2.17040, 2.13020/

        DATA(qext(2,6,ib,13),ib=1,17)  & ! Band 6, RH = 92%
        /0.00287, 0.00470, 0.00762, 0.01234, 0.01939, 0.03289  &
        ,0.06247, 0.14709, 0.38439, 0.99421, 2.07290, 3.21980  &
        ,2.84340, 2.20760, 2.19170, 2.17280, 2.12780/

        DATA(qext(2,6,ib,14),ib=1,17)  & ! Band 6, RH = 93%
        /0.00299, 0.00490, 0.00796, 0.01291, 0.02035, 0.03472  &
        ,0.06675, 0.15938, 0.41314, 1.05520, 2.16900, 3.25110  &
        ,2.74560, 2.24990, 2.18570, 2.17450, 2.12470/

        DATA(qext(2,6,ib,15),ib=1,17)  & ! Band 6, RH = 94%
        /0.00312, 0.00514, 0.00837, 0.01360, 0.02152, 0.03697  &
        ,0.07220, 0.17510, 0.44983, 1.12920, 2.28080, 3.29000  &
        ,2.62430, 2.29930, 2.19220, 2.17170, 2.12060/

        DATA(qext(2,6,ib,16),ib=1,17)  & ! Band 6, RH = 95%
        /0.00328, 0.00543, 0.00887, 0.01446, 0.02299, 0.03987  &
        ,0.07946, 0.19601, 0.49931, 1.22580, 2.40990, 3.30810  &
        ,2.48670, 2.35200, 2.21070, 2.16180, 2.11610/

        DATA(qext(2,6,ib,17),ib=1,17)  & ! Band 6, RH = 96%
        /0.00348, 0.00580, 0.00952, 0.01558, 0.02494, 0.04383  &
        ,0.08976, 0.22538, 0.57139, 1.36440, 2.57810, 3.29950  &
        ,2.33350, 2.38660, 2.22800, 2.14920, 2.11090/

        DATA(qext(2,6,ib,18),ib=1,17)  & ! Band 6, RH = 97%
        /0.00375, 0.00629, 0.01041, 0.01714, 0.02771, 0.04973  &
        ,0.10589, 0.27002, 0.68707, 1.57060, 2.79420, 3.22210  &
        ,2.18910, 2.36130, 2.21230, 2.14280, 2.10420/

        DATA(qext(2,6,ib,19),ib=1,17)  & ! Band 6, RH = 98%
        /0.00413, 0.00703, 0.01176, 0.01959, 0.03222, 0.06002  &
        ,0.13567, 0.34740, 0.88707, 1.87840, 3.06290, 2.99080  &
        ,2.15230, 2.23840, 2.16860, 2.12980, 2.09540/

        DATA(qext(2,6,ib,20),ib=1,17)  & ! Band 6, RH = 99%
        /0.00475, 0.00834, 0.01432, 0.02449, 0.04203, 0.08539  &
        ,0.21228, 0.53671, 1.29000, 2.47370, 3.27290, 2.40570  &
        ,2.36280, 2.21670, 2.15200, 2.11140, 2.08170/

        DATA(qext(2,7,ib,1),ib=1,17)  & ! Band 7, RH = 80%
        /0.00239, 0.00388, 0.00624, 0.01001, 0.01536, 0.02502  &
        ,0.04264, 0.08196, 0.19158, 0.47030, 1.13170, 2.22660  &
        ,3.17280, 2.62430, 2.28710, 2.21190, 2.16120/

        DATA(qext(2,7,ib,2),ib=1,17)  & ! Band 7, RH = 81%
        /0.00245, 0.00396, 0.00637, 0.01021, 0.01568, 0.02557  &
        ,0.04364, 0.08411, 0.19678, 0.48106, 1.15000, 2.24760  &
        ,3.17130, 2.59950, 2.29250, 2.21170, 2.15940/

        DATA(qext(2,7,ib,3),ib=1,17)  & ! Band 7, RH = 82%
        /0.00250, 0.00404, 0.00650, 0.01043, 0.01602, 0.02614  &
        ,0.04468, 0.08640, 0.20231, 0.49261, 1.16970, 2.26990  &
        ,3.17000, 2.57430, 2.29710, 2.21150, 2.15750/

        DATA(qext(2,7,ib,4),ib=1,17)  & ! Band 7, RH = 83%
        /0.00255, 0.00413, 0.00664, 0.01065, 0.01638, 0.02674  &
        ,0.04580, 0.08885, 0.20822, 0.50511, 1.19100, 2.29390  &
        ,3.16900, 2.54940, 2.30220, 2.21120, 2.15560/

        DATA(qext(2,7,ib,5),ib=1,17)  & ! Band 7, RH = 84%
        /0.00261, 0.00422, 0.00678, 0.01089, 0.01676, 0.02739  &
        ,0.04699, 0.09149, 0.21460, 0.51877, 1.21430, 2.32010  &
        ,3.16780, 2.52400, 2.30690, 2.21080, 2.15380/

        DATA(qext(2,7,ib,6),ib=1,17)  & ! Band 7, RH = 85%
        /0.00266, 0.00431, 0.00694, 0.01114, 0.01717, 0.02808  &
        ,0.04828, 0.09438, 0.22156, 0.53386, 1.24010, 2.34870  &
        ,3.16550, 2.49650, 2.31050, 2.21010, 2.15190/

        DATA(qext(2,7,ib,7),ib=1,17)  & ! Band 7, RH = 86%
        /0.00272, 0.00441, 0.00711, 0.01141, 0.01762, 0.02883  &
        ,0.04969, 0.09757, 0.22922, 0.55071, 1.26860, 2.38000  &
        ,3.16080, 2.46670, 2.31410, 2.20910, 2.15010/

        DATA(qext(2,7,ib,8),ib=1,17)  & ! Band 7, RH = 87%
        /0.00279, 0.00452, 0.00729, 0.01171, 0.01810, 0.02966  &
        ,0.05125, 0.10113, 0.23774, 0.56974, 1.30050, 2.41400  &
        ,3.15320, 2.43670, 2.31650, 2.20770, 2.14820/

        DATA(qext(2,7,ib,9),ib=1,17)  & ! Band 7, RH = 88%
        /0.00286, 0.00464, 0.00749, 0.01204, 0.01863, 0.03057  &
        ,0.05299, 0.10516, 0.24732, 0.59146, 1.33640, 2.45080  &
        ,3.14320, 2.40660, 2.31750, 2.20550, 2.14630/

        DATA(qext(2,7,ib,10),ib=1,17)  & ! Band 7, RH = 89%
        /0.00294, 0.00478, 0.00771, 0.01240, 0.01923, 0.03160  &
        ,0.05496, 0.10977, 0.25821, 0.61649, 1.37670, 2.49050  &
        ,3.13200, 2.37340, 2.31690, 2.20250, 2.14420/

        DATA(qext(2,7,ib,11),ib=1,17)  & ! Band 7, RH = 90%
        /0.00302, 0.00492, 0.00796, 0.01281, 0.01990, 0.03276  &
        ,0.05723, 0.11513, 0.27075, 0.64563, 1.42240, 2.53400  &
        ,3.11790, 2.34020, 2.31380, 2.19850, 2.14190/

        DATA(qext(2,7,ib,12),ib=1,17)  & ! Band 7, RH = 91%
        /0.00312, 0.00509, 0.00824, 0.01328, 0.02067, 0.03410  &
        ,0.05986, 0.12145, 0.28538, 0.67986, 1.47450, 2.58300  &
        ,3.09590, 2.30790, 2.30720, 2.19340, 2.13910/

        DATA(qext(2,7,ib,13),ib=1,17)  & ! Band 7, RH = 92%
        /0.00323, 0.00529, 0.00857, 0.01383, 0.02156, 0.03566  &
        ,0.06299, 0.12904, 0.30273, 0.72038, 1.53480, 2.63990  &
        ,3.06460, 2.27550, 2.29670, 2.18730, 2.13570/

        DATA(qext(2,7,ib,14),ib=1,17)  & ! Band 7, RH = 93%
        /0.00336, 0.00551, 0.00895, 0.01446, 0.02261, 0.03752  &
        ,0.06678, 0.13839, 0.32377, 0.76877, 1.60660, 2.70490  &
        ,3.02570, 2.24770, 2.28130, 2.18060, 2.13180/

        DATA(qext(2,7,ib,15),ib=1,17)  & ! Band 7, RH = 94%
        /0.00351, 0.00578, 0.00940, 0.01523, 0.02388, 0.03980  &
        ,0.07152, 0.15024, 0.35002, 0.82727, 1.69520, 2.77490  &
        ,2.96670, 2.22730, 2.26110, 2.17410, 2.12750/

        DATA(qext(2,7,ib,16),ib=1,17)  & ! Band 7, RH = 95%
        /0.00369, 0.00610, 0.00997, 0.01619, 0.02547, 0.04269  &
        ,0.07770, 0.16587, 0.38424, 0.89983, 1.80760, 2.85230  &
        ,2.88670, 2.21860, 2.23760, 2.16810, 2.12280/

        DATA(qext(2,7,ib,17),ib=1,17)  & ! Band 7, RH = 96%
        /0.00392, 0.00652, 0.01069, 0.01743, 0.02755, 0.04656  &
        ,0.08625, 0.18771, 0.43198, 0.99545, 1.95140, 2.93730  &
        ,2.76920, 2.22880, 2.21580, 2.16160, 2.11720/

        DATA(qext(2,7,ib,18),ib=1,17)  & ! Band 7, RH = 97%
        /0.00422, 0.00708, 0.01168, 0.01915, 0.03049, 0.05219  &
        ,0.09923, 0.22079, 0.50632, 1.13760, 2.14720, 3.01460  &
        ,2.60090, 2.26170, 2.20290, 2.15090, 2.11010/

        DATA(qext(2,7,ib,19),ib=1,17)  & ! Band 7, RH = 98%
        /0.00464, 0.00791, 0.01320, 0.02183, 0.03520, 0.06163  &
        ,0.12230, 0.27801, 0.64258, 1.37710, 2.43090, 3.03960  &
        ,2.37750, 2.29060, 2.19190, 2.13760, 2.10070/

        DATA(qext(2,7,ib,20),ib=1,17)  & ! Band 7, RH = 99%
        /0.00534, 0.00938, 0.01605, 0.02713, 0.04504, 0.08336  &
        ,0.17931, 0.41017, 0.94219, 1.85680, 2.85380, 2.80770  &
        ,2.22470, 2.22210, 2.16160, 2.11770, 2.08630/

        DATA(qext(2,8,ib,1),ib=1,17)  & ! Band 8, RH = 80%
        /0.00177, 0.00288, 0.00466, 0.00761, 0.01236, 0.02417  &
        ,0.06347, 0.21225, 0.66898, 1.71760, 3.12360, 3.25790  &
        ,2.27840, 2.29480, 2.24140, 2.15630, 2.11720/

        DATA(qext(2,8,ib,2),ib=1,17)  & ! Band 8, RH = 81%
        /0.00182, 0.00294, 0.00476, 0.00777, 0.01264, 0.02484  &
        ,0.06566, 0.21914, 0.68848, 1.74790, 3.14920, 3.22830  &
        ,2.28360, 2.29320, 2.23880, 2.15360, 2.11640/

        DATA(qext(2,8,ib,3),ib=1,17)  & ! Band 8, RH = 82%
        /0.00185, 0.00300, 0.00485, 0.00793, 0.01295, 0.02556  &
        ,0.06803, 0.22651, 0.70921, 1.77990, 3.17580, 3.19910  &
        ,2.28710, 2.29300, 2.23500, 2.15120, 2.11530/

        DATA(qext(2,8,ib,4),ib=1,17)  & ! Band 8, RH = 83%
        /0.00189, 0.00307, 0.00496, 0.00811, 0.01327, 0.02634  &
        ,0.07060, 0.23445, 0.73136, 1.81430, 3.20270, 3.16830  &
        ,2.28770, 2.29420, 2.23000, 2.14930, 2.11380/

        DATA(qext(2,8,ib,5),ib=1,17)  & ! Band 8, RH = 84%
        /0.00193, 0.00313, 0.00507, 0.00830, 0.01361, 0.02718  &
        ,0.07342, 0.24308, 0.75519, 1.85160, 3.22930, 3.13290  &
        ,2.29180, 2.29520, 2.22300, 2.14840, 2.11200/

        DATA(qext(2,8,ib,6),ib=1,17)  & ! Band 8, RH = 85%
        /0.00197, 0.00320, 0.00518, 0.00850, 0.01399, 0.02810  &
        ,0.07655, 0.25256, 0.78103, 1.89270, 3.25500, 3.09140  &
        ,2.30210, 2.29600, 2.21520, 2.14840, 2.11020/

        DATA(qext(2,8,ib,7),ib=1,17)  & ! Band 8, RH = 86%
        /0.00202, 0.00328, 0.00531, 0.00871, 0.01440, 0.02913  &
        ,0.08006, 0.26308, 0.80925, 1.93820, 3.27980, 3.04470  &
        ,2.31160, 2.29960, 2.20630, 2.14920, 2.10870/

        DATA(qext(2,8,ib,8),ib=1,17)  & ! Band 8, RH = 87%
        /0.00207, 0.00336, 0.00544, 0.00895, 0.01485, 0.03029  &
        ,0.08405, 0.27490, 0.84036, 1.98920, 3.30450, 2.99600  &
        ,2.31780, 2.29990, 2.19720, 2.15010, 2.10740/

        DATA(qext(2,8,ib,9),ib=1,17)  & ! Band 8, RH = 88%
        /0.00212, 0.00345, 0.00560, 0.00921, 0.01535, 0.03161  &
        ,0.08865, 0.28834, 0.87498, 2.04630, 3.33100, 2.94370  &
        ,2.32970, 2.30200, 2.18880, 2.15000, 2.10600/

        DATA(qext(2,8,ib,10),ib=1,17)  & ! Band 8, RH = 89%
        /0.00218, 0.00354, 0.00576, 0.00951, 0.01592, 0.03313  &
        ,0.09402, 0.30384, 0.91397, 2.10990, 3.36090, 2.88030  &
        ,2.34080, 2.30070, 2.18270, 2.14770, 2.10440/

        DATA(qext(2,8,ib,11),ib=1,17)  & ! Band 8, RH = 90%
        /0.00224, 0.00365, 0.00595, 0.00984, 0.01658, 0.03493  &
        ,0.10038, 0.32200, 0.95855, 2.18040, 3.39140, 2.80790  &
        ,2.35260, 2.29680, 2.18000, 2.14270, 2.10270/

        DATA(qext(2,8,ib,12),ib=1,17)  & ! Band 8, RH = 91%
        /0.00231, 0.00378, 0.00617, 0.01023, 0.01735, 0.03707  &
        ,0.10805, 0.34366, 1.01060, 2.25850, 3.41560, 2.73660  &
        ,2.35920, 2.28670, 2.18150, 2.13670, 2.10110/

        DATA(qext(2,8,ib,13),ib=1,17)  & ! Band 8, RH = 92%
        /0.00239, 0.00392, 0.00641, 0.01068, 0.01826, 0.03970  &
        ,0.11747, 0.37015, 1.07290, 2.34660, 3.42990, 2.64860  &
        ,2.36740, 2.27010, 2.18470, 2.13250, 2.09850/

        DATA(qext(2,8,ib,14),ib=1,17)  & ! Band 8, RH = 93%
        /0.00249, 0.00409, 0.00671, 0.01121, 0.01937, 0.04300  &
        ,0.12934, 0.40351, 1.14980, 2.45100, 3.44240, 2.55890  &
        ,2.36480, 2.24720, 2.18500, 2.13050, 2.09550/

        DATA(qext(2,8,ib,15),ib=1,17)  & ! Band 8, RH = 94%
        /0.00260, 0.00429, 0.00706, 0.01186, 0.02077, 0.04730  &
        ,0.14474, 0.44726, 1.24790, 2.57970, 3.44160, 2.46190  &
        ,2.35110, 2.22250, 2.17720, 2.12690, 2.09300/

        DATA(qext(2,8,ib,16),ib=1,17)  & ! Band 8, RH = 95%
        /0.00273, 0.00453, 0.00749, 0.01269, 0.02260, 0.05319  &
        ,0.16556, 0.50753, 1.37520, 2.73180, 3.40190, 2.37440  &
        ,2.32520, 2.20960, 2.16710, 2.12280, 2.08970/

        DATA(qext(2,8,ib,17),ib=1,17)  & ! Band 8, RH = 96%
        /0.00290, 0.00484, 0.00806, 0.01379, 0.02516, 0.06182  &
        ,0.19530, 0.59537, 1.54120, 2.90450, 3.32350, 2.30460  &
        ,2.29590, 2.22100, 2.16390, 2.11660, 2.08560/

        DATA(qext(2,8,ib,18),ib=1,17)  & ! Band 8, RH = 97%
        /0.00312, 0.00526, 0.00884, 0.01538, 0.02911, 0.07586  &
        ,0.24142, 0.73009, 1.77310, 3.12150, 3.14720, 2.29110  &
        ,2.28460, 2.22470, 2.14560, 2.11030, 2.08040/

        DATA(qext(2,8,ib,19),ib=1,17)  & ! Band 8, RH = 98%
        /0.00344, 0.00589, 0.01006, 0.01804, 0.03634, 0.10299  &
        ,0.32442, 0.95126, 2.14830, 3.34280, 2.80740, 2.34140  &
        ,2.29160, 2.17850, 2.13990, 2.10020, 2.07370/

        DATA(qext(2,8,ib,20),ib=1,17)  & ! Band 8, RH = 99%
        /0.00396, 0.00701, 0.01245, 0.02405, 0.05555, 0.17662  &
        ,0.54577, 1.45180, 2.80730, 3.34810, 2.33220, 2.30030  &
        ,2.21140, 2.16340, 2.11850, 2.08590, 2.06320/

! Qext for Mineral Dust:
        DATA(qext(3,1,ib,1),ib=1,17)  & ! Band 1, All RH's
        /0.00124, 0.00198, 0.00323, 0.00566, 0.01209, 0.03703  &
        ,0.13977, 0.46990, 1.22290, 2.26570, 2.61390, 2.64980  &
        ,2.54470, 2.19980, 2.20590, 2.13130, 2.09540/

        DATA(qext(3,2,ib,1),ib=1,17)  & ! Band 2, All RH's
        /0.00224, 0.00409, 0.00974, 0.03551, 0.16816, 0.69694  &
        ,2.02260, 3.42610, 3.04780, 2.63640, 2.33260, 2.27750  &
        ,2.19670, 2.14290, 2.10390, 2.07640, 2.05630/

        DATA(qext(3,3,ib,1),ib=1,17)  & ! Band 3, All RH's
        /0.00599, 0.01846, 0.08179, 0.36834, 1.28590, 2.91150  &
        ,3.72340, 2.62500, 2.42800, 2.26110, 2.21600, 2.16130  &
        ,2.12100, 2.09050, 2.06620, 2.04880, 2.03590/

        DATA(qext(3,4,ib,1),ib=1,17)  & ! Band 4, All RH's
        /0.00232, 0.00368, 0.00584, 0.00927, 0.01472, 0.02341  &
        ,0.03733, 0.05997, 0.09838, 0.17106, 0.34448, 0.88581  &
        ,2.31920, 3.26570, 2.89620, 2.60710, 2.45250/

        DATA(qext(3,5,ib,1),ib=1,17)  & ! Band 5, All RH's
        /0.00154, 0.00240, 0.00387, 0.00615, 0.00978, 0.01558  &
        ,0.02498, 0.04083, 0.07056, 0.14147, 0.36788, 1.06280  &
        ,2.46800, 3.53330, 2.71720, 2.50290, 2.33060/

        DATA(qext(3,6,ib,1),ib=1,17)  & ! Band 6, All RH's
        /0.00319, 0.00506, 0.00804, 0.01277, 0.02033, 0.03253  &
        ,0.05279, 0.08919, 0.16752, 0.38114, 0.93498, 1.87920  &
        ,2.68110, 2.74410, 2.47570, 2.30710, 2.22650/

        DATA(qext(3,7,ib,1),ib=1,17)  & ! Band 7, All RH's
        /0.00347, 0.00550, 0.00874, 0.01389, 0.02215, 0.03557  &
        ,0.05837, 0.10191, 0.20923, 0.57176, 1.73950, 3.22860  &
        ,3.29090, 2.47200, 2.45570, 2.32820, 2.24460/

        DATA(qext(3,8,ib,1),ib=1,17)  & ! Band 8, All RH's
        /0.00147, 0.00234, 0.00371, 0.00591, 0.00942, 0.01519  &
        ,0.02535, 0.04659, 0.10365, 0.25753, 0.65195, 1.44010  &
        ,2.51080, 2.84480, 2.33600, 2.25070, 2.18000/

! Returning value for qext:
        value = qext(aerotype,radband,bin,rh)

return
END SUBROUTINE aeroqext

!##############################################################################
Subroutine aeroqscat (aerotype,radband,bin,rh,value)

implicit none

        INTEGER aerotype,radband,bin,rh,ib
        REAL qscat(3,8,17,20),value

! Setting for aerosol with no growth (completely insoluble):
        IF(AEROTYPE.EQ.3) rh = 1

! Qscat for Ammonium Sulfate (with a 90% insoluble dust inclusion):
        DATA(qscat(1,1,ib,1),ib=1,17)  & ! Band 1, RH = 80%
        /0.00000, 0.00002, 0.00013, 0.00080, 0.00454, 0.02746  &
        ,0.13553, 0.49914, 1.31130, 2.38440, 2.54700, 2.37630  &
        ,1.94730, 1.59790, 1.42910, 1.27210, 1.16110/

        DATA(qscat(1,1,ib,2),ib=1,17)  & ! Band 1, RH = 81%
        /0.00000, 0.00002, 0.00013, 0.00081, 0.00459, 0.02775  &
        ,0.13664, 0.50267, 1.31790, 2.39340, 2.55060, 2.37430  &
        ,1.94120, 1.60120, 1.42960, 1.27260, 1.16180/

        DATA(qscat(1,1,ib,3),ib=1,17)  & ! Band 1, RH = 82%
        /0.00000, 0.00002, 0.00013, 0.00082, 0.00465, 0.02807  &
        ,0.13788, 0.50658, 1.32520, 2.40310, 2.55440, 2.37240  &
        ,1.93570, 1.60310, 1.43000, 1.27350, 1.16250/

        DATA(qscat(1,1,ib,4),ib=1,17)  & ! Band 1, RH = 83%
        /0.00000, 0.00002, 0.00013, 0.00083, 0.00471, 0.02843  &
        ,0.13925, 0.51092, 1.33340, 2.41380, 2.55820, 2.37100  &
        ,1.92910, 1.60270, 1.42970, 1.27440, 1.16330/

        DATA(qscat(1,1,ib,5),ib=1,17)  & ! Band 1, RH = 84%
        /0.00000, 0.00002, 0.00013, 0.00084, 0.00478, 0.02884  &
        ,0.14079, 0.51578, 1.34270, 2.42560, 2.56200, 2.37010  &
        ,1.91990, 1.60110, 1.43010, 1.27550, 1.16420/

        DATA(qscat(1,1,ib,6),ib=1,17)  & ! Band 1, RH = 85%
        /0.00000, 0.00002, 0.00013, 0.00086, 0.00486, 0.02930  &
        ,0.14253, 0.52125, 1.35310, 2.43850, 2.56570, 2.36890  &
        ,1.90880, 1.60400, 1.43220, 1.27660, 1.16510/

        DATA(qscat(1,1,ib,7),ib=1,17)  & ! Band 1, RH = 86%
        /0.00000, 0.00002, 0.00014, 0.00087, 0.00495, 0.02982  &
        ,0.14452, 0.52746, 1.36510, 2.45300, 2.56910, 2.36640  &
        ,1.89930, 1.61320, 1.43390, 1.27800, 1.16620/

        DATA(qscat(1,1,ib,8),ib=1,17)  & ! Band 1, RH = 87%
        /0.00000, 0.00002, 0.00014, 0.00089, 0.00506, 0.03043  &
        ,0.14679, 0.53456, 1.37880, 2.46900, 2.57220, 2.36120  &
        ,1.89220, 1.62010, 1.43540, 1.27940, 1.16740/

        DATA(qscat(1,1,ib,9),ib=1,17)  & ! Band 1, RH = 88%
        /0.00000, 0.00002, 0.00014, 0.00091, 0.00518, 0.03113  &
        ,0.14944, 0.54276, 1.39490, 2.48710, 2.57490, 2.35660  &
        ,1.88130, 1.61800, 1.43910, 1.28130, 1.16890/

        DATA(qscat(1,1,ib,10),ib=1,17)  & ! Band 1, RH = 89%
        /0.00000, 0.00002, 0.00015, 0.00093, 0.00532, 0.03197  &
        ,0.15255, 0.55234, 1.41380, 2.50770, 2.57720, 2.35330  &
        ,1.87260, 1.62400, 1.44170, 1.28330, 1.17070/

        DATA(qscat(1,1,ib,11),ib=1,17)  & ! Band 1, RH = 90%
        /0.00000, 0.00002, 0.00015, 0.00096, 0.00549, 0.03297  &
        ,0.15625, 0.56369, 1.43640, 2.53130, 2.57940, 2.34380  &
        ,1.86270, 1.63620, 1.44590, 1.28640, 1.17310/

        DATA(qscat(1,1,ib,12),ib=1,17)  & ! Band 1, RH = 91%
        /0.00000, 0.00002, 0.00015, 0.00100, 0.00570, 0.03420  &
        ,0.16075, 0.57732, 1.46390, 2.55900, 2.58290, 2.33190  &
        ,1.85060, 1.63320, 1.44980, 1.29070, 1.17640/

        DATA(qscat(1,1,ib,13),ib=1,17)  & ! Band 1, RH = 92%
        /0.00000, 0.00002, 0.00016, 0.00104, 0.00597, 0.03573  &
        ,0.16631, 0.59404, 1.49790, 2.59240, 2.58770, 2.32240  &
        ,1.84220, 1.64720, 1.45460, 1.29700, 1.18070/

        DATA(qscat(1,1,ib,14),ib=1,17)  & ! Band 1, RH = 93%
        /0.00000, 0.00003, 0.00017, 0.00110, 0.00630, 0.03769  &
        ,0.17340, 0.61501, 1.54080, 2.63430, 2.58920, 2.30550  &
        ,1.83310, 1.65330, 1.45950, 1.30520, 1.18570/

        DATA(qscat(1,1,ib,15),ib=1,17)  & ! Band 1, RH = 94%
        /0.00000, 0.00003, 0.00018, 0.00117, 0.00675, 0.04031  &
        ,0.18273, 0.64219, 1.59610, 2.68860, 2.58160, 2.28040  &
        ,1.82810, 1.65800, 1.46740, 1.31530, 1.19150/

        DATA(qscat(1,1,ib,16),ib=1,17)  & ! Band 1, RH = 95%
        /0.00000, 0.00003, 0.00019, 0.00128, 0.00738, 0.04396  &
        ,0.19560, 0.67898, 1.66920, 2.75720, 2.56080, 2.25160  &
        ,1.82730, 1.66620, 1.47870, 1.32690, 1.19970/

        DATA(qscat(1,1,ib,17),ib=1,17)  & ! Band 1, RH = 96%
        /0.00000, 0.00003, 0.00021, 0.00142, 0.00832, 0.04940  &
        ,0.21458, 0.73217, 1.76990, 2.83560, 2.53460, 2.21090  &
        ,1.82830, 1.67110, 1.49800, 1.34040, 1.21270/

        DATA(qscat(1,1,ib,18),ib=1,17)  & ! Band 1, RH = 97%
        /0.00000, 0.00003, 0.00024, 0.00166, 0.00987, 0.05834  &
        ,0.24568, 0.81783, 1.92050, 2.92820, 2.46480, 2.13970  &
        ,1.82830, 1.67450, 1.51900, 1.36490, 1.23070/

        DATA(qscat(1,1,ib,19),ib=1,17)  & ! Band 1, RH = 98%
        /0.00000, 0.00004, 0.00029, 0.00211, 0.01290, 0.07563  &
        ,0.30715, 0.98207, 2.18580, 3.02680, 2.34250, 2.02900  &
        ,1.82970, 1.69690, 1.55040, 1.39930, 1.26190/

        DATA(qscat(1,1,ib,20),ib=1,17)  & ! Band 1, RH = 99%
        /0.00001, 0.00005, 0.00040, 0.00323, 0.02131, 0.12172  &
        ,0.48233, 1.37510, 2.69820, 2.95220, 2.19410, 1.92460  &
        ,1.86720, 1.71550, 1.59730, 1.45850, 1.32080/

        DATA(qscat(1,2,ib,1),ib=1,17)  & ! Band 2, RH = 80%
        /0.00011, 0.00071, 0.00458, 0.02897, 0.15655, 0.67540  &
        ,1.95970, 3.31950, 2.88180, 2.33520, 1.97270, 1.79570  &
        ,1.59220, 1.41060, 1.26120, 1.16840, 1.12650/

        DATA(qscat(1,2,ib,2),ib=1,17)  & ! Band 2, RH = 81%
        /0.00011, 0.00072, 0.00461, 0.02916, 0.15731, 0.67728  &
        ,1.96030, 3.31930, 2.88280, 2.33570, 1.97470, 1.79770  &
        ,1.59440, 1.41280, 1.26260, 1.16900, 1.12650/

        DATA(qscat(1,2,ib,3),ib=1,17)  & ! Band 2, RH = 82%
        /0.00011, 0.00072, 0.00464, 0.02938, 0.15814, 0.67936  &
        ,1.96110, 3.31920, 2.88370, 2.33620, 1.97680, 1.79970  &
        ,1.59660, 1.41510, 1.26430, 1.16980, 1.12650/

        DATA(qscat(1,2,ib,4),ib=1,17)  & ! Band 2, RH = 83%
        /0.00011, 0.00073, 0.00468, 0.02962, 0.15906, 0.68169  &
        ,1.96200, 3.31930, 2.88460, 2.33660, 1.97900, 1.80160  &
        ,1.59910, 1.41770, 1.26620, 1.17050, 1.12650/

        DATA(qscat(1,2,ib,5),ib=1,17)  & ! Band 2, RH = 84%
        /0.00011, 0.00074, 0.00472, 0.02988, 0.16010, 0.68431  &
        ,1.96310, 3.31940, 2.88530, 2.33720, 1.98170, 1.80390  &
        ,1.60210, 1.42050, 1.26820, 1.17150, 1.12660/

        DATA(qscat(1,2,ib,6),ib=1,17)  & ! Band 2, RH = 85%
        /0.00012, 0.00074, 0.00477, 0.03019, 0.16127, 0.68728  &
        ,1.96450, 3.31990, 2.88570, 2.33740, 1.98440, 1.80680  &
        ,1.60540, 1.42370, 1.27070, 1.17240, 1.12660/

        DATA(qscat(1,2,ib,7),ib=1,17)  & ! Band 2, RH = 86%
        /0.00012, 0.00075, 0.00482, 0.03053, 0.16259, 0.69068  &
        ,1.96620, 3.32070, 2.88580, 2.33750, 1.98770, 1.80940  &
        ,1.60920, 1.42720, 1.27330, 1.17370, 1.12670/

        DATA(qscat(1,2,ib,8),ib=1,17)  & ! Band 2, RH = 87%
        /0.00012, 0.00076, 0.00488, 0.03093, 0.16410, 0.69460  &
        ,1.96830, 3.32190, 2.88540, 2.33770, 1.99110, 1.81280  &
        ,1.61340, 1.43150, 1.27640, 1.17510, 1.12680/

        DATA(qscat(1,2,ib,9),ib=1,17)  & ! Band 2, RH = 88%
        /0.00012, 0.00077, 0.00496, 0.03138, 0.16586, 0.69916  &
        ,1.97090, 3.32390, 2.88440, 2.33750, 1.99590, 1.81580  &
        ,1.61780, 1.43610, 1.28010, 1.17670, 1.12700/

        DATA(qscat(1,2,ib,10),ib=1,17)  & ! Band 2, RH = 89%
        /0.00012, 0.00078, 0.00504, 0.03192, 0.16791, 0.70454  &
        ,1.97430, 3.32680, 2.88300, 2.33790, 2.00150, 1.82010  &
        ,1.62330, 1.44160, 1.28450, 1.17870, 1.12730/

        DATA(qscat(1,2,ib,11),ib=1,17)  & ! Band 2, RH = 90%
        /0.00012, 0.00079, 0.00514, 0.03255, 0.17034, 0.71097  &
        ,1.97880, 3.33110, 2.88130, 2.33730, 2.00740, 1.82540  &
        ,1.62980, 1.44780, 1.28950, 1.18120, 1.12770/

        DATA(qscat(1,2,ib,12),ib=1,17)  & ! Band 2, RH = 91%
        /0.00012, 0.00081, 0.00525, 0.03332, 0.17327, 0.71878  &
        ,1.98480, 3.33690, 2.88030, 2.33650, 2.01350, 1.83110  &
        ,1.63780, 1.45510, 1.29570, 1.18440, 1.12830/

        DATA(qscat(1,2,ib,13),ib=1,17)  & ! Band 2, RH = 92%
        /0.00013, 0.00083, 0.00539, 0.03426, 0.17687, 0.72845  &
        ,1.99310, 3.34460, 2.87970, 2.33620, 2.02040, 1.83790  &
        ,1.64670, 1.46420, 1.30350, 1.18820, 1.12910/

        DATA(qscat(1,2,ib,14),ib=1,17)  & ! Band 2, RH = 93%
        /0.00013, 0.00085, 0.00557, 0.03545, 0.18141, 0.74072  &
        ,2.00510, 3.35410, 2.87720, 2.33430, 2.02940, 1.84500  &
        ,1.65830, 1.47550, 1.31260, 1.19350, 1.13040/

        DATA(qscat(1,2,ib,15),ib=1,17)  & ! Band 2, RH = 94%
        /0.00013, 0.00088, 0.00580, 0.03701, 0.18731, 0.75676  &
        ,2.02300, 3.36520, 2.87100, 2.33060, 2.04090, 1.85500  &
        ,1.67270, 1.48950, 1.32480, 1.20060, 1.13250/

        DATA(qscat(1,2,ib,16),ib=1,17)  & ! Band 2, RH = 95%
        /0.00014, 0.00092, 0.00611, 0.03912, 0.19529, 0.77856  &
        ,2.05140, 3.37710, 2.85500, 2.32270, 2.05760, 1.86700  &
        ,1.69140, 1.50760, 1.34150, 1.21130, 1.13620/

        DATA(qscat(1,2,ib,17),ib=1,17)  & ! Band 2, RH = 96%
        /0.00014, 0.00097, 0.00654, 0.04216, 0.20674, 0.80984  &
        ,2.09870, 3.39130, 2.82950, 2.31320, 2.07370, 1.88270  &
        ,1.71700, 1.53410, 1.36470, 1.22720, 1.14240/

        DATA(qscat(1,2,ib,18),ib=1,17)  & ! Band 2, RH = 97%
        /0.00015, 0.00104, 0.00719, 0.04693, 0.22465, 0.85868  &
        ,2.18040, 3.42460, 2.78720, 2.29430, 2.09670, 1.90580  &
        ,1.75050, 1.57290, 1.40170, 1.25330, 1.15470/

        DATA(qscat(1,2,ib,19),ib=1,17)  & ! Band 2, RH = 98%
        /0.00016, 0.00115, 0.00830, 0.05550, 0.25710, 0.94827  &
        ,2.32580, 3.46800, 2.69020, 2.24680, 2.12630, 1.94210  &
        ,1.80220, 1.63660, 1.46260, 1.30290, 1.18290/

        DATA(qscat(1,2,ib,20),ib=1,17)  & ! Band 2, RH = 99%
        /0.00017, 0.00135, 0.01067, 0.07559, 0.33861, 1.18810  &
        ,2.68340, 3.47430, 2.49650, 2.15760, 2.15960, 2.01730  &
        ,1.87270, 1.74060, 1.58500, 1.42040, 1.26840/

        DATA(qscat(1,3,ib,1),ib=1,17)  & ! Band 3, RH = 80%
        /0.00202, 0.01296, 0.07746, 0.36935, 1.25280, 2.81310  &
        ,3.56840, 2.40800, 2.09700, 1.85010, 1.68620, 1.50270  &
        ,1.33810, 1.21600, 1.14800, 1.12260, 1.11420/

        DATA(qscat(1,3,ib,2),ib=1,17)  & ! Band 3, RH = 81%
        /0.00203, 0.01305, 0.07791, 0.37075, 1.25460, 2.81260  &
        ,3.56790, 2.40940, 2.09830, 1.85140, 1.68820, 1.50500  &
        ,1.34000, 1.21710, 1.14830, 1.12240, 1.11390/

        DATA(qscat(1,3,ib,3),ib=1,17)  & ! Band 3, RH = 82%
        /0.00204, 0.01314, 0.07841, 0.37230, 1.25670, 2.81200  &
        ,3.56730, 2.41030, 2.09990, 1.85340, 1.69050, 1.50740  &
        ,1.34210, 1.21830, 1.14860, 1.12220, 1.11360/

        DATA(qscat(1,3,ib,4),ib=1,17)  & ! Band 3, RH = 83%
        /0.00206, 0.01324, 0.07896, 0.37402, 1.25900, 2.81160  &
        ,3.56670, 2.41100, 2.10140, 1.85550, 1.69310, 1.51030  &
        ,1.34450, 1.21970, 1.14890, 1.12200, 1.11330/

        DATA(qscat(1,3,ib,5),ib=1,17)  & ! Band 3, RH = 84%
        /0.00207, 0.01336, 0.07958, 0.37594, 1.26160, 2.81110  &
        ,3.56620, 2.41190, 2.10310, 1.85730, 1.69600, 1.51340  &
        ,1.34690, 1.22130, 1.14930, 1.12170, 1.11290/

        DATA(qscat(1,3,ib,6),ib=1,17)  & ! Band 3, RH = 85%
        /0.00209, 0.01348, 0.08027, 0.37810, 1.26450, 2.81080  &
        ,3.56580, 2.41280, 2.10580, 1.85970, 1.69940, 1.51670  &
        ,1.34980, 1.22300, 1.14980, 1.12150, 1.11250/

        DATA(qscat(1,3,ib,7),ib=1,17)  & ! Band 3, RH = 86%
        /0.00211, 0.01363, 0.08105, 0.38055, 1.26800, 2.81050  &
        ,3.56530, 2.41390, 2.10840, 1.86230, 1.70310, 1.52050  &
        ,1.35310, 1.22500, 1.15040, 1.12130, 1.11200/

        DATA(qscat(1,3,ib,8),ib=1,17)  & ! Band 3, RH = 87%
        /0.00213, 0.01379, 0.08193, 0.38335, 1.27200, 2.81050  &
        ,3.56370, 2.41440, 2.11160, 1.86550, 1.70750, 1.52480  &
        ,1.35670, 1.22720, 1.15100, 1.12100, 1.11150/

        DATA(qscat(1,3,ib,9),ib=1,17)  & ! Band 3, RH = 88%
        /0.00216, 0.01397, 0.08295, 0.38657, 1.27660, 2.81070  &
        ,3.56160, 2.41510, 2.11540, 1.86900, 1.71210, 1.52990  &
        ,1.36120, 1.23000, 1.15190, 1.12080, 1.11090/

        DATA(qscat(1,3,ib,10),ib=1,17)  & ! Band 3, RH = 89%
        /0.00219, 0.01418, 0.08412, 0.39034, 1.28220, 2.81140  &
        ,3.56020, 2.41410, 2.12040, 1.87270, 1.71740, 1.53530  &
        ,1.36620, 1.23330, 1.15300, 1.12050, 1.11030/

        DATA(qscat(1,3,ib,11),ib=1,17)  & ! Band 3, RH = 90%
        /0.00222, 0.01443, 0.08551, 0.39479, 1.28900, 2.81280  &
        ,3.55860, 2.41240, 2.12550, 1.87770, 1.72340, 1.54170  &
        ,1.37210, 1.23720, 1.15440, 1.12020, 1.10960/

        DATA(qscat(1,3,ib,12),ib=1,17)  & ! Band 3, RH = 91%
        /0.00225, 0.01472, 0.08715, 0.40014, 1.29720, 2.81530  &
        ,3.55690, 2.40990, 2.12980, 1.88390, 1.73120, 1.55000  &
        ,1.37920, 1.24210, 1.15610, 1.12010, 1.10870/

        DATA(qscat(1,3,ib,13),ib=1,17)  & ! Band 3, RH = 92%
        /0.00230, 0.01507, 0.08915, 0.40670, 1.30740, 2.81950  &
        ,3.55480, 2.40560, 2.13390, 1.88930, 1.74070, 1.55990  &
        ,1.38780, 1.24790, 1.15850, 1.12000, 1.10780/

        DATA(qscat(1,3,ib,14),ib=1,17)  & ! Band 3, RH = 93%
        /0.00235, 0.01550, 0.09162, 0.41494, 1.32050, 2.82600  &
        ,3.55210, 2.40080, 2.13970, 1.89680, 1.75090, 1.57220  &
        ,1.39900, 1.25600, 1.16170, 1.12000, 1.10670/

        DATA(qscat(1,3,ib,15),ib=1,17)  & ! Band 3, RH = 94%
        /0.00241, 0.01603, 0.09476, 0.42561, 1.33780, 2.83630  &
        ,3.54830, 2.39230, 2.14830, 1.90550, 1.76540, 1.58680  &
        ,1.41270, 1.26620, 1.16650, 1.12050, 1.10550/

        DATA(qscat(1,3,ib,16),ib=1,17)  & ! Band 3, RH = 95%
        /0.00249, 0.01672, 0.09891, 0.43996, 1.36130, 2.85320  &
        ,3.54130, 2.38190, 2.16300, 1.91810, 1.78270, 1.60670  &
        ,1.43190, 1.28040, 1.17360, 1.12160, 1.10420/

        DATA(qscat(1,3,ib,17),ib=1,17)  & ! Band 3, RH = 96%
        /0.00259, 0.01764, 0.10465, 0.46034, 1.39550, 2.88250  &
        ,3.52370, 2.36790, 2.18340, 1.93290, 1.80760, 1.63290  &
        ,1.45770, 1.30030, 1.18540, 1.12450, 1.10290/

        DATA(qscat(1,3,ib,18),ib=1,17)  & ! Band 3, RH = 97%
        /0.00272, 0.01894, 0.11316, 0.49152, 1.45080, 2.93570  &
        ,3.48900, 2.32910, 2.20870, 1.95730, 1.84510, 1.67120  &
        ,1.49660, 1.33330, 1.20570, 1.13130, 1.10220/

        DATA(qscat(1,3,ib,19),ib=1,17)  & ! Band 3, RH = 98%
        /0.00290, 0.02096, 0.12723, 0.54491, 1.55650, 3.03920  &
        ,3.42950, 2.25740, 2.24620, 2.00390, 1.89780, 1.73310  &
        ,1.56050, 1.39220, 1.24710, 1.15000, 1.10490/

        DATA(qscat(1,3,ib,20),ib=1,17)  & ! Band 3, RH = 99%
        /0.00317, 0.02451, 0.15542, 0.65754, 1.80500, 3.26390  &
        ,3.20780, 2.14190, 2.27970, 2.09880, 1.94200, 1.81610  &
        ,1.67600, 1.51750, 1.35540, 1.21840, 1.13060/

        DATA(qscat(1,4,ib,1),ib=1,17)  & ! Band 4, RH = 80%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00001  &
        ,0.00003, 0.00020, 0.00131, 0.00844, 0.05454, 0.31563  &
        ,1.12620, 1.76130, 1.40620, 1.22980, 1.23050/

        DATA(qscat(1,4,ib,2),ib=1,17)  & ! Band 4, RH = 81%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00001  &
        ,0.00003, 0.00021, 0.00132, 0.00853, 0.05505, 0.31755  &
        ,1.12530, 1.75470, 1.40390, 1.22910, 1.23000/

        DATA(qscat(1,4,ib,3),ib=1,17)  & ! Band 4, RH = 82%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00001  &
        ,0.00003, 0.00021, 0.00134, 0.00862, 0.05562, 0.31968  &
        ,1.12420, 1.74750, 1.40150, 1.22830, 1.22940/

        DATA(qscat(1,4,ib,4),ib=1,17)  & ! Band 4, RH = 83%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00001  &
        ,0.00003, 0.00021, 0.00135, 0.00873, 0.05626, 0.32204  &
        ,1.12310, 1.73970, 1.39880, 1.22750, 1.22880/

        DATA(qscat(1,4,ib,5),ib=1,17)  & ! Band 4, RH = 84%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00001  &
        ,0.00003, 0.00021, 0.00137, 0.00884, 0.05697, 0.32468  &
        ,1.12200, 1.73110, 1.39580, 1.22670, 1.22810/

        DATA(qscat(1,4,ib,6),ib=1,17)  & ! Band 4, RH = 85%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00001  &
        ,0.00003, 0.00022, 0.00139, 0.00898, 0.05779, 0.32766  &
        ,1.12070, 1.72180, 1.39250, 1.22580, 1.22730/

        DATA(qscat(1,4,ib,7),ib=1,17)  & ! Band 4, RH = 86%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00001  &
        ,0.00003, 0.00022, 0.00142, 0.00914, 0.05872, 0.33104  &
        ,1.11930, 1.71160, 1.38870, 1.22480, 1.22650/

        DATA(qscat(1,4,ib,8),ib=1,17)  & ! Band 4, RH = 87%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00001  &
        ,0.00004, 0.00023, 0.00145, 0.00932, 0.05981, 0.33492  &
        ,1.11790, 1.70030, 1.38460, 1.22370, 1.22550/

        DATA(qscat(1,4,ib,9),ib=1,17)  & ! Band 4, RH = 88%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00001  &
        ,0.00004, 0.00023, 0.00148, 0.00953, 0.06108, 0.33940  &
        ,1.11630, 1.68780, 1.37990, 1.22250, 1.22450/

        DATA(qscat(1,4,ib,10),ib=1,17)  & ! Band 4, RH = 89%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00001  &
        ,0.00004, 0.00024, 0.00152, 0.00978, 0.06259, 0.34463  &
        ,1.11470, 1.67390, 1.37450, 1.22130, 1.22330/

        DATA(qscat(1,4,ib,11),ib=1,17)  & ! Band 4, RH = 90%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00001  &
        ,0.00004, 0.00024, 0.00157, 0.01009, 0.06441, 0.35084  &
        ,1.11300, 1.65840, 1.36840, 1.21980, 1.22190/

        DATA(qscat(1,4,ib,12),ib=1,17)  & ! Band 4, RH = 91%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00001  &
        ,0.00004, 0.00025, 0.00162, 0.01047, 0.06665, 0.35832  &
        ,1.11130, 1.64090, 1.36130, 1.21820, 1.22030/

        DATA(qscat(1,4,ib,13),ib=1,17)  & ! Band 4, RH = 92%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00001  &
        ,0.00004, 0.00026, 0.00170, 0.01095, 0.06947, 0.36749  &
        ,1.10970, 1.62130, 1.35310, 1.21640, 1.21850/

        DATA(qscat(1,4,ib,14),ib=1,17)  & ! Band 4, RH = 93%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00001  &
        ,0.00004, 0.00028, 0.00180, 0.01157, 0.07314, 0.37902  &
        ,1.10830, 1.59900, 1.34330, 1.21430, 1.21630/

        DATA(qscat(1,4,ib,15),ib=1,17)  & ! Band 4, RH = 94%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00001  &
        ,0.00005, 0.00030, 0.00193, 0.01243, 0.07808, 0.39393  &
        ,1.10770, 1.57350, 1.33160, 1.21180, 1.21380/

        DATA(qscat(1,4,ib,16),ib=1,17)  & ! Band 4, RH = 95%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00001  &
        ,0.00005, 0.00033, 0.00212, 0.01366, 0.08508, 0.41397  &
        ,1.10860, 1.54430, 1.31710, 1.20890, 1.21060/

        DATA(qscat(1,4,ib,17),ib=1,17)  & ! Band 4, RH = 96%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00001  &
        ,0.00006, 0.00038, 0.00242, 0.01556, 0.09576, 0.44237  &
        ,1.11300, 1.51030, 1.29860, 1.20530, 1.20680/

        DATA(qscat(1,4,ib,18),ib=1,17)  & ! Band 4, RH = 97%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00001  &
        ,0.00007, 0.00046, 0.00295, 0.01891, 0.11385, 0.48578  &
        ,1.12530, 1.46890, 1.27450, 1.20070, 1.20200/

        DATA(qscat(1,4,ib,19),ib=1,17)  & ! Band 4, RH = 98%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00001  &
        ,0.00010, 0.00063, 0.00409, 0.02614, 0.15040, 0.56095  &
        ,1.15670, 1.41490, 1.24270, 1.19470, 1.19590/

        DATA(qscat(1,4,ib,20),ib=1,17)  & ! Band 4, RH = 99%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00003  &
        ,0.00018, 0.00121, 0.00804, 0.05065, 0.25404, 0.72700  &
        ,1.23010, 1.33440, 1.20190, 1.18720, 1.18810/

        DATA(qscat(1,5,ib,1),ib=1,17)  & ! Band 5, RH = 80%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00002  &
        ,0.00012, 0.00074, 0.00470, 0.02958, 0.16912, 0.66331  &
        ,1.65090, 2.22490, 1.39120, 1.22580, 1.17470/

        DATA(qscat(1,5,ib,2),ib=1,17)  & ! Band 5, RH = 81%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00002  &
        ,0.00012, 0.00074, 0.00473, 0.02977, 0.16963, 0.66110  &
        ,1.63450, 2.19630, 1.38460, 1.21940, 1.17280/

        DATA(qscat(1,5,ib,3),ib=1,17)  & ! Band 5, RH = 82%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00002  &
        ,0.00012, 0.00075, 0.00477, 0.02998, 0.17021, 0.65875  &
        ,1.61700, 2.16600, 1.37770, 1.21300, 1.17110/

        DATA(qscat(1,5,ib,4),ib=1,17)  & ! Band 5, RH = 83%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00002  &
        ,0.00012, 0.00076, 0.00481, 0.03021, 0.17087, 0.65624  &
        ,1.59820, 2.13360, 1.37040, 1.20660, 1.16930/

        DATA(qscat(1,5,ib,5),ib=1,17)  & ! Band 5, RH = 84%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00002  &
        ,0.00012, 0.00076, 0.00486, 0.03048, 0.17163, 0.65358  &
        ,1.57810, 2.09910, 1.36260, 1.20040, 1.16770/

        DATA(qscat(1,5,ib,6),ib=1,17)  & ! Band 5, RH = 85%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00002  &
        ,0.00012, 0.00077, 0.00492, 0.03079, 0.17250, 0.65074  &
        ,1.55650, 2.06200, 1.35440, 1.19420, 1.16610/

        DATA(qscat(1,5,ib,7),ib=1,17)  & ! Band 5, RH = 86%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00002  &
        ,0.00012, 0.00078, 0.00498, 0.03115, 0.17353, 0.64773  &
        ,1.53320, 2.02230, 1.34560, 1.18820, 1.16460/

        DATA(qscat(1,5,ib,8),ib=1,17)  & ! Band 5, RH = 87%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00002  &
        ,0.00013, 0.00080, 0.00506, 0.03156, 0.17474, 0.64454  &
        ,1.50800, 1.97960, 1.33630, 1.18250, 1.16310/

        DATA(qscat(1,5,ib,9),ib=1,17)  & ! Band 5, RH = 88%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00002  &
        ,0.00013, 0.00081, 0.00515, 0.03206, 0.17618, 0.64116  &
        ,1.48080, 1.93360, 1.32620, 1.17700, 1.16170/

        DATA(qscat(1,5,ib,10),ib=1,17)  & ! Band 5, RH = 89%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00002  &
        ,0.00013, 0.00083, 0.00526, 0.03266, 0.17794, 0.63763  &
        ,1.45130, 1.88380, 1.31520, 1.17180, 1.16020/

        DATA(qscat(1,5,ib,11),ib=1,17)  & ! Band 5, RH = 90%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00002  &
        ,0.00013, 0.00085, 0.00539, 0.03340, 0.18011, 0.63397  &
        ,1.41930, 1.82990, 1.30320, 1.16690, 1.15870/

        DATA(qscat(1,5,ib,12),ib=1,17)  & ! Band 5, RH = 91%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00002  &
        ,0.00014, 0.00088, 0.00556, 0.03433, 0.18284, 0.63028  &
        ,1.38460, 1.77130, 1.29010, 1.16240, 1.15710/

        DATA(qscat(1,5,ib,13),ib=1,17)  & ! Band 5, RH = 92%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00002  &
        ,0.00014, 0.00091, 0.00577, 0.03553, 0.18635, 0.62671  &
        ,1.34680, 1.70760, 1.27540, 1.15830, 1.15540/

        DATA(qscat(1,5,ib,14),ib=1,17)  & ! Band 5, RH = 93%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00002  &
        ,0.00015, 0.00096, 0.00606, 0.03712, 0.19102, 0.62354  &
        ,1.30580, 1.63810, 1.25900, 1.15460, 1.15340/

        DATA(qscat(1,5,ib,15),ib=1,17)  & ! Band 5, RH = 94%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00002  &
        ,0.00016, 0.00102, 0.00647, 0.03932, 0.19743, 0.62135  &
        ,1.26170, 1.56250, 1.24040, 1.15110, 1.15120/

        DATA(qscat(1,5,ib,16),ib=1,17)  & ! Band 5, RH = 95%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00003  &
        ,0.00017, 0.00112, 0.00706, 0.04254, 0.20666, 0.62120  &
        ,1.21480, 1.48050, 1.21920, 1.14770, 1.14860/

        DATA(qscat(1,5,ib,17),ib=1,17)  & ! Band 5, RH = 96%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00003  &
        ,0.00020, 0.00127, 0.00802, 0.04762, 0.22074, 0.62533  &
        ,1.16630, 1.39270, 1.19490, 1.14410, 1.14570/

        DATA(qscat(1,5,ib,18),ib=1,17)  & ! Band 5, RH = 97%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00001, 0.00004  &
        ,0.00024, 0.00155, 0.00974, 0.05659, 0.24414, 0.63873  &
        ,1.11930, 1.30100, 1.16810, 1.13990, 1.14260/

        DATA(qscat(1,5,ib,19),ib=1,17)  & ! Band 5, RH = 98%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00001, 0.00005  &
        ,0.00033, 0.00216, 0.01358, 0.07558, 0.28828, 0.67420  &
        ,1.08080, 1.21010, 1.14150, 1.13550, 1.14010/

        DATA(qscat(1,5,ib,20),ib=1,17)  & ! Band 5, RH = 99%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00001, 0.00009  &
        ,0.00065, 0.00435, 0.02689, 0.13280, 0.39431, 0.77241  &
        ,1.06860, 1.13330, 1.12400, 1.13420, 1.14010/

        DATA(qscat(1,6,ib,1),ib=1,17)  & ! Band 6, RH = 80%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00001, 0.00008  &
        ,0.00052, 0.00327, 0.02003, 0.10642, 0.35826, 0.81682  &
        ,1.26810, 1.40030, 1.24000, 1.12530, 1.11970/

        DATA(qscat(1,6,ib,2),ib=1,17)  & ! Band 6, RH = 81%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00001, 0.00008  &
        ,0.00052, 0.00330, 0.02017, 0.10695, 0.35922, 0.81974  &
        ,1.27370, 1.40510, 1.24010, 1.12480, 1.11950/

        DATA(qscat(1,6,ib,3),ib=1,17)  & ! Band 6, RH = 82%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00001, 0.00008  &
        ,0.00053, 0.00332, 0.02032, 0.10754, 0.36029, 0.82298  &
        ,1.28000, 1.41040, 1.24010, 1.12440, 1.11930/

        DATA(qscat(1,6,ib,4),ib=1,17)  & ! Band 6, RH = 83%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00001, 0.00008  &
        ,0.00053, 0.00335, 0.02049, 0.10820, 0.36151, 0.82662  &
        ,1.28690, 1.41620, 1.24010, 1.12380, 1.11900/

        DATA(qscat(1,6,ib,5),ib=1,17)  & ! Band 6, RH = 84%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00001, 0.00008  &
        ,0.00054, 0.00339, 0.02068, 0.10895, 0.36290, 0.83072  &
        ,1.29470, 1.42260, 1.24000, 1.12330, 1.11870/

        DATA(qscat(1,6,ib,6),ib=1,17)  & ! Band 6, RH = 85%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00001, 0.00009  &
        ,0.00054, 0.00343, 0.02090, 0.10980, 0.36451, 0.83537  &
        ,1.30350, 1.42980, 1.23970, 1.12270, 1.11840/

        DATA(qscat(1,6,ib,7),ib=1,17)  & ! Band 6, RH = 86%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00001, 0.00009  &
        ,0.00055, 0.00347, 0.02115, 0.11078, 0.36638, 0.84071  &
        ,1.31340, 1.43780, 1.23920, 1.12210, 1.11810/

        DATA(qscat(1,6,ib,8),ib=1,17)  & ! Band 6, RH = 87%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00001, 0.00009  &
        ,0.00056, 0.00352, 0.02145, 0.11191, 0.36859, 0.84687  &
        ,1.32470, 1.44680, 1.23850, 1.12150, 1.11770/

        DATA(qscat(1,6,ib,9),ib=1,17)  & ! Band 6, RH = 88%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00001, 0.00009  &
        ,0.00057, 0.00358, 0.02179, 0.11324, 0.37122, 0.85408  &
        ,1.33790, 1.45690, 1.23740, 1.12090, 1.11720/

        DATA(qscat(1,6,ib,10),ib=1,17)  & ! Band 6, RH = 89%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00001, 0.00009  &
        ,0.00058, 0.00366, 0.02220, 0.11482, 0.37441, 0.86263  &
        ,1.35320, 1.46850, 1.23580, 1.12020, 1.11660/

        DATA(qscat(1,6,ib,11),ib=1,17)  & ! Band 6, RH = 90%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00001, 0.00009  &
        ,0.00059, 0.00374, 0.02270, 0.11672, 0.37835, 0.87292  &
        ,1.37130, 1.48170, 1.23350, 1.11970, 1.11600/

        DATA(qscat(1,6,ib,12),ib=1,17)  & ! Band 6, RH = 91%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00002, 0.00010  &
        ,0.00061, 0.00385, 0.02332, 0.11907, 0.38330, 0.88553  &
        ,1.39310, 1.49700, 1.23030, 1.11940, 1.11530/

        DATA(qscat(1,6,ib,13),ib=1,17)  & ! Band 6, RH = 92%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00002, 0.00010  &
        ,0.00063, 0.00399, 0.02410, 0.12201, 0.38971, 0.90135  &
        ,1.41980, 1.51470, 1.22580, 1.11930, 1.11440/

        DATA(qscat(1,6,ib,14),ib=1,17)  & ! Band 6, RH = 93%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00002, 0.00010  &
        ,0.00066, 0.00417, 0.02512, 0.12583, 0.39827, 0.92175  &
        ,1.45330, 1.53550, 1.21890, 1.11980, 1.11340/

        DATA(qscat(1,6,ib,15),ib=1,17)  & ! Band 6, RH = 94%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00002, 0.00011  &
        ,0.00070, 0.00441, 0.02650, 0.13094, 0.41018, 0.94902  &
        ,1.49640, 1.56020, 1.20840, 1.12110, 1.11240/

        DATA(qscat(1,6,ib,16),ib=1,17)  & ! Band 6, RH = 95%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00002, 0.00012  &
        ,0.00076, 0.00477, 0.02848, 0.13812, 0.42765, 0.98717  &
        ,1.55410, 1.58920, 1.19320, 1.12370, 1.11150/

        DATA(qscat(1,6,ib,17),ib=1,17)  & ! Band 6, RH = 96%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00002, 0.00013  &
        ,0.00084, 0.00531, 0.03154, 0.14888, 0.45526, 1.04370  &
        ,1.63410, 1.61940, 1.16990, 1.12790, 1.11110/

        DATA(qscat(1,6,ib,18),ib=1,17)  & ! Band 6, RH = 97%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00002, 0.00015  &
        ,0.00099, 0.00626, 0.03682, 0.16662, 0.50348, 1.13400  &
        ,1.74790, 1.64320, 1.13760, 1.13120, 1.11060/

        DATA(qscat(1,6,ib,19),ib=1,17)  & ! Band 6, RH = 98%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00003, 0.00020  &
        ,0.00131, 0.00830, 0.04787, 0.20080, 0.59996, 1.29650  &
        ,1.92570, 1.62550, 1.10700, 1.12130, 1.10660/

        DATA(qscat(1,6,ib,20),ib=1,17)  & ! Band 6, RH = 99%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00005, 0.00034  &
        ,0.00234, 0.01505, 0.08212, 0.29489, 0.83761, 1.68570  &
        ,2.17570, 1.39960, 1.17480, 1.10880, 1.10420/

        DATA(qscat(1,7,ib,1),ib=1,17)  & ! Band 7, RH = 80%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00003, 0.00018  &
        ,0.00113, 0.00725, 0.04681, 0.27904, 1.07380, 1.99160  &
        ,1.71550, 1.15450, 1.22020, 1.20580, 1.19920/

        DATA(qscat(1,7,ib,2),ib=1,17)  & ! Band 7, RH = 81%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00003, 0.00018  &
        ,0.00113, 0.00727, 0.04691, 0.27900, 1.06960, 1.98870  &
        ,1.72160, 1.15300, 1.22000, 1.20480, 1.19820/

        DATA(qscat(1,7,ib,3),ib=1,17)  & ! Band 7, RH = 82%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00003, 0.00018  &
        ,0.00113, 0.00729, 0.04701, 0.27895, 1.06510, 1.98550  &
        ,1.72820, 1.15150, 1.21970, 1.20380, 1.19710/

        DATA(qscat(1,7,ib,4),ib=1,17)  & ! Band 7, RH = 83%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00003, 0.00018  &
        ,0.00114, 0.00732, 0.04713, 0.27889, 1.06020, 1.98180  &
        ,1.73540, 1.14990, 1.21950, 1.20270, 1.19600/

        DATA(qscat(1,7,ib,5),ib=1,17)  & ! Band 7, RH = 84%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00003, 0.00018  &
        ,0.00114, 0.00735, 0.04727, 0.27883, 1.05480, 1.97770  &
        ,1.74340, 1.14840, 1.21920, 1.20140, 1.19470/

        DATA(qscat(1,7,ib,6),ib=1,17)  & ! Band 7, RH = 85%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00003, 0.00018  &
        ,0.00115, 0.00738, 0.04741, 0.27876, 1.04880, 1.97310  &
        ,1.75220, 1.14680, 1.21900, 1.20010, 1.19330/

        DATA(qscat(1,7,ib,7),ib=1,17)  & ! Band 7, RH = 86%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00003, 0.00018  &
        ,0.00115, 0.00741, 0.04758, 0.27867, 1.04230, 1.96780  &
        ,1.76190, 1.14520, 1.21870, 1.19860, 1.19170/

        DATA(qscat(1,7,ib,8),ib=1,17)  & ! Band 7, RH = 87%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00003, 0.00018  &
        ,0.00116, 0.00745, 0.04778, 0.27857, 1.03490, 1.96170  &
        ,1.77280, 1.14380, 1.21840, 1.19700, 1.19000/

        DATA(qscat(1,7,ib,9),ib=1,17)  & ! Band 7, RH = 88%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00003, 0.00018  &
        ,0.00117, 0.00749, 0.04800, 0.27846, 1.02680, 1.95460  &
        ,1.78500, 1.14250, 1.21810, 1.19530, 1.18800/

        DATA(qscat(1,7,ib,10),ib=1,17)  & ! Band 7, RH = 89%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00003, 0.00018  &
        ,0.00118, 0.00755, 0.04826, 0.27833, 1.01760, 1.94630  &
        ,1.79880, 1.14160, 1.21770, 1.19330, 1.18570/

        DATA(qscat(1,7,ib,11),ib=1,17)  & ! Band 7, RH = 90%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00003, 0.00019  &
        ,0.00119, 0.00761, 0.04858, 0.27817, 1.00720, 1.93650  &
        ,1.81440, 1.14110, 1.21730, 1.19100, 1.18320/

        DATA(qscat(1,7,ib,12),ib=1,17)  & ! Band 7, RH = 91%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00003, 0.00019  &
        ,0.00120, 0.00769, 0.04896, 0.27800, 0.99543, 1.92480  &
        ,1.83230, 1.14150, 1.21670, 1.18850, 1.18020/

        DATA(qscat(1,7,ib,13),ib=1,17)  & ! Band 7, RH = 92%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00003, 0.00019  &
        ,0.00121, 0.00778, 0.04943, 0.27781, 0.98196, 1.91050  &
        ,1.85290, 1.14310, 1.21590, 1.18550, 1.17670/

        DATA(qscat(1,7,ib,14),ib=1,17)  & ! Band 7, RH = 93%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00003, 0.00019  &
        ,0.00123, 0.00790, 0.05003, 0.27762, 0.96649, 1.89300  &
        ,1.87670, 1.14700, 1.21470, 1.18200, 1.17260/

        DATA(qscat(1,7,ib,15),ib=1,17)  & ! Band 7, RH = 94%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00003, 0.00020  &
        ,0.00126, 0.00807, 0.05084, 0.27747, 0.94870, 1.87110  &
        ,1.90410, 1.15420, 1.21280, 1.17770, 1.16770/

        DATA(qscat(1,7,ib,16),ib=1,17)  & ! Band 7, RH = 95%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00003, 0.00020  &
        ,0.00130, 0.00830, 0.05197, 0.27745, 0.92833, 1.84320  &
        ,1.93510, 1.16690, 1.20960, 1.17220, 1.16170/

        DATA(qscat(1,7,ib,17),ib=1,17)  & ! Band 7, RH = 96%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00003, 0.00021  &
        ,0.00135, 0.00865, 0.05369, 0.27785, 0.90559, 1.80740  &
        ,1.96830, 1.18820, 1.20390, 1.16480, 1.15420/

        DATA(qscat(1,7,ib,18),ib=1,17)  & ! Band 7, RH = 97%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00003, 0.00023  &
        ,0.00145, 0.00924, 0.05664, 0.27950, 0.88222, 1.76160  &
        ,1.99960, 1.22220, 1.19440, 1.15460, 1.14470/

        DATA(qscat(1,7,ib,19),ib=1,17)  & ! Band 7, RH = 98%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00004, 0.00025  &
        ,0.00165, 0.01049, 0.06282, 0.28553, 0.86576, 1.71140  &
        ,2.02930, 1.26800, 1.17990, 1.14030, 1.13180/

        DATA(qscat(1,7,ib,20),ib=1,17)  & ! Band 7, RH = 99%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00005, 0.00034  &
        ,0.00227, 0.01454, 0.08216, 0.31465, 0.89350, 1.72680  &
        ,2.05730, 1.28860, 1.16810, 1.12220, 1.11360/

        DATA(qscat(1,8,ib,1),ib=1,17)  & ! Band 8, RH = 80%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00004, 0.00025  &
        ,0.00155, 0.00939, 0.05049, 0.19427, 0.60629, 1.41840  &
        ,2.31550, 2.08250, 1.57030, 1.27520, 1.17020/

        DATA(qscat(1,8,ib,2),ib=1,17)  & ! Band 8, RH = 81%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00004, 0.00025  &
        ,0.00158, 0.00955, 0.05129, 0.19688, 0.61344, 1.43180  &
        ,2.32740, 2.07790, 1.57040, 1.27600, 1.16880/

        DATA(qscat(1,8,ib,3),ib=1,17)  & ! Band 8, RH = 82%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00004, 0.00026  &
        ,0.00160, 0.00973, 0.05218, 0.19978, 0.62137, 1.44650  &
        ,2.34020, 2.07270, 1.57040, 1.27710, 1.16730/

        DATA(qscat(1,8,ib,4),ib=1,17)  & ! Band 8, RH = 83%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00004, 0.00026  &
        ,0.00164, 0.00993, 0.05318, 0.20303, 0.63022, 1.46290  &
        ,2.35420, 2.06690, 1.57020, 1.27850, 1.16570/

        DATA(qscat(1,8,ib,5),ib=1,17)  & ! Band 8, RH = 84%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00004, 0.00027  &
        ,0.00168, 0.01016, 0.05430, 0.20669, 0.64014, 1.48110  &
        ,2.36960, 2.06030, 1.56990, 1.28030, 1.16400/

        DATA(qscat(1,8,ib,6),ib=1,17)  & ! Band 8, RH = 85%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00004, 0.00027  &
        ,0.00172, 0.01042, 0.05557, 0.21084, 0.65136, 1.50160  &
        ,2.38640, 2.05270, 1.56930, 1.28260, 1.16240/

        DATA(qscat(1,8,ib,7),ib=1,17)  & ! Band 8, RH = 86%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00004, 0.00028  &
        ,0.00177, 0.01072, 0.05703, 0.21558, 0.66415, 1.52470  &
        ,2.40500, 2.04380, 1.56820, 1.28530, 1.16090/

        DATA(qscat(1,8,ib,8),ib=1,17)  & ! Band 8, RH = 87%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00005, 0.00029  &
        ,0.00183, 0.01107, 0.05872, 0.22107, 0.67885, 1.55090  &
        ,2.42540, 2.03300, 1.56590, 1.28840, 1.15970/

        DATA(qscat(1,8,ib,9),ib=1,17)  & ! Band 8, RH = 88%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00005, 0.00030  &
        ,0.00190, 0.01148, 0.06070, 0.22748, 0.69593, 1.58110  &
        ,2.44790, 2.01960, 1.56230, 1.29190, 1.15900/

        DATA(qscat(1,8,ib,10),ib=1,17)  & ! Band 8, RH = 89%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00005, 0.00031  &
        ,0.00198, 0.01197, 0.06304, 0.23508, 0.71603, 1.61610  &
        ,2.47260, 2.00280, 1.55720, 1.29560, 1.15910/

        DATA(qscat(1,8,ib,11),ib=1,17)  & ! Band 8, RH = 90%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00005, 0.00033  &
        ,0.00208, 0.01256, 0.06586, 0.24421, 0.74000, 1.65710  &
        ,2.49940, 1.98150, 1.55080, 1.29900, 1.16030/

        DATA(qscat(1,8,ib,12),ib=1,17)  & ! Band 8, RH = 91%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00005, 0.00035  &
        ,0.00220, 0.01329, 0.06931, 0.25542, 0.76907, 1.70580  &
        ,2.52800, 1.95550, 1.54190, 1.30130, 1.16270/

        DATA(qscat(1,8,ib,13),ib=1,17)  & ! Band 8, RH = 92%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00006, 0.00037  &
        ,0.00235, 0.01422, 0.07364, 0.26948, 0.80499, 1.76420  &
        ,2.55800, 1.92480, 1.52840, 1.30110, 1.16600/

        DATA(qscat(1,8,ib,14),ib=1,17)  & ! Band 8, RH = 93%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00006, 0.00040  &
        ,0.00256, 0.01543, 0.07921, 0.28764, 0.85037, 1.83530  &
        ,2.58880, 1.88630, 1.51060, 1.29640, 1.16870/

        DATA(qscat(1,8,ib,15),ib=1,17)  & ! Band 8, RH = 94%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00007, 0.00045  &
        ,0.00283, 0.01708, 0.08665, 0.31195, 0.90923, 1.92300  &
        ,2.62070, 1.83390, 1.48550, 1.28440, 1.16770/

        DATA(qscat(1,8,ib,16),ib=1,17)  & ! Band 8, RH = 95%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00008, 0.00051  &
        ,0.00323, 0.01944, 0.09704, 0.34604, 0.98798, 2.03400  &
        ,2.64840, 1.77250, 1.45270, 1.26450, 1.16000/

        DATA(qscat(1,8,ib,17),ib=1,17)  & ! Band 8, RH = 96%
        /0.00000, 0.00000, 0.00000, 0.00002, 0.00009, 0.00060  &
        ,0.00384, 0.02308, 0.11251, 0.39687, 1.09790, 2.18180  &
        ,2.65190, 1.69420, 1.41300, 1.24430, 1.15140/

        DATA(qscat(1,8,ib,18),ib=1,17)  & ! Band 8, RH = 97%
        /0.00000, 0.00000, 0.00000, 0.00002, 0.00011, 0.00076  &
        ,0.00491, 0.02938, 0.13787, 0.47899, 1.26320, 2.38650  &
        ,2.60680, 1.61500, 1.38070, 1.24570, 1.15620/

        DATA(qscat(1,8,ib,19),ib=1,17)  & ! Band 8, RH = 98%
        /0.00000, 0.00000, 0.00000, 0.00003, 0.00016, 0.00109  &
        ,0.00717, 0.04257, 0.18672, 0.62789, 1.55310, 2.65130  &
        ,2.41220, 1.56670, 1.39370, 1.24530, 1.14310/

        DATA(qscat(1,8,ib,20),ib=1,17)  & ! Band 8, RH = 99%
        /0.00000, 0.00000, 0.00000, 0.00004, 0.00029, 0.00214  &
        ,0.01460, 0.08332, 0.32481, 0.99785, 2.15220, 2.88160  &
        ,1.86390, 1.53510, 1.32100, 1.20400, 1.13080/

! Qscat for Sea Salt:
        DATA(qscat(2,1,ib,1),ib=1,17)  & ! Band 1, RH = 80%
        /0.00003, 0.00023, 0.00149, 0.00939, 0.04968, 0.22147  &
        ,0.77965, 1.94100, 3.20120, 2.66530, 2.25360, 2.08700  &
        ,1.94680, 1.87010, 1.79000, 1.71870, 1.63600/

        DATA(qscat(2,1,ib,2),ib=1,17)  & ! Band 1, RH = 81%
        /0.00004, 0.00024, 0.00156, 0.00984, 0.05183, 0.22852  &
        ,0.79833, 1.97240, 3.21400, 2.63550, 2.25630, 2.08880  &
        ,1.95830, 1.86100, 1.79370, 1.71480, 1.63280/

        DATA(qscat(2,1,ib,3),ib=1,17)  & ! Band 1, RH = 82%
        /0.00004, 0.00025, 0.00164, 0.01034, 0.05417, 0.23611  &
        ,0.81821, 2.00500, 3.22470, 2.60610, 2.24820, 2.07650  &
        ,1.95460, 1.85800, 1.78920, 1.71450, 1.63000/

        DATA(qscat(2,1,ib,4),ib=1,17)  & ! Band 1, RH = 83%
        /0.00004, 0.00027, 0.00173, 0.01088, 0.05672, 0.24437  &
        ,0.83960, 2.03930, 3.23330, 2.57200, 2.25660, 2.07190  &
        ,1.95030, 1.85180, 1.78940, 1.71120, 1.62450/

        DATA(qscat(2,1,ib,5),ib=1,17)  & ! Band 1, RH = 84%
        /0.00004, 0.00028, 0.00182, 0.01148, 0.05954, 0.25345  &
        ,0.86288, 2.07580, 3.24040, 2.53150, 2.25600, 2.06590  &
        ,1.96800, 1.85280, 1.78620, 1.70630, 1.62150/

        DATA(qscat(2,1,ib,6),ib=1,17)  & ! Band 1, RH = 85%
        /0.00004, 0.00030, 0.00193, 0.01216, 0.06268, 0.26355  &
        ,0.88858, 2.11520, 3.24720, 2.48740, 2.24190, 2.05630  &
        ,1.96270, 1.85210, 1.77880, 1.70090, 1.61750/

        DATA(qscat(2,1,ib,7),ib=1,17)  & ! Band 1, RH = 86%
        /0.00005, 0.00032, 0.00206, 0.01293, 0.06622, 0.27494  &
        ,0.91733, 2.15860, 3.25540, 2.44900, 2.23810, 2.05370  &
        ,1.95530, 1.85060, 1.77410, 1.69820, 1.61130/

        DATA(qscat(2,1,ib,8),ib=1,17)  & ! Band 1, RH = 87%
        /0.00005, 0.00034, 0.00220, 0.01382, 0.07027, 0.28797  &
        ,0.94994, 2.20710, 3.26490, 2.41120, 2.23060, 2.03180  &
        ,1.94810, 1.84910, 1.77440, 1.69360, 1.60810/

        DATA(qscat(2,1,ib,9),ib=1,17)  & ! Band 1, RH = 88%
        /0.00005, 0.00036, 0.00237, 0.01487, 0.07495, 0.30310  &
        ,0.98740, 2.26210, 3.27140, 2.36400, 2.21510, 2.01950  &
        ,1.93600, 1.84180, 1.77430, 1.68870, 1.60190/

        DATA(qscat(2,1,ib,10),ib=1,17)  & ! Band 1, RH = 89%
        /0.00006, 0.00039, 0.00256, 0.01611, 0.08042, 0.32095  &
        ,1.03090, 2.32420, 3.26950, 2.31350, 2.19120, 2.01160  &
        ,1.92730, 1.83760, 1.76450, 1.68410, 1.59630/

        DATA(qscat(2,1,ib,11),ib=1,17)  & ! Band 1, RH = 90%
        /0.00006, 0.00043, 0.00281, 0.01762, 0.08693, 0.34238  &
        ,1.08180, 2.39320, 3.25850, 2.27510, 2.15490, 1.98210  &
        ,1.89400, 1.82910, 1.75840, 1.67620, 1.59180/

        DATA(qscat(2,1,ib,12),ib=1,17)  & ! Band 1, RH = 91%
        /0.00007, 0.00047, 0.00310, 0.01948, 0.09478, 0.36857  &
        ,1.14150, 2.46850, 3.24460, 2.22580, 2.11910, 1.96730  &
        ,1.89710, 1.83400, 1.75170, 1.67160, 1.58420/

        DATA(qscat(2,1,ib,13),ib=1,17)  & ! Band 1, RH = 92%
        /0.00008, 0.00052, 0.00348, 0.02183, 0.10442, 0.40118  &
        ,1.21190, 2.55100, 3.22300, 2.19160, 2.07350, 1.95850  &
        ,1.88940, 1.82850, 1.74540, 1.66310, 1.57630/

        DATA(qscat(2,1,ib,14),ib=1,17)  & ! Band 1, RH = 93%
        /0.00009, 0.00060, 0.00398, 0.02488, 0.11652, 0.44251  &
        ,1.29540, 2.64610, 3.17240, 2.15710, 2.03900, 1.96080  &
        ,1.89320, 1.81750, 1.74010, 1.65410, 1.56660/

        DATA(qscat(2,1,ib,15),ib=1,17)  & ! Band 1, RH = 94%
        /0.00010, 0.00069, 0.00464, 0.02899, 0.13215, 0.49579  &
        ,1.39680, 2.76130, 3.10740, 2.14490, 1.99140, 1.97110  &
        ,1.88210, 1.80270, 1.72890, 1.64440, 1.55620/

        DATA(qscat(2,1,ib,16),ib=1,17)  & ! Band 1, RH = 95%
        /0.00011, 0.00082, 0.00559, 0.03478, 0.15313, 0.56556  &
        ,1.52720, 2.88630, 2.99490, 2.15750, 1.98940, 1.96780  &
        ,1.87020, 1.79730, 1.71910, 1.63360, 1.54360/

        DATA(qscat(2,1,ib,17),ib=1,17)  & ! Band 1, RH = 96%
        /0.00014, 0.00102, 0.00704, 0.04349, 0.18298, 0.65921  &
        ,1.70550, 3.02210, 2.83180, 2.19030, 2.03520, 1.93300  &
        ,1.86740, 1.78680, 1.70660, 1.61860, 1.52670/

        DATA(qscat(2,1,ib,18),ib=1,17)  & ! Band 1, RH = 97%
        /0.00018, 0.00135, 0.00949, 0.05791, 0.22988, 0.79447  &
        ,1.94980, 3.15690, 2.58510, 2.22550, 2.06930, 1.93750  &
        ,1.84350, 1.76800, 1.68880, 1.60050, 1.51170/

        DATA(qscat(2,1,ib,19),ib=1,17)  & ! Band 1, RH = 98%
        /0.00025, 0.00200, 0.01443, 0.08564, 0.31960, 1.03340  &
        ,2.32480, 3.22770, 2.27650, 2.14850, 1.98140, 1.89820  &
        ,1.81540, 1.74170, 1.66240, 1.57390, 1.48310/

        DATA(qscat(2,1,ib,20),ib=1,17)  & ! Band 1, RH = 99%
        /0.00041, 0.00372, 0.02853, 0.15654, 0.56295, 1.55470  &
        ,2.91620, 2.92000, 2.16930, 2.00640, 1.94160, 1.85940  &
        ,1.78410, 1.70410, 1.61750, 1.52570, 1.43480/

        DATA(qscat(2,2,ib,1),ib=1,17)  & ! Band 2, RH = 80%
        /0.00072, 0.00477, 0.02963, 0.15046, 0.53715, 1.54660  &
        ,3.05940, 3.42630, 2.44470, 2.26150, 2.22730, 2.15380  &
        ,2.12090, 2.08190, 2.05890, 2.03700, 2.01700/

        DATA(qscat(2,2,ib,2),ib=1,17)  & ! Band 2, RH = 81%
        /0.00076, 0.00501, 0.03099, 0.15589, 0.55370, 1.57750  &
        ,3.09000, 3.40460, 2.43820, 2.27160, 2.22420, 2.15200  &
        ,2.11780, 2.08800, 2.05630, 2.03580, 2.01610/

        DATA(qscat(2,2,ib,3),ib=1,17)  & ! Band 2, RH = 82%
        /0.00080, 0.00526, 0.03246, 0.16170, 0.57131, 1.61050  &
        ,3.12060, 3.37870, 2.43390, 2.27720, 2.21310, 2.15960  &
        ,2.12290, 2.09040, 2.05620, 2.03510, 2.01490/

        DATA(qscat(2,2,ib,4),ib=1,17)  & ! Band 2, RH = 83%
        /0.00084, 0.00553, 0.03407, 0.16799, 0.59015, 1.64610  &
        ,3.15130, 3.34920, 2.43020, 2.28990, 2.19800, 2.15550  &
        ,2.11540, 2.08230, 2.05730, 2.03420, 2.01430/

        DATA(qscat(2,2,ib,5),ib=1,17)  & ! Band 2, RH = 84%
        /0.00089, 0.00583, 0.03584, 0.17484, 0.61047, 1.68470  &
        ,3.18240, 3.31810, 2.42230, 2.30020, 2.19300, 2.15860  &
        ,2.11740, 2.08240, 2.05500, 2.03310, 2.01400/

        DATA(qscat(2,2,ib,6),ib=1,17)  & ! Band 2, RH = 85%
        /0.00094, 0.00617, 0.03782, 0.18240, 0.63257, 1.72710  &
        ,3.21490, 3.28710, 2.41990, 2.30660, 2.18880, 2.15850  &
        ,2.11360, 2.08300, 2.05380, 2.03290, 2.01300/

        DATA(qscat(2,2,ib,7),ib=1,17)  & ! Band 2, RH = 86%
        /0.00099, 0.00655, 0.04006, 0.19083, 0.65683, 1.77380  &
        ,3.24990, 3.25220, 2.41540, 2.31860, 2.18580, 2.16170  &
        ,2.11070, 2.07830, 2.05340, 2.03210, 2.01230/

        DATA(qscat(2,2,ib,8),ib=1,17)  & ! Band 2, RH = 87%
        /0.00106, 0.00698, 0.04262, 0.20034, 0.68374, 1.82540  &
        ,3.28840, 3.20860, 2.40670, 2.32910, 2.18470, 2.15730  &
        ,2.11170, 2.07770, 2.05120, 2.03050, 2.01070/

        DATA(qscat(2,2,ib,9),ib=1,17)  & ! Band 2, RH = 88%
        /0.00113, 0.00749, 0.04559, 0.21122, 0.71395, 1.88270  &
        ,3.33000, 3.15580, 2.40060, 2.33690, 2.18390, 2.15110  &
        ,2.10790, 2.07670, 2.05180, 2.03040, 2.00970/

        DATA(qscat(2,2,ib,10),ib=1,17)  & ! Band 2, RH = 89%
        /0.00122, 0.00809, 0.04908, 0.22387, 0.74833, 1.94650  &
        ,3.37200, 3.10070, 2.39400, 2.34940, 2.18590, 2.15040  &
        ,2.10630, 2.07600, 2.05060, 2.02860, 2.00890/

        DATA(qscat(2,2,ib,11),ib=1,17)  & ! Band 2, RH = 90%
        /0.00132, 0.00881, 0.05324, 0.23881, 0.78812, 2.01800  &
        ,3.41120, 3.03900, 2.38420, 2.34650, 2.19670, 2.15080  &
        ,2.10790, 2.07760, 2.05020, 2.02720, 2.00780/

        DATA(qscat(2,2,ib,12),ib=1,17)  & ! Band 2, RH = 91%
        /0.00144, 0.00969, 0.05830, 0.25683, 0.83512, 2.09920  &
        ,3.44820, 2.96140, 2.37110, 2.34970, 2.20440, 2.15580  &
        ,2.10690, 2.07290, 2.04620, 2.02620, 2.00630/

        DATA(qscat(2,2,ib,13),ib=1,17)  & ! Band 2, RH = 92%
        /0.00160, 0.01079, 0.06456, 0.27906, 0.89204, 2.19410  &
        ,3.48810, 2.88340, 2.35240, 2.33490, 2.20920, 2.15610  &
        ,2.10200, 2.06850, 2.04480, 2.02490, 2.00460/

        DATA(qscat(2,2,ib,14),ib=1,17)  & ! Band 2, RH = 93%
        /0.00179, 0.01221, 0.07250, 0.30731, 0.96288, 2.30790  &
        ,3.52540, 2.78770, 2.32800, 2.30970, 2.20170, 2.13850  &
        ,2.09820, 2.06730, 2.04360, 2.02430, 2.00360/

        DATA(qscat(2,2,ib,15),ib=1,17)  & ! Band 2, RH = 94%
        /0.00205, 0.01411, 0.08285, 0.34453, 1.05330, 2.44430  &
        ,3.54380, 2.68930, 2.29120, 2.28690, 2.19390, 2.13260  &
        ,2.09760, 2.06590, 2.04180, 2.02050, 1.99920/

        DATA(qscat(2,2,ib,16),ib=1,17)  & ! Band 2, RH = 95%
        /0.00241, 0.01674, 0.09687, 0.39566, 1.17080, 2.60430  &
        ,3.54910, 2.58460, 2.25250, 2.27190, 2.19490, 2.12740  &
        ,2.09050, 2.06060, 2.03810, 2.01890, 1.99730/

        DATA(qscat(2,2,ib,17),ib=1,17)  & ! Band 2, RH = 96%
        /0.00292, 0.02064, 0.11690, 0.46925, 1.32690, 2.80480  &
        ,3.49980, 2.49040, 2.22630, 2.25870, 2.17460, 2.12250  &
        ,2.08490, 2.05700, 2.03570, 2.01540, 1.99400/

        DATA(qscat(2,2,ib,18),ib=1,17)  & ! Band 2, RH = 97%
        /0.00373, 0.02695, 0.14781, 0.58045, 1.55090, 3.05660  &
        ,3.36790, 2.42100, 2.26320, 2.20750, 2.14610, 2.11300  &
        ,2.08050, 2.05380, 2.03080, 2.01050, 1.98990/

        DATA(qscat(2,2,ib,19),ib=1,17)  & ! Band 2, RH = 98%
        /0.00518, 0.03878, 0.20249, 0.76440, 1.91440, 3.35280  &
        ,3.06560, 2.38420, 2.34390, 2.18810, 2.14700, 2.10430  &
        ,2.07410, 2.04700, 2.02540, 2.00370, 1.98210/

        DATA(qscat(2,2,ib,20),ib=1,17)  & ! Band 2, RH = 99%
        /0.00852, 0.06786, 0.33636, 1.18000, 2.58570, 3.53280  &
        ,2.53990, 2.23140, 2.25990, 2.18940, 2.12290, 2.08720  &
        ,2.05920, 2.03590, 2.01470, 1.99480, 1.96930/

        DATA(qscat(2,3,ib,1),ib=1,17)  & ! Band 3, RH = 80%
        /0.01390, 0.07864, 0.33227, 1.04460, 2.28450, 3.49860  &
        ,2.98590, 2.19520, 2.34920, 2.22500, 2.14170, 2.11090  &
        ,2.08200, 2.06070, 2.04250, 2.03250, 2.02210/

        DATA(qscat(2,3,ib,2),ib=1,17)  & ! Band 3, RH = 81%
        /0.01464, 0.08201, 0.34275, 1.06840, 2.31620, 3.50960  &
        ,2.95390, 2.20050, 2.34220, 2.22620, 2.13850, 2.10850  &
        ,2.08410, 2.05840, 2.04400, 2.03110, 2.02240/

        DATA(qscat(2,3,ib,3),ib=1,17)  & ! Band 3, RH = 82%
        /0.01537, 0.08541, 0.35400, 1.09360, 2.34930, 3.52010  &
        ,2.92280, 2.20660, 2.32180, 2.21840, 2.14490, 2.10860  &
        ,2.07680, 2.05670, 2.04120, 2.03020, 2.02260/

        DATA(qscat(2,3,ib,4),ib=1,17)  & ! Band 3, RH = 83%
        /0.01612, 0.08908, 0.36614, 1.12070, 2.38480, 3.53090  &
        ,2.89170, 2.21720, 2.31400, 2.20760, 2.14590, 2.10120  &
        ,2.07690, 2.05470, 2.04240, 2.03030, 2.02300/

        DATA(qscat(2,3,ib,5),ib=1,17)  & ! Band 3, RH = 84%
        /0.01694, 0.09310, 0.37934, 1.14990, 2.42310, 3.53880  &
        ,2.85880, 2.22360, 2.29950, 2.20100, 2.14770, 2.09830  &
        ,2.07260, 2.05660, 2.04210, 2.02990, 2.02280/

        DATA(qscat(2,3,ib,6),ib=1,17)  & ! Band 3, RH = 85%
        /0.01784, 0.09754, 0.39382, 1.18200, 2.46470, 3.54480  &
        ,2.82230, 2.23140, 2.28320, 2.19710, 2.14640, 2.10400  &
        ,2.07280, 2.05520, 2.04100, 2.02960, 2.02140/

        DATA(qscat(2,3,ib,7),ib=1,17)  & ! Band 3, RH = 86%
        /0.01885, 0.10249, 0.40984, 1.21740, 2.50980, 3.54840  &
        ,2.78170, 2.24250, 2.26770, 2.18710, 2.15200, 2.09650  &
        ,2.07650, 2.05450, 2.03970, 2.02960, 2.02240/

        DATA(qscat(2,3,ib,8),ib=1,17)  & ! Band 3, RH = 87%
        /0.01999, 0.10807, 0.42771, 1.25700, 2.55930, 3.54920  &
        ,2.73530, 2.25490, 2.25220, 2.18260, 2.14880, 2.10060  &
        ,2.07430, 2.05380, 2.03970, 2.02990, 2.02110/

        DATA(qscat(2,3,ib,9),ib=1,17)  & ! Band 3, RH = 88%
        /0.02130, 0.11445, 0.44784, 1.30140, 2.61360, 3.55050  &
        ,2.68520, 2.27430, 2.23320, 2.16970, 2.14060, 2.10390  &
        ,2.07340, 2.05340, 2.03930, 2.02810, 2.02010/

        DATA(qscat(2,3,ib,10),ib=1,17)  & ! Band 3, RH = 89%
        /0.02282, 0.12184, 0.47072, 1.35140, 2.67260, 3.54930  &
        ,2.63090, 2.28640, 2.21800, 2.16130, 2.12830, 2.09900  &
        ,2.07080, 2.05150, 2.03840, 2.02900, 2.02020/

        DATA(qscat(2,3,ib,11),ib=1,17)  & ! Band 3, RH = 90%
        /0.02461, 0.13049, 0.49706, 1.40850, 2.73990, 3.54170  &
        ,2.57320, 2.30470, 2.20080, 2.15640, 2.13060, 2.09570  &
        ,2.06830, 2.05220, 2.03740, 2.02830, 2.02120/

        DATA(qscat(2,3,ib,12),ib=1,17)  & ! Band 3, RH = 91%
        /0.02676, 0.14077, 0.52789, 1.47490, 2.81590, 3.52420  &
        ,2.51200, 2.32930, 2.18840, 2.16050, 2.12120, 2.09040  &
        ,2.06970, 2.05100, 2.03810, 2.02760, 2.01930/

        DATA(qscat(2,3,ib,13),ib=1,17)  & ! Band 3, RH = 92%
        /0.02940, 0.15318, 0.56477, 1.55310, 2.90210, 3.49850  &
        ,2.44450, 2.34220, 2.18330, 2.16790, 2.11720, 2.08720  &
        ,2.06750, 2.05040, 2.03660, 2.02660, 2.02000/

        DATA(qscat(2,3,ib,14),ib=1,17)  & ! Band 3, RH = 93%
        /0.03270, 0.16840, 0.61014, 1.64650, 2.99490, 3.45990  &
        ,2.36890, 2.36290, 2.18670, 2.17930, 2.12060, 2.09150  &
        ,2.06570, 2.04780, 2.03540, 2.02620, 2.01870/

        DATA(qscat(2,3,ib,15),ib=1,17)  & ! Band 3, RH = 94%
        /0.03695, 0.18755, 0.66800, 1.75840, 3.09880, 3.39210  &
        ,2.30260, 2.37550, 2.20480, 2.17790, 2.12780, 2.08800  &
        ,2.06190, 2.04560, 2.03420, 2.02410, 2.01840/

        DATA(qscat(2,3,ib,16),ib=1,17)  & ! Band 3, RH = 95%
        /0.04262, 0.21252, 0.74478, 1.89610, 3.21790, 3.30210  &
        ,2.23270, 2.37720, 2.22460, 2.16150, 2.11120, 2.08170  &
        ,2.06330, 2.04430, 2.03360, 2.02530, 2.01840/

        DATA(qscat(2,3,ib,17),ib=1,17)  & ! Band 3, RH = 96%
        /0.05057, 0.24710, 0.84932, 2.07510, 3.34410, 3.16020  &
        ,2.18690, 2.36780, 2.23490, 2.14050, 2.10870, 2.07700  &
        ,2.05890, 2.04380, 2.03100, 2.02410, 2.01650/

        DATA(qscat(2,3,ib,18),ib=1,17)  & ! Band 3, RH = 97%
        /0.06256, 0.29951, 0.99765, 2.31660, 3.47200, 2.94570  &
        ,2.19220, 2.32250, 2.21080, 2.13720, 2.10380, 2.07590  &
        ,2.05420, 2.04030, 2.03000, 2.02130, 2.01610/

        DATA(qscat(2,3,ib,19),ib=1,17)  & ! Band 3, RH = 98%
        /0.08290, 0.38942, 1.23930, 2.67350, 3.53380, 2.61790  &
        ,2.28460, 2.20780, 2.15350, 2.12650, 2.09820, 2.06710  &
        ,2.04930, 2.03730, 2.02680, 2.01960, 2.01470/

        DATA(qscat(2,3,ib,20),ib=1,17)  & ! Band 3, RH = 99%
        /0.12530, 0.57315, 1.71420, 3.21310, 3.30290, 2.21790  &
        ,2.38380, 2.23300, 2.15730, 2.10010, 2.08180, 2.05830  &
        ,2.04370, 2.03130, 2.02260, 2.01610, 2.01310/

        DATA(qscat(2,4,ib,1),ib=1,17)  & ! Band 4, RH = 80%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00001, 0.00004  &
        ,0.00024, 0.00152, 0.00959, 0.05679, 0.25333, 0.66013  &
        ,1.11390, 1.27140, 1.18610, 1.16220, 1.16420/

        DATA(qscat(2,4,ib,2),ib=1,17)  & ! Band 4, RH = 81%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00001, 0.00004  &
        ,0.00025, 0.00162, 0.01019, 0.06010, 0.26415, 0.67730  &
        ,1.12620, 1.27060, 1.18430, 1.16280, 1.16480/

        DATA(qscat(2,4,ib,3),ib=1,17)  & ! Band 4, RH = 82%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00001, 0.00004  &
        ,0.00027, 0.00173, 0.01085, 0.06372, 0.27565, 0.69517  &
        ,1.13860, 1.26930, 1.18250, 1.16350, 1.16540/

        DATA(qscat(2,4,ib,4),ib=1,17)  & ! Band 4, RH = 83%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00001, 0.00005  &
        ,0.00029, 0.00184, 0.01159, 0.06771, 0.28796, 0.71387  &
        ,1.15090, 1.26770, 1.18080, 1.16410, 1.16610/

        DATA(qscat(2,4,ib,5),ib=1,17)  & ! Band 4, RH = 84%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00001, 0.00005  &
        ,0.00031, 0.00198, 0.01241, 0.07215, 0.30124, 0.73358  &
        ,1.16340, 1.26570, 1.17910, 1.16480, 1.16670/

        DATA(qscat(2,4,ib,6),ib=1,17)  & ! Band 4, RH = 85%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00001, 0.00005  &
        ,0.00033, 0.00213, 0.01335, 0.07714, 0.31569, 0.75450  &
        ,1.17590, 1.26320, 1.17760, 1.16550, 1.16730/

        DATA(qscat(2,4,ib,7),ib=1,17)  & ! Band 4, RH = 86%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00001, 0.00006  &
        ,0.00036, 0.00230, 0.01443, 0.08281, 0.33155, 0.77685  &
        ,1.18850, 1.26020, 1.17620, 1.16630, 1.16790/

        DATA(qscat(2,4,ib,8),ib=1,17)  & ! Band 4, RH = 87%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00001, 0.00006  &
        ,0.00039, 0.00251, 0.01569, 0.08933, 0.34912, 0.80090  &
        ,1.20110, 1.25670, 1.17480, 1.16710, 1.16860/

        DATA(qscat(2,4,ib,9),ib=1,17)  & ! Band 4, RH = 88%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00001, 0.00007  &
        ,0.00043, 0.00275, 0.01718, 0.09692, 0.36876, 0.82693  &
        ,1.21370, 1.25250, 1.17360, 1.16790, 1.16920/

        DATA(qscat(2,4,ib,10),ib=1,17)  & ! Band 4, RH = 89%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00001, 0.00007  &
        ,0.00048, 0.00304, 0.01898, 0.10588, 0.39091, 0.85525  &
        ,1.22620, 1.24760, 1.17260, 1.16870, 1.16980/

        DATA(qscat(2,4,ib,11),ib=1,17)  & ! Band 4, RH = 90%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00001, 0.00008  &
        ,0.00053, 0.00340, 0.02118, 0.11660, 0.41616, 0.88620  &
        ,1.23830, 1.24190, 1.17170, 1.16960, 1.17050/

        DATA(qscat(2,4,ib,12),ib=1,17)  & ! Band 4, RH = 91%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00001, 0.00009  &
        ,0.00060, 0.00385, 0.02393, 0.12962, 0.44525, 0.92013  &
        ,1.24990, 1.23530, 1.17110, 1.17060, 1.17110/

        DATA(qscat(2,4,ib,13),ib=1,17)  & ! Band 4, RH = 92%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00002, 0.00011  &
        ,0.00070, 0.00444, 0.02746, 0.14571, 0.47917, 0.95741  &
        ,1.26060, 1.22770, 1.17080, 1.17150, 1.17170/

        DATA(qscat(2,4,ib,14),ib=1,17)  & ! Band 4, RH = 93%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00002, 0.00013  &
        ,0.00082, 0.00521, 0.03212, 0.16600, 0.51934, 0.99844  &
        ,1.26970, 1.21930, 1.17080, 1.17250, 1.17230/

        DATA(qscat(2,4,ib,15),ib=1,17)  & ! Band 4, RH = 94%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00002, 0.00015  &
        ,0.00099, 0.00630, 0.03852, 0.19226, 0.56782, 1.04370  &
        ,1.27620, 1.21000, 1.17120, 1.17350, 1.17290/

        DATA(qscat(2,4,ib,16),ib=1,17)  & ! Band 4, RH = 95%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00003, 0.00019  &
        ,0.00124, 0.00788, 0.04775, 0.22733, 0.62776, 1.09350  &
        ,1.27870, 1.20000, 1.17200, 1.17450, 1.17330/

        DATA(qscat(2,4,ib,17),ib=1,17)  & ! Band 4, RH = 96%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00004, 0.00025  &
        ,0.00163, 0.01041, 0.06207, 0.27624, 0.70427, 1.14850  &
        ,1.27500, 1.19010, 1.17320, 1.17540, 1.17370/

        DATA(qscat(2,4,ib,18),ib=1,17)  & ! Band 4, RH = 97%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00005, 0.00036  &
        ,0.00235, 0.01496, 0.08672, 0.34865, 0.80602, 1.20730  &
        ,1.26180, 1.18110, 1.17490, 1.17630, 1.17370/

        DATA(qscat(2,4,ib,19),ib=1,17)  & ! Band 4, RH = 98%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00009, 0.00061  &
        ,0.00396, 0.02502, 0.13662, 0.46683, 0.94793, 1.26120  &
        ,1.23430, 1.17550, 1.17670, 1.17680, 1.17330/

        DATA(qscat(2,4,ib,20),ib=1,17)  & ! Band 4, RH = 99%
        /0.00000, 0.00000, 0.00000, 0.00003, 0.00021, 0.00146  &
        ,0.00967, 0.05933, 0.27152, 0.70240, 1.14990, 1.27720  &
        ,1.19220, 1.17600, 1.17830, 1.17630, 1.17150/

        DATA(qscat(2,5,ib,1),ib=1,17)  & ! Band 5, RH = 80%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00003, 0.00016  &
        ,0.00102, 0.00634, 0.03632, 0.15966, 0.42250, 0.78438  &
        ,1.05280, 1.11540, 1.11680, 1.12900, 1.13530/

        DATA(qscat(2,5,ib,2),ib=1,17)  & ! Band 5, RH = 81%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00003, 0.00017  &
        ,0.00109, 0.00671, 0.03828, 0.16598, 0.43241, 0.79364  &
        ,1.05510, 1.11410, 1.11710, 1.12950, 1.13570/

        DATA(qscat(2,5,ib,3),ib=1,17)  & ! Band 5, RH = 82%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00003, 0.00018  &
        ,0.00115, 0.00712, 0.04042, 0.17270, 0.44281, 0.80317  &
        ,1.05740, 1.11300, 1.11750, 1.13000, 1.13600/

        DATA(qscat(2,5,ib,4),ib=1,17)  & ! Band 5, RH = 83%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00003, 0.00019  &
        ,0.00123, 0.00758, 0.04278, 0.17991, 0.45380, 0.81303  &
        ,1.05970, 1.11190, 1.11790, 1.13050, 1.13640/

        DATA(qscat(2,5,ib,5),ib=1,17)  & ! Band 5, RH = 84%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00003, 0.00021  &
        ,0.00131, 0.00809, 0.04540, 0.18769, 0.46551, 0.82330  &
        ,1.06210, 1.11090, 1.11850, 1.13110, 1.13680/

        DATA(qscat(2,5,ib,6),ib=1,17)  & ! Band 5, RH = 85%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00003, 0.00022  &
        ,0.00141, 0.00867, 0.04833, 0.19615, 0.47810, 0.83407  &
        ,1.06460, 1.11010, 1.11900, 1.13170, 1.13730/

        DATA(qscat(2,5,ib,7),ib=1,17)  & ! Band 5, RH = 86%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00004, 0.00024  &
        ,0.00152, 0.00934, 0.05166, 0.20546, 0.49174, 0.84542  &
        ,1.06710, 1.10940, 1.11970, 1.13230, 1.13770/

        DATA(qscat(2,5,ib,8),ib=1,17)  & ! Band 5, RH = 87%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00004, 0.00026  &
        ,0.00165, 0.01011, 0.05548, 0.21576, 0.50663, 0.85747  &
        ,1.06960, 1.10880, 1.12050, 1.13300, 1.13810/

        DATA(qscat(2,5,ib,9),ib=1,17)  & ! Band 5, RH = 88%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00004, 0.00029  &
        ,0.00181, 0.01103, 0.05992, 0.22728, 0.52304, 0.87032  &
        ,1.07230, 1.10840, 1.12130, 1.13370, 1.13860/

        DATA(qscat(2,5,ib,10),ib=1,17)  & ! Band 5, RH = 89%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00005, 0.00032  &
        ,0.00199, 0.01213, 0.06515, 0.24028, 0.54125, 0.88408  &
        ,1.07500, 1.10820, 1.12230, 1.13440, 1.13910/

        DATA(qscat(2,5,ib,11),ib=1,17)  & ! Band 5, RH = 90%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00005, 0.00035  &
        ,0.00222, 0.01347, 0.07140, 0.25507, 0.56162, 0.89887  &
        ,1.07780, 1.10820, 1.12340, 1.13530, 1.13960/

        DATA(qscat(2,5,ib,12),ib=1,17)  & ! Band 5, RH = 91%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00006, 0.00040  &
        ,0.00250, 0.01515, 0.07898, 0.27209, 0.58459, 0.91482  &
        ,1.08070, 1.10850, 1.12470, 1.13610, 1.14010/

        DATA(qscat(2,5,ib,13),ib=1,17)  & ! Band 5, RH = 92%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00007, 0.00046  &
        ,0.00287, 0.01728, 0.08834, 0.29191, 0.61069, 0.93207  &
        ,1.08360, 1.10920, 1.12610, 1.13710, 1.14060/

        DATA(qscat(2,5,ib,14),ib=1,17)  & ! Band 5, RH = 93%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00008, 0.00053  &
        ,0.00336, 0.02008, 0.10014, 0.31535, 0.64066, 0.95078  &
        ,1.08660, 1.11030, 1.12760, 1.13810, 1.14110/

        DATA(qscat(2,5,ib,15),ib=1,17)  & ! Band 5, RH = 94%
        /0.00000, 0.00000, 0.00000, 0.00002, 0.00010, 0.00064  &
        ,0.00403, 0.02390, 0.11542, 0.34367, 0.67546, 0.97112  &
        ,1.08960, 1.11190, 1.12940, 1.13910, 1.14160/

        DATA(qscat(2,5,ib,16),ib=1,17)  & ! Band 5, RH = 95%
        /0.00000, 0.00000, 0.00000, 0.00002, 0.00012, 0.00080  &
        ,0.00502, 0.02938, 0.13589, 0.37894, 0.71650, 0.99317  &
        ,1.09270, 1.11410, 1.13150, 1.14030, 1.14210/

        DATA(qscat(2,5,ib,17),ib=1,17)  & ! Band 5, RH = 96%
        /0.00000, 0.00000, 0.00000, 0.00003, 0.00016, 0.00105  &
        ,0.00657, 0.03779, 0.16456, 0.42484, 0.76588, 1.01690  &
        ,1.09610, 1.11720, 1.13380, 1.14150, 1.14250/

        DATA(qscat(2,5,ib,18),ib=1,17)  & ! Band 5, RH = 97%
        /0.00000, 0.00000, 0.00001, 0.00004, 0.00023, 0.00150  &
        ,0.00934, 0.05212, 0.20727, 0.48860, 0.82691, 1.04210  &
        ,1.10020, 1.12140, 1.13650, 1.14270, 1.14280/

        DATA(qscat(2,5,ib,19),ib=1,17)  & ! Band 5, RH = 98%
        /0.00000, 0.00000, 0.00001, 0.00006, 0.00038, 0.00248  &
        ,0.01535, 0.08081, 0.27711, 0.58578, 0.90485, 1.06790  &
        ,1.10620, 1.12710, 1.13960, 1.14380, 1.14280/

        DATA(qscat(2,5,ib,20),ib=1,17)  & ! Band 5, RH = 99%
        /0.00000, 0.00000, 0.00002, 0.00014, 0.00088, 0.00586  &
        ,0.03514, 0.15840, 0.41680, 0.75471, 1.00690, 1.09230  &
        ,1.11740, 1.13510, 1.14320, 1.14450, 1.14170/

        DATA(qscat(2,6,ib,1),ib=1,17)  & ! Band 6, RH = 80%
        /0.00000, 0.00000, 0.00000, 0.00002, 0.00011, 0.00068  &
        ,0.00424, 0.02545, 0.12815, 0.43942, 1.22370, 2.36600  &
        ,2.65160, 1.26560, 1.38750, 1.18780, 1.11420/

        DATA(qscat(2,6,ib,2),ib=1,17)  & ! Band 6, RH = 81%
        /0.00000, 0.00000, 0.00000, 0.00002, 0.00011, 0.00071  &
        ,0.00444, 0.02659, 0.13253, 0.45157, 1.24550, 2.38380  &
        ,2.62430, 1.24700, 1.37230, 1.17960, 1.11400/

        DATA(qscat(2,6,ib,3),ib=1,17)  & ! Band 6, RH = 82%
        /0.00000, 0.00000, 0.00000, 0.00002, 0.00012, 0.00075  &
        ,0.00466, 0.02783, 0.13723, 0.46478, 1.26870, 2.40190  &
        ,2.59320, 1.23150, 1.35570, 1.17110, 1.11360/

        DATA(qscat(2,6,ib,4),ib=1,17)  & ! Band 6, RH = 83%
        /0.00000, 0.00000, 0.00000, 0.00002, 0.00012, 0.00078  &
        ,0.00491, 0.02919, 0.14231, 0.47924, 1.29350, 2.41980  &
        ,2.55870, 1.21730, 1.33760, 1.16220, 1.11280/

        DATA(qscat(2,6,ib,5),ib=1,17)  & ! Band 6, RH = 84%
        /0.00000, 0.00000, 0.00000, 0.00002, 0.00013, 0.00083  &
        ,0.00518, 0.03071, 0.14783, 0.49521, 1.31980, 2.43720  &
        ,2.52170, 1.20200, 1.31800, 1.15320, 1.11180/

        DATA(qscat(2,6,ib,6),ib=1,17)  & ! Band 6, RH = 85%
        /0.00000, 0.00000, 0.00000, 0.00002, 0.00014, 0.00088  &
        ,0.00549, 0.03241, 0.15391, 0.51297, 1.34800, 2.45390  &
        ,2.48310, 1.18690, 1.29660, 1.14420, 1.11040/

        DATA(qscat(2,6,ib,7),ib=1,17)  & ! Band 6, RH = 86%
        /0.00000, 0.00000, 0.00000, 0.00003, 0.00015, 0.00094  &
        ,0.00584, 0.03435, 0.16065, 0.53286, 1.37810, 2.47020  &
        ,2.44160, 1.17650, 1.27370, 1.13550, 1.10870/

        DATA(qscat(2,6,ib,8),ib=1,17)  & ! Band 6, RH = 87%
        /0.00000, 0.00000, 0.00000, 0.00003, 0.00016, 0.00100  &
        ,0.00625, 0.03659, 0.16822, 0.55527, 1.41060, 2.48740  &
        ,2.39430, 1.17070, 1.24900, 1.12770, 1.10680/

        DATA(qscat(2,6,ib,9),ib=1,17)  & ! Band 6, RH = 88%
        /0.00000, 0.00000, 0.00000, 0.00003, 0.00017, 0.00108  &
        ,0.00673, 0.03920, 0.17680, 0.58066, 1.44600, 2.50700  &
        ,2.33790, 1.16440, 1.22330, 1.12120, 1.10500/

        DATA(qscat(2,6,ib,10),ib=1,17)  & ! Band 6, RH = 89%
        /0.00000, 0.00000, 0.00000, 0.00003, 0.00018, 0.00118  &
        ,0.00731, 0.04230, 0.18665, 0.60951, 1.48530, 2.53050  &
        ,2.27220, 1.16350, 1.19670, 1.11660, 1.10330/

        DATA(qscat(2,6,ib,11),ib=1,17)  & ! Band 6, RH = 90%
        /0.00000, 0.00000, 0.00001, 0.00004, 0.00020, 0.00129  &
        ,0.00802, 0.04602, 0.19811, 0.64239, 1.53030, 2.55720  &
        ,2.19910, 1.17110, 1.17100, 1.11460, 1.10200/

        DATA(qscat(2,6,ib,12),ib=1,17)  & ! Band 6, RH = 91%
        /0.00000, 0.00000, 0.00001, 0.00004, 0.00023, 0.00144  &
        ,0.00890, 0.05059, 0.21164, 0.67993, 1.58350, 2.58260  &
        ,2.11420, 1.18160, 1.14790, 1.11540, 1.10130/

        DATA(qscat(2,6,ib,13),ib=1,17)  & ! Band 6, RH = 92%
        /0.00000, 0.00000, 0.00001, 0.00004, 0.00025, 0.00162  &
        ,0.01002, 0.05630, 0.22793, 0.72299, 1.64830, 2.60080  &
        ,2.01080, 1.20470, 1.13040, 1.11820, 1.10070/

        DATA(qscat(2,6,ib,14),ib=1,17)  & ! Band 6, RH = 93%
        /0.00000, 0.00000, 0.00001, 0.00005, 0.00029, 0.00186  &
        ,0.01150, 0.06362, 0.24811, 0.77306, 1.72740, 2.61330  &
        ,1.89310, 1.23470, 1.12230, 1.12080, 1.09980/

        DATA(qscat(2,6,ib,15),ib=1,17)  & ! Band 6, RH = 94%
        /0.00000, 0.00000, 0.00001, 0.00006, 0.00034, 0.00220  &
        ,0.01352, 0.07331, 0.27412, 0.83336, 1.82040, 2.62670  &
        ,1.75080, 1.27270, 1.12720, 1.11960, 1.09830/

        DATA(qscat(2,6,ib,16),ib=1,17)  & ! Band 6, RH = 95%
        /0.00000, 0.00000, 0.00001, 0.00007, 0.00042, 0.00269  &
        ,0.01642, 0.08668, 0.30970, 0.91153, 1.92660, 2.61830  &
        ,1.58800, 1.31280, 1.14470, 1.11200, 1.09690/

        DATA(qscat(2,6,ib,17),ib=1,17)  & ! Band 6, RH = 96%
        /0.00000, 0.00000, 0.00001, 0.00009, 0.00054, 0.00345  &
        ,0.02094, 0.10614, 0.36288, 1.02410, 2.06130, 2.57390  &
        ,1.40530, 1.33520, 1.16200, 1.10310, 1.09550/

        DATA(qscat(2,6,ib,18),ib=1,17)  & ! Band 6, RH = 97%
        /0.00000, 0.00000, 0.00002, 0.00013, 0.00075, 0.00481  &
        ,0.02879, 0.13681, 0.45204, 1.19610, 2.23820, 2.45290  &
        ,1.22660, 1.29830, 1.14790, 1.10240, 1.09390/

        DATA(qscat(2,6,ib,19),ib=1,17)  & ! Band 6, RH = 98%
        /0.00000, 0.00000, 0.00003, 0.00020, 0.00121, 0.00776  &
        ,0.04513, 0.19164, 0.61562, 1.45610, 2.44340, 2.15890  &
        ,1.14680, 1.16680, 1.10830, 1.09820, 1.09190/

        DATA(qscat(2,6,ib,20),ib=1,17)  & ! Band 6, RH = 99%
        /0.00000, 0.00001, 0.00006, 0.00044, 0.00273, 0.01760  &
        ,0.09380, 0.32978, 0.95046, 1.95710, 2.54660, 1.47480  &
        ,1.30660, 1.14630, 1.10370, 1.09440, 1.08890/

        DATA(qscat(2,7,ib,1),ib=1,17)  & ! Band 7, RH = 80%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00006, 0.00038  &
        ,0.00236, 0.01426, 0.07517, 0.27122, 0.80264, 1.72410  &
        ,2.46210, 1.70540, 1.23030, 1.13150, 1.10550/

        DATA(qscat(2,7,ib,2),ib=1,17)  & ! Band 7, RH = 81%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00006, 0.00039  &
        ,0.00247, 0.01488, 0.07782, 0.27758, 0.81510, 1.73760  &
        ,2.45240, 1.67410, 1.23210, 1.13080, 1.10410/

        DATA(qscat(2,7,ib,3),ib=1,17)  & ! Band 7, RH = 82%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00007, 0.00041  &
        ,0.00259, 0.01556, 0.08066, 0.28450, 0.82862, 1.75220  &
        ,2.44260, 1.64180, 1.23350, 1.13020, 1.10280/

        DATA(qscat(2,7,ib,4),ib=1,17)  & ! Band 7, RH = 83%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00007, 0.00043  &
        ,0.00272, 0.01631, 0.08375, 0.29209, 0.84342, 1.76800  &
        ,2.43250, 1.60920, 1.23490, 1.12960, 1.10160/

        DATA(qscat(2,7,ib,5),ib=1,17)  & ! Band 7, RH = 84%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00007, 0.00046  &
        ,0.00286, 0.01714, 0.08713, 0.30051, 0.85977, 1.78540  &
        ,2.42190, 1.57600, 1.23620, 1.12900, 1.10050/

        DATA(qscat(2,7,ib,6),ib=1,17)  & ! Band 7, RH = 85%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00008, 0.00048  &
        ,0.00303, 0.01807, 0.09087, 0.30996, 0.87799, 1.80470  &
        ,2.40990, 1.54100, 1.23670, 1.12820, 1.09960/

        DATA(qscat(2,7,ib,7),ib=1,17)  & ! Band 7, RH = 86%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00008, 0.00052  &
        ,0.00322, 0.01914, 0.09504, 0.32068, 0.89846, 1.82600  &
        ,2.39550, 1.50360, 1.23680, 1.12730, 1.09870/

        DATA(qscat(2,7,ib,8),ib=1,17)  & ! Band 7, RH = 87%
        /0.00000, 0.00000, 0.00000, 0.00002, 0.00009, 0.00055  &
        ,0.00344, 0.02037, 0.09975, 0.33301, 0.92161, 1.84960  &
        ,2.37760, 1.46490, 1.23600, 1.12600, 1.09790/

        DATA(qscat(2,7,ib,9),ib=1,17)  & ! Band 7, RH = 88%
        /0.00000, 0.00000, 0.00000, 0.00002, 0.00009, 0.00059  &
        ,0.00370, 0.02181, 0.10512, 0.34734, 0.94794, 1.87540  &
        ,2.35620, 1.42550, 1.23380, 1.12430, 1.09730/

        DATA(qscat(2,7,ib,10),ib=1,17)  & ! Band 7, RH = 89%
        /0.00000, 0.00000, 0.00000, 0.00002, 0.00010, 0.00064  &
        ,0.00401, 0.02352, 0.11132, 0.36421, 0.97796, 1.90340  &
        ,2.33180, 1.38350, 1.23010, 1.12200, 1.09660/

        DATA(qscat(2,7,ib,11),ib=1,17)  & ! Band 7, RH = 90%
        /0.00000, 0.00000, 0.00000, 0.00002, 0.00011, 0.00071  &
        ,0.00439, 0.02559, 0.11856, 0.38428, 1.01230, 1.93400  &
        ,2.30340, 1.34000, 1.22390, 1.11890, 1.09590/

        DATA(qscat(2,7,ib,12),ib=1,17)  & ! Band 7, RH = 91%
        /0.00000, 0.00000, 0.00000, 0.00002, 0.00012, 0.00078  &
        ,0.00486, 0.02813, 0.12713, 0.40838, 1.05190, 1.96800  &
        ,2.26680, 1.29710, 1.21470, 1.11510, 1.09500/

        DATA(qscat(2,7,ib,13),ib=1,17)  & ! Band 7, RH = 92%
        /0.00000, 0.00000, 0.00000, 0.00002, 0.00014, 0.00088  &
        ,0.00547, 0.03133, 0.13744, 0.43759, 1.09800, 2.00720  &
        ,2.21910, 1.25350, 1.20180, 1.11070, 1.09380/

        DATA(qscat(2,7,ib,14),ib=1,17)  & ! Band 7, RH = 93%
        /0.00000, 0.00000, 0.00000, 0.00003, 0.00016, 0.00101  &
        ,0.00626, 0.03546, 0.15012, 0.47326, 1.15310, 2.05200  &
        ,2.16040, 1.21390, 1.18450, 1.10610, 1.09230/

        DATA(qscat(2,7,ib,15),ib=1,17)  & ! Band 7, RH = 94%
        /0.00000, 0.00000, 0.00000, 0.00003, 0.00019, 0.00119  &
        ,0.00735, 0.04098, 0.16617, 0.51722, 1.22160, 2.09940  &
        ,2.08070, 1.18030, 1.16300, 1.10240, 1.09100/

        DATA(qscat(2,7,ib,16),ib=1,17)  & ! Band 7, RH = 95%
        /0.00000, 0.00000, 0.00001, 0.00004, 0.00023, 0.00145  &
        ,0.00892, 0.04868, 0.18739, 0.57241, 1.30980, 2.14940  &
        ,1.97550, 1.15820, 1.13940, 1.10010, 1.08990/

        DATA(qscat(2,7,ib,17),ib=1,17)  & ! Band 7, RH = 96%
        /0.00000, 0.00000, 0.00001, 0.00005, 0.00029, 0.00186  &
        ,0.01136, 0.06010, 0.21759, 0.64532, 1.42370, 2.20050  &
        ,1.82860, 1.15440, 1.11910, 1.09840, 1.08870/

        DATA(qscat(2,7,ib,18),ib=1,17)  & ! Band 7, RH = 97%
        /0.00000, 0.00000, 0.00001, 0.00007, 0.00040, 0.00259  &
        ,0.01563, 0.07859, 0.26627, 0.75393, 1.57760, 2.23340  &
        ,1.62510, 1.17320, 1.11030, 1.09430, 1.08710/

        DATA(qscat(2,7,ib,19),ib=1,17)  & ! Band 7, RH = 98%
        /0.00000, 0.00000, 0.00002, 0.00011, 0.00065, 0.00417  &
        ,0.02463, 0.11278, 0.36187, 0.94196, 1.80160, 2.19720  &
        ,1.35710, 1.19040, 1.10810, 1.09080, 1.08520/

        DATA(qscat(2,7,ib,20),ib=1,17)  & ! Band 7, RH = 99%
        /0.00000, 0.00000, 0.00003, 0.00023, 0.00146, 0.00946  &
        ,0.05239, 0.19656, 0.59175, 1.32850, 2.11640, 1.86710  &
        ,1.14870, 1.12170, 1.09590, 1.08710, 1.08250/

        DATA(qscat(2,8,ib,1),ib=1,17)  & ! Band 8, RH = 80%
        /0.00000, 0.00000, 0.00002, 0.00014, 0.00080, 0.00507  &
        ,0.03067, 0.15390, 0.56158, 1.53360, 2.83190, 2.82510  &
        ,1.68210, 1.55970, 1.40520, 1.23670, 1.15300/

        DATA(qscat(2,8,ib,2),ib=1,17)  & ! Band 8, RH = 81%
        /0.00000, 0.00000, 0.00002, 0.00015, 0.00084, 0.00532  &
        ,0.03209, 0.15943, 0.57857, 1.56040, 2.85230, 2.78920  &
        ,1.67970, 1.55290, 1.39780, 1.23060, 1.15060/

        DATA(qscat(2,8,ib,3),ib=1,17)  & ! Band 8, RH = 82%
        /0.00000, 0.00000, 0.00002, 0.00016, 0.00089, 0.00559  &
        ,0.03365, 0.16536, 0.59667, 1.58870, 2.87350, 2.75320  &
        ,1.67690, 1.54710, 1.38900, 1.22480, 1.14790/

        DATA(qscat(2,8,ib,4),ib=1,17)  & ! Band 8, RH = 83%
        /0.00000, 0.00000, 0.00003, 0.00016, 0.00094, 0.00590  &
        ,0.03536, 0.17177, 0.61606, 1.61920, 2.89480, 2.71520  &
        ,1.67180, 1.54220, 1.37890, 1.21970, 1.14500/

        DATA(qscat(2,8,ib,5),ib=1,17)  & ! Band 8, RH = 84%
        /0.00000, 0.00000, 0.00003, 0.00017, 0.00099, 0.00623  &
        ,0.03727, 0.17878, 0.63697, 1.65220, 2.91570, 2.67260  &
        ,1.66840, 1.53760, 1.36690, 1.21540, 1.14180/

        DATA(qscat(2,8,ib,6),ib=1,17)  & ! Band 8, RH = 85%
        /0.00000, 0.00000, 0.00003, 0.00018, 0.00105, 0.00662  &
        ,0.03941, 0.18650, 0.65968, 1.68870, 2.93550, 2.62400  &
        ,1.66980, 1.53250, 1.35380, 1.21200, 1.13860/

        DATA(qscat(2,8,ib,7),ib=1,17)  & ! Band 8, RH = 86%
        /0.00000, 0.00000, 0.00003, 0.00020, 0.00112, 0.00705  &
        ,0.04184, 0.19511, 0.68454, 1.72920, 2.95400, 2.56950  &
        ,1.67150, 1.52920, 1.33940, 1.20930, 1.13570/

        DATA(qscat(2,8,ib,8),ib=1,17)  & ! Band 8, RH = 87%
        /0.00000, 0.00000, 0.00003, 0.00021, 0.00120, 0.00756  &
        ,0.04465, 0.20482, 0.71198, 1.77480, 2.97200, 2.51210  &
        ,1.67000, 1.52340, 1.32450, 1.20670, 1.13310/

        DATA(qscat(2,8,ib,9),ib=1,17)  & ! Band 8, RH = 88%
        /0.00000, 0.00001, 0.00003, 0.00023, 0.00130, 0.00816  &
        ,0.04793, 0.21591, 0.74258, 1.82600, 2.99080, 2.45020  &
        ,1.67190, 1.51810, 1.31020, 1.20300, 1.13050/

        DATA(qscat(2,8,ib,10),ib=1,17)  & ! Band 8, RH = 89%
        /0.00000, 0.00001, 0.00004, 0.00025, 0.00141, 0.00888  &
        ,0.05181, 0.22875, 0.77709, 1.88330, 3.01190, 2.37720  &
        ,1.67350, 1.50950, 1.29770, 1.19700, 1.12770/

        DATA(qscat(2,8,ib,11),ib=1,17)  & ! Band 8, RH = 90%
        /0.00000, 0.00001, 0.00004, 0.00027, 0.00155, 0.00976  &
        ,0.05649, 0.24387, 0.81663, 1.94700, 3.03290, 2.29370  &
        ,1.67410, 1.49720, 1.28840, 1.18830, 1.12490/

        DATA(qscat(2,8,ib,12),ib=1,17)  & ! Band 8, RH = 91%
        /0.00000, 0.00001, 0.00005, 0.00030, 0.00173, 0.01085  &
        ,0.06221, 0.26202, 0.86288, 2.01740, 3.04690, 2.20850  &
        ,1.66990, 1.47860, 1.28270, 1.17830, 1.12210/

        DATA(qscat(2,8,ib,13),ib=1,17)  & ! Band 8, RH = 92%
        /0.00000, 0.00001, 0.00005, 0.00034, 0.00195, 0.01224  &
        ,0.06936, 0.28434, 0.91848, 2.09650, 3.04930, 2.10790  &
        ,1.66450, 1.45260, 1.27850, 1.17020, 1.11870/

        DATA(qscat(2,8,ib,14),ib=1,17)  & ! Band 8, RH = 93%
        /0.00000, 0.00001, 0.00006, 0.00039, 0.00225, 0.01407  &
        ,0.07852, 0.31268, 0.98759, 2.18970, 3.04670, 2.00080  &
        ,1.64720, 1.41920, 1.27050, 1.16420, 1.11500/

        DATA(qscat(2,8,ib,15),ib=1,17)  & ! Band 8, RH = 94%
        /0.00000, 0.00001, 0.00007, 0.00046, 0.00266, 0.01657  &
        ,0.09061, 0.35017, 1.07630, 2.30470, 3.02940, 1.88620  &
        ,1.61820, 1.38250, 1.25390, 1.15650, 1.11200/

        DATA(qscat(2,8,ib,16),ib=1,17)  & ! Band 8, RH = 95%
        /0.00000, 0.00001, 0.00008, 0.00056, 0.00325, 0.02017  &
        ,0.10725, 0.40242, 1.19210, 2.44110, 2.97010, 1.77580  &
        ,1.57480, 1.35600, 1.23420, 1.14830, 1.10860/

        DATA(qscat(2,8,ib,17),ib=1,17)  & ! Band 8, RH = 96%
        /0.00000, 0.00002, 0.00011, 0.00072, 0.00417, 0.02575  &
        ,0.13143, 0.47960, 1.34390, 2.59430, 2.86590, 1.67950  &
        ,1.52410, 1.35150, 1.22040, 1.13830, 1.10490/

        DATA(qscat(2,8,ib,18),ib=1,17)  & ! Band 8, RH = 97%
        /0.00000, 0.00002, 0.00014, 0.00099, 0.00579, 0.03540  &
        ,0.16962, 0.59953, 1.55510, 2.78250, 2.65520, 1.63000  &
        ,1.48860, 1.33560, 1.19030, 1.12860, 1.10110/

        DATA(qscat(2,8,ib,19),ib=1,17)  & ! Band 8, RH = 98%
        /0.00000, 0.00003, 0.00022, 0.00157, 0.00928, 0.05539  &
        ,0.23956, 0.79819, 1.89870, 2.96230, 2.26500, 1.63460  &
        ,1.46420, 1.26460, 1.17180, 1.11640, 1.09710/

        DATA(qscat(2,8,ib,20),ib=1,17)  & ! Band 8, RH = 99%
        /0.00001, 0.00006, 0.00046, 0.00342, 0.02065, 0.11394  &
        ,0.43157, 1.25490, 2.49680, 2.88860, 1.70360, 1.52330  &
        ,1.33570, 1.21560, 1.13740, 1.10350, 1.09250/

! Qscat for Mineral Dust:
        DATA(qscat(3,1,ib,1),ib=1,17)  & ! Band 1, All RH's
        /0.00000, 0.00002, 0.00010, 0.00060, 0.00379, 0.02291  &
        ,0.11493, 0.42401, 1.14410, 2.13480, 2.40680, 2.33500  &
        ,2.10380, 1.61510, 1.47100, 1.27120, 1.15620/

        DATA(qscat(3,2,ib,1),ib=1,17)  & ! Band 2, All RH's
        /0.00001, 0.00063, 0.00402, 0.02557, 0.14955, 0.65760  &
        ,1.94660, 3.29530, 2.83420, 2.32070, 1.91820, 1.74640  &
        ,1.53580, 1.35690, 1.22550, 1.15300, 1.12440/

        DATA(qscat(3,3,ib,1),ib=1,17)  & ! Band 3, All RH's
        /0.00177, 0.01133, 0.06890, 0.34280, 1.23190, 2.81180  &
        ,3.55540, 2.36150, 2.06270, 1.79890, 1.62930, 1.44460  &
        ,1.29170, 1.19010, 1.14040, 1.12370, 1.11720/

        DATA(qscat(3,4,ib,1),ib=1,17)  & ! Band 4, All RH's
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000  &
        ,0.00003, 0.00020, 0.00125, 0.00808, 0.05258, 0.31342  &
        ,1.17700, 1.80920, 1.38670, 1.24760, 1.24330/

        DATA(qscat(3,5,ib,1),ib=1,17)  & ! Band 5, All RH's
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00002  &
        ,0.00010, 0.00066, 0.00420, 0.02663, 0.15649, 0.65345  &
        ,1.77580, 2.55670, 1.53740, 1.30830, 1.20600/

        DATA(qscat(3,6,ib,1),ib=1,17)  & ! Band 6, All RH's
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00001, 0.00009  &
        ,0.00054, 0.00344, 0.02155, 0.12273, 0.47129, 1.15740  &
        ,1.71170, 1.59490, 1.27650, 1.14810, 1.13330/

        DATA(qscat(3,7,ib,1),ib=1,17)  & ! Band 7, All RH's
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00002, 0.00015  &
        ,0.00094, 0.00603, 0.03906, 0.24079, 1.02410, 2.13950  &
        ,1.96750, 1.11300, 1.22940, 1.20200, 1.19660/

        DATA(qscat(3,8,ib,1),ib=1,17)  & ! Band 8, All RH's
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00003, 0.00018  &
        ,0.00112, 0.00676, 0.03661, 0.14457, 0.46121, 1.13810  &
        ,2.05970, 2.21310, 1.51910, 1.29570, 1.16240/

! Returning value for qscat:
        value = qscat(aerotype,radband,bin,rh)

return
END SUBROUTINE aeroqscat

!##############################################################################
Subroutine aerogasym (aerotype,radband,bin,rh,value)

implicit none

        INTEGER aerotype,radband,bin,rh,ib
        REAL gasym(3,8,17,20),value

! Setting for aerosol with no growth (completely insoluble):
        IF(AEROTYPE.EQ.3) rh = 1

! gasym for Ammonium Sulfate (with a 90% insoluble dust inclusion):
        DATA(gasym(1,1,ib,1),ib=1,17)  & ! Band 1, RH = 80%
        /0.00117, 0.00297, 0.00750, 0.01891, 0.04501, 0.11364  &
        ,0.30041, 0.62353, 0.76482, 0.83189, 0.82722, 0.84191  &
        ,0.86233, 0.88780, 0.92578, 0.94922, 0.96426/

        DATA(gasym(1,1,ib,2),ib=1,17)  & ! Band 1, RH = 81%
        /0.00118, 0.00298, 0.00754, 0.01902, 0.04528, 0.11433  &
        ,0.30236, 0.62443, 0.76550, 0.83229, 0.82734, 0.84191  &
        ,0.86158, 0.88833, 0.92577, 0.94911, 0.96418/

        DATA(gasym(1,1,ib,3),ib=1,17)  & ! Band 1, RH = 82%
        /0.00118, 0.00300, 0.00759, 0.01914, 0.04557, 0.11510  &
        ,0.30451, 0.62540, 0.76626, 0.83272, 0.82747, 0.84173  &
        ,0.86062, 0.88873, 0.92584, 0.94902, 0.96409/

        DATA(gasym(1,1,ib,4),ib=1,17)  & ! Band 1, RH = 83%
        /0.00119, 0.00302, 0.00764, 0.01927, 0.04589, 0.11595  &
        ,0.30691, 0.62645, 0.76709, 0.83317, 0.82761, 0.84128  &
        ,0.86007, 0.88931, 0.92590, 0.94883, 0.96397/

        DATA(gasym(1,1,ib,5),ib=1,17)  & ! Band 1, RH = 84%
        /0.00119, 0.00304, 0.00769, 0.01942, 0.04625, 0.11690  &
        ,0.30959, 0.62761, 0.76803, 0.83367, 0.82776, 0.84069  &
        ,0.86013, 0.88994, 0.92575, 0.94868, 0.96384/

        DATA(gasym(1,1,ib,6),ib=1,17)  & ! Band 1, RH = 85%
        /0.00120, 0.00306, 0.00775, 0.01958, 0.04666, 0.11797  &
        ,0.31262, 0.62889, 0.76907, 0.83420, 0.82792, 0.84021  &
        ,0.86009, 0.88998, 0.92561, 0.94848, 0.96368/

        DATA(gasym(1,1,ib,7),ib=1,17)  & ! Band 1, RH = 86%
        /0.00121, 0.00308, 0.00782, 0.01977, 0.04712, 0.11919  &
        ,0.31605, 0.63030, 0.77026, 0.83478, 0.82805, 0.83996  &
        ,0.85896, 0.89010, 0.92573, 0.94821, 0.96351/

        DATA(gasym(1,1,ib,8),ib=1,17)  & ! Band 1, RH = 87%
        /0.00122, 0.00311, 0.00790, 0.01998, 0.04764, 0.12058  &
        ,0.31998, 0.63188, 0.77161, 0.83539, 0.82812, 0.83997  &
        ,0.85767, 0.89090, 0.92560, 0.94795, 0.96330/

        DATA(gasym(1,1,ib,9),ib=1,17)  & ! Band 1, RH = 88%
        /0.00123, 0.00314, 0.00798, 0.02022, 0.04824, 0.12219  &
        ,0.32452, 0.63367, 0.77317, 0.83605, 0.82808, 0.83918  &
        ,0.85796, 0.89187, 0.92527, 0.94763, 0.96305/

        DATA(gasym(1,1,ib,10),ib=1,17)  & ! Band 1, RH = 89%
        /0.00124, 0.00317, 0.00809, 0.02050, 0.04895, 0.12407  &
        ,0.32984, 0.63570, 0.77497, 0.83676, 0.82783, 0.83799  &
        ,0.85685, 0.89159, 0.92530, 0.94720, 0.96277/

        DATA(gasym(1,1,ib,11),ib=1,17)  & ! Band 1, RH = 90%
        /0.00125, 0.00321, 0.00820, 0.02083, 0.04978, 0.12631  &
        ,0.33613, 0.63806, 0.77709, 0.83752, 0.82724, 0.83802  &
        ,0.85651, 0.89251, 0.92482, 0.94675, 0.96243/

        DATA(gasym(1,1,ib,12),ib=1,17)  & ! Band 1, RH = 91%
        /0.00127, 0.00326, 0.00834, 0.02122, 0.05078, 0.12901  &
        ,0.34372, 0.64084, 0.77961, 0.83834, 0.82621, 0.83744  &
        ,0.85565, 0.89289, 0.92456, 0.94623, 0.96204/

        DATA(gasym(1,1,ib,13),ib=1,17)  & ! Band 1, RH = 92%
        /0.00128, 0.00332, 0.00851, 0.02170, 0.05200, 0.13232  &
        ,0.35303, 0.64421, 0.78262, 0.83925, 0.82515, 0.83570  &
        ,0.85597, 0.89394, 0.92415, 0.94573, 0.96155/

        DATA(gasym(1,1,ib,14),ib=1,17)  & ! Band 1, RH = 93%
        /0.00130, 0.00338, 0.00872, 0.02229, 0.05352, 0.13651  &
        ,0.36471, 0.64843, 0.78628, 0.84035, 0.82455, 0.83452  &
        ,0.85580, 0.89337, 0.92354, 0.94496, 0.96094/

        DATA(gasym(1,1,ib,15),ib=1,17)  & ! Band 1, RH = 94%
        /0.00133, 0.00347, 0.00898, 0.02305, 0.05550, 0.14196  &
        ,0.37980, 0.65397, 0.79078, 0.84189, 0.82410, 0.83402  &
        ,0.85608, 0.89339, 0.92216, 0.94413, 0.96012/

        DATA(gasym(1,1,ib,16),ib=1,17)  & ! Band 1, RH = 95%
        /0.00136, 0.00357, 0.00932, 0.02405, 0.05814, 0.14936  &
        ,0.39998, 0.66174, 0.79638, 0.84413, 0.82178, 0.83185  &
        ,0.85500, 0.89460, 0.92144, 0.94287, 0.95907/

        DATA(gasym(1,1,ib,17),ib=1,17)  & ! Band 1, RH = 96%
        /0.00139, 0.00371, 0.00977, 0.02545, 0.06190, 0.16004  &
        ,0.42820, 0.67353, 0.80355, 0.84662, 0.81739, 0.82839  &
        ,0.85644, 0.89306, 0.91969, 0.94117, 0.95762/

        DATA(gasym(1,1,ib,18),ib=1,17)  & ! Band 1, RH = 97%
        /0.00144, 0.00389, 0.01043, 0.02754, 0.06767, 0.17687  &
        ,0.46979, 0.69307, 0.81322, 0.84775, 0.81114, 0.82726  &
        ,0.85801, 0.89402, 0.91828, 0.93885, 0.95543/

        DATA(gasym(1,1,ib,19),ib=1,17)  & ! Band 1, RH = 98%
        /0.00150, 0.00417, 0.01148, 0.03107, 0.07786, 0.20775  &
        ,0.53350, 0.72611, 0.82632, 0.84856, 0.80037, 0.82507  &
        ,0.86271, 0.89448, 0.91520, 0.93532, 0.95194/

        DATA(gasym(1,1,ib,20),ib=1,17)  & ! Band 1, RH = 99%
        /0.00159, 0.00462, 0.01348, 0.03859, 0.10155, 0.28522  &
        ,0.61997, 0.77368, 0.84552, 0.83754, 0.79296, 0.83100  &
        ,0.86785, 0.89273, 0.91179, 0.92955, 0.94555/

        DATA(gasym(1,2,ib,1),ib=1,17)  & ! Band 2, RH = 80%
        /0.00537, 0.01357, 0.03414, 0.08571, 0.20982, 0.51975  &
        ,0.67945, 0.74366, 0.69649, 0.72276, 0.77858, 0.83411  &
        ,0.87045, 0.90423, 0.92975, 0.94564, 0.95243/

        DATA(gasym(1,2,ib,2),ib=1,17)  & ! Band 2, RH = 81%
        /0.00539, 0.01363, 0.03432, 0.08619, 0.21107, 0.52133  &
        ,0.68060, 0.74461, 0.69749, 0.72297, 0.77893, 0.83415  &
        ,0.87043, 0.90414, 0.92969, 0.94569, 0.95257/

        DATA(gasym(1,2,ib,3),ib=1,17)  & ! Band 2, RH = 82%
        /0.00541, 0.01370, 0.03451, 0.08671, 0.21245, 0.52304  &
        ,0.68185, 0.74565, 0.69858, 0.72328, 0.77939, 0.83428  &
        ,0.87043, 0.90399, 0.92964, 0.94574, 0.95272/

        DATA(gasym(1,2,ib,4),ib=1,17)  & ! Band 2, RH = 83%
        /0.00544, 0.01378, 0.03473, 0.08729, 0.21398, 0.52492  &
        ,0.68322, 0.74678, 0.69978, 0.72365, 0.77991, 0.83449  &
        ,0.87040, 0.90389, 0.92959, 0.94579, 0.95289/

        DATA(gasym(1,2,ib,5),ib=1,17)  & ! Band 2, RH = 84%
        /0.00547, 0.01387, 0.03497, 0.08793, 0.21570, 0.52698  &
        ,0.68473, 0.74802, 0.70110, 0.72409, 0.78041, 0.83462  &
        ,0.87037, 0.90386, 0.92949, 0.94583, 0.95306/

        DATA(gasym(1,2,ib,6),ib=1,17)  & ! Band 2, RH = 85%
        /0.00550, 0.01396, 0.03523, 0.08866, 0.21763, 0.52926  &
        ,0.68640, 0.74939, 0.70255, 0.72472, 0.78112, 0.83478  &
        ,0.87037, 0.90372, 0.92942, 0.94589, 0.95325/

        DATA(gasym(1,2,ib,7),ib=1,17)  & ! Band 2, RH = 86%
        /0.00554, 0.01407, 0.03553, 0.08948, 0.21982, 0.53179  &
        ,0.68825, 0.75092, 0.70415, 0.72547, 0.78177, 0.83508  &
        ,0.87038, 0.90354, 0.92931, 0.94593, 0.95346/

        DATA(gasym(1,2,ib,8),ib=1,17)  & ! Band 2, RH = 87%
        /0.00558, 0.01419, 0.03587, 0.09041, 0.22233, 0.53461  &
        ,0.69032, 0.75262, 0.70589, 0.72630, 0.78261, 0.83521  &
        ,0.87035, 0.90345, 0.92919, 0.94598, 0.95368/

        DATA(gasym(1,2,ib,9),ib=1,17)  & ! Band 2, RH = 88%
        /0.00562, 0.01432, 0.03626, 0.09148, 0.22522, 0.53779  &
        ,0.69265, 0.75455, 0.70778, 0.72729, 0.78326, 0.83565  &
        ,0.87037, 0.90321, 0.92901, 0.94600, 0.95393/

        DATA(gasym(1,2,ib,10),ib=1,17)  & ! Band 2, RH = 89%
        /0.00567, 0.01448, 0.03671, 0.09273, 0.22860, 0.54138  &
        ,0.69531, 0.75675, 0.70979, 0.72817, 0.78424, 0.83590  &
        ,0.87040, 0.90304, 0.92886, 0.94602, 0.95419/

        DATA(gasym(1,2,ib,11),ib=1,17)  & ! Band 2, RH = 90%
        /0.00573, 0.01466, 0.03723, 0.09419, 0.23261, 0.54547  &
        ,0.69835, 0.75930, 0.71184, 0.72949, 0.78557, 0.83612  &
        ,0.87030, 0.90273, 0.92863, 0.94602, 0.95448/

        DATA(gasym(1,2,ib,12),ib=1,17)  & ! Band 2, RH = 91%
        /0.00579, 0.01486, 0.03785, 0.09595, 0.23743, 0.55019  &
        ,0.70189, 0.76229, 0.71394, 0.73083, 0.78731, 0.83647  &
        ,0.87034, 0.90242, 0.92832, 0.94598, 0.95479/

        DATA(gasym(1,2,ib,13),ib=1,17)  & ! Band 2, RH = 92%
        /0.00587, 0.01511, 0.03859, 0.09808, 0.24335, 0.55569  &
        ,0.70606, 0.76584, 0.71640, 0.73203, 0.78905, 0.83695  &
        ,0.87035, 0.90206, 0.92794, 0.94588, 0.95511/

        DATA(gasym(1,2,ib,14),ib=1,17)  & ! Band 2, RH = 93%
        /0.00595, 0.01541, 0.03951, 0.10073, 0.25078, 0.56219  &
        ,0.71109, 0.77007, 0.71960, 0.73402, 0.79098, 0.83778  &
        ,0.87030, 0.90148, 0.92734, 0.94571, 0.95545/

        DATA(gasym(1,2,ib,15),ib=1,17)  & ! Band 2, RH = 94%
        /0.00606, 0.01578, 0.04065, 0.10412, 0.26041, 0.57001  &
        ,0.71730, 0.77509, 0.72358, 0.73632, 0.79331, 0.83847  &
        ,0.87039, 0.90095, 0.92660, 0.94537, 0.95576/

        DATA(gasym(1,2,ib,16),ib=1,17)  & ! Band 2, RH = 95%
        /0.00619, 0.01624, 0.04215, 0.10863, 0.27337, 0.57968  &
        ,0.72525, 0.78093, 0.72828, 0.73877, 0.79636, 0.83940  &
        ,0.87055, 0.90006, 0.92562, 0.94483, 0.95602/

        DATA(gasym(1,2,ib,17),ib=1,17)  & ! Band 2, RH = 96%
        /0.00635, 0.01685, 0.04418, 0.11494, 0.29181, 0.59216  &
        ,0.73581, 0.78747, 0.73166, 0.74242, 0.80100, 0.84085  &
        ,0.87058, 0.89893, 0.92409, 0.94380, 0.95606/

        DATA(gasym(1,2,ib,18),ib=1,17)  & ! Band 2, RH = 97%
        /0.00656, 0.01768, 0.04711, 0.12443, 0.32012, 0.60965  &
        ,0.75025, 0.79526, 0.73673, 0.74748, 0.80696, 0.84210  &
        ,0.87099, 0.89742, 0.92179, 0.94185, 0.95557/

        DATA(gasym(1,2,ib,19),ib=1,17)  & ! Band 2, RH = 98%
        /0.00684, 0.01891, 0.05178, 0.14060, 0.36901, 0.63865  &
        ,0.77047, 0.80597, 0.73957, 0.75422, 0.81475, 0.84571  &
        ,0.87128, 0.89539, 0.91761, 0.93788, 0.95349/

        DATA(gasym(1,2,ib,20),ib=1,17)  & ! Band 2, RH = 99%
        /0.00725, 0.02096, 0.06071, 0.17559, 0.46942, 0.69948  &
        ,0.80194, 0.81321, 0.73877, 0.76821, 0.82788, 0.85142  &
        ,0.87170, 0.89043, 0.90992, 0.92852, 0.94582/

        DATA(gasym(1,3,ib,1),ib=1,17)  & ! Band 3, RH = 80%
        /0.02652, 0.06689, 0.16982, 0.38963, 0.59803, 0.70274  &
        ,0.72603, 0.64355, 0.74197, 0.79950, 0.84232, 0.88036  &
        ,0.91173, 0.93410, 0.94606, 0.95018, 0.95093/

        DATA(gasym(1,3,ib,2),ib=1,17)  & ! Band 3, RH = 81%
        /0.02664, 0.06724, 0.17082, 0.39114, 0.59937, 0.70378  &
        ,0.72705, 0.64464, 0.74209, 0.79992, 0.84241, 0.88032  &
        ,0.91169, 0.93412, 0.94617, 0.95037, 0.95114/

        DATA(gasym(1,3,ib,3),ib=1,17)  & ! Band 3, RH = 82%
        /0.02678, 0.06763, 0.17192, 0.39279, 0.60083, 0.70491  &
        ,0.72817, 0.64591, 0.74225, 0.80011, 0.84247, 0.88028  &
        ,0.91165, 0.93412, 0.94628, 0.95057, 0.95136/

        DATA(gasym(1,3,ib,4),ib=1,17)  & ! Band 3, RH = 83%
        /0.02692, 0.06806, 0.17313, 0.39461, 0.60244, 0.70616  &
        ,0.72936, 0.64728, 0.74254, 0.80045, 0.84253, 0.88019  &
        ,0.91153, 0.93410, 0.94641, 0.95078, 0.95161/

        DATA(gasym(1,3,ib,5),ib=1,17)  & ! Band 3, RH = 84%
        /0.02708, 0.06853, 0.17448, 0.39663, 0.60422, 0.70755  &
        ,0.73064, 0.64868, 0.74289, 0.80099, 0.84257, 0.88022  &
        ,0.91143, 0.93410, 0.94654, 0.95102, 0.95188/

        DATA(gasym(1,3,ib,6),ib=1,17)  & ! Band 3, RH = 85%
        /0.02726, 0.06905, 0.17598, 0.39888, 0.60619, 0.70909  &
        ,0.73201, 0.65017, 0.74311, 0.80141, 0.84258, 0.88020  &
        ,0.91135, 0.93406, 0.94668, 0.95128, 0.95218/

        DATA(gasym(1,3,ib,7),ib=1,17)  & ! Band 3, RH = 86%
        /0.02745, 0.06964, 0.17767, 0.40140, 0.60839, 0.71081  &
        ,0.73359, 0.65177, 0.74361, 0.80186, 0.84266, 0.88011  &
        ,0.91121, 0.93404, 0.94683, 0.95156, 0.95250/

        DATA(gasym(1,3,ib,8),ib=1,17)  & ! Band 3, RH = 87%
        /0.02767, 0.07030, 0.17958, 0.40425, 0.61086, 0.71274  &
        ,0.73548, 0.65356, 0.74409, 0.80209, 0.84280, 0.88002  &
        ,0.91107, 0.93397, 0.94698, 0.95187, 0.95287/

        DATA(gasym(1,3,ib,9),ib=1,17)  & ! Band 3, RH = 88%
        /0.02791, 0.07104, 0.18177, 0.40750, 0.61365, 0.71492  &
        ,0.73754, 0.65535, 0.74491, 0.80274, 0.84313, 0.87996  &
        ,0.91088, 0.93393, 0.94715, 0.95222, 0.95328/

        DATA(gasym(1,3,ib,10),ib=1,17)  & ! Band 3, RH = 89%
        /0.02819, 0.07190, 0.18429, 0.41125, 0.61683, 0.71741  &
        ,0.73966, 0.65761, 0.74578, 0.80359, 0.84342, 0.87982  &
        ,0.91070, 0.93383, 0.94732, 0.95260, 0.95373/

        DATA(gasym(1,3,ib,11),ib=1,17)  & ! Band 3, RH = 90%
        /0.02850, 0.07288, 0.18724, 0.41560, 0.62048, 0.72028  &
        ,0.74212, 0.65986, 0.74711, 0.80424, 0.84366, 0.87978  &
        ,0.91041, 0.93367, 0.94749, 0.95302, 0.95425/

        DATA(gasym(1,3,ib,12),ib=1,17)  & ! Band 3, RH = 91%
        /0.02885, 0.07403, 0.19073, 0.42073, 0.62472, 0.72362  &
        ,0.74487, 0.66221, 0.74900, 0.80509, 0.84388, 0.87959  &
        ,0.91005, 0.93349, 0.94764, 0.95349, 0.95484/

        DATA(gasym(1,3,ib,13),ib=1,17)  & ! Band 3, RH = 92%
        /0.02927, 0.07540, 0.19491, 0.42687, 0.62971, 0.72758  &
        ,0.74806, 0.66500, 0.75110, 0.80661, 0.84425, 0.87948  &
        ,0.90966, 0.93325, 0.94779, 0.95401, 0.95552/

        DATA(gasym(1,3,ib,14),ib=1,17)  & ! Band 3, RH = 93%
        /0.02976, 0.07705, 0.20005, 0.43436, 0.63563, 0.73234  &
        ,0.75176, 0.66810, 0.75336, 0.80775, 0.84493, 0.87937  &
        ,0.90922, 0.93283, 0.94789, 0.95459, 0.95630/

        DATA(gasym(1,3,ib,15),ib=1,17)  & ! Band 3, RH = 94%
        /0.03034, 0.07908, 0.20649, 0.44370, 0.64280, 0.73821  &
        ,0.75610, 0.67233, 0.75594, 0.80926, 0.84533, 0.87916  &
        ,0.90845, 0.93229, 0.94791, 0.95521, 0.95722/

        DATA(gasym(1,3,ib,16),ib=1,17)  & ! Band 3, RH = 95%
        /0.03105, 0.08164, 0.21483, 0.45569, 0.65166, 0.74557  &
        ,0.76131, 0.67693, 0.75913, 0.81125, 0.84655, 0.87886  &
        ,0.90758, 0.93139, 0.94775, 0.95586, 0.95829/

        DATA(gasym(1,3,ib,17),ib=1,17)  & ! Band 3, RH = 96%
        /0.03194, 0.08501, 0.22605, 0.47164, 0.66294, 0.75515  &
        ,0.76732, 0.68256, 0.76585, 0.81410, 0.84773, 0.87844  &
        ,0.90638, 0.93002, 0.94723, 0.95645, 0.95954/

        DATA(gasym(1,3,ib,18),ib=1,17)  & ! Band 3, RH = 97%
        /0.03309, 0.08963, 0.24195, 0.49387, 0.67800, 0.76753  &
        ,0.77344, 0.68949, 0.77355, 0.81770, 0.84944, 0.87792  &
        ,0.90439, 0.92765, 0.94589, 0.95673, 0.96093/

        DATA(gasym(1,3,ib,19),ib=1,17)  & ! Band 3, RH = 98%
        /0.03464, 0.09648, 0.26626, 0.52682, 0.70007, 0.78416  &
        ,0.77991, 0.69494, 0.78672, 0.82441, 0.85278, 0.87796  &
        ,0.90114, 0.92354, 0.94257, 0.95579, 0.96215/

        DATA(gasym(1,3,ib,20),ib=1,17)  & ! Band 3, RH = 99%
        /0.03689, 0.10790, 0.30789, 0.58123, 0.73588, 0.80519  &
        ,0.77954, 0.71019, 0.80913, 0.83581, 0.85567, 0.87534  &
        ,0.89539, 0.91457, 0.93322, 0.94957, 0.96066/

        DATA(gasym(1,4,ib,1),ib=1,17)  & ! Band 4, RH = 80%
        /0.00001, 0.00002, 0.00005, 0.00013, 0.00030, 0.00076  &
        ,0.00191, 0.00482, 0.01212, 0.03042, 0.07694, 0.20715  &
        ,0.50616, 0.68219, 0.77144, 0.85518, 0.88175/

        DATA(gasym(1,4,ib,2),ib=1,17)  & ! Band 4, RH = 81%
        /0.00001, 0.00002, 0.00005, 0.00013, 0.00030, 0.00076  &
        ,0.00192, 0.00484, 0.01218, 0.03059, 0.07737, 0.20845  &
        ,0.50735, 0.68348, 0.77275, 0.85571, 0.88217/

        DATA(gasym(1,4,ib,3),ib=1,17)  & ! Band 4, RH = 82%
        /0.00001, 0.00002, 0.00005, 0.00013, 0.00030, 0.00077  &
        ,0.00194, 0.00487, 0.01225, 0.03077, 0.07785, 0.20990  &
        ,0.50866, 0.68489, 0.77416, 0.85629, 0.88263/

        DATA(gasym(1,4,ib,4),ib=1,17)  & ! Band 4, RH = 83%
        /0.00001, 0.00002, 0.00005, 0.00013, 0.00031, 0.00077  &
        ,0.00195, 0.00491, 0.01234, 0.03098, 0.07839, 0.21152  &
        ,0.51013, 0.68644, 0.77571, 0.85692, 0.88314/

        DATA(gasym(1,4,ib,5),ib=1,17)  & ! Band 4, RH = 84%
        /0.00001, 0.00002, 0.00005, 0.00013, 0.00031, 0.00078  &
        ,0.00196, 0.00494, 0.01243, 0.03121, 0.07899, 0.21334  &
        ,0.51177, 0.68816, 0.77740, 0.85762, 0.88369/

        DATA(gasym(1,4,ib,6),ib=1,17)  & ! Band 4, RH = 85%
        /0.00001, 0.00002, 0.00005, 0.00013, 0.00031, 0.00078  &
        ,0.00198, 0.00498, 0.01253, 0.03147, 0.07968, 0.21541  &
        ,0.51362, 0.69006, 0.77925, 0.85839, 0.88430/

        DATA(gasym(1,4,ib,7),ib=1,17)  & ! Band 4, RH = 86%
        /0.00001, 0.00002, 0.00005, 0.00013, 0.00031, 0.00079  &
        ,0.00200, 0.00503, 0.01265, 0.03177, 0.08046, 0.21777  &
        ,0.51572, 0.69220, 0.78130, 0.85925, 0.88497/

        DATA(gasym(1,4,ib,8),ib=1,17)  & ! Band 4, RH = 87%
        /0.00001, 0.00002, 0.00005, 0.00013, 0.00032, 0.00080  &
        ,0.00202, 0.00508, 0.01278, 0.03211, 0.08136, 0.22049  &
        ,0.51814, 0.69461, 0.78359, 0.86022, 0.88572/

        DATA(gasym(1,4,ib,9),ib=1,17)  & ! Band 4, RH = 88%
        /0.00001, 0.00002, 0.00005, 0.00013, 0.00032, 0.00081  &
        ,0.00204, 0.00515, 0.01294, 0.03251, 0.08241, 0.22366  &
        ,0.52093, 0.69735, 0.78614, 0.86130, 0.88655/

        DATA(gasym(1,4,ib,10),ib=1,17)  & ! Band 4, RH = 89%
        /0.00001, 0.00002, 0.00005, 0.00013, 0.00032, 0.00082  &
        ,0.00207, 0.00522, 0.01313, 0.03298, 0.08365, 0.22741  &
        ,0.52421, 0.70050, 0.78902, 0.86254, 0.88750/

        DATA(gasym(1,4,ib,11),ib=1,17)  & ! Band 4, RH = 90%
        /0.00001, 0.00002, 0.00005, 0.00014, 0.00033, 0.00083  &
        ,0.00210, 0.00531, 0.01335, 0.03354, 0.08514, 0.23189  &
        ,0.52810, 0.70416, 0.79230, 0.86395, 0.88857/

        DATA(gasym(1,4,ib,12),ib=1,17)  & ! Band 4, RH = 91%
        /0.00001, 0.00002, 0.00005, 0.00014, 0.00034, 0.00085  &
        ,0.00215, 0.00541, 0.01362, 0.03422, 0.08695, 0.23737  &
        ,0.53280, 0.70848, 0.79607, 0.86559, 0.88980/

        DATA(gasym(1,4,ib,13),ib=1,17)  & ! Band 4, RH = 92%
        /0.00001, 0.00002, 0.00006, 0.00014, 0.00034, 0.00087  &
        ,0.00220, 0.00554, 0.01395, 0.03507, 0.08921, 0.24418  &
        ,0.53858, 0.71364, 0.80043, 0.86750, 0.89122/

        DATA(gasym(1,4,ib,14),ib=1,17)  & ! Band 4, RH = 93%
        /0.00001, 0.00002, 0.00006, 0.00015, 0.00035, 0.00089  &
        ,0.00226, 0.00571, 0.01438, 0.03616, 0.09210, 0.25290  &
        ,0.54588, 0.71991, 0.80553, 0.86977, 0.89288/

        DATA(gasym(1,4,ib,15),ib=1,17)  & ! Band 4, RH = 94%
        /0.00001, 0.00002, 0.00006, 0.00015, 0.00036, 0.00093  &
        ,0.00235, 0.00593, 0.01494, 0.03759, 0.09595, 0.26441  &
        ,0.55537, 0.72768, 0.81156, 0.87251, 0.89487/

        DATA(gasym(1,4,ib,16),ib=1,17)  & ! Band 4, RH = 95%
        /0.00001, 0.00002, 0.00006, 0.00016, 0.00038, 0.00097  &
        ,0.00246, 0.00623, 0.01571, 0.03957, 0.10130, 0.28029  &
        ,0.56815, 0.73748, 0.81876, 0.87588, 0.89727/

        DATA(gasym(1,4,ib,17),ib=1,17)  & ! Band 4, RH = 96%
        /0.00001, 0.00002, 0.00006, 0.00017, 0.00040, 0.00103  &
        ,0.00263, 0.00667, 0.01685, 0.04249, 0.10926, 0.30344  &
        ,0.58622, 0.75008, 0.82756, 0.88017, 0.90024/

        DATA(gasym(1,4,ib,18),ib=1,17)  & ! Band 4, RH = 97%
        /0.00001, 0.00003, 0.00007, 0.00018, 0.00044, 0.00113  &
        ,0.00290, 0.00738, 0.01868, 0.04722, 0.12238, 0.33989  &
        ,0.61329, 0.76667, 0.83880, 0.88576, 0.90404/

        DATA(gasym(1,4,ib,19),ib=1,17)  & ! Band 4, RH = 98%
        /0.00001, 0.00003, 0.00007, 0.00020, 0.00050, 0.00132  &
        ,0.00341, 0.00872, 0.02215, 0.05630, 0.14816, 0.40341  &
        ,0.65647, 0.78999, 0.85407, 0.89340, 0.90910/

        DATA(gasym(1,4,ib,20),ib=1,17)  & ! Band 4, RH = 99%
        /0.00001, 0.00003, 0.00009, 0.00025, 0.00065, 0.00177  &
        ,0.00470, 0.01222, 0.03142, 0.08120, 0.22183, 0.52633  &
        ,0.72828, 0.82645, 0.87706, 0.90476, 0.91638/

        DATA(gasym(1,5,ib,1),ib=1,17)  & ! Band 5, RH = 80%
        /0.00002, 0.00004, 0.00011, 0.00029, 0.00069, 0.00175  &
        ,0.00440, 0.01106, 0.02775, 0.06976, 0.18122, 0.49305  &
        ,0.69945, 0.80568, 0.82301, 0.91174, 0.92981/

        DATA(gasym(1,5,ib,2),ib=1,17)  & ! Band 5, RH = 81%
        /0.00002, 0.00005, 0.00011, 0.00029, 0.00070, 0.00175  &
        ,0.00442, 0.01111, 0.02789, 0.07013, 0.18234, 0.49524  &
        ,0.70166, 0.80776, 0.82698, 0.91278, 0.93043/

        DATA(gasym(1,5,ib,3),ib=1,17)  & ! Band 5, RH = 82%
        /0.00002, 0.00005, 0.00011, 0.00029, 0.00070, 0.00176  &
        ,0.00444, 0.01118, 0.02805, 0.07055, 0.18358, 0.49764  &
        ,0.70405, 0.80999, 0.83110, 0.91382, 0.93104/

        DATA(gasym(1,5,ib,4),ib=1,17)  & ! Band 5, RH = 83%
        /0.00002, 0.00005, 0.00012, 0.00029, 0.00070, 0.00177  &
        ,0.00447, 0.01124, 0.02823, 0.07101, 0.18496, 0.50029  &
        ,0.70664, 0.81236, 0.83537, 0.91489, 0.93166/

        DATA(gasym(1,5,ib,5),ib=1,17)  & ! Band 5, RH = 84%
        /0.00002, 0.00005, 0.00012, 0.00029, 0.00071, 0.00179  &
        ,0.00450, 0.01132, 0.02843, 0.07153, 0.18651, 0.50322  &
        ,0.70944, 0.81490, 0.83980, 0.91597, 0.93229/

        DATA(gasym(1,5,ib,6),ib=1,17)  & ! Band 5, RH = 85%
        /0.00002, 0.00005, 0.00012, 0.00030, 0.00071, 0.00180  &
        ,0.00453, 0.01141, 0.02865, 0.07212, 0.18827, 0.50648  &
        ,0.71250, 0.81763, 0.84441, 0.91707, 0.93291/

        DATA(gasym(1,5,ib,7),ib=1,17)  & ! Band 5, RH = 86%
        /0.00002, 0.00005, 0.00012, 0.00030, 0.00072, 0.00181  &
        ,0.00457, 0.01151, 0.02891, 0.07279, 0.19028, 0.51013  &
        ,0.71584, 0.82056, 0.84919, 0.91819, 0.93355/

        DATA(gasym(1,5,ib,8),ib=1,17)  & ! Band 5, RH = 87%
        /0.00002, 0.00005, 0.00012, 0.00030, 0.00073, 0.00183  &
        ,0.00462, 0.01162, 0.02920, 0.07356, 0.19260, 0.51425  &
        ,0.71952, 0.82372, 0.85414, 0.91933, 0.93419/

        DATA(gasym(1,5,ib,9),ib=1,17)  & ! Band 5, RH = 88%
        /0.00002, 0.00005, 0.00012, 0.00030, 0.00073, 0.00185  &
        ,0.00467, 0.01176, 0.02954, 0.07446, 0.19529, 0.51893  &
        ,0.72358, 0.82713, 0.85929, 0.92050, 0.93486/

        DATA(gasym(1,5,ib,10),ib=1,17)  & ! Band 5, RH = 89%
        /0.00002, 0.00005, 0.00012, 0.00031, 0.00074, 0.00188  &
        ,0.00473, 0.01191, 0.02994, 0.07552, 0.19847, 0.52431  &
        ,0.72809, 0.83083, 0.86463, 0.92169, 0.93554/

        DATA(gasym(1,5,ib,11),ib=1,17)  & ! Band 5, RH = 90%
        /0.00002, 0.00005, 0.00012, 0.00031, 0.00075, 0.00190  &
        ,0.00480, 0.01210, 0.03042, 0.07679, 0.20228, 0.53053  &
        ,0.73314, 0.83485, 0.87016, 0.92292, 0.93625/

        DATA(gasym(1,5,ib,12),ib=1,17)  & ! Band 5, RH = 91%
        /0.00002, 0.00005, 0.00012, 0.00032, 0.00077, 0.00194  &
        ,0.00489, 0.01233, 0.03100, 0.07834, 0.20692, 0.53783  &
        ,0.73883, 0.83924, 0.87590, 0.92420, 0.93700/

        DATA(gasym(1,5,ib,13),ib=1,17)  & ! Band 5, RH = 92%
        /0.00002, 0.00005, 0.00013, 0.00032, 0.00078, 0.00198  &
        ,0.00500, 0.01261, 0.03172, 0.08026, 0.21271, 0.54652  &
        ,0.74533, 0.84404, 0.88183, 0.92555, 0.93779/

        DATA(gasym(1,5,ib,14),ib=1,17)  & ! Band 5, RH = 93%
        /0.00002, 0.00005, 0.00013, 0.00033, 0.00080, 0.00203  &
        ,0.00514, 0.01296, 0.03263, 0.08271, 0.22012, 0.55702  &
        ,0.75283, 0.84933, 0.88794, 0.92697, 0.93863/

        DATA(gasym(1,5,ib,15),ib=1,17)  & ! Band 5, RH = 94%
        /0.00002, 0.00005, 0.00013, 0.00034, 0.00083, 0.00210  &
        ,0.00532, 0.01343, 0.03385, 0.08597, 0.22994, 0.56997  &
        ,0.76162, 0.85521, 0.89421, 0.92851, 0.93955/

        DATA(gasym(1,5,ib,16),ib=1,17)  & ! Band 5, RH = 95%
        /0.00002, 0.00005, 0.00014, 0.00035, 0.00086, 0.00220  &
        ,0.00557, 0.01408, 0.03552, 0.09049, 0.24358, 0.58634  &
        ,0.77216, 0.86187, 0.90057, 0.93019, 0.94054/

        DATA(gasym(1,5,ib,17),ib=1,17)  & ! Band 5, RH = 96%
        /0.00002, 0.00005, 0.00014, 0.00037, 0.00091, 0.00233  &
        ,0.00593, 0.01502, 0.03798, 0.09721, 0.26377, 0.60767  &
        ,0.78515, 0.86958, 0.90698, 0.93206, 0.94162/

        DATA(gasym(1,5,ib,18),ib=1,17)  & ! Band 5, RH = 97%
        /0.00002, 0.00006, 0.00015, 0.00040, 0.00099, 0.00255  &
        ,0.00651, 0.01654, 0.04197, 0.10825, 0.29657, 0.63655  &
        ,0.80174, 0.87873, 0.91350, 0.93414, 0.94276/

        DATA(gasym(1,5,ib,19),ib=1,17)  & ! Band 5, RH = 98%
        /0.00002, 0.00006, 0.00017, 0.00045, 0.00112, 0.00294  &
        ,0.00759, 0.01942, 0.04962, 0.12986, 0.35812, 0.67803  &
        ,0.82399, 0.88994, 0.92023, 0.93648, 0.94390/

        DATA(gasym(1,5,ib,20),ib=1,17)  & ! Band 5, RH = 99%
        /0.00002, 0.00007, 0.00019, 0.00055, 0.00145, 0.00392  &
        ,0.01037, 0.02704, 0.07051, 0.19134, 0.50142, 0.74725  &
        ,0.85771, 0.90532, 0.92795, 0.93927, 0.94494/

        DATA(gasym(1,6,ib,1),ib=1,17)  & ! Band 6, RH = 80%
        /0.00004, 0.00011, 0.00029, 0.00072, 0.00173, 0.00436  &
        ,0.01099, 0.02767, 0.07017, 0.18631, 0.51579, 0.74556  &
        ,0.86027, 0.91201, 0.93629, 0.95085, 0.95889/

        DATA(gasym(1,6,ib,2),ib=1,17)  & ! Band 6, RH = 81%
        /0.00004, 0.00011, 0.00029, 0.00073, 0.00174, 0.00439  &
        ,0.01105, 0.02784, 0.07062, 0.18757, 0.51852, 0.74668  &
        ,0.86088, 0.91227, 0.93639, 0.95103, 0.95908/

        DATA(gasym(1,6,ib,3),ib=1,17)  & ! Band 6, RH = 82%
        /0.00005, 0.00011, 0.00029, 0.00073, 0.00175, 0.00442  &
        ,0.01113, 0.02803, 0.07111, 0.18896, 0.52153, 0.74791  &
        ,0.86155, 0.91255, 0.93650, 0.95122, 0.95928/

        DATA(gasym(1,6,ib,4),ib=1,17)  & ! Band 6, RH = 83%
        /0.00005, 0.00012, 0.00029, 0.00074, 0.00176, 0.00445  &
        ,0.01121, 0.02825, 0.07167, 0.19051, 0.52485, 0.74928  &
        ,0.86229, 0.91286, 0.93660, 0.95144, 0.95950/

        DATA(gasym(1,6,ib,5),ib=1,17)  & ! Band 6, RH = 84%
        /0.00005, 0.00012, 0.00029, 0.00074, 0.00178, 0.00449  &
        ,0.01131, 0.02849, 0.07229, 0.19226, 0.52853, 0.75081  &
        ,0.86311, 0.91319, 0.93671, 0.95168, 0.95975/

        DATA(gasym(1,6,ib,6),ib=1,17)  & ! Band 6, RH = 85%
        /0.00005, 0.00012, 0.00030, 0.00075, 0.00179, 0.00453  &
        ,0.01141, 0.02876, 0.07300, 0.19423, 0.53265, 0.75253  &
        ,0.86402, 0.91355, 0.93683, 0.95194, 0.96001/

        DATA(gasym(1,6,ib,7),ib=1,17)  & ! Band 6, RH = 86%
        /0.00005, 0.00012, 0.00030, 0.00076, 0.00181, 0.00458  &
        ,0.01154, 0.02907, 0.07380, 0.19648, 0.53728, 0.75447  &
        ,0.86503, 0.91395, 0.93694, 0.95224, 0.96031/

        DATA(gasym(1,6,ib,8),ib=1,17)  & ! Band 6, RH = 87%
        /0.00005, 0.00012, 0.00030, 0.00077, 0.00183, 0.00463  &
        ,0.01168, 0.02943, 0.07472, 0.19907, 0.54252, 0.75669  &
        ,0.86618, 0.91439, 0.93706, 0.95257, 0.96063/

        DATA(gasym(1,6,ib,9),ib=1,17)  & ! Band 6, RH = 88%
        /0.00005, 0.00012, 0.00031, 0.00078, 0.00186, 0.00469  &
        ,0.01184, 0.02984, 0.07580, 0.20209, 0.54850, 0.75924  &
        ,0.86748, 0.91487, 0.93718, 0.95294, 0.96099/

        DATA(gasym(1,6,ib,10),ib=1,17)  & ! Band 6, RH = 89%
        /0.00005, 0.00012, 0.00031, 0.00079, 0.00188, 0.00477  &
        ,0.01203, 0.03033, 0.07706, 0.20564, 0.55537, 0.76221  &
        ,0.86897, 0.91541, 0.93731, 0.95337, 0.96140/

        DATA(gasym(1,6,ib,11),ib=1,17)  & ! Band 6, RH = 90%
        /0.00005, 0.00012, 0.00031, 0.00080, 0.00192, 0.00485  &
        ,0.01225, 0.03091, 0.07856, 0.20989, 0.56333, 0.76571  &
        ,0.87070, 0.91601, 0.93745, 0.95388, 0.96186/

        DATA(gasym(1,6,ib,12),ib=1,17)  & ! Band 6, RH = 91%
        /0.00005, 0.00012, 0.00032, 0.00082, 0.00196, 0.00496  &
        ,0.01252, 0.03161, 0.08039, 0.21506, 0.57268, 0.76988  &
        ,0.87273, 0.91667, 0.93758, 0.95448, 0.96239/

        DATA(gasym(1,6,ib,13),ib=1,17)  & ! Band 6, RH = 92%
        /0.00005, 0.00013, 0.00033, 0.00084, 0.00201, 0.00509  &
        ,0.01286, 0.03248, 0.08265, 0.22150, 0.58374, 0.77493  &
        ,0.87511, 0.91740, 0.93766, 0.95518, 0.96299/

        DATA(gasym(1,6,ib,14),ib=1,17)  & ! Band 6, RH = 93%
        /0.00005, 0.00013, 0.00033, 0.00086, 0.00207, 0.00525  &
        ,0.01329, 0.03358, 0.08554, 0.22972, 0.59699, 0.78112  &
        ,0.87795, 0.91822, 0.93760, 0.95604, 0.96370/

        DATA(gasym(1,6,ib,15),ib=1,17)  & ! Band 6, RH = 94%
        /0.00005, 0.00013, 0.00035, 0.00089, 0.00215, 0.00546  &
        ,0.01385, 0.03503, 0.08935, 0.24061, 0.61300, 0.78881  &
        ,0.88128, 0.91911, 0.93741, 0.95712, 0.96453/

        DATA(gasym(1,6,ib,16),ib=1,17)  & ! Band 6, RH = 95%
        /0.00005, 0.00014, 0.00036, 0.00093, 0.00225, 0.00575  &
        ,0.01461, 0.03702, 0.09459, 0.25570, 0.63242, 0.79840  &
        ,0.88503, 0.91994, 0.93727, 0.95848, 0.96553/

        DATA(gasym(1,6,ib,17),ib=1,17)  & ! Band 6, RH = 96%
        /0.00005, 0.00014, 0.00038, 0.00099, 0.00241, 0.00617  &
        ,0.01572, 0.03992, 0.10229, 0.27802, 0.65566, 0.81017  &
        ,0.88885, 0.92037, 0.93673, 0.96024, 0.96674/

        DATA(gasym(1,6,ib,18),ib=1,17)  & ! Band 6, RH = 97%
        /0.00006, 0.00015, 0.00040, 0.00107, 0.00264, 0.00683  &
        ,0.01749, 0.04457, 0.11476, 0.31438, 0.68192, 0.82394  &
        ,0.89319, 0.92098, 0.93630, 0.96249, 0.96823/

        DATA(gasym(1,6,ib,19),ib=1,17)  & ! Band 6, RH = 98%
        /0.00006, 0.00016, 0.00045, 0.00122, 0.00306, 0.00802  &
        ,0.02075, 0.05329, 0.13856, 0.38365, 0.70776, 0.84055  &
        ,0.90334, 0.91878, 0.93739, 0.96510, 0.97000/

        DATA(gasym(1,6,ib,20),ib=1,17)  & ! Band 6, RH = 99%
        /0.00006, 0.00018, 0.00053, 0.00152, 0.00402, 0.01095  &
        ,0.02908, 0.07617, 0.20364, 0.55160, 0.75796, 0.87080  &
        ,0.91192, 0.90909, 0.94921, 0.96868, 0.97261/

        DATA(gasym(1,7,ib,1),ib=1,17)  & ! Band 7, RH = 80%
        /0.00004, 0.00011, 0.00028, 0.00071, 0.00171, 0.00431  &
        ,0.01083, 0.02720, 0.06855, 0.18178, 0.51389, 0.69553  &
        ,0.77061, 0.84653, 0.89565, 0.90847, 0.91468/

        DATA(gasym(1,7,ib,2),ib=1,17)  & ! Band 7, RH = 81%
        /0.00004, 0.00011, 0.00028, 0.00071, 0.00172, 0.00433  &
        ,0.01088, 0.02732, 0.06886, 0.18266, 0.51574, 0.69711  &
        ,0.77197, 0.84659, 0.89649, 0.90914, 0.91533/

        DATA(gasym(1,7,ib,3),ib=1,17)  & ! Band 7, RH = 82%
        /0.00004, 0.00011, 0.00028, 0.00071, 0.00172, 0.00435  &
        ,0.01093, 0.02745, 0.06921, 0.18364, 0.51778, 0.69884  &
        ,0.77346, 0.84667, 0.89741, 0.90988, 0.91604/

        DATA(gasym(1,7,ib,4),ib=1,17)  & ! Band 7, RH = 83%
        /0.00004, 0.00011, 0.00028, 0.00072, 0.00173, 0.00437  &
        ,0.01099, 0.02760, 0.06959, 0.18473, 0.52005, 0.70073  &
        ,0.77511, 0.84675, 0.89841, 0.91069, 0.91682/

        DATA(gasym(1,7,ib,5),ib=1,17)  & ! Band 7, RH = 84%
        /0.00004, 0.00011, 0.00028, 0.00072, 0.00174, 0.00439  &
        ,0.01106, 0.02777, 0.07003, 0.18597, 0.52256, 0.70282  &
        ,0.77693, 0.84685, 0.89951, 0.91157, 0.91767/

        DATA(gasym(1,7,ib,6),ib=1,17)  & ! Band 7, RH = 85%
        /0.00004, 0.00011, 0.00029, 0.00072, 0.00175, 0.00442  &
        ,0.01113, 0.02796, 0.07052, 0.18737, 0.52539, 0.70514  &
        ,0.77895, 0.84696, 0.90072, 0.91254, 0.91861/

        DATA(gasym(1,7,ib,7),ib=1,17)  & ! Band 7, RH = 86%
        /0.00004, 0.00011, 0.00029, 0.00073, 0.00177, 0.00446  &
        ,0.01122, 0.02818, 0.07108, 0.18897, 0.52857, 0.70772  &
        ,0.78121, 0.84711, 0.90206, 0.91362, 0.91965/

        DATA(gasym(1,7,ib,8),ib=1,17)  & ! Band 7, RH = 87%
        /0.00004, 0.00011, 0.00029, 0.00074, 0.00178, 0.00450  &
        ,0.01132, 0.02843, 0.07173, 0.19083, 0.53219, 0.71061  &
        ,0.78376, 0.84729, 0.90355, 0.91483, 0.92082/

        DATA(gasym(1,7,ib,9),ib=1,17)  & ! Band 7, RH = 88%
        /0.00005, 0.00012, 0.00029, 0.00074, 0.00180, 0.00454  &
        ,0.01143, 0.02873, 0.07249, 0.19300, 0.53633, 0.71388  &
        ,0.78664, 0.84752, 0.90521, 0.91618, 0.92212/

        DATA(gasym(1,7,ib,10),ib=1,17)  & ! Band 7, RH = 89%
        /0.00005, 0.00012, 0.00030, 0.00075, 0.00182, 0.00459  &
        ,0.01157, 0.02907, 0.07338, 0.19558, 0.54113, 0.71759  &
        ,0.78993, 0.84783, 0.90707, 0.91770, 0.92359/

        DATA(gasym(1,7,ib,11),ib=1,17)  & ! Band 7, RH = 90%
        /0.00005, 0.00012, 0.00030, 0.00076, 0.00184, 0.00466  &
        ,0.01173, 0.02949, 0.07446, 0.19867, 0.54674, 0.72185  &
        ,0.79372, 0.84825, 0.90916, 0.91944, 0.92527/

        DATA(gasym(1,7,ib,12),ib=1,17)  & ! Band 7, RH = 91%
        /0.00005, 0.00012, 0.00030, 0.00077, 0.00187, 0.00473  &
        ,0.01193, 0.02999, 0.07577, 0.20247, 0.55337, 0.72679  &
        ,0.79811, 0.84885, 0.91152, 0.92143, 0.92720/

        DATA(gasym(1,7,ib,13),ib=1,17)  & ! Band 7, RH = 92%
        /0.00005, 0.00012, 0.00031, 0.00079, 0.00191, 0.00483  &
        ,0.01217, 0.03062, 0.07741, 0.20723, 0.56133, 0.73259  &
        ,0.80326, 0.84969, 0.91420, 0.92375, 0.92943/

        DATA(gasym(1,7,ib,14),ib=1,17)  & ! Band 7, RH = 93%
        /0.00005, 0.00012, 0.00031, 0.00080, 0.00195, 0.00495  &
        ,0.01249, 0.03143, 0.07952, 0.21336, 0.57102, 0.73947  &
        ,0.80934, 0.85094, 0.91725, 0.92646, 0.93206/

        DATA(gasym(1,7,ib,15),ib=1,17)  & ! Band 7, RH = 94%
        /0.00005, 0.00012, 0.00032, 0.00083, 0.00201, 0.00510  &
        ,0.01290, 0.03250, 0.08232, 0.22156, 0.58302, 0.74777  &
        ,0.81658, 0.85283, 0.92072, 0.92969, 0.93518/

        DATA(gasym(1,7,ib,16),ib=1,17)  & ! Band 7, RH = 95%
        /0.00005, 0.00013, 0.00033, 0.00086, 0.00209, 0.00532  &
        ,0.01347, 0.03398, 0.08622, 0.23303, 0.59815, 0.75794  &
        ,0.82528, 0.85581, 0.92466, 0.93360, 0.93896/

        DATA(gasym(1,7,ib,17),ib=1,17)  & ! Band 7, RH = 96%
        /0.00005, 0.00013, 0.00035, 0.00090, 0.00221, 0.00563  &
        ,0.01431, 0.03616, 0.09201, 0.25020, 0.61750, 0.77066  &
        ,0.83583, 0.86073, 0.92905, 0.93844, 0.94362/

        DATA(gasym(1,7,ib,18),ib=1,17)  & ! Band 7, RH = 97%
        /0.00005, 0.00014, 0.00037, 0.00096, 0.00238, 0.00613  &
        ,0.01565, 0.03971, 0.10151, 0.27857, 0.64237, 0.78703  &
        ,0.84920, 0.86884, 0.93399, 0.94459, 0.94949/

        DATA(gasym(1,7,ib,19),ib=1,17)  & ! Band 7, RH = 98%
        /0.00005, 0.00015, 0.00040, 0.00108, 0.00270, 0.00706  &
        ,0.01818, 0.04646, 0.11991, 0.33373, 0.67384, 0.80995  &
        ,0.86888, 0.88109, 0.93985, 0.95262, 0.95704/

        DATA(gasym(1,7,ib,20),ib=1,17)  & ! Band 7, RH = 99%
        /0.00006, 0.00016, 0.00046, 0.00132, 0.00346, 0.00936  &
        ,0.02475, 0.06451, 0.17095, 0.47769, 0.72551, 0.85047  &
        ,0.89604, 0.90077, 0.94941, 0.96342, 0.96696/

        DATA(gasym(1,8,ib,1),ib=1,17)  & ! Band 8, RH = 80%
        /0.00012, 0.00031, 0.00078, 0.00198, 0.00470, 0.01184  &
        ,0.02979, 0.07506, 0.19257, 0.49463, 0.73032, 0.84980  &
        ,0.89773, 0.89442, 0.90417, 0.94218, 0.96453/

        DATA(gasym(1,8,ib,2),ib=1,17)  & ! Band 8, RH = 81%
        /0.00012, 0.00031, 0.00079, 0.00199, 0.00473, 0.01192  &
        ,0.02999, 0.07557, 0.19396, 0.49757, 0.73161, 0.85033  &
        ,0.89766, 0.89377, 0.90422, 0.94234, 0.96452/

        DATA(gasym(1,8,ib,3),ib=1,17)  & ! Band 8, RH = 82%
        /0.00012, 0.00031, 0.00079, 0.00200, 0.00476, 0.01200  &
        ,0.03022, 0.07615, 0.19550, 0.50079, 0.73304, 0.85090  &
        ,0.89758, 0.89305, 0.90424, 0.94252, 0.96452/

        DATA(gasym(1,8,ib,4),ib=1,17)  & ! Band 8, RH = 83%
        /0.00012, 0.00031, 0.00080, 0.00202, 0.00480, 0.01210  &
        ,0.03047, 0.07678, 0.19722, 0.50434, 0.73461, 0.85152  &
        ,0.89748, 0.89225, 0.90424, 0.94273, 0.96451/

        DATA(gasym(1,8,ib,5),ib=1,17)  & ! Band 8, RH = 84%
        /0.00012, 0.00032, 0.00080, 0.00204, 0.00484, 0.01221  &
        ,0.03075, 0.07750, 0.19914, 0.50826, 0.73636, 0.85219  &
        ,0.89738, 0.89136, 0.90425, 0.94295, 0.96451/

        DATA(gasym(1,8,ib,6),ib=1,17)  & ! Band 8, RH = 85%
        /0.00013, 0.00032, 0.00081, 0.00206, 0.00489, 0.01233  &
        ,0.03106, 0.07830, 0.20131, 0.51263, 0.73832, 0.85294  &
        ,0.89726, 0.89039, 0.90430, 0.94322, 0.96452/

        DATA(gasym(1,8,ib,7),ib=1,17)  & ! Band 8, RH = 86%
        /0.00013, 0.00032, 0.00082, 0.00208, 0.00494, 0.01247  &
        ,0.03142, 0.07922, 0.20378, 0.51751, 0.74051, 0.85376  &
        ,0.89713, 0.88932, 0.90443, 0.94354, 0.96454/

        DATA(gasym(1,8,ib,8),ib=1,17)  & ! Band 8, RH = 87%
        /0.00013, 0.00033, 0.00083, 0.00210, 0.00500, 0.01263  &
        ,0.03183, 0.08027, 0.20662, 0.52299, 0.74298, 0.85465  &
        ,0.89699, 0.88815, 0.90463, 0.94392, 0.96457/

        DATA(gasym(1,8,ib,9),ib=1,17)  & ! Band 8, RH = 88%
        /0.00013, 0.00033, 0.00084, 0.00213, 0.00507, 0.01282  &
        ,0.03230, 0.08148, 0.20992, 0.52919, 0.74579, 0.85564  &
        ,0.89683, 0.88682, 0.90484, 0.94434, 0.96463/

        DATA(gasym(1,8,ib,10),ib=1,17)  & ! Band 8, RH = 89%
        /0.00013, 0.00033, 0.00085, 0.00216, 0.00516, 0.01303  &
        ,0.03285, 0.08291, 0.21379, 0.53626, 0.74898, 0.85671  &
        ,0.89663, 0.88526, 0.90491, 0.94478, 0.96473/

        DATA(gasym(1,8,ib,11),ib=1,17)  & ! Band 8, RH = 90%
        /0.00013, 0.00034, 0.00086, 0.00220, 0.00525, 0.01329  &
        ,0.03351, 0.08460, 0.21841, 0.54437, 0.75263, 0.85788  &
        ,0.89638, 0.88330, 0.90494, 0.94527, 0.96486/

        DATA(gasym(1,8,ib,12),ib=1,17)  & ! Band 8, RH = 91%
        /0.00013, 0.00034, 0.00088, 0.00225, 0.00537, 0.01359  &
        ,0.03430, 0.08665, 0.22402, 0.55374, 0.75682, 0.85915  &
        ,0.89605, 0.88076, 0.90516, 0.94577, 0.96504/

        DATA(gasym(1,8,ib,13),ib=1,17)  & ! Band 8, RH = 92%
        /0.00013, 0.00035, 0.00090, 0.00230, 0.00551, 0.01397  &
        ,0.03528, 0.08918, 0.23096, 0.56464, 0.76164, 0.86050  &
        ,0.89554, 0.87763, 0.90546, 0.94625, 0.96525/

        DATA(gasym(1,8,ib,14),ib=1,17)  & ! Band 8, RH = 93%
        /0.00014, 0.00036, 0.00092, 0.00237, 0.00569, 0.01444  &
        ,0.03651, 0.09238, 0.23979, 0.57741, 0.76716, 0.86196  &
        ,0.89478, 0.87422, 0.90558, 0.94662, 0.96545/

        DATA(gasym(1,8,ib,15),ib=1,17)  & ! Band 8, RH = 94%
        /0.00014, 0.00037, 0.00095, 0.00246, 0.00592, 0.01506  &
        ,0.03812, 0.09657, 0.25141, 0.59239, 0.77349, 0.86363  &
        ,0.89369, 0.87011, 0.90613, 0.94682, 0.96555/

        DATA(gasym(1,8,ib,16),ib=1,17)  & ! Band 8, RH = 95%
        /0.00014, 0.00038, 0.00099, 0.00257, 0.00623, 0.01589  &
        ,0.04030, 0.10230, 0.26741, 0.60994, 0.78075, 0.86586  &
        ,0.89205, 0.86463, 0.90693, 0.94688, 0.96547/

        DATA(gasym(1,8,ib,17),ib=1,17)  & ! Band 8, RH = 96%
        /0.00015, 0.00039, 0.00104, 0.00274, 0.00666, 0.01708  &
        ,0.04346, 0.11065, 0.29084, 0.63023, 0.78934, 0.86939  &
        ,0.88866, 0.85896, 0.90863, 0.94720, 0.96543/

        DATA(gasym(1,8,ib,18),ib=1,17)  & ! Band 8, RH = 97%
        /0.00015, 0.00042, 0.00112, 0.00298, 0.00733, 0.01894  &
        ,0.04844, 0.12399, 0.32847, 0.65354, 0.80092, 0.87401  &
        ,0.88307, 0.85464, 0.91303, 0.94895, 0.96599/

        DATA(gasym(1,8,ib,19),ib=1,17)  & ! Band 8, RH = 98%
        /0.00016, 0.00045, 0.00124, 0.00338, 0.00851, 0.02229  &
        ,0.05760, 0.14906, 0.39807, 0.68417, 0.81977, 0.87745  &
        ,0.87024, 0.85966, 0.92344, 0.95183, 0.96660/

        DATA(gasym(1,8,ib,20),ib=1,17)  & ! Band 8, RH = 99%
        /0.00017, 0.00050, 0.00147, 0.00425, 0.01122, 0.03047  &
        ,0.08085, 0.21585, 0.54999, 0.75096, 0.84944, 0.87763  &
        ,0.83675, 0.88566, 0.93352, 0.95632, 0.96844/

! gasym for Sea Salt:
        DATA(gasym(2,1,ib,1),ib=1,17)  & ! Band 1, RH = 80%
        /0.00372, 0.00962, 0.02456, 0.06213, 0.14883, 0.39219  &
        ,0.65438, 0.78972, 0.83671, 0.78626, 0.78165, 0.83022  &
        ,0.85802, 0.86761, 0.88382, 0.89449, 0.90512/

        DATA(gasym(2,1,ib,2),ib=1,17)  & ! Band 1, RH = 81%
        /0.00383, 0.00989, 0.02521, 0.06379, 0.15294, 0.40271  &
        ,0.65885, 0.79259, 0.83781, 0.78397, 0.78391, 0.83075  &
        ,0.85349, 0.86989, 0.88407, 0.89752, 0.90520/

        DATA(gasym(2,1,ib,3),ib=1,17)  & ! Band 1, RH = 82%
        /0.00394, 0.01016, 0.02590, 0.06557, 0.15735, 0.41380  &
        ,0.66368, 0.79553, 0.83882, 0.78226, 0.78828, 0.83500  &
        ,0.85327, 0.87120, 0.88561, 0.89556, 0.90657/

        DATA(gasym(2,1,ib,4),ib=1,17)  & ! Band 1, RH = 83%
        /0.00405, 0.01044, 0.02665, 0.06747, 0.16210, 0.42554  &
        ,0.66893, 0.79854, 0.83966, 0.78107, 0.78868, 0.83865  &
        ,0.86204, 0.87445, 0.88581, 0.89632, 0.90781/

        DATA(gasym(2,1,ib,5),ib=1,17)  & ! Band 1, RH = 84%
        /0.00417, 0.01075, 0.02745, 0.06955, 0.16728, 0.43802  &
        ,0.67468, 0.80166, 0.84029, 0.77988, 0.79020, 0.83869  &
        ,0.85811, 0.87449, 0.88442, 0.89631, 0.90763/

        DATA(gasym(2,1,ib,6),ib=1,17)  & ! Band 1, RH = 85%
        /0.00429, 0.01109, 0.02833, 0.07182, 0.17299, 0.45137  &
        ,0.68102, 0.80492, 0.84070, 0.77778, 0.79472, 0.84174  &
        ,0.85738, 0.87524, 0.88720, 0.89890, 0.90901/

        DATA(gasym(2,1,ib,7),ib=1,17)  & ! Band 1, RH = 86%
        /0.00443, 0.01146, 0.02930, 0.07433, 0.17934, 0.46568  &
        ,0.68803, 0.80832, 0.84098, 0.77429, 0.79622, 0.83775  &
        ,0.85839, 0.87731, 0.88912, 0.90004, 0.90992/

        DATA(gasym(2,1,ib,8),ib=1,17)  & ! Band 1, RH = 87%
        /0.00458, 0.01187, 0.03038, 0.07715, 0.18650, 0.48106  &
        ,0.69577, 0.81187, 0.84136, 0.77204, 0.79778, 0.84188  &
        ,0.85927, 0.87715, 0.88897, 0.90059, 0.91080/

        DATA(gasym(2,1,ib,9),ib=1,17)  & ! Band 1, RH = 88%
        /0.00474, 0.01233, 0.03160, 0.08035, 0.19467, 0.49758  &
        ,0.70426, 0.81554, 0.84186, 0.77074, 0.79876, 0.84134  &
        ,0.86155, 0.87740, 0.88984, 0.90067, 0.91104/

        DATA(gasym(2,1,ib,10),ib=1,17)  & ! Band 1, RH = 89%
        /0.00493, 0.01285, 0.03300, 0.08403, 0.20411, 0.51526  &
        ,0.71348, 0.81924, 0.84218, 0.76760, 0.80246, 0.83760  &
        ,0.86210, 0.87828, 0.89180, 0.90254, 0.91217/

        DATA(gasym(2,1,ib,11),ib=1,17)  & ! Band 1, RH = 90%
        /0.00515, 0.01346, 0.03463, 0.08832, 0.21519, 0.53402  &
        ,0.72331, 0.82292, 0.84183, 0.76457, 0.80603, 0.84276  &
        ,0.86685, 0.88147, 0.89038, 0.90348, 0.91312/

        DATA(gasym(2,1,ib,12),ib=1,17)  & ! Band 1, RH = 91%
        /0.00540, 0.01417, 0.03655, 0.09340, 0.22839, 0.55361  &
        ,0.73358, 0.82659, 0.84074, 0.76428, 0.80602, 0.84499  &
        ,0.86246, 0.88025, 0.89539, 0.90508, 0.91456/

        DATA(gasym(2,1,ib,13),ib=1,17)  & ! Band 1, RH = 92%
        /0.00570, 0.01502, 0.03885, 0.09953, 0.24443, 0.57363  &
        ,0.74409, 0.83046, 0.83972, 0.76163, 0.81173, 0.84697  &
        ,0.86430, 0.88460, 0.89655, 0.90649, 0.91561/

        DATA(gasym(2,1,ib,14),ib=1,17)  & ! Band 1, RH = 93%
        /0.00606, 0.01605, 0.04167, 0.10709, 0.26435, 0.59349  &
        ,0.75476, 0.83499, 0.83755, 0.76377, 0.81183, 0.84940  &
        ,0.86876, 0.88303, 0.89586, 0.90712, 0.91671/

        DATA(gasym(2,1,ib,15),ib=1,17)  & ! Band 1, RH = 94%
        /0.00651, 0.01735, 0.04522, 0.11669, 0.28982, 0.61262  &
        ,0.76583, 0.84017, 0.83335, 0.76662, 0.81802, 0.85219  &
        ,0.87428, 0.88685, 0.89767, 0.90890, 0.91803/

        DATA(gasym(2,1,ib,16),ib=1,17)  & ! Band 1, RH = 95%
        /0.00707, 0.01902, 0.04987, 0.12941, 0.32355, 0.63113  &
        ,0.77813, 0.84456, 0.82795, 0.77394, 0.82352, 0.85770  &
        ,0.87246, 0.88918, 0.89986, 0.91026, 0.91965/

        DATA(gasym(2,1,ib,17),ib=1,17)  & ! Band 1, RH = 96%
        /0.00782, 0.02128, 0.05627, 0.14722, 0.37024, 0.65138  &
        ,0.79257, 0.84813, 0.81822, 0.78486, 0.82798, 0.86211  &
        ,0.87779, 0.89235, 0.90157, 0.91211, 0.92155/

        DATA(gasym(2,1,ib,18),ib=1,17)  & ! Band 1, RH = 97%
        /0.00888, 0.02459, 0.06582, 0.17448, 0.43773, 0.68085  &
        ,0.80980, 0.85182, 0.80225, 0.80110, 0.84028, 0.86074  &
        ,0.87985, 0.89430, 0.90479, 0.91487, 0.92395/

        DATA(gasym(2,1,ib,19),ib=1,17)  & ! Band 1, RH = 98%
        /0.01054, 0.03001, 0.08207, 0.22284, 0.53416, 0.72946  &
        ,0.82963, 0.85041, 0.77904, 0.81410, 0.84773, 0.86815  &
        ,0.88704, 0.89885, 0.90826, 0.91792, 0.92692/

        DATA(gasym(2,1,ib,20),ib=1,17)  & ! Band 1, RH = 99%
        /0.01363, 0.04120, 0.11833, 0.33813, 0.63653, 0.78644  &
        ,0.84970, 0.82810, 0.78299, 0.82584, 0.86127, 0.87790  &
        ,0.89090, 0.90419, 0.91355, 0.92296, 0.93167/

        DATA(gasym(2,2,ib,1),ib=1,17)  & ! Band 2, RH = 80%
        /0.01683, 0.04331, 0.11050, 0.28850, 0.58369, 0.74438  &
        ,0.81836, 0.80447, 0.73055, 0.77459, 0.81733, 0.83274  &
        ,0.84863, 0.85917, 0.86360, 0.86843, 0.87242/

        DATA(gasym(2,2,ib,2),ib=1,17)  & ! Band 2, RH = 81%
        /0.01733, 0.04452, 0.11345, 0.29658, 0.58950, 0.74770  &
        ,0.82000, 0.80393, 0.73114, 0.77681, 0.81572, 0.83395  &
        ,0.84796, 0.85772, 0.86497, 0.86979, 0.87284/

        DATA(gasym(2,2,ib,3),ib=1,17)  & ! Band 2, RH = 82%
        /0.01781, 0.04572, 0.11661, 0.30519, 0.59529, 0.75111  &
        ,0.82154, 0.80320, 0.73219, 0.78002, 0.81807, 0.83258  &
        ,0.84681, 0.85543, 0.86543, 0.86948, 0.87340/

        DATA(gasym(2,2,ib,4),ib=1,17)  & ! Band 2, RH = 83%
        /0.01830, 0.04701, 0.11999, 0.31444, 0.60112, 0.75462  &
        ,0.82297, 0.80204, 0.73364, 0.78163, 0.82102, 0.83552  &
        ,0.85075, 0.86037, 0.86432, 0.87073, 0.87386/

        DATA(gasym(2,2,ib,5),ib=1,17)  & ! Band 2, RH = 84%
        /0.01882, 0.04839, 0.12366, 0.32443, 0.60705, 0.75828  &
        ,0.82429, 0.80027, 0.73598, 0.78498, 0.81995, 0.83651  &
        ,0.84904, 0.86013, 0.86590, 0.87096, 0.87431/

        DATA(gasym(2,2,ib,6),ib=1,17)  & ! Band 2, RH = 85%
        /0.01938, 0.04989, 0.12767, 0.33533, 0.61318, 0.76211  &
        ,0.82555, 0.79818, 0.73709, 0.79020, 0.82072, 0.83782  &
        ,0.85120, 0.85892, 0.86730, 0.87137, 0.87490/

        DATA(gasym(2,2,ib,7),ib=1,17)  & ! Band 2, RH = 86%
        /0.01999, 0.05154, 0.13211, 0.34732, 0.61962, 0.76613  &
        ,0.82680, 0.79634, 0.73924, 0.79358, 0.82072, 0.83703  &
        ,0.85285, 0.86039, 0.86798, 0.87179, 0.87546/

        DATA(gasym(2,2,ib,8),ib=1,17)  & ! Band 2, RH = 87%
        /0.02067, 0.05338, 0.13709, 0.36063, 0.62650, 0.77037  &
        ,0.82814, 0.79455, 0.74193, 0.79697, 0.82168, 0.84015  &
        ,0.85140, 0.86180, 0.86834, 0.87270, 0.87592/

        DATA(gasym(2,2,ib,9),ib=1,17)  & ! Band 2, RH = 88%
        /0.02142, 0.05545, 0.14272, 0.37552, 0.63402, 0.77485  &
        ,0.82964, 0.79204, 0.74373, 0.80062, 0.82417, 0.84145  &
        ,0.85418, 0.86192, 0.86840, 0.87302, 0.87660/

        DATA(gasym(2,2,ib,10),ib=1,17)  & ! Band 2, RH = 89%
        /0.02228, 0.05780, 0.14920, 0.39230, 0.64240, 0.77962  &
        ,0.83120, 0.78817, 0.74632, 0.80321, 0.82706, 0.84247  &
        ,0.85506, 0.86377, 0.86879, 0.87357, 0.87728/

        DATA(gasym(2,2,ib,11),ib=1,17)  & ! Band 2, RH = 90%
        /0.02325, 0.06052, 0.15673, 0.41130, 0.65193, 0.78473  &
        ,0.83252, 0.78453, 0.74917, 0.80812, 0.82854, 0.84418  &
        ,0.85751, 0.86304, 0.86906, 0.87449, 0.87800/

        DATA(gasym(2,2,ib,12),ib=1,17)  & ! Band 2, RH = 91%
        /0.02439, 0.06372, 0.16566, 0.43288, 0.66296, 0.79028  &
        ,0.83332, 0.78048, 0.75208, 0.80862, 0.83105, 0.84615  &
        ,0.85732, 0.86441, 0.87112, 0.87513, 0.87846/

        DATA(gasym(2,2,ib,13),ib=1,17)  & ! Band 2, RH = 92%
        /0.02574, 0.06753, 0.17642, 0.45740, 0.67579, 0.79633  &
        ,0.83375, 0.77445, 0.75429, 0.81214, 0.83277, 0.84482  &
        ,0.85741, 0.86581, 0.87219, 0.87577, 0.87908/

        DATA(gasym(2,2,ib,14),ib=1,17)  & ! Band 2, RH = 93%
        /0.02736, 0.07218, 0.18969, 0.48504, 0.69060, 0.80288  &
        ,0.83432, 0.76893, 0.75687, 0.81547, 0.83550, 0.84979  &
        ,0.85766, 0.86814, 0.87204, 0.87633, 0.87988/

        DATA(gasym(2,2,ib,15),ib=1,17)  & ! Band 2, RH = 94%
        /0.02936, 0.07799, 0.20652, 0.51569, 0.70721, 0.80974  &
        ,0.83399, 0.76095, 0.76034, 0.81779, 0.83762, 0.85008  &
        ,0.85854, 0.86772, 0.87350, 0.87779, 0.88050/

        DATA(gasym(2,2,ib,16),ib=1,17)  & ! Band 2, RH = 95%
        /0.03190, 0.08551, 0.22871, 0.54865, 0.72499, 0.81699  &
        ,0.83222, 0.75349, 0.76463, 0.81971, 0.83948, 0.85158  &
        ,0.86297, 0.87029, 0.87469, 0.87837, 0.88158/

        DATA(gasym(2,2,ib,17),ib=1,17)  & ! Band 2, RH = 96%
        /0.03528, 0.09574, 0.25954, 0.58256, 0.74341, 0.82556  &
        ,0.82763, 0.74827, 0.77301, 0.82388, 0.84146, 0.85594  &
        ,0.86309, 0.87077, 0.87590, 0.87981, 0.88261/

        DATA(gasym(2,2,ib,18),ib=1,17)  & ! Band 2, RH = 97%
        /0.04004, 0.11070, 0.30572, 0.61709, 0.76374, 0.83377  &
        ,0.81822, 0.74813, 0.78878, 0.82730, 0.84369, 0.85811  &
        ,0.86754, 0.87232, 0.87738, 0.88085, 0.88399/

        DATA(gasym(2,2,ib,19),ib=1,17)  & ! Band 2, RH = 98%
        /0.04746, 0.13545, 0.38245, 0.66046, 0.78895, 0.84005  &
        ,0.79758, 0.75669, 0.81058, 0.83376, 0.84920, 0.86142  &
        ,0.86888, 0.87447, 0.87917, 0.88238, 0.88522/

        DATA(gasym(2,2,ib,20),ib=1,17)  & ! Band 2, RH = 99%
        /0.06134, 0.18747, 0.51997, 0.73291, 0.82151, 0.83610  &
        ,0.75735, 0.77150, 0.82576, 0.84278, 0.85340, 0.86711  &
        ,0.87323, 0.87765, 0.88149, 0.88446, 0.88713/

        DATA(gasym(2,3,ib,1),ib=1,17)  & ! Band 3, RH = 80%
        /0.09028, 0.23390, 0.46961, 0.68055, 0.77385, 0.81365  &
        ,0.76262, 0.72097, 0.80397, 0.82291, 0.83335, 0.84851  &
        ,0.85630, 0.86067, 0.86496, 0.86700, 0.86879/

        DATA(gasym(2,3,ib,2),ib=1,17)  & ! Band 3, RH = 81%
        /0.09305, 0.24036, 0.47702, 0.68471, 0.77665, 0.81460  &
        ,0.76139, 0.72276, 0.80423, 0.82000, 0.83609, 0.84798  &
        ,0.85289, 0.86159, 0.86407, 0.86804, 0.86939/

        DATA(gasym(2,3,ib,3),ib=1,17)  & ! Band 3, RH = 82%
        /0.09572, 0.24669, 0.48469, 0.68890, 0.77955, 0.81541  &
        ,0.75979, 0.72599, 0.80798, 0.82118, 0.83633, 0.84952  &
        ,0.85691, 0.86219, 0.86498, 0.86836, 0.86999/

        DATA(gasym(2,3,ib,4),ib=1,17)  & ! Band 3, RH = 83%
        /0.09841, 0.25336, 0.49265, 0.69317, 0.78244, 0.81583  &
        ,0.75812, 0.72909, 0.80724, 0.82354, 0.83865, 0.84906  &
        ,0.85590, 0.86321, 0.86591, 0.86857, 0.87033/

        DATA(gasym(2,3,ib,5),ib=1,17)  & ! Band 3, RH = 84%
        /0.10129, 0.26041, 0.50094, 0.69756, 0.78538, 0.81662  &
        ,0.75656, 0.73421, 0.80790, 0.82365, 0.84079, 0.85161  &
        ,0.85758, 0.86258, 0.86634, 0.86915, 0.87093/

        DATA(gasym(2,3,ib,6),ib=1,17)  & ! Band 3, RH = 85%
        /0.10440, 0.26793, 0.50964, 0.70212, 0.78843, 0.81709  &
        ,0.75477, 0.73853, 0.80949, 0.82281, 0.84233, 0.84842  &
        ,0.85731, 0.86385, 0.86704, 0.86945, 0.87143/

        DATA(gasym(2,3,ib,7),ib=1,17)  & ! Band 3, RH = 86%
        /0.10780, 0.27598, 0.51878, 0.70693, 0.79167, 0.81744  &
        ,0.75302, 0.74210, 0.80940, 0.82438, 0.84202, 0.84991  &
        ,0.85835, 0.86463, 0.86746, 0.86990, 0.87171/

        DATA(gasym(2,3,ib,8),ib=1,17)  & ! Band 3, RH = 87%
        /0.11154, 0.28468, 0.52844, 0.71205, 0.79497, 0.81768  &
        ,0.75058, 0.74726, 0.80966, 0.82402, 0.84153, 0.85015  &
        ,0.85866, 0.86446, 0.86818, 0.87067, 0.87227/

        DATA(gasym(2,3,ib,9),ib=1,17)  & ! Band 3, RH = 88%
        /0.11573, 0.29412, 0.53869, 0.71758, 0.79852, 0.81734  &
        ,0.74694, 0.75228, 0.81027, 0.82536, 0.84402, 0.85250  &
        ,0.86083, 0.86555, 0.86870, 0.87120, 0.87285/

        DATA(gasym(2,3,ib,10),ib=1,17)  & ! Band 3, RH = 89%
        /0.12046, 0.30444, 0.54964, 0.72354, 0.80224, 0.81720  &
        ,0.74357, 0.75964, 0.81048, 0.82748, 0.84634, 0.85461  &
        ,0.86036, 0.86708, 0.86947, 0.87178, 0.87323/

        DATA(gasym(2,3,ib,11),ib=1,17)  & ! Band 3, RH = 90%
        /0.12590, 0.31578, 0.56142, 0.72983, 0.80579, 0.81699  &
        ,0.73982, 0.76480, 0.81142, 0.82817, 0.84298, 0.85640  &
        ,0.86317, 0.86679, 0.87007, 0.87243, 0.87392/

        DATA(gasym(2,3,ib,12),ib=1,17)  & ! Band 3, RH = 91%
        /0.13222, 0.32835, 0.57429, 0.73649, 0.80936, 0.81594  &
        ,0.73664, 0.77042, 0.81320, 0.82963, 0.84684, 0.85505  &
        ,0.86363, 0.86739, 0.87048, 0.87279, 0.87444/

        DATA(gasym(2,3,ib,13),ib=1,17)  & ! Band 3, RH = 92%
        /0.13971, 0.34241, 0.58859, 0.74353, 0.81298, 0.81400  &
        ,0.73250, 0.77908, 0.81474, 0.83287, 0.84788, 0.85759  &
        ,0.86347, 0.86788, 0.87134, 0.87374, 0.87492/

        DATA(gasym(2,3,ib,14),ib=1,17)  & ! Band 3, RH = 93%
        /0.14874, 0.35842, 0.60473, 0.75111, 0.81661, 0.81155  &
        ,0.72798, 0.78519, 0.81856, 0.83556, 0.84886, 0.85778  &
        ,0.86367, 0.86874, 0.87203, 0.87436, 0.87543/

        DATA(gasym(2,3,ib,15),ib=1,17)  & ! Band 3, RH = 94%
        /0.15990, 0.37718, 0.62322, 0.75973, 0.82003, 0.80769  &
        ,0.72424, 0.79287, 0.82129, 0.83943, 0.85214, 0.86049  &
        ,0.86633, 0.87005, 0.87296, 0.87494, 0.87633/

        DATA(gasym(2,3,ib,16),ib=1,17)  & ! Band 3, RH = 95%
        /0.17410, 0.40022, 0.64446, 0.76933, 0.82368, 0.80243  &
        ,0.72191, 0.80250, 0.82666, 0.84324, 0.85335, 0.86023  &
        ,0.86645, 0.87120, 0.87388, 0.87565, 0.87709/

        DATA(gasym(2,3,ib,17),ib=1,17)  & ! Band 3, RH = 96%
        /0.19289, 0.43017, 0.66811, 0.78087, 0.82737, 0.79334  &
        ,0.72616, 0.81112, 0.83311, 0.84324, 0.85502, 0.86432  &
        ,0.86835, 0.87185, 0.87464, 0.87662, 0.87780/

        DATA(gasym(2,3,ib,18),ib=1,17)  & ! Band 3, RH = 97%
        /0.21911, 0.47083, 0.69366, 0.79592, 0.82979, 0.77907  &
        ,0.73832, 0.81622, 0.83326, 0.84596, 0.85760, 0.86539  &
        ,0.87008, 0.87354, 0.87586, 0.87724, 0.87857/

        DATA(gasym(2,3,ib,19),ib=1,17)  & ! Band 3, RH = 98%
        /0.25817, 0.52647, 0.72383, 0.81342, 0.82758, 0.75498  &
        ,0.76788, 0.81513, 0.83501, 0.84981, 0.85871, 0.86711  &
        ,0.87169, 0.87480, 0.87680, 0.87850, 0.87944/

        DATA(gasym(2,3,ib,20),ib=1,17)  & ! Band 3, RH = 99%
        /0.32054, 0.60325, 0.76351, 0.82916, 0.80845, 0.72768  &
        ,0.80522, 0.83202, 0.84412, 0.85774, 0.86561, 0.86980  &
        ,0.87372, 0.87640, 0.87834, 0.87969, 0.88053/

        DATA(gasym(2,4,ib,1),ib=1,17)  & ! Band 4, RH = 80%
        /0.00002, 0.00006, 0.00015, 0.00039, 0.00093, 0.00236  &
        ,0.00596, 0.01501, 0.03787, 0.09733, 0.26651, 0.58059  &
        ,0.76358, 0.85096, 0.89449, 0.91873, 0.92950/

        DATA(gasym(2,4,ib,2),ib=1,17)  & ! Band 4, RH = 81%
        /0.00002, 0.00006, 0.00016, 0.00040, 0.00096, 0.00243  &
        ,0.00612, 0.01543, 0.03893, 0.10018, 0.27471, 0.58770  &
        ,0.76704, 0.85241, 0.89521, 0.91890, 0.92941/

        DATA(gasym(2,4,ib,3),ib=1,17)  & ! Band 4, RH = 82%
        /0.00002, 0.00006, 0.00016, 0.00041, 0.00099, 0.00250  &
        ,0.00630, 0.01587, 0.04006, 0.10323, 0.28346, 0.59493  &
        ,0.77054, 0.85388, 0.89595, 0.91907, 0.92933/

        DATA(gasym(2,4,ib,4),ib=1,17)  & ! Band 4, RH = 83%
        /0.00003, 0.00007, 0.00017, 0.00043, 0.00102, 0.00257  &
        ,0.00649, 0.01635, 0.04128, 0.10654, 0.29287, 0.60231  &
        ,0.77411, 0.85540, 0.89672, 0.91925, 0.92925/

        DATA(gasym(2,4,ib,5),ib=1,17)  & ! Band 4, RH = 84%
        /0.00003, 0.00007, 0.00017, 0.00044, 0.00105, 0.00265  &
        ,0.00670, 0.01688, 0.04262, 0.11017, 0.30307, 0.60991  &
        ,0.77778, 0.85696, 0.89753, 0.91944, 0.92917/

        DATA(gasym(2,4,ib,6),ib=1,17)  & ! Band 4, RH = 85%
        /0.00003, 0.00007, 0.00018, 0.00045, 0.00108, 0.00274  &
        ,0.00693, 0.01745, 0.04409, 0.11418, 0.31420, 0.61781  &
        ,0.78157, 0.85860, 0.89837, 0.91964, 0.92910/

        DATA(gasym(2,4,ib,7),ib=1,17)  & ! Band 4, RH = 86%
        /0.00003, 0.00007, 0.00018, 0.00047, 0.00112, 0.00284  &
        ,0.00718, 0.01809, 0.04572, 0.11868, 0.32645, 0.62609  &
        ,0.78551, 0.86032, 0.89926, 0.91985, 0.92903/

        DATA(gasym(2,4,ib,8),ib=1,17)  & ! Band 4, RH = 87%
        /0.00003, 0.00007, 0.00019, 0.00049, 0.00117, 0.00296  &
        ,0.00746, 0.01881, 0.04756, 0.12376, 0.34003, 0.63484  &
        ,0.78964, 0.86216, 0.90021, 0.92007, 0.92896/

        DATA(gasym(2,4,ib,9),ib=1,17)  & ! Band 4, RH = 88%
        /0.00003, 0.00008, 0.00020, 0.00051, 0.00122, 0.00308  &
        ,0.00779, 0.01963, 0.04967, 0.12960, 0.35519, 0.64417  &
        ,0.79402, 0.86412, 0.90122, 0.92032, 0.92891/

        DATA(gasym(2,4,ib,10),ib=1,17)  & ! Band 4, RH = 89%
        /0.00003, 0.00008, 0.00021, 0.00053, 0.00128, 0.00323  &
        ,0.00816, 0.02057, 0.05210, 0.13638, 0.37221, 0.65418  &
        ,0.79868, 0.86623, 0.90230, 0.92058, 0.92886/

        DATA(gasym(2,4,ib,11),ib=1,17)  & ! Band 4, RH = 90%
        /0.00003, 0.00008, 0.00022, 0.00056, 0.00134, 0.00340  &
        ,0.00859, 0.02168, 0.05495, 0.14440, 0.39140, 0.66501  &
        ,0.80369, 0.86852, 0.90347, 0.92087, 0.92882/

        DATA(gasym(2,4,ib,12),ib=1,17)  & ! Band 4, RH = 91%
        /0.00003, 0.00009, 0.00023, 0.00059, 0.00142, 0.00360  &
        ,0.00911, 0.02299, 0.05834, 0.15402, 0.41310, 0.67679  &
        ,0.80911, 0.87101, 0.90473, 0.92119, 0.92880/

        DATA(gasym(2,4,ib,13),ib=1,17)  & ! Band 4, RH = 92%
        /0.00004, 0.00009, 0.00025, 0.00063, 0.00152, 0.00385  &
        ,0.00974, 0.02458, 0.06247, 0.16580, 0.43768, 0.68968  &
        ,0.81498, 0.87378, 0.90609, 0.92155, 0.92880/

        DATA(gasym(2,4,ib,14),ib=1,17)  & ! Band 4, RH = 93%
        /0.00004, 0.00010, 0.00026, 0.00068, 0.00164, 0.00415  &
        ,0.01051, 0.02655, 0.06761, 0.18061, 0.46549, 0.70386  &
        ,0.82136, 0.87688, 0.90758, 0.92196, 0.92881/

        DATA(gasym(2,4,ib,15),ib=1,17)  & ! Band 4, RH = 94%
        /0.00004, 0.00011, 0.00029, 0.00074, 0.00179, 0.00454  &
        ,0.01150, 0.02906, 0.07421, 0.19980, 0.49688, 0.71954  &
        ,0.82832, 0.88039, 0.90922, 0.92244, 0.92886/

        DATA(gasym(2,4,ib,16),ib=1,17)  & ! Band 4, RH = 95%
        /0.00004, 0.00012, 0.00032, 0.00082, 0.00198, 0.00505  &
        ,0.01281, 0.03241, 0.08309, 0.22574, 0.53221, 0.73703  &
        ,0.83608, 0.88442, 0.91102, 0.92300, 0.92895/

        DATA(gasym(2,4,ib,17),ib=1,17)  & ! Band 4, RH = 96%
        /0.00005, 0.00013, 0.00036, 0.00093, 0.00226, 0.00577  &
        ,0.01466, 0.03715, 0.09581, 0.26272, 0.57203, 0.75674  &
        ,0.84497, 0.88915, 0.91304, 0.92369, 0.92911/

        DATA(gasym(2,4,ib,18),ib=1,17)  & ! Band 4, RH = 97%
        /0.00006, 0.00016, 0.00042, 0.00110, 0.00269, 0.00688  &
        ,0.01750, 0.04450, 0.11590, 0.31919, 0.61788, 0.77927  &
        ,0.85536, 0.89489, 0.91538, 0.92457, 0.92935/

        DATA(gasym(2,4,ib,19),ib=1,17)  & ! Band 4, RH = 98%
        /0.00007, 0.00019, 0.00052, 0.00140, 0.00343, 0.00885  &
        ,0.02260, 0.05782, 0.15349, 0.41168, 0.67448, 0.80640  &
        ,0.86838, 0.90202, 0.91820, 0.92575, 0.92975/

        DATA(gasym(2,4,ib,20),ib=1,17)  & ! Band 4, RH = 99%
        /0.00009, 0.00026, 0.00075, 0.00208, 0.00522, 0.01365  &
        ,0.03525, 0.09188, 0.25351, 0.56302, 0.75191, 0.84210  &
        ,0.88706, 0.91128, 0.92206, 0.92754, 0.93046/

        DATA(gasym(2,5,ib,1),ib=1,17)  & ! Band 5, RH = 80%
        /0.00005, 0.00013, 0.00034, 0.00087, 0.00209, 0.00528  &
        ,0.01332, 0.03371, 0.08674, 0.23542, 0.56870, 0.77651  &
        ,0.87177, 0.91307, 0.93289, 0.94285, 0.94791/

        DATA(gasym(2,5,ib,2),ib=1,17)  & ! Band 5, RH = 81%
        /0.00005, 0.00014, 0.00035, 0.00090, 0.00214, 0.00542  &
        ,0.01368, 0.03462, 0.08919, 0.24252, 0.57811, 0.78044  &
        ,0.87338, 0.91371, 0.93313, 0.94287, 0.94784/

        DATA(gasym(2,5,ib,3),ib=1,17)  & ! Band 5, RH = 82%
        /0.00005, 0.00014, 0.00036, 0.00092, 0.00220, 0.00557  &
        ,0.01406, 0.03560, 0.09182, 0.25015, 0.58765, 0.78446  &
        ,0.87502, 0.91437, 0.93337, 0.94289, 0.94777/

        DATA(gasym(2,5,ib,4),ib=1,17)  & ! Band 5, RH = 83%
        /0.00006, 0.00015, 0.00037, 0.00095, 0.00227, 0.00573  &
        ,0.01447, 0.03665, 0.09467, 0.25839, 0.59738, 0.78858  &
        ,0.87669, 0.91504, 0.93362, 0.94292, 0.94770/

        DATA(gasym(2,5,ib,5),ib=1,17)  & ! Band 5, RH = 84%
        /0.00006, 0.00015, 0.00038, 0.00098, 0.00233, 0.00591  &
        ,0.01492, 0.03780, 0.09779, 0.26739, 0.60733, 0.79283  &
        ,0.87842, 0.91574, 0.93387, 0.94295, 0.94763/

        DATA(gasym(2,5,ib,6),ib=1,17)  & ! Band 5, RH = 85%
        /0.00006, 0.00015, 0.00040, 0.00101, 0.00241, 0.00610  &
        ,0.01541, 0.03907, 0.10124, 0.27729, 0.61755, 0.79725  &
        ,0.88021, 0.91646, 0.93414, 0.94298, 0.94756/

        DATA(gasym(2,5,ib,7),ib=1,17)  & ! Band 5, RH = 86%
        /0.00006, 0.00016, 0.00041, 0.00104, 0.00250, 0.00631  &
        ,0.01596, 0.04048, 0.10510, 0.28831, 0.62809, 0.80185  &
        ,0.88208, 0.91722, 0.93442, 0.94302, 0.94749/

        DATA(gasym(2,5,ib,8),ib=1,17)  & ! Band 5, RH = 87%
        /0.00006, 0.00017, 0.00042, 0.00108, 0.00259, 0.00656  &
        ,0.01658, 0.04208, 0.10947, 0.30069, 0.63900, 0.80669  &
        ,0.88403, 0.91803, 0.93472, 0.94307, 0.94742/

        DATA(gasym(2,5,ib,9),ib=1,17)  & ! Band 5, RH = 88%
        /0.00007, 0.00017, 0.00044, 0.00113, 0.00270, 0.00683  &
        ,0.01728, 0.04390, 0.11449, 0.31472, 0.65032, 0.81177  &
        ,0.88609, 0.91888, 0.93504, 0.94312, 0.94735/

        DATA(gasym(2,5,ib,10),ib=1,17)  & ! Band 5, RH = 89%
        /0.00007, 0.00018, 0.00046, 0.00118, 0.00282, 0.00715  &
        ,0.01809, 0.04600, 0.12031, 0.33079, 0.66211, 0.81714  &
        ,0.88827, 0.91979, 0.93538, 0.94319, 0.94729/

        DATA(gasym(2,5,ib,11),ib=1,17)  & ! Band 5, RH = 90%
        /0.00007, 0.00019, 0.00048, 0.00124, 0.00297, 0.00752  &
        ,0.01904, 0.04847, 0.12718, 0.34936, 0.67443, 0.82283  &
        ,0.89059, 0.92077, 0.93575, 0.94326, 0.94722/

        DATA(gasym(2,5,ib,12),ib=1,17)  & ! Band 5, RH = 91%
        /0.00007, 0.00020, 0.00051, 0.00131, 0.00314, 0.00796  &
        ,0.02017, 0.05141, 0.13542, 0.37103, 0.68737, 0.82888  &
        ,0.89306, 0.92184, 0.93616, 0.94336, 0.94716/

        DATA(gasym(2,5,ib,13),ib=1,17)  & ! Band 5, RH = 92%
        /0.00008, 0.00021, 0.00054, 0.00139, 0.00335, 0.00849  &
        ,0.02153, 0.05499, 0.14552, 0.39654, 0.70108, 0.83536  &
        ,0.89571, 0.92301, 0.93661, 0.94347, 0.94711/

        DATA(gasym(2,5,ib,14),ib=1,17)  & ! Band 5, RH = 93%
        /0.00008, 0.00022, 0.00058, 0.00150, 0.00360, 0.00915  &
        ,0.02323, 0.05944, 0.15818, 0.42683, 0.71583, 0.84233  &
        ,0.89861, 0.92429, 0.93711, 0.94360, 0.94706/

        DATA(gasym(2,5,ib,15),ib=1,17)  & ! Band 5, RH = 94%
        /0.00009, 0.00024, 0.00063, 0.00163, 0.00393, 0.01000  &
        ,0.02540, 0.06516, 0.17461, 0.46304, 0.73202, 0.84991  &
        ,0.90180, 0.92572, 0.93768, 0.94377, 0.94702/

        DATA(gasym(2,5,ib,16),ib=1,17)  & ! Band 5, RH = 95%
        /0.00010, 0.00026, 0.00070, 0.00180, 0.00436, 0.01111  &
        ,0.02828, 0.07285, 0.19685, 0.50646, 0.75031, 0.85825  &
        ,0.90536, 0.92733, 0.93833, 0.94398, 0.94701/

        DATA(gasym(2,5,ib,17),ib=1,17)  & ! Band 5, RH = 96%
        /0.00011, 0.00030, 0.00078, 0.00205, 0.00496, 0.01268  &
        ,0.03236, 0.08384, 0.22883, 0.55822, 0.77157, 0.86758  &
        ,0.90937, 0.92918, 0.93911, 0.94425, 0.94701/

        DATA(gasym(2,5,ib,18),ib=1,17)  & ! Band 5, RH = 97%
        /0.00012, 0.00034, 0.00092, 0.00241, 0.00588, 0.01510  &
        ,0.03868, 0.10114, 0.27893, 0.61835, 0.79685, 0.87836  &
        ,0.91411, 0.93139, 0.94007, 0.94462, 0.94706/

        DATA(gasym(2,5,ib,19),ib=1,17)  & ! Band 5, RH = 98%
        /0.00015, 0.00042, 0.00114, 0.00305, 0.00751, 0.01941  &
        ,0.05012, 0.13328, 0.36749, 0.68528, 0.82712, 0.89117  &
        ,0.91994, 0.93414, 0.94132, 0.94514, 0.94719/

        DATA(gasym(2,5,ib,20),ib=1,17)  & ! Band 5, RH = 99%
        /0.00019, 0.00057, 0.00165, 0.00454, 0.01142, 0.02999  &
        ,0.07911, 0.21829, 0.54490, 0.76654, 0.86496, 0.90770  &
        ,0.92785, 0.93794, 0.94316, 0.94596, 0.94746/

        DATA(gasym(2,6,ib,1),ib=1,17)  & ! Band 6, RH = 80%
        /0.00015, 0.00039, 0.00099, 0.00253, 0.00603, 0.01522  &
        ,0.03838, 0.09715, 0.25615, 0.63324, 0.78852, 0.87209  &
        ,0.88964, 0.84120, 0.93843, 0.96156, 0.97038/

        DATA(gasym(2,6,ib,2),ib=1,17)  & ! Band 6, RH = 81%
        /0.00015, 0.00040, 0.00102, 0.00259, 0.00619, 0.01564  &
        ,0.03942, 0.09985, 0.26384, 0.64208, 0.79358, 0.87364  &
        ,0.89017, 0.84265, 0.93975, 0.96205, 0.97071/

        DATA(gasym(2,6,ib,3),ib=1,17)  & ! Band 6, RH = 82%
        /0.00016, 0.00041, 0.00105, 0.00267, 0.00636, 0.01608  &
        ,0.04053, 0.10274, 0.27211, 0.65060, 0.79844, 0.87498  &
        ,0.89055, 0.84462, 0.94109, 0.96251, 0.97103/

        DATA(gasym(2,6,ib,4),ib=1,17)  & ! Band 6, RH = 83%
        /0.00016, 0.00042, 0.00108, 0.00274, 0.00655, 0.01655  &
        ,0.04174, 0.10586, 0.28109, 0.65875, 0.80305, 0.87617  &
        ,0.89083, 0.84727, 0.94241, 0.96298, 0.97134/

        DATA(gasym(2,6,ib,5),ib=1,17)  & ! Band 6, RH = 84%
        /0.00017, 0.00043, 0.00111, 0.00283, 0.00675, 0.01707  &
        ,0.04305, 0.10928, 0.29091, 0.66651, 0.80743, 0.87730  &
        ,0.89111, 0.85019, 0.94369, 0.96344, 0.97163/

        DATA(gasym(2,6,ib,6),ib=1,17)  & ! Band 6, RH = 85%
        /0.00017, 0.00045, 0.00115, 0.00292, 0.00698, 0.01764  &
        ,0.04449, 0.11304, 0.30177, 0.67382, 0.81156, 0.87854  &
        ,0.89151, 0.85308, 0.94491, 0.96392, 0.97192/

        DATA(gasym(2,6,ib,7),ib=1,17)  & ! Band 6, RH = 86%
        /0.00018, 0.00046, 0.00119, 0.00302, 0.00723, 0.01827  &
        ,0.04609, 0.11724, 0.31391, 0.68062, 0.81549, 0.88011  &
        ,0.89199, 0.85652, 0.94615, 0.96443, 0.97220/

        DATA(gasym(2,6,ib,8),ib=1,17)  & ! Band 6, RH = 87%
        /0.00018, 0.00048, 0.00123, 0.00314, 0.00750, 0.01898  &
        ,0.04790, 0.12198, 0.32762, 0.68687, 0.81928, 0.88226  &
        ,0.89225, 0.86119, 0.94743, 0.96499, 0.97248/

        DATA(gasym(2,6,ib,9),ib=1,17)  & ! Band 6, RH = 88%
        /0.00019, 0.00050, 0.00128, 0.00327, 0.00782, 0.01978  &
        ,0.04996, 0.12739, 0.34328, 0.69254, 0.82306, 0.88515  &
        ,0.89195, 0.86634, 0.94871, 0.96562, 0.97276/

        DATA(gasym(2,6,ib,10),ib=1,17)  & ! Band 6, RH = 89%
        /0.00020, 0.00052, 0.00134, 0.00342, 0.00818, 0.02071  &
        ,0.05233, 0.13366, 0.36137, 0.69770, 0.82703, 0.88868  &
        ,0.89111, 0.87203, 0.95005, 0.96634, 0.97304/

        DATA(gasym(2,6,ib,11),ib=1,17)  & ! Band 6, RH = 90%
        /0.00021, 0.00054, 0.00140, 0.00359, 0.00861, 0.02180  &
        ,0.05510, 0.14104, 0.38249, 0.70253, 0.83145, 0.89226  &
        ,0.89044, 0.87947, 0.95156, 0.96718, 0.97335/

        DATA(gasym(2,6,ib,12),ib=1,17)  & ! Band 6, RH = 91%
        /0.00022, 0.00057, 0.00148, 0.00380, 0.00911, 0.02308  &
        ,0.05840, 0.14985, 0.40742, 0.70753, 0.83661, 0.89495  &
        ,0.89006, 0.88727, 0.95321, 0.96813, 0.97368/

        DATA(gasym(2,6,ib,13),ib=1,17)  & ! Band 6, RH = 92%
        /0.00023, 0.00061, 0.00158, 0.00405, 0.00972, 0.02464  &
        ,0.06240, 0.16061, 0.43714, 0.71363, 0.84271, 0.89630  &
        ,0.88865, 0.89695, 0.95518, 0.96918, 0.97402/

        DATA(gasym(2,6,ib,14),ib=1,17)  & ! Band 6, RH = 93%
        /0.00024, 0.00065, 0.00169, 0.00435, 0.01046, 0.02657  &
        ,0.06735, 0.17405, 0.47285, 0.72241, 0.84971, 0.89759  &
        ,0.88628, 0.90700, 0.95756, 0.97026, 0.97436/

        DATA(gasym(2,6,ib,15),ib=1,17)  & ! Band 6, RH = 94%
        /0.00026, 0.00070, 0.00184, 0.00474, 0.01142, 0.02902  &
        ,0.07368, 0.19144, 0.51580, 0.73629, 0.85764, 0.90077  &
        ,0.88342, 0.91813, 0.96038, 0.97128, 0.97471/

        DATA(gasym(2,6,ib,16),ib=1,17)  & ! Band 6, RH = 95%
        /0.00029, 0.00077, 0.00202, 0.00525, 0.01268, 0.03227  &
        ,0.08213, 0.21495, 0.56668, 0.75791, 0.86722, 0.90382  &
        ,0.87974, 0.92950, 0.96351, 0.97218, 0.97507/

        DATA(gasym(2,6,ib,17),ib=1,17)  & ! Band 6, RH = 96%
        /0.00032, 0.00086, 0.00229, 0.00596, 0.01444, 0.03685  &
        ,0.09410, 0.24885, 0.62348, 0.78727, 0.87840, 0.90516  &
        ,0.87546, 0.94040, 0.96647, 0.97308, 0.97546/

        DATA(gasym(2,6,ib,18),ib=1,17)  & ! Band 6, RH = 97%
        /0.00036, 0.00099, 0.00268, 0.00703, 0.01713, 0.04389  &
        ,0.11272, 0.30262, 0.67653, 0.81631, 0.88584, 0.90556  &
        ,0.87435, 0.94957, 0.96848, 0.97415, 0.97586/

        DATA(gasym(2,6,ib,19),ib=1,17)  & ! Band 6, RH = 98%
        /0.00043, 0.00121, 0.00334, 0.00890, 0.02188, 0.05646  &
        ,0.14667, 0.40115, 0.71016, 0.83952, 0.89962, 0.90192  &
        ,0.89116, 0.95652, 0.97015, 0.97506, 0.97630/

        DATA(gasym(2,6,ib,20),ib=1,17)  & ! Band 6, RH = 99%
        /0.00055, 0.00167, 0.00480, 0.01323, 0.03325, 0.08738  &
        ,0.23464, 0.60652, 0.78084, 0.87799, 0.90836, 0.88433  &
        ,0.93891, 0.96699, 0.97369, 0.97604, 0.97682/

        DATA(gasym(2,7,ib,1),ib=1,17)  & ! Band 7, RH = 80%
        /0.00012, 0.00032, 0.00081, 0.00206, 0.00492, 0.01242  &
        ,0.03134, 0.07927, 0.20624, 0.54493, 0.75160, 0.86569  &
        ,0.90725, 0.89738, 0.92787, 0.96496, 0.97336/

        DATA(gasym(2,7,ib,2),ib=1,17)  & ! Band 7, RH = 81%
        /0.00013, 0.00033, 0.00083, 0.00212, 0.00505, 0.01276  &
        ,0.03219, 0.08146, 0.21230, 0.55717, 0.75707, 0.86831  &
        ,0.90831, 0.89791, 0.93035, 0.96570, 0.97367/

        DATA(gasym(2,7,ib,3),ib=1,17)  & ! Band 7, RH = 82%
        /0.00013, 0.00033, 0.00085, 0.00217, 0.00519, 0.01312  &
        ,0.03310, 0.08381, 0.21882, 0.56967, 0.76288, 0.87101  &
        ,0.90938, 0.89839, 0.93276, 0.96644, 0.97397/

        DATA(gasym(2,7,ib,4),ib=1,17)  & ! Band 7, RH = 83%
        /0.00013, 0.00034, 0.00088, 0.00224, 0.00534, 0.01351  &
        ,0.03409, 0.08634, 0.22589, 0.58246, 0.76904, 0.87377  &
        ,0.91049, 0.89892, 0.93513, 0.96716, 0.97427/

        DATA(gasym(2,7,ib,5),ib=1,17)  & ! Band 7, RH = 84%
        /0.00014, 0.00035, 0.00091, 0.00231, 0.00551, 0.01393  &
        ,0.03516, 0.08911, 0.23363, 0.59553, 0.77554, 0.87657  &
        ,0.91164, 0.89960, 0.93757, 0.96787, 0.97458/

        DATA(gasym(2,7,ib,6),ib=1,17)  & ! Band 7, RH = 85%
        /0.00014, 0.00036, 0.00094, 0.00238, 0.00569, 0.01439  &
        ,0.03633, 0.09216, 0.24219, 0.60888, 0.78236, 0.87937  &
        ,0.91278, 0.90033, 0.94001, 0.96859, 0.97489/

        DATA(gasym(2,7,ib,7),ib=1,17)  & ! Band 7, RH = 86%
        /0.00015, 0.00038, 0.00097, 0.00247, 0.00589, 0.01491  &
        ,0.03765, 0.09556, 0.25177, 0.62247, 0.78947, 0.88212  &
        ,0.91387, 0.90095, 0.94245, 0.96930, 0.97520/

        DATA(gasym(2,7,ib,8),ib=1,17)  & ! Band 7, RH = 87%
        /0.00015, 0.00039, 0.00100, 0.00256, 0.00612, 0.01549  &
        ,0.03912, 0.09939, 0.26261, 0.63619, 0.79678, 0.88476  &
        ,0.91489, 0.90152, 0.94495, 0.97000, 0.97552/

        DATA(gasym(2,7,ib,9),ib=1,17)  & ! Band 7, RH = 88%
        /0.00016, 0.00041, 0.00104, 0.00267, 0.00638, 0.01615  &
        ,0.04080, 0.10377, 0.27503, 0.64990, 0.80419, 0.88729  &
        ,0.91592, 0.90236, 0.94742, 0.97069, 0.97584/

        DATA(gasym(2,7,ib,10),ib=1,17)  & ! Band 7, RH = 89%
        /0.00016, 0.00042, 0.00109, 0.00279, 0.00667, 0.01691  &
        ,0.04274, 0.10884, 0.28944, 0.66337, 0.81160, 0.88978  &
        ,0.91707, 0.90341, 0.94999, 0.97138, 0.97618/

        DATA(gasym(2,7,ib,11),ib=1,17)  & ! Band 7, RH = 90%
        /0.00017, 0.00044, 0.00114, 0.00293, 0.00702, 0.01779  &
        ,0.04501, 0.11478, 0.30638, 0.67631, 0.81890, 0.89248  &
        ,0.91832, 0.90460, 0.95251, 0.97206, 0.97651/

        DATA(gasym(2,7,ib,12),ib=1,17)  & ! Band 7, RH = 91%
        /0.00018, 0.00047, 0.00121, 0.00310, 0.00743, 0.01884  &
        ,0.04770, 0.12188, 0.32660, 0.68839, 0.82604, 0.89574  &
        ,0.91935, 0.90627, 0.95509, 0.97273, 0.97685/

        DATA(gasym(2,7,ib,13),ib=1,17)  & ! Band 7, RH = 92%
        /0.00019, 0.00049, 0.00128, 0.00330, 0.00793, 0.02012  &
        ,0.05096, 0.13052, 0.35111, 0.69940, 0.83315, 0.89978  &
        ,0.92009, 0.90837, 0.95765, 0.97340, 0.97719/

        DATA(gasym(2,7,ib,14),ib=1,17)  & ! Band 7, RH = 93%
        /0.00020, 0.00053, 0.00138, 0.00355, 0.00854, 0.02169  &
        ,0.05500, 0.14129, 0.38140, 0.70944, 0.84056, 0.90424  &
        ,0.92108, 0.91157, 0.96021, 0.97409, 0.97754/

        DATA(gasym(2,7,ib,15),ib=1,17)  & ! Band 7, RH = 94%
        /0.00021, 0.00057, 0.00150, 0.00387, 0.00931, 0.02369  &
        ,0.06017, 0.15518, 0.41954, 0.71947, 0.84878, 0.90832  &
        ,0.92171, 0.91604, 0.96280, 0.97481, 0.97790/

        DATA(gasym(2,7,ib,16),ib=1,17)  & ! Band 7, RH = 95%
        /0.00023, 0.00063, 0.00165, 0.00428, 0.01034, 0.02635  &
        ,0.06706, 0.17388, 0.46850, 0.73214, 0.85832, 0.91235  &
        ,0.92184, 0.92243, 0.96547, 0.97559, 0.97827/

        DATA(gasym(2,7,ib,17),ib=1,17)  & ! Band 7, RH = 96%
        /0.00026, 0.00070, 0.00186, 0.00486, 0.01178, 0.03010  &
        ,0.07680, 0.20073, 0.53176, 0.75268, 0.87005, 0.91715  &
        ,0.92147, 0.93142, 0.96834, 0.97643, 0.97867/

        DATA(gasym(2,7,ib,18),ib=1,17)  & ! Band 7, RH = 97%
        /0.00029, 0.00081, 0.00218, 0.00574, 0.01398, 0.03585  &
        ,0.09192, 0.24320, 0.61010, 0.78691, 0.88521, 0.92196  &
        ,0.92012, 0.94321, 0.97149, 0.97729, 0.97910/

        DATA(gasym(2,7,ib,19),ib=1,17)  & ! Band 7, RH = 98%
        /0.00035, 0.00099, 0.00272, 0.00726, 0.01786, 0.04611  &
        ,0.11935, 0.32192, 0.68801, 0.82792, 0.89941, 0.92626  &
        ,0.91876, 0.95653, 0.97455, 0.97825, 0.97957/

        DATA(gasym(2,7,ib,20),ib=1,17)  & ! Band 7, RH = 99%
        /0.00045, 0.00136, 0.00392, 0.01080, 0.02716, 0.07134  &
        ,0.18956, 0.51157, 0.74833, 0.86898, 0.91861, 0.92646  &
        ,0.93225, 0.96876, 0.97709, 0.97928, 0.98014/

        DATA(gasym(2,8,ib,1),ib=1,17)  & ! Band 8, RH = 80%
        /0.00042, 0.00108, 0.00277, 0.00703, 0.01674, 0.04215  &
        ,0.10638, 0.27864, 0.61643, 0.76936, 0.84399, 0.83886  &
        ,0.79308, 0.88396, 0.92410, 0.94751, 0.96205/

        DATA(gasym(2,8,ib,2),ib=1,17)  & ! Band 8, RH = 81%
        /0.00043, 0.00111, 0.00284, 0.00722, 0.01719, 0.04329  &
        ,0.10931, 0.28684, 0.62246, 0.77267, 0.84581, 0.83798  &
        ,0.79638, 0.88633, 0.92539, 0.94838, 0.96257/

        DATA(gasym(2,8,ib,3),ib=1,17)  & ! Band 8, RH = 82%
        /0.00044, 0.00114, 0.00292, 0.00742, 0.01767, 0.04451  &
        ,0.11246, 0.29563, 0.62835, 0.77607, 0.84764, 0.83719  &
        ,0.80023, 0.88856, 0.92659, 0.94927, 0.96307/

        DATA(gasym(2,8,ib,4),ib=1,17)  & ! Band 8, RH = 83%
        /0.00045, 0.00117, 0.00300, 0.00764, 0.01819, 0.04583  &
        ,0.11586, 0.30513, 0.63414, 0.77959, 0.84941, 0.83661  &
        ,0.80440, 0.89104, 0.92778, 0.95020, 0.96356/

        DATA(gasym(2,8,ib,5),ib=1,17)  & ! Band 8, RH = 84%
        /0.00047, 0.00121, 0.00309, 0.00787, 0.01876, 0.04726  &
        ,0.11956, 0.31548, 0.63991, 0.78328, 0.85105, 0.83614  &
        ,0.80821, 0.89378, 0.92893, 0.95117, 0.96405/

        DATA(gasym(2,8,ib,6),ib=1,17)  & ! Band 8, RH = 85%
        /0.00048, 0.00125, 0.00319, 0.00813, 0.01938, 0.04884  &
        ,0.12365, 0.32687, 0.64572, 0.78721, 0.85251, 0.83539  &
        ,0.81239, 0.89631, 0.93011, 0.95221, 0.96454/

        DATA(gasym(2,8,ib,7),ib=1,17)  & ! Band 8, RH = 86%
        /0.00050, 0.00129, 0.00330, 0.00841, 0.02006, 0.05059  &
        ,0.12821, 0.33953, 0.65170, 0.79140, 0.85379, 0.83394  &
        ,0.81771, 0.89904, 0.93126, 0.95330, 0.96505/

        DATA(gasym(2,8,ib,8),ib=1,17)  & ! Band 8, RH = 87%
        /0.00051, 0.00133, 0.00343, 0.00873, 0.02084, 0.05255  &
        ,0.13334, 0.35372, 0.65797, 0.79588, 0.85494, 0.83181  &
        ,0.82321, 0.90209, 0.93253, 0.95442, 0.96558/

        DATA(gasym(2,8,ib,9),ib=1,17)  & ! Band 8, RH = 88%
        /0.00053, 0.00139, 0.00357, 0.00910, 0.02171, 0.05479  &
        ,0.13921, 0.36977, 0.66473, 0.80066, 0.85614, 0.82981  &
        ,0.82836, 0.90478, 0.93393, 0.95555, 0.96613/

        DATA(gasym(2,8,ib,10),ib=1,17)  & ! Band 8, RH = 89%
        /0.00055, 0.00145, 0.00372, 0.00951, 0.02272, 0.05737  &
        ,0.14600, 0.38809, 0.67226, 0.80574, 0.85762, 0.82801  &
        ,0.83483, 0.90789, 0.93556, 0.95663, 0.96668/

        DATA(gasym(2,8,ib,11),ib=1,17)  & ! Band 8, RH = 90%
        /0.00058, 0.00151, 0.00391, 0.01000, 0.02390, 0.06038  &
        ,0.15397, 0.40914, 0.68090, 0.81113, 0.85946, 0.82520  &
        ,0.84143, 0.91059, 0.93750, 0.95768, 0.96724/

        DATA(gasym(2,8,ib,12),ib=1,17)  & ! Band 8, RH = 91%
        /0.00061, 0.00159, 0.00413, 0.01057, 0.02530, 0.06395  &
        ,0.16349, 0.43345, 0.69111, 0.81690, 0.86129, 0.82161  &
        ,0.84824, 0.91334, 0.93976, 0.95874, 0.96783/

        DATA(gasym(2,8,ib,13),ib=1,17)  & ! Band 8, RH = 92%
        /0.00064, 0.00169, 0.00439, 0.01126, 0.02698, 0.06827  &
        ,0.17509, 0.46158, 0.70337, 0.82320, 0.86247, 0.81832  &
        ,0.85560, 0.91596, 0.94228, 0.95995, 0.96841/

        DATA(gasym(2,8,ib,14),ib=1,17)  & ! Band 8, RH = 93%
        /0.00068, 0.00181, 0.00471, 0.01211, 0.02905, 0.07363  &
        ,0.18956, 0.49400, 0.71804, 0.83009, 0.86306, 0.81357  &
        ,0.86282, 0.91842, 0.94481, 0.96132, 0.96900/

        DATA(gasym(2,8,ib,15),ib=1,17)  & ! Band 8, RH = 94%
        /0.00073, 0.00195, 0.00511, 0.01319, 0.03169, 0.08045  &
        ,0.20822, 0.53074, 0.73500, 0.83727, 0.86354, 0.81046  &
        ,0.87059, 0.92103, 0.94724, 0.96273, 0.96963/

        DATA(gasym(2,8,ib,16),ib=1,17)  & ! Band 8, RH = 95%
        /0.00079, 0.00214, 0.00564, 0.01461, 0.03518, 0.08953  &
        ,0.23336, 0.57082, 0.75331, 0.84405, 0.86235, 0.80778  &
        ,0.87889, 0.92489, 0.94989, 0.96424, 0.97027/

        DATA(gasym(2,8,ib,17),ib=1,17)  & ! Band 8, RH = 96%
        /0.00088, 0.00240, 0.00637, 0.01658, 0.04007, 0.10235  &
        ,0.26934, 0.61137, 0.77194, 0.85153, 0.86012, 0.81003  &
        ,0.88820, 0.93077, 0.95327, 0.96583, 0.97092/

        DATA(gasym(2,8,ib,18),ib=1,17)  & ! Band 8, RH = 97%
        /0.00100, 0.00277, 0.00745, 0.01956, 0.04752, 0.12218  &
        ,0.32554, 0.64888, 0.79263, 0.86063, 0.85335, 0.82235  &
        ,0.90136, 0.93724, 0.95670, 0.96759, 0.97157/

        DATA(gasym(2,8,ib,19),ib=1,17)  & ! Band 8, RH = 98%
        /0.00119, 0.00339, 0.00930, 0.02474, 0.06069, 0.15804  &
        ,0.42407, 0.69110, 0.81948, 0.86729, 0.83686, 0.85109  &
        ,0.91734, 0.94330, 0.96135, 0.96946, 0.97224/

        DATA(gasym(2,8,ib,20),ib=1,17)  & ! Band 8, RH = 99%
        /0.00154, 0.00466, 0.01338, 0.03675, 0.09234, 0.24911  &
        ,0.59662, 0.76793, 0.85173, 0.86554, 0.81500, 0.88824  &
        ,0.93140, 0.95415, 0.96660, 0.97150, 0.97288/

! gasym for Mineral Dust:
        DATA(gasym(3,1,ib,1),ib=1,17)  & ! Band 1, All RH's
        /0.00107, 0.00268, 0.00676, 0.01699, 0.04264, 0.10747  &
        ,0.28301, 0.61817, 0.76309, 0.83417, 0.83819, 0.85421  &
        ,0.88159, 0.89181, 0.92795, 0.95285, 0.96718/

        DATA(gasym(3,2,ib,1),ib=1,17)  & ! Band 2, All RH's
        /0.00490, 0.01231, 0.03085, 0.07719, 0.19842, 0.50460  &
        ,0.66831, 0.73464, 0.68896, 0.72559, 0.77863, 0.83901  &
        ,0.87691, 0.91052, 0.93413, 0.94697, 0.95155/

        DATA(gasym(3,3,ib,1),ib=1,17)  & ! Band 3, All RH's
        /0.02390, 0.05985, 0.15098, 0.36031, 0.58493, 0.69319  &
        ,0.71724, 0.63499, 0.74390, 0.80146, 0.84751, 0.88670  &
        ,0.91730, 0.93703, 0.94596, 0.94845, 0.94879/

        DATA(gasym(3,4,ib,1),ib=1,17)  & ! Band 4, All RH's
        /0.00001, 0.00002, 0.00005, 0.00011, 0.00029, 0.00073  &
        ,0.00183, 0.00461, 0.01159, 0.02909, 0.07352, 0.19723  &
        ,0.49168, 0.66569, 0.75801, 0.84784, 0.87389/

        DATA(gasym(3,5,ib,1),ib=1,17)  & ! Band 5, All RH's
        /0.00002, 0.00004, 0.00010, 0.00026, 0.00065, 0.00165  &
        ,0.00415, 0.01043, 0.02617, 0.06562, 0.16892, 0.46511  &
        ,0.67484, 0.78564, 0.77736, 0.89280, 0.91936/

        DATA(gasym(3,6,ib,1),ib=1,17)  & ! Band 6, All RH's
        /0.00004, 0.00011, 0.00026, 0.00067, 0.00169, 0.00426  &
        ,0.01071, 0.02690, 0.06781, 0.17736, 0.50452, 0.71737  &
        ,0.82994, 0.88285, 0.92003, 0.94431, 0.95420/

        DATA(gasym(3,7,ib,1),ib=1,17)  & ! Band 7, All RH's
        /0.00004, 0.00010, 0.00025, 0.00064, 0.00160, 0.00403  &
        ,0.01014, 0.02545, 0.06394, 0.16706, 0.49135, 0.67547  &
        ,0.75696, 0.82304, 0.89691, 0.91021, 0.91682/

        DATA(gasym(3,8,ib,1),ib=1,17)  & ! Band 8, All RH's
        /0.00011, 0.00027, 0.00069, 0.00173, 0.00436, 0.01087  &
        ,0.02761, 0.06953, 0.17772, 0.45807, 0.71820, 0.84603  &
        ,0.90347, 0.91302, 0.90597, 0.94409, 0.96733/

! Returning value for gasym:
        value = gasym(aerotype,radband,bin,rh)

return
END SUBROUTINE aerogasym
