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

Subroutine stars (ustar,tstar,rstar,ths,rvs,thetacan,can_rvap,zts,patch_rough  &
   ,vels_pat,dtllohcc,dens,dtll)

use rconstants

implicit none

real :: ustar,tstar,rstar,ths,rvs,thetacan,can_rvap,zts,patch_rough  &
       ,vels_pat,dtllohcc,dens,dtll
real :: b,csm,csh,d,a2,c1,ri,fm,fh,c2,cm,ch,c3
real :: ustaro,d_vel,d_veln,vel_new,delz,ubmin,tstaro,d_theta,d_thetan,theta_new
real, parameter :: ustmin = .1

! Routine to compute Louis (1981) surface layer parameterization.

b = 5.
csm = 7.5
csh = 5.
d = 5.

! a2 is the drag coefficient in neutral conditions, here same for h/m.
! ri is the bulk richardson numer, eq. 3.45 in Garratt.
! fm,fh are stability functions in various forms.
! However, a2*fm, a2*fh appear to be what Louis(1981) calls drag
! or transfer coefficients.

a2 = (vonk / log(zts / patch_rough)) ** 2
c1 = a2 * vels_pat
ri = 9.8 * zts * (ths - thetacan)  &
   / (.5 * (ths + thetacan) * vels_pat * vels_pat)

if (ths - thetacan > 0.) then   ! STABLE CASE
   fm = 1. / (1. + (2 * b * ri / sqrt(1 + d * ri)))
   fh = 1. / (1. + (3 * b * ri * sqrt(1 + d * ri)))
else                            ! UNSTABLE CASE
   c2 = b * a2 * sqrt(zts / patch_rough * (abs(ri)))
   cm = csm * c2
   ch = csh * c2
   fm = (1. - 2 * b * ri / (1. + 2 * cm))
   fh = (1. - 3 * b * ri / (1. + 3 * ch))
endif

ustar = max(ustmin,sqrt(c1 * vels_pat * fm))
c3 = c1 * fh / ustar
rstar = c3 * (rvs - can_rvap)
tstar = c3 * (ths - thetacan)

! Limit ustar so that the flux cannot take more than 1/2 velocity in a timestep

ustaro=ustar
delz = 2.*zts
d_vel = - ustar * ustar * dtll / delz
vel_new = vels_pat + d_vel
if(vel_new < .5 * vels_pat) then
   d_veln = .5 * vels_pat
   ustar = sqrt(d_veln * delz / dtll)
endif

! Limit tstar and rstar so that the direction of the gradients cannot change
!   sign due to the fluxes - also dampen solution for stability
tstaro = tstar
delz = 2.*zts
d_theta = dtllohcc * dens * cp * ustar * tstar
theta_new = thetacan + d_theta
if(thetacan < ths) then         ! STABLE CASE
   if(theta_new > ths) then
      d_thetan = ths - thetacan
      tstar = .75 * d_thetan / (dtllohcc * dens * cp * ustar)
   endif
else                            ! UNSTABLE CASE
   if(theta_new < ths) then
      d_thetan = ths - thetacan
      tstar = .75 * d_thetan / (dtllohcc * dens * cp * ustar)
   endif
endif

return
END SUBROUTINE stars

!##############################################################################
Subroutine sfclmcv (ustar,tstar,rstar,vels_pat,ups,vps,gzotheta,patch_area  &
     ,sflux_u,sflux_v,sflux_w,sflux_t,sflux_r)
!  +---------------------------------------------------------------+
!  \  This routine computes the turbulent fluxes of momentum,      \
!  \  heat and moisture from the surface layer using the           \
!  \  Manton-Cotton algebraic surface layer equations.             \
!  +---------------------------------------------------------------+

use rconstants

implicit none

real :: ustar,tstar,rstar,vels_pat,ups,vps,gzotheta,patch_area  &
       ,sflux_u,sflux_v,sflux_w,sflux_t,sflux_r,zoverl

real :: wtol,cosine1,sine1,vtscr,cx,psin

data wtol/1e-20/

cosine1 = ups / vels_pat
sine1   = vps / vels_pat

vtscr = ustar * patch_area

sflux_u = sflux_u - ustar * vtscr * cosine1 
sflux_v = sflux_v - ustar * vtscr * sine1  
sflux_t = sflux_t - tstar * vtscr
sflux_r = sflux_r - rstar * vtscr

zoverl = gzotheta * vonk * tstar / (ustar * ustar)

if (zoverl < 0.)then
   cx = zoverl * sqrt(sqrt(1. - 15. * zoverl))
else
   cx = zoverl / (1.0 + 4.7 * zoverl)
endif

psin = sqrt((1.-2.86 * cx) / (1. + cx * (-5.39 + cx * 6.998 )))
sflux_w = sflux_w + (0.27 * max(6.25 * (1. - cx) * psin,wtol)  &
   - 1.18 * cx * psin) * ustar * vtscr

return
END SUBROUTINE sfclmcv
