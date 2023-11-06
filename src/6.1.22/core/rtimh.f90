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

Subroutine timestep ()

use mem_basic
use node_mod
use mem_radiate
use mem_cuparm
use mem_varinit
use mem_turb
use mem_oda,   only:if_oda
use micphys,   only:level,icheckmic,itempnudge,iforceccn,fccnstart,iccnlev
use mem_grid

implicit none

integer :: callmassflux,massfluxfreq

if(iforceccn.eq.9) then
   if(time.lt.fccnstart) then
      iccnlev=0
   else
      iccnlev=2
   endif 
endif

 CALL acctimes ('INIT')

!        +-------------------------------------------------------------+
!        |   Timestep driver for the hybrid non-hydrostatic time-split |
!        |      model.                                                 |
!        +-------------------------------------------------------------+

!  Zero out all tendency arrays.   
!--------------------------------
 CALL tend0 ()          
 CALL acctimes ('TEND0')
! if(MOD(TIME,300.0).eq.0   ) then
!    print *, "Running temperature nudging" 
!    CALL temp_adj(mzp,mxp,myp,&
!       basic_g(ngrid)%thp(1,1,1),&
!       basic_g(ngrid)%theta(1,1,1))
! endif

! CALL temp_adj(mzp,mxp,myp,&
!    basic_g(ngrid)%thp(1,1,1),&
!    basic_g(ngrid)%theta(1,1,1))


! If temp nudging is on, and not aerosol forcing (ifoceccn=0) OR
! if temp nudging is on, and aerosol forcing but not yet FCCNSTART
if ((itempnudge.ne.0) .and. ((time < fccnstart) .or. (iforceccn.eq.0))) then
   ! print *, "Nudging temperature to initial sounding"
   
   if (time.eq.1) then
      call initial_temp_profile(mzp, mxp, myp, basic_g(ngrid)%thp(1,1,1))
   endif
   
   CALL temp_adj(mzp,mxp,myp, basic_g(ngrid)%thp(1,1,1))

endif

!  Thermodynamic diagnosis   
!--------------------------------
! if (level  /=  3) then
   CALL thermo (mzp,mxp,myp,ia,iz,ja,jz) 
! endif
 CALL acctimes ('THERMO')

!  Aerosol Initial & Sources of dust and salt
!--------------------------------------------
 CALL aerosols ()
 CALL acctimes ('AEROSOLS')

!  Radiation parameterization
!--------------------------------
 CALL radiate (mzp,mxp,myp,ia,iz,ja,jz)   
 CALL acctimes ('RADIATE')
   
!  Surface layer, soil, veggie, urban models
!--------------------------------------------
 CALL sfc_driver (mzp,mxp,myp,ia,iz,ja,jz) 
 CALL acctimes ('SFC_DRIVER')

!  Coriolis terms
!  ----------------------------------------
 CALL corlos (mzp,mxp,myp,i0,j0,ia,iz,ja,jz,izu,jzv) 
 CALL acctimes ('CORLOS')

!  Velocity advection
!----------------------------------------
 CALL advectc ('V',mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,mynum)
 CALL acctimes ('ADVECTv')

!  Cumulus parameterization
!----------------------------------------
 IF(NNQPARM(ngrid) == 1 ) CALL cuparm ()      
 IF(NNQPARM(ngrid) == 2 ) CALL kf_main ()      
 CALL acctimes ('CUPARM')

!  Analysis nudging and boundary condition
!------------------------------------------
 IF(NUD_TYPE > 0) CALL datassim ()  
 CALL acctimes ('DATASSIM')

!  Observation data assimilation 
!----------------------------------------
 IF(IF_ODA == 1) CALL oda_nudge ()  
 CALL acctimes ('DATASSIM')

!  Nested grid boundaries
!----------------------------------------
 if(nxtnest(ngrid) >= 1) CALL nstbdriv ()  
 CALL acctimes ('NSTBDRIV')

!  Rayleigh friction for theta
!----------------------------------------
 CALL rayft ()           
 CALL acctimes ('RAYFT')

!  Update the overlap region between parallel nodes
!---------------------------------------------------
 if(ipara == 1) then      
   CALL node_sendlbc (0)  
   CALL node_getlbc (0)  
 endif
 CALL acctimes ('UPDATElbc')

!  Sub-grid diffusion terms
!----------------------------------------
 CALL diffuse ()
 CALL acctimes ('DIFFUSE')

!  Scalar advection
!----------------------------------------
 CALL advectc ('T',mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,mynum)
 CALL acctimes ('ADVECTs')

!  Update scalars
!----------------------------------------
 CALL predtr ()          
 CALL acctimes ('PREDTR')

! Adele - added Thermodynamic diagnosis   
!--------------------------------
! if (level  /=  3) then
   CALL thermo (mzp,mxp,myp,ia,iz,ja,jz) 
! endif
 CALL acctimes ('THERMO')

!  Moisture variables positive definite
!----------------------------------------
 CALL negadj1 (mzp,mxp,myp) 
 CALL acctimes ('NEGADJ1')

!  Microphysics
!----------------------------------------
 if (level == 3) then
   !if (time<=5 .or. time>=3600) then
   CALL micro ()
   !endif
 endif
 CALL acctimes ('MICRO')

!  Check for negative micro and Nans
!----------------------------------------
 if(icheckmic == 1) then
  CALL checkmicro ('r')
 endif

!  Thermodynamic diagnosis
!----------------------------------------
! if (level /= 3) then
   CALL thermo (mzp,mxp,myp,1,mxp,1,myp) 
! endif
 CALL acctimes ('THERMO')

!  Apply scalar b.c.'s
!----------------------------------------
 CALL trsets ()          
 CALL acctimes ('TRSETS')

!  Lateral velocity boundaries - radiative
!-------------------------------------------
 CALL latbnd ()         
 CALL acctimes ('LATBND')

!  First stage Asselin filter
!----------------------------------------
 CALL hadvance (1)     
 CALL acctimes ('HADVANCE')

!  Buoyancy term for w equation
!----------------------------------------
 CALL buoyancy ()
 CALL acctimes ('BUOYANCY')

!  Acoustic small timesteps
!----------------------------------------
 CALL acoustic ()
 CALL acctimes ('ACOUSTIC')

!  Last stage of Asselin filter
!----------------------------------------
 CALL hadvance (2)      
 CALL acctimes ('HADVANCE')

!  Velocity/pressure boundary conditions
!----------------------------------------
 CALL vpsets ()          
 CALL acctimes ('VPSETS')

callmassflux=0    !flag for output BC mass flux: (0=off, 1==on)
massfluxfreq=300. !frequency of BC mass flux (seconds)
if(callmassflux==1) &
 CALL mass_flux_bc (mzp,mxp,myp &
      ,basic_g(ngrid)%up(1,1,1),basic_g(ngrid)%vp(1,1,1)  &
      ,basic_g(ngrid)%dn0(1,1,1) &
      ,basic_g(ngrid)%pp(1,1,1),basic_g(ngrid)%pi0(1,1,1) &
      ,grid_g(ngrid)%rtgu(1,1) ,grid_g(ngrid)%rtgv(1,1)    &
      ,grid_g(ngrid)%dyu(1,1)  ,grid_g(ngrid)%dyv(1,1))

return
END SUBROUTINE timestep

!##############################################################################
Subroutine acctimes (string)

use mem_all
use node_mod

implicit none

real, external :: valugp
real :: total,total2,dn,pitotal
real,save :: lasttotal
integer :: ip,jp,kp,i,j,k,patch
character(len=*) :: string
logical, external :: isnanr

kp=0
ip=20
jp=25
!  only here for debugging purposes
101    format ('DEBUG: NODE',i0,': LOC(',i0,',',i0,',',i0,') --> (' &
   ,i0,',',i0,',',i0,') ',a,': ',a10,f9.1,100e20.10)

!do ip=1,6
!do jp=1,6
!do kp=1,200
!if (isnanr(micro_g(1)%salt_film_np(kp,ip,jp) )) print*,string
!enddo
!enddo
!enddo
return
END SUBROUTINE acctimes

!##############################################################################
Subroutine output_points ()

use mem_all
use node_mod

implicit none

real :: a,b,c
integer :: k,i,j,ii,jj

do i=1,mxp
 do j=1,myp
  do k=1,mzp
    ii = i+mi0(ngrid)
    jj = j+mj0(ngrid)
    a=basic_g(ngrid)%wp(k,ii,jj)
    b=micro_g(ngrid)%rcp(k,ii,jj)
    if(k==50.and.ii==200.and.jj==1) print*,'output',ngrid,k,ii,jj,a,b
  enddo
 enddo
enddo

return
END SUBROUTINE output_points

!##############################################################################
Subroutine mass_flux_bc (m1,m2,m3,up,vp,dn0,pp,pi0,rtgu,rtgv,dyu,dxv)

use mem_grid
use mem_basic
use rconstants
use node_mod

implicit none

integer :: m1,m2,m3,i,j,k
real,dimension(m1,m2,m3) :: up,vp,dn0,pp,pi0
real,dimension(m2,m3) :: rtgu,rtgv,dyu,dxv
real :: wmass,emass,smass,nmass,prtot,tmass,ppp,area
real, save :: aintmass=0.

if(mod(time,300.).gt.0.1) return

if(ipara==1)then
 print*,'Cannot compute boundary condition mass flux in parallel.'
 stop 'See TIMESTEP and MASS_FLUX_BC in rtimh.f90'
endif

wmass=0.
emass=0.
smass=0.
nmass=0.
prtot=0.
area=0.
!************************************************************
if(jdim==1) then  !for 3D simulations
!  west/east bound
do j=2,nyp-1
   do k=2,nzp-1
      i=1
      wmass=wmass +  &
           up(k,i,j)*rtgu(i,j)/(dyu(i,j)*dzt(k))  &
           *(dn0(k,i,j)+dn0(k,i+1,j))*.5
      i=nxp-1
      emass=emass -  &
           up(k,i,j)*rtgu(i,j)/(dyu(i,j)*dzt(k))  &
           *(dn0(k,i,j)+dn0(k,i+1,j))*.5
   enddo
enddo
!  north/south bound
do i=2,nxp-1
   do k=2,nzp-1
      j=1
      smass=smass +  &
           vp(k,i,j)*rtgv(i,j)/(dxv(i,j)*dzt(k))  &
           *(dn0(k,i,j)+dn0(k,i,j+1))*.5
      j=nyp-1
      nmass=nmass -  &
           vp(k,i,j)*rtgv(i,j)/(dxv(i,j)*dzt(k))  &
           *(dn0(k,i,j)+dn0(k,i,j+1))*.5
   enddo
enddo
k=2
do j=2,nyp-1
   do i=2,nxp-1
      ppp= ( (pp(k,i,j)+pi0(k,i,j))/cp )**cpor*p00
      prtot=prtot+ppp/(dyu(i,j)*dxv(i,j))
   enddo
enddo
area=(nxp-2)*deltax*(nyp-2)*deltax
endif

!************************************************************
if(jdim==0) then  !for 2D simulations
!  west/east bound
do j=1,1
   do k=2,nzp-1
      i=1
      wmass=wmass +  &
           up(k,i,j)*rtgu(i,j)/(dyu(i,j)*dzt(k))  &
           *(dn0(k,i,j)+dn0(k,i+1,j))*.5
      i=nxp-1
      emass=emass -  &
           up(k,i,j)*rtgu(i,j)/(dyu(i,j)*dzt(k))  &
           *(dn0(k,i,j)+dn0(k,i+1,j))*.5
   enddo
enddo
k=2
do j=1,1
   do i=2,nxp-1
      ppp= ( (pp(k,i,j)+pi0(k,i,j))/cp )**cpor*p00
      prtot=prtot+ppp/dyu(i,j)
   enddo
enddo
area=(nxp-2)*deltax
endif

tmass=wmass+emass+smass+nmass
aintmass=aintmass+tmass*dtlong
print*,'==============================='
print*,' Mass flux - W, E, S, N'
print*,  wmass,emass,smass,nmass
print*, 'total (kg/(m2 s)):   ',tmass/area
print*, 'total (kg/m2):       ',aintmass/area
print*, 'total pr change (pa):',aintmass/area*9.8
print*, 'computed mean press: ',prtot/area
print*,'==============================='

return
END SUBROUTINE mass_flux_bc

