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

Subroutine inithh ()

!--------------------------------------------------------
! Initialization of model grid, topography, and fields
! for horizontally homogeneous fields.
!--------------------------------------------------------

use mem_basic
use mem_grid

implicit none

integer :: ifm,icm

!     Arrange the input sounding.

CALL arrsnd ()

!     For GRID 1, compute the 1-D reference state variables, the 3-D
!     reference state variables, the 3-D model fields, the surface
!     layer parameters, and the initial soil model fields.


do ifm = 1,ngrids
   icm = nxtnest(ifm)
   if (icm .eq. 0) then

      CALL newgrid (ifm)
      CALL refs1d ()

      CALL refs3d (nzp,nxp,nyp  &
         ,basic_g(ifm)%pi0  (1,1,1)  ,basic_g(ifm)%dn0  (1,1,1)  &
         ,basic_g(ifm)%dn0u (1,1,1)  ,basic_g(ifm)%dn0v (1,1,1)  &
         ,basic_g(ifm)%th0  (1,1,1)  ,grid_g(ifm)%topt  (1,1)    &
         ,grid_g(ifm)%rtgt  (1,1)                                )

      CALL flds3d (nzp,nxp,nyp  &
            ,basic_g(ifm)%uc  (1,1,1)  ,basic_g(ifm)%vc    (1,1,1)  &
            ,basic_g(ifm)%pi0 (1,1,1)  ,basic_g(ifm)%theta (1,1,1)  &
            ,basic_g(ifm)%thp (1,1,1)  ,basic_g(ifm)%rtp   (1,1,1)  &
            ,basic_g(ifm)%pc  (1,1,1)  ,basic_g(ifm)%rv    (1,1,1)  &
            ,grid_g(ifm)%topt (1,1)    ,grid_g(ifm)%topu   (1,1)    &
            ,grid_g(ifm)%topv (1,1)    ,grid_g(ifm)%rtgt   (1,1)    &
            ,grid_g(ifm)%rtgu (1,1)    ,grid_g(ifm)%rtgv   (1,1)    )

   endif
enddo

return
END SUBROUTINE inithh

!##############################################################################
Subroutine arrsnd ()

use mem_grid
use mem_scratch
use ref_sounding
use rconstants

implicit none

integer :: nnns,k,kk,kkk
real :: toffset,dir,spd,zold1,zold2,tavg,rtss,wt

!     Arrange the input sounding

if (ps(1) .eq. 0.) then

   do nsndg=1,maxsndg
      ps(nsndg) = 0.
      ts(nsndg) = 0.
      rts(nsndg) = 0.
      us(nsndg) = 0.
      vs(nsndg) = 0.
   enddo

   open(1,file='SOUND_IN',status='old',form='formatted')
   do nsndg=1,maxsndg
      read(1,*,end=1999) ps(nsndg),ts(nsndg),rts(nsndg),us(nsndg),vs(nsndg)
      if(ps(nsndg).le.0.) go to 1999
   enddo
1999    continue
   close(1)
endif

toffset=-1.e30 !Variable initialized
zold2=-1.e30   !Variable initialized

if(itsflg.eq.0)then
   toffset=273.15
elseif(itsflg.eq.1)then
   toffset=0.
endif

do nsndg=1,maxsndg
   nnns=nsndg
   if(ps(nsndg).eq.0.)go to 300
   if(us(nsndg).ne.9999.)then
      if(iusflg.ne.0)then
         dir=us(nsndg)
         spd=vs(nsndg)
         us(nsndg)=-spd*sin(pi180*dir)
         vs(nsndg)=-spd*cos(pi180*dir)
      endif
   endif

   if(ipsflg.eq.0)then
!     pressure given in millibars
      ps(nsndg)=ps(nsndg)*1.e2
   elseif(ipsflg.eq.1)then
!     pressure is height in meters with PS(1)=surface pressure

!  If sounding moisture is expressed as a mixing ratio (IRTSFLG=2),
!     take advantage of knowing virtual temperature effect when
!     integrating hydrostatically to get sounding pressure.
!
      if(irtsflg.eq.2)then
         vctr4(nsndg)=(1.+.61*rts(nsndg)*1.e-3)
      else
         vctr4(nsndg)=1.
      endif
      if(nsndg.eq.1)then
         ps(nsndg)=ps(nsndg)*100.
         zold2=0.
      else
         zold1=zold2
         zold2=ps(nsndg)
         if(itsflg.eq.0.or.itsflg.eq.1)then
            tavg = 0.5*( (ts(nsndg)+toffset)*vctr4(nsndg)  &
                 +        ts(nsndg-1)*vctr4(nsndg-1) )
            ps(nsndg)=ps(nsndg-1)*exp(-g*(zold2-zold1)/(rgas*tavg))
         elseif(itsflg.eq.2)then
            tavg=(ts(nsndg)*vctr4(nsndg)  &
                + ts(nsndg-1)*vctr4(nsndg-1)*p00k/ps(nsndg-1)**rocp)*.5
            ps(nsndg)=(ps(nsndg-1)**rocp-g*(zold2-zold1)*p00k  &
               /(cp*tavg))**cpor
         endif
      endif
   else
      write(6,131)ipsflg
131        format(' PRESSURE TYPE (IPSFLG=',I2,') NOT IMPLEMENTED')
      stop
   endif

   if(itsflg.eq.0)then
!     Temperature in degrees celsius
      ts(nsndg)=ts(nsndg)+273.15
   elseif(itsflg.eq.1)then
!     Temperature in degrees kelvin
   elseif(itsflg.eq.2)then
!     Temperature is potential temperature in kelvin
      ts(nsndg)=(ps(nsndg)*p00i)**rocp*ts(nsndg)
   else
      write(6,124)itsflg
124        format(' TEMPERATURE TYPE (ITSFLG=',I2,') NOT KNOWN')
      stop
   endif

   if(irtsflg.eq.0)then
!     Humidity given as dew point in degrees celsius
      CALL mrsl (1,ps(nsndg),rts(nsndg)+273.15,rts(nsndg))
   elseif(irtsflg.eq.1)then
!     Humidity given as dew point in degrees kelvin
      CALL mrsl (1,ps(nsndg),rts(nsndg),rts(nsndg))
   elseif(irtsflg.eq.2)then
!     Humidity given as mixing ratio in g/kg
      rts(nsndg)=rts(nsndg)*1.e-3
   elseif(irtsflg.eq.3)then
!     Humidity given as relative humidity in percent
      CALL mrsl (1,ps(nsndg),ts(nsndg),rtss)
      rts(nsndg)=rtss*rts(nsndg)*.01
   elseif(irtsflg.eq.4)then
!     Humidity given as dew point depression in kelvin
      CALL mrsl (1,ps(nsndg),ts(nsndg)-rts(nsndg),rts(nsndg))
   else
      write(6,125) irtsflg
125        format(' HUMIDITY TYPE (IRTSFLG=',I2,') NOT KNOWN')
      stop
   endif
enddo
300  continue

nsndg=nnns-1
do k=1,nsndg
   if(us(k).eq.9999.)then
      do kk=k,1,-1
         if(us(kk).ne.9999.)go to 17
      enddo
17         continue
      do kkk=k,nsndg,1
         if(us(kkk).ne.9999.) go to 18
      enddo
18         continue
      wt=(ps(k)-ps(kk))/(ps(kkk)-ps(kk))
      us(k)=us(kk)+wt*(us(kkk)-us(kk))
      vs(k)=vs(kk)+wt*(vs(kkk)-vs(kk))
   endif
enddo

!     compute height levels of input sounding.

do k=2,nsndg
   hs(k)=hs(k-1)-rgas*.5  &
      *(ts(k)*(1.+.61*rts(k))+ts(k-1)*(1.+.61*rts(k-1)))  &
      *(log(ps(k))-log(ps(k-1)))/g
enddo

if(hs(nsndg).lt.zt(nzp)) then
   write(6,1) hs(nsndg),zt(nzp)
1    format('    Input sounding needs to go higher ! !', /,  &
    '      Sounding top (m) = ',F12.2,'  Model top (m) = ',F12.2)
   stop 'ARRSND'
endif

do k=1,nsndg
   thds(k)=ts(k)*(p00/ps(k))**rocp
enddo

!Saleeby(2014): If simulation is 2D and sounding has V-winds, set to zero
if(jdim == 0) then
 print*,'Simulation is 2D, but sounding contained V-winds.'
 print*,'Will set V-winds to zero in this case.'
 do k=1,nsndg
    vs(k) = 0.0
 enddo
endif

return
END SUBROUTINE arrsnd

!##############################################################################
Subroutine refs1d ()

use mem_grid
use mem_scratch
use micphys
use ref_sounding
use rconstants

implicit none
! +---------------------------------------------------------------------
! \   This routine computes the reference state sounding on the model
! \     sigma-z levels from input sounding defined on pressure levels.
! +---------------------------------------------------------------------

integer :: k

if (ztn(nnzp(ngrid),ngrid) .gt. hs(nsndg)) then
   print*,' !!! Input sounding is not high enough !!!'
   print*,' !!! Sounding height: ',hs(nsndg)
   print*,' !!! Model top      : ',ztn(nnzp(ngrid),ngrid)
   stop 'refs1d'
endif

CALL htint (nsndg,thds,hs,nnzp(ngrid),vctr1,ztn(1,ngrid))
CALL htint (nsndg,us,hs,nnzp(ngrid),u01dn(1,ngrid),ztn(1,ngrid))
CALL htint (nsndg,vs,hs,nnzp(ngrid),v01dn(1,ngrid),ztn(1,ngrid))

if (level .ge. 1) then
   CALL htint (nsndg,rts,hs,nnzp(ngrid),rt01dn(1,ngrid),ztn(1,ngrid))
else
   do k = 1,nnzp(ngrid)
      rt01dn(k,ngrid) = 0.
   enddo
endif

do k = 1,nnzp(ngrid)
   th01dn(k,ngrid) = vctr1(k) * (1. + .61 * rt01dn(k,ngrid))
enddo
u01dn(1,ngrid) = u01dn(2,ngrid)
v01dn(1,ngrid) = v01dn(2,ngrid)
rt01dn(1,ngrid) = rt01dn(2,ngrid)
th01dn(1,ngrid) = th01dn(2,ngrid)

pi01dn(1,ngrid) = cp * (ps(1) * p00i) ** rocp  &
   + g * (hs(1) - ztn(1,ngrid))  &
   / (.5 * (th01dn(1,ngrid) + thds(1) * (1. + .61 * rts(1)) ) )
do k = 2,nnzp(ngrid)
   pi01dn(k,ngrid) = pi01dn(k-1,ngrid) - g / (dzmn(k-1,ngrid)  &
      * .5 * (th01dn(k,ngrid) + th01dn(k-1,ngrid)))
enddo

do k = 1,nnzp(ngrid)
   vctr4(k) = (pi01dn(k,ngrid) / cp) ** cpor * p00
   dn01dn(k,ngrid) = cp * vctr4(k)  &
      / (rgas * th01dn(k,ngrid) * pi01dn(k,ngrid))
enddo

return
END SUBROUTINE refs1d

!##############################################################################
Subroutine flds3d (n1,n2,n3,uc,vc,pi0,theta,thp,rtp,pc,rv  &
   ,topt,topu,topv,rtgt,rtgu,rtgv)

use mem_grid
use mem_scratch
use ref_sounding
use rconstants
use micphys

implicit none

integer :: n1,n2,n3
real :: pi0(n1,n2,n3),thp(n1,n2,n3),theta(n1,n2,n3)  &
   ,rtp(n1,n2,n3),pc(n1,n2,n3),rv(n1,n2,n3)  &
   ,uc(n1,n2,n3),vc(n1,n2,n3),topt(n2,n3),topu(n2,n3),topv(n2,n3)  &
   ,rtgt(n2,n3),rtgu(n2,n3),rtgv(n2,n3)

integer :: i,j,k
real :: qlatu,qlonu,qlatv,qlonv,dummy
real, dimension(nzpmax) :: p0,temp,rvls,rc

! +---------------------------------------------------------------------
!     This routine initializes the 3-D velocity and thermodynamic
!       fields from the 1-D reference state sounding.
! +---------------------------------------------------------------------

do j=1,nyp
   do i=1,nxp

      do k=1,nzp
         vctr11(k) = zt(k) * rtgt(i,j) + topt(i,j)
         vctr12(k) = zt(k) * rtgu(i,j) + topu(i,j)
         vctr13(k) = zt(k) * rtgv(i,j) + topv(i,j)
      enddo

      CALL htint (nzp, u01dn(1,ngrid),zt,nzp,vctr5,vctr12)
      CALL htint (nzp, v01dn(1,ngrid),zt,nzp,vctr6,vctr13)
      CALL htint (nzp,th01dn(1,ngrid),zt,nzp,vctr3,vctr11)
      if(level.ge.1)  &
      CALL htint (nzp,rt01dn(1,ngrid),zt,nzp,rtp(1,i,j),vctr11)

! If sounding winds are to be interpreted as eastward (U) and
! northward (V) components, rotate winds from geographic to
! polar stereographic orientation

      if (ihtran .eq. 1) then

         CALL xy_ll (qlatu,qlonu,polelat,polelon  &
            ,xm(i),yt(j))
         CALL xy_ll (qlatv,qlonv,polelat,polelon  &
            ,xt(i),ym(j))

         do k = 1,nzp
            CALL uevetouv (uc(k,i,j),dummy,vctr5(k),vctr6(k)  &
                   ,qlatu,qlonu,polelat,polelon)
            CALL uevetouv (dummy,vc(k,i,j),vctr5(k),vctr6(k)  &
                   ,qlatv,qlonv,polelat,polelon)
         enddo
      else
         do k = 1,nzp
            uc(k,i,j)=vctr5(k)
            vc(k,i,j)=vctr6(k)
         enddo
      endif

      if(level.ge.1)then
         do k=1,nzp
            thp(k,i,j)=vctr3(k)/(1.+.61*rtp(k,i,j))
         enddo
      else
         do k=1,nzp
            thp(k,i,j)=vctr3(k)
         enddo
      endif
      do k=1,nzp
         theta(k,i,j)=thp(k,i,j)
         pc(k,i,j)=0.
      enddo

      if(level.eq.1)then
         do k=1,nzp
            rv(k,i,j)=rtp(k,i,j)
         enddo
      endif

      if(level.ge.2)then
         do k=1,nzp
            p0(k)=(pi0(k,i,j)/cp)**cpor*p00
            temp(k)=pi0(k,i,j)*thp(k,i,j)/cp
         enddo
         CALL mrsl (nzp,p0(1),temp(1),rvls(1))
         do k=1,nzp
            rc(k)=max(0.,rtp(k,i,j)-rvls(k))
            thp(k,i,j)=theta(k,i,j)  &
                /(1.+(aklv*rc(k))/max(temp(k),253.))
            rv(k,i,j)=rtp(k,i,j)-rc(k)
         enddo
      endif

   enddo
enddo

return
END SUBROUTINE flds3d

