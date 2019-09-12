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

Subroutine varf_update (iswap,ifileok,initflag)

use mem_leaf
use mem_varinit
use mem_basic
use mem_grid
use mem_scratch
use micphys
use hdf5_utils

implicit none

integer :: ifileok,initflag,iswap

!---------------------------------------------------------------+
!    "Variable initialization"  initialization routines
!---------------------------------------------------------------+
logical :: there
integer :: iver_var,nc,iyearx,imonthx,idatex,ihourx  &
          ,nxpx,nypx,nzpx,kk,ii,jj
integer :: ndims,idims(4)          
real :: rlatx,wlon1x,deltaxx,deltazx,dzratx,dzmaxx
character(len=7) :: cgrid
character(len=strl1) :: flnm

!      Check and see what we are doing. If it is initial time, read
!        fields into regular arrays. If not, see if nudging will be done
!        on this grid if it is a nested grid.
if (ngrid > 1 .and. tnudcent+tnudtop < .001 .and. initflag == 0) return

! Put new fields into varinit future arrays. If iswap == 1, 
!     swap future into past first

if (iswap == 1) then

do jj=1,nnyp(ngrid)
 do ii=1,nnxp(ngrid)
  do kk=1,nnzp(ngrid)
   !The first 5 are needed for all nudged simulations
   varinit_g(ngrid)%varup(kk,ii,jj) = varinit_g(ngrid)%varuf(kk,ii,jj)
   varinit_g(ngrid)%varvp(kk,ii,jj) = varinit_g(ngrid)%varvf(kk,ii,jj)
   varinit_g(ngrid)%varpp(kk,ii,jj) = varinit_g(ngrid)%varpf(kk,ii,jj)
   varinit_g(ngrid)%vartp(kk,ii,jj) = varinit_g(ngrid)%vartf(kk,ii,jj)
   varinit_g(ngrid)%varrp(kk,ii,jj) = varinit_g(ngrid)%varrf(kk,ii,jj)
   !These next 2 are for condensate nudging from History-Varfiles
   if (nud_cond == 1) then
    varinit_g(ngrid)%varrph(kk,ii,jj) = varinit_g(ngrid)%varrfh(kk,ii,jj)
    varinit_g(ngrid)%varcph(kk,ii,jj) = varinit_g(ngrid)%varcfh(kk,ii,jj)
   endif
  enddo
 enddo
enddo

endif

write(cgrid,'(a2,i1,a3)') '-g',ngrid,'.h5'
nc=len_trim(fnames_varf(nvarffl))
flnm=fnames_varf(nvarffl)(1:nc-4)//trim(cgrid)
inquire(file=trim(flnm),exist=there)

! Gotta have grid 1...
if (.not.there .and. ngrid == 1) then
   print*
   print*,'No grid 1 varfile found: ',trim(flnm)
   print*
   stop 'no grid 1 varfile'
endif

if(there) then
   ifileok=1
else
   ifileok=0
   return
endif

! Read the varfile fields into the "future" varinit arrays. These will be 
!   swapped to the past arrays when needed.
   CALL shdf5_open (flnm,'R')
   ndims=1 ; idims(1)=1
   CALL shdf5_irec ('version',ivars=iver_var)
   CALL shdf5_irec ('year',ivars=iyearx)
   CALL shdf5_irec ('month',ivars=imonthx)
   CALL shdf5_irec ('day',ivars=idatex)
   CALL shdf5_irec ('hour',ivars=ihourx)
   CALL shdf5_irec ('nx',ivars=nxpx)
   CALL shdf5_irec ('ny',ivars=nypx)
   CALL shdf5_irec ('nz',ivars=nzpx)
   CALL shdf5_irec ('polelat',rvars=rlatx)
   CALL shdf5_irec ('polelon',rvars=wlon1x)
   CALL shdf5_irec ('dx',rvars=deltaxx)
   CALL shdf5_irec ('dz',rvars=deltazx)
   CALL shdf5_irec ('dzrat',rvars=dzratx)
   CALL shdf5_irec ('dzmax',rvars=dzmaxx)

if(nxp.ne.nxpx.or.  &
   nyp.ne.nypx.or.  &
   nzp.ne.nzpx.or.  &
   abs(deltax-deltaxx).gt..001.or.  &
   abs(deltaz-deltazx).gt..001.or.  &
   abs(dzrat-dzratx).gt..001.or.  & 
   abs(dzmax-dzmaxx).gt..001.or.  &
   abs(polelat-rlatx).gt..001.or.  &
   abs(polelon-wlon1x).gt..001) then
   
   print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   print*,'!!    GRID MISMATCH BETWEEN VARFILE AND NAMELIST !'
   print*,'!!          RUN IS STOPPED                       !'
   print*,'!!  File:',trim(flnm)
   print*,'!!  File, Namelist values for grid:',ngrid
   print*,'!!  nxp:',nxpx,nxp
   print*,'!!  nyp:',nypx,nyp
   print*,'!!  nzp:',nzpx,nzp
   print*,'!!  deltax:',deltaxx,deltax
   print*,'!!  deltaz:',deltazx,deltaz
   print*,'!!  dzrat:',dzratx,dzrat
   print*,'!!  dzmax:',dzmaxx,dzmax
   print*,'!!  polelat:',rlatx,polelat
   print*,'!!  polelon:',wlon1x,polelon
   PRINT*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   stop 'bad-vfile'
endif

   ndims=3 ; idims(1)=nnxp(ngrid); idims(2)=nnyp(ngrid); idims(3)=nnzp(ngrid)
   CALL shdf5_irec ('UP',rvara=scratch%scr1)
   CALL unarrange (nzp,nxp,nyp,scratch%scr1,varinit_g(ngrid)%varuf)
   CALL shdf5_irec ('VP',rvara=scratch%scr1)
   CALL unarrange (nzp,nxp,nyp,scratch%scr1,varinit_g(ngrid)%varvf)
   CALL shdf5_irec ('PI',rvara=scratch%scr1)
   CALL unarrange (nzp,nxp,nyp,scratch%scr1,varinit_g(ngrid)%varpf)
   CALL shdf5_irec ('THETA',rvara=scratch%scr1)
   CALL unarrange (nzp,nxp,nyp,scratch%scr1,varinit_g(ngrid)%vartf)
   CALL shdf5_irec ('RV',rvara=scratch%scr1)
   CALL unarrange (nzp,nxp,nyp,scratch%scr1,varinit_g(ngrid)%varrf)

!For Condensate nudging from History-Varfiles, RV = RTP so we load
!the RV variable into the nudging variable VARRFH
if(nud_cond == 1) then
   CALL shdf5_irec ('RV',rvara=scratch%scr1)
   CALL unarrange (nzp,nxp,nyp,scratch%scr1,varinit_g(ngrid)%varrfh)
   CALL shdf5_irec ('COND',rvara=scratch%scr1)
   CALL unarrange (nzp,nxp,nyp,scratch%scr1,varinit_g(ngrid)%varcfh)
endif

varinit_g(ngrid)%varrf(1:nzp,1:nxp,1:nyp)=  &
           max(1.e-8,varinit_g(ngrid)%varrf(1:nzp,1:nxp,1:nyp) )

!Extract soil/snow data from the varfile. Ignore other 2D fields for now.
!Do this for initialization only
if(initflag == 1 .and. iver_var == 3) then
      ndims=2 ; idims(1)=nnxp(ngrid); idims(2)=nnyp(ngrid)
      CALL shdf5_irec ('SOILMOIST1',rvara=leaf_g(ngrid)%soil_moist_bot)
      CALL shdf5_irec ('SOILMOIST2',rvara=leaf_g(ngrid)%soil_moist_top)
      CALL shdf5_irec ('SOILTEMP1', rvara=leaf_g(ngrid)%soil_temp_bot)
      CALL shdf5_irec ('SOILTEMP2', rvara=leaf_g(ngrid)%soil_temp_top)
      CALL shdf5_irec ('SNOWMASS',  rvara=leaf_g(ngrid)%snow_mass)
      CALL shdf5_irec ('SNOWDEPTH', rvara=leaf_g(ngrid)%snow_depth)    
endif

CALL shdf5_close ()

! Find the reference state

if(initflag == 1 .and. ngrid == 1)  &
     CALL varref (nzp,nxp,nyp &
         ,varinit_g(ngrid)%vartf(1,1,1) ,varinit_g(ngrid)%varpf(1,1,1)  &
         ,basic_g(ngrid)%pi0(1,1,1),     basic_g(ngrid)%th0(1,1,1)  &
         ,varinit_g(ngrid)%varrf(1,1,1), basic_g(ngrid)%dn0(1,1,1)  &
         ,basic_g(ngrid)%dn0u(1,1,1),    basic_g(ngrid)%dn0v(1,1,1)  &
         ,varinit_g(ngrid)%varuf(1,1,1), varinit_g(ngrid)%varvf(1,1,1)  &
         ,grid_g(ngrid)%topt(1,1),       grid_g(ngrid)%rtgt(1,1)  &
         ,level)

do jj=1,nnyp(ngrid)
 do ii=1,nnxp(ngrid)
  do kk=1,nnzp(ngrid)
    varinit_g(ngrid)%varpf(kk,ii,jj) =  &
        varinit_g(ngrid)%varpf(kk,ii,jj) - basic_g(ngrid)%pi0(kk,ii,jj)
  enddo
 enddo
enddo

! If this is an initialization, put data into regular arrays

if(initflag == 1 ) then
   CALL atob (nxyzp,varinit_g(ngrid)%varuf(1,1,1),basic_g(ngrid)%uc(1,1,1))
   CALL atob (nxyzp,varinit_g(ngrid)%varvf(1,1,1),basic_g(ngrid)%vc(1,1,1))
   CALL atob (nxyzp,varinit_g(ngrid)%varpf(1,1,1),basic_g(ngrid)%pc(1,1,1))
   CALL atob (nxyzp,varinit_g(ngrid)%vartf(1,1,1),basic_g(ngrid)%thp(1,1,1))
   CALL atob (nxyzp,varinit_g(ngrid)%varrf(1,1,1),basic_g(ngrid)%rtp(1,1,1))
endif

return
END SUBROUTINE varf_update

!##############################################################################
Subroutine varref (n1,n2,n3,thp,pc,pi0,th0,rtp,dn0,dn0u,dn0v,uc  &
                 ,vc,topt,rtgt,level)

use mem_grid
use ref_sounding
use mem_scratch
use rconstants
                 
implicit none

integer :: n1,n2,n3,level,i,j,k
real, dimension(n1,n2,n3) :: thp,pc,pi0,rtp,dn0,dn0u,dn0v,uc,vc,th0
real, dimension(n2,n3) :: topt,rtgt

!                Reference sounding is point with lowest topography
topref=1.e10
do j=1,nyp
   do i=1,nxp
      if(topt(i,j).lt.topref) then
         iref=i
         jref=j
         topref=topt(i,j)
      endif
   enddo
enddo

!  Set up 1-D reference state

do k=1,nzp
   vctr2(k)=ztn(k,ngrid)*(1.-topref/ztop)+topref
enddo
CALL htint2 (nzp,thp(1,iref,jref),vctr2,nzp,vctr1,zt)
CALL htint2 (nzp,uc(1,iref,jref),vctr2,nzp,u01dn(1,ngrid),zt)
CALL htint2 (nzp,vc(1,iref,jref),vctr2,nzp,v01dn(1,ngrid),zt)
if (level >= 1) then
   CALL htint2 (nzp,rtp(1,iref,jref),vctr2,nzp,rt01dn(1,ngrid),zt)
else
   rt01dn(1:nzp,ngrid) = 0.
endif

do k = 1,nzp
   th01dn(k,ngrid) = vctr1(k) * (1. + .61 * rt01dn(k,ngrid))
enddo
u01dn(1,ngrid) = u01dn(2,ngrid)
v01dn(1,ngrid) = v01dn(2,ngrid)
rt01dn(1,ngrid) = rt01dn(2,ngrid)
th01dn(1,ngrid) = th01dn(2,ngrid)

pi01dn(1,ngrid) = pc(1,iref,jref) + g * (vctr2(1) - zt(1))  &
   / (.5 * (th01dn(1,ngrid)  &
   + thp(1,iref,jref) * (1. + .61 * rtp(1,iref,jref))))
do k = 2,nzp
  pi01dn(k,ngrid) = pi01dn(k-1,ngrid) - g / (dzm(k-1) * .5  &
     * (th01dn(k,ngrid) + th01dn(k-1,ngrid)))
enddo

do k = 1,nzp
  vctr4(k) = (pi01dn(k,ngrid) / cp) ** cpor * p00
  dn01dn(k,ngrid) = cp * vctr4(k)  &
     / (rgas * th01dn(k,ngrid) * pi01dn(k,ngrid))
enddo

!        Compute 3-D reference state from 1-D reference state

CALL refs3d (nzp,nxp,nyp,pi0,dn0,dn0u,dn0v,th0,topt,rtgt)

return
END SUBROUTINE varref

