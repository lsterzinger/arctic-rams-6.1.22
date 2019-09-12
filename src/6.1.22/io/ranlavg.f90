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

Subroutine anlavg (n1,n2,n3)

use var_tables
use mem_scratch
use mem_turb
use mem_basic
use node_mod
use mem_grid
use io_params

implicit none

integer :: nv,izero,n1,n2,n3
real, pointer :: v_p, vm_p
integer, save :: ncall=0,navg
real, save :: avgtim1,frq,avgtim2,timem,timecent

! This routine averages all of the analysis variables over
! ANLAVG. It also accumulates precipitation for all of the
! microphysics variables since the last analysis write.

if(avgtim == 0.)return
if(frqmean == 0. .and. frqboth == 0.)return

if(ncall.eq.0) then
   ncall=1
!  Calculate the number of timestep to average over. Note that if
!    AVGTIM>0, then need to have an equal number of timesteps on
!    both sides of the analysis write, else can have an odd number
!    of timesteps if the write represent an average at the end
!    of the averaging period.
   avgtim1=avgtim
   if(avgtim.gt.0)then
     avgtim=max(avgtim,2.*dtlongn(1))
     avgtim=float(nint(avgtim/dtlongn(1)/2.))*dtlongn(1)*2.
   elseif(avgtim.lt.0)then
     avgtim=min(avgtim,-dtlongn(1))
     avgtim=float(nint(avgtim/dtlongn(1)))*dtlongn(1)
   endif
   if(avgtim.ne.avgtim1)then
     print*,' '
     print*,'***************************************************'
     print*,'***Changing AVGTIM so multiple of DTLONG',avgtim
     print*,'***************************************************'
     print*,' '
   endif
   navg=nint(abs(avgtim)/dtlongn(1))
   print*,' '
   print*,'****************************************************'
   print*,'**** AVERAGING OVER THESE # OF TIMESTEPS', navg
   print*,'**** AVERAGING TIME=',avgtim
   print*,'****************************************************'
   print*,' '

!  Choose a output frequency

   if(frqmean.gt.0.0.and.frqboth.gt.0.)then
      frq=min(frqmean,frqboth)
   elseif(frqmean.gt.0.0.and.frqboth.eq.0.)then
      frq=frqmean
   elseif(frqmean.eq.0.0.and.frqboth.gt.0.)then
      frq=frqboth
   endif

endif

izero=0
avgtim2=abs(avgtim)/2.
timem=time+0.01  
timecent=float(nint(timem/frq))*frq

! No need to execute if not in the averaging interval.
if(avgtim.gt.0.0)then
  if(timem.lt.avgtim2)return
  if(timem.lt.timecent-avgtim2.or.timem.gt.timecent+avgtim2) return
elseif(avgtim.lt.0.0)then
  if(mod(timem-avgtim,frq).gt.-avgtim)return
endif

! Zero out the averages before accumulating.
if(avgtim.gt.0.0.and.mod(timem+avgtim2,frq).lt.dtlongn(1)) izero=1
if(avgtim.lt.0.0.and.mod(timem-avgtim,frq).lt.dtlongn(1)) izero=1

! Implement call to NON_SCALAR_BC before updating averages.
! This sets all BC's for non-scalar (non-advected) quantities such
! as THETA, RV, and variables in LEAF, RADIATION, TURBULENCE.
CALL non_scalar_bc (n1,n2,n3,nzg,nzs)

! Loop through the main variable table

do nv = 1,num_var(ngrid)

   if (vtab_r(nv,ngrid)%imean == 1) then
   
      v_p => vtab_r(nv,ngrid)%var_p
      vm_p=> vtab_r(nv,ngrid)%var_m
      
      if(iprntstmt>=2)print*,'anlavg:',ngrid,nv,num_var(ngrid) &
                      ,trim(vtab_r(nv,ngrid)%name)
      if(izero == 1) CALL azero (vtab_r(nv,ngrid)%npts,vm_p)
      CALL average (vtab_r(nv,ngrid)%npts,vm_p,v_p,navg)
      
   endif

enddo

return
END SUBROUTINE anlavg

!##############################################################################
Subroutine average (m,av,v,navg)

implicit none

integer :: m, navg
real :: av(m),v(m)
integer :: i

do i=1,m
  av(i)=av(i)+v(i)/float(navg)
enddo

return
END SUBROUTINE average
