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

Subroutine tend0 ()

use mem_grid
use mem_tend
use var_tables
use node_mod
use micphys

implicit none

integer :: n,mxyzp

!     This routine simply sets all tendency arrays to zero.

!     First u,v tendencies

mxyzp = mxp * myp * mzp
CALL azero (mxyzp,tend%ut(1))
CALL azero (mxyzp,tend%vt(1))
CALL azero (mxyzp,tend%wt(1))
CALL azero (mxyzp,tend%pt(1))

!     Now sclrr tendencies

do n = 1,num_scalar(ngrid)
   CALL azero (mxyzp,scalar_tab(n,ngrid)%var_t)
enddo

!******************************************************************************
!Add tendency forcings here for current timestep 

!Convergence forcing
if (ipara.eq.0 .or. mchnum.eq.1) &
  print*,'Ngrid, Time, Domain Max W: ',ngrid,time,vertvel(ngrid)
if(ICONV.ge.1 .and. ICONV.le.5 .and. ICONGR.eq.ngrid) &
  CALL conv_forcing (tend%ut(1),tend%vt(1))

return
END SUBROUTINE tend0

!##############################################################################
Subroutine hadvance (iac)

use mem_grid
use mem_tend
use mem_basic
use mem_scratch
use node_mod

implicit none

integer :: iac
integer :: mxyzp

!     It is here that the Asselin filter is applied.  For the velocities
!     and pressure, this must be done in two stages, the first when
!     IAC=1 and the second when IAC=2.

mxyzp = mxp * myp * mzp
eps = .2

!     For both IAC=1 and IAC=2, call PREDICT for U, V, W, and P.

CALL predict (mxyzp,basic_g(ngrid)%uc(1,1,1)   &
   ,basic_g(ngrid)%up(1,1,1),tend%ut(1),scratch%vt3da(1),iac,dtlv)

if (icorflg .eq. 1 .or. jdim .eq. 1) then
   CALL predict (mxyzp,basic_g(ngrid)%vc(1,1,1)  &
      ,basic_g(ngrid)%vp(1,1,1),tend%vt(1),scratch%vt3da(1),iac,dtlv)
endif

CALL predict (mxyzp,basic_g(ngrid)%wc(1,1,1),basic_g(ngrid)%wp(1,1,1)  &
   ,tend%wt(1),scratch%vt3da(1),iac,dtlv)
CALL predict (mxyzp,basic_g(ngrid)%pc(1,1,1),basic_g(ngrid)%pp(1,1,1)  &
   ,tend%pt(1),scratch%vt3da(1),iac,dtlv)

return
END SUBROUTINE hadvance

!##############################################################################
Subroutine predict (npts,ac,ap,fa,af,iac,dtlp)

use mem_grid
use node_mod

implicit none

integer :: npts,iac,m
real :: epsu,dtlp
real, dimension(*) :: ac,ap,fa,af

!     This routine moves the arrays AC and AP forward by
!     1 time level by adding in the prescribed tendency. It also
!     applies the Asselin filter given by:

!              {AC} = AC + EPS * (AP - 2 * AC + AF)

!     where AP,AC,AF are the past, current and future time levels of A.
!     All IAC=1 does is to perform the {AC} calculation without the AF
!     term present.  IAC=2 completes the calculation of {AC} by adding
!     the AF term only, and advances AC by filling it with input AP
!     values which were already updated in ACOUSTC.
!
epsu = eps
if (ngbegun(ngrid) .eq. 0) epsu = 0.5

if (iac .eq. 1) then
   do m = 1,npts
      ac(m) = ac(m) + epsu * (ap(m) - 2. * ac(m))
   enddo
   return
elseif (iac .eq. 2) then
   do m = 1,npts
      af(m) = ap(m)
      ap(m) = ac(m) + epsu * af(m)
   enddo
endif

do m = 1,npts
  ac(m) = af(m)
enddo

return
END SUBROUTINE predict

!##############################################################################
Subroutine predtr ()

use mem_grid
use var_tables
use node_mod

implicit none

integer :: mxyzp,n

!   -  Step thermodynamic variables from  t  to  t+1.
!   -  Set top, lateral and bottom boundary conditions on some variables
!        if needed.
!   -  call adjustment to assure all positive definite quantities
!        remain positive.
!   -  Rediagnose some thermodynamic quantities for use on the small
!        timestep.

!     Update the scalars and apply lateral, top, and bottom boundary
!     conditions.

mxyzp = mxp * myp * mzp

do n = 1,num_scalar(ngrid)
   CALL update (mxyzp,scalar_tab(n,ngrid)%var_p  &
                    ,scalar_tab(n,ngrid)%var_t, dtlt)
enddo

return
END SUBROUTINE predtr






