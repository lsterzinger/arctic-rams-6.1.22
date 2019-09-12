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

Subroutine fmrefs1d_isan (ifm,icm,n0,n1  &
                   ,piref,thref,dnref,rtref)

use rconstants

implicit none

integer :: ifm,icm,n1,n0
real, dimension(n0,*) :: piref,thref,dnref,rtref
integer :: k
real :: c1,c2
real, allocatable :: vctr1(:),vctr2(:),vctr3(:),vctr4(:)

!     Interpolate the fine mesh 1-d reference state variables.

allocate(vctr1(n1),vctr2(n1),vctr3(n1),vctr4(n1))

c1 = rgas / (cp - rgas)
c2 = cp * (rgas / p00) ** c1
if (icm >= 1) then
   do k = 1,n1
      vctr1(k) = thref(k,icm) * dnref(k,icm)
      vctr2(k) = rtref(k,icm) * dnref(k,icm)
   enddo

  CALL eintp (dnref(1,icm),dnref(1,ifm),n1,1,1,n1  &
     ,1,1,ifm,1,'t',0,0)
  CALL eintp (vctr1,vctr3,n1,1,1,n1,1,1,ifm,1,'t',0,0)
  CALL eintp (vctr2,vctr4,n1,1,1,n1,1,1,ifm,1,'t',0,0)

   do k = 1,n1
      thref(k,ifm) = vctr3(k) / dnref(k,ifm)
      rtref(k,ifm) = vctr4(k) / dnref(k,ifm)
      piref(k,ifm) = c2 * (dnref(k,ifm) * thref(k,ifm)) ** c1
   enddo
endif

deallocate(vctr1,vctr2,vctr3,vctr4)

return
END SUBROUTINE fmrefs1d_isan

!##############################################################################
Subroutine fmrefs3d_isan (ifm,icm,n1f,n2f,n3f,n1c,n2c,n3c  &
          ,maxiz,maxix,maxiy,nbot,ntop,jd  &
          ,scr1,scr2,vt2da,toptf,toptc,dn0c,dn0f,th0c,th0f  &
          ,pi0f,dn0uf,dn0vf,zt,ztop)

use rconstants

implicit none

integer :: ifm,icm,n1f,n2f,n3f,n1c,n2c,n3c  &
          ,maxiz,maxix,maxiy,nbot,ntop,jd
real, dimension(*) :: scr1,scr2,vt2da,zt
real, dimension(n2f,n3f) :: toptf  
real, dimension(n2c,n3c) :: toptc     
real, dimension(n1f,n2f,n3f) :: dn0f,th0f,pi0f,dn0uf,dn0vf
real, dimension(n1c,n2c,n3c) :: dn0c,th0c   
real :: ztop  
real :: c1,c2,b(1)
integer :: i1,j1,i,j,k

!     Interpolate the fine mesh 3-D reference state variables.

if (icm == 0) return

CALL fmint3 (n1c,n2c,n3c,n1f,n2f,n3f,maxiz,maxix,maxiy  &
     ,ifm,icm,nbot,ntop,jd,1,0,0,'t'  &
     ,dn0c,dn0f,dn0c,dn0f,scr1,scr2,toptf,vt2da,b(1),b(1),b(1))

CALL fmint3 (n1c,n2c,n3c,n1f,n2f,n3f,maxiz,maxix,maxiy  &
     ,ifm,icm,nbot,ntop,jd,1,0,0,'t'  &
     ,th0c,th0f,dn0c,dn0f,scr1,scr2,toptf,vt2da,b(1),b(1),b(1))

c1 = rgas / (cp - rgas)
c2 = cp * (rgas / p00) ** c1
pi0f(1:n1f,1:n2f,1:n3f) = c2 * (dn0f(1:n1f,1:n2f,1:n3f)  &
                            *   th0f(1:n1f,1:n2f,1:n3f) ) ** c1


CALL fillscr (1,maxix,maxiy,1,n2c,n3c,1,1,scr1,toptc)
CALL eintp (scr1,scr2,1,maxix,maxiy,1,n2f,n3f,ifm,2,'t',0,0)
CALL fillvar (1,maxix,maxiy,1,n2f,n3f,1,1,scr2,scr1)

CALL rtgintrp_isan (n1f,n2f,n3f,th0f,scr1,toptf,zt,ztop)
CALL rtgintrp_isan (n1f,n2f,n3f,pi0f,scr1,toptf,zt,ztop)

! Define dn0u and dn0v

do j = 1,n3f
   j1 = min(j+1,n3f)
   do i = 1,n2f
      i1 = min(i+1,n2f)
      do k = 1,n1f
         dn0uf(k,i,j) = .5 * (dn0f(k,i,j) + dn0f(k,i1,j))
         dn0vf(k,i,j) = .5 * (dn0f(k,i,j) + dn0f(k,i,j1))
      enddo
   enddo
enddo

return
END SUBROUTINE fmrefs3d_isan

!##############################################################################
Subroutine fmdn0_isan (ifm,icm,n1f,n2f,n3f,n2c,n3c  &
          ,maxix,maxiy,scr1,scr2,toptf,toptc,dn0f,dn0uf,dn0vf,zt,ztop)

implicit none

integer :: ifm,icm,n1f,n2f,n3f,n2c,n3c,maxix,maxiy
real, dimension(*) :: scr1,scr2,zt
real, dimension(n2f,n3f) :: toptf 
real, dimension(n2c,n3c) :: toptc     
real, dimension(n1f,n2f,n3f) :: dn0f,dn0uf,dn0vf 
real :: ztop
integer :: i,j,i1,j1,k

!     Special vertical interpolation of DN0 must be done after all other
!     3-D reference state and prognostic variables are interpolated.

if (icm == 0) return

CALL fillscr (1,maxix,maxiy,1,n2c,n3c,1,1,scr1,toptc)
CALL eintp (scr1,scr2,1,maxix,maxiy,1,n2f,n3f,ifm,2,'t',0,0)
CALL fillvar (1,maxix,maxiy,1,n2f,n3f,1,1,scr2,scr1)

CALL rtgintrp_isan (n1f,n2f,n3f,dn0f,scr1,toptf,zt,ztop)

! Define dn0u and dn0v

do j = 1,n3f
   j1 = min(j+1,n3f)
   do i = 1,n2f
      i1 = min(i+1,n2f)
      do k = 1,n1f
         dn0uf(k,i,j) = .5 * (dn0f(k,i,j) + dn0f(k,i1,j))
         dn0vf(k,i,j) = .5 * (dn0f(k,i,j) + dn0f(k,i,j1))
      enddo
   enddo
enddo

return
END SUBROUTINE fmdn0_isan

!##############################################################################
Subroutine rtgintrp_isan (n1,n2,n3,fld,vt2da,topt,zt,ztop)

implicit none

integer :: n1,n2,n3
real, dimension(n1,n2,n3) :: fld
real, dimension(n2,n3) :: vt2da,topt
real, dimension(n1) :: zt
real :: ztop
integer :: i,j,k
real, allocatable, dimension(:) :: vctr1,vctr2,vctr3

!     Do special vertical interpolation in case terrain on this grid
!     (topt) is different from what would be interpolated from the
!     coarser grid (vt2da).

allocate(vctr1(n1),vctr2(n1),vctr3(n1))

do j = 1,n3
   do i = 1,n2
      do k = 1,n1
         vctr1(k) = zt(k) * (1. - vt2da(i,j) / ztop) + vt2da(i,j)
         vctr2(k) = zt(k) * (1. - topt(i,j) / ztop) + topt(i,j)
         vctr3(k) = fld(k,i,j)
      enddo
      CALL htint (n1,vctr3,vctr1,n1,fld(1,i,j),vctr2)
   enddo
enddo

deallocate(vctr1,vctr2,vctr3)

return
END SUBROUTINE rtgintrp_isan

!##############################################################################
Subroutine varfile_refstate (n1,n2,n3,thp,pc,pi0,th0,rtp,dn0  &
                 ,dn0u,dn0v,topt,rtgt,zt,ztop,piref,thref,dnref,rtref)

use rconstants
                 
implicit none

integer :: n1,n2,n3,ir,jr,i,j,k,i1,j1
real, dimension(n1,n2,n3) :: thp,pc,pi0,rtp,dn0,th0,dn0u,dn0v
real, dimension(n2,n3)    :: topt,rtgt
real, dimension(n1)       :: zt,piref,thref,rtref,dnref
real :: ztop,topref,c1,c2,c3
real, allocatable :: vctr1(:),vctr2(:)

!                Reference sounding is point with lowest topography

allocate(vctr1(n1),vctr2(n1))
ir=1
jr=1
topref=1.e10
do j=1,n3
   do i=1,n2
      if(topt(i,j).lt.topref) then
         ir=i
         jr=j
         topref=topt(i,j)
      endif
   enddo
enddo

do k=1,n1
   vctr2(k)=zt(k)*(1.-topref/ztop)+topref
enddo

CALL htint2 (n1,thp(1,ir,jr),vctr2,n1,vctr1,zt)
CALL htint2 (n1,rtp(1,ir,jr),vctr2,n1,rtref(1),zt)

do k = 1,n1
   thref(k) = vctr1(k) * (1. + .61 * rtref(k))
enddo
rtref(1) = rtref(2)
thref(1) = thref(2)

piref(1) = pc(1,ir,jr) + g * (vctr2(1) - zt(1))  &
                / (.5 * (thref(1)  &
                + thp(1,ir,jr) * (1. + .61 * rtp(1,ir,jr))))
do k = 2,n1
   piref(k) = piref(k-1) - g * (zt(k)-zt(k-1))  &
             / ( .5 * (thref(k) + thref(k-1)))
enddo

do k = 1,n1
  vctr1(k) = (piref(k) / cp) ** cpor * p00
  dnref(k) = cp * vctr1(k)  &
     / (rgas * thref(k) * piref(k))
enddo

!        Compute 3-D reference state from 1-D reference state

do j=1,n3
  do i=1,n2

    do k=1,n1
      vctr2(k)=zt(k)*rtgt(i,j)+topt(i,j)
    enddo
    CALL htint (n1,piref(1),zt,n1,pi0(1,i,j),vctr2)
    CALL htint (n1,thref(1),zt,n1,th0(1,i,j),vctr2)

    c1=g*2.*(1.-topt(i,j)/ztop)
    c2=(1-cpor)
    c3=cp**c2
    do k=n1-1,1,-1
      pi0(k,i,j)=pi0(k+1,i,j)  &
                +c1*(zt(k+1)-zt(k))/(th0(k,i,j)+th0(k+1,i,j))
    enddo

    do k=1,n1
      dn0(k,i,j)=(c3*p00)/(rgas*th0(k,i,j)*pi0(k,i,j)**c2)
    enddo

  enddo
enddo

do j = 1,n3
   j1 = min(j+1,n3)
   do i = 1,n2
      i1 = min(i+1,n2)
      do k = 1,n1
         dn0u(k,i,j) = .5 * (dn0(k,i,j) + dn0(k,i1,j))
         dn0v(k,i,j) = .5 * (dn0(k,i,j) + dn0(k,i,j1))
      enddo
   enddo
enddo

deallocate(vctr1,vctr2)

return
END SUBROUTINE varfile_refstate
