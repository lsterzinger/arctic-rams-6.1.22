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

Subroutine fldinit (initflg)

use mem_grid
use mem_basic
use mem_turb
use node_mod

implicit none

integer :: initflg

!         finish initializing past time level variables

CALL atob (nxyzp,basic_g(ngrid)%uc(1,1,1),basic_g(ngrid)%up(1,1,1))
CALL atob (nxyzp,basic_g(ngrid)%vc(1,1,1),basic_g(ngrid)%vp(1,1,1))
CALL atob (nxyzp,basic_g(ngrid)%wc(1,1,1),basic_g(ngrid)%wp(1,1,1))
CALL atob (nxyzp,basic_g(ngrid)%pc(1,1,1),basic_g(ngrid)%pp(1,1,1))

CALL tkeinit (nzp,nxp,nyp)

if(initflg == 1) then
   CALL azero (nxyzp,basic_g(ngrid)%wp(1,1,1))
endif

CALL dumset (nzp,nxp,nyp,15,basic_g(ngrid)%uc(1,1,1),'U')
CALL dumset (nzp,nxp,nyp,15,basic_g(ngrid)%vc(1,1,1),'V')
CALL dumset (nzp,nxp,nyp,15,basic_g(ngrid)%wc(1,1,1),'W')
CALL dumset (nzp,nxp,nyp,15,basic_g(ngrid)%up(1,1,1),'U')
CALL dumset (nzp,nxp,nyp,15,basic_g(ngrid)%vp(1,1,1),'V')
CALL dumset (nzp,nxp,nyp,15,basic_g(ngrid)%wp(1,1,1),'W')

return
END SUBROUTINE fldinit

!##############################################################################
Subroutine gridloc_prt ()

use mem_grid
use rconstants

implicit none

real :: gscr(maxgrds,19),centx,centy
integer :: ngr,midx,midy,i,j,ng

do ngr=1,ngrids
   CALL newgrid (ngr)
   gscr(ngr,1)=grid_g(ngr)%glat(1,1)
   gscr(ngr,2)=grid_g(ngr)%glon(1,1)
   gscr(ngr,3)=grid_g(ngr)%glat(nxp,1)
   gscr(ngr,4)=grid_g(ngr)%glon(nxp,1)
   gscr(ngr,5)=grid_g(ngr)%glat(nxp,nyp)
   gscr(ngr,6)=grid_g(ngr)%glon(nxp,nyp)
   gscr(ngr,7)=grid_g(ngr)%glat(1,nyp)
   gscr(ngr,8)=grid_g(ngr)%glon(1,nyp)
   midx = (nxp+1)/2
   midy = (nyp+1)/2
   if (nxp/2 .eq. midx) then
      centx = xm(nxp/2)
   else
      centx = xt(midx)
   endif
   if (nyp/2 .eq. midy) then
      centy = ym(nyp/2)
   else
      centy = yt(midy)
   endif
   CALL xy_ll (gscr(ngr,9),gscr(ngr,10),polelat,polelon  &
      ,centx,centy)

   gscr(ngr,12)=-1000.
   gscr(ngr,13)=1000.
   gscr(ngr,14)=1000.
   gscr(ngr,15)=-1000.
   gscr(ngr,16)=-1000.
   gscr(ngr,17)=1000.
   gscr(ngr,18)=-1000.
   gscr(ngr,19)=1000.
   do j=1,nyp
      do i=1,nxp
         if(j.eq.nyp)  &
            gscr(ngr,12)=max(gscr(ngr,12),grid_g(ngr)%glat(i,j))  ! outer s
         if(j.eq.1)  &
            gscr(ngr,13)=min(gscr(ngr,13),grid_g(ngr)%glat(i,j))  ! outer n
         if(i.eq.1)  &
            gscr(ngr,14)=min(gscr(ngr,14),grid_g(ngr)%glon(i,j))  ! outer w
         if(i.eq.nxp)  & 
            gscr(ngr,15)=max(gscr(ngr,15),grid_g(ngr)%glon(i,j))  ! outer e
         ! provide a conservative estimate so interpolation and
         ! round off error don't result in post process nans
         if(j.eq.3)  &
            gscr(ngr,16)=max(gscr(ngr,16),grid_g(ngr)%glat(i,j))  ! inner s
         if(j.eq.nyp-2)  &
            gscr(ngr,17)=min(gscr(ngr,17),grid_g(ngr)%glat(i,j))  ! inner n
         if(i.eq.3)  &
            gscr(ngr,18)=max(gscr(ngr,18),grid_g(ngr)%glon(i,j))  ! inner w
         if(i.eq.nxp-2)  &
            gscr(ngr,19)=min(gscr(ngr,19),grid_g(ngr)%glon(i,j))  ! inner e
       enddo
   enddo

ENDDO

do ng=1,ngrids
   write(6,100)ng
   write(6,101)gscr(ng,7),gscr(ng,8),gscr(ng,5),gscr(ng,6)
   write(6,102)gscr(ng,9),gscr(ng,10)
   write(6,105)gscr(ng,1),gscr(ng,2),gscr(ng,3),gscr(ng,4)
   print*, ' '
   write(6,106)xtn(1,ng),xtn(nnxp(ng),ng)
   write(6,107)ytn(1,ng),ytn(nnyp(ng),ng)
   write(6,109)deltaxn(ng)
   print*, ' '
   write(6,108)zmn(1,ng),zmn(nnzp(ng)-1,ng)
   write(6,110)zmn(2,ng)-zmn(1,ng)
   print*, ' '
   write(6,113) gscr(ng,12)
   write(6,114) gscr(ng,17)
   write(6,115) gscr(ng,14),gscr(ng,15)
   write(6,116) gscr(ng,18),gscr(ng,19)
   write(6,117) gscr(ng,16)
   write(6,118) gscr(ng,13)

enddo
!
100  FORMAT(/,'---------------Location and Dimensions of GRID',I4  &
        ,'---------------------------',/)
101  FORMAT('  NW lat/lon      NE lat/lon (deg) = '  &
   F8.4,',',F9.4,4X,F8.4,',',F9.4)
102  FORMAT('       Center lat/lon (deg)        =  '  &
   ,11X,F8.4,',',F9.4)
105  FORMAT('  SW lat/lon      SE lat/lon (deg) = '  &
   F8.4,',',F9.4,4X,F8.4,',',F9.4)
106  FORMAT('     West PS coord (km) =',-3PF10.3  &
     ,'               East PS coord (km) =',-3PF10.3)
107  FORMAT('    South PS coord (km) =',-3PF10.3  &
     ,'              North PS coord (km) =',-3PF10.3)
109  FORMAT('         Delta-X (m) =',F10.1)
108  FORMAT('    Bottom coord (m) =',F10.1  &
     ,'            Top coordinate (m) =',F10.1)
110  FORMAT('  Bottom Delta-Z (m) =',F10.1)
!111  format(/,' max N lat=',f9.3,'  min S lat=',f9.3  &
!      ,/,' min W lon=',f9.3,' max E lon=',f9.3 )
113  format('      Outer N lat (deg) =       ',F8.4)
114  format('      Inner N lat (deg) =       ',F8.4)
115  format('      Outer W lon (deg) = ',F9.4,2X,F9.4,'  = Outer E lon (deg)')
116  format('      Inner W lon (deg) = ',F9.4,2X,F9.4,'  = Inner E lon (deg)')
117  format('      Inner S lat (deg) =       ',F8.4)
118  format('      Outer S lat (deg) =       ',F8.4,/)
112  format(80('-'),//)

write(6,112)

return
END SUBROUTINE gridloc_prt

!##############################################################################
Subroutine refs3d (n1,n2,n3,pi0,dn0,dn0u,dn0v,th0,topt,rtgt)

use mem_grid
use mem_scratch
use ref_sounding
use rconstants

implicit none

integer :: n1,n2,n3
real :: pi0(n1,n2,n3),dn0(n1,n2,n3),dn0u(n1,n2,n3)  &
   ,dn0v(n1,n2,n3),th0(n1,n2,n3),topt(n2,n3),rtgt(n2,n3)

integer :: i,j,k
real :: c1,c2,c3

! +---------------------------------------------------------------------
! _    This routine initializes the 3-D reference state arrays from the
!        1-D reference state.
! +---------------------------------------------------------------------

do j=1,nyp
   do i=1,nxp
 
      do k = 1,n1
         vctr2(k) = zt(k) * rtgt(i,j) + topt(i,j)
      enddo
      CALL htint (nzp,pi01dn(1,ngrid),zt,nzp,pi0(1,i,j),vctr2)
      CALL htint (nzp,th01dn(1,ngrid),zt,nzp,th0(1,i,j),vctr2)
      c1 = g * 2. * (1. - topt(i,j) / ztop)

      c2 = 1. - cpor
      c3 = cp ** c2
      do k = n1-1,1,-1
         pi0(k,i,j) = pi0(k+1,i,j)  &
                    + c1 / ((th0(k,i,j) + th0(k+1,i,j)) * dzm(k))
       enddo

      do k = 1,n1
         dn0(k,i,j) = (c3 * p00) / (rgas * th0(k,i,j) * pi0(k,i,j) ** c2)
      enddo

   enddo
enddo

CALL fill_dn0uv (n1,n2,n3,dn0,dn0u,dn0v)

return
END SUBROUTINE refs3d

!##############################################################################
Subroutine fill_dn0uv (n1,n2,n3,dn0,dn0u,dn0v)

implicit none

integer :: n1,n2,n3
real :: dn0(n1,n2,n3),dn0u(n1,n2,n3),dn0v(n1,n2,n3)
integer :: i,j,k,i1,j1

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

return
END SUBROUTINE fill_dn0uv
