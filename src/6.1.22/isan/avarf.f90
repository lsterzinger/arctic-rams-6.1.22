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

Subroutine makevarf (ng)

use isan_coms
use mem_grid

implicit none

integer :: ng

!---------------------------------------------------------------+
!    Interpolate model grid from isentropic/sigmaz/surface data
!---------------------------------------------------------------+

!            Vertically interpolate isentropic data to sigma-z levels

CALL isnsig (nnzp(ng),nnxp(ng),nnyp(ng),is_grids(ng)%rr_u(1,1,1)  &
        ,is_grids(ng)%rr_v(1,1,1),is_grids(ng)%rr_t(1,1,1)  &
        ,is_grids(ng)%rr_r(1,1,1),grid_g(ng)%topt(1,1),ztn(1,ng),ztop)

!            Compute Exner func on model sigma-z surfaces
!              and change relative humidity to mixing ratio.

CALL vshyd (nnzp(ng),nnxp(ng),nnyp(ng),is_grids(ng)%rr_p(1,1,1)  &
     ,is_grids(ng)%rr_t(1,1,1)       ,is_grids(ng)%rr_r(1,1,1)  &
     ,grid_g(ng)%topt(1,1),grid_g(ng)%rtgt(1,1),ztn(1,ng))

!          Combine surface analysis with the upper air data.

CALL visurf (nnzp(ng),nnxp(ng),nnyp(ng) ,is_grids(ng)%rr_u(1,1,1)  &
      ,is_grids(ng)%rr_v(1,1,1)        ,is_grids(ng)%rr_t(1,1,1)  &
      ,is_grids(ng)%rr_r(1,1,1)        ,is_grids(ng)%rr_p(1,1,1)  &
      ,grid_g(ng)%topt(1,1),grid_g(ng)%rtgt(1,1),ztn(1,ng))

is_grids(ng)%rr_soilmoist1 (1:nnxp(ng),1:nnyp(ng)) =rs_soilmoist1 (1:nnxp(ng),1:nnyp(ng))
is_grids(ng)%rr_soilmoist2 (1:nnxp(ng),1:nnyp(ng)) =rs_soilmoist2 (1:nnxp(ng),1:nnyp(ng))
is_grids(ng)%rr_soiltemp1  (1:nnxp(ng),1:nnyp(ng)) =rs_soiltemp1  (1:nnxp(ng),1:nnyp(ng))
is_grids(ng)%rr_soiltemp2  (1:nnxp(ng),1:nnyp(ng)) =rs_soiltemp2  (1:nnxp(ng),1:nnyp(ng))
is_grids(ng)%rr_snowmass   (1:nnxp(ng),1:nnyp(ng)) =rs_snowmass   (1:nnxp(ng),1:nnyp(ng))
is_grids(ng)%rr_snowdepth  (1:nnxp(ng),1:nnyp(ng)) =rs_snowdepth  (1:nnxp(ng),1:nnyp(ng))

!          average the velocities to the correct points in the stagger
!             and rotate for polar stereographic transformation.

CALL varuv (nnzp(ng),nnxp(ng),nnyp(ng),is_grids(ng)%rr_u(1,1,1)  &
                      ,is_grids(ng)%rr_v(1,1,1))

return
END SUBROUTINE makevarf

!##############################################################################
Subroutine isnsig (n1,n2,n3,uu,vv,tt,rr,topt,zt,ztop)
   
use isan_coms
use rconstants

implicit none

integer :: n1,n2,n3
real, dimension(n1,n2,n3) :: uu,vv,tt,rr
real, dimension(n2,n3) :: topt
real, dimension(n1) :: zt
real :: ztop
real, allocatable :: v1(:),v2(:),v3(:),v4(:),v5(:),v6(:)
integer :: i,j,k,nki,ki,kbeg,kend
real :: rtg,wtsz

allocate (v1(n1),v2(nisn),v3(nisn),v4(nisn),v5(nisn),v6(nisn))

do j=1,n3
   do i=1,n2

      rtg=(1.-topt(i,j)/ztop)
      do k=1,n1
         v1(k)=topt(i,j)+zt(k)*rtg
      enddo

         nki=0
         do ki=1,nisn
            if(pi_p(i,j,ki).lt.1e20.and.pi_s(i,j,ki).lt.1e20  &
                 .and.pi_u(i,j,ki).lt.1e20.and.pi_v(i,j,ki).lt.1e20)  &
                 then
               nki=nki+1
               v2(nki)=(pi_s(i,j,ki)-cp*levth(ki)  &
                    *(pi_p(i,j,ki)*p00i)**rocp)/g                    
               v3(nki)=pi_u(i,j,ki)
               v4(nki)=pi_v(i,j,ki)
               v5(nki)=levth(ki)
               v6(nki)=pi_r(i,j,ki)
            endif
         enddo

         CALL htint (nki,v3,v2,n1,uu(1,i,j),v1)
         CALL htint (nki,v4,v2,n1,vv(1,i,j),v1)
         CALL htint (nki,v5,v2,n1,tt(1,i,j),v1)
         CALL htint (nki,v6,v2,n1,rr(1,i,j),v1)

      kbeg=n1+1
      kend=n1+1
      do k=1,nsigz
         if(zt(k).gt.hybbot) then
            kbeg=k-1
            exit
         endif
      enddo

      do k=1,nsigz
         if(zt(k).gt.hybtop) then
            kend=k
            exit
         endif
      enddo

      do k=1,nsigz

         if(k.lt.kbeg) then
            wtsz=sigzwt
         elseif(k.gt.kend) then
            wtsz=0.
         elseif(k.ge.kbeg.and.k.le.kend) then
            wtsz= (zt(kend)-zt(k))  &
                 /(zt(kend)-zt(kbeg))  &
                 *sigzwt
         endif

         uu(k,i,j)=(1.-wtsz)*uu(k,i,j)+wtsz*ps_u(i,j,k)
         vv(k,i,j)=(1.-wtsz)*vv(k,i,j)+wtsz*ps_v(i,j,k)
         tt(k,i,j)=(1.-wtsz)*tt(k,i,j)+wtsz*ps_t(i,j,k)
         rr(k,i,j)=(1.-wtsz)*rr(k,i,j)+wtsz*ps_r(i,j,k)

      enddo

   enddo
enddo

deallocate (v1,v2,v3,v4,v5,v6)

return
END SUBROUTINE isnsig

!##############################################################################
Subroutine visurf (n1,n2,n3,up,vp,thp,rtp,pp,topt,rtgt,zt)

use isan_coms
use rconstants

implicit none

integer :: n1,n2,n3
real, dimension(n1,n2,n3) :: up,vp,pp,thp,rtp
real, dimension(n2,n3) :: rtgt,topt
real, dimension(n1) :: zt

integer :: i,j,k 

real, allocatable :: v3(:)
real :: zgp,dzgobs,wt,ppp,ttt,rhm
real, external :: rsatmix


allocate(v3(n1))

do j=1,n3
   do i=1,n2
      if(rs_qual(i,j).gt..5) then

!          Quality of this surface analysis point is okay.

!          Extend the vertical influence of the surface analysis
!            up to a distance of SFCINF above the height of the
!            objectively analyzed surface.

         do k=1,n1
            zgp=topt(i,j)+zt(k)*rtgt(i,j)
            dzgobs=abs(zgp-rs_top(i,j))
            wt=dzgobs/sfcinf
            if(wt.le.1.) then
               if(rs_u(i,j).lt.1.e10)  &
                    up(k,i,j)=up(k,i,j)*wt+rs_u(i,j)*(1.-wt)
               if(rs_v(i,j).lt.1.e10)  &
                    vp(k,i,j)=vp(k,i,j)*wt+rs_v(i,j)*(1.-wt)
               ppp=(pp(k,i,j)/cp)**cpor*p00
               ttt=thp(k,i,j)*pp(k,i,j)/cp
               rhm=rtp(k,i,j)/rsatmix(ppp,ttt)
               if(rs_r(i,j).lt.1.e10)  &
                    rhm=rhm*wt+rs_r(i,j)*(1.-wt)
               if(rs_t(i,j).lt.1.e10)  &
                    thp(k,i,j)=thp(k,i,j)*wt+rs_t(i,j)*(1.-wt)
               ttt=thp(k,i,j)*pp(k,i,j)/cp
               rtp(k,i,j)=rhm*rsatmix(ppp,ttt)
            endif
         enddo

         do k=1,n1
            v3(k)=topt(i,j)+zt(k)*rtgt(i,j)
         enddo
         do k=n1-1,1,-1
            pp(k,i,j)=pp(k+1,i,j)+g*(v3(k+1)-v3(k))  &
                 /((thp(k,i,j)*(1.+.61*rtp(k,i,j))+thp(k+1,i,j)  &
                 *(1.+.61*rtp(k+1,i,j)))*.5)
         enddo

      endif

   enddo
enddo

deallocate(v3)

return
END SUBROUTINE visurf

!##############################################################################
Subroutine vshyd (n1,n2,n3,pp,tt,rr,topt,rtg,zt)
     
use isan_coms
use rconstants

implicit none
     
integer :: n1,n2,n3
real, dimension(n1,n2,n3) :: tt,rr,pp
real, dimension(n2,n3) :: topt,rtg
real, dimension(n1) :: zt

integer :: i,j,k,lbc,kabc
real :: thmin,ppp,ttt,tmpbc,rvibc,piibc,ziibc,gd2
real, external :: rsatmix

real, allocatable :: v3(:)

allocate(v3(n1))

!         Find internal boundary condition for the hydrostatic
!           equation

   thmin=tt(n1-1,1,1)
   do j=1,n3
      do i=1,n2
         thmin=min(thmin,tt(n1-1,i,j))
      enddo
   enddo

   do k=nisn,1,-1
      print*,'thmins-',k,levth(k),thmin
      if(float(levth(k)).le.thmin) then
         lbc=k
         go to 121
      endif
   enddo
   stop 'vi-nobc'
121     continue

   print 1212,lbc,levth(lbc)
1212    format(//' Hydrostatic boundary condition set at level'  &
        ,2I4,' K')

do j=1,n3
   do i=1,n2

!         Compute actual model heights

      do k=1,n1
         v3(k)=topt(i,j)+zt(k)*rtg(i,j)
      enddo

!         Find model level above boundary condition

         do k=1,n1
            if(tt(k,i,j).gt.levth(lbc))go to 141
         enddo
         stop 'vi-nbc2'
141      continue
         kabc=k

!         Integrate P to surface

         piibc=cp*(pi_p(i,j,lbc)*p00i)**rocp
         ziibc=(pi_s(i,j,lbc)-levth(lbc)*piibc)/g
         gd2=2.*g
         pp(kabc-1,i,j)=piibc+gd2*(ziibc-v3(kabc-1))  &
              /(levth(lbc)+tt(kabc-1,i,j))
         do k=kabc-2,1,-1
            pp(k,i,j)=pp(k+1,i,j)+gd2*(v3(k+1)-v3(k))  &
                 /(tt(k+1,i,j)+tt(k,i,j))
         enddo

!         Integrate P to model top

         do k=kabc,n1
            pp(k,i,j)=pp(k-1,i,j)-gd2*(v3(k)-v3(k-1))  &
                 /(tt(k,i,j)+tt(k-1,i,j))
         enddo
!
!         Compute mixing ratio from relative humidity and do final
!           hydrostatic integration.
!
         do k=1,n1
            ppp=(pp(k,i,j)/cp)**cpor*p00
            ttt=tt(k,i,j)*pp(k,i,j)/cp
            rr(k,i,j)=rr(k,i,j)*rsatmix(ppp,ttt)
         enddo

         tmpbc=levth(lbc)*piibc/cp
         rvibc=pi_r(i,j,lbc)*rsatmix(pi_p(i,j,lbc),tmpbc)

         pp(kabc-1,i,j)=piibc+gd2*(ziibc-v3(kabc-1))  &
              /(levth(lbc)*(1.+.61*rvibc)  &
               +tt(kabc-1,i,j)*(1.+.61*rr(kabc-1,i,j)))
         do k=kabc-2,1,-1
            pp(k,i,j)=pp(k+1,i,j)+gd2*(v3(k+1)-v3(k))  &
                 /(tt(k+1,i,j)*(1.+.61*rr(k+1,i,j))  &
                  +tt(k,i,j)*(1.+.61*rr(k,i,j)))
         enddo

         do k=kabc,n1
            pp(k,i,j)=pp(k-1,i,j)-gd2*(v3(k)-v3(k-1))  &
                 /((tt(k,i,j)*(1.+.61*rr(k,i,j))+tt(k-1,i,j)  &
                 *(1.+.61*rr(k-1,i,j))))
         enddo

   enddo
enddo

deallocate(v3)

return
END SUBROUTINE vshyd

!##############################################################################
Subroutine varuv (n1,n2,n3,up,vp)

implicit none

integer :: n1,n2,n3
real, dimension(n1,n2,n3) :: up,vp

integer :: i,j,k

!             Average to u,v points

do j=1,n3
   do i=1,n2-1
      do k=1,n1
         up(k,i,j)=(up(k,i,j)+up(k,i+1,j))*.5
      enddo
   enddo
enddo

do j=1,n3-1
   do i=1,n2
      do k=1,n1
         vp(k,i,j)=(vp(k,i,j)+vp(k,i,j+1))*.5
      enddo
   enddo
enddo

return
END SUBROUTINE varuv

!##############################################################################
Subroutine varfile_nstfeed (ifm,icm,n1f,n2f,n3f,n1c,n2c,n3c,nbot,ntop)

use isan_coms

implicit none

integer :: ifm,icm,n1f,n2f,n3f,n1c,n2c,n3c,nbot,ntop


!     Feed back the finer mesh to the coarser mesh.

CALL fdback (is_grids(icm)%rr_u   (1,1,1),is_grids(ifm)%rr_u   (1,1,1) &
           ,is_grids(icm)%rr_dn0u(1,1,1),is_grids(ifm)%rr_dn0u(1,1,1) &
           ,n1c,n2c,n3c,n1f,n2f,n3f,ifm,'u',rr_scr1(1))

CALL fdback (is_grids(icm)%rr_v   (1,1,1),is_grids(ifm)%rr_v   (1,1,1) &
           ,is_grids(icm)%rr_dn0v(1,1,1),is_grids(ifm)%rr_dn0v(1,1,1) &
           ,n1c,n2c,n3c,n1f,n2f,n3f,ifm,'v',rr_scr1(1))

CALL fdback (is_grids(icm)%rr_p  (1,1,1),is_grids(ifm)%rr_p  (1,1,1) &
           ,is_grids(icm)%rr_dn0(1,1,1),is_grids(ifm)%rr_dn0(1,1,1) &
           ,n1c,n2c,n3c,n1f,n2f,n3f,ifm,'p',rr_scr1(1))

CALL fdback (is_grids(icm)%rr_t  (1,1,1),is_grids(ifm)%rr_t  (1,1,1) &
           ,is_grids(icm)%rr_dn0(1,1,1),is_grids(ifm)%rr_dn0(1,1,1) &
           ,n1c,n2c,n3c,n1f,n2f,n3f,ifm,'t',rr_scr1(1))

CALL fdback (is_grids(icm)%rr_r  (1,1,1),is_grids(ifm)%rr_r  (1,1,1) &
           ,is_grids(icm)%rr_dn0(1,1,1),is_grids(ifm)%rr_dn0(1,1,1) &
           ,n1c,n2c,n3c,n1f,n2f,n3f,ifm,'t',rr_scr1(1))


if(nbot == 1) then
   CALL botset (n1c,n2c,n3c,is_grids(icm)%rr_u(1,1,1),'U')
   CALL botset (n1c,n2c,n3c,is_grids(icm)%rr_v(1,1,1),'V')
   CALL botset (n1c,n2c,n3c,is_grids(icm)%rr_p(1,1,1),'P')
   CALL botset (n1c,n2c,n3c,is_grids(icm)%rr_t(1,1,1),'T')
   CALL botset (n1c,n2c,n3c,is_grids(icm)%rr_r(1,1,1),'T')
endif

if(ntop == 1) then
   CALL topset (n1c,n2c,n3c,is_grids(icm)%rr_u(1,1,1),'U')
   CALL topset (n1c,n2c,n3c,is_grids(icm)%rr_v(1,1,1),'V')
   CALL topset (n1c,n2c,n3c,is_grids(icm)%rr_p(1,1,1),'P')
   CALL topset (n1c,n2c,n3c,is_grids(icm)%rr_t(1,1,1),'T')
   CALL topset (n1c,n2c,n3c,is_grids(icm)%rr_r(1,1,1),'T')
endif

return
END SUBROUTINE varfile_nstfeed
