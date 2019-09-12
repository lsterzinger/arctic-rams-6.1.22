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

Subroutine datassim ()

use mem_tend
use mem_basic
use mem_grid
use mem_varinit
use mem_scratch
use node_mod

implicit none

integer :: il,ir,jl,jr

!     Set bounds for nudging this sub-domain

il=ia  ;  ir=iz  ;  jl=ja  ;  jr=jz

!     West and east boundaries.
if(iand(ibcon,1) /= 0)  il=1
if(iand(ibcon,2) /= 0)  ir=mxp


!     South and north boundaries.
if(jdim == 1) then
   if(iand(ibcon,4) /= 0) jl=1
   if(iand(ibcon,8) /= 0) jr=myp
endif


! Basic boundary and analysis nudging scheme

 CALL nudge (mzp,mxp,myp,il,ir,jl,jr,varinit_g(ngrid)%varwts(1,1,1)  &

     ,varinit_g(ngrid)%varup(1,1,1),varinit_g(ngrid)%varvp(1,1,1)  &
     ,varinit_g(ngrid)%varpp(1,1,1),varinit_g(ngrid)%vartp(1,1,1)  &
     ,varinit_g(ngrid)%varrp(1,1,1) &
     
     ,varinit_g(ngrid)%varuf(1,1,1),varinit_g(ngrid)%varvf(1,1,1)  &
     ,varinit_g(ngrid)%varpf(1,1,1),varinit_g(ngrid)%vartf(1,1,1)  &
     ,varinit_g(ngrid)%varrf(1,1,1)  &
     
     ,basic_g(ngrid)%up(1,1,1)   ,basic_g(ngrid)%vp(1,1,1)  &
     ,basic_g(ngrid)%theta(1,1,1),basic_g(ngrid)%rtp(1,1,1)  &
     ,basic_g(ngrid)%pp(1,1,1)  &
     ,tend%ut(1),tend%vt(1),tend%tht(1),tend%rtt(1),tend%pt(1))


! Condensate nudging scheme

if (nud_cond == 1 .and. time >= tcond_beg .and. time <= tcond_end) &
 CALL nudge_cond (mzp,mxp,myp,il,ir,jl,jr,varinit_g(ngrid)%varwts(1,1,1)  &
        ,varinit_g(ngrid)%varrph(1,1,1),varinit_g(ngrid)%varcph(1,1,1) &
        ,varinit_g(ngrid)%varrfh(1,1,1),varinit_g(ngrid)%varcfh(1,1,1)  &
        ,basic_g(ngrid)%rtp(1,1,1),tend%rtt(1))

return
END SUBROUTINE datassim

!##############################################################################
Subroutine nudge (m1,m2,m3,ia,iz,ja,jz,varwts  &
     ,varup,varvp,varpp,vartp,varrp  &
     ,varuf,varvf,varpf,vartf,varrf  &
     ,up,vp,theta,rtp,pp,ut,vt,tht,rtt,pt)

use mem_grid
use mem_varinit
use mem_scratch

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz
real, dimension(m1,m2,m3) :: varup,varvp,vartp,varrp,varpp  &
                            ,varuf,varvf,vartf,varrf,varpf  &
                            ,varwts,up,vp,theta,rtp,pp  &
                            ,ut,vt,tht,rtt,pt             

integer :: i,j,k
real :: tfact, wt_uv, wt_th, wt_pi, wt_rt

!         Linearly interpolate values in time, then nudge.

if (nud_type == 1) then
   tfact=(time-vtime1)/(vtime2-vtime1)
else
   stop 'should not be nudging if nud_type=0'
endif

wt_uv=wt_nudge_g(ngrid)* wt_nudge_uv
wt_th=wt_nudge_g(ngrid)* wt_nudge_th
wt_pi=wt_nudge_g(ngrid)* wt_nudge_pi
wt_rt=wt_nudge_g(ngrid)* wt_nudge_rt

do j=ja,jz
   do i=ia,iz

      do k=1,m1
         vctr1(k)=varup(k,i,j)+(varuf(k,i,j)-varup(k,i,j))*tfact
         vctr2(k)=varvp(k,i,j)+(varvf(k,i,j)-varvp(k,i,j))*tfact
         vctr3(k)=vartp(k,i,j)+(vartf(k,i,j)-vartp(k,i,j))*tfact
         vctr4(k)=varpp(k,i,j)+(varpf(k,i,j)-varpp(k,i,j))*tfact
         vctr5(k)=varrp(k,i,j)+(varrf(k,i,j)-varrp(k,i,j))*tfact

         vctr10(k)=(varwts(k,i,j)+varwts(k,min(m2,i+1),j))*.5* wt_uv
         vctr11(k)=(varwts(k,i,j)+varwts(k,i,min(m3,j+jdim)))*.5* wt_uv
         vctr12(k)=varwts(k,i,j)* wt_th
         vctr13(k)=varwts(k,i,j)* wt_pi
         vctr14(k)=varwts(k,i,j)* wt_rt
      enddo

      do k=1,m1
         ut(k,i,j) = ut(k,i,j) + vctr10(k)*(vctr1(k)-up(k,i,j))
         vt(k,i,j) = vt(k,i,j) + vctr11(k)*(vctr2(k)-vp(k,i,j))
         tht(k,i,j)=tht(k,i,j) + vctr12(k)*(vctr3(k)-theta(k,i,j))
         pt(k,i,j) = pt(k,i,j) + vctr13(k)*(vctr4(k)-pp(k,i,j))
         rtt(k,i,j)=rtt(k,i,j) + vctr14(k)*(vctr5(k)-rtp(k,i,j))
      enddo

   enddo
enddo

return
END SUBROUTINE nudge

!##############################################################################
Subroutine nudge_cond (m1,m2,m3,ia,iz,ja,jz,varwts  &
                     ,varrph,varcph,varrfh,varcfh,rtp,rtt)

use mem_grid
use mem_varinit
use mem_scratch

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz
real, dimension(m1,m2,m3) :: varrph,varcph,varrfh,varcfh  &
                            ,varwts,rtp,rtt             

integer :: i,j,k
real :: tfact, wt_rc

tfact=(time-vtime1)/(vtime2-vtime1)

wt_rc=wt_nudgec(ngrid)/t_nudge_rc

do j=ja,jz
   do i=ia,iz

      do k=1,m1
         vctr1(k)=varrph(k,i,j)+(varrfh(k,i,j)-varrph(k,i,j))*tfact
         vctr2(k)=varcph(k,i,j)+(varcfh(k,i,j)-varcph(k,i,j))*tfact

         !Choose a nudging weight method below. By default we are only
         !applying the condensate nudging weight, but to the whole domain

         !This method sets nudging weights in conjunction with standard
         !analysis nudging weights in space and time
         vctr3(k)=varwts(k,i,j)*wt_rc

         !This method applies only condensate nudging weights only for 
         !all the domain. Currently use this to override option just above.
         vctr3(k)=wt_rc
      enddo

      ! Only nudging total water RTP where condensate mixing ratio exists.
      ! Condensate will have a value of -9999 for traditional varfiles using
      ! "MAKEVFILE" since we do not ingest condensate from "dp" files. We do 
      ! ingest condensate from History-Varfiles (MAKEHFILE), so values should
      ! be positive-definite.
      do k=1,m1
         if (vctr2(k) > 0.)  &
            rtt(k,i,j)=rtt(k,i,j) + vctr3(k)*(vctr1(k)-rtp(k,i,j))
      enddo

   enddo
enddo

return
END SUBROUTINE nudge_cond

!##############################################################################
Subroutine varweight (n1,n2,n3,varwts,topt,rtgt)

use mem_grid
use mem_varinit

implicit none

integer :: n1,n2,n3
real :: varwts(n1,n2,n3),topt(n2,n3),rtgt(n2,n3)

integer :: i,j,k
real :: tnudcenti,tnudtopi,tnudlati,rown,rows,rowe,roww,zloc,wttop &
       ,wtlat,delzi

!         Get weights for large scale and model tendencies

if (nudlat .le. 0) return

tnudcenti=0.
if(tnudcent.gt. .01) tnudcenti=1./tnudcent
tnudtopi=0.
if(tnudtop.gt. .01) tnudtopi=1./tnudtop-tnudcenti
tnudlati=0.
if(tnudlat.gt. .01) tnudlati=1./tnudlat-tnudcenti

delzi=0.0 !Variable initialized
if(ztop.gt.znudtop) then
   delzi=1./(ztop-znudtop)
elseif(tnudtop.gt. .01) then
   print*,'Incorrect specification of znudtop ! ! !'
   print*,' znudtop = ',znudtop
   print*,'    ztop = ',ztop
   stop 'varwt-znud'
endif


do j=1,n3
   do i=1,n2

!                       quadratic weight func for lateral boundaries

      if (n3 > 1) then
         rown=max(0.,float(j+nudlat-n3))
         rows=max(0.,float(nudlat+1-j))
      else
         rown=0.
         rows=0.
      endif

      rowe=max(0.,float(i+nudlat-n2))
      roww=max(0.,float(nudlat+1-i))
      wtlat=max(rown*rown,rows*rows,rowe*rowe,roww*roww)  &
           /float(nudlat*nudlat)

!                       linear weight func for top boundary

      do k=1,n1
         zloc=ztn(k,1)*rtgt(i,j)+topt(i,j)
         wttop=max(0.,(zloc-znudtop)*delzi)

!                       full 3-D weight func

         varwts(k,i,j)=tnudcenti  &
              +max(tnudlati*wtlat,tnudtopi*wttop)

      enddo

   enddo
enddo

return
END SUBROUTINE varweight

!##############################################################################
Subroutine vfintrpf (ifm,ifflag)

use mem_grid
use mem_scratch
use mem_varinit
use mem_basic

implicit none

integer :: ifm,ifflag,icm

icm = nxtnest(ifm)
if (icm == 0) return

!    Temporarily fill VT2DA with interpolated topography from coarser grid

CALL fillscr (1,maxnxp,maxnyp,1,nnxp(icm),nnyp(icm),1,1  &
   ,scratch%scr1(1),grid_g(icm)%topt(1,1))
CALL eintp (scratch%scr1(1),scratch%scr2(1)  &
   ,1,maxnxp,maxnyp,1,nnxp(ifm),nnyp(ifm),ifm,2,'t',0,0)
CALL fillvar (1,maxnxp,maxnyp,1,nnxp(ifm),nnyp(ifm),1,1  &
   ,scratch%scr2(1),scratch%vt2da(1))

if (ifflag == 1) then

!     Interpolate varwts

   CALL fmint4 (varinit_g(icm)%varwts(1,1,1)  &
              ,varinit_g(ifm)%varwts(1,1,1)  &
              ,basic_g(icm)%dn0(1,1,1),basic_g(ifm)%dn0(1,1,1)  &
              ,scratch%vt2da(1),ifm,icm,'t',0)

endif


if (ifflag == 2) then

!     Interpolate future level atmospheric variables

   CALL fmint4 (varinit_g(icm)%varuf(1,1,1),varinit_g(ifm)%varuf(1,1,1)  &
      ,basic_g(icm)%dn0u(1,1,1),basic_g(ifm)%dn0u(1,1,1)  &
      ,scratch%vt2da(1),ifm,icm,'u',1)
   CALL fmint4 (varinit_g(icm)%varvf(1,1,1),varinit_g(ifm)%varvf(1,1,1)  &
      ,basic_g(icm)%dn0v(1,1,1),basic_g(ifm)%dn0v(1,1,1)  &
      ,scratch%vt2da(1),ifm,icm,'v',1)
   CALL fmint4 (varinit_g(icm)%varpf(1,1,1),varinit_g(ifm)%varpf(1,1,1)  &
      ,basic_g(icm)%dn0v(1,1,1),basic_g(ifm)%dn0v(1,1,1)  &
      ,scratch%vt2da(1),ifm,icm,'t',1)
   CALL fmint4 (varinit_g(icm)%vartf(1,1,1),varinit_g(ifm)%vartf(1,1,1)  &
      ,basic_g(icm)%dn0(1,1,1),basic_g(ifm)%dn0(1,1,1)  &
      ,scratch%vt2da(1),ifm,icm,'t',1)
   CALL fmint4 (varinit_g(icm)%varrf(1,1,1),varinit_g(ifm)%varrf(1,1,1)  &
      ,basic_g(icm)%dn0(1,1,1),basic_g(ifm)%dn0(1,1,1)  &
      ,scratch%vt2da(1),ifm,icm,'t',1)

endif

return
END SUBROUTINE vfintrpf


