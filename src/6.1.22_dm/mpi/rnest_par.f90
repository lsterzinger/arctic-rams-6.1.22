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

Subroutine par_bintp (ac,as,dn0f,n1m,n2m,n3m,m1f  &
     ,m1,m2,m3,ifm,ivnam,i0,j0,ibcon,bx,by,bz)

! When calling par_bintp, pass in the following values for its local
!   index limits ia,iz,ja,jz, which we are calling ia0,iz0,ja0,jz0
!   in the call:
!      ia0=iz
!      if(ivnam.eq.1)ia0=iz-1
!      iz0=1
!      ja0=jz
!      if(ivnam.eq.2)ja0=jz-1
!      jz0=1
!      if(iand(ibcon,2).ne.0)iz0=ia0
!      if(iand(ibcon,1).ne.0)ia0=1
!      if(iand(ibcon,8).ne.0)jz0=ja0
!      if(iand(ibcon,4).ne.0)ja0=1

use mem_grid
     
implicit none

integer :: n1m,n2m,n3m,m1f,m1,m2,m3,ifm,ivnam,i0,j0,ibcon
real :: ac(n1m,n2m,n3m),as(n1m,n2m,n3m),dn0f(m1,m2,m3)  &
   ,bx(m1,m3,2),by(m1,m2,2),bz(m2,m3,2)
integer :: ia,iz,ja,jz,nc,kc,ic,jc,kf,if,jf,k1,k2,im,jm

!_____________________________________________________________________
!
!    Set various indices. DO NOT CHANGE THE ORDER OF THESE STATEMENTS!!

ia=m2
if(ivnam.eq.1)ia=m2-1
iz=1

ja=m3
if(ivnam.eq.2)ja=m3-1
jz=1

!      if(iand(ibcon,2).ne.0)iz=ia
!      if(iand(ibcon,1).ne.0)ia=1
!      if(iand(ibcon,8).ne.0)jz=ja
!      if(iand(ibcon,4).ne.0)ja=1

if(ibcon.ne.0)then
   if(ibcon.ne.1)iz=ia
   if(ibcon.ne.2)ia=1
   if(ibcon.ne.4)jz=ja
   if(ibcon.ne.8)ja=1
endif
!_____________________________________________________________________

! i0 and j0 now passed in as fm node offsets
! ia, iz, ja, and jz refer to locations where the b arrays need to
!   be filled; they are not the usual limits of prognostic points.
! The nstbot=0 and nsttop=0 options do not work with the current ipaths
!   array definitions.

nc=nxtnest(ifm)
k1=max(1,kpm(2,ifm)-2)
k2=min(nnzp(nc),kpm(m1f-1,ifm)+2)

if(ivnam.eq.1)then
   do jc=jpm(ja+j0,ifm)-2,jpm(jz+j0,ifm)+1
      do if=ia,iz
         im=if+i0
         ic=ipm(im,ifm)
         do kc=k1,k2
            as(kc,if,jc)=ei4(im,ifm)*ac(kc,ic-2,jc)  &
                        +ei5(im,ifm)*ac(kc,ic-1,jc)  &
                        +ei6(im,ifm)*ac(kc,ic  ,jc)  &
                        +ei7(im,ifm)*ac(kc,ic+1,jc)
         enddo
      enddo
   enddo
else
   do jc=jpm(ja+j0,ifm)-2,jpm(jz+j0,ifm)+1
      do if=ia,iz
         im=if+i0
         ic=ipm(im,ifm)
         do kc=k1,k2
            as(kc,if,jc)=ei1(im,ifm)*ac(kc,ic-1,jc)  &
                        +ei2(im,ifm)*ac(kc,ic  ,jc)  &
                        +ei3(im,ifm)*ac(kc,ic+1,jc)
         enddo
      enddo
   enddo
endif
!
if(ivnam.eq.2)then
   do jf=ja,jz
      jm=jf+j0
      jc=jpm(jm,ifm)
      do if=ia,iz
         do kc=k1,k2
            ac(kc,if,jf)=ej4(jm,ifm)*as(kc,if,jc-2)  &
                        +ej5(jm,ifm)*as(kc,if,jc-1)  &
                        +ej6(jm,ifm)*as(kc,if,jc  )  &
                        +ej7(jm,ifm)*as(kc,if,jc+1)
         enddo
      enddo
   enddo
else
   do jf=ja,jz
      jm=jf+j0
      jc=jpm(jm,ifm)
      do if=ia,iz
         do kc=k1,k2
            ac(kc,if,jf)=ej1(jm,ifm)*as(kc,if,jc-1)  &
                        +ej2(jm,ifm)*as(kc,if,jc  )  &
                        +ej3(jm,ifm)*as(kc,if,jc+1)
         enddo
      enddo
   enddo
endif
!
if(ivnam.eq.3)then

   if(iand(ibcon,1).ne.0) then
      if=ia
      do jf=ja,jz
         do kf=1,m1f-1
            kc=kpm(kf+1,ifm)
            bx(kf,jf,1)=  &
                (ek4(kf,ifm)*ac(max(1         ,kc-2),if,jf)  &
               + ek5(kf,ifm)*ac(               kc-1 ,if,jf)  &
               + ek6(kf,ifm)*ac(               kc   ,if,jf)  &
               + ek7(kf,ifm)*ac(min(nnzp(nc)-1,kc+1),if,jf))  &
               / (.5 * (dn0f(kf,if,jf) + dn0f(kf+1,if,jf)))
         enddo
      enddo
   endif

   if(iand(ibcon,2).ne.0) then
      if=iz
      do jf=ja,jz
         do kf=1,m1f-1
            kc=kpm(kf+1,ifm)
            bx(kf,jf,2)=  &
                (ek4(kf,ifm)*ac(max(1         ,kc-2),if,jf)  &
               + ek5(kf,ifm)*ac(               kc-1 ,if,jf)  &
               + ek6(kf,ifm)*ac(               kc   ,if,jf)  &
               + ek7(kf,ifm)*ac(min(nnzp(nc)-1,kc+1),if,jf))  &
               / (.5 * (dn0f(kf,if,jf) + dn0f(kf+1,if,jf)))
         enddo
      enddo
   endif

   if(iand(ibcon,4).ne.0) then
      jf=ja
      do if=ia,iz
         do kf=1,m1f-1
            kc=kpm(kf+1,ifm)
            by(kf,if,1)=  &
                (ek4(kf,ifm)*ac(max(1         ,kc-2),if,jf)  &
               + ek5(kf,ifm)*ac(               kc-1 ,if,jf)  &
               + ek6(kf,ifm)*ac(               kc   ,if,jf)  &
               + ek7(kf,ifm)*ac(min(nnzp(nc)-1,kc+1),if,jf))  &
               / (.5 * (dn0f(kf,if,jf) + dn0f(kf+1,if,jf)))
         enddo
      enddo
   endif

   if(iand(ibcon,8).ne.0) then
      jf=jz
      do if=ia,iz
         do kf=1,m1f-1
            kc=kpm(kf+1,ifm)
            by(kf,if,2)=  &
                (ek4(kf,ifm)*ac(max(1         ,kc-2),if,jf)  &
               + ek5(kf,ifm)*ac(               kc-1 ,if,jf)  &
               + ek6(kf,ifm)*ac(               kc   ,if,jf)  &
               + ek7(kf,ifm)*ac(min(nnzp(nc)-1,kc+1),if,jf))  &
               / (.5 * (dn0f(kf,if,jf) + dn0f(kf+1,if,jf)))
         enddo
      enddo
   endif

   if(nstbot.eq.0)then
      kf=1
      do jf=ja,jz
         do if=ia,iz
            kc=kpm(kf+1,ifm)
            bz(if,jf,1)=  &
                (ek4(kf,ifm)*ac(max(1         ,kc-2),if,jf)  &
               + ek5(kf,ifm)*ac(               kc-1 ,if,jf)  &
               + ek6(kf,ifm)*ac(               kc   ,if,jf)  &
               + ek7(kf,ifm)*ac(min(nnzp(nc)-1,kc+1),if,jf))  &
               / (.5 * (dn0f(kf,if,jf) + dn0f(kf+1,if,jf)))
         enddo
      enddo
   endif

   if(nsttop.eq.0)then
      kf=m1f-1
      do jf=ja,jz
         do if=ia,iz
            kc=kpm(kf+1,ifm)
            bz(if,jf,2)=  &
                (ek4(kf,ifm)*ac(max(1         ,kc-2),if,jf)  &
               + ek5(kf,ifm)*ac(               kc-1 ,if,jf)  &
               + ek6(kf,ifm)*ac(               kc   ,if,jf)  &
               + ek7(kf,ifm)*ac(min(nnzp(nc)-1,kc+1),if,jf))  &
               / (.5 * (dn0f(kf,if,jf) + dn0f(kf+1,if,jf)))
         enddo
      enddo
   endif

elseif(ivnam.eq.4)then

   if(iand(ibcon,1).ne.0) then
      if=ia
      do jf=ja,jz
         do kf=1,m1f
            bx(kf,jf,1)=  &
                (ek1(kf,ifm)*ac(kpm(kf,ifm)-1,if,jf)  &
               + ek2(kf,ifm)*ac(kpm(kf,ifm)  ,if,jf)  &
               + ek3(kf,ifm)*ac(kpm(kf,ifm)+1,if,jf))
         enddo
      enddo
   endif

   if(iand(ibcon,2).ne.0) then
      if=iz
      do jf=ja,jz
         do kf=1,m1f
            bx(kf,jf,2)=  &
                (ek1(kf,ifm)*ac(kpm(kf,ifm)-1,if,jf)  &
               + ek2(kf,ifm)*ac(kpm(kf,ifm)  ,if,jf)  &
               + ek3(kf,ifm)*ac(kpm(kf,ifm)+1,if,jf))
         enddo
      enddo
   endif

   if(iand(ibcon,4).ne.0) then
      jf=ja
      do if=ia,iz
         do kf=1,m1f
            by(kf,if,1)=  &
                (ek1(kf,ifm)*ac(kpm(kf,ifm)-1,if,jf)  &
               + ek2(kf,ifm)*ac(kpm(kf,ifm)  ,if,jf)  &
               + ek3(kf,ifm)*ac(kpm(kf,ifm)+1,if,jf))
         enddo
      enddo
   endif

   if(iand(ibcon,8).ne.0) then
      jf=jz
      do if=ia,iz
         do kf=1,m1f
            by(kf,if,2)=  &
                (ek1(kf,ifm)*ac(kpm(kf,ifm)-1,if,jf)  &
               + ek2(kf,ifm)*ac(kpm(kf,ifm)  ,if,jf)  &
               + ek3(kf,ifm)*ac(kpm(kf,ifm)+1,if,jf))
         enddo
      enddo
   endif

   if(nstbot.eq.0)then
      kf=1
      do jf=ja,jz
         do if=ia,iz
            bz(if,jf,1)=  &
                (ek1(kf,ifm)*ac(kpm(kf,ifm)-1,if,jf)  &
               + ek2(kf,ifm)*ac(kpm(kf,ifm)  ,if,jf)  &
               + ek3(kf,ifm)*ac(kpm(kf,ifm)+1,if,jf))
         enddo
      enddo
   endif

   if(nsttop.eq.0)then
      kf=m1f
      do jf=ja,jz
         do if=ia,iz
            bz(if,jf,2)=  &
                (ek1(kf,ifm)*ac(kpm(kf,ifm)-1,if,jf)  &
               + ek2(kf,ifm)*ac(kpm(kf,ifm)  ,if,jf)  &
               + ek3(kf,ifm)*ac(kpm(kf,ifm)+1,if,jf))
         enddo
      enddo
   endif

else

   if(iand(ibcon,1).ne.0) then
      if=ia
      do jf=ja,jz
         do kf=1,m1f
            bx(kf,jf,1)=  &
                (ek1(kf,ifm)*ac(kpm(kf,ifm)-1,if,jf)  &
               + ek2(kf,ifm)*ac(kpm(kf,ifm)  ,if,jf)  &
               + ek3(kf,ifm)*ac(kpm(kf,ifm)+1,if,jf))  &
               / dn0f(kf,if,jf)
         enddo
      enddo
   endif

   if(iand(ibcon,2).ne.0) then
      if=iz
      do jf=ja,jz
         do kf=1,m1f
            bx(kf,jf,2)=  &
                (ek1(kf,ifm)*ac(kpm(kf,ifm)-1,if,jf)  &
               + ek2(kf,ifm)*ac(kpm(kf,ifm)  ,if,jf)  &
               + ek3(kf,ifm)*ac(kpm(kf,ifm)+1,if,jf))  &
               / dn0f(kf,if,jf)
         enddo
      enddo
   endif

   if(iand(ibcon,4).ne.0) then
      jf=ja
      do if=ia,iz
         do kf=1,m1f
            by(kf,if,1)=  &
                (ek1(kf,ifm)*ac(kpm(kf,ifm)-1,if,jf)  &
               + ek2(kf,ifm)*ac(kpm(kf,ifm)  ,if,jf)  &
               + ek3(kf,ifm)*ac(kpm(kf,ifm)+1,if,jf))  &
               / dn0f(kf,if,jf)
         enddo
      enddo
   endif

   if(iand(ibcon,8).ne.0) then
      jf=jz
      do if=ia,iz
         do kf=1,m1f
            by(kf,if,2)=  &
                (ek1(kf,ifm)*ac(kpm(kf,ifm)-1,if,jf)  &
               + ek2(kf,ifm)*ac(kpm(kf,ifm)  ,if,jf)  &
               + ek3(kf,ifm)*ac(kpm(kf,ifm)+1,if,jf))  &
               / dn0f(kf,if,jf)
         enddo
      enddo
   endif

   if(nstbot.eq.0)then
      kf=1
      do jf=ja,jz
         do if=ia,iz
            bz(if,jf,1)=  &
                (ek1(kf,ifm)*ac(kpm(kf,ifm)-1,if,jf)  &
               + ek2(kf,ifm)*ac(kpm(kf,ifm)  ,if,jf)  &
               + ek3(kf,ifm)*ac(kpm(kf,ifm)+1,if,jf))  &
               / dn0f(kf,if,jf)
         enddo
      enddo
   endif

   if(nsttop.eq.0)then
      kf=m1f
      do jf=ja,jz
         do if=ia,iz
            bz(if,jf,2)=  &
                (ek1(kf,ifm)*ac(kpm(kf,ifm)-1,if,jf)  &
               + ek2(kf,ifm)*ac(kpm(kf,ifm)  ,if,jf)  &
               + ek3(kf,ifm)*ac(kpm(kf,ifm)+1,if,jf))  &
               / dn0f(kf,if,jf)
         enddo
      enddo
   endif

endif

return
END SUBROUTINE par_bintp



