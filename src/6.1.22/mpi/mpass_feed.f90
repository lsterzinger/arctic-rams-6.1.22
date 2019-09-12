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

Subroutine node_sendfeed (ngr)

use mem_grid
use node_mod
use mem_basic
use var_tables

implicit none

integer :: ngr
integer :: i1s,i2s,j1s,j2s,k1s,k2s,mtp,i1f,i2f,j1f,j2f
integer :: nm,icm,ifm,itype,itypef,nv,iptr,nvar
real, save, allocatable::pbuff(:)
integer, save :: nbuff_save=0

icm=nxtnest(ngr)
ifm=ngr

itype=6

!______________________
!
!   First, before we send anything, let's post the receives.

do nm=1,nmachs
   irecv_req(nm)=0
   if (iget_paths(itype,ifm,nm).ne.0) then
      CALL par_get_noblock (node_buffs(nm)%lbc_recv_buff(1)  &
          ,node_buffs(nm)%nrecv ,5500+icm,machs(nm),irecv_req(nm) )
          !  if(mynum==10) print*,mynum,'posted FEED receive',nm,irecv_req(nm)  &
          !        ,node_buffs(nm)%nrecv
   endif
enddo
!______________________

!     Allocate new temporary buffer if bigger than the old one.

if(nbuff_feed > nbuff_save) then
  if(iprntstmt>=1 .and. print_msg) &
     print*,'Allocating feed send buffer:',mynum,nbuff_feed,nbuff_save
  allocate (pbuff(nbuff_feed))
  nbuff_save=nbuff_feed
endif


!     Construct table of variables being fed back.

!ivar_feed(1)=2
!ivar_feed(2)=4
!ivar_feed(3)=6
!ivar_feed(4)=8

!ivar_type(1)=1
!ivar_type(2)=2
!ivar_type(3)=3
!ivar_type(4)=4
!nvar=4 + num_scalar(ifm)


!     Feed back this fine grid's portion of the each coarse grid node

k1s=kpm(2,ifm)
k2s=kpm(nnzp(ifm)-1,ifm)

do nm=1,nmachs
   isend_req(nm)=0
   if(ipaths(1,itype,ifm,nm).ne.0) then

!            mtp=total number of coarse grid points created
!              and sent from this fine grid node

      i1s=ipaths(1,itype,ifm,nm)
      i2s=ipaths(2,itype,ifm,nm)
      j1s=ipaths(3,itype,ifm,nm)
      j2s=ipaths(4,itype,ifm,nm)
      mtp=(i2s-i1s+1)*(j2s-j1s+1)*(k2s-k1s+1)

      itypef=7
      i1f=ipaths(1,itypef,ifm,nm)
      i2f=ipaths(2,itypef,ifm,nm)
      j1f=ipaths(3,itypef,ifm,nm)
      j2f=ipaths(4,itypef,ifm,nm)

      iptr=0
      CALL fdbackp (1,basic_g(ifm)%uc(1,1,1),pbuff(1+iptr),mtp  &
          ,basic_g(ifm)%dn0(1,1,1),basic_g(ifm)%dn0u(1,1,1)  &
          ,basic_g(ifm)%dn0v(1,1,1)  &
          ,mmzp(ifm),mmxp(ifm),mmyp(ifm)  &
          ,ifm,i1f-i0,i2f-i0,j1f-j0,j2f-j0  &
          ,i0,j0,nstratx(ifm),nstraty(ifm))
      iptr=iptr+mtp
      CALL fdbackp (2,basic_g(ifm)%vc(1,1,1),pbuff(1+iptr),mtp  &
          ,basic_g(ifm)%dn0(1,1,1),basic_g(ifm)%dn0u(1,1,1)  &
          ,basic_g(ifm)%dn0v(1,1,1)  &
          ,mmzp(ifm),mmxp(ifm),mmyp(ifm)  &
          ,ifm,i1f-i0,i2f-i0,j1f-j0,j2f-j0  &
          ,i0,j0,nstratx(ifm),nstraty(ifm))
      iptr=iptr+mtp
      CALL fdbackp (3,basic_g(ifm)%wc(1,1,1),pbuff(1+iptr),mtp  &
          ,basic_g(ifm)%dn0(1,1,1),basic_g(ifm)%dn0u(1,1,1)  &
          ,basic_g(ifm)%dn0v(1,1,1)  &
          ,mmzp(ifm),mmxp(ifm),mmyp(ifm)  &
          ,ifm,i1f-i0,i2f-i0,j1f-j0,j2f-j0  &
          ,i0,j0,nstratx(ifm),nstraty(ifm))
      iptr=iptr+mtp
      CALL fdbackp (4,basic_g(ifm)%pc(1,1,1),pbuff(1+iptr),mtp  &
          ,basic_g(ifm)%dn0(1,1,1),basic_g(ifm)%dn0u(1,1,1)  &
          ,basic_g(ifm)%dn0v(1,1,1)  &
          ,mmzp(ifm),mmxp(ifm),mmyp(ifm)  &
          ,ifm,i1f-i0,i2f-i0,j1f-j0,j2f-j0  &
          ,i0,j0,nstratx(ifm),nstraty(ifm))
      iptr=iptr+mtp

      do nv=1,num_scalar(ifm)
         CALL fdbackp (5,scalar_tab(nv,ifm)%var_p,pbuff(1+iptr),mtp  &
             ,basic_g(ifm)%dn0(1,1,1),basic_g(ifm)%dn0u(1,1,1)  &
             ,basic_g(ifm)%dn0v(1,1,1)  &
             ,mmzp(ifm),mmxp(ifm),mmyp(ifm)  &
             ,ifm,i1f-i0,i2f-i0,j1f-j0,j2f-j0  &
             ,i0,j0,nstratx(ifm),nstraty(ifm))
         iptr=iptr+mtp
      enddo

!     We will send master coarse grid indices to nodes.
      CALL par_init_put (node_buffs(nm)%lbc_send_buff(1)  &
                       ,node_buffs(nm)%nsend )
         CALL par_put_int (i1s,1)
         CALL par_put_int (i2s,1)
         CALL par_put_int (j1s,1)
         CALL par_put_int (j2s,1)
         CALL par_put_int (k1s,1)
         CALL par_put_int (k2s,1)
         CALL par_put_int (mynum,1)
         CALL par_put_int (nvar,1)
         CALL par_put_int (iptr,1)
         CALL par_put_float (pbuff(1),iptr)
      CALL par_send_noblock (ipaths(5,itype,ifm,nm),5500+icm  &
              ,isend_req(nm) )

   endif
enddo

return
END SUBROUTINE node_sendfeed

!##############################################################################
Subroutine node_getfeed (icm,ifm)

use mem_grid
use node_mod
use var_tables
use mem_basic

implicit none

integer :: icm,ifm
integer :: i1s,i2s,j1s,j2s,k1s,k2s,mtp,icall,i1z,i2z,j1z,j2z,k1z,k2z
integer :: i2zu,j2zv,k2zw,i2u,j2v,k2w,nfx,nfy,nfz
integer :: nm,itype,nv,iptr,nvar,ibytes,msgid,ihostnum,machf,nwds
real, save, allocatable::pbuff(:)
integer, save :: nbuff_save=0

itype=6

!_____________________________________________________________________
!
!  First, let's make sure our sends are all finished and de-allocated

do nm=1,nmachs
   if (ipaths(1,itype,ifm,nm).ne.0 ) then
      CALL par_wait (isend_req(nm),ibytes,msgid,ihostnum)
   endif
enddo
!      print*,mynum,'done FEED send wait'
!_____________________________________________________________________
!
!  Now, let's wait on our receives

do nm=1,nmachs
   if (iget_paths(itype,ifm,nm).ne.0) then
      CALL par_wait (irecv_req(nm),ibytes,msgid,ihostnum)
   endif
enddo
!      print*,mynum,'done FEED recv wait'
!_____________________________________________________________________

!
!     Can we use existing memory for the buffers? If not, allocate new
!       temporary buffer.

!     Allocate new temporary buffer if bigger than the old one.

!if(node_buffs(mynum)%nrecv > nbuff_save) then
!   print*,'Allocating feed recv buffer:',mynum,nbuff_feed,nbuff_save  &
!          ,node_buffs(mynum)%nrecv
!   if (allocated(pbuff)) deallocate (pbuff)
!  allocate (pbuff(node_buffs(mynum)%nrecv))
!  nbuff_save=node_buffs(mynum)%nrecv
!endif

icall=1
do nm=1,nmachs
   if(iget_paths(itype,ifm,nm).ne.0) then

      CALL par_assoc_buff (node_buffs(nm)%lbc_recv_buff(1)  &
                         ,node_buffs(nm)%nrecv)

      CALL par_get_int (i1s,1)
      CALL par_get_int (i2s,1)
      CALL par_get_int (j1s,1)
      CALL par_get_int (j2s,1)
      CALL par_get_int (k1s,1)
      CALL par_get_int (k2s,1)
      CALL par_get_int (machf,1)
      CALL par_get_int (nvar,1)
      CALL par_get_int (nwds,1)
      
      ! Make sure buffer for floating point info big enough
      
      if (nwds >  nbuff_save) then
         if (allocated(pbuff)) deallocate (pbuff)
         allocate (pbuff(nwds))
         nbuff_save=nwds
      endif      
      
      CALL par_get_float (pbuff,nwds)

      if(icall==1) then
         icall=0

! Set the portion of this coarse grid subdomain that will be filled
! with fine grid info to zero.

         i1z = max(ipm(2,ifm),1+mi0(icm))
         j1z = max(jpm(2,ifm),1+mj0(icm))
         k1z = kpm(2,ifm)

         i2z = min(ipm(nnxp(ifm)-1,ifm),mmxp(icm)+mi0(icm))
         j2z = min(jpm(nnyp(ifm)-1,ifm),mmyp(icm)+mj0(icm))
         k2z = kpm(nnzp(ifm)-1,ifm)

         i2zu = min(ipm(nnxp(ifm)-1-nstratx(ifm),ifm)  &
                   ,mmxp(icm)+mi0(icm))
         CALL zeroout (basic_g(icm)%uc(1,1,1),mmzp(icm),mmxp(icm),mmyp(icm) &
                 ,mi0(icm),mj0(icm),i1z,i2zu,j1z,j2z,k1z,k2z)
         j2zv = min(jpm(nnyp(ifm)-1-nstraty(ifm),ifm)  &
                   ,mmyp(icm)+mj0(icm))
         CALL zeroout (basic_g(icm)%vc(1,1,1),mmzp(icm),mmxp(icm),mmyp(icm) &
                 ,mi0(icm),mj0(icm),i1z,i2z,j1z,j2zv,k1z,k2z)
         k2zw = kpm(nnzp(ifm)-1-nrz(kpm(nnzp(ifm)-1,ifm)  &
                   ,ifm),ifm)
         CALL zeroout (basic_g(icm)%wc(1,1,1),mmzp(icm),mmxp(icm),mmyp(icm) &
                 ,mi0(icm),mj0(icm),i1z,i2z,j1z,j2z,k1z,k2zw)

         CALL zeroout (basic_g(icm)%pc(1,1,1),mmzp(icm),mmxp(icm),mmyp(icm) &
                 ,mi0(icm),mj0(icm),i1z,i2z,j1z,j2z,k1z,k2z)

         do nv=1,num_scalar(ifm)
            CALL zeroout (scalar_tab(nv,icm)%var_p  &
                 ,mmzp(icm),mmxp(icm),mmyp(icm)  &
                 ,mi0(icm),mj0(icm),i1z,i2z,j1z,j2z,k1z,k2z)
         enddo

      endif

      i2u = min(i2s,i2zu)
      j2v = min(j2s,j2zv)
      k2w = min(k2s,k2zw)

      nfx=i2s-i1s+1
      nfy=j2s-j1s+1
      nfz=k2s-k1s+1
      mtp=nfx*nfy*nfz

      iptr=0

      CALL unfdbackp (1,basic_g(icm)%uc(1,1,1),pbuff(1+iptr)  &
          ,basic_g(icm)%dn0(1,1,1),basic_g(icm)%dn0u(1,1,1) &
          ,basic_g(icm)%dn0v(1,1,1)  &
          ,mmzp(icm),mmxp(icm),mmyp(icm),mi0(icm),mj0(icm)  &
          ,i1s,i2s,j1s,j2s  &
          ,k1s,k2s,i2u,j2v,k2w,nfz,nfx,nfy)
      iptr=iptr+mtp
      CALL unfdbackp (2,basic_g(icm)%vc(1,1,1),pbuff(1+iptr)  &
          ,basic_g(icm)%dn0(1,1,1),basic_g(icm)%dn0u(1,1,1) &
          ,basic_g(icm)%dn0v(1,1,1)  &
          ,mmzp(icm),mmxp(icm),mmyp(icm),mi0(icm),mj0(icm)  &
          ,i1s,i2s,j1s,j2s  &
          ,k1s,k2s,i2u,j2v,k2w,nfz,nfx,nfy)
      iptr=iptr+mtp
      CALL unfdbackp (3,basic_g(icm)%wc(1,1,1),pbuff(1+iptr)  &
          ,basic_g(icm)%dn0(1,1,1),basic_g(icm)%dn0u(1,1,1) &
          ,basic_g(icm)%dn0v(1,1,1)  &
          ,mmzp(icm),mmxp(icm),mmyp(icm),mi0(icm),mj0(icm)  &
          ,i1s,i2s,j1s,j2s  &
          ,k1s,k2s,i2u,j2v,k2w,nfz,nfx,nfy)
      iptr=iptr+mtp
      CALL unfdbackp (4,basic_g(icm)%pc(1,1,1),pbuff(1+iptr)  &
          ,basic_g(icm)%dn0(1,1,1),basic_g(icm)%dn0u(1,1,1) &
          ,basic_g(icm)%dn0v(1,1,1)  &
          ,mmzp(icm),mmxp(icm),mmyp(icm),mi0(icm),mj0(icm)  &
          ,i1s,i2s,j1s,j2s  &
          ,k1s,k2s,i2u,j2v,k2w,nfz,nfx,nfy)
      iptr=iptr+mtp

      do nv=1,num_scalar(ifm)
         CALL unfdbackp (5,scalar_tab(nv,icm)%var_p,pbuff(1+iptr)  &
             ,basic_g(icm)%dn0(1,1,1),basic_g(icm)%dn0u(1,1,1) &
             ,basic_g(icm)%dn0v(1,1,1)  &
             ,mmzp(icm),mmxp(icm),mmyp(icm),mi0(icm),mj0(icm)  &
             ,i1s,i2s,j1s,j2s  &
             ,k1s,k2s,i2u,j2v,k2w,nfz,nfx,nfy)
         iptr=iptr+mtp
      enddo

   endif

enddo

if (nnstbot(icm) == 1) then
   CALL botset (mmzp(icm),mmxp(icm),mmyp(icm),basic_g(icm)%uc(1,1,1),'U')
   CALL botset (mmzp(icm),mmxp(icm),mmyp(icm),basic_g(icm)%vc(1,1,1),'V')
   CALL botset (mmzp(icm),mmxp(icm),mmyp(icm),basic_g(icm)%pc(1,1,1),'P')
endif

if (nnsttop(icm) == 1) then
   CALL topset (mmzp(icm),mmxp(icm),mmyp(icm),basic_g(icm)%uc(1,1,1),'U')
   CALL topset (mmzp(icm),mmxp(icm),mmyp(icm),basic_g(icm)%vc(1,1,1),'V')
   CALL topset (mmzp(icm),mmxp(icm),mmyp(icm),basic_g(icm)%pc(1,1,1),'P')
endif

do nv=1,num_scalar(ifm)
   if (nnstbot(icm) == 1) then   
      CALL botset (mmzp(icm),mmxp(icm),mmyp(icm),scalar_tab(nv,icm)%var_p,'T')
   endif

   if (nnsttop(icm) == 1) then
      CALL topset (mmzp(icm),mmxp(icm),mmyp(icm),scalar_tab(nv,icm)%var_p,'T')
   endif
enddo

return
END SUBROUTINE node_getfeed

!##############################################################################
Subroutine zeroout (ac,m1,m2,m3,i0,j0,i1z,i2z,j1z,j2z,k1z,k2z)

implicit none

integer :: m1,m2,m3,i0,j0,i1z,i2z,j1z,j2z,k1z,k2z
real :: ac(m1,m2,m3)
integer :: i,j,k

do j=j1z,j2z
   do i=i1z,i2z
      do k=k1z,k2z
         ac(k,i-i0,j-j0)=0.
      enddo
   enddo
enddo

return
END SUBROUTINE zeroout

!##############################################################################
Subroutine unfdbackp (ivarn,ac,acf,den,denu,denv,m1,m2,m3,i0,j0  &
   ,i1s,i2s,j1s,j2s,k1s,k2s,i2u,j2v,k2w,nf1,nf2,nf3)

implicit none

integer :: ivarn,m1,m2,m3,i0,j0  &
          ,i1s,i2s,j1s,j2s,k1s,k2s,i2u,j2v,k2w,nf1,nf2,nf3
real :: ac(m1,m2,m3),acf(nf1,nf2,nf3),den(m1,m2,m3)  &
       ,denu(m1,m2,m3),denv(m1,m2,m3)
integer :: i,j,k

!     ivarn = variable types 1- u
!                            2- v
!                            3- w
!                            4- p
!                            5- scalar

if(ivarn.ge.5) then
   do j=j1s,j2s
      do i=i1s,i2s
         do k=k1s,k2s
            ac(k,i-i0,j-j0)=ac(k,i-i0,j-j0)  &
                 +acf(k-k1s+1,i-i1s+1,j-j1s+1)/den(k,i-i0,j-j0)
         enddo
      enddo
   enddo
elseif(ivarn==1) then
   do j=j1s,j2s
      do i=i1s,i2u
         do k=k1s,k2s
            ac(k,i-i0,j-j0)=ac(k,i-i0,j-j0)  &
                 +acf(k-k1s+1,i-i1s+1,j-j1s+1)  &
                 /denu(k,i-i0,j-j0)
         enddo
      enddo
   enddo
elseif(ivarn==2) then
   do j=j1s,j2v
      do i=i1s,i2s
         do k=k1s,k2s
            ac(k,i-i0,j-j0)=ac(k,i-i0,j-j0)  &
                 +acf(k-k1s+1,i-i1s+1,j-j1s+1)  &
                 /denv(k,i-i0,j-j0)
         enddo
      enddo
   enddo
elseif(ivarn==3) then
   do j=j1s,j2s
      do i=i1s,i2s
         do k=k1s,k2w
            ac(k,i-i0,j-j0)=ac(k,i-i0,j-j0)  &
                 +acf(k-k1s+1,i-i1s+1,j-j1s+1)  &
                 /(den(k,i-i0,j-j0)+den(k+1,i-i0,j-j0))
         enddo
      enddo
   enddo
else
   do j=j1s,j2s
      do i=i1s,i2s
         do k=k1s,k2s
            ac(k,i-i0,j-j0)=ac(k,i-i0,j-j0)  &
                 +acf(k-k1s+1,i-i1s+1,j-j1s+1)
         enddo
      enddo
   enddo
endif

return
END SUBROUTINE unfdbackp
