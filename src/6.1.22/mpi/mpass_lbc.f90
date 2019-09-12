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

!############################################################################
Subroutine node_sendlbc (isflag)

! This routine along with node_getlbc() handles the updating of internal
! (internodal) overlap regions. These routines need to be run in sync with
! the node_sendcyclic() and node_getcyclic() routines. Run these routines
! right before calling the cyclic routines.
!
! Keep the coding for isflag in sync with that in the routine update_cyclic().
!
!   isflag                  var
!     0          all vars (from vtab_r)
!     1          all scalar vars (from scalar_tab)
!     2          past u velocity (up)
!     3          past v velocity (vp)
!     4          past p velocity (pp)
!     5          past w velocity (wp)
!
!

use mem_grid
use node_mod
use var_tables
use mem_scratch
use mem_basic

implicit none

integer :: isflag
integer :: itype,nm,i1,i2,j1,j2,nv,mtp
real, pointer :: scalarp

itype=1

!______________________
!
!   First, before we send anything, let's post the receives. Also, make sure
!     any pending sends are complete.

do nm=1,nmachs
   if (iget_paths(itype,ngrid,nm).ne.0) then
      CALL par_get_noblock (node_buffs(nm)%lbc_recv_buff(1)  &
          ,node_buffs(nm)%nrecv ,20000+ngrid,machs(nm),irecv_req(nm) )

   endif
enddo


!______________________
!
!   Now we can actually go on to sending the stuff

do nm=1,nmachs

   if(ipaths(1,itype,ngrid,nm).ne.0) then

      i1=ipaths(1,itype,ngrid,nm)
      i2=ipaths(2,itype,ngrid,nm)
      j1=ipaths(3,itype,ngrid,nm)
      j2=ipaths(4,itype,ngrid,nm)

      CALL par_init_put (node_buffs(nm)%lbc_send_buff(1)  &
                       ,node_buffs(nm)%nsend )
      CALL par_put_int (i1,1)
      CALL par_put_int (i2,1)
      CALL par_put_int (j1,1)
      CALL par_put_int (j2,1)
      CALL par_put_int (mynum,1)

      select case (isflag)
        case (0)
          do nv = 1,num_var(ngrid)
             if ( vtab_r(nv,ngrid)%impt1 == 1) then
                if ( vtab_r(nv,ngrid)%idim_type == 3) then
                   CALL mklbcbuff (mzp,mxp,myp,vtab_r(nv,ngrid)%var_p  &
                       ,scratch%vt3dp(1),i1-i0,i2-i0,j1-j0,j2-j0,mtp)
                   CALL par_put_float (scratch%vt3dp(1),mtp)
                elseif ( vtab_r(nv,ngrid)%idim_type == 2) then
                   CALL mklbcbuff (1,mxp,myp,vtab_r(nv,ngrid)%var_p  &
                       ,scratch%vt3dp(1),i1-i0,i2-i0,j1-j0,j2-j0,mtp)
                   CALL par_put_float (scratch%vt3dp(1),mtp)
                endif
             endif
          enddo
        case (1)
            ! scalar vars (all scalars are 3D)
            do nv = 1, num_scalar(ngrid)
              scalarp => scalar_tab(nv,ngrid)%var_p
              CALL mklbcbuff (mzp,mxp,myp,scalarp  &
                  ,scratch%scr1(1),i1-i0,i2-i0,j1-j0,j2-j0,mtp)
              CALL par_put_float (scratch%scr1(1),mtp)
            enddo
        case (2)
          CALL mklbcbuff (mzp,mxp,myp,basic_g(ngrid)%up(1,1,1)  &
              ,scratch%scr1(1),i1-i0,i2-i0,j1-j0,j2-j0,mtp)
          CALL par_put_float (scratch%scr1(1),mtp)
        case (3)
          CALL mklbcbuff (mzp,mxp,myp,basic_g(ngrid)%vp(1,1,1)  &
              ,scratch%scr1(1),i1-i0,i2-i0,j1-j0,j2-j0,mtp)
          CALL par_put_float (scratch%scr1(1),mtp)
        case (4)
          CALL mklbcbuff (mzp,mxp,myp,basic_g(ngrid)%pp(1,1,1)  &
              ,scratch%scr1(1),i1-i0,i2-i0,j1-j0,j2-j0,mtp)
          CALL par_put_float (scratch%scr1(1),mtp)
        case (5)
          CALL mklbcbuff (mzp,mxp,myp,basic_g(ngrid)%wp(1,1,1)  &
              ,scratch%scr1(1),i1-i0,i2-i0,j1-j0,j2-j0,mtp)
          CALL par_put_float (scratch%scr1(1),mtp)
      endselect

      CALL par_send_noblock (ipaths(5,itype,ngrid,nm)  &
           ,20000+ngrid,isend_req(nm))

   endif

enddo

return
END SUBROUTINE node_sendlbc

!############################################################################
Subroutine node_getlbc (isflag)

! This routine along with node_getlbc() handles the updating of internal
! (internodal) overlap regions. These routines need to be run in sync with
! the node_sendcyclic() and node_getcyclic() routines. Run these routines
! right before calling the cyclic routines.
!
! Keep the coding for isflag in sync with that in the routine update_cyclic().
!
!   isflag                  var
!     0          all vars (from vtab_r)
!     1          all scalar vars (from scalar_tab)
!     2          past u velocity (up)
!     3          past v velocity (vp)
!     4          past p velocity (pp)
!     5          past w velocity (wp)
!

use mem_grid
use node_mod
use var_tables
use mem_scratch
use mem_basic

implicit none

integer :: isflag
integer :: itype,nm,ibytes,msgid,ihostnum,i1,i2,j1,j2  &
          ,nv,node_src,mtc,mtp3,mtp2
real, pointer :: scalarp

itype=1

!_____________________________________________________________________
!
!  First, let's make sure our sends are all finished and de-allocated

do nm=1,nmachs
   if(ipaths(1,itype,ngrid,nm).ne.0) then
      CALL par_wait (isend_req(nm),ibytes,msgid,ihostnum)
   endif
enddo
!_____________________________________________________________________
!
!  Now, let's wait on our receives

do nm=1,nmachs
   if (iget_paths(itype,ngrid,nm).ne.0) then
      CALL par_wait (irecv_req(nm),ibytes,msgid,ihostnum)
   endif
enddo
!_____________________________________________________________________
!
!  We got all our stuff. Now unpack it into appropriate space.

do nm=1,nmachs

   if (iget_paths(itype,ngrid,nm).ne.0) then


      CALL par_assoc_buff (node_buffs(nm)%lbc_recv_buff(1)  &
                         ,node_buffs(nm)%nrecv) 


      CALL par_get_int (i1,1)
      CALL par_get_int (i2,1)
      CALL par_get_int (j1,1)
      CALL par_get_int (j2,1)
      CALL par_get_int (node_src,1)

      mtp2 = (i2-i1+1)*(j2-j1+1)
      mtp3 = nnzp(ngrid) * mtp2
      
      select case (isflag)
        case (0)
          do nv = 1,num_var(ngrid)
             if ( vtab_r(nv,ngrid)%impt1 == 1) then
                if ( vtab_r(nv,ngrid)%idim_type == 3) then
                   CALL par_get_float (scratch%vt3dp(1),mtp3)
                   CALL exlbcbuff (mzp,mxp,myp,vtab_r(nv,ngrid)%var_p  &
                       ,scratch%vt3dp(1),i1-i0,i2-i0,j1-j0,j2-j0,mtc)
                elseif ( vtab_r(nv,ngrid)%idim_type == 2) then
                   CALL par_get_float (scratch%vt3dp(1),mtp2)
                   CALL exlbcbuff (1,mxp,myp,vtab_r(nv,ngrid)%var_p  &
                       ,scratch%vt3dp(1),i1-i0,i2-i0,j1-j0,j2-j0,mtc)
                endif
             endif
          enddo
        case (1)
            ! scalar vars (all scalars are 3D)
            do nv = 1, num_scalar(ngrid)
              scalarp => scalar_tab(nv,ngrid)%var_p
              CALL par_get_float (scratch%scr1(1),mtp3)
              CALL exlbcbuff (mzp,mxp,myp,scalarp  &
                  ,scratch%scr1(1),i1-i0,i2-i0,j1-j0,j2-j0,mtc)
            enddo
        case (2)
          CALL par_get_float (scratch%scr1(1),mtp3)
          CALL exlbcbuff (mzp,mxp,myp,basic_g(ngrid)%up(1,1,1)  &
              ,scratch%scr1(1),i1-i0,i2-i0,j1-j0,j2-j0,mtc)
        case (3)
          CALL par_get_float (scratch%scr1(1),mtp3)
          CALL exlbcbuff (mzp,mxp,myp,basic_g(ngrid)%vp(1,1,1)  &
              ,scratch%scr1(1),i1-i0,i2-i0,j1-j0,j2-j0,mtc)
        case (4)
          CALL par_get_float (scratch%scr1(1),mtp3)
          CALL exlbcbuff (mzp,mxp,myp,basic_g(ngrid)%pp(1,1,1)  &
              ,scratch%scr1(1),i1-i0,i2-i0,j1-j0,j2-j0,mtc)
        case (5)
          CALL par_get_float (scratch%scr1(1),mtp3)
          CALL exlbcbuff (mzp,mxp,myp,basic_g(ngrid)%wp(1,1,1)  &
              ,scratch%scr1(1),i1-i0,i2-i0,j1-j0,j2-j0,mtc)
      endselect
   endif
enddo

return
END SUBROUTINE node_getlbc

!##############################################################################
Subroutine mklbcbuff (n1,n2,n3,a,b,il,ir,jb,jt,ind)

implicit none

integer :: n1,n2,n3,il,ir,jb,jt,ind
real :: a(n1,n2,n3),b(*)
integer :: i,j,k

ind=0
do j=jb,jt
   do i=il,ir
      do k=1,n1
         ind=ind+1
         b(ind)=a(k,i,j)
      enddo
   enddo
enddo

return
END SUBROUTINE mklbcbuff

!##############################################################################
Subroutine exlbcbuff (n1,n2,n3,a,b,il,ir,jb,jt,ind)

implicit none

integer :: n1,n2,n3,il,ir,jb,jt,ind
real :: a(n1,n2,n3),b(*)
integer :: i,j,k

ind=0
do j=jb,jt
   do i=il,ir
      do k=1,n1
         ind=ind+1
         a(k,i,j)=b(ind)
      enddo
   enddo
enddo

return
END SUBROUTINE exlbcbuff

