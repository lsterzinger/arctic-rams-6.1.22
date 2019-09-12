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

Subroutine master_sendinit ()

use var_tables
use mem_scratch
use mem_cuparm
use mem_varinit
use mem_grid
use io_params
use rpara

implicit none

integer :: nm,nwords,mxp,myp,mxyp,mxyzp,nv,npts,ng

do nm = 1,nmachs

   nwords = 4 * maxgrds + 6
   CALL par_init_put (vctr1,nwords)

   CALL par_put_float (vtime1,1)
   CALL par_put_float (vtime2,1)
   CALL par_put_float (htime1,1)
   CALL par_put_float (htime2,1)
   CALL par_put_float (ssttime1,maxgrds)
   CALL par_put_float (ssttime2,maxgrds)
   CALL par_put_float (ndvitime1,maxgrds)
   CALL par_put_float (ndvitime2,maxgrds)

   CALL par_send (machnum(nm),110)

   do ng = 1,ngrids

      mxp = nxend(nm,ng) - nxbeg(nm,ng) + 1
      myp = nyend(nm,ng) - nybeg(nm,ng) + 1
      mxyp = mxp * myp
      mxyzp = mxyp * nnzp(ng)

!      Send variables that the nodes will only have the subdomain portion of.

      do nv = 1,num_var(ng)
      
         if ( vtab_r(nv,ng)%impti == 1) then
      
         if ( vtab_r(nv,ng)%idim_type == 2) then

            CALL mk_2_buff (vtab_r(nv,ng)%var_p  &
                          ,scratch%scr2(1),nnxp(ng),nnyp(ng),mxp,myp  &
                          ,nxbeg(nm,ng),nxend(nm,ng)  &
                          ,nybeg(nm,ng),nyend(nm,ng))
            npts=mxyp

         elseif ( vtab_r(nv,ng)%idim_type == 3) then

            CALL mk_3_buff (vtab_r(nv,ng)%var_p  &
                          ,scratch%scr2(1),nnzp(ng),nnxp(ng),nnyp(ng)  &
                          ,nnzp(ng),mxp,myp  &
                          ,nxbeg(nm,ng),nxend(nm,ng)  &
                          ,nybeg(nm,ng),nyend(nm,ng))       
            npts=mxyzp

         elseif ( vtab_r(nv,ng)%idim_type == 4) then

            CALL mk_4_buff (vtab_r(nv,ng)%var_p  &
                          ,scratch%scr2(1),nzg,nnxp(ng),nnyp(ng)  &
                          ,npatch,nzg,mxp,myp,npatch &
                          ,nxbeg(nm,ng),nxend(nm,ng)  &
                          ,nybeg(nm,ng),nyend(nm,ng))
            npts=mxyp*nzg*npatch
         
         elseif ( vtab_r(nv,ng)%idim_type == 5) then

            CALL mk_4_buff (vtab_r(nv,ng)%var_p  &
                          ,scratch%scr2(1),nzs,nnxp(ng),nnyp(ng)  &
                          ,npatch,nzs,mxp,myp,npatch &
                          ,nxbeg(nm,ng),nxend(nm,ng)  &
                          ,nybeg(nm,ng),nyend(nm,ng))
            npts=mxyp*nzs*npatch
         
         elseif ( vtab_r(nv,ng)%idim_type == 6) then

            CALL mk_2p_buff (vtab_r(nv,ng)%var_p  &
                          ,scratch%scr2(1),nnxp(ng),nnyp(ng)  &
                          ,npatch,mxp,myp,npatch &
                          ,nxbeg(nm,ng),nxend(nm,ng)  &
                          ,nybeg(nm,ng),nyend(nm,ng))
            npts=mxyp*npatch
         
         endif
         
         CALL par_init_put (scratch%scr1(1),npts)
         CALL par_put_float (scratch%scr2(1),npts)
         CALL par_send (machnum(nm),110+ng)
         
         endif

      enddo

   enddo
enddo

return
END SUBROUTINE master_sendinit

!##############################################################################
Subroutine mk_2_buff (a,b,n1,n2,m1,m2,i1,i2,j1,j2)

implicit none

integer :: n1,n2,m1,m2,i1,i2,j1,j2
real :: a(n1,n2),b(m1,m2)
   
   b(1:m1,1:m2)=a(i1:i2,j1:j2)

return
END SUBROUTINE mk_2_buff

!##############################################################################
Subroutine mk_2p_buff (a,b,n1,n2,n3,m1,m2,m3,i1,i2,j1,j2)

implicit none

integer :: n1,n2,n3,m1,m2,m3,i1,i2,j1,j2
real :: a(n1,n2,n3),b(m1,m2,m3)
   
   b(1:m1,1:m2,1:m3)=a(i1:i2,j1:j2,1:n3)

return
END SUBROUTINE mk_2p_buff

!##############################################################################
Subroutine mk_3_buff (a,b,n1,n2,n3,m1,m2,m3,i1,i2,j1,j2)

implicit none

integer :: n1,n2,n3,m1,m2,m3,i1,i2,j1,j2
real :: a(n1,n2,n3),b(m1,m2,m3)
   
   b(1:m1,1:m2,1:m3)=a(1:n1,i1:i2,j1:j2)

return
END SUBROUTINE mk_3_buff

!##############################################################################
Subroutine mk_4_buff (a,b,n1,n2,n3,n4,m1,m2,m3,m4,i1,i2,j1,j2)

implicit none

integer :: n1,n2,n3,n4,m1,m2,m3,m4,i1,i2,j1,j2
real :: a(n1,n2,n3,m4),b(m1,m2,m3,m4)
   
   b(1:m1,1:m2,1:m3,1:m4)=a(1:n1,i1:i2,j1:j2,1:n4)

return
END SUBROUTINE mk_4_buff

!##############################################################################
Subroutine node_getinit ()

use var_tables
use mem_scratch
use node_mod
use mem_cuparm
use mem_varinit
use mem_grid
use io_params

implicit none

integer :: nwords,ibytes,msgtype,ihostnum,nv,ng

nwords = 4 * maxgrds + 6
CALL par_get_new (vctr1,nwords,110,ibytes,msgtype,ihostnum)

CALL par_get_float (vtime1,1)
CALL par_get_float (vtime2,1)
CALL par_get_float (htime1,1)
CALL par_get_float (htime2,1)
CALL par_get_float (ssttime1,maxgrds)
CALL par_get_float (ssttime2,maxgrds)
CALL par_get_float (ndvitime1,maxgrds)
CALL par_get_float (ndvitime2,maxgrds)

do ng = 1,ngrids

   do nv=1,num_var(ng)
      
      if (vtab_r(nv,ng)%impti == 1) then
         CALL par_get_new (scratch%scr2(1),vtab_r(nv,ng)%npts  &
              ,110+ng,ibytes,msgtype,ihostnum)
         CALL par_get_float (vtab_r(nv,ng)%var_p,vtab_r(nv,ng)%npts)
      endif
   
   enddo
   
enddo

return
END SUBROUTINE node_getinit

!##############################################################################
Subroutine node_sendall ()

use node_mod
use mem_grid
use var_tables
use mem_scratch

implicit none

integer, parameter :: ntags=5
integer :: msgtags(ntags),msgnum,nv,ng

msgnum=1001

do ng=1,ngrids

   msgtags(1)=mynum
   msgtags(2)=ng

   do nv=1,num_var(ng)

      if (vtab_r(nv,ng)%impt3 == 1) then
         msgtags(3)= vtab_r(nv,ng)%npts
         msgtags(4)= vtab_r(nv,ng)%idim_type
         msgtags(5)= nv
         CALL par_init_put (scratch%scr1(1),vtab_r(nv,ng)%npts+ntags)
         CALL par_put_int (msgtags,ntags)
         CALL par_put_float (vtab_r(nv,ng)%var_p,vtab_r(nv,ng)%npts)
         CALL par_send (master_num,msgnum)
      endif
      
   enddo

enddo

return
END SUBROUTINE node_sendall

!##############################################################################
Subroutine master_getall ()

use mem_grid
use rpara
use var_tables
use mem_scratch

implicit none

integer, parameter :: ntags=5
integer :: msgtags(ntags),numvars,nummess,nvvv,mxp,myp,mxyp
integer :: ibytes,msgtyp,ihostnum
integer :: nm,il1,ir2,jb1,jt2,nv,idim_type,npts,ng


numvars=0
do ng=1,ngrids
   do nv=1,num_var(ng)
      if (vtab_r(nv,ng)%impt3 == 1) numvars=numvars+1
   enddo
enddo
nummess=numvars*nmachs

do nvvv=1,nummess

   CALL par_get_new (scratch%scr1(1),ubound(scratch%scr1,1)  &
                   ,1001,ibytes,msgtyp,ihostnum)

   CALL par_get_int (msgtags,ntags)
   nm=msgtags(1)
   ng=msgtags(2)
   npts=msgtags(3)
   idim_type=msgtags(4)
   nv=msgtags(5)

   il1=nxbegc(nm,ng)
   ir2=nxendc(nm,ng)
   jb1=nybegc(nm,ng)
   jt2=nyendc(nm,ng)

   if(iand(ibcflg(nm,ng),1).ne.0) il1=il1-1
   if(iand(ibcflg(nm,ng),2).ne.0) ir2=ir2+1
   if(iand(ibcflg(nm,ng),4).ne.0) jb1=jb1-1
   if(iand(ibcflg(nm,ng),8).ne.0) jt2=jt2+1

   mxp = nxend(nm,ng) - nxbeg(nm,ng) + 1
   myp = nyend(nm,ng) - nybeg(nm,ng) + 1
   mxyp = mxp * myp

   CALL par_get_float (scratch%scr1(1),npts)
   
   if (idim_type == 2) then
      CALL ex_2_buff (vtab_r(nv,ng)%var_p,scratch%scr1(1)  &
          ,nnxp(ng),nnyp(ng) ,mxp,myp  &
          ,ixoff(nm,ng),iyoff(nm,ng),il1,ir2,jb1,jt2)
   elseif (idim_type == 3) then
      CALL ex_3_buff (vtab_r(nv,ng)%var_p,scratch%scr1(1)  &
          ,nnzp(ng),nnxp(ng),nnyp(ng),nnzp(ng),mxp,myp  &
          ,ixoff(nm,ng),iyoff(nm,ng),il1,ir2,jb1,jt2)
   elseif (idim_type == 4) then
      CALL ex_4_buff (vtab_r(nv,ng)%var_p,scratch%scr1(1)  &
          ,nzg,nnxp(ng),nnyp(ng),npatch,nzg,mxp,myp,npatch  &
          ,ixoff(nm,ng),iyoff(nm,ng),il1,ir2,jb1,jt2)
   elseif (idim_type == 5) then
      CALL ex_4_buff (vtab_r(nv,ng)%var_p,scratch%scr1(1)  &
          ,nzs,nnxp(ng),nnyp(ng),npatch,nzs,mxp,myp,npatch  &
          ,ixoff(nm,ng),iyoff(nm,ng),il1,ir2,jb1,jt2)
   elseif (idim_type == 6) then
      CALL ex_2p_buff (vtab_r(nv,ng)%var_p,scratch%scr1(1)  &
          ,nnxp(ng),nnyp(ng),npatch,mxp,myp,npatch  &
          ,ixoff(nm,ng),iyoff(nm,ng),il1,ir2,jb1,jt2)
   endif

enddo

return
END SUBROUTINE master_getall

!##############################################################################
Subroutine ex_2_buff (a,b,n1,n2,m1,m2,i0,j0,i1,i2,j1,j2)

implicit none

integer :: n1,n2,m1,m2,i0,j0,i1,i2,j1,j2
real :: a(n1,n2),b(m1,m2)

   a(i1+i0:i2+i0,j1+j0:j2+j0) = b(i1:i2,j1:j2)

return
END SUBROUTINE ex_2_buff

!##############################################################################
Subroutine ex_3_buff (a,b,n1,n2,n3,m1,m2,m3,i0,j0,i1,i2,j1,j2)

implicit none

integer :: n1,n2,n3,m1,m2,m3,i0,j0,i1,i2,j1,j2
real :: a(n1,n2,n3),b(m1,m2,m3)

   a(1:n1,i1+i0:i2+i0,j1+j0:j2+j0) = b(1:m1,i1:i2,j1:j2)

return
END SUBROUTINE ex_3_buff

!##############################################################################
Subroutine ex_4_buff (a,b,n1,n2,n3,n4,m1,m2,m3,m4,i0,j0,i1,i2,j1,j2)

implicit none

integer :: n1,n2,n3,n4,m1,m2,m3,m4,i0,j0,i1,i2,j1,j2
real :: a(n1,n2,n3,n4),b(m1,m2,m3,m4)

   a(1:n1,i1+i0:i2+i0,j1+j0:j2+j0,1:n4) = b(1:m1,i1:i2,j1:j2,1:m4)

return
END SUBROUTINE ex_4_buff

!##############################################################################
Subroutine ex_2p_buff (a,b,n1,n2,n3,m1,m2,m3,i0,j0,i1,i2,j1,j2)

implicit none

integer :: n1,n2,n3,m1,m2,m3,i0,j0,i1,i2,j1,j2
real :: a(n1,n2,n3),b(m1,m2,m3)

   a(i1+i0:i2+i0,j1+j0:j2+j0,1:n3) = b(i1:i2,j1:j2,1:m3)

return
END SUBROUTINE ex_2p_buff

!##############################################################################
Subroutine node_sendanl (vtype)

!Send just the 'vtype' variables back to the master process

use node_mod
use mem_grid
use var_tables
use mem_scratch

implicit none

character(len=*) :: vtype
integer, parameter :: ntags=5
integer :: msgtags(ntags),msgnum,nv,isend,ng
real, pointer :: v_pointer => NULL()

msgnum=1002
if(iprntstmt>=2) print*,'call node-send:',ngrids,trim(vtype)

do ng=1,ngrids

   msgtags(1)=mynum
   msgtags(2)=ng

   do nv=1,num_var(ng)

      if(iprntstmt>=2)print*,'node-send:',nv,num_var(ng) &
                       ,vtab_r(nv,ngrid)%npts,trim(vtab_r(nv,ng)%name)
      isend=0
      if(vtype == 'LITE' .and. vtab_r(nv,ng)%ilite == 1) then
         isend=1
         v_pointer => vtab_r(nv,ng)%var_p
      elseif( (vtype == 'MEAN' .or. vtype == 'BOTH')  &
                .and. vtab_r(nv,ng)%imean == 1) then
         isend=1
         v_pointer => vtab_r(nv,ng)%var_m
      endif

      if(isend == 1) then
         msgtags(3)= vtab_r(nv,ng)%npts
         msgtags(4)= vtab_r(nv,ng)%idim_type
         msgtags(5)= nv
         CALL par_init_put (scratch%scr1(1),vtab_r(nv,ng)%npts+ntags)
         CALL par_put_int (msgtags,ntags)
         CALL par_put_float (v_pointer, vtab_r(nv,ng)%npts)
         CALL par_send (master_num,msgnum)
      endif

   enddo

enddo

return
END SUBROUTINE node_sendanl

!##############################################################################
Subroutine master_getanl (vtype)

!Get just the 'vtype' variables from the nodes

use mem_grid
use rpara
use var_tables
use mem_scratch

implicit none

character(len=*) :: vtype
integer, parameter :: ntags=5
integer :: msgtags(ntags),numvars,nummess,nvvv,mxp,myp,mxyp &
          ,ibytes,msgtyp,ihostnum,nm,il1,ir2,jb1,jt2,nv,idim_type,npts,ng
real, pointer :: v_pointer => NULL()

numvars=0
do ng=1,ngrids
   do nv=1,num_var(ng)
      if( (vtype == 'LITE' .and. vtab_r(nv,ng)%ilite == 1) .or. &
          ( (vtype == 'MEAN' .or. vtype == 'BOTH') &
                           .and. vtab_r(nv,ng)%imean == 1) ) then
         numvars=numvars+1
      endif
   enddo
enddo
nummess=numvars*nmachs

do nvvv=1,nummess

   CALL par_get_new (scratch%scr1(1),ubound(scratch%scr1,1)  &
                   ,1002,ibytes,msgtyp,ihostnum)

   CALL par_get_int (msgtags,ntags)
   nm=msgtags(1)
   ng=msgtags(2)
   npts=msgtags(3)
   idim_type=msgtags(4)
   nv=msgtags(5)
   
   ! Find beginning of the variable array for this message
   if(vtype == 'LITE') then
      v_pointer => vtab_r(nv,ng)%var_p
   elseif( vtype == 'MEAN' .or. vtype == 'BOTH' ) then
      v_pointer => vtab_r(nv,ng)%var_m
   endif
   
   ! Boundary limits and number of points for this section of the variable
   il1=nxbegc(nm,ng)
   ir2=nxendc(nm,ng)
   jb1=nybegc(nm,ng)
   jt2=nyendc(nm,ng)

   if(iand(ibcflg(nm,ng),1) /= 0) il1=il1-1
   if(iand(ibcflg(nm,ng),2) /= 0) ir2=ir2+1
   if(iand(ibcflg(nm,ng),4) /= 0) jb1=jb1-1
   if(iand(ibcflg(nm,ng),8) /= 0) jt2=jt2+1

   mxp = nxend(nm,ng) - nxbeg(nm,ng) + 1
   myp = nyend(nm,ng) - nybeg(nm,ng) + 1
   mxyp = mxp * myp

   CALL par_get_float (scratch%scr1(1),npts)
   
   if (idim_type == 2) then
      CALL ex_2_buff (v_pointer,scratch%scr1(1)  &
          ,nnxp(ng),nnyp(ng) ,mxp,myp  &
          ,ixoff(nm,ng),iyoff(nm,ng),il1,ir2,jb1,jt2)
   elseif (idim_type == 3) then
      CALL ex_3_buff (v_pointer,scratch%scr1(1)  &
          ,nnzp(ng),nnxp(ng),nnyp(ng),nnzp(ng),mxp,myp  &
          ,ixoff(nm,ng),iyoff(nm,ng),il1,ir2,jb1,jt2)
   elseif (idim_type == 4) then
      CALL ex_4_buff (v_pointer,scratch%scr1(1)  &
          ,nzg,nnxp(ng),nnyp(ng),npatch,nzg,mxp,myp,npatch  &
          ,ixoff(nm,ng),iyoff(nm,ng),il1,ir2,jb1,jt2)
   elseif (idim_type == 5) then
      CALL ex_4_buff (v_pointer,scratch%scr1(1)  &
          ,nzs,nnxp(ng),nnyp(ng),npatch,nzs,mxp,myp,npatch  &
          ,ixoff(nm,ng),iyoff(nm,ng),il1,ir2,jb1,jt2)
   elseif (idim_type == 6) then
      CALL ex_2p_buff (v_pointer,scratch%scr1(1)  &
          ,nnxp(ng),nnyp(ng),npatch,mxp,myp,npatch  &
          ,ixoff(nm,ng),iyoff(nm,ng),il1,ir2,jb1,jt2)
   endif

enddo

return
END SUBROUTINE master_getanl
