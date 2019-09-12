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

Module mem_scratch

use grid_dims

implicit none

   Type scratch_vars

   real, allocatable, dimension(:) :: &
                         scr1,scr2 &
                        ,vt3da,vt3db,vt3dc,vt3dd,vt3de,vt3df,vt3dg  &
                        ,vt3dh,vt3di,vt3dj,vt3dk,vt3dl,vt3dm,vt3dn  &
                        ,vt3do,vt3dp
   real, allocatable, dimension(:) :: &
                        vt2da,vt2db,vt2dc,vt2dd,vt2de,vt2df
                         
   End Type
   
   type (scratch_vars) :: scratch

   real, dimension(maxdim) ::  vctr1 ,vctr2 ,vctr3 ,vctr4 ,vctr5 ,vctr6   &
                              ,vctr7 ,vctr8 ,vctr9 ,vctr10,vctr11,vctr12  &
                              ,vctr13,vctr14,vctr15,vctr16,vctr17,vctr18  &
                              ,vctr19,vctr20,vctr21,vctr22,vctr23,vctr24  &
                              ,vctr25,vctr26,vctr27,vctr28,vctr29,vctr30  &
                              ,vctr31,vctr32,vctr33,vctr34,vctr35,vctr36  &
                              ,vctr37,vctr38,vctr39,vctr40,vctr41

Contains

!##############################################################################
Subroutine alloc_scratch (nmzp,nmxp,nmyp,nnzp,nnxp,nnyp  &
              ,ngrs,nzg,nzs,npatch,proc_type,maxx,maxy,maxz)

implicit none
   
   integer, dimension (*) :: nmzp,nmxp,nmyp,nnzp,nnxp,nnyp
   integer :: ngrs,nzg,nzs,npatch,proc_type
   integer :: ng,ntpts,ntpts1,ntpts2,ntptsx,maxx,maxy,maxz

!print*, 'enter alloc_scratch'

!         Find the maximum number of grid points needed for any grid.
!           The max points in each direction are passed back for use by 
!           various nesting things.    

   maxx = 0
   maxy = 0
   maxz = 0
   ntpts=0
   ntpts2=0
   do ng=1,ngrs
      maxx = max(maxx,nnxp(ng))
      maxy = max(maxy,nnyp(ng))
      maxz = max(maxz,nnzp(ng))
      ntpts=max( nmxp(ng)*nmyp(ng)*nmzp(ng),ntpts )
      ntpts2=max( nmxp(ng)*nmyp(ng),ntpts2 )
   enddo
   ! scr1 and scr2 needs to be the max of a passed field
   ntptsx=max(maxx*maxy*maxz,ntpts2*nzg*npatch,ntpts2*nzs*npatch)+1000
   if(proc_type==1) then
      ntpts1=1
   else
      ntpts1=ntpts
   endif

! Allocate arrays based on options (if necessary).
!   -scr1 and scr2 need to be allocated to full domain (even on compute nodes)
!        to max(nx)*max(ny)*max(nz)
!   -do not need all these arrays if it is a master process in a parallel run,
!      so just allocate some to 1 word.

      allocate (scratch%scr1 (ntptsx))
      allocate (scratch%scr2 (ntptsx))
      allocate (scratch%vt3da(ntpts))
      allocate (scratch%vt3db(ntpts))     

      allocate (scratch%vt3dc(ntpts1))
      allocate (scratch%vt3dd(ntpts1))
      allocate (scratch%vt3de(ntpts1))
      allocate (scratch%vt3df(ntpts1))
      allocate (scratch%vt3dg(ntpts1))
      allocate (scratch%vt3dh(ntpts1))
      allocate (scratch%vt3di(ntpts1))
      allocate (scratch%vt3dj(ntpts1))
      allocate (scratch%vt3dk(ntpts1))
      allocate (scratch%vt3dl(ntpts1))
      allocate (scratch%vt3dm(ntpts1))
      allocate (scratch%vt3dn(ntpts1))
      allocate (scratch%vt3do(ntpts1))
      allocate (scratch%vt3dp(ntpts1))

      allocate (scratch%vt2da(ntpts2))
      allocate (scratch%vt2db(ntpts2))     
      allocate (scratch%vt2dc(ntpts2))
      allocate (scratch%vt2dd(ntpts2))
      allocate (scratch%vt2de(ntpts2))
      allocate (scratch%vt2df(ntpts2))

return
END SUBROUTINE alloc_scratch

!##############################################################################
Subroutine dealloc_scratch ()
   
implicit none

! Deallocate all scratch arrays

   if (allocated(scratch%scr1 ))  deallocate (scratch%scr1 )
   if (allocated(scratch%scr2 ))  deallocate (scratch%scr2 )
   if (allocated(scratch%vt3da))  deallocate (scratch%vt3da)
   if (allocated(scratch%vt3db))  deallocate (scratch%vt3db)
   if (allocated(scratch%vt3dc))  deallocate (scratch%vt3dc)
   if (allocated(scratch%vt3dd))  deallocate (scratch%vt3dd)
   if (allocated(scratch%vt3de))  deallocate (scratch%vt3de)
   if (allocated(scratch%vt3df))  deallocate (scratch%vt3df)
   if (allocated(scratch%vt3dg))  deallocate (scratch%vt3dg)
   if (allocated(scratch%vt3dh))  deallocate (scratch%vt3dh)
   if (allocated(scratch%vt3di))  deallocate (scratch%vt3di)
   if (allocated(scratch%vt3dj))  deallocate (scratch%vt3dj)
   if (allocated(scratch%vt3dk))  deallocate (scratch%vt3dk)
   if (allocated(scratch%vt3dl))  deallocate (scratch%vt3dl)
   if (allocated(scratch%vt3dm))  deallocate (scratch%vt3dm)
   if (allocated(scratch%vt3dn))  deallocate (scratch%vt3dn)
   if (allocated(scratch%vt3do))  deallocate (scratch%vt3do)
   if (allocated(scratch%vt3dp))  deallocate (scratch%vt3dp)
   if (allocated(scratch%vt2da))  deallocate (scratch%vt2da)
   if (allocated(scratch%vt2db))  deallocate (scratch%vt2db)
   if (allocated(scratch%vt2dc))  deallocate (scratch%vt2dc)
   if (allocated(scratch%vt2dd))  deallocate (scratch%vt2dd)
   if (allocated(scratch%vt2de))  deallocate (scratch%vt2de)
   if (allocated(scratch%vt2df))  deallocate (scratch%vt2df)
        
return
END SUBROUTINE dealloc_scratch

END MODULE mem_scratch
