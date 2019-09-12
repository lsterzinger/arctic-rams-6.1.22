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

Subroutine rams_mem_alloc (proc_type)

use mem_all
use node_mod

implicit none 

integer :: proc_type,ng,nv,imean,nsc
integer, pointer :: nmzp(:),nmxp(:),nmyp(:)

! First, depending on type of process, define grid point pointers correctly...

if (proc_type == 0 .or. proc_type == 1) then
   !  This is the call for either a single processor run or
   !    for the master process
   if(iprntstmt>=1 .and. print_msg) then
      print*,''
      print*,'MEMORY ALLOCATION ON MASTER NODE'
   endif
   nmzp => nnzp
   nmxp => nnxp
   nmyp => nnyp
elseif (proc_type == 2) then
   !  This is the call for an initial compute node process
   nmzp => mmzp
   nmxp => mmxp
   nmyp => mmyp
endif

!  If we are doing time-averaging for output, set flag ...
imean=0
if (avgtim /= 0.) imean=1

! Allocate universal variable tables

allocate (num_var(maxgrds))
allocate (vtab_r(maxvars,maxgrds))

num_var(1:maxgrds)=0
nvgrids=ngrids

! Allocate scalar table  

allocate(num_scalar(maxgrds))
allocate(scalar_tab(maxsclr,maxgrds))
num_scalar(1:maxgrds)=0

! Allocate Basic variables data type
if(iprntstmt>=1 .and. print_msg .and. proc_type /= 2) &
   print*,'start basic alloc'
allocate(basic_g(ngrids),basicm_g(ngrids))
do ng=1,ngrids
   CALL dealloc_basic (basic_g(ng))
   CALL dealloc_basic (basicm_g(ng))
   CALL alloc_basic (basic_g(ng),nmzp(ng),nmxp(ng),nmyp(ng)) 
   if (imean == 1) then  
      CALL alloc_basic (basicm_g(ng),nmzp(ng),nmxp(ng),nmyp(ng))
   elseif (imean == 0) then
      CALL alloc_basic (basicm_g(ng),1,1,1)
   endif
   
   CALL filltab_basic (basic_g(ng),basicm_g(ng),imean  &
          ,nmzp(ng),nmxp(ng),nmyp(ng),ng) 
enddo

! Allocate Cuparm variables data type
if(iprntstmt>=1 .and. print_msg .and. proc_type /= 2) &
   print*,'start cuparm alloc'
allocate(cuparm_g(ngrids),cuparmm_g(ngrids))
do ng=1,ngrids
   CALL dealloc_cuparm (cuparm_g(ng))
   CALL dealloc_cuparm (cuparmm_g(ng))
   CALL alloc_cuparm (cuparm_g(ng),nmzp(ng),nmxp(ng),nmyp(ng),ng) 
   if (imean == 1) then  
      CALL alloc_cuparm (cuparmm_g(ng),nmzp(ng),nmxp(ng),nmyp(ng),ng)
   elseif (imean == 0) then
      CALL alloc_cuparm (cuparmm_g(ng),1,1,1,ng)
   endif
   
   CALL filltab_cuparm (cuparm_g(ng),cuparmm_g(ng),imean  &
          ,nmzp(ng),nmxp(ng),nmyp(ng),ng) 
enddo

! Allocate Leaf type

if(iprntstmt>=1 .and. print_msg .and. proc_type /= 2) &
   print*,'start leaf alloc'
allocate(leaf_g(ngrids),leafm_g(ngrids))
do ng=1,ngrids
   CALL dealloc_leaf (leaf_g(ng))
   CALL dealloc_leaf (leafm_g(ng))
   CALL alloc_leaf (leaf_g(ng),nmxp(ng),nmyp(ng)  &
       ,nzg,nzs,npatch) 
   if (imean == 1) then  
      CALL alloc_leaf (leafm_g(ng),nmxp(ng),nmyp(ng)  &
          ,nzg,nzs,npatch)
   elseif (imean == 0) then
      CALL alloc_leaf (leafm_g(ng),1,1,1,1,1)
   endif
   
   CALL filltab_leaf (leaf_g(ng),leafm_g(ng),imean  &
          ,nmxp(ng),nmyp(ng),nzg,nzs,npatch,ng) 
enddo

! Allocate SiB-RAMS
if(iprntstmt>=1 .and. print_msg .and. proc_type /= 2) &
   print*,'start Sib alloc'
allocate(sib_g(ngrids),sibm_g(ngrids))
do ng=1,ngrids
   CALL dealloc_sib (sib_g(ng))
   CALL dealloc_sib (sibm_g(ng))
   CALL alloc_sib (sib_g(ng),nmzp(ng),nmxp(ng),nmyp(ng),npatch)
   if (imean==1) then
      CALL alloc_sib (sibm_g(ng),nmzp(ng),nmxp(ng),nmyp(ng),npatch)
   elseif (imean==0) then
      CALL alloc_sib (sibm_g(ng),1,1,1,1)
   endif

   CALL filltab_sib (sib_g(ng),sibm_g(ng),imean  &
          ,nmzp(ng),nmxp(ng),nmyp(ng),npatch,ng)
enddo

! Allocate Micro variables data type
if(iprntstmt>=1 .and. print_msg .and. proc_type /= 2)print*,'start micro alloc'
allocate(micro_g(ngrids),microm_g(ngrids))
do ng=1,ngrids
   CALL dealloc_micro (micro_g(ng))
   CALL dealloc_micro (microm_g(ng))
   CALL alloc_micro (micro_g(ng),nmzp(ng),nmxp(ng),nmyp(ng)) 
   if (imean == 1) then  
      CALL alloc_micro (microm_g(ng),nmzp(ng),nmxp(ng),nmyp(ng))
   elseif (imean == 0) then
      CALL alloc_micro (microm_g(ng),1,1,1)
   endif
   
   CALL filltab_micro (micro_g(ng),microm_g(ng),imean  &
          ,nmzp(ng),nmxp(ng),nmyp(ng),ng) 
enddo

! Allocate radiate variables data type
if(iprntstmt>=1 .and. print_msg .and. proc_type /= 2) &
   print*,'start radiate alloc'
allocate(radiate_g(ngrids),radiatem_g(ngrids))
do ng=1,ngrids
   CALL dealloc_radiate (radiate_g(ng))
   CALL dealloc_radiate (radiatem_g(ng))
   CALL alloc_radiate (radiate_g(ng),nmzp(ng),nmxp(ng),nmyp(ng)) 
   if (imean == 1) then  
      CALL alloc_radiate (radiatem_g(ng),nmzp(ng),nmxp(ng),nmyp(ng))
   elseif (imean == 0) then
      CALL alloc_radiate (radiatem_g(ng),1,1,1)
   endif
   
   CALL filltab_radiate (radiate_g(ng),radiatem_g(ng),imean  &
          ,nmzp(ng),nmxp(ng),nmyp(ng),ng) 
enddo

! Allocate turb variables data type
if(iprntstmt>=1 .and. print_msg .and. proc_type /= 2) &
   print*,'start turb alloc'
allocate(turb_g(ngrids),turbm_g(ngrids))
do ng=1,ngrids
   CALL dealloc_turb (turb_g(ng))
   CALL dealloc_turb (turbm_g(ng))
   CALL alloc_turb (turb_g(ng),nmzp(ng),nmxp(ng),nmyp(ng),ngrids) 
   if (imean == 1) then  
      CALL alloc_turb (turbm_g(ng),nmzp(ng),nmxp(ng),nmyp(ng),ngrids)
   elseif (imean == 0) then
      CALL alloc_turb (turbm_g(ng),1,1,1,ngrids)
   endif
   
   CALL filltab_turb (turb_g(ng),turbm_g(ng),imean  &
          ,nmzp(ng),nmxp(ng),nmyp(ng),ng) 
enddo

! Allocate varinit variables data type. 
!    These do not need "mean" type ever.
if(iprntstmt>=1 .and. print_msg .and. proc_type /= 2) &
   print*,'start varinit alloc'
allocate(varinit_g(ngrids),varinitm_g(ngrids))
do ng=1,ngrids
   CALL dealloc_varinit (varinit_g(ng))
   CALL dealloc_varinit (varinitm_g(ng))
   CALL alloc_varinit (varinit_g(ng),nmzp(ng),nmxp(ng),nmyp(ng)) 
   CALL alloc_varinit (varinitm_g(ng),1,1,1)
   
   CALL filltab_varinit (varinit_g(ng),varinitm_g(ng),0  &
          ,nmzp(ng),nmxp(ng),nmyp(ng),ng) 
enddo

! Allocate oda variables data type. 
!    These do not need "mean" type ever.
if(iprntstmt>=1 .and. print_msg .and. proc_type /= 2) &
   print*,'start oda alloc'
allocate(oda_g(ngrids),odam_g(ngrids))
do ng=1,ngrids
   CALL dealloc_oda (oda_g(ng))
   CALL dealloc_oda (odam_g(ng))
   CALL alloc_oda (oda_g(ng),nmzp(ng),nmxp(ng),nmyp(ng),proc_type) 
   CALL alloc_oda (odam_g(ng),1,1,1,proc_type)
   
   CALL filltab_oda (oda_g(ng),odam_g(ng),0  &
          ,nmzp(ng),nmxp(ng),nmyp(ng),ng) 
enddo

! Allocate grid variables data type. 

if(iprntstmt>=1 .and. print_msg .and. proc_type /= 2) &
   print*,'start grid alloc'
allocate(grid_g(ngrids),gridm_g(ngrids))
do ng=1,ngrids
   CALL dealloc_grid (grid_g(ng))
   CALL dealloc_grid (gridm_g(ng))
   CALL alloc_grid (grid_g(ng),nmxp(ng),nmyp(ng)) 
   if (imean == 1) then
      CALL alloc_grid (gridm_g(ng),nmxp(ng),nmyp(ng)) 
   elseif (imean == 0) then
      CALL alloc_grid (gridm_g(ng),1,1)
   endif

   CALL filltab_grid (grid_g(ng),gridm_g(ng),imean,nmxp(ng),nmyp(ng),ng) 
enddo


! Allocate any added Scalar types

if(iprntstmt>=1 .and. print_msg .and. proc_type /= 2) &
   print*,'start scalar alloc'
allocate(tracer_g(itracer,ngrids),tracerm_g(itracer,ngrids))
do ng=1,ngrids
   CALL dealloc_tracer (tracer_g(:,ng),ng)
   CALL dealloc_tracer (tracerm_g(:,ng),ng)
   CALL alloc_tracer (tracer_g(:,ng),nmzp(ng),nmxp(ng),nmyp(ng),ng)
   if (imean == 1) then  
    CALL alloc_tracer (tracerm_g(:,ng),nmzp(ng),nmxp(ng),nmyp(ng),ng)
   elseif (imean == 0) then
    CALL alloc_tracer (tracerm_g(:,ng),1,1,1,ng)
   endif
      
   CALL filltab_tracer (tracer_g(:,ng),tracerm_g(:,ng),imean  &
                       ,nmzp(ng),nmxp(ng),nmyp(ng),ng)
enddo

! Allocate Tendency data type,  filltab_tendency is responsible 
!   for filling the main atmospheric model variables in the scalar table,
!   so make sure to call any routines that define scalar variables first.

! Assuming same scalars on all grids!!!!!

if(iprntstmt>=1 .and. print_msg .and. proc_type /= 2) &
   print*,'start tendency alloc'
CALL dealloc_tend (ngrids)
CALL alloc_tend (nmzp,nmxp,nmyp,ngrids,proc_type)
do ng=1,ngrids
   CALL filltab_tend (basic_g(ng),micro_g(ng),turb_g(ng),sib_g(ng) &
      ,tracer_g(:,ng),ng)
enddo

! Allocate Scratch data type, This also fills the max's that are needed
!    by nesting stuff.

if(iprntstmt>=1 .and. print_msg .and. proc_type /= 2) &
   print*,'start scratch alloc'
CALL dealloc_scratch ()
CALL alloc_scratch (nmzp,nmxp,nmyp,nnzp,nnxp,nnyp,ngrids  &
                  ,nzg,nzs,npatch,proc_type,maxnxp,maxnyp,maxnzp)

! Allocate nested boundary interpolation arrays. All grids will be allocated.

if(iprntstmt>=1 .and. print_msg .and. proc_type /= 2) &
   print*,'start nestb alloc'
!if (proc_type == 0 .or. proc_type == 2) then
! We'll allocate this all the time for now, even though the master process in 
!   a parallel run doesn't really need this. The arrays appear in a call statement
!   during the initialization, but are never used.
   do ng=1,ngrids
      if(nxtnest(ng) == 0 ) then
         CALL alloc_nestb (ng,1,1,1)
      else
         CALL alloc_nestb (ng,nnxp(ng),nnyp(ng),nnzp(ng))
      endif
   enddo

! Set "Lite" variable flags according to namelist input LITE_VARS.

CALL lite_varset ()


! Set ALL variables in the vtab_r variable table to zero by default. These
!  are variables processed in the filltab_* routines with a call to vtables2.
!  This does NOT include scratch arrays, tendencies, or mean arrays.

do ng = 1, ngrids
   do nv = 1,num_var(ng)
      if(iprntstmt>=1 .and. print_msg .and. proc_type /= 2) &
         print*,'Zeroing out array:',ng,nv,vtab_r(nv,ng)%name
      CALL azero (vtab_r(nv,ng)%npts, vtab_r(nv,ng)%var_p)
   enddo
enddo

if(iprntstmt>=1 .and. print_msg .and. proc_type /= 2) &
   print*,'end alloc'

return
END SUBROUTINE rams_mem_alloc

