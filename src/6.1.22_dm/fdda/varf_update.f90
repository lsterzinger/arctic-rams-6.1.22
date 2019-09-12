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

Subroutine varf_update (iswap,ifileok,initflag)

use mem_leaf
use mem_varinit
use mem_basic
use mem_grid
use mem_scratch
use node_mod
use micphys
use hdf5_utils

implicit none

integer :: ifileok,initflag,iswap

!---------------------------------------------------------------+
!    "Variable initialization"  initialization routines
!---------------------------------------------------------------+
logical :: there
integer :: iver_var,nc,iyearx,imonthx,idatex,ihourx  &
          ,nxpx,nypx,nzpx,kk,ii,jj
real :: rlatx,wlon1x,deltaxx,deltazx,dzratx,dzmaxx
character(len=7) :: cgrid
character(len=strl1) :: flnm
integer :: h5_fid, iphdf5
type (hdf5_select_type) :: mem_select,file_select
integer, dimension(HDF5_MAX_DIMS) :: file_chunks

!      Check and see what we are doing. If it is initial time, read
!        fields into regular arrays. If not, see if nudging will be done
!        on this grid if it is a nested grid.
if (ngrid > 1 .and. tnudcent+tnudtop < .001 .and. initflag == 0) return

if (nmachs .gt. 1) then
  iphdf5 = 1
else
  iphdf5 = 0
endif

! Put new fields into varinit future arrays. If iswap == 1, 
!     swap future into past first

if (iswap == 1) then

do jj=1,mmyp(ngrid)
 do ii=1,mmxp(ngrid)
  do kk=1,mmzp(ngrid)
   !The first 5 are needed for all nudged simulations
   varinit_g(ngrid)%varup(kk,ii,jj) = varinit_g(ngrid)%varuf(kk,ii,jj)
   varinit_g(ngrid)%varvp(kk,ii,jj) = varinit_g(ngrid)%varvf(kk,ii,jj)
   varinit_g(ngrid)%varpp(kk,ii,jj) = varinit_g(ngrid)%varpf(kk,ii,jj)
   varinit_g(ngrid)%vartp(kk,ii,jj) = varinit_g(ngrid)%vartf(kk,ii,jj)
   varinit_g(ngrid)%varrp(kk,ii,jj) = varinit_g(ngrid)%varrf(kk,ii,jj)
   !These next 2 are for condensate nudging from History-Varfiles
   if (nud_cond == 1) then
    varinit_g(ngrid)%varrph(kk,ii,jj) = varinit_g(ngrid)%varrfh(kk,ii,jj)
    varinit_g(ngrid)%varcph(kk,ii,jj) = varinit_g(ngrid)%varcfh(kk,ii,jj)
   endif
  enddo
 enddo
enddo

endif

write(cgrid,'(a2,i1,a3)') '-g',ngrid,'.h5'
nc=len_trim(fnames_varf(nvarffl))
flnm=fnames_varf(nvarffl)(1:nc-4)//trim(cgrid)
inquire(file=trim(flnm),exist=there)

! Gotta have grid 1...
if (.not.there .and. ngrid == 1) then
   print*
   print*,'No grid 1 varfile found: ',trim(flnm)
   print*
   stop 'no grid 1 varfile'
endif

if(there) then
   ifileok=1
else
   ifileok=0
   return
endif

! Read the varfile fields into the "future" varinit arrays. These will be 
!   swapped to the past arrays when needed.
CALL shdf5_open (flnm,'R',iphdf5,h5_fid)

! scalar vars
CALL shdf5_set_hs_select (1,'R',ngrid,mem_select,file_select,file_chunks)
CALL shdf5_irec (h5_fid,iphdf5,'version',mem_select,file_select,ivars=iver_var)
CALL shdf5_irec (h5_fid,iphdf5,'year',mem_select,file_select,ivars=iyearx)
CALL shdf5_irec (h5_fid,iphdf5,'month',mem_select,file_select,ivars=imonthx)
CALL shdf5_irec (h5_fid,iphdf5,'day',mem_select,file_select,ivars=idatex)
CALL shdf5_irec (h5_fid,iphdf5,'hour',mem_select,file_select,ivars=ihourx)
CALL shdf5_irec (h5_fid,iphdf5,'nx',mem_select,file_select,ivars=nxpx)
CALL shdf5_irec (h5_fid,iphdf5,'ny',mem_select,file_select,ivars=nypx)
CALL shdf5_irec (h5_fid,iphdf5,'nz',mem_select,file_select,ivars=nzpx)
CALL shdf5_irec (h5_fid,iphdf5,'polelat',mem_select,file_select,rvars=rlatx)
CALL shdf5_irec (h5_fid,iphdf5,'polelon',mem_select,file_select,rvars=wlon1x)
CALL shdf5_irec (h5_fid,iphdf5,'dx',mem_select,file_select,rvars=deltaxx)
CALL shdf5_irec (h5_fid,iphdf5,'dz',mem_select,file_select,rvars=deltazx)
CALL shdf5_irec (h5_fid,iphdf5,'dzrat',mem_select,file_select,rvars=dzratx)
CALL shdf5_irec (h5_fid,iphdf5,'dzmax',mem_select,file_select,rvars=dzmaxx)

if(nxp.ne.nxpx.or.  &
   nyp.ne.nypx.or.  &
   nzp.ne.nzpx.or.  &
   abs(deltax-deltaxx).gt..001.or.  &
   abs(deltaz-deltazx).gt..001.or.  &
   abs(dzrat-dzratx).gt..001.or.  & 
   abs(dzmax-dzmaxx).gt..001.or.  &
   abs(polelat-rlatx).gt..001.or.  &
   abs(polelon-wlon1x).gt..001) then
   
   print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   print*,'!!    GRID MISMATCH BETWEEN VARFILE AND NAMELIST !'
   print*,'!!          RUN IS STOPPED                       !'
   print*,'!!  File:',trim(flnm)
   print*,'!!  File, Namelist values for grid:',ngrid
   print*,'!!  nxp:',nxpx,nxp
   print*,'!!  nyp:',nypx,nyp
   print*,'!!  nzp:',nzpx,nzp
   print*,'!!  deltax:',deltaxx,deltax
   print*,'!!  deltaz:',deltazx,deltaz
   print*,'!!  dzrat:',dzratx,dzrat
   print*,'!!  dzmax:',dzmaxx,dzmax
   print*,'!!  polelat:',rlatx,polelat
   print*,'!!  polelon:',wlon1x,polelon
   PRINT*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   stop 'bad-vfile'
endif

   ! Atmos 3D vars
   CALL shdf5_set_hs_select (3,'R',ngrid,mem_select,file_select,file_chunks)
   CALL shdf5_irec (h5_fid,iphdf5,'UP',mem_select,file_select &
                   ,rvara=scratch%scr1)
   CALL unarrange (mzp,mxp,myp,scratch%scr1,varinit_g(ngrid)%varuf)
   CALL shdf5_irec (h5_fid,iphdf5,'VP',mem_select,file_select &
                   ,rvara=scratch%scr1)
   CALL unarrange (mzp,mxp,myp,scratch%scr1,varinit_g(ngrid)%varvf)
   CALL shdf5_irec (h5_fid,iphdf5,'PI',mem_select,file_select &
                   ,rvara=scratch%scr1)
   CALL unarrange (mzp,mxp,myp,scratch%scr1,varinit_g(ngrid)%varpf)
   CALL shdf5_irec (h5_fid,iphdf5,'THETA',mem_select,file_select &
                   ,rvara=scratch%scr1)
   CALL unarrange (mzp,mxp,myp,scratch%scr1,varinit_g(ngrid)%vartf)
   CALL shdf5_irec (h5_fid,iphdf5,'RV',mem_select,file_select &
                   ,rvara=scratch%scr1)
   CALL unarrange (mzp,mxp,myp,scratch%scr1,varinit_g(ngrid)%varrf)

!For Condensate nudging from History-Varfiles, RV = RTP so we load
!the RV variable into the nudging variable VARRFH
if(nud_cond == 1) then
   CALL shdf5_irec (h5_fid,iphdf5,'RV',mem_select,file_select &
                   ,rvara=scratch%scr1)
   CALL unarrange (mzp,mxp,myp,scratch%scr1,varinit_g(ngrid)%varrfh)
   CALL shdf5_irec (h5_fid,iphdf5,'COND',mem_select,file_select &
                  ,rvara=scratch%scr1)
   CALL unarrange (mzp,mxp,myp,scratch%scr1,varinit_g(ngrid)%varcfh)
endif

varinit_g(ngrid)%varrf(1:mzp,1:mxp,1:myp)=  &
           max(1.e-8,varinit_g(ngrid)%varrf(1:mzp,1:mxp,1:myp) )

!Extract soil/snow data from the varfile. Ignore other 2D fields for now.
!Do this for initialization only
if(initflag == 1 .and. iver_var == 3) then
  CALL shdf5_set_hs_select (2,'R',ngrid,mem_select,file_select,file_chunks)
  CALL shdf5_irec (h5_fid,iphdf5,'SOILMOIST1',mem_select,file_select &
                  ,rvara=leaf_g(ngrid)%soil_moist_bot)
  CALL shdf5_irec (h5_fid,iphdf5,'SOILMOIST2',mem_select,file_select &
                  ,rvara=leaf_g(ngrid)%soil_moist_top)
  CALL shdf5_irec (h5_fid,iphdf5,'SOILTEMP1', mem_select,file_select &
                  ,rvara=leaf_g(ngrid)%soil_temp_bot)
  CALL shdf5_irec (h5_fid,iphdf5,'SOILTEMP2', mem_select,file_select &
                  ,rvara=leaf_g(ngrid)%soil_temp_top)
  CALL shdf5_irec (h5_fid,iphdf5,'SNOWMASS',  mem_select,file_select &
                  ,rvara=leaf_g(ngrid)%snow_mass)
  CALL shdf5_irec (h5_fid,iphdf5,'SNOWDEPTH', mem_select,file_select &
                  ,rvara=leaf_g(ngrid)%snow_depth)    
endif

CALL shdf5_close (h5_fid)

! Find the reference state

if(initflag == 1 .and. ngrid == 1)  &
     CALL varref (mzp,mxp,myp &
         ,varinit_g(ngrid)%vartf(1,1,1) ,varinit_g(ngrid)%varpf(1,1,1)  &
         ,basic_g(ngrid)%pi0(1,1,1),     basic_g(ngrid)%th0(1,1,1)  &
         ,varinit_g(ngrid)%varrf(1,1,1), basic_g(ngrid)%dn0(1,1,1)  &
         ,basic_g(ngrid)%dn0u(1,1,1),    basic_g(ngrid)%dn0v(1,1,1)  &
         ,varinit_g(ngrid)%varuf(1,1,1), varinit_g(ngrid)%varvf(1,1,1)  &
         ,grid_g(ngrid)%topt(1,1),       grid_g(ngrid)%rtgt(1,1)  &
         ,level)

do jj=1,mmyp(ngrid)
 do ii=1,mmxp(ngrid)
  do kk=1,mmzp(ngrid)
    varinit_g(ngrid)%varpf(kk,ii,jj) =  &
        varinit_g(ngrid)%varpf(kk,ii,jj) - basic_g(ngrid)%pi0(kk,ii,jj)
  enddo
 enddo
enddo

! If this is an initialization, put data into regular arrays

if(initflag == 1 ) then
   CALL atob (mxyzp,varinit_g(ngrid)%varuf(1,1,1),basic_g(ngrid)%uc(1,1,1))
   CALL atob (mxyzp,varinit_g(ngrid)%varvf(1,1,1),basic_g(ngrid)%vc(1,1,1))
   CALL atob (mxyzp,varinit_g(ngrid)%varpf(1,1,1),basic_g(ngrid)%pc(1,1,1))
   CALL atob (mxyzp,varinit_g(ngrid)%vartf(1,1,1),basic_g(ngrid)%thp(1,1,1))
   CALL atob (mxyzp,varinit_g(ngrid)%varrf(1,1,1),basic_g(ngrid)%rtp(1,1,1))
endif

return
END SUBROUTINE varf_update

!##############################################################################
Subroutine varref (n1,n2,n3,thp,pc,pi0,th0,rtp,dn0,dn0u,dn0v,uc  &
                 ,vc,topt,rtgt,level)

use mem_grid
use ref_sounding
use mem_scratch
use rconstants
use node_mod
                 
implicit none

integer :: n1,n2,n3,level,i,j,k
real, dimension(n1,n2,n3) :: thp,pc,pi0,rtp,dn0,dn0u,dn0v,uc,vc,th0
real, dimension(n2,n3) :: topt,rtgt

! No single node has the full domain version of the atmospheric
! variables (thp, pc, pi0,...). For this reason, have each node
! calculate a 1D reference state locally for their sub-domain.
!
! Then have mainnum gather these results and determine which column
! will be used for the global 1D reference state. Then have mainnum
! scatter this information back to all of the other nodes.

!                Reference sounding is point with lowest topography
topref=1.e10
do j=1,myp
   do i=1,mxp
      if(topt(i,j).lt.topref) then
         iref=i
         jref=j
         topref=topt(i,j)
      endif
   enddo
enddo

!  Set up 1-D reference state

do k=1,mzp
   vctr2(k)=ztn(k,ngrid)*(1.-topref/ztop)+topref
enddo
CALL htint2 (mzp,thp(1,iref,jref),vctr2,mzp,vctr1,zt)
CALL htint2 (mzp,uc(1,iref,jref),vctr2,mzp,u01dn(1,ngrid),zt)
CALL htint2 (mzp,vc(1,iref,jref),vctr2,mzp,v01dn(1,ngrid),zt)
if (level >= 1) then
   CALL htint2 (mzp,rtp(1,iref,jref),vctr2,mzp,rt01dn(1,ngrid),zt)
else
   rt01dn(1:mzp,ngrid) = 0.
endif

do k = 1,mzp
   th01dn(k,ngrid) = vctr1(k) * (1. + .61 * rt01dn(k,ngrid))
enddo
u01dn(1,ngrid) = u01dn(2,ngrid)
v01dn(1,ngrid) = v01dn(2,ngrid)
rt01dn(1,ngrid) = rt01dn(2,ngrid)
th01dn(1,ngrid) = th01dn(2,ngrid)

pi01dn(1,ngrid) = pc(1,iref,jref) + g * (vctr2(1) - zt(1))  &
   / (.5 * (th01dn(1,ngrid)  &
   + thp(1,iref,jref) * (1. + .61 * rtp(1,iref,jref))))
do k = 2,mzp
  pi01dn(k,ngrid) = pi01dn(k-1,ngrid) - g / (dzm(k-1) * .5  &
     * (th01dn(k,ngrid) + th01dn(k-1,ngrid)))
enddo

do k = 1,mzp
  vctr4(k) = (pi01dn(k,ngrid) / cp) ** cpor * p00
  dn01dn(k,ngrid) = cp * vctr4(k)  &
     / (rgas * th01dn(k,ngrid) * pi01dn(k,ngrid))
enddo

! Sync up all nodes with a single global 1D reference state.
if (nmachs .gt. 1) then
  CALL establish_1d_refstate ()
endif

!        Compute 3-D reference state from 1-D reference state

CALL refs3d (mzp,mxp,myp,pi0,dn0,dn0u,dn0v,th0,topt,rtgt)

! Sync up overlap regions between sub-domains on the complete
! var list which contains the 3D reference state
if (nmachs .gt. 1) then
  CALL node_sendlbc (0)
  CALL node_getlbc (0)
endif

return
END SUBROUTINE varref

!#######################################################
! This routine will establish a single global reference
! state to be used by all nodes in a parallel run.
!
! At this point, every node has calculated their own local
! reference state. We want to consider the full domain when
! selecting the global reference state, and have this
! parallel algorithm select the same column that the 
! sequential algorithm would have picked.
!
! Do this by having mainnum gather the local ref states,
! pick one for the global ref state, and then scatter this
! result back to the other nodes. When ties occur between
! local reference states (based on local minimum elevations)
! use the one that came from first the furthest south followed
! by the furthest west.

Subroutine establish_1d_refstate ()

use mem_grid
use node_mod
use ref_sounding

implicit none

  integer, dimension(:), allocatable :: all_iref,all_jref
  real, dimension(:), allocatable :: all_topref
  real, dimension(:,:), allocatable :: all_u01dn,all_v01dn,all_rt01dn &
                                      ,all_th01dn,all_pi01dn,all_dn01dn

  integer :: inode, sel_node, sel_iref, sel_jref
  real :: sel_topref
  logical :: sel_this_node

  ! Before gathering the 1D ref states, have each node translate
  ! (iref, jref) indicies to the full domain reference.
  iref = iref + i0
  jref = jref + j0

  ! Allocate the buffers for the MPI transfers
  allocate(all_iref(nmachs)) 
  allocate(all_jref(nmachs)) 
  allocate(all_topref(nmachs)) 
  allocate(all_u01dn(mzp,nmachs)) 
  allocate(all_v01dn(mzp,nmachs)) 
  allocate(all_rt01dn(mzp,nmachs)) 
  allocate(all_th01dn(mzp,nmachs)) 
  allocate(all_pi01dn(mzp,nmachs)) 
  allocate(all_dn01dn(mzp,nmachs)) 

  ! mainnum does the work
  if (my_rams_num .eq. mainnum) then
    ! Collect the ref state data from all the nodes
    CALL par_gather_ints (iref,all_iref,1,machnum(mainnum))
    CALL par_gather_ints (jref,all_jref,1,machnum(mainnum))
    CALL par_gather_floats (topref,all_topref,1,machnum(mainnum))
    CALL par_gather_floats (u01dn(1,ngrid),all_u01dn,mzp,machnum(mainnum))
    CALL par_gather_floats (v01dn(1,ngrid),all_v01dn,mzp,machnum(mainnum))
    CALL par_gather_floats (rt01dn(1,ngrid),all_rt01dn,mzp,machnum(mainnum))
    CALL par_gather_floats (th01dn(1,ngrid),all_th01dn,mzp,machnum(mainnum))
    CALL par_gather_floats (pi01dn(1,ngrid),all_pi01dn,mzp,machnum(mainnum))
    CALL par_gather_floats (dn01dn(1,ngrid),all_dn01dn,mzp,machnum(mainnum))

    ! Find the minumum elevation that would have been found in
    ! a sequential run.
    sel_node = 1
    sel_iref = all_iref(1)
    sel_jref = all_jref(1)
    sel_topref = all_topref(1)

    do inode = 1, nmachs
      if (inode .ne. mainnum) then
        sel_this_node = .false.
  
        if (all_topref(inode) .lt. sel_topref) then
          ! This node has a lesser elevation
          sel_this_node = .true.
        else if (all_topref(inode) .eq. sel_topref) then
          ! This node is in a tie with the currently selected node
          ! Make sure we pick the same one that the sequential
          ! algorithm would pick which is first the furthest
          ! south, then the furthest west.
          if (all_jref(inode) .lt. sel_jref) then
             ! This node is further south
             sel_this_node = .true.
          else if ((all_jref(inode) .eq. sel_jref) .and. &
                   (all_iref(inode) .lt. sel_iref)) then
             ! This node is in the same row as currently selected,
             ! and further west
             sel_this_node = .true.
          endif
        endif
  
        if (sel_this_node) then
          sel_node = inode
  
          sel_iref = all_iref(inode)
          sel_jref = all_jref(inode)
          sel_topref = all_topref(inode)
        endif
      endif
    enddo

    ! Broadcast the selected state to the other nodes
    iref = sel_iref
    jref = sel_jref
    topref = sel_topref
    u01dn(1:mzp,ngrid) = all_u01dn(1:mzp,sel_node)
    v01dn(1:mzp,ngrid) = all_v01dn(1:mzp,sel_node)
    rt01dn(1:mzp,ngrid) = all_rt01dn(1:mzp,sel_node)
    th01dn(1:mzp,ngrid) = all_th01dn(1:mzp,sel_node)
    pi01dn(1:mzp,ngrid) = all_pi01dn(1:mzp,sel_node)
    dn01dn(1:mzp,ngrid) = all_dn01dn(1:mzp,sel_node)

    CALL broadcast_vfile_refstate ()

  else
    ! Send the ref state data to mainnum
    CALL par_gather_ints (iref,all_iref,1,machnum(mainnum))
    CALL par_gather_ints (jref,all_jref,1,machnum(mainnum))
    CALL par_gather_floats (topref,all_topref,1,machnum(mainnum))
    CALL par_gather_floats (u01dn(1,ngrid),all_u01dn,mzp,machnum(mainnum))
    CALL par_gather_floats (v01dn(1,ngrid),all_v01dn,mzp,machnum(mainnum))
    CALL par_gather_floats (rt01dn(1,ngrid),all_rt01dn,mzp,machnum(mainnum))
    CALL par_gather_floats (th01dn(1,ngrid),all_th01dn,mzp,machnum(mainnum))
    CALL par_gather_floats (pi01dn(1,ngrid),all_pi01dn,mzp,machnum(mainnum))
    CALL par_gather_floats (dn01dn(1,ngrid),all_dn01dn,mzp,machnum(mainnum))

    ! Received the selected stated from mainnum
    CALL broadcast_vfile_refstate ()
  endif
  
  ! free up MPI buffers
  deallocate(all_iref) 
  deallocate(all_jref)
  deallocate(all_topref)
  deallocate(all_u01dn)
  deallocate(all_v01dn)
  deallocate(all_rt01dn)
  deallocate(all_th01dn)
  deallocate(all_pi01dn)
  deallocate(all_dn01dn)

return
END SUBROUTINE establish_1d_refstate
