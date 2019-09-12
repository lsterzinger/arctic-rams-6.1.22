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

Subroutine par_decomp (nxp,nyp  &
   ,nodes,work,workrow,nblocks,workload,workblock,workcol  &
   ,jrows,jrow,ixb,ixe,iyb,iye)

use grid_dims

implicit none

integer :: nxp,nyp,nodes
real :: work(nxp,*),workrow(*),workload(*),workblock(*),workcol(*)
integer :: jrows(*),jrow(*),nblocks(*),ixb(*),ixe(*),iyb(*),iye(*)

! default relspeed = 1.0 for nodes of uniform speed.
real :: relspeed(maxmach)
data relspeed/maxmach*1./

integer :: inode,i,j,islab,nslabs,min_blocks,nbigslabs,iblock &
          ,jnode,knode
real :: anodes,aslabs,totspeed,workdom,workaccum,worksofar &
       ,slabspeed,workslab

! This routine decomposes grid domains of size (nnxp,nnyp) into a number,
! specified by nodes, of rectangular subdomains.  The convention is followed
! that any internal boundaries (between subdomains) that are parallel to
! the x-axis run continuously across the full domain, while boundaries
! parallel to the y-axis may or may not run the full distance across the
! domain.  For convenience, regions of the domain bounded by adjacent
! east-west internal boundaries are termed "slabs", while smaller divisions
! within each slab are termed "blocks".  Each block is required to have
! a minimum dimension of 6 by 6 grid cells.  If this cannot be satisfied
! with the given input parameters, the routine stops.

! Estimate the number of slabs to be used (aslabs), and compute a final
! nearest integer value (nslabs) which is limited to allowable values.
! Zero out array for accumulating number of columns for each node.

   anodes = float(nodes)
   aslabs = sqrt(anodes * float(nyp) / float(nxp))
   nslabs = min(nodes,max(1,nint(aslabs)))

!          print*, 'nslabs',nslabs

   totspeed = 0.
   do inode = 1,nodes
      ixe(inode) = 0
      totspeed = totspeed + relspeed(inode)
   enddo

!          print*, 'totspeed',totspeed

! Compute total work load over each row and over entire domain.

   workdom = 0.
   do j = 1,nyp
      workrow(j) = 0.
      do i = 1,nxp
         workrow(j) = workrow(j) + work(i,j)
      enddo
      workdom = workdom + workrow(j)

!          print*, 'j,workdom,workrow(j)',j,workdom,workrow(j)

   enddo
   workrow(2) = workrow(2) + workrow(1)
   workrow(nyp-1) = workrow(nyp-1) + workrow(nyp)

! Determine number of blocks and the average workload for each slab.

   min_blocks = nodes / nslabs
   nbigslabs = nodes - min_blocks * nslabs
   inode = 0
   do islab = 1,nslabs
      workload(islab) = 0.
      nblocks(islab) = min_blocks
      if (islab .le. nbigslabs) nblocks(islab) = min_blocks + 1
      do iblock = 1,nblocks(islab)
         inode = inode + 1
         workload(islab) = workload(islab)  &
            + workdom * relspeed(inode) / totspeed

!           print*, 'islab,iblock,workload(islab),workdom,inode'
!           print*,  islab,iblock,workload(islab),workdom,inode

      enddo
   enddo

! Assign all j-rows to their respective slabs in a way that balances the work
! load among slabs according to their respective numbers of nodes (blocks).
! The array jrows counts the number of rows in each slab, and the array
! jrow is the index of the southernmost row in each slab.

   do islab = 1,nslabs
      jrows(islab) = 0
   enddo

   workaccum = 0.
   worksofar = 0.
   islab = 0

   do j = 2,nyp-1
      workaccum = workaccum + workrow(j)
      if (workaccum - .5 * workrow(j) .gt. worksofar .and.  &
              islab .lt. nslabs) then
         islab = islab + 1
         jrow(islab) = j
         worksofar = worksofar + workload(islab)
      endif
      jrows(islab) = jrows(islab) + 1
   enddo

   inode = 0
   jnode = 0
   knode = 0
   do islab = 1,nslabs

! Compute the total work load for each slab and for each i-column in the
! slab.

      slabspeed = 0.
      workslab = 0.
      do i = 1,nxp
         workcol(i) = 0.
         do j = jrow(islab),jrow(islab)+jrows(islab)-1
            workcol(i) = workcol(i) + work(i,j)
         enddo
         workslab = workslab + workcol(i)
      enddo
      workcol(2) = workcol(2) + workcol(1)
      workcol(nxp-1) = workcol(nxp-1) + workcol(nxp)

! Determine average workload for each block.

      do iblock = 1,nblocks(islab)
         jnode = jnode + 1
         slabspeed = slabspeed + relspeed(jnode)

!           print*, 'r1:iblock,jnode,slabspeed,relspeed(jnode)'
!           print*,     iblock,jnode,slabspeed,relspeed(jnode)

      enddo
      do iblock = 1,nblocks(islab)
         knode = knode + 1
         workblock(iblock) = workslab  &
            * relspeed(knode) / slabspeed

!       print*, 'islab,iblock,workblock,workslab,relspeed,slabspeed'
!       print*, islab,iblock,workblock(iblock),workslab,relspeed(knode)
!     +       ,slabspeed
!       print*, 'knode',knode

      enddo

! Assign the i-columns of each slab to their respective blocks in a way that
! balances the work load among the blocks.  The array ncols counts the number
! of i-columns on each node, and the array ncol is the index of the
! westernmost i-column on each node.

      workaccum = 0.
      worksofar = 0.

      iblock = 0
      do i = 2,nxp-1
         workaccum = workaccum + workcol(i)

!        print*, 'islab',islab
!        print*, 'i,workaccum,workcol(i),worksofar,iblock,nblocks'
!        print*, i,workaccum,workcol(i),worksofar,iblock,nblocks(islab)

         if (workaccum - .5 * workcol(i) .gt. worksofar .and.  &
                iblock .lt. nblocks(islab)) then
            iblock = iblock + 1

            !Defining node variables here
            !-------------------------------------------
            inode = inode + 1
            iyb(inode) = jrow(islab)
            ixb(inode) = i
            iye(inode) = iyb(inode) + jrows(islab) - 1
            !-------------------------------------------

            worksofar = worksofar + workblock(iblock)
         endif
         ixe(inode) = ixe(inode) + 1
      enddo
   enddo

   !Defining node variables here                                    
   !----------------------------------------------------------------
   do jnode = 1,nodes
      ixe(jnode) = ixb(jnode) + ixe(jnode) - 1
      !print*, 'jnode,ixb,ixe',jnode,ixb(jnode),ixe(jnode)  &
      !,(ixe(jnode)-ixb(jnode)+1),(iye(jnode)-iyb(jnode)+1) &
      !,(ixe(jnode)-ixb(jnode)+1)*(iye(jnode)-iyb(jnode)+1)
   enddo
   !----------------------------------------------------------------

! Check to make sure that each subdomain has at least 2 interior 
! rows and columns.

   do jnode = 1,nodes
      if (iye(jnode) - iyb(jnode) .lt. 1 .or.  &
          ixe(jnode) - ixb(jnode) .lt. 1) then
         print*, 'grid:',nxp,nyp,'  subdomain too small on node ',jnode
         print*, '(ixb,ixe,iyb,iye) = '  &
            ,ixb(jnode),ixe(jnode),iyb(jnode),iye(jnode)
         stop 'small_nodes'
      endif
   enddo

return
END SUBROUTINE par_decomp

!##############################################################################
Subroutine par_est_time (nxp,nyp,work)

implicit none

integer :: nxp,nyp
real :: work(nxp,nyp)
integer :: i,j
real :: bfact

! Sample routine to fill work elements with values proportional to the time
! required to perform model operations.

do j = 2,nyp-1
   do i = 2,nxp-1
      work(i,j) = 1.
   enddo
enddo

! Fill real boundaries with .2 of interior points

bfact=.2
do j = 1,nyp
   work(1,j) = bfact * work(2,j)
   work(nxp,j) = bfact * work(nxp-1,j)
enddo

do i = 1,nxp
   work(i,1) = bfact * work(i,2)
   work(i,nyp) = bfact * work(i,nyp-1)
enddo

return
END SUBROUTINE par_est_time

!##############################################################################
Subroutine par_node_paths (maxgrds,ngrids,nxpmax,nypmax,nnxp,nnyp  &
   ,nxtnest,maxmach,nodes,node_id,nxbeg,nxend,nybeg,nyend  &
   ,ipm,jpm,ixb,ixe,iyb,iye,ipaths,igetpaths,ibnd,jbnd)

use cyclic_mod

implicit none

integer :: maxgrds,ngrids,nxpmax,nypmax,maxmach,nodes

integer :: nnxp(*),nnyp(*),nxtnest(*),node_id(*)  &
     ,nxbeg(maxmach,*),nxend(maxmach,*),nybeg(maxmach,*)  &
     ,nyend(maxmach,*),ipm(nxpmax,*),jpm(nypmax,*),ixb(maxmach,*)  &
     ,ixe(maxmach,*),iyb(maxmach,*),iye(maxmach,*)  &
     ,ipaths(5,7,maxgrds,maxmach,maxmach)  &
     ,igetpaths(6,maxgrds,maxmach,maxmach)

integer :: ngr,isend_type,idn,isn,i,j,ibnd,jbnd,info,nnn &
          ,id,jd,is,js,is_t(2),js_t(2),is_u(2),js_u(2),is_v(2),js_v(2) &
          ,indt(2),indu(2),indv(2),nxp,nyp,avgflg_t(2),avgflg_u(2) &
          ,avgflg_v(2),dir,i1d,i2d,j1d,j2d,i1s,i2s,j1s,j2s

! if using cyclic boundary conditions, allocate ipaths_cyc array

if (ibnd .eq. 2 .or. jbnd .eq. 2) then
   CALL ipaths_cyc_alloc (nnxp(1),nnyp(1),ibnd,jbnd,maxmach)
endif

! Zero out ipaths array

do ngr=1,ngrids
   do isend_type = 1,6
      do idn=1,nodes
         do isn = 1,nodes
            igetpaths(isend_type,ngr,isn,idn) = 0
            do info = 1,5
               ipaths(info,isend_type,ngr,idn,isn) = 0
            enddo
         enddo
      enddo
   enddo
enddo


do ngr=1,ngrids
   nnn = nxtnest(ngr)

   do idn=1,nodes

      do isn = 1,nodes

         if (isn .eq. idn) go to 6

! Long timestep overlap regions

         if(nxbeg(idn,ngr) .gt. ixe(isn,ngr) .or.  &
            nybeg(idn,ngr) .gt. iye(isn,ngr) .or.  &
            nxend(idn,ngr) .lt. ixb(isn,ngr) .or.  &
            nyend(idn,ngr) .lt. iyb(isn,ngr)) go to 6
         igetpaths(1,ngr,isn,idn) = node_id(isn)

         ipaths(1,1,ngr,idn,isn)=max(ixb(isn,ngr),nxbeg(idn,ngr))
         ipaths(2,1,ngr,idn,isn)=min(ixe(isn,ngr),nxend(idn,ngr))
         ipaths(3,1,ngr,idn,isn)=max(iyb(isn,ngr),nybeg(idn,ngr))
         ipaths(4,1,ngr,idn,isn)=min(iye(isn,ngr),nyend(idn,ngr))
         ipaths(5,1,ngr,idn,isn)=node_id(idn)

! Expand ipaths to include [coarse] grid extern boundary points.

            if (ipaths(1,1,ngr,idn,isn) .eq. 2)  &
                ipaths(1,1,ngr,idn,isn) = 1
            if (ipaths(2,1,ngr,idn,isn) .eq. nnxp(ngr)-1)  &
                ipaths(2,1,ngr,idn,isn) = nnxp(ngr)
            if (ipaths(3,1,ngr,idn,isn) .eq. 2)  &
                ipaths(3,1,ngr,idn,isn) = 1
            if (ipaths(4,1,ngr,idn,isn) .eq. nnyp(ngr)-1)  &
                ipaths(4,1,ngr,idn,isn) = nnyp(ngr)

! Small timestep overlap regions for u

         if (ixb(idn,ngr)-1 .eq. ixe(isn,ngr) .and.  &
             iye(idn,ngr)   .ge. iyb(isn,ngr) .and.  &
             iyb(idn,ngr)   .le. iye(isn,ngr)) then

            igetpaths(2,ngr,isn,idn) = node_id(isn)

            ipaths(1,2,ngr,idn,isn)=ixe(isn,ngr)
            ipaths(2,2,ngr,idn,isn)=ixe(isn,ngr)
            ipaths(3,2,ngr,idn,isn)=max(iyb(isn,ngr),iyb(idn,ngr))
            ipaths(4,2,ngr,idn,isn)=min(iye(isn,ngr),iye(idn,ngr))
            ipaths(5,2,ngr,idn,isn)=node_id(idn)

         endif

! Small timestep overlap regions for v

         if (iyb(idn,ngr)-1 .eq. iye(isn,ngr) .and.  &
             ixe(idn,ngr)   .ge. ixb(isn,ngr) .and.  &
             ixb(idn,ngr)   .le. ixe(isn,ngr)) then

            igetpaths(3,ngr,isn,idn) = node_id(isn)

            ipaths(1,3,ngr,idn,isn)=max(ixb(isn,ngr),ixb(idn,ngr))
            ipaths(2,3,ngr,idn,isn)=min(ixe(isn,ngr),ixe(idn,ngr))
            ipaths(3,3,ngr,idn,isn)=iye(isn,ngr)
            ipaths(4,3,ngr,idn,isn)=iye(isn,ngr)
            ipaths(5,3,ngr,idn,isn)=node_id(idn)

         endif

! Small timestep overlap regions for pi'

         if (ixe(idn,ngr)+1 .eq. ixb(isn,ngr) .and.  &
             iye(idn,ngr)   .ge. iyb(isn,ngr) .and.  &
             iyb(idn,ngr)   .le. iye(isn,ngr)) then

            igetpaths(4,ngr,isn,idn) = node_id(isn)

            ipaths(1,4,ngr,idn,isn)=ixb(isn,ngr)
            ipaths(2,4,ngr,idn,isn)=ixb(isn,ngr)
            ipaths(3,4,ngr,idn,isn)=max(iyb(isn,ngr),iyb(idn,ngr))
            ipaths(4,4,ngr,idn,isn)=min(iye(isn,ngr),iye(idn,ngr))
            ipaths(5,4,ngr,idn,isn)=node_id(idn)

         elseif (iye(idn,ngr)+1 .eq. iyb(isn,ngr) .and.  &
                 ixe(idn,ngr)   .ge. ixb(isn,ngr) .and.  &
                 ixb(idn,ngr)   .le. ixe(isn,ngr)) then

            igetpaths(4,ngr,isn,idn) = node_id(isn)

            ipaths(1,4,ngr,idn,isn)=max(ixb(isn,ngr),ixb(idn,ngr))
            ipaths(2,4,ngr,idn,isn)=min(ixe(isn,ngr),ixe(idn,ngr))
            ipaths(3,4,ngr,idn,isn)=iyb(isn,ngr)
            ipaths(4,4,ngr,idn,isn)=iyb(isn,ngr)
            ipaths(5,4,ngr,idn,isn)=node_id(idn)

         endif

6             continue

! Coarse grid to fine grid interpolation communication

         if (nnn .eq. 0) go to 26

         if(ipm(ixb(idn,ngr)-1,ngr)-2 .gt. ixe(isn,nnn).or.  &
            jpm(iyb(idn,ngr)-1,ngr)-2 .gt. iye(isn,nnn).or.  &
            ipm(ixe(idn,ngr)+1,ngr)+1 .lt. ixb(isn,nnn).or.  &
            jpm(iye(idn,ngr)+1,ngr)+1 .lt. iyb(isn,nnn)) go to 16
         igetpaths(5,ngr,isn,idn) = node_id(isn)

         ipaths(1,5,ngr,idn,isn) = max(ixb(isn,nnn)  &
            ,ipm(ixb(idn,ngr)-1,ngr)-2)
         ipaths(2,5,ngr,idn,isn) = min(ixe(isn,nnn)  &
            ,ipm(ixe(idn,ngr)+1,ngr)+1)
         ipaths(3,5,ngr,idn,isn) = max(iyb(isn,nnn)  &
            ,jpm(iyb(idn,ngr)-1,ngr)-2)
         ipaths(4,5,ngr,idn,isn) = min(iye(isn,nnn)  &
            ,jpm(iye(idn,ngr)+1,ngr)+1)
         ipaths(5,5,ngr,idn,isn) = node_id(idn)

16            continue

! Fine grid to coarse grid averaging communication

! rtimh - micphys: in the following lines: change ixb, ixe, iyb, and iye
! to nxbeg, nxend, nybeg, and nyend on coarse grid destination nodes
! in order to keep internal CG boundaries up to date, since they have
! already been updated by their CG neighboring nodes??

         if(nxbeg(idn,nnn) .gt. ipm(ixe(isn,ngr),ngr).or.  &
            nybeg(idn,nnn) .gt. jpm(iye(isn,ngr),ngr).or.  &
            nxend(idn,nnn) .lt. ipm(ixb(isn,ngr),ngr).or.  &
            nyend(idn,nnn) .lt. jpm(iyb(isn,ngr),ngr)) go to 26
         igetpaths(6,ngr,isn,idn) = node_id(isn)

         ipaths(1,6,ngr,idn,isn) = max(ipm(ixb(isn,ngr),ngr)  &
            ,nxbeg(idn,nnn))
         ipaths(2,6,ngr,idn,isn) = min(ipm(ixe(isn,ngr),ngr)  &
            ,nxend(idn,nnn))
         ipaths(3,6,ngr,idn,isn) = max(jpm(iyb(isn,ngr),ngr)  &
            ,nybeg(idn,nnn))
         ipaths(4,6,ngr,idn,isn) = min(jpm(iye(isn,ngr),ngr)  &
            ,nyend(idn,nnn))
         ipaths(5,6,ngr,idn,isn) = node_id(idn)

! A second index value of 7 of the ipaths array is used to determine
! the loop limits in fdbackp for averaging the fm over the overlap
! between the cm node and fm node, rather than always over the full
! fm node.  It is not used for actually sending stuff.  The
! ipaths(*,6,*,*,*) part of the array is still used for sending the
! block of averaged cm points from the fm node to the cm node.

         do i = ixb(isn,ngr),ixe(isn,ngr)
            if (ipm(i,ngr) .eq. ipaths(1,6,ngr,idn,isn)) then
               ipaths(1,7,ngr,idn,isn) = i
               go to 21
            endif
         enddo
21            continue

         do i = ixe(isn,ngr),ixb(isn,ngr),-1
            if (ipm(i,ngr) .eq. ipaths(2,6,ngr,idn,isn)) then
               ipaths(2,7,ngr,idn,isn) = i
               go to 22
            endif
         enddo
22            continue

         do j = iyb(isn,ngr),iye(isn,ngr)
            if (jpm(j,ngr) .eq. ipaths(3,6,ngr,idn,isn)) then
               ipaths(3,7,ngr,idn,isn) = j
               go to 23
            endif
         enddo
23            continue

         do j = iye(isn,ngr),iyb(isn,ngr),-1
            if (jpm(j,ngr) .eq. ipaths(4,6,ngr,idn,isn)) then
               ipaths(4,7,ngr,idn,isn) = j
               go to 24
            endif
         enddo
24            continue

26            continue

      enddo
   enddo
enddo

! Cyclic boundary conditions (grid 1 only)

indt = 0
indu = 0
indv = 0

nxp = nnxp(1)
nyp = nnyp(1)

! Loop over all destination ij columns

do jd = 1,nyp
   do id = 1,nxp
      if (id <= 2 .or. id >= nxp-1 .or. jd <= 2 .or. jd >= nyp-1) then

         is_t = 0
         js_t = 0
         is_u = 0
         js_u = 0
         is_v = 0
         js_v = 0
         avgflg_t = 0
         avgflg_u = 0
         avgflg_v = 0

         ! Record source i,j for t, u, and v in the i direction
         ! (element 1 of is,js arrays)
         if (ibnd == 2) then
! t
            if (id <= 2) then
               is_t(1) = id+nxp-3
               js_t(1) = jd
            elseif (id >= nxp-1) then
               is_t(1) = id-nxp+3
               js_t(1) = jd
            endif
! u
            if (id == 1) then
               is_u(1) = id+nxp-3
               js_u(1) = jd
            elseif (id == nxp-1) then
               is_u(1) = id-nxp+3
               js_u(1) = jd
            endif
! v
            if (id <= 2) then
               is_v(1) = id+nxp-3
               js_v(1) = jd
            elseif (id >= nxp-1) then
               is_v(1) = id-nxp+3
               js_v(1) = jd
            endif

         endif
          
         ! Record source i,j for t, u, and v in the j direction
         ! (element 2 of is,js arrays)
         if (jbnd == 2) then
!t
            if (jd <= 2) then
               is_t(2) = id
               js_t(2) = jd+nyp-3
            elseif (jd >= nyp-1) then
               is_t(2) = id
               js_t(2) = jd-nyp+3
            endif
!u
            if (jd <= 2) then
               is_u(2) = id
               js_u(2) = jd+nyp-3
            elseif (jd >= nyp-1) then
               is_u(2) = id
               js_u(2) = jd-nyp+3
            endif
!v
            if (jd == 1) then
               is_v(2) = id
               js_v(2) = jd+nyp-3
            elseif (jd == nyp-1) then
               is_v(2) = id
               js_v(2) = jd-nyp+3
            endif

         endif

         ! Do averaging if on the "2" or "n-1" points
         if ((id .eq. 2) .or. (id .eq. nxp-1)) then
            avgflg_t(1) = 1
            avgflg_v(1) = 1
         endif
         if ((jd .eq. 2) .or. (jd .eq. nyp-1)) then
            avgflg_t(2) = 1
            avgflg_u(2) = 1
         endif
                  
! Loop over all destination nodes and check if current ij destination column 
! is in it

         do idn = 1,nodes
          ! grab end points (compute columns)
          i1d = ixb(idn,1)
          i2d = ixe(idn,1)
          j1d = iyb(idn,1)
          j2d = iye(idn,1)
          ! include domain boundaries
          if (i1d .eq. 2) i1d = 1
          if (i2d .eq. nxp-1) i2d = nxp
          if (j1d .eq. 2) j1d = 1
          if (j2d .eq. nyp-1) j2d = nyp

          if (id >= i1d .and. id <= i2d .and. jd >= j1d .and. jd <= j2d) then
                
! Loop over all source nodes and check if current ij source column is in it
! The t paths for each corresponding direction will always contain source nodes,
! however the u and v paths may not. (Eg, dir == 1, u paths only exist for 
! id == 1 and nxp-1 (and not for id == 2 and nxp).

           do dir = 1,2
            ! dir == 1 --> along the x-axis ("i" direction)
            ! dir == 2 --> along the y-axis ("j" direction)
            is = is_t(dir)
            js = js_t(dir)

            do isn = 1,nodes
             ! grab end points (compute columns)
             i1s = ixb(isn,1)
             i2s = ixe(isn,1)
             j1s = iyb(isn,1)
             j2s = iye(isn,1)
             ! include domain boundaries
             if (i1s .eq. 2) i1s = 1
             if (i2s .eq. nxp-1) i2s = nxp
             if (j1s .eq. 2) j1s = 1
             if (j2s .eq. nyp-1) j2s = nyp

             if(is>=i1s .and. is<=i2s .and. js>=j1s .and. js<=j2s) then
              indt(dir) = indt(dir) + 1
              ipathst_cyc(1,indt(dir),dir) = isn           ! source node #
              ipathst_cyc(2,indt(dir),dir) = is_t(dir)     ! source node i
              ipathst_cyc(3,indt(dir),dir) = js_t(dir)     ! source node j
              ipathst_cyc(4,indt(dir),dir) = idn           ! destination node #
              ipathst_cyc(5,indt(dir),dir) = id            ! destination node i
              ipathst_cyc(6,indt(dir),dir) = jd            ! destination node j
              ipathst_cyc(7,indt(dir),dir) = avgflg_t(dir) ! average flag
              if(is_u(dir) .gt. 0) then
               indu(dir) = indu(dir) + 1
               ipathsu_cyc(1,indu(dir),dir) = isn           ! source node #
               ipathsu_cyc(2,indu(dir),dir) = is_u(dir)     ! source node i
               ipathsu_cyc(3,indu(dir),dir) = js_u(dir)     ! source node j
               ipathsu_cyc(4,indu(dir),dir) = idn           ! destination node #
               ipathsu_cyc(5,indu(dir),dir) = id            ! destination node i
               ipathsu_cyc(6,indu(dir),dir) = jd            ! destination node j
               ipathsu_cyc(7,indu(dir),dir) = avgflg_u(dir) ! average flag
              endif
              if(is_v(dir) .gt. 0) then
               indv(dir) = indv(dir) + 1
               ipathsv_cyc(1,indv(dir),dir) = isn           ! source node #
               ipathsv_cyc(2,indv(dir),dir) = is_v(dir)     ! source node i
               ipathsv_cyc(3,indv(dir),dir) = js_v(dir)     ! source node j
               ipathsv_cyc(4,indv(dir),dir) = idn           ! destination node #
               ipathsv_cyc(5,indv(dir),dir) = id            ! destination node i
               ipathsv_cyc(6,indv(dir),dir) = jd            ! destination node j
               ipathsv_cyc(7,indv(dir),dir) = avgflg_v(dir) ! average flag
              endif

              exit ! since we found the source node
             endif
            enddo ! source nodes

           enddo  ! dir = 1, 2

           exit ! since we found the destination node
          endif
         enddo ! destination nodes

      endif ! id, jd along boundaries
   enddo ! id
enddo ! jd

return
END SUBROUTINE par_node_paths
