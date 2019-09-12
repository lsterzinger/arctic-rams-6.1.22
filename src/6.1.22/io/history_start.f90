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

Subroutine history_start ()

! This routine initializes the model from the history file

use io_params
use mem_grid
use mem_basic

implicit none


integer :: maxarr,ngr
integer :: ifm,icm

! Find maximum size of any array on history file. Allocate scratch array of
! this size.

maxarr=0
do ngr=1,ngridsh
   maxarr=max(maxarr,nnxp(ngr)*nnyp(ngr)*nnzp(ngr)  &
         ,nnxp(ngr)*nnyp(ngr)*nzg*npatch &
         ,nnxp(ngr)*nnyp(ngr)*nzs*npatch)
enddo

! read stuff here

CALL hist_read (maxarr,hfilin)

if(print_msg) print*,'back from read'

do ifm = 1,ngrids
   icm = nxtnest(ifm)
   if (icm  ==  0) then
      CALL newgrid (ifm)
      CALL refs3d (nzp,nxp,nyp  &
      ,basic_g(ifm)%pi0  (1,1,1),basic_g(ifm)%dn0  (1,1,1)  &
      ,basic_g(ifm)%dn0u (1,1,1),basic_g(ifm)%dn0v (1,1,1)  &
      ,basic_g(ifm)%th0  (1,1,1),grid_g(ifm)%topt  (1,1)    &
      ,grid_g(ifm)%rtgt  (1,1)  )
   endif
enddo

return
END SUBROUTINE history_start

!##############################################################################
Subroutine read_hheader (hnamein)

use an_header
use mem_grid
use ref_sounding
use node_mod

implicit none

  character (len=*) :: hnamein

  integer :: ngrids1,nnxp1(maxgrds),nnyp1(maxgrds),nnzp1(maxgrds) &
            ,nzg1,nzs1,npatch1

  integer :: ie,ngr,nv
  character(len=2) :: cng
  integer, external :: cio_i
  integer, external :: cio_f
  integer :: iunhd=11

! Read the input history header file and collect:
!   Environmental sounding
!   Simulation time
!   Number of grids and grid specs
!     In case RAMSIN is specifying more grids
!   List of variables contained in the history file
!
! Note that read_hheader() will set the global variable
! ngridsh which is used in the following loop and in hist_read().
! read_hheader() also sets the global variables nvbtab
! and anal_table of which read_hist uses.
!
! Note that read_hheader() will allocate anal_table, and
! hist_read() will deallocate anal_table.

    CALL rams_f_open (iunhd,trim(hnamein),'FORMATTED','OLD','READ',0)

    !Get history grid structure info so we can allocate space
    ie=cio_i(iunhd,1,'ngrids',ngrids1,1)
    ngridsh=ngrids1
    ie=cio_i(iunhd,1,'nnxp',nnxp1,ngrids1)
    ie=cio_i(iunhd,1,'nnyp',nnyp1,ngrids1)
    ie=cio_i(iunhd,1,'nnzp',nnzp1,ngrids1)
    ie=cio_i(iunhd,1,'npatch',npatch1,1)
    ie=cio_i(iunhd,1,'nzg',nzg1,1)
    ie=cio_i(iunhd,1,'nzs',nzs1,1)
    ie=cio_f(iunhd,1,'time',time,1)
    
    !Flag to determine if we are past first timesetp
    ie=cio_i(iunhd,1,'ngbegun',ngbegun,ngrids)
    
    !Get the sounding state (needed just for standard output)
    ie=cio_i(iunhd,1,'iref',iref,1)
    ie=cio_i(iunhd,1,'jref',jref,1)
    ie=cio_f(iunhd,1,'topref',topref,1)
    ie=cio_i(iunhd,1,'nsndg',nsndg,1)
    ie=cio_f(iunhd,1,'us',us,nsndg)
    ie=cio_f(iunhd,1,'vs',vs,nsndg)
    ie=cio_f(iunhd,1,'ts',ts,nsndg)
    ie=cio_f(iunhd,1,'thds',thds,nsndg)
    ie=cio_f(iunhd,1,'ps',ps,nsndg)
    ie=cio_f(iunhd,1,'hs',hs,nsndg)

    !Get original simulation type
    ie=cio_i(iunhd,1,'initorig',initorig,1)
    
    ! Get the 1-d reference state
    do ngr=1,ngridsh
       write(cng,1) ngr
1           format(i2.2)
       ie=cio_f(iunhd,1,'u01dn'//cng,u01dn(1,ngr),mmzp(ngr))
       ie=cio_f(iunhd,1,'v01dn'//cng,v01dn(1,ngr),mmzp(ngr))
       ie=cio_f(iunhd,1,'pi01dn'//cng,pi01dn(1,ngr),mmzp(ngr))
       ie=cio_f(iunhd,1,'th01dn'//cng,th01dn(1,ngr),mmzp(ngr))
       ie=cio_f(iunhd,1,'dn01dn'//cng,dn01dn(1,ngr),mmzp(ngr))
       ie=cio_f(iunhd,1,'rt01dn'//cng,rt01dn(1,ngr),mmzp(ngr))
    enddo

    !  Read variable header info

    rewind(iunhd)

    read(iunhd,*) nvbtab
    allocate (anal_table(nvbtab))
    do nv=1,nvbtab
       read(iunhd,*)  anal_table(nv)%string   &
                     ,anal_table(nv)%npointer  &
                     ,anal_table(nv)%idim_type  &
                     ,anal_table(nv)%ngrid  &
                     ,anal_table(nv)%nvalues
    enddo

    close(iunhd)
 
return
END SUBROUTINE read_hheader

!##############################################################################
Subroutine hist_read (maxarr,hnamein)

use an_header
use var_tables
use mem_grid
use hdf5_utils

implicit none

integer :: maxarr,checkhist
character(len=*) :: hnamein

integer :: ngr,npts,nc,nv,nvh,ndims,idims(4)
character(len=1) :: cgrid
character(len=strl1) :: hname
real, allocatable :: scr(:)
logical :: exists

allocate (scr(maxarr))

!Check to see that all history grids are present at this time
checkhist=1
do ngr=1,ngridsh
  write(cgrid,'(i1)') ngr
  nc=len_trim(hnamein)
  hname=hnamein(1:nc-9)//'-g'//cgrid//'.h5'
  inquire(file=hname,exist=exists)
  if(.not. exists)then
   checkhist=0
  endif
enddo
if(checkhist==0)then
 if(print_msg)then
 print*,'Not all original grids are present at this time: ',hnamein(1:nc-9)
 print*,'Choose a history restart time in which all original'
 print*,'simulation grids are available.'
 print*,''
 endif
 stop
endif

do ngr=1,ngridsh

   ! Open file
   write(cgrid,'(i1)') ngr
   nc=len_trim(hnamein)
   hname=hnamein(1:nc-9)//'-g'//cgrid//'.h5'

   CALL shdf5_open (hname,'R')

   ! Loop through all variables
   varloop: do nvh=1,nvbtab
      if(ngr /= anal_table(nvh)%ngrid) cycle varloop
      
      ! See if variable should be read and stored
      do nv = 1,num_var(ngr)
         if(anal_table(nvh)%string == vtab_r(nv,ngr)%name) then
            ! there is a match on this grid. see if hist flag is set...
            if(iprntstmt>=1 .and. print_msg) &
               print*,'found: ', trim(anal_table(nvh)%string)
            if (vtab_r(nv,ngr)%ianal /= 1) cycle varloop
            if(iprntstmt>=1 .and. print_msg) &
               print*,'read : ', trim(anal_table(nvh)%string)
            
            ! We want it...read, maybe rearrange, and store it

            ! call and output variable array names and size in hist file
            if(iprntstmt>=1 .and. print_msg) then
              CALL shdf5_info (trim(anal_table(nvh)%string),ndims,idims)
              print*,'name,ndims,idims: ', trim(anal_table(nvh)%string)
              print*,ndims,idims(1:ndims)
            endif

            CALL shdf5_irec (trim(anal_table(nvh)%string),rvara=scr)
            
            npts = vtab_r(nv,ngr)%npts
            
            select case(vtab_r(nv,ngr)%idim_type)
               case(2,6) ; CALL atob (npts,scr,vtab_r(nv,ngr)%var_p)
               case(3)
                  CALL unarrange (nnzp(ngr),nnxp(ngr),nnyp(ngr) &
                                 ,scr,vtab_r(nv,ngr)%var_p)
               case(4)
                  CALL unarrange_p (nnxp(ngr),nnyp(ngr),nzg,npatch &
                                 ,scr,vtab_r(nv,ngr)%var_p)
               case(5)
                  CALL unarrange_p (nnxp(ngr),nnyp(ngr),nzs,npatch &
                                 ,scr,vtab_r(nv,ngr)%var_p)
            end select
            cycle varloop
         endif
      enddo
   enddo varloop
   
   CALL shdf5_close ()

enddo

deallocate(scr,anal_table)

return
END SUBROUTINE hist_read
