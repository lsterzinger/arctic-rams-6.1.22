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

Subroutine sfc_read (ifm)

use mem_grid
use mem_leaf
use io_params
use hdf5_utils

implicit none

integer :: ifm,ip,i,j,k
logical :: there

character(len=strl1) :: flnm
character(len=2) :: cgrid
integer :: ndims,idims(4)

real, dimension(:), allocatable :: r_scratch

! read the "sfc" file

write(cgrid,'(a1,i1)') 'g',ifm
flnm=trim(sfcfiles)//'-S-'//cgrid//'.h5'

inquire(file=flnm,exist=there)

if(.not.there) then
   print*,'------------------------------------------------'
   print*,'SFC_read: file for grid ',ifm,' not there.'
   print*,'SFC_read: file:',trim(flnm)
   print*,'------------------------------------------------'
   stop 'sfc_read: no file'
endif

CALL shdf5_open (flnm,'R')
ndims=3; idims(1)=nnxp(ifm); idims(2)=nnyp(ifm); idims(3)=npatch
CALL shdf5_irec ('PATCH_AREA',rvara=leaf_g(ifm)%patch_area)
CALL shdf5_irec ('LEAF_CLASS',rvara=leaf_g(ifm)%leaf_class)
ndims=4; idims(1)=nnxp(ifm); idims(2)=nnyp(ifm); idims(3)=nzg; idims(4)=npatch
allocate(r_scratch(nnxysp(ifm)))
CALL shdf5_irec ('SOIL_TEXT',rvara=r_scratch)
CALL unarrange_p (nnxp(ifm),nnyp(ifm),nzg,npatch,r_scratch,leaf_g(ifm)%soil_text)
deallocate(r_scratch)
CALL shdf5_close ()

!Saleeby(2014): Set grid boundaries as in leaf3 for surface characteristics
!Do this here so that boundary assignments are made before the first 
!analysis write and first timestep.
do ip = 1,npatch
  do j = 1,nnyp(ifm)
     leaf_g(ifm)%leaf_class(1,j,ip)=leaf_g(ifm)%leaf_class(2,j,ip)
     leaf_g(ifm)%patch_area(1,j,ip)=leaf_g(ifm)%patch_area(2,j,ip)
     do k = 1,nzg
      leaf_g(ifm)%soil_text(k,1,j,ip)=leaf_g(ifm)%soil_text(k,2,j,ip)
     enddo
     leaf_g(ifm)%leaf_class(nnxp(ifm),j,ip)=leaf_g(ifm)%leaf_class(nnxp(ifm)-1,j,ip)
     leaf_g(ifm)%patch_area(nnxp(ifm),j,ip)=leaf_g(ifm)%patch_area(nnxp(ifm)-1,j,ip)
     do k = 1,nzg
      leaf_g(ifm)%soil_text(k,nnxp(ifm),j,ip)=leaf_g(ifm)%soil_text(k,nnxp(ifm)-1,j,ip)
     enddo
  enddo
  if (jdim == 1) then
   do i = 1,nnxp(ifm)
      leaf_g(ifm)%leaf_class(i,1,ip)=leaf_g(ifm)%leaf_class(i,2,ip)
      leaf_g(ifm)%patch_area(i,1,ip)=leaf_g(ifm)%patch_area(i,2,ip)
      do k = 1,nzg
       leaf_g(ifm)%soil_text(k,i,1,ip)=leaf_g(ifm)%soil_text(k,i,2,ip)
      enddo
      leaf_g(ifm)%leaf_class(i,nnyp(ifm),ip)=leaf_g(ifm)%leaf_class(i,nnyp(ifm)-1,ip)
      leaf_g(ifm)%patch_area(i,nnyp(ifm),ip)=leaf_g(ifm)%patch_area(i,nnyp(ifm)-1,ip)
      do k = 1,nzg
       leaf_g(ifm)%soil_text(k,i,nnyp(ifm),ip)=leaf_g(ifm)%soil_text(k,i,nnyp(ifm)-1,ip)
      enddo
   enddo
  endif
enddo

return
END SUBROUTINE sfc_read

!##############################################################################
Subroutine sfc_check (ifm,ierr)

! This routine checks for the existence of a surface file for
! grid number ifm, and if it exists, also checks for agreement of
! grid configuration between the file and the current model run.
! If the file does not exist or does not match grid configuration,
! the flag ifileok is returned with a value of 0.  If the file
! exists and is ok, ifileok is returned with a value of 1.

use mem_grid
use io_params
use hdf5_utils

implicit none

integer :: ifm,ierr
integer :: lc,nsfx,nsfy,nsfzg ,nsivegtflg,nsisoilflg,nspatch
real ::  sfdx,sfplat,sfplon,sflat,sflon,glatr,glonr

character(len=strl1) :: flnm
character(len=2) :: cgrid
logical there
integer :: ndims,idims(4)

lc=len_trim(sfcfiles)
write(cgrid,'(a1,i1)') 'g',ifm
flnm=trim(sfcfiles)//'-S-'//cgrid//'.h5'

print*,'------------------------------------------------'
print*,'---> Check grid:',ifm,' sfc file... '
print*,'--->   Filename:',trim(flnm)

inquire(file=flnm,exist=there)

if(.not.there) then
   ierr = 1
   print*,'SFCfile for grid ',ifm,' not there.'
   print*,'------------------------------------------------'
   return
endif

CALL xy_ll (glatr,glonr,polelat,polelon,xtn(1,ifm),ytn(1,ifm))

CALL shdf5_open (flnm,'R')
ndims=1 ; idims(1)=1
CALL shdf5_irec ('nx',ivars=nsfx)
CALL shdf5_irec ('ny',ivars=nsfy)
CALL shdf5_irec ('nzg',ivars=nsfzg)
CALL shdf5_irec ('npatch',ivars=nspatch)
CALL shdf5_irec ('dx',rvars=sfdx)
CALL shdf5_irec ('polelat',rvars=sfplat)
CALL shdf5_irec ('polelon',rvars=sfplon)
CALL shdf5_irec ('sw_lat',rvars=sflat)
CALL shdf5_irec ('sw_lon',rvars=sflon)
CALL shdf5_irec ('ivegtflg',ivars=nsivegtflg)
CALL shdf5_irec ('isoilflg',ivars=nsisoilflg)
CALL shdf5_close ()


if (nsfx                       .ne. nnxp(ifm)     .or.  &
    nsfy                       .ne. nnyp(ifm)     .or.  &
    nsfzg                      .ne. nzg           .or.  &
    nspatch                    .ne. npatch        .or.  &
    abs(sfdx-deltaxn(ifm))     .gt. .001          .or.  &
    abs(sfplat-polelat)        .gt. .001          .or.  &
    abs(sfplon-polelon)        .gt. .001          .or.  &
    abs(sflat-glatr)           .gt. .001          .or.  &
    abs(sflon-glonr)           .gt. .001          .or.  &
    nsivegtflg                 .ne. ivegtflg(ifm) .or.  &
    nsisoilflg                 .ne. isoilflg(ifm)) then

   ierr = 1

   print*,'SFCfile mismatch on grid:',ifm
   print*,'Values: model, file'
   print*,'-------------------'
   print*,'nnxp:',nnxp(ifm),nsfx
   print*,'nnyp:',nnyp(ifm),nsfy
   print*,'deltax:',deltaxn(ifm),sfdx
   print*,'polelat:',polelat,sfplat
   print*,'polelon:',polelon,sfplon
   print*,'SW lat:',glatr,sflat
   print*,'SW lon:',glonr,sflon
   print*,'ivegtflg:',ivegtflg(ifm),nsivegtflg
   print*,'isoilflg:',isoilflg(ifm),nsisoilflg
   print*,'-------------------'

else

   ierr = 0
   print*,'---> Grid:',ifm,' surface file data okay. '
   print*,'------------------------------------------------'

endif

return
END SUBROUTINE sfc_check

!##############################################################################
Subroutine sfc_write (ifm)

use mem_mksfc
use mem_grid
use io_params
use hdf5_utils

implicit none

integer :: ifm
real :: glatr,glonr
character(len=strl1) :: flnm
character(len=2) :: cgrid
integer :: ndims,idims(4)

real, dimension(:), allocatable :: r_scratch

!     write surface characteristics, one file for each grid


write(cgrid,'(a1,i1)') 'g',ifm

flnm=trim(sfcfiles)//'-S-'//cgrid//'.h5'

CALL xy_ll (glatr,glonr,polelat,polelon,xtn(1,ifm),ytn(1,ifm))


CALL shdf5_open (flnm,'W',iclobber)
ndims=1 ; idims(1)=1
CALL shdf5_orec (ndims,idims,'nx',ivars=nnxp(ifm))
CALL shdf5_orec (ndims,idims,'ny',ivars=nnyp(ifm))
CALL shdf5_orec (ndims,idims,'nzg',ivars=nzg)
CALL shdf5_orec (ndims,idims,'npatch',ivars=npatch)
CALL shdf5_orec (ndims,idims,'dx',rvars=deltaxn(ifm))
CALL shdf5_orec (ndims,idims,'polelat',rvars=polelat)
CALL shdf5_orec (ndims,idims,'polelon',rvars=polelon)
CALL shdf5_orec (ndims,idims,'sw_lat',rvars=glatr)
CALL shdf5_orec (ndims,idims,'sw_lon',rvars=glonr)
CALL shdf5_orec (ndims,idims,'ivegtflg',ivars=ivegtflg(ifm))
CALL shdf5_orec (ndims,idims,'isoilflg',ivars=isoilflg(ifm))
ndims=3; idims(1)=nnxp(ifm); idims(2)=nnyp(ifm); idims(3)=npatch
CALL shdf5_orec (ndims,idims,'PATCH_AREA',rvara=sfcfile_p(ifm)%patch_area)
CALL shdf5_orec (ndims,idims,'LEAF_CLASS',rvara=sfcfile_p(ifm)%leaf_class)
ndims=4; idims(1)=nnxp(ifm); idims(2)=nnyp(ifm); idims(3)=nzg; idims(4)=npatch
allocate(r_scratch(nnxysp(ifm)))
CALL rearrange_p (nnxp(ifm),nnyp(ifm),nzg,npatch,sfcfile_p(ifm)%soil_text,r_scratch)
CALL shdf5_orec (ndims,idims,'SOIL_TEXT',rvara=r_scratch)
deallocate(r_scratch)
CALL shdf5_close ()

return
END SUBROUTINE sfc_write
