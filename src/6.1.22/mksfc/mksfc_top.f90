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

Subroutine top_read (ifm)

use mem_grid
use io_params
use hdf5_utils

implicit none

integer :: ifm
character(len=strl1) :: flnm
character(len=2) :: cgrid
logical :: there
integer :: ndims,idims(4)

! read the "top" file

write(cgrid,'(a1,i1)') 'g',ifm
flnm=trim(topfiles)//'-S-'//cgrid//'.h5'

inquire(file=flnm,exist=there)

if(.not.there) then
   print*,'------------------------------------------------'
   print*,'TOP_read: file for grid ',ifm,' not there.'
   print*,'TOP_read: file:',trim(flnm)
   print*,'------------------------------------------------'
   stop 'top_read: no file'
endif

CALL shdf5_open (flnm,'R')
ndims=2 ; idims(1)=nnxp(ifm) ; idims(2)=nnyp(ifm)
CALL shdf5_irec ('TOPT',rvara=grid_g(ifm)%topt)
CALL shdf5_irec ('TOPZO',rvara=grid_g(ifm)%topzo)
CALL shdf5_close ()

return
END SUBROUTINE top_read

!##############################################################################
Subroutine top_check (ifm,ierr)

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

integer :: lc,nsfx,nsfy,nsitoptflg,nsitopsflg,nsiz0flg
real ::  sfdx,sfplat,sfplon,sflat,sflon,stoptenh,stoptwvl  &
   ,sz0max,sz0fact,glatr,glonr

character(len=strl1) :: flnm
character(len=2) :: cgrid
logical there
integer :: ndims,idims(4)

lc=len_trim(topfiles)
write(cgrid,'(a1,i1)') 'g',ifm
flnm=trim(topfiles)//'-S-'//cgrid//'.h5'

print*,'------------------------------------------------'
print*,'---> Check grid:',ifm,' top file... '
print*,'--->   Filename:',trim(flnm)

inquire(file=flnm,exist=there)

if(.not.there) then
   ierr = 1
   print*,'TOPfile for grid ',ifm,' not there.'
   print*,'------------------------------------------------'
   return
endif

CALL xy_ll (glatr,glonr,polelat,polelon,xtn(1,ifm),ytn(1,ifm))

CALL shdf5_open (flnm,'R')
ndims=1 ; idims(1)=1
CALL shdf5_irec ('nx',ivars=nsfx)
CALL shdf5_irec ('ny',ivars=nsfy)
CALL shdf5_irec ('dx',rvars=sfdx)
CALL shdf5_irec ('polelat',rvars=sfplat)
CALL shdf5_irec ('polelon',rvars=sfplon)
CALL shdf5_irec ('sw_lat',rvars=sflat)
CALL shdf5_irec ('sw_lon',rvars=sflon)
CALL shdf5_irec ('itoptflg',ivars=nsitoptflg)
CALL shdf5_irec ('itopsflg',ivars=nsitopsflg)
CALL shdf5_irec ('toptenh',rvars=stoptenh)
CALL shdf5_irec ('toptwvl',rvars=stoptwvl)
CALL shdf5_irec ('iz0flg',ivars=nsiz0flg)
CALL shdf5_irec ('z0max',rvars=sz0max)
CALL shdf5_irec ('z0fact',rvars=sz0fact)
CALL shdf5_close ()


if (nsfx                       .ne. nnxp(ifm)     .or.  &
    nsfy                       .ne. nnyp(ifm)     .or.  &
    abs(sfdx-deltaxn(ifm))     .gt. .001          .or.  &
    abs(sfplat-polelat)        .gt. .001          .or.  &
    abs(sfplon-polelon)        .gt. .001          .or.  &
    abs(sflat-glatr)           .gt. .001          .or.  &
    abs(sflon-glonr)           .gt. .001          .or.  &
    nsitoptflg                 .ne. itoptflg(ifm) .or.  &
    nsitopsflg                 .ne. itopsflg(ifm) .or.  &
    abs(stoptenh-toptenh(ifm)) .gt. .001          .or.  &
    abs(stoptwvl-toptwvl(ifm)) .gt. .001          .or.  &
    nsiz0flg                   .ne. iz0flg(ifm)   .or.  &
    abs(sz0max-z0max(ifm))     .gt. .001          .or.  &
    abs(sz0fact-z0fact)        .gt. .00001) then

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
   print*,'itoptflg:',itoptflg(ifm),nsitoptflg
   print*,'itopsflg:',itopsflg(ifm),nsitopsflg
   print*,'toptenh:',toptenh(ifm),stoptenh
   print*,'toptwvl:',toptwvl(ifm),stoptwvl
   print*,'iz0flg:',iz0flg(ifm),nsiz0flg
   print*,'z0max:',z0max(ifm),sz0max
   print*,'z0fact:',z0fact,sz0fact
   print*,'-------------------'

else

   ierr = 0
   print*,'---> Grid:',ifm,' topography file data okay. '
   print*,'------------------------------------------------'

endif

return
END SUBROUTINE top_check

!##############################################################################
Subroutine toptinit (n2,n3,topt,topzo)

implicit none

integer :: n2,n3,i,j
real, dimension(n2,n3) :: topt,topzo

! Fill the TOPT array with a default value of 0.  This default is used only
! when a standard RAMS topography dataset is not used and when no overrides
! to topography heights are defined in routine toptinit_user in the
! file ruser.f.

do j = 1,n3
   do i = 1,n2
      topt(i,j) = 0.
      topzo(i,j) = .0001
   enddo
enddo

return
END SUBROUTINE toptinit

!##############################################################################
Subroutine top_write (ifm)

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

!     write surface characteristics, one file for each grid

write(cgrid,'(a1,i1)') 'g',ifm

flnm=trim(topfiles)//'-S-'//cgrid//'.h5'

CALL xy_ll (glatr,glonr,polelat,polelon,xtn(1,ifm),ytn(1,ifm))

CALL shdf5_open (flnm,'W',iclobber)
ndims=1 ; idims(1)=1
CALL shdf5_orec (ndims,idims,'nx',ivars=nnxp(ifm))
CALL shdf5_orec (ndims,idims,'ny',ivars=nnyp(ifm))
CALL shdf5_orec (ndims,idims,'dx',rvars=deltaxn(ifm))
CALL shdf5_orec (ndims,idims,'polelat',rvars=polelat)
CALL shdf5_orec (ndims,idims,'polelon',rvars=polelon)
CALL shdf5_orec (ndims,idims,'sw_lat',rvars=glatr)
CALL shdf5_orec (ndims,idims,'sw_lon',rvars=glonr)
CALL shdf5_orec (ndims,idims,'itoptflg',ivars=itoptflg(ifm))
CALL shdf5_orec (ndims,idims,'itopsflg',ivars=itopsflg(ifm))
CALL shdf5_orec (ndims,idims,'toptenh',rvars=toptenh(ifm))
CALL shdf5_orec (ndims,idims,'toptwvl',rvars=toptwvl(ifm))
CALL shdf5_orec (ndims,idims,'iz0flg',ivars=iz0flg(ifm))
CALL shdf5_orec (ndims,idims,'z0max',rvars=z0max(ifm))
CALL shdf5_orec (ndims,idims,'z0fact',rvars=z0fact)
ndims=2 ; idims(1)=nnxp(ifm) ; idims(2)=nnyp(ifm)
CALL shdf5_orec (ndims,idims,'TOPT',rvara=sfcfile_p(ifm)%topt)
CALL shdf5_orec (ndims,idims,'TOPZO',rvara=sfcfile_p(ifm)%topzo)
CALL shdf5_close ()

return
END SUBROUTINE top_write
