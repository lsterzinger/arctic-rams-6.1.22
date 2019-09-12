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

Subroutine read_rams (action,ivtime,cgrid,ctime,iztrans,cvar,rh5_file,cvar_ok)

use an_header
use rcommons
use mem_grid

implicit none

character(len=*) :: action,cgrid,ctime
character(len=24) :: cvar,cdname,cdunits
character(len=strl1) :: flnm,cfile
integer :: ivtime,iztrans,rh5_file,cvar_ok,itime,idate,ng,ivtype &
          ,nf,lenhf,nv,nngd,idims(3)
integer, external :: lastchar
real, allocatable, save :: arra(:),arrb(:),arrbb(:)  &
                          ,topt(:),sfclat(:),sfclon(:)
real :: alat1(maxgrds),alon1(maxgrds),glats1,glatn2,glonw1,glone2,fcstsec

data alat1 /maxgrds*0./, alon1 /maxgrds*0./

!*************************************************************************
!Determine header file name to open
!*************************************************************************
read(ctime,*) nf
CALL rams_get_cdata (0,nf,flnm)
lenhf=len_trim(flnm)
write(cfile,*) flnm(1:lenhf-9)

!*************************************************************************
!Open and read header file and extra data
!*************************************************************************
open(10,file=flnm)
read(10,*) nvbtab
if(allocated(anal_table)) deallocate (anal_table)
allocate (anal_table(nvbtab))
do nv=1,nvbtab
   read(10,*)  anal_table(nv)%string   &
              ,anal_table(nv)%npointer  &
              ,anal_table(nv)%idim_type  &
              ,anal_table(nv)%ngrid  &
              ,anal_table(nv)%nvalues
enddo
CALL commio ('READ',10)
close(10)
   
CALL rams_get_fdata (1,nf,fcstsec)
CALL rams_get_fdata (2,nf,startutc)

!*************************************************************************
!Set up memory allocation based on max grid dimensions
!*************************************************************************
maxmem=0
do ng=1,ngrids
   maxmem=max(maxmem,(nnxp(ng)+1)*(nnyp(ng)+1)*  &
            max(nnzp(ng),(nzg+nzs+4)*npatch))
enddo

!print*,'compute memory',ngrids,maxmem

if(allocated(arra))    deallocate(arra);    allocate(arra(maxmem))
if(allocated(arrb))    deallocate(arrb);    allocate(arrb(maxmem))
if(allocated(arrbb))   deallocate(arrbb);   allocate(arrbb(maxmem))
if(allocated(topt))    deallocate(topt);    allocate(topt(maxmem))
if(allocated(sfclat))  deallocate(sfclat);  allocate(sfclat(maxmem))
if(allocated(sfclon))  deallocate(sfclon);  allocate(sfclon(maxmem))

!*************************************************************************
!Establish grid dimensions
!*************************************************************************
read(cgrid,*) nngd
idims(1)=nnxp(nngd)
idims(2)=nnyp(nngd)
idims(3)=nnzp(nngd)

!*************************************************************************
!Set up HORIZONTAL beginning, ending, and grid point increments to output
!*************************************************************************
if(ixbeg.le.0) nib=1-ixbeg
if(ixbeg.gt.0) nib=ixbeg
nib=max(0,min(nib,idims(1)))
if(ixend.le.0) nie=idims(1)+ixend
if(ixend.gt.0) nie=ixend
nie=max(nib,min(nie,idims(1)))

if(iybeg.le.0) njb=1-iybeg
if(iybeg.gt.0) njb=iybeg
njb=max(0,min(njb,idims(2)))
if(iyend.le.0) nje=idims(2)+iyend
if(iyend.gt.0) nje=iyend
nje=max(njb,min(nje,idims(2)))

niinc=ixstep
njinc=iystep

nii=nie-nib+1
njj=nje-njb+1

! Exit cleanly if grid too small to view
if((nii <= 1.and.idims(1) > 1).or.(nii < 1.and.idims(1) == 1)) then
   print*,'Warning: Selected grid too small in x direction'
   stop
endif
if((njj <= 1.and.idims(2) > 1).or.(njj < 1.and.idims(2) == 1)) then
   print*,'Warning: Selected grid too small in y direction'
   stop
endif

!*************************************************************************
!Finding lat/lon corners for grid
!*************************************************************************
if(action(1:4)=='TEXT'.or.action(1:4)=='HDF5') then

   CALL rams_fill_fld (idims(1),idims(2),idims(3)  &
                      ,arra,arrb,arrbb  &
                      ,flnm(1:lenhf-9),nngd,'lat',iztrans  &
                      ,cdname,cdunits,ivtype,topt)

   alat1(nngd)=arra(1) ! 1st point
   CALL thelatlon (idims(1),idims(2),idims(3),arra,alat1(nngd),nib,njb)
   CALL thelatlon (idims(1),idims(2),idims(3),arra,glats1,nib,njb)
   CALL thelatlon (idims(1),idims(2),idims(3),arra,glatn2,nie,nje)

   ! copy the surface latitude values from arra for the HDF5 format routine
   CALL copy_sfc_values (idims(1),idims(2),idims(3),arra,sfclat)

   !DEBUG
   ! look at latitude values
   !call textarray (idims(1),idims(2),idims(3),arra,sfclat,0)

   CALL rams_fill_fld (idims(1),idims(2),idims(3)  &
                      ,arra,arrb,arrbb  &
                      ,flnm(1:lenhf-9),nngd,'lon',iztrans  &
                      ,cdname,cdunits,ivtype,topt)

   alon1(nngd)=arra(1) ! 1st point
   CALL thelatlon (idims(1),idims(2),idims(3),arra,alon1(nngd),nib,njb)
   CALL thelatlon (idims(1),idims(2),idims(3),arra,glonw1,nib,njb)
   CALL thelatlon (idims(1),idims(2),idims(3),arra,glone2,nie,nje)

   ! copy the surface longitude values from arra for the HDF5 format routine
   CALL copy_sfc_values (idims(1),idims(2),idims(3),arra,sfclon)

   !DEBUG
   ! look at longitude values
   !call textarray (idims(1),idims(2),idims(3),arra,sfclon,1)
endif

!*************************************************************************
!Read primary field
!*************************************************************************
cvar_ok=0
if(cvar(1:4).ne.'none') then
   CALL rams_fill_fld (idims(1),idims(2),idims(3)  &
                      ,arra,arrb,arrbb  &
                      ,flnm(1:lenhf-9),nngd,cvar  &
                      ,iztrans,cdname,cdunits,ivtype,topt)
   if(ivtype==0.and.cvar(1:4).ne.'none') then
      print*,'Primary variable not found: ',cvar(1:len_trim(cvar))
      cdname=cvar(1:len_trim(cvar))//' not found'
      cdunits=';'
      return
   endif
   cvar_ok=1
endif

!*************************************************************************
!Set up beginning, ending, and grid point increments to output
!*************************************************************************
if(izbeg.le.0) nnb=1-izbeg
if(izbeg.gt.0) nnb=izbeg

if(ivtype.le.3.and.iztrans.eq.3) nnb=max(0,min(nnb,nplevs))
if(ivtype.le.3.and.iztrans.ne.3) nnb=max(0,min(nnb,idims(3)))
if(ivtype==6) nnb=max(0,min(nnb,npatch))
if(ivtype==5) nnb=max(0,min(nnb,nzg))
if(ivtype==4) nnb=max(0,min(nnb,nzs))
if(ivtype==2) nnb=1
!For now, override ivtype=2,4,5,6 so we output all physical levels
if(ivtype.ne.3)nnb=1

if(izend.le.0) then
   if(ivtype.le.3.and.iztrans.eq.3) nne=nplevs+izend
   if(ivtype.le.3.and.iztrans.ne.3) nne=idims(3)+izend
   if(ivtype==6) nne=npatch+izend
   if(ivtype==5) nne=nzg+izend
   if(ivtype==4) nne=nzs+izend
   if(ivtype==2) nne=1
endif
if(izend.gt.0) nne=izend

if(ivtype.le.3.and.iztrans.eq.3) nne=max(nnb,min(nne,nplevs))
if(ivtype.le.3.and.iztrans.ne.3) nne=max(nnb,min(nne,idims(3)))
if(ivtype==6) nne=max(nnb,min(nne,npatch))
if(ivtype==5) nne=max(nnb,min(nne,nzg))
if(ivtype==4) nne=max(nnb,min(nne,nzs))
if(ivtype==2) nne=1
!For now, override ivtype=2,4,5,6 so we output all physical levels
if(ivtype==6) nne=npatch
if(ivtype==5) nne=nzg
if(ivtype==4) nne=nzs
if(ivtype==2) nne=1

nninc=izstep

if(ivtype.eq.2 .or.ivtype.eq.4 .or. ivtype.eq.5 .or. ivtype.eq.6) nninc=1

!print*,'nib,nie,nii,niinc',nib,nie,nii,niinc
!print*,'njb,nje,njj,njinc',njb,nje,njj,njinc
!print*,'nnb,nne,nninc',nnb,nne,nninc
!print*,'idims(1),idims(2),idims(3)',idims(1),idims(2),idims(3)

!*************************************************************************
!Choose what to do with fields now that they are read and optionally
!interpolated
!*************************************************************************
if(action(1:4).eq.'TEXT') then

   ! if action is 'TEXT', call routine to output field

   CALL rams_text (arra,iztrans,ivtype,nngd,idims(1),idims(2)  &
                  ,idims(3),fcstsec,glats1,glonw1,glatn2,glone2 &
                  ,polelat,polelon,cvar)

elseif(action(1:4).eq.'HDF5') then

   ! if action is 'HDF5', call routine to output field

   CALL rams_get_idata (3,nf,itime)
   CALL rams_get_idata (2,nf,idate)
   CALL rams_hdf5 (rh5_file,arra,iztrans,ivtype,nngd,idims(1),idims(2)  &
        ,idims(3),ztn(1,nngd),nzg,idate,itime,ivtime  &
        ,fcstsec,cvar,cdname,cdunits,sfclat,sfclon)

endif

return
END SUBROUTINE read_rams

!##############################################################################
Subroutine copy_sfc_values (n1,n2,n3,a_3d,a_2d)

! This routine will copy the lowest level of a_3d (3D array) to
! a_2d (2D array). a_3d is organized as (i,j,k) and a_2d is organized
! as (i,j) where i, j and k step though spatial x, y, and z directions
! respectively.

implicit none

integer :: n1,n2,n3,i,j
real, dimension(n1,n2,n3) :: a_3d
real, dimension(n1,n2) :: a_2d

do i = 1, n1
  do j = 1, n2
    a_2d(i,j) = a_3d(i,j,1)
  enddo
enddo

return
END SUBROUTINE copy_sfc_values

!##############################################################################
Subroutine textarray (n1,n2,n3,array,array2d,TextCols)

implicit none

integer :: n1,n2,n3,TextCols,i
real, dimension(n1,n2,n3) :: array
real, dimension(n1,n2) :: array2d

print*,'DEBUG: TextArray:'
if (TextCols == 1) then
  print*,'DEBUG:   Columns: ', n1
  do i = 1,n1
    print*,'DEBUG:    ',i,array(i,1,1),array(i,n2/2,1),array(i,n2,1)
    print*,'DEBUG:    ',i,array2d(i,1),array2d(i,n2/2),array2d(i,n2)
    print*,''
  enddo
else
  print*,'DEBUG:   Rows: ', n2
  do i = 1,n2
    print*,'DEBUG:   ',i,array(1,i,1),array(n1/2,i,1),array(n1,i,1)
    print*,'DEBUG:   ',i,array2d(1,i),array2d(n1/2,i),array2d(n1,i)
    print*,''
  enddo
endif

return
END SUBROUTINE textarray
