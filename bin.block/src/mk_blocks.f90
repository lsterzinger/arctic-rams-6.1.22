!-------------------------------------------------------------------------
!      Subroutine to convert input lat-lon array
!      and reformat it into several smaller HDF
!      block files for input into RAMS.
!
!-------------------------------------------------------------------------
!         
Module mk_blocks_mod

integer :: nx, ny, iwlon, islat, iblsize, ibldim
real :: tres, offlat, offlon

character(len=256) :: hpath,fpref,subdir,field_name

character(len=80) :: icent, iyear,imonth,idate,itime

!      Before the mk_blocks routine is called, the preceeding variables 
!             must be filled
!
!        nx       - Number of x (longitude) gripoints in input array  
!        ny       - Number of y (latitude) gripoints in input array 
!        iwlon    - West longitude of data (-180 to 180 degrees)
!        islat    - South latitude of data (-90 to 90 degrees)
!        tres     - data resolution in degrees
!        iblsize  - Size in degrees of output files, they will have the
!                   same number of points in x and y.
!        offlat   - latitude offset - if data is not defined at integer latitudes,
!                                     but staggered by a bit
!        offlon   - longitude offset - same as offlat, but for longitude
!        hpath    - directory path where HEADER file will be written
!        fpref    - file name prefix for HEADER and data files
!        subdir   - directory path relative to hpath directory where data files 
!                   will be written
!        iyear    - year data is valid (yyyy) (0000 = climatology)
!        imonth   - month data is valid (mm)
!        idate    - date data is valid (dd)
!        itime    - time data is valid (hhmmss)
!        
!     Note that iwlon, islat, and iblsize are integers.

Contains

subroutine mk_blocks(source,datain,n1,n2)

use hdf5_utils

implicit none

integer :: n1,n2
real, dimension(n1,n2) :: datain

integer :: nsqx, nsqy, nsx, nsy
integer :: i, j, lat, lon, npt, jj, ii, nc
real, dimension(:),   allocatable :: scr
character(len=256) :: title1,title2,title3,fname,source

integer,external :: iargc

! hdf parameters
integer, parameter :: rank = 1,dfacc_create = 4,dfnt_float32 = 5
integer, parameter :: comp_code_deflate = 4 ,deflate_level = 6
character, parameter :: sds_name*9 = ' RAMS_data'
integer :: sd_id, sds_id, sds_index, status
integer, external :: sfstart, sfcreate, sfwdata, sfendacc, sfend, sfscompress
integer :: ndims, idims(4)

ibldim=int(float(iblsize)/tres+.001)+1

allocate (scr((ibldim+1)**2))

nsqx=nx/(ibldim-1)
nsqy=ny/(ibldim-1)

! Make sure hpath and subdir have slashes on the end...
nc=len_trim(hpath)  ; if(hpath(nc:nc) /= '/') hpath(nc+1:nc+1)='/'
nc=len_trim(subdir) ; if(subdir(nc:nc) /= '/') subdir(nc+1:nc+1)='/'

! Write the header file for this. If you are do more than one time, the header file
!   will have to be manually combined for now.
title3=trim(hpath)//trim(fpref)//'HEADER'
open(29,status='replace',file=title3,form='formatted')
print*, 'making file-',trim(title3)
write(29,'(4i5,2f10.6,1x,a)')iblsize,ibldim,islat,iwlon,offlat &
     ,offlon,trim(field_name)
write(29,'(i4)')1
icent='20'
iyear=source(1:2)
imonth=source(3:4)
idate=source(5:6)
itime='00'
write(29,'(1X,a7,1X,a2,a2,3x,a2,3x,a2,3x,a2)')source,icent,iyear,imonth,idate,itime
close(29)

title3=trim(fpref)//source
do nsx=1,nsqx
   do nsy=1,nsqy
      i=(nsx-1)*(ibldim-1)+1
      j=(nsy-1)*(ibldim-1)+1
      lat=islat+(nsy-1)*iblsize
      lon=iwlon+(nsx-1)*iblsize
      if(lat >= 0)then
          write(title1,'(i2.2,a1)')lat,'N'
      else
          write(title1,'(i2.2,a1)')abs(lat),'S'
      endif
      if(lon >= 0)then
          write(title2,'(i3.3,a1)')lon,'E'
      else
          write(title2,'(i3.3,a1)')abs(lon),'W'
      endif
      fname=trim(hpath)//trim(title3)//title1(1:3)//title2(1:4)//'.h5'
      !print*, nsx,nsy,i,j,lat,lon
      print*, 'making file-',trim(fname)
      npt=1
      do jj=j,(ibldim-1)+j
         do ii=i,(ibldim-1)+i
            scr(npt)=datain(ii,jj)
            npt=npt+1
         enddo
      enddo
      ! Write the data file block
      ! Write the HDF5 file
      call shdf5_open(fname,'W',1)
      ndims = 2 ; idims(1)=ibldim ; idims(2)=ibldim
      call shdf5_orec(ndims,idims,field_name,rvara=scr)
      call shdf5_close()
   enddo
enddo

end subroutine mk_blocks

end Module mk_blocks_mod
