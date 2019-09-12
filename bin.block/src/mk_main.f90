
program mk_main

use mk_blocks_mod

implicit none

character(len=256) :: source
integer, external :: iargc

real, allocatable :: datain(:,:)

if (iargc() < 1) then
  print*, 'usage: '
  print*, 'executable  input_file  output_directory  prefix'
  print*, '(ie. mk-blkfiles-6.1.0  131231_  output  S)'  
  print*, '(see SST example in directory example_modis_sst)'
  print*, '(file 131231_ is ASCII lat,lon,SST(C) MODIS data)'
  print*, '(files S131231_90S000E.h5 S131231_90S180W.h5 are the output)'
  print*, '(file SHEADER is the header file)'
  print*, ' '
  print*, 'where:'
  print*, '  input_file  - name of input data file'
  print*, '  header_path - directory where header file will be written'
  print*, '  prefix      - filename prefix for header and data files'
  stop ' usage'
endif

! Get the command line arguments
call getarg(1,source)  ! input file name
call getarg(2,hpath)   ! output file path (and location of header file)
call getarg(3,fpref)   ! output file name prefix for header and blocks

! Set the attributes of the input data. These could be set in an input namelist
!    but let's hardwire these in here for now:

!      Input parameters:
!
!         nx      - Number of x (longitude) gripoints in input array  
!         ny      - Number of y (latitude) gripoints in input array 
!         iwlon   - West longitude of data (-180 to 180 degrees)
!         islat   - South latitude of data (-90 to 90 degrees)
!         tres    - data resolution in degrees
!         iblsize - Size in degrees of output files, they will have the
!                   same number of points in x and y.
!         field_name - name given to field. Only relevant for hdf5 files
!
!     Note that iwlon, islat, and iblsize are all integers.
  
nx = 4321
ny = 2161
iwlon = -180
islat = -90
tres = .08333333
iblsize = 180
offlat=0.
offlon=0.                                             
field_name='MODIS_9km_sst'

! Allocate array to hold input data and scratch

allocate (datain(nx+1,ny+1))

call read_input(source,datain,nx,ny)

call mk_blocks(source,datain,nx,ny)

end

!#############################################################################
subroutine read_input(source,datain,n1,n2)

use mk_blocks_mod

implicit none

integer :: n1,n2
real :: datain(n1,n2),lat(n2),lon(n1)
character(len=256) :: source

integer :: i,j

! Read the input data. 
! This will have to be modified depending on actual files....

      open(30,status='old',file=source,form='formatted')
4     FORMAT(f8.2,f8.2,E11.3)
      DO I=1,4321
        DO J=1,2161
          read(30,4) LAT(J),LON(I),datain(I,J)             
          datain(I,J) = datain(I,J) + 273.16
        enddo
      enddo
      close(30)

! Any modifications to input? This example takes sst, overlaps a row in EW-NS,
!    then converts to K

!datain(4321,1:2161) = datain(1,1:2161)
!datain(1:4321,2161) = datain(1:4321,2160)
!datain = datain + 273.16

return
end
