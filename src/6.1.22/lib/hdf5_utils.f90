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

Module hdf5_utils

implicit none

Contains

Subroutine shdf5_open (locfn,access,idelete)
        
implicit none

character(len=*) :: locfn     ! file name
character(len=*) :: access    ! File access ('R','W','RW')
integer, optional :: idelete  ! If W, delete/overwrite file if exists? 1=yes, 0=no
                              ! Only needed when access='W'

integer :: hdferr ! Error flag
integer :: iaccess ! int access flag
character(len=2) :: caccess ! File access ('R ','W ','RW')

logical :: exists ! File existence

caccess = access

! Check for existence of RAMS file.

inquire(file=trim(locfn),exist=exists)

! Create a new file or open an existing RAMS file.
if (access(1:1) == 'R') then
   if (.not.exists) then
      print*,'shdf5_open:'
      print*,'   Attempt to open a file for reading that does not exist.'
      print*,'   Filename: ',trim(locfn)
      stop 'shdf5_open: no file'
   else
      if (caccess == 'R ') iaccess = 1
      if (caccess == 'RW') iaccess = 2
      CALL fh5f_open (trim(locfn)//char(0), iaccess, hdferr)
      
      if (hdferr < 0) then
         print*,'shdf5_open:'
         print*,'   Error opening hdf5 file - error -',hdferr
         print*,'   Filename: ',trim(locfn)
         stop 'shdf5_open: open error'      
      endif
   endif
elseif (access(1:1) == 'W') then
   if (.not.exists) then
      iaccess=2
      CALL fh5f_create (trim(locfn)//char(0), iaccess, hdferr)
   else
      if(.not.present(idelete) ) then
         print*,'shdf5_open: idelete not specified when access=W'
         stop 'shdf5_open: no idelete'
      endif
      
      if(idelete == 0) then
         print*,'In shdf5_open:'
         print*,'   Attempt to open an existing file for writing, '
         print*,'      but overwrite is disabled. idelete=',idelete
         print*,'   Filename: ',trim(locfn)
         stop 'shdf5_open'
      else
         CALL usystem ('rm -f '//trim(locfn)//char(0))
         iaccess=1
         CALL fh5f_create (trim(locfn)//char(0), iaccess, hdferr)
      endif
   endif
   if(hdferr < 0) then
      print*,'HDF5 file create failed:',hdferr
      print*,'file name:',trim(locfn),' ',trim(access), idelete
      stop 'shdf5_open: bad create'
   endif
endif

return
END SUBROUTINE shdf5_open

!##############################################################################
Subroutine shdf5_info (dsetname,ndims,dims)

implicit none

character(len=*) :: dsetname ! Dataset name
integer :: dims(*)
integer :: ndims ! Dataset rank (in file)
integer, dimension(4) :: dimsh
integer :: hdferr ! Error flag

! Open the dataset.
CALL fh5d_open (trim(dsetname)//char(0), hdferr)
if (hdferr < 0) then
   print*,'In h5info: call fh5d_open:',trim(dsetname),'hdf5 error =',hdferr
   stop   'shdf5_info'
endif

! Get dataset's dimensions

CALL fh5s_get_ndims (ndims)

! Need to reverse dims - see comments in the shdf5_orec routine
CALL fh5s_get_dims (dimsh)
CALL shdf5_reverse_dims (ndims, dimsh, dims)

print*,'name,ndims,dims: ',dsetname
print*,ndims,dims(1:ndims)

CALL fh5d_close (hdferr)

return
END SUBROUTINE shdf5_info

!##############################################################################
Subroutine shdf5_orec (ndims,dims,dsetname,ivara,rvara,cvara,dvara,lvara  &
                                          ,ivars,rvars,cvars,dvars,lvars)
    
implicit none

character(len=*) :: dsetname ! Variable label
integer :: ndims             ! Number of dimensions or rank
integer, dimension(*) :: dims ! Dataset dimensions.

! Array and scalar arguments for different types. Only specify one in each call.
integer, optional :: ivara(*),ivars
real,    optional :: rvara(*),rvars
character(len=*), optional :: cvara(*),cvars
double precision,    optional :: dvara(*),dvars
logical,   optional :: lvara(*),lvars

integer:: h5_type   ! Local type designator

integer, dimension(4) :: dimsh ! Dataset dimensions.

character(len=2) :: ctype    ! Variable type: int, real, char
integer :: hdferr ! Error flag

! Find which data type is input
    if(present(ivars)) then ; ctype='is'
elseif(present(rvars)) then ; ctype='rs'
elseif(present(cvars)) then ; ctype='cs'
elseif(present(dvars)) then ; ctype='ds'
elseif(present(lvars)) then ; ctype='ls'
elseif(present(ivara)) then ; ctype='ia'
elseif(present(rvara)) then ; ctype='ra'
elseif(present(cvara)) then ; ctype='ca'
elseif(present(dvara)) then ; ctype='da'
elseif(present(lvara)) then ; ctype='la'
else
   print*,'Incorrect or missing data field argument in shdf5_orec'
   stop 'shdf5_orec: bad data field'
endif

! Check dimensions and set compression chunk size

if (ndims <=0 .or. minval(dims(1:ndims)) <=0) then
   print*,'Dimension error in shdf5_orec:',ndims,dims(1:ndims)
   stop 'shdf5_orec: bad dims'
endif

! Set up hyperslab selection specs for memory and file spaces
!
! Since we are working in FORTRAN and switching to C for the HDF5 I/O,
! we need to reverse the specification of dimension sizes. This is
! necessary to keep the data in the HDF5 file consistent. This need
! comes about due to the column-major array storage order in FORTRAN
! versus the row-major storage in C. Eg., 2D array in FORTRAN that
! has 3 rows and 2 columns:
!    1 4
!    2 5
!    3 6
!
! It's typical to associate i with the rows and j with the columns,
! ie the array is indexed using (i,j).
!
! FORTRAN stores this in linear memory as:
!    1 2 3 4 5 6
!
! If we tell C that this linear storage represents an array that is
! 3 rows and 2 columns, C will see this as:
!    1 2
!    3 4
!    5 6
! Note that the association with (i,j) is now broken. Eg, A(i,j) where
! i = 1,2,3 and j is held at 1 gives 1 3 5, which in FORTRAN would have
! given 1 2 3.
!
! However, if we reverse the specs of the dimension sizes (tell C that
! the array is 2 rows by 3 columns), then C will see this as:
!    1 2 3
!    4 5 6
! Along with reversing the dimension sizes we also reverse the association
! with i and j: now the rows are associated with j and the columns are
! associated with i. With this scheme, when C is asked to give A(j,i)
! where i - 1,2,3 and j is held at 1, this results in 1 2 3 which matches
! what you would get in FORTRAN.
!
! Note: This could be also be accomplished by transposing the data in
! FORTRAN right before handing it to C, but it is way more efficient
! to reverse the dimension specs.

CALL shdf5_reverse_dims (ndims,dims,dimsh)
     
! Prepare memory and options for the write
CALL fh5_prepare_write (ndims,dimsh,hdferr)
if (hdferr /= 0) then
   print*,'shdf5_orec: cannot prepare requested field:',trim(dsetname)
   return
endif

if (ctype(1:1) == 'i') h5_type=1
if (ctype(1:1) == 'r') h5_type=2
if (ctype(1:1) == 'c') h5_type=3
if (ctype(1:1) == 'd') h5_type=4
if (ctype(1:1) == 'l') h5_type=5

! Write the dataset.
if (ctype == 'is') then
   CALL fh5_write (h5_type,ivars,trim(dsetname)//char(0),hdferr)
elseif (ctype == 'rs') then
   CALL fh5_write (h5_type,rvars,trim(dsetname)//char(0),hdferr)
elseif (ctype == 'cs') then
   CALL fh5_write (h5_type,cvars,trim(dsetname)//char(0),hdferr)
elseif (ctype == 'ds') then
   CALL fh5_write (h5_type,dvars,trim(dsetname)//char(0),hdferr)
elseif (ctype == 'ls') then
   CALL fh5_write (h5_type,lvars,trim(dsetname)//char(0),hdferr)
elseif (ctype == 'ia') then
   CALL fh5_write (h5_type,ivara,trim(dsetname)//char(0),hdferr)
elseif (ctype == 'ra') then
   CALL fh5_write (h5_type,rvara,trim(dsetname)//char(0),hdferr)
elseif (ctype == 'ca') then
   CALL fh5_write (h5_type,cvara,trim(dsetname)//char(0),hdferr)
elseif (ctype == 'da') then
   CALL fh5_write (h5_type,dvara,trim(dsetname)//char(0),hdferr)
elseif (ctype == 'la') then
   CALL fh5_write (h5_type,lvara,trim(dsetname)//char(0),hdferr)
endif

if (hdferr /= 0) then
   print*,'In shdf5_orec: hdf5 write error =',hdferr
   stop 'shdf5_orec: hdf5 write error'
endif


! Close the dataset, the dataspace for the dataset, and the dataspace properties.

CALL fh5_close_write (hdferr)

return
END SUBROUTINE shdf5_orec

!##############################################################################
Subroutine shdf5_irec (dsetname,ivara,rvara,cvara,dvara,lvara  &
                               ,ivars,rvars,cvars,dvars,lvars)
        
implicit none

character(len=*) :: dsetname ! Dataset name

! Array and scalar arguments for different types. Only specify one in each call.
integer, optional :: ivara(*),ivars
real,    optional :: rvara(*),rvars
character(len=*), optional :: cvara(*),cvars
double precision,    optional :: dvara(*),dvars
logical,   optional :: lvara(*),lvars

integer:: h5_type   ! Local type designator

integer :: hdferr ! Error flag

character(len=2) :: ctype

! Find which data type will be read
    if(present(ivars)) then ; ctype='is'
elseif(present(rvars)) then ; ctype='rs'
elseif(present(cvars)) then ; ctype='cs'
elseif(present(dvars)) then ; ctype='ds'
elseif(present(lvars)) then ; ctype='ls'
elseif(present(ivara)) then ; ctype='ia'
elseif(present(rvara)) then ; ctype='ra'
elseif(present(cvara)) then ; ctype='ca'
elseif(present(dvara)) then ; ctype='da'
elseif(present(lvara)) then ; ctype='la'
else
   print*,'Incorrect or missing data field argument in shdf5_irec'
   stop 'shdf5_irec: bad data field'
endif
    
! Prepare file and memory space for the read
CALL fh5_prepare_read (trim(dsetname)//char(0),hdferr)
if (hdferr < 0) then
   print*,'shdf5_irec: cannot prepare requested field:',trim(dsetname)
   return
endif

! Read data from hyperslab in the file into the hyperslab in memory.
if (ctype(1:1) == 'i') h5_type=1
if (ctype(1:1) == 'r') h5_type=2
if (ctype(1:1) == 'c') h5_type=3
if (ctype(1:1) == 'd') h5_type=4
if (ctype(1:1) == 'l') h5_type=5

if (ctype == 'is') then
   CALL fh5d_read (h5_type,ivars,hdferr)
elseif (ctype == 'rs') then
   CALL fh5d_read (h5_type,rvars,hdferr)
elseif (ctype == 'cs') then
   CALL fh5d_read (h5_type,cvars,hdferr)
elseif (ctype == 'ds') then
   CALL fh5d_read (h5_type,dvars,hdferr)
elseif (ctype == 'ls') then
   CALL fh5d_read (h5_type,lvars,hdferr)
elseif (ctype == 'ia') then
   CALL fh5d_read (h5_type,ivara,hdferr)
elseif (ctype == 'ra') then
   CALL fh5d_read (h5_type,rvara,hdferr)
elseif (ctype == 'ca') then
   CALL fh5d_read (h5_type,cvara,hdferr)
elseif (ctype == 'da') then
   CALL fh5d_read (h5_type,dvara,hdferr)
elseif (ctype == 'la') then
   CALL fh5d_read (h5_type,lvara,hdferr)
endif

if (hdferr /= 0) then
   print*,'shdf5_irec: call fh5d_read: hdf5 error =',hdferr
   stop
endif

! Close the dataset, the dataspace for the dataset, and the memory space.

CALL fh5_close_read (hdferr)

return
END SUBROUTINE shdf5_irec

!##############################################################################
Subroutine shdf5_close ()
        
implicit none

integer :: hdferr  ! Error flags

! Close RAMS hdf file.

CALL fh5f_close (hdferr)

return
END SUBROUTINE shdf5_close

!#############################################################################
Subroutine shdf5_reverse_dims (ndims,dims,rev_dims)

implicit none

  integer :: ndims
  integer, dimension(ndims) :: dims, rev_dims

  integer :: i, irev

  do i = 1,ndims
    irev = (ndims-i) + 1
    rev_dims(irev) = dims(i)
  enddo

return
END SUBROUTINE shdf5_reverse_dims

END MODULE hdf5_utils
