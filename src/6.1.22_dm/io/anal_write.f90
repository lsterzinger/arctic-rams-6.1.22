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

Subroutine anal_write (vtype)

use node_mod
use an_header
use var_tables
use mem_scratch
use mem_basic
use mem_turb
use mem_grid
use io_params
use anal_extra
use hdf5_utils
use micro_prm

implicit none

! This routine writes the chosen variables on the analysis file.

character(len=*) :: vtype
character(len=strl1) :: anamel,anamelh,command
character(len=2) :: cgrid
character(len=32) :: subaname,varn
character(len=1) :: vnam
logical :: exans,skip
integer, save :: ioaunt=10,ncall_head=0,nvtota=0,nvtotl=0,nvtotm=0  &
                ,nvtot
integer :: ngr,nv,nvcnt,lenl,npointer,iwrite,npts,idtype,nvl,litecheck
real :: timeold
logical :: first_call=.true.
character(len=strl1), save :: anameold
real, save :: time_save
integer :: ndims,idims(4)
real, pointer :: v_pointer
integer :: h5_fid, iphdf5
type (hdf5_select_type) :: mem_select,file_select
integer, dimension(HDF5_MAX_DIMS) :: file_chunks

type (head_table), allocatable,save :: aw_table(:)

if (ioutput == 0) return

if (nmachs .gt. 1) then
  iphdf5 = 1
else
  iphdf5 = 0
endif

if (ncall_head == 0) then
   ! Get "extra" variable count
   CALL anal_extra_init ()

   !  Find total number of fields to be written
   do ngr=1,ngrids
      do nv = 1,num_var(ngr)
         if ( vtab_r(nv,ngr)%ianal == 1) nvtota=nvtota+1
         if ( vtab_r(nv,ngr)%ilite == 1) nvtotl=nvtotl+1
         if ( vtab_r(nv,ngr)%imean == 1) nvtotm=nvtotm+1
      enddo
   enddo
   nvtot=max( nvtota + ngrids*num_extra_anal, nvtotl &
             ,nvtotm + ngrids*num_extra_anal)
   allocate (aw_table(nvtot))
   ncall_head=1
endif

if(print_msg) print '(a,4I6)',' NumVars(anal,extra,lite,mean):' &
   ,nvtota,num_extra_anal,nvtotl,nvtotm

timeold=time
if(vtype == 'MEAN'.or.vtype == 'BOTH') time=min(time,time-avgtim/2.)


! Construct header file name

if(vtype == 'INST') vnam='A'
if(vtype == 'LITE') vnam='L'
if(vtype == 'MEAN') vnam='M'
if(vtype == 'BOTH') vnam='B'
CALL makefnam (anamelh,afilepref,time,iyear1,imonth1,idate1,  &
     itime1*100,vnam,'head','txt')

! Loop through each nest

nvcnt=0

do ngr=1,ngrids
   
   ! If this routine is called, then 1 or more grids will get written.
   !     See if this isn't one of them...
   
   if(vtype == 'INST') then
      ! If this routine is called, then 1 or more grids will get written.
      !     See if this isn't one of them...
      if(mod(time,frqstate(ngr)) > dtlongn(1).and.  &
                         time  <  timmax - .01*dtlongn(1) .and.  &
                         iflag == 0) cycle 
   else
      ! If it isn't the instantaneous files, then we will be outputting all
      !  grids. We would need grid dependent frqmean, etc, and maybe averaging
      !  times also to vary the grids for other write types.
   endif
                        
      
   write(cgrid,'(a1,i1)') 'g',ngr
   CALL makefnam (anamel,afilepref,time,iyear1,imonth1,idate1,  &
           itime1*100,vnam,cgrid,'h5')

   lenl = len_trim(anamel)

   inquire(file=anamel,exist=exans)
   if(exans.and.iclobber == 0) then
      if(print_msg) then
        print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        print*,'!!!   trying to open file name :'
        print*,'!!!       ',anamel
        print*,'!!!   but it already exists. run is ended.'
        print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      endif
      stop 'anal_write'
   endif
   
   if(print_msg) then
     print*,'++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
     print*,'+ anal_write: open file:',trim(anamel)
     print*,'++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
   endif
   CALL shdf5_open (anamel,'W',iphdf5,h5_fid,iclobber)
   npointer=0

!  Loop through the main variable table and write those variables
!     with the correct flag set

   do nv = 1,num_var(ngr)
      !print*,'rio:',vtab_r(nv,ngr)%name,vtab_r(nv,ngr)%ilite

      iwrite=0
      if(vtype == 'INST' .and. vtab_r(nv,ngr)%ianal == 1) then
         iwrite=1
         v_pointer => vtab_r(nv,ngr)%var_p
      elseif(vtype == 'LITE' .and. vtab_r(nv,ngr)%ilite == 1) then
         iwrite=1
         v_pointer => vtab_r(nv,ngr)%var_p
      elseif(vtype == 'MEAN' .and. vtab_r(nv,ngr)%imean == 1 .and. &
                                   vtab_r(nv,ngr)%ianal == 1) then
         iwrite=1
         v_pointer => vtab_r(nv,ngr)%var_m
      elseif(vtype == 'BOTH' .and. vtab_r(nv,ngr)%ilite == 1 .and. &
                                   vtab_r(nv,ngr)%imean == 1) then
         iwrite=1
         v_pointer => vtab_r(nv,ngr)%var_m
      endif

      if(iwrite == 1) then
      
         varn= vtab_r(nv,ngr)%name
         if(iprntstmt>=2)print*,'anal_write: ',vtype,nv,varn,max_extra_anal

         ! All cases below in the following if..elseif..else statement
         ! will use the following assignment. The different cases just
         ! add on extra dims (beyond the first 2) as necessary.
         ndims=2 ; idims(1)=mmxp(ngr) ; idims(2)=mmyp(ngr)
         
         if(vtab_r(nv,ngr)%idim_type == 3) then
            !  Rearrange 3-d variables to (i,j,k)
            CALL rearrange (mmzp(ngr),mmxp(ngr),mmyp(ngr)  &
                          ,v_pointer,scratch%scr2(1))
            ndims=3 ; idims(3)=mmzp(ngr)
      
         elseif(vtab_r(nv,ngr)%idim_type == 4) then
            !  Rearrange 4-d leaf%soil variables to (i,j,k,ip)
            CALL rearrange_p (mmxp(ngr),mmyp(ngr),nzg,npatch  &
                            ,v_pointer,scratch%scr2(1))
            ndims=4 ; idims(3)=nzg ; idims(4)=npatch
      
         elseif(vtab_r(nv,ngr)%idim_type == 5) then
            !  Rearrange 4-d leaf%sfcwater variables to (i,j,k,ip)
            CALL rearrange_p (mmxp(ngr),mmyp(ngr),nzs,npatch  &
                            ,v_pointer,scratch%scr2(1))
            ndims=4 ; idims(3)=nzs ; idims(4)=npatch
         
         elseif(vtab_r(nv,ngr)%idim_type == 6) then
            ndims=3 ; idims(3)=npatch
            npts=product(idims(1:ndims))
            CALL atob (npts,v_pointer,scratch%scr2)

         elseif(vtab_r(nv,ngr)%idim_type == 7) then
            !  Rearrange 4-d micro bin variables to (i,j,k,kr)
            call rearrange_p(mmxp(ngr),mmyp(ngr),mmzp(ngr),nkr  &
                            ,v_pointer,scratch%scr2(1))
            ndims=4 ; idims(3)=mmzp(ngr); idims(4)=nkr
         else
            ! 2D var, copy into scratch, dims have already been set before this
            ! if..elseif..else statement.
            npts=product(idims(1:ndims))
            CALL atob (npts,v_pointer,scratch%scr2)
         endif

         nvcnt=nvcnt+1
         aw_table(nvcnt)%string=varn
         aw_table(nvcnt)%npointer=npointer
         aw_table(nvcnt)%idim_type=vtab_r(nv,ngr)%idim_type
         aw_table(nvcnt)%ngrid=ngr
         !Here npts contains number of points for each variable on sub-domains
         aw_table(nvcnt)%nvalues=vtab_r(nv,ngr)%npts
         !In the header file, output the total grid points for whole domain
         if(aw_table(nvcnt)%idim_type==2) &
            aw_table(nvcnt)%nvalues=nnxp(ngr)*nnyp(ngr)
         if(aw_table(nvcnt)%idim_type==3) &
            aw_table(nvcnt)%nvalues=nnxp(ngr)*nnyp(ngr)*nnzp(ngr)
         if(aw_table(nvcnt)%idim_type==4) &
            aw_table(nvcnt)%nvalues=nnxp(ngr)*nnyp(ngr)*nzg*npatch
         if(aw_table(nvcnt)%idim_type==5) &
            aw_table(nvcnt)%nvalues=nnxp(ngr)*nnyp(ngr)*nzs*npatch
         if(aw_table(nvcnt)%idim_type==6) &
            aw_table(nvcnt)%nvalues=nnxp(ngr)*nnyp(ngr)*npatch
         if(aw_table(nvcnt)%idim_type==7) &
            aw_table(nvcnt)%nvalues=nnxp(ngr)*nnyp(ngr)*nnzp(ngr)*nkr

         CALL shdf5_set_hs_select (vtab_r(nv,ngr)%idim_type,'W',ngr &
                ,mem_select,file_select,file_chunks)
         CALL shdf5_orec (h5_fid,iphdf5,varn &
                ,mem_select,file_select,file_chunks,rvara=scratch%scr2)
         !print*,nvcnt,varn,' : done'
      endif
   enddo
   
   ! Now loop through the "extra" variables, putting each in scratch%scr2...
   if(vtype == 'INST' .or. vtype == 'LITE') then  
    do nv = 1,max_extra_anal
      varn=an_extra(nv)%name
      if(iprntstmt>=2)print*,'anal_write_extra: ',vtype,nv,varn,max_extra_anal

      !See if the "extra" variables are in the list of LITE variables.
      !Only output "extra" variables to LITE files if they are in the list.
      if(vtype == 'LITE') then
       litecheck=0
       do nvl=1,nlite_vars
        if (varn == lite_vars(nvl)) litecheck=1
       enddo
       if(litecheck==0)cycle
      endif

      idtype = an_extra(nv)%idim_type
      CALL anal_extra_comp (nv,scratch%scr1,ngr,skip)
      
      if (.not.skip) then
         idims(1)=mmxp(ngr) ; idims(2)=mmyp(ngr)
      
         if(idtype == 2) then
            ndims=2 
            npts=product(idims(1:ndims))
            CALL atob (npts,scratch%scr1,scratch%scr2)

         elseif(idtype == 3) then
            ndims=3 ; idims(3)=mmzp(ngr)
            npts=product(idims(1:ndims))
            CALL rearrange (mmzp(ngr),mmxp(ngr),mmyp(ngr)  &
                           ,scratch%scr1,scratch%scr2)

         elseif(idtype == 4) then
            !  Rearrange 4-d leaf%soil variables to (i,j,k,ip)
            CALL rearrange_p (mmxp(ngr),mmyp(ngr),nzg,npatch  &
                           ,scratch%scr1,scratch%scr2)
            ndims=4 ; idims(3)=nzg ; idims(4)=npatch

         elseif(idtype == 5) then
            !  Rearrange 4-d leaf%sfcwater variables to (i,j,k,ip)
            CALL rearrange_p (mmxp(ngr),mmyp(ngr),nzs,npatch  &
                           ,scratch%scr1,scratch%scr2)
            ndims=4 ; idims(3)=nzs ; idims(4)=npatch

         elseif(idtype == 6) then
            ndims=3 ; idims(3)=npatch
            npts=product(idims(1:ndims))
            CALL atob (npts,scratch%scr1,scratch%scr2)

         elseif(idtype == 7) then
            ndims=4 ; idims(3)=mmzp(ngr); idims(4)=nkr
            !  Rearrange 4-d micro bin variables to (i,j,k,kr)
            CALL rearrange_p(mmxp(ngr),mmyp(ngr),mmzp(ngr),nkr  &
                           ,scratch%scr1,scratch%scr2)

         endif
      
         nvcnt=nvcnt+1
         aw_table(nvcnt)%string=an_extra(nv)%name
         aw_table(nvcnt)%npointer=1 ! not needed for hdf
         aw_table(nvcnt)%idim_type=idtype
         aw_table(nvcnt)%ngrid=ngr
         !Here npts contains number of points for each variable on sub-domains
         aw_table(nvcnt)%nvalues=npts
         !In the header file, output the total grid points for whole domain
         if(aw_table(nvcnt)%idim_type==2) &
            aw_table(nvcnt)%nvalues=nnxp(ngr)*nnyp(ngr)
         if(aw_table(nvcnt)%idim_type==3) &
            aw_table(nvcnt)%nvalues=nnxp(ngr)*nnyp(ngr)*nnzp(ngr)
         if(aw_table(nvcnt)%idim_type==4) &
            aw_table(nvcnt)%nvalues=nnxp(ngr)*nnyp(ngr)*nzg*npatch
         if(aw_table(nvcnt)%idim_type==5) &
            aw_table(nvcnt)%nvalues=nnxp(ngr)*nnyp(ngr)*nzs*npatch
         if(aw_table(nvcnt)%idim_type==6) &
            aw_table(nvcnt)%nvalues=nnxp(ngr)*nnyp(ngr)*npatch
         if(aw_table(nvcnt)%idim_type==7) &
            aw_table(nvcnt)%nvalues=nnxp(ngr)*nnyp(ngr)*nnzp(ngr)*nkr

         CALL shdf5_set_hs_select (idtype,'W',ngr &
                 ,mem_select,file_select,file_chunks)
         CALL shdf5_orec (h5_fid,iphdf5,varn &
                 ,mem_select,file_select,file_chunks,rvara=scratch%scr2)
         !print*,nvcnt,varn,' : done'
      endif
    enddo
   endif

   CALL shdf5_close (h5_fid)

enddo

! Write the header information out to the file. Have only
! the mainnum node do this so that you don't have many processors
! trying to write this file at the same time.

if ((my_rams_num .eq. mainnum) .or. (nmachs .eq. 1)) then
  CALL rams_f_open (ioaunt,anamelh,'FORMATTED','REPLACE','WRITE',iclobber)

  write(ioaunt,110) nvcnt
  do nv=1,nvcnt
     write(ioaunt,120) aw_table(nv)%string   &
                      ,aw_table(nv)%npointer  &
                      ,aw_table(nv)%idim_type  &
                      ,aw_table(nv)%ngrid  &
                      ,aw_table(nv)%nvalues
  enddo

110 format(i6)
120 format(a16,1x,i12,i3,i3,1x,i9)

  CALL commio ('WRITE',ioaunt)
  close(ioaunt)
endif

if(vtype == 'LITE')then
   subaname='  Analysis lite write          '
elseif(vtype == 'MEAN')then
   subaname='  Averaged analysis write      '
elseif(vtype == 'BOTH')then
   subaname='  Averaged analysis lite write '
else
   subaname='  Analysis write               '
endif

if(print_msg) then
  print 12,subaname,time,anamelh
12 format(/,1X,79('*'),/,  &
       A35,'  Time = ',F9.0,/,'      Header file name - ',A60  &
    ,/,1X,79('*'))
endif

! Reset the time back to the original
if(vtype == 'MEAN'.or.vtype == 'BOTH')  time=timeold


! See if previous state file should be deleted. This feature is for writing
!   frequent files for checkpointing, but only saving a subset of them.
!   DO NOT remove the old state file if IFLAG is set.

if(frqst_keep > 0.) then
   if(.not.first_call .and. iflag == 0 .and.  &
                            mod(time_save,frqst_keep) /= 0.) then
      command = 'rm -f '//anameold(1:len_trim(anameold)-5)//'*'
      if(print_msg) then
        print*,'*************************************************'
        print*,'*  Removing state file with command: '
        print*,'*    ',trim(command)
        print*,'*************************************************'
      endif
      CALL usystem (command)
   endif
   anameold = anamel
   time_save = time
   first_call = .false.
endif

return
END SUBROUTINE anal_write
