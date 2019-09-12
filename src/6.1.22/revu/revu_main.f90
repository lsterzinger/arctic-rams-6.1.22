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

Program revu

use rcommons

implicit none

integer :: nfl

integer :: i,numarg,iargc,bad=0
character(len=strl1) :: name_name,arg,cargs(0:maxargs)

! argument defaults
name_name='REVUIN'

! read arguments
do i=1,maxargs
 cargs(i)=''
enddo

numarg=iargc()

do i=0,numarg
   CALL ugetarg (i,arg)
   cargs(i)=trim(arg)//char(0) !Null terminate the string to pass to C
enddo

if(numarg > 0)then
 if(numarg == 2)then
   if(cargs(1)(1:1)=='-' .and. cargs(1)(1:2)=='-f') then
    if(len_trim(cargs(2)).gt.80) then
       print*,'Max filename length = 80 characters'
       bad=1
    endif
    name_name=cargs(2)(1:len_trim(cargs(2)))
   else
    bad=1
   endif
 else
   bad=1
 endif
endif

if (bad > 0) then
   print*,'REVU usage: ''exec name'' '
   print*,'  [-f ''Namelist file''] '
   stop 'bad command line arguments'
endif

print*,''
print*,'Starting REVU with: ',trim(name_name)
print*,''

! give some values to namelist parameters that may be missing
do i=1,maxrevu
  revuvar(i)=''
enddo

! read in the namelist
open(1,status='old',FILE=name_name)
CALL namein_revu (1,'$CONTROL')
close(1)

! set up the background configuration from input strings
if(len_trim(TVAR).eq.0.or.len_trim(ZVAR).eq.0.or.  &
   len_trim(XVAR).eq.0.or.len_trim(YVAR).eq.0) then
   print*,'TVAR: ',TVAR,' ZVAR: ',ZVAR
   print*,'XVAR: ',XVAR,' YVAR: ',YVAR
   stop 'plotspc: frame not defined'
endif

CALL backset ()

! fill the rams arrays
CALL rams_anal_init (nfl,anpref)

! This supports new version where all grids may not be available at all times.
!   Do grid inventory by looking at file names
CALL ra_grid_inv ()

! do the stuff
print*
write(*,'(2a)') 'mode: ',anatype

if(anatype.eq.'TEXT'.or.anatype.eq.'HDF5') then

   CALL plotspc (nfl,anatype)

else

   print*,'No action for: ',anatype
   print*,'  options= HDF5 or TEXT'

endif

print*
print*,'REVU finished normally'
print*

END PROGRAM revu

!##############################################################################
Subroutine ra_grid_inv ()

use rcommons

implicit none

integer :: nhfiles,n,nc,ngfiles,ng,ngr
character(len=strl1) :: hfiles(maxfiles),fpref,gfiles(10)

fpref=trim(anpref)//'*-head.txt'

CALL rams_filelist (hfiles,fpref,nhfiles)

do n=1,nhfiles
   nc=len_trim(hfiles(n))-9
   fpref=hfiles(n)(1:nc)//'-g*'
   CALL rams_filelist (gfiles,fpref,ngfiles)
   ! print*,'ngfiles:',ngfiles
   do ng=1,ngfiles
      ! print*,'gfiles:',ng,trim(gfiles(ng))
      ! Find last "." in file name. Grid number is character before it.
      do nc=len_trim(gfiles(ng)),1,-1
         if (gfiles(ng)(nc:nc) == '.') exit
      enddo
      nc=nc-1
      read(gfiles(ng)(nc:nc),*) ngr
   enddo
enddo

return
END SUBROUTINE ra_grid_inv

!##############################################################################
Subroutine nvfilla (GROUP,VA,ISUB1,ISUB2,IN,FV,CH,NVA)

use rcommons

implicit none

! Namelist reading routine, called by NAMEIN.

character(len=*) :: GROUP,VA,CH
integer, parameter :: NVCONTRL=10
integer :: ICONTROL(NVCONTRL)
character(len=8) :: CONTROL(NVCONTRL)
integer :: is1,is2,isub1,isub2,nva,in,inrflg
real :: fv

DATA ICONTROL/NVCONTRL*0/
DATA CONTROL/'ANATYPE','ANPREF','REVPREF','XVAR','YVAR','ZVAR','TVAR'  &
            ,'IGRID','IZTRAN','REVUVAR'/

IS1=ISUB1
IS2=ISUB2
IF(ISUB1.EQ.0)IS1=NVA

if(GROUP.eq.'$CONTROL') then
  CALL varchk (VA,GROUP,CONTROL,ICONTROL,NVCONTRL,INRFLG)
  IF(INRFLG.EQ.1) RETURN
  IF(VA.EQ.'ANPREF')   CALL varsetc (VA,ANPREF,IS1,1,CH,0,128)
  IF(VA.EQ.'REVPREF')  CALL varsetc (VA,REVPREF,IS1,1,CH,0,128)
  IF(VA.EQ.'ANATYPE')  CALL varsetc (VA,ANATYPE,IS1,1,CH,0,8)
  IF(VA.EQ.'XVAR')     CALL varsetc (VA,XVAR,IS1,1,CH,0,20)
  IF(VA.EQ.'YVAR')     CALL varsetc (VA,YVAR,IS1,1,CH,0,20)
  IF(VA.EQ.'ZVAR')     CALL varsetc (VA,ZVAR,IS1,1,CH,0,20)
  IF(VA.EQ.'TVAR')     CALL varsetc (VA,TVAR,IS1,1,CH,0,20)
  IF(VA.EQ.'IGRID')    CALL varseti (VA,IGRID,IS1,1,IN,0,9)
  IF(VA.EQ.'IZTRAN')   CALL varseti (VA,IZTRAN,IS1,1,IN,1,3)
  IF(VA.EQ.'REVUVAR')  CALL varsetc (VA,REVUVAR(IS1),IS1,MAXREVU,CH,0,64)
endif

return
END SUBROUTINE nvfilla

!##############################################################################
Subroutine namein_revu (IUNIT,GROUP)

use grid_dims

implicit none

! This routine is called by routines INITLZ, CONSTAT, and ARRSND
! to input the values for the NAMELIST specified by the character stri
! GROUP from input unit IUNIT.

integer :: iunit,nvalue,nvarn,nr,ncw,ntok,nt,isub1,isub2,int
real :: fnum
character(len=*) :: GROUP(*)
character(len=20) :: VARN
character(len=6) :: linelength
character(len=strl1) :: LINE,LINEW,VALUE(MAXVALUES),TOKENS(50)

integer, external :: letter
integer, external :: numberchk
integer, external :: letquo

REWIND IUNIT
NVALUE=0
NVARN=0

CALL findgr (IUNIT,GROUP)

DO NR=1,MAXREC
   write(linelength,'("(A",I3,")" )')  strl1
   READ(IUNIT,linelength,END=100,ERR=100) LINE
   CALL strip (LINE,LINEW,NCW)
   NCW=MAX(NCW,1)
   CALL toknze (LINEW,NCW,TOKENS,NTOK)
   IF(LINEW(1:NCW).EQ.'$END') THEN
      GOTO 100
   ENDIF

   NT=1
   20 CONTINUE
      IF(NT.GT.NTOK) GOTO 10
      IF(letter(TOKENS(NT)).EQ.1) THEN
         IF(NVARN.GT.0)  &
            CALL nvtran_revu (GROUP,VARN,ISUB1,ISUB2,VALUE,NVALUE)
         NVALUE=0
         ISUB1=0
         ISUB2=0
         VARN=TOKENS(NT)
         NVARN=NVARN+1
         NT=NT+1
         IF(TOKENS(NT).EQ.'(')THEN
            NT=NT+1
            CALL ch2int (TOKENS(NT),ISUB1)
            NT=NT+1
            IF(TOKENS(NT).EQ.',') THEN
               NT=NT+1
               CALL ch2int (TOKENS(NT),ISUB2)
               NT=NT+1
            ENDIF
         ENDIF
      ELSEIF(numberchk(TOKENS(NT)).EQ.1 .OR. letquo(TOKENS(NT)).EQ.1) THEN
         NVALUE=NVALUE+1
         VALUE(NVALUE)=TOKENS(NT)
      ENDIF
      NT=NT+1
   GOTO 20
   10 CONTINUE
ENDDO

100 CONTINUE
CALL nvtran_revu (GROUP,VARN,ISUB1,ISUB2,VALUE,NVALUE)
VARN='$END'
CALL nvfilla (GROUP,VARN,ISUB1,ISUB2,INT,FNUM,LINEW(1:NCW),NR)

return
END SUBROUTINE namein_revu

!##############################################################################
Subroutine nvtran_revu (GROUP,VARN,ISUB1,ISUB2,VALUE,NVALUE)

use grid_dims

implicit none

! This routine converts a variable value(s) (VALUE) read in as a
! character string(s) to the proper type (integer, real, or character)
! and assigns it (them) to the storage location(s) corresponding to
! the NAMELIST variable VARN(ISUB1,ISUB2).  GROUP is the NAMELIST
! group name.  NVALUE is the number of values to be assigned.

integer :: ISUB1,ISUB2,NVALUE
character(len=*) :: GROUP,VARN,VALUE(*)
character(len=strl1) :: CHVAL
integer :: nv,int,ncw
real :: fnum
integer, external :: letint
integer, external :: letquo

IF(letint(VARN).EQ.1) THEN
   DO NV=1,NVALUE
      IF(letquo(VALUE(NV)).EQ.0) THEN
         CALL ch2int (VALUE(NV),INT)
         CALL nvfilla (GROUP,VARN,ISUB1,ISUB2,INT,FNUM,CHVAL(1:NCW),NV)
      ELSE
         CALL ch2ch (VALUE(NV),CHVAL,NCW)
         CALL nvfilla (GROUP,VARN,ISUB1,ISUB2,INT,FNUM,CHVAL(1:NCW),NV)
      ENDIF
   ENDDO
ELSE
   DO NV=1,NVALUE
      IF(letquo(VALUE(NV)).EQ.0) THEN
         CALL ch2real (VALUE(NV),FNUM)
         CALL nvfilla (GROUP,VARN,ISUB1,ISUB2,INT,FNUM,CHVAL(1:NCW),NV)
      ELSE
         CALL ch2ch (VALUE(NV),CHVAL,NCW)
         CALL nvfilla (GROUP,VARN,ISUB1,ISUB2,INT,FNUM,CHVAL(1:NCW),NV)
      ENDIF
   ENDDO
ENDIF

return
END SUBROUTINE nvtran_revu
