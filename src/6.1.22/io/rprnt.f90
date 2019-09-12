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

Subroutine prtopt ()

use mem_grid
use mem_scratch
use ref_sounding
use rconstants
use mem_turb
use ref_sounding
use ke_coms

implicit none

integer :: k

IF(initial == 1 .or. (initial == 3 .and. initorig == 1))THEN
  CALL mrsl (nsndg,ps(1),ts(1),vctr5(1))
  do k=1,nsndg
     vctr1(k) = 100. * rts(k) / vctr5(k)
  enddo
  WRITE(*,41)
41      FORMAT(/,'------------------------------SOUNDING INPUT-------'  &
   ,'---------------------------',//,7X,'PS',9X,'HS',7X,'TS',6X  &
   ,'THDS',6X,'US',7X,'VS',7X,'RTS',5X,'REL HUM',/,6X  &
   ,'(Pa)',7X,'(m)',6X,'(K)',6X,'(K)',6X,'(m/s)',4X,'(m/s)'  &
   ,3X,'(kg/kg)',5X,'(%)',/)
  WRITE(*,42)(PS(K),HS(K),TS(K),THDS(K),US(K),VS(K),RTS(K)  &
             ,VCTR1(K),K=1,NSndg)
42      FORMAT(1X,F11.1,F10.1,2F9.2,2F9.2,F10.5,F9.1)
ENDIF
!
DO K=1,NNZP(1)
  VCTR1(K)=P00*(PI01DN(K,1)/CP)**CPOR
  VCTR2(K)=TH01DN(K,1)/(1.+.61*RT01DN(K,1))
ENDDO
WRITE(*,310)IREF,JREF,TOPREF,(ZTN(K,1),U01DN(K,1),V01DN(K,1)  &
  ,DN01DN(K,1),PI01DN(K,1),VCTR1(K),TH01DN(K,1),VCTR2(K)  &
  ,RT01DN(K,1),K=1,NNZP(1))
310   FORMAT(/,'--------REFERENCE STATE at I,J=(',I4,',',I4  &
      ,')   SFC ELEV (M)= ',F6.1,'-------------'  &
 ,//,4X,'Z',6X,'U01D',4X,'V01D',4X,'DN01D',4X  &
 ,'PI01D',4X,'PRESS',4X,'TH01D',4X,'THD',6X,'RT01D'  &
 ,/,3X,'(m)',5X,'(m/s)',3X,'(m/s)',2X,'(kg/m3)',2X  &
 ,'(J/kgK)',4X,'(Pa)',5X,'(K)',5X,'(K)',5X,'(kg/kg)'  &
 ,//,(1X,F7.1,2F8.2,F8.3,F10.2,F10.1,2F8.2,F10.5))
!
!STC..................................................................
!  Print out the values chosen for the empirical constants in
!  E-l and E-eps closures (they are the same for all the grids)
!  (S. Trini Castelli)
!STC..................................................................
!
WRITE(*,*) ' '
if (IDIFFK(1).eq.5) then
WRITE(*,*) 'Empirical constants and other parameters for E-l closure'
WRITE(*,105) ' ',C_MI,C_EPS,ALF_THT,ALF_TKE,ALF_EPS
if(IOPZL.eq.1) WRITE(*,*)  &
 'IOPZL=1 Constant asymptotic mixing length',AL0_CONST,' from Ying'
if(IOPZL.eq.2) WRITE(*,*)  &
 'IOPZL=2 Asymptotic mixing length from Mellor-Yamada'
if(IOPZL.eq.3) WRITE(*,*)  &
 'IOPZL=3 Asymptotic mixing length from Zilitinkevich'
endif
if (IDIFFK(1).eq.6) then
WRITE(*,*) 'Empirical constants and other parameters for E-eps closure'
WRITE(*,106) ' ',C_MI,C1_EPS,C2_EPS,ALF_THT,ALF_TKE,ALF_EPS
endif

105  FORMAT(A1,'C_MI=',f5.2,'  C_EPS=',f5.2,'  ALF_THT=',f5.2  &
      ,'  ALF_TKE=',f5.2,'  ALF_EPS=',f5.2,999(A1,/,I13,3I16))
106  FORMAT(A1,'C_MI=',f5.2,'  C1_EPS=',f5.2,'  C2_EPS=',f5.2,  &
     '  ALF_THT=',f5.2,'  ALF_TKE=',f5.2,'  ALF_EPS=', &
      f5.2,999(A1,/,I13,3I16))

return
END SUBROUTINE prtopt
