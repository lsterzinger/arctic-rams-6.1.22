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

Subroutine rams_text (a,iztrans,ivtype,nngd,n1,n2,n3,fcstsec &
                   ,swlat,swlon,nelat,nelon,polelatn,polelonn,cvar)

! Saleeby: Current get to GEMPAK text format
!-----------------------------------------------------------------
! Routine to write RAMS fields in text format as chosen through
! REVUIN. This is intended to be modified at will to produce the
! desired output. Basically, one 3 or 2-dimensional variable at a
! time will be sent here.

! Arguments:
! ----------
! cdname - name of variable
! cdunits - units of variable
! a - data
! n1,n2,n3 - actual dimensions of field
! iztrans - type of vertical transformation 1-sigma-z, 2-Cartesian, 3-pressure
! ivtype - type of variable 2-2d surface, 3-3d atmospheric
! nib,nie,njb,nje - horizontal or vertical "window" as chosen in namelist.
! nplevs - number of atmospheric coordinate levels (if pressure transformation)
! iplevs(nplevs)- atmospheric coordinate levels (if pressure transformation)
! zlev - atmospheric coordinate levels (if sigma_z or Cartesian)
! iyear1 - year of model start
! imonth1 - month of model start
! idate1 - date of model start
! itime1 - time of model start (hhmm)
! fcstsec - seconds into run
!-----------------------------------------------------------------

use an_header
use rcommons
use mem_grid

implicit none

integer :: iztrans,ivtype,n1,n2,n3,nngd,ihour2,imin2,idate2,imonth2,iyear2 &
          ,jcount,i,j,k,lastslash,iyr2,typlev,jj,kk
real :: fcstsec,swlat,swlon,nelat,nelon,polelatn,polelonn,fcstmin,fcsthrs &
          ,fcstday,fcstsec1,zlev
real :: a(n1,n2,n3)
character(len=4) :: coordn
character(len=5) :: cdnamegem
character(len=20) :: cvar
character(len=strl1) :: flnm,flnm1,out1,out2
integer, save :: ncall(maxgrds),iun,iuntag,numgrids &
                ,xbeg,xend,ybeg,yend,xinc,yinc,xpts,ypts,iyr1
real, save :: swlat1,swlon1,nelat1,nelon1

data ncall/maxgrds*0/

!print*,'===> text out =>',fcstsec
!print*,iyear1,imonth1,idate1,itime1,nngd
!print*,cdname(1:len_trim(cdname)-1),n1,n2,n3
!print*,iztrans,ivtype
!print*,nib,nie,njb,nje,nnb,nne

! If it is the first time into this routine, make and open a file.
if(ncall(nngd).eq.0) then

  if(niinc > 2 .or. njinc > 2) then
    print*,'Can only use increment of 1 or 2 right now'   
    stop
  endif

  CALL rams_get_cdata (0,1,flnm1)
  write(flnm,'(2a,2a1,i4.4,a1,i2.2,a1,i2.2,a1,i6.6,a2,i1)' )  &
      revpref(1:len_trim(revpref))  &
     ,flnm1(lastslash(flnm1)+1:len_trim(flnm1)-27),ftran,'-'  &
     ,iyear1,'-',imonth1,'-',idate1,'-',itime1*100,'-g',nngd

  write(out1,'(a,a4)'),trim(flnm),'.txt'
  iun=79
  open(iun,file=out1,status='unknown')  
  rewind iun

  write(out2,'(a,a4)')trim(flnm),'.tag'
  iuntag=78
  open(iuntag,file=out2,status='unknown')
  rewind iuntag

  swlat1 = swlat
  swlon1 = swlon
  nelat1 = nelat
  nelon1 = nelon
  xbeg = nib
  ybeg = njb
  xend = nie
  yend = nje
  xinc = niinc
  yinc = njinc
  xpts = xend/xinc-nib+1
  ypts = yend/yinc-njb+1
  
  if(nje==1) then
    ybeg = 1
    yend = 5
    yinc = 1
    ypts = 5
    nelat1 = swlat1 + 1.0
  endif

  if(iyear1 >= 2000) iyr1=iyear1-2000
  if(iyear1 <  2000) iyr1=iyear1-1900

  !write navigation file for creating gempak grid file
  write(iuntag,'(a7,a80)') &
   ,'GDEFIL=',out1
  write(iuntag,'(a7,3i2.2,i4.4,a1,i1,a9)') &
   ,'GDOUTF=',iyr1,imonth1,idate1,itime1,'g',nngd,'_rams.gem'
  write(iuntag,'(a7,3i2.2,i4.4,a1,i1,a9)') &
   ,'GDFILE=',iyr1,imonth1,idate1,itime1,'g',nngd,'_rams.gem'
  write(iuntag,'(a10,f8.2,a1,f8.2,a5)')'PROJ="str/',polelatn,';',polelonn,';0.0"'
  write(iuntag,'(a9,f8.2,a1,f8.2,a1,f8.2,a1,f8.2,a1)') &
   ,'GRDAREA="',swlat1,';',swlon1,';',nelat1,';',nelon1,'"'
  write(iuntag,'(a7,f8.2,a1,f8.2,a1,f8.2,a1,f8.2,a1)') &
   ,'GAREA="',swlat1,';',swlon1,';',nelat1,';',nelon1,'"'
  write(iuntag,'(a6,i5,a1,i5,a1)') 'KXKY="',xpts,';',ypts,'"'

  ncall(nngd)=1
  numgrids=0

endif

!************************************************************
!***** sets correct forecast time for gempak grid file *****
!************************************************************
ihour2=int(itime1/100.0)
imin2=itime1-ihour2*100.0
iyear2=iyear1
imonth2=imonth1
idate2=idate1
fcstsec1=fcstsec
fcstday=fcstsec/86400.0
if(fcstday >=1) then
 idate2=idate2+int(fcstday)
 fcstsec1=fcstsec1 - int(fcstday)*86400.0
endif
  
fcsthrs=fcstsec1/3600.0  
if(fcsthrs >=1) then
 ihour2=ihour2+int(fcsthrs)
 fcstsec1=fcstsec1 - int(fcsthrs)*3600.0
endif

fcstmin=fcstsec1/60.0
if(fcstmin >=1)then
 imin2=imin2+int(fcstmin)
 fcstsec1=fcstsec1 - int(fcstmin)*60.0
endif

if(imin2>=60) then
   imin2=imin2-60
   ihour2=ihour2+1
endif
if(ihour2>=24) then
   ihour2=ihour2-24
   idate2=idate2+1
endif  
if(imonth2 == 1 .and. idate2 > 31) then
   idate2=idate2-31
   imonth2=imonth2+1
endif
if(imonth2 == 2 .and. idate2 > 28 .and. mod(iyear2,4)>0) then
   idate2=idate2-28
   imonth2=imonth2+1
endif
if(imonth2 == 2 .and. idate2 > 29 .and. mod(iyear2,4)==0) then
   idate2=idate2-29
   imonth2=imonth2+1
endif
if(imonth2 == 3 .and. idate2 > 31) then
   idate2=idate2-31
   imonth2=imonth2+1
endif
if(imonth2 == 4 .and. idate2 > 30) then
   idate2=idate2-30
   imonth2=imonth2+1
endif
if(imonth2 == 5 .and. idate2 > 31) then
   idate2=idate2-31
   imonth2=imonth2+1
endif
if(imonth2 == 6 .and. idate2 > 30) then
   idate2=idate2-30
   imonth2=imonth2+1
endif
if(imonth2 == 7 .and. idate2 > 31) then
   idate2=idate2-31
   imonth2=imonth2+1
endif
if(imonth2 == 8 .and. idate2 > 31) then
   idate2=idate2-31
   imonth2=imonth2+1
endif
if(imonth2 == 9 .and. idate2 > 30) then
   idate2=idate2-30
   imonth2=imonth2+1
endif
if(imonth2 == 10 .and. idate2 > 31) then
   idate2=idate2-31
   imonth2=imonth2+1
endif
if(imonth2 == 11 .and. idate2 > 30) then
   idate2=idate2-30
   imonth2=imonth2+1
endif
if(imonth2 == 12 .and. idate2 > 31) then
   idate2=idate2-31
   imonth2=imonth2+1
   if(imonth2>12) then
     iyear2=iyear2+1
     imonth2=1
   endif
endif

if(fcstsec1 > 0.0) print*,'BAD TIME PROGRESSION: TEST CODE'

if(iyear2 >= 2000) iyr2=iyear2-2000
if(iyear2 <  2000) iyr2=iyear2-1900

!EMPTY VARIABLES - 1 variables
if(len_trim(cvar) .eq. 7  .and. cvar(1:7)  .eq. 'empty3d')             cdnamegem='EMT3'

!3D VELOCITY AND VORTICITY VARIABLES - 21 variables
if(len_trim(cvar) .eq. 1  .and. cvar(1:1)  .eq. 'u')                   cdnamegem='UWND'
if(len_trim(cvar) .eq. 1  .and. cvar(1:1)  .eq. 'v')                   cdnamegem='VWND'
if(len_trim(cvar) .eq. 5  .and. cvar(1:5)  .eq. 'u_avg')               cdnamegem='UWDA' 
if(len_trim(cvar) .eq. 5  .and. cvar(1:5)  .eq. 'v_avg')               cdnamegem='VWDA'
if(len_trim(cvar) .eq. 2  .and. cvar(1:2)  .eq. 'ue')                  cdnamegem='UEWD'
if(len_trim(cvar) .eq. 2  .and. cvar(1:2)  .eq. 've')                  cdnamegem='VEWD'
if(len_trim(cvar) .eq. 6  .and. cvar(1:6)  .eq. 'ue_avg')              cdnamegem='UEWA'
if(len_trim(cvar) .eq. 6  .and. cvar(1:6)  .eq. 've_avg')              cdnamegem='VEWA'
if(len_trim(cvar) .eq. 1  .and. cvar(1:1)  .eq. 'w')                   cdnamegem='WWND'
if(len_trim(cvar) .eq. 4  .and. cvar(1:4)  .eq. 'wcms')                cdnamegem='WCMS'
if(len_trim(cvar) .eq. 5  .and. cvar(1:5)  .eq. 'w_avg')               cdnamegem='WAVG'
if(len_trim(cvar) .eq. 5  .and. cvar(1:5)  .eq. 'speed')               cdnamegem='SPED'
if(len_trim(cvar) .eq. 9  .and. cvar(1:9)  .eq. 'speed_mph')           cdnamegem='SMPH'
if(len_trim(cvar) .eq. 8  .and. cvar(1:8)  .eq. 'speed10m')            cdnamegem='SP10'
if(len_trim(cvar) .eq. 9  .and. cvar(1:9)  .eq. 'direction')           cdnamegem='DRCT'
if(len_trim(cvar) .eq. 8  .and. cvar(1:8)  .eq. 'relvortx')            cdnamegem='XVOR'
if(len_trim(cvar) .eq. 8  .and. cvar(1:8)  .eq. 'relvorty')            cdnamegem='YVOR'
if(len_trim(cvar) .eq. 8  .and. cvar(1:8)  .eq. 'relvortz')            cdnamegem='ZVOR'
if(len_trim(cvar) .eq. 8  .and. cvar(1:8)  .eq. 'absvortz')            cdnamegem='AVOR'
if(len_trim(cvar) .eq. 8  .and. cvar(1:8)  .eq. 'potvortz')            cdnamegem='PVOR'
if(len_trim(cvar) .eq. 9  .and. cvar(1:9)  .eq. 'horiz_div')           cdnamegem='HDIV'

!3D THERMODYNAMIC PROPERTIES OF AIR - 18 variables
if(len_trim(cvar) .eq. 2  .and. cvar(1:2)  .eq. 'pi')                  cdnamegem='XNER'
if(len_trim(cvar) .eq. 5  .and. cvar(1:5)  .eq. 'press')               cdnamegem='PRES'
if(len_trim(cvar) .eq. 6  .and. cvar(1:6)  .eq. 'pprime')              cdnamegem='PPRM'
if(len_trim(cvar) .eq. 8  .and. cvar(1:8)  .eq. 'theta_il')            cdnamegem='THIL'
if(len_trim(cvar) .eq. 5  .and. cvar(1:5)  .eq. 'theta')               cdnamegem='THTA'
if(len_trim(cvar) .eq. 3  .and. cvar(1:3)  .eq. 'dn0')                 cdnamegem='DEN0'
if(len_trim(cvar) .eq. 3  .and. cvar(1:3)  .eq. 'pi0')                 cdnamegem='XNR0'
if(len_trim(cvar) .eq. 3  .and. cvar(1:3)  .eq. 'th0')                 cdnamegem='THV0'
if(len_trim(cvar) .eq. 13 .and. cvar(1:13) .eq. 'pert_pressure')       cdnamegem='PERT'
if(len_trim(cvar) .eq. 5  .and. cvar(1:5)  .eq. 'tempk')               cdnamegem='TMPK'
if(len_trim(cvar) .eq. 5  .and. cvar(1:5)  .eq. 'tempc')               cdnamegem='TMPC'
if(len_trim(cvar) .eq. 5  .and. cvar(1:5)  .eq. 'tempf')               cdnamegem='TMPF'
if(len_trim(cvar) .eq. 7  .and. cvar(1:7)  .eq. 'theta_e')             cdnamegem='THTE'
if(len_trim(cvar) .eq. 7  .and. cvar(1:7)  .eq. 'theta_v')             cdnamegem='THTV'
if(len_trim(cvar) .eq. 9  .and. cvar(1:9)  .eq. 'theta_rho')           cdnamegem='THTR'
if(len_trim(cvar) .eq. 15 .and. cvar(1:15) .eq. 'buoyancy_liquid')     cdnamegem='BOYL'
if(len_trim(cvar) .eq. 7  .and. cvar(1:7)  .eq. 'tempf2m')             cdnamegem='TMPF2'
if(len_trim(cvar) .eq. 7  .and. cvar(1:7)  .eq. 'tempc2m')             cdnamegem='TMPC2'

!3D MOISTURE MASS MIXING RATIOS AND HUMIDITY - 37 variables
if(len_trim(cvar) .eq. 10 .and. cvar(1:10) .eq. 'vapr_press')          cdnamegem='VPRS'
if(len_trim(cvar) .eq. 4  .and. cvar(1:4)  .eq. 'rslf')                cdnamegem='RSLF'
if(len_trim(cvar) .eq. 4  .and. cvar(1:4)  .eq. 'rsif')                cdnamegem='RSIF'
if(len_trim(cvar) .eq. 5  .and. cvar(1:5)  .eq. 'vapor')               cdnamegem='VMIX'
if(len_trim(cvar) .eq. 5  .and. cvar(1:5)  .eq. 'cloud')               cdnamegem='CMIX'
if(len_trim(cvar) .eq. 8  .and. cvar(1:8)  .eq. 'cloud_m3')            cdnamegem='CMXV'
if(len_trim(cvar) .eq. 4  .and. cvar(1:4)  .eq. 'rain')                cdnamegem='RMIX'
if(len_trim(cvar) .eq. 7  .and. cvar(1:7)  .eq. 'rain_m3')             cdnamegem='RMXV'
if(len_trim(cvar) .eq. 8  .and. cvar(1:8)  .eq. 'pristine')            cdnamegem='PMIX'
if(len_trim(cvar) .eq. 11 .and. cvar(1:11) .eq. 'pristine_m3')         cdnamegem='PMXV'
if(len_trim(cvar) .eq. 4  .and. cvar(1:4)  .eq. 'snow')                cdnamegem='SMIX'
if(len_trim(cvar) .eq. 7  .and. cvar(1:7)  .eq. 'snow_m3')             cdnamegem='SMXV'
if(len_trim(cvar) .eq. 10 .and. cvar(1:10) .eq. 'aggregates')          cdnamegem='AMIX'
if(len_trim(cvar) .eq. 13 .and. cvar(1:13) .eq. 'aggregates_m3')       cdnamegem='AMXV'
if(len_trim(cvar) .eq. 7  .and. cvar(1:7)  .eq. 'graupel')             cdnamegem='GMIX'
if(len_trim(cvar) .eq. 10 .and. cvar(1:10) .eq. 'graupel_m3')          cdnamegem='GMXV'
if(len_trim(cvar) .eq. 4  .and. cvar(1:4)  .eq. 'hail')                cdnamegem='HMIX'
if(len_trim(cvar) .eq. 7  .and. cvar(1:7)  .eq. 'hail_m3')             cdnamegem='HMXV'
if(len_trim(cvar) .eq. 7  .and. cvar(1:7)  .eq. 'drizzle')             cdnamegem='DMIX'
if(len_trim(cvar) .eq. 10 .and. cvar(1:10) .eq. 'drizzle_m3')          cdnamegem='DMXV'
if(len_trim(cvar) .eq. 11 .and. cvar(1:11) .eq. 'prissnowagg')         cdnamegem='PSAM'
if(len_trim(cvar) .eq. 9  .and. cvar(1:9)  .eq. 'grauphail')           cdnamegem='GHMX'
if(len_trim(cvar) .eq. 6  .and. cvar(1:6)  .eq. 'liquid')              cdnamegem='LMIX'
if(len_trim(cvar) .eq. 3  .and. cvar(1:3)  .eq. 'ice')                 cdnamegem='IMIX'
if(len_trim(cvar) .eq. 18 .and. cvar(1:18) .eq. 'ctop_tempc_sstbase')  cdnamegem='CTST'
if(len_trim(cvar) .eq. 17 .and. cvar(1:17) .eq. 'ctop_tempc_nobase')   cdnamegem='CTOP'
if(len_trim(cvar) .eq. 10 .and. cvar(1:10) .eq. 'total_cond')          cdnamegem='TMIX'
if(len_trim(cvar) .eq. 13 .and. cvar(1:13) .eq. 'total_cond_m3')       cdnamegem='TMXV'
if(len_trim(cvar) .eq. 7  .and. cvar(1:7)  .eq. 'r_total')             cdnamegem='MIXR'
if(len_trim(cvar) .eq. 11 .and. cvar(1:11) .eq. 'rtotal_orig')         cdnamegem='MIXR'
if(len_trim(cvar) .eq. 6  .and. cvar(1:6)  .eq. 'dewptk')              cdnamegem='DWPK'
if(len_trim(cvar) .eq. 6  .and. cvar(1:6)  .eq. 'dewptf')              cdnamegem='DWPF'
if(len_trim(cvar) .eq. 6  .and. cvar(1:6)  .eq. 'dewptc')              cdnamegem='DWPC'
if(len_trim(cvar) .eq. 6  .and. cvar(1:6)  .eq. 'relhum')              cdnamegem='RELH'
if(len_trim(cvar) .eq. 11 .and. cvar(1:11) .eq. 'relhum_frac')         cdnamegem='RHFR'
if(len_trim(cvar) .eq. 10 .and. cvar(1:10) .eq. 'clear_frac')          cdnamegem='CLRF'
if(len_trim(cvar) .eq. 10 .and. cvar(1:10) .eq. 'cloud_frac')          cdnamegem='CLDF'

!3D HYDROMETEOR NUMBER CONCENTRATIONS - 22 variables
if(len_trim(cvar) .eq. 15 .and. cvar(1:15) .eq. 'cloud_concen_mg')     cdnamegem='CNMG'
if(len_trim(cvar) .eq. 15 .and. cvar(1:15) .eq. 'cloud_concen_kg')     cdnamegem='CNKG'
if(len_trim(cvar) .eq. 14 .and. cvar(1:14) .eq. 'rain_concen_kg')      cdnamegem='RNKG'
if(len_trim(cvar) .eq. 14 .and. cvar(1:14) .eq. 'pris_concen_mg')      cdnamegem='PNMG'
if(len_trim(cvar) .eq. 14 .and. cvar(1:14) .eq. 'pris_concen_kg')      cdnamegem='PNKG'
if(len_trim(cvar) .eq. 14 .and. cvar(1:14) .eq. 'snow_concen_kg')      cdnamegem='SNKG'
if(len_trim(cvar) .eq. 13 .and. cvar(1:13) .eq. 'agg_concen_kg')       cdnamegem='ANKG'
if(len_trim(cvar) .eq. 15 .and. cvar(1:15) .eq. 'graup_concen_kg')     cdnamegem='GNKG'
if(len_trim(cvar) .eq. 14 .and. cvar(1:14) .eq. 'hail_concen_kg')      cdnamegem='HNKG'
if(len_trim(cvar) .eq. 17 .and. cvar(1:17) .eq. 'drizzle_concen_mg')   cdnamegem='DNMG'
if(len_trim(cvar) .eq. 17 .and. cvar(1:17) .eq. 'drizzle_concen_kg')   cdnamegem='DNKG'
if(len_trim(cvar) .eq. 16 .and. cvar(1:16) .eq. 'cloud_concen_cm3')    cdnamegem='CNC3'
if(len_trim(cvar) .eq. 14 .and. cvar(1:14) .eq. 'rain_concen_m3')      cdnamegem='RNM3'
if(len_trim(cvar) .eq. 15 .and. cvar(1:15) .eq. 'rain_concen_dm3')     cdnamegem='RND3'
if(len_trim(cvar) .eq. 14 .and. cvar(1:14) .eq. 'pris_concen_m3')      cdnamegem='PNM3'
if(len_trim(cvar) .eq. 15 .and. cvar(1:15) .eq. 'pris_concen_cm3')     cdnamegem='PNC3'
if(len_trim(cvar) .eq. 14 .and. cvar(1:14) .eq. 'snow_concen_m3')      cdnamegem='SNM3'
if(len_trim(cvar) .eq. 15 .and. cvar(1:15) .eq. 'snow_concen_cm3')     cdnamegem='SNC3'
if(len_trim(cvar) .eq. 13 .and. cvar(1:13) .eq. 'agg_concen_m3')       cdnamegem='ANM3'
if(len_trim(cvar) .eq. 15 .and. cvar(1:15) .eq. 'graup_concen_m3')     cdnamegem='GNM3'
if(len_trim(cvar) .eq. 14 .and. cvar(1:14) .eq. 'hail_concen_m3')      cdnamegem='HNM3'
if(len_trim(cvar) .eq. 18 .and. cvar(1:18) .eq. 'drizzle_concen_cm3')  cdnamegem='DNC3'

!HUCM-SBM SPECIFIC MICROPHYSICS – 18 variables
if(len_trim(cvar) .eq. 10 .and. cvar(1:10) .eq. 'ice_plates')          cdnamegem='IPMX'
if(len_trim(cvar) .eq. 11 .and. cvar(1:11) .eq. 'ice_columns')         cdnamegem='ICMX'
if(len_trim(cvar) .eq. 13 .and. cvar(1:13) .eq. 'ice_dendrites')       cdnamegem='IDMX'
if(len_trim(cvar) .eq. 16 .and. cvar(1:16) .eq. 'plates_concen_mg')    cdnamegem='PCMG'
if(len_trim(cvar) .eq. 16 .and. cvar(1:16) .eq. 'plates_concen_kg')    cdnamegem='PCKG'
if(len_trim(cvar) .eq. 17 .and. cvar(1:17) .eq. 'columns_concen_mg')   cdnamegem='CCMG'
if(len_trim(cvar) .eq. 17 .and. cvar(1:17) .eq. 'columns_concen_kg')   cdnamegem='CCKG'
if(len_trim(cvar) .eq. 19 .and. cvar(1:19) .eq. 'dendrites_concen_mg') cdnamegem='DCMG'
if(len_trim(cvar) .eq. 19 .and. cvar(1:19) .eq. 'dendrites_concen_kg') cdnamegem='DCKG'
if(len_trim(cvar) .eq. 6  .and. cvar(1:6)  .eq. 'pcpvip')              cdnamegem='PVIP'
if(len_trim(cvar) .eq. 6  .and. cvar(1:6)  .eq. 'pcpvic')              cdnamegem='PVIC'
if(len_trim(cvar) .eq. 6  .and. cvar(1:6)  .eq. 'pcpvid')              cdnamegem='PVID'
if(len_trim(cvar) .eq. 6  .and. cvar(1:6)  .eq. 'pcprip')              cdnamegem='PRIP'
if(len_trim(cvar) .eq. 6  .and. cvar(1:6)  .eq. 'pcpric')              cdnamegem='PRIC'
if(len_trim(cvar) .eq. 6  .and. cvar(1:6)  .eq. 'pcprid')              cdnamegem='PRID'
if(len_trim(cvar) .eq. 6  .and. cvar(1:6)  .eq. 'accpip')              cdnamegem='ACIP'
if(len_trim(cvar) .eq. 6  .and. cvar(1:6)  .eq. 'accpic')              cdnamegem='ACIC'
if(len_trim(cvar) .eq. 6  .and. cvar(1:6)  .eq. 'accpid')              cdnamegem='ACID'

!3D AEROSOLS NUMBER, MASS, SIZE, SOLUBILITY - 35 variables
if(len_trim(cvar) .eq. 13 .and. cvar(1:13) .eq. 'ifn_concen_mg')       cdnamegem='IFNM'
if(len_trim(cvar) .eq. 14 .and. cvar(1:14) .eq. 'ifn_concen_cm3')      cdnamegem='IFNC'
if(len_trim(cvar) .eq. 13 .and. cvar(1:13) .eq. 'ccn_concen_mg')       cdnamegem='CCNM'
if(len_trim(cvar) .eq. 14 .and. cvar(1:14) .eq. 'ccn_concen_cm3')      cdnamegem='CCNC'
if(len_trim(cvar) .eq. 14 .and. cvar(1:14) .eq. 'gccn_concen_mg')      cdnamegem='GCNM'
if(len_trim(cvar) .eq. 15 .and. cvar(1:15) .eq. 'gccn_concen_cm3')     cdnamegem='GCNC'
if(len_trim(cvar) .eq. 12 .and. cvar(1:12) .eq. 'dust1_concen')        cdnamegem='D1CN'
if(len_trim(cvar) .eq. 12 .and. cvar(1:12) .eq. 'dust2_concen')        cdnamegem='D2CN'
if(len_trim(cvar) .eq. 16 .and. cvar(1:16) .eq. 'salt_film_concen')    cdnamegem='SFCN'
if(len_trim(cvar) .eq. 15 .and. cvar(1:15) .eq. 'salt_jet_concen')     cdnamegem='SJCN'
if(len_trim(cvar) .eq. 17 .and. cvar(1:17) .eq. 'salt_spume_concen')   cdnamegem='SSCN'
if(len_trim(cvar) .eq. 18 .and. cvar(1:18) .eq. 'regen_aero1_concen')  cdnamegem='R1CN'
if(len_trim(cvar) .eq. 18 .and. cvar(1:18) .eq. 'regen_aero2_concen')  cdnamegem='R2CN'
if(len_trim(cvar) .eq. 8  .and. cvar(1:8)  .eq. 'ccn_mass')            cdnamegem='CCCM'
if(len_trim(cvar) .eq. 9  .and. cvar(1:9)  .eq. 'gccn_mass')           cdnamegem='GCCM'
if(len_trim(cvar) .eq. 10 .and. cvar(1:10) .eq. 'dust1_mass')          cdnamegem='D1CM'
if(len_trim(cvar) .eq. 10 .and. cvar(1:10) .eq. 'dust2_mass')          cdnamegem='D2CM'
if(len_trim(cvar) .eq. 14 .and. cvar(1:14) .eq. 'salt_film_mass')      cdnamegem='SFCM'
if(len_trim(cvar) .eq. 13 .and. cvar(1:13) .eq. 'salt_jet_mass')       cdnamegem='SJCM'
if(len_trim(cvar) .eq. 15 .and. cvar(1:15) .eq. 'salt_spume_mass')     cdnamegem='SSCM'
if(len_trim(cvar) .eq. 16 .and. cvar(1:16) .eq. 'regen_aero1_mass')    cdnamegem='R1CM'
if(len_trim(cvar) .eq. 16 .and. cvar(1:16) .eq. 'regen_aero2_mass')    cdnamegem='R2CM'
if(len_trim(cvar) .eq. 16 .and. cvar(1:16) .eq. 'resol_aero1_mass')    cdnamegem='R1SM'
if(len_trim(cvar) .eq. 16 .and. cvar(1:16) .eq. 'resol_aero2_mass')    cdnamegem='R2SM'
if(len_trim(cvar) .eq. 14 .and. cvar(1:14) .eq. 'regen1_epsilon')      cdnamegem='R1EP'
if(len_trim(cvar) .eq. 14 .and. cvar(1:14) .eq. 'regen2_epsilon')      cdnamegem='R2EP'
if(len_trim(cvar) .eq. 10 .and. cvar(1:10) .eq. 'ccn_medrad')          cdnamegem='CCCR'
if(len_trim(cvar) .eq. 11 .and. cvar(1:11) .eq. 'gccn_medrad')         cdnamegem='GCCR'
if(len_trim(cvar) .eq. 12 .and. cvar(1:12) .eq. 'dust1_medrad')        cdnamegem='D1CR'
if(len_trim(cvar) .eq. 12 .and. cvar(1:12) .eq. 'dust2_medrad')        cdnamegem='D2CR'
if(len_trim(cvar) .eq. 16 .and. cvar(1:16) .eq. 'salt_film_medrad')    cdnamegem='SFCR'
if(len_trim(cvar) .eq. 15 .and. cvar(1:15) .eq. 'salt_jet_medrad')     cdnamegem='SJCR'
if(len_trim(cvar) .eq. 17 .and. cvar(1:17) .eq. 'salt_spume_medrad')   cdnamegem='SSCR'
if(len_trim(cvar) .eq. 18 .and. cvar(1:18) .eq. 'regen_aero1_medrad')  cdnamegem='R1CR'
if(len_trim(cvar) .eq. 18 .and. cvar(1:18) .eq. 'regen_aero2_medrad')  cdnamegem='R2CR'

!3D AEROSOL TRACKING VARIABLES - 41 variables
if(len_trim(cvar) .eq. 18 .and. cvar(1:18) .eq. 'aerosol_cloud_mass')  cdnamegem='ARMC'
if(len_trim(cvar) .eq. 17 .and. cvar(1:17) .eq. 'aerosol_rain_mass')   cdnamegem='ARMR'
if(len_trim(cvar) .eq. 17 .and. cvar(1:17) .eq. 'aerosol_pris_mass')   cdnamegem='ARMP'
if(len_trim(cvar) .eq. 17 .and. cvar(1:17) .eq. 'aerosol_snow_mass')   cdnamegem='ARMS'
if(len_trim(cvar) .eq. 17 .and. cvar(1:17) .eq. 'aerosol_aggr_mass')   cdnamegem='ARMA'
if(len_trim(cvar) .eq. 17 .and. cvar(1:17) .eq. 'aerosol_grau_mass')   cdnamegem='ARMG'
if(len_trim(cvar) .eq. 17 .and. cvar(1:17) .eq. 'aerosol_hail_mass')   cdnamegem='ARMH'
if(len_trim(cvar) .eq. 17 .and. cvar(1:17) .eq. 'aerosol_driz_mass')   cdnamegem='ARMD'
if(len_trim(cvar) .eq. 18 .and. cvar(1:18) .eq. 'aerosol_hydro_mass')  cdnamegem='ARHY'
if(len_trim(cvar) .eq. 18 .and. cvar(1:18) .eq. 'soluble_cloud_mass')  cdnamegem='SLMC'
if(len_trim(cvar) .eq. 17 .and. cvar(1:17) .eq. 'soluble_rain_mass')   cdnamegem='SLMR'
if(len_trim(cvar) .eq. 17 .and. cvar(1:17) .eq. 'soluble_pris_mass')   cdnamegem='SLMP'
if(len_trim(cvar) .eq. 17 .and. cvar(1:17) .eq. 'soluble_snow_mass')   cdnamegem='SLMS'
if(len_trim(cvar) .eq. 17 .and. cvar(1:17) .eq. 'soluble_aggr_mass')   cdnamegem='SLMA'
if(len_trim(cvar) .eq. 17 .and. cvar(1:17) .eq. 'soluble_grau_mass')   cdnamegem='SLMG'
if(len_trim(cvar) .eq. 17 .and. cvar(1:17) .eq. 'soluble_hail_mass')   cdnamegem='SLMH'
if(len_trim(cvar) .eq. 17 .and. cvar(1:17) .eq. 'soluble_driz_mass')   cdnamegem='SLMD'
if(len_trim(cvar) .eq. 18 .and. cvar(1:18) .eq. 'soluble_hydro_mass')  cdnamegem='SLHY'
if(len_trim(cvar) .eq. 12 .and. cvar(1:12) .eq. 'aero_epsilon')        cdnamegem='EPSI'
if(len_trim(cvar) .eq. 15 .and. cvar(1:15) .eq. 'dust_cloud_mass')     cdnamegem='DUMC'
if(len_trim(cvar) .eq. 14 .and. cvar(1:14) .eq. 'dust_rain_mass')      cdnamegem='DUMR'
if(len_trim(cvar) .eq. 14 .and. cvar(1:14) .eq. 'dust_pris_mass')      cdnamegem='DUMP'
if(len_trim(cvar) .eq. 14 .and. cvar(1:14) .eq. 'dust_snow_mass')      cdnamegem='DUMS'
if(len_trim(cvar) .eq. 14 .and. cvar(1:14) .eq. 'dust_aggr_mass')      cdnamegem='DUMA'
if(len_trim(cvar) .eq. 14 .and. cvar(1:14) .eq. 'dust_grau_mass')      cdnamegem='DUMG'
if(len_trim(cvar) .eq. 14 .and. cvar(1:14) .eq. 'dust_hail_mass')      cdnamegem='DUMH'
if(len_trim(cvar) .eq. 14 .and. cvar(1:14) .eq. 'dust_driz_mass')      cdnamegem='DUMD'
if(len_trim(cvar) .eq. 15 .and. cvar(1:15) .eq. 'dust_hydro_mass')     cdnamegem='DUHY'
if(len_trim(cvar) .eq. 18 .and. cvar(1:18) .eq. 'dustifn_cloud_mass')  cdnamegem='DINC'
if(len_trim(cvar) .eq. 17 .and. cvar(1:17) .eq. 'dustifn_rain_mass')   cdnamegem='DINR'
if(len_trim(cvar) .eq. 17 .and. cvar(1:17) .eq. 'dustifn_pris_mass')   cdnamegem='DINP'
if(len_trim(cvar) .eq. 17 .and. cvar(1:17) .eq. 'dustifn_snow_mass')   cdnamegem='DINS'
if(len_trim(cvar) .eq. 17 .and. cvar(1:17) .eq. 'dustifn_aggr_mass')   cdnamegem='DINA'
if(len_trim(cvar) .eq. 17 .and. cvar(1:17) .eq. 'dustifn_grau_mass')   cdnamegem='DING'
if(len_trim(cvar) .eq. 17 .and. cvar(1:17) .eq. 'dustifn_hail_mass')   cdnamegem='DINH'
if(len_trim(cvar) .eq. 17 .and. cvar(1:17) .eq. 'dustifn_driz_mass')   cdnamegem='DIND'
if(len_trim(cvar) .eq. 18 .and. cvar(1:18) .eq. 'dustifn_hydro_mass')  cdnamegem='DIHY'
if(len_trim(cvar) .eq. 16 .and. cvar(1:16) .eq. 'ifn_nuc_numtrack')    cdnamegem='INTR'
if(len_trim(cvar) .eq. 11 .and. cvar(1:11) .eq. 'ifn_incloud')         cdnamegem='CICN'
if(len_trim(cvar) .eq. 10 .and. cvar(1:10) .eq. 'ifn_indriz')          cdnamegem='DICN'
if(len_trim(cvar) .eq. 10 .and. cvar(1:10) .eq. 'ifn_inrain')          cdnamegem='RICN'

!3D VERTICAL VELOCITY AND HYDROMETEOR BUDGETS - 50 variables
if(len_trim(cvar) .eq. 9  .and. cvar(1:9)  .eq. 'wp_advdif')           cdnamegem='WPAD'
if(len_trim(cvar) .eq. 13 .and. cvar(1:13) .eq. 'wp_buoy_theta')       cdnamegem='WPTH'
if(len_trim(cvar) .eq. 12 .and. cvar(1:12) .eq. 'wp_buoy_cond')        cdnamegem='WPCD'
if(len_trim(cvar) .eq. 8  .and. cvar(1:8)  .eq. 'nuccldrt')            cdnamegem='NUCRT'
if(len_trim(cvar) .eq. 9  .and. cvar(1:9)  .eq. 'cld2raint')           cdnamegem='CL2RT'
if(len_trim(cvar) .eq. 9  .and. cvar(1:9)  .eq. 'ice2raint')           cdnamegem='IC2RT'
if(len_trim(cvar) .eq. 8  .and. cvar(1:8)  .eq. 'nucicert')            cdnamegem='NUIRT'
if(len_trim(cvar) .eq. 7  .and. cvar(1:7)  .eq. 'vapliqt')             cdnamegem='VAPLT'
if(len_trim(cvar) .eq. 7  .and. cvar(1:7)  .eq. 'vapicet')             cdnamegem='VAPIT'
if(len_trim(cvar) .eq. 8  .and. cvar(1:8)  .eq. 'melticet')            cdnamegem='MELTT'
if(len_trim(cvar) .eq. 8  .and. cvar(1:8)  .eq. 'rimecldt')            cdnamegem='RIMCT'
if(len_trim(cvar) .eq. 9  .and. cvar(1:9)  .eq. 'rain2icet')           cdnamegem='R2ICT'
if(len_trim(cvar) .eq. 10 .and. cvar(1:10) .eq. 'aggregatet')          cdnamegem='AGGRT'
if(len_trim(cvar) .eq. 10 .and. cvar(1:10) .eq. 'latheatvap')          cdnamegem='LHVP'
if(len_trim(cvar) .eq. 11 .and. cvar(1:11) .eq. 'latheatvapt')         cdnamegem='LHVPT'
if(len_trim(cvar) .eq. 10 .and. cvar(1:10) .eq. 'latheatfrz')          cdnamegem='LHFZ'
if(len_trim(cvar) .eq. 11 .and. cvar(1:11) .eq. 'latheatfrzt')         cdnamegem='LHFZT'
if(len_trim(cvar) .eq. 9  .and. cvar(1:9)  .eq. 'inuchomrt')           cdnamegem='IHMRT'
if(len_trim(cvar) .eq. 10 .and. cvar(1:10) .eq. 'inuccontrt')          cdnamegem='ICORT'
if(len_trim(cvar) .eq. 9  .and. cvar(1:9)  .eq. 'inucifnrt')           cdnamegem='IINRT'
if(len_trim(cvar) .eq. 9  .and. cvar(1:9)  .eq. 'inuchazrt')           cdnamegem='IHZRT'
if(len_trim(cvar) .eq. 7  .and. cvar(1:7)  .eq. 'vapcldt')             cdnamegem='VAPCT'
if(len_trim(cvar) .eq. 8  .and. cvar(1:8)  .eq. 'vapraint')            cdnamegem='VAPRT'
if(len_trim(cvar) .eq. 8  .and. cvar(1:8)  .eq. 'vapprist')            cdnamegem='VAPPT'
if(len_trim(cvar) .eq. 8  .and. cvar(1:8)  .eq. 'vapsnowt')            cdnamegem='VAPST'
if(len_trim(cvar) .eq. 8  .and. cvar(1:8)  .eq. 'vapaggrt')            cdnamegem='VAPAT'
if(len_trim(cvar) .eq. 8  .and. cvar(1:8)  .eq. 'vapgraut')            cdnamegem='VAPGT'
if(len_trim(cvar) .eq. 8  .and. cvar(1:8)  .eq. 'vaphailt')            cdnamegem='VAPHT'
if(len_trim(cvar) .eq. 8  .and. cvar(1:8)  .eq. 'vapdrizt')            cdnamegem='VAPDT'
if(len_trim(cvar) .eq. 9  .and. cvar(1:9)  .eq. 'meltprist')           cdnamegem='MELPT'
if(len_trim(cvar) .eq. 9  .and. cvar(1:9)  .eq. 'meltsnowt')           cdnamegem='MELST'
if(len_trim(cvar) .eq. 9  .and. cvar(1:9)  .eq. 'meltaggrt')           cdnamegem='MELAT'
if(len_trim(cvar) .eq. 9  .and. cvar(1:9)  .eq. 'meltgraut')           cdnamegem='MELGT'
if(len_trim(cvar) .eq. 9  .and. cvar(1:9)  .eq. 'melthailt')           cdnamegem='MELHT'
if(len_trim(cvar) .eq. 12 .and. cvar(1:12) .eq. 'rimecldsnowt')        cdnamegem='RIMST'
if(len_trim(cvar) .eq. 12 .and. cvar(1:12) .eq. 'rimecldaggrt')        cdnamegem='RIMAT'
if(len_trim(cvar) .eq. 12 .and. cvar(1:12) .eq. 'rimecldgraut')        cdnamegem='RIMGT'
if(len_trim(cvar) .eq. 12 .and. cvar(1:12) .eq. 'rimecldhailt')        cdnamegem='RIMHT'
if(len_trim(cvar) .eq. 8  .and. cvar(1:8)  .eq. 'rain2prt')            cdnamegem='R2PRT'
if(len_trim(cvar) .eq. 8  .and. cvar(1:8)  .eq. 'rain2snt')            cdnamegem='R2SNT'
if(len_trim(cvar) .eq. 8  .and. cvar(1:8)  .eq. 'rain2agt')            cdnamegem='R2AGT'
if(len_trim(cvar) .eq. 8  .and. cvar(1:8)  .eq. 'rain2grt')            cdnamegem='R2GRT'
if(len_trim(cvar) .eq. 8  .and. cvar(1:8)  .eq. 'rain2hat')            cdnamegem='R2HAT'
if(len_trim(cvar) .eq. 13 .and. cvar(1:13) .eq. 'aggrselfprist')       cdnamegem='AGPPT'
if(len_trim(cvar) .eq. 13 .and. cvar(1:13) .eq. 'aggrselfsnowt')       cdnamegem='AGSST'
if(len_trim(cvar) .eq. 13 .and. cvar(1:13) .eq. 'aggrprissnowt')       cdnamegem='AGPST'
if(len_trim(cvar) .eq. 10 .and. cvar(1:10) .eq. 'dust1cldrt')          cdnamegem='D1CRT'
if(len_trim(cvar) .eq. 10 .and. cvar(1:10) .eq. 'dust2cldrt')          cdnamegem='D2CRT'
if(len_trim(cvar) .eq. 10 .and. cvar(1:10) .eq. 'dust1drzrt')          cdnamegem='D1DRT'
if(len_trim(cvar) .eq. 10 .and. cvar(1:10) .eq. 'dust2drzrt')          cdnamegem='D2DRT'

!3D HYDROMETEOR DIAMETERS - 9 variables
if(len_trim(cvar) .eq. 13 .and. cvar(1:13) .eq. 'cloudtop_diam')       cdnamegem='TDIM'
if(len_trim(cvar) .eq. 10 .and. cvar(1:10) .eq. 'cloud_diam')          cdnamegem='CDIM'
if(len_trim(cvar) .eq. 9  .and. cvar(1:9)  .eq. 'rain_diam')           cdnamegem='RDIM'
if(len_trim(cvar) .eq. 9  .and. cvar(1:9)  .eq. 'pris_diam')           cdnamegem='PDIM'
if(len_trim(cvar) .eq. 9  .and. cvar(1:9)  .eq. 'snow_diam')           cdnamegem='SDIM'
if(len_trim(cvar) .eq. 8  .and. cvar(1:8)  .eq. 'agg_diam')            cdnamegem='ADIM'
if(len_trim(cvar) .eq. 10 .and. cvar(1:10) .eq. 'graup_diam')          cdnamegem='GDIM'
if(len_trim(cvar) .eq. 9  .and. cvar(1:9)  .eq. 'hail_diam')           cdnamegem='HDIM'
if(len_trim(cvar) .eq. 12 .and. cvar(1:12) .eq. 'drizzle_diam')        cdnamegem='DDIM'

!3D HYDROMETEOR TEMPERATURE, ENERGY, LIQUID FRACTION - 11 variables
if(len_trim(cvar) .eq. 2  .and. cvar(1:2)  .eq. 'q2')                  cdnamegem='Q2RA'
if(len_trim(cvar) .eq. 2  .and. cvar(1:2)  .eq. 'q6')                  cdnamegem='Q6GR'
if(len_trim(cvar) .eq. 2  .and. cvar(1:2)  .eq. 'q7')                  cdnamegem='Q7HA'
if(len_trim(cvar) .eq. 9  .and. cvar(1:9)  .eq. 'rain_temp')           cdnamegem='RTMP'
if(len_trim(cvar) .eq. 10 .and. cvar(1:10) .eq. 'graup_temp')          cdnamegem='GTMP'
if(len_trim(cvar) .eq. 9  .and. cvar(1:9)  .eq. 'hail_temp')           cdnamegem='HTMP'
if(len_trim(cvar) .eq. 16 .and. cvar(1:16) .eq. 'rain_air_tempdif')    cdnamegem='RATD'
if(len_trim(cvar) .eq. 17 .and. cvar(1:17) .eq. 'graup_air_tempdif')   cdnamegem='GATD'
if(len_trim(cvar) .eq. 16 .and. cvar(1:16) .eq. 'hail_air_tempdif')    cdnamegem='HATD'
if(len_trim(cvar) .eq. 13 .and. cvar(1:13) .eq. 'graup_fracliq')       cdnamegem='GLIQ'
if(len_trim(cvar) .eq. 12 .and. cvar(1:12) .eq. 'hail_fracliq')        cdnamegem='HLIQ'

!3D MISCELLANEOUS FIELDS - 5 variables
if(len_trim(cvar) .eq. 3  .and. cvar(1:3)  .eq. 'geo')                 cdnamegem='HGHT'
if(len_trim(cvar) .eq. 3  .and. cvar(1:3)  .eq. 'tke')                 cdnamegem='TKET'
if(len_trim(cvar) .eq. 3  .and. cvar(1:3)  .eq. 'eps')                 cdnamegem='TKED'
if(len_trim(cvar) .eq. 6  .and. cvar(1:6)  .eq. 'pbl_ht')              cdnamegem='PBLH'
if(len_trim(cvar) .eq. 11 .and. cvar(1:11) .eq. 'reflect_all')         cdnamegem='DBZZ'

!CUMULUS PARM - RADIATION - TURBULENCE PARAMETERS - 11 variables
if(len_trim(cvar) .eq. 13 .and. cvar(1:13) .eq. 'cuparm_thetasrc')     cdnamegem='CVHR'
if(len_trim(cvar) .eq. 13 .and. cvar(1:13) .eq. 'cuparm_rtsrc')        cdnamegem='CVMR'
if(len_trim(cvar) .eq. 3  .and. cvar(1:3)  .eq. 'khh')                 cdnamegem='KHHC'
if(len_trim(cvar) .eq. 3  .and. cvar(1:3)  .eq. 'khv')                 cdnamegem='KHVC'
if(len_trim(cvar) .eq. 10 .and. cvar(1:10) .eq. 'visibility')          cdnamegem='VISB'
if(len_trim(cvar) .eq. 4  .and. cvar(1:4)  .eq. 'swup')                cdnamegem='SWUP'
if(len_trim(cvar) .eq. 4  .and. cvar(1:4)  .eq. 'swdn')                cdnamegem='SWDN'
if(len_trim(cvar) .eq. 4  .and. cvar(1:4)  .eq. 'lwup')                cdnamegem='LWUP'
if(len_trim(cvar) .eq. 4  .and. cvar(1:4)  .eq. 'lwdn')                cdnamegem='LWDN'
if(len_trim(cvar) .eq. 12 .and. cvar(1:12) .eq. 'rad_thetasrc')        cdnamegem='RAHR'
if(len_trim(cvar) .eq. 18 .and. cvar(1:18) .eq. 'column_net_rad_flx')  cdnamegem='NETR'

!2D SURFACE PRECIPITATION and VERTICALLY INTEGRATED FIELDS - 54 variables
if(len_trim(cvar) .eq. 5  .and. cvar(1:5)  .eq. 'accpr')               cdnamegem='ACCR'
if(len_trim(cvar) .eq. 5  .and. cvar(1:5)  .eq. 'accpp')               cdnamegem='ACCP'
if(len_trim(cvar) .eq. 5  .and. cvar(1:5)  .eq. 'accps')               cdnamegem='ACCS'
if(len_trim(cvar) .eq. 5  .and. cvar(1:5)  .eq. 'accpa')               cdnamegem='ACCA'
if(len_trim(cvar) .eq. 5  .and. cvar(1:5)  .eq. 'accpg')               cdnamegem='ACCG'
if(len_trim(cvar) .eq. 5  .and. cvar(1:5)  .eq. 'accph')               cdnamegem='ACCH'
if(len_trim(cvar) .eq. 5  .and. cvar(1:5)  .eq. 'accpd')               cdnamegem='ACCD'
if(len_trim(cvar) .eq. 8  .and. cvar(1:8)  .eq. 'accpaero')            cdnamegem='ACTA'
if(len_trim(cvar) .eq. 8  .and. cvar(1:8)  .eq. 'accpdust')            cdnamegem='ACDU'
if(len_trim(cvar) .eq. 6  .and. cvar(1:6)  .eq. 'totpcp')              cdnamegem='TRPM'
if(len_trim(cvar) .eq. 9  .and. cvar(1:9)  .eq. 'totpcp_in')           cdnamegem='TRPI'
if(len_trim(cvar) .eq. 6  .and. cvar(1:6)  .eq. 'precip')              cdnamegem='TAPM'
if(len_trim(cvar) .eq. 9  .and. cvar(1:9)  .eq. 'precip_in')           cdnamegem='TAPI'
if(len_trim(cvar) .eq. 5  .and. cvar(1:5)  .eq. 'pcprr')               cdnamegem='PCRR'
if(len_trim(cvar) .eq. 5  .and. cvar(1:5)  .eq. 'pcpvr')               cdnamegem='PCVR'
if(len_trim(cvar) .eq. 5  .and. cvar(1:5)  .eq. 'pcprp')               cdnamegem='PCRP'
if(len_trim(cvar) .eq. 5  .and. cvar(1:5)  .eq. 'pcpvp')               cdnamegem='PCVP'
if(len_trim(cvar) .eq. 5  .and. cvar(1:5)  .eq. 'pcprs')               cdnamegem='PCRS'
if(len_trim(cvar) .eq. 5  .and. cvar(1:5)  .eq. 'pcpvs')               cdnamegem='PCVS'
if(len_trim(cvar) .eq. 5  .and. cvar(1:5)  .eq. 'pcpra')               cdnamegem='PCRA'
if(len_trim(cvar) .eq. 5  .and. cvar(1:5)  .eq. 'pcpva')               cdnamegem='PCVA'
if(len_trim(cvar) .eq. 5  .and. cvar(1:5)  .eq. 'pcprg')               cdnamegem='PCRG'
if(len_trim(cvar) .eq. 5  .and. cvar(1:5)  .eq. 'pcpvg')               cdnamegem='PCVG'
if(len_trim(cvar) .eq. 5  .and. cvar(1:5)  .eq. 'pcprh')               cdnamegem='PCRH'
if(len_trim(cvar) .eq. 5  .and. cvar(1:5)  .eq. 'pcpvh')               cdnamegem='PCVH'
if(len_trim(cvar) .eq. 5  .and. cvar(1:5)  .eq. 'pcprd')               cdnamegem='PCRD'
if(len_trim(cvar) .eq. 5  .and. cvar(1:5)  .eq. 'pcpvd')               cdnamegem='PCVD'
if(len_trim(cvar) .eq. 4  .and. cvar(1:4)  .eq. 'pcpg')                cdnamegem='PCPG'
if(len_trim(cvar) .eq. 5  .and. cvar(1:5)  .eq. 'qpcpg')               cdnamegem='PCPQ'
if(len_trim(cvar) .eq. 5  .and. cvar(1:5)  .eq. 'dpcpg')               cdnamegem='PCPD'
if(len_trim(cvar) .eq. 7  .and. cvar(1:7)  .eq. 'pcprate')             cdnamegem='PRRM'
if(len_trim(cvar) .eq. 10 .and. cvar(1:10) .eq. 'pcprate_in')          cdnamegem='PRRI'
if(len_trim(cvar) .eq. 7  .and. cvar(1:7)  .eq. 'precipr')             cdnamegem='PRTM'
if(len_trim(cvar) .eq. 10 .and. cvar(1:10) .eq. 'precipr_in')          cdnamegem='PRTI'
if(len_trim(cvar) .eq. 6  .and. cvar(1:6)  .eq. 'conpcp')              cdnamegem='CNPR'
if(len_trim(cvar) .eq. 6  .and. cvar(1:6)  .eq. 'acccon')              cdnamegem='ACON'
if(len_trim(cvar) .eq. 9  .and. cvar(1:9)  .eq. 'vertmax_w')           cdnamegem='VMXW'
if(len_trim(cvar) .eq. 9  .and. cvar(1:9)  .eq. 'vertavg_w')           cdnamegem='VAVW'
if(len_trim(cvar) .eq. 12 .and. cvar(1:12) .eq. 'vertint_cond')        cdnamegem='COND'
if(len_trim(cvar) .eq. 10 .and. cvar(1:10) .eq. 'vertint_rt')          cdnamegem='WATR'
if(len_trim(cvar) .eq. 12 .and. cvar(1:12) .eq. 'vertint_orig')        cdnamegem='VERT'
if(len_trim(cvar) .eq. 13 .and. cvar(1:13) .eq. 'vertint_vapor')       cdnamegem='VRTV'
if(len_trim(cvar) .eq. 11 .and. cvar(1:11) .eq. 'vertint_liq')         cdnamegem='VRTL'
if(len_trim(cvar) .eq. 11 .and. cvar(1:11) .eq. 'vertint_ice')         cdnamegem='VRTI'
if(len_trim(cvar) .eq. 13 .and. cvar(1:13) .eq. 'vertint_cloud')       cdnamegem='VRTC'
if(len_trim(cvar) .eq. 12 .and. cvar(1:12) .eq. 'vertint_driz')        cdnamegem='VRTD'
if(len_trim(cvar) .eq. 12 .and. cvar(1:12) .eq. 'vertint_rain')        cdnamegem='VRTR'
if(len_trim(cvar) .eq. 12 .and. cvar(1:12) .eq. 'vertint_pris')        cdnamegem='VRTP'
if(len_trim(cvar) .eq. 12 .and. cvar(1:12) .eq. 'vertint_snow')        cdnamegem='VRTS'
if(len_trim(cvar) .eq. 12 .and. cvar(1:12) .eq. 'vertint_aggr')        cdnamegem='VRTA'
if(len_trim(cvar) .eq. 15 .and. cvar(1:15) .eq. 'vertint_graupel')     cdnamegem='VRTG'
if(len_trim(cvar) .eq. 12 .and. cvar(1:12) .eq. 'vertint_hail')        cdnamegem='VRTH'
if(len_trim(cvar) .eq. 12 .and. cvar(1:12) .eq. 'vertint_dust')        cdnamegem='VTDU'
if(len_trim(cvar) .eq. 18 .and. cvar(1:18) .eq. 'vertint_dust_hydro')  cdnamegem='VTDH'

!2D SEA ICE COVERAGE, DEPTH, ROUGHNESS, TEMP, SNOW COVER - 5 variables
if(len_trim(cvar) .eq. 14 .and. cvar(1:14) .eq. 'snowdepthonice')      cdnamegem='DEPS'
if(len_trim(cvar) .eq. 9  .and. cvar(1:9)  .eq. 'cicedepth')           cdnamegem='DEPI'
if(len_trim(cvar) .eq. 9  .and. cvar(1:9)  .eq. 'cicefract')           cdnamegem='ICEF'
if(len_trim(cvar) .eq. 8  .and. cvar(1:8)  .eq. 'cicetemp')            cdnamegem='ICET'
if(len_trim(cvar) .eq. 9  .and. cvar(1:9)  .eq. 'cicerough')           cdnamegem='ICER'

!2D SURFACE HEAT, MOISTURE, MOMENTUM AND RADIATIVE FLUXES - 12 variables
if(len_trim(cvar) .eq. 9  .and. cvar(1:9)  .eq. 'sens_flux')           cdnamegem='SFLX'
if(len_trim(cvar) .eq. 8  .and. cvar(1:8)  .eq. 'lat_flux')            cdnamegem='LFLX'
if(len_trim(cvar) .eq. 6  .and. cvar(1:6)  .eq. 'etrans')              cdnamegem='EVAP'
if(len_trim(cvar) .eq. 9  .and. cvar(1:9)  .eq. 'etrans_in')           cdnamegem='ETRI'
if(len_trim(cvar) .eq. 8  .and. cvar(1:8)  .eq. 'umom_flx')            cdnamegem='UFLX'
if(len_trim(cvar) .eq. 8  .and. cvar(1:8)  .eq. 'vmom_flx')            cdnamegem='VFLX'
if(len_trim(cvar) .eq. 8  .and. cvar(1:8)  .eq. 'wmom_flx')            cdnamegem='WFLX'
if(len_trim(cvar) .eq. 5  .and. cvar(1:5)  .eq. 'bowen')               cdnamegem='BOWN'
if(len_trim(cvar) .eq. 6  .and. cvar(1:6)  .eq. 'rshort')              cdnamegem='RSHT'
if(len_trim(cvar) .eq. 5  .and. cvar(1:5)  .eq. 'rlong')               cdnamegem='RLON'
if(len_trim(cvar) .eq. 7  .and. cvar(1:7)  .eq. 'rlongup')             cdnamegem='RLNU'
if(len_trim(cvar) .eq. 6  .and. cvar(1:6)  .eq. 'albedt')              cdnamegem='ALBE'

!2D TOPOGRAPHY AND GEOGRAPHIC VALUES - 3 variables
if(len_trim(cvar) .eq. 4  .and. cvar(1:4)  .eq. 'topt')                cdnamegem='TOPT'
if(len_trim(cvar) .eq. 3  .and. cvar(1:3)  .eq. 'lat')                 cdnamegem='LATI'
if(len_trim(cvar) .eq. 3  .and. cvar(1:3)  .eq. 'lon')                 cdnamegem='LONG'

!2D MISCELLANEOUS FIELDS - 3 variables
if(len_trim(cvar) .eq. 9  .and. cvar(1:9)  .eq. 'sea_press')           cdnamegem='MSLP'
if(len_trim(cvar) .eq. 7  .and. cvar(1:7)  .eq. 'sfc_div')             cdnamegem='SDIV'
if(len_trim(cvar) .eq. 3  .and. cvar(1:3)  .eq. 'sst')                 cdnamegem='SSTC'

!LEAF/SIB VARIABLES SECTION - 33 variables
if(len_trim(cvar) .eq. 10 .and. cvar(1:10) .eq. 'patch_area')          cdnamegem='PFRA'
if(len_trim(cvar) .eq. 4  .and. cvar(1:4)  .eq. 'land')                cdnamegem='LAND'
if(len_trim(cvar) .eq. 11 .and. cvar(1:11) .eq. 'snow_levels')         cdnamegem='SNOL'
if(len_trim(cvar) .eq. 13 .and. cvar(1:13) .eq. 'snow_depth_ps')       cdnamegem='SNOD'
if(len_trim(cvar) .eq. 12 .and. cvar(1:12) .eq. 'snow_mass_ps')        cdnamegem='SNOM'
if(len_trim(cvar) .eq. 12 .and. cvar(1:12) .eq. 'snow_temp_ps')        cdnamegem='SNOT'
if(len_trim(cvar) .eq. 10 .and. cvar(1:10) .eq. 'topo_z0_ps')          cdnamegem='TRUF'
if(len_trim(cvar) .eq. 9  .and. cvar(1:9)  .eq. 'net_z0_ps')           cdnamegem='NRUF'
if(len_trim(cvar) .eq. 10 .and. cvar(1:10) .eq. 'soil_z0_ps')          cdnamegem='SRUF'
if(len_trim(cvar) .eq. 9  .and. cvar(1:9)  .eq. 'veg_z0_ps')           cdnamegem='VRUF'
if(len_trim(cvar) .eq. 11 .and. cvar(1:11) .eq. 'veg_ndvi_ps')         cdnamegem='NDVI'
if(len_trim(cvar) .eq. 12 .and. cvar(1:12) .eq. 'veg_class_bp')        cdnamegem='VEGC'
if(len_trim(cvar) .eq. 13 .and. cvar(1:13) .eq. 'veg_albedo_ps')       cdnamegem='VEGA'
if(len_trim(cvar) .eq. 15 .and. cvar(1:15) .eq. 'veg_fracarea_ps')     cdnamegem='VEGF'
if(len_trim(cvar) .eq. 10 .and. cvar(1:10) .eq. 'veg_lai_ps')          cdnamegem='LAIF'
if(len_trim(cvar) .eq. 11 .and. cvar(1:11) .eq. 'veg_disp_ps')         cdnamegem='VDIS'
if(len_trim(cvar) .eq. 16 .and. cvar(1:16) .eq. 'canopy_mixrat_ps')    cdnamegem='CANM'
if(len_trim(cvar) .eq. 14 .and. cvar(1:14) .eq. 'grnd_mixrat_ps')      cdnamegem='GRDM'
if(len_trim(cvar) .eq. 14 .and. cvar(1:14) .eq. 'soil_mixrat_ps')      cdnamegem='SOIM'
if(len_trim(cvar) .eq. 12 .and. cvar(1:12) .eq. 'veg_moist_ps')        cdnamegem='VEGM'
if(len_trim(cvar) .eq. 11 .and. cvar(1:11) .eq. 'veg_temp_ps')         cdnamegem='VEGT'
if(len_trim(cvar) .eq. 15 .and. cvar(1:15) .eq. 'canopy_tempc_ps')     cdnamegem='CANC'
if(len_trim(cvar) .eq. 15 .and. cvar(1:15) .eq. 'canopy_tempf_ps')     cdnamegem='CANF'
if(len_trim(cvar) .eq. 8  .and. cvar(1:8)  .eq. 'ustar_ps')            cdnamegem='USTR'
if(len_trim(cvar) .eq. 8  .and. cvar(1:8)  .eq. 'tstar_ps')            cdnamegem='TSTR'
if(len_trim(cvar) .eq. 8  .and. cvar(1:8)  .eq. 'rstar_ps')            cdnamegem='RSTR'
if(len_trim(cvar) .eq. 8  .and. cvar(1:8)  .eq. 'sltex_bp')            cdnamegem='SLTX'
if(len_trim(cvar) .eq. 8  .and. cvar(1:8)  .eq. 'soilq_ps')            cdnamegem='SOIQ'
if(len_trim(cvar) .eq. 12 .and. cvar(1:12) .eq. 'soil_temp_ps')        cdnamegem='SOIT'
if(len_trim(cvar) .eq. 13 .and. cvar(1:13) .eq. 'soil_moist_ps')       cdnamegem='SLMS'
if(len_trim(cvar) .eq. 17 .and. cvar(1:17) .eq. 'soil_moistfrac_ps')   cdnamegem='SLMF'
if(len_trim(cvar) .eq. 13 .and. cvar(1:13) .eq. '5050_tempc_ps')       cdnamegem='50TC'
if(len_trim(cvar) .eq. 13 .and. cvar(1:13) .eq. '5050_tempf_ps')       cdnamegem='50TF'
!SIB VARIABLES SECTION - 40 variables
if(len_trim(cvar) .eq. 10 .and. cvar(1:10) .eq. 'co2_concen')          cdnamegem='CO2C'
if(len_trim(cvar) .eq. 8  .and. cvar(1:8)  .eq. 'snow1_ps')            cdnamegem='SNO1'
if(len_trim(cvar) .eq. 8  .and. cvar(1:8)  .eq. 'snow2_ps')            cdnamegem='SNO2'
if(len_trim(cvar) .eq. 9  .and. cvar(1:9)  .eq. 'capac1_ps')           cdnamegem='CAP1'
if(len_trim(cvar) .eq. 9  .and. cvar(1:9)  .eq. 'capca2_ps')           cdnamegem='CAP2'
if(len_trim(cvar) .eq. 9  .and. cvar(1:9)  .eq. 'pco2ap_ps')           cdnamegem='PCOA'
if(len_trim(cvar) .eq. 9  .and. cvar(1:9)  .eq. 'co2flx_ps')           cdnamegem='CO2F'
if(len_trim(cvar) .eq. 9  .and. cvar(1:9)  .eq. 'sfcswa_ps')           cdnamegem='SFAL'
if(len_trim(cvar) .eq. 9  .and. cvar(1:9)  .eq. 'uplwrf_ps')           cdnamegem='SFUP'
if(len_trim(cvar) .eq. 9  .and. cvar(1:9)  .eq. 'assimn_ps')           cdnamegem='ASSM'
if(len_trim(cvar) .eq. 8  .and. cvar(1:8)  .eq. 'respg_ps')            cdnamegem='RESP'
if(len_trim(cvar) .eq. 10 .and. cvar(1:10) .eq. 'rstfac1_ps')          cdnamegem='RST1'
if(len_trim(cvar) .eq. 10 .and. cvar(1:10) .eq. 'rstfac2_ps')          cdnamegem='RST2'
if(len_trim(cvar) .eq. 10 .and. cvar(1:10) .eq. 'rstfac3_ps')          cdnamegem='RST3'
if(len_trim(cvar) .eq. 6  .and. cvar(1:6)  .eq. 'ect_ps')              cdnamegem='ECTF'
if(len_trim(cvar) .eq. 6  .and. cvar(1:6)  .eq. 'eci_ps')              cdnamegem='ECIF'
if(len_trim(cvar) .eq. 6  .and. cvar(1:6)  .eq. 'egi_ps')              cdnamegem='EGIF'
if(len_trim(cvar) .eq. 6  .and. cvar(1:6)  .eq. 'egs_ps')              cdnamegem='EGSF'
if(len_trim(cvar) .eq. 5  .and. cvar(1:5)  .eq. 'hc_ps')               cdnamegem='HCFX'
if(len_trim(cvar) .eq. 5  .and. cvar(1:5)  .eq. 'hg_ps')               cdnamegem='HGFX'
if(len_trim(cvar) .eq. 5  .and. cvar(1:5)  .eq. 'ra_ps')               cdnamegem='RAST'
if(len_trim(cvar) .eq. 5  .and. cvar(1:5)  .eq. 'rb_ps')               cdnamegem='RBST'
if(len_trim(cvar) .eq. 5  .and. cvar(1:5)  .eq. 'rc_ps')               cdnamegem='RCST'
if(len_trim(cvar) .eq. 5  .and. cvar(1:5)  .eq. 'rd_ps')               cdnamegem='RDST'
if(len_trim(cvar) .eq. 7  .and. cvar(1:7)  .eq. 'roff_ps')             cdnamegem='ROFF'
if(len_trim(cvar) .eq. 8  .and. cvar(1:8)  .eq. 'green_ps')            cdnamegem='GREN'
if(len_trim(cvar) .eq. 7  .and. cvar(1:7)  .eq. 'apar_ps')             cdnamegem='APAR'
if(len_trim(cvar) .eq. 9  .and. cvar(1:9)  .eq. 'ventmf_ps')           cdnamegem='VENT'
if(len_trim(cvar) .eq. 8  .and. cvar(1:8)  .eq. 'pco2c_ps')            cdnamegem='PCOC'
if(len_trim(cvar) .eq. 8  .and. cvar(1:8)  .eq. 'pco2i_ps')            cdnamegem='PCOI'
if(len_trim(cvar) .eq. 8  .and. cvar(1:8)  .eq. 'pco2s_ps')            cdnamegem='PCOS'
if(len_trim(cvar) .eq. 8  .and. cvar(1:8)  .eq. 'pco2m_ps')            cdnamegem='PCOM'
if(len_trim(cvar) .eq. 5  .and. cvar(1:5)  .eq. 'ea_ps')               cdnamegem='EAPR'
if(len_trim(cvar) .eq. 5  .and. cvar(1:5)  .eq. 'em_ps')               cdnamegem='EMPR'
if(len_trim(cvar) .eq. 6  .and. cvar(1:6)  .eq. 'rha_ps')              cdnamegem='RHAC'
if(len_trim(cvar) .eq. 9  .and. cvar(1:9)  .eq. 'radvbc_ps')           cdnamegem='RVDR'
if(len_trim(cvar) .eq. 9  .and. cvar(1:9)  .eq. 'radvdc_ps')           cdnamegem='RVDF'
if(len_trim(cvar) .eq. 9  .and. cvar(1:9)  .eq. 'radnbc_ps')           cdnamegem='RNDR'
if(len_trim(cvar) .eq. 9  .and. cvar(1:9)  .eq. 'radndc_ps')           cdnamegem='RNDF'
if(len_trim(cvar) .eq. 6  .and. cvar(1:6)  .eq. 'psy_ps')              cdnamegem='PSYC'
!TRACER FIELDS SECTION - 6 variables
if(len_trim(cvar) .eq. 9  .and. cvar(1:9)  .eq. 'tracer001')           cdnamegem='T001'
if(len_trim(cvar) .eq. 9  .and. cvar(1:9)  .eq. 'tracer002')           cdnamegem='T002'
if(len_trim(cvar) .eq. 9  .and. cvar(1:9)  .eq. 'tracer003')           cdnamegem='T003'
if(len_trim(cvar) .eq. 9  .and. cvar(1:9)  .eq. 'tracer004')           cdnamegem='T004'
if(len_trim(cvar) .eq. 9  .and. cvar(1:9)  .eq. 'tracer005')           cdnamegem='T005'
if(len_trim(cvar) .eq. 9  .and. cvar(1:9)  .eq. 'tracer006')           cdnamegem='T006'
!*****************************************************************************
!**** Begin extracting the data **********************************************
!*****************************************************************************
if(ivtype.ge.2) then

 if(iztrans.eq.1.or.iztrans.eq.2) coordn='HGHT' !Sigma or cartesian coordinate
 if(iztrans.eq.3) coordn='PRES' !Pressure coordinate

 ! If this is a 3-dimensional atmospheric variable      
 do k=nnb,nne,nninc

  !2D variables
  if(ivtype==2) then
    typlev=0
    coordn='NONE'
  endif
  !LEAF surface model variables
  if(ivtype==4 .or. ivtype==5) then
    typlev=k
    coordn='SLEV'
  endif
  !3D variables
  if(ivtype==3) then
    if(iztrans.eq.1.or.iztrans.eq.2) typlev=int(ztn(k,nngd))
    if(iztrans.eq.3) typlev=int(iplevs(k))
  endif
  !Patch variables
  if(ivtype==6) then
    typlev=k
    coordn='PACH'
  endif

  !Advance the number of grid output thus far
  numgrids=numgrids+1

  !Write gempak grid header information
  write(iun,*)
  write(iun,*)
  write(iun,'(a12,3i2.2,i4.4,a9)') &
    ,' Grid file: ',iyr1,imonth1,idate1,itime1,'_rams.gem'
  write(iun,'(a17)'),' GRID IDENTIFIER:'
  write(iun,'(a60)') &
    ,'    TIME1             TIME2         LEVL1 LEVL2   VCORD PARM'
  write(iun,'(3i2.2,a1,2i2.2,24X,i6,10X,a4,1X,a5)') &
    ,iyr2,imonth2,idate2,'/',ihour2,imin2,typlev,coordn,cdnamegem
  write(iun,'(a6,f7.2,a1,f8.2,a1,f7.2,a1,f8.2,18X,a13,2i5)') &
    ,' AREA:',swlat1,';',swlon1,';',nelat1,';',nelon1 &
    ,'GRID SIZE: ',xpts,ypts
  write(iun,'(a17,i3,a18,i3)') &
    ,' COLUMNS:     1  ',xpts,'     ROWS:     1  ',ypts
  write(iun,*)
  write(iun,'(a21)') ' Scale factor: 10** 0'
  write(iun,*)
  write(iun,*)
  write(iun,'(a8,i7,2X,i7,2X,i7,2X,i7,2X,i7,2X,i7,2X,i7,2X,i7)') &
    ' COLUMN:',(i,i=1,8)
  write(iun,'(8X,i7,2X,i7,2X,i7,2X,i7,2X,i7,2X,i7,2X,i7,2X,i7)') &
    (i,i=9,xpts)

  !Determine how to write out arrays depending on vtype/dimension
  !and if simulation is 2D or 3D
  if(ivtype.eq.2) kk=1
  if(ivtype.ge.3) kk=k
  jj=1
  jcount=yend/yinc
  do j=yend,ybeg,-yinc
    if(nje>1)jj=j
    write(iun,'(a4,i3,a4,8e14.6)') &
        ' ROW',jcount,'    ',(a(i,jj,kk),i=xbeg,7*xinc+1,xinc)
    write(iun,'(11X,8e14.6)') (a(i,jj,kk),i=7*xinc+1+xinc,xend,xinc)
    jcount=jcount-1
  enddo

 enddo !loop k-levels

 write(iuntag,'(a9,i5,a1)'),'NUMGRDS="',numgrids,'"'
 backspace iuntag

else
   print*,'unknown ivtype',ivtype
   stop 'RAMS_text'
endif

return
END SUBROUTINE rams_text
