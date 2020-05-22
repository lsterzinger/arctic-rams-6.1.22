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

Subroutine micro_master ()

use micphys
use rconstants

implicit none

integer :: lhcat,khcat,lcat,nd1,nd2,nip,ilcat,ilhcat,idum
integer, dimension(8) :: lcat0

real, dimension(7,16) :: dstprms,dstprms1,dstprms2
real, dimension(16,16) :: jpairr,jpairc
character(len=strl1) :: dataline,cname

data lcat0 /1,2,3,4,5,6,7,16/ ! lcat corressponding to lhcat

data dstprms1/ &
!----------------------------------------------------------------------
! shape      cfmas   pwmas      cfvt    pwvt     dmb0      dmb1
!---------------------------------------------------------------------- 
    .5,      524.,     3.,    3173.e4,     2.,   2.e-6,   50.e-6,  & !cloud
    .5,      524.,     3.,     149.,     .5,   .1e-3,    5.e-3,  & !rain
  .179,     110.8,   2.91,  5.769e5,   1.88,  15.e-6,  125.e-6,  & !pris col
  .179,  2.739e-3,   1.74,  188.146,   .933,   .1e-3,   10.e-3,  & !snow col
    .5,      .496,    2.4,    3.084,     .2,   .1e-3,   10.e-3,  & !aggreg
    .5,      157.,     3.,     93.3,     .5,   .1e-3,    5.e-3,  & !graup
    .5,      471.,     3.,     161.,     .5,   .1e-3,   10.e-3,  & !hail 
  .429,     .8854,    2.5,     316.,   1.01,      00,       00,  & !pris hex
 .3183,   .377e-2,     2.,     316.,   1.01,      00,       00,  & !pris den
 .1803,   1.23e-3,    1.8,  5.769e5,   1.88,      00,       00,  & !pris ndl
    .5,     .1001,  2.256,   3.19e4,   1.66,      00,       00,  & !pris ros
  .429,     .8854,    2.5,    4.836,    .25,      00,       00,  & !snow hex
 .3183,   .377e-2,     2.,    4.836,    .25,      00,       00,  & !snow den
 .1803,   1.23e-3,    1.8,  188.146,   .933,      00,       00,  & !snow ndl
    .5,     .1001,  2.256,  1348.38,  1.241,      00,       00,  & !snow ros
    .5,      524.,     3.,    3173.e4,     2.,  65.e-6,  100.e-6/    !drizzle

data dstprms2/ &
!----------------------------------------------------------------------
! shape      cfmas   pwmas      cfvt    pwvt     dmb0      dmb1
!----------------------------------------------------------------------
    .5,      524.,     3.,    3173.e4,     2.,   2.e-6,   50.e-6,  & !cloud
    .5,      524.,     3.,     144.,   .497,   .1e-3,    5.e-3,  & !rain
  .179,     110.8,   2.91,    1538.,   1.00,  15.e-6,  125.e-6,  & !pris col
  .179,  2.739e-3,   1.74,     27.7,   .484,   .1e-3,   10.e-3,  & !snow col
    .5,      .496,    2.4,     16.1,   .416,   .1e-3,   10.e-3,  & !aggreg
    .5,      157.,     3.,     332.,   .786,   .1e-3,    5.e-3,  & !graup
    .5,      471.,     3.,    152.1,   .497,   .8e-3,   10.e-3,  & !hail
  .429,     .8854,    2.5,   20801.,  1.377,      00,       00,  & !pris hex
 .3183,   .377e-2,     2.,     56.4,   .695,      00,       00,  & !pris den
 .1803,   1.23e-3,    1.8,   1617.9,   .983,      00,       00,  & !pris ndl
    .5,     .1001,  2.256,    6239.,   1.24,      00,       00,  & !pris ros
  .429,     .8854,    2.5,    30.08,   .563,      00,       00,  & !snow hex
 .3183,   .377e-2,     2.,     3.39,   .302,      00,       00,  & !snow den
 .1803,   1.23e-3,    1.8,     44.6,   .522,      00,       00,  & !snow ndl
    .5,     .1001,  2.256,    125.7,   .716,      00,       00,  & !snow ros
    .5,      524.,     3.,   1.26e7,   1.91,  65.e-6,  100.e-6/    !drizzle

data jpairr/  &
     0,  0,  0,  1,  2,  3,  4,  0,  0,  0,  0,  5,  6,  7,  8,  0,  &
     0,  0,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21,  0,  &
     0, 22, 23, 24,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
    25, 26, 27, 28,  0,  0,  0, 29, 30, 31, 32,  0,  0,  0,  0, 33,  &
    34, 35, 36, 37,  0,  0,  0, 38, 39, 40, 41, 42, 43, 44, 45, 46,  &
    47, 48, 49, 50, 51,  0,  0, 52, 53, 54, 55, 56, 57, 58, 59, 60,  &
    61, 62, 63, 64, 65, 66,  0, 67, 68, 69, 70, 71, 72, 73, 74, 75,  &
     0, 76,  0, 77,  0,  0,  0, 78,  0,  0,  0, 79, 80, 81, 82,  0,  &
     0, 83,  0, 84,  0,  0,  0,  0, 85,  0,  0, 86, 87, 88, 89,  0,  &
     0, 90,  0, 91,  0,  0,  0,  0,  0, 92,  0, 93, 94, 95, 96,  0,  &
     0, 97,  0, 98,  0,  0,  0,  0,  0,  0, 99,100,101,102,103,  0,  &
   104,105,106,  0,  0,  0,  0,107,108,109,110,111,  0,  0,  0,112,  &
   113,114,115,  0,  0,  0,  0,116,117,118,119,  0,120,  0,  0,121,  &
   122,123,124,  0,  0,  0,  0,125,126,127,128,  0,  0,129,  0,130,  &
   131,132,133,  0,  0,  0,  0,134,135,136,137,  0,  0,  0,138,139,  &
     0,  0,  0,140,141,142,143,  0,  0,  0,  0,144,145,146,147,  0/

data jpairc/  &
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
     0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
     0,  2,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
     4,  5,  6,  7,  0,  0,  0,  8,  9, 10, 11,  0,  0,  0,  0, 12,  &
    13, 14, 15, 16, 17,  0,  0, 18, 19, 20, 21, 22, 23, 24, 25, 26,  &
    27, 28, 29, 30, 31, 32,  0, 33, 34, 35, 36, 37, 38, 39, 40, 41,  &
    42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57,  &
     0, 58,  0,  0,  0,  0,  0, 59,  0,  0,  0,  0,  0,  0,  0,  0,  &
     0, 60,  0,  0,  0,  0,  0,  0, 61,  0,  0,  0,  0,  0,  0,  0,  &
     0, 62,  0,  0,  0,  0,  0,  0,  0, 63,  0,  0,  0,  0,  0,  0,  &
     0, 64,  0,  0,  0,  0,  0,  0,  0,  0, 65,  0,  0,  0,  0,  0,  &
    66, 67, 68,  0,  0,  0,  0, 69, 70, 71, 72, 73,  0,  0,  0, 74,  &
    75, 76, 77,  0,  0,  0,  0, 78, 79, 80, 81,  0, 82,  0,  0, 83,  &
    84, 85, 86,  0,  0,  0,  0, 87, 88, 89, 90,  0,  0, 91,  0, 92,  &
    93, 94, 95,  0,  0,  0,  0, 96, 97, 98, 99,  0,  0,  0,100,101,  &
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0/

!  Define several parameters from above data list

do lhcat=1,nhcat
   !Using original RAMS 4.3 power laws
   if(iplaws==0) then
     dstprms(1,lhcat) = dstprms1(1,lhcat)
     dstprms(2,lhcat) = dstprms1(2,lhcat)
     dstprms(3,lhcat) = dstprms1(3,lhcat)
     dstprms(4,lhcat) = dstprms1(4,lhcat)
     dstprms(5,lhcat) = dstprms1(5,lhcat)
     dstprms(6,lhcat) = dstprms1(6,lhcat)
     dstprms(7,lhcat) = dstprms1(7,lhcat)
   !Using Carver/Mitchell 1996 power laws
   elseif(iplaws==1.or.iplaws==2) then
     dstprms(1,lhcat) = dstprms2(1,lhcat)
     dstprms(2,lhcat) = dstprms2(2,lhcat)
     dstprms(3,lhcat) = dstprms2(3,lhcat)
     dstprms(4,lhcat) = dstprms2(4,lhcat)
     dstprms(5,lhcat) = dstprms2(5,lhcat)
     dstprms(6,lhcat) = dstprms2(6,lhcat)
     dstprms(7,lhcat) = dstprms2(7,lhcat)
   endif

   shape(lhcat) = dstprms(1,lhcat)
   cfmas(lhcat) = dstprms(2,lhcat)
   pwmas(lhcat) = dstprms(3,lhcat)
   cfvt (lhcat) = dstprms(4,lhcat)
   pwvt (lhcat) = dstprms(5,lhcat)

   do khcat=1,nhcat
      ipairc(lhcat,khcat) = jpairc(lhcat,khcat)
      ipairr(lhcat,khcat) = jpairr(lhcat,khcat)
   enddo
enddo

do lcat=1,ncat
   lhcat = lcat0(lcat)
   emb0 (lcat) = cfmas(lhcat) * dstprms(6,lhcat) ** pwmas(lhcat)
   emb1 (lcat) = cfmas(lhcat) * dstprms(7,lhcat) ** pwmas(lhcat)
enddo

if (level .ne. 3) return

if(mkcoltab.lt.0.or.mkcoltab.gt.1)then
   print*, 'mkcoltab set to ',mkcoltab, 'which is out of bounds'
   stop 'mkcoltab'
endif

cname=coltabfn(1:len_trim(coltabfn))

if(mkcoltab.eq.1)then

! Make collection table and write to file

   CALL mkcoltb ()
   open(91,file=cname,form='formatted',status='unknown')
   rewind(91)
   write(91,181)
   do lcat = 1,ncat
      write(91,182)lcat,gnu(lcat),emb0(lcat),emb1(lcat)
   enddo
   write(91,180)
   write(91,183)
   do lhcat = 1,nhcat
      write(91,182)lhcat,cfmas(lhcat),pwmas(lhcat)  &
         ,cfvt(lhcat),pwvt(lhcat)
   enddo
   write(91,180)
   do nip=1,npairc
      write(91,186)nip
      write(91,184)(nd2,(coltabc(nd1,nd2,nip)  &
         ,nd1=1,nembc),nd2=1,nembc)
   enddo
   write(91,180)
   do nip=1,npairr
      write(91,187)nip
      write(91,184)(nd2,(coltabr(nd1,nd2,nip)  &
         ,nd1=1,nembc),nd2=1,nembc)
   enddo
   close(91)

endif

!  Read collection table regardless of whether or not your
!  just created it so that we keep the same variable precision
!  between initial computation and precision stored in collection
!  table.

   open(91,file=cname,form='formatted',status='old')
   read(91,185)dataline
   do ilcat = 1,ncat
      read(91,182)lcat,gnu(lcat),emb0(lcat),emb1(lcat)
   enddo
   read(91,185)dataline
   read(91,185)dataline
   do ilhcat = 1,nhcat
      read(91,182)lhcat,cfmas(lhcat),pwmas(lhcat)  &
         ,cfvt(lhcat),pwvt(lhcat)
   enddo
   read(91,185)dataline
   do nip=1,npairc
      read(91,185)dataline
      read(91,184)(idum,(coltabc(nd1,nd2,nip)  &
         ,nd1=1,nembc),nd2=1,nembc)
   enddo
   read(91,185)dataline
   do nip=1,npairr
      read(91,185)dataline
      read(91,184)(idum,(coltabr(nd1,nd2,nip)  &
         ,nd1=1,nembc),nd2=1,nembc)
   enddo
   close(91)

180  format(' ')
181  format(' lcat    gnu        emb0       emb1    ')
182  format(i4,7e11.4)
183  format(' lhcat  cfmas      pwmas       cfvt       pwvt')
184  format(i3,20f6.2)
185  format(a80)
186  format('ipairc',i4)
187  format('ipairr',i4)

return
END SUBROUTINE micro_master

!##############################################################################
Subroutine initqin (n1,n2,n3,q2,q6,q7,pi0,pp,theta)

use micphys
use rconstants

implicit none

integer :: n1,n2,n3,i,j,k
real, dimension(n1,n2,n3) :: q2,q6,q7,pi0,pp,theta

! Initialize Q2, Q6, Q7

do j = 1,n3
   do i = 1,n2
      do k = 1,n1
         pitot(k) = pi0(k,i,j) + pp(k,i,j)
         tair(k) = theta(k,i,j) * pitot(k) / cp

         if(irain .ge. 1) q2(k,i,j) = tair(k) - 193.16
         if(igraup .ge. 1) q6(k,i,j) = 0.5 * min(0.,tair(k) - 273.15)
         if(ihail .ge. 1) q7(k,i,j) = 0.5 * min(0.,tair(k) - 273.15)

      enddo
   enddo
enddo

return
END SUBROUTINE initqin

!##############################################################################
Subroutine init_ifn (n1,n2,n3,cifnp,dn0,ifm)

use micphys
use rconstants
use mem_grid

implicit none

integer :: n1,n2,n3,i,j,k,ifm
real, dimension(n1,n2,n3) :: cifnp,dn0
real :: cin_maxt

! Initialize IFN
if(iaeroprnt==1 .and. print_msg) then
 print*,' '
 print*,'Start Initializing Ice Nuclei concentration'
endif

do j = 1,n3
 do i = 1,n2
  do k = 1,n1

   !Convert RAMSIN #/mg to #/kg
   cin_maxt = cin_max * 1.e6

   !Set up Vertical profile 
   cifnp(k,i,j)=cin_maxt
   !if(k<=2) cifnp(k,i,j)=cin_maxt
   ! Exponential decrease that scales with pressure decrease
   !if(k>2)  cifnp(k,i,j)=cin_maxt*exp(-zt(k)/7000.)

   !Output initial sample profile
   if(iaeroprnt==1 .and. print_msg .and. i==1 .and. j==1) then
     if(k==1) print*,' Ice Nuclei - init (k,zt,ifn/kg,ifn/L) on Grid:',ifm
     print'(a9,i5,f11.1,2f17.7)',' IFN-init' &
        ,k,zt(k),cifnp(k,i,j),cifnp(k,i,j)/1.e3*dn0(k,i,j)
   endif


  enddo
 enddo
enddo

return
END SUBROUTINE init_ifn

!##############################################################################
Subroutine init_ccn (n1,n2,n3,cccnp,cccmp,dn0,ifm)

use micphys
use rconstants
use mem_grid

implicit none

integer :: n1,n2,n3,i,j,k,ifm
real, dimension(n1,n2,n3) :: cccnp,cccmp,dn0
real :: ccn_maxt

!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! Make sure that the profile here is consistent with the reset_ccn profile
!!!!!!!!!!!!!!!!!!!!!!!!!!

! Initialize CCN
if(iaeroprnt==1 .and. print_msg) then
 print*,' '
 print*,'Start Initializing CCN concentration'
endif

!Convert RAMSIN #/mg to #/kg
 ccn_maxt = ccn_max * 1.e6 

do j = 1,n3
 do i = 1,n2
  do k = 1,n1

   !Set up Vertical profile
   !if(k<=2) cccnp(k,i,j)=ccn_maxt
   cccnp(k,i,j)=ccn_maxt
   !if(zt(k)<=850) cccnp(k,i,j)=ccn_maxt
   !Exponential decrease that scales with pressure decrease
   !if(zt(k)>850) cccnp(k,i,j)=0.6*(zt(k)-850)+ccn_maxt ! 0.3 & 0.6
   !if(k>2)  cccnp(k,i,j)=ccn_maxt*exp(-zt(k)/7000.)

   !Output initial sample profile
   if(iaeroprnt==1 .and. print_msg .and. i==1 .and. j==1) then
     if(k==1) print*,' CCN-init (k,zt,ccn/mg,ccn/cc) on Grid:',ifm
     print'(a9,i5,f11.1,2f17.7)',' CCN-init' &
        ,k,zt(k),cccnp(k,i,j)/1.e6,cccnp(k,i,j)/1.e6*dn0(k,i,j)
   endif

   !Set up Field of CCN mass mixing ratio (kg/kg)
   cccmp(k,i,j) = ((aero_medrad(1)*aero_rg2rm(1))**3.) &
                *cccnp(k,i,j)/(0.23873/aero_rhosol(1))

  enddo
 enddo
enddo

return
END SUBROUTINE init_ccn

!##############################################################################
Subroutine reset_ccn (m1,time,rv,k1,k2)

use micphys
use rconstants
use mem_grid, only:zt

implicit none

integer :: k,m1,k1,k2
real :: ccn_maxt,time
real, dimension(m1) :: rv

!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! Make sure that the profile here is consistent with the init_ccn profile
!!!!!!!!!!!!!!!!!!!!!!!!!!

!Convert RAMSIN #/mg to #/kg
ccn_maxt = ccn_max * 1.e6 

!decay timescale
!fccnts = 28800. !8 hours
if (k1==0)k1=199
if (k2==0)k2=1

!print*, 'k1',k1,zt(k1)
do k = 1,m1

 if (rv(k)/rvlsair(k).lt.0.99) then !If subsaturated
  !Definitely force the aerosol profile where it is subsaturated
  !Allow decrease in forced value depending on value of iforceccn

  !1) decrease bl and zt(k)<blh or 
  !2) decrease ft and zt(k)>blh
  !3) decrease everywhere or
  !4) Hold reset everywhere
  select case (iforceccn)
  case(1) !decrease if below BLH/k1, reset elsewhere
    if(zt(k).lt.blh .or. zt(k).le.zt(k1)) then
      aerocon(k,1)=ccn_maxt * exp(-1.* (time-fccnstart)/fccnts)
    else
      aerocon(k,1)=ccn_maxt 
    endif
  
  case(2) !decrease if above BLH/k2, reset in bottom layer (k=2 is bottom layer of atmosphere)
    if (zt(k).ge.blh .or. zt(k).ge.zt(k2))  then
      aerocon(k,1)=ccn_maxt * exp(-1.* (time-fccnstart)/fccnts)
    else if (k.le.2) then
      aerocon(k,1)=ccn_maxt
    endif

  case(3) !decrease everywhere
    aerocon(k,1)=ccn_maxt * exp(-1.* (time-fccnstart)/fccnts)

  case(4) !hold constant above BLH/k2, and k<=2
    if (zt(k).ge.blh .or. zt(k).ge.zt(k2) .or. k.le.2) then
      aerocon(k,1)=ccn_maxt
    endif

  end select

  ! if (fccnts .gt. 0 .and. (iforceccn.eq.3 .or. &
  !     (iforceccn.eq.1 .and. (zt(k).lt.blh .or. zt(k).le.zt(k1))) .or. &
  !     (iforceccn.eq.2 .and. (zt(k).ge.blh .or. zt(k).ge.zt(k2))))) then
  !    aerocon(k,1)=ccn_maxt * exp(-1.* (time-fccnstart)/fccnts)
     
  ! ! If resetting in BL (IFORCECCN=2), only reset in bottom 10 layers of BL.
  ! else if (((iforceccn.eq.2).and.(k.le.2)).or.(iforceccn.ne.2)) then!Force the aerosol conc to be the same as at start
  !  ! else !Force the aerosol conc to be the same as at start
  !    aerocon(k,1)=ccn_maxt 
  ! endif

   !Set up Field of CCN mass mixing ratio (kg/kg)
   aeromas(k,1) = ((aero_medrad(1)*aero_rg2rm(1))**3.) &
                *aerocon(k,1)/(0.23873/aero_rhosol(1))
 endif
enddo

return
END SUBROUTINE reset_ccn

!##############################################################################
Subroutine init_gccn (n1,n2,n3,gccnp,gccmp,dn0,ifm)

use micphys
use rconstants
use mem_grid

implicit none

integer :: n1,n2,n3,i,j,k,ifm
real, dimension(n1,n2,n3) :: gccnp,gccmp,dn0
real :: gccn_maxt

! Initialize Giant-CCN
if(iaeroprnt==1 .and. print_msg) then
 print*,' '
 print*,'Start Initializing GCCN concentration'
endif

!Convert RAMSIN #/mg to #/kg
 gccn_maxt = gccn_max * 1.e6 

do j = 1,n3
 do i = 1,n2
  do k = 1,n1

   !Set up Vertical profile
   if(k<=2) gccnp(k,i,j)=gccn_maxt
   ! Exponential decrease that scales with pressure decrease
   if(k>2)  gccnp(k,i,j)=gccn_maxt*exp(-zt(k)/7000.)

   !Output initial sample profile
   if(iaeroprnt==1 .and. print_msg .and. i==1 .and. j==1) then
     if(k==1) print*,' GCCN-init (k,zt,gccn/mg,gccn/cc) on Grid:',ifm
     print'(a10,i5,f11.1,2f17.7)',' GCCN-init' &
        ,k,zt(k),gccnp(k,i,j)/1.e6,gccnp(k,i,j)/1.e6*dn0(k,i,j)
   endif

   !Set up Field of GCCN mass mixing ratio (kg/kg)
   gccmp(k,i,j) = ((aero_medrad(2)*aero_rg2rm(2))**3.) &
                *gccnp(k,i,j)/(0.23873/aero_rhosol(2))

  enddo
 enddo
enddo

return
END SUBROUTINE init_gccn

!##############################################################################
Subroutine init_dust (n1,n2,n3,md1np,md2np,md1mp,md2mp,ifm)

use micphys
use rconstants
use mem_grid

implicit none

integer :: n1,n2,n3,i,j,k,ifm
real, dimension(n1,n2,n3) :: md1np,md2np,md1mp,md2mp
real :: dust1_maxt,dust2_maxt

! Initialize Dust
if(iaeroprnt==1 .and. print_msg) then
 print*,' '
 print*,'Start Initializing DUST concentration'
 print*,'idust,iaerorad',idust,iaerorad
endif

!Convert RAMSIN #/mg to #/kg
 dust1_maxt = dust1_max * 1.e6
 dust2_maxt = dust2_max * 1.e6  

do j = 1,n3
 do i = 1,n2
  do k = 1,n1

   !If not using dust source model
   if(idust == 1) then
     !Set up concentration of SMALL MODE Mineral Dust (#/kg)
     if(k<=2) md1np(k,i,j)=dust1_maxt
     if(k>2)  md1np(k,i,j)=dust1_maxt*exp(-zt(k)/7000.)
     !Set up concentration of LARGE MODE Mineral Dust (#/kg)
     if(k<=2) md2np(k,i,j)=dust2_maxt
     if(k>2)  md2np(k,i,j)=dust2_maxt*exp(-zt(k)/7000.)

     !Set up Field of SMALL MODE DUST mass (kg/kg)
     md1mp(k,i,j) = ((aero_medrad(3)*aero_rg2rm(3))**3.) &
                    *md1np(k,i,j)/(0.23873/aero_rhosol(3))
     !Set up Field of LARGE MODE DUST mass (kg/kg)
     md2mp(k,i,j) = ((aero_medrad(4)*aero_rg2rm(4))**3.) &
                    *md2np(k,i,j)/(0.23873/aero_rhosol(4))

   !If using dust source model, do not initialize with background dust
   elseif(idust == 2) then
     md1np(k,i,j) = 0.
     md2np(k,i,j) = 0.
     md1mp(k,i,j) = 0.
     md2mp(k,i,j) = 0.
   endif

   !Output sample initial profile
   if(iaeroprnt==1 .and. print_msg .and. i==1 .and. j==1) then
     if(k==1) print*,' Dust-init (k,zt,dust1/mg,dust2/mg) on Grid:',ifm
     print'(a10,i5,f11.1,2f17.7)',' DUST-init' &
        ,k,zt(k),md1np(k,i,j)/1.e6,md2np(k,i,j)/1.e6
   endif

  enddo
 enddo
enddo

return
END SUBROUTINE init_dust

!##############################################################################
Subroutine init_salt (n1,n2,n3,salt_film_np,salt_jet_np,salt_spum_np &
                       ,salt_film_mp,salt_jet_mp,salt_spum_mp,ifm)

use micphys
use rconstants
use mem_grid

implicit none

integer :: n1,n2,n3,i,j,k,ifm
real, dimension(n1,n2,n3) :: salt_film_np,salt_jet_np,salt_spum_np
real, dimension(n1,n2,n3) :: salt_film_mp,salt_jet_mp,salt_spum_mp
real :: saltf_maxt,saltj_maxt,salts_maxt

! Initialize Sea-salt
if(iaeroprnt==1 .and. print_msg) then
 print*,' '
 print*,'Start Initializing SALT concentration'
 print*,'isalt,iaerorad',isalt,iaerorad
endif

!Convert RAMSIN #/mg to #/kg
 saltf_maxt = saltf_max * 1.e6
 saltj_maxt = saltj_max * 1.e6
 salts_maxt = salts_max * 1.e6

do j = 1,n3
 do i = 1,n2
  do k = 1,n1

   !If not using dust source model
   if(isalt == 1) then
     !Set up concentration of FILM MODE SALT (#/kg)
     if(k<=2) salt_film_np(k,i,j)=saltf_maxt
     if(k>2)  salt_film_np(k,i,j)=saltf_maxt*exp(-zt(k)/7000.)
     !Set up concentration of JET MODE Mineral Dust (#/kg)
     if(k<=2) salt_jet_np(k,i,j)=saltj_maxt
     if(k>2)  salt_jet_np(k,i,j)=saltj_maxt*exp(-zt(k)/7000.)
     !Set up concentration of SPUME MODE Mineral Dust (#/kg)
     if(k<=2) salt_spum_np(k,i,j)=salts_maxt
     if(k>2)  salt_spum_np(k,i,j)=salts_maxt*exp(-zt(k)/7000.)

     !Set up 3D Field of FILM MODE SALT mass (kg/kg)
     salt_film_mp(k,i,j) = ((aero_medrad(5)*aero_rg2rm(5))**3.) &
                           *salt_film_np(k,i,j)/(0.23873/aero_rhosol(5))
     !Set up 3D Field of JET MODE SALT mass (kg/kg)
     salt_jet_mp(k,i,j)  = ((aero_medrad(6)*aero_rg2rm(6))**3.) &
                           *salt_jet_np(k,i,j) /(0.23873/aero_rhosol(6))
     !Set up 3D Field of SPUME MODE SALT mass (kg/kg)
     salt_spum_mp(k,i,j) = ((aero_medrad(7)*aero_rg2rm(7))**3.) &
                           *salt_spum_np(k,i,j)/(0.23873/aero_rhosol(7))

   !If using salt source model, do not initialize with background salt
   elseif(isalt == 2) then
     salt_film_np(k,i,j) = 0.
     salt_jet_np(k,i,j)  = 0.
     salt_spum_np(k,i,j) = 0.
     salt_film_mp(k,i,j) = 0.
     salt_jet_mp(k,i,j)  = 0.
     salt_spum_mp(k,i,j) = 0.
   endif

   !Output initial sample profile
   if(iaeroprnt==1 .and. print_msg .and. i==1 .and. j==1) then
     if(k==1) print*,' Salt-init (k,zt,film/mg,jet/mg,spume/mg) on Grid:',ifm
     print'(a10,i5,f11.1,3f17.7)',' SALT-init',k,zt(k),salt_film_np(k,i,j)/1.e6 &
      ,salt_jet_np(k,i,j)/1.e6,salt_spum_np(k,i,j)/1.e6
   endif

  enddo
 enddo
enddo

return
END SUBROUTINE init_salt

!##############################################################################
Subroutine init_tracer (n1,n2,n3,tracerp,dn0,ifm,nsc)

use micphys
use rconstants
use mem_grid
use mem_tracer
use node_mod

implicit none

integer :: n1,n2,n3,i,j,k,ifm,nsc,ii,jj
real, dimension(n1,n2,n3) :: tracerp,dn0
real :: ccn_maxt

! Initialize Tracers
if(print_msg)then
 print*,' '
 print*,'Start Initializing Tracers, Grid:',ifm,' Tracer:',nsc
endif

!Convert RAMSIN #/mg to #/kg
 ccn_maxt = ccn_max * 1.e6 

do j = 1,n3
 do i = 1,n2
  do k = 1,n1

   !Get absolute grid points for parallel (& sequential) computation
   ii = i+mi0(ngrid)
   jj = j+mj0(ngrid)

   !Set up Vertical profile, Exponential decrease that scales with pressure
   if(nsc==1) then
    if(k<=2) tracerp(k,i,j)=ccn_maxt
    if(k>2)  tracerp(k,i,j)=ccn_maxt*exp(-zt(k)/7000.)
   endif
   !Set up Field of CCN mass mixing ratio (kg/kg)
   if(nsc==2) then
    tracerp(k,i,j) = ((aero_medrad(1)*aero_rg2rm(1))**3.) &
                *tracer_g(1,ifm)%tracerp(k,i,j)/(0.23873/aero_rhosol(1))
   endif

  enddo
 enddo
enddo

return
END SUBROUTINE init_tracer

!##############################################################################
Subroutine jnmbinit ()

use micphys

implicit none

if (level /= 3) then

   if (level <= 1) then
      jnmb(1) = 0
   else
      jnmb(1) = 4
   endif

   jnmb(2) = 0
   jnmb(3) = 0
   jnmb(4) = 0
   jnmb(5) = 0
   jnmb(6) = 0
   jnmb(7) = 0
   jnmb(8) = 0

else

   jnmb(1) = icloud
   jnmb(2) = irain
   jnmb(3) = ipris
   jnmb(4) = isnow
   jnmb(5) = iaggr
   jnmb(6) = igraup
   jnmb(7) = ihail
   jnmb(8) = idriz

   if (icloud .eq. 1) jnmb(1) = 4
   if (irain  .eq. 1) jnmb(2) = 2
   if (ipris  .ge. 1) jnmb(3) = 5
   if (isnow  .eq. 1) jnmb(4) = 2
   if (iaggr  .eq. 1) jnmb(5) = 2
   if (igraup .eq. 1) jnmb(6) = 2
   if (ihail  .eq. 1) jnmb(7) = 2
   if (idriz  .eq. 1) jnmb(8) = 4

   if (irain == 5 .or. isnow == 5 .or. iaggr == 5 .or.  &
      igraup == 5 .or. ihail == 5) then

      if (irain  >= 1) jnmb(2) = 5
      if (isnow  >= 1) jnmb(4) = 5
      if (iaggr  >= 1) jnmb(5) = 5
      if (igraup >= 1) jnmb(6) = 5
      if (ihail  >= 1) jnmb(7) = 5

   endif

endif

return
END SUBROUTINE jnmbinit

!##############################################################################
Subroutine micinit ()

use micphys

implicit none

integer :: lhcat,lcat,ia
integer, dimension(16) :: lcat0
real :: cfmasi,c1,glg,glg1,glg2,glgm,glgc,flngi,dpsi,embsip,dnsip
real, external :: gammln
real, external :: gammp
real, external :: gammq

data lcat0 /1,2,3,4,5,6,7,3,3,3,3,4,4,4,4,8/ ! lcat corressponding to lhcat

! Initialize arrays based on microphysics namelist parameters
! Note: pparm for pristine ice is obsolete since IPRIS = 0 or 5 only

parm(1) = cparm
parm(2) = rparm
parm(4) = sparm
parm(5) = aparm
parm(6) = gparm
parm(7) = hparm
parm(8) = dparm

if (icloud .le. 1) parm(1) = .3e9
if (irain  .eq. 1) parm(2) = .1e-2
if (isnow  .eq. 1) parm(4) = .1e-2
if (iaggr  .eq. 1) parm(5) = .1e-2
if (igraup .eq. 1) parm(6) = .1e-2
if (ihail  .eq. 1) parm(7) = .3e-2
if (idriz  .eq. 1) parm(8) = .1e6  !# per kg ~ m^3 
                                   !(mid-range avg from Feingold(99)

dps = 125.e-6
dps2 = dps ** 2
rictmin = 1.0001
rictmax = 0.9999 * float(nembc)

do lhcat = 1,nhcat
   lcat=lcat0(lhcat)

   cfden(lhcat) = cfmas(lhcat) * 6.0 / 3.14159
   pwden(lhcat) = pwmas(lhcat) - 3.
   emb0log(lcat) = log(emb0(lcat))
   emb1log(lcat) = log(emb1(lcat))

! Define coefficients [frefac1, frefac2] used for terminal velocity
! and Reynolds number

   cfmasi = 1. / cfmas(lhcat)
   pwmasi(lhcat) = 1. / pwmas(lhcat)
   pwen0(lhcat) = 1. / (pwmas(lhcat) + 1.)
   pwemb0(lhcat) = pwmas(lhcat) / (pwmas(lhcat) + 1.)
   c1 = 1.5 + .5 * pwvt(lhcat)

   glg   = gammln(gnu(lcat))
   glg1  = gammln(gnu(lcat) + 1.)
   glg2  = gammln(gnu(lcat) + 2.)
   glgm  = gammln(gnu(lcat) + pwmas(lhcat))
   glgc  = gammln(gnu(lcat) + c1)

   if (jnmb(lcat) .eq. 3) then
      cfemb0(lhcat) = cfmas(lhcat) * exp(glgm - glg)  &
         ** pwen0(lhcat) * (1. / parm(lcat)) ** pwemb0(lhcat)
      cfen0(lhcat) = parm(lcat) * (exp(glg - glgm) / parm(lcat))  &
         ** pwen0(lhcat)
   endif

   dnfac(lhcat) = (cfmasi * exp(glg - glgm)) ** pwmasi(lhcat)

   frefac1(lhcat) = shape(lhcat) * exp(glg1 - glg)  &
      * (cfmasi * exp(glg - glgm)) ** pwmasi(lhcat)

   frefac2(lhcat) = shape(lhcat) * 0.229 * sqrt(cfvt(lcat))  &
      * (cfmasi * exp(glg - glgm)) ** (pwmasi(lhcat) * c1)  &
      * exp(glgc - glg)

   sipfac(lhcat) = .785 * exp(glg2 - glg)  &
      * (cfmasi * exp(glg - glgm)) ** (2. * pwmasi(lhcat))

   dict(lcat) = float(nembc-1) / (emb1log(lcat) - emb0log(lcat))

   dpsmi(lhcat) = 1. / (cfmas(lhcat) * dps ** pwmas(lhcat))
   if (lhcat .le. 4) gamm(lhcat) = exp(glg)
   if (lhcat .le. 4) gamn1(lhcat) = exp(glg1)

! gam1   :  the integral of the pristine distribution from dps to infty
! gam2   :  the integral of the snow dist. from 0 to dps
! gam3   :  values of the exponential exp(-dps/dn)

enddo

!***********************************************************************
!************** Secondary Ice Production Arrays ************************
!***********************************************************************
flngi = 1. / float(ngam)
do ia=1,ngam
   dpsi = dps * 1.e6 / float(ia)

   gam(ia,1) = gammq(gnu(3) + 1., dpsi)
   gam(ia,2) = gammp(gnu(4) + 1., dpsi)
   gam(ia,3) = exp(-dpsi)

   GAMINC(IA,1)=gammq(GNU(3),dpsi)
   GAMINC(IA,2)=gammp(GNU(4),dpsi)

   embsip = emb1(1) * float(ia) * flngi
   dnsip = dnfac(1) * embsip ** pwmasi(1)
   gamsip13(1,ia)=gammp(gnu(1),13.e-6/dnsip)
   gamsip24(1,ia)=gammq(gnu(1),24.e-6/dnsip)

   embsip = emb1(8) * float(ia) * flngi
   dnsip = dnfac(16) * embsip ** pwmasi(16)
   gamsip13(2,ia)=gammp(gnu(8),13.e-6/dnsip)
   gamsip24(2,ia)=gammq(gnu(8),24.e-6/dnsip)
enddo

return
END SUBROUTINE micinit

!##############################################################################
Subroutine setupF94 ()

!Routine to set up rain size bin partitions for F94 collection
!routine for rain-pris ice collection forming hail

use micphys

implicit none

integer :: nbins,i
real :: pi314,mind,maxd,m,ii3,rmin,mmin,xjo,massdobl,tmp1

!--local constants
nbins=94
pi314=4.*ATAN(1.)
ii3=(1./3.)
massdobl=2. !mass doubling bin factor for rain

!--set up rain diameter bin partitions
mind=3.4668048E-07            !min rain diam [m] from subr 'mkcoltab'
maxd=1.7334025E-02            !max rain diam [m] from subr 'mkcoltab'
rmin=mind/2.                  !minimum radius
mmin=(4./3.)*pi314*(rmin)**3  !minimum mass
tmp1=massdobl*3.              !mass is doubled every 'massdobl' # of bins
xjo=tmp1/LOG(2.)
do i=1,nbins+1 !*NOTE: nbins+1 needed to numerically integrate Fgam
 m=mmin*exp(3*float(i-1)/(xjo))
 rdbF94(i)=2.*((3.*m/4./pi314)**ii3)  !bin diameter bounds
enddo

!--set up avg bin diams, term vel's, and mass
do i=1,nbins
 radF94(i)=0.5*(rdbF94(i)+rdbF94(i+1))    !avg bin diam
 avtF94(i)=cfvt(2)*(radF94(i))**pwvt(2)   !term vel of avg diam
 ramF94(i)=cfmas(2)*(radF94(i))**pwmas(2) !mass of drop w/ avg diam
enddo

return
END SUBROUTINE setupF94
