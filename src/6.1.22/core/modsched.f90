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

Subroutine modsched (ischtab,ntab,ngrids,nxtnest,nndtrat,nsub)

implicit none

integer :: ntab,ngrids,nsub,ischtab(ntab,*),nxtnest(*),nndtrat(*)  &
          ,ncnt(50),nnfm(50)
integer :: ng,nt,is,itopgrd,ntt,ngrid,isstp,n,k

do ng=1,ngrids
   nnfm(ng)=1
   ncnt(ng)=0
enddo
do nt=1,4
   do ntt=1,ntab
      ischtab(ntt,nt)=0
   enddo
enddo

is=0
itopgrd = 0
do ng=1,ngrids
   if(nxtnest(ng).eq.0)then
      ngrid=ng

30         ISSTP=MOD(NCNT(NGRID),NNDTRAT(NGRID))+1
      is=is+1

!        Table value 1 is the grid number to be stepped forward
      ischtab(is,1)=ngrid
!        Table value 3 is the sub-timestep counter
      ischtab(is,3)=ISSTP

      NCNT(NGRID)=NCNT(NGRID)+1

!        Table value 5 is the total sub-timestep counter
      ischtab(is,5)=NCNT(NGRID)

40         NNFM(NGRID)=NNFM(NGRID)+1
      IF (NNFM(NGRID).GT.NGRIDS) GO TO 50
      IF (NXTNEST(NNFM(NGRID)).NE.NGRID) GO TO 40
      NGRID=NNFM(NGRID)

!        Table value 2 is the grid number to interpolate. 0 = no interpolate.
      ischtab(is,2)=ngrid

      ISSTP=MOD(NCNT(NGRID),NNDTRAT(NGRID))+1

      GO TO 30

50         NNFM(NGRID)=1
      IF (NXTNEST(NGRID).EQ.0) GO TO 80
      IF(NCNT(NGRID).LT.NCNT(NXTNEST(NGRID))*NNDTRAT(NGRID))THEN
         GO TO 30
      ENDIF

!        Table value 4 is the number of grids to be fedback
      ischtab(is,4)=ischtab(is,4)+1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ngrid=nxtnest(ngrid)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      GO TO 40

80         CONTINUE

   endif
enddo

nsub=is

print*
print*,'=== Timestep Schedule ===='
print 333,(n,(ischtab(n,k),k=1,5),n=1,is)
333 format (6i5)

return
END SUBROUTINE modsched

!##############################################################################
Subroutine cfl (n1,n2,n3,i0,j0,mynum)

use mem_basic
use mem_grid

implicit none

integer :: n1,n2,n3,i0,j0,mynum

CALL cfll (n1,n2,n3,i0,j0,mynum  &
   ,basic_g(ngrid)%up   (1,1,1)  ,basic_g(ngrid)%vp   (1,1,1)  &
   ,basic_g(ngrid)%wp   (1,1,1)  ,grid_g(ngrid)%rtgt    (1,1)  &
   ,grid_g(ngrid)%f13t    (1,1)  ,grid_g(ngrid)%f23t    (1,1)  &
   ,grid_g(ngrid)%dxt     (1,1)  ,grid_g(ngrid)%dyt     (1,1)  )

return
END SUBROUTINE cfl

!##############################################################################
Subroutine cfll (n1,n2,n3,i0,j0,mynum,up,vp,wp,rtgt,f13t,f23t,dxt,dyt)

use mem_grid
use mem_scratch

implicit none

integer :: n1,n2,n3,i0,j0,mynum
real, dimension(n1,n2,n3) :: up,vp,wp
real, dimension(n2,n3)    :: rtgt,f13t,f23t,dxt,dyt
         
integer :: i,j,k,ifm,icm,innest,nprints
real :: c1x,c1y,c1z,cflnumh,cflnumv

!     This routine returns flags the model to bring itself down when the CFL
!     linear stability criteria on advection is exceeded.
!     (Actually check on 90% of CFL)

5    format('cflx,ngrid,k,i,j,mynum = ',f5.1,5i5)
6    format('cfly,ngrid,k,i,j,mynum = ',f5.1,5i5)
7    format('cflz,ngrid,k,i,j,mynum = ',f5.1,5i5)

nprints = 0
cflnumh = .90
cflnumv = .90
cflxy(ngrid) = 0.
cflz(ngrid) = 0.
vertvel(ngrid) = 0.

! Let's try a new thing... if we have a grid point that is on a 
!   coarse grid, but it is under a nested grid, we will ignore it
!   under the assumption that the fine grid values will overwrite
!   it, hence not allowing it to go numerically unstable.

jloop: do j = 1+jdim,n3-jdim
   iloop: do i = 2,n2-1
   
      ! See if this is under a fine grid horizontally... ignore vertical for now
      innest=0
      if (ngrids > ngrid) then
         do ifm=ngrid+1,ngrids
            icm=nxtnest(ifm)
            if(icm == ngrid .and. &
               i+i0 >= ipm(1,ifm) .and. i+i0 <= ipm(nnxp(ifm),ifm) .and. &
               j+j0 >= jpm(1,ifm) .and. j+j0 <= jpm(nnyp(ifm),ifm) ) then
               innest=1
               exit
            endif
         enddo
      endif
      
      if(innest == 1) then
         !print '(a,5i4)', 'cfl under grid-' &
         !         ,ngrid,ifm,i+i0,j+j0,mynum
         cycle iloop
      endif
   
      kloop: do k = 2,n1-1
      
         vctr1(k) = .5*(up(k,i,j)+up(k,i-1,j))*dtlt*dxt(i,j)
         vctr2(k) = .5*(vp(k,i,j)+vp(k,i,j-jdim))*dtlt*dyt(i,j)
         vctr3(k) = ((wp(k,i,j)+wp(k-1,i,j))  &
           +(up(k,i,j)+up(k,i-1,j))*f13t(i,j)*ht(k)*rtgt(i,j)  &
           +(vp(k,i,j)+vp(k,i,j-jdim))*f23t(i,j)*ht(k)*rtgt(i,j)  &
           )*.5*dtlt*dzt(k)
      enddo kloop
      
      do k = 2,n1-1
         c1x = abs(vctr1(k))
         c1y = abs(vctr2(k))
         c1z = abs(vctr3(k))

         if (nprints .le. 10) then
            if (c1x .gt. cflnumh) then
               nprints = nprints + 1
               print 5, c1x,ngrid,k,i+i0,j+j0,mynum
               print*,up(k,i,j),up(k,i-1,j),dtlt,dxt(i,j),vctr1(k)
            endif
            if (c1y .gt. cflnumh) then
               nprints = nprints + 1
               print 6, c1y,ngrid,k,i+i0,j+j0,mynum
            endif
            if (c1z .gt. cflnumv) then
               nprints = nprints + 1
               print 7, c1z,ngrid,k,i+i0,j+j0,mynum
            endif
         endif

         if (c1x .gt. cflxy(ngrid)) cflxy(ngrid) = c1x
         if (c1y .gt. cflxy(ngrid)) cflxy(ngrid) = c1y
         if (c1z .gt. cflz(ngrid)) cflz(ngrid) = c1z
         if (wp(k,i,j) .gt. vertvel(ngrid)) vertvel(ngrid) = wp(k,i,j)
      enddo
   enddo iloop
enddo jloop

return
END SUBROUTINE cfll

!##############################################################################
Subroutine dtset (nndtflg)

use mem_grid
use rconstants
use ref_sounding
use io_params

implicit none

integer :: nndtflg
integer, parameter  :: ndx=37,ndt=42
integer, dimension(maxgrds) :: idelt,nndtrat1,nnacoust1
real, dimension(maxgrds) :: sscourn,dtlongn1
real, dimension(nzpmax) :: vctr1

integer, save, dimension(ndx) :: idelx(maxgrds)
real, save :: cflnumh, cflnumv ,ssodx(maxgrds),delx(ndx),delt(ndt)
data delx/200000.,150000.,100000.,80000.,70000.,60000.,40000.,  &
           30000., 20000., 10000., 6000., 4000., 3000., 2000.,  &
            1000.,   800.,   600.,  500.,  400.,  300.,  200.,  &
             150.,   100.,    80.,   60.,   50.,   40.,   30.,  &
              20.,    10.,     8.,    6.,    5.,    4.,    3.,  &
               2.,     1. /
data delt/   300.,   240.,   180.,  150.,  120.,   90.,   60.,  &
              50.,    40.,    30.,   20.,   15.,   12.,   10.,  &
               6.,     5.,     4.,    3.,   2.5,   2.0,   1.5,  &
              1.2,    1.0,     .8,    .6,    .5,    .4,    .3,  &
               .2,     .1,    .08,   .06,   .05,   .04,   .03,  &
              .02,    .01,   .008,  .006,  .005,  .004,  .003/
data idelx(1)/0/

integer :: iabsdt,ifm,id,n2,n3,k,nn2,nn3,icm,ntf,ii
real :: ssmax,tmax,dxtmax,sscnmax,sspct0,cflxyz,timeleft

sscnmax=0 !Variable initialized
iabsdt = abs(ideltat)

! On the first call to this routine, initialize idelx, ssodx, dtlongn,
! nnacoust, and if required, nndtrat.

if (idelx(1) .eq. 0) then

   cflnumh = .90
   cflnumv = .90

   do ifm = 1,ngrids

      do id = ndx,1,-1
         if (delx(id) .le. deltaxn(ifm)) idelx(ifm) = id
      enddo

      n2 = nnxp(ifm)
      n3 = nnyp(ifm)
      do k = 1,nnzp(ifm)
         vctr1(k) = th01dn(k,1) * pi01dn(k,1) / cp
      enddo
      tmax = maxval(vctr1(1:nnzp(ifm)))
      ssmax = sqrt(cp / cv * rgas * tmax)

      nn2 = nnxp(ifm)
      nn3 = nnyp(ifm)

      dxtmax = max(grid_g(ifm)%dxt(1,1)      &
                  ,grid_g(ifm)%dxt(nn2,1)    &
                  ,grid_g(ifm)%dxt(nn2,nn3)  &
                  ,grid_g(ifm)%dxt(1,nn3)    )
      ssodx(ifm) = ssmax * dxtmax

   enddo

   if (ideltat .eq. 0) then

      sspct = 1.
      do ifm = 1,ngrids
         icm = nxtnest(ifm)
         if (icm .eq. 0) then
            dtlongn(ifm) = dtlong
         else
            dtlongn(ifm) = dtlongn(icm) / nndtrat(ifm)
         endif
         nnacoust(ifm) = nacoust
         sscourn(ifm) = 2. * ssodx(ifm) * dtlongn(ifm)
         sspct = min(sspct,  &
            .95*float(nnacoust(ifm))/(2.*sscourn(ifm)))            
      enddo

   else

      do ifm = 1,ngrids
         icm = nxtnest(ifm)
         dtlongn(ifm) = delt(idelx(ifm)+iabsdt-1)

! For coarse grid(s) adjust dtlongn so that it is an integer divisor of
! FRQSTATE.  For nested grids, compute nndtrat(ifm) as the first integer greater
! than or equal to the timestep ratio between a nested grid's parent and the
! nested grid.  Then adjust dtlongn(ifm) for the nested grid to be the
! parent grid timestep divided by nndtrat(ifm).

         if (icm .eq. 0) then
            ntf = nint(frqstate(1) / dtlongn(1))
            dtlongn(ifm) = frqstate(ifm) / ntf
         else
            nndtrat(ifm) = min(10,  &
               nint(dtlongn(icm) / dtlongn(ifm)))
            dtlongn(ifm) = dtlongn(icm) / nndtrat(ifm)
         endif

! Compute sst courant numbers (sstcourn(ifm)) for long timestep dtlongn.

         sscourn(ifm) = 2. * ssodx(ifm) * dtlongn(ifm)
         if (sscourn(ifm) .gt. sscnmax) sscnmax = sscourn(ifm)

      enddo

! Define trial sspct0 in terms of sscnmax using nonlinear formula intended to
! increase nnacoust as sspct decreases, but limit sspct0 to a minimum of .2.

      sspct0 = min(1.,(.95/sscnmax) ** .5)

      if (sspct0 .lt. .2) then
         print*, 'Sound speed percent is forced to be too low'
         stop 'low_sspct0'
      endif

      sspct = 1.
      do ifm = 1,ngrids
         nnacoust(ifm) = max(2,  &
            nint(sspct0 * sscourn(ifm) * 2. / .95))
         sspct = min(sspct,  &
            .95*float(nnacoust(ifm))/(2.*sscourn(ifm)))
      enddo

   endif

! Print out initial values of dtlongn, nndtrat, nnacoust, sscourn, and sspct


   write(6,120)
120        format(/,'Initial timestep info: ngrid, nndtrat, nnacoust,'  &
            ,' dtlongn, sscourn, sspct')

   do ifm = 1,ngrids

      write(6,121) ifm,nndtrat(ifm),nnacoust(ifm),dtlongn(ifm)  &
         ,sscourn(ifm),sspct
121        format(23x,i3,i8,i9,f13.3,2f8.3)
   enddo
endif

! check Courant numbers

nndtflg = 0

if (ideltat >= 0) then

   do ifm = 1,ngrids
      cflxyz = max(cflxy(ifm)/cflnumh,cflz(ifm)/cflnumv)
      if (cflxyz > 1.) then
         iflag = 1
         print*, 'Model will stop because CFL limit exceeded'
         print*, 'with ideltat .ge. 0.'
      endif
   enddo

   if (ideltat .eq. 0) then

! For the coarse grid(s), adjust dtlongn(1) to not jump past an analysis
! write time or the end of the model run time.

      sspct = 1.
      do ifm = 1,ngrids
         icm = nxtnest(ifm)
         if (ifm == 1) then
            dtlongn(ifm) = dtlong
            timeleft =  &
               min (timmax - time, frqstate(1) - mod(time,frqstate(1)))
            if (dtlongn(1) > timeleft) dtlongn(1) = timeleft
         else
            dtlongn(ifm) = dtlongn(icm) / nndtrat(ifm)
         endif
         nnacoust(ifm) = nacoust
         sscourn(ifm) = 2. * ssodx(ifm) * dtlongn(ifm)
         sspct = min(sspct,  &
            .95*float(nnacoust(ifm))/(2.*sscourn(ifm)))
      enddo
   endif

else

   do ifm = 1,ngrids
      icm = nxtnest(ifm)
      cflxyz = max(cflxy(ifm)/cflnumh,cflz(ifm)/cflnumv)
      nndtrat1(ifm) = nndtrat(ifm)
      dtlongn1(ifm) = dtlongn(ifm)
      nnacoust1(ifm) = nnacoust(ifm)

      do id = ndt,1,-1
         if (delt(id) * cflxyz .le. dtlongn(ifm)) idelt(ifm) = id
      enddo

      if (idelt(ifm) .gt. idelx(ifm)+4) then
         print*, 'Adjustable timestep forced to become too small'
         print*, 'on grid ',ifm,' idelx,idelt = '  &
            ,idelx(ifm),idelt(ifm)
         print*, 'Model will stop'
         iflag = 1
      else
         ii = max(idelx(ifm)+iabsdt-1,idelt(ifm))
         dtlongn(ifm) = delt(ii)

! For the coarse grid(s), adjust dtlongn(1) to not jump past an analysis
! write time or the end of the model run time.

         if (icm .eq. 0) then
            timeleft =  &
               min (timmax - time, frqstate(1) - mod(time,frqstate(1)))
            if (dtlongn(1) .gt. .95 * timeleft) then
               dtlongn(ifm) = timeleft
            endif
         else
            nndtrat(ifm) = min(10,  &
               nint(dtlongn(icm) / dtlongn(ifm)))
            dtlongn(ifm) = dtlongn(icm) / nndtrat(ifm)
         endif
      endif

! Compute sst courant numbers (sstcourn(ifm)) for long timestep dtlongn.

      sscourn(ifm) = 2. * ssodx(ifm) * dtlongn(ifm)
      if (sscourn(ifm) .gt. sscnmax) sscnmax = sscourn(ifm)
   enddo

! Define trial sspct0 in terms of sscnmax using nonlinear formula intended to
! increase nnacoust as sspct decreases, but limit sspct0 to a minimum of .2.

   sspct0 = min(1.,(.95/sscnmax) ** .5)
   if (sspct0 .lt. .2) then
      print*, 'Sound speed percent is forced to be too low'
      stop 'low_sspct0'
   endif

   sspct = 1.
   do ifm = 1,ngrids
      nnacoust(ifm) = max(2,  &
         nint(sspct0 * sscourn(ifm) * 2. / .95))
      sspct = min(sspct,  &
         .95*float(nnacoust(ifm))/(2.*sscourn(ifm)))

! If there are any updates to dtlongn(ifm), print out new values.

      if (abs(dtlongn(ifm) - dtlongn1(ifm)) .gt. 1.e-3) then
         write(6,122) ifm,nndtrat(ifm),nnacoust(ifm),dtlongn(ifm)
122           format('Timestep update: ngrid, nndtrat, nnacoust,'  &
            ,' dtlongn = ',3i3,f10.3)
      endif

! If there are any updates to nndtrat(ifm) or others, set nndtflg = 1 to flag
! new call to routine modsched and send new stuff to nodes.

         if (nndtrat(ifm) /= nndtrat1(ifm) .or. &
             nnacoust(ifm) /= nnacoust1(ifm) .or. &
             dtlongn(ifm) /= dtlongn1(ifm) ) nndtflg = 1

   enddo

endif

return
END SUBROUTINE dtset
