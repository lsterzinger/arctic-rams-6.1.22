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

Subroutine fdback (ac,af,dc,df,nzc,nxc,nyc,nzf,nxf,nyf,nf,vnam,sumflg)

use mem_grid

implicit none

integer :: nzc,nxc,nyc,nzf,nxf,nyf,nf,ibeg,jbeg,kbeg,iend,jend,kend  &
   ,nc,iinc,jinc,kv,ifbcf,ic,jc,kc,if,jf,kf
real, dimension(nzc,nxc,nyc) :: ac,dc,sumflg
real, dimension(nzf,nxf,nyf) :: af,df
character(len=*) :: vnam

nc = nxtnest(nf)
CALL azero (nzc*nxc*nyc,sumflg)
ibeg = 2
jbeg = 1 + jdim
kbeg = 2
iend = nxf - 1
jend = nyf - jdim
kend = nzf - 1
iinc = 1
jinc = 1
kv = 0

if (vnam .eq. 'u') then
   ibeg = 1 + nstratx(nf)
   iend = nxf - 1 - nstratx(nf)
   iinc = nstratx(nf)
   ifbcf = 1
elseif (vnam .eq. 'v') then
   jbeg = 1 + nstraty(nf) * jdim
   jend = nyf - (1 + nstraty(nf)) * jdim
   jinc = nstraty(nf)
   ifbcf = 2
elseif (vnam .eq. 'w' .or. vnam .eq. 'terr') then
   if (vnam .eq. 'w') then
      kbeg = 1 + nrz(kpm(2,nf),nf)
      kend = nzf - 1 - nrz(kpm(nzf-1,nf),nf)
      kv = 1
   else
      kbeg = 1
      kend = 1
   endif
   ifbcf = 3
else
   ifbcf = 4
endif

!print*,'fdback:',vnam,':',ibeg,iend,ipm(ibeg,nf),ipm(iend,nf)

kf = kbeg
1 continue
   kc = kpm(kf,nf)
   if (vnam .eq. 'terr') kc = 1

   do jf = jbeg,jend,jinc
      jc = jpm(jf,nf)
      if (vnam .eq. 'p' .or. vnam .eq. 'terr') then
         do if = ibeg,iend,iinc
            ic = ipm(if,nf)
            ac(kc,ic,jc) = ac(kc,ic,jc) * sumflg(kc,ic,jc)  &
                         + af(kf,if,jf) * fbcf(kf,nf,ifbcf)
            sumflg(kc,ic,jc) = 1.
         enddo
      elseif (vnam .eq. 'w') then
         do if = ibeg,iend,iinc
            ic = ipm(if,nf)
            ac(kc,ic,jc) = ac(kc,ic,jc) * sumflg(kc,ic,jc)  &
               + af(kf,if,jf) * fbcf(kf,nf,ifbcf)  &
               * (df(kf,if,jf) + df(kf+kv,if,jf))
            sumflg(kc,ic,jc) = 1.
         enddo
      else
         do if = ibeg,iend,iinc
            ic = ipm(if,nf)
            ac(kc,ic,jc) = ac(kc,ic,jc) * sumflg(kc,ic,jc)  &
               + af(kf,if,jf) * fbcf(kf,nf,ifbcf) * df(kf,if,jf)
            sumflg(kc,ic,jc) = 1.
         enddo
      endif
   enddo
   kf = kf + 1
   if (vnam .eq. 'w') kf = kf + nrz(kpm(kf,nf),nf) - 1
if (kf .le. kend) go to 1

do kc = kpm(kbeg,nf),kpm(kend,nf)
   do jc = jpm(jbeg,nf),jpm(jend,nf)

      if (vnam .eq. 'w') then

         do ic = ipm(ibeg,nf),ipm(iend,nf)
            ac(kc,ic,jc) = ac(kc,ic,jc) / (dc(kc,ic,jc) + dc(kc+kv,ic,jc))
         enddo

      elseif (vnam .ne. 'p' .and. vnam .ne. 'terr') then

         do ic = ipm(ibeg,nf),ipm(iend,nf)
            ac(kc,ic,jc) = ac(kc,ic,jc) / dc(kc,ic,jc)
         enddo

      endif
   enddo
enddo

return
END SUBROUTINE fdback

!##############################################################################
Subroutine fdbackp (ivarn,af,ac,mtp,df,dfu,dfv,m1,m2,m3,ifm  &
   ,ia,iz,ja,jz,i0,j0,nestratx,nestraty)

use mem_grid

implicit none

integer :: ivarn,mtp,m1,m2,m3,ifm,ia,iz,ja,jz,i0,j0,nestratx,nestraty
real, dimension(*) :: ac
real, dimension(m1,m2,m3) :: af,df,dfu,dfv   
integer :: kptscm,ikptscm,ibeg,jbeg,kbeg  &
   ,iend,jend,kend,iinc,jinc,kv,ifbcf,if,jf,jc,kf,kc,indcf


!     ivarn = variable types 1- u
!                            2- v
!                            3- w
!                            4- p
!                            5- scalar

! Local variables ia, iz, ja, and jz in this routine refer not to the
! limits of prognostic points but to the limits of points on this fm node
! to be averaged to a given cm node.


kptscm = kpm(m1-1,ifm) - kpm(2,ifm) + 1
ikptscm = kptscm * (ipm(iz+i0,ifm) - ipm(ia+i0,ifm)+1)
ibeg = ia
jbeg = ja
kbeg = 2
iend = iz
jend = jz
kend = m1 - 1
iinc = 1
jinc = 1
kv = 0

CALL azero (mtp,ac)

if (ivarn .eq. 1) then
   do ibeg = ia+i0,ia+i0+nestratx
      if (mod(ibeg-1,nestratx) .eq. 0) go to 11
   enddo
11      continue
   ibeg = ibeg - i0
   iinc = nestratx
   ifbcf = 1
elseif (ivarn .eq. 2) then
   do jbeg = ja+j0,ja+j0+nestraty
      if (mod(jbeg-1,nestraty) .eq. 0) go to 10
   enddo
10      continue
   jbeg = jbeg - j0
   jinc = nestraty
   ifbcf = 2
elseif (ivarn .eq. 3) then
   kbeg = 1 + nrz(kpm(2,ifm),ifm)
   kend = m1 - 2
   kv = 1
   ifbcf = 3
else
   ifbcf = 4
endif

kf = kbeg
1 continue
kc = kpm(kf,ifm) - kpm(2,ifm) + 1
if (ivarn .eq. 4) then
   do jf = jbeg,jend,jinc
      jc = (jpm(jf+j0,ifm) - jpm(ja+j0,ifm)) * ikptscm
      do if = ibeg,iend,iinc
         indcf = kc + jc + (ipm(if+i0,ifm) - ipm(ia+i0,ifm)) * kptscm
         ac(indcf) = ac(indcf) + af(kf,if,jf) * fbcf(kf,ifm,ifbcf)
      enddo
   enddo
elseif (ivarn .eq. 3) then
   do jf = jbeg,jend,jinc
      jc = (jpm(jf+j0,ifm) - jpm(ja+j0,ifm)) * ikptscm
      do if = ibeg,iend,iinc
         indcf = kc + jc + (ipm(if+i0,ifm) - ipm(ia+i0,ifm)) * kptscm
         ac(indcf) = ac(indcf) + af(kf,if,jf) * fbcf(kf,ifm,ifbcf)  &
            * (df(kf,if,jf) + df(kf+kv,if,jf))
      enddo
   enddo
elseif (ivarn .eq. 2) then
   do jf = jbeg,jend,jinc
      jc = (jpm(jf+j0,ifm) - jpm(ja+j0,ifm)) * ikptscm
      do if = ibeg,iend,iinc
         indcf = kc + jc + (ipm(if+i0,ifm) - ipm(ia+i0,ifm)) * kptscm
         ac(indcf) = ac(indcf) + af(kf,if,jf) * fbcf(kf,ifm,ifbcf)  &
            * dfv(kf,if,jf)
      enddo
   enddo
elseif (ivarn .eq. 1) then
   do jf = jbeg,jend,jinc
      jc = (jpm(jf+j0,ifm) - jpm(ja+j0,ifm)) * ikptscm
      do if = ibeg,iend,iinc
         indcf = kc + jc + (ipm(if+i0,ifm) - ipm(ia+i0,ifm)) * kptscm
         ac(indcf) = ac(indcf) + af(kf,if,jf) * fbcf(kf,ifm,ifbcf)  &
            * dfu(kf,if,jf)
      enddo
   enddo
else
   do jf = jbeg,jend,jinc
      jc = (jpm(jf+j0,ifm) - jpm(ja+j0,ifm)) * ikptscm
      do if = ibeg,iend,iinc
         indcf = kc + jc + (ipm(if+i0,ifm) - ipm(ia+i0,ifm)) * kptscm
         ac(indcf) = ac(indcf) + af(kf,if,jf) * fbcf(kf,ifm,ifbcf)  &
            * df(kf,if,jf)
      enddo
   enddo
endif
kf = kf + 1
if (ivarn .eq. 3) kf = kf + nrz(kpm(kf,ifm),ifm) - 1

if (kf .le. kend) go to 1

return
END SUBROUTINE fdbackp

