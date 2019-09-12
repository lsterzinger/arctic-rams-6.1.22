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

Subroutine rams_fill_fld (n1,n2,n3,a,b,bb,prefile,ngrd,cvar &
                         ,iztrans,cdname,cdunits,ivtype,topt)

use mem_grid

implicit none

integer :: n1,n2,n3,ngrd,iztrans,ivtype
real :: a(*),b(*),bb(*),topt(*)
character(len=*) :: prefile,cvar,cdname,cdunits
character(len=24) :: fdname,fdunits

CALL rams_varlib (cvar,n1,n2,n3,ngrd,a,b,prefile,cdname,cdunits)
CALL rams_varinfo (1,ivtype)
CALL rams_varlib ('topt',n1,n2,1,ngrd,topt,b,prefile,fdname,fdunits)

if(iztrans.eq.2.and.ivtype.eq.3) then
   CALL rams_ctrans (n1,n2,n3,a,b,topt,ztn(1,ngrd),zmn(nnzp(1)-1,1))
elseif(iztrans.eq.3.and.ivtype.eq.3) then
   CALL rams_varlib ('pi',n1,n2,n3,ngrd,bb,b,prefile,fdname,fdunits)
   CALL rams_ptrans (n1,n2,n3,a,bb,b)
endif

CALL rams_value_limit (n1,n2,n3,a)

return
END SUBROUTINE rams_fill_fld

!##############################################################################
Subroutine rams_value_limit (n1,n2,n3,a)

use rcommons

implicit none

!This routine keeps values from being too small since some output
!packages cannot handle really small values due to precision.

integer :: i,j,k,n1,n2,n3
real :: a(n1,n2,n3)

do j=1,n2
   do i=1,n1
      do k=1,n3
         if(a(i,j,k)>0.) a(i,j,k)=max(1.e-20,a(i,j,k))
         if(a(i,j,k)<0.) a(i,j,k)=min(-1.e-20,a(i,j,k))
      enddo
   enddo
enddo

return
END SUBROUTINE rams_value_limit

!##############################################################################
Subroutine rams_ptrans (n1,n2,n3,a,pi,b)

use rcommons

implicit none

integer :: i,j,k,kk,n1,n2,n3,np
real :: a(n1,n2,n3),pi(n1,n2,n3),b(n3,*)

do np=1,nplevs
   b(nplevs-np+1,4)=1004.*(float(iplevs(np))/1000.)**.286
enddo

do j=1,n2
   do i=1,n1
      do k=1,n3
         kk=n3-k+1
         b(kk,1)=a(i,j,k)
         b(kk,2)=pi(i,j,k)
      enddo
      CALL htint (n3,b(1,1),b(1,2),nplevs,b(1,3),b(1,4))
      do k=1,nplevs
         a(i,j,nplevs-k+1)=b(k,3)
      enddo
   enddo
enddo

return
END SUBROUTINE rams_ptrans

!##############################################################################
Subroutine rams_ctrans (n1,n2,n3,a,b,topt,zt,ztop)

implicit none

integer :: i,j,k,n1,n2,n3
real :: a(n1,n2,n3),b(n3,*),topt(n1,n2),zt(n3),ztop

do j=1,n2
   do i=1,n1
      do k=1,n3
         b(k,1)=a(i,j,k)
         b(k,2)=topt(i,j)+zt(k)*(1.-topt(i,j)/ztop)
      enddo
      CALL htint (n3,b(1,1),b(1,2),n3,b(1,3),zt)
      do k=1,n3
         a(i,j,k)=b(k,3)
      enddo
   enddo
enddo

return
END SUBROUTINE rams_ctrans

!##############################################################################
Subroutine thelatlon (n1,n2,n3,alatlon,all1,i1,j1)

implicit none

integer :: n1,n2,n3,i1,j1
real :: alatlon(n1,n2,n3),all1

all1=alatlon(i1,j1,1)

return
END SUBROUTINE thelatlon
