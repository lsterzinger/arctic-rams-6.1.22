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

Module mem_radiate

implicit none

   Type radiate_vars
   
      ! Variables to be dimensioned by (nzp,nxp,nyp)
   real, allocatable, dimension(:,:,:) :: &
                          fthrd,bext,swup,swdn,lwup,lwdn
                          
      ! Variables to be dimensioned by (nxp,nyp)
   real, allocatable, dimension(:,:) :: &
                          rshort,rlong,rlongup,albedt,cosz

   End Type
   
   type (radiate_vars), allocatable :: radiate_g(:), radiatem_g(:)
   
   integer :: lonrad,ilwrtyp,iswrtyp
   real    :: radfrq
  
Contains

!##############################################################################
Subroutine alloc_radiate (radiate,n1,n2,n3)

implicit none

   type (radiate_vars) :: radiate
   integer, intent(in) :: n1,n2,n3

! Allocate arrays based on options (if necessary)
      
      if(ilwrtyp+iswrtyp > 0)  then
                         allocate (radiate%fthrd(n1,n2,n3))
                         allocate (radiate%rshort(n2,n3))
                         allocate (radiate%rlong(n2,n3))
                         allocate (radiate%rlongup(n2,n3))
                         allocate (radiate%albedt(n2,n3))
                         allocate (radiate%cosz(n2,n3))
      endif
      if(ilwrtyp == 3 .or. iswrtyp == 3) then
         allocate (radiate%bext(n1,n2,n3))
      endif
      if(ilwrtyp == 3)  then
         allocate (radiate%lwup(n1,n2,n3))
         allocate (radiate%lwdn(n1,n2,n3))
      endif
      if(iswrtyp == 3)  then
         allocate (radiate%swup(n1,n2,n3))
         allocate (radiate%swdn(n1,n2,n3))
      endif
              
return
END SUBROUTINE alloc_radiate

!##############################################################################
Subroutine dealloc_radiate (radiate)

implicit none

   type (radiate_vars) :: radiate

   if (allocated(radiate%fthrd))    deallocate (radiate%fthrd)
   if (allocated(radiate%rshort))   deallocate (radiate%rshort)
   if (allocated(radiate%rlong))    deallocate (radiate%rlong)
   if (allocated(radiate%rlongup))  deallocate (radiate%rlongup)
   if (allocated(radiate%albedt))   deallocate (radiate%albedt)
   if (allocated(radiate%cosz))     deallocate (radiate%cosz)
   if (allocated(radiate%bext))     deallocate (radiate%bext)
   if (allocated(radiate%swup))     deallocate (radiate%swup)
   if (allocated(radiate%swdn))     deallocate (radiate%swdn)
   if (allocated(radiate%lwup))     deallocate (radiate%lwup)
   if (allocated(radiate%lwdn))     deallocate (radiate%lwdn)

return
END SUBROUTINE dealloc_radiate

!##############################################################################
Subroutine filltab_radiate (radiate,radiatem,imean,n1,n2,n3,ng)

use var_tables

implicit none

   type (radiate_vars) :: radiate,radiatem
   integer, intent(in) :: imean,n1,n2,n3,ng
   integer :: npts

! Fill arrays into variable tables

   npts=n1*n2*n3

   if (allocated(radiate%fthrd))  &
      CALL vtables2 (radiate%fthrd(1,1,1),radiatem%fthrd(1,1,1)  &
                 ,ng, npts, imean,  &
                 'FTHRD :3:anal:mpti:mpt3')
   if (allocated(radiate%bext))  &
      CALL vtables2 (radiate%bext(1,1,1),radiatem%bext(1,1,1)  &
                 ,ng, npts, imean,  &
                 'BEXT :3:anal:mpti:mpt3')
   if (allocated(radiate%swup))  &
      CALL vtables2 (radiate%swup(1,1,1),radiatem%swup(1,1,1)  &
                 ,ng, npts, imean,  &
                 'SWUP :3:anal:mpti:mpt3')
   if (allocated(radiate%swdn))  &
      CALL vtables2 (radiate%swdn(1,1,1),radiatem%swdn(1,1,1)  &
                 ,ng, npts, imean,  &
                 'SWDN :3:anal:mpti:mpt3')
   if (allocated(radiate%lwup))  &
      CALL vtables2 (radiate%lwup(1,1,1),radiatem%lwup(1,1,1)  &
                 ,ng, npts, imean,  &
                 'LWUP :3:anal:mpti:mpt3')
   if (allocated(radiate%lwdn))  &
      CALL vtables2 (radiate%lwdn(1,1,1),radiatem%lwdn(1,1,1)  &
                 ,ng, npts, imean,  &
                 'LWDN :3:anal:mpti:mpt3')


   npts=n2*n3
   if (allocated(radiate%rshort))  &
      CALL vtables2 (radiate%rshort(1,1),radiatem%rshort(1,1)  &
                 ,ng, npts, imean,  &
                 'RSHORT :2:anal:mpti:mpt3')
   if (allocated(radiate%rlong))  &
      CALL vtables2 (radiate%rlong(1,1),radiatem%rlong(1,1)  &
                 ,ng, npts, imean,  &
                 'RLONG :2:anal:mpti:mpt3')
   if (allocated(radiate%rlongup))  &
      CALL vtables2 (radiate%rlongup(1,1),radiatem%rlongup(1,1)  &
                 ,ng, npts, imean,  &
                 'RLONGUP :2:anal:mpti:mpt3')
   if (allocated(radiate%albedt))  &
      CALL vtables2 (radiate%albedt(1,1),radiatem%albedt(1,1)  &
                 ,ng, npts, imean,  &
                 'ALBEDT :2:anal:mpti:mpt3')
   if (allocated(radiate%cosz))  &
      CALL vtables2 (radiate%cosz(1,1),radiatem%cosz(1,1)  &
                 ,ng, npts, imean,  &
                 'COSZ :2:anal:mpti:mpt3')

return
END SUBROUTINE filltab_radiate

END MODULE mem_radiate
