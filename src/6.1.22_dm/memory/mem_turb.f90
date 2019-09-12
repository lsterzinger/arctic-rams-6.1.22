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

Module mem_turb

use grid_dims

implicit none

   Type turb_vars
   
      ! Variables to be dimensioned by (nzp,nxp,nyp)
   real, allocatable, dimension(:,:,:) :: &
                          tkep,epsp,hkm,vkm,vkh

      ! Variables to be dimensioned by (nxp,nyp)
   real, allocatable, dimension(:,:) :: &
                          sflux_u,sflux_v,sflux_w,sflux_t,sflux_r

   End Type

   type (turb_vars), allocatable :: turb_g(:), turbm_g(:)

   integer :: ihorgrad

   integer, dimension(maxgrds) :: idiffk

   real, dimension(maxgrds) :: zkhkm,xkhkm,csz,csx,akmin

Contains

!##############################################################################
Subroutine alloc_turb (turb,n1,n2,n3,ngrids)

implicit none

   type (turb_vars) :: turb
   integer, intent(in) :: n1,n2,n3,ngrids
   integer :: i,tke_alloc,epsp_alloc

! Allocate arrays based on options (if necessary)

! Saleeby(2009): For TKEP and EPSP, we need to allocate for all grids
! rather than being grid dependent since we allow grid dependent 
! variations in IDIFFK. There is also no grid dependent checking in
! mem_tend.f90. This can lead to segmentation faults if TKEP and EPSP
! not allocated on all grids. Allocate on all grids if even a single
! grid requires the allocation

      tke_alloc = 0
      do i=1,ngrids
        if(idiffk(i) == 1 .or. idiffk(i) == 4 .or.  &
           idiffk(i) == 5 .or. idiffk(i) == 6)  &
            tke_alloc = 1
      enddo   
      if(tke_alloc == 1) allocate (turb%tkep(n1,n2,n3))

      epsp_alloc = 0
      do i=1,ngrids
        if(idiffk(i) == 5 .or. idiffk(i) == 6)  &
            epsp_alloc = 1
      enddo   
      if(epsp_alloc == 1) allocate (turb%epsp(n1,n2,n3))

                         allocate (turb%hkm(n1,n2,n3))
                         allocate (turb%vkm(n1,n2,n3))
                         allocate (turb%vkh(n1,n2,n3))

                         allocate (turb%sflux_u(n2,n3))
                         allocate (turb%sflux_v(n2,n3))
                         allocate (turb%sflux_w(n2,n3))
                         allocate (turb%sflux_t(n2,n3))
                         allocate (turb%sflux_r(n2,n3))
                         
return
END SUBROUTINE alloc_turb

!##############################################################################
Subroutine dealloc_turb (turb)

implicit none

   type (turb_vars) :: turb

   if (allocated(turb%tkep))    deallocate (turb%tkep)
   if (allocated(turb%epsp))    deallocate (turb%epsp)
   if (allocated(turb%hkm))     deallocate (turb%hkm)
   if (allocated(turb%vkm))     deallocate (turb%vkm)
   if (allocated(turb%vkh))     deallocate (turb%vkh)
   if (allocated(turb%sflux_u)) deallocate (turb%sflux_u)
   if (allocated(turb%sflux_v)) deallocate (turb%sflux_v)
   if (allocated(turb%sflux_w)) deallocate (turb%sflux_w)
   if (allocated(turb%sflux_t)) deallocate (turb%sflux_t)
   if (allocated(turb%sflux_r)) deallocate (turb%sflux_r)

return
END SUBROUTINE dealloc_turb

!##############################################################################
Subroutine filltab_turb (turb,turbm,imean,n1,n2,n3,ng)

use var_tables

implicit none

   type (turb_vars) :: turb,turbm
   integer, intent(in) :: imean,n1,n2,n3,ng
   integer :: npts

! Fill arrays into variable tables

   npts=n1*n2*n3

   if (allocated(turb%tkep))  &
      CALL vtables2 (turb%tkep(1,1,1),turbm%tkep(1,1,1)  &
                 ,ng, npts, imean,  &
                 'TKEP :3:anal:mpti:mpt3:mpt1')
   if (allocated(turb%epsp))  &
      CALL vtables2 (turb%epsp(1,1,1),turbm%epsp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'EPSP :3:anal:mpti:mpt3:mpt1')

   if (allocated(turb%hkm))  &
      CALL vtables2 (turb%hkm(1,1,1),turbm%hkm(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RHKM :3:anal:mpti:mpt3:mpt1')
   if (allocated(turb%vkm))  &
      CALL vtables2 (turb%vkm(1,1,1),turbm%vkm(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RVKM :3:anal:mpti:mpt3:mpt1')
   if (allocated(turb%vkh))  &
      CALL vtables2 (turb%vkh(1,1,1),turbm%vkh(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RVKH :3:anal:mpti:mpt3:mpt1')

   npts=n2*n3
   if (allocated(turb%sflux_u))  &
      CALL vtables2 (turb%sflux_u(1,1),turbm%sflux_u(1,1)  &
                 ,ng, npts, imean,  &
                 'SFLUX_U :2:anal:mpt3:mpt1')
   if (allocated(turb%sflux_v))  &
      CALL vtables2 (turb%sflux_v(1,1),turbm%sflux_v(1,1)  &
                 ,ng, npts, imean,  &
                 'SFLUX_V :2:anal:mpt3:mpt1')
   if (allocated(turb%sflux_w))  &
      CALL vtables2 (turb%sflux_w(1,1),turbm%sflux_w(1,1)  &
                 ,ng, npts, imean,  &
                 'SFLUX_W :2:anal:mpt3')
   if (allocated(turb%sflux_t))  &
      CALL vtables2 (turb%sflux_t(1,1),turbm%sflux_t(1,1)  &
                 ,ng, npts, imean,  &
                 'SFLUX_T :2:anal:mpt3')
   if (allocated(turb%sflux_r))  &
      CALL vtables2 (turb%sflux_r(1,1),turbm%sflux_r(1,1)  &
                 ,ng, npts, imean,  &
                 'SFLUX_R :2:anal:mpt3')

return
END SUBROUTINE filltab_turb

END MODULE mem_turb
