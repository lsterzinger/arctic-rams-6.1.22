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

Module anal_extra

implicit none

integer, parameter :: max_anextra=10
type anextra
   character(len=16) :: name
   integer :: idim_type, npointer,ind_comp
end type

integer :: num_extra_anal

type(anextra), dimension(max_anextra) :: an_extra

Contains

!##############################################################################
Subroutine anal_extra_init ()

implicit none

num_extra_anal=7
an_extra(1)%name='PI'  ; an_extra(1)%idim_type=3 ; an_extra(1)%ind_comp=1
an_extra(2)%name='HKH' ; an_extra(2)%idim_type=3 ; an_extra(2)%ind_comp=2
an_extra(3)%name='VKH' ; an_extra(3)%idim_type=3 ; an_extra(3)%ind_comp=3
an_extra(4)%name='DN0' ; an_extra(4)%idim_type=3 ; an_extra(4)%ind_comp=4
an_extra(5)%name='SOIL_WATER-1';    an_extra(5)%idim_type=4;  an_extra(5)%ind_comp=5
an_extra(6)%name='SFCWATER_MASS-1'; an_extra(6)%idim_type=5;  an_extra(6)%ind_comp=6
an_extra(7)%name='STOM_RESIST-1';   an_extra(7)%idim_type=6;  an_extra(7)%ind_comp=7

return
END SUBROUTINE anal_extra_init

!##############################################################################
Subroutine anal_extra_comp (nv,a,ngr,skip)

! Compute the "extra" analysis file variables. Note 3-d arrays must be returned
!  as (i,j,k)

use mem_basic
use mem_grid
use mem_turb
use mem_leaf

implicit none

integer :: nv,idtype,ngr
real :: a(*)
logical :: skip

skip = .false.

select case (an_extra(nv)%ind_comp)
   case(1)
      ! Total exner func
      CALL ancomp_pi (nnxyzp(ngr)  &
                 ,basic_g(ngr)%pp,basic_g(ngr)%pi0,a)
   case(2)
      ! Horizontal heat eddy diffusivity
      CALL ancomp_hkh (nnxyzp(ngr),turb_g(ngr)%hkm  &
         ,turb_g(ngr)%vkh,basic_g(ngr)%dn0  &
         ,idiffk(ngr),xkhkm(ngr),a)
   case(3)
      ! Vertical heat eddy diffusivity
      CALL ancomp_vkh (nnxyzp(ngr),turb_g(ngr)%vkh  &
             ,basic_g(ngr)%dn0,a)
   case(4)
      ! Base state density
      CALL atob (nnxyzp(ngr),basic_g(ngr)%dn0,a)
   case(5)
      ! Check; add new here 4D leaf
      !CALL ancomp_check (nnxyzp(ngr),leaf_g(ngr)%soil_water,a)
      skip = .true.
   case(6)
      ! Check; add new here 5D leaf
      !CALL ancomp_check (nnxyzp(ngr),leaf_g(ngr)%sfcwater_mass,a)
      skip = .true.
   case(7)
      ! Check; add new here 6D leaf
      !CALL ancomp_check (nnxyzp(ngr),leaf_g(ngr)%stom_resist,a)
      skip = .true.

end select

return
END SUBROUTINE anal_extra_comp

!##############################################################################
Subroutine ancomp_check (n1,b,a)

implicit none

integer :: n1
real, dimension(n1) :: a,b
integer :: i

do i=1,n1
   a(i)=b(i)
enddo

return
END SUBROUTINE ancomp_check

!##############################################################################
Subroutine ancomp_pi (n1,b,c,a)

implicit none

integer :: n1
real, dimension(n1) :: a,b,c
integer :: i
 
do i=1,n1
   a(i)=b(i)+c(i)
enddo

return
END SUBROUTINE ancomp_pi

!##############################################################################
Subroutine ancomp_hkh (n1,hkm,vkh,dn0,idiffk,xkhkm,a)

implicit none

integer :: n1,idiffk
real :: xkhkm
real, dimension(n1) :: hkm,vkh,dn0,a
integer :: ind

! Convert to HKM to HKH (note that VKH is HKH for Deardorff)

if (idiffk <= 3) then
   do ind = 1,n1
      a(ind) = hkm(ind) * xkhkm / dn0(ind)
   enddo
elseif (idiffk >= 4) then
   do ind = 1,n1
      a(ind) = vkh(ind) / dn0(ind)
   enddo
endif

return
END SUBROUTINE ancomp_hkh

!##############################################################################
Subroutine ancomp_vkh (n1,vkh,dn0,a)

implicit none

integer :: n1
real :: vkh(n1),dn0(n1),a(n1)
integer :: ind

! Un-density weight VKH

do ind = 1,n1
   a(ind) = vkh(ind) / dn0(ind)
enddo

return
END SUBROUTINE ancomp_vkh

END MODULE anal_extra
