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

Module var_tables

use grid_dims

implicit none

    !    Define data type for main variable table

type var_tables_r
   
   real, pointer :: var_p,var_m
   integer :: npts, idim_type
   integer :: ianal,imean,ilite,impti,impt1,impt2,impt3,irecycle_sfc
   character(len=32) :: name
   
end type

    !    Main variable table
type(var_tables_r) :: vtab_r(maxvars,maxgrds)

    !    "nvgrids" is "ngrids", for convenience
integer :: nvgrids

    !    number of variables for each grid
integer :: num_var(maxgrds)



    !    Define data type for scalar variable table

type scalar_table
   
   real, pointer :: var_p,var_t
   character(len=32) :: name
   
end type

    !    Scalar variable table
type(scalar_table) :: scalar_tab(maxsclr,maxgrds)


    !    number of scalars for each grid
integer :: num_scalar(maxgrds)


END MODULE var_tables
