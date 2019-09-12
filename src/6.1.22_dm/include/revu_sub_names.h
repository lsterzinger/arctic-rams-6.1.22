/*
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
*/

#if defined (PC_LINUX1)

/* LINUX APPENDS AN UNDERSCORE TO C ROUTINES */ 
/* Other machines/OS might require no-underscore or all capitals */

#define rh5f_open            rh5f_open_
#define rh5f_close           rh5f_close_
#define rh5f_create          rh5f_create_
#define rh5l_exists          rh5l_exists_
#define rh5g_create          rh5g_create_
#define rh5g_open            rh5g_open_
#define rh5g_close           rh5g_close_
#define rh5d_open            rh5d_open_
#define rh5d_write           rh5d_write_
#define rh5d_setup_and_write rh5d_setup_and_write_
#define rh5d_read_get_dims   rh5d_read_get_dims_
#define rh5d_read            rh5d_read_
#define rh5d_close           rh5d_close_
#define rh5t_close           rh5t_close_
#define rh5s_close           rh5s_close_
#define rh5p_close           rh5p_close_
#define rh5a_write_anyscalar rh5a_write_anyscalar_
#define rh5a_read_anyscalar  rh5a_read_anyscalar_
#define rh5ds_attach_scale   rh5ds_attach_scale_
#define rh5ds_set_scale      rh5ds_set_scale_

#else

   print*,"You specified machine/OS other than PC_LINUX1"
   print*,"You need to modify filelist.F90 to add your machine/OS"
   stop

#endif
