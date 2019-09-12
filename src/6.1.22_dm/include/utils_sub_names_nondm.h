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

#define c_listfile         c_listfile_
#define irsleep            irsleep_
#define par_init_fortran   par_init_fortran_
#define par_init_put       par_init_put_
#define par_send           par_send_
#define par_put_int        par_put_int_
#define par_put_float      par_put_float_
#define par_put_char       par_put_char_
#define par_send_noblock   par_send_noblock_
#define par_get_noblock    par_get_noblock_
#define par_assoc_buff     par_assoc_buff_
#define par_wait           par_wait_
#define par_get_new        par_get_new_
#define par_get_int        par_get_int_
#define par_get_float      par_get_float_
#define par_get_char       par_get_char_
#define par_exit           par_exit_
#define par_pause          par_pause_
#define par_ready          par_ready_
#define fh5f_open          fh5f_open_
#define fh5f_create        fh5f_create_
#define fh5f_close         fh5f_close_
#define fh5d_open          fh5d_open_
#define fh5d_close         fh5d_close_
#define fh5s_get_ndims     fh5s_get_ndims_
#define fh5s_get_dims      fh5s_get_dims_
#define fh5_prepare_read   fh5_prepare_read_
#define fh5d_read          fh5d_read_
#define fh5_close_read     fh5_close_read_
#define fh5_prepare_write  fh5_prepare_write_
#define fh5_write          fh5_write_
#define fh5_close_write    fh5_close_write_

#else

   print*,"You specified machine/OS other than PC_LINUX1"
   print*,"You need to modify filelist.F90 to add your machine/OS"
   stop

#endif
