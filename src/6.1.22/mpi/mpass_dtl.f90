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

Subroutine master_getcflcpu ()

use mem_grid
use rpara
use mem_turb

implicit none

integer :: ngr,nm,ibytes,msgtyp,ihostnum
real, save, allocatable :: buff(:),buff1(:),buff2(:)
integer, save :: ncall=0,nwords

if (ncall==0) then
   nwords = 3 + 3 * maxgrds
   allocate (buff(nwords))
   allocate (buff1(maxgrds))
   allocate (buff2(maxgrds))
   ncall=1
endif

do ngr = 1,ngrids
   cflxy(ngr) = 0.
   cflz(ngr) = 0.
enddo

do nm=1,nmachs
   CALL par_get_new (buff,nwords,999,ibytes,msgtyp,ihostnum)
   CALL par_get_float (ptimes(ihostnum,1),1)
   CALL par_get_float (ptimes(ihostnum,2),1)
   CALL par_get_int   (ngbegun,maxgrds)
   CALL par_get_float (eps,1)
   CALL par_get_float (buff1,maxgrds)
   CALL par_get_float (buff2,maxgrds)
   do ngr = 1,ngrids
      cflxy(ngr) = max(cflxy(ngr),buff1(ngr))
      cflz(ngr) = max(cflz(ngr),buff2(ngr))
   enddo
   !print*,'master got:',nm,ptimes(ihostnum,1),ptimes(ihostnum,2)
enddo

return
END SUBROUTINE master_getcflcpu

!##############################################################################
Subroutine node_putcflcpu (totcpu,totwall)

use mem_grid
use node_mod
use mem_turb

implicit none

real :: totcpu,totwall
real, save, allocatable :: buff(:)
integer, save :: ncall=0,nwords

if (ncall==0) then
   nwords = 3 + 3 * maxgrds
   allocate (buff(nwords))
   ncall=1
endif

CALL par_init_put  (buff,nwords)
CALL par_put_float (totcpu,1)
CALL par_put_float (totwall,1)
CALL par_put_int   (ngbegun,maxgrds)
CALL par_put_float (eps,1)
CALL par_put_float (cflxy,maxgrds)
CALL par_put_float (cflz,maxgrds)
CALL par_send (master_num,999)

return
END SUBROUTINE node_putcflcpu

!##############################################################################
Subroutine master_putdtsched (isendflg,isendlite,isendmean  &
                            ,isendboth,ntsend)
use mem_grid
use rpara

implicit none

integer :: isendflg,isendlite,isendmean,isendboth,ntsend
real, save, allocatable :: buff(:)
integer, save :: ncall=0,nwords

integer :: nm

if (ncall==0) then
   nwords = 9 + 4 * maxgrds + maxsched * maxschent
   allocate (buff(nwords))
   ncall=1
endif

CALL par_init_put (buff,nwords)

CALL par_put_int (isendflg,1)
CALL par_put_int (isendlite,1)
CALL par_put_int (isendmean,1)
CALL par_put_int (isendboth,1)

CALL par_put_int (ntsend,1)
if(ntsend == 1) then
   CALL par_put_int   (nnacoust,maxgrds)
   CALL par_put_int   (nndtrat,maxgrds)
   CALL par_put_int   (ngbegun,maxgrds)
   CALL par_put_int   (isched,maxsched*maxschent)
   CALL par_put_int   (nsubs,1)
   CALL par_put_float (dtlongn,maxgrds)
   CALL par_put_float (sspct,1)
   CALL par_put_float (eps,1)
   CALL par_put_int   (initorig,1)
endif

do nm=1,nmachs
   CALL par_send (machnum(nm),44)
enddo

return
END SUBROUTINE master_putdtsched

!##############################################################################
Subroutine node_getdtsched (isendflg,isendlite,isendmean,isendboth)

use mem_grid
use node_mod

implicit none

integer :: isendflg,isendlite,isendmean,isendboth,ntsend
real, save, allocatable :: buff(:)
integer, save :: ncall=0,nwords
integer :: ibytes,msgtype,ihostnum

if (ncall==0) then
   nwords = 9 + 4 * maxgrds + maxsched * maxschent
   allocate (buff(nwords))
   ncall=1
endif

CALL par_get_new (buff,nwords,44,ibytes,msgtype,ihostnum)

CALL par_get_int (isendflg,1)
CALL par_get_int (isendlite,1)
CALL par_get_int (isendmean,1)
CALL par_get_int (isendboth,1)
CALL par_get_int (ntsend,1)

if(ntsend == 1) then
   CALL par_get_int   (nnacoust,maxgrds)
   CALL par_get_int   (nndtrat,maxgrds)
   CALL par_get_int   (ngbegun,maxgrds)
   CALL par_get_int   (isched,maxsched*maxschent)
   CALL par_get_int   (nsubs,1)
   CALL par_get_float (dtlongn,maxgrds)
   CALL par_get_float (sspct,1)
   CALL par_get_float (eps,1)
   CALL par_get_int   (initorig,1)
endif

return
END SUBROUTINE node_getdtsched

!##############################################################################
Subroutine master_putvertvel ()

use mem_grid
use rpara

implicit none

real, save, allocatable :: buff(:)
integer, save :: ncall=0,nwords

integer :: nm

if (ncall==0) then
   nwords = 10 + 1 * maxgrds
   allocate (buff(nwords))
   ncall=1
endif

CALL par_init_put (buff,nwords)

CALL par_put_float (vertvel,maxgrds)

do nm=1,nmachs
   CALL par_send (machnum(nm),45)
enddo

return
END SUBROUTINE master_putvertvel

!##############################################################################
Subroutine node_getvertvel ()

use mem_grid
use node_mod

implicit none

real, save, allocatable :: buff(:)
integer, save :: ncall=0,nwords
integer :: ibytes,msgtype,ihostnum

if (ncall==0) then
   nwords = 10 + 1 * maxgrds
   allocate (buff(nwords))
   ncall=1
endif

CALL par_get_new (buff,nwords,45,ibytes,msgtype,ihostnum)

CALL par_get_float (vertvel,maxgrds)

return
END SUBROUTINE node_getvertvel

!##############################################################################
Subroutine master_getvertvel ()

use mem_grid
use rpara
use mem_turb

implicit none

integer :: ngr,nm,ibytes,msgtyp,ihostnum
real, save, allocatable :: buff(:),buff1(:)
integer, save :: ncall=0,nwords

if (ncall==0) then
   nwords = 3 + 1 * maxgrds
   allocate (buff(nwords))
   allocate (buff1(maxgrds))
   ncall=1
endif

do ngr = 1,ngrids
   vertvel(ngr) = 0.
enddo

do nm=1,nmachs
   CALL par_get_new (buff,nwords,998,ibytes,msgtyp,ihostnum)
   CALL par_get_float (buff1,maxgrds)
   do ngr = 1,ngrids
      vertvel(ngr) = max(vertvel(ngr),buff1(ngr))
   enddo
enddo

return
END SUBROUTINE master_getvertvel

!##############################################################################
Subroutine node_putvertvel ()

use mem_grid
use node_mod

implicit none

real, save, allocatable :: buff(:)
integer, save :: ncall=0,nwords

if (ncall==0) then
   nwords = 3 + 1 * maxgrds
   allocate (buff(nwords))
   ncall=1
endif

CALL par_init_put  (buff,nwords)
CALL par_put_float (vertvel,maxgrds)
CALL par_send (master_num,998)

return
END SUBROUTINE node_putvertvel
