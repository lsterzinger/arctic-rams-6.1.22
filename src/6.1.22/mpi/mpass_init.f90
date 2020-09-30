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

Subroutine masterput_processid (nproc,taskids,master_num)

use rpara

implicit none

integer :: taskids(*),master_num,nproc
!   +------------------------------------------------------------------
!   ! This routine gives basic processor ID info to the nodes.
!   +------------------------------------------------------------------

real, allocatable :: buff(:)
integer :: nm,nwords

mainnum=master_num
nmachs=nproc

nwords= 4 + nmachs
allocate(buff(nwords))

do nm=1,nmachs
   machnum(nm)=taskids(nm)
enddo

do nm=1,nmachs
   CALL par_init_put (buff,nwords)
   CALL par_put_int (mainnum,1)
   CALL par_put_int (machnum(nm),1)
   CALL par_put_int (nm,1)
   CALL par_put_int (nmachs,1)
   CALL par_put_int (machnum(1),nmachs)
   CALL par_send (machnum(nm),1)
enddo

deallocate(buff)

return
END SUBROUTINE masterput_processid

!##############################################################################
Subroutine masterput_nl ()

use mem_all
use rpara

implicit none

real, allocatable :: buff(:)
integer :: nwords, nm

!Saleeby(2016)
!Increment memory buffer size here if you add RAMSIN Namelist variables.
!Add to the appropriate section below as (#-of-them * arraysize).
nwords = 133 * 1                 & !single values
       +   2 * 8                 & !micro (8-hydromet types for gnu)
       +   5 * 9                 & !micro (9-aerosol species)
       +  28 * maxgrds           & !grid-dependent (max grids)
       +   4 * nzpmax            & !max vertical levels
       +   3 * nzgmax            & !max soil levels
       +   1 * maxlite * 32      & !lite variables 32 char length strings
       + 103                       !extras so we have enough buffer

allocate (buff(nwords))

CALL par_init_put (buff,nwords)

CALL par_put_float (TIMMAX,1)
CALL par_put_int   (IMONTH1,1)
CALL par_put_int   (IDATE1,1)
CALL par_put_int   (IYEAR1,1)
CALL par_put_int   (ITIME1,1)
CALL par_put_int   (NGRIDS,1)
CALL par_put_int   (NNXP,MAXGRDS)
CALL par_put_int   (NNYP,MAXGRDS)
CALL par_put_int   (NNZP,MAXGRDS)
CALL par_put_int   (NZG,1)
CALL par_put_int   (NZS,1)
CALL par_put_int   (NXTNEST,MAXGRDS)
CALL par_put_int   (INESTING,1)
CALL par_put_int   (IPRNTSTMT,1)
CALL par_put_int   (IHTRAN,1)
CALL par_put_float (DELTAX,1)
CALL par_put_float (DELTAZ,1)
CALL par_put_float (DZRAT,1)
CALL par_put_float (DZMAX,1)
CALL par_put_float (ZZ,NZPMAX)
CALL par_put_int   (NACOUST,1)
CALL par_put_int   (IDELTAT,1)
CALL par_put_int   (NSTRATX,MAXGRDS)
CALL par_put_int   (NSTRATY,MAXGRDS)
CALL par_put_int   (NESTZ,1)
CALL par_put_int   (NSTRATZ,NZPMAX)
CALL par_put_float (POLELAT,1)
CALL par_put_float (POLELON,1)
CALL par_put_int   (NINEST,MAXGRDS)
CALL par_put_int   (NJNEST,MAXGRDS)
CALL par_put_int   (NKNEST,MAXGRDS)
CALL par_put_float (CENTLAT,MAXGRDS)
CALL par_put_float (CENTLON,MAXGRDS)
CALL par_put_int   (NNSTTOP,MAXGRDS)
CALL par_put_int   (NNSTBOT,MAXGRDS)
CALL par_put_int   (INITIAL,1)
CALL par_put_int   (INORAINTIME,1)
CALL par_put_int   (ITEMPNUDGE,1)
CALL par_put_float (TEMPNUDGEVALS, NZPMAX)
CALL par_put_int   (ITNTS,1)
CALL par_put_int   (NUD_TYPE,1)
CALL par_put_int   (NUDLAT,1)
CALL par_put_float (TNUDLAT,1)
CALL par_put_float (TNUDCENT,1)
CALL par_put_float (TNUDTOP,1)
CALL par_put_float (ZNUDTOP,1)
CALL par_put_float (WT_NUDGE_G,maxgrds)
CALL par_put_float (WT_NUDGE_UV,1)
CALL par_put_float (WT_NUDGE_TH,1)
CALL par_put_float (WT_NUDGE_PI,1)
CALL par_put_float (WT_NUDGE_RT,1)
CALL par_put_int   (NUD_COND,1)
CALL par_put_float (TCOND_BEG,1)
CALL par_put_float (TCOND_END,1)
CALL par_put_float (T_NUDGE_RC,1)
CALL par_put_float (WT_NUDGEC,maxgrds)
CALL par_put_int   (IF_ODA,1)
CALL par_put_int   (IOUTPUT,1)
CALL par_put_float (FRQSTATE,MAXGRDS)
CALL par_put_float (FRQLITE,1)
CALL par_put_float (NLITE_VARS,1)
do nm = 1, nlite_vars
   print*,'lite pack:',nm,trim(LITE_VARS(nm))
   CALL par_put_char (LITE_VARS(nm),32)
enddo
CALL par_put_float (AVGTIM,1)
CALL par_put_float (FRQMEAN,1)
CALL par_put_float (FRQBOTH,1)
CALL par_put_int   (IUPDNDVI,1)
CALL par_put_int   (IUPDSST,1)
CALL par_put_int   (ICORFLG,1)
CALL par_put_int   (IBND,1)
CALL par_put_int   (JBND,1)
CALL par_put_float (CPHAS,1)
CALL par_put_int   (LSFLG,1)
CALL par_put_int   (NFPT,1)
CALL par_put_float (DISTIM,1)
CALL par_put_int   (ILWRTYP,1)
CALL par_put_int   (ISWRTYP,1)
CALL par_put_float (RADFRQ,1)
CALL par_put_int   (LONRAD,1)
CALL par_put_int   (NNQPARM,MAXGRDS)
CALL par_put_float (CONFRQ,1)
CALL par_put_float (WCLDBS,1)
CALL par_put_int   (NPATCH,1)
CALL par_put_int   (NVEGPAT,1)
CALL par_put_int   (ISFCL,1)
CALL par_put_int   (ISOILDAT,1)
CALL par_put_int   (ISNOWDAT,1)
CALL par_put_int   (NVGCON,1)
CALL par_put_float (PCTLCON,1)
CALL par_put_int   (NSLCON,1)
CALL par_put_float (ZROUGH,1)
CALL par_put_float (ALBEDO,1)
CALL par_put_float (SEATMP,1)
CALL par_put_float (DTHCON,1)
CALL par_put_float (DRTCON,1)
CALL par_put_float (SLZ,NZGMAX)
CALL par_put_float (SLMSTR,NZGMAX)
CALL par_put_float (STGOFF,NZGMAX)
CALL par_put_float (CO2_INIT,NZPMAX)
CALL par_put_float (DIVLS,1)
CALL par_put_int   (IDIFFK,MAXGRDS)
CALL par_put_int   (IHORGRAD,1)
CALL par_put_float (CSX,MAXGRDS)
CALL par_put_float (CSZ,MAXGRDS)
CALL par_put_float (XKHKM,MAXGRDS)
CALL par_put_float (ZKHKM,MAXGRDS)
CALL par_put_float (AKMIN,MAXGRDS)
CALL par_put_int   (ICONV,1)
CALL par_put_int   (ICONGR,1)
CALL par_put_int   (ICICENT,1)
CALL par_put_int   (ICJCENT,1)
CALL par_put_float (CXRAD,1)
CALL par_put_float (CYRAD,1)
CALL par_put_int   (ICVERT,1)
CALL par_put_int   (ICKMAX,1)
CALL par_put_float (CZRAD,1)
CALL par_put_int   (ICKCENT,1)
CALL par_put_float (CDIVMAX,1)
CALL par_put_float (CTAU,1)
CALL par_put_float (CTMAX,1)
CALL par_put_int   (ITRACER,1)
CALL par_put_int   (ITRACHIST,1)
CALL par_put_int   (LEVEL,1)
CALL par_put_int   (ICHECKMIC,1)
CALL par_put_int   (IMBUDGET,1)
CALL par_put_int   (IRIME,1)
CALL par_put_int   (IPLAWS,1)
CALL par_put_int   (ISEDIM,1)
CALL par_put_int   (ICLOUD,1)
CALL par_put_int   (IDRIZ,1)
CALL par_put_int   (IRAIN,1)
CALL par_put_int   (IPRIS,1)
CALL par_put_int   (ISNOW,1)
CALL par_put_int   (IAGGR,1)
CALL par_put_int   (IGRAUP,1)
CALL par_put_int   (IHAIL,1)
CALL par_put_float (CPARM,1)
CALL par_put_float (DPARM,1)
CALL par_put_float (RPARM,1)
CALL par_put_float (PPARM,1)
CALL par_put_float (SPARM,1)
CALL par_put_float (APARM,1)
CALL par_put_float (GPARM,1)
CALL par_put_float (HPARM,1)
CALL par_put_float (GNU,8)
CALL par_put_int   (IAEROSOL,1)
CALL par_put_int   (ISALT,1)
CALL par_put_int   (IDUST,1)
CALL par_put_int   (ICCNLEV,1)
CALL par_put_int   (IFORCECCN,1)
CALL par_put_float (BLH,1)
CALL par_put_float (FCCNTS,1)
CALL par_put_float (FCCNSTART,1)
CALL par_put_int   (IIFN,1)
CALL par_put_int   (IAERORAD,1)
CALL par_put_int   (IAERODEP,1)
CALL par_put_int   (IAEROPRNT,1)
CALL par_put_int   (IAEROHIST,1)
CALL par_put_float (CIN_MAX,1)
CALL par_put_float (CCN_MAX,1)
CALL par_put_float (GCCN_MAX,1)
CALL par_put_float (DUST1_MAX,1)
CALL par_put_float (DUST2_MAX,1)
CALL par_put_float (SALTF_MAX,1)
CALL par_put_float (SALTJ_MAX,1)
CALL par_put_float (SALTS_MAX,1)
CALL par_put_int   (IAEROLBC,MAXGRDS)
CALL par_put_int   (ICO2LBC,MAXGRDS)
CALL par_put_float (BCTAU,MAXGRDS)
CALL par_put_int   (IAERO_CHEM,9)
CALL par_put_float (AERO_EPSILON,9)
CALL par_put_float (AERO_MEDRAD,9)
CALL par_put_int   (ITRKEPSILON,1)
CALL par_put_int   (ITRKDUST,1)
CALL par_put_int   (ITRKDUSTIFN,1)
CALL par_put_int   (jnmb,8)
CALL par_put_int   (aero_vanthoff,9)
CALL par_put_float (aero_rhosol,9)
do nm=1,nmachs
   CALL par_send (machnum(nm),2)
enddo

deallocate (buff)

return
END SUBROUTINE masterput_nl

!##############################################################################
Subroutine masterput_gridinit ()

use mem_grid
use rpara

implicit none

real, allocatable :: buff(:)
integer :: nwords,nm

nwords = 1 + 10 * maxgrds
allocate (buff(nwords))

CALL par_init_put (buff,nwords)

CALL par_put_int (NNX,MAXGRDS)
CALL par_put_int (NNX1,MAXGRDS)
CALL par_put_int (NNX2,MAXGRDS)
CALL par_put_int (NNY,MAXGRDS)
CALL par_put_int (NNY1,MAXGRDS)

CALL par_put_int (NNY2,MAXGRDS)
CALL par_put_int (NNZ,MAXGRDS)
CALL par_put_int (NNXYZP,MAXGRDS)
CALL par_put_int (NNXYSP,MAXGRDS)
CALL par_put_int (NNXYP,MAXGRDS)
CALL par_put_int (JDIM,1)

do nm=1,nmachs

   CALL par_send (machnum(nm),12)

enddo

deallocate (buff)

return
END SUBROUTINE masterput_gridinit

!##############################################################################
Subroutine masterput_grid_dimens ()

use mem_grid
use cyclic_mod
use rpara

implicit none

real, allocatable :: buff(:)
integer :: nwords,nm,nxpts,nypts,nzpts

nwords=5*7*maxgrds*maxmach+6*maxgrds*maxmach+2+2*maxmach+42*npts_cyc
allocate (buff(nwords))

CALL par_init_put (buff,nwords)

do nm=1,nmachs
   do ngrid=1,ngrids
      nxpts=nxend(nm,ngrid)-nxbeg(nm,ngrid)+1
      nypts=nyend(nm,ngrid)-nybeg(nm,ngrid)+1
      nzpts=nnzp(ngrid)
      CALL par_put_int (nxpts,1)
      CALL par_put_int (nypts,1)
      CALL par_put_int (nzpts,1)
      CALL par_put_int (nxbegc(nm,ngrid),1)
      CALL par_put_int (nxendc(nm,ngrid),1)
      CALL par_put_int (nybegc(nm,ngrid),1)
      CALL par_put_int (nyendc(nm,ngrid),1)
      CALL par_put_int (ixoff(nm,ngrid),1)
      CALL par_put_int (iyoff(nm,ngrid),1)
      CALL par_put_int (ibcflg(nm,ngrid),1)
   enddo
   CALL par_put_int (machnum(nm),1)
enddo

do nm=1,nmachs
   CALL par_send (machnum(nm),23)
enddo

do nm=1,nmachs
   CALL par_init_put (buff,nwords)
   CALL par_put_int (inode_paths_master(1,1,1,1,nm),5*7*maxgrds*maxmach)
   CALL par_put_int (iget_paths_master(1,1,1,nm),6*maxgrds*maxmach)
   if (npts_cyc > 0) then
      CALL par_put_int (ipathst_cyc,14*npts_cyc)
      CALL par_put_int (ipathsu_cyc,14*npts_cyc)
      CALL par_put_int (ipathsv_cyc,14*npts_cyc)
   endif
   CALL par_put_int (lbc_buffs(1,1,nm),2*maxmach)
   CALL par_put_int (newbuff_nest1(nm),1)
   CALL par_put_int (nbuff_nest1(nm),1)
   CALL par_send (machnum(nm),24)
enddo

deallocate (buff)

return
END SUBROUTINE masterput_grid_dimens

!##############################################################################
Subroutine masterput_gridset ()

use mem_grid
use rpara

implicit none

real, allocatable :: buff(:)
integer :: nwords,nm

nwords=1+maxgrds*(2+8*nzpmax+3*(nxpmax+nypmax))
allocate (buff(nwords))

CALL par_init_put (buff,nwords)

CALL par_put_int   (nrz,nzpmax*maxgrds)
CALL par_put_int   (ipm,nxpmax*maxgrds)
CALL par_put_int   (jpm,nypmax*maxgrds)
CALL par_put_int   (kpm,nzpmax*maxgrds)
CALL par_put_float (xmn,nxpmax*maxgrds)
CALL par_put_float (ymn,nypmax*maxgrds)
CALL par_put_float (zmn,nzpmax*maxgrds)
CALL par_put_float (xtn,nxpmax*maxgrds)
CALL par_put_float (ytn,nypmax*maxgrds)
CALL par_put_float (ztn,nzpmax*maxgrds)
CALL par_put_float (dzmn,nzpmax*maxgrds)
CALL par_put_float (dzm2n,nzpmax*maxgrds)
CALL par_put_float (dztn,nzpmax*maxgrds)
CALL par_put_float (dzt2n,nzpmax*maxgrds)

CALL par_put_float (deltaxn,maxgrds)
CALL par_put_float (deltazn,maxgrds)

CALL par_put_float (ztop,1)

do nm=1,nmachs
   CALL par_send (machnum(nm),35)
enddo

deallocate (buff)

return
END SUBROUTINE masterput_gridset

!##############################################################################
Subroutine masterput_cofnest ()

use mem_grid
use rpara

implicit none

real, allocatable :: buff(:)
integer :: nwords,nm


nwords=maxgrds*(7*(nxpmax+nypmax)+11*nzpmax)
allocate (buff(nwords))

CALL par_init_put (buff,nwords)

CALL par_put_float (ei1,nxpmax*maxgrds)
CALL par_put_float (ei2,nxpmax*maxgrds)
CALL par_put_float (ei3,nxpmax*maxgrds)
CALL par_put_float (ei4,nxpmax*maxgrds)
CALL par_put_float (ei5,nxpmax*maxgrds)
CALL par_put_float (ei6,nxpmax*maxgrds)
CALL par_put_float (ei7,nxpmax*maxgrds)

CALL par_put_float (ej1,nypmax*maxgrds)
CALL par_put_float (ej2,nypmax*maxgrds)
CALL par_put_float (ej3,nypmax*maxgrds)
CALL par_put_float (ej4,nypmax*maxgrds)
CALL par_put_float (ej5,nypmax*maxgrds)
CALL par_put_float (ej6,nypmax*maxgrds)
CALL par_put_float (ej7,nypmax*maxgrds)

CALL par_put_float (ek1,nzpmax*maxgrds)
CALL par_put_float (ek2,nzpmax*maxgrds)
CALL par_put_float (ek3,nzpmax*maxgrds)
CALL par_put_float (ek4,nzpmax*maxgrds)
CALL par_put_float (ek5,nzpmax*maxgrds)
CALL par_put_float (ek6,nzpmax*maxgrds)
CALL par_put_float (ek7,nzpmax*maxgrds)
CALL par_put_float (fbcf,nzpmax*maxgrds*4)

do nm=1,nmachs
   CALL par_send (machnum(nm),36)
enddo

deallocate (buff)

return
END SUBROUTINE masterput_cofnest

!##############################################################################
Subroutine masterput_misc ()

use mem_grid
use rpara
use mem_cuparm
use ref_sounding

implicit none

integer :: nwords, nm
real, allocatable :: buff(:)

nwords=4+13*nzpmax*maxgrds+6*nzpmax
allocate (buff(nwords))

CALL par_init_put (buff,nwords)

CALL par_put_int   (nsubs,1)
CALL par_put_int   (itopo,1)
CALL par_put_int   (impl,1)
CALL par_put_float (time,1)

CALL par_put_float (wsub,nzpmax*maxgrds)
CALL par_put_float (u01dn,nzpmax*maxgrds)
CALL par_put_float (v01dn,nzpmax*maxgrds)
CALL par_put_float (pi01dn,nzpmax*maxgrds)
CALL par_put_float (th01dn,nzpmax*maxgrds)
CALL par_put_float (dn01dn,nzpmax*maxgrds)
CALL par_put_float (rt01dn,nzpmax*maxgrds)

CALL par_put_float (htn,nzpmax*maxgrds)
CALL par_put_float (hwn,nzpmax*maxgrds)
CALL par_put_float (ht2n,nzpmax*maxgrds)
CALL par_put_float (ht4n,nzpmax*maxgrds)
CALL par_put_float (hw2n,nzpmax*maxgrds)
CALL par_put_float (hw4n,nzpmax*maxgrds)
CALL par_put_float (ht,nzpmax)
CALL par_put_float (hw,nzpmax)
CALL par_put_float (ht2,nzpmax)
CALL par_put_float (ht4,nzpmax)
CALL par_put_float (hw2,nzpmax)
CALL par_put_float (hw4,nzpmax)

do nm=1,nmachs
   CALL par_send (machnum(nm),37)
enddo

deallocate (buff)

return
END SUBROUTINE masterput_misc

!##############################################################################
Subroutine masterput_micphys ()

use micphys
use rpara

implicit none

real, allocatable :: buff(:)
integer :: nwords,nm

nwords=2*ncat+nhcat*(5+2*nhcat)+nembc*nembc*(npairc+npairr)
allocate (buff(nwords))

CALL par_init_put (buff,nwords)

CALL par_put_float (shape,nhcat)
CALL par_put_float (cfmas,nhcat)
CALL par_put_float (pwmas,nhcat)
CALL par_put_float (cfvt,nhcat)
CALL par_put_float (pwvt,nhcat)
CALL par_put_float (emb0,ncat)
CALL par_put_float (emb1,ncat)
CALL par_put_float (coltabc,nembc*nembc*npairc)
CALL par_put_float (coltabr,nembc*nembc*npairr)

CALL par_put_int (ipairc,nhcat*nhcat)
CALL par_put_int (ipairr,nhcat*nhcat)


do nm=1,nmachs
   CALL par_send (machnum(nm),41)
enddo

deallocate (buff)

return
END SUBROUTINE masterput_micphys

!##############################################################################
Subroutine nodeget_processid (init)

use grid_dims
use node_mod

implicit none

integer :: init
integer :: ibytes,msgtype,ihostnum
real, allocatable :: buff(:)
integer :: nwords

if(init.eq.1) then
!          get process identifying info from the master process
!          -------------------------------------------------------
   nwords=4+maxmach
   allocate (buff(nwords))
   CALL par_get_new (buff,nwords,1,ibytes,msgtype,ihostnum)
   CALL par_get_int (master_num,1)
   CALL par_get_int (mchnum,1)
   CALL par_get_int (mynum,1)

   CALL par_get_int (nmachs,1)
   CALL par_get_int (machs,nmachs)
   
   deallocate(buff)

endif
!print*,mynum,' ==== got first message',ibytes,msgtype,ihostnum

return
END SUBROUTINE nodeget_processid

!##############################################################################
Subroutine nodeget_nl ()

use mem_all
use node_mod

implicit none

real, allocatable :: buff(:)
integer :: nwords,ibytes,msgtype,ihostnum,nm

!Saleeby(2016)
!Increment memory buffer size here if you add RAMSIN Namelist variables.
!Add to the appropriate section below as (#-of-them * arraysize).
nwords = 133 * 1                 & !single values
       +   2 * 8                 & !micro (8-hydromet types for gnu)
       +   5 * 9                 & !micro (9-aerosol species)
       +  28 * maxgrds    	 & !grid-dependent (max grids)
       +   4 * nzpmax            & !max vertical levels
       +   3 * nzgmax            & !max soil levels
       +   1 * maxlite * 32      & !lite variables 32 char length strings
       + 103                       !extras so we have enough buffer

allocate (buff(nwords))

CALL par_get_new (buff,nwords,2,ibytes,msgtype,ihostnum)

CALL par_get_float (TIMMAX,1)
CALL par_get_int   (IMONTH1,1)
CALL par_get_int   (IDATE1,1)
CALL par_get_int   (IYEAR1,1)
CALL par_get_int   (ITIME1,1)
CALL par_get_int   (NGRIDS,1)
CALL par_get_int   (NNXP,MAXGRDS)
CALL par_get_int   (NNYP,MAXGRDS)
CALL par_get_int   (NNZP,MAXGRDS)
CALL par_get_int   (NZG,1)
CALL par_get_int   (NZS,1)
CALL par_get_int   (NXTNEST,MAXGRDS)
CALL par_get_int   (INESTING,1)
CALL par_get_int   (IPRNTSTMT,1)
CALL par_get_int   (IHTRAN,1)
CALL par_get_float (DELTAX,1)
CALL par_get_float (DELTAZ,1)
CALL par_get_float (DZRAT,1)
CALL par_get_float (DZMAX,1)
CALL par_get_float (ZZ,NZPMAX)
CALL par_get_int   (NACOUST,1)
CALL par_get_int   (IDELTAT,1)
CALL par_get_int   (NSTRATX,MAXGRDS)
CALL par_get_int   (NSTRATY,MAXGRDS)
CALL par_get_int   (NESTZ,1)
CALL par_get_int   (NSTRATZ,NZPMAX)
CALL par_get_float (POLELAT,1)
CALL par_get_float (POLELON,1)
CALL par_get_int   (NINEST,MAXGRDS)
CALL par_get_int   (NJNEST,MAXGRDS)
CALL par_get_int   (NKNEST,MAXGRDS)
CALL par_get_float (CENTLAT,MAXGRDS)
CALL par_get_float (CENTLON,MAXGRDS)
CALL par_get_int   (NNSTTOP,MAXGRDS)
CALL par_get_int   (NNSTBOT,MAXGRDS)
CALL par_get_int   (INITIAL,1)
CALL par_get_int   (INORAINTIME,1)
CALL par_get_int   (ITEMPNUDGE,1)
CALL par_get_float (TEMPNUDGEVALS, NZPMAX)
CALL par_get_int   (ITNTS,1)
CALL par_get_int   (NUD_TYPE,1)
CALL par_get_int   (NUDLAT,1)
CALL par_get_float (TNUDLAT,1)
CALL par_get_float (TNUDCENT,1)
CALL par_get_float (TNUDTOP,1)
CALL par_get_float (ZNUDTOP,1)
CALL par_get_float (WT_NUDGE_G,maxgrds)
CALL par_get_float (WT_NUDGE_UV,1)
CALL par_get_float (WT_NUDGE_TH,1)
CALL par_get_float (WT_NUDGE_PI,1)
CALL par_get_float (WT_NUDGE_RT,1)
CALL par_get_int   (NUD_COND,1)
CALL par_get_float (TCOND_BEG,1)
CALL par_get_float (TCOND_END,1)
CALL par_get_float (T_NUDGE_RC,1)
CALL par_get_float (WT_NUDGEC,maxgrds)
CALL par_get_int   (IF_ODA,1)
CALL par_get_int   (IOUTPUT,1)
CALL par_get_float (FRQSTATE,MAXGRDS)
CALL par_get_float (FRQLITE,1)
CALL par_get_float (NLITE_VARS,1)
do nm = 1, nlite_vars
   CALL par_get_char (LITE_VARS(nm),32)
enddo
CALL par_get_float (AVGTIM,1)
CALL par_get_float (FRQMEAN,1)
CALL par_get_float (FRQBOTH,1)
CALL par_get_int   (IUPDNDVI,1)
CALL par_get_int   (IUPDSST,1)
CALL par_get_int   (ICORFLG,1)
CALL par_get_int   (IBND,1)
CALL par_get_int   (JBND,1)
CALL par_get_float (CPHAS,1)
CALL par_get_int   (LSFLG,1)
CALL par_get_int   (NFPT,1)
CALL par_get_float (DISTIM,1)
CALL par_get_int   (ILWRTYP,1)
CALL par_get_int   (ISWRTYP,1)
CALL par_get_float (RADFRQ,1)
CALL par_get_int   (LONRAD,1)
CALL par_get_int   (NNQPARM,MAXGRDS)
CALL par_get_float (CONFRQ,1)
CALL par_get_float (WCLDBS,1)
CALL par_get_int   (NPATCH,1)
CALL par_get_int   (NVEGPAT,1)
CALL par_get_int   (ISFCL,1)
CALL par_get_int   (ISOILDAT,1)
CALL par_get_int   (ISNOWDAT,1)
CALL par_get_int   (NVGCON,1)
CALL par_get_float (PCTLCON,1)
CALL par_get_int   (NSLCON,1)
CALL par_get_float (ZROUGH,1)
CALL par_get_float (ALBEDO,1)
CALL par_get_float (SEATMP,1)
CALL par_get_float (DTHCON,1)
CALL par_get_float (DRTCON,1)
CALL par_get_float (SLZ,NZGMAX)
CALL par_get_float (SLMSTR,NZGMAX)
CALL par_get_float (STGOFF,NZGMAX)
CALL par_get_float (CO2_INIT,NZPMAX)
CALL par_get_float (DIVLS,1)
CALL par_get_int   (IDIFFK,MAXGRDS)
CALL par_get_int   (IHORGRAD,1)
CALL par_get_float (CSX,MAXGRDS)
CALL par_get_float (CSZ,MAXGRDS)
CALL par_get_float (XKHKM,MAXGRDS)
CALL par_get_float (ZKHKM,MAXGRDS)
CALL par_get_float (AKMIN,MAXGRDS)
CALL par_get_int   (ICONV,1)
CALL par_get_int   (ICONGR,1)
CALL par_get_int   (ICICENT,1)
CALL par_get_int   (ICJCENT,1)
CALL par_get_float (CXRAD,1)
CALL par_get_float (CYRAD,1)
CALL par_get_int   (ICVERT,1)
CALL par_get_int   (ICKMAX,1)
CALL par_get_float (CZRAD,1)
CALL par_get_int   (ICKCENT,1)
CALL par_get_float (CDIVMAX,1)
CALL par_get_float (CTAU,1)
CALL par_get_float (CTMAX,1)
CALL par_get_int   (ITRACER,1)
CALL par_get_int   (ITRACHIST,1)
CALL par_get_int   (LEVEL,1)
CALL par_get_int   (ICHECKMIC,1)
CALL par_get_int   (IMBUDGET,1)
CALL par_get_int   (IRIME,1)
CALL par_get_int   (IPLAWS,1)
CALL par_get_int   (ISEDIM,1)
CALL par_get_int   (ICLOUD,1)
CALL par_get_int   (IDRIZ,1)
CALL par_get_int   (IRAIN,1)
CALL par_get_int   (IPRIS,1)
CALL par_get_int   (ISNOW,1)
CALL par_get_int   (IAGGR,1)
CALL par_get_int   (IGRAUP,1)
CALL par_get_int   (IHAIL,1)
CALL par_get_float (CPARM,1)
CALL par_get_float (DPARM,1)
CALL par_get_float (RPARM,1)
CALL par_get_float (PPARM,1)
CALL par_get_float (SPARM,1)
CALL par_get_float (APARM,1)
CALL par_get_float (GPARM,1)
CALL par_get_float (HPARM,1)
CALL par_get_float (GNU,8)
CALL par_get_int   (IAEROSOL,1)
CALL par_get_int   (ISALT,1)
CALL par_get_int   (IDUST,1)
CALL par_get_int   (ICCNLEV,1)
CALL par_get_int   (IFORCECCN,1)
CALL par_get_float (BLH,1)
CALL par_get_float (FCCNTS,1)
CALL par_get_float (FCCNSTART,1)
CALL par_get_int   (IIFN,1)
CALL par_get_int   (IAERORAD,1)
CALL par_get_int   (IAERODEP,1)
CALL par_get_int   (IAEROPRNT,1)
CALL par_get_int   (IAEROHIST,1)
CALL par_get_float (CIN_MAX,1)
CALL par_get_float (CCN_MAX,1)
CALL par_get_float (GCCN_MAX,1)
CALL par_get_float (DUST1_MAX,1)
CALL par_get_float (DUST2_MAX,1)
CALL par_get_float (SALTF_MAX,1)
CALL par_get_float (SALTJ_MAX,1)
CALL par_get_float (SALTS_MAX,1)
CALL par_get_int   (IAEROLBC,MAXGRDS)
CALL par_get_int   (ICO2LBC,MAXGRDS)
CALL par_get_float (BCTAU,MAXGRDS)
CALL par_get_int   (IAERO_CHEM,9)
CALL par_get_float (AERO_EPSILON,9)
CALL par_get_float (AERO_MEDRAD,9)
CALL par_get_int   (ITRKEPSILON,1)
CALL par_get_int   (ITRKDUST,1)
CALL par_get_int   (ITRKDUSTIFN,1)
CALL par_get_int   (jnmb,8)
CALL par_get_int   (aero_vanthoff,9)
CALL par_get_float (aero_rhosol,9)

deallocate (buff)

return
END SUBROUTINE nodeget_nl

!##############################################################################
Subroutine nodeget_gridinit ()

use mem_grid
use node_mod

implicit none

real, allocatable :: buff(:)
integer :: nwords,ibytes,msgtype,ihostnum

nwords = 1 + 10 * maxgrds
allocate (buff(nwords))

CALL par_get_new (buff,nwords,12,ibytes,msgtype,ihostnum)

CALL par_get_int (NNX,MAXGRDS)
CALL par_get_int (NNX1,MAXGRDS)
CALL par_get_int (NNX2,MAXGRDS)
CALL par_get_int (NNY,MAXGRDS)
CALL par_get_int (NNY1,MAXGRDS)
CALL par_get_int (NNY2,MAXGRDS)
CALL par_get_int (NNZ,MAXGRDS)
CALL par_get_int (NNXYZP,MAXGRDS)
CALL par_get_int (NNXYSP,MAXGRDS)
CALL par_get_int (NNXYP,MAXGRDS)
CALL par_get_int (JDIM,1)

deallocate (buff)

return
END SUBROUTINE nodeget_gridinit

!##############################################################################
Subroutine nodeget_grid_dimens ()

use mem_grid
use node_mod
use cyclic_mod 

implicit none

real, allocatable :: buff(:)
integer :: nwords,ibytes,msgtype,ihostnum,nm,ng

!        get gridpoint information from the master process

nwords=5*7*maxgrds*maxmach+6*maxgrds*maxmach+2+2*maxmach+42*npts_cyc
allocate (buff(nwords))
CALL par_get_new (buff,nwords,23,ibytes,msgtype,ihostnum)

do nm=1,nmachs
   do ng=1,ngrids
      CALL par_get_int (nodemxp(nm,ng),1)
      CALL par_get_int (nodemyp(nm,ng),1)
      CALL par_get_int (nodemzp(nm,ng),1)
      CALL par_get_int (nodeia(nm,ng),1)
      CALL par_get_int (nodeiz(nm,ng),1)
      CALL par_get_int (nodeja(nm,ng),1)
      CALL par_get_int (nodejz(nm,ng),1)
      CALL par_get_int (nodei0(nm,ng),1)
      CALL par_get_int (nodej0(nm,ng),1)
      CALL par_get_int (nodeibcon(nm,ng),1)
   enddo
   CALL par_get_int (machs(nm),1)
enddo

CALL par_get_new (buff,nwords,24,ibytes,msgtype,ihostnum)

CALL par_get_int (ipaths,5*7*maxgrds*maxmach)
CALL par_get_int (iget_paths,6*maxgrds*maxmach)
if(npts_cyc > 0) then
   CALL par_get_int (ipathst_cyc,14*npts_cyc)
   CALL par_get_int (ipathsu_cyc,14*npts_cyc)
   CALL par_get_int (ipathsv_cyc,14*npts_cyc)
endif

do nm=1,maxmach
   CALL par_get_int (node_buffs(nm)%nsend,1) 
   CALL par_get_int (node_buffs(nm)%nrecv,1)
enddo

CALL par_get_int (newbuff_nest,1)
CALL par_get_int (nbuff_nest,1)

deallocate (buff)

do ng=1,ngrids
   mmxp(ng)=nodemxp(mynum,ng)
   mmyp(ng)=nodemyp(mynum,ng)
   mmzp(ng)=nodemzp(mynum,ng)
   mia(ng)=nodeia(mynum,ng)
   miz(ng)=nodeiz(mynum,ng)
   mja(ng)=nodeja(mynum,ng)
   mjz(ng)=nodejz(mynum,ng)
   mi0(ng)=nodei0(mynum,ng)
   mj0(ng)=nodej0(mynum,ng)
   mibcon(ng)=nodeibcon(mynum,ng)
enddo

return
END SUBROUTINE nodeget_grid_dimens

!##############################################################################
Subroutine nodeget_gridset ()

use mem_grid
use node_mod

implicit none

real, allocatable :: buff(:)
integer :: nwords,ibytes,msgtype,ihostnum

nwords=1+maxgrds*(2+8*nzpmax+3*(nxpmax+nypmax))
allocate (buff(nwords))

CALL par_get_new (buff,nwords,35,ibytes,msgtype,ihostnum)

CALL par_get_int   (nrz,nzpmax*maxgrds)
CALL par_get_int   (ipm,nxpmax*maxgrds)
CALL par_get_int   (jpm,nypmax*maxgrds)
CALL par_get_int   (kpm,nzpmax*maxgrds)
CALL par_get_float (xmn,nxpmax*maxgrds)
CALL par_get_float (ymn,nypmax*maxgrds)
CALL par_get_float (zmn,nzpmax*maxgrds)
CALL par_get_float (xtn,nxpmax*maxgrds)
CALL par_get_float (ytn,nypmax*maxgrds)
CALL par_get_float (ztn,nzpmax*maxgrds)
CALL par_get_float (dzmn,nzpmax*maxgrds)
CALL par_get_float (dzm2n,nzpmax*maxgrds)
CALL par_get_float (dztn,nzpmax*maxgrds)
CALL par_get_float (dzt2n,nzpmax*maxgrds)

CALL par_get_float (deltaxn,maxgrds)
CALL par_get_float (deltazn,maxgrds)

CALL par_get_float (ztop,1)

deallocate (buff)

return
END SUBROUTINE nodeget_gridset

!##############################################################################
Subroutine nodeget_cofnest ()

use mem_grid

implicit none

real, allocatable :: buff(:)
integer :: nwords,ibytes,msgtype,ihostnum

nwords=maxgrds*(7*(nxpmax+nypmax)+11*nzpmax)
allocate (buff(nwords))

CALL par_get_new (buff,nwords,36,ibytes,msgtype,ihostnum)

CALL par_get_float (ei1,nxpmax*maxgrds)
CALL par_get_float (ei2,nxpmax*maxgrds)
CALL par_get_float (ei3,nxpmax*maxgrds)
CALL par_get_float (ei4,nxpmax*maxgrds)
CALL par_get_float (ei5,nxpmax*maxgrds)
CALL par_get_float (ei6,nxpmax*maxgrds)
CALL par_get_float (ei7,nxpmax*maxgrds)

CALL par_get_float (ej1,nypmax*maxgrds)
CALL par_get_float (ej2,nypmax*maxgrds)
CALL par_get_float (ej3,nypmax*maxgrds)
CALL par_get_float (ej4,nypmax*maxgrds)
CALL par_get_float (ej5,nypmax*maxgrds)
CALL par_get_float (ej6,nypmax*maxgrds)
CALL par_get_float (ej7,nypmax*maxgrds)

CALL par_get_float (ek1,nzpmax*maxgrds)
CALL par_get_float (ek2,nzpmax*maxgrds)
CALL par_get_float (ek3,nzpmax*maxgrds)
CALL par_get_float (ek4,nzpmax*maxgrds)
CALL par_get_float (ek5,nzpmax*maxgrds)
CALL par_get_float (ek6,nzpmax*maxgrds)
CALL par_get_float (ek7,nzpmax*maxgrds)
CALL par_get_float (fbcf,nzpmax*maxgrds*4)

deallocate (buff)

return
END SUBROUTINE nodeget_cofnest

!##############################################################################
Subroutine nodeget_misc ()

use mem_grid
use mem_cuparm
use ref_sounding

implicit none

integer :: nwords,ibytes,msgtype,ihostnum
real, allocatable :: buff(:)

nwords=4 + 13*nzpmax*maxgrds + 6*nzpmax
allocate (buff(nwords))

CALL par_get_new (buff,nwords,37,ibytes,msgtype,ihostnum)

CALL par_get_int   (nsubs,1)
CALL par_get_int   (itopo,1)
CALL par_get_int   (impl,1)
CALL par_get_float (time,1)

CALL par_get_float (wsub,nzpmax*maxgrds)
CALL par_get_float (u01dn,nzpmax*maxgrds)
CALL par_get_float (v01dn,nzpmax*maxgrds)
CALL par_get_float (pi01dn,nzpmax*maxgrds)
CALL par_get_float (th01dn,nzpmax*maxgrds)
CALL par_get_float (dn01dn,nzpmax*maxgrds)
CALL par_get_float (rt01dn,nzpmax*maxgrds)

CALL par_get_float (htn,nzpmax*maxgrds)
CALL par_get_float (hwn,nzpmax*maxgrds)
CALL par_get_float (ht2n,nzpmax*maxgrds)
CALL par_get_float (ht4n,nzpmax*maxgrds)
CALL par_get_float (hw2n,nzpmax*maxgrds)
CALL par_get_float (hw4n,nzpmax*maxgrds)
CALL par_get_float (ht,nzpmax)
CALL par_get_float (hw,nzpmax)
CALL par_get_float (ht2,nzpmax)
CALL par_get_float (ht4,nzpmax)
CALL par_get_float (hw2,nzpmax)
CALL par_get_float (hw4,nzpmax)

deallocate (buff)

return
END SUBROUTINE nodeget_misc

!##############################################################################
Subroutine nodeget_micphys ()

use micphys

implicit none

real, allocatable :: buff(:)
integer :: nwords,ibytes,msgtype,ihostnum

nwords=2*ncat+nhcat*(5+2*nhcat)+nembc*nembc*(npairc+npairr)
allocate (buff(nwords))

CALL par_get_new (buff,nwords,41,ibytes,msgtype,ihostnum)

CALL par_get_float (shape,nhcat)
CALL par_get_float (cfmas,nhcat)
CALL par_get_float (pwmas,nhcat)
CALL par_get_float (cfvt,nhcat)
CALL par_get_float (pwvt,nhcat)

CALL par_get_float (emb0,ncat)
CALL par_get_float (emb1,ncat)
CALL par_get_float (coltabc,nembc*nembc*npairc)
CALL par_get_float (coltabr,nembc*nembc*npairr)
CALL par_get_int   (ipairc,nhcat*nhcat)
CALL par_get_int   (ipairr,nhcat*nhcat)

deallocate (buff)

return
END SUBROUTINE nodeget_micphys
