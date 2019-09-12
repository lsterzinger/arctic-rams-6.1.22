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

Module mem_tend

implicit none

   Type tend_vars

   real, allocatable, dimension(:) :: &
         ut, vt, wt, pt, tht, rtt  &
        ,rct, rdt, rrt, rpt, rst, rat, rgt ,rht      &
        ,cct, cdt, crt, cpt, cst, cat, cgt ,cht      &
        ,cifnt, tket, epst                           &
        ,cccmt, gccmt, cccnt, gccnt                  &
        ,md1nt, md2nt, md1mt, md2mt                  &
        ,salt_film_nt,salt_jet_nt,salt_spum_nt       &
        ,salt_film_mt,salt_jet_mt,salt_spum_mt       &
        ,regen_aero1_nt,regen_aero1_mt               &
        ,regen_aero2_nt,regen_aero2_mt               &
        ,immerct, immerdt, immerrt ,ifnnuct          &
        ,cnmct, cnmdt, cnmrt, cnmpt, cnmst           &
        ,cnmat, cnmgt, cnmht                         &
        ,dnmct, dnmdt, dnmrt, dnmpt, dnmst           &
        ,dnmat, dnmgt, dnmht                         &
        ,dinct, dindt, dinrt, dinpt, dinst           &
        ,dinat, dingt, dinht                         &
        ,snmct, snmdt, snmrt, snmpt, snmst           &
        ,snmat, snmgt, snmht                         &
        ,resol_aero1_mt,resol_aero2_mt               &
        ,rco2t

   End Type
   
   type (tend_vars) :: tend

Contains

!##############################################################################
Subroutine alloc_tend (nmzp,nmxp,nmyp,ngrs,proc_type)

use mem_basic
use mem_tracer
use mem_sib
use mem_micro
use mem_turb

implicit none
   
   integer, dimension (*) :: nmzp,nmxp,nmyp
   integer, intent(in) :: ngrs,proc_type
   integer :: ng,ntpts,nsc

!print*, 'enter alloc_tend'

!         Find the maximum number of grid points needed for any grid.

   if(proc_type==1) then
      ntpts=1
   else
      ntpts=0
      do ng=1,ngrs
         ntpts=max( nmxp(ng)*nmyp(ng)*nmzp(ng),ntpts )
      enddo
   endif

! Allocate arrays based on options (if necessary).
! Do not need these arrays if it is a master process in a parallel run,
! so just allocate to 1 word.

!!!!!  WE ARE ONLY CHECKING GRID 1 !!!!!!!!!
!!!!!  All grids must have same scalars defined !!!!!!!

   if (allocated(basic_g(1)%up))      allocate (tend%ut(ntpts))
   if (allocated(basic_g(1)%vp))      allocate (tend%vt(ntpts))
   if (allocated(basic_g(1)%wp))      allocate (tend%wt(ntpts))
   if (allocated(basic_g(1)%pp))      allocate (tend%pt(ntpts))
      
   if (allocated(basic_g(1)%thp))     allocate (tend%tht(ntpts))
   if (allocated(basic_g(1)%rtp))     allocate (tend%rtt(ntpts))
   if (allocated(micro_g(1)%rcp))     allocate (tend%rct(ntpts))
   if (allocated(micro_g(1)%rdp))     allocate (tend%rdt(ntpts))
   if (allocated(micro_g(1)%rrp))     allocate (tend%rrt(ntpts))
   if (allocated(micro_g(1)%rpp))     allocate (tend%rpt(ntpts))
   if (allocated(micro_g(1)%rsp))     allocate (tend%rst(ntpts))
   if (allocated(micro_g(1)%rap))     allocate (tend%rat(ntpts))
   if (allocated(micro_g(1)%rgp))     allocate (tend%rgt(ntpts))
   if (allocated(micro_g(1)%rhp))     allocate (tend%rht(ntpts))
   if (allocated(micro_g(1)%ccp))     allocate (tend%cct(ntpts))
   if (allocated(micro_g(1)%cdp))     allocate (tend%cdt(ntpts))
   if (allocated(micro_g(1)%crp))     allocate (tend%crt(ntpts))
   if (allocated(micro_g(1)%cpp))     allocate (tend%cpt(ntpts))
   if (allocated(micro_g(1)%csp))     allocate (tend%cst(ntpts))
   if (allocated(micro_g(1)%cap))     allocate (tend%cat(ntpts))
   if (allocated(micro_g(1)%cgp))     allocate (tend%cgt(ntpts))
   if (allocated(micro_g(1)%chp))     allocate (tend%cht(ntpts))
   if (allocated(micro_g(1)%cifnp))   allocate (tend%cifnt(ntpts))

   if (allocated(turb_g(1)%tkep))     allocate (tend%tket(ntpts))
   if (allocated(turb_g(1)%epsp))     allocate (tend%epst(ntpts))

   if (allocated(micro_g(1)%cccnp))   allocate (tend%cccnt(ntpts))
   if (allocated(micro_g(1)%gccnp))   allocate (tend%gccnt(ntpts))
   if (allocated(micro_g(1)%cccmp))   allocate (tend%cccmt(ntpts))
   if (allocated(micro_g(1)%gccmp))   allocate (tend%gccmt(ntpts))
   if (allocated(micro_g(1)%md1np))   allocate (tend%md1nt(ntpts))
   if (allocated(micro_g(1)%md2np))   allocate (tend%md2nt(ntpts))
   if (allocated(micro_g(1)%md1mp))   allocate (tend%md1mt(ntpts))
   if (allocated(micro_g(1)%md2mp))   allocate (tend%md2mt(ntpts))
   if (allocated(micro_g(1)%salt_film_np))  allocate (tend%salt_film_nt(ntpts))
   if (allocated(micro_g(1)%salt_jet_np))   allocate (tend%salt_jet_nt(ntpts))
   if (allocated(micro_g(1)%salt_spum_np))  allocate (tend%salt_spum_nt(ntpts))
   if (allocated(micro_g(1)%salt_film_mp))  allocate (tend%salt_film_mt(ntpts))
   if (allocated(micro_g(1)%salt_jet_mp))   allocate (tend%salt_jet_mt(ntpts))
   if (allocated(micro_g(1)%salt_spum_mp))  allocate (tend%salt_spum_mt(ntpts))
   if (allocated(micro_g(1)%regen_aero1_np)) allocate (tend%regen_aero1_nt(ntpts))
   if (allocated(micro_g(1)%regen_aero1_mp)) allocate (tend%regen_aero1_mt(ntpts))
   if (allocated(micro_g(1)%regen_aero2_np)) allocate (tend%regen_aero2_nt(ntpts))
   if (allocated(micro_g(1)%regen_aero2_mp)) allocate (tend%regen_aero2_mt(ntpts))

   if (allocated(micro_g(1)%immercp)) allocate (tend%immerct(ntpts))
   if (allocated(micro_g(1)%immerdp)) allocate (tend%immerdt(ntpts))
   if (allocated(micro_g(1)%immerrp)) allocate (tend%immerrt(ntpts))
   if (allocated(micro_g(1)%ifnnucp)) allocate (tend%ifnnuct(ntpts))

   if (allocated(micro_g(1)%cnmcp))   allocate (tend%cnmct(ntpts))
   if (allocated(micro_g(1)%cnmdp))   allocate (tend%cnmdt(ntpts))
   if (allocated(micro_g(1)%cnmrp))   allocate (tend%cnmrt(ntpts))
   if (allocated(micro_g(1)%cnmpp))   allocate (tend%cnmpt(ntpts))
   if (allocated(micro_g(1)%cnmsp))   allocate (tend%cnmst(ntpts))
   if (allocated(micro_g(1)%cnmap))   allocate (tend%cnmat(ntpts))
   if (allocated(micro_g(1)%cnmgp))   allocate (tend%cnmgt(ntpts))
   if (allocated(micro_g(1)%cnmhp))   allocate (tend%cnmht(ntpts))

   if (allocated(micro_g(1)%dnmcp))   allocate (tend%dnmct(ntpts))
   if (allocated(micro_g(1)%dnmdp))   allocate (tend%dnmdt(ntpts))
   if (allocated(micro_g(1)%dnmrp))   allocate (tend%dnmrt(ntpts))
   if (allocated(micro_g(1)%dnmpp))   allocate (tend%dnmpt(ntpts))
   if (allocated(micro_g(1)%dnmsp))   allocate (tend%dnmst(ntpts))
   if (allocated(micro_g(1)%dnmap))   allocate (tend%dnmat(ntpts))
   if (allocated(micro_g(1)%dnmgp))   allocate (tend%dnmgt(ntpts))
   if (allocated(micro_g(1)%dnmhp))   allocate (tend%dnmht(ntpts))

   if (allocated(micro_g(1)%dincp))   allocate (tend%dinct(ntpts))
   if (allocated(micro_g(1)%dindp))   allocate (tend%dindt(ntpts))
   if (allocated(micro_g(1)%dinrp))   allocate (tend%dinrt(ntpts))
   if (allocated(micro_g(1)%dinpp))   allocate (tend%dinpt(ntpts))
   if (allocated(micro_g(1)%dinsp))   allocate (tend%dinst(ntpts))
   if (allocated(micro_g(1)%dinap))   allocate (tend%dinat(ntpts))
   if (allocated(micro_g(1)%dingp))   allocate (tend%dingt(ntpts))
   if (allocated(micro_g(1)%dinhp))   allocate (tend%dinht(ntpts))

   if (allocated(micro_g(1)%snmcp))   allocate (tend%snmct(ntpts))
   if (allocated(micro_g(1)%snmdp))   allocate (tend%snmdt(ntpts))
   if (allocated(micro_g(1)%snmrp))   allocate (tend%snmrt(ntpts))
   if (allocated(micro_g(1)%snmpp))   allocate (tend%snmpt(ntpts))
   if (allocated(micro_g(1)%snmsp))   allocate (tend%snmst(ntpts))
   if (allocated(micro_g(1)%snmap))   allocate (tend%snmat(ntpts))
   if (allocated(micro_g(1)%snmgp))   allocate (tend%snmgt(ntpts))
   if (allocated(micro_g(1)%snmhp))   allocate (tend%snmht(ntpts))
   if (allocated(micro_g(1)%resol_aero1_mp)) allocate (tend%resol_aero1_mt(ntpts))
   if (allocated(micro_g(1)%resol_aero2_mp)) allocate (tend%resol_aero2_mt(ntpts))

   if (allocated(sib_g(1)%rco2p))     allocate (tend%rco2t(ntpts))

   do ng=1,ngrs
    do nsc=1,itracer
      if (allocated(tracer_g(nsc,ng)%tracerp).and.  &
             (.not.allocated(tracer_g(nsc,ng)%tracert)))  &
              allocate (tracer_g(nsc,ng)%tracert(ntpts))
    enddo
   enddo

return
END SUBROUTINE alloc_tend

!##############################################################################
Subroutine dealloc_tend (ngrs)

use mem_tracer

implicit none

integer, intent(in) :: ngrs
integer :: nsc,ng

! Deallocate all tendency arrays

   if (allocated(tend%ut))   deallocate (tend%ut)
   if (allocated(tend%vt))   deallocate (tend%vt)
   if (allocated(tend%wt))   deallocate (tend%wt)
   if (allocated(tend%pt))   deallocate (tend%pt)
   if (allocated(tend%tht))  deallocate (tend%tht)
   if (allocated(tend%rtt))  deallocate (tend%rtt)
   if (allocated(tend%rct))  deallocate (tend%rct)
   if (allocated(tend%rdt))  deallocate (tend%rdt)
   if (allocated(tend%rrt))  deallocate (tend%rrt)
   if (allocated(tend%rpt))  deallocate (tend%rpt)
   if (allocated(tend%rst))  deallocate (tend%rst)
   if (allocated(tend%rat))  deallocate (tend%rat)
   if (allocated(tend%rgt))  deallocate (tend%rgt)
   if (allocated(tend%rht))  deallocate (tend%rht)
   if (allocated(tend%cct))  deallocate (tend%cct)
   if (allocated(tend%cdt))  deallocate (tend%cdt)
   if (allocated(tend%crt))  deallocate (tend%crt)
   if (allocated(tend%cpt))  deallocate (tend%cpt)
   if (allocated(tend%cst))  deallocate (tend%cst)
   if (allocated(tend%cat))  deallocate (tend%cat)
   if (allocated(tend%cgt))  deallocate (tend%cgt)
   if (allocated(tend%cht))  deallocate (tend%cht)
   if (allocated(tend%cifnt))deallocate (tend%cifnt)

   if (allocated(tend%tket)) deallocate (tend%tket)
   if (allocated(tend%epst)) deallocate (tend%epst)

   if (allocated(tend%cccnt))deallocate (tend%cccnt)
   if (allocated(tend%gccnt))deallocate (tend%gccnt)
   if (allocated(tend%cccmt)) deallocate (tend%cccmt)
   if (allocated(tend%gccmt)) deallocate (tend%gccmt)
   if (allocated(tend%md1nt)) deallocate (tend%md1nt)
   if (allocated(tend%md2nt)) deallocate (tend%md2nt)
   if (allocated(tend%md1mt)) deallocate (tend%md1mt)
   if (allocated(tend%md2mt)) deallocate (tend%md2mt)
   if (allocated(tend%salt_film_nt))   deallocate (tend%salt_film_nt)
   if (allocated(tend%salt_jet_nt))    deallocate (tend%salt_jet_nt)
   if (allocated(tend%salt_spum_nt))   deallocate (tend%salt_spum_nt)
   if (allocated(tend%salt_film_mt))   deallocate (tend%salt_film_mt)
   if (allocated(tend%salt_jet_mt))    deallocate (tend%salt_jet_mt)
   if (allocated(tend%salt_spum_mt))   deallocate (tend%salt_spum_mt)
   if (allocated(tend%regen_aero1_nt)) deallocate (tend%regen_aero1_nt)
   if (allocated(tend%regen_aero1_mt)) deallocate (tend%regen_aero1_mt)
   if (allocated(tend%regen_aero2_nt)) deallocate (tend%regen_aero2_nt)
   if (allocated(tend%regen_aero2_mt)) deallocate (tend%regen_aero2_mt)

   if (allocated(tend%immerct))        deallocate (tend%immerct)
   if (allocated(tend%immerdt))        deallocate (tend%immerdt)
   if (allocated(tend%immerrt))        deallocate (tend%immerrt)
   if (allocated(tend%ifnnuct))        deallocate (tend%ifnnuct)

   if (allocated(tend%cnmct))   deallocate (tend%cnmct)
   if (allocated(tend%cnmdt))   deallocate (tend%cnmdt)
   if (allocated(tend%cnmrt))   deallocate (tend%cnmrt)
   if (allocated(tend%cnmpt))   deallocate (tend%cnmpt)
   if (allocated(tend%cnmst))   deallocate (tend%cnmst)
   if (allocated(tend%cnmat))   deallocate (tend%cnmat)
   if (allocated(tend%cnmgt))   deallocate (tend%cnmgt)
   if (allocated(tend%cnmht))   deallocate (tend%cnmht)

   if (allocated(tend%dnmct)) deallocate (tend%dnmct)
   if (allocated(tend%dnmdt)) deallocate (tend%dnmdt)
   if (allocated(tend%dnmrt)) deallocate (tend%dnmrt)
   if (allocated(tend%dnmpt)) deallocate (tend%dnmpt)
   if (allocated(tend%dnmst)) deallocate (tend%dnmst)
   if (allocated(tend%dnmat)) deallocate (tend%dnmat)
   if (allocated(tend%dnmgt)) deallocate (tend%dnmgt)
   if (allocated(tend%dnmht)) deallocate (tend%dnmht)

   if (allocated(tend%dinct)) deallocate (tend%dinct)
   if (allocated(tend%dindt)) deallocate (tend%dindt)
   if (allocated(tend%dinrt)) deallocate (tend%dinrt)
   if (allocated(tend%dinpt)) deallocate (tend%dinpt)
   if (allocated(tend%dinst)) deallocate (tend%dinst)
   if (allocated(tend%dinat)) deallocate (tend%dinat)
   if (allocated(tend%dingt)) deallocate (tend%dingt)
   if (allocated(tend%dinht)) deallocate (tend%dinht)

   if (allocated(tend%snmct)) deallocate (tend%snmct)
   if (allocated(tend%snmdt)) deallocate (tend%snmdt)
   if (allocated(tend%snmrt)) deallocate (tend%snmrt)
   if (allocated(tend%snmpt)) deallocate (tend%snmpt)
   if (allocated(tend%snmst)) deallocate (tend%snmst)
   if (allocated(tend%snmat)) deallocate (tend%snmat)
   if (allocated(tend%snmgt)) deallocate (tend%snmgt)
   if (allocated(tend%snmht)) deallocate (tend%snmht)
   if (allocated(tend%resol_aero1_mt)) deallocate (tend%resol_aero1_mt)
   if (allocated(tend%resol_aero2_mt)) deallocate (tend%resol_aero2_mt)

   if (allocated(tend%rco2t)) deallocate (tend%rco2t)

   do ng=1,ngrs
    do nsc=1,itracer
      if (allocated(tracer_g(nsc,ng)%tracert)) &
        deallocate (tracer_g(nsc,ng)%tracert)
    enddo
   enddo

return
END SUBROUTINE dealloc_tend

!##############################################################################
Subroutine filltab_tend (basic,micro,turb,sib,tracer,ng)

use mem_basic
use mem_micro
use mem_turb
use mem_sib
use mem_tracer
use var_tables

implicit none

   type (basic_vars) :: basic
   type (micro_vars) :: micro
   type (turb_vars)  :: turb
   type (sib_vars)   :: sib
   type (tracer_vars) :: tracer(*)
   integer, intent(in) :: ng
   integer :: nsc
   character (len=10) :: sname

! Fill scalar arrays into scalar tables

   if (allocated(tend%tht))  &
      CALL vtables_scalar (basic%thp(1,1,1),tend%tht(1),ng,'THP')
   if (allocated(tend%rtt))  &
      CALL vtables_scalar (basic%rtp(1,1,1),tend%rtt(1),ng,'RTP')
   if (allocated(tend%rct))  &
      CALL vtables_scalar (micro%rcp(1,1,1),tend%rct(1),ng,'RCP')
   if (allocated(tend%rdt))  &
      CALL vtables_scalar (micro%rdp(1,1,1),tend%rdt(1),ng,'RDP')
   if (allocated(tend%rrt))  &
      CALL vtables_scalar (micro%rrp(1,1,1),tend%rrt(1),ng,'RRP')
   if (allocated(tend%rpt))  &
      CALL vtables_scalar (micro%rpp(1,1,1),tend%rpt(1),ng,'RPP')
   if (allocated(tend%rst))  &
      CALL vtables_scalar (micro%rsp(1,1,1),tend%rst(1),ng,'RSP')
   if (allocated(tend%rat))  &
      CALL vtables_scalar (micro%rap(1,1,1),tend%rat(1),ng,'RAP')
   if (allocated(tend%rgt))  &
      CALL vtables_scalar (micro%rgp(1,1,1),tend%rgt(1),ng,'RGP')
   if (allocated(tend%rht))  &
      CALL vtables_scalar (micro%rhp(1,1,1),tend%rht(1),ng,'RHP')
   if (allocated(tend%cct))  &
      CALL vtables_scalar (micro%ccp(1,1,1),tend%cct(1),ng,'CCP')
   if (allocated(tend%cdt))  &
      CALL vtables_scalar (micro%cdp(1,1,1),tend%cdt(1),ng,'CDP')
   if (allocated(tend%crt))  &
      CALL vtables_scalar (micro%crp(1,1,1),tend%crt(1),ng,'CRP')
   if (allocated(tend%cpt))  &
      CALL vtables_scalar (micro%cpp(1,1,1),tend%cpt(1),ng,'CPP')
   if (allocated(tend%cst))  &
      CALL vtables_scalar (micro%csp(1,1,1),tend%cst(1),ng,'CSP')
   if (allocated(tend%cat))  &
      CALL vtables_scalar (micro%cap(1,1,1),tend%cat(1),ng,'CAP')
   if (allocated(tend%cgt))  &
      CALL vtables_scalar (micro%cgp(1,1,1),tend%cgt(1),ng,'CGP')
   if (allocated(tend%cht))  &
      CALL vtables_scalar (micro%chp(1,1,1),tend%cht(1),ng,'CHP')
   if (allocated(tend%cifnt))  &
      CALL vtables_scalar (micro%cifnp(1,1,1),tend%cifnt(1),ng,'CIFNP')

   if( allocated(tend%tket))  &
      CALL vtables_scalar (turb%tkep(1,1,1),tend%tket(1),ng,'TKEP')
   if( allocated(tend%epst))  &
      CALL vtables_scalar (turb%epsp(1,1,1),tend%epst(1),ng,'EPSP')

   if (allocated(tend%cccnt))  &
      CALL vtables_scalar (micro%cccnp(1,1,1),tend%cccnt(1),ng,'CCCNP')
   if (allocated(tend%gccnt))  &
      CALL vtables_scalar (micro%gccnp(1,1,1),tend%gccnt(1),ng,'GCCNP')
   if (allocated(tend%cccmt))  &
      CALL vtables_scalar (micro%cccmp(1,1,1),tend%cccmt(1),ng,'CCCMP')
   if (allocated(tend%gccmt))  &
      CALL vtables_scalar (micro%gccmp(1,1,1),tend%gccmt(1),ng,'GCCMP')
   if (allocated(tend%md1nt))  &
      CALL vtables_scalar (micro%md1np(1,1,1),tend%md1nt(1),ng,'MD1NP')
   if (allocated(tend%md2nt))  &
      CALL vtables_scalar (micro%md2np(1,1,1),tend%md2nt(1),ng,'MD2NP')
   if (allocated(tend%md1mt))  &
      CALL vtables_scalar (micro%md1mp(1,1,1),tend%md1mt(1),ng,'MD1MP')
   if (allocated(tend%md2mt))  &
      CALL vtables_scalar (micro%md2mp(1,1,1),tend%md2mt(1),ng,'MD2MP')
   if (allocated(tend%salt_film_nt))  &
      CALL vtables_scalar (micro%salt_film_np(1,1,1) &
                          ,tend%salt_film_nt(1),ng,'SALT_FILM_NP')
   if (allocated(tend%salt_jet_nt))  &
      CALL vtables_scalar (micro%salt_jet_np(1,1,1) &
                          ,tend%salt_jet_nt(1),ng,'SALT_JET_NP')
   if (allocated(tend%salt_spum_nt))  &
      CALL vtables_scalar (micro%salt_spum_np(1,1,1) &
                          ,tend%salt_spum_nt(1),ng,'SALT_SPUM_NP')
   if (allocated(tend%salt_film_mt))  &
      CALL vtables_scalar (micro%salt_film_mp(1,1,1) &
                          ,tend%salt_film_mt(1),ng,'SALT_FILM_MP')
   if (allocated(tend%salt_jet_mt))  &
      CALL vtables_scalar (micro%salt_jet_mp(1,1,1) &
                          ,tend%salt_jet_mt(1),ng,'SALT_JET_MP')
   if (allocated(tend%salt_spum_mt))  &
      CALL vtables_scalar (micro%salt_spum_mp(1,1,1) &
                          ,tend%salt_spum_mt(1),ng,'SALT_SPUM_MP')
   if (allocated(tend%regen_aero1_nt))  &
      CALL vtables_scalar (micro%regen_aero1_np(1,1,1) &
                          ,tend%regen_aero1_nt(1),ng,'REGEN_AERO1_NP')
   if (allocated(tend%regen_aero1_mt))  &
      CALL vtables_scalar (micro%regen_aero1_mp(1,1,1) &
                          ,tend%regen_aero1_mt(1),ng,'REGEN_AERO1_MP')
   if (allocated(tend%regen_aero2_nt))  &
      CALL vtables_scalar (micro%regen_aero2_np(1,1,1) &
                          ,tend%regen_aero2_nt(1),ng,'REGEN_AERO2_NP')
   if (allocated(tend%regen_aero2_mt))  &
      CALL vtables_scalar (micro%regen_aero2_mp(1,1,1) &
                          ,tend%regen_aero2_mt(1),ng,'REGEN_AERO2_MP')


   if (allocated(tend%immerct))  &
      CALL vtables_scalar (micro%immercp(1,1,1),tend%immerct(1),ng,'IMMERCP')
   if (allocated(tend%immerdt))  &
      CALL vtables_scalar (micro%immerdp(1,1,1),tend%immerdt(1),ng,'IMMERDP')
   if (allocated(tend%immerrt))  &
      CALL vtables_scalar (micro%immerrp(1,1,1),tend%immerrt(1),ng,'IMMERRP')
   if (allocated(tend%ifnnuct))  &
      CALL vtables_scalar (micro%ifnnucp(1,1,1),tend%ifnnuct(1),ng,'IFNNUCP')

   if (allocated(tend%cnmct))  &
      CALL vtables_scalar (micro%cnmcp(1,1,1),tend%cnmct(1),ng,'CNMCP')
   if (allocated(tend%cnmdt))  &
      CALL vtables_scalar (micro%cnmdp(1,1,1),tend%cnmdt(1),ng,'CNMDP')
   if (allocated(tend%cnmrt))  &
      CALL vtables_scalar (micro%cnmrp(1,1,1),tend%cnmrt(1),ng,'CNMRP')
   if (allocated(tend%cnmpt))  &
      CALL vtables_scalar (micro%cnmpp(1,1,1),tend%cnmpt(1),ng,'CNMPP')
   if (allocated(tend%cnmst))  &
      CALL vtables_scalar (micro%cnmsp(1,1,1),tend%cnmst(1),ng,'CNMSP')
   if (allocated(tend%cnmat))  &
      CALL vtables_scalar (micro%cnmap(1,1,1),tend%cnmat(1),ng,'CNMAP')
   if (allocated(tend%cnmgt))  &
      CALL vtables_scalar (micro%cnmgp(1,1,1),tend%cnmgt(1),ng,'CNMGP')
   if (allocated(tend%cnmht))  &
      CALL vtables_scalar (micro%cnmhp(1,1,1),tend%cnmht(1),ng,'CNMHP')

   if (allocated(tend%dnmct))  &
      CALL vtables_scalar (micro%dnmcp(1,1,1),tend%dnmct(1),ng,'DNMCP')
   if (allocated(tend%dnmdt))  &
      CALL vtables_scalar (micro%dnmdp(1,1,1),tend%dnmdt(1),ng,'DNMDP')
   if (allocated(tend%dnmrt))  &
      CALL vtables_scalar (micro%dnmrp(1,1,1),tend%dnmrt(1),ng,'DNMRP')
   if (allocated(tend%dnmpt))  &
      CALL vtables_scalar (micro%dnmpp(1,1,1),tend%dnmpt(1),ng,'DNMPP')
   if (allocated(tend%dnmst))  &
      CALL vtables_scalar (micro%dnmsp(1,1,1),tend%dnmst(1),ng,'DNMSP')
   if (allocated(tend%dnmat))  &
      CALL vtables_scalar (micro%dnmap(1,1,1),tend%dnmat(1),ng,'DNMAP')
   if (allocated(tend%dnmgt))  &
      CALL vtables_scalar (micro%dnmgp(1,1,1),tend%dnmgt(1),ng,'DNMGP')
   if (allocated(tend%dnmht))  &
      CALL vtables_scalar (micro%dnmhp(1,1,1),tend%dnmht(1),ng,'DNMHP')

   if (allocated(tend%dinct))  &
      CALL vtables_scalar (micro%dincp(1,1,1),tend%dinct(1),ng,'DINCP')
   if (allocated(tend%dindt))  &
      CALL vtables_scalar (micro%dindp(1,1,1),tend%dindt(1),ng,'DINDP')
   if (allocated(tend%dinrt))  &
      CALL vtables_scalar (micro%dinrp(1,1,1),tend%dinrt(1),ng,'DINRP')
   if (allocated(tend%dinpt))  &
      CALL vtables_scalar (micro%dinpp(1,1,1),tend%dinpt(1),ng,'DINPP')
   if (allocated(tend%dinst))  &
      CALL vtables_scalar (micro%dinsp(1,1,1),tend%dinst(1),ng,'DINSP')
   if (allocated(tend%dinat))  &
      CALL vtables_scalar (micro%dinap(1,1,1),tend%dinat(1),ng,'DINAP')
   if (allocated(tend%dingt))  &
      CALL vtables_scalar (micro%dingp(1,1,1),tend%dingt(1),ng,'DINGP')
   if (allocated(tend%dinht))  &
      CALL vtables_scalar (micro%dinhp(1,1,1),tend%dinht(1),ng,'DINHP')

   if (allocated(tend%snmct))  &
      CALL vtables_scalar (micro%snmcp(1,1,1),tend%snmct(1),ng,'SNMCP')
   if (allocated(tend%snmdt))  &
      CALL vtables_scalar (micro%snmdp(1,1,1),tend%snmdt(1),ng,'SNMDP')
   if (allocated(tend%snmrt))  &
      CALL vtables_scalar (micro%snmrp(1,1,1),tend%snmrt(1),ng,'SNMRP')
   if (allocated(tend%snmpt))  &
      CALL vtables_scalar (micro%snmpp(1,1,1),tend%snmpt(1),ng,'SNMPP')
   if (allocated(tend%snmst))  &
      CALL vtables_scalar (micro%snmsp(1,1,1),tend%snmst(1),ng,'SNMSP')
   if (allocated(tend%snmat))  &
      CALL vtables_scalar (micro%snmap(1,1,1),tend%snmat(1),ng,'SNMAP')
   if (allocated(tend%snmgt))  &
      CALL vtables_scalar (micro%snmgp(1,1,1),tend%snmgt(1),ng,'SNMGP')
   if (allocated(tend%snmht))  &
      CALL vtables_scalar (micro%snmhp(1,1,1),tend%snmht(1),ng,'SNMHP')
   if (allocated(tend%resol_aero1_mt))  &
      CALL vtables_scalar (micro%resol_aero1_mp(1,1,1) &
                          ,tend%resol_aero1_mt(1),ng,'RESOL_AERO1_MP')
   if (allocated(tend%resol_aero2_mt))  &
      CALL vtables_scalar (micro%resol_aero2_mp(1,1,1) &
                          ,tend%resol_aero2_mt(1),ng,'RESOL_AERO2_MP')

   if (allocated(tend%rco2t))  &
      CALL vtables_scalar (sib%rco2p(1,1,1),tend%rco2t(1),ng,'RCO2P')

   do nsc=1,itracer
      write(sname,'(a7,i3.3)') 'TRACERP',nsc
      if (allocated(tracer(nsc)%tracert))  &
         CALL vtables_scalar (tracer(nsc)%tracerp(1,1,1)  &
                             ,tracer(nsc)%tracert(1),ng,sname)
   enddo

return
END SUBROUTINE filltab_tend

END MODULE mem_tend
