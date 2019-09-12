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

Module mem_micro

implicit none

   Type micro_vars
   
      ! Variables to be dimensioned by (nzp,nxp,nyp)
   real, allocatable, dimension(:,:,:) :: &
                          rcp,rdp,rrp,rpp,rsp,rap,rgp,rhp           &
                         ,ccp,cdp,crp,cpp,csp,cap,cgp,chp           &
                         ,cifnp,q2,q6,q7                            &
                         ,pcpvd,pcpvr,pcpvp,pcpvs,pcpva,pcpvg,pcpvh &
                         !Aerosol categories mass and number
                         ,cccnp,gccnp,cccmp,gccmp                   &
                         ,md1np,md2np,md1mp,md2mp                   &
                         ,salt_film_np,salt_jet_np,salt_spum_np     &
                         ,salt_film_mp,salt_jet_mp,salt_spum_mp     &
                         ,regen_aero1_np,regen_aero1_mp             &
                         ,regen_aero2_np,regen_aero2_mp             &
                         !Immersion freezing nuclei tracking
                         ,immercp,immerdp,immerrp,ifnnucp           &
                         !Aerosol tracking variables
                         ,cnmcp,cnmdp,cnmrp,cnmpp,cnmsp             &
                         ,cnmap,cnmgp,cnmhp                         &
                         ,dnmcp,dnmdp,dnmrp,dnmpp,dnmsp             &
                         ,dnmap,dnmgp,dnmhp                         &
                         ,dincp,dindp,dinrp,dinpp,dinsp             &
                         ,dinap,dingp,dinhp                         &
                         ,snmcp,snmdp,snmrp,snmpp,snmsp             &
                         ,snmap,snmgp,snmhp                         &
                         ,resol_aero1_mp,resol_aero2_mp &
           ! MICRO BUDGET PROCESSES (imbudget >=1)
           ,latheatvap,latheatfrz,nuccldrt,cld2raint,ice2raint,nucicert    &
           ,vapliqt,vapicet,melticet,rimecldt,rain2icet,aggregatet         &
           ,latheatvapt,latheatfrzt,nuccldct,nucicect,sipct,siprt,collpsct &
           ! MICRO BUDGET PROCESSES (imbudget >=2)
           ,inuchomrt,inuccontrt,inucifnrt,inuchazrt,vapcldt,vapraint      &
           ,inuchomct,inuccontct,inucifnct,inuchazct                       &
           ,vapprist,vapsnowt,vapaggrt,vapgraut,vaphailt,vapdrizt          &
           ,meltprist,meltsnowt,meltaggrt,meltgraut,melthailt              &
           ,rimecldsnowt,rimecldaggrt,rimecldgraut,rimecldhailt            &
           ,rain2prt,rain2snt,rain2agt,rain2grt,rain2hat                   &
           ,aggrselfprist,aggrselfsnowt,aggrprissnowt                      &
           ! MICRO BUDGET PROCESSES (imbudget >=3)
           ,dust1cldrt,dust2cldrt,dust1drzrt,dust2drzrt

   ! Variables to be dimensioned by (nxp,nyp)
   real, allocatable, dimension(:,:) :: &
                          accpr,accpp,accps,accpa,accpg,accph,accpd &
                         ,pcprr,pcprp,pcprs,pcpra,pcprg,pcprh,pcprd &
                         ,pcpg,qpcpg,dpcpg                          &
                         !Accumulated aerosols and accumulation rate
                         ,accpdust,pcprdust,accpaero,pcpraero
                          
   End Type               
                          
   type (micro_vars), allocatable :: micro_g(:), microm_g(:)
                          
Contains                  
                          
!##############################################################################
Subroutine alloc_micro (micro,n1,n2,n3)

use micphys

implicit none          

   type (micro_vars) :: micro
   integer, intent(in) :: n1,n2,n3

! Allocate arrays based on options (if necessary)
      if (level >= 0) then
         if(iaerosol > 0) then
            allocate (micro%cccnp(n1,n2,n3))
            allocate (micro%cccmp(n1,n2,n3))
            allocate (micro%gccnp(n1,n2,n3))
            allocate (micro%gccmp(n1,n2,n3))
         endif
         if(idust > 0) then
            allocate (micro%md1np(n1,n2,n3))
            allocate (micro%md2np(n1,n2,n3))
            allocate (micro%md1mp(n1,n2,n3))
            allocate (micro%md2mp(n1,n2,n3))
         endif
         if(isalt > 0) then
            allocate (micro%salt_film_np(n1,n2,n3))
            allocate (micro%salt_jet_np(n1,n2,n3))
            allocate (micro%salt_spum_np(n1,n2,n3))
            allocate (micro%salt_film_mp(n1,n2,n3))
            allocate (micro%salt_jet_mp(n1,n2,n3))
            allocate (micro%salt_spum_mp(n1,n2,n3))
         endif
      endif

      if (level >= 2) then
         allocate (micro%rcp(n1,n2,n3))
      endif

      if (level >= 3) then
         if(ipris>=5 .and. (iifn==1.or.iifn==2)) allocate (micro%cifnp(n1,n2,n3))
         if(idriz >= 1)  then
            allocate (micro%rdp(n1,n2,n3))
            allocate (micro%accpd(n2,n3))
            allocate (micro%pcprd(n2,n3))
            allocate (micro%pcpvd(n1,n2,n3))
         endif
         if(irain >= 1)  then
            allocate (micro%rrp(n1,n2,n3))
            allocate (micro%accpr(n2,n3))
            allocate (micro%pcprr(n2,n3))
            allocate (micro%pcpvr(n1,n2,n3))
            allocate (micro%q2(n1,n2,n3))
         endif
         if(ipris >= 1)  then
            allocate (micro%rpp(n1,n2,n3))
            allocate (micro%accpp(n2,n3))
            allocate (micro%pcprp(n2,n3))
            allocate (micro%pcpvp(n1,n2,n3))
         endif
         if(isnow >= 1)  then
            allocate (micro%rsp(n1,n2,n3))
            allocate (micro%accps(n2,n3))
            allocate (micro%pcprs(n2,n3))
            allocate (micro%pcpvs(n1,n2,n3))
         endif
         if(iaggr >= 1)  then
            allocate (micro%rap(n1,n2,n3))
            allocate (micro%accpa(n2,n3))
            allocate (micro%pcpra(n2,n3))
            allocate (micro%pcpva(n1,n2,n3))
         endif
         if(igraup >= 1) then
            allocate (micro%rgp(n1,n2,n3))
            allocate (micro%accpg(n2,n3))
            allocate (micro%pcprg(n2,n3))
            allocate (micro%pcpvg(n1,n2,n3))
            allocate (micro%q6(n1,n2,n3))
         endif
         if(ihail >= 1)  then
            allocate (micro%rhp(n1,n2,n3))
            allocate (micro%accph(n2,n3))
            allocate (micro%pcprh(n2,n3))
            allocate (micro%pcpvh(n1,n2,n3))
            allocate (micro%q7(n1,n2,n3))
         endif

         if(jnmb(1) >= 5) allocate (micro%ccp(n1,n2,n3))
         if(jnmb(8) >= 5) allocate (micro%cdp(n1,n2,n3))
         if(jnmb(2) >= 5) allocate (micro%crp(n1,n2,n3))
         if(jnmb(3) >= 5) allocate (micro%cpp(n1,n2,n3))
         if(jnmb(4) >= 5) allocate (micro%csp(n1,n2,n3))
         if(jnmb(5) >= 5) allocate (micro%cap(n1,n2,n3))
         if(jnmb(6) >= 5) allocate (micro%cgp(n1,n2,n3))
         if(jnmb(7) >= 5) allocate (micro%chp(n1,n2,n3))

         allocate (micro%pcpg(n2,n3))
         allocate (micro%qpcpg(n2,n3))
         allocate (micro%dpcpg(n2,n3))

         if(iccnlev >= 2) then
           allocate (micro%regen_aero1_np(n1,n2,n3))
           allocate (micro%regen_aero1_mp(n1,n2,n3))
           allocate (micro%regen_aero2_np(n1,n2,n3))
           allocate (micro%regen_aero2_mp(n1,n2,n3))
           if(jnmb(1) >= 1) allocate (micro%cnmcp(n1,n2,n3))
           if(jnmb(8) >= 1) allocate (micro%cnmdp(n1,n2,n3))
           if(jnmb(2) >= 1) allocate (micro%cnmrp(n1,n2,n3))
           if(jnmb(3) >= 1) allocate (micro%cnmpp(n1,n2,n3))
           if(jnmb(4) >= 1) allocate (micro%cnmsp(n1,n2,n3))
           if(jnmb(5) >= 1) allocate (micro%cnmap(n1,n2,n3))
           if(jnmb(6) >= 1) allocate (micro%cnmgp(n1,n2,n3))
           if(jnmb(7) >= 1) allocate (micro%cnmhp(n1,n2,n3))
           allocate (micro%accpaero(n2,n3))
           allocate (micro%pcpraero(n2,n3))
           if(itrkdust==1) then
            allocate (micro%accpdust(n2,n3))
            allocate (micro%pcprdust(n2,n3))
            if(jnmb(1) >= 1) allocate (micro%dnmcp(n1,n2,n3))
            if(jnmb(8) >= 1) allocate (micro%dnmdp(n1,n2,n3))
            if(jnmb(2) >= 1) allocate (micro%dnmrp(n1,n2,n3))
            if(jnmb(3) >= 1) allocate (micro%dnmpp(n1,n2,n3))
            if(jnmb(4) >= 1) allocate (micro%dnmsp(n1,n2,n3))
            if(jnmb(5) >= 1) allocate (micro%dnmap(n1,n2,n3))
            if(jnmb(6) >= 1) allocate (micro%dnmgp(n1,n2,n3))
            if(jnmb(7) >= 1) allocate (micro%dnmhp(n1,n2,n3))
           endif
           if(itrkdustifn==1) then
            if(jnmb(1) >= 1) allocate (micro%dincp(n1,n2,n3))
            if(jnmb(8) >= 1) allocate (micro%dindp(n1,n2,n3))
            if(jnmb(2) >= 1) allocate (micro%dinrp(n1,n2,n3))
            if(jnmb(3) >= 1) allocate (micro%dinpp(n1,n2,n3))
            if(jnmb(4) >= 1) allocate (micro%dinsp(n1,n2,n3))
            if(jnmb(5) >= 1) allocate (micro%dinap(n1,n2,n3))
            if(jnmb(6) >= 1) allocate (micro%dingp(n1,n2,n3))
            if(jnmb(7) >= 1) allocate (micro%dinhp(n1,n2,n3))
           endif
           if(itrkepsilon==1) then
            allocate (micro%resol_aero1_mp(n1,n2,n3))
            allocate (micro%resol_aero2_mp(n1,n2,n3))
            if(jnmb(1) >= 1) allocate (micro%snmcp(n1,n2,n3))
            if(jnmb(8) >= 1) allocate (micro%snmdp(n1,n2,n3))
            if(jnmb(2) >= 1) allocate (micro%snmrp(n1,n2,n3))
            if(jnmb(3) >= 1) allocate (micro%snmpp(n1,n2,n3))
            if(jnmb(4) >= 1) allocate (micro%snmsp(n1,n2,n3))
            if(jnmb(5) >= 1) allocate (micro%snmap(n1,n2,n3))
            if(jnmb(6) >= 1) allocate (micro%snmgp(n1,n2,n3))
            if(jnmb(7) >= 1) allocate (micro%snmhp(n1,n2,n3))
           endif
         endif

         if(iifn==3 .and. iccnlev>=1) then
           if(jnmb(1) >= 5) allocate (micro%ifnnucp(n1,n2,n3))
           if(jnmb(1) >= 5) allocate (micro%immercp(n1,n2,n3))
           if(jnmb(8) >= 5) allocate (micro%immerdp(n1,n2,n3))
           if(jnmb(2) >= 5) allocate (micro%immerrp(n1,n2,n3))
         endif

         !COMPUTE AND OUTPUT MICRO BUDGET PROCESSES
         if(imbudget>=1) then
           allocate (micro%latheatvap(n1,n2,n3))
           allocate (micro%latheatfrz(n1,n2,n3))
           allocate (micro%nuccldrt(n1,n2,n3))
           allocate (micro%nuccldct(n1,n2,n3))
           allocate (micro%cld2raint(n1,n2,n3))
           allocate (micro%ice2raint(n1,n2,n3))
           allocate (micro%nucicert(n1,n2,n3))
           allocate (micro%nucicect(n1,n2,n3))
           allocate (micro%siprt(n1,n2,n3))
           allocate (micro%sipct(n1,n2,n3))
           allocate (micro%collpsct(n1,n2,n3))
           allocate (micro%vapliqt(n1,n2,n3))
           allocate (micro%vapicet(n1,n2,n3))
           allocate (micro%melticet(n1,n2,n3))
           allocate (micro%rimecldt(n1,n2,n3))
           allocate (micro%rain2icet(n1,n2,n3))
           allocate (micro%aggregatet(n1,n2,n3))
           allocate (micro%latheatvapt(n1,n2,n3))
           allocate (micro%latheatfrzt(n1,n2,n3))
         endif
         if(imbudget>=2) then
           allocate (micro%inuchomrt(n1,n2,n3))
           allocate (micro%inuccontrt(n1,n2,n3))
           allocate (micro%inucifnrt(n1,n2,n3))
           allocate (micro%inuchazrt(n1,n2,n3))
           allocate (micro%inuchomct(n1,n2,n3))
           allocate (micro%inuccontct(n1,n2,n3))
           allocate (micro%inucifnct(n1,n2,n3))
           allocate (micro%inuchazct(n1,n2,n3))
           allocate (micro%vapcldt(n1,n2,n3))
           allocate (micro%vapraint(n1,n2,n3))
           allocate (micro%vapprist(n1,n2,n3))
           allocate (micro%vapsnowt(n1,n2,n3))
           allocate (micro%vapaggrt(n1,n2,n3))
           allocate (micro%vapgraut(n1,n2,n3))
           allocate (micro%vaphailt(n1,n2,n3))
           allocate (micro%vapdrizt(n1,n2,n3))
           allocate (micro%meltprist(n1,n2,n3))
           allocate (micro%meltsnowt(n1,n2,n3))
           allocate (micro%meltaggrt(n1,n2,n3))
           allocate (micro%meltgraut(n1,n2,n3))
           allocate (micro%melthailt(n1,n2,n3))
           allocate (micro%rimecldsnowt(n1,n2,n3))
           allocate (micro%rimecldaggrt(n1,n2,n3))
           allocate (micro%rimecldgraut(n1,n2,n3))
           allocate (micro%rimecldhailt(n1,n2,n3))
           allocate (micro%rain2prt(n1,n2,n3))
           allocate (micro%rain2snt(n1,n2,n3))
           allocate (micro%rain2agt(n1,n2,n3))
           allocate (micro%rain2grt(n1,n2,n3))
           allocate (micro%rain2hat(n1,n2,n3))
           allocate (micro%aggrselfprist(n1,n2,n3))
           allocate (micro%aggrselfsnowt(n1,n2,n3))
           allocate (micro%aggrprissnowt(n1,n2,n3))
         endif
         if(imbudget==3 .and. idust>=1) then
           allocate (micro%dust1cldrt(n1,n2,n3))
           allocate (micro%dust2cldrt(n1,n2,n3))
           allocate (micro%dust1drzrt(n1,n2,n3))
           allocate (micro%dust2drzrt(n1,n2,n3))
         endif
      endif

return
END SUBROUTINE alloc_micro

!##############################################################################
Subroutine dealloc_micro (micro)

implicit none

   type (micro_vars) :: micro
   
   if (allocated(micro%rcp))     deallocate (micro%rcp)
   if (allocated(micro%rdp))     deallocate (micro%rdp)
   if (allocated(micro%rrp))     deallocate (micro%rrp)
   if (allocated(micro%rpp))     deallocate (micro%rpp)
   if (allocated(micro%rsp))     deallocate (micro%rsp)
   if (allocated(micro%rap))     deallocate (micro%rap)
   if (allocated(micro%rgp))     deallocate (micro%rgp)
   if (allocated(micro%rhp))     deallocate (micro%rhp)
   if (allocated(micro%ccp))     deallocate (micro%ccp)
   if (allocated(micro%cdp))     deallocate (micro%cdp)
   if (allocated(micro%crp))     deallocate (micro%crp)
   if (allocated(micro%cpp))     deallocate (micro%cpp)
   if (allocated(micro%csp))     deallocate (micro%csp)
   if (allocated(micro%cap))     deallocate (micro%cap)
   if (allocated(micro%cgp))     deallocate (micro%cgp)
   if (allocated(micro%chp))     deallocate (micro%chp)
   if (allocated(micro%cifnp))   deallocate (micro%cifnp)
   if (allocated(micro%q2))      deallocate (micro%q2)
   if (allocated(micro%q6))      deallocate (micro%q6)
   if (allocated(micro%q7))      deallocate (micro%q7)

   if (allocated(micro%cccnp))   deallocate (micro%cccnp)
   if (allocated(micro%gccnp))   deallocate (micro%gccnp)
   if (allocated(micro%cccmp))   deallocate (micro%cccmp)
   if (allocated(micro%gccmp))   deallocate (micro%gccmp)
   if (allocated(micro%md1np))   deallocate (micro%md1np)
   if (allocated(micro%md2np))   deallocate (micro%md2np)
   if (allocated(micro%md1mp))   deallocate (micro%md1mp)
   if (allocated(micro%md2mp))   deallocate (micro%md2mp)
   if (allocated(micro%salt_film_np))  deallocate (micro%salt_film_np)
   if (allocated(micro%salt_jet_np))   deallocate (micro%salt_jet_np)
   if (allocated(micro%salt_spum_np))  deallocate (micro%salt_spum_np)
   if (allocated(micro%salt_film_mp))  deallocate (micro%salt_film_mp)
   if (allocated(micro%salt_jet_mp))   deallocate (micro%salt_jet_mp)
   if (allocated(micro%salt_spum_mp))  deallocate (micro%salt_spum_mp)
   if (allocated(micro%regen_aero1_np)) deallocate (micro%regen_aero1_np)
   if (allocated(micro%regen_aero1_mp)) deallocate (micro%regen_aero1_mp)
   if (allocated(micro%regen_aero2_np)) deallocate (micro%regen_aero2_np)
   if (allocated(micro%regen_aero2_mp)) deallocate (micro%regen_aero2_mp)

   if (allocated(micro%immercp)) deallocate (micro%immercp)
   if (allocated(micro%immerdp)) deallocate (micro%immerdp)
   if (allocated(micro%immerrp)) deallocate (micro%immerrp)
   if (allocated(micro%ifnnucp)) deallocate (micro%ifnnucp)

   if (allocated(micro%cnmcp))   deallocate (micro%cnmcp)
   if (allocated(micro%cnmdp))   deallocate (micro%cnmdp)
   if (allocated(micro%cnmrp))   deallocate (micro%cnmrp)
   if (allocated(micro%cnmpp))   deallocate (micro%cnmpp)
   if (allocated(micro%cnmsp))   deallocate (micro%cnmsp)
   if (allocated(micro%cnmap))   deallocate (micro%cnmap)
   if (allocated(micro%cnmgp))   deallocate (micro%cnmgp)
   if (allocated(micro%cnmhp))   deallocate (micro%cnmhp)
   if (allocated(micro%accpaero))deallocate (micro%accpaero)
   if (allocated(micro%pcpraero))deallocate (micro%pcpraero)

   if (allocated(micro%dnmcp))   deallocate (micro%dnmcp)
   if (allocated(micro%dnmdp))   deallocate (micro%dnmdp)
   if (allocated(micro%dnmrp))   deallocate (micro%dnmrp)
   if (allocated(micro%dnmpp))   deallocate (micro%dnmpp)
   if (allocated(micro%dnmsp))   deallocate (micro%dnmsp)
   if (allocated(micro%dnmap))   deallocate (micro%dnmap)
   if (allocated(micro%dnmgp))   deallocate (micro%dnmgp)
   if (allocated(micro%dnmhp))   deallocate (micro%dnmhp)
   if (allocated(micro%accpdust))deallocate (micro%accpdust)
   if (allocated(micro%pcprdust))deallocate (micro%pcprdust)

   if (allocated(micro%dincp))   deallocate (micro%dincp)
   if (allocated(micro%dindp))   deallocate (micro%dindp)
   if (allocated(micro%dinrp))   deallocate (micro%dinrp)
   if (allocated(micro%dinpp))   deallocate (micro%dinpp)
   if (allocated(micro%dinsp))   deallocate (micro%dinsp)
   if (allocated(micro%dinap))   deallocate (micro%dinap)
   if (allocated(micro%dingp))   deallocate (micro%dingp)
   if (allocated(micro%dinhp))   deallocate (micro%dinhp)

   if (allocated(micro%snmcp))   deallocate (micro%snmcp)
   if (allocated(micro%snmdp))   deallocate (micro%snmdp)
   if (allocated(micro%snmrp))   deallocate (micro%snmrp)
   if (allocated(micro%snmpp))   deallocate (micro%snmpp)
   if (allocated(micro%snmsp))   deallocate (micro%snmsp)
   if (allocated(micro%snmap))   deallocate (micro%snmap)
   if (allocated(micro%snmgp))   deallocate (micro%snmgp)
   if (allocated(micro%snmhp))   deallocate (micro%snmhp)
   if (allocated(micro%resol_aero1_mp)) deallocate (micro%resol_aero1_mp)
   if (allocated(micro%resol_aero2_mp)) deallocate (micro%resol_aero2_mp)

   if (allocated(micro%pcpvr))   deallocate (micro%pcpvr)
   if (allocated(micro%pcpvp))   deallocate (micro%pcpvp)
   if (allocated(micro%pcpvs))   deallocate (micro%pcpvs)
   if (allocated(micro%pcpva))   deallocate (micro%pcpva)
   if (allocated(micro%pcpvg))   deallocate (micro%pcpvg)
   if (allocated(micro%pcpvh))   deallocate (micro%pcpvh)
   if (allocated(micro%pcpvd))   deallocate (micro%pcpvd)

   if (allocated(micro%accpr))   deallocate (micro%accpr)
   if (allocated(micro%accpp))   deallocate (micro%accpp)
   if (allocated(micro%accps))   deallocate (micro%accps)
   if (allocated(micro%accpa))   deallocate (micro%accpa)
   if (allocated(micro%accpg))   deallocate (micro%accpg)
   if (allocated(micro%accph))   deallocate (micro%accph)
   if (allocated(micro%accpd))   deallocate (micro%accpd)
   if (allocated(micro%pcprr))   deallocate (micro%pcprr)
   if (allocated(micro%pcprp))   deallocate (micro%pcprp)
   if (allocated(micro%pcprs))   deallocate (micro%pcprs)
   if (allocated(micro%pcpra))   deallocate (micro%pcpra)
   if (allocated(micro%pcprg))   deallocate (micro%pcprg)
   if (allocated(micro%pcprh))   deallocate (micro%pcprh)
   if (allocated(micro%pcprd))   deallocate (micro%pcprd)
   if (allocated(micro%pcpg))    deallocate (micro%pcpg)
   if (allocated(micro%qpcpg))   deallocate (micro%qpcpg)
   if (allocated(micro%dpcpg))   deallocate (micro%dpcpg)

    !COMPUTE AND OUTPUT MICRO BUDGET PROCESSES
    if (allocated(micro%latheatvap))    deallocate (micro%latheatvap)
    if (allocated(micro%latheatfrz))    deallocate (micro%latheatfrz)
    if (allocated(micro%nuccldrt))      deallocate (micro%nuccldrt)
    if (allocated(micro%nuccldct))      deallocate (micro%nuccldct)
    if (allocated(micro%cld2raint))     deallocate (micro%cld2raint)
    if (allocated(micro%ice2raint))     deallocate (micro%ice2raint)
    if (allocated(micro%nucicert))      deallocate (micro%nucicert)
    if (allocated(micro%nucicect))      deallocate (micro%nucicect)
    if (allocated(micro%siprt))         deallocate (micro%siprt)
    if (allocated(micro%sipct))         deallocate (micro%sipct)
    if (allocated(micro%collpsct))      deallocate (micro%collpsct)
    if (allocated(micro%vapliqt))       deallocate (micro%vapliqt)
    if (allocated(micro%vapicet))       deallocate (micro%vapicet)
    if (allocated(micro%melticet))      deallocate (micro%melticet)
    if (allocated(micro%rimecldt))      deallocate (micro%rimecldt)
    if (allocated(micro%rain2icet))     deallocate (micro%rain2icet)
    if (allocated(micro%aggregatet))    deallocate (micro%aggregatet)
    if (allocated(micro%latheatvapt))   deallocate (micro%latheatvapt)
    if (allocated(micro%latheatfrzt))   deallocate (micro%latheatfrzt)

    if (allocated(micro%inuchomrt))     deallocate (micro%inuchomrt)
    if (allocated(micro%inuccontrt))    deallocate (micro%inuccontrt)
    if (allocated(micro%inucifnrt))     deallocate (micro%inucifnrt)
    if (allocated(micro%inuchazrt))     deallocate (micro%inuchazrt)
    if (allocated(micro%inuchomct))     deallocate (micro%inuchomct)
    if (allocated(micro%inuccontct))    deallocate (micro%inuccontct)
    if (allocated(micro%inucifnct))     deallocate (micro%inucifnct)
    if (allocated(micro%inuchazct))     deallocate (micro%inuchazct)
    if (allocated(micro%vapcldt))       deallocate (micro%vapcldt)
    if (allocated(micro%vapraint))      deallocate (micro%vapraint)
    if (allocated(micro%vapprist))      deallocate (micro%vapprist)
    if (allocated(micro%vapsnowt))      deallocate (micro%vapsnowt)
    if (allocated(micro%vapaggrt))      deallocate (micro%vapaggrt)
    if (allocated(micro%vapgraut))      deallocate (micro%vapgraut)
    if (allocated(micro%vaphailt))      deallocate (micro%vaphailt)
    if (allocated(micro%vapdrizt))      deallocate (micro%vapdrizt)
    if (allocated(micro%meltprist))     deallocate (micro%meltprist)
    if (allocated(micro%meltsnowt))     deallocate (micro%meltsnowt)
    if (allocated(micro%meltaggrt))     deallocate (micro%meltaggrt)
    if (allocated(micro%meltgraut))     deallocate (micro%meltgraut)
    if (allocated(micro%melthailt))     deallocate (micro%melthailt)
    if (allocated(micro%rimecldsnowt))  deallocate (micro%rimecldsnowt)
    if (allocated(micro%rimecldaggrt))  deallocate (micro%rimecldaggrt)
    if (allocated(micro%rimecldgraut))  deallocate (micro%rimecldgraut)
    if (allocated(micro%rimecldhailt))  deallocate (micro%rimecldhailt)
    if (allocated(micro%rain2prt))      deallocate (micro%rain2prt)
    if (allocated(micro%rain2snt))      deallocate (micro%rain2snt)
    if (allocated(micro%rain2agt))      deallocate (micro%rain2agt)
    if (allocated(micro%rain2grt))      deallocate (micro%rain2grt)
    if (allocated(micro%rain2hat))      deallocate (micro%rain2hat)
    if (allocated(micro%aggrselfprist)) deallocate (micro%aggrselfprist)
    if (allocated(micro%aggrselfsnowt)) deallocate (micro%aggrselfsnowt)
    if (allocated(micro%aggrprissnowt)) deallocate (micro%aggrprissnowt)

    if (allocated(micro%dust1cldrt))        deallocate (micro%dust1cldrt)
    if (allocated(micro%dust2cldrt))        deallocate (micro%dust2cldrt)
    if (allocated(micro%dust1drzrt))        deallocate (micro%dust1drzrt)
    if (allocated(micro%dust2drzrt))        deallocate (micro%dust2drzrt)

return
END SUBROUTINE dealloc_micro

!##############################################################################
Subroutine filltab_micro (micro,microm,imean,n1,n2,n3,ng)

use var_tables

implicit none

   type (micro_vars) :: micro,microm
   integer, intent(in) :: imean,n1,n2,n3,ng
   integer :: npts

! Fill arrays into variable tables

   npts=n1*n2*n3
   if (allocated(micro%rcp))   &
      CALL vtables2 (micro%rcp(1,1,1),microm%rcp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RCP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%rdp))   &
      CALL vtables2 (micro%rdp(1,1,1),microm%rdp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RDP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%rrp))   &
      CALL vtables2 (micro%rrp(1,1,1),microm%rrp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RRP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%rpp))   &
      CALL vtables2 (micro%rpp(1,1,1),microm%rpp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RPP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%rsp))   &
      CALL vtables2 (micro%rsp(1,1,1),microm%rsp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RSP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%rap))   &
      CALL vtables2 (micro%rap(1,1,1),microm%rap(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RAP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%rgp))   &
      CALL vtables2 (micro%rgp(1,1,1),microm%rgp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RGP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%rhp))   &
      CALL vtables2 (micro%rhp(1,1,1),microm%rhp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RHP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%ccp))   &
      CALL vtables2 (micro%ccp(1,1,1),microm%ccp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'CCP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%cdp))   &
      CALL vtables2 (micro%cdp(1,1,1),microm%cdp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'CDP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%crp))   &
      CALL vtables2 (micro%crp(1,1,1),microm%crp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'CRP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%cpp))   &
      CALL vtables2 (micro%cpp(1,1,1),microm%cpp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'CPP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%csp))   &
      CALL vtables2 (micro%csp(1,1,1),microm%csp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'CSP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%cap))   &
      CALL vtables2 (micro%cap(1,1,1),microm%cap(1,1,1)  &
                 ,ng, npts, imean,  &
                 'CAP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%cgp))   &
      CALL vtables2 (micro%cgp(1,1,1),microm%cgp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'CGP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%chp))   &
      CALL vtables2 (micro%chp(1,1,1),microm%chp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'CHP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%cifnp)) &
      CALL vtables2 (micro%cifnp(1,1,1),microm%cifnp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'CIFNP :3:anal:mpti:mpt3:mpt1')

   if (allocated(micro%q2))   &
      CALL vtables2 (micro%q2(1,1,1),microm%q2(1,1,1)  &
                 ,ng, npts, imean,  &
                 'Q2 :3:anal:mpti:mpt3')
   if (allocated(micro%q6)) &
      CALL vtables2 (micro%q6(1,1,1),microm%q6(1,1,1)  &
                 ,ng, npts, imean,  &
                 'Q6 :3:anal:mpti:mpt3')
   if (allocated(micro%q7)) &
      CALL vtables2 (micro%q7(1,1,1),microm%q7(1,1,1)  &
                 ,ng, npts, imean,  &
                 'Q7 :3:anal:mpti:mpt3')

!Aerosol categories mass and number
   if (allocated(micro%cccnp)) &
      CALL vtables2 (micro%cccnp(1,1,1),microm%cccnp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'CCCNP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%gccnp)) &
      CALL vtables2 (micro%gccnp(1,1,1),microm%gccnp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'GCCNP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%cccmp)) &
      CALL vtables2 (micro%cccmp(1,1,1),microm%cccmp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'CCCMP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%gccmp)) &
      CALL vtables2 (micro%gccmp(1,1,1),microm%gccmp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'GCCMP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%md1np)) &
      CALL vtables2 (micro%md1np(1,1,1),microm%md1np(1,1,1)  &
                 ,ng, npts, imean,  &
                 'MD1NP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%md2np)) &
      CALL vtables2 (micro%md2np(1,1,1),microm%md2np(1,1,1)  &
                 ,ng, npts, imean,  &
                 'MD2NP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%md1mp)) &
      CALL vtables2 (micro%md1mp(1,1,1),microm%md1mp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'MD1MP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%md2mp)) &
      CALL vtables2 (micro%md2mp(1,1,1),microm%md2mp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'MD2MP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%salt_film_np)) &
      CALL vtables2 (micro%salt_film_np(1,1,1),microm%salt_film_np(1,1,1)  &
                 ,ng, npts, imean,  &
                 'SALT_FILM_NP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%salt_jet_np)) &
      CALL vtables2 (micro%salt_jet_np(1,1,1),microm%salt_jet_np(1,1,1)  &
                 ,ng, npts, imean,  &
                 'SALT_JET_NP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%salt_spum_np)) &
      CALL vtables2 (micro%salt_spum_np(1,1,1),microm%salt_spum_np(1,1,1)  &
                 ,ng, npts, imean,  &
                 'SALT_SPUM_NP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%salt_film_mp)) &
      CALL vtables2 (micro%salt_film_mp(1,1,1),microm%salt_film_mp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'SALT_FILM_MP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%salt_jet_mp)) &
      CALL vtables2 (micro%salt_jet_mp(1,1,1),microm%salt_jet_mp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'SALT_JET_MP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%salt_spum_mp)) &
      CALL vtables2 (micro%salt_spum_mp(1,1,1),microm%salt_spum_mp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'SALT_SPUM_MP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%regen_aero1_np)) &
      CALL vtables2 (micro%regen_aero1_np(1,1,1),microm%regen_aero1_np(1,1,1)  &
                 ,ng, npts, imean,  &
                 'REGEN_AERO1_NP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%regen_aero1_mp)) &
      CALL vtables2 (micro%regen_aero1_mp(1,1,1),microm%regen_aero1_mp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'REGEN_AERO1_MP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%regen_aero2_np)) &
      CALL vtables2 (micro%regen_aero2_np(1,1,1),microm%regen_aero2_np(1,1,1)  &
                 ,ng, npts, imean,  &
                 'REGEN_AERO2_NP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%regen_aero2_mp)) &
      CALL vtables2 (micro%regen_aero2_mp(1,1,1),microm%regen_aero2_mp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'REGEN_AERO2_MP :3:anal:mpti:mpt3:mpt1')

!Immersion freezing nuclei number tracking variables
   if (allocated(micro%immercp)) &
      CALL vtables2 (micro%immercp(1,1,1),microm%immercp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'IMMERCP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%immerdp)) &
      CALL vtables2 (micro%immerdp(1,1,1),microm%immerdp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'IMMERDP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%immerrp)) &
      CALL vtables2 (micro%immerrp(1,1,1),microm%immerrp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'IMMERRP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%ifnnucp)) &
      CALL vtables2 (micro%ifnnucp(1,1,1),microm%ifnnucp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'IFNNUCP :3:anal:mpti:mpt3:mpt1')

!Total aerosol mass-in-hydrometeors tracking variables
   if (allocated(micro%cnmcp)) &
      CALL vtables2 (micro%cnmcp(1,1,1),microm%cnmcp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'CNMCP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%cnmdp)) &
      CALL vtables2 (micro%cnmdp(1,1,1),microm%cnmdp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'CNMDP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%cnmrp)) &
      CALL vtables2 (micro%cnmrp(1,1,1),microm%cnmrp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'CNMRP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%cnmpp)) &
      CALL vtables2 (micro%cnmpp(1,1,1),microm%cnmpp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'CNMPP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%cnmsp)) &
      CALL vtables2 (micro%cnmsp(1,1,1),microm%cnmsp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'CNMSP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%cnmap)) &
      CALL vtables2 (micro%cnmap(1,1,1),microm%cnmap(1,1,1)  &
                 ,ng, npts, imean,  &
                 'CNMAP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%cnmgp)) &
      CALL vtables2 (micro%cnmgp(1,1,1),microm%cnmgp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'CNMGP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%cnmhp)) &
      CALL vtables2 (micro%cnmhp(1,1,1),microm%cnmhp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'CNMHP :3:anal:mpti:mpt3:mpt1')

!Dust aerosol mass-in-hydrometeors tracking variables
   if (allocated(micro%dnmcp)) &
      CALL vtables2 (micro%dnmcp(1,1,1),microm%dnmcp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'DNMCP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%dnmdp)) &
      CALL vtables2 (micro%dnmdp(1,1,1),microm%dnmdp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'DNMDP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%dnmrp)) &
      CALL vtables2 (micro%dnmrp(1,1,1),microm%dnmrp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'DNMRP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%dnmpp)) &
      CALL vtables2 (micro%dnmpp(1,1,1),microm%dnmpp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'DNMPP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%dnmsp)) &
      CALL vtables2 (micro%dnmsp(1,1,1),microm%dnmsp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'DNMSP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%dnmap)) &
      CALL vtables2 (micro%dnmap(1,1,1),microm%dnmap(1,1,1)  &
                 ,ng, npts, imean,  &
                 'DNMAP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%dnmgp)) &
      CALL vtables2 (micro%dnmgp(1,1,1),microm%dnmgp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'DNMGP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%dnmhp)) &
      CALL vtables2 (micro%dnmhp(1,1,1),microm%dnmhp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'DNMHP :3:anal:mpti:mpt3:mpt1')

!Dust as Ice nuclei mass-in-hydrometeors tracking variables
   if (allocated(micro%dincp)) &
      CALL vtables2 (micro%dincp(1,1,1),microm%dincp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'DINCP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%dindp)) &
      CALL vtables2 (micro%dindp(1,1,1),microm%dindp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'DINDP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%dinrp)) &
      CALL vtables2 (micro%dinrp(1,1,1),microm%dinrp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'DINRP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%dinpp)) &
      CALL vtables2 (micro%dinpp(1,1,1),microm%dinpp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'DINPP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%dinsp)) &
      CALL vtables2 (micro%dinsp(1,1,1),microm%dinsp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'DINSP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%dinap)) &
      CALL vtables2 (micro%dinap(1,1,1),microm%dinap(1,1,1)  &
                 ,ng, npts, imean,  &
                 'DINAP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%dingp)) &
      CALL vtables2 (micro%dingp(1,1,1),microm%dingp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'DINGP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%dinhp)) &
      CALL vtables2 (micro%dinhp(1,1,1),microm%dinhp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'DINHP :3:anal:mpti:mpt3:mpt1')

!Total aerosol soluble mass-in-hydrometeors for aerosol tracking
   if (allocated(micro%snmcp)) &
      CALL vtables2 (micro%snmcp(1,1,1),microm%snmcp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'SNMCP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%snmdp)) &
      CALL vtables2 (micro%snmdp(1,1,1),microm%snmdp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'SNMDP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%snmrp)) &
      CALL vtables2 (micro%snmrp(1,1,1),microm%snmrp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'SNMRP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%snmpp)) &
      CALL vtables2 (micro%snmpp(1,1,1),microm%snmpp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'SNMPP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%snmsp)) &
      CALL vtables2 (micro%snmsp(1,1,1),microm%snmsp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'SNMSP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%snmap)) &
      CALL vtables2 (micro%snmap(1,1,1),microm%snmap(1,1,1)  &
                 ,ng, npts, imean,  &
                 'SNMAP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%snmgp)) &
      CALL vtables2 (micro%snmgp(1,1,1),microm%snmgp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'SNMGP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%snmhp)) &
      CALL vtables2 (micro%snmhp(1,1,1),microm%snmhp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'SNMHP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%resol_aero1_mp)) &
      CALL vtables2 (micro%resol_aero1_mp(1,1,1),microm%resol_aero1_mp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RESOL_AERO1_MP :3:anal:mpti:mpt3:mpt1')
   if (allocated(micro%resol_aero2_mp)) &
      CALL vtables2 (micro%resol_aero2_mp(1,1,1),microm%resol_aero2_mp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RESOL_AERO2_MP :3:anal:mpti:mpt3:mpt1')

!Vertical precipitation rates
   if (allocated(micro%pcpvr)) &
      CALL vtables2 (micro%pcpvr(1,1,1),microm%pcpvr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'PCPVR :3:anal:mpti:mpt3')
   if (allocated(micro%pcpvp)) &
      CALL vtables2 (micro%pcpvp(1,1,1),microm%pcpvp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'PCPVP :3:anal:mpti:mpt3')
   if (allocated(micro%pcpvs)) &
      CALL vtables2 (micro%pcpvs(1,1,1),microm%pcpvs(1,1,1)  &
                 ,ng, npts, imean,  &
                 'PCPVS :3:anal:mpti:mpt3')
   if (allocated(micro%pcpva)) &
      CALL vtables2 (micro%pcpva(1,1,1),microm%pcpva(1,1,1)  &
                 ,ng, npts, imean,  &
                 'PCPVA :3:anal:mpti:mpt3')
   if (allocated(micro%pcpvg)) &
      CALL vtables2 (micro%pcpvg(1,1,1),microm%pcpvg(1,1,1)  &
                 ,ng, npts, imean,  &
                 'PCPVG :3:anal:mpti:mpt3')
   if (allocated(micro%pcpvh)) &
      CALL vtables2 (micro%pcpvh(1,1,1),microm%pcpvh(1,1,1)  &
                 ,ng, npts, imean,  &
                 'PCPVH :3:anal:mpti:mpt3')
   if (allocated(micro%pcpvd)) &
      CALL vtables2 (micro%pcpvd(1,1,1),microm%pcpvd(1,1,1)  &
                 ,ng, npts, imean,  &
                 'PCPVD :3:anal:mpti:mpt3')

   !COMPUTE AND OUTPUT MICRO BUDGET PROCESSES
   if (allocated(micro%latheatvap)) &
      CALL vtables2 (micro%latheatvap(1,1,1),microm%latheatvap(1,1,1)  &
                 ,ng, npts, imean,  &
                 'LATHEATVAP :3:anal:mpt3')
   if (allocated(micro%latheatfrz)) &
      CALL vtables2 (micro%latheatfrz(1,1,1),microm%latheatfrz(1,1,1)  &
                 ,ng, npts, imean,  &
                 'LATHEATFRZ :3:anal:mpt3')
   if (allocated(micro%nuccldrt)) &
      CALL vtables2 (micro%nuccldrt(1,1,1),microm%nuccldrt(1,1,1)  &
                 ,ng, npts, imean,  &
                 'NUCCLDRT :3:anal:mpti:mpt3')
   if (allocated(micro%nuccldct)) &
      CALL vtables2 (micro%nuccldct(1,1,1),microm%nuccldct(1,1,1)  &
                 ,ng, npts, imean,  &
                 'NUCCLDCT :3:anal:mpti:mpt3')
   if (allocated(micro%cld2raint)) &
      CALL vtables2 (micro%cld2raint(1,1,1),microm%cld2raint(1,1,1)  &
                 ,ng, npts, imean,  &
                 'CLD2RAINT :3:anal:mpti:mpt3')
   if (allocated(micro%ice2raint)) &
      CALL vtables2 (micro%ice2raint(1,1,1),microm%ice2raint(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ICE2RAINT :3:anal:mpti:mpt3')
   if (allocated(micro%nucicert)) &
      CALL vtables2 (micro%nucicert(1,1,1),microm%nucicert(1,1,1)  &
                 ,ng, npts, imean,  &
                 'NUCICERT :3:anal:mpti:mpt3')
   if (allocated(micro%nucicect)) &
      CALL vtables2 (micro%nucicect(1,1,1),microm%nucicect(1,1,1)  &
                 ,ng, npts, imean,  &
                 'NUCICECT :3:anal:mpti:mpt3')
   if (allocated(micro%siprt)) &
      CALL vtables2 (micro%siprt(1,1,1),microm%siprt(1,1,1)  &
                 ,ng, npts, imean,  &
                 'SIPRT :3:anal:mpti:mpt3')
   if (allocated(micro%sipct)) &
      CALL vtables2 (micro%sipct(1,1,1),microm%sipct(1,1,1)  &
                 ,ng, npts, imean,  &
                 'SIPCT :3:anal:mpti:mpt3')
   if (allocated(micro%collpsct)) &
      CALL vtables2 (micro%collpsct(1,1,1),microm%collpsct(1,1,1)  &
                 ,ng, npts, imean,  &
                 'COLLPSCT :3:anal:mpti:mpt3')
   if (allocated(micro%vapliqt)) &
      CALL vtables2 (micro%vapliqt(1,1,1),microm%vapliqt(1,1,1)  &
                 ,ng, npts, imean,  &
                 'VAPLIQT :3:anal:mpti:mpt3')
   if (allocated(micro%vapicet)) &
      CALL vtables2 (micro%vapicet(1,1,1),microm%vapicet(1,1,1)  &
                 ,ng, npts, imean,  &
                 'VAPICET :3:anal:mpti:mpt3')
   if (allocated(micro%melticet)) &
      CALL vtables2 (micro%melticet(1,1,1),microm%melticet(1,1,1)  &
                 ,ng, npts, imean,  &
                 'MELTICET :3:anal:mpti:mpt3')
   if (allocated(micro%rimecldt)) &
      CALL vtables2 (micro%rimecldt(1,1,1),microm%rimecldt(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RIMECLDT :3:anal:mpti:mpt3')
   if (allocated(micro%rain2icet)) &
      CALL vtables2 (micro%rain2icet(1,1,1),microm%rain2icet(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RAIN2ICET :3:anal:mpti:mpt3')
   if (allocated(micro%aggregatet)) &
      CALL vtables2 (micro%aggregatet(1,1,1),microm%aggregatet(1,1,1)  &
                 ,ng, npts, imean,  &
                 'AGGREGATET :3:anal:mpti:mpt3')
   if (allocated(micro%latheatvapt)) &
      CALL vtables2 (micro%latheatvapt(1,1,1),microm%latheatvapt(1,1,1)  &
                 ,ng, npts, imean,  &
                 'LATHEATVAPT :3:anal:mpti:mpt3')
   if (allocated(micro%latheatfrzt)) &
      CALL vtables2 (micro%latheatfrzt(1,1,1),microm%latheatfrzt(1,1,1)  &
                 ,ng, npts, imean,  &
                 'LATHEATFRZT :3:anal:mpti:mpt3')

   if (allocated(micro%inuchomrt)) &
      CALL vtables2 (micro%inuchomrt(1,1,1),microm%inuchomrt(1,1,1)  &
                 ,ng, npts, imean,  &
                 'INUCHOMRT :3:anal:mpti:mpt3')
   if (allocated(micro%inuccontrt)) &
      CALL vtables2 (micro%inuccontrt(1,1,1),microm%inuccontrt(1,1,1)  &
                 ,ng, npts, imean,  &
                 'INUCCONTRT :3:anal:mpti:mpt3')
   if (allocated(micro%inucifnrt)) &
      CALL vtables2 (micro%inucifnrt(1,1,1),microm%inucifnrt(1,1,1)  &
                 ,ng, npts, imean,  &
                 'INUCIFNRT :3:anal:mpti:mpt3')
   if (allocated(micro%inuchazrt)) &
      CALL vtables2 (micro%inuchazrt(1,1,1),microm%inuchazrt(1,1,1)  &
                 ,ng, npts, imean,  &
                 'INUCHAZRT :3:anal:mpti:mpt3')
   if (allocated(micro%inuchomct)) &
      CALL vtables2 (micro%inuchomct(1,1,1),microm%inuchomct(1,1,1)  &
                 ,ng, npts, imean,  &
                 'INUCHOMCT :3:anal:mpti:mpt3')
   if (allocated(micro%inuccontct)) &
      CALL vtables2 (micro%inuccontct(1,1,1),microm%inuccontct(1,1,1)  &
                 ,ng, npts, imean,  &
                 'INUCCONTCT :3:anal:mpti:mpt3')
   if (allocated(micro%inucifnct)) &
      CALL vtables2 (micro%inucifnct(1,1,1),microm%inucifnct(1,1,1)  &
                 ,ng, npts, imean,  &
                 'INUCIFNCT :3:anal:mpti:mpt3')
   if (allocated(micro%inuchazct)) &
      CALL vtables2 (micro%inuchazct(1,1,1),microm%inuchazct(1,1,1)  &
                 ,ng, npts, imean,  &
                 'INUCHAZCT :3:anal:mpti:mpt3')
   if (allocated(micro%vapcldt)) &
      CALL vtables2 (micro%vapcldt(1,1,1),microm%vapcldt(1,1,1)  &
                 ,ng, npts, imean,  &
                 'VAPCLDT :3:anal:mpti:mpt3')
   if (allocated(micro%vapraint)) &
      CALL vtables2 (micro%vapraint(1,1,1),microm%vapraint(1,1,1)  &
                 ,ng, npts, imean,  &
                 'VAPRAINT :3:anal:mpti:mpt3')
   if (allocated(micro%vapprist)) &
      CALL vtables2 (micro%vapprist(1,1,1),microm%vapprist(1,1,1)  &
                 ,ng, npts, imean,  &
                 'VAPPRIST :3:anal:mpti:mpt3')
   if (allocated(micro%vapsnowt)) &
      CALL vtables2 (micro%vapsnowt(1,1,1),microm%vapsnowt(1,1,1)  &
                 ,ng, npts, imean,  &
                 'VAPSNOWT :3:anal:mpti:mpt3')
   if (allocated(micro%vapaggrt)) &
      CALL vtables2 (micro%vapaggrt(1,1,1),microm%vapaggrt(1,1,1)  &
                 ,ng, npts, imean,  &
                 'VAPAGGRT :3:anal:mpti:mpt3')
   if (allocated(micro%vapgraut)) &
      CALL vtables2 (micro%vapgraut(1,1,1),microm%vapgraut(1,1,1)  &
                 ,ng, npts, imean,  &
                 'VAPGRAUT :3:anal:mpti:mpt3')
   if (allocated(micro%vaphailt)) &
      CALL vtables2 (micro%vaphailt(1,1,1),microm%vaphailt(1,1,1)  &
                 ,ng, npts, imean,  &
                 'VAPHAILT :3:anal:mpti:mpt3')
   if (allocated(micro%vapdrizt)) &
      CALL vtables2 (micro%vapdrizt(1,1,1),microm%vapdrizt(1,1,1)  &
                 ,ng, npts, imean,  &
                 'VAPDRIZT :3:anal:mpti:mpt3') 
   if (allocated(micro%meltprist)) &
      CALL vtables2 (micro%meltprist(1,1,1),microm%meltprist(1,1,1)  &
                 ,ng, npts, imean,  &
                 'MELTPRIST :3:anal:mpti:mpt3')
   if (allocated(micro%meltsnowt)) &
      CALL vtables2 (micro%meltsnowt(1,1,1),microm%meltsnowt(1,1,1)  &
                 ,ng, npts, imean,  &
                 'MELTSNOWT :3:anal:mpti:mpt3')
   if (allocated(micro%meltaggrt)) &
      CALL vtables2 (micro%meltaggrt(1,1,1),microm%meltaggrt(1,1,1)  &
                 ,ng, npts, imean,  &
                 'MELTAGGRT :3:anal:mpti:mpt3')
   if (allocated(micro%meltgraut)) &
      CALL vtables2 (micro%meltgraut(1,1,1),microm%meltgraut(1,1,1)  &
                 ,ng, npts, imean,  &
                 'MELTGRAUT :3:anal:mpti:mpt3')
   if (allocated(micro%melthailt)) &
      CALL vtables2 (micro%melthailt(1,1,1),microm%melthailt(1,1,1)  &
                 ,ng, npts, imean,  &
                 'MELTHAILT :3:anal:mpti:mpt3')
   if (allocated(micro%rimecldsnowt)) &
      CALL vtables2 (micro%rimecldsnowt(1,1,1),microm%rimecldsnowt(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RIMECLDSNOWT :3:anal:mpti:mpt3')
   if (allocated(micro%rimecldaggrt)) &
      CALL vtables2 (micro%rimecldaggrt(1,1,1),microm%rimecldaggrt(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RIMECLDAGGRT :3:anal:mpti:mpt3')
   if (allocated(micro%rimecldgraut)) &
      CALL vtables2 (micro%rimecldgraut(1,1,1),microm%rimecldgraut(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RIMECLDGRAUT :3:anal:mpti:mpt3')
   if (allocated(micro%rimecldhailt)) &
      CALL vtables2 (micro%rimecldhailt(1,1,1),microm%rimecldhailt(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RIMECLDHAILT :3:anal:mpti:mpt3')
   if (allocated(micro%rain2prt)) &
      CALL vtables2 (micro%rain2prt(1,1,1),microm%rain2prt(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RAIN2PRT :3:anal:mpti:mpt3')
   if (allocated(micro%rain2snt)) &
      CALL vtables2 (micro%rain2snt(1,1,1),microm%rain2snt(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RAIN2SNT :3:anal:mpti:mpt3')
   if (allocated(micro%rain2agt)) &
      CALL vtables2 (micro%rain2agt(1,1,1),microm%rain2agt(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RAIN2AGT :3:anal:mpti:mpt3')
   if (allocated(micro%rain2grt)) &
      CALL vtables2 (micro%rain2grt(1,1,1),microm%rain2grt(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RAIN2GRT :3:anal:mpti:mpt3')
   if (allocated(micro%rain2hat)) &
      CALL vtables2 (micro%rain2hat(1,1,1),microm%rain2hat(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RAIN2HAT :3:anal:mpti:mpt3')
   if (allocated(micro%aggrselfprist)) &
      CALL vtables2 (micro%aggrselfprist(1,1,1),microm%aggrselfprist(1,1,1)  &
                 ,ng, npts, imean,  &
                 'AGGRSELFPRIST :3:anal:mpti:mpt3')
   if (allocated(micro%aggrselfsnowt)) &
      CALL vtables2 (micro%aggrselfsnowt(1,1,1),microm%aggrselfsnowt(1,1,1)  &
                 ,ng, npts, imean,  &
                 'AGGRSELFSNOWT :3:anal:mpti:mpt3')
   if (allocated(micro%aggrprissnowt)) &
      CALL vtables2 (micro%aggrprissnowt(1,1,1),microm%aggrprissnowt(1,1,1)  &
                 ,ng, npts, imean,  &
                 'AGGRPRISSNOWT :3:anal:mpti:mpt3')

   if (allocated(micro%dust1cldrt)) &
        CALL vtables2 (micro%dust1cldrt(1,1,1),microm%dust1cldrt(1,1,1)  &
                       ,ng, npts, imean,  &
                       'DUST1CLDRT :3:anal:mpti:mpt3')
   if (allocated(micro%dust2cldrt)) &
        CALL vtables2 (micro%dust2cldrt(1,1,1),microm%dust2cldrt(1,1,1)  &
                       ,ng, npts, imean,  &
                       'DUST2CLDRT :3:anal:mpti:mpt3')
   if (allocated(micro%dust1drzrt)) &
        CALL vtables2 (micro%dust1drzrt(1,1,1),microm%dust1drzrt(1,1,1)  &
                       ,ng, npts, imean,  &
                       'DUST1DRZRT :3:anal:mpti:mpt3')
   if (allocated(micro%dust2drzrt)) &
        CALL vtables2 (micro%dust2drzrt(1,1,1),microm%dust2drzrt(1,1,1)  &
                       ,ng, npts, imean,  &
                       'DUST2DRZRT :3:anal:mpti:mpt3')

   npts=n2*n3
   if (allocated(micro%accpr)) &
      CALL vtables2 (micro%accpr(1,1),microm%accpr(1,1)  &
                 ,ng, npts, imean,  &
                 'ACCPR :2:anal:mpti:mpt3')
   if (allocated(micro%accpp)) &
      CALL vtables2 (micro%accpp(1,1),microm%accpp(1,1)  &
                 ,ng, npts, imean,  &
                 'ACCPP :2:anal:mpti:mpt3')
   if (allocated(micro%accps)) &
      CALL vtables2 (micro%accps(1,1),microm%accps(1,1)  &
                 ,ng, npts, imean,  &
                 'ACCPS :2:anal:mpti:mpt3')
   if (allocated(micro%accpa)) &
      CALL vtables2 (micro%accpa(1,1),microm%accpa(1,1)  &
                 ,ng, npts, imean,  &
                 'ACCPA :2:anal:mpti:mpt3')
   if (allocated(micro%accpg)) &
      CALL vtables2 (micro%accpg(1,1),microm%accpg(1,1)  &
                 ,ng, npts, imean,  &
                 'ACCPG :2:anal:mpti:mpt3')
   if (allocated(micro%accph)) &
      CALL vtables2 (micro%accph(1,1),microm%accph(1,1)  &
                 ,ng, npts, imean,  &
                 'ACCPH :2:anal:mpti:mpt3')
   if (allocated(micro%accpd)) &
      CALL vtables2 (micro%accpd(1,1),microm%accpd(1,1)  &
                 ,ng, npts, imean,  &
                 'ACCPD :2:anal:mpti:mpt3')
   if (allocated(micro%pcprr)) &
      CALL vtables2 (micro%pcprr(1,1),microm%pcprr(1,1)  &
                 ,ng, npts, imean,  &
                 'PCPRR :2:anal:mpt3')
   if (allocated(micro%pcprp)) &
      CALL vtables2 (micro%pcprp(1,1),microm%pcprp(1,1)  &
                 ,ng, npts, imean,  &
                 'PCPRP :2:anal:mpt3')
   if (allocated(micro%pcprs)) &
      CALL vtables2 (micro%pcprs(1,1),microm%pcprs(1,1)  &
                 ,ng, npts, imean,  &
                 'PCPRS :2:anal:mpt3')
   if (allocated(micro%pcpra)) &
      CALL vtables2 (micro%pcpra(1,1),microm%pcpra(1,1)  &
                 ,ng, npts, imean,  &
                 'PCPRA :2:anal:mpt3')
   if (allocated(micro%pcprg)) &
      CALL vtables2 (micro%pcprg(1,1),microm%pcprg(1,1)  &
                 ,ng, npts, imean,  &
                 'PCPRG :2:anal:mpt3')
   if (allocated(micro%pcprh)) &
      CALL vtables2 (micro%pcprh(1,1),microm%pcprh(1,1)  &
                 ,ng, npts, imean,  &
                 'PCPRH :2:anal:mpt3')
   if (allocated(micro%pcprd)) &
      CALL vtables2 (micro%pcprd(1,1),microm%pcprd(1,1)  &
                 ,ng, npts, imean,  &
                 'PCPRD :2:anal:mpt3')
   if (allocated(micro%pcpg)) &
      CALL vtables2 (micro%pcpg(1,1),microm%pcpg(1,1)  &
                 ,ng, npts, imean,  &
                 'PCPG :2:anal:mpti:mpt3')
   if (allocated(micro%qpcpg)) &
      CALL vtables2 (micro%qpcpg(1,1),microm%qpcpg(1,1)  &
                 ,ng, npts, imean,  &
                 'QPCPG :2:anal:mpti:mpt3')
   if (allocated(micro%dpcpg)) &
      CALL vtables2 (micro%dpcpg(1,1),microm%dpcpg(1,1)  &
                 ,ng, npts, imean,  &
                 'DPCPG :2:anal:mpti:mpt3')

!Total aerosol mass precipitation to the surface
   if (allocated(micro%accpaero)) &
      CALL vtables2 (micro%accpaero(1,1),microm%accpaero(1,1)  &
                 ,ng, npts, imean,  &
                 'ACCPAERO :2:anal:mpti:mpt3')
   if (allocated(micro%pcpraero)) &
      CALL vtables2 (micro%pcpraero(1,1),microm%pcpraero(1,1)  &
                 ,ng, npts, imean,  &
                 'PCPRAERO :2:anal:mpt3')

!Dust mass precipitation to the surface
   if (allocated(micro%accpdust)) &
      CALL vtables2 (micro%accpdust(1,1),microm%accpdust(1,1)  &
                 ,ng, npts, imean,  &
                 'ACCPDUST :2:anal:mpti:mpt3')
   if (allocated(micro%pcprdust)) &
      CALL vtables2 (micro%pcprdust(1,1),microm%pcprdust(1,1)  &
                 ,ng, npts, imean,  &
                 'PCPRDUST :2:anal:mpt3')

return
END SUBROUTINE filltab_micro

END MODULE mem_micro
