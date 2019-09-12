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

Subroutine aerosols ()

use mem_basic
use mem_micro
use mem_grid
use mem_leaf
use node_mod
use micphys

implicit none

integer :: i,j

!Run the SEASALT and DUST Source model before the call to Micro
if(idust==2) then
  CALL dust_sources (mzp,mxp,myp,ia,iz,ja,jz                  &
    ,leaf_g(ngrid)%leaf_class  ,grid_g(ngrid)%rtgm           &
    ,leaf_g(ngrid)%patch_area  ,leaf_g(ngrid)%veg_rough      &
    ,basic_g(ngrid)%up         ,basic_g(ngrid)%vp            &
    ,micro_g(ngrid)%md1np      ,micro_g(ngrid)%md2np         &
    ,micro_g(ngrid)%md1mp      ,micro_g(ngrid)%md2mp         &
    ,leaf_g(ngrid)%soil_water  ,leaf_g(ngrid)%soil_text      &
    ,grid_g(ngrid)%glat        ,grid_g(ngrid)%glon           &
    ,basic_g(ngrid)%dn0)
endif
if(isalt==2) then
  CALL salt_sources (mzp,mxp,myp,ia,iz,ja,jz                   &
    ,leaf_g(ngrid)%leaf_class  ,grid_g(ngrid)%rtgm            &
    ,leaf_g(ngrid)%patch_area                                 &
    ,basic_g(ngrid)%up           ,basic_g(ngrid)%vp           &
    ,micro_g(ngrid)%salt_film_np ,micro_g(ngrid)%salt_jet_np  &
    ,micro_g(ngrid)%salt_spum_np ,micro_g(ngrid)%salt_film_mp &
    ,micro_g(ngrid)%salt_jet_mp  ,micro_g(ngrid)%salt_spum_mp &
    ,basic_g(ngrid)%dn0)
endif

! Aerosol dry and wet deposition call when micro LEVEL < 3
if(iaerodep==1 .and. level<3) then
 do j = ja,jz
  do i = ia,iz

   CALL aero_copy (1,mzp &
    ,micro_g(ngrid)%cccnp(1,i,j),micro_g(ngrid)%cccmp(1,i,j) &
    ,micro_g(ngrid)%gccnp(1,i,j),micro_g(ngrid)%gccmp(1,i,j) &
    ,micro_g(ngrid)%md1np(1,i,j),micro_g(ngrid)%md1mp(1,i,j) &
    ,micro_g(ngrid)%md2np(1,i,j),micro_g(ngrid)%md2mp(1,i,j) &
    ,micro_g(ngrid)%salt_film_np(1,i,j),micro_g(ngrid)%salt_film_mp(1,i,j) &
    ,micro_g(ngrid)%salt_jet_np(1,i,j) ,micro_g(ngrid)%salt_jet_mp(1,i,j)  &
    ,micro_g(ngrid)%salt_spum_np(1,i,j),micro_g(ngrid)%salt_spum_mp(1,i,j))

   CALL deposition_driver (i,j,mzp,imonth1,zm &
    ,basic_g(ngrid)%rv(1,i,j) &
    ,basic_g(ngrid)%pi0(1,i,j) &
    ,basic_g(ngrid)%pp(1,i,j) &
    ,basic_g(ngrid)%theta(1,i,j) & 
    ,basic_g(ngrid)%up(1,i,j) &
    ,basic_g(ngrid)%vp(1,i,j) &
    ,basic_g(ngrid)%dn0(1,i,j) &
    ,leaf_g(ngrid)%ustar(i,j,1:npatch)      &
    ,leaf_g(ngrid)%leaf_class(i,j,1:npatch) &
    ,leaf_g(ngrid)%patch_area(i,j,1:npatch) &
    ,leaf_g(ngrid)%veg_rough(i,j,1:npatch)  &
    ,leaf_g(ngrid)%soil_rough(i,j,1:npatch) &
    ,grid_g(ngrid)%rtgt(i,j))

   CALL aero_copy (2,mzp &
    ,micro_g(ngrid)%cccnp(1,i,j),micro_g(ngrid)%cccmp(1,i,j) &
    ,micro_g(ngrid)%gccnp(1,i,j),micro_g(ngrid)%gccmp(1,i,j) &
    ,micro_g(ngrid)%md1np(1,i,j),micro_g(ngrid)%md1mp(1,i,j) &
    ,micro_g(ngrid)%md2np(1,i,j),micro_g(ngrid)%md2mp(1,i,j) &
    ,micro_g(ngrid)%salt_film_np(1,i,j),micro_g(ngrid)%salt_film_mp(1,i,j) &
    ,micro_g(ngrid)%salt_jet_np(1,i,j) ,micro_g(ngrid)%salt_jet_mp(1,i,j)  &
    ,micro_g(ngrid)%salt_spum_np(1,i,j),micro_g(ngrid)%salt_spum_mp(1,i,j))

  enddo
 enddo
endif

return
END SUBROUTINE aerosols

!##############################################################################
Subroutine aerosol_init ()

use micphys
use mem_grid, only:iprntstmt,print_msg

implicit none

real :: weightfac

! Set aerosol density depending on chemistry and soluble fraction
! Pure quantity densities (kg/m3) are:
! NH42S04 = 1769. (ammonium sulfate)
! Clay Dust (smaller) = 2500.
! Silt Dust (larger) = 2650.
! NaCl = 2165. (sodium chloride)

! Set Aerosol density (kg/m3) based on weighted mixture of soluble
! and insoluble material. Assume insoluble core to be like that of 
! silt dust with density = 2650 kg/m3, except for acat=3 which is
! already set to small sized clay dust
! Also set vanthoff factors for given chemistry

if(iprntstmt>=1 .and. print_msg) print*,''
if(iprntstmt>=1 .and. print_msg) print*,'Setting up default aerosol densities:'

do acat=1,aerocat
 
 aero_rhosol(acat)   = 0.0
 aero_vanthoff(acat) = 0.0

 if(acat==3) then !if small dust (clay)
   weightfac = 2500. * (1.0-aero_epsilon(acat)) !clay dust core
 else !all other
   weightfac = 2650. * (1.0-aero_epsilon(acat)) !silt dust core
 endif

 if(iaero_chem(acat)==1) then !NH42S04
   aero_rhosol(acat) = 1769. * aero_epsilon(acat) + weightfac
   aero_vanthoff(acat) = 3
 elseif(iaero_chem(acat)==2) then !NaCl
   aero_rhosol(acat) = 2165. * aero_epsilon(acat) + weightfac
   aero_vanthoff(acat) = 2
 endif

 if(iprntstmt>=1 .and. print_msg) print*,'acat,rg,rho,i:',acat &
     ,aero_medrad(acat),aero_rhosol(acat),aero_vanthoff(acat)

enddo

if(iprntstmt>=1 .and. print_msg) print*,''

return
END SUBROUTINE aerosol_init

!##############################################################################
Subroutine aero_copy (aflag,m1,cccnp,cccmp,gccnp,gccmp,md1np,md1mp &
                    ,md2np,md2mp,salt_film_np,salt_film_mp,salt_jet_np &
                    ,salt_jet_mp,salt_spum_np,salt_spum_mp)

!This routine is called in the event that MICRO LEVEL=1,2 so that
!aerosols can still be allowed to impact radiation.

use micphys

implicit none

integer :: m1,k,aflag
real, dimension(m1) :: cccnp,cccmp,gccnp,gccmp,md1np,md1mp &
                    ,md2np,md2mp,salt_film_np,salt_film_mp,salt_jet_np &
                    ,salt_jet_mp,salt_spum_np,salt_spum_mp

if(aflag==1)then
 !Zero out aerosol scratch arrays
 do acat = 1,aerocat
     do k = 1,m1
       aerocon(k,acat) = 0.0
       aeromas(k,acat) = 0.0
     enddo
 enddo
 !Fill scratch arrays for aerosol modes for level=1,2
 do k = 1,m1-1
   if (iaerosol > 0) then
     aerocon(k,1) = cccnp(k)
     aeromas(k,1) = cccmp(k)
     aerocon(k,2) = gccnp(k)
     aeromas(k,2) = gccmp(k)
   endif
   if (idust > 0) then
     aerocon(k,3) = md1np(k)
     aeromas(k,3) = md1mp(k)
     aerocon(k,4) = md2np(k)
     aeromas(k,4) = md2mp(k)
   endif
   if (isalt > 0) then
     aerocon(k,5) = salt_film_np(k)
     aeromas(k,5) = salt_film_mp(k)
     aerocon(k,6) = salt_jet_np(k)
     aeromas(k,6) = salt_jet_mp(k)
     aerocon(k,7) = salt_spum_np(k)
     aeromas(k,7) = salt_spum_mp(k)
   endif
 enddo

elseif(aflag==2)then
 !Copy back scratch arrays to aerosol modes for level=1,2
 do k = 1,m1-1
   if (iaerosol > 0) then
    cccnp(k) = aerocon(k,1)
    cccmp(k) = aeromas(k,1)
    gccnp(k) = aerocon(k,2)
    gccmp(k) = aeromas(k,2)
   endif
   if (idust > 0) then
    md1np(k) = aerocon(k,3)
    md1mp(k) = aeromas(k,3)
    md2np(k) = aerocon(k,4)
    md2mp(k) = aeromas(k,4)
   endif
   if (isalt > 0) then
    salt_film_np(k) = aerocon(k,5)
    salt_film_mp(k) = aeromas(k,5)
    salt_jet_np(k)  = aerocon(k,6)
    salt_jet_mp(k)  = aeromas(k,6)
    salt_spum_np(k) = aerocon(k,7)
    salt_spum_mp(k) = aeromas(k,7)
   endif
 enddo
endif

return
END SUBROUTINE aero_copy

!##############################################################################
Subroutine checkmicro ()

use mem_basic
use mem_micro
use mem_grid
use node_mod
use micphys
use rconstants

implicit none

integer :: k,i,j,prtflg,mprtflg
logical, external :: isnanr

prtflg=0
mprtflg=0

do j = ja,jz
do i = ia,iz
do k = 1,mzp

   !CHECK FOR NEGATIVE VAPOR OR MIXING RATIO
   if(level>=1)then
     if(basic_g(ngrid)%rv(k,i,j)<0.0.or.basic_g(ngrid)%rtp(k,i,j)<0.0) mprtflg=1
   endif
   if(level>=2)then
     if(micro_g(ngrid)%rcp(k,i,j)<0.0) mprtflg=1
   endif
   if(level==3)then
     if((idriz  >= 1 .and. micro_g(ngrid)%rdp(k,i,j)<0.0) .or. &
        (irain  >= 1 .and. micro_g(ngrid)%rrp(k,i,j)<0.0) .or. &
        (ipris  >= 1 .and. micro_g(ngrid)%rpp(k,i,j)<0.0) .or. &
        (isnow  >= 1 .and. micro_g(ngrid)%rsp(k,i,j)<0.0) .or. &
        (iaggr  >= 1 .and. micro_g(ngrid)%rap(k,i,j)<0.0) .or. &
        (igraup >= 1 .and. micro_g(ngrid)%rgp(k,i,j)<0.0) .or. &
        (ihail  >= 1 .and. micro_g(ngrid)%rhp(k,i,j)<0.0)) mprtflg=1
   endif

   !CHECK DUST MODES
   if(idust > 0)then
     if(isnanr(micro_g(ngrid)%md1np(k,i,j)) .or. &
        isnanr(micro_g(ngrid)%md1mp(k,i,j)) .or. &
        isnanr(micro_g(ngrid)%md2np(k,i,j)) .or. &
        isnanr(micro_g(ngrid)%md2mp(k,i,j)) ) prtflg=1
   endif
   !CHECK SEA SALT MODES
   if(isalt > 0)then
     if(isnanr(micro_g(ngrid)%salt_film_np(k,i,j)) .or. &
        isnanr(micro_g(ngrid)%salt_jet_np(k,i,j))  .or. &
        isnanr(micro_g(ngrid)%salt_spum_np(k,i,j)) .or. &
        isnanr(micro_g(ngrid)%salt_film_mp(k,i,j)) .or. &
        isnanr(micro_g(ngrid)%salt_jet_mp(k,i,j))  .or. &
        isnanr(micro_g(ngrid)%salt_spum_mp(k,i,j)) ) prtflg=1
   endif
   !CHECK CCN AND GCCN
   if(iaerosol > 0)then
     if(isnanr(micro_g(ngrid)%cccnp(k,i,j)) .or. &
        isnanr(micro_g(ngrid)%cccmp(k,i,j)) .or. &
        isnanr(micro_g(ngrid)%gccnp(k,i,j)) .or. &
        isnanr(micro_g(ngrid)%gccmp(k,i,j)) ) prtflg=1
   endif
   !CHECK REGENERATED AEROSOL MODES
   if(iccnlev>=2)then
     if(isnanr(micro_g(ngrid)%regen_aero1_np(k,i,j)) .or. &
        isnanr(micro_g(ngrid)%regen_aero2_np(k,i,j)) .or. &
        isnanr(micro_g(ngrid)%regen_aero1_mp(k,i,j)) .or. &
        isnanr(micro_g(ngrid)%regen_aero2_mp(k,i,j)) ) prtflg=1
   endif

   !CHECK TOTAL AEROSOL TRACKING VARIABLES
   if(iccnlev>=2)then
     if(icloud>=1)then 
       if(isnanr(micro_g(ngrid)%cnmcp(k,i,j))) prtflg=1 
     endif
     if(irain>=1)then
       if(isnanr(micro_g(ngrid)%cnmrp(k,i,j))) prtflg=1
     endif
     if(ipris>=1)then
       if(isnanr(micro_g(ngrid)%cnmpp(k,i,j))) prtflg=1
     endif
     if(isnow>=1)then
       if(isnanr(micro_g(ngrid)%cnmsp(k,i,j))) prtflg=1
     endif
     if(iaggr>=1)then
       if(isnanr(micro_g(ngrid)%cnmap(k,i,j))) prtflg=1
     endif
     if(igraup>=1)then
       if(isnanr(micro_g(ngrid)%cnmgp(k,i,j))) prtflg=1
     endif
     if(ihail>=1)then
       if(isnanr(micro_g(ngrid)%cnmhp(k,i,j))) prtflg=1
     endif
     if(idriz>=1)then
       if(isnanr(micro_g(ngrid)%cnmdp(k,i,j))) prtflg=1
     endif
     if(itrkepsilon==1) then
       if(icloud>=1)then
         if(isnanr(micro_g(ngrid)%snmcp(k,i,j))) prtflg=1
       endif
       if(irain>=1)then
         if(isnanr(micro_g(ngrid)%snmrp(k,i,j))) prtflg=1
       endif
       if(ipris>=1)then
         if(isnanr(micro_g(ngrid)%snmpp(k,i,j))) prtflg=1
       endif
       if(isnow>=1)then
         if(isnanr(micro_g(ngrid)%snmsp(k,i,j))) prtflg=1
       endif
       if(iaggr>=1)then
         if(isnanr(micro_g(ngrid)%snmap(k,i,j))) prtflg=1
       endif
       if(igraup>=1)then
         if(isnanr(micro_g(ngrid)%snmgp(k,i,j))) prtflg=1
       endif
       if(ihail>=1)then
         if(isnanr(micro_g(ngrid)%snmhp(k,i,j))) prtflg=1
       endif
       if(idriz>=1)then
         if(isnanr(micro_g(ngrid)%snmdp(k,i,j))) prtflg=1
       endif
       if(isnanr(micro_g(ngrid)%resol_aero1_mp(k,i,j)) .or. &
          isnanr(micro_g(ngrid)%resol_aero2_mp(k,i,j)) ) prtflg=1
     endif
     if(itrkdust==1)then
       if(icloud>=1)then
         if(isnanr(micro_g(ngrid)%dnmcp(k,i,j))) prtflg=1
       endif
       if(irain>=1)then
         if(isnanr(micro_g(ngrid)%dnmrp(k,i,j))) prtflg=1
       endif
       if(ipris>=1)then
         if(isnanr(micro_g(ngrid)%dnmpp(k,i,j))) prtflg=1
       endif
       if(isnow>=1)then
         if(isnanr(micro_g(ngrid)%dnmsp(k,i,j))) prtflg=1
       endif
       if(iaggr>=1)then
         if(isnanr(micro_g(ngrid)%dnmap(k,i,j))) prtflg=1
       endif
       if(igraup>=1)then
         if(isnanr(micro_g(ngrid)%dnmgp(k,i,j))) prtflg=1
       endif
       if(ihail>=1)then
         if(isnanr(micro_g(ngrid)%dnmhp(k,i,j))) prtflg=1
       endif
       if(idriz>=1)then
         if(isnanr(micro_g(ngrid)%dnmdp(k,i,j))) prtflg=1
       endif
     endif
     if(itrkdustifn==1)then
       if(icloud>=1)then
         if(isnanr(micro_g(ngrid)%dincp(k,i,j))) prtflg=1
       endif
       if(irain>=1)then
         if(isnanr(micro_g(ngrid)%dinrp(k,i,j))) prtflg=1
       endif
       if(ipris>=1)then
         if(isnanr(micro_g(ngrid)%dinpp(k,i,j))) prtflg=1
       endif
       if(isnow>=1)then
         if(isnanr(micro_g(ngrid)%dinsp(k,i,j))) prtflg=1
       endif
       if(iaggr>=1)then
         if(isnanr(micro_g(ngrid)%dinap(k,i,j))) prtflg=1
       endif
       if(igraup>=1)then
         if(isnanr(micro_g(ngrid)%dingp(k,i,j))) prtflg=1
       endif
       if(ihail>=1)then
         if(isnanr(micro_g(ngrid)%dinhp(k,i,j))) prtflg=1
       endif
       if(idriz>=1)then
         if(isnanr(micro_g(ngrid)%dindp(k,i,j))) prtflg=1
       endif
     endif
   endif

   !CHECK IMMERSION FREEZING NUCLEI TRACKING VARIABLES
   if(iifn==3 .and. iccnlev>=1)then
     if(jnmb(1) >= 5)then
       if(isnanr(micro_g(ngrid)%ifnnucp(k,i,j))) prtflg=1
     endif
     if(jnmb(1) >= 5)then
       if(isnanr(micro_g(ngrid)%immercp(k,i,j))) prtflg=1
     endif
     if(jnmb(8) >= 5)then
       if(isnanr(micro_g(ngrid)%immerdp(k,i,j))) prtflg=1
     endif
     if(jnmb(2) >= 5)then
       if(isnanr(micro_g(ngrid)%immerrp(k,i,j))) prtflg=1
     endif
   endif

   !CHECK HYDROMETEOR NUMBER CONCENTRATIONS
   if(jnmb(1) >= 5)then
     if(isnanr(micro_g(ngrid)%ccp(k,i,j))) prtflg=1
   endif
   if(jnmb(8) >= 5)then
     if(isnanr(micro_g(ngrid)%cdp(k,i,j))) prtflg=1
   endif
   if(jnmb(2) >= 5)then
     if(isnanr(micro_g(ngrid)%crp(k,i,j))) prtflg=1
   endif
   if(jnmb(3) >= 5)then
     if(isnanr(micro_g(ngrid)%cpp(k,i,j))) prtflg=1
   endif
   if(jnmb(4) >= 5)then
     if(isnanr(micro_g(ngrid)%csp(k,i,j))) prtflg=1
   endif
   if(jnmb(5) >= 5)then
     if(isnanr(micro_g(ngrid)%cap(k,i,j))) prtflg=1
   endif
   if(jnmb(6) >= 5)then
     if(isnanr(micro_g(ngrid)%cgp(k,i,j))) prtflg=1
   endif
   if(jnmb(7) >= 5)then
     if(isnanr(micro_g(ngrid)%chp(k,i,j))) prtflg=1
   endif
   !CHECK HYDROMETEOR MIXING RATIOS
   if(level>=2)then
     if(isnanr(micro_g(ngrid)%rcp(k,i,j))) prtflg=1
   endif
   if(idriz>=1)then
     if(isnanr(micro_g(ngrid)%rdp(k,i,j))) prtflg=1
   endif
   if(irain>=1)then
     if(isnanr(micro_g(ngrid)%rrp(k,i,j))) prtflg=1
   endif   
   if(ipris>=1)then
     if(isnanr(micro_g(ngrid)%rpp(k,i,j))) prtflg=1
   endif
   if(isnow>=1)then
     if(isnanr(micro_g(ngrid)%rsp(k,i,j))) prtflg=1
   endif
   if(iaggr>=1)then
     if(isnanr(micro_g(ngrid)%rap(k,i,j))) prtflg=1
   endif
   if(igraup>=1)then
     if(isnanr(micro_g(ngrid)%rgp(k,i,j))) prtflg=1
   endif
   if(ihail>=1)then
     if(isnanr(micro_g(ngrid)%rhp(k,i,j))) prtflg=1
   endif
   !CHECK 3D PRECIPITATION RATES
   if(idriz>=1)then
     if(isnanr(micro_g(ngrid)%pcpvd(k,i,j))) prtflg=1
   endif
   if(irain>=1)then
     if(isnanr(micro_g(ngrid)%pcpvr(k,i,j))) prtflg=1
   endif   
   if(ipris>=1)then
     if(isnanr(micro_g(ngrid)%pcpvp(k,i,j))) prtflg=1
   endif
   if(isnow>=1)then
     if(isnanr(micro_g(ngrid)%pcpvs(k,i,j))) prtflg=1
   endif
   if(iaggr>=1)then
     if(isnanr(micro_g(ngrid)%pcpva(k,i,j))) prtflg=1
   endif
   if(igraup>=1)then
     if(isnanr(micro_g(ngrid)%pcpvg(k,i,j))) prtflg=1
   endif
   if(ihail>=1)then
     if(isnanr(micro_g(ngrid)%pcpvh(k,i,j))) prtflg=1
   endif

 !*******************************************************************
 !IF NEGATIVE MIXING RATIO or NAN EXISTS THEN PRINT MICRO INFORMATION
 !*******************************************************************
 if(mprtflg==1)then
    print*,'Negative Condensate MICRO (ngrid,k,i,j):',ngrid,k,i,j
    if(level>=1)then
      print*,'vapor,rtp:',basic_g(ngrid)%rv(k,i,j),basic_g(ngrid)%rtp(k,i,j)
    endif
    if(level>=2)then 
      print*,'cloud:  ',micro_g(ngrid)%rcp(k,i,j)
    endif
    if(level==3)then
      if(idriz  >= 1) print*,'driz:   ',micro_g(ngrid)%rdp(k,i,j)
      if(irain  >= 1) print*,'rain:   ',micro_g(ngrid)%rrp(k,i,j)
      if(ipris  >= 1) print*,'ice:    ',micro_g(ngrid)%rpp(k,i,j)
      if(isnow  >= 1) print*,'snow:   ',micro_g(ngrid)%rsp(k,i,j)
      if(iaggr  >= 1) print*,'aggr:   ',micro_g(ngrid)%rap(k,i,j)
      if(igraup >= 1) print*,'graup:  ',micro_g(ngrid)%rgp(k,i,j)
      if(ihail  >= 1) print*,'hail:   ',micro_g(ngrid)%rhp(k,i,j)
    endif
    print*,'Try shorter timestep. These can become negative due to'
    print*,' diffusion, advection, cumulus parameterization, etc.'
 endif

 if(prtflg==1)then
    print*,'NAN k,i,j',k,i+i0,j+j0,ngrid
    if(iaerosol > 0)then
      print*,'cccnp',micro_g(ngrid)%cccnp(k,i,j)
      print*,'cccmp',micro_g(ngrid)%cccmp(k,i,j)
      print*,'gccnp',micro_g(ngrid)%gccnp(k,i,j)
      print*,'gccmp',micro_g(ngrid)%gccmp(k,i,j)
    endif
    if(idust > 0)then
      print*,'md1np',micro_g(ngrid)%md1np(k,i,j)
      print*,'md1mp',micro_g(ngrid)%md1mp(k,i,j)
      print*,'md2np',micro_g(ngrid)%md2np(k,i,j)
      print*,'md2mp',micro_g(ngrid)%md2mp(k,i,j)
    endif
    if(isalt > 0)then
      print*,'salt_film_np',micro_g(ngrid)%salt_film_np(k,i,j)
      print*,'salt_jet_np' ,micro_g(ngrid)%salt_jet_np(k,i,j)
      print*,'salt_spum_np',micro_g(ngrid)%salt_spum_np(k,i,j)
      print*,'salt_film_mp',micro_g(ngrid)%salt_film_mp(k,i,j)
      print*,'salt_jet_mp' ,micro_g(ngrid)%salt_jet_mp(k,i,j)
      print*,'salt_spum_mp',micro_g(ngrid)%salt_spum_mp(k,i,j)
    endif
    if(iccnlev>=2)then
      print*,'regen1',micro_g(ngrid)%regen_aero1_np(k,i,j)
      print*,'regen2',micro_g(ngrid)%regen_aero2_np(k,i,j)
      print*,'regem1',micro_g(ngrid)%regen_aero1_mp(k,i,j)
      print*,'regem2',micro_g(ngrid)%regen_aero2_mp(k,i,j)
      if(icloud>=1)print*,'cnmcp',micro_g(ngrid)%cnmcp(k,i,j)
      if(irain>=1) print*,'cnmrp',micro_g(ngrid)%cnmrp(k,i,j)
      if(ipris>=1) print*,'cnmpp',micro_g(ngrid)%cnmpp(k,i,j)
      if(isnow>=1) print*,'cnmsp',micro_g(ngrid)%cnmsp(k,i,j)
      if(iaggr>=1) print*,'cnmap',micro_g(ngrid)%cnmap(k,i,j)
      if(igraup>=1)print*,'cnmgp',micro_g(ngrid)%cnmgp(k,i,j)
      if(ihail>=1) print*,'cnmhp',micro_g(ngrid)%cnmhp(k,i,j)
      if(idriz>=1) print*,'cnmdp',micro_g(ngrid)%cnmdp(k,i,j)
      if(itrkdust==1)then
       if(icloud>=1)print*,'dnmcp',micro_g(ngrid)%dnmcp(k,i,j)
       if(irain>=1) print*,'dnmrp',micro_g(ngrid)%dnmrp(k,i,j)
       if(ipris>=1) print*,'dnmpp',micro_g(ngrid)%dnmpp(k,i,j)
       if(isnow>=1) print*,'dnmsp',micro_g(ngrid)%dnmsp(k,i,j)
       if(iaggr>=1) print*,'dnmap',micro_g(ngrid)%dnmap(k,i,j)
       if(igraup>=1)print*,'dnmgp',micro_g(ngrid)%dnmgp(k,i,j)
       if(ihail>=1) print*,'dnmhp',micro_g(ngrid)%dnmhp(k,i,j)
       if(idriz>=1) print*,'dnmdp',micro_g(ngrid)%dnmdp(k,i,j)
      endif
      if(itrkdustifn==1)then
       if(icloud>=1)print*,'dincp',micro_g(ngrid)%dincp(k,i,j)
       if(irain>=1) print*,'dinrp',micro_g(ngrid)%dinrp(k,i,j)
       if(ipris>=1) print*,'dinpp',micro_g(ngrid)%dinpp(k,i,j)
       if(isnow>=1) print*,'dinsp',micro_g(ngrid)%dinsp(k,i,j)
       if(iaggr>=1) print*,'dinap',micro_g(ngrid)%dinap(k,i,j)
       if(igraup>=1)print*,'dingp',micro_g(ngrid)%dingp(k,i,j)
       if(ihail>=1) print*,'dinhp',micro_g(ngrid)%dinhp(k,i,j)
       if(idriz>=1) print*,'dindp',micro_g(ngrid)%dindp(k,i,j)
      endif
      if(itrkepsilon==1)then
       print*,'regemsol1',micro_g(ngrid)%resol_aero1_mp(k,i,j)
       print*,'regemsol2',micro_g(ngrid)%resol_aero2_mp(k,i,j)
       if(icloud>=1)print*,'snmcp',micro_g(ngrid)%snmcp(k,i,j)
       if(irain>=1) print*,'snmrp',micro_g(ngrid)%snmrp(k,i,j)
       if(ipris>=1) print*,'snmpp',micro_g(ngrid)%snmpp(k,i,j)
       if(isnow>=1) print*,'snmsp',micro_g(ngrid)%snmsp(k,i,j)
       if(iaggr>=1) print*,'snmap',micro_g(ngrid)%snmap(k,i,j)
       if(igraup>=1)print*,'snmgp',micro_g(ngrid)%snmgp(k,i,j)
       if(ihail>=1) print*,'snmhp',micro_g(ngrid)%snmhp(k,i,j)
       if(idriz>=1) print*,'snmdp',micro_g(ngrid)%snmdp(k,i,j)
      endif
    endif
    if(iifn==3 .and. iccnlev>=1)then
      if(icloud>=5) print*,'ifnnucp',micro_g(ngrid)%ifnnucp(k,i,j)
      if(icloud>=5) print*,'immercp',micro_g(ngrid)%immercp(k,i,j)
      if(idriz>=5)  print*,'immerdp',micro_g(ngrid)%immerdp(k,i,j)
      if(irain>=5)  print*,'immerrp',micro_g(ngrid)%immerrp(k,i,j)
    endif
    if(icloud>=5)print*,'ccp',micro_g(ngrid)%ccp(k,i,j)
    if(idriz>=5) print*,'cdp',micro_g(ngrid)%cdp(k,i,j)
    if(irain>=5) print*,'crp',micro_g(ngrid)%crp(k,i,j)
    if(ipris>=5) print*,'cpp',micro_g(ngrid)%cpp(k,i,j)
    if(isnow>=5) print*,'csp',micro_g(ngrid)%csp(k,i,j)
    if(iaggr>=5) print*,'cap',micro_g(ngrid)%cap(k,i,j)
    if(igraup>=5)print*,'cgp',micro_g(ngrid)%cgp(k,i,j)
    if(ihail>=5) print*,'chp',micro_g(ngrid)%chp(k,i,j)
    if(icloud>=1)print*,'rcp',micro_g(ngrid)%rcp(k,i,j)
    if(idriz>=1) print*,'rdp',micro_g(ngrid)%rdp(k,i,j)
    if(irain>=1) print*,'rrp',micro_g(ngrid)%rrp(k,i,j)
    if(ipris>=1) print*,'rpp',micro_g(ngrid)%rpp(k,i,j)
    if(isnow>=1) print*,'rsp',micro_g(ngrid)%rsp(k,i,j)
    if(iaggr>=1) print*,'rap',micro_g(ngrid)%rap(k,i,j)
    if(igraup>=1)print*,'rgp',micro_g(ngrid)%rgp(k,i,j)
    if(ihail>=1) print*,'rhp',micro_g(ngrid)%rhp(k,i,j)
    if(idriz>=1) print*,'pcpvd',micro_g(ngrid)%pcpvd(k,i,j)
    if(irain>=1) print*,'pcpvr',micro_g(ngrid)%pcpvr(k,i,j)
    if(ipris>=1) print*,'pcpvp',micro_g(ngrid)%pcpvp(k,i,j)
    if(isnow>=1) print*,'pcpvs',micro_g(ngrid)%pcpvs(k,i,j)
    if(iaggr>=1) print*,'pcpva',micro_g(ngrid)%pcpva(k,i,j)
    if(igraup>=1)print*,'pcpvg',micro_g(ngrid)%pcpvg(k,i,j)
    if(ihail>=1) print*,'pcpvh',micro_g(ngrid)%pcpvh(k,i,j)
 endif

 if(mprtflg==1 .or. prtflg==1) stop

enddo
enddo
enddo

return
END SUBROUTINE checkmicro

!##############################################################################
logical Function isnanr (x)

implicit none

real :: x

isnanr = .false. 
if( (x .ne. x) .or. (x .ne. 0.0 .and. x/x .ne. 1) .or. (x*0.0 .ne. 0.0))then
   isnanr = .true.
endif
    
return
END FUNCTION isnanr

