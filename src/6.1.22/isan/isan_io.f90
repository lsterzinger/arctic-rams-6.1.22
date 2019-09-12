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

Subroutine isenio (inout,n1,n2)

use isan_coms
use hdf5_utils

implicit none

integer :: n1,n2,npts,nx3,ny3,ninn,l,levnn(maxisn),ndims,idims(4)
character(len=*) :: inout

if(inout == 'IN') THEN

   ndims=1 ; idims(1)=1
   CALL shdf5_irec ('isen_year',ivars=iyy)
   CALL shdf5_irec ('isen_month',ivars=imm)
   CALL shdf5_irec ('isen_date',ivars=idd)
   CALL shdf5_irec ('isen_hour',ivars=ihh)
   CALL shdf5_irec ('isen_nx',ivars=nx3)
   CALL shdf5_irec ('isen_ny',ivars=ny3)
   CALL shdf5_irec ('isen_nisn',ivars=ninn)
   ndims=1 ; idims(1)=ninn
   CALL shdf5_irec ('isen_levth',ivara=levnn)

   if(nx3.ne.n1.or.ny3.ne.n2.or.ninn.ne.nisn) then
      print*,'Isentropic stage grid dimensions do not match'
      print*,'   configuration file on read !'
      print*,' File dimens - ',nx3,ny3,ninn
      print*,' Run  dimens - ',n1,n2,nisn
      stop 'IO3-2'
   endif

   npts=n1*n2*nisn
   ndims=3 ; idims(1)=n1 ; idims(2)=n2 ; idims(3)=nisn
   CALL shdf5_irec ('isen_u',rvara=pi_u)
   CALL vmissr (pi_u,npts,1e30,-9998.)
   CALL shdf5_irec ('isen_v',rvara=pi_v)
   CALL vmissr (pi_v,npts,1e30,-9998.)
   CALL shdf5_irec ('isen_s',rvara=pi_s)
   CALL vmissr (pi_p,npts,1e30,-9998.)
   CALL shdf5_irec ('isen_p',rvara=pi_p)
   CALL vmissr (pi_s,npts,1e30,-9998.)
   CALL shdf5_irec ('isen_r',rvara=pi_r)
   CALL vmissr (pi_r,npts,1e30,-9998.)

   npts=n1*n2
   ndims=2 ; idims(1)=n1 ; idims(2)=n2
   CALL shdf5_irec ('sfc_u',rvara=rs_u)
   CALL vmissr (rs_u,npts,1e30,-9998.)
   CALL shdf5_irec ('sfc_v',rvara=rs_v)
   CALL vmissr (rs_v,npts,1e30,-9998.)
   CALL shdf5_irec ('sfc_p',rvara=rs_p)
   CALL vmissr (rs_p,npts,1e30,-9998.)
   CALL shdf5_irec ('sfc_t',rvara=rs_t)
   CALL vmissr (rs_t,npts,1e30,-9998.)
   CALL shdf5_irec ('sfc_r',rvara=rs_r)
   CALL vmissr (rs_r,npts,1e30,-9998.)
   CALL shdf5_irec ('sfc_s',rvara=rs_s)
   CALL vmissr (rs_s,npts,1e30,-9998.)
   CALL shdf5_irec ('sfc_topo',rvara=rs_top)
   CALL vmissr (rs_top,npts,1e30,-9998.)
   CALL shdf5_irec ('sfc_qual',rvara=rs_qual)
   CALL vmissr (rs_qual,npts,1e30,-9998.)

   CALL shdf5_irec ('sfc_soilmoist1',rvara=rs_soilmoist1)
   CALL vmissr (rs_soilmoist1,npts,1e30,-9998.)
   CALL shdf5_irec ('sfc_soilmoist2',rvara=rs_soilmoist2)
   CALL vmissr (rs_soilmoist2,npts,1e30,-9998.)
   CALL shdf5_irec ('sfc_soiltemp1',rvara=rs_soiltemp1)
   CALL vmissr (rs_soiltemp1,npts,1e30,-9998.)
   CALL shdf5_irec ('sfc_soiltemp2',rvara=rs_soiltemp2)
   CALL vmissr (rs_soiltemp2,npts,1e30,-9998.)
   CALL shdf5_irec ('sfc_snowmass',rvara=rs_snowmass)
   CALL vmissr (rs_snowmass,npts,1e30,-9998.)
   CALL shdf5_irec ('sfc_snowdepth',rvara=rs_snowdepth)
   CALL vmissr (rs_snowdepth,npts,1e30,-9998.)

   print 201,' *****  Isentropic file input *****************'  &
        ,iyear,imonth,idate,ihour,n1,n2,nisn  &
        ,(levth(l),l=1,nisn)
   201 format(//,a,//  &
        ,' *',7X,' Date (year,month,day,hour)  - ',4I5,/  &
        ,' *',7X,' Number of X,Y points        - ',2I5,/  &
        ,' *',7X,' Number of isentropic levels - ',I5,/  &
        ,' *',7X,' Isentropic levels (K)       - '/,(32X,8I5))
   print '(a)',' **********************************************'

endif

if(inout == 'OUT') then

   ndims=1 ; idims(1)=1
   CALL shdf5_orec (ndims,idims,'isen_year',ivars=iyear)
   CALL shdf5_orec (ndims,idims,'isen_month',ivars=imonth)
   CALL shdf5_orec (ndims,idims,'isen_date',ivars=idate)
   CALL shdf5_orec (ndims,idims,'isen_hour',ivars=ihour)
   CALL shdf5_orec (ndims,idims,'isen_nx',ivars=n1)
   CALL shdf5_orec (ndims,idims,'isen_ny',ivars=n2)
   CALL shdf5_orec (ndims,idims,'isen_nisn',ivars=nisn)
   ndims=1 ; idims(1)=nisn
   CALL shdf5_orec (ndims,idims,'isen_levth',ivara=levth)

   npts=n1*n2*nisn
   ndims=3 ; idims(1)=n1 ; idims(2)=n2 ; idims(3)=nisn
   CALL vmissw (pi_u,npts,pi_scra,1E30,-9999.)
   CALL shdf5_orec (ndims,idims,'isen_u',rvara=pi_scra)
   CALL vmissw (pi_v,npts,pi_scra,1E30,-9999.)
   CALL shdf5_orec (ndims,idims,'isen_v',rvara=pi_scra)
   CALL vmissw (pi_p,npts,pi_scra,1E30,-9999.)
   CALL shdf5_orec (ndims,idims,'isen_p',rvara=pi_scra)
   CALL vmissw (pi_s,npts,pi_scra,1E30,-9999.)
   CALL shdf5_orec (ndims,idims,'isen_s',rvara=pi_scra)
   CALL vmissw (pi_r,npts,pi_scra,1E30,-9999.)
   CALL shdf5_orec (ndims,idims,'isen_r',rvara=pi_scra)

   npts=n1*n2
   ndims=2 ; idims(1)=n1 ; idims(2)=n2
   CALL vmissw (rs_u,npts,pi_scra,1E30,-9999.)
   CALL shdf5_orec (ndims,idims,'sfc_u',rvara=pi_scra)
   CALL vmissw (rs_v,npts,pi_scra,1E30,-9999.)
   CALL shdf5_orec (ndims,idims,'sfc_v',rvara=pi_scra)
   CALL vmissw (rs_p,npts,pi_scra,1E30,-9999.)
   CALL shdf5_orec (ndims,idims,'sfc_p',rvara=pi_scra)
   CALL vmissw (rs_t,npts,pi_scra,1E30,-9999.)
   CALL shdf5_orec (ndims,idims,'sfc_t',rvara=pi_scra)
   CALL vmissw (rs_r,npts,pi_scra,1E30,-9999.)
   CALL shdf5_orec (ndims,idims,'sfc_r',rvara=pi_scra)
   CALL vmissw (rs_s,npts,pi_scra,1E30,-9999.)
   CALL shdf5_orec (ndims,idims,'sfc_s',rvara=pi_scra)
   CALL vmissw (rs_top,npts,pi_scra,1E30,-9999.)
   CALL shdf5_orec (ndims,idims,'sfc_topo',rvara=pi_scra)
   CALL vmissw (rs_qual,npts,pi_scra,1E30,-9999.)
   CALL shdf5_orec (ndims,idims,'sfc_qual',rvara=pi_scra)
   
   CALL vmissw (rs_soilmoist1,npts,pi_scra,1E30,-9999.)
   CALL shdf5_orec (ndims,idims,'sfc_soilmoist1',rvara=pi_scra)
   CALL vmissw (rs_soilmoist2,npts,pi_scra,1E30,-9999.)
   CALL shdf5_orec (ndims,idims,'sfc_soilmoist2',rvara=pi_scra)
   CALL vmissw (rs_soiltemp1,npts,pi_scra,1E30,-9999.)
   CALL shdf5_orec (ndims,idims,'sfc_soiltemp1',rvara=pi_scra)
   CALL vmissw (rs_soiltemp2,npts,pi_scra,1E30,-9999.)
   CALL shdf5_orec (ndims,idims,'sfc_soiltemp2',rvara=pi_scra)
   CALL vmissw (rs_snowmass,npts,pi_scra,1E30,-9999.)
   CALL shdf5_orec (ndims,idims,'sfc_snowmass',rvara=pi_scra)
   CALL vmissw (rs_snowdepth,npts,pi_scra,1E30,-9999.)
   CALL shdf5_orec (ndims,idims,'sfc_snowdepth',rvara=pi_scra)

   print 201,' *****  Isentropic file written *************'  &
        ,iyear,imonth,idate,ihour,n1,n2,nisn  &
        ,(levth(l),l=1,nisn)

   print 303,igridfl,gobsep,gobrad
   303 format(/,  &
         ' Grid flag (IGRIDFL)               -',I4,/  &
        ,' Grid-obs separation in degrees    -',F5.2,/  &
        ,' Grid-obs radius influence degrees -',F5.2)

endif

return
END SUBROUTINE isenio

!##############################################################################
Subroutine sigzio (inout,n1,n2)

use isan_coms
use hdf5_utils

implicit none

integer :: n1,n2,npts,l,ninn,nx3,ny3,ndims,idims(4)
character(len=*) :: inout

if(inout == 'IN') then

   ndims=1 ; idims(1)=1
   CALL shdf5_irec ('isen_year',ivars=iyy)
   CALL shdf5_irec ('isen_month',ivars=imm)
   CALL shdf5_irec ('isen_date',ivars=idd)
   CALL shdf5_irec ('isen_hour',ivars=ihh)
   CALL shdf5_irec ('isen_nx',ivars=nx3)
   CALL shdf5_irec ('isen_ny',ivars=ny3)
   CALL shdf5_irec ('sigz_nsigz',ivars=ninn)
   ndims=1 ; idims(1)=ninn
   CALL shdf5_irec ('sigz_sigz',rvara=sigz)

   if(nx3.ne.n1.or.ny3.ne.n2.or.ninn.ne.nsigz)then
      print*,'Sigma-z grid dimensions do not match'
      print*,'   input data on read !'
      print*,' File  dimensions - ',nx3,ny3,ninn
      print*,' Input dimensions - ',n1,n2,nsigz
      stop 'iO3-2'
   endif

   npts=n1*n2*nsigz
   ndims=3 ; idims(1)=n1 ; idims(2)=n2 ; idims(3)=nsigz
   CALL shdf5_irec ('sigz_u',rvara=ps_u)
   CALL vmissr (ps_u,npts,1e30,-9998.)
   CALL shdf5_irec ('sigz_v',rvara=ps_v)
   CALL vmissr (ps_v,npts,1e30,-9998.)
   CALL shdf5_irec ('sigz_p',rvara=ps_p)
   CALL vmissr (ps_p,npts,1e30,-9998.)
   CALL shdf5_irec ('sigz_t',rvara=ps_t)
   CALL vmissr (ps_t,npts,1e30,-9998.)
   CALL shdf5_irec ('sigz_r',rvara=ps_r)
   CALL vmissr (ps_r,npts,1e30,-9998.)

   print 201,' *****  Sigma-z file input *****************'  &
        ,iyear,imonth,idate,ihour,n1,n2,nsigz  &
        ,(sigz(l),l=1,nsigz)
   201 format(//,a,//  &
        ,' *',7X,' Date (year,month,day,hour)  - ',4I5,/  &
        ,' *',7X,' Number of X,Y points        - ',2I5,/  &
        ,' *',7X,' Number of sigma-z levels    - ',I5,/  &
        ,' *',7X,' Sigma-z levels (m)          - '/,(32X,7F8.1))
   print '(a)',' **********************************************'

endif

if(inout == 'OUT') then
   
   ndims=1 ; idims(1)=1
   CALL shdf5_orec (ndims,idims,'sigz_nsigz',ivars=nsigz)
   ndims=1 ; idims(1)=nsigz
   CALL shdf5_orec (ndims,idims,'sigz_sigz',rvara=sigz)

   npts=n1*n2*nsigz
   ndims=3 ; idims(1)=n1 ; idims(2)=n2 ; idims(3)=nsigz
   CALL vmissw (ps_u,npts,ps_scra,1E30,-9999.)
   CALL shdf5_orec (ndims,idims,'sigz_u',rvara=ps_scra)
   CALL vmissw (ps_v,npts,ps_scra,1E30,-9999.)
   CALL shdf5_orec (ndims,idims,'sigz_v',rvara=ps_scra)
   CALL vmissw (ps_p,npts,ps_scra,1E30,-9999.)
   CALL shdf5_orec (ndims,idims,'sigz_p',rvara=ps_scra)
   CALL vmissw (ps_t,npts,ps_scra,1E30,-9999.)
   CALL shdf5_orec (ndims,idims,'sigz_t',rvara=ps_scra)
   CALL vmissw (ps_r,npts,ps_scra,1E30,-9999.)
   CALL shdf5_orec (ndims,idims,'sigz_r',rvara=ps_scra)

   print 201,' *****  Sigma-z file written *************'  &
        ,iyear,imonth,idate,ihour,n1,n2,nsigz   &
        ,(sigz(l),l=1,nsigz)

endif

return
END SUBROUTINE sigzio

!##############################################################################
Subroutine vmissw (af,n,as,fm,fx)

implicit none

integer :: n
real :: af(*),as(*),fm,fx
integer :: i

do i=1,n
   as(i)=af(i)
   if(af(i).ge.fm) as(i)=fx
enddo

return
END SUBROUTINE vmissw

!##############################################################################
Subroutine vmissr (af,n,fm,fx)

implicit none

integer :: n
real :: af(*),fm,fx
integer :: i

do i=1,n
   if(af(i).le.fx) af(i)=fm
enddo

return
END SUBROUTINE vmissr

