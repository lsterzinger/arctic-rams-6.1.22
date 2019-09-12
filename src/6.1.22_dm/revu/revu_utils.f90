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

Subroutine error_mess (msg)

!Used only for REVU

implicit none

character(len=*) :: msg

print*,msg

return
END SUBROUTINE error_mess

!##############################################################################
integer Function lastchar (str)

!Used only for REVU

implicit none

character(len=*) :: str
integer :: n,ln

! returns last non-blank character position from a string

ln=len(str)
do n=ln,1,-1
   if(str(n:n).ne.' ') then
      lastchar=n
      return
   endif
enddo
lastchar=0

return
END FUNCTION lastchar

!##############################################################################
Subroutine tokenize (str1,tokens,ntok,toksep,nsep)

!Used only for REVU

use grid_dims

implicit none

integer :: nsep,ntok,npt,nch,nc,ns
character(len=*) :: str1,tokens(*)
character(len=1) :: toksep(nsep)
character(len=strl1) :: str

! this routine "parses" character string str into different pieces
! or tokens by looking for possible token separators (toks)
! str contains nch characters. The number of tokens identified is nto
! the character string tokens are stored in tokens.

ntok=0
npt=1
CALL deblank (str1,str,nch)
do nc=1,nch
   do ns=1,nsep
      if(str(nc:nc).eq.toksep(ns).or.nc.eq.nch) then
         if(nc-npt.ge.1)then
            ntok=ntok+1
            tokens(ntok)=str(npt:nc-1)
            if(nc.eq.nch.and.str(nc:nc).ne.toksep(ns)) then
               tokens(ntok)=str(npt:nc)
               goto 10
            endif
         endif
         ntok=ntok+1
         tokens(ntok)=str(nc:nc)
         npt=nc+1
         goto 10
      endif
   enddo
10      continue
enddo

return
END SUBROUTINE tokenize
