!Sample program for solving Lid-driven cavity flow test using SIMPLE-algorithm
! covective terms approximation -  HLPA scheme implementation modul
!Copyright (C) 2010  Michail Kiričkov

!This program is free software; you can redistribute it and/or
!modify it under the terms of the GNU General Public License
!as published by the Free Software Foundation; either version 2
!of the License, or (at your option) any later version.

!This program is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU General Public License for more details.

!You should have received a copy of the GNU General Public License
!along with this program; if not, write to the Free Software
!Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

!**********************************************************************
! HLPA scheme

Subroutine HLPA(Uw,Fww,Fw,Fp,Fe,Delta_f)


If (Uw.GE.0.) then

	if ( ABS( Fp - 2.* Fw + Fww ).LT.ABS( Fp - Fww ) ) then

		Alpha_pl = 1.
												 else

		Alpha_pl = 0.

	End If

	If ( (Fp - Fww).NE.0. ) then

		Delta_f = Alpha_pl * (Fp - Fw)* (Fw - Fww) / (Fp - Fww)

	End If

End If

!-------------------------------------------------------------------------

If (Uw.LT.0.) then

	if ( ABS( Fw - 2.* Fp + Fe ).LT.ABS( Fw - Fe ) ) then

		Alpha_mn = 1.
												 else

		Alpha_mn = 0.

	End If

	If ( (Fw - Fe).NE.0. ) then

		Delta_f = Alpha_mn * (Fw - Fp)* (Fp - Fe) / (Fw - Fe)

	End If

End If

100 continue

Return

End
