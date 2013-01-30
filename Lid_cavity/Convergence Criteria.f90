!Sample program for solving Lid-driven cavity flow test using SIMPLE-algorithm
! calculation of covergence criteria modul
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
Subroutine Convergence_Criteria(NF,res_sum)

include 'icomm_1.f90'
Dimension Res(nx,ny)

res =  0.


Res_Sum = 0.
Res_vol = 0.

	Do 20 I=2,NXmaxC-1
	Do 20 J=2,NYmaxC-1

	  Res_vol =	 Ap(i,j,nf) * F(i  ,j  ,nf) -   &

		 (	    As(i,j) * F(i  ,j-1,nf) +   &
			    An(i,j) * F(i  ,j+1,nf) +   &
			    Aw(i,j) * F(i-1,j  ,nf) +   &
			    Ae(i,j) * F(i+1,j  ,nf)  )  &
													- Sp(i,j,nf)

	  res(i,j) =  Res_vol

	  Res_Sum = Res_Sum + Res_vol

	20 continue

write(*,*)'covergence',Res_Sum

Return
End
