!Sample program for solving Lid-driven cavity flow test using SIMPLE-algorithm
! solution of linear system of equations by Thomas algorithm modul
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
Subroutine TDMA_1(NF)

include 'icomm_1.f90'

DOUBLE PRECISION  P(nx),Q(nx)

!--------------------------------------------------------------------------------------------------------------------------------------------------------
	Do 101 J =  2, NYmax

	 	    P(1) =  0.
  		    Q(1) =  F(1,j,nf)

		    P(NXmaxC) = 0.
		    Q(NXmaxC) = F(NXmaxC,j,nf)

 !		Forward Elimination

	      Do 10 i = 2,NXmaxC-1

			temp =  Ap(i,j,nf) - Aw(i,j) * P(i-1)

 			Spp= Sp(i,j,nf) + As(i,j) * F(i,j-1,nf) + &
		 	  	          An(i,j) * F(i,j+1,nf)

					P(i) = Ae(i,j) / temp

					Q(i) = (Spp + Aw(i,j)*Q(i-1)) / temp

		10    continue

!		 Back Substitution

	  Do 20 i = NXmaxC-1,1,-1

          F(i,j,nf) = P(i)*F(i+1,j,nf) + Q(i)

		20    continue

     101 continue
!--------------------------------------------------------------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------------------------------------------------------------
	Do 301 J =  NYmax,2,-1

	 	    P(1) =  0.
 		    Q(1) =  F(1,j,nf)

		    P(NXmaxC) = 0.
		    Q(NXmaxC) = F(NXmaxC,j,nf)

 !		Forward Elimination

	      Do 32 i = 2,NXmaxC-1

		    temp =  Ap(i,j,nf) - Aw(i,j) * P(i-1)

 		    Spp= Sp(i,j,nf) + As(i,j) * F(i,j-1,nf) + &
				      An(i,j) * F(i,j+1,nf)

					P(i) = Ae(i,j) / temp

					Q(i) = (Spp + Aw(i,j)*Q(i-1)) / temp

		32    continue

!		 Back Substitution

	  Do 30 i = NXmaxC-1,2,-1

         F(i,j,nf) = P(i)*F(i+1,j,nf) + Q(i)

		30    continue

     301 continue
!--------------------------------------------------------------------------------------------------------------------------------------------------------


Return
End