!Sample program for solving Lid-driven cavity flow test using SIMPLE-algorithm
! calculation of grid modul
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
Subroutine Grid_rectangular

include 'icomm_1.f90'

x=0.
y=0.

 SLx = 1.
 SLy = 1.

 Xbeg =  0.
 Ybeg =  0.


NXmax =  121
NYmax =  121

NXmaxC =  NXmax + 1
NYmaxC =  NYmax + 1

!-------------------------------------------------------------------------------------------
	DO 2 I=1, NXmax
	DO 2 J=1, NYmax

	X(I,J) = Xbeg + (SLx /(NXmax-1))* (i-1)

  2     continue

 	DO 3 J=1, NYmax
		DO 3 I=1, NXmax

	Y(I,J) = Ybeg + (SLy /(NYmax-1))* (j-1)

  3     continue

!----------------------------------------------------------

	open (22,file='GRID.dat')

	WRITE(22,*)'VARIABLES = "X", "Y" '

      WRITE (22,*)' ZONE I=' ,NXmax, ', J=', NYmax, ', F=POINT'

	DO 1 J=1, NYmax
	DO 1 I=1, NXmax

	 WRITE (22,*) X(I,J), Y(I,J)

  1     continue


	open (21,file='GRID.txt')

	write(21,*) NYmax,NXmax

	DO 4 J=1, NYmax
	DO 4 I=1, NXmax

	 WRITE (21,*) X(I,J), Y(I,J)

  4     continue

	close(21)

Return
End
