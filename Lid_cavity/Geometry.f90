
!Sample program for solving Lid-driven cavity flow test using SIMPLE-algorithm
! Calculation of Xc and Yc with possibility for further development modul
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
Subroutine Geom

include 'icomm_1.f90'

! calculation Xc,Yc

! ------------------------------------------------------------------------

	do  2 I=2,NXmax
        do  2 J=2,NYmax

     	    Xc(I,J)=(  X(i-1,j-1) + X(i-1,j  ) + &
     	               X(i  ,j  ) + X(i  ,j-1)      ) * 0.25

	    Yc(I,J)=(  Y(i-1,j-1) + Y(i-1,j  ) + &
         	       Y(i  ,j  ) + Y(i  ,j-1)     ) * 0.25

	2 continue

! ------------------------------------------------------------------------
	do 4 I=2,NXmax

		Xc(i,1      ) = ( X(i  ,1    ) + X(i-1,1    ) ) * 0.5
		Xc(i,NYmax+1) = ( X(i  ,NYmax) + X(i-1,NYmax) ) * 0.5


		Yc(i,1      ) = ( Y(i  ,1    ) + Y(i-1,1    ) ) * 0.5
		Yc(i,NYmax+1) = ( Y(i  ,NYmax) + Y(i-1,NYmax) ) * 0.5

	4 continue

 ! ------------------------------------------------------------------------

	Xc(1      ,      1) = X(    1,    1)
	Xc(NXmax+1,      1) = X(NXmax,    1)
	Xc(      1,NYmax+1) = X(    1,NYmax)
	Xc(NXmax+1,NYmax+1) = X(NXmax,NYmax)

	Yc(1      ,      1) = Y(    1,    1)
	Yc(NXmax+1,      1) = Y(NXmax,    1)
	Yc(      1,NYmax+1) = Y(    1,NYmax)
	Yc(NXmax+1,NYmax+1) = Y(NXmax,NYmax)

!--------------------------------------------------------------------------
! ------------------------------------------------------------------------
	do 5 J=2,NYmax

		Yc(1      ,j ) = ( Y(1     ,j) + Y(1    ,j-1) ) * 0.5
		Yc(NXmax+1,j ) = ( Y(NXmax ,j) + Y(NXmax,j-1) ) * 0.5
		Xc(1      ,j ) = ( X(1     ,j) + X(1    ,j-1) ) * 0.5
		Xc(NXmax+1,j ) = ( X(NXmax ,j) + X(NXmax,j-1) ) * 0.5

	5 continue
!--------------------------------------------------------------------------
! ------------------------------------------------------------------------
! Xi (vertical)

	Do 101 I=1,NXmax
	Do 101 J=1,NYmax-1

		X_xi(I,J) = X(i  ,j+1) - X(i  ,j  )
		Y_xi(I,J) = Y(i  ,j+1) - Y(i  ,j  )

	101 continue

! Eta (horisontal)

	Do 102 I=1,NXmax-1
	Do 102 J=1,NYmax

		X_et(I,J) = X(i+1,j  ) - X(i  ,j  )
		Y_et(I,J) = Y(i+1,j  ) - Y(i  ,j  )

	102 continue

!--------------------------------------------------------------------------
! ------------------------------------------------------------------------

! Xi (vertical)

	Do 201 I=1,NXmaxC
	Do 201 J=1,NYmax

		Del_X_xi(i  ,j  ) =  Xc(i  ,j+1) - Xc(i  ,j  )
		Del_Y_xi(i  ,j  ) =  Yc(i  ,j+1) - Yc(i  ,j  )

	201 continue


! Eta (horisontal)

	Do 202 I=1,NXmax
	Do 202 J=1,NYmaxC

		Del_X_et(i  ,j  ) =  Xc(i+1,j  ) - Xc(i  ,j  )
		Del_Y_et(i  ,j  ) =  Yc(i+1,j  ) - Yc(i  ,j  )

	202 continue

!--------------------------------------------------------------------------

	Do 303 I=2,NXmaxC-1
	Do 303 J=2,NYmaxC-1

	     Dx_c(i,j) = X(i,j) - X(i-1,j)
 	     Dy_c(i,j) = Y(i,j) - Y(i,j-1)


	303 continue

Return
End
