
!Sample program for solving Lid-driven cavity flow test using SIMPLE-algorithm
! include modul
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
parameter nx=500, ny=500

 DOUBLE PRECISION  U(nx,ny),V(nx,ny),F(nx,ny,10), &
			Xc(nx,ny),Yc(nx,ny),X(nx,ny),Y(nx,ny),Gam(nx,ny),Ro(nx,ny)

 DOUBLE PRECISION  X_xi(nx,ny),    &
				  Y_xi(nx,ny),    &
				  X_et(nx,ny),    &
				  Y_et(nx,ny),	  &

                          Del_X_xi(nx,ny),    &
			  Del_Y_xi(nx,ny),    &

			  Del_X_et(nx,ny),    &
			  Del_Y_et(nx,ny),    &


				  Dx_c(nx,ny),    &
				  Dy_c(nx,ny)

 DOUBLE PRECISION  DPx_c(nx,ny),  &
	           DPy_c(nx,ny),  &

		   Con_e(nx,ny),  &
		   Con_n(nx,ny),  &
		   CheckFlux(nx,ny),  &

				 DpU(nx,ny),   &
				 DpV(nx,ny),   &

			      Ap(nx,ny,10),   &
			      As(nx,ny),   &
			      An(nx,ny),   &
			      Aw(nx,ny),   &
			      Ae(nx,ny),   &
			      Sp(nx,ny,10)



!-----------------------------------------------


Common /var/  U,V,F, &
			Xc,Yc,X,Y,Gam,Ro

Common /var_Geom/  X_xi,    &
		   Y_xi,    &
		   X_et,    &
		   Y_et,    &

              Del_X_xi,    &
 	      Del_Y_xi,    &

	      Del_X_et,    &
	      Del_Y_et,    &

		  Dx_c,    &
		  Dy_c

Common /koeff/ Con_e,  &
 	       Con_n,  &

	       DPx_c,  &
	       DPy_c,  &

		   CheckFlux,  &

				 DpU,   &
				 DpV,   &

			      Ap,   &
			      As,   &
			      An,   &
			      Aw,   &
			      Ae,   &
			      Sp


Common /geomt/ NXmax,NYmax,NXmaxC,NYmaxC  