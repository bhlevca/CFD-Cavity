!Sample program for solving Lid-driven cavity flow test using SIMPLE-algorithm
! solution of pressure-correction quation and correction U,V, and P
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
Subroutine Solve_Pressure_Correction
include 'icomm_1.f90'

      DpU(2:NXmax,2:NYmax) = 1./Ap(2:NXmax,2:NYmax,1)
      DpV(2:NXmax,2:NYmax) = 1./Ap(2:NXmax,2:NYmax,2)


!****************************************************************************************
! vertical  faces East - e
   Do 101 I=1,NXmax
   Do 101 J=1,NYmax-1

    If((i.ne.1).and.(i.ne.NXmax))then


 	     DXc_e = Xc(i+1,j+1) - Xc(i,j+1)
	     S_e   =  Y(i  ,j+1) -  Y(i,j  )

		 VOL_e = DXc_e * S_e

	!-------------------------------------------------------------
         APU_e =       0.5 * (   DpU(i,j+1)   +   DpU(i+1,j+1)   )
         Ul_e  =       0.5 * (     F(i,j+1,1) +     F(i+1,j+1,1) )
         DPl_e =       0.5 * ( DPx_c(i,j+1)   + DPx_c(i+1,j+1)   )

         DPx_e =             ( F(i+1,j+1,4) - F(i,j+1,4) ) / DXc_e

         U_e   = Ul_e  + APU_e * VOL_e * ( DPl_e - DPx_e)

        !--------------------------------------------------------------
		 Con_e(i,j) =   U_e  * S_e


	End If


 	     Con_e(1,j) = 0.
  	     Con_e(NXmax,j) = 0.

 101 continue

!****************************************************************************************
! horisontal  faces
   Do 102 I=1,NXmax-1
   Do 102 J=1,NYmax

    If((j.ne.1).and.(j.ne.NYmax))then

 	     DYc_n = Yc(i+1,j+1) - Yc(i+1,j)
	     S_n   =  X(i+1,j  ) -  X(i  ,j)

	     VOL_n = DYc_n * S_n

	!-----------------------------------------------------------
         APV_n =       0.5 * (   DpV(i+1,j)   +   DpV(i+1,j+1)   )
         Vl_n  =       0.5 * (     F(i+1,j,2) +     F(i+1,j+1,2) )
         DPl_n =       0.5 * ( DPy_c(i+1,j)   + DPy_c(i+1,j+1)   )


         DPy_n = ( F(i+1,j+1,4) - F(i+1,j,4) ) / DYc_n

  		 V_n   = Vl_n  + APV_n * VOL_n * ( DPl_n - DPy_n )
        !-----------------------------------------------------------

		 Con_n(i,j) =   V_n  * S_n


	End If


	     Con_n(i,1) = 0.
 	     Con_n(i,NYmax) = 0.

 102 continue

!****************************************************************************************
  Summ = 0.

	Do 200 I=2,NXmax
	Do 200 J=2,NYmax

		 S_e   =  Y(i  ,j) -  Y(i,j-1)
		 S_n   =  X(i  ,j) -    X(i-1,j)


 	 An(i,j) = 0.5 * (   DpV(i,j)   +   DpV(i,j+1)   )  * S_n * S_n
	 As(i,j) = 0.5 * (   DpV(i,j)   +   DpV(i,j-1)   )  * S_n * S_n
	 Ae(i,j) = 0.5 * (   DpU(i,j)   +   DpU(i+1,j)   )  * S_e * S_e
	 Aw(i,j) = 0.5 * (   DpU(i,j)   +   DpU(i-1,j)   )  * S_e * S_e

	 Ap(i,j,3) = (An(i,j) + As(i,j) + Ae(i,j) +	Aw(i,j))

    Sp(i,j,3) = -1. * ( Con_n(i-1,j) - Con_n(i-1,j-1) + Con_e(i,j-1) - Con_e(i-1,j-1) )

	 Summ = Summ + Sp(i,j,3) !
 200 continue

			write(*,*)'Sum',summ
!***********************************************************************************
!***********************************************************************************
		write(*,*) Summ

write(*,*)'solve PP'
 niter = 0


	20	call Convergence_Criteria(3,Res_sum_before)
         niter= niter + 1

		Call TDMA_1(3)
		call Convergence_Criteria(3,Res_sum_After)
	If((abs(Res_sum_before-Res_sum_After).Ge.0.0000001).and.(niter.le.500).and.(abs(Res_sum_After).ge.0.0001))go to 20
	write(*,*)'Res_sum_before-Res_sum_After',Res_sum_before-Res_sum_After,niter




!***********************************************************************************
!***********************************************************************************
!  velocities and pressure correction

	Do 300 I=2,NXmax
	Do 300 J=2,NYmax

		DY = Y(i,j)-Y(i,j-1)
		DX = X(i,j)-X(i-1,j)

		PPe = 0.5 * ( F(i,j,3) + F(i+1,j,3) )
		PPw = 0.5 * ( F(i,j,3) + F(i-1,j,3) )
		PPn = 0.5 * ( F(i,j,3) + F(i,j+1,3) )
		PPs = 0.5 * ( F(i,j,3) + F(i,j-1,3) )

		F(i,j,1) = F(i,j,1) + (PPw - PPe)  * DpU(i,j)  * Dy
		F(i,j,2) = F(i,j,2) + (PPs - PPn)  * DpV(i,j)  * Dx

		F(i,j,4) = F(i,j,4) + 0.1 *   F(i,j,3)

 300 continue

!***********************************************************************************
	! correct convective fluxes  through vertical East faces
   Do 701 I=1,NXmax
   Do 701 J=1,NYmax-1

    If((i.ne.1).and.(i.ne.NXmax))then

	  Con_e(i,j) =  Con_e(i,j) + Ae(i,j+1) * (F(i+1,j+1,3)-F(i,j+1,3) )

     end if

	 if((i.eq.1).or.(i.eq.NXmax))Con_e(i,j)=0.

	701 continue

! correct  convective fluxes through horisontal  North faces
   Do 702 I=1,NXmax-1
   Do 702 J=1,NYmax

    If((j.ne.1).and.(j.ne.NYmax))then

		Con_n(i,j) = Con_n(i,j) + An(i+1,j) * (F(i+1,j+1,3)-F(i+1,j,3) )

     end if

	if((j.eq.1).or.(j.eq.NYmax))Con_n(i,j)=0.

	702 continue

!***********************************************************************************
	Do 501 I=1,NXmaxC

       F(i,1,4) = F(i,2,4)            !
       F(i,NYmaxC,4) = F(i,NYmaxC-1,4)!


 501 continue

 	Do 502 J=1,NYmaxC

       F(1,j,4) = F(2,j,4)             !
       F(NXmaxC,j,4) = F(NXmaxC-1,j,4) !


 502 continue

 F(:,:,4) = F(:,:,4) - F(3,4,4)
!***********************************************************************************



Return
End
