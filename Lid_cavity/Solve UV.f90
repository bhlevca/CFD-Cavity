!Sample program for solving Lid-driven cavity flow test using SIMPLE-algorithm
! solution of momentum equation for U and V modul
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
Subroutine Solve_UV
include 'icomm_1.f90'


!  calculation of fluxes
!  all geometry has rectangular 2D notation

	Do 100 I= 2,NXmax
	Do 100 J= 2,NYmax

		Gam_e = ( Gam(i+1,j  ) + Gam(i  ,j  ) ) * 0.5
		Gam_w = ( Gam(i-1,j  ) + Gam(i  ,j  ) ) * 0.5
		Gam_s = ( Gam(i  ,j-1) + Gam(i  ,j  ) ) * 0.5
		Gam_n = ( Gam(i  ,j+1) + Gam(i  ,j  ) ) * 0.5


	   !----------------------------------------------------
		Area_w = Y(i-1,j)-Y(i-1,j-1)
		Area_e = Y(i  ,j)-Y(i  ,j-1)

		Area_s = X(i,j-1)-X(i-1,j-1)
		Area_n = X(i,j  )-X(i-1,j  )
	   !----------------------------------------------------

		Del_w  = Xc(i  ,j)-Xc(i-1,j)
		Del_e  = Xc(i+1,j)-Xc(i  ,j)

		Del_s  = Yc(i,j  )-Yc(i,j-1)
		Del_n  = Yc(i,j+1)-Yc(i,j  )
	   !----------------------------------------------------

! upwind differencing (all other will be included into the source term)

 	    Conv_w = Area_w *  ( F(i,j,1) + F(i-1,j  ,1) ) * 0.5
            Conv_e = Area_e *  ( F(i,j,1) + F(i+1,j  ,1) ) * 0.5
 	    Conv_s = Area_s *  ( F(i,j,2) + F(i  ,j-1,2) ) * 0.5
 	    Conv_n = Area_n *  ( F(i,j,2) + F(i  ,j+1,2) ) * 0.5

		if(i.eq.2    )Conv_w = 0.
		if(i.eq.NXmax)Conv_e = 0.
                if(j.eq.2    )Conv_s = 0.
		if(j.eq.NYmax)Conv_n = 0.

		Diff_e = Area_e * Gam_e / Del_e
 		Diff_w = Area_w * Gam_w / Del_w
		Diff_s = Area_s * Gam_s / Del_s
		Diff_n = Area_n * Gam_n / Del_n


		Aw(i,j) = Diff_w + max(     Conv_w,0.)
		Ae(i,j) = Diff_e + max(-1.* Conv_e,0.)
		As(i,j) = Diff_s + max(     Conv_s,0.)
		An(i,j) = Diff_n + max(-1.* Conv_n,0.)

	Ap(i,j,1:2)= Aw(i,j) + Ae(i,j) + An(i,j) + As(i,j)

		Sp(i,j,1:2)= 0.

!-------------------------------- HLPA SCHEME----------------------------
 !  go to 600 ! (now HLPA is "off")

	DO 500 nf=1,2

! Subroutine HLPA(Uw,Fww,Fw,Fp,Fe,Delta_f)

 if( (i.GT.2).AND.(i.LT.NXmax-0).and.(j.GT.2).AND.(j.LT.NYmax-0) ) then

   !------------------ w face -------------------
	Fww = F(i-2,j,nf)
	Fw  = F(i-1,j,nf)
	Fp  = F(i  ,j,nf)
	Fe  = F(i+1,j,nf)

	call  HLPA(Conv_w,Fww,Fw,Fp,Fe,Delta_f)

	Sp(i,j,nf) = Sp(i,j,nf) + Conv_w * Delta_f

  !------------------ e face--------------------

	Fww = F(i-1,j,nf)
	Fw  = F(i  ,j,nf)
	Fp  = F(i+1,j,nf)
	Fe  = F(i+2,j,nf)

	call  HLPA(Conv_e,Fww,Fw,Fp,Fe,Delta_f)

    Sp(i,j,nf) = Sp(i,j,nf) + Conv_e * Delta_f * (-1.)

  !------------------ s face--------------------
	Fww = F(i  ,j-2,nf)
	Fw  = F(i  ,j-1,nf)
	Fp  = F(i  ,j  ,nf)
	Fe  = F(i  ,j+1,nf)

	call  HLPA(Conv_s,Fww,Fw,Fp,Fe,Delta_f)

   Sp(i,j,nf) = Sp(i,j,nf) + Conv_s * Delta_f

  !------------------ n face--------------------

	Fww = F(i  ,j-1,nf)
	Fw  = F(i  ,j  ,nf)
	Fp  = F(i  ,j+1,nf)
	Fe  = F(i  ,j+2,nf)

	call  HLPA(Conv_n,Fww,Fw,Fp,Fe,Delta_f)

   Sp(i,j,nf) = Sp(i,j,nf) + Conv_n * Delta_f *(-1.)


 end if

   500 continue


   600 continue


!------------------------------------------------------------------------

100 continue ! coefficient cycle

!----------------------------- pressure gradient ------------------------

  	Do 200 I= 2,NXmax
	Do 200 J= 2,NYmax

      DX = X(i,j) - X(i-1,j)
      DY = Y(i,j) - Y(i,j-1)

	  VOL = DX * DY

	  PE = ( F(i,j,4) + F(i+1,j,4) ) * 0.5
	  PW = ( F(i,j,4) + F(i-1,j,4) ) * 0.5
	  PN = ( F(i,j,4) + F(i,j+1,4) ) * 0.5
	  PS = ( F(i,j,4) + F(i,j-1,4) ) * 0.5

	  DPx_c(i,j) = (PE-PW)/DX
	  DPy_c(i,j) = (PN-PS)/DY

	  Sp(i,j,1) = Sp(i,j,1) - DPx_c(i,j) * VOL
	  Sp(i,j,2) = Sp(i,j,2) - DPy_c(i,j) * VOL

	200 continue

!---------------------------- under-relaxation ---------------------------------

Alfa = 0.85
Urf  = 1. / Alfa


    Ap(1:NXmaxC,1:NYmaxC,1) = Ap(1:NXmaxC,1:NYmaxC,1)  * Urf
    Sp(1:NXmaxC,1:NYmaxC,1) = Sp(1:NXmaxC,1:NYmaxC,1)  + (1. - Alfa )* Ap(1:NXmaxC,1:NYmaxC,1)*F(1:NXmaxC,1:NYmaxC,1) ! / Alfa

    Ap(1:NXmaxC,1:NYmaxC,2) = Ap(1:NXmaxC,1:NYmaxC,2)  * Urf
    Sp(1:NXmaxC,1:NYmaxC,2) = Sp(1:NXmaxC,1:NYmaxC,2)  + (1. - Alfa )* Ap(1:NXmaxC,1:NYmaxC,2)*F(1:NXmaxC,1:NYmaxC,2) ! / Alfa

!---------------------------------------------------------------------------------

!******************************************************************

niter = 0
write(*,*)'solve U'
	call Convergence_Criteria(1,Res_sum_before)

	10 continue
	     niter= niter + 1
		Call TDMA_1(1)
		call Convergence_Criteria(1,Res_sum_After)
	If((abs(Res_sum_before-Res_sum_After).Ge.0.00000000001).and.niter.le.20)then

		         Res_sum_before = Res_sum_After
				     go to 10
    End if



niter = 0
write(*,*)'solve V'
	call Convergence_Criteria(2,Res_sum_before)

	20	continue

     niter= niter + 1
		Call TDMA_1(2)
		call Convergence_Criteria(2,Res_sum_After)
	If((abs(Res_sum_before-Res_sum_After).Ge.0.00000000001).and.niter.le.20)then

		        Res_sum_before = Res_sum_After
				go to 20
    End if

!***********************************************************************

!---------------------------------------------------------------------------------
Return
End
