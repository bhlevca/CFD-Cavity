!Sample program for solving Lid-Driven cavity test using SIMPLE-algorithm
! Main modul
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

!-----------------------------------------------------------------------
Program Main

include 'icomm_1.f90'

Dimension F_out(nx,ny),x_exct(30), v_exct(30)
Dimension              y_exct(30), u_exct(30)
Character  Filename*10

	Call Ini
	Call Grid_rectangular

	Call Geom
	Call Init_all_cavity

!

! ------------------------------------------------
! 	Open(10,file='binary_dat_all.dat')

! 	 	 read(10,*)F(:,:,1:5)

! 	close(10)

! ------------------------------------------------
Niter_last = 0
Niter_print = 0

	Do 100 Niter=1,1500


write(*,*) '------------------------',Niter,'------------------------',Niter


 	Call Solve_UV
    Call Solve_Pressure_Correction


!--------------------------------------------------------------------------
	if((Niter-Niter_last).eq.10)then

		Niter_last = Niter

	open (23,file='Domain_all.dat')

		WRITE(23,*)'VARIABLES = "Xp", "Yp" , "Up" , "Vp" , "Pp" '

        WRITE (23,*)' ZONE I=' ,NXmaxC, ', J=', NYmaxC, ', F=POINT'

	DO 88 J=1, NYmaxC
	DO 88 I=1, NXmaxC

	 WRITE (23,*) Xc(I,J), Yc(I,J) , F(i,j,1) , F(i,j,2) , F(i,j,4)

  88     continue

	close(23)
!------------------------------------------------------
	open (23,file='Profiles_V.dat')


	    WRITE(23,*)'VARIABLES = "Xcalc", "Vcalc" '

        WRITE (23,*)' ZONE F=POINT, I=', NYmaxC


	DO 881 i=1,NYmaxC

	 WRITE (23,*) Xc(i,61), F(i,61,2)

  881     continue


x_exct(1)=0.00282901 ; v_exct(1)=-0.07028612
x_exct(2)=6.42882E-2 ; v_exct(2)=2.71771E+0
x_exct(3)=7.35542E-2 ; v_exct(3)=2.86215E+0
x_exct(4)=8.09753E-2 ; v_exct(4)=3.02728E+0
x_exct(5)=9.76486E-2 ; v_exct(5)=3.25421E+0
x_exct(6)=1.58722E-1 ; v_exct(6)=3.70750E+0
x_exct(7)=2.28897E-1 ; v_exct(7)=3.31348E+0
x_exct(8)=2.34432E-1 ; v_exct(8)=3.25139E+0
x_exct(9)=5.01948E-1 ; v_exct(9)=1.87997E-1
x_exct(10)=8.04508E-1; v_exct(10)=-3.33066E+0
x_exct(11)=8.59780E-1; v_exct(11)=-4.42685E+0
x_exct(12)=9.07688E-1; v_exct(12)=-5.33693E+0
x_exct(13)=9.44865E-1; v_exct(13)=-4.07737E+0
x_exct(14)=9.54199E-1; v_exct(14)=-3.51971E+0
x_exct(15)=9.61692E-1; v_exct(15)=-2.92069E+0
x_exct(16)=9.69195E-1; v_exct(16)=-2.25968E+0
x_exct(17)=1.00098E+0; v_exct(17)=-9.09091E-2



        WRITE (23,*)' ZONE F=POINT, I=', 17


	DO 882 i=1,17

	 WRITE (23,*) x_exct(i), v_exct(i)/10.

  882     continue


	close(23)
!--------------------------------------------------------
	open (24,file='Profiles_U.dat')

   WRITE(24,*)'VARIABLES = "Ycalc", "Ucalc" '

        WRITE (24,*)' ZONE F=POINT, I=', NYmaxC


	DO 871 j=1,NYmaxC

	 WRITE (24,*) Yc(61,j), F(61,j,1)

  871     continue


y_exct(1)=-1.85185E-3 ; u_exct(1)=0.00000E+0
y_exct(2)=5.00000E-2 ; u_exct(2)=-1.84615E+0
y_exct(3)=5.92593E-2 ; u_exct(3)=-2.03077E+0
y_exct(4)=6.66667E-2 ; u_exct(4)=-2.21538E+0
y_exct(5)=1.00000E-1 ; u_exct(5)=-2.98462E+0
y_exct(6)=1.68519E-1 ; u_exct(6)=-3.81538E+0
y_exct(7)=2.77778E-1 ; u_exct(7)=-2.76923E+0
y_exct(8)=4.48148E-1 ; u_exct(8)=-1.07692E+0
y_exct(9)=4.96296E-1 ; u_exct(9)=-6.15385E-1
y_exct(10)=6.09259E-1 ; u_exct(10)=5.84615E-1
y_exct(11)=7.31481E-1 ; u_exct(11)=1.84615E+0
y_exct(12)=8.50000E-1 ; u_exct(12)=3.32308E+0
y_exct(13)=9.50000E-1 ; u_exct(13)=4.64615E+0
y_exct(14)=9.57407E-1 ; u_exct(14)=5.07692E+0
y_exct(15)=9.64815E-1 ; u_exct(15)=5.72308E+0
y_exct(16)=9.74074E-1 ; u_exct(16)=6.58462E+0
y_exct(17)=9.96296E-1 ; u_exct(17)=1.00000E+1



        WRITE (24,*)' ZONE F=POINT, I=', 17

	DO 872 j=1,17

	 WRITE (24,*) y_exct(j), u_exct(j)/10.

  872     continue


	close(24)
	end if
!--------------------------------------------------------------------------
if((Niter-Niter_print).eq.50)then

    Niter_print = Niter

	Open(10,file='binary_dat_all.dat')

	  write(10,*)F(:,:,1:5)
	close(10)

end if
!--------------------------------------------------------------------------

	100 continue

!********************************************************************************
!********************************************************************************
!	Open(10,file='binary_dat_all.dat')
!
!	  write(10,*)F(:,:,1:5)
!	close(10)

!----------------------------------------------------------------
	NImax = NXmaxC
	NJmax = NYmaxC

	F_out     = F(:,:,1)

!		        1234567890
	Filename  ='1_U_s.txt'

    Call  Out_array(F_out,NImax,NJmax,Filename)
!-------------------------------------------------------------------
!----------------------------------------------------------------
	NImax = NXmaxC
	NJmax = NYmaxC

	F_out     = F(:,:,2)

!		        1234567890
	Filename  ='1_V_s.txt'

    Call  Out_array(F_out,NImax,NJmax,Filename)
!-------------------------------------------------------------------
!----------------------------------------------------------------
	NImax = NXmaxC
	NJmax = NYmaxC

	F_out     = F(:,:,5)

!		        1234567890
	Filename  ='1_T_s.txt'

!    Call  Out_array(F_out,NImax,NJmax,Filename)
!-------------------------------------------------------------------

!--------------------------------------------------------------------------
	open (23,file='Domain_all.dat')

	WRITE(23,*)'VARIABLES = "Xp", "Yp" , "Up" , "Vp" , "Pp" '

        WRITE (23,*)' ZONE I=' ,NXmaxC, ', J=', NYmaxC, ', F=POINT'

	DO 4 J=1, NYmaxC
	DO 4 I=1, NXmaxC

	 WRITE (23,*) Xc(I,J), Yc(I,J) , F(i,j,1) , F(i,j,2) , F(i,j,4)

  4     continue

	close(23)
!--------------------------------------------------------------------------

	WRITE(*,*) 'PRIVET'
	STOP
	END
