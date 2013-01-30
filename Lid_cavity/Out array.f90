!Sample program for solving Lid-driven cavity flow test using SIMPLE-algorithm
! array output into the .txt file modul
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
Subroutine Out_array(F_out,NImax,NJmax,Filename)

Dimension F_out(500,500)

Character  Filename*10

	Open(10,file=Filename)

		Do 12 J=1,NJmax

!write(10,*)'ssss'

	 		write(10,11)j,(F_out(i,j),i=1,NImax)

		11 format(I5,15F6.3)

		12 continue


	close(10)

Return
End
