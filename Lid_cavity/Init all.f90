!Sample program for solving Lid-driven cavity flow test using SIMPLE-algorithm
! initiation of arrays modul
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

Subroutine Ini

include 'icomm_1.f90'


  U(:,:) = 0.
  V(:,:) = 0.

  F(:,:,:) = 0.
 Xc(:,:) = 0.
 Yc(:,:) = 0.
  X(:,:) = 0.
  Y(:,:) = 0.
Gam(:,:) = 0.
 Ro(:,:) = 1.

Con_e(:,:) = 0.
Con_n(:,:) = 0.

DPx_c(:,:) = 0.
DPy_c(:,:) = 0.

 Dx_c(:,:) = 0.
 Dy_c(:,:) = 0.

 Ap(:,:,:) = 0.
 As(:,:) = 0.
 An(:,:) = 0.
 Aw(:,:) = 0.
 Ae(:,:) = 0.
 Sp(:,:,:) = 0.


 Return
 End
