!Sample program for solving Lid-driven cavity flow test using SIMPLE-algorithm
! initiation data for lid-driven cavity flow Re=1000 modul
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

Subroutine Init_all_cavity

include 'icomm_1.f90'

!------------------------------------------------

     F(:,:,1) = 0.0
     F(:,:,2) = 0.0
     F(:,:,3) = 0.0
     F(:,:,4) = 0.0


     Gam(:,:) = 1.  /1000.
     Ro(:,:)  = 1.

     F(1     ,:     ,1) = 0.
     F(NXmaxC,:     ,1) = 0.
     F(:     ,1     ,1) = 0.
     F(:     ,NYmaxC,1) = 1.

     F(:     ,1     ,2) = 0.
     F(:     ,NYmaxC,2) = 0.
     F(1     ,:     ,2) = 0.
     F(NXmaxC,:     ,2) = 0.


!------------------------------------------------

Return
End
