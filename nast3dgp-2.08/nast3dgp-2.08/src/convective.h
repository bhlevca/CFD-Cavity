/* NaSt3DGP - The Parallel 3D Navier-Stokes Solver
 * Copyright (C) 2003 Institute for Numerical Simulation
 *                    University of Bonn
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 */

#define DDP(L,M,R,DA,ii) (S.ddPstar[DA][0][ii]*L+S.ddPstar[DA][1][ii]*M+S.ddPstar[DA][2][ii]*R)
#define DD(L,M,R,DA,ii)  (S.ddstar[DA][0][ii]*L+S.ddstar[DA][1][ii]*M+S.ddstar[DA][2][ii]*R)
#define DDS(L,M,R,DA,ii) (S.ddSstar[DA][0][ii]*L+S.ddSstar[DA][1][ii]*M+S.ddSstar[DA][2][ii]*R)
#define cmp(i,j,k,a,b) if(a!=b) printf("%2d %2d %2d %g %g\n",i,j,k,a,b);

