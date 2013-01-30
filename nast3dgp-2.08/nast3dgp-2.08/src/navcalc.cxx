/* NaSt3DGP - The Parallel 3D Navier-Stokes Solver
 * Copyright (C) 2003 Institute for Numerical Simulation
 *                    University of Bonn
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 */

#include "typen.h"
#include "matrix.h"
#include "parallel.h"
#include "navier.h"
#include "convective.h"
#include <math.h>
#include <errno.h>
#include <stdlib.h>
#include <cstring>

//##########################################################
//## Compute auxiliary-velocities without free-boundaries ##
//##########################################################
void NavierCalc::CompFGH() {
    int oldtimer=timer->Start(4,Par);
    int          i,j,k,border;
    NS_REAL      du2dx, duvdy, duwdz, d2udx2, d2udy2, d2udz2;
    NS_REAL      dv2dy, dvudx, dvwdz, d2vdx2, d2vdy2, d2vdz2;
    NS_REAL      dw2dz, dwudx, dwvdy, d2wdx2, d2wdy2, d2wdz2,a=0.0,b=0.0;
    NS_REAL      Gx, Gy, Gz;
    Gx = S.g[0]/S.froude; 
    Gy = S.g[1]/S.froude; 
    Gz = S.g[2]/S.froude; 



    
    // ############################################################################################################
    // ###                 <=5 point stencils (VONOS, SMART, HLPA, QUICK, DC)                                    ##               
    // ############################################################################################################
    assert(S.GC<=2) ;
    
    PJLOOP PKLOOP for (i =(((Par->iug==1) && !(S.periodbound&1)) ? 1:(Par->iug-1)); i<=Par->iog; i++)        
	IFFLUID(flag[i][j][k]) IFFLUID(flag[i+1][j][k]) { 
	//
	// Attention: current test for fluid cells assumes that number of ghostcells S.GC <=2 !!!! 
	// for larger S.GC a loop has to be implemented
	//
	if(ISFLUID(flag[(i-(S.GC-1))][j][k]) && ISFLUID(flag[(i+S.GC)][j][k])) 
	{border=FALSE;a=U[i-S.GC][j][k]; b=U[i+S.GC][j][k];} else {border=TRUE;a=0;b=0;} 
	du2dx= DUU(a,U[i-1][j][k],U[i][j][k],U[i+1][j][k],b,DX[i-1],DX[i],DX[i+1],DX[i+2],S.alpha,border);
	
	if((ISFLUID(flag[i][j-1][k]) && ISFLUID(flag[i][j+1][k])))
	{border=FALSE;a=U[i][j-S.GC][k]; b=U[i][j+S.GC][k];} else {border=TRUE; a=0;b=0;}
	duvdy= DUV(a,U[i][j-1][k],U[i][j][k],U[i][j+1][k],b,
		   (DX[i+1]*V[i][j-1][k]+DX[i]*V[i+1][j-1][k])/(DX[i]+DX[i+1]),
		   (DX[i+1]*V[i][j][k]+DX[i]*V[i+1][j][k])/(DX[i]+DX[i+1]),
		   DY[j-2],DY[j-1],DY[j],DY[j+1],DY[j+2],S.alpha,border);
	
	if((ISFLUID(flag[i][j][k-1]) && ISFLUID(flag[i][j][k+1])))
	{border=FALSE;a=U[i][j][k-S.GC]; b=U[i][j][k+S.GC];} else {border=TRUE;a=0;b=0;}
	duwdz= DUV(a,U[i][j][k-1],U[i][j][k],U[i][j][k+1],b,
		   (DX[i+1]*W[i][j][k-1]+DX[i]*W[i+1][j][k-1])/(DX[i]+DX[i+1]),
		   (DX[i+1]*W[i][j][k]+DX[i]*W[i+1][j][k])/(DX[i]+DX[i+1]),
		   DZ[k-2],DZ[k-1],DZ[k],DZ[k+1],DZ[k+2],S.alpha,border);  
	
	d2udx2=DD (U[i-1][j][k],U[i][j][k],U[i+1][j][k],0,i);
	d2udy2=DDS(U[i][j-1][k],U[i][j][k],U[i][j+1][k],1,j);
	d2udz2=DDS(U[i][j][k-1],U[i][j][k],U[i][j][k+1],2,k);
	F[actual_fgh][i][j][k] = S.nu*(d2udx2+d2udy2+d2udz2) - du2dx-duvdy-duwdz + Gx;
    }
	 
	 
    PILOOP PKLOOP for(j =(((Par->jug==1) && !(S.periodbound&2)) ? 1:(Par->jug-1)); j<=Par->jog; j++)
	IFFLUID(flag[i][j][k]) IFFLUID(flag[i][j+1][k]) {
	if(ISFLUID(flag[i][(j-(S.GC-1))][k]) && ISFLUID(flag[i][(j+S.GC)][k]))
	{border=FALSE;a=V[i][j-S.GC][k]; b=V[i][j+S.GC][k];} else {border=TRUE;a=0;b=0;}
	dv2dy= DUU(a,V[i][j-1][k],V[i][j][k],V[i][j+1][k],b,DY[j-1],DY[j],DY[j+1],DY[j+2],S.alpha,border);
	
	if((ISFLUID(flag[i-1][j][k]) && ISFLUID(flag[i+1][j][k])))
	{border=FALSE;a=V[i-S.GC][j][k]; b=V[i+S.GC][j][k];}  else {border=TRUE;a=0;b=0;}
	dvudx= DUV(a,V[i-1][j][k],V[i][j][k],V[i+1][j][k],b,
		   (DY[j+1]*U[i-1][j][k]+DY[j]*U[i-1][j+1][k])/(DY[j]+DY[j+1]),
		   (DY[j+1]*U[i][j][k]+DY[j]*U[i][j+1][k])/(DY[j]+DY[j+1]),
		   DX[i-2],DX[i-1],DX[i],DX[i+1],DX[i+2],S.alpha,border);
	
	if((ISFLUID(flag[i][j][k-1]) && ISFLUID(flag[i][j][k+1]))) 
	{border=FALSE;a=V[i][j][k-S.GC]; b=V[i][j][k+S.GC];} else {border=TRUE; a=0;b=0;}
	dvwdz= DUV(a,V[i][j][k-1],V[i][j][k],V[i][j][k+1],b,
		   (DY[j+1]*W[i][j][k-1]+DY[j]*W[i][j+1][k-1])/(DY[j]+DY[j+1]),
		   (DY[j+1]*W[i][j][k]+DY[j]*W[i][j+1][k])/(DY[j]+DY[j+1]),
		   DZ[k-2],DZ[k-1],DZ[k],DZ[k+1],DZ[k+2],S.alpha,border);
	
	d2vdx2=DDS(V[i-1][j][k],V[i][j][k],V[i+1][j][k],0,i);
	d2vdy2=DD (V[i][j-1][k],V[i][j][k],V[i][j+1][k],1,j);
	d2vdz2=DDS(V[i][j][k-1],V[i][j][k],V[i][j][k+1],2,k);
	G[actual_fgh][i][j][k] = S.nu*(d2vdx2+d2vdy2+d2vdz2) - dv2dy-dvudx-dvwdz + Gy;
    }
    
    
    PILOOP PJLOOP for(k =(((Par->kug==1) && !(S.periodbound&4)) ? 1:(Par->kug-1)); k<=Par->kog; k++)
	IFFLUID(flag[i][j][k]) IFFLUID(flag[i][j][k+1]) {
	if(ISFLUID(flag[i][j][(k-(S.GC-1))]) && ISFLUID(flag[i][j][(k+S.GC)]))
	{border=FALSE;a=W[i][j][k-S.GC]; b=W[i][j][k+S.GC];}  else {border=TRUE; a=0;b=0;}      
	dw2dz= DUU(a,W[i][j][k-1],W[i][j][k],W[i][j][k+1],b,DZ[k-1],DZ[k],DZ[k+1],DZ[k+2],S.alpha,border);
	
	if((ISFLUID(flag[i-1][j][k]) && ISFLUID(flag[i+1][j][k])))
	{border=FALSE;a=W[i-S.GC][j][k]; b=W[i+S.GC][j][k];} else {border=TRUE; a=0;b=0;}
	dwudx= DUV(a,W[i-1][j][k],W[i][j][k],W[i+1][j][k],b,
		   (DZ[k+1]*U[i-1][j][k]+DZ[k]*U[i-1][j][k+1])/(DZ[k]+DZ[k+1]),
		   (DZ[k+1]*U[i][j][k]+DZ[k]*U[i][j][k+1])/(DZ[k]+DZ[k+1]),
		   DX[i-2],DX[i-1],DX[i],DX[i+1],DX[i+2],S.alpha,border);
	
	if((ISFLUID(flag[i][j-1][k]) && ISFLUID(flag[i][j+1][k])))
	{border=FALSE;a=W[i][j-S.GC][k]; b=W[i][j+S.GC][k];} else {border=TRUE;a=0;b=0;}
	dwvdy= DUV(a,W[i][j-1][k],W[i][j][k],W[i][j+1][k],b,
		   (DZ[k+1]*V[i][j-1][k]+DZ[k]*V[i][j-1][k+1])/(DZ[k]+DZ[k+1]),
		   (DZ[k+1]*V[i][j][k]+DZ[k]*V[i][j][k+1])/(DZ[k]+DZ[k+1]),
		   DY[j-2],DY[j-1],DY[j],DY[j+1],DY[j+2],S.alpha,border);
	
	d2wdx2=DDS(W[i-1][j][k],W[i][j][k],W[i+1][j][k],0,i);
	d2wdy2=DDS(W[i][j-1][k],W[i][j][k],W[i][j+1][k],1,j);
	d2wdz2=DD (W[i][j][k-1],W[i][j][k],W[i][j][k+1],2,k);
	H[actual_fgh][i][j][k] = S.nu*(d2wdx2+d2wdy2+d2wdz2) - dw2dz-dwudx-dwvdy + Gz;
    }
    
    
    if(S.CompTemp) {
	double g0=S.g[0]*S.beta, g1=S.g[1]*S.beta, g2=S.g[2]*S.beta ;
	
	if(S.g[0]!=0.0) PJLOOP PKLOOP
			    for(i =(((Par->iug==1) && !(S.periodbound&1)) ? 1:(Par->iug-1)); i<=Par->iog; i++)
				IFFLUID(flag[i][j][k]) IFFLUID(flag[i+1][j][k]) 
				    F[actual_fgh][i][j][k] -= g0*( (DX[i+1]*T[i][j][k]+DX[i]*T[i+1][j][k])/(DX[i+1]+DX[i]) - S.TempRef);
	
	if(S.g[1]!=0.0) PILOOP PKLOOP 
			    for(j =(((Par->jug==1) && !(S.periodbound&2)) ? 1:(Par->jug-1)); j<=Par->jog; j++) 
				IFFLUID(flag[i][j][k]) IFFLUID(flag[i][j+1][k]) 
				    G[actual_fgh][i][j][k] -= g1*( (DY[j+1]*T[i][j][k]+DY[j]*T[i][j+1][k])/(DY[j+1]+DY[j]) - S.TempRef);
	
	if(S.g[2]!=0.0) PILOOP PJLOOP 
			    for(k =(((Par->kug==1) && !(S.periodbound&4)) ? 1:(Par->kug-1)); k<=Par->kog; k++)
				IFFLUID(flag[i][j][k]) IFFLUID(flag[i][j][k+1]) 
				    H[actual_fgh][i][j][k] -= g2*( (DZ[k+1]*T[i][j][k]+DZ[k]*T[i][j][k+1])/(DZ[k+1]+DZ[k]) - S.TempRef);
    }
    timer->Stop(oldtimer,Par);
}



//############################################
//## time handling for auxiliary-velocities ##
//############################################
void NavierCalc::CompTUVWfromFGH(int timestepmethod) {  
    int i,j,k;

    if(S.TimeDis==0){
	//##############################
	//##     EULER 1st ORDER      ##
	//##############################
	PJLOOP PKLOOP for(i =(((Par->iug==1) && !(S.periodbound&1)) ? 1:(Par->iug-1)); i<=Par->iog; i++)
	    U[i][j][k]+=S.delt*F[actual_fgh][i][j][k];
	PILOOP PKLOOP for(j =(((Par->jug==1) && !(S.periodbound&2)) ? 1:(Par->jug-1)); j<=Par->jog; j++)
	    V[i][j][k]+=S.delt*G[actual_fgh][i][j][k];
	PILOOP PJLOOP for(k =(((Par->kug==1) && !(S.periodbound&4)) ? 1:(Par->kug-1)); k<=Par->kog; k++)
	    W[i][j][k]+=S.delt*H[actual_fgh][i][j][k];
    } else {
	//##############################
	//##   Adams-Bash. 2nd ORDER  ##
	//##############################
	switch(timestepmethod) {
	case 0:  // Euler 1st order
	    PJLOOP PKLOOP for(i =(((Par->iug==1) && !(S.periodbound&1)) ? 1:(Par->iug-1)); i<=Par->iog; i++) 
		U[i][j][k]+=S.delt*F[actual_fgh][i][j][k];
	    PILOOP PKLOOP for(j =(((Par->jug==1) && !(S.periodbound&2)) ? 1:(Par->jug-1)); j<=Par->jog; j++) 
		V[i][j][k]+=S.delt*G[actual_fgh][i][j][k];
	    PILOOP PJLOOP for(k =(((Par->kug==1) && !(S.periodbound&4)) ? 1:(Par->kug-1)); k<=Par->kog; k++) 
		W[i][j][k]+=S.delt*H[actual_fgh][i][j][k];
	    break;  
	case 1:  // Adams-Bashforth for non-equidistant time-steps
	    PJLOOP PKLOOP for(i =(((Par->iug==1) && !(S.periodbound&1)) ? 1:(Par->iug-1)); i<=Par->iog; i++)
		U[i][j][k]+=0.5*S.delt*((S.delt+2*S.delt_old)/S.delt_old*F[actual_fgh][i][j][k]-S.delt/S.delt_old*F[actual_fgh^1][i][j][k]);
	    PILOOP PKLOOP for(j =(((Par->jug==1) && !(S.periodbound&2)) ? 1:(Par->jug-1)); j<=Par->jog; j++)  
		V[i][j][k]+=0.5*S.delt*((S.delt+2*S.delt_old)/S.delt_old*G[actual_fgh][i][j][k]-S.delt/S.delt_old*G[actual_fgh^1][i][j][k]);
	    PILOOP PJLOOP for(k =(((Par->kug==1) && !(S.periodbound&4)) ? 1:(Par->kug-1)); k<=Par->kog; k++) 
		W[i][j][k]+=0.5*S.delt*((S.delt+2*S.delt_old)/S.delt_old*H[actual_fgh][i][j][k]-S.delt/S.delt_old*H[actual_fgh^1][i][j][k]);
	    break;
	default:
	    if (Par->IsParent()) printf("Error in Time Scheme! Use EU1 or AB2 in *.nav configuration file\n");  
	}
    }
}


//##################################################
//## Compute right hand side for poisson equation ##
//##################################################
void NavierCalc::CompRHS() {
    int oldtimer=timer->Start(5,Par);
    int i, j, k;

    PIJKLOOP IFFLUID(flag[i][j][k]) { 
            RHS[i][j][k]=
	      ((U[i][j][k]-U[i-1][j][k])/(DX[i])+
               (V[i][j][k]-V[i][j-1][k])/(DY[j])+
               (W[i][j][k]-W[i][j][k-1])/(DZ[k])) /* /S.delt */ ;
    }
    
    RHSbalance=0.0;                      //   calculate mass balance
    PIJKLOOP IFFLUID(flag[i][j][k]) RHSbalance+=RHS[i][j][k]*RHS[i][j][k];
    RHSbalance=sqrt(Par->CommSum(RHSbalance)/(S.gridp[0]*S.gridp[1]*S.gridp[2]-S.ObstacleCount));

    timer->Stop(oldtimer,Par);
}


//#########################################################
//## calc new velocities by pressure-gradient correction ##
//#########################################################
void NavierCalc::AdapUVW() {
    int oldtimer=timer->Start(2,Par);
    int i,j,k;

    PJLOOP PKLOOP for(i =(((Par->iug==1) && !(S.periodbound&1)) ? 1:(Par->iug-1)); i<=Par->iog; i++) 
	IFFLUID(flag[i][j][k]) IFFLUID(flag[i+1][j][k])
	    U[i][j][k]-=/*S.delt* */ (P[i+1][j][k]-P[i][j][k])/DXM[i+1];      
    PILOOP PKLOOP for(j =(((Par->jug==1) && !(S.periodbound&2)) ? 1:(Par->jug-1)); j<=Par->jog; j++) 
	IFFLUID(flag[i][j][k]) IFFLUID(flag[i][j+1][k])
	    V[i][j][k]-=/*S.delt* */ (P[i][j+1][k]-P[i][j][k])/DYM[j+1];
    PILOOP PJLOOP for(k =(((Par->kug==1) && !(S.periodbound&4)) ? 1:(Par->kug-1)); k<=Par->kog; k++) 
	IFFLUID(flag[i][j][k]) IFFLUID(flag[i][j][k+1])
	    W[i][j][k]-=/*S.delt* */ (P[i][j][k+1]-P[i][j][k])/DZM[k+1];
        
    timer->Stop(oldtimer,Par);
}


//########################################
//## 'solve' transport equation (space) ##
//########################################
void NavierCalc::CompTG(Matrix<NS_REAL>& M, NS_REAL alpha, NS_REAL diffconst){
    int oldtimer=timer->Start(3,Par);
    int     i, j, k, border;
    NS_REAL    duTdx, dvTdy, dwTdz, d2Tdx2, d2Tdy2, d2Tdz2, a, b;


    PIJKLOOP IFFLUID(flag[i][j][k]) {       // T2=F(T^n)
	if(ISFLUID(flag[(i-(S.GC-1))][j][k]) && ISFLUID(flag[(i+(S.GC-1))][j][k]))
	{border=FALSE; a=M[i-2][j][k]; b=M[i+2][j][k];} else {border=TRUE; a=0.0; b=0.0;}      
	duTdx=DUV(a,M[i-1][j][k],M[i][j][k],M[i+1][j][k],b,U[i-1][j][k],U[i][j][k],
		  DX[i-2],DX[i-1],DX[i],DX[i+1],DX[i+2],alpha,border);
	
	if(ISFLUID(flag[i][(j-(S.GC-1))][k]) &&  ISFLUID(flag[i][(j+(S.GC-1))][k]))
	{border=FALSE;a=M[i][j-2][k]; b=M[i][j+2][k];} else {border=TRUE;a=0.0; b=0.0;} 
	dvTdy=DUV(a,M[i][j-1][k],M[i][j][k],M[i][j+1][k],b,V[i][j-1][k],V[i][j][k],
		  DY[j-2],DY[j-1],DY[j],DY[j+1],DY[j+2],alpha,border);
	
	if(ISFLUID(flag[i][j][(k-(S.GC-1))]) && ISFLUID(flag[i][j][(k+(S.GC-1))]))
	{border=FALSE;a=M[i][j][k-2]; b=M[i][j][k+2];} else {border=TRUE;a=0.0; b=0.0;} 
	dwTdz=DUV(a,M[i][j][k-1],M[i][j][k],M[i][j][k+1],b,W[i][j][k-1],W[i][j][k],
		  DZ[k-2],DZ[k-1],DZ[k],DZ[k+1],DZ[k+2],alpha,border);
	
	d2Tdx2=DDS(M[i-1][j][k],M[i][j][k],M[i+1][j][k],0,i);
	d2Tdy2=DDS(M[i][j-1][k],M[i][j][k],M[i][j+1][k],1,j);
	d2Tdz2=DDS(M[i][j][k-1],M[i][j][k],M[i][j][k+1],2,k); 
	T2[i][j][k]=-duTdx-dvTdy-dwTdz+diffconst*(d2Tdx2+d2Tdy2+d2Tdz2);
    }

    timer->Stop(oldtimer,Par);
}


//#########################################################
//## 'solve' transport equation (time) for scalar values ##
//#########################################################
void NavierCalc::CompTT(Matrix<NS_REAL>& M, Matrix<NS_REAL>& FT, NS_REAL alpha, NS_REAL diffconst, int n, int timestepmethod){
    int     i, j, k;
   
    switch(S.TimeDis){
    case 0:
	//##############################
	//##     Euler 1st ORDER      ##
	//##############################
	PIJKLOOP IFFLUID(flag[i][j][k]) {
	    M[i][j][k]+=S.delt*T2[i][j][k];
	    FT[i][j][k]=T2[i][j][k];
	}
	break;
    case 1:
	//##############################
	//##   Adams-Bash. 2nd ORDER  ##
	//##############################
	switch(timestepmethod) {
	case 0:
	    PIJKLOOP IFFLUID(flag[i][j][k]) { // Euler
		M[i][j][k]+=S.delt*T2[i][j][k];
		FT[i][j][k]=T2[i][j][k];
	    }
	    break;
	case 1:
	    PIJKLOOP IFFLUID(flag[i][j][k]) { // Adams-Bashforth for non-equidistant time steps
		M[i][j][k]+=0.5*S.delt*((S.delt+2*S.delt_old)/S.delt_old*T2[i][j][k]-S.delt/S.delt_old*FT[i][j][k]);
		FT[i][j][k]=T2[i][j][k];
	    }
	    break;
	}
	break;	    
    default:
	if (Par->IsParent()) printf("Error in Time Scheme! Use EU1 or AB2 in *.nav configuration file\n");  
    }
}


//#######################
//## Time-Step Control ##
//#######################
NS_REAL NavierCalc::TimeStep() {
    int             oldtimer=timer->Start(6,Par)  ;
    NS_REAL         deltneu=S.deltmax ; 
    static NS_REAL  xmin=1e15,ymin=1e15,zmin=1e15 ;
    int             i,j,k,l;

    if(xmin==1e15 && ymin==1e15 && zmin==1e15) {                   
      PIJKLOOP {xmin=std::min(xmin,DX[i]); ymin=std::min(ymin,DY[j]); zmin=std::min(zmin,DZ[k]);} 
      xmin=Par->CommMin(xmin);
      ymin=Par->CommMin(ymin);
      zmin=Par->CommMin(zmin);
    }

    // local CFL-condition including volume forces
    PIJKLOOP {
      deltneu=std::min(deltneu, S.tfconv*2/(fabs(U[i][j][k])/DXM[i]+
				       sqrt( (U[i][j][k]*U[i][j][k])/(DXM[i]*DXM[i])+4*fabs(S.g[0])/(S.froude*DXM[i])) +1E-15));
      deltneu=std::min(deltneu, S.tfconv*2/(fabs(V[i][j][k])/DYM[j]+
				       sqrt( (V[i][j][k]*V[i][j][k])/(DYM[j]*DYM[j])+4*fabs(S.g[1])/(S.froude*DYM[j])) +1E-15));
      deltneu=std::min(deltneu, S.tfconv*2/(fabs(W[i][j][k])/DZM[k]+
				       sqrt( (W[i][j][k]*W[i][j][k])/(DZM[k]*DZM[k])+4*fabs(S.g[2])/(S.froude*DZM[k])) +1E-15));
    }
    
    // Viscosity
    deltneu=std::min(deltneu,S.tfdiff*S.re*0.5/(1./(xmin*xmin)+1./(ymin*ymin)+1./(zmin*zmin)));
    
 
    if(S.CompTemp) if(S.nu!=0.0 && S.prandtl!=0) // temperature
      deltneu=std::min(deltneu,S.tfdiff*S.re*S.prandtl*0.5*1./(1./(xmin*xmin)+1./(ymin*ymin)+1./(zmin*zmin)));
    if(S.CompChem) for(l=0; l<S.nchem; l++) if(S.chemc[l]!=0.0)  // scalars viscosity
      deltneu=std::min(deltneu,S.tfdiff/S.chemc[l]*0.5*1./(1./(xmin*xmin)+1./(ymin*ymin)+1./(zmin*zmin)));
 
    deltneu=Par->CommMin(deltneu); 
    if(deltneu<1E-15) {
      printf("Maybe ERROR in parameters: Time-Step = 0.0 ! \n");exit(1);    
    } 

    if(S.TimeDis==1) {      // Compute new time step (8 steps smoothing through parabola for AB)
      deltneu=(S.timesteps<9)?((1-((S.timesteps+1)*0.1-1)*((S.timesteps+1)*0.1-1))*deltneu):deltneu;
    }

    timer->Stop(oldtimer,Par);
    return deltneu;
}


//
// Computes mean value
//
NS_REAL NavierCalc::CalcMeanVal(Matrix<NS_REAL>& mat) {
    NS_REAL   sum=0,count=0;
    int       i,j,k;

    PIJKLOOP IFFLUID(flag[i][j][k]) {sum+=mat[i][j][k];count+=1.0;}
    sum=Par->CommSum(sum);
    count=Par->CommSum(count);
    return(sum/count);
}


void NavierCalc::InfoOutputs(NS_REAL res) {
    NS_REAL        APSbalance;
    int               i, j, k;

    if (Par->IsParent() && (printinfo>0)) {  // Some Info-outputs
	if(printinfo==1) printf("%d %f\n",S.timesteps,res); else {
	    printf("Iter: %7d  Res:%15.8e     all It: %8d\n",iterations,res,alliterations);
	    printf("Mass Balance (rhs): %2.8e\n",RHSbalance);
	    fflush(stdout);
	}
    }
    
    if(printinfo>=3 /*&& (S.iobd==0 || S.ObstacleCount==0)*/){
	PIJKLOOP IFFLUID(flag[i][j][k]) { 
	    RHS[i][j][k]=
		((U[i][j][k]-U[i-1][j][k])/(DX[i])+
		 (V[i][j][k]-V[i][j-1][k])/(DY[j])+
		 (W[i][j][k]-W[i][j][k-1])/(DZ[k])) /* /S.delt */ ;
	}
	
	APSbalance=0.0;             //   calculate mass balance after proj. step
	PIJKLOOP IFFLUID(flag[i][j][k]) APSbalance+=RHS[i][j][k]*RHS[i][j][k] /*DX[i]*DY[j]*DZ[k])*/;
	APSbalance=sqrt(Par->CommSum(APSbalance)/(S.gridp[0]*S.gridp[1]*S.gridp[2]-S.ObstacleCount)); 
	if(Par->IsParent()) printf("Mass Balance (aps): %2.8e\n",APSbalance); fflush(stdout);
    }
        
    if(Par->IsParent() && (printinfo>2) && S.iobd>=1) {
	printf("flow_in=%g  flow_out=%g  mass_diff=%g",flow_in, flow_out, mass_diff);fflush(stdout);
	if(S.iobd==2) printf("  corrected_mass_diff=%g",corrected_mass_diff);
	printf("\n");fflush(stdout);
    }
    
    if(printinfo>=4) { // info output routines
	int i,j,k;
	NS_REAL mval=0.0;
	PIJKALLLOOP IFFLUID(flag[i][j][k]) mval=std::max(mval,fabs(U[i][j][k]));
	mval=Par->CommMax(mval);
	if(Par->IsParent()) printf("u_max=%9.6f ",mval);
	mval=0.0;
	PIJKALLLOOP IFFLUID(flag[i][j][k]) mval=std::max(mval,fabs(V[i][j][k]));
	mval=Par->CommMax(mval);
	if(Par->IsParent()) printf("v_max=%9.6f ",mval);
	mval=0.0;
	PIJKALLLOOP IFFLUID(flag[i][j][k]) mval=std::max(mval,fabs(W[i][j][k]));
	mval=Par->CommMax(mval);
	if(Par->IsParent()) printf("w_max=%9.6f ",mval);
	mval=0.0;
	PIJKALLLOOP IFFLUID(flag[i][j][k]) mval=std::max(mval,fabs(P[i][j][k]));
	mval=Par->CommMax(mval);
	if(Par->IsParent()) printf("p_max=%9.6f\n",mval);
	mval=CalcMeanVal(U);
	if(Par->IsParent()) printf("u_mid=%9.6f ",mval);
	mval=CalcMeanVal(V);
	if(Par->IsParent()) printf("v_mid=%9.6f ",mval);
	mval=CalcMeanVal(W);
	if(Par->IsParent()) printf("w_mid=%9.6f ",mval);
	mval=CalcMeanVal(P);
	if(Par->IsParent()) printf("p_mid=%9.6f\n",mval);
	fflush(stdout);
    }
}


void NavierCalc::OutputData(){
    char* name = new char[255];
    char* targetfilename = new char[255];
    static int        q=1;
    static double     p=0;
    
    // File writing routine depending on prstep (printing after x timesteps) and
    // prtTime (printing after x physical time) check routines
    if(((S.timesteps%(S.prstep+1)==0)&&(S.timesteps>0)) || (S.t==S.Tfin) || (S.prtTime>0.0 && S.t>=p)) { // File writing
	if(Par->IsParent() && (printinfo>=2)) printf("Output...\n"); fflush(stdout);	    
	if(!WriteFile(binfile)) {fprintf(stderr,"Error writing file %s.\n",binfile);exit(1);}
	if(S.prtTime>0.0 && S.t>=p){
	    p+=S.prtTime;
	    if(Par->IsParent()){ 
		if (strrchr(binfile, '/') == NULL) 
		    strcpy(targetfilename, binfile);
		else
		    strcpy(targetfilename, strrchr(binfile, '/') + 1);		// don't copy the '/'
		sprintf(name,"cp %s %s/%03d.%s",binfile, S.targetdir, q, targetfilename);
		printf("Writing file \"%s\" to %s/%03d.%s ...\n",binfile, S.targetdir, q++, targetfilename);
		system(name);
	    } 
	}
    }
    delete [] name;
    delete [] targetfilename;
}
