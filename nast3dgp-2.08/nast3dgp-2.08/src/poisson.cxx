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
#include <math.h>
#include <stdlib.h>

#define DDP(L,M,R,DA,i) (S.ddPstar[DA][0][i]*L+S.ddPstar[DA][1][i]*M+S.ddPstar[DA][2][i]*R)


NS_REAL  NavierCalc::CompPoisson(double starttime){
 
  if (Par->IsParent() && (printinfo>=2)) {
	printf("\nStep:  %6d Time: %4.6f TimeStep: %2.14f ", S.timesteps, S.t,S.delt);fflush(stdout); 
	double actualtime=Par->GetTime();
	double timerd=actualtime-starttime ;
	int    hours =((int)timerd)/3600 ; timerd -=3600*hours ;
	int    mints =((int)timerd)/60   ; timerd -=  60*mints ;
	printf(" Time: %3dh %2dm %10.2fs\n",hours,mints,timerd) ; fflush(stdout);
  }
  
  return Poisson();                            // Solve Pressure equation  
}


//
// compute residual \| \Delta u-RHS \|_* (discrete)
//
 
NS_REAL NavierCalc::CalcRes(NS_REAL ijkroot) {
    int oldtimer=timer->Start(10,Par);
    NS_REAL    res=0.0, sq, balance=0.0;
    NS_REAL rhoXR,rhoXL,rhoYR,rhoYL,rhoZR,rhoZL;
    int     i,j,k;

    PIJKLOOP IFFLUID(flag[i][j][k]){                                                
	SetPLocalBorder(P,i,j,k);
	
	sq=DDP(P[i-1][j][k],P[i][j][k],P[i+1][j][k],0,i) +                            
	    DDP(P[i][j-1][k],P[i][j][k],P[i][j+1][k],1,j) +                            
	    DDP(P[i][j][k-1],P[i][j][k],P[i][j][k+1],2,k) - RHS[i][j][k];
	res+=sq*sq;
	if(printinfo>=4) balance += sq*DX[i]*DY[j]*DZ[k] ;
    }
    
    
    if (printinfo>=4) {                 //    Calculate mass balance 
	balance=Par->CommSum(balance);
	if (Par->IsParent()) {printf ("Mass Balance (poisson) %2.10e   \n", balance);fflush(stdout);}
    }
    
    timer->Stop(oldtimer,Par) ;
    return (sqrt(Par->CommSum(res)*ijkroot)/*S.delt*/);
}

inline void NavierCalc::PSOR(int i,int j,int k,NS_REAL omg1,NS_REAL lomg) {
     
    // Relaxiertes Gaus-Seidel-Verfahren
    //
    IFFLUID(flag[i][j][k]) {
	SetPLocalBorder(P,i,j,k);
	P[i][j][k]=omg1*P[i][j][k] - lomg/(S.ddPstar[0][1][i]+S.ddPstar[1][1][j]+S.ddPstar[2][1][k]) *  
	    (S.ddPstar[0][2][i]*P[i+1][j][k] + S.ddPstar[0][0][i]*P[i-1][j][k] + 
	     S.ddPstar[1][2][j]*P[i][j+1][k] + S.ddPstar[1][0][j]*P[i][j-1][k] +
	     S.ddPstar[2][2][k]*P[i][j][k+1] + S.ddPstar[2][0][k]*P[i][j][k-1] 
	     - RHS[i][j][k]);
    }
}

//
// solve linear system "\Delta p=RHS" using SOR iteration in a 8-color scheme
// boundary values are set before each SOR iteration 
// iteration is terminated when 
//   - the residue is too small 
//   - a maximum value of iterations is exceeded
// if SSOR is defined, symmetric SOR iteration is used. 
//

NS_REAL NavierCalc::Poisson() {
    int oldtimer=timer->Start(11,Par);
    NS_REAL    res;
    // Use BICGStab-Solver
    if (S.Solver == 5 || S.Solver == 6){
      res=BiCGStab();
    }else{
      // Use SOR-Solver-Types
      int     i, j, k, it=0, rescount=0, printcount=0;
      NS_REAL    ijkroot =1./(S.gridp[0]*S.gridp[1]*S.gridp[2]-S.ObstacleCount);
      NS_REAL    omg1    =1-S.omg;
      //    NS_REAL    oldres  =100*S.eps;
      NS_REAL    lomg    =S.omg;
      
      if(S.periodbound) {
	  if (S.Solver == 0 || S.Solver == 1 || S.Solver == 2)// Simple SOR or Red-Black Scheme
	      Par->CommMatrix(P, 0, 0, 0, PCOMM, S.GC); 
		
	  if (S.Solver == 3 || S.Solver == 4) {  // 7 Color-Scheme
	      Par->SendColour(P,0,0,0);
	      Par->SendColour(P,0,0,1);
	      Par->SendColour(P,0,1,0);
	      Par->SendColour(P,0,1,1);
	      Par->SendColour(P,1,0,0);
	      Par->SendColour(P,1,0,1);
	      Par->SendColour(P,1,1,0);
	      Par->SendColour(P,1,1,1);
	      Par->ReceiveColour(P,0,0,0);
	      Par->ReceiveColour(P,0,0,1);
	      Par->ReceiveColour(P,0,1,0);
	      Par->ReceiveColour(P,0,1,1);
	      Par->ReceiveColour(P,1,0,0);
	      Par->ReceiveColour(P,1,0,1);
	      Par->ReceiveColour(P,1,1,0);
	      Par->ReceiveColour(P,1,1,1);
	  } 
      }
      
      if((Par->IsParent()) && (printinfo>=4)) printf("  it:     %3d  ",it);
      res=CalcRes(ijkroot);
      
      //###################################################################
      //#   Begin of Poisson Loop for SOR,SSOR,RED BLACK or 8-COLOR SOR   #
      //###################################################################
      do {                
	if(S.Solver == 0 || S.Solver == 1 ){ // Simple SOR
	  //
	  //  SIMPLE PARALLELISATION WITHOUT COLORS
	  //
	    it++;
	    PIJKLOOP PSOR(i,j,k,omg1,lomg);
	    Par->CommMatrix(P, 0, 0, 0, PCOMM, S.GC);
	}
	
	if(S.Solver == 2){	   
	  //
	  //  RED-BLACK PARALLELISATION
	  //
	    it++;
	    // RED
	    PILOOP PJLOOP for(k=(Par->kug)+((i+j+(Par->kug))%2); k<=Par->kog; k+=2) PSOR(i,j,k,omg1,lomg);
	    Par->CommMatrix(P, 0, 0, 0, PCOMM, S.GC);	    
	    // BLACK
	    PILOOP PJLOOP for(k=(Par->kug)+((i+j+(Par->kug)+1)%2); k<=Par->kog; k+=2) PSOR(i,j,k,omg1,lomg);
	    Par->CommMatrix(P, 0, 0, 0, PCOMM, S.GC);
	}
	
	if(S.Solver == 3 || S.Solver == 4){
	  //
	  //    PARALLELISATION WITH COLOURS
	  //
	  //    ColorGroups:   - X    | Y    / Z
	  //      7   3   7
	  //    1 4 5 8 1 4
	  //    6 7 2 3 6 7
	  //    1   5   
	    it++;
	    if(rescount>0) Par->ReceiveColour(P,1,0,0);                         // 5
	    if(rescount>0) Par->ReceiveColour(P,0,1,0);                         // 6
	    if(rescount>0) Par->ReceiveColour(P,0,0,1);                         // 7
	    if(rescount>0) Par->ReceiveColour(P,1,1,1);                         // 8
//	    SetPBorder(P,0|0|0);
	    IALB JALB KALB PSOR(i,j,k,omg1,lomg);
	    Par->SendColour(P,0,0,0);                                           // 1
//	    SetPBorder(P,4|2|0);
	    IBLB JBLB KALB PSOR(i,j,k,omg1,lomg);
	    Par->SendColour(P,1,1,0);                                           // 2
//	    SetPBorder(P,4|0|1);
	    IBLB JALB KBLB PSOR(i,j,k,omg1,lomg);
	    Par->SendColour(P,1,0,1);                                           // 3
//	    SetPBorder(P,0|2|1);
	    IALB JBLB KBLB PSOR(i,j,k,omg1,lomg);
	    Par->SendColour(P,0,1,1);                                           // 4
	    
	    Par->ReceiveColour(P,0,0,0);                                        // 1
	    Par->ReceiveColour(P,1,1,0);                                        // 2
	    Par->ReceiveColour(P,1,0,1);                                        // 3
	    Par->ReceiveColour(P,0,1,1);                                        // 4
	    
//	    SetPBorder(P,4|0|0);
	    IBLB JALB KALB PSOR(i,j,k,omg1,lomg);
	    Par->SendColour(P,1,0,0);                                           // 5
//	    SetPBorder(P,0|2|0);
	    IALB JBLB KALB PSOR(i,j,k,omg1,lomg);
	    Par->SendColour(P,0,1,0);                                           // 6
//	    SetPBorder(P,0|0|1);
	    IALB JALB KBLB PSOR(i,j,k,omg1,lomg);
	    Par->SendColour(P,0,0,1);                                           // 7
//	    SetPBorder(P,4|2|1);
	    IBLB JBLB KBLB PSOR(i,j,k,omg1,lomg);
	    Par->SendColour(P,1,1,1);                                           // 8
	}
	
	
	if(S.Solver == 1){
	  //
	  //  SIMPLE SSOR PARALLELISATION WITHOUT COLOURS
	  //
	    for(i=Par->iog; i<=Par->iug; i--)
		for(j=Par->jog; j<=Par->jug; j--)
		    for(k=Par->kog; k<=Par->kug; k--) PSOR(i,j,k,omg1,lomg);
	    
	    Par->CommMatrix(P, 0, 0, 0, PCOMM, S.GC); 
	}
	
	if(S.Solver == 4){
	  // 
	  // 8 COLOR SSOR PARALLELISATION
	  //
	  Par->ReceiveColour(P,1,0,0);                                        // 5
	  Par->ReceiveColour(P,0,1,0);                                        // 6
	  Par->ReceiveColour(P,0,0,1);                                        // 7
	  Par->ReceiveColour(P,1,1,1);                                        // 8
//	  SetPBorder(P,0|0|0);
	  IAL JAL KAL PSOR(i,j,k,omg1,lomg);
	  Par->SendColour(P,0,0,0);                                           // 1
//	  SetPBorder(P,4|2|0);
	  IBL JBL KAL PSOR(i,j,k,omg1,lomg);
	  Par->SendColour(P,1,1,0);                                           // 2
//	  SetPBorder(P,4|0|1);
	  IBL JAL KBL PSOR(i,j,k,omg1,lomg);
	  Par->SendColour(P,1,0,1);                                           // 3
//	  SetPBorder(P,0|2|1);
	  IAL JBL KBL PSOR(i,j,k,omg1,lomg);
	  Par->SendColour(P,0,1,1);                                           // 4
	  
	  Par->ReceiveColour(P,0,0,0);                                        // 1
	  Par->ReceiveColour(P,1,1,0);                                        // 2
	  Par->ReceiveColour(P,1,0,1);                                        // 3
	  Par->ReceiveColour(P,0,1,1);                                        // 4
	  
//	  SetPBorder(P,4|0|0);
	  IBL JAL KAL PSOR(i,j,k,omg1,lomg);
	  Par->SendColour(P,1,0,0);                                           // 5
//	  SetPBorder(P,0|2|0);
	  IAL JBL KAL PSOR(i,j,k,omg1,lomg);
	  Par->SendColour(P,0,1,0);                                           // 6
//	  SetPBorder(P,0|0|1);
	  IAL JAL KBL PSOR(i,j,k,omg1,lomg);
	  Par->SendColour(P,0,0,1);                                           // 7
//	  SetPBorder(P,4|2|1);
	  IBL JBL KBL PSOR(i,j,k,omg1,lomg);
	  Par->SendColour(P,1,1,1);                                           // 8
	}

	if((--rescount==0) || (it>=S.itermax) || ((it<10) && (iterations<25))) {
	    
	    if(S.Solver == 4){
		Par->ReceiveColour(P, 1, 0, 0);
		Par->ReceiveColour(P, 0, 1, 0);
		Par->ReceiveColour(P, 0, 0, 1);
		Par->ReceiveColour(P, 1, 1, 1); 
	    }
	    //	    oldres=res;
	    if((Par->IsParent()) && (printinfo>=4)){printf("  it:     %3d  ",it);}
	    
	    res=CalcRes(ijkroot);
	} 
	if (rescount==-1) rescount+=10 ; 
	
        if((Par->IsParent())) if((printinfo>=3) && (printcount==0)) {
	    printf("  it: %7d  res: %e  omg: %1.2f  itmax: %d\n", it, res,lomg,S.itermax);fflush(stdout);
	    printcount=50 ; 
        }
        if((Par->IsParent())) printcount-- ;
	
      } while ((res>=S.eps) && (it<S.itermax)) ;      // End of  Poisson Loop
      
      iterations=it;

      if(S.periodbound) {
	if(S.Solver == 0 || S.Solver == 1 || S.Solver == 2){
	  Par->CommMatrix(P, 0, 0, 0, PCOMM, S.GC);  
	}
	if( S.Solver ==  3 || S.Solver ==  4){
	  Par->SendColour(P,0,0,0);
	  Par->SendColour(P,0,0,1);
	  Par->SendColour(P,0,1,0);
	  Par->SendColour(P,0,1,1);
	  Par->SendColour(P,1,0,0);
	  Par->SendColour(P,1,0,1);
	  Par->SendColour(P,1,1,0);
	  Par->SendColour(P,1,1,1);
	  Par->ReceiveColour(P,0,0,0);
	  Par->ReceiveColour(P,0,0,1);
	  Par->ReceiveColour(P,0,1,0);
	  Par->ReceiveColour(P,0,1,1);
	  Par->ReceiveColour(P,1,0,0);
	  Par->ReceiveColour(P,1,0,1);
	  Par->ReceiveColour(P,1,1,0);
	  Par->ReceiveColour(P,1,1,1);  
	}
      }
      alliterations+=iterations;
    }
    timer->Stop(oldtimer,Par);
    return (res);
}


void  NavierCalc::CalcKoeffPreCond(){
    int i,j,k;
    
    // Jacobi Preconditioner coeff
    PIJKLOOP IFFLUID(flag[i][j][k]){
	KoeffPreCond[i][j][k]= -1.0/(S.ddPstar[0][1][i]+S.ddPstar[1][1][j]+S.ddPstar[2][1][k]);
    }
}



inline NS_REAL NavierCalc::Precondition(NS_REAL P, int i, int j, int k){
  return KoeffPreCond[i][j][k] * P;
}

inline NS_REAL NavierCalc::PreconditionBack(NS_REAL P, int i, int j, int k){
  return 1.0/KoeffPreCond[i][j][k] * P;
}

inline NS_REAL  NavierCalc::MatVecMult(Matrix<NS_REAL>& A,int i, int j, int k, int RES){
  NS_REAL result;

  result= S.ddPstar[0][2][i]*A[i+1][j][k] + S.ddPstar[0][0][i]*A[i-1][j][k] + 
      S.ddPstar[1][2][j]*A[i][j+1][k] + S.ddPstar[1][0][j]*A[i][j-1][k] +
      S.ddPstar[2][2][k]*A[i][j][k+1] + S.ddPstar[2][0][k]*A[i][j][k-1] +
      (S.ddPstar[0][1][i]+S.ddPstar[1][1][j]+S.ddPstar[2][1][k]) * A[i][j][k];
  
  if(RES) result=RHS[i][j][k]-result;
  
  return result;
}


NS_REAL NavierCalc::BiCGStab() {
    NS_REAL a,b,w1,w2,w,res=1E9,norm_vj,norm_r0,norm_sj,norm_rj ;
    NS_REAL rj1r0, rjr0, sig ;
    int i,j,k;
    int it=0;
    NS_REAL ijkroot = 1./(S.gridp[0]*S.gridp[1]*S.gridp[2]-S.ObstacleCount);
    
    CalcKoeffPreCond();

 restart:
    rjr0=norm_r0=0.0 ;
    PIJKLOOP IFFLUID(flag[i][j][k]) {
	SetPLocalBorder(P,i,j,k);
	BiCG_r0[i][j][k]=BiCG_pj[i][j][k]=BiCG_rj[i][j][k]=MatVecMult(P,i,j,k,1);
	rjr0 += BiCG_rj[i][j][k]*BiCG_r0[i][j][k] ;
    }
    rjr0   =Par->CommSum(rjr0) ;
    norm_r0=sqrt(rjr0);
    PIJKLOOP IFFLUID(flag[i][j][k]) P[i][j][k]=PreconditionBack(P[i][j][k],i,j,k);  //  PRECONDITIONBACK P
    Par->CommMatrix(P, 0, 0, 0, PCOMM, S.GC);
    
    if((res>=S.eps) && (it<S.itermax)){
	
	do{
	    sig=0.0;
	    norm_vj=0.0;
	    norm_rj=0.0;
  
	    // APpj      
	    PIJKLOOP IFFLUID(flag[i][j][k]) MatBuf[i][j][k]=Precondition(BiCG_pj[i][j][k],i,j,k);  //  PRECONDITION pj 
	    Par->CommMatrix(MatBuf, 0, 0, 0, PCOMM, S.GC);
	    PIJKLOOP IFFLUID(flag[i][j][k]){
		SetPLocalBorder(MatBuf,i,j,k);
		BiCG_vj[i][j][k]=MatVecMult(MatBuf,i,j,k,0);  	
		sig   += BiCG_vj[i][j][k]*BiCG_r0[i][j][k];
		norm_vj += BiCG_vj[i][j][k]*BiCG_vj[i][j][k];
		norm_rj += BiCG_rj[i][j][k]*BiCG_rj[i][j][k]; // NEU
	    }
	    sig   = Par->CommSum(sig);
	    norm_vj = sqrt(Par->CommSum(norm_vj));
	    norm_rj = sqrt(Par->CommSum(norm_rj));  // NEU
	    a=rjr0 / sig ;
	    
	    if(fabs(sig) <= (1e-12*(norm_vj*norm_r0))){
		PIJKLOOP IFFLUID(flag[i][j][k]) P[i][j][k]=Precondition(P[i][j][k],i,j,k);  //  PRECONDITION P
		Par->CommMatrix(P, 0, 0, 0, PCOMM, S.GC);
		res=CalcRes(ijkroot); it++;
		if((Par->IsParent()) && (printinfo>=3) && (it%50==1)) {
		    printf("  it: %7d  res: %e    itmax: %d\n", it, res, S.itermax);fflush(stdout);
		}
		goto restart ;
	    }
	    if((fabs(a)*norm_vj/(norm_rj==0?1E-15:norm_rj))<=0.08){
		PIJKLOOP IFFLUID(flag[i][j][k]) P[i][j][k]=Precondition(P[i][j][k],i,j,k);  //  PRECONDITION P
		Par->CommMatrix(P, 0, 0, 0, PCOMM, S.GC);
		res=CalcRes(ijkroot); it++;
		if((Par->IsParent()) && (printinfo>=3) && (it%50==1)) {
		    printf("  it: %7d  res: %e    itmax: %d\n", it, res, S.itermax);fflush(stdout);
		}
		goto restart;  // NEU
	    }
	    
	    // sj =rj -a*vj
	    PIJKLOOP IFFLUID(flag[i][j][k]) 
		BiCG_sj[i][j][k] = BiCG_rj[i][j][k]-a*BiCG_vj[i][j][k] ;
      
	    norm_sj=0.0;
	    PIJKLOOP IFFLUID(flag[i][j][k]) norm_sj += BiCG_sj[i][j][k]*BiCG_sj[i][j][k];
	    norm_sj=sqrt(Par->CommSum(norm_sj));
	    
	    if(norm_sj>S.eps){
		// tj = Asj 
		PIJKLOOP IFFLUID(flag[i][j][k]) MatBuf[i][j][k]=Precondition(BiCG_sj[i][j][k],i,j,k);  //  PRECONDITION sj
		Par->CommMatrix(MatBuf, 0, 0, 0, PCOMM, S.GC);
		PIJKLOOP IFFLUID(flag[i][j][k]){
		    SetPLocalBorder(MatBuf,i,j,k);		    
		    BiCG_tj[i][j][k]=MatVecMult(MatBuf,i,j,k,0); 	
		}
		w1=w2=0.0;
		PIJKLOOP IFFLUID(flag[i][j][k]){
		    w1 += BiCG_tj[i][j][k]*BiCG_sj[i][j][k];
		    w2 += BiCG_tj[i][j][k]*BiCG_tj[i][j][k];
		}
		w1=Par->CommSum(w1);
		w2=Par->CommSum(w2);
		w=w1/(w2==0?1E-15:w2);
		
		// p += a*pj + w*sj 
		PIJKLOOP IFFLUID(flag[i][j][k]) 
		    P[i][j][k] += a*BiCG_pj[i][j][k] + w*BiCG_sj[i][j][k];
		
		// rj+1 = sj -w*tj
		PIJKLOOP IFFLUID(flag[i][j][k]) 
		    BiCG_rj[i][j][k]=BiCG_sj[i][j][k]-w*BiCG_tj[i][j][k] ;
		
		// < rj+1, r0 >  
		rj1r0=0.0;
		PIJKLOOP IFFLUID(flag[i][j][k])
		    rj1r0 += BiCG_rj[i][j][k]*BiCG_r0[i][j][k] ;
		
		rj1r0 = Par->CommSum(rj1r0);
		b=a*rj1r0/(w*rjr0==0?1E-15:(w*rjr0));
		
		PIJKLOOP IFFLUID(flag[i][j][k]) 
		    BiCG_pj[i][j][k] = BiCG_rj[i][j][k] + b*(BiCG_pj[i][j][k]-w*BiCG_vj[i][j][k]);
		
	    }else{
		PIJKLOOP IFFLUID(flag[i][j][k]) P[i][j][k] += a*BiCG_pj[i][j][k];
		
		rj1r0=0.0 ;
		PIJKLOOP IFFLUID(flag[i][j][k]) {
		    BiCG_rj[i][j][k]=BiCG_sj[i][j][k];
		    rj1r0 += BiCG_rj[i][j][k]*BiCG_r0[i][j][k] ;
		}
		rj1r0=Par->CommSum(rj1r0);
	    }
	    
	    // j=j+1 
	    rjr0 = rj1r0 ;
	    
	    res=0.0;
	    PIJKLOOP IFFLUID(flag[i][j][k]) res += BiCG_rj[i][j][k]*BiCG_rj[i][j][k]; 
	    res=sqrt(Par->CommSum(res)*ijkroot) /*S.delt*/;
	    
	    it++;
	    if((Par->IsParent()) && (printinfo>=3) && (it%50==1)) {
		printf("  it: %7d  res: %e    itmax: %d\n", it, res, S.itermax);fflush(stdout);
	    } 
	    
	}while( (res>=S.eps) && (it<S.itermax) ); 
	
    } // End If
    
    
    PIJKLOOP IFFLUID(flag[i][j][k]) P[i][j][k]=Precondition(P[i][j][k],i,j,k);  //  PRECONDITION P
    Par->CommMatrix(P, 0, 0, 0, PCOMM, S.GC);
    res=CalcRes(ijkroot); 

    iterations = it ;   
    alliterations+=iterations;
    
    return res;
}


