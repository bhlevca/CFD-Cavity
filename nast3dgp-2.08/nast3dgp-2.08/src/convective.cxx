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


//####################################################################
//#                                                                  #
//# Discretization of convective Terms:                              #
//#                                                                  #
//#       e.g. discretization of \partial(uu)/\partial x             #
//#                                                                  #
//#    ------------- -------------- -------------- --------------    #
//#   |             |              |              |              |   #
//#   |             |              |              |              |   #
//#   |             |              |              |              |   #
//#   LL            L       KL     M      KR      R              RR  #
//#   |             |              |              |              |   #
//#   |             |              |              |              |   #
//#   |             |              |              |              |   #
//#    ------------- -------------- -------------- --------------    #
//#    \____  _____/ \_____  _____/ \______  ____/ \_____  _____/    #
//#         \/             \/              \/            \/          #
//#         DD             D               DP            PP          #
//#                                                                  #
//####################################################################

NS_REAL NavierCalc::DUU(NS_REAL LL, NS_REAL L, NS_REAL M, NS_REAL R, NS_REAL RR, NS_REAL DD, NS_REAL D, NS_REAL DP, NS_REAL PP, NS_REAL alpha, int border) {
    NS_REAL d, KL=(L+M)/2, KR=(R+M)/2, UR, UL, KRABS, KLABS, Phi_R, Phi_L, UUR, UUL, ksiR, ksiL, rR, rL, q1R, q2R, q1L, q2L; 
    
    if(S.Convective != DonorCell){ // if donor-cell is not defined
      UR=(KR>=0) ? 0:1;
      UL=(KL>=0) ? 0:1;
    }

    switch (S.Convective){
    case DonorCell:
      KRABS=fabs(KR);           // Donor-Cell
      KLABS=fabs(KL);
      d=(KR*(R+M) -KL*(L+M))/(D +DP); // central differences
      d=(1-alpha)*d+alpha*( (KR-KRABS)*R + ( (KR+KRABS)-(KL-KLABS) )*M - (KL+KLABS)*L )/(D+DP); 
      break;
    case QUICK:   
      if(border==FALSE){     // Quadratic upwind interpolation for convective kinematics (QUICK)       
	UR=(1-UR)*((3.*R+6.*M-L)/8.) + UR*((3.*M+6.*R-RR)/8.);
	UL=(1-UL)*((3.*M+6.*L-LL)/8.) + UL*((3.*L+6.*M-R)/8.);
	d=(KR*UR-KL*UL)/(0.5*(D +DP));
      }else{// Donor-Cell full upwind
	KRABS=fabs(KR);       
	KLABS=fabs(KL);
	d=(KR*(R+M) -KL*(L+M))/(D +DP); // central differences
	// alpha=0 => central differences,  alpha=1 => full upwind
	d=(1-alpha)*d + alpha*( (KR-KRABS)*R + ( (KR+KRABS)-(KL-KLABS) )*M - (KL+KLABS)*L )/(D+DP); 
      }  
      break;
    case HLPA:     
      if(border==FALSE){               // Hybrid-Linear Parabolic Appoximation (HLPA) 2nd-Order
	Phi_R=(1-UR)*(M-L)/(((R-L)!=0.0)?(R-L):1E-15)  + UR*(R-RR)/(((M-RR)!=0.0)?(M-RR):1E-15);
	Phi_L=(1-UL)*(L-LL)/(((M-LL)!=0.0)?(M-LL):1E-15) + UL*(M-R)/(((L-R)!=0.0)?(L-R):1E-15);
	
	if((Phi_R>1) || (Phi_R<0)) 
	  UR=(1-UR)*M+UR*R;
	else
	  UR=(1-UR)*(M+(R-M)*Phi_R) + UR*(R+(M-R)*Phi_R);
	
	if((Phi_L>1) || (Phi_L<0)) 
	  UL=(1-UL)*L+UL*M;
	else
	  UL=(1-UL)*(L+(M-L)*Phi_L) + UL*(M+(L-M)*Phi_L);
	
	d=(KR*UR-KL*UL)/(0.5*(D +DP));
      }else{// Donor-Cell full upwind
	KRABS=fabs(KR);        
	KLABS=fabs(KL);
	d=(KR*(R+M) -KL*(L+M))/(D +DP); // central differences
	// alpha=0 => central differences,  alpha=1 => full upwind
	d=(1-alpha)*d + alpha*( (KR-KRABS)*R + ( (KR+KRABS)-(KL-KLABS) )*M - (KL+KLABS)*L )/(D+DP); 
      }  
      break;
    case SMART:
      if(border==FALSE){         // Sharp and Monotonic Algorithm for Realistic Transport (SMART) 2nd-Order
	Phi_R=(1-UR)*(M-L)/(((R-L)!=0.0)?(R-L):1E-15)  + UR*(R-RR)/(((M-RR)!=0.0)?(M-RR):1E-15);
	Phi_L=(1-UL)*(L-LL)/(((M-LL)!=0.0)?(M-LL):1E-15) + UL*(M-R)/(((L-R)!=0.0)?(L-R):1E-15);
	
	if((Phi_R>=0) && (Phi_R<(3./74)))
	  UR=(1-UR)*(10*M-9*L)+UR*(10*R-9*RR);
	else if((Phi_R>=(3./74)) && (Phi_R<(5./6)))
	    UR=(1-UR)*((3.*R+6.*M-L)/8.) + UR*((3.*M+6.*R-RR)/8.);
//	  UR=(1-UR)*((3./8)*R+(6./8)*M-(1./8)*L)+UR*((3./8)*M+(6./8)*R-(1./8)*RR);
	else if((Phi_R>=(5./6)) && (Phi_R<=1))
	  UR=(1-UR)*R+UR*M;
	else 
	  UR=(1-UR)*M+UR*R;
	
	if((Phi_L>=0) && (Phi_L<(3./74)))
	  UL=(1-UL)*(10*L-9*LL)+UL*(10*M-9*R);
	else if((Phi_L>=(3./74)) && (Phi_L<(5./6)))
          UL=(1-UL)*((3.*M+6.*L-LL)/8.) + UL*((3.*L+6.*M-R)/8.);
          //UL=(1-UL)*((3./8)*M+(6./8)*L-(1./8)*LL)+UL*((3./8)*L+(6./8)*M-(1./8)*R);
	else if((Phi_L>=(5./6)) && (Phi_L<=1))
	  UL=(1-UL)*M+UL*L;
	else 
	  UL=(1-UL)*L+UL*M;
	
	d=(KR*UR-KL*UL)/(0.5*(D +DP));
      }else{// Donor-Cell 
	KRABS=fabs(KR);       
	KLABS=fabs(KL);
	d=(KR*(R+M) -KL*(L+M))/(D +DP); // central differences
	// alpha=0 => central differences,  alpha=1 => full upwind
	d=(1-alpha)*d + alpha*( (KR-KRABS)*R + ( (KR+KRABS)-(KL-KLABS) )*M - (KL+KLABS)*L )/(D+DP); 
      }  
      break;
    case VONOS:
      if(border==FALSE){             // Variable-Order Non-Oscillatory Scheme (VONOS) 2nd/3rd-Order for non-uniform grids	
        ksiR=(KR>=0)?(DP/(2*D)):(DP/(2*PP)); 
	rR=(KR>=0)?(DP/D):(DP/PP);
        ksiL=(KL>=0)?(D/(2*DD)):(D/(2*DP));  
	rL=(KL>=0)?(D/DD):(D/DP);
	q1R=(2+rR)/4, q2R=(2+rR)/(4*(1+rR));
	q1L=(2+rL)/4, q2L=(2+rL)/(4*(1+rL));
	Phi_R=(1-UR)*(M-L)/(((R-L)!=0.0)?(R-L):1E-15) + UR*(R-RR)/(((M-RR)!=0.0)?(M-RR):1E-15);
	Phi_L=(1-UL)*(L-LL)/(((M-LL)!=0.0)?(M-LL):1E-15) + UL*(M-R)/(((L-R)!=0.0)?(L-R):1E-15);

	if((Phi_R>=0) && (Phi_R<(q2R/(10-q1R))))
	  UUR=(1-UR)*(10*M-9*L)+ UR*(10*R-9*RR);
	else if((Phi_R>=(q2R/(10-q1R))) && (Phi_R<(q2R/(1+ksiR-q1R))))
	  UUR=(1-UR)*(q1R*(M-L)+q2R*(R-L)+L)+ UR*(q1R*(R-RR)+q2R*(M-RR)+RR);
	else if((Phi_R>=(q2R/(1+ksiR-q1R))) && (Phi_R<(1./(1+ksiR))))
	  UUR=(1-UR)*((1+ksiR)*(M-L)+L)+ UR*((1+ksiR)*(R-RR)+RR);
	else if((Phi_R>=(1./(1+ksiR))) && (Phi_R<=1))
	  UUR=(1-UR)*R+ UR*M;
	else 
	  UUR=(1-UR)*M+ UR*R;
	
	if((Phi_L>=0) && (Phi_L<(q2L/(10-q1L))))
	  UUL=(1-UL)*(10*L-9*LL)+UL*(10*M-9*R);
	else if((Phi_L>=(q2L/(10-q1L))) && (Phi_L<(q2L/(1+ksiL-q1L))))
	  UUL=(1-UL)*(q1L*(L-LL)+q2L*(M-LL)+LL) +UL*(q1L*(M-R)+q2L*(L-R)+R);
	else if((Phi_L>=(q2L/(1+ksiL-q1L))) && (Phi_L<(1./(1+ksiL))))
	  UUL=(1-UL)*((1+ksiL)*(L-LL)+LL)+ UL*((1+ksiL)*(M-R)+R);
	else if((Phi_L>=(1./(1+ksiL))) && (Phi_L<=1))
	  UUL=(1-UL)*M+UL*L;
	else 
	  UUL=(1-UL)*L+UL*M;
	
	d=(KR*UUR-KL*UUL)/(0.5*(D +DP));
      }else{      //   Hybrid scheme 1st/2nd order
	KRABS=fabs(KR);         
	KLABS=fabs(KL);
	d=(KR*(R+M) -KL*(L+M))/(D +DP); // central differences
	// alpha=0 => central differences,  alpha=1 => full upwind
	d=(1-alpha)*d+alpha*( (KR-KRABS)*R + ( (KR+KRABS)-(KL-KLABS) )*M - (KL+KLABS)*L )/(D+DP); 
	/*
	  NS_REAL UUR, UUL;
	  UUR= alpha * ((1-UR)*M + UR*R) + (1-alpha) * (M+R)/2; // blending upwind 1st order  or  central differences
	  UUL= alpha * ((1-UL)*L + UL*M) + (1-alpha) * (L+M)/2;
	  d=(KR*UUR-KL*UUL)/(0.5*(D +DP));
	*/
      }
      break;
    default:	
      break;
    }

      
    return (d) ;
}


//#################################################################################
//#                                                                               #
//#             e.g. discretization of \partial(uv)/\partial x                    #                 
//#                                                                               #
//#  ------- ------------- -------------- -------------- -------------- -------   #
//#  .      |      .      |      .       |      .       |      .       |      .   #
//#  .      |      .      |      .       |      .       |      .       |      .   #
//#  .      |      .      |      .       |      .       |      .       |      .   #
//#  .      LL     .      L      KL      M      KR      R      .       RR     .   #
//#  .      |      .      |      .       |      .       |      .       |      .   #
//#  .      |      .      |      .       |      .       |      .       |      .   #
//#  .      |      .      |      .       |      .       |      .       |      .   #
//#  --------------------- -------------- -------------- -------------- -------   #
//#   \_____  ____/ \____  _____/ \_____  _____/ \______  ____/ \_____  _____/    #
//#         \/           \/             \/              \/            \/          #
//#         DD           DM             D               DP            PP          #
//#                                                                               #
//#################################################################################

//##########################################################################
//## The discretization for 'mixed' convective term is somewhat different ##
//##########################################################################

NS_REAL NavierCalc::DUV(NS_REAL LL, NS_REAL L, NS_REAL M, NS_REAL R, NS_REAL RR, NS_REAL KL, NS_REAL KR, NS_REAL DD, NS_REAL DM, NS_REAL D, NS_REAL DP, NS_REAL PP, NS_REAL alpha, int border) {
    NS_REAL d, KRABS,KLABS, Phi_L, Phi_R, UUR, UUL, ksiR, ksiL, rL, rR, q1R, q1L, q2R, q2L;
    NS_REAL UR=(KR>=0) ? 0:1;
    NS_REAL UL=(KL>=0) ? 0:1;

    switch(S.Convective){
    case DonorCell:
      KRABS=fabs(KR);          // Donor-Cell
      KLABS=fabs(KL);
      UR=(D*R +DP*M)/(D +DP);// weighted values at cross points of stagg. gr.
      UL=(DM*M +D*L)/(DM +D);
      d=(KR*UR -KL*UL)/D; 
      // alpha=0 => central differences,  alpha=1 => full upwind
      d=(1-alpha)*d+alpha*( (KR-KRABS)*R + ( (KR+KRABS)-(KL-KLABS) )*M - (KL+KLABS)*L )/(2*D); 
      break;
    case QUICK:   
      if(border==FALSE){     // Quadratic upwind interpolation for convective kinematics (QUICK)       
	UR=(1-UR)*((3.*R+6.*M-L)/8.) + UR*((3.*M+6.*R-RR)/8.);
	UL=(1-UL)*((3.*M+6.*L-LL)/8.) + UL*((3.*L+6.*M-R)/8.);
	d=(KR*UR-KL*UL)/D;
      }else{// Donor-Cell full upwind
	KRABS=fabs(KR);       
	KLABS=fabs(KL);
	UR=(D*R +DP*M)/(D +DP);// weighted values at cross points of stagg. gr.
	UL=(DM*M +D*L)/(DM +D);
	d=(KR*UR -KL*UL)/D; 
	// alpha=0 => central differences,  alpha=1 => full upwind
	d=(1-alpha)*d+alpha*( (KR-KRABS)*R + ( (KR+KRABS)-(KL-KLABS) )*M - (KL+KLABS)*L )/(2*D); 
      }  
      break;
    case HLPA:
      if(border==FALSE){              // Hybrid-Linear Parabolic Appoximation (HLPA) 2nd-Order
	Phi_R=(1-UR)*(M-L)/(((R-L)!=0.0)?(R-L):1E-15)  + UR*(R-RR)/(((M-RR)!=0.0)?(M-RR):1E-15);
	Phi_L=(1-UL)*(L-LL)/(((M-LL)!=0.0)?(M-LL):1E-15) + UL*(M-R)/(((L-R)!=0.0)?(L-R):1E-15);
	
	if((Phi_R>1) || (Phi_R<0)) 
	  UR=(1-UR)*M+UR*R;
	else
	  UR=(1-UR)*(M+(R-M)*Phi_R) + UR*(R+(M-R)*Phi_R);
	
	if((Phi_L>1) || (Phi_L<0)) 
	  UL=(1-UL)*L+UL*M;
	else
	  UL=(1-UL)*(L+(M-L)*Phi_L) + UL*(M+(L-M)*Phi_L);
	
	d=(KR*UR-KL*UL)/D ;
      }else{// Donor-Cell full upwind
	KRABS=fabs(KR);       
	KLABS=fabs(KL);
	UR=(D*R +DP*M)/(D +DP);// weighted values at cross points of stagg. gr.
	UL=(DM*M +D*L)/(DM +D);
	d=(KR*UR -KL*UL)/D; 
	// alpha=0 => central differences,  alpha=1 => full upwind
	d=(1-alpha)*d+alpha*( (KR-KRABS)*R + ( (KR+KRABS)-(KL-KLABS) )*M - (KL+KLABS)*L )/(2*D); 
      }
      break;
    case SMART:
      if(border==FALSE){         // Sharp and Monotonic Algorithm for Realistic Transport (SMART) 2nd-Order
	Phi_R=(1-UR)*(M-L)/(((R-L)!=0.0)?(R-L):1E-15)  + UR*(R-RR)/(((M-RR)!=0.0)?(M-RR):1E-15);
        Phi_L=(1-UL)*(L-LL)/(((M-LL)!=0.0)?(M-LL):1E-15) + UL*(M-R)/(((L-R)!=0.0)?(L-R):1E-15);
	
	if((Phi_R>=0) && (Phi_R<(3./74)))
	  UR=(1-UR)*(10*M-9*L)+UR*(10*R-9*RR);
	else if((Phi_R>=(3./74)) && (Phi_R<(5./6)))
	    UR=(1-UR)*((3.*R+6.*M-L)/8.) + UR*((3.*M+6.*R-RR)/8.);
	    //UR=(1-UR)*((3./8)*R+(6./8)*M-(1./8)*L)+UR*((3./8)*M+(6./8)*R-(1./8)*RR);
	else if((Phi_R>=(5./6)) && (Phi_R<=1))
	  UR=(1-UR)*R+UR*M;
	else 
	  UR=(1-UR)*M+UR*R;
	
	if((Phi_L>=0) && (Phi_L<(3./74)))
	  UL=(1-UL)*(10*L-9*LL)+UL*(10*M-9*R);
	else if((Phi_L>=(3./74)) && (Phi_L<(5./6)))
	  UL=(1-UL)*((3.*M+6.*L-LL)/8.) + UL*((3.*L+6.*M-R)/8.);
	  // UL=(1-UL)*((3./8)*M+(6./8)*L-(1./8)*LL)+UL*((3./8)*L+(6./8)*M-(1./8)*R);
	else if((Phi_L>=(5./6)) && (Phi_L<=1))
	  UL=(1-UL)*M+UL*L;
	else 
	  UL=(1-UL)*L+UL*M;
	
	d=(KR*UR-KL*UL)/D;
      }else{// Donor-Cell full upwind
	KRABS=fabs(KR);       
	KLABS=fabs(KL);
	UR=(D*R +DP*M)/(D +DP);// weighted values at cross points of stagg. gr.
	UL=(DM*M +D*L)/(DM +D);
	d=(KR*UR -KL*UL)/D; 
	// alpha=0 => central differences,  alpha=1 => full upwind
	d=(1-alpha)*d+alpha*( (KR-KRABS)*R + ( (KR+KRABS)-(KL-KLABS) )*M - (KL+KLABS)*L )/(2*D); 
      }
      break;
    case VONOS:     
      if(border==FALSE){             // Variable-Order Non-Oscillatory Scheme (VONOS) 2nd/3rd-Order for non-uniform grids
	ksiR = (KR>=0)?((D+DP)/(2*(D+DM))):((D+DP)/(2*(PP+DP)));  
	rR   = (KR>=0)?((D+DP)/(D+DM)):((DP+D)/(DP+PP));
	ksiL = (KL>=0)?((DM+D)/(2*(DM+DD))):((D+DM)/(2*(D+DP)));
	rL   = (KL>=0)?((D+DM)/(DD+DM)):((D+DM)/(D+DP));
        q1R  =  (2+rR)/4,  q2R=(2+rR)/(4*(1+rR));
        q1L  =  (2+rL)/4,  q2L=(2+rL)/(4*(1+rL));
        Phi_R=(1-UR)*(M-L)/(((R-L)!=0.0)?(R-L):1E-15) + UR*(R-RR)/(((M-RR)!=0.0)?(M-RR):1E-15);
        Phi_L=(1-UL)*(L-LL)/(((M-LL)!=0.0)?(M-LL):1E-15) + UL*(M-R)/(((L-R)!=0.0)?(L-R):1E-15);
	
	if((Phi_R>=0) && (Phi_R<(q2R/(10-q1R))))
	  UUR=(1-UR)*(10*M-9*L)+ UR*(10*R-9*RR);
	else if((Phi_R>=(q2R/(10-q1R))) && (Phi_R<(q2R/(1+ksiR-q1R))))
	  UUR=(1-UR)*(q1R*(M-L)+q2R*(R-L)+L)+ UR*(q1R*(R-RR)+q2R*(M-RR)+RR);
	else if((Phi_R>=(q2R/(1+ksiR-q1R))) && (Phi_R<(1./(1+ksiR))))
	  UUR=(1-UR)*((1+ksiR)*(M-L)+L)+ UR*((1+ksiR)*(R-RR)+RR);
	else if((Phi_R>=(1./(1+ksiR))) && (Phi_R<=1))
	  UUR=(1-UR)*R+ UR*M;
	else 
	  UUR=(1-UR)*M+ UR*R;
	
	if((Phi_L>=0) && (Phi_L<(q2L/(10-q1L))))
	  UUL=(1-UL)*(10*L-9*LL)+UL*(10*M-9*R);
	else if((Phi_L>=(q2L/(10-q1L))) && (Phi_L<(q2L/(1+ksiL-q1L))))
	  UUL=(1-UL)*(q1L*(L-LL)+q2L*(M-LL)+LL) +UL*(q1L*(M-R)+q2L*(L-R)+R);
	else if((Phi_L>=(q2L/(1+ksiL-q1L))) && (Phi_L<(1./(1+ksiL))))
	  UUL=(1-UL)*((1+ksiL)*(L-LL)+LL)+ UL*((1+ksiL)*(M-R)+R);
	else if((Phi_L>=(1./(1+ksiL))) && (Phi_L<=1))
	  UUL=(1-UL)*M+UL*L;
	else 
	  UUL=(1-UL)*L+UL*M;
	
	d=(KR*UUR-KL*UUL)/D;
      }else{      // Hybrid-scheme 1st/2nd order
	KRABS=fabs(KR);       
	KLABS=fabs(KL);
	UR=(D*R +DP*M)/(D +DP);// weighted values at cross points of stagg. gr.
	UL=(DM*M +D*L)/(DM +D);
	d=(KR*UR -KL*UL)/D; 
	// alpha=0 => central differences,  alpha=1 => full upwind
	d=(1-alpha)*d+alpha*( (KR-KRABS)*R + ( (KR+KRABS)-(KL-KLABS) )*M - (KL+KLABS)*L )/(2*D); 
	/* 
	   NS_REAL UUR, UUL;
	   UUR= alpha * ((1-UR)*M + UR*R) + (1-alpha) * (D*R +DP*M)/(D +DP); //  blending upwind 1st order  or  central differences
	   UUL= alpha * ((1-UL)*L + UL*M) + (1-alpha) * (DM*M +D*L)/(DM +D);
	   d=(KR*UUR-KL*UUL)/D;
	*/
      }
      break;
    default: 
      break;
    }
    return (d) ;
}
