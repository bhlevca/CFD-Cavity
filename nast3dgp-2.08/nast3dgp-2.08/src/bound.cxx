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
#include <errno.h>
#include <stdlib.h>

#ifdef TIMED
#define STARTTIMER(ii) int oldtimer=timer.Start(ii,Par);
#define STOPTIMER timer.Stop(oldtimer,Par);
#else
#define STARTTIMER(ii)
#define STOPTIMER
#endif


void NavierCalc::GENENTRY(int i,int j,int k,NS_List<CELLREC>& list) {
    NS_List<CELLREC>* p2=new NS_List<CELLREC>;
    p2->data.SetCoord(i,j,k);
    list.InsertFirst(p2);
}

void NavierCalc::GENSPECENTRY(int i,int j,int k, NS_List<SPECCELLREC>& list) {
    NS_List<SPECCELLREC>* p3=new NS_List<SPECCELLREC>;
    p3->data.SetCoord(i,j,k);      
    p3->data.gweight=1;                                                                 
    if((i>=Par->iug) && (flag[i][j][k]&SOUTH )) p3->data.gweight*=DX[i];                  
    if((i<=Par->iog) && (flag[i][j][k]&NORTH )) p3->data.gweight*=DX[i];                  
    if((j>=Par->jug) && (flag[i][j][k]&EAST  )) p3->data.gweight*=DY[j];                  
    if((j<=Par->jog) && (flag[i][j][k]&WEST  )) p3->data.gweight*=DY[j];                  
    if((k>=Par->kug) && (flag[i][j][k]&BOTTOM)) p3->data.gweight*=DZ[k];                  
    if((k<=Par->kog) && (flag[i][j][k]&TOP   )) p3->data.gweight*=DZ[k];                  
    p3->data.den=0;                                                                     
    if((i>=Par->iug) && (flag[i][j][k]&SOUTH )) p3->data.den+=p3->data.gweight/DX[i];     
    if((i<=Par->iog) && (flag[i][j][k]&NORTH )) p3->data.den+=p3->data.gweight/DX[i];     
    if((j>=Par->jug) && (flag[i][j][k]&EAST  )) p3->data.den+=p3->data.gweight/DY[j];     
    if((j<=Par->jog) && (flag[i][j][k]&WEST  )) p3->data.den+=p3->data.gweight/DY[j];     
    if((k>=Par->kug) && (flag[i][j][k]&BOTTOM)) p3->data.den+=p3->data.gweight/DZ[k];     
    if((k<=Par->kog) && (flag[i][j][k]&TOP   )) p3->data.den+=p3->data.gweight/DZ[k];     
    list.InsertFirst(p3);
}

void NavierCalc::GENPENTRY(int i,int j,int k,int cg,int listnumber) {
    if(flag[i][j][k]&INOUT) {
        switch(INOUTTYP(flag[i][j][k])) {
        case 0: // simple neumann zero bound. cond.    
	    GENENTRY(i,j,k,PInOut1List[cg][listnumber]);break;
        case 1: // convective time dep. bound. cond.    
	    GENENTRY(i,j,k,PInOut2List[cg][listnumber]);break;
        case 2:// this bound. cond. are under construction
	    if(((listnumber==0) || (listnumber==1)) && 
	       !( (j>=Par->jug)&&(j<=Par->jog)&&(k>=Par->kug)&&(k<=Par->kog) )) break; 
	    if(((listnumber==2) || (listnumber==3)) && 
	       !((i>=Par->iug)&&(i<=Par->iog)&& (k>=Par->kug)&&(k<=Par->kog))) break;
	    if(((listnumber==4) || (listnumber==5)) && 
	       !((j>=Par->jug)&&(j<=Par->jog)&& (i>=Par->iug)&&(i<=Par->iog))) break;
	    GENENTRY(i,j,k,PInOut2List[cg][listnumber]);break; 
        default:    
	    GENENTRY(i,j,k,PdpdnList[cg][listnumber]);break;
        }
    } else GENENTRY(i,j,k,PdpdnList[cg][listnumber]);
}


#define TEMPINFLOWFLAG 0x80000000

//
// all OBSTACLE and BOUNDARY (pressure) cells are listed in various lists 
// each of these lists  collects cells which share some particular property
// e.g. this can be all cells where inflow boundary conditions are prescribed.
// other examples are 
// 
//   InflowList                           collect all cells, where Inflow data is prescribed
// 
//   SlipList                             cells, where Slip-BC are prescribed
//
//   NoSlipList                           cells, where NoSlip-BC are prescribed
//
//   InOutList                            cells, where any outflow-boundary cond. is assumed
//
//   The colored SOR iteration for the pressure poisson equation suggests  to
//   split the lists above further into cells with same color 0..7
//
// 
//   PdpdnList[color][face]               cells where prescribed boundary conditions (BC) lead to homogeneous
//                                        Neumann BC for the pressure and which have exactly 
//                                        one fluid cell neighbouring at face 'face'
//                                        The cells have color 'color' (color in 0...7)
//   
//   PSpecdpdnList[color]                 cells where prescribed boundary conditions (BC) lead to homogeneous
//                                        Neumann BC for the pressure but have more than one 
//                                        neighbouring fluid cell. This requires some more
//                                        involved boundary treatment
//
//   PInOut2List[color]                   cells with the prescribed outflow condition 2, 
//                                        it is assumed that no cells with more than one neighbouring fluid 
//                                        cell exist
//                                        Some special boundary treatment is required.
//   
//
//  extract local lists (InflowList,..) or each processor of global lists (S.InlowList,..)
//  the global lists are traversed if one element is responsible for a cell whithin the 
//  local (process dependent) computational domain this element is insert in the local 
//  counterpart of the global list.
//

void NavierCalc::BuildLocalLists() {                             // create local lists of boundary cells
    NS_List<INFLOW>   *p ,*q  ;
    NS_List<INOUTREC> *p1,*q1 ;
    int i,j,k,n;

    p=S.InflowList.GetNext() ;                
    while(p) {
        i=p->data.i , j=p->data.j , k=p->data.k ;
        if(PPINALLBORDER(i,j,k)) IFBORDER(flag[i][j][k]) {
            flag[i][j][k] |= TEMPINFLOWFLAG  ;       // mark INFLOW-Cells in the flag field to have 
            q              = new NS_List<INFLOW>;    // later easy access to the information whether a cell
            q->data        = p->data         ;       // has inflow boundary conditions or not
            InflowList.InsertFirst(q)        ;
        }
        p=p->GetNext();
    }
   
    p1=S.InOutList.GetNext();                        
    while(p1) {
        i=p1->data.i, j=p1->data.j, k=p1->data.k;
        if(PPINALLBORDER(i,j,k)) IFBORDER(flag[i][j][k]) {
            if(p1->data.type!=3) {                    // collect all cells with inout 1 or 2 b.c. 
                q1       = new NS_List<INOUTREC> ;
                q1->data = p1->data           ;       // contains type of b.c. and addit. parameters
                InOutList.InsertFirst(q1)     ;
            } else {                                  // collect all cells with inout 3 b.c.
		printf(" Wrong Boundary Conditions (only b.c. type 1 and b.c. type 2) ! \n"); exit(1);
// --- UNDER CONSTRUCTION ---
//                 NS_List<INOUTRECDATA>* q2=new NS_List<INOUTRECDATA>;
//                 q2->data=p1->data;
//                 InOut3List.InsertFirst(q2);
// --------------------------
            }
        }
        p1=p1->GetNext();
    }

    //
    // Set local lists, depending only on the value of the flag field
    //
    NS_List<CELLREC> *p2;
    PIJKALLLOOP IFBORDER(flag[i][j][k]) {
      if ((flag[i][j][k]&(INOUT|TEMPINFLOWFLAG))==0) {   // no INOUT/INFLOW cell -> Slip or NoSlip cell
        p2=new NS_List<CELLREC>;
        p2->data.SetCoord(i,j,k);
        if (flag[i][j][k]&SLIP) 
          SlipList.InsertFirst(p2)  ; 
        else 
          NoslipList.InsertFirst(p2);
      }

      // All boundary cells with TEMP/CHEMB have homog.
      // Neumann-Conditions for T/CH these cells are collected
      // here 
      if(S.CompTemp) if(flag[i][j][k]&TEMPB) GENSPECENTRY(i,j,k,dTempdnList); 
      if(S.CompChem) if(flag[i][j][k]&CHEMB) GENSPECENTRY(i,j,k,dChemdnList);
    }

 
    PIJKALLLOOP {                                        // generate list for p boundaries
         n=(4*(i&1)+2*(j&1)+1*(k&1))^7;                  // color group
                                                         // Except the cells with the second outflow b.c. 
                                                         // all pressure cells have dp/dn=0

        if(i<=Par->iog) IFBORDER(flag[i+1][j][k]) if(FREECELLS(flag[i+1][j][k])==SOUTH)  GENPENTRY(i+1,j,k,n,0);
        if(i>=Par->iug) IFBORDER(flag[i-1][j][k]) if(FREECELLS(flag[i-1][j][k])==NORTH)  GENPENTRY(i-1,j,k,n,1);
        if(j<=Par->jog) IFBORDER(flag[i][j+1][k]) if(FREECELLS(flag[i][j+1][k])==EAST)   GENPENTRY(i,j+1,k,n,2); 
        if(j>=Par->jug) IFBORDER(flag[i][j-1][k]) if(FREECELLS(flag[i][j-1][k])==WEST)   GENPENTRY(i,j-1,k,n,3);
        if(k<=Par->kog) IFBORDER(flag[i][j][k+1]) if(FREECELLS(flag[i][j][k+1])==BOTTOM) GENPENTRY(i,j,k+1,n,4);
        if(k>=Par->kug) IFBORDER(flag[i][j][k-1]) if(FREECELLS(flag[i][j][k-1])==TOP)    GENPENTRY(i,j,k-1,n,5);

        // treat boundary cells with more than one 
        // neighbourting fluid cell:
        // !!! at such cells NO inout BC is 
        // allowed !!!
        IFBORDER(flag[i][j][k]) if((FREECELLS(flag[i][j][k])!=NORTH)  && (FREECELLS(flag[i][j][k])!=SOUTH) &&
                                   (FREECELLS(flag[i][j][k])!=EAST)   && (FREECELLS(flag[i][j][k])!=WEST)  &&
                                   (FREECELLS(flag[i][j][k])!=BOTTOM) && (FREECELLS(flag[i][j][k])!=TOP))
          GENSPECENTRY(i,j,k,PSpecdpdnList[n^7]);
    }
    PIJKALLLOOP flag[i][j][k]&=~TEMPINFLOWFLAG;                 // clear INFLOW flag  

    //
    // Calculate inflow-surface 
    //
    inflow_surface=0.0;
    flow_in=0.0;
    unsigned flg;   
    for(p=InflowList.GetNext();p;p=p->GetNext()) {
      i=p->data.i; j=p->data.j; k=p->data.k ;
      flg=flag[i][j][k];       
      if((flg&SOUTH) && (i >Par->iug) && (j>=Par->jug) && (j<=Par->jog) && (k>=Par->kug) && (k<=Par->kog)) 
	{inflow_surface+=DZ[k]*DY[j]; flow_in+=DZ[k]*DY[j]*p->data.u;}   
      if((flg&NORTH) && (i<=Par->iog) && (j>=Par->jug) && (j<=Par->jog) && (k>=Par->kug) && (k<=Par->kog)) 
	{inflow_surface+=DZ[k]*DY[j]; flow_in+=DZ[k]*DY[j]*p->data.u;}
      
      if((flg&EAST)  && (j >Par->jug) && (i>=Par->iug) && (i<=Par->iog) && (k>=Par->kug) && (k<=Par->kog)) 
	{inflow_surface+=DZ[k]*DX[i]; flow_in+=DZ[k]*DX[i]*p->data.v;}  
      if((flg&WEST)  && (j<=Par->jog) && (i>=Par->iug) && (i<=Par->iog) && (k>=Par->kug) && (k<=Par->kog)) 
	{inflow_surface+=DZ[k]*DX[i]; flow_in+=DZ[k]*DX[i]*p->data.v;}
      
      if((flg&BOTTOM)&& (k >Par->kug) && (i>=Par->iug) && (i<=Par->iog) && (j>=Par->jug) && (j<=Par->jog)) 
	{inflow_surface+=DY[j]*DX[i]; flow_in+=DY[j]*DX[i]*p->data.w;}  
      if((flg&TOP)   && (k<=Par->kog) && (i>=Par->iug) && (i<=Par->iog) && (j>=Par->jug) && (j<=Par->jog)) 
	{inflow_surface+=DY[j]*DX[i]; flow_in+=DY[j]*DX[i]*p->data.w;}	
    }
    inflow_surface=Par->CommSum(inflow_surface);
    flow_in=Par->CommSum(flow_in),

    outflow_surface=0.0;
    for(p1=InOutList.GetNext();p1;p1=p1->GetNext()) {
	i=p1->data.i; j=p1->data.j; k=p1->data.k;
	flg=flag[i][j][k];      

      if((flg&SOUTH) && (i >Par->iug) && (j>=Par->jug) && (j<=Par->jog) && (k>=Par->kug) && (k<=Par->kog)) 
	outflow_surface+=DZ[k]*DY[j];   
      if((flg&NORTH) && (i<=Par->iog) && (j>=Par->jug) && (j<=Par->jog) && (k>=Par->kug) && (k<=Par->kog)) 
	outflow_surface+=DZ[k]*DY[j];
      
      if((flg&EAST)  && (j >Par->jug) && (i>=Par->iug) && (i<=Par->iog) && (k>=Par->kug) && (k<=Par->kog)) 
	outflow_surface+=DZ[k]*DX[i];   
      if((flg&WEST)  && (j<=Par->jog) && (i>=Par->iug) && (i<=Par->iog) && (k>=Par->kug) && (k<=Par->kog)) 
	outflow_surface+=DZ[k]*DX[i];
      
      if((flg&BOTTOM)&& (k >Par->kug) && (i>=Par->iug) && (i<=Par->iog) && (j>=Par->jug) && (j<=Par->jog)) 
	outflow_surface+=DY[j]*DX[i];   
      if((flg&TOP)   && (k<=Par->kog) && (i>=Par->iug) && (i<=Par->iog) && (j>=Par->jug) && (j<=Par->jog)) 
	outflow_surface+=DY[j]*DX[i];	
    }
    outflow_surface=Par->CommSum(outflow_surface);
    if(Par->IsParent()) if(printinfo>=2) 
      printf("inflow_surface = %g  outflow_surface = %g  flow_in = %g\n",inflow_surface, outflow_surface, flow_in); 
}

//
//   Set pressure boundary conditions whithin the SOR iteration
//

// Difference Stars required for implem. of second outflow b.c.
#define PXX2(i,j,k) (S.ddstar[0][2][i]*P[i+1][j][k]  +  S.ddstar[0][0][i]*P[i-1][j][k])
#define PYY2(i,j,k) (S.ddstar[1][2][j]*P[i][j+1][k]  +  S.ddstar[1][0][j]*P[i][j-1][k])
#define PZZ2(i,j,k) (S.ddstar[2][2][k]*P[i][j][k+1]  +  S.ddstar[2][0][k]*P[i][j][k-1])

void NavierCalc::SetPBorder(Matrix<NS_REAL>& B,int cg) {
    int oldtimer=timer->Start(7,Par);
    NS_List<CELLREC>* p;
    NS_List<SPECCELLREC>* p1;

    for(p=PdpdnList[cg][0].GetNext();p;p=p->GetNext()) B[p->data.i][p->data.j][p->data.k] =B[p->data.i-1][p->data.j][p->data.k]; // NORD
    for(p=PdpdnList[cg][1].GetNext();p;p=p->GetNext()) B[p->data.i][p->data.j][p->data.k] =B[p->data.i+1][p->data.j][p->data.k]; // SOUTH
    for(p=PdpdnList[cg][2].GetNext();p;p=p->GetNext()) B[p->data.i][p->data.j][p->data.k] =B[p->data.i][p->data.j-1][p->data.k]; // WEST
    for(p=PdpdnList[cg][3].GetNext();p;p=p->GetNext()) B[p->data.i][p->data.j][p->data.k] =B[p->data.i][p->data.j+1][p->data.k]; // EAST
    for(p=PdpdnList[cg][4].GetNext();p;p=p->GetNext()) B[p->data.i][p->data.j][p->data.k] =B[p->data.i][p->data.j][p->data.k-1]; // TOP
    for(p=PdpdnList[cg][5].GetNext();p;p=p->GetNext()) B[p->data.i][p->data.j][p->data.k] =B[p->data.i][p->data.j][p->data.k+1]; // BOTTOM
  
    for(p1=PSpecdpdnList[cg].GetNext();p1;p1=p1->GetNext()) {   
	// for cells with more than one neighbouring fluid cell
	// we use weighted copies of the neighb. pressure cells
	// this agrees with dp/dn=0 where the normal direction  
	// is a weighted from the two or more normals to the fluid cells
	int i=p1->data.i, j=p1->data.j, k=p1->data.k;
	B[i][j][k]=0;
	unsigned flg=flag[i][j][k];
	if((flg&SOUTH ) && (i>=Par->iug)) {
	    B[i][j][k]+=(B[i-1][j][k]*(p1->data.gweight/DX[i]));
	}
	if((flg&NORTH ) && (i<=Par->iog)) {
	    B[i][j][k]+=(B[i+1][j][k]*(p1->data.gweight/DX[i]));
	}
	if((flg&EAST  ) && (j>=Par->jug)) {
	    B[i][j][k]+=(B[i][j-1][k]*(p1->data.gweight/DY[j]));
	}
	if((flg&WEST  ) && (j<=Par->jog)) {
	    B[i][j][k]+=(B[i][j+1][k]*(p1->data.gweight/DY[j]));
	}
	if((flg&BOTTOM) && (k>=Par->kug)) {
	    B[i][j][k]+=(B[i][j][k-1]*(p1->data.gweight/DZ[k]));
	}
	if((flg&TOP   ) && (k<=Par->kog)) {
	    B[i][j][k]+=(B[i][j][k+1]*(p1->data.gweight/DZ[k]));
	}
	B[i][j][k]/=p1->data.den; 
    }
        
    //
    // Treatment of InOutBC1 for pressure
    //
    for(p=PInOut1List[cg][0].GetNext();p;p=p->GetNext()) // NORD
	B[p->data.i][p->data.j][p->data.k] = B[p->data.i-1][p->data.j][p->data.k];
    for(p=PInOut1List[cg][1].GetNext();p;p=p->GetNext()) // SOUTH
	B[p->data.i][p->data.j][p->data.k] = B[p->data.i+1][p->data.j][p->data.k];
    for(p=PInOut1List[cg][2].GetNext();p;p=p->GetNext()) // WEST
	B[p->data.i][p->data.j][p->data.k] = B[p->data.i][p->data.j-1][p->data.k];
    for(p=PInOut1List[cg][3].GetNext();p;p=p->GetNext()) // EAST
	B[p->data.i][p->data.j][p->data.k] = B[p->data.i][p->data.j+1][p->data.k];
    for(p=PInOut1List[cg][4].GetNext();p;p=p->GetNext()) // TOP
	B[p->data.i][p->data.j][p->data.k] = B[p->data.i][p->data.j][p->data.k-1];
    for(p=PInOut1List[cg][5].GetNext();p;p=p->GetNext()) // BOTTOM
	B[p->data.i][p->data.j][p->data.k] = B[p->data.i][p->data.j][p->data.k+1];

    //
    // Treatment of InOutBC2 for pressure 
    //
    for(p=PInOut2List[cg][0].GetNext();p;p=p->GetNext()) // NORD
	B[p->data.i][p->data.j][p->data.k] = B[p->data.i-1][p->data.j][p->data.k];
    for(p=PInOut2List[cg][1].GetNext();p;p=p->GetNext()) // SOUTH
	B[p->data.i][p->data.j][p->data.k] = B[p->data.i+1][p->data.j][p->data.k];
    for(p=PInOut2List[cg][2].GetNext();p;p=p->GetNext()) // WEST
	B[p->data.i][p->data.j][p->data.k] = B[p->data.i][p->data.j-1][p->data.k];
    for(p=PInOut2List[cg][3].GetNext();p;p=p->GetNext()) // EAST
	B[p->data.i][p->data.j][p->data.k] = B[p->data.i][p->data.j+1][p->data.k];
    for(p=PInOut2List[cg][4].GetNext();p;p=p->GetNext()) // TOP
	B[p->data.i][p->data.j][p->data.k] = B[p->data.i][p->data.j][p->data.k-1];
    for(p=PInOut2List[cg][5].GetNext();p;p=p->GetNext()) // BOTTOM
	B[p->data.i][p->data.j][p->data.k] = B[p->data.i][p->data.j][p->data.k+1];

    timer->Stop(oldtimer,Par);
}

void NavierCalc::SetPALLBorder(Matrix <NS_REAL> &X){
    SetPBorder(X,0|0|0);
    SetPBorder(X,4|2|0);
    SetPBorder(X,4|0|1);
    SetPBorder(X,0|2|1);
    SetPBorder(X,4|0|0);
    SetPBorder(X,0|2|0);
    SetPBorder(X,0|0|1);
    SetPBorder(X,4|2|1);
}


void NavierCalc::SetPLocalBorder(Matrix <NS_REAL> &X, int i, int j, int k){
    IFOBSTACLE(flag[i-1][j][k]) X[i-1][j][k]=X[i][j][k];
    IFOBSTACLE(flag[i+1][j][k]) X[i+1][j][k]=X[i][j][k];
    IFOBSTACLE(flag[i][j-1][k]) X[i][j-1][k]=X[i][j][k];
    IFOBSTACLE(flag[i][j+1][k]) X[i][j+1][k]=X[i][j][k];
    IFOBSTACLE(flag[i][j][k-1]) X[i][j][k-1]=X[i][j][k];
    IFOBSTACLE(flag[i][j][k+1]) X[i][j][k+1]=X[i][j][k];
}



/***********************************************************************
  
 In SetObstacleCondTilde Boundary values for the auxiliary velocity fields are set,
 in general at this stage this is only required for the normal velocity components !

 Temperature and Species are also treated here

 ***********************************************************************/  

inline void NavierCalc::SetHomogeneousNeumannNormal(Matrix<NS_REAL>& U,Matrix<NS_REAL>& V,Matrix<NS_REAL>& W,int i,int j,int k) {
    unsigned flg=flag[i][j][k];
        
    // Sets one sided neumann zero values only for normal velocity components along the boundaey
    if((flg&SOUTH) && (i >Par->iug)) U[i-1][j][k]=U[i-2][j][k];   
    if((flg&NORTH) && (i<=Par->iog)) U[i  ][j][k]=U[i+1][j][k];

    if((flg&EAST)  && (j >Par->jug)) V[i][j-1][k]=V[i][j-2][k];   
    if((flg&WEST)  && (j<=Par->jog)) V[i][j  ][k]=V[i][j+1][k];

    if((flg&BOTTOM)&& (k >Par->kug)) W[i][j][k-1]=W[i][j][k-2];   
    if((flg&TOP)   && (k<=Par->kog)) W[i][j][k  ]=W[i][j][k+1];
} 

inline void NavierCalc::SetDirichletNormal(Matrix<NS_REAL>& U,Matrix<NS_REAL>& V,Matrix<NS_REAL>& W,
					   int i,int j,int k,NS_REAL u,NS_REAL v,NS_REAL w) {
    unsigned flg=flag[i][j][k];

    // Sets dirichlet values only for normal velocity components along the boundaey
    if ((flg&SOUTH) && (i>=Par->iug)) U[i-1][j][k]=u;
    if ((flg&NORTH) && (i<=Par->iog)) U[i][j][k]  =u;
    if ((flg&EAST)  && (j>=Par->jug)) V[i][j-1][k]=v;
    if ((flg&WEST)  && (j<=Par->jog)) V[i][j][k]  =v;
    if ((flg&BOTTOM)&& (k>=Par->kug)) W[i][j][k-1]=w;
    if ((flg&TOP)   && (k<=Par->kog)) W[i][j][k]  =w;                                                                        
} 

inline void NavierCalc::SetVelocityExtrapolation(Matrix<NS_REAL>& U,Matrix<NS_REAL>& V,Matrix<NS_REAL>& W,int i,int j,int k) {
    unsigned flg=flag[i][j][k];
//
//  time-dep. convective outflow boundary cond. 2nd order 
//
    if(S.timestepmethod){
	if((S.TimeDis==2 || S.TimeDis==3) && !actual_fgh) S.delt/=2;  // Only for Runge-Kutta schemes
	if((flg&SOUTH) && (i >Par->iug)) {
	    U[i-1][j][k]  = U[i-1][j][k] - S.delt * U[i-1][j][k] * (3*U[i-1][j][k] - 4*U[i-2][j][k] + U[i-3][j][k])/(DX[i-1]+DX[i-2]);
	    V[i  ][j][k]  = V[i  ][j][k] - S.delt * V[i  ][j][k] * (3*V[i  ][j][k] - 4*V[i-1][j][k] + V[i-2][j][k])/(DXM[i-1]+DXM[i-2]);
	    W[i  ][j][k]  = W[i  ][j][k] - S.delt * W[i  ][j][k] * (3*W[i  ][j][k] - 4*W[i-1][j][k] + W[i-2][j][k])/(DXM[i-1]+DXM[i-2]);
	}
	if((flg&NORTH) && (i<=Par->iog)) {
	    U[i  ][j][k]  = U[i  ][j][k] - S.delt * U[i  ][j][k] * (-U[i+2][j][k] + 4*U[i+1][j][k] - 3*U[i  ][j][k])/(DX[i+1]+DX[i+2]);  
	    V[i  ][j][k]  = V[i  ][j][k] - S.delt * V[i  ][j][k] * (-V[i+2][j][k] + 4*V[i+1][j][k] - 3*V[i  ][j][k])/(DXM[i+1]+DXM[i+2]);
	    W[i  ][j][k]  = W[i  ][j][k] - S.delt * W[i  ][j][k] * (-W[i+2][j][k] + 4*W[i+1][j][k] - 3*W[i  ][j][k])/(DXM[i+1]+DXM[i+2]);
	}
	
	if((flg&EAST)  && (j >Par->jug)) {
	    U[i][j  ][k]  = U[i][j  ][k] - S.delt * U[i][j  ][k] * (3*U[i][j  ][k] - 4*U[i][j-1][k] + U[i][j-2][k])/(DYM[j-1]+DYM[j-2]);
	    V[i][j-1][k]  = V[i][j-1][k] - S.delt * V[i][j-1][k] * (3*V[i][j-1][k] - 4*V[i][j-2][k] + V[i][j-3][k])/(DY[j-1]+DY[j-2]);
	    W[i][j  ][k]  = W[i][j  ][k] - S.delt * W[i][j  ][k] * (3*W[i][j  ][k] - 4*W[i][j-1][k] + W[i][j-2][k])/(DYM[j-1]+DYM[j-2]);
	}
	if((flg&WEST)  && (j<=Par->jog)) {
	    U[i][j  ][k]  = U[i][j  ][k] - S.delt * U[i][j  ][k] * (-U[i][j+2][k] + 4*U[i][j+1][k] - 3*U[i][j  ][k])/(DYM[j+1]+DYM[j+2]);
	    V[i][j  ][k]  = V[i][j  ][k] - S.delt * V[i][j  ][k] * (-V[i][j+2][k] + 4*V[i][j+1][k] - 3*V[i][j  ][k])/(DY[j+1]+DY[j+2]);
	    W[i][j  ][k]  = W[i][j  ][k] - S.delt * W[i][j  ][k] * (-W[i][j+2][k] + 4*W[i][j+1][k] - 3*W[i][j  ][k])/(DYM[j+1]+DYM[j+2]);
	}
	
	if((flg&BOTTOM)&& (k >Par->kug)) {
	    U[i][j][k  ]  = U[i][j][k  ] - S.delt * U[i][j][k  ] * (3*U[i][j][k  ] - 4*U[i][j][k-1] + U[i][j][k-2])/(DZM[k-1]+DZM[k-2]);
	    V[i][j][k  ]  = V[i][j][k  ] - S.delt * V[i][j][k  ] * (3*V[i][j][k  ] - 4*V[i][j][k-1] + V[i][j][k-2])/(DZM[k-1]+DZM[k-2]);    
	    W[i][j][k-1]  = W[i][j][k-1] - S.delt * W[i][j][k-1] * (3*W[i][j][k-1] - 4*W[i][j][k-2] + W[i][j][k-3])/(DZ[k-1]+DZ[k-2]);
	}
	if((flg&TOP)   && (k<=Par->kog)) {
	    U[i][j][k  ]  = U[i][j][k  ] - S.delt * U[i][j][k  ] * (-U[i][j][k+2] + 4*U[i][j][k+1] - 3*U[i][j][k  ])/(DZM[k+1]+DZM[k+2]);
	    V[i][j][k  ]  = V[i][j][k  ] - S.delt * V[i][j][k  ] * (-V[i][j][k+2] + 4*V[i][j][k+1] - 3*V[i][j][k  ])/(DZM[k+1]+DZM[k+2]);
	    W[i][j][k  ]  = W[i][j][k  ] - S.delt * W[i][j][k  ] * (-W[i][j][k+2] + 4*W[i][j][k+1] - 3*W[i][j][k  ])/(DZ[k+1]+DZ[k+2]);
	}
	if((S.TimeDis==2 || S.TimeDis==3) && !actual_fgh) S.delt*=2;  // Only for Runge-Kutta schemes
	if(S.TimeDis==2 || S.TimeDis==3) actual_fgh^=1;              // Only for Runge-Kutta schemes
    }       
}


inline void NavierCalc::SetVelocityCorrection(Matrix<NS_REAL>& U,Matrix<NS_REAL>& V,Matrix<NS_REAL>& W,int i,int j,int k) {
  //
  // Outflow boundary correction for normal velocity, 
  // so that \int_{\partial \Omega} (velocityfield * normalvector) dx = 0 
  //       
  unsigned flg=flag[i][j][k];

  if((flg&SOUTH) && (i >Par->iug)) U[i-1][j][k]+=mass_diff/outflow_surface;   
  if((flg&NORTH) && (i<=Par->iog)) U[i  ][j][k]+=mass_diff/outflow_surface;
  
  if((flg&EAST)  && (j >Par->jug)) V[i][j-1][k]+=mass_diff/outflow_surface;   
  if((flg&WEST)  && (j<=Par->jog)) V[i][j  ][k]+=mass_diff/outflow_surface;
  
  if((flg&BOTTOM)&& (k >Par->kug)) W[i][j][k-1]+=mass_diff/outflow_surface;   
  if((flg&TOP)   && (k<=Par->kog)) W[i][j][k  ]+=mass_diff/outflow_surface;
} 


void NavierCalc::OutflowVelocityCorrection(Matrix<NS_REAL>& U,Matrix<NS_REAL>& V,Matrix<NS_REAL>& W, int method) {
    //
    //  OUTFLOW velocity correction
    //
    int                  i,j,k;
    unsigned             flg;
    NS_List<INOUTREC>*   p1;
    
    if(printinfo>2 || S.iobd){
	flow_out=0.0;
	for(p1=InOutList.GetNext();p1;p1=p1->GetNext()) {
	    if(method!=-1){
		i=p1->data.i; j=p1->data.j; k=p1->data.k;
		flg=flag[i][j][k];      
		if((flg&SOUTH) && (i >Par->iug) && (j>=Par->jug) && (j<=Par->jog) && (k>=Par->kug) && (k<=Par->kog)) 
		    flow_out+=DZ[k]*DY[j]*U[i-1][j][k]; 
		if((flg&NORTH) && (i<=Par->iog) && (j>=Par->jug) && (j<=Par->jog) && (k>=Par->kug) && (k<=Par->kog)) 
		    flow_out+=DZ[k]*DY[j]*U[i][j][k];
		
		if((flg&EAST)  && (j >Par->jug) && (i>=Par->iug) && (i<=Par->iog) && (k>=Par->kug) && (k<=Par->kog)) 
		    flow_out+=DZ[k]*DX[i]*V[i][j-1][k];   
		if((flg&WEST)  && (j<=Par->jog) && (i>=Par->iug) && (i<=Par->iog) && (k>=Par->kug) && (k<=Par->kog)) 
		    flow_out+=DZ[k]*DX[i]*V[i][j][k];
		
		if((flg&BOTTOM)&& (k >Par->kug) && (i>=Par->iug) && (i<=Par->iog) && (j>=Par->jug) && (j<=Par->jog)) 
		    flow_out+=DY[j]*DX[i]*W[i][j][k-1];   
		if((flg&TOP)   && (k<=Par->kog) && (i>=Par->iug) && (i<=Par->iog) && (j>=Par->jug) && (j<=Par->jog)) 
		    flow_out+=DY[j]*DX[i]*W[i][j][k];
	    }
	}	
	flow_out = Par->CommSum(flow_out);
	mass_diff= flow_in-flow_out; 
    }
    
    if(S.iobd==2 && method!=-1){	
	for(p1=InOutList.GetNext();p1;p1=p1->GetNext()) {
	    i=p1->data.i; j=p1->data.j; k=p1->data.k;
	    SetVelocityCorrection(U,V,W,i,j,k);
	}
	corrected_mass_diff=0.0;
	for(p1=InOutList.GetNext();p1;p1=p1->GetNext()) {
	    i=p1->data.i; j=p1->data.j; k=p1->data.k;
	    flg=flag[i][j][k];
	    if((flg&SOUTH) && (i >Par->iug) && (j>=Par->jug) && (j<=Par->jog) && (k>=Par->kug) && (k<=Par->kog))
		corrected_mass_diff+=DZ[k]*DY[j]*U[i-1][j][k];
	    if((flg&NORTH) && (i<=Par->iog) && (j>=Par->jug) && (j<=Par->jog) && (k>=Par->kug) && (k<=Par->kog))
		corrected_mass_diff+=DZ[k]*DY[j]*U[i][j][k];
	    
	    if((flg&EAST)  && (j >Par->jug) && (i>=Par->iug) && (i<=Par->iog) && (k>=Par->kug) && (k<=Par->kog))
		corrected_mass_diff+=DZ[k]*DX[i]*V[i][j-1][k];
	    if((flg&WEST)  && (j<=Par->jog) && (i>=Par->iug) && (i<=Par->iog) && (k>=Par->kug) && (k<=Par->kog))
		corrected_mass_diff+=DZ[k]*DX[i]*V[i][j][k];
	    
	    if((flg&BOTTOM)&& (k >Par->kug) && (i>=Par->iug) && (i<=Par->iog) && (j>=Par->jug) && (j<=Par->jog))
		corrected_mass_diff+=DY[j]*DX[i]*W[i][j][k-1];
	    if((flg&TOP)   && (k<=Par->kog) && (i>=Par->iug) && (i<=Par->iog) && (j>=Par->jug) && (j<=Par->jog))
		corrected_mass_diff+=DY[j]*DX[i]*W[i][j][k];	    
	}
	corrected_mass_diff = Par->CommSum(corrected_mass_diff);
	corrected_mass_diff = flow_in-corrected_mass_diff;
    }    
}


// --- UNDER CONSTRUCTION ---
// Set normal boundary cond. for auxilliary fields UVWTC 
//
void NavierCalc::SetObstacleCondTilde(int method) {
//     int oldtimer=timer->Start(8,Par);
//     int                 i,j,k,n;
//     unsigned            flg;
//     NS_List<CELLREC>*      p;
//     NS_List<INOUTREC>*     p1;
//     NS_List<INFLOW>*       p2;
//     NS_List<SPECCELLREC>*  p3;
//     NS_List<INOUTRECDATA>* io3;
//    
//    timer->Stop(oldtimer,Par);
}
// --------------------------

/***************************************************************

    Set boundary Conditions for velocity after projection step

    we have to treat all tangential velocity components and 
    sometimes the normal velocity components

****************************************************************/                 

inline NS_REAL NavierCalc::STAGGVALUESLIP(Matrix<NS_REAL>& M,int DIR,int i,int j,int k) {
    switch(DIR) {
    case NORTH:
      return (flag[i][j][k]&NORTH) ? ((i<=Par->iog) ? M[i+1][j][k] : M[i][j][k]):
                                     ((i>=Par->iug) ? M[i-1][j][k] : M[i][j][k]);break;
    case WEST:
      return (flag[i][j][k]&WEST ) ? ((j<=Par->jog) ? M[i][j+1][k] : M[i][j][k]):
                                     ((j>=Par->jug) ? M[i][j-1][k] : M[i][j][k]);break;
    case TOP: 
      return (flag[i][j][k]&TOP  ) ? ((k<=Par->kog) ? M[i][j][k+1] : M[i][j][k]):
                                     ((k>=Par->kug) ? M[i][j][k-1] : M[i][j][k]);break;
    default:
      assert(0);break;
    }
    return 0;
}


inline NS_REAL NavierCalc::STAGGVALUENOSLIP(Matrix<NS_REAL>& M,int DIR,int i,int j,int k,NS_REAL value) {
    switch(DIR) {
    case NORTH:
      return (flag[i][j][k]&NORTH)?((i<=Par->iog) ? S.ddiv[0][0][i]*M[i+1][j][k]+value*(1-S.ddiv[0][0][i]) : M[i][j][k]):
	((i>=Par->iug) ? S.ddiv[1][0][i]*M[i-1][j][k]+value*(1-S.ddiv[1][0][i]) : M[i][j][k]);break;
    case WEST:
      return (flag[i][j][k]&WEST )?((j<=Par->jog) ? S.ddiv[0][1][j]*M[i][j+1][k]+value*(1-S.ddiv[0][1][j]) : M[i][j][k]):
	((j>=Par->jug) ? S.ddiv[1][1][j]*M[i][j-1][k]+value*(1-S.ddiv[1][1][j]) : M[i][j][k]);break;
    case TOP:
      return (flag[i][j][k]&TOP  )?((k<=Par->kog) ? S.ddiv[0][2][k]*M[i][j][k+1]+value*(1-S.ddiv[0][2][k]) : M[i][j][k]):
	((k>=Par->kug) ? S.ddiv[1][2][k]*M[i][j][k-1]+value*(1-S.ddiv[1][2][k]) : M[i][j][k]);break;
    default:
      assert(0);break;
    }
    return 0;
}


inline void NavierCalc::SetHomogeneousNeumannTangential(Matrix<NS_REAL>& U,Matrix<NS_REAL>& V,Matrix<NS_REAL>& W,
                                                        int i,int j,int k) {
    unsigned flg=flag[i][j][k];

    // du/dt=0
    //
    if(!(flg&(NORTH|SOUTH)) && (i<=Par->iog)) {                         // U is tangential
        if(!((flg&(TOP|BOTTOM)) && (flg&(EAST|WEST)))) {                // mixing not necessary
            if(flg&(TOP|BOTTOM))  U[i][j][k]=STAGGVALUESLIP(U,TOP ,i,j,k); else
                                  U[i][j][k]=STAGGVALUESLIP(U,WEST,i,j,k);
        } else {                                                        // two 'fluid faces' => mixing necessary
            U[i][j][k]=(DY[j]*STAGGVALUESLIP(U,TOP ,i,j,k)+
                        DZ[k]*STAGGVALUESLIP(U,WEST,i,j,k))/(DY[j]+DZ[k]);
        }
    }
    if(!(flg&(WEST|EAST)) && (j<=Par->jog)) {                           // V is tangential
        if(!((flg&(TOP|BOTTOM)) && (flg&(NORTH|SOUTH)))) {              // mixing not necessary
            if(flg&(TOP|BOTTOM))  V[i][j][k]=STAGGVALUESLIP(V,TOP  ,i,j,k); else
                                  V[i][j][k]=STAGGVALUESLIP(V,NORTH,i,j,k);
        } else {                                                        // mixing necessary
            V[i][j][k]=(DX[i]*STAGGVALUESLIP(V,TOP  ,i,j,k)+
                        DZ[k]*STAGGVALUESLIP(V,NORTH,i,j,k))/(DX[i]+DZ[k]);
        }
    }
    if(!(flg&(TOP|BOTTOM)) && (k<=Par->kog)) {                          // W is tangential
        if(!((flg&(NORTH|SOUTH)) && (flg&(EAST|WEST)))) {               // mixing not necessary
            if(flg&(NORTH|SOUTH)) W[i][j][k]=STAGGVALUESLIP(W,NORTH,i,j,k); else
                                  W[i][j][k]=STAGGVALUESLIP(W,WEST ,i,j,k);
        } else {                                                        // mixing necessary
            W[i][j][k]=(DY[j]*STAGGVALUESLIP(W,NORTH,i,j,k)+
                        DX[i]*STAGGVALUESLIP(W,WEST ,i,j,k))/(DX[i]+DY[j]);
        }
    }
}


inline void NavierCalc::SetDirichletTangential(Matrix<NS_REAL>& U,Matrix<NS_REAL>& V,Matrix<NS_REAL>& W,
                                                        int i,int j,int k,NS_REAL value[3]) {
        unsigned flg=flag[i][j][k];
        if(!(flg&(NORTH|SOUTH)) && (i<=Par->iog)) {                         // U is tangential
            if(!((flg&(TOP|BOTTOM)) && (flg&(EAST|WEST)))) {                // mixing not necessary
                if(flg&(TOP|BOTTOM))  U[i][j][k]=STAGGVALUENOSLIP(U,TOP ,i,j,k,value[0]); else
                                      U[i][j][k]=STAGGVALUENOSLIP(U,WEST,i,j,k,value[0]);
            } else {                                                        // mixing necessary
                U[i][j][k]=(DY[j]*STAGGVALUENOSLIP(U,TOP ,i,j,k,value[0])+
                            DZ[k]*STAGGVALUENOSLIP(U,WEST,i,j,k,value[0]))/(DY[j]+DZ[k]);
            }
        }
        if(!(flg&(WEST|EAST)) && (j<=Par->jog)) {                           // V is tangential
            if(!((flg&(TOP|BOTTOM)) && (flg&(NORTH|SOUTH)))) {              // mixing not necessary
                if(flg&(TOP|BOTTOM))  V[i][j][k]=STAGGVALUENOSLIP(V,TOP  ,i,j,k,value[1]); else
                                      V[i][j][k]=STAGGVALUENOSLIP(V,NORTH,i,j,k,value[1]);
            } else {                                                        // mixing necessary
                V[i][j][k]=(DX[i]*STAGGVALUENOSLIP(V,TOP  ,i,j,k,value[1])+
                            DZ[k]*STAGGVALUENOSLIP(V,NORTH,i,j,k,value[1]))/(DX[i]+DZ[k]);
            }
        }
        if(!(flg&(TOP|BOTTOM)) && (k<=Par->kog)) {                          // W is tangential
            if(!((flg&(NORTH|SOUTH)) && (flg&(EAST|WEST)))) {               // mixing not necessary
                if(flg&(NORTH|SOUTH)) W[i][j][k]=STAGGVALUENOSLIP(W,NORTH,i,j,k,value[2]); else
                                      W[i][j][k]=STAGGVALUENOSLIP(W,WEST ,i,j,k,value[2]);
            } else {                                                        // mixing necessary
                W[i][j][k]=(DY[j]*STAGGVALUENOSLIP(W,NORTH,i,j,k,value[2])+
                            DX[i]*STAGGVALUENOSLIP(W,WEST ,i,j,k,value[2]))/(DX[i]+DY[j]);
            }
        }
}


void NavierCalc::SetObstacleForTG(Matrix<NS_REAL>& M, int n) {
    int                 i,j,k;
    unsigned            flg;
    NS_List<INFLOW>*       p2;
    NS_List<SPECCELLREC>*  p3;
   
    // Treat Dirichlet cond. for temperature and chemicals 
    //
    if((n==-2) || (n>=0)){
      for(p2=InflowList.GetNext();p2;p2=p2->GetNext()) {
	// temperature and scalars
	// determine neighbouring fluid cell for interpolation of Dirichlet value
	// DirichletValueT == w*T[i][j][k] + w1*T[i1][j1][k1] , where w+w1=1 (linear interpolatation)
	i=p2->data.i; j=p2->data.j; k=p2->data.k ;
	int i1=-10000, j1=-10000, k1=-10000 ;
	double w=0.0,w1=0.0,ww1=0.0 ;
	switch(FREECELLS(flag[i][j][k])) {
        case SOUTH  :if(i>=Par->iug) {i1=i-1, j1=j  , k1=k  , w=S.d[0][i1], w1=S.d[0][i];}
	  break;
        case NORTH  :if(i<=Par->iog) {i1=i+1, j1=j  , k1=k  , w=S.d[0][i1], w1=S.d[0][i];}
	  break;
        case EAST   :if(j>=Par->jug) {i1=i  , j1=j-1, k1=k  , w=S.d[1][j1], w1=S.d[1][j];}
        break;
        case WEST   :if(j<=Par->jog) {i1=i  , j1=j+1, k1=k  , w=S.d[1][j1], w1=S.d[1][j];}
	  break;
        case BOTTOM :if(k>=Par->kug) {i1=i  , j1=j  , k1=k-1, w=S.d[2][k1], w1=S.d[2][k];}
	  break;
        case TOP    :if(k<=Par->kog) {i1=i  , j1=j  , k1=k+1, w=S.d[2][k1], w1=S.d[2][k];}
	  break;
        default: {i1=i, j1=j, k1=k; w=w1=0.5;}
	}
	ww1=w+w1 ;
	w =w /ww1 ;
	w1=w1/ww1 ;
	
	if(n==-2)   M[i][j][k]=(p2->data.t        -  M[i1][j1][k1]*w1)/w;   // Temperature
	if(n>=0)    M[i][j][k]=(p2->data.chem[n]  -  M[i1][j1][k1]*w1)/w;   // Chemicals
      }
    }
    

    // Neumann conditions for temparature and chemicals 
    //
    if(n==-2) for(p3=dTempdnList.GetNext();p3;p3=p3->GetNext()) {       // Temperature
        i=p3->data.i;j=p3->data.j;k=p3->data.k;
        flg=flag[i][j][k];
        switch(FREECELLS(flg)) {
        case SOUTH  :if(i>=Par->iug) M[i][j][k]=M[i-1][j][k];break;
        case NORTH  :if(i<=Par->iog) M[i][j][k]=M[i+1][j][k];break;
        case EAST   :if(j>=Par->jug) M[i][j][k]=M[i][j-1][k];break;
        case WEST   :if(j<=Par->jog) M[i][j][k]=M[i][j+1][k];break;
        case BOTTOM :if(k>=Par->kug) M[i][j][k]=M[i][j][k-1];break;
        case TOP    :if(k<=Par->kog) M[i][j][k]=M[i][j][k+1];break;
        default:    // treatment for more than one neighbouring cells (weighting)
            M[i][j][k]=0;
            if((flg&SOUTH ) && (i>=Par->iug)) M[i][j][k]+=M[i-1][j][k]*(p3->data.gweight/DX[i]);
            if((flg&NORTH ) && (i<=Par->iog)) M[i][j][k]+=M[i+1][j][k]*(p3->data.gweight/DX[i]);
            if((flg&EAST  ) && (j>=Par->jug)) M[i][j][k]+=M[i][j-1][k]*(p3->data.gweight/DY[j]);
            if((flg&WEST  ) && (j<=Par->jog)) M[i][j][k]+=M[i][j+1][k]*(p3->data.gweight/DY[j]);
            if((flg&BOTTOM) && (k>=Par->kug)) M[i][j][k]+=M[i][j][k-1]*(p3->data.gweight/DZ[k]);
            if((flg&TOP   ) && (k<=Par->kog)) M[i][j][k]+=M[i][j][k+1]*(p3->data.gweight/DZ[k]);
            M[i][j][k]/=p3->data.den;break;
        }
    } 
    if(n>=0) for(p3=dChemdnList.GetNext();p3;p3=p3->GetNext()) {      // Chemicals
        i=p3->data.i;j=p3->data.j;k=p3->data.k;
        flg=flag[i][j][k];
        switch(FREECELLS(flg)) {
        case SOUTH : if(i>=Par->iug) M[i][j][k]=M[i-1][j][k];break;
        case NORTH : if(i<=Par->iog) M[i][j][k]=M[i+1][j][k];break;
        case EAST  : if(j>=Par->jug) M[i][j][k]=M[i][j-1][k];break;
        case WEST  : if(j<=Par->jog) M[i][j][k]=M[i][j+1][k];break;
        case BOTTOM: if(k>=Par->kug) M[i][j][k]=M[i][j][k-1];break;
        case TOP   : if(k<=Par->kog) M[i][j][k]=M[i][j][k+1];break;
        default:    // treatment for more than one neighbouring cells (weighting)
	    M[i][j][k]=0;
	    if((flg&SOUTH ) && (i>=Par->iug))  M[i][j][k]+=M[i-1][j][k]*(p3->data.gweight/DX[i]);
	    if((flg&NORTH ) && (i<=Par->iog))  M[i][j][k]+=M[i+1][j][k]*(p3->data.gweight/DX[i]);
	    if((flg&EAST  ) && (j>=Par->jug))  M[i][j][k]+=M[i][j-1][k]*(p3->data.gweight/DY[j]);
	    if((flg&WEST  ) && (j<=Par->jog))  M[i][j][k]+=M[i][j+1][k]*(p3->data.gweight/DY[j]);
	    if((flg&BOTTOM) && (k>=Par->kug))  M[i][j][k]+=M[i][j][k-1]*(p3->data.gweight/DZ[k]);
	    if((flg&TOP   ) && (k<=Par->kog))  M[i][j][k]+=M[i][j][k+1]*(p3->data.gweight/DZ[k]);
	    M[i][j][k]/=p3->data.den;
	    break;
        }
    }
}


void NavierCalc::SetObstacleCond(int method) {
    int oldtimer=timer->Start(8,Par);
    int                 i,j,k;
    unsigned            flg;
    NS_List<CELLREC>*      p;
    NS_List<INOUTREC>*     p1;
    NS_List<INFLOW>*       p2;
    NS_List<SPECCELLREC>*  p3;
    NS_List<INOUTRECDATA>* io3;
    
    // #######################################
    // # Set normal boundary cond. for UVWTC #
    // #######################################
    
    // Treat Dirichlet Cond. 
    //
    for(p2=InflowList.GetNext();p2;p2=p2->GetNext()) {
	// Set dirichlet normal boundary cond. for velocity u,v,w
	i=p2->data.i; j=p2->data.j; k=p2->data.k;
	SetDirichletNormal(U,V,W, i,j,k, p2->data.u, p2->data.v, p2->data.w);
    }
    

    // #####################################
    // # Set normal boundary cond. for UVW #
    // #####################################

    // Set Slip dirichlet normal boundary cond.
    for(p=SlipList.GetNext();p;p=p->GetNext()) 
	SetDirichletNormal(U,V,W,p->data.i,p->data.j,p->data.k,0,0,0);
    
    // Set NoSlip dirichlet normal boundary cond.
    for(p=NoslipList.GetNext();p;p=p->GetNext()) 
	SetDirichletNormal(U,V,W,p->data.i,p->data.j,p->data.k,0,0,0);
    
    //
    // Treat InOut boundary cond. in normal direction 
    //
    for(p1=InOutList.GetNext();p1;p1=p1->GetNext()) {
	i=p1->data.i;j=p1->data.j;k=p1->data.k; 
	switch(p1->data.type) {
	case 1:       // du/dn=0; simple Neuman zero boundary cond.
	    SetHomogeneousNeumannNormal(U,V,W,i,j,k) ;
	    break;
	case 2:	      // time-dependent convective boundary cond. 
	    SetHomogeneousNeumannNormal(U,V,W,i,j,k) ;
	    break;
	} 
    }    
   
//     --- UNDER CONSTRUCTION --- 
//     BoundaryCondition 3 is implemented somewhat different 
//     here just for 1st order Euler        
//     for(io3=InOut3List.GetNext();io3;io3=io3->GetNext()) {
//       i=io3->data.i;j=io3->data.j;k=io3->data.k;
//       flg=flag[i][j][k];
//       if((flg&SOUTH) && (i >Par->iug)) U[i-1][j][k]+=S.delt*(S.g[0] - io3->data.param*io3->data.data[0]);
//       if((flg&NORTH) && (i<=Par->iog)) U[i+1][j][k]+=S.delt*(S.g[0] - io3->data.param*io3->data.data[0]);
      
//       if((flg&EAST ) && (j >Par->jug)) V[i][j-1][k]+=S.delt*(S.g[1] - io3->data.param*io3->data.data[1]);
//       if((flg&WEST ) && (j<=Par->jog)) V[i][j+1][k]+=S.delt*(S.g[1] - io3->data.param*io3->data.data[1]);

//       if((flg&BOTTOM)&& (k >Par->kug)) W[i][j][k-1]+=S.delt*(S.g[2] - io3->data.param*io3->data.data[2]);
//       if((flg&TOP  ) && (k<=Par->kog)) W[i][j][k+1]+=S.delt*(S.g[2] - io3->data.param*io3->data.data[2]);
//     } 
//     --------------------------
    
    // #########################################
    // # Set tangential boundary cond. for UVW #
    // #########################################
 
    //
    // Treat Dirichlet BC 
    // 
    for(p2=InflowList.GetNext();p2;p2=p2->GetNext()) {  
        NS_REAL uvwvec[3]={p2->data.u,p2->data.v,p2->data.w};
        SetDirichletTangential(U,V,W,p2->data.i,p2->data.j,p2->data.k,uvwvec);
    }

    NS_REAL uvwvec[3]={0,0,0};
    for(p=NoslipList.GetNext();p;p=p->GetNext()) 
      SetDirichletTangential(U,V,W,p->data.i,p->data.j,p->data.k,uvwvec);
  
    //
    // Treat Neumann BC for Slip BC ; slightly simplified discretization ...
    //
    for(p=SlipList.GetNext();p;p=p->GetNext()) 
      SetHomogeneousNeumannTangential(U,V,W,p->data.i,p->data.j,p->data.k);  

    //
    //  Treat InOut BC 1 and BC 2
    //
    for(p1=InOutList.GetNext(); p1; p1=p1->GetNext()) {
        i=p1->data.i; j=p1->data.j; k=p1->data.k; flg=flag[i][j][k];
        switch(p1->data.type) {
        case 1: // simple neumann zero condition
	    SetHomogeneousNeumannTangential(U,V,W,i,j,k);
	    break;
        case 2: // time-dependent convective boundary cond.
	    SetHomogeneousNeumannTangential(U,V,W,i,j,k);
	    SetVelocityExtrapolation(U,V,W,i,j,k) ;
	    break;
        }
    }

    // correction of boundary conditions (optional)
    OutflowVelocityCorrection(U,V,W,method);

   
//     --- UNDER CONSTRUCTION ---
    
//     initialize boundary values for tangential velocity components at places where BC 3 is prescribed 
     
//     if (method==-1) for(io3=InOut3List.GetNext(); io3; io3=io3->GetNext()) {
//         i=io3->data.i; j=io3->data.j; k=io3->data.k;
//         flg=flag[i][j][k];
//         if((flg&SOUTH) && (i >Par->iug)) {
//             V[i][j][k]=V[i-1][j][k];  W[i][j][k]=W[i-1][j][k];
//             io3->data.data[1]=io3->data.data[2]=0;
//         }
//         if((flg&NORTH) && (i<=Par->iog)) {
//             V[i][j][k]=V[i+1][j][k];  W[i][j][k]=W[i+1][j][k];
//             io3->data.data[1]=io3->data.data[2]=0;
//         }
//         if((flg&EAST ) && (j >Par->jug)) {
//             U[i][j][k]=U[i][j-1][k];  W[i][j][k]=W[i][j-1][k];
//             io3->data.data[0]=io3->data.data[2]=0;
//         }
//         if((flg&WEST ) && (j<=Par->jog)) {
//             U[i][j][k]=U[i][j+1][k];  W[i][j][k]=W[i][j+1][k];
//             io3->data.data[0]=io3->data.data[2]=0;
//         }
//         if((flg&BOTTOM)&& (k >Par->kug)) {
//             U[i][j][k]=U[i][j][k-1];  V[i][j][k]=V[i][j][k-1];
//             io3->data.data[0]=io3->data.data[1]=0;
//         }
//         if((flg&TOP  ) && (k<=Par->kog)) {
//             U[i][j][k]=U[i][j][k+1];  V[i][j][k]=V[i][j][k+1];
//             io3->data.data[0]=io3->data.data[1]=0;
//         }
//     }

//     --- UNDER CONSTRUCTION --- 
//     set BC 3 for all time-steps, this amounts to store the normal gradients for all
//     velocity components and to set the BC for the tangential vel. components the BC 
//     for the normal vel. comp. is implicitely set by SetObstacleCondTilde and  AdapUVW 
    
//     for(io3=InOut3List.GetNext();io3;io3=io3->GetNext()) {
// 	i=io3->data.i; j=io3->data.j; k=io3->data.k;
// 	flg=flag[i][j][k];
	
//         if((flg&SOUTH) && (i >Par->iug)) {
//             V[i  ][j][k] += S.delt*(S.g[1]-io3->data.param*io3->data.data[1]);
//             W[i  ][j][k] += S.delt*(S.g[2]-io3->data.param*io3->data.data[2]);
//             io3->data.data[0]=(U[i-1][j][k]-U[i-2][j][k])/DX [i-1];
//             io3->data.data[1]=(V[i  ][j][k]-V[i-1][j][k])/DXM[i-1];
//             io3->data.data[2]=(W[i  ][j][k]-W[i-1][j][k])/DXM[i-1];
//         }
//         if((flg&NORTH) && (i<=Par->iog)) {
//             V[i  ][j][k] += S.delt*(S.g[1]-io3->data.param*io3->data.data[1]);
//             W[i  ][j][k] += S.delt*(S.g[2]-io3->data.param*io3->data.data[2]);
//             io3->data.data[0]=(U[i  ][j][k]-U[i+1][j][k])/DX [i+1];
//             io3->data.data[1]=(V[i+1][j][k]-V[i  ][j][k])/DXM[i  ];
//             io3->data.data[2]=(W[i+1][j][k]-W[i  ][j][k])/DXM[i  ];
//         }
//         if((flg&EAST ) && (j >Par->jug)) {
//             U[i  ][j][k] += S.delt*(S.g[0]-io3->data.param*io3->data.data[0]);
//             W[i  ][j][k] += S.delt*(S.g[2]-io3->data.param*io3->data.data[2]);
//             io3->data.data[0]=(U[i  ][j][k]-U[i][j-1][k])/DYM[j-1];
//             io3->data.data[1]=(V[i][j-1][k]-U[i][j-2][k])/DY [j-1];
//             io3->data.data[2]=(W[i  ][j][k]-W[i][j-1][k])/DYM[j-1];
//         }
//         if((flg&WEST ) && (j<=Par->jog)) {
//             U[i  ][j][k] += S.delt*(S.g[0]-io3->data.param*io3->data.data[0]);
//             W[i  ][j][k] += S.delt*(S.g[2]-io3->data.param*io3->data.data[2]);
//             io3->data.data[0]=(U[i][j+1][k]-U[i  ][j][k])/DYM[j  ];
//             io3->data.data[1]=(V[i  ][j][k]-V[i][j+1][k])/DY [j+1];
//             io3->data.data[2]=(W[i][j+1][k]-W[i  ][j][k])/DYM[j  ];
//         }
//         if((flg&BOTTOM)&& (k >Par->kug)) {
//             U[i  ][j][k] += S.delt*(S.g[0]-io3->data.param*io3->data.data[0]);
//             V[i  ][j][k] += S.delt*(S.g[1]-io3->data.param*io3->data.data[1]);
//             io3->data.data[0]=(U[i  ][j][k]-U[i][j][k-1])/DZM[k-1];
//             io3->data.data[1]=(V[i  ][j][k]-V[i][j][k-1])/DZM[k-1];
//             io3->data.data[2]=(W[i][j][k-1]-W[i][j][k-2])/DZ [k-1];
//         }
//         if((flg&TOP  ) && (k<=Par->kog)) {
//             U[i  ][j][k] += S.delt*(S.g[0]-io3->data.param*io3->data.data[0]);
//             V[i  ][j][k] += S.delt*(S.g[1]-io3->data.param*io3->data.data[1]);
//             io3->data.data[0]=(U[i][j][k+1]-U[i  ][j][k])/DZM[k  ];
//             io3->data.data[1]=(V[i][j][k+1]-V[i  ][j][k])/DZM[k  ];
//             io3->data.data[2]=(W[i  ][j][k]-W[i][j][k+1])/DZ [k+1];
//         }
//     }         
//     -----------------------------------

    timer->Stop(oldtimer,Par);
}

