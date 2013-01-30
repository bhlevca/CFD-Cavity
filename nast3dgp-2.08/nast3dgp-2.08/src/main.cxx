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
#include <stdlib.h>
#include <string.h>


const char* procnames[]={"",
    "NavierCalc::DoIt()",
    "NavierCalc::AdapUVW()",
    "NavierCalc::CompTG()",
    "NavierCalc::CompFGH()",
    "NavierCalc::CompTUVWfromFGH()",
    "NavierCalc::SetObstacleCondTilde()", 
    "NavierCalc::CompRHS()",
    "NavierCalc::TimeStep()",
    "NavierCalc::SetPBorder()",
    "NavierCalc::SetPALLBorder()",
    "NavierCalc::SetObstacleCond()",
    "NavierCalc::WriteFile()",
    "NavierCalc::CalcMeanVal()", 
    "NavierCalc::CalcRes()",
    "NavierCalc::Poisson()",
    "NavierCalc::BiCGStab()",
    "ParParams::CommMin()",
    "ParParams::CommMax()",
    "ParParams::CommSum()",
    "ParParams::SendColour()",
    "ParParams::ReceiveColour()",
    "ParParams::CommMatrix()",
    "ParParams::CommMatrixFlag()" 
    };

// 
// parse command line
//

void NavierCalc::ParseArgs(int argc,char** argv) {
    int      error=FALSE;
    char*    env_str ;

    while (--argc>0)
        if (argv[argc][0]=='-')
            switch (argv[argc][1]) {
            case 'b':                  
                binfile[0]=0 ; 
                env_str=getenv("NAV_BIN");
                if(env_str) {if(argv[argc+1][0]!='/' && argv[argc+1][0]!='\\') strcpy(binfile,env_str);}
                strcat(binfile, argv[argc+1]) ;
                break;
            case 'p':
                printinfo=atoi(argv[argc]+2);

                if(printinfo<0) printinfo=0;
                if(printinfo>4) printinfo=4;
                break;
            case 't':
                time_to_run=atoi(argv[argc+1]);
                break;
            case 'h':
            case '?':
                error=TRUE;
                break;
            default:
                break;
            }
    if (error==TRUE) {
        printf("This is NAVCALC version %s\nParameter:\n -b <file> binary data file\n"
               " -p <n>    info level\n"
               " -t <n>    number of seconds to run\n", VERSION);
        exit(1);
    }
}    

NavierCalc::~NavierCalc() {  // Destructor
    if(CH)  delete[] CH;
    if(Fch) delete[] Fch;
    delete Par; 
    delete timer;
}

NavierCalc::NavierCalc(int* pargc,char*** pargv):Navier() {
    int i,j,k;
    time_to_run=0;
    printinfo=1  ;
    actual_fgh=0 ;
    Fch=NULL ;

    Par=new ParParams(pargc,pargv);                         // Initialize parallelization
    timer=new Timer               ;                          
    Par->SetTimer(timer)          ;
    
    ParseArgs(*pargc,*pargv)      ;                         // Read Command Line Args  

    iterations=alliterations=n=0;
    /*if(Par->IsParent())*/ {                               // Read parts of data file for global info 
        if(!ReadScene(binfile)) { fprintf(stderr,"Cannot read file %s.\n",binfile);exit(1);}
        if(printinfo>=2) if(Par->IsParent()) {
            printf("Reading file %s...\n",binfile);
            S.DumpInfo(stdout);
        }
    }

    // Distribute global info such as number of gridpoints, filename and number of chemicals to be computed    
    // Par->ShareInfo(S,binfile);
    Par->SetBounds(S);                                          // Domain decomposition for parallelization

    if (S.Solver == 3 || S.Solver == 4) Par->CalcColourCount(); // Prepare for coloured send/receive

    if(printinfo>=2) {
        Par->PrintInitInfo();
        Par->Queue();
        printf("Process %d (%d-%d, %d-%d, %d-%d)  CompTemp:%d  CompChem:%d  nChem:%d\n",
            Par->Me(),Par->iug,Par->iog,Par->jug,Par->jog,Par->kug,Par->kog,S.CompTemp,S.CompChem,S.nchem);
        fflush(stdout);
        Par->QueueNext();
    }
    
    // Initialize arrays   
    //
    FILE* f=fopen(binfile,"rb");
    if(f==0) { fprintf(stderr,"Cannot open file %s.\n",binfile); exit(1);}
    S.ReadFile(f);

    if((printinfo>=2) && Par->IsParent()) { printf ("Init. Arrays\n") ; fflush(stdout) ;}

    int  dims[6]={Par->iug-1,Par->iog+1,Par->jug-1,Par->jog+1,Par->kug-1,Par->kog+1}; 
    int dims1[6]={Par->iug,Par->iog,Par->jug,Par->jog,Par->kug,Par->kog};
    int dims2[6]={Par->iug-(Par->iug>1?S.GC:1),
                  Par->iog+(Par->iog<S.gridp[0]?S.GC:1),
                  Par->jug-(Par->jug>1?S.GC:1),
                  Par->jog+(Par->jog<S.gridp[1]?S.GC:1),
                  Par->kug-(Par->kug>1?S.GC:1),
                  Par->kog+(Par->kog<S.gridp[2]?S.GC:1)};
    int dims4[6]={Par->iug-4, Par->iog+4,   // ghostcells for level-set (ENO 4th order) 
	          Par->jug-4, Par->jog+4,
	          Par->kug-4, Par->kog+4};

    if((Par->iug==1)          && (S.periodbound&1)) {dims2[0]-=(S.GC-1);}
    if((Par->iog==S.gridp[0]) && (S.periodbound&1)) {dims2[1]+=(S.GC-1);}
    if((Par->jug==1)          && (S.periodbound&2)) {dims2[2]-=(S.GC-1);}
    if((Par->jog==S.gridp[1]) && (S.periodbound&2)) {dims2[3]+=(S.GC-1);}
    if((Par->kug==1)          && (S.periodbound&4)) {dims2[4]-=(S.GC-1);}
    if((Par->kog==S.gridp[2]) && (S.periodbound&4)) {dims2[5]+=(S.GC-1);}

    if((dims2[0]>0) || (S.periodbound&1)) {dims2[0]--;U.Init(dims2);dims2[0]++;} else U.Init(dims2);
    if((dims2[2]>0) || (S.periodbound&2)) {dims2[2]--;V.Init(dims2);dims2[2]++;} else V.Init(dims2);
    if((dims2[4]>0) || (S.periodbound&4)) {dims2[4]--;W.Init(dims2);dims2[4]++;} else W.Init(dims2);

    flag.Init(dims2);

    P.Init(dims);
    F[0].Init(dims); 
    G[0].Init(dims);
    H[0].Init(dims);
    F[1].Init(dims); 
    G[1].Init(dims);
    H[1].Init(dims);
    RHS.Init(dims1);
    Tmp1.Init(dims4); // temporary field

    // arrays for BiCGSTAB and BiCGSTAB+JPC Solvers
    //
    if (S.Solver == 5 || S.Solver == 6) { 	
	BiCG_r0.Init(dims1); BiCG_rj.Init(dims1); 
	BiCG_pj.Init(dims1); BiCG_vj.Init(dims1);
	BiCG_sj.Init(dims1); BiCG_tj.Init(dims1); 
	MatBuf.Init(dims); KoeffPreCond.Init(dims1);
    }

    if(S.CompChem || S.CompTemp) T2.Init(dims);
    if(S.CompTemp) {T.Init(dims2); Ft.Init(dims);}
    if(S.CompChem) {
      CH =new Matrix<NS_REAL>[S.nchem];
      for(i=0;i<S.nchem;i++) CH[i].Init(dims2);
      Fch=new Matrix<NS_REAL>[S.nchem];for(i=0;i<S.nchem;i++) Fch[i].Init(dims);
    }

    if((Par->iug==1)          && (S.periodbound&1)) {dims2[0]+=(S.GC-1);}
    if((Par->iog==S.gridp[0]) && (S.periodbound&1)) {dims2[1]-=(S.GC-1);}
    if((Par->jug==1)          && (S.periodbound&2)) {dims2[2]+=(S.GC-1);}
    if((Par->jog==S.gridp[1]) && (S.periodbound&2)) {dims2[3]-=(S.GC-1);}
    if((Par->kug==1)          && (S.periodbound&4)) {dims2[4]+=(S.GC-1);}
    if((Par->kog==S.gridp[2]) && (S.periodbound&4)) {dims2[5]-=(S.GC-1);} 

    // Read values flag/U/V/W/... from file
    //
    flag.ReadPartial(f,dims2);                                  
    U.ReadPartial(f,dims2);
    V.ReadPartial(f,dims2);
    W.ReadPartial(f,dims2);
    P.ReadPartial(f);
    if(S.CompTemp) T.ReadPartial(f,dims);
    if(S.CompChem) for(i=0;i<S.nchem;i++) CH[i].ReadPartial(f,dims);
    fclose(f); 

    // the only sequential file access 
    //
    S.vol_of_iobc1_cells=0;   // calc parameter for IOBC1 
    NS_List<INOUTREC>* p=S.InOutList.GetNext();
    while(p) {
        i=p->data.i;j=p->data.j;k=p->data.k;
        if(p->data.type==1 || p->data.type==2) if(PPINALLBORDER(i,j,k)) switch(FREECELLS(flag[i][j][k])) {
            case NORTH: if(PINBORDER(i+1,j,k)) S.vol_of_iobc1_cells+=DZ[k]*DY[j];break;
            case SOUTH: if(PINBORDER(i-1,j,k)) S.vol_of_iobc1_cells+=DZ[k]*DY[j];break;
            case WEST:  if(PINBORDER(i,j+1,k)) S.vol_of_iobc1_cells+=DX[i]*DZ[k];break;
            case EAST:  if(PINBORDER(i,j-1,k)) S.vol_of_iobc1_cells+=DX[i]*DZ[k];break;
            case BOTTOM:if(PINBORDER(i,j,k-1)) S.vol_of_iobc1_cells+=DX[i]*DY[j];break;
            case TOP:   if(PINBORDER(i,j,k+1)) S.vol_of_iobc1_cells+=DX[i]*DY[j];break;
        }
        p=p->GetNext();
    }

    // calculate volumes
    S.vol_of_iobc1_cells=Par->CommSum(S.vol_of_iobc1_cells);
    printf("iobc1_vol: %e\n",S.vol_of_iobc1_cells);fflush(stdout);

    vol_domain=(S.kabs[0][S.gridp[0]] - S.kabs[0][0])*
               (S.kabs[1][S.gridp[1]] - S.kabs[1][0])*
               (S.kabs[2][S.gridp[2]] - S.kabs[2][0]);

    if(Par->IsParent()) printf("volume of domain = %e\n",vol_domain); fflush(stdout);
   
    // Generate local lists with cells of different types of boundary conditions
    //
    if(printinfo>=2) {printf("Build local Lists  %d\n",Par->Me()); fflush(stdout);}
    BuildLocalLists();
    if(printinfo>=2) {printf ("Initialization done %d\n",Par->Me()); fflush(stdout);}
}


int NavierCalc::WriteMatrix(FILE* f,Matrix<NS_REAL>& mat,const int* dims) {
#ifdef NOPARALLEL
    return mat.WritePartial(f,dims);
#else
    if(Par->IsParent()) {
        int ret=mat.WriteParallel(f,ParParams::writevect,Par,dims); 
        Par->EndWriteMatrix();
        return ret;
    } else Par->ChildWriteMatrix(mat,dims);
    return TRUE;
#endif
}

int NavierCalc::WriteFile(char* filename) { // Every processor recognizes, at which position in the file 
    int oldtimer=timer->Start(9,Par);       // the data has to be written
    int i,j,k;
    FILE* f=NULL ;

    if(Par->IsParent()) { // only master processor 
        f=fopen(filename,"rb+");
        if(!f) return FALSE;
        if(!S.WriteFile(f)) return FALSE;
        if(!flag.Skip(f)) return FALSE;
    }

    int dims[6]={(Par->iug==1) ? 0:Par->iug, (Par->iog==S.gridp[0]) ? Par->iog+1:Par->iog,
                 (Par->jug==1) ? 0:Par->jug, (Par->jog==S.gridp[1]) ? Par->jog+1:Par->jog,
                 (Par->kug==1) ? 0:Par->kug, (Par->kog==S.gridp[2]) ? Par->kog+1:Par->kog};

    if(!WriteMatrix(f,U,dims)) return FALSE;
    if(!WriteMatrix(f,V,dims)) return FALSE;
    if(!WriteMatrix(f,W,dims)) return FALSE;
    if(S.delt_old > 0.0) {
	PIJKALLLOOP Tmp1[i][j][k]=P[i][j][k]/S.delt_old; // calc original values from delt-scaled pressure
	if(!WriteMatrix(f,Tmp1,dims)) return FALSE;  // write original pressure field
    } else
	if(!WriteMatrix(f,P,dims)) return FALSE;  // write pressure field
    if(S.CompTemp) if(!WriteMatrix(f,T,dims)) return FALSE;
    if(S.CompChem) for(int i=0;i<S.nchem;i++) if(!WriteMatrix(f,CH[i],dims)) return FALSE;
    if(Par->IsParent()) fclose(f);
    
    timer->Stop(oldtimer,Par);
    return TRUE;
}


int NavierCalc::ReadScene(char* filename) {
    FILE* f=fopen(filename,"rb");
    if(f==0) return FALSE ;
    if(!S.ReadFile(f)) return FALSE;
    fclose(f);
    return TRUE;
}


#define cmp(mat,ii,jj,kk,iii,jjj,kkk) if(mat[ii][jj][kk]!=mat[iii][jjj][kkk]) printf("%2d %2d %2d %f %2d %2d %2d %f %f\n",ii,jj,kk,mat[ii][jj][kk],iii,jjj,kkk,mat[iii][jjj][kkk],mat[ii][jj][kk]-mat[iii][jjj][kkk]);


void NavierCalc::DoIt() {      // Main Loop          
    timer->Start(1,Par);
    NS_REAL    res=0;
    double     starttime,actualtime,globalstarttime ;
    double     time_in_prev_step,time_in_this_step,timediff;
    int        i,j,k,l,n;

    S.delt=0.0;
    S.timestepmethod=0;
    
    // Initialize auxiliary fields
    for(l=0; l<2; l++) PIJKALLLOOP F[l][i][j][k]=G[l][i][j][k]=H[l][i][j][k]=0;
    PIJKALLLOOP P[i][j][k]*=S.delt_old;

    // One time only flag-field communication
    Par->CommMatrixFlag(flag,0,0,0,flagCOMM, S.GC); 

    SetObstacleCond(-1);   

    Par->CommMatrix(U, 1, 0, 0, UCOMM, S.GC);  
    Par->CommMatrix(V, 0, 1, 0, VCOMM, S.GC);
    Par->CommMatrix(W, 0, 0, 1, WCOMM, S.GC);
    
    time_in_prev_step=starttime=actualtime=Par->GetTime();
    globalstarttime=Par->CommMax(starttime);
    printf("Clock skews(%d) %e\n",Par->Me() , globalstarttime-starttime) ;fflush(stdout);

    //###################################################################################
    //##############################               ######################################
    //##############################   MAIN LOOP   ######################################
    //##############################               ######################################
    //###################################################################################
    do {        	                       
        			
	S.delt=TimeStep();                         // Compute new time step 

	SetObstacleCond(S.timestepmethod);         // Set Boundary conditions for new velocity including time dep. bc	

	Par->CommMatrix(U, 1, 0, 0, UCOMM, S.GC);  // Parallel communication for U-velocity 
	Par->CommMatrix(V, 0, 1, 0, VCOMM, S.GC);  // Parallel communication for V-velocity 
	Par->CommMatrix(W, 0, 0, 1, WCOMM, S.GC);  // Parallel communication for W-velocity 

	if(S.CompTemp) {SetObstacleForTG(T,-2); Par->CommMatrix(T, 0, 0, 0, TCOMM, S.GC);}
	if(S.CompChem) for(n=0; n<S.nchem; n++) {SetObstacleForTG(CH[n],n); Par->CommMatrix(CH[n],0,0,0,CHCOMM+n,S.GC);}

	OutputData(); S.t+=S.delt;          // Data writing routine and calculation of new physical time

	if(S.CompTemp) {                                            // if required advect+diffuse temperature explicitely
	    CompTG(T,S.alphatg,S.nu/S.prandtl);                     // Transport Temperature (space)
	    CompTT(T,Ft,S.alphatg,S.nu/S.prandtl,-2,S.timestepmethod);// Transport Temperature (time)
	    Par->CommMatrix(T, 0, 0, 0, TCOMM, S.GC);               // Parallel communication for Temperature
	}
	if(S.CompChem) {                                   // if required advect+diffuse chemicals explicitely
	    for(n=0; n<S.nchem; n++) {                     // Loop for all chemicals
		CompTG(CH[n],S.alphatg,S.chemc[n]);        // Transport of chemicals (space)
		CompTT(CH[n],Fch[n],S.alphatg,S.chemc[n],n,S.timestepmethod);// Transport chemicals (time)
		Par->CommMatrix(CH[n],0,0,0,CHCOMM+n,S.GC);// Parallel communication for chemicals
 	    }    
	}
	
	CompFGH();                               // Compute predicted velocities
	
	CompTUVWfromFGH(S.timestepmethod);       // Compute predicted velocities with 2nd order Adams-Bashfort or Runge-Kutta 
	CompRHS();                               // Compute right hand side for poisson equation
	res=CompPoisson(starttime);              // Solve poisson equation  
	AdapUVW();                               // Compute new velocity (velocity correction by pressure-gradient)        
	InfoOutputs(res);                        // Some Info outputs

	
	S.delt_old=S.delt;                       // old time step needed for adams-bashfort          
	S.timestepmethod|=1;                     // used to mark the first time step for AB discretization and data output 
	if(S.TimeDis==1) actual_fgh^=1;          // store old values for Adams-Bashforth       
	S.timesteps++;

	time_in_this_step=Par->GetTime();
	time_in_this_step=Par->CommMax(time_in_this_step);
	timediff=time_in_this_step-time_in_prev_step;
	time_in_prev_step=time_in_this_step;
	timediff=Par->CommMax(timediff);
	
    } while (((S.Tfin < 0) ? (res>S.eps):(S.Tfin>S.t)) && 
	     (time_to_run ? (time_in_this_step+2*timediff-starttime<time_to_run):TRUE));
    
    if(Par->IsParent() && (printinfo>=2)) printf("----------------------------------------\n");

    S.delt=TimeStep();       // Compute last time step for outflow boundary cond.
    SetObstacleCond(S.timestepmethod);  // Set boundary conditions
    
    Par->CommMatrix(U, 1, 0, 0, UCOMM, S.GC);
    Par->CommMatrix(V, 0, 1, 0, VCOMM, S.GC);
    Par->CommMatrix(W, 0, 0, 1, WCOMM, S.GC);
    
    if(!WriteFile(binfile)) {fprintf(stderr,"Error writing file %s.\n",binfile);exit(1);}
    
    timer->Stop(0,Par); 
    if(printinfo>=2) timer->Dump(procnames);
}

int main(int argc, char *argv[]) { 
    NavierCalc NaSt3DGPF(&argc,&argv);
    NaSt3DGPF.DoIt();
    return(0);
}
