/* NaSt3DGP - The Parallel 3D Navier-Stokes Solver
 * Copyright (C) 2003 Institute for Numerical Simulation
 *                    University of Bonn
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 */

#ifndef NOPARALLEL
#include "typen.h"
#include "matrix.h"
#include "parallel.h"
#include <string.h>
#include <stdlib.h>

//
// Queue-Operations for quasi sequential operations
//
// Each process gets the ID of the (0,0,0) in the cartesian proc-topology
// then, in a loop process (0,0,0) sends to all others prosesses a message,
//                                 but continues in the loop only after it has received an answer
//                                 from actual addressee   
// the other processes wait to receive their message and then leave Queue. 
// The answer message is send by running QueueNext       

void ParParams::Queue() {
      int        i,j,k,l,m;
      MPI_Status status;
      int p[3]={0,0,0};
      MPI_Cart_rank(MYCOMM,p,&m) ; 
      if (me==m) {
         for (i=pproc[0]-1; i>=0; i--) for(j=pproc[1]-1;j>=0;j--) for(k=pproc[2]-1;k>=0;k--) {
             int pos[3]={i,j,k};
             if((i!=0) || (j!=0) || (k!=0)) {   
             MPI_Cart_rank(MYCOMM,pos,&m);
             MPI_Send(&i,1,MPI_INT,m,QUEUECOMM,MYCOMM);
             MPI_Recv(&l,1,MPI_INT,m,QUEUECOMM,MYCOMM,&status);
             }
         }
      } else {
         MPI_Recv(&j,1,MPI_INT,0,QUEUECOMM,MYCOMM,&status);
      }
}

void ParParams::QueueNext() {
      int p[3]={0,0,0};
      int m;
      MPI_Cart_rank(MYCOMM,p,&m);  
      if (m!=me) {
         int i;
         MPI_Send(&i,1,MPI_INT,m,QUEUECOMM,MYCOMM);
      }
}

//
// Global Operations Min/Max/Sum using Allreduce Operation
//

NS_REAL ParParams::CommMin(NS_REAL lmax) {
    int oldtimer=timer->Start(12,this);
    NS_REAL llmax=lmax;
    MPI_Allreduce(&lmax,&llmax,1,MPI_MYNS_REAL,MPI_MIN,MYCOMM);
    timer->Stop(oldtimer,this);
    return llmax;
}

NS_REAL ParParams::CommMax(NS_REAL lmax) {
    int oldtimer=timer->Start(13,this);
    NS_REAL llmax=lmax;
    MPI_Allreduce(&lmax,&llmax,1,MPI_MYNS_REAL,MPI_MAX,MYCOMM);
    timer->Stop(oldtimer,this);
    return llmax;
}

NS_REAL ParParams::CommSum(NS_REAL lres) {
    int oldtimer=timer->Start(14,this);
    NS_REAL llres=lres;
    MPI_Allreduce(&lres,&llres,1,MPI_MYNS_REAL,MPI_SUM,MYCOMM);
    timer->Stop(oldtimer,this);
    return llres;
}

//
// Miscelaneous parallel operations 
//

// ParParams:
// initializes MPI and the message buffers, which are dynamically allocated in SetBounds

ParParams::ParParams(int* pargc,char*** pargv) {
    int i,j; 
    for(i=0;i<6;i++) for(j=0;j<8;j++) sbuf[i][j]=rbuf[i][j]=NULL;
    MPI_Init(pargc,pargv);
    MPI_Comm_rank(MPI_COMM_WORLD,&me);
}

// ShareInfo:
// parent has already read the most important paramters, such as number of gridpoints  
// number of chemicals to be computed,... and now sends this information to all 
// other processes  
// This information is required for the dynamic allocation of the 3D arrays 
// periodbound is required to initialize a torus proc-topology or a cartesian topology in SetBounds

void ParParams::ShareInfo(Scene &S,char* filename) {
    MPI_Bcast(S.gridp         ,3  ,MPI_INT ,0,MPI_COMM_WORLD) ;
    MPI_Bcast(filename        ,256,MPI_CHAR,0,MPI_COMM_WORLD) ;
    MPI_Bcast(&(S.nchem)      ,1  ,MPI_INT ,0,MPI_COMM_WORLD) ;
    MPI_Bcast(&(S.periodbound),1  ,MPI_INT ,0,MPI_COMM_WORLD) ;
    MPI_Bcast(&(S.CompChem)   ,1  ,MPI_INT ,0,MPI_COMM_WORLD) ;  
    MPI_Bcast(&(S.CompTemp)   ,1  ,MPI_INT ,0,MPI_COMM_WORLD) ;       
}

// ~ParParams:
// Finalizes MPI  

ParParams::~ParParams() {
    int i,j;
    for(i=0;i<6;i++) for(j=0;j<8;j++) {
        if(sbuf[i][j]) delete[] sbuf[i][j];
        if(rbuf[i][j]) delete[] rbuf[i][j];
    }
    MPI_Finalize(); 
}

//
// Coloured communication for color0-color7 SOR iteration 
// this allows to run most of the communication in background
//

#define ILP for(i=((iug+1+cx)&~1)-1-cx;i<=iog+1;i+=2)
#define JLP for(j=((jug+1+cy)&~1)-1-cy;j<=jog+1;j+=2)
#define KLP for(k=((kug+1+cz)&~1)-1-cz;k<=kog+1;k+=2)

void ParParams::CalcColourCount()  {
    int i,j,k,l,cx,cy,cz;
    for(i=0;i<8;i++) for(j=0;j<3;j++) colourcount[i][j]=0;
    for(l=0;l<8;l++) {
        cx=l & 1;cy=(l >> 1) & 1; cz=(l >> 2) & 1;
        JLP KLP colourcount[l][0]++;
        ILP KLP colourcount[l][1]++;
        ILP JLP colourcount[l][2]++;
    }
}

void ParParams::SendColour(Matrix<NS_REAL>&matrix, int cx, int cy, int cz) {
    int oldtimer=timer->Start(15,this);
    int i, j, k, z;
    int num=cx+2*cy+4*cz;

    if ((ts!=-1) && ((iug&1)!=cx)) {          // SOUTH
      z=0;
      JLP KLP sbuf[1][num][z++]=matrix[iug][j][k];
      assert(z==colourcount[num][0]);
      MPI_Isend(sbuf[1][num],z,MPI_MYNS_REAL,ts,PCOMMS+num,MYCOMM,&srequest[1][num]);
    }
    if ((tn!=-1) && ((iog&1)!=cx)) {           // NORTH
      z=0;
      JLP KLP sbuf[0][num][z++]=matrix[iog][j][k];
      assert(z==colourcount[num][0]);
      MPI_Isend(sbuf[0][num],z,MPI_MYNS_REAL,tn,PCOMMN+num,MYCOMM,&srequest[0][num]);
    }                                   
    if ((te!=-1) && ((jug&1)!=cy)) {          // EAST
      z=0;
      ILP KLP sbuf[3][num][z++]=matrix[i][jug][k];
      assert(z==colourcount[num][1]);
      MPI_Isend(sbuf[3][num],z,MPI_MYNS_REAL,te,PCOMME+num,MYCOMM,&srequest[3][num]);
    }
    if ((tw!=-1) && ((jog&1)!=cy)) {          // WEST
      z=0;
      ILP KLP sbuf[2][num][z++]=matrix[i][jog][k]; 
      assert(z==colourcount[num][1]);
      MPI_Isend(sbuf[2][num],z,MPI_MYNS_REAL,tw,PCOMMW+num,MYCOMM,&srequest[2][num]);
    }
    if ((tb!=-1) && ((kug&1)!=cz)) {          // BOTTOM
      z=0;
      ILP JLP sbuf[5][num][z++]=matrix[i][j][kug];
      assert(z==colourcount[num][2]);
      MPI_Isend(sbuf[5][num],z,MPI_MYNS_REAL,tb,PCOMMB+num,MYCOMM,&srequest[5][num]);
    }
    if ((tt!=-1) && ((kog&1)!=cz)) {          // TOP
      z=0;
      ILP JLP sbuf[4][num][z++]=matrix[i][j][kog];
      assert(z==colourcount[num][2]);   
      MPI_Isend(sbuf[4][num],z,MPI_MYNS_REAL,tt,PCOMMT+num,MYCOMM,&srequest[4][num]);
    }

    if ((tt!=-1) && ((kog&1)==cz)) {              // TOP      
       MPI_Irecv(rbuf[4][num],colourcount[num][2],MPI_MYNS_REAL,tt,PCOMMB+num,MYCOMM,&rrequest[4][num]);
    }
    if ((tb!=-1) && ((kug&1)==cz)) {              // BOTTOM
       MPI_Irecv(rbuf[5][num],colourcount[num][2],MPI_MYNS_REAL,tb,PCOMMT+num,MYCOMM,&rrequest[5][num]);
    }
    if ((tw!=-1) && ((jog&1)==cy)) {              // WEST
       MPI_Irecv(rbuf[2][num],colourcount[num][1],MPI_MYNS_REAL,tw,PCOMME+num,MYCOMM,&rrequest[2][num]);
    }
    if ((te!=-1) && ((jug&1)==cy)) {              // EAST
       MPI_Irecv(rbuf[3][num],colourcount[num][1],MPI_MYNS_REAL,te,PCOMMW+num,MYCOMM,&rrequest[3][num]);
    }
    if ((tn!=-1) && ((iog&1)==cx)) {              // NORTH
       MPI_Irecv(rbuf[0][num],colourcount[num][0],MPI_MYNS_REAL,tn,PCOMMS+num,MYCOMM,&rrequest[0][num]);
    } 
    if ((ts!=-1) && ((iug&1)==cx)) {              // SOUTH
       MPI_Irecv(rbuf[1][num],colourcount[num][0],MPI_MYNS_REAL,ts,PCOMMN+num,MYCOMM,&rrequest[1][num]);
    }
    timer->Stop(oldtimer,this);
}

#define rcb(buf,i,j,k) matrix[i][j][k]=buf[z++];

void ParParams::ReceiveColour(Matrix<NS_REAL>&matrix, int cx, int cy, int cz) {
    int oldtimer=timer->Start(16,this);
    int         i, j, k, z;
    MPI_Status  status;
    int         num=cx+2*cy+4*cz;

    if ((ts!=-1) && ((iug&1)!=cx)) {              // WAIT FOR SOUTH DATA
       MPI_Wait(&srequest[1][num],&status);
    }
    if ((tn!=-1) && ((iog&1)!=cx)) {              // NORTH
       MPI_Wait(&srequest[0][num],&status);
    }                                   
    if ((te!=-1) && ((jug&1)!=cy)) {              // EAST
       MPI_Wait(&srequest[3][num],&status);
    }
    if ((tw!=-1) && ((jog&1)!=cy)) {              // WEST
       MPI_Wait(&srequest[2][num],&status);
    }
    if ((tb!=-1) && ((kug&1)!=cz)) {              // BOTTOM
       MPI_Wait(&srequest[5][num],&status);
    }
    if ((tt!=-1) && ((kog&1)!=cz)) {              // TOP
       MPI_Wait(&srequest[4][num],&status);
    }

    if ((tt!=-1) && ((kog&1)==cz)) {              // WRITE IN TOP      
       MPI_Wait(&rrequest[4][num],&status);
       z=0;
       ILP JLP rcb(rbuf[4][num],i, j, kog+1);
     }
    if ((tb!=-1) && ((kug&1)==cz)) {              // BOTTOM
       MPI_Wait(&rrequest[5][num],&status);
       z=0;
       ILP JLP rcb(rbuf[5][num],i, j, kug-1);
    }
    if ((tw!=-1) && ((jog&1)==cy)) {              // WEST
       MPI_Wait(&rrequest[2][num],&status);
       z=0;
       ILP KLP rcb(rbuf[2][num],i, jog+1, k);
    }
    if ((te!=-1) && ((jug&1)==cy)) {              // EAST
       MPI_Wait(&rrequest[3][num],&status);
       z=0;
       ILP KLP rcb(rbuf[3][num],i, jug-1, k);
    }
    if ((tn!=-1) && ((iog&1)==cx)) {              // NORTH
       MPI_Wait(&rrequest[0][num],&status);
       z=0;
       JLP KLP rcb(rbuf[0][num],iog+1, j, k);
    }
    if ((ts!=-1) && ((iug&1)==cx)) {              // SOUTH
       MPI_Wait(&rrequest[1][num],&status);
       z=0;
       JLP KLP rcb(rbuf[1][num],iug-1, j, k);
    }

    timer->Stop(oldtimer,this);
}

//
//  Send / Receive for matrix boundary values  
//
void ParParams::CommMatrix(Matrix<NS_REAL>&matrix, int vi, int vj, int vk, int tag, int GC) {
    int oldtimer=timer->Start(17,this);

    int         i, j, k, z, gc, row;
    MPI_Status  status;

    if((tag==UCOMM) || (tag==VCOMM) || (tag==WCOMM) || (tag==TCOMM) || ((CHCOMM-tag)<=0)) row=GC; else row=1; 

    if(tn!=-1) {
      z=0;
      for (gc=0; gc<row; gc++)
	for (j=jug-1; j<jog+2; j++) for (k=kug-1; k<kog+2; k++)
	  sbuf[0][0][z++]=matrix[iog-vi-gc][j][k];
      MPI_Isend(sbuf[0][0],z,MPI_MYNS_REAL,tn,tag+0,MYCOMM,&srequest[0][0]); 
      MPI_Irecv(rbuf[0][0],z,MPI_MYNS_REAL,tn,tag+1,MYCOMM,&rrequest[0][0]);
    }
    if(ts!=-1) {
      z=0;
      for (gc=0; gc<row; gc++) 
	for (j=jug-1; j<jog+2; j++) for (k=kug-1; k<kog+2; k++)
	  sbuf[1][0][z++]=matrix[iug+gc][j][k];
      MPI_Isend(sbuf[1][0],z,MPI_MYNS_REAL,ts,tag+1,MYCOMM,&srequest[1][0]);
      MPI_Irecv(rbuf[1][0],z,MPI_MYNS_REAL,ts,tag+0,MYCOMM,&rrequest[1][0]);
    }
    if(tw!=-1) {
      z=0;
      for (gc=0; gc<row; gc++) 
	for (i=iug-1; i<iog+2; i++) for (k=kug-1; k<kog+2; k++)
	  sbuf[2][0][z++]=matrix[i][jog-vj-gc][k];
      MPI_Isend(sbuf[2][0],z,MPI_MYNS_REAL,tw,tag+2,MYCOMM,&srequest[2][0]);
      MPI_Irecv(rbuf[2][0],z,MPI_MYNS_REAL,tw,tag+3,MYCOMM,&rrequest[2][0]);
    }
    if(te!=-1) {
      z=0;
      for (gc=0; gc<row; gc++) 
	for (i=iug-1; i<iog+2; i++) for (k=kug-1; k<kog+2; k++)
	  sbuf[3][0][z++]=matrix[i][jug+gc][k];
      MPI_Isend(sbuf[3][0],z,MPI_MYNS_REAL,te,tag+3,MYCOMM,&srequest[3][0]);
      MPI_Irecv(rbuf[3][0],z,MPI_MYNS_REAL,te,tag+2,MYCOMM,&rrequest[3][0]);
    }
    if(tt!=-1) {
      z=0;
      for (gc=0; gc<row; gc++)
	for (i=iug-1; i<iog+2; i++) for (j=jug-1; j<jog+2; j++)
	  sbuf[4][0][z++]=matrix[i][j][kog-vk-gc];
      MPI_Isend(sbuf[4][0],z,MPI_MYNS_REAL,tt,tag+4,MYCOMM,&srequest[4][0]);
      MPI_Irecv(rbuf[4][0],z,MPI_MYNS_REAL,tt,tag+5,MYCOMM,&rrequest[4][0]);
    }
    if(tb!=-1) {
      z=0;
      for (gc=0; gc<row; gc++) 
	for (i=iug-1; i<iog+2; i++) for (j=jug-1; j<jog+2; j++)
	  sbuf[5][0][z++]=matrix[i][j][kug+gc];
      MPI_Isend(sbuf[5][0],z,MPI_MYNS_REAL,tb,tag+5,MYCOMM,&srequest[5][0]);
      MPI_Irecv(rbuf[5][0],z,MPI_MYNS_REAL,tb,tag+4,MYCOMM,&rrequest[5][0]);
    }
    if (tn!=-1) {
      MPI_Wait(&srequest[0][0],&status);
      MPI_Wait(&rrequest[0][0],&status);
      z=0;
      for (gc=0; gc<row; gc++)
	for (j=jug-1; j<jog+2; j++) for (k=kug-1; k<kog+2; k++)
	  matrix[iog+1+gc][j][k]=rbuf[0][0][z++];
    }
    if(ts!=-1) {
      MPI_Wait(&srequest[1][0],&status);
      MPI_Wait(&rrequest[1][0],&status);  
      z=0;
      for (gc=0; gc<row; gc++) 
	for (j=jug-1; j<jog+2; j++) for (k=kug-1; k<kog+2; k++)
	  matrix[iug-vi-1-gc][j][k]=rbuf[1][0][z++];
    }
    if(tw!=-1) {
      MPI_Wait(&srequest[2][0],&status);
      MPI_Wait(&rrequest[2][0],&status);  
      z=0;
      for (gc=0; gc<row; gc++) 
	for (i=iug-1; i<iog+2; i++) for (k=kug-1; k<kog+2; k++)
	  matrix[i][jog+1+gc][k]=rbuf[2][0][z++];
    }
    if(te!=-1) {
      MPI_Wait(&srequest[3][0],&status);
      MPI_Wait(&rrequest[3][0],&status);  
      z=0;
      for (gc=0; gc<row; gc++) 
	for (i=iug-1; i<iog+2; i++) for (k=kug-1; k<kog+2; k++)
	  matrix[i][jug-vj-1-gc][k]=rbuf[3][0][z++];
    }
    if(tt!=-1) {
      MPI_Wait(&srequest[4][0],&status);
      MPI_Wait(&rrequest[4][0],&status);  
      z=0;
      for (gc=0; gc<row; gc++) 
	for (i=iug-1; i<iog+2; i++) for (j=jug-1; j<jog+2; j++)
	  matrix[i][j][kog+1+gc]=rbuf[4][0][z++];
    }
    if(tb!=-1) {
      MPI_Wait(&srequest[5][0],&status);
      MPI_Wait(&rrequest[5][0],&status);
      z=0;
      for (gc=0; gc<row; gc++) 
	for (i=iug-1; i<iog+2; i++) for (j=jug-1; j<jog+2; j++)
	  matrix[i][j][kug-vk-1-gc]=rbuf[5][0][z++];
    }
    timer->Stop(oldtimer,this);
}

void ParParams::CommMatrixFlag(Matrix<unsigned int>&matrix, int vi, int vj, int vk, int tag, int GC) {
    int oldtimer=timer->Start(17,this);

    int         i, j, k, z, gc, row=GC;
    MPI_Status  status;

    if(tn!=-1) {
      z=0;
      for (gc=0; gc<row; gc++)
	for (j=jug-1; j<jog+2; j++) for (k=kug-1; k<kog+2; k++)
	  sbuf[0][0][z++]=(NS_REAL)matrix[iog-vi-gc][j][k];
      MPI_Isend(sbuf[0][0],z,MPI_MYNS_REAL,tn,tag+0,MYCOMM,&srequest[0][0]); 
      MPI_Irecv(rbuf[0][0],z,MPI_MYNS_REAL,tn,tag+1,MYCOMM,&rrequest[0][0]);
    }
    if(ts!=-1) {
      z=0;
      for (gc=0; gc<row; gc++) 
	for (j=jug-1; j<jog+2; j++) for (k=kug-1; k<kog+2; k++)
	  sbuf[1][0][z++]=(NS_REAL)matrix[iug+gc][j][k];
      MPI_Isend(sbuf[1][0],z,MPI_MYNS_REAL,ts,tag+1,MYCOMM,&srequest[1][0]);
      MPI_Irecv(rbuf[1][0],z,MPI_MYNS_REAL,ts,tag+0,MYCOMM,&rrequest[1][0]);
    }
    if(tw!=-1) {
      z=0;
      for (gc=0; gc<row; gc++) 
	for (i=iug-1; i<iog+2; i++) for (k=kug-1; k<kog+2; k++)
	  sbuf[2][0][z++]=(NS_REAL)matrix[i][jog-vj-gc][k];
      MPI_Isend(sbuf[2][0],z,MPI_MYNS_REAL,tw,tag+2,MYCOMM,&srequest[2][0]);
      MPI_Irecv(rbuf[2][0],z,MPI_MYNS_REAL,tw,tag+3,MYCOMM,&rrequest[2][0]);
    }
    if(te!=-1) {
      z=0;
      for (gc=0; gc<row; gc++) 
	for (i=iug-1; i<iog+2; i++) for (k=kug-1; k<kog+2; k++)
	  sbuf[3][0][z++]=(NS_REAL)matrix[i][jug+gc][k];
      MPI_Isend(sbuf[3][0],z,MPI_MYNS_REAL,te,tag+3,MYCOMM,&srequest[3][0]);
      MPI_Irecv(rbuf[3][0],z,MPI_MYNS_REAL,te,tag+2,MYCOMM,&rrequest[3][0]);
    }
    if(tt!=-1) {
      z=0;
      for (gc=0; gc<row; gc++)
	for (i=iug-1; i<iog+2; i++) for (j=jug-1; j<jog+2; j++)
	  sbuf[4][0][z++]=(NS_REAL)matrix[i][j][kog-vk-gc];
      MPI_Isend(sbuf[4][0],z,MPI_MYNS_REAL,tt,tag+4,MYCOMM,&srequest[4][0]);
      MPI_Irecv(rbuf[4][0],z,MPI_MYNS_REAL,tt,tag+5,MYCOMM,&rrequest[4][0]);
    }
    if(tb!=-1) {
      z=0;
      for (gc=0; gc<row; gc++) 
	for (i=iug-1; i<iog+2; i++) for (j=jug-1; j<jog+2; j++)
	  sbuf[5][0][z++]=(NS_REAL)matrix[i][j][kug+gc];
      MPI_Isend(sbuf[5][0],z,MPI_MYNS_REAL,tb,tag+5,MYCOMM,&srequest[5][0]);
      MPI_Irecv(rbuf[5][0],z,MPI_MYNS_REAL,tb,tag+4,MYCOMM,&rrequest[5][0]);
    }
    if (tn!=-1) {
      MPI_Wait(&srequest[0][0],&status);
      MPI_Wait(&rrequest[0][0],&status);
      z=0;
      for (gc=0; gc<row; gc++)
	for (j=jug-1; j<jog+2; j++) for (k=kug-1; k<kog+2; k++)
	  matrix[iog+1+gc][j][k]=(unsigned int)rbuf[0][0][z++];
    }  
    if(ts!=-1) {
      MPI_Wait(&srequest[1][0],&status);
      MPI_Wait(&rrequest[1][0],&status);  
      z=0;
      for (gc=0; gc<row; gc++) 
	for (j=jug-1; j<jog+2; j++) for (k=kug-1; k<kog+2; k++)
	  matrix[iug-vi-1-gc][j][k]=(unsigned int)rbuf[1][0][z++];
    }
    if(tw!=-1) {
      MPI_Wait(&srequest[2][0],&status);
      MPI_Wait(&rrequest[2][0],&status);  
      z=0;
      for (gc=0; gc<row; gc++) 
	for (i=iug-1; i<iog+2; i++) for (k=kug-1; k<kog+2; k++)
	  matrix[i][jog+1+gc][k]=(unsigned int)rbuf[2][0][z++];
    }
    if(te!=-1) {
      MPI_Wait(&srequest[3][0],&status);
      MPI_Wait(&rrequest[3][0],&status);  
      z=0;
      for (gc=0; gc<row; gc++) 
	for (i=iug-1; i<iog+2; i++) for (k=kug-1; k<kog+2; k++)
	  matrix[i][jug-vj-1-gc][k]=(unsigned int)rbuf[3][0][z++];
    }
    if(tt!=-1) {
      MPI_Wait(&srequest[4][0],&status);
      MPI_Wait(&rrequest[4][0],&status);  
      z=0;
      for (gc=0; gc<row; gc++) 
	for (i=iug-1; i<iog+2; i++) for (j=jug-1; j<jog+2; j++)
	  matrix[i][j][kog+1+gc]=(unsigned int)rbuf[4][0][z++];
    }
    if(tb!=-1) {
      MPI_Wait(&srequest[5][0],&status);
      MPI_Wait(&rrequest[5][0],&status);
      z=0;
      for (gc=0; gc<row; gc++) 
	for (i=iug-1; i<iog+2; i++) for (j=jug-1; j<jog+2; j++)
	  matrix[i][j][kug-vk-1-gc]=(unsigned int)rbuf[5][0][z++];
    }
    timer->Stop(oldtimer,this);
}


//
//  parallel domain decomosition
//
// Note, the current implementation uses the (known !) numbers of gridpoints in each
// coordinate direction for a more elaborate domain decomposition

// compute all factorizations of a positive number n into NumberFaktors factors

# define MAXNUMBERFACTORIZATIONS 1000

void factorize(int n,int NumberFactors,int **Factors,int *NumberFactorizations) {
  int i,j,nfalt ;
  
  if ((*NumberFactorizations) == MAXNUMBERFACTORIZATIONS) {
    printf ("too much factorizations\n") ;fflush(stdout); MPI_Finalize() ;exit(0) ;
  }
  
  if(NumberFactors==1) {Factors[*NumberFactorizations][0]=n;(*NumberFactorizations)++;return;}
  
  for (i=1 ;i<=n ;i++)  if(n%i==0) { 
    nfalt= *NumberFactorizations ;
    factorize(n/i,NumberFactors-1,Factors,NumberFactorizations) ;
    for (j=nfalt ;j< *NumberFactorizations; j++)
      Factors[j][NumberFactors-1] = i ;
  }
}

void ParParams::PrintInitInfo() {
  if(me==0) {
    printf("Number of started processes %d\n",procs);
  }   
  int nlen;
  char pname[2048];
  int coords[3];
  MPI_Get_processor_name(pname,&nlen);
  MPI_Cart_coords(MYCOMM,me,3,coords);
  printf("%d (%d/%d/%d) procs, I am %d (%d/%d/%d) (%s),\n neighbours N:%2d S:%2d W:%2d E:%2d T:%2d B:%2d.\n",
	 procs,Iproc(),Jproc(),Kproc(),me,coords[0],coords[1],coords[2],pname,tn,ts,tw,te,tt,tb);
}

void ParParams::SetBounds(Scene &S) {
    int  coords[3],period[3] ;
    int  i,j;
   
    int  NumberFactorizations=0 ;
    int  *Factors[MAXNUMBERFACTORIZATIONS]       ;
    int  mem[3*MAXNUMBERFACTORIZATIONS]          ;
  
    for (i=0 ;i<MAXNUMBERFACTORIZATIONS ;i++) Factors[i]=&mem[3*i] ;
    double best,cost                            ;

    MPI_Comm_size(MPI_COMM_WORLD,&procs);                       // Get total number of parallel processes

    factorize(procs,3,Factors,&NumberFactorizations);           //  Domain Decomposition

    best=1e+20;                                                 // search best decomposition with respect to
    for (i=0; i<NumberFactorizations; i++) {                    // the induced amount of communication
        cost=S.gridp[0]/Factors[i][0]*S.gridp[1]/Factors[i][1]+ // simple cost function
             S.gridp[0]/Factors[i][0]*S.gridp[2]/Factors[i][2] +
             S.gridp[1]/Factors[i][1]*S.gridp[2]/Factors[i][2] ;
        if(best >= cost) {
            best=cost ;
            pproc[0]=Factors[i][0] ;  
            pproc[1]=Factors[i][1] ;
            pproc[2]=Factors[i][2] ;
        } 
    }

    if((pproc[0]>S.gridp[0]) || (pproc[1]>S.gridp[1]) || (pproc[2]>S.gridp[2])){
      if(IsParent()){printf("\nNo domain decomposition possible !!!\n\n");
                     fflush(stdout);MPI_Finalize();exit(0);}
    }

    period[0]= (S.periodbound&1) > 0 ;                       // let MPI distribute the procs processes in 3D 
    period[1]= (S.periodbound&2) > 0 ;
    period[2]= (S.periodbound&4) > 0 ; 
    MPI_Cart_create(MPI_COMM_WORLD,3,pproc,period,TRUE,&MYCOMM);

    MPI_Comm_rank(MYCOMM,&me);                                  // who am i ?
    MPI_Cart_coords(MYCOMM,me,3,coords);                        // and in the 3D process topology ?

    if ((coords[0]>0) || period[0] ) {                          // get neighbours in north/south
        coords[0]=(coords[0]+pproc[0]-1)%pproc[0]; 
        MPI_Cart_rank(MYCOMM,coords,&ts)         ;
        coords[0]=(coords[0]+1         )%pproc[0];
    } else ts=-1 ; 

    if ((coords[0]<pproc[0]-1) || period[0] ) { 
        coords[0]=(coords[0]+1)%pproc[0]; 
        MPI_Cart_rank(MYCOMM,coords,&tn)         ;
        coords[0]=(coords[0]+pproc[0]-1)%pproc[0];
    } else tn=-1 ;     

    if ((coords[1]>0) || period[1] ) {                          // get neighbours in east/west
        coords[1]=(coords[1]+pproc[1]-1)%pproc[1]; 
        MPI_Cart_rank(MYCOMM,coords,&te)         ;
        coords[1]=(coords[1]+1         )%pproc[1];
    } else te=-1 ; 

    if ((coords[1]<pproc[1]-1) || period[1] ) { 
        coords[1]=(coords[1]+1)%pproc[1]; 
        MPI_Cart_rank(MYCOMM,coords,&tw)         ;
        coords[1]=(coords[1]+pproc[1]-1)%pproc[1];
    } else tw=-1 ;     

    if ((coords[2]>0) || period[2] ) {                          // get neighbours in top/bottom
        coords[2]=(coords[2]+pproc[2]-1)%pproc[2]; 
        MPI_Cart_rank(MYCOMM,coords,&tb)         ;
        coords[2]=(coords[2]+1         )%pproc[2];
    } else tb=-1 ; 

    if ((coords[2]<pproc[2]-1) || period[2] ) { 
        coords[2]=(coords[2]+1)%pproc[2]; 
        MPI_Cart_rank(MYCOMM,coords,&tt)         ;
        coords[2]=(coords[2]+pproc[2]-1)%pproc[2];
    } else tt=-1 ;     
   
    biug=(coords[0]*S.gridp[0])/Iproc()+1   ;                   // set local lower and upper indices 
    biog=(coords[0]*S.gridp[0]+S.gridp[0])/Iproc();             // for computational domain 
    bjug=(coords[1]*S.gridp[1])/Jproc()+1   ;
    bjog=(coords[1]*S.gridp[1]+S.gridp[1])/Jproc();
    bkug=(coords[2]*S.gridp[2])/Kproc()+1   ;
    bkog=(coords[2]*S.gridp[2]+S.gridp[2])/Kproc();

    for(i=0;i<6;i++) for(j=0;j<8;j++) {                         // create send/receive buffers
        if(rbuf[i][j]) delete rbuf[i][j];
        if(sbuf[i][j]) delete sbuf[i][j];
    }

    for(j=0;j<8;j++) {
        rbuf[0][j]=new NS_REAL[4*(bjog-bjug+3)*(bkog-bkug+3)];
        rbuf[1][j]=new NS_REAL[4*(bjog-bjug+3)*(bkog-bkug+3)];
        rbuf[2][j]=new NS_REAL[4*(biog-biug+3)*(bkog-bkug+3)];
        rbuf[3][j]=new NS_REAL[4*(biog-biug+3)*(bkog-bkug+3)];
        rbuf[4][j]=new NS_REAL[4*(biog-biug+3)*(bjog-bjug+3)];
        rbuf[5][j]=new NS_REAL[4*(biog-biug+3)*(bjog-bjug+3)];
        sbuf[0][j]=new NS_REAL[4*(bjog-bjug+3)*(bkog-bkug+3)];
        sbuf[1][j]=new NS_REAL[4*(bjog-bjug+3)*(bkog-bkug+3)];
        sbuf[2][j]=new NS_REAL[4*(biog-biug+3)*(bkog-bkug+3)];
        sbuf[3][j]=new NS_REAL[4*(biog-biug+3)*(bkog-bkug+3)];
        sbuf[4][j]=new NS_REAL[4*(biog-biug+3)*(bjog-bjug+3)];
        sbuf[5][j]=new NS_REAL[4*(biog-biug+3)*(bjog-bjug+3)];
    }
    
}


int ParParams::writevect(FILE* f,int ip,int jp,int* kdim,const ParParams* par,const Matrix<NS_REAL>& mat,const int* dims) {
    NS_REAL* buffer=new NS_REAL[kdim[1]-kdim[0]+1];
    int dims1[4]={ip,jp,kdim[0],kdim[1]};
    MPI_Bcast(dims1,4,MPI_INT,0,par->MYCOMM);
    int matrixdims[2]={1,0};
    if(ip>=dims[0] && ip<=dims[1] && jp>=dims[2] && jp<=dims[3]) memcpy(matrixdims,dims+4,sizeof(int)*2);
    int* allmatrixdims=new int[2*par->procs];
    int* recvcount=new int[par->procs];
    int* recvdispl=new int[par->procs];
    MPI_Gather(matrixdims,2,MPI_INT,allmatrixdims,2,MPI_INT,0,par->MYCOMM);
    int i;

//    printf("(%d %d %d-%d [%x]) ",ip,jp,kdim[0],kdim[1],ftell(f)); 
//    for(i=0;i<par->procs;i++) printf("%d %d  ",allmatrixdims[2*i],allmatrixdims[2*i+1]);
//    printf("\n");

    for(i=0;i<par->procs;i++) {
        recvcount[i]=allmatrixdims[2*i+1]-allmatrixdims[2*i]+1;
        recvdispl[i]=allmatrixdims[2*i]-kdim[0];
//        printf("%d %d  ",recvcount[i],recvdispl[i]);
    }
//    printf("\n");
   
    int count=matrixdims[1]-matrixdims[0]+1;
//    printf("%d: %d\n",par->me,count);
    MPI_Gatherv((count!=0)?&(mat[dims1[0]][dims1[1]][matrixdims[0]]):0,
                count,MPI_MYNS_REAL,buffer,recvcount,recvdispl,MPI_MYNS_REAL,0,par->MYCOMM);

    // printf ("fseek\n") ;
    // fseek(f,0,SEEK_CUR);
    if(fwrite(buffer,kdim[1]-kdim[0]+1,sizeof(NS_REAL),f)!=sizeof(NS_REAL)) return FALSE;
    // fseek(f,0,SEEK_CUR);
    
    delete[] recvdispl;
    delete[] recvcount;
    delete[] allmatrixdims;
    delete[] buffer;
    return TRUE;    
}

void ParParams::ChildWriteMatrix(Matrix<NS_REAL>& mat,const int* dims) {
    int dims1[4];

    while(TRUE) {
        MPI_Bcast(dims1,4,MPI_INT,0,MYCOMM);
        if(dims1[0]==dims1[1] && dims1[1]==dims1[2] && dims1[2]==dims1[3] && dims1[3]==-1) break;
        int matrixdims[2]={1,0};
        if(dims1[0]>=dims[0] && dims1[0]<=dims[1] && dims1[1]>=dims[2] && dims1[1]<=dims[3]) 
          memcpy(matrixdims,dims+4,sizeof(int)*2);
        MPI_Gather(matrixdims,2,MPI_INT,0,0,MPI_INT,0,MYCOMM);
        int count=matrixdims[1]-matrixdims[0]+1;
        MPI_Gatherv((count!=0)?&(mat[dims1[0]][dims1[1]][matrixdims[0]]):0,
                    count,MPI_MYNS_REAL,0,0,0,MPI_MYNS_REAL,0,MYCOMM);
    } 
}

void ParParams::EndWriteMatrix() {
    int dims[4]={-1,-1,-1,-1};
    MPI_Bcast(dims,4,MPI_INT,0,MYCOMM);
}


// --------------------------------------------------------------------------------------------------------------------------
#else // --------------------------------------- NOPARALLEL -----------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------
#include "typen.h"
#include "matrix.h"
#include "parallel.h"
#include <string.h>
#include <stdlib.h>

//
// Queue-Operations for quasi sequential operations
//
// Each process gets the ID of the (0,0,0) in the cartesian proc-topology
// then, in a loop process (0,0,0) sends to all others prosesses a message,
//                                 but continues in the loop only after it has received an answer
//                                 from actual addressee   

// the other processes wait to receive their message and then leave Queue. 
// The answer message is send by running QueueNext       

void ParParams::Queue() {

}

void ParParams::QueueNext() {

}


//
// Global Operations Min/Max/Sum using Allreduce Operation
//

NS_REAL ParParams::CommMin(NS_REAL lmax) {
    return lmax;
}

NS_REAL ParParams::CommMax(NS_REAL lmax) {
    return lmax;
}

NS_REAL ParParams::CommSum(NS_REAL lres) {
    return lres;
}


//
// Miscelaneous parallel operations 
//

// ParParams:
// initializes MPI and the message buffers, which are dynamically allocated in SetBounds

ParParams::ParParams(int* pargc,char*** pargv) {
    int i,j; 
    for(i=0;i<6;i++) for(j=0;j<8;j++) sbuf[i][j]=rbuf[i][j]=NULL;
    ts=tn=te=tw=tt=tb=-1;
    me=0;
}

// ShareInfo:
// parent has already read the most important paramters, such as number of gridpoints  
// number of chemicals to be computed,... and now sends this information to all 
// other processes  
// This information is required for the dynamic allocation of the 3D arrays 
// periodbound is required to initialize a torus proc-topology or a cartesian topology in SetBounds

void ParParams::ShareInfo(Scene &S,char* filename) {

}

// ~ParParams:
// Finalizes MPI  

ParParams::~ParParams() {
    int i,j;
    for(i=0;i<6;i++) for(j=0;j<8;j++) {
        if(sbuf[i][j]) delete[] sbuf[i][j];
        if(rbuf[i][j]) delete[] rbuf[i][j];
    }
}

//
// Coloured communication for color0-color7 SOR iteration 
// this allows to run most of the communication in background
//

#define ILP for(i=((iug+1+cx)&~1)-1-cx;i<=iog+1;i+=2)
#define JLP for(j=((jug+1+cy)&~1)-1-cy;j<=jog+1;j+=2)
#define KLP for(k=((kug+1+cz)&~1)-1-cz;k<=kog+1;k+=2)

void ParParams::CalcColourCount()  {
    int i,j,k,l,cx,cy,cz;
    for(i=0;i<8;i++) for(j=0;j<3;j++) colourcount[i][j]=0;
    for(l=0;l<8;l++) {
        cx=l & 1;cy=(l >> 1) & 1; cz=(l >> 2) & 1;
        JLP KLP colourcount[l][0]++;
        ILP KLP colourcount[l][1]++;
        ILP JLP colourcount[l][2]++;
    }
}

void ParParams::SendColour(Matrix<NS_REAL>&matrix, int cx, int cy, int cz) {
    int i, j, k, z;
    int num=cx+2*cy+4*cz;

    if ((ts!=-1) || (tn!=-1)) {
        if ((ts!=-1) && ((iug&1)!=cx)) {          // SOUTH
            z=0;
            JLP KLP sbuf[1][num][z++]=matrix[iug][j][k];
            assert(z==colourcount[num][0]);
        }
        if ((tn!=-1) && ((iog&1)!=cx)) {           // NORTH
            z=0;
            JLP KLP sbuf[0][num][z++]=matrix[iog][j][k];
            assert(z==colourcount[num][0]);
        }                                   
    }
    if ((te!=-1) || (tw!=-1)) {
        if ((te!=-1) && ((jug&1)!=cy)) {          // EAST
            z=0;
            ILP KLP sbuf[3][num][z++]=matrix[i][jug][k];
            assert(z==colourcount[num][1]);
        }
        if ((tw!=-1) && ((jog&1)!=cy)) {          // WEST
            z=0;
            ILP KLP sbuf[2][num][z++]=matrix[i][jog][k];
            assert(z==colourcount[num][1]);
        }
    }
    if ((tt!=-1)||(tb!=-1)) {
        if ((tb!=-1) && ((kug&1)!=cz)) {      // BOTTOM
            z=0;
            ILP JLP sbuf[5][num][z++]=matrix[i][j][kug];
            assert(z==colourcount[num][2]);
        }
        if ((tt!=-1) && ((kog&1)!=cz)) {      // TOP
            z=0;
            ILP JLP sbuf[4][num][z++]=matrix[i][j][kog];
            assert(z==colourcount[num][2]);   
        }
    }
}

#define rcb(buf,i,j,k) matrix[i][j][k]=buf[z++];

void ParParams::ReceiveColour(Matrix<NS_REAL>&matrix, int cx, int cy, int cz) {
    int oldtimer=timer->Start(16,this);
    int         i, j, k, z;
    int         num=cx+2*cy+4*cz;

    if ((tb!=-1) || (tt!=-1)) {
        if ((tt!=-1) && ((kog&1)==cz)) {              // TOP      
            z=0;
            ILP JLP rcb(sbuf[5][num],i, j, kog+1);
        }
        if ((tb!=-1) && ((kug&1)==cz)) {              // BOTTOM
            z=0;
            ILP JLP rcb(sbuf[4][num],i, j, kug-1);
        }
    }
    if ((tw!=-1) || (te!=-1)) {
        if ((tw!=-1) && ((jog&1)==cy)) {              // WEST
            z=0;
            ILP KLP rcb(sbuf[3][num],i, jog+1, k);
        }
        if ((te!=-1) && ((jug&1)==cy)) {              // EAST
            z=0;
            ILP KLP rcb(sbuf[2][num],i, jug-1, k);
        }
    }
    if ((ts!=-1) || (tn!=-1)) {
        if ((tn!=-1) && ((iog&1)==cx)) {              // NORTH
            z=0;
            JLP KLP rcb(sbuf[1][num],iog+1, j, k);
        }
        if ((ts!=-1) && ((iug&1)==cx)) {              // SOUTH
            z=0;
            JLP KLP rcb(sbuf[0][num],iug-1, j, k);
        }
    }
    timer->Stop(oldtimer,this);
}

//
//  Send / Receive for matrix boundary values  
//

void ParParams::CommMatrix(Matrix<NS_REAL>&matrix, int vi, int vj, int vk, int tag, int GC) {
    int         i, j, k,gc,row;

    if((tag==UCOMM) || (tag==VCOMM) || (tag==WCOMM) || (tag==TCOMM) || ((CHCOMM-tag)<=0)) row=GC; else row=1;

    if(tn!=-1) {
      for (gc=0; gc<row; gc++)        
        for (j=jug-1; j<jog+2; j++) for (k=kug-1; k<kog+2; k++)
	  matrix[iug-vi-1-gc][j][k]=matrix[iog-vi-gc][j][k];     
    }
    if(ts!=-1) {
      for (gc=0; gc<row; gc++) 
        for (j=jug-1; j<jog+2; j++) for (k=kug-1; k<kog+2; k++)
	  matrix[iog+1+gc][j][k]=matrix[iug+gc][j][k];
    }
    if(tw!=-1) {
      for (gc=0; gc<row; gc++)
        for (i=iug-1; i<iog+2; i++) for (k=kug-1; k<kog+2; k++)
	  matrix[i][jug-vj-1-gc][k]=matrix[i][jog-vj-gc][k];
    }
    if(te!=-1) {
      for (gc=0; gc<row; gc++)      
        for (i=iug-1; i<iog+2; i++) for (k=kug-1; k<kog+2; k++)
	  matrix[i][jog+1+gc][k]=matrix[i][jug+gc][k];
    }
    if(tt!=-1) {
      for (gc=0; gc<row; gc++) 
        for (i=iug-1; i<iog+2; i++) for (j=jug-1; j<jog+2; j++)
	  matrix[i][j][kug-vk-1-gc]=matrix[i][j][kog-vk-gc];
    }
    if(tb!=-1) {
      for (gc=0; gc<row; gc++)
        for (i=iug-1; i<iog+2; i++) for (j=jug-1; j<jog+2; j++)
	  matrix[i][j][kog+1+gc]=matrix[i][j][kug+gc];
    }
    
}


void ParParams::CommMatrixFlag(Matrix<unsigned int>&matrix, int vi, int vj, int vk, int tag, int GC) {
    int         i, j, k,gc,row=GC;
    
    if(tn!=-1) {
      for (gc=0; gc<row; gc++)        
        for (j=jug-1; j<jog+2; j++) for (k=kug-1; k<kog+2; k++)
	  matrix[iug-vi-1-gc][j][k]=matrix[iog-vi-gc][j][k];     
    }
    if(ts!=-1) {
      for (gc=0; gc<row; gc++) 
        for (j=jug-1; j<jog+2; j++) for (k=kug-1; k<kog+2; k++)
	  matrix[iog+1+gc][j][k]=matrix[iug+gc][j][k];
    }
    if(tw!=-1) {
      for (gc=0; gc<row; gc++)
        for (i=iug-1; i<iog+2; i++) for (k=kug-1; k<kog+2; k++)
	  matrix[i][jug-vj-1-gc][k]=matrix[i][jog-vj-gc][k];
    }
    if(te!=-1) {
      for (gc=0; gc<row; gc++)      
        for (i=iug-1; i<iog+2; i++) for (k=kug-1; k<kog+2; k++)
	  matrix[i][jog+1+gc][k]=matrix[i][jug+gc][k];
    }
    if(tt!=-1) {
      for (gc=0; gc<row; gc++) 
        for (i=iug-1; i<iog+2; i++) for (j=jug-1; j<jog+2; j++)
	  matrix[i][j][kug-vk-1-gc]=matrix[i][j][kog-vk-gc];
    }
    if(tb!=-1) {
        for (gc=0; gc<row; gc++)
        for (i=iug-1; i<iog+2; i++) for (j=jug-1; j<jog+2; j++)
	  matrix[i][j][kog+1+gc]=matrix[i][j][kug+gc];
    }
    
}  


void ParParams::SetBounds(Scene &S) {
    int i,j;
    procs=1;
    pproc[0]=pproc[1]=pproc[2]=1;
    me=0;

    biug=1;
    biog=S.gridp[0];
    bjug=1;
    bjog=S.gridp[1];
    bkug=1;
    bkog=S.gridp[2];
    if(S.periodbound&1) {ts=tn=0;}
    if(S.periodbound&2) {tw=te=0;}
    if(S.periodbound&4) {tt=tb=0;}
    for(i=0;i<6;i++) for(j=0;j<8;j++) {   // create send/receive buffers
        if(rbuf[i][j]) delete[] rbuf[i][j];
        if(sbuf[i][j]) delete[] sbuf[i][j];
    }
    for(j=0;j<8;j++) {                    // only needed for pressure calculation
        rbuf[0][j]=new NS_REAL[(bjog-bjug+3)*(bkog-bkug+3)];
        rbuf[1][j]=new NS_REAL[(bjog-bjug+3)*(bkog-bkug+3)];
        rbuf[2][j]=new NS_REAL[(biog-biug+3)*(bkog-bkug+3)];
        rbuf[3][j]=new NS_REAL[(biog-biug+3)*(bkog-bkug+3)];
        rbuf[4][j]=new NS_REAL[(biog-biug+3)*(bjog-bjug+3)];
        rbuf[5][j]=new NS_REAL[(biog-biug+3)*(bjog-bjug+3)];
        sbuf[0][j]=new NS_REAL[(bjog-bjug+3)*(bkog-bkug+3)];
        sbuf[1][j]=new NS_REAL[(bjog-bjug+3)*(bkog-bkug+3)];
        sbuf[2][j]=new NS_REAL[(biog-biug+3)*(bkog-bkug+3)];
        sbuf[3][j]=new NS_REAL[(biog-biug+3)*(bkog-bkug+3)];
        sbuf[4][j]=new NS_REAL[(biog-biug+3)*(bjog-bjug+3)];
        sbuf[5][j]=new NS_REAL[(biog-biug+3)*(bjog-bjug+3)];
    }
}

void ParParams::PrintInitInfo() {
    printf("\nRunning non-parallel.\n");
}

#endif

void Timer::Dump(const char** names) {
#ifdef TIMED
    double gtime=0;
    int i;
    printf("\n");
    for(i=1;i<100;i++) gtime+=times[i];
    for(i=1;i<100;i++) if(times[i]>0) printf("%7.4f %6.3f%c  %s\n",
            times[i],100*times[i]/gtime,37,names?names[i]?names[i]:"":"");
    printf("\n");
#endif
}
