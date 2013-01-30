/* NaSt3DGP - The Parallel 3D Navier-Stokes Solver
 * Copyright (C) 2003 Institute for Numerical Simulation
 *                    University of Bonn
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 */

#ifndef PARALLEL_INCLUDED
#define PARALLEL_INCLUDED
 #include <sys/types.h>
 #include <time.h>
#ifndef NOPARALLEL
 #include "mpi.h"
#endif
#ifdef WIN32
#include <time.h>
#endif

#include "typen.h"

class Timer;

//! various parameters and method for parallel computation
/*!
  This class provides various parameters and methods for parallel computation.
 */
class ParParams {
protected: 
     //! the local, rectangular, computational domain are the cells [biug,biog]x[bjug,bjog]x[bkug,bkog]                                       
    int         biug, biog, bjug, bjog, bkug, bkog;   
    //! number of processes
    int         procs    ;                
    //! The processes are distributed in 3D cartesian topology: pproc[0]*pproc[1]*pproc[2]=procs    
    int         pproc[3] ;                
    //! ID of 'my' task
    int         me;                       
    //! IDs of  neighbours in north,south,west,...; -1 means no neighbour
    int         tn, ts, tw, te, tb, tt;   

    //! buffer sizes
    int         colourcount[8][3];        
    //! send buffers for every color/direction
    NS_REAL*    sbuf[6][8];               
    //! receive buffers for every color/direction
    NS_REAL*    rbuf[6][8];
    //! pointer to timer class (for speedup measurements etc.)
    Timer*      timer;

#ifndef NOPARALLEL
    MPI_Comm    MYCOMM;                   // MPI specific stuff
    MPI_Request srequest[6][8]; 
    MPI_Request rrequest[6][8];  
#endif

public:
                ParParams(int* pargc,char*** pargv);
               ~ParParams();

    void        SetTimer(Timer* t) {timer=t;}
    //! prints initialization info to stdout
    void        PrintInitInfo();          
    //! return bounds of process domain
    int         Iug() {return biug;}   
    //! return bounds of process domain
    int         Iog() {return biog;}
    //! return bounds of process domain
    int         Jug() {return bjug;}
    //! return bounds of process domain
    int         Jog() {return bjog;}
    //! return bounds of process domain
    int         Kug() {return bkug;}
    //! return bounds of process domain
    int         Kog() {return bkog;}
    void        SetBounds(Scene& S);
    //! detect the parent process (rank 0)
    int         IsParent() {return(me==0);}

    //! communication of minimum value using Allreduce
    NS_REAL     CommMin (NS_REAL lmax);      
    //! communication of maximum value using Allreduce
    NS_REAL     CommMax (NS_REAL lmax);
    //! communication of sum value using Allreduce
    NS_REAL     CommSum (NS_REAL lres);
    //! communication routine (send) for 8-Color SSOR method
    void        SendColour (Matrix < NS_REAL > &matrix, int cx, int cy, int cz);
    //! communication routine (receive) for 8-Color SSOR method
    void        ReceiveColour (Matrix < NS_REAL > &matrix, int cx, int cy, int cz);
    void        CalcColourCount();  
    //! communication routine to exchange boundary Planes of Matrix for pressure
    /*!
      this method is used for communication of boundary planes in field of type Matrix
      \param tag describes type of communication: UCOMM(VCOMM, WCOMM) for velocity component u(v,w),
      PCOMM for pressure, 
      LCOMM for level set
     */
    void        CommMatrix (Matrix < NS_REAL > &matrix, int vi, int vj, int vk, int tag, int GC);
    //! communication routine to exchange boundary Planes of Matrix for flag field
    void        CommMatrixFlag (Matrix < unsigned int > &matrix, int vi, int vj, int vk, int tag, int GC);
    void        Queue();
    void        QueueNext();
    void        ShareInfo(Scene& S,char *filename);
    //! method to synchronize processes based on MPI_Barrier
    void        Barrier() {

#ifndef NOPARALLEL        
        MPI_Barrier(MYCOMM);
#endif    
            }
    static int  writevect(FILE* f,int ip,int jp,int* kdim,const ParParams* par,const Matrix<NS_REAL>&,const int* dims);
    void        ChildWriteMatrix(Matrix<NS_REAL>& mat,const int* dims);
    void        EndWriteMatrix();

    int         Me() {return me;}
    int         Iproc() {return pproc[0];}
    int         Jproc() {return pproc[1];}
    int         Kproc() {return pproc[2];}
    double      GetTime() {
#ifdef NOPARALLEL
#ifdef WIN32
        return (double)clock()/(double)CLOCKS_PER_SEC;   
#else
        return (double)time(NULL) ;
#endif
#else
        return MPI_Wtime();
#endif
                }  
};

//! simple timer class based on MPI_Wtime
/*!
  this class provides a facility to measure elapsed time in parallel calculations. this is useful
  for measuring speedup etc. For parallel code it is based on MPI_Wtime, for serial code it uses
  standard system interface functions.
  
 */
class Timer {
protected:
    double          times[100];
    int             actual;
    double          starttime;
public:
                    Timer() {for(int i=0;i<100;i++) times[i]=0;actual=0;}
    inline int      Start(int tmr,ParParams* p) {
#ifdef TIMED
                        int oldact=actual;
                        double atime=p->GetTime();
                        times[actual]+=atime-starttime;
                        starttime=atime;actual=tmr;
                        return oldact;
#else 
                        return 0;
#endif
                    }
    inline void     Stop(int oldtmr,ParParams* p) {
#ifdef TIMED
                        double atime=p->GetTime();
                        times[actual]+=atime-starttime;
                        starttime=atime;actual=oldtmr;
#endif
                    }
                        
    void            Dump(const char** names);        
};

#define iug Iug()
#define iog Iog()
#define jug Jug()
#define jog Jog()
#define kug Kug()
#define kog Kog()

#define UCOMM           0*6+1           // several communication IDs
#define VCOMM           1*6+1
#define WCOMM           2*6+1
#define PCOMM           3*6+1           // for simple pressure communication
#define TCOMM           4*6+1
#define LCOMM           5*6+1
#define flagCOMM        6*6+1
#define RESCOMM         7*6+1+6
#define MAXCOMM         7*6+1+7
#define QUEUECOMM       7*6+1+8      
#define PCOMMN          7*6+1+9+0*8
#define PCOMMS          7*6+1+9+1*8
#define PCOMMW          7*6+1+9+2*8
#define PCOMME          7*6+1+9+3*8
#define PCOMMT          7*6+1+9+4*8
#define PCOMMB          7*6+1+9+5*8
#define CHCOMM          7*6+1+12+6*8       
#define FILECOMM1       10000
#define FILECOMM2       10001

#endif
