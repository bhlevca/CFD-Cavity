/* NaSt3DGP - The Parallel 3D Navier-Stokes Solver
 * Copyright (C) 2003 Institute for Numerical Simulation
 *                    University of Bonn
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 */

#ifndef NAVIER_H
#define NAVIER_H

#include "typen.h"
#include "matrix.h"
#include <limits.h>
#include <fstream>

// change the following definition if you have too much objects in your NAV file
typedef unsigned int   flagtype;
const   flagtype       flagtypemsb=UINT_MAX & ~(UINT_MAX >> 1);

enum Errors {ERR_OPEN,ERR_READ,ERR_WRITE,ERR_DEFAULT};

class Timer;
class Object;

//! This class contains all data arrays for physical values and mathematical computations 
class Navier {
protected:    
     //! U is the array of the velocity component in x-direction
     Matrix <NS_REAL>     U;
     //! V is the array of the velocity component in y-direction
     Matrix <NS_REAL>     V;
     //! W is the array of the velocity component in z-direction
     Matrix <NS_REAL>     W; 
     //! P is the array for pressure
     Matrix <NS_REAL>     P;
     //! T is the array for temperature
     Matrix <NS_REAL>     T;
     //! Tmp1 is a temporary array needed for time discretization
     Matrix <NS_REAL>    Tmp1;
     //! Tmp2 is a temporary array needed for time discretization
     Matrix <NS_REAL>    Tmp2;
     //! Tmp3 is a temporary array needed for time discretization
     Matrix <NS_REAL>    Tmp3;
     //! Temporary array needed for the BiCGSTAB solver scheme  
     Matrix <NS_REAL>    BiCG_r0, BiCG_rj, BiCG_pj, BiCG_vj, BiCG_sj, BiCG_tj, MatBuf;
     //! KoeffPreCond array contains the coefficients for the jacobi preconditioner included in the BiCGStab scheme 
     Matrix <NS_REAL>    KoeffPreCond;
     //! flag field, this array contains the geometry and tags all cells types (FLUID, OBSTACLE,...)  
     Matrix <flagtype>   flag;               
     //! CH is a vector [1...n] of arrays needed to calculate n additional chemicals  
     Matrix <NS_REAL>*   CH;                           
     //! S contains all physical and computational parameters, the so called scene description 
     Scene               S;                                // physical+computational parameters
     //! the name of the binary file
     char                binfile[2048];          
    
public:
                        Navier();
                       ~Navier();
    static void         Error(Errors ErrNo,...);

protected:
    int                 ReadFile(char* filename);
    int                 WriteFile(char* filename);
};

//! Class NavierSetup contains all methods used to convert the binary output file   
class NavierSetup:public Navier {

public:
 	static unsigned     swapEndian;

protected:
    char                scenefile[2048],outfile[2048],loadfile[2048];
    //! tag for what kind of readable files to write
    unsigned            output;           
    unsigned            output_mode;
    //! Change endian
    Endian              systemEndian;
    Endian              desiredEndian;
    //! tag for loading flow data from file(s)
    unsigned            loadmode;            
    //! tag for writing the velocity field
    unsigned            write_veloc_field;     
    unsigned            save_mem;
    //! tag for adding the scalar value to the velocity field (only for visualization with vtk)
    unsigned            addU,addV,addW,addP,addT,addL,addCF,addC[99]; 
    //! parse command line
    void                ParseArgs(int argc,char** argv);     
    
    int                 ConvertFile(char* filename,int(*writehead)(const char*,const char*,FILE**,const Scene& S,int pos),
                             int(*writeonefield)(FILE**,Matrix<NS_REAL>&,const Scene&,int,int),
                             int(*writevelocfield)(Matrix<NS_REAL>&,Matrix<NS_REAL>&,Matrix<NS_REAL>&,const Scene&,FILE**,int));
    int                 ConvertOneFile(FILE* f,const char* extension,int(*writehead)(const char*,const char*,FILE**,const Scene& S,int pos),
                             int(*writeonefield)(FILE**,Matrix<NS_REAL>&,const Scene&,int,int),int pos);
    int                 ConvertFileVTK(char* filename,int(*writehead)(const char*,const char*,FILE**,const Scene& S,int pos),
                             int(*writeonefield)(FILE**,Matrix<NS_REAL>&,const Scene&,int,int),
                             int(*writevelocfield)(Matrix<NS_REAL>&,Matrix<NS_REAL>&,Matrix<NS_REAL>&,const Scene&,FILE**,int));
    int                 ConvertOneFileVTK(FILE* f,const char* extension,int(*writehead)(const char*,const char*,FILE**,const Scene& S, int pos),
					  int(*writeonefield)(FILE**,Matrix<NS_REAL>&,const Scene&,int,int),int pos);
    int                 ConvertTecplotFile(char* filename, int(*writehead)(const char*,const char*,FILE**,const Scene& S,int pos));  
    void                ReadFieldsFromFile(FILE* f,int(*readhead)(int*,double**,int,const char*,FILE**),
                             int(*readnextslice)(Matrix<NS_REAL>&,FILE*,int,int,int));
    void                ReadFieldFromFile(FILE* f,int(*readhead)(int*,double**,int,const char*,FILE**),
                             int(*readnextslice)(Matrix<NS_REAL>&,FILE*,int,int,int),int which,int ipmode);


    static int          WriteExplorerHead(const char* filename,const char* extension,FILE** f,const Scene& S,int pos); 
    static int          WriteExplorerOneField(FILE** f,Matrix<NS_REAL>& mat,const Scene& S,int slicenum,int);
    static int          WriteExplorerVelocField(Matrix<NS_REAL>& U,Matrix<NS_REAL>& V,Matrix<NS_REAL>& W,const Scene& S,FILE** f,int slicenum);
    static int          ReadExplorerHead(int* dims,double** spaces,int which,const char* name,FILE** f);
    static int          ReadExplorerNextSlice(Matrix<NS_REAL>& mat,FILE* f,int dimi,int dimj,int slice);

    static int          WriteTecplotHead(const char* filename,const char* extension,FILE** f,const Scene& S,int pos); 
    static int          WriteTecplotHead2(const char* filename,const char* extension,FILE** f,const Scene& S,int pos); 
    static int          WriteTecplotOneField(FILE** f,Matrix<NS_REAL>& mat,const Scene& S,int slicenum,int joker);
    static int          WriteTecplotVelocField(Matrix<NS_REAL>& U,Matrix<NS_REAL>& V,Matrix<NS_REAL>& W,const Scene& S,FILE** f,int slicenum); 
 
    static int          WriteVTKHead(const char* filename,const char* extension,FILE** f,const Scene& S,int pos); 
    static int          WriteVTKOneField(FILE** f,Matrix<NS_REAL>& mat,const Scene& S,int slicenum,int);
    static int          WriteVTKVelocField(Matrix<NS_REAL>& U,Matrix<NS_REAL>& V,Matrix<NS_REAL>& W,const Scene& S,FILE** f,int slicenum);

    static int          WriteVTKSPHead(const char* filename,const char* extension,FILE** f,const Scene& S,int pos); 
    static int          WriteVTKSPOneField(FILE** f,Matrix<NS_REAL>& mat,const Scene& S,int slicenum,int);

    static int          ReadBinHead(int* dims,double** spaces,int which,const char* name,FILE** f);
    static int          ReadBinNextSlice(Matrix<NS_REAL>& mat,FILE* f,int dimi,int dimj,int slice);

    static int          WriteTextHead(const char* filename,const char* extension,FILE** f,const Scene& S,int pos); 
    static int          WriteTextOneField(FILE** f,Matrix<NS_REAL>& mat,const Scene& S,int slicenum,int);
    static int          WriteTextVelocField(Matrix<NS_REAL>& U,Matrix<NS_REAL>& V,Matrix<NS_REAL>& W,const Scene& S,FILE** f,int slicenum);
    static int          ReadTextHead(int* dims,double** spaces,int which,const char* name,FILE** f);
    static int          ReadTextNextSlice(Matrix<NS_REAL>& mat,FILE* f,int dimi,int dimj,int slice);
    NS_REAL             GetValueAt(float x,float y,float z,Matrix<NS_REAL>& mat,const Scene& S);

    Object*             Parse(NavTokenizer& t);             // parse scene description (NAV) file 
    void                MarkBound(Matrix<flagtype>& flag);
    void                CheckFlag(Matrix<flagtype>& flag);  // check frag field 

public:
                        NavierSetup();
    void                DoIt(int argc,char* argv[]);
};

 
class ParParams;

/*! Class NavierCalc contains all kind of lists needed for 
    different boundary treatments for the velocities, pressure,
    temperature and chemicals. In addition there are all
    methods included which are important for the numerics.  
*/    
class NavierCalc:public Navier {
protected:
     //! RHS is the right hand side array of the poisson equation
    Matrix <NS_REAL>       RHS;
    //! some temporarly used arrays for velocity, temperature and level-set
    Matrix <NS_REAL>       F[2], G[2], H[2], T2, Ft, Fl;
    //! some temporarly used arrays for chemicals
    Matrix <NS_REAL>*      Fch;       
    //! parallelization parameters
    ParParams*             Par;                        
    //! list of inflow cells
    NS_List<INFLOW>        InflowList;              
    //! list of outflow cells
    NS_List<INOUTREC>      InOutList;                 
    NS_List<INOUTRECDATA>  InOut3List;      // list of inout 3 cells, not included in InOutList
    //! list of slip boundary condition cells
    NS_List<CELLREC>       SlipList;           
    //! list of noslip boundary condition cells
    NS_List<CELLREC>       NoslipList;                 
    //! list of Temperature/Chemicals d/dn=0 boundary condition cells
    NS_List<SPECCELLREC>   dTempdnList,dChemdnList;    
    //! pressure dp/dn=0 lists for each color group and each single direction
    NS_List<CELLREC>       PdpdnList[8][6];        
    //! pressure dp/dn=0 lists for each color group and multiple directions
    NS_List<SPECCELLREC>   PSpecdpdnList[8];          
    //! pressure boundary cells for outflow boundary condition II
    NS_List<CELLREC>       PInOut2List[8][6];         
    //! pressure boundary cells for outflow boundary condition I
    NS_List<CELLREC>       PInOut1List[8][6];            

    //! number of iterations for the poisson solver
    int                    iterations, alliterations; 
    //! number of steps
    int                    n;
    //! defines amount of information printed to stdout
    int                    printinfo;           
    //! tag which defines if F,G,H arrays are from actual or previous time step
    int                    actual_fgh;

    //! stdout variables during calculation
    NS_REAL             RHSbalance, vol_domain;
    //! stdout variables during calculation
    NS_REAL             flow_in, flow_out, outflow_surface; 
    //! some variables needed for calculation of outflow boundary values
    NS_REAL             inflow_surface, mass_diff, corrected_mass_diff; 
    //! number of seconds to run program (=0 -> infinite loop)
    int                 time_to_run;            
  
    //! velocity correction with pressure gradient due to chorin projection method (see User's Guide)   
    void                AdapUVW();                
    //! calculate temporary velocity fields (in space) due to chorin projection method (see User's Guide)
    void                CompFGH();
    //! BiCGStab poisson solver method 
    NS_REAL             BiCGStab();
    //! calculate temporary velocity fields (in time) due to chorin projection method (see User's Guide)
    void                CompTUVWfromFGH(int timestepmethod); // computes \tilde u 
    //! compute right hand side of poisson equation
    void                CompRHS();
    //! some info outputs on stdout
    void                InfoOutputs(NS_REAL res);
    //! write binary data output 
    void                OutputData();
    //! calculate new time step
    NS_REAL             TimeStep();
    //! compute transport equation in space
    void                CompTG(Matrix<NS_REAL>& CH, NS_REAL alpha, NS_REAL diffconst);
    //! compute transport equation in time
    void                CompTT(Matrix<NS_REAL>& CH,Matrix<NS_REAL>& FT,NS_REAL alpha,NS_REAL diffconst, int n, int timestepmethod);
    //! compute poisson equation
    NS_REAL             CompPoisson(double starttime);
    //! compute poisson equation
    NS_REAL             Poisson();
    //! compute coefficients for jacobi precondidioner
    void                CalcKoeffPreCond();

    //! set color dependent pressure boundary conditions (done in every SOR iteration)
    void                SetPBorder(Matrix<NS_REAL>& B,int cg);
    //! set global pressure boundary conditions
    void                SetPALLBorder(Matrix <NS_REAL> &X);
    //! set local pressure boundary conditions
    void                SetPLocalBorder(Matrix <NS_REAL> &X, int i, int j, int k);    
    //! set velocity/temperature/chemical boundary conditions
    void                SetObstacleCond(int method);		  
    void                SetObstacleCondTilde(int method);
    //! set boundary conditions for transport equation
    void                SetObstacleForTG(Matrix<NS_REAL>& M, int n);
    void                OutflowVelocityCorrection(Matrix<NS_REAL>& U,Matrix<NS_REAL>& V,Matrix<NS_REAL>& W, int method);
    NS_REAL             CalcRes(NS_REAL ijkroot);              
    NS_REAL             CalcMeanVal(Matrix<NS_REAL>& mat);        // compute \|mat\|_{l^2}
    //! generate lists described above
    void                BuildLocalLists();                   
    int                 ReadScene(char* filename);
    int                 WriteFile(char* filename);
    int                 WriteMatrix(FILE* f,Matrix<NS_REAL>& mat,const int* dims);
    void                ParseArgs(int argc,char** argv);

    //  the following functions are used locally
    //
    //! interpolate values by third order lagrange interpolation
    NS_REAL             LagrangeIntp(NS_REAL LL, NS_REAL L, NS_REAL R, NS_REAL RR, NS_REAL DD, NS_REAL D, NS_REAL P, NS_REAL PP);
    //! compute the preconditioning 
    NS_REAL             Precondition(NS_REAL P, int i, int j, int k);
    //! reverse computation of the preconditioning
    NS_REAL             PreconditionBack(NS_REAL P, int i, int j, int k);
    //! matrix vector multiplication
    NS_REAL             MatVecMult(Matrix<NS_REAL>& A,int i, int j, int k, int RES);
    //! compute the mixed partial derivative
    NS_REAL             DUV(NS_REAL LL,NS_REAL L,NS_REAL M,NS_REAL R,NS_REAL RR,NS_REAL KL,NS_REAL KR,NS_REAL DD,NS_REAL DM,NS_REAL D,NS_REAL DP,NS_REAL PP,NS_REAL alpha,int border);
    //! compute the partial derivative
    NS_REAL             DUU(NS_REAL LL,NS_REAL L,NS_REAL M,NS_REAL R,NS_REAL RR,NS_REAL DD,NS_REAL D,NS_REAL DP,NS_REAL PP,NS_REAL alpha, int border);
    NS_REAL             DL(NS_REAL LLL, NS_REAL LL, NS_REAL L, NS_REAL M, NS_REAL R, NS_REAL RR, NS_REAL RRR, 
			   NS_REAL DDD, NS_REAL DD, NS_REAL D, NS_REAL DP, NS_REAL PP, NS_REAL PPP, 
			   NS_REAL sign, int reinit);  
    void                GENSPECENTRY(int i,int j,int k,NS_List<SPECCELLREC>& list);
    void                GENENTRY(int i,int j,int k,NS_List<CELLREC>& list);
    void                GENPENTRY(int i,int j,int k,int cg,int listnumber);
    void                PSOR(int i,int j,int k,NS_REAL omg1,NS_REAL lomg);
    NS_REAL             STAGGVALUENOSLIP(Matrix<NS_REAL>& M,int DIR,int i,int j,int k,NS_REAL value);    
    NS_REAL             STAGGVALUESLIP(Matrix<NS_REAL>& M,int DIR,int i,int j,int k);    

    void                SetHomogeneousNeumannTangential(Matrix<NS_REAL>& U,Matrix<NS_REAL>& V,Matrix<NS_REAL>& W,
                                                        int i,int j,int k) ;
    void                SetDirichletTangential(Matrix<NS_REAL>& U,Matrix<NS_REAL>& V,Matrix<NS_REAL>& W,
					       int i,int j,int k,NS_REAL value[3]);
    void                SetHomogeneousNeumannNormal(Matrix<NS_REAL>& U,Matrix<NS_REAL>& V,Matrix<NS_REAL>& W,
						    int i,int j,int k) ;
    void                SetDirichletNormal(Matrix<NS_REAL>& U,Matrix<NS_REAL>& V,Matrix<NS_REAL>& W,
                                           int i,int j,int k,NS_REAL u,NS_REAL v,NS_REAL w) ;
    void                SetVelocityExtrapolation(Matrix<NS_REAL>& U,Matrix<NS_REAL>& V,Matrix<NS_REAL>& W,
					      int i,int j,int k) ;
    void                SetVelocityCorrection(Matrix<NS_REAL>& U,Matrix<NS_REAL>& V,Matrix<NS_REAL>& W,
					      int i,int j,int k) ;
    Timer*              timer;
   
public:
    //! this method is the whole numeric engine for the computation of the simulation   
    void                DoIt();

                        NavierCalc(int* pargc,char*** pargv);
                       ~NavierCalc();
};


#endif









