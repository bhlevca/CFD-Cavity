/* NaSt3DGP - The Parallel 3D Navier-Stokes Solver
 * Copyright (C) 2003 Institute for Numerical Simulation
 *                    University of Bonn
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 */

#ifndef TYPEN_INCLUDED
#define TYPEN_INCLUDED
#include <math.h>
#include <stdio.h>
#include <assert.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#else
#error "config.h not found but needed for compilation. Run configure!"
#endif

#ifdef STDC_HEADERS
#include <float.h>
#else
#error "Need ANSI C headers. Sorry, given up."
#endif

#define NS_REAL       double
#define MAXNS_REAL    DBL_MAX
#define MINNS_REAL    DBL_MIN
#define MPI_MYNS_REAL MPI_DOUBLE
// NS_REAL data type

#ifndef TRUE
#define TRUE  1
#endif
#ifndef FALSE
#define FALSE 0
#endif

// for better portability, DBL_MAX from float.h is used
//#ifndef MAXDOUBLE
//#define MAXDOUBLE 1.79769313486231500e+308
//#endif

// for better portability, FLT_MAX from float.h is used
//#ifndef MAXFLOAT
//#define MAXFLOAT (float)3.402823466e38
//#endif

#define sqr(x)      ((x)*(x))
#define min3(x,y,z) ((((x) < (y)) && ((x) < (z)))?(x):(((y) < (z))?(y):(z)))
#define max(x,y)    ((x) > (y)?(x):(y))
#define min(x,y)    ((x) < (y)?(x):(y))
#define max3(x,y,z) (max(max((x),(y)),(z)))
#define ffabs(a) (((a)>0?(a):(-a)))

#ifdef NOPARALLEL
#ifdef TIMED
#undef TIMED
#endif
#endif

enum		Endian	{Big, little};

void		SwapEndian(char* x, int size, int n, unsigned swap);
int			CheckEndian(void);

enum {
  OBSTACLE  = 0x0001,        // =1;    Zelle in Hindernis ohne Fluid-Nachbarn 
  SLIP      = 0x0002,        // =2;    Rutschbedingung 
  INOUT     = 0x0004,        // =4;    Ein- bzw. Ausstroemzelle 
  NORTH     = 0x0008,        // =8;    noerdl. Zelle Fluid 
  SOUTH     = 0x0010,        // =16;   suedl.    "     "   
  EAST      = 0x0020,        // =32;   oestl.    "     "   
  WEST      = 0x0040,        // =64;   westl.    "     "   
  TOP       = 0x0080,        // =128;  obere     "     "   
  BOTTOM    = 0x0100,        // =256;  untere    "     "   
  TEMPB     = 0x0200,        // =512;  Temperatur-Randwerte berechnen
  CHEMB     = 0x0400,        // =1024; Chemie-Randwerte berechnen  
  INOUTTYP2 = 0x0800,        // =2048; two bits for type of INOUT condition
  INOUTTYP3 = 0x1000         // =4096;
};

enum ConvectiveType {DonorCell, QUICK, HLPA, SMART, VONOS};

class Scene;

#include "list.h"

//! Class for inflow boundary conditions
class Inflow {
public:
     //! cartesian coordinate indices (i,j,k) of the computational grid
    int        i,j,k;
    //! u,v,w are the velocity values in x,y and z direction, respectively
    NS_REAL    u,v,w;
    //! t,l are the temperature and level-set values
    NS_REAL    t,l;
    //! chem is a vector [1...n] containing all viscosity values of chemicals 
    NS_REAL*   chem;
    //! nchem is the total number of chemicals 
    int        nchem;
    void       SetCoord(int ii,int jj,int kk) ;
    int        WriteFile(FILE* f) ;
    int        ReadFile (FILE* f) ;
               Inflow() ;
               ~Inflow() ;
    void operator=(const Inflow& inflw) ;
    void       Print(char *head) ;
};

typedef Inflow INFLOW;

//! Class for outflow boundary conditions
class InOutRec {
public:
    //! cartesian coordinate indices (i,j,k) of the computational grid
    int        i,j,k;
    //! type defines the type of outflow boundary conditions (1=neumann zero, 2=time dependent convective) 
    int        type;
    //! param is an additional parameter to enforce the compatibility condition for outflow boundaries (0=OFF, 1=ON) 
    NS_REAL    param;
    void       SetCoord(int ii,int jj,int kk) ;
    int        WriteFile(FILE* f) ;
    int        ReadFile (FILE* f)  ;
               InOutRec() ;
    void operator=(const InOutRec& ior) ;
};
typedef InOutRec INOUTREC;

//! Includes some data operating routines for class InOutRec 
class InOutRecData:public InOutRec {
 public:
     NS_REAL    data[3];
     
     int     WriteFile(FILE* f) ;
     int     ReadFile (FILE* f) ;
     InOutRecData():InOutRec() {for(int i=0;i<3;i++) data[i]=0;}
     void operator=(const InOutRec& ior) {i=ior.i;j=ior.j;k=ior.k;type=ior.type;param=ior.param;}
     void operator=(const InOutRecData& ior) {i=ior.i;j=ior.j;k=ior.k;type=ior.type;
     param=ior.param;
     for(int i=0;i<3;i++) data[i]=ior.data[i];}
};
typedef InOutRecData INOUTRECDATA;

//! Recording of boundary cells 
class CellRec {
public:
    int     i,j,k;
    void    SetCoord(int ii,int jj,int kk) {i=ii;j=jj;k=kk;}
    int     WriteFile(FILE* f) ;
    int     ReadFile (FILE* f) ;
            CellRec() {i=j=k=0;}
    void operator=(const CellRec& cr) {i=cr.i;j=cr.j;k=cr.k;}        
};
typedef CellRec CELLREC;

//! Recording of special boundary cells, like cells which have more than one neigbouring fluid cell
class SpecCellRec {
public:
    int        i,j,k;
    NS_REAL    den;
    NS_REAL    gweight;
    void       SetCoord(int ii,int jj,int kk) {i=ii;j=jj;k=kk;}
    int        WriteFile(FILE* f) ;
    int        ReadFile (FILE* f) ;
               SpecCellRec() {i=j=k=0;den=gweight=0;}
    void operator=(const CellRec& cr) {i=cr.i;j=cr.j;k=cr.k;}        
};
typedef SpecCellRec SPECCELLREC;

class NavTokenizer;

//! Class Scene includes all variables and methods important for complete physical and geometrical scene description 
class Scene {
 public:
     //! version number
     char 	          version[7];       
     //! time
     NS_REAL              t;               
     //! number of gridpoints in each co-ordinate direction
     int                  gridp[3];           
     //! physical domain length in each co-ordinate direction  
     NS_REAL              dimension[3];
     //! number of time steps between output
     int                  prstep; 
     //! maximal number of iterations for solving the poisson equation
     int   	         itermax;    
     //! number of calculated timesteps
     int                  timesteps;
     //! timestepmethod needed for Adams-Bashfort time discretization, since first step is Euler
     int                  timestepmethod;    
     //! total number of obstacle cells
     int                  ObstacleCount;   
     //! single actual time step 
     NS_REAL              delt;
     //! single previous time step
     NS_REAL          delt_old;
     //! final time at which the calculation finishes 
     NS_REAL              Tfin;
     //! accuracy threshold for residual of the poisson iteration 
     NS_REAL              eps;
     //! relaxation paramameter for SOR, SSOR, Red-Black and 7color poisson solver scheme 
     NS_REAL omg;         
     //! donor-cell upwind parameter for velocity in momentum equations
     NS_REAL              alpha;
     //! donor-cell upwind parameter for velocity in transport equation
     NS_REAL alphatg;      
     //! CFL-number for diffusion  
     NS_REAL              tfdiff;
     //! CFL-number for convection
     NS_REAL              tfconv;      
     //! volume expansion coefficient for temperature calculation
     NS_REAL              beta;
     //! volume forces in x,y- and z-direction
     NS_REAL              g[3];
     //! kinematic viscosity (nu=1/re)
     NS_REAL              nu;
     //! reynolds number (re=1/nu)
     NS_REAL              re;
     //! prandtl number
     NS_REAL         prandtl;
     //! in/outfow boundaries tag
     unsigned int         iobd;        
     //! time diskretisation type (Euler or Adams-Bashfort or Runge-Kutta)
     unsigned int         TimeDis;            
     //! tag for dimensionless calculation  
     unsigned int         dimless;        
     //! order of spacial discretization for the reinitialisation scheme (1,2,3 or 5 are possible) 
     unsigned int         reinitOrder;        
     //! tag for output file to  include closure of level-set
     unsigned int         LevelSetClosure;    
     //! solver type for poisson equation (SOR, Red-Black or BiCGStab with Jacobi preconditioneer)
     unsigned int         Solver;             
     //! convective terms handling (VONOS, SMART, HLPA, QUICK or Donor Cell (DC))
     ConvectiveType       Convective;    
     //! number of ghostcells for convective terms 1 or 2 gc
     unsigned int         GC;                 
     //! maximal timestep size
     NS_REAL              deltmax;            
     //! reference temperature for Boussinesq-Approximation  
     NS_REAL              TempRef;     
     //NS_REAL            TempCold, TempHot; no longer supported
     //! initial velocity value in x-direction
     NS_REAL              ui;
     //! initial velocity value in y-direction;
     NS_REAL              vi;
     //! initial velocity value in z-direction;
     NS_REAL              wi;
     //! initial value for pressure
     NS_REAL		  pi;
     //! total surface-volume of outflow boundaries  
     NS_REAL              vol_of_iobc1_cells;
     //! total number of chemicals 
     int                  nchem;          
     //! vector [1-n] containing viscosity constants for chemicals
     NS_REAL*             chemc;    
     //! physical time step when writing the actual data in 'targetdir'
     NS_REAL              prtTime;      
     //! name of the 'targetdir' directory in which a sequence of output-files should be written
     char                 targetdir[2000];    

     //! tag for computing temperature 
     unsigned int         CompTemp;  
     //! tag for computing chemicals 
     unsigned int         CompChem;
          
     //! dimensionless froude number
     NS_REAL 		  froude;
     
     //! gridspacings for pressure cells  
     NS_REAL*             d[3];               
     //! (DX[i]+DX[i-1])/2 
     NS_REAL*             dm[3];              
     //! absolute coordinates
     NS_REAL*             kabs[3];            
     //! 3 pt stars second derivative for coordinates x,y,z; weights l(eft),m(iddle),r(ight)  
     NS_REAL*             ddstar[3][3];       
     //! 3 pt stars second derivative for coordinates x,y,z; weights l(eft),m(iddle),r(ight)  
     NS_REAL*             ddPstar[3][3];      
     //! 3 pt stars second derivative for coordinates x,y,z; weights l(eft),m(iddle),r(ight)  
     NS_REAL*             ddSstar[3][3];      
     //! ddiv[0][i][j]=d[i][j+1]/d[i][j], ddiv[1][i][j]=d[i][j-1]/d[i][j]
     NS_REAL*             ddiv[2][3];         
     //! periodic boundaries (000000...000ZYX bitwise)
     unsigned int         periodbound;        
     //! list of cells with inflow b.c.
     NS_List<Inflow>      InflowList;         
     //! list of cells with outflow b.c.
     NS_List<InOutRec>    InOutList;          
     
     //! default constructor
     Scene();
     //! default destructor
     ~Scene ();
     int               WriteFile (FILE * f);
     int               ReadFile (FILE * f);
     //! prints object info to out
     void              DumpInfo(FILE* out);    
     //! init kabs,...  
     void              Init();                 
     void              ParseDimension(NavTokenizer& t);
     void              ParseParameter(NavTokenizer& t);
     void              DeleteAll();
};

#define PI 3.14159265358979323846264338327950288

#define ISFLUID(flag)     (((flag)&OBSTACLE) == 0)
#define IFFLUID(flag)     if(ISFLUID(flag))
#define IFNOFLUID(flag)   if(((flag)&OBSTACLE)==OBSTACLE)
#define IFOBSTACLE(flag)  if(((flag)&OBSTACLE)==OBSTACLE)
#define IFBORDER(flag)    if(FREECELLS(flag)!=0)
#define OBSTACLETYPE(flag)((flag)&(SLIP|INOUT))
#define FREECELLS(flag)   ((flag)&(NORTH|SOUTH|WEST|EAST|TOP|BOTTOM)) 

#define INOUTTYP(flag)    (((flag)&(INOUTTYP2|INOUTTYP3))/INOUTTYP2)
#define IFINOUTTYP1(flag)   if((((flag)/INOUTTYP2)&3) ==0)
#define IFINOUTTYP2(flag)   if((((flag)/INOUTTYP2)&3) ==1)

#define DX S.d[0]
#define DY S.d[1]
#define DZ S.d[2]
#define DXM S.dm[0]
#define DYM S.dm[1]
#define DZM S.dm[2]

#define ILOOP       for(i=1;i<=S.gridp[0];i++)
#define JLOOP       for(j=1;j<=S.gridp[1];j++)
#define KLOOP       for(k=1;k<=S.gridp[2];k++)
#define IJKLOOP ILOOP JLOOP KLOOP
#define IALLLOOP    for(i=0;i<=S.gridp[0]+1;i++)
#define JALLLOOP    for(j=0;j<=S.gridp[1]+1;j++)
#define KALLLOOP    for(k=0;k<=S.gridp[2]+1;k++)
#define IJKALLLOOP IALLLOOP JALLLOOP KALLLOOP

#define INBORDER(i,j,k) (((i)<=S.gridp[0]) && ((i)>0) && ((j)<=S.gridp[1]) && ((j)>0) && ((k)<=S.gridp[2]) && ((k)>0))
#define INABORDER(i,j,k) (((i)<=S.gridp[0]+1) && ((i)>=0) && ((j)<=S.gridp[1]+1) && ((j)>=0) && ((k)<=S.gridp[2]+1) && ((k)>=0))

#define PILOOP      for(i=Par->iug; i<=Par->iog; i++)
#define PIALLLOOP   for(i=Par->iug-1; i<Par->iog+2; i++)
#define PJLOOP      for(j=Par->jug; j<=Par->jog; j++)
#define PJALLLOOP   for(j=Par->jug-1; j<Par->jog+2; j++)
#define PKLOOP      for(k=Par->kug; k<=Par->kog; k++)
#define PKALLLOOP   for(k=Par->kug-1; k<Par->kog+2; k++)
#define PIJKLOOP PILOOP PJLOOP PKLOOP
#define PIJKALLLOOP PIALLLOOP PJALLLOOP PKALLLOOP

#define NOBORDERX(i,j,k) (((i)>=(((Par->iug==1)&& !(S.periodbound&1))? 2:0))  &&  ((i)<=(((Par->iog==S.gridp[0])&& !(S.periodbound&1))? (S.gridp[0]-2):(S.gridp[0])))  &&  ((j)>1)  &&  ((j)<S.gridp[1])  &&  ((k)>1)  &&  ((k)<S.gridp[2]))

#define NOBORDERY(i,j,k) (((i)>1)  &&  ((i)<S.gridp[0])  &&  ((j)>=(((Par->jug==1)&& !(S.periodbound&2))? 2:0))  &&  ((j)<=(((Par->jog==S.gridp[1])&& !(S.periodbound&2))?(S.gridp[1]-2):(S.gridp[1])))  &&  ((k)>1)  &&  ((k)<S.gridp[2]))

#define NOBORDERZ(i,j,k) (((i)>1)  &&  ((i)<S.gridp[0])  &&  ((j)>1)  &&  ((j)<S.gridp[1])  &&  ((k)>=(((Par->kug==1)&& !(S.periodbound&4))? 2:0))  &&  ((k)<=(((Par->kog==S.gridp[2])&& !(S.periodbound&4))?(S.gridp[2]-2):(S.gridp[2]))))

#define NOBORDER(i,j,k)  (((i)>=(((Par->iug==1)&& !(S.periodbound&1))? 2:1))  &&  ((i)<=(((Par->iog==S.gridp[0])&& !(S.periodbound&1))?(S.gridp[0]-1):S.gridp[0]))  &&  ((j)>=(((Par->jug==1)&& !(S.periodbound&2))? 2:1))  &&  ((j)<=(((Par->jog==S.gridp[1])&& !(S.periodbound&2))?(S.gridp[1]-1):S.gridp[1]))  &&  ((k)>=(((Par->kug==1)&& !(S.periodbound&4))? 2:1))  &&  ((k)<=(((Par->kog==S.gridp[2])&& !(S.periodbound&4))?(S.gridp[2]-1):S.gridp[2])))

/*
#define PIBALLLOOP  for(i=((Par->iug==1)?0:Par->iug); i<=((Par->iog==S.gridp[0])?Par->iog+1:Par->iog); i++)
#define PJBALLLOOP  for(j=((Par->jug==1)?0:Par->jug); j<=((Par->jog==S.gridp[1])?Par->jog+1:Par->jog); j++)
#define PKBALLLOOP  for(k=((Par->kug==1)?0:Par->kug); k<=((Par->kog==S.gridp[2])?Par->kog+1:Par->kog); k++)
#define PIJKBALLLOOP PIBALLLOOP PJBALLLOOP PKBALLLOOP

#define PIBLOOP     for(i=((Par->iug==1)?1:Par->iug-1); i<=((Par->iog==S.gridp[0])?Par->iog-1:Par->iog); i++)
#define PJBLOOP     for(j=((Par->jug==1)?1:Par->jug-1); j<=((Par->jog==S.gridp[1])?Par->jog-1:Par->jog); j++)
#define PKBLOOP     for(k=((Par->kug==1)?1:Par->kug-1); k<=((Par->kog==S.gridp[2])?Par->kog-1:Par->kog); k++)
#define PIJKBLOOP PIBLOOP PJBLOOP PKBLOOP
*/

#define PINBORDER(i,j,k) (((i)>=Par->iug)&&((i)<=Par->iog)&&\
                ((j)>=Par->jug)&&((j)<=Par->jog)&&\
                ((k)>=Par->kug)&&((k)<=Par->kog))

#define PINALLBORDER(i,j,k) (((i)>(Par->iug-3))&&((i)<(Par->iog+2))&&\
                ((j)>(Par->jug-3))&&((j)<(Par->jog+2))&&\
                ((k)>(Par->kug-3))&&((k)<(Par->kog+2)))
#define PPINALLBORDER(i,j,k)    (((i)>(Par->iug-2))&&((i)<(Par->iog+2))&&\
                ((j)>(Par->jug-2))&&((j)<(Par->jog+2))&&\
                ((k)>(Par->kug-2))&&((k)<(Par->kog+2)))

#define IAL for(i=(((Par->iug)&~1)+1);i<=Par->iog;i+=2)
#define JAL for(j=(((Par->jug)&~1)+1);j<=Par->jog;j+=2)
#define KAL for(k=(((Par->kug)&~1)+1);k<=Par->kog;k+=2)
#define IBL for(i=((Par->iug+1)&~1);i<=Par->iog;i+=2)
#define JBL for(j=((Par->jug+1)&~1);j<=Par->jog;j+=2)
#define KBL for(k=((Par->kug+1)&~1);k<=Par->kog;k+=2)


#define IALB for(i=(((Par->iog+1)&~1)-1);i>=Par->iug;i-=2)
#define JALB for(j=(((Par->jog+1)&~1)-1);j>=Par->jug;j-=2)
#define KALB for(k=(((Par->kog+1)&~1)-1);k>=Par->kug;k-=2)
#define IBLB for(i=((Par->iog)&~1);i>=Par->iug;i-=2)
#define JBLB for(j=((Par->jog)&~1);j>=Par->jug;j-=2)
#define KBLB for(k=((Par->kog)&~1);k>=Par->kug;k-=2)


// Debugging macros
#define TP(a) printf("%s\n",a);
#define TPN(n) printf(" TPN %d\n",n);
#define NFILE char fn[256];sprintf(fn,"D:\\c\\s\\d.%d",Par->Me());FILE* debugfile=fopen(fn,"a");
#define TPF(i,j,k,e) fprintf(debugfile,"%3d %3d %3d %e\n",i,j,k,e);fflush(debugfile);
#define CFILE fclose(debugfile);
#define NFILEP char fn[256];sprintf(fn,"D:\\c\\s\\d.%d",me);FILE* debugfile=fopen(fn,"a");
#define PRINTMATRIX(mat,i) {int j,k;for(j=Par->jug-1;j<=Par->jog+1;j++) {printf("%2d ",j);for(k=Par->kug-1;k<=Par->kog+1;k++) printf("%-13.10g ",mat[i][j][k]);printf("\n");}printf("\n");}
#define PRINTFLAGENTRY(flag,U) printf("mat[%d][%d][%d]=%e (%s%s%s%s%s%s) %s%s\n",i,j,k,U[i][j][k],((flag[i][j][k]&NORTH)?"N":"."),    \
                        ((flag[i][j][k]&SOUTH)?"S":"."),((flag[i][j][k]&EAST)?"E":"."),((flag[i][j][k]&WEST)?"W":"."),((flag[i][j][k]&TOP)?"T":"."),((flag[i][j][k]&BOTTOM)?"B":"."),((flag[i][j][k]&TEMPB)?"t":"."),((flag[i][j][k]&CHEMB)?"c":"."));

#endif
