/* NaSt3DGP - The Parallel 3D Navier-Stokes Solver
 * Copyright (C) 2003 Institute for Numerical Simulation
 *                    University of Bonn
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 */

#include "typen.h"
#include <cstdio>
#include <cstdlib>
//FIXME: use C++ strings instead?
#include <cstring>


Scene::Scene() {
  // FIXME: here it is assumed that the VERSION macro defined in config.h
  // is a 4-char string
  sprintf(version, "NV%s", VERSION);

    int i,j;
    t=0.0;
    timesteps=1;
    timestepmethod=0;
    gridp[0]=gridp[1]=gridp[2]=15;
    dimension[0]=dimension[1]=dimension[2]=1.0;
    prstep          =20;    // write data output after 20 time step
    itermax         =100;   // maximal number of iterations for poisson-solver
    ObstacleCount   =0;      
    Tfin            =1.0;
    delt            =0.0;
    delt_old        =0.0;
    deltmax         =1.0;
    eps             =0.001;
    omg             =1.7;
    alpha=alphatg   =1.0;
    tfdiff          =0.3;
    tfconv          =0.3;   
    beta            =1.e-4; // volume expansion coefficient for temperature
    g[0]=g[1]=g[2]  =0.0;   // volume forces
    re              =10.0;  // reynolds number
    prandtl         =1.0;  
    froude          =1.0;
    dimless         =1;
    TimeDis         =0;   
    prtTime         =0;
    targetdir[0]    ='.';
    // TempCold        =0.0; TempHot         =1.0;
    TempRef         =273.0;
    ui=vi=wi=pi     =0.0;
    nchem           =0;
    chemc           =0;
    vol_of_iobc1_cells=0;
    Solver = 6;
    Convective = VONOS; // VONOS scheme as default scheme
    GC = 2; // number of ghostcells needed
    CompTemp=CompChem=FALSE;
    d[0]=d[1]=d[2]=0;
    dm[0]=dm[1]=dm[2]=0;
    kabs[0]=kabs[1]=kabs[2]=0;
    for(i=0;i<3;i++) for(j=0;j<3;j++) 
        ddSstar[i][j]=ddPstar[i][j]=ddstar[i][j]=ddiv[0][i]=ddiv[1][i]=0; 
    periodbound=0;
    iobd=0;
}


Scene::~Scene() {
  DeleteAll();
}


void Scene::DeleteAll() {
    InflowList.Delete();
    InOutList.Delete();
    for(int i=0;i<3;i++) {
        if(d[i]) delete[] --(--d[i]);
        if(dm[i]) {dm[i]-=2;delete[] dm[i];}
        if(kabs[i]) delete[] kabs[i];
        for(int j=0;j<3;j++) {
            if(ddstar[i][j])  delete[] ddstar[i][j];
            if(ddSstar[i][j]) delete[] ddSstar[i][j];
            if(ddPstar[i][j]) delete[] ddPstar[i][j];
        }
        if(ddiv[0][i]) delete[] ddiv[0][i];
        if(ddiv[1][i]) delete[] ddiv[1][i];
    }
    if(chemc!=NULL) { // delete viscosities for chemicals
	delete[] chemc;
    }
}


int Scene::WriteFile(FILE* f) {
    int i;
    if(fwrite(this, sizeof(Scene), 1, f)!=1) return FALSE;
    if(chemc) if((int)fwrite(chemc,sizeof(NS_REAL),nchem,f)!=nchem) return FALSE; 
    for(i=0;i<3;i++) if(fwrite(d[i],gridp[i]+2,sizeof(NS_REAL),f)!=sizeof(NS_REAL)) return FALSE;
    if(!InflowList.WriteFile(f)) return FALSE;
    if(!InOutList.WriteFile(f)) return FALSE;
    fflush(f);
    return TRUE;
}


int Scene::ReadFile(FILE*f) {  // Reads total numbers of gridpoints, number of chemicals
    int i,j;                   // gridspacings, and Inflow/InOutList

    DeleteAll();
    if(fread(this, sizeof(Scene), 1, f)!=1) return FALSE;
    if(version[0]!='N' || version[1]!='V') {
        printf("\nError: This is not a valid NAV file.\n");exit(1);
    }
    char * helper_string = new char[7];
    char * read_version = new char[5];
    //FIXME: numerischer Vergleich der Versionsnummern besser?
    read_version[0] = version[2];
    read_version[1] = version[3];
    read_version[2] = version[4];
    read_version[3] = version[5];
    read_version[4] = version[6];
    sprintf(helper_string, "NV%s", VERSION);
    if(strcmp(helper_string, version) != 0){
      fprintf(stderr, "\nWARNING: This is not a NAV file created by NaSt3DGPF version %s ! Incompatibilities may occur!\n",
    	      VERSION, read_version);
    //exit(1);
    //FIXME: try alternative ReadFile method instead of terminating.
    }
    delete[] helper_string;
    delete[] read_version;
    InflowList.Init();
    InOutList.Init();
    for(i=0;i<3;i++) for(j=0;j<3;j++) ddSstar[i][j]=ddPstar[i][j]=ddstar[i][j]=ddiv[0][i]=ddiv[1][i]=0;
    for(i=0;i<3;i++) kabs[i]=dm[i]=0;

    chemc=new NS_REAL[nchem]; assert(chemc);

    if((int)fread(chemc,sizeof(NS_REAL),nchem,f)!=nchem) return FALSE;

    for(i=0;i<3;i++) {
	d[i]=new NS_REAL[gridp[i]+5];assert(d[i]);d[i]+=2;
	if(fread(d[i],gridp[i]+2,sizeof(NS_REAL),f)!=sizeof(NS_REAL)) return FALSE;
    }

    if(!InflowList.ReadFile(f)) return FALSE;
    if(!InOutList.ReadFile(f)) return FALSE;
    
    Init();
    return TRUE;
}


void Scene::DumpInfo(FILE* out) {                       // print some information to file
  fprintf(out,"Time: %g, finish at %g\n"                            
	  "Gridpoints: %d (%dx%dx%d), containing %d obstacle cells\n"
	  "   eps=%g   omg=%g   alpha=%g   alphatg=%g\n"  
	  "   tfdiff=%g   tfconv=%g\n",
	  (float)t,(float)Tfin,gridp[0]*gridp[1]*gridp[2],gridp[0],gridp[1],
          gridp[2],ObstacleCount,
          (float)eps,(float)omg,(float)alpha,(float)alphatg,
          (float)tfdiff,(float)tfconv);
  fprintf(out,"   reynolds=%g   froude=%g\n",(float)re,(float)froude);
  fprintf(out,"   gx=%g   gy=%g   gz=%g\n"
	  "temperature is%s calculated\n"       
	  "chemicals are%s calculated\n",
	  g[0],g[1],g[2], 
	  (CompTemp?"":" not"),(CompChem?"":" not"));
  if(CompTemp){
    fprintf(out,"Temperature list:\n"                         \
                "   TmpRef=%g   prandtl=%g   beta=%g\n",
	    (float)TempRef,(float)prandtl,(float)beta);
    fprintf(out,"\n");
  }
  if(CompChem){ 
    fprintf(out,"Chemical list:\n");
    for(int i=0;i<nchem;i++) {
      fprintf(out,"   chemc[%2d]=%g",i,(float)chemc[i]);
      if(i%4==3) fprintf(out,"\n");
    }
  fprintf(out,"\n");
  }

  fprintf(out,"\n");
  if(periodbound!=0){ 
    fprintf(out,"periodic bound in %s%s%s-direction.\n",(periodbound&1)?"x":"",
	    (periodbound&2)?"y":"",(periodbound&4)?"z":"");
  }

  fprintf(out,"            time handling:  ");
  switch (TimeDis){
  case 0: fprintf(out,"Euler 1st \n");break;
  case 1: fprintf(out,"Adams-Bashfort 2nd \n");break;
  default: fprintf(out,"ERROR: wrong time scheme defined!\n");break;
  }
    
  fprintf(out,"         iteration scheme:  ");
  switch (Solver){
  case 0: fprintf(out,"SOR\n"); break;
  case 1: fprintf(out,"SSOR (fw/bw)\n"); break;
  case 2: fprintf(out,"Red-Black\n"); break;
  case 3: fprintf(out,"8-color SOR\n"); break;
  case 4: fprintf(out,"8-color SSOR (fw/bw)\n"); break;
  case 5: fprintf(out,"BiCGStab\n"); break;
  case 6: fprintf(out,"BiCGStab with Jacobi PC\n"); break;
  default: fprintf(out,"ERROR: wrong solver defined!\n"); break;
  }

  fprintf(out,"convective terms handling:  ");
  switch (Convective){
  case DonorCell: fprintf(out,"Donor-Cell \n\n"); break;
  case QUICK:     fprintf(out,"QUICK \n\n");      break;
  case HLPA:      fprintf(out,"HLPA \n\n");       break;
  case SMART:     fprintf(out,"SMART \n\n");      break;
  case VONOS:     fprintf(out,"VONOS \n\n");      break;
  default: fprintf(out,"ERROR: wrong convective terms scheme defined!\n"); break;
  }
  fflush(out);
}                        


//
// initialize some precalculated values:
// kabs     absolute grid positions
// d        cell widths
// dm       middle-of-cell distances
// dstar    constant factors for discretization of first derivative
// ddstar   constant factors for discretization of second derivative
// ddiv     quotient of adjacent cell widths
//
// timestepping parameters are calculated
//

void Scene::Init() {
    int i,j,k;
  
    for(i=0;i<3;i++) if(d[i]==0) { // if d[i] undef., generate  aequidistant grid  UEBERFLUESSIG!!! -> Siehe parse.cxx
        NS_REAL del=dimension[i]/gridp[i];
        d[i]=new NS_REAL[gridp[i]+5];assert(d[i]);d[i]+=2;
        for(j=-2;j<=gridp[i]+2;j++) d[i][j]=del;
    }
    if(periodbound&1) {
      d[0][-2]=d[0][gridp[0]-2];
      d[0][-1]=d[0][gridp[0]-1];
      d[0][0]=d[0][gridp[0]];
      d[0][gridp[0]+1]=d[0][1];
      d[0][gridp[0]+2]=d[0][2];
    }else{
      d[0][-2]=d[0][1];
      d[0][-1]=d[0][1];
      d[0][0]=d[0][1];
      d[0][gridp[0]+1]=d[0][gridp[0]];
      d[0][gridp[0]+2]=d[0][gridp[0]];
    }
    if(periodbound&2) {
      d[1][-2]=d[1][gridp[1]-2];
      d[1][-1]=d[1][gridp[1]-1];
      d[1][0]=d[1][gridp[1]];
      d[1][gridp[1]+1]=d[1][1];
      d[1][gridp[1]+2]=d[1][2];
    }else{
      d[1][-2]=d[1][1];
      d[1][-1]=d[1][1];
      d[1][0]=d[1][1];
      d[1][gridp[1]+1]=d[1][gridp[1]];
      d[1][gridp[1]+2]=d[1][gridp[1]];
    }
    if(periodbound&4) {
      d[2][-2]=d[2][gridp[2]-2];
      d[2][-1]=d[2][gridp[2]-1];
      d[2][0]=d[2][gridp[2]];
      d[2][gridp[2]+1]=d[2][1];
      d[2][gridp[2]+2]=d[2][2];
    }else{
      d[2][-2]=d[2][1];
      d[2][-1]=d[2][1];
      d[2][0]=d[2][1];
      d[2][gridp[2]+1]=d[2][gridp[2]];
      d[2][gridp[2]+2]=d[2][gridp[2]];
    }

//      d[0][-1]=d[0][gridp[0]-1];
//      d[1][-1]=d[1][gridp[1]-1];
//      d[2][-1]=d[2][gridp[2]-1];
    
    for(i=0;i<3;i++) {
        if(kabs[i]) delete kabs[i];
        if(dm[i]) {dm[i]-=2; delete dm[i];}
        kabs[i]=new NS_REAL[gridp[i]+3];assert(kabs[i]);
        dm[i]=new NS_REAL[gridp[i]+5];assert(dm[i]); 
	dm[i]+=2;
        dm[i][-2]=dm[i][-1]=dm[i][0]=d[i][0];                             
        for(j=1;j<=gridp[i]+1;j++) dm[i][j]=0.5*(d[i][j]+d[i][j-1]);
	dm[i][gridp[i]+2]=dm[i][gridp[i]+1];
        kabs[i][0]=-d[i][0];                                         
        for(j=1;j<gridp[i]+3;j++) kabs[i][j]=kabs[i][j-1]+d[i][j-1]; 
    }
//     for(i=0;i<3;i++) for(j=-2;j<=gridp[i]+2;j++)
// 	printf("dm[%i][%i]=%g \n",i,j,d[i][j]);
    
    for(k=0;k<3;k++) for(j=0;j<3;j++) {
        if(ddPstar[j][k]) delete ddPstar[j][k];
        ddPstar[j][k]=new NS_REAL[gridp[j]+1];assert(ddPstar[j][k]);
        if(ddSstar[j][k]) delete ddSstar[j][k];
        ddSstar[j][k]=new NS_REAL[gridp[j]+1];assert(ddSstar[j][k]);
        if(ddstar[j][k]) delete ddstar[j][k];
        ddstar[j][k]=new NS_REAL[gridp[j]+1];assert(ddstar[j][k]);
    }
    for(j=0;j<3;j++) for(i=0;i<=gridp[j];i++) {
        NS_REAL h,H;
        h=0.5*(d[j][i-1]+d[j][i]);
        H=0.5*(d[j][i]+d[j][i+1]); 
/*
        ddSstar[j][0][i]=2.0/(h*(h+H));
        ddSstar[j][1][i]=-2.0/(h*H);
        ddSstar[j][2][i]=2.0/(H*(h+H));
*/      
        ddPstar[j][0][i]=1.0/(d[j][i]*dm[j][i]); 
        ddPstar[j][1][i]=-(1.0/dm[j][i+1]+1.0/dm[j][i])/d[j][i]; 
        ddPstar[j][2][i]=1.0/(d[j][i]*dm[j][i+1]);

        ddstar[j][0][i]=1.0/(dm[j][i+1]*d[j][i]);
        ddstar[j][1][i]=-(1/d[j][i]+1/d[j][i+1])/dm[j][i+1];
        ddstar[j][2][i]=1.0/(dm[j][i+1]*d[j][i+1]);

        ddSstar[j][0][i]=1/(d[j][i]*dm[j][i]) ;
        ddSstar[j][1][i]=-1/(d[j][i]*dm[j][i+1])-1/(d[j][i]*dm[j][i]) ;
        ddSstar[j][2][i]=1/(d[j][i]*dm[j][i+1]) ;
    }

    for(j=0;j<3;j++) for(i=0;i<3;i++) {
        ddSstar[j][i][0]=ddSstar[j][i][gridp[j]];
    }

    if((chemc==0) && (nchem>0)) {
        chemc=new NS_REAL[nchem];
        for(i=0;i<nchem;i++) chemc[i]=1.0;
    }
    for(i=0;i<3;i++) {
        if(ddiv[0][i]) delete ddiv[0][i];
        if(ddiv[1][i]) delete ddiv[1][i];
        ddiv[0][i]=new NS_REAL[gridp[i]+2];ddiv[1][i]=new NS_REAL[gridp[i]+2];
        for(j=0;j<=gridp[i];j++) ddiv[0][i][j]=-d[i][j]/d[i][j+1];
        for(j=1;j<=gridp[i]+1;j++) ddiv[1][i][j]=-d[i][j]/d[i][j-1];
    }

   if(re!=0.0) nu=1./re; else nu=DBL_MAX;  // corresponding kinematic viscosity
 
   /*
    double dxyzmax=DBL_MIN/2, numax=nu ;
    for(i=1;i<=gridp[0];i++) for(j=1;j<=gridp[1];j++) for(k=1;k<=gridp[2];k++) 
    dxyzmax = max(dxyzmax,1./sqr(d[0][i])+1./sqr(d[1][j])+1./sqr(d[2][k])); 
    if(CompTemp)                          numax=max(nu/prandtl,numax) ;
    if(CompChem) for(int l=0;l<nchem;l++) numax=max(chemc[l]  ,numax) ; 
    deltmax=min(deltmax,tfdiff/(numax*dxyzmax)) ;
    */

}


Inflow::Inflow () {i=j=k=nchem=0;chem=0;u=v=w=t=0;}
Inflow::~Inflow() {if(chem) delete[] chem;}


void Inflow::operator=(const Inflow& inflw) {
  i=inflw.i; j=inflw.j; k=inflw.k;
  u=inflw.u; v=inflw.v; w=inflw.w;
  t=inflw.t;

  nchem=inflw.nchem;
  if(chem) delete chem;
  if(nchem>0) {
    chem=new NS_REAL[nchem];
    for(int l=0; l<nchem; l++) chem[l]=inflw.chem[l];
   }
}


void Inflow::Print(char *head) {
  printf ("%s u(%d,%d,%d)=<%e,%e,%e>\n",(head) ? head : "", i,j,k,u,v,w) ;
  printf ("t=%e\n",t) ;
  for (int l=0; l<nchem; l++) printf ("%e ",chem[l]);
  printf("\n") ;
}


int Inflow::WriteFile(FILE* f) {
    NS_REAL* temp=chem;
    chem=0;
    if(fwrite(this,sizeof(Inflow),1,f)!=1) return FALSE;
    chem=temp;
    if(nchem>0) if(fwrite(chem,nchem,sizeof(NS_REAL),f)!=sizeof(NS_REAL)) return FALSE;
    return TRUE;
}


int Inflow::ReadFile(FILE* f) {
  if(fread(this,sizeof(Inflow),1,f)!=1) return FALSE;
  if(nchem>0) {
    chem=new NS_REAL[nchem];
    if(fread(chem,nchem,sizeof(NS_REAL),f)!=sizeof(NS_REAL)) return FALSE;
   }
  return TRUE;
}


void Inflow::SetCoord(int ii,int jj,int kk) {i=ii;j=jj;k=kk;}


InOutRec::InOutRec() {i=j=k=type=0;param=0;}
int InOutRec::WriteFile(FILE* f) {
  if(fwrite(this,sizeof(InOutRec),1,f)!=1) return FALSE;
  return TRUE;
}


int InOutRec::ReadFile(FILE* f) {
  if(fread(this,sizeof(InOutRec),1,f)!=1) return FALSE;
  return TRUE;
}

//
void InOutRec::operator=(const InOutRec& ior) {i=ior.i;j=ior.j;k=ior.k;type=ior.type;param=ior.param;}

void InOutRec::SetCoord(int ii,int jj,int kk) {i=ii;j=jj;k=kk;}

//
int InOutRecData::WriteFile(FILE* f) {
  if(fwrite(this,sizeof(InOutRecData),1,f)!=1) return FALSE;
  return TRUE;
}

int InOutRecData::ReadFile(FILE* f) {
  if(fread(this,sizeof(InOutRecData),1,f)!=1) return FALSE;
  return TRUE;
}

//
int CellRec::WriteFile(FILE* f) {
  if(fwrite(this,sizeof(CellRec),1,f)!=1) return FALSE;
  return TRUE;
}

int CellRec::ReadFile(FILE* f) {
  if(fread(this,sizeof(CellRec),1,f)!=1) return FALSE;
  return TRUE;
}

int SpecCellRec::WriteFile(FILE* f) {
  if(fwrite(this,sizeof(CellRec),1,f)!=1) return FALSE;
  return TRUE;
}

int SpecCellRec::ReadFile(FILE* f) {
  if(fread(this,sizeof(CellRec),1,f)!=1) return FALSE;
  return TRUE;
}

void SwapEndian(char* x, int size, int n, unsigned swap)
{
    if(!swap)
      return;

    char	one_byte;
    char	*pos;
    int		i, j;
    int		num= size;

    pos=	x;
    for(j= 0; j< n; j++)
    {
      for (i = 0; i < num/2; i++)
      {
          one_byte= pos[i];
	  pos[i]= pos[num-i-1];
          pos[num-i-1]= one_byte;
      }
      pos+= num;
    }
}

int  CheckEndian(void)
{
	unsigned char	EndianTest[2] = { 1, 0 };
	short x;

	x = *(short *) EndianTest;
	if(x > 1)
		return Big;
	else
		return little;
}

