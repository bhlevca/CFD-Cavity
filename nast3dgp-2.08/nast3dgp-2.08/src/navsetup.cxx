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
#include "navier.h"
#include "object.h"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define TEXTFILEDOUBLEFORMAT "%.10g "

unsigned NavierSetup::swapEndian(0);

NavierSetup::NavierSetup():Navier() {
    char* envstr;
    envstr=getenv("NAV_SCENES");
    if(envstr) strcpy(scenefile,envstr); else scenefile[0]=0;
    strcat(scenefile,"default.nav");
    envstr=getenv("NAV_VISUAL");
    if(envstr) strcpy(outfile,envstr); else outfile[0]=0;
    strcat(outfile,"default");
    output=0;
    desiredEndian= Big;
    systemEndian= (Endian)CheckEndian();
    swapEndian= (systemEndian == desiredEndian) ? FALSE : TRUE;
    loadfile[0]=0;loadmode=0;write_veloc_field=save_mem=addU=addV=addW=addT=addP=addCF=FALSE;
    for(int i=0; i<=98; i++) addC[i]=FALSE;
    output_mode=0;
}

void NavierSetup::ParseArgs(int argc, char** argv) {
    int      error=FALSE,n=0, i= 0, cur= 0;
    char*    env_str;   
    char     aux[256];
 
    while (--argc>0)
        if (argv[argc][0]=='-')
            switch (argv[argc][1]) {
            case 'b': // binary data file
                binfile[0]  =0 ; 
                env_str=getenv("NAV_BIN");
                if(env_str) {if(argv[argc+1][0]!='/' && argv[argc+1][0]!='\\') strcpy(binfile,env_str); } 
                strcat(binfile, argv[argc+1]) ;
                break;
            case 'g': // read binary data file and write as MatLab readable files
                binfile[0]  =0 ;
                env_str=getenv("NAV_BIN");
                if(env_str) {if(argv[argc+1][0]!='/' && argv[argc+1][0]!='\\') strcpy(binfile,env_str); }
                strcat(binfile,argv[argc+1]);
                output=1;
                break;
            case 'a': // read UVWPTC from text file
                binfile[0]  =0 ;
                env_str=getenv("NAV_BIN");
                if(env_str) {if(argv[argc+1][0]!='/' && argv[argc+1][0]!='\\') strcpy(binfile,env_str); }
                strcat(binfile,argv[argc+1]);
                output=3;
                break;
            case 'G': // read binary data file and write as VTK rectilinear grid files
                if(argv[argc][2]=='S') output_mode=1;      // write as stuctured points
                binfile[0]  =0 ;
                env_str=getenv("NAV_BIN");
                if(env_str) {if(argv[argc+1][0]!='/' && argv[argc+1][0]!='\\') strcpy(binfile,env_str); }
                strcat(binfile,argv[argc+1]);
                output=2;
                break;
	    case 'T': // read binary data file and write as Tecplot ASCII file
                if(argv[argc][2]=='C') output_mode=1;      // write as stuctured points
                binfile[0]  =0 ;
                env_str=getenv("NAV_BIN");
                if(env_str) {if(argv[argc+1][0]!='/' && argv[argc+1][0]!='\\') strcpy(binfile,env_str); }
                strcat(binfile,argv[argc+1]);
                output=4;
                break;
            case 'o': // output filename for -g, -G, -GS, -T and -a
                outfile[0]  =0 ;
                env_str=getenv("NAV_VISUAL");
                if(env_str) {if(argv[argc+1][0]!='/' && argv[argc+1][0]!='\\') strcpy(outfile,env_str); } 
                strcat(outfile,argv[argc+1]);
                break; 
            case 's': // scene (.nav) file
                scenefile[0]=0 ;
                env_str=getenv("NAV_SCENES");
                if(env_str) {if(argv[argc+1][0]!='/' && argv[argc+1][0]!='\\') strcpy(scenefile,env_str); }
                strcat(scenefile, argv[argc+1]);
                break;
            case 'l': // read UVWPTC from binary data file
                loadfile[0] =0 ;
                env_str=getenv("NAV_BIN");
                if(env_str) {if(argv[argc+1][0]!='/' && argv[argc+1][0]!='\\') strcpy(loadfile,env_str); }
                strcat(loadfile, argv[argc+1]);loadmode=0;
                break;
            case 'L': // read UVWPTC from MatLab readable files file.u, file.v,...
                loadfile[0] =0 ;
                env_str=getenv("NAV_BIN");
                if(env_str) {if(argv[argc+1][0]!='/' && argv[argc+1][0]!='\\') strcpy(loadfile,env_str); }
                strcat(loadfile, argv[argc+1]);loadmode=1;
                break;
            case 'A': // read UVWPTC from text file
                loadfile[0] =0 ;
                env_str=getenv("NAV_BIN");
                if(env_str) {if(argv[argc+1][0]!='/' && argv[argc+1][0]!='\\') strcpy(loadfile,env_str); }
                strcat(loadfile, argv[argc+1]);loadmode=2;
                break;
            case 'O': // Use certain byte order instead of default
            	n=strlen(argv[argc+1]);
            	cur= 0;
            	while(argv[argc+1][cur] != 0 && cur < n) {
            		switch (argv[argc+1][cur]) {
            		case 'B':
            			desiredEndian= Big;
            			swapEndian= (systemEndian == desiredEndian) ? FALSE : TRUE;
            			cur++;
            			break;
            		case 'l':
            			desiredEndian= little;
            			swapEndian= (systemEndian == desiredEndian) ? FALSE : TRUE;
            			cur++;
            			break;
            		}
            	}
            	break;
 		case 'f': // create VTK-file wich contains the specified data
                n=strlen(argv[argc+1]);
		cur= 0;
	        while(argv[argc+1][cur] != 0 && cur < n) {
                switch (argv[argc+1][cur]) {
		case 'V':
		    write_veloc_field=TRUE;
		    cur++;
                    break;
		case 'u':
		    addU=TRUE;
		    cur++;
		    break;
		case 'v':
		    addV=TRUE;
		    cur++;
		    break;
		case 'w':
		    addW=TRUE;                
		    cur++;
                    break;
		case 'p':
		    addP=TRUE;
		    cur++;
		    break;
		case 't':
		    addT=TRUE;
		    cur++;
		    break;
		case 'c':
		    addCF=TRUE;
		    cur++;
		    strcpy(aux, (const char*) &(argv[argc+1][cur]));
		    i=0;
		    while(isalnum(aux[i]) && !(isalpha(aux[i])))	// while numeric
                      i++;
		    aux[i]= 0;
		    addC[atoi(aux)]=TRUE;
		    cur+= i;
                    break;
		  default:
		      error=TRUE;
		}
	      }
		break;
            case 'm': 
                save_mem=TRUE;
                break;
            case 'h':
            case '?':
                error=TRUE;
                break;
            default:
                error=TRUE;
                break;
            }
    if (error==TRUE) {
	printf("This is NAVSETUP version %s\n", VERSION) ;
	printf("Parameter:\n");
	printf(" -b <file>     binary data file\n") ;
	printf(" -O  <B|l>     Determine byte order of binary data output\n");
	printf("               Use Option \"B\" for bigendian or \"l\" for littleendian\n");
	printf("               The default is big endian, which is the most common for use with VTK\n");
	printf(" -g <file>     read binary data file and write as MatLab readable files\n");
	printf(" -T <file>     read binary data file and write as tecplot files\n");
	printf(" -TC <file>    read binary data file and write all data in one tecplot file\n");
	printf(" -G <file>     read binary data file and write as VTK rectilinear grid files\n");
	printf(" -GS <file>    read binary data file and write as VTK structured points files\n");
	printf(" -o <file>     filename for -g, -G, -GS and -a\n");
	printf(" -s <file>     scene (.nav) file\n" );
	printf(" -l <file>     read UVWPTC from binary data file\n");
	printf(" -L <file>     read UVWPTC from MatLab readable files file.u, file.v,...\n");
	printf(" -f <uvwptVc0c1...> (you must  -G) create VTK-file wich contains the specified data\n");
	printf("                    u,v,w: the velocity components, p: pressure, t: temperature, V: velocity vector, c<i>: scalar\n") ;
	printf("                    for example: navsetup -f utc1 -G <file> writes u, the temperature and scalar 1 in a vtk-file\n") ;
	printf("                                 navsetup -f Vp   -G <file> writes u,v,w as vectors and p as a scalar in a vtk-file\n") ; 
	printf(" -a <file>     write Text files\n");
	printf(" -A <file>     read UVWPTC from text file\n\n");
	printf("Binary data files are read from directory specified in the environment variable NAV_BIN\n");
	printf("scene files from NAV_SCENES, and MatLab readable files are written to NAV_VISUAL.\n");
	exit(1);
    }
}

// --------------------------------------------------
// Explorer file functions (Matlab?)
// --------------------------------------------------

int NavierSetup::WriteExplorerHead(const char* filename,const char* extension,FILE** f,const Scene& S,int pos) {
    char* name=new char[strlen(filename)+strlen(extension)+2];assert(name);
    int i;
    float fld;
    sprintf(name,"%s.%s",filename,extension);
    *f=fopen(name,"wb");
    printf("Writing %s...\n",name);
    delete name;
    if(!(*f)) return FALSE;
    if(pos==-1) {                                       // write head for UVW file
        i=S.gridp[0]+1;if(fwrite(&i,sizeof(int),1,*f)!=1) return FALSE;
        i=S.gridp[1]+1;if(fwrite(&i,sizeof(int),1,*f)!=1) return FALSE;
        i=S.gridp[2]+1;if(fwrite(&i,sizeof(int),1,*f)!=1) return FALSE;
        for(i=0;i<S.gridp[0]+1; i++) {
            fld=(float)0.5*(float)(S.kabs[0][i]+S.kabs[0][i+1]);
            if(fwrite(&fld,sizeof(float),1,*f)!=1) return FALSE;
        }
        for(i=0;i<S.gridp[1]+1;i++) {
            fld=(float)0.5*(float)(S.kabs[1][i]+S.kabs[1][i+1]);
            if(fwrite(&fld,sizeof(float),1,*f)!=1) return FALSE;
        }
        for(i=0;i<S.gridp[2]+1;i++) {
            fld=(float)0.5*(float)(S.kabs[2][i]+S.kabs[2][i+1]);
            if(fwrite(&fld,sizeof(float),1,*f)!=1) return FALSE;
        }
        return TRUE;
    } else {                                            // write normal head
        int i;
        float fld;
        i=S.gridp[0]+2;fwrite(&i,sizeof(int),1,*f);
        i=S.gridp[1]+2;fwrite(&i,sizeof(int),1,*f);
        i=S.gridp[2]+2;fwrite(&i,sizeof(int),1,*f);
        for(i=0;i<=S.gridp[0]+1; i++) {
            fld=(float)((pos&1)?S.kabs[0][i+1]:0.5*(S.kabs[0][i]+S.kabs[0][i+1]));
            if(fwrite(&fld,sizeof(float),1,*f)!=1) return FALSE;
        }
        for(i=0;i<=S.gridp[1]+1;i++) {
            fld=(float)((pos&2)?S.kabs[1][i+1]:0.5*(S.kabs[1][i]+S.kabs[1][i+1]));
            if(fwrite(&fld,sizeof(float),1,*f)!=1) return FALSE;
        }
        for(i=0;i<=S.gridp[2]+1;i++) {
            fld=(float)((pos&4)?S.kabs[2][i+1]:0.5*(S.kabs[2][i]+S.kabs[2][i+1]));
            if(fwrite(&fld,sizeof(float),1,*f)!=1) return FALSE;
        }
    }
    return TRUE;
}    

int NavierSetup::WriteExplorerOneField(FILE** f,Matrix<NS_REAL>& mat,const Scene& S,int slicenum,int joker) {
// joker wird nur fuer VTK benoetigt
    int  i,j;
    float* fl=new float[S.gridp[0]+2];assert(fl);
    JALLLOOP {
        IALLLOOP fl[i]=(float)mat[i][j][slicenum];
        if((int)fwrite(fl,sizeof(float),S.gridp[0]+2,*f)!=S.gridp[0]+2) {delete fl;return FALSE;}
    }
    delete fl;
    return TRUE;
}

int NavierSetup::WriteExplorerVelocField(Matrix<NS_REAL>& U,Matrix<NS_REAL>& V,Matrix<NS_REAL>& W,const Scene& S,FILE** f,int slicenum) {
    int i,j;
    float* fl=new float[3*S.gridp[0]];assert(fl);
    for(j=1;j<=S.gridp[1];j++) {
        for(i=1;i<=S.gridp[0];i++) {
            fl[3*i-3]=(float)0.5*(float)(U[i-1][j][slicenum]+U[i][j][slicenum]);
            fl[3*i-2]=(float)0.5*(float)(V[i][j-1][slicenum]+V[i][j][slicenum]);
            fl[3*i-1]=(float)0.5*(float)(W[i][j][slicenum-1]+W[i][j][slicenum]);
        }
        if((int)fwrite(fl,sizeof(float),3*(S.gridp[0]+1),*f)!=3*(S.gridp[0]+1)) {delete fl;return FALSE;}
    }
    return TRUE;
}

int NavierSetup::ReadExplorerHead(int* dims,double** spaces,int which,const char* name,FILE** f) {
    char* filename=new char[strlen(name)+10];
    switch(which&7) {
    case 0:     sprintf(filename,"%s.u",name);break;
    case 1:     sprintf(filename,"%s.v",name);break;
    case 2:     sprintf(filename,"%s.w",name);break;
    case 3:     sprintf(filename,"%s.p",name);break;
    case 4:     sprintf(filename,"%s.t",name);break;
    case 5:     sprintf(filename,"%s.l",name);break;
    case 6:     sprintf(filename,"%s.c%d",name,which >> 3);break;
    default:    return FALSE;
    }
    printf("Reading MatLab file %s ...",filename);fflush(stdout);
    *f=fopen(filename,"rb");
    delete filename;
    if(!(*f)) {printf("not found.\n");return FALSE;}
    if(fread(dims,sizeof(int),3,*f)!=3) {printf("read error.\n");return FALSE;}       // read dimensions
    float* fl=new float[std::max(dims[0],std::max(dims[1],dims[2]))];
    int i;
    for(i=0;i<3;i++) {
        spaces[i]=new double[dims[i]];
        if((int)fread(fl,sizeof(float),dims[i],*f)!=dims[i]) {
            for(int j=0;j<3;j++) delete spaces[i];delete fl;
            printf("read error.\n");
            return FALSE;
        }
        for(int j=0;j<dims[i];j++) spaces[i][j]=(double)fl[j];
    }
    delete fl;
    return TRUE;
}

int NavierSetup::ReadExplorerNextSlice(Matrix<NS_REAL>& mat,FILE* f,int dimi,int dimj,int slice) {
    int i,j;
    float* fl=new float[dimi];
    for(j=0;j<dimj;j++) {
        if((int)fread(fl,sizeof(float),dimi,f)!=dimi) {delete fl;return FALSE;}
        for(i=0;i<dimi;i++) mat[i][j][1]=(NS_REAL)fl[i];
    }
    delete fl;
    return TRUE;
}

// --------------------------------------------------
// VTK file functions
// --------------------------------------------------

#ifdef VTK_WORDS_BIGENDIAN
void Swap4BERange(char *,int){}
#else
void Swap4BERange(char *mem_ptr1,int num) {
    char one_byte;
    char *pos;
    int i;
    pos = mem_ptr1;
    for (i = 0; i < num; i++) {
        one_byte = pos[0];pos[0] = pos[3];pos[3] = one_byte;
        one_byte = pos[1];pos[1] = pos[2];pos[2] = one_byte;
        pos+=4;
    }
}
#endif

int NavierSetup::WriteVTKHead(const char* filename,const char* extension,FILE** f,const Scene& S,int pos) {
    char* name=new char[strlen(filename)+strlen(extension)+6];assert(name);
    int i;
    float fld;
    sprintf(name,"%s.%s.vtk",filename,extension);
    *f=fopen(name,"wb");
    printf("Writing %s...\n",name);
    delete[] name;
    if(!(*f)) return FALSE;
    if(pos==-1) {                                       // write head for UVW file
        fprintf(*f,"# vtk DataFile Version 2.0\nnav\nBINARY\nDATASET RECTILINEAR_GRID\nDIMENSIONS %d %d %d\nX_COORDINATE %d float\n",
            S.gridp[0],S.gridp[1],S.gridp[2],S.gridp[0]);
        for(i=1;i<=S.gridp[0]; i++) {
            fld=0.5*(float)(S.kabs[0][i]+S.kabs[0][i+1]);
//            Swap4BERange((char*)&fld,1);
            SwapEndian((char*)&fld, sizeof(fld), 1, NavierSetup::swapEndian);
            if(fwrite(&fld,sizeof(float),1,*f)!=1) return FALSE;
        }
        fprintf(*f,"\nY_COORDINATE %d float\n",S.gridp[1]);
        for(i=1;i<=S.gridp[1];i++) {
            fld=0.5*(float)(S.kabs[1][i]+S.kabs[1][i+1]);
//            Swap4BERange((char*)&fld,1);
            SwapEndian((char*)&fld, sizeof(fld), 1, NavierSetup::swapEndian);
            if(fwrite(&fld,sizeof(float),1,*f)!=1) return FALSE;
        }
        fprintf(*f,"\nZ_COORDINATE %d float\n",S.gridp[2]);
        for(i=1;i<=S.gridp[2];i++) {
            fld=0.5*(float)(S.kabs[2][i]+S.kabs[2][i+1]);
//            Swap4BERange((char*)&fld,1);
            SwapEndian((char*)&fld, sizeof(fld), 1, NavierSetup::swapEndian);
            if(fwrite(&fld,sizeof(float),1,*f)!=1) return FALSE;
        }
        fprintf(*f,"\nPOINT_DATA %d\nVECTORS vectors float\n",
            (S.gridp[0])*(S.gridp[1])*(S.gridp[2]));
    } else {                                            // write normal head
        int i;
        float fld;
        fprintf(*f,"# vtk DataFile Version 2.0\nnav\nBINARY\nDATASET RECTILINEAR_GRID\nDIMENSIONS %d %d %d\nX_COORDINATE %d float\n",
            S.gridp[0],S.gridp[1],S.gridp[2],S.gridp[0]);
        for(i=1;i<=S.gridp[0]; i++) {
            fld=0.5*(float)(S.kabs[0][i]+S.kabs[0][i+1]);
//            Swap4BERange((char*)&fld,1);
            SwapEndian((char*)&fld, sizeof(fld), 1, NavierSetup::swapEndian);
            if(fwrite(&fld,sizeof(float),1,*f)!=1) return FALSE;
        }
        fprintf(*f,"\nY_COORDINATE %d float\n",S.gridp[1]);
        for(i=1;i<=S.gridp[1];i++) {
            fld=0.5*(float)(S.kabs[1][i]+S.kabs[1][i+1]);
//            Swap4BERange((char*)&fld,1);
            SwapEndian((char*)&fld, sizeof(fld), 1, NavierSetup::swapEndian);
            if(fwrite(&fld,sizeof(float),1,*f)!=1) return FALSE;
        }
        fprintf(*f,"\nZ_COORDINATE %d float\n",S.gridp[2]);
        for(i=1;i<=S.gridp[2];i++) {
            fld=0.5*(float)(S.kabs[2][i]+S.kabs[2][i+1]);
//            Swap4BERange((char*)&fld,1);
            SwapEndian((char*)&fld, sizeof(fld), 1, NavierSetup::swapEndian);
            if(fwrite(&fld,sizeof(float),1,*f)!=1) return FALSE;
        }
	fprintf(*f,"\nPOINT_DATA %d",(S.gridp[0])*(S.gridp[1])*(S.gridp[2]));
	if(strlen(extension)!=3) fprintf(*f,"\nSCALARS %s float\nlookup_table default\n",extension);
    }    
    return TRUE;
}    

int NavierSetup::WriteVTKOneField(FILE** f,Matrix<NS_REAL>& mat,const Scene& S,int slicenum,int pos) {
    int  i,j;
    float* fl=new float[S.gridp[0]];assert(fl);
    for(j=1; j<=S.gridp[1]; j++) {
        for(i=1; i<=S.gridp[0]; i++) {
	  switch(pos){
	  case 1:     fl[i-1]=0.5*((float)mat[i][j][slicenum]+(float)mat[i-1][j][slicenum]);break;    
	  case 2:     fl[i-1]=0.5*((float)mat[i][j][slicenum]+(float)mat[i][j-1][slicenum]);break;
	  case 4:     fl[i-1]=0.5*((float)mat[i][j][slicenum]+(float)mat[i][j][slicenum-1]);break;
	  case 0:     fl[i-1]=(float)mat[i][j][slicenum];break;  
	  }                                                    
	}
//        Swap4BERange((char*)fl,S.gridp[0]);
        SwapEndian((char*)fl, sizeof((*fl)), S.gridp[0], NavierSetup::swapEndian);
        if((int)fwrite(fl,sizeof(float),S.gridp[0],*f)!=S.gridp[0]) {delete fl;return FALSE;}
    }
    delete fl;
    return TRUE;
}

int NavierSetup::WriteVTKVelocField(Matrix<NS_REAL>& U,Matrix<NS_REAL>& V,Matrix<NS_REAL>& W,const Scene& S,FILE** f,int slicenum) {
    int i,j;
    float* fl=new float[3*S.gridp[0]];assert(fl);
    for(j=1;j<=S.gridp[1];j++) {
        for(i=1;i<=S.gridp[0];i++) {
            fl[3*i-3]=0.5*(float)(U[i-1][j][slicenum]+U[i][j][slicenum]);
            fl[3*i-2]=0.5*(float)(V[i][j-1][slicenum]+V[i][j][slicenum]);
            fl[3*i-1]=0.5*(float)(W[i][j][slicenum-1]+W[i][j][slicenum]);
        }
//        Swap4BERange((char*)fl,3*(S.gridp[0]));
        SwapEndian((char*)fl, sizeof((*fl)), 3*(S.gridp[0]), NavierSetup::swapEndian);
        if((int)fwrite(fl,sizeof(float),3*(S.gridp[0]),*f)!=3*(S.gridp[0])) {delete fl;return FALSE;}
    }
    return TRUE;
}

int NavierSetup::WriteVTKSPHead(const char* filename,const char* extension,FILE** f,const Scene& S,int pos) {
  if(pos==-1) {printf("Can't create velocity-field in structured points format\n\n"); exit(1);} // Routine muss noch geschrieben werden
    int i,j;
    for(j=0;j<3;j++) for(i=0;i<=S.gridp[j];i++) if(S.d[j][i]!=S.d[j][i+1])      // check for equidistant grid 
      NavierSetup::Error(ERR_DEFAULT,"Cannot write structured points for not equidistant grid.\n");
    
    char* name=new char[strlen(filename)+strlen(extension)+8];assert(name);
    sprintf(name,"%s.vtksp.%s",filename,extension);
    *f=fopen(name,"wb");
    printf("Writing %s...\n",name);
    delete name;
    if(!(*f)) return FALSE;

    fprintf(*f,"# vtk DataFile Version 2.0\nnav\nBINARY\ndataset STRUCTURED_POINTS\nDIMENSIONS %d %d %d\nSPACING %g %g %g\n",
        S.gridp[0],S.gridp[1],S.gridp[2],S.d[0][0],S.d[1][0],S.d[2][0]);
    fprintf(*f,"ORIGIN %g %g %g\n",
        0.5*(S.kabs[0][1]+S.kabs[0][2]),
        0.5*(S.kabs[1][1]+S.kabs[1][2]), 
        0.5*(S.kabs[2][1]+S.kabs[2][2]));
    fprintf(*f,"point_data %d\n",(S.gridp[0])*(S.gridp[1])*(S.gridp[2]));
    if(strlen(extension)!=3) fprintf(*f,"SCALARS %s unsigned_char\nlookup_table default\n",extension);
    return TRUE;
}    

int NavierSetup::WriteVTKSPOneField(FILE** f,Matrix<NS_REAL>& mat,const Scene& S,int slicenum,int pos) {
  int  i,j;
  unsigned char* fl=new unsigned char[S.gridp[0]];assert(fl);
  NS_REAL min=mat.Min(),max=mat.Max();
  for(j=1; j<=S.gridp[1]; j++) {
    for(i=1; i<=S.gridp[0]; i++) {
      switch(pos){ 
      case 1:     fl[i-1]=(unsigned char)(255*(max-0.5*(mat[i][j][slicenum]+mat[i-1][j][slicenum]))/(max-min));
      case 2:     fl[i-1]=(unsigned char)(255*(max-0.5*(mat[i][j][slicenum]+mat[i][j-1][slicenum]))/(max-min));
      case 4:     fl[i-1]=(unsigned char)(255*(max-0.5*(mat[i][j][slicenum]+mat[i][j][slicenum-1]))/(max-min));
      case 0:     fl[i-1]=(unsigned char)(255*(max-mat[i][j][slicenum])/(max-min));
      }
    } 
    if((int)fwrite(fl,sizeof(unsigned char),S.gridp[0],*f)!=S.gridp[0]) {delete fl;return FALSE;}
  }
  delete fl;
  return TRUE;
}


// --------------------------------------------------
// TECPLOT file functions
// --------------------------------------------------
// --- writes all data in one file ----
int NavierSetup::ConvertTecplotFile(char* filename, 
				    int(*writehead)(const char*,const char*,FILE**,const Scene& S,int pos)) 
{
    int i,j,k,n;
    long filepos;
    FILE* f=fopen(filename,"rb");
    if(!f) return FALSE;
    FILE* wf;
    
    
    if(!S.ReadFile(f)) return FALSE;
    S.DumpInfo(stdout);
    int dims[6]={0,S.gridp[0]+1,0,S.gridp[1]+1,0,S.gridp[2]+1};

    if(!writehead(outfile,"dat",&wf,S,0)) return FALSE;

    // reading flag-field from binary file
    if(!flag.Init(dims)) return FALSE;
    if(!U.Init(dims)) return FALSE;
    if(!flag.ReadPartial(f,dims)) return FALSE;
    IALLLOOP JALLLOOP KALLLOOP {
      U[i][j][k]=(NS_REAL)(flag[i][j][k]&(SLIP|INOUT|OBSTACLE));
      if(i==0 || i==S.gridp[0]+1 || j==0 || j==S.gridp[1]+1 || k==0 || k==S.gridp[2]+1) U[i][j][k]-=1024;
    }

    IJKALLLOOP flag[i][j][k]=(int) U[i][j][k];

    if(!U.Init(dims)) return FALSE;
    if(!U.ReadPartial(f,dims)) return FALSE;

    if(!V.Init(dims)) return FALSE;
    if(!V.ReadPartial(f,dims)) return FALSE;

    if(!W.Init(dims)) return FALSE;
    if(!W.ReadPartial(f,dims)) return FALSE;

    if(!P.Init(dims)) return FALSE;
    if(!P.ReadPartial(f,dims)) return FALSE;

    if(S.CompTemp) {
      if(!T.Init(dims)) return FALSE;
      if(!T.ReadPartial(f,dims)) return FALSE;
    }
 
    if(S.CompChem){
	CH =new Matrix<NS_REAL>[S.nchem];
	for(i=0;i<S.nchem;i++) {
	    //char fe[5];
	    //sprintf(fe,"c%d",i);      
	    if(!CH[i].Init(dims)) return FALSE;
	    if(!CH[i].ReadPartial(f,dims)) return FALSE;
	}
    }

    float u,v,w;
    for(int slicenum=1; slicenum<=S.gridp[2]; slicenum++) {
	for(j=1; j<=S.gridp[1]; j++) {
	    for(i=1; i<=S.gridp[0]; i++) {
		fprintf(wf,"\n%f %f %f ",
			0.5*(S.kabs[0][i]+S.kabs[0][i+1]),
			0.5*(S.kabs[1][j]+S.kabs[1][j+1]),
			0.5*(S.kabs[2][slicenum]+S.kabs[2][slicenum+1]));
		
		// cell-middle-point-evaluation 
		u=(float)0.5*(float)(U[i-1][j][slicenum]+U[i][j][slicenum]);
		v=(float)0.5*(float)(V[i][j-1][slicenum]+V[i][j][slicenum]);
		w=(float)0.5*(float)(W[i][j][slicenum-1]+W[i][j][slicenum]);
		
		fprintf(wf,"%e %e %e %e %i ",u , v, w, (float) P[i][j][slicenum], flag[i][j][slicenum]);
		if(S.CompTemp) { fprintf(wf,"%e ",(float) T[i][j][slicenum]);}
		if(S.CompChem) for(int n=0;n<S.nchem;n++) { fprintf(wf,"%e ",(float) CH[n][i][j][slicenum]);}
	    }
	}
    }
    fclose(wf);
    
    printf("\n");
    fclose(f);
    return TRUE;
}

// --- Tecplot head for file containing all data --- 
int NavierSetup::WriteTecplotHead2(const char* filename,const char* extension,FILE** wf,const Scene& S,int pos) {
  char* name=new char[strlen(filename)+strlen(extension)+6];assert(name);
  int i;
  float fld;
  sprintf(name,"%s.%s",filename,extension);
  *wf=fopen(name,"w");
  printf("Writing all data in one tecplot-file %s...\n",name);
  delete[] name;
  if(!(*wf)) return FALSE;
  fprintf(*wf,"VARIABLES = \"x\" \"y\" \"z\" \"u\" \"v\" \"w\" \"p\" \"flg\" ");
  if(S.CompTemp) fprintf(*wf,"\"t\" ");
  if(S.CompChem) for(int i=0;i<S.nchem;i++) fprintf(*wf,"\"Ch%i\" ",i+1);
  fprintf(*wf,"\nZONE I=%i, J=%i, K=%i F=POINT \n\n",S.gridp[0],S.gridp[1],S.gridp[2]);
 
  return TRUE;
}

// --- Tecplot head for splitted writing of data files ---  
int NavierSetup::WriteTecplotHead(const char* filename,const char* extension,FILE** f,const Scene& S,int pos) {
    char* name=new char[strlen(filename)+strlen(extension)+6];assert(name);
    sprintf(name,"%s.%s.dat",filename,extension);
    *f=fopen(name,"w");
    if(extension!="uvw"){
	fprintf(*f,"VARIABLES = \"x\" \"y\" \"z\" \"%s\" \n",extension);
	fprintf(*f,"ZONE I=%4d, J=%4d, K=%4d, F=POINT\n",S.gridp[0]+2,S.gridp[1]+2,S.gridp[2]+2);
    }else{
	fprintf(*f,"VARIABLES = \"x\" \"y\" \"z\" \"u\" \"v\" \"w\" \n");
	fprintf(*f,"ZONE I=%4d, J=%4d, K=%4d, F=POINT\n",S.gridp[0],S.gridp[1],S.gridp[2]);
    }
    printf("Writing %s...\n",name);
    delete[] name;
    if(!(*f)) return FALSE;
    return TRUE;
}
                                            
// --- Tecplot routine for splitted writing of data files ---                                                                                                                        
int NavierSetup::WriteTecplotOneField(FILE** f,Matrix<NS_REAL>& mat,const Scene& S,int slicenum,int pos) {
    int  i,j;
                                                                                                                                                               
    for(j=0;j<=S.gridp[1]+1;j++)
	for(i=0;i<=S.gridp[0]+1;i++) {
	    fprintf(*f,"%f %f %f ",
		    ((pos&1) ? S.kabs[0][i+1] : 0.5*(S.kabs[0][i]+S.kabs[0][i+1])),
		    ((pos&2) ? S.kabs[1][j+1] : 0.5*(S.kabs[1][j]+S.kabs[1][j+1])),
		    ((pos&4) ? S.kabs[2][slicenum+1] : 0.5*(S.kabs[2][slicenum]+S.kabs[2][slicenum+1])));
	    fprintf(*f,TEXTFILEDOUBLEFORMAT,mat[i][j][slicenum]);
	    fprintf(*f,"\n");
	}
    
    return TRUE;
}
               
// --- Tecplot routine for writing velocity vectors ---                                                                                      
int NavierSetup::WriteTecplotVelocField(Matrix<NS_REAL>& U,Matrix<NS_REAL>& V,Matrix<NS_REAL>& W,const Scene& S,FILE** f,int slicenum) {
    int i,j,k;
    float u,v,w;
    
    if (slicenum<=S.gridp[2]) {
	for(j=1;j<=S.gridp[1];j++) {
	    for(i=1;i<=S.gridp[0];i++) {
		// cell-middle-point-evaluation 
		u=(float)0.5*(float)(U[i-1][j][slicenum]+U[i][j][slicenum]);
		v=(float)0.5*(float)(V[i][j-1][slicenum]+V[i][j][slicenum]);
		w=(float)0.5*(float)(W[i][j][slicenum-1]+W[i][j][slicenum]);

		fprintf(*f,"%f %f %f ", 
			0.5*(S.kabs[0][i]+S.kabs[0][i+1]), 
			0.5*(S.kabs[1][j]+S.kabs[1][j+1]), 
			0.5*(S.kabs[2][slicenum]+S.kabs[2][slicenum+1]));
		fprintf(*f,TEXTFILEDOUBLEFORMAT,u);
		fprintf(*f,TEXTFILEDOUBLEFORMAT,v);
		fprintf(*f,TEXTFILEDOUBLEFORMAT,w);
		fprintf(*f,"\n");
	    }
	}
    }
    
    return TRUE;
}


// --------------------------------------------------
//        BIN file functions
// --------------------------------------------------

int NavierSetup::ReadBinHead(int* dims,double** spaces,int which,const char* name,FILE** f) {
    int ipmode;
    printf("Reading BIN file %s ",name);fflush(stdout);
    *f=fopen(name,"rb");
    if(!(*f)) {printf("not found.\n");return FALSE;}
    Scene myscene;
    if(!myscene.ReadFile(*f)) {printf("read error.\n");return FALSE;}
    int i,j;
    for(i=0;i<3;i++) dims[i]=myscene.gridp[i]+2;
    switch(which) {
    case 0: ipmode=1;break;
    case 1: ipmode=2;break;
    case 2: ipmode=4;break;
    default:ipmode=0;break;
    }
    for(i=0;i<3;i++) {
        spaces[i]=new double[dims[i]];
        for(j=0;j<dims[i];j++) spaces[i][j]=(ipmode & (1 << i))?myscene.kabs[i][j+1]:0.5*(myscene.kabs[i][j]+myscene.kabs[i][j+1]);
    }
    int matinit[6]={1,1,1,1,1,1};
    Matrix<NS_REAL> U;
    Matrix<unsigned> fl;
    fl.Init(matinit);
    fl.ReadPartial(*f);
    if((which&7)==0) {printf("(u)...");return TRUE;}
    U.Init(matinit);
    U.ReadPartial(*f);
    if((which&7)==1) {printf("(v)...");return TRUE;}
    U.ReadPartial(*f);
    if((which&7)==2) {printf("(w)...");return TRUE;}
    U.ReadPartial(*f);
    if((which&7)==3) {printf("(p)...");return TRUE;}
    U.ReadPartial(*f);
    if(myscene.CompTemp) {
        if((which&7)==4) {printf("(t)...");return TRUE;}
        U.ReadPartial(*f);
    }
    if(myscene.CompChem && (which&7)==6 && myscene.nchem>(which>>3)) {
        for(i=0;i<(which>>3);i++) U.ReadPartial(*f);
        printf("(ch%d)...",(which>>3));return TRUE;
    }
    return FALSE;
}

int NavierSetup::ReadBinNextSlice(Matrix<NS_REAL>& mat,FILE* f,int dimi,int dimj,int slice) {
    int i,j;
    Matrix<NS_REAL> U;
    int matinit[6]={0,dimi-1,0,dimj-1,slice,slice};
    U.Init(matinit);
    long filepos=ftell(f);
    U.ReadPartial(f);
    fseek(f,filepos,SEEK_SET);
    for(j=0;j<dimj;j++) for(i=0;i<dimi;i++) mat[i][j][1]=U[i][j][slice];
    return TRUE;
}

// --------------------------------------------------
// Text file functions
// --------------------------------------------------

int NavierSetup::WriteTextHead(const char* filename,const char* extension,FILE** f,const Scene& S,int pos) {
    char* name=new char[strlen(filename)+strlen(extension)+6];assert(name);
    int i;
    double fld;
    sprintf(name,"%s.txt.%s",filename,extension);
    *f=fopen(name,"w");
    printf("Writing %s...\n",name);
    delete name;
    if(!(*f)) return FALSE;
    if(pos==-1) {                                       // write head for UVW file
        i=S.gridp[0]+1;if(fprintf(*f,"%d\n",i)<0) return FALSE;
        i=S.gridp[1]+1;if(fprintf(*f,"%d\n",i)<0) return FALSE;
        i=S.gridp[2]+1;if(fprintf(*f,"%d\n",i)<0) return FALSE;
        for(i=0;i<S.gridp[0]+1; i++) {
            fld=0.5*(S.kabs[0][i]+S.kabs[0][i+1]);
            if(fprintf(*f,TEXTFILEDOUBLEFORMAT,fld)<0) return FALSE;
        }
        fprintf(*f,"\n");
        for(i=0;i<S.gridp[1]+1;i++) {
            fld=0.5*(S.kabs[1][i]+S.kabs[1][i+1]);
            if(fprintf(*f,TEXTFILEDOUBLEFORMAT,fld)<0) return FALSE;
        }
        fprintf(*f,"\n");
        for(i=0;i<S.gridp[2]+1;i++) {
            fld=0.5*(S.kabs[2][i]+S.kabs[2][i+1]);
            if(fprintf(*f,TEXTFILEDOUBLEFORMAT,fld)<0) return FALSE;
        }
        fprintf(*f,"\n");
        return TRUE;
    } else {                                            // write normal head
        int i;
        double fld;
        i=S.gridp[0]+2;if(fprintf(*f,"%d\n",i)<0) return FALSE;
        i=S.gridp[1]+2;if(fprintf(*f,"%d\n",i)<0) return FALSE;
        i=S.gridp[2]+2;if(fprintf(*f,"%d\n",i)<0) return FALSE;
        for(i=0;i<=S.gridp[0]+1; i++) {
            fld=((pos&1)?S.kabs[0][i+1]:0.5*(S.kabs[0][i]+S.kabs[0][i+1]));
            if(fprintf(*f,TEXTFILEDOUBLEFORMAT,fld)<0) return FALSE;
        }
        fprintf(*f,"\n");
        for(i=0;i<=S.gridp[1]+1;i++) {
            fld=((pos&2)?S.kabs[1][i+1]:0.5*(S.kabs[1][i]+S.kabs[1][i+1]));
            if(fprintf(*f,TEXTFILEDOUBLEFORMAT,fld)<0) return FALSE;
        }
        fprintf(*f,"\n");
        for(i=0;i<=S.gridp[2]+1;i++) {
            fld=((pos&4)?S.kabs[2][i+1]:0.5*(S.kabs[2][i]+S.kabs[2][i+1]));
            if(fprintf(*f,TEXTFILEDOUBLEFORMAT,fld)<0) return FALSE;
        }
        fprintf(*f,"\n");
    }
    return TRUE;
}    

int NavierSetup::WriteTextOneField(FILE** f,Matrix<NS_REAL>& mat,const Scene& S,int slicenum,int joker) {
// joker wird nur fuer VTK benoetigt
    int  i,j;
    JALLLOOP {
        fprintf(*f,"%4d %4d ",j,slicenum);
        IALLLOOP fprintf(*f,TEXTFILEDOUBLEFORMAT,mat[i][j][slicenum]);
        fprintf(*f,"\n");
    }
    return TRUE;
}

int NavierSetup::WriteTextVelocField(Matrix<NS_REAL>& U,Matrix<NS_REAL>& V,Matrix<NS_REAL>& W,const Scene& S,FILE** f,int slicenum) {
    int i,j;
    for(j=1;j<=S.gridp[1];j++) {
        fprintf(*f,"%4d %4d ",j,slicenum);
        for(i=1;i<=S.gridp[0];i++) {
            fprintf(*f,TEXTFILEDOUBLEFORMAT,0.5*(U[i-1][j][slicenum]+U[i][j][slicenum]));
            fprintf(*f,TEXTFILEDOUBLEFORMAT,0.5*(V[i][j-1][slicenum]+V[i][j][slicenum]));
            fprintf(*f,TEXTFILEDOUBLEFORMAT,0.5*(W[i][j][slicenum-1]+W[i][j][slicenum]));
        }
        fprintf(*f,"\n");
    }
    return TRUE;
}

int NavierSetup::ReadTextHead(int* dims,double** spaces,int which,const char* name,FILE** f) {
    char* filename=new char[strlen(name)+16];
    switch(which&7) {
    case 0:     sprintf(filename,"%s.txt.u",name);break;
    case 1:     sprintf(filename,"%s.txt.v",name);break;
    case 2:     sprintf(filename,"%s.txt.w",name);break;
    case 3:     sprintf(filename,"%s.txt.p",name);break;
    case 4:     sprintf(filename,"%s.txt.t",name);break;
    case 5:     sprintf(filename,"%s.txt.l",name);break;
    case 6:     sprintf(filename,"%s.txt.c%d",name,which >> 3);break;
    default:    return FALSE;
    }
    printf("Reading text file %s ...",filename);fflush(stdout);
    *f=fopen(filename,"r");
    delete filename;
    if(!(*f)) {printf("not found.\n");return FALSE;}
    if(fscanf(*f,"%d\n%d\n%d\n",&dims[0],&dims[1],&dims[2])<1) {printf("read error.\n");return FALSE;}       // read dimensions
    int i;
    double d;
    for(i=0;i<3;i++) {
        spaces[i]=new double[dims[i]];
        for(int j=0;j<dims[i];j++) {
            if(fscanf(*f,"%le ",&d)<1) {
                for(int j=0;j<3;j++) delete spaces[i];
                printf("read error.\n");
                return FALSE;
            }
            spaces[i][j]=d;
        }
        fscanf(*f,"\n");
    }
    return TRUE;
}

int NavierSetup::ReadTextNextSlice(Matrix<NS_REAL>& mat,FILE* f,int dimi,int dimj,int slice) {
    int i,j,k,l;
    double d;
    for(j=0;j<dimj;j++) {
        if(fscanf(f,"%d %d ",&k,&l)<1) return FALSE;
        assert(k==j && l==slice);
        for(i=0;i<dimi;i++) {
            if(fscanf(f,"%le ",&d)<1) return FALSE;
            mat[i][j][1]=(NS_REAL)d;
        }
        fscanf(f,"\n");
    }
    return TRUE;
}

// --------------------------------------------------
// common file functions
// --------------------------------------------------

int NavierSetup::ConvertOneFile(FILE* f,const char* extension,int(*writehead)(const char*,const char*,FILE**,const Scene& S,int pos),
                                int(*writeonefield)(FILE**,Matrix<NS_REAL>&,const Scene&,int,int),int pos) {
    long filepos=ftell(f);
    int k;
    FILE* wf;
    NS_REAL Min=MAXNS_REAL,Max=MINNS_REAL/2;
    int dims[6]={0,S.gridp[0]+1,0,S.gridp[1]+1,0,S.gridp[2]+1};
    if(!writehead(outfile,extension,&wf,S,pos)) return FALSE;
    //Marke
    if(save_mem) {
        for(k=0;k<=S.gridp[2]+1;k++) {
            dims[4]=dims[5]=k;
            if(!U.Init(dims)) return FALSE;
            fseek(f,filepos,SEEK_SET);
            if(!U.ReadPartial(f,dims)) return FALSE;
            Min=std::min(Min,U.Min());
            Max=std::max(Max,U.Max());
            if(!writeonefield(&wf,U,S,k,pos)) return FALSE;
        }
    } else {
        if(!U.Init(dims)) return FALSE;
        if(!U.ReadPartial(f,dims)) return FALSE;
        Min=U.Min();
        Max=U.Max();
        
        for(k=0;k<=S.gridp[2]+1;k++) if(!writeonefield(&wf,U,S,k,pos)) return FALSE;        
    }
    fclose(wf);
    printf("%s: min  %10.5e   max  %10.5e\n",extension,Min,Max);
    return TRUE;
}

int NavierSetup::ConvertFile(char* filename,int(*writehead)(const char*,const char*,FILE**,const Scene& S,int pos),
                             int(*writeonefield)(FILE**,Matrix<NS_REAL>&,const Scene&,int,int),
                             int(*writevelocfield)(Matrix<NS_REAL>&,Matrix<NS_REAL>&,Matrix<NS_REAL>&,const Scene&,FILE**,int)) {
    int i,j,k;
    long filepos, filepos_l;
    FILE* f=fopen(filename,"rb");
    if(!f) return FALSE;
    FILE* wf;
    if(!S.ReadFile(f)) return FALSE;
    S.DumpInfo(stdout);
    int dims[6]={0,S.gridp[0]+1,0,S.gridp[1]+1,0,S.gridp[2]+1};

    filepos=ftell(f);
    if(!writehead(outfile,"fl",&wf,S,0)) return FALSE;
    if(save_mem) {
        for(k=0;k<=S.gridp[2]+1;k++) {
            dims[4]=dims[5]=k;
            if(!flag.Init(dims)) return FALSE;
            if(!U.Init(dims)) return FALSE;
            fseek(f,filepos,SEEK_SET);
            if(!flag.ReadPartial(f,dims)) return FALSE;
            IALLLOOP JALLLOOP {
                U[i][j][k]=(NS_REAL)(flag[i][j][k]&(SLIP|INOUT|OBSTACLE)); 
                if(i==0 || i==S.gridp[0]+1 || j==0 || j==S.gridp[1]+1 || k==0 || k==S.gridp[2]+1) U[i][j][k]-=1024;
            }
            if(!writeonefield(&wf,U,S,k,0)) return FALSE;
        }
    } else {
        if(!flag.Init(dims)) return FALSE;
        if(!U.Init(dims)) return FALSE;
        if(!flag.ReadPartial(f,dims)) return FALSE;
        IALLLOOP JALLLOOP KALLLOOP {
            U[i][j][k]=(NS_REAL)(flag[i][j][k]&(SLIP|INOUT|OBSTACLE));
            if(i==0 || i==S.gridp[0]+1 || j==0 || j==S.gridp[1]+1 || k==0 || k==S.gridp[2]+1) U[i][j][k]-=1024;
        }
        for(k=0;k<=S.gridp[2]+1;k++) if(!writeonefield(&wf,U,S,k,0)) return FALSE; 
    }
    fclose(wf);
    filepos=ftell(f);  // remember U position
    if(!ConvertOneFile(f,"u",writehead,writeonefield,1)) return FALSE;
    if(!ConvertOneFile(f,"v",writehead,writeonefield,2)) return FALSE;
    if(!ConvertOneFile(f,"w",writehead,writeonefield,4)) return FALSE;
    if(!ConvertOneFile(f,"p",writehead,writeonefield,0)) return FALSE;
    if(S.CompTemp) if(!ConvertOneFile(f,"t",writehead,writeonefield,0)) return FALSE;
    filepos_l=ftell(f);
    if(S.CompChem) for(i=0;i<S.nchem;i++) {
        char fe[5];
        sprintf(fe,"c%d",i);
        if(!ConvertOneFile(f,fe,writehead,writeonefield,0)) return FALSE;
    }
    if(write_veloc_field) {
        if(!writehead(outfile,"uvw",&wf,S,-1)) return FALSE;
        for(k=1;k<=S.gridp[2];k++) {
            dims[4]=dims[5]=k;
            fseek(f,filepos,SEEK_SET);
            if(!U.Init(dims)) return FALSE;
            if(!U.ReadPartial(f,dims)) return FALSE;
            if(!V.Init(dims)) return FALSE;
            if(!V.ReadPartial(f,dims)) return FALSE;
            dims[4]--;
            if(!W.Init(dims)) return FALSE;
            if(!W.ReadPartial(f,dims)) return FALSE;
            if(!writevelocfield(U,V,W,S,&wf,k)) return FALSE;
        }
        fclose(wf);
    }
    printf("\n");
    fclose(f);
    return TRUE;
}


int NavierSetup::ConvertOneFileVTK(FILE* f,const char* extension,
                                   int(*writehead)(const char*,const char*,FILE**,const Scene& S,int pos),
                                   int(*writeonefield)(FILE**,Matrix<NS_REAL>&,const Scene&,int,int),
                                   int pos) {
    long filepos=ftell(f);
    int k;
    FILE* wf;
    NS_REAL Min=MAXNS_REAL,Max=MINNS_REAL/2;
    int dims[6]={0,S.gridp[0]+1,0,S.gridp[1]+1,0,S.gridp[2]+1};
    if(!writehead(outfile,extension,&wf,S,pos)) return FALSE;
    //Marke
    if(save_mem) {
        for(k=1;k<=S.gridp[2];k++) {
            dims[4]=dims[5]=k;
            if(!U.Init(dims)) return FALSE;
            fseek(f,filepos,SEEK_SET);
            if(!U.ReadPartial(f,dims)) return FALSE;
            Min=std::min(Min,U.Min());
            Max=std::max(Max,U.Max());
            if(!writeonefield(&wf,U,S,k,pos)) return FALSE;
        }
    } else {
        if(!U.Init(dims)) return FALSE;
        if(!U.ReadPartial(f,dims)) return FALSE;
        Min=U.Min();
        Max=U.Max();        
        for(k=1;k<=S.gridp[2];k++) if(!writeonefield(&wf,U,S,k,pos)) return FALSE;        
    }
    fclose(wf);
    printf("%s: min  %10.5e   max  %10.5e\n",extension,Min,Max);
    return TRUE;
}


int NavierSetup::ConvertFileVTK(char* filename,int(*writehead)(const char*,const char*,FILE**,const Scene& S,int pos),
                                int(*writeonefield)(FILE**,Matrix<NS_REAL>&,const Scene&,int,int),
                                int(*writevelocfield)(Matrix<NS_REAL>&,Matrix<NS_REAL>&,Matrix<NS_REAL>&,const Scene&,FILE**,int)) {
    int i,j,k;
    long filepos,filepos_l;
    FILE* f=fopen(filename,"rb");
    if(!f) return FALSE;
    FILE* wf;
    if(!S.ReadFile(f)) return FALSE;
    S.DumpInfo(stdout);
    int dims[6]={0,S.gridp[0]+1,0,S.gridp[1]+1,0,S.gridp[2]+1};

    filepos=ftell(f);
    printf("\n");
    if(!writehead(outfile,"fl",&wf,S,0)) return FALSE;
    if(save_mem) {
        for(k=1;k<=S.gridp[2];k++) {
            dims[4]=dims[5]=k;
            if(!flag.Init(dims)) return FALSE;
            if(!U.Init(dims)) return FALSE;
            fseek(f,filepos,SEEK_SET);
            if(!flag.ReadPartial(f)) return FALSE;
            IALLLOOP JALLLOOP {
                U[i][j][k]=(NS_REAL)(flag[i][j][k]&(SLIP|INOUT|OBSTACLE));
                if(i==0 || i==S.gridp[0]+1 || j==0 || j==S.gridp[1]+1 || k==0 || k==S.gridp[2]+1) U[i][j][k]-=1024;
            }
            if(!writeonefield(&wf,U,S,k,0)) return FALSE;
        }
    } else {
        if(!flag.Init(dims)) return FALSE;
        if(!U.Init(dims)) return FALSE;
        if(!flag.ReadPartial(f)) return FALSE;
        IALLLOOP JALLLOOP KALLLOOP {
            U[i][j][k]=(NS_REAL)(flag[i][j][k]&(SLIP|INOUT|OBSTACLE));
            if(i==0 || i==S.gridp[0]+1 || j==0 || j==S.gridp[1]+1 || k==0 || k==S.gridp[2]+1) U[i][j][k]-=1024;
        }
        for(k=1;k<=S.gridp[2];k++) if(!writeonefield(&wf,U,S,k,0)) return FALSE; 
    }
    fclose(wf);
    filepos=ftell(f);                                   // remember U position
    if(!ConvertOneFileVTK(f,"u",writehead,writeonefield,1)) return FALSE;
    if(!ConvertOneFileVTK(f,"v",writehead,writeonefield,2)) return FALSE;
    if(!ConvertOneFileVTK(f,"w",writehead,writeonefield,4)) return FALSE;
    if(!ConvertOneFileVTK(f,"p",writehead,writeonefield,0)) return FALSE;
    if(S.CompTemp) if(!ConvertOneFileVTK(f,"t",writehead,writeonefield,0)) return FALSE;   
    filepos_l=ftell(f);
    if(S.CompChem) for(i=0;i<S.nchem;i++) {
        char fe[5];
        sprintf(fe,"c%d",i);
        if(!ConvertOneFileVTK(f,fe,writehead,writeonefield,0)) return FALSE;
    }
    if(write_veloc_field) {
        if(!writehead(outfile,"mix",&wf,S,-1)) return FALSE;
        printf("   ...add velocity-field to mix-file...");
        for(k=1;k<=S.gridp[2];k++) {
            dims[4]=dims[5]=k;
            fseek(f,filepos,SEEK_SET);
            if(!U.Init(dims)) return FALSE;
            if(!U.ReadPartial(f,dims)) return FALSE;
            if(!V.Init(dims)) return FALSE;
            if(!V.ReadPartial(f,dims)) return FALSE;
            dims[4]--;
            if(!W.Init(dims)) return FALSE;
            if(!W.ReadPartial(f,dims)) return FALSE;
            if(!writevelocfield(U,V,W,S,&wf,k)) return FALSE;          
        }
	printf(" o.k.\n");
	dims[4]=0;
	dims[5]=S.gridp[2]+1;
    }
   
    if(addU || addV || addW || addP || addT || addCF){
      if(!write_veloc_field) if(!writehead(outfile,"mix",&wf,S,0)) return FALSE;
      fseek(f,filepos,SEEK_SET);
      if(!U.Init(dims)) return FALSE;
      if(!U.ReadPartial(f,dims)) return FALSE;
      if(addU) {
	printf("   ...add u-scalars to mix-file........");
	fprintf(wf,"\nSCALARS u float\nLOOKUP_TABLE default\n");
	for(k=1;k<=S.gridp[2];k++) if(!writeonefield(&wf,U,S,k,1)) return FALSE;
	printf(" o.k.\n");
      }
      if(!V.Init(dims)) return FALSE;
      if(!V.ReadPartial(f,dims)) return FALSE;
      if(addV) {
	printf("   ...add v-scalars to mix-file........");
	fprintf(wf,"\nSCALARS v float\nLOOKUP_TABLE default\n");
	for(k=1;k<=S.gridp[2];k++) if(!writeonefield(&wf,V,S,k,2)) return FALSE;
	printf(" o.k.\n");
      }
      if(!W.Init(dims)) return FALSE;
      if(!W.ReadPartial(f,dims)) return FALSE;
      if(addW) {
	printf("   ...add w-scalars to mix-file........");
	fprintf(wf,"\nSCALARS w float\nLOOKUP_TABLE default\n");
	for(k=1;k<=S.gridp[2];k++) if(!writeonefield(&wf,W,S,k,4)) return FALSE;
	printf(" o.k.\n");
      }
      if(!P.Init(dims)) return FALSE;
      if(!P.ReadPartial(f)) return FALSE;
      if(addP) {
	printf("   ...add p-scalars to mix-file........");
	fprintf(wf,"\nSCALARS p float\nLOOKUP_TABLE default\n");
	for(k=1;k<=S.gridp[2];k++) if(!writeonefield(&wf,P,S,k,0)) return FALSE;
	printf(" o.k.\n");
      }
      if(S.CompTemp){
	if(!T.Init(dims)) return FALSE;
	if(!T.ReadPartial(f)) return FALSE;
	if(addT) {
	  printf("   ...add t-scalars to mix-file........");
	  fprintf(wf,"\nSCALARS t float\nLOOKUP_TABLE default\n");
	  for(k=1;k<=S.gridp[2];k++) if(!writeonefield(&wf,T,S,k,0)) return FALSE;
	  printf(" o.k.\n");
	}
      }
      if(S.CompChem){     
	for(i=0;i<S.nchem;i++) {
	  char ec[5];                 
	  sprintf(ec,"c%d",i);
	  if(!T.Init(dims)) return FALSE;
	  if(!T.ReadPartial(f)) return FALSE; 
	  if(addC[i]) {
	    printf("   ...add %s-scalars to mix-file.......",ec);
	    fprintf(wf,"\nSCALARS %s float\nLOOKUP_TABLE default\n",ec);
	    for(k=1;k<=S.gridp[2];k++) if(!writeonefield(&wf,T,S,k,0)) return FALSE;
	    printf(" o.k.\n");
	  }
	}
      }
    }
    printf("\n");
    
    if(write_veloc_field || addU || addV || addW || addP || addT || addCF) fclose(wf);
    fclose(f);
    return TRUE;
}


#define LFA(ii,j) ( spaces[j][ii##p]-wspaces[j][ii])
#define LFB(ii,j) (wspaces[j][ii   ]- spaces[j][ii##p-1])
void NavierSetup::ReadFieldFromFile(FILE* f,int(*readhead)(int*,double**,int,const char*,FILE**),
                                     int(*readnextslice)(Matrix<NS_REAL>&,FILE*,int,int,int),int which,int ipmode) {
    double* spaces[3];
    double* wspaces[3];
    int     dims[3];
    FILE*   rf;
    int     i,j,k,ip,jp,kp;
    long filepos;

    filepos=ftell(f);
    if(!readhead(dims,spaces,which,loadfile,&rf)) return;
    if(spaces[0][dims[0]-1]<=S.dimension[0] || 
       spaces[1][dims[1]-1]<=S.dimension[1] || 
       spaces[2][dims[2]-1]<=S.dimension[2]) {
         printf("domain too small.\n");
         fclose(rf);
         for(i=0;i<3;i++) delete spaces[i];
         return;
    }
    int rmat_dims[6]={0,dims[0]-1,0,dims[1]-1,0,1};
    int wmat_dims[6]={0,S.gridp[0]+1,0,S.gridp[1]+1,0,S.gridp[2]+1};
    Matrix<NS_REAL> rmat;
    rmat.Init(rmat_dims);
    readnextslice(rmat,rf,dims[0],dims[1],0);
    for(i=0;i<dims[0];i++) for(j=0;j<dims[1];j++) rmat[i][j][0]=rmat[i][j][1];
    readnextslice(rmat,rf,dims[0],dims[1],1);
    kp=1;
    for(i=0;i<3;i++) {
        wspaces[i]=new double[S.gridp[i]+2];
        for(j=0;j<=S.gridp[i]+1;j++) {
            wspaces[i][j]=(ipmode & (1 << i))?S.kabs[i][j+1]:0.5*(S.kabs[i][j]+S.kabs[i][j+1]);
        }
    }
    if(save_mem) {
        for(k=0;k<=S.gridp[2]+1;k++) {
            wmat_dims[4]=wmat_dims[5]=k;
            U.Init(wmat_dims);
            while(wspaces[2][k]>=spaces[2][kp] && kp<dims[2]-1) {         // load new plane if necessary
                for(i=0;i<dims[0];i++) for(j=0;j<dims[1];j++) rmat[i][j][0]=rmat[i][j][1];
                readnextslice(rmat,rf,dims[0],dims[1],++kp);
            }
            for(i=0,ip=0;i<=S.gridp[0]+1;i++) {
                while(wspaces[0][i]>=spaces[0][ip]) ip++;
                if(ip>=dims[0]) ip=dims[0]-1;
                for(j=0,jp=0;j<=S.gridp[1]+1;j++) {
                    while(wspaces[1][j]>=spaces[1][jp]) jp++;
                    if(jp>=dims[1]) jp=dims[1]-1;
                    U[i][j][k]= rmat[ip-1][jp-1][0]*LFA(i,0)*LFA(j,1)*LFA(k,2)+
                                rmat[ip  ][jp-1][0]*LFB(i,0)*LFA(j,1)*LFA(k,2)+
                                rmat[ip-1][jp  ][0]*LFA(i,0)*LFB(j,1)*LFA(k,2)+
                                rmat[ip  ][jp  ][0]*LFB(i,0)*LFB(j,1)*LFA(k,2)+
                                rmat[ip-1][jp-1][1]*LFA(i,0)*LFA(j,1)*LFB(k,2)+
                                rmat[ip  ][jp-1][1]*LFB(i,0)*LFA(j,1)*LFB(k,2)+
                                rmat[ip-1][jp  ][1]*LFA(i,0)*LFB(j,1)*LFB(k,2)+
                                rmat[ip  ][jp  ][1]*LFB(i,0)*LFB(j,1)*LFB(k,2);
                    U[i][j][k]/=(spaces[0][ip]-spaces[0][ip-1])*(spaces[1][jp]-spaces[1][jp-1])*(spaces[2][kp]-spaces[2][kp-1]);
                }
            }
            fseek(f,filepos,SEEK_SET);
            if(!U.WritePartial(f)) NavierSetup::Error(ERR_WRITE,binfile);
        }
    } else {
        U.Init(wmat_dims);
        for(k=0;k<=S.gridp[2]+1;k++) {
            while(wspaces[2][k]>=spaces[2][kp] && kp<dims[2]-1) {         // load new plane if necessary
                for(i=0;i<dims[0];i++) for(j=0;j<dims[1];j++) rmat[i][j][0]=rmat[i][j][1];
                readnextslice(rmat,rf,dims[0],dims[1],++kp);
            }
            for(i=0,ip=1;i<=S.gridp[0]+1;i++) {
                while(wspaces[0][i]>=spaces[0][ip] && ip<dims[0]-1) ip++;
                for(j=0,jp=1;j<=S.gridp[1]+1;j++) {
                    while(wspaces[1][j]>=spaces[1][jp] && jp<dims[1]-1) jp++;
                    U[i][j][k]= rmat[ip-1][jp-1][0]*LFA(i,0)*LFA(j,1)*LFA(k,2)+
                                rmat[ip  ][jp-1][0]*LFB(i,0)*LFA(j,1)*LFA(k,2)+
                                rmat[ip-1][jp  ][0]*LFA(i,0)*LFB(j,1)*LFA(k,2)+
                                rmat[ip  ][jp  ][0]*LFB(i,0)*LFB(j,1)*LFA(k,2)+
                                rmat[ip-1][jp-1][1]*LFA(i,0)*LFA(j,1)*LFB(k,2)+
                                rmat[ip  ][jp-1][1]*LFB(i,0)*LFA(j,1)*LFB(k,2)+
                                rmat[ip-1][jp  ][1]*LFA(i,0)*LFB(j,1)*LFB(k,2)+
                                rmat[ip  ][jp  ][1]*LFB(i,0)*LFB(j,1)*LFB(k,2);
                    U[i][j][k]/=(spaces[0][ip]-spaces[0][ip-1])*(spaces[1][jp]-spaces[1][jp-1])*(spaces[2][kp]-spaces[2][kp-1]);
                }
            }
	}
	if(!U.WritePartial(f)) NavierSetup::Error(ERR_WRITE,binfile);
    }
    printf("ok.\n");fclose(rf);
    for(i=0;i<3;i++) {
        delete spaces[i];
        delete wspaces[i];
    }
}

#define LFA1(ii,j,v) (S.kabs[j][ii]-v)
#define LFB1(ii,j,v) (v-S.kabs[j][ii-1])

NS_REAL NavierSetup::GetValueAt(float x,float y,float z,Matrix<NS_REAL>& mat,const Scene& S) {
    int i=1,j=1,k=1;
    while(x>=S.kabs[0][i] && i<S.gridp[0]) i++;
    while(y>=S.kabs[1][j] && j<S.gridp[1]) j++;
    while(z>=S.kabs[2][k] && k<S.gridp[2]) k++;
    NS_REAL f=         mat[i-1][j-1][k-1]*LFA1(i,0,x)*LFA1(j,1,y)*LFA1(k,2,z)+
                    mat[i  ][j-1][k-1]*LFB1(i,0,x)*LFA1(j,1,y)*LFA1(k,2,z)+
                    mat[i-1][j  ][k-1]*LFA1(i,0,x)*LFB1(j,1,y)*LFA1(k,2,z)+
                    mat[i  ][j  ][k-1]*LFB1(i,0,x)*LFB1(j,1,y)*LFA1(k,2,z)+
                    mat[i-1][j-1][k  ]*LFA1(i,0,x)*LFA1(j,1,y)*LFB1(k,2,z)+
                    mat[i  ][j-1][k  ]*LFB1(i,0,x)*LFA1(j,1,y)*LFB1(k,2,z)+
                    mat[i-1][j  ][k  ]*LFA1(i,0,x)*LFB1(j,1,y)*LFB1(k,2,z)+
                    mat[i  ][j  ][k  ]*LFB1(i,0,x)*LFB1(j,1,y)*LFB1(k,2,z);
    return f/(DX[i]*DY[j]*DZ[k]);
}


void NavierSetup::ReadFieldsFromFile(FILE* f,int(*readhead)(int*,double**,int,const char*,FILE**),
                                     int(*readnextslice)(Matrix<NS_REAL>&,FILE*,int,int,int)) {
    ReadFieldFromFile(f,readhead,readnextslice,0,1);
    ReadFieldFromFile(f,readhead,readnextslice,1,2);
    ReadFieldFromFile(f,readhead,readnextslice,2,4);
    ReadFieldFromFile(f,readhead,readnextslice,3,0);
    if(S.CompTemp) ReadFieldFromFile(f,readhead,readnextslice,4,0);
    if(S.CompChem) for(int n=0;n<S.nchem;n++) ReadFieldFromFile(f,readhead,readnextslice,6+(n<<3),0);
}

void NavierSetup::DoIt(int argc,char* argv[]) {
    ParseArgs(argc,argv);
    switch(output) {
    case 1: // Explorer file routines
      if(!ConvertFile(binfile, NavierSetup::WriteExplorerHead, NavierSetup::WriteExplorerOneField, WriteExplorerVelocField)) 
	NavierSetup::Error(ERR_DEFAULT,"Error converting file.");
      break;
    case 2: // VTK file routines
      switch(output_mode) {
      case 0: // VTK rectilinear grid
	if(!ConvertFileVTK(binfile, NavierSetup::WriteVTKHead, NavierSetup::WriteVTKOneField, WriteVTKVelocField)) 
	  NavierSetup::Error(ERR_DEFAULT,"Error converting file.");
	break;
      case 1: // VTK structured points
	if(!ConvertFileVTK(binfile, NavierSetup::WriteVTKSPHead, NavierSetup::WriteVTKSPOneField, WriteVTKVelocField)) 
	  NavierSetup::Error(ERR_DEFAULT,"Error converting file.");
	break;
      }
      break;
    case 3: // ASCII file routines
      if(!ConvertFile(binfile, NavierSetup::WriteTextHead, NavierSetup::WriteTextOneField, WriteTextVelocField)) 
	NavierSetup::Error(ERR_DEFAULT,"Error converting file.");
      break;
    case 4: // Tecplot file routines
      switch(output_mode) {
      case 0: // write all fields in own files
	if(!ConvertFile(binfile,NavierSetup::WriteTecplotHead,NavierSetup::WriteTecplotOneField,WriteTecplotVelocField))
	  NavierSetup::Error(ERR_DEFAULT,"Error converting file.");
	break;
      case 1: // write all data in one file
	if(!ConvertTecplotFile(binfile, NavierSetup::WriteTecplotHead2)) 
	  NavierSetup::Error(ERR_DEFAULT,"Error converting file.");
	break;
      }
      break;
    case 0:  // Parser and binary file routines
      Object* objlist;
      FILE* f;
      printf("\nCreating new scene: %s -> %s\n",scenefile,binfile);
      
      // parse objects from file
      NavTokenizer t(scenefile);
      objlist=Parse(t);
      if(!objlist) NavierSetup::Error(ERR_DEFAULT,"Error parsing file.");
      
      if(S.alphatg==-1) S.alphatg=S.alpha;
      S.Init();
      f=fopen(binfile,"wb+");
      if(!f) NavierSetup::Error(ERR_OPEN,binfile);
      int matinit[6]={0,S.gridp[0]+1,0,S.gridp[1]+1,0,S.gridp[2]+1};
      
      printf("Generating flag field...\n");
      Matrix<flagtype> tempflagfield;
      tempflagfield.Init(matinit);
      
      objlist->CreatePreFlag(tempflagfield,S);
      CheckFlag(tempflagfield);    
      objlist->CreateInOutFlowList(tempflagfield,S);
      
      MarkBound(tempflagfield);
      S.ObstacleCount=objlist->CountInnerObstacles(tempflagfield,S);    
      S.DumpInfo(stdout);
      int i,k;                                                  
      
      if(!S.WriteFile(f)) NavierSetup::Error(ERR_WRITE,binfile);
      
      // 
      // Allocate Pointers to matrices
      Matrix <NS_REAL>* matrix[100];
      int n_matrices;   
      for (i=0;i<100;i++) matrix[i]=NULL;    
      
      matrix[0]=NULL; // place holder for flag flied
      matrix[1]=&U;
      matrix[2]=&V;
      matrix[3]=&W;
      matrix[4]=&P;
      n_matrices=5;
      if (S.CompTemp) matrix[n_matrices++]=&T;
      if (S.CompChem) {
	CH=new Matrix<NS_REAL>[S.nchem];
	for (i=0;i<S.nchem;i++) matrix[n_matrices++]=&(CH[i]);
      }
      
      if (n_matrices >100) NavierSetup::Error(ERR_WRITE,binfile);
      
      // initializes bounds for slices
      int slice_init[6];
      for (i=0;i<6;i++) slice_init[i]=matinit[i];
      
      long fileposu=-1;
      
      // Write all matrices:
      for (int n=0;n<n_matrices;n++) {
	if (n>0) assert(matrix[n]!=NULL) ;
	
	printf("\nWriting matrix %d \n",n);
	
	if (n==0){
	  if (!flag.WriteHead(f,matinit)) NavierSetup::Error(ERR_WRITE,binfile);
	}
	else
	  if (!matrix[n]->WriteHead(f,matinit)) NavierSetup::Error(ERR_WRITE,binfile);
	
	// Write Slices [i,:,:]
	for (i=0;i<=S.gridp[0]+1;i++)  {
	  slice_init[0]=slice_init[1]=i;
	  
	  // init U,V,W,P (without objects) with default values (zero)
	  flag.Init(slice_init);
	  for (k=1;k<5;k++) matrix[k]->Init(slice_init);
	  
	  // init temp
	  int pch=5 ;
	  if(S.CompTemp) {
	    matrix[pch]->Init(slice_init);
	    // use TempRef instead of zero as default initial value
	    for(int j=0;j<=S.gridp[1]+1;j++) 
	      for(int k=0;k<=S.gridp[2]+1;k++) (*matrix[pch])[i][j][k]=S.TempRef ;
	    // printf("Temp=%e\n",matrix[pch]->Min()) ;
	    pch++ ;
	  }
	  // init chem.
	  if(S.CompChem) {
	    for(k=0;k<S.nchem;k++) matrix[pch+k]->Init(slice_init);
	  }
	  
	  // CreateFlag initializes the slice for all fields U,V,W,.... 
	  // uncritical relict of an older code version
	  objlist->CreateFlag(tempflagfield,flag,U,V,W,P,T,CH,S,i);       
	  
	  printf("\rWriting slice %d ",i);fflush(stdout);
	  if(n==0){ 
	    if(!flag.WriteSlice(f,slice_init)) NavierSetup::Error(ERR_WRITE,binfile);
	    fileposu=ftell(f);  // we need the Position of matrix U later
	  }else
	    if(!matrix[n]->WriteSlice(f,slice_init)) NavierSetup::Error(ERR_WRITE,binfile);
	}
      }
    
      objlist->DeleteAll();
      delete objlist;
      printf("\n");       
      fseek(f,fileposu,SEEK_SET);
      if(loadfile[0]!=0) {
	switch(loadmode) {
	case 1: ReadFieldsFromFile(f,ReadExplorerHead,ReadExplorerNextSlice);break;
	case 0: ReadFieldsFromFile(f,ReadBinHead,ReadBinNextSlice);break;
	case 2: ReadFieldsFromFile(f,ReadTextHead,ReadTextNextSlice);break;
	}
      }
      fclose(f);
      break;
    }
}

int main(int argc,char* argv[]) {
    NavierSetup n;
    n.DoIt(argc,argv);
 return 0 ;
}


