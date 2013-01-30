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
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdarg.h>

int Navier::ReadFile(char* filename) {
    int i;
    FILE* f=fopen(filename,"rb");
    if(f==0) return FALSE;
    if(!S.ReadFile(f)) return FALSE;
    int dims[6]={0,S.gridp[0]+1,0,S.gridp[1]+1,0,S.gridp[2]+1};
    if(!flag.Init(dims)) return FALSE;
    if(!U.Init(dims)) return FALSE;
    if(!V.Init(dims)) return FALSE;
    if(!W.Init(dims)) return FALSE;
    if(!P.Init(dims)) return FALSE;
    if(S.CompTemp) if(!T.Init(dims)) return FALSE;
    if(S.CompChem) {
        if(CH) delete[] CH;
        CH=new Matrix<NS_REAL>[S.nchem];
        for(i=0;i<S.nchem;i++) if(!CH[i].Init(dims)) return FALSE;
    }
    if(!flag.Read(f)) return FALSE;
    if(!U.Read(f)) return FALSE;
    if(!V.Read(f)) return FALSE;
    if(!W.Read(f)) return FALSE;
    if(!P.Read(f)) return FALSE;
    if(S.CompTemp) if(!T.Read(f)) return FALSE;
    if(S.CompChem) for(i=0;i<S.nchem;i++) if(!CH[i].Read(f)) return FALSE;
    fclose(f);
    return TRUE;
}

int Navier::WriteFile(char* filename) {
    FILE* f=fopen(filename,"wb");
    if(f==0) return FALSE;
    if(!S.WriteFile(f)) return FALSE;
    if(!flag.Write(f)) return FALSE;
    if(!U.Write(f)) return FALSE;
    if(!V.Write(f)) return FALSE;
    if(!W.Write(f)) return FALSE;
    if(!P.Write(f)) return FALSE;
    if(S.CompTemp) if(!T.Write(f)) return FALSE;
    if(S.CompChem) for(int i=0;i<S.nchem;i++) if(!CH[i].Write(f)) return FALSE;
    fclose(f);
    return TRUE;
}

Navier::Navier() {
    CH=NULL;
    char* envstr=getenv("NAV_BIN");
    if(envstr) strcpy(binfile,envstr); else binfile[0]=0;
    strcat(binfile,"default.bin");
}

Navier::~Navier() {
  //    if(CH) delete[] CH;
}

void Navier::Error(Errors ErrNo,...) {
    va_list vlist;
    va_start(vlist,ErrNo);
    switch(ErrNo) {
        case ERR_OPEN:  fprintf(stderr,"Cannot open file %s.\n",va_arg(vlist,char*));break;
        case ERR_READ:  fprintf(stderr,"Cannot read from file %s.\n",va_arg(vlist,char*));break;
        case ERR_WRITE: fprintf(stderr,"Cannot write to file %s.\n",va_arg(vlist,char*));break;
        case ERR_DEFAULT: fprintf(stderr,va_arg(vlist,char*));break;        
        default:        fprintf(stderr,"Unknown error.\n");break;
    }
    va_end(vlist);
    exit(1);
}

