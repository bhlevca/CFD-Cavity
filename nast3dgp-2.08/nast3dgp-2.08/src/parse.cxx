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
#include "tokenizer.h"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define BETA_TOKEN          1
#define BOTTOM_TOKEN        2
#define BOX_TOKEN           3
#define COMMA_TOKEN         4
#define COORDS_TOKEN        5
#define DASH_TOKEN          6
#define DELT_TOKEN          7
#define DIFFERENCE_TOKEN    8
#define DIMENSION_TOKEN     9
#define EAST_TOKEN          10
#define END_OF_FILE_TOKEN   11
#define EPS_TOKEN           12
#define FLOAT_TOKEN         14
#define FLUID_TOKEN         15
#define GX_TOKEN            17
#define GY_TOKEN            18
#define GZ_TOKEN            19
#define INFLOW_TOKEN        20
#define INIT_TOKEN          21
#define INOUT_TOKEN         22
#define INTERSECTION_TOKEN  23
#define ITERMAX_TOKEN       24
#define LEFT_ANGLE_TOKEN    25
#define LEFT_CURLY_TOKEN    26
#define LEFT_PAREN_TOKEN    27
#define LENGTH_TOKEN        28
#define LENGTHREF_TOKEN     29
#define MAP_TOKEN           30
#define MAX_TOKEN           31
#define NORTH_TOKEN         32
#define NOSLIP_TOKEN        33
#define RE_TOKEN            34
#define OMEGA_TOKEN         35
#define PARAMS_TOKEN        36
#define PLUS_TOKEN          38
#define PRANDTL_TOKEN       39
#define PRSTEP_TOKEN        40
#define RELATIVE_TOKEN      41
#define RES_TOKEN           42
#define RIGHT_ANGLE_TOKEN   43
#define RIGHT_CURLY_TOKEN   44
#define RIGHT_PAREN_TOKEN   45
#define RHOREF_TOKEN        46
#define SLASH_TOKEN         47
#define SLIP_TOKEN          48
#define SOUTH_TOKEN         49
#define SPHERE_TOKEN        50
#define STREAK_TOKEN        51
#define STRING_TOKEN        52
#define TEMPERATURE_TOKEN   53
#define TEMPCOLD_TOKEN      54
#define TEMPHOT_TOKEN       55
#define TFIN_TOKEN          56
#define TIME_TOKEN          57
#define TOP_TOKEN           58
#define TRACE_TOKEN         59
#define UNION_TOKEN         60
#define VELOCREF_TOKEN      61
#define WEST_TOKEN          62
#define X_TOKEN             63
#define Y_TOKEN             64
#define Z_TOKEN             65
#define ZYLINDER_TOKEN      66
#define CWCACOMP_TOKEN      67
#define CHEMNUM_TOKEN       68
#define CHEMINIT_TOKEN      69
#define TFDIFF_TOKEN        73
#define TFCONV_TOKEN        74
#define PERIODX_TOKEN       75
#define PERIODY_TOKEN       76
#define PERIODZ_TOKEN       77
#define HESSE_TOKEN         78
#define HALFSPACE_TOKEN     79
#define ALPHA_TOKEN         80
#define COMMENT_TOKEN       81
#define ALPHATC_TOKEN       82
#define POINTS_TOKEN        83
#define VERTICES_TOKEN      84
#define POLY_TOKEN          85
#define INFLOWDATA_TOKEN    86
#define INITDATA_TOKEN      87

#define FROUDE_TOKEN        93
#define TIME_DISKR_TOKEN    97
#define DIM_LESS_TOKEN      98 

#define TEMPREF_TOKEN       101

#define EU1_TOKEN           103
#define AB2_TOKEN           104

#define ON_TOKEN            107
#define OFF_TOKEN           108

#define SOLVER_TOKEN        109

#define SOR_TOKEN           110
#define SSOR_TOKEN          111
#define REDBLACK_TOKEN      112
#define COLORSOR_TOKEN      113
#define COLORSSOR_TOKEN     114
#define BICGSTAB_TOKEN      115
#define BICGSTAB_JVK_TOKEN  116

#define CONVECTIVE_TOKEN    117

#define DC_TOKEN            118
#define HLPA_TOKEN          119
#define QUICK_TOKEN         120
#define SMART_TOKEN         121
#define VONOS_TOKEN         122

#define PRTTIME_TOKEN       123
#define TARGETDIR_TOKEN     124


Token tok[]={
    {"alpha",ALPHA_TOKEN},
    {"beta",BETA_TOKEN},
    {"bottom",BOTTOM_TOKEN},
    {"box",BOX_TOKEN},
    {"cheminit",CHEMINIT_TOKEN},
    {"nuC",CHEMNUM_TOKEN},
    {"coords",COORDS_TOKEN},
    {"deltmax",DELT_TOKEN},
    {"difference",DIFFERENCE_TOKEN},
    {"dimension",DIMENSION_TOKEN},
    {"east",EAST_TOKEN},
    {"eps",EPS_TOKEN},
    {"fluid",FLUID_TOKEN},
    {"gx",GX_TOKEN},
    {"gy",GY_TOKEN},
    {"gz",GZ_TOKEN},
    {"halfspace",HALFSPACE_TOKEN},
    {"hesse",HESSE_TOKEN},
    {"inflow",INFLOW_TOKEN},
    {"init",INIT_TOKEN},
    {"inout",INOUT_TOKEN},
    {"intersection",INTERSECTION_TOKEN},
    {"itermax",ITERMAX_TOKEN},
    {"length",LENGTH_TOKEN},
    {"north",NORTH_TOKEN},
    {"noslip",NOSLIP_TOKEN},
    {"reynolds",RE_TOKEN},
    {"omega",OMEGA_TOKEN},
    {"parameter",PARAMS_TOKEN},
    {"periodboundx",PERIODX_TOKEN},
    {"periodboundy",PERIODY_TOKEN},
    {"periodboundz",PERIODZ_TOKEN},
    {"prandtl",PRANDTL_TOKEN},
    {"froude",FROUDE_TOKEN},
    {"prstep",PRSTEP_TOKEN},
    {"absolute",RELATIVE_TOKEN},
    {"resolution",RES_TOKEN},
    {"slip",SLIP_TOKEN},
    {"south",SOUTH_TOKEN},
    {"sphere",SPHERE_TOKEN},
    {"temperature",TEMPERATURE_TOKEN},
    {"TempCold",TEMPCOLD_TOKEN},
    {"TempHot",TEMPHOT_TOKEN},
    {"TempRef",TEMPREF_TOKEN},
    {"tfconv",TFCONV_TOKEN},
    {"tfdiff",TFDIFF_TOKEN},
    {"Tfin",TFIN_TOKEN},
    {"TimeDis",TIME_DISKR_TOKEN},
    {"dimensionless",DIM_LESS_TOKEN},
    {"top",TOP_TOKEN},
    {"union",UNION_TOKEN},
    {"west",WEST_TOKEN},
    {"cylinder",ZYLINDER_TOKEN},
    {"x",X_TOKEN},
    {"y",Y_TOKEN},
    {"z",Z_TOKEN},
    {"<",LEFT_ANGLE_TOKEN},
    {",",COMMA_TOKEN},
    {">",RIGHT_ANGLE_TOKEN},
    {"}",RIGHT_CURLY_TOKEN},
    {"{",LEFT_CURLY_TOKEN},
    {"//",COMMENT_TOKEN},
    {"alphaTC",ALPHATC_TOKEN},
    {"points",POINTS_TOKEN},
    {"vertices",VERTICES_TOKEN},
    {"poly",POLY_TOKEN},
    {"InflowData",INFLOWDATA_TOKEN},
    {"InitData",INITDATA_TOKEN},
    {"ON", ON_TOKEN},
    {"OFF", OFF_TOKEN},
    {"EU1", EU1_TOKEN},
    {"AB2", AB2_TOKEN},
    {"PoissonSolver", SOLVER_TOKEN},
    {"SOR", SOR_TOKEN},
    {"SSOR", SSOR_TOKEN},
    {"RedBlack", REDBLACK_TOKEN},
    {"8ColorSOR",COLORSOR_TOKEN},
    {"8ColorSSOR",COLORSSOR_TOKEN},
    {"BiCGStab", BICGSTAB_TOKEN},
    {"BiCGStabJPC", BICGSTAB_JVK_TOKEN},
    {"ConvectiveTerms", CONVECTIVE_TOKEN},
    {"DC", DC_TOKEN},
    {"HLPA", HLPA_TOKEN},
    {"QUICK", QUICK_TOKEN},
    {"SMART", SMART_TOKEN},
    {"VONOS", VONOS_TOKEN},
    {"TimePrintStep", PRTTIME_TOKEN},
    {"TargetDirectory", TARGETDIR_TOKEN},
    {0,0}             
};


const char* syn_err="syntax error";
char* end_token_chars="{/<";

NavTokenizer::NavTokenizer(const char* c):Tokenizer(c) {
    SetEndTokenChars(end_token_chars);
}

void NavTokenizer::ParseVector(double* v,int size) {
    int i,Exit_Flag=FALSE;
    while(!Exit_Flag) {
        switch(GetToken()) {
            case LEFT_ANGLE_TOKEN:     
                v[0]=GetDouble();
                for(i=1;i<size;i++) {
                    ParseToken(COMMA_TOKEN);
                    v[i]=GetDouble();
                }
                if(GetToken()!=RIGHT_ANGLE_TOKEN) Error("'>' expected");
                Exit_Flag=TRUE;
                break;
            default:
                Error("vector read fail");
        }
    }
}

void NavTokenizer::ParseToken(int tokenid) {
    if(GetToken()!=tokenid) Error(syn_err);
}

int NavTokenizer::GetToken() {
    int token;
    do {
        token=Tokenizer::GetToken(tok);
        if(token==COMMENT_TOKEN) NewLine();
    } while(token==COMMENT_TOKEN);
    return token;
}

int NavTokenizer::ShowToken() {
    int token;
    do {
        token=Tokenizer::ShowToken(tok);
        if(token==COMMENT_TOKEN) NewLine();
    } while(token==COMMENT_TOKEN);
    return token;
}

enum {
    OBJTOKEN_INFLOW =0x0001,
    OBJTOKEN_TEMP   =0x0002,
    OBJTOKEN_CHEM   =0x0004
};

int Object::Parse(NavTokenizer& t,Scene& S) {
    int n;
    NS_REAL vector[3];
    switch(t.GetToken()) {
        case COORDS_TOKEN:
            t.ParseVector(vector,3);
            xmin=vector[0];ymin=vector[1];zmin=vector[2];
            t.ParseToken(COMMA_TOKEN);
            t.ParseVector(vector,3);
            xmax=vector[0],ymax=vector[1],zmax=vector[2];
            if(xmin>xmax) {
                NS_REAL x=xmin;xmin=xmax;xmax=x;
            }
            if(ymin>ymax) {
                NS_REAL y=ymin;ymin=ymax;ymax=y;
            }
            if(zmin>zmax) {
                NS_REAL z=zmin;zmin=zmax;zmax=z;
            }
            break;
        case FLUID_TOKEN:
            flag=0;
            break;    
        case NOSLIP_TOKEN:
            flag&=~SLIP;
            break;
        case SLIP_TOKEN:
            flag|=SLIP;
            break;
        case INOUT_TOKEN:
            inouttype=t.GetLong();
            if(inouttype<1 || inouttype>3) t.Error("Wrong INOUT type");
            t.ParseToken(COMMA_TOKEN);
            inoutparam=t.GetDouble();
	    if(inoutparam<0 || inoutparam>1) t.Error("Wrong INOUT type");
	    if(inoutparam==0) S.iobd=1; else S.iobd=2;
            flag|=INOUT;
            flag|=(inouttype-1)*INOUTTYP2;
            break;
        case INIT_TOKEN:
            t.ParseVector(initv,3);
            t.ParseToken(COMMA_TOKEN);
            initp=t.GetDouble();
            break;
        case INFLOW_TOKEN:
            t.ParseVector(inflv,3);
            add_inflow=TRUE;
            tokens_given|=OBJTOKEN_INFLOW;
	    if(S.iobd==0) S.iobd=1;
            break;
        case RELATIVE_TOKEN:
            is_relative=TRUE;
            break;
        case TEMPERATURE_TOKEN:
            initt=t.GetDouble();
            flag &=~TEMPB ;
            tokens_given|=OBJTOKEN_TEMP;
            break;
        case INITDATA_TOKEN:{
                int which_ifd=0;
                switch(t.GetToken()) {
                case X_TOKEN: which_ifd=0;break;
                case Y_TOKEN: which_ifd=1;break;
                case Z_TOKEN: which_ifd=2;break;
                default: t.Error("Illegal direction.");
                }
                t.ParseToken(COMMA_TOKEN);
                char* filename=t.GetString();
                if(initdata[which_ifd]) delete initdata[which_ifd];
                initdata[which_ifd]=new InflowData(filename); 
                delete filename;
                }
                break;
        case CHEMINIT_TOKEN:
            if(S.nchem>0) {
                if(!initc) initc=new NS_REAL[S.nchem];
                initc[0]=t.GetDouble();
                for(n=1;n<S.nchem;n++) {
                    t.ParseToken(COMMA_TOKEN);
                    initc[n]=t.GetDouble();
                }
            }
            tokens_given|=OBJTOKEN_CHEM;
            flag&=~CHEMB;
            break;  
        case RIGHT_CURLY_TOKEN:
            if(tokens_given&OBJTOKEN_INFLOW) {
                if(S.CompTemp && !(tokens_given&OBJTOKEN_TEMP)) 
                    t.Error("Inflow temperature has to be specified.");
                if(S.CompChem && !(tokens_given&OBJTOKEN_CHEM)) 
                    t.Error("Inflow chemical concentrations have to be specified.");
            }            
            return TRUE;
        default:
            t.Error(syn_err);
    }
    return FALSE;
}

int Box::Parse(NavTokenizer& t,Scene& S) {
    int inflowdata_set=0;
    initt=S.TempRef ; //FK
    t.ParseToken(LEFT_CURLY_TOKEN);
    while(TRUE) {
        switch(t.ShowToken()) {
            case NORTH_TOKEN:
                t.GetToken();
                wall=NORTH;
                break;
            case SOUTH_TOKEN:
                t.GetToken();
                wall=SOUTH;
                break;
            case WEST_TOKEN:
                t.GetToken();
                wall=WEST;
                break;
            case EAST_TOKEN:
                t.GetToken();
                wall=EAST;
                break;
            case TOP_TOKEN:
                t.GetToken();
                wall=TOP;
                break;
            case BOTTOM_TOKEN:
                t.GetToken();
                wall=BOTTOM;
                break;
            case INFLOWDATA_TOKEN:{
                t.GetToken();
                inflowdata_set=1;
                int which_ifd=0;
                switch(t.GetToken()) {
                case X_TOKEN: which_ifd=0;break;
                case Y_TOKEN: which_ifd=1;break;
                case Z_TOKEN: which_ifd=2;break;
                default: t.Error("Illegal direction.");
                }
                t.ParseToken(COMMA_TOKEN);
                char* filename=t.GetString();
                if(inflowdata[which_ifd]) delete inflowdata[which_ifd];
                inflowdata[which_ifd]=new InflowData(filename); 
                delete filename;
                }
                break;
            default:
                if(Object::Parse(t,S)) {
                    if(inflowdata_set && wall==0) t.Error("InflowData can only be specified for wall objects.");
                    return TRUE;
                }
        }
    }
}

int Sphere::Parse(NavTokenizer& t,Scene& S) {
  //initt=S.TempRef ; //FK
    t.ParseToken(LEFT_CURLY_TOKEN);
    while(TRUE) {
        switch(t.ShowToken()) {
            default:
                if(Object::Parse(t,S)) return TRUE;    
        }
    }
}

int Poly::Parse(NavTokenizer& t,Scene& S) {
    int i;
    //initt=S.TempRef ; //FK
    t.ParseToken(LEFT_CURLY_TOKEN);
    while(TRUE) {
        switch(t.ShowToken()) {
            case POINTS_TOKEN:
                t.GetToken();
                npoint=(int)t.GetLong();
                if(point) delete point;
                point=new double[npoint*3];assert(point);
                for(i=0;i<npoint;i++) {
                    NS_REAL vector[3];
                    t.ParseToken(COMMA_TOKEN);
                    t.ParseVector(vector,3);
                    point[3*i]=vector[0];
                    point[3*i+1]=vector[1];
                    point[3*i+2]=vector[2];
                }
                break;
            case VERTICES_TOKEN:
                t.GetToken();
                nvert=(int)t.GetLong();
                if(vert) delete vert;
                vert=new int[nvert];assert(vert);
                for(i=0;i<nvert;i++) {
                    t.ParseToken(COMMA_TOKEN);
                    vert[i]=(int)t.GetLong();
                }
                break;
            default:
                if(Object::Parse(t,S)) return TRUE;    
        }
    }
}

int Zylinder::Parse(NavTokenizer& t,Scene& S) {
    initt=S.TempRef ; //FK
    t.ParseToken(LEFT_CURLY_TOKEN);
    while(TRUE) {
        switch(t.ShowToken()) {
            case X_TOKEN:
                t.GetToken();
                heading=0;
                break;
            case Y_TOKEN:
                t.GetToken();
                heading=1;
                break;
            case Z_TOKEN:
                t.GetToken();
                heading=2;
                break;
            default:
                if(Object::Parse(t,S)) return TRUE;    
        }
    }
}

const char* parse_csg_err="CSG - Only two objects allowed";
int CSG::Parse(NavTokenizer& t,Scene& S) {
    t.ParseToken(LEFT_CURLY_TOKEN);
    int have=0;
    Object* objnew;
    while(TRUE) {
        switch(t.ShowToken()) {
            case BOX_TOKEN:
                t.GetToken();
                if(have==2) t.Error(parse_csg_err);
                objnew=new Box() ;
                objnew->Parse(t,S);
                if(have==0) csg1=objnew; else csg2=objnew;
                have++;
                break;
            case SPHERE_TOKEN:
                t.GetToken();
                if(have==2) t.Error(parse_csg_err);
                objnew=new Sphere();
                objnew->Parse(t,S);
                if(have==0) csg1=objnew; else csg2=objnew;
                have++;
                break;
            case ZYLINDER_TOKEN:
                t.GetToken();
                if(have==2) t.Error(parse_csg_err);
                objnew=new Zylinder();
                objnew->Parse(t,S);
                if(have==0) csg1=objnew; else csg2=objnew;
                have++;
                break;                    
            case UNION_TOKEN:
                t.GetToken();
                if(have==2) t.Error(parse_csg_err);
                objnew=new CSG(2);
                objnew->Parse(t,S);
                if(have==0) csg1=objnew; else csg2=objnew;
                have++;
                break;                    
            case DIFFERENCE_TOKEN:
                t.GetToken();
                if(have==2) t.Error(parse_csg_err);
                objnew=new CSG(0);
                objnew->Parse(t,S);
                if(have==0) csg1=objnew; else csg2=objnew;
                have++;
                break;                    
            case INTERSECTION_TOKEN:
                t.GetToken();
                if(have==2) t.Error(parse_csg_err);
                objnew=new CSG(1);
                objnew->Parse(t,S);
                if(have==0) csg1=objnew; else csg2=objnew;
                have++;
                break;                    
            case HALFSPACE_TOKEN:
                t.GetToken();
                if(have==2) t.Error(parse_csg_err);
                objnew=new HalfSpace();
                objnew->Parse(t,S);
                if(have==0) csg1=objnew; else csg2=objnew;
                have++;
                break;                    
            case POLY_TOKEN:
                t.GetToken();
                if(have==2) t.Error(parse_csg_err);
                objnew=new Poly();
                objnew->Parse(t,S);
                ((Poly*)objnew)->Init();
                if(have==0) csg1=objnew; else csg2=objnew;
                have++;
                break;                    
            default:
                if(Object::Parse(t,S)) {
                    if(have==2) return TRUE; else
                        t.Error("CSG - not enough objects");
                }
        }
    }
}

int HalfSpace::Parse(NavTokenizer& t,Scene& S) {
    t.ParseToken(LEFT_CURLY_TOKEN);
    NS_REAL vector[3];
    while(TRUE) {
         switch(t.ShowToken()) {
            case HESSE_TOKEN:
                t.GetToken();
                t.ParseVector(vector,3);
                a=vector[0];b=vector[1];c=vector[2];
                t.ParseToken(COMMA_TOKEN);
                d=t.GetDouble();
                break;
            default:
                if(Object::Parse(t,S)) return TRUE;    
        }
    }
}


void Scene::ParseDimension(NavTokenizer& t) {
    NS_REAL v[3];
    int i;
    t.ParseToken(LEFT_CURLY_TOKEN);
    int have_length=FALSE;
    int have_dxdydz=0;    
    while(TRUE) {
        switch(t.GetToken()) {
            case LENGTH_TOKEN:
                t.ParseVector(dimension,3);
                break;
            case RES_TOKEN:
                t.ParseVector(v,3);
                gridp[0]=(int)v[0];gridp[1]=(int)v[1];gridp[2]=(int)v[2];
                have_length=TRUE;
                break;
            case X_TOKEN:
                if(!have_length) t.Error("Define 'resolution' before 'x'");
                if(d[0]) delete[] d[0];
                d[0]=new NS_REAL[gridp[0]+5];assert(d[0]);d[0]+=2;
                for(i=1;i<=gridp[0]+1;i++) {
                    t.GetDouble();
                    d[0][i]=t.GetDouble();
                }
                for(i=1;i<=gridp[0];i++) d[0][i]=d[0][i+1]-d[0][i];
                d[0][0]=d[0][1];d[0][gridp[0]+1]=d[0][gridp[0]];
                have_dxdydz|=1;
                break;
            case Y_TOKEN:
                if(!have_length) t.Error("Define 'resolution'  before 'y'");
                if(d[1]) delete[] d[1];
                d[1]=new NS_REAL[gridp[1]+5];assert(d[1]);d[1]+=2;
                for(i=1;i<=gridp[1]+1;i++) {
                    t.GetDouble();
                    d[1][i]=t.GetDouble();
                }
                for(i=1;i<=gridp[1];i++) d[1][i]=d[1][i+1]-d[1][i];
                d[1][0]=d[1][1];d[1][gridp[1]+1]=d[1][gridp[1]];        
                have_dxdydz|=2;
                break;
            case Z_TOKEN:
                if(!have_length) t.Error("'Define 'resolution' before 'z'");
                if(d[2]) delete[] d[2];
                d[2]=new NS_REAL[gridp[2]+5];assert(d[2]);d[2]+=2;
                for(i=1;i<=gridp[2]+1;i++) {
                    t.GetDouble();
                    d[2][i]=t.GetDouble();
                }
                for(i=1;i<=gridp[2];i++) d[2][i]=d[2][i+1]-d[2][i];
                d[2][0]=d[2][1];d[2][gridp[2]+1]=d[2][gridp[2]];
                have_dxdydz|=4;
                break;
	    case RIGHT_CURLY_TOKEN:
	        if((have_dxdydz&1)==0) {                    // falls keine Gitterweiten eingelesen, aequidistantes Gitter setzen
                    d[0]=new NS_REAL[gridp[0]+5];assert(d[0]);d[0]+=2;
                    for(i=0;i<gridp[0]+2;i++) d[0][i]=dimension[0]/gridp[0];
                }
                if((have_dxdydz&2)==0) {
                    d[1]=new NS_REAL[gridp[1]+5];assert(d[1]);d[1]+=2;
                    for(i=0;i<gridp[1]+2;i++) d[1][i]=dimension[1]/gridp[1];
                }
                if((have_dxdydz&4)==0) {
                    d[2]=new NS_REAL[gridp[2]+5];assert(d[2]);d[2]+=2;
                    for(i=0;i<gridp[2]+2;i++) d[2][i]=dimension[2]/gridp[2];
		}
		return;
	   default:
                t.Error(syn_err);
        }
    }
}

void Scene::ParseParameter(NavTokenizer& t) {
    int n;
    int tempparam_given=0;
    NS_REAL value;
    t.ParseToken(LEFT_CURLY_TOKEN);
    while(TRUE) {
	switch(t.GetToken()) {
	case FROUDE_TOKEN:
	    froude=t.GetDouble();
	    break;
	case TIME_DISKR_TOKEN:
	    switch(t.GetToken()){
	    case EU1_TOKEN:
		TimeDis=0;
		break;
	    case AB2_TOKEN:
		TimeDis=1;
		break;
	    default:
		t.Error("ERROR: wrong TimeDis argument, use:\n EU1 (euler), AB2 (adams-bashfort)\n");
	    }
	    break;
	case DIM_LESS_TOKEN:
	    switch(t.GetToken()){
	    case ON_TOKEN:
		dimless=1;
		break;
	    case OFF_TOKEN:
		dimless=0;
		break;
	    default:
		t.Error("ERROR: wrong dimless argument!\n");
	    }
	    break;
	case SOLVER_TOKEN:
	    switch(t.GetToken()) {
	    case SOR_TOKEN:
		Solver = 0;
		break;
	    case SSOR_TOKEN:
		Solver = 1;
		break;
	    case REDBLACK_TOKEN:
		Solver = 2;
		break;
	    case COLORSOR_TOKEN:
		Solver = 3;
		break;
	    case COLORSSOR_TOKEN:
		Solver = 4;
		break;
	    case BICGSTAB_TOKEN:
		Solver = 5;
		break;
	    case BICGSTAB_JVK_TOKEN:
		Solver = 6;
		break;
	    default:
		t.Error(syn_err);
	    }
	    break;
	case CONVECTIVE_TOKEN:
	    switch(t.GetToken()) {
	    case DC_TOKEN:
		Convective = DonorCell;
		GC = 1;
		break;
	    case QUICK_TOKEN:
		Convective = QUICK;
		GC = 2;
		break;
	    case HLPA_TOKEN:
		Convective = HLPA;
		GC = 2;
		break;
	    case SMART_TOKEN:
		Convective = SMART;
		GC = 2;
		break;
	    case VONOS_TOKEN:
		Convective = VONOS;
		GC = 2;
		break;
	    default:
		t.Error(syn_err);
	    }
	    break;
	case PRSTEP_TOKEN:
	    prstep=(int)t.GetLong();
	    break;
	case ITERMAX_TOKEN:
	    itermax=(int)t.GetLong();
	    break;
	case DELT_TOKEN:
	    deltmax=t.GetDouble();
	    break;
	case TFIN_TOKEN:
	    Tfin=t.GetDouble();
	    break;
	case EPS_TOKEN:
	    eps=t.GetDouble();
	    break;
	case OMEGA_TOKEN:
	    omg=t.GetDouble();
	    break;
	case ALPHA_TOKEN:
	    alpha=t.GetDouble();
	    break;
	case ALPHATC_TOKEN:
	    alphatg=t.GetDouble();
	    break;
	case RE_TOKEN:
	    re=t.GetDouble();
	    if(re<0.0) t.Error("ERROR: negative Reynolds number!");
	    break;
	case GX_TOKEN:
	    g[0]=t.GetDouble();
	    break;
	case GY_TOKEN:
	    g[1]=t.GetDouble();
	    break;
	case GZ_TOKEN:
	    g[2]=t.GetDouble();
	    break;
	case PRANDTL_TOKEN:
	    prandtl=t.GetDouble();
	    break;
	case BETA_TOKEN:
	    beta=t.GetDouble();
	    break;
	case TEMPCOLD_TOKEN:
	    printf("WARNING: TempCold no longer supported, set TempRef only\n") ;
	    TempRef=t.GetDouble();
	    CompTemp=TRUE;
	    tempparam_given|=1;
	    break;
	case TEMPHOT_TOKEN:
	    printf("WARNING: TempHot no longer supported, set TempRef only\n") ;
	    TempRef=t.GetDouble();
	    CompTemp=TRUE;
	    tempparam_given|=2;
	    break;
	case TEMPREF_TOKEN:
	    TempRef=t.GetDouble();
	    CompTemp=TRUE;
	    //tempparam_given|=2;
	    break;
	case TFDIFF_TOKEN:
	    tfdiff=t.GetDouble();
	    break;
	case TFCONV_TOKEN:
	    tfconv=t.GetDouble();
	    break;    
	case TARGETDIR_TOKEN:
	    strcpy(targetdir,t.GetString());
	    break;
	case PRTTIME_TOKEN:
	    prtTime=t.GetDouble();
	    break;
	case CHEMNUM_TOKEN:
	    nchem=(int)t.GetDouble();
	    if(chemc) delete chemc;
	    chemc=new NS_REAL[nchem];
	    assert(chemc);
	    for(n=0;n<nchem;n++) {
		t.ParseToken(COMMA_TOKEN);
		chemc[n]=t.GetDouble();
	    }
	    CompChem=TRUE;
	    break;
	case PERIODX_TOKEN:
	    periodbound|=1;
	    break;
	case PERIODY_TOKEN:
	    periodbound|=2;
	    break;
	case PERIODZ_TOKEN:
	    periodbound|=4;
	    break;           
	case RIGHT_CURLY_TOKEN:
	    //if (CompTemp==TRUE) if (tempparam_given!=3) t.Error("Missing initializers for temperature references");
	    //TempRef=0.5*(TempHot+TempCold);
	    return;             
	default:
	    t.Error(syn_err);
	}
    }
}

Object* NavierSetup::Parse(NavTokenizer& t) {
    Object* objlist=new Box(),*objnew;
    int have=0 ;
    while(TRUE) {
        switch (t.GetToken()) {
            case DIMENSION_TOKEN:
                if(have&1) t.Error("Multiple definition of 'dimension'");
                S.ParseDimension(t);
                have|=1;
                break;
            case PARAMS_TOKEN:
                if(have&2) t.Error("Multiple definition of 'parameter'");
                S.ParseParameter(t);
                have|=2;
                break;
            case BOX_TOKEN:
                objnew=new Box();
                objnew->Parse(t,S);
		//objnew->initt=S.TempRef ;
                objlist->Insert(objnew);
                break;
            case SPHERE_TOKEN:
                objnew=new Sphere();
                objnew->Parse(t,S);
                objlist->Insert(objnew);
                break;
            case ZYLINDER_TOKEN:
                objnew=new Zylinder();
                objnew->Parse(t,S);
                objlist->Insert(objnew);
                break;                    
            case UNION_TOKEN:
                objnew=new CSG(2);
                objnew->Parse(t,S);
                objlist->Insert(objnew);
                break;
            case DIFFERENCE_TOKEN:
                objnew=new CSG(0);
                objnew->Parse(t,S);
                objlist->Insert(objnew);
                break;
            case INTERSECTION_TOKEN:
                objnew=new CSG(1);
                objnew->Parse(t,S);
                objlist->Insert(objnew);
                break;
            case HALFSPACE_TOKEN:
                objnew=new HalfSpace();
                objnew->Parse(t,S);
                objlist->Insert(objnew);
                break;
            case POLY_TOKEN:
                objnew=new Poly();
                objnew->Parse(t,S);
                ((Poly*)objnew)->Init();
                objlist->Insert(objnew);
                break;
            case TOKEN_EOF: 
                return objlist;         
            default:
                t.Error(syn_err);
        }
    }
    return 0;
}












