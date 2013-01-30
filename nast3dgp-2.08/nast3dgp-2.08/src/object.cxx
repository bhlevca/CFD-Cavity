/* NaSt3DGP - The Parallel 3D Navier-Stokes Solver
 * Copyright (C) 2003 Institute for Numerical Simulation
 *                    University of Bonn
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 */

#include "object.h"
#include "typen.h"
#include <stdlib.h>
#include <stdio.h>


Object::Object() {
    int i;
    xmin=xmax=ymin=ymax=zmin=zmax=0.0;
    is_relative=FALSE;
    add_inflow=FALSE;
    for(i=0;i<3;i++) initv[i]=inflv[i]=0.0;
    for(i=0;i<3;i++) initdata[i]=0;
    initp=0.0;
    initt=0.0;
    initc=0;
    flag=TEMPB|CHEMB|OBSTACLE;
    next=prev=this;
    inouttype=0;
    inoutparam=0.0;
    tokens_given=0;
}


Object::~Object() {
    assert(next);assert(prev);
    for(int i=0;i<3;i++) if(initdata[i]) delete initdata[i];
    prev->next=next;next->prev=prev;
    if(initc) delete[] initc;
}

//
// insert *obj  in List *this right before *this (i.a. at the beginning)
void Object::Insert(Object* obj) {
    assert(prev);
    obj->next=this;
    obj->prev=prev;
    prev->next=obj;
    prev=obj; // prev==obj  end mark ..
}


void Object::DeleteAll() {
    while(prev!=this) {
      assert(prev);assert(next);
      delete next;
    }
}


NS_REAL Object::GetInflowValue(int i,int j,int k,int num,const Scene& S) {
    switch(num) {
    case 0: return inflv[0];
    case 1: return inflv[1];
    case 2: return inflv[2];
    default: return 0;
    }
}


#define XMIN (is_relative?xmin:0.5*(S.kabs[0][(int)xmin]+S.kabs[0][(int)xmin+1]))
#define XMAX (is_relative?xmax:0.5*(S.kabs[0][(int)xmax]+S.kabs[0][(int)xmax+1]))
#define YMIN (is_relative?ymin:0.5*(S.kabs[1][(int)ymin]+S.kabs[1][(int)ymin+1]))
#define YMAX (is_relative?ymax:0.5*(S.kabs[1][(int)ymax]+S.kabs[1][(int)ymax+1]))
#define ZMIN (is_relative?zmin:0.5*(S.kabs[2][(int)zmin]+S.kabs[2][(int)zmin+1]))
#define ZMAX (is_relative?zmax:0.5*(S.kabs[2][(int)zmax]+S.kabs[2][(int)zmax+1]))


int Box::IsInside(int i,int j,int k,const Scene& S) {
    switch(wall) {
        case NORTH:   return(i==S.gridp[0]+1);
        case SOUTH:   return(i==0);
        case WEST:    return(j==S.gridp[1]+1);
        case EAST:    return(j==0);
        case TOP:     return(k==S.gridp[2]+1);
        case BOTTOM:  return(k==0);
    }
    if(is_relative) {
        NS_REAL    x=0.5*(S.kabs[0][i]+S.kabs[0][i+1]);        // get position of point
        NS_REAL    y=0.5*(S.kabs[1][j]+S.kabs[1][j+1]);
        NS_REAL    z=0.5*(S.kabs[2][k]+S.kabs[2][k+1]);
        TransformBack(x,y,z);
        return( x>=xmin && x<=xmax && y>=ymin && y<=ymax && z>=zmin && z<=zmax );
    }
    return( i>=xmin && i<=xmax && j>=ymin && j<=ymax && k>=zmin && k<=zmax);
}


int Box::IsInside(NS_REAL x,NS_REAL y,NS_REAL z,const Scene& S) {
    switch(wall) {
        case NORTH: return(x>=S.dimension[0]);
        case SOUTH: return(x<=0);
        case WEST:  return(y>=S.dimension[1]);
        case EAST:  return(y<=0);
        case TOP:   return(z>=S.dimension[2]);
        case BOTTOM:return(z<=0);
    }
    TransformBack(x,y,z);
    return(x>=XMIN && x<=XMAX && y>=YMIN && y<=YMAX && z>=ZMIN && z<=ZMAX);
}


NS_REAL Box::GetInflowValue(int i,int j,int k,int num,const Scene& S) {
    assert(num>=0 && num<3);
    if(!inflowdata[num]) switch(num) {
    case 0: return inflv[0];
    case 1: return inflv[1];
    case 2: return inflv[2];
    default: return 0;
    } else {
        NS_REAL        x=(num==0)?S.kabs[0][i+1]:0.5*(S.kabs[0][i]+S.kabs[0][i+1]);
        NS_REAL        y=(num==1)?S.kabs[1][j+1]:0.5*(S.kabs[1][j]+S.kabs[1][j+1]);
        NS_REAL        z=(num==2)?S.kabs[2][k+1]:0.5*(S.kabs[2][k]+S.kabs[2][k+1]);
        switch(wall) {
        case SOUTH:
        case NORTH: return inflowdata[num]->GetValueAt(y,z);
        case WEST:
        case EAST:  return inflowdata[num]->GetValueAt(x,z);
        case TOP:
        case BOTTOM:return inflowdata[num]->GetValueAt(x,y);
        }
    }
    return 0;
}


int Sphere::IsInside(int i,int j,int k,const Scene& S) {
    NS_REAL        xm=XMIN+0.5*(XMAX-XMIN),xr=XMAX-xm;          // Mittelpunkt+Radien bestimmen
    NS_REAL        ym=YMIN+0.5*(YMAX-YMIN),yr=YMAX-ym;
    NS_REAL        zm=ZMIN+0.5*(ZMAX-ZMIN),zr=ZMAX-zm;
    NS_REAL        x=0.5*(S.kabs[0][i]+S.kabs[0][i+1]);         // Position des Punktes ermitteln
    NS_REAL        y=0.5*(S.kabs[1][j]+S.kabs[1][j+1]);
    NS_REAL        z=0.5*(S.kabs[2][k]+S.kabs[2][k+1]);
    TransformBack(x,y,z);
    return(sqr((xm-x)/xr)+sqr((ym-y)/yr)+sqr((zm-z)/zr)<=1.0);
}


int Sphere::IsInside(NS_REAL x,NS_REAL y,NS_REAL z,const Scene& S) {
    NS_REAL        xm=XMIN+0.5*(XMAX-XMIN),xr=XMAX-xm;          // Mittelpunkt+Radien bestimmen
    NS_REAL        ym=YMIN+0.5*(YMAX-YMIN),yr=YMAX-ym;
    NS_REAL        zm=ZMIN+0.5*(ZMAX-ZMIN),zr=ZMAX-zm;
    TransformBack(x,y,z);
    return(sqr((xm-x)/xr)+sqr((ym-y)/yr)+sqr((zm-z)/zr)<=1.0);
}


int Zylinder::IsInside(int i,int j,int k,const Scene& S) {
    NS_REAL        xm=XMIN+0.5*(XMAX-XMIN),xr=XMAX-xm;          // Mittelpunkt+Radien bestimmen
    NS_REAL        ym=YMIN+0.5*(YMAX-YMIN),yr=YMAX-ym;
    NS_REAL        zm=ZMIN+0.5*(ZMAX-ZMIN),zr=ZMAX-zm;
    NS_REAL        x=0.5*(S.kabs[0][i]+S.kabs[0][i+1]);         // Position des Punktes ermitteln
    NS_REAL        y=0.5*(S.kabs[1][j]+S.kabs[1][j+1]);
    NS_REAL        z=0.5*(S.kabs[2][k]+S.kabs[2][k+1]);
    TransformBack(x,y,z);
    assert(heading<3 && heading>=0);
    switch(heading) {
        case 0: if(is_relative) {if(x>xmax || x<xmin) return FALSE;}
                   else {if(i>xmax || i<xmin) return FALSE;}
                return(sqr((ym-y)/yr)+sqr((zm-z)/zr)<=1.0);
        case 1: if(is_relative) {if(y>ymax || y<ymin) return FALSE;}
                   else {if(j>ymax || j<ymin) return FALSE;}
                return(sqr((xm-x)/xr)+sqr((zm-z)/zr)<=1.0);
        case 2: if(is_relative) {if(z>zmax || z<zmin) return FALSE;}
                   else {if(k>zmax || k<zmin) return FALSE;}
                return(sqr((ym-y)/yr)+sqr((xm-x)/xr)<=1.0);
    }
    return FALSE;
}


int Zylinder::IsInside(NS_REAL x,NS_REAL y,NS_REAL z,const Scene& S) {
    NS_REAL        xm=XMIN+0.5*(XMAX-XMIN),xr=XMAX-xm;          // Mittelpunkt+Radien bestimmen
    NS_REAL        ym=YMIN+0.5*(YMAX-YMIN),yr=YMAX-ym;
    NS_REAL        zm=ZMIN+0.5*(ZMAX-ZMIN),zr=ZMAX-zm;
    assert(heading<3 && heading>=0);
    switch(heading) {
        case 0: if(x>XMAX || x<XMIN) return FALSE;
                return(sqr((ym-y)/yr)+sqr((zm-z)/zr)<=1.0);
        case 1: if(y>YMAX || y<YMIN) return FALSE;
                return(sqr((xm-x)/xr)+sqr((zm-z)/zr)<=1.0);
        case 2: if(z>ZMAX || z<ZMIN) return FALSE;
                return(sqr((ym-y)/yr)+sqr((xm-x)/xr)<=1.0);
    }
    return FALSE;
}


int HalfSpace::IsInside(NS_REAL x,NS_REAL y,NS_REAL z,const Scene& S) {
   TransformBack(x,y,z);
   return(a*x+b*y+c*z+d>=0);
}


int HalfSpace::IsInside(int i,int j,int k,const Scene& S) {
    NS_REAL        x=0.5*(S.kabs[0][i]+S.kabs[0][i+1]);    // Position des Punktes ermitteln
    NS_REAL        y=0.5*(S.kabs[1][j]+S.kabs[1][j+1]);
    NS_REAL        z=0.5*(S.kabs[2][k]+S.kabs[2][k+1]);
    return IsInside(x,y,z,S);
}


int Poly::IsInside(int i,int j,int k,const Scene& S) {
    NS_REAL        x=0.5*(S.kabs[0][i]+S.kabs[0][i+1]);    // Position des Punktes ermitteln
    NS_REAL        y=0.5*(S.kabs[1][j]+S.kabs[1][j+1]);
    NS_REAL        z=0.5*(S.kabs[2][k]+S.kabs[2][k+1]);
    return IsInside(x,y,z,S);
}


inline double Poly::DCLineLine(const double* l1,const double* l2) {
// line (l1 & l2):     x,y,z,      l,m,n
//                     0,1,2,      3,4,5
//                  St�tzvektor  Richtungsvektor

    double result=( (l1[0]-l2[0])*l1[4]*l2[5] + (l1[1]-l2[1])*l1[5]*l2[3] + (l1[2]-l2[2])*l1[3]*l2[4] -
		    (l1[2]-l2[2])*l1[4]*l2[3] - (l1[0]-l2[0])*l1[5]*l2[4] - (l1[1]-l2[1])*l1[3]*l2[5] );

    if(result!=0) return result;
    if(l2[0]-l1[0]==0 && l1[3]==0) return 1;
    if(l2[1]-l1[1]==0 && l1[4]==0) return 1;
    if(l2[2]-l1[2]==0 && l1[5]==0) return 1;
    double x;
    if(l1[3]!=0) {
        x=(l2[0]-l1[0])/l1[3];
        if(l1[4]*x!=(l2[1]-l1[1])) return 1;
        if(l1[5]*x!=(l2[2]-l1[2])) return 1;
        return 0;
    }
    if(l1[4]!=0) {
        x=(l2[1]-l1[1])/l1[4];
        if(l1[5]*x!=(l2[2]-l1[2])) return 1;
        return 0;
    }
    return 1;
}


inline double Poly::DPoint(const double* l,const double* p) {   //  ||(p-l) x (z)||^2 / ||(z)||^2
  // l[0-2] = (x,y,z)    asked point
  // p[0-2]              fixed point
  return ( ( sqr( (p[0]-l[0])*l[4]-(p[1]-l[1])*l[3] )+
	     sqr( (p[1]-l[1])*l[5]-(p[2]-l[2])*l[4] )+
	     sqr( (p[2]-l[2])*l[3]-(p[0]-l[0])*l[5] ) ) / (sqr(l[3])+sqr(l[4])+sqr(l[5])) );
}


inline int Poly::Pnpoly(const int* v, const double* xp, const double* yp, double x, double y) {
  int i=1,j=0,c=FALSE;
  while(v[i]!=-1) {
      if ((((yp[3*v[i]]<=y) && (y<yp[3*v[j]])) ||
	   ((yp[3*v[j]]<=y) && (y<yp[3*v[i]]))) &&
	  (x < (xp[3*v[j]] - xp[3*v[i]]) * (y - yp[3*v[i]]) / (yp[3*v[j]] - yp[3*v[i]]) + xp[3*v[i]]))
	c = !c;
      i++;j++;
  }
  return c;
}


inline double Poly::PlaneIntersect(const int* v,const double* l,const double* plane) {
  // returns z coordinate of plane intersection
  // FLT_MAX for no intersection
  double p[3];
  double r;
  r=plane[0]*l[3]+plane[1]*l[4]+plane[2]*l[5];   // Skalarprodukt von (0,0,1) mit plane=Normalenvektor der Ebene 
  if(r==0) return DBL_MAX;                     // parallel
  r=(plane[0]*l[0] + plane[1]*l[1] + plane[2]*l[2] + plane[3])/r; // (SP von (x,y,z) mit plane + plane[3] )/ (SP von (0,0,1) mit plane) 
  p[0]=l[0]-l[3]*r; // p0 = x - 0*r
  p[1]=l[1]-l[4]*r; // p1 = y - 0*r
  p[2]=l[2]-l[5]*r; // p2 = z - 1*r
  if(Pnpoly(v,point,point+1,p[0],p[1])) return p[2];
  return DBL_MAX;
}


int Poly::IsInside(NS_REAL x,NS_REAL y,NS_REAL z,const Scene& S) {
  TransformBack(x,y,z);
  double strahl[6]={x,y,z,0,0,1};
  double line[6];
  unsigned intersectok;
  int i,j,k=0,l;
  double* aplane=planes;
  unsigned number_of_intersections=0;
  double zc;
  
  if(x<boundbox[0] || x>boundbox[1] ||            // check if (x,y,z) is inside bounding box
     y<boundbox[2] || y>boundbox[3] || 
     z<boundbox[4] || z>boundbox[5]) return FALSE; 
  
  if(x==precalcx && y==precalcy) {                // use precalcbuffer    
    i=0;
    while(z>precalcbuf[i]) i++;
    return i&1;
  } else {                                        // create precalcbuffer
    precalcx=x;precalcy=y;
    do {
      intersectok=TRUE;                       // point and line intersections are not allowed, move test point a little bit
      for(i=0;i<npoint;i++)                   // point intersection 
	if(DPoint(strahl,point+3*i)==0) {     // check ob Strecke durch (x,y,z) und point nicht kolinear zu (0,0,1) ist 
	  intersectok=FALSE; 
	  i=npoint;
	}  
      for(i=0;i<nvert-1;i++)                  // line intersection 
	if(vert[i]!=-1 && vert[i+1]!=-1) {    // Strecke zwischen vert[i] und vert[i+1]
	  for(j=0;j<3;j++) line[j]=line[j+3]= point[3*vert[i]+j];   // St�tzvektor
	  for(j=0;j<3;j++) line[j+3]       -= point[3*vert[i+1]+j]; // Richtungsvektor
	  if(DCLineLine(strahl,line)==0) { 
	    intersectok=FALSE; 
	    i=nvert;
	  }
	}
      if(!intersectok) {if(k==1) strahl[0]+=DX[1]*1e-6; else strahl[1]+=DY[1]*1e-6; k^=1;}  // move test point a little bit
    } while(!intersectok);
    j=0;i=0;
    while(j<nvert) {                            // store z coords of intersections in precalcbuf
      zc=PlaneIntersect(vert+j,strahl,aplane);
      aplane+=4;
      j++;
      while(vert[j-1]!=-1) j++;
      if(zc!=DBL_MAX) {precalcbuf[i++]=zc; /*number_of_intersections++;*/}
    }
    precalcbuf[i++]=DBL_MAX;
    for(j=0;j<i-1;j++) {                          // sort precalcbuf    
      l=j;
      for(k=j+1;k<i;k++) if(precalcbuf[l]>precalcbuf[k]) l=k;
      zc=precalcbuf[j]; precalcbuf[j]=precalcbuf[l]; precalcbuf[l]=zc;
    }
    return (number_of_intersections & 1);
  }
  return 0;
}


void Poly::Init() {
    int i,j;
    int* v;
    boundbox[0]=boundbox[2]=boundbox[4]=FLT_MAX;
    boundbox[1]=boundbox[3]=boundbox[5]=FLT_MIN;
    for(i=0;i<npoint;i++) {                        // bounding box of point-cloud
        boundbox[0]=std::min(boundbox[0],point[3*i]);   // std::min x coordinate
        boundbox[1]=std::max(boundbox[1],point[3*i]);   // std::max x coordinate
        boundbox[2]=std::min(boundbox[2],point[3*i+1]); // std::min y coordinate
        boundbox[3]=std::max(boundbox[3],point[3*i+1]); // std::max y coordinate
        boundbox[4]=std::min(boundbox[4],point[3*i+2]); // std::min z coordinate
        boundbox[5]=std::max(boundbox[5],point[3*i+2]); // std::max z coordinate
    }
    i=0;
    for(j=0;j<nvert;j++) if(vert[j]==-1) i++;      // i = number of polygons 
    planes=new double[4*i];       assert(planes);  // one ploygon defines one plane
    precalcbuf=new double[nvert]; assert(precalcbuf);
    precalcx=FLT_MAX;
    v=vert;
    for(j=0;j<i;j++) {    // planes[0-2] = (v2-v0) x (v1-v0) Normalenvektor zur Ebene durch v0 v1 v2 
        planes[4*j+0]=(point[3*v[2]+1]-point[3*v[0]+1]) * (point[3*v[1]+2]-point[3*v[0]+2]) -
                      (point[3*v[2]+2]-point[3*v[0]+2]) * (point[3*v[1]+1]-point[3*v[0]+1]);
        planes[4*j+1]=(point[3*v[2]+2]-point[3*v[0]+2]) * (point[3*v[1]+0]-point[3*v[0]+0])-
                      (point[3*v[2]+0]-point[3*v[0]+0]) * (point[3*v[1]+2]-point[3*v[0]+2]);
        planes[4*j+2]=(point[3*v[2]+0]-point[3*v[0]+0]) * (point[3*v[1]+1]-point[3*v[0]+1])-
                      (point[3*v[2]+1]-point[3*v[0]+1]) * (point[3*v[1]+0]-point[3*v[0]+0]);
        planes[4*j+3]=-(planes[4*j+0]*point[3*v[0]] + planes[4*j+1]*point[3*v[0]+1] + planes[4*j+2]*point[3*v[0]+2]); // -<n,v0> 
        while((*v)!=-1) v++;   // begin new polygon
        v++;   
    }
}


int CSG::IsInside(int i,int j,int k,const Scene& S) {
    assert(type>=0 && type<3);
    assert(csg1);assert(csg2);
    switch(type){
        case 0: return((csg1->IsInside(i,j,k,S)) && !(csg2->IsInside(i,j,k,S)));
        case 1: return((csg1->IsInside(i,j,k,S)) &&  (csg2->IsInside(i,j,k,S)));
        case 2: return((csg1->IsInside(i,j,k,S)) ||  (csg2->IsInside(i,j,k,S)));
    }
    return FALSE;
}


int CSG::IsInside(NS_REAL x,NS_REAL y,NS_REAL z,const Scene& S) {
    assert(type>=0 && type<3);
    assert(csg1);assert(csg2);
    switch(type){
        case 0: return((csg1->IsInside(x,y,z,S)) && !(csg2->IsInside(x,y,z,S)));
        case 1: return((csg1->IsInside(x,y,z,S)) &&  (csg2->IsInside(x,y,z,S)));
        case 2: return((csg1->IsInside(x,y,z,S)) ||  (csg2->IsInside(x,y,z,S)));
    }
    return FALSE;
}

CSG::~CSG() {
    if(csg1) {csg1->DeleteAll();delete csg1;}
    if(csg2) {csg2->DeleteAll();delete csg2;}
}

//------------------------------------------------------------------------------------------------------------

//
// flag(i,j,k)  =   number of the last object O such that (i,j,k) \in O
// if (i,j,k) is fluid, set also the flagtypemsb bit (bit 31)
//
void Object::CreatePreFlag(Matrix<flagtype>& flag,const Scene& S) {
    int i,j,k,c;
    Object* obj;
    for(i=0;i<=S.gridp[0]+1;i++) {
        printf("\r%4d",i);fflush(stdout);
        for(j=0;j<=S.gridp[1]+1;j++) 
	  for(k=0;k<=S.gridp[2]+1;k++) {
            c=0; 
            obj=this; 
            flag[i][j][k]=0;
            // find last object O (number c) such that (i,j,k) \in O 
            while(obj->next!=this) {
	      c++; 
              obj=obj->next;
	      if(obj->IsInside(i,j,k,S)) flag[i][j][k]=c | ((obj->GetFlag()!=0)?flagtypemsb:0);
            }
        }
    }
    printf("\r        \r");
}


#define isfluid(object_id) (object_id & flagtypemsb?0:1)
#define nofluid(object_id) (object_id & flagtypemsb?1:0)

void NavierSetup::CheckFlag(Matrix<flagtype>& flag) {
    int i,j,k;
    int error;
    int errcount=0;
    do {
        error=0;
        IJKLOOP if(nofluid(flag[i][j][k])) {
            if(isfluid(flag[i-1][j][k]) && isfluid(flag[i+1][j][k])) {
	        printf("flag error x at (%2d,%2d,%2d)  ...corrected \n",i,j,k);fflush(stdout);
	        error=1; errcount++;
		flag[i][j][k]=flag[i-1][j][k]; 
            }
            if(isfluid(flag[i][j-1][k]) && isfluid(flag[i][j+1][k])) {
	        printf("flag error y at (%2d,%2d,%2d)  ...corrected \n",i,j,k);fflush(stdout);
	        error=1; errcount++;
                flag[i][j][k]=flag[i][j-1][k];
            }
            if(isfluid(flag[i][j][k-1]) && isfluid(flag[i][j][k+1])) {
                printf("flag error z at (%2d,%2d,%2d)  ...corrected \n",i,j,k);fflush(stdout);
	        error=1; errcount++; 
                flag[i][j][k]=flag[i][j][k-1]; 
            }
        }
    } while(error);
    if(errcount>0) printf("\n %d flag field errors corrected.\n",errcount);

    IJKLOOP if(nofluid(flag[i][j][k])) {    
        if(isfluid(flag[i-1][j][k]) && isfluid(flag[i+1][j][k])) printf("flag error x at (%2d,%2d,%2d)\n",i,j,k);fflush(stdout);
        if(isfluid(flag[i][j-1][k]) && isfluid(flag[i][j+1][k])) printf("flag error y at (%2d,%2d,%2d)\n",i,j,k);fflush(stdout);
        if(isfluid(flag[i][j][k-1]) && isfluid(flag[i][j][k+1])) printf("flag error z at (%2d,%2d,%2d)\n",i,j,k);fflush(stdout);
    }
}


void Object::CreateInOutFlowList(Matrix<flagtype>& flag,Scene& S) {     // create inout and inflow lists
    int     i,j,k,l;
    Object* obj;
  
    IJKALLLOOP if(flag[i][j][k]) {
        obj=this;
        // find (number of last) object in which (i,j,k) is located
        int objc=(int)(flag[i][j][k]& ~flagtypemsb) ;
        for(l=0; l<objc; l++) obj=obj->next;

        if(obj->add_inflow) {// create INFLOW list
            NS_List<INFLOW>* p;
            p=new NS_List<INFLOW>;
            p->data.i=i;p->data.j=j;p->data.k=k;
            p->data.u=obj->GetInflowValue(i,j,k,0,S);
            p->data.v=obj->GetInflowValue(i,j,k,1,S);
            p->data.w=obj->GetInflowValue(i,j,k,2,S);
            if(S.CompTemp) p->data.t=obj->initt;
            if(S.CompChem && S.nchem>0) {
                p->data.nchem=S.nchem;
                p->data.chem=new NS_REAL[S.nchem];
                for(int l=0;l<S.nchem;l++) p->data.chem[l]=(obj->initc)?obj->initc[l]:0.0;
            }
            S.InflowList.InsertFirst(p);
        }
        if(obj->inouttype) {// create INOUT list
            NS_List<INOUTREC>* p;
            p=new NS_List<INOUTREC>;
            p->data.i=i;p->data.j=j;p->data.k=k;
            p->data.type=obj->inouttype;  // type of boundary cond. 1=neumann zero, 2=convective time dep. bound. cond.
            p->data.param=obj->inoutparam;// extra parameter for outflow bound. cond. (under construction) 
            S.InOutList.InsertFirst(p);
        }
    }
}


int Object::CountInnerObstacles(Matrix<flagtype>& flag,const Scene& S) {
    int i,j,k,count=0;
    IJKLOOP if(flag[i][j][k]&flagtypemsb) count++;
    return count;
}

  
void Object::CreateFlag(Matrix<flagtype>& tempflag,Matrix<flagtype>& flag,Matrix<NS_REAL>& U,
                        Matrix<NS_REAL>& V,Matrix<NS_REAL>& W,Matrix<NS_REAL>& P,Matrix<NS_REAL>& T,
                        Matrix<NS_REAL>* CH,const Scene& S,int slice) {
    int i=slice,j,k,c,n;
    Object* obj;
    for(j=0;j<=S.gridp[1]+1;j++) 
      for(k=0;k<=S.gridp[2]+1;k++) {
	if(tempflag[i][j][k]) {
          // find last object O such that (i,j,k) \in O 
	  obj=this;
          int objc=(int)(tempflag[i][j][k]&~flagtypemsb) ;
          for(c=0; c<objc; c++) obj=obj->next;

	  flag[i][j][k]=obj->GetFlag();
	  if((i==0 || i==S.gridp[0]+1) && !(S.periodbound&1) && !(flag[i][j][k]&OBSTACLE)) 
	    NavierSetup::Error(ERR_DEFAULT,"\nError: Fluid object on ghost cells while non-periodic boundary condition in x-direction!\n");
	  if((j==0 || j==S.gridp[1]+1) && !(S.periodbound&2) && !(flag[i][j][k]&OBSTACLE)) 
	    NavierSetup::Error(ERR_DEFAULT,"\nError: Fluid object on ghost cells while non-periodic boundary condition in y-direction!\n");
	  if((k==0 || k==S.gridp[2]+1) && !(S.periodbound&4) && !(flag[i][j][k]&OBSTACLE)) 
	    NavierSetup::Error(ERR_DEFAULT,"\nError: Fluid object on ghost cells while non-periodic boundary condition in z-direction!\n");
	  if(tempflag[i][j][k]&flagtypemsb) {
	    if(i<=S.gridp[0])   if(!(tempflag[i+1][j][k]&flagtypemsb)) flag[i][j][k]|=NORTH;
	    if(i>=1)            if(!(tempflag[i-1][j][k]&flagtypemsb)) flag[i][j][k]|=SOUTH;
	    if(j<=S.gridp[1])   if(!(tempflag[i][j+1][k]&flagtypemsb)) flag[i][j][k]|=WEST;
	    if(j>=1)            if(!(tempflag[i][j-1][k]&flagtypemsb)) flag[i][j][k]|=EAST;
	    if(k<=S.gridp[2])   if(!(tempflag[i][j][k+1]&flagtypemsb)) flag[i][j][k]|=TOP;
	    if(k>=1)            if(!(tempflag[i][j][k-1]&flagtypemsb)) flag[i][j][k]|=BOTTOM;
	  }
	  unsigned flg=FREECELLS(flag[i][j][k]);
	  
	  if((flag[i][j][k]&INOUT) && (flg!=SOUTH) && (flg!=NORTH) && (flg!=BOTTOM) && (flg!=TOP) && (flg!=EAST) && (flg!=WEST) && (flg!=0)) 
	      NavierSetup::Error(ERR_DEFAULT,"\nError: INOUT b.c. and cells with faces in multiple directions not allowed.\n");
	  
	  U[i][j][k]=(!obj->add_inflow)?obj->GetInitValue(0,0.5*(S.kabs[1][j]+S.kabs[1][j+1]),0.5*(S.kabs[2][k]+S.kabs[2][k+1])):obj->GetInflowValue(i,j,k,0,S); 
	  V[i][j][k]=(!obj->add_inflow)?obj->GetInitValue(1,0.5*(S.kabs[0][i]+S.kabs[0][i+1]),0.5*(S.kabs[2][k]+S.kabs[2][k+1])):obj->GetInflowValue(i,j,k,1,S);
	  W[i][j][k]=(!obj->add_inflow)?obj->GetInitValue(2,0.5*(S.kabs[0][i]+S.kabs[0][i+1]),0.5*(S.kabs[1][j]+S.kabs[1][j+1])):obj->GetInflowValue(i,j,k,2,S);
	  P[i][j][k]=obj->initp;

	  if(S.CompTemp){
            T[i][j][k]=obj->initt ;
	    // if((obj->initt != S.TempHot) && (obj->initt != S.TempCold)) obj->initt=S.TempRef; //  NEW  !!!
	    // T[i][j][k]=(obj->initt-S.TempRef)/(S.TempHot-S.TempCold);
	  }	  	  
	  if(S.CompChem) for(n=0;n<S.nchem;n++) CH[n][i][j][k]=(obj->initc)?obj->initc[n]:0;
	} // if (tempflag) end
      }
}


void NavierSetup::MarkBound(Matrix<flagtype>& flag) {
    int i,j,k;
    
    if((S.periodbound&4)==0) IALLLOOP JALLLOOP {
        if(!flag[i][j][0])              flag[i][j][0]=flagtypemsb;
        if(!flag[i][j][S.gridp[2]+1])   flag[i][j][S.gridp[2]+1]=flagtypemsb;
    }
    if((S.periodbound&2)==0) IALLLOOP KALLLOOP {
        if(!flag[i][0][k])              flag[i][0][k]=flagtypemsb;
        if(!flag[i][S.gridp[1]+1][k])   flag[i][S.gridp[1]+1][k]=flagtypemsb;
    }
    if((S.periodbound&1)==0) JALLLOOP KALLLOOP {
        if(!flag[0][j][k])              flag[0][j][k]=flagtypemsb;
        if(!flag[S.gridp[0]+1][j][k])   flag[S.gridp[0]+1][j][k]=flagtypemsb;
    }
}


NS_REAL InflowData::GetValueAt(NS_REAL x,NS_REAL y) {
    assert(spaces[0]);assert(spaces[1]);
    if(x<spaces[0][0] || x>spaces[0][dims[0]-1] || y<spaces[1][0] || y>spaces[1][dims[1]-1]) return 0;
    int ip=0,jp=0;
    while(spaces[0][ip]<=x) ip++;
    while(spaces[1][jp]<=y) jp++;
    NS_REAL d=data[ip-1][jp-1]*(spaces[0][ip  ]-x)*(spaces[1][jp  ]-y)+
              data[ip  ][jp-1]*(x-spaces[0][ip-1])*(spaces[1][jp  ]-y)+
              data[ip  ][jp  ]*(x-spaces[0][ip-1])*(y-spaces[1][jp-1])+
              data[ip-1][jp  ]*(spaces[0][ip  ]-x)*(y-spaces[1][jp-1]);
    d/=(spaces[0][ip]-spaces[0][ip-1])*(spaces[1][jp]-spaces[1][jp-1]);
    return d;
}


InflowData::InflowData(char* filename) {
    static const char* errstring="Error in inflow data file.\n";
    InflowData();
    FILE* f=fopen(filename,"r");
    if(!f) NavierSetup::Error(ERR_OPEN,filename);
    printf("Reading inflow data file %s...",filename);fflush(stdout);
    if(fscanf(f,"%d\n%d\n",&dims[0],&dims[1])<0) NavierSetup::Error(ERR_DEFAULT,errstring);
    if(dims[0]<=0 || dims[1]<=0) NavierSetup::Error(ERR_DEFAULT,errstring);
    int planeinit[4]={0,dims[0]-1,0,dims[1]-1};
    data.Init(planeinit);
    spaces[0]=new double[dims[0]];spaces[1]=new double[dims[1]];
    int i,j;
    for(i=0;i<dims[0];i++) if(fscanf(f,"%le",&spaces[0][i])<0) NavierSetup::Error(ERR_DEFAULT,errstring);
    fscanf(f,"\n");
    for(i=0;i<dims[1];i++) if(fscanf(f,"%le",&spaces[1][i])<0) NavierSetup::Error(ERR_DEFAULT,errstring);
    fscanf(f,"\n");
    for(i=0;i<dims[0];i++) {
        for(j=0;j<dims[1];j++) if(fscanf(f,"%le",&data[i][j])<0) NavierSetup::Error(ERR_DEFAULT,errstring);
        if(i<dims[0]-1) fscanf(f,"\n");
    }
    fclose(f);
    printf("ok.\n");
}
