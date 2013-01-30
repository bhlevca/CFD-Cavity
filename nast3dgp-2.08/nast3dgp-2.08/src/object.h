/* NaSt3DGP - The Parallel 3D Navier-Stokes Solver
 * Copyright (C) 2003 Institute for Numerical Simulation
 *                    University of Bonn
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 */

#ifndef OBJECT_HPP
#define OBJECT_HPP
#include "navier.h"
#include "tokenizer.h"

//! provide inflow data at boundaries
class InflowData {
protected:
    int             dims[2];
    double*         spaces[2];
    Plane<double>     data;
public:
                    InflowData() {spaces[0]=spaces[1]=0;}
                    InflowData(char* filename);
                   ~InflowData() {if(spaces[0]) delete spaces[0];if(spaces[1]) delete spaces[1];} 
    NS_REAL            GetValueAt(NS_REAL x,NS_REAL y);
};

//! parser for nav scene files

class NavTokenizer:public Tokenizer {
 public:
     //! read a vector included in angles (<,>)
     /*!
       parses a vector which is included in angles.
       \param int num number of vector elements to parse
       \param double* v pointer to the vector which contains the result
      */
    void    ParseVector(double* val,int num);
    //! try to parse a token of type tokenid
    void    ParseToken(int tokenid);
    NavTokenizer(const char* c);
    //! get the next token. If token is a comment, continue parsing, else return token
    int     GetToken();
    int     ShowToken();
};

//! provides basic interface for geometric objects
class Object {
protected:    
     //! coordinates
    NS_REAL       xmin,xmax,ymin,ymax,zmin,zmax;      
    //! relative coordinates
    unsigned      is_relative;                        
    //! Inflow?
    unsigned      add_inflow;                         
    //! u/v/w init-values
    NS_REAL       initv[3];                           
    //! p init-value
    NS_REAL       initp;                              
    //! t init-value
    NS_REAL       initt;                              
    //! c init-values
    NS_REAL*      initc;             
    //! u/v/w inflow-values
    NS_REAL       inflv[3];                           
    NS_REAL       inoutparam;
    int           inouttype;
    //! flags for computation of boundary-values, TEMPB,CHEMB,OBSTACLE  
    unsigned      flag;                               
    Object*       next;
    Object*       prev;
    unsigned      tokens_given;  
    InflowData*   initdata[3];

    void          Transform(NS_REAL& x,NS_REAL& y,NS_REAL& z) {}   // TODO: rotation
    void          TransformBack(NS_REAL& x,NS_REAL& y,NS_REAL& z) {}

public:
                     Object();
    virtual         ~Object();
    //! insert at end of list
    void             Insert(Object* obj);                         
    //! delete whole list, except myself
    void             DeleteAll();       
    //! is the point inside
    virtual int      IsInside(int i,int j,int k,const Scene& S)=0;
    virtual int      IsInside(NS_REAL x,NS_REAL y,NS_REAL z,const Scene& S)=0;
    unsigned         GetFlag() {return flag;}
    //! parse the input file
    virtual int      Parse(NavTokenizer& t,Scene& S);              
    virtual NS_REAL  GetInflowValue(int i,int j,int k,int num,const Scene& S);
    void             CreatePreFlag(Matrix<flagtype>& flag,const Scene& S);  
    void             CreateFlag   (Matrix<flagtype>& tempflag,Matrix<unsigned>& flag,Matrix<NS_REAL>& U,
                                   Matrix<NS_REAL>& V,Matrix<NS_REAL>& W,Matrix<NS_REAL>& P,Matrix<NS_REAL>& T,
                                   Matrix<NS_REAL>* CH,const Scene& S,int slice);

    void             CreateInOutFlowList(Matrix<flagtype>& flag,Scene& S);
    int              CountInnerObstacles(Matrix<flagtype>& flag,const Scene& S);  
    int              AddInflow() {return add_inflow;}
    NS_REAL*         Inflow() {return inflv;}
    NS_REAL          GetInitValue(int i,NS_REAL px,NS_REAL py) {
                      if(initdata[i]) return initdata[i]->GetValueAt(px,py); 
                        else return initv[i];
                  }
};

//! rectangular box
class Box:public Object {
protected:
  int           wall;    
  InflowData*   inflowdata[3];

public:
                Box():Object() {wall=0;for(int i=0;i<3;i++) inflowdata[i]=0;}
               ~Box() {for(int i=0;i<3;i++) if(inflowdata[i]) delete inflowdata[i];}
  int           IsInside(int i,int j,int k,const Scene& S);
  NS_REAL       GetInflowValue(int i,int j,int k,int num,const Scene& S);
  int           IsInside(NS_REAL x,NS_REAL y,NS_REAL z,const Scene& S);
  int           Parse(NavTokenizer& t,Scene& S);
};

class Sphere:public Object {
public:
  int           IsInside(int i,int j,int k,const Scene& S);
  int           IsInside(NS_REAL x,NS_REAL y,NS_REAL z,const Scene& S);
  int           Parse(NavTokenizer& t,Scene& S);
};  

class Zylinder:public Object {
protected:
  int           heading;                    // Richtung des Zylinders: X=0,Y=1,Z=2
public:
                Zylinder():Object() {heading=0;}
  int           IsInside(int i,int j,int k,const Scene& S);
  int           IsInside(NS_REAL x,NS_REAL y,NS_REAL z,const Scene& S);
  int           Parse(NavTokenizer& t,Scene& S);
};

class HalfSpace:public Object {
protected:
  NS_REAL          a,b,c,d;
public:
                HalfSpace():Object() {a=b=c=d=0;}
  int           IsInside(int i,int j,int k,const Scene& S);
  int           IsInside(NS_REAL x,NS_REAL y,NS_REAL z,const Scene& S);
  int           Parse(NavTokenizer& t,Scene& S);
};

//! poly object
class Poly:public Object {
protected:
     //! field for points
     double*       point;              
     //! bounding box of point-cloud
     double        boundbox[6];        
     //! number of points
     int           npoint;             
     //! field for vertices
     int*          vert;              
     //! number of vertices
     int           nvert;        
     //! field for planes
     double*       planes;       
     double*       precalcbuf;
  double        precalcx,precalcy;

  double        DCLineLine(const double* l1,const double* l2);       
  double        DPoint(const double* l,const double* p);          
  //! returns z coordinate of intersection point
  double        PlaneIntersect(const int* v,const double* l,const double* plane); 
                                                
  int           Pnpoly(const int* v, const double* xp, const double* yp, double x, double y);

public:
                Poly():Object() {point=0;vert=0;npoint=nvert=0;planes=0;}
               ~Poly() {if(point) delete point;if(vert) delete vert;if(planes) delete planes;if(precalcbuf) delete precalcbuf;}
  int           IsInside(int i,int j,int k,const Scene& S);
  int           IsInside(NS_REAL x,NS_REAL y,NS_REAL z,const Scene& S);
  int           Parse(NavTokenizer& t,Scene& S);
  void          Init();
};

//! CSG operations on objects
class CSG:public Object {
protected:
     //! difference=0, intersection=1, union=2
  int           type;         
  Object*       csg1;
  Object*       csg2;
public:
                CSG(int typ):Object() {type=typ;csg1=csg2=0;}
                ~CSG();
  int           IsInside(int i,int j,int k,const Scene& S);
  int           IsInside(NS_REAL x,NS_REAL y,NS_REAL z,const Scene& S);
  int           Parse(NavTokenizer& t,Scene& S);
};    

#endif
