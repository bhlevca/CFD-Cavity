/* NaSt3DGP - The Parallel 3D Navier-Stokes Solver
 * Copyright (C) 2003 Institute for Numerical Simulation
 *                    University of Bonn
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 */

#ifndef PARSE_HPP
#define PARSE_HPP
#include "navier.h"
#include "tokenizer.h"

//! This class contains variables and methods which are important for geometry object generation     
class Object {
protected:    
     //! minimal x-coordinate value   
     NS_REAL      xmin;
     //! maximal x-coordinate value
     NS_REAL      xmax;
     //! minimal y-coordinate value
     NS_REAL      ymin;
     //! maximal y-coordinate value
     NS_REAL      ymax;
     //! minimal z-coordinate value
     NS_REAL      zmin;
     //! maximal z-coordinate value
     NS_REAL      zmax;   
     //! tag for switching to relative coordinates
     int          is_relative;   
     //! tag for additional inflow boundaries
     int           add_inflow;        
     //! initial values for u/v/w-velocity field
     NS_REAL       initv[3];                       
     //! initial value for pressure
     NS_REAL       initp;                 
     //! initial value for temperature
     NS_REAL       initt;                     
     //! initial value for all chemicals
     NS_REAL*      initc;          
     //! fixed prescibed inflow values for u/v/w-velocity field  
     NS_REAL       inflv[3];  
     //! additional parameter for outflow boundaries to enforce compatibility equation  
     NS_REAL       inoutparam;
     //! tag for type of outflow boundary
     int           inouttype;
     int           add_cwca;                 
     //! flag helps to calculate the boundary values, TEMPB,CHEMB,OBSTACLE
     int           flag;
     //! pointer to next object
     Object*       next;
     //! pointer to previous object
     Object*       prev;

  void			Transform(NS_REAL& x,NS_REAL& y,NS_REAL& z) {}
  void			TransformBack(NS_REAL& x,NS_REAL& y,NS_REAL& z) {}

public:
                Object();
  virtual       ~Object();
  //! add to the end of the list 
  void          Insert(Object* obj);       
  //! delete all list entries
  void          DeleteAll(); 
  //! is point (int i, int j, int k) inside an given geometry
  virtual int   IsInside(int i,int j,int k,const Scene& S)=0;        
  //! is point (NS_REAL x,NS_REAL y,NS_REAL z) inside an given geometry
  virtual int	IsInside(NS_REAL x,NS_REAL y,NS_REAL z,const Scene& S)=0;
  //! gives the flag value 
  int           GetFlag() {return flag;}
  //! parse input file
  virtual int   Parse(Tokenizer& t,Scene& S); 
  void          CreatePreFlag(Matrix<unsigned>& flag,Scene& S);
  //! creates the flag field which describes the geometry
  void          CreateFlag(Matrix<unsigned>& flag,Matrix<NS_REAL>& U,Matrix<NS_REAL>& V,Matrix<NS_REAL>& W,Matrix<NS_REAL>& P,Matrix<NS_REAL>& T,
                           Matrix<NS_REAL>* CH,Scene& S);
  //void          Create_VListElt(Matrix<unsigned>& flag,Scene& S);
  int           AddInflow() {return add_inflow;}
  NS_REAL*      Inflow() {return inflv;}

};

class Box:public Object {
protected:
  int           wall;    
public:
                Box():Object() {wall=0;}
  int           IsInside(int i,int j,int k,const Scene& S);
  int		IsInside(NS_REAL x,NS_REAL y,NS_REAL z,const Scene& S);
  int           Parse(Tokenizer& t,Scene& S);
};

class Sphere:public Object {
public:
  int			IsInside(int i,int j,int k,const Scene& S);
  int			IsInside(NS_REAL x,NS_REAL y,NS_REAL z,const Scene& S);
  int			Parse(Tokenizer& t,Scene& S);
};  

class Zylinder:public Object {
protected:
  int           heading;                    // Richtung des Zylinders: X=0,Y=1,Z=2
public:
                Zylinder():Object() {heading=0;}
  int           IsInside(int i,int j,int k,const Scene& S);
  int		IsInside(NS_REAL x,NS_REAL y,NS_REAL z,const Scene& S);
  int           Parse(Tokenizer& t,Scene& S);
};

class HalfSpace:public Object {
protected:
  NS_REAL       a,b,c,d;
public:
                HalfSpace():Object() {a=b=c=d=0;}
  int           IsInside(int i,int j,int k,const Scene& S);
  int		IsInside(NS_REAL x,NS_REAL y,NS_REAL z,const Scene& S);
  int           Parse(Tokenizer& t,Scene& S);
};

class CSG:public Object {
protected:
  int           type;                       // difference=0,intersection=1,union=2
  Object*       csg1;
  Object*       csg2;
public:
                CSG(int typ):Object() {type=typ;csg1=csg2=0;}
                ~CSG();
  int           IsInside(int i,int j,int k,const Scene& S);
  int			IsInside(NS_REAL x,NS_REAL y,NS_REAL z,const Scene& S);
  int           Parse(Tokenizer& t,Scene& S);
};    



struct Reserved_Word_Struct {
   TOKEN               Token_Number;
   char               *Token_Name;
};

//! this describes a token for the parser
class               Token_Struct {
public:
     //! type of the token
   TOKEN Token_Id;
   int                 Token_Line_No;
   char               *Token_String;
   NS_REAL             Token_Float;
   TOKEN               Begin_Id;
   int                 Unget_Token, End_Of_File;
   char               *Filename;

                       Token_Struct();
   char               *Get_Token_String(TOKEN Token_Id);
   void                Write_Token(TOKEN Token_Id, DATA_FILE * Data_File, char *String);

   void                UngetToken() {Unget_Token=TRUE;}
};

//! simple wrapper around a file which serves as datafile for the parser
class               Data_File_Struct {
public:
     //! file handle
   FILE * File;
   char               *Filename;
   int                 Line_Number;

                       Data_File_Struct(char *filename);
                      ~Data_File_Struct();

		      //! get next char from file
   int                 Getc() {
      return getc(File);
   } void              Ungetc(int c) {
      ungetc(c, File);
   } int               Skip_Spaces();
   int                 Parse_C_Comments();
   void                Token_Error(char *str);

};

//! central part of the parser
class               Tokenizer {
 public:
     Token_Struct          Token;
     
 protected:
     //! the file to be parsed
     Data_File_Struct*   Data_File;
     char**              Symbol_Table;
     int                 String_Index;
     int                 Number_Of_Symbols;
     char                String[MAX_STRING_INDEX];
     int                 token_count;
     TOKEN*              Brace_Stack;
     int                 Brace_Index;
     
     void                Begin_String() {String_Index = 0;}
     int                 Find_Symbol();
     void                Stuff_Character(int c);
     void                End_String() {Stuff_Character((int)'\0');}
     int                 Read_Float();
     void                Parse_String();
     int                 Read_Symbol();
     int                 Find_Reserved();
     void                Parse_Error_Str(char *str);
     void                Parse_Error(TOKEN Token_Id);
     void                Found_Instead();
 public:
     Tokenizer(char *filename);
     ~Tokenizer();
     /*! read a vector of the form <v0, v1, ..., vn> with v0 .. vn size
      * double numbers and save them to v.
      */
     void                 ParseVector(NS_REAL* v,int size);
     void                 Warning(char* s); 
     void                 Parse_Begin();
     void                 Parse_End();
     void                 Parse_Comma();
     NS_REAL                 Parse_Float();
     void                 Get_Token();
     void                 Where_Error();
     
};

#endif
