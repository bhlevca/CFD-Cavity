/* NaSt3DGP - The Parallel 3D Navier-Stokes Solver
 * Copyright (C) 2003 Institute for Numerical Simulation
 *                    University of Bonn
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 */

#ifndef MATRIX_INCLUDED
#define MATRIX_INCLUDED
#include <stdio.h>
#include <assert.h>

#ifndef FALSE
#define FALSE (1==0)
#endif
#ifndef TRUE
#define TRUE (0==0)
#endif

// 
// vector/plane/matrix with distributed file access functionality
//

#define INBOUNDS(a,b,i) ((a[i]>=b[i]) && (a[i]<=b[i+1]) && (a[i+1]>=b[i]) && (a[i+1]<=b[i+1]) && (a[i]<=a[i+1]) && (b[i]<=b[i+1]))


class ParParams;
template<class T> class Matrix;

//! A vector class
/*!
  This is a template class for a vector-type. This vector is also used to build the Plane and Matrix datatypes.
  The class also provides I/O-routines.
 */
template < class T >
class Vector {
 public:
     Vector() {data=0;size[0]=0;}
     int         Init(int* sz) {
	  data+=size[0];
	  if(data) delete[] data;
	  assert(sz[0]<=sz[1]);
	  size[0]=sz[0];size[1]=sz[1];
	  data=new T[size[1]-size[0]+1];
	  if(data==0) return FALSE;
	  data-=size[0];
	  for(int i=size[0];i<=size[1];i++) data[i]=0;
	  return TRUE;
     }
     Vector(int* sz) {
	  data=0;
	  size[0]=0;
	  Init(sz);
     }
     //! access operator
     /*!
       This operator is the standard way to access a vector element.
       \param int i index into the vector
       \return reference to the data type T
      */
     inline T&   operator[] (int i) const {
#ifndef NDEBUG
	  if(!((data + size[0]) && (i >= size[0]) && (i <= size[1]))) {
	       printf("Matrix access at ?,?,%d (%d/%d)\n",i,size[0],size[1]);
	       //assert((data + size[0]) && (i >= size[0]) && (i <= size[1]));
	  }
#endif
	  return data[i];
     }
     //! output routine
     /*!
       This method writes the vector to a file
       \param FILE*f file pointer
       \return bool indicating success
      */
     int         Write(FILE* f) {
	  if(fwrite(size,sizeof(int),2,f)!=2) return FALSE;
	  if(fwrite(data+size[0],size[1]-size[0]+1,sizeof(T),f)!=sizeof(T)) return FALSE;
	  return TRUE;
     }
     
     int         WriteParallel(FILE* f,
			       int(*writevect)(FILE* f,int,int,int*,const ParParams*,const Matrix<T>&,const int*),
			       const ParParams* par,int ip,int jp,const Matrix<T>& mat,
			       const int* dims , int *bigsize) {
	  if(fwrite(bigsize,sizeof(int),2,f)!=2) return FALSE;
	  if(!writevect(f,ip,jp, bigsize ,par,mat,dims)) return FALSE;
	  return TRUE;
     }
     
     static int  WriteFrame(FILE* f,int* sz) {
	  if(fwrite(sz,sizeof(int),2,f)!=2) return FALSE;
	  T* mydata=new T[sz[1]-sz[0]+1];
	  //assert(mydata);
	  for(int i=0;i<sz[1]-sz[0]+1;i++) mydata[i]=0;
	  if(fwrite(mydata,sz[1]-sz[0]+1,sizeof(T),f)!=sizeof(T)) {delete mydata;return FALSE;}
	  delete mydata;
	  return TRUE;
     }
     int         Read(FILE* f) {
	  int mysize[2];
	  if(fread(mysize,sizeof(int),2,f)!=2) return FALSE;
	  //assert(mysize[0]==size[0] && mysize[1]==size[1]);
	  if(fread(data+size[0],size[1]-size[0]+1,sizeof(T),f)!=sizeof(T)) return FALSE;
	  return TRUE;
     }
     //! write part of a big vector with size bigsize in file
     int         WritePartial(FILE* f) {                     
	  int bigsize[2];
	  if(fread(bigsize,sizeof(int),2,f)!=2) return FALSE;
	  //assert(INBOUNDS(size,bigsize,0));
	  long filepos=ftell(f);     
	  if(fseek(f,(long)sizeof(T)*(size[0]-bigsize[0]),SEEK_CUR)!=0) return FALSE;
	  if(fwrite(data+size[0],size[1]-size[0]+1,sizeof(T),f)!=sizeof(T)) return FALSE;
	  if(fseek(f,filepos+(long)sizeof(T)*(bigsize[1]-bigsize[0]+1),SEEK_SET)!=0) return FALSE;
	  return TRUE;
     }
     int         WritePartial(FILE* f,const int* sz) {       // write part of a big vector with size bigsize (in file)
	  int bigsize[2];
	  if(fread(bigsize,sizeof(int),2,f)!=2) return FALSE;
	  //assert(INBOUNDS(sz,bigsize,0));
	  long filepos=ftell(f);     
	  if(fseek(f,(long)sizeof(T)*(sz[0]-bigsize[0]),SEEK_CUR)!=0) return FALSE;
	  if(fwrite(data+sz[0],sz[1]-sz[0]+1,sizeof(T),f)!=sizeof(T)) return FALSE;
	  if(fseek(f,filepos+(long)sizeof(T)*(bigsize[1]-bigsize[0]+1),SEEK_SET)!=0) return FALSE;
	  return TRUE;
     }
     int         ReadPartial(FILE* f) {
	  int bigsize[2];
	  if(fread(bigsize,sizeof(int),2,f)!=2) return FALSE;
	  //assert(INBOUNDS(size,bigsize,0));
	  long filepos=ftell(f);
	  if(fseek(f,(long)sizeof(T)*(size[0]-bigsize[0]),SEEK_CUR)!=0) return FALSE;
	  if(fread(data+size[0],size[1]-size[0]+1,sizeof(T),f)!=sizeof(T)) return FALSE;
	  if(fseek(f,filepos+(long)sizeof(T)*(bigsize[1]-bigsize[0]+1),SEEK_SET)!=0) return FALSE;
	  return TRUE;
     }           
     int         ReadPartial(FILE* f,int* sz) {                     // write part of a big vector with size bigsize (in file)
	  int bigsize[2];
	  if(fread(bigsize,sizeof(int),2,f)!=2) return FALSE;
	  //assert(INBOUNDS(sz,size,0));
	  //assert(INBOUNDS(sz,bigsize,0));
	  long filepos=ftell(f);     
	  if(fseek(f,(long)sizeof(T)*(sz[0]-bigsize[0]),SEEK_CUR)!=0) return FALSE;
	  if(fread(data+sz[0],sz[1]-sz[0]+1,sizeof(T),f)!=sizeof(T)) return FALSE;
	  if(fseek(f,filepos+(long)sizeof(T)*(bigsize[1]-bigsize[0]+1),SEEK_SET)!=0) return FALSE;
	  return TRUE;
     }
     //! get the minimum entry
     /*!
       This method returns the minimum entry of the vector
       \return T minimum value of the vector
      */
     T           Min() {
	  T min=data[size[1]],minc;
	  for(int i=size[0];i<size[1];i++) if(min>(minc=data[i])) min=minc;
	  return min;
     }
     //! get the maximum entry
     /*!
       This method returns the maximum entry of the vector
       \return T maximal value of the vector
     */
     T           Max() {
	  T max=data[size[1]],maxc;
	  for(int i=size[0];i<size[1];i++) if(max<(maxc=data[i])) max=maxc;
	  return max;
     }
     //! determine the length of the vector
     /*!
       This method gives the length of the vector
       \return int length of vector
      */
     int         Size(int* bigsize) {return(2*sizeof(int)+sizeof(T)*(bigsize[1]-bigsize[0]+1));}
     
     ~Vector() {
	  data+=size[0];
	  if(data) delete[] data;
     }
 protected:
     //! pointer to data of type T
     T*          data;
     //! array to store index range
     int         size[2];
     int         kl, kr;
};

//! Plane datatype
/*!
  This class is an array with two indices. It is built from the vector class.
  It provides I/O-routines and a standard C-style access-operator.
  This dataype is useful for handling of boundary data of a 3D-field.
 */
template < class T >
class Plane {
 public:
     //! default constructor
     Plane() {
	  data=0;
	  size[0]=0;
     }
     int         Init(int* sz) {
	  data+=size[0];
	  if(data) delete[] data;
	  //assert(sz[0]<=sz[1] && sz[2]<=sz[3]);
	  size[0]=sz[0];size[1]=sz[1];
	  size[2]=sz[2];size[3]=sz[3];
	  data=new Vector<T>[size[1]-size[0]+1];
	  if (data==0) return FALSE;
	  data-=size[0];
	  for(int i=size[0];i<=size[1];i++) data[i].Init(&(size[2]));
	  return TRUE;
     }
     Plane(int* sz) {
	  data=0;
	  size[0]=0;
	  Init(sz);
     }
     //! access operator
     /*!
       This access operator provides access to the vectors of which the Plane is built.
       Together with the access operator of the vector class this gives a C-style indexing into the 2D-field.
       \param int i index to vector
       \return reference to vector
      */
     inline Vector<T>& operator[] (int i) const{
#ifndef NDEBUG
	  if(!((data + size[0]) && (i >= size[0]) && (i <= size[1]))) {
	       printf("Matrix access at ?,%d,? (%d/%d)\n",i,size[0],size[1]);
	       //assert((data + size[0]) && (i >= size[0]) && (i <= size[1]));
	  }
#endif
	  return data[i];
     }
     int         Write(FILE* f) {
	  if(fwrite(size,sizeof(int),4,f)!=4) return FALSE;
	  for(int i=size[0];i<=size[1];i++) if(!data[i].Write(f)) return FALSE;
	  return TRUE;
     }
     
     int         WriteParallel(FILE* f,
			       int(*writevect)(FILE* f,int,int,int*,const ParParams*,const Matrix<T>&,const int*),
			       const ParParams* par,int ip,const Matrix<T>& mat,
			       const int* dims , int *bigsize) {
	  if(fwrite(bigsize,sizeof(int),4,f)!=4) return FALSE;
	  for(int i=bigsize[0];i<=bigsize[1];i++) 
	       if(!data[i].WriteParallel(f,writevect,par,ip,i,mat,dims , &bigsize[2])) return FALSE;
	  return TRUE;
     }
     
     
     static int  WriteFrame(FILE* f,int* sz) {
	  if(fwrite(sz,sizeof(int),4,f)!=4) return FALSE;
	  for(int i=sz[0];i<=sz[1];i++) if(!Vector<T>::WriteFrame(f,sz+2)) return FALSE;               
	  return TRUE;
     }
     int         Read(FILE* f) {
	  int mysize[4];
	  if(fread(mysize,sizeof(int),4,f)!=4) return FALSE;
	  //assert(mysize[0]==size[0] && mysize[1]==size[1]);
	  //assert(mysize[2]==size[2] && mysize[3]==size[3]);
	  for(int i=size[0];i<=size[1];i++) if(!data[i].Read(f)) return FALSE;
	  return TRUE;
     }
     int         WritePartial(FILE* f) {                     // write part of a big plane with size bigsize (in file)
	  int bigsize[4];
	  if(fread(bigsize,sizeof(int),4,f)!=4) return FALSE;
	  //assert(INBOUNDS(size,bigsize,0) && INBOUNDS(size,bigsize,2));
	  long int filepos=ftell(f);     
	  if(fseek(f,data[size[0]].Size(&(bigsize[2]))*(size[0]-bigsize[0]),SEEK_CUR)!=0) return FALSE;
	  for(int i=size[0];i<=size[1];i++) if(!data[i].WritePartial(f)) return FALSE;
	  if(fseek(f,filepos+data[size[0]].Size(&(bigsize[2]))*(bigsize[1]-bigsize[0]+1),SEEK_SET)!=0) return FALSE;
	  return TRUE;
     }
     int         WritePartial(FILE* f,const int* sz) {                     // write part of a big plane with size bigsize (in file)
	  int bigsize[4];
	  if(fread(bigsize,sizeof(int),4,f)!=4) return FALSE;
	  //assert(INBOUNDS(sz,size,0) && INBOUNDS(sz,size,2));
	  //assert(INBOUNDS(sz,bigsize,0) && INBOUNDS(sz,bigsize,2));
	  long int filepos=ftell(f);     
	  if(fseek(f,data[sz[0]].Size(&(bigsize[2]))*(sz[0]-bigsize[0]),SEEK_CUR)!=0) return FALSE;
	  for(int i=sz[0];i<=sz[1];i++) if(!data[i].WritePartial(f,sz+2)) return FALSE;
	  if(fseek(f,filepos+data[sz[0]].Size(&(bigsize[2]))*(bigsize[1]-bigsize[0]+1),SEEK_SET)!=0) return FALSE;
	  return TRUE;
     }
     int         ReadPartial(FILE* f) {
	  int bigsize[4];
	  if(fread(bigsize,sizeof(int),4,f)!=4) return FALSE;
	  //assert(INBOUNDS(size,bigsize,0) && INBOUNDS(size,bigsize,2));
	  long int filepos=ftell(f);     
	  if(fseek(f,data[size[0]].Size(&(bigsize[2]))*(size[0]-bigsize[0]),SEEK_CUR)!=0) return FALSE;
	  for(int i=size[0];i<=size[1];i++) if(!data[i].ReadPartial(f)) return FALSE;
	  if(fseek(f,filepos+data[size[0]].Size(&(bigsize[2]))*(bigsize[1]-bigsize[0]+1),SEEK_SET)!=0) return FALSE;
	  return TRUE;
     }           
     int         ReadPartial(FILE* f,int* sz) {                     // write part of a big plane with size bigsize (in file)
	  int bigsize[4];
	  if(fread(bigsize,sizeof(int),4,f)!=4) return FALSE;
	  //assert(INBOUNDS(sz,size,0) && INBOUNDS(sz,size,2));
	  //assert(INBOUNDS(sz,bigsize,0) && INBOUNDS(sz,bigsize,2));
	  long int filepos=ftell(f);     
	  if(fseek(f,(long)data[sz[0]].Size(&(bigsize[2]))*(sz[0]-bigsize[0]),SEEK_CUR)!=0) return FALSE;
	  for(int i=sz[0];i<=sz[1];i++) if(!data[i].ReadPartial(f,sz+2)) return FALSE;
	  if(fseek(f,(long)filepos+data[sz[0]].Size(&(bigsize[2]))*(bigsize[1]-bigsize[0]+1),SEEK_SET)!=0) return FALSE;
	  return TRUE;
     }
     T           Min() {
	  T min=data[size[1]].Min(),minc;
	  for(int i=size[0];i<size[1];i++) if(min>(minc=data[i].Min())) min=minc;
	  return min;
     }
     T           Max() {
	  T max=data[size[1]].Max(),maxc;
	  for(int i=size[0];i<size[1];i++) if(max<(maxc=data[i].Max())) max=maxc;
	  return max;
     }
     int         Size(int* bigsize) {return(4*sizeof(int)+data[size[0]].Size(&(bigsize[2]))*(bigsize[1]-bigsize[0]+1));}
     ~Plane() {
	  data+=size[0];
	  if(data) delete[] data;
     }
 protected:
     //! array to store the index ranges
     int         size[4];
     //! pointer to the actual data
     Vector<T>*  data;
};

//! Matrix datatype
/*!
  This class provides a datatype which is a 3D-field. It is built from the Plane datatype.
  It provides an index-operator in the same way as Plane and vector do. This results in a C-style
  indexing scheme into the 3D-field.
  The Matrix datatype is used to store the field data like pressure or velocity components.
 */ 
template < class T >
class Matrix {
     
 public:
     //! default constructor
     Matrix() {
	  data=0;
	  size[0]=0;
     }
     int         Init(int* sz) {
	  data+=size[0];
	  if(data) delete[] data;
	  //assert(sz[0]<=sz[1] && sz[2]<=sz[3] && sz[4]<=sz[5]);
	  int i;
	  for(i=0;i<6;i++) size[i]=sz[i];
	  data=new Plane<T>[size[1]-size[0]+1];
	  if(data==0) return FALSE;
	  data-=size[0];
	  for(i=size[0];i<=size[1];i++) data[i].Init(&(size[2]));
	  return TRUE;
     }
     Matrix(int* sz) {
	  data=0;
	  size[0]=0;
	  Init(sz);
     }
     //! acces operator
     /*!
       This access operator provides C-style indexing into the 3D-field.
       \param int i index into the Matrix
       \return reference to Plane
      */
     inline Plane<T>& operator[] (int i) const {
#ifndef NDEBUG
	  if(!((data + size[0]) && (i >= size[0]) && (i <= size[1]))) {
	       printf("Matrix access at %d,?,? (%d/%d)\n",i,size[0],size[1]);
	       //assert((data + size[0]) && (i >= size[0]) && (i <= size[1]));
	  }
#endif
	  return data[i];
     }
     
     int         Write(FILE * f) {
	  if(fwrite(size,sizeof(int),6,f)!=6) return FALSE;
	  for(int i=size[0];i<=size[1];i++) if(!data[i].Write(f)) return FALSE;
	  return TRUE;
     }
     
     int         WriteParallel(FILE* f,
			       int(*writevect)(FILE* f,int,int,int*,const ParParams*,const Matrix<T>&,const int*),
			       const ParParams* par,const int* dims) {
	  
	  int bigsize[6] ;
	  fseek(f,0,SEEK_CUR);
	  if(fread(bigsize,sizeof(int),6,f)!=6) return FALSE;
	  fseek(f,0,SEEK_CUR);
	  
	  for(int i=bigsize[0];i<=bigsize[1];i++) {
	       if ((i%10)==0){ printf("\rWrite %d / %d",i,bigsize[1]) ;fflush(stdout);}
	       if (!data[i].WriteParallel(f,writevect,par,i,*this,dims , &bigsize[2])) return FALSE;
	  }
	  printf("\n");
	  return TRUE;
     }
     
     static int  WriteFrame(FILE* f,int* sz) {
	  if(fwrite(sz,sizeof(int),6,f)!=6) return FALSE;
	  for(int i=sz[0];i<=sz[1];i++) if(!Plane<T>::WriteFrame(f,sz+2)) return FALSE;
	  return TRUE;
     }
     int         Skip(FILE* f) {
	  int bigsize[6];
	  if(fread(bigsize,sizeof(int),6,f)!=6) return FALSE;
	  //for(int i=0; i<6;i++) assert(INBOUNDS(size,bigsize,0) && INBOUNDS(size,bigsize,2) && INBOUNDS(size,bigsize,4));
	  long int filepos=ftell(f);     
	  if(fseek(f,filepos+data[size[0]].Size(&(bigsize[2]))*(bigsize[1]-bigsize[0]+1),SEEK_SET)!=0) return FALSE;
	  return TRUE;
     }
     int         Read(FILE* f) {
	  int mysize[6];
	  if(fread(mysize,sizeof(int),6,f)!=6) return FALSE;
	  //assert(mysize[0]==size[0] && mysize[1]==size[1]);
	  //assert(mysize[2]==size[2] && mysize[3]==size[3]);
	  //assert(mysize[4]==size[4] && mysize[5]==size[5]);
	  for(int i=size[0];i<=size[1];i++) if(!data[i].Read(f)) return FALSE;
	  return TRUE;
     }
     int         WritePartial(FILE* f) {                     // write part of a big plane with size bigsize (in file)
	  int bigsize[6];
	  if(fread(bigsize,sizeof(int),6,f)!=6) return FALSE;
	  //assert(INBOUNDS(size,bigsize,0) && INBOUNDS(size,bigsize,2) && INBOUNDS(size,bigsize,4));
	  long int filepos=ftell(f);     
	  if(fseek(f,data[size[0]].Size(&(bigsize[2]))*(size[0]-bigsize[0]),SEEK_CUR)!=0) return FALSE;
	  for(int i=size[0];i<=size[1];i++) if(!data[i].WritePartial(f)) return FALSE;
	  if(fseek(f,filepos+data[size[0]].Size(&(bigsize[2]))*(bigsize[1]-bigsize[0]+1),SEEK_SET)!=0) return FALSE;
	  return TRUE;
     }
     int         WritePartial(FILE* f,const int* sz) {                     // write part of a big plane with size bigsize (in file)
	  int bigsize[6];int i;
	  if(fread(bigsize,sizeof(int),6,f)!=6) return FALSE;
	  //assert(INBOUNDS(sz,size,0) && INBOUNDS(sz,size,2) && INBOUNDS(sz,size,4));
	  //assert(INBOUNDS(sz,bigsize,0) && INBOUNDS(sz,bigsize,2) && INBOUNDS(sz,bigsize,4));
	  long int filepos=ftell(f);     
	  if(fseek(f,data[sz[0]].Size(&(bigsize[2]))*(sz[0]-bigsize[0]),SEEK_CUR)!=0) return FALSE;
	  for(i=sz[0];i<=sz[1];i++) if(!data[i].WritePartial(f,sz+2)) return FALSE;
	  if(fseek(f,filepos+data[sz[0]].Size(&(bigsize[2]))*(bigsize[1]-bigsize[0]+1),SEEK_SET)!=0) return FALSE;
	  return TRUE;
     }
     int         ReadPartial(FILE* f) {
	  int bigsize[6];
	  if(fread(bigsize,sizeof(int),6,f)!=6) return FALSE;
	  //assert(INBOUNDS(size,bigsize,0) && INBOUNDS(size,bigsize,2) && INBOUNDS(size,bigsize,4));
	  long int filepos=ftell(f);     
	  if(fseek(f,data[size[0]].Size(&(bigsize[2]))*(size[0]-bigsize[0]),SEEK_CUR)!=0) return FALSE;
	  for(int i=size[0];i<=size[1];i++) if(!data[i].ReadPartial(f)) return FALSE;
	  if(fseek(f,filepos+data[size[0]].Size(&(bigsize[2]))*(bigsize[1]-bigsize[0]+1),SEEK_SET)!=0) return FALSE;
	  return TRUE;
     }           
     int         ReadPartial(FILE* f,int* sz) {                     // write part of a big plane with size bigsize (in file)
	  int bigsize[6];int i;
	  if(fread(bigsize,sizeof(int),6,f)!=6) return FALSE;
	  //assert(INBOUNDS(sz,size,0) && INBOUNDS(sz,size,2) && INBOUNDS(sz,size,4));
	  //assert(INBOUNDS(sz,bigsize,0) && INBOUNDS(sz,bigsize,2) && INBOUNDS(sz,bigsize,4));
	  long int filepos=ftell(f);     
	  if(fseek(f,data[sz[0]].Size(&(bigsize[2]))*(sz[0]-bigsize[0]),SEEK_CUR)!=0) return FALSE;
	  for(i=sz[0];i<=sz[1];i++) if(!data[i].ReadPartial(f,sz+2)) return FALSE;
	  if(fseek(f,filepos+data[sz[0]].Size(&(bigsize[2]))*(bigsize[1]-bigsize[0]+1),SEEK_SET)!=0) return FALSE;
	  return TRUE;
     }
// ------- June99: New Routines to speed up writing of binfile: WriteHead,WriteSlice
  
     int         WriteHead(FILE *f, int *sz) {
	  // Write bounds of the matrix
	  if(fwrite(sz,sizeof(int),6,f)!=6) return FALSE;
	  return TRUE;
     }

     int         WriteSlice(FILE* f,int* sz) {
		
	  if (sz[0]!=sz[1])
	  {
	       printf("\nERROR in Proced. WriteSlice: Matrix is not a Slice\n");  
	       return FALSE;
	  }
	  // Write one Slice
	  if(!data[sz[0]].Write(f)) return FALSE;
	  return TRUE;
     }
// -------

     T           Min() {
	  T min=data[size[1]].Min(),minc;
	  for(int i=size[0];i<size[1];i++) if(min>(minc=data[i].Min())) min=minc;
	  return min;
     }
     T           Max() {
	  T max=data[size[1]].Max(),maxc;
	  for(int i=size[0];i<size[1];i++) if(max<(maxc=data[i].Max())) max=maxc;
	  return max;
     }
     void        Dump(char* str) {
	  int i,j,k;
	  for(i=size[0];i<=size[1];i++) for(j=size[2];j<=size[3];j++) for(k=size[4];k<=size[5];k++) 
	       printf("%s[%3d][%3d][%3d]=%le\n",str,i,j,k,data[i][j][k]);
     }
     ~Matrix() {
	  data+=size[0];
	  if(data) delete[] data;
     }
 protected:
     //! array to store the index range
     int         size[6];
     //! pointer to actual data(Planes)
     Plane<T>*   data;
};
#endif
