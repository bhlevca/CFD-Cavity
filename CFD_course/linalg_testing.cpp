/*
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <ostream>

//#include <boost/python.hpp>   // if you include this, then you don't
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/module.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/python/list.hpp>
#include <boost/python/numeric.hpp>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include "gauss_seidel.h"

#include <boost/python/stl_iterator.hpp>

using namespace boost::python;




// Class Declarations
class Data {
    public:
    double value;

    // Commented-out methods are not strictly necessary, but nice to
    // implement so that your object will act Pythonic

    Data():value(0.0) {}
    Data(double value):value(value) {}
    double repr() const { return value; }
    void reset() { value = 0.0; }
    bool operator==(Data const& n) const { return value == n.value; }
    bool operator!=(Data const& n) const { return value != n.value; }
};

template <typename T> class CppBlas {
    public:
	
	double **  initialA;

	//define the initial values
	double* initialX;

	//define the right hand side
	double* initialB;


	//allocate the matrix and the vectors
	matrix<double>* A;
	vector<double>* b ;
	vector<double>* x0;
    
	std::vector<T> nodes;

    std::vector<T> getData (void);
    int getSize (void);
    void setB (object o);
    double sum(list o);
    double sum(numeric::array y, int n);
    double sum(double *x, int n);
    double sum(numeric::array y, int n, int m);
};

// Class definitions
template <typename T> std::vector<T> CppBlas<T>::getData () {
    return this->nodes;
}
// Class definitions
template <typename T> int CppBlas<T>::getSize () {
    return this->nodes.size();
}


template<typename T> void CppBlas<T>::setB(object b) {

    try
      {
        object iter_obj = object( handle<>( PyObject_GetIter( b.ptr() ) ));
        while( 1 ) {
        	object obj = extract<object>( iter_obj.attr( "next" )() ); // Should   always work
        	double val = extract<double>( obj ); // Should launch an exception if you wannot extract an double ...
        	nodes.push_back(val);
        }
      }catch( error_already_set )
      {
        PyErr_Clear(); // If there is an exception (no iterator, extract failed or end of the list reached), clear it and exit the function
        return;
      }
}

template<typename T> void CppBlas<T>::setX0(object x0) {

    try
      {
        object iter_obj = object( handle<>( PyObject_GetIter( x0.ptr() ) ));
        while( 1 ) {
        	object obj = extract<object>( iter_obj.attr( "next" )() ); // Should   always work
        	double val = extract<double>( obj ); // Should launch an exception if you wannot extract an double ...
        	nodes.push_back(val);
        }
      }catch( error_already_set )
      {
        PyErr_Clear(); // If there is an exception (no iterator, extract failed or end of the list reached), clear it and exit the function
        return;
      }
}


template<typename T> void CppBlas<T>::setA(object A) {

    try
      {
        object row_obj = object( handle<>( PyObject_GetIter( A.ptr() ) ));
        while ( 1 ) {

        	object iter_obj = object( handle<>( PyObject_GetIter( row.ptr() ) ));
        	object obj = extract<object>( row_obj.attr( "next" )() ); // Should   always work

			while( 1 ) {
				object col_obj = extract<object>( obj.attr( "next" )() ); // Should   always work
				double val = extract<double>( obj ); // Should launch an exception if you wannot extract an double ...
				nodes.push_back(val);
			}
			catch( error_already_set ){
			     PyErr_Clear(); // If there is an exception (no iterator, extract failed or end of the list reached),
			     // GO TO THE NEXT LINE
			}
        }catch( error_already_set ) {
        	PyErr_Clear(); // If there is an exception (no iterator, extract failed or end of the list reached), clear it and exit the function
        	return;
      }
}


template<typename T> double CppBlas<T>::sum(double *x, int n)
{
  double thesum = 0.;
  for(int i = 0; i < n; i++)
    {
      thesum += x[i];
    }
  return(thesum);
}

template<typename T> double CppBlas<T>::sum(list o) {
      std::size_t n = len(o);
      double *tmp = new double[n];
      for(unsigned int i = 0; i < n; i++)
        {
          tmp[i] = extract<double>(o[i]);
        }
      double thesum = CppBlas<T>::sum(tmp, n);
      delete tmp;
      return(thesum);
}

template<typename T> double CppBlas<T>::setX0(numeric::array y, int n) {
      double *tmp = new double[n];
      for(int i = 0; i < n; i++)
        {
          tmp[i] = extract<double>(y[i]);
        }
      double thesum = CppBlas<T>::sum(tmp, n);
      delete tmp;
      return(thesum);
}

template<typename T> double CppBlas<T>::setB(numeric::array y, int n) {
      double *tmp = new double[n];
      for(int i = 0; i < n; i++)
        {
          tmp[i] = extract<double>(y[i]);
        }
      double thesum = CppBlas<T>::sum(tmp, n);
      delete tmp;
      return(thesum);
}


template<typename T> double CppBlas<T>::setA(numeric::array y, int n, int m) {
      double *tmp = new double[m];
      double * sum = new double[n];
      for(int i = 0; i < n; i++) {
    	  for (int j = 0; j < m; j++) {
    		  tmp[j] = extract<double>(y[i][j]);
    	  }
    	  sum[i] = CppBlas<T>::sum(tmp, m);
      }
      double thesum = CppBlas<T>::sum(sum, n);
      delete tmp;
      delete sum;
      return(thesum);
}



// Boost.python definitions to expose classes to Python
BOOST_PYTHON_MODULE(linalg) {
	numeric::array::set_module_and_type("numpy", "ndarray");
    class_<Data> ("Data")
        .def_readwrite("value", &Data::value)
    ;
    // Line below is necessary to expose the vector "nodes" to Python
    // and have it function as expected in Python
    class_< std::vector<Data> > ("DataArray")
        .def(vector_indexing_suite< std::vector<Data> >())
    ;

    class_<CppBlas<Data> > ("CppBlas")
        .def_readwrite ("nodes", &CppBlas<Data>::nodes)
        .def("getData", &CppBlas<Data>::getData)
        .def("getSize", &CppBlas<Data>::getSize)
        .def("setList": &CppBlas<Data>::setList)
        .def("sum",(double ( ::CppBlas<Data>::* )( list ) )(&CppBlas<Data>::sum ), ( arg("o") ) )
        .def("setB", (double ( ::CppBlas<Data>::* )( numeric::array,int ) )( &CppBlas<Data>::setB ), ( arg("y"), arg("n") ) )
        .def("setX0", (double ( ::CppBlas<Data>::* )( numeric::array,int ) )( &CppBlas<Data>::sumX0 ), ( arg("y"), arg("n") ) )
        .def("setA", (double ( ::CppBlas<Data>::* )( numeric::array,int,int ) )( &CppBlas<Data>::setA ), ( arg("y"), arg("n"), arg("m") ) );
    ;


}
*/
