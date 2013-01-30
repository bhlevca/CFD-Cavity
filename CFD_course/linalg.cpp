#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <ostream>

//#include <boost/python.hpp>   // if you include this, then you don't
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/module.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/python/list.hpp>
#include <boost/python/dict.hpp>
#include <boost/python/numeric.hpp>
#include <boost/python/tuple.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/python.hpp>
#include <numpy/noprefix.h>
#include "gauss_seidel.h"
#include "thomas.h"



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
	
	//allocate the matrix and the vectors
	matrix<double> A;
	vector<double> b ;   // free term for GS or subdiagonal for TDMA
	vector<double> x0;
	double omega;
	int    kmax;
	double precision;
	std::vector<T> nodes;

	//For TDMA
	vector<double> a ;      //main diagonal
	vector<double> c ;      //supradiagonal
	vector<double> q ;      //free term for TDMA
	vector<double> x ;      //solution


	std::vector<T> getData (void);
    int getSize (void);
    double sum(list o);
    double sum(double *x, int n);

    //Gauss Seidel Methods
    void setB(numeric::array y, int n);
    void setX0(numeric::array y, int n);
    void setA(numeric::array y, int n, int m);
    void setParams(double omega, int kmax, double ACC );
    boost::python::dict solve(numeric::array x, int n);

    //Thomas Methods
    void setTDMA(numeric::array bb,
      	         numeric::array aa,
    		     numeric::array cc,
    			 numeric::array qq,
    			 int nn);

    boost::python::dict solveTDMA (numeric::array xx, int nn);

    void   clean();

};

// Class definitions
template <typename T> std::vector<T> CppBlas<T>::getData () {
    return this->nodes;
}
// Class definitions
template <typename T> int CppBlas<T>::getSize () {
    return this->nodes.size();
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

template<typename T> void CppBlas<T>::setX0(numeric::array y, int n) {
      x0.resize(n);
      for(int i = 0; i < n; i++)         {
    	  x0(i) = extract<double>(y[i]);
      }
}

template<typename T> void CppBlas<T>::setB(numeric::array y, int n) {
	  b.resize(n);

	  for(int i = 0; i < n; i++)         {
		  b(i) = extract<double>(y[i]);
	  }
}


template<typename T> void CppBlas<T>::setA(numeric::array y, int n, int m) {

	 A.resize(n,m);
	 for(int i = 0; i < n; i++)    {
		  for(int j = 0; j < m; j++)    {
			  A(i,j) = extract<double>(y[i][j]);
			  //std::cout<< "A[i][j]:" << A(i,j)<< std::endl;
		  }
	  }
}

template<typename T> void CppBlas<T>::setParams(double omega, int kmax, double ACC ){
	this->omega= omega;
	this->precision = ACC;
	this->kmax = kmax;
}


template <typename T> boost::python::dict CppBlas<T>::solve (numeric::array x, int n) {
	//std::cout<< "omega:" << omega<< " precision:" << precision << " kmax:" << kmax << std::endl;
	//std::cout<< "A" << A<< " b:" << b << " x0" << x0 << std::endl;

	struct gsResult* r=gauss_seidel(A,b,x0,omega,kmax,precision);
	//std::cout<< "r.x" << r->x << " r->kiter:" << r->kiter << std::endl;
	for (unsigned int i=0; i< r->x.size() ; i++) {
			x[i]= r->x(i);
	}
	boost::python::dict retval;
	retval["solution"] = x;
	retval["kiter"]=r->kiter;
	return retval;
}

template<typename T> void CppBlas<T>::setTDMA(numeric::array bb,
											  numeric::array aa,
											  numeric::array cc,
											  numeric::array qq,
											  int nn){


	 a.resize(nn);
	 b.resize(nn-1);
	 c.resize(nn-1);
	 q.resize(nn);

	 for(int i = 0; i < nn; i++)    {
		  a(i) = extract<double>(aa[i]);
		  q(i) = extract<double>(qq[i]);
	 }

	 for(int i = 0; i < nn-1; i++)    {

		 b(i) = extract<double>(bb[i]);
		 c(i) = extract<double>(cc[i]);
		  //std::cout<< "A[i][j]:" << A(i,j)<< std::endl;
     }

}

template <typename T> boost::python::dict CppBlas<T>::solveTDMA (numeric::array xx, int nn) {
	//std::cout<< "A" << A<< " b:" << b << " x0" << x0 << std::endl;
	x.resize(nn);
	ThomasAlgorithm(4, b, a, c, x, q);
	for (unsigned int i=0; i< x.size() ; i++) {
			xx[i]= x(i);
	}
	boost::python::dict retval;
	retval["solution"] = xx;
	return retval;
}


template<typename T> void CppBlas<T>::clean() {
	//delete alocations

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
        .def("clean", &CppBlas<Data>::clean)
        .def("solve",(boost::python::dict ( ::CppBlas<Data>::* )( numeric::array,int) )(&CppBlas<Data>::solve ), ( arg("x"), arg("n")))
        .def("sum",(double ( ::CppBlas<Data>::* )( list ) )(&CppBlas<Data>::sum ), ( arg("o") ) )
        .def("setB", (void ( ::CppBlas<Data>::* )( numeric::array,int ) )( &CppBlas<Data>::setB ), ( arg("y"), arg("n") ) )
        .def("setX0", (void ( ::CppBlas<Data>::* )( numeric::array,int ) )( &CppBlas<Data>::setX0 ), ( arg("y"), arg("n") ) )
        .def("setA", (void ( ::CppBlas<Data>::* )( numeric::array,int,int ) )( &CppBlas<Data>::setA ), ( arg("y"), arg("n"), arg("m") ) )
    	.def("setParams", (void ( ::CppBlas<Data>::* )( double, int, double ) )( &CppBlas<Data>::setParams ), ( arg("omega"), arg("kmax"), arg("ACC") ) )
    	.def("setTDMA", (void ( ::CppBlas<Data>::* )( numeric::array, numeric::array, numeric::array, numeric::array, int) )( &CppBlas<Data>::setTDMA ), ( arg("bb"), arg("aa"), arg("cc"),arg("qq"), arg("nn") ) )
    	.def("solveTDMA",(boost::python::dict ( ::CppBlas<Data>::* )( numeric::array, int ) )(&CppBlas<Data>::solveTDMA ), ( arg("xx"), arg("nn") ) )
    	;
    ;


}

