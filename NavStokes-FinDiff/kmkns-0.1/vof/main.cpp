#include <Python.h>
#include <numpy/arrayobject.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

//#define DEBUG
//#undef NDEBUG
#include <cassert>

#include "vof.h"

using namespace vof;

const double EPSILON = 1e-12;

template<class T>
inline T& index1D(PyArrayObject* v, int i)
{
	assert(v && "Trying to index NULL.");
	assert(v->nd == 1 && "Incorrect array dimensions.");
	assert(i >= 0 && i < v->dimensions[0] && "Index out of bounds.");
	return *((T*)(v->data + i * v->strides[0]));
}

template<class T>
inline T& index2D(PyArrayObject* v, int i, int j)
{
	assert(v && "Trying to index NULL.");
	assert(v->nd == 2 && "Incorrect array dimensions.");
	assert(i >= 0 && i < v->dimensions[0] && "Index out of bounds.");
	assert(j >= 0 && j < v->dimensions[1] && "Index out of bounds.");
	return *((T*)(v->data + i * v->strides[0] + j * v->strides[1]));
}

// Pick a 3x3 block from 'phi' for easier indexing.
static void local_grid(PyArrayObject* phi, double* f, int i, int j)
{
	for(int l = 0; l < 3; ++l)
		for(int k = 0; k < 3; ++k)
			f[3*l+k] = index2D<double>(phi, i+k-1, j+l-1);
}

static PyArrayObject* create_similar_array(PyArrayObject* v)
{
	int* dims = new int[v->nd];
	for(int i = 0; i < v->nd; ++i)
		dims[i] = (int)v->dimensions[i];
	PyObject* w = PyArray_FromDims(v->nd, dims, v->descr->type_num);
	delete dims;
	return (PyArrayObject*)w;
}

static void gradient_main(PyArrayObject* phi, PyArrayObject* mask, PyArrayObject* delta, 
						  PyArrayObject* gx, PyArrayObject* gy)
{
	double dx = index1D<double>(delta, 0);
	double dy = index1D<double>(delta, 1);
	double r2p1 = 1.0 + sqr(dx / dy);
	double invr2p1 = 1.0 + sqr(dx / dy);
	double factor_x = 1.0/(dx * (4.0 + 2.0 * invr2p1));
	double factor_y = 1.0/(dy * (4.0 + 2.0 * r2p1));
	double f[9];
	
	// Stencil w/weights for x-component:
	//        -1    0    1
	//          +---+---+
	//          |\  |  /|
	//          | \ | / |
	//          |  \|/  |
	// -invr2p1 +---+---+ invr2p1
	//          |  /|\  |
	//          | / | \ |
	//          |/  |  \|
	//          +---+---+
	//        -1    0    1

	for(int j = 1; j < phi->dimensions[1] - 1; ++j)
	{
		for(int i = 1; i < phi->dimensions[0] - 1; ++i)
		{
			if(index2D<long>(mask, i, j) & 1)
			{
				local_grid(phi, f, i, j);
				index2D<double>(gx, i, j) = factor_x * (f[8] + f[2] - f[6] - f[0] + invr2p1 * (f[5] - f[3]));
				index2D<double>(gy, i, j) = factor_y * (f[8] + f[6] - f[2] - f[0] + r2p1 * (f[7] - f[1]));
			}
		}
	}
}

static char gradient_doc[] = "gx, gy = gradient(phi, mask, delta)\n"
	"\n"
	"Calculate the gradient of 'phi' using a 3x3 stencil.\n"
	"'phi' and 'mask' must have ghost cells around the domain.\n"
	"'gx' and 'gy' will also have ghost cells, but the gradients in the ghost\n"
	"cells are undefined.\n"
	"Ghost cells and obstacle cells must have suitable values.\n"
	"\n"
	"  Input:\n"
	"  - phi = array of volume-of-fluid values\n"
	"  - mask = array of flags, 1 for fluid cell, 0 for obstacle cell\n"
	"  - delta = grid spacing dx and dy\n"
	"  Output:\n"
	"  - gx = array of gradient x-components\n"
	"  - gy = array of gradient y-components\n";

static PyObject* gradient(PyObject* self, PyObject* args)
{
	PyObject* result = NULL;
	PyObject* phi_;
	PyObject* mask_;
	PyObject* delta_;
	PyArrayObject* phi = NULL;
	PyArrayObject* mask = NULL;
	PyArrayObject* delta = NULL;
	PyArrayObject* gx = NULL;
	PyArrayObject* gy = NULL;

	if(!PyArg_ParseTuple(args, "OOO:gradient", &phi_, &mask_, &delta_))
		return NULL;

	bool failed = false;
	if(!failed && !(phi = (PyArrayObject*)PyArray_FromObject(phi_, PyArray_DOUBLE, 2, 2)))
		failed = true;
	if(!failed && !(mask = (PyArrayObject*)PyArray_FromObject(mask_, PyArray_LONG, 2, 2)))
		failed = true;
	if(!failed && (phi->dimensions[0] != mask->dimensions[0] || phi->dimensions[1] != mask->dimensions[1])) 
	{
		PyErr_Format(PyExc_TypeError, "'phi' and 'mask' must have the same size.");
		failed = true;
	}

	if(!failed && !(delta = (PyArrayObject*)PyArray_FromObject(delta_, PyArray_DOUBLE, 1, 1)))
		failed = true;
	if(!failed && delta->dimensions[0] != 2)
	{
		PyErr_Format(PyExc_TypeError, "'delta' must consist of two elements.");
		failed = true;
	}

	if(!failed && !(gx = create_similar_array(phi)))
		failed = true;
	if(!failed && !(gy = create_similar_array(phi)))
		failed = true;

	if(!failed)
	{
		gradient_main(phi, mask, delta, gx, gy);
		result = Py_BuildValue("OO", gx, gy);
	}
	Py_XDECREF(phi);
	Py_XDECREF(mask);
	Py_XDECREF(delta);
	Py_XDECREF(gx);
	Py_XDECREF(gy);
	return result;
}

static void reconstruct_main(PyArrayObject* phi, PyArrayObject* mask, PyArrayObject* delta, 
							 PyArrayObject* nx, PyArrayObject* ny, PyArrayObject* c)
{
	double dx = index1D<double>(delta, 0);
	double dy = index1D<double>(delta, 1);

	double f[9];
	for(int j = 1; j < phi->dimensions[1] - 1; ++j)
	{
		for(int i = 1; i < phi->dimensions[0] - 1; ++i)
		{
			if(index2D<long>(mask, i, j) & 1)
			{
				local_grid(phi, f, i, j);
				reconstruct(index2D<double>(nx, i, j), index2D<double>(ny, i, j), index2D<double>(c, i, j), f, dx, dy);
			}
			else
			{
				index2D<double>(nx, i, j) = 0.0;
				index2D<double>(ny, i, j) = 0.0;
				index2D<double>(c, i, j) = 0.0;
			}
		}
	}
}

static char reconstruct_doc[] = "nx, ny, c = reconstruct(phi, mask, delta)\n"
	"\n"
	"Reconstruct the interface in the volume-of-fluid method.\n"
	"'phi' defines how much of each cell is filled with fluid.\n"
	"Zero indicates an empty cell, one a full cell.\n"
	"If the value is in the range (0, 1), an interface will pass\n"
	"through the cell. (nx, ny) is the normal of the interface in\n"
	"such cells. The returned value 'c' is the signed distance\n"
	"from the centre of a cell to its interface. 'c' is undefined for\n"
	"cells with 'phi'=0 or 'phi'=1.\n"
	"\n"
	"  Input:\n"
	"  - phi   = array of volume-of-fluid values\n"
	"  - delta = grid spacing dx and dy\n"
	"  Output:\n"
	"  - nx    = array of normal x-components\n"
	"  - ny    = array of normal y-components\n"
	"  - c     = array of gradient x-components\n";

static PyObject* reconstruct(PyObject* self, PyObject* args)
{
	PyObject* result = NULL;
	PyObject* phi_;
	PyObject* mask_;
	PyObject* delta_;
	PyArrayObject* phi = NULL;
	PyArrayObject* mask = NULL;
	PyArrayObject* delta = NULL;
	PyArrayObject* nx = NULL;
	PyArrayObject* ny = NULL;
	PyArrayObject* c = NULL;

	if(!PyArg_ParseTuple(args, "OOO:reconstruct", &phi_, &mask_, &delta_))
		return NULL;

	bool failed = false;
	if(!failed && !(phi = (PyArrayObject*)PyArray_FromObject(phi_, PyArray_DOUBLE, 2, 2)))
		failed = true;
	if(!failed && !(mask = (PyArrayObject*)PyArray_FromObject(mask_, PyArray_DOUBLE, 2, 2)))
		failed = true;
	if(!failed && (phi->dimensions[0] != mask->dimensions[0] || phi->dimensions[1] != mask->dimensions[1]))
	{
		PyErr_Format(PyExc_TypeError, "'mask' and 'phi' must have the same size.");
		failed = true;
	}

	if(!failed && !(delta = (PyArrayObject*)PyArray_FromObject(delta_, PyArray_DOUBLE, 1, 1)))
		failed = true;
	if(!failed && delta->dimensions[0] != 2) 
	{
		PyErr_Format(PyExc_TypeError, "'delta' must have two elements.");
		failed = true;
	}

	if(!failed && !(nx = create_similar_array(phi)))
		failed = true;
	if(!failed && !(ny = create_similar_array(phi)))
		failed = true;
	if(!failed && !(c = create_similar_array(phi)))
		failed = true;

	if(!failed)
	{
		reconstruct_main(phi, mask, delta, nx, ny, c);
		result = Py_BuildValue("OOO", nx, ny, c);
	}
	Py_XDECREF(phi);
	Py_XDECREF(mask);
	Py_XDECREF(delta);
	Py_XDECREF(nx);
	Py_XDECREF(ny);
	Py_XDECREF(c);
	return result;
}

// semicircle function.
static double f(double x, double cx, double cy, double r)
{
	return sqrt(vof::max(r * r - sqr(x - cx), 0.0)) + cy;
}

// integral of f().
static double F(double x, double cx, double cy, double r)
{
	double ymy0 = sqrt(vof::max(r * r - sqr(x - cx), 0.0));
	return 0.5 * (r * r * acos(clamp((cx - x) / r, -1.0, 1.0)) - r * ymy0 + (x - cx + r) * (2 * cy + ymy0));
}

// partial derivative of F() with respect to cx.
static double dFdx(double x, double cx, double cy, double r)
{
	return -f(x, cx, cy, r);
}

// partial derivative of F() with respect to cy.
static double dFdy(double x, double cx, double cy, double r)
{
	return x - cx + r;
}

// partial derivative of F() with respect to r.
static double dFdr(double x, double cx, double cy, double r)
{
	return cy + r * acos(clamp((cx - x) / r, -1.0, 1.0));
}

// determinant of 3x3 matrix, stored row by row
static double det3x3(double* m)
{
	double det = 0.0;
	for(int i = 0; i < 3; ++i)
		det += (m[i] * m[3+(i+1)%3] * m[6+(i+2)%3]) - (m[i] * m[3+(i+2)%3] * m[6+(i+1)%3]);
	return det;
}

// Use Newton's method to find the circle that corresponds to the given heights.
static bool refine_curvature(double dx, double* heights, double& cx, double& cy, double& r)
{
	// Solve the following equation system:
	// F(-dx/2)  - F(-3*dx/2) = heights[0]
	// F(dx/2)   - F(-dx/2)   = heights[1]
	// F(3*dx/2) - F(dx/2)    = heights[2]
	//
	// Let:
	//                [ F(-dx/2)  - F(-3*dx/2) - heights[0] ]
	// G(cx, cy, r) = [ F(dx/2)   - F(-dx/2)   - heights[1] ]
	//                [ F(3*dx/2) - F(dx/2)    - heights[2] ]
	// x_0 = [cx, cy, r]
	//
	// One iteration is given by:
	// x_{n+1} = x_n - J(x_n)^{-1}G(x_n)
	// where J is the Jacobian of G.
	
	double G[3], J[9], invJ[9];
	double residual, init_res;

	for(int i = 0; i < 3; ++i)
	{
		G[i] = F((i - 0.5) * dx, cx, cy, r) - F((i - 1.5) * dx, cx, cy, r) - heights[i] * dx;
		J[3*i+0] = dFdx((i - 0.5) * dx, cx, cy, r) - dFdx((i - 1.5) * dx, cx, cy, r);
		J[3*i+1] = dFdy((i - 0.5) * dx, cx, cy, r) - dFdy((i - 1.5) * dx, cx, cy, r);
		J[3*i+2] = dFdr((i - 0.5) * dx, cx, cy, r) - dFdr((i - 1.5) * dx, cx, cy, r);
	}
	residual = init_res = sqr(G[0]) + sqr(G[1]) + sqr(G[2]);

	for(int it = 0; it < 32; ++it)
	{
		// Invert J, only 3x3, so it should be fast enough.
		double det = det3x3(J);

		if(det == 0.0)
		{
			//fprintf(stderr, "FAILED at it=%i, cx=%f, cy=%f, r=%f\n", it, cx, cy, r);
			//fprintf(stderr, "  G=[%f, %f, %f]\n", G[0], G[1], G[2]);
			return false; // failed. cannot calculate next step.
		}

		invJ[0] = (J[4] * J[8] - J[5] * J[7]) / det;
		invJ[3] = (J[5] * J[6] - J[3] * J[8]) / det;
		invJ[6] = (J[3] * J[7] - J[4] * J[6]) / det;
		invJ[1] = (J[2] * J[7] - J[1] * J[8]) / det;
		invJ[4] = (J[0] * J[8] - J[2] * J[6]) / det;
		invJ[7] = (J[1] * J[6] - J[0] * J[7]) / det;
		invJ[2] = (J[1] * J[5] - J[2] * J[4]) / det;
		invJ[5] = (J[2] * J[3] - J[0] * J[5]) / det;
		invJ[8] = (J[0] * J[4] - J[1] * J[3]) / det;

		//for(int i = 0; i < 3; ++i)
		//	for(int j = 0; j < 3; ++j)
		//		fprintf(stderr, "%f, ", invJ[3*i+0] * J[3*0+j] + invJ[3*i+1] * J[3*1+j] + invJ[3*i+2] * J[3*2+j]);
		//fprintf(stderr, "\n");

		// Multiply J^{-1} with G to find x_{n+1}.
		cx -= 1.0 * (invJ[0] * G[0] + invJ[1] * G[1] + invJ[2] * G[2]);
		cy -= 1.0 * (invJ[3] * G[0] + invJ[4] * G[1] + invJ[5] * G[2]);
		r  -= 1.0 * (invJ[6] * G[0] + invJ[7] * G[1] + invJ[8] * G[2]);

		// Force to valid values.
		// If left side of the circle does not span the left interval (-3*dx/2, -dx/2), 
		// enlarge the circle so that it does.
		if(cx - r > -1.5 * dx)
		{
			double delta = 0.5 * ((cx - r) + 1.5 * dx);
			cx -= delta;
			r += delta;
		}
		// If right side of the circle does not span the right interval (dx/2, 3*dx/2), 
		// enlarge the circle so that it does.
		if(cx + r < 1.5 * dx)
		{
			double delta = 0.5 * (1.5 * dx - (cx + r));
			cx += delta;
			r += delta;
		}

		for(int i = 0; i < 3; ++i)
		{
			G[i] = F((i - 0.5) * dx, cx, cy, r) - F((i - 1.5) * dx, cx, cy, r) - heights[i] * dx;
			J[3*i+0] = dFdx((i - 0.5) * dx, cx, cy, r) - dFdx((i - 1.5) * dx, cx, cy, r);
			J[3*i+1] = dFdy((i - 0.5) * dx, cx, cy, r) - dFdy((i - 1.5) * dx, cx, cy, r);
			J[3*i+2] = dFdr((i - 0.5) * dx, cx, cy, r) - dFdr((i - 1.5) * dx, cx, cy, r);
		}
		residual = sqr(G[0]) + sqr(G[1]) + sqr(G[2]);
	}

	/*
	if(residual <= init_res)
	{
		fprintf(stderr, "it=it_max, cx=%f, cy=%f, r=%f\n", cx, cy, r);
		fprintf(stderr, "  G=[%f, %f, %f]\n", G[0], G[1], G[2]);
		fprintf(stderr, "  heights * dx=[%f, %f, %f]\n", heights[0] * dx, heights[1] * dx, heights[2] * dx);
	}
	else
	{
		fprintf(stderr, "FAILED at it=it_max, cx=%f, cy=%f, r=%f\n", cx, cy, r);
		fprintf(stderr, "  G=[%f, %f, %f]\n", G[0], G[1], G[2]);
	}
	*/

	return residual <= init_res;
}

static void curvature_main(PyArrayObject* phi, PyArrayObject* mask, PyArrayObject* delta, 
						   bool refine, PyArrayObject* nx, PyArrayObject* ny, PyArrayObject* kappa)
{
	const int HALF_STENCIL_HEIGHT = 4; // (H+1)/2 where H is actual stencil height
	// If a radius is less than REFINE_RADIUS * dx (or dy), it is iteratively refined 
	// for better accuracy. If zero, the radius is never refined. 
	const double REFINE_RADIUS = 100.0;

	double dx = index1D<double>(delta, 0);
	double dy = index1D<double>(delta, 1);

	for(int j = 1; j < phi->dimensions[1] - 1; ++j)
	{
		for(int i = 1; i < phi->dimensions[0] - 1; ++i)
		{
			double phi_c = index2D<double>(phi, i, j);
			if((index2D<long>(mask, i, j) & 1) && 
				isInterfaceCell(phi_c, index2D<double>(phi, i, j-1), index2D<double>(phi, i-1, j),
				index2D<double>(phi, i, j+1), index2D<double>(phi, i+1, j)))
			{
				if(vof::abs(index2D<double>(nx, i, j) * dx) > vof::abs(index2D<double>(ny, i, j) * dy))
				{
					// vertical interface
					double r[3];
					int dir = (index2D<double>(nx, i, j) > 0.0 ? 1 : -1);
					for(int l = 0; l < 3; ++l)
					{
						r[l] = index2D<double>(phi, i, j+l-1);
						for(int k = i + dir; k != i + HALF_STENCIL_HEIGHT * dir; k += dir)
						{
							if(k < 0 || k >= phi->dimensions[0] || (index2D<long>(mask, k, j+l-1) & 1) == 0 ||
								(index2D<double>(phi, k, j+l-1) <= EPSILON))
								break;

							r[l] += index2D<double>(phi, k, j+l-1);
						}
						for(int k = i - dir; k != i - HALF_STENCIL_HEIGHT * dir; k -= dir)
						{
							if(k < 0 || k >= phi->dimensions[0] || (index2D<long>(mask, k, j+l-1) & 1) == 0 ||
								(index2D<double>(phi, k, j+l-1) >= 1.0 - EPSILON)
								|| (dir * index2D<double>(nx, k, j+l-1) < 0.0))
								break;

							r[l] += index2D<double>(phi, k, j+l-1) - 1.0;
						}
						r[l] *= dx;
					}
					if(j == 1) r[0] = r[1];
					if(j == phi->dimensions[1] - 2) r[2] = r[1];

					double a = (r[0] - 2.0 * r[1] + r[2]) / (dy * dy);
					double b = (r[2] - r[0]) / (2.0 * dy);

					// dx, heights, cx, cy, r)
					double len = sqrt(1.0 + b * b);
					// Initial approximation. Overwrite if better is found.
					index2D<double>(kappa, i, j) = -a / (len * len * len);
					// If radius is less than REFINE_RADIUS times greater than the 
					// grid spacing, refine the curvature.
					if(refine && (len * len * len < REFINE_RADIUS * dy * vof::abs(a)))
					{
						double s = vof::sign(a);
						if(s > 0.0)
							for(int i = 0; i < 3; ++i)
								r[i] = -r[i];

						double n_y = -s * b / len, n_x = 1.0 / len;
						double radius = (len * len * len) / vof::abs(a);
						double cy = -radius * n_y;
						double cx = r[1] - radius * n_x;

						if(refine_curvature(dy, r, cy, cx, radius))
							index2D<double>(kappa, i, j) = -s / radius;
					}
				}
				else
				{
					// horizontal interface
					double c[3];
					int dir = (index2D<double>(ny, i, j) > 0.0 ? 1 : -1);
					for(int k = 0; k < 3; ++k)
					{
						c[k] = index2D<double>(phi, i+k-1, j);
						for(int l = j + dir; l != j + HALF_STENCIL_HEIGHT * dir; l += dir)
						{
							if(l < 0 || l >= phi->dimensions[1] || (index2D<long>(mask, i+k-1, l) & 1) == 0 ||
								(index2D<double>(phi, i+k-1, l) <= EPSILON))
								break;

							c[k] += index2D<double>(phi, i+k-1, l);
						}
						for(int l = j - dir; l != j - HALF_STENCIL_HEIGHT * dir; l -= dir)
						{
							if(l < 0 || l >= phi->dimensions[1] || (index2D<long>(mask, i+k-1, l) & 1) == 0 ||
								(index2D<double>(phi, i+k-1, l) >= 1.0 - EPSILON))
								break;

							c[k] += index2D<double>(phi, i+k-1, l) - 1.0;
						}
						c[k] *= dy;
					}
					if(i == 1) c[0] = c[1];
					if(i == phi->dimensions[0] - 2) c[2] = c[1];

					double a = (c[0] - 2.0 * c[1] + c[2]) / (dx * dx);
					double b = (c[2] - c[0]) / (2.0 * dx);

					// dx, heights, cx, cy, r)
					double len = sqrt(1.0 + b * b);
					// Initial approximation. Overwrite if better is found.
					index2D<double>(kappa, i, j) = -a / (len * len * len);
					// If radius is less than REFINE_RADIUS times greater than the 
					// grid spacing, refine the curvature.
					if(refine && (len * len * len < REFINE_RADIUS * dx * vof::abs(a)))
					{
						double s = vof::sign(a);
						if(s > 0.0)
							for(int i = 0; i < 3; ++i)
								c[i] = -c[i];

						double n_x = -s * b / len, n_y = 1.0 / len;
						double radius = (len * len * len) / vof::abs(a);
						double cx = -radius * n_x;
						double cy = c[1] - radius * n_y;

						if(refine_curvature(dx, c, cx, cy, radius))
							index2D<double>(kappa, i, j) = -s / radius;
					}
				}
			}
			else
				index2D<double>(kappa, i, j) = 0.0;
		}
	}
}

static char curvature_doc[] = "kappa = curvature(phi, mask, delta[, refine, nx, ny])\n"
	"\n"
	"Calculates the curvature of the interface in the volume-of-fluid\n"
	"method. 'phi' defines how much of each cell is filled with fluid.\n"
	"Zero indicates an empty cell, one a full cell.\n"
	"If the value is in the range (0, 1), an interface will pass\n"
	"through the cell. (nx, ny) is the normal of the interface in\n"
	"such cells. The returned value 'kappa' is 1e308 in cells with\n"
	"'phi'=0 or 'phi'=1 except when the interface passes along the edge\n"
	"between two cells. In this case, the curvature is calculated in at\n"
	"least one of the two bordering cells. 'kappa' is undefined in ghost\n"
	"cells.\n"
	"\n"
	"  Input:\n"
	"  - phi    = array of volume-of-fluid values\n"
	"  - mask   = array of flags, 1 for fluid cells, 0 for obstacle cells\n"
	"  - delta  = grid spacing dx and dy\n"
	"  - refine = if true, the curvature is iteratively refined\n"
	"  - nx     = array of normal x-components (optional)\n"
	"  - ny     = array of normal y-components (optional)\n"
	"  Output:\n"
	"  - kappa  = array of curvature values\n";

static PyObject* curvature(PyObject* self, PyObject* args)
{
	PyObject* result = NULL;
	PyObject* phi_ = NULL;
	PyObject* mask_ = NULL;
	PyObject* delta_ = NULL;
	PyObject* nx_ = NULL;
	PyObject* ny_ = NULL;
	PyArrayObject* phi = NULL;
	PyArrayObject* mask = NULL;
	PyArrayObject* delta = NULL;
	PyArrayObject* nx = NULL;
	PyArrayObject* ny = NULL;
	PyArrayObject* kappa = NULL;
	int refine;

	if(!PyArg_ParseTuple(args, "OOO|iOO:curvature", &phi_, &mask_, &delta_, &refine, &nx_, &ny_))
		return NULL;

	bool failed = false;
	if(!failed && !(phi = (PyArrayObject*)PyArray_FromObject(phi_, PyArray_DOUBLE, 2, 2)))
		failed = true;
	if(!failed && !(mask = (PyArrayObject*)PyArray_FromObject(mask_, PyArray_LONG, 2, 2)))
		failed = true;
	if(!failed && (phi->dimensions[0] != mask->dimensions[0] || phi->dimensions[1] != mask->dimensions[1]))
	{
		PyErr_Format(PyExc_TypeError, "'mask' and 'phi' must have the same size.");
		failed = true;
	}

	if(!failed && !(delta = (PyArrayObject*)PyArray_FromObject(delta_, PyArray_DOUBLE, 1, 1)))
		failed = true;
	if(!failed && delta->dimensions[0] != 2) 
	{
		PyErr_Format(PyExc_TypeError, "'delta' must have two elements.");
		failed = true;
	}

	if(!failed && !(kappa = create_similar_array(phi)))
		failed = true;

	if(!failed)
	{
		if(nx_ != NULL && ny_ != NULL)
		{
			failed |= !(nx = (PyArrayObject*)PyArray_FromObject(nx_, PyArray_DOUBLE, 2, 2));
			failed |= !(ny = (PyArrayObject*)PyArray_FromObject(ny_, PyArray_DOUBLE, 2, 2));
			if(!failed && (phi->dimensions[0] != nx->dimensions[0] || phi->dimensions[1] != nx->dimensions[1] ||
				phi->dimensions[0] != ny->dimensions[0] || phi->dimensions[1] != ny->dimensions[1]))
			{
				PyErr_Format(PyExc_TypeError, "'phi', 'nx' and 'ny' must have the same size.");
				failed = true;
			}
		}
		else
		{
			failed |= !(nx = create_similar_array(phi));
			failed |= !(ny = create_similar_array(phi));
			if(!failed)
			// gradient() is cheaper than reconstruct() and for this purpose, it is just as good.
				//gradient_main(phi, mask, delta, nx, ny);
				reconstruct_main(phi, mask, delta, nx, ny, kappa); // kappa is used as a dummy here
		}
	}

	if(!failed)
	{
		curvature_main(phi, mask, delta, refine != 0, nx, ny, kappa);
		result = Py_BuildValue("O", kappa);
	}
	Py_XDECREF(phi);
	Py_XDECREF(mask);
	Py_XDECREF(delta);
	Py_XDECREF(nx);
	Py_XDECREF(ny);
	Py_XDECREF(kappa);
	return result;
}

// phi_temp is assumed to be a copy of phi on entry
static void advect_horizontally(PyArrayObject* phi, PyArrayObject* u, PyArrayObject* v, 
								PyArrayObject* nx, PyArrayObject* ny, PyArrayObject* c, 
								PyArrayObject* phi_temp, double dt, double dx, double dy)
{
	double vel, flux, phi_ij, nx_ij, ny_ij;
	for(int j = 1; j < u->dimensions[1] - 1; ++j)
		for(int i = 0; i < u->dimensions[0]; ++i)
		{
			vel = index2D<double>(u, i, j);
			// If the velocity is greater than zero, calculate the flux from the 
			// left cell, otherwise from the right cell.
			if(vel > 0.0)
			{
				phi_ij = index2D<double>(phi, i, j);
				nx_ij = index2D<double>(nx, i, j);
				ny_ij = index2D<double>(ny, i, j);
				if(phi_ij >= 1.0 || phi_ij <= 0.0 || (nx_ij == 0.0 && ny_ij == 0.0))
					flux = phi_ij * vel * dt / dx;
				else
					flux = rightFlux(nx_ij, ny_ij, index2D<double>(c, i, j), dx, dy, vel * dt);
			}
			else
			{
				phi_ij = index2D<double>(phi, i+1, j);
				nx_ij = index2D<double>(nx, i+1, j);
				ny_ij = index2D<double>(ny, i+1, j);
				if(phi_ij >= 1.0 || phi_ij <= 0.0 || (nx_ij == 0.0 && ny_ij == 0.0))
					flux = phi_ij * vel * dt / dx;
				else
					flux = -leftFlux(nx_ij, ny_ij, index2D<double>(c, i+1, j), dx, dy, -vel * dt);
			}

			// Do not change the value in ghost cells
			if(i > 0)
				index2D<double>(phi_temp, i, j) -= flux;
			if(i < u->dimensions[0] - 1)
				index2D<double>(phi_temp, i+1, j) += flux;
		}
}

// phi_temp is assumed to be a copy of phi on entry
static void advect_vertically(PyArrayObject* phi, PyArrayObject* u, PyArrayObject* v, 
							  PyArrayObject* nx, PyArrayObject* ny, PyArrayObject* c, 
							  PyArrayObject* phi_temp, double dt, double dx, double dy)
{
	double vel, flux, phi_ij, nx_ij, ny_ij;
	for(int j = 0; j < v->dimensions[1]; ++j)
		for(int i = 1; i < v->dimensions[0] - 1; ++i)
		{
			vel = index2D<double>(v, i, j);
			// If the velocity is greater than zero, calculate the flux from the 
			// bottom cell, otherwise from the top cell.
			if(vel > 0.0)
			{
				phi_ij = index2D<double>(phi, i, j);
				nx_ij = index2D<double>(nx, i, j);
				ny_ij = index2D<double>(ny, i, j);
				if(phi_ij >= 1.0 || phi_ij <= 0.0 || (nx_ij == 0.0 && ny_ij == 0.0))
					flux = phi_ij * vel * dt / dy;
				else
					flux = topFlux(nx_ij, ny_ij, index2D<double>(c, i, j), dx, dy, vel * dt);
			}
			else
			{
				phi_ij = index2D<double>(phi, i, j+1);
				nx_ij = index2D<double>(nx, i, j+1);
				ny_ij = index2D<double>(ny, i, j+1);
				if(phi_ij >= 1.0 || phi_ij <= 0.0 || (nx_ij == 0.0 && ny_ij == 0.0))
					flux = phi_ij * vel * dt / dy;
				else
					flux = -bottomFlux(nx_ij, ny_ij, index2D<double>(c, i, j+1), dx, dy, -vel * dt);
			}

			// Do not change the value in ghost cells
			if(j > 0)
				index2D<double>(phi_temp, i, j) -= flux;
			if(j < v->dimensions[1] - 1)
				index2D<double>(phi_temp, i, j+1) += flux;
		}
}

// phi is overwritten, no other arrays are modified
static void advect_main(PyArrayObject* phi, PyArrayObject* u, PyArrayObject* v, 
						PyArrayObject* mask, PyArrayObject* delta, 
						PyArrayObject* nx, PyArrayObject* ny, PyArrayObject* c, double dt, int dir)
{
	PyArrayObject* phi_temp = (PyArrayObject*)PyArray_Copy(phi);
	PyArrayObject* nx_temp = create_similar_array(nx);
	PyArrayObject* ny_temp = create_similar_array(ny);
	PyArrayObject* c_temp = create_similar_array(c);

	double dx = index1D<double>(delta, 0), dy = index1D<double>(delta, 1);
	if((dir & 1) == 0)
	{
		// first advect horizontally
		advect_horizontally(phi, u, v, nx, ny, c, phi_temp, dt, dx, dy);

		// scale according to sussman's article
		for(int j = 1; j < phi_temp->dimensions[1] - 1; ++j)
			for(int i = 1; i < phi_temp->dimensions[0] - 1; ++i)
			{
				index2D<double>(phi_temp, i, j) /= (1.0 - dt / dx * (index2D<double>(u, i, j) - index2D<double>(u, i-1, j)));
				index2D<double>(phi, i, j) = index2D<double>(phi_temp, i, j);
			}
	}
	else
	{
		// first advect vertically
		advect_vertically(phi, u, v, nx, ny, c, phi_temp, dt, dx, dy);

		// scale according to sussman's article
		for(int j = 1; j < phi_temp->dimensions[1] - 1; ++j)
			for(int i = 1; i < phi_temp->dimensions[0] - 1; ++i)
			{
				index2D<double>(phi_temp, i, j) /= (1.0 - dt / dy * (index2D<double>(v, i, j) - index2D<double>(v, i, j-1)));
				index2D<double>(phi, i, j) = index2D<double>(phi_temp, i, j);
			}
	}

	// reconstruct interface
	reconstruct_main(phi_temp, mask, delta, nx_temp, ny_temp, c_temp);

	if((dir & 1) == 0)
	{
		// second advect vertically
		advect_vertically(phi_temp, u, v, nx_temp, ny_temp, c_temp, phi, dt, dx, dy);

		// adjust according to sussman's article
		for(int j = 1; j < phi_temp->dimensions[1] - 1; ++j)
			for(int i = 1; i < phi_temp->dimensions[0] - 1; ++i)
			{
				index2D<double>(phi, i, j) += index2D<double>(phi_temp, i, j) * 
					dt / dy * (index2D<double>(v, i, j) - index2D<double>(v, i, j-1));
				if(index2D<double>(phi, i, j) < EPSILON)
					index2D<double>(phi, i, j) = 0.0;
				else if(index2D<double>(phi, i, j) > 1.0 - EPSILON)
					index2D<double>(phi, i, j) = 1.0;
			}
	}
	else
	{
		// second advect horizontally
		advect_horizontally(phi_temp, u, v, nx_temp, ny_temp, c_temp, phi, dt, dx, dy);

		// adjust according to sussman's article
		for(int j = 1; j < phi_temp->dimensions[1] - 1; ++j)
			for(int i = 1; i < phi_temp->dimensions[0] - 1; ++i)
			{
				index2D<double>(phi, i, j) += index2D<double>(phi_temp, i, j) * 
					dt / dx * (index2D<double>(u, i, j) - index2D<double>(u, i-1, j));
				if(index2D<double>(phi, i, j) < EPSILON)
					index2D<double>(phi, i, j) = 0.0;
				else if(index2D<double>(phi, i, j) > 1.0 - EPSILON)
					index2D<double>(phi, i, j) = 1.0;
			}
	}

	Py_XDECREF(phi_temp);
	Py_XDECREF(nx_temp);
	Py_XDECREF(ny_temp);
	Py_XDECREF(c_temp);
}

static char advect_doc[] = "phi = advect(phi, u, v, mask, delta, dt, dir[, nx, ny, c])\n"
	"\n"
	"Advects the volume-of-fluid field in place. The field is advected in two\n"
	"steps. If 'dir' is zero, the field is first advected horizontally, then\n"
	"vertically. If 'dir' is one, the order is opposite. 'dir' should alternate\n"
	"between zero and one for each iteration for higher accuracy.\n"
	"\n"
	"The normal (nx, ny) is assumed to be normalized (or zero).\n"
	"\n"
	"  Input/output:\n"
	"  - phi = array of volume-of-fluid values\n"
	"  Input:\n"
	"  - u     = array of cell face velocities, x-component\n"
	"  - v     = array of cell face velocities, y-component\n"
	"  - mask  = array of flags, 1 for fluid cells, 0 for obstacles\n"
	"  - delta = grid spacing, dx and dy\n"
	"  - dt    = time step value\n"
	"  - dir   = which direction to advect first\n"
	"  - nx    = array of normal x-components (optional)\n"
	"  - ny    = array of normal y-components (optional)\n"
	"  - c     = array of interface coefficients (optional)\n";

static PyObject* advect(PyObject* self, PyObject* args)
{
	PyObject* phi_ = NULL;
	PyObject* u_ = NULL;
	PyObject* v_ = NULL;
	PyObject* mask_ = NULL;
	PyObject* delta_ = NULL;
	PyObject* nx_ = NULL;
	PyObject* ny_ = NULL;
	PyObject* c_ = NULL;

	double dt = 0.0;
	int dir = 0;
	PyObject* result = NULL;
	PyArrayObject* phi = NULL;
	PyArrayObject* u = NULL;
	PyArrayObject* v = NULL;
	PyArrayObject* mask = NULL;
	PyArrayObject* delta = NULL;
	PyArrayObject* nx = NULL;
	PyArrayObject* ny = NULL;
	PyArrayObject* c = NULL;

	//phi, u, v, mask, delta, dt, dir[, nx, ny, c]
	if(!PyArg_ParseTuple(args, "OOOOOdi|OOO:advect", &phi_, &u_, &v_, &mask_, &delta_, &dt, &dir, &nx_, &ny_, &c_))
		return NULL;

	bool failed = false;
	if(!failed && !(phi = (PyArrayObject*)PyArray_FromObject(phi_, PyArray_DOUBLE, 2, 2)))
		failed = true;
	if(!failed && !(u = (PyArrayObject*)PyArray_FromObject(u_, PyArray_DOUBLE, 2, 2)))
		failed = true;
	if(!failed && (phi->dimensions[0] != u->dimensions[0] + 1 || phi->dimensions[1] != u->dimensions[1]))
	{
		PyErr_Format(PyExc_TypeError, "'u' must have the same height as 'phi' and width one less than 'phi'.");
		failed = true;
	}

	if(!failed && !(v = (PyArrayObject*)PyArray_FromObject(v_, PyArray_DOUBLE, 2, 2)))
		failed = true;
	if(!failed && (phi->dimensions[0] != v->dimensions[0] || phi->dimensions[1] != v->dimensions[1] + 1))
	{
		PyErr_Format(PyExc_TypeError, "'v' must have the same width as 'phi' and height one less than 'phi'.");
		failed = true;
	}

	if(!failed && !(mask = (PyArrayObject*)PyArray_FromObject(mask_, PyArray_LONG, 2, 2)))
		failed = true;
	if(!failed && (phi->dimensions[0] != mask->dimensions[0] || phi->dimensions[1] != mask->dimensions[1]))
	{
		PyErr_Format(PyExc_TypeError, "'mask' and 'phi' must have the same size.");
		failed = true;
	}

	if(!failed && !(delta = (PyArrayObject*)PyArray_FromObject(delta_, PyArray_DOUBLE, 1, 1)))
		failed = true;
	if(!failed && delta->dimensions[0] != 2) 
	{
		PyErr_Format(PyExc_TypeError, "'delta' must have two elements.");
		failed = true;
	}

	if(!failed)
	{
		if(nx_ != NULL && ny_ != NULL && c_ != NULL)
		{
			failed |= !(nx = (PyArrayObject*)PyArray_FromObject(nx_, PyArray_DOUBLE, 2, 2));
			failed |= !(ny = (PyArrayObject*)PyArray_FromObject(ny_, PyArray_DOUBLE, 2, 2));
			failed |= !(c = (PyArrayObject*)PyArray_FromObject(c_, PyArray_DOUBLE, 2, 2));
			if(!failed && (phi->dimensions[0] != nx->dimensions[0] || phi->dimensions[1] != nx->dimensions[1] ||
				phi->dimensions[0] != ny->dimensions[0] || phi->dimensions[1] != ny->dimensions[1] ||
				phi->dimensions[0] != c->dimensions[0] || phi->dimensions[1] != c->dimensions[1])) 
			{
				PyErr_Format(PyExc_TypeError, "'phi', 'nx', 'ny' and 'c' must have the same size.");
				failed = true;
			}
		}
		else
		{
			failed |= !(nx = create_similar_array(phi));
			failed |= !(ny = create_similar_array(phi));
			failed |= !(c = create_similar_array(phi));
			if(!failed)
				reconstruct_main(phi, mask, delta, nx, ny, c);
		}
	}

	if(!failed)
	{
		advect_main(phi, u, v, mask, delta, nx, ny, c, dt, dir);
		result = Py_BuildValue("O", phi);
	}

	//phi, u, v, mask, delta, dt, dir[, nx, ny, c]
	Py_XDECREF(phi);
	Py_XDECREF(u);
	Py_XDECREF(v);
	Py_XDECREF(mask);
	Py_XDECREF(delta);
	Py_XDECREF(nx);
	Py_XDECREF(ny);
	Py_XDECREF(c);
	return result;
}

static void circle_fraction_main(double cx, double cy, double r, 
								 PyArrayObject* x, PyArrayObject* y, PyArrayObject* f)
{
	for(int j = 0; j < y->dimensions[0] - 1; ++j)
	{
		double y1 = index1D<double>(y, j);
		double y2 = index1D<double>(y, j + 1);
		for(int i = 0; i < x->dimensions[0] - 1; ++i)
		{
			double x1 = index1D<double>(x, i);
			double x2 = index1D<double>(x, i + 1);
			if(x1 >= cx + r || x2 <= cx - r || y1 >= cy + r || y2 <= cy - r)
				index2D<double>(f, i, j) = 0.0;
			else
				index2D<double>(f, i, j) = rectCircleIntersection(x1, y1, x2, y2, cx, cy, r) / 
					((x2 - x1) * (y2 - y1));
		}
	}
}

static char circle_fraction_doc[] = "f = circle_fraction(cx, cy, r, x, y)\n"
	"\n"
	"Calculate the fraction of a cell inside the circle with centre\n"
	"(cx, cy) and radius 'r'. 'x' and 'y' are 1D arrays with the x- and\n"
	"y-coordinate of the cell corners respectively.\n"
	"\n"
	"  Input:\n"
	"  - cx = x-coordinate of the circle centre\n"
	"  - cy = y-coordinate of the circle centre\n"
	"  - r  = radius of the circle\n"
	"  - x  = array of grid node x-coordinates\n"
	"  - y  = array of grid node y-coordinates\n"
	"  Output:\n"
	"  - f  = array of fractions for each cell.\n";

static PyObject* circle_fraction(PyObject* self, PyObject* args)
{
	PyObject* result = NULL;
	PyObject* x_;
	PyObject* y_;
	PyArrayObject* x = NULL;
	PyArrayObject* y = NULL;
	PyArrayObject* f = NULL;
	double cx, cy, r;

	if(!PyArg_ParseTuple(args, "dddOO:circle_fraction", &cx, &cy, &r, &x_, &y_))
		return NULL;

	bool failed = false;
	if(!failed && !(x = (PyArrayObject*)PyArray_FromObject(x_, PyArray_DOUBLE, 1, 1)))
		failed = true;
	if(!failed && !(y = (PyArrayObject*)PyArray_FromObject(y_, PyArray_DOUBLE, 1, 1)))
		failed = true;
	
	if(!failed)
	{
		int dims[] = {(int)x->dimensions[0] - 1, (int)y->dimensions[0] - 1};
		failed = !(f = (PyArrayObject*)PyArray_FromDims(2, dims, PyArray_DOUBLE));
	}

	if(!failed)
	{
		circle_fraction_main(cx, cy, r, x, y, f);
		result = Py_BuildValue("O", f);
	}

	Py_XDECREF(x);
	Py_XDECREF(y);
	Py_XDECREF(f);
	return result;
}

static void rectangle_fraction_main(double left, double bottom, double right, double top,
									PyArrayObject* x, PyArrayObject* y, PyArrayObject* f)
{
	for(int j = 0; j < y->dimensions[0] - 1; ++j)
	{
		double y1 = index1D<double>(y, j);
		double y2 = index1D<double>(y, j + 1);
		for(int i = 0; i < x->dimensions[0] - 1; ++i)
		{
			double x1 = index1D<double>(x, i);
			double x2 = index1D<double>(x, i + 1);
			if(x1 >= right || x2 <= left || y1 >= top || y2 <= bottom)
				index2D<double>(f, i, j) = 0.0;
			else
				index2D<double>(f, i, j) = rectRectIntersection(x1, y1, x2, y2, 
					left, bottom, right, top) / ((x2 - x1) * (y2 - y1));
		}
	}
}

static char rectangle_fraction_doc[] = "f = rectangle_fraction(cx, cy, r, x, y)\n"
	"\n"
	"Calculate the fraction of a cell inside the rectangle with lower left\n"
	"corner (x1, y1) and upper right corner (x2, y2). 'x' and 'y' are 1D \n"
	"arrays with the x- and y-coordinate of the cell corners respectively.\n"
	"\n"
	"  Input:\n"
	"  - x1 = x-coordinate of the lower left corner of the rectangle\n"
	"  - y1 = y-coordinate of the lower left corner of the rectangle\n"
	"  - x2 = x-coordinate of the upper right corner of the rectangle\n"
	"  - y2 = y-coordinate of the upper right corner of the rectangle\n"
	"  - x  = array of grid node x-coordinates\n"
	"  - y  = array of grid node y-coordinates\n"
	"  Output:\n"
	"  - f  = array of fractions for each cell.\n";

static PyObject* rectangle_fraction(PyObject* self, PyObject* args)
{
	PyObject* result = NULL;
	PyObject* x_;
	PyObject* y_;
	PyArrayObject* x = NULL;
	PyArrayObject* y = NULL;
	PyArrayObject* f = NULL;
	double x1, y1, x2, y2;

	if(!PyArg_ParseTuple(args, "ddddOO:rectangle_fraction", &x1, &y1, &x2, &y2, &x_, &y_))
		return NULL;

	bool failed = false;
	if(!failed && !(x = (PyArrayObject*)PyArray_FromObject(x_, PyArray_DOUBLE, 1, 1)))
		failed = true;
	if(!failed && !(y = (PyArrayObject*)PyArray_FromObject(y_, PyArray_DOUBLE, 1, 1)))
		failed = true;
	
	if(!failed)
	{
		int dims[] = {(int)x->dimensions[0] - 1, (int)y->dimensions[0] - 1};
		failed = !(f = (PyArrayObject*)PyArray_FromDims(2, dims, PyArray_DOUBLE));
	}

	if(!failed)
	{
		rectangle_fraction_main(x1, y1, x2, y2, x, y, f);
		result = Py_BuildValue("O", f);
	}

	Py_XDECREF(x);
	Py_XDECREF(y);
	Py_XDECREF(f);
	return result;
}



static PyMethodDef module_methods[] = {
	{"gradient", gradient, METH_VARARGS, gradient_doc},
	{"reconstruct", reconstruct, METH_VARARGS, reconstruct_doc},
	{"curvature", curvature, METH_VARARGS, curvature_doc},
	{"advect", advect, METH_VARARGS, advect_doc},
	{"circle_fraction", circle_fraction, METH_VARARGS, circle_fraction_doc},
	{"rectangle_fraction", rectangle_fraction, METH_VARARGS, rectangle_fraction_doc},
	{NULL, NULL}
};

static char module_doc[] = "module vof:\n"
	"  gx, gy = gradient(phi, mask, delta)\n"
	"  nx, ny, c = reconstruct(phi, delta)\n"
	"  kappa = curvature(phi, mask, delta[, nx, ny])\n"
	"  phi = advect(phi, u, v, mask, delta, dt, dir[, nx, ny, c])\n"
	"  f = circle_fraction(cx, cy, r, x, y)\n"
	"  f = rectangle_fraction(x1, y1, x2, y2, x, y)\n";

#ifdef __cplusplus
extern "C" PyMODINIT_FUNC initvof();
#endif

PyMODINIT_FUNC initvof()
{
	Py_InitModule3("vof", module_methods, module_doc);
	import_array();
}
