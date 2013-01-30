#include <Python.h>
#include <numpy/arrayobject.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

//#define DEBUG
//#undef NDEBUG
#include <cassert>

#include "sparse.h"

// Include SuperLU. The next three files must be included in this order.
#include "slu_ddefs.h"
#include "supermatrix.h"
#include "slu_util.h"

template<class T>
T& index(PyArrayObject* v, const array<int>& idx, int dim = 0, int offs = 0)
{
	assert(v && "Trying to index NULL.");
	assert((int)idx.size() == v->nd && "Incorrect index size.");

	char* p = v->data;
	for(int i = 0; i < (int)idx.size(); ++i)
	{
		assert(idx[i] >= 0 && idx[i] < v->dimensions[i] && "Index out of bounds.");
		p += idx[i] * v->strides[i];
	}
	p += offs * v->strides[dim];
	return *((T*)p);
}

// A cast from void* to long causes a warning in Visual Studio 2005, so
// disable the warning first.
inline long ptr2long(void* p)
{
#pragma warning(disable:4311)
	return (long)p;
#pragma warning(default:4311)
}

static void build_matrix(SparseMatrix<double>& A, int dims, npy_intp* sizes, 
	double* deltas, PyArrayObject* mask, PyArrayObject* rho)
{
	array<int> idx(dims, 1);
	array<int> offsets(dims);
	offsets.back() = 1;
	for(int i = dims - 1; i > 0; --i)
		offsets[i-1] = offsets[i] * (sizes[i]-2);
	int prod = offsets[0] * (sizes[0]-2);
	A.clear(prod, prod, prod * (2 * dims + 1));
	double rho_avg;

	// for each line...
	int row = 0;
	while(idx[0] < sizes[0]-1)
	{
		assert(row < prod && "Row number has exceeded matrix size.");
		if(index<long>(mask, idx, 0, 0) == 0)
			A.append(row, row, 1.0); // inactive cell
		else
		{
			for(int j = 0; j < dims; ++j)
				if(idx[j] > 1 && index<long>(mask, idx, j, -1) != 0 &&
					(idx[j] < sizes[j]-2 || index<long>(mask, idx, j, 1) != 0))
				{
					rho_avg = 0.5 * (index<double>(rho, idx) + index<double>(rho, idx, j, -1));
					A.append(row, row - offsets[j], -1.0/(rho_avg*deltas[j]*deltas[j]));
				}

			double diag = 0.0;
			for(int j = 0; j < dims; ++j)
			{
				// dirichlet boundary condition:
				// constant on boundary => (phi[0] + phi[1])/2 = c
				// phi[0] - 2phi[1] + phi[2] = -f*h^2 => 3phi[1] - phi[2] = f*h^2 + 2c
				double val = 0.0;
				if(!(idx[j] == 1 && index<long>(mask, idx, j, -1) == 0) &&
				   !(idx[j] == sizes[j]-2 && index<long>(mask, idx, j, 1) == 0))
				{
					rho_avg = 0.5 * (index<double>(rho, idx) + index<double>(rho, idx, j, -1));
					if(idx[j] == 1)
						val = index<long>(mask, idx, j, -1) - 1;
					else
						val = index<long>(mask, idx, j, -1) & 1;
					diag += val/(rho_avg*deltas[j]*deltas[j]);

					rho_avg = 0.5 * (index<double>(rho, idx) + index<double>(rho, idx, j, 1));
					if(idx[j] == sizes[j]-2)
						val = index<long>(mask, idx, j, 1) - 1;
					else
						val = index<long>(mask, idx, j, 1) & 1;
					diag += val/(rho_avg*deltas[j]*deltas[j]);

					if(index<long>(mask, idx) == 3)
						diag += 1.0/(index<double>(rho, idx)*deltas[j]*deltas[j]); // locked cell
				}
			}
			A.append(row, row, diag);

			for(int j = dims - 1; j >= 0; --j)
				if(idx[j] < sizes[j]-2 && index<long>(mask, idx, j, 1) != 0 &&
				   (idx[j] > 1 || index<long>(mask, idx, j, -1) != 0))
				{
					rho_avg = 0.5 * (index<double>(rho, idx) + index<double>(rho, idx, j, 1));
					A.append(row, row + offsets[j], -1.0/(rho_avg*deltas[j]*deltas[j]));
				}
		}

		int i = (int)idx.size() - 1;
		++idx[i];
		while(i > 0 && idx[i] >= sizes[i]-1)
		{
			idx[i--] = 1;
			++idx[i];
		}
		++row;
	}
}

struct SUPERLU_DATA
{
	SUPERLU_DATA(int height, int width) : Pc(height), Pr(width), 
		etree(height), Dc(height), Dr(width), equed(0)
	{
		memset(&L, 0, sizeof(L));
		memset(&U, 0, sizeof(U));
	}

	~SUPERLU_DATA()
	{
		freeLU();
	}

	void freeLU()
	{
		if(L.Store)
			Destroy_SuperNode_Matrix(&L);
		if(U.Store)
			Destroy_CompCol_Matrix(&U);
		memset(&L, 0, sizeof(L));
		memset(&U, 0, sizeof(U));
		equed = 0;
	}

	void resize(int height, int width)
	{
		Pc.resize(height);
		Pr.resize(width);
		Dc.resize(height);
		Dr.resize(width);
		etree.resize(height);
	}

	array<int> Pc, Pr, etree;
	array<double> Dc, Dr;
	SuperMatrix L, U;
	char equed;
};

// Python object for SUPERLU_DATA
struct poisson_data
{
	PyObject_HEAD
	SUPERLU_DATA* data;
};

static PyObject* poisson_data_new(PyObject* self, PyObject* args);
static void poisson_data_dealloc(PyObject* self);
static int poisson_data_compare(PyObject* self, PyObject* other);
static long poisson_data_hash(PyObject* self);
static PyObject* poisson_data_repr(PyObject *self);

static PyTypeObject poisson_data_Type = {
	PyObject_HEAD_INIT(&PyType_Type)
	0,
	"poisson_data",
	sizeof(poisson_data),
	0,
	poisson_data_dealloc,
	NULL,
	NULL,
	NULL,
	poisson_data_compare,
	poisson_data_repr,
	NULL,
	NULL,
	NULL,
	poisson_data_hash,
	NULL,
	NULL,
};

static char poisson_data_doc[] = "data = poisson_data()\n"
	"\n"
	"Create a data structure for holding permutations and factors calculated\n"
	"in poisson(). These can be reused in other calls to poisson().\n";

static PyObject* poisson_data_new(PyObject* self, PyObject* args)
{
	poisson_data* obj = PyObject_New(poisson_data, &poisson_data_Type);
	if(obj)
		obj->data = new SUPERLU_DATA(0, 0);
	return (PyObject*)obj;
}

static void poisson_data_dealloc(PyObject* self)
{
	poisson_data* obj = (poisson_data*)self;
	if(obj)
		delete obj->data;
	PyObject_Del(self);
}

static int poisson_data_compare(PyObject* self, PyObject* other)
{
	poisson_data* obj = (poisson_data*)self;
	PyObject* tmp = PyInt_FromLong(ptr2long(obj->data));
	int result = (tmp ? PyObject_Compare(tmp, other) : 0);
	Py_XDECREF(tmp);
	return result;
}

static long poisson_data_hash(PyObject* self)
{
	poisson_data* obj = (poisson_data*)self;
	PyObject *tmp = PyInt_FromLong(ptr2long(obj->data));
	long result = (tmp ? PyObject_Hash(tmp) : -1);
	Py_XDECREF(tmp);
	return result;
}

static PyObject* poisson_data_repr(PyObject *self)
{
	return PyString_FromString("poisson_data()");
}

// The SuperLU library is not well explained, so some of this code is more or
// less copied from examples.
//
// reuseLevel:
// - 0 = solve from scratch.
// - 1 = use same column permutation (same A pattern).
// - 2 = use same permutations and LU storage (similar A).
// - 3 = use same factorization (same A).
static void superLU(SparseMatrix<double>& A, array<double>& b, array<double>& x, 
	SUPERLU_DATA* previousData = NULL, int reuseLevel = 0)
{
	// 'A' is stored in row major order while SuperLU operates on column major matrices.
	// This is why width and height is swapped in many places.
	SUPERLU_DATA* data;
	if(previousData == NULL)
	{
		data = new SUPERLU_DATA(A.height(), A.width());
		reuseLevel = 0;
	}
	else
	{
		data = previousData;
		if((int)data->Pr.size() != A.width() || (int)data->Pc.size() != A.height())
		{
			data->resize(A.height(), A.width());
			reuseLevel = 0;
		}
	}

	int info;
	double rpg, rcond, ferr, berr; // ?, condition number, forward error bound,  backward error bound
	SuperMatrix slu_A, slu_b, slu_x; // SuperLU versions of A, b and x
	superlu_options_t options;
	mem_usage_t mem_usage;
	SuperLUStat_t stat;
	//memcpy(x.c_array(), b.c_array(), x.size() * sizeof(double));

	dCreate_CompRow_Matrix(&slu_A, A.height(), A.width(), A.length(), 
		A.values(), A.columns(), A.rows(), SLU_NR, SLU_D, SLU_GE);
	dCreate_Dense_Matrix(&slu_b, (int)b.size(), 1, b.c_array(), (int)b.size(), 
		SLU_DN, SLU_D, SLU_GE);
	dCreate_Dense_Matrix(&slu_x, (int)x.size(), 1, x.c_array(), (int)x.size(), 
		SLU_DN, SLU_D, SLU_GE);

    set_default_options(&options);
	options.Equil = YES; // must be tested for time vs. accuracy

	switch(reuseLevel)
	{
	case(1):
		options.Fact = SamePattern;
		data->freeLU();
		break;
	case(2):
		options.Fact = SamePattern_SameRowPerm;
		// use same LU storage
		break;
	case(3):
		options.Fact = FACTORED;
		// use same LU
		break;
	default:
		options.Fact = DOFACT;
		data->freeLU();
		break;
	}

	StatInit(&stat);

	dgssvx(&options, &slu_A, data->Pc.c_array(), data->Pr.c_array(), 
		data->etree.c_array(), &data->equed, data->Dr.c_array(), 
		data->Dc.c_array(), &data->L, &data->U, NULL, 0, &slu_b, &slu_x, &rpg, 
		&rcond, &ferr, &berr, &mem_usage, &stat, &info);

	Destroy_SuperMatrix_Store(&slu_A);
	Destroy_SuperMatrix_Store(&slu_b);
	Destroy_SuperMatrix_Store(&slu_x);
	StatFree(&stat);

	if(previousData == NULL)
		delete data;
}

static bool next_cell(PyArrayObject* phi, array<int>& idx)
{
	++idx.back();
	int dim = (int)idx.size() - 1;
	while(idx[dim] >= phi->dimensions[dim] - 1)
	{
		idx[dim] = 1;
		if(--dim < 0)
			break;
		++idx[dim];
	}
	return (dim < 0);
}

static void poisson_main(PyArrayObject* phi, PyArrayObject* f, PyArrayObject* delta, 
	PyArrayObject* mask, PyArrayObject* rho, SUPERLU_DATA* data = NULL, int reuseLevel = 0)
{
	int dims = phi->nd;
	for(int i = 0; i < dims; ++i)
		if(phi->dimensions[i] <= 2)
			return;
	array<int> idx(phi->nd, 1); // start by indexing cell (1,1,1)
	array<double> x, b;
	SparseMatrix<double> A;
	double* deltas = new double[dims];
	for(int i = 0; i < dims; ++i)
		deltas[i] = *((double*)(delta->data + i * delta->strides[0]));

	// vectorize 'phi' and 'f' (-> 'x' and 'b')
	for(;;)
	{
		if(index<long>(mask, idx) == 0)
		{
			// inactive. keep current pressure.
			b.push_back(0.0);
			x.push_back(0.0);
		}
		else
		{
			b.push_back(index<double>(f, idx));
			x.push_back(index<double>(phi, idx)); // use the given phi as initial guess

			// correct 'b' to eliminate boundary values
			for(int i = 0; i < dims; ++i)
			{
				if(idx[i] == 1 && index<long>(mask, idx, i, -1) != 0)
				{
					double rho_avg = 0.5 * (index<double>(rho, idx) +
						index<double>(rho, idx, i, -1));
					--idx[i];
					b.back() += index<double>(phi, idx) / (rho_avg * deltas[i] * deltas[i]);
					++idx[i];
				}
				if(idx[i] == phi->dimensions[i]-2 && index<long>(mask, idx, i, 1) != 0)
				{
					double rho_avg = 0.5 * (index<double>(rho, idx) + 
						index<double>(rho, idx, i, 1));
					++idx[i];
					b.back() += index<double>(phi, idx) / (rho_avg * deltas[i] * deltas[i]);
					--idx[i];
				}
			}
		}

		// update index
		//int dim = phi->nd - 1;
		//++idx[dim];
		//while(dim > 0 && idx[dim] >= phi->dimensions[dim] - 1)
		//{
		//	idx[dim--] = 1;
		//	++idx[dim];
		//}
		if(next_cell(phi, idx))
			break;
	}

	build_matrix(A, dims, phi->dimensions, deltas, mask, rho);

	// for debugging:
	//static bool printed = false;
	//if(!printed)
	//	std::cout << A << std::endl;
	//printed = true;


	// solve linear equation system
	superLU(A, b, x, data, reuseLevel);

	//idx.clear();
	//idx.resize(dims, 1);
	int element = 0;
	for(;;)
	{
		index<double>(phi, idx) = x[element];
		// correct ghost cells
		for(int i = 0; i < dims; ++i)
		{
			if(idx[i] == 1 && index<long>(mask, idx, i, -1) != 0)
			{
				--idx[i];
				index<double>(phi, idx) -= (index<long>(mask, idx, i, 0)-2)*x[element];
				++idx[i];
			}
			if(idx[i] == phi->dimensions[i]-2 && index<long>(mask, idx, i, 1) != 0)
			{
				++idx[i];
				index<double>(phi, idx) -= (index<long>(mask, idx, i, 0)-2)*x[element];
				--idx[i];
			}
		}
		++element;
		//int dim = phi->nd - 1;
		//++idx[dim];
		//while(dim > 0 && idx[dim] >= phi->dimensions[dim] - 1)
		//{
		//	idx[dim--] = 1;
		//	++idx[dim];
		//}
		if(next_cell(phi, idx))
			break;
	}
}

static char poisson_doc[] = "phi = poisson(phi, f, delta, mask, rho[, data, reuse])\n"
	"\n"
	"Solve the Poisson equation \"div(grad(phi)/rho) = -f\".\n"
	"\n"
	"The ghost nodes of 'phi' contains information about the pressure boundary\n"
	"condition value. The ghost nodes of 'mask' contains information about the\n"
	"pressure boundary condition type:\n"
	"  mask = 1  :  phi_input = (phi_0 - phi_1)\n"
	"  mask = 2  :  phi_input = (phi_0)\n"
	"  mask = 3  :  phi_input = (phi_0 + phi_1)\n"
	"phi_0 is the ghost cell node, phi_1 is the node inside the domain.\n"
	"phi_input is the ghost cell node of 'phi' when passed to poisson().\n"
	"phi_0 and phi_1 are the values of 'phi' when returning from the function.\n"
	"\n"
	"For interior cells, mask has the following interpretation:\n"
	"  mask = 0  :  obstacle cell, Neumann boundary condition\n"
	"  mask = 1  :  fluid cell, no boundaries\n"
	"  mask = 3  :  fluid cell, phi fixed to zero\n"
	"\n"
	"  Input/output:\n"
	"  - phi   = solution of equation. The ghost cell values must have been set\n"
	"            on entry.\n"
	"  Input:\n"
	"  - f     = right hand side of equation.\n"
	"  - delta = grid spacing.\n"
	"  - mask  = array of boundary conditions.\n"
	"  - rho   = left hand side coefficient.\n"
	"  - data  = data structure holding permutations and factors that may be\n"
	"            reused (optional).\n"
	"  - reuse = a number in the range [0, 3] indicating how much of previously\n"
	"            calculated data should be reused (optional):\n"
	"            - 0 = calculate from scratch\n"
	"            - 1 = reuse sparsity pattern\n"
	"            - 2 = reuse permutations (similar equation system)\n"
	"            - 3 = reuse all (same equation system)\n";

static PyObject* poisson(PyObject* _self, PyObject* _args)
{
	PyObject* _result = NULL;
	PyObject* _phi;
	PyObject* _f;
	PyObject* _delta;
	PyObject* _mask;
	PyObject* _rho;
	poisson_data* data = NULL;
	PyArrayObject* phi = NULL;
	PyArrayObject* f = NULL;
	PyArrayObject* delta = NULL;
	PyArrayObject* mask = NULL;
	PyArrayObject* rho = NULL;
	int reuse = 0;
	if(!PyArg_ParseTuple(_args, "OOOOO|O!i:poisson", &_phi, &_f, &_delta, &_mask, &_rho, 
		&poisson_data_Type, &data, &reuse))
		return NULL;
	int failed = FALSE;
	if(!failed && !(phi = (PyArrayObject*)PyArray_FromObject(_phi, PyArray_DOUBLE, 0, 0)))
		failed = TRUE;
	if(!failed && !(f = (PyArrayObject*)PyArray_FromObject(_f, PyArray_DOUBLE, 0, 0)))
		failed = TRUE;
	if(!failed && !(delta = (PyArrayObject*)PyArray_FromObject(_delta, PyArray_DOUBLE, 1, 1)))
		failed = TRUE;
	if(!failed && !(mask = (PyArrayObject*)PyArray_FromObject(_mask, PyArray_LONG, 0, 0)))
		failed = TRUE;
	if(!failed && !(rho = (PyArrayObject*)PyArray_FromObject(_rho, PyArray_DOUBLE, 0, 0)))
		failed = TRUE;
	if(!failed)
	{
		failed |= (f->nd != phi->nd) || (f->nd != delta->dimensions[0]) || 
			(f->nd != mask->nd) || (f->nd != rho->nd);
		for(int i = 0; i < f->nd && i < phi->nd; ++i)
		{
			failed |= (f->dimensions[i] != phi->dimensions[i]) || 
				(f->dimensions[i] != mask->dimensions[i]) ||
				(f->dimensions[i] != rho->dimensions[i]);
		}
		if(failed)
			PyErr_Format(PyExc_TypeError, "Incorrect array size");
	}
	if(!failed)
	{
		poisson_main(phi, f, delta, mask, rho, (data ? data->data : NULL), reuse);
		_result = Py_BuildValue("O", phi);
	}
	Py_XDECREF(phi);
	Py_XDECREF(f);
	Py_XDECREF(delta);
	Py_XDECREF(mask);
	Py_XDECREF(rho);
	return _result;
}

static const int FLOOD_FILL_BIT = 0x00010000;

// moves to the next cell unless an obstacle cell or ghost cell has been reached.
static bool remove_singularity_flood_fill_next_cell(PyArrayObject* mask, array<int>& idx, array<int>& dir, array<int>& origin)
{
	idx.back() += dir.back(); // move to next cell
	// wrap index
	int dim = (int)idx.size() - 1;
	while((idx[dim] >= mask->dimensions[dim] - 1) || (idx[dim] < 1) || 
		((index<long>(mask, idx) & 1) == 0))
	{
		if(dir[dim] > 0)
		{
			dir[dim] = -1;
			idx[dim] = origin[dim] - 1;
		}
		else
		{
			dir[dim] = 1;
			idx[dim] = origin[dim];
			if(--dim < 0)
				break;
			idx[dim] += dir[dim];
		}
	}
	return dim < 0;
}

static void remove_singularity_flood_fill(PyArrayObject* mask, array<int>& origin)
{
	// if cell is already visited or is not a fluid cell, return
	if((index<long>(mask, origin) & FLOOD_FILL_BIT) != 0 || (index<long>(mask, origin) & 1) == 0)
		return;

	array<int> idx(origin);
	array<int> dir(idx.size(), 1); // used to keep track of which direction the index is moving
	// Fill order:
	// +---+---+---+
	// | 8 | 2 | 5 |
	// +---+---+---+
	// | 7 | 1 | 4 |
	// +---+---+---+
	// | 9 | 3 | 6 |
	// +---+---+---+

	for(;;)
	{
		index<long>(mask, idx) |= FLOOD_FILL_BIT;
		if(remove_singularity_flood_fill_next_cell(mask, idx, dir, origin))
			break;
	}

	// 'idx' should now be back to 'origin'

	// Iterate through the same cells again, but this time, call this function recursively
	// on each cell's neighbour such that the no connected cell is missed.
	// Most calls will terminate immediately because most cells were marked in the first pass
	// already.
	for(;;)
	{
		// Skip neighbours along the last dimension axis. 
		// There are no missed cells in that direction.
		// Call this function recursively on all neighbouring cell.
		for(int dim = 0; dim < (int)idx.size() - 1; ++dim)
		{
			++idx[dim];
			if(idx[dim] < mask->dimensions[dim] - 1)
				remove_singularity_flood_fill(mask, idx);
			idx[dim] -= 2;
			if(idx[dim] > 0)
				remove_singularity_flood_fill(mask, idx);
			++idx[dim];
		}

		// Find the next cell that was marked in the first pass.
		if(remove_singularity_flood_fill_next_cell(mask, idx, dir, origin))
			break;
	}
}

static bool next_boundary_cell(PyArrayObject* phi, array<int>& idx, int fixed_dim)
{
	int dim = (int)idx.size() - 1;
	if(dim == fixed_dim)
		if(--dim < 0)
			return true;			
	++idx[dim];
	while(idx[dim] >= phi->dimensions[dim] - 1)
	{
		idx[dim] = 1;
		if(--dim == fixed_dim)
			--dim;
		if(dim < 0)
			break;
		++idx[dim];
	}
	return (dim < 0);
}

static void remove_singularity_main(PyArrayObject* mask)
{
	array<int> idx(mask->nd, 1); // N-dimensional index

	// find all pressure nodes that can be determined uniquely,
	// that is, all nodes that are connected to a node with fixed pressure.

	// for each boundary (east, west, north, south etc.)...
	for(int dim = 0; dim < mask->nd; ++dim)
	{
		// iterate through all ghost cells at the high end
		idx[dim] = (int)mask->dimensions[dim] - 2;
		for(;;)
		{
			// if a node has fixed pressure (mask = 2 or mask = 3),
			// mark all connected nodes.
			if((index<long>(mask, idx, dim, 1) & 3) > 1)
				remove_singularity_flood_fill(mask, idx);

			// move to next node on the boundary
			if(next_boundary_cell(mask, idx, dim))
				break;
		}

		// iterate through all ghost cells at the low end
		idx[dim] = 1;
		for(;;)
		{
			if((index<long>(mask, idx, dim, -1) & 3) > 1)
				remove_singularity_flood_fill(mask, idx);
			if(next_boundary_cell(mask, idx, dim))
				break;
		}
	}

	// 'idx' should now be [1, 1, ..., 1]
	// fix pressure nodes that cannot be uniquely determined,
	// that is, those nodes that have not yet been marked.
	for(;;)
	{
		// if there is a fluid cell that has not been marked, fix the pressure.
		if((index<long>(mask, idx) & FLOOD_FILL_BIT) == 0 && (index<long>(mask, idx) & 1) != 0)
		{
			index<long>(mask, idx) |= 2; // fix pressure
			// then mark all nodes connected to it.
			remove_singularity_flood_fill(mask, idx);
		}

		// move to the next cell in the grid
		if(next_cell(mask, idx))
			break;
	}

	// clear visited flag.
	for(;;)
	{
		index<long>(mask, idx) &= ~FLOOD_FILL_BIT; // clear flag
		if(next_cell(mask, idx))
			break;
	}
}

static char remove_singularity_doc[] = "mask = remove_singularity(mask)\n"
	"\n"
	"Fix pressure nodes such that there is a unique solution to the Poisson\n"
	"equation. This will work in any dimension and with any number of closed\n"
	"fluid areas. No pressure nodes must be fixed before calling this function.\n"
	"\n"
	"Input/output:\n"
	"- mask = array of flags (see poisson() for more info).\n";

static PyObject* remove_singularity(PyObject* _self, PyObject* _args)
{
	PyObject* result = NULL;
	PyObject* _mask;
	PyArrayObject* mask = NULL;
	if(!PyArg_ParseTuple(_args, "O:remove_singularity", &_mask))
		return NULL;
	int failed = FALSE;
	if(!failed && !(mask = (PyArrayObject*)PyArray_FromObject(_mask, PyArray_LONG, 0, 0)))
		failed = TRUE;
	if(!failed)
	{
		remove_singularity_main(mask);
		result = Py_BuildValue("O", mask);
	}
	Py_XDECREF(mask);
	return result;
}

static PyMethodDef module_methods[] = {
	{"poisson", poisson, METH_VARARGS, poisson_doc},
	{"remove_singularity", remove_singularity, METH_VARARGS, remove_singularity_doc},
	{"poisson_data", poisson_data_new, METH_VARARGS, poisson_data_doc},
	{NULL, NULL}
};

static char module_doc[] = "module poisson:\n"
	"  phi = poisson(phi, f, delta, mask, rho[, data, reuse])\n"
	"  mask = remove_singularity(mask)\n"
	"  data = poisson_data()\n";

#ifdef __cplusplus
extern "C" PyMODINIT_FUNC initpoisson();
#endif

PyMODINIT_FUNC initpoisson()
{
	Py_InitModule3("poisson", module_methods, module_doc);
	import_array();
}
