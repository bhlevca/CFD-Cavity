#ifndef SPARSE_H
#define SPARSE_H

#include <iostream>
#include <cassert>
#include <cmath>
#include <cstdlib>

// This class should only be used for simple types: 
// float, double, int, bool, char etc.
template<class T>
class array
{
public:
	array(size_t size = 0, const T& value = 0) : m_v(NULL), m_size(size), m_capacity(size)
	{
		if(m_capacity == 0)
			return;
		m_v = new T[m_capacity];
		assert(m_v);
		for(size_t i = 0; i < size; ++i)
			m_v[i] = value;
	}

	array(const array& a) : m_v(NULL), m_size(a.m_size), m_capacity(a.m_size)
	{
		if(m_capacity == 0)
			return;
		m_v = new T[m_capacity];
		assert(m_v);
		for(size_t i = 0; i < m_size; ++i)
			m_v[i] = a.m_v[i];
	}

	~array()
	{
		if(m_v) delete[] m_v;
	}

	T* c_array() {return m_v;}

	size_t size() const {return m_size;}

	bool empty() const {return m_size == 0;}

	void push_back(const T& value)
	{
		if(m_size == m_capacity)
			reserve(m_capacity < MIN_CAPACITY ? MIN_CAPACITY : 2 * m_capacity);
		m_v[m_size] = value;
		++m_size;
	}

	void pop_back()
	{
		assert(m_size > 0);
		--m_size;
	}

	void reserve(size_t capacity)
	{
		if(capacity <= m_capacity)
			return;
		T* v_new = new T[capacity];
		assert(v_new);
		if(m_v)
		{
			for(size_t i = 0; i < m_size; ++i)
				v_new[i] = m_v[i];
			delete[] m_v;
		}
		m_v = v_new;
		m_capacity = capacity;
	}

	void resize(size_t size, const T& value = 0)
	{
		if(size > m_capacity)
			reserve(size);
		for(size_t i = m_size; i < size; ++i)
			m_v[i] = value;
		m_size = size;
	}

	void clear()
	{
		m_size = 0;
	}

	T& operator [] (size_t index) 
	{
		assert(index < m_size); 
		return m_v[index];
	}

	const T& operator [] (size_t index) const
	{
		assert(index < m_size); 
		return m_v[index];
	}

	T& back()
	{
		assert(m_size > 0);
		return m_v[m_size-1];
	}

	const T& back() const
	{
		assert(m_size > 0);
		return m_v[m_size-1];
	}

	T& front()
	{
		assert(m_size > 0);
		return m_v[0];
	}

	const T& front() const
	{
		assert(m_size > 0);
		return m_v[0];
	}

private:
	static const size_t MIN_CAPACITY = 16;
	T* m_v;
	size_t m_size, m_capacity;
};


template<class T>
class SparseMatrix
{
public:
	SparseMatrix(int height = 0, int width = 0, int size = 0) : m_rows(1, 0), m_height(height), m_width(width)
	{
		m_values.reserve(size);
		m_columns.reserve(size);
	}

	void clear(int height = 0, int width = 0, int size = 0)
	{
		if(height > 0) m_height = height;
		if(width > 0) m_width = width;
		m_values.clear();
		m_columns.clear();
		m_rows.clear();
		m_values.reserve(size);
		m_columns.reserve(size);
		m_rows.resize(1, 0);
	}

	int width() const {return m_width;}
	int height() const {return m_height;}
	int length() const {return (int)m_values.size();}

	// compressed row storage
	// returns all values row by row
	T* values() {return m_values.c_array();}
	// returns the columns corresponding to the values
	int * columns() {return m_columns.c_array();}
	// returns an the indices to the first element of each row, 
	// plus an element indicating total length
	int * rows() {return m_rows.c_array();}

	T operator () (int row, int column) const
	{
		if(row > (int)m_rows.size() - 2)
			return T(0);
		int k = search(row, column);
		if(k < m_rows[row+1] && m_columns[k] == column)
			return m_values[k];
		return T(0);
	}

	void append(int row, int column, T value)
	{
		assert(row >= 0 && row < m_height && column >= 0 && column < m_width);
		assert((row == (int)m_rows.size() - 2 && (m_columns.empty() || column > m_columns.back())) || (row > (int)m_rows.size() - 2));
		while(row > (int)m_rows.size() - 2)
			m_rows.push_back((int)m_values.size());
		++m_rows.back();
		m_values.push_back(value);
		m_columns.push_back(column);
	}

	T& operator () (int row, int column)
	{
		if(row > (int)m_rows.size() - 2)
		{
			append(row, column, T(0));
			return m_values.back();
		}

		int k = search(row, column);
		if(k < m_rows[row+1] && m_columns[k] == column)
			return m_values[k];
/*
		To insert a new value...

		[ 0 1 2 ]    m_values = {1, 2, 3, 4}
		[ 0 0 0 ] -> m_columns = {1, 2, 0, 2}
		[ 3 0 4 ]    m_rows = {0, 2, 2, 4}

		[ 0 1 2 ]    m_values = {1, 2, x, 3, 4}
		[ 0 x 0 ] -> m_columns = {1, 2, 1, 0, 2}
		[ 3 0 4 ]    m_rows = {0, 2, 3, 5}
*/
		// shift elements to make room for the new one
		// if a zero is found, it can be overwritten
		int zeroIdx = k;
		while(zeroIdx < (int)m_values.size() && m_values[zeroIdx] != T(0))
			++zeroIdx;

		// if there is no zero to overwrite, append a zero
		if(zeroIdx == m_values.size())
		{
			m_values.push_back(T(0));
			m_columns.push_back(0);
		}
		for(int l = zeroIdx; l > k; --l)
		{
			m_values[l] = m_values[l-1];
			m_columns[l] = m_columns[l-1];
		}

		for(int l = row+1; l < (int)m_rows.size() && m_rows[l] <= zeroIdx; ++l)
			++m_rows[l];

		m_columns[k] = column;
		return (m_values[k] = T(0));
	}

	// 'in' and 'out' must not be the same object.
	friend void transpose(const SparseMatrix<T>& in, SparseMatrix<T>& out)
	{
		out.clear(in.m_width, in.m_height, (int)in.m_values.size());
		out.m_values.resize(in.m_values.size());
		out.m_columns.resize(in.m_values.size());
		out.m_rows.resize(out.m_height + 1, 0);

		// "counting sort"
		for(size_t k = 0; k < in.m_columns.size(); ++k)
			++out.m_rows[in.m_columns[k]];
		// add up frequencies to get indices to first element on each row.
		for(size_t k = 1; k < out.m_rows.size(); ++k)
			out.m_rows[k] += out.m_rows[k-1];
		// copy from 'in' to 'out'
		for(int i = in.m_rows.size() - 2; i >= 0; --i)
			for(int k = in.m_rows[i+1] - 1; k >= in.m_rows[i]; --k)
			{
				int idx = --out.m_rows[in.m_columns[k]];
				out.m_values[idx] = in.m_values[k];
				out.m_columns[idx] = i;
			}
		// remove empty rows
		while(out.m_rows.size() > 1 && (out.m_rows.back() == out.m_rows[out.m_rows.size() - 2]))
			out.m_rows.pop_back();
	}

	// Matrix-vector multiplication: out = m*v. 'v' and 'out' are treated as column vectors.
	// Sizes must match. The size of 'out' will be set.
	// 'v' and 'out' must not be the same object.
	// From "http://www.cs.utk.edu/~dongarra/etemplates/node382.html".
	friend void mul(const SparseMatrix<T>& m, const array<T>& v, array<T>& out)
	{
		assert(v.size() == m.m_width);
		out.clear();
		out.resize(m.m_height, T(0));
		for(int i = 0; i < (int)m.m_rows.size() - 1; ++i)
			for(int j = m.m_rows[i]; j < m.m_rows[i+1]; ++j)
				out[i] += m.m_values[j] * v[m.m_columns[j]];
	}

	// Vector-matrix multiplication: out = v*m. 'v' and 'out' are treated as row vectors.
	// Sizes must match. The size of 'out' will be set.
	// 'v' and 'out' must not be the same object.
	// From "http://www.cs.utk.edu/~dongarra/etemplates/node382.html".
	friend void mul(const array<T>& v, const SparseMatrix<T>& m, array<T>& out) 
	{
		assert(v.size() == m.m_height);
		out.clear();
		out.resize(m.m_width, T(0));
		for(int i = 0; i < m.m_rows.size() - 1; ++i)
			for(int j = m.m_rows[i]; j < m.m_rows[i+1]; ++j)
				out[m.m_columns[j]] += m.m_values[j] * v[i];
	}

	// Matrix-matrix multiplication: out = mat1*mat2.
	// Sizes must match. The size of 'out' will be set.
	// 'out' must not be the same object as 'mat1' or 'mat2'.
	friend void mul(const SparseMatrix<T>& mat1, const SparseMatrix<T>& mat2, SparseMatrix<T>& out)
	{
		assert(mat1.m_width == mat2.m_height);
		out.clear(mat1.m_height, mat2.m_width);
		for(int i = 0; i < (int)mat1.m_rows.size() - 1; ++i)
			for(int j = mat1.m_rows[i]; j < mat1.m_rows[i+1]; ++j)
				if(mat1.m_columns[j] < (int)mat2.m_rows.size() - 1)
					for(int k = mat2.m_rows[mat1.m_columns[j]]; k < mat2.m_rows[mat1.m_columns[j] + 1]; ++k)
						out(i, mat2.m_columns[k]) += mat2.m_values[k]*mat1.m_values[j];
	}

	// Matrix addition: out = mat1+mat2.
	// Sizes must match. The size of 'out' will be set.
	// 'out' must not be the same object as 'mat1' or 'mat2'.
	friend void add(const SparseMatrix<T>& mat1, const SparseMatrix<T>& mat2, SparseMatrix<T>& out)
	{
		assert(mat1.m_width == mat2.m_width && mat1.m_height == mat2.m_height);
		out.clear(mat1.m_height, mat2.m_width);
		int i = 0;
		while(i < (int)mat1.m_rows.size() - 1 && i < (int)mat2.m_rows.size() - 1)
		{
			int j1 = mat1.m_rows[i];
			int j2 = mat2.m_rows[i];
			while(j1 < mat1.m_rows[i+1] && j2 < mat2.m_rows[i+1])
			{
				while(j1 < mat1.m_rows[i+1] && mat1.m_columns[j1] < mat2.m_columns[j2])
				{
					out.m_values.push_back(mat1.m_values[j1]);
					out.m_columns.push_back(mat1.m_columns[j1]);
					++j1;
				}
				while(j2 < mat2.m_rows[i+1] && mat2.m_columns[j2] < mat1.m_columns[j1])
				{
					out.m_values.push_back(mat2.m_values[j2]);
					out.m_columns.push_back(mat2.m_columns[j2]);
					++j2;
				}
				while(j1 < mat1.m_rows[i+1] && j2 < mat2.m_rows[i+1] && mat1.m_columns[j1] == mat2.m_columns[j2])
				{
					out.m_values.push_back(mat1.m_values[j1] + mat2.m_values[j2]);
					out.m_columns.push_back(mat1.m_columns[j1]);
					++j1; ++j2;
				}
			}
			while(j1 < mat1.m_rows[i+1])
			{
				out.m_values.push_back(mat1.m_values[j1]);
				out.m_columns.push_back(mat1.m_columns[j1]);
				++j1;
			}
			while(j2 < mat2.m_rows[i+1])
			{
				out.m_values.push_back(mat2.m_values[j2]);
				out.m_columns.push_back(mat2.m_columns[j2]);
				++j2;
			}
			out.m_rows.push_back((int)out.m_values.size());
			++i;
		}
		while(i < (int)mat1.m_rows.size() - 1)
		{
			for(int j = mat1.m_rows[i]; j < mat1.m_rows[i+1]; ++j)
			{
				out.m_values.push_back(mat1.m_values[j]);
				out.m_columns.push_back(mat1.m_columns[j]);
			}
			out.m_rows.push_back((int)out.m_values.size());
			++i;
		}
		while(i < (int)mat2.m_rows.size() - 1)
		{
			for(int j = mat2.m_rows[i]; j < mat2.m_rows[i+1]; ++j)
			{
				out.m_values.push_back(mat2.m_values[j]);
				out.m_columns.push_back(mat2.m_columns[j]);
			}
			out.m_rows.push_back((int)out.m_values.size());
			++i;
		}
	}

	// Matrix subtraction: out = mat1-mat2.
	// Sizes must match. The size of 'out' will be set.
	// 'out' must not be the same object as 'mat1' or 'mat2'.
	friend void sub(const SparseMatrix<T>& mat1, const SparseMatrix<T>& mat2, SparseMatrix<T>& out)
	{
		assert(mat1.m_width == mat2.m_width && mat1.m_height == mat2.m_height);
		out.clear(mat1.m_height, mat2.m_width);
		int i = 0;
		while(i < (int)mat1.m_rows.size() - 1 && i < (int)mat2.m_rows.size() - 1)
		{
			int j1 = mat1.m_rows[i];
			int j2 = mat2.m_rows[i];
			while(j1 < mat1.m_rows[i+1] && j2 < mat2.m_rows[i+1])
			{
				while(j1 < mat1.m_rows[i+1] && mat1.m_columns[j1] < mat2.m_columns[j2])
				{
					out.m_values.push_back(mat1.m_values[j1]);
					out.m_columns.push_back(mat1.m_columns[j1]);
					++j1;
				}
				while(j2 < mat2.m_rows[i+1] && mat2.m_columns[j2] < mat1.m_columns[j1])
				{
					out.m_values.push_back(-mat2.m_values[j2]);
					out.m_columns.push_back(mat2.m_columns[j2]);
					++j2;
				}
				while(j1 < mat1.m_rows[i+1] && j2 < mat2.m_rows[i+1] && mat1.m_columns[j1] == mat2.m_columns[j2])
				{
					out.m_values.push_back(mat1.m_values[j1] - mat2.m_values[j2]);
					out.m_columns.push_back(mat1.m_columns[j1]);
					++j1; ++j2;
				}
			}
			while(j1 < mat1.m_rows[i+1])
			{
				out.m_values.push_back(mat1.m_values[j1]);
				out.m_columns.push_back(mat1.m_columns[j1]);
				++j1;
			}
			while(j2 < mat2.m_rows[i+1])
			{
				out.m_values.push_back(-mat2.m_values[j2]);
				out.m_columns.push_back(mat2.m_columns[j2]);
				++j2;
			}
			out.m_rows.push_back((int)out.m_values.size());
			++i;
		}
		while(i < (int)mat1.m_rows.size() - 1)
		{
			for(int j = mat1.m_rows[i]; j < mat1.m_rows[i+1]; ++j)
			{
				out.m_values.push_back(mat1.m_values[j]);
				out.m_columns.push_back(mat1.m_columns[j]);
			}
			out.m_rows.push_back((int)out.m_values.size());
			++i;
		}
		while(i < (int)mat2.m_rows.size() - 1)
		{
			for(int j = mat2.m_rows[i]; j < mat2.m_rows[i+1]; ++j)
			{
				out.m_values.push_back(-mat2.m_values[j]);
				out.m_columns.push_back(mat2.m_columns[j]);
			}
			out.m_rows.push_back((int)out.m_values.size());
			++i;
		}
	}

	void negate()
	{
		for(int i = 0; i < (int)m_values.size(); ++i)
			m_values[i] = -m_values[i];
	}

	// Kronecker multiplication: out = mat1*mat2.
	// The size of 'out' will be set.
	// 'out' must not be the same object as 'mat1' or 'mat2'.
	friend void kron(const SparseMatrix<T>& mat1, const SparseMatrix<T>& mat2, SparseMatrix<T>& out)
	{
		out.clear(mat1.m_height*mat2.m_height, mat1.m_width*mat2.m_width, (int)(mat1.m_values.size()*mat2.m_values.size()));

		// out.m_rows[0] is already zero
		for(int i = 0; i < (int)mat1.m_rows.size() - 1; ++i)
			for(int j = 0; j < mat2.m_height; ++j)
			{
				if(j < (int)mat2.m_rows.size() - 1)
					for(int k = mat1.m_rows[i]; k < mat1.m_rows[i+1]; ++k)
						for(int l = mat2.m_rows[j]; l < mat2.m_rows[j+1]; ++l)
						{
							out.m_values.push_back(mat1.m_values[k]*mat2.m_values[l]);
							out.m_columns.push_back(mat1.m_columns[k]*mat2.m_width + mat2.m_columns[l]);
						}
				out.m_rows.push_back((int)out.m_values.size());
			}
	}

	void loadIdentity(T diagonal = T(1), int size = 0)
	{
		clear(size, size, size);
		if(m_height != m_width) m_width = m_height;
		for(int i = 0; i < m_height; ++i)
		{
			m_rows.push_back(i+1);
			m_values.push_back(diagonal);
			m_columns.push_back(i);
		}
	}

	void loadTridiag(T sub, T diag, T super, int size = 0)
	{
		clear(size, size, size);
		if(m_height != m_width) m_width = m_height;
		m_values.push_back(diag);
		m_columns.push_back(0);
		for(int i = 1; i < m_height; ++i)
		{
			m_values.push_back(super);
			m_values.push_back(sub);
			m_values.push_back(diag);
			m_columns.push_back(i);
			m_columns.push_back(i-1);
			m_columns.push_back(i);
			m_rows.push_back(3 * i - 1);
		}
		m_rows.push_back(m_height * 3 - 2);
	}

	void loadZeros(int elementsPerRow, int height = 0, int width = 0)
	{
		clear(height, width);
		assert(elementsPerRow <= m_width);
		m_values.resize(elementsPerRow * m_height, T(0));
		m_columns.reserve(elementsPerRow * m_height);
		m_rows.reserve(m_height + 1);
		for(int i = 0; i < m_height; ++i)
		{
			for(int j = 0; j < elementsPerRow; ++j)
				m_columns.push_back(m_width - elementsPerRow + j);
			m_rows.push_back(elementsPerRow * (i+1));
		}
	}

	void loadLaplace(int dims, int* sizes, T* deltas, bool dirichlet = true)
	{
		array<int> index(dims, 0);
		array<int> offsets(dims);
		offsets.back() = 1;
		for(int i = dims - 1; i > 0; --i)
			offsets[i-1] = offsets[i] * sizes[i];
		int prod = offsets[0] * sizes[0];

		clear(prod, prod, prod * (2 * dims + 1));

		// for each line...
		int row = 0;
		while(index[0] < sizes[0])
		{
			for(int j = 0; j < dims; ++j)
				if(index[j] != 0)
				{
					m_values.push_back(T(-1)/(deltas[j]*deltas[j]));
					m_columns.push_back(row - offsets[j]);
				}
			m_values.push_back(T(0));
			m_columns.push_back(row);
			for(int j = 0; j < dims; ++j)
			{
				T val = T(0);
				if((index[j] != 0) || dirichlet) val += T(1);
				if((index[j] != sizes[j] - 1) || dirichlet) val += T(1);
				m_values.back() += val/(deltas[j]*deltas[j]);
			}
			for(int j = dims - 1; j >= 0; --j)
				if(index[j] != sizes[j] - 1)
				{
					m_values.push_back(T(-1)/(deltas[j]*deltas[j]));
					m_columns.push_back(row + offsets[j]);
				}
			m_rows.push_back((int)m_values.size());

			int idx = (int)index.size() - 1;
			++index[idx];
			while(idx > 0 && index[idx] >= sizes[idx])
			{
				index[idx--] = 0;
				++index[idx];
			}
			++row;
		}
	}

	bool operator == (const SparseMatrix<T>& mat) const
	{
		if(m_height != mat.m_height || m_width != mat.m_width)
			return false;

		int i = 0;
		while(i < (int)m_rows.size() - 1 && i < (int)mat.m_rows.size() - 1)
		{
			int j = m_rows[i];
			int mat_j = mat.m_rows[i];

			while(j < m_rows[i+1] && mat_j < mat.m_rows[i+1])
			{
				while(j < m_rows[i+1] && m_columns[j] < mat.m_columns[mat_j])
					if(m_values[j++] != T(0))
						return false;
				while(mat_j < mat.m_rows[i+1] && mat.m_columns[mat_j] < m_columns[j])
					if(mat.m_values[mat_j++] != T(0))
						return false;
				while(j < m_rows[i+1] && mat_j < mat.m_rows[i+1] && m_columns[j] == mat.m_columns[mat_j])
					if(m_values[j++] != mat.m_values[mat_j++])
						return false;
			}
			while(j < m_rows[i+1])
				if(m_values[j++] != T(0))
					return false;
			while(mat_j < mat.m_rows[i+1])
				if(mat.m_values[mat_j++] != T(0))
					return false;
			++i;
		}
		while(i < (int)m_rows.size() - 1)
		{
			int j = m_rows[i];
			while(j < m_rows[i+1])
				if(m_values[j++] != T(0))
					return false;
			++i;
		}
		while(i < (int)mat.m_rows.size() - 1)
		{
			int j = mat.m_rows[i];
			while(j < mat.m_rows[i+1])
				if(mat.m_values[j++] != T(0))
					return false;
			++i;
		}

		return true;
	}

	bool operator != (const SparseMatrix<T>& mat) const
	{
		return !(*this == mat);
	}

	friend std::ostream& operator << (std::ostream& os, const SparseMatrix<T>& mat)
	{
		os << "[";
		for(int i = 0; i < mat.m_height; ++i)
		{
//			if(i != 0) os << ", ";
			if(i != 0) os << std::endl << ' ';
			os << "[";
			int col = 0;
			if(i < (int)mat.m_rows.size() - 1)
			{
				for(int j = mat.m_rows[i]; j < mat.m_rows[i+1]; ++j)
				{
					while(col < mat.m_columns[j])
					{
						if(col != 0) os << ", ";
						os << T(0);
						++col;
					}
					if(col != 0) os << ", ";
					os << mat.m_values[j];
					++col;
				}
			}
			while(col < mat.m_width)
			{
				if(col != 0) os << ", ";
				os << T(0);
				++col;
			}
			os << "]";
		}
		os << "]";
		return os;
	}

private:
	// Returns the index of the value at position (i, j) in the matrix.
	// If not found, returns the index where the value should have been.
	// (O(log n) where n is the number of elements on one row)
	int search(int i, int j) const
	{
		assert(i >= 0 && i < m_height && j >= 0 && j < m_width);
		// binary search for element

		// search for 5:
		//   0 1 2 3[4]6 7 8 9
		// ^         ^         ^
		//   0 1 2 3 4 6 7[8]9
		//           ^     ^   ^
		//   0 1 2 3 4 6[7]8 9
		//           ^   ^ ^
		//   0 1 2 3 4[6]7 8 9
		//           ^ ^ ^
		//   0 1 2 3 4[6]7 8 9
		//           ^ ^

		// search for 0:
		//   0 1 2 3[4]6 7 8 9
		// ^         ^         ^
		//   0 1[2]3 4 6 7 8 9
		// ^     ^   ^
		//   0[1]2 3 4 6 7 8 9
		// ^   ^ ^
		//  [0]1 2 3 4 6 7 8 9
		// ^ ^ ^
		//  [0]1 2 3 4 6 7 8 9
		// ^ ^

		if(i > (int)m_rows.size() - 2)
			return (int)m_values.size();

		int k0 = m_rows[i] - 1;
		int k1 = m_rows[i+1];
		int k = (k1 + k0 + 1) >> 1;
		while(k1 - k0 > 1)
		{
			if(j <= m_columns[k])
				k1 = k;
			else
				k0 = k;
			k = (k1 + k0 + 1) >> 1;
		}
		return k;
	}

	// Compressed Row Storage Format as described on
	// "http://www.cs.utk.edu/~dongarra/etemplates/node373.html"
	// (Slightly modified: 'm_rows' can be shorter if the last rows are zero)
	array<T> m_values;	// all non-zero elements
	array<int> m_columns;	// column index
	array<int> m_rows;	// index of first element in a row
	int m_height, m_width;
};

// Sizes must match.
template<class T>
T dot(const array<T>& x, const array<T>& y)
{
	assert(x.size() == y.size());
	T sum = T(0);
	for(size_t i = 0; i < x.size(); ++i)
		sum += x[i]*y[i];
	return sum;
}

// Solves the equation Ax=b by using conjugate gradients.
// A must be symmetric and positive definite.
// Sizes must match. The size of 'x' will be set.
// From "http://www.math-linux.com/spip.php?article55"
template<class T>
void conjGrad(const SparseMatrix<T>& A, const array<T>& b, array<T>& x, T eps = 1e-14)
{
	assert(A.height() == b.size());
/*
	x.clear();
	x.resize(b.size(), T(0));
	array<T> rk, pk, Apk;
	rk = b;
	pk = rk;
	T alpha, beta, dot_rk, dot_rkp1;
	dot_rkp1 = dot(rk, rk);
	if(dot_rkp1 <= eps*eps)
		return;
	for(int k = 0; k < (int)b.size(); ++k)
	{
		dot_rk = dot_rkp1;
		mul(A, pk, Apk);
		alpha = dot_rk/dot(pk, Apk);
		for(int i = 0; i < (int)x.size(); ++i)
		{
			x[i] += alpha*pk[i];
			rk[i] -= alpha*Apk[i];
		}
		// test 'rk'
		dot_rkp1 = dot(rk, rk);
		if(dot_rkp1 <= eps*eps)
			return;
		beta = dot_rkp1/dot_rk;
		for(int i = 0; i < (int)pk.size(); ++i)
			pk[i] = rk[i] + beta*pk[i];
	}
/*/
	array<T> r, z(b.size()), Ad, invC(b.size());
	T alpha, beta, dot_r, dot_rp1;

	mul(A, x, r);
	for(int i = 0; i < (int)x.size(); ++i)
	{
		r[i] = b[i] - r[i]; // r_0 = b - A x_0
		invC[i] = T(1)/A(i, i);
		z[i] = invC[i] * r[i];
	}
	if(dot(r, r) <= eps*eps)
		return;

	array<T> d(z);
	dot_rp1 = dot(z, r);
	for(int k = 0; k < (int)x.size(); ++k)
	{
		dot_r = dot_rp1;
		mul(A, d, Ad);

		alpha = dot_r/dot(d, Ad);
		for(size_t i = 0; i < x.size(); ++i)
		{
			x[i] += alpha*d[i];
			r[i] -= alpha*Ad[i];
			z[i] = invC[i]*r[i];
		}
		if(dot(r, r) <= eps*eps)
			return;
		dot_rp1 = dot(z, r);
		beta = dot_rp1/dot_r;
		for(size_t i = 0; i < x.size(); ++i)
			d[i] = z[i] + beta*d[i];
	}

/**/
}

template<class T>
std::ostream& operator << (std::ostream& os, const array<T>& v)
{
	os << "[";
	if(v.size() != 0)
		os << v[0];
	for(size_t i = 1; i < v.size(); ++i)
		os << ", " << v[i];
	return os << "]";
}

#endif
