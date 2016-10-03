// This code is in the public domain -- castanyo@yahoo.es

#pragma once
#ifndef NV_MATH_SPARSE_H
#define NV_MATH_SPARSE_H

#include <vector>
#include "nvmath.h"


// Full and sparse vector and matrix classes. BLAS subset.

namespace nv
{
class FullVector;
class SparseMatrix;


/// Fixed size vector class.
class FullVector
{
public:

	FullVector(uint dim);
	FullVector(const FullVector &v);

	const FullVector &operator=(const FullVector &v);

	uint dimension() const
	{
		return m_array.size();
	}

	const float &operator[]( uint index ) const
	{
		return m_array[index];
	}
	float &operator[] ( uint index )
	{
		return m_array[index];
	}

	void fill(float f);

	void operator+= (const FullVector &v);
	void operator-= (const FullVector &v);
	void operator*= (const FullVector &v);

	void operator+= (float f);
	void operator-= (float f);
	void operator*= (float f);


private:

	std::vector<float> m_array;

};

// Pseudo-BLAS interface.
void saxpy(float a, const FullVector &x, FullVector &y);   // y = a * x + y
void copy(const FullVector &x, FullVector &y);
void scal(float a, FullVector &x);
float dot(const FullVector &x, const FullVector &y);


enum Transpose {
	NoTransposed = 0,
	Transposed = 1
};

/**
* Sparse matrix class. The matrix is assumed to be sparse and to have
* very few non-zero elements, for this reason it's stored in indexed
* format. To multiply column vectors efficiently, the matrix stores
* the elements in indexed-column order, there is a list of indexed
* elements for each row of the matrix. As with the FullVector the
* dimension of the matrix is constant.
**/
class SparseMatrix
{
public:

	// An element of the sparse array.
	struct Coefficient {
		uint x;  // column
		float v; // value
	};


public:

	SparseMatrix(uint d);
	SparseMatrix(uint w, uint h);
	SparseMatrix(const SparseMatrix &m);

	const SparseMatrix &operator=(const SparseMatrix &m);


	uint width() const
	{
		return m_width;
	}
	uint height() const
	{
		return m_array.size();
	}
	bool isSquare() const
	{
		return width() == height();
	}

	float getCoefficient(uint x, uint y) const; // x is column, y is row

	void setCoefficient(uint x, uint y, float f);
	void addCoefficient(uint x, uint y, float f);
	void mulCoefficient(uint x, uint y, float f);

	float sumRow(uint y) const;
	float dotRow(uint y, const FullVector &v) const;
	void madRow(uint y, float alpha, FullVector &v) const;

	void clearRow(uint y);
	void scaleRow(uint y, float f);
	void normalizeRow(uint y);

	void clearColumn(uint x);
	void scaleColumn(uint x, float f);

	const std::vector<Coefficient> &getRow(uint y) const;

	bool isSymmetric() const;

private:

	/// Number of columns.
	const uint m_width;

	/// Array of matrix elements.
	std::vector< std::vector<Coefficient> > m_array;

};

void transpose(const SparseMatrix &A, SparseMatrix &B);

void mult(const SparseMatrix &M, const FullVector &x, FullVector &y);
void mult(Transpose TM, const SparseMatrix &M, const FullVector &x, FullVector &y);

// y = alpha*A*x + beta*y
void sgemv(float alpha, const SparseMatrix &A, const FullVector &x, float beta, FullVector &y);
void sgemv(float alpha, Transpose TA, const SparseMatrix &A, const FullVector &x, float beta, FullVector &y);

void mult(const SparseMatrix &A, const SparseMatrix &B, SparseMatrix &C);
void mult(Transpose TA, const SparseMatrix &A, Transpose TB, const SparseMatrix &B, SparseMatrix &C);

// C = alpha*A*B + beta*C
void sgemm(float alpha, const SparseMatrix &A, const SparseMatrix &B, float beta, SparseMatrix &C);
void sgemm(float alpha, Transpose TA, const SparseMatrix &A, Transpose TB, const SparseMatrix &B, float beta, SparseMatrix &C);

// C = At * A
void sqm(const SparseMatrix &A, SparseMatrix &C);

} // nv namespace


#endif // NV_MATH_SPARSE_H
