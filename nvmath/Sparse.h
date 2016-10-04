// This code is in the public domain -- castanyo@yahoo.es

#pragma once
#ifndef NV_MATH_SPARSE_H
#define NV_MATH_SPARSE_H

#include <vector>
#include "nvmath.h"


// Full and sparse vector and matrix classes. BLAS subset.

namespace nv
{
class SparseMatrix;

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
		uint32_t x;  // column
		float v; // value
	};


public:

	SparseMatrix(uint32_t d);
	SparseMatrix(uint32_t w, uint32_t h);
	SparseMatrix(const SparseMatrix &m);

	const SparseMatrix &operator=(const SparseMatrix &m);


	uint32_t width() const
	{
		return m_width;
	}
	uint32_t height() const
	{
		return m_array.size();
	}
	bool isSquare() const
	{
		return width() == height();
	}

	float getCoefficient(uint32_t x, uint32_t y) const; // x is column, y is row

	void setCoefficient(uint32_t x, uint32_t y, float f);
	void addCoefficient(uint32_t x, uint32_t y, float f);
	void mulCoefficient(uint32_t x, uint32_t y, float f);

	float sumRow(uint32_t y) const;
	float dotRow(uint32_t y, const FullVector &v) const;
	void madRow(uint32_t y, float alpha, FullVector &v) const;

	void clearRow(uint32_t y);
	void scaleRow(uint32_t y, float f);
	void normalizeRow(uint32_t y);

	void clearColumn(uint32_t x);
	void scaleColumn(uint32_t x, float f);

	const std::vector<Coefficient> &getRow(uint32_t y) const;

	bool isSymmetric() const;

private:

	/// Number of columns.
	const uint32_t m_width;

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
