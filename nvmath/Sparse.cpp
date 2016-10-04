// This code is in the public domain -- Ignacio Castaño <castanyo@yahoo.es>

#include "nvmath.h"

namespace nv {
namespace sparse {

void saxpy(float a, const FullVector &x, FullVector &y)
{
	nvDebugCheck(x.dimension() == y.dimension());
	const uint32_t dim = x.dimension();
	for (uint32_t i = 0; i < dim; i++) {
		y[i] += a * x[i];
	}
}

void copy(const FullVector &x, FullVector &y)
{
	nvDebugCheck(x.dimension() == y.dimension());
	const uint32_t dim = x.dimension();
	for (uint32_t i = 0; i < dim; i++) {
		y[i] = x[i];
	}
}

void scal(float a, FullVector &x)
{
	const uint32_t dim = x.dimension();
	for (uint32_t i = 0; i < dim; i++) {
		x[i] *= a;
	}
}

float dot(const FullVector &x, const FullVector &y)
{
	nvDebugCheck(x.dimension() == y.dimension());
	const uint32_t dim = x.dimension();
	float sum = 0;
	for (uint32_t i = 0; i < dim; i++) {
		sum += x[i] * y[i];
	}
	return sum;
}

// y = M * x
void mult(const Matrix &M, const FullVector &x, FullVector &y)
{
	mult(NoTransposed, M, x, y);
}

void mult(Transpose TM, const Matrix &M, const FullVector &x, FullVector &y)
{
	const uint32_t w = M.width();
	const uint32_t h = M.height();
	if (TM == Transposed) {
		nvDebugCheck( h == x.dimension() );
		nvDebugCheck( w == y.dimension() );
		y.fill(0.0f);
		for (uint32_t i = 0; i < h; i++) {
			M.madRow(i, x[i], y);
		}
	} else {
		nvDebugCheck( w == x.dimension() );
		nvDebugCheck( h == y.dimension() );
		for (uint32_t i = 0; i < h; i++) {
			y[i] = M.dotRow(i, x);
		}
	}
}

// y = alpha*A*x + beta*y
void sgemv(float alpha, const Matrix &A, const FullVector &x, float beta, FullVector &y)
{
	sgemv(alpha, NoTransposed, A, x, beta, y);
}

void sgemv(float alpha, Transpose TA, const Matrix &A, const FullVector &x, float beta, FullVector &y)
{
	const uint32_t w = A.width();
	const uint32_t h = A.height();
	if (TA == Transposed) {
		nvDebugCheck( h == x.dimension() );
		nvDebugCheck( w == y.dimension() );
		for (uint32_t i = 0; i < h; i++) {
			A.madRow(i, alpha * x[i], y);
		}
	} else {
		nvDebugCheck( w == x.dimension() );
		nvDebugCheck( h == y.dimension() );
		for (uint32_t i = 0; i < h; i++) {
			y[i] = alpha * A.dotRow(i, x) + beta * y[i];
		}
	}
}


// dot y-row of A by x-column of B
static float dotRowColumn(int y, const Matrix &A, int x, const Matrix &B)
{
	const std::vector<Matrix::Coefficient> &row = A.getRow(y);
	const uint32_t count = row.size();
	float sum = 0.0f;
	for (uint32_t i = 0; i < count; i++) {
		const Matrix::Coefficient &c = row[i];
		sum += c.v * B.getCoefficient(x, c.x);
	}
	return sum;
}

// dot y-row of A by x-row of B
static float dotRowRow(int y, const Matrix &A, int x, const Matrix &B)
{
	const std::vector<Matrix::Coefficient> &row = A.getRow(y);
	const uint32_t count = row.size();
	float sum = 0.0f;
	for (uint32_t i = 0; i < count; i++) {
		const Matrix::Coefficient &c = row[i];
		sum += c.v * B.getCoefficient(c.x, x);
	}
	return sum;
}

// dot y-column of A by x-column of B
static float dotColumnColumn(int y, const Matrix &A, int x, const Matrix &B)
{
	nvDebugCheck(A.height() == B.height());
	const uint32_t h = A.height();
	float sum = 0.0f;
	for (uint32_t i = 0; i < h; i++) {
		sum += A.getCoefficient(y, i) * B.getCoefficient(x, i);
	}
	return sum;
}


void transpose(const Matrix &A, Matrix &B)
{
	nvDebugCheck(A.width() == B.height());
	nvDebugCheck(B.width() == A.height());
	const uint32_t w = A.width();
	for (uint32_t x = 0; x < w; x++) {
		B.clearRow(x);
	}
	const uint32_t h = A.height();
	for (uint32_t y = 0; y < h; y++) {
		const std::vector<Matrix::Coefficient> &row = A.getRow(y);
		const uint32_t count = row.size();
		for (uint32_t i = 0; i < count; i++) {
			const Matrix::Coefficient &c = row[i];
			nvDebugCheck(c.x < w);
			B.setCoefficient(y, c.x, c.v);
		}
	}
}

// C = A * B
void mult(const Matrix &A, const Matrix &B, Matrix &C)
{
	mult(NoTransposed, A, NoTransposed, B, C);
}

void mult(Transpose TA, const Matrix &A, Transpose TB, const Matrix &B, Matrix &C)
{
	sgemm(1.0f, TA, A, TB, B, 0.0f, C);
}

// C = alpha*A*B + beta*C
void sgemm(float alpha, const Matrix &A, const Matrix &B, float beta, Matrix &C)
{
	sgemm(alpha, NoTransposed, A, NoTransposed, B, beta, C);
}

void sgemm(float alpha, Transpose TA, const Matrix &A, Transpose TB, const Matrix &B, float beta, Matrix &C)
{
	const uint32_t w = C.width();
	const uint32_t h = C.height();
	uint32_t aw = (TA == NoTransposed) ? A.width() : A.height();
	uint32_t ah = (TA == NoTransposed) ? A.height() : A.width();
	uint32_t bw = (TB == NoTransposed) ? B.width() : B.height();
	uint32_t bh = (TB == NoTransposed) ? B.height() : B.width();
	nvDebugCheck(aw == bh);
	nvDebugCheck(bw == ah);
	nvDebugCheck(w == bw);
	nvDebugCheck(h == ah);
	for (uint32_t y = 0; y < h; y++) {
		for (uint32_t x = 0; x < w; x++) {
			float c = beta * C.getCoefficient(x, y);
			if (TA == NoTransposed && TB == NoTransposed) {
				// dot y-row of A by x-column of B.
				c += alpha * dotRowColumn(y, A, x, B);
			} else if (TA == Transposed && TB == Transposed) {
				// dot y-column of A by x-row of B.
				c += alpha * dotRowColumn(x, B, y, A);
			} else if (TA == Transposed && TB == NoTransposed) {
				// dot y-column of A by x-column of B.
				c += alpha * dotColumnColumn(y, A, x, B);
			} else if (TA == NoTransposed && TB == Transposed) {
				// dot y-row of A by x-row of B.
				c += alpha * dotRowRow(y, A, x, B);
			}
			C.setCoefficient(x, y, c);
		}
	}
}

} // namespace sparse
} // namespace nv

