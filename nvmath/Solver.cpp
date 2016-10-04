// This code is in the public domain -- castanyo@yahoo.es

#include "nvmath.h"


using namespace nv;

namespace nv {
namespace solver {

class JacobiPreconditioner
{
public:

	JacobiPreconditioner(const sparse::Matrix &M, bool symmetric) : m_inverseDiagonal(M.width())
	{
		nvCheck(M.isSquare());
		for (uint32_t x = 0; x < M.width(); x++) {
			float elem = M.getCoefficient(x, x);
			//nvDebugCheck( elem != 0.0f ); // This can be zero in the presence of zero area triangles.
			if (symmetric) {
				m_inverseDiagonal[x] = (elem != 0) ? 1.0f / sqrtf(fabsf(elem)) : 1.0f;
			} else {
				m_inverseDiagonal[x] = (elem != 0) ? 1.0f / elem : 1.0f;
			}
		}
	}

	void apply(const FullVector &x, FullVector &y) const
	{
		nvDebugCheck(x.dimension() == m_inverseDiagonal.dimension());
		nvDebugCheck(y.dimension() == m_inverseDiagonal.dimension());
		// @@ Wrap vector component-wise product into a separate function.
		const uint32_t D = x.dimension();
		for (uint32_t i = 0; i < D; i++) {
			y[i] = m_inverseDiagonal[i] * x[i];
		}
	}

private:

	FullVector m_inverseDiagonal;

};

/**
* Compute the solution of the sparse linear system Ab=x using the Conjugate
* Gradient method.
*
* Solving sparse linear systems:
* (1)		A·x = b
*
* The conjugate gradient algorithm solves (1) only in the case that A is
* symmetric and positive definite. It is based on the idea of minimizing the
* function
*
* (2)		f(x) = 1/2·x·A·x - b·x
*
* This function is minimized when its gradient
*
* (3)		df = A·x - b
*
* is zero, which is equivalent to (1). The minimization is carried out by
* generating a succession of search directions p.k and improved minimizers x.k.
* At each stage a quantity alfa.k is found that minimizes f(x.k + alfa.k·p.k),
* and x.k+1 is set equal to the new point x.k + alfa.k·p.k. The p.k and x.k are
* built up in such a way that x.k+1 is also the minimizer of f over the whole
* vector space of directions already taken, {p.1, p.2, . . . , p.k}. After N
* iterations you arrive at the minimizer over the entire vector space, i.e., the
* solution to (1).
*
* For a really good explanation of the method see:
*
* "An Introduction to the Conjugate Gradient Method Without the Agonizing Pain",
* Jonhathan Richard Shewchuk.
*
**/
/*static*/ bool ConjugateGradientSolver(const sparse::Matrix &A, const FullVector &b, FullVector &x, float epsilon)
{
	nvDebugCheck( A.isSquare() );
	nvDebugCheck( A.width() == b.dimension() );
	nvDebugCheck( A.width() == x.dimension() );
	int i = 0;
	const int D = A.width();
	const int i_max = 4 * D;   // Convergence should be linear, but in some cases, it's not.
	FullVector r(D);   // residual
	FullVector p(D);   // search direction
	FullVector q(D);   //
	float delta_0;
	float delta_old;
	float delta_new;
	float alpha;
	float beta;
	// r = b - A·x;
	sparse::copy(b, r);
	sparse::sgemv(-1, A, x, 1, r);
	// p = r;
	sparse::copy(r, p);
	delta_new = sparse::dot( r, r );
	delta_0 = delta_new;
	while (i < i_max && delta_new > epsilon * epsilon * delta_0) {
		i++;
		// q = A·p
		mult(A, p, q);
		// alpha = delta_new / p·q
		alpha = delta_new / sparse::dot( p, q );
		// x = alfa·p + x
		sparse::saxpy(alpha, p, x);
		if ((i & 31) == 0) { // recompute r after 32 steps
			// r = b - A·x
			sparse::copy(b, r);
			sparse::sgemv(-1, A, x, 1, r);
		} else {
			// r = r - alpha·q
			sparse::saxpy(-alpha, q, r);
		}
		delta_old = delta_new;
		delta_new = sparse::dot( r, r );
		beta = delta_new / delta_old;
		// p = beta·p + r
		sparse::scal(beta, p);
		sparse::saxpy(1, r, p);
	}
	return delta_new <= epsilon * epsilon * delta_0;
}


// Conjugate gradient with preconditioner.
/*static*/ bool ConjugateGradientSolver(const JacobiPreconditioner &preconditioner, const sparse::Matrix &A, const FullVector &b, FullVector &x, float epsilon)
{
	nvDebugCheck( A.isSquare() );
	nvDebugCheck( A.width() == b.dimension() );
	nvDebugCheck( A.width() == x.dimension() );
	int i = 0;
	const int D = A.width();
	const int i_max = 4 * D;   // Convergence should be linear, but in some cases, it's not.
	FullVector r(D);    // residual
	FullVector p(D);    // search direction
	FullVector q(D);    //
	FullVector s(D);    // preconditioned
	float delta_0;
	float delta_old;
	float delta_new;
	float alpha;
	float beta;
	// r = b - A·x
	sparse::copy(b, r);
	sparse::sgemv(-1, A, x, 1, r);
	// p = M^-1 · r
	preconditioner.apply(r, p);
	delta_new = sparse::dot(r, p);
	delta_0 = delta_new;
	while (i < i_max && delta_new > epsilon * epsilon * delta_0) {
		i++;
		// q = A·p
		mult(A, p, q);
		// alpha = delta_new / p·q
		alpha = delta_new / sparse::dot(p, q);
		// x = alfa·p + x
		sparse::saxpy(alpha, p, x);
		if ((i & 31) == 0) { // recompute r after 32 steps
			// r = b - A·x
			sparse::copy(b, r);
			sparse::sgemv(-1, A, x, 1, r);
		} else {
			// r = r - alfa·q
			sparse::saxpy(-alpha, q, r);
		}
		// s = M^-1 · r
		preconditioner.apply(r, s);
		delta_old = delta_new;
		delta_new = sparse::dot( r, s );
		beta = delta_new / delta_old;
		// p = s + beta·p
		sparse::scal(beta, p);
		sparse::saxpy(1, s, p);
	}
	return delta_new <= epsilon * epsilon * delta_0;
}

static bool SymmetricSolver(const sparse::Matrix &A, const FullVector &b, FullVector &x, float epsilon = 1e-5f)
{
	nvDebugCheck(A.height() == A.width());
	nvDebugCheck(A.height() == b.dimension());
	nvDebugCheck(b.dimension() == x.dimension());
	JacobiPreconditioner jacobi(A, true);
	return ConjugateGradientSolver(jacobi, A, b, x, epsilon);
}

// Solve the symmetric system: At·A·x = At·b
bool LeastSquaresSolver(const sparse::Matrix &A, const FullVector &b, FullVector &x, float epsilon/*1e-5f*/)
{
	nvDebugCheck(A.width() == x.dimension());
	nvDebugCheck(A.height() == b.dimension());
	nvDebugCheck(A.height() >= A.width()); // @@ If height == width we could solve it directly...
	const uint32_t D = A.width();
	sparse::Matrix At(A.height(), A.width());
	sparse::transpose(A, At);
	FullVector Atb(D);
	sparse::mult(At, b, Atb);
	sparse::Matrix AtA(D);
	sparse::mult(At, A, AtA);
	return SymmetricSolver(AtA, Atb, x, epsilon);
}


// See section 10.4.3 in: Mesh Parameterization: Theory and Practice, Siggraph Course Notes, August 2007
bool LeastSquaresSolver(const sparse::Matrix &A, const FullVector &b, FullVector &x, const uint32_t *lockedParameters, uint32_t lockedCount, float epsilon/*= 1e-5f*/)
{
	nvDebugCheck(A.width() == x.dimension());
	nvDebugCheck(A.height() == b.dimension());
	nvDebugCheck(A.height() >= A.width() - lockedCount);
	// @@ This is not the most efficient way of building a system with reduced degrees of freedom. It would be faster to do it on the fly.
	const uint32_t D = A.width() - lockedCount;
	nvDebugCheck(D > 0);
	// Compute: b - Al * xl
	FullVector b_Alxl(b);
	for (uint32_t y = 0; y < A.height(); y++) {
		const uint32_t count = A.getRow(y).size();
		for (uint32_t e = 0; e < count; e++) {
			uint32_t column = A.getRow(y)[e].x;
			bool isFree = true;
			for (uint32_t i = 0; i < lockedCount; i++) {
				isFree &= (lockedParameters[i] != column);
			}
			if (!isFree) {
				b_Alxl[y] -= x[column] * A.getRow(y)[e].v;
			}
		}
	}
	// Remove locked columns from A.
	sparse::Matrix Af(D, A.height());
	for (uint32_t y = 0; y < A.height(); y++) {
		const uint32_t count = A.getRow(y).size();
		for (uint32_t e = 0; e < count; e++) {
			uint32_t column = A.getRow(y)[e].x;
			uint32_t ix = column;
			bool isFree = true;
			for (uint32_t i = 0; i < lockedCount; i++) {
				isFree &= (lockedParameters[i] != column);
				if (column > lockedParameters[i]) ix--; // shift columns
			}
			if (isFree) {
				Af.setCoefficient(ix, y, A.getRow(y)[e].v);
			}
		}
	}
	// Remove elements from x
	FullVector xf(D);
	for (uint32_t i = 0, j = 0; i < A.width(); i++) {
		bool isFree = true;
		for (uint32_t l = 0; l < lockedCount; l++) {
			isFree &= (lockedParameters[l] != i);
		}
		if (isFree) {
			xf[j++] = x[i];
		}
	}
	// Solve reduced system.
	bool result = LeastSquaresSolver(Af, b_Alxl, xf, epsilon);
	// Copy results back to x.
	for (uint32_t i = 0, j = 0; i < A.width(); i++) {
		bool isFree = true;
		for (uint32_t l = 0; l < lockedCount; l++) {
			isFree &= (lockedParameters[l] != i);
		}
		if (isFree) {
			x[i] = xf[j++];
		}
	}
	return result;
}
} // namespace solver
} // namespace nv
