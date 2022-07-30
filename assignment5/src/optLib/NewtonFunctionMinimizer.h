#pragma once

#include "ObjectiveFunction.h"
#include "GradientDescentFunctionMinimizer.h"

/**
	use Newton's method to optimize a function. p will store the final value that minimizes the function, and its initial value
	is used to start the optimization method.

	Task: find p that minimize f(p). This means that df/dp(p) = 0.
	df/dp(p+dp) ~ df/dp(p) + d/dp(df/dp) * dp = 0 ==> -df/dp(p) = d/dp(df/dp) * dp
	Iterating the above, will hopefully get p that minimizes f.
*/

typedef Eigen::SparseMatrix<double> SparseXd;

class NewtonFunctionMinimizer : public GradientDescentFunctionMinimizer
{
public:
	NewtonFunctionMinimizer(int maxIterations = 100, double solveResidual = 0.0001, int maxLineSearchIterations = 15)
		: GradientDescentFunctionMinimizer(maxIterations, solveResidual, maxLineSearchIterations) {	}

	virtual ~NewtonFunctionMinimizer() { }

	SparseMatrixd H;
	std::vector<Tripletd> hessianEntries;

protected:
	// The search direction is given by -Hinv * g
	virtual void computeSearchDirection(ObjectiveFunction* function, const VectorXd& x, VectorXd& dx)
	{
		hessianEntries.clear();
		function->addHessianEntriesTo(hessianEntries, x);

		SparseXd hessian(x.rows(), x.rows());
		hessian.setFromTriplets(hessianEntries.begin(), hessianEntries.end());

		VectorXd gradient(x.rows());
		gradient.setZero();
		function->addGradientTo(gradient, x);

		Eigen::SimplicialLDLT<SparseXd> solver;
		solver.compute(hessian);
		dx = solver.solve(-gradient.sparseView());
	}
};
