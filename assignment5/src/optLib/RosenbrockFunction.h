#pragma once

#include "ObjectiveFunction.h"

class RosenbrockFunction : public ObjectiveFunction
{
public:
	double a, b;

	RosenbrockFunction() : a(1), b(100) { }

	virtual double computeValue(const VectorXd& x)
	{
		return pow(a - x[0], 2) + b * pow(x[1] - pow(x[0], 2), 2);
	}

	virtual void addGradientTo(VectorXd& grad, const VectorXd& x)
	{
		auto temp = pow(x[0], 2) - x[1];
		grad[0] += -2 * a + 2 * x[0] * (1 + 2 * b * temp);
		grad[1] += -2 * b * temp;
	}

	virtual void addHessianEntriesTo(std::vector<Tripletd>& hessianEntries, const VectorXd& x)
	{
		auto temp = -4 * b * x[0];
		hessianEntries.push_back(Tripletd(0, 0, -4 * b * x[1] - 3 * x[0] * temp + 2));
		hessianEntries.push_back(Tripletd(0, 1, temp));
		hessianEntries.push_back(Tripletd(1, 0, temp));
		hessianEntries.push_back(Tripletd(1, 1, 2 * b));
	}
};
