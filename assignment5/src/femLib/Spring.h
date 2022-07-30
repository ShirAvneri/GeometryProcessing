#pragma once

#include "Element.h"

/**
	This class implements the interface for an elementary energy unit. As a function of deformed, undeformed,
	and other parameters, such as boundary conditions, each class that extends this one will define a potential energy.
	The deformed energy depends on a number of nodes.
*/
class Spring : public Element
{
public:
	Spring(const std::array<int, 2>& nodeIndices, const VectorXd& X) : nodeIndices(nodeIndices) { }
	virtual ~Spring() { }

	// Returns the number of nodes this unit depends on
	virtual int getNumNodes() const { return 2; }

	// Returns the global index of node `i`
	virtual int getNodeIndex(int i) const { return nodeIndices[i]; }

	// Returns the element's mass
	virtual double getMass() const { return 0; }

	// Returns the energy value given deformed `x` and undeformed `X` state
	virtual double getEnergy(const VectorXd& x, const VectorXd& X)
	{
		double l = (getNodePos(1, x) - getNodePos(0, x)).norm(), L = (getNodePos(1, X) - getNodePos(0, X)).norm();
		return 0.5 * k * pow(l / L - 1, 2) * L;
	}

	// Adds the gradient to `grad` given deformed `x` and undeformed `X` state
	virtual void addEnergyGradientTo(const VectorXd& x, const VectorXd& X, VectorXd& grad)
	{
		Vector2d l_vec = getNodePos(1, x) - getNodePos(0, x);
		double l = l_vec.norm(), L = (getNodePos(1, X) - getNodePos(0, X)).norm();
		Vector2d f = -k * (l / L - 1) * l_vec / l;

		size_t i = static_cast<size_t>(2) * getNodeIndex(0);
		grad(i) += f(0);
		grad(i + 1) += f(1);
		i = static_cast<size_t>(2) * getNodeIndex(1);
		grad(i) += -f(0);
		grad(i + 1) += -f(1);
	}

	// Adds the hessian entries to `hesEntries` given deformed `x` and undeformed `X` state
	virtual void addEnergyHessianTo(const VectorXd& x, const VectorXd& X, std::vector<Tripletd>& hesEntries)
	{
		Vector2d l_vec = getNodePos(1, x) - getNodePos(0, x);
		double l = l_vec.norm(), L = (getNodePos(1, X) - getNodePos(0, X)).norm();

		Matrix2d temp = (l_vec * l_vec.transpose()) / (l_vec.transpose() * l_vec);
		Matrix2d H = -k * (1 / L * temp + (l / L - 1) / l * (Matrix2d::Identity() - temp));

		size_t i1 = static_cast<size_t>(2) * getNodeIndex(0), i2 = static_cast<size_t>(2) * getNodeIndex(1);

		hesEntries.push_back(Tripletd(i1, i1, -H(0, 0)));
		hesEntries.push_back(Tripletd(i1, i1 + 1, -H(0, 1)));
		hesEntries.push_back(Tripletd(i1 + 1, i1, -H(1, 0)));
		hesEntries.push_back(Tripletd(i1 + 1, i1 + 1, -H(1, 1)));

		hesEntries.push_back(Tripletd(i1, i2, H(0, 0)));
		hesEntries.push_back(Tripletd(i1, i2 + 1, H(0, 1)));
		hesEntries.push_back(Tripletd(i1 + 1, i2, H(1, 0)));
		hesEntries.push_back(Tripletd(i1 + 1, i2 + 1, H(1, 1)));

		hesEntries.push_back(Tripletd(i2, i1, H(0, 0)));
		hesEntries.push_back(Tripletd(i2, i1 + 1, H(0, 1)));
		hesEntries.push_back(Tripletd(i2 + 1, i1, H(1, 0)));
		hesEntries.push_back(Tripletd(i2 + 1, i1 + 1, H(1, 1)));

		hesEntries.push_back(Tripletd(i2, i2, -H(0, 0)));
		hesEntries.push_back(Tripletd(i2, i2 + 1, -H(0, 1)));
		hesEntries.push_back(Tripletd(i2 + 1, i2, -H(1, 0)));
		hesEntries.push_back(Tripletd(i2 + 1, i2 + 1, -H(1, 1)));
	}

protected:
	// the collection of nodes that define the triangle element
	std::array<int, 2> nodeIndices;
	// spring stiffness
	double k = 20.0;
};
