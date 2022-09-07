//Sydney Loats SML4267

#include "matrix.hpp"
#include "vector.hpp"
#pragma once

using namespace std;

class Beam{

	private:
		Matrix K_;            //system stiffness matrix
		Vector f_;            //right hand forcing vector
		int n_;
		double L_;
		double EI_;
		double q_;
		bool turbo_ = false;
		int iters_;
		double time_;

	public:
		Beam(int n, double L, double EI, double q); // constructor
		Matrix getStiffnessMatrix();                // return system stiffness matrix
		Vector getSystemRHS();                      // return system forcing vector
		Vector getExactSoln();                      // return exact solution
		Vector getApproxSoln();                     // return finit-difference solution
		Vector getCoordValues();                    // retun discretized x-coord values
		double l2norm(Vector exact, Vector approx); // return l2 error norm between exact and approx. solutions

		double compute(double x);
		
		int getSolveIters();
		double getSolveTime();
		void setTurbo(bool mode);
		double getStepSize();

};
