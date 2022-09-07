//Sydney Loats SML4267

#include "beam.hpp"
#include <math.h>
#include <iostream>
#include "matrix.hpp"
#include "vector.hpp"
#include "timer.hpp"

using namespace std;

Beam:: Beam(int n, double L, double EI, double q){
	//cout << "entering beam class" << endl;
	n_ = n;
	L_ = L;
	EI_ = EI;
	q_ = q;
}

Matrix Beam:: getStiffnessMatrix(){
	Matrix ret;
	ret.initSize(n_, n_);
	
	for(int r=0; r<n_; r++){
		for(int c=0; c<n_; c++){
			ret.setVal(r, c, 0);
		}
	}

	ret.setVal(0, 0, 1);
	ret.setVal(n_-1, n_-1, 1);

	for(int i=1; i< n_-1; i++){
		ret.setVal(i, i, -2);
		ret.setVal(i, i-1, 1);
		ret.setVal(i, i+1, 1);
	}

	K_ = ret;

	return K_;
}

Vector Beam:: getSystemRHS(){
	Vector ret;
	ret.allocateData(n_);
	ret.setVal(0, 0);
	ret.setVal(n_-1, 0);
	double step = L_/(n_-1);
	double xval = step;
	for(int i=1; i<n_-1; i++){
		double val = (1/EI_) * (-q_/2) * (pow(xval, 2) - (xval*L_)) * pow(step, 2);
		ret.setVal(i, val);
		xval += step;
	}
	f_ = ret;
	return f_;
}

Vector Beam:: getExactSoln(){
	Vector ret;
	ret.allocateData(n_);
	
	ret.setVal(0, 0);
	ret.setVal(n_-1, 0);
	double step = L_/(n_-1);
	double xval = step;
	for(int i=1; i<n_-1; i++){
		ret.setVal(i, compute(xval));
		xval += step;
	}

	return ret;
}

double Beam:: compute(double x){
	double step = L_/(n_-1);
	double term1 = (-q_*pow(x, 4))/24;
	double term2 = (q_*L_*pow(x, 3))/12;
	double term3 = ((-q_*pow(L_, 3)*x)/24);
	double ret = (1/EI_) * (term1 + term2 + term3);
	return ret;
}

Vector Beam:: getApproxSoln(){
	if(turbo_) {K_.setMatrixTurbo(true);}
	Vector ret = K_.Solve(f_, GAUSS_SEIDEL);
	iters_ = K_.getSolveIters();
	time_ = K_.getSolveTime();
	return ret;
}


double Beam:: getStepSize(){
	return L_/(n_-1);
}


Vector Beam:: getCoordValues(){
	Vector ret;
	ret.allocateData(n_);
	ret.setVal(0, 0);
	ret.setVal(n_-1, L_);
	double step = L_/(n_-1);
	double val = step;
	for(int i=1; i<n_-1; i++){
		ret.setVal(i, val);
		val += step;	
	}
	return ret;
}


double Beam:: l2norm(Vector exact, Vector approx){
	Vector sub = exact.subtractFrom(approx);
	double val = 0;
	for(int i=0; i<sub.numElems(); i++){
		val += pow(sub.getVal(i), 2);
	}
	val = pow(val, 0.5);
	double denom = 0;
	for(int j=0; j<exact.numElems(); j++){
		denom += pow(exact.getVal(j), 2);
	}
	denom = pow(denom, 0.5);
	if(denom == 0){
		cout << "error: the l2 norm is diving by zero, please use valid values" << endl;
		exit(1);
	}

	return val/denom;
}

int Beam:: getSolveIters(){
	return iters_;
}

double Beam:: getSolveTime(){
	return time_;
}

void Beam:: setTurbo(bool mode){
	turbo_ = mode;
}
