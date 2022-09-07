//Sydney Loats SML4267

#include "matrix.hpp"
#include "vector.hpp"
#include "timer.hpp"
#include "beam.hpp"
#include <vector>

using namespace std;

int main(int argc, char *argv[]){

	if(argc < 5){
		cout << "error: wrong number of inputs" << endl;
		cout << endl;
		cout << "usage: ./beamer [n] [L] [EI] [q] <turbo>" << endl;
		cout << endl;
		cout << "where:" << endl;
		cout << "[n]          number of points" << endl;
		cout << "[L]          length of the beam" << endl;
		cout << "[EI]         (Young's modulus)*(2nd moment of intertia)" << endl;
		cout << "[q]          distributed load per unit length" << endl;
		cout << "<turbo>      optional argument (1=turbo mode)" << endl;
		exit(1);
	}

	int n = stoi(argv[1]);
	if(n<3) {
		cout << "error: your value of n is invald, please enter a number 3 or greater" << endl;
		exit(1);
	}
	double L = stod(argv[2]);
	double EI = stod(argv[3]);
	if(EI == 0){
		cout << "error: EI must be an non-zero value" << endl;
		exit(1);
	}
	double q = stod(argv[4]);
	if(q == 0){
		cout << "error: q must be an non-zero value" << endl;
		exit(1);
	}

	Beam b(n, L, EI, q);

	if(argc == 6){ 
		int t = stoi(argv[5]);
		if(t==1) {b.setTurbo(true);}
		else if(t==0) {b.setTurbo(false);}
		else{
			cout << "error: only values of 1 and 0 are valid for turbo mode" << endl;
			exit(1);
		}
	}

	Vector coord = b.getCoordValues();
	Matrix m = b.getStiffnessMatrix();
	Vector rhs = b.getSystemRHS();
	Vector exact = b.getExactSoln();
	Vector approx = b.getApproxSoln();

	cout << "Beamer Solver for Analytic Solution of a Beam" << endl;
	cout << "---------------------------------------------" << endl;

	cout << "n: " << n << endl;
	cout << "L (m): " << L << endl;
	cout << "EI (Nm^2): " << EI << endl;
	cout << "q (N/m): " << q << endl;
	cout << endl;
	cout << "Number of iterations required to find the solution: " << b.getSolveIters() << endl;
	cout << "Elapsed wallclock time (s): " << b.getSolveTime() << endl;
	cout << "Constant mesh interval, h (step size): " << b.getStepSize() << endl;
	cout << "l2 error: " << b.l2norm(exact, approx) << endl;


	return 0;
}
