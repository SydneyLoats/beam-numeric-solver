
//Sydney Loats SML4267

#include <iostream>
#include <math.h>
#include <fstream>
#include "matrix.hpp"
#include "timer.hpp"
#include "vector.hpp"

using namespace std;

Matrix:: Matrix(){
//	cout << "constructor: entering the program" << endl;
}

Matrix::~Matrix(){
//	cout << "destructor: leaving the program" << endl;
}

//initialization error
void Matrix:: checkInit(){
	if(numRows_ < 0) {
		cout << "error: matrix has not been initialized properly" << endl;
		exit(1);
	}
	else if(numCols_ < 0){
		cout << "error: matrix has not been initialized properly" << endl;
		exit(1);
	}
}

//initializes identity matrix to size nxn
void Matrix:: initIdentity(int n){

	numRows_ = n;
	numCols_ = n;
	initSize(numRows_, numCols_);
	for(int r=0; r<numRows_; r++){
		for(int c=0; c<numCols_; c++){
			if(r==c){
				matrix_[r][c] = 1;
			}
		}
	}

}

//initializes matrix by reading in from a file
void Matrix:: initFromFile(string fileName){

	ifstream fin(fileName);
        int r, c;
	double d;
	fin >> r >> c;
	numRows_ = r;
	numCols_ = c;
	initSize(numRows_, numCols_);
	for(int i=0; i<numRows_; i++){
		for(int j=0; j<numCols_; j++){
			fin >> d;
			matrix_[i][j] = d;
		}
	}

}

//initializing matrix to size rxc
void Matrix:: initSize(int r, int c){
	matrix_.resize(r);
	for (int j = 0; j < r; j++){
		matrix_[j].resize(c);
        }
	numRows_ = r;
	numCols_ = c;
}

//checks for if matrix is square
bool Matrix:: isSquare(){
	checkInit();
	if(numRows_ == numCols_){ return true;}
	else {return false;}
}

//returns number of rows
int Matrix:: numRows(){
	checkInit();
	return numRows_;
}

//returns number of columns
int Matrix:: numCols(){
	checkInit();
	return numCols_;
}

//sets number of rows to r
void Matrix:: setRows(int r){
	numRows_ = r;
}

//sets number of columns to c
void Matrix:: setCols(int c){
	numCols_ = c;
}

//returns value at row,col location
double Matrix:: getVal(int row, int col){
	checkInit();

	if(row>numRows_-1){
		cout << "error: not a valid row/col" << endl;
		exit(1);
	}
	else if(col>numCols_-1){
		cout << "error: not a valid row/col" << endl;
		exit(1);
	}
	else{
		return matrix_[row][col];
	}
}

//sets value to val at row,col location
void Matrix:: setVal(int row, int col, double val){
	checkInit();

	if(row>numRows_-1){
		cout << "error: not a valid row/col" << endl;
		exit(1);
	}
	else if(col>numCols_-1){
		cout << "error: not a valid row/col" << endl;
		exit(1);
	}
	else{
		matrix_[row][col] = val;
	}
}

//returns matrix multiplied by Matrix B
Matrix Matrix:: Multiply(Matrix B){
	checkInit();

	Matrix multiply;
	int rA = numRows_; 
	int cA = numCols_;
	int rB = B.numRows();
	int cB = B.numCols();

	if(!(cA==rB)){
		cout << "error: matrices are not compatible sizes" << endl;
		exit(1);
	}

	else{
		multiply.initSize(rA, cB);
		multiply.setRows(rA);
		multiply.setCols(cB);

		for(int r=0; r<rA; r++){
			for(int c=0; c<cB; c++){
				double helpVal = helpMultiply(r, c, B);
				multiply.setVal(r, c, helpVal);
			}
		}
	}

	return multiply;

}

//helper for the Multiply(Matrix B) function
double Matrix:: helpMultiply(int r, int c, Matrix B){
	checkInit();

	double ret;

	if(!(numCols()==B.numRows())){
		cout << "error: matrices are not a compatible sizes" << endl;
		exit(1);
	}

	else{

		for(int i=0; i<numCols(); i++){
			double num1 = getVal(r, i);
			double num2 = B.getVal(i, c);
			double mul = num1*num2;
			ret += mul;
		}
		
	}
	return ret;

}

//returns matrix multiplied by scalar A
Matrix Matrix:: Multiply(double A){
	checkInit();

	Matrix scalar;
	scalar.initSize(numRows_, numCols_);
	scalar.setRows(numRows_);
	scalar.setCols(numCols_);
	for(int r=0; r<numRows_; r++){
		for(int c=0; c<numCols_; c++){
			scalar.setVal(r, c, A*matrix_[r][c]);
		}
	}
	return scalar;

}

//returns transposed matrix
Matrix Matrix:: Transpose(){
	checkInit();

	Matrix tran;
	tran.initSize(numCols_, numRows_);
	tran.setRows(numCols_);
	tran.setCols(numRows_);
	for(int r=0; r<tran.numRows(); r++){
		for(int c=0; c<tran.numCols(); c++){
			tran.setVal(r, c, matrix_[c][r]);
		}
	}
	return tran;

}

//returns vector of diagonal matrix values
vector<double> Matrix:: Diagonal(){
	checkInit();

	vector<double> diag;
	for(int r=0; r<numRows_; r++){
		for(int c=0; c<numCols_; c++){
			if(r==c){
				diag.push_back(matrix_[r][c]);
			}
		}
	}

	return diag;

}

//prints matrix
void Matrix:: Print(){
	checkInit();

	for(int i=0; i<numRows_; i++){
		cout << "|   " ;
		for(int j=0; j<numCols_; j++){
			cout << matrix_[i][j] << "   ";
		}
		cout << "|" << endl;
	}
}

//prints matrix with matrix name
void Matrix:: Print(string name){
	checkInit();
	
	int length = name.size();
	cout << name << " = " << endl;

	Print();
}

Vector Matrix:: createVector(int row){
	Vector ret;
	ret.allocateData(numCols_);
	for(int c=0; c<numCols_; c++){
		double val = matrix_[row][c];
		ret.setVal(c, val);
	}
	return ret;
}

void Matrix:: replaceRow(int r, Vector b){
	for(int i=0; i< b.numElems(); i++){
		matrix_[r][i] = b.getVal(i);
	}
}

void Matrix:: addColumn(Vector b){
	numCols_ += 1;
	for (int r = 0; r < numRows_; r++){
		matrix_[r].resize(numCols_);
		matrix_[r][numCols_-1] = b.getVal(r);
	}

}

void Matrix:: setMatrix(vector< vector<double>> m){
	matrix_ = m;
}


	
Vector Matrix:: Solve(Vector b, int mode){
	if(mode == GAUSS_ELIM) {

		Timer timer;
		timer.Start();

		Matrix reduced = reducedMatrix(b);
		Vector ret = reduced.solveMatrix();

		timer.Stop();
		solveTime_ = timer.ElapsedTime();

		return ret;
	}
	else if(mode == JACOBI){

		Timer timer;
		timer.Start();

		Vector recent;

		recent.allocateData(numRows_);

		iters_ = 1;

		for(int i=0; i<numRows_; i++){
			recent.setVal(i, 0);
		}
		
		Vector ret;
		ret.allocateData(numRows_);

		ret.setPrevious(recent.getVector());

		for(int iter=0; iter<maxIters_; iter++){
			ret = solveJacobi(b, recent);
			double norm = ret.l2norm();
			if(getDebugMode()){
				cout << "--> Iteration:    " << iters_ << "  norm = " << norm << endl; 
			}
			recent = ret;
			ret.setPrevious(recent.getVector());
			if(norm < tolerance_) {
				if(getDebugMode()){cout << "Converged!" << endl;}
				break;
			}
			iters_ ++;
		}	

		timer.Stop();
		solveTime_ = timer.ElapsedTime();

		return ret;
	}
	else if(mode == GAUSS_SEIDEL){
		
		Timer timer;
		timer.Start();

		Vector recent;
		recent.allocateData(numRows_);
		
		iters_ = 1;

		for(int i=0; i<numRows_; i++){
			recent.setVal(i, 0);
		}

		Vector ret;
		ret.allocateData(numRows_);
		ret.setPrevious(recent.getVector());

		for(int iter=0; iter<maxIters_; iter++){
			ret = solveGaussSeidel(b, recent);
			double norm = ret.l2norm();
			if(getDebugMode()){
				cout << "--> Iteration:    " << iters_ << "  norm = " << norm << endl; 
			}
			recent = ret;
			ret.setPrevious(recent.getVector());
			if(norm< tolerance_) {
				if(getDebugMode()){cout << "Converged!" << endl;}
				break;

			}
			iters_ ++;
		}	

		timer.Stop();
		solveTime_ = timer.ElapsedTime();

		return ret;
	
	}
	else{
		cout << "error: not a valid mode" << endl;
		exit(1);
	}

}

Matrix Matrix:: separateMatrix(){
	Matrix ret;
	ret.initSize(numRows_, numCols_-1);
	ret.setRows(numRows_);
	ret.setCols(numCols_-1);
	for(int r=0; r< numRows_; r++){
		for(int c=0; c< numCols_-1; c++){
			ret.setVal(r, c, matrix_[r][c]);
		}
	}

	return ret;
}

int Matrix:: getMaxIters(){
	return maxIters_;
}

Vector Matrix:: separateVector(){
	Vector ret;
	ret.allocateData(numRows_);
	for(int i=0; i< numRows_; i++){
		ret.setVal(i, matrix_[i][numCols_-1]);
	}
	return ret;
}

Matrix Matrix:: reducedMatrix(Vector b){

	int rowIters = numRows_ - 1;
	int colIters = numCols_ - 1;
	
	Matrix gauss;
	gauss.initSize(numRows_, numCols_);
	gauss.setMatrix(matrix_);
	gauss.setRows(numRows_);
	gauss.setCols(numCols_);
	gauss.addColumn(b);

	int initialIter = 1;
	for(int c=0; c<colIters; c++){
		for(int r=initialIter; r<rowIters+1; r++){
			Vector eqI = gauss.createVector(c);
			Vector eqX = gauss.createVector(r);
			if(gauss.getVal(c, c) == 0) {
				cout << "error: diving by zero, please flip one of the rows" << endl;
				exit(1);
			}
			double multiplier = gauss.getVal(r, c)/gauss.getVal(c, c);
			Vector eqITemp = eqI.multiply(multiplier);
			Vector nEqX = eqITemp.subtractFrom(eqX); 
			gauss.replaceRow(r, nEqX);
		}
		initialIter++;
	}	

	return gauss;

	
}

Vector Matrix:: solveMatrix(){
	
	Vector ret;
	ret.allocateData(numRows_);
	for(int i=0; i<numRows_; i++){
		ret.setVal(i, 0);
	}
	for(int r = numRows_-1; r>=0; r--){
		double val = matrix_[r][numCols_-1];
		for(int c=0; c<numCols_-1; c++){
			if(r==c) {c++;}
			double temp = ret.getVal(c)*matrix_[r][c];
			val -= ret.getVal(c) * matrix_[r][c];
		}
		double set = val/matrix_[r][r];
		ret.setVal(r, set);
	}
	return ret;
}

Vector Matrix:: solveJacobi(Vector b, Vector recent){
	Vector ret;
	ret.allocateData(numRows_);

	for(int r=0; r<numRows_; r++){

	       double retVal = b.getVal(r);

	       for(int c=0; c<numCols_; c++){
		        if(r==c) {c++;}
			retVal -= matrix_[r][c]*recent.getVal(c);
	       }
	       if(matrix_[r][r] == 0){
			cout << "error: dividing by zero, cannot be solved" << endl;
			exit(1);
	       }
	       ret.setVal(r,retVal/matrix_[r][r]);
	}

	return ret;
}

void Matrix:: setMatrixTurbo(bool t){
	turbo_ = t;
}

Vector Matrix:: solveGaussSeidel(Vector b, Vector recent){
	Vector ret;
	ret.allocateData(numRows_);

	for(int i=0; i<numRows_; i++){
		ret.setVal(i, 0);
	}
	
	for(int r=1; r<numRows_-1; r++){
		
		double retVal = b.getVal(r);
		
		if(turbo_){

			for(int c=0; c<2; c++){
				if(c==0){
					retVal -= matrix_[r][r-1]*ret.getVal(r-1);
				}
				else if(c==1){
					retVal -= matrix_[r][r+1]*recent.getVal(r+1);
				}
			
			}

		}
		else{
			for(int c=0; c<numCols_; c++){
				if(c<r){
					retVal -= matrix_[r][c]*ret.getVal(c);
				}
				else if(c>r){
					retVal -= matrix_[r][c]*recent.getVal(c);
				}
			}
		}

		if(matrix_[r][r] == 0){
			cout << "error: diving by zero, matrix cannot be solved" << endl;
			exit(1);
		}
		ret.setVal(r, retVal/matrix_[r][r]);

	}

	return ret;

}



double Matrix:: getSolveTime(){
	return solveTime_;
}

void Matrix:: setSolveDebugMode(bool flag){
	debug_ = flag;
}

bool Matrix:: getDebugMode(){
	return debug_;
}

void Matrix:: setSolveMaxIters(int iters){
	maxIters_ = iters;
}

void Matrix:: setSolveTolerance(double tol){
	tolerance_ = tol;
}

int Matrix:: getSolveIters(){
	return iters_;
}


//hw05


Vector Matrix:: linFit(int mode){

	Vector ret;

	if(mode == LINEAR){
		//compute a1
		ret.allocateData(2);
		Vector x = extract(0);
		Vector y = extract(1);
		double num = (numRows_*x.sum(y))-(x.sum()*y.sum());
		double denom = (numRows_*x.sumSquared())-(x.sum()*x.sum());
		double a1 = num/denom;
	//	double a1 = ((numRows_*x.sum(y))-(x.sum()*y.sum()))/(numRows_*x.sumSquared()-(x.sum()*x.sum()));
		double a0 = (1/double(numRows_))*(y.sum()-(a1*x.sum()));
		cout << "a0: " << a0 << endl;
		ret.setVal(0, a0);
		ret.setVal(1, a1);

	}

	else if(mode == EXPONENTIAL){
		ret.allocateData(2);
		Vector x = extract(0);
		Vector y = extract(1);
		for(int i=0; i<y.numElems(); i++){
			double yprev = y.getVal(i);
			y.setVal(i, log(yprev));
		}
		double num = (numRows_*x.sum(y))-(x.sum()*y.sum());
		double denom = (numRows_*x.sumSquared())-(x.sum()*x.sum());
		double a1 = num/denom;
	//	double a1 = ((numRows_*x.sum(y))-(x.sum()*y.sum()))/(numRows_*x.sumSquared()-(x.sum()*x.sum()));
		double a0 = (1/double(numRows_))*(y.sum()-(a1*x.sum()));
		ret.setVal(0, a0);
		ret.setVal(1, a1);

	}
	else if(mode == POWER){
		ret.allocateData(2);
		Vector x = extract(0);
		for(int i=0; i<x.numElems(); i++){
			double xprev = x.getVal(i);
			x.setVal(i, log10(xprev));
		}
		Vector y = extract(1);
		for(int j=0; j<y.numElems(); j++){
			double yprev = y.getVal(j);
			y.setVal(j, log10(yprev));
		}
		double num = (numRows_*x.sum(y))-(x.sum()*y.sum());
		double denom = (numRows_*x.sumSquared())-(x.sum()*x.sum());
		double a1 = num/denom;
	//	double a1 = ((numRows_*x.sum(y))-(x.sum()*y.sum()))/(numRows_*x.sumSquared()-(x.sum()*x.sum()));
		double a0 = (1/double(numRows_))*(y.sum()-(a1*x.sum()));
		ret.setVal(0, a0);
		ret.setVal(1, a1);
		
	}

	else{

		cout << "error: not a valid mode" << endl;
		exit(1);

	}

	return ret;

}

Matrix Matrix:: evalLinFit(int mode, Vector fit){
	Matrix ret;
	ret.initSize(numRows_, numCols_);
	
	if(mode == LINEAR){
		for(int i=0; i<numRows_; i++){
			ret.setVal(i, 0, matrix_[i][0]);
			ret.setVal(i, 1, (matrix_[i][1]*fit.getVal(1) + fit.getVal(0)));
		}
	}

	else if(mode == EXPONENTIAL){
		for(int i=0; i< numRows_; i++){
			ret.setVal(i, 0, matrix_[i][0]);
	//		cout << "e^a0: " << exp(fit.getVal(0)) << endl;
	//		cout << "a1: " << fit.getVal(1) << endl;
	//		cout << "x: " << matrix_[i][1] << endl;
	//		cout << "e^a1x: " << exp(fit.getVal(1)*matrix_[i][1]) << endl;
			double yval = exp(fit.getVal(0)) * exp(fit.getVal(1)*matrix_[i][0]);
			ret.setVal(i, 1, yval);
		}
	}

	else if(mode == POWER){
		for(int i=0; i< numRows_; i++){
			ret.setVal(i, 0, matrix_[i][0]);
			double coeff = pow(10, fit.getVal(0));
			double mult = pow(matrix_[i][0], fit.getVal(1));
			ret.setVal(i, 1, coeff*mult);
		}
	}
	else{
		cout << "error: not a valid mode" << endl;
		exit(1);
	}
	
	return ret;
}


Vector Matrix:: extract(int col){
	Vector ret;
	ret.allocateData(numRows_);
	for(int r=0; r<numRows_; r++){
		ret.setVal(r, matrix_[r][col]);
	}
	return ret;
}

void Matrix:: saveToFile(string fileName){
	ofstream myfile;
	myfile.open(fileName);
	myfile << numRows_ << " " << numCols_ << endl;
	for(int r=0; r<numRows_; r++){
		for(int c=0; c<numCols_; c++){
			myfile << matrix_[r][c] << " ";
		}
		myfile << endl;
	}
}
