
//Sydney Loats SML4267

#include <iostream>
#include <vector>
#include <math.h>
#include <string>
#include <limits.h>
#pragma once
#include "vector.hpp"

using namespace std;

enum solveModes {GAUSS_ELIM,JACOBI,GAUSS_SEIDEL};
enum fitModes {LINEAR=1, EXPONENTIAL=2, POWER=3};

class Matrix

{

  private:
	int numRows_ = -1;
	int numCols_ = -1;
	vector<vector<double>> matrix_;
	int iters_ = 1;
	double tolerance_ = 1*pow(10, -8);
	double solveTime_;
	int maxIters_ = INT_MAX;
	bool debug_ = false;
	bool turbo_ = false;

  public:
	Matrix();                                          //constructor
	~Matrix();                                         //destructor
	void checkInit();                                  //ensures the matrix has been initialized before the function is called
	void initIdentity(int n);                          //initialize as an identity matrix of size nxn
	void initFromFile(string fileName);                //read the matrix in from a file filename
	void initSize(int r, int c);                       //initialize matrix to size rxc
	bool isSquare();                                   //test wether matrix is square
	int numRows();                                     //return number of rows
	int numCols();                                     //return number of columns
	void setRows(int r);                               //set number of rows to r
	void setCols(int c);                               //set number of columns to c
	double getVal(int row, int col);                   //return the matrix value at given row/col location (0-index based)
	void setVal(int row, int col, double val);         //set the matrix value at given row/col location (0-index based)
	Matrix Multiply(Matrix B);                         //post-multiply by B and return resulting matrix
	double helpMultiply(int r, int c, Matrix B);       //function to help the Multiply(Matrix B) function
	Matrix Multiply(double A);                         //multiply by a scalar and return resulting matrix
	Matrix Transpose();                                //return transpose of the matrix
	vector<double> Diagonal();                         //returns a vector containing diagonal elements of the matrix
	void Print();                                      //print the matrix
	void Print(string name);                           //print the matrix with the name name
	int getMaxIters();                                 //returns maximum iterations


	//linear solver support
	
	Vector Solve(Vector b, int mode);                  //return solution of [A]x = b using desired solver mode
	double getSolveTime();                             //wall clock time (in secs) required for last solve
	bool getDebugMode();                               //return debug mode
	void setSolveDebugMode(bool flag);                 //set flag to toggle debug output mode for linear solves

	//iterative methods solver support
	
	void setSolveMaxIters(int iters);                  //set cap on max # of iterations
	void setSolveTolerance(double tol);                // set desired stopping tolerance
	int  getSolveIters();                              // return number of iters completed from last solve

	//helper methods
	Vector createVector(int row);                      //creates a vector from row row
	void replaceRow(int r, Vector b);                  //replaces a row of the matrix with vector b
	void addColumn(Vector b);                          //adds column vector b to the end of the matrix
	Matrix reducedMatrix(Vector b);                    //returns reduced matrix
	void setMatrix(vector< vector<double>> m);         //sets matrix
	Vector solveMatrix();                              //helper for gaussian elimination
	Vector solveJacobi(Vector b, Vector recent);       //helper for jacobi
	Vector solveGaussSeidel(Vector b, Vector recent);  //helper for gauss seidel
	Matrix separateMatrix();                           //separates matrix from the reduced matrix
	Vector separateVector();                           //separates vector from the reduced matrix



	//homework 5 new functions
	//regression support (for Nx2 matrices)
	Vector linFit(int mode);                           //linear regression fit (mode=LINEAR, EXPONENTIAL, or POWER)
	Matrix evalLinFit(int mode, Vector fit);     // evaluate linear fit - ouput is an Nx2 matrix

	// additional utilities to aid in regression
	Vector extract(int col);                     // extract vector from matrix corresponding to the col index
	void saveToFile(string fileName);            // save matrix to file (same file format as initFromFile)

	//proj02
	void setMatrixTurbo(bool t);

};
