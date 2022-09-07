//Sydney Loats SML4267

#include <vector>
#include <string>
#pragma once
using namespace std;


class Vector
{
private:
	int length_ = -1;
	vector<double> vector_;
	vector<double> previous_;
public:
	Vector();                                  // constructor
	~Vector();                                 // destructor
	void allocateData(int nelems);             // allocate space for nelem entries (of double type)
        void initFromFile(string fileName);        // read the vector from a file
        int numElems();                            // return number of elements
        double getVal(int i);                      // return the ith element
        double l2norm();                           // return l2 norm
        void setVal(int i, double val);            // set the ith element to val
        void setAllVals(double val);               // set all elements of the vector to val
        void Print();                              // print the vector contents to stdout
	void Print(string name);                   // print the vector contents to stdout with a name prefix
	Vector subtractFrom(Vector b);             //subtract current vector from vector b
	Vector multiply(double val);               //multiplies whole vector by constant
	void setPrevious(vector<double>);          //assigns previous vector value for interative methods
	vector<double> getVector();                //returns vector
	double norm();                             //helper for l2 norm

	//regression support
	double sum();                              //return the sum of a vector
	double sumSquared();                       //return the sum of (vector elements)^2
	double sum(Vector y);                      //return the sum of (current vector)*(y vector) elements
};
