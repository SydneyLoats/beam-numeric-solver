//Sydney Loats SML4267

#include <iostream>
#include <math.h>
#include <fstream>
#include "vector.hpp"

using namespace std;

Vector:: Vector(){
//	cout << "constructor: entering the program" << endl;
}

Vector::~Vector(){
//	cout << "destructor: leaving the program" << endl;
}

void Vector:: allocateData(int nelems){
	vector_.resize(nelems);
	length_ = nelems;
}

void Vector:: initFromFile(string filename){
	
	ifstream fin(filename);
	int l;
	fin >> l;
	length_ = l;
	double d;
	allocateData(l);
	for(int i=0; i< l; i++){
		fin >> d;
		vector_[i] = d;
	}
}

Vector Vector:: subtractFrom(Vector b){
	Vector ret;
	ret.allocateData(length_);
	for(int i=0; i<length_; i++){
		double val = b.getVal(i) - vector_[i];
		ret.setVal(i,val); 
	}
	return ret;
}

int Vector:: numElems(){
	return length_;
}

vector<double> Vector:: getVector(){
	return vector_;
}

void Vector:: setPrevious(vector<double> b){
	previous_ = b;
}

double Vector:: getVal(int i){
	return vector_[i];
}

double Vector:: norm(){
	double root_of = 0;
	for(int i=0; i<length_; i++){
		root_of += pow(vector_[i], 2);
	}
	double ret = pow(root_of, 0.5);
	return ret;
}

double Vector:: l2norm(){

	Vector sub;
        sub.allocateData(length_);
	for(int i=0; i<length_; i++){
		sub.setVal(i, vector_[i]-previous_[i]);
	}

	double num = sub.norm();
	double denom = norm();
	return num/denom;
}

void Vector:: setVal(int i, double val){
	
	if(i > length_-1){
		cout << "error: invalid index" << endl;
		exit(1);
	}
	else{
		vector_[i] = val;
	}
}

void Vector:: setAllVals(double val){
	
	for(int i=0; i< length_; i++){
		vector_[i] = val;
	}
}

Vector Vector:: multiply(double val){
	Vector ret;
	ret.allocateData(length_);
	for(int i=0; i< length_; i++){
		double newVal = vector_[i]*val;
		ret.setVal(i, newVal);
	}
	return ret;
}

void Vector:: Print(){
	
	for(int i=0; i<length_; i++){
		cout << "| " << vector_[i] << " |" << endl;
	}
}

void Vector:: Print(string name){

	cout << name << " = " << endl;
	for(int i=0; i<length_; i++){
		cout << "| " << vector_[i] << " |" << endl;
	}
}

//hw05 additions

double Vector:: sum(){
	double ret = vector_[0];
	for(int i=1; i<length_; i++){
		ret += vector_[i];
	}
	return ret;
}

double Vector:: sumSquared(){
	double ret = pow(vector_[0], 2);
	for(int i=1; i<length_; i++){
		ret+= pow(vector_[i], 2);
	}
	return ret;
}

double Vector:: sum(Vector y){
	double ret = vector_[0] * y.getVal(0);
	for(int i=1; i<length_; i++){
		ret+= vector_[i] * y.getVal(i);
	}
	return ret;
}
