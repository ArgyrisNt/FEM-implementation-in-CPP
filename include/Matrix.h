#pragma once

#include <iostream>
#include <vector>

class Matrix
{
public:
	// Constructors
	Matrix();
	Matrix(int my_rows, int my_cols);
	Matrix(const Matrix&);

	// Overloaded operators
	Matrix& operator=(Matrix);
	Matrix operator+(Matrix&);
	Matrix operator-(Matrix&);
	Matrix operator*(Matrix&);
	Matrix operator*(const double&);
	double operator()(int row, int col);

	//Member functions
	void setValue(int row, int col, double val);
	void calcDet();
	void print();
	Matrix transpose();
	Matrix inverse();
	double norm(std::vector<double>& v);
	std::vector<double> linsolve(std::vector<double> b, int iters);

	// Destructor
	~Matrix();

	//Member variables
	int rows;
	int cols;
	double** mat;
	double detVal;
};