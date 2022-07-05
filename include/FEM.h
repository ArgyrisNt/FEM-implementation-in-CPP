#pragma once
#include "C:\Users\argir\FEM-implementation-in-CPP\include\Matrix.h"
#include "C:\Users\argir\FEM-implementation-in-CPP\include\Mesh.h"

class FEM
{
public:
	// Constructors
	FEM();
	FEM(Mesh&);
	FEM(FEM&);

	//Destructor
	~FEM();

	// Member variables
	Mesh my_mesh;
	Matrix stiff;
	std::vector<double> rhs;
	std::vector<double> solution;

	// Member functions
	void discretize();
	void print(std::vector<double>&);
	void print(std::vector<int>&);
	std::vector<double> solve(int);
	void write_solution(std::string&);

private:
	// Member functions
	void calcStiff();
	void calcRhs();
	void calcLocalStiff(std::vector<double>, std::vector<double>, std::vector<double>);
	void apply_boundary();
	void expand_sol(std::vector<double>&);

	Matrix localStiff;
};