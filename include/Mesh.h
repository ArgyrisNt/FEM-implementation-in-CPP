#pragma once
#include <iostream>
#include <vector>

class Mesh
{
public:
	// Constructors
	Mesh();
	Mesh(std::vector<std::vector<double>>&, std::vector<std::vector<int>>&);
	Mesh(const Mesh&);

	// Destructor
	~Mesh();

	// Member functions
	void calcInteriorNodes();
	void refine(int);
	void write_mesh(std::string&);

	// Member variables
	std::vector<std::vector<double>> nodes;
	std::vector<std::vector<int>> connectivity;
	int numNodes;
	int numElements;
	std::vector<int> interior_nodes;
};