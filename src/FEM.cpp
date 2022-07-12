#include <iostream>
#include <cmath>
#include "C:\Users\argir\FEM-implementation-in-CPP\include\Matrix.h"
#include "C:\Users\argir\FEM-implementation-in-CPP\include\Mesh.h"
#include "C:\Users\argir\FEM-implementation-in-CPP\include\FEM.h"

FEM::FEM() {}

FEM::FEM(Mesh& m)
{
	my_mesh = m;
}

FEM::FEM(FEM& b)
{
	my_mesh = b.my_mesh;
	stiff = b.stiff;
	rhs = b.rhs;
}

FEM::~FEM() {}

void FEM::calcStiff()
{
	Matrix A(my_mesh.numNodes, my_mesh.numNodes);
	for (int k = 0; k < my_mesh.numElements; k++)
	{
		int n1 = my_mesh.connectivity[k][0];
		int n2 = my_mesh.connectivity[k][1];
		int n3 = my_mesh.connectivity[k][2];
		calcLocalStiff(my_mesh.nodes[n1], my_mesh.nodes[n2], my_mesh.nodes[n3]);
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				double temp = A(my_mesh.connectivity[k][i], my_mesh.connectivity[k][j]) + localStiff(i, j);
				A.setValue(my_mesh.connectivity[k][i], my_mesh.connectivity[k][j], temp);
			}
		}
	}
	stiff = A;
}

void FEM::calcRhs()
{
	double v1x;
	double v1y;
	double v2x;
	double v2y;
	double v3x;
	double v3y;
	std::vector<double> F(my_mesh.numNodes);

	for (int i = 0; i < my_mesh.numElements; i++)
	{
		int n1 = my_mesh.connectivity[i][0];
		int n2 = my_mesh.connectivity[i][1];
		int n3 = my_mesh.connectivity[i][2];

		v3x = (my_mesh.nodes[n2][0] - my_mesh.nodes[n1][0]);
		v3y = (my_mesh.nodes[n2][1] - my_mesh.nodes[n1][1]);
		v1x = (my_mesh.nodes[n3][0] - my_mesh.nodes[n2][0]);
		v1y = (my_mesh.nodes[n3][1] - my_mesh.nodes[n2][1]);
		v2x = (my_mesh.nodes[n1][0] - my_mesh.nodes[n3][0]);
		v2y = (my_mesh.nodes[n1][1] - my_mesh.nodes[n3][1]);

		double val = 0.5 * abs(-v3x * v2y + v3y * v2x) * 10;

		F[n1] += val;
		F[n2] += val;
		F[n3] += val;
	}

	rhs = F;
}

void FEM::print(std::vector<double>& vec)
{
	for (auto el : vec)
	{
		std::cout << el << std::endl;
	}
}

void FEM::print(std::vector<int>& vec)
{
	for (auto el : vec)
	{
		std::cout << el << std::endl;
	}
}

void FEM::calcLocalStiff(std::vector<double> coords1, std::vector<double> coords2, std::vector<double> coords3)
{
	Matrix Alocal(3, 3);

	// B is the Jacobian
	Matrix B(2, 2);
	B.setValue(0, 0, coords1[0] - coords3[0]);
	B.setValue(0, 1, coords2[0] - coords3[0]);
	B.setValue(1, 0, coords1[1] - coords3[1]);
	B.setValue(1, 1, coords2[1] - coords3[1]);
	B.calcDet();
	double area = 0.5 * B.detVal;
	if (area < 0) area = -area;
	Matrix invTransB = B.transpose().inverse();

	// G is the phij(FT(x))
	Matrix G(2, 3);
	G.setValue(0, 0, 1);
	G.setValue(0, 1, 0);
	G.setValue(0, 2, -1);
	G.setValue(1, 0, 0);
	G.setValue(1, 1, 1);
	G.setValue(1, 2, -1);

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			Matrix G_temp_j(2, 1);
			G_temp_j.setValue(0, 0, G(0, j));
			G_temp_j.setValue(1, 0, G(1, j));
			Matrix BG_j = invTransB * G_temp_j;

			Matrix G_temp_i(2, 1);
			G_temp_i.setValue(0, 0, G(0, i));
			G_temp_i.setValue(1, 0, G(1, i));
			Matrix BG_i = invTransB * G_temp_i;

			double val = area * (BG_j(0, 0) * BG_i(0, 0) + BG_j(1, 0) * BG_i(1, 0));
			Alocal.setValue(i, j, val);
		}
	}
	localStiff = Alocal;
}

void FEM::apply_boundary()
{
	Matrix new_stiff(my_mesh.interior_nodes.size(), my_mesh.interior_nodes.size());
	std::vector<double> new_rhs;
	int i = 0;
	for (auto it1 = my_mesh.interior_nodes.begin(); it1 != my_mesh.interior_nodes.end(); it1++)
	{
		int j = 0;
		for (auto it2 = my_mesh.interior_nodes.begin(); it2 != my_mesh.interior_nodes.end(); it2++)
		{
			new_stiff.setValue(i, j, stiff(*it1, *it2));
			j++;
		}
		new_rhs.push_back(rhs[*it1]);
		i++;
	}	
	stiff = new_stiff;
	rhs = new_rhs;
}

void FEM::discretize()
{
	calcStiff();
	calcRhs();
	apply_boundary();
}

void FEM::expand_sol(std::vector<double>& init_sol)
{
	std::vector<double> final_sol;
	int j = 0;
	for (int i = 0; i < my_mesh.numNodes; i++)
	{
		if (std::count(my_mesh.interior_nodes.begin(), my_mesh.interior_nodes.end(), i))
		{
			final_sol.push_back(init_sol[j]);
			j++;
		}
		else
		{
			final_sol.push_back(0.0);
		}
	}
	solution = final_sol;
}

std::vector<double> FEM::solve(int iters = 10)
{
	solution = stiff.Jacobi_iterator(rhs, iters);
	expand_sol(solution);

	return solution;
}

void FEM::write_solution(std::string& filename)
{
	std::ofstream my_file(filename);
	if (my_file.is_open())
	{
		my_file << "VARIABLES= " << "\"x\"" << "," << "\"y\"" << "," << "\"sol\"" << "\n";
		my_file << "ZONE T= " << "\"Frame 0\"" << ",I=" << my_mesh.numNodes << ",J=" << my_mesh.numNodes << "\n";
		for (int i = 0; i < my_mesh.numNodes; i++)
		{
			my_file << my_mesh.nodes[i][0] << " " << my_mesh.nodes[i][1] << " " << solution[i] << "\n";
		}
		my_file.close();
	}
}