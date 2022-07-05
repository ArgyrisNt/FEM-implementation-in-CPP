#include <iostream>
#include <algorithm>
#include <fstream>
#include "C:\Users\argir\FEM-implementation-in-CPP\include\Mesh.h"

Mesh::Mesh()
{
	std::vector<std::vector<double>> my_nodes;
	nodes = my_nodes;
	numNodes = 0;

	std::vector<std::vector<int>> my_connectivity;
	connectivity = my_connectivity;
	numElements = 0;
}

Mesh::Mesh(std::vector<std::vector<double>>& my_nodes, std::vector<std::vector<int>>& my_connectivity)
{
	nodes = my_nodes;
	connectivity = my_connectivity;
	numNodes = nodes.size();
	numElements = connectivity.size();
	calcInteriorNodes();
}

Mesh::Mesh(const Mesh& m)
{
	nodes = m.nodes;
	connectivity = m.connectivity;
	numNodes = m.numNodes;
	numElements = m.numElements;
}

Mesh::~Mesh() {}

void Mesh::calcInteriorNodes()
{
	std::vector<std::vector<int>> edges;
	for (int i = 0; i < numElements; i++)
	{
		int n1 = connectivity[i][0];
		int n2 = connectivity[i][1];
		int n3 = connectivity[i][2];
		int count1 = 0;
		int count2 = 0;
		int count3 = 0;
		for (int j = 0; j < edges.size(); j++)
		{
			if ((edges[j][1] == n1 && edges[j][2] == n2) || (edges[j][1] == n2 && edges[j][2] == n1))
			{
				edges[j][0] += 1;
				count1++;
			}
			if ((edges[j][1] == n2 && edges[j][2] == n3) || (edges[j][1] == n3 && edges[j][2] == n2))
			{
				edges[j][0] += 1;
				count2++;
			}
			if ((edges[j][1] == n3 && edges[j][2] == n1) || (edges[j][1] == n1 && edges[j][2] == n3))
			{
				edges[j][0] += 1;
				count3++;
			}
		}
		if (!count1)
		{
			edges.push_back({ 1,n1,n2 });
		}
		if (!count2)
		{
			edges.push_back({ 1,n2,n3 });
		}
		if (!count3)
		{
			edges.push_back({ 1,n3,n1 });
		}
	}

	std::vector<int> temp_bound;
	for (auto el1 : edges)
	{
		if (el1[0] == 1)
		{
			temp_bound.push_back(el1[1]);
			temp_bound.push_back(el1[2]);
		}
	}
	sort(temp_bound.begin(), temp_bound.end());
	auto it1 = unique(temp_bound.begin(), temp_bound.end());
	temp_bound.resize(distance(temp_bound.begin(), it1));
	std::vector<int> boundary_nodes = temp_bound;

	interior_nodes = {};
	for (int i = 0; i < numNodes; i++)
	{
		if (!std::count(boundary_nodes.begin(), boundary_nodes.end(), i))
		{
			interior_nodes.push_back(i);
		}
	}
}

void Mesh::refine(int iters)
{
	for (int iter = 0; iter < iters; iter++)
	{
		std::vector<std::vector<int>> new_connectivity;
		int temp = numElements;
		for (int i = 0; i < temp; i++)
		{
			// Compute nodes of the element
			int n1 = connectivity[i][0];
			int n2 = connectivity[i][1];
			int n3 = connectivity[i][2];

			// Compute mid points on each edge
			std::vector<double> mid12{ (nodes[n1][0] + nodes[n2][0]) / 2,(nodes[n1][1] + nodes[n2][1]) / 2 };
			std::vector<double> mid23{ (nodes[n2][0] + nodes[n3][0]) / 2,(nodes[n2][1] + nodes[n3][1]) / 2 };
			std::vector<double> mid31{ (nodes[n3][0] + nodes[n1][0]) / 2,(nodes[n3][1] + nodes[n1][1]) / 2 };
			std::vector<std::vector<double>> mids{ mid12,mid23,mid31 };

			// Check if mid points already exist as new nodes
			std::vector<int> ids{ 0,0,0 };
			for (int i = 0; i < mids.size(); i++)
			{
				auto it = find(nodes.begin(), nodes.end(), mids[i]);
				if (it != nodes.end())
				{
					ids[i] = it - nodes.begin();
				}
				else
				{
					nodes.push_back(mids[i]);
					ids[i] = nodes.size() - 1;
				}
			}

			// Insert new elements
			std::vector<int> elem1{ ids[0],n2,ids[1] };
			std::vector<int> elem2{ ids[1],n3,ids[2] };
			std::vector<int> elem3{ ids[2],n1,ids[0] };
			std::vector<int> elem4{ ids[0],ids[1],ids[2] };
			new_connectivity.push_back(elem1);
			new_connectivity.push_back(elem2);
			new_connectivity.push_back(elem3);
			new_connectivity.push_back(elem4);

			numNodes = nodes.size();
			numElements = new_connectivity.size();
		}
		connectivity = new_connectivity;
		calcInteriorNodes();
	}
}

void Mesh::write_mesh(std::string& filename)
{
	std::ofstream my_file(filename);
	if (my_file.is_open())
	{
		my_file << "# " << numNodes << "\n";
		my_file << "# " << numElements << "\n";
		for (int i = 0; i < numNodes; i++)
		{
			my_file << "v " << nodes[i][0] << " " << nodes[i][1] << "\n";
		}
		for (int i = 0; i < numElements; i++)
		{
			my_file << "f " << connectivity[i][0] + 1 << " " << connectivity[i][1] + 1 << " " << connectivity[i][2] + 1 << "\n";
		}
		my_file.close();
	}
}