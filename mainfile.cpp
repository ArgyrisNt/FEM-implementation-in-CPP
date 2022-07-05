#include <iostream>
#include <vector>
#include <algorithm>
#include "C:\Users\argir\FEM-implementation-in-CPP\include\Matrix.h"
#include "C:\Users\argir\FEM-implementation-in-CPP\src\Matrix.cpp"
#include "C:\Users\argir\FEM-implementation-in-CPP\include\Mesh.h"
#include "C:\Users\argir\FEM-implementation-in-CPP\src\Mesh.cpp"
#include "C:\Users\argir\FEM-implementation-in-CPP\include\FEM.h"
#include "C:\Users\argir\FEM-implementation-in-CPP\src\FEM.cpp"

// - - - - - - - - - - TO DO - - - - - - - - - -
// 
// Define the calcError function
// Define the refinement class
// Define several options for shape functions

int main()
{
    std::vector<std::vector<double>> my_nodes{ {-1.0, -1.0},
                        {-0.5, -1.0},
                        {0.0, -1.0},
                        {0.0, -0.5},
                        {0.0, 0.0},
                        {0.5, 0.0},
                        {1.0, 0.0},
                        {1.0, 0.5},
                        {1.0, 1.0},
                        {0.5, 1.0},
                        {0.0, 1.0},
                        {-0.5, 1.0},
                        {-1.0, 1.0},
                        {-1.0, 0.5},
                        {-1.0, 0.0},
                        {-1.0, -0.5},
                        {-0.5, -0.5},
                        {-0.5, 0.0},
                        {-0.5, 0.5},
                        {0.0, 0.5},
                        {0.5, 0.5} };

    std::vector<std::vector<int>> my_connectivity{ {1, 15, 0},
                                {15, 1, 16},
                                {2, 16, 1},
                                {16, 2, 3},
                                {16, 14, 15},
                                {14, 16, 17},
                                {3, 17, 16},
                                {17, 3, 4},
                                {18, 14, 17},
                                {14, 18, 13},
                                {19, 17, 4},
                                {17, 19, 18},
                                {11, 13, 18},
                                {13, 11, 12},
                                {10, 18, 19},
                                {18, 10, 11},
                                {5, 19, 4},
                                {19, 5, 20},
                                {6, 20, 5},
                                {20, 6, 7},
                                {20, 10, 19},
                                {10, 20, 9},
                                {7, 9, 20},
                                {9, 7, 8} };

    // - - - - - Construct mesh - - - - -
    Mesh my_mesh(my_nodes, my_connectivity);

    // - - - - - Write mesh data - - - - -
    my_mesh.write_mesh(std::string("mesh_init.obj"));

    // - - - - - Define basis - - - - -
    FEM my_FEM(my_mesh);

    // Compute stiffness matrix and rhs
    my_FEM.discretize();

    // - - - - - Solve - - - - -
    my_FEM.solve();

    // - - - - - Write solution data - - - - -
    my_FEM.write_solution(std::string("sol_init.dat"));

    // - - - - - Uniform refinement - - - - -
    int iters = 1;
    my_mesh.refine(iters);

    // - - - - - Write new mesh data - - - - -
    my_mesh.write_mesh(std::string("mesh_ref.obj"));

    // - - - - - Update basis - - - - -
    my_FEM.my_mesh = my_mesh;

    // Re-compute stiffness matrix and rhs
    my_FEM.discretize();

    // - - - - - Solve again - - - - -
    my_FEM.solve();

    // - - - - - Write new solution data - - - - -
    my_FEM.write_solution(std::string("sol_ref.dat"));

	return 0;
}