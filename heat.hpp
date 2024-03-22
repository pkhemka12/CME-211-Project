#ifndef HEAT_HPP
#define HEAT_HPP

#include <string>
#include <vector>

#include "sparse.hpp"

class HeatEquation2D
{
  private:
    SparseMatrix A;
    std::vector<double> b, x, T_x;
    double Th; 
  public:
    /* Method to setup Ax=b system */
    int Setup(std::string inputfile);

    /* Method to solve system using CGsolver */
    int Solve(std::string soln_prefix);

    /* Method to compute temperature at cold boundary */ 
    double Tx(double x_pt, double len, double Tc); 
};

#endif /* HEAT_HPP */
