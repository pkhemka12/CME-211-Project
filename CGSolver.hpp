#ifndef CGSOLVER_HPP
#define CGSOLVER_HPP

#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include "matvecops.hpp"
#include "sparse.hpp"


/* Function that implements the CG algorithm for a linear system
 *
 * Ax = b
 *
 * where A is in CSR format.  The starting guess for the solution
 * is provided in x, and the solver runs a maximum number of iterations
 * equal to the size of the linear system.  Function returns the
 * number of iterations to converge the solution to the specified
 * tolerance, or -1 if the solver did not converge.
 */
int CGSolver(SparseMatrix &A,
	     const std::vector<double> &b,
	     std::vector<double> &x,
	     double              tol,
	     std::string soln_prefix,
	     const std::vector<double> &T_x,
	     double Th);

/* Function that implements writing solution to output files
   where soln_prefix is received from argument line, niter is
   the number of iterations for the name of the file, Tx is 
   the vector of temps, x is the solution file, and Th is the
   temp of the hot boundary */ 
void WriteSoln(std::string soln_prefix,
	       int niter,
	       const std::vector<double> &T_x,
	       const std::vector<double> &x,
	       double Th);

#endif /* CGSOLVER_HPP */
