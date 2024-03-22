#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "CGSolver.hpp"
#include "matvecops.hpp"
#include "sparse.hpp"

//Finds iterations to convergence
int CGSolver(SparseMatrix &A,
	     const std::vector<double> &b,
	     std::vector<double> &x,
	     double              tol,
	     std::string soln_prefix,
	     const std::vector<double> &T_x,
	     double Th){

  //Initializing variable types
  int niter, nitermax;
  std::vector<double> r, p, Ax;
  double alpha, beta;

  //Assigning/calculating variables needed for loop calc
  nitermax = int(x.size());
  Ax = A.MulVec(x);
  //Pass -1 as the last argument to subtract
  r = add(b, Ax, -1.);
  double normr0 = L2norm(r);
  p = r;
  niter = 0;

  //Write first guess solution to file
  WriteSoln(soln_prefix, niter, T_x, x, Th);

  //Iterating until solution converges/niter exceeds matrix dim
  while(niter < nitermax){
    niter++;
    //Calculating vec dot vec and mat dot vec for req vars
    double rr = dotprod(r, r);
    std::vector<double> Ap = A.MulVec(p);
    double pAp = dotprod(p, Ap);

    //Computing alpha and updating x, r
    alpha = rr/pAp;
    //Pass (-)alpha into add to perform scalar product and addition
    x = add(x, p, alpha);
    r = add(r, Ap, -alpha);

    //Computing the new norm
    double normr = L2norm(r);

    //Checks if niter is a multile of 10 and writes soln to file
    if (niter % 10 == 0){
      WriteSoln(soln_prefix, niter, T_x, x, Th);
    }

    //Checks if convergence found
    if ((normr/normr0) < tol) break;
    
    //Calculating beta and updating p for next iteration
    beta = dotprod(r,r)/rr;
    //Pass beta into add to perform scalar product and addition
    p = add(r, p, beta);
  }

  //Checks whether solution actually converged or not
  if (niter <= nitermax){
    //Writes to file
    WriteSoln(soln_prefix, niter, T_x, x, Th);
    return niter;
  }
  else{
    return -1;
  }
}

//Writes to solution files & outputs
void WriteSoln(std::string soln_prefix,
	       int niter,
	       const std::vector<double> &T_x,
	       const std::vector<double> &x,
	       double Th){

  //Open string stream for output filenames
  std::stringstream os;

  //Creating the correct name for files
  os << std::setw(3) << std::setfill('0') << niter;
  soln_prefix += os.str();
  soln_prefix += ".txt";

  //Open output stream
  std::ofstream outf(soln_prefix);

  if (outf.is_open()){
    //Write hot isothermal boundary
    for (int i = 0; i < (int)T_x.size(); i++){
      outf << Th << std::endl;
    }

    //Write Solution
    for (int i = 0; i < (int)x.size(); i++){
      outf << x[i] << std::endl;
    }

    //Write hot iosthermal boundary
    for (int i = 0; i < (int)T_x.size(); i++){
      outf << T_x[i] << std::endl;
    }

    outf.close();
  }
  //Output file error
  else{
    std::cerr << "ERROR: Cannot open solution file" << std::endl;
    exit(1);
  }
}

