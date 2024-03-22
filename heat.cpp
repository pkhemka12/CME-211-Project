#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "CGSolver.hpp"
#include "heat.hpp"
#include "sparse.hpp"

/* Method to setup Ax=b system */
int HeatEquation2D::Setup(std::string inputfile){
 
  int nx, ny; 
  double len, width, h, Tc;

  /* Creating input stream */
  std::ifstream inputf(inputfile);

  /* Opening file and reading in data */
  if(inputf.is_open()){
    inputf >> len >> width >> h;
    inputf >> Tc >> Th;

    inputf.close(); 
  }
  /* Error message when file fails to open */
  else{
    std::cerr << "Input file can't be read/opened!" << std::endl;
    return -1;
  }
  
  /* Check to make sure length and width for matrix dim are ints */
  if (floor(len/h) == len/h && floor(width/h) == width/h){
    nx = int(len/h);
    ny = int(width/h - 1);
  }
  /* Error message to prevent population of matrix with irrregular dims */
  else{
    std::cerr << "Length or width is not an integer" << std::endl;
    return -1;
  }

  /* Resize matrix and soln vector to correct dims */
  A.Resize(nx*ny, nx*ny);
  x.resize(nx*ny, 1.);
  
  /* Pushes Tx for x = 0 into the T_x vector */
  T_x.push_back(Tx(0., Tc, len));

  /* Filling matrix A */
  for(int j = 0; j < ny; j++){
    for(int i = 0; i < nx; i++){
      /* Var to relate loop vars to positions in stencil */ 
      int idxA = j * nx + i;
      /* Add self point */
      A.AddEntry(idxA, idxA, 4.);
      /* Check if pos not on right boundary */
      if (i < nx - 1){
	/* Add right point */
	A.AddEntry(idxA, idxA + 1, -1.);
	/* Check if pos on left boundary */
        if (i == 0){
	  /* Add "left" point, which wraps to right side */
	  A.AddEntry(idxA, idxA + nx - 1, -1.);
	}
	else{
	  /* Add left point */
	  A.AddEntry(idxA, idxA - 1, -1.);
	}
      }
      /* Pos on right boundary */
      else{
	/* Add "right" point aka self point */
	A.AddEntry(idxA, idxA, -1.);
	/* Add left point */
	A.AddEntry(idxA, idxA - 1, -1.);
      }

      /* Check if pos not on top boundary */
      if (j < ny - 1){
	/* Add point above */
	A.AddEntry(idxA, idxA + nx, -1.);
	/* Pos on bottom boundary */
	if (j == 0){
	  /* Calc x based on i & T(x) for cold boundary */
	  double x_pt = h * i;
	  double tx = Tx(x_pt, Tc, len);
	  b.push_back(tx);
	  T_x.push_back(tx);
	}
	else{
	  /* Add point below */
	  A.AddEntry(idxA, idxA - nx, -1.);
	  b.push_back(0);
	}
      }
      /* Pos on top boundary */ 
      else{
	/* Add point below */
	A.AddEntry(idxA, idxA - nx, -1.);
	b.push_back(Th);
      }
    }
  }
  /* Convert A to CSR */
  A.ConvertToCSR();
  
  return 0; 
}

/* Method to solve system using CGsolver */
int HeatEquation2D::Solve(std::string soln_prefix){

  /* Initialize tolerance */
  double tol = 1.e-5;
  /* Find number of iterations for convergence */
  int niter = CGSolver(A, b, x, tol, soln_prefix, T_x, Th);

  /* Checks if convergence and outputs */
  if (niter == -1){
    std::cout << "Solution did not converge!" << std::endl;
    return -1;
  }
  else{
    std::cout << "SUCCESS: CG solver converged in " << niter;
    std::cout << " iterations." << std::endl;
    return 0;
  }
}

/* Compute temperature at cold boundary */
double HeatEquation2D::Tx(double x_pt, double Tc, double len){
  /* Finds exponent and calculates Tx */
  double expon = -10. * pow((x_pt - len/2), 2);
  double tx = -Tc * (exp(expon) - 2); 

  return tx;
}	
 
