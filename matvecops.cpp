#include <cmath>
#include <iostream>
#include <vector>

#include "matvecops.hpp"

//Calculates matrix dot vector
std::vector<double> sp_matdotvec(const std::vector<double> &val,
				 const std::vector<int> &row_ptr,
				 const std::vector<int> &col_idx,
				 const std::vector<double> &x){

  //Initializing output vector
  std::vector<double> b;

  //Performing the product for a sparse matrix
  for (unsigned int i = 0; i < row_ptr.size() - 1; i++){
    double count = 0;
    for (int j = row_ptr[i]; j < row_ptr[i+1]; j++){
      count += val[j] * x[col_idx[j]];
    }
    b.push_back(count);
  }
  return b;
}

//Calculates addition of vectors with scalar product
std::vector<double> add(const std::vector<double> &a,
			const std::vector<double> &b, double c){
  
  //Initializing output vector
  std::vector<double> sumvec;

  //Calculation of the sum
  for (unsigned int i = 0; i < a.size(); i++){
    double sum = a[i] + c * b[i];
    sumvec.push_back(sum);
  }
  return(sumvec);
}

//Calculates dotproduct of vectors
double dotprod(const std::vector<double> &a, const std::vector<double> &b){

  //Initializing the output variable
  double prodTot = 0;

  //Calculating the dotproduct
  for(unsigned int i = 0; i < a.size(); i++){
    prodTot += a[i] * b[i];
  }
  return prodTot; 
}

//Calculates the L2 norm of a vector
double L2norm(const std::vector<double> &a){

  //Computes square root of the dot product
  return(std::sqrt(dotprod(a, a)));
}
