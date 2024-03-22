#ifndef MATVECOPS_HPP
#define MATVECOPS_HPP

#include <vector>

/* Calculates the matrix vector product under the assumption that the 
dimensions match 
Arguments: 
val: Values of matrix (vector of doubles)
row_ptr: CSR row pointers (vector of ints)
col_idx: Col indices (vector of ints)
x : vector for multiplication (vector of doubles)
Output: b: product vector (vector of doubles) */  
std::vector<double> sp_matdotvec(const std::vector<double> &val,
				 const std::vector<int> &row_ptr,
				 const std::vector<int> &col_idx,
				 const std::vector<double> &x);

/* Calculates addition or subtraction of 2 vectors by using scalar
product as part of the calculation, assuming size of vectors matches
Arguments:
a,b: 2 input vectors for operation (vectors of doubles)
c: a constant for scalar product with b (double)
Output: summed vector (vector of doubles) */
std::vector<double> add(const std::vector<double> &a,                                   
			const std::vector<double> &b, double c);

/* Calculates the dot product of 2 vectors assuming their sizes match
Arguments: 
a, b: 2 input vectors for operation (vectors of doubles)
Output: dotproduct, a scalar quantity (double)  */
double dotprod(const std::vector<double> &a, const std::vector<double> &b);

/* Calculates the norm of a vector, using the dot product of that vector
and then taking its square root 
Argument:
a: 2 input vectors for operation (vectors of doubles)                                       
Output: norm, a scalar quantity (double) */
double L2norm(const std::vector<double> &a);

#endif /* MATVECOPS_HPP */
