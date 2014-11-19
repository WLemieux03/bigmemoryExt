#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <math.h>
#include <cstddef>

#include <boost/lexical_cast.hpp>

#include "bigmemory/BigMatrix.h"
#include "bigmemory/MatrixAccessor.hpp"
#include "bigmemory/util.h"
#include "bigmemory/isna.hpp"

#include <stdio.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <stdlib.h>
#include <sys/types.h>

//#include <Rcpp.h>
//
//using namespace Rcpp;
using namespace std;

/* 
 * ===  TEMPLATE  ======================================================================
 *         Name:  transposeBM
 *  Description:  This produces a transpose of a big.matrix in a new object.  Expected 
 to be used when you need the original structure and transpose in two separate objects
 for simultaneous use.
 * =====================================================================================
 */
template<typename in_CType, typename in_BMAccessorType, 
  typename out_CType, typename out_BMAccessorType>
void transposeBM(BigMatrix *pInMat, BigMatrix *pOutMat, SEXP rowInds, SEXP colInds)
{
  in_BMAccessorType inMat( *pInMat );
  out_BMAccessorType outMat( *pOutMat );
  
  index_type nRows = GET_LENGTH(rowInds);
  index_type nCols = GET_LENGTH(colInds);
  
  if (nRows != pOutMat->ncol())
    Rf_error("length of row indices does not equal # of rows in new matrix");
  if (nCols != pOutMat->nrow())
    Rf_error("length of col indices does not equal # of cols in new matrix");
  
  index_type i = 0;
  index_type j = 0;
  
  for(i = 0; i < nRows; i++) {
    for(j=0; j < nCols; j++){
      outMat[i][j] = inMat[j][i];
    }
  }
  
  return;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  CtransposeMatrix
 *  Description:  This provides the SEXP for the transposeBM template.  This allows the
 function to be called in R with '.Call'.
 * =====================================================================================
 */
extern "C"
{
  #define CALL_transpose_2(IN_CTYPE, IN_ACCESSOR, OUT_ACCESSOR) \
    switch(pOutMat->matrix_type()) \
    { \
      case 1: \
        transposeBM<IN_CTYPE, IN_ACCESSOR<IN_CTYPE>, char, OUT_ACCESSOR<char> >( \
          pInMat, pOutMat, rowInds, colInds); \
        break; \
      case 2: \
        transposeBM<IN_CTYPE, IN_ACCESSOR<IN_CTYPE>, short, OUT_ACCESSOR<short> >( \
          pInMat, pOutMat, rowInds, colInds); \
        break; \
      case 4: \
        transposeBM<IN_CTYPE, IN_ACCESSOR<IN_CTYPE>, int, OUT_ACCESSOR<int> >( \
          pInMat, pOutMat, rowInds, colInds); \
        break; \
      case 8: \
        transposeBM<IN_CTYPE, IN_ACCESSOR<IN_CTYPE>, double, OUT_ACCESSOR<double> >( \
          pInMat, pOutMat, rowInds, colInds); \
        break; \
    }

  #define CALL_transpose_1(IN_ACCESSOR, OUT_ACCESSOR) \
    switch(pInMat->matrix_type()) \
    { \
      case 1: \
        CALL_transpose_2(char, IN_ACCESSOR, OUT_ACCESSOR) \
        break; \
      case 2: \
        CALL_transpose_2(short, IN_ACCESSOR, OUT_ACCESSOR) \
        break; \
      case 4: \
        CALL_transpose_2(int, IN_ACCESSOR, OUT_ACCESSOR) \
        break; \
      case 8: \
        CALL_transpose_2(double, IN_ACCESSOR, OUT_ACCESSOR) \
        break; \
    }
      
  SEXP CtransposeMatrix(SEXP inAddr, SEXP outAddr, SEXP rowInds, SEXP colInds, 
    SEXP typecast_warning)
  {
    BigMatrix *pInMat = reinterpret_cast<BigMatrix*>(
      R_ExternalPtrAddr(inAddr));
    BigMatrix *pOutMat = reinterpret_cast<BigMatrix*>(
      R_ExternalPtrAddr(outAddr));
    
    if ((pOutMat->matrix_type() < pInMat->matrix_type()) & 
      (LOGICAL_VALUE(typecast_warning) == (Rboolean)TRUE))
    {
      string type_names[9] = {
        "", "char", "short", "", "integer", "", "", "", "double"};
      
      std::string warnMsg = string("Assignment will down cast from ") + 
        type_names[pInMat->matrix_type()] + string(" to ") + 
        type_names[pOutMat->matrix_type()] + string("\n") + 
        string("Hint: To remove this warning type: ") + 
        string("options(bigmemory.typecast.warning=FALSE)");
      Rf_warning(warnMsg.c_str());
    }
    
    // Not sure if there is a better way to do these function calls
    if (pInMat->separated_columns() && pOutMat->separated_columns()) {
      CALL_transpose_1(SepMatrixAccessor, SepMatrixAccessor)
    }
    else if(pInMat->separated_columns() && !(pOutMat->separated_columns()))
    {
      CALL_transpose_1(SepMatrixAccessor, MatrixAccessor)
    }
    else if(!(pInMat->separated_columns()) && pOutMat->separated_columns())
    {
      CALL_transpose_1(MatrixAccessor, SepMatrixAccessor)
    }
    else
    {
      CALL_transpose_1(MatrixAccessor, MatrixAccessor)
    }

    return R_NilValue;
  }
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  gcd
 *  Description:  This calculates the greatest common denominator between two integers. 
 * =====================================================================================
 */

int gcd(int a, int b) 
{
	  int c;

		while (a != 0) {
				c = a;
				a = b%a;
				b = c;
		}
		return b;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  r
 *  Description:  This calculates the indices to remove conflicts.  This was developed 
 from the code algorithm described in:
 "Catanzaro, Bryan, Alexander Keller, and Michael Garland.  "A decomposition for 
 in-place matrix transposition." Proceedings of the 19th ACM SIGPLAN symposium on 
 Principles and practice of parallel programming. ACM, 2014.
 * =====================================================================================
 */
int r(int i, int j, int m, int n) 
{
  int c = gcd(m,n);
  int b = n/c;
  int tmp = 1;
  int out;
  int tmp_in = static_cast<int>(tmp);
  
  tmp_in = i + floor(j/b);
  out = tmp_in % m;
  return out;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  d_prime
 *  Description:  This calculates the indicies for row shuffling.  This was developed 
 from the code algorithm described in:
 "Catanzaro, Bryan, Alexander Keller, and Michael Garland.  "A decomposition for 
 in-place matrix transposition." Proceedings of the 19th ACM SIGPLAN symposium on 
 Principles and practice of parallel programming. ACM, 2014.
 * =====================================================================================
 */
int d_prime(int i, int j, int m, int n) 
{
  int c = gcd(m,n);
  int b = n/c;
  int out;
  int tmp = 1;
  int tmp_in = static_cast<int>(tmp);
  
  tmp_in = i + floor(j/b);
  out = ((tmp_in % m) + j*m) % n;
          return out;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  s_prime
 *  Description:  This calculates the indices for column shuffling.  This was developed 
 from the code algorithm described in:
 "Catanzaro, Bryan, Alexander Keller, and Michael Garland.  "A decomposition for 
 in-place matrix transposition." Proceedings of the 19th ACM SIGPLAN symposium on 
 Principles and practice of parallel programming. ACM, 2014.
 * =====================================================================================
 */
int s_prime(int i, int j, int m, int n) 
{
  int c = gcd(m,n);
  int a = m/c;
  int out;
  int tmp = 1;
  int tmp_in = static_cast<int>(tmp);
  
  tmp_in = j + i*n - floor(i/a);
  out = tmp_in % m;
  return out;
}
    
/* 
 * ===  TEMPLATE  ======================================================================
 *         Name:  ipt
 *  Description:  This applys the logic for in-place transposition of matrix A. This 
 was developed from the code algorithm described in:
 "Catanzaro, Bryan, Alexander Keller, and Michael Garland.  "A decomposition for 
 in-place matrix transposition." Proceedings of the 19th ACM SIGPLAN symposium on 
 Principles and practice of parallel programming. ACM, 2014.
 * =====================================================================================
 */
template<typename T>
void ipt( MatrixAccessor<T> A, int m, int n) {
  
  // Allocate temporary space
  int *tmp = new int[max(m, n)];
  
  
  if(gcd(m, n) > 1) { // Removing conflicts
                      for(int j=0; j<n; j++){
                        //Gather with r(i,j)
                        for(int i=0; i<m; i++) {
                          tmp[i] = A[r(i,j,m,n)][j];
                        }
                        for(int i=0; i<m; i++) {
                          A[i][j] = tmp[i];
                        }
                      }
  }
  
  for(int i=0; i<m; i++) { // Row shuffles
                           // Scatter with d_prime(i,j)
                           for(int j=0; j<n; j++) {
                             tmp[d_prime(i,j,m,n)] = A[i][j];
                           } 
                           for(int j=0; j<n; j++) {
                             A[i][j]=tmp[j];
                           }
  }
  
  for(int j=0; j<n; j++) { // Column shuffles
                           // Gather with s_prime(i,j)
                           for(int i=0; i<m; i++) {
                             tmp[i] = A[s_prime(i,j,m,n)][j];
                           }
                           for(int i=0; i<m; i++) {
                             A[i][j] = tmp[i];
                           }
  }
  delete[] tmp;
}

/* 
 * ===  TEMPLATE  ======================================================================
 *         Name:  IPTBM
 *  Description:  This provides the logic for applying the ipt function, namely
 accessing the big.matrix elements. 
 * =====================================================================================
 */
template<typename in_CType, typename in_BMAccessorType>
void IPTBM(BigMatrix *pInMat)
{
  in_BMAccessorType inMat( *pInMat);
  
  // run in-place-transpose transpose
  ipt( inMat, pInMat->ncol(), pInMat->nrow() );
  
  return;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  CIPTMatrix
 *  Description:  This provides the SEXP for the IPTBM template.  This allows the
 function to be called in R with '.Call'. 
 * =====================================================================================
 */
extern "C"
{
  #define CALL_IPT_1(IN_CTYPE, IN_ACCESSOR) \
    switch(pInMat->matrix_type()) \
    { \
      case 1: \
        IPTBM<IN_CTYPE, IN_ACCESSOR<IN_CTYPE> >( \
          pInMat); \
        break; \
      case 2: \
        IPTBM<IN_CTYPE, IN_ACCESSOR<IN_CTYPE> >( \
          pInMat); \
        break; \
      case 4: \
        IPTBM<IN_CTYPE, IN_ACCESSOR<IN_CTYPE> >( \
          pInMat); \
        break; \
      case 8: \
        IPTBM<IN_CTYPE, IN_ACCESSOR<IN_CTYPE> >( \
          pInMat); \
        break; \
    }

  #define CALL_IPT_2(IN_ACCESSOR) \
    switch(pInMat->matrix_type()) \
    { \
      case 1: \
        CALL_IPT_1(char, IN_ACCESSOR) \
        break; \
      case 2: \
        CALL_IPT_1(short, IN_ACCESSOR) \
        break; \
      case 4: \
        CALL_IPT_1(int, IN_ACCESSOR) \
        break; \
      case 8: \
        CALL_IPT_1(double, IN_ACCESSOR) \
        break; \
    }
      
  SEXP CIPTMatrix(SEXP inAddr)
  {
    BigMatrix *pInMat = reinterpret_cast<BigMatrix*>(
      R_ExternalPtrAddr(inAddr));
    
    // Not sure if there is a better way to do these function calls
    if (pInMat->separated_columns()) {
      // Need method for separated_columns
      //CALL_IPT_2(SepMatrixAccessor)
    }
    else
    {
      CALL_IPT_2(MatrixAccessor)
    }

    return R_NilValue;
  }
}