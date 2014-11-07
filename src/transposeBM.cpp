#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

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

#include <Rcpp.h>

using namespace Rcpp;


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
