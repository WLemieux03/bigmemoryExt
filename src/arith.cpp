#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cmath>

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
using namespace std;


template<typename in_CType, typename in_BMAccessorType, 
  typename out_CType, typename out_BMAccessorType>
void IntSubtractBM(BigMatrix *pInMat, BigMatrix *pOutMat, SEXP rowInds, SEXP colInds, int value_)
{
  in_BMAccessorType inMat( *pInMat );
  out_BMAccessorType outMat( *pOutMat );
  
  index_type nRows = GET_LENGTH(rowInds);
  index_type nCols = GET_LENGTH(colInds);
  
  if (nRows != pOutMat->nrow())
    Rf_error("length of row indices does not equal # of rows in new matrix");
  if (nCols != pOutMat->ncol())
    Rf_error("length of col indices does not equal # of cols in new matrix");
  
  index_type i = 0;
  index_type j = 0;
  
  for(i=0; i < nRows; i++) {
    for(j=0; j < nCols; j++){
      outMat[j][i] = value_ - inMat[j][i];
    }
  }
  
  return;
}

extern "C"
{
  #define CALL_intSBM_2(IN_CTYPE, IN_ACCESSOR, OUT_ACCESSOR, value) \
    switch(pOutMat->matrix_type()) \
    { \
      case 4: \
        IntSubtractBM<IN_CTYPE, IN_ACCESSOR<IN_CTYPE>, int, OUT_ACCESSOR<int> >( \
          pInMat, pOutMat, rowInds, colInds, value); \
        break; \
      case 8: \
        IntSubtractBM<IN_CTYPE, IN_ACCESSOR<IN_CTYPE>, double, OUT_ACCESSOR<double> >( \
          pInMat, pOutMat, rowInds, colInds, value); \
        break; \
    }

  #define CALL_intSBM_1(IN_ACCESSOR, OUT_ACCESSOR, value) \
    switch(pInMat->matrix_type()) \
    { \
      case 4: \
        CALL_intSBM_2(int, IN_ACCESSOR, OUT_ACCESSOR, value) \
        break; \
      case 8: \
        CALL_intSBM_2(double, IN_ACCESSOR, OUT_ACCESSOR, value) \
        break; \
    }
      
  SEXP CsubtIntBM(SEXP inAddr, SEXP outAddr, SEXP rowInds, SEXP colInds, 
    SEXP typecast_warning, SEXP value_)
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
    
    int value = as<int>(value_);
    
    // Not sure if there is a better way to do these function calls
    if (pInMat->separated_columns() && pOutMat->separated_columns()) {
      CALL_intSBM_1(SepMatrixAccessor, SepMatrixAccessor, value)
    }
    else if(pInMat->separated_columns() && !(pOutMat->separated_columns()))
    {
      CALL_intSBM_1(SepMatrixAccessor, MatrixAccessor, value)
    }
    else if(!(pInMat->separated_columns()) && pOutMat->separated_columns())
    {
      CALL_intSBM_1(MatrixAccessor, SepMatrixAccessor, value)
    }
    else
    {
      CALL_intSBM_1(MatrixAccessor, MatrixAccessor, value)
    }

    return R_NilValue;
  }
}

template<typename in_CType, typename in_BMAccessorType,
  typename out_CType, typename out_BMAccessorType>
void powBM(BigMatrix *pInMat, BigMatrix *pOutMat, SEXP rowInds, SEXP colInds, int value)
{
  in_BMAccessorType inMat( *pInMat );
  out_BMAccessorType outMat( *pOutMat );
  
  index_type nRows = GET_LENGTH(rowInds);
  index_type nCols = GET_LENGTH(colInds);
  
  for(index_type i=0; i < nRows; i++) {
    for(index_type j=0; j < nCols; j++){
      outMat[j][i] = pow(inMat[j][i], value);
    }
  }
  
  return;
}

template<typename in_CType, typename in_BMAccessorType>
void powBMIP(BigMatrix *pInMat, int value)
{
  in_BMAccessorType inMat( *pInMat );
  
  for(index_type i=0; i<pInMat->nrow(); i++ ){
    for(index_type j=0; j<pInMat->ncol(); j++){
      inMat[j][i] = pow(inMat[j][i], value);
    }
  }
  
  return;
}

extern "C"
{
  #define CALL_powBM_2(IN_CTYPE, IN_ACCESSOR, OUT_ACCESSOR, value) \
    switch(pOutMat->matrix_type()) \
    { \
      case 4: \
        powBM<IN_CTYPE, IN_ACCESSOR<IN_CTYPE>, int, OUT_ACCESSOR<int> >( \
          pInMat, pOutMat, rowInds, colInds, value); \
        break; \
      case 8: \
        powBM<IN_CTYPE, IN_ACCESSOR<IN_CTYPE>, double, OUT_ACCESSOR<double> >( \
          pInMat, pOutMat, rowInds, colInds, value); \
        break; \
    }

  #define CALL_powBM_1(IN_ACCESSOR, OUT_ACCESSOR, value) \
    switch(pInMat->matrix_type()) \
    { \
      case 4: \
        CALL_powBM_2(int, IN_ACCESSOR, OUT_ACCESSOR, value) \
        break; \
      case 8: \
        CALL_powBM_2(double, IN_ACCESSOR, OUT_ACCESSOR, value) \
        break; \
    }
      
  SEXP CpowBM(SEXP inAddr, SEXP outAddr, SEXP rowInds, SEXP colInds,
    SEXP typecast_warning, SEXP value)
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
    
    int val = as<int>(value);
    
    // Not sure if there is a better way to do these function calls
    if (pInMat->separated_columns() && pOutMat->separated_columns()) {
      CALL_powBM_1(SepMatrixAccessor, SepMatrixAccessor, val)
    }
    else if(pInMat->separated_columns() && !(pOutMat->separated_columns()))
    {
      CALL_powBM_1(SepMatrixAccessor, MatrixAccessor, val)
    }
    else if(!(pInMat->separated_columns()) && pOutMat->separated_columns())
    {
      CALL_powBM_1(MatrixAccessor, SepMatrixAccessor, val)
    }
    else
    {
      CALL_powBM_1(MatrixAccessor, MatrixAccessor, val)
    }

    return R_NilValue;
  }
  
  #define CALL_powBMIP_2(IN_CTYPE, IN_ACCESSOR, value) \
    switch(pInMat->matrix_type()) \
    { \
      case 4: \
        powBMIP<IN_CTYPE, IN_ACCESSOR<IN_CTYPE> >( \
          pInMat, value); \
        break; \
      case 8: \
        powBMIP<IN_CTYPE, IN_ACCESSOR<IN_CTYPE> >( \
          pInMat, value); \
        break; \
    }

  #define CALL_powBMIP_1(IN_ACCESSOR, value) \
    switch(pInMat->matrix_type()) \
    { \
      case 4: \
        CALL_powBMIP_2(int, IN_ACCESSOR, value) \
        break; \
      case 8: \
        CALL_powBMIP_2(double, IN_ACCESSOR, value) \
        break; \
    }
      
  SEXP CpowBMIP(SEXP inAddr, SEXP value, SEXP typecast_warning)
  {
    BigMatrix *pInMat = reinterpret_cast<BigMatrix*>(
      R_ExternalPtrAddr(inAddr));
    
    // Convert R numeric to C++ integer
    int val = as<int>(value);
    
    // Not sure if there is a better way to do these function calls
    if (pInMat->separated_columns() ) {
      CALL_powBMIP_1(SepMatrixAccessor, val)
    }
    else
    {
      CALL_powBMIP_1(MatrixAccessor, val)
    }

    return R_NilValue;
  }
}