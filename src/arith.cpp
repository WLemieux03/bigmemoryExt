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

/* 
 * ===  NEW-BM-TEMPLATES  ============================================================
 *  Description:  The following templates are for creating new big.matrix objects
 following the called mathematical function.
 * =====================================================================================
 */
template<typename in_CType, typename in_BMAccessorType, 
  typename out_CType, typename out_BMAccessorType>
void IntSubtractBM(BigMatrix *pInMat, BigMatrix *pOutMat, int value_)
{
  in_BMAccessorType inMat( *pInMat );
  out_BMAccessorType outMat( *pOutMat );
  
  for(index_type i=0; i < pInMat->nrow(); i++) {
    for(index_type j=0; j < pInMat->ncol(); j++){
      outMat[j][i] = value_ - inMat[j][i];
    }
  }
  
  return;
}

template<typename in_CType, typename in_BMAccessorType,
  typename out_CType, typename out_BMAccessorType>
void powBM(BigMatrix *pInMat, BigMatrix *pOutMat, double value)
{
  in_BMAccessorType inMat( *pInMat );
  out_BMAccessorType outMat( *pOutMat );
  
  for(index_type i=0; i < pInMat->nrow(); i++) {
    for(index_type j=0; j < pInMat->ncol(); j++){
      outMat[j][i] = pow(inMat[j][i], value);
    }
  }
  
  return;
}

template<typename in_CType, typename in_BMAccessorType,
  typename out_CType, typename out_BMAccessorType>
void expBM(BigMatrix *pInMat, BigMatrix *pOutMat)
{
  in_BMAccessorType inMat( *pInMat );
  out_BMAccessorType outMat( *pOutMat );
  
  for(index_type i=0; i < pInMat->nrow(); i++) {
    for(index_type j=0; j < pInMat->ncol(); j++){
      outMat[j][i] = exp(inMat[j][i]);
    }
  }
  
  return;
}

template<typename in_CType, typename in_BMAccessorType,
  typename out_CType, typename out_BMAccessorType>
void log10BM(BigMatrix *pInMat, BigMatrix *pOutMat)
{
  in_BMAccessorType inMat( *pInMat );
  out_BMAccessorType outMat( *pOutMat );
  
  for(index_type i=0; i < pInMat->nrow(); i++) {
    for(index_type j=0; j < pInMat->ncol(); j++){
      outMat[j][i] = log10(inMat[j][i]);
    }
  }
  
  return;
}

template<typename in_CType, typename in_BMAccessorType,
  typename out_CType, typename out_BMAccessorType>
void logBaseBM(BigMatrix *pInMat, BigMatrix *pOutMat, double base)
{
  in_BMAccessorType inMat( *pInMat );
  out_BMAccessorType outMat( *pOutMat );
    
  for(index_type i=0; i < pInMat->nrow(); i++) {
    for(index_type j=0; j < pInMat->ncol(); j++){
      outMat[j][i] = log10(inMat[j][i])/log10(base);
    }
  }
  
  return;
}

template<typename in_CType, typename in_BMAccessorType,
  typename out_CType, typename out_BMAccessorType>
void tanhBM(BigMatrix *pInMat, BigMatrix *pOutMat)
{
  in_BMAccessorType inMat( *pInMat );
  out_BMAccessorType outMat( *pOutMat );
  
  for(index_type i=0; i < pInMat->nrow(); i++) {
    for(index_type j=0; j < pInMat->ncol(); j++){
      outMat[j][i] = tanh(inMat[j][i]);
    }
  }
  
  return;
}

template<typename in_CType, typename in_BMAccessorType,
  typename out_CType, typename out_BMAccessorType>
void coshBM(BigMatrix *pInMat, BigMatrix *pOutMat)
{
  in_BMAccessorType inMat( *pInMat );
  out_BMAccessorType outMat( *pOutMat );
  
  for(index_type i=0; i < pInMat->nrow(); i++) {
    for(index_type j=0; j < pInMat->ncol(); j++){
      outMat[j][i] = cosh(inMat[j][i]);
    }
  }
  
  return;
}

template<typename in_CType, typename in_BMAccessorType,
  typename out_CType, typename out_BMAccessorType>
void sinhBM(BigMatrix *pInMat, BigMatrix *pOutMat)
{
  in_BMAccessorType inMat( *pInMat );
  out_BMAccessorType outMat( *pOutMat );
  
  for(index_type i=0; i < pInMat->nrow(); i++) {
    for(index_type j=0; j < pInMat->ncol(); j++){
      outMat[j][i] = sinh(inMat[j][i]);
    }
  }
  
  return;
}

/* 
 * ===  IN-PLACE-TEMPLATES  ============================================================
 *  Description:  The following templates are for the in-place modifications of
 big.matrix objects
 * =====================================================================================
 */
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

template<typename in_CType, typename in_BMAccessorType>
void expBMIP(BigMatrix *pInMat)
{
  in_BMAccessorType inMat( *pInMat );
  
  for(index_type i=0; i<pInMat->nrow(); i++ ){
    for(index_type j=0; j<pInMat->ncol(); j++){
      inMat[j][i] = exp(inMat[j][i]);
    }
  }
  
  return;
}

template<typename in_CType, typename in_BMAccessorType>
void logBMIP(BigMatrix *pInMat, double base)
{
  in_BMAccessorType inMat( *pInMat );
  
  for(index_type i=0; i<pInMat->nrow(); i++ ){
    for(index_type j=0; j<pInMat->ncol(); j++){
      inMat[j][i] = log10(inMat[j][i])/log10(base);
    }
  }
  
  return;
}

template<typename in_CType, typename in_BMAccessorType>
void tanhBMIP(BigMatrix *pInMat)
{
  in_BMAccessorType inMat( *pInMat );
  
  for(index_type i=0; i<pInMat->nrow(); i++ ){
    for(index_type j=0; j<pInMat->ncol(); j++){
      inMat[j][i] = tanh(inMat[j][i]);
    }
  }
  
  return;
}

template<typename in_CType, typename in_BMAccessorType>
void coshBMIP(BigMatrix *pInMat)
{
  in_BMAccessorType inMat( *pInMat );
  
  for(index_type i=0; i<pInMat->nrow(); i++ ){
    for(index_type j=0; j<pInMat->ncol(); j++){
      inMat[j][i] = cosh(inMat[j][i]);
    }
  }
  
  return;
}

template<typename in_CType, typename in_BMAccessorType>
void sinhBMIP(BigMatrix *pInMat)
{
  in_BMAccessorType inMat( *pInMat );
  
  for(index_type i=0; i<pInMat->nrow(); i++ ){
    for(index_type j=0; j<pInMat->ncol(); j++){
      inMat[j][i] = sinh(inMat[j][i]);
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
          pInMat, pOutMat, value); \
        break; \
      case 8: \
        IntSubtractBM<IN_CTYPE, IN_ACCESSOR<IN_CTYPE>, double, OUT_ACCESSOR<double> >( \
          pInMat, pOutMat, value); \
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
      
  SEXP CsubtIntBM(SEXP inAddr, SEXP outAddr, SEXP typecast_warning, SEXP value_)
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
  
  #define CALL_powBM_2(IN_CTYPE, IN_ACCESSOR, OUT_ACCESSOR, value) \
    switch(pOutMat->matrix_type()) \
    { \
      case 4: \
        powBM<IN_CTYPE, IN_ACCESSOR<IN_CTYPE>, int, OUT_ACCESSOR<int> >( \
          pInMat, pOutMat, value); \
        break; \
      case 8: \
        powBM<IN_CTYPE, IN_ACCESSOR<IN_CTYPE>, double, OUT_ACCESSOR<double> >( \
          pInMat, pOutMat, value); \
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
      
  SEXP CpowBM(SEXP inAddr, SEXP outAddr, SEXP value, SEXP typecast_warning)
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
    
    double val = as<double>(value);
    
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
    int val = as<double>(value);
    
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
  
  #define CALL_expBM_2(IN_CTYPE, IN_ACCESSOR, OUT_ACCESSOR) \
    switch(pOutMat->matrix_type()) \
    { \
      case 4: \
        expBM<IN_CTYPE, IN_ACCESSOR<IN_CTYPE>, int, OUT_ACCESSOR<int> >( \
          pInMat, pOutMat); \
        break; \
      case 8: \
        expBM<IN_CTYPE, IN_ACCESSOR<IN_CTYPE>, double, OUT_ACCESSOR<double> >( \
          pInMat, pOutMat); \
        break; \
    }

  #define CALL_expBM_1(IN_ACCESSOR, OUT_ACCESSOR) \
    switch(pInMat->matrix_type()) \
    { \
      case 4: \
        CALL_expBM_2(int, IN_ACCESSOR, OUT_ACCESSOR) \
        break; \
      case 8: \
        CALL_expBM_2(double, IN_ACCESSOR, OUT_ACCESSOR) \
        break; \
    }
      
  SEXP CexpBM(SEXP inAddr, SEXP outAddr, SEXP typecast_warning)
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
      CALL_expBM_1(SepMatrixAccessor, SepMatrixAccessor)
    }
    else if(pInMat->separated_columns() && !(pOutMat->separated_columns()))
    {
      CALL_expBM_1(SepMatrixAccessor, MatrixAccessor)
    }
    else if(!(pInMat->separated_columns()) && pOutMat->separated_columns())
    {
      CALL_expBM_1(MatrixAccessor, SepMatrixAccessor)
    }
    else
    {
      CALL_expBM_1(MatrixAccessor, MatrixAccessor)
    }

    return R_NilValue;
  }
  
  #define CALL_expBMIP_2(IN_CTYPE, IN_ACCESSOR) \
    switch(pInMat->matrix_type()) \
    { \
      case 4: \
        expBMIP<IN_CTYPE, IN_ACCESSOR<IN_CTYPE> >( \
          pInMat); \
        break; \
      case 8: \
        expBMIP<IN_CTYPE, IN_ACCESSOR<IN_CTYPE> >( \
          pInMat); \
        break; \
    }

  #define CALL_expBMIP_1(IN_ACCESSOR) \
    switch(pInMat->matrix_type()) \
    { \
      case 4: \
        CALL_expBMIP_2(int, IN_ACCESSOR) \
        break; \
      case 8: \
        CALL_expBMIP_2(double, IN_ACCESSOR) \
        break; \
    }
      
  SEXP CexpBMIP(SEXP inAddr, SEXP typecast_warning)
  {
    BigMatrix *pInMat = reinterpret_cast<BigMatrix*>(
      R_ExternalPtrAddr(inAddr));
    
    // Not sure if there is a better way to do these function calls
    if (pInMat->separated_columns() ) {
      CALL_expBMIP_1(SepMatrixAccessor)
    }
    else
    {
      CALL_expBMIP_1(MatrixAccessor)
    }

    return R_NilValue;
  }
  
  #define CALL_log10BM_2(IN_CTYPE, IN_ACCESSOR, OUT_ACCESSOR) \
    switch(pOutMat->matrix_type()) \
    { \
      case 4: \
        log10BM<IN_CTYPE, IN_ACCESSOR<IN_CTYPE>, int, OUT_ACCESSOR<int> >( \
          pInMat, pOutMat); \
        break; \
      case 8: \
        log10BM<IN_CTYPE, IN_ACCESSOR<IN_CTYPE>, double, OUT_ACCESSOR<double> >( \
          pInMat, pOutMat); \
        break; \
    }

  #define CALL_log10BM_1(IN_ACCESSOR, OUT_ACCESSOR) \
    switch(pInMat->matrix_type()) \
    { \
      case 4: \
        CALL_log10BM_2(int, IN_ACCESSOR, OUT_ACCESSOR) \
        break; \
      case 8: \
        CALL_log10BM_2(double, IN_ACCESSOR, OUT_ACCESSOR) \
        break; \
    }
      
  SEXP Clog10BM(SEXP inAddr, SEXP outAddr, SEXP typecast_warning)
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
      CALL_log10BM_1(SepMatrixAccessor, SepMatrixAccessor)
    }
    else if(pInMat->separated_columns() && !(pOutMat->separated_columns()))
    {
      CALL_log10BM_1(SepMatrixAccessor, MatrixAccessor)
    }
    else if(!(pInMat->separated_columns()) && pOutMat->separated_columns())
    {
      CALL_log10BM_1(MatrixAccessor, SepMatrixAccessor)
    }
    else
    {
      CALL_log10BM_1(MatrixAccessor, MatrixAccessor)
    }

    return R_NilValue;
  }
  
  #define CALL_logBaseBM_2(IN_CTYPE, IN_ACCESSOR, OUT_ACCESSOR, base) \
    switch(pOutMat->matrix_type()) \
    { \
      case 4: \
        logBaseBM<IN_CTYPE, IN_ACCESSOR<IN_CTYPE>, int, OUT_ACCESSOR<int> >( \
          pInMat, pOutMat, base); \
        break; \
      case 8: \
        logBaseBM<IN_CTYPE, IN_ACCESSOR<IN_CTYPE>, double, OUT_ACCESSOR<double> >( \
          pInMat, pOutMat, base); \
        break; \
    }

  #define CALL_logBaseBM_1(IN_ACCESSOR, OUT_ACCESSOR, base) \
    switch(pInMat->matrix_type()) \
    { \
      case 4: \
        CALL_logBaseBM_2(int, IN_ACCESSOR, OUT_ACCESSOR, base) \
        break; \
      case 8: \
        CALL_logBaseBM_2(double, IN_ACCESSOR, OUT_ACCESSOR, base) \
        break; \
    }
      
  SEXP ClogBaseBM(SEXP inAddr, SEXP outAddr, SEXP _base, SEXP typecast_warning)
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
    
    double base = as<double>(_base);
    
    // Not sure if there is a better way to do these function calls
    if (pInMat->separated_columns() && pOutMat->separated_columns()) {
      CALL_logBaseBM_1(SepMatrixAccessor, SepMatrixAccessor, base)
    }
    else if(pInMat->separated_columns() && !(pOutMat->separated_columns()))
    {
      CALL_logBaseBM_1(SepMatrixAccessor, MatrixAccessor, base)
    }
    else if(!(pInMat->separated_columns()) && pOutMat->separated_columns())
    {
      CALL_logBaseBM_1(MatrixAccessor, SepMatrixAccessor, base)
    }
    else
    {
      CALL_logBaseBM_1(MatrixAccessor, MatrixAccessor, base)
    }

    return R_NilValue;
  }
  
  #define CALL_logBMIP_2(IN_CTYPE, IN_ACCESSOR, base) \
    switch(pInMat->matrix_type()) \
    { \
      case 4: \
        logBMIP<IN_CTYPE, IN_ACCESSOR<IN_CTYPE> >( \
          pInMat, base); \
        break; \
      case 8: \
        logBMIP<IN_CTYPE, IN_ACCESSOR<IN_CTYPE> >( \
          pInMat, base); \
        break; \
    }

  #define CALL_logBMIP_1(IN_ACCESSOR, base) \
    switch(pInMat->matrix_type()) \
    { \
      case 4: \
        CALL_logBMIP_2(int, IN_ACCESSOR, base) \
        break; \
      case 8: \
        CALL_logBMIP_2(double, IN_ACCESSOR, base) \
        break; \
    }
      
  SEXP ClogBMIP(SEXP inAddr, SEXP _base, SEXP typecast_warning)
  {
    BigMatrix *pInMat = reinterpret_cast<BigMatrix*>(
      R_ExternalPtrAddr(inAddr));
    
    double base = as<double>(_base);
    
    // Not sure if there is a better way to do these function calls
    if (pInMat->separated_columns() ) {
      CALL_logBMIP_1(SepMatrixAccessor, base)
    }
    else
    {
      CALL_logBMIP_1(MatrixAccessor, base)
    }

    return R_NilValue;
  }
  
  #define CALL_tanhBM_2(IN_CTYPE, IN_ACCESSOR, OUT_ACCESSOR) \
    switch(pOutMat->matrix_type()) \
    { \
      case 4: \
        tanhBM<IN_CTYPE, IN_ACCESSOR<IN_CTYPE>, int, OUT_ACCESSOR<int> >( \
          pInMat, pOutMat); \
        break; \
      case 8: \
        tanhBM<IN_CTYPE, IN_ACCESSOR<IN_CTYPE>, double, OUT_ACCESSOR<double> >( \
          pInMat, pOutMat); \
        break; \
    }

  #define CALL_tanhBM_1(IN_ACCESSOR, OUT_ACCESSOR) \
    switch(pInMat->matrix_type()) \
    { \
      case 4: \
        CALL_tanhBM_2(int, IN_ACCESSOR, OUT_ACCESSOR) \
        break; \
      case 8: \
        CALL_tanhBM_2(double, IN_ACCESSOR, OUT_ACCESSOR) \
        break; \
    }
      
  SEXP CtanhBM(SEXP inAddr, SEXP outAddr, SEXP typecast_warning)
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
      CALL_tanhBM_1(SepMatrixAccessor, SepMatrixAccessor)
    }
    else if(pInMat->separated_columns() && !(pOutMat->separated_columns()))
    {
      CALL_tanhBM_1(SepMatrixAccessor, MatrixAccessor)
    }
    else if(!(pInMat->separated_columns()) && pOutMat->separated_columns())
    {
      CALL_tanhBM_1(MatrixAccessor, SepMatrixAccessor)
    }
    else
    {
      CALL_tanhBM_1(MatrixAccessor, MatrixAccessor)
    }

    return R_NilValue;
  }
  
  #define CALL_tanhBMIP_2(IN_CTYPE, IN_ACCESSOR) \
    switch(pInMat->matrix_type()) \
    { \
      case 4: \
        tanhBMIP<IN_CTYPE, IN_ACCESSOR<IN_CTYPE> >( \
          pInMat); \
        break; \
      case 8: \
        tanhBMIP<IN_CTYPE, IN_ACCESSOR<IN_CTYPE> >( \
          pInMat); \
        break; \
    }

  #define CALL_tanhBMIP_1(IN_ACCESSOR) \
    switch(pInMat->matrix_type()) \
    { \
      case 4: \
        CALL_tanhBMIP_2(int, IN_ACCESSOR) \
        break; \
      case 8: \
        CALL_tanhBMIP_2(double, IN_ACCESSOR) \
        break; \
    }
      
  SEXP CtanhBMIP(SEXP inAddr, SEXP typecast_warning)
  {
    BigMatrix *pInMat = reinterpret_cast<BigMatrix*>(
      R_ExternalPtrAddr(inAddr));
    
    // Not sure if there is a better way to do these function calls
    if (pInMat->separated_columns() ) {
      CALL_tanhBMIP_1(SepMatrixAccessor)
    }
    else
    {
      CALL_tanhBMIP_1(MatrixAccessor)
    }

    return R_NilValue;
  }
  
  #define CALL_coshBM_2(IN_CTYPE, IN_ACCESSOR, OUT_ACCESSOR) \
    switch(pOutMat->matrix_type()) \
    { \
      case 4: \
        coshBM<IN_CTYPE, IN_ACCESSOR<IN_CTYPE>, int, OUT_ACCESSOR<int> >( \
          pInMat, pOutMat); \
        break; \
      case 8: \
        coshBM<IN_CTYPE, IN_ACCESSOR<IN_CTYPE>, double, OUT_ACCESSOR<double> >( \
          pInMat, pOutMat); \
        break; \
    }

  #define CALL_coshBM_1(IN_ACCESSOR, OUT_ACCESSOR) \
    switch(pInMat->matrix_type()) \
    { \
      case 4: \
        CALL_coshBM_2(int, IN_ACCESSOR, OUT_ACCESSOR) \
        break; \
      case 8: \
        CALL_coshBM_2(double, IN_ACCESSOR, OUT_ACCESSOR) \
        break; \
    }
      
  SEXP CcoshBM(SEXP inAddr, SEXP outAddr, SEXP typecast_warning)
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
      CALL_coshBM_1(SepMatrixAccessor, SepMatrixAccessor)
    }
    else if(pInMat->separated_columns() && !(pOutMat->separated_columns()))
    {
      CALL_coshBM_1(SepMatrixAccessor, MatrixAccessor)
    }
    else if(!(pInMat->separated_columns()) && pOutMat->separated_columns())
    {
      CALL_coshBM_1(MatrixAccessor, SepMatrixAccessor)
    }
    else
    {
      CALL_coshBM_1(MatrixAccessor, MatrixAccessor)
    }

    return R_NilValue;
  }
  
  #define CALL_coshBMIP_2(IN_CTYPE, IN_ACCESSOR) \
    switch(pInMat->matrix_type()) \
    { \
      case 4: \
        coshBMIP<IN_CTYPE, IN_ACCESSOR<IN_CTYPE> >( \
          pInMat); \
        break; \
      case 8: \
        coshBMIP<IN_CTYPE, IN_ACCESSOR<IN_CTYPE> >( \
          pInMat); \
        break; \
    }

  #define CALL_coshBMIP_1(IN_ACCESSOR) \
    switch(pInMat->matrix_type()) \
    { \
      case 4: \
        CALL_coshBMIP_2(int, IN_ACCESSOR) \
        break; \
      case 8: \
        CALL_coshBMIP_2(double, IN_ACCESSOR) \
        break; \
    }
      
  SEXP CcoshBMIP(SEXP inAddr, SEXP typecast_warning)
  {
    BigMatrix *pInMat = reinterpret_cast<BigMatrix*>(
      R_ExternalPtrAddr(inAddr));
    
    // Not sure if there is a better way to do these function calls
    if (pInMat->separated_columns() ) {
      CALL_coshBMIP_1(SepMatrixAccessor)
    }
    else
    {
      CALL_coshBMIP_1(MatrixAccessor)
    }

    return R_NilValue;
  }
  
  #define CALL_sinhBM_2(IN_CTYPE, IN_ACCESSOR, OUT_ACCESSOR) \
    switch(pOutMat->matrix_type()) \
    { \
      case 4: \
        sinhBM<IN_CTYPE, IN_ACCESSOR<IN_CTYPE>, int, OUT_ACCESSOR<int> >( \
          pInMat, pOutMat); \
        break; \
      case 8: \
        sinhBM<IN_CTYPE, IN_ACCESSOR<IN_CTYPE>, double, OUT_ACCESSOR<double> >( \
          pInMat, pOutMat); \
        break; \
    }

  #define CALL_sinhBM_1(IN_ACCESSOR, OUT_ACCESSOR) \
    switch(pInMat->matrix_type()) \
    { \
      case 4: \
        CALL_sinhBM_2(int, IN_ACCESSOR, OUT_ACCESSOR) \
        break; \
      case 8: \
        CALL_sinhBM_2(double, IN_ACCESSOR, OUT_ACCESSOR) \
        break; \
    }
      
  SEXP CsinhBM(SEXP inAddr, SEXP outAddr, SEXP typecast_warning)
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
      CALL_sinhBM_1(SepMatrixAccessor, SepMatrixAccessor)
    }
    else if(pInMat->separated_columns() && !(pOutMat->separated_columns()))
    {
      CALL_sinhBM_1(SepMatrixAccessor, MatrixAccessor)
    }
    else if(!(pInMat->separated_columns()) && pOutMat->separated_columns())
    {
      CALL_sinhBM_1(MatrixAccessor, SepMatrixAccessor)
    }
    else
    {
      CALL_sinhBM_1(MatrixAccessor, MatrixAccessor)
    }

    return R_NilValue;
  }
  
  #define CALL_sinhBMIP_2(IN_CTYPE, IN_ACCESSOR) \
    switch(pInMat->matrix_type()) \
    { \
      case 4: \
        sinhBMIP<IN_CTYPE, IN_ACCESSOR<IN_CTYPE> >( \
          pInMat); \
        break; \
      case 8: \
        sinhBMIP<IN_CTYPE, IN_ACCESSOR<IN_CTYPE> >( \
          pInMat); \
        break; \
    }

  #define CALL_sinhBMIP_1(IN_ACCESSOR) \
    switch(pInMat->matrix_type()) \
    { \
      case 4: \
        CALL_sinhBMIP_2(int, IN_ACCESSOR) \
        break; \
      case 8: \
        CALL_sinhBMIP_2(double, IN_ACCESSOR) \
        break; \
    }
      
  SEXP CsinhBMIP(SEXP inAddr, SEXP typecast_warning)
  {
    BigMatrix *pInMat = reinterpret_cast<BigMatrix*>(
      R_ExternalPtrAddr(inAddr));
    
    // Not sure if there is a better way to do these function calls
    if (pInMat->separated_columns() ) {
      CALL_sinhBMIP_1(SepMatrixAccessor)
    }
    else
    {
      CALL_sinhBMIP_1(MatrixAccessor)
    }

    return R_NilValue;
  }
}