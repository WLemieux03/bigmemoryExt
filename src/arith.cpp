#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cmath>

#include <boost/lexical_cast.hpp>

#include <Rcpp.h>

#include "bigmemory/BigMatrix.h"
#include "bigmemory/MatrixAccessor.hpp"
#include "bigmemory/util.h"
#include "bigmemory/isna.hpp"

#include <stdio.h>
#include <R.h>
#include <Rinternals.h>
// #include <Rdefines.h>
#include <stdlib.h>
#include <sys/types.h>


using namespace std;

/* 
 * ===  NEW-BM-TEMPLATES  ============================================================
 *  Description:  The following templates are for creating new big.matrix objects
 *  following the called mathematical function.
 * =====================================================================================
 */



/* 
 * ===  IN-PLACE-TEMPLATES  ============================================================
 *  Description:  The following templates are for the in-place modifications of
 *  big.matrix objects
 * =====================================================================================
 */

extern "C"
{

}