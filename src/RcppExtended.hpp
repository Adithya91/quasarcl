#ifndef CONVERT_H_
#define CONVERT_H_

#include <RcppCommon.h>
#define CL_STD   "-cl-std=CL1.1"
#define CL_USE_DEPRECATED_OPENCL_1_1_APIS
#include <CL/cl.h>
#undef CL_VERSION_1_2
#include <CL/cl.hpp>

namespace Rcpp 
{
		
		 template <> cl_double2 as (SEXP);
		 template <> cl_double4 as (SEXP);
		 template <> cl_double8 as (SEXP);
		
		 template <> SEXP wrap(const cl_double2 & obj);
		 template <> SEXP wrap(const cl_double4 & obj);
		 template <> SEXP wrap(const cl_double8 & obj);
		
}


#endif
