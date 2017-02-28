#include "RcppExtended.hpp"
#include <Rcpp.h>


namespace Rcpp {

	
	template <> cl_double2 as (SEXP vector_)
	{
		NumericVector vector(vector_);
		cl_double2 data;
		if (vector.size() != sizeof(data.s)/sizeof(*data.s)) {
			stop("No known conversion to cl_double2");
		}
		for(size_t i = 0; i < vector.size(); ++i)
		{
			data.s[i] = vector[i];
		}
		return data;
	};
	
	template <> cl_double4 as (SEXP vector_)
	{
		NumericVector vector(vector_);
		cl_double4 data;
		if (vector.size() != sizeof(data.s)/sizeof(*data.s)) {
			stop("No known conversion to cl_double4");
		}
		for(size_t i = 0; i < vector.size(); ++i)
		{
			data.s[i] = vector[i];
		}
		return data;
	};
	
	template <> cl_double8 as (SEXP vector_)
	{
		NumericVector vector(vector_);
		cl_double8 data;
		if (vector.size() != sizeof(data.s)/sizeof(*data.s)) {
			stop("No known conversion to cl_double8");
		}
		for(size_t i = 0; i < vector.size(); ++i)
		{
			data.s[i] = vector[i];
		}
		return data;
	};
	
	template <> SEXP wrap(const cl_double2 & obj)
	{
		NumericVector out(obj.s, obj.s+sizeof(obj.s)/sizeof(*obj.s));
		return out;
	};


	template <> SEXP wrap(const cl_double4 & obj)
	{
		NumericVector out(obj.s, obj.s+sizeof(obj.s)/sizeof(*obj.s));
		return out;
	};
	
	template <> SEXP wrap(const cl_double8 & obj)
	{
		NumericVector out(obj.s, obj.s+sizeof(obj.s)/sizeof(*obj.s));
		return out;
	};
}