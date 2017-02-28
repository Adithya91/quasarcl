#include "QuasarCL.hpp"

	
//[[Rcpp::export]]
SEXP cppInitialize(SEXP sources_) {
	std::vector<std::string> sources = Rcpp::as<std::vector<std::string> >(sources_);
	quasarcl::QuasarCL* quasarCL = new quasarcl::QuasarCL();
	quasarCL->initialize(sources);
	Rcpp::XPtr<quasarcl::QuasarCL> ptr(quasarCL);
	return ptr;
}



//[[Rcpp::export]]
SEXP cppIsInitialized(SEXP quasarclPtr_) {
	Rcpp::XPtr<quasarcl::QuasarCL> quasarclPtr(quasarclPtr_);
	return Rcpp::wrap(quasarclPtr->isInitialized());
}