#include "QuasarCL.hpp"

using namespace quasarcl;


std::vector<const char*> kernelNames = {
	"matrix_log10",
	"matrix_transpose",
	"matrix_minus_scalar",
	"matrix_minus_matrix",
	"matrix_divide_matrix",
	"matrix_multiply_col_vector",
	"simple_mavg",
	"centered_mavg",
	"generateWavelengthsMatrix",
	"addSpectrum",
	"filterWithWavelengthWindows",
	"filterNonpositive",
	"filterZeros",
	"filterInfs",
	"convolve",
	"copyIfNotInf",
	"countIfNotInf",
	"reglin",
	"reglin_yax",
	"chisq",
	"integrate_trapz",
	"fix_reglin_results",
	"calc_cfun_dcfun",
	"calc_cw",
	"reduce_continuum_chisqs",
	"reduce_fe_chisqs",
	"fit_gaussian",
	"calc_gaussian",
	"calc_gaussian_chisq",
	"calc_gaussian_fwhm"
};



QuasarCL::QuasarCL() {
	this->initialized = false;
}


void QuasarCL::initialize(std::vector<std::string>& sourceFilenames) {
	try {
		std::vector<cl::Device> devices = getDevices();
		this->context =  cl::Context(devices);
		this->queue = cl::CommandQueue(context, devices[0], CL_QUEUE_PROFILING_ENABLE);
		
		std::vector<std::string> programSources = getSources(sourceFilenames);
		cl::Program::Sources sources;
		for (auto& source : programSources)
		{
			sources.push_back({ source.c_str(), source.length() });
		}
		this->program = cl::Program(context, sources);
		this->program.build(devices);
	
		for (const auto& kernelName : kernelNames) {
			cl::Kernel kernel = cl::Kernel(this->program, kernelName);
			addKernel(kernelName, kernel);
		}
		this->initialized = true;
	}
	catch(const std::exception& e) {
		Rcpp::stop("Initialization has failed: " + std::string(e.what()));
	}
}

bool QuasarCL::isInitialized() {
	return this->initialized;
}


cl::CommandQueue& QuasarCL::getQueue() {
	return this->queue;
}

cl::Context& QuasarCL::getContext() {
	return this->context;
}

cl::Program& QuasarCL::getProgram(){
	return this->program;
}

cl::Kernel& QuasarCL::getKernelByName(const char* name) {
	auto it = this->kernels.find(name);
	if (it == this->kernels.end()) {
		Rcpp::stop("Kernel missing: " + std::string(name));
	}
	return it->second;
}


void QuasarCL::addKernel(const char* name, cl::Kernel& kernel) {
	this->kernels.insert({ name, kernel });
}


std::vector<cl::Device> QuasarCL::getDevices() {
	std::vector<cl::Platform> platforms;
	cl::Platform::get(&platforms);

	if (platforms.empty())
	{
		Rcpp::stop("No OpenCL platform");
	}

	cl::Platform platform = platforms[0];

	std::vector<cl::Device> devices;
	platform.getDevices(CL_DEVICE_TYPE_GPU, &devices);
	if (devices.empty())
	{
		Rcpp::stop("No OpenCL GPU devices");
	}
	return devices;
}


