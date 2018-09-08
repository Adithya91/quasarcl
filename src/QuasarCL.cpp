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
		cl_int err = CL_SUCCESS;
		this->program = cl::Program(context, sources, &err);
		std::cout << "Error: " << err << std::endl;
		try{
			this->program.build(devices,CL_STD);
		}
		catch(const cl::Error &e) {
		 if (e.err() == CL_BUILD_PROGRAM_FAILURE)
  		 {
    			for (cl::Device dev : devices)
    			{
      			// Check the build status
      			cl_build_status status = program.getBuildInfo<CL_PROGRAM_BUILD_STATUS>(dev);
      			if (status != CL_BUILD_ERROR)
        			continue;

      			// Get the build log
      			std::string name     = dev.getInfo<CL_DEVICE_NAME>();
      			std::string buildlog = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(dev);
      			std::cerr << "Build log for " << name << ":" << std::endl << buildlog << std::endl;
  			}
		 }
  		 else
  		 {
    		   	throw e;
  		 }		
		}

		for (const auto& kernelName : kernelNames) {
			cl::Kernel kernel = cl::Kernel(this->program, kernelName);
			addKernel(kernelName, kernel);
		}
		this->initialized = true;
	}
	//catch(const std::exception& e) {
	catch(const cl::Error &e) {
		Rcpp::stop("Initialization has failed: " + std::string(e.what()) + std::string(": ") + std::to_string(e.err()));
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

	int device_id = 0;

	std::cout << "Number of Platforms: " << platforms.size() << std::endl;

	if (platforms.empty())
	{
		Rcpp::stop("No OpenCL platform");
	}

	cl::Platform platform = platforms[0];

	std::vector<cl::Device> devices;
	platform.getDevices(CL_DEVICE_TYPE_GPU, &devices);
	for(std::vector<cl::Device>::iterator it2 = devices.begin(); it2 != devices.end(); ++it2){
            cl::Device device(*it2);

            std::cout << "\tDevice " << device_id++ << ": " << std::endl;
            std::cout << "\t\tDevice Name: " << device.getInfo<CL_DEVICE_NAME>() << std::endl;
            std::cout << "\t\tDevice Type: " << device.getInfo<CL_DEVICE_TYPE>();
            std::cout << " (GPU: " << CL_DEVICE_TYPE_GPU << ", CPU: " << CL_DEVICE_TYPE_CPU << ")" << std::endl;
            std::cout << "\t\tDevice Vendor: " << device.getInfo<CL_DEVICE_VENDOR>() << std::endl;
            std::cout << "\t\tDevice Max Compute Units: " << device.getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>() << std::endl;
            std::cout << "\t\tDevice Global Memory: " << device.getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>() << std::endl;
            std::cout << "\t\tDevice Max Clock Frequency: " << device.getInfo<CL_DEVICE_MAX_CLOCK_FREQUENCY>() << std::endl;
            std::cout << "\t\tDevice Max Allocateable Memory: " << device.getInfo<CL_DEVICE_MAX_MEM_ALLOC_SIZE>() << std::endl;
            std::cout << "\t\tDevice Local Memory: " << device.getInfo<CL_DEVICE_LOCAL_MEM_SIZE>() << std::endl;
            std::cout << "\t\tDevice Available: " << device.getInfo< CL_DEVICE_AVAILABLE>() << std::endl;
        }
        std::cout<< std::endl;

	if (devices.empty())
	{
		Rcpp::stop("No OpenCL GPU devices");
	}
	return devices;
}


