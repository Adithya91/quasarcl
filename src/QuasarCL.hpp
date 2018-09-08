#ifndef QUASAR_CL_H_
#define QUASAR_CL_H_

#define ASTRO_OBJ_SPEC_SIZE 4096

#ifndef __CL_ENABLE_EXCEPTIONS
#define __CL_ENABLE_EXCEPTIONS
#endif

#include <CL/cl.h>
#undef CL_VERSION_1_2
#include <CL/cl.hpp>
#include <map>

#include "utils.hpp"



namespace quasarcl {

	class QuasarCL {
	public:
		QuasarCL();
		//~QuasarCL();
		
		void initialize(std::vector<std::string>& sourcesPaths);
		bool isInitialized();
		
		cl::CommandQueue& getQueue();
		cl::Context& getContext();
		cl::Program& getProgram();
		cl::Kernel& getKernelByName(const char* name);
		
	private:
		void addKernel(const char* name, cl::Kernel& kernel);
		std::vector<cl::Device> getDevices();

		cl::CommandQueue queue;
		cl::Context context;
		cl::Program program;
		std::map<const char*, cl::Kernel> kernels;	
		bool initialized;
	};

}

#endif
