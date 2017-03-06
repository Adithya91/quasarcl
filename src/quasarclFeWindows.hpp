#ifndef QUASARCL_FE_WIN_H_
#define QUASARCL_FE_WIN_H_

#include "QuasarCL.hpp"


namespace quasarcl {
	


	inline void reduceFeChisqs(QuasarCL& quasarcl, cl::Buffer& chisqs, cl::Buffer& sizes, const size_t size)
	{
		cl::Kernel reduceFeChisqsKernel = quasarcl.getKernelByName("reduce_fe_chisqs");
		
		unsigned int arg = 0;
		reduceFeChisqsKernel.setArg(arg++, sizeof(cl_double *), &chisqs);
		reduceFeChisqsKernel.setArg(arg++, sizeof(cl_uint *), &sizes);
		reduceFeChisqsKernel.setArg(arg++, static_cast<cl_uint>(size));

		cl::CommandQueue queue = quasarcl.getQueue();
		cl::Device device = queue.getInfo<CL_QUEUE_DEVICE>();
		
		size_t maxWorkGroupSize = device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>();
		size_t globalSize = calcGlobalSize(maxWorkGroupSize, size);

		cl::NDRange global(globalSize);
		cl::NDRange local(maxWorkGroupSize);

		queue.enqueueNDRangeKernel(
				reduceFeChisqsKernel,
				cl::NullRange,
				global,
				local);
		queue.finish();
	}

}

#endif