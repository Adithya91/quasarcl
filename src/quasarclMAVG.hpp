#ifndef QUASARCL_MAVG_H_
#define QUASARCL_MAVG_H_


#include "QuasarCL.hpp"


namespace quasarcl {
	

	inline void simpleMAVG(QuasarCL& quasarcl, cl::Buffer& in, const size_t width, const size_t height, cl::Buffer& out,
			const unsigned int window)
	{
		cl::Kernel smavgKernel = quasarcl.getKernelByName("simple_mavg");
		
		unsigned int arg = 0;
		smavgKernel.setArg(arg++, sizeof(cl_double *), &in);
		smavgKernel.setArg(arg++, static_cast<cl_uint>(width));
		smavgKernel.setArg(arg++, static_cast<cl_uint>(height));
		smavgKernel.setArg(arg++, sizeof(cl_double *), &out);
		smavgKernel.setArg(arg++, window);

		cl::CommandQueue queue = quasarcl.getQueue();
		cl::Device device = queue.getInfo<CL_QUEUE_DEVICE>();
		
		size_t maxWorkGroupSize = device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>();
		size_t globalSize = calcGlobalSize(maxWorkGroupSize, width);

		cl::NDRange global(globalSize);
		cl::NDRange local(maxWorkGroupSize);


		queue.enqueueNDRangeKernel(
				smavgKernel,
				cl::NullRange,
				global,
				local);
	}


	
	inline void centeredMAVG(QuasarCL& quasarcl, cl::Buffer& in, const size_t width, const size_t height, cl::Buffer& sizes, cl::Buffer& out,
			const unsigned int window)
	{
		cl::Kernel cmavgKernel = quasarcl.getKernelByName("centered_mavg");
		
		unsigned int arg = 0;
		cmavgKernel.setArg(arg++, sizeof(cl_double *), &in);
		cmavgKernel.setArg(arg++, static_cast<cl_uint>(width));
		cmavgKernel.setArg(arg++, static_cast<cl_uint>(height));
		cmavgKernel.setArg(arg++, sizeof(cl_uint *), &sizes);
		cmavgKernel.setArg(arg++, sizeof(cl_double *), &out);
		cmavgKernel.setArg(arg++, window);

		cl::CommandQueue queue = quasarcl.getQueue();
		cl::Device device = queue.getInfo<CL_QUEUE_DEVICE>();
		
		size_t maxWorkGroupSize = device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>();
		size_t globalSize = calcGlobalSize(maxWorkGroupSize, width);

		cl::NDRange global(globalSize);
		cl::NDRange local(maxWorkGroupSize);

		queue.enqueueNDRangeKernel(
				cmavgKernel,
				cl::NullRange,
				global,
				local);
	}
}

#endif