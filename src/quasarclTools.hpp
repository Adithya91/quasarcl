#ifndef QUASARCL_TOOLS_H_
#define QUASARCL_TOOLS_H_


#include "QuasarCL.hpp"



namespace quasarcl {
	

	inline void convolve(QuasarCL& quasarcl, cl::Buffer& input, const size_t width, const size_t height, cl::Buffer& sizes, 
				  cl::Buffer& filter, const size_t filterSize, cl::Buffer& output)
	{
		cl::Kernel convolveKernel = quasarcl.getKernelByName("convolve");
		
		unsigned int arg = 0;
		convolveKernel.setArg(arg++, sizeof(cl_double *), &input);
		convolveKernel.setArg(arg++, sizeof(cl_uint *), &sizes);
		convolveKernel.setArg(arg++, sizeof(cl_double *), &filter);
		convolveKernel.setArg(arg++, static_cast<cl_uint>(filterSize));
		convolveKernel.setArg(arg++, sizeof(cl_double *), &output);

		cl::CommandQueue queue = quasarcl.getQueue();
		cl::Device device = queue.getInfo<CL_QUEUE_DEVICE>();

		cl::NDRange global(width, height);
		cl::NDRange local(1, device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>());

		queue.enqueueNDRangeKernel(
				convolveKernel,
				cl::NullRange,
				global,
				local);
		queue.finish();
	}


	
	inline void copyIfNotInf(QuasarCL& quasarcl, cl::Buffer& input, const size_t width, const size_t height, 
					  cl::Buffer& output, const size_t outputHeight)
	{
		cl::Kernel copyIfNotInfKernel = quasarcl.getKernelByName("copyIfNotInf");
		
		unsigned int arg = 0;
		copyIfNotInfKernel.setArg(arg++, sizeof(cl_double *), &input);
		copyIfNotInfKernel.setArg(arg++, static_cast<cl_uint>(width));
		copyIfNotInfKernel.setArg(arg++, static_cast<cl_uint>(height));
		copyIfNotInfKernel.setArg(arg++, sizeof(cl_double *), &output);
		copyIfNotInfKernel.setArg(arg++, static_cast<cl_uint>(outputHeight));

		cl::CommandQueue queue = quasarcl.getQueue();
		cl::Device device = queue.getInfo<CL_QUEUE_DEVICE>();
		
		size_t maxWorkGroupSize = device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>();
		size_t workGroupMultiple = copyIfNotInfKernel.getWorkGroupInfo<CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE>(device);
		size_t globalSize = calcGlobalSize(maxWorkGroupSize, width);

		cl::NDRange global(globalSize);
		cl::NDRange local(workGroupMultiple);

		queue.enqueueNDRangeKernel(
				copyIfNotInfKernel,
				cl::NullRange,
				global,
				local);
		queue.finish();
	}


	
	inline void countIfNotInf(QuasarCL& quasarcl, cl::Buffer& input, const size_t width, const size_t height, 
					   cl::Buffer& sizes)
	{
		cl::Kernel countIfNotInfKernel = quasarcl.getKernelByName("countIfNotInf");
		
		unsigned int arg = 0;
		countIfNotInfKernel.setArg(arg++, sizeof(cl_double *), &input);
		countIfNotInfKernel.setArg(arg++, static_cast<cl_uint>(width));
		countIfNotInfKernel.setArg(arg++, static_cast<cl_uint>(height));
		countIfNotInfKernel.setArg(arg++, sizeof(cl_uint *), &sizes);

		cl::CommandQueue queue = quasarcl.getQueue();
		cl::Device device = queue.getInfo<CL_QUEUE_DEVICE>();
		
		size_t maxWorkGroupSize = device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>();
		size_t workGroupMultiple = countIfNotInfKernel.getWorkGroupInfo<CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE>(device);
		size_t globalSize = calcGlobalSize(maxWorkGroupSize, width);

		cl::NDRange global(globalSize);
		cl::NDRange local(workGroupMultiple);

		queue.enqueueNDRangeKernel(
				countIfNotInfKernel,
				cl::NullRange,
				global,
				local);
		queue.finish();
	}

	

	inline void reglin(QuasarCL& quasarcl, cl::Buffer& xs, cl::Buffer& ys, const size_t width, const size_t height, cl::Buffer& sizes, cl::Buffer& output)
	{
		cl::Kernel reglinKernel = quasarcl.getKernelByName("reglin");
		
		unsigned int arg = 0;
		reglinKernel.setArg(arg++, sizeof(cl_double *), &xs);
		reglinKernel.setArg(arg++, sizeof(cl_double *), &ys);
		reglinKernel.setArg(arg++, static_cast<cl_uint>(width));
		reglinKernel.setArg(arg++, static_cast<cl_uint>(height));
		reglinKernel.setArg(arg++, sizeof(cl_uint *), &sizes);
		reglinKernel.setArg(arg++, sizeof(cl_double8 *), &output);

		cl::CommandQueue queue = quasarcl.getQueue();
		cl::Device device = queue.getInfo<CL_QUEUE_DEVICE>();
		
		size_t maxWorkGroupSize = device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>();
		size_t workGroupMultiple = reglinKernel.getWorkGroupInfo<CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE>(device);
		size_t globalSize = calcGlobalSize(maxWorkGroupSize, width);

		cl::NDRange global(globalSize);
		cl::NDRange local(workGroupMultiple);

		queue.enqueueNDRangeKernel(
				reglinKernel,
				cl::NullRange,
				global,
				local);
		queue.finish();
	}

	

	inline void reglinYax(QuasarCL& quasarcl, cl::Buffer& xs, cl::Buffer& ys, const size_t width, const size_t height, cl::Buffer& sizes, 
				   cl::Buffer& output)
	{
		cl::Kernel reglinYaxKernel = quasarcl.getKernelByName("reglin_yax");
		
		unsigned int arg = 0;
		reglinYaxKernel.setArg(arg++, sizeof(cl_double *), &xs);
		reglinYaxKernel.setArg(arg++, sizeof(cl_double *), &ys);
		reglinYaxKernel.setArg(arg++, static_cast<cl_uint>(width));
		reglinYaxKernel.setArg(arg++, static_cast<cl_uint>(height));
		reglinYaxKernel.setArg(arg++, sizeof(cl_uint *), &sizes);
		reglinYaxKernel.setArg(arg++, sizeof(cl_double8 *), &output);
		
		cl::CommandQueue queue = quasarcl.getQueue();
		cl::Device device = queue.getInfo<CL_QUEUE_DEVICE>();
		
		size_t maxWorkGroupSize = device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>();
		size_t workGroupMultiple = reglinYaxKernel.getWorkGroupInfo<CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE>(device);
		size_t globalSize = calcGlobalSize(maxWorkGroupSize, width);

		cl::NDRange global(globalSize);
		cl::NDRange local(workGroupMultiple);

		queue.enqueueNDRangeKernel(
				reglinYaxKernel,
				cl::NullRange,
				global,
				local);
		queue.finish();
	}
	
	

	inline void chisq(QuasarCL& quasarcl, cl::Buffer& fs, cl::Buffer& ys, cl::Buffer& errors, const size_t width, const size_t height, 
			   cl::Buffer& sizes, cl::Buffer& output)
	{
		cl::Kernel chisqKernel = quasarcl.getKernelByName("chisq");
		
		unsigned int arg = 0;
		chisqKernel.setArg(arg++, sizeof(cl_double *), &fs);
		chisqKernel.setArg(arg++, sizeof(cl_double *), &ys);
		chisqKernel.setArg(arg++, sizeof(cl_double *), &errors);
		chisqKernel.setArg(arg++, static_cast<cl_uint>(width));
		chisqKernel.setArg(arg++, static_cast<cl_uint>(height));
		chisqKernel.setArg(arg++, sizeof(cl_uint *), &sizes);
		chisqKernel.setArg(arg++, sizeof(cl_double *), &output);

		cl::CommandQueue queue = quasarcl.getQueue();
		cl::Device device = queue.getInfo<CL_QUEUE_DEVICE>();
		
		size_t maxWorkGroupSize = device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>();
		size_t workGroupMultiple = chisqKernel.getWorkGroupInfo<CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE>(device);
		size_t globalSize = calcGlobalSize(maxWorkGroupSize, width);

		cl::NDRange global(globalSize);
		cl::NDRange local(workGroupMultiple);

		queue.enqueueNDRangeKernel(
				chisqKernel,
				cl::NullRange,
				global,
				local);
		queue.finish();
	}


	
	inline void trapz(QuasarCL& quasarcl, cl::Buffer& ys, cl::Buffer& xs, const size_t width, const size_t height, cl::Buffer& sizes, 
			   cl::Buffer& output)
	{
		cl::Kernel trapzKernel = quasarcl.getKernelByName("integrate_trapz");
		
		unsigned int arg = 0;
		trapzKernel.setArg(arg++, sizeof(cl_double *), &ys);
		trapzKernel.setArg(arg++, sizeof(cl_double *), &xs);
		trapzKernel.setArg(arg++, static_cast<cl_uint>(width));
		trapzKernel.setArg(arg++, static_cast<cl_uint>(height));
		trapzKernel.setArg(arg++, sizeof(cl_uint *), &sizes);
		trapzKernel.setArg(arg++, sizeof(cl_double8 *), &output);

		cl::CommandQueue queue = quasarcl.getQueue();
		cl::Device device = queue.getInfo<CL_QUEUE_DEVICE>();
		
		size_t maxWorkGroupSize = device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>();
		size_t workGroupMultiple = trapzKernel.getWorkGroupInfo<CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE>(device);
		size_t globalSize = calcGlobalSize(maxWorkGroupSize, width);

		cl::NDRange global(globalSize);
		cl::NDRange local(workGroupMultiple);

		queue.enqueueNDRangeKernel(
				trapzKernel,
				cl::NullRange,
				global,
				local);
		queue.finish();
	}

}

#endif