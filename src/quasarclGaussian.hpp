#ifndef QUASARCL_GAUSSIAN_H_
#define QUASARCL_GAUSSIAN_H_

#ifndef __CL_ENABLE_EXCEPTIONS
#define __CL_ENABLE_EXCEPTIONS
#endif

#define MAX_FITGAUSSIAN_LM_ITERS 500

#include "QuasarCL.hpp"


namespace quasarcl {

	inline void fitGaussian(QuasarCL& quasarcl, cl::Buffer& ys, cl::Buffer& xs, cl::Buffer& sizes, const size_t size, cl::Buffer& results)
	{
		cl::Kernel fitGaussianKernel = quasarcl.getKernelByName("fit_gaussian");
		
		unsigned int arg = 0;
		fitGaussianKernel.setArg(arg++, sizeof(cl_double *), &ys);
		fitGaussianKernel.setArg(arg++, sizeof(cl_double *), &xs);
		fitGaussianKernel.setArg(arg++, static_cast<cl_uint>(size));
		fitGaussianKernel.setArg(arg++, sizeof(cl_uint *), &sizes);
		fitGaussianKernel.setArg(arg++, static_cast<cl_uint>(MAX_FITGAUSSIAN_LM_ITERS));
		fitGaussianKernel.setArg(arg++, sizeof(cl_double4 *), &results);

		cl::CommandQueue queue = quasarcl.getQueue();
		cl::Device device = queue.getInfo<CL_QUEUE_DEVICE>();
		
		size_t workGroupMultiple = fitGaussianKernel.getWorkGroupInfo<CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE>(device);
		size_t globalSize = calcGlobalSize(workGroupMultiple, size);

		cl::NDRange global(globalSize);

		queue.enqueueNDRangeKernel(
				fitGaussianKernel,
				cl::NullRange,
				global,
				cl::NullRange);
	}


	
	inline void calcGaussian(QuasarCL& quasarcl, cl::Buffer& xs, cl::Buffer& gaussianParams, cl::Buffer& sizes, const size_t size, 
					  const size_t maxSpectrumSize, cl::Buffer& fxs)
	{
		cl::Kernel calcGaussianKernel = quasarcl.getKernelByName("calc_gaussian");
		
		unsigned int arg = 0;
		calcGaussianKernel.setArg(arg++, sizeof(cl_double *), &xs);
		calcGaussianKernel.setArg(arg++, sizeof(cl_double4 *), &gaussianParams);
		calcGaussianKernel.setArg(arg++, static_cast<cl_uint>(size));
		calcGaussianKernel.setArg(arg++, sizeof(cl_uint *), &sizes);
		calcGaussianKernel.setArg(arg++, sizeof(cl_double *), &fxs);

		cl::CommandQueue queue = quasarcl.getQueue();
		cl::Device device = queue.getInfo<CL_QUEUE_DEVICE>();
		
		size_t maxWorkGroupSize = device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>();
		size_t workGroupMultiple = calcGaussianKernel.getWorkGroupInfo<CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE>(device);
		size_t globalSize = calcGlobalSize(maxWorkGroupSize, maxSpectrumSize);

		cl::NDRange global(size, globalSize);
		cl::NDRange local(1, workGroupMultiple);

		queue.enqueueNDRangeKernel(
				calcGaussianKernel,
				cl::NullRange,
				global,
				local);
	}

	
	
	inline void calcGaussianChisq(QuasarCL& quasarcl, cl::Buffer& xs, cl::Buffer& ys, cl::Buffer& errors, cl::Buffer& gaussianParams,
			cl::Buffer& sizes, const size_t size, cl::Buffer& results)
	{
		cl::Kernel calcGaussianChisqKernel = quasarcl.getKernelByName("calc_gaussian_chisq");
		
		unsigned int arg = 0;
		calcGaussianChisqKernel.setArg(arg++, sizeof(cl_double *), &xs);
		calcGaussianChisqKernel.setArg(arg++, sizeof(cl_double *), &ys);
		calcGaussianChisqKernel.setArg(arg++, sizeof(cl_double *), &errors);
		calcGaussianChisqKernel.setArg(arg++, sizeof(cl_double4 *), &gaussianParams);
		calcGaussianChisqKernel.setArg(arg++, static_cast<cl_uint>(size));
		calcGaussianChisqKernel.setArg(arg++, sizeof(cl_uint *), &sizes);
		calcGaussianChisqKernel.setArg(arg++, sizeof(cl_double *), &results);

		cl::CommandQueue queue = quasarcl.getQueue();
		cl::Device device = queue.getInfo<CL_QUEUE_DEVICE>();
		
		size_t workGroupMultiple = calcGaussianChisqKernel.getWorkGroupInfo<CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE>(device);
		size_t globalSize = calcGlobalSize(workGroupMultiple, size);

		cl::NDRange global(globalSize);
		cl::NDRange local(workGroupMultiple);
		
		queue.enqueueNDRangeKernel(
				calcGaussianChisqKernel,
				cl::NullRange,
				global,
				local);
	}


	
	inline void calcGaussianFWHM(QuasarCL& quasarcl, cl::Buffer& gaussianParams, cl::Buffer& gaussianFWHMs, const size_t size)
	{
		cl::Kernel calcGaussianFWHMKernel = quasarcl.getKernelByName("calc_gaussian_fwhm");
		unsigned int arg = 0;
		calcGaussianFWHMKernel.setArg(arg++, sizeof(cl_double *), &gaussianParams);
		calcGaussianFWHMKernel.setArg(arg++, sizeof(cl_double *), &gaussianFWHMs);
		calcGaussianFWHMKernel.setArg(arg++, static_cast<cl_uint>(size));

		cl::CommandQueue queue = quasarcl.getQueue();
		cl::Device device = queue.getInfo<CL_QUEUE_DEVICE>();
		
		size_t maxWorkGroupSize = device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>();
		size_t workGroupMultiple = calcGaussianFWHMKernel.getWorkGroupInfo<CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE>(device);
		size_t globalSize = calcGlobalSize(maxWorkGroupSize, size);

		cl::NDRange global(globalSize);
		cl::NDRange local(workGroupMultiple);

		queue.enqueueNDRangeKernel(
				calcGaussianFWHMKernel,
				cl::NullRange,
				global,
				local);
	}

}

#endif