#ifndef QUASARCL_CONTINUUM_H_
#define QUASARCL_CONTINUUM_H_

#define ASTRO_OBJ_SPEC_SIZE 4096

#include "QuasarCL.hpp"


namespace quasarcl {
	

	inline void fixReglinResults(QuasarCL& quasarcl, cl::Buffer& cReglinResults, cl::Buffer& reglinResults, const size_t size)
	{
		cl::Kernel fixReglinResultsKernel = quasarcl.getKernelByName("fix_reglin_results");
		
		unsigned int arg = 0;
		fixReglinResultsKernel.setArg(arg++, sizeof(cl_double8 *), &cReglinResults);
		fixReglinResultsKernel.setArg(arg++, sizeof(cl_double8 *), &reglinResults);
		fixReglinResultsKernel.setArg(arg++, static_cast<cl_uint>(size));

		cl::CommandQueue queue = quasarcl.getQueue();
		cl::Device device = quasarcl.getQueue().getInfo<CL_QUEUE_DEVICE>();
		
		size_t maxWorkGroupSize = device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>();
		size_t globalSize = calcGlobalSize(maxWorkGroupSize, size);

		cl::NDRange global(globalSize);
		cl::NDRange local(maxWorkGroupSize);

		quasarcl.getQueue().enqueueNDRangeKernel(
				fixReglinResultsKernel,
				cl::NullRange,
				global,
				local);
		queue.finish();
	}

	

	inline void calcCfunDcfun(QuasarCL& quasarcl, cl::Buffer& wavelengthsMatrix, cl::Buffer& dcontinuumsMatrix, cl::Buffer& continuumsMatrix,
					   const size_t size, cl::Buffer& cReglinResults, cl::Buffer& reglinResults)
	{
		cl::Kernel calcCfunDcfunKernel = quasarcl.getKernelByName("calc_cfun_dcfun");
		
		unsigned int arg = 0;
		calcCfunDcfunKernel.setArg(arg++, sizeof(cl_double *), &wavelengthsMatrix);
		calcCfunDcfunKernel.setArg(arg++, sizeof(cl_double *), &dcontinuumsMatrix);
		calcCfunDcfunKernel.setArg(arg++, sizeof(cl_double *), &continuumsMatrix);
		calcCfunDcfunKernel.setArg(arg++, sizeof(cl_double8 *), &cReglinResults);
		calcCfunDcfunKernel.setArg(arg++, sizeof(cl_double8 *), &reglinResults);

		cl::CommandQueue queue = quasarcl.getQueue();
		cl::Device device = queue.getInfo<CL_QUEUE_DEVICE>();
		
		cl::NDRange global(size, ASTRO_OBJ_SPEC_SIZE);
		cl::NDRange local(1, device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>());

		queue.enqueueNDRangeKernel(
				calcCfunDcfunKernel,
				cl::NullRange,
				global,
				local);
		queue.finish();
	}


	
	inline void calcCw(QuasarCL& quasarcl, cl::Buffer& wavelengthsMatrix, cl::Buffer& continuumsMatrix, const size_t maxSpectrumSize, 
				const size_t size, cl::Buffer& cReglinResults)
	{
		cl::Kernel calcCwKernel = quasarcl.getKernelByName("calc_cw");
		
		unsigned int arg = 0;
		calcCwKernel.setArg(arg++, sizeof(cl_double *), &wavelengthsMatrix);
		calcCwKernel.setArg(arg++, sizeof(cl_double *), &continuumsMatrix);
		calcCwKernel.setArg(arg++, static_cast<cl_uint>(maxSpectrumSize));
		calcCwKernel.setArg(arg++, sizeof(cl_double8 *), &cReglinResults);

		cl::CommandQueue queue = quasarcl.getQueue();
		cl::Device device = queue.getInfo<CL_QUEUE_DEVICE>();
		
		size_t maxWorkGroupSize = device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>();
		size_t globalSize = calcGlobalSize(maxWorkGroupSize, maxSpectrumSize);

		cl::NDRange global(size, globalSize);
		cl::NDRange local(1, maxWorkGroupSize);

		queue.enqueueNDRangeKernel(
				calcCwKernel,
				cl::NullRange,
				global,
				local);
		queue.finish();
	}

	

	inline void reduceContinuumChisqs(QuasarCL& quasarcl, cl::Buffer& chisqs, cl::Buffer& sizes, const size_t size)
	{
		cl::Kernel reduceContinuumChisqsKernel = quasarcl.getKernelByName("reduce_continuum_chisqs");
		
		unsigned int arg = 0;
		reduceContinuumChisqsKernel.setArg(arg++, sizeof(cl_double *), &chisqs);
		reduceContinuumChisqsKernel.setArg(arg++, sizeof(cl_uint *), &sizes);
		reduceContinuumChisqsKernel.setArg(arg++, static_cast<cl_uint>(size));

		cl::CommandQueue queue = quasarcl.getQueue();
		cl::Device device = queue.getInfo<CL_QUEUE_DEVICE>();
		
		size_t maxWorkGroupSize = device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>();
		size_t globalSize = calcGlobalSize(maxWorkGroupSize, size);

		cl::NDRange global(globalSize);
		cl::NDRange local(maxWorkGroupSize);

		queue.enqueueNDRangeKernel(
				reduceContinuumChisqsKernel,
				cl::NullRange,
				global,
				local);
		queue.finish();
	}

}

#endif