#ifndef QUASARCL_SPECTRUMS_H_
#define QUASARCL_SPECTRUMS_H_


#include "QuasarCL.hpp"


namespace quasarcl {
	

	inline void generateWavelengthsMatrix(QuasarCL& quasarcl, cl::Buffer& abz_, const size_t size, cl::Buffer& wavelengthsMatrix)
	{
		cl::Kernel generateWavelengthsMatrixKernel = quasarcl.getKernelByName("generateWavelengthsMatrix");
		
		unsigned int arg = 0;
		generateWavelengthsMatrixKernel.setArg(arg++, sizeof(cl_double4 *), &abz_);
		generateWavelengthsMatrixKernel.setArg(arg++, sizeof(cl_double *), &wavelengthsMatrix);

		cl::CommandQueue queue = quasarcl.getQueue();
		cl::Device device = queue.getInfo<CL_QUEUE_DEVICE>();
		
		cl::NDRange global(size, ASTRO_OBJ_SPEC_SIZE);
		cl::NDRange local(1, device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>());
		
		queue.enqueueNDRangeKernel(
				generateWavelengthsMatrixKernel,
				cl::NullRange,
				global,
				local);
	}

	

	inline void addSpectrum(QuasarCL& quasarcl, cl::Buffer& spectrumsMatrix, cl::Buffer& wavelengthsMatrix, cl::Buffer& sizes, const size_t size, 
					 cl::Buffer& toAddWavelengths, cl::Buffer& toAddSpectrum, const size_t toAddSize, cl::Buffer& outputSpectrumsMatrix)
	{
		cl::Kernel addSpectrumKernel = quasarcl.getKernelByName("addSpectrum");
		
		unsigned int arg = 0;
		addSpectrumKernel.setArg(arg++, sizeof(cl_double *), &wavelengthsMatrix);
		addSpectrumKernel.setArg(arg++, sizeof(cl_double *), &spectrumsMatrix);
		addSpectrumKernel.setArg(arg++, sizeof(cl_uint *), &sizes);
		addSpectrumKernel.setArg(arg++, static_cast<cl_uint>(size));
		addSpectrumKernel.setArg(arg++, sizeof(cl_double *), &toAddWavelengths);
		addSpectrumKernel.setArg(arg++, sizeof(cl_double *), &toAddSpectrum);
		addSpectrumKernel.setArg(arg++, static_cast<cl_uint>(toAddSize));
		addSpectrumKernel.setArg(arg++, sizeof(cl_double *), &outputSpectrumsMatrix);

		cl::CommandQueue queue = quasarcl.getQueue();
		cl::Device device = queue.getInfo<CL_QUEUE_DEVICE>();
		
		size_t maxWorkGroupSize = device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>();
		size_t globalSize = calcGlobalSize(maxWorkGroupSize, size);

		cl::NDRange global(globalSize);
		cl::NDRange local(maxWorkGroupSize);

		queue.enqueueNDRangeKernel(
				addSpectrumKernel,
				cl::NullRange,
				global,
				local);
	}


	
	inline void filterWithWavelengthWindows(QuasarCL& quasarcl, cl::Buffer& spectrumsMatrix, cl::Buffer& wavelengthsMatrix, cl::Buffer& errorsMatrix, 
									 cl::Buffer& sizes, const size_t size, cl::Buffer& windows, const size_t windowsSize)
	{
		cl::Kernel filterWithWavelengthWindowsKernel = quasarcl.getKernelByName("filterWithWavelengthWindows");
		
		unsigned int arg = 0;
		filterWithWavelengthWindowsKernel.setArg(arg++, sizeof(cl_double *), &wavelengthsMatrix);
		filterWithWavelengthWindowsKernel.setArg(arg++, sizeof(cl_double *), &spectrumsMatrix);
		filterWithWavelengthWindowsKernel.setArg(arg++, sizeof(cl_double *), &errorsMatrix);
		filterWithWavelengthWindowsKernel.setArg(arg++, sizeof(cl_uint *), &sizes);
		filterWithWavelengthWindowsKernel.setArg(arg++, sizeof(cl_double2 *), &windows);
		filterWithWavelengthWindowsKernel.setArg(arg++, static_cast<cl_uint>(windowsSize));

		cl::CommandQueue queue = quasarcl.getQueue();
		cl::Device device = queue.getInfo<CL_QUEUE_DEVICE>();

		cl::NDRange global(size, ASTRO_OBJ_SPEC_SIZE);
		cl::NDRange local(1, device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>());

		queue.enqueueNDRangeKernel(
				filterWithWavelengthWindowsKernel,
				cl::NullRange,
				global,
				local);
	}

	
	
	inline void filterNonpositive(QuasarCL& quasarcl, cl::Buffer& spectrumsMatrix, cl::Buffer& matrixA, cl::Buffer& matrixB, 
						   cl::Buffer& sizes, const size_t size)
	{
		cl::Kernel filterNonpositiveKernel = quasarcl.getKernelByName("filterNonpositive");
		
		unsigned int arg = 0;
		filterNonpositiveKernel.setArg(arg++, sizeof(cl_double *), &spectrumsMatrix);
		filterNonpositiveKernel.setArg(arg++, sizeof(cl_double *), &matrixA);
		filterNonpositiveKernel.setArg(arg++, sizeof(cl_double *), &matrixB);
		filterNonpositiveKernel.setArg(arg++, sizeof(cl_uint *), &sizes);

		cl::CommandQueue queue = quasarcl.getQueue();
		cl::Device device = queue.getInfo<CL_QUEUE_DEVICE>();
		
		cl::NDRange global(size, ASTRO_OBJ_SPEC_SIZE);
		cl::NDRange local(1, device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>());

		queue.enqueueNDRangeKernel(
				filterNonpositiveKernel,
				cl::NullRange,
				global,
				local);
	}



	inline void filterZeros(QuasarCL& quasarcl, cl::Buffer& spectrumsMatrix, cl::Buffer& matrixA, cl::Buffer& matrixB, 
					 cl::Buffer& sizes, const size_t size)
	{
		cl::Kernel filterZerosKernel = quasarcl.getKernelByName("filterZeros");

		unsigned int arg = 0;
		filterZerosKernel.setArg(arg++, sizeof(cl_double *), &spectrumsMatrix);
		filterZerosKernel.setArg(arg++, sizeof(cl_double *), &matrixA);
		filterZerosKernel.setArg(arg++, sizeof(cl_double *), &matrixB);
		filterZerosKernel.setArg(arg++, sizeof(cl_uint *), &sizes);

		cl::CommandQueue queue = quasarcl.getQueue();
		cl::Device device = queue.getInfo<CL_QUEUE_DEVICE>();

		cl::NDRange global(size, ASTRO_OBJ_SPEC_SIZE);
		cl::NDRange local(1, device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>());

		queue.enqueueNDRangeKernel(
				filterZerosKernel,
				cl::NullRange,
				global,
				local);
	}


	
	inline void filterInfs(QuasarCL& quasarcl, cl::Buffer& spectrumsMatrix, cl::Buffer& matrixA, cl::Buffer& matrixB, 
					cl::Buffer& sizes, const size_t size)
	{
		cl::Kernel filterInfsKernel = quasarcl.getKernelByName("filterInfs");
		
		unsigned int arg = 0;
		filterInfsKernel.setArg(arg++, sizeof(cl_double *), &spectrumsMatrix);
		filterInfsKernel.setArg(arg++, sizeof(cl_double *), &matrixA);
		filterInfsKernel.setArg(arg++, sizeof(cl_double *), &matrixB);
		filterInfsKernel.setArg(arg++, sizeof(cl_uint *), &sizes);

		cl::CommandQueue queue = quasarcl.getQueue();
		cl::Device device = queue.getInfo<CL_QUEUE_DEVICE>();

		cl::NDRange global(size, ASTRO_OBJ_SPEC_SIZE);
		cl::NDRange local(1, device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>());

		queue.enqueueNDRangeKernel(
				filterInfsKernel,
				cl::NullRange,
				global,
				local);
	}


	
	inline void filterInfs(QuasarCL& quasarcl, cl::Buffer& spectrumsMatrix, cl::Buffer& matrixA, 
					cl::Buffer& sizes, const size_t size)
	{
		cl::Kernel filterInfsKernel = quasarcl.getKernelByName("filterInfs");
		
		unsigned int arg = 0;
		filterInfsKernel.setArg(arg++, sizeof(cl_double *), &spectrumsMatrix);
		filterInfsKernel.setArg(arg++, sizeof(cl_double *), &matrixA);
		filterInfsKernel.setArg(arg++, sizeof(cl_double *), nullptr);
		filterInfsKernel.setArg(arg++, sizeof(cl_uint *), &sizes);

		cl::CommandQueue queue = quasarcl.getQueue();
		cl::Device device = queue.getInfo<CL_QUEUE_DEVICE>();

		cl::NDRange global(size, ASTRO_OBJ_SPEC_SIZE);
		cl::NDRange local(1, device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>());
		
		queue.enqueueNDRangeKernel(
				filterInfsKernel,
				cl::NullRange,
				global,
				local);
	}

}


#endif