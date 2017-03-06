#ifndef QUASARCL_BASICS_H_
#define QUASARCL_BASICS_H_

#include "QuasarCL.hpp"


namespace quasarcl {

	inline void log10(QuasarCL& quasarcl, cl::Buffer& input, const size_t width, const size_t height)
	{
		cl::Kernel log10Kernel = quasarcl.getKernelByName("matrix_log10");
		
		unsigned int arg = 0;
		log10Kernel.setArg(arg++, sizeof(cl_double *), &input);
		log10Kernel.setArg(arg++, static_cast<cl_uint>(width));

		cl::CommandQueue queue = quasarcl.getQueue();
		cl::Device device = queue.getInfo<CL_QUEUE_DEVICE>();
		
		size_t maxWorkGroupSize = device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>();
		size_t globalSize = calcGlobalSize(maxWorkGroupSize, width);

		cl::NDRange global(height, globalSize);
		cl::NDRange local(1, maxWorkGroupSize);

		queue.enqueueNDRangeKernel(
				log10Kernel,
				cl::NullRange,
				global,
				local);
		queue.finish();
	}


	
	inline void minus(QuasarCL& quasarcl, cl::Buffer& input, const size_t width, const size_t height,
			cl::Buffer& subtrahend_matrix, cl::Buffer& output)
	{
		cl::Kernel minusMatrixKernel = quasarcl.getKernelByName("matrix_minus_matrix");
		
		unsigned int arg = 0;
		minusMatrixKernel.setArg(arg++, sizeof(cl_double *), &input);
		minusMatrixKernel.setArg(arg++, static_cast<cl_uint>(width));
		minusMatrixKernel.setArg(arg++, sizeof(cl_double *), &subtrahend_matrix);
		minusMatrixKernel.setArg(arg++, sizeof(cl_double *), &output);

		cl::CommandQueue queue = quasarcl.getQueue();
		cl::Device device = queue.getInfo<CL_QUEUE_DEVICE>();
		
		size_t maxWorkGroupSize = device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>();
		size_t globalSize = calcGlobalSize(maxWorkGroupSize, width);

		cl::NDRange global(height, globalSize);
		cl::NDRange local(1, maxWorkGroupSize);

		queue.enqueueNDRangeKernel(
				minusMatrixKernel,
				cl::NullRange,
				global,
				local);
		queue.finish();
	}



	inline void divide(QuasarCL& quasarcl, cl::Buffer& dividend_matrix, const size_t width, const size_t height,
				cl::Buffer& divisor_matrix, cl::Buffer& output)
	{
		cl::Kernel divideMatrixKernel = quasarcl.getKernelByName("matrix_divide_matrix");
		
		unsigned int arg = 0;
		divideMatrixKernel.setArg(arg++, sizeof(cl_double *), &dividend_matrix);
		divideMatrixKernel.setArg(arg++, static_cast<cl_uint>(width));
		divideMatrixKernel.setArg(arg++, sizeof(cl_double *), &divisor_matrix);
		divideMatrixKernel.setArg(arg++, sizeof(cl_double *), &output);
		
		cl::CommandQueue queue = quasarcl.getQueue();
		cl::Device device = queue.getInfo<CL_QUEUE_DEVICE>();
		
		size_t maxWorkGroupSize = device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>();
		size_t globalSize = calcGlobalSize(maxWorkGroupSize, width);

		cl::NDRange global(height, globalSize);
		cl::NDRange local(1, maxWorkGroupSize);
		
		queue.enqueueNDRangeKernel(
				divideMatrixKernel,
				cl::NullRange,
				global,
				local);
		queue.finish();
	}



	inline void minus(QuasarCL& quasarcl, cl::Buffer& input, const size_t width, const size_t height,
			cl_double subtrahend)
	{
		cl::Kernel minusScalarKernel = quasarcl.getKernelByName("matrix_minus_scalar");
		
		unsigned int arg = 0;
		minusScalarKernel.setArg(arg++, sizeof(cl_double *), &input);
		minusScalarKernel.setArg(arg++, static_cast<cl_uint>(width));
		minusScalarKernel.setArg(arg++, subtrahend);

		cl::CommandQueue queue = quasarcl.getQueue();
		cl::Device device = queue.getInfo<CL_QUEUE_DEVICE>();
		
		size_t maxWorkGroupSize = device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>();
		size_t globalSize = calcGlobalSize(maxWorkGroupSize, width);

		cl::NDRange global(height, globalSize);
		cl::NDRange local(1, maxWorkGroupSize);

		queue.enqueueNDRangeKernel(
				minusScalarKernel,
				cl::NullRange,
				global,
				local);
		queue.finish();
	}



	inline void multiplyCol(QuasarCL& quasarcl, cl::Buffer& input, const size_t width, const size_t height,
					cl::Buffer& vector, cl::Buffer& output)
	{
		cl::Kernel multiplyColKernel = quasarcl.getKernelByName("matrix_multiply_col_vector");
		
		unsigned int arg = 0;
		multiplyColKernel.setArg(arg++, sizeof(cl_double *), &input);
		multiplyColKernel.setArg(arg++, static_cast<cl_uint>(width));
		multiplyColKernel.setArg(arg++, sizeof(cl_double *), &vector);
		multiplyColKernel.setArg(arg++, sizeof(cl_double *), &output);
		
		cl::CommandQueue queue = quasarcl.getQueue();
		cl::Device device = queue.getInfo<CL_QUEUE_DEVICE>();
		
		size_t maxWorkGroupSize = device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>();
		size_t globalSize = calcGlobalSize(maxWorkGroupSize, width);

		cl::NDRange global(height, globalSize);
		cl::NDRange local(1, maxWorkGroupSize);

		queue.enqueueNDRangeKernel(
				multiplyColKernel,
				cl::NullRange,
				global,
				local);
		queue.finish();
	}




	inline void transpose(QuasarCL& quasarcl, cl::Buffer& matrix, const size_t width, const size_t height, 
				cl::Buffer& tmatrix)
	{

		cl::Kernel transposeKernel = quasarcl.getKernelByName("matrix_transpose");
		size_t BLOCK_DIM = 16;

		unsigned int arg = 0;
		transposeKernel.setArg(arg++, sizeof(cl_double *), &matrix);
		transposeKernel.setArg(arg++, sizeof(cl_double*), &tmatrix);
		transposeKernel.setArg(arg++, static_cast<cl_uint>(width));
		transposeKernel.setArg(arg++, static_cast<cl_uint>(height));
		transposeKernel.setArg(arg++, cl::__local(sizeof(cl_double) * BLOCK_DIM * (BLOCK_DIM + 1)));
		
		cl::CommandQueue queue = quasarcl.getQueue();

		size_t globalSize1 = calcGlobalSize(BLOCK_DIM, width);
		size_t globalSize2 = calcGlobalSize(BLOCK_DIM, height);

		cl::NDRange global(globalSize1, globalSize2);
		cl::NDRange local(BLOCK_DIM, BLOCK_DIM);

		queue.enqueueNDRangeKernel(
				transposeKernel,
				cl::NullRange,
				global,
				local);
		queue.finish();
	}

}

#endif