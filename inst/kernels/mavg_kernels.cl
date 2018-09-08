#if __OPENCL_VERSION__ < 120
  #if cl_khr_fp64
     #pragma OPENCL EXTENSION cl_khr_fp64 : enable
  #elif cl_amd_fp64
     #pragma OPENCL EXTENSION cl_amd_fp64 : enable
  #else
     #error Missing double precision extension
  #endif
#endif


__kernel void simple_mavg
	(
		__global double * input,
		uint width,
		uint height,
		__global double * output,
		uint window_width
	)
{
	// kolumna 
	uint gid = get_global_id(0);
	uint idx = gid;

	if(gid >= width)
	{
		return;
	}
		
	unsigned long end = idx + height * width;
	uint i = 0;
	double lastSum = 0;
	double result = 0;

	while(i < window_width && idx < end)
	{
		lastSum += input[idx];		
		result = lastSum / ((double)(i+1));
		output[idx] = result;
		
		idx += width;	
		i++;
	}
	
	double fwindow_width = (double)(window_width);
	while(idx < end)
	{
		double new = input[idx];
		double old = input[idx - (window_width*width)];

		lastSum = lastSum - old + new;	
		result = lastSum / fwindow_width;
		output[idx] = result;

		idx += width;	
	}

	return;
}

__kernel void centered_mavg
	(
		__global double * input,
		uint width,
		uint height,
		__constant uint * cols_heights,
		__global double * output,
		uint window_width
	)
{
	// kolumna 
	uint gid = get_global_id(0);
	uint idx = gid;

	if(gid >= width)
	{
		return;
	}

	uint col_height = cols_heights[gid];
		
	unsigned long end = idx + col_height * width;
	uint i = 0;
	double lastSum = 0;
	double result = 0;

	//
	// Dla pierwszych window_width - 1 elementów [0; window_width - 1)
	//

	// Wartość średniej kroczącej dla pierwsztch window_width elementów jest
	// wyznaczana tak, że dla i-tego elementu jest równe śrędniej arytmetycznej
	// od elementu o indeksie 0 do elementu o indeksie i+window_width
	while(i < window_width && idx < end)
	{
		lastSum += input[idx];			
		idx += width;	
		i++;
	}

	i = 0;
	while(i < window_width && idx < end)
	{	
		result = lastSum / ((double)(window_width + i));
		output[idx - (window_width * width)] = result;
		
		double new = input[idx];
		lastSum = lastSum + new;

		idx += width;	
		i++;
	}

	//
	// Dla elementów o indeksie z przedziału [window_width; ilość_elementów - window_width]
	//

	// Wartość średniej kroczącej dla i-tego elementu jest średnią arytmetyczną elementów
	// o indeksach w przedziale (i - window_width; i + window_width)
	//idx = gid;
	//lastSum = 0;
	//i = 0;
	//while(i < 2 * window_width && idx < end)
	//{
	//	lastSum += input[idx];			
	//	idx += width;	
	//	i++;
	//}
	
	double fwindow_width = (double)(2 * window_width);
	while(idx < end)
	{			
		result = lastSum / fwindow_width;
		output[idx - (window_width * width)] = result;

		double new = input[idx];
		double old = input[idx - (2 * window_width * width)];
		lastSum = lastSum - old + new;	

		idx += width;	
	}
	
	//
	// Dla elementów o indeksie z przedziału (ilość_elementów - window_width; ilość_elementów]
	//

	// Wartość średniej kroczącej dla i-tego elementu jest średnią arytmetyczną elementów
	// o indeksach w przedziale (i; ilość_elementów)

	lastSum = 0.0f;
	idx -= 2 * window_width * width;
	while(idx < end)
	{		
		lastSum += input[idx];
		idx += width;	
	}

	idx -= window_width * width;	
	i = 2 * window_width;
	while(idx < end)
	{			
		result = lastSum / ((double)(i));
		output[idx] = result;
		
		double old = input[idx - window_width * width];
		lastSum = lastSum - old;	

		idx += width;	
		i--;
	}

	return;
}
