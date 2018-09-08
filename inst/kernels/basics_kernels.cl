#if __OPENCL_VERSION__ < 120
  #if cl_khr_fp64
     #pragma OPENCL EXTENSION cl_khr_fp64 : enable
  #elif cl_amd_fp64
     #pragma OPENCL EXTENSION cl_amd_fp64 : enable
  #else
     #error Missing double precision extension
  #endif
#endif


__kernel void matrix_log10
	(
		__global double * matrix,	// Macierz
		uint row_size			// Rozmiar wiersza macierzy.
	)
{
	// gid0 - numer wiersza macierzy input
	uint gid0 = get_global_id(0);
	// gid1 - numer elementu w wierszu.
	uint gid1 = get_global_id(1);

	if(gid1 >= row_size)
	{
		return;
	}
	
	uint idx = gid0 * row_size + gid1;

	double m = matrix[idx];
	m = log10(m);
	matrix[idx] = m;
}

__kernel void matrix_minus_scalar
	(
		__global double * matrix,	// Macierz
		uint row_size,			// Rozmiar wiersza macierzy.
		double subtrahend		// Liczba do odjęcia
	)
{
	// gid0 - numer wiersza macierzy input
	uint gid0 = get_global_id(0);
	// gid1 - numer elementu w wierszu.
	uint gid1 = get_global_id(1);

	if(gid1 >= row_size)
	{
		return;
	}
	
	uint idx = gid0 * row_size + gid1;

	double m = matrix[idx];
	m -= subtrahend;
	matrix[idx] = m;
}

__kernel void matrix_minus_matrix
	(
		__global double * matrix,	// Macierz
		uint row_size,			// Rozmiar wiersza macierzy.
		__global double * subtrahend_matrix,		// Macierz do odjęcia
		__global double * output_matrix	// Wynik
	)
{
	// gid0 - numer wiersza macierzy input
	uint gid0 = get_global_id(0);
	// gid1 - numer elementu w wierszu.
	uint gid1 = get_global_id(1);

	if(gid1 >= row_size)
	{
		return;
	}
	
	uint idx = gid0 * row_size + gid1;

	double value = matrix[idx];
	double subtrahend = subtrahend_matrix[idx];

	value -= subtrahend;

	output_matrix[idx] = value;
}

__kernel void matrix_divide_matrix
	(
		__global double * dividend_matrix,	// Macierz
		uint row_size,				// Rozmiar wiersza macierzy.
		__global double * divisor_matrix,	// Macierz 
		__global double * output_matrix		// Wynik
	)
{
	// gid0 - numer wiersza macierzy input
	uint gid0 = get_global_id(0);
	// gid1 - numer elementu w wierszu.
	uint gid1 = get_global_id(1);

	if(gid1 >= row_size)
	{
		return;
	}
	
	uint idx = gid0 * row_size + gid1;

	double dividend = dividend_matrix[idx];
	double divisor = divisor_matrix[idx];

	output_matrix[idx] = dividend /= divisor;
}

// Mnoży każdy element z i-tej kolumny przez i-ty element
// podanego wektora.
//
__kernel void matrix_multiply_col_vector
	(
		__global double * matrix,	// Macierz
		uint row_size,			// Rozmiar wiersza macierzy.
		__constant double * vector,	// Wektor, których zawiera co najmniej
						// tyle elementów ile kolumn ma matrix.
		__global double * output_matrix	// Wynik
	)
{
	// gid0 - numer wiersza macierzy input
	uint gid0 = get_global_id(0);
	// gid1 - numer elementu w wierszu (numer kolumny).
	uint gid1 = get_global_id(1);

	if(gid1 >= row_size)
	{
		return;
	}
	
	uint idx = gid0 * row_size + gid1;
	uint col_idx = gid1;

	double value = matrix[idx];
	double m = vector[col_idx];

	value *= m;

	output_matrix[idx] = value;
}

#define M_TRANSPOSE_BLOCK_DIM 16

__kernel void matrix_transpose
	(
		__global double * matrix,	// Macierz wejściowa
		__global double * tmatrix,	// Macierz transponowana
		uint width,			// Ilość kolumn, szerokość macierzy
		uint height,			// Ilość wierszy, wysokość macierzy
		__local double * scratch
	)
{
	// gid0 - numer wiersza
	uint x_idx = get_global_id(0);
	// gid1 - numer kolumny
	uint y_idx = get_global_id(1);
	
	uint idx;

	// Pobieranie wartości z matrix	
	if((x_idx < width) && (y_idx < height))
	{	
		idx = y_idx * width + x_idx;
		scratch[get_local_id(1)*(M_TRANSPOSE_BLOCK_DIM+1)+get_local_id(0)] = matrix[idx];
	}
	barrier(CLK_LOCAL_MEM_FENCE);


	// Pobieranie wartości z matrix	
	x_idx = get_group_id(1) * M_TRANSPOSE_BLOCK_DIM + get_local_id(0);
	y_idx = get_group_id(0) * M_TRANSPOSE_BLOCK_DIM + get_local_id(1);
	if((x_idx < height) && (y_idx < width))
	{	
		idx = y_idx * height + x_idx;
		tmatrix[idx] = scratch[get_local_id(0)*(M_TRANSPOSE_BLOCK_DIM+1)+get_local_id(1)];
	}
}



























