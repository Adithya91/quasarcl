#pragma OPENCL EXTENSION cl_khr_fp64 : enable

#define ASTRO_OBJ_SPEC_SIZE 4096

__kernel void generateWavelengthsMatrix
	(
		__global double4 * abz_buffer,
		__global double * wavelengths_matrix
	)
{
	// gid0 - numer widma kwazaru
	uint gid0 = get_global_id(0);	
	
	// parametry a, b i z kwazaru
	double4 abz_;
	__local double4 local_abz_;

	if (get_local_id(0) == 0)
	{
  		local_abz_ = abz_buffer[gid0];
	}
	barrier(CLK_LOCAL_MEM_FENCE);
	abz_ = local_abz_;
	
	// gid1 - numer elementu widma (indeks od 0 do 4095)
	uint gid1 = get_global_id(1);
	uint idx = gid0 * get_global_size(1) + gid1;
	
	// Wyliczenie lambdy dla tego gid1
	double wl = (abz_.x * (double)(gid1)) + abz_.y;
	wl = pow((double)(10),wl);
	
	// Uwzglednienie przesuniecia kosmologicznego - przejscie do ukladu emitujacego z wykorzystaniem redshiftu
	wavelengths_matrix[idx] = wl / (abz_.z + (double)(1));
	
	return;
}


// Dodaje do każdego z widm w spectrums_matrix jedno wskazane widmo
// z pojedynczą interpolacją.
//
// Widma, długości widm są tu ustawione kolumnami, a nie wierszami.
// dzięki temu otrzymujemy duży skok wydajności.
//
__kernel void addSpectrum
	(
		__global double * wavelengths_matrix,	// Długości dla widm
		__global double * spectrums_matrix,	// Widma po kolumnach
		__global uint  * sizes,			// Rozmiary widm
		uint size,				// Ilość widm (ilość kolumn)
		__constant double * to_add_wavelengths,	// Długości fal widma, które dodajemy
		__constant double * to_add_spectrum,	// Widmo, które dodajemy
		uint to_add_spectrum_size,	
		__global double * output_matrix		// Macierz zapisu wyniku sumy
	)
{
	// gid0 - numer widma kwazaru
	uint gid0 = get_global_id(0);
	if(gid0 >= size)
	{
		return;
	}	
	
	uint spectrum_size = sizes[gid0];

	// Indek elementu ze wavelengths_matrix/spectrums_matrix
	uint idx = gid0;
	uint idx_max = idx + spectrum_size * size;	
	
	// Indeks elementu z to_add_wavelengths/to_add_spectrum
	uint to_add_idx = 0;
	uint to_add_idx_max = to_add_spectrum_size - 1;
	
	double wavelength = wavelengths_matrix[idx];
	double value = spectrums_matrix[idx];
	
	double to_add_wl 	= to_add_wavelengths[to_add_idx];
	double to_add_wl_next 	= to_add_wavelengths[to_add_idx+1];		
	double to_add_value 	= to_add_spectrum[to_add_idx];
	double to_add_value_next = to_add_spectrum[to_add_idx+1];	
	

	int idx_max_flag = 0;
	int to_add_idx_max_flag = 0;

	while(1)
	{				
		if(wavelength >= to_add_wl)
		{			
			if(wavelength <= to_add_wl_next)
			{	
				double a = (to_add_value_next - to_add_value)/(to_add_wl_next - to_add_wl);
				double b = to_add_value - (a * to_add_wl);
				value = value + (a * wavelength + b);		
				output_matrix[idx] = value;		
				
				idx += size;
				// Sprawdzanie czy idx nie przekroczył idx_max
				// Zanim odczytamy dane.
				idx_max_flag = select(1, 0, idx < idx_max);
				if(idx_max_flag)
				{
					break;
				}
				wavelength = wavelengths_matrix[idx];
				value = spectrums_matrix[idx];
			}
			else
			{
				to_add_idx++;
				// Sprawdzanie czy to_add_idx nie przekroczył to_add_idx_max
				// Zanim odczytamy dane.
				to_add_idx_max_flag = select(1, 0, to_add_idx < to_add_idx_max);
				if(to_add_idx_max_flag)
				{
					break;
				}
				to_add_wl = to_add_wl_next;				
				to_add_wl_next = to_add_wavelengths[to_add_idx+1];
				
				to_add_value = to_add_value_next;
				to_add_value_next = to_add_spectrum[to_add_idx+1];
			}
		}
		else
		{	
			output_matrix[idx] = value;
		
			idx += size;
			// Sprawdzanie czy idx nie przekroczył idx_max
			// Zanim odczytamy dane.
			idx_max_flag = select(1, 0, idx < idx_max);
			if(idx_max_flag)
			{
				break;
			}
			wavelength = wavelengths_matrix[idx];
			value = spectrums_matrix[idx];
		}	
	}
	
	while(idx < idx_max)
	{
		value = spectrums_matrix[idx];
		output_matrix[idx] = value; 
		idx += size;
	}	
}




// Filtruje tylko te dane, których odpowiadająca wartość
// długości fali z wavelengths_matrix mieści się, w jakimś oknie widmowym.
//
__kernel void filterWithWavelengthWindows
	(
		__global double * wavelengths_matrix,	// Długości fal widm
		__global double * spectrums_matrix,	// Widma
		__global double * errors_matrix,		// Błędy pomiaru widm
		__global uint  * sizes,		 	// Rozmiary widm w spectrums_matrix
		__constant double2 * windows,	// Okna	długości fal
		uint windows_size 		// Ilość okien
	)	
{
	// gid0 - numer widma kwazaru
	uint gid0 = get_global_id(0);
	// gid1 - numer elementu w widmie
	uint gid1 = get_global_id(1);
	// indeks 
	uint idx = (gid0 * ASTRO_OBJ_SPEC_SIZE) + gid1;

	uint size = sizes[idx % get_global_size(0)];

	// Zero wskazuje na brak dopasownia do któregokolwiek okna
	int global_flag = 0;

	double wl_result = wavelengths_matrix[idx];
	double spec_result = spectrums_matrix[idx];
	double error_result = errors_matrix[idx];

	uint window_idx = 0;
	double2 window;
	__local double2 window_local;
	
	// Pętla po wszystkich oknach
	int window_flag;
	for(; window_idx < windows_size; window_idx++)
	{
		window_flag = 0;
		if (get_local_id(0) == 0)
		{
  			window_local = windows[window_idx];
		}
		//barrier(CLK_LOCAL_MEM_FENCE);

		window = window_local;

		window_flag = select(0, 1, wl_result >= window.x);
		window_flag *= select(0, 1, wl_result <= window.y);
		// Jeżeli oba warunki spełnione to mamy dopasowania
		// W przeciwnym przypadku zostawiamy to co było.

		global_flag = select(global_flag, 1, window_flag == 1);
	}
	
	spec_result = select((double)INFINITY, spec_result, (long)(global_flag == 1));
	wl_result = select((double)INFINITY, wl_result, (long)(global_flag == 1));
	error_result = select((double)INFINITY, error_result, (long)(global_flag == 1));	

	uint row_idx = idx / get_global_size(0);
	spec_result = select((double)INFINITY, spec_result, (long)(row_idx < size));
	wl_result = select((double)INFINITY, wl_result, (long)(row_idx < size));
	error_result = select((double)INFINITY, error_result, (long)(row_idx < size));

	wavelengths_matrix[idx] = wl_result;
	//spectrums_matrix[idx] = (double)(global_flag==1);
	spectrums_matrix[idx] = spec_result;
	errors_matrix[idx] = error_result;
}

/*__kernel void filterWithWavelengthWindows
	(
		__global double * wavelengths_matrix,	// Długości fal widm
		__global double * spectrums_matrix,	// Widma
		__global double * errors_matrix,		// Błędy pomiaru widm
		__global uint  * sizes,		 	// Rozmiary widm w spectrums_matrix
		__constant double2 * windows,	// Okna	długości fal
		uint windows_size 		// Ilość okien
	)	
{
	// gid0 - numer widma kwazaru
	uint gid0 = get_global_id(0);
	// gid1 - numer elementu w widmie
	uint gid1 = get_global_id(1);
	// indeks 
	uint idx = (gid0 * ASTRO_OBJ_SPEC_SIZE) + gid1;

	uint size = sizes[idx % get_global_size(0)];

	// Zero wskazuje na brak dopasownia do któregokolwiek okna
	int global_flag = 0;

	double wl_result = wavelengths_matrix[idx];
	double spec_result = spectrums_matrix[idx];
	double error_result = errors_matrix[idx];

	uint window_idx = 0;
	double2 window;
	__local double2 window_local;
	
	// Pętla po wszystkich oknach
	int window_flag;
	for(; window_idx < windows_size; window_idx++)
	{
		window_flag = 0;
		if (get_local_id(0) == 0)
		{
  			window_local = windows[window_idx];
		}
		//barrier(CLK_LOCAL_MEM_FENCE);

		window = window_local;

		if (wl_result >= window.x && wl_result <= window.y) {
			global_flag = 1;
		}
	}
	int row_idx = idx / get_global_size(0);
	if (!(global_flag == 1 && row_idx < size)) {
		wavelengths_matrix[idx] = (double)INFINITY;
		spectrums_matrix[idx] = (double)INFINITY;
		errors_matrix[idx] = (double)INFINITY;
	}
}*/

// Filtruje tylko te dane, których odpowiadające
// wartości ze spectrum_matrix są większe od zera.
//
__kernel void filterNonpositive
	(		
		__global double * spectrums_matrix,	// Widma
		__global double * a_matrix,	// Dowolna macierz wielkości spectrums_matrix
		__global double * b_matrix,	// Dowolna macierz wielkości spectrums_matrix
		__constant uint * sizes 	// Rozmiary widm w spectrums_matrix
	)
{
	uint gid0 = get_global_id(0);
	uint gid1 = get_global_id(1);
	// Indeks elementu
	uint idx = (gid0 * ASTRO_OBJ_SPEC_SIZE) + gid1;
	
	uint size = sizes[idx % get_global_size(0)];
	uint row_idx = idx / get_global_size(0);

	double spectrum = spectrums_matrix[idx];

	// Obliczenie flagi, czy spec_result > 0.
	// Zero wskazuje na liczbe mniejszą równą zero
	long flag = select(0, 1, spectrum >= (double)FLT_MIN);
	// Dodatkowe sprawdzenie, czy element jest znaczący.
	flag *= select(0, 1, row_idx < size);

	spectrum = select((double)INFINITY, spectrum, flag);
	spectrums_matrix[idx] = spectrum;

	double a = a_matrix[idx];
	a = select((double)INFINITY, a, flag);
	a_matrix[idx] = a;	

	double b = b_matrix[idx];
	b = select((double)INFINITY, b, flag);	
	b_matrix[idx] = b;
}


// Filtruje tylko te dane, których odpowiadające
// wartości ze spectrums_matrix, są różne od zera.
//
__kernel void filterZeros
	(
		__global double * spectrums_matrix,	// Widma
		__global double * a_matrix,	// Dowolna macierz wielkości spectrums_matrix
		__global double * b_matrix,	// Dowolna macierz wielkości spectrums_matrix
		__constant uint * sizes 	// Rozmiary widm w spectrums_matrix
	)
{
	//numer wiersza
	uint gid0 = get_global_id(0);
	//numer kolumny
	uint gid1 = get_global_id(1);
	// Indeks elementu
	uint idx = (gid0 * ASTRO_OBJ_SPEC_SIZE) + gid1;
	
	uint size = sizes[idx % get_global_size(0)];
	uint row_idx = idx / get_global_size(0);

	double spectrum = spectrums_matrix[idx];

	// Obliczenie flagi, czy spec_result != 0.
	long flag = select(0, 1, ((spectrum != 0.0) && (row_idx < size)));
	// Dodatkowe sprawdzenie, czy element jest znaczący.
	//flag* = select(0, 1, row_idx < size);
	//double inf = (double)INFINITY;
	
	spectrum = select((double)INFINITY, spectrum, flag);
	spectrums_matrix[idx] = spectrum;

	double a = a_matrix[idx];
	a = select((double)INFINITY, a, flag);
	a_matrix[idx] = a;	

	double b = b_matrix[idx];
	b = select((double)INFINITY, b, flag);	
	b_matrix[idx] = b;
}


// Filtruje tylko te dane, których odpowiadające
// wartości ze spectrums_matrix są różne od +/- INFINITY.
//

__kernel void filterInfs
	(
		__global double * spectrums_matrix,	// Widma
		__global double * a_matrix,	// Dowolna macierz wielkości spectrums_matrix
		__global double * b_matrix,	// Dowolna macierz wielkości spectrums_matrix
		__constant uint * sizes 	// Rozmiary widm w spectrums_matrix
	)
{
	uint gid0 = get_global_id(0);
	uint gid1 = get_global_id(1);
	// Indeks elementu
	uint idx = (gid0 * ASTRO_OBJ_SPEC_SIZE) + gid1;

	uint size = sizes[idx % get_global_size(0)];
	uint row_idx = idx / get_global_size(0);

	double spectrum = spectrums_matrix[idx];

	// Obliczenie flagi, czy spec_result != +/- INFINITY.
	long flag = select(0, 1, spectrum != (double)INFINITY);
	flag *= select(0, 1, spectrum != -(double)INFINITY);
	// Dodatkowe sprawdzenie, czy element jest znaczący.
	flag *= select(0, 1, row_idx < size);

	spectrum = select((double)INFINITY, spectrum, flag);
	spectrums_matrix[idx] = spectrum;

	double a = a_matrix[idx];
	a = select((double)INFINITY, a, flag);
	a_matrix[idx] = a;
		
	if(b_matrix != 0)
	{
		double b = b_matrix[idx];
		b = select((double)INFINITY, b, flag);	
		b_matrix[idx] = b;
	}
}
