#if __OPENCL_VERSION__ < 120
  #if cl_khr_fp64
     #pragma OPENCL EXTENSION cl_khr_fp64 : enable
  #elif cl_amd_fp64
     #pragma OPENCL EXTENSION cl_amd_fp64 : enable
  #else
     #error Missing double precision extension
  #endif
#endif

typedef double real_t;
typedef double2 real2_t;
typedef double4 real4_t;
#define REAL_MIN DBL_MIN

real_t max_in_vector(real4_t vector)
{
	return max(max(vector.x, vector.y), max(vector.z, vector.w));
}

// Gaussian
real_t f(real_t x, real_t a, real_t b, real_t c)
{
	return a * exp( -0.5 * pown((x - b), 2) / pown(c, 2));
}

// dfda
real_t dfda(real_t x, real_t b, real_t c)
{
	return exp( -0.5 * pown((x - b), 2) / pown(c, 2));
}

// dfdb
real_t dfdb(real_t x, real_t a, real_t b, real_t c)
{
	return ((a * (x - b)) / pown(c, 2)) * exp( -0.5 * pown((x - b), 2) / pown(c, 2));
}

// dfdc
real_t dfdc(real_t x, real_t a, real_t b, real_t c)
{
	return ((a * pown((x - b), 2)) / pown(c, 3)) * exp( -0.5 * pown((x - b), 2) / pown(c, 2));
}

//
// Dopasowuje do podanych wartości (y'ów i x'ów) gaussiana (funkcje Gaussa),
// tj. znajduję współczynniki a, b i c funkcji 
//
// f(x) = a * exp( - 0.5 * (x - b)^2 / c^2)
//
// , dla których funkcja najlepiej przybliża zadane wartości.
//
// Dane są ułożone kolumnami w macierzach ys i xs.
//
// Kernel korzysta z metody Levenberg–Marquardt dla nieliniowego problemu
// najmniejszych kwadratów. 
// Wikipedia: http://en.wikipedia.org/wiki/Levenberg%E2%80%93Marquardt_algorithm
//
// Potrzebne definicje
#define MAX_ALFA 1.0e+10
#define MIN_ALFA 1.0e-7
#define ALFA_CHANGE 10
//
__kernel void fit_gaussian
	(
		__global double * ys,	// Wartości, które dopasowujemy
		__global double * xs,	// Argumenty (x'y) dla których dopasowujemy
		uint width,
		__constant uint * cols_sizes, 	// Tablica ilośći znaczących elememntów 
					   	// z kolejnych kolumn macierzy ys
		uint max_iterations,
		__global double4 * results	// Wyniki dopasowania,
						// na początku znajdują się tam początkowe
						// wartości dla współczynników a, b i c.
	)
{
	uint gid = get_global_id(0);		
	if(gid >= width)
	{
		return;
	}

	// W metodzie GN iteracyjnie dfrom vectorochodzimy do "najlepszych"
	// współczynników rozwiązując za każdym razem równanie
	// (J^T * J + alfa * diag(JT * J)) * dB = J^T * R  dla dB
	//
	// (Można rozwiązać równienie dB = (J^T * J + alfa * diag(JT * J))^-1 * J^T * R
	// ale odwracanie macierzy jest problematyczne, a tu są tylko
	// 3 parametry więc łatwiej użyć wzorów Cramera.)
	//
	// B - współczynniki f(x)
	// dB - zmiana wartości współczynników.
	// R - różnice y - f(x, B) 
	// J - jakobian dla funkcji f(x)

	// J^T * J
	real_t JTJ[3][3];
	// Współczynniki
   	real4_t B;
	// J^T * R	
    	real_t JTR[3]; 

	// Początkowe wartości współczyników
	// uzyskane metodą "wyszukanego zgadnięcia" lub zdrowego rozsądku.
	{		
		//#if defined(cl_khr_fp64)
			//B = convert_double4(results[gid]);
		//#else
			B = results[gid];
		//#endif
	}
	// Flaga, że dopasowanie nie doszło do skutku
	B.w = 0.0;

	// Indeks
	uint idx = gid;
	uint col_height = cols_sizes[gid];
	if(col_height == 0)
	{
		return;
	}
	uint idx_max = idx + col_height * width;

	//
	real_t alfa = 100;

	// Dopasowanie dla aktualnych współczynników
	real_t chisq_B = 0.0;
	// Dopasowanie dla nowych współczynników
	real_t chisq_newB = 0.0;
	
	for(uint i = 0; i < max_iterations; i++)
	{		
		// Zerowanie danych
		JTJ[0][0] = 0;
		JTJ[0][1] = 0;
		JTJ[0][2] = 0;
		JTJ[1][0] = 0;
		JTJ[1][1] = 0;
		JTJ[1][2] = 0;
		JTJ[2][0] = 0;
		JTJ[2][1] = 0;
		JTJ[2][2] = 0;

		JTR[0] = 0;
		JTR[1] = 0;
		JTR[2] = 0;

		chisq_B = 0.0;

		idx = gid;
		while(idx < idx_max)
		{
			real_t x = xs[idx];

			// dfda
			real_t dfda_x = dfda(x, B.y, B.z);
			//jacobian[idx] = dfd_;
			// dfdb
			real_t dfdb_x = dfdb(x, B.x, B.y, B.z);
			//jacobian[idx + width] = dfd_;
			// dfdc
			real_t dfdc_x = dfdc(x, B.x, B.y, B.z);
			//jacobian[idx + 2 * width] = dfd_

			real_t y = ys[idx];
		
        		JTJ[0][0] += dfda_x * dfda_x;
        		JTJ[0][1] += dfda_x * dfdb_x;
        		JTJ[0][2] += dfda_x * dfdc_x;

        		JTJ[1][0] += dfdb_x * dfda_x;
        		JTJ[1][1] += dfdb_x * dfdb_x;
        		JTJ[1][2] += dfdb_x * dfdc_x;

        		JTJ[2][0] += dfdc_x * dfda_x;
        		JTJ[2][1] += dfdc_x * dfdb_x;
        		JTJ[2][2] += dfdc_x * dfdc_x;

       			// R[idx]
			real_t r = y - f(x, B.x, B.y, B.z);
			// chisq_B
			chisq_B += pown(r, 2);	
			//JT * R
        		JTR[0] += dfda_x * r;
        		JTR[1] += dfdb_x * r;
        		JTR[2] += dfdc_x * r;   

			idx += width;  
		}
		chisq_B /= 2.0;

		real_t diagJTJ[3];
		diagJTJ[0] = JTJ[0][0];
		diagJTJ[1] = JTJ[1][1];
		diagJTJ[2] = JTJ[2][2];

		// Metoda największego spadku w LM
		// (modyfikacja względem algorytmu Gaussa Newtona)
		JTJ[0][0] += alfa * diagJTJ[0];
		JTJ[1][1] += alfa * diagJTJ[1];
		JTJ[2][2] += alfa * diagJTJ[2];

    		// (JT * J + alfa * diag(JT * J) ) * dB = JT * R jest równaniem typu Ax = b
    		// A = (JT * J), dB = x, JT * R = b
    		// Rozwiązanie za pomocą wzorów Cramera, wikipedia: http://en.wikipedia.org/wiki/Cramer%27s_rule
    		// x_i = det(A_i)/det(A)

    		real_t detA = 
   		JTJ[0][0] * (JTJ[1][1] * JTJ[2][2] - JTJ[1][2] * JTJ[2][1]) -
    		JTJ[0][1] * (JTJ[1][0] * JTJ[2][2] - JTJ[1][2] * JTJ[2][0]) + 
    		JTJ[0][2] * (JTJ[1][0] * JTJ[2][1] - JTJ[1][1] * JTJ[2][0]) ;

    		real_t detA1 =
    		JTR[0]	  * (JTJ[1][1] * JTJ[2][2] - JTJ[1][2] * JTJ[2][1]) -
    		JTJ[0][1] * (  JTR[1]  * JTJ[2][2] - JTJ[1][2] * JTR[2]   ) + 
    		JTJ[0][2] * (  JTR[1]  * JTJ[2][1] - JTJ[1][1] * JTR[2]   ) ;

    		real_t detA2 = 
    		JTJ[0][0] * (JTR[1]    * JTJ[2][2] - JTJ[1][2] * JTR[2]   ) -
    		JTR[0]	  * (JTJ[1][0] * JTJ[2][2] - JTJ[1][2] * JTJ[2][0]) + 
    		JTJ[0][2] * (JTJ[1][0] * JTR[2]    - JTR[1]    * JTJ[2][0]) ;

    		real_t detA3 = 
    		JTJ[0][0] * (JTJ[1][1] * JTR[2]    - JTR[1]    * JTJ[2][1]) -
    		JTJ[0][1] * (JTJ[1][0] * JTR[2]    - JTR[1]    * JTJ[2][0]) + 
    		JTR[0]	  * (JTJ[1][0] * JTJ[2][1] - JTJ[1][1] * JTJ[2][0]) ;

		if(fabs(detA) < REAL_MIN)
		{			
			break;
		}				

		// Zmiana i sprawdzenie warunków stopu
		{    		
			real4_t dB = (real4_t)(detA1/detA, detA2/detA, detA3/detA, 0.0f);
			
			// B(k+1) = B(k) + dB	
			real4_t newB = B + dB;

					
			// Obliczenie dopasowanie dla nowych współczynników
			// jeżeli następuje pierwsza zmiana.
			chisq_newB = 0.0;						
			idx = gid;			
			while(idx < idx_max)
			{
				real_t x = xs[idx];
				real_t fx = f(x, newB.x, newB.y, newB.z);
				real_t y = ys[idx];
			
				// 
				chisq_newB += pown(y - fx, 2);
				idx += width;  
			}
			chisq_newB /= 2.0;
			
			// Sprawdzenie, czy nowe współczynniki są lepsze					
			if(chisq_newB < chisq_B)
			{
				// B(k+1) = B(k)+ dB	
    				B = newB;

				// Modyfikacja w stronę metody Gaussa-Newtonwa
				alfa = max(alfa/ALFA_CHANGE, MIN_ALFA);	
			}
			else
			{
				// Zwiększamy udział metody największego spadku
				// aż dojdziemy do maksymalnego wpływu.
				while(alfa != MAX_ALFA && i < max_iterations)
				{
					i++;

					// Modyfikacja w stronę metody największego spadku
					alfa = min(alfa*ALFA_CHANGE, MAX_ALFA);	

					// Metoda największego spadku w LM
					// (modyfikacja względem algorytmu Gaussa Newtona)
					JTJ[0][0] += (alfa - alfa/ALFA_CHANGE) * diagJTJ[0];
					JTJ[1][1] += (alfa - alfa/ALFA_CHANGE) * diagJTJ[1];
					JTJ[2][2] += (alfa - alfa/ALFA_CHANGE) * diagJTJ[2];

    					detA = 
   					JTJ[0][0] * (JTJ[1][1] * JTJ[2][2] - JTJ[1][2] * JTJ[2][1]) -
    					JTJ[0][1] * (JTJ[1][0] * JTJ[2][2] - JTJ[1][2] * JTJ[2][0]) + 
    					JTJ[0][2] * (JTJ[1][0] * JTJ[2][1] - JTJ[1][1] * JTJ[2][0]) ;

    					detA1 =
    					JTR[0]	  * (JTJ[1][1] * JTJ[2][2] - JTJ[1][2] * JTJ[2][1]) -
    					JTJ[0][1] * (  JTR[1]  * JTJ[2][2] - JTJ[1][2] * JTR[2]   ) + 
    					JTJ[0][2] * (  JTR[1]  * JTJ[2][1] - JTJ[1][1] * JTR[2]   ) ;

    					detA2 = 
    					JTJ[0][0] * (JTR[1]    * JTJ[2][2] - JTJ[1][2] * JTR[2]   ) -
    					JTR[0]	  * (JTJ[1][0] * JTJ[2][2] - JTJ[1][2] * JTJ[2][0]) + 
    					JTJ[0][2] * (JTJ[1][0] * JTR[2]    - JTR[1]    * JTJ[2][0]) ;

    					detA3 = 
    					JTJ[0][0] * (JTJ[1][1] * JTR[2]    - JTR[1]    * JTJ[2][1]) -
    					JTJ[0][1] * (JTJ[1][0] * JTR[2]    - JTR[1]    * JTJ[2][0]) + 
    					JTR[0]	  * (JTJ[1][0] * JTJ[2][1] - JTJ[1][1] * JTJ[2][0]) ;

					if(fabs(detA) < REAL_MIN)
					{			
						break;
					}	

					dB = (real4_t)(detA1/detA, detA2/detA, detA3/detA, 0.0);
			
					// B(k+1) = B(k) + dB	
					newB = B + dB;

					
					// Obliczenie dopasowanie dla nowych współczynników
					// jeżeli następuje pierwsza zmiana.
					chisq_newB = 0.0;						
					idx = gid;			
					while(idx < idx_max)
					{
						real_t x = xs[idx];
						real_t fx = f(x, newB.x, newB.y, newB.z);
						real_t y = ys[idx];
			
						// 
						chisq_newB += pown(y - fx, 2);
						idx += width;  
					}
					chisq_newB /= 2.0;	

					if(chisq_newB < chisq_B)
					{
						// B(k+1) = B(k)+ dB	
    						B = newB;

						// Modyfikacja w stronę metody Gaussa-Newtonwa
						alfa = max(alfa/ALFA_CHANGE, MIN_ALFA);	
						break;
					}					
				}
				
				// Nie udało się osiągnąć lepszego wyniku dla
				// największego alfa, więc kończymy cały algorytm.
				if(alfa == MAX_ALFA)
				{
					break;
				}									
			}
			
		}
	};

	// Zapisanie flagi, ze dopasowanie doszło do skutku.
	B.w = 1.0;
	// Zapisanie wyników	
	//results[gid] = convert_double4(B);
	results[gid] = B;
}

//
// Oblicza funkcje Gaussa, gdzie xs jest macierzą z argumentami funkcji.
//
__kernel void calc_gaussian
	(
		__global double * xs,	// Argumenty dla których obliczamy gaussiana
		__global double4 * gparams,	// Współczynniki funkcji Gaussa
		uint width,		
		__constant uint * cols_sizes, 	// Tablica ilośći znaczących elememntów 
					   	// z kolejnych kolumn macierzy xs		
		__global double * ys	// Wyniki funkcji
	)
{
	uint gid0 = get_global_id(0);		
	if(gid0 >= width)
	{
		return;
	}

	// Ilość elementów w kolumnie
	uint col_height = cols_sizes[gid0];

	uint gid1 = get_global_id(1);	
	if(gid1 >= col_height)
	{
		return;
	}

	// Indeks
	uint idx = gid0 + width * gid1;	
	
	// Pobranie x
	real_t x = xs[idx];

	// Pobranie parametrów dla tej kolumny x'ów
	__local double4 abc_local;
	real4_t abc;
	if(gid1 == 0)
	{
		abc_local = gparams[gid0];
	}
	barrier(CLK_LOCAL_MEM_FENCE);
	#if defined(cl_khr_fp64)
		abc = convert_double4(abc_local);
	#else
		abc = abc_local;
	#endif

	// Obliczenie
	real_t fx = abc.x * exp( (real_t)(-0.5) * pown((x - abc.y), 2) / pown(abc.z, 2));

	// Zapis
	ys[idx] = convert_double(fx);
}

//
// 
//
__kernel void calc_gaussian_chisq
	(
		__global double * xs,	// Argumenty dla których obliczamy gaussiana
		__global double * ys,	// Argumenty dla których obliczamy gaussiana
		__global double * errors,	// 
		__global double4 * gparams,	// Współczynniki funkcji Gaussa
		uint width,		
		__constant uint * cols_sizes, 	// Tablica ilośći znaczących elememntów 
					   	// z kolejnych kolumn macierzy xs		
		__global double * chisqs	// Wyniki chi kwadrat
	)
{
	uint gid = get_global_id(0);		
	if(gid >= width)
	{
		return;
	}

	// Ilość elementów w kolumnie
	uint col_height = cols_sizes[gid];

	// Indeks
	uint idx = gid;	
	uint idx_max = gid + col_height * width;

	real4_t abc;
	#if defined(cl_khr_fp64)
		abc = convert_double4(gparams[gid]);
	#else
		abc = gparams[gid];
	#endif
	
	real_t chisq = (real_t)(0.0);
	real_t c = (real_t)(0.0);
	while (idx < idx_max) 
	{
		real_t f, y, e, t, u, x;

		y = ys[idx];
		x = xs[idx];
		f = abc.x * exp( (real_t)(-0.5) * pown((x - abc.y), 2) / pown(abc.z, 2));
		e = errors[idx];

		u = pown((y-f),2) / pown(e, 2) - c;
		t = chisq + u;
		c = (t - chisq) - u;

		chisq = t;

		idx += width;
	}

	// Zapis
	chisqs[gid] = convert_double(chisq);
}

#define C (real_t)(299792458.0)

//
// 
//
__kernel void calc_gaussian_fwhm
	(
		__global double4 * gparams,	// Współczynniki funkcji Gaussa
		__global double * fwhms,	// 
		uint width
	)
{
	uint idx = get_global_id(0);		
	if(idx >= width)
	{
		return;
	}

	real4_t abc;
	#if defined(cl_khr_fp64)
		abc = convert_double4(gparams[idx]);
	#else
		abc = gparams[idx];
	#endif

	real_t b = abc.y;
	real_t c = abc.z;

	// przeliczenie jednoski współczynnika c (sigma): A -> km/s
	real_t c_kms = C;
 	c_kms *= (pown(((real_t)(1.0) + c/b), 2) - (real_t)(1.0));
	c_kms /= (pown(((real_t)(1.0) + c/b), 2) + (real_t)(1.0));
	c_kms *= (real_t)(1.0e-3);
	
	// c (sigma) -> FWHM
	real_t fwhm = c_kms * (real_t)(2.0);
	fwhm *= pow( ( (real_t)(2.0) * log((real_t)(2.0)) ), (real_t)(0.5) );

/*
	// Tak jest w skrypcie, ale z komentarzy w tym samym skrypcie wynika
	// , że to jest błęda wersja.

	c = c * (real_t)(2.0);
	c *= pow( ( (real_t)(2.0) * log((real_t)(2.0)) ), (real_t)(0.5) );

	real_t fwhm = C;
 	fwhm *= (pown(((real_t)(1.0) + c/b), 2) - (real_t)(1.0));
	fwhm /= (pown(((real_t)(1.0) + c/b), 2) + (real_t)(1.0));
	fwhm *= (real_t)(1.0e-3);
*/	

	// Zapis
	fwhms[idx] = convert_double(fwhm);
}

#if defined(cl_khr_fp64)
	#pragma OPENCL EXTENSION cl_khr_fp64 : disable
#endif
