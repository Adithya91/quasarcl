
__kernel void reduce_fe_chisqs
	(
		__global double * chisqs, 		// Bufor z chisq
		__constant uint * filtered_sizes, 	// Ilość znaczących elementów po filtracji
		uint size				// Ilość chisq i filtered_sizes
	)
{
	// gid0 - numer elementu z chisqs (czyli jednego chisq)
	uint gid0 = get_global_id(0);

	if(gid0 >= size)
	{
		return;
	}
		
	double filtered_size, chisq;

	filtered_size = (double)filtered_sizes[gid0];
	chisq = chisqs[gid0];

	chisq /= (filtered_size - 1);

	chisqs[gid0] = chisq;	
}
