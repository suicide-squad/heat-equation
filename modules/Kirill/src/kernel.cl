void __kernel csr_mult(const int num_rows,
		__global int * Ap,
		__global int * Aj,
		__global double * Ax,
		__global double * x,
		__global double * y)
{
	int row = get_global_id(0);

	if(row < num_rows) {
		double sum = y[row];

		const int row_start = Ap[row];
		const int row_end = Ap[row+1];

		int jj = 0;
		for (jj = row_start; jj < row_end; jj++)
			sum += Ax[jj] * x[Aj[jj]];

		y[row] = sum;
	}
}
