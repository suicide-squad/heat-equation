void __kernel csr_mult(const int        num_rows,
                      __global int *    Ap,
                      __global int *    Aj,
                      __global double * Au,
                      __global double * u,
                      __global double * un,
                        const int       size_time)
{
    __global double* tmp;
    int row = get_global_id(0);
    for (int i = 0; i < size_time; i++) {
        if(row < num_rows) {
            double sum = 0.;

            const int row_start = Ap[row];
            const int row_end = Ap[row+1];

            int jj = 0;
            for (jj = row_start; jj < row_end; jj++)
                sum += Au[jj] * u[Aj[jj]];

            un[row] = sum;

            // barrier(CLK_GLOBAL_MEM_FENCE);
            tmp = u;
            u = un;
            un = tmp;
        }
    }
}
