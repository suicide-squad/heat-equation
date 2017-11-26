//#pragma OPENCL EXTENSION cl_khr_fp64 : enable



void __kernel csr_mult_f(const int        num_rows,
                      __global int *    Ap,
                      __global int *    Aj,
                      __global float *  Au,
                      __global float *  u,
                      __global float *  un,
                      const int         size_time)
{
    int row = get_global_id(0);
   // printf("%d\n", row);

    if(row < num_rows) {
        float sum = 0.f;

        const int row_start = Ap[row];
        const int row_end = Ap[row+1];

        int jj = 0;
        for (jj = row_start; jj < row_end; jj++)
            sum += Au[jj] * u[Aj[jj]];

        un[row] = sum;
    }
}