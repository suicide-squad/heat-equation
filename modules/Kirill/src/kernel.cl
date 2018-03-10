//#pragma OPENCL EXTENSION cl_khr_fp64 : enable


void __kernel csr_mult_d( const int               num_rows,
                          __global int *          Ap,
                          __global int *          Aj,
                          __global const double * Au,
                          __global const double * u,
                          __global double *       un) {
    int row = get_global_id(0);

    double sum = 0.0;

    const int row_start = Ap[row];
    const int row_end = Ap[row+1];

    int jj = 0;
    for (jj = row_start; jj < row_end; jj++)
        sum += Au[jj] * u[Aj[jj]];

    un[row] = sum;
}

#define I(x,y,z) ((x) + (y)*nx + (z)*ny*nx)

void __kernel naive_formula_d(const int               nx,
                              const int               ny,
                              const int               nz,
                              const double            coeff0,
                              const double            coeff1,
                              const double            coeff2,
                              const double            coeff3,
                              __global const double*  u,
                              __global double*        un) {

    const int ox = get_global_id(0);
    const int oy = get_global_id(1);
    const int oz = get_global_id(2);

    int ix = ox;
    int iy = oy;
    int iz = oz;

    if (ix == 0) {
      ix++;
    } else if (ix == nx - 1) {
      ix--;
    }

    if (iy == 0) {
      iy++;
    } else if (iy == ny - 1) {
      iy--;
    }

    if (iz == 0) {
      iz++;
    } else if (iz == nz - 1) {
      iz--;
    }

    un[I(ox, oy, oz)] = coeff3*(u[I(ix, iy, iz - 1)] + u[I(ix, iy, iz + 1)]) + coeff2*(u[I(ix, iy - 1, iz)] +
        u[I(ix, iy + 1, iz)]) + coeff1*(u[I(ix - 1, iy, iz)] + u[I(ix + 1, iy, iz)]) + coeff0*u[I(ix, iy, iz)];
}

void __kernel csr_mult_f(const int     num_rows,
                      __global int *   Ap,
                      __global int *   Aj,
                      __global float * Au,
                      __global float * u,
                      __global float * un)
{
    int row = get_global_id(0);

    if(row < num_rows) {
        float sum = 0.0;

        const int row_start = Ap[row];
        const int row_end = Ap[row+1];

        int jj = 0;
        for (jj = row_start; jj < row_end; jj++)
            sum += Au[jj] * u[Aj[jj]];

        un[row] = sum;
    }
}
