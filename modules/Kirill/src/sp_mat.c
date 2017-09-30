#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>
#include <sp_mat.h>

#include "sp_mat.h"

void initSpMat(SpMatrix *mat, int nz, int nRows) {
  mat->nz = nz;
  mat->nRows = nRows;
  mat->value = (TYPE *)aligned_alloc(32, sizeof(TYPE) * nz);
  mat->col = (int *)aligned_alloc(32, sizeof(int) * nz);
  mat->rowIndex = (int *)calloc(nRows + 1, sizeof(int));
}

void freeSpMat(SpMatrix* mat) {
  free(mat->value);
  free(mat->col);
  free(mat->rowIndex);
}

inline void multMV_default(TYPE* result, SpMatrix mat, TYPE* vec) {
  TYPE localSum;
#pragma omp parallel private(localSum) if (ENABLE_PARALLEL)
  {
#pragma omp for nowait
    for (int i = 0; i < mat.nRows; i++) {
      localSum = 0.0;
      for (int j = mat.rowIndex[i]; j < mat.rowIndex[i + 1]; j++)
        localSum += mat.value[j] * vec[mat.col[j]];
      result[i] = localSum;
    }
  }
}

#define IND2(x,y,z) ((x) + (y)*nx + (z)*ny*nx)


inline void multMV_AVX_1(TYPE* result, SpMatrix mat, TYPE* vec, TYPE* coeff) {
  TYPE * val = mat.value;
  TYPE *vec2;
  int nx = 66;
  int ny = 66;
  int nz = 66;

    for (int z = 0; z < nz; z++) {
    for (int y = 0; y < ny; y++) {
      for (int x = 0; x < nx; x++) {
        if (x==0)
          vec[IND2(x,y,z)]=vec[IND2(x+1,y,z)];
        if (x==nx-1)
          vec[IND2(x,y,z)]=vec[IND2(x-1,y,z)];
        if (y==0)
          vec[IND2(x,y,z)]=vec[IND2(x,y+1,z)];
        if (y==nx-1)
          vec[IND2(x,y,z)]=vec[IND2(x,y-1,z)];
        if (z==0)
          vec[IND2(x,y,z)]=vec[IND2(x,y,z+1)];
        if (z==nx-1)
          vec[IND2(x,y,z)]=vec[IND2(x,y,z-1)];
      }
    }
  }

  for (int z = 0; z < nz; z++) {
    if (z == 0 || z == nz - 1) {
      val += nx*ny * NR;  result += nx*ny;
      continue;
    }
    for (int y = 0; y < ny; y++) {
      if (y == 0 || y == ny - 1) {
        val += nx * NR;   result += nx;
        continue;
      }

      for (int x = 0; x < nx; x ++) {
        if (x == 0 || x == nx - 1) {
          val += 1 * NR;   result += 1;
          continue;
        }

        vec2 = &vec[IND2(x,y,z-1)];
        m_ind index4 = mm_set_epi32(0);
        m_real val4 = mm_i32gather(val, index4);
        m_real vec4 = mm_load(vec2);
        m_real sum4 = _mm256_mul_pd(val4, vec4);

        index4 = mm_set_epi32(1);
        val4 = mm_i32gather(val, index4);
        vec2 += nx*ny - nx;
        vec4 = mm_load(vec2);
        sum4 = mm_fmadd(val4, vec4, sum4);

        index4 = mm_set_epi32(2);
        val4 = mm_i32gather(val, index4);
        vec2 += nx -1;
        vec4 = mm_load(vec2);
        sum4 = mm_fmadd(val4, vec4, sum4);

        index4 = mm_set_epi32(3);
        val4 = mm_i32gather(val, index4);
        vec2 ++;
        vec4 = mm_load(vec2);
        sum4 = mm_fmadd(val4, vec4, sum4);

        index4 = mm_set_epi32(4);
        val4 = mm_i32gather(val, index4);
        vec2 ++;
        vec4 = mm_load(vec2);
        sum4 = mm_fmadd(val4, vec4, sum4);

        index4 = mm_set_epi32(5);
        val4 = mm_i32gather(val, index4);
        vec2 += -1 + nx;
        vec4 = mm_load(vec2);
        sum4 = mm_fmadd(val4, vec4, sum4);

        index4 = mm_set_epi32(6);
        val4 = mm_i32gather(val, index4);
        vec2 += -nx + nx*ny;
        vec4 = mm_load(vec2);
        sum4 = mm_fmadd(val4, vec4, sum4);
        mm_stream(result, sum4);

        x+=3;
        result+=4;
        val += 4 * NR;
      }
    }
  }
}



inline void multMV_AVX(TYPE* result, SpMatrix mat, TYPE* vec, TYPE* coeff) {
  TYPE * val = mat.value;
  TYPE *vec2;
  int nx = 66;
  int ny = 66;
  int nz = 66;
//  vec+=67;
  //  for (int z = 0; z < nz; z++) {
//    for (int y = 0; y < ny; y++) {
//      for (int x = 0; x < nx; x++) {
//        if (x==0)
//          vec[IND2(x,y,z)]=vec[IND2(x+1,y,z)];
//        if (x==nx-1)
//          vec[IND2(x,y,z)]=vec[IND2(x-1,y,z)];
//        if (y==0)
//          vec[IND2(x,y,z)]=vec[IND2(x,y+1,z)];
//        if (y==nx-1)
//          vec[IND2(x,y,z)]=vec[IND2(x,y-1,z)];
//        if (z==0)
//          vec[IND2(x,y,z)]=vec[IND2(x,y,z+1)];
//        if (z==nx-1)
//          vec[IND2(x,y,z)]=vec[IND2(x,y,z-1)];
//      }
//    }
//  }

  m_real val41 = _mm256_set1_pd(coeff[0]);
  m_real val42 = _mm256_set1_pd(coeff[1]);
  m_real val43 = _mm256_set1_pd(coeff[2]);
  m_real val44 = _mm256_set1_pd(coeff[3]);

  for (int z = 0; z < nz; z++) {
    if (z == 0 || z == nz - 1) {
      val += nx*ny * NR;  result += nx*ny;
      continue;
    }
//    vec2 = vec + 2*nx;
    for (int y = 0; y < ny; y++) {
      if (y == 0 || y == ny - 1) {
        val += nx * NR;   result += nx;
        continue;
      }
//      vec2 = vec + 12;

      for (int x = 0; x < nx; x ++) {
        if (x == 0 || x == nx - 1) {
          val += 1 * NR;
          result += 1;
          continue;
        }
//              printf("x=%d y=%d z=%d\n", x, y ,z);
//        m_real sum4 = _mm256_setzero_pd();

//        for (int j = 0; j < NR; j++) {
//          m_ind index4 = mm_set_epi32(j);
//          m_real val4 = mm_i32gather(val, index4);

        TYPE * tmp = vec2;
        m_real vec4 = mm_load(tmp);
        m_real sum4 = _mm256_mul_pd(val44, vec4);
        tmp+=(nx)*ny;
        vec4 = mm_load(tmp);
        sum4 = mm_fmadd(val43, vec4, sum4);
        tmp+=nx;
        vec4 = mm_load(tmp);
        sum4 = mm_fmadd(val42, vec4, sum4);
        tmp+=1;
        vec4 = mm_load(tmp);
        sum4 = mm_fmadd(val41, vec4, sum4);
        tmp+=1;
        vec4 = mm_load(tmp);
        sum4 = mm_fmadd(val42, vec4, sum4);
        tmp+=nx;
        vec4 = mm_load(tmp);
        sum4 = mm_fmadd(val43, vec4, sum4);
        tmp+=(nx)*ny;
        vec4 = mm_load(tmp);
        sum4 = mm_fmadd(val44, vec4, sum4);
//        }
        mm_stream(result, sum4);
        x+=3;
//        val +=4*NR;
        result+=4;
        vec2+=4;
      }
    }
  }
}

inline void multMV_AVX_v2(TYPE* result, SpMatrix mat, TYPE* vec) {
//    __m128i mask = _mm_setr_epi32(-1, -1, -1, 0);
//  __m256d mask2 = _mm256_setr_pd(-1,-1,-1,0);
//  __m256d zero =_mm256_setzero_pd();
//  #pragma omp parallel for if (ENABLE_PARALLEL)
  for (int row = 0; row < mat.nRows; row ++) {
    m_ind v_col1 = mm_load_si(mat.col);
    m_real v_val1 = mm_load(mat.value);
    m_real v_vec1 = mm_i32gather(vec, v_col1);

    m_real rowsum = mm_mul(v_vec1, v_val1);

#ifdef DOUBLE_TYPE
//    __m128i v_col2 = _mm_maskload_epi32(mat.col+4, mask);
    m_ind v_col2 = mm_load_si(mat.col+LENVEC);
    m_real v_val2 = mm_load(mat.value+LENVEC);
    m_real v_vec2 = mm_i32gather(vec, v_col2);
//    __m256d v_vec2 = _mm256_mask_i32gather_pd(zero, vec, v_col2, mask2, 8);
    rowsum = mm_fmadd(v_vec2, v_val2, rowsum);
    rowsum = mm_hadd(rowsum, rowsum);
#endif
    result[row] = ((TYPE*)&rowsum)[0] + ((TYPE*)&rowsum)[2];
#ifdef FLOAT_TYPE
    result[row] += ((TYPE*)&rowsum)[4] + ((TYPE*)&rowsum)[6];
#endif
    mat.value += NR;
    mat.col += NR;
  }
}

# define multMV_mkl_f(result, mat, vec)                                                 \
  mkl_cspblas_scsrgemv("N", &mat.nRows, mat.value, mat.rowIndex, mat.col, vec, result); \

# define multMV_mkl_d(result, mat, vec)                                                 \
  mkl_cspblas_dcsrgemv("N", &mat.nRows, mat.value, mat.rowIndex, mat.col, vec, result); \


void multMV(TYPE* result, SpMatrix mat, TYPE* vec, TYPE* coeff) {
#if defined(MKL_RUN)
  #if defined(FLOAT_TYPE)
    multMV_mkl_f(result, mat, vec);
  #elif defined(DOUBLE_TYPE)
    multMV_mkl_d(result, mat, vec);
  #endif
#elif defined(AVX2_RUN)
  multMV_AVX_1(result,mat, vec, coeff);
#else
  multMV_default(result, mat, vec);
#endif
}

void sumV(TYPE **result, TYPE *U, TYPE *k1, TYPE *k2, TYPE *k3, TYPE *k4, int N, TYPE h) {
  #pragma omp parallel for if (ENABLE_PARALLEL)
  for (int i = 0; i < N; i++)
    (*result)[i] = U[i] + h*(k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
}


void printSpMat(SpMatrix mat) {
  for (int i = 0; i < mat.nRows; i++) {
    for (int j = 0; j < mat.nRows; j++)
      printf("%.0lf", procedure(mat, i, j));
    printf("\n");
  }
}

TYPE procedure(SpMatrix mat, int i, int j) {
  TYPE result = 0;
  int N1 = mat.rowIndex[i];
  int N2 = mat.rowIndex[i+1];
  for(int k = N1; k < N2; k++) {
    if (mat.col[k] == j) {
      result = mat.value[k];
      break;
    }
  }
  return result;
}

void denseMult(TYPE **result, TYPE **mat, TYPE *vec, int dim) {
  memset(*result, 0, dim*sizeof(TYPE));
  for (int x = 0; x < dim; x++) {
    for (int i = 0;i < dim;i++)
      (*result)[x]+=mat[x][i]*vec[i];
  }
}

void createExplicitSpMat(SpMatrix *mat, TYPE coeffs[4], int dim, int NX, int NXY) {
  int index = 0, j, k, shiftIndex;
  mat->rowIndex[0] = 0;
  for (int i = 0; i < dim; i++) {
    if (i % NX == 0)
      shiftIndex = 1;
    else if (i % NX == NX - 1)
      shiftIndex = -1;
    else
      shiftIndex = 0;

    mat->col[index] = i + shiftIndex;
    mat->value[index] = coeffs[1];
    index++;

    // ***************************************
    //      Смещение на x - 1
    // ***************************************
    mat->col[index] = i + shiftIndex - 1;
    mat->value[index] = coeffs[0];
    index++;

    // ***************************************
    //      Смещение на x + 1
    // ***************************************
    mat->col[index] = i + shiftIndex + 1;
    mat->value[index] = coeffs[0];
    index++;
    // ***************************************
    //      Смещение на y - 1
    // ***************************************
    j = (i + shiftIndex) % NXY - NX;
    if (j < 0) {
//      mat->col[index] = DIM + j;
      mat->col[index] = NXY + j;
      mat->value[index] = coeffs[2];

      index++;
    } else {
      mat->col[index] = j;
      mat->value[index] = coeffs[2];

      index++;
    }
    // ***************************************
    //      Смещение на y + 1
    // ***************************************
    mat->col[index] = ((i + shiftIndex) / NXY) * NXY + (i + shiftIndex + NX) % NXY;
    //mat->col[index] = (i + shiftIndex + NX) % NXY;
    //mat->col[index] = (i + NX) % NXY;
    mat->value[index] = coeffs[2];

    index++;
    // ***************************************
    //    Смещение на  z - 1
    // ***************************************
    k = i + shiftIndex - NXY;
    if (k <= 0) {
      mat->col[index] = dim + k;
      mat->value[index] = coeffs[3];
      index++;
    } else {
      mat->col[index] = k;
      mat->value[index] = coeffs[3];
      index++;
    }
    // ***************************************
    //    Смещение на  z + 1
    // ***************************************
    mat->col[index] = (i + shiftIndex + NXY) % dim;
    mat->value[index] = coeffs[3];
    index++;
    // ***************************************

//    mat->col[index] = 0;
//    mat->value[index] = 0.0;
//    index++;

    mat->rowIndex[i + 1] = mat->rowIndex[i] + NR;
  }
}

void createExplicitSpMatV2(SpMatrix *mat, TYPE coeffs[4], int nx, int ny, int nz) {
  int index = 0, k = 0;
  int shiftIndexX, shiftIndexY, shiftIndexZ;
  int shift;
  mat->rowIndex[0] = 0;

  int realIndex;
  for (int z = 0; z < nz; z++) {
    if (z == 0)
      shiftIndexZ = nx*ny;
    else if (z == nz - 1)
      shiftIndexZ = -nx*ny;
    else
      shiftIndexZ = 0;
    for (int y = 0; y < ny; y++) {
      if (y == 0)
        shiftIndexY = nx;
      else if (y == ny - 1)
        shiftIndexY = -nx;
      else
        shiftIndexY = 0;

      for (int x = 0; x < nx; x++) {
        if (x == 0)
          shiftIndexX = 1;
        else if (x == nx - 1)
          shiftIndexX = -1;
        else
          shiftIndexX = 0;

        shift = shiftIndexZ + shiftIndexY + shiftIndexX;
        realIndex = x + y*nx + z*nx*ny + shift;

        // ***************************************
        //      Смещение на z - 1
        // ***************************************
        mat->col[index] = realIndex - nx*ny;
        mat->value[index] = coeffs[3];
        index++;

        // ***************************************
        //      Смещение на y - 1
        // ***************************************
        mat->col[index] = realIndex - nx;
        mat->value[index] = coeffs[2];
        index++;

        // ***************************************
        //      Смещение на x - 1
        // ***************************************
        mat->col[index] = realIndex - 1;
        mat->value[index] = coeffs[1];
        index++;

// Отсутствие смещения
        mat->col[index] = realIndex;
        mat->value[index] = coeffs[0];
        index++;


        // ***************************************
        //      Смещение на x + 1
        // ***************************************
        mat->col[index] = realIndex + 1;
        mat->value[index] = coeffs[1];
        index++;


        // ***************************************
        //      Смещение на y + 1
        // ***************************************
        mat->col[index] = realIndex + nx;
        mat->value[index] = coeffs[2];
        index++;

        // ***************************************
        //      Смещение на z + 1
        // ***************************************
        mat->col[index] = realIndex + nx*ny;
        mat->value[index] = coeffs[3];
        index++;

//        mat->col[index] = 0;
//        mat->value[index] = 0.0;
//        index++;

        k++;
        mat->rowIndex[k] = mat->rowIndex[k-1] + NR ;

      }

    }
  }
}

void createExplicitSpMatV2R(SpMatrix* mat, TYPE* coeffs, int nx, int ny, int nz, MPI_Comm comm) {
  int index = 0, k = 0;
  int shiftIndexX=0, shiftIndexY=0, shiftIndexZ=0;
  int shift;
  mat->rowIndex[0] = 0;

  // Определения соседних ранков в декардовой решётке
  int rank_left, rank_right, rank_down, rank_top;
  MPI_Cart_shift(comm, 1, 1, &rank_left, &rank_right);
  MPI_Cart_shift(comm, 0, 1, &rank_down, &rank_top);

  int realIndex;
  for (int z = 0; z < nz; z++) {
    if (z == 0 && rank_down==-2)
      shiftIndexZ = 4*nx*ny;
    else if (z==1&& rank_down==-2)
      shiftIndexZ = 3*nx*ny;
    else if (z==2&& rank_down==-2)
      shiftIndexZ = 2*nx*ny;
    else if (z==3&& rank_down==-2)
      shiftIndexZ = 1*nx*ny;
    else if (z == nz - 1&& rank_top==-2)
      shiftIndexZ = -1*nx*ny;
    else if (z == nz - 2&& rank_top==-2)
      shiftIndexZ = -2*nx*ny;
    else if (z == nz - 3&& rank_top==-2)
      shiftIndexZ = -3*nx*ny;
    else if (z == nz - 4&& rank_top==-2)
      shiftIndexZ = -4*nx*ny;
    else
      shiftIndexZ = 0;
    for (int y = 0; y < ny; y++) {
      if (y == 0&&rank_left==-2)
        shiftIndexY = 4*nx;
      else if (y==1&&rank_left==-2)
        shiftIndexY = 3*nx;
      else if (y==2&&rank_left==-2)
        shiftIndexY = 2*nx;
      else if (y==3&&rank_left==-2)
        shiftIndexY = 1*nx;
      else if (y == ny - 1 &&rank_right==-2)
        shiftIndexY = -1*nx;
      else if (y == ny - 2&&rank_right==-2)
        shiftIndexY = -2*nx;
      else if (y == ny - 3&&rank_right==-2)
        shiftIndexY = -3*nx;
      else if (y == ny - 4&&rank_right==-2)
        shiftIndexY = -4*nx;
      else
        shiftIndexY = 0;

      for (int x = 0; x < nx; x++) {
        if (x == 0)
          shiftIndexX = 4;
        else if (x==1)
          shiftIndexX = 3;
        else if (x==2)
          shiftIndexX = 2;
        else if (x==3)
          shiftIndexX = 1;
        else if (x == nx - 1)
          shiftIndexX = -1;
        else if (x == nx - 2)
          shiftIndexX = -2;
        else if (x == nx - 3)
          shiftIndexX = -3;
        else if (x == nx - 4)
          shiftIndexX = -4;
        else
          shiftIndexX = 0;

        shift = shiftIndexZ + shiftIndexY + shiftIndexX;
        realIndex = x + y*nx + z*nx*ny + shift;

        // ***************************************
        //      Смещение на z - 1
        // ***************************************
        mat->col[index] = realIndex - nx*ny;
        mat->value[index] = coeffs[3];
        index++;

        // ***************************************
        //      Смещение на y - 1
        // ***************************************
        mat->col[index] = realIndex - nx;
        mat->value[index] = coeffs[2];
        index++;

        // ***************************************
        //      Смещение на x - 1
        // ***************************************
        mat->col[index] = realIndex - 1;
        mat->value[index] = coeffs[1];
        index++;

        // Отсутствие смещения
        mat->col[index] = realIndex;
         mat->value[index] = coeffs[0];
        index++;

        // ***************************************
        //      Смещение на x + 1
        // ***************************************
        mat->col[index] = realIndex + 1;
        mat->value[index] = coeffs[1];
        index++;


        // ***************************************
        //      Смещение на y + 1
        // ***************************************
        mat->col[index] = realIndex + nx;
        mat->value[index] = coeffs[2];
        index++;

        // ***************************************
        //      Смещение на z + 1
        // ***************************************
        mat->col[index] = realIndex + nx*ny;
        mat->value[index] = coeffs[3];
        index++;

        k++;
        mat->rowIndex[k] = mat->rowIndex[k-1] + 7;

      }

    }
  }
}

void createImplicitSpMat(SpMatrix* mat, TYPE* coeffs, int nx, int ny, int nz) {
  int index = 0, k = 0;
  int shiftIndexX, shiftIndexY, shiftIndexZ;
  int shift;
  mat->rowIndex[0] = 0;

  int realIndex;
  for (int z = 0; z < nz; z++) {
    if (z == 0)
      shiftIndexZ = nx*ny;
    else if (z == nz - 1)
      shiftIndexZ = -nx*ny;
    else
      shiftIndexZ = 0;
    for (int y = 0; y < ny; y++) {
      if (y == 0)
        shiftIndexY = nx;
      else if (y == ny - 1)
        shiftIndexY = -nx;
      else
        shiftIndexY = 0;

      for (int x = 0; x < nx; x++) {
        if (x == 0)
          shiftIndexX = 1;
        else if (x == nx - 1)
          shiftIndexX = -1;
        else
          shiftIndexX = 0;

        shift = shiftIndexZ + shiftIndexY + shiftIndexX;
        realIndex = x + y*nx + z*nx*ny + shift;

        // ***************************************
        //      Смещение на z - 1
        // ***************************************
        mat->col[index] = realIndex - nx*ny;
        mat->value[index] = coeffs[2];
        index++;

        // ***************************************
        //      Смещение на y - 1
        // ***************************************
        mat->col[index] = realIndex - nx;
        mat->value[index] = coeffs[1];
        index++;

        // ***************************************
        //      Смещение на x - 1
        // ***************************************
        mat->col[index] = realIndex - 1;
        mat->value[index] = coeffs[0];
        index++;

        // ***************************************
        //      Смещение на x + 1
        // ***************************************
        mat->col[index] = realIndex + 1;
        mat->value[index] = coeffs[0];
        index++;


        // ***************************************
        //      Смещение на y + 1
        // ***************************************
        mat->col[index] = realIndex + nx;
        mat->value[index] = coeffs[1];
        index++;

        // ***************************************
        //      Смещение на z + 1
        // ***************************************
        mat->col[index] = realIndex + nx*ny;
        mat->value[index] = coeffs[2];
        index++;

        k++;
        mat->rowIndex[k] = mat->rowIndex[k - 1] + 6;

      }

    }
  }
}