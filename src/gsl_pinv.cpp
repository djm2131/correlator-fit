/**********************************************************************
    
    CorrelatorFit (www.github.com/djm2131/correlator-fit)

    Source file: src/gsl_pinv.cpp

    Author: David Murphy (djmurphy@mit.edu)

    Copyright (c) 2019 MIT. All rights reserved.

**********************************************************************/

#include "gsl_pinv.h"

gsl_matrix* moore_penrose_pinv(gsl_matrix* A, const double& svd_cut)
{
  gsl_matrix* V;
  gsl_matrix* Sigma_pinv;
  gsl_matrix* U;
  gsl_matrix* A_pinv;

  gsl_matrix* _tmp_mat = nullptr;

  gsl_vector* _tmp_vec;
  gsl_vector* u;

  double x, cutoff;

  size_t i, j;
  size_t n = A -> size1;
  size_t m = A -> size2;
  bool transposed(false);

  // GSL only handles m <= n, so transpose if necessary
  if(m > n)
  {
    transposed = true;
    _tmp_mat = gsl_matrix_alloc(m, n);
    gsl_matrix_transpose_memcpy(_tmp_mat, A);
    A = _tmp_mat;
    i = m;
    m = n;
    n = i;
  }

  // Compute SVD
  V = gsl_matrix_alloc(m, m);
  u = gsl_vector_alloc(m);
  _tmp_vec = gsl_vector_alloc(m);
  gsl_linalg_SV_decomp(A, V, u, _tmp_vec);
  gsl_vector_free(_tmp_vec);

  // Compute Sigma_pinv
  Sigma_pinv = gsl_matrix_alloc(m, n);
  gsl_matrix_set_zero(Sigma_pinv);
  cutoff = svd_cut * gsl_vector_max(u);

  for(i=0; i<m; i++)
  {
    if(gsl_vector_get(u,i) > cutoff){ x = 1.0 / gsl_vector_get(u,i); }
    else{ x = 0.0; }
    gsl_matrix_set(Sigma_pinv, i, i, x);
  }

  U = gsl_matrix_alloc(n, n);
  gsl_matrix_set_zero(U);
  for(i=0; i<n; i++){
  for(j=0; j<m; j++){
    gsl_matrix_set(U, i, j, gsl_matrix_get(A, i, j));
  }}
  if(_tmp_mat != nullptr){ gsl_matrix_free(_tmp_mat); }

  // Assemble pseudoinverse
  _tmp_mat = gsl_matrix_alloc(m, n);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, V, Sigma_pinv, 0.0, _tmp_mat);

  if(transposed){
    A_pinv = gsl_matrix_alloc(n, m);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, U, _tmp_mat, 0.0, A_pinv);
  } else {
    A_pinv = gsl_matrix_alloc(m, n);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, _tmp_mat, U, 0.0, A_pinv);
  }

  // Clean up allocated memory
  gsl_vector_free(u);
  gsl_matrix_free(_tmp_mat);
  gsl_matrix_free(U);
  gsl_matrix_free(Sigma_pinv);
  gsl_matrix_free(V);

  return A_pinv;
}

void gsl_pinv(std::vector<std::vector<double>>& M, const double& svd_cut)
{
  gsl_matrix* A = gsl_matrix_alloc(M.size(), M[0].size());

  for(int i=0; i<M.size(); i++){
  for(int j=0; j<M[0].size(); j++){
    gsl_matrix_set(A, i, j, M[i][j]);
  }}

  gsl_matrix* A_pinv = moore_penrose_pinv(A, svd_cut);

  for(int i=0; i<M.size(); i++){
  for(int j=0; j<M[0].size(); j++){
    M[i][j] = gsl_matrix_get(A_pinv, i, j);
  }}

  gsl_matrix_free(A);
  gsl_matrix_free(A_pinv);
}