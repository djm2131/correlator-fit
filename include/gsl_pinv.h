/**********************************************************************
    
    CorrelatorFit (www.github.com/djm2131/correlator-fit)

    Source file: include/gsl_pinv.h

    Author: David Murphy (djmurphy@mit.edu)

    Copyright (c) 2019 MIT. All rights reserved.

**********************************************************************/

// Uses GSL to compute Moore-Penrose pseudoinverse of a matrix.
// Adapted from gist.github.com/turingbirds/5e99656e08dbe1324c99/moore_penrose_pseudoinverse.c

#ifndef __GSL_PINV_H_INCLUDED__
#define __GSL_PINV_H_INCLUDED__

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <vector>

gsl_matrix* moore_penrose_pinv(gsl_matrix* A, const double& svd_cut);

void gsl_pinv(std::vector<std::vector<double>>& M, const double& svd_cut);

#endif 
