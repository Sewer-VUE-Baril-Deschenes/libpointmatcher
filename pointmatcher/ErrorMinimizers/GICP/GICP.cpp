/*************************************************************
  Generalized-ICP Copyright (c) 2009 Aleksandr Segal.
  All rights reserved.

  Redistribution and use in source and binary forms, with
  or without modification, are permitted provided that the
  following conditions are met:

* Redistributions of source code must retain the above
  copyright notice, this list of conditions and the
  following disclaimer.
* Redistributions in binary form must reproduce the above
  copyright notice, this list of conditions and the
  following disclaimer in the documentation and/or other
  materials provided with the distribution.
* The names of the contributors may not be used to endorse
  or promote products derived from this software
  without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
  CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
  WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
  INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
  OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
  DAMAGE.
*************************************************************/

#include "GICP.h"
#include "transform.h"
#include "optimize.h"

template<typename T>
GICPErrorMinimizer<T>::GICPErrorMinimizer(const Parameters& params)
{
	// initialize some stuff
}

template<typename T>
GICPErrorMinimizer<T>::GICPErrorMinimizer(const ParametersDoc paramsDoc, const Parameters& params)
{
	// initialize some stuff
}

template<typename T>
typename PointMatcher<T>::TransformationParameters GICPErrorMinimizer<T>::compute(const ErrorElements& mPtsIn)
{
	ErrorElements mPts = mPtsIn;
	int num_matches = 0;
	int n = mPts.reading.getNbPoints();
	double delta = 0.;
	dgc_transform_t t;
	dgc_transform_t t_last;
	dgc_transform_t identity; // TODO: check if initialization is fine
	dgc_transform_identity(identity);
	bool debug_ = false;
	double query_point[3];

	dgc::gicp::gicp_mat_t* mahalanobis = new dgc::gicp::gicp_mat_t[n];

	gsl_matrix* gsl_R = gsl_matrix_alloc(3, 3);
	gsl_matrix* gsl_temp = gsl_matrix_alloc(3, 3);

	bool converged = false;
	int iteration = 0;
	bool opt_status = false;


	/* set up the optimization parameters */
	dgc::gicp::GICPOptData<T> opt_data;
	opt_data.nn_indecies = &mPts.matches.ids;
	opt_data.p1 = &mPts.reading;
	opt_data.p2 = &mPts.reference;
	opt_data.M = mahalanobis;
	opt_data.solve_rotation = true;
	dgc_transform_copy(opt_data.base_t, identity);

	dgc::gicp::GICPOptimizer<T> opt;
	opt.SetDebug(debug_);
	opt.SetMaxIterations(1); // TODO: make sure 1 is the right value
	/* set up the mahalanobis matricies */
	/* these are identity for now to ease debugging */
	for(int i = 0; i < n; i++)
	{
		for(int k = 0; k < 3; k++)
		{
			for(int l = 0; l < 3; l++)
			{
				mahalanobis[i][k][l] = (k == l) ? 1 : 0.;
			}
		}
	}

	while(!converged)
	{
		dgc_transform_t transform_R;
		dgc_transform_copy(transform_R, identity);
		dgc_transform_left_multiply(transform_R, t);
		// copy the rotation component of the current total transformation (including base), into a gsl matrix
		for(int i = 0; i < 3; i++)
		{
			for(int j = 0; j < 3; j++)
			{
				gsl_matrix_set(gsl_R, i, j, transform_R[i][j]);
			}
		}
		/* find correpondences */
		num_matches = 0;
		for(int i = 0; i < n; i++)
		{
			query_point[0] = mPts.reading.features(0, i);
			query_point[1] = mPts.reading.features(1, i);
			query_point[2] = mPts.reading.features(2, i);

			dgc_transform_point(&query_point[0], &query_point[1], &query_point[2], identity);
			dgc_transform_point(&query_point[0], &query_point[1], &query_point[2], t);

			if(nn_dist_sq < max_d_sq)
			{
				// set up the updated mahalanobis matrix here
				gsl_matrix_view C1 = gsl_matrix_view_array(&scan->point_[i].C[0][0], 3, 3);
				gsl_matrix_view C2 = gsl_matrix_view_array(&point_[nn_indecies[i]].C[0][0], 3, 3);
				gsl_matrix_view M = gsl_matrix_view_array(&mahalanobis[i][0][0], 3, 3);
				gsl_matrix_set_zero(&M.matrix);
				gsl_matrix_set_zero(gsl_temp);

				// M = R*C1  // using M as a temp variable here
				gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., gsl_R, &C1.matrix, 1., &M.matrix);

				// temp = M*R' // move the temp value to 'temp' here
				gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1., &M.matrix, gsl_R, 0., gsl_temp);

				// temp += C2
				gsl_matrix_add(gsl_temp, &C2.matrix);
				// at this point temp = C2 + R*C1*R'

				// now invert temp to get the mahalanobis distance metric for gicp
				// M = temp^-1
				gsl_matrix_set_identity(&M.matrix);
				gsl_linalg_cholesky_decomp(gsl_temp);
				for(int k = 0; k < 3; k++)
				{
					gsl_vector_view row_view = gsl_matrix_row(&M.matrix, k);
					gsl_linalg_cholesky_svx(gsl_temp, &row_view.vector);
				}
				num_matches++;
			}
			else
			{
				nn_indecies[i] = -1; // no match
			}
		}
		opt_data.num_matches = num_matches;

		/* optimize transformation using the current assignment and Mahalanobis metrics*/
		dgc_transform_copy(t_last, t);
		opt_status = opt.Optimize(t, opt_data);

		/* compute the delta from this iteration */
		delta = 0.;
		for(int k = 0; k < 4; k++)
		{
			for(int l = 0; l < 4; l++)
			{
				double ratio = 1;
				if(k < 3 && l < 3)
				{ // rotation part of the transform
					ratio = 1. / epsilon_rot_;
				}
				else
				{
					ratio = 1. / epsilon_;
				}
				double c_delta = ratio * fabs(t_last[k][l] - t[k][l]);

				if(c_delta > delta)
				{
					delta = c_delta;
				}
			}
		}

		/* check convergence */
		iteration++;
		if(iteration >= max_iteration_ || delta < 1)
		{
			converged = true;
		}
	}
	if(mahalanobis != NULL)
	{
		delete[] mahalanobis;
	}
	if(gsl_R != NULL)
	{
		gsl_matrix_free(gsl_R);
	}
	if(gsl_temp != NULL)
	{
		gsl_matrix_free(gsl_temp);
	}

	return iteration;
}

template
struct GICPErrorMinimizer<float>;
template
struct GICPErrorMinimizer<double>;