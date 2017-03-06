// MATLAB MEX function for tvl1optflow estimation. The code is provided by
// Javier Sanchez Peres <jsanchez@dis.ulpgc.es>
// The MEX file was created by Mauricio Delbracio <mdelbra@gmail.com>
//
//
// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2011, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// All rights reserved.

#include "mex.h"
#include "matrix.h"


#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

//COMPUTE OPTICAL FLOW IN 1/3 Resolution of the input image
//Then upsample to original resolution
#define PAR_DEFAULT_ZOUTINI 0.3

#define DISABLE_OMP 1

#ifndef DISABLE_OMP
#include <omp.h>
#endif//DISABLE_OMP

#include "tvl1flow_lib.c"

#define PAR_DEFAULT_OUTFLOW "flow.flo"
#define PAR_DEFAULT_NPROC   0
#define PAR_DEFAULT_TAU     0.25
#define PAR_DEFAULT_LAMBDA  0.15
#define PAR_DEFAULT_THETA   0.3
#define PAR_DEFAULT_NSCALES 100
#define PAR_DEFAULT_ZFACTOR 0.5
#define PAR_DEFAULT_NWARPS  5
#define PAR_DEFAULT_EPSILON 0.01
#define PAR_DEFAULT_VERBOSE 0


/**
 *
 *  Main program:
 *   This program reads the following parameters from the console and
 *   then computes the optical flow:
 *   -nprocs      number of threads to use (OpenMP library)
 *   -I0          first image
 *   -I1          second image
 *   -tau         time step in the numerical scheme
 *   -lambda      data term weight parameter
 *   -theta       tightness parameter
 *   -nscales     number of scales in the pyramidal structure
 *   -zfactor     downsampling factor for creating the scales
 *   -nwarps      number of warps per scales
 *   -epsilon     stopping criterion threshold for the iterative process
 *   -out         name of the output flow field
 *   -verbose     switch on/off messages
 *
 */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    //nlhs: number output
    //plhs: array of pointers
    //nrhs: number of input
    //prhs: array of input
    
    
    if (nrhs < 10) {
        mexErrMsgTxt("Incorrect number of parameters");
        return;
    }
    
    //Input Images
    int ndims1 = mxGetNumberOfDimensions(prhs[0]);
    int ndims2 = mxGetNumberOfDimensions(prhs[1]);
    
    if (ndims1 !=2 || ndims2 !=2)
        mexErrMsgTxt("Images must be grayscaled");
    
    // Make sure that image is double precision
    if (mxGetElementSize(prhs[0]) != 4 || mxGetElementSize(prhs[1]) != 4 )
        mexErrMsgTxt("image must be single precision.");
    
    int ny = mxGetM(prhs[0]);
    int nx = mxGetN(prhs[0]);
    
    int ny2 = mxGetM(prhs[1]);
    int nx2 = mxGetN(prhs[1]);
    
    
    if (nx != nx2 && ny != ny2)
    {
        char buffer[400];
        sprintf(buffer,"ERROR: input images size mismatch %dx%d != %dx%d\n", nx, ny, nx2, ny2);
        mexErrMsgTxt(buffer);
    }

    float *I0in = (float*)mxGetData(prhs[0]);
    float *I1in = (float*)mxGetData(prhs[1]);
    
    float *I0 = xmalloc(nx * ny * sizeof(float));
    float *I1 = xmalloc(nx * ny * sizeof(float));

    for (int i = 0; i < ny; i++)
        for (int j = 0; j < nx; j++)
            {
                I0[nx*i + j] =  (float) I0in[i + j*ny];
                I1[nx*i + j] =  (float) I1in[i + j*ny];

            }
    
    
    //read the rest of the parameters
    int i = 2;
    
    //char* outfile = (nrhs>i)? plhs[i]: PAR_DEFAULT_OUTFLOW;       i++;
#ifndef DISABLE_OMP
    int   nproc   = (nrhs>i)? mxGetScalar(prhs[i]): PAR_DEFAULT_NPROC;   i++;
#endif//DISABLE_OMP

    float zfini   = (nrhs>i)? mxGetScalar(prhs[i]): PAR_DEFAULT_ZOUTINI; i++; //Added by mdelbra. 19Jun2015
    float tau     = (nrhs>i)? mxGetScalar(prhs[i]): PAR_DEFAULT_TAU;     i++;
    float lambda  = (nrhs>i)? mxGetScalar(prhs[i]): PAR_DEFAULT_LAMBDA;  i++;
    float theta   = (nrhs>i)? mxGetScalar(prhs[i]): PAR_DEFAULT_THETA;   i++;
    int   nscales = (nrhs>i)? mxGetScalar(prhs[i]): PAR_DEFAULT_NSCALES; i++;
    float zfactor = (nrhs>i)? mxGetScalar(prhs[i]): PAR_DEFAULT_ZFACTOR; i++;
    int   nwarps  = (nrhs>i)? mxGetScalar(prhs[i]): PAR_DEFAULT_NWARPS;  i++;
    float epsilon = (nrhs>i)? mxGetScalar(prhs[i]): PAR_DEFAULT_EPSILON; i++;
    int   verbose = (nrhs>i)? mxGetScalar(prhs[i]): PAR_DEFAULT_VERBOSE; i++;
    
    //	//check parameters
    //	if (nproc < 0) {
    //		nproc = PAR_DEFAULT_NPROC;
    //		if (verbose) mexErrMsgTxt(sprinttf("warning: "
    //				"nproc changed to %d\n", nproc));
    //	}
    //	if (tau <= 0 || tau > 0.25) {
    //		tau = PAR_DEFAULT_TAU;
    //		if (verbose) mexErrMsgTxt(sprintf("warning: "
    //				"tau changed to %g\n", tau);
    //	}
    //	if (lambda <= 0) {
    //		lambda = PAR_DEFAULT_LAMBDA;
    //		if (verbose) fprintf(stderr, "warning: "
    //				"lambda changed to %g\n", lambda);
    //	}
    //	if (theta <= 0) {
    //		theta = PAR_DEFAULT_THETA;
    //		if (verbose) fprintf(stderr, "warning: "
    //				"theta changed to %g\n", theta);
    //	}
    //	if (nscales <= 0) {
    //		nscales = PAR_DEFAULT_NSCALES;
    //		if (verbose) fprintf(stderr, "warning: "
    //				"nscales changed to %d\n", nscales);
    //	}
    //	if (zfactor <= 0 || zfactor >= 1) {
    //		zfactor = PAR_DEFAULT_ZFACTOR;
    //		if (verbose) fprintf(stderr, "warning: "
    //				"zfactor changed to %g\n", zfactor);
    //	}
    //	if (nwarps <= 0) {
    //		nwarps = PAR_DEFAULT_NWARPS;
    //		if (verbose) fprintf(stderr, "warning: "
    //				"nwarps changed to %d\n", nwarps);
    //	}
    //	if (epsilon <= 0) {
    //		epsilon = PAR_DEFAULT_EPSILON;
    //		if (verbose) fprintf(stderr, "warning: "
    //				"epsilon changed to %f\n", epsilon);
    //	}
    
#ifndef DISABLE_OMP
    if (nproc > 0)
        omp_set_num_threads(nproc);
#endif//DISABLE_OMP
    
    //	// read the input images
    //	int    nx, ny, nx2, ny2;
    //	float *I0 = read_image(image1_name, &nx, &ny);
    //	float *I1 = read_image(image2_name, &nx2, &ny2);
    
    //read the images and compute the optical flow
    //Set the number of scales according to the size of the
    //images.  The value N is computed to assure that the smaller
    //images of the pyramid don't have a size smaller than 16x16
    const float N = 1 + log(hypot(nx, ny)/16.0) / log(1/zfactor);
    if (N < nscales)
        nscales = N;
    
    if (verbose)
    {
        char buffer[400];
        sprintf(buffer,"zfini=%f tau=%f lambda=%f theta=%f nscales=%d zfactor=%f nwarps=%d epsilon=%g\n",
                zfini,tau, lambda, theta, nscales,
                zfactor, nwarps, epsilon);
        
        mexPrintf(buffer);
    }
    
    //shuochen: output flow size might not match in the input images due to the disable or upscale
    int nx_out, ny_out;
    zoom_size(nx,ny,&nx_out,&ny_out,zfini);

    //allocate memory for the flow
    float *u = xmalloc(2 * nx_out * ny_out * sizeof*u);
    float *v = u + nx_out*ny_out;
    
    //compute the optical flow
    Dual_TVL1_optic_flow_multiscale(
                                    I0, I1, u, v, nx, ny, zfini, tau, lambda, theta,
                                    nscales, zfactor, nwarps, epsilon, verbose
                                    );
    
    

    int dims[3];
    dims[0] = ny_out;
    dims[1] = nx_out;
    dims[2] = 2;
    
    plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    
    double *oFinal = mxGetPr(plhs[0]);

    for (int i = 0; i < ny_out; i++)
        for (int j = 0; j < nx_out; j++)
            for (int k = 0; k < 2; k++)
            {
                oFinal[ny_out*nx_out*k + ny_out*j + i] =  (double) u[j  + i*nx_out  + nx_out*ny_out*k];
                // oFinal[ny*nx*k + nx*i + j] =  (double) u[j  + i*nx  + nx*ny*k];
                
            }
    
    //delete allocated memory
    free(I0);
    free(I1);
    free(u);
    
    
    return;
}
