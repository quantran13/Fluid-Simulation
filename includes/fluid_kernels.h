//
// Created by quan on 4/16/17.
//

#ifndef FLUID_SIMULATION_FLUID_KERNELS_H
#define FLUID_SIMULATION_FLUID_KERNELS_H

#include <math.h>
#include <utility.h>
#include <stdio.h>

__global__ void advect_kernel(double *d, double *d0, double *vecX, double *velocY,
                              double *lo, double dt, int N, int k);
__global__ void set_bnd_kernel(int b, double *x, int N);
__global__ void project_kernel(double *velocX, double *velocY, double *velocZ,
                               double *p, double *div, int iter, int N,
                               double N_recip, int k);

#endif //FLUID_SIMULATION_FLUID_KERNELS_H
