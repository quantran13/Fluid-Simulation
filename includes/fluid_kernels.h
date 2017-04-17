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

#endif //FLUID_SIMULATION_FLUID_KERNELS_H
