//
// Created by quan on 4/16/17.
//

#ifndef FLUID_SIMULATION_FLUID_KERNELS_H
#define FLUID_SIMULATION_FLUID_KERNELS_H

__device__ void advect_kernel(double *d, double *d0, double *vecX, double *vecY,
                              double *vecZ, double dt, int N);

#endif //FLUID_SIMULATION_FLUID_KERNELS_H
