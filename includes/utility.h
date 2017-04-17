//
// Created by quan on 4/13/17.
//

#ifndef FLUID_SIMULATION_UTILITY_H
#define FLUID_SIMULATION_UTILITY_H

#include <stdlib.h>
#include <sys/time.h>

#define IX(x, y, z) ((x) + (y) * N + (z) * N * N)

typedef struct {
    double timeDrawSquare;
    double timeDrawing;
    double timeDiffuse;
    double timeProject;
    double timeAdvect;
    int totalDiffuse;
    int totalProject;
    int totalAdvect;
} perf_t;

extern "C"
{
double get_time();
}

#endif //FLUID_SIMULATION_UTILITY_H
