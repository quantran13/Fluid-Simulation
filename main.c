#include <stdio.h>
#include <sys/time.h>
#include "fluid.h"

double get_time()
{
    struct timeval t;
    double retval;
    
    gettimeofday(&t, NULL);
    retval = t.tv_sec;
    retval += t.tv_usec / 1000000.0;
    return retval;
}

int main(int argc, char **argv)
{
    int n = 100;
    double start = get_time();

    FluidCube* cube = FluidCubeCreate(n, 1, 1, 1);
    FluidCubeStep(cube);
    FluidCubeFree(cube);

    double end = get_time();
    double elapsed = end - start;
    printf("Elapsed time: %fs.\n", elapsed);

    return 0;
}
