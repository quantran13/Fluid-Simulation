//
// Created by quan on 4/16/17.
//

#include <fluid_kernels.h>
#include <device_launch_parameters.h>

__global__ void advect_kernel(double *d, double *d0, double *velocX, double *velocY,
                              double *velocZ, double dt, int N, int k)
{
    double Ndouble = (double) N;
    double dtx = dt * (N - 2);
    double dty = dt * (N - 2);
    double dtz = dt * (N - 2);

    int j = blockIdx.x;
    int i = threadIdx.x;
    double idouble = (double) i;
    double jdouble = (double) j;
    double kdouble = (double) k;

    double s0, s1, t0, t1, u0, u1;
    double tmp1, tmp2, tmp3, x, y, z;
    double i0, i1, j0, j1, k0, k1;

    tmp1 = dtx * velocX[IX(i, j, k)];
    tmp2 = dty * velocY[IX(i, j, k)];
    tmp3 = dtz * velocZ[IX(i, j, k)];
    x    = idouble - tmp1;
    y    = jdouble - tmp2;
    z    = kdouble - tmp3;

    if(x < 0.5f) x = 0.5f;
    if(x > Ndouble + 0.5f) x = Ndouble + 0.5f;
    i0 = floor(x);
    i1 = i0 + 1.0f;
    if(y < 0.5f) y = 0.5f;
    if(y > Ndouble + 0.5f) y = Ndouble + 0.5f;
    j0 = floor(y);
    j1 = j0 + 1.0f;
    if(z < 0.5f) z = 0.5f;
    if(z > Ndouble + 0.5f) z = Ndouble + 0.5f;
    k0 = floor(z);
    k1 = k0 + 1.0f;

    s1 = x - i0;
    s0 = 1.0f - s1;
    t1 = y - j0;
    t0 = 1.0f - t1;
    u1 = z - k0;
    u0 = 1.0f - u1;

    int i0i = (int) i0;
    int i1i = (int) i1;
    int j0i = (int) j0;
    int j1i = (int) j1;
    int k0i = (int) k0;
    int k1i = (int) k1;

    d[IX(i, j, k)] =

            s0 * ( t0 * (u0 * d0[IX(i0i, j0i, k0i)]
                         +u1 * d0[IX(i0i, j0i, k1i)])
                   +( t1 * (u0 * d0[IX(i0i, j1i, k0i)]
                            +u1 * d0[IX(i0i, j1i, k1i)])))
            +s1 * ( t0 * (u0 * d0[IX(i1i, j0i, k0i)]
                          +u1 * d0[IX(i1i, j0i, k1i)])
                    +( t1 * (u0 * d0[IX(i1i, j1i, k0i)]
                             +u1 * d0[IX(i1i, j1i, k1i)])));
}
