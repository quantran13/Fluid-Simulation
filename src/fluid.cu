#include <fluid.h>

#define SERIAL_SET_BND 0
#define SERIAL_ADVECT 0
#define SERIAL_PROJECT 0
#define SERIAL_LIN_SOLVE 0

FluidCube *FluidCubeCreate(int size, int diffusion, int viscosity, double dt)
{
    FluidCube *cube = (FluidCube *) malloc(sizeof(*cube));
    size_t N = (size_t) size;

    cube->size = size;
    cube->dt = dt;
    cube->diff = diffusion;
    cube->visc = viscosity;

    cudaMallocManaged((void **) &cube->s, N * N * N * sizeof(double));
    cudaMallocManaged((void **) &cube->density, N * N * N * sizeof(double));

    cudaMallocManaged((void **) &cube->Vx, N * N * N * sizeof(double));
    cudaMallocManaged((void **) &cube->Vy, N * N * N * sizeof(double));
    cudaMallocManaged((void **) &cube->Vz, N * N * N * sizeof(double));

    cudaMallocManaged((void **) &cube->Vx0, N * N * N * sizeof(double));
    cudaMallocManaged((void **) &cube->Vy0, N * N * N * sizeof(double));
    cudaMallocManaged((void **) &cube->Vz0, N * N * N * sizeof(double));

    return cube;
}

void FluidCubeFree(FluidCube *cube)
{
    cudaFree(cube->s);
    cudaFree(cube->density);

    cudaFree(cube->Vx);
    cudaFree(cube->Vy);
    cudaFree(cube->Vz);

    cudaFree(cube->Vx0);
    cudaFree(cube->Vy0);
    cudaFree(cube->Vz0);

    free(cube);
}

static void set_bnd_serial(int b, double *x, int N)
{
    for(int j = 1; j < N - 1; j++) {
        for(int i = 1; i < N - 1; i++) {
            x[IX(i, j, 0  )] = b == 3 ? -x[IX(i, j, 1  )] : x[IX(i, j, 1  )];
            x[IX(i, j, N-1)] = b == 3 ? -x[IX(i, j, N-2)] : x[IX(i, j, N-2)];
        }
    }
    for(int k = 1; k < N - 1; k++) {
        for(int i = 1; i < N - 1; i++) {
            x[IX(i, 0  , k)] = b == 2 ? -x[IX(i, 1  , k)] : x[IX(i, 1  , k)];
            x[IX(i, N-1, k)] = b == 2 ? -x[IX(i, N-2, k)] : x[IX(i, N-2, k)];
        }
    }
    for(int k = 1; k < N - 1; k++) {
        for(int j = 1; j < N - 1; j++) {
            x[IX(0  , j, k)] = b == 1 ? -x[IX(1  , j, k)] : x[IX(1  , j, k)];
            x[IX(N-1, j, k)] = b == 1 ? -x[IX(N-2, j, k)] : x[IX(N-2, j, k)];
        }
    }

    x[IX(0, 0, 0)]       = 0.33f * (x[IX(1, 0, 0)]
                                    + x[IX(0, 1, 0)]
                                    + x[IX(0, 0, 1)]);
    x[IX(0, N-1, 0)]     = 0.33f * (x[IX(1, N-1, 0)]
                                    + x[IX(0, N-2, 0)]
                                    + x[IX(0, N-1, 1)]);
    x[IX(0, 0, N-1)]     = 0.33f * (x[IX(1, 0, N-1)]
                                    + x[IX(0, 1, N-1)]
                                    + x[IX(0, 0, N)]);
    x[IX(0, N-1, N-1)]   = 0.33f * (x[IX(1, N-1, N-1)]
                                    + x[IX(0, N-2, N-1)]
                                    + x[IX(0, N-1, N-2)]);
    x[IX(N-1, 0, 0)]     = 0.33f * (x[IX(N-2, 0, 0)]
                                    + x[IX(N-1, 1, 0)]
                                    + x[IX(N-1, 0, 1)]);
    x[IX(N-1, N-1, 0)]   = 0.33f * (x[IX(N-2, N-1, 0)]
                                    + x[IX(N-1, N-2, 0)]
                                    + x[IX(N-1, N-1, 1)]);
    x[IX(N-1, 0, N-1)]   = 0.33f * (x[IX(N-2, 0, N-1)]
                                    + x[IX(N-1, 1, N-1)]
                                    + x[IX(N-1, 0, N-2)]);
    x[IX(N-1, N-1, N-1)] = 0.33f * (x[IX(N-2, N-1, N-1)]
                                    + x[IX(N-1, N-2, N-1)]
                                    + x[IX(N-1, N-1, N-2)]);
}

static void set_bnd(int b, double *x, int N)
{
#if SERIAL_SET_BND
    for(int j = 1; j < N - 1; j++) {
        for(int i = 1; i < N - 1; i++) {
            x[IX(i, j, 0  )] = b == 3 ? -x[IX(i, j, 1  )] : x[IX(i, j, 1  )];
            x[IX(i, j, N-1)] = b == 3 ? -x[IX(i, j, N-2)] : x[IX(i, j, N-2)];
        }
    }
    for(int k = 1; k < N - 1; k++) {
        for(int i = 1; i < N - 1; i++) {
            x[IX(i, 0  , k)] = b == 2 ? -x[IX(i, 1  , k)] : x[IX(i, 1  , k)];
            x[IX(i, N-1, k)] = b == 2 ? -x[IX(i, N-2, k)] : x[IX(i, N-2, k)];
        }
    }
    for(int k = 1; k < N - 1; k++) {
        for(int j = 1; j < N - 1; j++) {
            x[IX(0  , j, k)] = b == 1 ? -x[IX(1  , j, k)] : x[IX(1  , j, k)];
            x[IX(N-1, j, k)] = b == 1 ? -x[IX(N-2, j, k)] : x[IX(N-2, j, k)];
        }
    }

    x[IX(0, 0, 0)]       = 0.33f * (x[IX(1, 0, 0)]
                                  + x[IX(0, 1, 0)]
                                  + x[IX(0, 0, 1)]);
    x[IX(0, N-1, 0)]     = 0.33f * (x[IX(1, N-1, 0)]
                                  + x[IX(0, N-2, 0)]
                                  + x[IX(0, N-1, 1)]);
    x[IX(0, 0, N-1)]     = 0.33f * (x[IX(1, 0, N-1)]
                                  + x[IX(0, 1, N-1)]
                                  + x[IX(0, 0, N)]);
    x[IX(0, N-1, N-1)]   = 0.33f * (x[IX(1, N-1, N-1)]
                                  + x[IX(0, N-2, N-1)]
                                  + x[IX(0, N-1, N-2)]);
    x[IX(N-1, 0, 0)]     = 0.33f * (x[IX(N-2, 0, 0)]
                                  + x[IX(N-1, 1, 0)]
                                  + x[IX(N-1, 0, 1)]);
    x[IX(N-1, N-1, 0)]   = 0.33f * (x[IX(N-2, N-1, 0)]
                                  + x[IX(N-1, N-2, 0)]
                                  + x[IX(N-1, N-1, 1)]);
    x[IX(N-1, 0, N-1)]   = 0.33f * (x[IX(N-2, 0, N-1)]
                                  + x[IX(N-1, 1, N-1)]
                                  + x[IX(N-1, 0, N-2)]);
    x[IX(N-1, N-1, N-1)] = 0.33f * (x[IX(N-2, N-1, N-1)]
                                  + x[IX(N-1, N-2, N-1)]
                                  + x[IX(N-1, N-1, N-2)]);
#endif

#if not SERIAL_SET_BND
    set_bnd_kernel1 <<< N-2, N-2 >>> (b, x, N);
    set_bnd_kernel2 <<< 1, 1 >>> (x, N);
    cudaDeviceSynchronize();
#endif
}

static void lin_solve(int b, double *x, double *x0, double a, double c, int N)
{
    double cRecip = 1.0 / c;
    int iter = 4;

#if not SERIAL_LIN_SOLVE
    double *x_next;
    cudaMallocManaged((void **) &x_next, N * N * N * sizeof(double));

    for (int k = 0; k < iter; k++) {
        for (int m = 1; m < N - 1; m++) {
            lin_solve_kernel <<< N-2, N-2 >>> (x_next, x, x0, a, cRecip, N, m);
        }

        for (int m = 1; m < N - 1; m++) {
            set_values_kernel <<< N-2, N-2 >>> (x_next, x, m, N);
        }

        set_bnd(b, x, N);
    }

    cudaFree(x_next);
#endif

#if SERIAL_LIN_SOLVE
    for (int k = 0; k < iter; k++) {
        for (int m = 1; m < N - 1; m++) {
            for (int j = 1; j < N - 1; j++) {
                for (int i = 1; i < N - 1; i++) {
                    x[IX(i, j, m)] =
                            (x0[IX(i, j, m)]
                             + a * (x[IX(i + 1, j, m)]
                                    + x[IX(i - 1, j, m)]
                                    + x[IX(i, j + 1, m)]
                                    + x[IX(i, j - 1, m)]
                                    + x[IX(i, j, m + 1)]
                                    + x[IX(i, j, m - 1)]
                            )) * cRecip;
                }
            }
        }

        set_bnd_serial(b, x, N);
        //set_bnd(b, x, N);
    }
#endif
}

static void diffuse(int b, double *x, double *x0, double diff, double dt, int N)
{
    double a = dt * diff * (N - 2) * (N - 2);
    lin_solve(b, x, x0, a, 1 + 6 * a, N);
}

static void advect(int b, double *d, double *d0, double *velocX,
                   double *velocY, double *velocZ, double dt, int N)
{
#if SERIAL_ADVECT
    double i0, i1, j0, j1, k0, k1;

    double dtx = dt * (N - 2);
    double dty = dt * (N - 2);
    double dtz = dt * (N - 2);

    double s0, s1, t0, t1, u0, u1;
    double tmp1, tmp2, tmp3, x, y, z;

    double Ndouble = N;
    double idouble, jdouble, kdouble;
    int i, j, k;

    for(k = 1, kdouble = 1; k < N - 1; k++, kdouble++) {
        for(j = 1, jdouble = 1; j < N - 1; j++, jdouble++) {
            for(i = 1, idouble = 1; i < N - 1; i++, idouble++) {
                tmp1 = dtx * velocX[IX(i, j, k)];
                tmp2 = dty * velocY[IX(i, j, k)];
                tmp3 = dtz * velocZ[IX(i, j, k)];
                x    = idouble - tmp1;
                y    = jdouble - tmp2;
                z    = kdouble - tmp3;

                if(x < 0.5f) x = 0.5f;
                if(x > Ndouble + 0.5f) x = Ndouble + 0.5f;
                i0 = floorf(x);
                i1 = i0 + 1.0f;
                if(y < 0.5f) y = 0.5f;
                if(y > Ndouble + 0.5f) y = Ndouble + 0.5f;
                j0 = floorf(y);
                j1 = j0 + 1.0f;
                if(z < 0.5f) z = 0.5f;
                if(z > Ndouble + 0.5f) z = Ndouble + 0.5f;
                k0 = floorf(z);
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
        }
    }
#endif

#if not SERIAL_ADVECT
    for (int k = 1; k < N - 1; k++) {
        advect_kernel <<< N-2, N-2 >>> (d, d0, velocX, velocY, velocZ, dt, N, k);
    }
    cudaDeviceSynchronize();
#endif

    set_bnd(b, d, N);
}

static void project(double *velocX, double *velocY, double *velocZ,
                    double *p, double *div, int N)
{
    double N_recip = 1 / N;
    for (int k = 1; k < N - 1; k++) {
#if not SERIAL_PROJECT
        project_kernel1 <<< N-2, N-2 >>> (velocX, velocY, velocZ, p, div,
                                          N, N_recip, k);
#endif

#if SERIAL_PROJECT
        for (int j = 1; j < N - 1; j++) {
            for (int i = 1; i < N - 1; i++) {
                div[IX(i, j, k)] = -0.5f*(
                        velocX[IX(i+1, j  , k  )]
                       -velocX[IX(i-1, j  , k  )]
                       +velocY[IX(i  , j+1, k  )]
                       -velocY[IX(i  , j-1, k  )]
                       +velocZ[IX(i  , j  , k+1)]
                       -velocZ[IX(i  , j  , k-1)]
                   ) * N_recip;
                p[IX(i, j, k)] = 0;
            }
        }
#endif
    }

#if not SERIAL_PROJECT
    cudaDeviceSynchronize();
#endif

    set_bnd(0, div, N);
    set_bnd(0, p, N);
    lin_solve(0, p, div, 1, 6, N);

    for (int k = 1; k < N - 1; k++) {
#if SERIAL_PROJECT
        for (int j = 1; j < N - 1; j++) {
            for (int i = 1; i < N - 1; i++) {
                velocX[IX(i, j, k)] -= 0.5f * (p[IX(i+1, j, k)]
                                               - p[IX(i-1, j, k)]) * N;
                velocY[IX(i, j, k)] -= 0.5f * (p[IX(i, j+1, k)]
                                               - p[IX(i, j-1, k)]) * N;
                velocZ[IX(i, j, k)] -= 0.5f * (p[IX(i, j, k+1)]
                                               - p[IX(i, j, k-1)]) * N;
            }
        }
#endif

#if not SERIAL_PROJECT
        project_kernel2 <<< N-2, N-2 >>> (velocX, velocY, velocZ, p, N, k);
#endif
    }

#if not SERIAL_PROJECT
    cudaDeviceSynchronize();
#endif

    set_bnd(1, velocX, N);
    set_bnd(2, velocY, N);
    set_bnd(3, velocZ, N);
}

void FluidCubeStep(FluidCube *cube, perf_t *perf_struct)
{
    int N = cube->size;
    double visc = cube->visc;
    double diff = cube->diff;
    double dt = cube->dt;
    double *Vx = cube->Vx;
    double *Vy = cube->Vy;
    double *Vz = cube->Vz;
    double *Vx0 = cube->Vx0;
    double *Vy0 = cube->Vy0;
    double *Vz0 = cube->Vz0;
    double *s = cube->s;
    double *density = cube->density;

    double start = 0, end = 0;

    start = get_time();
    diffuse(1, Vx0, Vx, visc, dt, N);
    end = get_time();
    perf_struct->timeDiffuse += end - start;

    start = get_time();
    diffuse(2, Vy0, Vy, visc, dt, N);
    end = get_time();
    perf_struct->timeDiffuse += end - start;

    start = get_time();
    diffuse(3, Vz0, Vz, visc, dt, N);
    end = get_time();
    perf_struct->timeDiffuse += end - start;

    start = get_time();
    project(Vx0, Vy0, Vz0, Vx, Vy, N);
    end = get_time();
    perf_struct->timeProject += end - start;

    start = get_time();
    advect(1, Vx, Vx0, Vx0, Vy0, Vz0, dt, N);
    end = get_time();
    perf_struct->timeAdvect += end - start;

    start = get_time();
    advect(2, Vy, Vy0, Vx0, Vy0, Vz0, dt, N);
    end = get_time();
    perf_struct->timeAdvect += end - start;

    start = get_time();
    advect(3, Vz, Vz0, Vx0, Vy0, Vz0, dt, N);
    end = get_time();
    perf_struct->timeAdvect += end - start;

    start = get_time();
    project(Vx, Vy, Vz, Vx0, Vy0, N);
    end = get_time();
    perf_struct->timeProject += end - start;

    start = get_time();
    diffuse(0, s, density, diff, dt, N);
    end = get_time();
    perf_struct->timeDiffuse += end - start;

    start = get_time();
    advect(0, density, s, Vx, Vy, Vz, dt, N);
    end = get_time();
    perf_struct->timeAdvect += end - start;

    perf_struct->totalDiffuse += 4;
    perf_struct->totalAdvect += 4;
    perf_struct->totalProject += 2;
}

void FluidCubeAddDensity(FluidCube *cube, int x, int y, int z, double amount)
{
    int N = cube->size;
    cube->density[IX(x, y, z)] += amount;
}

void FluidCubeAddVelocity(FluidCube *cube, int x, int y, int z,
                          double amountX, double amountY, double amountZ)
{
    int N = cube->size;
    int index = IX(x, y, z);

    cube->Vx[index] += amountX;
    cube->Vy[index] += amountY;
    cube->Vz[index] += amountZ;
}
