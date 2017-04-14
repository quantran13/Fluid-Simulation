#include <fluid.h>
#include <utility.h>

FluidCube *FluidCubeCreate(int size, int diffusion, int viscosity, double dt)
{
    FluidCube *cube = malloc(sizeof(*cube));
    size_t N = (size_t) size;
    
    cube->size = size;
    cube->dt = dt;
    cube->diff = diffusion;
    cube->visc = viscosity;
    
    cube->s = calloc(N * N * N, sizeof(double));
    cube->density = calloc(N * N * N, sizeof(double));
    
    cube->Vx = calloc(N * N * N, sizeof(double));
    cube->Vy = calloc(N * N * N, sizeof(double));
    cube->Vz = calloc(N * N * N, sizeof(double));
    
    cube->Vx0 = calloc(N * N * N, sizeof(double));
    cube->Vy0 = calloc(N * N * N, sizeof(double));
    cube->Vz0 = calloc(N * N * N, sizeof(double));
    
    return cube;
}

void FluidCubeFree(FluidCube *cube)
{
    free(cube->s);
    free(cube->density);
    
    free(cube->Vx);
    free(cube->Vy);
    free(cube->Vz);
    
    free(cube->Vx0);
    free(cube->Vy0);
    free(cube->Vz0);
    
    free(cube);
}

static void set_bnd(int b, double *x, int N)
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

static void lin_solve(int b, double *x, double *x0, double a, double c, int iter, int N)
{
    double cRecip = 1.0 / c;
    for (int k = 0; k < iter; k++) {
        for (int m = 1; m < N - 1; m++) {
            for (int j = 1; j < N - 1; j++) {
                for (int i = 1; i < N - 1; i++) {
                    x[IX(i, j, m)] =
                        (x0[IX(i, j, m)]
                            + a*(    x[IX(i+1, j  , m  )]
                                    +x[IX(i-1, j  , m  )]
                                    +x[IX(i  , j+1, m  )]
                                    +x[IX(i  , j-1, m  )]
                                    +x[IX(i  , j  , m+1)]
                                    +x[IX(i  , j  , m-1)]
                           )) * cRecip;
                }
            }
        }
        set_bnd(b, x, N);
    }
}

static void diffuse (int b, double *x, double *x0, double diff, double dt, int iter, int N)
{
    double a = dt * diff * (N - 2) * (N - 2);
    lin_solve(b, x, x0, a, 1 + 6 * a, iter, N);
}

static void advect(int b, double *d, double *d0,  double *velocX,
                   double *velocY, double *velocZ, double dt, int N)
{
    double i0, i1, j0, j1, k0, k1;
    
    double dtx = dt * (N - 2);
    double dty = dt * (N - 2);
    double dtz = dt * (N - 2);
    

    double Ndouble = N;
    double idouble, jdouble, kdouble;
    int i, j, k;

    for(k = 1, kdouble = 1; k < N - 1; k++, kdouble++) {
        for(j = 1, jdouble = 1; j < N - 1; j++, jdouble++) {
            for(i = 1, idouble = 1; i < N - 1; i++, idouble++) {
                double s0, s1, t0, t1, u0, u1;
                double tmp1, tmp2, tmp3, x, y, z;
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
        }
    }
    set_bnd(b, d, N);
}

static void project(double *velocX, double *velocY, double *velocZ,
                    double *p, double *div, int iter, int N)
{
    double N_recip = 1 / N;
    for (int k = 1; k < N - 1; k++) {
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
    }
    set_bnd(0, div, N); 
    set_bnd(0, p, N);
    lin_solve(0, p, div, 1, 6, iter, N);

    for (int k = 1; k < N - 1; k++) {
        for (int j = 1; j < N - 1; j++) {
            for (int i = 1; i < N - 1; i++) {
                velocX[IX(i, j, k)] -= 0.5f * (  p[IX(i+1, j, k)]
                                                -p[IX(i-1, j, k)]) * N;
                velocY[IX(i, j, k)] -= 0.5f * (  p[IX(i, j+1, k)]
                                                -p[IX(i, j-1, k)]) * N;
                velocZ[IX(i, j, k)] -= 0.5f * (  p[IX(i, j, k+1)]
                                                -p[IX(i, j, k-1)]) * N;
            }
        }
    }
    set_bnd(1, velocX, N);
    set_bnd(2, velocY, N);
    set_bnd(3, velocZ, N);
}

void FluidCubeStep(FluidCube *cube, perf_t *perf_struct)
{
    int N          = cube->size;
    double visc     = cube->visc;
    double diff     = cube->diff;
    double dt       = cube->dt;
    double *Vx      = cube->Vx;
    double *Vy      = cube->Vy;
    double *Vz      = cube->Vz;
    double *Vx0     = cube->Vx0;
    double *Vy0     = cube->Vy0;
    double *Vz0     = cube->Vz0;
    double *s       = cube->s;
    double *density = cube->density;

    double start = 0, end = 0;

    start = get_time();
    diffuse(1, Vx0, Vx, visc, dt, 4, N);
    end = get_time();
    perf_struct->timeDiffuse += end - start;

    start = get_time();
    diffuse(2, Vy0, Vy, visc, dt, 4, N);
    end = get_time();
    perf_struct->timeDiffuse += end - start;

    start = get_time();
    diffuse(3, Vz0, Vz, visc, dt, 4, N);
    end = get_time();
    perf_struct->timeDiffuse += end - start;

    start = get_time();
    project(Vx0, Vy0, Vz0, Vx, Vy, 4, N);
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
    project(Vx, Vy, Vz, Vx0, Vy0, 4, N);
    end = get_time();
    perf_struct->timeProject += end - start;

    start = get_time();
    diffuse(0, s, density, diff, dt, 4, N);
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

