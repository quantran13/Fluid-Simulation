#include <GL/glut.h>
#include <stdio.h>
#include <sys/time.h>
#include "fluid.h"
#include "graphic.h"

double get_time();
void printUsage(char *program_name);

int main(int argc, char **argv)
{
    if (argc != 4)
        printUsage(argv[0]);

    int n = atoi(argv[1]);
    int steps = atoi(argv[2]);
    int display_graphic = atoi(argv[3]);

    // Init graphic
    int width = 1000;
    int height = 1000;

    if (display_graphic) {
        glutInit(&argc, argv);
        glutInitWindowPosition(300, 0);
        glutInitWindowSize(width, height);
        glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
        glutCreateWindow("Fluid Simulation");
        glutDisplayFunc(init_render);
    }

    // Init the cube
    FluidCube* cube = FluidCubeCreate(n, 1, 1, 1);

    for (int i = n/4; i < 3*n/4; i++) {
        for (int j = n/4; j < 3*n/4; j++) {
            for (int k = n/4; k < 3*n/4; k++) {
                FluidCubeAddDensity(cube, i, j, k, 2000);
                FluidCubeAddVelocity(cube, i, j, k, 100, 100, 100);
            }
        }
    }

    /*
    for (int z = 0; z < n; z++) {
        FluidCubeAddDensity(cube, n-2, n-2, z, 200);
        //FluidCubeAddVelocity(cube, 0, 0, z, 100, 100, 100);
    }

    for (int x = 0; x < n-1; x++) {
        for (int y = 0; y < n-1; y++) {
            for (int z = 0; z < n-1; z++) {
                FluidCubeAddVelocity(cube, x, y, z, 100, 100, 100);
            }
        }
    }
     */

    // Start the simulation
    double start = get_time();
    for (int step = 0; step < steps; step++) {
        FluidCubeStep(cube);
        draw_cube(cube);
        printf("---------------------- done step -----------------------\n");
    }

    FluidCubeFree(cube);

    double end = get_time();
    double elapsed = end - start;

    printf("Elapsed time: %fs.\n", elapsed);

    return 0;
}

void printUsage(char *program_name)
{
    printf("Usage: %s <n> <steps> <display graphic>\n\n", program_name);
    exit(0);
}

double get_time()
{
    struct timeval t;
    double retval;

    gettimeofday(&t, NULL);
    retval = t.tv_sec;
    retval += t.tv_usec / 1000000.0;
    return retval;
}

