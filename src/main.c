#include <fluid.h>
#include <graphic.h>
#include <string.h>
#include <getopt.h>
#include <utility.h>

const int initial_density = 2000;
const int initial_velocity = 10;

void print_usage(char *program_name);
FluidCube* get_input_from_file(char *file_name);
void get_command_line_args(int argc, char **argv, int *steps,
                           int *display_graphic, char *file_name);

int main(int argc, char **argv)
{
    if (argc != 7)
        print_usage(argv[0]);

    // Get command line arguments
    int steps, display_graphic;
    char file_name[128];
    get_command_line_args(argc, argv, &steps, &display_graphic, file_name);

    // Init graphic
    int width = 1000;
    int height = 1000;

    glutInit(&argc, argv);
    glutInitWindowPosition(300, 0);
    glutInitWindowSize(width, height);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
    glutCreateWindow("Fluid Simulation");
    glutDisplayFunc(init_render);

    // Init the cube
    FluidCube* cube = get_input_from_file(file_name);
    int n = cube->size;

    // Start the simulation
    perf_t perf_struct;
    perf_struct.timeDrawSquare = 0;
    perf_struct.timeDrawing = 0;
    perf_struct.timeDiffuse = 0;
    perf_struct.timeAdvect = 0;
    perf_struct.timeProject = 0;
    perf_struct.totalDiffuse = 0;
    perf_struct.totalAdvect = 0;
    perf_struct.totalProject = 0;

    double start = get_time();

    if (display_graphic)
        draw_cube(cube, &perf_struct);
    for (int step = 0; step < steps; step++) {
        printf("---------- done step -----------\n");
        FluidCubeStep(cube, &perf_struct);

        if (display_graphic)
            draw_cube(cube, &perf_struct);
    }

    FluidCubeFree(cube);

    double end = get_time();
    double elapsed = end - start;

    printf("Elapsed time: %fs.\n", elapsed);
    printf("Average - diffuse: %f\n", perf_struct.timeDiffuse / perf_struct.totalDiffuse);
    printf("Average - advect : %f\n", perf_struct.timeAdvect / perf_struct.totalAdvect);
    printf("Average - project: %f\n", perf_struct.timeProject / perf_struct.totalProject);
    printf("Average - drawing: %f\n", perf_struct.timeDrawing / (steps + 1));
    printf("Average - draw square: %f\n", perf_struct.timeDrawSquare / (steps + 1));

    return 0;
}

void get_command_line_args(int argc, char **argv, int *steps,
                           int *display_graphic, char *file_name)
{
    char c;
    while ((c = (char) getopt(argc, argv, "s:g:f:")) != -1) {
        switch (c) {
            case 's':
                *steps = atoi(optarg);
                break;
            case 'g':
                *display_graphic = atoi(optarg);
                break;
            case 'f':
                strncpy(file_name, optarg, strlen(optarg));
                break;
            default:
                print_usage(argv[0]);
        }
    }
}

void print_usage(char *program_name)
{
    printf("Usage: %s -s <steps> -g <display graphic> -f <file name>\n\n", program_name);
    exit(0);
}

FluidCube* get_input_from_file(char *file_name)
{
    FluidCube *cube;
    FILE *fin = fopen(file_name, "r");
    char n_str[10], diffusion_str[10], viscosity_str[10];
    int n, diffusion, viscosity;

    fgets(n_str, 9, fin);
    fgets(diffusion_str, 10, fin);
    fgets(viscosity_str, 10, fin);
    n = atoi(n_str);
    diffusion = atoi(diffusion_str);
    viscosity = atoi(viscosity_str);
    cube = FluidCubeCreate(n, diffusion, viscosity, 1);

    char line[64];
    for (int x = 0; x < n; x++) {
        for (int y = 0; y < n; y++) {
            for (int z = 0; z < n; z++) {
                fgets(line, 63, fin);
                char c = line[strlen(line) - 2];
                if (c == '1') {
                    FluidCubeAddDensity(cube, x, y, z, initial_density);
                    FluidCubeAddVelocity(cube, x, y, z, initial_velocity,
                                         initial_velocity, initial_velocity);
                }
            }
        }
    }

    fclose(fin);

    return cube;
}

