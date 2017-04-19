#include <fluid.h>
#include <graphic.h>
#include <string.h>
#include <getopt.h>

const int initial_density = 2000;
const int initial_velocity = 10;

void print_usage(char *program_name);
FluidCube* get_input_from_file(char *file_name);
void get_command_line_args(int argc, char **argv, int *steps, int *write_result,
                           int *display_graphic, char *file_name);
char *get_result_file_name(char *input_file);
void print_result(FluidCube *cube, FILE *result_file, int print_full);
char *double_to_string(double x);
int is_boundary(int i, int j, int k, int N);

int main(int argc, char **argv)
{
    if (argc != 9)
        print_usage(argv[0]);

    // Get command line arguments
    int steps, display_graphic, write_result;
    char file_name[128];
    get_command_line_args(argc, argv, &steps, &write_result, &display_graphic,
                          file_name);

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

    // Delete the output file's content
    char *result_file_name = get_result_file_name(file_name);
    FILE *result_file = NULL;
    if (write_result) {
        result_file = fopen(result_file_name, "w");
        if (result_file == NULL) {
            fprintf(stderr, "WARNING: cannot open result file '%s' to write!",
                    result_file_name);
        }
    }

    // Set up the performance info struct
    perf_t perf_struct;
    perf_struct.timeDrawSquare = 0;
    perf_struct.timeDrawing = 0;
    perf_struct.timeDiffuse = 0;
    perf_struct.timeAdvect = 0;
    perf_struct.timeProject = 0;
    perf_struct.totalDiffuse = 0;
    perf_struct.totalAdvect = 0;
    perf_struct.totalProject = 0;

    // Start the simulation
    double start = get_time();

    printf("---------- starting  ----------\n");
    if (display_graphic)
        draw_cube(cube, &perf_struct);

    for (int step = 0; step < steps; step++) {
        FluidCubeStep(cube, &perf_struct);
        print_result(cube, result_file, cube->size <= 66);

        if (display_graphic)
            draw_cube(cube, &perf_struct);
        printf("---------- done step ----------\n");
    }

    printf("---------- finished  ----------\n\n");
    double end = get_time();
    double elapsed = end - start;

    // End the simulation and print the profiling result
    FluidCubeFree(cube);
    free(result_file_name);
    if (result_file != NULL)
        fclose(result_file);

    printf("Elapsed time: %fs.\n", elapsed);
    printf("Average - diffuse: %f\n", perf_struct.timeDiffuse / perf_struct.totalDiffuse);
    printf("Average - advect : %f\n", perf_struct.timeAdvect / perf_struct.totalAdvect);
    printf("Average - project: %f\n", perf_struct.timeProject / perf_struct.totalProject);
    printf("Average - drawing: %f\n", perf_struct.timeDrawing / (steps + 1));
    printf("Average - draw square: %f\n", perf_struct.timeDrawSquare / (steps + 1));

    return 0;
}

void get_command_line_args(int argc, char **argv, int *steps, int *write_result,
                           int *display_graphic, char *file_name)
{
    char c;
    while ((c = (char) getopt(argc, argv, "s:g:f:r:")) != -1) {
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
            case 'r':
                *write_result = atoi(optarg);
                break;
            default:
                print_usage(argv[0]);
        }
    }
}

void print_usage(char *program_name)
{
    printf("Usage: %s -s <steps> -g <display graphic> -f <file name>"
                   " -r <write result>\n\n", program_name);
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

char *get_result_file_name(char *input_file)
{
    char *result = (char *) malloc(sizeof(char) * 152);
    int i = 0;
    for (i = (int) strlen(input_file) - 1; i >= 0; i--)
        if (input_file[i] == '/') break;

    char *file_name = input_file + i + 1;
    strcpy(result, "results/result_");
    strcat(result, file_name);
    return result;
}

void print_result(FluidCube *cube, FILE *result_file, int print_full)
{
    if (result_file == NULL)
        return;

    int N = cube->size;
    char space[2] = " ";
    char newline[2] = "\n";
    char linebreak[6] = "----\n";

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                if (print_full || is_boundary(i, j, k, N)) {
                    char *num = double_to_string(cube->s[IX(k, j, i)]);
                    fputs(num, result_file);
                    fputs(space, result_file);
                    free(num);

                    num = double_to_string(cube->density[IX(k, j, i)]);
                    fputs(num, result_file);
                    fputs(space, result_file);
                    free(num);

                    num = double_to_string(cube->Vx[IX(k, j, i)]);
                    fputs(num, result_file);
                    fputs(space, result_file);
                    free(num);

                    num = double_to_string(cube->Vy[IX(k, j, i)]);
                    fputs(num, result_file);
                    fputs(space, result_file);
                    free(num);

                    num = double_to_string(cube->Vz[IX(k, j, i)]);
                    fputs(num, result_file);
                    fputs(space, result_file);
                    free(num);

                    fputs(newline, result_file);
                }
            }
        }
    }

    fputs(linebreak, result_file);
}

int is_boundary(int i, int j, int k, int N)
{
    int bool1 = i > 1 && i < N - 2;
    int bool2 = j > 1 && j < N - 2;
    int bool3 = k > 1 && k < N - 2;
    int count = 0;

    if (bool1) count++;
    if (bool2) count++;
    if (bool3) count++;

    return count <= 1;
}

char *double_to_string(double x)
{
    char *num = (char *) malloc(sizeof(char) * 9);
    sprintf(num, "%.4f", x);

    return num;
}
