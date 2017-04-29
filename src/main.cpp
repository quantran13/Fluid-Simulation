#include <fluid.h>
#include <graphic.h>
#include <string.h>
#include <getopt.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct {
    perf_t *perf;
    FluidCube** cubes;
    int steps;
} info_t;

// MPI
#define MAX_N 1000000
#define MAX_N_PROCS 100
#define DATA_MSG 0
#define NEWDATA_MSG 1

// Fluid stuff
const int initial_density = 2000;
const int initial_velocity = 10;
const int print_full_threshold = 66;

int nnodes;
int n;
int me;
int list[MAX_N];
int count[MAX_N_PROCS];

void init(int argc, char **argv);
void managernode(int argc, char **argv);
void workernode();
void print_usage(char *program_name);
FluidCube* get_input_from_file(char *file_name);
FluidCube** get_input_parallel(char *file_name, int p);
void get_command_line_args(int argc, char **argv, int *steps, int *write_result,
                           int *display_graphic, char *file_name);
char *get_result_file_name(char *input_file);
void print_result(FluidCube *cube, FILE *result_file, int print_full);
char *double_to_string(double x);
int is_boundary(int i, int j, int k, int N);
FluidCube* combineCubes(FluidCube** cubes);

int main(int argc, char **argv)
{
  init(argc, argv);
  if(me == 0){
    managernode(argc, argv);
  }
  else{
    workernode();
  }
  MPI_Finalize();
  return 0;
}

void managernode(int argc, char **argv){
    MPI_Status status;
    int i;

    if(argc != 9)
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
    FluidCube** cubes;
    if (cubes == NULL) {
        fprintf(stderr, "Cannot open input file %s!\n", file_name);
        print_usage(argv[0]);
    }

    // Delete the output file's content
    char *result_file_name = get_result_file_name(file_name);
    FILE *result_file = NULL;
    if (write_result) {
        result_file = fopen(result_file_name, "w");
        if (result_file == NULL) {
            fprintf(stderr, "WARNING: cannot open result file '%s' to write!\n\n",
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
    int slice = 20;

    // MPI division of work
    info_t *info;
    info->cubes = cubes;
    info->perf = &perf_struct;
    info->steps = steps;

    for(i = 1; i<nnodes; i++){
      MPI_Send(info, 8, MPI_INT, i, DATA_MSG, MPI_COMM_WORLD);
    }

    FluidCube* myCube1 = cubes[0];
    FluidCube* myCube2 = cubes[1];
    printf("---------- starting  ----------\n");
    //if (display_graphic)
        //draw_cube(cube, &perf_struct, slice);
    //print_result(cube, result_file, cube->size < print_full_threshold);
    //print_result(cube, result_file, cube->size < print_full_threshold);

    for (int step = 0; step < steps; step++) {
        FluidCubeStep(myCube1, &perf_struct);
        //print_result(myCube1, result_file, myCube1->size < print_full_threshold);
        FluidCubeStep(myCube2, &perf_struct);
        //print_result(myCube2, result_file, myCube2->size < print_full_threshold);
        //if (display_graphic)
            //draw_cube(cube, &perf_struct, slice);
        printf("---------- done step ----------\n");
    }
    for(i = 1; i<nnodes; i++){
      MPI_Recv(info, 8, MPI_INT, i, DATA_MSG, MPI_COMM_WORLD, &status);
    }

    FluidCube* finalCube = combineCubes(cubes);
    printf("---------- finished  ----------\n\n");
    double end = get_time();
    double elapsed = end - start;

    print_result(finalCube, result_file, myCube2->size < print_full_threshold);
    if (display_graphic)
        draw_cube(finalCube, &perf_struct, slice);
    // End the simulation and print the profiling result
    FluidCubeFree(finalCube);
    FluidCubeFree(finalCube);
    FluidCubeFree(finalCube);

    free(result_file_name);
    if (result_file != NULL)
        fclose(result_file);

    printf("Elapsed time: %fs.\n", elapsed);
    printf("Average - diffuse: %f\n", perf_struct.timeDiffuse / perf_struct.totalDiffuse);
    printf("Average - advect : %f\n", perf_struct.timeAdvect / perf_struct.totalAdvect);
    printf("Average - project: %f\n", perf_struct.timeProject / perf_struct.totalProject);
    printf("Average - drawing: %f\n", perf_struct.timeDrawing / (steps + 1));
    printf("Average - draw square: %f\n", perf_struct.timeDrawSquare / (steps + 1));

    return;
}

void workernode(){
  MPI_Status status;
  info_t *info;
  MPI_Recv(info, 8, MPI_INT, 0, NEWDATA_MSG, MPI_COMM_WORLD, &status);
  FluidCube* mycube1 = info->cubes[me*2];
  FluidCube* mycube2 = info->cubes[me*2-1];
  int steps = info->steps;
  for (int step = 0; step < steps; step++) {
      FluidCubeStep(mycube1, info->perf);
      //print_result(mycube1, result_file, mycube1->size < print_full_threshold);
      FluidCubeStep(mycube2, info->perf);
      //print_result(mycube2, result_file, mycube2->size < print_full_threshold);
      printf("---------- done step ----------\n");
  }
  info->cubes[me*2] = mycube1;
  info->cubes[me*2+1] = mycube2;
  MPI_Send(info, 8, MPI_INT, 0, NEWDATA_MSG, MPI_COMM_WORLD);
}
FluidCube* combineCubes(FluidCube** cubes){


}


void init(int argc, char **argv){
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nnodes);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
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
                file_name[strlen(optarg)] = 0;
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
FluidCube** get_input_parallel(char *file_name, int p)
{
    FluidCube **cubes;
    FILE *fin = fopen(file_name, "r");
    char n_str[10], diffusion_str[10], viscosity_str[10];
    int n, diffusion, viscosity;

    if (fin == NULL)
        return NULL;

    fgets(n_str, 9, fin);
    fgets(diffusion_str, 10, fin);
    fgets(viscosity_str, 10, fin);
    n = atoi(n_str);
    diffusion = atoi(diffusion_str);
    viscosity = atoi(viscosity_str);
    for(int i=0; i<p; i++){
      cubes[i] = FluidCubeCreate(n/2, diffusion, viscosity, 1);

    }
    char line[5];
    for (int x = 0; x < n; x++) {
        for (int y = 0; y < n; y++) {
            for (int z = 0; z < n; z++) {
                fgets(line, 4, fin);
                char c = line[0];
                int subCube;
                if(x>n/2){
                  subCube += 4;
                }
                if(y>n/2){
                  subCube += 2;
                }
                if(z>n/2){
                  subCube += 1;
                }
                if (c == '1') {
                    FluidCubeAddDensity(cubes[subCube], x, y, z, initial_density);
                    FluidCubeAddVelocity(cubes[subCube], x, y, z, initial_velocity,
                                         initial_velocity, initial_velocity);
                }
            }
        }
    }

    fclose(fin);

    return cubes;
}

FluidCube* get_input_from_file(char *file_name)
{
    FluidCube *cube;
    FILE *fin = fopen(file_name, "r");
    char n_str[10], diffusion_str[10], viscosity_str[10];
    int n, diffusion, viscosity;

    if (fin == NULL)
        return NULL;

    fgets(n_str, 9, fin);
    fgets(diffusion_str, 10, fin);
    fgets(viscosity_str, 10, fin);
    n = atoi(n_str);
    diffusion = atoi(diffusion_str);
    viscosity = atoi(viscosity_str);
    cube = FluidCubeCreate(n, diffusion, viscosity, 1);

    char line[5];
    for (int x = 0; x < n; x++) {
        for (int y = 0; y < n; y++) {
            for (int z = 0; z < n; z++) {
                fgets(line, 4, fin);
                char c = line[0];
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
    strcpy(result, "output/result_");
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
    // round x to 5 decimal points
    long longX = (long) x * 100000;
    x = longX / 100000.0;
    if (x == -0.0) x = 0.0;

    char *num = (char *) malloc(sizeof(char) * 9);
    sprintf(num, "%.5f", x);

    return num;
}
