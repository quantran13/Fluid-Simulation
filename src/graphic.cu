//
// Created by quan on 4/11/17.
//

#include "graphic.h"

void init_render()
{
    glClear(GL_COLOR_BUFFER_BIT);
    glClearColor(1.0, 1.0, 1.0, 1.0);
    glPointSize(3.0);

    glutSwapBuffers();
}

void draw_cube(FluidCube *cube, perf_t *perf_struct)
{
    double start_drawing = get_time();

    int N = cube->size;
    double vsize = 2.0 / N;
    double hsize = vsize;
    double start = 0, end = 0;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            color_t color;
            double density = cube->density[IX(i, j, 25)];
            color.r = density;
            color.g = density;
            color.b = density;

            start = get_time();
            draw_square(i, j, vsize, hsize, &color);
            end = get_time();
            perf_struct->timeDrawSquare += end - start;
            //printf("%d %d: %f\n", i, j, cube->density[IX(i, j, 8)]);
        }
    }

    double end_drawing = get_time();
    perf_struct->timeDrawing += end_drawing - start_drawing;
}

void draw_square(int x, int y, double vsize, double hsize, color_t* color)
{
    double left_x = -1.0 + (double) x * hsize;
    double right_x = -1.0 + (double) (x + 1) * hsize;
    double up_y = 1.0 - (double) y * vsize;
    double down_y = 1.0 - (double) (y + 1) * vsize;

    glColor3d(color->r, color->g, color->b);
    glBegin(GL_POLYGON);
    glVertex2d(left_x, up_y);
    glVertex2d(right_x, up_y);
    glVertex2d(right_x, down_y);
    glVertex2d(left_x, down_y);
    glVertex2d(left_x, up_y);
    glEnd();

    glutSwapBuffers();
}
