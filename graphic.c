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

void draw_cube(FluidCube *cube, int n)
{
    double vsize = 2.0 / n;
    double hsize = vsize;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            draw_square(i, j, vsize, hsize, 1);
        }
    }
}

void draw_square(int x, int y, double vsize, double hsize, double color)
{
    double left_x = -1.0 + (double) x * hsize;
    double right_x = -1.0 + (double) (x + 1) * hsize;
    double up_y = 1.0 - (double) y * vsize;
    double down_y = 1.0 - (double) (y + 1) * vsize;

    glColor3d(color, color, color);
    glBegin(GL_POLYGON);
    glVertex2d(left_x, up_y);
    glVertex2d(right_x, up_y);
    glVertex2d(right_x, down_y);
    glVertex2d(left_x, down_y);
    glVertex2d(left_x, up_y);
    glEnd();

    glColor3f(0, 0, 0);
    glBegin(GL_LINES);
    glVertex2d(left_x, up_y);
    glVertex2d(left_x, down_y);
    glEnd();
    glBegin(GL_LINES);
    glVertex2d(left_x, up_y);
    glVertex2d(right_x, up_y);
    glEnd();
    glBegin(GL_LINES);
    glVertex2d(right_x, up_y);
    glVertex2d(right_x, down_y);
    glEnd();
    glBegin(GL_LINES);
    glVertex2d(left_x, down_y);
    glVertex2d(right_x, down_y);
    glEnd();

    glutSwapBuffers();
}
