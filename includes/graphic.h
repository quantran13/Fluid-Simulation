//
// Created by quan on 4/11/17.
//

#ifndef PARTICLE_SIM_GRAPHIC_H
#define PARTICLE_SIM_GRAPHIC_H

#include <stdio.h>
#include <GL/glut.h>
#include <fluid.h>
#include <utility.h>

typedef struct {
    double r;
    double b;
    double g;
} color_t;

void init_render();
void draw_cube(FluidCube *cube, perf_t *perf_struct);
void draw_square(int x, int y, double vsize, double hsize, color_t* color);

#endif //PARTICLE_SIM_GRAPHIC_H
