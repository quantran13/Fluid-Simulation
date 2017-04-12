//
// Created by quan on 4/11/17.
//

#ifndef PARTICLE_SIM_GRAPHIC_H
#define PARTICLE_SIM_GRAPHIC_H

#include <GL/glut.h>
#include "fluid.h"

void init_render();
void draw_cube(FluidCube *cube, int n);
void draw_square(int x, int y, double vsize, double hsize, double color);

#endif //PARTICLE_SIM_GRAPHIC_H
