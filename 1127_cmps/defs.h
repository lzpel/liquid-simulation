//
// Created by misumi on 2021/11/27.
//

#pragma once

enum class Type { GHOST=0, FLUID, WALL, DEM };
#define VectorScalar(x,y,z,scalar) {(x)*(scalar),(y)*(scalar),(z)*(scalar)}
#define InnerProduct(a,b) ((a)[0]*(b)[0]+(a)[1]*(b)[1]+(a)[2]*(b)[2])
#define Kernel(distance, radius) ((radius/distance) - 1.0)

#define clockStart() clock_t start_clock = clock();
#define clockPrint() printf("clock %s:%f\n",__func__,(double)(clock() - start_clock) / CLOCKS_PER_SEC);