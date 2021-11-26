//
// Created by misumi on 2021/11/27.
//

#pragma once

enum class Type { GHOST=0, FLUID, WALL, DEM };
struct Vector{
    double x,y,z;
    Vector(double x,double y,double z){
        x=x;
        y=y;
        z=z;
    }
    inline Vector scalar(double scalar){
        return Vector(x*scalar,y*scalar,z*scalar);
    }
};
