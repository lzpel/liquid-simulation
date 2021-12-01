//
// Created by misumi on 2021/11/27.
//

#pragma once

#include <dirent.h>
#include <string>

enum class Type { GHOST=0, FLUID, WALL, DEM };
inline double distance2(double x,double y,double z){return x*x+y*y+z*z;}


#define VectorScalar(x,y,z,scalar) {(x)*(scalar),(y)*(scalar),(z)*(scalar)}
#define InnerProduct(a,b) ((a)[0]*(b)[0]+(a)[1]*(b)[1]+(a)[2]*(b)[2])
#define Kernel(distance, radius) ((radius/distance) - 1.0)

#define clockStart() clock_t start_clock = clock();
#define clockPrint() printf("clock %s:%f\n",__func__,(double)(clock() - start_clock) / CLOCKS_PER_SEC);

char removeProf(void){
    struct dirent *dir;
    DIR *di = opendir("."); //specify the directory name
    if (!di)return -1;
    while ((dir = readdir(di)) != NULL){
        char tmp[256];
        strcpy(tmp,dir->d_name);
        const char* ptr1=strtok(tmp,".");
        const char* ptr2=strtok(NULL,".");
        if(ptr2!=NULL&&strcmp(ptr2,"prof")==0)remove(dir->d_name);
    }
    closedir(di);
    return 0;
}