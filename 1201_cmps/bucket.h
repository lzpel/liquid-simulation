#pragma once

#include "math.h"
#include "vector"
#include "cfloat"

#define BUCKETINDEXXYZ(P) const int ix = int((P[0] - min[0]) / length) + 1, iy = int((P[1] - min[1]) / length) + 1,iz = int((P[2] - min[2]) / length) + 1
#define BUCKETINDEX(dx,dy,dz) ((dz) * sizeX * sizeY + (dy) * sizeX + (dx))

struct Bucket {
    double max[3], min[3];
    signed sizeX, sizeY, sizeZ;
    double length;
    std::vector<int> first;
    std::vector<int> nxt;
    std::vector<int> near;

    Bucket() {
        //(0.046*3)^3/(0.02^3)=328
        near=std::vector<int>(1024);
        rangeInit();
    }

    void rangeInit(){
        min[0] = min[1] = min[2] = DBL_MAX;
        max[0] = max[1] = max[2] = DBL_MIN;
    }
    void rangeUpdate(const double *pos) {
        min[0] = fmin(min[0], pos[0]);
        min[1] = fmin(min[1], pos[1]);
        min[2] = fmin(min[2], pos[2]);
        max[0] = fmax(max[0], pos[0]);
        max[1] = fmax(max[1], pos[1]);
        max[2] = fmax(max[2], pos[2]);
    }

    void rangeUpdateFinish(const int size, const double len) {
        length = len;
        sizeX = int((max[0] - min[0]) / length) + 1 + 2;
        sizeY = int((max[1] - min[1]) / length) + 1 + 2;
        sizeZ = int((max[2] - min[2]) / length) + 1 + 2;
        nxt.resize(size);
        first.resize(sizeX * sizeY * sizeZ);
        printf("%lf <= x <= %lf, %lf <= y <= %lf, %lf <= z <= %lf\n", min[0], max[0], min[1], max[1], min[2], max[2]);
        printf("bucketLength = %lf, bucketSizeX = %d, bucketSizeY = %d, bucketSizeZ = %d, bucketSizeXYZ = %d\n", length, sizeX, sizeY, sizeZ, sizeX * sizeY * sizeZ);
        fflush(stdout);
    }

    void clear() {
        for (int i = 0; i < nxt.size(); i++) { nxt[i] = -1; }
        for (int i = 0; i < first.size(); i++) { first[i] = -1; }
    }

    void add(const double pos[3], int index) {
        BUCKETINDEXXYZ(pos);
        const int ib=BUCKETINDEX(ix,iy,iz);
        nxt[index] = first[ib];
        first[ib] = index;
    }

    inline const int* iterator(const double *pos) {
        BUCKETINDEXXYZ(pos);
        near.clear();
        for(int z=0;z<3;z++){
            for(int y=0;y<3;y++){
                for(int x=0;x<3;x++){
                    for(int r=first[BUCKETINDEX(ix-1+x,iy-1+y,iz-1+z)];r!=-1;r=nxt[r]){
                        near.emplace_back(r);
                    }
                }
            }
        }
        near.emplace_back(-1);
        return &near[0];
    }
};