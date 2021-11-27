#pragma once

#include "math.h"
#include "vector"
#include "cfloat"

#define BUCKETINDEXXYZ(P) const int ix = int((P[0] - min[0]) / length) + 1, iy = int((P[1] - min[1]) / length) + 1,iz = int((P[2] - min[2]) / length) + 1
#define BUCKETINDEX(dx,dy,dz) ((iz+(dz)) * sizeX * sizeY + (iy+(dy)) * sizeX + (ix+(dx)))

struct Bucket {
    double max[3], min[3];
    signed sizeX, sizeY, sizeZ;
    double length;
    std::vector<int> first;
    std::vector<int> nxt;
    signed list333index;
    signed list333[27];

    Bucket() {
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
        const int ib=BUCKETINDEX(0,0,0);
        nxt[index] = first[ib];
        first[ib] = index;
    }

    int iterator(const double *pos) {
        BUCKETINDEXXYZ(pos);
        list333index=0;
        list333[0]=BUCKETINDEX(0,0,0);
        list333[1]=BUCKETINDEX(0,0,+1);
        list333[2]=BUCKETINDEX(0,0,-1);
        list333[3]=BUCKETINDEX(0,+1,0);
        list333[4]=BUCKETINDEX(0,+1,+1);
        list333[5]=BUCKETINDEX(0,+1,-1);
        list333[6]=BUCKETINDEX(0,-1,0);
        list333[7]=BUCKETINDEX(0,-1,+1);
        list333[8]=BUCKETINDEX(0,-1,-1);
        list333[9]=BUCKETINDEX(+1,0,0);
        list333[10]=BUCKETINDEX(+1,0,+1);
        list333[11]=BUCKETINDEX(+1,0,-1);
        list333[12]=BUCKETINDEX(+1,+1,0);
        list333[13]=BUCKETINDEX(+1,+1,+1);
        list333[14]=BUCKETINDEX(+1,+1,-1);
        list333[15]=BUCKETINDEX(+1,-1,0);
        list333[16]=BUCKETINDEX(+1,-1,+1);
        list333[17]=BUCKETINDEX(+1,-1,-1);
        list333[18]=BUCKETINDEX(-1,0,0);
        list333[19]=BUCKETINDEX(-1,0,+1);
        list333[20]=BUCKETINDEX(-1,0,-1);
        list333[21]=BUCKETINDEX(-1,+1,0);
        list333[22]=BUCKETINDEX(-1,+1,+1);
        list333[23]=BUCKETINDEX(-1,+1,-1);
        list333[24]=BUCKETINDEX(-1,-1,0);
        list333[25]=BUCKETINDEX(-1,-1,+1);
        list333[26]=BUCKETINDEX(-1,-1,-1);
        return first[list333[list333index]];
    }
    inline int next(int n){
        int i=nxt[n];
        if(i!=-1)return i;
        list333index++;
        if(list333index<27)return first[list333[list333index]];
        return -1;
    }
};