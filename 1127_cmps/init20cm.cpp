#include "stdio.h"
#include "line.h"

int main(void){
    const int height=200;//mm
    const int width=200;//mm
    const int depth=20;//mm
    const int particleDistance=4;//mm
    const int thickNum=3;//個
    const int density=1000;//kg/m^3
    FILE *fp = fopen("init20cm.prof", "w");
    fprintf(fp, "%lf %lf\n",particleDistance,density);
    for(int ix=-thickNum;ix<width/particleDistance+thickNum;ix++){
        for(int iy=-thickNum;iy<depth/particleDistance+thickNum;iy++){
            for(int iz=-thickNum;iz<height/particleDistance;iz++){
                Type type=Type::WALL;
                if((ix>=0)&&(iy>=0)&&(iz>=0)&&(ix<width/particleDistance)&&(iy<depth/particleDistance)&&(iz<height/particleDistance)){
                    type=Type::FLUID;
                }else{
                    type=Type::WALL;
                }
                fprintf(fp, "%d %lf %lf %lf\n",type,ix*particleDistance/1000.0,iy*particleDistance/1000.0,iz*particleDistance/1000.0);
            }
        }
    }
    fclose(fp);
    return 0;
}