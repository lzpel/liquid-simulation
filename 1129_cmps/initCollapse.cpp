#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "defs.h"

#define OUTPUT_FILE "initCollapse.prof"

#define MIN_X  0.0
#define MIN_Y  0.0
#define MIN_Z  0.0
#define MAX_X  1.0
#define MAX_Y  0.2
#define MAX_Z  0.6

#define WAVE_HEIGHT 0.5
#define WAVE_WIDTH 0.25
#define PARTICLE_DISTANCE 0.02

int nx = (int)((MAX_X - MIN_X)/PARTICLE_DISTANCE)+6;
int ny = (int)((MAX_Y - MIN_Y)/PARTICLE_DISTANCE)+6;
int nz = (int)((MAX_Z - MIN_Z)/PARTICLE_DISTANCE)+6;
int nxy = nx*ny;
int nxyz = nx*ny*nz;
int NumberOfParticle;
Type* ParticleType;
float *Position;

int main(int argc, char** argv){

	printf("start mk_particle\n");
	printf("nx:%d ny:%d nz:%d nxyz:%d\n", nx, ny, nz,nxyz);

	ParticleType = (Type*)malloc(sizeof(Type)*nxyz);
	Position = (float*)malloc(sizeof(float)*nxyz*3);
	int NumberOfParticle = 0;

	for(int iz=0;iz<nz;iz++){
	for(int iy=0;iy<ny;iy++){
	for(int ix=0;ix<nx;ix++){
		int ip = iz*nxy + iy*nx + ix;
		ParticleType[ip] = Type::GHOST;
		Position[ip*3  ] = MIN_X + PARTICLE_DISTANCE*(ix-3);
		Position[ip*3+1] = MIN_Y + PARTICLE_DISTANCE*(iy-3);
		Position[ip*3+2] = MIN_Z + PARTICLE_DISTANCE*(iz-3);
	}}}

	for(int iz=0;iz<nz;iz++){
	for(int iy=0;iy<ny;iy++){
	for(int ix=0;ix<nx;ix++){
		int ip = iz*nxy + iy*nx + ix;
		if(ix<3 || ix >=nx-3 || iy<3 || iy>=ny-3 || iz<3){
			ParticleType[ip] = Type::WALL;
			NumberOfParticle++;
		}
		else if(Position[ip*3+2] <= WAVE_HEIGHT && Position[ip*3] <= WAVE_WIDTH){
			ParticleType[ip] = Type::FLUID;
			NumberOfParticle++;
		}
	}}}


	FILE* fp;
	fp = fopen(OUTPUT_FILE, "w");
    fprintf(fp, "%lf %d\n",PARTICLE_DISTANCE, 1000);
	int k=0;
	for(int iz=0;iz<nz;iz++){
	for(int iy=0;iy<ny;iy++){
	for(int ix=0;ix<nx;ix++){
		int ip = iz*nxy + iy*nx + ix;
		if(ParticleType[ip]==Type::GHOST)continue;
		fprintf(fp," %d %lf %lf %lf\n",ParticleType[ip],Position[ip*3],Position[ip*3+1],Position[ip*3+2]);
		k++;
	}}}
	fclose(fp);

	free(ParticleType);	free(Position);
	printf("end mk_particle\n");
	return 0;
}
