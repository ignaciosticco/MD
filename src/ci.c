#include "ci.h"
//#include "time.h"
#include "stdio.h"
#include <stdlib.h>
#include <math.h>

int posiciones_iniciales(int n, float densidad,float *vector_x,float *vector_y,float *vector_z){
	float cociente = (float)n/densidad;
	float lado = pow(cociente,1/3.);
	int partxlado = pow(n,1/3.);
	float a = lado/(float)partxlado;

	for(int k=0; k<partxlado;k++){
		for(int j=0; j<partxlado;j++){
			for(int i=0; i<partxlado;i++){
				vector_x[i+partxlado*j+partxlado*partxlado*k] = i*a+a/2;
				vector_y[i+partxlado*j+partxlado*partxlado*k] = j*a+a/2;
				vector_z[i+partxlado*j+partxlado*partxlado*k] = k*a+a/2;
			}
		}
	}
	for(int i=0;i<n;i++){
		//printf("%f\t%f\t%f\n",vector_x[i],vector_y[i],vector_z[i]);
	}
	return 0;
}

int velocidades_iniciales(int n,float *vector_x,float *vector_y,float *vector_z){

		float vx_medio=0.0;
		float vy_medio=0.0;
		float vz_medio=0.0;

	for(int i=0; i<n;i++){
		float ex = (float)rand()/RAND_MAX;
		float ey = (float)rand()/RAND_MAX;
		float ez = (float)rand()/RAND_MAX;

		vector_x[i]=ex;
		vector_y[i]=ey;
		vector_z[i]=ez;

		vx_medio+=ex;
		vy_medio+=ey;
		vz_medio+=ez;
	}
	vx_medio=vx_medio/n;
	vy_medio=vy_medio/n;
	vz_medio=vz_medio/n;

	for(int i=0; i<n;i++){
	
		vector_x[i]=vector_x[i]-vx_medio;
		vector_y[i]=vector_y[i]-vy_medio;
		vector_z[i]=vector_z[i]-vz_medio;
	}	

	for(int i=0;i<n;i++){
		printf("%f\t%f\t%f\n",vector_x[i],vector_y[i],vector_z[i]);
	}
	return 0;
}


