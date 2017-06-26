#include "ci.h"
//#include "time.h"
#include "stdio.h"
#include <stdlib.h>
#include <math.h>

double posiciones_iniciales(int n, float densidad,double *vector_x,double *vector_y,double *vector_z){
	float cociente = (float)n/densidad;
	double lado = pow(cociente,1/3.);
	int partxlado = pow(n,1/3.);
	double a = lado/(double)partxlado;

	for(int k=0; k<partxlado;k++){
		for(int j=0; j<partxlado;j++){
			for(int i=0; i<partxlado;i++){
				vector_x[i+partxlado*j+partxlado*partxlado*k] = i*a+a/2;
				vector_y[i+partxlado*j+partxlado*partxlado*k] = j*a+a/2;
				vector_z[i+partxlado*j+partxlado*partxlado*k] = k*a+a/2;
			}
		}
	}

return lado;
}

int velocidades_iniciales(int n,double *vector_x,double *vector_y,double *vector_z){

		double vx_medio=0.0;
		double vy_medio=0.0;
		double vz_medio=0.0;

	for(int i=0; i<n;i++){
		double ex = (double)rand()/RAND_MAX-0.5;
		double ey = (double)rand()/RAND_MAX-0.5;
		double ez = (double)rand()/RAND_MAX-0.5;

		vector_x[i]=ex;
		vector_y[i]=ey;
		vector_z[i]=ez;

		vx_medio+=ex;
		vy_medio+=ey;
		vz_medio+=ez;
	}
	vx_medio=vx_medio/(double)n;
	vy_medio=vy_medio/(double)n;
	vz_medio=vz_medio/(double)n;

	for(int i=0; i<n;i++){
	
		vector_x[i]=vector_x[i]-vx_medio;
		vector_y[i]=vector_y[i]-vy_medio;
		vector_z[i]=vector_z[i]-vz_medio;
	}	

return 0;
}

int fuerzas_iniciales(int n, double *f_x_t, double *f_y_t, double *f_z_t) {
	
	for (int ii=0; ii < n; ii++){
		f_x_t[ii] = 0.0;
		f_y_t[ii] = 0.0;
		f_z_t[ii] = 0.0;
	}

return 0;
}
