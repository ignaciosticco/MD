#include "ci.h"
#include "stdio.h"
#include <stdlib.h>
#include <math.h>

double posiciones_iniciales(int n, float densidad,double *vector_x,double *vector_y,double *vector_z){
	int partxlado;
	float cociente = (float)n/densidad;
	double lado = pow(cociente,1/3.);
	partxlado = pow(n,1/3.);
	if(n % partxlado>0) partxlado = partxlado + 1;
	double a = lado/(double)partxlado;
	//printf("lado: %g\n",lado);
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

int velocidades_iniciales(float T_actual,int n,double *vector_x,double *vector_y,double *vector_z){

	double vx_medio=0.0;
	double vy_medio=0.0;
	double vz_medio=0.0;
	float sigma = sqrt(T_actual);
	float media = 0.0;

	for(int i=0; i<n;i++){
		double ex = rand_gausiano(media,sigma);
		double ey = rand_gausiano(media,sigma);
		double ez = rand_gausiano(media,sigma);

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
		//printf("velocidad: %f \n:",vector_x[i]*vector_x[i]+vector_y[i]*vector_y[i]+vector_z[i]*vector_z[i]);
	
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


//////////////////// Codigo que genera un random gausiano  //////////////////// 
double drand(){   
/* uniform distribution, (0..1] */
  return (rand()+1.0)/(RAND_MAX+1.0);
}
double random_normal(){
/* normal distribution, centered on 0, std dev 1 */
  return sqrt(-2*log(drand())) * cos(2*M_PI*drand());
}
double rand_gausiano(float media,float sigma){
  return media + sigma*random_normal();
}
// Fuente: https://stackoverflow.com/questions/7034930/how-to-generate-gaussian-pseudo-random-numbers-in-c-for-a-given-mean-and-varianc
/////////// //////////////////////////  ///////////
