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

int fuerzas_iniciales(int n, double d_corte, double *pos_x, double *pos_y, double *pos_z, double *f_x_t, double *f_y_t, double *f_z_t, double *fuerzas, int nf) {
	
	double distancia_cuadrado, distancia;

	for (int ii=0; ii < n; ii++){
		f_x_t[ii] = 0.0;
		f_y_t[ii] = 0.0;
		f_z_t[ii] = 0.0;
	}

	//calculo las fuerzas sobre cada una de las partículas
	for (int i=0; i < (n-1); i++){
		for (int j=i+1; j < n; j++){
			
			distancia_cuadrado = pow(pos_x[i]-pos_x[j],2.)+pow(pos_y[i]-pos_y[j],2.)+pow(pos_z[i]-pos_z[j],2.);
			distancia = pow(distancia_cuadrado,1/2.);	
			
			if (distancia <= d_corte){

				int a = (int) (distancia*nf/4.0) - 1.0;	
						
				f_x_t[i] += fuerzas[a]*(pos_x[j]-pos_x[i])/distancia;
				f_x_t[j] -= fuerzas[a]*(pos_x[j]-pos_x[i])/distancia;
				f_y_t[i] += fuerzas[a]*(pos_y[j]-pos_y[i])/distancia;
				f_y_t[j] -= fuerzas[a]*(pos_y[j]-pos_y[i])/distancia;
				f_z_t[i] += fuerzas[a]*(pos_z[j]-pos_z[i])/distancia;
				f_z_t[j] -= fuerzas[a]*(pos_z[j]-pos_z[i])/distancia;
				
				//if (i==0) printf("distancia=%g\ta=%i\tf_x_t=%g\tfuerzas[a]=%g\tj=%i\n",distancia,a,f_x_t[i],fuerzas[a],j+1);
			}
		}
	}

return 0;
}
