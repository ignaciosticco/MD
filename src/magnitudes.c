#include "magnitudes.h"
#include "stdio.h"
#include <stdlib.h>
#include <math.h>

int potencial(int n, double d_corte, double *pos_x, double *pos_y, double *pos_z, double *vector_potencial, double *potenciales, int nf){
	double distancia_cuadrado, distancia;

	for (int i = 0; i < n; i++){
		vector_potencial[i]=0.0;
	}

	//calculo la energia potencial sobre cada partícula
	for (int i=0; i < (n-1); i++){
		for (int j=i+1; j < n; j++){
			
			distancia_cuadrado = pow(pos_x[i]-pos_x[j],2.)+pow(pos_y[i]-pos_y[j],2.)+pow(pos_z[i]-pos_z[j],2.);
			distancia = pow(distancia_cuadrado,1/2.);	
			
			if (distancia <= d_corte){

				int a = (int) (distancia*nf/4.0) - 1.0;	
						
				vector_potencial[i] += potenciales[a];
				vector_potencial[j] += potenciales[a];
			}
		}
	}
    return 0;
}

int cinetica(int n, double *vel_x, double *vel_y, double *vel_z, double *vector_cinetico){
	
	double k;
	for (int i = 0; i < n; i++){
		k = vel_x[i]*vel_x[i]+vel_y[i]*vel_y[i]+vel_z[i]*vel_z[i];
		vector_cinetico[i]=pow(k,1/2.);
	}
    return 0;
}
