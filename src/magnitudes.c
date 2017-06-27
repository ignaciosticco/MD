#include "magnitudes.h"
#include "stdio.h"
#include <stdlib.h>
#include <math.h>

int potencial(double *vector_potencial, double va, int a,int i, int j){
	
	vector_potencial[i] += va;
	vector_potencial[j] += va;
	
    return 0;
}

double cinetica(int n, double *vel_x, double *vel_y, double *vel_z, double *vector_cinetico){
	
	double k;
	double energia_cinetica = 0.0;	
	for (int i = 0; i < n; i++){
		k = vel_x[i]*vel_x[i]+vel_y[i]*vel_y[i]+vel_z[i]*vel_z[i];
		//printf("k: %f\n", k );
		vector_cinetico[i]=k; // pow(k,1/2.);
		energia_cinetica+=vector_cinetico[i];
	}
    return energia_cinetica;
}

double energia_potencial(int n,double *vector_potencial){
	double energia_potencial = 0.0;
	for (int i = 0; i < n; ++i){
		energia_potencial+=vector_potencial[i];
	}
	return energia_potencial;
}
