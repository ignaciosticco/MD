#include "reescaleo.h"
#include "stdio.h"
#include <stdlib.h>
#include <math.h>

int reescaleo_velocidades(float T_actual, float T_deseada, double *vel_x, double *vel_y, double *vel_z, int n){

	double raiz;
	raiz=sqrt(T_deseada/T_actual);
	
	for (int i = 0; i < n; i++){
		vel_x[i] =raiz * vel_x[i]; 
		vel_y[i] =raiz * vel_y[i]; 
		vel_z[i] =raiz * vel_z[i]; 
	}

return 0;
}
