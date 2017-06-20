#include "ci.h"
#include "potencial.h"
#include "stdio.h"
#include <stdlib.h>
#include <math.h>


int main(){
	int n = 27;
	float densidad = 0.8442;
	float   *posicion_x = malloc(n * sizeof(float));
	float   *posicion_y = malloc(n * sizeof(float));
	float   *posicion_z = malloc(n * sizeof(float));
	float   *velocidad_x = malloc(n * sizeof(float));
	float   *velocidad_y = malloc(n * sizeof(float));
	float   *velocidad_z = malloc(n * sizeof(float));
	float   *vector_potencial = malloc(n * sizeof(float));


	// Inicializo vector de potencial
	for(int i=0;i<n;i++){
		vector_potencial[i] = 0.0;
	}

	posiciones_iniciales(n,densidad,posicion_x,posicion_y,posicion_z);
	//velocidades_iniciales(n,velocidad_x,velocidad_y,velocidad_z);
	potencial(posicion_x,posicion_y,posicion_z,vector_potencial,n);

	printf("\n Tabla de potenciales:\n \n");
	for(int i=0;i<n;i++){
		printf("%f\t%f\t%f\t%f\n",posicion_x[i],posicion_y[i],posicion_z[i],vector_potencial[i]);
	}

	free(posicion_x);
	free(posicion_y);
	free(posicion_z);
	free(velocidad_x);
	free(velocidad_y);
	free(velocidad_z);
	free(vector_potencial);

	return 0;
}
