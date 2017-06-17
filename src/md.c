#include "ci.h"
#include "stdio.h"
#include <stdlib.h>
#include <math.h>


int main(){
	int n = 27;
	float densidad = 0.333;
	float   *posicion_x = malloc(n * sizeof(float));
	float   *posicion_y = malloc(n * sizeof(float));
	float   *posicion_z = malloc(n * sizeof(float));

	float   *velocidad_x = malloc(n * sizeof(float));
	float   *velocidad_y = malloc(n * sizeof(float));
	float   *velocidad_z = malloc(n * sizeof(float));

	posiciones_iniciales(n,densidad,posicion_x,posicion_y,posicion_z);
	velocidades_iniciales(n,velocidad_x,velocidad_y,velocidad_z);

	free(posicion_x);
	free(posicion_y);
	free(posicion_z);
	free(velocidad_x);
	free(velocidad_y);
	free(velocidad_z);

	return 0;
}
