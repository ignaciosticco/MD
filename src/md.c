#include "ci.h"
#include "verlet.h"
#include "tablas.h"
#include "magnitudes.h"
#include "stdio.h"
#include <stdlib.h>
#include <math.h>

void escribir(double vector1[],double vector2[],double vector3[],int niter);

int main(){
	int n = 512; //cantidad de particulas
	int nf = 4001; //cantidad de bins del potencial
	float densidad = 0.8442; //dada por el problema
	double d_corte = 0.5*pow((float)n/densidad,1/3.)-0.1; //distancia de corte (potencial)
	int np = 1;
	float T = 0.728;

	double   *pos_x = malloc(n * sizeof(double));
	double   *pos_y = malloc(n * sizeof(double));
	double   *pos_z = malloc(n * sizeof(double));
	double   *vel_x = malloc(n * sizeof(double));
	double   *vel_y = malloc(n * sizeof(double));
	double   *vel_z = malloc(n * sizeof(double));
	double   *f_x_t = malloc(n * sizeof(double));
	double   *f_y_t = malloc(n * sizeof(double));
	double   *f_z_t = malloc(n * sizeof(double));

	double   *potenciales = malloc(nf * sizeof(double));
	double   *fuerzas = malloc(nf * sizeof(double));

	double   *vector_potencial = malloc(n * sizeof(double));
	double   *vector_cinetico = malloc(n * sizeof(double));
	double   *energia_potencial_total = malloc(np * sizeof(double));
	double   *energia_cinetica = malloc(np * sizeof(double));

	//Armo la tabla de potenciales y de fuerzas
	tabla_potenciales(nf,potenciales); 
	tabla_fuerzas(nf,fuerzas); 
	
	//Condiciones iniciales (usamos un lattice de tipo simple cubic)
	double lado = posiciones_iniciales(n,densidad,pos_x,pos_y,pos_z);
	velocidades_iniciales(T,n,vel_x,vel_y,vel_z);
	escribir(pos_x,pos_y,pos_z,n);

	//Calculo las fuerzas en t=0
	fuerzas_iniciales(n,f_x_t,f_y_t,f_z_t);

	printf("Cinetica 	potencial 	total\n");
	for (int p = 0; p < np; p++) {   		
		algoritmo_verlet(n,d_corte,pos_x,pos_y,pos_z,vel_x,vel_y,vel_z,f_x_t,f_y_t,f_z_t,fuerzas,nf,lado,vector_potencial,potenciales);
		energia_cinetica[p] = cinetica(n,vel_x,vel_y,vel_z,vector_cinetico);
		energia_potencial_total[p] = energia_potencial(n,vector_potencial);

		printf("%g\t%g\t%g\n", energia_cinetica[p],energia_potencial_total[p], energia_cinetica[p]+energia_potencial_total[p]);
	}
    	//escribir(vector_potencial, vector_potencial, vector_cinetico, n); //El primer parametro no sé para que sirve
	
	free(pos_x);
	free(pos_y);
	free(pos_z);
	free(vel_x);
	free(vel_y);
	free(vel_z);
	free(f_x_t);
	free(f_y_t);
	free(f_z_t);
	free(fuerzas);
	free(potenciales);
	free(vector_potencial);
	free(vector_cinetico);
	free(energia_potencial_total);
	free(energia_cinetica);

	return 0;
}

void escribir(double vector1[],double vector2[],double vector3[],int niter){
  int i;
  FILE *fp;
  fp = fopen("output.txt","w");
  for(i=0;i<niter;i++){
    fprintf(fp, "%.6f \t %.6f \t %.6f \n",vector1[i],vector2[i],vector3[i]);
  }  

  fclose(fp);
}
