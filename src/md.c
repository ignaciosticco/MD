#include "ci.h"
#include "verlet.h"
#include "tablas.h"
#include "magnitudes.h"
#include "stdio.h"
#include <stdlib.h>
#include <math.h>

void escribir(double *f, double *m, double *e, int n);

int main(){
	int n = 27; //cantidad de particulas
	int nf = 4001; //cantidad de bins del potencial
	float densidad = 0.8442; //dada por el problema
	double d_corte = 0.5*pow((float)n/densidad,1/3.); //distancia de corte (potencial)
	int np = 5000;

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

	//Armo la tabla de potenciales y de fuerzas
	tabla_potenciales(nf,potenciales); 
	tabla_fuerzas(nf,fuerzas); 
	
	//Condiciones iniciales (usamos un lattice de tipo simple cubic)
	double lado = posiciones_iniciales(n,densidad,pos_x,pos_y,pos_z);
	velocidades_iniciales(n,vel_x,vel_y,vel_z);

	//Calculo las fuerzas en t=0
	fuerzas_iniciales(n,d_corte,pos_x,pos_y,pos_z,f_x_t,f_y_t,f_z_t,fuerzas,nf);

	//calculo la energia potencial y cinetica en t=0
	potencial(n,d_corte,pos_x,pos_y,pos_z,vector_potencial,potenciales,nf);
	cinetica(n,vel_x,vel_y,vel_z,vector_cinetico);

	for (int p = 0; p < np; p++) {   		
		algoritmo_verlet(n,d_corte,pos_x,pos_y,pos_z,vel_x,vel_y,vel_z,f_x_t,f_y_t,f_z_t,fuerzas,nf,lado);
		potencial(n,d_corte,pos_x,pos_y,pos_z,vector_potencial,potenciales,nf);
		cinetica(n,vel_x,vel_y,vel_z,vector_cinetico);
	}

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

	return 0;
}

void escribir(double *f, double *m, double *e, int n)
{
	int i;
	FILE *fp;
	fp=fopen("r0.txt","a");
	for (i=0;i<n;i++) fprintf(fp,"%.6f\n",f[i]);
	fclose(fp);
}
