#include "ci.h"
#include "verlet.h"
#include "tablas.h"
#include "magnitudes.h"
#include "reescaleo.h"
#include "stdio.h"
#include <stdlib.h>
#include <math.h>

void escribir(double vector1[], double vector2[], double vector3[], int nt_maximo);

int main(int argc,char *argv[]){

	int densidad_entera;
	if (argc==2) 
     	{
       		sscanf(argv[1],"%i",&densidad_entera); //primer argumento: densidad
     	}

	int n = 125; //cantidad de particulas
	int nf = 4001; //cantidad de bins del potencial
	float densidad = 0.2*densidad_entera; //dada por el problema
	double d_corte = 0.5*pow((float)n/densidad,1/3.)-0.1; //distancia de corte (potencial)
	float T_actual = 2.0; //temperatura actual
	float T_deseada; //temperatura deseada
	int nt_maximo = 1001; //cantidad de temperaturas del barrido
	int n_term; //cantidad de pasos para termalizar
	int paso_descorrelacion = 500; //cantidad de pasos para descorrelacionar
	int muestreo_total = 1000; // cantidad de mediciones realizadas para cada temperatura
	double p;
	int pares_interaccion = 0;

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
	double   *energia_potencial_total = malloc(nt_maximo * sizeof(double));
	double   *energia_cinetica = malloc(nt_maximo * sizeof(double));
	double   *presion_de_interaccion = malloc(nt_maximo * sizeof(double));

	//Armo la tabla de potenciales y de fuerzas
	tabla_potenciales(nf,potenciales,d_corte); 
	tabla_fuerzas(nf,fuerzas,d_corte); 
	
	//Condiciones iniciales (usamos un lattice de tipo simple cubic)
	double lado = posiciones_iniciales(n,densidad,pos_x,pos_y,pos_z);
	velocidades_iniciales(T_actual,n,vel_x,vel_y,vel_z);
	
	//Calculo las fuerzas en t=0
	fuerzas_iniciales(n,f_x_t,f_y_t,f_z_t);

	for (int i = 0; i < n; i++) pares_interaccion += i;

	//hace un barrido en la temperatura
	for (int nt = 0; nt < nt_maximo; nt++) { 
		T_actual=2.0*(1.0-0.8*nt/(float)nt_maximo);
		
		energia_cinetica[nt] = 0.0;
		energia_potencial_total[nt] = 0.0;
		presion_de_interaccion[nt] = 0.0;

		//Termalización
		if(nt==0) {
			n_term = 500;
			T_deseada = T_actual;
		}
		else {
			n_term = 200;
			T_deseada = 2.0*(1.0-0.8*(nt+1)/(float)nt_maximo);	
		}
		
		for (int paso_term = 0; paso_term < n_term; paso_term++) {   		
		algoritmo_verlet(n,d_corte,pos_x,pos_y,pos_z,vel_x,vel_y,vel_z,f_x_t,f_y_t,f_z_t,fuerzas,nf,lado,vector_potencial,potenciales);
		}

		//Mido la energía cinetica y potencial cada "paso_descorrelacion" pasos
		for (int paso_muestreo = 0; paso_muestreo < muestreo_total; paso_muestreo++) {   		
			
			for (int paso_correlacion = 0; paso_correlacion < paso_descorrelacion; paso_correlacion++){
		p = algoritmo_verlet(n,d_corte,pos_x,pos_y,pos_z,vel_x,vel_y,vel_z,f_x_t,f_y_t,f_z_t,fuerzas,nf,lado,vector_potencial,potenciales);
			}			
		energia_cinetica[nt] += cinetica(n,vel_x,vel_y,vel_z,vector_cinetico)/(double)n;
		energia_potencial_total[nt] += energia_potencial(n,vector_potencial)/(double)n;
		presion_de_interaccion[nt] += p/(double)(pares_interaccion*3*lado);
		}
	//al final reescaleo las velocidaes para modificar la temperatura del sistema
	reescaleo_velocidades(T_actual,T_deseada,vel_x,vel_y,vel_z,n);
	if (nt==200 || nt==400 || nt==600 || nt==800) paso_descorrelacion = 2*paso_descorrelacion;

	}
    	//escribir(energia_cinetica, energia_potencial_total, presion_de_interaccion, nt_maximo); 
	
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
	free(presion_de_interaccion);

	return 0;
}

void escribir(double vector1[],double vector2[], double vector3[], int nt_maximo){
  int i;
  FILE *fp;
  fp = fopen("cinetica_potencial_presion_barrido_T.txt","w");
  for(i = 0; i < nt_maximo; i++){
    fprintf(fp, "%.6f \t %.6f \t %.6f \n",vector1[i],vector2[i],vector3[i]);
  }  

  fclose(fp);
}
