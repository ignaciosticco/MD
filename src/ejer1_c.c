#include "ci.h"
#include "verlet.h"
#include "tablas.h"
#include "magnitudes.h"
#include "stdio.h"
#include <stdlib.h>
#include <math.h>

void 	escribir(double vector1[],double vector2[],double vector3[],int niter);
double 	mean(double vector[],double size);
double 	std(double vector[],double size,double mean);
int		rescaleo(double *vel_x,double *vel_y,double *vel_z,double T_actual,double T_deseada,int n);


int main(){
	int 	n = 512; //cantidad de particulas
	int 	nf = 4001; //cantidad de bins del potencial
	float 	densidad = 0.8442; //dada por el problema
	double 	d_corte = 0.5*pow((float)n/densidad,1/3.)-0.1; //distancia de corte (potencial)
	int 	t_termaliz = 1; //500;
	int 	t_correlacion = 2;// 250;
	int 	sample_max = 1;
	int   	cant_T = 100;
  	double 	T_min = 0.728;
  	double	T_max = 4.0;
  	double 	delta_T = (T_max - T_min)/(cant_T-1);
  	double 	vector_T[cant_T];
  	double 	energia_total[cant_T];
  	double 	std_energia[cant_T];
  	double  mean_energy;
  	double  T, T_actual, T_deseada;

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

	// Inicializo los vectores
	for(int itera_T = 0;itera_T<cant_T;itera_T++){
	    vector_T[itera_T] = T_max - delta_T*itera_T;
	    energia_total[itera_T] = 0.0;
	    std_energia[itera_T] = 0.0;
  	}
  	//Armo la tabla de potenciales y de fuerzas
	tabla_potenciales(nf,potenciales,d_corte); 
	tabla_fuerzas(nf,fuerzas,d_corte); 

	

	//Condiciones iniciales (usamos un lattice de tipo simple cubic)
	double lado = posiciones_iniciales(n,densidad,pos_x,pos_y,pos_z);
	velocidades_iniciales(T,n,vel_x,vel_y,vel_z);
	//Calculo las fuerzas en t=0
	fuerzas_iniciales(n,f_x_t,f_y_t,f_z_t);

	// Loop de temperatura
	for (int itera_T = 0;itera_T<cant_T-1;itera_T++){
		T_actual = vector_T[itera_T];
		T_deseada = vector_T[itera_T+1];

		// Loop de termalizacion
		for (int t = 0; t < t_termaliz; t++){
			algoritmo_verlet(n,d_corte,pos_x,pos_y,pos_z,vel_x,vel_y,vel_z,f_x_t,f_y_t,f_z_t,fuerzas,nf,lado,vector_potencial,potenciales);
		}
		// Loop de sampleo 
		for (int sample = 0; sample < sample_max; sample++){
			for (int t = 0; t < t_correlacion; t++) {   		
				algoritmo_verlet(n,d_corte,pos_x,pos_y,pos_z,vel_x,vel_y,vel_z,f_x_t,f_y_t,f_z_t,fuerzas,nf,lado,vector_potencial,potenciales);
			}	
			energia_total[sample] = cinetica(n,vel_x,vel_y,vel_z,vector_cinetico) + energia_potencial(n,vector_potencial);
		}
		mean_energy = mean(energia_total,sample_max);
		energia_total[itera_T] = mean_energy ;
	    std_energia[itera_T] = std(energia_total,sample_max,mean_energy);

	    rescaleo(vel_x,vel_y,vel_z,T_actual,T_deseada,n);
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

    escribir(vector_T, energia_total, std_energia , cant_T); 
	

	return 0;
}

void escribir(double vector1[],double vector2[],double vector3[],int niter){
  int i;
  FILE *fp;
  fp = fopen("energias_1c.txt","w");
  for(i=0;i<niter;i++){
    fprintf(fp, "%.6f \t %.6f \t %.6f \n",vector1[i],vector2[i],vector3[i]);
  }  

  fclose(fp);
}

double mean(double vector[],double size){
	double mean=0.0;
	for (int i = 0; i < size; i++){
		mean+=vector[i];
	}
	return mean/size;
}

double std(double vector[],double size,double mean){
	double std=0.0;
	for (int i = 0; i < size; i++){
		std+= (mean-vector[i])*(mean-vector[i]);
	}
	return sqrt(std/size);
}

int rescaleo(double *vel_x,double *vel_y,double *vel_z,double T_actual,double T_deseada,int n){
	double raiz = sqrt(T_deseada/T_actual);
	for (int i = 0; i < n; i++){
		vel_x[i] = raiz*vel_x[i];
		vel_y[i] = raiz*vel_y[i]; 
		vel_z[i] = raiz*vel_z[i];
	}
	return 0;

}