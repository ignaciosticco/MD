#include "verlet.h"
#include "stdio.h"
#include <stdlib.h>
#include <math.h>

int algoritmo_verlet(int n, double d_corte, double *pos_x, double *pos_y, double *pos_z, double *vel_x, double *vel_y, double *vel_z, double *f_x_t, double *f_y_t, double *f_z_t, double *fuerzas, int nf) {

	double h = 0.001;
	double distancia_cuadrado, distancia;

	double   *f_x_t_h = malloc(n * sizeof(double));
	double   *f_y_t_h = malloc(n * sizeof(double));
	double   *f_z_t_h = malloc(n * sizeof(double));	

	//actualizo las posiciones
	for (int i=0; i < n; i++){
		pos_x[i] = pos_x[i] + vel_x[i]*h+(0.5)*f_x_t[i]*h*h;
		pos_y[i] = pos_y[i] + vel_y[i]*h+(0.5)*f_y_t[i]*h*h;
		pos_z[i] = pos_z[i] + vel_z[i]*h+(0.5)*f_z_t[i]*h*h;
	}
	
	for (int i=0; i < n; i++){
		f_x_t_h[i] = 0.0;
		f_y_t_h[i] = 0.0;
		f_z_t_h[i] = 0.0;
	}

	//calculo las nuevas fuerzas sobre cada una de las partículas (en el paso t+h)
	for (int i = 0; i < (n-1); i++){
		for (int j=i+1; j < n; j++){
			
			distancia_cuadrado = pow(pos_x[i]-pos_x[j],2.)+pow(pos_y[i]-pos_y[j],2.)+pow(pos_z[i]-pos_z[j],2.);
			distancia = pow(distancia_cuadrado,1/2.);
			
			if (distancia <= d_corte){
	
				int a = (int) (distancia*nf/4.0) - 1.0;		
		
				f_x_t_h[i] += fuerzas[a]*(pos_x[j]-pos_x[i])/distancia;
				f_x_t_h[j] -= fuerzas[a]*(pos_x[j]-pos_x[i])/distancia;
				f_y_t_h[i] += fuerzas[a]*(pos_y[j]-pos_y[i])/distancia;
				f_y_t_h[j] -= fuerzas[a]*(pos_y[j]-pos_y[i])/distancia;
				f_z_t_h[i] += fuerzas[a]*(pos_z[j]-pos_z[i])/distancia;
				f_z_t_h[j] -= fuerzas[a]*(pos_z[j]-pos_z[i])/distancia;
			}
		}
	}
	
	//actualizo las velocidades
	for (int i = 0; i < n; i++){
		vel_x[i] = vel_x[i] + (0.5)*(f_x_t[i]+f_x_t_h[i])*h;
		vel_y[i] = vel_y[i] + (0.5)*(f_y_t[i]+f_y_t_h[i])*h;
		vel_z[i] = vel_z[i] + (0.5)*(f_z_t[i]+f_z_t_h[i])*h;
	}

	//las nuevas fuerzas me sirven para calcular las nuevas posiciones (cuando entre de nuevo al "verlet")
	for (int i = 0; i < n; i++){
		f_x_t[i] = f_x_t_h[i];
		f_y_t[i] = f_y_t_h[i];
		f_z_t[i] = f_z_t_h[i];
	}

	free(f_x_t_h);
	free(f_y_t_h);
	free(f_z_t_h);

	return 0;
}
