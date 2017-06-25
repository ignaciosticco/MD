#include "verlet.h"
#include "stdio.h"
#include <stdlib.h>
#include <math.h>

int algoritmo_verlet(int n, double d_corte, double *pos_x, double *pos_y, double *pos_z, double *vel_x, double *vel_y, double *vel_z, double *f_x_t, double *f_y_t, double *f_z_t, double *fuerzas, int nf,double lado) {

	double h = 0.001;
	double distancia_cuadrado, distancia, actualizacion, d_corte_cuadrado, dx, dy, dz;

	double   *f_x_t_h = malloc(n * sizeof(double));
	double   *f_y_t_h = malloc(n * sizeof(double));
	double   *f_z_t_h = malloc(n * sizeof(double));	

	//actualizo las posiciones
	for (int i=0; i < n; i++){
		pos_x[i] = pos_x[i] + vel_x[i]*h+(0.5)*f_x_t[i]*h*h;
		pos_x[i] = pos_x[i] - lado* (int)(pos_x[i]/lado);  // Efecto mario bros
		pos_y[i] = pos_y[i] + vel_y[i]*h+(0.5)*f_y_t[i]*h*h;
		pos_y[i] = pos_y[i] -  lado* (int)(pos_y[i]/lado); // Efecto mario bros
		pos_z[i] = pos_z[i] + vel_z[i]*h+(0.5)*f_z_t[i]*h*h;
		pos_z[i] = pos_z[i] -  lado* (int)(pos_z[i]/lado);  // Efecto mario bros
	}
	
	for (int i=0; i < n; i++){
		f_x_t_h[i] = 0.0;
		f_y_t_h[i] = 0.0;
		f_z_t_h[i] = 0.0;
	}

	//calculo las nuevas fuerzas sobre cada una de las partículas (en el paso t+h)
	d_corte_cuadrado = d_corte*d_corte;
	for (int i = 0; i < (n-1); i++){
		for (int j=i+1; j < n; j++){
			
			
			distancia = 0;
            dx = pos_x[i]-pos_x[j];
            dy = pos_y[i]-pos_y[j];
            dz = pos_z[i]-pos_z[j];
            distancia_cuadrado = dx*dx+dy*dy+dz*dz;
            
            //Decide si la partícula i interactúa con la j de forma real o virtual o no interactúa directamente
			if (distancia_cuadrado <= d_corte_cuadrado){
                distancia = sqrt(distancia_cuadrado);
			}else{
                //i está muy lejos de j real, calcula entonces en qué región debería hubicarse la partícula virtual.
                //Divide por cuadrantes la caja real y se fija en cada una de las direcciones a que cuadrante va a ir a parar la virtual
                if(abs(dx) > d_corte) dx -= lado*( (dx > 0) - (dx < 0) ); //Esto es el x_j = x_j + L*signo(dx)
                if(abs(dy) > d_corte) dy -= lado*( (dy > 0) - (dy < 0) ); 
                if(abs(dz) > d_corte) dz -= lado*( (dz > 0) - (dz < 0) );
                distancia_cuadrado = dx*dx+dy*dy+dz*dz;
                if (distancia_cuadrado <= d_corte_cuadrado) distancia = sqrt(distancia_cuadrado);                
			}
			
			if(distancia > 0){
				int a = (int) (distancia*nf/4.0) - 1.0;		
				actualizacion = fuerzas[a]*(pos_x[j]-pos_x[i])/distancia;
				f_x_t_h[i] += actualizacion;
				f_x_t_h[j] -= actualizacion;
				actualizacion = fuerzas[a]*(pos_y[j]-pos_y[i])/distancia;
				f_y_t_h[i] += actualizacion;
				f_y_t_h[j] -= actualizacion;
				actualizacion = fuerzas[a]*(pos_z[j]-pos_z[i])/distancia;
				f_z_t_h[i] += actualizacion;
				f_z_t_h[j] -= actualizacion;
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
