#include "potencial.h"
#include "stdio.h"
#include <stdlib.h>
#include <math.h>

float rin = 1.8;
float rout = 2.0;
float sigma = 1;


void potencial(float *vector_x,float *vector_y,float *vector_z, float *vector_potencial,int n){
	float deltax,deltay,deltaz,r,s;
	for (int i = 0; i < n; ++i){
		for (int j = 0; j < n; ++i){
			deltax = vector_x[i] - vector_x[j];
			deltay = vector_y[i] - vector_y[j];
			deltaz = vector_z[i] - vector_z[j];
			r = sqrt(deltax*deltax + deltay*deltay + deltaz*deltaz);
			if (r<rin && r!=0){
				vector_potencial[i]+=4*(pow(sigma/r,12) - pow(sigma/r,6));
			}
			else if(rin<r && r<rout){  
				 s = switching_function(r,rin,rout);
				 vector_potencial[i]+=s*4*(pow(sigma/r,12) - pow(sigma/r,6));
			}
			else{
				vector_potencial[i]+= 0;	
			}
		} 
	}
	return;
}

// Esta es una funcion de spline.
float switching_function(float r,float rin,float rout){
	 float numerador1, numerador2,denominador,resultado;
	 numerador1 = (rout*rout-r*r)*(rout*rout-r*r);
	 numerador2 = rout*rout+2*r*r-3*rin*rin;
	 denominador = (rout*rout - rin*rin)*(rout*rout - rin*rin)*(rout*rout - rin*rin);   
	 resultado = numerador1*numerador2/denominador;
	return resultado;
}