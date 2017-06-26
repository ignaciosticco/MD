#include "tablas.h"
#include "stdio.h"
#include <stdlib.h>
#include <math.h>

double rin = 2.0;  //  r de corte interno
double rout = 4.0; //  r de corte interno (V=0 si r>r0ut)
double sigma = 1.0;

int tabla_potenciales(int nf, double *potenciales){

	double r, numerador1, numerador2, denominador, s;
	   
	for (int i=0; i < nf; i++){
			r = 4.0*(i+1)/(double)nf;
			if (r <= rin) {
				potenciales[i] = 4.0*(pow(sigma/r,12.)-pow(sigma/r,6.));
			}
			else if (rin < r && r < rout) {

/*Esta es una funcion de spline. Fuente: http://www4.ncsu.edu/~franzen/public_html/CH795N/lecture/III/trunc/trunc.html
Suaviza al potencial desde rin hasta rout.
*/
				numerador1 = (rout*rout-r*r)*(rout*rout-r*r);
				numerador2 = rout*rout+2.0*r*r-3.0*rin*rin;
				denominador = (rout*rout - rin*rin)*(rout*rout - rin*rin)*(rout*rout - rin*rin);
				s = (numerador1*numerador2)/denominador;
				potenciales[i] = s*4.0*(pow(sigma/r,12.) - pow(sigma/r,6.));
			}
			else {
				potenciales[i] = 0.0;	
			}
	}

return 0;
}


int tabla_fuerzas(int nf, double *fuerzas){

	double r, numerador1, numerador2, denominador, s, nuevo_numerador1, nuevo_numerador2, nuevo_denominador, nuevo_s;

	for (int i=0; i < nf; i++){
			r = 4.0*(i+1)/(double)nf;
			if (r <= rin) {
				fuerzas[i] = 48.0*(pow(sigma/r,12.)-0.5*pow(sigma/r,6.))/r;
			}
			else if (rin < r && r < rout) {
				//primero derivo s
				nuevo_numerador1 = -4.0*r*(rout*rout-r*r)*(rout*rout+2.0*r*r-3.0*rin*rin);
				nuevo_numerador2 = pow(rout*rout-r*r,2.)*4.0*r;
				nuevo_denominador = pow(rout*rout - rin*rin,3.);
				nuevo_s = (nuevo_numerador1 + nuevo_numerador2)/nuevo_denominador;

				//segundo derivo la componente de potencias
				numerador1 = pow(rout*rout-r*r,2.);
				numerador2 = rout*rout+2*r*r-3*rin*rin;
				denominador = pow(rout*rout - rin*rin,3.);
				s = (numerador1*numerador2)/denominador;

			       	fuerzas[i] = nuevo_s*4.0*(pow(sigma/r,12.)-pow(sigma/r,6.))+s*48.0*(pow(sigma/r,12.)-0.5*pow(sigma/r,6.))/r;
			}
			else {
				fuerzas[i] = 0.0;	
			}
	}

return 0;
}
