#ifndef MAGNITUDES_H
#define MAGNITUDES_H
int potencial(double *vector_potencial, double va, int a,int i, int j);
double cinetica(int n, double *vel_x, double *vel_y, double *vel_z, double *vector_cinetico);
double energia_potencial(int n,double *vector_potencial);
#endif

