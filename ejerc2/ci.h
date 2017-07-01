#ifndef CI_H
#define CI_H
double posiciones_iniciales(int n, float densidad,double *vector_x,double *vector_y,double *vector_z);
int velocidades_iniciales(float T_actual,int n,double *vector_x,double *vector_y,double *vector_z);
int fuerzas_iniciales(int n, double *f_x_t, double *f_y_t, double *f_z_t);
double rand_gausiano(float media,float sigma);
#endif
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
