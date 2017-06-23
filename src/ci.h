#ifndef CI_H
#define CI_H
double posiciones_iniciales(int n, float densidad,double *vector_x,double *vector_y,double *vector_z);
int velocidades_iniciales(int n,double *vector_x,double *vector_y,double *vector_z);
int fuerzas_iniciales(int n, double d_corte, double *pos_x, double *pos_y, double *pos_z, double *f_x_t, double *f_y_t, double *f_z_t, double *fuerzas, int nf);
#endif

