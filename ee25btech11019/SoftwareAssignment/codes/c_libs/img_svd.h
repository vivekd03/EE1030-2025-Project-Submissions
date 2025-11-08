#ifndef img_svd_h
#define img_svd_h
//helper
double **creat_mat(int m, int n);
void free_mat(double ** M, int m);
void mult_mat(double **A, double **B, double** C, int m, int p, int n);
void trans_mat(double **A, double **At, int m, int n);
double vec_dot(double *a, double*b, int n);
void vec_devide(double *a, int n, double k);
void vecmatMult(double **A, int m, int n , double*x, double *y);
//improvement in power iteration
void gramSchmidt(double *v,double **old_V, int n, int count);
//svd
void power_iteration(double**A, int m, int n, double *v, double *sing, double** old_V, int count, int maxI, double tolerence);
void topkSingular(double **A, int m, int n, int k, double **U, double **V, double *sing);
//image input output
double** read_image_grayscale(const char *filename, int *width, int *height);
void write_image_grayscale(const char *filename, double **matrix, int width,int height, int quality);
//frobenius error
double frobn(double **A,double **B, int m,int n);
#endif
