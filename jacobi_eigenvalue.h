#include "etching.h"
#pragma acc routine seq
void jacobi_eigenvalue( int n, double a[], int it_max, double v[], double d[], int &it_num, int &rot_num );
#pragma acc routine seq
void r8mat_diag_get_vector ( int n, double a[], double v[] );
#pragma acc routine seq
void r8mat_identity ( int n, double a[] );
#pragma acc routine seq
double r8mat_is_eigen_right ( int n, int k, double a[], double x[], double lambda[] );
#pragma acc routine seq
double r8mat_norm_fro ( int m, int n, double a[] );
#pragma acc routine seq
void r8mat_print ( int m, int n, double a[], std::string title );
#pragma acc routine seq
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi, int jhi, std::string title );
#pragma acc routine seq
void r8vec_print ( int n, double a[], std::string title );
#pragma acc routine seq
void timestamp ( );
