// function.h
#ifndef FUNCTION_H
#define FUNCTION_H

// mathematical constant
const double pi = M_PI;
const complex < double > I = complex<double>(0, 1.0);
// system constant
const int int_max = 2147483647;
const int int_min = 0-2147483647-1;
// capacity
const double over_capacity = 1.0e15;
const double under_capacity = 1.0e-15;
// illegal constant
#define nan numeric_limits < double > ::quiet_NaN()
#define inf numeric_limits < double > ::infinity()
// declare function
template < typename T > string to_str ( const T &t );
template < typename T > T *alloc1d ( int n, T *init_array=NULL );
template < typename T > T **alloc2d ( int m, int n, T **init_array=NULL );
template < typename T > T ***alloc3d ( int l, int m, int n, T ***init_array=NULL );
template < typename T > T *free1d ( T *arr );
template < typename T > T **free2d ( T **arr );
template < typename T > T ***free3d ( T ***arr );
int grep_c ( string filename, string keyword );
double fixed_rand ( double min, double max );
int fit_int ( double x );
double sgn ( double x );
double Kronecker_delta ( int i, int j );
double fact ( int k );
double comb ( int n, int k );
int min ( int x, int y, int z );
int max ( int x, int y, int z );
double min ( double x, double y, double z );
double max ( double x, double y, double z );
double Lambert_W0 ( double x );
void cart2shpr ( double *r, double *theta, double *phi,
		    double x, double y, double z );
double diff_pow ( double z, int m, int n );
double Legendre_Pl ( int l, double z );
double Legendre_Plm ( int l, int m, double z );
double Laguerre_Lnl ( int n, int l, double z );
complex<double> spherical_Ylm ( int l, int m, double theta, double phi );
double spherical_Bessel_jl ( int l, double z );
double Wigner_3j ( int j1, int j2, int j3, int m1, int m2, int m3 );
double LeviCivita_eps ( vector < int > ids );
// mpi function
void mpi_local_part ( int *local_start, int *local_size, int size, int rank, int nmpi );

#endif
