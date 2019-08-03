// function.cpp

#include "header.h"
#include "function.h"

// num to string
template < typename T > string to_str ( const T &t ) {
	ostringstream os;
	os << t;
	return os.str();
}
// memory alloc and free 1d/2d/3d
template < typename T > T *alloc1d ( int n, T *init_array=NULL ) {
	T *arr;
	try {
		arr = new T [n];
		for ( int i=0; i<n; i++ ) {
			if ( init_array ) arr[i] = init_array[i];
			else arr[i] = 0;
		}
	} catch ( ... ) { arr = NULL; }
	return arr;
}
template < typename T > T **alloc2d ( int m, int n, T **init_array=NULL ) {
	T **arr;
	try {
		arr = new T* [m];
		arr[0] = new T [m*n];
		for ( int i=0; i<m; i++ ) { arr[i] = arr[0] + i*n;
		for ( int j=0; j<n; j++ ) {
			if ( init_array ) arr[i][j] = init_array[i][j];
			else arr[i][j] = 0;
		}
		}
	} catch ( ... ) { arr = NULL; }
	return arr;
}
template < typename T > T ***alloc3d ( int l, int m, int n, T ***init_array=NULL ) {
	T ***arr;
	try {
		arr = new T** [l];
		arr[0] = new T* [l*m];
		arr[0][0] = new T [l*m*n];
		for ( int i=0; i<l; i++ ) { arr[i] = arr[0] + i*m;
		for ( int j=0; j<m; j++ ) { arr[i][j] = arr[0][0] + i*m*n + j*n;
		for ( int k=0; k<n; k++ ) {
			if ( init_array ) arr[i][j][k] = init_array[i][j][k];
			else arr[i][j][k] = 0;
		}
		}
		}
	} catch ( ... ) { arr = NULL; }
	return arr;
}
template < typename T > T ****alloc4d ( int n1, int n2, int n3, int n4, T ****init_array=NULL ) {
	T ****arr;
	try {
		arr = new T*** [n1];
		arr[0] = new T** [n1*n2];
		arr[0][0] = new T* [n1*n2*n3];
		arr[0][0][0] = new T [n1*n2*n3*n4];
		for ( int i=0; i<n1; i++ ) { arr[i] = arr[0] + i*n2;
		for ( int j=0; j<n2; j++ ) { arr[i][j] = arr[0][0] + i*n2*n3 + j*n3;
		for ( int k=0; k<n3; k++ ) {
			arr[i][j][k] = arr[0][0][0] + i*n2*n3*n4 + j*n3*n4 + k*n4;
			for ( int l=0; l<n4; l++ ) {
				if ( init_array ) arr[i][j][k] = init_array[i][j][k];
				else arr[i][j][k] = 0;
			}
		}
		}
		}
	} catch ( ... ) { arr = NULL; }
	return arr;
}
template < typename T > T *free1d ( T *arr ) {
	if ( arr ) {
		delete [] arr;
		arr = NULL;
	}
	return NULL;
}
template < typename T > T **free2d ( T **arr ) {
	if ( arr ) {
		delete [] arr[0];
		delete [] arr;
		arr = NULL;
	}
	return NULL;
}
template < typename T > T ***free3d ( T ***arr ) {
	if ( arr ) {
		delete [] arr[0][0];
		delete [] arr[0];
		delete [] arr;
		arr = NULL;
	}
	return NULL;
}
template < typename T > T ****free4d ( T ****arr ) {
	if ( arr ) {
		delete [] arr[0][0][0];
		delete [] arr[0][0];
		delete [] arr[0];
		delete [] arr;
		arr = NULL;
	}
	return NULL;
}
// grep -c command
int grep_c ( string filename, string keyword ) {
	int cnt=0;
	string str;
	ifstream ifs(filename.c_str(),ios::in);
	if ( ifs.is_open() ) {
	while ( getline(ifs,str) ) {
		if ( !str.find(keyword) ) cnt ++;
	}
	}
	return cnt;
}
// random function in [min,max]
double fixed_rand ( double min, double max ) {
	return min+(double)rand()*(max-min)/(1.0+(double)RAND_MAX);
}
// fitting integer
int fit_int ( double x ) {
	int xp,xm,fit_x;
	double dp,dm;
	xm = (int)x; xp = xm + 1;
	dm = abs(x-xm); dp = abs(xp-x);
	if ( dp < dm ) fit_x = xp;
	else if ( dp == dm ) fit_x = xp;
	else fit_x = xm;
	return fit_x;
}
// sign function
double sgn ( double x ) {
	if ( x != 0 ) return x / abs(x);
	else return 0;
}

// Kronecker delta
double Kronecker_delta ( int i, int j ) {
	if ( i == j ) return 1.0;
	else return 0;
}
// factorial
double fact ( int k ) {
	if ( k == 0 ) return 1.0;
	else return (double)k * fact(k-1);
}
// combinatorial
double comb ( int n, int k ) {
	return fact(n) / ( fact(n-k)*fact(k) );
}
// minimum and maximum
int min ( int x, int y, int z ) {
	int xymin = min(x,y);
	if ( z < xymin ) return z;
	else return xymin;
}
int max ( int x, int y, int z ) {
	int xymax = max(x,y);
	if ( z > xymax ) return z;
	else return xymax;
}
double min ( double x, double y, double z ) {
	double xymin = min(x,y);
	if ( z < xymin ) return z;
	else return xymin;
}
double max ( double x, double y, double z ) {
	double xymax = max(x,y);
	if ( z > xymax ) return z;
	else return xymax;
}
// Lambert W function ( real )
double Lambert_W0 ( double x ) {
	double w, D, eps;
	w = 0.0; D = 1.0; eps = 0.0001;
	while ( abs(D) > eps ) {
		D = 2.0*w+2.0;
		D = (w+2)*(w*exp(w)-x)/D;
		D = (w+1)*exp(w)-D;
		D = (w*exp(w)-x)/D;
		w -= D;
	}
	return w;
}
// cartesian to spherical
void cart2shpr ( double *r, double *theta, double *phi,
		    double x, double y, double z ){
	*r = sqrt( x*x + y*y + z*z );
	if ( *r == 0 ) *theta = 0;
	else *theta = acos( z/(*r) );
	*phi = atan2(y,x);
}
// n-order differential of z^m
double diff_pow ( double z, int m, int n ) {
	if ( m < n ) return 0;
	else if ( m == n ) return fact(m);
	else return fact(m)*pow(z,(double)m-n)/fact(m-n);
}
// Legendre polynomials
double Legendre_Pl ( int l, double z ) {
	int m;
	double p0,p1,pl;
	p0 = 1.0; p1 = z; pl = 0;
	if ( l == 0 ) return p0;
	else if ( l == 1 ) return p1;
	else {
		for ( m=1; m<l; m++ ) {
			pl = ((double)(2*m+1)*z*p1-(double)m*p0)/(double)(m+1);
			p0 = p1; p1 = pl;
		}
		return pl;
	}
}
// associated Legendre function
double Legendre_Plm ( int l, int m, double z ) {
	int k; double p;
	p = 0;
	for ( k=0; k<=l; k++ ) p += pow(-1.0,k)*comb(l,k)*diff_pow(z,2*l-2*k,l+abs(m));
	p *= pow(1.0-z*z,(double)abs(m)*0.5)/(pow(2.0,(double)l)*fact(l));
	if ( m < 0 ) return pow(-1.0,(double)abs(m))*fact(l-abs(m))*p/fact(l+abs(m));
	else return p;
}
// associated Laguerre polynomials
double Laguerre_Lnl ( int n, int l, double z ) {
	int k; double lg;
	lg = 0;
	for ( k=0; k<=n-l; k++ ) {
		lg += pow(-1.0,(double)l+k)*pow(fact(n),2.0)*pow(z,(double)k)/( fact(k)*fact(l+k)*fact(n-l-k) );
	}
	return lg;
}
// spherical harmonics
complex<double> spherical_Ylm ( int l, int m, double theta, double phi ) {
	return sqrt( (double)(2*l+1)*fact(l-m)/(4.0*M_PI*fact(l+m)) )*Legendre_Plm(l,m,cos(theta))*polar(1.0,(double)m*phi);
}
// spherical Bessel function
double spherical_Bessel_jl ( int l, double z ) {
	int m;
	double j0,j1,jl;
	m = j0 = j1 = jl = 0;
	if ( z == 0 ) return pow(2.0*z,(double)l)*fact(l)/fact(2*l+1);
	else {
		j0 = sin(z)/z; j1 = (sin(z)-z*cos(z))/(z*z);
		if ( l == 0 ) return j0;
		else if ( l == 1 ) return j1;
		else {
			for ( m=1; m<l; m++ ) {
				jl = (double)(2*m+1)*j1/z-j0;
				j0 = j1; j1 = jl;
			}
			return jl;
		}
	}
}

// Wigner-3j symbol
double Wigner_3j ( int j1, int j2, int j3, int m1, int m2, int m3 ) {
	int k; double sum = 0, den, prod_part;
	if ( m1+m2+m3 == 0 && abs(j1-j2) <= j3 && j3 <= j1+j2 ) {
		prod_part = fact(j1+j2-j3)*fact(j1-j2+j3)*fact(-j1+j2+j3);
		prod_part *= fact(j1-m1)*fact(j1+m1)*fact(j2-m2)*fact(j2+m2)*fact(j3-m3)*fact(j3+m3);
		prod_part /= fact(j1+j2+j3+1);
		for ( k=max(0,j2-j3-m1,j1-j3+m2); k<=min(j1+j2-j3,j1-m1,j2+m2); k++ ) {
			den = fact(k)*fact(j1+j2-j3-k)*fact(j1-m1-k)*fact(j2+m2-k)*fact(j3-j2+m1+k)*fact(j3-j1-m2+k);
			sum += pow(-1.0,(double)k)/den;
		}
		return pow(-1.0,(double)(j1-j2-m3))*sqrt(prod_part)*sum;
	}
	else return 0;
}

// Levi-Civita symbol ( double )
double LeviCivita_eps ( vector < int > ids ) {
	int idx,len;
	map < int,int > sgm;
	vector < int > pmt;
	vector < vector < int > > pmt2;
	for ( int i=0; i<(int)ids.size(); i++ ) sgm[i] = ids[i];
	while ( !sgm.empty() ) {
		idx = sgm.begin()->first;
		pmt.push_back(idx);
		while ( sgm.begin()->first != sgm[idx] ) {
			pmt.push_back(sgm[idx]);
			idx = sgm[idx];
		}
		for ( int i=0; i<(int)pmt.size(); i++ ) sgm.erase(pmt[i]); 
		pmt2.push_back(pmt);
		pmt.clear();
	}
	len = 0;
	for ( int i=0; i<(int)pmt2.size(); i++ ) {
		len += pmt2.size()-1;
	}
	vector < int > ().swap(pmt);
	vector < vector < int > > ().swap(pmt2);
	return pow(-1.0,(double)len);
}

// mpi function
void mpi_local_part ( int *local_start, int *local_size, int size, int rank, int nmpi ) {
	int ls0 = (int)((double)size/(double)nmpi);
	int res = size-nmpi*ls0;
	if ( 0 < res ) {
		if ( rank < res ) {
			*local_start = (ls0+1) * rank;
			*local_size = ls0 + 1;
		} else {
			*local_start = (ls0+1)*res+ls0*(rank-res);
			*local_size = ls0;
		}
	} else {
		*local_start = ls0 * rank;
		*local_size = ls0;
	}
}
// re-declare function explicitly
// int
template int *alloc1d < int > ( int n, int *init_array );
template int **alloc2d < int > ( int m, int n, int **init_array );
template int ***alloc3d < int > ( int l, int m, int n, int ***init_array );
template int ****alloc4d < int > ( int n1, int n2, int n3, int n4, int ****init_array );
template int *free1d < int > ( int *arr );
template int **free2d < int > ( int **arr );
template int ***free3d < int > ( int ***arr );
template int ****free4d < int > ( int ****arr );
// double
template double *alloc1d < double > ( int n, double *init_array );
template double **alloc2d < double > ( int m, int n, double **init_array );
template double ***alloc3d < double > ( int l, int m, int n, double ***init_array );
template double ****alloc4d < double > ( int n1, int n2, int n3, int n4, double ****init_array );
template double *free1d < double > ( double *arr );
template double **free2d < double > ( double **arr );
template double ***free3d < double > ( double ***arr );
template double ****free4d < double > ( double ****arr );
// complex < int >
template complex < int > *alloc1d < complex < int > > ( int n, complex < int > *init_array );
template complex < int > **alloc2d < complex < int > > ( int m, int n, complex < int > **init_array );
template complex < int > ***alloc3d < complex < int > > ( int l, int m, int n, complex < int > ***init_array );
template complex < int > ****alloc4d < complex < int > > ( int n1, int n2, int n3, int n4, complex < int > ****init_array );
template complex < int > *free1d < complex < int > > ( complex < int > *arr );
template complex < int > **free2d < complex < int > > ( complex < int > **arr );
template complex < int > ***free3d < complex < int > > ( complex < int > ***arr );
template complex < int > ****free4d < complex < int > > ( complex < int > ****arr );
// complex < double >
template complex < double > *alloc1d < complex < double > > ( int n, complex < double > *init_array );
template complex < double > **alloc2d < complex < double > > ( int m, int n, complex < double > **init_array );
template complex < double > ***alloc3d < complex < double > > ( int l, int m, int n, complex < double > ***init_array );
template complex < double > ****alloc4d < complex < double > > ( int n1, int n2, int n3, int n4, complex < double > ****init_array );
template complex < double > *free1d < complex < double > > ( complex < double > *arr );
template complex < double > **free2d < complex < double > > ( complex < double > **arr );
template complex < double > ***free3d < complex < double > > ( complex < double > ***arr );
template complex < double > ****free4d < complex < double > > ( complex < double > ****arr );
