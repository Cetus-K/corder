//mat.cpp

#include <header.h>
#include "vec.h"
#include "mat.h"

#define alloc_err "ERROR: Bad allocation has occured."
#define mdim_eq_err "ERROR: The matrices dimension must be equivalent."
#define mdim_def_err "ERROR: The matrix dimension must be larger than 0."
#define mdim_uni_err "ERROR: The matrix must be set square."

#define nan numeric_limits < double > ::quiet_NaN()
#define inf numeric_limits < double > ::infinity()

// constructor and destructor
mat::mat ( int vrt, int hrz ) {
	// check dimension
	if ( vrt < 1 || hrz < 1 ) {
		cerr << mdim_def_err << endl;
		exit(1);
	}
	// set dimension
	this->row = vrt;
	this->col = hrz;
	try {
		this->x = new double* [this->row];
		this->x[0] = new double [(this->row)*(this->col)];
		for ( int i=0; i<(this->row); i++ ) {
			this->x[i] = this->x[0] + i*(this->col);
			for ( int j=0; j<(this->col); j++ ) this->x[i][j] = 0;
		}
	} catch ( bad_alloc ) {
		// error code = bad_alloc
		cerr << alloc_err << endl;
		exit(1);
	}
}
mat::mat ( const mat &m ) {
	// set dimension
	this->row = m.row;
	this->col = m.col;
	try {
		this->x = new double* [this->row];
		this->x[0] = new double [(this->row)*(this->col)];
		for ( int i=0; i<(this->row); i++ ) {
			this->x[i] = this->x[0] + i*(this->col);
			for ( int j=0; j<(this->col); j++ ) this->x[i][j] = m.x[i][j];
		}
	} catch ( bad_alloc ) {
		// error code = bad_alloc
		cerr << alloc_err << endl;
		exit(1);
	}
}
mat::~mat() {
	if ( this->x ) {
		delete [] this->x[0];
		delete [] this->x;
		this->x = NULL;
	}
}
// size
void mat::resize ( int vrt, int hrz ) {
	// check dimension
	if ( vrt < 1 || hrz < 1 ) {
		cerr << mdim_def_err << endl;
		exit(1);
	}
	// free
	if ( this->x ) {
		delete [] this->x[0];
		delete [] this->x;
	}
	// set dimension
	this->row = vrt;
	this->col = hrz;
	try {
		this->x = new double* [this->row];
		this->x[0] = new double [(this->row)*(this->col)];
		for ( int i=0; i<(this->row); i++ ) {
			this->x[i] = this->x[0] + i*(this->col);
			for ( int j=0; j<(this->col); j++ ) this->x[i][j] = 0;
		}
	} catch ( bad_alloc ) {
		// error code = bad_alloc
		cerr << alloc_err << endl;
		exit(1);
	}
}
// one-matrix operator
void mat::zero ( void ) {
	for ( int i=0; i<(this->row); i++ ) {
	for ( int j=0; j<(this->col); j++ ) this->x[i][j] = 0;
	}
}
void mat::disp ( void ) {
	for ( int i=0; i<(this->row); i++ ) {
		for ( int j=0; j<(this->col); j++ ) {
			if ( 0 < this->x[i][j] ) cout << "\t " << this->x[i][j];
			else cout << "\t" << this->x[i][j];
		}
		cout << endl;
	}
	cout << endl;
}
void mat::unit ( void ) {
	if ( this->row != this->col ) {
		cerr << mdim_uni_err << endl;
		exit(0);
	}
	resize(this->row,this->col);
	for ( int i=0; i<(this->row); i++ ) this->x[i][i] = 1.0;
}
double mat::det2d ( void ) { return x[0][0]*x[1][1]-x[0][1]*x[1][0]; }
double mat::det3d ( void ) {
	double det3=0;
	det3 += x[0][0]*(x[1][1]*x[2][2]-x[1][2]*x[2][1]);
	det3 -= x[0][1]*(x[1][0]*x[2][2]-x[1][2]*x[2][0]);
	det3 += x[0][2]*(x[1][0]*x[2][1]-x[1][1]*x[2][0]);
	return det3;
}
mat mat::transpose ( void ) {
	mat m(this->col,this->row);
	for ( int i=0; i<(this->col); i++ ) {
	for ( int j=0; j<(this->row); j++ ) m.x[i][j] = this->x[j][i];
	}
	return m;
}
mat mat::inv2d ( void ) {
	double det2=det2d();
	mat inv2(2,2);
	if ( det2 ) {
		inv2.x[0][0] = x[1][1]; inv2.x[0][1] = -x[0][1];
		inv2.x[1][0] = -x[1][0]; inv2.x[1][1] = x[0][0];
		inv2 *= 1.0/det2;
	} else {
		inv2.x[0][0] = inf; inv2.x[0][1] = inf;
		inv2.x[1][0] = inf; inv2.x[1][1] = inf;
	}
	return inv2;
}
mat mat::inv3d ( void ) {
	double det3=det3d();
	mat inv3(3,3);
	if ( det3 ) {
		inv3.x[0][0] = x[1][1]*x[2][2]-x[1][2]*x[2][1];
		inv3.x[0][1] = x[0][2]*x[2][1]-x[0][1]*x[2][2];
		inv3.x[0][2] = x[0][1]*x[1][2]-x[0][2]*x[1][1];
		inv3.x[1][0] = x[1][2]*x[2][0]-x[1][0]*x[2][2];
		inv3.x[1][1] = x[0][0]*x[2][2]-x[0][2]*x[2][0];
		inv3.x[1][2] = x[0][2]*x[1][0]-x[0][0]*x[1][2];
		inv3.x[2][0] = x[1][0]*x[2][1]-x[1][1]*x[2][0];
		inv3.x[2][1] = x[0][1]*x[2][0]-x[0][0]*x[2][1];
		inv3.x[2][2] = x[0][0]*x[1][1]-x[0][1]*x[1][0];
		inv3 *= 1.0/det3;
	} else {
		inv3.x[0][0] = inf; inv3.x[0][1] = inf; inv3.x[0][2] = inf;
		inv3.x[1][0] = inf; inv3.x[1][1] = inf; inv3.x[1][2] = inf;
		inv3.x[2][0] = inf; inv3.x[2][1] = inf; inv3.x[2][2] = inf;
	}
	return inv3;
}
// operator
// substitution
mat mat::operator= ( const mat &m ) {
	if ( this->row != m.row || this->col != m.col ) resize(m.row,m.col);
	for ( int i=0; i<m.row; i++ ) {
	for ( int j=0; j<m.col; j++ ) this->x[i][j] = m.x[i][j];
	}
	return *this;
}
// addition
mat mat::operator+ ( const mat &m ) {
	// check dimension
	if ( this->row != m.row || this->col != m.col ) {
		cerr << mdim_eq_err << endl;
		exit(0);
	}
	// add matrix v
	mat n(this->row,this->col);
	for ( int i=0; i<(this->row); i++ ) {
	for ( int j=0; j<(this->col); j++ ) n.x[i][j] = this->x[i][j] + m.x[i][j];
	}
	return n;
}
mat mat::operator+= ( const mat &m ) {
	// check dimension
	if ( this->row != m.row || this->col != m.col ) {
		cerr << mdim_eq_err << endl;
		exit(0);
	}
	// add matrix v
	for ( int i=0; i<(this->row); i++ ) {
	for ( int j=0; j<(this->col); j++ ) this->x[i][j] += m.x[i][j];
	}
	return *this;
}
// subtraction
mat mat::operator- ( const mat &m ) {
	// check dimension
	if ( this->row != m.row || this->col != m.col ) {
		cerr << mdim_eq_err << endl;
		exit(0);
	}
	// subtract matrix v
	mat n(this->row,this->col);
	for ( int i=0; i<(this->row); i++ ) {
	for ( int j=0; j<(this->col); j++ ) n.x[i][j] = this->x[i][j] - m.x[i][j];
	}
	return n;
}
mat mat::operator-= ( const mat &m ) {
	// check dimension
	if ( this->row != m.row || this->col != m.col ) {
		cerr << mdim_eq_err << endl;
		exit(0);
	}
	// subtract matrix v
	for ( int i=0; i<(this->row); i++ ) {
	for ( int j=0; j<(this->col); j++ ) this->x[i][j] -= m.x[i][j];
	}
	return *this;
}
// ( inner ) production
vec mat::operator* ( const vec &v ) {
	// check dimension
	if ( this->col != v.dim ) {
		cerr << mdim_eq_err << endl;
		exit(0);
	}
	// product vector v
	vec w(this->row);
	for ( int i=0; i<(this->row); i++ ) {
		w[i] = 0;
		for ( int k=0; k<(this->col); k++ ) w[i] += this->x[i][k]*v.x[k];
	}
	return w;
}
vec mat::dot ( const vec &v ) {
	// check dimension
	if ( this->col != v.dim ) {
		cerr << mdim_eq_err << endl;
		exit(0);
	}
	// product vector v
	vec w(this->row);
	for ( int i=0; i<(this->row); i++ ) {
		w[i] = 0;
		for ( int k=0; k<(this->col); k++ ) w[i] += this->x[i][k]*v.x[k];
	}
	return w;
}
// ( outer ) production
mat mat::operator* ( const mat &m ) {
	// check dimension
	if ( this->col != m.row ) {
		cerr << mdim_eq_err << endl;
		exit(0);
	}
	// product matrix m
	mat n(this->row,m.col);
	for ( int i=0; i<(this->row); i++ ) {
	for ( int j=0; j<m.col; j++ ) {
		n[i][j] = 0;
		for ( int k=0; k<(this->col); k++ ) n[i][j] += this->x[i][k]*m.x[k][j];
	}
	}
	return n;
}
mat mat::dot ( const mat &m ) {
	// check dimension
	if ( this->col != m.row ) {
		cerr << mdim_eq_err << endl;
		exit(0);
	}
	// product matrix m
	mat n(this->row,m.col);
	for ( int i=0; i<(this->row); i++ ) {
	for ( int j=0; j<m.col; j++ ) {
		n[i][j] = 0;
		for ( int k=0; k<(this->col); k++ ) n[i][j] += this->x[i][k]*m.x[k][j];
	}
	}
	return n;
}
// ( const ) production
mat mat::operator* ( double c ) {
	mat n(this->row,this->col);
	for ( int i=0; i<(this->row); i++ ) {
	for ( int j=0; j<(this->col); j++ ) n[i][j] = c * this->x[i][j];
	}
	return n;
}
mat mat::operator*= ( double c ) {
	for ( int i=0; i<(this->row); i++ ) {
	for ( int j=0; j<(this->col); j++ ) this->x[i][j] *= c;
	}
	return *this;
}
mat operator* ( const mat &m, double c ) {
	mat n(m.row,m.col);
	for ( int i=0; i<(m.row); i++ ) {
	for ( int j=0; j<(m.col); j++ ) n[i][j] = c * m.x[i][j];
	}
	return n;
}
mat operator* ( double c, const mat &m ) {
	mat n(m.row,m.col);
	for ( int i=0; i<(m.row); i++ ) {
	for ( int j=0; j<(m.col); j++ ) n[i][j] = c * m.x[i][j];
	}
	return n;
}
