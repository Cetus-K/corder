// vec.cpp

#include <header.h>
#include "vec.h"

#define alloc_err "ERROR: Bad allocation has occured."
#define vdim_eq_err "ERROR: The vectors dimension must be equivalent."
#define vdim_def_err "ERROR: The vector dimension must be larger than 0."

// constructor and destructor
vec::vec ( int len ) {
	// check dimension
	if ( len < 1 ) {
		cerr << vdim_def_err << endl;
		exit(1);
	}
	// set dimension
	this->dim = len;
	try {
		this->x = new double [this->dim];
		for ( int i=0; i<(this->dim); i++ ) this->x[i] = 0;
	} catch ( ... ) {
		// error code = bad_alloc
		cerr << alloc_err << endl;
		exit(1);
	}
}
vec::vec ( const vec &v ) {
	// set dimension
	this->dim = v.dim;
	try {
		this->x = new double [this->dim];
		for ( int i=0; i<(this->dim); i++ ) this->x[i] = v.x[i];
	} catch ( ... ) {
		// error code = bad_alloc
		cerr << alloc_err << endl;
		exit(1);
	}
}
vec::~vec() {
	if ( this->x ) {
		delete [] this->x;
		this->x = NULL;
	}
}
// size
void vec::resize ( int len ) {
	// check dimension
	if ( len < 1 ) {
		cerr << vdim_def_err << endl;
		exit(1);
	}
	// free
	if ( this->x ) delete [] this->x;
	// set dimension
	this->dim = len;
	try {
		this->x = new double [this->dim];
		for ( int i=0; i<(this->dim); i++ ) this->x[i] = 0;
	} catch ( ... ) {
		// error code = bad_alloc
		cerr << alloc_err << endl;
		exit(1);
	}
}
// one-vector operator
void vec::zero ( void ) {
	for ( int i=0; i<(this->dim); i++ ) this->x[i] = 0;
}
void vec::disp ( void ) {
	for ( int i=0; i<(this->dim); i++ ) {
		if ( 0 < this->x[i] ) cout << "\t " << this->x[i];
		else cout << "\t" << this->x[i];
	}
	cout << endl << endl;
}
double vec::norm ( void ) {
	double c = 0;
	for ( int i=0; i<(this->dim); i++ ) c += pow(this->x[i],2.0);
	return sqrt(c);
}
// operator
// substitution
vec vec::operator= ( const vec &v ) {
	if ( this->dim != v.dim ) resize(v.dim);
	for ( int i=0; i<v.dim; i++ ) this->x[i] = v.x[i];
	return *this;
}
// addition
vec vec::operator+ ( const vec &v ) {
	// check dimension
	if ( this->dim != v.dim ) {
		cerr << vdim_eq_err << endl;
		exit(0);
	}
	// add vector v
	vec w(this->dim);
	for ( int i=0; i<(this->dim); i++ ) w.x[i] = this->x[i] + v.x[i];
	return w;
}
vec vec::operator+= ( const vec &v ) {
	// check dimension
	if ( this->dim != v.dim ) {
		cerr << vdim_eq_err << endl;
		exit(0);
	}
	// add vector v
	for ( int i=0; i<(this->dim); i++ ) this->x[i] += v.x[i];
	return *this;
}
// subtraction
vec vec::operator- ( const vec &v ) {
	// check dimension
	if ( this->dim != v.dim ) {
		cerr << vdim_eq_err << endl;
		exit(0);
	}
	// subtract vector v
	vec w(this->dim);
	for ( int i=0; i<(this->dim); i++ ) w.x[i] = this->x[i] - v.x[i];
	return w;
}
vec vec::operator-= ( const vec &v ) {
	// check dimension
	if ( this->dim != v.dim ) {
		cerr << vdim_eq_err << endl;
		exit(0);
	}
	// subtract vector v
	for ( int i=0; i<(this->dim); i++ ) this->x[i] -= v.x[i];
	return *this;
}
// ( inner ) production
double vec::operator* ( const vec &v ) {
	// check dimension
	if ( this->dim != v.dim ) {
		cerr << vdim_eq_err << endl;
		exit(0);
	}
	// product scalar c
	double c = 0;
	for ( int i=0; i<(this->dim); i++ ) c += this->x[i] * v.x[i];
	return c;
}
double vec::dot ( const vec &v ) {
	// check dimension
	if ( this->dim != v.dim ) {
		cerr << vdim_eq_err << endl;
		exit(0);
	}
	// product scalar c
	double c = 0;
	for ( int i=0; i<(this->dim); i++ ) c += this->x[i] * v.x[i];
	return c;
}
// ( outer ) production
// assuming: cartesian coordination ( det(g) = 1 )
// assuming: 3-dimensional vector
vec vec::cross ( const vec &v ) {
	// check dimension
	if ( this->dim != v.dim ) {
	if ( this->dim == 3 ) {
		cerr << vdim_eq_err << endl;
		exit(0);
	}
	}
	// product vector v
	vec w(this->dim);
	w.x[0] = this->x[1]*v.x[2]-this->x[2]*v.x[1];
	w.x[1] = this->x[2]*v.x[0]-this->x[0]*v.x[2];
	w.x[2] = this->x[0]*v.x[1]-this->x[1]*v.x[0];
	return w;
}
// ( const ) production
vec vec::operator* ( double c ) {
	vec w(this->dim);
	for ( int i=0; i<(this->dim); i++ ) w[i] = c * this->x[i];
	return w;
}
vec vec::operator*= ( double c ) {
	for ( int i=0; i<(this->dim); i++ ) this->x[i] *= c;
	return *this;
}
vec operator* ( const vec &v, double c ) {
	vec w(v.dim);
	for ( int i=0; i<v.dim; i++ ) w[i] = c * v.x[i];
	return w;
}
vec operator* ( double c, const vec &v ) {
	vec w(v.dim);
	for ( int i=0; i<v.dim; i++ ) w[i] = c * v.x[i];
	return w;
}
