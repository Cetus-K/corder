// vec.h
#ifndef VEC_H
#define VEC_H

#include <header.h>

#define alloc_err "ERROR: Bad allocation has occured."
#define vdim_eq_err "ERROR: The vectors dimension must be equivalent."
#define vdim_def_err "ERROR: The vector dimension must be larger than 0."

// vector class
class vec {

	public:
		int dim;
		double *x;

	public:	// size
		int size(){ return dim; }
		void resize ( int len );

	public:	// a one-vector operator
		void zero ( void );
		void disp ( void );
		double norm ( void );

	public:
		double &operator[]( int idx ){ return x[idx]; }
		double &operator()( int idx ){ return x[idx]; }
		double operator* ( const vec &v );
		double dot ( const vec &v );
		vec cross ( const vec &v );
		vec operator= ( const vec &v );
		vec operator+ ( const vec &v );
		vec operator- ( const vec &v );
		vec operator* ( double c );
		vec operator+= ( const vec &v );
		vec operator-= ( const vec &v );
		vec operator*= ( double c );
		friend vec operator* ( const vec &v, double c );
		friend vec operator* ( double c, const vec &v );

	public:
		vec ( int len=1 );		// default constructor
		vec ( const vec &v );	// copy constructor
		~vec();			// destructor

};

#endif
