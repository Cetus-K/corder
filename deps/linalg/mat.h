// mat.h
#ifndef MAT_H
#define MAT_H

#include <header.h>
#include "vec.h"

#define alloc_err "ERROR: Bad allocation has occured."
#define mdim_eq_err "ERROR: The matrices dimension must be equivalent."
#define mdim_def_err "ERROR: The matrix dimension must be larger than 0."
#define mdim_uni_err "ERROR: The matrix must be set square."

// matrix class
class mat {

	public:
		int row,col;
		double **x;

	public:	// size
		int row_size(){ return row; }
		int col_size(){ return col; }
		void resize ( int vrt, int hrz );

	public:	// a one-matrix operator
		void zero ( void );
		void disp ( void );
		void unit ( void );
		double det2d ( void );
		double det3d ( void );
		mat transpose ( void );
		mat inv2d ( void );
		mat inv3d ( void );

	public:
		double* &operator[]( int i ){ return x[i]; }
		double &operator()( int i, int j ){ return x[i][j]; }
		vec operator* ( const vec &v );
		mat operator* ( const mat &m );
		vec dot ( const vec &v );
		mat dot ( const mat &m );
		mat operator= ( const mat &m );
		mat operator+ ( const mat &m );
		mat operator- ( const mat &m );
		mat operator* ( double c );
		mat operator+= ( const mat &m );
		mat operator-= ( const mat &m );
		mat operator*= ( double c );
		friend mat operator* ( const mat &m, double c );
		friend mat operator* ( double c, const mat &m );

	public:
		mat ( int vrt=1, int hrz=1 );	// default constructor
		mat ( const mat &m );		// copy constructor
		~mat();				// destructor

};

#endif
