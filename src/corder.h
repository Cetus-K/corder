// corder.h
#ifndef CORDER_H
#define CORDER_H

#include <header.h>
#include <function.h>
#include <vec.h>
#include <mat.h>
#include <mpi.h>

// using voro++ code
#include <voro++.hh>
using namespace voro;

// using qhull code
#include <libqhullcpp/Qhull.h>
#include <libqhull_r/qhull_ra.h>
using namespace orgQhull;

// system class ( basis class )
class corder {

	public:
		// constructor and destructor
		corder() { ; }
		~corder();

		/*--------------------------
			declare variable
		--------------------------*/

		// read param
		vector < string > mode;	// calc mode
		vector < int > l;		// angular l-number array
		int f0,f1,df;			// sampling flames
		double r0,r1,dr;		// r-range
		double q0,q1,dq;		// q-range
		double ql0,ql1,dql;		// ql-range
		double wl0,wl1,dwl;		// wl-range
		double cs0,cs1,dcs;		// cs-range
		double dt;			// time step
		int trans;			// super cell
		int lcut;			// cut-off angular l-number
		int dump;			// dumping interval
		vector < vector < string > > stype;	// specified types
		vector < vector < int > > stypeindex;	// and corresponding indices
		vector < vector < int > > stypelist;	// list of specified types
		vector < int > sitem;	// number of specified types
		int nstypes;			// number of specified types

		// atomic information ( read from XDATCAR )
		bool latdyn;			// switch lattice dynamics
		int nions,ntypes;		// number of ions/types/items
		int ncells;			// number of super cells
		vector < string > type;	// type; species string
		vector < int > item;		// item; number of thistype

		// iterative variable
		int itr;
		vec a1,a2,a3;
		vector < vec > config;
		mat latmat;

		/*------------------
			function
		------------------*/

		// read param.in
		void read_input ( void );
		// calculate fdtd parameter
		void read_system ( void );
		// set above 2 function
		void set_param ( void );
		// write system information
		void info ( void );

		// the type of i-th index ion
		int mytype ( int idx );
		// the address that start itype ion indices
		int mystart ( int itype );
		// update value in the system
		double update ( double prev, double add );
		double update_delay ( double prev, double add, int delay );

};

// gofr class
class gofr {

	private:
		corder *co;

	public:
		// constructor and destructor
		gofr() { co=NULL; }
		gofr ( corder &cdr ) { co = &cdr; }
		void set ( corder &cdr ) { co = &cdr; }
		~gofr();

		/*--------------------------
			declare variable
		--------------------------*/

		// average density: rho, rhoab, rhosb
		// total radial distribution function: g(r)
		// partial radial distribution function: gab(r)
		// specified partial radial distribution function: gsab(r)
		// total radial density distribution function: rho*g(r)
		// partial radial density distribution function: rhoab*gab(r)
		// specified partial radial density distribution function: rhosab*gsab(r)
		double rho,*rhob,*rhosb;
		double *gr,*rhor;
		double **gabr,**rhoabr;
		double **gsabr,**rhosabr;

		// radial grid
		int nr;

		/*------------------
			function
		------------------*/

		// initialize
		void init ( void );

		// update gofr
		void update ( int rank, int nmpi );

		// write out: g(r)
		void write ( string filename );

};

// sofq class
class sofq {

	private:
		corder *co;

	public:
		// constructor and destructor
		sofq() { co=NULL; }
		sofq ( corder &cdr ) { co = &cdr; }
		void set ( corder &cdr ) { co = &cdr; }
		~sofq() { ; }

		/*--------------------------
			declare variable
		--------------------------*/

		// reciprocal radial grid
		int nq;

		/*------------------
			function
		------------------*/

		// initialize
		void init ( void );

		// write out: s(q)
		void write ( string filename, gofr *gf );

};

// dfc class
class dfc {

	private:
		corder *co;

	public:
		// constructor and destructor
		dfc() { co=NULL; }
		dfc ( corder &cdr ) { co = &cdr; }
		void set ( corder &cdr );
		~dfc();

		/*--------------------------
			declare variable
		--------------------------*/

		// total mean-square displacement: MSD
		// partial mean-square displacement: MSDa
		// specified partial mean-square displacement: MSDsa
		// initial and previous configuration
		// time-averaged velocity: <v>, <va>, <vsa>
		double msd,*msda,*msdsa;
		vector < vec > config_init,config_prev;
		vec veloc;
		vector < vec > veloca,velocsa;

		/*------------------
			function
		------------------*/

		// initialize
		void init ( void );

		// update dfc
		void update ( void );

		// write out: dfc
		void write ( string filename );

};

// covoro class for lboo,bacf,qwdist,bofree,povoro
class covoro {

	private:
		corder *co;
		mat tfmat;	// input information into voro++

	public:
		// constructor and destructor
		covoro() { co=NULL; }
		covoro ( corder &cdr ) { co = &cdr; }
		~covoro();

		/*--------------------------
			declare variable
		--------------------------*/

		// output information from voro++
		vector < vec > cell;
		// total information
		vector < vec > config;
		vector < vector < vec > > nv;
		vector < vector < vector < vec > > > eg;
		vector < vector < int > > nnid;
		vector < vector < double > > surf;
		// partial information
		vector < vector < vector < vec > > > nvab;
		vector < vector < vector < vector < vec > > > > egab;
		vector < vector < vector < int > > > nnidab;
		vector < vector < vector < double > > > surfab;
		// specified partial information
		vector < vector < vector < vec > > > nvsab;
		vector < vector < vector < vector < vec > > > > egsab;
		vector < vector < vector < int > > > nnidsab;
		vector < vector < vector < double > > > surfsab;

		/*------------------
			function
		------------------*/

		// initialize
		void init ( void );

		// set corder copy constructor
		// and match lattice vector format to voro++
		void set ( corder &cdr );

		// get informations that voro++ calculated
		// get values: surface area, normal vector and near-neighbor ids
		void get ( void );
		void pget ( int a, int b );
		void spget ( int spa, int spb );

};

// lboo class
class lboo {

	private:
		corder *co;
		covoro *cv;

	public:
		// constructor and destructor
		lboo() { co=NULL; cv=NULL; }
		lboo ( corder &cdr ) { co=&cdr; cv=NULL; }
		lboo ( covoro &cvr ) { co=NULL; cv=&cvr; }
		lboo ( corder &cdr, covoro &cvr ) { co=&cdr; cv=&cvr; }
		void set ( corder &cdr ) { co = &cdr; }
		void set ( covoro &cvr ) { cv = &cvr; }
		void set ( corder &cdr, covoro &cvr ) { co=&cdr; cv=&cvr; }
		~lboo() { ; }

		/*--------------------------
			declare variable
		--------------------------*/

		// number of angular l-number
		int nl;

		/*------------------
			function
		------------------*/

		// initialize
		void init ( void );

		// write out: lboo
		void write ( string filename );

};

// bacf class
class bacf {

	private:
		corder *co;
		covoro *cv;

	public:
		// constructor and destructor
		bacf() { co=NULL; cv=NULL; }
		bacf ( corder &cdr ) { co=&cdr; cv=NULL; }
		bacf ( covoro &cvr ) { co=NULL; cv=&cvr; }
		bacf ( corder &cdr, covoro &cvr ) { co=&cdr; cv=&cvr; }
		void set ( corder &cdr ) { co = &cdr; }
		void set ( covoro &cvr ) { cv = &cvr; }
		void set ( corder &cdr, covoro &cvr ) { co=&cdr; cv=&cvr; }
		~bacf();

		/*--------------------------
			declare variable
		--------------------------*/

		// number of angular l-number
		int nl;
		// radial grid
		int nr;
		// G0(r)*Gl(r) data
		double **gglr,***gglababr,***gglsababr;

		/*------------------
			function
		------------------*/

		// initialize
		void init ( void );

		// update bacf
		void update ( int rank, int nmpi );

		// write out: bacf
		void write ( string filename );

};

// pofqw class
class pofqw {

	private:
		corder *co;
		covoro *cv;

	public:
		// constructor and destructor
		pofqw() { co=NULL; cv=NULL; }
		pofqw ( corder &cdr ) { co=&cdr; cv=NULL; }
		pofqw ( covoro &cvr ) { co=NULL; cv=&cvr; }
		pofqw ( corder &cdr, covoro &cvr ) { co=&cdr; cv=&cvr; }
		void set ( corder &cdr ) { co = &cdr; }
		void set ( covoro &cvr ) { cv = &cvr; }
		void set ( corder &cdr, covoro &cvr ) { co=&cdr; cv=&cvr; }
		~pofqw();

		/*--------------------------
			declare variable
		--------------------------*/

		// number of angular l-number
		int nl,nql,nwl;
		// p(ql),p(wl),p(qlab),p(wlab) data
		double **pql,**pwl,***pqlab,***pwlab,***pqlsab,***pwlsab;

		/*------------------
			function
		------------------*/

		// initialize
		void init ( void );

		// update pofqw
		void update ( void );

		// write out: pofqw
		void write ( string filename );

};

// bofree class
class bofree {

	private:
		corder *co;
		covoro *cv;

	public:
		// constructor and destructor
		bofree() { co=NULL; cv=NULL; }
		bofree ( corder &cdr ) { co=&cdr; cv=NULL; }
		bofree ( covoro &cvr ) { co=NULL; cv=&cvr; }
		bofree ( corder &cdr, covoro &cvr ) { co=&cdr; cv=&cvr; }
		void set ( corder &cdr ) { co = &cdr; }
		void set ( covoro &cvr ) { cv = &cvr; }
		void set ( corder &cdr, covoro &cvr ) { co=&cdr; cv=&cvr; }
		~bofree();

		/*--------------------------
			declare variable
		--------------------------*/

		// number of angular l-number
		// from 0 to lcut
		int nl;
		double *il,**jll,**ilab,***jllab,**ilsab,***jllsab;

		/*------------------
			function
		------------------*/

		// initialize
		void init ( void );

		// update bofree
		void update ( int rank, int nmpi );

		// write out: bofree
		void write ( string filename );

};

// povoro class
class povoro {

	private:
		corder *co;
		covoro *cv;

	public:
		// constructor and destructor
		povoro() { co=NULL; cv=NULL; }
		povoro ( corder &cdr ) { co=&cdr; cv=NULL; }
		povoro ( covoro &cvr ) { co=NULL; cv=&cvr; }
		povoro ( corder &cdr, covoro &cvr ) { co=&cdr; cv=&cvr; }
		void set ( corder &cdr ) { co = &cdr; }
		void set ( covoro &cvr ) { cv = &cvr; }
		void set ( corder &cdr, covoro &cvr ) { co=&cdr; cv=&cvr; }
		~povoro() { ; }

		/*------------------
			function
		------------------*/

		// write out: povoro
		void write ( string filename );
		void pwrite ( string filename, int a, int b );
		void spwrite ( string filename, int spa, int spb );

};

// coqhull class
class coqhull {

	private:
		corder *co;
		covoro *cv;

	public:
		// constructor and destructor
		coqhull() { co=NULL; cv=NULL; }
		coqhull ( corder &cdr ) { co=&cdr; cv=NULL; }
		coqhull ( covoro &cvr ) { co=NULL; cv=&cvr; }
		coqhull ( corder &cdr, covoro &cvr ) { co=&cdr; cv=&cvr; }
		void set ( corder &cdr ) { co = &cdr; }
		void set ( covoro &cvr ) { cv = &cvr; }
		void set ( corder &cdr, covoro &cvr ) { co=&cdr; cv=&cvr; }
		~coqhull();

		/*--------------------------
			declare variable
		--------------------------*/

		// nv: normal vectors of facet
		// eg: edges belonging facet
		// vt: vertices of facet
		// vnv: vertices which construct nv
		// veg: vertices which construct eg
		// total information
		vector < vector < vec > > nv;
		vector < vector < vec > > eg;
		vector < vector < vec > > vt;
		vector < vector < vector < int > > > vnv;
		vector < vector < vector < int > > > veg;
		// partial information
		vector < vector < vector < vec > > > nvab;
		vector < vector < vector < vec > > > egab;
		vector < vector < vector < vec > > > vtab;
		vector < vector < vector < vector < int > > > > vnvab;
		vector < vector < vector < vector < int > > > > vegab;
		// specified partial information
		vector < vector < vector < vec > > > nvsab;
		vector < vector < vector < vec > > > egsab;
		vector < vector < vector < vec > > > vtsab;
		vector < vector < vector < vector < int > > > > vnvsab;
		vector < vector < vector < vector < int > > > > vegsab;

		/*------------------
			function
		------------------*/

		// initialize
		void init ( void );

		// set corder copy constructor and
		// get informations that qhull calculated
		// get values: normal vectors of facets, edges and vertices
		void get_one ( vector < vec > points, int index );
		void pget_one ( vector < vec > points, int index, int a, int b );
		void spget_one ( vector < vec > points, int index, int spa, int spb );
		void get ( vector < vector < vec > > points );
		void pget ( vector < vector < vec > > points, int a, int b );
		void spget ( vector < vector < vec > > points, int spa, int spb );

};

// csform class
class csform {

	private:
		corder *co;
		covoro *cv;
		coqhull *cq;

	public:
		// constructor and destructor
		csform() { co=NULL; cv=NULL; cq=NULL; }
		csform ( corder &cdr ) { co=&cdr; cv=NULL; cq=NULL; }
		csform ( covoro &cvr ) { co=NULL; cv=&cvr; cq=NULL; }
		csform ( coqhull &cqh ) { co=NULL; cv=NULL; cq=&cqh; }
		csform ( corder &cdr, covoro &cvr ) { co=&cdr; cv=&cvr; cq=NULL; }
		csform ( corder &cdr, coqhull &cqh ) { co=&cdr; cv=NULL; cq=&cqh; }
		csform ( covoro &cvr, coqhull &cqh ) { co=NULL; cv=&cvr; cq=&cqh; }
		csform ( corder &cdr, covoro &cvr, coqhull &cqh ) { co=&cdr; cv=&cvr; cq=&cqh; }
		void set ( corder &cdr ) { co = &cdr; }
		void set ( covoro &cvr ) { cv = &cvr; }
		void set ( coqhull &cqh ) { cq = &cqh; }
		void set ( corder &cdr, covoro &cvr ) { co=&cdr; cv=&cvr; }
		void set ( corder &cdr, coqhull &cqh ) { co=&cdr; cq=&cqh; }
		void set ( covoro &cvr, coqhull &cqh ) { cv=&cvr; cq=&cqh; }
		void set ( corder &cdr, covoro &cvr, coqhull &cqh ) { co=&cdr; cv=&cvr; cq=&cqh; }
		~csform();

		/*--------------------------
			declare variable
		--------------------------*/

		// number of countable range
		int ncs,*ddim,**ddimab,**ddimsab;
		// total information
		double *pb,*pf,*pe,*pv;
		// partial information
		double **pbab,**pfab,**peab,**pvab;
		// specified partial information
		double **pbsab,**pfsab,**pesab,**pvsab;

		/*------------------
			function
		------------------*/

		// initialize
		void init ( void );

		// flag for calc qhull
		bool qhflag ( vector < vec > points );
		// check degeneration for convex hull
		int dim_convex ( vector < vec > points );
		// calc cyclic vertex list and center position of 2d-polygon
		void calc_flat ( vector < vec > points, vector < int > *list, vec *g );
		// for calculating the collision
		bool separation ( vector < vec > vti, vector < vec > vtj, vec gij, vec eij );
		// check whether the two x-dimensional convex hulls are in collision
		bool collision ( int i, int j, int iab, int isab );
		// get convex hull
		void get ( void );
		void pget ( void );
		void spget ( void );

		// update csform
		// count sharing constructures
		void update ( int rank, int nmpi );

		// write out: csform
		void write ( string filename );

};


#endif
