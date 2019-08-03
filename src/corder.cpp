#include "corder.h"

// corder destructor
corder::~corder() {
	vector < int > ().swap(l);
	vector < int > ().swap(item);
	vector < int > ().swap(sitem);
	vector < vector < int > > ().swap(stypeindex);
	vector < vector < int > > ().swap(stypelist);
	vector < string > ().swap(mode);
	vector < string > ().swap(type);
	vector < vector < string > > ().swap(stype);
	stypelist.clear();
}

// the type of i-th index ion
int corder::mytype ( int idx ) {
	int itype,cnt = 0;
	for ( itype=0; itype<ntypes; itype++ ) {
		cnt += item[itype];
		if ( idx < cnt ) break;
	}
	return itype;
}

// the address that start itype ion indices
// range: [ mystart(a), mystart(a)+item[a] )
int corder::mystart ( int itype ) {
	int cnt = 0;
	for ( int a=0; a<itype; a++ ) cnt += item[a];
	return cnt;
}

// update
double corder::update ( double prev, double add ) {
	return ( (double)(itr-f0)/df*prev+add )/((double)(itr-f0)/df+1.0);
}

// main control
int main(int argc, char* argv[]){

	// declare
	int rank,nmpi;		// mpi rank ans its size
	double scale;			// lattice scale
	double beg,time,prog;	// progress
	map < string,bool > mode;	// calc mode
	vec cfg(3);			// vector and matrix class
	ifstream ifs;			// fstream
	istringstream iss;
	string str;
	corder co;			// corder basis class

	// set corder parameter
	co.set_param();

	// calculation mode
	mode["gofr"] = false;	// g(r),gab(r)
	mode["sofq"] = false;	// s(q),sab(q)
	mode["dfc"] = false;		// msd,msda,<v>,<va>
	mode["lboo"] = false;	// ql,wl,qlab,wlab
	mode["bacf"] = false;	// gl(r),glabab(r)
	mode["pofqw"] = false;	// p(ql),p(wl),p(qlab),p(wlab)
	mode["bofree"] = false;	// f,fl,fll,fab,flab,fllab
	mode["povoro"] = false;	// dump trajectory as pov-ray format
	mode["csform"] = false;	// counting each sharing formations
	// dependence parameter
	mode["covoro"] = false;	// engine for lboo/bacf/pofqw/bofree/voro
	mode["coqhull"] = false;	// engine for ...
	mode["instant"] = false;	// instantsneous calculation for write out
	for ( int i=0; i<(int)co.mode.size(); i++ ) mode[co.mode[i]] = true;

	// dependency
	mode["sofq"] = mode["sofq"] and mode["gofr"];
	mode["covoro"] = mode["lboo"] or mode["bacf"] or mode["pofqw"] or mode["bofree"] or mode["csform"] or mode["povoro"] or mode["csform"];
	mode["instant"] = !( mode["bacf"] or mode["pofqw"] or mode["bofree"] or mode["csform"] );

	// check read
	if ( co.f1 < co.f0 ) co.f1 = grep_c("XDATCAR","Direct");

	// memory allocation:
	co.config.resize(co.nions);
	for ( int i=0; i<co.nions; i++ ) co.config[i].resize(3);

	// initialize
	// each class related to mode
	// create copy constructor by refference handling
	gofr co_gofr(co); co_gofr.init();
	sofq co_sofq(co); co_sofq.init();
	dfc co_dfc(co); co_dfc.init();
	covoro co_covoro(co); co_covoro.init();
	lboo co_lboo(co,co_covoro); co_lboo.init();
	bacf co_bacf(co,co_covoro); co_bacf.init();
	pofqw co_pofqw(co,co_covoro); co_pofqw.init();
	bofree co_bofree(co,co_covoro); co_bofree.init();
	povoro co_povoro(co,co_covoro);
	coqhull co_coqhull(co); co_coqhull.init();
	coqhull co_cv_coqhull(co,co_covoro); co_cv_coqhull.init();
	csform co_csform(co,co_covoro,co_cv_coqhull); co_csform.init();

	// MPI environment
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&nmpi);
	MPI_Barrier(MPI_COMM_WORLD);
	beg = prog = MPI_Wtime();

	// read XDATCAR and skip initial lines without axis at latdyn = false
	ifs.open("XDATCAR",ios::in);
	if ( !co.latdyn ) for ( int i=0; i<7; i++ ) getline(ifs,str);

	// dump system information
	if ( rank == 0 ) {
		co.info();
		cout << endl;
		cout << "------------------" << endl;
		cout << "     Progress" << endl;
		cout << "------------------" << endl;
		cout << endl;
		cout << "Iteration ( completion [%] ) : 1 dump-step cpu time ( total ) [s]" << endl;
		cout << "-----------------------------------------------------------------" << endl;
	}
	MPI_Barrier(MPI_COMM_WORLD);

	// iterator init
	co.itr = 0;
	while ( co.itr<=co.f1 && !ifs.eof() ) {
		co.itr ++;

		if ( ( co.f0<=co.itr && co.itr<=co.f1 && (co.itr-co.f0)%co.df==0 ) || co.itr==co.f1 ) {
			// if lattice dynamics is true
			// read lattice system
			if ( co.latdyn ) {
				// label
				getline(ifs,str);
				// scale
				getline(ifs,str); iss.str(str);
				iss >> scale; iss.clear();
				// lattice vector >> transform matrix
				getline(ifs,str); iss.str(str);
				iss >> co.a1(0) >> co.a1(1) >> co.a1(2);
				iss.clear();  co.a1 *= scale;
				getline(ifs,str); iss.str(str);
				iss >> co.a2(0) >> co.a2(1) >> co.a2(2);
				iss.clear(); co.a2 *= scale;
				getline(ifs,str); iss.str(str);
				iss >> co.a3(0) >> co.a3(1) >> co.a3(2);
				iss.clear(); co.a3 *= scale;
				for ( int i=0; i<3; i++ ) co.latmat(i,0) = co.a1(i);
				for ( int i=0; i<3; i++ ) co.latmat(i,1) = co.a2(i);
				for ( int i=0; i<3; i++ ) co.latmat(i,2) = co.a3(i);
				// type ands item
				for ( int i=0; i<2; i++ ) getline(ifs,str);
			}
			// axis line
			getline(ifs,str);
			// configuration
			for ( int i=0; i<co.nions; i++ ) {
				getline(ifs,str); iss.str(str);
				for ( int j=0; j<3; j++ ) iss >> cfg(j);
				co.config[i] = co.latmat*cfg;
				iss.clear();
			}
			MPI_Barrier(MPI_COMM_WORLD);
			// update mode
			if ( mode["gofr"] ) {
				co_gofr.set(co);
				co_gofr.update(rank,nmpi);
			}
			MPI_Barrier(MPI_COMM_WORLD);
//			not using so far
//			if ( mode["sofq"] ) {
//				co_sofq.set(co);
//				MPI_Barrier(MPI_COMM_WORLD);
//			}
			if ( mode["dfc"] && rank == 0 ) {
				co_dfc.set(co);
				co_dfc.update();
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if ( mode["covoro"] && !mode["instant"] ) {
				co_covoro.set(co);
				co_covoro.get();
				for ( int a=0; a<co.ntypes; a++ ) {
				for ( int b=0; b<co.ntypes; b++ ) {
					co_covoro.pget(a,b);
				}
				}
				for ( int spa=0; spa<co.nstypes; spa++ ) {
				for ( int spb=0; spb<co.nstypes; spb++ ) {
					co_covoro.spget(spa,spb);
				}
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
//			not using so far
//			if ( mode["lboo"] && rank == 0 ) {
//				co_lboo.set(co,co_covoro);
//				MPI_Barrier(MPI_COMM_WORLD);
//			}
			if ( mode["bacf"] ) {
				co_bacf.set(co,co_covoro);
				co_bacf.update(rank,nmpi);
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if ( mode["pofqw"] && rank == 0 ) {
				co_pofqw.set(co,co_covoro);
				co_pofqw.update();
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if ( mode["bofree"] && ( (co.itr-co.f0) % (co.df*co.dump) == 0 || co.itr==co.f1 ) ) {
				co_bofree.set(co,co_covoro);
				co_bofree.update(rank,nmpi);
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if ( mode["csform"] ) {
				co_cv_coqhull.set(co,co_covoro);
				co_csform.get();
				co_csform.pget();
				co_csform.spget();
				co_csform.update(rank,nmpi);
			}
			MPI_Barrier(MPI_COMM_WORLD);
			// write out
			if ( ( (co.itr-co.f0) % (co.df*co.dump) == 0 || co.itr==co.f1 ) && rank == 0 ) {

				// write out
				if ( mode["gofr"] ) co_gofr.write("gofr");
				if ( mode["sofq"] ) {
					co_sofq.set(co);
					co_sofq.write("sofq",&co_gofr);
				}
				if ( mode["dfc"] ) co_dfc.write("dfc");
				if ( mode["lboo"] ) {
					if ( mode["covoro"] && mode["instant"] ) {
						co_covoro.set(co);
						co_covoro.get();
						for ( int a=0; a<co.ntypes; a++ ) {
						for ( int b=0; b<co.ntypes; b++ ) {
							co_covoro.pget(a,b);
						}
						}
						for ( int spa=0; spa<co.nstypes; spa++ ) {
						for ( int spb=0; spb<co.nstypes; spb++ ) {
							co_covoro.spget(spa,spb);
						}
						}
					}
					co_lboo.set(co,co_covoro);
					co_lboo.write("lboo");
				}
				if ( mode["bacf"] ) co_bacf.write("bacf");
				if ( mode["pofqw"] ) co_pofqw.write("pofqw");
				if ( mode["bofree"] ) co_bofree.write("bofree");
				if ( mode["povoro"] ) {
					if ( mode["covoro"] && mode["instant"] ) {
						co_covoro.set(co);
						co_covoro.get();
						for ( int a=0; a<co.ntypes; a++ ) {
						for ( int b=0; b<co.ntypes; b++ ) {
							co_covoro.pget(a,b);
						}
						}
						for ( int spa=0; spa<co.nstypes; spa++ ) {
						for ( int spb=0; spb<co.nstypes; spb++ ) {
							co_covoro.spget(spa,spb);
						}
						}
					}
					co_povoro.set(co,co_covoro);
					co_povoro.write("povoro");
					for ( int a=0; a<co.ntypes; a++ ) {
					for ( int b=0; b<co.ntypes; b++ ) {
						co_povoro.pwrite("povoro",a,b);
					}
					}
				}
				if ( mode["csform"] ) {
					co_csform.write("csform");
				}

				// dump progress
				time = MPI_Wtime();
				cout << co.itr << " ( " << (double)(co.itr-co.f0+1)/(double)(co.f1-co.f0+1)*100.0 << " )";
				cout << " : " << time-prog << " ( " << time-beg << " )" << endl;
				prog = time;
			}
			MPI_Barrier(MPI_COMM_WORLD);
		} else {
			// skip
			if ( co.latdyn ) for ( int i=0; i<7; i++ ) getline(ifs,str);
			getline(ifs,str);	// axis
			for ( int i=0; i<co.nions; i++ ) getline(ifs,str);
		}

	}
	MPI_Barrier(MPI_COMM_WORLD);
	ifs.close();
	MPI_Finalize();

	// free
	mode.clear();

	return 0;
}

