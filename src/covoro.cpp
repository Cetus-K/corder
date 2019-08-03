#include "corder.h"

// destructor
covoro::~covoro() {
	vector < vec > ().swap(cell);
	vector < vec > ().swap(config);
	vector < vector < vec > > ().swap(nv);
	vector < vector < vector < vec > > > ().swap(eg);
	vector < vector < int > > ().swap(nnid);
	vector < vector < double > > ().swap(surf);
	vector < vector < vector < vec > > > ().swap(nvab);
	vector < vector < vector < vector < vec > > > > ().swap(egab);
	vector < vector < vector < int > > > ().swap(nnidab);
	vector < vector < vector < double > > > ().swap(surfab);
	vector < vector < vector < vec > > > ().swap(nvsab);
	vector < vector < vector < vector < vec > > > > ().swap(egsab);
	vector < vector < vector < int > > > ().swap(nnidsab);
	vector < vector < vector < double > > > ().swap(surfsab);
}

void covoro::init ( void ) {

	// memory allocation
	cell.resize(3);
	for ( int i=0; i<3; i++ ) cell[i].resize(3);
	// total information
	config.resize(co->nions);
	for ( int i=0; i<(co->nions); i++ ) config[i].resize(3);
	nv.resize(co->nions);
	eg.resize(co->nions);
	nnid.resize(co->nions);
	surf.resize(co->nions);
	// partial information
	nvab.resize(co->ntypes*co->ntypes);
	egab.resize(co->ntypes*co->ntypes);
	nnidab.resize(co->ntypes*co->ntypes);
	surfab.resize(co->ntypes*co->ntypes);
	for ( int iab=0; iab<(co->ntypes*co->ntypes); iab++ ) {
		nvab[iab].resize(co->nions);
		egab[iab].resize(co->nions);
		nnidab[iab].resize(co->nions);
		surfab[iab].resize(co->nions);
	}
	// specofoed partial information
	nvsab.resize(co->nstypes*co->nstypes);
	egsab.resize(co->nstypes*co->nstypes);
	nnidsab.resize(co->nstypes*co->nstypes);
	surfsab.resize(co->nstypes*co->nstypes);
	for ( int isab=0; isab<(co->nstypes*co->nstypes); isab++ ) {
		nvsab[isab].resize(co->nions);
		egsab[isab].resize(co->nions);
		nnidsab[isab].resize(co->nions);
		surfsab[isab].resize(co->nions);
	}

}

// set corder copy constructor
// and match lattice vector format to voro++
// a1,a2,a3 to the form of
//	a1 ,  0,  0
//	a2x,a2y,  0
//	a3x,a3y,a3z
void covoro::set ( corder &cdr ) {
	co = &cdr;

	// temporary variable
	double ya,za,xb;
	vec v(3);
	mat x(3,3),y(3,3),z(3,3);

	// trasform matrix z for a1 to (a1x,0,a1z)
	v = co->a1;
	za = -atan2(v(1),v(0));
	z(0,0) = cos(za); z(0,1) = -sin(za); z(0,2) = 0;
	z(1,0) = sin(za); z(1,1) = cos(za); z(1,2) = 0;
	z(2,0) = 0; z(2,1) = 0; z(2,2) = 1.0;

	// transform matrix y for z*a1 to (a1,0,0)
	v = z*co->a1;
	ya = atan2(v(2),v(0));
	y(0,0) = cos(ya); y(0,1) = 0; y(0,2) = sin(ya);
	y(1,0) = 0; y(1,1) = 1.0; y(1,2) = 0;
	y(2,0) = -sin(ya); y(2,1) = 0; y(2,2) = cos(ya);

	// trasform matrix x for y*z*a2 to (a2x,a2y,0)
	v = (y*z)*co->a2;
	xb = -atan2(v(2),v(1));
	x(0,0) = 1.0; x(0,1) = 0; x(0,2) = 0;
	x(1,0) = 0; x(1,1) = cos(xb); x(1,2) = -sin(xb);
	x(2,0) = 0; x(2,1) = sin(xb); x(2,2) = cos(xb);

	// cell vector which is matched format for voro++
	// transform matrix which match atomic configurations format for voro++
	tfmat = x*y*z;
	cell[0] = tfmat*co->a1;
	cell[1] = tfmat*co->a2;
	cell[2] = tfmat*co->a3;

}

// get informations that voro++ calculated
// get values: surface area, normal vector and near-neighbor ids
// total information
void covoro::get ( void ) {
	int gx=3,gy=3,gz=3,mem=16;		// input
	double ax,axy,ay,axz,ayz,az;
	vec p(3);
	int index;				// looping
	double nabs,*pp;
	vector < int > rmindex;
	vector < double > normal,edge;

	// set up the geometry of constructor
	ax = cell[0](0);
	axy = cell[1](0); ay = cell[1](1);
	axz = cell[2](0); ayz = cell[2](1); az = cell[2](2);
	container_periodic con(ax,axy,ay,axz,ayz,az,gx,gy,gz,mem);

	// put particles
	for ( int i=0; i<(co->nions); i++ ) {
		p = tfmat*co->config[i];
		con.put(i,p(0),p(1),p(2));
	}

	// set looping over all particles
	c_loop_all_periodic vl(con);
	voronoicell_neighbor c;

	// start looping
	if ( vl.start() ) {
		do {
		if ( con.compute_cell(c,vl) ) {

			// this index and the position
			index = con.id[vl.ijk][vl.q];
			pp = con.p[vl.ijk]+con.ps*vl.q;
			config[index](0) = *pp;
			config[index](1) = pp[1];
			config[index](2) = pp[2];

			// pre-saving voro++ information
			c.neighbors(nnid[index]);
			c.face_areas(surf[index]);
			c.normals_with_edges(normal,edge);

			// eliminate illegal normal vector = (0,0,0)
			// counting
			for ( int ij=0; ij<(int)nnid[index].size(); ij++ ) {
				nabs = 0;
				for ( int j=0; j<3; j++ ) nabs += pow(normal[ij*3+j],2.0);
				if ( nabs < under_capacity ) rmindex.push_back(ij);
			}
			// eliminating
			for ( int i=0; i<(int)rmindex.size(); i++ ) {
				nnid[index].erase(nnid[index].begin()+rmindex[i]-i);
				surf[index].erase(surf[index].begin()+rmindex[i]-i);
				normal.erase(normal.begin()+(rmindex[i]-i)*3,normal.begin()+(rmindex[i]-i)*3+3);
				edge.erase(edge.begin()+(rmindex[i]-i)*9,edge.begin()+(rmindex[i]-i)*9+9);
			}
			rmindex.clear();

			// save normal vector to nv[index]
			// save edge and its cross point to
			// 	eg[index][0], eg[index][1] and eg[index][2]
			// memory allocation; jag array allocated
			// due to eliminate illegal normal vectors
			if ( !nv[index].empty() ) vector < vec > ().swap(nv[index]);
			if ( !eg[index].empty() ) {
				vector < vector < vec > > ().swap(eg[index]);
			}
			nv[index].resize(nnid[index].size());
			eg[index].resize(nnid[index].size());
			for ( int ij=0; ij<(int)nnid[index].size(); ij++ ) {
				nv[index][ij].resize(3);
				eg[index][ij].resize(3);
				for ( int j=0; j<3; j++ ) eg[index][ij][j].resize(3);
				for ( int j=0; j<3; j++ ) {
					nv[index][ij](j) = normal[ij*3+j];
					eg[index][ij][0](j) = edge[ij*9+j];
					eg[index][ij][1](j) = edge[ij*9+3+j];
					eg[index][ij][2](j) = pp[j]+0.5*edge[ij*9+6+j];
				}
			}

		}
		} while ( vl.inc() );
	}

	// check
//	c.extract_edges_neighbors();
//	exit(0);

	// free
	vector < int > ().swap(rmindex);
	vector < double > ().swap(normal);
	vector < double > ().swap(edge);

}

// partial information
void covoro::pget ( int a, int b ) {
	int gx=3,gy=3,gz=3,mem=16;		// input
	double ax,axy,ay,axz,ayz,az;
	vec p(3);
	int index;				// looping
	int a_start=co->mystart(a),b_start=co->mystart(b);
	int iab=a*co->ntypes+b;
	double nabs,*pp;
	vector < int > rmindex;
	vector < double > normal,edge;

	// set up the geometry of constructor
	ax = cell[0](0);
	axy = cell[1](0); ay = cell[1](1);
	axz = cell[2](0); ayz = cell[2](1); az = cell[2](2);
	container_periodic con(ax,axy,ay,axz,ayz,az,gx,gy,gz,mem);

	// put particles
	// item[a] ions
	for ( int i=a_start; i<(a_start+co->item[a]); i++ ) {
		p = tfmat*co->config[i];
		con.put(i,p(0),p(1),p(2));
	}
	// item[a]+item[b] ions
	if ( a != b ) {
		for ( int i=b_start; i<(b_start+co->item[b]); i++ ) {
			p = tfmat*co->config[i];
			con.put(i,p(0),p(1),p(2));
		}
	}

	// set looping over all particles
	c_loop_all_periodic vl(con);
	voronoicell_neighbor c;

	// start looping
	if ( vl.start() ) {
		do {
		if ( con.compute_cell(c,vl) ) {

			// this index and the position
			index = con.id[vl.ijk][vl.q];
			pp = con.p[vl.ijk]+con.ps*vl.q;

			// pre-saving voro++ information
			c.neighbors(nnidab[iab][index]);
			c.face_areas(surfab[iab][index]);
			c.normals_with_edges(normal,edge);

			// eliminate illegal normal vector = (0,0,0)
			// counting
			for ( int ij=0; ij<(int)nnidab[iab][index].size(); ij++ ) {
				nabs = 0;
				for ( int j=0; j<3; j++ ) nabs += pow(normal[ij*3+j],2.0);
				if ( nabs < under_capacity ) rmindex.push_back(ij);
			}
			// eliminating
			for ( int i=0; i<(int)rmindex.size(); i++ ) {
				nnidab[iab][index].erase(nnidab[iab][index].begin()+rmindex[i]-i);
				surfab[iab][index].erase(surfab[iab][index].begin()+rmindex[i]-i);
				normal.erase(normal.begin()+(rmindex[i]-i)*3,normal.begin()+(rmindex[i]-i)*3+3);
				edge.erase(edge.begin()+(rmindex[i]-i)*9,edge.begin()+(rmindex[i]-i)*9+9);
			}
			rmindex.clear();

			// save normal vector to nv[index]
			// save edge and its cross point to
			// 	eg[index][0], eg[index][1] and eg[index][2]
			// memory allocation; jag array allocated
			// due to eliminate illegal normal vectors
			if ( !nvab[iab][index].empty() ) vector < vec > ().swap(nvab[iab][index]);
			if ( !egab[iab][index].empty() ) {
				vector < vector < vec > > ().swap(egab[iab][index]);
			}
			nvab[iab][index].resize(nnidab[iab][index].size());
			egab[iab][index].resize(nnidab[iab][index].size());
			for ( int ij=0; ij<(int)nnidab[iab][index].size(); ij++ ) {
				nvab[iab][index][ij].resize(3);
				egab[iab][index][ij].resize(3);
				for ( int j=0; j<3; j++ ) egab[iab][index][ij][j].resize(3);
				for ( int j=0; j<3; j++ ) {
					nvab[iab][index][ij](j) = normal[ij*3+j];
					egab[iab][index][ij][0](j) = edge[ij*9+j];
					egab[iab][index][ij][1](j) = edge[ij*9+3+j];
					egab[iab][index][ij][2](j) = pp[j]+0.5*edge[ij*9+6+j];
				}
			}

		}
		} while ( vl.inc() );
	}

	// free
	vector < int > ().swap(rmindex);
	vector < double > ().swap(normal);
	vector < double > ().swap(edge);

}

// specified partial information
void covoro::spget ( int spa, int spb ) {
	int gx=3,gy=3,gz=3,mem=16;		// input
	double ax,axy,ay,axz,ayz,az;
	vec p(3);
	int index;				// looping
	int a_start,b_start;
	int isab=spa*co->nstypes+spb;
	double nabs,*pp;
	map < int,bool > sbflag;
	vector < int > rmindex;
	vector < double > normal,edge;

	// set up the geometry of constructor
	ax = cell[0](0);
	axy = cell[1](0); ay = cell[1](1);
	axz = cell[2](0); ayz = cell[2](1); az = cell[2](2);
	container_periodic con(ax,axy,ay,axz,ayz,az,gx,gy,gz,mem);

	// set sbflag; whether put b in sb
	for ( int b=0; b<(co->ntypes); b++ ) sbflag[b] = true;

	// put particles
	// sitem[sa] ions
	for ( int ia=0; ia<(int)co->stypeindex[spa].size(); ia++ ) {
		a_start=co->mystart(co->stypeindex[spa][ia]);
		for ( int i=a_start; i<(a_start+co->item[co->stypeindex[spa][ia]]); i++ ) {
			p = tfmat*co->config[i];
			con.put(i,p(0),p(1),p(2));
		}
		sbflag[co->stypeindex[spa][ia]] = false;
	}
	// sitem[sa]+sitem[sb] ions
	for ( int ib=0; ib<(int)co->stypeindex[spb].size(); ib++ ) {
	if ( sbflag[co->stypeindex[spb][ib]] ) {
		b_start=co->mystart(co->stypeindex[spb][ib]);
		for ( int i=b_start; i<(b_start+co->item[co->stypeindex[spb][ib]]); i++ ) {
			p = tfmat*co->config[i];
			con.put(i,p(0),p(1),p(2));
		}
	}
	}
	sbflag.clear();

	// set looping over all particles
	c_loop_all_periodic vl(con);
	voronoicell_neighbor c;

	// start looping
	if ( vl.start() ) {
		do {
		if ( con.compute_cell(c,vl) ) {

			// this index and the position
			index = con.id[vl.ijk][vl.q];
			pp = con.p[vl.ijk]+con.ps*vl.q;

			// pre-saving voro++ information
			c.neighbors(nnidsab[isab][index]);
			c.face_areas(surfsab[isab][index]);
			c.normals_with_edges(normal,edge);

			// eliminate illegal normal vector = (0,0,0)
			// counting
			for ( int ij=0; ij<(int)nnidsab[isab][index].size(); ij++ ) {
				nabs = 0;
				for ( int j=0; j<3; j++ ) nabs += pow(normal[ij*3+j],2.0);
				if ( nabs < under_capacity ) rmindex.push_back(ij);
			}
			// eliminating
			for ( int i=0; i<(int)rmindex.size(); i++ ) {
				nnidsab[isab][index].erase(nnidsab[isab][index].begin()+rmindex[i]-i);
				surfsab[isab][index].erase(surfsab[isab][index].begin()+rmindex[i]-i);
				normal.erase(normal.begin()+(rmindex[i]-i)*3,normal.begin()+(rmindex[i]-i)*3+3);
				edge.erase(edge.begin()+(rmindex[i]-i)*9,edge.begin()+(rmindex[i]-i)*9+9);
			}
			rmindex.clear();

			// save normal vector to nv[index]
			// save edge and its cross point to
			// 	eg[index][0], eg[index][1] and eg[index][2]
			// memory allocation; jag array allocated
			// due to eliminate illegal normal vectors
			if ( !nvsab[isab][index].empty() ) vector < vec > ().swap(nvsab[isab][index]);
			if ( !egsab[isab][index].empty() ) {
				vector < vector < vec > > ().swap(egsab[isab][index]);
			}
			nvsab[isab][index].resize(nnidsab[isab][index].size());
			egsab[isab][index].resize(nnidsab[isab][index].size());
			for ( int ij=0; ij<(int)nnidsab[isab][index].size(); ij++ ) {
				nvsab[isab][index][ij].resize(3);
				egsab[isab][index][ij].resize(3);
				for ( int j=0; j<3; j++ ) egsab[isab][index][ij][j].resize(3);
				for ( int j=0; j<3; j++ ) {
					nvsab[isab][index][ij](j) = normal[ij*3+j];
					egsab[isab][index][ij][0](j) = edge[ij*9+j];
					egsab[isab][index][ij][1](j) = edge[ij*9+3+j];
					egsab[isab][index][ij][2](j) = pp[j]+0.5*edge[ij*9+6+j];
				}
			}

		}
		} while ( vl.inc() );
	}

	// free
	vector < int > ().swap(rmindex);
	vector < double > ().swap(normal);
	vector < double > ().swap(edge);

}
