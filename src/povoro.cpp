#include "corder.h"

void povoro::write ( string filename ) {
	int gx=3,gy=3,gz=3,mem=16;
	double ax,axy,ay,axz,ayz,az;
	vec p(3);

	// set up the geometry of constructor
	ax = cv->cell[0](0);
	axy = cv->cell[1](0); ay = cv->cell[1](1);
	axz = cv->cell[2](0); ayz = cv->cell[2](1); az = cv->cell[2](2);
	container_periodic con(ax,axy,ay,axz,ayz,az,gx,gy,gz,mem);

	// put particles
	for ( int i=0; i<(co->nions); i++ ) {
		p = cv->config[i];
		con.put(i,p(0),p(1),p(2));
	}

	// dump voronoi cell trajectory as pov-ray format
	con.draw_particles_povs((filename+"-particle.pov").c_str(),co->itr,co->f0,co->ntypes,co->item,co->f1-co->f0+1);
	con.draw_cells_povs((filename+"-cell.pov").c_str(),co->itr,co->f0,co->f1-co->f0+1);

}

void povoro::pwrite ( string filename, int a, int b ) {
	int gx=3,gy=3,gz=3,mem=16;
	double ax,axy,ay,axz,ayz,az;
	int a_start=co->mystart(a),b_start=co->mystart(b);
	vec p(3);
	string pair=co->type[a]+co->type[b];

	// set up the geometry of constructor
	ax = cv->cell[0](0);
	axy = cv->cell[1](0); ay = cv->cell[1](1);
	axz = cv->cell[2](0); ayz = cv->cell[2](1); az = cv->cell[2](2);
	container_periodic con(ax,axy,ay,axz,ayz,az,gx,gy,gz,mem);

	// put particles
	// item[a] ions
	for ( int i=a_start; i<(a_start+co->item[a]); i++ ) {
		p = cv->config[i];
		con.put(i,p(0),p(1),p(2));
	}
	// item[a]+item[b] ions
	if ( a != b ) {
		for ( int i=b_start; i<(b_start+co->item[b]); i++ ) {
			p = cv->config[i];
			con.put(i,p(0),p(1),p(2));
		}
	}

	// dump voronoi cell trajectory as pov-ray format
	con.draw_particles_povs((filename+"-"+pair+"-particle.pov").c_str(),co->itr,co->f0,co->ntypes,co->item,co->f1-co->f0+1);
	con.draw_cells_povs((filename+"-"+pair+"-cell.pov").c_str(),co->itr,co->f0,co->f1-co->f0+1);

}

void povoro::spwrite ( string filename, int spa, int spb ) {
	int gx=3,gy=3,gz=3,mem=16;
	double ax,axy,ay,axz,ayz,az;
	int a_start,b_start;
	int isab=spa*co->nstypes+spb;
	map < int,bool > sbflag;
	vec p(3);
	string pair="";

	// set name
	for ( int ia=0; ia<(int)co->stypeindex[spa].size(); ia++ )
		pair += co->type[co->stypeindex[spa][ia]];
	for ( int ib=0; ib<(int)co->stypeindex[spb].size(); ib++ )
		pair += co->type[co->stypeindex[spb][ib]];

	// set up the geometry of constructor
	ax = cv->cell[0](0);
	axy = cv->cell[1](0); ay = cv->cell[1](1);
	axz = cv->cell[2](0); ayz = cv->cell[2](1); az = cv->cell[2](2);
	container_periodic con(ax,axy,ay,axz,ayz,az,gx,gy,gz,mem);

	// set sbflag; whether put b in sb
	for ( int b=0; b<(co->ntypes); b++ ) sbflag[b] = true;

	// put particles
	// sitem[sa] ions
	for ( int ia=0; ia<(int)co->stypeindex[spa].size(); ia++ ) {
		a_start=co->mystart(co->stypeindex[spa][ia]);
		for ( int i=a_start; i<(a_start+co->item[co->stypeindex[spa][ia]]); i++ ) {
			p = cv->config[i];
			con.put(i,p(0),p(1),p(2));
		}
		sbflag[co->stypeindex[spa][ia]] = false;
	}
	// sitem[sa]+sitem[sb] ions
	for ( int ib=0; ib<(int)co->stypeindex[spb].size(); ib++ ) {
	if ( sbflag[co->stypeindex[spb][ib]] ) {
		b_start=co->mystart(co->stypeindex[spb][ib]);
		for ( int i=b_start; i<(b_start+co->item[co->stypeindex[spb][ib]]); i++ ) {
			p = cv->config[i];
			con.put(i,p(0),p(1),p(2));
		}
	}
	}
	sbflag.clear();

	// dump voronoi cell trajectory as pov-ray format
	con.draw_particles_povs((filename+"-"+pair+"-particle.pov").c_str(),co->itr,co->f0,co->nstypes,co->sitem,co->f1-co->f0+1);
	con.draw_cells_povs((filename+"-"+pair+"-cell.pov").c_str(),co->itr,co->f0,co->f1-co->f0+1);

}
