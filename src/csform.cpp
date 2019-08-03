#include "corder.h"

// destructor
csform::~csform() {
	free1d < int > (ddim);
	free2d < int > (ddimab);
	free2d < int > (ddimsab);
	free1d < double > (pb);
	free1d < double > (pf);
	free1d < double > (pe);
	free1d < double > (pv);
	free2d < double > (pbab);
	free2d < double > (pfab);
	free2d < double > (peab);
	free2d < double > (pvab);
	free2d < double > (pbsab);
	free2d < double > (pfsab);
	free2d < double > (pesab);
	free2d < double > (pvsab);
}

void csform::init ( void ) {

	// temporary variable
	ncs = (int)((co->cs1-co->cs0)/co->dcs);

	// memory allocation
	ddim = alloc1d < int > (co->nions);
	ddimab = alloc2d < int > (co->ntypes*co->ntypes,co->nions);
	ddimsab = alloc2d < int > (co->nstypes*co->nstypes,co->nions);
	pb = alloc1d < double > (ncs);
	pf = alloc1d < double > (ncs);
	pe = alloc1d < double > (ncs);
	pv = alloc1d < double > (ncs);
	pbab = alloc2d < double > (co->ntypes*co->ntypes,ncs);
	pfab = alloc2d < double > (co->ntypes*co->ntypes,ncs);
	peab = alloc2d < double > (co->ntypes*co->ntypes,ncs);
	pvab = alloc2d < double > (co->ntypes*co->ntypes,ncs);
	pbsab = alloc2d < double > (co->nstypes*co->nstypes,ncs);
	pfsab = alloc2d < double > (co->nstypes*co->nstypes,ncs);
	pesab = alloc2d < double > (co->nstypes*co->nstypes,ncs);
	pvsab = alloc2d < double > (co->nstypes*co->nstypes,ncs);

}

// flag for calc qhull
// false: all vertices are in a coplanar
//        or 2 < number of vertices
//  true: otherwise
bool csform::qhflag ( vector < vec > points ) {
	// temporary variable
	bool flag=false;
	double y1,z1,x2;
	vec v(3);
	mat x(3,3),y(3,3),z(3,3),tfmat(3,3);

	if ( 2 < (int)points.size() ) {
		// trasform matrix z for r1-r0 to (rx,0,rz)
		v = points[1]-points[0];
		z1 = -atan2(v(1),v(0));
		z(0,0) = cos(z1); z(0,1) = -sin(z1); z(0,2) = 0;
		z(1,0) = sin(z1); z(1,1) = cos(z1); z(1,2) = 0;
		z(2,0) = 0; z(2,1) = 0; z(2,2) = 1.0;
		// transform matrix y for z*(r1-r0) to (r,0,0)
		v = z*(points[1]-points[0]);
		y1 = atan2(v(2),v(0));
		y(0,0) = cos(y1); y(0,1) = 0; y(0,2) = sin(y1);
		y(1,0) = 0; y(1,1) = 1.0; y(1,2) = 0;
		y(2,0) = -sin(y1); y(2,1) = 0; y(2,2) = cos(y1);
		// trasform matrix x for y*z*(r2-r0) to (rx,ry,0)
		v = (y*z)*(points[2]-points[0]);
		x2 = -atan2(v(2),v(1));
		x(0,0) = 1.0; x(0,1) = 0; x(0,2) = 0;
		x(1,0) = 0; x(1,1) = cos(x2); x(1,2) = -sin(x2);
		x(2,0) = 0; x(2,1) = sin(x2); x(2,2) = cos(x2);
		// transform matrix
		tfmat = x*y*z;

		// transform ri -> x*y*z*(ri-r0)
		for ( int i=0; i<(int)points.size(); i++ ) {
			v = tfmat*(points[i]-points[0]);
			if ( under_capacity < pow(v(2),2.0) ) {
				flag = true; break;
			}
		}
	}

	return flag;

}

// calc center position of 2d-polygon
// assuming: all points are in coplanar
//           , 2 < (int)points.size()
//           and the face is convex
// return: vertex's indices cyclically on the face
//         center position of the face
void csform::calc_flat ( vector < vec > points, vector < int > *list, vec *g ) {
	// temporary variable
	int maxindex;
	double y1,z1,x2;
	double cosw,cosw0=0.0,si,stot=0.0;
	vector < double > phi;
	vector < vec > pts(points);
	vec v(3),v1(3),v2(3);
	mat x(3,3),y(3,3),z(3,3),tfmat(3,3);

	// find cyclical index alignment
	// search minimum angular from r1-r0 for starting
	(*list).push_back(1);
	v1 = points[1]-points[0];
	for ( int i=2; i<(int)points.size(); i++ ) {
		v2 = points[i]-points[0];
		cosw = 1.0/(v1.norm()*v2.norm())*v1*v2;
		if ( cosw0 < cosw ) {
			cosw0 = cosw; maxindex = i;
		}
	}
	// start from r2-r1
	if ( (int)points.size() == 3 ) {
		(*list).push_back(2);
		(*list).push_back(0);
	} else {
		int k=0;
		while ( maxindex != (int)(*list)[0] ) {
			cosw0 = 0.0;
			(*list).push_back(maxindex);
			v1 = points[(*list)[k+1]]-points[(*list)[k]];
			for ( int j=0; j<(int)points.size(); j++ ) {
			if ( j != (*list)[k] and j != (*list)[k+1] ) {
				v2 = points[j]-points[(*list)[k]];
				cosw = 1.0/(v1.norm()*v2.norm())*v1*v2;
				if ( cosw0 < cosw ) {
					cosw0 = cosw; maxindex = j;
				}
			}
			}
			k ++;
		}
	}

	// maping into xy-plane
	// trasform matrix z for r1-r0 to (rx,0,rz)
	v = points[(*list)[1]]-points[(*list)[0]];
	z1 = -atan2(v(1),v(0));
	z(0,0) = cos(z1); z(0,1) = -sin(z1); z(0,2) = 0;
	z(1,0) = sin(z1); z(1,1) = cos(z1); z(1,2) = 0;
	z(2,0) = 0; z(2,1) = 0; z(2,2) = 1.0;
	// transform matrix y for z*(r1-r0) to (r,0,0)
	v = z*(points[(*list)[1]]-points[(*list)[0]]);
	y1 = atan2(v(2),v(0));
	y(0,0) = cos(y1); y(0,1) = 0; y(0,2) = sin(y1);
	y(1,0) = 0; y(1,1) = 1.0; y(1,2) = 0;
	y(2,0) = -sin(y1); y(2,1) = 0; y(2,2) = cos(y1);
	// trasform matrix x for y*z*(r2-r0) to (rx,ry,0)
	v = (y*z)*(points[(*list)[2]]-points[(*list)[0]]);
	x2 = -atan2(v(2),v(1));
	x(0,0) = 1.0; x(0,1) = 0; x(0,2) = 0;
	x(1,0) = 0; x(1,1) = cos(x2); x(1,2) = -sin(x2);
	x(2,0) = 0; x(2,1) = sin(x2); x(2,2) = cos(x2);
	// transform matrix
	tfmat = x*y*z;

	// calculate center point of the face
	(*g).zero();
	for ( int i=0; i<(int)(*list).size(); i++ )
		pts[i] = tfmat*(points[(*list)[i]]-points[(*list)[0]]);
	for ( int i=1; i<(int)(*list).size()-1; i++ ) {
		v1=pts[i];v2=pts[i+1];
		si = 0.5*(v1.cross(v2))(2);
		stot += si;
		(*g) += si*1.0/3.0*(v1+v2);
	}
	(*g) = 1.0/stot*tfmat.inv3d()*(*g)+points[(*list)[0]];
	vector < vec > ().swap(pts);

}

/* --- check degeneration for convex hull
| ddim=3: all vertices are not degenerated
   | ddim=2(3d->2d): all vertices are in a coplanar
      | ddim=1(2d->1d): all vertices are in a colinear
         | ddim=0(1d->0d): all vertices are in the same point (?)
            ( * this case is occured at (int)points.size() = 1 )
            | ddim=-1: the ion has no neighbor
               ( * this case may not be occured )
----- */
int csform::dim_convex ( vector < vec > points ) {
	int di=3;
	vec v1(3),v2(3),vijk(3);
	bool d0=true,d1=true;
	// check dimension
	// 3d -> 2d ?
	if ( !qhflag(points) and !points.empty() ) {
		di = 2;
		// 2d -> 1d ?
		if ( (int)points.size() < 3 ) {
			di = 1;
			// 1d -> 0d ?
			if ( pow((v2-v1).norm(),2.0) < under_capacity ) di = 0;
		} else {
			v1=points[0];v2=points[1];
			// check all atoms are in a colinear
			for ( int ij=2; ij<(int)points.size(); ij++ ) {
				vijk = (v2-v1).cross(points[ij]-v1);
				if ( under_capacity < pow(vijk.norm(),2.0) ) {
					d0=false;d1=false; break;
				} else {
					// 1d -> 0d ?
					if ( under_capacity < pow((v1-v2).norm(),2.0) ) d0=false;
				}
				if ( !d1 ) break;
			}
		}
		if ( d1 ) di = 1;
		if ( d0 ) di = 0;
	// no neighbor ?
	} else if ( points.empty() ) di=-1;
	// return
	return di;
}

// get convex hull
void csform::get ( void ) {
	// temporary variables
	double s0,tol=0.05;
	vec pair(3); mat nwe(3,3);
	vec edge(3);
	vector < int > veijk(2);
	

	// set vertices
	for ( int i=0; i<(co->nions); i++ ) {
		// maximum surface area
		s0 = 0.0;
		for ( int ij=0; ij<(int)cv->surf[i].size(); ij++ )
			if ( s0 < cv->surf[i][ij] ) s0 = cv->surf[i][ij];
		// free if vt[i] is remaining
		if ( !cq->vt[i].empty() ) 
			vector < vec > ().swap(cq->vt[i]);
		// loop for all neghbors around i-th ion
		// set neghboring atoms from i-th ion
		for ( int ij=0; ij<(int)cv->nnid[i].size(); ij++ ) {
			// calc r(i)'s pair position of r(nnid[i][ij])
			// from normal vector of voronoi surface
			// which direction is r(i) to r(nnid[i][ij])
			for ( int j=0; j<3; j++ ) {
				nwe(0,j) = cv->eg[i][ij][0](j);
				nwe(1,j) = cv->eg[i][ij][1](j);
				nwe(2,j) = cv->nv[i][ij](j);
			}
			pair(0) = cv->config[i]*cv->eg[i][ij][0];
			pair(1) = cv->config[i]*cv->eg[i][ij][1];
			pair(2) = (2.0*cv->eg[i][ij][2]-cv->config[i])*cv->nv[i][ij];
			pair = nwe.inv3d()*pair;
			// round x to 0.0 for x^2 < under_capacity
			for ( int ijk=0; ijk<3; ijk++ ) {
				if ( pow(pair(ijk),2.0) < under_capacity )
					pair(ijk) = 0.0;
			}
			// set point
			if ( tol < cv->surf[i][ij]/s0 )
				cq->vt[i].push_back(pair);
		}
		// check whether the vertices are in coplanar or colinear
		ddim[i] = dim_convex(cq->vt[i]);
		if ( ddim[i] == 3 ) {
			// calc qhull for one point sets
			cq->get_one(cq->vt[i],i);
		} else {
			// free if nv[i] and vnv[i] are remaining
			if ( !cq->nv[i].empty() ) {
				vector < vec > ().swap(cq->nv[i]);
				vector < vector < int > > ().swap(cq->vnv[i]);
			}
			// free if eg[i] and veg[i] are remaining
			if ( !cq->eg[i].empty() ) {
				vector < vec > ().swap(cq->eg[i]);
				vector < vector < int > > ().swap(cq->veg[i]);
			}
			// dealing cases of few neighbor
			// ddim = -1: dump error
			// ddim = 0: remove vt[i]'s ions
			// ddim = 1: do nothing ( only vertex information )
			// ddim = 2: save edge information
			if ( ddim[i] == 0 ) vector < vec > ().swap(cq->vt[i]);
			if ( ddim[i] == 2 ) {
				vec v1(3),v2(3),g(3);
				vector < int > list;
				// get cyclical vertex list and center position
				calc_flat(cq->vt[i],&list,&g);
				// save edge vector and whose vertices
				for ( int ij=0; ij<(int)list.size()-1; ij++ ) {
					v1 = cq->vt[i][list[ij]];
					v2 = cq->vt[i][list[ij+1]];
					edge = 1.0/(v1-v2).norm()*(v1-v2);
					veijk[0] = list[ij];
					veijk[1] = list[ij+1];
					cq->eg[i].push_back(edge);
					cq->veg[i].push_back(veijk);
				}
				// free
				vector < int > ().swap(list);
			}
		}
	}

	// check
//	double sf;
//	for ( int i=0; i<(co->nions); i++ ) {
//		cout << "---------------------------" << endl;
//		cout << "atom(" << i+1 << ") = " << endl;
//		cv->config[i].disp();
//		sf = 0.0;
//		for ( int ij=0; ij<(int)cq->vt[i].size(); ij++ ) sf += cv->surf[i][ij];
//		cout << "sf(" << i+1 << ") = " << sf << endl;
//		cout << (int)cq->vt[i].size() << " vertices" << endl;
//		cout << cq->eg[i].size() << " edges" << endl;
//		cout << cq->nv[i].size() << " facets" << endl;
//		cout << "---------------------------" << endl;
//		for ( int ij=0; ij<(int)cq->vt[i].size(); ij++ ) {
//			cout << "surf(" << i+1 << "," << ij+1 << ") = ";
//			cout << cv->surf[i][ij] << endl;
//			cout << "nn(" << i+1 << "," << ij+1 << ") = ";
//			cout << "atom(" << (int)cv->nnid[i][ij]+1 << ")" << endl;
//			cq->vt[i][ij].disp();
//		}
//		cout << "---------------------------" << endl;
//	}

}

// get convex hull as species pair
void csform::pget ( void ) {
	// temporary variables
	int iab,a_start;
	double s0,tol=0.05;
	vec pair(3); mat nwe(3,3);
	vec edge(3);
	vector < int > veijk(2);

	for ( int a=0; a<(co->ntypes); a++ ) { a_start = co->mystart(a);
	for ( int b=0; b<(co->ntypes); b++ ) { iab = a*co->ntypes+b;
	for ( int i=a_start; i<(a_start+co->item[a]); i++ ) {
		s0 = 0.0;
		for ( int ij=0; ij<(int)cv->surfab[iab][i].size(); ij++ )
			if ( s0 < cv->surfab[iab][i][ij] ) s0 = cv->surfab[iab][i][ij];
		if ( !cq->vtab[iab][i].empty() ) 
			vector < vec > ().swap(cq->vtab[iab][i]);
		for ( int ij=0; ij<(int)cv->nnidab[iab][i].size(); ij++ ) {
			for ( int j=0; j<3; j++ ) {
				nwe(0,j) = cv->egab[iab][i][ij][0](j);
				nwe(1,j) = cv->egab[iab][i][ij][1](j);
				nwe(2,j) = cv->nvab[iab][i][ij](j);
			}
			pair(0) = cv->config[i]*cv->egab[iab][i][ij][0];
			pair(1) = cv->config[i]*cv->egab[iab][i][ij][1];
			pair(2) = (2.0*cv->egab[iab][i][ij][2]-cv->config[i])*cv->nvab[iab][i][ij];
			pair = nwe.inv3d()*pair;
			for ( int ijk=0; ijk<3; ijk++ ) {
				if ( pow(pair(ijk),2.0) < under_capacity )
					pair(ijk) = 0.0;
			}
			// set point ab-pair only
			// a-center, b-coordinate
			if ( tol < cv->surfab[iab][i][ij]/s0 ) {
			if ( a == b ) cq->vtab[iab][i].push_back(pair);
			else {
				if ( b == co->mytype(cv->nnidab[iab][i][ij]) )
					cq->vtab[iab][i].push_back(pair);
			}
			}
		}

		ddimab[iab][i] = dim_convex(cq->vtab[iab][i]);
		if ( ddimab[iab][i] == 3 ) {
			cq->pget_one(cq->vtab[iab][i],i,a,b);
		} else {
			if ( !cq->nvab[iab][i].empty() ) {
				vector < vec > ().swap(cq->nvab[iab][i]);
				vector < vector < int > > ().swap(cq->vnvab[iab][i]);
			}
			if ( !cq->egab[iab][i].empty() ) {
				vector < vec > ().swap(cq->egab[iab][i]);
				vector < vector < int > > ().swap(cq->vegab[iab][i]);
			}
			if ( ddimab[iab][i] == 0 ) vector < vec > ().swap(cq->vtab[iab][i]);
			if ( ddimab[iab][i] == 2 ) {
				vec v1(3),v2(3),g(3);
				vector < int > list;
				// get cyclical vertex list and center position
				calc_flat(cq->vtab[iab][i],&list,&g);
				// save edge vector and whose vertices
				for ( int ij=0; ij<(int)list.size()-1; ij++ ) {
					v1 = cq->vtab[iab][i][list[ij]];
					v2 = cq->vtab[iab][i][list[ij+1]];
					edge = 1.0/(v1-v2).norm()*(v1-v2);
					veijk[0] = list[ij];
					veijk[1] = list[ij+1];
					cq->egab[iab][i].push_back(edge);
					cq->vegab[iab][i].push_back(veijk);
				}
				vector < int > ().swap(list);
			}
		}

	}
	}
	}

	// check
//	double sf;
//	for ( int a=0; a<(co->ntypes); a++ ) { a_start = co->mystart(a);
//	for ( int b=0; b<(co->ntypes); b++ ) { iab = a*co->ntypes+b;
//	for ( int i=a_start; i<(a_start+co->item[a]); i++ ) {
//		cout << "---------------------------" << endl;
//		cout << "type(" << a+1 << "," << b+1 << "):" << endl;
//		cout << "atom(" << i+1 << ") = " << endl;
//		cv->config[i].disp();
//		sf = 0.0;
//		for ( int ij=0; ij<(int)cq->vtab[iab][i].size(); ij++ ) sf += cv->surfab[iab][i][ij];
//		cout << "sf(" << i+1 << ") = " << sf << endl;
//		cout << (int)cq->vtab[iab][i].size() << " vertices" << endl;
//		cout << cq->egab[iab][i].size() << " edges" << endl;
//		cout << cq->nvab[iab][i].size() << " facets" << endl;
//		cout << "---------------------------" << endl;
//		for ( int ij=0; ij<(int)cq->vtab[iab][i].size(); ij++ ) {
//			cout << "surf(" << i+1 << "," << ij+1 << ") = ";
//			cout << cv->surfab[iab][i][ij] << endl;
//			cout << "nn(" << i+1 << "," << ij+1 << ") = ";
//			cout << "atom(" << (int)cv->nnidab[iab][i][ij]+1 << ")" << endl;
//			cq->vtab[iab][i][ij].disp();
//		}
//		cout << "---------------------------" << endl;
//	}
//	}
//	}

}

// get convex hull as species pair
void csform::spget ( void ) {
	// temporary variables
	int isab,a_start;
	double s0,tol=0.05;
	vec pair(3); mat nwe(3,3);
	vec edge(3);
	vector < int > veijk(2);

	for ( int spa=0; spa<(co->nstypes); spa++ ) {
	for ( int spb=0; spb<(co->nstypes); spb++ ) { isab = spa*co->nstypes+spb;
	for ( int ia=0; ia<(int)co->stypeindex[spa].size(); ia++ ) {
		int a=co->stypeindex[spa][ia];a_start=co->mystart(a);
		for ( int i=a_start; i<(a_start+co->item[a]); i++ ) {
			s0 = 0.0;
			for ( int ij=0; ij<(int)cv->surfsab[isab][i].size(); ij++ )
				if ( s0 < cv->surfsab[isab][i][ij] ) s0 = cv->surfsab[isab][i][ij];
			if ( !cq->vtsab[isab][i].empty() ) 
				vector < vec > ().swap(cq->vtsab[isab][i]);
			for ( int ij=0; ij<(int)cv->nnidsab[isab][i].size(); ij++ ) {
				for ( int j=0; j<3; j++ ) {
					nwe(0,j) = cv->egsab[isab][i][ij][0](j);
					nwe(1,j) = cv->egsab[isab][i][ij][1](j);
					nwe(2,j) = cv->nvsab[isab][i][ij](j);
				}
				pair(0) = cv->config[i]*cv->egsab[isab][i][ij][0];
				pair(1) = cv->config[i]*cv->egsab[isab][i][ij][1];
				pair(2) = (2.0*cv->egsab[isab][i][ij][2]-cv->config[i])*cv->nvsab[isab][i][ij];
				pair = nwe.inv3d()*pair;
				for ( int ijk=0; ijk<3; ijk++ ) {
					if ( pow(pair(ijk),2.0) < under_capacity )
						pair(ijk) = 0.0;
				}
				// set point ab-pair only
				// a-center, b-coordinate
				if ( tol < cv->surfsab[isab][i][ij]/s0 ) {
				if ( spa == spb ) cq->vtsab[isab][i].push_back(pair);
				else {
					int b=co->mytype(cv->nnidsab[isab][i][ij]);
					vector<int>::iterator ib=find(co->stypelist[b].begin(),co->stypelist[b].end(),spb);
					if ( ib != co->stypelist[b].end() )
						cq->vtsab[isab][i].push_back(pair);
				}
				}
			}
	
			ddimsab[isab][i] = dim_convex(cq->vtsab[isab][i]);
			if ( ddimsab[isab][i] == 3 ) {
				cq->spget_one(cq->vtsab[isab][i],i,spa,spb);
			} else {
				if ( !cq->nvsab[isab][i].empty() ) {
					vector < vec > ().swap(cq->nvsab[isab][i]);
					vector < vector < int > > ().swap(cq->vnvsab[isab][i]);
				}
				if ( !cq->egsab[isab][i].empty() ) {
					vector < vec > ().swap(cq->egsab[isab][i]);
					vector < vector < int > > ().swap(cq->vegsab[isab][i]);
				}
				if ( ddimsab[isab][i] == 0 ) vector < vec > ().swap(cq->vtsab[isab][i]);
				if ( ddimsab[isab][i] == 2 ) {
					vec v1(3),v2(3),g(3);
					vector < int > list;
					// get cyclical vertex list and center position
					calc_flat(cq->vtsab[isab][i],&list,&g);
					// save edge vector and whose vertices
					for ( int ij=0; ij<(int)list.size()-1; ij++ ) {
						v1 = cq->vtsab[isab][i][list[ij]];
						v2 = cq->vtsab[isab][i][list[ij+1]];
						edge = 1.0/(v1-v2).norm()*(v1-v2);
						veijk[0] = list[ij];
						veijk[1] = list[ij+1];
						cq->egsab[isab][i].push_back(edge);
						cq->vegsab[isab][i].push_back(veijk);
					}
					vector < int > ().swap(list);
				}
			}
		}
	}
	}
	}

}

// find separation axis
// true: no collision
// false: collision
// --- separation(cq->vt[i],cq->vt[j],gij,eij);
// search the range of the separating axis [tij0:tij1]
// line [ e = eij*t + gij ] and r(i)'s point [ rij ]
// required: (rij-e).eij = 0 => tij = (rij-gij).eij/eij.eij
// search: tij0 = min(tij) and tij1 = max(tij)
bool csform::separation ( vector < vec > vti, vector < vec > vtj, vec gij, vec eij ) {
	bool sep=false,bdr0,bdr1;
	double tc,tij0,tij1,tijk0,tijk1;
	// search range from rij to eij-axis for detect
	tij0=1000.0;tij1=-1000.0;
	for ( int ij=0; ij<(int)vti.size(); ij++ ) {
		tc = 1.0/pow(eij.norm(),2.0)*(vti[ij]-gij)*eij;
		if ( tc < tij0 ) tij0 = tc;
		if ( tij1 < tc ) tij1 = tc;
	}
	// search range from rjk to eij-axis for detect
	tijk0=1000.0;tijk1=-1000.0;
	for ( int jk=0; jk<(int)vtj.size(); jk++ ) {
		tc = 1.0/pow(eij.norm(),2.0)*(vtj[jk]-gij)*eij;
		if ( tc < tijk0 ) tijk0 = tc;
		if ( tijk1 < tc ) tijk1 = tc;
	}
	// check whether the two convex hulls are in collision
	bdr0 = ( pow(tij0-tijk1,4.0) < under_capacity );
	bdr1 = ( pow(tijk0-tij1,4.0) < under_capacity );
	if ( bdr0 or bdr1 or tijk1 < tij0 or tij1 < tijk0 ) sep = true;
	// return
	return sep;
}

// check whether the two x-dimensional convex hulls are in collision
// iab = -1 and isab = -1: total information
// iab > -1 and isab = -1: partial information
// iab = -1 and isab > -1: specified partial information
// ----
// searching the separating axis for all edges ( m+n+mn )
//   I: loop for r(i)'s edges eij
//  II: loop for r(j)'s edges ejk
// III: loop for r(i)'s and r(j)'s edges eij x ejk
bool csform::collision ( int i, int j, int iab, int isab ) {
	double tij0,tij1,tijk0,tijk1;
	bool coll=false;
	vec eijk(3),g(3);
	vec ri0(3),ri1(3),rj0(3),rj1(3);
	// save bounding box
	for ( int ijk=0; ijk<3; ijk++ ) {
		ri0(ijk)=1000.0;ri1(ijk)=-1000.0;
	}
	for ( int jkl=0; jkl<3; jkl++ ) {
		rj0(jkl)=1000.0;rj1(jkl)=-1000.0;
	}
	if ( iab < 0 and isab < 0 ) {

		// --- simple comparing --- //
		// search whole enclosing bounding boxes
		for ( int ij=0; ij<(int)cq->vt[i].size(); ij++ ) {
			for ( int ijk=0; ijk<3; ijk++ ) {
				if ( cq->vt[i][ij](ijk) < ri0(ijk) ) ri0(ijk) = cq->vt[i][ij](ijk);
				if ( cq->vt[i][ij](ijk) > ri1(ijk) ) ri1(ijk) = cq->vt[i][ij](ijk);
			}
		}
		for ( int jk=0; jk<(int)cq->vt[j].size(); jk++ ) {
			for ( int jkl=0; jkl<3; jkl++ ) {
				if ( cq->vt[j][jk](jkl) < rj0(jkl) ) rj0(jkl) = cq->vt[j][jk](jkl);
				if ( cq->vt[j][jk](jkl) > rj1(jkl) ) rj1(jkl) = cq->vt[j][jk](jkl);
			}
		}
		// check whether the two bounding boxes are in collision
		for ( int p=0; p<3; p++ ) {
			if ( rj1(p) < ri0(p) or ri1(p) < rj0(p) ) goto END;
		}

		// --- search separating axis --- //
		// xd vs yd ( 1 < x )
		if ( !cq->eg[i].empty() ) {

			// loop for r(i)'s edges
			for ( int ij=0; ij<(int)cq->eg[i].size(); ij++ ) {
				g = 0.5*(cq->vt[i][cq->veg[i][ij][1]]+cq->vt[i][cq->veg[i][ij][0]]);
				if ( separation(cq->vt[i],cq->vt[j],g,cq->eg[i][ij]) ) goto END;
			}

			// loop for r(j)'s edges
			if ( !cq->eg[j].empty() ) { // xd vs yd ( 1 < x and 1 < y )
				// loop for r(j)'s edges
				for ( int jk=0; jk<(int)cq->eg[j].size(); jk++ ) {
					// search range from rjl to ejk-axis for detect
					g = 0.5*(cq->vt[j][cq->veg[j][jk][1]]+cq->vt[j][cq->veg[j][jk][0]]);
					if ( separation(cq->vt[j],cq->vt[i],g,cq->eg[j][jk]) ) goto END;
				}
			} else if ( !cq->vt[j].empty() ) { // end of xd vs yd ( 1 < x and 1 < y ) and begin xd vs 1d
				g = 0.5*(cq->vt[j][1]+cq->vt[j][0]);
				if ( separation(cq->vt[j],cq->vt[i],g,cq->vt[j][1]-cq->vt[j][0]) ) goto END;
			} // xd vs 1d

			// search more detail if the cross flag is still true
			// loop for r(i)'s edges
			for ( int ij=0; ij<(int)cq->eg[i].size(); ij++ ) {
				g = 0.5*(cq->vt[i][cq->veg[i][ij][1]]+cq->vt[i][cq->veg[i][ij][0]]);
				if ( !cq->eg[j].empty() ) { // xd vs yd ( 1 < x and 1 < y )
					// loop for r(j)'s edges
					for ( int jk=0; jk<(int)cq->eg[j].size(); jk++ ) {
						eijk = cq->eg[i][ij].cross(cq->eg[j][jk]);
						if ( under_capacity < pow(eijk.norm(),4.0) )
							if ( separation(cq->vt[i],cq->vt[j],g,eijk) ) goto END;
					} 
				} else if ( !cq->vt[j].empty() ) { // end of xd vs yd ( 1 < x and 1 < y ) and begin of xd vs 1d
					// loop for r(j)'s edges
					eijk = cq->eg[i][ij].cross(cq->vt[j][1]-cq->vt[j][0]);
					if ( under_capacity < pow(eijk.norm(),4.0) )
						if ( separation(cq->vt[i],cq->vt[j],g,eijk) ) goto END;
				} // xd vs 1d
			}

		// end of xd vs yd ( 1 < x ) and begin of 1d vs yd
		} else if ( !cq->vt[i].empty() ) {

			// for r(i)'s edge
			g = 0.5*(cq->vt[i][1]+cq->vt[i][0]);
			coll = separation(cq->vt[i],cq->vt[j],g,cq->vt[i][1]-cq->vt[i][0]);

			// loop for r(j)'s edges
			if ( !cq->eg[j].empty() ) { // 1d vs yd ( 1 < y )
				for ( int jk=0; jk<(int)cq->eg[j].size(); jk++ ) {
					g = 0.5*(cq->vt[j][cq->veg[j][jk][1]]+cq->vt[j][cq->veg[j][jk][0]]);
					if ( separation(cq->vt[j],cq->vt[i],g,cq->eg[j][jk]) ) goto END;
				}
			} else if ( !cq->vt[j].empty() ) { // end of 1d vs yd ( 1 < y ) and begin 1d vs 1d
				g = 0.5*(cq->vt[j][1]+cq->vt[j][0]);
				if ( separation(cq->vt[j],cq->vt[i],g,cq->vt[j][1]-cq->vt[j][0]) ) goto END;
			} // 1d vs 1d

			// search more detail if the cross flag is still true
			g = 0.5*(cq->vt[i][1]+cq->vt[i][0]);
			if ( !cq->eg[j].empty() ) { // 1d vs yd ( 1 < y )
				// loop for r(j)'s edges
				for ( int jk=0; jk<(int)cq->eg[j].size(); jk++ ) {
					eijk = (cq->vt[i][1]-cq->vt[i][0]).cross(cq->eg[j][jk]);
					if ( under_capacity < pow(eijk.norm(),4.0) )
						if ( separation(cq->vt[i],cq->vt[j],g,eijk) ) goto END;
				}
			} else if ( !cq->vt[j].empty() ) { // end of 1d vs yd ( 1 < y ) and begin of 1d vs 1d
				eijk = (cq->vt[i][1]-cq->vt[i][0]).cross(cq->vt[j][1]-cq->vt[j][0]);
				if ( under_capacity < pow(eijk.norm(),4.0) )
					if ( separation(cq->vt[i],cq->vt[j],g,eijk) ) goto END;
			} // 1d vs 1d

		} // 1d vs yd

	// end of total and begin of partial
	} else if ( -1 < iab and isab < 0 ) {

		// simple comparing
		// search whole enclosing bounding boxes
		for ( int ij=0; ij<(int)cq->vtab[iab][i].size(); ij++ ) {
			for ( int ijk=0; ijk<3; ijk++ ) {
				if ( cq->vtab[iab][i][ij](ijk) < ri0(ijk) ) ri0(ijk) = cq->vtab[iab][i][ij](ijk);
				if ( cq->vtab[iab][i][ij](ijk) > ri1(ijk) ) ri1(ijk) = cq->vtab[iab][i][ij](ijk);
			}
		}
		for ( int jk=0; jk<(int)cq->vtab[iab][j].size(); jk++ ) {
			for ( int jkl=0; jkl<3; jkl++ ) {
				if ( cq->vtab[iab][j][jk](jkl) < rj0(jkl) ) rj0(jkl) = cq->vtab[iab][j][jk](jkl);
				if ( cq->vtab[iab][j][jk](jkl) > rj1(jkl) ) rj1(jkl) = cq->vtab[iab][j][jk](jkl);
			}
		}
		// check whether the two bounding boxes are in collision
		for ( int p=0; p<3; p++ ) {
			if ( rj1(p) < ri0(p) or ri1(p) < rj0(p) ) goto END;
		}

		// search separating axis
		// xd vs yd ( 1 < x )
		if ( !cq->egab[iab][i].empty() ) {

			// loop for r(i)'s edges
			for ( int ij=0; ij<(int)cq->egab[iab][i].size(); ij++ ) {
				g = 0.5*(cq->vtab[iab][i][cq->vegab[iab][i][ij][1]]+cq->vtab[iab][i][cq->vegab[iab][i][ij][0]]);
				if ( separation(cq->vtab[iab][i],cq->vtab[iab][j],g,cq->egab[iab][i][ij]) ) goto END;
			}

			// loop for r(j)'s edges
			if ( !cq->egab[iab][j].empty() ) { // xd vs yd ( 1 < x and 1 < y )
				// loop for r(j)'s edges
				for ( int jk=0; jk<(int)cq->egab[iab][j].size(); jk++ ) {
					// search range from rjl to ejk-axis for detect
					g = 0.5*(cq->vtab[iab][j][cq->vegab[iab][j][jk][1]]+cq->vtab[iab][j][cq->vegab[iab][j][jk][0]]);
					if ( separation(cq->vtab[iab][j],cq->vtab[iab][i],g,cq->egab[iab][j][jk]) ) goto END;
				}
			} else if ( !cq->vtab[iab][j].empty() ) { // end of xd vs yd ( 1 < x and 1 < y ) and begin xd vs 1d
				g = 0.5*(cq->vtab[iab][j][1]+cq->vtab[iab][j][0]);
				if ( separation(cq->vtab[iab][j],cq->vtab[iab][i],g,cq->vtab[iab][j][1]-cq->vtab[iab][j][0]) ) goto END;
			} // xd vs 1d

			// search more detail if the cross flag is still true
			// loop for r(i)'s edges
			for ( int ij=0; ij<(int)cq->egab[iab][i].size(); ij++ ) {
				g = 0.5*(cq->vtab[iab][i][cq->vegab[iab][i][ij][1]]+cq->vtab[iab][i][cq->vegab[iab][i][ij][0]]);
				if ( !cq->egab[iab][j].empty() ) { // xd vs yd ( 1 < x and 1 < y )
					// loop for r(j)'s edges
					for ( int jk=0; jk<(int)cq->egab[iab][j].size(); jk++ ) {
						eijk = cq->egab[iab][i][ij].cross(cq->egab[iab][j][jk]);
						if ( under_capacity < pow(eijk.norm(),4.0) ) {
							if ( separation(cq->vtab[iab][i],cq->vtab[iab][j],g,eijk) ) goto END;
						}
					} 
				} else if ( !cq->vtab[iab][j].empty() ) { // end of xd vs yd ( 1 < x and 1 < y ) and begin of xd vs 1d
					// loop for r(j)'s edges
					eijk = cq->egab[iab][i][ij].cross(cq->vtab[iab][j][1]-cq->vtab[iab][j][0]);
					if ( under_capacity < pow(eijk.norm(),4.0) )
						if ( separation(cq->vtab[iab][i],cq->vtab[iab][j],g,eijk) ) goto END;
				} // xd vs 1d
			}

		// end of xd vs yd ( 1 < x ) and begin of 1d vs yd
		} else if ( !cq->vtab[iab][i].empty() ) {

			// for r(i)'s edge
			g = 0.5*(cq->vtab[iab][i][1]+cq->vtab[iab][i][0]);
			if ( separation(cq->vtab[iab][i],cq->vtab[iab][j],g,cq->vtab[iab][i][1]-cq->vtab[iab][i][0]) ) goto END;

			// loop for r(j)'s edges
			if ( !cq->egab[iab][j].empty() ) { // 1d vs yd ( 1 < y )
				for ( int jk=0; jk<(int)cq->egab[iab][j].size(); jk++ ) {
					g = 0.5*(cq->vtab[iab][j][cq->vegab[iab][j][jk][1]]+cq->vtab[iab][j][cq->vegab[iab][j][jk][0]]);
					if ( separation(cq->vtab[iab][j],cq->vtab[iab][i],g,cq->egab[iab][j][jk]) ) goto END;
				}
			} else if ( !cq->vtab[iab][j].empty() ) { // end of 1d vs yd ( 1 < y ) and begin 1d vs 1d
				g = 0.5*(cq->vtab[iab][j][1]+cq->vtab[iab][j][0]);
				if ( separation(cq->vtab[iab][j],cq->vtab[iab][i],g,cq->vtab[iab][j][1]-cq->vtab[iab][j][0]) ) goto END;
			} // 1d vs 1d

			// search more detail if the cross flag is still true
			g = 0.5*(cq->vtab[iab][i][1]+cq->vtab[iab][i][0]);
			if ( !cq->egab[iab][j].empty() ) { // 1d vs yd ( 1 < y )
				// loop for r(j)'s edges
				for ( int jk=0; jk<(int)cq->egab[iab][j].size(); jk++ ) {
					eijk = (cq->vtab[iab][i][1]-cq->vtab[iab][i][0]).cross(cq->egab[iab][j][jk]);
					if ( under_capacity < pow(eijk.norm(),4.0) )
						if ( separation(cq->vtab[iab][i],cq->vtab[iab][j],g,eijk) ) goto END;
				}
			} else if ( !cq->vtab[iab][j].empty() ) { // end of 1d vs yd ( 1 < y ) and begin of 1d vs 1d
				eijk = (cq->vtab[iab][i][1]-cq->vtab[iab][i][0]).cross(cq->vtab[iab][j][1]-cq->vtab[iab][j][0]);
				if ( under_capacity < pow(eijk.norm(),4.0) )
					if ( separation(cq->vtab[iab][i],cq->vtab[iab][j],g,eijk) ) goto END;
			} // 1d vs 1d

		} // 1d vs yd

	} else if ( iab < 0 and -1 < isab ) { // end of partial and begin of specified partial

		// simple comparing
		// search whole enclosing bounding boxes
		for ( int ij=0; ij<(int)cq->vtsab[isab][i].size(); ij++ ) {
			for ( int ijk=0; ijk<3; ijk++ ) {
				if ( cq->vtsab[isab][i][ij](ijk) < ri0(ijk) ) ri0(ijk) = cq->vtsab[isab][i][ij](ijk);
				if ( cq->vtsab[isab][i][ij](ijk) > ri1(ijk) ) ri1(ijk) = cq->vtsab[isab][i][ij](ijk);
			}
		}
		for ( int jk=0; jk<(int)cq->vtsab[isab][j].size(); jk++ ) {
			for ( int jkl=0; jkl<3; jkl++ ) {
				if ( cq->vtsab[isab][j][jk](jkl) < rj0(jkl) ) rj0(jkl) = cq->vtsab[isab][j][jk](jkl);
				if ( cq->vtsab[isab][j][jk](jkl) > rj1(jkl) ) rj1(jkl) = cq->vtsab[isab][j][jk](jkl);
			}
		}
		// check whether the two bounding boxes are in collision
		for ( int p=0; p<3; p++ ) {
			if ( rj1(p) < ri0(p) or ri1(p) < rj0(p) ) goto END;
		}

		// search separating axis
		// xd vs yd ( 1 < x )
		if ( !cq->egsab[isab][i].empty() ) {

			// loop for r(i)'s edges
			for ( int ij=0; ij<(int)cq->egsab[isab][i].size(); ij++ ) {
				g = 0.5*(cq->vtsab[isab][i][cq->vegsab[isab][i][ij][1]]+cq->vtsab[isab][i][cq->vegsab[isab][i][ij][0]]);
				if ( separation(cq->vtsab[isab][i],cq->vtsab[isab][j],g,cq->egsab[isab][i][ij]) ) goto END;
			}

			// loop for r(j)'s edges
			if ( !cq->egsab[isab][j].empty() ) { // xd vs yd ( 1 < x and 1 < y )
				// loop for r(j)'s edges
				for ( int jk=0; jk<(int)cq->egsab[isab][j].size(); jk++ ) {
					// search range from rjl to ejk-axis for detect
					g = 0.5*(cq->vtsab[isab][j][cq->vegsab[isab][j][jk][1]]+cq->vtsab[isab][j][cq->vegsab[isab][j][jk][0]]);
					if ( separation(cq->vtsab[isab][j],cq->vtsab[isab][i],g,cq->egsab[isab][j][jk]) ) goto END;
				}
			} else if ( !cq->vtsab[isab][j].empty() ) { // end of xd vs yd ( 1 < x and 1 < y ) and begin xd vs 1d
				g = 0.5*(cq->vtsab[isab][j][1]+cq->vtsab[isab][j][0]);
				if ( separation(cq->vtsab[isab][j],cq->vtsab[isab][i],g,cq->vtsab[isab][j][1]-cq->vtsab[isab][j][0]) ) goto END;
			} // xd vs 1d

			// search more detail if the cross flag is still true
			// loop for r(i)'s edges
			for ( int ij=0; ij<(int)cq->egsab[isab][i].size(); ij++ ) {
				g = 0.5*(cq->vtsab[isab][i][cq->vegsab[isab][i][ij][1]]+cq->vtsab[isab][i][cq->vegsab[isab][i][ij][0]]);
				if ( !cq->egsab[isab][j].empty() ) { // xd vs yd ( 1 < x and 1 < y )
					// loop for r(j)'s edges
					for ( int jk=0; jk<(int)cq->egsab[isab][j].size(); jk++ ) {
						eijk = cq->egsab[isab][i][ij].cross(cq->egsab[isab][j][jk]);
						if ( under_capacity < pow(eijk.norm(),4.0) ) {
							if ( separation(cq->vtsab[isab][i],cq->vtsab[isab][j],g,eijk) ) goto END;
						}
					} 
				} else if ( !cq->vtsab[isab][j].empty() ) { // end of xd vs yd ( 1 < x and 1 < y ) and begin of xd vs 1d
					// loop for r(j)'s edges
					eijk = cq->egsab[isab][i][ij].cross(cq->vtsab[isab][j][1]-cq->vtsab[isab][j][0]);
					if ( under_capacity < pow(eijk.norm(),4.0) )
						if ( separation(cq->vtsab[isab][i],cq->vtsab[isab][j],g,eijk) ) goto END;
				} // xd vs 1d
			}

		// end of xd vs yd ( 1 < x ) and begin of 1d vs yd
		} else if ( !cq->vtsab[isab][i].empty() ) {

			// for r(i)'s edge
			g = 0.5*(cq->vtsab[isab][i][1]+cq->vtsab[isab][i][0]);
			if ( separation(cq->vtsab[isab][i],cq->vtsab[isab][j],g,cq->vtsab[isab][i][1]-cq->vtsab[isab][i][0]) ) goto END;

			// loop for r(j)'s edges
			if ( !cq->egsab[isab][j].empty() ) { // 1d vs yd ( 1 < y )
				for ( int jk=0; jk<(int)cq->egsab[isab][j].size(); jk++ ) {
					g = 0.5*(cq->vtsab[isab][j][cq->vegsab[isab][j][jk][1]]+cq->vtsab[isab][j][cq->vegsab[isab][j][jk][0]]);
					if ( separation(cq->vtsab[isab][j],cq->vtsab[isab][i],g,cq->egsab[isab][j][jk]) ) goto END;
				}
			} else if ( !cq->vtsab[isab][j].empty() ) { // end of 1d vs yd ( 1 < y ) and begin 1d vs 1d
				g = 0.5*(cq->vtsab[isab][j][1]+cq->vtsab[isab][j][0]);
				if ( separation(cq->vtsab[isab][j],cq->vtsab[isab][i],g,cq->vtsab[isab][j][1]-cq->vtsab[isab][j][0]) ) goto END;
			} // 1d vs 1d

			// search more detail if the cross flag is still true
			g = 0.5*(cq->vtsab[isab][i][1]+cq->vtsab[isab][i][0]);
			if ( !cq->egsab[isab][j].empty() ) { // 1d vs yd ( 1 < y )
				// loop for r(j)'s edges
				for ( int jk=0; jk<(int)cq->egsab[isab][j].size(); jk++ ) {
					eijk = (cq->vtsab[isab][i][1]-cq->vtsab[isab][i][0]).cross(cq->egsab[isab][j][jk]);
					if ( under_capacity < pow(eijk.norm(),4.0) )
						if ( separation(cq->vtsab[isab][i],cq->vtsab[isab][j],g,eijk) ) goto END;
				}
			} else if ( !cq->vtsab[isab][j].empty() ) { // end of 1d vs yd ( 1 < y ) and begin of 1d vs 1d
				eijk = (cq->vtsab[isab][i][1]-cq->vtsab[isab][i][0]).cross(cq->vtsab[isab][j][1]-cq->vtsab[isab][j][0]);
				if ( under_capacity < pow(eijk.norm(),4.0) )
					if ( separation(cq->vtsab[isab][i],cq->vtsab[isab][j],g,eijk) ) goto END;
			} // 1d vs 1d

		} // 1d vs yd

	} // specified partial

	// two convex hulls are in collision
	coll = true;

END:

	return coll;
}

// count sharing forms
void csform::update ( int rank, int nmpi ) {
	// temporary variables
	bool face,edge,corner;
	int iab,isab,a_start;
	int local_start,local_size;		// mpi: temporary
	int jb,jf,je,jv;
	int *nb,*nf,*ne,*nv,**nbab,**nfab,**neab,**nvab,**nbsab,**nfsab,**nesab,**nvsab;
	int *nb0,*nf0,*ne0,*nv0,**nbab0,**nfab0,**neab0,**nvab0,**nbsab0,**nfsab0,**nesab0,**nvsab0;
	double *s,**sab,**ssab;
	double *eb,*ef,*ee,*ev,**ebab,**efab,**eeab,**evab,**ebsab,**efsab,**eesab,**evsab;
	vec v1(3),v2(3),g(3);
	vector < vector < vec > > geg,gnv;
	vector < vector < vector < vec > > > gegab,gnvab;
	vector < vector < vector < vec > > > gegsab,gnvsab;

	// memory allocation
	// *ab: partial information of *
	// n*: distribution of *s which calculated from e* 
	// e*: expected value of *s constructing its convex hull
	// nb,eb: * = bicap
	// nf,ef: * = facet
	// ne,ee: * = edge
	// nv,ev: * = vertice
	nb = alloc1d < int > (ncs);
	nf = alloc1d < int > (ncs);
	ne = alloc1d < int > (ncs);
	nv = alloc1d < int > (ncs);
	eb = alloc1d < double > (co->nions);
	ef = alloc1d < double > (co->nions);
	ee = alloc1d < double > (co->nions);
	ev = alloc1d < double > (co->nions);
	nbab = alloc2d < int > (co->ntypes*co->ntypes,ncs);
	nfab = alloc2d < int > (co->ntypes*co->ntypes,ncs);
	neab = alloc2d < int > (co->ntypes*co->ntypes,ncs);
	nvab = alloc2d < int > (co->ntypes*co->ntypes,ncs);
	ebab = alloc2d < double > (co->ntypes*co->ntypes,co->nions);
	efab = alloc2d < double > (co->ntypes*co->ntypes,co->nions);
	eeab = alloc2d < double > (co->ntypes*co->ntypes,co->nions);
	evab = alloc2d < double > (co->ntypes*co->ntypes,co->nions);
	nbsab = alloc2d < int > (co->nstypes*co->nstypes,ncs);
	nfsab = alloc2d < int > (co->nstypes*co->nstypes,ncs);
	nesab = alloc2d < int > (co->nstypes*co->nstypes,ncs);
	nvsab = alloc2d < int > (co->nstypes*co->nstypes,ncs);
	ebsab = alloc2d < double > (co->nstypes*co->nstypes,co->nions);
	efsab = alloc2d < double > (co->nstypes*co->nstypes,co->nions);
	eesab = alloc2d < double > (co->nstypes*co->nstypes,co->nions);
	evsab = alloc2d < double > (co->nstypes*co->nstypes,co->nions);
	// gathering memory
	nb0 = alloc1d < int > (ncs);
	nf0 = alloc1d < int > (ncs);
	ne0 = alloc1d < int > (ncs);
	nv0 = alloc1d < int > (ncs);
	nbab0 = alloc2d < int > (co->ntypes*co->ntypes,ncs);
	nfab0 = alloc2d < int > (co->ntypes*co->ntypes,ncs);
	neab0 = alloc2d < int > (co->ntypes*co->ntypes,ncs);
	nvab0 = alloc2d < int > (co->ntypes*co->ntypes,ncs);
	nbsab0 = alloc2d < int > (co->nstypes*co->nstypes,ncs);
	nfsab0 = alloc2d < int > (co->nstypes*co->nstypes,ncs);
	nesab0 = alloc2d < int > (co->nstypes*co->nstypes,ncs);
	nvsab0 = alloc2d < int > (co->nstypes*co->nstypes,ncs);
	// normalization factor: Voronoi surface area
	s = alloc1d < double > (co->nions);
	sab = alloc2d < double > (co->ntypes*co->ntypes,co->nions);
	ssab = alloc2d < double > (co->nstypes*co->nstypes,co->nions);
	// midpoint of edge and center point of the facet
	geg.resize(co->nions);
	gnv.resize(co->nions);
	gegab.resize(co->ntypes*co->ntypes);
	gnvab.resize(co->ntypes*co->ntypes);
	for ( iab=0; iab<(int)(co->ntypes*co->ntypes); iab++ ) {
		gegab[iab].resize(co->nions);
		gnvab[iab].resize(co->nions);
	}
	gegsab.resize(co->nstypes*co->nstypes);
	gnvsab.resize(co->nstypes*co->nstypes);
	for ( isab=0; isab<(int)(co->nstypes*co->nstypes); isab++ ) {
		gegsab[isab].resize(co->nions);
		gnvsab[isab].resize(co->nions);
	}

	// calc midpoint of edge and center point of facet
	for ( int i=0; i<(co->nions); i++ ) {
		// midpoint of edge
		for ( int ij=0; ij<(int)cq->eg[i].size(); ij++ ) {
			v1 = cq->vt[i][cq->veg[i][ij][0]];
			v2 = cq->vt[i][cq->veg[i][ij][1]];
			geg[i].push_back(0.5*(v1+v2));
		}
		// center of facet
		for ( int ij=0; ij<(int)cq->nv[i].size(); ij++ ) {
			vector < int > list;
			vector < vec > points;
			for ( int ijk=0; ijk<(int)cq->vnv[i][ij].size(); ijk++ )
				points.push_back(cq->vt[i][(int)cq->vnv[i][ij][ijk]]);
			if ( dim_convex(points) == 2 ) {
				calc_flat(points,&list,&g);
				gnv[i].push_back(g);
				vector < int > ().swap(list);
			}
			vector < vec > ().swap(points);
		}
	}
	// partial
	for ( int a=0; a<(co->ntypes); a++ ) { a_start = co->mystart(a);
	for ( int b=0; b<(co->ntypes); b++ ) { iab = a*co->ntypes+b;
	for ( int i=a_start; i<(a_start+co->item[a]); i++ ) {
		// midpoint of edge
		for ( int ij=0; ij<(int)cq->egab[iab][i].size(); ij++ ) {
			v1 = cq->vtab[iab][i][cq->vegab[iab][i][ij][0]];
			v2 = cq->vtab[iab][i][cq->vegab[iab][i][ij][1]];
			gegab[iab][i].push_back(0.5*(v1+v2));
		}
		// center of facet
		for ( int ij=0; ij<(int)cq->nvab[iab][i].size(); ij++ ) {
			vector < int > list;
			vector < vec > points;
			for ( int ijk=0; ijk<(int)cq->vnvab[iab][i][ij].size(); ijk++ )
				points.push_back(cq->vtab[iab][i][cq->vnvab[iab][i][ij][ijk]]);
			if ( dim_convex(points) == 2 ) {
				calc_flat(points,&list,&g);
				gnvab[iab][i].push_back(g);
//				// check
//				vec ev(3);
//				for ( int ijk=0; ijk<(int)list.size()+1; ijk++ ) {
//					ev = points[list[ijk%(int)list.size()]];
//					cout << ev(0) << " " << ev(1) << " " << ev(2) << " ";
//					cout << g(0) << " " << g(1) << " " << g(2) << endl;
//				}
				vector < int > ().swap(list);
			}
//			cout << endl;
			vector < vec > ().swap(points);
		}
//		cout << endl;
	}
	}
	}
	// specified partial
	for ( int spa=0; spa<(co->nstypes); spa++ ) {
	for ( int spb=0; spb<(co->nstypes); spb++ ) { isab = spa*co->nstypes+spb;
	for ( int ia=0; ia<(int)co->stypeindex[spa].size(); ia++ ) {
		int a=co->stypeindex[spa][ia];a_start=co->mystart(a);
		for ( int i=a_start; i<(a_start+co->item[a]); i++ ) {
			// midpoint of edge
			for ( int ij=0; ij<(int)cq->egsab[isab][i].size(); ij++ ) {
				v1 = cq->vtsab[isab][i][cq->vegsab[isab][i][ij][0]];
				v2 = cq->vtsab[isab][i][cq->vegsab[isab][i][ij][1]];
				gegsab[isab][i].push_back(0.5*(v1+v2));
			}
			// center of facet
			for ( int ij=0; ij<(int)cq->nvsab[isab][i].size(); ij++ ) {
				vector < int > list;
				vector < vec > points;
				for ( int ijk=0; ijk<(int)cq->vnvsab[isab][i][ij].size(); ijk++ )
					points.push_back(cq->vtsab[isab][i][cq->vnvsab[isab][i][ij][ijk]]);
				if ( dim_convex(points) == 2 ) {
					calc_flat(points,&list,&g);
					gnvsab[isab][i].push_back(g);
					vector < int > ().swap(list);
				}
				vector < vec > ().swap(points);
			}
		}
	}
	}
	}
	MPI_Barrier(MPI_COMM_WORLD);

	/*--- search sharing form
	   r(i)'s neighboring ions construct 3d convex hull
	   ( it corresponds to the condition of "ddim = 3" )
	   priority: bicap- ( > non- ) > face- > edge- > corner- > non-sharing
		| if: two clusters are in collision ( bicap-sharing )
		|	possible sharing form is only bicap
		|	and only count sharing vertices
		| else if: no collision ( non-sharing )
		|	do nothing;
		|	two clusters are not contacted
		| else: contacting at its boundary
			| if: two facets are same ( face-sharing )
			|	count it with average Voronoi weightings
			|	related to vertices that construct r(i)'s facet
			| else: no facets are same
				| if: two edges are same ( edge-sharing )
				|	count it with average Voronoi weightings
				|	related to vertices that construct r(i)'s edge
				| else: no edges are same
					| if: two vertices are same ( corner-sharing )
					|	count it with average Voronoi weightings
					| else: no vertices are same
					|	do nothing
					|	two cluster are not contacted
	   r(i)'s neighboring ions construct 2d convex hull
	   ( it corresponds to the condition of "ddim = 2" )
	   priority: bicap- ( > non- ) > edge- > corner- > non-sharing
		| if: two clusters are in collision ( bicap-sharing )
		|     ( * this case may be occurred at few probability )
		|	possible sharing form is only bicap
		|	and only count sharing vertices
		| else if: no collision ( non-sharing )
		|	do nothing;
		|	two clusters are not contacted
		| else: contacting at its boundary
			| if: two edges are same ( edge-sharing )
			|	count it with average Voronoi weightings
			|	related to vertices that construct r(i)'s edge
			| else: no edges are same
				| if: two vertices are same ( corner-sharing )
				|	count it with average Voronoi weightings
				| else: no vertices are same
				|	do nothing
				|	two cluster are not contacted
	   r(i)'s neighboring ions construct 1d convex hull
	   ( it corresponds to the condition of "ddim = 1" )
	   priority: bicap- ( > non- ) > corner- > non-sharing
		| if: two clusters are in collision ( bicap-sharing )
		|     ( * this case may be occurred at few probability )
		|	possible sharing form is only bicap
		|	and only count sharing vertices
		| else if: no collision ( non-sharing )
		|	do nothing;
		|	two clusters are not contacted
		| else: contacting at its boundary
			| if: two vertices are same ( corner-sharing )
			|	count it with average Voronoi weightings
			| else: no vertices are same
			|	do nothing
			|	two cluster are not contacted
	-----*/

	// mpi parameter
	mpi_local_part(&local_start,&local_size,co->nions,rank,nmpi);

	// total information
	for ( int i=local_start; i<local_start+local_size; i++ ) {
		for ( int j=0; j<(co->nions); j++ ) {
		if ( i != j ) {
			// reset boolean
			face=false;edge=false;corner=false;
			// check whether two convex hulls are in collision
			if ( collision (i,j,-1,-1) ) {
				for ( int ij=0; ij<(int)cq->vt[i].size(); ij++ ) {
				for ( int jk=0; jk<(int)cq->vt[j].size(); jk++ ) {
					if ( pow((cq->vt[i][ij]-cq->vt[j][jk]).norm(),2.0) < under_capacity ) {
							eb[i]+=cv->surf[i][ij];s[i]+=cv->surf[i][ij];
					}
				}
				}
			// othercase: face/edge/corner
			} else {
				// r(i)'s and r(j)'s facets
				for ( int ij=0; ij<(int)gnv[i].size(); ij++ ) {
				for ( int jk=0; jk<(int)gnv[j].size(); jk++ ) {
				if ( pow((gnv[i][ij]-gnv[j][jk]).norm(),2.0) < under_capacity ) {
					// accumulate surface area
					face = true;
					for ( int ijk=0; ijk<(int)cq->vnv[i][ij].size(); ijk++ ) {
						ef[i] += cv->surf[i][cq->vnv[i][ij][ijk]];
						s[i] += cv->surf[i][cq->vnv[i][ij][ijk]];
					}
				}
				}
				}
				// othercase: edge/corner/none
				if ( !face ) {
					// r(i)'s and r(j)'s edges
					for ( int ij=0; ij<(int)geg[i].size(); ij++ ) {
					for ( int jk=0; jk<(int)geg[j].size(); jk++ ) {
					if ( pow((geg[i][ij]-geg[j][jk]).norm(),2.0) < under_capacity ) {
						// accumulate surface area
						edge = true;
						for ( int ijk=0; ijk<(int)cq->veg[i][ij].size(); ijk++ ) {
							ee[i] += cv->surf[i][cq->veg[i][ij][ijk]];
							s[i] += cv->surf[i][cq->veg[i][ij][ijk]];
						}
					}
					}
					}
					// othercase: corner/none
					if ( !edge ) {
						// r(i)'s and r(j)'s vertices
						for ( int ij=0; ij<(int)cq->vt[i].size(); ij++ ) {
						for ( int jk=0; jk<(int)cq->vt[j].size(); jk++ ) {
						if ( pow((cq->vt[i][ij]-cq->vt[j][jk]).norm(),2.0) < under_capacity ) {
							// accumulate surface area
							corner = true;
							ev[i]+=cv->surf[i][ij];s[i]+=cv->surf[i][ij];
						}
						}
						}
					} // corner
				} // edge
			} // bicap or face

		} // if ( i != j )
		} // j-th atom

		// expected value
		if ( under_capacity < pow(s[i],4.0) ) {
			eb[i]/=s[i];ef[i]/=s[i];ee[i]/=s[i];ev[i]/=s[i];
			// counting values
			jb = (int)((eb[i]-co->cs0)/co->dcs);
			jf = (int)((ef[i]-co->cs0)/co->dcs);
			je = (int)((ee[i]-co->cs0)/co->dcs);
			jv = (int)((ev[i]-co->cs0)/co->dcs);
			if ( 0 <= jb && jb < ncs ) nb[jb] ++;
			if ( 0 <= jf && jf < ncs ) nf[jf] ++;
			if ( 0 <= je && je < ncs ) ne[je] ++;
			if ( 0 <= jv && jv < ncs ) nv[jv] ++;
		}

	}
	MPI_Barrier(MPI_COMM_WORLD);

	// partial
	for ( int a=0; a<(co->ntypes); a++ ) { a_start = co->mystart(a);
	for ( int b=0; b<(co->ntypes); b++ ) { iab = a*co->ntypes+b;

		// mpi parameter
		mpi_local_part(&local_start,&local_size,co->item[a],rank,nmpi);
	
		for ( int i=a_start+local_start; i<a_start+local_start+local_size; i++ ) {
			for ( int j=a_start; j<(a_start+co->item[a]); j++ ) {
			if ( i != j ) {
				// reset boolean
				face=false;edge=false;corner=false;
				// check whether two convex hulls are in collision
				if ( collision (i,j,iab,-1) ) {
					for ( int ij=0; ij<(int)cq->vtab[iab][i].size(); ij++ ) {
					for ( int jk=0; jk<(int)cq->vtab[iab][j].size(); jk++ ) {
						if ( pow((cq->vtab[iab][i][ij]-cq->vtab[iab][j][jk]).norm(),4.0) < under_capacity ) {
							ebab[iab][i] += cv->surfab[iab][i][ij];
							sab[iab][i] += cv->surfab[iab][i][ij];
						}
					}
					}
				// othercase: face/edge/corner
				} else {
					// r(i)'s and r(j)'s facets
					for ( int ij=0; ij<(int)gnvab[iab][i].size(); ij++ ) {
					for ( int jk=0; jk<(int)gnvab[iab][j].size(); jk++ ) {
					if ( pow((gnvab[iab][i][ij]-gnvab[iab][j][jk]).norm(),4.0) < under_capacity ) {
						// accumulate surface area
						face = true;
						for ( int ijk=0; ijk<(int)cq->vnvab[iab][i][ij].size(); ijk++ ) {
							efab[iab][i] += cv->surfab[iab][i][cq->vnvab[iab][i][ij][ijk]];
							sab[iab][i] += cv->surfab[iab][i][cq->vnvab[iab][i][ij][ijk]];
						}
					}
					}
					}
					// othercase: edge/corner/none
					if ( !face ) {
						// r(i)'s and r(j)'s edges
						for ( int ij=0; ij<(int)gegab[iab][i].size(); ij++ ) {
						for ( int jk=0; jk<(int)gegab[iab][j].size(); jk++ ) {
						if ( pow((gegab[iab][i][ij]-gegab[iab][j][jk]).norm(),4.0) < under_capacity ) {
							// accumulate surface area
							edge = true;
							for ( int ijk=0; ijk<(int)cq->vegab[iab][i][ij].size(); ijk++ ) {
								eeab[iab][i] += cv->surfab[iab][i][cq->vegab[iab][i][ij][ijk]];
								sab[iab][i] += cv->surfab[iab][i][cq->vegab[iab][i][ij][ijk]];
							}
						}
						}
						}
						// othercase: corner/none
						if ( !edge ) {
							// r(i)'s and r(j)'s vertices
							for ( int ij=0; ij<(int)cq->vtab[iab][i].size(); ij++ ) {
							for ( int jk=0; jk<(int)cq->vtab[iab][j].size(); jk++ ) {
							if ( pow((cq->vtab[iab][i][ij]-cq->vtab[iab][j][jk]).norm(),4.0) < under_capacity ) {
								// accumulate surface area
								corner = true;
								evab[iab][i] += cv->surfab[iab][i][ij];
								sab[iab][i] += cv->surfab[iab][i][ij];
							}
							}
							}
						} // corner
					} // edge
				} // bicap or face
	
			} // if ( i != j )
			} // j-th atom

//			if ( rank == 0 ) {
//				cout << "pair: " << a+1 << "-" << b+1 << endl;
//				cout << "atom(" << i+1 << "):"<< endl;
//				cout << "\t bicap -> " << ebab[iab][i] << endl;
//				cout << "\t face  -> " << efab[iab][i] << endl;
//				cout << "\t edge  -> " << eeab[iab][i] << endl;
//				cout << "\tcorner -> " << evab[iab][i] << endl;
//			}

			// expected value
			if ( under_capacity < pow(sab[iab][i],2.0) ) {
				ebab[iab][i]/=sab[iab][i];efab[iab][i]/=sab[iab][i];
				eeab[iab][i]/=sab[iab][i];evab[iab][i]/=sab[iab][i];
				// counting values
				jb = (int)((ebab[iab][i]-co->cs0)/co->dcs);
				jf = (int)((efab[iab][i]-co->cs0)/co->dcs);
				je = (int)((eeab[iab][i]-co->cs0)/co->dcs);
				jv = (int)((evab[iab][i]-co->cs0)/co->dcs);
				if ( 0 <= jb && jb < ncs ) nbab[iab][jb] ++;
				if ( 0 <= jf && jf < ncs ) nfab[iab][jf] ++;
				if ( 0 <= je && je < ncs ) neab[iab][je] ++;
				if ( 0 <= jv && jv < ncs ) nvab[iab][jv] ++;
			}

		}
	}
	}
	MPI_Barrier(MPI_COMM_WORLD);

	// specified partial
	for ( int spa=0; spa<(co->nstypes); spa++ ) {
	for ( int spb=0; spb<(co->nstypes); spb++ ) { isab = spa*co->nstypes+spb;
	for ( int ia=0; ia<(int)co->stypeindex[spa].size(); ia++ ) {
		int a=co->stypeindex[spa][ia];a_start=co->mystart(a);

		// mpi parameter
		mpi_local_part(&local_start,&local_size,co->item[a],rank,nmpi);

		for ( int i=a_start+local_start; i<a_start+local_start+local_size; i++ ) {
			for ( int j=a_start; j<(a_start+co->item[a]); j++ ) {
			if ( i != j ) {
				// reset boolean
				face=false;edge=false;corner=false;
				// check whether two convex hulls are in collision
				if ( collision (i,j,-1,isab) ) {
					for ( int ij=0; ij<(int)cq->vtsab[isab][i].size(); ij++ ) {
					for ( int jk=0; jk<(int)cq->vtsab[isab][j].size(); jk++ ) {
						if ( pow((cq->vtsab[isab][i][ij]-cq->vtsab[isab][j][jk]).norm(),4.0) < under_capacity ) {
							ebsab[isab][i] += cv->surfsab[isab][i][ij];
							ssab[isab][i] += cv->surfsab[isab][i][ij];
						}
					}
					}
				// othercase: face/edge/corner
				} else {
					// r(i)'s and r(j)'s facets
					for ( int ij=0; ij<(int)gnvsab[isab][i].size(); ij++ ) {
					for ( int jk=0; jk<(int)gnvsab[isab][j].size(); jk++ ) {
					if ( pow((gnvsab[isab][i][ij]-gnvsab[isab][j][jk]).norm(),4.0) < under_capacity ) {
						// accumulate surface area
						face = true;
						for ( int ijk=0; ijk<(int)cq->vnvsab[isab][i][ij].size(); ijk++ ) {
							efsab[isab][i] += cv->surfsab[isab][i][cq->vnvsab[isab][i][ij][ijk]];
							ssab[isab][i] += cv->surfsab[isab][i][cq->vnvsab[isab][i][ij][ijk]];
						}
					}
					}
					}
					// othercase: edge/corner/none
					if ( !face ) {
						// r(i)'s and r(j)'s edges
						for ( int ij=0; ij<(int)gegsab[isab][i].size(); ij++ ) {
						for ( int jk=0; jk<(int)gegsab[isab][j].size(); jk++ ) {
						if ( pow((gegsab[isab][i][ij]-gegsab[isab][j][jk]).norm(),4.0) < under_capacity ) {
							// accumulate surface area
							edge = true;
							for ( int ijk=0; ijk<(int)cq->vegsab[isab][i][ij].size(); ijk++ ) {
								eesab[isab][i] += cv->surfsab[isab][i][cq->vegsab[isab][i][ij][ijk]];
								ssab[isab][i] += cv->surfsab[isab][i][cq->vegsab[isab][i][ij][ijk]];
							}
						}
						}
						}
						// othercase: corner/none
						if ( !edge ) {
							// r(i)'s and r(j)'s vertices
							for ( int ij=0; ij<(int)cq->vtsab[isab][i].size(); ij++ ) {
							for ( int jk=0; jk<(int)cq->vtsab[isab][j].size(); jk++ ) {
							if ( pow((cq->vtsab[isab][i][ij]-cq->vtsab[isab][j][jk]).norm(),4.0) < under_capacity ) {
								// accumulate surface area
								corner = true;
								evsab[isab][i] += cv->surfsab[isab][i][ij];
								ssab[isab][i] += cv->surfsab[isab][i][ij];
							}
							}
							}
						} // corner
					} // edge
				} // bicap or face
	
			} // if ( i != j )
			} // j-th atom
		
//			if ( rank == 0 ) {
//				cout << "pair: " << spa+1 << "-" << spb+1 << endl;
//				cout << "atom(" << i+1 << "):"<< endl;
//				cout << "\t bicap -> " << ebsab[isab][i] << endl;
//				cout << "\t face  -> " << efsab[isab][i] << endl;
//				cout << "\t edge  -> " << eesab[isab][i] << endl;
//				cout << "\tcorner -> " << evsab[isab][i] << endl;
//			}

			// expected value
			if ( under_capacity < pow(ssab[isab][i],2.0) ) {
				ebsab[isab][i]/=ssab[isab][i];efsab[isab][i]/=ssab[isab][i];
				eesab[isab][i]/=ssab[isab][i];evsab[isab][i]/=ssab[isab][i];
				// counting values
				jb = (int)((ebsab[isab][i]-co->cs0)/co->dcs);
				jf = (int)((efsab[isab][i]-co->cs0)/co->dcs);
				je = (int)((eesab[isab][i]-co->cs0)/co->dcs);
				jv = (int)((evsab[isab][i]-co->cs0)/co->dcs);
				if ( 0 <= jb && jb < ncs ) nbsab[isab][jb] ++;
				if ( 0 <= jf && jf < ncs ) nfsab[isab][jf] ++;
				if ( 0 <= je && je < ncs ) nesab[isab][je] ++;
				if ( 0 <= jv && jv < ncs ) nvsab[isab][jv] ++;
			}

		}
	}
	}
	}
	MPI_Barrier(MPI_COMM_WORLD);

	// data gathering
	MPI_Allreduce(&(nb[0]),&(nb0[0]),ncs,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(&(nf[0]),&(nf0[0]),ncs,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(&(ne[0]),&(ne0[0]),ncs,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(&(nv[0]),&(nv0[0]),ncs,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(&(nbab[0][0]),&(nbab0[0][0]),co->ntypes*co->ntypes*ncs,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(&(nfab[0][0]),&(nfab0[0][0]),co->ntypes*co->ntypes*ncs,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(&(neab[0][0]),&(neab0[0][0]),co->ntypes*co->ntypes*ncs,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(&(nvab[0][0]),&(nvab0[0][0]),co->ntypes*co->ntypes*ncs,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(&(nbsab[0][0]),&(nbsab0[0][0]),co->nstypes*co->nstypes*ncs,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(&(nfsab[0][0]),&(nfsab0[0][0]),co->nstypes*co->nstypes*ncs,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(&(nesab[0][0]),&(nesab0[0][0]),co->nstypes*co->nstypes*ncs,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(&(nvsab[0][0]),&(nvsab0[0][0]),co->nstypes*co->nstypes*ncs,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	// update
	for ( jb=0; jb<ncs; jb++ )
		pb[jb] = co->update(pb[jb],(double)nb0[jb]/(double)(co->nions*co->dcs));
	for ( jf=0; jf<ncs; jf++ )
		pf[jf] = co->update(pf[jf],(double)nf0[jf]/(double)(co->nions*co->dcs));
	for ( je=0; je<ncs; je++ )
		pe[je] = co->update(pe[je],(double)ne0[je]/(double)(co->nions*co->dcs));
	for ( jv=0; jv<ncs; jv++ )
		pv[jv] = co->update(pv[jv],(double)nv0[jv]/(double)(co->nions*co->dcs));
	for ( int a=0; a<(co->ntypes); a++ ) {
	for ( int b=0; b<(co->ntypes); b++ ) { iab = a*co->ntypes+b;
		for ( jb=0; jb<ncs; jb++ )
			pbab[iab][jb] = co->update(pbab[iab][jb],(double)nbab0[iab][jb]/(double)(co->item[a]*co->dcs));
		for ( jf=0; jf<ncs; jf++ )
			pfab[iab][jf] = co->update(pfab[iab][jf],(double)nfab0[iab][jf]/(double)(co->item[a]*co->dcs));
		for ( je=0; je<ncs; je++ )
			peab[iab][je] = co->update(peab[iab][je],(double)neab0[iab][je]/(double)(co->item[a]*co->dcs));
		for ( jv=0; jv<ncs; jv++ )
			pvab[iab][jv] = co->update(pvab[iab][jv],(double)nvab0[iab][jv]/(double)(co->item[a]*co->dcs));
	}
	}
	for ( int spa=0; spa<(co->nstypes); spa++ ) {
	for ( int spb=0; spb<(co->nstypes); spb++ ) { isab = spa*co->nstypes+spb;
		for ( jb=0; jb<ncs; jb++ )
			pbsab[isab][jb] = co->update(pbsab[isab][jb],(double)nbsab0[isab][jb]/(double)(co->sitem[spa]*co->dcs));
		for ( jf=0; jf<ncs; jf++ )
			pfsab[isab][jf] = co->update(pfsab[isab][jf],(double)nfsab0[isab][jf]/(double)(co->sitem[spa]*co->dcs));
		for ( je=0; je<ncs; je++ )
			pesab[isab][je] = co->update(pesab[isab][je],(double)nesab0[isab][je]/(double)(co->sitem[spa]*co->dcs));
		for ( jv=0; jv<ncs; jv++ )
			pvsab[isab][jv] = co->update(pvsab[isab][jv],(double)nvsab0[isab][jv]/(double)(co->sitem[spa]*co->dcs));
	}
	}
	MPI_Barrier(MPI_COMM_WORLD);

	// free
	free1d < int > (nb);
	free1d < int > (nf);
	free1d < int > (ne);
	free1d < int > (nv);
	free1d < double > (eb);
	free1d < double > (ef);
	free1d < double > (ee);
	free1d < double > (ev);
	free2d < int > (nbab);
	free2d < int > (nfab);
	free2d < int > (neab);
	free2d < int > (nvab);
	free2d < double > (ebab);
	free2d < double > (efab);
	free2d < double > (eeab);
	free2d < double > (evab);
	free2d < int > (nbsab);
	free2d < int > (nfsab);
	free2d < int > (nesab);
	free2d < int > (nvsab);
	free2d < double > (ebsab);
	free2d < double > (efsab);
	free2d < double > (eesab);
	free2d < double > (evsab);

	free1d < int > (nb0);
	free1d < int > (nf0);
	free1d < int > (ne0);
	free1d < int > (nv0);
	free2d < int > (nbab0);
	free2d < int > (nfab0);
	free2d < int > (neab0);
	free2d < int > (nvab0);
	free2d < int > (nbsab0);
	free2d < int > (nfsab0);
	free2d < int > (nesab0);
	free2d < int > (nvsab0);

	free1d < double > (s);
	free2d < double > (sab);
	free2d < double > (ssab);
	vector < vector < vec > > ().swap(geg);
	vector < vector < vec > > ().swap(gnv);
	vector < vector < vector < vec > > > ().swap(gegab);
	vector < vector < vector < vec > > > ().swap(gnvab);
	vector < vector < vector < vec > > > ().swap(gegsab);
	vector < vector < vector < vec > > > ().swap(gnvsab);

}

void csform::write ( string filename ) {
	double cs;
	ofstream ofs,ofsab,ofssab;

	// write out
	ofs.open((filename+"-tot.dat").c_str(),ios::out);
	ofsab.open((filename+".dat").c_str(),ios::out);
	ofssab.open((filename+"-spec.dat").c_str(),ios::out);
	ofs << scientific << setprecision(5);
	ofsab << scientific << setprecision(5);
	ofssab << scientific << setprecision(5);
	// total
	for ( int ics=0; ics<ncs; ics++ ) {
		cs = (double)ics*co->dcs + co->cs0;
		ofs << cs;
		ofs << "\t" << pb[ics] << "\t" << pf[ics];
		ofs << "\t" << pe[ics] << "\t" << pv[ics];
		ofs << endl;
	}
	// partial
	for ( int iab=0; iab<(co->ntypes*co->ntypes); iab++ ) {
		for ( int ics=0; ics<ncs; ics++ ) {
			cs = (double)ics*co->dcs + co->cs0;
			ofsab << cs;
			ofsab << "\t" << pbab[iab][ics] << "\t" << pfab[iab][ics];
			ofsab << "\t" << peab[iab][ics] << "\t" << pvab[iab][ics];
			ofsab << endl;
		}
		ofsab << endl;
	}
	// specified partial
	for ( int isab=0; isab<(co->nstypes*co->nstypes); isab++ ) {
		for ( int ics=0; ics<ncs; ics++ ) {
			cs = (double)ics*co->dcs + co->cs0;
			ofssab << cs;
			ofssab << "\t" << pbsab[isab][ics] << "\t" << pfsab[isab][ics];
			ofssab << "\t" << pesab[isab][ics] << "\t" << pvsab[isab][ics];
			ofssab << endl;
		}
		ofssab << endl;
	}
	ofs.close();
	ofsab.close();
	ofssab.close();

}
