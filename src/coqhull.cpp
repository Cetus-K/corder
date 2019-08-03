#include "corder.h"

// destructor
coqhull::~coqhull() {
	vector < vector < vec > > ().swap(nv);
	vector < vector < vec > > ().swap(eg);
	vector < vector < vec > > ().swap(vt);
	vector < vector < vector < int > > > ().swap(vnv);
	vector < vector < vector < int > > > ().swap(veg);
	vector < vector < vector < vec > > > ().swap(nvab);
	vector < vector < vector < vec > > > ().swap(egab);
	vector < vector < vector < vec > > > ().swap(vtab);
	vector < vector < vector < vector < int > > > > ().swap(vnvab);
	vector < vector < vector < vector < int > > > > ().swap(vegab);
}

void coqhull::init ( void ) {

	// memory allocation
	// total information
	nv.resize(co->nions);
	eg.resize(co->nions);
	vt.resize(co->nions);
	vnv.resize(co->nions);
	veg.resize(co->nions);
	// partial information
	nvab.resize(co->ntypes*co->ntypes);
	egab.resize(co->ntypes*co->ntypes);
	vtab.resize(co->ntypes*co->ntypes);
	vnvab.resize(co->ntypes*co->ntypes);
	vegab.resize(co->ntypes*co->ntypes);
	for ( int iab=0; iab<(co->ntypes*co->ntypes); iab++ ) {
		nvab[iab].resize(co->nions);
		egab[iab].resize(co->nions);
		vtab[iab].resize(co->nions);
		vnvab[iab].resize(co->nions);
		vegab[iab].resize(co->nions);
	}
	// specified partial information
	nvsab.resize(co->nstypes*co->nstypes);
	egsab.resize(co->nstypes*co->nstypes);
	vtsab.resize(co->nstypes*co->nstypes);
	vnvsab.resize(co->nstypes*co->nstypes);
	vegsab.resize(co->nstypes*co->nstypes);
	for ( int isab=0; isab<(co->nstypes*co->nstypes); isab++ ) {
		nvsab[isab].resize(co->nions);
		egsab[isab].resize(co->nions);
		vtsab[isab].resize(co->nions);
		vnvsab[isab].resize(co->nions);
		vegsab[isab].resize(co->nions);
	}

}

// set corder copy constructor and
// get informations that qhull calculated
// get values: normal vectors of facets, edges and vertices
void coqhull::get_one ( vector < vec > points, int index ) {

	// qhull inputs
	int dim=3;
	realT pts[256];
	facetT *facet,*neighbor;
	ridgeT *ridge,**ridgep;
	vertexT *vertex,**vertexp;
	setT *vertices;
	// dealing variables
	int ecnt;
	vector < int > vnijk,veijk(2);
	vector < vec > eijk(2);
	vec normal(3),edge(3);

	// free if nv[index] and eg[index] are remaining
	if ( !nv[index].empty() ) {
		vector < vec > ().swap(nv[index]);
		vector < vector < int > > ().swap(vnv[index]);
	}
	if ( !eg[index].empty() ) {
		vector < vec > ().swap(eg[index]);
		vector < vector < int > > ().swap(veg[index]);
	}
	// set points
	for ( int ij=0; ij<(int)points.size(); ij++ ) {
	for ( int ijk=0; ijk<dim; ijk++ ) {
		pts[ij*dim+ijk] = static_cast<realT>(points[ij](ijk));
	}
	}
	// set qhull and execute
	Qhull qhl("",dim,(int)points.size(),pts,"Q2 Pp");
	// gets:
	// 	normal vectors of facets;	{f1,f2,f3,...}
	// 	edges of convex hull;	{e1,e2,e3,...}
	//	vertices of convex hull;	{v1,v2,v3,...}
	// 	vertices constructing nv;	{vf11,vf12,...},{vf21,vf22,...},...
	// 	vertices constructing eg;	{ve11,ve12,...},{ve21,ve22,...},...
	qhl.qh()->visit_id++;
	FORALLfacet_(qhl.qh()->facet_list) {
		// vertices constructing the facet in ordering
		vertices = qh_facet3vertex(qhl.qh(),facet);
		FOREACHvertex_(vertices)
			vnijk.push_back(qh_pointid(qhl.qh(),vertex->point));
		if ( !vnijk.empty() ) {
			vnv[index].push_back(vnijk);
			vector < int > ().swap(vnijk);
		}
		// normals
		for ( int fi=0; fi<3; fi++ ) normal(fi) = facet->normal[fi];
		nv[index].push_back(normal);
		// ridges constructing the facet
		qh_makeridges(qhl.qh(),facet);
		facet->visitid = qhl.qh()->visit_id;
		FOREACHridge_(facet->ridges) {
			neighbor = otherfacet_(ridge,facet);
			if ( neighbor->visitid != qhl.qh()->visit_id ) {
				ecnt=0;
				FOREACHvertex_(ridge->vertices) {
					eijk[ecnt] = points[qh_pointid(qhl.qh(),vertex->point)];
					veijk[ecnt] = qh_pointid(qhl.qh(),vertex->point);
					ecnt ++;
				}
				if ( !veijk.empty() ) {
					edge = 1.0/(eijk[0]-eijk[1]).norm()*(eijk[0]-eijk[1]);
					eg[index].push_back(edge);
					veg[index].push_back(veijk);
				}
			}
		}
	}
	// free
	qh_settempfree(qhl.qh(),&vertices);

}

void coqhull::pget_one ( vector < vec > points, int index, int a, int b ) {

	// qhull inputs
	int dim=3;
	realT pts[256];
	facetT *facet,*neighbor;
	ridgeT *ridge,**ridgep;
	vertexT *vertex,**vertexp;
	setT *vertices;
	// dealing variables
	int iab=a*co->ntypes+b;
	int ecnt;
	vector < int > vnijk,veijk(2);
	vector < vec > eijk(2);
	vec normal(3),edge(3);

	// free if nv[index] and eg[index] are remaining
	if ( !nvab[iab][index].empty() ) {
		vector < vec > ().swap(nvab[iab][index]);
		vector < vector < int > > ().swap(vnvab[iab][index]);
	}
	if ( !egab[iab][index].empty() ) {
		vector < vec > ().swap(egab[iab][index]);
		vector < vector < int > > ().swap(vegab[iab][index]);
	}
	// set points
	for ( int ij=0; ij<(int)points.size(); ij++ ) {
	for ( int ijk=0; ijk<dim; ijk++ ) {
		pts[ij*dim+ijk] = static_cast<realT>(points[ij](ijk));
	}
	}
	// set qhull and execute
	Qhull qhl("",dim,(int)points.size(),pts,"Q2 Pp");
	// gets:
	// 	normal vectors of facets;	{f1,f2,f3,...}
	// 	edges of convex hull;	{e1,e2,e3,...}
	//	vertices of convex hull;	{v1,v2,v3,...}
	// 	vertices constructing nv;	{vf11,vf12,...},{vf21,vf22,...},...
	// 	vertices constructing eg;	{ve11,ve12,...},{ve21,ve22,...},...
	qhl.qh()->visit_id++;
	FORALLfacet_(qhl.qh()->facet_list) {
		// vertices constructing the facet in ordering
		vertices = qh_facet3vertex(qhl.qh(),facet);
		FOREACHvertex_(vertices)
			vnijk.push_back(qh_pointid(qhl.qh(),vertex->point));
		if ( !vnijk.empty() ) {
			vnvab[iab][index].push_back(vnijk);
			vector < int > ().swap(vnijk);
		}
		// normals
		for ( int fi=0; fi<3; fi++ ) normal(fi) = facet->normal[fi];
		nvab[iab][index].push_back(normal);
		// ridges constructing the facet
		qh_makeridges(qhl.qh(),facet);
		facet->visitid = qhl.qh()->visit_id;
		FOREACHridge_(facet->ridges) {
			neighbor = otherfacet_(ridge,facet);
			if ( neighbor->visitid != qhl.qh()->visit_id ) {
				ecnt=0;
				FOREACHvertex_(ridge->vertices) {
					eijk[ecnt] = points[qh_pointid(qhl.qh(),vertex->point)];
					veijk[ecnt] = qh_pointid(qhl.qh(),vertex->point);
					ecnt ++;
				}
				if ( !veijk.empty() ) {
					edge = 1.0/(eijk[0]-eijk[1]).norm()*(eijk[0]-eijk[1]);
					egab[iab][index].push_back(edge);
					vegab[iab][index].push_back(veijk);
				}
			}
		}
	}
	// free
	qh_settempfree(qhl.qh(),&vertices);

}

void coqhull::spget_one ( vector < vec > points, int index, int spa, int spb ) {

	// qhull inputs
	int dim=3;
	realT pts[256];
	facetT *facet,*neighbor;
	ridgeT *ridge,**ridgep;
	vertexT *vertex,**vertexp;
	setT *vertices;
	// dealing variables
	int isab=spa*co->nstypes+spb;
	int ecnt;
	vector < int > vnijk,veijk(2);
	vector < vec > eijk(2);
	vec normal(3),edge(3);

	// free if nv[index] and eg[index] are remaining
	if ( !nvsab[isab][index].empty() ) {
		vector < vec > ().swap(nvsab[isab][index]);
		vector < vector < int > > ().swap(vnvsab[isab][index]);
	}
	if ( !egsab[isab][index].empty() ) {
		vector < vec > ().swap(egsab[isab][index]);
		vector < vector < int > > ().swap(vegsab[isab][index]);
	}
	// set points
	for ( int ij=0; ij<(int)points.size(); ij++ ) {
	for ( int ijk=0; ijk<dim; ijk++ ) {
		pts[ij*dim+ijk] = static_cast<realT>(points[ij](ijk));
	}
	}
	// set qhull and execute
	Qhull qhl("",dim,(int)points.size(),pts,"Q2 Pp");
	// gets:
	// 	normal vectors of facets;	{f1,f2,f3,...}
	// 	edges of convex hull;	{e1,e2,e3,...}
	//	vertices of convex hull;	{v1,v2,v3,...}
	// 	vertices constructing nv;	{vf11,vf12,...},{vf21,vf22,...},...
	// 	vertices constructing eg;	{ve11,ve12,...},{ve21,ve22,...},...
	qhl.qh()->visit_id++;
	FORALLfacet_(qhl.qh()->facet_list) {
		// vertices constructing the facet in ordering
		vertices = qh_facet3vertex(qhl.qh(),facet);
		FOREACHvertex_(vertices)
			vnijk.push_back(qh_pointid(qhl.qh(),vertex->point));
		if ( !vnijk.empty() ) {
			vnvsab[isab][index].push_back(vnijk);
			vector < int > ().swap(vnijk);
		}
		// normals
		for ( int fi=0; fi<3; fi++ ) normal(fi) = facet->normal[fi];
		nvsab[isab][index].push_back(normal);
		// ridges constructing the facet
		qh_makeridges(qhl.qh(),facet);
		facet->visitid = qhl.qh()->visit_id;
		FOREACHridge_(facet->ridges) {
			neighbor = otherfacet_(ridge,facet);
			if ( neighbor->visitid != qhl.qh()->visit_id ) {
				ecnt=0;
				FOREACHvertex_(ridge->vertices) {
					eijk[ecnt] = points[qh_pointid(qhl.qh(),vertex->point)];
					veijk[ecnt] = qh_pointid(qhl.qh(),vertex->point);
					ecnt ++;
				}
				if ( !veijk.empty() ) {
					edge = 1.0/(eijk[0]-eijk[1]).norm()*(eijk[0]-eijk[1]);
					egsab[isab][index].push_back(edge);
					vegsab[isab][index].push_back(veijk);
				}
			}
		}
	}
	// free
	qh_settempfree(qhl.qh(),&vertices);

}

// get each geometric parameters
// vt[i]: i-th convex's vertices
// eg[i]: i-th convex's edges
// nv[i]: i-th convex's facet normal vectors
void coqhull::get ( vector < vector < vec > > points ) {
	// loop for all atoms
	for ( int i=0; i<(co->nions); i++ ) {
		if ( 2<(int)points[i].size() ) get_one(points[i],i);
	}
}
void coqhull::pget ( vector < vector < vec > > points, int a, int b ) {
	int a_start=co->mystart(a),b_start=co->mystart(b);
	// loop for all atoms
	for ( int a=0; a<(co->ntypes); a++ ) { a_start = co->mystart(a);
	for ( int i=a_start; i<(a_start+co->item[a]); i++ ) {
		if ( 2<(int)points[i].size() ) pget_one(points[i],i,a,b);
	}
	}
}
void coqhull::spget ( vector < vector < vec > > points, int spa, int spb ) {
	int a_start,b_start;
	// loop for all atoms
	for ( int ia=0; ia<(int)co->stypeindex[spa].size(); ia++ ) {
		int a=co->stypeindex[spa][ia];a_start=co->mystart(a);
		for ( int i=a_start; i<(a_start+co->item[a]); i++ ) {
			if ( 2<(int)points[i].size() ) spget_one(points[i],i,spa,spb);
		}
	}
}