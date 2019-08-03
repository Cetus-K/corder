#include "corder.h"

// set
void dfc::set ( corder &cdr ) {
	co = &cdr;
	if ( co->itr == co->f0 ) {
		for ( int i=0; i<(co->nions); i++ ) {
			config_init[i] = co->config[i];
			config_prev[i] = co->config[i];
		}
	}
}

// destructor
dfc::~dfc() {
	free1d < double > (msda);
	free1d < double > (msdsa);
	vector < vec > ().swap(config_init);
	vector < vec > ().swap(config_prev);
	vector < vec > ().swap(veloca);
	vector < vec > ().swap(velocsa);
}

void dfc::init ( void ) {

	// memory allocation
	msd = 0;
	msda = alloc1d < double > (co->ntypes);
	msdsa = alloc1d < double > (co->nstypes);
	config_init.resize(co->nions);
	config_prev.resize(co->nions);
	veloc.resize(3);
	veloca.resize(co->ntypes);
	velocsa.resize(co->nstypes);
	for ( int i=0; i<(co->nions); i++ ) {
		config_init[i].resize(3);
		config_prev[i].resize(3);
	}
	for ( int a=0; a<(co->ntypes); a++ )
		veloca[a].resize(3);
	for ( int spa=0; spa<(co->nstypes); spa++ )
		velocsa[spa].resize(3);

}

void dfc::update ( void ) {

	// temporary variable
	int a;
	double msdj=0,*msdaj,*msdsaj;
	double r;
	vec diff(3),crct(3),tol(3);	// check whether applied pb
	mat invmat(3,3);
	vec vlc(3);			// velocity
	vector < vec > vlca,vlcsa;

	// memory allocation
	msdaj = alloc1d < double > (co->ntypes);
	msdsaj = alloc1d < double > (co->nstypes);
	vlca.resize(co->ntypes);
	vlcsa.resize(co->nstypes);
	for ( a=0; a<(co->ntypes); a++ ) vlca[a].resize(3);
	for ( int spa=0; spa<(co->nstypes); spa++ ) vlcsa[spa].resize(3);

	// inverse matrix
	invmat = co->latmat.inv3d();

	// tolerrance to pb correction
	tol(0) = 0.7;
	tol(1) = 0.7;
	tol(2) = 0.7;

	if ( co->itr == co->f0 ) {	
		vlc.zero();
		for ( a=0; a<(co->ntypes); a++ ) vlca[a].zero();
		for ( int spa=0; spa<(co->nstypes); spa++ ) vlcsa[spa].zero();
	} else {

		// calc msd and msd*t
		// simultaneously correct the configuration
		for ( int i=0; i<(co->nions); i++ ) {
			a = co->mytype(i);
	
			// check pb
			diff = invmat*(co->config[i]-config_prev[i]);
			for ( int j=0; j<3; j++ ) {
				if ( diff(j) < -tol(j) ) crct(j) = 1.0;
				else if ( diff(j) > tol(j) ) crct(j) = -1.0;
				else crct(j) = 0;
			}
			crct = co->latmat*crct;
	
			// calc velocity
			vlc += (co->config[i]+crct-config_prev[i])*(1.0/co->dt);
			vlca[a] += (co->config[i]+crct-config_prev[i])*(1.0/co->dt);
			for ( int spa=0; spa<(int)co->stypelist[a].size(); spa++ )
				vlcsa[co->stypelist[a][spa]] += (co->config[i]+crct-config_prev[i])*(1.0/co->dt);

			// calc msd and msd*t
			r = (co->config[i]+crct-config_init[i]).norm();
			msdj += pow(r,2.0);
			msdaj[a] += pow(r,2.0);
			for ( int spa=0; spa<(int)co->stypelist[a].size(); spa++ )
				msdsaj[co->stypelist[a][spa]] += pow(r,2.0);
	
			// update previous
			config_prev[i] = co->config[i] + crct;
	
		}
	
		// update velocity
		for ( int i=0; i<3; i++ ) {
			veloc(i) = co->update(veloc(i),vlc(i)/(double)co->nions);
			for ( a=0; a<(co->ntypes); a++ )
				veloca[a](i) = co->update(veloca[a](i),vlca[a](i)/(double)co->item[a]);
			for ( int spa=0; spa<(co->nstypes); spa++ )
				velocsa[spa](i) = co->update(velocsa[spa](i),vlcsa[spa](i)/(double)co->sitem[spa]);
		}

		// update msd
		msd = msdj/(double)co->nions;
		for ( a=0; a<(co->ntypes); a++ )
			msda[a] = msdaj[a]/(double)co->item[a];
		for ( int spa=0; spa<(co->nstypes); spa++ )
			msdsa[spa] = msdsaj[spa]/(double)co->sitem[spa];

	}

	// free
	free1d < double > (msdaj);
	free1d < double > (msdsaj);
	vector < vec > ().swap(vlca);
	vector < vec > ().swap(vlcsa);

}

void dfc::write ( string filename ) {
	double t;
	ofstream ofs,ofsa,ofssa;

	// write out
	if ( co->itr == co->f0 ) {
		ofs.open((filename+"-tot.dat").c_str(),ios::out);
		ofsa.open((filename+".dat").c_str(),ios::out);
		ofssa.open((filename+"-spec.dat").c_str(),ios::out);
	} else {
		ofs.open((filename+"-tot.dat").c_str(),ios::app);
		ofsa.open((filename+".dat").c_str(),ios::app);
		ofssa.open((filename+"-spec.dat").c_str(),ios::app);
	}
	ofs << scientific << setprecision(5);
	ofsa << scientific << setprecision(5);
	ofssa << scientific << setprecision(5);

	t = (double)co->itr*co->dt;
	ofs << t << "\t" << msd << "\t" << t << veloc.norm() << endl;
	ofsa << t; ofssa << t;
	for ( int a=0; a<(co->ntypes); a++ ) ofsa << "\t" << msda[a];
	for ( int spa=0; spa<(co->nstypes); spa++ ) ofssa << "\t" << msdsa[spa];
	ofsa << "\t" << t; ofssa << "\t" << t;
	for ( int a=0; a<(co->ntypes); a++ ) ofsa << "\t" << veloca[a].norm();
	for ( int spa=0; spa<(co->nstypes); spa++ ) ofssa << "\t" << velocsa[spa].norm();
	ofsa << endl; ofssa << endl;

	ofs.close();
	ofsa.close();
	ofssa.close();

}
