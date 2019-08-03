#include "corder.h"

// destructor
gofr::~gofr() {
	free1d < double > (rhob);
	free1d < double > (rhosb);
	free1d < double > (gr);
	free2d < double > (gabr);
	free2d < double > (gsabr);
	free1d < double > (rhor);
	free2d < double > (rhoabr);
	free2d < double > (rhosabr);
}

void gofr::init ( void ) {

	// temporary variable
	nr = (int)((co->r1-co->r0)/co->dr);

	// memory allocation
	rho = 0;
	rhob = alloc1d < double > (co->ntypes);
	rhosb = alloc1d < double > (co->nstypes);
	gr = alloc1d < double > (nr);
	gabr = alloc2d < double > (co->ntypes*co->ntypes,nr);
	gsabr = alloc2d < double > (co->nstypes*co->nstypes,nr);
	rhor = alloc1d < double > (nr);
	rhoabr = alloc2d < double > (co->ntypes*co->ntypes,nr);
	rhosabr = alloc2d < double > (co->nstypes*co->nstypes,nr);

}

void gofr::update ( int rank, int nmpi ) {

	// temporary variable
	int ir=0,a=0,b=0;
	int local_start,local_size;
	double r,vol,fr,fabr,fsabr;
	int *ngr,*ngr0;
	int **ngabr,**ngabr0;
	int **ngsabr,**ngsabr0;
	vec T(3);

	// mpi parameter
	mpi_local_part(&local_start,&local_size,co->nions,rank,nmpi);

	// volume
	vol = co->a1*(co->a2.cross(co->a3));
	rho = co->update(rho,(double)co->nions/vol);
	for ( b=0; b<(co->ntypes); b++ ) {
		rhob[b] = co->update(rhob[b],(double)co->item[b]/vol);
	}
	for ( int spb=0; spb<(co->nstypes); spb++ ) {
		rhosb[spb] = co->update(rhosb[spb],(double)co->sitem[spb]/vol);
	}

	// memory allocation
	ngr = alloc1d < int > (nr);
	ngabr = alloc2d < int > (co->ntypes*co->ntypes,nr);
	ngsabr = alloc2d < int > (co->nstypes*co->nstypes,nr);
	ngr0 = alloc1d < int > (nr);
	ngabr0 = alloc2d < int > (co->ntypes*co->ntypes,nr);
	ngsabr0 = alloc2d < int > (co->nstypes*co->nstypes,nr);

	// calc gofr
	for ( int i=local_start; i<local_start+local_size; i++ ) {
	for ( int j=0; j<(co->nions); j++ ) {
		a = co->mytype(i);
		b = co->mytype(j);
		if ( i != j ) {
			for ( int t1=-(co->trans); t1<=(co->trans); t1++ ) {
			for ( int t2=-(co->trans); t2<=(co->trans); t2++ ) {
			for ( int t3=-(co->trans); t3<=(co->trans); t3++ ) {
				T = (double)t1*co->a1+(double)t2*co->a2+(double)t3*co->a3;
				r = (co->config[i]-co->config[j]+T).norm();
				ir = (int)((r-co->r0)/co->dr);
				if ( ir < nr ) {
					ngr[ir] ++;
					ngabr[a*co->ntypes+b][ir] ++;
					for ( int spa=0; spa<(int)co->stypelist[a].size(); spa++ ) {
					for ( int spb=0; spb<(int)co->stypelist[b].size(); spb++ ) {
						ngsabr[co->stypelist[a][spa]*co->nstypes+co->stypelist[b][spb]][ir] ++;
					}
					}
				}
			}
			}
			}
		}
	}
	}
	MPI_Barrier(MPI_COMM_WORLD);

	// data gathering
	MPI_Allreduce(&(ngr[0]),&(ngr0[0]),nr,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(&(ngabr[0][0]),&(ngabr0[0][0]),co->ntypes*co->ntypes*nr,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(&(ngsabr[0][0]),&(ngsabr0[0][0]),co->nstypes*co->nstypes*nr,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	// temporary check: multi-dimensional array that same as ngabr
//	MPI_Datatype stype;
//	MPI_Type_vector(co->ntypes*co->ntypes,1,nr,MPI_INT,&stype);
//	MPI_Type_commit(&stype);
//	for ( ir=0; ir<nr; ir++ ) {
//		MPI_Reduce(&(ngabr[0][ir]),&(ngabr0[0][ir]),1,stype,MPI_SUM,0,MPI_COMM_WORLD);
//	}

	// free
	free1d < int > (ngr);
	free2d < int > (ngabr);
	free2d < int > (ngsabr);

	// update g(r) and rho*g(r)
	// update gab(r) and rhob*gab(r)
	gr[0] = 0;
	rhor[0] = 0;
	for ( int iab=0; iab<(co->ntypes*co->ntypes); iab++ ) {
		gabr[iab][0] = 0;
		rhoabr[iab][0] = 0;
	}
	for ( int isab=0; isab<(co->nstypes*co->nstypes); isab++ ) {
		gsabr[isab][0] = 0;
		rhosabr[isab][0] = 0;
	}
	for ( ir=1; ir<nr; ir++ ) {
		r = (double)ir*co->dr + co->r0;
		fr = 4.0*pi*r*r*co->dr*rho*(double)co->nions;
		gr[ir] = co->update(gr[ir],(double)ngr0[ir]/fr);
		rhor[ir] = co->update(rhor[ir],rho*(double)ngr0[ir]/fr);
		for ( a=0; a<(co->ntypes); a++ ) {
		for ( b=0; b<(co->ntypes); b++ ) {
			fabr = 4.0*pi*r*r*co->dr*rhob[b]*(double)co->item[a];
			gabr[a*co->ntypes+b][ir] = co->update(gabr[a*co->ntypes+b][ir],(double)ngabr0[a*co->ntypes+b][ir]/fabr);
			rhoabr[a*co->ntypes+b][ir] = co->update(rhoabr[a*co->ntypes+b][ir],rhob[b]*(double)ngabr0[a*co->ntypes+b][ir]/fabr);
		}
		}
		for ( int spa=0; spa<(co->nstypes); spa++ ) {
		for ( int spb=0; spb<(co->nstypes); spb++ ) {
			fsabr = 4.0*pi*r*r*co->dr*rhosb[spb]*(double)co->sitem[spa];
			gsabr[spa*co->nstypes+spb][ir] = co->update(gsabr[spa*co->nstypes+spb][ir],(double)ngsabr0[spa*co->nstypes+spb][ir]/fsabr);
			rhosabr[spa*co->nstypes+spb][ir] = co->update(rhosabr[spa*co->nstypes+spb][ir],rhosb[spb]*(double)ngsabr0[spa*co->nstypes+spb][ir]/fsabr);
		}
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);

	// free
	free1d < int > (ngr0);
	free2d < int > (ngabr0);
	free2d < int > (ngsabr0);

}

void gofr::write ( string filename ) {
	double r,cr=0,*cabr,*csabr;
	ofstream ofs,ofsab,ofssab;

	// memory allocation
	cabr = alloc1d < double > (co->ntypes*co->ntypes);
	csabr = alloc1d < double > (co->nstypes*co->nstypes);

	// write out
	ofs.open((filename+"-tot.dat").c_str(),ios::out);
	ofsab.open((filename+".dat").c_str(),ios::out);
	ofssab.open((filename+"-spec.dat").c_str(),ios::out);
	ofs << scientific << setprecision(5);
	ofsab << scientific << setprecision(5);
	ofssab << scientific << setprecision(5);
	// g(r)
	for ( int ir=0; ir<nr; ir++ ) {
		r = (double)ir*co->dr + co->r0;
		ofs << r << "\t" << gr[ir] << endl;
		ofsab << r; ofssab << r;
		for ( int iab=0; iab<(co->ntypes*co->ntypes); iab++ )
			ofsab << "\t" << gabr[iab][ir];
		for ( int isab=0; isab<(co->nstypes*co->nstypes); isab++ )
			ofssab << "\t" << gsabr[isab][ir];
		ofsab << endl; ofssab << endl;
	}
	ofs << endl; ofsab << endl; ofssab << endl;
	ofs << endl; ofsab << endl; ofssab << endl;
	// c(r)
	for ( int ir=0; ir<nr; ir++ ) {
		r = (double)ir*co->dr + co->r0;

		// integrate 4*pi*r*r*rho*g(r) >> coordination number c(r)
		if ( ir == 0 ) {
			cr = rhor[0]*co->dr;
			for ( int iab=0; iab<(co->ntypes*co->ntypes); iab++ )
				cabr[iab] = r*r*rhoabr[iab][0];
			for ( int isab=0; isab<(co->nstypes*co->nstypes); isab++ )
				csabr[isab] = r*r*rhosabr[isab][0];
		} else if ( ir == 1 ) {
			cr = 0.5*(rhor[0]+rhor[1])*co->dr;
			for ( int iab=0; iab<(co->ntypes*co->ntypes); iab++ )
				cabr[iab] = 0.5*r*r*(rhoabr[iab][0]+rhoabr[iab][1]);
			for ( int isab=0; isab<(co->nstypes*co->nstypes); isab++ )
				csabr[isab] = 0.5*r*r*(rhosabr[isab][0]+rhosabr[isab][1]);
		} else {
			cr += 0.5*(rhor[ir-1]+rhor[ir])*co->dr;
			for ( int iab=0; iab<(co->ntypes*co->ntypes); iab++ )
				cabr[iab] += 0.5*r*r*(rhoabr[iab][ir-1]+rhoabr[iab][ir]);
			for ( int isab=0; isab<(co->nstypes*co->nstypes); isab++ )
				csabr[isab] += 0.5*r*r*(rhosabr[isab][ir-1]+rhosabr[isab][ir]);
		}
		ofs << r << "\t" << 4.0*pi*cr*co->dr << endl;
		ofsab << r; ofssab << r;
		for ( int iab=0; iab<(co->ntypes*co->ntypes); iab++ )
			ofsab << "\t" << 4.0*pi*cabr[iab]*co->dr;
		for ( int isab=0; isab<(co->nstypes*co->nstypes); isab++ )
			ofssab << "\t" << 4.0*pi*csabr[isab]*co->dr;
		ofsab << endl; ofssab << endl;
	}
	ofs.close();
	ofsab.close();
	ofssab.close();

	// free
	free1d < double > (cabr);
	free1d < double > (csabr);
}
