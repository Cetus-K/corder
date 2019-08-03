#include "corder.h"

// destructor
bofree::~bofree() {
	free1d < double > (il);
	free2d < double > (jll);
	free2d < double > (ilab);
	free3d < double > (jllab);
	free2d < double > (ilsab);
	free3d < double > (jllsab);
}

void bofree::init ( void ) {

	// temporary variable
	nl = co->lcut+1;

	// memory allocation
	// free energy density coefficient
	// fl = rl*il+sum_{ l'=0...l }[ wl'l*jl'l ]
	il = alloc1d < double > (nl);
	jll = alloc2d < double > (nl,nl);
	ilab = alloc2d < double > (nl,co->ntypes*co->ntypes);
	jllab = alloc3d < double > (nl,nl,co->ntypes*co->ntypes);
	ilsab = alloc2d < double > (nl,co->nstypes*co->nstypes);
	jllsab = alloc3d < double > (nl,nl,co->nstypes*co->nstypes);

}

void bofree::update ( int rank, int nmpi ) {

	// temporary variable
	int a,b,a_start,iab,isab;			// voro++
	double *s,**sb,**ssb;
	int local_start,local_size;		// mpi: temporary
	int l,m,m1,m2,m3;			// bofree: temporary
	double r,theta,phi;
	complex < double > Ylm(0,0),wl_buf(0,0);
	complex < double > ***qlm,****qlmb,****qlmsb;	// boofree: saving
	double *flcf,**fllcf,**flabcf,***fllabcf,**flsabcf,***fllsabcf;

	// total/partial surface area
	s = alloc1d < double > (co->nions);
	sb = alloc2d < double > (co->nions,co->ntypes);
	ssb = alloc2d < double > (co->nions,co->nstypes);
	for ( int i=0; i<(co->nions); i++ ) {
	for ( int ij=0; ij<(int)cv->surf[i].size(); ij++ ) {
		s[i] += cv->surf[i][ij];
	}
	}
	// partial surface area
	for ( int a=0; a<(co->ntypes); a++ ) { a_start = co->mystart(a);
	for ( int b=0; b<(co->ntypes); b++ ) { iab = a*co->ntypes+b;
	for ( int i=a_start; i<(a_start+co->item[a]); i++ ) {
	for ( int ij=0; ij<(int)cv->surfab[iab][i].size(); ij++ ) {
		// for na
		if ( a == b ) { 
			sb[i][b] += cv->surfab[iab][i][ij];
		// for nb of na+nb
		} else {
			if ( b == co->mytype(cv->nnidab[iab][i][ij]) ) {
				sb[i][b] += cv->surfab[iab][i][ij];
			}
		}
	}
	}
	}
	}
	// specified partial surface area
	for ( int spa=0; spa<(co->nstypes); spa++ ) {
	for ( int spb=0; spb<(co->nstypes); spb++ ) { isab = spa*co->nstypes+spb;
	for ( int ia=0; ia<(int)co->stypeindex[spa].size(); ia++ ) {
		int a=co->stypeindex[spa][ia];a_start=co->mystart(a);
		for ( int i=a_start; i<(a_start+co->item[a]); i++ ) {
		for ( int ij=0; ij<(int)cv->surfsab[isab][i].size(); ij++ ) {
			// for nsa
			if ( spa == spb ) { 
				ssb[i][spb] += cv->surfsab[isab][i][ij];
			// for nsb of nsa+nsb
			} else {
				int b=co->mytype(cv->nnidsab[isab][i][ij]);
				vector<int>::iterator ib=find(co->stypelist[b].begin(),co->stypelist[b].end(),spb);
				if ( ib != co->stypelist[b].end() )
					ssb[i][spb] += cv->surfsab[isab][i][ij];
			}
		}
		}
	}
	}
	}

	// local free energy density coefficient
	// fl = rl*il+sum_{ l'=0...l }[ wl'l*jl'l ]
	flcf = alloc1d < double > (nl);
	fllcf = alloc2d < double > (nl,nl);
	flabcf = alloc2d < double > (nl,co->ntypes*co->ntypes);
	fllabcf = alloc3d < double > (nl,nl,co->ntypes*co->ntypes);
	flsabcf = alloc2d < double > (nl,co->nstypes*co->nstypes);
	fllsabcf = alloc3d < double > (nl,nl,co->nstypes*co->nstypes);

	// jag array
	qlm = new complex < double >** [nl];
	qlmb = new complex < double >*** [nl];
	qlmsb = new complex < double >*** [nl];
	for ( l=0; l<nl; l++ ) {
		qlm[l] = new complex < double >* [2*l+1];
		qlmb[l] = new complex < double >** [2*l+1];
		qlmsb[l] = new complex < double >** [2*l+1];
		for ( m=-l; m<=l; m++ ) {
			qlm[l][m+l] = alloc1d < complex < double > > (co->nions);
			qlmb[l][m+l] = alloc2d < complex < double > > (co->nions,co->ntypes);
			qlmsb[l][m+l] = alloc2d < complex < double > > (co->nions,co->nstypes);
		}
	}

	// calc qlm,qlmb,qlmsb
	for ( l=0; l<nl; l++ ) {
		// total information
		for ( int i=0; i<(co->nions); i++ ) {
		for ( m=-l; m<=l; m++ ) {
			if ( !cv->nnid[i].empty() ) {
				for ( int ij=0; ij<(int)cv->nnid[i].size(); ij++ ) {
					cart2shpr(&r,&theta,&phi,cv->nv[i][ij](0),cv->nv[i][ij](1),cv->nv[i][ij](2));
					Ylm = spherical_Ylm(l,m,theta,phi);
					// qlm
					qlm[l][m+l][i] += cv->surf[i][ij]*conj(Ylm);
				}
				qlm[l][m+l][i] /= s[i];
			}
		}
		}
		// partial information
		for ( int a=0; a<(co->ntypes); a++ ) { a_start = co->mystart(a);
		for ( int b=0; b<(co->ntypes); b++ ) { iab = a*co->ntypes+b;
		for ( int i=a_start; i<(a_start+co->item[a]); i++ ) {
		for ( m=-l; m<=l; m++ ) {
			if ( !cv->nnid[i].empty() ) {
				for ( int ij=0; ij<(int)cv->nnidab[iab][i].size(); ij++ ) {
					cart2shpr(&r,&theta,&phi,cv->nvab[iab][i][ij](0),cv->nvab[iab][i][ij](1),cv->nvab[iab][i][ij](2));
					Ylm = spherical_Ylm(l,m,theta,phi);
					// qlmb
					// for na
					if ( a == b ) {
						qlmb[l][m+l][i][b] += cv->surfab[iab][i][ij]*conj(Ylm);
					// for nb of na+nb
					} else {
						if ( b == co->mytype(cv->nnidab[iab][i][ij]) ) {
							qlmb[l][m+l][i][b] += cv->surfab[iab][i][ij]*conj(Ylm);
						}
					}
				}
				qlmb[l][m+l][i][b] /= sb[i][b];
			}
		}
		}
		}
		}
		// specified partial information
		for ( int spa=0; spa<(co->nstypes); spa++ ) {
		for ( int spb=0; spb<(co->nstypes); spb++ ) { isab = spa*co->nstypes+spb;
		for ( int ia=0; ia<(int)co->stypeindex[spa].size(); ia++ ) {
			int a=co->stypeindex[spa][ia];a_start=co->mystart(a);
			for ( int i=a_start; i<(a_start+co->item[a]); i++ ) {
			for ( m=-l; m<=l; m++ ) {
				if ( !cv->nnidsab[isab][i].empty() ) {
					for ( int ij=0; ij<(int)cv->nnidsab[isab][i].size(); ij++ ) {
						cart2shpr(&r,&theta,&phi,cv->nvsab[isab][i][ij](0),cv->nvsab[isab][i][ij](1),cv->nvsab[isab][i][ij](2));
						Ylm = spherical_Ylm(l,m,theta,phi);
						// qlmb
						// for na
						if ( spa == spb ) {
							qlmsb[l][m+l][i][spb] += cv->surfsab[isab][i][ij]*conj(Ylm);
						// for nb of na+nb
						} else {
							int b=co->mytype(cv->nnidsab[isab][i][ij]);
							vector<int>::iterator ib=find(co->stypelist[b].begin(),co->stypelist[b].end(),spb);
							if ( ib != co->stypelist[b].end() )
								qlmsb[l][m+l][i][spb] += cv->surfsab[isab][i][ij]*conj(Ylm);
						}
					}
					qlmsb[l][m+l][i][spb] /= ssb[i][spb];
				}
			}
			}
		}
		}
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);

	// mpi parameter
	mpi_local_part(&local_start,&local_size,co->nions,rank,nmpi);

	// calculate bond-orientational Landau free energy coefficient
	for ( l=0; l<nl; l++ ) {
	for ( int i=local_start; i<local_start+local_size; i++ ) {
		a = co->mytype(i);

		// total bofree coefficient
		for ( m=-l; m<=l; m++ ) flcf[l] += norm(qlm[l][m+l][i]);
		for( int lj=0; lj<nl; lj++ ){
			wl_buf = complex < double > (0,0);
			for ( m1=-l; m1<=l; m1++ ) {
			for ( m2=max(-l,-(lj+m1)); m2<=min(l,lj-m1); m2++ ) {
				m3 = -m1-m2;
				wl_buf += Wigner_3j(l,l,lj,m1,m2,m3)*qlm[l][m1+l][i]*qlm[l][m2+l][i]*qlm[lj][m3+lj][i];
			}
			}
			fllcf[l][lj] += wl_buf.real();
		}

		// partial bofree coefficient
		for ( b=0; b<(co->ntypes); b++ ) {
			for ( m=-l; m<=l; m++ ) flabcf[l][a*co->ntypes+b] += norm(qlmb[l][m+l][i][b]);
			for ( int lj=0; lj<nl; lj++ ) {
				wl_buf = complex < double > (0,0);
				for ( m1=-l; m1<=l; m1++ ) {
				for ( m2=max(-l,-(lj+m1)); m2<=min(l,lj-m1); m2++ ) {
					m3 = -m1-m2;
					wl_buf += Wigner_3j(l,l,lj,m1,m2,m3)*qlmb[l][m1+l][i][b]*qlmb[l][m2+l][i][b]*qlmb[lj][m3+lj][i][b];
				}
				}
				fllabcf[l][lj][a*co->ntypes+b] += wl_buf.real();
			}
		}

		// specified partial bofree coefficient
		for ( int spa=0; spa<(int)co->stypelist[a].size(); spa++ ) {
		for ( int spb=0; spb<(co->nstypes); spb++ ) {
			for ( m=-l; m<=l; m++ ) flsabcf[l][co->stypelist[a][spa]*co->nstypes+spb] += norm(qlmsb[l][m+l][i][spb]);
			for ( int lj=0; lj<nl; lj++ ) {
				wl_buf = complex < double > (0,0);
				for ( m1=-l; m1<=l; m1++ ) {
				for ( m2=max(-l,-(lj+m1)); m2<=min(l,lj-m1); m2++ ) {
					m3 = -m1-m2;
					wl_buf += Wigner_3j(l,l,lj,m1,m2,m3)*qlmsb[l][m1+l][i][spb]*qlmb[l][m2+l][i][spb]*qlmsb[lj][m3+lj][i][spb];
				}
				}
				fllsabcf[l][lj][co->stypelist[a][spa]*co->nstypes+spb] += wl_buf.real();
			}
		}
		}

	}
	}
	MPI_Barrier(MPI_COMM_WORLD);

	// data gathering
	MPI_Allreduce(&(flcf[0]),&(il[0]),nl,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(&(fllcf[0][0]),&(jll[0][0]),nl*nl,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(&(flabcf[0][0]),&(ilab[0][0]),nl*co->ntypes*co->ntypes,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(&(fllabcf[0][0][0]),&(jllab[0][0][0]),nl*nl*co->ntypes*co->ntypes,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(&(flsabcf[0][0]),&(ilsab[0][0]),nl*co->nstypes*co->nstypes,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(&(fllsabcf[0][0][0]),&(jllsab[0][0][0]),nl*nl*co->nstypes*co->nstypes,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	// free
	free1d < double > (s);
	free2d < double > (sb);
	free2d < double > (ssb);
	free1d < double > (flcf);
	free2d < double > (fllcf);
	free2d < double > (flabcf);
	free3d < double > (fllabcf);
	free2d < double > (flsabcf);
	free3d < double > (fllsabcf);
	for ( l=0; l<nl; l++ ) {
		for ( m=-l; m<=l; m++ ) {
			free1d < complex < double > > (qlm[l][m+l]);
			free2d < complex < double > > (qlmb[l][m+l]);
			free2d < complex < double > > (qlmsb[l][m+l]);
		}
		delete [] qlm[l]; qlm[l] = NULL;
		delete [] qlmb[l]; qlmb[l] = NULL;
		delete [] qlmsb[l]; qlmsb[l] = NULL;
	}
	delete [] qlm; qlm = NULL;
	delete [] qlmb; qlmb = NULL;
	delete [] qlmsb; qlmsb = NULL;

}

void bofree::write ( string filename ) {
	double t=0;
	ofstream ofsi,ofsj,ofsiab,ofsjab,ofsisab,ofsjsab;

	// write out
	if ( co->itr == co->f0 ) {
		ofsi.open((filename+"-i-tot.dat").c_str(),ios::out);
		ofsj.open((filename+"-j-tot.dat").c_str(),ios::out);
		ofsiab.open((filename+"-i.dat").c_str(),ios::out);
		ofsjab.open((filename+"-j.dat").c_str(),ios::out);
		ofsisab.open((filename+"-i-spec.dat").c_str(),ios::out);
		ofsjsab.open((filename+"-j-spec.dat").c_str(),ios::out);
	} else {
		ofsi.open((filename+"-i-tot.dat").c_str(),ios::app);
		ofsj.open((filename+"-j-tot.dat").c_str(),ios::app);
		ofsiab.open((filename+"-i.dat").c_str(),ios::app);
		ofsjab.open((filename+"-j.dat").c_str(),ios::app);
		ofsisab.open((filename+"-i-spec.dat").c_str(),ios::app);
		ofsjsab.open((filename+"-j-spec.dat").c_str(),ios::app);
	}
	ofsi << scientific << setprecision(5);
	ofsj << scientific << setprecision(5);
	ofsiab << scientific << setprecision(5);
	ofsjab << scientific << setprecision(5);
	ofsisab << scientific << setprecision(5);
	ofsjsab << scientific << setprecision(5);

	if ( nl ) {
		t = (double)co->itr*co->dt;

		// il,jll
		ofsi << t;
		for ( int l=0; l<nl; l++ ) {
			// il
			ofsi << "\t" << il[l];
			// jll
			ofsj << jll[l][0];
			for ( int lj=1; lj<nl; lj++ ) ofsj << "\t" << jll[l][lj];
			ofsj << endl;
		}
		ofsi << endl;
		ofsj << endl << endl;

		// ilab,jllab
		for ( int iab=0; iab<(co->ntypes*co->ntypes); iab++ ) {
			// ilab
			ofsiab << ilab[0][iab];
			for ( int l=1; l<nl; l++ ) ofsiab << "\t" << ilab[l][iab];
			ofsiab << endl << endl;
			// jllab
			for ( int l=0; l<nl; l++ ) {
				ofsjab << jllab[l][0][iab];
				for ( int lj=1; lj<nl; lj++ ) ofsjab << "\t" << jllab[l][lj][iab];
				ofsjab << endl;
			}
			ofsjab << endl << endl;
		}
		ofsiab << endl;
		ofsjab << endl << endl;

		// ilsab,jllsab
		for ( int isab=0; isab<(co->nstypes*co->nstypes); isab++ ) {
			// ilsab
			ofsisab << ilsab[0][isab];
			for ( int l=1; l<nl; l++ ) ofsisab << "\t" << ilsab[l][isab];
			ofsisab << endl << endl;
			// jllsab
			for ( int l=0; l<nl; l++ ) {
				ofsjsab << jllsab[l][0][isab];
				for ( int lj=1; lj<nl; lj++ ) ofsjsab << "\t" << jllsab[l][lj][isab];
				ofsjsab << endl;
			}
			ofsjsab << endl << endl;
		}
		ofsisab << endl;
		ofsjsab << endl << endl;

	}
	ofsi.close();
	ofsj.close();
	ofsiab.close();
	ofsjab.close();
	ofsisab.close();
	ofsjsab.close();

}
