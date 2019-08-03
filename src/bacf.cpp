#include "corder.h"

// destructor
bacf::~bacf() {
	free2d < double > (gglr);
	free3d < double > (gglababr);
	free3d < double > (gglsababr);
}

void bacf::init ( void ) {

	// temporary variable
	nl = (int)co->l.size();
	nr = (int)((co->r1-co->r0)/co->dr);

	// memory allocation
	gglr = alloc2d < double > (nl,nr);
	gglababr = alloc3d < double > (nl,co->ntypes*co->ntypes*co->ntypes*co->ntypes,nr);
	gglsababr = alloc3d < double > (nl,co->nstypes*co->nstypes*co->nstypes*co->nstypes,nr);

}

void bacf::update ( int rank, int nmpi ) {

	// temporary variable
	int a,b,iab,isab,a_start=0;			// voro++
	vec pair(3),bmp(3),T(3);
	mat nwe(3,3);
	vector < vector < int > > iabs,isabs;
	vector < double > s;
	vector < vector < double > > sab,ssab;
	vector < vec > bond,dirc;
	vector < vector < vec > > bondab,dircab;
	vector < vector < vec > > bondsab,dircsab;
	string bmpstr;				// for utilizing map function
	stringstream ss;
	map < string,bool > duple;
	int local_start,local_size;			// mpi: temporary
	int ir,l;					// bacf: temporary
	double r,ra,rb,thetaa,thetab,phia,phib,cosw;
	double *sr,**sababr,**ssababr,*sr0,**sababr0,**ssababr0;	// bacf: saving
	double **ssplr,***ssplababr,***ssplsababr,**ssplr0,***ssplababr0,***ssplsababr0;

	// memory allocation
	// middle point of a bond and its direction
	// total/partial surface area
	sr = alloc1d < double > (nr);
	sababr = alloc2d < double > (co->ntypes*co->ntypes*co->ntypes*co->ntypes,nr);
	ssababr = alloc2d < double > (co->nstypes*co->nstypes*co->nstypes*co->nstypes,nr);
	ssplr = alloc2d < double > (nl,nr);
	ssplababr = alloc3d < double > (nl,co->ntypes*co->ntypes*co->ntypes*co->ntypes,nr);
	ssplsababr = alloc3d < double > (nl,co->nstypes*co->nstypes*co->nstypes*co->nstypes,nr);
	iabs.resize(co->ntypes*co->ntypes);
	isabs.resize(co->nstypes*co->nstypes);
	sab.resize(co->ntypes*co->ntypes);
	ssab.resize(co->nstypes*co->nstypes);
	bondab.resize(co->ntypes*co->ntypes);
	bondsab.resize(co->nstypes*co->nstypes);
	dircab.resize(co->ntypes*co->ntypes);
	dircsab.resize(co->nstypes*co->nstypes);

	// save bonds which are not duplicate to each other points
	ss.precision(5);
	// total information
	for ( int i=0; i<(co->nions); i++ ) {
		if ( !cv->nnid[i].empty() ) {
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
		
				// save middle point of the bond not duplicately
				// and its direction; normal vector
				bmp = 0.5*(cv->config[i]+pair);
				// erase uncertainty
				for ( int ijk=0; ijk<3; ijk++ ) {
					if ( pow(bmp(ijk),2.0) < under_capacity )
						bmp(ijk) = 0.0;
				}
				ss << bmp(0) << "," << bmp(1) << "," << bmp(2);
				bmpstr = ss.str();
				ss.str(""); ss.clear();
				if ( !duple[bmpstr] ) {
					duple[bmpstr] = true;
					s.push_back(cv->surf[i][ij]);
					bond.push_back(bmp);
					dirc.push_back(cv->nv[i][ij]);
				}
		
			}
		}
	}
	// free boolean duple array
	duple.clear();

	// partial information
	for ( int a=0; a<(co->ntypes); a++ ) { a_start = co->mystart(a);
	for ( int b=0; b<(co->ntypes); b++ ) { iab = a*co->ntypes+b;
		for ( int i=a_start; i<(a_start+co->item[a]); i++ ) {
			if ( !cv->nnidab[iab][i].empty() ) {
				for ( int ij=0; ij<(int)cv->nnidab[iab][i].size(); ij++ ) {
					// calc r(i)'s pair position of r(nnid[i][ij])
					// from normal vector of voronoi surface
					// which direction is r(i) to r(nnid[i][ij])
					for ( int j=0; j<3; j++ ) {
						nwe(0,j) = cv->egab[iab][i][ij][0](j);
						nwe(1,j) = cv->egab[iab][i][ij][1](j);
						nwe(2,j) = cv->nvab[iab][i][ij](j);
					}
					pair(0) = cv->config[i]*cv->egab[iab][i][ij][0];
					pair(1) = cv->config[i]*cv->egab[iab][i][ij][1];
					pair(2) = (2.0*cv->egab[iab][i][ij][2]-cv->config[i])*cv->nvab[iab][i][ij];
					pair = nwe.inv3d()*pair;
	
					// save middle point of the bond not duplicately
					// and its direction; normal vector
					if ( a == b ) {
						bmp = 0.5*(cv->config[i]+pair);
						// erase uncertainty
						for ( int ijk=0; ijk<3; ijk++ ) {
							if ( pow(bmp(ijk),2.0) < under_capacity )
								bmp(ijk) = 0.0;
						}
						ss << bmp(0) << "," << bmp(1) << "," << bmp(2);
						bmpstr = ss.str();
						ss.str(""); ss.clear();
						if ( !duple[bmpstr] ) {
							duple[bmpstr] = true;
							iabs[iab].push_back(iab);
							sab[iab].push_back(cv->surfab[iab][i][ij]);
							bondab[iab].push_back(bmp);
							dircab[iab].push_back(cv->nvab[iab][i][ij]);
						}
					} else {
						if ( b == co->mytype(cv->nnidab[iab][i][ij]) ) {
							bmp = 0.5*(cv->config[i]+pair);
							// erase uncertainty
							for ( int ijk=0; ijk<3; ijk++ ) {
								if ( pow(bmp(ijk),2.0) < under_capacity )
									bmp(ijk) = 0.0;
							}
							ss << bmp(0) << "," << bmp(1) << "," << bmp(2);
							bmpstr = ss.str();
							ss.str(""); ss.clear();
							if ( !duple[bmpstr] ) {
								duple[bmpstr] = true;
								iabs[iab].push_back(iab);
								sab[iab].push_back(cv->surfab[iab][i][ij]);
								bondab[iab].push_back(bmp);
								dircab[iab].push_back(cv->nvab[iab][i][ij]);
							}
						}
					}
	
				}
			}
		}
		// free boolean duple array
		duple.clear();
	}
	}

	// specified partial information
	for ( int sa=0; sa<(co->nstypes); sa++ ) {
	for ( int sb=0; sb<(co->nstypes); sb++ ) { isab = sa*co->nstypes+sb;
	for ( int ia=0; ia<(int)co->stypeindex[sa].size(); ia++ ) {
		int a=co->stypeindex[sa][ia];a_start=co->mystart(a);
		for ( int i=a_start; i<(a_start+co->item[a]); i++ ) {
			if ( !cv->nnidsab[isab][i].empty() ) {
				for ( int ij=0; ij<(int)cv->nnidsab[isab][i].size(); ij++ ) {
					// calc r(i)'s pair position of r(nnid[i][ij])
					// from normal vector of voronoi surface
					// which direction is r(i) to r(nnid[i][ij])
					for ( int j=0; j<3; j++ ) {
						nwe(0,j) = cv->egsab[isab][i][ij][0](j);
						nwe(1,j) = cv->egsab[isab][i][ij][1](j);
						nwe(2,j) = cv->nvsab[isab][i][ij](j);
					}
					pair(0) = cv->config[i]*cv->egsab[isab][i][ij][0];
					pair(1) = cv->config[i]*cv->egsab[isab][i][ij][1];
					pair(2) = (2.0*cv->egsab[isab][i][ij][2]-cv->config[i])*cv->nvsab[isab][i][ij];
					pair = nwe.inv3d()*pair;
	
					// save middle point of the bond not duplicately
					// and its direction; normal vector
					if ( sa == sb ) {
						bmp = 0.5*(cv->config[i]+pair);
						// erase uncertainty
						for ( int ijk=0; ijk<3; ijk++ ) {
							if ( pow(bmp(ijk),2.0) < under_capacity )
								bmp(ijk) = 0.0;
						}
						ss << bmp(0) << "," << bmp(1) << "," << bmp(2);
						bmpstr = ss.str();
						ss.str(""); ss.clear();
						if ( !duple[bmpstr] ) {
							duple[bmpstr] = true;
							isabs[isab].push_back(isab);
							ssab[isab].push_back(cv->surfsab[isab][i][ij]);
							bondsab[isab].push_back(bmp);
							dircsab[isab].push_back(cv->nvsab[isab][i][ij]);
						}
					} else {
						int b=co->mytype(cv->nnidsab[isab][i][ij]);
						vector<int>::iterator ib=find(co->stypelist[b].begin(),co->stypelist[b].end(),sb);
						if ( ib != co->stypelist[b].end() ) {
							bmp = 0.5*(cv->config[i]+pair);
							// erase uncertainty
							for ( int ijk=0; ijk<3; ijk++ ) {
								if ( pow(bmp(ijk),2.0) < under_capacity )
									bmp(ijk) = 0.0;
							}
							ss << bmp(0) << "," << bmp(1) << "," << bmp(2);
							bmpstr = ss.str();
							ss.str(""); ss.clear();
							if ( !duple[bmpstr] ) {
								duple[bmpstr] = true;
								isabs[isab].push_back(isab);
								ssab[isab].push_back(cv->surfsab[isab][i][ij]);
								bondsab[isab].push_back(bmp);
								dircsab[isab].push_back(cv->nvsab[isab][i][ij]);
							}
						}
					}
	
				}
			}
		}
		// free boolean duple array
		duple.clear();
	}
	}
	}

	// mpi parameter
	mpi_local_part(&local_start,&local_size,(int)bond.size(),rank,nmpi);
	// calculate G0(r)*Gl(r)
	for ( int i=local_start; i<local_start+local_size; i++ ) {
	for ( int j=0; j<(int)bond.size(); j++ ) {
	if ( i != j ) {
		for ( int t1=-(co->trans); t1<=(co->trans); t1++ ) {
		for ( int t2=-(co->trans); t2<=(co->trans); t2++ ) {
		for ( int t3=-(co->trans); t3<=(co->trans); t3++ ) {
			T = (double)t1*cv->cell[0]+(double)t2*cv->cell[1]+(double)t3*cv->cell[2];
			r = (bond[i]-bond[j]+T).norm();
			ir = (int)((r-co->r0)/co->dr);
			if ( ir < nr ) {
				sr[ir] += s[i]*s[j];
				cart2shpr(&ra,&thetaa,&phia,dirc[i](0),dirc[i](1),dirc[i](2));
				cart2shpr(&rb,&thetab,&phib,dirc[j](0),dirc[j](1),dirc[j](2));
				cosw = cos(thetaa)*cos(thetab)+sin(thetaa)*sin(thetab)*cos(phia-phib);
				for ( int il=0; il<nl; il++ ) {
					l = co->l[il];
					ssplr[il][ir] += s[i]*s[j]*Legendre_Pl(l,cosw);
				}
			}
		}
		}
		}
	}
	}
	}
	MPI_Barrier(MPI_COMM_WORLD);

	// calculate G0abab(r)*Glabab(r)
	for ( iab=0; iab<(co->ntypes*co->ntypes); iab++ ) {
		// mpi parameter
		mpi_local_part(&local_start,&local_size,(int)bondab[iab].size(),rank,nmpi);
		for ( int i=local_start; i<local_start+local_size; i++ ) {
		for ( int icd=0; icd<(co->ntypes*co->ntypes); icd++ ) {
		for ( int j=0; j<(int)bondab[icd].size(); j++ ) {
			for ( int t1=-(co->trans); t1<=(co->trans); t1++ ) {
			for ( int t2=-(co->trans); t2<=(co->trans); t2++ ) {
			for ( int t3=-(co->trans); t3<=(co->trans); t3++ ) {
				T = (double)t1*cv->cell[0]+(double)t2*cv->cell[1]+(double)t3*cv->cell[2];
				r = (bondab[icd][j]-bondab[iab][i]+T).norm();
				ir = (int)((r-co->r0)/co->dr);
				if ( under_capacity < ir && ir < nr ) {
					sababr[iabs[iab][i]*co->ntypes*co->ntypes+iabs[icd][j]][ir] += sab[iab][i]*sab[icd][j];
					cart2shpr(&ra,&thetaa,&phia,dircab[iab][i](0),dircab[iab][i](1),dircab[iab][i](2));
					cart2shpr(&rb,&thetab,&phib,dircab[icd][j](0),dircab[icd][j](1),dircab[icd][j](2));
					cosw = cos(thetaa)*cos(thetab)+sin(thetaa)*sin(thetab)*cos(phia-phib);
					for ( int il=0; il<nl; il++ ) {
						l = co->l[il];
						ssplababr[il][iabs[iab][i]*co->ntypes*co->ntypes+iabs[icd][j]][ir] += sab[iab][i]*sab[icd][j]*Legendre_Pl(l,cosw);
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

	// calculate G0sabab(r)*Glsabab(r)
	for ( isab=0; isab<(co->nstypes*co->nstypes); isab++ ) {
		// mpi parameter
		mpi_local_part(&local_start,&local_size,(int)bondsab[isab].size(),rank,nmpi);
		for ( int i=local_start; i<local_start+local_size; i++ ) {
		for ( int iscd=0; iscd<(co->nstypes*co->nstypes); iscd++ ) {
		for ( int j=0; j<(int)bondsab[iscd].size(); j++ ) {
			for ( int t1=-(co->trans); t1<=(co->trans); t1++ ) {
			for ( int t2=-(co->trans); t2<=(co->trans); t2++ ) {
			for ( int t3=-(co->trans); t3<=(co->trans); t3++ ) {
				T = (double)t1*cv->cell[0]+(double)t2*cv->cell[1]+(double)t3*cv->cell[2];
				r = (bondsab[iscd][j]-bondsab[isab][i]+T).norm();
				ir = (int)((r-co->r0)/co->dr);
				if ( under_capacity < ir && ir < nr ) {
					ssababr[isabs[isab][i]*co->nstypes*co->nstypes+isabs[iscd][j]][ir] += ssab[isab][i]*ssab[iscd][j];
					cart2shpr(&ra,&thetaa,&phia,dircsab[isab][i](0),dircsab[isab][i](1),dircsab[isab][i](2));
					cart2shpr(&rb,&thetab,&phib,dircsab[iscd][j](0),dircsab[iscd][j](1),dircsab[iscd][j](2));
					cosw = cos(thetaa)*cos(thetab)+sin(thetaa)*sin(thetab)*cos(phia-phib);
					for ( int il=0; il<nl; il++ ) {
						l = co->l[il];
						ssplsababr[il][isabs[isab][i]*co->nstypes*co->nstypes+isabs[iscd][j]][ir] += ssab[isab][i]*ssab[iscd][j]*Legendre_Pl(l,cosw);
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

	// free voro++ output
	vector < double > ().swap(s);
	vector < vec > ().swap(bond);
	vector < vec > ().swap(dirc);
	vector < vector < int > > ().swap(iabs);
	vector < vector < double > > ().swap(sab);
	vector < vector < vec > > ().swap(bondab);
	vector < vector < vec > > ().swap(dircab);
	vector < vector < int > > ().swap(isabs);
	vector < vector < double > > ().swap(ssab);
	vector < vector < vec > > ().swap(bondsab);
	vector < vector < vec > > ().swap(dircsab);

	// memory allocation for data gathering
	// normalization surface area s12(r)
	// radial correlation function s1(r)*s2(r)*Pl(cos(w12))
	sr0 = alloc1d < double > (nr);
	sababr0 = alloc2d < double > (co->ntypes*co->ntypes*co->ntypes*co->ntypes,nr);
	ssababr0 = alloc2d < double > (co->nstypes*co->nstypes*co->nstypes*co->nstypes,nr);
	ssplr0 = alloc2d < double > (nl,nr);
	ssplababr0 = alloc3d < double > (nl,co->ntypes*co->ntypes*co->ntypes*co->ntypes,nr);
	ssplsababr0 = alloc3d < double > (nl,co->nstypes*co->nstypes*co->nstypes*co->nstypes,nr);
	MPI_Barrier(MPI_COMM_WORLD);

	// data gathering
	MPI_Allreduce(&(sr[0]),&(sr0[0]),nr,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(&(sababr[0][0]),&(sababr0[0][0]),co->ntypes*co->ntypes*co->ntypes*co->ntypes*nr,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(&(ssababr[0][0]),&(ssababr0[0][0]),co->nstypes*co->nstypes*co->nstypes*co->nstypes*nr,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(&(ssplr[0][0]),&(ssplr0[0][0]),nl*nr,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(&(ssplababr[0][0][0]),&(ssplababr0[0][0][0]),nl*co->ntypes*co->ntypes*co->ntypes*co->ntypes*nr,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(&(ssplsababr[0][0][0]),&(ssplsababr0[0][0][0]),nl*co->nstypes*co->nstypes*co->nstypes*co->nstypes*nr,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	// free local data
	free1d < double > (sr);
	free2d < double > (sababr);
	free2d < double > (ssababr);
	free2d < double > (ssplr);
	free3d < double > (ssplababr);
	free3d < double > (ssplsababr);

	// update s1(r)*s2(r)*Pl(cos(w12))/s12(r)
	for ( int il=0; il<nl; il++ ) {
	for ( ir=0; ir<nr; ir++ ) {

		// total bacf
		if ( sr0[ir] ) gglr[il][ir] = co->update(gglr[il][ir],ssplr0[il][ir]/sr0[ir]);
		else gglr[il][ir] = co->update(gglr[il][ir],0);

		// partial bacf
		for ( int iabab=0; iabab<(co->ntypes*co->ntypes*co->ntypes*co->ntypes); iabab++ ) {
			if ( sababr0[iabab][ir] ) {
				gglababr[il][iabab][ir] = co->update(gglababr[il][iabab][ir],ssplababr0[il][iabab][ir]/sababr0[iabab][ir]);
			} else {
				gglababr[il][iabab][ir] = co->update(gglababr[il][iabab][ir],0);
			}
		}

		// specified partial bacf
		for ( int isabab=0; isabab<(co->nstypes*co->nstypes*co->nstypes*co->nstypes); isabab++ ) {
			if ( ssababr0[isabab][ir] ) {
				gglsababr[il][isabab][ir] = co->update(gglsababr[il][isabab][ir],ssplsababr0[il][isabab][ir]/ssababr0[isabab][ir]);
			} else {
				gglsababr[il][isabab][ir] = co->update(gglsababr[il][isabab][ir],0);
			}
		}

	}
	}
	MPI_Barrier(MPI_COMM_WORLD);

	// free integrated data
	free1d < double > (sr0);
	free2d < double > (sababr0);
	free2d < double > (ssababr0);
	free2d < double > (ssplr0);
	free3d < double > (ssplababr0);
	free3d < double > (ssplsababr0);

}

void bacf::write ( string filename ) {
	double r;
	ofstream ofs,ofsabab,ofssabab;

	// write out
	ofs.open((filename+"-tot.dat").c_str(),ios::out);
	ofsabab.open((filename+".dat").c_str(),ios::out);
	ofssabab.open((filename+"-spec.dat").c_str(),ios::out);
	ofs << scientific << setprecision(5);
	ofsabab << scientific << setprecision(5);
	ofssabab << scientific << setprecision(5);
	// total G0(r)*Gl(r)
	for ( int ir=0; ir<nr; ir++ ) {
		r = (double)ir*co->dr + co->r0;
		ofs << r;
		for ( int il=0; il<nl; il++ ) ofs << "\t" << gglr[il][ir];
		ofs << endl;
	}
	// partial G0(r)*Gl(r0
	for ( int iabab=0; iabab<(co->ntypes*co->ntypes*co->ntypes*co->ntypes); iabab++ ) {
		for ( int ir=0; ir<nr; ir++ ) {
			r = (double)ir*co->dr + co->r0;
			ofsabab << r;
			for ( int il=0; il<nl; il++ ) ofsabab << "\t" << gglababr[il][iabab][ir];
			ofsabab << endl;
		}
		ofsabab << endl;
	}
	// specified partial G0(r)*Gl(r0
	for ( int isabab=0; isabab<(co->nstypes*co->nstypes*co->nstypes*co->nstypes); isabab++ ) {
		for ( int ir=0; ir<nr; ir++ ) {
			r = (double)ir*co->dr + co->r0;
			ofssabab << r;
			for ( int il=0; il<nl; il++ ) ofssabab << "\t" << gglsababr[il][isabab][ir];
			ofssabab << endl;
		}
		ofssabab << endl;
	}
	ofs.close();
	ofsabab.close();
	ofssabab.close();

}
