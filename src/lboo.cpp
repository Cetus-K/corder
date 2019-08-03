#include "corder.h"

void lboo::init ( void ) {

	// temporary variable
	nl = (int)co->l.size();

}

void lboo::write ( string filename ) {

	// temporary variable
	int a_start=0,iab,isab;	// voro++
	double *s,**sb,**ssb;
	int l=0,m,m1,m2,m3;		// lboo: temporary
	double r,theta,phi;
	complex < double > Ylm(0,0),wl_buf(0,0);
	double **ql,**wl,***qlb,***wlb,***qlsb,***wlsb;	// lboo: saving
	complex < double > ***qlm,****qlmb,****qlmsb;
	ofstream ofsq,ofsw,ofsqab,ofswab,ofsqsab,ofswsab;

	// total surface area
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
			if ( b == co->mytype(cv->nnidab[iab][i][ij]) )
				sb[i][b] += cv->surfab[iab][i][ij];
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

	// check
//	cout << "partial summation surfsce area:" << endl;
//	for ( int i=0; i<(co->nions); i++ ) {
//		cout << "ion-" << i+1 << ":" << endl;
//		for ( int b=0; b<(co->ntypes); b++ ) {
//			cout << "\tS(" << i+1 << "," << b << ") = " << sb[i][b] << endl;
//		}
//	}
//	cout << endl;
//	cout << "partial surface area:" << endl;
//	for ( int a=0; a<(co->ntypes); a++ ) {
//	for ( int b=0; b<(co->ntypes); b++ ) {
//		a_start = co->mystart(a);
//		iab = a*co->ntypes+b;
//		// for na
//		if ( a == b ) {
//			for ( int i=a_start; i<(a_start+co->item[a]); i++ ) {
//				cout << "ion-" << i+1 << ":" << endl;
//				for ( int ij=0; ij<(int)cv->surfab[iab][i].size(); ij++ ) {
//					cout << "\tbond-" << cv->nnidab[iab][i][ij]+1;
//					cout << "; S(" << ij+1 << ") = " <<  cv->surfab[iab][i][ij] << endl;
//				}
//			}
//		} else {
//			for ( int i=a_start; i<(a_start+co->item[a]); i++ ) {
//				cout << "ion-" << i+1 << ":" << endl;
//				for ( int ij=0; ij<(int)cv->surfab[iab][i].size(); ij++ ) {
//					if ( b == co->mytype(cv->nnidab[iab][i][ij]) ) {
//						cout << "\tbond-" << cv->nnidab[iab][i][ij]+1;
//						cout << "; S(" << ij+1 << ") = " <<  cv->surfab[iab][i][ij] << endl;
//					}
//				}
//			}
//		}
//	}
//	}
//	exit(0);

	// total/partial local-boo parameter
	ql = alloc2d < double > (nl,co->nions);
	wl = alloc2d < double > (nl,co->nions);
	qlb = alloc3d < double > (nl,co->nions,co->ntypes);
	wlb = alloc3d < double > (nl,co->nions,co->ntypes);
	qlsb = alloc3d < double > (nl,co->nions,co->nstypes);
	wlsb = alloc3d < double > (nl,co->nions,co->nstypes);

	// jag array
	qlm = new complex < double >** [nl];
	qlmb = new complex < double >*** [nl];
	qlmsb = new complex < double >*** [nl];
	for ( int il=0; il<nl; il++ ) {
		l = co->l[il];
		qlm[il] = new complex < double >* [2*l+1];
		qlmb[il] = new complex < double >** [2*l+1];
		qlmsb[il] = new complex < double >** [2*l+1];
		for ( m=-l; m<=l; m++ ) {
			qlm[il][m+l] = alloc1d < complex < double > > (co->nions);
			qlmb[il][m+l] = alloc2d < complex < double > > (co->nions,co->ntypes);
			qlmsb[il][m+l] = alloc2d < complex < double > > (co->nions,co->nstypes);
		}
	}

	// calc ql,qlb,wl,wlb,qlsb,wlsb
	for ( int il=0; il<nl; il++ ) {
		l = co->l[il];
		// total information
		for ( int i=0; i<(co->nions); i++ ) {

			// calc qlm
			for ( m=-l; m<=l; m++ ) {
				if ( !cv->nnid[i].empty() ) {
					for ( int ij=0; ij<(int)cv->nnid[i].size(); ij++ ) {
						cart2shpr(&r,&theta,&phi,cv->nv[i][ij](0),cv->nv[i][ij](1),cv->nv[i][ij](2));
						Ylm = spherical_Ylm(l,m,theta,phi);
						// qlm
						qlm[il][m+l][i] += cv->surf[i][ij]*conj(Ylm);
					}
					qlm[il][m+l][i] /= s[i];
				}
				// sum_{ |m|<=l }[ | qlm(i) |^2 ] = (2l+1)/4pi*ql^2
				// needs to multiply 4pi/(2l+1) and square-root it
				ql[il][i] += norm(qlm[il][m+l][i]);
			}
	
			// modulation to ql
			if ( ql[il][i] < under_capacity ) ql[il][i] = 0;
	
			// calc wl
			if ( ql[il][i] ) {
				wl_buf = complex < double > (0,0);
				for ( m1=-l; m1<=l; m1++ ) {
				for ( m2=max(-l,-(l+m1)); m2<=min(l,l-m1); m2++ ) {
					m3 = -m1-m2;
					wl_buf += Wigner_3j(l,l,l,m1,m2,m3)*qlm[il][m1+l][i]*qlm[il][m2+l][i]*qlm[il][m3+l][i];
				}
				}
				wl[il][i] = wl_buf.real()/pow(ql[il][i],1.5);
			} else { wl[il][i] = nan; }
			// calc ql
			ql[il][i] = sqrt(4.0*pi*ql[il][i]/(double)(2*l+1));
	
		}
		// partial information
		for ( int a=0; a<(co->ntypes); a++ ) { a_start = co->mystart(a);
		for ( int b=0; b<(co->ntypes); b++ ) { iab = a*co->ntypes+b;
		for ( int i=a_start; i<(a_start+co->item[a]); i++ ) {

			// calc qlmb
			for ( m=-l; m<=l; m++ ) {
				if ( !cv->nnidab[iab][i].empty() ) {
					for ( int ij=0; ij<(int)cv->nnidab[iab][i].size(); ij++ ) {
						cart2shpr(&r,&theta,&phi,cv->nvab[iab][i][ij](0),cv->nvab[iab][i][ij](1),cv->nvab[iab][i][ij](2));
						Ylm = spherical_Ylm(l,m,theta,phi);
						// qlmb
						// for na
						if ( a == b ) {
							qlmb[il][m+l][i][b] += cv->surfab[iab][i][ij]*conj(Ylm);
						// for nb of na+nb
						} else {
							if ( b == co->mytype(cv->nnidab[iab][i][ij]) ) {
								qlmb[il][m+l][i][b] += cv->surfab[iab][i][ij]*conj(Ylm);
							}
						}
					}
					qlmb[il][m+l][i][b] /= sb[i][b];
				}
				// sum_{ |m|<=l }[ | qlm(i) |^2 ] = (2l+1)/4pi*ql^2
				// needs to multiply 4pi/(2l+1) and square-root it
				qlb[il][i][b] += norm(qlmb[il][m+l][i][b]);
			}
	
			// modulation to qlb
			if ( qlb[il][i][b] < under_capacity ) qlb[il][i][b] = 0;
	
			// calc qlb,wlb
			if ( qlb[il][i][b] ) {
				wl_buf = complex < double > (0,0);
				for ( m1=-l; m1<=l; m1++ ) {
				for ( m2=max(-l,-(l+m1)); m2<=min(l,l-m1); m2++ ) {
					m3 = -m1-m2;
					wl_buf += Wigner_3j(l,l,l,m1,m2,m3)*qlmb[il][m1+l][i][b]*qlmb[il][m2+l][i][b]*qlmb[il][m3+l][i][b];
				}
				}
				wlb[il][i][b] = wl_buf.real()/pow(qlb[il][i][b],1.5);
			} else { wlb[il][i][b] = nan; }
			// calc qlb
			qlb[il][i][b] = sqrt(4.0*pi*qlb[il][i][b]/(double)(2*l+1));
	
		}
		}
		}
		// specified partial surface area
		for ( int spa=0; spa<(co->nstypes); spa++ ) {
		for ( int spb=0; spb<(co->nstypes); spb++ ) { isab = spa*co->nstypes+spb;
		for ( int ia=0; ia<(int)co->stypeindex[spa].size(); ia++ ) {
			int a=co->stypeindex[spa][ia];a_start=co->mystart(a);
			for ( int i=a_start; i<(a_start+co->item[a]); i++ ) {
				// calc qlmsb
				for ( m=-l; m<=l; m++ ) {
					if ( !cv->nnidsab[isab][i].empty() ) {
						for ( int ij=0; ij<(int)cv->nnidsab[isab][i].size(); ij++ ) {
							cart2shpr(&r,&theta,&phi,cv->nvsab[isab][i][ij](0),cv->nvsab[isab][i][ij](1),cv->nvsab[isab][i][ij](2));
							Ylm = spherical_Ylm(l,m,theta,phi);
							// qlmb
							// for na
							if ( spa == spb ) {
								qlmsb[il][m+l][i][spb] += cv->surfsab[isab][i][ij]*conj(Ylm);
							// for nb of na+nb
							} else {
								int b=co->mytype(cv->nnidsab[isab][i][ij]);
								vector<int>::iterator ib=find(co->stypelist[b].begin(),co->stypelist[b].end(),spb);
								if ( ib != co->stypelist[b].end() )
									qlmsb[il][m+l][i][spb] += cv->surfsab[isab][i][ij]*conj(Ylm);
							}
						}
						qlmsb[il][m+l][i][spb] /= ssb[i][spb];
					}
					// sum_{ |m|<=l }[ | qlm(i) |^2 ] = (2l+1)/4pi*ql^2
					// needs to multiply 4pi/(2l+1) and square-root it
					qlsb[il][i][spb] += norm(qlmsb[il][m+l][i][spb]);
				}
		
				// modulation to qlsb
				if ( qlsb[il][i][spb] < under_capacity ) qlsb[il][i][spb] = 0;
		
				// calc qlsb,wlsb
				if ( qlsb[il][i][spb] ) {
					wl_buf = complex < double > (0,0);
					for ( m1=-l; m1<=l; m1++ ) {
					for ( m2=max(-l,-(l+m1)); m2<=min(l,l-m1); m2++ ) {
						m3 = -m1-m2;
						wl_buf += Wigner_3j(l,l,l,m1,m2,m3)*qlmsb[il][m1+l][i][spb]*qlmsb[il][m2+l][i][spb]*qlmsb[il][m3+l][i][spb];
					}
					}
					wlsb[il][i][spb] = wl_buf.real()/pow(qlsb[il][i][spb],1.5);
				} else { wlsb[il][i][spb] = nan; }
				// calc qlsb
				qlsb[il][i][spb] = sqrt(4.0*pi*qlsb[il][i][spb]/(double)(2*l+1));
			}
		}
		}
		}
	}

	// write out
	if ( co->itr == co->f0 ) {
		ofsq.open((filename+"-q-tot.dat").c_str(),ios::out);
		ofsw.open((filename+"-w-tot.dat").c_str(),ios::out);
		ofsqab.open((filename+"-q.dat").c_str(),ios::out);
		ofswab.open((filename+"-w.dat").c_str(),ios::out);
		ofsqsab.open((filename+"-q-spec.dat").c_str(),ios::out);
		ofswsab.open((filename+"-w-spec.dat").c_str(),ios::out);
	} else {
		ofsq.open((filename+"-q-tot.dat").c_str(),ios::app);
		ofsw.open((filename+"-w-tot.dat").c_str(),ios::app);
		ofsqab.open((filename+"-q.dat").c_str(),ios::app);
		ofswab.open((filename+"-w.dat").c_str(),ios::app);
		ofsqsab.open((filename+"-q-spec.dat").c_str(),ios::app);
		ofswsab.open((filename+"-w-spec.dat").c_str(),ios::app);
	}
	ofsq << scientific << setprecision(5);
	ofsw << scientific << setprecision(5);
	ofsqab << scientific << setprecision(5);
	ofswab << scientific << setprecision(5);
	ofsqsab << scientific << setprecision(5);
	ofswsab << scientific << setprecision(5);

	if ( nl ) {

		// ql,wl
		for ( int i=0; i<(co->nions); i++ ) {
			ofsq << ql[0][i]; ofsw << wl[0][i];
			for ( int il=1; il<nl; il++ ) {
				ofsq << "\t" << ql[il][i];
				ofsw << "\t" << wl[il][i];
			}
			ofsq << endl; ofsw << endl;
		}
		ofsq << endl; ofsw << endl;
		ofsq << endl; ofsw << endl;

		// qlab,wlab
		for ( int a=0; a<(co->ntypes); a++ ) {
			a_start = co->mystart(a);
			for ( int b=0; b<(co->ntypes); b++ ) {
				for ( int i=a_start; i<(a_start+co->item[a]); i++ ) {
					ofsqab << qlb[0][i][b]; ofswab << wlb[0][i][b];
					for ( int il=1; il<nl; il++ ) {
						ofsqab << "\t" << qlb[il][i][b];
						ofswab << "\t" << wlb[il][i][b];
					}
					ofsqab << endl; ofswab << endl;
				}
				ofsqab << endl; ofswab << endl;
			}
		}
		ofsqab << endl; ofswab << endl;

		// qlsab,wlsab
		for ( int spa=0; spa<(co->nstypes); spa++ ) {
		for ( int spb=0; spb<(co->nstypes); spb++ ) {
			for ( int ia=0; ia<(int)co->stypeindex[spa].size(); ia++ ) {
				int a=co->stypeindex[spa][ia];a_start=co->mystart(a);
				for ( int i=a_start; i<(a_start+co->item[a]); i++ ) {
					ofsqsab << qlsb[0][i][spb]; ofswsab << wlsb[0][i][spb];
					for ( int il=1; il<nl; il++ ) {
						ofsqsab << "\t" << qlsb[il][i][spb];
						ofswsab << "\t" << wlsb[il][i][spb];
					}
					ofsqsab << endl; ofswsab << endl;
				}
			}
			ofsqsab << endl; ofswsab << endl;
		}
		}
		ofsqsab << endl; ofswsab << endl;

	}
	ofsq.close();
	ofsw.close();
	ofsqab.close();
	ofswab.close();
	ofsqsab.close();
	ofswsab.close();

	// free
	free1d < double > (s);
	free2d < double > (sb);
	free2d < double > (ssb);
	free2d < double > (ql);
	free2d < double > (wl);
	free3d < double > (qlb);
	free3d < double > (wlb);
	free3d < double > (qlsb);
	free3d < double > (wlsb);
	for ( int il=0; il<nl; il++ ) {
		l = co->l[il];
		for ( m=-l; m<=l; m++ ) {
			free1d < complex < double > > (qlm[il][m+l]);
			free2d < complex < double > > (qlmb[il][m+l]);
			free2d < complex < double > > (qlmsb[il][m+l]);
		}
		delete [] qlm[il]; qlm[il] = NULL;
		delete [] qlmb[il]; qlmb[il] = NULL;
		delete [] qlmsb[il]; qlmsb[il] = NULL;
	}
	delete [] qlm; qlm = NULL;
	delete [] qlmb; qlmb = NULL;
	delete [] qlmsb; qlmsb = NULL;

}
