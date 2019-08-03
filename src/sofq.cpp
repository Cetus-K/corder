#include "corder.h"

void sofq::init ( void ) {

	// temporary variable
	nq = (int)((co->q1-co->q0)/co->dq);

}

void sofq::write ( string filename, gofr *gf ) {

	// temporary variable
	int iab,isab;
	double r,q,sincqr;
	double *sq,**sabq,**ssabq;
	ofstream ofs,ofsab,ofssab;

	// memory allocation
	sq = alloc1d < double > (nq);
	sabq = alloc2d < double > (co->ntypes*co->ntypes,nq);
	ssabq = alloc2d < double > (co->nstypes*co->nstypes,nq);

	// calc [s(q)-1]/(4*pi)/(dr/2)
	for ( int iq=0; iq<nq; iq++ ) {
		q = (double)iq*co->dq + co->q0;
		sq[iq] = 0;
		for ( iab=0; iab<(co->ntypes*co->ntypes); iab++ ) sabq[iab][iq] = 0;
		for ( isab=0; isab<(co->nstypes*co->nstypes); isab++ ) ssabq[isab][iq] = 0;
		for ( int ir=0; ir<(gf->nr); ir++ ) {
			r = (double)ir*co->dr + co->r0;

			// sinc function
			if ( q*r == 0 ) sincqr = 1.0;
			else sincqr = sin(q*r)/q*r;

			// trapezoidal integrate g(r)
			if ( ir == 0 || ir == gf->nr-1 ) {
				sq[iq] += (gf->rhor[ir]-gf->rho)*sincqr*r*r;
				for ( int a=0; a<(co->ntypes); a++ ) {
				for ( int b=0; b<(co->ntypes); b++ ) {
					iab = a*co->ntypes+b;
					sabq[iab][iq] += (gf->rhoabr[iab][ir]-gf->rhob[b])*sincqr*r*r;
				}
				}
				for ( int spa=0; spa<(co->nstypes); spa++ ) {
				for ( int spb=0; spb<(co->nstypes); spb++ ) {
					isab = spa*co->nstypes+spb;
					ssabq[isab][iq] += (gf->rhosabr[isab][ir]-gf->rhosb[spb])*sincqr*r*r;
				}
				}
			} else {
				sq[iq] += 2.0*(gf->rhor[ir]-gf->rho)*sincqr*r*r;
				for ( int a=0; a<co->ntypes; a++ ) {
				for ( int b=0; b<co->ntypes; b++ ) {
					iab = a*co->ntypes+b;
					sabq[iab][iq] += 2.0*(gf->rhoabr[iab][ir]-gf->rhob[b])*sincqr*r*r;
				}
				}
				for ( int spa=0; spa<(co->nstypes); spa++ ) {
				for ( int spb=0; spb<(co->nstypes); spb++ ) {
					isab = spa*co->nstypes+spb;
					ssabq[isab][iq] += 2.0*(gf->rhosabr[isab][ir]-gf->rhosb[spb])*sincqr*r*r;
				}
				}
			}

		}
	}

	// write out
	ofs.open((filename+"-tot.dat").c_str(),ios::out);
	ofsab.open((filename+".dat").c_str(),ios::out);
	ofssab.open((filename+"-spec.dat").c_str(),ios::out);
	ofs << scientific << setprecision(5);
	ofsab << scientific << setprecision(5);
	ofssab << scientific << setprecision(5);
	for ( int iq=0; iq<nq; iq++ ) {
		q = (double)iq*co->dq + co->q0;
		ofs << q << "\t" << 1.0+4.0*pi*sq[iq]*(0.5*co->dr) << endl;
		ofsab << q; ofssab << q;
		for ( iab=0; iab<(co->ntypes*co->ntypes); iab++ )
			ofsab << "\t" << 1.0+4.0*pi*sabq[iab][iq]*(0.5*co->dr);
		for ( isab=0; isab<(co->nstypes*co->nstypes); isab++ )
			ofssab << "\t" << 1.0+4.0*pi*ssabq[isab][iq]*(0.5*co->dr);
		ofsab << endl; ofssab << endl;
	}
	ofs.close();
	ofsab.close();
	ofssab.close();

	// free
	free1d < double > (sq);
	free2d < double > (sabq);
	free2d < double > (ssabq);

}
