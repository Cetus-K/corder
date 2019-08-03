// lceps.h
#ifndef LCEPS_H
#define LCEPS_H

#include "header.h"

// Levi-Civita symbol ( double )
// sgn(sgm) = -1,+1,0
double LC_eps ( vector < int > ids ) {
	int idx,len;
	map < int,int > sgm;
	vector < int > pmt;
	vector < vector < int > > pmt2;
	for ( int i=0; i<ids.size(); i++ ) sgm[i] = ids[i];
	while ( !sgm.empty() ) {
		idx = sgm.begin()->first;
		pmt.push_back(idx);
		while ( sgm.begin()->first != sgm[idx] ) {
			idx = sgm[idx];
			pmt.push_back(idx);
		}
		for ( int i=0; i<pmt.size(); i++ ) sgm.erase(pmt[i]); 
		pmt2.push_back(pmt);
		pmt.clear();
	}
	len = 0;
	for ( int i=0; i<pmt2.size(); i++ ) {
		len += pmt2[i].size()-1;
	}
	vector < int > ().swap(pmt);
	vector < vector < int > > ().swap(pmt2);
	return pow(-1.0,(double)len);
}

#endif
