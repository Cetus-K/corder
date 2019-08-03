#include "corder.h"

void corder::read_input ( void ) {

	// fstream
	ifstream ifs;
	istringstream iss;
	string str;

	// temporary variable
	int num=0;
	vector < string > ist;

	ifs.open("param.in",ios::in);
	while ( !ifs.eof() ) {

		// get line in param.in
		getline(ifs,str); iss.str(str); iss >> str;

		// reference key instead of string
		if ( str == "mode:" ) {
			while ( iss.rdbuf() -> in_avail() > 1 ) {
				iss >> str;
				mode.push_back(str);
			}
		}
		if ( str == "lnum:" ) {
			while ( iss.rdbuf() -> in_avail() > 1 ) {
				iss >> num;
				l.push_back(num);
			}
			sort(l.begin(),l.end());
		}
		if ( str == "spec:" ) {
			// presave
			while ( iss.rdbuf() -> in_avail() > 1 ) {
				iss >> str;
				ist.push_back(str);
			}
			// save
			stype.push_back(ist);
		}
		if ( str == "flame:" ) iss >> f0 >> f1 >> df;
		if ( str == "rrange:" ) iss >> r0 >> r1 >> dr;
		if ( str == "qrange:" ) iss >> q0 >> q1 >> dq;
		if ( str == "qlrange:" ) iss >> ql0 >> ql1 >> dql;
		if ( str == "wlrange:" ) iss >> wl0 >> wl1 >> dwl;
		if ( str == "csrange:" ) iss >> cs0 >> cs1 >> dcs;
		if ( str == "tstep:" ) iss >> dt;
		if ( str == "trans:" ) iss >> trans;
		if ( str == "lcut:" ) iss >> lcut;
		if ( str == "dump:" ) iss >> dump;

		// clear stream
		iss.clear();
		if ( !ist.empty() ) vector < string > ().swap(ist);
	}
	ifs.close();

}

void corder::read_system ( void ) {

	// fstream
	ifstream ifs;
	istringstream iss;
	string str;

	// temporary variable
	int thistype;
	double scale;
	string label;
	vector < int > typecount;
	map < string,int > typeindex;
	map < string, map < int,bool > > typexist;

	// number of super cells
	ncells = (2*trans+1)*(2*trans+1)*(2*trans+1);

	// set vec and mat
	a1.resize(3);
	a2.resize(3);
	a3.resize(3);
	latmat.resize(3,3);

	// read XDATCAR
	ifs.open("XDATCAR",ios::in);
	// label
	getline(ifs,label);
	// scale
	getline(ifs,str); iss.str(str);
	iss >> scale; iss.clear();
	// lattice vectors
	getline(ifs,str); iss.str(str);
	iss >> a1(0) >> a1(1) >> a1(2);
	iss.clear(); a1 *= scale;
	getline(ifs,str); iss.str(str);
	iss >> a2(0) >> a2(1) >> a2(2);
	iss.clear(); a2 *= scale;
	getline(ifs,str); iss.str(str);
	iss >> a3(0) >> a3(1) >> a3(2);
	iss.clear(); a3 *= scale;
	for ( int i=0; i<3; i++ ) latmat(i,0) = a1(i);
	for ( int i=0; i<3; i++ ) latmat(i,1) = a2(i);
	for ( int i=0; i<3; i++ ) latmat(i,2) = a3(i);
	// species
	getline(ifs,str); iss.str(str);
	while ( iss.rdbuf() -> in_avail() > 1 ) {
		iss >> str;
		type.push_back(str);
	}
	iss.clear();
	ntypes = type.size();
	// items
	item.resize(ntypes);
	nions = 0;
	getline(ifs,str); iss.str(str);
	for ( int a=0; a<ntypes; a++ ) {
		iss >> item[a];
		nions += item[a];
		typeindex[type[a]] = a;
		typecount.push_back(0);
	}
	iss.clear();
	// axis
	getline(ifs,str);
	// configuration at first step
	for ( int i=0; i<nions; i++ ) getline(ifs,str);
	// second step
	getline(ifs,str);
	if ( str == label ) latdyn = true;
	else latdyn = false;
	// close XDATCAR
	ifs.close();

	// make spectypes array
	// stypeindex[a][j]= [j-th ai]
	// stypelist[ai][j-th ai] = a
	// ai: an atomic specie
	// j-th ai: j-th array of ai
	// a: partial total component including ai
	// example:
	// 	=============
	// 	  a | 0 1 2
	// 	-------------
	// 	    | 0 3 2
	// 	 ai | 1   3
	// 	    | 2
	// 	=============
	// 	stypeindex[0] = [0,1,2]
	// 	stypeindex[1] = [3]
	// 	stypeindex[2] = [2,3]
	// 	stypelist[0] = [0]
	// 	stypelist[1] = [0]
	// 	stypelist[2] = [0,2]
	// 	stypelist[3] = [1,2]
	// usage: when some indices v and w are given,
	// 	for iv=[0,stypeindex[v].size())
	// 	v_start = co->mystart(stypeindex[v][iv]);
	// usage: when some species v and w are given,
	// 	for iv=[0,stypelist[v].size()) and
	// 	    iw=[0,stypelist[w].size())
	// 	arr = stypelist[v][iv]*nstypes+stypelist[w][iw]
	if ( !stype.empty() ) {
		nstypes = (int)stype.size();
		stypeindex.resize(nstypes);
		stypelist.resize(ntypes);
		sitem.resize(nstypes);
		// set duple-cheking boolean true
		for ( int spa=0; spa<nstypes; spa++ ) {
		for ( int aj=0; aj<(int)stype[spa].size(); aj++ ) {
			thistype = typeindex[stype[spa][aj]];
			typexist[stype[spa][aj]][typecount[thistype]] = true;
			typecount[thistype] ++;
		}
		}
		// reset typecount
		for ( int a=0; a<ntypes; a++ ) typecount[a] = 0;
		for ( int spa=0; spa<nstypes; spa++ ) {
		for ( int aj=0; aj<(int)stype[spa].size(); aj++ ) {
			thistype = typeindex[stype[spa][aj]];
			if ( typexist[stype[spa][aj]][typecount[thistype]] ) {
				stypeindex[spa].push_back(thistype);
				stypelist[thistype].push_back(spa);
				sitem[spa] += item[thistype];
				typexist[stype[spa][aj]][typecount[thistype]] = false;
				typecount[thistype] ++;
			} else {
				cout << "ERROR: Specified specie " << stype[spa][aj];
				cout << " is not found in this system." << endl;
			}
		}
		}
	}
	vector < int > ().swap(typecount);
	typeindex.clear();
	typexist.clear();
}

void corder::set_param ( void ) {

	read_input();
	read_system();

}

void corder::info ( void ) {

	map < string,bool > flag;

	// set flag
	flag["rrange"] = true;
	flag["lnum"] = true;

	cout << endl;
	cout << "------------------------" << endl;
	cout << "     Input settings" << endl;
	cout << "------------------------" << endl;
	cout << endl;
	cout << "Calclation mode:" << endl;
	cout << "\t" << mode[0];
	for ( int i=1; i<(int)mode.size(); i++ ) cout << " " << mode[i];
	cout << endl << endl;
	cout << "Sampling flame:" << endl;
	cout << "\t" << f0 << " ~ " << f1 << " by " << df << endl;
	cout << endl;
	cout << "Super cell:" << endl;
	cout << "\t" << trans << " ( " << (int)(pow((double)(2*trans+1),3.0)) << " cells )" << endl;
	cout << endl;
	cout << "Write out inverval:" << endl;
	cout << "\t" << df*dump << " flames" << endl;
	cout << endl;
	for ( int i=0; i<(int)mode.size(); i++ ) {
		if ( mode[i] == "gofr" || mode[i] == "bacf" ) {
		if ( flag["rrange"] ) {
			cout << "r-range [ang]:" << endl;
			cout << "\t" << r0 << " ~ " << r1 << " by " << dr << endl;
			cout << endl;
			flag["rrange"] = false;
		}
		}
		if ( mode[i] == "sofq" ) {
			cout << "q-range:" << endl;
			cout << "\t" << q0 << " ~ " << q1 << " by " << dq << endl;
			cout << endl;
		}
		if ( mode[i] == "lboo" || mode[i] == "bacf" || mode[i] == "pofqw" ) {
		if ( flag["lnum"] ) {
			cout << "Angular l-number:" << endl;
			cout << "\t" << l[0];
			for ( int j=1; j<(int)l.size(); j++ ) cout << " " << l[j];
			cout << endl << endl;
			flag["lnum"] = false;
		}
		}
		if ( mode[i] == "pofqw" ) {
			cout << "ql-range:" << endl;
			cout << "\t" << ql0 << " ~ " << ql1 << " by " << dql << endl;
			cout << endl;
			cout << "wl-range:" << endl;
			cout << "\t" << wl0 << " ~ " << wl1 << " by " << dwl << endl;
			cout << endl;
		}
		if ( mode[i] == "dfc" ) {
			cout << "Time step [ps]:" << endl;
			cout << "\t" << dt << endl;
			cout << endl;
		}
		if ( mode[i] == "bofree" ) {
			cout << "Cut-off angular l-number:" << endl;
			cout << "\t" << lcut << endl;
			cout << endl;
		}
		if ( mode[i] == "csform" ) {
			cout << "cs-range:" << endl;
			cout << "\t" << cs0 << " ~ " << cs1 << " by " << dcs << endl;
			cout << endl;
		}
	}

	cout << endl;
	cout << "--------------------------" << endl;
	cout << "     system parameter" << endl;
	cout << "--------------------------" << endl;
	cout << endl;

	cout << "Number of ions / types:" << endl;
	cout << "\t" << nions << " / " << ntypes << endl;
	cout << endl;
	cout << "Species and its items:" << endl;
	cout << "\t" << item[0] << " " << type[0];
	for ( int i=1; i<ntypes; i++ ) {
		cout << " + " << item[i] << " " << type[i];
	}
	cout << endl;
	cout << endl;
	cout << "Lattice dynamics:" << endl;
	if ( latdyn ) cout << "\t" << "on" << endl;
	else cout << "\t" << "off" << endl;
	cout << endl;

}
