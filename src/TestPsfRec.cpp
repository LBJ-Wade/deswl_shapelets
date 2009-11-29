#include <iostream>
#include <sstream>
#include <vector>
#include <string>

#include "dbg.h"

#include "PsfCatalog.h"
#include "FittedPsf.h"

#include "ConfigFile.h"
#include "fp.h" // Generated with xxd -i fitsparams.config fp.h

// for Position
#include "Bounds.h"

using namespace std;

std::ostream* dbgout = 0;
bool XDEBUG = true;

#define xdbgout (XDEBUG ? dbgout : 0)


int main(int argc, char **argv) 
{

	if (argc < 4) {
		std::cout<<"usage: test_psfrec config psf_file fitpsf_file\n";
		exit(45);
	}

	//dbgout = &std::cerr;

	string config_file=argv[1];
	string psf_file=argv[2];
	string fitpsf_file=argv[3];


	cerr<<"\n";
	cerr<<"config_file: "<<config_file<<"\n";
	cerr<<"psf_file: "<<psf_file<<"\n";
	cerr<<"fitpsf_file: "<<fitpsf_file<<"\n";

	cerr<<"Loading config...\n";
	ConfigFile params(config_file);
	string fp((const char*)fitsparams_config,fitsparams_config_len);
	istringstream is(fp);
	params.read(is);


	cerr<<"Loading a PSFCatalog\n";
	PsfCatalog psfcat(params, psf_file);

	cerr<<"Loading a FittedPSF\n";
	FittedPsf fitpsf(params, fitpsf_file);

	BVec ipsf(fitpsf.getPsfOrder(), fitpsf.getSigma());

	// Test the reconstructions
	cout<<"NROWS = "<<psfcat.size()<<"\n";
	cout<<"{'_DELIM': ' ',\n";
	cout<<" '_DTYPE': [('id','i4'),\n";
	cout<<"            ('psf_flags','i4'),\n";
	cout<<"            ('x','f4'),('y','f4'),\n";
	cout<<"            ('e1','f4'),('e2','f4'),\n";
	cout<<"            ('e1interp','f4'),('e2interp','f4')],\n";
	cout<<" 'config_file': '"<<config_file<<"',\n";
	cout<<" 'psf_file': '"<<psf_file<<"',\n";
	cout<<" 'fitpsf_file': '"<<fitpsf_file<<"'}\n";
	cout<<"END\n";
	cout<<"\n";
	const int nStars = psfcat.size();
	for (int i=0; i<nStars; i++) {

		double e1 = sqrt(2)*psfcat.getPsf(i)[3];
		double e2 = sqrt(2)*psfcat.getPsf(i)[4];

		fitpsf.interpolate(psfcat.getPos(i), ipsf);

		double ie1 = sqrt(2)*ipsf[3];
		double ie2 = sqrt(2)*ipsf[4];
		cout<<psfcat.getId(i)
			<<" "<<psfcat.getFlags(i)
			<<" "<<psfcat.getPos(i).getX()
			<<" "<<psfcat.getPos(i).getY()
			<<" "<<e1
			<<" "<<e2
			<<" "<<ie1
			<<" "<<ie2
			<<"\n";
	}

	exit(0);

}
