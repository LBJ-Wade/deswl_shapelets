#include <iostream>
#include <sstream>
#include <vector>
#include <string>

#include "TMV.h"
#include "TMV_Small.h"

#include "dbg.h"

#include "PSFCatalog.h"
#include "FittedPSF.h"

#include "ConfigFile.h"
#include "fp.h" // Generated with xxd -i fitsparams.config fp.h

// for Position
#include "Bounds.h"

#include <CCfits/CCfits>

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
	
	dbgout = &std::cout;

	string config_file=argv[1];
	string psf_file=argv[2];
	string fitpsf_file=argv[3];


	cout<<"config_file: "<<config_file<<"\n";
	cout<<"psf_file: "<<psf_file<<"\n";
	cout<<"fitpsf_file: "<<fitpsf_file<<"\n";

	cout<<"Loading config...\n";
	ConfigFile params(config_file);
	std::string fp((const char*)fitsparams_config,fitsparams_config_len);
	std::istringstream is(fp);
	params.Read(is);


	cout<<"Loading a PSFCatalog\n";
	PSFCatalog psf(params, psf_file);

	cout<<"Loading a FittedPSF\n";
	FittedPSF fitpsf(params, fitpsf_file);

	// Test the reconstructions
	for (int i=0; i<psf.pos.size(); i++) {
		double b00 = psf.psf[i][0];
		cout<<"b00 = "<<b00<<"\n";
		cout<<"1 = "<<psf.psf[i][1]<<"\n";
		cout<<"2 = "<<psf.psf[i][2]<<"\n";
		double e1 = psf.psf[i][3]; // e1/sqrt(2)
		double e2 = psf.psf[i][4]; // e2/sqrt(2)
		cout<<"e1: "<<e1<<" e2: "<<e2<<"\n";
	}

	exit(0);

}
