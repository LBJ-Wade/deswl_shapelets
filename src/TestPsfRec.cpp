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

std::ostream* dbgout = 0;
bool XDEBUG = true;


int main(int argc, char **argv) 
{

    if (argc < 4) {
        std::cout<<"usage: test_psfrec config psf_file fitpsf_file\n";
        exit(45);
    }

    std::string configFile=argv[1];
    std::string psfFile=argv[2];
    std::string fitPsfFile=argv[3];

    std::cerr<<"\n";
    std::cerr<<"config_file: "<<configFile<<"\n";
    std::cerr<<"psf_file: "<<psfFile<<"\n";
    std::cerr<<"fitpsf_file: "<<fitPsfFile<<"\n";

    std::cerr<<"Loading config...\n";
    ConfigFile params(configFile);
    std::string fp((const char*)fitsparams_config,fitsparams_config_len);
    std::istringstream is(fp);
    params.read(is);


    std::cerr<<"Loading a PSFCatalog\n";
    PsfCatalog psfCat(params, psfFile);

    std::cerr<<"Loading a FittedPSF\n";
    FittedPsf fitPsf(params);
    fitPsf.read(fitPsfFile);

    BVec iPsf(fitPsf.getPsfOrder(), fitPsf.getSigma());

    // Test the reconstructions
    std::cout<<"NROWS = "<<psfCat.size()<<"\n";
    std::cout<<"{'_DELIM': ' ',\n";
    std::cout<<" '_DTYPE': [('id','i4'),\n";
    std::cout<<"            ('psf_flags','i4'),\n";
    std::cout<<"            ('x','f4'),('y','f4'),\n";
    std::cout<<"            ('e1','f4'),('e2','f4'),\n";
    std::cout<<"            ('e1interp','f4'),('e2interp','f4')],\n";
    std::cout<<" 'config_file': '"<<configFile<<"',\n";
    std::cout<<" 'psf_file': '"<<psfFile<<"',\n";
    std::cout<<" 'fitpsf_file': '"<<fitPsfFile<<"'}\n";
    std::cout<<"END\n";
    std::cout<<"\n";
    const int nStars = psfCat.size();
    for (int i=0; i<nStars; ++i) {

        double e1 = sqrt(2)*psfCat.getPsf(i)(3);
        double e2 = sqrt(2)*psfCat.getPsf(i)(4);

        fitPsf.interpolate(psfCat.getPos(i), iPsf);

        double ie1 = sqrt(2)*iPsf(3);
        double ie2 = sqrt(2)*iPsf(4);
        std::cout<<psfCat.getId(i)
            <<" "<<psfCat.getFlags(i)
            <<" "<<psfCat.getPos(i).getX()
            <<" "<<psfCat.getPos(i).getY()
            <<" "<<e1
            <<" "<<e2
            <<" "<<ie1
            <<" "<<ie2
            <<"\n";
    }

    exit(0);

}
