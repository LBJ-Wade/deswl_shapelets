#include <iostream>
#include <sstream>
#include <vector>
#include <string>

#include "PsfCatalog.h"
#include "FittedPsf.h"
#include "ConfigFile.h"
#include "fp.h" // Generated with xxd -i fitsparams.config fp.h

// for Position
#include "Bounds.h"

#if defined (__INTEL_COMPILER) && defined(OPENMP_LINK)
__thread std::ostream* dbgout = 0;
__thread bool XDEBUG = false;
#else
std::ostream* dbgout = 0;
bool XDEBUG = false;
#endif

int main(int argc, char **argv) 
{

    if (argc < 4) {
        std::cout<<"usage: test_psfrec config psf_file fitpsf_file\n";
        exit(45);
    }

    std::string config_file=argv[1];
    std::string psf_file=argv[2];
    std::string fitpsf_file=argv[3];

    std::cerr<<"\n";
    std::cerr<<"config_file: "<<config_file<<"\n";
    std::cerr<<"psf_file: "<<psf_file<<"\n";
    std::cerr<<"fitpsf_file: "<<fitpsf_file<<"\n";

    std::cerr<<"Loading config...\n";
    ConfigFile params(config_file);
    std::string fp((const char*)fitsparams_config,fitsparams_config_len);
    std::istringstream is(fp);
    params.read(is);


    std::cerr<<"Loading a PSFCatalog\n";
    PsfCatalog psfcat(params);
    psfcat.read(psf_file);

    std::cerr<<"Loading a FittedPSF\n";
    FittedPsf fitpsf(params);
    fitpsf.read(fitpsf_file);

    BVec ipsf(fitpsf.getPsfOrder(), fitpsf.getSigma());

    // Test the reconstructions
    std::cout<<"NROWS = "<<psfcat.size()<<"\n";
    std::cout<<"{'_DELIM': ' ',\n";
    std::cout<<" '_DTYPE': [('id','i4'),\n";
    std::cout<<"            ('psf_flags','i4'),\n";
    std::cout<<"            ('x','f4'),('y','f4'),\n";
    std::cout<<"            ('e1','f4'),('e2','f4'),\n";
    std::cout<<"            ('e1interp','f4'),('e2interp','f4')],\n";
    std::cout<<" 'config_file': '"<<config_file<<"',\n";
    std::cout<<" 'psf_file': '"<<psf_file<<"',\n";
    std::cout<<" 'fitpsf_file': '"<<fitpsf_file<<"'}\n";
    std::cout<<"END\n";
    std::cout<<"\n";
    const int nstars = psfcat.size();
    for (int i=0; i<nstars; ++i) {

        double e1 = sqrt(2.)*psfcat.getPsf(i)(3);
        double e2 = sqrt(2.)*psfcat.getPsf(i)(4);

        fitpsf.interpolate(psfcat.getPos(i), ipsf);

        double ie1 = sqrt(2.)*ipsf(3);
        double ie2 = sqrt(2.)*ipsf(4);
        std::cout<<psfcat.getId(i)
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
