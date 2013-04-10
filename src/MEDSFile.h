#ifndef MEDSFile_H
#define MEDSFile_H

#ifdef __INTEL_COMPILER
#pragma warning (disable : 1418)
#endif
#include "boost/shared_ptr.hpp"
#ifdef __INTEL_COMPILER
#pragma warning (default : 1418)
#endif

#include <vector>
#include <string>
extern "C" {
#include "meds.h"
}
#include "dbg.h"
#include "ConfigFile.h"
#include "FittedPsf.h"
#include "Pixel.h"

class MEDSFile 
{

public :

    // Just sets up the parameters.  Normally followed by read().
    MEDSFile(const ConfigFile& _params);
    ~MEDSFile();

    long size() const;

    void getPixels(std::vector<PixelList>& pix_list, long i, long& flag) const;
    void getPSFs(std::vector<BVec>& psf_list, long i, long& flag) const;

    std::string getCoaddCatFile() const;

private :

    struct meds* _medsPtr;

    // Data read from other files
    std::vector<boost::shared_ptr<FittedPsf> > _fitpsf;

    ConfigFile _params;

};

#endif
