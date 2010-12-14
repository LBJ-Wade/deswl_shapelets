
#include <valarray>
#include "Image.h"
#include "FittedPsf.h"
#include "Log.h"
#include "CoaddCatalog.h"
#include "MultiShearCatalog.h"
#include "BasicSetup.h"
#include <sys/time.h>

#include <iostream>
#include <fstream>

#ifdef VG
#ifdef __INTEL_COMPILER
#pragma warning (disable : 593)
#pragma warning (disable : 1469)
#include <valgrind/memcheck.h>
#endif
#endif

static void doMeasureMultiShear(ConfigFile& params, ShearLog& log) 
{
    const double ARCSEC_PER_RAD = 206264.806247;

    bool isTiming = params.read("timing",false);
    timeval tp;
    double t1=0.,t2=0.;

    if (isTiming) {
        gettimeofday(&tp,0);
        t1 = tp.tv_sec + tp.tv_usec/1.e6;
    }

    bool shouldOutputDots = params.read("output_dots",false);

    // Read the coadd catalog.
    CoaddCatalog coaddCat(params);
    coaddCat.read();
    dbg<<"Made coaddcat\n";

    if (isTiming) {
        gettimeofday(&tp,0);
        t2 = tp.tv_sec + tp.tv_usec/1.e6;
        std::cout<<"Time: Read CoaddCatalog = "<<t2-t1<<std::endl;
        t1 = t2;
    }

    // Make the multishear catalog.
    MultiShearCatalog shearCat(coaddCat,params);
    dbg<<"Made multishearcat\n";
    dbg<<"Total bounds are "<<shearCat.getSkyBounds()<<std::endl;
    double area = shearCat.getSkyBounds().getArea();
    area /= 3600.;
    dbg<<"Raw area = "<<area<<" square arcmin\n";
    double dec = shearCat.getSkyBounds().getCenter().getY();
    double cosdec = cos(dec / ARCSEC_PER_RAD);
    area *= cosdec;
    dbg<<"Corrected area = "<<area<<" square arcmin";
    area /= 3600.;
    dbg<<" = "<<area<<" square degrees\n";

    log._nGals = shearCat.size();
    xdbg<<"nGals = "<<log._nGals<<std::endl;

    //On abe, this next line is giving a weird error and then aborting.
    //I have no idea what is going on here, so just avoid the stl count.
    //log._nGoodIn = std::count(
        //shearCat.getFlagsList().begin(),shearCat.getFlagsList().end(),0);
    int nGood = 0;
    for(size_t i=0;i<shearCat.getFlagsList().size();++i)
        if (shearCat.getFlagsList()[i] == 0) ++nGood;
    log._nGoodIn = nGood;
    xdbg<<"nGood = "<<log._nGoodIn<<std::endl;
    dbg<<log._nGoodIn<<"/"<<log._nGals<<" galaxies with no input flags\n";

    if (isTiming) {
        gettimeofday(&tp,0);
        t2 = tp.tv_sec + tp.tv_usec/1.e6;
        std::cout<<"Time: Make MultiShearCatalog = "<<t2-t1<<std::endl;
        t1 = t2;
    }

    double loadTime=0.;
    double calcTime=0.;

    // Break the area up into sections and load each section separately.
    // It is most efficient to load everything at once, but it can take 
    // a lot of memory.  So this bit keeps the memory requirement for each
    // section manageable.  So long as each section does a significant amount
    // of work, the extra I/O time won't be much of an issue. 
    // 20'-30' armin seems to be a good size to keep the memory under around
    // 4-5 GB, and only add about 5-10% to the running time.
    long nShear = 0;
    std::vector<Bounds> sectionBounds = shearCat.splitBounds();
    for(size_t i=0;i<sectionBounds.size();++i) {
        dbg<<"Starting section "<<(i+1)<<"/"<<sectionBounds.size()<<std::endl;
        dbg<<sectionBounds[i]<<std::endl;
        if (shouldOutputDots) {
            std::cerr<<"Starting section ";
            std::cerr<<(i+1)<<"/"<<sectionBounds.size()<<std::endl;
        }
        // Load the pixel information for each galaxy in the section.
        int nPix = shearCat.getPixels(sectionBounds[i]);
        if (shouldOutputDots)
            std::cerr<<nPix<<" galaxies in this section.\n";

        if (isTiming) {
            gettimeofday(&tp,0);
            t2 = tp.tv_sec + tp.tv_usec/1.e6;
            std::cout<<"Time: GetPixels = "<<t2-t1<<std::endl;
            loadTime += t2-t1;
            t1 = t2;
        }

        // Measure the shears.
        long nShear1 = shearCat.measureMultiShears(sectionBounds[i],log);

        nShear += nShear1;
        dbg<<"After MeasureShears: nshear = "<<nShear1<<"  "<<nShear<<std::endl;
        if (shouldOutputDots) {
            std::cerr<<std::endl;
            std::cerr<<nShear1<<" successful shear measurements ";
            std::cerr<<"(total so far: "<<nShear<<")\n";
        }

        if (isTiming) {
            gettimeofday(&tp,0);
            t2 = tp.tv_sec + tp.tv_usec/1.e6;
            std::cout<<"Time: MeasureMultiShears = "<<t2-t1<<std::endl;
            calcTime += t2-t1;
            t1 = t2;
        }

#ifdef VG
        dbg<<"Valgrind Leak Check:\n";
        VALGRIND_DO_LEAK_CHECK;
#endif
    }

    if (isTiming) {
        std::cout<<"Time: Total load time = "<<loadTime<<std::endl;
        std::cout<<"Time: Total calc time = "<<calcTime<<std::endl;
    }

    dbg<<"Done: "<<log._nsGamma<<" successful shear measurements, ";
    dbg<<(log._nGoodIn-log._nsGamma)<<" unsuccessful, ";
    dbg<<(log._nGals-log._nGoodIn)<<" with input flags\n";
    //log._nGood = std::count(
        //shearCat.getFlagsList().begin(),shearCat.getFlagsList().end(),0);
    nGood = 0;
    for(size_t i=0;i<shearCat.getFlagsList().size();++i)
        if (shearCat.getFlagsList()[i] == 0) ++nGood;
    log._nGood = nGood;
    dbg<<log._nGood<<" successful measurements with no measurement flags.\n";

    if (shouldOutputDots) {
        std::cerr
            <<std::endl
            <<"Success rate: "<<log._nsGamma<<"/"<<log._nGoodIn
            <<"  # with no flags: "<<log._nGood
            <<std::endl;
    }

    // Write the catalog.
    shearCat.write();
    dbg<<"After Write\n";

    if (isTiming) {
        gettimeofday(&tp,0);
        t2 = tp.tv_sec + tp.tv_usec/1.e6;
        std::cout<<"Time: Write = "<<t2-t1<<std::endl;
        t1 = t2;
    }

    if (nShear == 0) {
        throw ProcessingException(
            "No successful shear measurements");
    }

    xdbg<<"Log: \n"<<log<<std::endl;
}

int main(int argc, char **argv) try 
{
    ConfigFile params;
    if (basicSetup(argc,argv,params,"multishear")) return EXIT_FAILURE;

    // Setup Log
    std::string logFile = ""; // Default is to stdout
    if (params.keyExists("log_file") || params.keyExists("log_ext")) 
        logFile = makeName(params,"log",false,false);
    std::string multishear_file = makeName(params,"multishear",false,false);
    std::auto_ptr<ShearLog> log (
        new ShearLog(params,logFile,multishear_file)); 

    try {
        doMeasureMultiShear(params,*log);
    }
#if 0
    // Change to 1 to let gdb see where the program bombed out.
    catch(int) {}
#else
    CATCHALL;
#endif

    if (dbgout && dbgout != &std::cout) {delete dbgout; dbgout=0;}
    return EXIT_SUCCESS; // = 0 typically.  Defined in <cstdlib>
} catch (std::exception& e) {
    std::cerr<<"Fatal error: Caught \n"<<e.what()<<std::endl;
    std::cout<<"STATUS5BEG Fatal error: "<<e.what()<<" STATUS5END\n";
    return EXIT_FAILURE;
} catch (...) {
    std::cerr<<"Fatal error: Caught an exception.\n";
    std::cout<<"STATUS5BEG Fatal error: unknown exception STATUS5END\n";
    return EXIT_FAILURE;
}
