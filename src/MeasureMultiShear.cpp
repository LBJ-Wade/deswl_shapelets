
#include <valarray>
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

#include "Image.h"
#include "FittedPsf.h"
#include "Log.h"
#include "CoaddCatalog.h"
#include "MultiShearCatalog.h"
#include "BasicSetup.h"

// If desired stop after this many sections.
//#define ENDAT 2

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

    bool shouldOutputInfo = params.read("output_info",true);
    bool shouldOutputDesQa = params.read("des_qa",true);

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
    int nFlags = shearCat.getFlagsList().size();
    for(int i=0;i<nFlags;++i) if (shearCat.getFlagsList()[i] == 0) ++nGood;
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
    long nShear = 0;
    std::vector<Bounds> sectionBounds = shearCat.splitBounds();
    int nSection = sectionBounds.size();
    int nresplit = params.read("multishear_max_resplits",1);
    for(int i=0;i<nSection;++i) {
#ifdef ENDAT
        if (i == ENDAT) break;
#endif
        dbg<<"Starting section "<<(i+1)<<"/"<<nSection<<std::endl;
        dbg<<sectionBounds[i]<<std::endl;
        if (shouldOutputInfo) {
            std::cerr<<"Starting section ";
            std::cerr<<(i+1)<<"/"<<nSection<<std::endl;
        }
        // Load the pixel information for each galaxy in the section.
        if (!shearCat.getPixels(sectionBounds[i])) {
            dbg<<"Section exceeded maximum memory usage.\n";
            if (shouldOutputInfo) {
                std::cerr<<"Section exceeded maximum memory usage.\n";
            }
            if (nresplit > 0) {
                --nresplit;
                std::vector<Bounds> split = sectionBounds[i].quarter();
                dbg<<"Split bounds into four new bounds:\n";
                dbg<<"b1 = "<<split[0];
                dbg<<"b2 = "<<split[1];
                dbg<<"b3 = "<<split[2];
                dbg<<"b4 = "<<split[3];
                dbg<<nresplit<<" more resplits allowed.\n";
                sectionBounds.insert(
                    sectionBounds.end(),
                    split.begin(),split.end());
                if (shouldOutputInfo) {
                    std::cerr<<"Will try splitting it up and continuing.\n";
                    std::cerr<<nresplit<<" more resplits allowed.\n";
                }
                if (shouldOutputDesQa) {
                    std::cout
                        << "STATUS3BEG Warning: "
                        << "Memory usage exceeded maximum allowed.  "
                        << "Trying to recover by resplitting bounds "
                        << "for section "<<i<<".  STATUS3END"<<std::endl;
                }
                continue;
            } else {
                dbg<<"No more allowed resplits.\nAborting.\n";
                if (shouldOutputInfo) {
                    std::cerr<<"No more allowed resplits.\n";
                    std::cerr<<"Either reduce multishear_section_size, or\n";
                    std::cerr<<"increase max_vmem or multishear_max_resplits\n";
                    std::cerr<<"Current values = "<<
                        params["multishear_section_size"]<<" , "<<
                        params["max_vmem"]<<" , "<<
                        params["multishear_max_resplits"]<<std::endl;
                    throw ProcessingException(
                        "Memory exceeded maximum allowed, and unable to "
                        "recover by resplitting bounds");
                }
            }
        }
        if (shouldOutputInfo) {
            std::cerr<<shearCat.getNGalsWithPixels()<<
                " galaxies in this section.\n";
        }

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
        if (shouldOutputInfo) {
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
    nFlags = shearCat.getFlagsList().size();
    for(int i=0;i<nFlags;++i) if (shearCat.getFlagsList()[i] == 0) ++nGood;
    log._nGood = nGood;
    dbg<<log._nGood<<" successful measurements with no measurement flags.\n";
    dbg<<"Breakdown of flags:\n";
    if (dbgout) PrintFlags(shearCat.getFlagsList(),*dbgout);

    if (shouldOutputInfo) {
        std::cerr
            <<std::endl
            <<"Success rate: "<<log._nsGamma<<"/"<<log._nGoodIn
            <<"  # with no flags: "<<log._nGood
            <<std::endl;
        std::cerr<<"Breakdown of flags:\n";
        PrintFlags(shearCat.getFlagsList(),std::cerr);
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
