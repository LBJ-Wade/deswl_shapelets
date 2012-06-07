
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

static void DoMeasureMultiShear(ConfigFile& params, ShearLog& log) 
{
    const double ARCSEC_PER_RAD = 206264.806247;

    bool timing = params.read("timing",false);
    timeval tp;
    double t1=0.,t2=0.;

    if (timing) {
        gettimeofday(&tp,0);
        t1 = tp.tv_sec + tp.tv_usec/1.e6;
    }

    bool output_info = params.read("output_info",true);
    bool des_qa = params.read("des_qa",true);

    // Read the coadd catalog.
    CoaddCatalog coaddcat(params);
    coaddcat.read();
    dbg<<"Made coaddcat\n";

    if (timing) {
        gettimeofday(&tp,0);
        t2 = tp.tv_sec + tp.tv_usec/1.e6;
        std::cout<<"Time: Read CoaddCatalog = "<<t2-t1<<std::endl;
        t1 = t2;
    }

    // Make the multishear catalog.
    MultiShearCatalog shearcat(coaddcat,params);
    dbg<<"Made multishearcat\n";
    dbg<<"Total bounds are "<<shearcat.getSkyBounds()<<std::endl;
    double area = shearcat.getSkyBounds().getArea();
    area /= 3600.;
    dbg<<"Raw area = "<<area<<" square arcmin\n";
    double dec = shearcat.getSkyBounds().getCenter().getY();
    double cosdec = cos(dec / ARCSEC_PER_RAD);
    area *= cosdec;
    dbg<<"Corrected area = "<<area<<" square arcmin";
    area /= 3600.;
    dbg<<" = "<<area<<" square degrees\n";

    log._ngals = shearcat.size();
    xdbg<<"ngals = "<<log._ngals<<std::endl;

    //On abe, this next line is giving a weird error and then aborting.
    //I have no idea what is going on here, so just avoid the stl count.
    //log._ngoodin = std::count(
        //shearcat.getFlagsList().begin(),shearcat.getFlagsList().end(),0);
    int ngood = 0;
    int nFlags = shearcat.getFlagsList().size();
    for(int i=0;i<nFlags;++i) if (shearcat.getFlagsList()[i] == 0) ++ngood;
    log._ngoodin = ngood;
    xdbg<<"ngood = "<<log._ngoodin<<std::endl;
    dbg<<log._ngoodin<<"/"<<log._ngals<<" galaxies with no input flags\n";

    if (timing) {
        gettimeofday(&tp,0);
        t2 = tp.tv_sec + tp.tv_usec/1.e6;
        std::cout<<"Time: Make MultiShearCatalog = "<<t2-t1<<std::endl;
        t1 = t2;
    }

    double load_time=0.;
    double calc_time=0.;

    // Break the area up into sections and load each section separately.
    // It is most efficient to load everything at once, but it can take 
    // a lot of memory.  So this bit keeps the memory requirement for each
    // section manageable.  So long as each section does a significant amount
    // of work, the extra I/O time won't be much of an issue. 
    long nshear = 0;
    std::vector<Bounds> section_bounds = shearcat.splitBounds();
    int nSection = section_bounds.size();
    int nresplit = params.read("multishear_max_resplits",1);
    for(int i=0;i<nSection;++i) {
#ifdef ENDAT
        if (i == ENDAT) break;
#endif
        dbg<<"Starting section "<<(i+1)<<"/"<<nSection<<std::endl;
        dbg<<section_bounds[i]<<std::endl;
        if (output_info) {
            std::cerr<<"Starting section ";
            std::cerr<<(i+1)<<"/"<<nSection<<std::endl;
        }
        // Load the pixel information for each galaxy in the section.
        if (!shearcat.getPixels(section_bounds[i])) {
            dbg<<"Section exceeded maximum memory usage.\n";
            if (output_info) {
                std::cerr<<"Section exceeded maximum memory usage.\n";
            }
            if (nresplit > 0) {
                --nresplit;
                std::vector<Bounds> split = section_bounds[i].quarter();
                dbg<<"Split bounds into four new bounds:\n";
                dbg<<"b1 = "<<split[0];
                dbg<<"b2 = "<<split[1];
                dbg<<"b3 = "<<split[2];
                dbg<<"b4 = "<<split[3];
                dbg<<nresplit<<" more resplits allowed.\n";
                section_bounds.insert(
                    section_bounds.end(),
                    split.begin(),split.end());
                if (output_info) {
                    std::cerr<<"Will try splitting it up and continuing.\n";
                    std::cerr<<nresplit<<" more resplits allowed.\n";
                }
                if (des_qa) {
                    std::cout
                        << "STATUS3BEG Warning: "
                        << "Memory usage exceeded maximum allowed.  "
                        << "Trying to recover by resplitting bounds "
                        << "for section "<<i<<".  STATUS3END"<<std::endl;
                }
                continue;
            } else {
                dbg<<"No more allowed resplits.\nAborting.\n";
                if (output_info) {
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
        if (output_info) {
            std::cerr<<shearcat.getNGalsWithPixels()<<
                " galaxies in this section.\n";
        }

        if (timing) {
            gettimeofday(&tp,0);
            t2 = tp.tv_sec + tp.tv_usec/1.e6;
            std::cout<<"Time: GetPixels = "<<t2-t1<<std::endl;
            load_time += t2-t1;
            t1 = t2;
        }

        // Measure the shears.
        long nshear1 = shearcat.measureMultiShears(section_bounds[i],log);

        nshear += nshear1;
        dbg<<"After MeasureShears: nshear = "<<nshear1<<"  "<<nshear<<std::endl;
        if (output_info) {
            std::cerr<<std::endl;
            std::cerr<<nshear1<<" successful shear measurements ";
            std::cerr<<"(total so far: "<<nshear<<")\n";
        }

        if (timing) {
            gettimeofday(&tp,0);
            t2 = tp.tv_sec + tp.tv_usec/1.e6;
            std::cout<<"Time: MeasureMultiShears = "<<t2-t1<<std::endl;
            calc_time += t2-t1;
            t1 = t2;
        }

#ifdef VG
        dbg<<"Valgrind Leak Check:\n";
        VALGRIND_DO_LEAK_CHECK;
#endif
    }

    if (timing) {
        std::cout<<"Time: Total load time = "<<load_time<<std::endl;
        std::cout<<"Time: Total calc time = "<<calc_time<<std::endl;
    }

    dbg<<"Done: "<<log._ns_gamma<<" successful shear measurements, ";
    dbg<<(log._ngoodin-log._ns_gamma)<<" unsuccessful, ";
    dbg<<(log._ngals-log._ngoodin)<<" with input flags\n";
    //log._ngood = std::count(
        //shearcat.getFlagsList().begin(),shearcat.getFlagsList().end(),0);
    ngood = 0;
    nFlags = shearcat.getFlagsList().size();
    for(int i=0;i<nFlags;++i) if (shearcat.getFlagsList()[i] == 0) ++ngood;
    log._ngood = ngood;
    dbg<<log._ngood<<" successful measurements with no measurement flags.\n";
    dbg<<"Breakdown of flags:\n";
    if (dbgout) PrintFlags(shearcat.getFlagsList(),*dbgout);

    if (output_info) {
        std::cerr
            <<std::endl
            <<"Success rate: "<<log._ns_gamma<<"/"<<log._ngoodin
            <<"  # with no flags: "<<log._ngood
            <<std::endl;
        std::cerr<<"Breakdown of flags:\n";
        PrintFlags(shearcat.getFlagsList(),std::cerr);
    }

    // Write the catalog.
    shearcat.write();
    dbg<<"After Write\n";

    if (timing) {
        gettimeofday(&tp,0);
        t2 = tp.tv_sec + tp.tv_usec/1.e6;
        std::cout<<"Time: Write = "<<t2-t1<<std::endl;
        t1 = t2;
    }

    if (nshear == 0) {
        throw ProcessingException(
            "No successful shear measurements");
    }

    xdbg<<"Log: \n"<<log<<std::endl;
}

int main(int argc, char **argv) try 
{
    ConfigFile params;
    if (BasicSetup(argc,argv,params,"multishear")) return EXIT_FAILURE;

    // Setup Log
    std::string log_file = ""; // Default is to stdout
    if (params.keyExists("log_file") || params.keyExists("log_ext")) 
        log_file = MakeName(params,"log",false,false);
    std::string multishear_file = MakeName(params,"multishear",false,false);
    std::auto_ptr<ShearLog> log (
        new ShearLog(params,log_file,multishear_file)); 

    try {
        DoMeasureMultiShear(params,*log);
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
