
#include <sys/time.h>
#include "fitsio.h"
#include "BVec.h"
#include "InputCatalog.h"
#include "FittedPsf.h"
#include "Image.h"
#include "BasicSetup.h"

int main(int argc, char **argv) try 
{
    //const double PI = 3.141592653589793;

    ConfigFile params;
    if (BasicSetup(argc,argv,params,"g10star")) return EXIT_FAILURE;

    // Setup Log
    std::string log_file = ""; // Default is to stdout
    if (params.keyExists("log_file") || params.keyExists("log_ext"))
        log_file = MakeName(params,"log",false,false);
    std::string psf_file=MakeName(params,"psf",false,false);
    std::auto_ptr<PsfLog> log (
        new PsfLog(params,log_file,psf_file));

    try {
        bool isTiming = params.read("timing",false);
        timeval tp;
        double t1=0.,t2=0.;

        if (isTiming) {
            gettimeofday(&tp,0);
            t1 = tp.tv_sec + tp.tv_usec/1.e6;
        }

        // Read distortion function
        Transformation trans(params);

        if (isTiming) {
            gettimeofday(&tp,0);
            t2 = tp.tv_sec + tp.tv_usec/1.e6;
            std::cout<<"Time: Read Transformation = "<<t2-t1<<std::endl;
            t1 = t2;
        }

        // Read input catalog
        // Note: this is the target list for where to interpolate.
        InputCatalog incat(params);
        incat.read();

        if (isTiming) {
            gettimeofday(&tp,0);
            t2 = tp.tv_sec + tp.tv_usec/1.e6;
            std::cout<<"Time: Read InputCatalog = "<<t2-t1<<std::endl;
            t1 = t2;
        }

        // Read the fitted psf file
        FittedPsf fitpsf(params);
        fitpsf.read();

        if (isTiming) {
            gettimeofday(&tp,0);
            t2 = tp.tv_sec + tp.tv_usec/1.e6;
            std::cout<<"Time: Read FittedPSF = "<<t2-t1<<std::endl;
            t1 = t2;
        }

        const int ntarget = incat.size();
        const int pwidth = 48;
        const int pheight = 48;
        const double xOffset = params.read("cat_x_offset",0.);
        const double yOffset = params.read("cat_y_offset",0.);
        const Position cen(24.,24.);
        std::vector<Image<double> > im_list(
            ntarget, Image<double>(pwidth,pheight));

#ifdef _OPENMP
#pragma omp for
#endif
        for(int n=0; n<ntarget; ++n) {
            Position pos = incat.getPos(n);
            BVec b = fitpsf(pos);
            b.makeImage(im_list[n],cen,0.,trans,xOffset,yOffset);
        }

        // Write data to fits file:
        // (This is largely copied from Tom's example C code.)
        std::string output_file = 
            std::string("!") + MakeFitsName(params, "star_image");
        // ! at the start means overwrite existing file if any.
        fitsfile *fptr;
        int fits_err=0;
        fits_create_file(&fptr, output_file.c_str(), &fits_err);
        if (fits_err != 0) {
            fits_report_error(stderr,fits_err);
            throw WriteException(
                "Error opening fits file " + output_file);
        }
        dbg<<"opened output file: "<<output_file<<std::endl;

        /* (-)number of bits per pixel*/
        int bitpix = -32;

        /*define the size of the 3d data cube*/
        int anaxis = 3;
        long anaxes_out[3];
        anaxes_out[0] = pwidth;
        anaxes_out[1] = pheight;
        anaxes_out[2] = ntarget;
        int size3d = ntarget * pwidth * pheight;
        long fpixel_out[3] = {1,1,1};

        /* Need to fill the 3D data cube correctly*/
        std::vector<float> array3d(size3d);
        const int postage_size = pwidth*pheight;
        for (int n=0; n<ntarget; n++)
        {
            /* write postage stamps into 3D array */
            for (int i=0; i<postage_size; i++)
            {
                int j = n*postage_size + i;
                array3d[j] = (&im_list[n](0,0))[i];
            }
        }

        /*first create the image, bits, and dimensionsm, to be filled with data*/
        fits_create_img(fptr, bitpix, anaxis, anaxes_out, &fits_err);
        if (fits_err != 0) {
            fits_report_error(stderr,fits_err);
            throw WriteException(
                "Error creating image " + output_file);
        }

        /*then fill the define image with data from array3d*/
        if (fits_write_pix(
                fptr, TFLOAT, fpixel_out, size3d, &array3d[0], &fits_err))
        {
            fits_report_error(stderr,fits_err);
            throw WriteException(
                "Error writing postage stamps to " + output_file);
        }
        /* close the fits file */
        fits_close_file(fptr, &fits_err);
        if (fits_err != 0) {
            fits_report_error(stderr,fits_err);
            throw WriteException(
                "Error closing fits file " + output_file);
        }
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
    std::cerr<<"Fatal error: Cought an exception.\n";
    std::cout<<"STATUS5BEG Fatal error: unknown exception STATUS5END\n";
    return EXIT_FAILURE;
}
