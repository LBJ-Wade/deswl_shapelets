
#include "MultiShearCatalog.h"
#include "ConfigFile.h"
#include "Params.h"
#include "Log.h"
#include "MeasureShearAlgo.h"

int MultiShearCatalog::measureMEDS(const MEDSFile& meds, ShearLog& log)
{
    dbg<<"Start MeasureMEDS\n";
    int ngals = size();
    dbg<<"ngals = "<<ngals<<std::endl;

    // Read some needed parameters
    bool output_dots = _params.read("output_dots",false);
#ifdef _OPENMP
    bool des_qa = _params.read("des_qa",false); 
#endif

    int nSuccess = 0;

#ifdef ENDAT
    ngals = ENDAT;
#endif

    // Main loop to measure shears
#ifdef _OPENMP
#pragma omp parallel 
    {
        try {
#endif
            ShearLog log1(_params); // just for this thread
            log1.noWriteLog();
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
            for(int i=0;i<ngals;++i) {
                if (_flags[i]) continue;
#ifdef STARTAT
                if (i < STARTAT) continue;
#endif
#ifdef SINGLEGAL
                if (i < SINGLEGAL) continue;
                if (i > SINGLEGAL) break;
#endif

                if (output_dots) {
#ifdef _OPENMP
#pragma omp critical (output)
#endif
                    {
                        std::cerr<<"."; std::cerr.flush(); 
                    }
                }
                dbg<<"galaxy "<<i<<":\n";
                dbg<<"id = "<<_id[i]<<std::endl;
                dbg<<"chippos = "<<_chippos[i]<<std::endl;
                dbg<<"skypos = "<<_skypos[i]<<std::endl;

                // Get the pixels and PSFs from the MEDS file.
                std::vector<PixelList> pix_list;
                meds.getPixels(pix_list,i,_flags[i]);
                std::vector<BVec> psf_list;
                meds.getPSFs(psf_list,i,_flags[i]);
                dbg<<"Using "<<pix_list.size()<<" epochs\n";
                Assert(pix_list.size() == psf_list.size());

#if 1
                MeasureSingleShear(
                    // Input data:
                    pix_list, psf_list,
                    // Parameters:
                    _meas_galorder[i], _params,
                    // Log information
                    log1,
                    // Ouput values:
                    _shape[i], _shear[i], _cov[i], _nu[i], _flags[i]);
#else
                _meas_galorder[i] = 0;
                _shear[i] = std::complex<double>(0.1,0.2);
                _cov[i] << 1., 0., 0., 1.;
                _shape[i].vec().setZero();
                _nu[i] = 10.;
#endif

                if (!_flags[i]) {
                    dbg<<"Successful shear measurements: \n";
                    dbg<<"shape = "<<_shape[i]<<std::endl;
                    dbg<<"shear = "<<_shear[i]<<std::endl;
                    dbg<<"cov = "<<_cov[i]<<std::endl;
                    dbg<<"nu = "<<_nu[i]<<std::endl;
                    ++nSuccess;
                } else {
                    dbg<<"Unsuccessful shear measurement\n"; 
                    dbg<<"flag = "<<_flags[i]<<std::endl;
                }
            }
#ifdef _OPENMP
#pragma omp critical (add_log)
#endif
            {
                log += log1;
            }
#ifdef _OPENMP
        } catch (std::exception& e) {
            // This isn't supposed to happen.
            if (des_qa) {
                std::cerr<<"STATUS5BEG Caught error in parallel region STATUS5END\n";
            } 
            std::cerr<<"Caught "<<e.what()<<std::endl;
            std::cerr<<"Caught error in parallel region.  Aborting.\n";
            exit(1);
        } catch (...) {
            if (des_qa) {
                std::cerr<<"STATUS5BEG Caught error in parallel region STATUS5END\n";
            }
            std::cerr<<"Caught error in parallel region.  Aborting.\n";
            exit(1);
        }
    }
#endif

    dbg<<nSuccess<<" successful shear measurements in this pass.\n";
    dbg<<log._ns_gamma<<" successful shear measurements so far.\n";
    xdbg<<log<<std::endl;

    return nSuccess;
}

