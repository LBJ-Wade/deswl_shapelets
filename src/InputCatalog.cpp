
#include <valarray>
#include <CCfits/CCfits>

#include "InputCatalog.h"
#include "StarCatalog.h"
#include "dbg.h"
#include "Name.h"
#include "Params.h"
#include "Image.h"

void InputCatalog::flagStars(const StarCatalog& starCat)
{
    const int nGals = _id.size();
    for (int i=0; i<nGals; ++i) {
        xdbg<<"i = "<<i<<", id[i] = "<<_id[i]<<std::endl;
        // It doesn't seem worth making a separate flag for this.
        // If it's considered a star, let's just admit that it's too small 
        // to bother trying to measure a shear for it.
        if (starCat.isAStar(_id[i])) {
            xdbg<<"Flag this one as a star\n";
            _flags[i] |= TOO_SMALL;
            xdbg<<i<<" is a star: flag -> "<<_flags[i]<<std::endl;
        }
    }
}


static void readGain(const std::string& file, int hdu, ConfigFile& params)
{
    if (!doesFileExist(file)) {
        throw FileNotFoundException(file);
    }
    xdbg<<"ReadGain: from fits file "<<file<<std::endl;
    CCfits::FITS fits(file, CCfits::Read);
    if (hdu > 1) fits.read(hdu-1);

    double gain, readNoise;

    std::vector<std::string> gainKey = 
        params.read<std::vector<std::string> >("image_gain_key");
    std::vector<std::string> readNoiseKey =
        params.read<std::vector<std::string> >("image_readnoise_key");

    const int nGainKeys = gainKey.size();
    const int nReadNoiseKeys = readNoiseKey.size();
    for(int k=0;k<nGainKeys;++k) {
        xdbg<<"try "<<gainKey[k]<<std::endl;
        try {
            fits.pHDU().readKey(gainKey[k],gain);
            break;
        } catch (CCfits::FitsException& e) {
            xdbg<<"Caught FitsException: \n"<<e.message()<<std::endl;
            if (k == nGainKeys-1) {
                throw ReadException(
                    "Error reading gain from Fits file "+file+
                    "\nCCfits error message: \n"+e.message());
            }
        }
    }

    for(int k=0;k<nReadNoiseKeys;++k) {
        xdbg<<"try "<<readNoiseKey[k]<<std::endl;
        try {
            fits.pHDU().readKey(readNoiseKey[k], readNoise);
            break;
        } catch (CCfits::FitsException& e) { 
            xdbg<<"Caught FitsException: \n"<<e.message()<<std::endl;
            if (k == nReadNoiseKeys-1) {
                throw ReadException(
                    "Error reading read_noise from Fits file "+file+
                    "\nCCfits error message: \n"+e.message());
            }
        }
    }

    params["image_gain"] = gain;
    params["image_readnoise"] = readNoise;
}

InputCatalog::InputCatalog(ConfigFile& params, const Image<double>* im) : 
    _params(params), _im(im),
    _noiseValue(0), _gain(0), _readNoise(0)
{
  this->determineNoiseMethod();
}

void InputCatalog::init(ConfigFile& params, const Image<double>* im)
{
  _noiseValue=0;
  _gain=0;
  _readNoise=0;
  _im=im;
  this->loadParams(params);
  this->determineNoiseMethod();
}

void InputCatalog::determineNoiseMethod() {
  // Setup noise calculation:
  std::string noiseMethod = _params.get("noise_method");
  if (noiseMethod == "VALUE") {
    _nm = VALUE;
    _noiseValue = _params.read<double>("noise");
  } else if (noiseMethod == "CATALOG") {
    _nm = CATALOG;
    Assert(_params.keyExists("cat_noise_col"));
  } else if (noiseMethod == "CATALOG_SIGMA") {
    _nm = CATALOG_SIGMA;
    Assert(_params.keyExists("cat_noise_col"));
  } else if (noiseMethod == "GAIN_VALUE") {
    _nm = GAIN_VALUE;
    _gain = _params.read<double>("image_gain");
    _readNoise = _params.read<double>("image_readnoise");
    dbg<<"gain, readnoise = "<<_gain<<"  "<<_readNoise<<std::endl;
  } else if (noiseMethod == "GAIN_FITS") {
    std::string imageName = makeName(_params,"image",true,true);
    int hdu = getHdu(_params,"image",imageName,1);
    readGain(imageName,hdu,_params);
    xdbg<<"Read gain = "<<_params["image_gain"]<<
      ", rdn = "<<_params["image_readnoise"]<<std::endl;
    _gain = _params.read<double>("image_gain");
    _readNoise = _params.read<double>("image_readnoise");
    _nm = GAIN_VALUE;
    dbg<<"gain, readnoise = "<<_gain<<"  "<<_readNoise<<std::endl;
  } else if (noiseMethod == "WEIGHT_IMAGE") {
    dbg<<"using WEIGHT_IMAGE as noise method"<<std::endl;
    _nm = WEIGHT_IMAGE;
    Assert(_params.keyExists("weight_ext") || 
        _params.keyExists("weight_file"));
  } else {
    throw ParameterException("Unknown noise method "+noiseMethod);
  }

}


void InputCatalog::read()
{
    std::string file = makeName(_params,"cat",true,true);
    dbg<< "Reading input cat from file: " << file << std::endl;

    bool isFitsIo = false;
    if (_params.keyExists("cat_io")) 
        isFitsIo = (_params["cat_io"] == "FITS");
    else if (file.find("fits") != std::string::npos) 
        isFitsIo = true;

    if (!doesFileExist(file)) {
        throw FileNotFoundException(file);
    }

    try {
        if (isFitsIo) {
            readFits(file);
        } else {
            std::string delim = "  ";
            if (_params.keyExists("cat_delim")) delim = _params["cat_delim"];
            else if (file.find("csv") != std::string::npos) delim = ",";
            readAscii(file,delim);
        }
    } catch (CCfits::FitsException& e) {
        xdbg<<"Caught FitsException: \n"<<e.message()<<std::endl;
        throw ReadException(
            "Error reading from "+file+" -- caught error\n" + e.message());
    } catch (std::exception& e) {
        xdbg<<"Caught std::exception: \n"<<e.what()<<std::endl;
        throw ReadException(
            "Error reading from "+file+" -- caught error\n" + e.what());
    } catch (...) {
        xdbg<<"Caught unknown exception: "<<std::endl;
        throw ReadException(
            "Error reading from "+file+" -- caught unknown error");
    }
    // Update the bounds:
    const int nPos = _pos.size();
    for(int i=0;i<nPos;++i) _bounds += _pos[i];

    // These are the only two fields guaranteed to be set at this point.
    Assert(_id.size() == _pos.size());
    int nRows = _id.size();

    // Fix sky if necessary
    double badSkyVal = _params.read("cat_bad_sky",-999.);
    if (_sky.size() == 0) _sky.resize(_id.size(),badSkyVal);
    if (std::find(_sky.begin(), _sky.end(), badSkyVal) != _sky.end()) {
        double globSky = 0.;
        if (_params.keyExists("cat_global_sky")) {
            globSky = _params["cat_global_sky"];
        } else {
            if (_im) {
                globSky = _im->median();
            } else {
                Image<double> im1(_params);
                globSky = im1.median();
            }
            dbg<<"Found global sky from image median value.\n";
        }
        dbg<<"Set global value of sky to "<<globSky<<std::endl;
        for(int i=0;i<nRows;++i) {
            if (_sky[i] == badSkyVal) _sky[i] = globSky;
        }
    }
    Assert(_sky.size() == _pos.size());

    // MJ: <100 have basically no chance to find the stars
    if (_params.read("des_qa",false)) {
        int minRows = _params.read("cat_nrows",0);
        if (nRows <= minRows) {
            try {
                std::string name = makeName(_params,"cat",true,false);
                std::cout<<"STATUS3BEG Warning: Input catalog only has "
                    <<nRows<<" rows for Name="<<name
                    <<". STATUS3END"<<std::endl;
            } catch (AssertFailureException& ) {
                std::cout<<"STATUS3BEG Warning: Input catalog only has "
                    <<nRows<<" rows. STATUS3END"<<std::endl;
            }
        }
    }

    // Update noise calculation if necessary
    Assert(_nm == VALUE || _nm == CATALOG || _nm == CATALOG_SIGMA ||
           _nm == GAIN_VALUE || _nm == WEIGHT_IMAGE);
    if (_nm == VALUE) {
        _noise.resize(nRows,0);
        for(int i=0;i<nRows;++i) {
            _noise[i] = _noiseValue;
        }
        xdbg<<"Set all noise to "<<_noiseValue<<std::endl;
    } else if (_nm == CATALOG) {
        Assert(_noise.size() == _id.size());
    } else if (_nm == CATALOG_SIGMA) {
        Assert(_noise.size() == _id.size());
        for(int i=0;i<nRows;++i) {
            _noise[i] = _noise[i]*_noise[i];
        }
        xdbg<<"Squared noise values from catalog\n";
    } else if (_nm == GAIN_VALUE) {
        _noise.resize(nRows,0);
        double extraSky=_params.read("image_extra_sky",0.);
        xdbg<<"Calculate noise from sky, gain, readnoise\n";
        for(int i=0;i<nRows;++i) {
            _noise[i] = (_sky[i]+extraSky)/_gain + _readNoise;
        }
    } else if (_nm == WEIGHT_IMAGE) {
        // Then we don't need the noise vector, but it is easier
        // to just fill it with 0's.
        _noise.resize(nRows,0);
        xdbg<<"Set noise values to 0, since using weight_im, not noise vector.\n";
    }

    // Convert input flags into our flag schema
    if (_flags.size() == 0) {
        dbg<<"No flags read in -- starting all with 0\n";
        _flags.resize(_id.size(),0);
    } else {
        long ignoreFlags = ~0L;
        dbg<<std::hex<<std::showbase;
        if (_params.keyExists("cat_ignore_flags")) {
            ignoreFlags = _params["cat_ignore_flags"];
            dbg<<"Using ignore flag parameter = "<<ignoreFlags<<std::endl;
        } else if (_params.keyExists("cat_ok_flags")) {
            ignoreFlags = _params["cat_ok_flags"];
            dbg<<"Using ok flag parameter = "<<ignoreFlags<<std::endl;
            ignoreFlags = ~ignoreFlags;
            dbg<<"ignore flag = "<<ignoreFlags<<std::endl;
        } else {
            dbg<<"No ok or ignore parameter: use ignore flag = "<<
                ignoreFlags<<std::endl;
        }
        Assert(_flags.size() == _id.size());
        dbg<<"INPUT_FLAG: "<<INPUT_FLAG<<"\n";
        for(int i=0;i<nRows;++i) {
            _flags[i] = (_flags[i] & ignoreFlags) ? INPUT_FLAG : 0;
        }
        dbg<<std::dec<<std::noshowbase;
    }

    // Calculate the skyBounds using only the good objects:
    if (_skyPos.size() == _id.size()) {
        Bounds skyBounds2; // In degrees, rather than arcsec
        for(int i=0;i<nRows;++i) if (!_flags[i]) {
            _skyBounds += _skyPos[i];
            Position temp = _skyPos[i];
            temp /= 3600.;
            skyBounds2 += temp;
        }
        dbg<<"skyBounds = "<<_skyBounds<<std::endl;
        dbg<<"in degrees: "<<skyBounds2<<std::endl;
    }

    // At this point the vectors that are guaranteed to be filled are:
    // id, pos, sky, noise, flags
    Assert(_pos.size() == _id.size());
    Assert(_sky.size() == _id.size());
    Assert(_noise.size() == _id.size());
    Assert(_flags.size() == _id.size());

    dbg<<"Done Read InputCatalog\n";
}

void InputCatalog::readFits(std::string file) {
    int hdu = getHdu(_params,"cat",file,1);
    readFits(file, hdu);
}
void InputCatalog::readFits(std::string file, int hdu)
{
    dbg<< "Reading cat from FITS file: " << file << std::endl;

    dbg<<"Opening InputCatalog file "<<file<<" at hdu "<<hdu<<std::endl;
    CCfits::FITS fits(file, CCfits::Read);
    if (hdu > 1) fits.read(hdu-1);

    CCfits::ExtHDU& table=fits.extension(hdu-1);

    long nRows=table.rows();

    dbg<<"  nrows = "<<nRows<<std::endl;
    if (nRows <= 0) {
        throw ReadException(
            "InputCatalog found to have 0 rows.  Must have > 0 rows.");
    }

    // Read each column in turn:
    //
    long start=1;
    long end=nRows;

    dbg<<"Reading columns"<<std::endl;


    if (_params.keyExists("cat_id_col")) {
        std::string idCol = _params["cat_id_col"];
        dbg<<"  "<<idCol<<std::endl;
        table.column(idCol).read(_id, start, end);
    } else {
        _id.resize(nRows);
        for (int i=0;i<nRows;++i) _id[i] = i+1;
    }
    Assert(int(_id.size()) == nRows);

    // Position (on chip)
    std::vector<double> posX;
    std::vector<double> posY;
    _pos.resize(nRows);
    std::string xCol=_params.get("cat_x_col");
    std::string yCol=_params.get("cat_y_col");
    dbg<<"  "<<xCol<<std::endl;
    table.column(xCol).read(posX, start, end);
    dbg<<"  "<<yCol<<std::endl;
    table.column(yCol).read(posY, start, end);

    for (long i=0; i< nRows; ++i) {
        _pos[i] = Position(posX[i], posY[i]);
    }

    // Local sky
    if (_params.keyExists("cat_sky_col")) {
        std::string skyCol=_params["cat_sky_col"];
        dbg<<"  "<<skyCol<<std::endl;
        table.column(skyCol).read(_sky, start, end);
    }

    // Magnitude
    if (_params.keyExists("cat_mag_col")) {
        std::string magCol=_params["cat_mag_col"];
        dbg<<"  "<<magCol<<std::endl;
        table.column(magCol).read(_mag, start, end);
    }

    // Magnitude Error
    if (_params.keyExists("cat_mag_err_col")) {
        std::string magErrCol=_params["cat_mag_err_col"];
        dbg<<"  "<<magErrCol<<std::endl;
        table.column(magErrCol).read(_magErr, start, end);
    }

    // Star Galaxy separation
    if (_params.keyExists("cat_sg_col")) {
        std::string sgCol=_params["cat_sg_col"];
        dbg<<"  "<<sgCol<<std::endl;
        table.column(sgCol).read(_sg, start, end);
    } 

    // Size
    if (_params.keyExists("cat_size_col")) {
        std::string sizeCol=_params["cat_size_col"];
        dbg<<"  "<<sizeCol<<std::endl;
        table.column(sizeCol).read(_objSize, start, end);
        if (_params.keyExists("cat_size2_col")) {
            std::string size2Col=_params["cat_size2_col"];
            dbg<<"  "<<size2Col<<std::endl;
            std::vector<double> objSize2(_objSize.size());
            table.column(size2Col).read(objSize2, start, end);
            Assert(_objSize.size() == _id.size());
            Assert(objSize2.size() == _id.size());
            for(int i=0;i<nRows;++i) _objSize[i] += objSize2[i];
        }
    }

    // Error flags
    if (_params.keyExists("cat_flag_col")) {
        std::string flagCol=_params["cat_flag_col"];
        dbg<<"  "<<flagCol<<std::endl;
        table.column(flagCol).read(_flags, start, end);
    }


    // RA
    if (_params.keyExists("cat_ra_col")) {
        if (!_params.keyExists("cat_dec_col"))
            throw ParameterException("cat_ra_col, but no cat_dec_col");

        std::string raCol=_params["cat_ra_col"];
        std::string declCol=_params["cat_dec_col"];

        std::vector<double> ra;
        std::vector<double> decl;
        dbg<<"  "<<raCol<<std::endl;
        table.column(raCol).read(ra, start, end);
        dbg<<"  "<<declCol<<std::endl;
        table.column(declCol).read(decl, start, end);
        
        _skyPos.resize(nRows);
        for(long i=0;i<nRows;++i) {
            _skyPos[i] = Position(ra[i],decl[i]);
            
            // The convention for Position is to use arcsec for everything.
            // ra and dec come in as degrees.  So wee need to convert to arcsec.
            _skyPos[i] *= 3600.;  // deg -> arcsec
        }
    } else if (_params.keyExists("cat_dec_col")) {
        throw ParameterException("cat_dec_col, but no cat_ra_col");
    }

    // Noise
    if (_params.keyExists("cat_noise_col")) {
        std::string noiseCol=_params["cat_noise_col"];
        dbg<<"  "<<noiseCol<<std::endl;
        table.column(noiseCol).read(_noise, start, end);
    }
}

// Helper functino to convert ASCII line into vector of tokens
static void getTokens(
    std::string delim, std::string line,
    std::vector<ConvertibleString>& tokens)
{
    std::istringstream lineIn(line);
    if (delim == "  ") {
        std::string temp;
        while (lineIn >> temp) tokens.push_back(temp);
    } else {
        if (delim.size() > 1) {
            // getline only works with single character delimiters.
            // Since I don't really expect a multicharacter delimiter to
            // be used ever, I'm just going to throw an exception here 
            // if we do need it, and I can write the workaround then.
            throw ParameterException(
                "ReadAscii delimiter must be a single character");
        }
        char d = delim[0];
        std::string temp;
        while (getline(lineIn,temp,d)) {
            tokens.push_back(temp);
        }
    }
}

void InputCatalog::readAscii(std::string file, std::string delim)
{
    std::ifstream catIn(file.c_str());
    if (!catIn) {
        throw ReadException("Error opening input catalog file"+file);
    }
    xdbg<<"Opened catalog "<<file<<std::endl;

    // x,y is required
    std::string line;
    int xCol = _params.read<int>("cat_x_col");
    int yCol = _params.read<int>("cat_y_col");

    // Set column numbers for optional columns
    int idCol = _params.read("cat_id_col",0);

    int magCol = _params.read("cat_mag_col",0);
    int magErrCol = _params.read("cat_mag_err_col",0);

    int sgCol = _params.read("cat_sg_col",0);

    int sizeCol = _params.read("cat_size_col",0);
    int size2Col = _params.read("cat_size2_col",0);

    int flagCol = _params.read("cat_flag_col",0);

    int skyCol = _params.read("cat_sky_col",0);

    int raCol = _params.read("cat_ra_col",0);
    int declCol = _params.read("cat_dec_col",0);

    int noiseCol = _params.read("cat_noise_col",0);

    // Set up allowed comment markers
    std::vector<std::string> commentMarker = 
        _params.read("cat_comment_marker",std::vector<std::string>(1,"#"));

    // Keep running id value when idCol = 0
    int idVal = 0;

    // Read each line from catalog file
    while (getline(catIn,line)) {

        // Skip if this is a comment.
        bool shouldSkip = false;
        const int nCommentMarkers = commentMarker.size();
        for(int k=0;k<nCommentMarkers;++k) {
            if (std::string(line,0,commentMarker[k].size()) == 
                commentMarker[k]) {
                shouldSkip = true;
            }
        }
        if (shouldSkip) continue;

        // Convert line into vector of tokens
        std::vector<ConvertibleString> tokens;
        getTokens(delim,line,tokens);

        // ID
        if (idCol) {
            Assert(idCol <= int(tokens.size()));
            idVal = tokens[idCol-1];
        } else {
            // if not reading id, then just increment to get sequential values
            ++idVal;
        }
        _id.push_back(idVal);

        // Position
        Assert(xCol <= int(tokens.size()));
        Assert(yCol <= int(tokens.size()));
        double x = tokens[xCol-1];
        double y = tokens[yCol-1];
        _pos.push_back(Position(x,y));

        // Sky
        if (skyCol) {
            // Note: if sky not read in, then sky.size() is still 0
            // This is indicator to update with global given or median 
            // value later.
            double skyVal = 0.;
            Assert(skyCol <= int(tokens.size()));
            skyVal = tokens[skyCol-1];
            _sky.push_back(skyVal);
        } 

        // Magnitude
        if (magCol) {
            double magVal = 0.;
            Assert(magCol <= int(tokens.size()));
            magVal = tokens[magCol-1];
            _mag.push_back(magVal);
        } 

        // Magnitude error
        if (magErrCol) {
            double magErrVal = 0.;
            Assert(magErrCol <= int(tokens.size()));
            magErrVal = tokens[magErrCol-1];
            _magErr.push_back(magErrVal);
        } 

        // Star-galaxy
        if (sgCol) {
            double sgVal = 0.;
            Assert(sgCol <= int(tokens.size()));
            sgVal = tokens[sgCol-1];
            _sg.push_back(sgVal);
        } 

        // Size
        if (sizeCol) {
            double sizeVal = 0.;
            Assert(sizeCol <= int(tokens.size()));
            sizeVal = tokens[sizeCol-1];
            if (size2Col) {
                Assert(size2Col <= int(tokens.size()));
                sizeVal += double(tokens[size2Col-1]);
            } 
            _objSize.push_back(sizeVal);
        } 

        // Flags
        if (flagCol) {
            long flagVal = 0;
            Assert(flagCol <= int(tokens.size()));
            flagVal = tokens[flagCol-1];
            _flags.push_back(flagVal);
        } 

        // RA
        if (raCol) {
            if (!declCol) 
                throw ParameterException("cat_ra_col, but no cat_dec_col");
            double raVal = 0.;
            double declVal = 0.;
            Assert(raCol <= int(tokens.size()));
            Assert(declCol <= int(tokens.size()));
            raVal = tokens[raCol-1];
            declVal = tokens[declCol-1];
            Position skyPosVal(raVal,declVal);
            skyPosVal *= 3600.; // deg -> arcsec
            _skyPos.push_back(skyPosVal);
        } else if (declCol) {
            throw ParameterException("cat_dec_col, but no cat_ra_col");
        }

        // Noise
        if (noiseCol) {
            double noiseVal = 0.;
            Assert(noiseCol <= int(tokens.size()));
            noiseVal = tokens[noiseCol-1];
            _noise.push_back(noiseVal);
        }
    }
}

void InputCatalog::printall(int i) {
  if (int(_id.size()) > i) {
    std::cout<<"  InputCatalog::printall id["<<i<<"]: "<<_id[i]<<"\n";
  } else {
    std::cout<<"  InputCatalog::printall index "<<i<<" is larger than id.size\n";
  }

  if (int(_pos.size()) > i) {
    std::cout<<"  InputCatalog::printall pos["<<i<<"]: "<<_pos[i]<<"\n";
  } else {
    std::cout<<"  InputCatalog::printall index "<<i<<" is larger than pos.size\n";
  }

  if (int(_sky.size()) > i) {
    std::cout<<"  InputCatalog::printall sky["<<i<<"]: "<<_sky[i]<<"\n";
  } else {
    std::cout<<"  InputCatalog::printall index "<<i<<" is larger than sky.size\n";
  }

  if (int(_mag.size()) > i) {
    std::cout<<"  InputCatalog::printall mag["<<i<<"]: "<<_mag[i]<<"\n";
  } else {
    std::cout<<"  InputCatalog::printall index "<<i<<" is larger than mag.size\n";
  }

  if (int(_magErr.size()) > i) {
    std::cout<<"  InputCatalog::printall magErr["<<i<<"]: "<<_magErr[i]<<"\n";
  } else {
    std::cout<<"  InputCatalog::printall index "<<i<<" is larger than magErr.size\n";
  }

  if (int(_sg.size()) > i) {
    std::cout<<"  InputCatalog::printall sg["<<i<<"]: "<<_sg[i]<<"\n";
  } else {
    std::cout<<"  InputCatalog::printall index "<<i<<" is larger than sg.size\n";
  }

  if (int(_objSize.size()) > i) {
    std::cout<<"  InputCatalog::printall objSize["<<i<<"]: "<<_objSize[i]<<"\n";
  } else {
    std::cout<<"  InputCatalog::printall index "<<i<<" is larger than objSize.size\n";
  }

  if (int(_flags.size()) > i) {
    std::cout<<"  InputCatalog::printall flags["<<i<<"]: "<<this->_flags[i]<<"\n";
  } else {
    std::cout<<"  InputCatalog::printall index "<<i<<" is larger than flags.size\n";
  }

  if (int(_noise.size()) > i) {
    std::cout<<"  InputCatalog::printall noise["<<i<<"]: "<<this->_noise[i]<<"\n";
  } else {
    std::cout<<"  InputCatalog::printall index "<<i<<" is larger than noise.size\n";
  }
}
