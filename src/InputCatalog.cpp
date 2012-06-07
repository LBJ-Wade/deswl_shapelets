
#include <valarray>
#include <CCfits/CCfits>

#include "dbg.h"
#include "InputCatalog.h"
#include "StarCatalog.h"
#include "Name.h"
#include "Params.h"
#include "Image.h"

void InputCatalog::flagStars(const StarCatalog& starcat)
{
    const int ngals = _id.size();
    for (int i=0; i<ngals; ++i) {
        xdbg<<"i = "<<i<<", id[i] = "<<_id[i]<<std::endl;
        // It doesn't seem worth making a separate flag for this.
        // If it's considered a star, let's just admit that it's too small 
        // to bother trying to measure a shear for it.
        if (starcat.isStar(_id[i])) {
            xdbg<<"Flag this one as a star\n";
            _flags[i] |= TOO_SMALL;
            xdbg<<i<<" is a star: flag -> "<<_flags[i]<<std::endl;
        }
    }
}

static void readGain(const std::string& file, int hdu, ConfigFile& params)
{
    if (!DoesFileExist(file)) {
        throw FileNotFoundException(file);
    }
    xdbg<<"ReadGain: from fits file "<<file<<std::endl;
    CCfits::FITS fits(file, CCfits::Read);
    if (hdu > 1) fits.read(hdu-1);

    double gain, read_noise;

    std::vector<std::string> gain_key = 
        params.read<std::vector<std::string> >("image_gain_key");
    std::vector<std::string> read_noise_key =
        params.read<std::vector<std::string> >("image_readnoise_key");

    const int ngain_keys = gain_key.size();
    const int nread_noise_keys = read_noise_key.size();
    for(int k=0;k<ngain_keys;++k) {
        xdbg<<"try "<<gain_key[k]<<std::endl;
        try {
            fits.pHDU().readKey(gain_key[k],gain);
            break;
        } catch (CCfits::FitsException& e) {
            xdbg<<"Caught FitsException: \n"<<e.message()<<std::endl;
            if (k == ngain_keys-1) {
                throw ReadException(
                    "Error reading gain from Fits file "+file+
                    "\nCCfits error message: \n"+e.message());
            }
        }
    }

    for(int k=0;k<nread_noise_keys;++k) {
        xdbg<<"try "<<read_noise_key[k]<<std::endl;
        try {
            fits.pHDU().readKey(read_noise_key[k], read_noise);
            break;
        } catch (CCfits::FitsException& e) { 
            xdbg<<"Caught FitsException: \n"<<e.message()<<std::endl;
            if (k == nread_noise_keys-1) {
                throw ReadException(
                    "Error reading read_noise from Fits file "+file+
                    "\nCCfits error message: \n"+e.message());
            }
        }
    }

    params["image_gain"] = gain;
    params["image_readnoise"] = read_noise;
}

InputCatalog::InputCatalog(ConfigFile& params, const Image<double>* im) : 
    _params(params), _im(im),
    _noise_value(0), _gain(0), _read_noise(0)
{
    this->determineNoiseMethod();
}

void InputCatalog::init(ConfigFile& params, const Image<double>* im)
{
    _noise_value=0;
    _gain=0;
    _read_noise=0;
    _im=im;
    this->loadParams(params);
    this->determineNoiseMethod();
}

void InputCatalog::determineNoiseMethod() {
    // Setup noise calculation:
    std::string noiseMethod = _params.get("noise_method");
    if (noiseMethod == "VALUE") {
        _nm = VALUE;
        _noise_value = _params.read<double>("noise");
    } else if (noiseMethod == "CATALOG") {
        _nm = CATALOG;
        Assert(_params.keyExists("cat_noise_col"));
    } else if (noiseMethod == "CATALOG_SIGMA") {
        _nm = CATALOG_SIGMA;
        Assert(_params.keyExists("cat_noise_col"));
    } else if (noiseMethod == "GAIN_VALUE") {
        _nm = GAIN_VALUE;
        _gain = _params.read<double>("image_gain");
        _read_noise = _params.read<double>("image_readnoise");
        dbg<<"gain, readnoise = "<<_gain<<"  "<<_read_noise<<std::endl;
    } else if (noiseMethod == "GAIN_FITS") {
        std::string imageName = MakeName(_params,"image",true,true);
        int hdu = GetHdu(_params,"image",imageName,1);
        readGain(imageName,hdu,_params);
        xdbg<<"Read gain = "<<_params["image_gain"]<<
            ", rdn = "<<_params["image_readnoise"]<<std::endl;
        _gain = _params.read<double>("image_gain");
        _read_noise = _params.read<double>("image_readnoise");
        _nm = GAIN_VALUE;
        dbg<<"gain, readnoise = "<<_gain<<"  "<<_read_noise<<std::endl;
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
    std::string file = MakeName(_params,"cat",true,true);
    dbg<< "Reading input cat from file: " << file << std::endl;

    bool is_fitsio = false;
    if (_params.keyExists("cat_io")) 
        is_fitsio = (_params["cat_io"] == "FITS");
    else if (file.find("fits") != std::string::npos) 
        is_fitsio = true;

    if (!DoesFileExist(file)) {
        throw FileNotFoundException(file);
    }

    try {
        if (is_fitsio) {
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

    // These are the only two fields guaranteed to be set at this point.
    Assert(_id.size() == _pos.size());
    int nrows = _id.size();

    // Fix sky if necessary
    double bad_sky_val = _params.read("cat_bad_sky",-999.);
    if (_sky.size() == 0) _sky.resize(_id.size(),bad_sky_val);
    if (std::find(_sky.begin(), _sky.end(), bad_sky_val) != _sky.end()) {
        double glob_sky = 0.;
        if (_params.keyExists("cat_global_sky")) {
            glob_sky = _params["cat_global_sky"];
        } else {
            if (_im) {
                glob_sky = _im->median();
            } else {
                Image<double> im1(_params);
                glob_sky = im1.median();
            }
            dbg<<"Found global sky from image median value.\n";
        }
        dbg<<"Set global value of sky to "<<glob_sky<<std::endl;
        for(int i=0;i<nrows;++i) {
            if (_sky[i] == bad_sky_val) _sky[i] = glob_sky;
        }
    }
    Assert(_sky.size() == _pos.size());

    // MJ: <100 have basically no chance to find the stars
    if (_params.read("des_qa",false)) {
        int minrows = _params.read("cat_nrows",0);
        if (nrows <= minrows) {
            try {
                std::string name = MakeName(_params,"cat",true,false);
                std::cout<<"STATUS3BEG Warning: Input catalog only has "
                    <<nrows<<" rows for Name="<<name
                    <<". STATUS3END"<<std::endl;
            } catch (AssertFailureException& ) {
                std::cout<<"STATUS3BEG Warning: Input catalog only has "
                    <<nrows<<" rows. STATUS3END"<<std::endl;
            }
        }
    }

    // Update noise calculation if necessary
    Assert(_nm == VALUE || _nm == CATALOG || _nm == CATALOG_SIGMA ||
           _nm == GAIN_VALUE || _nm == WEIGHT_IMAGE);
    if (_nm == VALUE) {
        _noise.resize(nrows,0);
        for(int i=0;i<nrows;++i) {
            _noise[i] = _noise_value;
        }
        xdbg<<"Set all noise to "<<_noise_value<<std::endl;
    } else if (_nm == CATALOG) {
        Assert(_noise.size() == _id.size());
    } else if (_nm == CATALOG_SIGMA) {
        Assert(_noise.size() == _id.size());
        for(int i=0;i<nrows;++i) {
            _noise[i] = _noise[i]*_noise[i];
        }
        xdbg<<"Squared noise values from catalog\n";
    } else if (_nm == GAIN_VALUE) {
        _noise.resize(nrows,0);
        double extra_sky=_params.read("image_extra_sky",0.);
        xdbg<<"Calculate noise from sky, gain, readnoise\n";
        for(int i=0;i<nrows;++i) {
            _noise[i] = (_sky[i]+extra_sky)/_gain + _read_noise;
        }
    } else if (_nm == WEIGHT_IMAGE) {
        // Then we don't need the noise vector, but it is easier
        // to just fill it with 0's.
        _noise.resize(nrows,0);
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
        for(int i=0;i<nrows;++i) {
            if (_flags[i]) {
                dbg<<std::dec<<std::noshowbase;
                dbg<<"Object "<<i;
                dbg<<std::hex<<std::showbase;
                dbg<<" has sextractor flag "<<_flags[i];
                if (_flags[i] & ignoreFlags) {
                    dbg<<"    Marking with flag INPUT_FLAG\n";
                    _flags[i] = INPUT_FLAG;
                } else {
                    dbg<<"\n";
                }
            } else {
                _flags[i] = 0;
            }
        }
        dbg<<std::dec<<std::noshowbase;
    }

    const int nPos = _pos.size();
    // If we have an image, then also flag any objects whose input position
    // is outside the bounds of the image.
    if (_im) {
        for(int i=0;i<nPos;++i) {
            Bounds imb = _im->getBounds();
            if (!imb.includes(_pos[i])) {
                dbg<<"Object "<<i<<" has invalid position "<<_pos[i]<<std::endl;
                _flags[i] = INPUT_FLAG;
            }
        }
    }

    // Update the bounds with valid input positions.
    for(int i=0;i<nPos;++i) if (!_flags[i]) _bounds += _pos[i];

    // Calculate the skybounds using only the good objects:
    if (_skypos.size() == _id.size()) {
        Bounds skybounds2; // In degrees, rather than arcsec
        for(int i=0;i<nrows;++i) if (!_flags[i]) {
            _skybounds += _skypos[i];
            Position temp = _skypos[i];
            temp /= 3600.;
            skybounds2 += temp;
        }
        dbg<<"skybounds = "<<_skybounds<<std::endl;
        dbg<<"in degrees: "<<skybounds2<<std::endl;
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
    int hdu = GetHdu(_params,"cat",file,1);
    readFits(file, hdu);
}

void InputCatalog::readFits(std::string file, int hdu)
{
    dbg<< "Reading cat from FITS file: " << file << std::endl;

    dbg<<"Opening InputCatalog file "<<file<<" at hdu "<<hdu<<std::endl;
    CCfits::FITS fits(file, CCfits::Read);
    if (hdu > 1) fits.read(hdu-1);

    CCfits::ExtHDU& table=fits.extension(hdu-1);

    long nrows=table.rows();

    dbg<<"  nrows = "<<nrows<<std::endl;
    if (nrows <= 0) {
        throw ReadException(
            "InputCatalog found to have 0 rows.  Must have > 0 rows.");
    }

    // Read each column in turn:
    //
    long start=1;
    long end=nrows;

    dbg<<"Reading columns"<<std::endl;


    if (_params.keyExists("cat_id_col")) {
        std::string id_col = _params["cat_id_col"];
        dbg<<"  "<<id_col<<std::endl;
        table.column(id_col).read(_id, start, end);
    } else {
        _id.resize(nrows);
        for (int i=0;i<nrows;++i) _id[i] = i+1;
    }
    Assert(int(_id.size()) == nrows);

    // Position (on chip)
    std::vector<double> posX;
    std::vector<double> posY;
    _pos.resize(nrows);
    std::string x_col=_params.get("cat_x_col");
    std::string y_col=_params.get("cat_y_col");
    dbg<<"  "<<x_col<<std::endl;
    table.column(x_col).read(posX, start, end);
    dbg<<"  "<<y_col<<std::endl;
    table.column(y_col).read(posY, start, end);

    for (long i=0; i< nrows; ++i) {
        _pos[i] = Position(posX[i], posY[i]);
    }

    // Local sky
    if (_params.keyExists("cat_sky_col")) {
        std::string sky_col=_params["cat_sky_col"];
        dbg<<"  "<<sky_col<<std::endl;
        table.column(sky_col).read(_sky, start, end);
    }

    // Magnitude
    if (_params.keyExists("cat_mag_col")) {
        std::string mag_col=_params["cat_mag_col"];
        dbg<<"  "<<mag_col<<std::endl;
        table.column(mag_col).read(_mag, start, end);
    }

    // Magnitude Error
    if (_params.keyExists("cat_mag_err_col")) {
        std::string mag_err_col=_params["cat_mag_err_col"];
        dbg<<"  "<<mag_err_col<<std::endl;
        table.column(mag_err_col).read(_mag_err, start, end);
    }

    // Star Galaxy separation
    if (_params.keyExists("cat_sg_col")) {
        std::string sg_col=_params["cat_sg_col"];
        dbg<<"  "<<sg_col<<std::endl;
        table.column(sg_col).read(_sg, start, end);
    } 

    // Size
    if (_params.keyExists("cat_size_col")) {
        std::string size_col=_params["cat_size_col"];
        dbg<<"  "<<size_col<<std::endl;
        table.column(size_col).read(_obj_size, start, end);
        if (_params.keyExists("cat_size2_col")) {
            std::string size2_col=_params["cat_size2_col"];
            dbg<<"  "<<size2_col<<std::endl;
            std::vector<double> obj_size2(_obj_size.size());
            table.column(size2_col).read(obj_size2, start, end);
            Assert(_obj_size.size() == _id.size());
            Assert(obj_size2.size() == _id.size());
            for(int i=0;i<nrows;++i) _obj_size[i] += obj_size2[i];
        }
    }

    // Error flags
    if (_params.keyExists("cat_flag_col")) {
        std::string flag_col=_params["cat_flag_col"];
        dbg<<"  "<<flag_col<<std::endl;
        table.column(flag_col).read(_flags, start, end);
    }


    // RA
    if (_params.keyExists("cat_ra_col")) {
        if (!_params.keyExists("cat_dec_col"))
            throw ParameterException("cat_ra_col, but no cat_dec_col");

        std::string ra_col=_params["cat_ra_col"];
        std::string decl_col=_params["cat_dec_col"];

        std::vector<double> ra;
        std::vector<double> decl;
        dbg<<"  "<<ra_col<<std::endl;
        table.column(ra_col).read(ra, start, end);
        dbg<<"  "<<decl_col<<std::endl;
        table.column(decl_col).read(decl, start, end);

        _skypos.resize(nrows);
        for(long i=0;i<nrows;++i) {
            _skypos[i] = Position(ra[i],decl[i]);

            // The convention for Position is to use arcsec for everything.
            // ra and dec come in as degrees.  So wee need to convert to arcsec.
            _skypos[i] *= 3600.;  // deg -> arcsec
        }
    } else if (_params.keyExists("cat_dec_col")) {
        throw ParameterException("cat_dec_col, but no cat_ra_col");
    }

    // Noise
    if (_params.keyExists("cat_noise_col")) {
        std::string noise_col=_params["cat_noise_col"];
        dbg<<"  "<<noise_col<<std::endl;
        table.column(noise_col).read(_noise, start, end);
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
    int x_col = _params.read<int>("cat_x_col");
    int y_col = _params.read<int>("cat_y_col");

    // Set column numbers for optional columns
    int id_col = _params.read("cat_id_col",0);

    int mag_col = _params.read("cat_mag_col",0);
    int mag_err_col = _params.read("cat_mag_err_col",0);

    int sg_col = _params.read("cat_sg_col",0);

    int size_col = _params.read("cat_size_col",0);
    int size2_col = _params.read("cat_size2_col",0);

    int flag_col = _params.read("cat_flag_col",0);

    int sky_col = _params.read("cat_sky_col",0);

    int ra_col = _params.read("cat_ra_col",0);
    int decl_col = _params.read("cat_dec_col",0);

    int noise_col = _params.read("cat_noise_col",0);

    // Set up allowed comment markers
    std::vector<std::string> commentMarker = 
        _params.read("cat_comment_marker",std::vector<std::string>(1,"#"));

    // Keep running id value when id_col = 0
    int id_val = 0;

    // Read each line from catalog file
    while (getline(catIn,line)) {

        // Skip if this is a comment.
        bool skip = false;
        const int nCommentMarkers = commentMarker.size();
        for(int k=0;k<nCommentMarkers;++k) {
            if (std::string(line,0,commentMarker[k].size()) == 
                commentMarker[k]) {
                skip = true;
            }
        }
        if (skip) continue;

        // Convert line into vector of tokens
        std::vector<ConvertibleString> tokens;
        getTokens(delim,line,tokens);

        // ID
        if (id_col) {
            Assert(id_col <= int(tokens.size()));
            id_val = tokens[id_col-1];
        } else {
            // if not reading id, then just increment to get sequential values
            ++id_val;
        }
        _id.push_back(id_val);

        // Position
        Assert(x_col <= int(tokens.size()));
        Assert(y_col <= int(tokens.size()));
        double x = tokens[x_col-1];
        double y = tokens[y_col-1];
        _pos.push_back(Position(x,y));

        // Sky
        if (sky_col) {
            // Note: if sky not read in, then sky.size() is still 0
            // This is indicator to update with global given or median 
            // value later.
            double sky_val = 0.;
            Assert(sky_col <= int(tokens.size()));
            sky_val = tokens[sky_col-1];
            _sky.push_back(sky_val);
        } 

        // Magnitude
        if (mag_col) {
            double mag_val = 0.;
            Assert(mag_col <= int(tokens.size()));
            mag_val = tokens[mag_col-1];
            _mag.push_back(mag_val);
        } 

        // Magnitude error
        if (mag_err_col) {
            double mag_err_val = 0.;
            Assert(mag_err_col <= int(tokens.size()));
            mag_err_val = tokens[mag_err_col-1];
            _mag_err.push_back(mag_err_val);
        } 

        // Star-galaxy
        if (sg_col) {
            double sg_val = 0.;
            Assert(sg_col <= int(tokens.size()));
            sg_val = tokens[sg_col-1];
            _sg.push_back(sg_val);
        } 

        // Size
        if (size_col) {
            double size_val = 0.;
            Assert(size_col <= int(tokens.size()));
            size_val = tokens[size_col-1];
            if (size2_col) {
                Assert(size2_col <= int(tokens.size()));
                size_val += double(tokens[size2_col-1]);
            } 
            _obj_size.push_back(size_val);
        } 

        // Flags
        if (flag_col) {
            long flag_val = 0;
            Assert(flag_col <= int(tokens.size()));
            flag_val = tokens[flag_col-1];
            _flags.push_back(flag_val);
        } 

        // RA
        if (ra_col) {
            if (!decl_col) 
                throw ParameterException("cat_ra_col, but no cat_dec_col");
            double ra_val = 0.;
            double decl_val = 0.;
            Assert(ra_col <= int(tokens.size()));
            Assert(decl_col <= int(tokens.size()));
            ra_val = tokens[ra_col-1];
            decl_val = tokens[decl_col-1];
            Position skypos_val(ra_val,decl_val);
            skypos_val *= 3600.; // deg -> arcsec
            _skypos.push_back(skypos_val);
        } else if (decl_col) {
            throw ParameterException("cat_dec_col, but no cat_ra_col");
        }

        // Noise
        if (noise_col) {
            double noise_val = 0.;
            Assert(noise_col <= int(tokens.size()));
            noise_val = tokens[noise_col-1];
            _noise.push_back(noise_val);
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

    if (int(_mag_err.size()) > i) {
        std::cout<<"  InputCatalog::printall mag_err["<<i<<"]: "<<_mag_err[i]<<"\n";
    } else {
        std::cout<<"  InputCatalog::printall index "<<i<<" is larger than mag_err.size\n";
    }

    if (int(_sg.size()) > i) {
        std::cout<<"  InputCatalog::printall sg["<<i<<"]: "<<_sg[i]<<"\n";
    } else {
        std::cout<<"  InputCatalog::printall index "<<i<<" is larger than sg.size\n";
    }

    if (int(_obj_size.size()) > i) {
        std::cout<<"  InputCatalog::printall obj_size["<<i<<"]: "<<_obj_size[i]<<"\n";
    } else {
        std::cout<<"  InputCatalog::printall index "<<i<<" is larger than obj_size.size\n";
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
