import os
import re
from esutil import xmltools
from esutil.ostools import path_join, getenv_check
import pprint


class Runconfig(object):
    def __init__(self, run=None):
        
        self.fileclass='wlbnl'

        self.run_types={}
        self.run_types['wlse'] = {'name':'wlse',
                                  'fileclass': self.fileclass, 
                                  'filetype':'se_shapelet'}
        self.run_types['wlme'] = {'name':'wlme',
                                  'fileclass': self.fileclass, 
                                  'filetype':'me_shapelet'}


        self.fileclasses = {}
        self.fileclasses['wlse'] = self.fileclass


        self.se_executables = ['findstars','measurepsf','measureshear']
        self.se_filetypes = ['stars','psf','fitpsf','shear',
                             'findstars_log', 
                             'measurepsf_log', 
                             'measureshear_log']
        self.me_executables = ['multishear']
        self.me_filetypes = ['multishear','multishear_log']

        self.config={}
        if run is not None:
            self.load(run)

    def getbasedir(self):
        DESDATA=getenv_check('DESDATA')
        return os.path.join(DESDATA, 'runconfig')

    def getdir(self, run):
        basedir=self.getbasedir()
        return os.path.join(basedir, run)

    def getname(self, run):
        rdir = self.getdir(run)
        return os.path.join(rdir, run+'-config.xml')

    def load(self, run):
        self.config = self.read(run)

    def read(self, run):
        name=self.getname(run)
        if not os.path.exists(name):
            mess="runconfig '%s' not found at %s\n" % (run, name)
            raise RuntimeError(mess)
        runconfig = xmltools.xml2dict(name, noroot=True)
        return runconfig

    def asdict(self):
        return self.config

    def __getitem__(self, name):
        if name in self.config:
            return self.config[name]
        else:
            raise KeyError("%s" % name)

    def __iter__(self):
        return self.config.__iter__()

    def __repr__(self):
        return pprint.pformat(self.config)


def fileclass_dir(fileclass, rootdir=None):
    if rootdir is not None:
        DESDATA=rootdir
    else:
        DESDATA=getenv_check('DESDATA')
    return os.path.join(DESDATA, fileclass)
    
def run_dir(fileclass, run, rootdir=None):
    fcdir = fileclass_dir(fileclass, rootdir=rootdir)
    return os.path.join(fcdir, run)

def filetype_dir(fileclass, run, filetype, rootdir=None):
    rundir=run_dir(fileclass, run, rootdir=rootdir)
    return os.path.join(rundir, filetype)

def exposure_dir(fileclass, run, filetype, exposurename, rootdir=None):
    ftdir=filetype_dir(fileclass, run, filetype, rootdir=rootdir)
    return os.path.join(ftdir, exposurename)

def red_image_name(run, exposurename, ccdnum, fz=True, rootdir=None):
    fileclass='red'
    filetype='red'
    expdir = exposure_dir(fileclass, run, filetype, exposurename, 
                         rootdir=rootdir)
    imagename = '%s_%02i.fits' % (exposurename, ccdnum)
    if fz:
        imagename=imagename+'.fz'
    imagename = os.path.join(expdir, imagename)
    return imagename

def extract_image_exposure_names(flist):
    allinfo={}
    for imfile in flist:
        info=get_info_from_path(imfile, 'red')
        allinfo[info['exposurename']] = info['exposurename']

    return allinfo.keys()


def tile_dir(coaddrun, rootdir=None):
    """
    Coadds are different than other output files:  The fileclass dir is 
    where the files are!
    """
    fileclass='coadd'
    filetype='coadd'
    tiledir=filetype_dir(fileclass, coaddrun, filetype, rootdir=rootdir)
    return tiledir

def coadd_image_name(coaddrun, tilename, band, fz=True):
    """
    More frustration:  the path element for "coaddrun"
    contains the tilename but in principle it can be arbitrary, so 
    we can't simply pass in pieces:  coaddrun is treated independently
    """
    tiledir=tile_dir(coaddrun)
    fname=tilename+'_'+band+'.fits'
    if fz:
        fname+='.fz'
    return os.path.join(tiledir, fname)

def coadd_cat_name(catalogrun, tilename, band):
    """
    More frustration:  the path element for "catalogrun"
    contains the tilename but in principle it can be arbitrary, so 
    we can't simply pass in pieces:  catalogrun is treated independently

    Note often catalogrun is the same as the coaddrun
    """
    tiledir=tile_dir(catalogrun)
    fname=tilename+'_'+band+'_cat.fits'
    return fname


_MOHR_MAP={'stars':'shpltall', 
           'psf':'shpltpsf', 
           'fitpsf':'psfmodel', 
           'shear':'shear'}
def mohrify_name(name):
    if name in _MOHR_MAP:
        return _MOHR_MAP[name]
    else:
        return name



# names and directories for  weak lensing processing
# se means single epoch
# me means multi-epoch

def wlse_dir(run, exposurename, rootdir=None):
    rc=Runconfig()

    fileclass=rc.run_types['wlse']['fileclass']
    filetype=rc.run_types['wlse']['filetype']

    ftdir=filetype_dir(fileclass, run, filetype, rootdir=rootdir)
    wldir=os.path.join(ftdir, exposurename)
 
    return wldir

def wlme_dir(merun, tilename, rootdir=None):
    rc=Runconfig()

    fileclass=rc.run_types['wlme']['fileclass']
    filetype=rc.run_types['wlme']['filetype']

    ftdir=filetype_dir(fileclass, merun, filetype, rootdir=rootdir)
    wldir=os.path.join(ftdir, tilename)
 
    return wldir

def wlse_name(run, redrun, exposurename, ccdnum, wltype, 
              mohrify=True, rootdir=None, ext='.fits'):
    """
    name=wlse_name(run, redrun, exposurename, ccdnum, wltype, 
                  mohrify=True, rootdir=None, ext='.fits')
    Return the SE output file name for the given inputs

    Because this is designed for single epoch image processing, the
    resulting directory structure is different than uiuc one running
    from the tiling:

    fileclass/run/filetype/redrun_exposurename/
                                  redrun_exposurename_ccd_wltype.fits
    """
    if mohrify:
        wltype_use=mohrify_name(wltype)
    else:
        wltype_use=wltype

    wldir = wlse_dir(run, exposurename, rootdir=rootdir)
    name_front = redrun+'_'+os.path.basename(wldir)
    name = '%s_%s_%02i_%s.%s' % (run,name_front, int(ccdnum), wltype_use, ext)
    fullpath = os.path.join(wldir, name)
    return fullpath

def wlme_basename(tilename, band, extra=None, ext='.fits'):
    """
    basename for multi-epoch weak lensing output files
    """

    name=[tilename,band,'multishear']
    if extra is not None:
        name.append(extra)

    name='_'.join(name)+ext
    return name


def wlme_fullpath(merun, tilename, band, extra=None, 
                  dir=None, rootdir=None, ext='.fits'):
    """
    Return a multi-epoch shear output file name

    currently just for multishear and logs, etc 
    """

    name = wlme_basename(tilename, band, extra=extra,ext=ext)

    if dir is None:
        dir=wlme_dir(merun, tilename, rootdir=rootdir)

    name = os.path.join(dir, name)
    return name

def wlme_srclist_name(tilename, band):
    srclist=[tilename,band,'srclist']
    srclist='-'.join(srclist)+'.dat'
    return srclist
    

def wlme_srclist_fullpath(dataset, host, tilename, band):
    srclist=wlme_srclist_name(tilename,band)
    desfiles_dir=getenv_check('DESFILES_DIR')
    srclist=path_join(desfiles_dir,dataset,host,'multishear-srclists',srclist)
    return srclist


def wlme_pbs_name(merun,tilename, band):
    pbsfile=[tilename,band,merun]
    pbsfile='-'.join(pbsfile)+'.pbs'
    return pbsfile
    
def wlme_pbs_fullpath(dataset, host, merun, tilename, band):
    pbsfile=wlme_pbsfile_name(tilename,band)
    desfiles_dir=getenv_check('DESFILES_DIR')
    pbsfile=path_join(desfiles_dir,dataset,host,merun,'multishear-pbs', pbsfile)
    return pbsfile




# run configurations
def runconfig_basedir():
    DESDATA=getenv_check('DESDATA')
    return os.path.join(DESDATA, 'runconfig')

def runconfig_dir(run):
    basedir=runconfig_basedir()
    return os.path.join(basedir, run)

def runconfig_name(run):
    rdir=runconfig_dir(run)
    return os.path.join(rdir, run+'-config.xml')


def get_info_from_path(filepath, fileclass):
    """
    This is a more extensive version get_info_from_image_path.  Also works for
    the cat files.  If fileclass='coadd' then treated differently
    """

    odict = {}
    # remove the extension
    dotloc=filepath.find('.')
    
    odict['ext'] = filepath[dotloc:]
    root=filepath[0:dotloc]



    DESDATA=getenv_check('DESDATA')
    if DESDATA[-1] == os.sep:
        DESDATA=DESDATA[0:len(DESDATA)-1]

    endf = filepath.replace(DESDATA, '')


    odict['filename'] = filepath
    odict['root'] = root
    odict['DESDATA'] = DESDATA


    fsplit = endf.split(os.sep)


    if fileclass == 'red':
        nparts = 6
        if len(fsplit) < nparts:
            raise ValueError,\
                'Found too few file elements, '+\
                '%s instead of %s: "%s"\n' % (len(fsplit), nparts, filepath)

        if fsplit[1] != 'red' or fsplit[3] != 'red':
            mess="In file path %s expected 'red' in positions "+\
                "1 and 3 after DESDATA=%s" % (filepath,DESDATA)
            raise ValueError,mess

        odict['redrun'] = fsplit[2]
        odict['exposurename'] = fsplit[4]
        odict['basename'] = fsplit[5]

        base=os.path.basename(odict['root'])


        # last 3 bytes are _[ccd number]
        odict['ccd']=int(base[-2:])
        base=base[0:len(base)-3]

        # next -[repeatnum].  not fixed length so we must go by the dash
        dashloc=base.rfind('-')
        odict['repeatnum'] = base[dashloc+1:]
        base=base[0:dashloc]

        # now is -[band]
        dashloc=base.rfind('-')
        band=base[-1]
        odict['band'] = band
        base=base[0:dashloc]
        
        # now we have something like decam--24--11 or decam--24--9, so
        # just remove the decam-
        # might have to alter this some time if they change the names
        if base[0:6] != 'decam-':
            raise ValueError,'Expected "decam-" at beginning'
        odict['pointing'] = base.replace('decam-','')

        
    elif fileclass == 'coadd':
        nparts = 5
        if len(fsplit) < nparts:
            raise ValueError,\
                'Found too few file elements, '+\
                '%s instead of %s: "%s"\n' % (len(fsplit), nparts, filepath)

        if fsplit[1] != 'coadd' or fsplit[3] != 'coadd':
            mess="In file path %s expected 'red' in positions "+\
                "1 and 3 after DESDATA=%s" % (filepath,DESDATA)
            raise ValueError,mess


        odict['coaddrun'] = fsplit[2]
        odict['basename'] = fsplit[4]

        tmp = os.path.basename(root)
        if tmp[-3:] == 'cat':
            odict['coaddtype'] = 'cat'
            tmp = tmp[0:-4]
        else:
            odict['coaddtype'] = 'image'

        odict['band'] = tmp[-1]

        last_underscore = tmp.rfind('_')
        if last_underscore == -1:
            raise ValueError,\
                'Expected t find an underscore in basename of %s' % filepath
        odict['tilename'] = tmp[0:last_underscore]

    
    return odict


# utilities
def replacedir(oldfile, newdir):
    if oldfile is not None:
        fbase=os.path.basename(oldfile)
        newf = os.path.join(newdir, fbase)
        return newf
    else:
        return None


_FITS_EXTSTRIP=re.compile( '\.fits(\.fz)?$' )
def remove_fits_extension(fname):
    """
    See if the pattern exists and replace it with '' if it does
    """
    return _FITS_EXTSTRIP.sub('', fname)




