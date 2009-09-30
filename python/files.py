from sys import stdout,stderr
import os
import re
import esutil
from esutil import xmltools
from esutil.ostools import path_join, getenv_check
import pprint
import deswl


class Runconfig(object):
    def __init__(self, run=None):
        
        self.run=run
        self.fileclass='wlbnl'

        self.run_types={}
        self.run_types['wlse'] = {'name':'wlse',
                                  'fileclass': self.fileclass, 
                                  'filetype':'se_shapelet',
                                  'runvarname':'serun'}
        self.run_types['wlme'] = {'name':'wlme',
                                  'fileclass': self.fileclass, 
                                  'filetype':'me_shapelet',
                                  'runvarname': 'merun'}


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

    def get_basedir(self):
        basedir=getenv_check('DESFILES_DIR')
        return path_join(basedir, 'runconfig')

    def getpath(self, run=None):
        if run is None:
            if self.run is None:
                raise ValueError("Either send run= keyword or load a runconfig")
            else:
                run=self.run
        rdir = self.get_basedir()
        return path_join(rdir, run+'-config.xml')
    

    def generate_new_runconfig_name(self, run_type, test=False):
        """
        generates a new name like run_type000002, checking against names in
        the runconfig directory
        """

        if run_type not in self.run_types:
            raise ValueError("Unknown run type: '%s'.  Must be one "
                             "of: (%s)" % (run_type, ', '.join(self.run_types)))


        i=0
        run_name=self._run_name_from_type_number(run_type, i, test=test)
        fullpath=self.getpath(run_name)
        while os.path.exists(fullpath):
            i+=1
            run_name=self._run_name_from_type_number(run_type, i, test=test)
            fullpath=self.getpath(run_name)

        return run_name

    def _run_name_from_type_number(self, run_type, number, test=False):
        if test:
            name='%stest%06i' % (run_type, number)
        else:
            name='%s%06i' % (run_type, number)
        return name


    def generate_new_runconfig(self, run_type, dataset, localid, 
                               test=False, **extra):

        # wlme runs depend on wlse runs
        if run_type == 'wlme':
            err="You must send serun=something for wlme"
            if 'serun' not in extra:
                raise RuntimeError(mess)


        run_name=self.generate_new_runconfig_name(run_type, test=test)
        config={}

        fileclass = self.run_types[run_type]['fileclass']
        filetype = self.run_types[run_type]['filetype']

        # software versions
        pyvers=deswl.get_python_version()
        esutilvers=esutil.version()
        wlvers=deswl.version()
        tmvvers=deswl.get_tmv_version()

        runconfig={'run':run_name, 
                   'run_type':run_type,
                   'fileclass': fileclass,
                   'filetype':filetype,
                   'pyvers':pyvers, 
                   'esutilvers': esutilvers,
                   'wlvers':wlvers,
                   'tmvvers':tmvvers,
                   'dataset':dataset,
                   'localid':localid}
        for e in extra:
            if e not in runconfig:
                runconfig[e] = extra[e]

        fullpath=self.getpath(run_name)
        stdout.write('Writing to file: %s\n' % fullpath)
        xmltools.dict2xml(runconfig, fullpath, roottag='runconfig')

        stdout.write("\n * Don't forget to check it into SVN!!!\n\n")
        return runconfig

    def verify(self):
        """

        verify current setup against the run config, make sure all versions and
        such are the same.  This does not check the data set or localid or
        serun for merun runconfigs.

        """

        wlvers=deswl.version()
        esutilvers=esutil.version()
        tmvvers=deswl.get_tmv_version()
        pyvers=deswl.get_python_version()

        if self.config['esutilvers'] != esutilvers:
            raise ValueError("current esutil '%s' does not match "
                    "runconfig '%s'" % (esutilvers,self.config['esutilvers']))

        if self.config['wlvers'] != wlvers:
            raise ValueError("current wlvers '%s' does not match "
                             "runconfig '%s'" % (wlvers,self.config['wlvers']))

        if self.config['tmvvers'] != tmvvers:
            raise ValueError("current tmvvers '%s' does not match "
                        "runconfig '%s'" % (tmvvers,self.config['tmvvers']))

        if self.config['pyvers'] != pyvers:
            raise ValueError("current pyvers '%s' does not match "
                        "runconfig '%s'" % (pyvers,self.config['pyvers']))



    def load(self, run):
        self.config = self.read(run)
        self.run=run

    def read(self, run):
        name=self.getpath(run)
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

def red_image_path(run, exposurename, ccdnum, fz=True, rootdir=None, check=False):
    fileclass='red'
    filetype='red'
    expdir = exposure_dir(fileclass, run, filetype, exposurename, 
                          rootdir=rootdir)
    imagename = '%s_%02i.fits' % (exposurename, int(ccdnum))
    basic_path = os.path.join(expdir, imagename)

    if check:
        # check both with and without .fz
        path=basic_path
        if not os.path.exists(path):
            path += '.fz'
            if not os.path.exists(path):
                stdout.write("SE image not found: %s(.fz)\n" % basic_path)
                return None
    else:
        if fz:
            path=basic_path + '.fz'
        
    return path

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

def coadd_image_path(coaddrun, tilename, band, fz=True, rootdir=None, check=False):
    """
    More frustration:  the path element for "coaddrun"
    contains the tilename but in principle it can be arbitrary, so 
    we can't simply pass in pieces:  coaddrun is treated independently

    If check=True, fz will be also be tried if not already set to True
    """
    tiledir=tile_dir(coaddrun, rootdir=rootdir)
    fname=tilename+'_'+band+'.fits'
    basic_path=path_join(tiledir, fname)
    if check:
        # check both with and without .fz
        path=basic_path
        if not os.path.exists(path):
            path += '.fz'
            if not os.path.exists(path):
                stdout.write("coadd image not found: %s(.fz)\n" % basic_path)
                return None
    else:
        if fz:
            path=basic_path + '.fz'
    return path

def coadd_cat_path(catalogrun, tilename, band, rootdir=None, check=False):
    """
    More frustration:  the path element for "catalogrun"
    contains the tilename but in principle it can be arbitrary, so 
    we can't simply pass in pieces:  catalogrun is treated independently

    Note often catalogrun is the same as the coaddrun
    """
    tiledir=tile_dir(catalogrun, rootdir=rootdir)
    fname=tilename+'_'+band+'_cat.fits'
    path=path_join(tiledir, fname)
    if check:
        if not os.path.exists(path):
            stdout.write("coadd catalog not found: %s(.fz)\n" % path)
            return None
    return path


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

def wlse_path(run, redrun, exposurename, ccdnum, wltype, 
              mohrify=True, rootdir=None, ext='.fits'):
    """
    name=wlse_path(run, redrun, exposurename, ccdnum, wltype, 
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
    name = '%s_%s_%02i_%s%s' % (run,name_front, int(ccdnum), wltype_use, ext)
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


def wlme_path(merun, tilename, band, extra=None, 
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


def collated_coaddfiles_dir(dataset, localid=None):
    desfiles_dir=getenv_check('DESFILES_DIR')
    desfiles_dir = path_join(desfiles_dir,dataset)
    if localid is not None:
        desfiles_dir=path_join(desfiles_dir, localid)
    return desfiles_dir


# need to rationalize these names
def collated_coaddfiles_name(dataset, band, 
                             withsrc=True, 
                             serun=None,
                             localid=None, 
                             newest=True):
    name=[dataset,'images','catalogs']
    if withsrc:
        name.append('withsrc')

    if localid is not None:
        name.append(localid)
        if newest:
            name.append('newest')

    if serun is not None:
        name.append(serun)
    name.append(band)

    return '-'.join(name)+'.xml'

def collated_coaddfiles_path(dataset, band, 
                             withsrc=True, 
                             serun=None,
                             localid=None, 
                             newest=True):
    tdir=collated_coaddfiles_dir(dataset, localid=localid)
    name = collated_coaddfiles_name(dataset, band, 
                                    withsrc=withsrc, 
                                    serun=serun,
                                    localid=localid, 
                                    newest=newest)
    return path_join(tdir, name)

def collated_coaddfiles_read(dataset, band, 
                             withsrc=True, 
                             serun=None,
                             localid=None, 
                             newest=True, 
                             getname=False):
    f=collated_coaddfiles_path(dataset, band, 
                               withsrc=withsrc, 
                               serun=serun,
                               localid=localid, 
                               newest=newest)
    stdout.write('Reading tileinfo file: %s\n' % f)
    tileinfo=xmltools.xml2dict(f,noroot=True)
    if getname:
        return tileinfo, f
    else:
        return tileinfo



def wlme_srclist_dir(dataset, localid):
    sdir=collated_coaddfiles_dir(dataset, localid=localid)
    sdir=path_join(sdir, 'multishear-srclists')
    return sdir
def wlme_srclist_name(tilename, band, serun):
    srclist=[tilename,band,serun,'srclist']
    srclist='-'.join(srclist)+'.dat'
    return srclist
    

def wlme_srclist_path(dataset, localid, tilename, band, serun):
    dir=wlme_srclist_dir(dataset, localid)
    srclist=wlme_srclist_name(tilename,band,serun)
    return path_join(dir, srclist)




def wlme_pbs_dir(dataset, localid, tilename, band, merun):
    desfiles_dir=getenv_check('DESFILES_DIR')
    pbsdir=path_join(desfiles_dir,dataset,localid,'multishear-pbs', merun)
    return pbsdir

def wlme_pbs_name(tilename, band, merun):
    pbsfile=[tilename,band,merun]
    pbsfile='-'.join(pbsfile)+'.pbs'
    return pbsfile
    
def wlme_pbs_path(dataset, localid, tilename, band, merun):
    pbsfile=wlme_pbsfile_name(tilename,band, merun)
    pbsdir=wlme_pbs_dir(tilename, band, merun)
    pbsfile=path_join(pbsdir, pbsfile)
    return pbsfile




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




