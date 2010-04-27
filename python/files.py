from sys import stdout,stderr
import platform
import os
import re
import esutil
from esutil import json_util
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

        self.se_filetypes = ['stars', 'psf',  'fitpsf','shear','qa',  'stat']
        self.se_fext       = {'stars':  '.fits', 
                              'psf':    '.fits',
                              'fitpsf': '.fits', 
                              'shear':  '.fits',
                              'qa':     '.dat',
                              'stat':   '.json',
                              'checkpsf':'.rec'}
        self.se_collated_type = ['badlist','goodlist','collated']
        self.se_collated_filetypes = {'badlist':'json',
                                      'goodlist':'json',
                                      'gal':'fits'}

        self.me_filetypes = ['multishear','qa',  'stat']
        self.me_fext       = {'multishear':'.fits',
                              'qa':'.dat',
                              'stat':'.json'}


        self.me_executables = ['multishear']

        self.config={}
        if run is not None:
            self.load(run)

    def get_basedir(self):
        basedir=getenv_check('DESFILES_DIR')
        return path_join(basedir, 'runconfig')

    def getpath(self, run=None):
        if run is None:
            if self.run is None:
                raise ValueError("Either send run= keyword or "
                                 "load a runconfig")
            else:
                run=self.run
        rdir = self.get_basedir()
        return path_join(rdir, run+'-config.json')
    

    def generate_new_runconfig_name(self, run_type, test=False):
        """
        generates a new name like run_type000002, checking against names in
        the runconfig directory
        """

        if run_type not in self.run_types:
            mess="Unknown run type: '%s'.  Must be one "+\
                 "of: (%s)" % (run_type, ', '.join(self.run_types))
            raise ValueError(mess)

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
            name='%stest%04i' % (run_type, number)
        else:
            name='%s%04i' % (run_type, number)
        return name


    def generate_new_runconfig(self, 
                               run_type, 
                               dataset, 
                               wl_config,
                               localid=None, 
                               test=False, 
                               dryrun=False, 
                               pyvers=None,
                               esutilvers=None,
                               wlvers=None,
                               tmvvers=None,
                               comment=None,
                               **extra):

        # wlme runs depend on wlse runs
        if run_type == 'wlme':
            if 'serun' not in extra:
                raise RuntimeError("You must send serun=something for wlme")
            if localid is None:
                raise RuntimeError("You must send localid= for wlme runs")


        run_name=self.generate_new_runconfig_name(run_type, test=test)

        fileclass = self.run_types[run_type]['fileclass']
        filetype = self.run_types[run_type]['filetype']

        # software versions.  Default to whatever is in our environment
        if pyvers is None:
            pyvers=deswl.get_python_version()
        if esutilvers is None:
            esutilvers=esutil.version()
        if wlvers is None:
            wlvers=deswl.version()
        if tmvvers is None:
            tmvvers=deswl.get_tmv_version()

        runconfig={'run':run_name, 
                   'run_type':run_type,
                   'fileclass': fileclass,
                   'pyvers':pyvers, 
                   'esutilvers': esutilvers,
                   'wlvers':wlvers,
                   'tmvvers':tmvvers,
                   'dataset':dataset,
                   'wl_config':wl_config}

        if comment is not None:
            runconfig['comment'] = comment

        if localid is not None:
            runconfig['localid'] = localid

        for e in extra:
            if e not in runconfig:
                runconfig[e] = extra[e]

        fullpath=self.getpath(run_name)
        stdout.write('Writing to file: %s\n' % fullpath)
        if not dryrun:
            json_util.write(runconfig, fullpath)
            stdout.write("\n * Don't forget to check it into SVN!!!\n\n")
        else:
            stdout.write(" .... dry run, skipping file write\n")

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
        runconfig = json_util.read(name)
        return runconfig

    def asdict(self):
        from copy import copy
        conf=copy(self.config)
        return conf

    def __getitem__(self, name):
        if name in self.config:
            return self.config[name]
        else:
            raise KeyError("%s" % name)

    def __iter__(self):
        return self.config.__iter__()

    def __repr__(self):
        return pprint.pformat(self.config)


def convert_to_degrees(data):
    try:
        if 'ra' in data.dtype.names:
            data['ra'] /= 3600

        if 'dec' in data.dtype.names:
            data['dec'] /= 3600
    except:
        # probably names is None
        pass





def des_rootdir():
    return getenv_check('DESDATA')


def fileclass_dir(fileclass, rootdir=None):
    if rootdir is None:
        rootdir=des_rootdir()
    return os.path.join(rootdir, fileclass)


def run_dir(fileclass, run, rootdir=None):
    if rootdir is None:
        rootdir=des_rootdir()
    dir=path_join(rootdir, fileclass, run)
    return dir
 

def filetype_dir(fileclass, run, filetype, rootdir=None):
    rundir=run_dir(fileclass, run, rootdir=rootdir)
    return os.path.join(rundir, filetype)

def exposure_dir(fileclass, run, filetype, exposurename, rootdir=None):
    ftdir=filetype_dir(fileclass, run, filetype, rootdir=rootdir)
    return os.path.join(ftdir, exposurename)

def red_image_path(run, exposurename, ccd, fz=True, rootdir=None, check=False):
    fileclass='red'
    filetype='red'
    expdir = exposure_dir(fileclass, run, filetype, exposurename, 
                          rootdir=rootdir)
    imagename = '%s_%02i.fits' % (exposurename, int(ccd))
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

    # list() of for py3k
    return list( allinfo.keys() )


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
    rundir=run_dir(fileclass, run, rootdir=rootdir)
    dir=path_join(rundir, exposurename)
    return dir

def wlse_basename(exposurename, ccd, ftype,
                  serun=None,
                  mohrify=False, fext=None):

    if fext is None:
        rc=Runconfig()
        if ftype in rc.se_fext:
            fext=rc.se_fext[ftype]
        else:
            fext='.fits'

    el=[exposurename, '%02i' % int(ccd), ftype]

    if serun is not None:
        el=[serun]+el

    basename='_'.join(el) + fext
    return basename



def wlse_path(exposurename, ccd, ftype, 
              serun=None,
              mohrify=False, 
              dir=None, 
              rootdir=None, 
              fext=None):
    """
    name=wlse_path(exposurename, ccd, ftype, serun=None,
                   mohrify=False, rootdir=None, fext=None)
    Return the SE output file name for the given inputs
    """
    if mohrify:
        ftype_use=mohrify_name(ftype)
    else:
        ftype_use=ftype

    if dir is None:
        if serun is not None:
            dir = wlse_dir(serun, exposurename, rootdir=rootdir)
        else:
            dir='.'

    name = wlse_basename(exposurename, ccd, ftype_use, 
                         serun=serun,
                         mohrify=mohrify, fext=fext)
    path=path_join(dir,name)
    return path


def wlse_read(exposurename, ccd, ftype, 
              serun=None,
              mohrify=False, 
              dir=None, 
              rootdir=None, 
              fext=None,
              header=False, 
              ext=0):

    fpath=wlse_path(exposurename, ccd, ftype, 
                    serun=serun,
                    mohrify=mohrify, 
                    dir=dir, 
                    rootdir=rootdir, 
                    fext=fext)

    return esutil.io.read(fpath, header=header, ext=ext)

def generate_se_filenames(exposurename, ccd, serun=None,
                          rootdir=None, dir=None):
    fdict={}

    rc=deswl.files.Runconfig()

    # output file names
    for ftype in rc.se_filetypes:

        name= wlse_path(exposurename, 
                        ccd, 
                        ftype,
                        serun=serun,
                        dir=dir,
                        rootdir=rootdir)
        fdict[ftype] = name


    return fdict


def wlse_coldir(serun):
    #dir=run_dir('wlbnl',serun)
    dir = wlse_collated_dir(serun)
    dir=os.path.join(dir,serun+'.cols')
    return dir

def wlse_collated_dir(serun):
    dir=run_dir('wlbnl',serun)
    dir = path_join(dir, 'collated')
    return dir


def wlse_collated_path(serun, 
                       objclass, 
                       ftype=None, 
                       delim=None,
                       region=None,
                       dir=None):
    """
    Can add more functionality later
    """
    fname=[serun,objclass]

    # add a region identifier
    if region is not None:
        if isinstance(region, list):
            regstr=[str(r) for r in region]
            regstr = '-'.join(regstr)
        else:
            regstr=str(region)
        addstr='region%s' % regstr
        fname.append(addstr)

    fname='-'.join(fname)


    # determine the file type
    if ftype is None:
        rc=Runconfig()
        if objclass in rc.se_collated_filetypes:
            ftype=rc.se_collated_filetypes[objclass]
        else:
            typ='fits'
        
    # determine the extension
    fext = esutil.io.ftype2fext(ftype)


    # add delimiter info to the file name
    if fext == 'rec' and delim is not None:
        if delim == '\t':
            fname += '-tab'
        elif delim == ',':
            fname += '-csv'
        else:
            raise ValueError,'delim should be , or tab'

    fname += '.'+fext

    if dir is None:
        dir = wlse_collated_dir(serun)
    outpath = path_join(dir,fname)

    return outpath


def wlse_collated_read(serun, 
                       objclass, 
                       ftype=None, 
                       delim=None,
                       region=None,
                       dir=None,
                       ext=0,
                       header=False,
                       rows=None,columns=None,fields=None,
                       norecfile=False, verbose=False):

    fpath=wlse_collated_path(serun, objclass, ftype=ftype, delim=delim, 
                             region=region, dir=dir)

    return esutil.io.read(fpath, header=header, 
                          rows=rows, columns=columns, fields=fields,
                          norecfile=norecfile, verbose=verbose, 
                          ext=ext) 


def wlme_dir(merun, tilename, rootdir=None):
    rc=Runconfig()

    fileclass=rc.run_types['wlme']['fileclass']
    rundir=run_dir(fileclass, merun, rootdir=rootdir)
    wldir=os.path.join(rundir, tilename)
    return wldir

def wlme_basename(tilename, band, ftype, merun=None, extra=None, ext=None):
    """
    basename for multi-epoch weak lensing output files
    """


    if ext is None:
        rc=Runconfig()
        if ftype in rc.me_fext:
            ext=rc.me_fext[ftype]
        else:
            ext='.fits'

    name=[tilename,band,ftype]
    
    if merun is not None:
        name=[merun] + name

    if extra is not None:
        name.append(extra)

    name='_'.join(name)+ext
    return name


def wlme_path(tilename, band, ftype, merun=None, 
              extra=None, 
              dir=None, rootdir=None, ext=None):
    """
    Return a multi-epoch shear output file name

    currently just for multishear and logs, etc 
    """


    if dir is None:
        if merun is not None:
            dir = wlme_dir(merun, tilename, rootdir=rootdir)
        else:
            dir='.'

    name = wlme_basename(tilename, band, ftype, 
                         merun=merun, extra=extra, ext=ext)
    name = path_join(dir, name)
    return name


def generate_me_filenames(tilename, band, 
                          merun=None, dir=None, rootdir=None):

    fdict={}

    rc=deswl.files.Runconfig()

    # output file names
    for ftype in rc.me_filetypes:
        name= wlme_path(tilename,
                        band,
                        ftype,
                        merun=merun,
                        dir=dir,
                        rootdir=rootdir)
        fdict[ftype] = name


    return fdict





def collated_redfiles_dir(dataset):
    desfiles_dir=getenv_check('DESFILES_DIR')
    desfiles_dir = path_join(desfiles_dir,dataset)
    return desfiles_dir

def collated_redfiles_name(dataset, band):
    name=[dataset,'images','catalogs']
    name.append(band)
    return '-'.join(name)+'.json'

def collated_redfiles_path(dataset, band):
    tdir=collated_redfiles_dir(dataset)
    name = collated_redfiles_name(dataset, band)
    return path_join(tdir, name)

_redfiles_cache={'dataset':None,'band':None,'data':None,'filename':None}
def collated_redfiles_read(dataset, band, getpath=False):
    if (_redfiles_cache['dataset'] == dataset 
            and _redfiles_cache['band'] == band):
        stdout.write('Re-using redfiles cache\n')
        f=_redfiles_cache['filename']
        tileinfo=_redfiles_cache['data']
    else:
        f=collated_redfiles_path(dataset, band)
        stdout.write('Reading tileinfo file: %s\n' % f)
        tileinfo=json_util.read(f)

        _redfiles_cache['filename'] = f
        _redfiles_cache['data'] = tileinfo
        _redfiles_cache['dataset'] = dataset
        _redfiles_cache['band'] = band

    if getpath:
        return tileinfo, f
    else:
        return tileinfo

def collated_redfiles_write(dataset, band, flist):
    f=collated_redfiles_path(dataset, band)
    fdir=os.path.dirname(f)
    if not os.path.exists(fdir):
        stdout.write("Creating directory tree: %s\n" % fdir)
        os.makedirs(fdir)

    output={'rootdir':deswl.files.des_rootdir(),
            'hostname': platform.node(),
            'flist':flist}

    stdout.write('Writing tileinfo file: %s\n' % f)
    json_util.write(output, f)
    stdout.write("    Don't forget to check it into SVN!!!!!")






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

    return '-'.join(name)+'.json'

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
                             getpath=False):
    f=collated_coaddfiles_path(dataset, band, 
                               withsrc=withsrc, 
                               serun=serun,
                               localid=localid, 
                               newest=newest)
    stdout.write('Reading tileinfo file: %s\n' % f)
    tileinfo=json_util.read(f)
    if getpath:
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



def pbs_dir(run, subdir=None):
    outdir=path_join('~','pbs',run)
    if subdir is not None:
        outdir=path_join(outdir, subdir)
    outdir=os.path.expanduser(outdir)
    return outdir


def wlme_pbs_name(tilename, band):
    pbsfile=[tilename,band]
    pbsfile='-'.join(pbsfile)+'.pbs'
    return pbsfile
    
def wlme_pbs_path(merun, tilename, band):
    pbsdir=pbs_dir(merun)
    pbsfile=wlme_pbs_name(tilename,band)
    pbsfile=path_join(pbsdir, pbsfile)
    return pbsfile


def wlse_pbs_name(exposurename, typ='fullpipe', ccd=None):
    pbsfile=[exposurename]

    if typ != 'fullpipe':
        pbsfile.append(typ)

    if ccd is not None:
        pbsfile.append('%02i' % int(ccd))
    pbsfile='-'.join(pbsfile)+'.pbs'
    return pbsfile
    
def wlse_pbs_path(serun, exposurename, typ='fullpipe', ccd=None):
    if ccd is not None:
        subdir='byccd'
    else:
        subdir=None

    pbsdir=pbs_dir(serun, subdir=subdir)
    pbsfile=wlse_pbs_name(exposurename, typ=typ, ccd=ccd)
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
            tmp=(len(fsplit), nparts, filepath)
            mess="Found too few file elements, %s instead of %s: '%s'\n" % tmp
            raise ValueError(mess)

        if fsplit[1] != 'red' or fsplit[3] != 'red':
            mess="In file path %s expected 'red' in positions "+\
                "1 and 3 after DESDATA=%s" % (filepath,DESDATA)
            raise ValueError(mess)

        odict['redrun'] = fsplit[2]
        odict['exposurename'] = fsplit[4]
        odict['basename'] = fsplit[5]

        base=os.path.basename(odict['root'])

        if base[-4:] == '_cat':
            base=base[:-4]

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
            raise ValueError('Expected "decam-" at beginning')
        odict['pointing'] = base.replace('decam-','')

        
    elif fileclass == 'coadd':
        nparts = 5
        if len(fsplit) < nparts:
            mess= 'Found too few file elements, '+\
                '%s instead of %s: "%s"\n' % (len(fsplit), nparts, filepath)
            raise ValueError(mess)

        if fsplit[1] != 'coadd' or fsplit[3] != 'coadd':
            mess="In file path %s expected 'red' in positions "+\
                "1 and 3 after DESDATA=%s" % (filepath,DESDATA)
            raise ValueError(mess)


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
            mess='Expected t find an underscore in basename of %s' % filepath
            raise ValueError(mess)
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




