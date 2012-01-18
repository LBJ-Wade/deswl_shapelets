from sys import stdout,stderr
import platform
import os
import re
import esutil as eu
from esutil import json_util
from esutil.ostools import path_join, getenv_check

import pprint
import deswl

# for separating file elements
file_el_sep = '-'

def get_proc_environ():
    e={}
    e['TMV_DIR']=getenv_check('TMV_DIR')
    e['WL_DIR']=getenv_check('WL_DIR')
    e['ESUTIL_DIR']=getenv_check('ESUTIL_DIR')
    e['DESFILES_DIR']=getenv_check('DESFILES_DIR')
    e['pyvers']=deswl.get_python_version()
    e['esutilvers']=eu.version()
    e['wlvers']=deswl.version()
    e['tmvvers']=deswl.get_tmv_version()
    return e

class Runconfig(dict):
    def __init__(self, run=None):
        
        self.run=run
        self.fileclass='wlbnl'

        self.run_types={}
        self.run_types['se'] = {'name':'se',
                                'fileclass': self.fileclass, 
                                'filetype':'se_shapelet',
                                'runvarname':'serun'}
        self.run_types['me'] = {'name':'me',
                                'fileclass': self.fileclass, 
                                'filetype':'me_shapelet',
                                'runvarname': 'merun'}


        self.fileclasses = {}
        self.fileclasses['se'] = self.fileclass


        self.se_executables = ['findstars','measurepsf','measureshear']

        self.se_filetypes = ['stars', 'stars1','stars2',
                             'psf','psf1','psf2',
                             'fitpsf','fitpsf1','fitpsf2',
                             'shear','shear1','shear2',
                             'qa','stat','debug']
        self.se_fext       = {'stars':  '.fits', 
                              'stars1':  '.fits', 
                              'stars2':  '.fits', 
                              'psf':    '.fits',
                              'psf1':    '.fits',
                              'psf2':    '.fits',
                              'fitpsf': '.fits', 
                              'fitpsf1': '.fits', 
                              'fitpsf2': '.fits', 
                              'shear':  '.fits',
                              'shear1':  '.fits',
                              'shear2':  '.fits',
                              'qa':     '.dat',
                              'debug':'.dat',
                              'stat':   '.json',
                              'checkpsf':'.rec'}
        self.se_collated_filetypes = {'badlist':'json',
                                      'goodlist':'json',
                                      'gal':'fits'}
        self.me_collated_filetypes = {'badlist':'json', 'goodlist':'json'}


        self.me_filetypes = ['multishear','qa','stat']
        self.me_fext       = {'multishear':'.fits',
                              'qa':'.dat',
                              'stat':'.json'}


        self.me_executables = ['multishear']

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
    

    def generate_new_runconfig_name(self, run_type, band, test=False, old=False):
        """
        generates a new name like se013i, checking against names in
        the runconfig directory
        """

        if run_type not in self.run_types:
            mess="Unknown run type: '%s'.  Must be one "+\
                 "of: (%s)" % (run_type, ', '.join(self.run_types))
            raise ValueError(mess)

        if run_type == 'me':
            i=3
        else:
            i=13

        run_name=self._run_name_from_type_number(run_type, band, i, test=test)

        fullpath=self.getpath(run_name)
        while os.path.exists(fullpath):
            i+=1
            run_name=self._run_name_from_type_number(run_type, band, i, test=test)
            fullpath=self.getpath(run_name)

        return run_name

    def _run_name_from_type_number(self, run_type, band, number, test=False, old=False):
        if old:
            name='%(type)s%(num)04i' % {'type':run_type, 'num':number}
        else:
            name='%(type)s%(num)03i%(band)s' % {'type':run_type,
                                                'num':number,
                                                'band':band}
        if test:
            name += 't'
        return name


    def generate_new_runconfig(self, 
                               run_type, 
                               dataset, 
                               band,
                               wl_config,
                               run_name=None,
                               test=False, 
                               dryrun=False, 
                               pyvers=None,
                               esutilvers=None,
                               wlvers=None,
                               tmvvers=None,
                               comment=None,
                               old=False,
                               **extra):
        """
        Class:
            Runconfig
        Name:
            generate_new_runconfig
        Purpose:
            Generate a new run configuration and write into the usual
            place.
        Calling Sequence:
            rc=deswl.files.Runconfig()
            rc.generate_new_runconfig(run_type, dataset, wl_config, keywords)
        Inputs:
            run_type: 
                'se' or 'me'
            dataset: 
                e.g. 'dc5b' or 'dr012'
            band:
                'g','r','i','z','Y'
            wl_config: 
                Location of the weak lensing config you want to use for this
                run. E.g. '$DESFILES_DIR/wl.config/wl05.config' Such
                environment vars like $DESFILES_DIR will be expanded as needed.

        Keywords:
            run_name: 
                e.g. se003i or me007zt.  If not sent, will be generated.
            test: 
                If true, this is a test run. Generated names will be like
                se008t
            dryrun: 
                If true, just show what would have been written to the file.

            pyvers, esutilvers, wlvers, tmvvers: 
                The versions to use.  Default is the current version.

            comment: Add an additional comment.

        """

        # me runs depend on se runs
        if run_type == 'me':
            if 'serun' not in extra:
                raise RuntimeError("You must send serun=something for me")


        if run_name is None:
            run_name=self.generate_new_runconfig_name(run_type, band, test=test,old=old)

        fileclass = self.run_types[run_type]['fileclass']
        filetype = self.run_types[run_type]['filetype']

        # software versions.  Default to whatever is in our environment
        if pyvers is None:
            pyvers=deswl.get_python_version()
        if esutilvers is None:
            esutilvers=eu.version()
        if wlvers is None:
            wlvers=deswl.version()
        if tmvvers is None:
            tmvvers=deswl.get_tmv_version()

        runconfig={'run':run_name, 
                   'run_type':run_type,
                   'band':band,
                   'fileclass': fileclass,
                   'pyvers':pyvers, 
                   'esutilvers': esutilvers,
                   'wlvers':wlvers,
                   'tmvvers':tmvvers,
                   'dataset':dataset,
                   'wl_config':wl_config}

        if comment is not None:
            runconfig['comment'] = comment

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
        such are the same.  This does not check the data set or or
        serun for merun runconfigs.

        """

        vdict = get_proc_environ()

        #for type in vdict:
        for type in ['wlvers','esutilvers','tmvvers','pyvers']:
            vers=self[type]
            cvers=vdict[type]
            if vers != cvers:
                # check for possiblity we are using a local install of trunk
                # declared with -r
                # note, when making the pbs files we automaticaly expand
                # trunk to ~/exports/{prodname}-work so putting in the -r
                # still is usually not necessary
                if cvers == 'trunk' and vers.find('-r') == -1:
                    raise ValueError("current %s '%s' does not match "
                                     "runconfig '%s'" % (type,cvers,vers))


    def load(self, run):
        config = self.read(run)
        for k in config:
            self[k] = config[k]
        self.run=run

    def read(self, run):
        name=self.getpath(run)
        if not os.path.exists(name):
            mess="runconfig '%s' not found for %s\n" % (run, name)
            raise RuntimeError(mess)
        runconfig = json_util.read(name)
        return runconfig


def ftype2fext(ftype_input):
    ftype=ftype_input.lower()

    if ftype == 'fits' or ftype == 'fit':
        return 'fits'
    elif ftype == 'rec' or ftype == 'pya':
        return 'rec'
    elif ftype == 'json':
        return 'json'
    elif ftype == 'yaml':
        return 'yaml'
    elif ftype == 'xml':
        return 'xml'
    elif ftype == 'eps':
        return 'eps'
    elif ftype == 'ps':
        return 'ps'
    elif ftype == 'png':
        return 'png'
    else:
        raise ValueError("Don't know about '%s' files" % ftype)



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

def hdfs_rootdir():
    return 'hdfs:///user/esheldon/DES'



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

def exposure_dir(fileclass, run, filetype, expname, rootdir=None):
    ftdir=filetype_dir(fileclass, run, filetype, rootdir=rootdir)
    return os.path.join(ftdir, expname)

def red_image_path(run, expname, ccd, fz=True, rootdir=None, check=False):
    fileclass='red'
    filetype='red'
    expdir = exposure_dir(fileclass, run, filetype, expname, 
                          rootdir=rootdir)
    imagename = '%s_%02i.fits' % (expname, int(ccd))
    basic_path = os.path.join(expdir, imagename)

    if check:
        # check both with and without .fz
        path=basic_path
        if not os.path.exists(path):
            path += '.fz'
            if not os.path.exists(path):
                raise RuntimeError("SE image not found: %s(.fz)\n" % basic_path)
    else:
        if fz:
            path=basic_path + '.fz'
        
    return path

def red_cat_path(run, expname, ccd, rootdir=None, check=False):
    fileclass='red'
    filetype='red'
    expdir = exposure_dir(fileclass, run, filetype, expname, 
                          rootdir=rootdir)
    imagename = '%s_%02i_cat.fits' % (expname, int(ccd))
    path = os.path.join(expdir, imagename)

    if check:
        if not os.path.exists(path):
            raise RuntimeError("SE catalog not found: %s\n" % path)
        
    return path




def extract_image_exposure_names(flist):
    allinfo={}
    for imfile in flist:
        info=get_info_from_path(imfile, 'red')
        allinfo[info['expname']] = info['expname']

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

def se_dir(run, expname, rootdir=None):
    rc=Runconfig()

    fileclass=rc.run_types['se']['fileclass']
    rundir=run_dir(fileclass, run, rootdir=rootdir)
    dir=path_join(rundir, expname)
    return dir

def se_basename(expname, ccd, ftype,
                  serun=None,
                  mohrify=False, fext=None):

    if fext is None:
        rc=Runconfig()
        if ftype in rc.se_fext:
            fext=rc.se_fext[ftype]
        else:
            fext='.fits'

    el=[expname, '%02i' % int(ccd), ftype]

    if serun is not None:
        el=[serun]+el

    basename=file_el_sep.join(el) + fext
    return basename



def se_path(expname, ccd, ftype, 
              serun=None,
              mohrify=False, 
              dir=None, 
              rootdir=None, 
              fext=None):
    """
    name=se_path(expname, ccd, ftype, serun=None,
                   mohrify=False, rootdir=None, fext=None)
    Return the SE output file name for the given inputs
    """
    if mohrify:
        ftype_use=mohrify_name(ftype)
    else:
        ftype_use=ftype

    if dir is None:
        if serun is not None:
            dir = se_dir(serun, expname, rootdir=rootdir)
        else:
            dir=os.path.abspath('.')

    name = se_basename(expname, ccd, ftype_use, 
                         serun=serun,
                         mohrify=mohrify, fext=fext)
    path=path_join(dir,name)
    return path


def se_read(expname, ccd, ftype, 
              serun=None,
              mohrify=False, 
              dir=None, 
              rootdir=None, 
              fext=None,
              **keys):
              #ext=0, 
              #verbose=False):

    fpath=se_path(expname, ccd, ftype, 
                    serun=serun,
                    mohrify=mohrify, 
                    dir=dir, 
                    rootdir=rootdir, 
                    fext=fext)

    return eu.io.read(fpath, **keys)

def generate_se_filenames(expname, ccd, serun=None,
                          rootdir=None, dir=None,
                          split=False):
    fdict={}

    rc=deswl.files.Runconfig()

    # output file names
    for ftype in rc.se_filetypes:
        if not split:
            if ftype[-1] == '1' or ftype[-1] == '2':
                continue
        name= se_path(expname, 
                      ccd, 
                      ftype,
                      serun=serun,
                      dir=dir,
                      rootdir=rootdir)
        fdict[ftype] = name


    return fdict

def coldir_open(serun):
    import columns
    coldir=coldir(serun)

    cols = columns.Columns(coldir)
    return cols
    

def coldir(run, fits=False):
    dir = collated_dir(run)
    if fits:
        dir=os.path.join(dir,run+'-fits')
    else:
        dir=os.path.join(dir,run+'.cols')
    return dir



def se_test_dir(serun, subdir=None):
    dir=run_dir('wlbnl',serun)
    dir = path_join(dir, 'test')
    if subdir is not None:
        dir = path_join(dir, subdir)

    return dir

def se_test_path(serun, subdir=None, extra=None, ext='fits'):
    """
    e.g. se_test_path('se011it', subdir=['checksg', 'decam--22--44-i-11'],
                        extra='%02d' % ccd, ext='eps')
    """

    fname = [serun]
    if subdir is not None:
        if isinstance(subdir,(list,tuple)):
            fname += list(subdir)
        else:
            fname += [subdir]

    # extra only goes on the name
    if extra is not None:
        if isinstance(extra,(list,tuple)):
            fname += list(extra)
        else:
            fname += [str(extra)]

    fname='-'.join(fname)
    fname += '.'+ext

    dir = se_test_dir(serun, subdir=subdir)
    outpath = path_join(dir,fname)

    return outpath



def collated_dir(run):
    dir=run_dir('wlbnl',run)
    dir = path_join(dir, 'collated')
    return dir




def se_collated_path(serun, 
                     objclass, 
                     ftype=None, 
                     delim=None):

    fname=[serun,objclass]
    fname='-'.join(fname)

    # determine the file type
    if ftype is None:
        rc=Runconfig()
        if objclass in rc.se_collated_filetypes:
            ftype=rc.se_collated_filetypes[objclass]
        else:
            ftype='fits'
        
    # determine the extension
    fext = eu.io.ftype2fext(ftype)


    # add delimiter info to the file name
    if fext == 'rec' and delim is not None:
        if delim == '\t':
            fname += '-tab'
        elif delim == ',':
            fname += '-csv'
        else:
            raise ValueError,'delim should be , or tab'

    fname += '.'+fext

    dir = collated_dir(serun)
    outpath = path_join(dir,fname)

    return outpath

def me_collated_path(merun, 
                     objclass, 
                     ftype=None, 
                     delim=None):

    fname=[merun,objclass]
    fname='-'.join(fname)

    # determine the file type
    if ftype is None:
        rc=Runconfig()
        if objclass in rc.me_collated_filetypes:
            ftype=rc.me_collated_filetypes[objclass]
        else:
            ftype='fits'
        
    # determine the extension
    fext = eu.io.ftype2fext(ftype)


    # add delimiter info to the file name
    if fext == 'rec' and delim is not None:
        if delim == '\t':
            fname += '-tab'
        elif delim == ',':
            fname += '-csv'
        else:
            raise ValueError,'delim should be , or tab'

    fname += '.'+fext

    dir = collated_dir(merun)
    outpath = path_join(dir,fname)

    return outpath




def se_collated_read(serun, 
                       objclass, 
                       ftype=None, 
                       delim=None,
                       dir=None,
                       ext=0,
                       header=False,
                       rows=None,columns=None,fields=None,
                       norecfile=False, verbose=False):

    fpath=se_collated_path(serun, objclass, ftype=ftype, delim=delim)

    return eu.io.read(fpath, header=header, 
                          rows=rows, columns=columns, fields=fields,
                          norecfile=norecfile, verbose=verbose, 
                          ext=ext) 




def collated_redfiles_dir(dataset, html=False):
    desfiles_dir=getenv_check('DESFILES_DIR')
    if html:
        desfiles_dir=path_join(desfiles_dir,'html')
    desfiles_dir = path_join(desfiles_dir,dataset)
    return desfiles_dir

def collated_redfiles_name(dataset, band):
    #name=[dataset,'images','catalogs']
    #name.append(band)
    name=[dataset,'red',band,'info']
    name = '-'.join(name)
    name += '.json'
    return name

def collated_redfiles_path(dataset, band, html=False):
    tdir=collated_redfiles_dir(dataset, html=html)
    name = collated_redfiles_name(dataset, band)
    return path_join(tdir, name)

_redfiles_cache={'dataset':None,'band':None,'data':None,'url':None}
def collated_redfiles_read(dataset, band, getpath=False):
    if (_redfiles_cache['dataset'] == dataset 
            and _redfiles_cache['band'] == band):
        stdout.write('Re-using redfiles cache\n')
        f=_redfiles_cache['url']
        tileinfo=_redfiles_cache['data']
    else:
        f=collated_redfiles_path(dataset, band)
        stdout.write('Reading red info file: %s\n' % f)
        tileinfo=json_util.read(f)

        _redfiles_cache['url'] = f
        _redfiles_cache['data'] = tileinfo
        _redfiles_cache['dataset'] = dataset
        _redfiles_cache['band'] = band

    if getpath:
        return tileinfo, f
    else:
        return tileinfo

def collated_redfiles_write_web(dataset, band):
    json_file=collated_redfiles_path(dataset, band)
    dat_file=collated_redfiles_path(dataset, band, html=True)
    dat_file = dat_file.replace('.json','-pointings.json')

    html_dir = collated_redfiles_dir(dataset, html=True)

    stdout.write("Writing to: %s\n" % dat_file)
    if not os.path.exists(html_dir):
        os.makedirs(html_dir)

    data = eu.io.read(json_file)
    flist = data['flist']

    #DESDATA=data['rootdir']


    # just for getting unique ones
    phist = {}
    # for json output
    pdata = {}
    pdata['band'] = band
    pdata['dataset'] = dataset
    #pdata['DESDATA'] = DESDATA
    pdata['hostname'] = data['hostname']
    pdata['flist'] = []
    for f in flist:

        imfile=f['image_path']#.replace(DESDATA,'$DESDATA')
        catfile=f['cat_path']#.replace(DESDATA,'$DESDATA')

        imdir=os.path.dirname(imfile)
        catdir=os.path.dirname(catfile)

        imid = '{redrun}/red/{expname}'.format(redrun=f['redrun'],
                                               expname=f['expname'])
        if imid not in phist:
            phist[imid] = 1
            pdata['flist'].append( {'imdir':imdir,'catdir':catdir} )


    eu.io.write(dat_file,pdata,verbose=True)
                                    




# 
# This is a coadd info list that also has 'srclist'
def coadd_info_dir(dataset):
    desfiles_dir=getenv_check('DESFILES_DIR')
    desfiles_dir = path_join(desfiles_dir,dataset)
    return desfiles_dir


def coadd_info_url(dataset, srclist=False):
    if srclist:
        name=[dataset,'coadd','srclist']
    else:
        name=[dataset,'coadd','info']
    name = '-'.join(name)+'.pyobj'
    tdir=coadd_info_dir(dataset)
    return path_join(tdir, name)

_coadd_info_cache={'dataset':None,
                   'data':None,
                   'url':None,
                   'has_srclist':None}

def coadd_info_read(dataset, srclist=False, geturl=False):

    if (_coadd_info_cache['dataset'] == dataset 
            and _coadd_info_cache['has_srclist']==srclist):
        stdout.write('Re-using coadd info cache\n')
        url=_coadd_info_cache['url']
        data=_coadd_info_cache['data']
    else:
        url=coadd_info_url(dataset,srclist=srclist)
        print 'Reading coadd info:',url
        data=eu.io.read(url)

        _coadd_info_cache['dataset'] = dataset
        _coadd_info_cache['data'] = data
        _coadd_info_cache['url'] = url
        _coadd_info_cache['has_srclist'] = srclist

    if geturl:
        return data, url
    else:
        return data

def coadd_info_select(dataset, ids=None, tiles=None, bands=None, srclist=False):
    """
    Select from the input id or tile list.  Id and tile can be scalar
    or list/tupe.
    
    id is faster since that is they key, but keys by tile will be added as you
    go.

    """

    if ids is not None:
        return coadd_info_select_ids(dataset, ds, srclist=srclist)
    elif tiles is not None:
        return coadd_info_select_tile(dataset, tiles, bands=None, srclist=srclist)
    else:
        raise ValueError("send either id= or tiles=")

def coadd_info_select_ids(dataset, ids, srclist=False):

    cinfo=coadd_info_read(dataset, srclist=srclist)
    if isinstance(ids,(list,tuple)): 
        output=[]
        for id in ids:
            ci=_extract_coadd_info_id(cinfo, id)
            output.append(ci)
    else:
        output = _extract_coadd_info_id(cinfo, ids)

    return output

def _extract_coadd_info_id(cinfo, id):
    try:
        ci = cinfo[id]
    except KeyError:
        msg='id %s not found in coadd info for dataset %s'
        msg=msg % (id,_coadd_info_cache['dataset'])
        raise KeyError(msg)

    return ci

def coadd_info_select_tile(dataset, tiles, bands=None, srclist=False):

    cinfo=coadd_info_read(dataset, srclist=srclist)
    if isinstance(tiles,(list,tuple)): 
        output=[]
        for tile in tiles:
            ci=_extract_coadd_info_tile(cinfo, tile, bands=bands)
            output.append(ci)
    else:
        output = _extract_coadd_info_tile(cinfo, tiles)

    return output

def _extract_coadd_info_tile(cinfo, tile):
    for key in cinfo:
        ci=cinfo[key]
        if ci['tilename'] == tile:
            # also cache this tilename as key
            cinfo[tile] = ci
            return ci

    msg='tile %s not found in coadd info for dataset %s'
    msg=msg % (tile,_coadd_info_cache['dataset'])
    raise ValueError(msg)

def _extract_coadd_info_tile_band(cinfo, tile, bands):
    if not isinstance(bands,(type,list)):
        bands=[bands]

    for key in cinfo:
        ci=cinfo[key]
        for band in bands:
            if ci['tilename'] == tile and ci['band'] == band:
                # also cache this tilename-band as key
                tk='%s-%s' % (tile,band)
                cinfo[tk] = ci
                return ci

    msg='tile %s band %s not found in coadd info for dataset %s band %s'
    msg=msg % (tile,band,_coadd_info_cache['dataset'])
    raise ValueError(msg)

#
# multishear output files
#

def me_dir(merun, tilename, rootdir=None):
    rc=Runconfig()

    fileclass=rc.run_types['me']['fileclass']
    rundir=run_dir(fileclass, merun, rootdir=rootdir)
    wldir=os.path.join(rundir, tilename)
    return wldir

def me_basename(tilename, band, ftype, merun=None, extra=None, ext=None):
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

    name=file_el_sep.join(name)+ext
    return name

def me_url(tilename, band, ftype, merun=None, 
           extra=None, 
           dir=None, rootdir=None, ext=None):
    """
    Return a multi-epoch shear output file name

    currently just for multishear and logs, etc 
    """


    if dir is None:
        if merun is not None:
            dir = me_dir(merun, tilename, rootdir=rootdir)
        else:
            dir='.'

    name = me_basename(tilename, band, ftype, 
                         merun=merun, extra=extra, ext=ext)
    name = path_join(dir, name)
    return name

def generate_me_output_urls(tilename, band, 
                            merun=None, dir=None, rootdir=None):
    """
    This is the files output by multishear
    """

    fdict={}

    rc=deswl.files.Runconfig()

    # output file names
    for ftype in rc.me_filetypes:
        name= me_url(tilename,
                     band,
                     ftype,
                     merun=merun,
                     dir=dir,
                     rootdir=rootdir)
        fdict[ftype] = name


    return fdict


#
# generates the input image/catalog urls using the database. Generates
# the single epoch # list url using me_seinputs_url
#
# gets the wl config from the runconfig
#
# gets the output files from generate_me_output_urls
#
class MultishearFiles:
    """
    Generate the input file names, output file names, and condor file names
    """
    def __init__(self, merun, conn=None):
        self.merun=merun
        self.rc=Runconfig(self.merun)

        if conn is None:
            import desdb
            self.conn=desdb.Connection()
        else:
            self.conn=conn


    def get_flist(self):
        query="""
        select
            id,tilename,band
        from
            %(release)s_files
        where
            filetype='coadd'
            and band='%(band)s'\n""" % {'release':self.rc['dataset'],
                                        'band':self.rc['band']}
        print query
        res=self.conn.quick(query)
        if len(res) == 1:
            raise ValueError("Found no coadds for release %s "
                             "band %s" % (self.rc['dataset'],self.rc['band']))

        print 'getting all the file sets'
        all_fdicts=[]
        for r in res:
            print '%s-%s' % (r['tilename'],r['band'])
            id=r['id']
            fdict=self.get_files(id=id)
            all_fdicts.append(fdict)
        return all_fdicts

    def get_files(self, 
                  id=None, 
                  tilename=None,
                  dir=None):
        """
        Call with
            f=get_files(id=)
        or
            f=get_files(tilename=)

        Send dir= to pick a different output directory
        """
        import desdb

        # first the coadd image and catalog
        c=desdb.files.Coadd(id=id, 
                            band=self.rc['band'], 
                            dataset=self.rc['dataset'], 
                            tilename=tilename, 
                            conn=self.conn)
        c.load()

        # now get the output files
        outfiles=generate_me_output_urls(c['tilename'], 
                                         self.rc['band'], 
                                         merun=self.merun, 
                                         dir=dir)

        files={}
        files['wl_config'] = os.path.expandvars(self.rc['wl_config'])
        files['image'] = c['image_url']
        files['cat'] = os.path.expandvars(c['cat_url'])
        files['id'] = os.path.expandvars(c['image_id'])
        files['tilename'] = c['tilename']

        srclist = me_seinputs_url(self.merun, c['tilename'], self.rc['band'])
        files['srclist'] = os.path.expandvars(srclist)

        for k in outfiles:
            files[k] = outfiles[k]

        return files

#
#
# these are lists of srcfile,shearfile,fitpsffile
# for input to multishear
#

def me_seinputs_dir(merun):
    rc=Runconfig(merun)
    d=os.path.join('$DESFILES_DIR',
                   rc['dataset'],
                   'multishear-inputs-%s' % merun)
    return d

def me_seinputs_url(merun, tilename, band):
    d=me_seinputs_dir(merun)

    inputs=[tilename,band,merun,'inputs']
    inputs='-'.join(inputs)+'.dat'
    return path_join(d, inputs)

class MultishearSEInputs:
    """
    This is the single epoch inputs, use MultishearInputs to just generate
    the input file names for multishear.
    """
    def __init__(self, merun, conn=None):
        self.merun=merun
        self.rc=Runconfig(self.merun)

        if conn is None:
            import desdb
            self.conn=desdb.Connection()
        else:
            self.conn=conn


    def generate_all_inputs(self):
        """
        Generate all inputs for the input merun
        """
        query="""
        select
            id
        from
            %(release)s_files
        where
            filetype='coadd'
            and band='%(band)s'\n""" % {'release':self.rc['dataset'],
                                       'band':self.rc['band']}
        res=self.conn.quick(query)
        if len(res) == 1:
            raise ValueError("Found no coadds for release %s "
                             "band %s" % (self.rc['dataset'],self.rc['band']))

        for r in res:
            id=r['id']
            self.generate_inputs(id=id)
       
    def generate_inputs(self, id=None, tilename=None):
        import desdb
        if tilename is None and id is None:
            raise ValueError("Send tilename or id")
        if id is None:
            id=self.get_id(tilename)
        elif tilename is None:
            tilename=self.get_tilename(id)

        c=desdb.files.Coadd(id=id, conn=self.conn)
        c.load(srclist=True)

        url=self.get_url(tilename=tilename)
        url=os.path.expandvars(url)
        d=os.path.dirname(url)
        if not os.path.exists(d):
            os.makedirs(d)

        print 'writing me inputs:',url
        with open(url,'w') as fobj:
            for s in c.srclist:
                names=generate_se_filenames(s['expname'],
                                            s['ccd'],
                                            serun=self.rc['serun'])

                line = '%s %s %s\n' % (s['url'],names['shear'],names['fitpsf'])
                fobj.write(line)

    def get_url(self, id=None, tilename=None):
        """
        Get the location of the single-epoch input list
        for the input id or tilename.

        If the url does not exist, and generate=True, the
        data are generated and written.
        """
        
        if tilename is None and id is None:
            raise ValueError("Send tilename or id")
        if tilename is None:
            tilename=self.get_tilename(id)
        url=me_seinputs_url(self.merun, tilename, self.rc['band'])
        return url

    def get_tilename(self, id):
        query="""
        select
            tilename
        from
            coadd
        where
            id=%(id)d\n""" % {'id':id}

        res=self.conn.quick(query)
        if len(res) != 1:
            raise ValueError("Expected a single match for id "
                             "%s, found %s" % (id,len(res)))
        return res[0]['tilename']

    def get_id(self, tile):
        query="""
        select
            id
        from
            %(release)s_files
        where
            filetype='coadd'
            and tilename='%(tile)s'
            and band='%(band)s'\n""" % {'tile':tile,
                                        'band':self.rc['band'],
                                        'release':self.rc['dataset']}
        res=self.conn.quick(query)
        if len(res) != 1:
            raise ValueError("Expected a single match for tilename "
                             "%s, found %s" % (tile,len(res)))
        return res[0]['id']


class MultishearCondorJob(dict):
    """
    You most certainly want to run the script
        generate-me-condor merun
    Instead of use this directly
    """
    def __init__(self, merun, 
                 files=None,
                 id=None, 
                 tilename=None, 
                 nthread=None, 
                 dir=None,
                 conn=None,
                 verbose=False,
                 dryrun=False):
        """
        Note merun implies a dataset which, combined with tilename,
        implies the band.

        Or send files= (get_files from MultishearFiles) for quicker 
        results, no sql calls required.
        """
        if id is None and tilename is None and files is None:
            raise ValueError("Send id= or tilename= or files=")

        self.rc=Runconfig(merun)
        self['dataset'] = self.rc['dataset']
        self['band'] = self.rc['band']

        self['run'] = merun
        self['nthread'] = nthread
        self['verbose']=verbose
        self['dryrun']=dryrun


        if files is not None:
            self['config'] = files
        else:
            mfobj=MultishearFiles(merun, conn=conn)
            self['config'] = mfobj.get_files(id=id, tilename=tilename,dir=dir)

        self['id'] = self['config']['id']
        self['tilename'] = self['config']['tilename']

        #self['config']['nodots']=True
        self['config']['merun']=merun

        self['config_file']=me_config_path(self['run'],self['tilename'],self['band'])
        self['script_file']=me_script_path(self['run'],self['tilename'],self['band'])
        self['condor_file']=me_condor_path(self['run'],self['tilename'],self['band'])

    def write_all(self):
        self.write_config()
        self.write_script()
        self.write_condor()

    def write_condor(self):
        self.make_condor_dir()

        print "writing to condor submit file:",self['condor_file']
        text = self.condor_text()
        if self['verbose'] or self['dryrun']:
            print text 
        if not self['dryrun']:
            with open(self['condor_file'],'w') as fobj:
                fobj.write(text)
        else:
            print "this is just a dry run" 

    def write_script(self):
        self.make_condor_dir()

        print "writing to condor script file:",self['script_file']
        text = self.script_text()
        if self['verbose'] or self['dryrun']:
            print text 
        if not self['dryrun']:
            with open(self['script_file'],'w') as fobj:
                fobj.write(text)

            print "changing mode of file to executable"
            os.popen('chmod 755 '+self['script_file'])
        else:
            print "this is just a dry run" 


    def write_config(self):
        """
        This is the config meaning the list of files and parameters
        for input to multishear, not the wl config
        """
        self.make_condor_dir()

        print "writing to config file:",self['config_file']
        if self['verbose'] or self['dryrun']:
            pprint.pprint(self.mf)

        if not self['dryrun']:
            with open(self['config_file'],'w') as fobj:
                eu.json_util.write(self['config'], fobj)

    def condor_text(self):
        condor_dir = self.condor_dir()
        script_base = me_script_base(self['tilename'],self['band'])
        condor_text="""
Universe        = parallel
Notification    = Error
GetEnv          = True
Notify_user     = esheldon@bnl.gov
+Experiment     = "astro"
Requirements    = (CPU_Experiment == "astro") && (TotalSlots == 12 || TotalSlots == 8)
Initialdir      = {condor_dir}

Executable      = {script_base}.sh
Output          = {script_base}.out
Error           = {script_base}.err
Log             = {script_base}.log

Queue\n""".format(condor_dir=condor_dir, script_base=script_base)

        return condor_text


    def script_text(self):
        rc=self.rc
        wl_load = _make_load_command('wl',rc['wlvers'])
        tmv_load = _make_load_command('tmv',rc['tmvvers'])
        esutil_load = _make_load_command('esutil', rc['esutilvers'])
     
        thread_text=''
        if self['nthread'] is not None:
            thread_text='export OMP_NUM_THREADS=%s' % self['nthread']

        script_text="""#!/bin/bash
source ~astrodat/setup/setup.sh
source ~astrodat/setup/setup-modules.sh
module load use.own
{tmv_load}
{wl_load}
{esutil_load}
module load desfiles

{thread_text}
multishear-run -c {config_file}
""".format(wl_load=wl_load, 
           tmv_load=tmv_load, 
           esutil_load=esutil_load,
           thread_text=thread_text,
           config_file=os.path.basename(self['config_file']))

        return script_text

    def condor_dir(self):
        return condor_dir(self['run'])

    def make_condor_dir(self):
        condor_dir = self.condor_dir()
        if not os.path.exists(condor_dir):
            print 'making output dir',condor_dir
            os.makedirs(condor_dir)


class MultishearWQJob(dict):
    """
    You most certainly want to run the script
        generate-me-wq merun
    Instead of use this directly
    """
    def __init__(self, merun, 
                 files=None,
                 id=None, 
                 tilename=None, 
                 nthread=None, 
                 dir=None,
                 conn=None,
                 groups=None,
                 verbose=False,
                 dryrun=False):
        """
        Note merun implies a dataset which, combined with tilename,
        implies the band.

        Or send files= (get_files from MultishearFiles) for quicker 
        results, no sql calls required.
        """
        if id is None and tilename is None and files is None:
            raise ValueError("Send id= or tilename= or files=")

        self.rc=Runconfig(merun)
        self['dataset'] = self.rc['dataset']
        self['band'] = self.rc['band']

        self['run'] = merun
        self['nthread'] = nthread
        self['verbose']=verbose
        self['dryrun']=dryrun


        if files is not None:
            self['config'] = files
        else:
            mfobj=MultishearFiles(merun, conn=conn)
            self['config'] = mfobj.get_files(id=id, tilename=tilename,dir=dir)

        self['id'] = self['config']['id']
        self['tilename'] = self['config']['tilename']

        #self['config']['nodots']=True
        self['config']['merun']=merun

        self['job_file']=me_wq_path(self['run'],self['tilename'])
        self['config_file']=me_config_path(self['run'],self['tilename'])
        self['log_file']=os.path.basename(self['job_file']).replace('yaml','out')


        if groups is None:
            groups = '[gen3,gen4,gen5]'
        else:
            groups = '['+groups+']'
        self['groups'] = groups

    def write_all(self):
        self.write_config()
        self.write_job_file()

    def write_job_file(self):
        make_wq_dir(self['run'])

        print "writing to wq job file:",self['job_file']
        text = self.job_file_text()
        if self['verbose'] or self['dryrun']:
            print text 
        if not self['dryrun']:
            with open(self['job_file'],'w') as fobj:
                fobj.write(text)
        else:
            print "this is just a dry run" 

    def write_config(self):
        """
        This is the config meaning the list of files and parameters
        for input to multishear, not the wl config
        """
        import yaml
        make_wq_dir(self['run'])

        print "writing to config file:",self['config_file']
        if self['verbose'] or self['dryrun']:
            print yaml.dump(self.mf)

        if not self['dryrun']:
            with open(self['config_file'],'w') as fobj:
                for k in self['config']:
                    # I want it a bit prettier so writing my own, but will have
                    # to be careful
                    kk = k+': '
                    val = self['config'][k]
                    if isinstance(val,bool):
                        if val:
                            val='true'
                        else:
                            val='false'
                    fobj.write('%-15s %s\n' % (kk,val))
                #yaml.dump(self['config'], fobj)


    def job_file_text(self):
        job_name=self['tilename'] + '-'+self['band']
        groups = self['groups']

        rc=self.rc
        wl_load = _make_load_command('wl',rc['wlvers'])
        tmv_load = _make_load_command('tmv',rc['tmvvers'])
        esutil_load = _make_load_command('esutil', rc['esutilvers'])
     
        thread_text='\n'
        if self['nthread'] is not None:
            thread_text='\nexport OMP_NUM_THREADS=%s' % self['nthread']

        text="""
command: |
    hostname
    source ~astrodat/setup/setup.sh
    source ~/.dotfiles/bash/astro.bnl.gov/modules.sh
    {esutil_load}
    {tmv_load}
    {wl_load}
    {thread_text}
    multishear-run -c {config_file} &> {log_file}

group: {groups}
mode: bynode
priority: low
job_name: {job_name}\n""".format(wl_load=wl_load, 
                                 esutil_load=esutil_load,
                                 tmv_load=tmv_load, 
                                 thread_text=thread_text,
                                 config_file=os.path.basename(self['config_file']),
                                 log_file=self['log_file'],
                                 groups=groups, 
                                 job_name=job_name)

        return text

class ShearFiles(dict):
    """
    Currently, hdfs is only used for reading images and catalogs
    Writing is still done to nfs, which should be fine
    """
    def __init__(self, serun, conn=None, fs='nfs'):
        import desdb
        self['run'] = serun
        self['fs'] = fs

        self.rc = Runconfig(self['run'])
        self['dataset'] = self.rc['dataset']

        self.conn=conn
        if self.conn is None:
            self.conn=desdb.Connection()

        self.expnames = None

    def get_expname_query(self):
        query = """
        select 
            distinct(file_exposure_name) 
        from 
            %(release)s_files
        where 
            filetype = 'red'
            and band = '%(band)s'
        order by
            file_exposure_name\n""" % {'release':self.rc['dataset'],
                                       'band':self.rc['band']}
        return query

    def get_expnames(self):
        if self.expnames is None:
            query = self.get_expname_query()
            print query
            curs = self.conn.cursor()
            curs.execute(query)

            self.expnames = [r[0] for r in curs]

            curs.close()

        return self.expnames

    def get_flist(self):
        """
        Get the red info from the database and add in
        the file info for the wl outputs
        """

        import desdb
        infolist=desdb.files.get_red_info(self.rc['dataset'],self.rc['band'])
        for info in infolist:
            fdict=deswl.files.generate_se_filenames(info['expname'],
                                                    info['ccd'],
                                                    serun=self['run'])
            info['wl_config'] = self.rc['wl_config']
            for k,v in fdict.iteritems():
                info[k] = v

            # we often leave $DESDATA etc. in names
            for k,v in info.iteritems():
                if isinstance(v,basestring):
                    info[k] = os.path.expandvars(v)


            info['serun'] = self['run']

        if self['fs'] == 'hdfs':
            self.convert_to_hdfs(infolist)
        return infolist

    def convert_to_hdfs(self, infolist):
        """
        Modifies in place
        """
        nfs_root=des_rootdir()
        hdfs_root=hdfs_rootdir()
        for info in infolist:
            for k,v in info.iteritems():
                if isinstance(v,basestring):
                    info[k] = v.replace(nfs_root,hdfs_root)
            #for k in ['image_url','cat_url','stars','fitpsf','psf','shear','qa','stat','debug']:
            #    info[k] = info[k].replace(nfs_root,hdfs_root)

    def get_flist_old(self):
        expnames=self.get_expnames()
        flist = []
        for expname in expnames:
            for ccd in xrange(1,62+1):
                fdict=deswl.files.generate_se_filenames(expname,ccd,serun=self['run'])
                fdict['expname'] = expname
                fdict['ccd'] = ccd
                flist.append(fdict)
        return flist


class SEWQJob(dict):
    def __init__(self, 
                 serun, 
                 files=None, 
                 groups=None, 
                 verbose=False,
                 debug_level=-1,
                 type='fullpipe'):
        """
        For now just take files list for a given ccd.

        Files should be an element from ShearFiles get_flist(), which returns a
        list of dicts
        """
        if files is None:
            raise ValueError("Send files=")

        self['config'] = files
        self['config']['serun'] = serun
        self['config']['type'] = type
        self['config']['debug_level'] = debug_level
        self['verbose']=verbose

        self['job_file']=se_wq_path(self['config']['serun'],
                                    self['config']['expname'],
                                    type=self['config']['type'],
                                    ccd=self['config']['ccd'])
        self['log_file']=os.path.basename(self['job_file']).replace('yaml','out')
        self['config_file']=se_config_path(self['config']['serun'],
                                           self['config']['expname'],
                                           type=self['config']['type'],
                                           ccd=self['config']['ccd'])


        if groups is None:
            groups = '[new,new2]'
        else:
            groups = '['+groups+']'
        self['groups'] = groups


    def write_all(self):
        self.write_config()
        self.write_job_file()

 
    def write_job_file(self):
        make_wq_dir(self['config']['serun'])

        f=se_wq_path(self['config']['serun'], 
                     self['config']['expname'], 
                     type=self['config']['type'], 
                     ccd=self['config']['ccd'])

        text = self.job_file_text()
        if self['verbose']:
            print "writing to wq job file:",f 

        with open(f,'w') as fobj:
            fobj.write(text)


    def job_file_text(self):

        groups = self['groups']

        ccd=self['config']['ccd']
        expname=self['config']['expname']

        tname = se_wq_path(self['config']['serun'], expname,ccd=ccd,
                           type=self['config']['type'])
        tname = os.path.basename(tname)
        log_name = tname.replace('.yaml','.out')

        rc=Runconfig(self['config']['serun'])
        wl_load = _make_load_command('wl',rc['wlvers'])
        tmv_load = _make_load_command('tmv',rc['tmvvers'])
        esutil_load = _make_load_command('esutil', rc['esutilvers'])
     
        if ccd == None:
            ccd=''
        job_name='%s-%02i' % (expname,ccd)

        config_file=os.path.basename(self['config_file'])
        text = """
command: |
    hostname
    source /opt/astro/SL53/bin/setup.hadoop.sh
    source ~astrodat/setup/setup.sh
    source ~/.dotfiles/bash/astro.bnl.gov/modules.sh
    %(esutil_load)s
    %(tmv_load)s
    %(wl_load)s

    export OMP_NUM_THREADS=1
    shear-run %(config_file)s &> %(log_file)s
    echo Done

group: %(groups)s
priority: low
job_name: %(job_name)s\n""" % {'esutil_load':esutil_load,
                               'tmv_load':tmv_load,
                               'wl_load':wl_load,
                               'config_file':config_file,
                               'log_file':self['log_file'],
                               'groups':groups,
                               'job_name':job_name}


        return text

    def write_config(self):
        """
        This is the config meaning the list of files and parameters
        for input to multishear, not the wl config
        """
        import yaml
        make_wq_dir(self['config']['serun'], subdir='byccd')

        if self['verbose']:
            print "writing to config file:",self['config_file']

        with open(self['config_file'],'w') as fobj:
            for k in self['config']:
                # I want it a bit prettier so writing my own, but will have
                # to be careful
                kk = k+': '
                val = self['config'][k]
                if isinstance(val,bool):
                    if val:
                        val='true'
                    else:
                        val='false'
                fobj.write('%-15s %s\n' % (kk,val))
            #yaml.dump(self['config'], fobj)



class SEWQJobOld(dict):
    def __init__(self, serun, expname, ccd=None, type='fullpipe'):
        self['run'] = serun
        self['expname'] = expname
        self['ccd'] = ccd
        self['type'] = type

    def write_jobfile(self, verbose=False, dryrun=False):
        make_wq_dir(self['run'])

        f=se_wq_path(self['run'], 
                     self['expname'], 
                     type=self['type'], 
                     ccd=self['ccd'])

        print "writing to wq job file:",f 
        text = self.job_file_text()
        if verbose or dryrun:
            print text 
        if not dryrun:
            with open(f,'w') as fobj:
                fobj.write(text)
        else:
            print "this is just a dry run" 


    def job_file_text(self):

        groups = '[gen3,gen4,gen5]'

        ccd=self['ccd']
        tname = se_wq_path(self['run'], self['expname'],ccd=ccd,type=self['type'])
        tname = os.path.basename(tname)
        log_name = tname.replace('.yaml','.out')

        rc=Runconfig(self['run'])
        wl_load = _make_load_command('wl',rc['wlvers'])
        tmv_load = _make_load_command('tmv',rc['tmvvers'])
        esutil_load = _make_load_command('esutil', rc['esutilvers'])
     
        if ccd == None:
            ccd=''

        text = """
command: |
    hostname
    source ~astrodat/setup/setup.sh
    source ~/.dotfiles/bash/astro.bnl.gov/modules.sh
    %(esutil_load)s
    %(tmv_load)s
    %(wl_load)s

    export OMP_NUM_THREADS=1
    shear-run --serun=%(serun)s --nodots %(expname)s %(ccd)s &> %(log)s

group: %(groups)s
priority: low
job_name: %(job_name)s\n""" % {'esutil_load':esutil_load,
                               'tmv_load':tmv_load,
                               'wl_load':wl_load,
                               'serun':self['run'],
                               'expname':self['expname'],
                               'ccd':ccd,
                               'log':log_name,
                               'groups':groups,
                               'job_name':self['expname']}


        return text





def _make_load_command(modname, vers):
    """
    convention is that trunk is exported to ~/exports/modname-work
    """
    load_command = "module unload {modname} && module load {modname}".format(modname=modname)
    if vers == 'trunk' and modname != 'desfiles':
        load_command += '/work'
    else:
        load_command += '/%s' % vers

    return load_command




#
# old pbs scripts for running a tile. These will be used by condor
#



def pbs_dir(run, subdir=None):
    outdir=path_join('~','pbs',run)
    if subdir is not None:
        outdir=path_join(outdir, subdir)
    outdir=os.path.expanduser(outdir)
    return outdir


def me_pbs_name(tilename, band):
    pbsfile=[tilename,band]
    pbsfile='-'.join(pbsfile)+'.pbs'
    return pbsfile
    
def me_pbs_path(merun, tilename, band):
    pbsdir=pbs_dir(merun)
    pbsfile=me_pbs_name(tilename,band)
    pbsfile=path_join(pbsdir, pbsfile)
    return pbsfile


def se_pbs_name(expname, typ='fullpipe', ccd=None):
    pbsfile=[expname]

    if typ != 'fullpipe':
        pbsfile.append(typ)

    if ccd is not None:
        pbsfile.append('%02i' % int(ccd))
    pbsfile='-'.join(pbsfile)+'.pbs'
    return pbsfile
    
def se_pbs_path(serun, expname, typ='fullpipe', ccd=None):
    if ccd is not None:
        subdir='byccd'
    else:
        subdir=None

    pbsdir=pbs_dir(serun, subdir=subdir)
    pbsfile=se_pbs_name(expname, typ=typ, ccd=ccd)
    pbsfile=path_join(pbsdir, pbsfile)
    return pbsfile






#
# new wq stuff
#

def wq_dir(run, subdir=None):
    outdir=path_join('~','des-wq',run)
    if subdir is not None:
        outdir=path_join(outdir, subdir)
    outdir=os.path.expanduser(outdir)
    return outdir

def make_wq_dir(run, subdir=None):

    d = wq_dir(run, subdir=subdir)

    if not os.path.exists(d):
        print 'making output dir',d
        os.makedirs(d)



def se_wq_path(run, expname, type='fullpipe', ccd=None):
    wqfile=[expname]

    if type != 'fullpipe':
        wqfile.append(type)

    if ccd is not None:
        wqfile.append('%02i' % int(ccd))

    wqfile='-'.join(wqfile)+'.yaml'

    if ccd is not None:
        subdir='byccd'
    else:
        subdir=None

    d = wq_dir(run, subdir=subdir)

    wqpath=path_join(d, wqfile)
    return wqpath

def se_config_path(run, expname, type='fullpipe', ccd=None):
    """
    This is the file with all the paths and parameters,
    not a wl_config
    """
    
    f=se_wq_path(run,expname,type=type,ccd=ccd)
    f=f[0:f.rfind('.')]+'-config.yaml'
    return f



def me_wq_path(run, tilename):
    rc = Runconfig(run)
    band = rc['band']


    f=[tilename,band]
    f='-'.join(f)+'.yaml'

    d = wq_dir(run)

    wqpath=path_join(d, f)
    return wqpath

def me_config_path(merun, tilename):
    """
    This is the file with all the paths and parameters,
    not a wl_config
    """
    
    f=me_wq_path(merun, tilename)
    f=f[0:f.rfind('.')]+'-config.yaml'
    return f



#
# condor submit files and bash scripts

def condor_dir(run, subdir=None):
    outdir=path_join('~','condor','wl',run)
    if subdir is not None:
        outdir=path_join(outdir, subdir)
    outdir=os.path.expanduser(outdir)
    return outdir


def me_condor_name(tilename, band):
    condorfile=[tilename,band]
    condorfile='-'.join(condorfile)+'.condor'
    return condorfile
    
def me_condor_path(merun, tilename, band):
    condordir=condor_dir(merun)
    condorfile=me_condor_name(tilename,band)
    condorfile=path_join(condordir, condorfile)
    return condorfile


def se_condor_name(expname, typ='fullpipe', ccd=None):
    condorfile=[expname]

    if typ != 'fullpipe':
        condorfile.append(typ)

    if ccd is not None:
        condorfile.append('%02i' % int(ccd))
    condorfile='-'.join(condorfile)+'.condor'
    return condorfile
    
def se_condor_path(serun, expname, typ='fullpipe', ccd=None):
    if ccd is not None:
        subdir='byccd'
    else:
        subdir=None

    condordir=condor_dir(serun, subdir=subdir)
    condorfile=se_condor_name(expname, typ=typ, ccd=ccd)
    condorfile=path_join(condordir, condorfile)
    return condorfile



def me_script_base(tilename, band):
    script_parts=[tilename,band]
    script_base='-'.join(script_parts)
    return script_base

def me_script_name(tilename, band):
    script_base = me_script_base(tilename, band)
    script_name = script_base+'.sh'
    return script_name
    
def me_script_path(merun, tilename, band):
    scriptdir=condor_dir(merun)
    scriptfile=me_script_name(tilename,band)
    scriptfile=path_join(scriptdir, scriptfile)
    return scriptfile




def se_script_base(expname, typ='fullpipe', ccd=None):
    script_parts=[expname]

    if typ != 'fullpipe':
        script_parts.append(typ)

    if ccd is not None:
        script_parts.append('%02i' % int(ccd))
    script_base='-'.join(script_parts)
    return script_base

def se_script_name(expname, typ='fullpipe', ccd=None):
    script_base = se_script_base(expname, typ=typ, ccd=ccd)
    script_name = script_base+'.sh'
    return script_name
    
def se_script_path(serun, expname, typ='fullpipe', ccd=None):
    if ccd is not None:
        subdir='byccd'
    else:
        subdir=None

    scriptdir=condor_dir(serun, subdir=subdir)
    scriptfile=se_script_name(expname, typ=typ, ccd=ccd)
    scriptfile=path_join(scriptdir, scriptfile)
    return scriptfile












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
        odict['expname'] = fsplit[4]
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




