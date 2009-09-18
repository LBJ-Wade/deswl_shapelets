"""
    Order of doing things for SE image processing:
        # create a run name and its configuration
        GenerateSeRunConfig()

        # For a given run, band generate an input file list.  Also split
        # it up into chunks for each processor
        # also write out a python script to run each chunk.  
        # these go under $DESDATA/runconfig/runname/images-cats-$band-split 
        
        GenerateSeRunFileListAndScripts(run, band, rootdir=default_rootdir, badlist=None)

        # each script calls This
        # run through image list and execute the code
        ProcessSeImageList(serun, executable, 
                          config_file=default_config, 
                          imlist_file=None,
                          flist=None,
                          badlist_file=None,
                          outdir=None,
                          rootdir=None,
                          verbose=True,
                          mohrify=True):

        # you can run a set of scripts for the 8 processors on each machine
        # with this shell script. Run in the directory with the scripts.
        run-wlpipe machine_number executable
        (need to convert to batch system)


    For coadds, using multishear:
        On a machine that can directly access the desdm database:

            import deswl
            deswl.oracle.CollateCoaddCatalogsAndImages(band,collated_file,
                                                       getsrc=True)

        To process data on your local machien, you need to find the files on
        local disk:

            des.wlpipe.FindCollatedCoaddFiles(collated_file, toxml=outfile,
                                              wlserun=run_for_fitpsf_shear)

        For dc4 I'm putting all this stuff into the nyu repo 'desfiles':
            deswl/desfiles/trunk/dc4-coadd

        This data is not really versioned.  The product should be installed and
        set up so that DESFILES_DIR is set (this is also done by setup wl).
        Then the files are under:

            $DESFILES_DIR/dataset/
                (dataset)-images-catalogs-withsrc-(host)-newest-(band).xml

        Then send this file to 

            des.wlpipe.GenerateMeSrclists(tileinfo,outdir='.')

        Which calls GenerateMeSrclist.  This will put files under the
        given output directory.  Results are usually checked into desfiles under
            $DESFILES_DIR/dataset/multishear-srclist


        This calls the function RunMultishear(tilename,band,...) defined
        in this module.  This requires only the tilename and band to
        be set, and figures out all the input files and srclist files.


        
        Running the code is easiest from the stand-alone wrapper
        multishear-run. If you use eups to manager your code, all you should
        need to run the code is to type the following command, or put it into
        your startup file:

            setup wl

        This sets up wl/desfiles/tmv/ccfits/cfitsio as well.  Then from the
        command line:

            multishear-run [options] tilename band


    Using your own local install of wl:
        Say you have an svn checkout in ~/svn/wl and you want the local
        install to be in ~/exports/wl-work.  Then the following are required
        to run the code:

        Put this into your startup file (e.g. bashrc):
            setup wl -r ~/svn/wl-work

        Each time you re-compile of the C++ code or change the python code,
        also do an install to your local install directory.

            scons install PREFIX=~/svn/wl-work
            multishear-run tilename band


"""

import sys
from sys import stdout,stderr
import os
import subprocess
import time
import signal
import re
import datetime

from deswl import GetWlVersion
from deswl import version
from deswl import oracle

from esutil.ostools import path_join, getenv_check
from esutil import xmltools
from esutil.misc import ptime

# default timeout 5 minutes
default_timeout_se = 300
# default timeout 1 hour
default_timeout_me = 60*60

# weak lensing config file
default_config=os.path.join(os.getenv('WL_DIR'),'etc/wl_dc4.config')

_fileclass='wlbnl'

run_types={}
run_types['wlse'] = \
        {'name':'wlse','fileclass': _fileclass, 'filetype':'se_shapelet'}
run_types['wlme'] = \
        {'name':'wlme','fileclass': _fileclass, 'filetype':'me_shapelet'}


fileclasses = {}
fileclasses['wlse'] = _fileclass


se_executables = ['findstars','measurepsf','measureshear']
se_filetypes = ['stars','psf','fitpsf','shear',
                'findstars_log', 
                'measurepsf_log', 
                'measureshear_log']

def FileclassDir(fileclass, rootdir=None):
    if rootdir is not None:
        DESDATA=rootdir
    else:
        DESDATA=os.getenv('DESDATA')
        if DESDATA is None:
            raise ValueError,'DESDATA environment variable not set'
    return os.path.join(DESDATA, fileclass)
    
def RunDir(fileclass, run, rootdir=None):
    fcdir = FileclassDir(fileclass, rootdir=rootdir)
    return os.path.join(fcdir, run)

def FiletypeDir(fileclass, run, filetype, rootdir=None):
    rundir=RunDir(fileclass, run, rootdir=rootdir)
    return os.path.join(rundir, filetype)

def ExposureDir(fileclass, run, filetype, exposurename, rootdir=None):
    ftdir=FiletypeDir(fileclass, run, filetype, rootdir=rootdir)
    return os.path.join(ftdir, exposurename)

def RedImageName(run, exposurename, ccdnum, fz=True, rootdir=None):
    fileclass='red'
    filetype='red'
    expdir = ExposureDir(fileclass, run, filetype, exposurename, 
                         rootdir=rootdir)
    imagename = '%s_%02i.fits' % (exposurename, ccdnum)
    if fz:
        imagename=imagename+'.fz'
    imagename = os.path.join(expdir, imagename)
    return imagename

def TileDir(coaddrun, rootdir=None):
    """
    Coadds are different than other output files:  The fileclass dir is 
    where the files are!
    """
    fileclass='coadd'
    filetype='coadd'
    tiledir=FiletypeDir(fileclass, coaddrun, filetype, rootdir=rootdir)
    return tiledir

def CoaddImageName(coaddrun, tilename, band, fz=True):
    """
    More frustration:  the path element for "coaddrun"
    contains the tilename but in principle it can be arbitrary, so 
    we can't simply pass in pieces:  coaddrun is treated independently
    """
    tiledir=TileDir(coaddrun)
    fname=tilename+'_'+band+'.fits'
    if fz:
        fname+='.fz'
    return fname

def CoaddCatName(catalogrun, tilename, band):
    """
    More frustration:  the path element for "catalogrun"
    contains the tilename but in principle it can be arbitrary, so 
    we can't simply pass in pieces:  catalogrun is treated independently

    Note often catalogrun is the same as the coaddrun
    """
    tiledir=TileDir(catalogrun)
    fname=tilename+'_'+band+'_cat.fits'
    return fname



def RedFileList(band=None, extra=None, latest=False, return_dict=False):
    """
    Send extra="_cat" to get the catalog list
    """
    fileclass='red'
    filetype='red'

    command="des-list-files -c red -e ''"
    if band is not None:
        command = command + " -b %s" % band
    if extra is not None:
        command = command + " -e %s" % extra

    exit_status, flist, serr = ExecuteWlCommand(command)

    if exit_status != 0:
        raise RuntimeError,\
          'des-list-files: exit status: %s stderr: %s' % (exit_status, serr)
    
    flist=flist.split('\n')
    while flist[-1] == '':
        flist.pop()

    if latest:
        return SelectLatestRedRun(flist, return_dict=return_dict)
    else:
        return flist

def SelectLatestRedRun(image_list, return_dict=False):
    """
    Return the latest versions in file list based on run.  By default
    returns the list, but can also return the info dict keyed by
    exposurename-ccdnum

    seems to work with cat files too
    """
    kinfo={}
    for imfile in image_list:
        #info=GetInfoFromImagePath(imfile)
        info=GetInfoFromPath(imfile, 'red')
        key=info['exposurename']+'-'+str(info['ccdnum'])
        key = '%s-%02i' % (info['exposurename'], info['ccdnum'])
        if key in kinfo:
            stdout.write('found\n')
            if info['redrun'] > kinfo[key]['redrun']:
                s='replacing %s > %s\n' % (info['redrun'],kinfo[key]['redrun'])
                stdout.write(s)
                kinfo[key]['redrun'] = info['redrun']
                kinfo[key]['file'] = imfile
        else:
            kinfo[key]={}
            kinfo[key]['redrun'] = info['redrun']
            kinfo[key]['file'] = imfile

    if return_dict:
        return kinfo
    else:
        keeplist=[]
        for k in kinfo:
            keeplist.append(kinfo[k]['file'])
        return keeplist




def GetMatchedRedImageCatList(band=None, combine=True):
    """
    Find the latest catalogs and images and match them up
    """
    stdout.write('Getting latest image list\n')
    imdict = RedFileList(band=band, latest=True, return_dict=True)
    stdout.write('Getting latest cat list\n')
    catdict = RedFileList(band=band, latest=True, extra='_cat', 
                          return_dict=True)


    stdout.write('Matching %s images with %s catalogs\n' % 
                 (len(imdict),len(catdict)) )


    if not combine:
        imlist=[]
        catlist=[]
        for key in sorted(imdict):
            if key in catdict:
                imlist.append(imdict[key]['file'])
                catlist.append(catdict[key]['file'])
        return imlist, catlist
    else:
        outlist=[]
        for key in sorted(imdict):
            if key in catdict:
                t={'imfile':imdict[key]['file'],'catfile':catdict[key]['file']}
                outlist.append(t)
        stdout.write('Matched %s\n' % len(outlist))
        return outlist



_mohr_map={'stars':'shpltall', 
           'psf':'shpltpsf', 
           'fitpsf':'psfmodel', 
           'shear':'shear'}
def Mohrify(name):
    if name in _mohr_map:
        return _mohr_map[name]
    else:
        return name

def WlSeDir(run, exposurename, rootdir=None):
    fileclass=run_types['wlse']['fileclass']
    filetype=run_types['wlse']['filetype']

    ftdir=FiletypeDir(fileclass, run, filetype, rootdir=rootdir)
    wldir=os.path.join(ftdir, exposurename)
 
    return wldir

def WlMeDir(run, tilename, rootdir=None):
    fileclass=run_types['wlme']['fileclass']
    filetype=run_types['wlme']['filetype']

    ftdir=FiletypeDir(fileclass, run, filetype, rootdir=rootdir)
    wldir=os.path.join(ftdir, tilename)
 
    return wldir


def WlSeName(run, redrun, exposurename, ccdnum, wltype, 
             mohrify=True, rootdir=None, ext='fits'):
    """
    name=WlSeName(run, redrun, exposurename, ccdnum, wltype, 
                  mohrify=True, rootdir=None, ext='fits')
    Return the SE output file name for the given inputs

    Because this is designed for single epoch image processing, the
    resulting directory structure is different than uiuc one running
    from the tiling:

    fileclass/run/filetype/redrun_exposurename/
                                  redrun_exposurename_ccd_wltype.fits
    """
    if mohrify:
        wltype_use=Mohrify(wltype)
    else:
        wltype_use=wltype

    wldir = WlSeDir(run, exposurename, rootdir=rootdir)
    name_front = redrun+'_'+os.path.basename(wldir)
    name = '%s_%s_%02i_%s.%s' % (run,name_front, int(ccdnum), wltype_use, ext)
    fullpath = os.path.join(wldir, name)
    return fullpath

def WlMeName(run, coaddrun, tilename, band, wltype):
    """
    Return a multi-epoch shear output file name
    """
    pass

def RunConfigBaseDir():
    DESDATA=os.getenv('DESDATA')
    return os.path.join(DESDATA, 'runconfig')

def RunConfigDir(run):
    basedir=RunConfigBaseDir()
    return os.path.join(basedir, run)

def RunConfigName(run):
    rdir=RunConfigDir(run)
    return os.path.join(rdir, run+'-config.xml')

def ReadRunConfig(run):
    name=RunConfigName(run)
    runconfig = xmltools.xml2dict(name, noroot=True)
    return runconfig




def GenerateRunName(run_type):
    """
    Not used
    """

    if run_type not in run_types:
        raise ValueError,'Unknown run type: %s.  Must be one of: (%s)' % (run_type, ', '.join(run_types))


    now=datetime.datetime.now()
    run = now.strftime('%Y%m%d%H%M%S')
    run = run_type+run
    return run

def GenerateRunNameOrdered(run_type):

    if run_type not in run_types:
        raise ValueError,"Unknown run type: '%s'.  Must be one of: (%s)" % (run_type, ', '.join(run_types))


    dir=RunConfigBaseDir()
    i=0
    run_name = '%s%06i' % (run_type, i)
    rundir=os.path.join(dir, run_name)
    while os.path.exists(rundir):
        i+=1
        run_name = '%s%06i' % (run_type,i)
        rundir=os.path.join(dir, run_name) 

    return run_name



def GenerateSeRunConfig():
    """
    Generate single epoch wl run configuration
    """
    run_type = 'wlse'
    return GenerateRunConfig(run_type)

def GenerateMeRunConfig():
    """
    Generate multi epoch wl run configuration
    """
    run_type = 'wlme'
    return GenerateRunConfig(run_type)


def GenerateRunConfig(run_type):
    """
    Generate single epoch wl run configuration
    """
    # need to add a tmv version here

    run = GenerateRunNameOrdered(run_type)
    fileclass = run_types[run_type]['fileclass']
    filetype = run_types[run_type]['filetype']


    runconfig_dir=RunConfigDir(run)
    runconfig_name = RunConfigName(run)

    if not os.path.exists(runconfig_dir):
        os.mkdir(runconfig_dir)

    # software versions
    pyvers=GetPythonVersion()
    wlvers=GetWlVersion()
    tmvvers=GetTmvVersion()

    runconfig={'serun':run, 
               'run_type':run_type,
               'fileclass': fileclass,
               'filetype':filetype,
               'pyvers':pyvers, 
               'wlvers':wlvers,
               'tmvvers':tmvvers}

    stdout.write('Writing to file: %s\n' % runconfig_name)
    xmltools.dict2xml(runconfig, runconfig_name, roottag='runconfig')
    return runconfig



def VerifyRunConfig(runconfig):
    """
    Make sure the current values for versions and such match that expected
    by the run configuration
    """

    wlvers=GetWlVersion()
    tmvvers=GetTmvVersion()
    pyvers=GetPythonVersion()

    if runconfig['wlvers'] != wlvers:
        raise ValueError,'wlvers "%s" does not match runconfig "%s"' % (wlvers,runconfig['wlvers'])

    if runconfig['tmvvers'] != tmvvers:
        raise ValueError,'tmvvers "%s" does not match runconfig "%s"' % (tmvvers,runconfig['tmvvers'])

    if runconfig['pyvers'] != pyvers:
        raise ValueError,'pyvers "%s" does not match runconfig "%s"' % (pyvers,runconfig['pyvers'])




default_rootdir='/scratch/esheldon/DES'
def GenerateSeRunFileListAndScripts(run, band, rootdir=default_rootdir, badlist=None):
    """
    This is very generic, you will also want to do do more finely grained 
    selections, especially when debugging.  In that case it may make more
    sense to run individual files than work in the full framework.

    Note badlist is used as an input to the *processing* of this filelist
    later, not in building the file list.
    """

    runconfig_dir=RunConfigDir(run)
    runconfig_file=RunConfigName(run)
    runconfig=xmltools.xml2dict(runconfig_file)

    allfile = os.path.join(runconfig_dir, 'images-cats-'+band+'.xml')
    splitdir=os.path.join(runconfig_dir, 'images-cats-'+band+'-bachsplit')
    
    stdout.write('Getting matched image and catalog lists\n')
    flist = GetMatchedRedImageCatList(band=band)

    stdout.write('Found %s images\n' % len(flist))

    fdict={'files':flist}
    stdout.write('Writing all to: %s\n' % allfile)
    xmltools.dict2xml(fdict, allfile, roottag='filelist')

    stdout.write('Splitting files\n')
    # this splits the file list and writes out new xml files.  Also returns
    # the names of those files
    splitfiles = BachSplitImageCatFileList(allfile)

    # Now write out python scripts
    conf_file=default_config

    for ex in se_executables:
        for f in splitfiles:
            pyfilename=f.replace('.xml','-'+ex+'.py')
            stdout.write('Writing python script: %s\n' % pyfilename)
            pyfile=open(pyfilename,'w')
            pyfile.write('from deswl import wlpipe\n')
            pyfile.write('\n')
            pyfile.write('run="%s"\n' % run)
            pyfile.write('ex="%s"\n' % ex)
            pyfile.write('conf="%s"\n' % conf_file)
            pyfile.write('flist="%s"\n' % f)
            pyfile.write('\n')
            extra=''
            if rootdir is not None:
                extra += ', rootdir="%s"' % rootdir
            if badlist is not None:
                extra +=', badlist_file="%s"' % badlist

            command=\
                'wlpipe.ProcessSeRunImageList(run, ex, conf, flist'+extra+')'
            pyfile.write(command+'\n')
            pyfile.close()


def ExtractImageExposureNames(flist):
    allinfo={}
    for imfile in flist:
        #info=GetInfoFromImagePath(imfile)
        info=GetInfoFromPath(imfile, 'red')
        allinfo[info['exposurename']] = info['exposurename']

    return allinfo.keys()

def _PollSubProcess(pobj, timeout_in):
    # minimum of 2 seconds
    timeout = timeout_in
    if timeout <= 1:
        timeout=2

    i=0
    while i < timeout:
        exit_status = pobj.poll()
        if exit_status is not None:
            break
        time.sleep(1)
        i += 1

    if exit_status is None:
        stderr.write('Process is taking longer than %s seconds.  Ending process\n' % timeout)
        os.kill(pobj.pid, signal.SIGTERM)
        exit_status = 1024
        stdout_ret, stderr_ret = None, None
    else:
        stdout_ret, stderr_ret = pobj.communicate()

    return exit_status, stdout_ret, stderr_ret


bitshift = 8
def ExecuteWlCommand(command, timeout=None, 
                     stdout_file=subprocess.PIPE, 
                     stderr_file=subprocess.PIPE, 
                     shell=True,
                     verbose=False):
    """
    exit_status, stdout_returned, stderr_returned = \
       ExecuteWlCommand(command, 
                        timeout=None, 
                        stdout=subprocess.PIPE, 
                        stderr=subprocess.PIPE, 
                        verbose=False)

    Execute the command with a timeout in seconds
    """

    # the user can send file names, PIPE, or a file object
    if isinstance(stdout_file, str):
        stdout_fileobj = open(stdout_file, 'w')
    else:
        stdout_fileobj=stdout_file

    if isinstance(stderr_file, str):
        stderr_fileobj = open(stderr_file, 'w')
    else:
        stderr_fileobj = stderr_file


    # if a list was entered, convert to a string.  Also print the command
    # if requested
    if isinstance(command, list):
        cmd = ' '.join(command)
        if verbose:
            stderr.write(command[0] + '    \\\n')
            for c in command[1:]:
                stderr.write('    '+c+'    \\\n')
        #print 'actual cmd:',cmd
    else:
        cmd=command
        if verbose:
            stderr.write('%s\n' % cmd)



    pobj = subprocess.Popen(cmd, 
                            stdout=stdout_fileobj, 
                            stderr=stderr_fileobj, 
                            shell=shell)

    if timeout is not None:
        exit_status, stdout_ret, stderr_ret = _PollSubProcess(pobj, timeout)
    else:
        # this just waits for the process to end
        stdout_ret, stderr_ret = pobj.communicate()
        # this is not set until we call pobj.communicate()
        exit_status = pobj.returncode

    if exit_status > 100:
        if exit_status > 100:
            mess='This system is probably returning bit-shifted exit codes. Old value: %d\n' % exit_status

            if verbose:
                stderr.write(mess)
            exit_status = exit_status >> bitshift

    # If they were opened files, close them
    if isinstance(stdout_fileobj, file):
        stdout_fileobj.close()
    if isinstance(stderr_fileobj, file):
        stderr_fileobj.close()

    return exit_status, stdout_ret, stderr_ret




def ReplaceDirectory(oldfile, newdir):
    if oldfile is not None:
        fbase=os.path.basename(oldfile)
        newf = os.path.join(newdir, fbase)
        return newf
    else:
        return None




_fits_extstrip=re.compile( '\.fits(\.fz)?$' )
def RemoveFitsExtension(fname):
    """
    See if the pattern exists and replace it with '' if it does
    """
    return _fits_extstrip.sub('', fname)





def MakeSeCommandList(executable, conf_file, image_file, 
                      cat_file=None,
                      stars_file=None, 
                      psf_file=None, 
                      fitpsf_file=None, 
                      shear_file=None):

    root = RemoveFitsExtension(image_file)

    command = \
        [executable,
         conf_file,
         'root='+root]

    if cat_file is not None:
        command.append('cat_file='+cat_file)
    if stars_file is not None:
        command.append('stars_file='+stars_file)
    if psf_file is not None:
        command.append('psf_file='+psf_file)
    if fitpsf_file is not None:
        command.append('fitpsf_file='+fitpsf_file)
    if shear_file is not None:
        command.append('shear_file='+shear_file)

    return command

def ProcessSeImage(executable, config_file, image_filename, 
                   cat_file=None,
                   stars_file=None, 
                   psf_file=None, 
                   fitpsf_file=None, 
                   shear_file=None,
                   stdout_file=None,
                   stderr_file=None,
                   timeout=default_timeout_se,
                   verbose=None):

    command = MakeSeCommandList(executable, config_file, image_filename,
                                cat_file=cat_file,
                                stars_file=stars_file,
                                psf_file=psf_file,
                                fitpsf_file=fitpsf_file,
                                shear_file=shear_file)

    exit_status, stdout_ret, stderr_ret = \
        ExecuteWlCommand(command, timeout=timeout, 
                         stdout_file=stdout_file, 
                         stderr_file=stderr_file,
                         verbose=verbose)

    if verbose:
        stderr.write('    exit_status: %s\n' % exit_status)
    return exit_status







def ProcessSeRunImageList(serun, executable, 
                          config_file=default_config, 
                          imlist_file=None,
                          flist=None,
                          badlist_file=None,
                          outdir=None,
                          rootdir=None,
                          verbose=True,
                          mohrify=True):
    """
    This is for processing one of the big lists such as generated from
    GenerateSeRunFileList()

    The run is required for generating file names and for verifying the
    configuration is as expected
    """

    import time

    tm0 = time.time()

    runconfig = ReadRunConfig(serun)
    VerifyRunConfig(runconfig)

    # flist takes precedence
    if flist is None:
        if imlist_file is None:
            raise ValueError,'Either enter flist or imlist_file'
        imdict = xmltools.xml2dict(imlist_file, noroot=True)
        flist = imdict['files']

    if type(flist) != type([]):
        flist=[flist]

    if badlist_file is not None:
        badlist=xmltools.xml2dict(badlist_file)
        badlist=badlist['filelist']['files']
        badlist = [fd['imfile'] for fd in badlist]
    else:
        badlist=[]

    if verbose:
        stderr.write('Run: %s\n' % serun)

    execbase=os.path.basename(executable)
    for fdict in flist:

        imfile=fdict['imfile']
        catfile=fdict['catfile']

        tm1=time.time()
        if verbose:
            m='Running %s on Image file: %s\n               cat file: %s\n' % \
                    (executable, imfile, catfile)
            stderr.write(m)

        sefdict=GenerateSeFileNames(serun, imfile, 
                                    mohrify=mohrify, 
                                    rootdir=rootdir,
                                    outdir=outdir)

        outd = os.path.dirname(sefdict['stars'])
        if not os.path.exists(outd):
            estatus, stdo, stde = ExecuteWlCommand('mkdir -p '+outd)
            if estatus != 0:
                raise RuntimeError,'Could not make directory: '+outd
            #os.mkdir(outd)
        stdout_file=sefdict[execbase+'_stdout']
        stderr_file=sefdict[execbase+'_stderr']
        log_file=sefdict[execbase+'_log']

        if imfile in badlist:
            if verbose:
                stderr.write('Skipping known bad image: %s\n' % imfile)
            exit_status=-9999
            skipped=True
        else:
            skipped=False
            exit_status = ProcessSeImage(executable, config_file, imfile, 
                                         cat_file=catfile,
                                         stars_file=sefdict['stars'],
                                         psf_file=sefdict['psf'],
                                         fitpsf_file=sefdict['fitpsf'],
                                         shear_file=sefdict['shear'],
                                         stdout_file=stdout_file,
                                         stderr_file=stderr_file,
                                         verbose=verbose)
        # add config info to the log
        log=sefdict
        log['image'] = imfile
        log['catalog'] = catfile
        log['exit_status'] = exit_status
        log['executable'] = executable
        log['skipped'] = skipped
        for k in runconfig:
            log[k] = runconfig[k]

        xmltools.dict2xml(log, log_file, roottag=execbase+'_log')
        tm2=time.time()
        ptime(tm2-tm1, fobj=stderr, format='    %s\n')

    tmfinal=time.time()
    ptime(tmfinal-tm0, fobj=stderr)






def GetExternalVersion(version_command):
    """
    Run a command to get a version string and return it through the
    standard output
    """
    pobj=subprocess.Popen(version_command, shell=True, 
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # this will wait until process terminates
    out=pobj.communicate()
    sout=out[0]
    serr=out[1]
    estatus=pobj.poll()
    if sout == '' or estatus != 0:
        raise RuntimeError, \
         'Could not get version from command: %s: %s' % (version_command,serr)
    return sout.strip()

def GetTmvVersion():
    return GetExternalVersion('tmv-version')

def GetPythonVersion():
    pyvers='v%s.%s.%s' % sys.version_info[0:3]
    return pyvers




def GenerateSeFileNames(serun, imfile, 
                        rootdir=None, outdir=None, mohrify=True):
    #info = GetInfoFromImagePath(imfile)
    info=GetInfoFromPath(imfile, 'red')
    fdict={}

    # output file names
    for wltype in se_filetypes:
        name= WlSeName(serun, 
                       info['redrun'], 
                       info['exposurename'], 
                       info['ccdnum'], 
                       wltype,
                       rootdir=rootdir,
                       mohrify=mohrify)
        if outdir is not None:
            name=ReplaceDirectory(name, outdir)
        fdict[wltype] = name

    # stdout, stderr, and log names
    for executable in se_executables:
        key=executable+'_stdout'
        name = WlSeName(serun, 
                        info['redrun'], 
                        info['exposurename'], 
                        info['ccdnum'], 
                        executable,
                        rootdir=rootdir,
                        ext='stdout')

        if outdir is not None:
            name=ReplaceDirectory(name, outdir)
        fdict[key] = name

        key=executable+'_stderr'
        name = WlSeName(serun, 
                        info['redrun'], 
                        info['exposurename'], 
                        info['ccdnum'], 
                        executable,
                        rootdir=rootdir,
                        ext='stderr')
        if outdir is not None:
            name=ReplaceDirectory(name, outdir)
        fdict[key] = name

        key=executable+'_log'
        name = WlSeName(serun, 
                        info['redrun'], 
                        info['exposurename'], 
                        info['ccdnum'], 
                        executable+'_log',
                        rootdir=rootdir,
                        ext='xml')
        if outdir is not None:
            name=ReplaceDirectory(name, outdir)
        fdict[key] = name



    return fdict




def GetInfoFromImagePath(imfile):
    """
    This will probably also work for cat files
    """

    root = RemoveFitsExtension(imfile)

    DESDATA=os.getenv('DESDATA')
    if DESDATA is None:
        raise ValueError,'Environment variable DESDATA does not exist'
    if DESDATA[-1] == os.sep:
        DESDATA=DESDATA[0:len(DESDATA)-1]

    endf = imfile.replace(DESDATA, '')

    fsplit = endf.split(os.sep)
    if len(fsplit) < 6:
        raise ValueError,'Found too few file elements: "%s"\n' % imfile
    if fsplit[1] != 'red' or fsplit[3] != 'red':
        mess="""In image path %s expected "red" in positions 1 and 3 after DESDATA=%s""" % (imfile,DESDATA)
        raise ValueError,mess

    odict = {}
    odict['filename'] = imfile
    odict['root'] = root
    odict['DESDATA'] = DESDATA
    odict['redrun'] = fsplit[2]
    odict['exposurename'] = fsplit[4]
    odict['basename'] = fsplit[5]

    namefront=RemoveFitsExtension(odict['basename'])
    #odict['ccdnum'] = int(namefront.split('_')[-1])
    odict['ccdnum'] = int(namefront.split('_')[1])
    
    return odict


def GetInfoFromPath(filepath, fileclass):
    """
    This is a more extensive version GetInfoFromImagePath.  Also works for
    the cat files.  If fileclass='coadd' then treated differently
    """

    odict = {}
    # remove the extension
    dotloc=filepath.find('.')
    
    odict['ext'] = filepath[dotloc:]
    root=filepath[0:dotloc]



    DESDATA=os.getenv('DESDATA')
    if DESDATA is None:
        raise ValueError,'Environment variable DESDATA does not exist'
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



def BachSplitImageCatFileList(input_fname, roottag='bachsplit'):
    """

    Split list into files for input to the pipeline, one for each core on each
    bachNN machine

    """

    flistdict, roottag = xmltools.xml2dict(input_fname, seproot=True)
    flist=flistdict['files']

    ntot = len(flist)

    # number of processors
    # number of machines
    nmachine = 4
    cores_per_machine = 8

    ncores = nmachine*cores_per_machine

    npercore = ntot/ncores
    nleft = ntot % ncores

    machine_front = 'bach'

    # get basename without extension
    dname=os.path.dirname(input_fname)
    ifb=os.path.basename(input_fname)
    doti=ifb.find('.')
    if doti == -1:
        raise ValueError,'filehas no extension'

    ifb=ifb[0:doti]

    outdir=os.path.join(dname, ifb+'-split')
    if not os.path.exists(outdir):
        os.mkdir(outdir)


    i=0
    minlist=[]
    maxlist=[]
    filenames=[]
    for imach in range(nmachine):
        for icore in range(cores_per_machine):
            name=ifb+'-'+machine_front + '%02i-%02i.xml' % (imach, icore)
            name = os.path.join(outdir, name)
            filenames.append(name)
            stdout.write(name+'\n')
            stdout.write('    %i : %i\n' % (i*npercore,(i+1)*npercore))

            minlist.append(i*npercore)
            maxlist.append((i+1)*npercore)
            i+=1
    if nleft > 0:
        stdout.write('    and with remainder: %i\n' % (i*npercore+nleft) )
        maxlist[-1] = i*npercore+nleft

    for i in range(ncores):
        flist_sub = flist[minlist[i]:maxlist[i]]
        fdict={'input_file':os.path.abspath(input_fname),
               'files':flist_sub}
        xmltools.dict2xml(fdict, filenames[i], roottag=roottag)

    return filenames









def _find_wl_output(wlserun, src, typename, mohrify=True, verbosity=1):
    fname=WlSeName(wlserun, 
                   src['run'],
                   src['exposurename'],
                   int(src['ccd']),
                   typename,
                   mohrify=mohrify)
    if not os.path.exists(fname):
        if verbosity > 1:
            stdout.write("file not found: %s\n" % fname)
        return False
    else:
        src[typename+'file'] = fname
        return True


def FindCollatedCoaddFiles(infofile, xmlfile, 
                           wlserun=None, mohrify=True,
                           verbosity=1):
    """

    Find coadd files on the local disk.

    Also if srclist is in the dictionaries these files are also searched for
    as well as the corresponding fitpsf file needed for input into the
    multishear code, and the shear catalog for SE.  For the fitpsf files you
    must enter wlserun.

    """
    infodict = xmltools.xml2dict(infofile, noroot=True)
    allinfo=infodict['info']

    DESDATA=os.getenv('DESDATA')


    n_srcmissing=0
    n_srcshort=0
    n_catmissing=0
    n_imagemissing=0
    output = []

    if 'srclist' in allinfo[0]:
        stdout.write("Checking for sources\n")
        checksources=True
    else:
        checksources=False

    for info in allinfo:
        
        catfile=os.sep.join([DESDATA,'coadd',info['catalogrun'],
                             'coadd',info['catalogname']])
        imfile=os.sep.join([DESDATA,'coadd',info['coaddrun'],
                            'coadd',info['coaddname']])
        #catfile = os.path.join(DESDATA, catfile)
        #imfile = os.path.join(DESDATA, imfile)

        srcfound=True
        catfound=True
        imfound=True
        if not os.path.exists(imfile):
            imfile = imfile + '.fz'
            if not os.path.exists(imfile):
                stdout.write('image file not found: %s\n' % imfile)
                imfound=False
                n_imagemissing += 1

        if not os.path.exists(catfile):
            stdout.write('catalog file not found: %s\n' % catfile)
            catfound=False
            n_catmissing += 1

        srclist=[]
        if checksources and imfound and catfound:
            for src in info['srclist']:
                keepsrc = True
                sparts = [DESDATA,
                          'red',
                          src['run'],
                          'red',
                          src['exposurename'],
                          src['filename']]
                srcfile=os.sep.join(sparts)
                if not os.path.exists(srcfile):
                    srcfile += '.fz'
                    if not os.path.exists(srcfile):
                        keepsrc=False
                        stdout.write("srcfile not found: %s\n" % srcfile)
                    else:
                        src['srcfile'] = srcfile
                else:
                    src['srcfile'] = srcfile


                # check for the fitted psf file and shear file
                if wlserun is not None:
                    for ftype in ['fitpsf','shear']:
                        res=_find_wl_output(wlserun,src,ftype,mohrify=mohrify,
                                            verbosity=verbosity)
                        if not res:
                            keepsrc=False
                
                if keepsrc:
                    srclist.append(src)

            if len(srclist) < len(info['srclist']):
                stdout.write("Found fewer sources than listed %s/%s\n" % (len(srclist), len(info['srclist'])))
                n_srcshort += 1
                if len(srclist) == 0:
                    stdout.write("skipping the rest for this coadd %s\n" % imfile)
                    srcfound=False
                    n_srcmissing += 1

        if catfound and imfound and \
                ( (checksources and len(srclist) > 0) or (not checksources)):
            info['catfile'] = catfile
            info['imfile'] = imfile
            if checksources:
                info['srclist'] = srclist
            output.append(info)

    stdout.write("\nout of %s coadds\n" % len(allinfo))
    stdout.write("  n_imagemissing = %s\n" % n_imagemissing)
    stdout.write("  n_catmissing = %s\n" % n_catmissing)
    stdout.write("  n_srcmissing = %s\n" % n_srcmissing)
    stdout.write("  Of the good ones, n_srcshort = %s\n" % n_srcshort)
    stdout.write("  (may not be unique)\n")

    stdout.write("Writing to xml file: %s\n" % xmlfile)
    xmltools.dict2xml({'info':output}, xmlfile, roottag='tileinfo_local')


def GetLatestCollatedCoaddFiles(infofile, latest_file):
    """
    Take an output file from FindCollatedCoaddFiles and return only the
    latest version of each tile (by coaddrun-catalogrun) 
    """


    infodict,root = xmltools.xml2dict(infofile, seproot=True)
    # get the list of info
    allinfo=infodict['info']

    tmpdict={}
    for info in allinfo:
        tilename=info['tilename']
        if tilename not in tmpdict:
            tmpdict[tilename] = info
        else:
            # If this is newer, copy it over the existing one
            oldinfo=tmpdict[tilename]
            oldkey = oldinfo['coaddrun']+'-'+oldinfo['catalogrun']
            newkey = info['coaddrun']+'-'+info['catalogrun']
            if newkey > oldkey:
                stdout.write("Replacing '%s' with '%s'\n" % (oldkey,newkey))
                tmpdict[tilename] = info

    output = {'info': tmpdict.values()}

    stdout.write('Finally keeping %s/%s\n' % (len(tmpdict), len(allinfo)))
    stdout.write('Writing to file: %s\n' % latest_file)
    xmltools.dict2xml(output, latest_file, roottag=root)

def GenerateMeSrclist(tileinfo, outfile):
    """
    Generate an input file for multishear for the input tileinfo dict, which
    should describe a single tile version.  

    This just writes out what is stored in the tile info structure
    (found_coaddinfo).
    """


    fobj=open(outfile, 'w')
    for s in tileinfo['srclist']:
        line = '%s %s %s\n' % (s['srcfile'],s['shearfile'],s['fitpsffile'])
        fobj.write(line)

    fobj.close()

def GenerateMeSrclists(tileinfo_found_latest, outdir='.'):
    """
    The tileinfo_found_latest is the output file generated by
        GetLatestCollatedCoaddFiles
    """
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    tileinfo_dict = xmltools.xml2dict(tileinfo_found_latest)
    tileinfo = tileinfo_dict['tileinfo_local']['info']

    for ti in tileinfo:
        name=[ti['tilename'],ti['band'],'srclist']
        name='-'.join(name)+'.dat'

        outfile = os.path.join(outdir, name)
        stdout.write('Writing multishear input file: %s\n' % outfile)
        GenerateMeSrclist(ti, outfile)

    
    readme_fname=os.path.join(outdir, 'README')
    stdout.write('Creating %s\n' % readme_fname)
    readme=open(readme_fname, 'w')
    mess="""
    Input file list: 
        %s
    Which should be a unique list of tiles.  When duplicate tiles are present
    the latest based on coaddrun-catalogrun was chosen.
    """ % tileinfo_found_latest
    readme.write(mess)
    readme.close()



def MakeMeCommandList(executable, 
                      config_file, 
                      coadd_srclist, 
                      coaddimage_file, 
                      coaddcat_file, 
                      multishear_file):
    command = [executable,
               config_file,
               'coadd_srclist='+coadd_srclist,
               'coaddimage_file='+coaddimage_file,
               'coaddcat_file='+coaddcat_file,
               'multishear_file='+multishear_file]
    return command

def ProcessMeTile(executable, 
                  config_file, 
                  coadd_srclist,
                  coaddimage_file, 
                  coaddcat_file,
                  multishear_file,
                  stdout_file=None,
                  stderr_file=None,
                  timeout=default_timeout_me,
                  verbose=False):
    """
    All paths must be sent explicitly.  Use the wrapper RunMultishear to
    specify a smaller subset
    """

    command = MakeMeCommandList(executable, 
                                config_file, 
                                coadd_srclist, 
                                coaddimage_file, 
                                coaddcat_file, 
                                multishear_file)

    exit_status, stdout_ret, stderr_ret = \
        ExecuteWlCommand(command, timeout=timeout, 
                         stdout_file=stdout_file, 
                         stderr_file=stderr_file,
                         verbose=verbose)

    if verbose:
        stderr.write('    exit_status: %s\n' % exit_status)
    return exit_status



def RunMultishear(tilename, band, 
                  outdir='.', 
                  dataset='dc4coadd', 
                  host='tutti',
                  srclist=None):
    """

    Wrapper function requiring minimial info to process a tile, namely the
    tilename and the band.  The latest version of the tile is found and used.
    All paths are generated.  The output file is placed in '.' unless
    specified via the optional 'outdir' parameter.

    For this code to work, you need to setup the products WL and DESFILES,
    which sets paths and environment variables such as $WL_DIR.  Note WL
    automaticaly sets up DESFILES.  The information for the tiles is kept in
    the file

        $DESFILES_DIR/(dataset)/
           (dataset)-images-catalogs-withsrc-(host)-newest-(band).xml

    This contains a list of the newest versions of all the unique tiles.  The
    exececutable file for multishear and the config file are under $WL_DIR/bin
    and $WL_DIR/etc.  Note the variables 'dataset' and 'host', which default
    to dataset='dc4coadd' and host='tutti'.  The host may be relevant since
    different sites may contain different newest versions.

    The input source list is in 
        $DESFILES_DIR/(dataset)/multishear-srclists/
           (tilename)-(band)-srclist.dat

    You can override the srclist file with the srclist= keyword.

    Example using the library:
        deswl.wlpipe.RunMultishear('DES2215-2853', 'i')\n\n"""

    tm1=time.time()
    wl_dir=getenv_check('WL_DIR')
    desfiles_dir=getenv_check('DESFILES_DIR')

    # the executable
    executable=path_join(wl_dir, 'bin','multishear')

    # config file
    config=path_join(wl_dir, 'etc','wl.config')

    # source list can be determined from tilename/band
    if srclist is None:
        srclist=[tilename,band,'srclist']
        srclist='-'.join(srclist)+'.dat'
        srclist=path_join(desfiles_dir, dataset,'multishear-srclists', srclist)

    # info on all the unique, latest tiles and catalogs
    tileinfo_file=\
        '%s-images-catalogs-withsrc-%s-newest-%s.xml' % (dataset,host,band)

    tileinfo_file=path_join(desfiles_dir,dataset,tileinfo_file)


    # Print some info
    stdout.write('Configuration:\n')
    stdout.write('    dataset: %s\n' % dataset)
    stdout.write('    host: %s\n' % host)
    stdout.write('    WL_DIR: %s\n' % wl_dir)
    stdout.write('    DESFILES_DIR: %s\n' % desfiles_dir)
    stdout.write('    executable: %s\n' % executable)
    stdout.write('    tileinfo file: %s\n' % tileinfo_file)
    stdout.write('    tilename: %s\n' % tilename)
    stdout.write('    band: %s\n' % band)
    stdout.write('    srclist: %s\n' % srclist)

    stdout.write('\nReading tileinfo file\n')
    tileinfo=xmltools.xml2dict(tileinfo_file,noroot=True)
    tileinfo=tileinfo['info']

    # Here we rely on the fact that we picked out the unique, latest
    # version for each tile
    coaddimage=None
    for ti in tileinfo:
        if ti['tilename'] == tilename:
            coaddimage=ti['imfile']
            coaddcat=ti['catfile']
            break

    if coaddimage is None:
        raise ValueError,"Tilename '%s' not found in '%s'" % \
                (tilename, tileinfo_file)

    multishear_file=\
        RemoveFitsExtension(coaddimage)+'_multishear.fits'
    multishear_file=\
        path_join(outdir,os.path.basename(multishear_file))
    stdout_file=multishear_file.replace('.fits','.stdout')
    stderr_file=multishear_file.replace('.fits','.stderr')


    ProcessMeTile(executable, 
                  config, 
                  srclist, 
                  coaddimage, 
                  coaddcat, 
                  multishear_file, 
                  stdout_file=stdout_file,
                  stderr_file=stderr_file,
                  verbose=True)
    tm2=time.time()
    ptime(tm2-tm1, format='RunMultishear Execution time: %s\n')


