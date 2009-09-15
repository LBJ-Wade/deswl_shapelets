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
        ProcessSeRunImageList(serun, executable, 
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
        On a machine that can directly access the desdm database, run (from
        ES des.util module:

            deswl.oracle.CollateCoaddCatalogsAndImages(band,getsrc=True,
                                                     toxml=collated_file)
        Then we need to find the files on local disk, filename from above:
            des.wlpipe.FindCollatedCoaddFiles(collated_file, toxml=outfile,
                                              wlserun=run_for_fitpsf_shear)

        For dc4 I'm putting all this stuff into the nyu repo:
            deswl/desfiles/trunk/dc4-coadd
        This is data it is not really versioned.  The product should
        be installed and set up so that DESFILES_DIR is set.  Then the
        files are under, e.g. 
            $DESFILES_DIR/
              dc4-coadd/dc4-collate-coadd-catalogs-images-withsrc-i-found.xml
        Then send this file to 
            des.util.GenerateMeInputFiles(tileinfo,merun,outdir)
        Which calls GenerateMeInputFile.  This will put files under the
        given output directory.
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

# xmltools is a basic requirement, but I'm not putting it into the DES svn
# instead it must be exported from ES's svn and placed into the path somewhere

try:
    import xmltools
except:
    stderr.write('Failed to import xmltools\n')

# default timeout 5 minutes
default_timeout = 300

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
        info=GetInfoFromImagePath(imfile)
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
    runconfig, roottag = xmltools.xml2dict(name, seproot=True)
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
        info=GetInfoFromImagePath(imfile)
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
        stdout, stderr = None, None
    else:
        stdout, stderr = pobj.communicate()

    return exit_status, stdout, stderr


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
    pobj = subprocess.Popen(command, 
                            stdout=stdout_file, 
                            stderr=stderr_file, 
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





def MakeSeCommand(executable, conf_file, image_file, 
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
                   timeout=default_timeout,
                   verbose=None):

    cmd = MakeSeCommand(executable, config_file, image_filename,
                        cat_file=cat_file,
                        stars_file=stars_file,
                        psf_file=psf_file,
                        fitpsf_file=fitpsf_file,
                        shear_file=shear_file)
    if verbose:
        stderr.write(cmd[0] + '    \\\n')
        for c in cmd[1:]:
            stderr.write('    '+c+'    \\\n')
        stderr.write('    1> '+stdout_file+'\n    2> '+stderr_file+'\n') 
    command = ' '.join(cmd)

    if stdout_file is not None:
        stdo=open(stdout_file,'w')
    else:
        stdo=None
    if stderr_file is not None:
        stde=open(stderr_file,'w')
    else:
        stde=None

    exit_status, stdout_ret, stderr_ret = \
        ExecuteWlCommand(command, timeout=timeout, 
                         stdout_file=stdo, 
                         stderr_file=stde,
                         verbose=verbose)

    if stdo is not None:
        stdo.close()
    if stde is not None:
        stde.close()
    if verbose:
        stderr.write('    exit_status: %s\n' % exit_status)
    return exit_status






default_config=os.path.join(os.getenv('WL_DIR'),'etc/wl_dc4.config')

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
        imdict,roottag = xmltools.xml2dict(imlist_file, seproot=True)
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
    info = GetInfoFromImagePath(imfile)
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





def ptime(seconds, fobj=None, format=None):
    from sys import stdout

    min, sec = divmod(seconds, 60.0)
    hr, min = divmod(min, 60.0)
    days, hr = divmod(hr, 24.0)
    yrs,days = divmod(days, 365.0)

    if yrs > 0:
        tstr="%d years %d days %d hours %d min %f sec" % (yrs,days,hr,min,sec)
    elif days > 0:
        tstr="%d days %d hours %d min %f sec" % (days,hr,min,sec)
    elif hr > 0:
        tstr="%d hours %d min %f sec" % (hr,min,sec)
    elif min > 0:
        tstr="%d min %f sec" % (min,sec)
    else:
        tstr="%f sec" % sec

    if format is None:
        format='%s\n'

    if fobj is None:
        stdout.write(format % tstr)
    else:
        fobj.write(format % tstr)


