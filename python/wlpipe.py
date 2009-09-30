"""
    Order of doing things for SE image processing:
        (note this is the old pre-batch system notes)

        # create a run name and its configuration
        rc=deswl.files.Runconfig()
        rc.generate_new_runconfig('wlse')

        # For a given run, band generate an input file list.  Also split
        # it up into chunks for each processor
        # also write out a python script to run each chunk.  
        # these go under $DESDATA/runconfig/runname/images-cats-$band-split 
        
        generate_serun_filelists_scripts(run, band, rootdir=default_scratch_rootdir, badlist=None)

        # each script calls This
        # run through image list and execute the code
        process_serun_imagelist(serun, executable, 
                          config_file=None, 
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
            deswl.oracle.collate_coadd_catim(band,collated_file,
                                             getsrc=True)


        For dc4 I'm putting all this stuff into the nyu repo 'desfiles':
            $DESFILES_DIR/dc4coadd
        Generally
            $DESFILES_DIR/(dataset)/(dataset)-images-catalogs-(band).xml
            $DESFILES_DIR/(dataset)/(dataset)-images-catalogs-withsrc-(band).xml

        This data is not really versioned.  The product should be installed and
        set up so that DESFILES_DIR is set (this is also done by setup wl).

        To process data on your local machine, you need to find the files on
        local disk:

            deswl.wlpipe.find_collated_coaddfiles(collated_file, toxml=outfile,
                                                wlserun=run_for_fitpsf_shear)

        which I'm keeping here:
            $DESFILES_DIR/(dataset)/(localid)/(dataset)-images-catalogs-withsrc-(localid)-(band).xml

        And to get the newest versions, send that file to:
            deswl.wlpipe.get_newest_coaddfiles(infofile, output_file)
        Stored in
            $DESFILES_DIR/(dataset)/(localid)/(dataset)-images-catalogs-withsrc-(localid)-newest-(band).xml

        Then this file is used to generate source lists (SE images that contribute to the 
        coadd tiles) for input to multishear.

            deswl.wlpipe.generate_me_srclists(dataset, band, localid, outdir='.')

        Which calls generate_me_srclist.  This will put files under the
        given output directory.  Results are usually checked into desfiles under
            $DESFILES_DIR/(dataset)/(localid)/multishear-srclist/(tilename)-(band)-srclist.dat

        This calls the function run_multishear(tilename,band,...) defined
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
import platform
import os
import subprocess
import time
import signal
import re
import datetime
import shutil

import deswl

import esutil
from esutil.ostools import path_join, getenv_check
from esutil import xmltools
from esutil.misc import ptime


def default_dataset():
    return 'dc4coadd'
def default_localid():
    return 'tutti2'
def default_serun():
    return 'wlse000001'
def default_scratch_rootdir():
    return '/scratch/esheldon/DES'


def get_red_filelist(band=None, extra=None, latest=False, return_dict=False):
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

    exit_status, flist, serr = execute_command(command)

    if exit_status != 0:
        raise RuntimeError,\
          'des-list-files: exit status: %s stderr: %s' % (exit_status, serr)
    
    flist=flist.split('\n')
    while flist[-1] == '':
        flist.pop()

    if latest:
        return select_latest_red_run(flist, return_dict=return_dict)
    else:
        return flist

def select_latest_red_run(image_list, return_dict=False):
    """
    Return the latest versions in file list based on run.  By default
    returns the list, but can also return the info dict keyed by
    exposurename-ccdnum

    seems to work with cat files too
    """
    kinfo={}
    for imfile in image_list:
        info=deswl.files.get_info_from_path(imfile, 'red')
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




def get_matched_red_image_catlist(band=None, combine=True):
    """
    Find the latest catalogs and images and match them up
    """
    stdout.write('Getting latest image list\n')
    imdict = get_red_filelist(band=band, latest=True, return_dict=True)
    stdout.write('Getting latest cat list\n')
    catdict = get_red_filelist(band=band, latest=True, extra='_cat', 
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

def _run_name_from_type_number(run_type, number, test=False):
    if test:
        name='%stest%06i' % (run_type, number)
    else:
        name='%s%06i' % (run_type, number)
    return name

def generate_run_name(run_type, test=False):
    """
    generates a new name like run_type000002, checking against names in
    runconfig_basedir()
    """

    rc=deswl.files.Runconfig()
    if run_type not in rc.run_types:
        raise ValueError,"Unknown run type: '%s'.  Must be one of: (%s)" % (run_type, ', '.join(rc.run_types))


    dir=deswl.files.runconfig_basedir()
    i=0
    #run_name = '%s%06i' % (run_type, i)
    run_name=_run_name_from_type_number(run_type, i, test=test)
    rundir=os.path.join(dir, run_name)
    while os.path.exists(rundir):
        i+=1
        #run_name = '%s%06i' % (run_type,i)
        run_name=_run_name_from_type_number(run_type, i, test=test)
        rundir=os.path.join(dir, run_name) 

    return run_name



def generate_runconfig(run_type, dataset, localid, test=False):
    """
    Generate wl run configuration. run_types supported currently are
        'wlse': single epoch
        'wlme': multi epoch
    """

    rc=deswl.files.Runconfig()
    run = rc.generate_new_runconfig_path(run_type, test=test)

    fileclass = rc.run_types[run_type]['fileclass']
    filetype = rc.run_types[run_type]['filetype']

    # software versions
    pyvers=deswl.get_python_version()
    esutilvers=esutil.version()
    wlvers=deswl.version()
    tmvvers=deswl.get_tmv_version()

    runconfig={'run':run, 
               'run_type':run_type,
               'fileclass': fileclass,
               'filetype':filetype,
               'pyvers':pyvers, 
               'esutilvers': esutilvers,
               'wlvers':wlvers,
               'tmvvers':tmvvers,
               'dataset':dataset,
               'localid':localid}

    stdout.write('Writing to file: %s\n' % runconfig_name)
    xmltools.dict2xml(runconfig, runconfig_name, roottag='runconfig')
    stdout.write("Dont' forget to check it into SVN!\n")
    return runconfig



def generate_serun_filelists_scripts(run, band, 
                                     rootdir=None, badlist=None):
    """
    This is partially no longer needed since we have a batch system in place


    This is very generic, you will also want to do do more finely grained 
    selections, especially when debugging.  In that case it may make more
    sense to run individual files than work in the full framework.

    Note badlist is used as an input to the *processing* of this filelist
    later, not in building the file list.
    """

    if rootdir is None:
        rootdir=default_scratch_rootdir()
    runconfig_dir=deswl.files.runconfig_dir(run)
    runconfig_file=deswl.files.runconfig_name(run)
    runconfig=xmltools.xml2dict(runconfig_file)

    allfile = os.path.join(runconfig_dir, 'images-cats-'+band+'.xml')
    splitdir=os.path.join(runconfig_dir, 'images-cats-'+band+'-bachsplit')
    
    stdout.write('Getting matched image and catalog lists\n')
    flist = get_matched_red_image_catlist(band=band)

    stdout.write('Found %s images\n' % len(flist))

    fdict={'files':flist}
    stdout.write('Writing all to: %s\n' % allfile)
    xmltools.dict2xml(fdict, allfile, roottag='filelist')

    stdout.write('Splitting files\n')
    # this splits the file list and writes out new xml files.  Also returns
    # the names of those files
    splitfiles = bach_split_imagecat_filelist(allfile)

    # Now write out python scripts
    conf_file=os.path.join(getenv_check('WL_DIR'),'etc/wl_dc4.config')
    rc=deswl.files.Runconfig()
    for ex in rc.se_executables:
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
                'wlpipe.process_serun_imagelist(run, ex, conf, flist'+extra+')'
            pyfile.write(command+'\n')
            pyfile.close()


def _poll_subprocess(pobj, timeout_in):
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
        stdout.write('Process is taking longer than %s seconds.  Ending process\n' % timeout)
        os.kill(pobj.pid, signal.SIGTERM)
        exit_status = 1024
        stdout_ret, stderr_ret = None, None
    else:
        stdout_ret, stderr_ret = pobj.communicate()

    return exit_status, stdout_ret, stderr_ret


bitshift = 8
def execute_command(command, timeout=None, 
                     stdout_file=subprocess.PIPE, 
                     stderr_file=subprocess.PIPE, 
                     shell=True,
                     verbose=False):
    """
    exit_status, stdout_returned, stderr_returned = \
       execute_command(command, 
                        timeout=None, 
                        stdout=subprocess.PIPE, 
                        stderr=subprocess.PIPE, 
                        verbose=False)

    Execute the command with a timeout in seconds
    """

    # the user can send file names, PIPE, or a file object
    if isinstance(stdout_file, str):
        stdout.write("redirecting stdout to %s\n" % stdout_file)
        stdout_fileobj = open(stdout_file, 'w')
    else:
        stdout_fileobj=stdout_file

    if isinstance(stderr_file, str):
        stdout.write("redirecting stderr to %s\n" % stderr_file)
        stderr_fileobj = open(stderr_file, 'w')
    else:
        stderr_fileobj = stderr_file


    # if a list was entered, convert to a string.  Also print the command
    # if requested
    stdout.write('Executing command: \n')
    if isinstance(command, list):
        cmd = ' '.join(command)
        if verbose:
            stdout.write(command[0] + '    \\\n')
            for c in command[1:]:
                stdout.write('    '+c+'    \\\n')
        #print 'actual cmd:',cmd
    else:
        cmd=command
        if verbose:
            stdout.write('%s\n' % cmd)



    stdout.flush()
    stderr.flush()
    pobj = subprocess.Popen(cmd, 
                            stdout=stdout_fileobj, 
                            stderr=stderr_fileobj, 
                            shell=shell)

    if timeout is not None:
        exit_status, stdout_ret, stderr_ret = _poll_subprocess(pobj, timeout)
    else:
        # this just waits for the process to end
        stdout_ret, stderr_ret = pobj.communicate()
        # this is not set until we call pobj.communicate()
        exit_status = pobj.returncode

    if exit_status > 100:
        if exit_status > 100:
            mess='This system is probably returning bit-shifted exit codes. Old value: %d\n' % exit_status

            if verbose:
                stdout.write(mess)
            exit_status = exit_status >> bitshift

    # If they were opened files, close them
    if isinstance(stdout_fileobj, file):
        stdout_fileobj.close()
    if isinstance(stderr_fileobj, file):
        stderr_fileobj.close()

    return exit_status, stdout_ret, stderr_ret








def make_se_commandlist(executable, conf_file, image_file, 
                      cat_file=None,
                      stars_file=None, 
                      psf_file=None, 
                      fitpsf_file=None, 
                      shear_file=None):

    root = deswl.files.remove_fits_extension(image_file)

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

def process_se_image(executable, config_file, image_filename, 
                     cat_file=None,
                     stars_file=None, 
                     psf_file=None, 
                     fitpsf_file=None, 
                     shear_file=None,
                     stdout_file=None,
                     stderr_file=None,
                     timeout=5*60,
                     verbose=None):

    command = make_se_commandlist(executable, config_file, image_filename,
                                  cat_file=cat_file,
                                  stars_file=stars_file,
                                  psf_file=psf_file,
                                  fitpsf_file=fitpsf_file,
                                  shear_file=shear_file)

    exit_status, stdout_ret, stderr_ret = \
        execute_command(command, timeout=timeout, 
                         stdout_file=stdout_file, 
                         stderr_file=stderr_file,
                         verbose=verbose)

    if verbose:
        stdout.write('exit_status: %s\n' % exit_status)
    return exit_status







def process_serun_imagelist(serun, executable, 
                            config_file=None, 
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

    if config_file is None:
        config_file=os.path.join(getenv_check('WL_DIR'),'etc/wl_dc4.config')

    tm0 = time.time()

    runconfig = deswl.files.Runconfig(serun)
    runconfig.verify()

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
        stdout.write('Run: %s\n' % serun)

    execbase=os.path.basename(executable)
    for fdict in flist:

        imfile=fdict['imfile']
        catfile=fdict['catfile']

        tm1=time.time()
        if verbose:
            m='Running %s on Image file: %s\n               cat file: %s\n' % \
                    (executable, imfile, catfile)
            stdout.write(m)

        sefdict=generate_se_filenames(serun, imfile, 
                                      mohrify=mohrify, 
                                      rootdir=rootdir,
                                      outdir=outdir)

        outd = os.path.dirname(sefdict['stars'])
        if not os.path.exists(outd):
            if verbose:
                stdout.write('    Making output directory path: %s\n' % outd)
            os.makedirs(outd)

        stdout_file=sefdict[execbase+'_stdout']
        stderr_file=sefdict[execbase+'_stderr']
        log_file=sefdict[execbase+'_log']

        if imfile in badlist:
            if verbose:
                stdout.write('Skipping known bad image: %s\n' % imfile)
            exit_status=-9999
            skipped=True
        else:
            skipped=False
            exit_status = process_se_image(executable, config_file, imfile, 
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








def generate_se_filenames(serun, imfile, 
                          rootdir=None, outdir=None, mohrify=True):
    info=deswl.files.get_info_from_path(imfile, 'red')
    fdict={}

    rc=deswl.files.Runconfig()

    # output file names
    for wltype in rc.se_filetypes:
        name= deswl.files.wlse_path(serun, 
                                    info['redrun'], 
                                    info['exposurename'], 
                                    info['ccdnum'], 
                                    wltype,
                                    rootdir=rootdir,
                                    mohrify=mohrify)
        if outdir is not None:
            name=deswl.files.replacedir(name, outdir)
        fdict[wltype] = name

    # stdout, stderr, and log names
    for executable in rc.se_executables:
        key=executable+'_stdout'
        name = deswl.files.wlse_path(serun, 
                                     info['redrun'], 
                                     info['exposurename'], 
                                     info['ccdnum'], 
                                     executable,
                                     rootdir=rootdir,
                                     ext='.stdout')

        if outdir is not None:
            name=deswl.files.replacedir(name, outdir)
        fdict[key] = name

        key=executable+'_stderr'
        name = deswl.files.wlse_path(serun, 
                                     info['redrun'], 
                                     info['exposurename'], 
                                     info['ccdnum'], 
                                     executable,
                                     rootdir=rootdir,
                                     ext='.stderr')
        if outdir is not None:
            name=deswl.files.replacedir(name, outdir)
        fdict[key] = name

        key=executable+'_log'
        name = deswl.files.wlse_path(serun, 
                                     info['redrun'], 
                                     info['exposurename'], 
                                     info['ccdnum'], 
                                     executable+'_log',
                                     rootdir=rootdir,
                                     ext='.xml')
        if outdir is not None:
            name=deswl.files.replacedir(name, outdir)
        fdict[key] = name



    return fdict




def get_info_from_image_path(imfile):
    """
    This will probably also work for cat files
    """

    root = deswl.files.remove_fits_extension(imfile)

    DESDATA=getenv_check('DESDATA')
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

    namefront=deswl.files.remove_fits_extension(odict['basename'])
    #odict['ccdnum'] = int(namefront.split('_')[-1])
    odict['ccdnum'] = int(namefront.split('_')[1])
    
    return odict




def bach_split_imagecat_filelist(input_fname, roottag='bachsplit'):
    """
    We don't need this any more since we have a batch system

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









def _find_wl_output(wlserun, src, typename, mohrify=True):
    fname=deswl.files.wlse_path(wlserun, 
                                src['run'],
                                src['exposurename'],
                                int(src['ccd']),
                                typename,
                                mohrify=mohrify)
    if not os.path.exists(fname):
        stdout.write("WL file not found: %s\n" % fname)
        return False
    else:
        src[typename+'file'] = fname
        return True






def find_collated_coaddfiles(dataset, band, localid, 
                             withsrc=True, 
                             serun=None, mohrify=True,
                             rootdir=None,
                             verbosity=1, dowrite=True):
    """

    Find coadd files on the local disk.  The collated_coaddfiles, created
    from the db, is used as the starting point.

    If withsrc=True then the SE images that when into the coadd are also
    searched for.

    If serun is not None then the corresponding WL outputs for the SE images
    are searched for.
    
    and collate with the input source
    SE images and the weak lensing SE outputs.

    Also if srclist is in the dictionaries these files are also searched for as
    well as the corresponding fitpsf file needed for input into the multishear
    code, and the shear catalog for SE.  For the fitpsf files you must enter
    serun.

    """


    # note newest=False since we haven't discovered the newest yet
    output_path = deswl.files.collated_coaddfiles_path(dataset, band, 
                                                       withsrc=withsrc,
                                                       serun=serun,
                                                       newest=False,
                                                       localid=localid)
    if os.path.exists(output_path) and dowrite:
        raise ValueError("Output file already exists: %s" % output_path)
    stdout.write('Output file: %s\n' % output_path)

    infodict = deswl.files.collated_coaddfiles_read(dataset, band,
                                                    withsrc=withsrc)
    allinfo=infodict['info']

    # where to look
    if rootdir is None:
        rootdir=getenv_check('DESDATA')


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
        
        imfile=deswl.files.coadd_image_path(info['catalogrun'], 
                                            info['tilename'], 
                                            band, 
                                            rootdir=rootdir, 
                                            check=True)
        catfile=deswl.files.coadd_cat_path(info['catalogrun'], 
                                           info['tilename'], 
                                           band, 
                                           rootdir=rootdir, 
                                           check=True)

        srcfound=True
        catfound=True
        imfound=True
        if imfile is None:
            imfound=False
            n_imagemissing += 1

        if catfile is None:
            catfound=False
            n_catmissing += 1

        srclist=[]
        if checksources and imfound and catfound:
            for src in info['srclist']:
                keepsrc = True
                srcfile=deswl.files.red_image_path(src['run'], 
                                                   src['exposurename'], 
                                                   src['ccd'],
                                                   check=True, 
                                                   rootdir=rootdir)
                if srcfile is None:
                    keepsrc=False
                else:
                    src['srcfile'] = srcfile


                # check for the fitted psf file and shear file
                if serun is not None:
                    for ftype in ['fitpsf','shear']:
                        res=_find_wl_output(serun,src,ftype,mohrify=mohrify)
                        if not res:
                            keepsrc=False
                
                if keepsrc:
                    srclist.append(src)

            if len(srclist) < len(info['srclist']):
                stdout.write("Found fewer sources than listed %s/%s\n" \
                             % (len(srclist), len(info['srclist'])))
                n_srcshort += 1
                if len(srclist) == 0:
                    stdout.write("skipping the rest for this coadd %s\n" \
                                 % imfile)
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
    stdout.write("Keeping %s\n" % len(output))

    stdout.write("Output file: %s\n" % output_path)

    if dowrite:
        outdir = os.path.dirname(output_path)
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        outdict={'rootdir':rootdir, 'info': output}
        xmltools.dict2xml(outdict, output_path, roottag='tileinfo_local')
        stdout.write("\n  * Don't forget to check it into SVN!!!\n\n")


def get_newest_coaddfiles(dataset, band, localid, withsrc=True, serun=None, 
                          dowrite=True):
    """
    Take an output file from find_collated_coaddfiles() and return only the
    newest version of each tile (by coaddrun-catalogrun) 
    """


    output_path = deswl.files.collated_coaddfiles_path(dataset, band, 
                                                       withsrc=withsrc,
                                                       serun=serun,
                                                       newest=True,
                                                       localid=localid)

    if os.path.exists(output_path) and dowrite:
        raise ValueError("Output file already exists: %s" % output_path)
    stdout.write('Output file: %s\n' % output_path)


    infodict=deswl.files.collated_coaddfiles_read(dataset, band, 
                                                  localid=localid,
                                                  withsrc=withsrc,
                                                  newest=False,
                                                  serun=serun)

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

    output = {'info': tmpdict.values(), 'rootdir':infodict['rootdir']}

    stdout.write('Finally keeping %s/%s\n' % (len(tmpdict), len(allinfo)))
    stdout.write('Output file: %s\n' % output_path)

    if dowrite:
        xmltools.dict2xml(output, output_path, roottag='tileinfo_local')

def generate_me_srclist(dataset, band, localid, serun, tileinfo):
    """
    Generate an input file for multishear for the input tileinfo dict, which
    should describe a single tile version.  

    This just writes out what is stored in the tile info structure
    (found_coaddinfo).
    """

    tilename=tileinfo['tilename']
    band=tileinfo['band']
    fpath=deswl.files.wlme_srclist_path(dataset,localid,tilename,band,serun)
    outdir=os.path.dirname(fpath)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    stdout.write('Writing multishear input file: %s\n' % fpath)


    fobj=open(fpath, 'w')
    for s in tileinfo['srclist']:
        line = '%s %s %s\n' % (s['srcfile'],s['shearfile'],s['fitpsffile'])
        fobj.write(line)

    fobj.close()

def generate_me_srclists(dataset, band, localid, serun):
    """
    Generate multiepoch srclists for the given dataset and band.
    """

    infodict, tileinfo_file=\
            deswl.files.collated_coaddfiles_read(dataset, 
                                                 band, 
                                                 localid=localid,
                                                 serun=serun,
                                                 getname=True)
    tileinfo=infodict['info']


    for ti in tileinfo:
        generate_me_srclist(dataset, band, localid, serun, ti)

def _make_setup_command(prodname, vers):
    """
    convention is that trunk is exported to ~/exports/prodname-work
    """
    setup_command = "setup %s" % prodname
    if vers == 'trunk' and prodname != 'desfiles':
        setup_command += ' -r ~/exports/%s-work' % prodname
    else:
        setup_command += ' %s' % vers

    return setup_command

def generate_me_pbsfile(merun, tilename, band, outfile,
                        nodes=1, ppn=8, walltime='2:00:00'):

    # the job name
    jobname=tilename+'-'+band

    # The useful log file we redirect from the actual script call
    logf=os.path.basename(outfile).replace('.pbs','.log')
    logf=path_join('${HOME}', 'desprocess', 'pbs', merun, logf)

    # the generally useless pbslog
    pbslog_file = os.path.basename(logf.replace('.log','.pbslog'))

    scratch_rootdir=default_scratch_rootdir()

    header="""#!/bin/bash
#PBS -S /bin/bash
#PBS -N %s
#PBS -j oe
#PBS -o %s
#PBS -m a
#PBS -V
#PBS -r n
#PBS -W umask=0022
#PBS -l walltime=%s
#PBS -l nodes=%s:ppn=%s

# explanation of above directives
# see also: http://wiki.hpc.ufl.edu/index.php/PBS_Directives
#
# -N name - the name of the job
# -j oe   - merge stdout and stderr in the output log (see -o)
# -o name - put the output here.  Not usually useful for logging that
#                   needs to be watched since it isn't actually written till the
#                       script finishes.
# -m a    - send mail when there is an abort
# -V      - All environment variables in the qsub command's environment
#                       are to be exported to the batch job.
# -r n    - Don't try to rerun the job
# -W umask=0022 - you can specify lots of things with -W. In this case make
#                 sure the umask is correct.  I often put this in the 
#                 script itself too
# -l walltime=24:00:00 - "-l" is used to specify resources.  walltime gives
#                        the max wall time the script can take before it
#                        gets killed.
#
# some stuff not in this header
#
# -l nodes=number_of_nodes:ppn=process_per_node
#    Note, in modern parlance, ppn is really number of cores to use.
# -q queue  - Run in this queue

umask 0022

# The main log file, updated as the script runs
logf="%s"

# make sure our log directory exists
d=$(dirname $logf)
if [ ! -e "$d" ]; then
    mkdir -p "$d"
fi

# get eups ready
source /global/data/products/eups/bin/setups.sh
""" % (jobname, pbslog_file, walltime, nodes, ppn, logf)

    rc=deswl.files.Runconfig(merun)
    esutil_setup = _make_setup_command('esutil', rc['esutilvers'])
    tmv_setup = _make_setup_command('tmv',rc['tmvvers'])
    wl_setup = _make_setup_command('wl',rc['wlvers'])

    
    fobj=open(outfile, 'w')

    fobj.write(header)
    fobj.write(tmv_setup+'\n')
    fobj.write(wl_setup+'\n')
    # must come after if it is not the current version
    fobj.write(esutil_setup+'\n')

    fobj.write("multishear-run    \\\n")
    fobj.write('    --merun=%s    \\\n' % merun)
    fobj.write('    --rootdir=%s  \\\n' % scratch_rootdir)
    fobj.write('    --copyroot    \\\n')
    fobj.write('    --cleanup     \\\n')
    fobj.write('    --redirect    \\\n')
    fobj.write('    %s %s &> "$logf"\n' % (tilename, band))

    fobj.write('\n')
    fobj.close()

def generate_me_pbsfiles(merun, dataset, localid, band, serun, outdir=None): 

    if outdir is None:
        outdir=path_join('~','desprocess','pbs',merun)

    outdir=os.path.expanduser(outdir)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    infodict,tileinfo_file = \
        deswl.files.collated_coaddfiles_read(dataset, band, 
                                             localid=localid,
                                             serun=serun,getname=True)
    tileinfo = infodict['info']

    stdout.write("Creating pbs files\n")

    submit_fname=os.path.join(outdir, 'submit.sh')
    submit=open(submit_fname, 'w')

    for ti in tileinfo:
        tilename=ti['tilename']
        band=ti['band']
        name=deswl.files.wlme_pbs_name(tilename,band,merun)

        outfile = os.path.join(outdir, name)
        stdout.write('Writing pbs file: %s\n' % outfile)
        generate_me_pbsfile(merun, tilename, band, outfile)

        submit.write('qsub %s\n' % os.path.basename(outfile))

    
    submit.close()

    readme_fname=os.path.join(outdir, 'README')
    stdout.write('Creating %s\n' % readme_fname)
    readme=open(readme_fname, 'w')
    mess="""
    Input file list: 
        %s
    Which should be a unique list of tiles.  When duplicate tiles are present
    the latest based on coaddrun-catalogrun was chosen.
    \n""" % tileinfo_file

    readme.write(mess)
    readme.close()




def make_me_commandlist(executable, 
                        config_file, 
                        coadd_srclist, 
                        coaddimage_file, 
                        coaddcat_file, 
                        multishear_file,
                        shear_dc4_input_format=True,
                        merun=None,
                        debug=0):

    # eventually change multishear_sky_method to MAP when we get that data
    command = [executable,
               config_file,
               'coadd_srclist='+coadd_srclist,
               'coaddimage_file='+coaddimage_file,
               'coaddcat_file='+coaddcat_file,
               'multishear_file='+multishear_file,
               'multishear_sky_method=NEAREST']
    if shear_dc4_input_format:
        command.append('shear_dc4_input_format=true')

    if merun is not None:
        command.append('merun=%s' % merun)


    if debug:
        debug_file=multishear_file.replace('.fits','.debug')
        command += ['verbose=%s' % debug, 'debug_file='+debug_file]
    return command

def process_me_tile(executable, 
                    config_file, 
                    coadd_srclist,
                    coaddimage_file, 
                    coaddcat_file,
                    multishear_file,
                    stdout_file=None,
                    stderr_file=None,
                    merun=None,
                    timeout=60*60,
                    debug=0):
    """
    All paths must be sent explicitly.  Use the wrapper run_multishear to
    specify a smaller subset
    """

    command = make_me_commandlist(executable, 
                                  config_file, 
                                  coadd_srclist, 
                                  coaddimage_file, 
                                  coaddcat_file, 
                                  multishear_file,
                                  merun=merun,
                                  debug=debug)

    exit_status, stdout_ret, stderr_ret = \
        execute_command(command, timeout=timeout, 
                         stdout_file=stdout_file, 
                         stderr_file=stderr_file,
                         verbose=True)

    stdout.write('exit_status: %s\n' % exit_status)
    return exit_status


def generate_me_filenames(tilename, band, merun=None, dir=None, rootdir=None):
    if merun is not None:
        multishear_file=deswl.files.wlme_path(merun, tilename, band, 
                                              dir=dir, rootdir=rootdir)
        stdout_file=deswl.files.wlme_path(merun, tilename, band, 
                                          ext='.stdout',
                                          dir=dir, rootdir=rootdir)
        stderr_file=deswl.files.wlme_path(merun, tilename, band, 
                                          ext='.stderr',
                                          dir=dir, rootdir=rootdir)
        log_file=deswl.files.wlme_path(merun, tilename, band, 
                                       extra='log', ext='.xml',
                                       dir=dir, rootdir=rootdir)
    else:
        if dir is None:
            dir='.'

        multishear_file=deswl.files.wlme_basename(tilename, band)
        stdout_file=deswl.files.wlme_basename(tilename, band, ext='.stdout')
        stderr_file=deswl.files.wlme_basename(tilename, band, ext='.stderr')
        log_file=deswl.files.wlme_basename(tilename, band, 
                                           extra='log', ext='.xml')

    named={'multishear':multishear_file,
           'stdout':stdout_file,
           'stderr':stderr_file,
           'log':log_file}
    return named

def run_multishear(tilename, band, 
                   outdir=None, 
                   rootdir=None,
                   copyroot=False,
                   cleanup=False,
                   merun=None,
                   dataset=None, 
                   localid=None,
                   serun=None,
                   srclist=None,
                   redirect=False,
                   debug=0):
    """

    deswl.wlpipe.run_multishear(tilename,band
                                outdir=None, 
                                rootdir=None,
                                copyroot=False,
                                cleanup=False,
                                merun=None,
                                dataset=None, 
                                localid=None,
                                serun=None,
                                srclist=None,
                                redirect=True,
                                debug=0)                               

    Wrapper function requiring minimial info to process a tile, namely the
    tilename and the band.  The latest version of the tile is found and used.
    All paths are generated.  The output file is placed in '.' unless
    specified via the optional outdir/rootdir parameters or the merun is sent.

    For this code to work, you need to setup the products WL and DESFILES,
    which sets paths and environment variables such as $WL_DIR.  Note WL
    automaticaly sets up DESFILES.  The information for the tiles is kept in
    the file

        $DESFILES_DIR/(dataset)/
           (dataset)-images-catalogs-withsrc-(localid)-newest-(band).xml

    This contains a list of the newest versions of all the unique tiles.  The
    exececutable file for multishear and the config file are under $WL_DIR/bin
    and $WL_DIR/etc.  Note the variables 'dataset' and 'localid', which default
    to dataset='dc4coadd' and localid='tutti'.  The localid may be relevant since
    different sites may contain different newest versions.

    The input source list is in 
        $DESFILES_DIR/(dataset)/multishear-srclists/
           (tilename)-(band)-srclist.dat

    You can override the srclist file with the srclist= keyword.

    Example using the library:
        deswl.wlpipe.run_multishear('DES2215-2853', 'i')
        
   
    Parameters
        tilename: e.g. DES2215-2853  
        band: e.g. 'i'

    Optional parameters:
        merun: The multi-epoch run identifier.  If sent then the default output 
            directory is $DESDATA/wlbnl/run/me_shapelet/tilename.
            The merun= keyword is sent to multishear and added to the
            header.  An merun implies and overrides dataset and localid.
        outdir: The output directory.  Default '.'.  If merun is sent the
            default directory is the "standard" place.  See below.
        rootdir: The root directory when merun sent.  Default is $DESDATA.
        copyroot: Copy the output files to the default path under $DESDATA.
            Useful when writing to local scratch disk and you want to copy
            the data over to shared disk upon completion.
        cleanup: If copyroot=True and cleanup=True, then the local copyis
            deleted.
        dataset:  The dataset to use.  This tells the code how to get the
            list of available tiles (see above).  Default is 
            is used if merun is not sent.  see default_dataset()
        localid:  The localid holding the files.  The available tile list may
            be different on different machines.  Default is used if
            merun is not sent. see default_localid()
        serun:  The serun used to generate the SE shear input for the
            srclist.  Default defined in default_serun()
        srclist: The list of input SE images used to create the coadd.  If
            not sent the default is as listed above.
        redirect: If True redirect the stdout/stderr to files
            with the default names.
        debug: The level of debugging information printed.  Default 0.  
            1 gives "dbg" level and 2 gives "xdbg" level.\n\n"""

    # for making full directory trees

    if merun is not None:
        runconfig=deswl.files.Runconfig(merun)
        # make sure we have consistency in the software versions
        stdout.write('Verifying runconfig: ...')
        runconfig.verify()
        stdout.write('OK\n')
        if outdir is None:
            outdir = deswl.files.wlme_dir(merun, tilename, rootdir=rootdir)

        dataset=runconfig['dataset']
        localid=runconfig['localid']
        serun=runconfig['serun']
    else:
        runconfig=None
        if outdir is None:
            outdir='.'

        if dataset is None:
            dataset = default_dataset()
        if localid is None:
            localid = default_localid()
        if serun is None:
            serun = default_serun()



    tm1=time.time()
    wl_dir=getenv_check('WL_DIR')
    desfiles_dir=getenv_check('DESFILES_DIR')

    pyvers=deswl.get_python_version()
    esutilvers=esutil.version()
    wlvers=deswl.version()
    tmvvers=deswl.get_tmv_version()


    # the executable
    executable=path_join(wl_dir, 'bin','multishear')

    # config file
    config=path_join(wl_dir, 'etc','wl.config')

    # source list can be determined from tilename/band
    if srclist is None:
        srclist=\
            deswl.files.wlme_srclist_path(dataset,localid,tilename,band,serun)

    # info on all the unique, newest tiles and catalogs
    stdout.write('\n')
    infodict, tileinfo_file=\
        deswl.files.collated_coaddfiles_read(dataset, band, 
                                             localid=localid,
                                             serun=serun,
                                             getname=True)
    tileinfo=infodict['info']
    stdout.write('\n')

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

    fdict=generate_me_filenames(tilename, band, 
                                merun=merun, dir=outdir, rootdir=rootdir)
    multishear_file=fdict['multishear']
    if redirect:
        stdout_file=fdict['stdout']
        stderr_file=fdict['stderr']
    else:
        stdout_file=None
        stderr_file=None
    log_file=fdict['log']


    log={}
    if merun is not None:
        log['merun'] = merun
        #log['runconfig'] = runconfig.asdict()
    else:
        log['merun'] = 'unspecified'
        #log['runconfig'] = 'unspecified'

    log['dataset'] = dataset
    log['localid'] = localid
    log['serun'] = serun
    log['pyvers'] = pyvers
    log['esutilvers'] = esutilvers
    log['wlvers'] = wlvers
    log['tmvvers'] = tmvvers
    log['WL_DIR'] = wl_dir
    log['DESFILES_DIR'] = desfiles_dir
    log['executable'] = executable
    log['tileinfo_file'] = tileinfo_file
    log['tilename'] = tilename
    log['band'] = band
    log['srclist'] = srclist
    log['outdir'] = outdir
    if rootdir is None:
        rootdir=getenv_check('DESDATA')
    log['rootdir'] = rootdir
    log['copyroot'] = copyroot
    log['cleanup'] = cleanup
    log['debug'] = debug

    log['coaddimage'] = coaddimage
    log['coaddcat'] = coaddcat
    log['multishear_file'] = multishear_file
    log['hostname'] = platform.node()
    log['date'] = datetime.datetime.now().strftime("%Y-%m-%d-%X")

    # Print some info
    stdout.write('Configuration:\n')
    if merun is not None:
        stdout.write('    merun: %s\n' % merun)
        for key in runconfig:
            stdout.write("        rc['%s']: %s\n" % (key, runconfig[key]))
    stdout.write('    dataset: %s\n' % dataset)
    stdout.write('    data localid: %s\n' % localid)
    stdout.write('    serun: %s\n' % serun)
    stdout.write('    python version: %s\n' % pyvers)
    stdout.write('    tmv version: %s\n' % tmvvers)
    stdout.write('    wl version: %s\n' % wlvers)
    stdout.write('    WL_DIR: %s\n' % wl_dir)
    stdout.write('    DESFILES_DIR: %s\n' % desfiles_dir)
    stdout.write('    executable: %s\n' % executable)
    stdout.write('    tileinfo file: %s\n' % tileinfo_file)
    stdout.write('    tilename: %s\n' % tilename)
    stdout.write('    band: %s\n' % band)
    stdout.write('    srclist: %s\n' % srclist)
    stdout.write('    outdir: %s\n' % outdir)
    stdout.write('    rootdir: %s\n' % rootdir)
    stdout.write('    copyroot: %s\n' % copyroot)
    stdout.write('    cleanup: %s\n' % cleanup)
    stdout.write('    debug: %s\n' % debug)
    stdout.write('    coaddimage: %s\n' % coaddimage)
    stdout.write('    coaddcat: %s\n' % coaddcat)
    stdout.write('    multishear_file: %s\n' % multishear_file)
    stdout.write('    date: %s\n' % log['date'])
    stdout.write('    running on: %s\n' % log['hostname'])

    # make sure outdir exists
    if not os.path.exists(outdir):
        stdout.write('\nMaking output directory path: %s\n\n' % outdir)
        os.makedirs(outdir)


    exit_status = process_me_tile(executable, 
                                  config, 
                                  srclist, 
                                  coaddimage, 
                                  coaddcat, 
                                  multishear_file, 
                                  stdout_file=stdout_file,
                                  stderr_file=stderr_file,
                                  merun=merun,
                                  debug=debug)

    
    log['exit_status'] = exit_status
    stdout.write('Writing log file: %s\n' % log_file)
    execbase=os.path.basename(executable)
    xmltools.dict2xml(log, log_file, roottag=execbase+'_log')


    # If requested, and if rootdir != the default, then copy the data files
    # into the default.  Also, if cleanup=True, then remove the local copy

    if copyroot:
        fdict_def=generate_me_filenames(tilename, band, merun=merun)
        # see if they point to the same thing.
        if fdict['log'] != fdict_def['log']:
            dirname=os.path.dirname(fdict_def['log'])
            if not os.path.exists(dirname):
                os.makedirs(dirname)
            # copy the files
            for key in ['log','stdout','stderr','multishear']:
                f=fdict[key]
                df=fdict_def[key]
                stdout.write(' * copying %s to %s\n' % (f,df) )
                if not os.path.exists(f):
                    stdout.write('    Source file does not exist\n')
                else:
                    shutil.copy2(f, df)
                    if cleanup:
                        stdout.write(' * deleting %s\n' % (f,) )
                        os.remove(f)




    tm2=time.time()
    ptime(tm2-tm1, format='multishear execution time: %s\n')


def check_multishear(merun, band,
                     rootdir=None, 
                     tilelist=None,
                     badlist_file=None):
    """
    Currently only works by merun identifier
    """

    rc=deswl.files.Runconfig('merun')
    dataset=rc['dataset']
    localid=rc['localid']
    if tilelist is None:
        infodict=deswl.files.collated_coaddfiles_read(dataset, band, localid=localid)
        tileinfo=infodict['info']
    else:
        # fake up a tileinfo structure
        tileinfo=[]
        for tile in tilelist:
            ti={'tilename': tile, 'band': band}
            tileinfo.append(ti)

    badlist=[]
    for ti in tileinfo:
        tilename=ti['tilename']
        fdict=generate_me_filenames(tilename, band, 
                                    merun=merun, rootdir=rootdir)

        
        logfile=fdict['log']
        stdofile=fdict['stdout']
        stdefile=fdict['stderr']
        msfile=fdict['multishear']

        error_found=False
        missing=[]
        if not os.path.exists(logfile):
            missing.append('logfile')
            error_found=True
        if not os.path.exists(stdofile):
            missing.append('stdout')
            error_found=True
        if not os.path.exists(stdefile):
            missing.append('stderr')
            error_found=True
        if not os.path.exists(msfile):
            missing.append('multishear')
            error_found=True

        if len(missing) > 0:
            missmess=','.join(missing)
            stdout.write('%s-%s: %s missing\n' % (tilename,band,missmess))

        if 'logfile' not in missing:
            log=xmltools.xml2dict(logfile,noroot=True)
            exit_status=int( log['exit_status'] )
            if exit_status != 0:
                stdout.write("%s-%s: Found non-zero exit "
                             "status %s in logfile %s\n" % \
                                (tilename,band,exit_status,logfile))
                error_found=True
            
        if error_found:
            badlist.append(ti)

    if badlist_file is not None:
        xmltools.dict2xml({'files':badlist}, badlist_file, roottag='filelist')
    return badlist
        
