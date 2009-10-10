"""
    Order of doing things for SE image processing:


    Find the collated, newest versions of SE "red" images and catalogs
        deswl.wlpipe.find_collated_redfiles(dataset, band)
    dataset is just a label, this command currently will look for all
    such files. Files get put under
        $DESFILES_DIR/(dataset)/(dataset)-images-catalogs-(band).xml



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
                          mohrify=False):

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


def default_dataset(run_type):
    if run_type == 'wlme':
        return 'dc4coadd'
    if run_type == 'wlse':
        return 'dc4'

def default_localid():
    return 'tutti2'
def default_serun():
    return 'wlse000001'
def default_scratch_rootdir():
    return '/scratch/esheldon/DES'


def get_red_filelist(band=None, extra=None, newest=False, return_dict=False):
    """
    Send extra="_cat" to get the catalog list
    """

    command="des-list-files -c red -e ''"
    if band is not None:
        command = command + " -b %s" % band
    if extra is not None:
        command = command + " -e %s" % extra

    exit_status, flist, serr = execute_command(command, verbose=True)

    if exit_status != 0:
        raise RuntimeError,\
          'des-list-files: exit status: %s stderr: %s' % (exit_status, serr)
    
    flist=flist.split('\n')
    while flist[-1] == '':
        flist.pop()

    if newest:
        return select_newest_red_run(flist, return_dict=return_dict)
    else:
        return flist

def select_newest_red_run(image_list, return_dict=False):
    """
    Return the newest versions in file list based on run.  By default
    returns the list, but can also return the info dict keyed by
    exposurename-ccd

    seems to work with cat files too
    """
    kinfo={}
    for imfile in image_list:
        info=deswl.files.get_info_from_path(imfile, 'red')
        key=info['exposurename']+'-'+str(info['ccd'])
        key = '%s-%02i' % (info['exposurename'], info['ccd'])
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
    Find the newest catalogs and images and match them up

    This is *everything on disk*!
    """
    stdout.write('Getting newest image list\n')
    imdict = get_red_filelist(band=band, newest=True, return_dict=True)
    stdout.write('Getting newest cat list\n')
    catdict = get_red_filelist(band=band, newest=True, extra='_cat', 
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

def find_collated_redfiles(dataset, band):
    """
    Get the matched SE image and catalog lists and write to a file

    dataset is just a label, in fact everything from the disk is lumped
    together.  Need a way to distinguish the old from the new.
    """

    flist = get_matched_red_image_catlist(band)

    infolist=[]
    for fpair in flist:
        imfile=fpair['imfile']
        info=deswl.files.get_info_from_path(imfile, 'red') 
        info['imfile'] = imfile
        info['catfile'] = fpair['catfile']
        infolist.append(info)

    deswl.files.collated_redfiles_write(dataset, band, infolist)


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
    if verbose:
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

    #tm0=time.time()

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

    #ptime(time.time()-tm0,format='Execution time: %s\n')

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





def generate_se_pbsfile_byccd(serun, exposurename, ccd, outfile,
                              nodes=1, ppn=8, walltime='1:00:00'):

    """
    """
    # the job name
    jobname=exposurename+'-%02i' % int(ccd)

    # The useful log file we redirect from the actual script call
    logf=outfile+'.log'

    # the generally useless pbslog
    pbslog_file = os.path.basename(outfile.replace('.pbs','.pbslog'))

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

    rc=deswl.files.Runconfig(serun)
    esutil_setup = _make_setup_command('esutil', rc['esutilvers'])
    tmv_setup = _make_setup_command('tmv',rc['tmvvers'])
    wl_setup = _make_setup_command('wl',rc['wlvers'])

    
    fobj=open(outfile, 'w')

    fobj.write(header)
    fobj.write(tmv_setup+'\n')
    fobj.write(wl_setup+'\n')
    # must come after if it is not the current version
    fobj.write(esutil_setup+'\n')

    fobj.write("shear-run    \\\n")
    fobj.write('    --serun=%s    \\\n' % serun)
    fobj.write('    --rootdir=%s  \\\n' % scratch_rootdir)
    fobj.write('    --copyroot    \\\n')
    fobj.write('    --cleanup     \\\n')
    fobj.write('    --nodots      \\\n')
    fobj.write('    --writelog    \\\n')
    fobj.write('    %s %s &> "$logf"\n' % (exposurename, ccd))

    fobj.write('\n')
    fobj.close()


def generate_se_pbsfile(serun, exposurename, outfile,
                        nodes=1, ppn=8, walltime='4:00:00'):

    """
    These jobs can run on a single core/processor.  So we set
    nodes=1 and ppn=1
    """

    # the job name
    jobname=exposurename

    # The useful log file we redirect from the actual script call
    logf=outfile+'.log'

    # the generally useless pbslog
    pbslog_file = os.path.basename(outfile.replace('.pbs','.pbslog'))

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

    rc=deswl.files.Runconfig(serun)
    esutil_setup = _make_setup_command('esutil', rc['esutilvers'])
    tmv_setup = _make_setup_command('tmv',rc['tmvvers'])
    wl_setup = _make_setup_command('wl',rc['wlvers'])

    
    fobj=open(outfile, 'w')

    fobj.write(header)
    fobj.write(tmv_setup+'\n')
    fobj.write(wl_setup+'\n')
    # must come after if it is not the current version
    fobj.write(esutil_setup+'\n')

    fobj.write("shear-run    \\\n")
    fobj.write('    --serun=%s    \\\n' % serun)
    fobj.write('    --rootdir=%s  \\\n' % scratch_rootdir)
    fobj.write('    --copyroot    \\\n')
    fobj.write('    --cleanup     \\\n')
    fobj.write('    --nodots      \\\n')
    fobj.write('    --writelog    \\\n')
    fobj.write('    %s &> "$logf"\n' % (exposurename,))

    fobj.write('\n')
    fobj.close()




def generate_se_pbsfiles(serun, dataset, band, byccd=False): 

    stdout.write("Creating pbs files\n")


    # get the file lists
    infodict,infodict_file = \
        deswl.files.collated_redfiles_read(dataset, band, getpath=True)
    fileinfo = infodict['flist']


    if byccd:
        outdir=deswl.files.pbs_dir(serun, subdir='byccd')
    else:
        outdir=deswl.files.pbs_dir(serun)
    if not os.path.exists(outdir):
        os.makedirs(outdir)


    submit_fname=os.path.join(outdir, 'submit.sh')
    if byccd:
        submit_fname=submit_fname.replace('submit','submit-byccd')

    submit=open(submit_fname, 'w')

    if not byccd:
        # get unique exposure names
        edict={}
        for fi in fileinfo:
            edict[fi['exposurename']] = fi
        
        for exposurename in edict:
            outfile = deswl.files.wlse_pbs_path(serun,exposurename)
            stdout.write('Writing pbs file: %s\n' % outfile)
            generate_se_pbsfile(serun, exposurename, outfile)
            submit.write('qsub %s\n' % os.path.basename(outfile))
    else:
        for fi in fileinfo:
            exposurename=fi['exposurename']
            ccd=fi['ccd']

            outfile = deswl.files.wlse_pbs_path(serun,exposurename,ccd)
            stdout.write('Writing pbs file: %s\n' % outfile)
            generate_se_pbsfile_byccd(serun, exposurename, ccd, outfile)

            submit.write('qsub %s\n' % os.path.basename(outfile))

    
    submit.close()

    readme_fname=os.path.join(outdir, 'README')
    stdout.write('Creating %s\n' % readme_fname)
    readme=open(readme_fname, 'w')
    mess="""
    Input file list: 
        %s
    Which should be a unique list of files.  When duplicateis are present
    the newest based on run ids was chosen.
    \n""" % infodict_file

    readme.write(mess)
    readme.close()







def make_se_commandlist(fdict, debug=0):
    command=[fdict['executable'],
             fdict['config'],
             'image_file='+fdict['image'],
             'cat_file='+fdict['cat'],
             'stars_file='+fdict['stars'],
             'fitpsf_file='+fdict['fitpsf'],
             'psf_file='+fdict['psf'],
             'shear_file='+fdict['shear'] ]

    if fdict['nodots']:
        command.append('output_dots=false')

    if 'serun' in fdict:
        command.append('serun=%s' % fdict['serun'])


    if debug:
        debug_file=fdict['qa'].replace('qa','debug')
        command += ['verbose=%s' % debug, 'debug_file='+debug_file]

    return command

def process_se_image(fdict, writelog=False, debug=0, timeout=5*60):

    # should we re-direct stdout (which is only QA and STATUS messages)
    # to the QA log file?
    if writelog:
        stdout.write('\nWill write QA/STATUS to %s\n\n' % fdict['qa'])
        stdout_file=open(fdict['qa'], 'w')
    else:
        stdout_file=None

    command = make_se_commandlist(fdict, debug=debug)

    # stdout_file=None because we want to allow the caller to redirect that
    exit_status, stdout_ret, stderr_ret = \
        execute_command(command, timeout=timeout,
                        stdout_file=stdout_file,stderr_file=None)

    stdout.write('exit_status: %s\n' % exit_status)
    return exit_status


def getband_from_exposurename(exposurename):
    for band in ['g','r','i','z','Y']:
        pattern='-%s-' % band
        if exposurename.find(pattern) != -1:
            return band
    raise ValueError("No band found in exposurename: %s\n" % exposurename)

def run_shear(exposurename, ccd=None,
              outdir=None, 
              rootdir=None,
              copyroot=False,
              cleanup=False,
              serun=None,
              dataset=None, 
              writelog=False,
              nodots=False,
              debug=0):
    """

    deswl.wlpipe.run_shear(tilename,band
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
    tilename and the band.  The newest version of the tile is found and used.
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
        deswl.wlpipe.run_shear('DES2215-2853', 'i')
        
   
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




    # if ccd not sent, process the whole exposure
    if ccd is None:
        tmall0=time.time()
        for ccd in range(1,62+1):
            run_shear(exposurename, ccd,
                      outdir=outdir, 
                      rootdir=rootdir,
                      copyroot=copyroot,
                      cleanup=cleanup,
                      serun=serun,
                      dataset=dataset, 
                      writelog=writelog,
                      nodots=nodots,
                      debug=debug)
        ptime(time.time() - tmall0, format='Time for all 62 ccds: %s\n')
        return



    tm1=time.time()


    band=getband_from_exposurename(exposurename)

    if serun is not None:
        runconfig=deswl.files.Runconfig(serun)
        # make sure we have consistency in the software versions
        stdout.write('Verifying runconfig: ...')
        runconfig.verify()
        stdout.write('OK\n')
        if outdir is None:
            outdir = deswl.files.wlse_dir(serun, exposurename, rootdir=rootdir)

        dataset=runconfig['dataset']
    else:
        runconfig=None
        if outdir is None:
            outdir='.'

        if dataset is None:
            dataset = default_dataset('wlse')



    wl_dir=getenv_check('WL_DIR')
    desfiles_dir=getenv_check('DESFILES_DIR')

    pyvers=deswl.get_python_version()
    esutilvers=esutil.version()
    wlvers=deswl.version()
    tmvvers=deswl.get_tmv_version()


    # the executable
    executable=path_join(wl_dir, 'bin','fullpipe')

    # config file
    config=path_join(wl_dir, 'etc','wl.config')

    # info on all the unique, newest tiles and catalogs
    stdout.write('\n')
    infodict, infodict_file=\
        deswl.files.collated_redfiles_read(dataset, band, getpath=True)

    infolist=infodict['flist']
    stdout.write('\n')

    # Here we rely on the fact that we picked out the unique, newest
    # version for each tile
    image=None
    bname = exposurename+'_%02i' % int(ccd)
    for ti in infolist:
        tbname=deswl.files.remove_fits_extension(ti['basename'])
        if tbname == bname:
            image=ti['imfile']
            cat=ti['catfile']
            break


    if image is None:
        raise ValueError("Exposure ccd '%s' not found in '%s'" % \
                         (bname, infodict_file))

    fdict=deswl.files.generate_se_filenames(exposurename,ccd,
                                            serun=serun, 
                                            dir=outdir, 
                                            rootdir=rootdir)

    fdict['config'] = config
    fdict['executable'] = executable
    fdict['image'] = image
    fdict['cat'] = cat
    fdict['nodots'] = nodots
    if serun is not None:
        fdict['serun'] = serun



    stat={}
    if serun is not None:
        stat['serun'] = serun
    else:
        stat['serun'] = 'unspecified'

    stat['dataset'] = dataset
    stat['pyvers'] = pyvers
    stat['esutilvers'] = esutilvers
    stat['wlvers'] = wlvers
    stat['tmvvers'] = tmvvers
    stat['WL_DIR'] = wl_dir
    stat['DESFILES_DIR'] = desfiles_dir
    stat['infodict_file'] = infodict_file
    stat['exposurename'] = exposurename
    stat['ccd'] = ccd
    stat['outdir'] = outdir
    if rootdir is None:
        rootdir=getenv_check('DESDATA')
    stat['rootdir'] = rootdir
    stat['copyroot'] = copyroot
    stat['cleanup'] = cleanup
    stat['debug'] = debug

    stat['imfile'] = image
    stat['catfile'] = cat
    stat['hostname'] = platform.node()
    stat['date'] = datetime.datetime.now().strftime("%Y-%m-%d-%X")
    for f in fdict:
        stat[f] = fdict[f]

    hostshort=stat['hostname'].split('.')[0]

    # Print some info
    stdout.write('Configuration:\n')
    if serun is not None:
        stdout.write('    serun: %s\n' % serun)
        for key in runconfig:
            stdout.write("        rc['%s']: %s\n" % (key, runconfig[key]))
    stdout.write('    dataset: %s\n' % dataset)
    stdout.write('    python version: %s\n' % pyvers)
    stdout.write('    tmv version: %s\n' % tmvvers)
    stdout.write('    wl version: %s\n' % wlvers)
    stdout.write('    WL_DIR: %s\n' % wl_dir)
    stdout.write('    DESFILES_DIR: %s\n' % desfiles_dir)
    stdout.write('    infodict_file: %s\n' % infodict_file)
    stdout.write('    exposurename: %s\n' % exposurename)
    stdout.write('    ccd: %s\n' % ccd)
    stdout.write('    outdir: %s\n' % outdir)
    stdout.write('    rootdir: %s\n' % rootdir)
    stdout.write('    copyroot: %s\n' % copyroot)
    stdout.write('    cleanup: %s\n' % cleanup)
    stdout.write('    debug: %s\n' % debug)
    for f in fdict:
        stdout.write("        fdict['%s']: %s\n" % (f,fdict[f]))
    stdout.write('    date: %s\n' % stat['date'])
    stdout.write('    running on: %s\n' % stat['hostname'])

    # make sure outdir exists
    if not os.path.exists(outdir):
        stdout.write('\nMaking output directory path: %s\n\n' % outdir)
        os.makedirs(outdir)


    # need to write one of these for se
    exit_status = process_se_image(fdict, writelog=writelog, debug=debug)

    
    stat['exit_status'] = exit_status
    stdout.write('\nWriting stat file: %s\n' % fdict['stat'])
    execbase=os.path.basename(executable)
    xmltools.dict2xml(stat, fdict['stat'], roottag=execbase+'_stats')


    # If requested, and if rootdir != the default, then copy the data files
    # into the default.  Also, if cleanup=True, then remove the local copy

    if copyroot:
        fdict_def=\
            deswl.files.generate_se_filenames(exposurename, ccd, serun=serun)
        # see if they point to the same thing.
        if fdict['stat'] != fdict_def['stat']:
            dirname=os.path.dirname(fdict_def['stat'])
            if not os.path.exists(dirname):
                os.makedirs(dirname)
            # copy the files
            rc=deswl.files.Runconfig()
            for ftype in rc.se_filetypes:
                f=fdict[ftype]
                df=fdict_def[ftype]
                stdout.write(' * copying %s:%s to %s\n' % (hostshort,f,df) )
                if not os.path.exists(f):
                    stdout.write('    Source file does not exist\n')
                else:
                    shutil.copy2(f, df)
                    if cleanup:
                        stdout.write(' * deleting %s:%s\n' % (hostshort,f,) )
                        os.remove(f)




    tm2=time.time()
    ptime(tm2-tm1, format='fullpipe execution time: %s\n')










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









def _find_wl_output(wlserun, src, typename, mohrify=False):
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
                             serun=None, mohrify=False,
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
                                                 getpath=True)
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
    logf=outfile+'.log'

    # the generally useless pbslog
    pbslog_file = os.path.basename(outfile.replace('.pbs','.pbslog'))

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
    fobj.write('    --nodots      \\\n')
    fobj.write('    --writelog    \\\n')
    fobj.write('    %s %s &> "$logf"\n' % (tilename, band))

    fobj.write('\n')
    fobj.close()

def generate_me_pbsfiles(merun, dataset, localid, band, serun): 

    stdout.write("Creating pbs files\n")


    # get the file lists
    infodict,tileinfo_file = \
        deswl.files.collated_coaddfiles_read(dataset, band, 
                                             localid=localid,
                                             serun=serun,getpath=True)
    tileinfo = infodict['info']



    outdir=deswl.files.pbs_dir(merun)
    if not os.path.exists(outdir):
        os.makedirs(outdir)


    submit_fname=os.path.join(outdir, 'submit.sh')
    submit=open(submit_fname, 'w')



    for ti in tileinfo:
        tilename=ti['tilename']
        band=ti['band']

        outfile = deswl.files.wlme_pbs_path(merun,tilename,band)
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
    the newest based on coaddrun-catalogrun was chosen.
    \n""" % tileinfo_file

    readme.write(mess)
    readme.close()




def make_me_commandlist(fdict,shear_dc4_input_format=True, debug=0):

    # eventually change multishear_sky_method to MAP when we get that data
    command=[fdict['executable'],
             fdict['config'],
             'coadd_srclist='+fdict['srclist'],
             'coaddimage_file='+fdict['coaddimage'],
             'coaddcat_file='+fdict['coaddcat'],
             'multishear_file='+fdict['multishear'],
             'multishear_sky_method=NEAREST']


    if fdict['nodots']:
        command.append('output_dots=false')

    if shear_dc4_input_format:
        command.append('shear_dc4_input_format=true')

    if 'merun' in fdict:
        command.append('merun=%s' % fdict['merun'])

    if debug:
        debug_file=fdict['multishear'].replace('.fits','.debug')
        command += ['verbose=%s' % debug, 'debug_file='+debug_file]
    return command

def process_me_tile(fdict, 
                    writelog=False,
                    timeout=2*60*60,  # two hours
                    debug=0):
    """
    All paths must be sent explicitly.  Use the wrapper run_multishear to
    specify a smaller subset
    """

    # should we re-direct stdout (which is only QA and STATUS messages)
    # to the QA log file?
    if writelog:
        stdout.write('\nWill write QA/STATUS to %s\n\n' % fdict['qa'])
        stdout_file=open(fdict['qa'], 'w')
    else:
        stdout_file=None


    command = make_me_commandlist(fdict,debug=debug)

    exit_status, stdout_ret, stderr_ret = \
        execute_command(command, timeout=timeout,
                        stdout_file=stdout_file,stderr_file=None)

    stdout.write('exit_status: %s\n' % exit_status)
    return exit_status



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
                   writelog=False,
                   nodots=False,
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
                                writelog=False,
                                debug=0)                               

    Wrapper function requiring minimial info to process a tile, namely the
    tilename and the band.  The newest version of the tile is found and used.
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
            with the default names. (currently ignored)
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
            dataset = default_dataset('wlme')
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
                                             getpath=True)
    tileinfo=infodict['info']
    stdout.write('\n')

    # Here we rely on the fact that we picked out the unique, newest
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

    fdict=deswl.files.generate_me_filenames(tilename, band, 
                                            merun=merun, dir=outdir, 
                                            rootdir=rootdir)
    multishear_file=fdict['multishear']


    fdict['config'] = config
    fdict['executable'] = executable
    fdict['coaddimage'] = coaddimage
    fdict['coaddcat'] = coaddcat
    fdict['srclist'] = srclist
    fdict['nodots'] = nodots
    if merun is not None:
        fdict['merun'] = merun


    stat={}
    if merun is not None:
        stat['merun'] = merun
        #stat['runconfig'] = runconfig.asdict()
    else:
        stat['merun'] = 'unspecified'
        #stat['runconfig'] = 'unspecified'

    stat['dataset'] = dataset
    stat['localid'] = localid
    stat['serun'] = serun
    stat['pyvers'] = pyvers
    stat['esutilvers'] = esutilvers
    stat['wlvers'] = wlvers
    stat['tmvvers'] = tmvvers
    stat['WL_DIR'] = wl_dir
    stat['DESFILES_DIR'] = desfiles_dir
    stat['executable'] = executable
    stat['tileinfo_file'] = tileinfo_file
    stat['tilename'] = tilename
    stat['band'] = band
    stat['srclist'] = srclist
    stat['outdir'] = outdir
    if rootdir is None:
        rootdir=getenv_check('DESDATA')
    stat['rootdir'] = rootdir
    stat['copyroot'] = copyroot
    stat['cleanup'] = cleanup
    stat['debug'] = debug

    stat['coaddimage'] = coaddimage
    stat['coaddcat'] = coaddcat
    stat['multishear_file'] = multishear_file
    stat['hostname'] = platform.node()
    stat['date'] = datetime.datetime.now().strftime("%Y-%m-%d-%X")

    hostshort=stat['hostname'].split('.')[0]

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
    stdout.write('    outdir: %s\n' % outdir)
    stdout.write('    rootdir: %s\n' % rootdir)
    stdout.write('    copyroot: %s\n' % copyroot)
    stdout.write('    cleanup: %s\n' % cleanup)
    stdout.write('    debug: %s\n' % debug)
    stdout.write('    date: %s\n' % stat['date'])
    stdout.write('    running on: %s\n' % stat['hostname'])
    for f in fdict:
        stdout.write("        fdict['%s']: %s\n" % (f,fdict[f]))

    # make sure outdir exists
    if not os.path.exists(outdir):
        stdout.write('\nMaking output directory path: %s\n\n' % outdir)
        os.makedirs(outdir)


    exit_status = process_me_tile(fdict, writelog=writelog, debug=debug)
    
    stat['exit_status'] = exit_status
    stdout.write('Writing stat file: %s\n' % fdict['stat'])
    execbase=os.path.basename(executable)
    xmltools.dict2xml(stat, fdict['stat'], roottag=execbase+'_stats')


    # If requested, and if rootdir != the default, then copy the data files
    # into the default.  Also, if cleanup=True, then remove the local copy

    if copyroot:
        fdict_def=\
            deswl.files.generate_me_filenames(tilename, band, merun=merun)
        # see if they point to the same thing.
        if fdict['stat'] != fdict_def['stat']:
            dirname=os.path.dirname(fdict_def['stat'])
            if not os.path.exists(dirname):
                os.makedirs(dirname)
            # copy the files
            rc=deswl.files.Runconfig()
            for ftype in rc.me_filetypes:
                f=fdict[ftype]
                df=fdict_def[ftype]
                stdout.write(' * copying %s:%s to %s\n' % (hostshort,f,df) )
                if not os.path.exists(f):
                    stdout.write('    Source file does not exist\n')
                else:
                    shutil.copy2(f, df)
                    if cleanup:
                        stdout.write(' * deleting %s:%s\n' % (hostshort,f,) )
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
        
