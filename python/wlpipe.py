"""
    vim: set ft=python:
    Order of doing things for SE image processing

        I now use the database to figure out what SE files are available
        for a given *release*.  The release dr012 corresponds to DC6b.

        Use ~/python/desdb/bin/get-release-runs.py to get a list of runs
        and then use wget-des/wget-des-parallel to download them.

        Generate md5sums using ~/shell_scripts/checksum-files.  Then verify against 
          https://desar.cosmology.illinois.edu/DESFiles/desardata/Release/{release}/MD5SUM
        Here release is all caps.  Use 
            ~/python/des/bin/checksum-compare.py 
            or 
            ~/python/des/src/checksum-compare (compile it)
        for the comparison.

        Then, use ~/python/desdb/bin/get-red-info.py to produce a .json file
        with the image and cat info, including path relative to $DESDATA You
        should output this to $DESFILES/{release} with a name like
        {release}-red-{band}-info.json

        # create a run name and its configuration
            rc=deswl.files.Runconfig()
            rc.generate_new_runconfig('wlse', dataset, wl_config)

        Note I've started setting wl_config here, instead of using
        what is under etc. This allows us to use the same version
        of code but run with a different configuration.
        I keep these in $DESFILES_DIR/wl.config, e.g. wl01.config
        
        You can send test=True and dryrun=True also

        # create the condor scripts
        sec=deswl.wlpipe.SECondorJobs(wlserun,band,type='fullpipe', nthread=1)
        sec.write()
        OR
        sec.write_byccd()

        OLDER PBS STUFF
            # create pbs scripts, either by pointing or by ccd number
            # these files execute the shear-run stand alone code

            generate_se_pbsfiles(serun, band, byccd=False): 

        # after running you can check the results.  This will write
        # out lists of successful and unsuccessful jobs

        check_shear(serun, band)

        # the files written are these.  Need to make it bandpass 
        # dependent
        deswl.files.wlse_collated_path(serun, 'goodlist')
        deswl.files.wlse_collated_path(serun, 'badlist')

    For coadds, using multishear:
        On a machine that can directly access the desdm database:

            import deswl
            deswl.oracle.collate_coadd_catim(band,collated_file,
                                             getsrc=True)


        For DC4 I'm putting all this stuff into the nyu repo 'desfiles':
            $DESFILES_DIR/dc4coadd
        Generally
          $DESFILES_DIR/(dataset)/(dataset)-images-catalogs-(band).json
          $DESFILES_DIR/(dataset)/(dataset)-images-catalogs-withsrc-(band).json

        This data is not really versioned.  The product should be installed and
        set up so that DESFILES_DIR is set (this is also done by setup wl).

        To process data on your local machine, you need to find the files on
        local disk:

            deswl.wlpipe.find_collated_coaddfiles(dataset, band, localid, 
                                                  withsrc=True,
                                                  serun=None,
                                                  rootdir=None,
                                                  verbosity=1,
                                                  dowrite=True)
        most common usage will be
            find_collated_coaddfiles(dataset, band, localid, serun=something)
        because this will be the input to the coadd stuff later.

        I'm keeping these here:
            $DESFILES_DIR/(dataset)/(localid)/
               (dataset)-images-catalogs-withsrc-(localid)-(band).json

        And to get the newest versions, send that file to:
            deswl.wlpipe.get_newest_coaddfiles(infofile, output_file)


        Stored in
            $DESFILES_DIR/(dataset)/(localid)/
              (dataset)-images-catalogs-withsrc-(localid)-newest-(band).json


        Then this file is used to generate source lists (SE images that 
        contribute to the coadd tiles) for input to multishear.

            deswl.wlpipe.generate_me_srclists(dataset,band,localid,serun)

        Which calls generate_me_srclist.  This will put files under the
        given output directory.  Results are usually checked into desfiles 
        under
            $DESFILES_DIR/(dataset)/(localid)/multishear-srclist/
                          (tilename)-(band)-srclist.dat





        create a run name and its configuration
            rc=deswl.files.Runconfig()
            rc.generate_new_runconfig('wlme',dataset,localid=,serrun=)

        You can send test=True and dryrun=True also

        


        CHECK THE ABOVE FILES INTO SVN!!!


        To write out the pbs files needed to process the new run:
            deswl.wlpipe.generate_me_pbsfiles(merun,band)


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
import traceback
import platform
import os
import glob
import subprocess
import time
import copy
import signal
import re
import datetime
import shutil
import pprint
import logging

import deswl

import esutil
from esutil.ostools import path_join, getenv_check
from esutil import json_util
from esutil.misc import ptime

if deswl.get_python_version(numerical=True) >= 3:
    inputfunc = input
else:
    inputfunc = raw_input


def default_dataset(run_type):
    if run_type == 'wlme':
        return 'dc4coadd'
    if run_type == 'wlse':
        return 'dc4new'

def default_localid():
    return 'tutti2'
def default_serun():
    return 'wlse000001'
def default_scratch_rootdir():
    return '/scratch/users/esheldon/DES'

ERROR_SE_UNKNOWN=2**0
ERROR_SE_MISC=2**1
ERROR_SE_MISSING_FILE=2**2
ERROR_SE_LOAD_CONFIG=2**3
ERROR_SE_LOAD_FILES=2**4
ERROR_SE_FINDSTARS=2**5
ERROR_SE_MEASURE_PSF=2**6
ERROR_SE_MEASURE_SHEAR=2**7
ERROR_SE_SPLIT_STARS=2**8
ERROR_SE_IO=2**9
ERROR_SE_SET_LOG=2**10


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


def generate_se_checkpsf_pbsfile(serun, exposurename, outfile,
                                 nodes=1, ppn=1, walltime='24:00:00'):

    """
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

# get eups, DESDATA, etc ready
source /home/users/dataman/.bashrc
SCRATCH=/scratch/esheldon/DES

""" % (jobname, pbslog_file, walltime, nodes, ppn, logf)

    rc=deswl.files.Runconfig(serun)
    esutil_setup = _make_setup_command('esutil', rc['esutilvers'])
    tmv_setup = _make_setup_command('tmv',rc['tmvvers'])
    wl_setup = _make_setup_command('wl',rc['wlvers'])

    
    fobj=open(outfile, 'w')

    fobj.write(header)
    fobj.write(wl_setup+'\n')

    # must come after if it is not the current version
    fobj.write(tmv_setup+'\n')
    fobj.write(esutil_setup+'\n')

    script="""
config=${WL_DIR}/etc/wl.config

serun=%s
exposurename=%s
input_dir=${DESDATA}/wlbnl/${serun}/${exposurename}
outdir=${SCRATCH}/wlbnl/${serun}/${exposurename}

if [ ! -e $outdir ]; then
    mkdir -p $outdir
fi

ccds=$(seq -w 1 62)

for ccd in $ccds; do

    froot=${serun}_${exposurename}_${ccd}

    psf_file=${input_dir}/${froot}_psf.fits
    fitpsf_file=${input_dir}/${froot}_fitpsf.fits
    out_file=${outdir}/${froot}_checkpsf.rec

    test-psfrec $config $psf_file $fitpsf_file > $out_file
    echo "output file: $out_file"
done
    """ % (serun, exposurename)

    fobj.write(script)

    fobj.close()

class SECondorJobs(dict):
    def __init__(self, wlserun, band, 
                 type='fullpipe',
                 nodes=1, ppn=1, nthread=1):
        """
        Need to implement nodes, ppn
        """
        self['run'] = wlserun
        self['band'] = band
        self['type'] = type
        self['nodes'] = nodes
        self['ppn'] = ppn
        self['nthread'] = nthread

        self.rc = deswl.files.Runconfig(self['run'])
        self['dataset'] = self.rc['dataset']

        # get the file lists
        self.fileinfo = deswl.files.collated_redfiles_read(self['dataset'], band)


    def write(self, dryrun=False):
        # get unique exposure names
        edict={}
        for fi in self.fileinfo:
            edict[fi['exposurename']] = fi
        
        for exposurename in edict:
            sejob = SECondorJob(self['run'], exposurename, 
                                type=self['type'],
                                nodes=self['nodes'], 
                                ppn=self['ppn'], 
                                nthread=self['nthread'])
            sejob.write_submit(dryrun=dryrun)
            sejob.write_script(dryrun=dryrun)

    def write_byccd(self, dryrun=False):
        for fi in self.fileinfo:
            exposurename=fi['exposurename']
            ccd=fi['ccd']

            sejob = SECondorJob(self['run'], exposurename, ccd=ccd,
                                type=self['type'],
                                nodes=self['nodes'], 
                                ppn=self['ppn'], 
                                nthread=self['nthread'])
            sejob.write_submit(dryrun=dryrun)
            sejob.write_script(dryrun=dryrun)



 
class SECondorJob(dict):
    def __init__(self, wlserun, exposurename, 
                 ccd=None, type='fullpipe',
                 nodes=1, ppn=1, nthread=1):
        """
        Need to implement nodes, ppn
        """
        self['run'] = wlserun
        self['exposurename'] = exposurename
        self['ccd'] = ccd
        self['type'] = type
        self['nodes'] = nodes
        self['ppn'] = ppn
        self['nthread'] = nthread

    def write_submit(self, verbose=False, dryrun=False):
        self.make_condor_dir()
        f=deswl.files.wlse_condor_path(self['run'], 
                                       self['exposurename'], 
                                       typ=self['type'], 
                                       ccd=self['ccd'])

        print "writing to condor submit file:",f 
        text = self.submit_text()
        if verbose or dryrun:
            print text 
        if not dryrun:
            with open(f,'w') as fobj:
                fobj.write(text)
        else:
            print "this is just a dry run" 

    def write_script(self, verbose=False, dryrun=False):
        self.make_condor_dir()
        f=deswl.files.wlse_script_path(self['run'], 
                                       self['exposurename'], 
                                       typ=self['type'], 
                                       ccd=self['ccd'])

        print "writing to condor script file:",f 
        text = self.script_text()
        if verbose or dryrun:
            print text 
        if not dryrun:
            with open(f,'w') as fobj:
                fobj.write(text)

            print "changing mode of file to executable"
            os.popen('chmod 755 '+f)
        else:
            print "this is just a dry run" 


    def submit_text(self):
        condor_dir = self.condor_dir()
        script_base = deswl.files.wlse_script_base(self['exposurename'],
                                                   ccd=self['ccd'],
                                                   typ=self['type'])
        submit_text="""
Universe        = vanilla
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

        return submit_text

    def script_text(self):
        rc=deswl.files.Runconfig(self['run'])
        wl_setup = _make_setup_command('wl',rc['wlvers'])
        tmv_setup = _make_setup_command('tmv',rc['tmvvers'])
        esutil_setup = _make_setup_command('esutil', rc['esutilvers'])
     
        if self['ccd'] == None:
            ccd=''
        else:
            ccd=self['ccd']

        script_text="""#!/bin/bash
source ~astrodat/setup/setup.sh
source ~astrodat/setup/setup-modules.sh
module load use.own
module unload tmv
module load tmv/work
module unload wl
module load wl/work
module unload esutil
module load esutil/work
module load desfiles

export OMP_NUM_THREADS={nthread}
shear-run                     \\
     --serun={serun}          \\
     --nodots                 \\
     {expname} {ccd} 2>&1
""".format(wl_setup=wl_setup, 
           tmv_setup=tmv_setup, 
           esutil_setup=esutil_setup,
           nthread=self['nthread'],
           serun=self['run'],
           expname=self['exposurename'],
           ccd=ccd)

        return script_text

    def condor_dir(self):
        if self['ccd'] is None:
            return deswl.files.condor_dir(self['run'])
        else:
            return deswl.files.condor_dir(self['run'], subdir='byccd')


    def make_condor_dir(self):
        condor_dir = self.condor_dir()
        if not os.path.exists(condor_dir):
            print 'making output dir',condor_dir
            os.makedirs(condor_dir)






def generate_se_pbsfile(serun, exposurename, outfile, ccd=None,
                        nodes=1, ppn=1, walltime=None, queue='fast', 
                        nthread=1):

    """
    """

    # the job name
    if ccd is None:
        jobname=exposurename
        ccd = ''
    else:
        jobname=exposurename+'-%02i' % int(ccd)
    jobname = jobname.replace('decam-','')

    # The useful log file we redirect from the actual script call
    logf=outfile+'.log'

    # the generally useless pbslog
    pbslog_file = os.path.basename(outfile.replace('.pbs','.pbslog'))

    scratch_rootdir=default_scratch_rootdir()

    if walltime is not None:
        walltime='#PBS -l walltime=%s' % walltime
    else:
        walltime=''

    header="""#!/bin/bash
#PBS -S /bin/bash
#PBS -N {jobname}
#PBS -j oe
#PBS -o {pbslog_file}
#PBS -m a
#PBS -V
#PBS -r n
#PBS -W umask=0022
{walltime}
#PBS -l nodes={nodes}:ppn={ppn}
#PBS -q {queue}

umask 0022

# The main log file, updated as the script runs
logf="{logf}"

# get eups ready
source /global/data/products/eups/bin/setups.sh
""".format(jobname=jobname, 
           pbslog_file=pbslog_file, 
           walltime=walltime, 
           nodes=nodes, 
           ppn=ppn, 
           logf=logf,
           queue=queue)

    rc=deswl.files.Runconfig(serun)
    wl_setup = _make_setup_command('wl',rc['wlvers'])
    tmv_setup = _make_setup_command('tmv',rc['tmvvers'])
    esutil_setup = _make_setup_command('esutil', rc['esutilvers'])
    
    fobj=open(outfile, 'w')

    fobj.write(header)

    fobj.write('%s\n' % wl_setup)
    fobj.write('%s\n' % tmv_setup)
    fobj.write('%s\n' % esutil_setup)
    fobj.write('export OMP_NUM_THREADS=%d\n' % nthread)

    fobj.write("\nshear-run       \\\n")
    fobj.write('    --serun=%s    \\\n' % serun)
    #fobj.write('    --rootdir=%s  \\\n' % scratch_rootdir)
    #fobj.write('    --copyroot    \\\n')
    fobj.write('    --nodots      \\\n')
    fobj.write('    %s %s &> "$logf"\n' % (exposurename,ccd))

    fobj.write('\n')
    fobj.close()




def generate_se_pbsfiles(serun, band, typ='fullpipe', byccd=False, nthread=1): 

    stdout.write("Creating pbs files\n")

    rc = deswl.files.Runconfig(serun)
    dataset = rc['dataset']

    # get the file lists
    fileinfo,info_file = \
        deswl.files.collated_redfiles_read(dataset, band, getpath=True)


    if byccd:
        outdir=deswl.files.pbs_dir(serun, subdir='byccd')
    else:
        outdir=deswl.files.pbs_dir(serun)
    if not os.path.exists(outdir):
        os.makedirs(outdir)


    if typ == 'fullpipe':
        submit_fname = 'submit.sh'
    else:
        submit_fname = 'submit-%s.sh' % typ
    submit_fname=os.path.join(outdir, submit_fname)

    if byccd:
        submit_fname=submit_fname.replace('submit','submit-byccd')

    submit=open(submit_fname, 'w')

    if not byccd:
        # get unique exposure names
        edict={}
        for fi in fileinfo:
            edict[fi['exposurename']] = fi
        
        for exposurename in edict:
            outfile = deswl.files.wlse_pbs_path(serun,exposurename,typ=typ)
            stdout.write('Writing pbs file: %s\n' % outfile)
            if typ == 'checkpsf':
                generate_se_checkpsf_pbsfile(serun, exposurename, outfile)
            else:
                generate_se_pbsfile(serun, exposurename, outfile, nthread=nthread)

            bfile = os.path.basename(outfile)
            submit.write('echo -n "%s "; qsub %s\n' % (bfile,bfile))
    else:
        for fi in fileinfo:
            exposurename=fi['exposurename']
            ccd=fi['ccd']

            outfile = deswl.files.wlse_pbs_path(serun,exposurename,
                                                ccd=ccd,typ=typ)
            stdout.write('Writing pbs file: %s\n' % outfile)
            generate_se_pbsfile(serun, exposurename, outfile, ccd=ccd)

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
    \n""" % info_file

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
        command.append('wlserun=%s' % fdict['serun'])


    if debug:
        debug_file=fdict['qa'].replace('qa','debug')
        command += ['verbose=%s' % debug, 'debug_file='+debug_file]

    return command


def run_se_pipe(types, fdict, writelog=False, debug=0, dosplit=True):
    import deswl

    config_fname = fdict['config']

    image = fdict['image']
    cat = fdict['cat']

    # usually is saved as $DESFILES_DIR/..
    config_fname = os.path.expandvars(config_fname)

    # usually the names are saves as $DESDATA/...
    image = os.path.expandvars(image)
    cat = os.path.expandvars(cat)

    if not os.path.exists(config_fname):
        return ERROR_SE_MISSING_FILE
    if not os.path.exists(image):
        return ERROR_SE_MISSING_FILE
    if not os.path.exists(cat):
        return ERROR_SE_MISSING_FILE
    fdict['image']=image
    fdict['cat'] = cat

    try:
        stdout.write("Loading config file '%s'\n" % config_fname)
        t=deswl.WL(config_fname)
    except RuntimeError as e:
        stdout.write('%s\n' % e)
        return ERROR_SE_LOAD_CONFIG

    if writelog:
        # setting log_file is really the QA file.
        stdout.write('\nWill write QA/STATUS to %s\n\n' % fdict['qa'])
        #t.set_param("log_file",fdict['qa'])
        t.set_log(fdict['qa'])

    if debug:
        # setting log_file is really the QA file.
        stdout.write('\nWill write debug to %s\n\n' % fdict['debug'])
        t.set_param('debug',str(debug))
        t.set_param("debug_file",fdict['debug'])

    try:

        stdout.write("loading images from: '%s'\n" % fdict['image'])
        t.load_images(fdict['image'])
        stdout.write("loading catalog from '%s'\n" % fdict['cat'])
        t.load_catalog(fdict['cat'])
    except RuntimeError as e:
        stdout.write('%s\n' % e)
        return ERROR_SE_LOAD_FILES

    if 'stars' in types:
        try:
            stdout.write("finding stars; writing to '%s'\n" % fdict['stars'])
            t.find_stars(fdict['stars'])
        except RuntimeError as e:
            stdout.write('%s\n' % e)
            return ERROR_SE_FINDSTARS

    if 'psf' in types:
        try:
            stdout.write("measuring psf; writing to \n"
                         "    '%s'\n"
                         "    '%s'\n" % (fdict['psf'],fdict['fitpsf']))
            t.measure_psf(fdict['psf'],fdict['fitpsf'])
        except RuntimeError as e:
            stdout.write('%s\n' % e)
            return ERROR_SE_MEASURE_PSF

    if 'shear' in types:
        try:
            stdout.write("measuring shear; writing to '%s'\n" % fdict['shear'])
            t.measure_shear(fdict['shear'])
        except RuntimeError as e:
            stdout.write('%s\n' % e)
            return ERROR_SE_MEASURE_SHEAR


    if 'split' in types:
        try:
            stdout.write("splitting catalog: \n"
                         "    '%s'\n"
                         "    '%s'\n" % (fdict['stars1'],fdict['stars2']))
            t.split_starcat(fdict["stars1"],fdict["stars2"])

            stdout.write('loading stars1\n')
            t.load_starcat(fdict['stars1'])
            stdout.write('measuring psf1\n')
            t.measure_psf(fdict['psf1'],fdict['fitpsf1'])
            stdout.write('measuring shear1\n')
            t.measure_shear(fdict['shear1'])

            stdout.write('loading stars2\n')
            t.load_starcat(fdict['stars2'])
            stdout.write('measuring psf2\n')
            t.measure_psf(fdict['psf2'],fdict['fitpsf2'])
            stdout.write('measuring shear2\n')
            t.measure_shear(fdict['shear2'])

        except RuntimeError as e:
            stdout.write('%s\n' % e)
            return ERROR_SE_SPLIT_STARS

    return 0



def process_se_image_old(fdict, writelog=False, debug=0, timeout=5*60):

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

def process_se_image(config, image, cat, **args):
    ip = ImageProcessor(config, image, cat, **args)
    types = args.get('types','all')
    ip.process(types)

class ImageProcessor(deswl.WL):
    """
    Much simpler than ExposureProcessor for processing a single image/catalog pair.
    Good for interactive use.
    """
    def __init__(self, config, image, cat, **args):
        deswl.WL.__init__(self, config)

        self.config_url = config
        self.image_url = image
        self.cat_url = cat

        outdir = args.get('outdir',None)
        if outdir is None:
            outdir = os.path.dirname(self.cat_url)
        self.outdir = outdir
        if self.outdir == '':
            self.outdir = '.'

        self.set_output_names()

        self.imcat_loaded = False

    def set_output_names(self):
        out_base = os.path.basename(self.cat_url)
        out_base = out_base.replace('_cat.fits','')
        out_base = out_base.replace('.fits','')
        out_base = out_base.replace('.fit','')

        out_base = path_join(self.outdir, out_base)

        rc = deswl.files.Runconfig()
        self.names = {}
        for type in rc.se_filetypes:
            self.names[type] = out_base+'-'+type+'.fits'

    def load_data(self):
        stdout.write('Loading %s\n' % self.image_url)
        self.load_images(self.image_url)
        stdout.write('Loading %s\n' % self.cat_url)
        self.load_catalog(self.cat_url)
        self.imcat_loaded = True

    def process(self, types=['stars','psf','shear','split']):
        if not self.imcat_loaded:
            self.load_data()
        if 'stars' in types:
            self.find_stars()
        if 'psf' in types:
            self.measure_psf()
        if 'shear' in types:
            self.measure_shear()
        if 'split' in types:
            self.split_starcat()
            # note only 'psf' and 'shear' will actually get processed
            self.process_split(types=types)


    def process_full(self):
        if not self.imcat_loaded:
            self.load_data()
        self.find_stars()
        self.measure_psf()
        self.measure_shear()

    def find_stars(self):
        print("Finding stars.  Will write to: %s" % self.names['stars'])
        deswl.WL.find_stars(self, self.names['stars'])
    def measure_psf(self):
        print("Measuring PSF. Will write to: "
              "\n    %s\n    %s" % (self.names['psf'],self.names['fitpsf']))
        deswl.WL.measure_psf(self,self.names['psf'],self.names['fitpsf'])
    def measure_shear(self):
        print("Measuring shear.  Will write to: %s" % self.names['shear'])
        deswl.WL.measure_shear(self,self.names['shear'])

    def split_starcat(self):
        stdout.write("splitting catalog: \n"
                     "    '%s'\n"
                     "    '%s'\n" % (self.names['stars1'],self.names['stars2']))
        deswl.WL.split_starcat(self, self.names["stars1"],self.names["stars2"])

        """
        stdout.write('loading stars1\n')
        self.load_starcat(fdict['stars1'])
        stdout.write('measuring psf1\n')
        deswl.WL.measure_psf(self, self.names['psf1'], self.names['fitpsf1'])
        stdout.write('measuring shear1\n')
        deswl.WL.measure_shear(self,self.names['shear1'])

        stdout.write('loading stars2\n')
        self.load_starcat(fdict['stars2'])
        stdout.write('measuring psf2\n')
        deswl.WL.measure_psf(self, self.names['psf2'], self.names['fitpsf2'])
        stdout.write('measuring shear2\n')
        deswl.WL.measure_shear(self,self.names['shear2'])
        """

    def process_split(self, splitnum=[1,2], types=['psf','shear']):
        for num in splitnum:
            nstr = str(num)

            t = ImageProcessor(self.config_url, self.image_url, self.cat_url, 
                               outdir=self.outdir) 

            t.names['psf'] = t.names['psf'+nstr]
            t.names['fitpsf'] = t.names['fitpsf'+nstr]
            t.names['shear'] = t.names['shear'+nstr]

            t.load_starcat(t.names['stars'+nstr])

            # need to do this explicitly just in case 'stars' or 'split' 
            # were accidentally sent
            if 'psf' in types:
                t.process('psf')
            if 'shear' in types:
                t.process('shear')



class ExposureProcessor:
    def __init__(self, exposurename, **keys):
        """
        We don't load anything in the constructor because we
        want to have state available and always write the
        status file.
        """
        self.all_types = ['stars','psf','shear','split']

        # things in stat will get written to the status QA
        # file
        self.stat = {}
        
        self.stat['exposurename'] = exposurename
        self.stat['serun']     = keys.get('serun',None)
        self.stat['config']    = keys.get('config',None)
        self.stat['outdir']    = keys.get('outdir',None)
        self.stat['copyroot']  = keys.get('copyroot',False)
        self.stat['nodots']    = keys.get('nodots',False)

        DESDATA = deswl.files.des_rootdir()
        self.stat['DESDATA'] = DESDATA
        self.stat['rootdir']   = keys.get('rootdir',DESDATA)

        self.stat['error']        = 0
        self.stat['error_string'] = ''

        self.stat['hostname'] = platform.node()
        self.stat['date'] = datetime.datetime.now().strftime("%Y-%m-%d-%X")

        self.runconfig = None

        self.logger = logging.getLogger('ExposureProcessor')
        log_level = keys.get('log_level',logging.ERROR)
        self.logger.setLevel(log_level)

        self.verbose=0


    def set_config(self):
        # config file
        defconfig = path_join(self.stat['WL_DIR'], 'etc','wl.config')

        if self.stat['config'] is None:
            stdout.write("Determining config...\n")
            if self.runconfig is not None:
                if 'wl_config' in self.runconfig:
                    self.stat['config'] = self.runconfig['wl_config']
                else:
                    self.stat['config'] = defconfig
            else:
                self.stat['config'] = defconfig
            stdout.write("    using: '%s'\n" % self.stat['config'])
        else:
            stdout.write("Using input config: '%s'\n" % self.stat['config'])

    def load_serun(self):

        if self.stat['serun'] is not None:
            stdout.write("    Loading runconfig for serun: '%s'\n" % self.stat['serun'])
            stdout.flush()
            try:
                self.runconfig=deswl.files.Runconfig(self.stat['serun'])
            except RuntimeError as e:
                self.logger.debug("caught RuntimeError load_serun loading Runconfig")
                self.stat['error'] = ERROR_SE_IO
                self.stat['error_string'] = "Error loading run config: '%s'" % e
                stdout.write(self.stat['error_string']+'\n')
                raise e
            except IOError as e:
                self.stat['error'] = ERROR_SE_IO
                self.stat['error_string'] = "Error loading run config: '%s'" % e
                raise e

                
            # make sure we have consistency in the software versions
            stdout.write('    Verifying runconfig: ...'); stdout.flush()
            try:
                self.runconfig.verify()
            except ValueError as e:
                self.logger.debug("caught RuntimeError in load_serun verifying runconfig")
                self.stat['error'] = ERROR_SE_MISC
                self.stat['error_string'] = "Error verifying runconfig: '%s'" % e
                stdout.write(self.stat['error_string']+'\n')
                raise e

            stdout.write('OK\n'); stdout.flush()

            self.stat['dataset'] = self.runconfig['dataset']
        else:
            self.runconfig=None

    def set_outdir(self):
        """
        This method is unused
        """
        stdout.write("Setting outdir:\n")
        if self.stat['outdir'] is None:
            if self.stat['serun'] is not None:
                self.stat['outdir'] = \
                        deswl.files.wlse_dir(self.stat['serun'], 
                                             self.stat['exposurename'], 
                                             rootdir=self.stat['rootdir'])
            else:
                self.stat['outdir']='.'

        self.stat['outdir'] = os.path.expandvars(self.stat['outdir'])
        if self.stat['outdir'] == '.':
            self.stat['outdir']=os.path.abspath('.')

        stdout.write("    outdir: '%s'\n" % self.stat['outdir'])


    def get_environ(self):
        try:
            self.stat['TMV_DIR']=getenv_check('TMV_DIR')
            self.stat['WL_DIR']=getenv_check('WL_DIR')
            self.stat['ESUTIL_DIR']=getenv_check('ESUTIL_DIR')
            self.stat['DESFILES_DIR']=getenv_check('DESFILES_DIR')
            self.stat['pyvers']=deswl.get_python_version()
            self.stat['esutilvers']=esutil.version()
            self.stat['wlvers']=deswl.version()
            self.stat['tmvvers']=deswl.get_tmv_version()

        except RuntimeError as e:
            self.logger.debug("caught RuntimeError in get_environ")
            self.stat['error'] = ERROR_SE_MISC
            self.stat['error_string'] = unicode(str(e),errors='replace')
            stdout.write(self.stat['error_string']+'\n')
            raise e

    def load_wl(self):
        # usually is saved as $DESFILES_DIR/..
        config_fname = os.path.expandvars(self.stat['config'])
        stdout.write("Loading config file '%s'\n" % config_fname)
        stdout.flush()

        if not os.path.exists(config_fname):
            self.stat['error'] = ERROR_SE_MISSING_FILE
            self.stat['error_string'] = \
                "Config File does not exist: '%s'" % config_fname
            raise IOError(self.stat['error_string'])
        
        try:
            self.wl = deswl.WL(str(config_fname))
            if self.stat['nodots']:
                self.wl.set_param("output_dots","false");
            if self.verbose > 0:
                self.wl.set_verbose(self.verbose)
            
            # this will go in the header
            if self.stat['serun']  is not None:
                self.wl.set_param("wlserun",self.stat['serun'])

        except RuntimeError as e:
            self.stat['error'] = ERROR_SE_LOAD_CONFIG
            self.stat['error_string'] = 'Error loading config: %s' % e
            stdout.write(self.stat['error_string']+'\n')
            raise e

    def set_verbose(self, verbose):
        self.verbose=int(verbose)

    def load_file_lists(self):
        if not hasattr(self,'infolist'):
            try:
                info, ifile=\
                        deswl.files.collated_redfiles_read(self.stat['dataset'], 
                                                           self.stat['band'], 
                                                           getpath=True)
                self.infolist = info
                self.stat['info_url'] = ifile
            except IOError as e:
                self.stat['error'] = ERROR_SE_IO
                self.stat['error_string'] = 'Error loading file lists: %s' % e
                stdout.write(self.stat['error_string']+'\n')
                raise e

    def get_band(self):
        stdout.write("    Getting band for '%s'..." % self.stat['exposurename'])
        stdout.flush()
        try:
            self.stat['band']=\
                getband_from_exposurename(self.stat['exposurename'])
            stdout.write(" '%s'\n" % self.stat['band'])
        except RuntimeError as e:
            stdout.write('\n')
            self.logger.debug("caught RuntimeError in get_band")
            self.stat['error'] = ERROR_SE_MISC
            self.stat['error_string'] = unicode(str(e),errors='replace')
            stdout.write(self.stat['error_string']+'\n')
            raise e

    def setup(self):
        """
        Don't call this, let the processing code call it!
        """

        self.logger.debug("In setup")

        # this stuff only needs to be loaded once.
        if not hasattr(self, 'tmv_dir'):

            stdout.write('running setup\n'); stdout.flush()

            self.get_band()
            self.get_environ()
            self.stat['dataset']=default_dataset('wlse')

            # only loads if not None.  Will override the runconfig
            # and dataset
            self.load_serun()

            #self.set_outdir()
            self.set_config()

            self.load_wl()

            self.load_file_lists()

    def get_image_cat(self, ccd):
        # Here we rely on the fact that we picked out the unique, newest
        # version for each

        self.logger.debug("Entering get_image_cat")

        image=None
        expname=self.stat['exposurename']
        for ti in self.infolist:
            if ti['exposurename'] == expname and ti['ccd'] == ccd:
                image=ti['image_url']
                cat=ti['cat_url']
                break

        if image is None:
            self.logger.debug("caught error loading image/cat")
            self.stat['error'] = ERROR_SE_MISC
            self.stat['error_string'] = \
                "Exposure ccd '%s' not found in '%s'" % \
                             (bname, self.info_file)
            raise ValueError(self.stat['error_string'])
        image = os.path.expandvars(image)
        cat = os.path.expandvars(cat)
        return image,cat


    def get_output_filenames(self, ccd):

        try:
            fdict=deswl.files.generate_se_filenames(self.stat['exposurename'],
                                                    ccd,
                                                    serun=self.stat['serun'], 
                                                    dir=self.stat['outdir'], 
                                                    rootdir=self.stat['rootdir'])
            for k in fdict:
                fdict[k] = os.path.expandvars(fdict[k])
            self.stat['output_files'] = fdict

            # if outdir was none, we picked a proper dir
            if self.stat['outdir'] is None:
                self.stat['outdir'] = os.path.dirname(fdict['stars'])

        except RuntimeError as e:
            self.stat['error'] = ERROR_SE_MISC
            self.stat['error_string'] = unicode(str(e),errors='replace')
            stdout.write(self.stat['error_string']+'\n')
            # since we cannot write to a file, just convert it
            # to a string
            stdout.write('%s\n' % pprint.pprint(self.stat))
            raise e
        except AttributeError as e:
            self.stat['error'] = ERROR_SE_MISC
            self.stat['error_string'] = unicode(str(e),errors='replace')
            stdout.write(self.stat['error_string']+'\n')
            stdout.write('%s\n' % pprint.pprint(self.stat))
            raise e
        except NameError as e:
            self.stat['error'] = ERROR_SE_MISC
            self.stat['error_string'] = unicode(str(e),errors='replace')
            stdout.write(self.stat['error_string']+'\n')
            stdout.write('%s\n' % pprint.pprint(self.stat))
            raise e
        except:
            self.stat['error'] = ERROR_SE_UNKNOWN
            self.stat['error_string'] = "Unexpected error: '%s'" % traceback.format_exc()
            stdout.write(self.stat['error_string']+'\n')
            stdout.write('%s\n' % self.stat)
            raise RuntimeError(self.stat['error_string'])
    def load_images(self):

        stdout.write("Loading images from '%s'\n" % self.stat['image'])
        stdout.flush()
        try:
            self.wl.load_images(str(self.stat['image']))
        except RuntimeError as e:
            # the internal routines currently throw const char* which
            # swig is converting to runtime
            self.logger.debug("caught error running load_images")
            self.stat['error'] = ERROR_SE_IO
            self.stat['error_string'] = unicode(str(e),errors='replace')
            stdout.write(self.stat['error_string']+'\n')
            raise e
        except TypeError as e:
            self.logger.debug("caught type error running load_images, probably string conversion")
            self.stat['error'] = ERROR_SE_IO
            self.stat['error_string'] = unicode(str(e),errors='replace')
            stdout.write(self.stat['error_string']+'\n')
            raise e


    def load_catalog(self):

        stdout.write("Loading catalog from '%s'\n" % self.stat['cat'])
        stdout.flush()
        try:
            self.wl.load_catalog(str(self.stat['cat']))
        except RuntimeError as e:
            # the internal routines currently throw const char* which
            # swig is converting to runtime
            self.logger.debug("caught error running load_catalog")
            self.stat['error'] = ERROR_SE_IO
            self.stat['error_string'] = unicode(str(e),errors='replace')
            stdout.write(self.stat['error_string']+'\n')
            raise e
        except TypeError as e:
            self.logger.debug("caught type error running load_catalog, probably string conversion")
            self.stat['error'] = ERROR_SE_IO
            self.stat['error_string'] = unicode(str(e),errors='replace')
            stdout.write(self.stat['error_string']+'\n')
            raise e



    def write_status(self):
        statfile = self.stat['output_files']['stat']
        d=os.path.dirname(statfile)
        if not os.path.exists(d):
            os.makedirs(d)

        stdout.write("Writing status file: '%s'\n" % statfile)
        stdout.write("  staterr: %s\n" % self.stat['error'])
        stdout.write("  staterr string: '%s'\n" % self.stat['error_string'])
        json_util.write(self.stat, statfile)


    def run_findstars(self):
        starsfile = self.stat['output_files']['stars']
        try:
            stdout.write("\nFinding stars. Will write to \n    '%s'\n" \
                         % starsfile)
            stdout.flush()
            self.wl.find_stars(starsfile)
            stars_loaded=True
        except RuntimeError as e:
            self.stat['error'] = ERROR_SE_FINDSTARS
            self.stat['error_string'] = unicode(str(e),errors='replace')
            stdout.write(self.stat['error_string']+'\n')
            raise e

    def ensure_stars_loaded(self):
        if not self.stars_loaded:
            starsfile = self.stat['output_files']['stars']
            stdout.write("Loading star cat: '%s'\n" % starsfile)
            self.wl.load_starcat(starsfile)
            self.stars_loaded=True

    def run_measurepsf(self):
        try:
            self.ensure_stars_loaded() 

            files=self.stat['output_files']
            stdout.write("\nMeasuring psf; Will write to: \n"
                         "    '%s'\n"
                         "    '%s'\n" % (files['psf'],files['fitpsf']))
            stdout.flush()
            self.wl.measure_psf(files['psf'],files['fitpsf'])
            self.psf_loaded=True
        except RuntimeError as e:
            self.stat['error'] = ERROR_SE_MEASURE_PSF
            self.stat['error_string'] = unicode(str(e),errors='replace')
            stdout.write(self.stat['error_string']+'\n')
            raise e
    def ensure_psf_loaded(self):
        if not self.psf_loaded:
            files=self.stat['output_files']
            stdout.write("Loading psf cat: '%s'\n" % files['psf'])
            self.wl.load_psfcat(files['psf'])
            stdout.write("Loading fitpsf: '%s'\n" % files['fitpsf'])
            self.wl.load_fitpsf(files['fitpsf'])
            self.psf_loaded=True

    def run_shear(self):
        try:
            self.ensure_stars_loaded()
            self.ensure_psf_loaded()

            files=self.stat['output_files']
            stdout.write("\nMeasuring shear. Will write to \n    '%s'\n" \
                         % files['shear'])
            stdout.flush()
            self.wl.measure_shear(files['shear'])
            self.shear_loaded=True
        except RuntimeError as e:
            self.stat['error'] = ERROR_SE_MEASURE_SHEAR
            self.stat['error_string'] = unicode(str(e),errors='replace')
            stdout.write(self.stat['error_string']+'\n')
            raise e

    def run_split(self):
        try:
            files=self.stat['output_files']

            self.ensure_stars_loaded()
            
            stdout.write("\nSplitting catalog: \n"
                         "    '%s'\n"
                         "    '%s'\n" % (files['stars1'],files['stars2']))
            stdout.flush()
            self.wl.split_starcat(files["stars1"],files["stars2"])


            stdout.write("\nLoading stars1 \n'%s'\n" % files['stars1'])
            stdout.flush()
            self.wl.load_starcat(files['stars1'])

            stdout.write("\nMeasuring psf1; Will write to: \n"
                         "    '%s'\n"
                         "    '%s'\n" % (files['psf1'],files['fitpsf1']))
            stdout.flush()
            self.wl.measure_psf(files['psf1'],files['fitpsf1'])

            stdout.write("\nMeasuring shear1. Will write to \n'%s'\n" \
                         % files['shear1'])
            stdout.flush()
            self.wl.measure_shear(files['shear1'])

            stdout.write("\n\nLoading stars2 \n'%s'\n" % files['stars2']); stdout.flush()
            self.wl.load_starcat(files['stars2'])

            stdout.write("measuring psf2; Will write to: \n"
                         "    '%s'\n"
                         "    '%s'\n" % (files['psf2'],files['fitpsf2']))
            stdout.flush()
            self.wl.measure_psf(files['psf2'],files['fitpsf2'])

            stdout.write("Measuring shear2. Will write to \n'%s'\n" \
                         % files['shear2'])
            stdout.flush()
            self.wl.measure_shear(files['shear2'])
        except RuntimeError as e:
            self.stat['error'] = ERROR_SE_SPLIT_STARS
            self.stat['error_string'] = unicode(str(e),errors='replace')
            stdout.write(self.stat['error_string']+'\n')
            raise e


    def set_log(self):
        self.logger.debug("In set_log")
        try:
            qafile=self.stat['output_files']['qa']
            stdout.write("    Setting qa file: '%s'\n" % qafile)
            stdout.flush()
            self.wl.set_log(str(qafile))
        except RuntimeError as e:
            self.logger.debug("caught internal RuntimeError running set_log")
            self.stat['error'] = ERROR_SE_SET_LOG
            self.stat['error_string'] = unicode(str(e),errors='replace')
        except TypeError as e:
            self.logger.debug("caught type error running set_log, probably string conversion")
            self.stat['error'] = ERROR_SE_SET_LOG
            self.stat['error_string'] = unicode(str(e),errors='replace')
            stdout.write(self.stat['error_string']+'\n')
            raise e



    def _process_ccd_types(self, types):

        self.stars_loaded=False
        self.psf_loaded=False
        self.shear_loaded=False

        if 'stars' in types:
            self.run_findstars()

        if 'psf' in types:
            self.run_measurepsf()

        if 'shear' in types:
            self.run_shear() 

        if 'split' in types:
            self.run_split()

    def copy_root(self):
        """

        This method does *not* change the self.stat dictionary.  This is
        important because we don't want things like the 'outdir' to change.
        We *might* want to eventually write a new stat file though.

        """
        # this assumes the default root, to which we will copy
        fdict_def=\
            deswl.files.generate_se_filenames(self.stat['exposurename'], 
                                              self.stat['ccd'], 
                                              serun=self.stat['serun'])


        # see if the file is already in the destination directory.  
        if self.stat['output_files']['stat'] != fdict_def['stat']:

            dirname=os.path.dirname(fdict_def['stat'])
            old_dirname = os.path.dirname(self.stat['output_files']['stat'])

            if not os.path.exists(dirname):
                os.makedirs(dirname)
            # copy the files
            rc=deswl.files.Runconfig()
            for ftype in rc.se_filetypes:
                f=self.stat['output_files'][ftype]
                df=fdict_def[ftype]
                hostshort=self.stat['hostname'].split('.')[0]
                stdout.write(' * moving %s:%s to\n    %s\n' % (hostshort,f,df))
                if not os.path.exists(f):
                    stdout.write('    Source file does not exist\n')
                    stdout.flush()
                else:
                    shutil.copy2(f, df)
                    os.remove(f)



    def process_ccd(self, ccd, **keys):
        t0=time.time()

        self.stat['error']        = 0
        self.stat['error_string'] = ''

        self.stat['types'] = keys.get('types',self.all_types)
        # set file names in self.stat
        # if this fails we can't even write the status file, so let the
        # exceptions fall through
        self.get_output_filenames(ccd)

        self.stat['ccd'] = ccd


        # ok, now we at least have the status file name to write to
        try:
            stdout.write('-'*72 + '\n')
            stdout.write("host: %s\n" % platform.node())
            stdout.write("ccd:  %02d\n" % ccd)

            self.setup()

            image,cat = self.get_image_cat(ccd)

            self.stat['image'] = image
            self.stat['cat'] = cat


            self.make_output_dir()
            self.load_images()
            self.load_catalog()
            self.set_log()
            self._process_ccd_types(self.stat['types'])
        except:
            if self.stat['error'] == 0:
                # this was an unexpected error
                self.stat['error'] = ERROR_SE_UNKNOWN
                self.stat['error_string'] = "Unexpected error: '%s'" % traceback.format_exc()
            else:
                # we caught this exception already and set an
                # error state.  Just proceed and write the status file
                pass
        self.write_status()
        if self.stat['copyroot']:
            self.copy_root()
        ptime(time.time() - t0, format='Time to process ccd: %s\n')

    def make_output_dir(self):
        outdir=self.stat['outdir']
        if not os.path.exists(outdir):
            stdout.write("Creating output dir: '%s'\n" % outdir)
            os.makedirs(outdir)

    def process_all_ccds(self, **keys):
        t0=time.time()
        for ccd in xrange(1,62+1):
            self.process_ccd(ccd, **keys)
            stdout.flush()
        ptime(time.time() - t0, format='Time for all 62 ccds: %s\n')
        stdout.flush()


def run_shear(exposurename, 
              ccd=None,
              types=None,
              config_file=None,
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

    deswl.wlpipe.run_shear(exposurename, 
                           ccd=None,
                           types=None,
                           config_file=None,
                           outdir=None, 
                           rootdir=None,
                           copyroot=False,
                           cleanup=False,
                           serun=None,
                           dataset=None, 
                           writelog=False,
                           nodots=False,
                           debug=0):


    Wrapper function requiring minimial info to process an exposure or ccd,
    namely the exposurename and possile ccd.  The output file is placed in '.'
    unless specified via the optional outdir/rootdir parameters or the serun is
    sent.

    For this code to work, you need to setup the products WL and DESFILES,
    which sets paths and environment variables such as $WL_DIR.  Note WL
    automaticaly sets up DESFILES.  

    Parameters
        exposurename: e.g. decam--25--20-i-12 

    Optional parameters:
        ccd: The ccd number from 1-62.
        types: A list of the types of processing to perform.  Default is
            all: stars,psf,shear,split
        config_file: The path to the config file.  Defaults to 
            $WL_DIR/etc/wl.config
        outdir: The output directory.  Default '.'.  If serun is sent the
            default directory is the "standard" place.  See below.
        rootdir: The root directory when serun sent.  Default is $DESDATA.
        copyroot: Copy the output files to the default path under $DESDATA.
            Useful when writing to local scratch disk and you want to copy
            the data over to shared disk upon completion.
        cleanup: If copyroot=True and cleanup=True, then the local copyis
            deleted.
        dataset:  The dataset to use.  This tells the code how to get the
            list of available tiles (see above).  Default is 
            is used if serun is not sent.  see default_dataset()
        serun:  The single epoch shear run identifier.  If sent, is used to
            generate file names and paths.
        writelog: Write out the QA/Status file.
        nodots: Don't print dots as we process to indicate progress.
        debug: The level of debugging information printed.  Default 0.  
            1 gives "dbg" level and 2 gives "xdbg" level.\n\n"""




    # if ccd not sent, process the whole exposure
    if ccd is None:
        tmall0=time.time()
        for ccd in range(1,62+1):
            run_shear(exposurename, ccd,
                      types=types,
                      config_file=config_file,
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

    if types is None:
        types=['stars','psf','shear','split']

    band=getband_from_exposurename(exposurename)

    if serun is not None:
        runconfig=deswl.files.Runconfig(serun)
        # make sure we have consistency in the software versions
        stdout.write('Verifying runconfig: ...')
        runconfig.verify()
        stdout.write('OK\n')
        stdout.flush()
        if outdir is None:
            outdir=deswl.files.wlse_dir(serun, exposurename, rootdir=rootdir)

        dataset=runconfig['dataset']
    else:
        runconfig=None
        if outdir is None:
            outdir='.'

        if dataset is None:
            dataset = default_dataset('wlse')



    tmv_dir=getenv_check('TMV_DIR')
    wl_dir=getenv_check('WL_DIR')
    esutil_dir=getenv_check('ESUTIL_DIR')
    desfiles_dir=getenv_check('DESFILES_DIR')

    pyvers=deswl.get_python_version()
    esutilvers=esutil.version()
    wlvers=deswl.version()
    tmvvers=deswl.get_tmv_version()


    # the executable
    executable=path_join(wl_dir, 'bin','fullpipe')

    # config file
    if config_file is not None:
        config = config_file
    elif runconfig is not None:
        if 'wl_config' in runconfig:
            config = runconfig['wl_config']
    else:   
        config=path_join(wl_dir, 'etc','wl.config')


    # info on all the unique, newest tiles and catalogs
    stdout.write('\n')
    infolist, info_file=\
        deswl.files.collated_redfiles_read(dataset, band, getpath=True)

    stdout.write('\n')

    # Here we rely on the fact that we picked out the unique, newest
    # version for each
    image=None
    for ti in infolist:
        if exposurename == ti['exposurename'] and ti['ccd'] == ccd:
            image=ti['image_url']
            cat=ti['cat_url']
            break


    if image is None:
        raise ValueError("Exposure ccd '%s' not found in '%s'" % \
                         (bname, info_file))

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
    stat['TMV_DIR'] = tmv_dir
    stat['ESUTIL_DIR'] = esutil_dir
    stat['WL_DIR'] = wl_dir
    stat['DESFILES_DIR'] = desfiles_dir
    stat['info_url'] = info_file
    stat['exposurename'] = exposurename
    stat['ccd'] = ccd
    stat['outdir'] = outdir
    stat['types'] = types
    if rootdir is None:
        rootdir=getenv_check('DESDATA')
    stat['rootdir'] = rootdir
    stat['copyroot'] = copyroot
    stat['cleanup'] = cleanup
    stat['debug'] = debug

    stat['image_url'] = image
    stat['cat_url'] = cat
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
    stdout.write('    TMV_DIR: %s\n' % tmv_dir)
    stdout.write('    ESUTIL_DIR: %s\n' % esutil_dir)
    stdout.write('    WL_DIR: %s\n' % wl_dir)
    stdout.write('    DESFILES_DIR: %s\n' % desfiles_dir)
    stdout.write('    info_file: %s\n' % info_file)
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
    #exit_status = process_se_image_old(fdict, writelog=writelog, debug=debug)
    exit_status = run_se_pipe(types, fdict, writelog=writelog, debug=debug)

    
    stat['exit_status'] = exit_status
    stdout.write('\nWriting stat file: %s\n' % fdict['stat'])
    execbase=os.path.basename(executable)
    json_util.write(stat, fdict['stat'])


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
                stdout.write(' * copying %s:%s to\n    %s\n' % (hostshort,f,df))
                if not os.path.exists(f):
                    stdout.write('    Source file does not exist\n')
                else:
                    shutil.copy2(f, df)
                    if cleanup:
                        stdout.write(' * deleting %s:%s\n' % (hostshort,f,) )
                        os.remove(f)




    tm2=time.time()
    ptime(tm2-tm1, format='fullpipe execution time: %s\n')



def check_shear(serun, band, rootdir=None, outdir=None):
    """
    rootdir= and outdir= refer to the directories where the shear outputs
    are located, not outputs of this program
    """
    import copy

    rc=deswl.files.Runconfig(serun)
    dataset = rc['dataset']

    infolist=deswl.files.collated_redfiles_read(dataset, band)

    badlist=[]
    goodlist=[]
    for info in infolist:

        problem_found=False

        exposurename=info['exposurename']
        ccd=info['ccd']

        fdict=deswl.files.generate_se_filenames(exposurename,ccd,
                                                serun=serun, 
                                                dir=outdir, 
                                                rootdir=rootdir)
        info['outputs'] = fdict
        stars_file  = fdict['stars'],
        fitpsf_file = fdict['fitpsf']
        psf_file    = fdict['psf']
        shear_file  = fdict['shear']

        status_read = False
        if not os.path.exists(fdict['stat']):
            print 'Stat file missing:',fdict['stat']
            info = copy.deepcopy(fdict)
            info['error'] = 9999
            info['error_string'] = 'stat file missing'
            info['failtype'] = 'stat file missing'
            badlist.append(info)
        else:
            try:
                stat=json_util.read(fdict['stat'])
                status_read = True
            except ValueError as e:
                print "Error reading file:",fdict['stat']
                info = copy.deepcopy(fdict)
                info['error'] = 9998
                info['error_string'] = 'Error reading stat file'
                info['failtype'] = 'Error reading stat file'
                badlist.append(info)

        if status_read:
            error_code = int( stat['error'] )
            estring = stat['error_string']
            info['error'] = error_code
            info['error_string'] = estring

            if error_code != 0:
                stdout.write("'%s': error: %s: '%s'\n" % (fdict['stat'],error_code,estring))
                info['failtype'] = 'Processing error'
                problem_found=True
                badlist.append(info)

        """
        for ftype in ['stars','fitpsf','psf','shear','qa','stat']:
            fpath=fdict[ftype]
            if not os.path.exists(fpath):
                stdout.write("%s file missing: %s\n" % (ftype,fpath))
                info['failtype'] = ftype
                problem_found=True
                badlist.append(info)
                # break out of file type loop
            elif ftype == 'stat':
                try:
                    stat=json_util.read(fpath)
                except ValueError as e:
                    stdout.write("Error reading file: %s\n" % fpath)
                    raise ValueError(e)

                error_code = int( stat['error'] )
                estring = stat['error_string']
                info['error'] = error_code
                info['error_string'] = estring

                if error_code != 0:
                    stdout.write("'%s': error: %s: '%s'\n" % (fpath,error_code,estring))
                    info['failtype'] = 'error'
                    problem_found=True
                    badlist.append(info)
        """

        if not problem_found:
            goodlist.append(info)

    stdout.write('Found %s/%s problems\n' % (len(badlist),len(infolist)))

    goodfile=deswl.files.wlse_collated_path(serun, 'goodlist')
    badfile=deswl.files.wlse_collated_path(serun, 'badlist')
    dir=deswl.files.wlse_collated_dir(serun)
    if not os.path.exists(dir):
        stdout.write("Creating output dir: %s\n" % dir)
        os.makedirs(dir)
    stdout.write("Writing goodlist: %s\n" % goodfile)
    json_util.write(goodlist, goodfile)
    stdout.write("Writing badlist: %s\n" % badfile)
    json_util.write(badlist, badfile)
    

def check_shear_qa(badlist):
    for b in badlist:
        qa=b['outputs']['qa']
        data = open(qa,'r').read()

        stdout.write('from: %s\n' % qa)
        stdout.write(data)
        stdout.write('-'*50 + '\n')
        key=inputfunc('hit a key ')
        if key.strip().lower() == 'q': 
            return

def check_shear_input_images(badlist, ext=1):
    import pyfits
    for b in badlist:
        imfile=b['image_url']
        stdout.write('\nReading %s[%s]\n' % (imfile,ext))
        try:
            data=pyfits.getdata(imfile,ext=ext)
        except:
            stdout.write("Error: %s\n" % repr(sys.exc_info()))

        key=inputfunc('hit a key ')
        if key.strip().lower() == 'q': 
            return










def _find_wl_output(wlserun, src, typename, mohrify=False):
    fname=deswl.files.wlse_path(src['exposurename'],
                                int(src['ccd']),
                                typename,
                                serun=wlserun,
                                mohrify=mohrify)
    if not os.path.exists(fname):
        stdout.write("WL file not found: %s\n" % fname)
        return False
    else:
        src[typename+'_path'] = fname
        return True




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

    infodict, tileinfo_url=\
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
                        nodes=1, ppn=8, walltime=None,
                        queue='fast'):

    # the job name
    jobname=tilename+'-'+band

    # The useful log file we redirect from the actual script call
    logf=outfile+'.log'

    # the generally useless pbslog
    pbslog_file = os.path.basename(outfile.replace('.pbs','.pbslog'))

    scratch_rootdir=default_scratch_rootdir()

    if walltime is not None:
        walltime='#PBS -l walltime=%s' % walltime
    else:
        walltime=''

    header="""#!/bin/bash
#PBS -S /bin/bash
#PBS -N {jobname}
#PBS -j oe
#PBS -o {pbslog_file}
#PBS -m a
#PBS -V
#PBS -r n
#PBS -W umask=0022
{walltime}
#PBS -l nodes={nodes}:ppn={ppn}
#PBS -q {queue}

umask 0022

# The main log file, updated as the script runs
logf="{logf}"

# get eups ready
source /global/data/products/eups/bin/setups.sh
""".format(jobname=jobname, 
           pbslog_file=pbslog_file, 
           walltime=walltime, 
           nodes=nodes, 
           ppn=ppn, 
           logf=logf,
           queue=queue)
    rc=deswl.files.Runconfig(merun)

    # tmv,esutil must come after wl if wl is not current
    wl_setup = _make_setup_command('wl',rc['wlvers'])
    tmv_setup = _make_setup_command('tmv',rc['tmvvers'])
    esutil_setup = _make_setup_command('esutil', rc['esutilvers'])
    
    fobj=open(outfile, 'w')

    fobj.write(header)
    fobj.write('%s\n' % wl_setup)
    fobj.write('%s\n' % tmv_setup)
    fobj.write('%s\n' % esutil_setup)


    fobj.write("\nmultishear-run    \\\n")
    fobj.write('    --merun=%s    \\\n' % merun)
    fobj.write('    --rootdir=%s  \\\n' % scratch_rootdir)
    fobj.write('    --copyroot    \\\n')
    fobj.write('    --cleanup     \\\n')
    #fobj.write('    --nodots      \\\n')
    fobj.write('    --writelog    \\\n')
    fobj.write('    %s %s &> "$logf"\n' % (tilename, band))

    fobj.write('\n')
    fobj.close()

def generate_me_pbsfiles(merun, band): 
    """
    generate_me_pbsfiles(merun, band)
    """

    stdout.write("Creating pbs files\n")

    rc = deswl.files.Runconfig(merun)
    dataset = rc['dataset']
    localid = rc['localid']
    serun   = rc['serun']

    # get the file lists
    infodict,tileinfo_url = \
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
    \n""" % tileinfo_url

    readme.write(mess)
    readme.close()




def make_me_commandlist(fdict,shear_dc4_input_format=False, debug=0):

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
        command.append('wlmerun=%s' % fdict['merun'])

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
                        stdout_file=stdout_file,stderr_file=None,
                        verbose=True)

    stdout.write('exit_status: %s\n' % exit_status)
    return exit_status



def run_multishear(tilename, band, 
                   config_file=None,
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
                                config_file=None,
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
           (dataset)-images-catalogs-withsrc-(localid)-newest-(band).json

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
        config_file: The path to the config file.  Defaults to 
            $WL_DIR/etc/wl.config
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
    tmv_dir=getenv_check('TMV_DIR')
    wl_dir=getenv_check('WL_DIR')
    esutil_dir=getenv_check('ESUTIL_DIR')
    desfiles_dir=getenv_check('DESFILES_DIR')

    pyvers=deswl.get_python_version()
    esutilvers=esutil.version()
    wlvers=deswl.version()
    tmvvers=deswl.get_tmv_version()


    # the executable
    executable=path_join(wl_dir, 'bin','multishear')

    # config file
    if config_file is not None:
        config = config_file
    elif runconfig is not None:
        if 'wl_config' in runconfig:
            config = runconfig['wl_config']
    else:   
        config=path_join(wl_dir, 'etc','wl.config')


    config = os.path.expandvars(config)

    # source list can be determined from tilename/band
    if srclist is None:
        srclist=\
            deswl.files.wlme_srclist_path(dataset,localid,tilename,band,serun)

    # info on all the unique, newest tiles and catalogs
    stdout.write('\n')
    infodict, tileinfo_url=\
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
            coaddimage=ti['image_url']
            coaddcat=ti['cat_url']
            break

    if coaddimage is None:
        raise ValueError("Tilename '%s' not found in '%s'" % \
                         (tilename, tileinfo_url))

    fdict=deswl.files.generate_me_filenames(tilename, band, 
                                            merun=merun, dir=outdir, 
                                            rootdir=rootdir)
    multishear_url=fdict['multishear']


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
    stat['TMV_DIR'] = tmv_dir
    stat['ESUTIL_DIR'] = esutil_dir
    stat['WL_DIR'] = wl_dir
    stat['DESFILES_DIR'] = desfiles_dir
    stat['executable'] = executable
    stat['tileinfo_url'] = tileinfo_url
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
    stat['multishear_url'] = multishear_url
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
    stdout.write('    TMV_DIR: %s\n' % tmv_dir)
    stdout.write('    ESUTIL_DIR: %s\n' % esutil_dir)
    stdout.write('    WL_DIR: %s\n' % wl_dir)
    stdout.write('    DESFILES_DIR: %s\n' % desfiles_dir)
    stdout.write('    executable: %s\n' % executable)
    stdout.write('    tileinfo file: %s\n' % tileinfo_url)
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
    json_util.write(stat, fdict['stat'])


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
                stdout.write(' * copying %s:%s to\n    %s\n' % (hostshort,f,df))
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

    rc=deswl.files.Runconfig(merun)
    dataset=rc['dataset']
    localid=rc['localid']
    serun=rc['serun']

    if tilelist is None:
        infodict=deswl.files.collated_coaddfiles_read(dataset, band, 
                                                      localid=localid,
                                                      serun=serun)
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
        fdict=deswl.files.generate_me_filenames(tilename, band, 
                                                merun=merun, rootdir=rootdir)

        
        statfile=fdict['stat']
        qafile=fdict['qa']
        msfile=fdict['multishear']

        error_found=False
        missing=[]
        if not os.path.exists(statfile):
            missing.append('statfile')
            error_found=True
        if not os.path.exists(qafile):
            missing.append('qa')
            error_found=True
        if not os.path.exists(msfile):
            missing.append('multishear')
            error_found=True

        if len(missing) > 0:
            missmess=','.join(missing)
            stdout.write('%s-%s: %s missing\n' % (tilename,band,missmess))
            stdout.write('    %s\n' % fdict['stat'])

        if 'statfile' not in missing:
            stat=json_util.read(statfile)
            # left over from old xml format
            exit_status=int( stat['exit_status'] )
            if exit_status != 0:
                stdout.write("%s-%s: Found non-zero exit "
                             "status %s in stat file %s\n" % \
                                (tilename,band,exit_status,statfile))
                error_found=True
            
        if error_found:
            badlist.append(ti)

    if badlist_file is not None:
        json_util.write(badlist, badlist_file)
    return badlist
        
