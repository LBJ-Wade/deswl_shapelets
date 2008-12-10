#!/usr/bin/env python
"""
Run a single image through either findstars, measurepsf or measureshear
"""

import os
import sys
import subprocess
import time
import signal

# How long to wait for sub process to finish.  Make configurable
timeout = 120

finfo={}

crap="""
decam--24--29-i-2_52_star.csv
decam--24--29-i-2_52_findstars_log.csv
decam--24--29-i-2_52_all.csv
decam--24--29-i-2_52_psf.csv
decam--24--29-i-2_52_measurepsf_log.csv
decam--24--29-i-2_52_fitpsf.dat
decam--24--29-i-2_52_shear.csv
decam--24--29-i-2_52_measureshear_log.csv
"""

finfo['cat']         = {'suffix':'_cat','ext':'fits'}

finfo['allcat']      = {'suffix':'_all','ext':'csv'}
finfo['starcat']     = {'suffix':'_star','ext':'csv'}

finfo['psf']         = {'suffix':'_psf','ext':'csv'}
finfo['psf_log']     = {'suffix':'_psf_log','ext':'csv'}
finfo['psf_debug']   = {'suffix':'_psf_debug','ext':'txt'}
finfo['fitpsf']      = {'suffix':'_fitpsf','ext':'dat'}

finfo['shear']       = {'suffix':'_shear','ext':'csv'}
finfo['shear_log']   = {'suffix':'_shear_log','ext':'csv'}
finfo['shear_debug'] = {'suffix':'_shear_debug','ext':'txt'}

# these should be configurable
prefix="/home/esheldon/local"
data_root = "/data/Archive"
output_root="/home/esheldon/data"

bitshift = 8

def MakeOutputFroot(root):
    return root.replace(data_root, output_root)

def MakeFileName(root, ftype):
    if not finfo.has_key(ftype):
        raise ValueError,'No such file type: '+str(ftype)
    info=finfo[ftype]
    return root+info['suffix']+'.'+info['ext']

def MakeAllFileNames(image_file):
    fn = {}

    # root of file name without various extensions
    froot = image_file.replace('.fz','')
    froot = froot.replace('.fits','')

    # root in the output directory
    output_froot = MakeOutputFroot(froot)
    output_dir = os.path.dirname(output_froot)

    fn['froot'] = froot
    fn['output_froot'] = output_froot
    fn['output_dir'] = output_dir

    # cat is in the froot always
    fn['cat'] = MakeFileName(froot, 'cat')

    # these may be in an alternative directory
    for ftype in finfo.keys():
        fn[ftype] = MakeFileName(output_froot, ftype)

    return fn

def GetExecutableArgs(executable, fnames):
    cmd = []
    if executable == 'measurepsf':
        cmd.append('starcat_file='+fnames['starcat'])
        cmd.append('debug_file='+fnames['psf-debug'])
        cmd.append('log_file='+fnames['psf-log'])
        cmd.append('psf_file='+fnames['psf'])
        cmd.append('fitpsf_file='+fnames['fitpsf'])
    elif executable == 'measureshear':
        cmd.append('allcat_file='+fnames['allcat'])
        cmd.append('fitpsf_file='+fnames['fitpsf'])
        cmd.append('log_file='+fnames['shear-log'])
        cmd.append('debug_file='+fnames['shear-debug'])
        cmd.append('shear_file='+fnames['shear'])
    else:
        raise ValueError,'Unkown executable: '+str(executable)
    return cmd


def MakeCommandOld(image_file, executable):

    fnames = MakeAllFileNames(image_file)

    if not os.path.exists(fnames['output_dir']):
        # create all intermediate components of path
        os.makedirs(fnames['output_dir'])

    # Just run measurepsf right now
    cmd = []
    cmd.append(executable)
    cmd.append(GetConf())
    cmd.append('root='+fnames['froot'])
    cmd.append('image_file='+image_file)
    extra_cmd = GetExecutableArgs(executable, fnames)
    cmd += extra_cmd

    return cmd



def GetConf():
    wlconf=prefix+"/share/wl/wl.config"
    return wlconf

def MakeNamePieces(image_file):
    base=os.path.basename(image_file)

    # root of file name without various extensions
    froot = base.replace('.fz','')
    froot = froot.replace('.fits','')

    input_prefix = os.path.dirname(image_file)
    output_prefix = input_prefix.replace(data_root, output_root)

    return froot, input_prefix, output_prefix


def MakeCommand(executable, image_file):
    root, input_prefix, output_prefix = MakeNamePieces(image_file)

    cmd = []

    cmd.append(executable)
    cmd.append(GetConf())
    cmd.append('root='+root)
    cmd.append('log_ext=_'+executable+'_log.csv')

    if input_prefix != '':
        input_prefix+=os.sep
        cmd.append('input_prefix='+input_prefix)

    if output_prefix != '':
        if not os.path.exists(output_prefix):
            # create all intermediate components of path
            os.makedirs(output_prefix)
        output_prefix+=os.sep
        cmd.append('output_prefix='+output_prefix)
    return cmd


def ExecuteCommand(command):
    pobj = subprocess.Popen(command, shell=True)

    # Create a new sub process
    # end the process after timeout seconds
    # poll interval is 1 second for simplicity
    # but for faster processes may want to change
    # this
    i=0
    while i < timeout:
        res = pobj.poll()
        if res is not None:
            break
        time.sleep(1)
        i += 1

    if res is None:
        print 'Process is taking longer than',timeout,'seconds.  Ending process'
        os.kill(pobj.pid, signal.SIGTERM)
    else:
        if res > 100:
            print 'This system is probably returning bit-shifted '+\
                    'exit codes. Old value:',res
            res = res >> bitshift
    print 'exit code:',res


if len(sys.argv) < 3:
    print 'usage: '+os.path.basename(sys.argv[0])+' executable image_file'
    sys.exit(45)

executable = sys.argv[1]
image_file = sys.argv[2]

cmd = MakeCommand(executable, image_file)

print cmd[0]
for c in cmd[1:]:
    print '\t'+c
command = ' '.join(cmd)

ExecuteCommand(command)
