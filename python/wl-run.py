"""
Run a single image through either measure_stars, measurepsf or measureshear
"""

import os
import sys
import subprocess
import time
import signal

# How long to wait for sub process to finish.  Make configurable
timeout = 120

finfo={}
finfo['cat'] = {'suffix':'_cat','ext':'fits'}

finfo['allcat'] = {'suffix':'-all','ext':'csv'}
finfo['starcat'] = {'suffix':'-star','ext':'csv'}

finfo['psf'] = {'suffix':'-psf','ext':'csv'}
finfo['psf-log'] = {'suffix':'-psf-log','ext':'csv'}
finfo['psf-debug'] = {'suffix':'-psf-debug','ext':'txt'}
finfo['fitpsf'] = {'suffix':'-fitpsf','ext':'csv'}

finfo['shear'] = {'suffix':'-shear','ext':'csv'}
finfo['shear-log'] = {'suffix':'-shear-log','ext':'csv'}
finfo['shear-debug'] = {'suffix':'-shear-debug','ext':'txt'}

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

def GetConf():
    wlconf=prefix+"/share/wl/wl.config"
    return wlconf

def MakeCommand(image_file, executable):

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
    print 'usage: '+sys.argv[0]+' image_file executable'
    sys.exit(45)

image_file = sys.argv[1]
executable = sys.argv[2]

cmd = MakeCommand(image_file, executable)

print cmd[0]
for c in cmd[1:]:
    print '\t'+c
command = ' '.join(cmd)

ExecuteCommand(command)
