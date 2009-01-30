#!/usr/bin/env python
"""
    wl-run.py executable image_file

Run a single image through either findstars, measurepsf or measureshear
"""

import os
import sys
import subprocess
import time
import signal
import re


from optparse import OptionParser

parser = OptionParser(__doc__)

parser.add_option("--cat-file", default=None,
                  help="input catalog file")
parser.add_option("--stars-file", default=None,
                  help="Findstars output file")
parser.add_option("--psf-file", default=None,
                  help="measurepsf PSF output file")
parser.add_option("--fitpsf-file", default=None,
                  help="measurepsf fitted PSF output file")
parser.add_option("--shear-file", default=None,
                  help="measureshear output file")

parser.add_option("-t","--timeout", default=300,
                  help="timeout in seconds. Default is %default")

#parser.add_option("-c","--clobber",action="store_true", dest="verbose",
#                  help="Clobber, or overwrite, existing files")

parser.add_option("-v","--verbose",action="store_true", dest="verbose",
                  help="print messages.  By default "+\
                    "only the stdout/stderr of the process is printed")

extstrip=re.compile( '\.fits(\.fz)?$' )

bitshift = 8

def ExecuteCommand(command, timeout, verbose=False):
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
        res = 1024
    else:
        if res > 100:
            if verbose:
                print 'This system is probably returning bit-shifted '+\
                        'exit codes. Old value:',res
            res = res >> bitshift
    if verbose:
        print 'exit code:',res

    return 1024

def MakeRootPrefix(image_file):
    image_root = extstrip.sub('', image_file)
    root = os.path.basename(image_root)
    input_prefix = os.path.dirname(image_root)
    if input_prefix != '':
        input_prefix += os.sep
    return root, input_prefix

def MakeRoot(image_file):
    image_root = extstrip.sub('', image_file)
    return image_root

def MakeCommand(executable, conf_file, image_file, 
                cat_file=None,
                stars_file=None, 
                psf_file=None, fitpsf_file=None, 
                shear_file=None):

    #root, input_prefix = MakeRoot(image_file)
    root = MakeRoot(image_file)

    command = \
        [executable,
         conf_file,
         'root='+root]

    #if input_prefix != '':
        #     command.append('input_prefix='+input_prefix)

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

options, args = parser.parse_args(sys.argv[1:])

if len(args) < 2:
    parser.print_help()
    sys.exit(45)

executable = args[0]
conf_file = args[1]
image_file = args[2]
verbose=options.verbose
timeout=int(options.timeout)

cmd = MakeCommand(executable, conf_file, image_file,
                  cat_file=options.cat_file,
                  stars_file=options.stars_file,
                  psf_file=options.psf_file,
                  fitpsf_file=options.fitpsf_file,
                  shear_file=options.shear_file)

if verbose:
    print cmd[0] + '    \\'
    for c in cmd[1:]:
        print '    '+c+'    \\'
command = ' '.join(cmd)

res=ExecuteCommand(command, timeout, verbose=verbose)
sys.exit(res)
