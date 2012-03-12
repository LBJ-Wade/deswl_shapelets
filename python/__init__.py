import sys
from sys import stdout

def version():
    from sys import stderr

    # version info
    # You need to run 'svn propset svn:keywords HeadURL' on the file and commit
    # before this works.
    #
    # Don't edit these svn properties by hand
    _property_headurl='$HeadURL$'



    thisname='/python/__init__.py'
    badvers="NOTAG: unparseable"

    psplit=_property_headurl.split()
    if len(psplit) != 3:
        mess="headurl did not split into 3: '%s'\n" % _property_headurl
        stderr.write(mess)
        return badvers

    url=psplit[1]

    if url.find(thisname) == -1:
        mess="url '%s' does not contain string '%s'\n" % \
                (_property_headurl, thisname)
        stderr.write(mess)
        return badvers

    urlfront = url.replace(thisname, '')

    tag=urlfront.split('/')[-1]
    return tag

def get_external_version(version_command):
    """
    Run a command to get a version string and return it through the
    standard output
    """
    import subprocess
    pobj=subprocess.Popen(version_command, shell=True, 
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # this will wait until process terminates
    out=pobj.communicate()
    sout=out[0]
    serr=out[1]
    estatus=pobj.poll()
    if sout == '' or estatus != 0:
        raise RuntimeError("Could not get version from "
                           "command: %s: %s" % (version_command,serr))
    return sout.strip()

def get_tmv_version():
    return get_external_version('tmv-version')

def get_python_version(numerical=False):
    if numerical:
        v=sys.version_info[0:3]
        pyvers=v[0] + 0.1*v[1] + 0.01*v[2]
    else:
        pyvers='v%s.%s.%s' % sys.version_info[0:3]
    return pyvers


try:
    import cwl
    from cwl import CWL as WL
except:
    stdout.write('Could not import cwl\n')
    pass


import wlpipe
import files

import generic
import impyp


