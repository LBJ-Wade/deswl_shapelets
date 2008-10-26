# vim: set filetype=python :
# to do:  
#   Add more support for other compilers
#   will support g++ and icc at least
# 
# Always run scons from the root directory of tmv

import os
import sys


# Subdirectories containing SConscript files.  We always process src but
# there are some other optional ones
src_dir = 'src'
subdirs=[]

# Configurations will be saved here so command line options don't
# have to be sent more than once
config_file = 'wl.conf'

# The following path pieces are the only UNIX specific things I am aware
# of in the script.  On the other hand, these are not required for the
# script to work since prefix can be set on the command line and the
# extra paths are not needed.

# Some extra path places to look.  For extra_paths it is assumed that /bin
# /lib and /include come after the specified path.  For extra_include_paths
# it is assumed this is the full path to the /include directory (or whatever).
# Similarly for extra_library_paths and extra_bin_paths
# Note there is also the option to import paths from the shell environment
#extra_paths = ['/sw.unstable','/sw','/opt/local']
extra_paths = []
extra_include_paths = []
extra_library_paths = []
extra_bin_paths = []

# Default directory for installation.  
default_prefix = '/usr/local'


# Now set up the environment
initial_env = Environment()

# first check for a saved conf file
opts = Options(config_file)

# Now set up options for the command line
opts.Add(PathOption('PREFIX',
                    'prefix for installation',
                    default_prefix,
                    PathOption.PathAccept))
opts.Add('CXX','Name of c++ compiler',None)
opts.Add(BoolOption('DBG','Turn on debugging output',False))
opts.Add(BoolOption('WITH_ATLAS',
                    'Look for atlas libraries and link if found.',
                    True))
opts.Add(BoolOption('IMPORT_PATHS',
                    'Import PATH, C_INCLUDE_PATH and LIBRARY_PATH/LD_INCLUDE_PATH environment variables',
                    False))

opts.Update(initial_env)
opts.Save(config_file,initial_env)
Help(opts.GenerateHelpText(initial_env))

# Extra libraries to link in
initial_env['_extralibs'] = []


# This helps us determine of openmp is available
openmp_mingcc_vers = 4.2
openmp_minicc_vers = 8.0



def CCFlags(env):
    """
    Make this check for icc or other compilers in addition to these gcc
    specific settings
    """
    # making this a list makes it easier to extend
    # dynamically
    compiler = os.path.basename(env['CXX'])
    if compiler[0] == 'g':
        gflags=['-O2','-fno-strict-aliasing']
    else:
        raise ValueError,'Need to add CCFLAGS for compilers other than gcc'
    return gflags

def AddOpenMPFlag(env):
    """
    Make sure you do this after you have determined the version of
    the compiler.

    gcc uses -fopemnp and an extra library, -lgomp
    icc uses -openmp, library for linking?  
    
    Other compilers?
    """
    compiler = os.path.basename(env['CXX'])
    version = env['CXXVERSION_NUMERICAL']
    if compiler[0] == 'g':
        # g++
        if version < openmp_mingcc_vers: 
            return
        flag = '-fopenmp'
        ldflag = '-fopenmp'
    elif compiler[0] == 'i':
        # icc
        if version < openmp_minicc_vers:
            return
        flag = '-openmp'
        # library?
    else:
        raise ValueError,'compiler must be g++ or icc'

    print 'Adding openmp support:',flag
    env.Append(CXXFLAGS=[flag])
    env.Append(LINKFLAGS=[ldflag])

def NDebugFlag(compiler):
    """
    need to make work with other compilers
    """
    return '-DNDEBUG'

def GetCompilerVersion(compiler):
    """
    I think this should work for both gcc and icc, but need to check icc
    """
    import re
    line = os.popen(compiler + ' --version').readline()
    match = re.search(r'[0-9]+(\.[0-9]+)+', line)
    if match:
        return match.group(0) 

def GetNumericalVersion(version):
    # Get the version up to the first decimal, e.g. for 4.3.1 we only keep 4.3
    vnum = version[0:version.find('.')+2]
    return float(vnum)

def AddPath(pathlist, newpath):
    """
    Add path(s) to a list of paths.  Check the path exists and that it is
    not already in the list
    """
    if type(newpath) == list:
        for l in newpath:
            AddPath(pathlist, l)
    else:
        # to deal with expansions and possible end / which 
        # messes up uniqueness test
        p = os.path.abspath(newpath) 
        if os.path.exists(p):
            if pathlist.count(p) == 0:
                pathlist.append(p)

def AddExtraPaths(env):
    """
    Add some include and library paths.
    Also merge in $PATH, $C_INCLUDE_PATH and $LIBRARY_PATH/$LD_LIBRARY_PATH 
    environment variables if requested.
    
    The set itself is created in order of appearance here, but then this 
    whole set is prepended.  The order within this list is:

        paths from the user's environment
        paths in extra_bin_paths, extra_library_paths, and extra_include_paths
        paths under extra_paths

    Only paths that actually exists are kept.
    """
    bin_paths = []
    cpp_paths = []
    lib_paths = []

    # These come first
    if env['IMPORT_PATHS'] and os.environ.has_key('PATH'):
        paths=os.environ['PATH']
        paths=paths.split(os.pathsep)
        AddPath(bin_paths, paths)

    if env['IMPORT_PATHS'] and os.environ.has_key('C_INCLUDE_PATH'):
        paths=os.environ['C_INCLUDE_PATH']
        paths=paths.split(os.pathsep)
        AddPath(cpp_paths, paths)

    if env['IMPORT_PATHS'] and os.environ.has_key('LIBRARY_PATH'):
        paths=os.environ['LIBRARY_PATH']
        paths=paths.split(os.pathsep)
        AddPath(lib_paths, paths)

    if env['IMPORT_PATHS'] and os.environ.has_key('LD_LIBRARY_PATH'):
        paths=os.environ['LD_LIBRARY_PATH']
        paths=paths.split(os.pathsep)
        AddPath(lib_paths, paths)

    for p in extra_include_paths:
        AddPath(cpp_paths, p)
    for p in extra_library_paths:
        AddPath(lib_paths, p)
    for p in extra_bin_paths:
        AddPath(bin_paths, p)

    for expath in extra_paths:
        AddPath(bin_paths,os.path.join(expath,'bin')) 
        AddPath(cpp_paths,os.path.join(expath,'include')) 
        AddPath(lib_paths,os.path.join(expath,'lib')) 

    #env.AppendENVPath('PATH', bin_paths)
    #env.Append(LIBPATH= lib_paths)
    #env.Append(CPPPATH= cpp_paths)
    env.PrependENVPath('PATH', bin_paths)
    env.Prepend(LIBPATH= lib_paths)
    env.Prepend(CPPPATH= cpp_paths)

def ReadFileList(fname):
    """
    This reads a list of whitespace separated values from the input file fname
    and stores it as a list.  We will make this part of the environment so
    other SConscripts can use it
    """
    try:
        files=open(fname).read().split()
    except:
        print 'Could not open file:',fname
        sys.exit(45)
    files = [f.strip() for f in files]
    return files


def DoLibraryAndHeaderChecks(config):
    """
    Check for some headers.  In some cases we may just sat a flag during
    compilation, such as if no blas is found we set -DNOBLAS.  In other
    cases we may raise an exception or just exit
    """
    # Add cfitsio here
    if not config.CheckLibWithHeader('cfitsio','fitsio.h','C'):
        print 'cfitsio not found'
        Exit(1)
    if config.env['WITH_ATLAS'] and config.CheckLib('atlas'):
        if config.CheckLibWithHeader('cblas','cblas.h','C'):
            config.env.Append(CXXFLAGS=['-DATLAS'])
            # I've had lots of issues trying to get lapack support to work.
            # MJ says his routines are nearly as fast, so let's just skip it
            config.env.Append(CXXFLAGS=['-DNOLAP'])
            # The order is important here
            #config.env.Append(LINKFLAGS='-lcblas -latlas')
            #config.env.Append(LINKFLAGS='/sw.unstable/lib/libcblas.a /sw.unstable/lib/libatlas.a')
            config.env['_extralibs'] += ['cblas','atlas']

    # here is a silly example of exiting using the Exit() scons function
    #if not config.CheckHeader('stdio.h'):
    #    print 'stdio.h not found'
    #    Exit(1)

def DoConfig(env):
    """
    Configure the system
    """
    config = env.Configure()

    # local includes and lib paths
    cpppath=['#include']
    libpath=['#lib']

    config.env.Append(CPPPATH=cpppath)
    config.env.Append(LIBPATH=libpath)

    # add some extra paths 
    AddExtraPaths(config.env)

    # The basic flags for this compiler
    cxx=config.env['CXX']
    config.env['CXXVERSION'] = GetCompilerVersion(cxx)
    config.env['CXXVERSION_NUMERICAL'] = \
            GetNumericalVersion(config.env['CXXVERSION'])
    v=config.env['CXXVERSION']
    print '\nUsing compiler:',cxx
    print 'compiler version:',v

    # Get the version up to the first decimal, e.g. for 4.3.1 we only keep 4.3
    v = v[0:v.find('.')+2]

    cxxflags=CCFlags(config.env)
    config.env.Append(CXXFLAGS=[cxxflags])

    # Not supporting openmp in tmv at this moment
    #AddOpenMPFlag(config.env)

    DoLibraryAndHeaderChecks(config)

    if not env['DBG']:
        print 'Debugging turned off'
        flag = NDebugFlag(cxx)
        config.env.Append(CXXFLAGS=[flag])

    print

    # Now finish configuration and get the new environment
    env = config.Finish()
    return env


#
# main program
#

if not GetOption('help'):

    # Set up the configuration
    env = DoConfig(initial_env)

    # subdirectory SConscript files can use this function
    env['__readfunc'] = ReadFileList

    # subdirectores to process.  We process src by default
    script_files = [os.path.join(src_dir,'SConscript')]
    for d in subdirs:
        if d in COMMAND_LINE_TARGETS:
            script_files.append(os.path.join(d,'SConscript'))

    SConscript(script_files, exports='env')


