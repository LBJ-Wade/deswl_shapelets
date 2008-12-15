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
#subdirs=['python','example_config']
subdirs=['python','example_config']

# Configurations will be saved here so command line options don't
# have to be sent more than once
config_file = 'wl_scons.conf'

# Some extra path places to look.  For extra_paths it is assumed that /bin
# /lib and /include come after the specified path.  For extra_include_paths
# it is assumed this is the full path to the /include directory (or whatever).
# Similarly for extra_library_paths and extra_bin_paths
# Note there is also the option to import paths from the shell environment
extra_paths = []
extra_include_paths = []
extra_library_paths = []
extra_bin_paths = []

# Default directory for installation.  
# This is the only UNIX specific things I am aware
# of in the script.  On the other hand, these are not required for the
# script to work since prefix can be set on the command line and the
# extra paths are not needed, but I wish I knew how to get the default 
# prefix for the system so I didn't have to set this.

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
opts.Add(BoolOption('DBG','Turn on debugging',True))
opts.Add(BoolOption('WITH_LAPACK',
                    'Look for lapack libraries and link if found.',
                    True))
opts.Add(BoolOption('WITH_BLAS',
                    'Look for blas libraries and link if found.',
                    True))
opts.Add(BoolOption('WITH_OPENMP',
                    'Look for openmp and use if found.',
                    True))
opts.Add(BoolOption('STATIC',
                    'Use static linkage',
                    False))
opts.Add(BoolOption('IMPORT_PATHS',
                    'Import PATH, C_INCLUDE_PATH and LIBRARY_PATH/LD_LIBRARY_PATH environment variables',
                    False))
opts.Add(BoolOption('IMPORT_ENV',
                    'Import full environment from calling shell',
                    False))
opts.Add(PathOption('EXTRA_PATH',
                    'Extra paths for executables (separated by : if more than 1)',
                    '',PathOption.PathAccept))
opts.Add(PathOption('EXTRA_LIB_PATH',
                    'Extra paths for linking (separated by : if more than 1)',
                    '',PathOption.PathAccept))
opts.Add(PathOption('EXTRA_INCLUDE_PATH',
                    'Extra paths for header files (separated by : if more than 1)',
                    '',PathOption.PathAccept))


opts.Update(initial_env)
opts.Save(config_file,initial_env)
Help(opts.GenerateHelpText(initial_env))

# Extra libraries to link in
initial_env['_extralibs'] = []


# This helps us determine of openmp is available
openmp_mingcc_vers = 4.2
openmp_minicc_vers = 9.0
openmp_minpgcc_vers = 6.0

def RunInstall(env, targets, subdir):
    install_dir = os.path.join(env['PREFIX'], subdir)
    env.Alias(target='install',
              source=env.Install(dir=install_dir, source=targets))

def RunUninstall(env, targets, subdir):
    # There is no env.Uninstall method, we must build our own
    install_dir = os.path.join(env['PREFIX'], subdir)
    deltarget = Delete("$TARGET")

    # delete from $prefix/bin/
    files = []
    for t in targets:
        ifile = os.path.join(install_dir, os.path.basename(str(t))) 
        files.append(ifile)

    for f in files:
        env.Alias('uninstall', env.Command(f, None, deltarget))


def CCFlags(env):
    """
    """
    # making this a list makes it easier to extend
    # dynamically
    compiler = os.path.basename(env['CXX'])
    version = env['CXXVERSION_NUMERICAL']
    if compiler[0] == 'g':
        cflags=['-O2','-Wall','-Werror','-ggdb']
        if version <= 4.2:
            cflags += ['-fno-strict-aliasing']
    elif compiler[0] == 'i':
        cflags=['-O2','-Wall','-Werror','-g','-wd383,810,981']
        if version >= 10:
            cflags += ['-wd1572']
    elif compiler[0] == 'p':
        cflags=['-O2','-fast','-Mcache_align','-g']
        #env.Append(CPPDEFINES=['USE_STEDC'])
    else:
        raise ValueError,'Need to add CCFLAGS for compilers other than g++, icpc, pgCC'
    return cflags

def AddOpenMPFlag(config):
    """
    Make sure you do this after you have determined the version of
    the compiler.

    gcc uses -fopemnp and an extra library
    icc uses -openmp
    pgCC uses -mp
    
    Other compilers?
    """
    compiler = os.path.basename(config.env['CXX'])
    version = config.env['CXXVERSION_NUMERICAL']
    if compiler[0] == 'g':
        # g++
        # In principle we should be able to look for libgomp, which is the
        # underlying library, but this isn't in the paths often.  All I can
        # figure out is a version check.
        #if config.CheckLib('gomp'):
        if version < openmp_mingcc_vers: 
            return
	if version == openmp_mingcc_vers:
	    print 'Warning: OpenMP works with g++ 4.2.3, but not with 4.2.1.'
	    print 'Currently, we only check for the first decimal, so if you'
	    print 'have 4.2.1, this might not compile.'
	    print '(The status of 4.2.2 is unknown.)'
        flag = '-fopenmp'
        ldflag = '-fopenmp'
    elif compiler[0] == 'i':
        # icpc
        if version < openmp_minicc_vers:
            return
        flag = '-openmp'
        ldflag = '-openmp'
    elif compiler[0] == 'p':
        # pgCC
        if version < openmp_minicc_vers:
            return
        flag = '-mp'
        ldflag = '-mp'
    else:
        print 'No OpenMP support for compiler ',compiler

    #print '\tAdding openmp flag:',flag
    print 'Using OpenMP'
    config.env.Append(CXXFLAGS=[flag])
    config.env.Append(LINKFLAGS=[ldflag])
    config.env['_extralibs'] += ['pthread']

def NDebugFlag(compiler):
    """
    """
    return 'NDEBUG'

def GetCompilerVersion(compiler):
    """
    """
    import re
    lines = os.popen(compiler + ' --version').readlines()
    # pgCC puts the version number on the second line of output.
    if compiler[0] == 'p':
        line = lines[1]
    else:
        line = lines[0]
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
    bin_paths = env['EXTRA_PATH'].split(':')
    lib_paths = env['EXTRA_LIB_PATH'].split(':')
    cpp_paths = env['EXTRA_INCLUDE_PATH'].split(':')

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

    AddPath(lib_paths, os.path.join(env['PREFIX'], 'lib'))
    AddPath(cpp_paths, os.path.join(env['PREFIX'], 'include'))

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

    #print 'LIBPATH = ',env['LIBPATH']
    #print 'CPPPATH = ',env['CPPPATH']

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

    # Mike Jarvis' matrix libraries
    if not config.CheckLibWithHeader('tmv','TMV.h',language='C++'):
        print 'tmv library or TMV.h not found'
        Exit(1)
    if not config.CheckLibWithHeader('tmv_symband','TMV.h',language='C++'):
        print 'tmv_symband library not found'
        Exit(1)

    # We need cfitsio in the search path
    if not config.CheckLibWithHeader('cfitsio','fitsio.h',language='C'):
        print 'cfitsio not found'
        Exit(1)

    # This bit needs to be the same as in the TMV SConstruct file
    # to make sure the correct libraries are linked
    # We don't need the CPPDEFINES here though.
    compiler = os.path.basename(config.env['CXX'])
    if config.env['WITH_BLAS']:
        # Look for ACML, MKL before more generic (and probably less
        # optimized) ATLAS library.
        if (compiler[0] == 'i' and not config.env['STATIC'] and
                (config.CheckLib('mkl',language='C++') or
                 config.CheckLib('mkl',language='C') )) :
            if config.env['WITH_LAPACK']:
                print 'Using MKL LAPACK'
                config.env['_extralibs'] += ['mkl_lapack']
            else:
	        pass
                #config.env.Append(CPPDEFINES=['NOLAP'])
            print 'Using MKL BLAS'
            #config.env.Append(CPPDEFINES=['MKL'])
            config.env['_extralibs'] += ['mkl']

        elif (compiler[0] == 'i' and 
                 (config.CheckLib('mkl_ia32',language='C++') or 
                  config.CheckLib('mkl_ia32',language='C') )) :
            if config.env['WITH_LAPACK']:
                print 'Using MKL LAPACK'
                config.env['_extralibs'] += ['mkl_lapack']
            else:
	        pass
                #config.env.Append(CPPDEFINES=['NOLAP'])
            print 'Using MKL BLAS'
            #config.env.Append(CPPDEFINES=['MKL'])
            config.env['_extralibs'] += ['mkl_ia32','mkl_sequential',
	        'mkl_core','guide','pthread']

        elif (compiler[0] == 'i' and 
                 (config.CheckLib('mkl_intel_lp64',language='C++')  or
                  config.CheckLib('mkl_intel_lp64',language='C') )):
            if config.env['WITH_LAPACK']:
                print 'Using MKL LAPACK'
                config.env['_extralibs'] += ['mkl_lapack']
            else:
	        pass
                #config.env.Append(CPPDEFINES=['NOLAP'])
            print 'Using MKL BLAS'
            #config.env.Append(CPPDEFINES=['MKL'])
	    # MJ -- This is untested.  I don't have a 64 bit machine with
	    #       icpc to test this on.
            config.env['_extralibs'] += ['mkl_intel_lp64','mkl_sequential',
	        'mkl_core','guide','pthread']

        elif (compiler[0] == 'p' and 
                (config.CheckLibWithHeader('acml','acml.h',language='C++') or
                 config.CheckLibWithHeader('acml','acml.h',language='C')) ):
            if config.env['WITH_LAPACK']:
                print 'Using ACML LAPACK'
            else:
	        pass
                #config.env.Append(CPPDEFINES=['NOLAP'])
            print 'Using ACML BLAS'
            #config.env.Append(CPPDEFINES=['ACML'])
            config.env['_extralibs'] += ['acml','pgftnrtl']

        elif ( (config.CheckLib('atlas',language='C++') or
	        config.CheckLib('atlas',language='C')) and
	       (config.CheckLibWithHeader('cblas','cblas.h',language='C++') or
		config.CheckLibWithHeader('cblas','cblas.h',language='C'))) :
            if config.env['WITH_LAPACK']:
                print 'No LAPACK support for non MKL, ACML option.'
                    #config.env.Append(CPPDEFINES=['NOLAP'])
            #config.env.Append(CPPDEFINES=['NOLAP'])
            print 'Using ATLAS BLAS'
            #config.env.Append(CPPDEFINES=['ATLAS'])
            config.env['_extralibs'] += ['cblas','atlas']
        else:
            print 'No BLAS libraries found'
            #config.env.Append(CPPDEFINES=['NOBLAS'])
    else:
	pass
        #config.env.Append(CPPDEFINES=['NOBLAS'])

    if config.env['STATIC']:
        config.env.Append(LINKFLAGS=['-static'])



def DoConfig(env):
    """
    Configure the system
    """
    if env['IMPORT_ENV']:
        # I couldn't figure out how to get this option before the 
        # initial constructor.  So this seems a bit inefficient to me.
        # But I think it works, so good enough for now.
        env = Environment(ENV=os.environ)
        # Now repeat the stuff that has already been done to env
        opts.Update(env)
        opts.Save(config_file,env)
        Help(opts.GenerateHelpText(env))
        env['_extralibs'] = []

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

    # work with multiple processors
    if config.env['WITH_OPENMP']:
        AddOpenMPFlag(config)

    DoLibraryAndHeaderChecks(config)

    if not env['DBG']:
        print 'Debugging turned off'
        flag = NDebugFlag(cxx)
        config.env.Append(CPPDEFINES=[flag])

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
    env['_InstallProgram'] = RunInstall
    env['_UninstallProgram'] = RunUninstall

    # subdirectores to process.  We process src by default
    script_files = [os.path.join(src_dir,'SConscript')]
    for d in subdirs:
        script_files.append(os.path.join(d,'SConscript'))

    SConscript(script_files, exports='env')


