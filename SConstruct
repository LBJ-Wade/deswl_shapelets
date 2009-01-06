# vim: set filetype=python :
# to do:  
#   Add more support for other compilers
#   will support g++ and icpc at least
# 
# Always run scons from the root directory of tmv

import os
import sys


# Subdirectories containing SConscript files.  We always process src but
# there are some other optional ones
src_dir = 'src'
subdirs=['python','example_config']

# Configurations will be saved here so command line options don't
# have to be sent more than once
config_file = 'wl_scons.conf'

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
            'prefix for installation', default_prefix, PathOption.PathAccept))
opts.Add('CXX','Name of c++ compiler',None)
opts.Add('FLAGS','Compile flags to send to the compiler','')
opts.Add(BoolOption('DBG','Turn on debugging',True))
opts.Add(BoolOption('WITH_LAPACK',
            'Look for lapack libraries and link if found.', True))
opts.Add(BoolOption('WITH_BLAS',
            'Look for blas libraries and link if found.', True))
opts.Add(BoolOption('WITH_OPENMP',
            'Look for openmp and use if found.', True))
opts.Add(BoolOption('STATIC',
            'Use static linkage', False))
opts.Add(BoolOption('IMPORT_PATHS',
            'Import PATH, C_INCLUDE_PATH and LIBRARY_PATH/LD_LIBRARY_PATH environment variables',
            False))
opts.Add(BoolOption('IMPORT_ENV',
            'Import full environment from calling shell', False))
opts.Add(PathOption('EXTRA_PATH',
            'Extra paths for executables (separated by : if more than 1)',
            '',PathOption.PathAccept))
opts.Add(PathOption('EXTRA_LIB_PATH',
            'Extra paths for linking (separated by : if more than 1)',
            '',PathOption.PathAccept))
opts.Add(PathOption('EXTRA_INCLUDE_PATH',
            'Extra paths for header files (separated by : if more than 1)',
            '',PathOption.PathAccept))
opts.Add(BoolOption('MEM_TEST',
            'Test for memory leaks',False))


opts.Update(initial_env)
opts.Save(config_file,initial_env)
Help(opts.GenerateHelpText(initial_env))

# This helps us determine of openmp is available
openmp_mingcc_vers = 4.2
openmp_minicpc_vers = 9.0
openmp_minpgcc_vers = 7.0
# MJ -- pgCC 6.1 used to work with openmp, but I'm getting seg faults now.
#       I need to look into it a bit more.

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


def BasicCCFlags(env):
    """
    """
    if env['FLAGS'] == '':
        compiler = env['CXXTYPE']
        version = env['CXXVERSION_NUMERICAL']

        if compiler == 'g++':
            env.Replace(CCFLAGS=['-O2','-Wall','-Werror','-ggdb'])
            if version <= 4.2:
                  env['CCFLAGS'] += ['-fno-strict-aliasing']

        elif compiler == 'icpc':
            env.Replace(CCFLAGS=['-O2','-Wall','-Werror','-g','-wd383,810,981'])
            if version >= 9:
                env.Append(CCFLAGS=['-wd1572'])

        elif compiler == 'pgCC':
            env.Replace(CCFLAGS=['-O2','-fast','-Mcache_align','-g'])

        else:
            print 'Warning: Unknown compiler.  You should set FLAGS directly.'
            print 'Trying simple -O flag '
            env.Replace(CCFLAGS=['-O'])
    else :
        # If flags are specified as an option use them:
        cxx_flags = env['FLAGS'].split(' ')
        env.Replace(CCFLAGS=cxx_flags)


def AddOpenMPFlag(env):
    """
    Make sure you do this after you have determined the version of
    the compiler.

    gcc uses -fopemnp and an extra library
    icpc uses -openmp
    pgCC uses -mp
    
    Other compilers?
    """
    compiler = env['CXXTYPE']
    version = env['CXXVERSION_NUMERICAL']
    if compiler == 'g++':
        # In principle we should be able to look for libgomp, which is the
        # underlying library, but this isn't in the paths often.  All I can
        # figure out is a version check.
        #if config.CheckLib('gomp'):
        if version < openmp_mingcc_vers: 
	    print 'No OpenMP support for g++ versions before ',openmp_mingcc_vers
            return
        if version == openmp_mingcc_vers:
            print 'Warning: OpenMP works with g++ 4.2.3, but not with 4.2.1.'
            print 'Currently, we only check for the first decimal, so if you'
            print 'have 4.2.1, this might not compile.'
            print '(The status of 4.2.2 is unknown.)'
        flag = ['-fopenmp']
        ldflag = ['-fopenmp']
        xlib = ['pthread']
    elif compiler == 'icpc':
        if version < openmp_minicpc_vers:
	    print 'No OpenMP support for icpc versions before ',openmp_minicpc_vers
            return
        flag = ['-openmp']
        ldflag = ['-openmp']
        xlib = ['pthread']
    elif compiler == 'pgCC':
        if version < openmp_minpgcc_vers:
	    print 'No OpenMP support for pgCC versions before ',openmp_minpgcc_vers
            return
        flag = ['-mp','--exceptions']
        ldflag = ['-mp','-lpthread']
        xlib = ['pthread']
    else:
        print 'No OpenMP support for compiler ',compiler

    #print '\tAdding openmp flag:',flag
    print 'Using OpenMP'
    env.Append(CXXFLAGS=flag)
    env.Append(LINKFLAGS=ldflag)
    env.Append(LIBS=xlib)


def GetCompilerVersion(env):
    """
    """
    compiler = env['CXX']

    # Get the compiler type without suffix or path.  
    # e.g. /sw/bin/g++-4 -> g++
    if 'icpc' in compiler :
        compilertype = 'icpc'
    elif 'pgCC' in compiler :
        compilertype = 'pgCC'
    elif 'g++' in compiler :
        compilertype = 'g++'
    else :
        compilertype = 'unknown'

    lines = os.popen(compiler + ' --version').readlines()
    # pgCC puts the version number on the second line of output.
    if compilertype == 'pgCC':
        line = lines[1]
    else:
        line = lines[0]

    import re
    match = re.search(r'[0-9]+(\.[0-9]+)+', line)

    if match:
        version = match.group(0)
        # Get the version up to the first decimal
        # e.g. for 4.3.1 we only keep 4.3
        vnum = version[0:version.find('.')+2]
    else:
        version = 0
        vnum = 0

    print '\nUsing compiler:',compiler
    print 'compiler version:',version

    env['CXXTYPE'] = compilertype
    env['CXXVERSION'] = version
    env['CXXVERSION_NUMERICAL'] = float(vnum)


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

        paths in EXTRA_*PATH parameters
        paths from the user's environment
	local lib and include paths
        paths in PREFIX directory

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

    # local includes and lib paths
    # The # symbol means to interpret these from the top-level scons
    # directory even when we are in a sub-directory (src,test,etc.)
    cpp_paths += ['#include']
    lib_paths += ['#lib']

    # PREFIX directory
    AddPath(lib_paths, os.path.join(env['PREFIX'], 'lib'))
    AddPath(cpp_paths, os.path.join(env['PREFIX'], 'include'))

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


def CheckLibs(context,try_libs,source_file):
    init_libs = context.env['LIBS']
    context.env.Prepend(LIBS=try_libs)
    result = context.TryLink(source_file,'.cpp')
    if not result :
        context.env.Replace(LIBS=init_libs)
    return result
      

def CheckMKL(context):
    mkl_source_file = """
#include "mkl.h"
int main()
{
  char uplo='U', compq='I';
  int n=1,ldu=1,ldv=1,*iq=0,*iwork=0,*info=0;
  double *d=0, *e=0, *u=0, *v=0, *q=0, *work=0;
  dbdsdc(uplo,compq,n,d,e,u,ldu,v,ldv,q,iq,work,iwork,info);
  return 0;
}
"""

    context.Message('Checking for MKL... ')
    result = (
        (not context.env['STATIC'] and 
                CheckLibs(context,['mkl'],mkl_source_file)) or
        CheckLibs(context,['mkl_ia32','mkl_sequential','mkl_core',
	            'guide','pthread'],mkl_source_file) or
        CheckLibs(context,['mkl_intel_lp64','mkl_sequential','mkl_core',
	            'guide','pthread'],mkl_source_file) or
        CheckLibs(context,['mkl_ipf','guide','pthread'],mkl_source_file) or
        CheckLibs(context,['mkl_em64t','guide','pthread'],mkl_source_file) )
    context.Result(result)
    return result


def CheckACML(context):
    acml_source_file = """
#include "acml.h"
int main()
{
  char uplo='U', compq='I';
  int n=1,ldu=1,ldv=1,*iq=0,*iwork=0,*info=0;
  double *d=0, *e=0, *u=0, *v=0, *q=0, *work=0;
  dbdsdc(uplo,compq,n,d,e,u,ldu,v,ldv,q,iq,work,iwork,info);
  return 0;
}
"""

    context.Message('Checking for ACML... ')
    result = CheckLibs(context,['acml','pgftnrtl'],acml_source_file)
    context.Result(result)
    return result


def CheckGOTO(context):
    goto_source_file = """
extern "C" {
#include "fblas.h"
}
int main()
{
  char ta='N', tb='N';
  int M=1,N=1,K=1,lda=1,ldb=1,ldc=1;
  double alpha=1.,beta=1., *A=0, *B=0, *C=0;
  dgemm_(ta,tb,M,N,K,alpha,A,lda,B,ldb,beta,C,ldc,1,1);
  return 0;
}
"""

    context.Message('Checking for GotoBLAS... ')
    result = CheckLibs(context,['goto'],goto_source_file)
    context.Result(result)
    return result

def CheckATLAS(context):
    atlas_source_file = """
extern "C" {
#include "cblas.h"
}
int main()
{
  int M=1,N=1,K=1,lda=1,ldb=1,ldc=1;
  double alpha=1,beta=1, *A=0, *B=0, *C=0;
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
      M,N,K,alpha,A,lda,B,ldb,beta,C,ldc);
  return 0;
}
"""

    context.Message('Checking for ATLAS... ')
    result = CheckLibs(context,['cblas','atlas'],atlas_source_file)
    context.Result(result)
    return result


def CheckCBLAS(context):
    cblas_source_file = """
extern "C" {
#include "cblas.h"
}
int main()
{
  int M=1,N=1,K=1,lda=1,ldb=1,ldc=1;
  double alpha=1.,beta=1., *A=0, *B=0, *C=0;
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
          M,N,K,alpha,A,lda,B,ldb,beta,C,ldc);
  return 0;
}
"""

    context.Message('Checking for CBLAS... ')
    result = CheckLibs(context,['cblas'],cblas_source_file)
    context.Result(result)
    return result

def CheckFBLAS(context):
    fblas_source_file = """
extern "C" {
#include "fblas.h"
}
int main()
{
  char ta='N', tb='N';
  int M=1,N=1,K=1,lda=1,ldb=1,ldc=1;
  double alpha=1.,beta=1., *A=0, *B=0, *C=0;
  dgemm_(ta,tb,M,N,K,alpha,A,lda,B,ldb,beta,C,ldc,1,1);
  return 0;
}
"""

    context.Message('Checking for Fotran BLAS... ')
    result = CheckLibs(context,['blas'],fblas_source_file)
    context.Result(result)
    return result


def CheckATLAS_LAP(context):
    atlas_lapack_source_file = """
extern "C" {
#include "clapack.h"
}
int main()
{
  int N=1,lda=1;
  double* A=0;
  clapack_dlauum(CblasRowMajor,CblasUpper,N,A,lda);
  return 0;
}
"""

    context.Message('Checking for ATLAS version of LAPACK... ')
    result = CheckLibs(context,['lapack'],atlas_lapack_source_file)
    context.Result(result)
    return result


def CheckCLAPACK(context):
    clapack_source_file = """
extern "C" {
#include "f2c.h"
#include "clapack.h"
}
int main()
{
  char uplo='U', compq='I';
  int n=1,ldu=1,ldv=1,*iq=0,*iwork=0,*info=0;
  double *d=0, *e=0, *u=0, *v=0, *q=0, *work=0;
  dbdsdc_(&uplo,&compq,&n,d,e,u,&ldu,v,&ldv,q,iq,work,iwork,info);
  return 0;
}
"""

    context.Message('Checking for CLAPACK... ')
    result = CheckLibs(context,['lapack','cblaswr','f2c'],clapack_source_file)
    context.Result(result)
    return result


def CheckFLAPACK(context):
    flapack_source_file = """
extern "C" {
#include "flapack.h"
}
int main()
{
  char uplo='U', compq='I';
  int n=1,ldu=1,ldv=1,*iq=0,*iwork=0,*info=0;
  double *d=0, *e=0, *u=0, *v=0, *q=0, *work=0;
  dbdsdc_(uplo,compq,n,d,e,u,ldu,v,ldv,q,iq,work,iwork,info);
  return 0;
}
"""

    context.Message('Checking for Fortran LAPACK... ')
    result = CheckLibs(context,['lapack'],flapack_source_file)
    context.Result(result)
    return result


def DoLibraryAndHeaderChecks(config):
    """
    Check for some headers.  In some cases we may just set a flag during
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
    if not config.CheckLibWithHeader('cfitsio','fitsio.h',language='C++'):
        print 'cfitsio not found'
        Exit(1)

    # The rest needs to be the same as in the TMV SConstruct file
    # to make sure the correct libraries are linked
    # We don't need the CPPDEFINES or env['LAP'] here though.
    compiler = config.env['CXXTYPE']
    if not (config.env.has_key('LIBS')) :
      config.env['LIBS'] = []

    def DoCheckLAPACK(config) :
        if config.CheckCLAPACK() :
	    print 'Using CLAPACK'
	elif config.CheckFLAPACK() :
	    print 'Using Fortran LAPACK'
	else :
	    return 0
	return 1


    if config.env['WITH_BLAS']:
        if compiler == 'icpc' and config.CheckMKL() :
            if (config.env['WITH_LAPACK']) :
                print 'Using MKL LAPACK'
            print 'Using MKL BLAS'
 
        elif config.CheckACML() :
            if config.env['WITH_LAPACK']:
                print 'Using ACML LAPACK'
            print 'Using ACML BLAS'

        elif config.CheckGOTO() :
            if config.env['WITH_LAPACK'] and DoCheckLAPACK(config) :
	        pass
            print 'Using GotoBLAS'

        elif config.CheckCBLAS() :
            if config.env['WITH_LAPACK'] and DoCheckLAPACK(config) :
	        pass
            print 'Using CBLAS'

        elif config.CheckATLAS() :
            if config.env['WITH_LAPACK'] and DoCheckLAPACK(config) :
	        pass
            elif config.env['WITH_LAPACK'] and config.CheckATLAS_LAP() :
	        print 'Using ATLAS LAPACK'
            print 'Using ATLAS'

        elif config.CheckFBLAS() :
            if config.env['WITH_LAPACK'] and DoCheckLAPACK(config) :
	        pass
            print 'Using Fortran BLAS'

        else:
            print 'No BLAS libraries found'

    if config.env['STATIC'] :
        config.env.Append(LINKFLAGS=['static'])


def DoConfig(env):
    """
    Configure the system
    """


    # add some extra paths 
    AddExtraPaths(env)

    # Figure out what kind of compiler we are dealing with
    GetCompilerVersion(env)

    # The basic flags for this compiler if not explicitly specified
    BasicCCFlags(env)

    # Some extra flags depending on the options:
    if env['WITH_OPENMP']:
        AddOpenMPFlag(env)
    if not env['DBG']:
        print 'Debugging turned off'
        env.Append(CPPDEFINES='NDEBUG')
    if env['MEM_TEST']:
        env.Append(CPPDEFINES=['MEM_TEST'])

    import SCons.SConf

    # Figure out what BLAS and/or LAPACK libraries are on the system
    # MJ: I have had bad luck with scons figuring out when the cache
    #     is invalid.  This just forces a check every time.
    SCons.SConf.SetCacheMode('force')
    config = env.Configure(custom_tests = {
        'CheckMKL' : CheckMKL ,
        'CheckACML' : CheckACML ,
        'CheckGOTO' : CheckGOTO ,
        'CheckATLAS' : CheckATLAS ,
        'CheckCBLAS' : CheckCBLAS ,
        'CheckFBLAS' : CheckFBLAS ,
        'CheckATLAS_LAP' : CheckATLAS_LAP ,
        'CheckCLAPACK' : CheckCLAPACK ,
        'CheckFLAPACK' : CheckFLAPACK })
    DoLibraryAndHeaderChecks(config)
    env = config.Finish()
    SCons.SConf.SetCacheMode('auto')



#
# main program
#

if not GetOption('help'):

    env = initial_env

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

    # Set up the configuration
    DoConfig(env)

    # subdirectory SConscript files can use this function
    env['__readfunc'] = ReadFileList
    env['_InstallProgram'] = RunInstall
    env['_UninstallProgram'] = RunUninstall

    # subdirectores to process.  We process src by default
    script_files = [os.path.join(src_dir,'SConscript')]
    for d in subdirs:
        script_files.append(os.path.join(d,'SConscript'))

    SConscript(script_files, exports='env')


