# vim: set filetype=python et ts=4 sw=4:
# to do:  
#   Add more support for other compilers
#   will support g++ and icpc at least
# 
# Always run scons from the root directory of wl

import os
import sys
from sys import stdout,stderr


# Subdirectories containing SConscript files.  We always process src but
# there are some other optional ones
src_dir = 'src'
subdirs=['python','example_config','bin']

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

# This helps us determine of openmp is available
openmp_mingcc_vers = 4.1
openmp_minicpc_vers = 9.0
openmp_minpgcc_vers = 6.0

# Now set up the environment
initial_env = Environment()

# first check for a saved conf file
opts = Variables(config_file)

# Now set up options for the command line
opts.Add('CXX','Name of c++ compiler')
opts.Add('FLAGS','Compile flags to send to the compiler','')
opts.Add(BoolVariable('DEBUG','Turn on debugging statements',True))
opts.Add(PathVariable('PREFIX','prefix for installation','', PathVariable.PathAccept))

opts.Add(PathVariable('EXTRA_PATH',
            'Extra paths for executables (separated by : if more than 1)',
            '',PathVariable.PathAccept))
opts.Add(PathVariable('EXTRA_LIB_PATH',
            'Extra paths for linking (separated by : if more than 1)',
            '',PathVariable.PathAccept))
opts.Add(PathVariable('EXTRA_INCLUDE_PATH',
            'Extra paths for header files (separated by : if more than 1)',
            '',PathVariable.PathAccept))
opts.Add(BoolVariable('IMPORT_PATHS',
            'Import PATH, C_INCLUDE_PATH and LIBRARY_PATH/LD_LIBRARY_PATH environment variables',
            False))
opts.Add(BoolVariable('IMPORT_ENV','Import full environment from calling shell', False))

opts.Add('TMV_LINK','File that contains the linking instructions for TMV','')
opts.Add('LIBS','Libraries to send to the linker','')
opts.Add(BoolVariable('CACHE_LIB','Cache the results of the library checks',True))

opts.Add(BoolVariable('WITH_OPENMP','Look for openmp and use if found.', True))
opts.Add(BoolVariable('WITH_TMV','Use TMV for Matrix/Vector (rather than Eigen)', True))
opts.Add(BoolVariable('STATIC','Use static linkage', False))
opts.Add(BoolVariable('MEM_TEST','Test for memory leaks',False))
opts.Add(BoolVariable('WARN','Add warning compiler flags, like -Wall', False))

opts.Add(BoolVariable('WITH_UPS',
            'Create ups/wl.table.  Install the ups directory under PREFIX/ups',False))
opts.Add(BoolVariable('WITH_PROF',
            'Use the compiler flag -pg to include profiling info for gprof',False))


opts.Update(initial_env)
opts.Save(config_file,initial_env)
Help(opts.GenerateHelpText(initial_env))

def RunInstall(env, targets, subdir):
    install_dir = os.path.join(env['INSTALL_PREFIX'], subdir)
    env.Alias(target='install',
              source=env.Install(dir=install_dir, source=targets))

def RunUninstall(env, targets, subdir):
    # There is no env.Uninstall method, we must build our own
    install_dir = os.path.join(env['INSTALL_PREFIX'], subdir)
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
            env.Replace(CCFLAGS=['-O2'])
            if version <= 4.2:
                env.Append(CCFLAGS=['-fno-strict-aliasing'])
            if env['WARN']:
                #env.Append(CCFLAGS=['-ggdb','-Wall','-Werror'])
                env.Append(CCFLAGS=['-g','-Wall','-Werror'])
    
        elif compiler == 'icpc':
            env.Replace(CCFLAGS=['-O2'])
            if version >= 10:
                env.Append(CCFLAGS=['-vec-report0'])
            if env['WARN']:
                env.Append(CCFLAGS=['-Wall','-Werror','-wd383,810,981'])
                if version >= 9:
                    env.Append(CCFLAGS=['-wd1572'])
                if version >= 11:
                    env.Append(CCFLAGS=['-wd2259'])

        elif compiler == 'pgCC':
            env.Replace(CCFLAGS=['-O2','-fast','-Mcache_align'])

        elif compiler == 'cl':
            env.Replace(CCFLAGS=['/EHsc','/nologo','/O2','/Oi'])
            if env['WARN']:
                env.Append(CCFLAGS=['/W2','/WX'])

        else:
            print 'Warning: Unknown compiler.  You should set FLAGS directly.'
            env.Replace(CCFLAGS=[])

    else :
        # If flags are specified as an option use them:
        cxx_flags = env['FLAGS'].split(' ')
        env.Replace(CCFLAGS=cxx_flags)

    # Also parse the LIBS options if present
    if env['LIBS'] == '':
        env.Replace(LIBS=[])
    else:
        libs = env['LIBS'].split(' ')
        env.Replace(LIBS=libs)

    if env['WITH_PROF']:
        env.Append(CCFLAGS=['-pg'])
        env.Append(LINKFLAGS=['-pg'])


def AddOpenMPFlag(env):
    """
    Make sure you do this after you have determined the version of
    the compiler.

    g++ uses -fopemnp
    icpc uses -openmp
    pgCC uses -mp
    
    Other compilers?
    """
    compiler = env['CXXTYPE']
    version = env['CXXVERSION_NUMERICAL']
    if compiler == 'g++':
        if version < openmp_mingcc_vers: 
            print 'No OpenMP support for g++ versions before ',openmp_mingcc_vers
            return
        flag = ['-fopenmp']
        #ldflag = []
        #xlib = ['gomp','gfortran','pthread']
        ldflag = ['-fopenmp']
        # Note: gcc_eh is required on MacOs, but not linux
        # However, it seems to be safe for both (???), so I always include it.
        # If it turns out to break something, I'll have to revisit this.
        xlib = ['pthread','gcc_eh']
        # Also add a #define to let program know that linking will use openmp
        # even for object files that don't.
        env.Append(CPPDEFINES=['OPENMP_LINK'])
    elif compiler == 'icpc':
        if version < openmp_minicpc_vers:
            print 'No OpenMP support for icpc versions before ',openmp_minicpc_vers
            return
        #env.Append(CCFLAGS=['-openmp'])
        #flag = []
        flag = ['-openmp']
        ldflag = ['-openmp']
        xlib = ['pthread']
        env.Append(CPPDEFINES=['OPENMP_LINK'])
    elif compiler == 'pgCC':
        if version < openmp_minpgcc_vers:
            print 'No OpenMP support for pgCC versions before ',openmp_minpgcc_vers
            return
        flag = ['-mp','--exceptions']
        ldflag = ['-mp']
        xlib = ['pthread']
    elif compiler == 'cl':
        flag = ['/openmp']
        ldflag = ['/openmp']
        xlib = []
    else:
        print 'Warning: No OpenMP support for compiler ',compiler

    #print 'Adding openmp support:',flag
    env['OMP_FLAGS'] = flag
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
        versionflag = '--version'
        linenum=0
    elif 'pgCC' in compiler :
        compilertype = 'pgCC'
        versionflag = '--version'
        linenum=1
        # pgCC puts the version number on the second line of output.
    elif 'g++' in compiler :
        compilertype = 'g++'
        versionflag = '--version'
        linenum=0
    elif 'cl' in compiler :
        compilertype = 'cl'
        versionflag = ''
        linenum=0
        # With cl, the version seems to be printed in the first line,
        # but the lines read out with popen below seem to skip the
        # first two lines.  So the code below ends up with version = 0.
        # It doesn't really matter though, since we don't use the cl 
        # version for anything. 
    else :
        compilertype = 'unknown'
        version = 0
        vnum = 0

    if compilertype != 'unknown':
        cmd = compiler + ' ' + versionflag
        lines = os.popen(cmd).readlines()
        line = lines[linenum]
    
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

        local lib and include paths
        paths in EXTRA_*PATH parameters
        paths from the user's environment
        paths in PREFIX directory

    Only paths that actually exists are kept.
    """
    # local includes and lib paths
    # The # symbol means to interpret these from the top-level scons
    # directory even when we are in a sub-directory (src,test,etc.)
    bin_paths = []
    cpp_paths = ['#include']
    lib_paths = ['#lib']

    # Paths specified in EXTRA_*
    bin_paths += env['EXTRA_PATH'].split(':')
    lib_paths += env['EXTRA_LIB_PATH'].split(':')
    cpp_paths += env['EXTRA_INCLUDE_PATH'].split(':')

    # Paths found in environment paths
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

    # PREFIX directory
    # If none given, then don't add them to the -L and -I directories.
    # But still use the default /usr/local for installation
    if env['PREFIX'] == '':
        env['INSTALL_PREFIX'] = default_prefix
    else:
        AddPath(bin_paths, os.path.join(env['PREFIX'], 'bin'))
        AddPath(lib_paths, os.path.join(env['PREFIX'], 'lib'))
        AddPath(cpp_paths, os.path.join(env['PREFIX'], 'include'))
        env['INSTALL_PREFIX'] = env['PREFIX']
    

    #print 'bin paths = ',bin_paths
    #print 'cpp paths = ',cpp_paths
    #print 'lib paths = ',lib_paths

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


def CheckLibs(context,try_libs,source_file):
    init_libs = context.env['LIBS']
    context.env.Prepend(LIBS=try_libs)
    result = context.TryLink(source_file,'.cpp')
    if not result :
        context.env.Replace(LIBS=init_libs)
    return result
      

def CheckTMV(context):
    tmv_source_file = """
#include "TMV_Sym.h"
int main()
{
  tmv::SymMatrix<double> S(10,4.);
  tmv::Matrix<double> m(10,3,2.);
  tmv::Matrix<double> m2 = m / S;
  return 0;
}
"""

    print 'Checking for correct TMV linkage... (this may take a little while)'
    context.Message('Checking for correct TMV linkage... ')

    if not context.TryCompile(tmv_source_file,'.cpp'):
        context.Result(0)
        print 'Error: TMV file failed to compile.'
        Exit(1)

    if not CheckLibs(context,['tmv_symband','tmv'],tmv_source_file):
        # If that didn't work, we might need to add the openmp flag to the 
        # linking step.
        if not env['WITH_OPENMP']:
            env1 = context.env.Clone()
            AddOpenMPFlag(env1)
            context.env['LINKFLAGS'] = env1['LINKFLAGS']
            if not CheckLibs(context,['tmv_symband','tmv'],tmv_source_file):
                context.Result(0)
                print 'Error: TMV file failed to link correctly'
                Exit(1)
        else:
            context.Result(0)
            print 'Error: TMV file failed to link correctly'
            Exit(1)
    context.Result(1)
    return 1


def CheckEigen(context):
    eigen_source_file = """
#include "Eigen/Core"
#include "Eigen/Cholesky"
int main()
{
  Eigen::MatrixXd S(10,10);
  S.setConstant(4.);
  Eigen::MatrixXd m(10,3);
  m.setConstant(2.);
  Eigen::MatrixXd m2(10,3);
  S.llt().solve(m,&m2);
  return 0;
}
"""

    print 'Checking for correct Eigen linkage... (this may take a little while)'
    context.Message('Checking for correct Eigen linkage... ')

    if not context.TryCompile(eigen_source_file,'.cpp'):
        context.Result(0)
        print 'Error: Eigen file failed to compile.'
        Exit(1)

    # Eigen doesn't actually need any linking, but I didn't bother to 
    # change the format of this check.  So [] for the linking libraries.
    if not CheckLibs(context,[],eigen_source_file):
        context.Result(0)
        print 'Error: Eigen file failed to link correctly'
        Exit(1)

    context.Result(1)
    return 1



def DoLibraryAndHeaderChecks(config):
    """
    Check for some headers.  
    """

    # We need cfitsio in the search path
    if not config.CheckLibWithHeader('cfitsio','fitsio.h',language='C++'):
        print 'cfitsio not found'
        Exit(1)

    if not config.CheckLibWithHeader('CCfits','CCfits/CCfits.h',
                                     language='C++'):
        stdout.write('CCfits library or header not found\n')
        Exit(1)

    if config.env['WITH_TMV'] :
        # First do a simple check that the library and header are in the path.
        # We check the linking with the BLAS library below.
        if not config.CheckLibWithHeader('tmv','TMV.h',language='C++',autoadd=0):
            print 'tmv library or TMV.h not found'
            Exit(1)
        if not config.CheckLibWithHeader('tmv_symband','TMV_Sym.h',language='C++',
                autoadd=0):
            print 'tmv_symband library not found'
            Exit(1)

        compiler = config.env['CXXTYPE']
        version = config.env['CXXVERSION_NUMERICAL']
    
        if not (config.env.has_key('LIBS')) :
            config.env['LIBS'] = []
    
        if (config.env['TMV_LINK'] == '') :
            tmv_link_file = config.env['INSTALL_PREFIX'] + '/share/tmv-link'
        else :
            tmv_link_file = config.env['TMV_LINK']
    
        print 'Using TMV_LINK file:',tmv_link_file
        try:
            tmv_link = open(tmv_link_file).read().strip()
        except:
            print 'Could not open TMV link file: ',tmv_link_file
            Exit(1)
        print '    ',tmv_link
        config.env.Append(LINKFLAGS=tmv_link)

        config.CheckTMV()
        config.env.Append(CPPDEFINES=['USE_TMV'])

    else :
        # No need for the BLAS libraries if we are using Eigen.

        if not config.CheckHeader('Eigen/Core',language='C++') :
            print 'Eigen/Core not found'
            Exit(1)

        config.CheckEigen()


def DoConfig(env):
    """
    Configure the system
    """


    # Add some extra paths 
    AddExtraPaths(env)

    # Figure out what kind of compiler we are dealing with
    GetCompilerVersion(env)
   
    # The basic flags for this compiler if not explicitly specified
    BasicCCFlags(env)

    # Some extra flags depending on the options:
    if env['WITH_OPENMP']:
        print 'Using OpenMP'
        AddOpenMPFlag(env)
    if not env['DEBUG']:
        print 'Debugging turned off'
        env.Append(CPPDEFINES=['NDEBUG'])
    if env['MEM_TEST']:
        env.Append(CPPDEFINES=['MEMTEST'])
    if env['STATIC'] :
        if env['CXXTYPE'] == 'pgCC':
            env.Append(LINKFLAGS=['-Bstatic'])
        else:
            env.Append(LINKFLAGS=['-static'])

    import SCons.SConf

    # Figure out what BLAS and/or LAPACK libraries are on the system
    # MJ: I have had bad luck with scons figuring out when the cache
    #     is invalid.  This just forces a check every time.
    if not env['CACHE_LIB']:
        SCons.SConf.SetCacheMode('force')
    config = env.Configure(custom_tests = {
        'CheckTMV' : CheckTMV ,
        'CheckEigen' : CheckEigen
        })
    DoLibraryAndHeaderChecks(config)
    env = config.Finish()
    # MJ: Turn the cache back on now, since we want it for the
    #     main compilation steps.
    if not env['CACHE_LIB']:
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

    # Set up the configuration
    DoConfig(env)

    # subdirectory SConscript files can use this function
    env['__readfunc'] = ReadFileList
    env['_InstallProgram'] = RunInstall
    env['_UninstallProgram'] = RunUninstall

    if env['WITH_UPS']:
        subdirs += ['ups']

    # subdirectores to process.  We process src by default
    script_files = [os.path.join(src_dir,'SConscript')]
    for d in subdirs:
        script_files.append(os.path.join(d,'SConscript'))

    SConscript(script_files, exports='env')


