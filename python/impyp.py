import os
import deswl
import esutil as eu

_impyp_files={'shear':'%(run)s-%(expname)s-%(ccd)02d-shear.dat'}

def impyp_dir(run, expname, **keys):
    rc=deswl.files.Runconfig()

    fileclass=rc.run_types['impyp']['fileclass']
    rundir=deswl.files.run_dir(fileclass, run, **keys)
    dir=os.path.join(rundir, expname)
    return dir

def impyp_url(run, expname, ccd, ftype, **keys):
    if ftype not in _impyp_files:
        raise ValueError("bad impyp ftype: '%s'" % ftype)

    basename=_impyp_files[ftype] % {'run':run,
                                    'expname':expname,
                                    'ccd':int(ccd)}

    dir = impyp_dir(run, expname, **keys)

    url = os.path.join(dir, basename)
    return url




def generate_impyp_filenames(run, expname, ccd, **keys):
    """
    Output filenames for impyp
    """
    fdict={}

    rc=deswl.files.Runconfig()

    # output file names
    for ftype in _impyp_files:
        name= impyp_url(run, expname, ccd, ftype, **keys)
        fdict[ftype] = name

    return fdict


