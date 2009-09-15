from sys import stdout,stderr
try:
    import oracle_util
    import cx_Oracle as cxo
except:
    stderr.write('Could not import oracle utilities\n')

def GetTileRun(tilename, verbose=False, multi=False):
    """
    Try to get the appropriate run identifier for the input tilename
    Currently choses the latest one
    """

    #runreg='DES%'
    #query="""
    #    SELECT 
    #        distinct(run) 
    #    FROM 
    #        location
    #    WHERE
    #        tilename like '%s'
    #        AND run like '%s'
    #    ORDER BY
    #        run DESC
    #""" % (tilename, runreg)

    query="""
        SELECT 
            distinct(run) 
        FROM 
            location
        WHERE
            tilename like '%s'
        ORDER BY
            run DESC
    """ % (tilename,)
    

    if verbose:
        sys.stdout.write('%s\n' % query)

    o=oracle_util.Connection()
    res = o.Execute(query)
    if res is None:
        return None

    if multi:
        return res['run']
    else:
        run = res['run'][0]
        return run

def GetCoaddId(tilename, band, run=None, verbose=False):
    if run is None:
        run = GetTileRun(tilename, verbose=verbose)
        if run is None:
            return None

    query="""
        SELECT
            id
        FROM
            coadd
        WHERE
            tilename like '%s'
            AND run like '%s'
            AND band like '%s'
    """ % (tilename, run, band)

    if verbose:
        sys.stdout.write('%s\n' % query)
    o=oracle_util.Connection()
    res = o.Execute(query)
    if res is None:
        return None
    return res['id']





def CoaddCatalogInfo(band=None, dictlist=False):
    """
    we demand the band since band is not set in the database for most
    tables and we want to set it explicitly

    get all the DES coadd catalog info
    """

    q="""
    select
        *
    from
        catalog
    where
        catalogtype = 'coadd_cat' 
        and project = 'DES'
        and catalogname like '%%_%s_cat%%'
    order by
        catalog.run
    """ % band

    stdout.write(q)

    conn = oracle_util.Connection()
    res = conn.Execute(q, dictlist=dictlist, lower=True)

    if res is not None:
        if dictlist:
            for i in range(len(res)):
                res[i]['band'] = band
        else:
            res['band'] = band
    return res


def GetCoaddSrcLocations(coadd_id=None, tilename=None, band=None, run=None, 
                         dictlist=False, verbose=False):
    """
    Get info about original SE 'red' images associated with a given coadd
    from the locations table.  This can be a complicated process.

    depends:
        (GetTileRun)
        (GetCoaddId)

    """

    if coadd_id is None:
        coadd_id = GetCoaddId(tilename, band, run=run, verbose=verbose)
        coadd_id=coadd_id[0]


    # Get a list of coadd_src images for the coadd_id
    query="""
        SELECT
            image.parentid
        FROM
            image,coadd_src
        WHERE
            coadd_src.coadd_imageid = %d
            AND coadd_src.src_imageid = image.id
    """ % coadd_id

    if verbose:
        sys.stdout.write('%s\n' % query)
    o = oracle_util.Connection()
    res = o.Execute(query,lower=True)

    if res is None:
        return None

    # now we need to find the single epoch red images we are looking for
    # typically we now have the ids of the remap images, which means the very
    # next query for parents will get 'red' but we might have to iterate

    idstrings = [str(id) for id in res['parentid']]
    ftype = 'none'
    i=0
    itmax=5

    while (ftype != 'red') and (i < itmax):
        idlist = ', '.join(idstrings)
        query="""
            SELECT
                image.id,
                image.imagetype,
                image.parentid
            FROM
                image
            WHERE
                image.id in (%s)
        """ % idlist

        if verbose:
            sys.stdout.write('%s\n' % query)
        res = o.Execute(query, lower=True)
        if res is None:
            return None

        idstrings = [str(id) for id in res['parentid']]
        ftype = res['imagetype'][0]
        i += 1

    if ftype != 'red':
        raise ValueError,"Reach itmax=%s before finding 'red' images. last is %s" % (itmax, ftype)

    # now the idlist comes from id instead of parentid
    idstrings = [str(id) for id in res['id']]
    idlist = ', '.join(idstrings)
    query="""
        select
            location.id,
            location.run,
            location.exposurename,
            location.filename,
            location.band,
            location.filetype,
            location.ccd
        from
            location
        where
            location.id in (%s)
    """ % idlist

    if verbose:
        sys.stdout.write('%s\n' % query)
    res = o.Execute(query, dictlist=dictlist, lower=True)

    if verbose:
        stdout.write("# iterations: %s\n" % i) 
    return res



def CollateCoaddCatalogsAndImages(band, dictlist=False,
                                  toxml=None, 
                                  getsrc=False):
    """

    join the catalog table (catalogtype='coadd_cat' and project='DES') to the
    coadd table, using parentid as the connection

    Note this join is useful even if we didn't need the info since there are
    entries in the database for which the id info is nonsense and we have to
    skip them.

    we demand the band since band is not set in the database for many entries
    so we must set it explicitly

    if getsrc=True add a list of the source 'red' images that made up the
    coadd.  This depends on
        DBGetCoaddSrcLocations
            (GetCoaddId)
                (GetTileRun)

    """

    q="""
    select
        catalog.id as catalogid,
        catalog.run as catalogrun,
        catalog.band,
        catalog.tilename,
        catalog.catalogname,
        catalog.parentid,
        location.run as coaddrun,
        location.filename as coaddname
    from
        catalog, location
    where
        catalog.catalogtype = 'coadd_cat' 
        and catalog.project = 'DES'
        and location.id = catalog.parentid
        and catalogname like '%%_%s_cat%%'
    order by
        catalog.run
    """ % band

    stdout.write(q+'\n')

    conn = oracle_util.Connection()

    if toxml or getsrc:
        dictlist=True

    res = conn.Execute(q, dictlist=dictlist, lower=True)

    if res is not None:
        if dictlist:
            for i in range(len(res)):
                res[i]['band'] = band
        else:
            res['band'] = band


    if getsrc:
        # add the single epoch coadd source images
        stdout.write('Getting src lists\n')
        nbad=0
        for i in range(len(res)):
            coadd_id = res[i]['parentid']
            srclist = DBGetCoaddSrcLocations(coadd_id,dictlist=True)
            if srclist is None:
                raise ValueError,\
                    'Could not get source list for id=%s\n' % coadd_id
            ftype_bad = \
                [sl['filetype'] for sl in srclist if sl['filetype'] != 'red']
            if len(ftype_bad) > 0:
                raise ValueError,"Expected 'red' type got '%s'\n" % ftype_bad[0]

            # fix up the bands
            for sl in srclist:
                if sl['band'] != band:
                    stdout.write('Fixing band\n')
                    sl['band'] = band
            res[i]['srclist'] = srclist


    if toxml is not None:
        stdout.write('Writing to xml file: %s\n' % toxml)
        xmltools.dict2xml({'info':res}, toxml, roottag='coaddinfo')
    return res


def DumpQuery(q, fname):
    import sfile
    conn = oracle_util.Connection()
    stdout.write(q+'\n')
    res = conn.Execute(q)
    stdout.write("Writing to file: %s\n" % fname)
    sfile.write(res, fname)
    res=0
    return

def DumpQueryCurs(q, sep=',', fobj=None, fname=None):
    """
    currently only good for numbers since not fixed length fields
    """
    import sys
    import cx_Oracle as cxo

    if fobj is None and fname is None:
        fobj = sys.stdout
    elif fname is not None:
        stdout.write("Writing to file: %s\n" % fname)
        fobj = open(fname, 'w')

    defconn=os.environ['ORACLE_CONNINFO']
    conn = cxo.connect(defconn)
    curs = conn.cursor()

    stdout.write(q+'\n')
    curs.execute(q)
    
    for res in curs:
        line = [str(r) for r in res]
        line = sep.join(line) +'\n'
        fobj.write(line)

    conn.close()
    return




def DumpTableInChunks(tables_in, columns_in, fname, chunksize=1000000,
                      constraints='',
                      dryrun=False):
    import sfile
    conn = oracle_util.Connection()

    if isinstance(tables_in, str):
        tables = [tables_in]
    else:
        tables = tables_in
    if isinstance(columns_in, str):
        columns = [columns_in]
    else:
        columns = columns_in



    q="select count(*) from %s" % tables[0]
    stdout.write(q+'\n')
    res = conn.Execute(q)
    nrows = res['count(*)'][0]

    nchunk = nrows/chunksize
    nleft = nrows % chunksize

    stdout.write('Total number of rows: %s\n' % nrows)
    stdout.write('Chunksize: %s\n' % chunksize)
    stdout.write('    Main chunks: %s\n' % nchunk)
    stdout.write('    leftover:    %s\n' % nleft)
    stdout.write("Output file: '%s'\n" % fname)

    nget = nchunk
    if nleft > 0:
        nget += 1

    colstr = ','.join(columns)
    tabstr = ','.join(tables)
    for i in range(nget):

        rowmin = i*chunksize + 1
        rowmax = (i+1)*chunksize + 1
        q = """
        select 
            %s
        from
            %s
        where
            rownum >= %s
            AND rownum < %s
            %s
        """ % (colstr, tabstr, rowmin, rowmax, constraints)

        stdout.write(q+'\n')
        res = conn.Execute(q)

        if res is not None:
            stdout.write("Found: %s\n" % len(res))
        else:
            stderr.write("Found no rows\n");

        if not dryrun and res is not None:
            if i == 0:
                sfile.write(res, fname)
            else:
                sfile.write(res, fname, append=True)

    conn.Close()


def GetIndexNorun(loc):
    """
    This is darren's index that does not include the run.
    """
    ind = []
    keys = ['project','fileclass','filetype', 'exposurename', 'filename', 
            'tilename', 'band', 'ccd']
    for i in range(len(loc)):
        tind = []
        for key in  keys:
            tind.append( str(loc[i][key]).strip() )
        tind = '-'.join(tind)
        ind.append(tind)
    return ind



def QueryByValue(table, colname, values):
    """
    Assumes these are unique
    """

    conn = oracle_util.Connection()

    n = len(values)

    chunksize = 1000
    nchunk = n/chunksize
    nleft = n % chunksize

    res = None
    if nleft > 0:
        nget = nchunk+1
    else:
        nget = nchunk
    for i in range(nget):
        vals = [str(val) for val in values[i*chunksize:(i+1)*chunksize]]
        valstring = ','.join(vals)
        q = """
        select 
            *
        from
            %s
        where %s in (%s)
        """ % (table, colname, valstring)
        
        tres = conn.Execute(q)
        if res is None:
            res = numpy.zeros(n, dtype=tres.dtype)

        res[i*chunksize:i*chunksize + len(tres)] = tres
        
    return res




def GetLatestDC4Locations(limit=None, band=None):
    """
    At this stage we just need the files with nite like '2008%'.

    Then we just get the lastest version of the file (not sure
    at this point why there are multiple). That is the step we
    do for removing dups but keeping the latest.
    """

    conn = oracle_util.Connection()

    nite = '2008%'


    stdout.write('Getting images corresponding to catalogs\n')

    cat_loc_join = """
    SELECT
        catalog.id as cat_id, 
        catalog.parentid,
        location.id as image_id,
        location.project,
        location.fileclass,
        location.filetype,
        location.exposurename,
        location.filename,
        location.tilename,
        location.band,
        location.ccd,
        location.run
    FROM
        catalog, location
    WHERE
        catalog.project = 'DES'
        AND catalog.catalogtype = 'red_cat'
        AND catalog.nite LIKE '%s'
        AND catalog.parentid = location.id
    """ % (nite,)

    if band is not None:
        cat_loc_join += "    AND location.band = '%s'" % band

    stdout.write(cat_loc_join+'\n')
    imloc = conn.Execute(cat_loc_join)

    stdout.write("Returned %d rows\n" % imloc.size)

    stdout.write('Getting latest runs\n') 
    ind = GetLatestRuns(imloc)
    imloc = imloc[ind]


    catloc = QueryByValue('location','id',imloc['cat_id'])

    return imloc, catloc
    

    # First get the images
    image_query = """
    SELECT
        *
    FROM
        location
    WHERE
        project = 'DES'
        AND filetype='red'
        AND nite LIKE '%s'
    """ % (nite,)


    if limit is not None:
        image_query += ' ' + limit

     
    stdout.write(image_query+'\n')
    imloc = conn.Execute(image_query)

    stdout.write('Getting latest runs\n')
    ind = GetLatestRuns(imloc)
    imloc = imloc[ind]

        
    return imloc, catloc




