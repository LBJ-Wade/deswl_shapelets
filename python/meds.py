"""
Defines the MEDS class to work with MEDS (Multi Epoch Data Structures)

See docs for the MEDS class for more info
"""
import fitsio

class MEDS(object):
    """
    Class to work with MEDS (Multi Epoch Data Structures)

    One can extract cutouts using get_cutout() and get_mosaic() and
    get_cutout_list()

    One can access fields using [field_name] notation. The number of entries is
    in the .size attribute. Access the .cat catalog for further functionality
    with regards the recarray

    Fields

     id                 i4       id from coadd catalog
     ncutout            i4       number of cutouts for this object
     box_size           i4       box size for each cutout
     file_id            i4[NMAX] zero-offset id into the file names in the 
                                 second extension
     start_row          i4[NMAX] zero-offset, points to start of each cutout.
     orig_row           f8[NMAX] zero-offset position in original image
     orig_col           f8[NMAX] zero-offset position in original image
     orig_start_row     i4[NMAX] zero-offset start corner in original image
     orig_start_col     i4[NMAX] zero-offset start corner in original image
     cutout_row         f8[NMAX] zero-offset position in cutout imag
     cutout_col         f8[NMAX] zero-offset position in cutout image

    The first cutout for an object is always from the coadd.

    examples
    --------
    import meds
    m=meds.MEDS(filename)

    # number of coadd objects
    num=m.size

    # number of cutouts for object 35
    m['ncutout'][35]

    # get cutout 3 for object 35
    im=m.get_cutout(35,3)

    # get the all the cutouts for object 35 as a single image
    im=m.get_mosaic(35)

    # get a list of all the cutout images for object 35
    im=m.get_cutout_list(35)

    # the last two share the same underlying storage for the images

    # you can access any of the columns in the
    # catalog (stored as a recarray) directly

    # e.g. get the center in the cutout for use in image processing
    row = m['row_cutout'][35]
    col = m['col_cutout'][35]

    """
    def __init__(self, filename):
        self._filename=filename
        
        self._fits=fitsio.FITS(filename)

        self._cat=self._fits[1][:]

    def get_cutout(self, iobj, icutout):
        """
        Get a single cutout for the indicated entry

        parameters
        ----------
        iobj:
            Index of the object
        icutout:
            Index of the cutout for this object.

        returns
        -------
        The cutout image
        """
        self._check_indices(iobj, icutout=icutout)

        box_size=self._cat['box_size'][iobj]
        start_row = self._cat['start_row'][iobj,icutout]
        row_end = start_row+box_size

        return self._fits["image_cutouts"][start_row:row_end,:]

    def get_mosaic(self, iobj):
        """
        Get a mosaic of all cutouts associated with this coadd object

        parameters
        ----------
        iobj:
            Index of the object

        returns
        -------
        An image holding all cutouts
        """

        self._check_indices(iobj)

        ncutout=self._cat['ncutout'][iobj]
        box_size=self._cat['box_size'][iobj]

        start_row = self._cat['start_row'][iobj,0]
        row_end = start_row+box_size*ncutout

        mosaic=self._fits["image_cutouts"][start_row:row_end,:]
        return mosaic

    def get_cutout_list(self, iobj):
        """
        Get an image list with all cutouts associated with this coadd object

        parameters
        ----------
        iobj:
            Index of the object

        returns
        -------
        A list of images hold all cutouts.
        """

        mosaic=this.get_mosaic(iobj)
        ncutout=self._cat['ncutout'][iobj]
        box_size=self._cat['box_size'][iobj]
        return self._split_mosaic(mosaic, box_size, ncutout)

    def _split_mosaic(self, mosaic, box_size, ncutout):
        imlist=[]
        for i in xrange(ncutout):
            r1=i*box_size
            r2=(i+1)*box_size
            imlist.append( mosaic[r1:r2, :] )

        return imlist


    def _check_indices(self, iobj, icutout=None):
        if iobj >= self._cat.size:
            raise ValueError("object index should be within "
                             "[0,%s)" % self._cat.size)

        ncutout=self._cat['ncutout'][iobj]
        if ncutout==0:
            raise ValueError("object %s has no cutouts" % iobj)

        if icutout is not None:
            if icutout >= ncutout:
                raise ValueError("requested cutout index %s for "
                                 "object %s should be in bounds "
                                 "[0,%s)" % (iobj,ncutout))

    def __repr__(self):
        return repr(self._fits[1])
    def __getitem__(self, item):
        return self._cat[item]

    @property
    def size(self):
        return self._cat.size
