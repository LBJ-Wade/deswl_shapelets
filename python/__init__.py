
# version info
# You need to run 'svn propset svn:keywords HeadURL' on the file and commit
# before this works.
#
# Don't edit these svn properties by hand
_property_headurl='$HeadURL$'

def version():
    import os
    from sys import stderr

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

    tag = os.path.basename(urlfront)
    return tag


import wlpipe
import files
import oracle

