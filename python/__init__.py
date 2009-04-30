# Don't edit these svn properties by hand
_property_headurl='$HeadURL$'

#_property_rev='$Rev: 1339 $'

#def GetRevision():
#    revision = _property_rev.split()[1]
#    return revision

def GetWlVersion():
    import re

    thisname='/python/__init__.py'
    headstrip=re.compile( '^\$HeadURL: .*/trunk/' )
    stripped_headurl = headstrip.sub('trunk/', _property_headurl)
    tag = stripped_headurl.replace(' $','')
    tag = tag.replace(thisname,'')

    return tag

version=GetWlVersion()
import wlpipe
