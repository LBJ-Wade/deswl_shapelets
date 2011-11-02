"""
    %prog merun

Generate the input list for multishear. Each row is
    seimage  shear  fitpsf

Requires desdb and cx_Oracle
"""
import sys
import deswl
from optparse import OptionParser
parser = OptionParser(__doc__)
options, args = parser.parse_args(sys.argv[1:])

if len(args) < 1:
    parser.print_help()
    sys.exit(45)

merun=args[0]

mi=deswl.files.MultishearSEInputs(merun)
mi.generate_all_inputs()
