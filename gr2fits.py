from optparse import OptionParser
import sdrutils
    
##################################################
def main():
    usage="Usage: %prog <file>\n"
    parser = OptionParser(usage=usage)
    parser.add_option('-s','--separate',action='store_false',
                      dest='merge',
                      default=True,
                      help='Keep blocks separate?')
    parser.add_option('-v','--verbose',
                  action="store_true",dest="verbose",default=False,
                  help="Increase verbosity of output")

    (options, args) = parser.parse_args()
    for file in args:
        sdrutils.gr2fits(file, options.merge, options.verbose)
        
######################################################################
# Running as executable
if __name__=='__main__':
    main()
