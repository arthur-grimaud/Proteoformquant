#!/usr/bin/env python
import 

def doArgs(argList, name):
    parser = argparse.ArgumentParser(description=name)

    parser.add_argument('-v', "--verbose", action="store_true", help="Enable verbose debugging", default=False)
    parser.add_argument('--input', action="store", dest="inputFn", type=str, help="Input file name", required=True)
    parser.add_argument('--output', action="store", dest="outputFn", type=str, help="Output file name", required=True)

    return parser.parse_args(argList)

def main():
    progName = "Template"
    args = doArgs(sys.argv[1:], progName)

    verbose = args.verbose
    inputFn = args.inputFn
    outputFn = args.outputFn

    print "Starting %s" % (progName)
    startTime = float(time.time())

    if not os.path.isfile(inputFn):
        print "Input doesn't exist, exiting"
        return

    outputBase = os.path.dirname(outputFn)
    if outputBase!='' and not os.path.exists(outputBase):
        print "Output directory doesn't exist, making output dirs: %s" % (outputBase)
        os.makedirs(outputBase)


    print "Finished in %0.4f seconds" % (time.time() - startTime)
    return

if __name__ == '__main__':
    #sys.argv = ["programName.py","--input","test.txt","--output","tmp/test.txt"]
    main()
