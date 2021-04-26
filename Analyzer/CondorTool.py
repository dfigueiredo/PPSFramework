from cmd import Cmd
from condorutil import colors
from condorutil import dataset_parser
from condorutil import condor_library
from optparse import OptionParser
import os

def getOptions():

    """
    Parse and return the arguments provided by the user.
    """

    usage = "usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-f", "--filename",
                  metavar="FILE", help="XML mapping file", default='condor_sd_elec.xml')
    parser.add_option("-p", "--parsing",
                  help="parsing: commands which can be passed from SHELL directly. [parsing: --p \"submit --file filename.xml\"]")

    parser.add_option("-v", "--verbose",
                  action="store_true", dest="verbose", default=False,
                  help="make lots of noise [default]")

    (options, arguments) = parser.parse_args()
    return options

def execute(fromXML, condorsub):

     try:
                for i in fromXML:
				# dataset, era, mode, xangle, datatype, output, params, enable
                                condorsample = condorsub.doSubmit(i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7])

     except:
                print color.FAIL+color.BOLD+'\tFailed! Please try again!'+color.ENDC+color.HEADER+color.ENDC
                exit(0)

class MyPrompt(Cmd):
    prompt = 'condor_submission> '
    intro = "Type ? to list commands"
 
    def do_exit(self, inp):
        print("\n Exiting...\n")
        return True
    
    def help_exit(self):
        print('exit the application. Shorthand: x q Ctrl-D.')
 
    def do_submit(self, arg):

	os.system("voms-proxy-init --voms cms")

	arcmd = arg.split()
	filename = options.filename

	if ("--file" or "--f") in arg: 
                if len(arcmd)>1:
                        for i in arcmd[1:]:
                                filename=str(i)

		if os.path.exists(filename):
	        	xml = dataset_parser.Parser(filename, options.verbose)
	        	fromXML = xml.GetFromXML(options.verbose)
		        condorsub = condor_library.CondorLibrary()
			execute(fromXML, condorsub)
		
		else:
        	        print color.FAIL+color.BOLD+'\tXML file does not exist or wrong path! Please use the option --file filename.xml or run the application again with the option --f filename.xml'+color.ENDC+color.HEADER+color.ENDC

	else:
                xml = dataset_parser.Parser(options.filename, options.verbose)
                fromXML = xml.GetFromXML(options.verbose)
                condorsub = condor_library.CondorLibrary()
                execute(fromXML, condorsub)			


    def help_submit(self):
        print("\n\nUse the options --file filename.xml\n\tEx: submit --file filename.xml <press enter>\n\n")

    def do_status(self, inp):
	os.system("_CONDOR_SCHEDD_HOST=bigbird17.cern.ch _CONDOR_CREDD_HOST=bigbird17.cern.ch condor_q")

    def help_status(self):
        print("\n\nCheck the condor jobs status.\n\n")

    def do_kill(self, inp):
	os.system("_CONDOR_SCHEDD_HOST=bigbird17.cern.ch _CONDOR_CREDD_HOST=bigbird17.cern.ch condor_rm -all")

    def help_kill(self):
        print("\n\nKill all the condor jobs!\n\n")
    
    def do_holdreason(self, inp):
        os.system("_CONDOR_SCHEDD_HOST=bigbird17.cern.ch _CONDOR_CREDD_HOST=bigbird17.cern.ch condor_q -hold -af HoldReason")

    def help_holdreason(self):
        print("\n\nRetrieve the reason why the job was put on hold state\n\n")

    def default(self, inp):
        if inp == 'x' or inp == 'q':
            return self.do_exit(inp)
 
        print("Default: {}".format(inp))

    def emptyline(self):
        pass

    do_EOF = do_exit
    help_EOF = help_exit
 
if __name__ == '__main__':

    color = colors.Paint()
    options = getOptions()

    # Creating Output Folder
    # submit from (CTPPSAnalysisCode/Analyzer) folder
    path = os.getcwd()

    # Checking if cmsenv has been done
    os.system("cd "+path)
    try:
        ReleaseBase = os.path.join(os.environ['CMSSW_BASE'], "src")
        ReleaseVersion = os.environ['CMSSW_VERSION']
    except KeyError:
        print "CMSSW enviroment not set, please run cmsenv!"
        sys.exit()

    if options.parsing:
	MyPrompt().onecmd(''.join(options.parsing))
    else:
        MyPrompt().cmdloop()
