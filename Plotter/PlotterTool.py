from cmd import Cmd
from plotterlib import colors
from plotterlib import dataset_parser
from plotterlib import plotter_library
from optparse import OptionParser
import os

def getOptions():

    """
    Parse and return the arguments provided by the user.
    """

    usage = "usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-f", "--filename",
                  metavar="FILE", help="XML mapping file", default='plots.xml')
    parser.add_option("-p", "--parsing",
                  help="parsing: commands which can be passed from SHELL directly. [parsing: --p \"submit --file filename.xml\"]")

    parser.add_option("-v", "--verbose",
                  action="store_true", dest="verbose", default=False,
                  help="make lots of noise [default]")

    (options, arguments) = parser.parse_args()
    return options

def execute(fromXML, plotsub):

     try:
                for i in fromXML:

                                #(signal file, bkg file, variable name, cosmetics, bin size, first bin, last bin, signal cuts, bkg cuts, tagname)
                                plotdraw = plotsub.Drawing(i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7], i[8], i[9], i[10])

     except:
                print color.FAIL+color.BOLD+'\tFailed! Please try again!'+color.ENDC+color.HEADER+color.ENDC
                exit(0)

class MyPrompt(Cmd):
    prompt = 'plotter_cmd> '
    intro = "Type ? to list commands"
 
    def do_exit(self, inp):
        print("\n Exiting...\n")
        return True
    
    def help_exit(self):
        print('exit the application. Shorthand: x q Ctrl-D.')
 
    def do_draw(self, arg):

	arcmd = arg.split()
	filename = options.filename

	print filename

	if ("--file" or "--f") in arg: 
                if len(arcmd)>1:
                        for i in arcmd[1:]:
                                filename=str(i)

		if os.path.exists(filename):
	        	xml = dataset_parser.Parser(filename, options.verbose)
	        	fromXML = xml.GetFromXML(options.verbose)
		        plotsub = plotter_library.PlotterLibrary()
			execute(fromXML, plotsub)
		
		else:
        	        print color.FAIL+color.BOLD+'\tXML file does not exist or wrong path! Please use the option --file filename.xml or run the application again with the option --f filename.xml'+color.ENDC+color.HEADER+color.ENDC

	else:
                xml = dataset_parser.Parser(options.filename, options.verbose)
                fromXML = xml.GetFromXML(options.verbose)
                plotsub = plotter_library.PlotterLibrary()
                execute(fromXML, plotsub)			
	
    def help_draw(self):
        print("\n\nUse the options --file filename.xml\n\tEx: submit --file filename.xml <press enter>\n\n")

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

    if options.parsing:
	MyPrompt().onecmd(''.join(options.parsing))
    else:
        MyPrompt().cmdloop()

