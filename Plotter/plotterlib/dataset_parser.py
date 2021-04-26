#!/usr/bin/env python

from xml.dom import minidom
from plotterlib import colors
import numpy as numpy
import sys

color = colors.Paint()

class Parser():

    def __init__(self, filename, verbose):

        # parse an xml file by name
        self.mydoc = minidom.parse(filename)
        cNodes = self.mydoc.childNodes

	if verbose:
	        print "\nReading file..."
        	print color.BOLD + cNodes[0].toxml() + color.ENDC
	        print "\n"

    def GetFromXML(self, verbose):

        command = []
        samples_node = self.mydoc.getElementsByTagName('plots')
        for samples in samples_node:
            dataset_node = samples.getElementsByTagName('drawing')
            for dataset in dataset_node:

		#try:
                	enable = dataset.getElementsByTagName("enable")[0]
	                signal = dataset.getElementsByTagName("signal_file")[0]
	                bkg = dataset.getElementsByTagName("bkg_file")[0]
	                variable = dataset.getElementsByTagName("variable_name")[0]
	                cosmetics = dataset.getElementsByTagName("title_and_axis_name")[0]
	                binsize = dataset.getElementsByTagName("bin_size")[0]
	                firstbin = dataset.getElementsByTagName("first_bin")[0]
	                lastbin = dataset.getElementsByTagName("last_bin")[0]
	                signal_cut = dataset.getElementsByTagName("signal_cuts")[0]
	                bkg_cut = dataset.getElementsByTagName("bkg_cuts")[0]
	                tagname = dataset.getElementsByTagName("tagname")[0]
			command.append([str(enable.firstChild.data), str(signal.firstChild.data), str(bkg.firstChild.data), str(variable.firstChild.data),
					str(cosmetics.firstChild.data), str(binsize.firstChild.data), str(firstbin.firstChild.data), str(lastbin.firstChild.data),
					str(signal_cut.firstChild.data), str(bkg_cut.firstChild.data), str(tagname.firstChild.data)])
               		if verbose:
                	   print color.OKBLUE + "\tEnable: " + enable.firstChild.data + color.ENDC
	                   print color.OKBLUE + "\tSignal File: " + signal.firstChild.data + color.ENDC
	                   print color.OKBLUE + "\tBackground File: " + bkg.firstChild.data + color.ENDC
                	   print color.OKBLUE + "\tVariable: " + variable.firstChild.data + color.ENDC
                	   print color.OKBLUE + "\tCosmetics: " + cosmetics.firstChild.data + color.ENDC
                	   print color.OKBLUE + "\tBin Size: " + binsize.firstChild.data + color.ENDC
                	   print color.OKBLUE + "\tFirst Bin: " + firstbin.firstChild.data + color.ENDC
                	   print color.OKBLUE + "\tLast Bin: " + lastbin.firstChild.data + color.ENDC
                	   print color.OKBLUE + "\tSignal Cut: " + signal_cut.firstChild.data + color.ENDC
                	   print color.OKBLUE + "\tBkg Cut: " + bkg_cut.firstChild.data + color.ENDC
                	   print color.OKBLUE + "\tTagname: " + tagname.firstChild.data + color.ENDC
		
		#except:
                #        print color.FAIL+color.BOLD+'\tFailed to get all the parameters from XML file! Please, check your XML file, there is(are) some error(s)!'+color.ENDC+color.HEADER+color.ENDC
                #        exit(0)

        return command

