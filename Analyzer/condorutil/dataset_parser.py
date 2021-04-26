#!/usr/bin/env python

from xml.dom import minidom
from condorutil import colors
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
        samples_node = self.mydoc.getElementsByTagName('samples')
        for samples in samples_node:
            dataset_node = samples.getElementsByTagName('dataset')
            for dataset in dataset_node:

		try:
	                samplename = dataset.getElementsByTagName("name")[0]
                	era = dataset.getElementsByTagName("era")[0]
                	mode = dataset.getElementsByTagName("mode")[0]
 			xa = dataset.getElementsByTagName("xa")[0]
                	output = dataset.getElementsByTagName("output")[0]
			datatype = dataset.getElementsByTagName("datatype")[0]
			params = dataset.getElementsByTagName("parameters")[0]
                	enable = dataset.getElementsByTagName("enable")[0]
			command.append([str(samplename.firstChild.data), str(era.firstChild.data), str(mode.firstChild.data), str(xa.firstChild.data), str(output.firstChild.data), str(datatype.firstChild.data), str(params.firstChild.data), str(enable.firstChild.data)])
               		if verbose:
	                   print color.OKBLUE + "\tName: " + samplename.firstChild.data + color.ENDC
        	           print color.OKBLUE + "\tEra: " + era.firstChild.data + color.ENDC
                	   print color.OKBLUE + "\tMode: " + mode.firstChild.data + color.ENDC
			   print color.OKBLUE + "\tX-Angle: " + xa.firstChild.data + color.ENDC
	                   print color.OKBLUE + "\tOutput: " + output.firstChild.data + color.ENDC
			   print color.OKBLUE + "\tDatatype: " + datatype.firstChild.data + color.ENDC
        	           print color.OKBLUE + "\tParamaters: " + params.firstChild.data + color.ENDC
        	           print color.OKBLUE + "\tEnable: " + enable.firstChild.data + color.ENDC
	
		except:
                        print color.FAIL+color.BOLD+'\tFailed to get all the parameters from XML file! Please, check your XML file, there is(are) some error(s)!'+color.ENDC+color.HEADER+color.ENDC
                        exit(0)

        return command

