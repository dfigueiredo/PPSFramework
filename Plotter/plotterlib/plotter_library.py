#!/usr/bin/env python

from plotterlib import colors
import numpy as numpy

import time
import os
import sys

color = colors.Paint()

class PlotterLibrary():

  def Drawing(self, enable, signal_file, bkg_file, variable_name, title_and_axis_name, bin_size, first_bin, last_bin, signal_cuts, bkg_cuts, tagname):

    print color.BOLD + color.OKBLUE + "[Checking parameters] " + color.ENDC 
    print "\t" + color.BOLD + "Enable: " + color.ENDC,
    print "\t" + color.OKGREEN + enable + color.ENDC

    print "\t" + color.BOLD + "Signal File: " + color.ENDC,
    print "\t" + color.OKGREEN + signal_file + color.ENDC

    print "\t" + color.BOLD + "Background File: " + color.ENDC,
    print "\t" + color.OKGREEN + bkg_file + color.ENDC

    print "\t" + color.BOLD + "Variable Name: " + color.ENDC,
    print "\t" + color.OKGREEN + variable_name + color.ENDC

    print "\t" + color.BOLD + "Title and Axis Name: " + color.ENDC,
    print "\t" + color.OKGREEN + title_and_axis_name + color.ENDC

    print "\t" + color.BOLD + "Bin size: " + color.ENDC,
    print "\t" + color.OKGREEN + str(bin_size) + color.ENDC

    print "\t" + color.BOLD + "First bin: " + color.ENDC,
    print "\t" + color.OKGREEN + str(first_bin) + color.ENDC

    print "\t" + color.BOLD + "Last Bin: " + color.ENDC,
    print "\t" + color.OKGREEN + str(last_bin) + color.ENDC

    print "\t" + color.BOLD + "Signal Cuts: " + color.ENDC,
    print "\t" + color.OKGREEN + str(signal_cuts) + color.ENDC

    print "\t" + color.BOLD + "Bkg Cuts: " + color.ENDC,
    print "\t" + color.OKGREEN + str(bkg_cuts) + color.ENDC

    print "\t" + color.BOLD + "Tag Name: " + color.ENDC,
    print "\t" + color.OKGREEN + str(tagname) + color.ENDC

    if int(enable):
    	print "\t" + color.BOLD + color.HEADER + "-- Drawing enabled --" + color.ENDC
	command = "./Plotter --fsignal " + str(signal_file) + " --fbkg " + str(bkg_file) + " --signalcuts " + str(signal_cuts) + " --bkgcuts " + str(bkg_cuts) + " --title " + str(title_and_axis_name) + " --variable " + str(variable_name) + " --tagname " + str(tagname) + " --first_bin " + str(first_bin) + " --last_bin " + str(last_bin) + " --binsize " + str(bin_size)
	os.system(command);
    else:
    	print "\t" + color.BOLD + color.HEADER + "-- Drawing not enabled --" + color.ENDC
    print "\n"	
