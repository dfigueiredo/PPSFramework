#!/usr/bin/env python

from condorutil import colors
import numpy as numpy

import time
import os
import sys

color = colors.Paint()

class CondorLibrary():

  def doSubmit(self, name, era, mode, xa, output, datatype, params, enable):

    print color.BOLD + color.OKBLUE + "[Submitting the following ntuples on condor] " + color.ENDC 
    print "\t" + color.BOLD + "Name: " + color.ENDC,
    print "\t" + color.OKGREEN + name + color.ENDC

    print "\t" + color.BOLD + "Era: " + color.ENDC,
    print "\t" + color.OKGREEN + era + color.ENDC

    print "\t" + color.BOLD + "Mode: " + color.ENDC,
    print "\t" + color.OKGREEN + mode + color.ENDC

    print "\t" + color.BOLD + "X-Angle: " + color.ENDC,
    print "\t" + color.OKGREEN + xa + color.ENDC

    print "\t" + color.BOLD + "Output Folder: " + color.ENDC,
    print "\t" + color.OKGREEN + output + color.ENDC

    print "\t" + color.BOLD + "Datatype: " + color.ENDC,
    print "\t" + color.OKGREEN + datatype + color.ENDC

    print "\t" + color.BOLD + "Parameters: " + color.ENDC,
    print "\t" + color.OKGREEN + params + color.ENDC

    print "\t" + color.BOLD + "Enable: " + color.ENDC,
    print "\t" + color.OKGREEN + enable + color.ENDC

    if int(enable):	
	    print "\t" + color.BOLD + color.HEADER + "-- Submittion enabled --" + color.ENDC

    	# Creating Output Folder
	    path = os.getcwd()
	    folderout = str(path)+"/"+output
	    #folderout = output

	    if not os.path.exists(folderout):
        	os.makedirs(folderout)
	        os.system("cp proton*.txt "+folderout+"/.")
	    else:
	        os.system("rm -rf "+folderout)

	    files_input_condor = "transfer_input_files\t\t\t= ../proton_reco_rphorizontal_Electron_era_B_xa_120.txt, ../proton_reco_rphorizontal_Electron_era_B_xa_130.txt, ../proton_reco_rphorizontal_Electron_era_B_xa_140.txt, ../proton_reco_rphorizontal_Electron_era_B_xa_150.txt, ../proton_reco_rphorizontal_Electron_era_C_xa_120.txt, ../proton_reco_rphorizontal_Electron_era_C_xa_130.txt, ../proton_reco_rphorizontal_Electron_era_C_xa_140.txt, ../proton_reco_rphorizontal_Electron_era_C_xa_150.txt, ../proton_reco_rphorizontal_Electron_era_D_xa_120.txt, ../proton_reco_rphorizontal_Electron_era_D_xa_130.txt, ../proton_reco_rphorizontal_Electron_era_D_xa_140.txt, ../proton_reco_rphorizontal_Electron_era_D_xa_150.txt, ../proton_reco_rphorizontal_Electron_era_E_xa_120.txt, ../proton_reco_rphorizontal_Electron_era_E_xa_130.txt, ../proton_reco_rphorizontal_Electron_era_E_xa_140.txt, ../proton_reco_rphorizontal_Electron_era_E_xa_150.txt, ../proton_reco_rphorizontal_Electron_era_F_xa_120.txt, ../proton_reco_rphorizontal_Electron_era_F_xa_130.txt, ../proton_reco_rphorizontal_Electron_era_F_xa_140.txt, ../proton_reco_rphorizontal_Electron_era_F_xa_150.txt, ../proton_reco_rphorizontal_Electron_era_preTS2_xa_120.txt, ../proton_reco_rphorizontal_Electron_era_preTS2_xa_130.txt, ../proton_reco_rphorizontal_Electron_era_preTS2_xa_140.txt, ../proton_reco_rphorizontal_Electron_era_preTS2_xa_150.txt, ../proton_reco_rphorizontal_Electron_era_postTS2_xa_120.txt, ../proton_reco_rphorizontal_Electron_era_postTS2_xa_130.txt, ../proton_reco_rphorizontal_Electron_era_postTS2_xa_140.txt, ../proton_reco_rphorizontal_Electron_era_postTS2_xa_150.txt, ../proton_reco_rphorizontal_Muon_era_B_xa_120.txt, ../proton_reco_rphorizontal_Muon_era_B_xa_130.txt, ../proton_reco_rphorizontal_Muon_era_B_xa_140.txt, ../proton_reco_rphorizontal_Muon_era_B_xa_150.txt, ../proton_reco_rphorizontal_Muon_era_C_xa_120.txt, ../proton_reco_rphorizontal_Muon_era_C_xa_130.txt, ../proton_reco_rphorizontal_Muon_era_C_xa_140.txt, ../proton_reco_rphorizontal_Muon_era_C_xa_150.txt, ../proton_reco_rphorizontal_Muon_era_D_xa_120.txt, ../proton_reco_rphorizontal_Muon_era_D_xa_130.txt, ../proton_reco_rphorizontal_Muon_era_D_xa_140.txt, ../proton_reco_rphorizontal_Muon_era_D_xa_150.txt, ../proton_reco_rphorizontal_Muon_era_E_xa_120.txt, ../proton_reco_rphorizontal_Muon_era_E_xa_130.txt, ../proton_reco_rphorizontal_Muon_era_E_xa_140.txt, ../proton_reco_rphorizontal_Muon_era_E_xa_150.txt, ../proton_reco_rphorizontal_Muon_era_F_xa_120.txt, ../proton_reco_rphorizontal_Muon_era_F_xa_130.txt, ../proton_reco_rphorizontal_Muon_era_F_xa_140.txt, ../proton_reco_rphorizontal_Muon_era_F_xa_150.txt, ../proton_reco_rphorizontal_Muon_era_preTS2_xa_120.txt, ../proton_reco_rphorizontal_Muon_era_preTS2_xa_130.txt, ../proton_reco_rphorizontal_Muon_era_preTS2_xa_140.txt, ../proton_reco_rphorizontal_Muon_era_preTS2_xa_150.txt, ../proton_reco_rphorizontal_Muon_era_postTS2_xa_120.txt, ../proton_reco_rphorizontal_Muon_era_postTS2_xa_130.txt, ../proton_reco_rphorizontal_Muon_era_postTS2_xa_140.txt, ../proton_reco_rphorizontal_Muon_era_postTS2_xa_150.txt, ../PreliminaryEfficiencies_October92019_1D2DMultiTrack.root\n"

	#../proton_reco_rphorizontal_Electron_era_B_xa_120.txt, ../proton_reco_rphorizontal_Electron_era_B_xa_130.txt, ../proton_reco_rphorizontal_Electron_era_B_xa_140.txt, ../proton_reco_rphorizontal_Electron_era_B_xa_150.txt, ../proton_reco_rphorizontal_Electron_era_C_xa_120.txt, ../proton_reco_rphorizontal_Electron_era_C_xa_130.txt, ../proton_reco_rphorizontal_Electron_era_C_xa_140.txt, ../proton_reco_rphorizontal_Electron_era_C_xa_150.txt, ../proton_reco_rphorizontal_Electron_era_D_xa_120.txt, ../proton_reco_rphorizontal_Electron_era_D_xa_130.txt, ../proton_reco_rphorizontal_Electron_era_D_xa_140.txt, ../proton_reco_rphorizontal_Electron_era_D_xa_150.txt, ../proton_reco_rphorizontal_Electron_era_E_xa_120.txt, ../proton_reco_rphorizontal_Electron_era_E_xa_130.txt, ../proton_reco_rphorizontal_Electron_era_E_xa_140.txt, ../proton_reco_rphorizontal_Electron_era_E_xa_150.txt, ../proton_reco_rphorizontal_Electron_era_F_xa_120.txt, ../proton_reco_rphorizontal_Electron_era_F_xa_130.txt, ../proton_reco_rphorizontal_Electron_era_F_xa_140.txt, ../proton_reco_rphorizontal_Electron_era_F_xa_150.txt, ../proton_reco_rphorizontal_Electron_era_preTS2_xa_120.txt, ../proton_reco_rphorizontal_Electron_era_preTS2_xa_130.txt, ../proton_reco_rphorizontal_Electron_era_preTS2_xa_140.txt, ../proton_reco_rphorizontal_Electron_era_preTS2_xa_150.txt, ../proton_reco_rphorizontal_Electron_era_postTS2_xa_120.txt, ../proton_reco_rphorizontal_Electron_era_postTS2_xa_130.txt, ../proton_reco_rphorizontal_Electron_era_postTS2_xa_140.txt, ../proton_reco_rphorizontal_Electron_era_postTS2_xa_150.txt, ../proton_reco_rphorizontal_Muon_era_B_xa_120.txt, ../proton_reco_rphorizontal_Muon_era_B_xa_130.txt, ../proton_reco_rphorizontal_Muon_era_B_xa_140.txt, ../proton_reco_rphorizontal_Muon_era_B_xa_150.txt, ../proton_reco_rphorizontal_Muon_era_C_xa_120.txt, ../proton_reco_rphorizontal_Muon_era_C_xa_130.txt, ../proton_reco_rphorizontal_Muon_era_C_xa_140.txt, ../proton_reco_rphorizontal_Muon_era_C_xa_150.txt, ../proton_reco_rphorizontal_Muon_era_D_xa_120.txt, ../proton_reco_rphorizontal_Muon_era_D_xa_130.txt, ../proton_reco_rphorizontal_Muon_era_D_xa_140.txt, ../proton_reco_rphorizontal_Muon_era_D_xa_150.txt, ../proton_reco_rphorizontal_Muon_era_E_xa_120.txt, ../proton_reco_rphorizontal_Muon_era_E_xa_130.txt, ../proton_reco_rphorizontal_Muon_era_E_xa_140.txt, ../proton_reco_rphorizontal_Muon_era_E_xa_150.txt, ../proton_reco_rphorizontal_Muon_era_F_xa_120.txt, ../proton_reco_rphorizontal_Muon_era_F_xa_130.txt, ../proton_reco_rphorizontal_Muon_era_F_xa_140.txt, ../proton_reco_rphorizontal_Muon_era_F_xa_150.txt, ../proton_reco_rphorizontal_Muon_era_preTS2_xa_120.txt, ../proton_reco_rphorizontal_Muon_era_preTS2_xa_130.txt, ../proton_reco_rphorizontal_Muon_era_preTS2_xa_140.txt, ../proton_reco_rphorizontal_Muon_era_preTS2_xa_150.txt, ../proton_reco_rphorizontal_Muon_era_postTS2_xa_120.txt, ../proton_reco_rphorizontal_Muon_era_postTS2_xa_130.txt, ../proton_reco_rphorizontal_Muon_era_postTS2_xa_140.txt, ../proton_reco_rphorizontal_Muon_era_postTS2_xa_150.txt,

	    with open('job_condor_tmp.sub', 'w') as fout:
        	fout.write("initialdir\t\t\t= "+output+"\n")
	        fout.write("executable\t\t\t= ./MissingMassNtupleAnalyzer\n")
        	fout.write("arguments\t\t\t= --f $(filename) --era "+era+" --mode "+mode+" --xa "+xa+" --datatype "+datatype+" --jobid $(ProcId) "+params+"\n")
	        fout.write(files_input_condor)
	        fout.write("output\t\t\t= missing_mass.$(ClusterId).$(ProcId).out\n")
	        fout.write("error\t\t\t= missing_mass.$(ClusterId).$(ProcId).err\n")
	        fout.write("log\t\t\t= missing_mass.$(ClusterId).$(ProcId).log\n")
	        fout.write("getenv\t\t\t= True\n\n")
                #fout.write("Should_Transfer_Files = YES\n")
                #fout.write("output_destination = /eos/cms/store/group/phys_exotica/PPS-Exo/Condor/\n")
                #fout.write("WhenToTransferOutput = ON_EXIT\n")

        	###########################
	        # espresso     = 20 minutes
        	# microcentury = 1 hour
	        # longlunch    = 2 hours
	        # workday      = 8 hours
	        # tomorrow     = 1 day
	        # testmatch    = 3 days
	        # nextweek     = 1 week
	        ##########################

        	fout.write("+JobFlavour\t\t\t= \"espresso\"\n")
        	#fout.write("+JobFlavour\t\t\t= \"workday\"\n")
	        #fout.write("RequestCpus\t\t\t = 4\n")
	        #fout.write("+MaxRuntime\t\t\t = 7200\n")
		#fout.write("max_transfer_input_mb = 2048\n")
		#fout.write("max_transfer_output_mb = 2048\n")
        	fout.write("requirements\t\t\t = (OpSysAndVer =?= \"CentOS7\")\n\n")
	        fout.write("queue filename matching("+name+"*.root)\n")

	    os.system('_CONDOR_SCHEDD_HOST=bigbird17.cern.ch _CONDOR_CREDD_HOST=bigbird17.cern.ch condor_submit job_condor_tmp.sub')
	    os.system('rm job_condor_tmp.sub')
	
    else:
        print "\t" + color.BOLD + color.HEADER + "-- Submittion not enabled --" + color.ENDC
 
    print "\n"


