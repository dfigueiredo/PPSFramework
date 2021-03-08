##########!/usr/bin/env python

from crabutil import colors
import numpy as numpy

import time
import os
import sys

import CRABClient
from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
#from CRABClient.UserUtilities import config, 
from CRABClient.UserUtilities import config
from httplib import HTTPException

color = colors.Paint()

class CrabLibrary():

  def doSubmit(self, dataset, mode, era, year, json, tier, configfile, filesPerJob, tagname, enable, lfndir, config):

    timestr = time.strftime("%Y-%m-%d_UTC%H-%M-%S")
    pathfull = str(lfndir) + '/' + str(tagname) + '_' + str(timestr)

    print(color.BOLD + color.OKBLUE + "\n[Submitting the following dataset on grid] " + color.ENDC)
    print(color.BOLD +  "\tSample: " + color.ENDC + color.OKGREEN + "\t" + dataset + color.ENDC)
    print(color.BOLD + "\tConfig file: " + color.ENDC + color.OKGREEN + "\t" + configfile + color.ENDC)
    print(color.BOLD + "\tTagName: " + color.ENDC + color.OKGREEN + "\t" + tagname + color.ENDC)
    print(color.BOLD + "\tEra: " + color.ENDC + color.OKGREEN + "\t" + era + color.ENDC)
    print(color.BOLD + "\tJSON: " + color.ENDC + color.OKGREEN + "\t" + json + color.ENDC)
    print(color.BOLD + "\tMode: " + color.ENDC + color.OKGREEN + "\t" + mode + color.ENDC)
    print(color.BOLD + "\tYear: " + color.ENDC + color.OKGREEN + "\t" + str(year) + color.ENDC)
    print(color.BOLD + "\tTier: " + color.ENDC + color.OKGREEN + "\t" + str(tier) + color.ENDC)
    print(color.BOLD + "\tFiles per Job: " + color.ENDC + color.OKGREEN + "\t" + str(filesPerJob) + color.ENDC)
    print(color.BOLD + "\tEnable: " + color.ENDC + color.OKGREEN + "\t" + str(enable) + color.ENDC)
    print(color.BOLD + "\tLFN output dir: " + color.ENDC + color.OKGREEN + "\t" + str(pathfull) + color.ENDC)

    timestr = time.strftime("%Y-%m-%d_UTC%H-%M-%S")

    config.General.transferLogs = True
    config.General.transferOutputs = True
    config.JobType.pluginName = 'Analysis'
    config.Data.inputDBS = 'global' #phys03

    config.Data.splitting = 'LumiBased' # or 'LumiBased', ''FileBased','Automatic'
    #config.Data.unitsPerJob = int(filesPerJob)
    config.Data.unitsPerJob = 100
    config.Data.lumiMask = json # need to be removed when it is a MC
    
    #NJOBS = 1000
    #config.Data.totalUnits = 600000

    config.Data.publication = False
    config.JobType.psetName = configfile
    config.JobType.outputFiles = ['output.root']
    config.General.workArea = 'crab_ProtonReco_' + str(tagname) + '_' + str(timestr)
    config.General.requestName = str(tagname)
    config.Data.inputDataset = str(dataset)
    ##config.JobType.allowUndistributedCMSSW = True
    #config.JobType.maxMemoryMB = 2500
    config.Data.outputDatasetTag = str(tagname)
    config.Site.storageSite = 'T2_CH_CERN' #T2_CH_CERNBOX
    config.Data.outLFNDirBase = pathfull
    config.JobType.pyCfgParams = ["Mode="+str(mode),"Year="+str(year)]

    if int(enable):
        res = crabCommand('submit', config = config)
    else:
        print("\t" + color.BOLD + color.HEADER + "-- Submittion not enabled --" + color.ENDC)
    print("\n")	
