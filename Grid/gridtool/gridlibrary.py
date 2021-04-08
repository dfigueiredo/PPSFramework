#!/usr/bin/env python

#!interpreter [optional-arg]
# -*- coding: utf-8 -*-

""" Module used to submit crab jobs on grid.
This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.
This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with
this program. If not, see <http://www.gnu.org/licenses/>.
"""
__author__ = "Diego Figueiredo"
__contact__ = "dmf@cern.ch"
__copyright__ = "Copyright 2021, INFN-Pisa"
__credits__ = ["Diego Figueiredo"]
__date__ = "2021/04/05"
__deprecated__ = False
__email__ =  "dmf@cern.ch"
__license__ = "GPLv3"
__maintainer__ = "developer"
__status__ = "Validation"
__version__ = "0.0.1"

import re
import json
from xml.dom import minidom
from gridtool import colors
import numpy as numpy
import sys
import time
import os
import getpass

from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
from CRABClient.UserUtilities import config
from httplib import HTTPException
from multiprocessing import Process

color = colors.Paint()
NBLOCKS = 12 # number of keys per dataset block.

class Parser():

    def __init__(self, filename, verbose):

        # Count number of JSON inputs (datasets)
        self.count_data = 0
        self.count_key = 0
        self.parameters = []        
 
        # parse a JSON file by name
        with open(filename) as json_file:
         try:
          self.data = json.load(json_file)
         except ValueError as err:
          print("\n"+color.FAIL+"[gridtool] The file {} is not valid.\nPlease, check it."+color.ENDC+"\n").format(json_file)
          exit()

        if "datasets" in self.data:
         for p in self.data["datasets"]:
          self.count_data = self.count_data + 1
          if "id" not in p:
            print("\n"+color.FAIL+"[gridtool] The key \"id\" must be defined in each dataset bracket.\nPlease, check it."+color.ENDC+"\n")
            exit()
          else:
            check_id = isinstance(p["id"], int)
            if not check_id:
             print("\n"+color.FAIL+"[gridtool] The key \"id\" must be integer. Id is defined as {}.\nPlease, check it."+color.ENDC+"\n").format(type(p["id"]))
             exit()
          for key in p:
           if "id" == key or\
              "enable" == key or\
              "localpath" == key or\
              "eospath" == key or\
              "sample" == key or\
              "mode" == key or\
              "lumimask" == key or\
              "config" == key or\
              "parameters" == key or\
              "output" == key or\
              "events" == key or\
              "site" == key:
            self.count_key = self.count_key + 1
           else:
            print("\n"+color.FAIL+"[gridtool] {} is _not_ a valid key for the dataset id {}.\nPlease, check your JSON file."+color.ENDC+"\n").format(str(key), str(self.count_key))
            exit()
        else:
            print("\n"+color.FAIL+"[gridtool] JSON file is not correct.\nPlease, check it."+color.ENDC+"\n")
            exit()

        check = int(self.count_key)/int(self.count_data)
        if not (check==NBLOCKS):
            error_message = '\n'+color.FAIL+'[gridtool] Defined {} keys instead of a total of {} keys in your JSON file.\n'\
                            'In case you are sure the file has the exactly number of parameters you need, you should re-define NBLOCKS in the code accordinly.'\
                            '\nPlease, correct it instead.'+color.ENDC+'\n'
            print(error_message).format(check, NBLOCKS)
            exit()
 
        if verbose:
         print "\n[gridtool] Reading file...\n"
         print(color.BOLD + json.dumps(self.data, indent=4, sort_keys=True) + color.ENDC)
         print "\n"

    def submit(self, config):
        try:
            crabCommand('submit', config = config)
        except HTTPException as hte:
            print "Failed submitting task: %s" % (hte.headers)
        except ClientException as cle:
            print "Failed submitting task: %s" % (cle)

    def prepareSubmission(self, config):

     if "datasets" in self.data:
      for p in self.data["datasets"]:

       timestr = time.strftime("%Y-%m-%d_UTC%H-%M-%S")
       pathname = "PATHNAME"
       pathfull = '/store/user/dmf/%s_%s/' % (pathname, timestr)

       print "\n\t" + color.BOLD + "id: " + color.ENDC,
       print "\t" + color.OKGREEN + str(p["id"]) + color.ENDC

       print "\t" + color.BOLD + "enable: " + color.ENDC,
       print "\t" + color.OKGREEN + str(p["enable"]) + color.ENDC

       print "\t" + color.BOLD + "localpath: " + color.ENDC,
       print "\t" + color.OKGREEN + str(p["localpath"]) + color.ENDC
 
       print "\t" + color.BOLD + "eospath: " + color.ENDC,
       print "\t" + color.OKGREEN + str(p["eospath"]) + color.ENDC

       print "\t" + color.BOLD + "sample: " + color.ENDC,
       print "\t" + color.OKGREEN + str(p["sample"]) + color.ENDC

       print "\t" + color.BOLD + "mode: " + color.ENDC,
       print "\t" + color.OKGREEN + str(p["mode"]) + color.ENDC

       print "\t" + color.BOLD + "lumimask: " + color.ENDC,
       print "\t" + color.OKGREEN + str(p["lumimask"]) + color.ENDC

       print "\t" + color.BOLD + "config: " + color.ENDC,
       print "\t" + color.OKGREEN + str(p["config"]) + color.ENDC

       print "\t" + color.BOLD + "parameters: " + color.ENDC,
       print "\t" + color.OKGREEN + str(p["parameters"]) + color.ENDC

       print "\t" + color.BOLD + "output: " + color.ENDC,
       print "\t" + color.OKGREEN + str(p["output"]) + color.ENDC

       print "\t" + color.BOLD + "events: " + color.ENDC,
       print "\t" + color.OKGREEN + str(p["events"]) + color.ENDC

       print "\t" + color.BOLD + "site: " + color.ENDC,
       print "\t" + color.OKGREEN + str(p["site"]) + color.ENDC

       tagname = 'crab_%s_%s' % (getpass.getuser(), timestr)
       localpath = '/%s/%s/' % (p["localpath"], tagname)
       eospath = '/%s/%s/' % (p["eospath"], tagname)

       # Crab common paratemers 
       config.JobType.psetName = p["config"]
       config.JobType.pyCfgParams = [p["parameters"]] # check how to make it
       config.JobType.outputFiles = [p["output"]]
       config.Data.outLFNDirBase = eospath
       config.General.transferOutputs = True
       config.General.requestName = tagname
       config.General.workArea = localpath
       config.Site.storageSite = p["site"]

       # Crab parameters for each submission case
       if p["mode"] == "data_analysis":
        config.JobType.pluginName = 'Analysis'
        config.Data.inputDataset = p["sample"]
        config.Data.inputDBS = 'global'
        config.Data.splitting = 'LumiBased'
        config.Data.unitsPerJob = 20
        config.Data.publication = False
        config.Data.lumiMask = p["lumimask"]
        config.General.transferLogs = True
       elif p["mode"] == "mc_analysis":
        config.JobType.pluginName = 'Analysis'
        config.Data.inputDataset = p["sample"]
        config.Data.inputDBS = 'global'
        config.Data.splitting = 'LumiBased'
        config.Data.unitsPerJob = 20
        config.Data.publication = False
        config.Data.lumiMask = p["lumimask"]
        config.General.transferLogs = True
       elif p["mode"] == "mc_private_production":
        config.JobType.pluginName = 'PrivateMC'
        config.Data.outputPrimaryDataset = p["sample"]
        config.Data.inputDBS = 'phys03'
        config.Data.splitting = 'EventBased'
        config.Data.unitsPerJob = 500
        config.Data.totalUnits = 6000 * config.Data.unitsPerJob
        config.Data.publication = True
        config.General.transferLogs = False
       else:
        print("\n"+color.WARNING+"[gridtool] mode {} is not defined.\nJobs have not been configured for the task id {}."+color.ENDC+"\n").format(p["mode"], p["id"])
        continue

       if int(p["enable"]):
        print("\n")
        p = Process(target=self.submit, args=(config,))
        p.start()
        p.join()
       else:
        print "\t" + color.BOLD + color.HEADER + "-- Submittion not enabled --" + color.ENDC
       print("\n")
