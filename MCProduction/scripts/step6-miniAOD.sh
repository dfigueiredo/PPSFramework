#!/bin/bash

#############################################################
#   This script is used by McM when it performs automatic   #
#  validation in HTCondor or submits requests to computing  #
#                                                           #
#      !!! THIS FILE IS NOT MEANT TO BE RUN BY YOU !!!      #
# If you want to run validation script yourself you need to #
#     get a "Get test" script which can be retrieved by     #
#  clicking a button next to one you just clicked. It will  #
# say "Get test command" when you hover your mouse over it  #
#      If you try to run this, you will have a bad time     #
#############################################################

cd /afs/cern.ch/cms/PPD/PdmV/work/McM/submit/TOP-RunIISummer20UL17MiniAOD-00006/

# Make voms proxy
voms-proxy-init --voms cms --out $(pwd)/voms_proxy.txt --hours 4
export X509_USER_PROXY=$(pwd)/voms_proxy.txt

export SCRAM_ARCH=slc7_amd64_gcc700

source /cvmfs/cms.cern.ch/cmsset_default.sh
if [ -r CMSSW_10_6_17_patch1/src ] ; then
  echo release CMSSW_10_6_17_patch1 already exists
else
  scram p CMSSW CMSSW_10_6_17_patch1
fi
cd CMSSW_10_6_17_patch1/src
eval `scram runtime -sh`

scram b
cd ../..

# Maximum validation duration: 28800s
# Margin for validation duration: 20%
# Validation duration with margin: 28800 * (1 - 0.20) = 23040s
# Time per event for each sequence: 0.3391s
# Threads for each sequence: 8
# Time per event for single thread for each sequence: 8 * 0.3391s = 2.7128s
# Which adds up to 2.7128s per event
# Single core events that fit in validation duration: 23040s / 2.7128s = 8493
# Produced events limit in McM is 10000
# According to 1.0000 efficiency, up to 10000 / 1.0000 = 10000 events should run
# Clamp (put value) 8493 within 1 and 10000 -> 8493
# It is estimated that this validation will produce: 8493 * 1.0000 = 8493 events
EVENTS=8493


# cmsDriver command
cmsDriver.py  --python_filename TOP-RunIISummer20UL17MiniAOD-00006_1_cfg.py --eventcontent MINIAODSIM --customise Configuration/DataProcessing/Utils.addMonitoring --datatier MINIAODSIM --fileout file:TOP-RunIISummer20UL17MiniAOD-00006.root --conditions 106X_mc2017_realistic_v6 --step PAT --geometry DB:Extended --filein "dbs:/ExclusiveTTbar_LepHad_TuneCP5_13TeV-FPMC-madspin-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v1/AODSIM" --era Run2_2017 --runUnscheduled --no_exec --mc -n $EVENTS || exit $? ;

