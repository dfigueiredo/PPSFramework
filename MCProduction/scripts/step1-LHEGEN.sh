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

cd /afs/cern.ch/cms/PPD/PdmV/work/McM/submit/TOP-RunIISummer20UL17pLHEGEN-00002/

# Make voms proxy
voms-proxy-init --voms cms --out $(pwd)/voms_proxy.txt --hours 4
export X509_USER_PROXY=$(pwd)/voms_proxy.txt

export SCRAM_ARCH=slc7_amd64_gcc700

source /cvmfs/cms.cern.ch/cmsset_default.sh
if [ -r CMSSW_10_6_18/src ] ; then
  echo release CMSSW_10_6_18 already exists
else
  scram p CMSSW CMSSW_10_6_18
fi
cd CMSSW_10_6_18/src
eval `scram runtime -sh`

# Download fragment from McM
curl -s -k https://cms-pdmv.cern.ch/mcm/public/restapi/requests/get_fragment/TOP-RunIISummer20UL17pLHEGEN-00002 --retry 3 --create-dirs -o Configuration/GenProduction/python/TOP-RunIISummer20UL17pLHEGEN-00002-fragment.py
[ -s Configuration/GenProduction/python/TOP-RunIISummer20UL17pLHEGEN-00002-fragment.py ] || exit $?;
scram b
cd ../..

# Maximum validation duration: 57600s
# Margin for validation duration: 20%
# Validation duration with margin: 57600 * (1 - 0.20) = 46080s
# Time per event for each sequence: 0.0229s, 0.0229s
# Threads for each sequence: 1, 1
# Time per event for single thread for each sequence: 1 * 0.0229s = 0.0229s, 1 * 0.0229s = 0.0229s
# Which adds up to 0.0458s per event
# Single core events that fit in validation duration: 46080s / 0.0458s = 1006113
# Produced events limit in McM is 10000
# According to 1.0000 efficiency, up to 10000 / 1.0000 = 10000 events should run
# Clamp (put value) 1006113 within 1 and 10000 -> 10000
# It is estimated that this validation will produce: 10000 * 1.0000 = 10000 events
EVENTS=1000


# cmsDriver command
#cmsDriver.py Configuration/GenProduction/python/TOP-RunIISummer20UL17pLHEGEN-00002-fragment.py --python_filename TOP-RunIISummer20UL17pLHEGEN-00002_1_cfg.py --eventcontent LHE --customise Configuration/DataProcessing/Utils.addMonitoring --datatier LHE --fileout file:TOP-RunIISummer20UL17pLHEGEN-00002_0.root --conditions 106X_mc2017_realistic_v6 --step NONE --filein "lhe:19332" --era Run2_2017 --no_exec --mc -n $EVENTS || exit $? ;

# cmsDriver command
#cmsDriver.py Configuration/GenProduction/python/TOP-RunIISummer20UL17pLHEGEN-00002-fragment.py --python_filename TOP-RunIISummer20UL17pLHEGEN-00002_2_cfg.py --eventcontent RAWSIM --customise Configuration/DataProcessing/Utils.addMonitoring --datatier GEN --fileout file:TOP-RunIISummer20UL17pLHEGEN-00002.root --conditions 106X_mc2017_realistic_v6 --beamspot Realistic25ns13TeVEarly2017Collision --step GEN --geometry DB:Extended --filein file:TOP-RunIISummer20UL17pLHEGEN-00002_0.root --era Run2_2017 --no_exec --mc -n $EVENTS || exit $? ;


# cmsDriver command
cmsDriver.py Configuration/GenProduction/python/PYTHIA8_SingleDiffractiveTop_13TeV_cff.py --python_filename SD-TOP-RunIISummer20UL17_cfg.py --eventcontent RAWSIM --customise Configuration/DataProcessing/Utils.addMonitoring --datatier GEN --fileout file:SD-TOP-RunIISummer20UL17.root --conditions 106X_mc2017_realistic_v6 --beamspot Realistic25ns13TeVEarly2017Collision --step GEN --geometry DB:Extended --era Run2_2017 --no_exec --mc -n $EVENTS || exit $? ;
