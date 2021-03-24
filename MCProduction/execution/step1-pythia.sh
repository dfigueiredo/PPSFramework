#!/bin/bash

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
curl -s -k https://raw.githubusercontent.com/dfigueiredo/PPSFramework/main/MCProduction/Configuration/PYTHIA8_SingleDiffractiveTop_13TeV_cff.py --retry 3 --create-dirs -o Configuration/GenProduction/python/pythia8-singlediffraction-top-fragment.py
[ -s Configuration/GenProduction/python/pythia8-singlediffraction-top-fragment.py ] || exit $?;
scram b
cd ../..

EVENTS=100

# cmsDriver command
cmsDriver.py Configuration/GenProduction/python/pythia8-singlediffraction-top-fragment.py --python_filename SD-TOP-PYTHIA8_cfg.py --eventcontent RAWSIM --customise Configuration/DataProcessing/Utils.addMonitoring --datatier GEN --fileout file:RunIISummer20UL17GEN.root --conditions 106X_mc2017_realistic_v6 --beamspot Realistic25ns13TeVEarly2017Collision --step GEN --geometry DB:Extended --era Run2_2017 --no_exec --mc -n $EVENTS || exit $? ;
