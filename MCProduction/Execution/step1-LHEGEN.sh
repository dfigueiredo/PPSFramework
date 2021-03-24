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
curl -s -k https://cms-pdmv.cern.ch/mcm/public/restapi/requests/get_fragment/TOP-RunIISummer20UL17pLHEGEN-00002 --retry 3 --create-dirs -o Configuration/GenProduction/python/RunIISummer20UL17LHEGEN-fragment.py
[ -s Configuration/GenProduction/python/RunIISummer20UL17LHEGEN-fragment.py ] || exit $?;
scram b
cd ../..

EVENTS=1000

# cmsDriver command
cmsDriver.py Configuration/GenProduction/python/RunIISummer20UL17LHEGEN-fragment.py --python_filename RunIISummer20UL17LHE_cfg.py --eventcontent LHE --customise Configuration/DataProcessing/Utils.addMonitoring --datatier LHE --fileout file:RunIISummer20UL17LHE.root --conditions 106X_mc2017_realistic_v6 --step NONE --filein "lhe_file.lhe" --era Run2_2017 --no_exec --mc -n $EVENTS || exit $? ;

# cmsDriver command
cmsDriver.py Configuration/GenProduction/python/RunIISummer20UL17LHEGEN-fragment.py --python_filename RunIISummer20UL17GEN_cfg.py --eventcontent RAWSIM --customise Configuration/DataProcessing/Utils.addMonitoring --datatier GEN --fileout file:RunIISummer20UL17GEN.root --conditions 106X_mc2017_realistic_v6 --beamspot Realistic25ns13TeVEarly2017Collision --step GEN --geometry DB:Extended --filein file:RunIISummer20UL17LHE.root --era Run2_2017 --no_exec --mc -n $EVENTS || exit $? ;
