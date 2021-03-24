#!/bin/bash

# Make voms proxy
voms-proxy-init --voms cms --out $(pwd)/voms_proxy.txt --hours 4
export X509_USER_PROXY=$(pwd)/voms_proxy.txt

export SCRAM_ARCH=slc7_amd64_gcc700

source /cvmfs/cms.cern.ch/cmsset_default.sh
if [ -r CMSSW_10_6_17/src ] ; then
  echo release CMSSW_10_6_17_patch1 already exists
else
  scram p CMSSW CMSSW_10_6_17
fi
cd CMSSW_10_6_17/src
eval `scram runtime -sh`

git cms-init
git cms-addpkg IOMC/EventVertexGenerators RecoCTPPS Validation/CTPPS CalibPPS/ESProducers
cd CalibPPS/ESProducers
git clone https://github.com/cms-data/CalibPPS-ESProducers.git data
cd ../..
git clone https://github.com/jjhollar/PPtoPPWWjets.git
cd PPtoPPWWjets/
cd ..
git clone https://github.com/cms-jet/JetToolbox.git JMEAnalysis/JetToolbox -b jetToolbox_102X_v3

curl -s -k https://raw.githubusercontent.com/dfigueiredo/PPSFramework/main/FrameworkScripts/jetToolbox_cff.py --retry 3 --create-dirs -o JMEAnalysis/JetToolbox/python/jetToolbox_cff.py
[ -s JMEAnalysis/JetToolbox/python/jetToolbox_cff.py ] || exit $?;
git clone https://github.com/AndreaBellora/protonPreMix.git
scram b -j8

cd ../..
