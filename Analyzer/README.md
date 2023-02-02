# Missing Mass Searches Skimmer
Skimmer code includes Forward Proton POG containers. Furthermore, it includes jets, leptons, particle flow collections used for Missing Mass Searches.
Instructions how to produce your ntuples:

## Running Locally

```sh
cmsrel CMSSW_10_6_0
cd CMSSW_10_6_0/src
cmsenv
cd MissingMass/Skimmer/test
cmsRun RunMissingMassSearchesAOD.py Mode=Muon (or Electron) Year=2017 (or 2018) [AOD]
cmsRun RunMissingMassSearches.py Mode=Muon (or Electron) Year=2017 (or 2018) [miniAOD]
```

## Running Locally

```sh
cmsRun RunMissingMassSearches.py year=2017 physics=bjet trigger=True era=D mode=data
```

## Running on Condor

This option is to submit condor jobs which will execute all the N output files from crab (outputs of your skimmer) in a paralelized way.

```sh
source /cvmfs/cms.cern.ch/crab3/crab.sh
```

Then, setup the tool JobSubmitter accordingly the instructions (https://github.com/dfigueiredo/JobSubmitter). As an example:

```sh
python SubmitterTool.py --f samples_condor_doublemuon.json
mode condor
submit
```

Where the samples_condor_doublemuon.json has the following contents: https://github.com/dfigueiredo/JobSubmitter/blob/master/samples_condor_doublemuon.json
The same idea can be follow-up for all the samples_condor files.
