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

## Running on Grid

First, you need to set your crab environment doing:

```sh
source /cvmfs/cms.cern.ch/crab3/crab.sh
```

This tool will submit on grid job tasks for a set of CMS datasets, which are defined in a XML file. As an example, check the file grid_samples_2017_aod.xml. In order to submit, please do [i.e]:

```sh
python GridTool.py -p "submit --file grid_samples_2017_timing_aod.xml" 
```

or

```sh
python GridTool.py
submit --file grid_samples_2017_timing_aod.xml [press enter]
```

For more options:

```sh
python GridTool.py --h
Usage: GridTool.py [options]

Options:
  -h, --help            show this help message and exit
  -f FILE, --filename=FILE
                        XML mapping file
  -p PARSING, --parsing=PARSING
                        parsing: commands which can be passed from SHELL
                        directly. [parsing: --p "submit --file filename.xml"]
  -v, --verbose         make lots of noise [default]
```
