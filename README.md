# PPS Framework

This package is used for the Missing Mass Search analysis and acts on the CMSSW official datasets, producing a root NTuple which has pre-selected CMSSW collections (and physics-computing objects). Secondly, a tool which saves a root Ttree for final analysis. 

## Description

| Folder Name       | Comments |
| ------------- |:-------------:|
| Skimmer | to store CMSSW pre-selected objects in a NTuple (step 1)|
| Analyzer | produce the final analysis Ttree (step 2) |
| Plotter | produce automatic plots for your final analysis (using Ttree format) |
| Utils | helpful scripts and commands |

## Installation

Download the package from git. It is also important to set your voms-certificate as well as CMSSW (cmsenv) where your plugins are set. 

```bash
SCRAM_ARCH=slc7_amd64_gcc700; export SCRAM_ARCH
cmsrel CMSSW_10_6_17_patch1
cd CMSSW_10_6_17_patch1/src
cmsenv
git clone https://github.com/dfigueiredo/PPSFramework.git -b CMSSW_10_6_17_patch1_2023
cd PPSFramework/
scram b -j 8
source set_environment.sh
cd working
```

Once CMSSW and your package is installed, you need to setup everytime you connect in a new terminal:

```bash
SCRAM_ARCH=slc7_amd64_gcc700; export SCRAM_ARCH (or setenv SCRAM_ARCH slc7_amd64_gcc700 for tcsh) 
cd CMSSW_10_6_17_patch1/src
cmsenv
```

# Skimmer

Used to create the NTuple (Skimmer) from a CMSSW dataset type format (AOD or miniAOD).

## Usage

As an example, edit the file working/RunMissingMassSearches.py and substitute the line:

```bash
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
```

for


```bash
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5000) )
```

After that, run the following command:

```bash
cmsRun RunMissingMassSearches.py year=2017 physics=bjet trigger=True era=D mode=data
```

Please, pay attention that the options mode (data or mc), year (2017 or 2018) must be set correctly accordingly to the dataset which has been choosen. If you are not using crab, remember to change the input file name to run locally (in general for a test or to debug the code).

## Options

| Options       | Explanation | Comments |
| ------------- |:-------------:|-------------:|
| physics      | (string) muon, electron, bjet, displacedjet, emu, zerobias | to set the correct trigger path for a specific dataset stream |
| mode      | (string) data, mc | If a data sample, it must be set accordingly (same for monte carlo) |
| era      | (string) A, B, C, D | It depends of the year. Used to load PPS proton reconstruction parameters, cutoff or for simulation |
| year   | (string) 2017, 2018 | To deal with HLT names and conditions |
| prescales   | (bool) False or True | enable/disable trigger prescales values. Default is False (only for data) |
| trigger   | (bool) False or True | enable/disable the trigger. Default is True. |
| unmatching   | (bool) False or True | True: leptons and jets not associated in the same cone. Default is false |
| ppstagging   | (bool) False or True | True: selecting events with at least one proton per arm. Default is true |
| debugging   | (bool) False or True | True: enabling skimmer debugger. Default is false |

# Analyzer

The analyzer is used to produce a TTree after some analysis basic cuts and selecting events with one proton per arm. As an example, you can now run the Analyzer to process the NTuple produced by the previous skimmer (RunMissinMassSearches.py):

```bash
cd ../Analyzer
./MissingMassNtupleAnalyzer --f ../working/output.root --year 2017 --era D --mode data --physics bjet
```

Where all the options are following the Skimmer options (same year, era, mode, physics). The output.root file was produced by the Skimmer (RunMissingMassSearches.py). Now, the output file (NTuple_data_bjet.root) can be opened and the analysis variables can be plotted direclty from the Ttree (extra cuts can also be applied).

## Compilation

Everytime the source code is changed, you need to recompile it to run (MissingMassNtupleAnalyzer.cc):

```bash
make clean
make
```

## Options

| Options       | Comments |
| ------------- | -------------:|
| --protonfile | Used to create, from the detector data per physics (muon, electron, bjets), a root file with Multi and Single Xi histograms. This file is used for the --random option |
| --eventfile | Used to create, from the detector data, a text data format with golded selected events (to be visualized later on Fireworks) |
| --short | Option to process only the proton selection, skipping the jets or leptons selection |
| --notrigger | Skip the trigger selection |
| --random | If there is no proton in one arm, it adds randomly from the Xi data distribution (single and multi). Used for MC cases |
| --zerobias | Only for zerobias datasets to study background |
| --noppstagging | Without proton selection (one proton per arm) |
| --debugging | For debugging the selections. It prints some outputs on the screen |
| --help | It prints all the options accepted by the code |
| --mode | mc or data |
| --era | A,B,C,D,E or F |
| --physics | electron, muon, bjet or displacedjet (to enable each trigger menu) |
| --jobid | It adds a string in the end of the output filename. Used in condor submission to avoid file overwriting (string) |
| --output | It creates an output folder name (string) |
| --f | Input file name |

*** Important ***: the options "--protonfile" and "--random" _can not_ be used together! In case you do not have the random file with the protons xi, first run --protonfile to create the root file histogram used for the --random option.

## Examples

### Pythia8 SD Top

```bash
./MissingMassNtupleAnalyzer --f $(filename) --year 2017 --era C --mode mc --physics muon --random --jobid 0
```

### Toy MC (Higgs or Z)

```bash
./MissingMassNtupleAnalyzer --f $(filename) --year 2017 --era C --mode mc --physics bjet --jobid 0
```

### Data 

```bash
./MissingMassNtupleAnalyzer --f $(filename) --year 2017 --era C --mode data --physics displacedjet --jobid 0
```

### Creating Random File

```bash
./MissingMassNtupleAnalyzer --f $(filename) --year 2017 --era D --mode data --physics bjet --protonfile
```

the file RandomProtons_bjet_eraD.root will be created.
