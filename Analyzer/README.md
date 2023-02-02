# Missing Mass Searches Analyzer
Tool to produce the final ttree for analysis.

## Compiling the source 

Before running, you should compile the code:

```sh
make clean
make
```

## Running Locally

Example to run with detector data (produced by the Missing Mass Search Skimmer):

```sh
./MissingMassNtupleAnalyzer --f skimmer_output.root --year 2017 --era C --mode data --physics displacedjet --jobid 0
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
