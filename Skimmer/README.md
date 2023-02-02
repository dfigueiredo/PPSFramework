# Missing Mass Searches Skimmer
Skimmer code includes Forward Proton POG containers. Furthermore, it includes jets, leptons, particle flow collections used for Missing Mass Searches.
Instructions how to produce your ntuples:

## Running Locally

```sh
cmsRun RunMissingMassSearches.py year=2017 physics=bjet trigger=True era=D mode=data
```

## Running on Grid

First, you need to set your crab environment doing:

```sh
source /cvmfs/cms.cern.ch/crab3/crab.sh
```

Then, setup the tool JobSubmitter accordingly the instructions (https://github.com/dfigueiredo/JobSubmitter). As an example, to generate a skimmer from a CMS dataset:

```sh
python SubmitterTool.py --f samples_skimmer2018_bjets.json 
mode crab
submit
```

Where the samples_skimmer2018_bjets.json has the following contents: https://github.com/dfigueiredo/JobSubmitter/blob/master/samples_skimmer2018_bjets.json
The same idea can be follow-up for all the samples_skimmer files.

# Missing Mass Efficiency Studies

Tool to check analysis cuts and efficiencies.

```sh
cmsRun RunMissingMassEfficiency.py trigger=False mode=mc physics=muon debugging=True verbosity=0 discriminator=False
```
