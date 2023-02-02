# Utilities

Simple and helpful scripts.

## Plotter

Root script example to make plots applying specific cuts from selected branches (from PPS Framework TTree). Please, insert your filename accordingly.

```bash
root -l PlotterExample.C
```

## CMS DAS API

Python script to automatically get the total number of events from a CMS dataset. Edit the script in order to insert the dataset name you are interested to check.

```bash
python GettingDASData.py
```

## Crab Check Job Status

Place this script in the same directory where all your crab submitter folders are located. This script will check the status for all the crab submissions.

```bash
source StatusJobsCrab.sh
```

## Useful commands

### Copying a file saved in any T2 to your local area. Example:

You can check the file path in the CMS DAS website (https://cmsweb.cern.ch/das/).

```bash
xrdcp root://cmsxrootd.fnal.gov//store/user/dmf/crab_dmf_2021-04-28_UTC11-16-01/PYTHIA8-SD-TOP-GEN/PYTHIA8-SD-TOP-MINIAOD-13TEV/210428_091727/0000/RunIISummer20UL17MiniAOD_141.root test_pythia_miniaod.root
```

### Getting files from a specific run and lumisection in cms DAS (for detector data):

```bash
 file dataset=/EGamma/Run2022E-v1/RAW run=359699 lumi=265
```
