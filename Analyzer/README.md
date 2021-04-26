# Missing Mass Searches Analyzer
Tool which produces the final samples used for Missing Mass Search analysis.

### Compilation

```sh
cd MissingMass/Analyzer
make
```

### Running Locally

You can use test the following command:

```sh
./MissingMassNtupleAnalyzer --f output2017.root --era B --mode Muon --jobid 1
```

Where the file 'output2017.root' has been produced by the Skimmer.

In order to check more options, run the option --help

```sh
./MissingMassNtupleAnalyzer --help
```

**Other options**

```
	 --f filename.root (input file)
	 --era B (B, C, D, E or F)
	 --mode Muon (or Electron)
	 --jobid job1 (tag to be added in the outputfile, _option_)
	 --outdir Output (Output dir for jobs, _option_)
	 --protonfile (it generates a text file with info from protons)
	 --eventfile (it generates a text file with run:ls:event_number format)
	 --random (performing analysis using random protons)
	 --single (fixing arm 45 with random protons.)
	 --zerobias (option to run with zerobias triggers. It includes the option --noprotonsfilter)
	 --noprotonsfilter (option to run with proton selection)

```

### Running on HTC Condor

In order to produce your own sample for the final physics analysis, you can use a script which automatically produces condor jobs per each crab output file (produced by the Skimmer). @CERN, HTC condor is very recommended. 

First, you need to generate a template xml file which includes automatically all the EOS folders output from your crab jobs. Please, note that in our case, the crab outputs are in "/eos/cms/store/group/phys_exotica/PPS-Exo/":

```sh
source CreateTemplateXML.sh
```

After that, you must fill the template file with the needed inputs. Please, check the example file "condor_samples.xml". If you would like to use extra parameters for MissingMassNtupleAnalyzer, you should include them in &lt;parameters&gt;--random --single &lt;/parameters&gt; for instance.

```sh
python CondorTool.py

Type ? to list commands
condor_submission> help

Documented commands (type help <topic>):
========================================
EOF  exit  help  kill  status  submit
```

#### Condor Commands

```sh
condor_submit job_condor.sub (Submit with automatic output transfer)
watch condor_q (Check job)
condor_q (Check job)
condor_submit -spool job_condor.sub (Submit, but not automatic output transfer)
condor_transfer_data $LOGNAME -const 'JobStatus == 4' (Retrieve data when submitted with -spool option)
condor_rm -all (kill all jobs)
```

**More information here:** [https://twiki.cern.ch/twiki/bin/view/ABPComputing/LxbatchHTCondor](https://twiki.cern.ch/twiki/bin/view/ABPComputing/LxbatchHTCondor)

**Gui for Monitoring:** [https://monit-grafana.cern.ch/d/000000869/user-batch-jobs?orgId=5&refresh=5m&var-cluster=cernprod](https://monit-grafana.cern.ch/d/000000869/user-batch-jobs?orgId=5&refresh=5m&var-cluster=cernprod)


### Debugging

```sh
gdb MissingMassNtupleAnalyzer
run --f output2018.root --era C --mode Muon
```
