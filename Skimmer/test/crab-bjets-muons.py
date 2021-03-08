import time
from CRABClient.UserUtilities import config

config = config()

#Set Run Type
stream = 'DoubleMuon' # DoubleMuon or BTagCSV
run_name = 'F'

timestr = time.strftime("%Y-%m-%d_UTC%H-%M-%S")
lfndir = '/store/group/phys_exotica/PPS-Exo'
mode = 'Muon'
year = '2017'

dataset_name = '/'+stream+'/Run2017'+run_name+'-09Aug2019_UL2017-v1/MINIAOD'

tagname = 'crab_' + str(stream) + "_" + str(run_name) + '_' + str(timestr)
pathfull = str(lfndir) + '/' + str(tagname) + '_' + str(timestr)

print dataset_name
print pathfull

config.General.transferLogs = True
config.General.transferOutputs = True
config.JobType.pluginName = 'Analysis'
config.JobType.maxMemoryMB = 2500
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 20
config.Data.publication = False
config.JobType.psetName = 'RunMissingMassSearches.py'
config.JobType.outputFiles = ['output.root']
config.General.workArea = str(tagname)
config.General.requestName = str(tagname)
config.Data.inputDataset = str(dataset_name)
config.Data.lumiMask = 'combined_RPIN_CMS_2017.json'
config.Data.outputDatasetTag = str(tagname)
config.Site.storageSite = 'T2_CH_CERN'
config.Data.outLFNDirBase = str(pathfull)
config.JobType.pyCfgParams = ["Mode="+str(mode),"Year="+str(year)]

