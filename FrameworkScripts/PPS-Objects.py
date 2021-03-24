import FWCore.ParameterSet.Config as cms
import sys

import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing ()
options.register('year',
  '',
  VarParsing.VarParsing.multiplicity.singleton,
  VarParsing.VarParsing.varType.int,
  "year")
options.register('era',
  '',
  VarParsing.VarParsing.multiplicity.singleton,
  VarParsing.VarParsing.varType.string,
  "CMS era")
options.register('sampleTag',
  '',
  VarParsing.VarParsing.multiplicity.singleton,
  VarParsing.VarParsing.varType.string,
  "Tag of the sample to analyze")
options.register('outputDir',
  '',
  VarParsing.VarParsing.multiplicity.singleton,
  VarParsing.VarParsing.varType.string,
  "Output directory")
options.register('partNumber',
  '',
  VarParsing.VarParsing.multiplicity.singleton,
  VarParsing.VarParsing.varType.int,
  "Number of the part")
options.register('nParts',
  '',
  VarParsing.VarParsing.multiplicity.singleton,
  VarParsing.VarParsing.varType.int,
  "Number of parts")
options.partNumber = 0
options.nParts = 0
options.parseArguments() 

if ((options.partNumber == 0) ^ (options.nParts == 0)):
  print "ERROR: either none or both partNumber and nParts must be set"
  sys.exit(1)
elif ((options.partNumber == 0) & (options.nParts == 0)):
  options.partNumber = 1
  options.nParts = 1

###########################################################
# possible sampleTags:                                    #
# WW_A0W1e-6, WW_A0W2e-6, WW_A0W5e-6                      #
# WW_ACW2e-5, WW_ACW5e-6, WW_ACW8e-6                      #
# ZZ_A0Z-1e-5, ZZ_A0Z1e-5, ZZ_A0Z5e-5                     #
# ZZ_ACZ-1e-5, ZZ_ACZ1e-5, ZZ_ACZ5e-5                     #
###########################################################

MC=True
YEAR=options.year
ERA=options.era
MINIAOD=False
SAMPLETAG=options.sampleTag
OUTPUTDIR=options.outputDir
PARTNUMBER=options.partNumber
NPARTS=options.nParts
DoTheorySystematics=True
UseMCProtons=True

print "\nRunning with year = "+str(YEAR)+" era = "+ERA+" sampleTag = "+SAMPLETAG+"\n"

from Configuration.Eras.Modifier_ctpps_2016_cff import ctpps_2016
from Configuration.Eras.Modifier_ctpps_2017_cff import ctpps_2017
from Configuration.Eras.Modifier_ctpps_2018_cff import ctpps_2018
from Configuration.ProcessModifiers.run2_miniAOD_UL_cff import run2_miniAOD_UL

# load common code
if MC == True and YEAR == 2016 and (ERA == "B" or ERA == "C" or ERA == "G"):
   process = cms.Process('CTPPSTestAcceptance', ctpps_2016,run2_miniAOD_UL)
   from Validation.CTPPS.simu_config.year_2016_preTS2_cff import *
if MC == True and YEAR == 2016 and (ERA == "H"):
   process = cms.Process('CTPPSTestAcceptance', ctpps_2016,run2_miniAOD_UL)
   from Validation.CTPPS.simu_config.year_2016_postTS2_cff import *
if MC == True and YEAR == 2017 and (ERA == "B" or ERA == "C" or ERA == "D"):
   process = cms.Process('CTPPSTestAcceptance', ctpps_2017,run2_miniAOD_UL)
   from Validation.CTPPS.simu_config.year_2017_preTS2_cff import *
if MC == True and YEAR == 2017 and (ERA == "E" or ERA == "F"):
   process = cms.Process('CTPPSTestAcceptance', ctpps_2017,run2_miniAOD_UL)
   from Validation.CTPPS.simu_config.year_2017_postTS2_cff import *
if MC == True and YEAR == 2018 and (ERA == "A" or ERA == "B1"):
   process = cms.Process('CTPPSTestAcceptance', ctpps_2018,run2_miniAOD_UL)
   from Validation.CTPPS.simu_config.year_2018_cff import *
if MC == True and YEAR == 2018 and (ERA == "B2" or ERA == "C" or ERA == "D1"):
   process = cms.Process('CTPPSTestAcceptance', ctpps_2018,run2_miniAOD_UL)
   from Validation.CTPPS.simu_config.year_2018_cff import *
if MC == True and YEAR == 2018 and (ERA == "D2"):
   process = cms.Process('CTPPSTestAcceptance', ctpps_2018,run2_miniAOD_UL)
   from Validation.CTPPS.simu_config.year_2018_cff import *

if MC == True and YEAR == 2016 and (ERA == "B" or ERA == "C" or ERA == "G"):
   process.load("Validation.CTPPS.simu_config.year_2016_preTS2_cff")
if MC == True and YEAR == 2016 and (ERA == "H"):
   process.load("Validation.CTPPS.simu_config.year_2016_postTS2_cff")
if MC == True and YEAR == 2017 and (ERA == "B" or ERA == "C" or ERA == "D"):
   process.load("Validation.CTPPS.simu_config.year_2017_preTS2_cff")
if MC == True and YEAR == 2017 and (ERA == "E" or ERA == "F"):
   process.load("Validation.CTPPS.simu_config.year_2017_postTS2_cff")
if MC == True and YEAR == 2018 and (ERA == "A" or ERA == "B1"):
   process.load("Validation.CTPPS.simu_config.year_2018_cff")
   process.ctppsRPAlignmentCorrectionsDataESSourceXML.MisalignedFiles = ["PPtoPPWWjets/PPtoPPWWjets/python/PPS_2018_Alignments/2018_preTS1.xml"]
   process.ctppsRPAlignmentCorrectionsDataESSourceXML.RealFiles = ["PPtoPPWWjets/PPtoPPWWjets/python/PPS_2018_Alignments/2018_preTS1.xml"]
if MC == True and YEAR == 2018 and (ERA == "B2" or ERA == "C" or ERA == "D1"):
   process.load("Validation.CTPPS.simu_config.year_2018_cff")
   process.ctppsRPAlignmentCorrectionsDataESSourceXML.MisalignedFiles = ["PPtoPPWWjets/PPtoPPWWjets/python/PPS_2018_Alignments/2018_TS1_TS2.xml"]
   process.ctppsRPAlignmentCorrectionsDataESSourceXML.RealFiles = ["PPtoPPWWjets/PPtoPPWWjets/python/PPS_2018_Alignments/2018_TS1_TS2.xml"]
if MC == True and YEAR == 2018 and (ERA == "D2"):
   process.load("Validation.CTPPS.simu_config.year_2018_cff")
   process.ctppsRPAlignmentCorrectionsDataESSourceXML.MisalignedFiles = ["PPtoPPWWjets/PPtoPPWWjets/python/PPS_2018_Alignments/2018_postTS2.xml"]
   process.ctppsRPAlignmentCorrectionsDataESSourceXML.RealFiles = ["PPtoPPWWjets/PPtoPPWWjets/python/PPS_2018_Alignments/2018_postTS2.xml"]

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag

process.load('Configuration.EventContent.EventContent_cff')

##KS added GT splitting by year
if MC == True and YEAR == 2016:
   # JH - change this for AOD instead of MiniAOD, according to PDMV twiki                                                                                                         
   process.GlobalTag.globaltag ='80X_mcRun2_asymptotic_2016_TrancheIV_v8'
   #   process.GlobalTag.globaltag ='94X_mcRun2_asymptotic_v3,' # This one is for Summer16 MC campaign
if MC == True and YEAR == 2017:
   process.GlobalTag.globaltag ='94X_mc2017_realistic_v17' # this one for 94X MC
if MC == True and YEAR == 2018:
   process.GlobalTag.globaltag ='102X_upgrade2018_realistic_v20' # this one is for Autumn18 MC campaign

# load the right optical functions with python acrobatics
if MC == True and YEAR == 2017:
  process.ctppsOpticalFunctionsESSource.configuration[0].opticalFunctions = cms.VPSet(
    cms.PSet( xangle = cms.double(120), fileName = cms.FileInPath("CalibPPS/ESProducers/data/optical_functions/2017/version5tim/120urad.root") ),
    cms.PSet( xangle = cms.double(130), fileName = cms.FileInPath("CalibPPS/ESProducers/data/optical_functions/2017/version5tim/130urad.root") ),
    cms.PSet( xangle = cms.double(140), fileName = cms.FileInPath("CalibPPS/ESProducers/data/optical_functions/2017/version5tim/140urad.root") )
  )
if MC == True and YEAR == 2018:
  process.ctppsOpticalFunctionsESSource.configuration[0].opticalFunctions = cms.VPSet(
    cms.PSet( xangle = cms.double(120), fileName = cms.FileInPath("CalibPPS/ESProducers/data/optical_functions/2018/version6/120urad.root") ),
    cms.PSet( xangle = cms.double(130), fileName = cms.FileInPath("CalibPPS/ESProducers/data/optical_functions/2018/version6/130urad.root") ),
    cms.PSet( xangle = cms.double(140), fileName = cms.FileInPath("CalibPPS/ESProducers/data/optical_functions/2018/version6/140urad.root") )
  )

from signalSamples_cfi import *
print "Signal sample file list: "
print "\n".join(signalSamples(YEAR,ERA,SAMPLETAG,PARTNUMBER,NPARTS))
print "\n"
process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring(*signalSamples(YEAR,ERA,SAMPLETAG,PARTNUMBER,NPARTS))
    fileNames = cms.untracked.vstring("file:../../../../../../TOP-RunIISummer20UL17MiniAOD-00006.root")
)

# add pre-mixing of recHits
# remove the sub-era tag for 2018, to comply with the PU samples
pu_ERA = ERA
if len(ERA) == 2 and YEAR == 2018:
  pu_ERA = ERA[0]

from pileupSamples_cfi import *
print "Mixing PU samples: "
print "\n".join(pileupSamples(YEAR,pu_ERA)[:3])
print "...\n" 

process.load("protonPreMix.protonPreMix.ctppsPreMixProducer_cfi")
process.ctppsPreMixProducer.PUFilesList = cms.vstring(*pileupSamples(YEAR,pu_ERA))
process.ctppsPreMixProducer.Verbosity = 0

import FWCore.PythonUtilities.LumiList as LumiList
# use JSON for 2018 PU in order to mix the right events with the right period
if MC == True and YEAR == 2018:
  pu_jsonFile = "/eos/project-c/ctpps/Operations/DataExternalConditions/2018/CMSgolden_2RPGood_anyarms_Era"+ERA+".json"
  process.ctppsPreMixProducer.lumisToProcess = LumiList.LumiList(filename = pu_jsonFile).getVLuminosityBlockRange()
  print "Using JSON file for PU: "+pu_jsonFile+"\n"
  process.ctppsPreMixProducer.includeStrips = False

# rng service for premixing
process.RandomNumberGeneratorService.ctppsPreMixProducer = cms.PSet(initialSeed = cms.untracked.uint32(42))

# number of events
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport = cms.untracked.PSet(
  optionalPSet = cms.untracked.bool(True),
  reportEvery = cms.untracked.int32(100),
  limit = cms.untracked.int32(2000)
)

process.MessageLogger.destinations.append('logFile')
process.MessageLogger.logFile = cms.untracked.PSet(
  output = cms.untracked.string(OUTPUTDIR+"ExclWWjets_"+SAMPLETAG+"_Part"+str(PARTNUMBER)+"of"+str(NPARTS)+".out"),
  optionalPSet = cms.untracked.bool(True),
    INFO = cms.untracked.PSet(
      limit = cms.untracked.int32(0)
     ),
  noTimeStamps = cms.untracked.bool(False),
  FwkReport = cms.untracked.PSet(
    optionalPSet = cms.untracked.bool(True),
    reportEvery = cms.untracked.int32(1000),
    limit = cms.untracked.int32(10000000)
    ),
  default = cms.untracked.PSet(
    limit = cms.untracked.int32(10000000)
    ),
  Root_NoDictionary = cms.untracked.PSet(
    optionalPSet = cms.untracked.bool(True),
    limit = cms.untracked.int32(0)
    ),
  FwkSummary = cms.untracked.PSet(
    optionalPSet = cms.untracked.bool(True),
    reportEvery = cms.untracked.int32(1),
    limit = cms.untracked.int32(10000000)
    ),
  threshold = cms.untracked.string('INFO')
)

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

'''
### ADD SOME NEW JET COLLECTIONS                                                                                                                                           
# New (March 8, 2019) - to recover ak8 CHS jets with 2017 MiniAOD                                                                                                          
#################################                                                                                                                                          
###  JET TOOLBOX FOR CHS ###                                                                                                                                               
#################################                                                                                                                                          
# AK R=0.8 jets from CHS inputs with basic grooming, W tagging, and top tagging                                                                                            
from JMEAnalysis.JetToolbox.jetToolbox_cff import *

jetToolbox( process, 'ak8', 'ak8JetSubs', 'noOutput',
  PUMethod='CHS',
  addPruning=True, addSoftDrop=False ,           # add basic grooming                                                                                            
  addTrimming=False, addFiltering=False,
  addPrunedSubjets=False, addSoftDropSubjets=False,
  addNsub=True, maxTau=4,                       # add Nsubjettiness tau1, tau2, tau3, tau4                                                                       
  #miniAOD = MINIAOD,                                                                                                                                            
  dataTier = 'AOD',
  runOnMC=MC,
  bTagDiscriminators = None,  # blank means default list of discriminators, None means none                                                                      
  bTagInfos = None,
  subjetBTagDiscriminators = None,
  subjetBTagInfos = None,
  # added L1FastJet on top of the example config file                                                                                                            
  JETCorrPayload = 'AK8PFchs', JETCorrLevels = ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']
)


#################################
###  PT CUT 
#################################
process.goodJets = cms.EDFilter("CandViewSelector",
  src = cms.InputTag("selectedPatJetsAK8PFCHS"),
  cut = cms.string("pt > 30.0"),
  filter = cms.bool(True)
)


#################################                                                                                                                                                   
###  JER SMEARING ###                                                                                                                                                          
#################################                                                                                                                                                    
 
from RecoMET.METProducers.METSigParams_cfi import *
process.slimmedAK8JetsSmeared = cms.EDProducer('SmearedPATJetProducer',
  src = cms.InputTag("slimmedJetsAK8JetId"),
                                         enabled = cms.bool(True),
  rho = cms.InputTag("fixedGridRhoFastjetAll"),
  algo = cms.string('AK8PFchs'),
  algopt = cms.string('AK8PFchs_pt'),
  genJets = cms.InputTag('ak8GenJets'),
  # genJets = cms.InputTag('slimmedGenJetsAK8'),                                                                                
  dRMax = cms.double(0.4),
  dPtMaxFactor = cms.double(3),
  seed = cms.uint32(37428479),
  debug = cms.untracked.bool(False),
  # Systematic variation                                                                                                                                                 
  # 0: Nominal                                                                                                                                                           
  # -1: -1 sigma (down variation)                                                                                                                                        
  # 1: +1 sigma (up variation)                                                                                                                                           
  variation = cms.int32(0)  # If not specified, default to 0                                                                                                             
)

#################################                                                                                                                                          
###  JET ID  ###                                                                                                                                                           
#################################                                                                                                                                          
from PhysicsTools.SelectorUtils.pfJetIDSelector_cfi import pfJetIDSelector
process.slimmedJetsAK8JetId = cms.EDFilter("PFJetIDSelectionFunctorFilter",
  filterParams = pfJetIDSelector.clone(),
  src = cms.InputTag("selectedPatJetsAK8PFCHS"),
  filter = cms.bool(True)
)
'''


#########################                                                                                                                                                  
# override LHCInfo source
# for direct sim                                                                                                                                                           
#########################                   
process.load("CalibPPS.ESProducers.ctppsLHCInfoRandomXangleESSource_cfi")

#JH - debugging small shift                                                                                                                                                           
process.load("CalibPPS.ESProducers.ctppsBeamParametersFromLHCInfoESSource_cfi")
process.ctppsBeamParametersFromLHCInfoESSource.beamDivX45 = process.ctppsBeamParametersESSource.beamDivX45
process.ctppsBeamParametersFromLHCInfoESSource.beamDivX56 = process.ctppsBeamParametersESSource.beamDivX56
process.ctppsBeamParametersFromLHCInfoESSource.beamDivY45 = process.ctppsBeamParametersESSource.beamDivY45
process.ctppsBeamParametersFromLHCInfoESSource.beamDivY56 = process.ctppsBeamParametersESSource.beamDivY56

process.ctppsBeamParametersFromLHCInfoESSource.vtxOffsetX45 = process.ctppsBeamParametersESSource.vtxOffsetX45
process.ctppsBeamParametersFromLHCInfoESSource.vtxOffsetX56 = process.ctppsBeamParametersESSource.vtxOffsetX56
process.ctppsBeamParametersFromLHCInfoESSource.vtxOffsetY45 = process.ctppsBeamParametersESSource.vtxOffsetY45
process.ctppsBeamParametersFromLHCInfoESSource.vtxOffsetY56 = process.ctppsBeamParametersESSource.vtxOffsetY56
process.ctppsBeamParametersFromLHCInfoESSource.vtxOffsetZ45 = process.ctppsBeamParametersESSource.vtxOffsetZ45
process.ctppsBeamParametersFromLHCInfoESSource.vtxOffsetZ56 = process.ctppsBeamParametersESSource.vtxOffsetZ56

process.ctppsBeamParametersFromLHCInfoESSource.vtxStddevX = process.ctppsBeamParametersESSource.vtxStddevX
process.ctppsBeamParametersFromLHCInfoESSource.vtxStddevY = process.ctppsBeamParametersESSource.vtxStddevY
process.ctppsBeamParametersFromLHCInfoESSource.vtxStddevZ = process.ctppsBeamParametersESSource.vtxStddevZ

# do not apply vertex smearing again
process.ctppsBeamParametersESSource.vtxStddevX = 0
process.ctppsBeamParametersESSource.vtxStddevY = 0
process.ctppsBeamParametersESSource.vtxStddevZ = 0

# undo CMS vertex shift
#process.ctppsBeamParametersESSource.vtxOffsetX45 = +0.2475 * 1E-1
#process.ctppsBeamParametersESSource.vtxOffsetY45 = -0.6924 * 1E-1
#process.ctppsBeamParametersESSource.vtxOffsetZ45 = -8.1100 * 1E-1
# These should be the vtx shift values for 2017...
process.ctppsBeamParametersESSource.vtxOffsetX45 = +0.024793
process.ctppsBeamParametersESSource.vtxOffsetY45 = -0.0692861
process.ctppsBeamParametersESSource.vtxOffsetZ45 = -7.89895
# end JH

process.ctppsLHCInfoRandomXangleESSource.generateEveryNEvents = 1
if MC == True and YEAR == 2016:
    process.ctppsLHCInfoRandomXangleESSource.xangleHistogramFile = "CrossingAngles2016.root"
if MC == True and YEAR == 2017:
    process.ctppsLHCInfoRandomXangleESSource.xangleHistogramFile = "CrossingAngles2017.root"
if MC == True and YEAR == 2018:
    process.ctppsLHCInfoRandomXangleESSource.xangleHistogramFile = "CrossingAngles2018.root"

process.ctppsLHCInfoRandomXangleESSource.xangleHistogramObject = "hxang"
process.ctppsLHCInfoRandomXangleESSource.beamEnergy = 6500.
process.ctppsLHCInfoRandomXangleESSource.betaStar = 0.40

process.esPreferLHCInfo = cms.ESPrefer("CTPPSLHCInfoRandomXangleESSource", "ctppsLHCInfoRandomXangleESSource")

process.beamDivergenceVtxGenerator.src = cms.InputTag("")
process.beamDivergenceVtxGenerator.srcGenParticle = cms.VInputTag(
  # Don't propagate PU protons
  # cms.InputTag("genPUProtons","genPUProtons")
  cms.InputTag("prunedGenParticles")
)

# point the track producers to the correct products
#process.ctppsPixelLocalTracks.label = "ctppsPreMixProducer"
#process.ctppsPixelLocalTracks.label = "ctppsPixelRecHits"
#process.totemRPUVPatternFinder.tagRecHit = cms.InputTag("ctppsPreMixProducer")
#process.totemRPUVPatternFinder.tagRecHit = cms.InputTag("totemRPRecHitProducer","","CTPPSTestAcceptance")

# add PPS efficiency simulation
# remove timing tracks from trackLites, they are not produced by protonPreMix
#process.ctppsLocalTrackLiteProducer.includeDiamonds = False

# reconstruction (if 2016 data, remove modules for RPs which did not exist at that time)
def RemoveModules(pr):
  pr.reco_local.remove(pr.ctppsPixelLocalTracks)
  pr.ctppsLocalTrackLiteProducer.includePixels = False
  pr.ctppsPreMixProducer.includePixels = False

(ctpps_2016 & ~ctpps_2017 & ~ctpps_2018).toModify(process, RemoveModules)

# rng service for efficiency
process.load("protonPreMix.protonPreMix.ppsEfficiencyProducer_cfi")

process.RandomNumberGeneratorService.ppsEfficiencyProducer = cms.PSet(initialSeed = cms.untracked.uint32(43))
process.ppsEfficiencyProducer.year = cms.int32(YEAR)
process.ppsEfficiencyProducer.era = cms.string(ERA)

if YEAR == 2016:
  process.ppsEfficiencyProducer.efficiencyFileName_Near = cms.string("/eos/project/c/ctpps/subsystems/Strips/StripsTracking/PreliminaryEfficiencies_July132020_1D2DMultiTrack.root")
  process.ppsEfficiencyProducer.efficiencyFileName_Far = cms.string("/eos/project/c/ctpps/subsystems/Strips/StripsTracking/PreliminaryEfficiencies_July132020_1D2DMultiTrack.root")
if YEAR == 2017:
  process.ppsEfficiencyProducer.efficiencyFileName_Near = cms.string("/eos/project/c/ctpps/subsystems/Strips/StripsTracking/PreliminaryEfficiencies_July132020_1D2DMultiTrack.root")
  process.ppsEfficiencyProducer.efficiencyFileName_Far = cms.string("/eos/project/c/ctpps/subsystems/Pixel/RPixTracking/pixelEfficiencies_singleRP.root")
  #process.ppsEfficiencyProducer.efficiencyFileName_Far = cms.string("/afs/cern.ch/user/a/abellora/cernbox/PPS/Efficiency_reMiniAOD/pixelEfficiencies_singleRP_reMiniAOD.root")

if YEAR == 2018:
  process.ppsEfficiencyProducer.efficiencyFileName_Near = cms.string("/eos/project/c/ctpps/subsystems/Pixel/RPixTracking/pixelEfficiencies_radiation.root")
  #process.ppsEfficiencyProducer.efficiencyFileName_Near = cms.string("/afs/cern.ch/user/a/abellora/cernbox/PPS/Efficiency_reMiniAOD/pixelEfficiencies_radiation_reMiniAOD.root")  
  process.ppsEfficiencyProducer.efficiencyFileName_Far = cms.string("/eos/project/c/ctpps/subsystems/Pixel/RPixTracking/pixelEfficiencies_radiation.root")
  #process.ppsEfficiencyProducer.efficiencyFileName_Far = cms.string("/afs/cern.ch/user/a/abellora/cernbox/PPS/Efficiency_reMiniAOD/pixelEfficiencies_singleRP_reMiniAOD.root")

#########################                                                                                                                                                           
# Configure Analyzer
#########################                                                                                                                                                           
process.load("PPtoPPWWjets.PPtoPPWWjets.CfiFile_cfi")
process.TFileService = cms.Service("TFileService", fileName = cms.string(OUTPUTDIR+"ExclWWjets_"+SAMPLETAG+"_Part"+str(PARTNUMBER)+"of"+str(NPARTS)+".root"))
print "Output will be saved in: "+OUTPUTDIR+"ExclWWjets_"+SAMPLETAG+"_Part"+str(PARTNUMBER)+"of"+str(NPARTS)+".root"
print "Log will be saved in: "+OUTPUTDIR+"ExclWWjets_"+SAMPLETAG+"_Part"+str(PARTNUMBER)+"of"+str(NPARTS)+".out"
print ""

process.demo = cms.EDAnalyzer('PPtoPPWWjets')

if YEAR == 2016:
  process.load("PPtoPPWWjets.PPtoPPWWjets.HLTFilter2016_cfi")
  process.hltFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
if YEAR == 2017:
  process.load("PPtoPPWWjets.PPtoPPWWjets.HLTFilter_cfi")
  process.hltFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
if YEAR == 2018:
  process.load("PPtoPPWWjets.PPtoPPWWjets.HLTFilter2018_cfi")
  process.hltFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")


if MC == True and YEAR == 2016:
  process.demo.dataPileupFile = cms.string("PUHistos_data_2016.root")
  process.demo.mcPileupFile = cms.string("PUHistos_mc_2016.root")
if MC == True and YEAR == 2017:
  process.demo.dataPileupFile = cms.string("PUHistos_data_2017.root")
  process.demo.mcPileupFile = cms.string("PUHistos_mc_2017.root")
  # Special cases for buggy datasets in 2017 MC                                                                                                
  #    process.demo.mcPileupFile = cms.string("PUHistos_mc_2017_QCDPt300to470.root")
  #    process.demo.mcPileupFile = cms.string("PUHistos_mc_2017_QCDPt600to800.root")
  #    process.demo.mcPileupFile = cms.string("PUHistos_mc_2017_QCDPt1000to1400.root")                                                         
  #    process.demo.mcPileupFile = cms.string("PUHistos_mc_2017_Wjets.root")
  #    process.demo.mcPileupFile = cms.string("PUHistos_mc_2017_Zjets.root")     
if MC == True and YEAR == 2018:
  process.demo.dataPileupFile = cms.string("PUHistos_data_2018.root")
  process.demo.mcPileupFile = cms.string("PUHistos_mc_2018.root")

process.demo.jetAK8CHSCollection = cms.InputTag("slimmedJetsAK8JetId")
process.demo.ppsRecoProtonSingleRPTag = cms.InputTag("ctppsProtons", "singleRP")
# process.demo.ppsRecoProtonMultiRPTag = cms.InputTag("ctppsProtons", "multiRP")
process.demo.ppsRecoProtonMultiRPTag = cms.InputTag("ppsEfficiencyProducer", "multiRP")
process.demo.isMC = cms.bool(MC)
process.demo.doMCTheorySystematics = cms.bool(DoTheorySystematics)
process.demo.useMCProtons = cms.bool(UseMCProtons)
process.demo.year = cms.int32(YEAR)
process.demo.era = cms.string(ERA)

process.maxEvents = cms.untracked.PSet(
  # input = cms.untracked.int32(20)
  input = cms.untracked.int32(-1)
)

# reconstruction plotter (analysis example)
# process.ctppsProtonReconstructionPlotter = cms.EDAnalyzer("CTPPSProtonReconstructionPlotter",
#   tagTracks = cms.InputTag("ctppsLocalTrackLiteProducer"),
#   tagRecoProtonsSingleRP = cms.InputTag("ctppsProtons", "singleRP"),
#   # tagRecoProtonsMultiRP = cms.InputTag("ctppsProtons", "multiRP"),
#   tagRecoProtonsMultiRP = cms.InputTag("ppsEfficiencyProducer", "multiRP"),
#   rpId_45_F = cms.uint32(23),
#   rpId_45_N = cms.uint32(3),
#   rpId_56_N = cms.uint32(103),
#   rpId_56_F = cms.uint32(123),

#   association_cuts_45 = process.ctppsProtons.association_cuts_45,
#   association_cuts_56 = process.ctppsProtons.association_cuts_56,

#   outputFile = cms.string(OUTPUTDIR+"output_protons.root")
# )


#process.output = cms.OutputModule("PoolOutputModule",
# fileName = cms.untracked.string("output_data.root"),
# outputCommands = cms.untracked.vstring("keep *_ctpps*_*_*")
#)

#process.out = cms.EndPath(process.output)

#output file
process.out = cms.OutputModule('PoolOutputModule',
    fileName = cms.untracked.string('output.root'),
    outputCommands = process.MINIAODEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

process.out.outputCommands.extend(
 cms.untracked.vstring(
        'keep *_*_*_*',
        #'keep *_*_*_CTPPSTestAcceptance'
        'drop *_ctpps*_*_RECO'
    )
)

print(process.reco_local)

process.outpath = cms.EndPath(process.out)

# processing path
process.p = cms.Path(
  #process.hltFilter
  process.beamDivergenceVtxGenerator
  * process.ctppsDirectProtonSimulation
  #* process.ctppsPreMixProducer
  #* process.ctppsDiamondLocalReconstruction
  * process.reco_local
  * process.ctppsProtons
  #* process.ppsEfficiencyProducer

  # * process.ctppsProtonReconstructionPlotter

  #* process.genParticlesForJetsNoNu 
  #* process.goodOfflinePrimaryVertices  
  #* process.slimmedJetsAK8JetId 
  #* process.demo
)
