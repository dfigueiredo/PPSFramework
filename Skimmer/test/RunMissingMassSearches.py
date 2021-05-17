#flake8: noqa

#
# Forcing DB changes for Proton. The new optics is not yet into the GT
#
forceDB = False

import sys
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

# Setting Input Parameters from Line Command
options = VarParsing ('analysis')
options.register('physics','muon',VarParsing.multiplicity.singleton, VarParsing.varType.string,"physics(string): muon, electron, zerobias or emu. Default is muon.")
options.register('mode','data',VarParsing.multiplicity.singleton, VarParsing.varType.string,"mode(string): data or mc. Default is data.")
options.register('year','2017',VarParsing.multiplicity.singleton, VarParsing.varType.string,"year(string): 2017 or 2018. Default is 2017.")
options.register('prescale',False,VarParsing.multiplicity.singleton, VarParsing.varType.bool,"prescale(bool): enable/disable trigger prescales values. Default is False (it works only for data).")
options.register('trigger',True,VarParsing.multiplicity.singleton, VarParsing.varType.bool,"trigger(bool): enable/disable the trigger. Default is True.")
options.register('unmatching',False,VarParsing.multiplicity.singleton, VarParsing.varType.bool,"unmatching(bool): \n\t\tTrue: leptons and jets not associated in the same cone.\n\t\tFalse: otherwise\n. Default is False.")
options.register('ppstagging',True,VarParsing.multiplicity.singleton, VarParsing.varType.bool,"ppstagging(bool): \n\t\tTrue: selecting events with at least one proton per arm.\n\t\tFalse: no filter.\n. Default is True.")
options.register('debugging',False,VarParsing.multiplicity.singleton, VarParsing.varType.bool,"debugging(bool): \n\t\tTrue: enabling skimmer debugger.\n\t\tFalse: otherwise\n. Default is False.")
options.parseArguments()

#fileinput = '/store/data/Run2017B/DoubleMuon/MINIAOD/09Aug2019_UL2017-v1/260000/F1B17EC1-36AF-6340-B359-774C7D88948D.root'
#fileinput = '/store/data/Run2017D/BTagCSV/MINIAOD/09Aug2019_UL2017-v1/60000/CAEC3BF9-E990-AD41-A55B-B28269E242C1.root'
#fileinput = '/store/data/Run2018C/DoubleMuon/MINIAOD/ForValUL2018-v1/230000/1A9817AE-B1F4-A746-9460-34B48037B519.root'
#fileinput = '/store/user/dmf/crab_dmf_2021-05-07_UTC14-15-29/PYTHIA8-SD-TOP-GEN/PYTHIA8-SD-TOP-MINIAOD-withpps-13TEV/210507_121551/0000/output_10.root'
fileinput = '/store/user/dmf/PYTHIA8-SD-top-13TeV-miniAOD-xangle140-PostTS2-miniAOD_PostTS2_140_2020-05-08_UTC15-08-15/PYTHIA8-SD-top-13TeV/PYTHIA8-SD-top-13TeV-miniAOD-xangle140-PostTS2-miniAOD_PostTS2_140_2020-05-08_UTC15-08-15/200508_130830/0000/output_123.root' #testing old pythia8
#fileinput = 'root://xrootd-cms.infn.it//store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/90000/FEB3954C-4942-E811-8A09-008CFAC91A38.root'
#fileinput = '/store/data/Run2017C/DoubleEG/MINIAOD/09Aug2019_UL2017-v1/50000/F696E155-EF11-CB48-82CA-96EA48359D54.root'
#fileinput = '/store/data/Run2018C/EGamma/MINIAOD/ForValUL2018-v2/230000/5C5DBBE7-A5E1-E746-835A-D33E0C37071F.root'
#fileinput = '/store/data/Run2017E/MuonEG/MINIAOD/09Aug2019_UL2017-v1/30000/D68314AC-E5B8-BB43-A0F9-86F599A8DA82.root'
#fileinput = '/store/data/Run2018A/MuonEG/MINIAOD/12Nov2019_UL2018-v1/00000/51E69F77-C0F4-B942-B704-3826C9BA0F26.root'

process = cms.Process("analysis")

if options.year != "2017" and options.year !="2018":
  	  print("\n Please, use only years 2017 or 2018!\n")
	  sys.exit(0)

if options.year == "2018" and options.mode =="zerobias":
          print("\n Please, zerobias is not yet configured to run over 2018 data!\n")
          sys.exit(0)

if options.physics == "muon":
  print("")
  print("####")
  print("Muon")
  print("####")
  print("")
  triggerlist = 'HLT_IsoMu27_v*', 'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*', 'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v*', 'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*', 'HLT_DoubleMu43NoFiltersNoVtx_v*', 'HLT_DoubleMu48NoFiltersNoVtx_v*', 'HLT_BTagMu_AK4DiJet20_Mu5_v*', 'HLT_BTagMu_AK4DiJet40_Mu5_v*', 'HLT_BTagMu_AK4DiJet70_Mu5_v*', 'HLT_BTagMu_AK4DiJet110_Mu5_v*', 'HLT_BTagMu_AK4DiJet170_Mu5_v*', 'HLT_BTagMu_AK4Jet300_Mu5_v*', 'HLT_BTagMu_AK8DiJet170_Mu5_v*', 'HLT_BTagMu_AK8Jet300_Mu5_v*', 'HLT_PFHT380_SixJet32_DoubleBTagCSV_p075_v*', 'HLT_PFHT430_SixJet40_BTagCSV_p080_v*'
elif options.physics == "electron":
  print("")
  print("########")
  print("Electron")
  print("########")
  print("")
  triggerlist = 'HLT_Ele27_WPTight_Gsf_v*', 'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*', 'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*', 'HLT_DoubleEle33_CaloIdL_MW_v*', 'HLT_BTagMu_AK4DiJet20_Mu5_v*', 'HLT_BTagMu_AK4DiJet40_Mu5_v*', 'HLT_BTagMu_AK4DiJet70_Mu5_v*', 'HLT_BTagMu_AK4DiJet110_Mu5_v*', 'HLT_BTagMu_AK4DiJet170_Mu5_v*', 'HLT_BTagMu_AK4Jet300_Mu5_v*', 'HLT_BTagMu_AK8DiJet170_Mu5_v*', 'HLT_BTagMu_AK8Jet300_Mu5_v*', 'HLT_PFHT380_SixJet32_DoubleBTagCSV_p075_v*', 'HLT_PFHT430_SixJet40_BTagCSV_p080_v*','HLT_PFMET100_PFMHT100_IDTight_CaloBTagCSV_3p1_v*', 'HLT_PFMET110_PFMHT110_IDTight_CaloBTagCSV_3p1_v*', 'HLT_PFMET120_PFMHT120_IDTight_CaloBTagCSV_3p1_v*', 'HLT_PFMET130_PFMHT130_IDTight_CaloBTagCSV_3p1_v*','HLT_PFMET140_PFMHT140_IDTight_CaloBTagCSV_3p1_v*', 'HLT_Ele15_IsoVVVL_PFHT450_CaloBTagCSV_4p5_v*', 'DST_HT250_CaloBTagScouting_v*'
elif options.physics == "zerobias":
  print("")
  print("########")
  print("ZeroBias")
  print("########")
  print("")
  triggerlist = 'HLT_ZeroBias_v*', 'HLT_ZeroBias_FirstBXAfterTrain_v*','HLT_ZeroBias_IsolatedBunches_v*','HLT_ZeroBias_Alignment_v*','HLT_ZeroBias_Beamspot_v*','HLT_L1UnpairedBunchBptxMinus_v*','HLT_L1UnpairedBunchBptxPlus_v*','HLT_Physics_v*','HLT_L1SingleMu_v*'
elif options.physics == "emu":
  print("")
  print("###")
  print("EMu")
  print("###")
  print("")
  triggerlist = 'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*', 'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*', 'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*'
else:
  print("")
  print("########################################################################################################")
  print("Please, try physics==muon or physics=electron or physics=zerobias or physics=emu")
  print("########################################################################################################")
  print("")
  sys.exit(0)

mc = False
if options.mode == "data":
  print("################")
  print("Running for Data")
  print("################")
  print("")
elif options.mode == "mc":
  print("#######################")
  print("Running for Monte Carlo")
  print("#######################")
  print("")
  mc = True
else:
  print("")
  print("################################")
  print("Please, try mode=data or mode=mc")
  print("################################")
  print("")
  sys.exit(0)

if options.prescale:
  print("#################################")
  print("Including prescales in the Ntuple")
  print("#################################")
  print("")
else:
  print("#################################")
  print("_without_ prescales in the Ntuple")
  print("#################################")
  print("")

if options.trigger:
  print("################################")
  print("Including triggers in the Ntuple")
  print("################################")
  print("")
else:
  print("#################################")
  print("_without_ triggers in the Ntuple")
  print("#################################")
  print("")

if options.unmatching:
  print("#################################################################")
  print("leading leptons not associated with jets cone (unmatching = true)")
  print("#################################################################")
  print("")
else:
  print("##################################################################")
  print("leading leptons associated with the jets cone (unmatching = false)")
  print("##################################################################")
  print("")

if options.debugging:
  print("################")
  print("debugging active")
  print("################")
  print("")

#########################
#    General options    #
#########################

process.load("FWCore.MessageService.MessageLogger_cfi")
process.options   = cms.untracked.PSet(
    #wantSummary = cms.untracked.bool(True),
    #SkipEvent = cms.untracked.vstring('ProductNotFound'),
    allowUnscheduled = cms.untracked.bool(True),
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#########################
#      Input files      #
#########################

process.source = cms.Source("PoolSource",
				fileNames = cms.untracked.vstring(fileinput),
				#firstEvent = cms.untracked.uint32(0)
)

#########################
#        Triggers       #
#########################
process.load("PPSFramework.Skimmer.HLTFilter_cfi")
process.hltFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.hltFilter.HLTPaths = (triggerlist)

#########################
#      Preskimming      #
#########################
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag

'''
if mc:
 if options.year==2017:
  process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mc2017_realistic_v17') # 94X MC campaing
 if options.year==2018:
  process.GlobalTag = GlobalTag(process.GlobalTag, '102X_upgrade2018_realistic_v20') # Autumn18 MC campaing
else:
 process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v22') # global tag which includes PPS alignment and optics. Default: auto:run2_data
'''

process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v22') # global tag which includes PPS alignment and optics. Default: auto:run2_data

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

#########################
#     Proton Filter     #
#########################
process.protontagFilter = cms.EDFilter("ProtonTagFilter",
					debugging = cms.bool(options.debugging),
					singlePotMode = cms.bool(True), #if False, filter requires at least one proton per arm.
					includeProtonsReco = cms.bool(True),
				        protonSingleTag = cms.InputTag("ctppsProtons", "singleRP"),
					protonMultiTag = cms.InputTag("ctppsProtons", "multiRP"),
					pixelsppsTag = cms.InputTag('ctppsPixelLocalTracks'),
					timingppsTag = cms.InputTag('ctppsDiamondLocalTracks'),
					stripstotemTag = cms.InputTag('totemRPLocalTrackFitter')
					)

#########################
#       Analysis        #
#########################
process.load("PPSFramework.Skimmer.countsAnalyzer_cfi")
process.CounterNoCuts = process.countsAnalyzer.clone()
process.CounterAfterTrigger = process.countsAnalyzer.clone()
process.CounterTriggered = process.countsAnalyzer.clone()
process.CounterPPSTagged = process.countsAnalyzer.clone()

process.load("PPSFramework.Skimmer.MissingMassSearches_cfi")
process.missing_mass.physics = cms.string(options.physics) 
process.missing_mass.year = cms.string(options.year) 
process.missing_mass.enableMC = cms.bool(mc)
process.missing_mass.enableTrigger = cms.bool(options.trigger)
process.missing_mass.enablePrescales = cms.bool(options.prescale)
process.missing_mass.unmatching = cms.bool(options.unmatching) #if True, force lepton/jets to be separated
process.missing_mass.debugging = cms.bool(options.debugging)
process.missing_mass.includeMuons = cms.bool(True)
process.missing_mass.includeElectrons = cms.bool(True)
process.missing_mass.includeJets = cms.bool(True)
process.missing_mass.includePhotons = cms.bool(True)
process.missing_mass.includeMET = cms.bool(True)
process.missing_mass.includeVertices = cms.bool(True)
process.missing_mass.includeParticleFlow = cms.bool(True)
process.missing_mass.includeProtonsReco = cms.bool(True)
process.missing_mass.includePPSInfo = cms.bool(True)
process.missing_mass.triggersList = process.hltFilter.HLTPaths
process.missing_mass.tier = cms.string('miniAOD')
process.missing_mass.ppslite = cms.bool(True)
process.missing_mass.GenPartTag = cms.InputTag('prunedGenParticles')
process.missing_mass.GenMETTag = cms.InputTag('genMetTrue')
process.missing_mass.GenJetAlgoA = cms.InputTag('slimmedGenJets')
process.missing_mass.GenJetAlgoB = cms.InputTag('slimmedGenJetsAK8')
process.missing_mass.patJetAlgoA = cms.InputTag('slimmedJets')
process.missing_mass.patJetAlgoB = cms.InputTag('slimmedJetsAK8')
process.missing_mass.electronTag = cms.InputTag("slimmedElectrons")
process.missing_mass.muonTag = cms.InputTag("slimmedMuons")
process.missing_mass.pfTag = cms.InputTag('particleFlow')
process.missing_mass.packedTag = cms.InputTag('packedPFCandidates')
process.missing_mass.photonTag = cms.InputTag('slimmedPhotons')
process.missing_mass.patmetTag = cms.InputTag('slimmedMETs')
process.missing_mass.vertexTag = cms.InputTag('offlineSlimmedPrimaryVertices')
process.missing_mass.protonTag = cms.InputTag('ctppsLocalTrackLiteProducer')
process.missing_mass.protonSingleTag = cms.InputTag("ctppsProtons", "singleRP")
process.missing_mass.protonMultiTag = cms.InputTag("ctppsProtons", "multiRP")
process.missing_mass.pixelsppsTag = cms.InputTag('ctppsPixelLocalTracks')
process.missing_mass.timingppsTag = cms.InputTag('ctppsDiamondLocalTracks')
process.missing_mass.stripstotemTag = cms.InputTag('totemRPLocalTrackFitter')
process.missing_mass.caloTowerTag = cms.InputTag('towerMaker')
process.missing_mass.lhcInfoLabel = cms.string('')

# Proton, check optics
process.ctppsOpticsPlotter = cms.EDAnalyzer("CTPPSOpticsPlotter",
    outputFile = cms.string("output_optics.root")
)

# Proton Reconstruction Plotter
process.ctppsProtonReconstructionPlotter = cms.EDAnalyzer("CTPPSProtonReconstructionPlotter",
    tagTracks = cms.InputTag("ctppsLocalTrackLiteProducer"),
    tagRecoProtonsSingleRP = cms.InputTag("ctppsProtons", "singleRP"),
    tagRecoProtonsMultiRP = cms.InputTag("ctppsProtons", "multiRP"),

    rpId_45_F = cms.uint32(23),
    rpId_45_N = cms.uint32(3),
    rpId_56_N = cms.uint32(103),
    rpId_56_F = cms.uint32(123),

    outputFile = cms.string("outputproton.root")
)


# b-tagging for RECO
process.load('RecoBTag/Configuration/RecoBTag_cff')
process.load("JetMETCorrections.Type1MET.correctedMet_cff")

# prepare the output file
process.TFileService = cms.Service('TFileService',
    fileName = cms.string("output.root"),
    closeFileFast = cms.untracked.bool(True)
)

#####################################
# Proton Reconstruction, only for AOD
####################################

process.load("RecoCTPPS.Configuration.recoCTPPS_cff")

if forceDB:
	from CondCore.CondDB.CondDB_cfi import *
	process.CondDBOptics = CondDB.clone( connect = 'frontier://FrontierProd/CMS_CONDITIONS' )
	process.PoolDBESSourceOptics = cms.ESSource("PoolDBESSource",
	    process.CondDBOptics,
	    DumpStat = cms.untracked.bool(False),
	    toGet = cms.VPSet(cms.PSet(
        	record = cms.string('CTPPSOpticsRcd'),
	        tag = cms.string("PPSOpticalFunctions_offline_v5") # always check on https://twiki.cern.ch/twiki/bin/view/CMS/TaggedProtonsRecommendations
    	)),
	)
	process.esPreferDBFileOptics = cms.ESPrefer("PoolDBESSource","PoolDBESSourceOptics")

	from CondCore.CondDB.CondDB_cfi import *
	process.CondDBAlignment = CondDB.clone( connect = 'frontier://FrontierProd/CMS_CONDITIONS' )
	process.PoolDBESSourceAlignment = cms.ESSource("PoolDBESSource",
	    process.CondDBAlignment,
	    #timetype = cms.untracked.string('runnumber'),
	    toGet = cms.VPSet(cms.PSet(
        	record = cms.string('RPRealAlignmentRecord'),
	        tag = cms.string('CTPPSRPAlignment_real_offline_v4') # always check on https://twiki.cern.ch/twiki/bin/view/CMS/TaggedProtonsRecommendations
    	))
	)

	process.esPreferDBFileAlignment = cms.ESPrefer("PoolDBESSource",
							"PoolDBESSourceAlignment")

process.MiniAOD = cms.Sequence(
					process.CounterNoCuts
					*process.hltFilter
					*process.CounterTriggered
					*process.protontagFilter
					*process.CounterPPSTagged
					*process.missing_mass
                                        *process.ctppsOpticsPlotter
                                )


if not options.trigger:
   process.MiniAOD.remove(process.CounterTriggered)
   process.MiniAOD.remove(process.hltFilter)

if not options.debugging:
   process.MiniAOD.remove(process.ctppsOpticsPlotter)

if not options.ppstagging:
   process.MiniAOD.remove(process.CounterPPSTagged)
   process.MiniAOD.remove(process.protontagFilter)

process.p = cms.Path(process.MiniAOD)
