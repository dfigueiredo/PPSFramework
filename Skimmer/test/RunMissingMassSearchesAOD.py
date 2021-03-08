#flake8: noqa

#
# Forcing DB changes for Proton. The new optics is not yet into the GT
#
forceDB = True

import sys
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

# Setting Input Parameters from Line Command
options = VarParsing ('analysis')
options.register('Mode','Muon',VarParsing.multiplicity.singleton, VarParsing.varType.string,"Option to Run: Muon or Electron")
options.register('Year','2018',VarParsing.multiplicity.singleton, VarParsing.varType.string,"Year: by default, 2018")
options.register('MC',False,VarParsing.multiplicity.singleton, VarParsing.varType.bool,"If running over MC, this option must be True")
options.parseArguments()

# For Proton Reconstruction
from Configuration.StandardSequences.Eras import eras
process = cms.Process("MissMass",eras.ctpps_2016)

if options.Year != "2017" and options.Year !="2018":
  	  print("\n Please, use only years 2017 or 2018!\n")
	  sys.exit(0)

if options.Year == "2018" and options.Mode =="ZeroBias":
          print("\n Please, ZeroBias is not yet configured to run over 2018 data!\n")
          sys.exit(0)

if options.Mode == "Muon":
  print("")
  print("#########")
  print("Data Muon")
  print("#########")
  print("")
  # 2017, 2018, same path!
  triggerlist = 'HLT_IsoMu27_v*', 'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*', 'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v*', 'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*', 'HLT_DoubleMu43NoFiltersNoVtx_v*', 'HLT_DoubleMu48NoFiltersNoVtx_v*' 
  if options.Year == "2017":
	testfile = '/store/data/Run2017C/DoubleMuon/AOD/17Nov2017-v1/50001/2005DA54-85D3-E711-BF10-02163E011E35.root'
  if options.Year == "2018":
	testfile = '/store/data/Run2018D/DoubleMuon/AOD/PromptReco-v2/000/325/172/00000/583B14CE-796E-C54D-866C-9312930623ED.root'
elif options.Mode == "MC":
  options.MC=True
  options.Tier="AOD"
  print("")
  print("############")
  print("Data MC Muon")
  print("############")
  print("")
  # 2017, 2018, same path!
  triggerlist = 'HLT_IsoMu27_v*', 'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*', 'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v*', 'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*', 'HLT_DoubleMu43NoFiltersNoVtx_v*', 'HLT_DoubleMu48NoFiltersNoVtx_v*'
  print('\t\t-->Tier is AOD\n')
  testfile = '/store/mc/RunIISpring16DR80/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/PUFlat0to50_80X_mcRun2_asymptotic_2016_v3-v1/80000/F4BEDC9D-84F2-E511-8E98-001E67A40451.root'
elif options.Mode == "Electron":
  print("")
  print("#############")
  print("Data Electron")
  print("#############")
  print("")
  # 2017, 2018, same path!
  triggerlist = 'HLT_Ele27_WPTight_Gsf_v*', 'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*', 'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*', 'HLT_DoubleEle33_CaloIdL_MW_v*'
  if options.Year == "2017":
	testfile = '/store/data/Run2017B/DoubleEG/AOD/17Nov2017-v1/50002/F8BF871A-ECD9-E711-893C-14187763B750.root'
  if options.Year == "2018":
	testfile = '/store/data/Run2018A/EGamma/AOD/PromptReco-v1/000/315/322/00000/7A86BEAE-604C-E811-9A63-FA163EB308B9.root'
elif options.Mode == "Zerobias" or options.Mode == "ZEROBIAS":
  print("")
  print("########")
  print("ZeroBias")
  print("########")
  print("")
  triggerlist = 'HLT_ZeroBias_v*', 'HLT_ZeroBias_FirstBXAfterTrain_v*','HLT_ZeroBias_IsolatedBunches_v*','HLT_ZeroBias_Alignment_v*','HLT_ZeroBias_Beamspot_v*','HLT_L1UnpairedBunchBptxMinus_v*','HLT_L1UnpairedBunchBptxPlus_v*','HLT_Physics_v*','HLT_L1SingleMu_v*'
  if options.Year == "2017":
	testfile = '/store/data/Run2017E/ZeroBias/AOD/17Nov2017-v1/20000/103F179D-55D1-E711-AA3F-001E67792518.root'
  if options.Year == "2018":
        testfile = '/store/data/Run2018C/ZeroBias/AOD/17Sep2018-v1/00000/1101F169-58F3-3442-8432-6ACF16F56A47.root'
else:
  print("")
  print("#######################################################")
  print("Please, try Mode=Muon or Mode=Electron or Mode=ZeroBias")
  print("#######################################################")
  print("")
  sys.exit(0)

if options.MC:
  print("##############")
  print("Running for MC")
  print("##############")
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
				fileNames = cms.untracked.vstring(testfile),
				#firstEvent = cms.untracked.uint32(0)
)

#########################
#        Triggers       #
#########################
process.load("MissingMass.Skimmer.HLTFilter_cfi")
process.hltFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.hltFilter.HLTPaths = (triggerlist)

#########################
#      Preskimming      #
#########################
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data')
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v22') # global tag which includes PPS alignment and optics!
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")


#########################
#     PAT-ification     #
#########################
## Look at https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATTools#Core_Tools for more information

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('PATuple.root'),
    outputCommands = cms.untracked.vstring(
        'drop *',
       	'keep *_offline*PrimaryVertices*_*_*',
        'keep *_selectedPatMuons*_*_*',
       	'keep *_*lectron*_*_*',
        'keep *_selectedPatElectrons*_*_*',
       	'keep *_selectedPat*Photons*_*_*',
       	'keep *_selectedPatPhotons*_*_*',
        'keep *_selectedPatJets*_*_*',
        'keep *_*Photons*_*_*',
        'keep *_*MET*_*_*',
       	'keep *_*particleFlow*_*_*',
       	'keep *_*patJets*_*_*',
       	'keep *_*patPhotons*_*_*',
	'keep patPhotons*_*_*_*',
	'keep patMETs_patMETs__PAT'
    ),
)
       
from PhysicsTools.PatAlgos.tools.helpers import getPatAlgosToolsTask
patAlgosToolsTask = getPatAlgosToolsTask(process)

process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
patAlgosToolsTask.add(process.patCandidatesTask)

process.load("PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff")
patAlgosToolsTask.add(process.selectedPatCandidatesTask)

## DeepCSV meta discriminators (simple arithmethic on output probabilities)
process.load('RecoBTag.Combined.deepFlavour_cff')
patAlgosToolsTask.add(process.pfDeepCSVDiscriminatorsJetTags)
process.patJets.discriminatorSources.extend([
	cms.InputTag('pfDeepCSVDiscriminatorsJetTags:BvsAll' ),
	cms.InputTag('pfDeepCSVDiscriminatorsJetTags:CvsB'   ),
	cms.InputTag('pfDeepCSVDiscriminatorsJetTags:CvsL'   ),
])

noDeepFlavourDiscriminators = [x.value() for x in process.patJets.discriminatorSources if not "DeepFlavour" in x.value()]
from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection
addJetCollection(process, postfix   = "", labelName = 'AK8PFCHS', jetSource = cms.InputTag('ak8PFJetsCHS'), 
          jetCorrections = ('AK8PFchs', ['L1FastJet','L2Relative', 'L3Absolute'], 'None'),
          algo= 'ak8', rParam = 0.8, btagDiscriminators = noDeepFlavourDiscriminators
	)

from PhysicsTools.PatAlgos.tools.coreTools import runOnData
if not options.MC:
	runOnData( process )


#########################
#      E/Gamma ID       #
#########################
from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq,makeEgammaPATWithUserData
setupEgammaPostRecoSeq(process,isMiniAOD=False,era='2018-Prompt')

makeEgammaPATWithUserData(process,eleTag=cms.InputTag("selectedPatElectrons"),
                          phoTag=cms.InputTag("selectedPatPhotons"),
                          runVID=True,
                          runEnergyCorrections=True,
                          suffex="WithUserData",
                          era='2018-Prompt', force=True)

patAlgosToolsTask.add(process.egammaPostRecoPatUpdatorTask)

#########################
# Jet Box               #
#########################
from JMEAnalysis.JetToolbox.jetToolbox_cff import jetToolbox
jetToolbox( process, 'ak4', 'jetSequence', 'noOutput', PUMethod='CHS', addPruning=False,  dataTier='AOD', runOnMC=False, addQGTagger=True, addPUJetID=True)
jetToolbox( process, 'ak8', 'jetSequence', 'noOutput', PUMethod='CHS', addPruning=False,  dataTier='AOD', runOnMC=False)


#########################
#     Proton Filter     #
#########################
process.protontagFilter = cms.EDFilter("ProtonTagFilter",
                                        debugging = cms.bool(False),
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
process.load("MissingMass.Skimmer.countsAnalyzer_cfi")
process.CounterBeforeTrigger = process.countsAnalyzer.clone()
process.CounterAfterTrigger = process.countsAnalyzer.clone()
process.CounterBeforePPSTagging = process.countsAnalyzer.clone()
process.CounterAfterPPSTagging = process.countsAnalyzer.clone()

process.load("MissingMass.Skimmer.MissingMassSearches_cfi")
process.missing_mass.mode = cms.string(options.Mode) #Muons or Electrons or Zerobias
process.missing_mass.year = cms.string(options.Year)
process.missing_mass.debugging = cms.bool(False)
process.missing_mass.includeMuons = cms.bool(True)
process.missing_mass.includeElectrons = cms.bool(True)
process.missing_mass.includeJets = cms.bool(True)
process.missing_mass.includePhotons = cms.bool(True)
process.missing_mass.includeMET = cms.bool(True)
process.missing_mass.includeVertices = cms.bool(True)
process.missing_mass.includeParticleFlow = cms.bool(False)
process.missing_mass.includeProtonsReco = cms.bool(True)
process.missing_mass.includePPSInfo = cms.bool(True)
process.missing_mass.triggersList = process.hltFilter.HLTPaths
process.missing_mass.tier = cms.string('AOD')
process.missing_mass.ppslite = cms.bool(False)
process.missing_mass.patJetAlgoA = cms.InputTag('selectedPatJetsAK4PFCHS')
process.missing_mass.patJetAlgoB = cms.InputTag('selectedPatJetsAK8PFCHS')
process.missing_mass.electronTag = cms.InputTag("selectedPatElectronsWithUserData")
process.missing_mass.muonTag = cms.InputTag("patMuons")
process.missing_mass.pfTag = cms.InputTag('particleFlow')
process.missing_mass.packedTag = cms.InputTag('packedPFCandidates')
process.missing_mass.photonTag = cms.InputTag('selectedPatPhotonsWithUserData')
process.missing_mass.patmetTag = cms.InputTag('patMETs')
process.missing_mass.vertexTag = cms.InputTag('offlinePrimaryVertices')
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

# prepare the output file
process.TFileService = cms.Service('TFileService',
    fileName = cms.string('output.root'),
    closeFileFast = cms.untracked.bool(True)
)

process.printContent = cms.EDAnalyzer("EventContentAnalyzer",
    #should we print data? (sets to 'true' if verboseForModuleLabels has entries)
    verbose = cms.untracked.bool(False),
    #how much to indent when printing verbosely
    verboseIndentation = cms.untracked.string('  '),
    #string used at the beginning of all output of this module
    indentation = cms.untracked.string('++'),
    #data from which modules to print (all if empty)
    verboseForModuleLabels = cms.untracked.vstring(),
    # which data from which module should we get without printing
    getDataForModuleLabels = cms.untracked.vstring(),
    #should we get data? (sets to 'true' if getDataFormModuleLabels has entries)
    getData = cms.untracked.bool(False)
)

process.load("JetMETCorrections.Type1MET.correctedMet_cff")

import HLTrigger.HLTcore.hltEventAnalyzerAOD_cfi
process.hltAOD = HLTrigger.HLTcore.hltEventAnalyzerAOD_cfi.hltEventAnalyzerAOD.clone()
process.hltAOD.processName = cms.string("HLT")
process.hltAOD.triggerResults = cms.InputTag("TriggerResults","","HLT")
process.hltAOD.triggerEvent = cms.InputTag("hltTriggerSummaryAOD","","HLT")

#######################
# Proton Reconstruction
#######################

process.load("RecoCTPPS.Configuration.recoCTPPS_cff")
if forceDB:
	from CondCore.CondDB.CondDB_cfi import *
	process.CondDBOptics = CondDB.clone( connect = 'frontier://FrontierProd/CMS_CONDITIONS' )
	process.PoolDBESSourceOptics = cms.ESSource("PoolDBESSource",
	    process.CondDBOptics,
	    DumpStat = cms.untracked.bool(False),
	    toGet = cms.VPSet(cms.PSet(
        	record = cms.string('CTPPSOpticsRcd'),
	        tag = cms.string("PPSOpticalFunctions_offline_v5")
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
	        tag = cms.string('CTPPSRPAlignment_real_offline_v4')
    	))
	)

	process.esPreferDBFileAlignment = cms.ESPrefer("PoolDBESSource",
							"PoolDBESSourceAlignment")



#
# CMSSW PATHS
#

#
# OPTIONAL: for debugging purposes
# 
# process.ctppsProtonReconstructionPlotter: Analyzer for proton reconstruction.
# process.ctppsOpticsPlotter: Analyzer for optics.
# process.printContent: It shows all the modules loaded.
#
process.calibratedDiamondRecHits = process.ctppsDiamondRecHits.clone()
process.ctppsDiamondLocalTracks.recHitsTag = cms.InputTag('calibratedDiamondRecHits')

process.AOD = cms.Sequence(
		 	process.totemRPUVPatternFinder*
                    	process.totemRPLocalTrackFitter*
                   	process.calibratedDiamondRecHits*
                    	process.ctppsDiamondLocalTracks*
                    	process.ctppsPixelLocalTracks*
                    	process.ctppsLocalTrackLiteProducer*
                    	process.ctppsProtons*
                    	process.CounterBeforeTrigger*
                    	process.hltFilter*
                    	process.CounterBeforePPSTagging*
                    	process.protontagFilter*
                    	process.CounterAfterPPSTagging*
                    	process.btagging*
                    	process.egammaPostRecoSeq*
                    	process.missing_mass
                    	#process.ctppsProtonReconstructionPlotter*
                    	#process.ctppsOpticsPlotter
                    	#*process.hltAOD
                    	#*process.printContent
                        )

process.p = cms.Path(process.AOD)
process.outpath = cms.EndPath(patAlgosToolsTask)
