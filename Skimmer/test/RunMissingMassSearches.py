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
options.register('Mode','MC_Muon',VarParsing.multiplicity.singleton, VarParsing.varType.string,"Option to Run: Muon, MC_Muon, Electron, MC_Electron, ZeroBias or Emu")
options.register('Year','2017',VarParsing.multiplicity.singleton, VarParsing.varType.string,"Year: by default, 2018")
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
  triggerlist = 'HLT_IsoMu27_v*', 'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*', 'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v*', 'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*', 'HLT_DoubleMu43NoFiltersNoVtx_v*', 'HLT_DoubleMu48NoFiltersNoVtx_v*', 'HLT_BTagMu_AK4DiJet20_Mu5_v*', 'HLT_BTagMu_AK4DiJet40_Mu5_v*', 'HLT_BTagMu_AK4DiJet70_Mu5_v*', 'HLT_BTagMu_AK4DiJet110_Mu5_v*', 'HLT_BTagMu_AK4DiJet170_Mu5_v*', 'HLT_BTagMu_AK4Jet300_Mu5_v*', 'HLT_BTagMu_AK8DiJet170_Mu5_v*', 'HLT_BTagMu_AK8Jet300_Mu5_v*', 'HLT_PFHT380_SixJet32_DoubleBTagCSV_p075_v*', 'HLT_PFHT430_SixJet40_BTagCSV_p080_v*'
  print('\t\t-->Tier is miniAOD\n')
  if options.Year == "2017":
	  #testfile = '/store/data/Run2017B/DoubleMuon/MINIAOD/09Aug2019_UL2017-v1/260000/F1B17EC1-36AF-6340-B359-774C7D88948D.root'
	  testfile = '/store/data/Run2017D/BTagCSV/MINIAOD/09Aug2019_UL2017-v1/60000/CAEC3BF9-E990-AD41-A55B-B28269E242C1.root'
  if options.Year == "2018":
	  testfile = '/store/data/Run2018C/DoubleMuon/MINIAOD/ForValUL2018-v1/230000/1A9817AE-B1F4-A746-9460-34B48037B519.root'

elif options.Mode == "MC_Muon":
  options.MC=True
  print("")
  print("#######")
  print("MC Muon")
  print("#######")
  print("")
  # 2017, 2018, same path!
  triggerlist = 'HLT_IsoMu27_v*', 'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*', 'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v*', 'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*', 'HLT_DoubleMu43NoFiltersNoVtx_v*', 'HLT_DoubleMu48NoFiltersNoVtx_v*' 
  print('\t\t-->Tier is miniAOD\n')
  #testfile = 'root://cms-xrd-global.cern.ch////store/user/dmf/HXToyMC-13TeV-miniAOD_2021-01-05_UTC12-28-09/HXToyMC_NONE_PostTS2_850_130/HXToyMC-13TeV-miniAOD_2021-01-05_UTC12-28-09/210105_112841/0000/output_130.root'
  #testfile = 'file:////afs/cern.ch/user/d/dmf/private/work/private/CMSPhysicsAnalysis/PrivateMCProduction/2019_production_13TeV/EXO-19-009/PPS_MC_Production/CMSSW_10_6_8_patch1/src/ProtonMonteCarloProduction/Grid/step4.root'
  testfile = 'file:////afs/cern.ch/user/d/dmf/private/work/private/CMSPhysicsAnalysis/PrivateMCProduction/2019_production_13TeV/EXO-19-009/PPS_MC_Production/CMSSW_10_6_8_patch1/src/ProtonMonteCarloProduction/Grid/step4_not_stable_b.root'
  #testfile = 'root://cms-xrd-global.cern.ch///store/user/dmf/PYTHIA8-SD-top-pt70-13TeV-miniAOD-xangle130-PreTS2-miniAOD_PreTS2_130_2020-12-31_UTC22-35-07/Pythia-sd-top-pt70/PYTHIA8-SD-top-pt70-13TeV-miniAOD-xangle130-PreTS2-miniAOD_PreTS2_130_2020-12-31_UTC22-35-07/201231_213522/0000/output_127.root'
  #testfile = 'root://cms-xrd-global.cern.ch///store/user/dmf/Pomwig-Zmumu-Plus-13TeV-xangle150-PostTS2-miniAOD_PostTS2_150_2020-03-12_UTC19-46-43/Pomwig-Zmumu-Plus-GEN/Pomwig-Zmumu-Plus-13TeV-xangle150-PostTS2-miniAOD_PostTS2_150_2020-03-12_UTC19-46-43/200312_184700/0000/output_15.root'
  #testfile = 'root://cmsxrootd.fnal.gov//store/user/dmf/PYTHIA8-SD-top-13TeV-miniAOD-xangle120-PreTS2-miniAOD_PreTS2_120_2020-05-08_UTC15-06-33/PYTHIA8-SD-top-13TeV/PYTHIA8-SD-top-13TeV-miniAOD-xangle120-PreTS2-miniAOD_PreTS2_120_2020-05-08_UTC15-06-33/200508_130659/0004/output_4214.root'
#'file:/eos/user/l/lpagliai/PPZX_Simulation/MiniAOD/Muons/versionLorenzo/Z/m_X_950/xangle_150/2017_postTS2/output.root'
#'root://xrootd-cms.infn.it//store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/90000/FEB3954C-4942-E811-8A09-008CFAC91A38.root'
#'file:/eos/user/l/lpagliai/PPZX_Simulation/MiniAOD/Muons/versionLorenzo/Z/m_X_950/xangle_140/2017_preTS2/output.root' 
#'root://xrootd-cms.infn.it//store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/90000/FEB3954C-4942-E811-8A09-008CFAC91A38.root'
#'file:/afs/cern.ch/work/l/lpagliai/public/PUmixing/Single_diffraction/output/versionLorenzo/Z/Minus/xangle_150/2017_postTS2/miniAOD.root'
#'file:/eos/user/l/lpagliai/PPZX_Simulation/MiniAOD/Muons/versionLorenzo/Z/m_X_950/xangle_140/2017_preTS2/output.root' 
#'file:/afs/cern.ch/work/l/lpagliai/public/DYSimulator/CMSSW_10_6_0/src/MissingMass/Skimmer/test/miniAOD_withProtons.root'
#'/store/user/lviliani/POMWIG_SingleDiffractiveZmumu_13TeV_GEN_JanCuts/POMWIG_SingleDiffractiveZmumu_13TeV_GEN_JanCuts/190318_163531/0000/POMWIG_SingleDiffractiveZmumu_13TeV_JanCuts_15.root'
#'file:/afs/cern.ch/work/l/lpagliai/public/PUmixing/Ntuples/Muons/versionLorenzo/Z/m_X_950/xangle_120/2017_postTS2/miniAOD.root'
 #'file:/afs/cern.ch/work/l/lpagliai/public/PUmixing/production_106xSummer19UL17/output/versionLorenzo/Z/m_X_950/xangle_120/2017_postTS2/miniAOD.root'
  #'file:/eos/user/l/lpagliai/PPZX_Simulation/Ntuples/Muons/versionLorenzo/Z/m_X_950/xangle_120/2017_preTS2/ntuple.root'

elif options.Mode == "Electron":
  print("")
  print("#############")
  print("Data Electron")
  print("#############")
  print("")
  # 2017, 2018, same path!
  triggerlist = 'HLT_Ele27_WPTight_Gsf_v*', 'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*', 'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*', 'HLT_DoubleEle33_CaloIdL_MW_v*', 'HLT_BTagMu_AK4DiJet20_Mu5_v*', 'HLT_BTagMu_AK4DiJet40_Mu5_v*', 'HLT_BTagMu_AK4DiJet70_Mu5_v*', 'HLT_BTagMu_AK4DiJet110_Mu5_v*', 'HLT_BTagMu_AK4DiJet170_Mu5_v*', 'HLT_BTagMu_AK4Jet300_Mu5_v*', 'HLT_BTagMu_AK8DiJet170_Mu5_v*', 'HLT_BTagMu_AK8Jet300_Mu5_v*', 'HLT_PFHT380_SixJet32_DoubleBTagCSV_p075_v*', 'HLT_PFHT430_SixJet40_BTagCSV_p080_v*','HLT_PFMET100_PFMHT100_IDTight_CaloBTagCSV_3p1_v*', 'HLT_PFMET110_PFMHT110_IDTight_CaloBTagCSV_3p1_v*', 'HLT_PFMET120_PFMHT120_IDTight_CaloBTagCSV_3p1_v*', 'HLT_PFMET130_PFMHT130_IDTight_CaloBTagCSV_3p1_v*','HLT_PFMET140_PFMHT140_IDTight_CaloBTagCSV_3p1_v*', 'HLT_Ele15_IsoVVVL_PFHT450_CaloBTagCSV_4p5_v*', 'DST_HT250_CaloBTagScouting_v*'

  print('\t\t-->Tier is miniAOD\n')
  if options.Year == "2017":
  	testfile = '/store/data/Run2017C/DoubleEG/MINIAOD/09Aug2019_UL2017-v1/50000/F696E155-EF11-CB48-82CA-96EA48359D54.root'
  if options.Year == "2018":
	testfile = '/store/data/Run2018C/EGamma/MINIAOD/ForValUL2018-v2/230000/5C5DBBE7-A5E1-E746-835A-D33E0C37071F.root'

elif options.Mode == "MC_Electron":
  options.MC=True
  print("")
  print("###########")
  print("MC Electron")
  print("###########")
  print("")
  # 2017, 2018, same path!
  triggerlist = 'HLT_Ele27_WPTight_Gsf_v*', 'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*', 'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*', 'HLT_DoubleEle33_CaloIdL_MW_v*' 
  print('\t\t-->Tier is miniAOD\n')
  #testfile = 'root://cms-xrd-global.cern.ch///store/user/lpagliai/ZXToyMC-MINIAOD_Electron_PreTS2_950_120_2020-03-13_UTC14-30-11/ZXToyMC-GEN_Electron_PreTS2_950_120/ZXToyMC-MINIAOD_Electron_PreTS2_950_120_2020-03-13_UTC14-30-11/200313_133059/0000/output_99.root'
  testfile = 'root://cms-xrd-global.cern.ch////store/user/dmf/HXToyMC-13TeV-miniAOD_2021-01-05_UTC12-28-09/HXToyMC_NONE_PostTS2_850_130/HXToyMC-13TeV-miniAOD_2021-01-05_UTC12-28-09/210105_112841/0000/output_130.root'
#'root://xrootd-cms.infn.it//store/user/dmf/Pomwig-Zee-Plus-13TeV-xangle140-PreTS2-miniAOD_PreTS2_140_2020-03-12_UTC19-39-40/Pomwig-Zee-Plus-GEN/Pomwig-Zee-Plus-13TeV-xangle140-PreTS2-miniAOD_PreTS2_140_2020-03-12_UTC19-39-40/200312_183958/0000/output_98.root'
#'file:/eos/user/l/lpagliai/PPZX_Simulation/MiniAOD/Electrons/versionLorenzo/Z/m_X_950/xangle_150/2017_postTS2/output.root'
#'root://xrootd-cms.infn.it//store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/90000/FEB3954C-4942-E811-8A09-008CFAC91A38.root'


elif options.Mode == "ZeroBias" or options.Mode == "ZEROBIAS":
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

elif options.Mode == "EMu":
  print("")
  print("#########")
  print("Data EMu")
  print("#########")
  print("")
  # 2017, 2018, same path!
  triggerlist = 'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*', 'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*', 'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*' #OTHER TO BE INSERTED?
  print('\t\t-->Tier is miniAOD\n')
  if options.Year == "2017":
          testfile = '/store/data/Run2017E/MuonEG/MINIAOD/09Aug2019_UL2017-v1/30000/D68314AC-E5B8-BB43-A0F9-86F599A8DA82.root'
  if options.Year == "2018":
          testfile = '/store/data/Run2018A/MuonEG/MINIAOD/12Nov2019_UL2018-v1/00000/51E69F77-C0F4-B942-B704-3826C9BA0F26.root'

elif options.Mode == "MC_EMu":
  options.MC=True
  print("")
  print("################")
  print("MC Electron/Muon")
  print("################")
  print("")
  # 2017, 2018, same path!
  triggerlist = 'HLT_Ele27_WPTight_Gsf_v*', 'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*', 'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*', 'HLT_DoubleEle33_CaloIdL_MW_v*'
  print('\t\t-->Tier is miniAOD\n')
  testfile = 'root://cms-xrd-global.cern.ch///store/user/dmf/PYTHIA8-SD-top-pt70-13TeV-miniAOD-xangle130-PreTS2-miniAOD_PreTS2_130_2020-12-31_UTC22-35-07/Pythia-sd-top-pt70/PYTHIA8-SD-top-pt70-13TeV-miniAOD-xangle130-PreTS2-miniAOD_PreTS2_130_2020-12-31_UTC22-35-07/201231_213522/0000/output_127.root'
  #testfile = 'root://cmsxrootd.fnal.gov//store/user/dmf/PYTHIA8-SD-top-13TeV-miniAOD-xangle120-PreTS2-miniAOD_PreTS2_120_2020-05-08_UTC15-06-33/PYTHIA8-SD-top-13TeV/PYTHIA8-SD-top-13TeV-miniAOD-xangle120-PreTS2-miniAOD_PreTS2_120_2020-05-08_UTC15-06-33/200508_130659/0004/output_4214.root'
  #testfile = '/store/group/phys_exotica/PPS-Exo/Test/miniAOD_test.root'
#'root://xrootd-cms.infn.it//store/user/dmf/Pomwig-Zee-Plus-13TeV-xangle140-PreTS2-miniAOD_PreTS2_140_2020-03-12_UTC19-39-40/Pomwig-Zee-Plus-GEN/Pomwig-Zee-Plus-13TeV-xangle140-PreTS2-miniAOD_PreTS2_140_2020-03-12_UTC19-39-40/200312_183958/0000/output_98.root'
#'file:/eos/user/l/lpagliai/PPZX_Simulation/MiniAOD/Electrons/versionLorenzo/Z/m_X_950/xangle_150/2017_postTS2/output.root'
#'root://xrootd-cms.infn.it//store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/90000/FEB3954C-4942-E811-8A09-008CFAC91A38.root'

else:
  print("")
  print("########################################################################################################")
  print("Please, try Mode=Muon or Mode=MC_Muon or Mode=Electron or Mode=MC_Electron or Mode=ZeroBias or Mode =EMu")
  print("########################################################################################################")
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
# Jet Box               #
#########################
from JMEAnalysis.JetToolbox.jetToolbox_cff import jetToolbox
jetToolbox( process, 'ak4', 'jetSequence', 'noOutput', PUMethod='CHS', addPruning=False,  dataTier='miniAOD', runOnMC=options.MC, addQGTagger=True, addPUJetID=True)
jetToolbox( process, 'ak8', 'jetSequence', 'noOutput', PUMethod='CHS', addPruning=False,  dataTier='miniAOD', runOnMC=options.MC)

#########################
#     Proton Filter     #
#########################
process.protontagFilter = cms.EDFilter("ProtonTagFilter",
					debugging = cms.bool(False),
					singlePotMode = cms.bool(True),
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
process.missing_mass.mode = cms.string(options.Mode) 
process.missing_mass.year = cms.string(options.Year) 
process.missing_mass.enableMC = cms.bool(options.MC)
process.missing_mass.enableTrigger = cms.bool(True)
process.missing_mass.matching = cms.bool(False) #if True, force lepton/jets to be separated
process.missing_mass.debugging = cms.bool(False)
process.missing_mass.includeMuons = cms.bool(True)
process.missing_mass.includeElectrons = cms.bool(True)
process.missing_mass.includeJets = cms.bool(True)
process.missing_mass.includePhotons = cms.bool(True)
process.missing_mass.includeMET = cms.bool(True)
process.missing_mass.includeVertices = cms.bool(True)
process.missing_mass.includeParticleFlow = cms.bool(True)
process.missing_mass.includeProtonsReco = cms.bool(True)
process.missing_mass.includePPSInfo = cms.bool(True)
print(process.hltFilter.HLTPaths)
process.missing_mass.triggersList = process.hltFilter.HLTPaths
#process.missing_mass.triggersList = (triggerlist)
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

# prepare the output file
process.TFileService = cms.Service('TFileService',
    fileName = cms.string("output.root"),
    closeFileFast = cms.untracked.bool(True)
)

process.printContent = cms.EDAnalyzer("EventContentAnalyzer",
    #should we print data? (sets to 'true' if verboseForModuleLabels has entries)
    verbose = cms.untracked.bool(True),
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

#######################
# Proton Reconstruction
#######################

#
# Only for AOD. Miniaod has Proton Reconstruction already done!
#

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


process.MiniAOD = cms.Sequence(
					process.CounterBeforeTrigger
					*process.hltFilter
					*process.CounterBeforePPSTagging
					*process.protontagFilter
					*process.CounterAfterPPSTagging
					*process.missing_mass
					#*process.ctppsOpticsPlotter
                                )

process.p = cms.Path(process.MiniAOD)
