#flake8: noqa

import sys
import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.Eras.Modifier_ctpps_2016_cff import ctpps_2016
from Configuration.Eras.Modifier_ctpps_2017_cff import ctpps_2017
from Configuration.Eras.Modifier_ctpps_2018_cff import ctpps_2018
from Configuration.ProcessModifiers.run2_miniAOD_UL_cff import run2_miniAOD_UL

# Setting Input Parameters from Line Command
options = VarParsing ('analysis')
options.register('physics','muon',VarParsing.multiplicity.singleton, VarParsing.varType.string,"physics(string): muon, electron, zerobias or emu. Default is muon.")
options.register('mode','mc',VarParsing.multiplicity.singleton, VarParsing.varType.string,"mode(string): data or mc. Default is data.")
options.register('era','C',VarParsing.multiplicity.singleton, VarParsing.varType.string,"era(string): B, C, D, E or F.")
options.register('year','2018',VarParsing.multiplicity.singleton, VarParsing.varType.string,"year(string): 2017 or 2018. Default is 2017.")
options.register('prescale',False,VarParsing.multiplicity.singleton, VarParsing.varType.bool,"prescale(bool): enable/disable trigger prescales values. Default is False (it works only for data).")
options.register('trigger',False,VarParsing.multiplicity.singleton, VarParsing.varType.bool,"trigger(bool): enable/disable the trigger. Default is True.")
options.register('discriminator',True,VarParsing.multiplicity.singleton, VarParsing.varType.bool,"discriminator(bool): enable/disable the ordering of the jets by btagging discriminator. Default is True.")
options.register('unmatching',False,VarParsing.multiplicity.singleton, VarParsing.varType.bool,"unmatching(bool): \n\t\tTrue: leptons and jets not associated in the same cone.\n\t\tFalse: otherwise\n. Default is False.")
options.register('debugging',False,VarParsing.multiplicity.singleton, VarParsing.varType.bool,"debugging(bool): \n\t\tTrue: enabling debugging.\n\t\t: Default is False.")
options.register('verbosity',0,VarParsing.multiplicity.singleton, VarParsing.varType.int,"verbosity(int): \n\t\tDebugging verbosity. Default is 0.")
options.parseArguments()

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

if options.year != "2017" and options.year !="2018":
          print("\n Please, use only years 2017 or 2018!\n")
          sys.exit(0)

if options.year == "2018" and options.mode =="zerobias":
          print("\n Please, zerobias is not yet configured to run over 2018 data!\n")
          sys.exit(0)

# Input file
#-----------
#fileinput = [
#	'/store/user/dmf/TOY-HX900-MC-eraC-2018-miniAOD-PPS-13TEV/crab_dmf_2022-01-26_UTC15-50-31/TOY-HX-MC-LHE-2018-13TEV_2022-01-13_UTC16-41-34/TOY-HX900-MC-eraC-2018-miniAOD-PPS-13TEV_2022-01-26_UTC15-50-31/220126_145046/0000/output_110.root',
#	'/store/user/dmf/TOY-HX900-MC-eraC-2018-miniAOD-PPS-13TEV/crab_dmf_2022-01-26_UTC15-50-31/TOY-HX-MC-LHE-2018-13TEV_2022-01-13_UTC16-41-34/TOY-HX900-MC-eraC-2018-miniAOD-PPS-13TEV_2022-01-26_UTC15-50-31/220126_145046/0000/output_111.root',
#'/store/user/dmf/TOY-HX900-MC-eraC-2018-miniAOD-PPS-13TEV/crab_dmf_2022-01-26_UTC15-50-31/TOY-HX-MC-LHE-2018-13TEV_2022-01-13_UTC16-41-34/TOY-HX900-MC-eraC-2018-miniAOD-PPS-13TEV_2022-01-26_UTC15-50-31/220126_145046/0000/output_112.root'
#	]

#fileinput = '/store/user/dmf/PYTHIA8-SD-TOP-MINIAOD-eraB-2017-withpps-13TEV/crab_dmf_2021-06-25_UTC09-48-16/PYTHIA8-SD-TOP-GEN/PYTHIA8-SD-TOP-MINIAOD-eraB-2017-withpps-13TEV/210625_074911/0000/output_101.root'
#fileinput = '/store/user/dmf/PYTHIA8-SD-TOP-MINIAOD-eraD-2017-withpps-13TEV/crab_dmf_2021-07-12_UTC18-52-49/PYTHIA8-SD-TOP-GEN/PYTHIA8-SD-TOP-MINIAOD-eraD-2017-withpps-13TEV/210712_165318/0000/output_1.root'
#fileinput = '/store/data/Run2017D/BTagCSV/MINIAOD/09Aug2019_UL2017-v1/60000/CAEC3BF9-E990-AD41-A55B-B28269E242C1.root'
fileinput = '/store/data/Run2018C/DoubleMuon/MINIAOD/ForValUL2018-v1/230000/1A9817AE-B1F4-A746-9460-34B48037B519.root'
#fileinput = '/store/user/dmf/PYTHIA8-SD-TOP-MINIAOD-eraD-2017-withpps-13TEV/crab_dmf_2021-08-09_UTC16-52-30/PYTHIA8-SD-TOP-GEN/PYTHIA8-SD-TOP-MINIAOD-eraD-2017-withpps-13TEV_2021-08-09_UTC16-52-30/210809_145244/0000/output_138.root'
#fileinput = '/store/user/dmf/ZXToyMC-13TeV-miniAOD_2020-03-23_UTC13-55-30/ZXToyMC-GEN_Muon_PostTS2_950_150/ZXToyMC-13TeV-miniAOD_2020-03-23_UTC13-55-30/200323_125619/0000/output_14.root'
#fileinput = '/store/user/dmf/ZXToyMC-13TeV-miniAOD_2020-03-23_UTC13-55-30/ZXToyMC-GEN_Muon_PostTS2_950_150/ZXToyMC-13TeV-miniAOD_2020-03-23_UTC13-55-30/200323_125619/0000/output_53.root'


#-----------------------------
# PPS Direct Sim, only for MC
#-----------------------------

if mc:
   if options.year == "2016" and (options.era == "B" or options.era == "C" or options.era == "G"):
      process = cms.Process('Analysis', ctpps_2016,run2_miniAOD_UL)
      from Validation.CTPPS.simu_config.year_2016_preTS2_cff import *
   if options.year == "2016" and (options.era == "H"):
      process = cms.Process('Analysis', ctpps_2016,run2_miniAOD_UL)
      from Validation.CTPPS.simu_config.year_2016_postTS2_cff import *
   if options.year == "2017" and (options.era == "B" or options.era == "C" or options.era == "D"):
      process = cms.Process('Analysis', ctpps_2017,run2_miniAOD_UL)
      from Validation.CTPPS.simu_config.year_2017_preTS2_cff import *
   if options.year == "2017" and (options.era == "E" or options.era == "F"):
      process = cms.Process('Analysis', ctpps_2017,run2_miniAOD_UL)
      from Validation.CTPPS.simu_config.year_2017_postTS2_cff import *
   if options.year == "2018" and (options.era == "A" or options.era == "B1"):
      process = cms.Process('Analysis', ctpps_2018,run2_miniAOD_UL)
      from Validation.CTPPS.simu_config.year_2018_cff import *
   if options.year == "2018" and (options.era == "B2" or options.era == "C" or options.era == "D1"):
      process = cms.Process('Analysis', ctpps_2018,run2_miniAOD_UL)
      from Validation.CTPPS.simu_config.year_2018_cff import *
   if options.year == "2018" and (options.era == "D2"):
      process = cms.Process('Analysis', ctpps_2018,run2_miniAOD_UL)
      from Validation.CTPPS.simu_config.year_2018_cff import *

   if options.year == "2016" and (options.era == "B" or options.era == "C" or options.era == "G"):
      process.load("Validation.CTPPS.simu_config.year_2016_preTS2_cff")
   if options.year == "2016" and (options.era == "H"):
      process.load("Validation.CTPPS.simu_config.year_2016_postTS2_cff")
   if options.year == "2017" and (options.era == "B" or options.era == "C" or options.era == "D"):
      process.load("Validation.CTPPS.simu_config.year_2017_preTS2_cff")
   if options.year == "2017" and (options.era == "E" or options.era == "F"):
      process.load("Validation.CTPPS.simu_config.year_2017_postTS2_cff")
   if options.year == "2018" and (options.era == "A" or options.era == "B1"):
      process.load("Validation.CTPPS.simu_config.year_2018_cff")
      #process.ctppsRPAlignmentCorrectionsDataESSourceXML.MisalignedFiles = ["../data/2018_preTS1.xml"]
      process.ctppsRPAlignmentCorrectionsDataESSourceXML.MisalignedFiles = ["PPSFramework/Skimmer/data/2018_preTS1.xml"]
      #process.ctppsRPAlignmentCorrectionsDataESSourceXML.RealFiles = ["../data/2018_preTS1.xml"]
      process.ctppsRPAlignmentCorrectionsDataESSourceXML.RealFiles = ["PPSFramework/Skimmer/data/2018_preTS1.xml"]
   if options.year == "2018" and (options.era == "B2" or options.era == "C" or options.era == "D1"):
      process.load("Validation.CTPPS.simu_config.year_2018_cff")
      #process.ctppsRPAlignmentCorrectionsDataESSourceXML.MisalignedFiles = ["../2018_TS1_TS2.xml"]
      process.ctppsRPAlignmentCorrectionsDataESSourceXML.MisalignedFiles = ["PPSFramework/Skimmer/data/2018_TS1_TS2.xml"]
      #process.ctppsRPAlignmentCorrectionsDataESSourceXML.RealFiles = ["..data/2018_TS1_TS2.xml"]
      process.ctppsRPAlignmentCorrectionsDataESSourceXML.RealFiles = ["PPSFramework/Skimmer/data/2018_TS1_TS2.xml"]
   if options.year == "2018" and (options.era == "D2"):
      process.load("Validation.CTPPS.simu_config.year_2018_cff")
      #process.ctppsRPAlignmentCorrectionsDataESSourceXML.MisalignedFiles = ["..data/2018_postTS2.xml"]
      process.ctppsRPAlignmentCorrectionsDataESSourceXML.MisalignedFiles = ["PPSFramework/Skimmer/data/2018_postTS2.xml"]
      #process.ctppsRPAlignmentCorrectionsDataESSourceXML.RealFiles = ["..data/2018_postTS2.xml"]
      process.ctppsRPAlignmentCorrectionsDataESSourceXML.RealFiles = ["PPSFramework/Skimmer/data/2018_postTS2.xml"]

   if options.year == "2017":
      process.ctppsOpticalFunctionsESSource.configuration[0].opticalFunctions = cms.VPSet(
        cms.PSet( xangle = cms.double(120), fileName = cms.FileInPath("CalibPPS/ESProducers/data/optical_functions/2017/version5tim/120urad.root") ),
        cms.PSet( xangle = cms.double(130), fileName = cms.FileInPath("CalibPPS/ESProducers/data/optical_functions/2017/version5tim/130urad.root") ),
        cms.PSet( xangle = cms.double(140), fileName = cms.FileInPath("CalibPPS/ESProducers/data/optical_functions/2017/version5tim/140urad.root") )
      )
   if options.year == "2018":
      process.ctppsOpticalFunctionsESSource.configuration[0].opticalFunctions = cms.VPSet(
        cms.PSet( xangle = cms.double(120), fileName = cms.FileInPath("CalibPPS/ESProducers/data/optical_functions/2018/version6/120urad.root") ),
        cms.PSet( xangle = cms.double(130), fileName = cms.FileInPath("CalibPPS/ESProducers/data/optical_functions/2018/version6/130urad.root") ),
        cms.PSet( xangle = cms.double(140), fileName = cms.FileInPath("CalibPPS/ESProducers/data/optical_functions/2018/version6/140urad.root") )
      )

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

   process.ctppsBeamParametersESSource.vtxStddevX = 0
   process.ctppsBeamParametersESSource.vtxStddevY = 0
   process.ctppsBeamParametersESSource.vtxStddevZ = 0

   process.load("CalibPPS.ESProducers.ctppsLHCInfoRandomXangleESSource_cfi")
   process.ctppsLHCInfoRandomXangleESSource.generateEveryNEvents = 1
   if options.year == "2016":
      process.ctppsBeamParametersESSource.vtxOffsetX45 = -1.048 * 1E-1
      process.ctppsBeamParametersESSource.vtxOffsetY45 = -1.686 * 1E-1
      process.ctppsBeamParametersESSource.vtxOffsetZ45 = +10.04 * 1E-1
      #process.ctppsLHCInfoRandomXangleESSource.xangleHistogramFile = "../data/CrossingAngles2016.root"
      process.ctppsLHCInfoRandomXangleESSource.xangleHistogramFile = "CrossingAngles2016.root"
   if options.year == "2017":
      process.ctppsBeamParametersESSource.vtxOffsetX45 = +0.24793  * 1E-1
      process.ctppsBeamParametersESSource.vtxOffsetY45 = -0.692861 * 1E-1
      process.ctppsBeamParametersESSource.vtxOffsetZ45 = -7.89895  * 1E-1
      #process.ctppsLHCInfoRandomXangleESSource.xangleHistogramFile = "../data/CrossingAngles2017.root"
      process.ctppsLHCInfoRandomXangleESSource.xangleHistogramFile = "CrossingAngles2017.root"
   if options.year == "2018":
      process.ctppsBeamParametersESSource.vtxOffsetX45 = -0.1078 * 1E-1
      process.ctppsBeamParametersESSource.vtxOffsetY45 = -0.4189 * 1E-1
      process.ctppsBeamParametersESSource.vtxOffsetZ45 = -0.2488 * 1E-1
      #process.ctppsLHCInfoRandomXangleESSource.xangleHistogramFile = "../data/CrossingAngles2018.root"
      process.ctppsLHCInfoRandomXangleESSource.xangleHistogramFile = "CrossingAngles2018.root"

   process.ctppsLHCInfoRandomXangleESSource.xangleHistogramObject = "hxang"
   process.ctppsLHCInfoRandomXangleESSource.beamEnergy = 6500.
   process.ctppsLHCInfoRandomXangleESSource.betaStar = 0.40
   process.esPreferLHCInfo = cms.ESPrefer("CTPPSLHCInfoRandomXangleESSource", "ctppsLHCInfoRandomXangleESSource")

   process.beamDivergenceVtxGenerator.src = cms.InputTag("")
   process.beamDivergenceVtxGenerator.srcGenParticle = cms.VInputTag(
      cms.InputTag("prunedGenParticles")
   )

else:
   process = cms.Process("analysis")

#---------
# Triggers
#---------

if options.physics == "muon":
  print("")
  print("####")
  print("Muon")
  print("####")
  print("")
  if(options.year == "2017"):
    triggerlist = 'HLT_IsoMu27_v*', 'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*', 'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v*', 'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*', 'HLT_DoubleMu43NoFiltersNoVtx_v*', 'HLT_DoubleMu48NoFiltersNoVtx_v*'
  if(options.year == "2018"):
    triggerlist = 'HLT_TkMu100_v*','HLT_Mu50_v*','HLT_IsoMu24_eta2p1_v*','HLT_IsoMu24_v*','HLT_IsoMu27_v*','HLT_IsoMu30_v*','HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_v*','HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1_v*','HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx_v*','HLT_TrkMu16_DoubleTrkMu6NoFiltersNoVtx_v*','HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx_v*','HLT_Mu37_TkMu27_v*','HLT_DoubleL2Mu50_v*','HLT_DoubleMu3_DZ_PFMET90_PFMHT90_v*','HLT_DoubleMu48NoFiltersNoVtx_v*','HLT_DoubleMu40NoFiltersNoVtxDisplaced_v*','HLT_Mu18_Mu9_v*','HLT_TripleMu_12_10_5_v*'
elif options.physics == "jet":
  print("")
  print("####")
  print("BJet")
  print("####")
  print("")
  if(options.year == "2017"):
    triggerlist = 'HLT_BTagMu_AK4DiJet20_Mu5_v*', 'HLT_BTagMu_AK4DiJet40_Mu5_v*', 'HLT_BTagMu_AK4DiJet70_Mu5_v*', 'HLT_BTagMu_AK4DiJet110_Mu5_v*', 'HLT_BTagMu_AK4DiJet170_Mu5_v*', 'HLT_BTagMu_AK4Jet300_Mu5_v*', 'HLT_BTagMu_AK8DiJet170_Mu5_v*', 'HLT_BTagMu_AK8Jet300_Mu5_v*', 'HLT_PFHT380_SixJet32_DoubleBTagCSV_p075_v*', 'HLT_PFHT430_SixJet40_BTagCSV_p080_v*'
  if(options.year == "2018"):  
    triggerlist = 'HLT_BTagMu_AK8Jet300_Mu5_v*',' HLT_BTagMu_AK8Jet300_Mu5_noalgo_v*',' HLT_BTagMu_AK8Jet170_DoubleMu5_v*',' HLT_BTagMu_AK8Jet170_DoubleMu5_noalgo_v*','HLT_BTagMu_AK8DiJet170_Mu5_v*','HLT_BTagMu_AK8DiJet170_Mu5_noalgo_v*','HLT_BTagMu_AK4Jet300_Mu5_v*','HLT_BTagMu_AK4Jet300_Mu5_noalgo_v*','HLT_BTagMu_AK4DiJet70_Mu5_v*','HLT_BTagMu_AK4DiJet70_Mu5_noalgo_v*','HLT_BTagMu_AK4DiJet40_Mu5_v*','HLT_BTagMu_AK4DiJet40_Mu5_noalgo_v*','HLT_BTagMu_AK4DiJet20_Mu5_v*','HLT_BTagMu_AK4DiJet20_Mu5_noalgo_v*','HLT_BTagMu_AK4DiJet170_Mu5_v*','HLT_BTagMu_AK4DiJet170_Mu5_noalgo_v*','HLT_BTagMu_AK4DiJet110_Mu5_v*','HLT_BTagMu_AK4DiJet110_Mu5_noalgo_v*'
elif options.physics == "electron":
  print("")
  print("########")
  print("Electron")
  print("########")
  print("")
  if(options.year == "2017"):
    triggerlist = 'HLT_Ele27_WPTight_Gsf_v*', 'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*', 'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*', 'HLT_DoubleEle33_CaloIdL_MW_v*', 'HLT_BTagMu_AK4DiJet20_Mu5_v*', 'HLT_BTagMu_AK4DiJet40_Mu5_v*', 'HLT_BTagMu_AK4DiJet70_Mu5_v*', 'HLT_BTagMu_AK4DiJet110_Mu5_v*', 'HLT_BTagMu_AK4DiJet170_Mu5_v*', 'HLT_BTagMu_AK4Jet300_Mu5_v*', 'HLT_BTagMu_AK8DiJet170_Mu5_v*', 'HLT_BTagMu_AK8Jet300_Mu5_v*', 'HLT_PFHT380_SixJet32_DoubleBTagCSV_p075_v*', 'HLT_PFHT430_SixJet40_BTagCSV_p080_v*','HLT_PFMET100_PFMHT100_IDTight_CaloBTagCSV_3p1_v*', 'HLT_PFMET110_PFMHT110_IDTight_CaloBTagCSV_3p1_v*', 'HLT_PFMET120_PFMHT120_IDTight_CaloBTagCSV_3p1_v*', 'HLT_PFMET130_PFMHT130_IDTight_CaloBTagCSV_3p1_v*','HLT_PFMET140_PFMHT140_IDTight_CaloBTagCSV_3p1_v*', 'HLT_Ele15_IsoVVVL_PFHT450_CaloBTagCSV_4p5_v*', 'DST_HT250_CaloBTagScouting_v*','HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_v*'
  if(options.year == "2018"):
    triggerlist = 'HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_v*','HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55_v*','HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto_Mass55_v*','HLT_DoubleEle25_CaloIdL_MW_v*','HLT_DoubleEle27_CaloIdL_MW_v*','HLT_DoubleEle33_CaloIdL_MW_v*','HLT_DoublePhoton70_v*','HLT_DoublePhoton85_v*','HLT_Ele135_CaloIdVT_GsfTrkIdT_v*','HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*','HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_v*','HLT_Ele27_Ele37_CaloIdL_MW_v*','HLT_Ele27_WPTight_Gsf_v*','HLT_Photon110EB_TightID_TightIso_v*','HLT_Photon200_v*','HLT_TriplePhoton_20_20_20_CaloIdLV2_R9IdVL_v*','HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350_v*'
elif options.physics == "displacedjet":
  print("")
  print("##############")
  print("Displaced Jets")
  print("##############")
  print("")
  if(options.year == "2017"):
    print("\n Please, there is no displaced jets for 2017!\n")
    sys.exit(0)
  if(options.year == "2018"):
    triggerlist='HLT_HT400_DisplacedDijet40_DisplacedTrack_v*','HLT_HT425_v*','HLT_HT430_DisplacedDijet40_DisplacedTrack_v*','HLT_HT430_DisplacedDijet60_DisplacedTrack_v*','HLT_HT500_DisplacedDijet40_DisplacedTrack_v*','HLT_HT550_DisplacedDijet60_Inclusive_v*','HLT_HT650_DisplacedDijet60_Inclusive_v*'
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
  print("leading leptons not associated with jets cone (unmatching = True)")
  print("#################################################################")
  print("")
else:
  print("##################################################################")
  print("leading leptons associated with the jets cone (unmatching = False)")
  print("##################################################################")
  print("")

if options.discriminator:
  print("#####################################################################################")
  print("selecting jets ordered by the two biggests b-tag discriminator (discriminator = True)")
  print("#####################################################################################")
  print("")
else:
  print("################################################################")
  print("selecting two leading jets ordered by pT (discriminator = False)")
  print("################################################################")
  print("")

if options.debugging:
  print("################")
  print("debugging active")
  print("################")
  print("")

#-----------------
# General options
#-----------------

process.load("FWCore.MessageService.MessageLogger_cfi")
process.options   = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    allowUnscheduled = cms.untracked.bool(True),
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#-------------
# Input files    
#-------------
process.source = cms.Source("PoolSource",
				fileNames = cms.untracked.vstring(fileinput),
)

#----------
# Triggers       
#----------
process.load("PPSFramework.Skimmer.HLTFilter_cfi")
process.hltFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.hltFilter.HLTPaths = (triggerlist)

#-------------
# Preskimming      
#-------------
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag

if mc:
 if options.year=="2017":
  process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mc2017_realistic_v17') # 94X MC campaing
 if options.year=="2018":
  # 106X_upgrade2018_realistic_v20
  process.GlobalTag = GlobalTag(process.GlobalTag, '106X_upgrade2018_realistic_v15_L1v1')
else:
 if options.year=="2017":
  process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v22') # global tag which includes PPS alignment and optics. Default: auto:run2_data
 if options.year=="2018":
  process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v32') # global tag which includes PPS alignment and optics. Default: auto:run2_data

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

#---------------
# Proton Filter     
#---------------
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

#----------
# Analysis        
#----------
process.load("PPSFramework.Skimmer.countsAnalyzer_cfi")
process.CounterNoCuts = process.countsAnalyzer.clone()
process.CounterAfterTrigger = process.countsAnalyzer.clone()
process.CounterTriggered = process.countsAnalyzer.clone()
process.CounterPPSTagged = process.countsAnalyzer.clone()

process.load("PPSFramework.Skimmer.MissingMassEfficiency_cfi")
process.missing_mass_efficiency.debugging = cms.bool(options.debugging)
process.missing_mass_efficiency.verbosity = cms.int32(options.verbosity)
process.missing_mass_efficiency.physics = cms.string(options.physics) 
process.missing_mass_efficiency.unmatching = cms.bool(options.unmatching) #if True, force lepton/jets to be separated
process.missing_mass_efficiency.orderDiscriminator = cms.bool(options.discriminator)
process.missing_mass_efficiency.enableMC = cms.bool(mc)
process.missing_mass_efficiency.enableTrigger = cms.bool(options.trigger)
process.missing_mass_efficiency.enablePrescales = cms.bool(options.prescale)
process.missing_mass_efficiency.triggersList = process.hltFilter.HLTPaths
process.missing_mass_efficiency.GenPartTag = cms.InputTag('prunedGenParticles')
process.missing_mass_efficiency.GenMETTag = cms.InputTag('genMetTrue')
process.missing_mass_efficiency.GenJetAlgoA = cms.InputTag('slimmedGenJets')
process.missing_mass_efficiency.GenJetAlgoB = cms.InputTag('slimmedGenJetsAK8')
process.missing_mass_efficiency.patJetAlgoA = cms.InputTag('slimmedJets')
process.missing_mass_efficiency.patJetAlgoB = cms.InputTag('slimmedJetsAK8')
process.missing_mass_efficiency.electronTag = cms.InputTag("slimmedElectrons")
process.missing_mass_efficiency.muonTag = cms.InputTag("slimmedMuons")
process.missing_mass_efficiency.photonTag = cms.InputTag('slimmedPhotons')
process.missing_mass_efficiency.patmetTag = cms.InputTag('slimmedMETs')
process.missing_mass_efficiency.vertexTag = cms.InputTag('offlineSlimmedPrimaryVertices')
process.missing_mass_efficiency.protonSingleTag = cms.InputTag("ctppsProtons", "singleRP")
process.missing_mass_efficiency.protonMultiTag = cms.InputTag("ctppsProtons", "multiRP")

# b-tagging for RECO
process.load('RecoBTag/Configuration/RecoBTag_cff')
process.load("JetMETCorrections.Type1MET.correctedMet_cff")

# prepare the output file
process.TFileService = cms.Service('TFileService',
    fileName = cms.string("output_efficiency.root"),
    closeFileFast = cms.untracked.bool(True)
)

process.MiniAOD = cms.Sequence(     
                                        process.hltFilter
                                        *process.missing_mass_efficiency
                                )

if mc:
   process.MiniAOD.insert(0, process.beamDivergenceVtxGenerator)
   process.MiniAOD.insert(1, process.ctppsDirectProtonSimulation)
   process.MiniAOD.insert(2, process.reco_local)
   process.MiniAOD.insert(3, process.ctppsProtons)

if not options.trigger:
   process.MiniAOD.remove(process.hltFilter)

process.p = cms.Path(process.MiniAOD)
