import FWCore.ParameterSet.Config as cms

missing_mass = cms.EDAnalyzer('MissingMassSearches',
    triggerResults = cms.InputTag('TriggerResults', '', 'HLT'),
    energyThresholdHB = cms.double(1.5),
    energyThresholdHE = cms.double(2.0),
    energyThresholdHF = cms.double(4.0),
    energyThresholdEB = cms.double(0.6),
    energyThresholdEE = cms.double(1.5)
)
