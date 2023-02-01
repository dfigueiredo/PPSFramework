import FWCore.ParameterSet.Config as cms

missing_mass_efficiency = cms.EDAnalyzer('MissingMassEfficiency',
    triggerResults = cms.InputTag('TriggerResults', '', 'HLT')
)
