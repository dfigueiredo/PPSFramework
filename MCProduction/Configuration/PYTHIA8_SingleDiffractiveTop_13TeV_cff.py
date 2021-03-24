import FWCore.ParameterSet.Config as cms

generator = cms.EDFilter("Pythia8GeneratorFilter",
    maxEventsToPrint = cms.untracked.int32(1),
    pythiaPylistVerbosity = cms.untracked.int32(1),
    filterEfficiency = cms.untracked.double(1.0),
    pythiaHepMCVerbosity = cms.untracked.bool(True),
    comEnergy = cms.double(13000.),

    PythiaParameters = cms.PSet(
        py8DiffSettings = cms.vstring(
                        'HardQCD:all = on',
                        'PhaseSpace:pTHatMin = 80.' # changing top pT
                        'Diffraction:doHard = on',
                        'Diffraction:sampleType = 4',
                        'Diffraction:sampleType = 3',
                        'SigmaDiffractive:PomFlux = 7',
                        'PDF:PomSet = 6',
                        'PDF:PomSet = 4'
        ),
        py8ProcessSettings = cms.vstring(
                                        'Top:gg2ttbar = on',
                                        'Top:qqbar2ttbar = on',
                                        'Top:ffbar2ttbar(s:gmZ) = on',
                                        'Top:gmgm2ttbar = on'
        ),
        parameterSets = cms.vstring( 'py8DiffSettings',
                                     'py8ProcessSettings'
                                   )
    )
)

