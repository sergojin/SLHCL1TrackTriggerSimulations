import FWCore.ParameterSet.Config as cms

process = cms.Process('NTUPLE')
runOnMC = True

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')
defaultInputFiles = [
'/store/mc/TTI2023Upg14D/PYTHIA6_Tauola_TTbar_TuneZ2star_14TeV/GEN-SIM-DIGI-RAW/PU200_DES23_62_V1-v1/110000/004C20AB-4D9E-E611-AE77-00266CFFBDAC.root',
#'/store/mc/TTI2023Upg14D/PYTHIA6_Tauola_TTbar_TuneZ2star_14TeV/GEN-SIM-DIGI-RAW/PU200_DES23_62_V1-v1/110000/007F6C2E-5A9E-E611-B0E2-C4346BC8D390.root',
#'/store/mc/TTI2023Upg14D/PYTHIA6_Tauola_TTbar_TuneZ2star_14TeV/GEN-SIM-DIGI-RAW/PU200_DES23_62_V1-v1/110000/0208AEA5-4D9E-E611-ACE8-00266CFFCCC8.root',
#'/store/mc/TTI2023Upg14D/PYTHIA6_Tauola_TTbar_TuneZ2star_14TeV/GEN-SIM-DIGI-RAW/PU200_DES23_62_V1-v1/110000/02A8431A-4E9E-E611-B3D0-00266CFFBF34.root',
]
options.setDefault('inputFiles', defaultInputFiles)
options.setDefault('outputFile', 'edmhack.root')
options.parseArguments()


## MessageLogger
process.load('FWCore.MessageService.MessageLogger_cfi')

## Options
process.options = cms.untracked.PSet(

)

## Input Source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles)
)

## Maximal Number of Events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

## Geometry and Global Tags
process.load('Configuration.Geometry.GeometryExtended2023TTIReco_cff')
process.load('Configuration.Geometry.GeometryExtended2023TTI_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')

process.load('Geometry.TrackerGeometryBuilder.StackedTrackerGeometry_cfi')
process.load('Configuration.StandardSequences.L1TrackTrigger_cff')
process.L1TrackTrigger_custom=cms.Sequence(process.BeamSpotFromSim+process.TTTracksFromPixelDigis+process.TTTrackAssociatorFromPixelDigis)

# Output definition

process.load('Configuration.EventContent.EventContent_cff')

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    fileName = cms.untracked.string(options.outputFile),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW')
    ),
)


## Attach the ntuple back to EDM
process.load("SLHCL1TrackTriggerSimulations.NTupleTools.edmHacker_cfi")

## Paths and schedule
process.p = cms.Path(process.edmHacker)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)
process.schedule = cms.Schedule(process.p,process.RAWSIMoutput_step)

## filter all path with the skim
#process.load("SLHCL1TrackTriggerSimulations.NTupleTools.simpleSkimmer_cfi")
#for path in process.paths:
#    getattr(process,path)._seq = process.simpleSkimmer * getattr(process,path)._seq

# customisation of the process.

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.combinedCustoms
from SLHCUpgradeSimulations.Configuration.combinedCustoms import cust_2023TTI

#call to customisation function cust_2023TTI imported from SLHCUpgradeSimulations.Configuration.combinedCustoms
process = cust_2023TTI(process)

# End of customisation functions


# Configure framework report and summary
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.MessageLogger.cerr.FwkReport.reportEvery = 1

