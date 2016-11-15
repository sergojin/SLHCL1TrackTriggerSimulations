import FWCore.ParameterSet.Config as cms

process = cms.Process('NTUPLE')
runOnMC = True

class Options:
    pass

options = Options()
options.inputFiles = ['root://cmsxrootd.fnal.gov///store/mc/TTI2023Upg14D/PYTHIA6_Tauola_TTbar_TuneZ2star_14TeV/GEN-SIM-DIGI-RAW/PU200_DES23_62_V1-v1/110000/004C20AB-4D9E-E611-AE77-00266CFFBDAC.root']
options.outputFile = 'ntuple.root'
options.maxEvents = -1


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

## Make the ntuple
process.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.outputFile)
)
process.load("SLHCL1TrackTriggerSimulations.NTupleTools.sequences_cff")

## Paths and schedule
process.L1TrackTrigger_step=cms.Path(process.L1TrackTrigger_custom)
process.p = cms.Path(process.ntupleSequence_TTI)
process.schedule = cms.Schedule(process.L1TrackTrigger_step,process.p)


# Configure framework report and summary
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.MessageLogger.cerr.FwkReport.reportEvery = 1

