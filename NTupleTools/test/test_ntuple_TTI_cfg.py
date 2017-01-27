import FWCore.ParameterSet.Config as cms

process = cms.Process('NTUPLE')
runOnMC = True

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')
options.setDefault('inputFiles', ['file:rawsim_numEvent100.root'])
options.setDefault('outputFile', 'ntuple.root')
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

## Make the ntuple
process.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.outputFile)
)
process.load("SLHCL1TrackTriggerSimulations.NTupleTools.sequences_cff")

## Paths and schedule
process.p = cms.Path(process.ntupleSequence_TTI)
process.schedule = cms.Schedule(process.p)


# Configure framework report and summary
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

