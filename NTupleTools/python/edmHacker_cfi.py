import FWCore.ParameterSet.Config as cms

edmHacker = cms.EDProducer("EDMHacker",
    inputTag = cms.InputTag("TTStubsFromPixelDigis", "StubAccepted"),
    ntupleFile = cms.string("tracks.root"),
    verbosity = cms.untracked.int32(0),
)

