import FWCore.ParameterSet.Config as cms

simpleSkimmer = cms.EDFilter("SimpleSkimmer",
    inputTag = cms.InputTag('mix', 'MergedTrackTruth'),
    cut = cms.string(
        "(abs(pdgId) = 13 & status = 1 & pt > 3 & abs(eta) < 2.4) && " # muon
        "(eventId().event() == 0) && "                                 # signal
        "(eventId().bunchCrossing() == 0) && "                         # intime
        "(sqrt(vx*vx+vy*vy) < 0.5 && abs(vz) < 30.0) "                 # primary
    ),
    minN = cms.uint32(1),
)

