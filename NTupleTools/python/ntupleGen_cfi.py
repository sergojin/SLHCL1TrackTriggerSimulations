import FWCore.ParameterSet.Config as cms

#from SLHCL1TrackTriggerSimulations.NTupleTools.prunedGenParticles_cfi import prunedGenParticles

ntupleGenParticles = cms.EDProducer('NTupleGenParticles',
    #inputTag = cms.InputTag('prunedGenParticles'),
    inputTag = cms.InputTag('genParticles'),
    prefix = cms.string('genParts@'),
    suffix = cms.string(''),
    cut = cms.string(''),
    maxN = cms.uint32(999999)
)

from RecoJets.Configuration.GenJetParticles_cff import genParticlesForJetsNoMuNoNu
from RecoJets.Configuration.RecoGenJets_cff import ak5GenJetsNoMuNoNu
ak4GenJetsNoMuNoNu = ak5GenJetsNoMuNoNu.clone( rParam = 0.4 )

ntupleGenJets = cms.EDProducer('NTupleGenJets',
    inputTag = cms.InputTag('ak4GenJetsNoMuNoNu'),
    prefix = cms.string('genJets@'),
    suffix = cms.string(''),
    cut = cms.string(''),
    maxN = cms.uint32(999999)
)

from RecoMET.Configuration.GenMETParticles_cff import genParticlesForMETAllVisible
from RecoMET.Configuration.RecoGenMET_cff import genMetTrue

ntupleGenMET = cms.EDProducer('NTupleGenMET',
    inputTag = cms.InputTag('genMetTrue'),
    prefix = cms.string('genMET@'),
    suffix = cms.string(''),
    cut = cms.string(''),
    maxN = cms.uint32(999999)
)

from SLHCL1TrackTriggerSimulations.NTupleTools.simBeamSpot_cfi import simBeamSpot

ntupleBeamSpot = cms.EDProducer('NTupleBeamSpot',
    inputTag = cms.InputTag('simBeamSpot', 'BeamSpot'),
    prefix = cms.string('beamSpot@'),
    suffix = cms.string(''),
    cut = cms.string(''),
    maxN = cms.uint32(999999)
)

ntupleGenEventInfo = cms.EDProducer('NTupleGenEventInfo',
    genEventInfo = cms.InputTag('generator'),
    pileupInfo = cms.InputTag('addPileupInfo'),
    pileupWeight = cms.InputTag(''),
    pdfWeight = cms.InputTag(''),
    #randomSeed = cms.InputTag('randomEngineSeedKeeper'),
    prefix = cms.string('gen@'),
    suffix = cms.string(''),
)

#ntupleGen = cms.Sequence((prunedGenParticles * ntupleGenParticles) * (genParticlesForJetsNoMuNoNu * ak4GenJetsNoMuNoNu * ntupleGenJets) * (genParticlesForMETAllVisible * genMetTrue * ntupleGenMET) * (simBeamSpot * ntupleBeamSpot) * ntupleGenEventInfo)
ntupleGen = cms.Sequence(ntupleGenParticles * (genParticlesForJetsNoMuNoNu * ak4GenJetsNoMuNoNu * ntupleGenJets) * (genParticlesForMETAllVisible * genMetTrue * ntupleGenMET) * (simBeamSpot * ntupleBeamSpot) * ntupleGenEventInfo)
#ntupleGen = cms.Sequence(ntupleGenParticles * (simBeamSpot * ntupleBeamSpot) * ntupleGenEventInfo)

