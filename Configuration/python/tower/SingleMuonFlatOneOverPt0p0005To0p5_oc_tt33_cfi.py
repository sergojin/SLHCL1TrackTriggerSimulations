import FWCore.ParameterSet.Config as cms

# reference: https://github.com/cms-sw/genproductions/blob/master/python/EightTeV/SingleMuMinusFlatPt0p2To100_cff.py

# Calculate tower eta, phi boundaries
_tower = 33
_pi = 3.141592653589793
_min_eta = (_tower/8)*2.2/3 - 2.2
_max_eta = _min_eta + 2.2/3
_min_phi = (((_tower%8)+4)%8)*_pi/4 - _pi
_max_phi = _min_phi + _pi/4
#print _min_eta, _max_eta, _min_phi, _max_phi

generator = cms.EDProducer("FlatRandomOneOverPtGunProducer2",
    PGunParameters = cms.PSet(
        MaxOneOverPt = cms.double(0.5),
        MinOneOverPt = cms.double(0.0005),
        PartID = cms.vint32(13),
        MinEta = cms.double(_min_eta),
        MaxEta = cms.double(_max_eta),
        MinPhi = cms.double(_min_phi),
        MaxPhi = cms.double(_max_phi),
        #XFlatSpread = cms.double(1.5),  ## in mm
        #YFlatSpread = cms.double(1.5),  ## in mm
        #ZFlatSpread = cms.double(150.), ## in mm
        RandomCharge = cms.bool(True),
        ReallyFlat = cms.bool(True),
        VtxSmeared = cms.bool(True),
        UseRStar = cms.bool(True),
    ),
    Verbosity = cms.untracked.int32(0),
    psethack = cms.string('single muon+/- 1/pt 0.0005 to 0.5 tower %i' % _tower),
    AddAntiParticle = cms.bool(False),
    firstRun = cms.untracked.uint32(1)
)
