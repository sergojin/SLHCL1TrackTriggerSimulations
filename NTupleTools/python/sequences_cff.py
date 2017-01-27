import FWCore.ParameterSet.Config as cms

from SLHCL1TrackTriggerSimulations.NTupleTools.ntupleGen_cfi import *
from SLHCL1TrackTriggerSimulations.NTupleTools.ntupleGenExtra_cfi import *
from SLHCL1TrackTriggerSimulations.NTupleTools.ntupleSim_cfi import *
from SLHCL1TrackTriggerSimulations.NTupleTools.ntupleDigi_cfi import *
from SLHCL1TrackTriggerSimulations.NTupleTools.ntupleL1TrackTrigger_cfi import *
from SLHCL1TrackTriggerSimulations.NTupleTools.ntupleMaker_cfi import *

ntupleSequence = cms.Sequence(ntupleGen * ntupleGenExtra * ntupleSim * ntupleDigi * ntupleL1TrackTrigger * (ntupleEventInfo * ntupler))

ntupleSequence_TTI = cms.Sequence(ntupleGen * ntupleSim * ntupleDigi_TTI * ntupleL1TrackTrigger_TTI * (ntupleEventInfo * ntupler))

ntupleSequence_GENSIM = cms.Sequence(ntupleGen * ntupleSim * (ntupleEventInfo * ntupler))

