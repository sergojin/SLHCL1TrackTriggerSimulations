#include "FWCore/Framework/interface/MakerMacros.h"
// #include "IOMC/Input/interface/MCFileSource.h"

// Julia Yarba : related to particle gun prototypes
//

#include "SLHCL1TrackTriggerSimulations/ParticleGuns/interface/FlatRandomPtGunProducer2.h"
#include "SLHCL1TrackTriggerSimulations/ParticleGuns/interface/FlatRandomOneOverPtGunProducer2.h"

// particle gun prototypes
//
  
using edm::FlatRandomPtGunProducer2;
DEFINE_FWK_MODULE(FlatRandomPtGunProducer2);
using edm::FlatRandomOneOverPtGunProducer2;
DEFINE_FWK_MODULE(FlatRandomOneOverPtGunProducer2);
