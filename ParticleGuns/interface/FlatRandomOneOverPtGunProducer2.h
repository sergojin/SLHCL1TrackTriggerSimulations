#ifndef FlatRandomOneOverPtGunProducer2_H
#define FlatRandomOneOverPtGunProducer2_H

/** \class FlatRandomOneOverPtGunProducer2
 *
 * Generates single particle gun flat in (1/pt) in HepMC format
 **************************************************************/

#include "SLHCL1TrackTriggerSimulations/ParticleGuns/interface/BaseFlatGunProducer2.h"

namespace edm
{
  
  class FlatRandomOneOverPtGunProducer2 : public BaseFlatGunProducer2
  {
  
  public:
    FlatRandomOneOverPtGunProducer2(const ParameterSet & pset);
    virtual ~FlatRandomOneOverPtGunProducer2();
   
    virtual void produce(Event & e, const EventSetup& es) override;

  private:
    
    // data members
    
    double            fMinOneOverPt   ;
    double            fMaxOneOverPt   ;
    double            fXFlatSpread   ;
    double            fYFlatSpread   ;
    double            fZFlatSpread   ;
    bool              fRandomCharge  ;
    bool              fReallyFlat    ;
    bool              fVtxSmeared    ;
    bool              fUseRStar      ;

  };
} 

#endif
