#ifndef AMSimulationDataFormats_PairAssignment_h_
#define AMSimulationDataFormats_PairAssignment_h_

namespace slhcl1tt{

struct PairAssignment //structure for saving pairs for the pairwise Delta Delta S combination cleanup
{
  unsigned Stub1Address;
  unsigned Stub2Address;
  bool SurvivingPair;
  PairAssignment(unsigned Stub1Address_, unsigned Stub2Address_, bool SurvivingPair_){
    Stub1Address=Stub1Address_;
    Stub2Address=Stub2Address_;
    SurvivingPair=SurvivingPair_;
  }
};

}

#endif
