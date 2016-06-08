#ifndef AMSimulation_PairCombinationFactory_h_
#define AMSimulation_PairCombinationFactory_h_

#include <vector>
#include <iostream>
#include <cmath>
#include "SLHCL1TrackTriggerSimulations/AMSimulationDataFormats/interface/PairAssignment.h"

namespace slhcl1tt{

class PairCombinationFactory{
  public:
    PairCombinationFactory() : verbose_(true) {}
    
    ~PairCombinationFactory() {}
    
    // Enum
    enum Flag { BAD=999999999 };
    
    std::vector<std::vector<unsigned> > combine(const std::vector<std::vector<unsigned> >& groups, std::vector<std::vector<float> > DeltaS, bool FiveOfSix);
    
    //debug
    void print();
    
  private:
    int verbose_;
    
    //DDS value coarsening based on layer number
    float Coarsener(float DeltaS, int precision, int layer);

    //central functions
    std::vector<PairAssignment> Stage1Pair(std::vector<std::pair<unsigned, float> > first, std::vector<std::pair<unsigned, float> > second, float cut);
    std::vector<std::vector<unsigned> > CombinationBuilder(std::vector<PairAssignment> Layer10, std::vector<PairAssignment> Layer12, std::vector<PairAssignment> Layer32, std::vector<PairAssignment> Layer34, std::vector<PairAssignment> Layer54, unsigned l0, unsigned l2, unsigned l4);
};

    
}

#endif
