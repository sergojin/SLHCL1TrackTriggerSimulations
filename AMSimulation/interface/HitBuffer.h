#ifndef AMSimulation_HitBuffer_h_
#define AMSimulation_HitBuffer_h_

#include "SLHCL1TrackTriggerSimulations/AMSimulationDataFormats/interface/Pattern.h"
#include <map>
#include <vector>

namespace slhcl1tt {

class HitBuffer {
  public:
    // Constructor
    HitBuffer() : nlayers_(0), nss_(0), frozen_(false) {}

    // Destructor
    ~HitBuffer() {}

    // Functions
    int init(unsigned nlayers, unsigned nss);

    void reset();

    void insert(unsigned layer, superstrip_type ss, unsigned stubRef);

    void freeze(unsigned maxStubs);

    bool isHit(unsigned layer, superstrip_type ss) const;

    std::vector<unsigned> getHits(unsigned layer, superstrip_type ss) const;

    // Debug
    void print();

  private:
    // Member data
    std::map<superstrip_type, std::vector<unsigned> > superstripHits_;   // superstrip --> stubRefs (std::map)
    std::vector<bool>                                 superstripBools_;  // superstrip --> hit or empty (hash table)
    unsigned nlayers_;
    unsigned nss_;
    bool frozen_;
};

}

#endif
