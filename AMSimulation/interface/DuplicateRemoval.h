#ifndef AMSimulation_DuplicateRemoval_h_
#define AMSimulation_DuplicateRemoval_h_

#include "SLHCL1TrackTriggerSimulations/AMSimulationDataFormats/interface/TTTrack2.h"
#include <vector>
#include <algorithm>


namespace slhcl1tt {

class DuplicateRemoval {
  public:
    // Constructor
    DuplicateRemoval() {}

    // Destructor
    ~DuplicateRemoval() {}

    // dupRm is a parameter 0..6, if num of shared stubs > dupRm, it is a duplicate
    void checkTracks(std::vector<TTTrack2>& all_tracks, int dupRm);
};

}

#endif
