#ifndef AMSimulation_RoadMerging_h_
#define AMSimulation_RoadMerging_h_

#include "SLHCL1TrackTriggerSimulations/AMSimulationIO/interface/PatternBankReader.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulationIO/interface/TTTrackReader.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulationIO/interface/TTTrackWriter.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/Helper.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/ProgramOption.h"
using namespace slhcl1tt;

// Codes originally written by Roberto Rossin (Florida, now Padova)
// Modified for inclusion into AMSimulation


class RoadMerging {
public:
  typedef PatternBankReaderT<kRoadMerging> PatternBankReader;
  typedef TTTrackReaderT<kRoadMerging> TTRoadReader;
  typedef TTTrackWriterT<kRoadMerging> TTRoadWriter;

  // ___________________________________________________________________________
  // TTRoad

  struct RMTTRoad {
    unsigned patternRef;
    unsigned tower;
    unsigned nstubs;
    float    patternInvPt;
    unsigned patternFreq;

    std::vector<unsigned> superstripIds;
    std::vector<std::vector<unsigned> > stubRefs;  // [layer i][stub j]

    std::vector<std::vector<unsigned> > superstripIdsUnited;  // [layer i][superstrip j]
  };

  // ___________________________________________________________________________
  // Pattern

  struct RMPattern {
    // pattern proper (vector of superstrips)
    std::vector<unsigned> superstripIds;

    // popularity
    unsigned frequency;

    // attributes
    float invPt_mean;
    //float invPt_sigma;
    //float cotTheta_mean;
    //float cotTheta_sigma;
    //float phi_mean;
    //float phi_sigma;
    //float z0_mean;
    //float z0_sigma;

    // index in original frequency-sorted list
    unsigned index;

    // index in merged pattern bank
    unsigned indToMerged;

    // pattern proper big-league (vector of vector of superstrips)
    std::vector<std::vector<unsigned> > superstripIdsUnited;

    // sibling indices in original frequency-sorted list
    std::vector<unsigned> indFromMerged;
  };

  // ___________________________________________________________________________
  // RoadMerging

  // Constructor
  RoadMerging(const ProgramOption& po) :
      po_(po),
      nEvents_(po.maxEvents), verbose_(po.verbose) {}

  // Destructor
  ~RoadMerging() {}

  // Main driver
  int run();

  int processEvents(TString bank, TString src, TString out, float targetCoverage) const;

  void mergeRoads(
      const std::vector<RMPattern>& patterns,
      const std::vector<RMPattern>& merged_patterns,
      const std::vector<RMTTRoad>& roads,
      std::vector<RMTTRoad>& merged_roads
  ) const;

private:
  static const int nLayers = 6;

  // Program options
  const ProgramOption po_;
  long long nEvents_;
  int verbose_;
};

#endif
