#ifndef AMSimulation_PatternMerging_h_
#define AMSimulation_PatternMerging_h_

#include "SLHCL1TrackTriggerSimulations/AMSimulationIO/interface/PatternBankReader.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulationIO/interface/PatternBankWriter.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/Helper.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/ProgramOption.h"
using namespace slhcl1tt;

// Codes originally written by Luciano Ristori (FNAL) and developed by
// Roberto Rossin (Florida, now Padova)
// Modified for inclusion into AMSimulation


class PatternMerging {
public:
  typedef PatternBankReaderT<kPatternMerging> PatternBankReader;
  //typedef PatternBankWriterT<kPatternMerging> PatternBankWriter;

  // ___________________________________________________________________________
  // Sibling

  struct PMSibling {
    //unsigned patternInd;
    //unsigned siblingInd;
    unsigned index;  // siblingInd
    unsigned frequency;
    int layer;
    int delta;
  };

  // ___________________________________________________________________________
  // Pattern

  struct PMPattern {
    // pattern proper (vector of superstrips)
    std::vector<unsigned> superstripIds;

    // popularity
    unsigned frequency;

    // attributes
    //float invPt_mean;
    //float invPt_sigma;
    //float cotTheta_mean;
    //float cotTheta_sigma;
    float phi_mean;
    //float phi_sigma;
    //float z0_mean;
    //float z0_sigma;

    // index in original frequency-sorted list
    unsigned index;

    // Definition of Sibling
    // Requires two patterns to differ in one and only one layer
    // In case of success returns true and assigns layer number and delta
    // where delta is the distance of non-matching superstrip
    // layer and delta are undefined in case of false
    bool isSibling(const PMPattern& p, int& layer, int& delta) const {
      int nDiff = 0;
      for (int i = 0; i < nLayers && nDiff <= 1; ++i) {
        if (superstripIds[i] != p.superstripIds[i]) {
          ++nDiff;
          layer = i;
          delta = (int) p.superstripIds[i] - (int) superstripIds[i];
        }
      }
      return (nDiff == 1); // one and only one layer
    }

  };

  // ___________________________________________________________________________
  // PatternMerging

  // Constructor
  PatternMerging(const ProgramOption& po) :
      po_(po),
      nEvents_(po.maxEvents), verbose_(po.verbose) {}

  // Destructor
  ~PatternMerging() {}

  // Main driver
  int run();

  int mergePatterns(TString src, TString out, unsigned deltaN, float targetCoverage) const;

  void selectSiblings(
      unsigned patternInd,
      const std::vector<PMSibling>& siblings,
      const std::vector<PMPattern>& patternList,
      const std::vector<bool>& merged,
      const std::map<std::vector<unsigned>, unsigned>& patternMap,
      std::vector<unsigned>& selectedSiblings,
      unsigned& selectedFrequency
  ) const;

private:
  static const int nLayers = 6;

  // Program options
  const ProgramOption po_;
  long long nEvents_;
  int verbose_;
};

#endif
