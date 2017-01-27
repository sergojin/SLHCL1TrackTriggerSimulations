#ifndef AMSimulation_StubCleaner_h_
#define AMSimulation_StubCleaner_h_

#include "SLHCL1TrackTriggerSimulations/AMSimulationIO/interface/TTStubReader.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulationIO/interface/TTStubWriter.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/Helper.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/HelperMath.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/ProgramOption.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/Picky.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/ModuleOverlapMap.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/TrackParametersToTT.h"
using namespace slhcl1tt;


class StubCleaner {
  public:
    typedef TTStubReaderT<kStubCleaner> TTStubReader;
    typedef TTStubWriterT<kStubCleaner> TTStubWriter;

    // Constructor
    StubCleaner(const ProgramOption& po)
    : po_(po),
      nEvents_(po.maxEvents), verbose_(po.verbose) {

        // Initialize
        picky_ = new Picky();

        momap_ = 0;
        if (po_.removeOverlap) {
            momap_   = new ModuleOverlapMap();
            momap_->readModuleOverlapMap(po_.datadir);
        }
    }

    // Destructor
    ~StubCleaner() {
        if (picky_)  delete picky_;
        if (momap_)  delete momap_;
    }

    // Main driver
    int run();


  private:
    // Member functions
    // Select one unique stub per layer
    int cleanStubs(TString src, TString out);

    // Program options
    const ProgramOption po_;
    long long nEvents_;
    int verbose_;

    // Operators
    Picky * picky_;
    ModuleOverlapMap  * momap_;
};

#endif
