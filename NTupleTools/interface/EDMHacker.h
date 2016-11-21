#ifndef NTupleTools_EDMHacker_h_
#define NTupleTools_EDMHacker_h_

#include "SLHCL1TrackTriggerSimulations/NTupleTools/interface/NTupleCommon.h"
#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"

#include "Geometry/Records/interface/StackedTrackerGeometryRecord.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/StackedTrackerGeometry.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "TFile.h"
#include "TTree.h"


class EDMHacker : public edm::EDProducer {
  public:
    explicit EDMHacker(const edm::ParameterSet&);
    virtual ~EDMHacker();

  private:
    virtual void beginJob();
    virtual void produce(edm::Event&, const edm::EventSetup&);
    virtual void endJob();

    virtual void beginRun(const edm::Run&, const edm::EventSetup&);
    virtual void endRun(const edm::Run&, const edm::EventSetup&);

    // For TTree branches
    std::vector<unsigned> *               vb_modId;       // stub module ID
    std::vector<float> *                  vt_pt;          // track momentum [GeV]
    std::vector<float> *                  vt_eta;         // track momentum [GeV]
    std::vector<float> *                  vt_rinv;        // = 0.003 * 3.8 * track fitted q/pT
    std::vector<float> *                  vt_invPt;       // track q/pT [1/GeV]
    std::vector<float> *                  vt_phi0;        // track phi at r=0
    std::vector<float> *                  vt_cottheta;    // track cot theta at r=0
    std::vector<float> *                  vt_z0;          // track vertex z [cm]
    std::vector<float> *                  vt_d0;          // track impact parameter [cm]
    std::vector<float> *                  vt_chi2;        // track fit chi-squared
    std::vector<int> *                    vt_ndof;        // track fit num of degrees of freedom
    std::vector<unsigned> *               vt_tower;       // track tower
    std::vector<std::vector<unsigned> > * vt_stubRefs;    // track stub ids

    // For event setup
    const TrackerGeometry * theGeometry;
    const StackedTrackerGeometry * theStackedGeometry;
    const MagneticField* theMagneticField;

    // Member data
    const edm::InputTag inputTag_;  // Input TTStub
    const std::string ntupleFile_;
    int verbose_;

    TFile* tfile_;
    TTree* ttree_;
    long long ttreeEntry_;
};

#endif
