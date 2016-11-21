#include "SLHCL1TrackTriggerSimulations/NTupleTools/interface/EDMHacker.h"

#include "SLHCL1TrackTriggerSimulations/NTupleTools/interface/ModuleIdFunctor.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTClusterAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTStubAssociationMap.h"

#include <cmath>


EDMHacker::EDMHacker(const edm::ParameterSet& iConfig) :
  inputTag_    (iConfig.getParameter<edm::InputTag>("inputTag")),
  ntupleFile_  (iConfig.getParameter<std::string>("ntupleFile")),
  verbose_     (iConfig.getUntrackedParameter<int>("verbosity", 0)),
  tfile_(0), ttree_(0), ttreeEntry_(0)
{
    // Output product
    produces<std::vector<TTTrack<Ref_PixelDigi_> > >("AMTTTracks");
}

EDMHacker::~EDMHacker() {
    if (ttree_)  delete ttree_;
    if (tfile_)  delete tfile_;
}

void EDMHacker::beginJob() {
    tfile_ = TFile::Open(ntupleFile_.c_str());
    if (verbose_)  std::cout << "[INFO] opened file: " << ntupleFile_ << std::endl;

    ttree_ = (TTree*) tfile_->Get("ntupler/tree");
    assert(ttree_ != 0);

    // For TTree branches
    vb_modId    = 0;
    vt_pt       = 0;
    vt_eta      = 0;
    vt_rinv     = 0;
    vt_invPt    = 0;
    vt_phi0     = 0;
    vt_cottheta = 0;
    vt_z0       = 0;
    vt_d0       = 0;
    vt_chi2     = 0;
    vt_ndof     = 0;
    vt_tower    = 0;
    vt_stubRefs = 0;

    // Stub information
    ttree_->SetBranchAddress("TTStubs_modId"     , &(vb_modId));

    // Track information
    TString prefix = "AMTTTracks_", suffix = "";
    ttree_->SetBranchAddress(prefix + "pt"           + suffix, &(vt_pt));
    ttree_->SetBranchAddress(prefix + "eta"          + suffix, &(vt_eta));
    ttree_->SetBranchAddress(prefix + "rinv"         + suffix, &(vt_rinv));
    ttree_->SetBranchAddress(prefix + "invPt"        + suffix, &(vt_invPt));
    ttree_->SetBranchAddress(prefix + "phi0"         + suffix, &(vt_phi0));
    ttree_->SetBranchAddress(prefix + "cottheta"     + suffix, &(vt_cottheta));
    ttree_->SetBranchAddress(prefix + "z0"           + suffix, &(vt_z0));
    ttree_->SetBranchAddress(prefix + "d0"           + suffix, &(vt_d0));
    ttree_->SetBranchAddress(prefix + "chi2"         + suffix, &(vt_chi2));
    ttree_->SetBranchAddress(prefix + "ndof"         + suffix, &(vt_ndof));
    ttree_->SetBranchAddress(prefix + "tower"        + suffix, &(vt_tower));
    ttree_->SetBranchAddress(prefix + "stubRefs"     + suffix, &(vt_stubRefs));

    if (verbose_)  std::cout << "[INFO] got ttree with num entries: " << ttree_->GetEntries() << std::endl;
}

void EDMHacker::endJob() {
    assert(ttreeEntry_ == ttree_->GetEntries());
}

void EDMHacker::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup) {
    /// Geometry setup
    edm::ESHandle<TrackerGeometry> geometryHandle;
    iSetup.get<TrackerDigiGeometryRecord>().get(geometryHandle);
    theGeometry = geometryHandle.product();

    edm::ESHandle<StackedTrackerGeometry> stackedGeometryHandle;
    iSetup.get<StackedTrackerGeometryRecord>().get(stackedGeometryHandle);
    theStackedGeometry = stackedGeometryHandle.product();
    assert(theStackedGeometry->getCBC3MaxStubs() == 0);

    /// Magnetic field setup
    edm::ESHandle<MagneticField> magneticFieldHandle;
    iSetup.get<IdealMagneticFieldRecord>().get(magneticFieldHandle);
    theMagneticField = magneticFieldHandle.product();
}

void EDMHacker::endRun(const edm::Run& iRun, const edm::EventSetup& iSetup) {

}

void EDMHacker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

    // _________________________________________________________________________
    // Prepare output
    std::auto_ptr<std::vector<TTTrack<Ref_PixelDigi_> > > out_tracks(new std::vector<TTTrack<Ref_PixelDigi_> >());

    // _________________________________________________________________________
    // Load TTree

    ttree_->GetEntry(ttreeEntry_);

    // _________________________________________________________________________
    // Get handles

    typedef edmNew::DetSetVector<TTStub<Ref_PixelDigi_> >    dsv_stub;
    typedef edmNew::DetSet<TTStub<Ref_PixelDigi_> >          ds_stub;
    typedef edm::Ref<dsv_stub, TTStub<Ref_PixelDigi_> >      ref_stub;
    typedef TTStubAssociationMap<Ref_PixelDigi_>             assocmap_stub;

    edm::Handle<dsv_stub> pixelDigiTTStubs;
    iEvent.getByLabel(inputTag_, pixelDigiTTStubs);

    /// Prepare detId -> moduleId
    ModuleIdFunctor getModuleId;

    // _________________________________________________________________________
    // Collect TTStubs

    const unsigned nstubs = vb_modId->size();

    if (verbose_)  std::cout << "Processing event: " << ttreeEntry_ << " expect # stubs: " << nstubs << std::endl;

    std::vector<ref_stub> allStubRefs;

    if (pixelDigiTTStubs.isValid()) {
        edm::LogInfo("EDMHacker") << "Size: " << pixelDigiTTStubs->size();

        unsigned n = 0;

        for (dsv_stub::const_iterator itv = pixelDigiTTStubs->begin(); itv != pixelDigiTTStubs->end(); ++itv) {
            for (ds_stub::const_iterator it = itv->begin(); it != itv->end(); ++it) {
                const StackedTrackerDetId stackDetId(it->getDetId());

                if (verbose_)  std::cout << ".. Processing stub: " << n << std::endl;

                ref_stub iref = edmNew::makeRefTo(pixelDigiTTStubs, it);
                allStubRefs.push_back(iref);

                // Sanity check
                const DetId geoId0 = theStackedGeometry->idToDet(stackDetId, 0)->geographicalId();
                const DetId geoId1 = theStackedGeometry->idToDet(stackDetId, 1)->geographicalId();
                const unsigned moduleId0 = getModuleId(geoId0);
                const unsigned moduleId1 = getModuleId(geoId1);
                assert(moduleId0 == moduleId1);
                edm::LogInfo("EDMHacker") << "geoId0: " << geoId0.rawId() << " geoId1: " << geoId1.rawId() << " modId0: " << moduleId0 << " modId1: " << moduleId1;

                unsigned stub_moduleId = vb_modId->at(n);
                assert(stub_moduleId == moduleId0);

                ++n;
            }
        }

    } else {
        edm::LogError("EDMHacker") << "Cannot get the product: " << inputTag_;
    }

    // Sanity check
    assert(allStubRefs.size() == nstubs);

    // _________________________________________________________________________
    // Build TTTrack

    const unsigned ntracks = vt_pt->size();

    for (unsigned itrack = 0; itrack < ntracks; ++itrack) {

        if (verbose_)  std::cout << ".. Processing track: " << itrack << std::endl;

        const auto& track_pt       = vt_pt      ->at(itrack);
        const auto& track_eta      = vt_eta     ->at(itrack);
        const auto& track_rinv     = vt_rinv    ->at(itrack);
        //const auto& track_invPt    = vt_invPt   ->at(itrack);
        const auto& track_phi0     = vt_phi0    ->at(itrack);
        //const auto& track_cottheta = vt_cottheta->at(itrack);
        const auto& track_z0       = vt_z0      ->at(itrack);
        //const auto& track_d0       = vt_d0      ->at(itrack);
        const auto& track_chi2     = vt_chi2    ->at(itrack);
        //const auto& track_ndof     = vt_ndof    ->at(itrack);
        const auto& track_tower    = vt_tower   ->at(itrack);
        const auto& track_stubRefs = vt_stubRefs->at(itrack);

        // Track stubs
        std::vector<ref_stub> theStubRefs;

        unsigned track_nstubs = track_stubRefs.size();

        for (unsigned istub = 0; istub < track_nstubs; ++istub) {
            unsigned myStubRef = track_stubRefs.at(istub);
            if (myStubRef == 999999999)  continue;
            assert(myStubRef < nstubs);
            ref_stub iref = allStubRefs.at(myStubRef);
            theStubRefs.push_back(iref);
        }

        // Track parameters
        GlobalVector mom(GlobalVector::Cylindrical(track_pt, track_phi0, track_pt * std::sinh(track_eta)));
        GlobalPoint poca(0.0, 0.0, track_z0);
        int nPar = 4;

        // Set track properties
        TTTrack<Ref_PixelDigi_> out_track(theStubRefs);
        out_track.setMomentum(mom, nPar);
        out_track.setRInv(track_rinv, nPar);
        out_track.setPOCA(poca, nPar);
        out_track.setSector(track_tower);
        out_track.setWedge(999999);
        out_track.setChi2(track_chi2, nPar);
        out_track.setStubPtConsistency(999999., nPar);

        // Done!
        out_tracks->push_back(out_track);
    }

    assert(out_tracks->size() == ntracks);

    // _________________________________________________________________________
    // Increment ttree entry number
    ++ttreeEntry_;

    // _________________________________________________________________________
    // Save output
    iEvent.put(out_tracks, "AMTTTracks");
}
