#ifndef __TTTrackWriter_h__
#define __TTTrackWriter_h__

#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include <memory>
#include <vector>

#include "BasicReader.h"

// _____________________________________________________________________________
// This is a simple wrapper around TTree. It clones a tree and adds new
// branches if needed

namespace slhcl1tt {
class TTRoad;
class TTTrack2;
}

template<size_t T>
class TTTrackWriterT {
public:
    TTTrackWriterT(int verbose=1);
    ~TTTrackWriterT();

    void init(TTree* intree, TString out);

    void fill() { ttree_->Fill(); }

    void fillRoads(const std::vector<slhcl1tt::TTRoad>& roads);

    void fillRoads(const std::vector<slhcl1tt::TTRoad>& roads, const std::vector<std::string>& stubs_bitString, const std::vector<unsigned>& stubs_superstripId);

    void fillTracks(const std::vector<slhcl1tt::TTTrack2>& tracks);

    void write() { tfile_->Write(); }

    TTree * getTree() { return ttree_; }

    // Stub information
    std::unique_ptr<std::vector<std::string> >    vb_bitString;     // stub info encoded as in firmware
    std::unique_ptr<std::vector<unsigned> >       vb_superstripId;  // stub superstrip ID

    // Road information
    std::unique_ptr<std::vector<unsigned> >                             vr_patternRef;    // road pattern id
    std::unique_ptr<std::vector<unsigned> >                             vr_tower;         // road tower
    std::unique_ptr<std::vector<unsigned> >                             vr_nstubs;        // road total number of stubs
    std::unique_ptr<std::vector<float> >                                vr_patternInvPt;  // road pattern q/pT [1/GeV]
    std::unique_ptr<std::vector<unsigned> >                             vr_patternFreq;   // road pattern frequency
    std::unique_ptr<std::vector<std::vector<unsigned> > >               vr_superstripIds; // road superstrip ids
    std::unique_ptr<std::vector<std::vector<std::vector<unsigned> > > > vr_stubRefs;      // road stub ids

    // Track information
    //std::unique_ptr<std::vector<float> >                  vt_px;          // track momentum [GeV]
    //std::unique_ptr<std::vector<float> >                  vt_py;          // track momentum [GeV]
    //std::unique_ptr<std::vector<float> >                  vt_pz;          // track momentum [GeV]
    std::unique_ptr<std::vector<float> >                  vt_pt;          // track momentum [GeV]
    std::unique_ptr<std::vector<float> >                  vt_eta;         // track momentum [GeV]
    //std::unique_ptr<std::vector<float> >                  vt_phi;         // track momentum [GeV]
    //std::unique_ptr<std::vector<float> >                  vt_vx;          // track vertex [cm]
    //std::unique_ptr<std::vector<float> >                  vt_vy;          // track vertex [cm]
    //std::unique_ptr<std::vector<float> >                  vt_vz;          // track vertex [cm]
    std::unique_ptr<std::vector<float> >                  vt_rinv;        // = 0.003 * 3.8 * track fitted q/pT
    std::unique_ptr<std::vector<float> >                  vt_invPt;       // track q/pT [1/GeV]
    std::unique_ptr<std::vector<float> >                  vt_phi0;        // track phi at r=0
    std::unique_ptr<std::vector<float> >                  vt_cottheta;    // track cot theta at r=0
    std::unique_ptr<std::vector<float> >                  vt_z0;          // track vertex z [cm]
    std::unique_ptr<std::vector<float> >                  vt_d0;          // track impact parameter [cm]
    std::unique_ptr<std::vector<float> >                  vt_chi2;        // track fit chi-squared
    std::unique_ptr<std::vector<int> >                    vt_ndof;        // track fit num of degrees of freedom
    std::unique_ptr<std::vector<float> >                  vt_synMatchChi2;// track MC match chi-squared
    std::unique_ptr<std::vector<int> >                    vt_synMatchCat; // track MC match category
    std::unique_ptr<std::vector<int> >                    vt_tpId;        // associated tracking particle id
    std::unique_ptr<std::vector<int> >                    vt_synTpId;     // associated tracking particle id using 'synthetic' approach
    std::unique_ptr<std::vector<unsigned> >               vt_tower;       // track tower
    std::unique_ptr<std::vector<unsigned> >               vt_hitBits;     // track stub composition by layers
    std::unique_ptr<std::vector<unsigned> >               vt_ptSegment;   // obsolete
    std::unique_ptr<std::vector<unsigned> >               vt_roadRef;     // track road id
    std::unique_ptr<std::vector<unsigned> >               vt_combRef;     // track combination id
    std::unique_ptr<std::vector<unsigned> >               vt_patternRef;  // track pattern id
    std::unique_ptr<std::vector<std::vector<unsigned> > > vt_stubRefs;    // track stub ids
    std::unique_ptr<std::vector<std::vector<float> > >    vt_principals;  // track principle components

protected:
    TFile* tfile_;
    TTree* ttree_;
    const int verbose_;
};

// _____________________________________________________________________________
// Implementation is included in the header file to simplify ROOT library generation

#include <cassert>
#include <stdexcept>

#include "SLHCL1TrackTriggerSimulations/AMSimulationDataFormats/interface/TTRoad.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulationDataFormats/interface/TTTrack2.h"
using namespace slhcl1tt;


template<size_t T>
TTTrackWriterT<T>::TTTrackWriterT(int verbose) :
    vb_bitString      (new std::vector<std::string>()),
    vb_superstripId   (new std::vector<unsigned>()),
    //
    vr_patternRef     (new std::vector<unsigned>()),
    vr_tower          (new std::vector<unsigned>()),
    vr_nstubs         (new std::vector<unsigned>()),
    vr_patternInvPt   (new std::vector<float>()),
    vr_patternFreq    (new std::vector<unsigned>()),
    vr_superstripIds  (new std::vector<std::vector<unsigned> >()),
    vr_stubRefs       (new std::vector<std::vector<std::vector<unsigned> > >()),
    //
    //vt_px           (new std::vector<float>()),
    //vt_py           (new std::vector<float>()),
    //vt_pz           (new std::vector<float>()),
    vt_pt           (new std::vector<float>()),
    vt_eta          (new std::vector<float>()),
    //vt_phi          (new std::vector<float>()),
    //vt_vx           (new std::vector<float>()),
    //vt_vy           (new std::vector<float>()),
    //vt_vz           (new std::vector<float>()),
    vt_rinv         (new std::vector<float>()),
    vt_invPt        (new std::vector<float>()),
    vt_phi0         (new std::vector<float>()),
    vt_cottheta     (new std::vector<float>()),
    vt_z0           (new std::vector<float>()),
    vt_d0           (new std::vector<float>()),
    vt_chi2         (new std::vector<float>()),
    vt_ndof         (new std::vector<int>()),
    vt_synMatchChi2 (new std::vector<float>()),
    vt_synMatchCat  (new std::vector<int>()),
    vt_tpId         (new std::vector<int>()),
    vt_synTpId      (new std::vector<int>()),
    vt_tower        (new std::vector<unsigned>()),
    vt_hitBits      (new std::vector<unsigned>()),
    vt_ptSegment    (new std::vector<unsigned>()),
    vt_roadRef      (new std::vector<unsigned>()),
    vt_combRef      (new std::vector<unsigned>()),
    vt_patternRef   (new std::vector<unsigned>()),
    vt_stubRefs     (new std::vector<std::vector<unsigned> >()),
    vt_principals   (new std::vector<std::vector<float> >()),
    //
    tfile_(0), ttree_(0), verbose_(verbose) {}

template<size_t T>
TTTrackWriterT<T>::~TTTrackWriterT() {
    if (ttree_)  delete ttree_;
    if (tfile_)  delete tfile_;
}

template<size_t T>
void TTTrackWriterT<T>::init(TTree* intree, TString out) {
    if (!out.EndsWith(".root")) {
        TString msg = "Output filename must be .root";
        throw std::invalid_argument(msg.Data());
    }

    //if (verbose_)  std::cout << "Opening " << out << std::endl;
    tfile_ = TFile::Open(out, "RECREATE");

    if (tfile_) {
        //if (verbose_)  std::cout << "Successfully opened " << out << std::endl;
    } else {
        TString msg = "Failed to open " + out;
        throw std::invalid_argument(msg.Data());
    }

    tfile_->mkdir("ntupler")->cd();
    ttree_ = (TTree*) intree->CloneTree(0); // Do not copy the data yet

    TString prefix = "", suffix = "";

    if (T == kPatternMatcher) {
        prefix = "AMTTRoads_";
        ttree_->Branch("TTStubs_bitString"   , &(*vb_bitString));
        ttree_->Branch("TTStubs_superstripId", &(*vb_superstripId));
        //
        ttree_->Branch(prefix + "patternRef"    + suffix, &(*vr_patternRef));
        ttree_->Branch(prefix + "tower"         + suffix, &(*vr_tower));
        ttree_->Branch(prefix + "nstubs"        + suffix, &(*vr_nstubs));
        ttree_->Branch(prefix + "patternInvPt"  + suffix, &(*vr_patternInvPt));
        ttree_->Branch(prefix + "patternFreq"   + suffix, &(*vr_patternFreq));
        ttree_->Branch(prefix + "superstripIds" + suffix, &(*vr_superstripIds));
        ttree_->Branch(prefix + "stubRefs"      + suffix, &(*vr_stubRefs));
    }

    if (T == kTrackFitter) {
        prefix = "AMTTTracks_";
        //ttree_->Branch(prefix + "px"             + suffix, &(*vt_px));
        //ttree_->Branch(prefix + "py"             + suffix, &(*vt_py));
        //ttree_->Branch(prefix + "pz"             + suffix, &(*vt_pz));
        ttree_->Branch(prefix + "pt"             + suffix, &(*vt_pt));
        ttree_->Branch(prefix + "eta"            + suffix, &(*vt_eta));
        //ttree_->Branch(prefix + "phi"            + suffix, &(*vt_phi));
        //ttree_->Branch(prefix + "vx"             + suffix, &(*vt_vx));
        //ttree_->Branch(prefix + "vy"             + suffix, &(*vt_vy));
        //ttree_->Branch(prefix + "vz"             + suffix, &(*vt_vz));
        ttree_->Branch(prefix + "rinv"           + suffix, &(*vt_rinv));
        ttree_->Branch(prefix + "invPt"          + suffix, &(*vt_invPt));
        ttree_->Branch(prefix + "phi0"           + suffix, &(*vt_phi0));
        ttree_->Branch(prefix + "cottheta"       + suffix, &(*vt_cottheta));
        ttree_->Branch(prefix + "z0"             + suffix, &(*vt_z0));
        ttree_->Branch(prefix + "d0"             + suffix, &(*vt_d0));
        ttree_->Branch(prefix + "chi2"           + suffix, &(*vt_chi2));
        ttree_->Branch(prefix + "ndof"           + suffix, &(*vt_ndof));
        ttree_->Branch(prefix + "synMatchChi2"   + suffix, &(*vt_synMatchChi2));
        ttree_->Branch(prefix + "synMatchCat"    + suffix, &(*vt_synMatchCat));
        ttree_->Branch(prefix + "tpId"           + suffix, &(*vt_tpId));
        ttree_->Branch(prefix + "synTpId"        + suffix, &(*vt_synTpId));
        ttree_->Branch(prefix + "tower"          + suffix, &(*vt_tower));
        ttree_->Branch(prefix + "hitBits"        + suffix, &(*vt_hitBits));
        ttree_->Branch(prefix + "ptSegment"      + suffix, &(*vt_ptSegment));
        ttree_->Branch(prefix + "roadRef"        + suffix, &(*vt_roadRef));
        ttree_->Branch(prefix + "combRef"        + suffix, &(*vt_combRef));
        ttree_->Branch(prefix + "patternRef"     + suffix, &(*vt_patternRef));
        ttree_->Branch(prefix + "stubRefs"       + suffix, &(*vt_stubRefs));
        ttree_->Branch(prefix + "principals"     + suffix, &(*vt_principals));
    }
}

template<size_t T>
void TTTrackWriterT<T>::fillRoads(const std::vector<TTRoad>& roads) {
    if (T == kPatternMatcher) {
        vr_patternRef   ->clear();
        vr_tower        ->clear();
        vr_nstubs       ->clear();
        vr_patternInvPt ->clear();
        vr_patternFreq  ->clear();
        vr_superstripIds->clear();
        vr_stubRefs     ->clear();

        const unsigned nroads = roads.size();
        for (unsigned i=0; i<nroads; ++i) {
            const TTRoad& road = roads.at(i);
            vr_patternRef   ->push_back(road.patternRef);
            vr_tower        ->push_back(road.tower);
            vr_nstubs       ->push_back(road.nstubs);
            vr_patternInvPt ->push_back(road.patternInvPt);
            vr_patternFreq  ->push_back(road.patternFreq);
            vr_superstripIds->push_back(road.superstripIds);
            vr_stubRefs     ->push_back(road.stubRefs);
        }
        assert(vr_patternRef->size() == nroads);

        ttree_->Fill();
    }
}

template<size_t T>
void TTTrackWriterT<T>::fillRoads(const std::vector<TTRoad>& roads, const std::vector<std::string>& stubs_bitString, const std::vector<unsigned>& stubs_superstripId) {
    if (T == kPatternMatcher) {
        *vb_bitString = stubs_bitString;
        *vb_superstripId = stubs_superstripId;
        fillRoads(roads);
    }
}

template<size_t T>
void TTTrackWriterT<T>::fillTracks(const std::vector<TTTrack2>& tracks) {
    if (T == kTrackFitter) {
        //vt_px              ->clear();
        //vt_py              ->clear();
        //vt_pz              ->clear();
        vt_pt              ->clear();
        vt_eta             ->clear();
        //vt_phi             ->clear();
        //vt_vx              ->clear();
        //vt_vy              ->clear();
        //vt_vz              ->clear();
        vt_rinv            ->clear();
        vt_invPt           ->clear();
        vt_phi0            ->clear();
        vt_cottheta        ->clear();
        vt_z0              ->clear();
        vt_d0              ->clear();
        vt_chi2            ->clear();
        vt_ndof            ->clear();
        vt_synMatchChi2    ->clear();
        vt_synMatchCat     ->clear();
        vt_tpId            ->clear();
        vt_synTpId         ->clear();
        vt_tower           ->clear();
        vt_hitBits         ->clear();
        vt_ptSegment       ->clear();
        vt_roadRef         ->clear();
        vt_combRef         ->clear();
        vt_patternRef      ->clear();
        vt_stubRefs        ->clear();
        vt_principals      ->clear();

        const unsigned ntracks = tracks.size();
        for (unsigned i=0; i<ntracks; ++i) {
            const TTTrack2& track = tracks.at(i);

            //vt_px              ->push_back(track.px());
            //vt_py              ->push_back(track.py());
            //vt_pz              ->push_back(track.pz());
            vt_pt              ->push_back(track.pt());
            vt_eta             ->push_back(track.eta());
            //vt_phi             ->push_back(track.phi());
            //vt_vx              ->push_back(track.vx());
            //vt_vy              ->push_back(track.vy());
            //vt_vz              ->push_back(track.vz());
            vt_rinv            ->push_back(track.rinv());
            vt_invPt           ->push_back(track.invPt());
            vt_phi0            ->push_back(track.phi0());
            vt_cottheta        ->push_back(track.cottheta());
            vt_z0              ->push_back(track.z0());
            vt_d0              ->push_back(track.d0());
            vt_chi2            ->push_back(track.chi2());
            vt_ndof            ->push_back(track.ndof());
            vt_synMatchChi2    ->push_back(track.synMatchChi2());
            vt_synMatchCat     ->push_back(track.synMatchCat());
            vt_tpId            ->push_back(track.tpId());
            vt_synTpId         ->push_back(track.synTpId());
            vt_tower           ->push_back(track.tower());
            vt_hitBits         ->push_back(track.hitBits());
            vt_ptSegment       ->push_back(track.ptSegment());
            vt_roadRef         ->push_back(track.roadRef());
            vt_combRef         ->push_back(track.combRef());
            vt_patternRef      ->push_back(track.patternRef());
            vt_stubRefs        ->push_back(track.stubRefs());
            vt_principals      ->push_back(track.principals());
        }
        assert(vt_pt->size() == ntracks);

        ttree_->Fill();
    }
}

#endif  // __TTTrackWriter_h__
