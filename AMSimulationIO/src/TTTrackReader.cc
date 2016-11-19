#include "SLHCL1TrackTriggerSimulations/AMSimulationIO/interface/TTTrackReader.h"

#include "SLHCL1TrackTriggerSimulations/AMSimulationIO/interface/Helper.h"
using namespace slhcl1tt;


// _____________________________________________________________________________
TTTrackReader::TTTrackReader(int verbose)
: TTRoadReader(verbose) {}

TTTrackReader::~TTTrackReader() {}

int TTTrackReader::init(TString src, TString prefixRoad, TString prefixTrack, TString suffix) {
    if (TTRoadReader::init(src, prefixRoad, suffix))
        return 1;

    return 0;
}


// _____________________________________________________________________________
TTTrackWriter::TTTrackWriter(int verbose)
: BasicWriter(verbose),

  vt_px           (new std::vector<float>()),
  vt_py           (new std::vector<float>()),
  vt_pz           (new std::vector<float>()),
  vt_pt           (new std::vector<float>()),
  vt_eta          (new std::vector<float>()),
  vt_phi          (new std::vector<float>()),
  vt_vx           (new std::vector<float>()),
  vt_vy           (new std::vector<float>()),
  vt_vz           (new std::vector<float>()),
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
  vt_principals   (new std::vector<std::vector<float> >()) {}


TTTrackWriter::~TTTrackWriter() {}

int TTTrackWriter::init(TChain* tchain, TString out, TString prefix, TString suffix) {
    if (BasicWriter::init(tchain, out))
        return 1;

  //ttree->Branch(prefix + "px"             + suffix, &(*vt_px));
  //ttree->Branch(prefix + "py"             + suffix, &(*vt_py));
  //ttree->Branch(prefix + "pz"             + suffix, &(*vt_pz));
    ttree->Branch(prefix + "pt"             + suffix, &(*vt_pt));
    ttree->Branch(prefix + "eta"            + suffix, &(*vt_eta));
  //ttree->Branch(prefix + "phi"            + suffix, &(*vt_phi));
  //ttree->Branch(prefix + "vx"             + suffix, &(*vt_vx));
  //ttree->Branch(prefix + "vy"             + suffix, &(*vt_vy));
  //ttree->Branch(prefix + "vz"             + suffix, &(*vt_vz));
    ttree->Branch(prefix + "rinv"           + suffix, &(*vt_rinv));
    ttree->Branch(prefix + "invPt"          + suffix, &(*vt_invPt));
    ttree->Branch(prefix + "phi0"           + suffix, &(*vt_phi0));
    ttree->Branch(prefix + "cottheta"       + suffix, &(*vt_cottheta));
    ttree->Branch(prefix + "z0"             + suffix, &(*vt_z0));
    ttree->Branch(prefix + "d0"             + suffix, &(*vt_d0));
    ttree->Branch(prefix + "chi2"           + suffix, &(*vt_chi2));
    ttree->Branch(prefix + "ndof"           + suffix, &(*vt_ndof));
    ttree->Branch(prefix + "synMatchChi2"   + suffix, &(*vt_synMatchChi2));
    ttree->Branch(prefix + "synMatchCat"    + suffix, &(*vt_synMatchCat));
    ttree->Branch(prefix + "tpId"           + suffix, &(*vt_tpId));
    ttree->Branch(prefix + "synTpId"        + suffix, &(*vt_synTpId));
    ttree->Branch(prefix + "tower"          + suffix, &(*vt_tower));
    ttree->Branch(prefix + "hitBits"        + suffix, &(*vt_hitBits));
    ttree->Branch(prefix + "ptSegment"      + suffix, &(*vt_ptSegment));
    ttree->Branch(prefix + "roadRef"        + suffix, &(*vt_roadRef));
    ttree->Branch(prefix + "combRef"        + suffix, &(*vt_combRef));
    ttree->Branch(prefix + "patternRef"     + suffix, &(*vt_patternRef));
    ttree->Branch(prefix + "stubRefs"       + suffix, &(*vt_stubRefs));
    ttree->Branch(prefix + "principals"     + suffix, &(*vt_principals));
    return 0;
}

void TTTrackWriter::fill(const std::vector<TTTrack2>& tracks) {
    vt_px              ->clear();
    vt_py              ->clear();
    vt_pz              ->clear();
    vt_pt              ->clear();
    vt_eta             ->clear();
    vt_phi             ->clear();
    vt_vx              ->clear();
    vt_vy              ->clear();
    vt_vz              ->clear();
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

    ttree->Fill();
    assert(vt_pt->size() == ntracks);
}
