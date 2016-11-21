#ifndef __PatternBankWriter_h__
#define __PatternBankWriter_h__

#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include <memory>
#include <vector>

#include "BasicReader.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulationDataFormats/interface/Pattern.h"

// _____________________________________________________________________________
// This is a simple wrapper around TTree.

template<size_t T>
class PatternBankWriterT {
public:
    // Typedefs
    typedef slhcl1tt::superstrip_type superstrip_type;
    typedef slhcl1tt::frequency_type frequency_type;
    //typedef uint32_t superstrip_type;
    //typedef uint16_t frequency_type;

    PatternBankWriterT(int verbose=1);
    ~PatternBankWriterT();

    void init(TString out);

    void fillPatternMerging() { ttree4_->Fill(); ttree5_->Fill(); }

    void fillPatternAttributes() { ttree3_->Fill(); }

    void fillPatternBankInfo() { ttree2_->Fill(); }

    void fillPatternBank() { ttree1_->Fill(); }

    void write() { tfile_->Write(); }

    TTree * getTree()     { return ttree1_; }
    TTree * getInfoTree() { return ttree2_; }
    TTree * getAttrTree() { return ttree3_; }

    TTree * getTree_toMerged()   { return ttree4_; }
    TTree * getTree_fromMerged() { return ttree5_; }

    // Pattern merging indices
    std::unique_ptr<std::vector<unsigned> >        pb_indFromMerged;
    std::unique_ptr<std::vector<unsigned> >        pb_indToMerged;

     // Pattern attributes
    std::unique_ptr<float>                         pb_invPt_mean;
    std::unique_ptr<float>                         pb_invPt_sigma;
    std::unique_ptr<float>                         pb_cotTheta_mean;
    std::unique_ptr<float>                         pb_cotTheta_sigma;
    std::unique_ptr<float>                         pb_phi_mean;
    std::unique_ptr<float>                         pb_phi_sigma;
    std::unique_ptr<float>                         pb_z0_mean;
    std::unique_ptr<float>                         pb_z0_sigma;

    // Pattern bank statistics
    std::unique_ptr<float>                         pb_coverage;
    std::unique_ptr<unsigned>                      pb_count;
    std::unique_ptr<unsigned>                      pb_tower;
    std::unique_ptr<std::string>                   pb_superstrip;
    //std::unique_ptr<unsigned>                      pb_superstrip_nx;
    //std::unique_ptr<unsigned>                      pb_superstrip_nz;

    // Pattern bank
    std::unique_ptr<frequency_type>                pb_frequency;
    std::unique_ptr<std::vector<superstrip_type> > pb_superstripIds;

protected:
    TFile* tfile_;
    TTree* ttree1_;  // for pattern bank
    TTree* ttree2_;  // for pattern bank statistics
    TTree* ttree3_;  // for pattern attributes
    TTree* ttree4_;  // for merged pattern bank: toMerged
    TTree* ttree5_;  // for merged pattern bank: fromMerged
    const int verbose_;
};

// _____________________________________________________________________________
// Implementation is included in the header file to simplify ROOT library generation

#include <cassert>
#include <stdexcept>

template<size_t T>
PatternBankWriterT<T>::PatternBankWriterT(int verbose) :
    pb_indFromMerged  (new std::vector<unsigned>()),
    pb_indToMerged    (new std::vector<unsigned>()),
    //
    pb_invPt_mean     (new float(0.)),
    pb_invPt_sigma    (new float(0.)),
    pb_cotTheta_mean  (new float(0.)),
    pb_cotTheta_sigma (new float(0.)),
    pb_phi_mean       (new float(0.)),
    pb_phi_sigma      (new float(0.)),
    pb_z0_mean        (new float(0.)),
    pb_z0_sigma       (new float(0.)),
    //
    pb_coverage       (new float(0.)),
    pb_count          (new unsigned(0)),
    pb_tower          (new unsigned(0)),
    pb_superstrip     (new std::string("")),
    //pb_superstrip_nx  (new unsigned(0)),
    //pb_superstrip_nz  (new unsigned(0)),
    //
    pb_frequency      (new frequency_type(0)),
    pb_superstripIds  (new std::vector<superstrip_type>()),
    //
    tfile_(0), ttree1_(0), ttree2_(0), ttree3_(0), ttree4_(0), ttree5_(0), verbose_(verbose) {}

template<size_t T>
PatternBankWriterT<T>::~PatternBankWriterT() {
    if (ttree5_) delete ttree5_;
    if (ttree4_) delete ttree4_;
    if (ttree3_) delete ttree3_;
    if (ttree2_) delete ttree2_;
    if (ttree1_) delete ttree1_;
    if (tfile_)  delete tfile_;
}

template<size_t T>
void PatternBankWriterT<T>::init(TString out) {
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

    if (T == kPatternMerging) {
        // Pattern merging indices
        ttree5_ = new TTree("fromMerged", "");
        ttree5_->Branch("indFromMerged" , &(*pb_indFromMerged));

        ttree4_ = new TTree("toMerged", "");
        ttree4_->Branch("indToMerged"   , &(*pb_indToMerged));
    }

    if (T == kPatternGenerator) {
        // Pattern attributes
        ttree3_ = new TTree("patternAttributes", "");
        ttree3_->Branch("invPt_mean"    , &(*pb_invPt_mean));
        ttree3_->Branch("invPt_sigma"   , &(*pb_invPt_sigma));
        ttree3_->Branch("cotTheta_mean" , &(*pb_cotTheta_mean));
        ttree3_->Branch("cotTheta_sigma", &(*pb_cotTheta_sigma));
        ttree3_->Branch("phi_mean"      , &(*pb_phi_mean));
        ttree3_->Branch("phi_sigma"     , &(*pb_phi_sigma));
        ttree3_->Branch("z0_mean"       , &(*pb_z0_mean));
        ttree3_->Branch("z0_sigma"      , &(*pb_z0_sigma));

        // Pattern bank statistics
        ttree2_ = new TTree("patternBankInfo", "");
        ttree2_->Branch("coverage"      , &(*pb_coverage));
        ttree2_->Branch("count"         , &(*pb_count));
        ttree2_->Branch("tower"         , &(*pb_tower));
        ttree2_->Branch("superstrip"    , &(*pb_superstrip));

        // Pattern bank
        ttree1_ = new TTree("patternBank", "");
        ttree1_->Branch("frequency"      , &(*pb_frequency));
        ttree1_->Branch("superstripIds"  , &(*pb_superstripIds));
    }
}

#endif  // __PatternBankWriter_h__
