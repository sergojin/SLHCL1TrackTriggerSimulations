#ifndef __PatternBankReader_h__
#define __PatternBankReader_h__

#include "TFile.h"
#include "TFileCollection.h"
#include "TFileInfo.h"
#include "TString.h"
#include "TTree.h"
#include <memory>
#include <vector>

#include "BasicReader.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulationDataFormats/interface/Pattern.h"

// _____________________________________________________________________________
// This is a simple wrapper around multiple TTrees. It sets the branch names, addresses
// etc. Its functions are essentially the same as the functions of TTree.

template<size_t T>
class PatternBankReaderT {
public:
    // Typedefs
    typedef slhcl1tt::superstrip_type superstrip_type;
    typedef slhcl1tt::frequency_type frequency_type;
    //typedef uint32_t superstrip_type;
    //typedef uint16_t frequency_type;

    PatternBankReaderT(int verbose=1);
    ~PatternBankReaderT();

    void init(TString src);

    Int_t getPattern    (Long64_t entry)   { return ttree1_->GetEntry(entry); }
    Int_t getPatternInfo(Long64_t entry=0) { return ttree2_->GetEntry(0);     }  // only 1 entry
    Int_t getPatternAttr(Long64_t entry)   { return ttree3_->GetEntry(entry); }

    Int_t getEntry_toMerged  (Long64_t entry=0) { return ttree4_->GetEntry(0);     }  // only 1 entry
    Int_t getEntry_fromMerged(Long64_t entry)   { return ttree5_->GetEntry(entry); }

    void getPatternInvPt(Long64_t entry, float& invPt_mean) {
        ttree3_->GetEntry(entry);
        invPt_mean = pb_invPt_mean;
    }

    Long64_t getEntries() const { return ttree1_->GetEntries(); }

    TTree * getTree()     { return ttree1_; }
    TTree * getInfoTree() { return ttree2_; }
    TTree * getAttrTree() { return ttree3_; }

    TTree * getTree_toMerged()   { return ttree4_; }
    TTree * getTree_fromMerged() { return ttree5_; }

    // Pattern merging indices
    std::vector<unsigned> *        pb_indFromMerged;
    std::vector<unsigned> *        pb_indToMerged;

    // Pattern attributes
    float                          pb_invPt_mean;
    float                          pb_invPt_sigma;
    float                          pb_cotTheta_mean;
    float                          pb_cotTheta_sigma;
    float                          pb_phi_mean;
    float                          pb_phi_sigma;
    float                          pb_z0_mean;
    float                          pb_z0_sigma;

    // Pattern bank statistics
    float                          pb_coverage;
    unsigned                       pb_count;
    unsigned                       pb_tower;
    std::string *                  pb_superstrip;
    //unsigned                       pb_superstrip_nx;
    //unsigned                       pb_superstrip_nz;

    // Pattern bank
    frequency_type                 pb_frequency;
    std::vector<superstrip_type> * pb_superstripIds;

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
PatternBankReaderT<T>::PatternBankReaderT(int verbose) :
    pb_indFromMerged  (0),
    pb_indToMerged    (0),
    //
    pb_invPt_mean     (0.),
    pb_invPt_sigma    (0.),
    pb_cotTheta_mean  (0.),
    pb_cotTheta_sigma (0.),
    pb_phi_mean       (0.),
    pb_phi_sigma      (0.),
    pb_z0_mean        (0.),
    pb_z0_sigma       (0.),
    //
    pb_coverage       (0.),
    pb_count          (0),
    pb_tower          (0),
    pb_superstrip     (0),
    //pb_superstrip_nx  (0),
    //pb_superstrip_nz  (0),
    //
    pb_frequency      (0),
    pb_superstripIds  (0),
    //
    tfile_(0), ttree1_(0), ttree2_(0), ttree3_(0), ttree4_(0), ttree5_(0), verbose_(verbose) {}

template<size_t T>
PatternBankReaderT<T>::~PatternBankReaderT() {
    if (ttree5_) delete ttree5_;
    if (ttree4_) delete ttree4_;
    if (ttree3_) delete ttree3_;
    if (ttree2_) delete ttree2_;
    if (ttree1_) delete ttree1_;
    if (tfile_)  delete tfile_;
}

template<size_t T>
void PatternBankReaderT<T>::init(TString src) {
    if (!src.EndsWith(".root") && !src.EndsWith(".txt")) {
        TString msg = "Input source must be either .root or .txt";
        throw std::invalid_argument(msg.Data());
    }

    //if (verbose_)  std::cout << "Opening " << src << std::endl;

    if (src.EndsWith(".txt")) {
        TFileCollection fc("fileinfolist", "", src);

        TCollection* filelist = (TCollection*) fc.GetList();
        TIter next(filelist);
        TFileInfo* fi = (TFileInfo*) next();
        src = (fi->GetCurrentUrl()) ? fi->GetCurrentUrl()->GetUrl() : 0;
        if (!src.EndsWith(".root")) {
            TString msg = "Input source must be .root";
            throw std::invalid_argument(msg.Data());
        }
    }

    tfile_ = TFile::Open(src);

    if (tfile_) {
        //if (verbose_)  std::cout << "Successfully read " << src << std::endl;
    } else {
        TString msg = "Failed to read " + src;
        throw std::invalid_argument(msg.Data());
    }

    ttree5_ = (TTree*) tfile_->Get("fromMerged");
    if (ttree5_ != 0) {
        ttree5_->SetBranchAddress("indFromMerged", &(pb_indFromMerged));
    } else {
        //std::cout << "Cannot find TTree \"fromMerged\"" << std::endl;
    }

    ttree4_ = (TTree*) tfile_->Get("toMerged");
    if (ttree4_ != 0) {
        ttree4_->SetBranchAddress("indToMerged", &(pb_indToMerged));
    } else {
        //std::cout << "Cannot find TTree \"toMerged\"" << std::endl;
    }

    ttree3_ = (TTree*) tfile_->Get("patternAttributes");
    if (ttree3_ != 0) {
        if (T == kPatternMatcher) {
            ttree3_->SetBranchStatus("*"              , 0);
            ttree3_->SetBranchStatus("invPt_mean"     , 1);
        }
        ttree3_->SetBranchAddress("invPt_mean"    , &(pb_invPt_mean));
        ttree3_->SetBranchAddress("invPt_sigma"   , &(pb_invPt_sigma));
        ttree3_->SetBranchAddress("cotTheta_mean" , &(pb_cotTheta_mean));
        ttree3_->SetBranchAddress("cotTheta_sigma", &(pb_cotTheta_sigma));
        ttree3_->SetBranchAddress("phi_mean"      , &(pb_phi_mean));
        ttree3_->SetBranchAddress("phi_sigma"     , &(pb_phi_sigma));
        ttree3_->SetBranchAddress("z0_mean"       , &(pb_z0_mean));
        ttree3_->SetBranchAddress("z0_sigma"      , &(pb_z0_sigma));
    } else {
        //std::cout << "Cannot find TTree \"patternAttributes\"" << std::endl;
    }

    ttree2_ = (TTree*) tfile_->Get("patternBankInfo");
    assert(ttree2_ != 0);

    ttree2_->SetBranchAddress("coverage"      , &(pb_coverage));
    ttree2_->SetBranchAddress("count"         , &(pb_count));
    ttree2_->SetBranchAddress("tower"         , &(pb_tower));
    ttree2_->SetBranchAddress("superstrip"    , &(pb_superstrip));
    //ttree2_->SetBranchAddress("superstrip_nx" , &(pb_superstrip_nx));
    //ttree2_->SetBranchAddress("superstrip_nz" , &(pb_superstrip_nz));

    ttree1_ = (TTree*) tfile_->Get("patternBank");
    assert(ttree1_ != 0);

    ttree1_->SetBranchAddress("frequency"     , &(pb_frequency));
    ttree1_->SetBranchAddress("superstripIds" , &(pb_superstripIds));
}

#endif  // __PatternBankReader_h__
