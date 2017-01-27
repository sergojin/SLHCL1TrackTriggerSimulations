#ifndef __TTStubWriter_h__
#define __TTStubWriter_h__

#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include <memory>
#include <vector>

#include "BasicReader.h"

// _____________________________________________________________________________
// This is a simple wrapper around TTree. It clones a tree and adds new
// branches if needed

template<size_t T>
class TTStubWriterT {
public:
    TTStubWriterT(int verbose=1);
    ~TTStubWriterT();

    void init(TTree* intree, TString out);

    void fill() { ttree_->Fill(); }

    void write() { tfile_->Write(); }

    TTree * getTree() { return ttree_; }

protected:
    TFile* tfile_;
    TTree* ttree_;
    const int verbose_;
};

// _____________________________________________________________________________
// Implementation is included in the header file to simplify ROOT library generation

#include <cassert>
#include <stdexcept>

template<size_t T>
TTStubWriterT<T>::TTStubWriterT(int verbose) :
    tfile_(0), ttree_(0), verbose_(verbose) {}

template<size_t T>
TTStubWriterT<T>::~TTStubWriterT() {
    if (ttree_)  delete ttree_;
    if (tfile_)  delete tfile_;
}

template<size_t T>
void TTStubWriterT<T>::init(TTree* intree, TString out) {
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
    //ttree_ = (TTree*) intree->GetTree()->CloneTree(0); // Do not copy the data yet
    if (ttree_ == 0) {
        TString msg = "Failed to clone the tree";
        throw std::invalid_argument(msg.Data());
    }

}

#endif  // __TTStubWriter_h__
