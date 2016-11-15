#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/AssociativeMemory.h"
using namespace slhcl1tt;

#include <cassert>
#include <iostream>
#include <stdexcept>


// _____________________________________________________________________________
int AssociativeMemory::init(unsigned npatterns) {
    patternBank_.clear();
    patternBank_.reserve(npatterns);

    patternAttributes_invPt_.clear();
    patternAttributes_invPt_.reserve(npatterns);

    patternAttributes_freq_.clear();
    patternAttributes_freq_.reserve(npatterns);

    return 0;
}

// _____________________________________________________________________________
void AssociativeMemory::insert(std::vector<superstrip_type>::const_iterator begin, std::vector<superstrip_type>::const_iterator end, const float invPt, const unsigned freq) {
    //patternBank_.insert(patternBank_.end(), begin, end);

    pattern_type patt;
    unsigned i = 0;
    for (std::vector<superstrip_type>::const_iterator it = begin; it != end; ++it, ++i) {
        patt.at(i) = *it;
    }
    patternBank_.push_back(patt);
    patternAttributes_invPt_.push_back(invPt);
    patternAttributes_freq_.push_back(freq);
}

void AssociativeMemory::insert(const pattern_type& patt, const float invPt, const unsigned freq) {
    patternBank_.push_back(patt);
    patternAttributes_invPt_.push_back(invPt);
    patternAttributes_freq_.push_back(freq);
}

// _____________________________________________________________________________
void AssociativeMemory::freeze() {
    assert(patternBank_.size() == patternAttributes_invPt_.size());
    assert(patternBank_.size() == patternAttributes_freq_.size());
    frozen_ = true;
}

// _____________________________________________________________________________
std::vector<unsigned> AssociativeMemory::lookup(const HitBuffer& hitBuffer, const unsigned nLayers, const unsigned maxMisses) {
    std::vector<unsigned> firedPatterns;

    for (std::vector<pattern_type>::const_iterator itpatt = patternBank_.begin();
         itpatt != patternBank_.end(); ++itpatt) {
        unsigned nMisses = 0;

        for (pattern_type::const_reverse_iterator itlayer = itpatt->rend() - nLayers;
             itlayer != itpatt->rend(); ++itlayer) {

            const superstrip_type ss = *itlayer;
            if (!hitBuffer.isHit(ss))
                ++nMisses;

            // Skip if more misses than allowed
            if (nMisses > maxMisses)
                break;
        }
        if (nMisses <= maxMisses)
            firedPatterns.push_back(itpatt - patternBank_.begin());
    }
    return firedPatterns;
}

// _____________________________________________________________________________
void AssociativeMemory::retrieve(const unsigned patternRef, pattern_type& superstripIds, float& invPt, unsigned& freq) {
    superstripIds = patternBank_            .at(patternRef);
    invPt         = patternAttributes_invPt_.at(patternRef);
    freq          = patternAttributes_freq_ .at(patternRef);
}

// _____________________________________________________________________________
void AssociativeMemory::print() {
    std::cout << "npatterns: " << patternBank_.size() << std::endl;
}
