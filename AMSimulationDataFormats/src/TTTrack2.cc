#include "SLHCL1TrackTriggerSimulations/AMSimulationDataFormats/interface/TTTrack2.h"
#include <iostream>


namespace slhcl1tt {

std::ostream& operator<<(std::ostream& o, const TTTrack2& track) {
    o << "rinv: " << track.rinv() << " phi0: " << track.phi0() << " cottheta: " << track.cottheta() << " z0: " << track.z0() << " d0: " << track.d0()
      << " chi2: " << track.chi2() << " ndof: " << track.ndof() << " synMatchChi2: " << track.synMatchChi2() << " synMatchCat: " << track.synMatchCat()
      << " tpId: " << track.tpId() << " tower: " << track.tower() << " hitBits: " << track.hitBits() << " ptSegment: " << track.ptSegment()
      << " roadRef: " << track.roadRef() << " combRef: " << track.combRef() << " patternRef: " << track.patternRef()
      << " # stubs: " << track.stubRefs().size() << " # principals: " << track.principals().size()
      << std::endl;
    return o;
}

}  // namespace slhcl1tt
