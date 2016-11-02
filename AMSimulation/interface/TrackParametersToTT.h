#ifndef AMSimulation_TrackParametersToTT_h_
#define AMSimulation_TrackParametersToTT_h_

// Codes originally written by Olmo Cerri (SNS Pisa)
// See https://github.com/ocerri/SLHCL1TrackTriggerSimulations/blob/dev_tt25/AMSimulation/src/TrackParametersToTT.h
// Modified for better coherence with the rest of AMSimulation

namespace slhcl1tt {

  class TrackParametersToTT {
  public:
    int get_tt(double phi, double invPt, double eta, double vz);

    double get_phiStar_from_phi(double phi, double invPt, double rStar=90.);
    double get_phi_from_phiStar(double phiStar, double invPt, double rStar=90.);

    double get_etaStar_from_eta(double eta, double vz, double invPt, double rStar=60.);
    double get_eta_from_etaStar(double etaStar, double vz, double invPt, double rStar=60.);
  };

}  // namespace slhcl1tt

#endif
