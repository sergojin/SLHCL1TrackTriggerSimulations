#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/TrackParametersToTT.h"
using namespace slhcl1tt;

#include <cmath>
#include <cassert>


int TrackParametersToTT::get_tt(double phi, double invPt, double eta, double vz) {
  const double max_eta = 2.2;
  const double max_vz = 15.;
  const double max_invPt = 1./3;

  double etaStar = get_etaStar_from_eta(eta, vz, invPt);
  double phiStar = get_phiStar_from_phi(phi, invPt);

  if (std::abs(etaStar) > max_eta || std::abs(vz) > max_vz || std::abs(invPt) > max_invPt) {
    return -1;
  }

  // Bring phi -> [0, 2pi]
  while (phiStar > M_PI*2.) {
    phiStar -= M_PI*2.;
  }
  while (phiStar <= 0.) {
    phiStar += M_PI*2.;
  }

  // Assign trigger tower
  int tt_eta = (etaStar + max_eta) / (max_eta*2./6);
  int tt_phi = phiStar / (M_PI*2./8);
  int tt = tt_eta * 8 + tt_phi;
  assert(tt<48);
  return tt;
}

double TrackParametersToTT::get_phiStar_from_phi(double phi, double invPt, double rStar) {
  constexpr double mPtFactor = 0.3*3.8*1e-2/2.0;
  double dphi = - asin(mPtFactor * rStar * invPt);
  return phi + dphi;
}

double TrackParametersToTT::get_phi_from_phiStar(double phiStar, double invPt, double rStar) {
  constexpr double mPtFactor = 0.3*3.8*1e-2/2.0;
  double dphi = - asin(mPtFactor * rStar * invPt);
  return phiStar - dphi;
}

double TrackParametersToTT::get_etaStar_from_eta(double eta, double vz, double invPt, double rStar) {
  constexpr double mPtFactor = 0.3*3.8*1e-2/2.0;
  if (std::abs(invPt) < 1e-10 && invPt <  0.) invPt = -1e-10;
  if (std::abs(invPt) < 1e-10 && invPt >= 0.) invPt = +1e-10;
  double cot = sinh(eta);
  double cotStar = (cot * (asin(mPtFactor * rStar * invPt)/(mPtFactor * invPt)) + vz) / rStar;
  return asinh(cotStar);
}

double TrackParametersToTT::get_eta_from_etaStar(double etaStar, double vz, double invPt, double rStar) {
  constexpr double mPtFactor = 0.3*3.8*1e-2/2.0;
  if (std::abs(invPt) < 1e-10 && invPt <  0.) invPt = -1e-10;
  if (std::abs(invPt) < 1e-10 && invPt >= 0.) invPt = +1e-10;
  double cotStar = sinh(etaStar);
  double cot = (rStar * cotStar - vz) / (asin(mPtFactor * rStar * invPt)/(mPtFactor * invPt));
  return asinh(cot);
}
