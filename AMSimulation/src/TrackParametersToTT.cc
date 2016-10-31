#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/TrackParametersToTT.h"
using namespace slhcl1tt;

#include <cmath>
#include <cassert>


int TrackParametersToTT::get_tt(double phi, double invPt, double eta, double z0) {
  const double etaStar_max = 2.2;
  const double z0_max = 15.;
  const double invPt_max = 1./3;

  double etaStar = get_etaStar_from_eta(eta, z0, invPt);
  double phiStar = get_phiStar_from_phi(phi, invPt);

  if (std::abs(etaStar) > etaStar_max || std::abs(z0) > z0_max || std::abs(invPt) > invPt_max) {
    return -1;
  }

  // Bring phi -> [0, 2pi]
  while (phiStar > M_PI*2.) {
    phiStar -= M_PI*2.;
  }
  while (phiStar <= 0.) {
    phiStar += M_PI*2.;
  }

  int tt_phi = phiStar / (M_PI*2./8);
  int tt_eta = (etaStar + etaStar_max) / (etaStar_max*2./6);
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

double TrackParametersToTT::get_etaStar_from_eta(double eta, double z0, double invPt, double rStar) {
  constexpr double mPtFactor = 0.3*3.8*1e-2/2.0;
  if (std::abs(invPt) < 1e-10 && invPt <  0.) invPt = -1e-10;
  if (std::abs(invPt) < 1e-10 && invPt >= 0.) invPt = +1e-10;
  double cot = sinh(eta);
  double cotStar = (cot * (asin(mPtFactor * rStar * invPt)/(mPtFactor * invPt)) + z0) / rStar;
  return asinh(cotStar);
}

double TrackParametersToTT::get_eta_from_etaStar(double etaStar, double z0, double invPt, double rStar) {
  constexpr double mPtFactor = 0.3*3.8*1e-2/2.0;
  if (std::abs(invPt) < 1e-10 && invPt <  0.) invPt = -1e-10;
  if (std::abs(invPt) < 1e-10 && invPt >= 0.) invPt = +1e-10;
  double cotStar = sinh(etaStar);
  double cot = (rStar * cotStar - z0) / (asin(mPtFactor * rStar * invPt)/(mPtFactor * invPt));
  return asinh(cot);
}
