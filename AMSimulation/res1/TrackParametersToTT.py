from math import pi, asin, sinh, asinh

def get_phiStar_from_phi(phi, invPt, rStar=90.):
    mPtFactor = 0.3*3.8*1e-2/2.0
    dphi = - asin(mPtFactor * rStar * invPt)
    phiStar = phi + dphi
    while phiStar < -pi:
      phiStar += pi*2.
    while phiStar >= pi:
      phiStar -= pi*2.
    return phiStar

def get_etaStar_from_eta(eta, z0, invPt, rStar=60.):
    mPtFactor = 0.3*3.8*1e-2/2.0
    if abs(rStar) < 1e-10: return eta
    if abs(invPt) < 1e-10 and invPt <  0.: invPt = -1e-10
    if abs(invPt) < 1e-10 and invPt >= 0.: invPt = +1e-10
    cot = sinh(eta)
    cotStar = (cot * (asin(mPtFactor * rStar * invPt)/(mPtFactor * invPt)) + z0) / rStar
    return asinh(cotStar)

def TrackParametersToTT(phi, invPt, eta, z0):
    max_eta = 2.2
    max_z0 = 15.      # [cm]
    max_invPt = 1./3  # [1/GeV]
    etaStar = get_etaStar_from_eta(eta, z0, invPt)
    phiStar = get_phiStar_from_phi(phi, invPt)
    if abs(etaStar) > max_eta or abs(z0) > max_z0:
      return -1

    tt_eta = int((etaStar + max_eta) / (max_eta*2./6))
    tt_phi = int((phiStar + pi) / (pi*2./8))
    tt = tt_eta * 8 + tt_phi
    assert(0<=tt and tt<48)
    return tt
