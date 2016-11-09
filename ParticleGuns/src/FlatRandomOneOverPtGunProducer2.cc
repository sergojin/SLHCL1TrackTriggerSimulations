#include <ostream>

#include "SLHCL1TrackTriggerSimulations/ParticleGuns/interface/FlatRandomOneOverPtGunProducer2.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "CLHEP/Units/GlobalSystemOfUnits.h"
using namespace edm;

namespace {
  double get_phi_from_phiStar(double phiStar, double invPt, double rStar=90.) {
    constexpr double mPtFactor = 0.3*3.8*1e-2/2.0;
    double dphi = - asin(mPtFactor * rStar * invPt);
    return phiStar - dphi;
  }

  double get_eta_from_etaStar(double etaStar, double vz, double invPt, double rStar=60.) {
    constexpr double mPtFactor = 0.3*3.8*1e-2/2.0;
    if (std::abs(invPt) < 1e-10 && invPt <  0.) invPt = -1e-10;
    if (std::abs(invPt) < 1e-10 && invPt >= 0.) invPt = +1e-10;
    double cotStar = sinh(etaStar);
    double cot = (rStar * cotStar - vz) / (asin(mPtFactor * rStar * invPt)/(mPtFactor * invPt));
    return asinh(cot);
  }
}

FlatRandomOneOverPtGunProducer2::FlatRandomOneOverPtGunProducer2(const edm::ParameterSet& pset) : 
  BaseFlatGunProducer2(pset) {


  edm::ParameterSet defpset ;
  edm::ParameterSet pgun_params = 
    pset.getParameter<ParameterSet>("PGunParameters") ;
  
  fMinOneOverPt = pgun_params.getParameter<double>("MinOneOverPt");
  fMaxOneOverPt = pgun_params.getParameter<double>("MaxOneOverPt");
  fXFlatSpread  = pgun_params.exists("XFlatSpread") ? pgun_params.getParameter<double>("XFlatSpread") : 0.;
  fYFlatSpread  = pgun_params.exists("YFlatSpread") ? pgun_params.getParameter<double>("YFlatSpread") : 0.;
  fZFlatSpread  = pgun_params.exists("ZFlatSpread") ? pgun_params.getParameter<double>("ZFlatSpread") : 0.;
  fRandomCharge = pgun_params.exists("RandomCharge")? pgun_params.getParameter<bool>  ("RandomCharge"): false;
  fReallyFlat   = pgun_params.exists("ReallyFlat")  ? pgun_params.getParameter<bool>  ("ReallyFlat")  : false;
  fVtxSmeared   = pgun_params.exists("VtxSmeared")  ? pgun_params.getParameter<bool>  ("VtxSmeared")  : false;
  fUseRStar     = pgun_params.exists("UseRStar")    ? pgun_params.getParameter<bool>  ("UseRStar")    : false;

  
  produces<HepMCProduct>();
  produces<GenEventInfoProduct>();

  edm::LogInfo("ParticleGun") << "FlatRandomOneOverPtGunProducer2: initialized with minimum and maximum 1/pt " << fMinOneOverPt << ":" << fMaxOneOverPt;
}

FlatRandomOneOverPtGunProducer2::~FlatRandomOneOverPtGunProducer2() {
  // no need to cleanup GenEvent memory - done in HepMCProduct
}

void FlatRandomOneOverPtGunProducer2::produce(Event &e, const EventSetup& es) {

  LogDebug("ParticleGun") << " FlatRandomOneOverPtGunProducer2 : Begin New Event Generation"; 

  // event loop (well, another step in it...)
          
  // no need to clean up GenEvent memory - done in HepMCProduct
  // 
   
  // here re-create fEvt (memory)
  //
  fEvt = new HepMC::GenEvent() ;
   
  // now actualy, cook up the event from PDGTable and gun parameters
  //
  // 1st, primary vertex
  //
  double xpos = fXFlatSpread > 0. ? fRandomGenerator->fire(-fXFlatSpread,fXFlatSpread) : 0.;
  double ypos = fYFlatSpread > 0. ? fRandomGenerator->fire(-fYFlatSpread,fYFlatSpread) : 0.;
  double zpos = fZFlatSpread > 0. ? fRandomGenerator->fire(-fZFlatSpread,fZFlatSpread) : 0.;
  if (fVtxSmeared)
         zpos = fRandomGaussGenerator->fire(0,4.5*10.);
  HepMC::GenVertex* Vtx = new HepMC::GenVertex(HepMC::FourVector(xpos,ypos,zpos));

  // loop over particles
  //
  int barcode = 1 ;
  for (unsigned int ip=0; ip<fPartIDs.size(); ++ip) {

    double xx     = fRandomGenerator->fire(0.0,1.0);
    double pt     = std::exp((1.-xx)*std::log(fMinOneOverPt)+
			     xx*std::log(fMaxOneOverPt)) ;
    if (fReallyFlat)
           pt     = fMinOneOverPt + xx * (fMaxOneOverPt - fMinOneOverPt);
    double eta    = fRandomGenerator->fire(fMinEta, fMaxEta) ;
    double phi    = fRandomGenerator->fire(fMinPhi, fMaxPhi) ;
    if (pt != 0) pt = 1./pt;
    int PartID = fPartIDs[ip] ;
    if (fRandomCharge && (fRandomGenerator->fire(0, 1) < 0.5))
        PartID = - PartID;
    const HepPDT::ParticleData* 
      PData = fPDGTable->particle(HepPDT::ParticleID(PartID)) ;
    double mass   = PData->mass().value() ;

    if (fUseRStar)
           phi    = get_phi_from_phiStar(phi, PData->charge()/pt, 90.);  // fRStarForPhi = 90.
    if (fUseRStar)
           eta    = get_eta_from_etaStar(eta, zpos/10., PData->charge()/pt, 60.);  // fRStarForEta = 60.

    double theta  = 2.*atan(exp(-eta)) ;
    double mom    = pt/sin(theta) ;
    double px     = pt*cos(phi) ;
    double py     = pt*sin(phi) ;
    double pz     = mom*cos(theta) ;
    double energy2= mom*mom + mass*mass ;
    double energy = sqrt(energy2) ; 
    HepMC::FourVector p(px,py,pz,energy) ;
    HepMC::GenParticle* Part = new HepMC::GenParticle(p,PartID,1);
    Part->suggest_barcode( barcode ) ;
    barcode++ ;
    Vtx->add_particle_out(Part);
    LogDebug("ParticleGun") << "FlatRandomOneOverPtGunProducer2: Event generated with pt:eta:phi " << pt << " " << eta << " " << phi << " (" << theta/CLHEP::deg << ":" << phi/CLHEP::deg << ")";

    if ( fAddAntiParticle ) {
      HepMC::FourVector ap(-px,-py,-pz,energy) ;
      int APartID = -PartID ;
      if ( PartID == 22 || PartID == 23 ) {
	APartID = PartID ;
      }	  
      HepMC::GenParticle* APart = new HepMC::GenParticle(ap,APartID,1);
      APart->suggest_barcode( barcode ) ;
      barcode++ ;
      Vtx->add_particle_out(APart) ;
    }

  }

  fEvt->add_vertex(Vtx) ;
  fEvt->set_event_number(e.id().event()) ;
  fEvt->set_signal_process_id(20) ; 
        
  if ( fVerbosity > 0 ) {
    fEvt->print() ;  
  }

  std::auto_ptr<HepMCProduct> BProduct(new HepMCProduct()) ;
  BProduct->addHepMCData( fEvt );
  e.put(BProduct);

  std::auto_ptr<GenEventInfoProduct> genEventInfo(new GenEventInfoProduct(fEvt));
  e.put(genEventInfo);
    
  LogDebug("ParticleGun") << " FlatRandomOneOverPtGunProducer2 : Event Generation Done ";
}
