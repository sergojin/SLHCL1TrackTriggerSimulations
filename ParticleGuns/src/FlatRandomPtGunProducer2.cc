/*
 *  $Date: 2009/02/19 21:52:40 $
 *  $Revision: 1.4 $
 *  \author Julia Yarba
 */

#include <ostream>

#include "SLHCL1TrackTriggerSimulations/ParticleGuns/interface/FlatRandomPtGunProducer2.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


using namespace edm;
using namespace std;

FlatRandomPtGunProducer2::FlatRandomPtGunProducer2(const ParameterSet& pset) : 
   BaseFlatGunProducer2(pset)
{


   ParameterSet defpset ;
   ParameterSet pgun_params = 
      pset.getParameter<ParameterSet>("PGunParameters") ;
  
   fMinPt = pgun_params.getParameter<double>("MinPt");
   fMaxPt = pgun_params.getParameter<double>("MaxPt");
   fXFlatSpread  = pgun_params.exists("XFlatSpread") ? pgun_params.getParameter<double>("XFlatSpread") : 0.;
   fYFlatSpread  = pgun_params.exists("YFlatSpread") ? pgun_params.getParameter<double>("YFlatSpread") : 0.;
   fZFlatSpread  = pgun_params.exists("ZFlatSpread") ? pgun_params.getParameter<double>("ZFlatSpread") : 0.;
   fRandomCharge = pgun_params.exists("RandomCharge")? pgun_params.getParameter<bool>  ("RandomCharge"): false;
  
  produces<HepMCProduct>();
  produces<GenEventInfoProduct>();
}

FlatRandomPtGunProducer2::~FlatRandomPtGunProducer2()
{
   // no need to cleanup GenEvent memory - done in HepMCProduct
}

void FlatRandomPtGunProducer2::produce(Event &e, const EventSetup& es) 
{

   if ( fVerbosity > 0 )
   {
      cout << " FlatRandomPtGunProducer2 : Begin New Event Generation" << endl ; 
   }
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
   HepMC::GenVertex* Vtx = new HepMC::GenVertex(HepMC::FourVector(xpos,ypos,zpos));

   // loop over particles
   //
   int barcode = 1 ;
   for (unsigned int ip=0; ip<fPartIDs.size(); ++ip)
   {

       double pt     = fRandomGenerator->fire(fMinPt, fMaxPt) ;
       double eta    = fRandomGenerator->fire(fMinEta, fMaxEta) ;
       double phi    = fRandomGenerator->fire(fMinPhi, fMaxPhi) ;
       int PartID = fPartIDs[ip] ;
       if (fRandomCharge && (fRandomGenerator->fire(0, 1) < 0.5))
           PartID = - PartID;
       const HepPDT::ParticleData* 
          PData = fPDGTable->particle(HepPDT::ParticleID(abs(PartID))) ;
       double mass   = PData->mass().value() ;
       double theta  = 2.*atan(exp(-eta)) ;
       double mom    = pt/sin(theta) ;
       double px     = pt*cos(phi) ;
       double py     = pt*sin(phi) ;
       double pz     = mom*cos(theta) ;
       double energy2= mom*mom + mass*mass ;
       double energy = sqrt(energy2) ; 
       HepMC::FourVector p(px,py,pz,energy) ;
       HepMC::GenParticle* Part = 
           new HepMC::GenParticle(p,PartID,1);
       Part->suggest_barcode( barcode ) ;
       barcode++ ;
       Vtx->add_particle_out(Part);

       if ( fAddAntiParticle )
       {
          HepMC::FourVector ap(-px,-py,-pz,energy) ;
	  int APartID = -PartID ;
	  if ( PartID == 22 || PartID == 23 )
	  {
	     APartID = PartID ;
	  }	  
	  HepMC::GenParticle* APart =
	     new HepMC::GenParticle(ap,APartID,1);
	  APart->suggest_barcode( barcode ) ;
	  barcode++ ;
	  Vtx->add_particle_out(APart) ;
       }

   }

   fEvt->add_vertex(Vtx) ;
   fEvt->set_event_number(e.id().event()) ;
   fEvt->set_signal_process_id(20) ; 
        
   if ( fVerbosity > 0 )
   {
      fEvt->print() ;  
   }

   auto_ptr<HepMCProduct> BProduct(new HepMCProduct()) ;
   BProduct->addHepMCData( fEvt );
   e.put(BProduct);

   auto_ptr<GenEventInfoProduct> genEventInfo(new GenEventInfoProduct(fEvt));
   e.put(genEventInfo);
    
   if ( fVerbosity > 0 )
   {
      // for testing purpose only
      // fEvt->print() ; // prints empty info after it's made into edm::Event
      cout << " FlatRandomPtGunProducer2 : Event Generation Done " << endl;
   }
}
//#include "FWCore/Framework/interface/MakerMacros.h"
//DEFINE_FWK_MODULE(FlatRandomPtGunProducer2);
