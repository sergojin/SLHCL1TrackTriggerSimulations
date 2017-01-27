#ifndef NTupleTools_SimpleSkimmer_h_
#define NTupleTools_SimpleSkimmer_h_

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"

#include "TrackParametersToTT.h"


class SimpleSkimmer : public edm::EDFilter {
  public:
    explicit SimpleSkimmer(const edm::ParameterSet&);

  private:
    //virtual void beginJob();
    virtual bool filter(edm::Event&, const edm::EventSetup&);
    //virtual void endJob();

    const edm::InputTag inputTag_;

    StringCutObjectSelector<TrackingParticle> selector_;
    const unsigned minN_;

    TH1F * h_denom, * h_num;
};
#endif


// _____________________________________________________________________________
#include "TH1F.h"

SimpleSkimmer::SimpleSkimmer(const edm::ParameterSet& iConfig) :
  inputTag_(iConfig.getParameter<edm::InputTag>("inputTag")),
  selector_(iConfig.getParameter<std::string>("cut")),
  minN_    (iConfig.getParameter<unsigned>("minN"))
{
    edm::Service<TFileService> fs;
    h_denom = fs->make<TH1F>("denom", "", 2, 0, 2);
    h_num   = fs->make<TH1F>("num", "", 2, 0, 2);
}

bool SimpleSkimmer::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    bool result = false;

    slhcl1tt::TrackParametersToTT tt_operator;
    const int the_tt = 25;

    if (!iEvent.isRealData()) {
        edm::Handle<TrackingParticleCollection> parts;
        iEvent.getByLabel(inputTag_, parts);

        if (parts.isValid()) {
            edm::LogInfo("SimpleSkimmer") << "Size: " << parts->size();

            unsigned n = 0;
            for (TrackingParticleCollection::const_iterator it = parts->begin(); it != parts->end(); ++it) {
                if (n >= minN_)
                    break;
                if (!selector_(*it))
                    continue;

                double phi   = it->phi();
                double invPt = double(it->charge())/it->pt();
                double eta   = it->eta();
                double z0    = it->vz();
                int tt = tt_operator.get_tt(phi, invPt, eta, z0);
                if (tt != the_tt)
                    continue;

                ++n;
            }

            result = (n >= minN_);

        } else {
           edm::LogError("SimpleSkimmer") << "Cannot get the product: " << inputTag_;
        }
    }

    h_denom->Fill(1);
    if (result)  h_num->Fill(1);

    return result;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SimpleSkimmer);
