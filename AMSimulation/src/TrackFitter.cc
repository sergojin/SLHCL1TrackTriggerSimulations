#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/TrackFitter.h"

#include "SLHCL1TrackTriggerSimulations/AMSimulationIO/interface/TTRoadReader.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulationIO/interface/TTTrackReader.h"


namespace {
unsigned getPtSegment(float invPt) {  // for PCA
    return (invPt - PCA_MIN_INVPT) / (PCA_MAX_INVPT - PCA_MIN_INVPT) * PCA_NSEGMENTS;
}

unsigned getHitBits(const std::vector<bool>& stubs_bool) {
    unsigned bitset = 0;

    for (unsigned i=0; i<stubs_bool.size(); ++i) {
        bitset |= (stubs_bool.at(i) << i);
    }

    switch (bitset) {
    case 0b111111:  return 0;
    case 0b111110:  return 1;
    case 0b111101:  return 2;
    case 0b111011:  return 3;
    case 0b110111:  return 4;
    case 0b101111:  return 5;
    case 0b011111:  return 6;
    default      :  return 7;
    }
}

// Comparator
bool sortByPt(const TTTrack2& lhs, const TTTrack2& rhs) {
    return lhs.pt() > rhs.pt();
}

// principal component cut function
bool PrincipalCuts(std::vector<float> princes){
  bool pass=true;
  const float Cuts6[12]={58.,33.,22.,12.,-1.,-1.,9.,3.,3.,3.,-1.,-1.};
  const float Cuts5[10]={33.,22.,12.,-1.,-1.,3.,3.,3.,-1.,-1.};
  if(princes.size()==12) for(unsigned i=0; i<princes.size(); ++i){
    if(Cuts6[i]==-1) continue;
    if(fabs(princes[i])>Cuts6[i]){
      pass=false;
      break;
    }
  }
  else for(unsigned i=0; i<princes.size(); ++i){
    if(Cuts5[i]==-1) continue;
    if(fabs(princes[i])>Cuts5[i]){
      pass=false;
      break;
    }
  }
  return pass;
}

  struct TrackStubConsistency{
    TFile *input;
    std::vector<TH1F*> extrapolators;
    void init(){
      input=new TFile("SLHCL1TrackTriggerSimulations/AMSimulation/data/DeltaSbands_AllExtrapolated_AllWidthsAverage.root","READ");
      for(unsigned i=0; i<6; ++i) extrapolators.push_back((TH1F*)input->Get(TString::Format("Extrapolator%d",i)));   
    }
    float DeltaSchi2(std::vector<std::vector<float> > DeltaS_, float CpT_){
      std::vector<float> chi2;
      float ndof=0;
      float Chi2_DeltaS=0;
      float cuts[6]={12.75,14.75,15.85,11.55,13.45,15.84};
      for(unsigned i=0; i<DeltaS_.size(); ++i) for(unsigned j=0; j<DeltaS_[i].size(); ++j){
      	++ndof;
      	const int pos=extrapolators[i]->FindBin(DeltaS_[i][j]);
	const float result=pow((CpT_-extrapolators[i]->GetBinContent(pos))/extrapolators[i]->GetBinError(pos),2);
      	if(result<cuts[i]) Chi2_DeltaS+=result;
	else Chi2_DeltaS+=999;
      }
      return Chi2_DeltaS/ndof;
    }
    void cleanup(){
      extrapolators.clear();
      input->Close();
    }
  };

}


// _____________________________________________________________________________
// Do track fitting
int TrackFitter::makeTracks(TString src, TString out) {
    if (verbose_)  std::cout << Info() << "Reading " << nEvents_ << " events and fitting tracks." << std::endl;

    // _________________________________________________________________________
    // For reading
    TTRoadReader reader(verbose_);

    if (reader.init(src, prefixRoad_, suffix_)) {
        std::cout << Error() << "Failed to initialize TTRoadReader." << std::endl;
        return 1;
    }

    // _________________________________________________________________________
    // For writing
    TTTrackWriter writer(verbose_);
    if (writer.init(reader.getChain(), out, prefixTrack_, suffix_)) {
        std::cout << Error() << "Failed to initialize TTTrackWriter." << std::endl;
        return 1;
    }

    //test consistency of stubs within track to track fit pT
    TrackStubConsistency StubTest;
    StubTest.init();

    // _________________________________________________________________________
    // Loop over all events

    // Containers
    std::vector<TTTrack2> tracks;
    tracks.reserve(300);

    // Bookkeepers
    long int nRead = 0, nKept = 0;

    //input histogram file


    for (long long ievt=0; ievt<nEvents_; ++ievt) {
        if (reader.loadTree(ievt) < 0)  break;
        reader.getEntry(ievt);

        const unsigned nroads = reader.vr_patternRef->size();
        if (verbose_>1 && ievt%100==0)  std::cout << Debug() << Form("... Processing event: %7lld, fitting: %7ld", ievt, nKept) << std::endl;
        if (verbose_>2)  std::cout << Debug() << "... evt: " << ievt << " # roads: " << nroads << std::endl;

        if (!nroads) {  // skip if no road
            writer.fill(std::vector<TTTrack2>());
            ++nRead;
            continue;
        }

        tracks.clear();
        int fitstatus = 0;


        // _____________________________________________________________________
        // Track fitters taking fit combinations

        // Loop over the roads
        for (unsigned iroad=0; iroad<nroads; ++iroad) {
            if (iroad >= (unsigned) po_.maxRoads)  break;

            const unsigned patternRef = reader.vr_patternRef->at(iroad);
            if (patternRef >= (unsigned) po_.maxPatterns)  continue;

            // Get combinations of stubRefs
            std::vector<std::vector<unsigned> > stubRefs = reader.vr_stubRefs->at(iroad);
	    //clean duplicate stubs
	    if(po_.removeClones) for(unsigned ilayer=0; ilayer<stubRefs.size(); ++ilayer){
	      std::vector<bool> Remove;
	      for(unsigned stubs=0; stubs<stubRefs[ilayer].size(); ++stubs) Remove.push_back(false);
	      for(int stub=0; stub<(int)stubRefs[ilayer].size()-1; ++stub) for(unsigned j=stub+1; j<stubRefs[ilayer].size(); ++j){
		const unsigned stub1=stubRefs[ilayer][stub];
		const unsigned stub2=stubRefs[ilayer][j];
		if(reader.vb_r->at(stub1)==reader.vb_r->at(stub2) && reader.vb_phi->at(stub1)==reader.vb_phi->at(stub2) && reader.vb_z->at(stub1)==reader.vb_z->at(stub2)){
		  if(reader.vb_tpId->at(stub1) <0) Remove[stub]=true;
		  else if(reader.vb_tpId->at(stub2) <0) Remove[j]=true;
		}
	      }//end adjacent stub comparison
	      unsigned offset=0;
	      for(unsigned i=0; i<Remove.size(); ++i){
		if(Remove[i]){
		  stubRefs[ilayer].erase(stubRefs[ilayer].begin()+i-offset); //accounts for previous deletions
		  ++offset;
		}
	      }//end removal procedure
	    }//end layer loop


	    std::vector<std::vector<float> > stubDeltaS; //pass DeltaS information for each stub to the PDDS
            for (unsigned ilayer=0; ilayer<stubRefs.size(); ++ilayer) {
	        std::vector<float> placeholderTemp;
	        stubDeltaS.push_back(placeholderTemp);
		if(po_.PDDS) for(unsigned istub=0; istub<stubRefs[ilayer].size(); ++istub) stubDeltaS[ilayer].push_back(reader.vb_trigBend->at(stubRefs[ilayer][istub]));
		else for(unsigned istub=0; istub<stubRefs[ilayer].size(); ++istub) stubDeltaS[ilayer].push_back(0.); //default DDS is 0 to disable PDDS cleaning
                if (stubRefs.at(ilayer).size() > (unsigned) po_.maxStubs){
                    stubRefs.at(ilayer).resize(po_.maxStubs);
		    stubDeltaS.at(ilayer).resize(po_.maxStubs);
		}
            }
	    
	    //choose either the normal combination building or the 5/6 permutations per 6/6 road in addition and/or pairwise Delta Delta S cleaning (PDDS)
	    std::vector<std::vector<unsigned> > combinations;
	    if (po_.oldCB) combinations = combinationFactory_.combine(stubRefs);
            else if(po_.PDDS) combinations = pairCombinationFactory_.combine(stubRefs, stubDeltaS, po_.FiveOfSix, po_.PDDS);
	    else combinations = combinationBuilderFactory_->combine(stubRefs);

	    // std::cout << "combinations = " << combinations.size() << std::endl;

            for (unsigned icomb=0; icomb<combinations.size(); ++icomb)
                assert(combinations.at(icomb).size() == reader.vr_stubRefs->at(iroad).size());

            if (verbose_>2) {
                std::cout << Debug() << "... ... road: " << iroad << " # combinations: " << combinations.size() << std::endl;
            }

            // Loop over the combinations
            for (unsigned icomb=0; icomb<combinations.size(); ++icomb) {
                if (icomb >= (unsigned) po_.maxCombs)  break;

                // Create and set TTRoadComb
                TTRoadComb acomb;
                acomb.roadRef    = iroad;
                acomb.combRef    = icomb;
                acomb.patternRef = patternRef;
                acomb.ptSegment  = getPtSegment(reader.vr_patternInvPt->at(iroad));
                acomb.stubRefs   = combinations.at(icomb);

                acomb.stubs_r   .clear();
                acomb.stubs_phi .clear();
                acomb.stubs_z   .clear();
                acomb.stubs_bool.clear();

		std::vector<std::vector<float> > DeltaSvector;

                for (unsigned istub=0; istub<acomb.stubRefs.size(); ++istub) {
                    const unsigned stubRef = acomb.stubRefs.at(istub);
                    if (stubRef != CombinationFactory::BAD) {
                        acomb.stubs_r   .push_back(reader.vb_r   ->at(stubRef));
                        acomb.stubs_phi .push_back(reader.vb_phi ->at(stubRef));
                        acomb.stubs_z   .push_back(reader.vb_z   ->at(stubRef));
                        acomb.stubs_bool.push_back(true);
			DeltaSvector.push_back({reader.vb_trigBend->at(stubRef)});
                    } else {
                        acomb.stubs_r   .push_back(0.);
                        acomb.stubs_phi .push_back(0.);
                        acomb.stubs_z   .push_back(0.);
                        acomb.stubs_bool.push_back(false);
			DeltaSvector.push_back({});
                    }
                }

                acomb.hitBits = getHitBits(acomb.stubs_bool);

                if (verbose_>2) {
                    std::cout << Debug() << "... ... ... comb: " << icomb << " " << acomb;
                    std::cout << std::endl;
                }

                // _____________________________________________________________
                // Fit
                TTTrack2 atrack;
                fitstatus = fitter_->fit(acomb, atrack);

                atrack.setTower     (po_.tower);
                atrack.setRoadRef   (acomb.roadRef);
                atrack.setCombRef   (acomb.combRef);
                atrack.setPatternRef(acomb.patternRef);
                atrack.setPtSegment (acomb.ptSegment);
                atrack.setHitBits   (acomb.hitBits);
                atrack.setStubRefs  (acomb.stubRefs);

                if(!po_.CutPrincipals){
		  if (atrack.chi2Red() < po_.maxChi2)  // reduced chi^2 = chi^2 / ndof
 			tracks.push_back(atrack);
		}
		else if(PrincipalCuts(atrack.principals())){
		  float StubChi2=StubTest.DeltaSchi2(DeltaSvector, atrack.invPt());
		  if(StubChi2<po_.maxDeltaSChi2) tracks.push_back(atrack);
		}

                if (verbose_>2)  std::cout << Debug() << "... ... ... track: " << icomb << " status: " << fitstatus << " reduced chi2: " << atrack.chi2Red() << " invPt: " << atrack.invPt() << " phi0: " << atrack.phi0() << " cottheta: " << atrack.cottheta() << " z0: " << atrack.z0() << std::endl;
            }
        }  // loop over the roads

        if (verbose_>2)  std::cout << Debug() << "... evt: " << ievt << " # tracks: " << tracks.size() << std::endl;
        if (verbose_>3) {
            for (unsigned itrack=0; itrack!=tracks.size(); ++itrack) {
                std::cout << "... ... track: " << itrack << " " << tracks.at(itrack) << std::endl;
            }
        }

        // _____________________________________________________________________
        // Find ghosts

        for (unsigned itrack=0; itrack<tracks.size(); ++itrack) {  // all tracks
            for (unsigned jtrack=0; jtrack<itrack; ++jtrack) {  // only non ghost tracks
                if (tracks.at(jtrack).isGhost())  continue;

                bool isGhost = ghostBuster_.isGhostTrack(tracks.at(jtrack).stubRefs(), tracks.at(itrack).stubRefs());
                if (isGhost) {
                    tracks.at(itrack).setAsGhost();
                }
            }
        }

        if (! tracks.empty())
            ++nKept;

        if (tracks.size() > (unsigned) po_.maxTracks)
            tracks.resize(po_.maxTracks);



	// ---------------------------------------------------------------------
	// Classify tracks as duplicates or not (for duplicate removal)
	// In the algorithm tracking particles are sorted by pT
	// And AM tracks are sorted by logic and pT
	// ---------------------------------------------------------------------
	DuplicateRemoval flagDuplicates;
	flagDuplicates.CheckTracks(tracks, po_.rmDuplicate);

	//----------------------------------------------------------------------
	// Identify and flag duplicates by defining a track-parameter space 
	// inside of which anything is considered to be a single track
	//----------------------------------------------------------------------
	ParameterDuplicateRemoval RemoveParameterDuplicates;
	if(po_.rmParDuplicate) RemoveParameterDuplicates.ReduceTracks(tracks);

        std::sort(tracks.begin(), tracks.end(), sortByPt);


        // _____________________________________________________________________        // Track categorization

        if (po_.speedup<1) {
            const unsigned nparts = reader.vp2_primary->size();
            if (verbose_>2)  std::cout << Debug() << "... evt: " << ievt << " # particles: " << nparts << std::endl;

            std::vector<TrackingParticle> trkParts;
            for (unsigned ipart=0; ipart<nparts; ++ipart) {
                bool  primary         = reader.vp2_primary->at(ipart);
                int   simCharge       = reader.vp2_charge->at(ipart);

                if (simCharge!=0 && primary) {
                    float simPt           = reader.vp2_pt->at(ipart);
                    float simEta          = reader.vp2_eta->at(ipart);
                    float simPhi          = reader.vp2_phi->at(ipart);
                    //float simVx           = reader.vp2_vx->at(ipart);
                    //float simVy           = reader.vp2_vy->at(ipart);
                    float simVz           = reader.vp2_vz->at(ipart);
                    int   simCharge       = reader.vp2_charge->at(ipart);
                    int   simPdgId        = reader.vp2_pdgId->at(ipart);

                    float simCotTheta     = std::sinh(simEta);
                    float simChargeOverPt = float(simCharge)/simPt;

                    trkParts.emplace_back(TrackingParticle{  // using POD type constructor
                        (int) ipart,
                        simPdgId,
                        simChargeOverPt,
                        simPhi,
                        simCotTheta,
                        simVz,
                        0.
                    });

                    if (verbose_>3)  std::cout << Debug() << "... ... part: " << ipart << " primary: " << primary << " " << trkParts.back();
                }
            }
            truthAssociator_.associate(trkParts, tracks);
        }

        writer.fill(tracks);
        ++nRead;
    }

    if (nRead == 0) {
        std::cout << Error() << "Failed to read any event." << std::endl;
        return 1;
    }

    if (verbose_)  std::cout << Info() << Form("Read: %7ld, triggered: %7ld", nRead, nKept) << std::endl;

    StubTest.cleanup();

    // _________________________________________________________________________
    // Write histograms

    for (std::map<TString, TH1F *>::const_iterator it=fitter_->histograms_.begin();
         it!=fitter_->histograms_.end(); ++it) {
        if (it->second)  it->second->SetDirectory(gDirectory);
    }

    long long nentries = writer.writeTree();
    assert(nentries == nRead);

    return 0;
}


// _____________________________________________________________________________
// Main driver
int TrackFitter::run() {
    int exitcode = 0;
    Timing(1);

    exitcode = makeTracks(po_.input, po_.output);
    if (exitcode)  return exitcode;
    Timing();

    return exitcode;
}
