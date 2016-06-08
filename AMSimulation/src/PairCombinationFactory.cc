#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/PairCombinationFactory.h"
using namespace slhcl1tt;

#include <iostream>

float PairCombinationFactory::Coarsener(float DeltaS, int precision, int layer){ //maximum precision is 4, minimum 2
  float Coarsed=0.;
  if(layer<3) Coarsed = DeltaS-((int)(2*DeltaS)%(int)pow(2,3-(precision-1)))/2.;
  else        Coarsed = DeltaS-((int)(2*DeltaS)%(int)pow(2,4-precision))/2.;
  return Coarsed;
}

std::vector<PairAssignment> PairCombinationFactory::Stage1Pair(std::vector<std::pair<unsigned, float> > first, std::vector<std::pair<unsigned, float> > second, float cut)
{  
  //initialize pair assignment struct, loop and filter pairs
  std::vector<PairAssignment> PairedOnes;
  for(unsigned i=0; i<first.size(); ++i){
    for(unsigned j=0; j<second.size(); ++j){
      bool Survivor=true;
      if(first[i].first==PairCombinationFactory::BAD && second[j].first==PairCombinationFactory::BAD) Survivor=false; //catch double dummy cases
      else if(fabs(first[i].second-second[j].second)>cut && !(first[i].first==PairCombinationFactory::BAD || second[j].first==PairCombinationFactory::BAD) ) Survivor=false; //pairwise DDS cut, spares any pair with only a single dummy
      PairedOnes.push_back(PairAssignment(first[i].first,second[j].first,Survivor));
    }
  }
  //give back assignment map
  return PairedOnes;
}

//combination(stubAddresses,pass)
std::vector<std::vector<unsigned> > PairCombinationFactory::CombinationBuilder(std::vector<PairAssignment> Layer10, std::vector<PairAssignment> Layer12, std::vector<PairAssignment> Layer32, std::vector<PairAssignment> Layer34, std::vector<PairAssignment> Layer54, unsigned l0, unsigned l2, unsigned l4)
{
  //initialize return object
  std::vector<std::vector<unsigned> >  Combinations;
  //loop over all permutations
  for(unsigned i=0; i<Layer10.size(); ++i){
    for(unsigned j=0; j<Layer32.size(); ++j){
      for(unsigned k=0; k<Layer54.size(); ++k){
	if(!(Layer10[i].SurvivingPair * Layer32[j].SurvivingPair * Layer54[k].SurvivingPair)) continue; //3-logic part
	//5-logic stub position lookup
	int pos12=i/l0*l2+j%l2;
	int pos34=j/l2*l4+k%l4;
	if(!(Layer12[pos12].SurvivingPair * Layer34[pos34].SurvivingPair)) continue; //5-logic extension
	std::vector<unsigned> layers={Layer10[i].Stub2Address, Layer10[i].Stub1Address, Layer32[j].Stub2Address, Layer32[j].Stub1Address, Layer54[k].Stub2Address, Layer54[k].Stub1Address};

	//catch >1 dummy cases, always put dummy in the first position per layer!
	unsigned dummyCounter=0;
	for(unsigned l=0; l<layers.size(); ++l) if(layers[l]==PairCombinationFactory::BAD) ++dummyCounter;
	if(dummyCounter>1) continue;

	//build combination object
	Combinations.push_back(layers);
      }
    }
  }
  return Combinations;
}

std::vector<std::vector<unsigned> > PairCombinationFactory::combine(const std::vector<std::vector<unsigned> >& groups, std::vector<std::vector<float> > DeltaS, bool FiveOfSix) {
        std::vector<std::vector<unsigned> > combinations;
        std::vector<std::vector<unsigned> > CombinationBase=groups;

        unsigned ngroups = groups.size(); //number of layers
	std::vector<std::vector<std::pair<unsigned, float> > > Layer; //0 is innermost, 5 is outermost

	//check for empty layers, add 1 dummy to nonempty layers
	for(unsigned i=0; i<ngroups; ++i){
	  if(!CombinationBase[i].size()){
	    CombinationBase[i].push_back(PairCombinationFactory::BAD);
	    DeltaS[i].push_back(999.);
	  }
	  else if(FiveOfSix){
	    CombinationBase[i].insert(CombinationBase[i].begin(),PairCombinationFactory::BAD); //only insert dummies in full layers, if 5/6 permutations are desired
	    DeltaS[i].insert(DeltaS[i].begin(),999.);
	  }
	}

	if(ngroups!=6) return combinations; //catch less than 6 layers present
	for(unsigned layer=0; layer<6; ++layer){ //compose all the layer information
	  std::vector<std::pair<unsigned, float> > OneLayer;
	  Layer.push_back(OneLayer);
	  for(unsigned j=0; j<CombinationBase[layer].size(); ++j) Layer[layer].push_back(std::make_pair(CombinationBase[layer][j], Coarsener(DeltaS[layer][j],4,layer))); //pushes back stub iterator and coarsened DeltaS value
	}
	
	//build layer pairs & pass cut values
	std::vector<PairAssignment> Layer01=Stage1Pair(Layer[1],Layer[0],1.5);
	std::vector<PairAssignment> Layer12=Stage1Pair(Layer[1],Layer[2],2.5);
	std::vector<PairAssignment> Layer23=Stage1Pair(Layer[3],Layer[2],3.0);
	std::vector<PairAssignment> Layer34=Stage1Pair(Layer[3],Layer[4],2.0);
	std::vector<PairAssignment> Layer45=Stage1Pair(Layer[5],Layer[4],2.0);

	//build the actual combinations
	combinations = CombinationBuilder(Layer01,Layer12,Layer23,Layer34,Layer45,CombinationBase[0].size(),CombinationBase[2].size(),CombinationBase[4].size());
	//for(unsigned c=0; c<combinations.size(); ++c) std::cout<<combinations[c][0]<<","<<combinations[c][1]<<","<<combinations[c][2]<<","<<combinations[c][3]<<","<<combinations[c][4]<<","<<combinations[c][5]<<std::endl;  //combination verbose
        return combinations;
}

// _____________________________________________________________________________
void PairCombinationFactory::print() {
    std::cout << std::endl;
}


