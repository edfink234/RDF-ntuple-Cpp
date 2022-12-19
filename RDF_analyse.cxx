#include <iostream>
#include <vector>
#include <chrono>
#include <algorithm>
#include <string>
#include <cstdlib>
#include <memory>

#include <ROOT/RLogger.hxx>
#include "Math/VectorUtil.h"
#include "ROOT/RDF/RResultMap.hxx"
#include "ROOT/RDataFrame.hxx"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH1D.h"

#include "RDFObjects.h"
#include "MakeRDF.h"
#include "RDFevent.h"

using namespace ROOT::VecOps; // RVec
using namespace ROOT::Math::VectorUtil; // DeltaR
using namespace ROOT::Math; // PtEtaPhiEVector

using ROOT::RDF::Experimental::VariationsFor;
//using ROOT::RDF::Experimental::RResultMap;
//using ct = typename ROOT::RDF::RDisplay;
//using ct = typename ROOT::RDF::Experimental::RResultMap<RDisplay>;

using Clock = std::chrono::high_resolution_clock;

//https://root-forum.cern.ch/t/rdataframe-count-and-report-re-looping-over-whole-dataframe/46592
void RDF_analyse()
{
    auto start_time = Clock::now();
    
    Event::systematics = {"EG_RESOLUTION_ALL__1down"};
    std::vector<std::string> input_filenames = {"/Users/edwardfinkelstein/ATLAS_axion/user.kschmied.28655874._000025.LGNTuple.root"};
    
    SchottDataFrame df(MakeRDF(input_filenames));
    
//    std::cout << '\n' << df.Display<RVec<Photon>> ({"photons"},100)->AsString() << '\n';
    

    
    
    auto disp = df.Display<RVec<Photon>> ({"photons"},100);
//    auto Varied = VariationsFor(disp);
    
//    df.Describe().Print();
    
    auto end_time = Clock::now();
    std::cout << "\nTime difference:"
       << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count() << " nanoseconds" << std::endl;
    
    
}


int main()
{
    RDF_analyse();
}


