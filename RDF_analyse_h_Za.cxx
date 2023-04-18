#include <iostream>
#include <vector>
#include <chrono>
#include <algorithm>
#include <string>
#include <cstdlib>

#include <ROOT/RLogger.hxx>
#include "Math/VectorUtil.h"
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

using Clock = std::chrono::high_resolution_clock;

//https://root-forum.cern.ch/t/rdataframe-count-and-report-re-looping-over-whole-dataframe/46592
void RDF_analyse_h_Za()
{
    auto start_time = Clock::now();
    
    std::vector<std::vector<std::string>> input_filenames =
//    {{"/home/common/Za/NTuples/Ntuple_MC_Za_mA5p0_v4.root"}};
    {{"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root"}};
    SchottDataFrame df(MakeRDF(input_filenames[0],8));
    
   // df.Describe().Print();
   // exit(1);
    auto verbosity = ROOT::Experimental::RLogScopedVerbosity(ROOT::Detail::RDF::RDFLogChannel(), ROOT::Experimental::ELogLevel::kInfo);
    
    auto stable_truth_photons = df.Define("stable_truth_photons", [&](RVec<TruthParticle> truth_particles)
    {
        truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
        [](const TruthParticle& x)
        {
            return (x.mc_pdg_id!=22 || x.mc_status!=1);
            
        }), truth_particles.end());
        
        return truth_particles;
    }, {"truth_particles"});
    
//    std::cout << '\n' << stable_truth_photons.Display<RVec<TruthParticle>,RVec<TruthParticle>> ({"truth_particles", "stable_truth_photons"},100)->AsString() << '\n';
    
    auto stable_truth_photons_pt = stable_truth_photons.Define("stable_truth_photons_pt", [&](const RVec<TruthParticle>& stable_truth_photons)
    {
        RVec<float> pt;
        pt.reserve(stable_truth_photons.size());
        for (auto &i: stable_truth_photons)
        {
            pt.push_back(i.mc_pt/1e3);
        }
        return pt;

    }, {"stable_truth_photons"});
    
//    std::cout << '\n' << stable_truth_photons_pt.Display<RVec<TruthParticle>,RVec<TruthParticle>, RVec<float>> ({"truth_particles", "stable_truth_photons", "stable_truth_photons_pt"},100)->AsString() << '\n';
    
    auto stable_truth_photons_eta = stable_truth_photons.Define("stable_truth_photons_eta", [&](const RVec<TruthParticle>& stable_truth_photons)
    {
        RVec<float> eta;
        eta.reserve(stable_truth_photons.size());
        for (auto &i: stable_truth_photons)
        {
            eta.push_back(i.mc_eta);
        }
        return eta;

    }, {"stable_truth_photons"});
    
//    std::cout << '\n' << stable_truth_photons_eta.Display<RVec<TruthParticle>,RVec<TruthParticle>, RVec<float>> ({"truth_particles", "stable_truth_photons", "stable_truth_photons_eta"},100)->AsString() << '\n';
    
    auto stable_truth_leptons = df.Define("stable_truth_leptons", [&](RVec<TruthParticle> truth_particles)
    {
        truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
        [](const TruthParticle& x)
        {
              return (((abs(x.mc_pdg_id)!=11)
                  && (abs(x.mc_pdg_id)!=12)
                  && (abs(x.mc_pdg_id)!=13)
                  && (abs(x.mc_pdg_id)!=14)
                  && (abs(x.mc_pdg_id)!=15)
                  && (abs(x.mc_pdg_id)!=16)
                  && (abs(x.mc_pdg_id)!=17)
                  && (abs(x.mc_pdg_id)!=18))
                  || (x.mc_status!=1));

        }), truth_particles.end());
        
        return truth_particles;

    }, {"truth_particles"});
    
    auto passed_truth_leptons = stable_truth_leptons.Filter(
    [&](const RVec<TruthParticle>& stable_truth_leptons)
    {
        return
        (
            (stable_truth_leptons.size()==2)
         
            &&
         
            (DeltaR(
                PtEtaPhiEVector(stable_truth_leptons[0].mc_pt,
                                stable_truth_leptons[0].mc_eta,
                                stable_truth_leptons[0].mc_phi,
                                stable_truth_leptons[0].mc_e),

                PtEtaPhiEVector(stable_truth_leptons[1].mc_pt,
                                stable_truth_leptons[1].mc_eta,
                                stable_truth_leptons[1].mc_phi,
                                stable_truth_leptons[1].mc_e)
               ) > 0.01)
         
            &&
         
            (stable_truth_leptons[0].mc_charge == -1*stable_truth_leptons[1].mc_charge)

            &&

            ((((stable_truth_leptons[0].mc_e + stable_truth_leptons[1].mc_e)
            * (stable_truth_leptons[0].mc_e + stable_truth_leptons[1].mc_e))
            -
            ((stable_truth_leptons[0].mc_pt + stable_truth_leptons[1].mc_pt)
            * (stable_truth_leptons[0].mc_pt + stable_truth_leptons[1].mc_pt)))
            >= 6561e6)

            &&

            ((((stable_truth_leptons[0].mc_e + stable_truth_leptons[1].mc_e)
            * (stable_truth_leptons[0].mc_e + stable_truth_leptons[1].mc_e))
            -
            ((stable_truth_leptons[0].mc_pt + stable_truth_leptons[1].mc_pt)
            * (stable_truth_leptons[0].mc_pt + stable_truth_leptons[1].mc_pt)))
            <= 10201e6)

            &&

            (((stable_truth_leptons[0].mc_pt > 20e3) && (stable_truth_leptons[1].mc_pt > 27e3))
            ||
            ((stable_truth_leptons[1].mc_pt > 20e3) && (stable_truth_leptons[0].mc_pt > 27e3)))

        );
        
    } , {"stable_truth_leptons"}, "passed_truth_leptons");
    
    //    auto cutReport = passed_truth_leptons.Report();
    //    cutReport->Print();

    auto passed_truth_leptons_pt = passed_truth_leptons.Define("passed_truth_leptons_pt", [&](const RVec<TruthParticle>& passed_truth_leptons)
    {
        RVec<float> pt;
        pt.reserve(passed_truth_leptons.size());
        for (auto &i: passed_truth_leptons)
        {
            pt.push_back(i.mc_pt/1e3);
        }
        return pt;

    }, {"stable_truth_leptons"});
    
    auto passed_truth_leptons_eta = passed_truth_leptons.Define("passed_truth_leptons_eta", [&](RVec<TruthParticle>& passed_truth_leptons)
    {
        RVec<float> eta;
        eta.reserve(passed_truth_leptons.size());
        for (auto &i: passed_truth_leptons)
        {
            eta.push_back(i.mc_eta);
        }
        return eta;

    }, {"stable_truth_leptons"});
    
    auto selected_photons = df.Define("selected_photons", [&](RVec<Photon> photons)
    {
        photons.erase(std::remove_if(photons.begin(),photons.end(),
        [](Photon& x)
        {
            return ((!x.photon_id)|| (x.photon_pt<=10000) || (x.photon_eta >= 2.37)
                    || (abs(x.photon_eta)>1.37 && abs(x.photon_eta)<1.52));
            
        }), photons.end());
        
        return photons;

    }, {"photons"});
    
//    std::cout << '\n' << selected_photons.Display<RVec<Photon>, RVec<Photon>> ({"photons", "selected_photons"},100)->AsString() << '\n';
    
    auto selected_photons_pt = selected_photons.Define("selected_photons_pt", [&](const RVec<Photon>& selected_photons)
    {
        RVec<float> pt;
        pt.reserve(selected_photons.size());
        for (auto &i: selected_photons)
        {
            pt.push_back(i.photon_pt/1e3);
        }
        return pt;

    }, {"selected_photons"});

    auto selected_photons_eta = selected_photons.Define("selected_photons_eta", [&](const RVec<Photon>& selected_photons)
    {
        RVec<float> eta;
        eta.reserve(selected_photons.size());
        for (auto &i: selected_photons)
        {
            eta.push_back(i.photon_eta);
        }
        return eta;

    }, {"selected_photons"});
    
//    std::cout << '\n' << selected_photons_eta.Display<RVec<Photon>, RVec<Photon>, RVec<float>> ({"photons", "selected_photons", "selected_photons_eta"},100)->AsString() << '\n';
    
    auto passed_selected_photons = selected_photons.Filter(
    [&](const RVec<Photon>& selected_photons)
    {
        return (selected_photons.size()==2);

    } , {"selected_photons"}, "passed_selected_photons");
    
    auto passed_selected_photons_pt = passed_selected_photons.Define("passed_selected_photons_pt", [&](const RVec<Photon>& passed_selected_photons)
    {
        RVec<float> pt;
        pt.reserve(passed_selected_photons.size());
        for (auto &i: passed_selected_photons)
        {
            pt.push_back(i.photon_pt/1e3);
        }
        return pt;

    }, {"selected_photons"});

    auto passed_selected_photons_eta = passed_selected_photons.Define("passed_selected_photons_eta", [&](const RVec<Photon>& passed_selected_photons)
    {
        RVec<float> eta;
        eta.reserve(passed_selected_photons.size());
        for (auto &i: passed_selected_photons)
        {
            eta.push_back(i.photon_eta);
        }
        return eta;

    }, {"selected_photons"});

//    std::cout << '\n' << passed_selected_photons_eta.Display<RVec<Photon>,RVec<Photon>,RVec<float>>({"photons", "selected_photons", "passed_selected_photons_eta"},100)->AsString() << '\n';
    
    auto nEntriesAfterCuts = passed_selected_photons_eta.Count();
    std::cout << "# events for nominal = " << *nEntriesAfterCuts << '\n';
    
    std::vector<ROOT::RDF::RResultPtr<TH1D>> histos;
    histos.reserve(12);
    
    histos = {stable_truth_photons_pt.Histo1D<RVec<float>>({"stable_truth_photons_pt", "histTitle", 20u, 0, 200}, "stable_truth_photons_pt"),
        stable_truth_photons_eta.Histo1D<RVec<float>>({"stable_truth_photons_eta", "histTitle", 40u, -6, 6}, "stable_truth_photons_eta"),
        stable_truth_photons_pt.Histo1D<RVec<float>>({"stable_truth_photons_pt_tight", "histTitle", 100u, 0, 100}, "stable_truth_photons_pt"),
        passed_truth_leptons_pt.Histo1D<RVec<float>>({"passed_truth_leptons_pt", "histTitle", 20u, 0, 200}, "passed_truth_leptons_pt"),
        passed_truth_leptons_eta.Histo1D<RVec<float>>({"passed_truth_leptons_eta", "histTitle", 40u, -6, 6}, "passed_truth_leptons_eta"),
        passed_truth_leptons_pt.Histo1D<RVec<float>>({"passed_truth_leptons_pt_tight", "histTitle", 100u, 0, 100}, "passed_truth_leptons_pt"),
        selected_photons_pt.Histo1D<RVec<float>>({"selected_photons_pt", "histTitle", 20u, 0, 200}, "selected_photons_pt"),
        selected_photons_pt.Histo1D<RVec<float>>({"selected_photons_pt_tight", "histTitle", 100u, 0, 100}, "selected_photons_pt"),
        selected_photons_eta.Histo1D<RVec<float>>({"selected_photons_eta", "histTitle", 40u, -6, 6}, "selected_photons_eta"),
        passed_selected_photons_pt.Histo1D<RVec<float>>({"passed_selected_photons_pt_tight", "histTitle", 100u, 0, 100}, "passed_selected_photons_pt"),
        passed_selected_photons_pt.Histo1D<RVec<float>>({"passed_selected_photons_pt", "histTitle", 20u, 0, 200}, "passed_selected_photons_pt"),
        passed_selected_photons_eta.Histo1D<RVec<float>>({"passed_selected_photons_eta", "histTitle", 40u, -6, 6}, "passed_selected_photons_eta"),
    };
    
    TCanvas *c1;
    for (auto& h: histos)
    {
        c1 = new TCanvas("","",800, 700);
        h->Draw("same");
        c1->SaveAs((h->GetName()+std::string(".png")).c_str());
    }
    
//    system("convert *pdf -quality 100 file.pdf"); //Only if imagemagick is installed.
//    system(R"--(ls *pdf | grep -xv "file.pdf" | parallel rm)--");

    auto end_time = Clock::now();
    std::cout << "Time difference:"
       << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count() << " nanoseconds" << std::endl;
    
}


int main()
{
    RDF_analyse_h_Za();
}

