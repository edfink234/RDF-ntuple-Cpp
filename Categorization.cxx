#include <iostream>
#include <vector>
#include <chrono>
#include <algorithm>
#include <unordered_map>
#include <iomanip>
#include <fstream>

#include <string>
#include <array>
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
#include "THStack.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TLatex.h"
#include "TLegend.h"
#include "Rtypes.h"
#include "ROOT/RDFHelpers.hxx"

#include "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/RDFObjects.h"
#include "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/MakeRDF.h"
#include "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/RDFevent.h"

using namespace ROOT::VecOps; // RVec
using namespace ROOT::Math::VectorUtil; // DeltaR
using namespace ROOT::Math; // PtEtaPhiEVector

using Clock = std::chrono::high_resolution_clock;

constexpr std::array<const char*,35> triggers =
{
    "HLT_e26_lhtight_nod0_ivarloose",
    "HLT_e60_lhmedium_nod0",
    "HLT_e140_lhloose_nod0",
    "HLT_mu26_ivarmedium",
    "HLT_mu50",
    
    "HLT_e26_lhtight_nod0_ivarloose",
    "HLT_e60_lhmedium_nod0",
    "HLT_e140_lhloose_nod0",
    "HLT_mu26_ivarmedium",
    "HLT_mu50",
    
    "HLT_e24_lhtight_nod0_ivarloose",
    "HLT_e26_lhtight_nod0_ivarloose",
    "HLT_e60_lhmedium_nod0",
    "HLT_e60_medium",
    "HLT_e140_lhloose_nod0",
    "HLT_mu26_ivarmedium",
    "HLT_mu50",
    
    "HLT_e24_lhmedium_L1EM20VH",
    "HLT_e60_lhmedium",
    "HLT_e120_lhloose",
    "HLT_mu20_iloose_L1MU15",
    "HLT_mu50",
    
    "HLT_2e17_lhvloose_nod0_L12EM15VHI",
    "HLT_2e17_lhvloose_nod0",
    "HLT_2e24_lhvloose_nod0",
    "HLT_mu22_mu8noL1",
    
    "HLT_2e17_lhvloose_nod0_L12EM15VHI",
    "HLT_2e17_lhvloose_nod0",
    "HLT_2e24_lhvloose_nod0",
    "HLT_mu22_mu8noL1",
    
    "HLT_2e15_lhvloose_nod0_L12EM13VHI",
    "HLT_2e17_lhvloose_nod0",
    "HLT_mu22_mu8noL1",
    
    "HLT_2e12_lhvloose_L12EM10VH",
    "HLT_mu18_mu8noL1",
};
/*
void Table4()
{
    std::vector<std::string> input_filenames = {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root", "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root",
    };
    
    constexpr std::array<float, 2> prefixes = {5.0, 1.0};
    
    auto findParentInChain = [](int targetBarcode, RVec<TruthParticle>& startParticles, RVec<TruthParticle>& truthChain)
    {
        RVec<TruthParticle> truthSelected;
        bool foundParent;
        if (truthChain.size() >= 1)
        {
            TruthParticle tp;
            for (auto& tpe: startParticles)
            {
                tp = tpe;
                while (true)
                {
                    if (tp.mc_parent_barcode == targetBarcode)
                    {
                        truthSelected.push_back(tp);
                        break;
                    }
                    else
                    {
                        foundParent = false;
                        for (auto& tmp: truthChain)
                        {
                            if (tp.mc_parent_barcode == tmp.mc_barcode)
                            {
                                tp = tmp;
                                foundParent = true;
                                break;
                            }
                        }
                        if (foundParent == false)
                        {
                            break;
                        }
                    }
                }
            }
        }
        return truthSelected;
    };
    std::vector<ROOT::RDF::RResultHandle> Nodes;
    std::cout << R"--(\section*{Table 4})--" << '\n';
    std::cout << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
    std::cout << R"--(\begin{tabular}{|c|c|c|c|c|})--" << '\n';
    std::cout << R"--(\hline)--" << '\n';
    std::cout << R"--($m_A$ \, (GeV) & Both photons matched (\%) & 1 photon matched (\%) & At least 1 photon matched (\%) & 0 photons matched (\%) \\ \hline )--" << '\n';
    for (auto& i: input_filenames)
    {
        SchottDataFrame df(MakeRDF({i}, 8));
        
        auto preselection = df.Filter(
        [&](const RVec<std::string>& trigger_passed_triggers, RVec<TruthParticle> truth_particles)
        {
            bool trigger_found = (std::find_first_of(trigger_passed_triggers.begin(), trigger_passed_triggers.end(), triggers.begin(), triggers.end()) != trigger_passed_triggers.end());

            if (!trigger_found)
            {
                return false;
            }
            
            truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
            [](TruthParticle& x)
            {
                return (abs(x.mc_pdg_id) != 11 && abs(x.mc_pdg_id) != 12 && abs(x.mc_pdg_id) != 13 &&
                        abs(x.mc_pdg_id) != 14 && abs(x.mc_pdg_id) != 15 && abs(x.mc_pdg_id) != 16 &&
                        abs(x.mc_pdg_id) != 17 && abs(x.mc_pdg_id) != 18);
                
            }), truth_particles.end());
            
            if (truth_particles.size() != 2)
            {
                return false;
            }
            
            if (DeltaR(truth_particles[0].Vector(), truth_particles[1].Vector()) <= 0.01)
            {
                return false;
            }
            
            if (truth_particles[0].mc_charge*truth_particles[1].mc_charge >= 0)
            {
                return false;
            }
            
            if (!((truth_particles[0].mc_pt > 27e3 && truth_particles[1].mc_pt > 20e3)
                                                   ||
                  (truth_particles[1].mc_pt > 27e3 && truth_particles[0].mc_pt > 20e3)))
            {
                return false;
            }
            
            PtEtaPhiEVector dilepton = truth_particles[0].Vector() + truth_particles[1].Vector();
            
            if ((dilepton.M() < 81e3) || (dilepton.M() > 101e3))
            {
                return false;
            }
            
            if ((truth_particles[0].Vector() + truth_particles[1].Vector()).Pt() <= 10e3)
            {
                return false;
            }
            
            return true;
            
        }, {"trigger_passed_triggers", "truth_particles"});
                
        auto truth_photons_from_axions = preselection.Define("truth_photons",[&](RVec<TruthParticle> truth_particles)
        {
            truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
            [](TruthParticle& x)
            {
                return (abs(x.mc_pdg_id) != 22);
                
            }), truth_particles.end());
            
            return truth_particles;
            
        }, {"truth_particles"})
        .Define("truth_axions", [&](RVec<TruthParticle> truth_particles)
        {
            truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
            [](TruthParticle& x)
            {
                return (abs(x.mc_pdg_id) != 35);
                
            }), truth_particles.end());
            
            return truth_particles;
            
        }, {"truth_particles"}).Define("truth_photons_from_axions",
        [&](RVec<TruthParticle>& truth_photons, RVec<TruthParticle>& truth_axions)
        {
            return findParentInChain(truth_axions[0].mc_barcode, truth_photons, truth_axions);
            
        }, {"truth_photons", "truth_axions"})
        .Filter([&] (RVec<TruthParticle>& truth_photons_from_axions)
        {
            return truth_photons_from_axions.size()==2;
        }, {"truth_photons_from_axions"})
        .Define("photons_pass_cuts",
        [&](RVec<Photon> photons)
        {
            photons.erase(std::remove_if(photons.begin(),photons.end(),
            [](Photon& x)
            {
                return ((abs(x.photon_eta) >= 2.37) || (abs(x.photon_eta) > 1.37 && abs(x.photon_eta) < 1.52));
            }), photons.end());
            
            return photons;
        }, {"photons"}).Filter([](RVec<Photon>& reco_photons_matched)
        {
            if (reco_photons_matched.size() < 2)
            {
                return false;
            }
            auto combs = Combinations(reco_photons_matched, 2);
            size_t length = combs[0].size();
            double delta_r, m, pt, X, best_X, pt1, pt2, chosen_delta_r;

            for (size_t i=0; i<length; i++)
            {
                delta_r = DeltaR(reco_photons_matched[combs[0][i]].Vector(), reco_photons_matched[combs[1][i]].Vector());
                m = (reco_photons_matched[combs[0][i]].Vector() + reco_photons_matched[combs[1][i]].Vector()).M();
                pt = (reco_photons_matched[combs[0][i]].Vector() + reco_photons_matched[combs[1][i]].Vector()).Pt();
                X = delta_r*(pt/(2.0*m));
                if (i==0 || abs(1-X) < abs(1-best_X))
                {
                    best_X = X;
                    pt1 = reco_photons_matched[combs[0][i]].photon_pt;
                    pt2 = reco_photons_matched[combs[1][i]].photon_pt;
                    chosen_delta_r = delta_r;
                }
            }
            if (pt1 > 10e3 && pt2 > 10e3 && best_X > 0.96 && best_X < 1.2 && chosen_delta_r < 1.5)
            {
                return true;
            }
            return false;
            
        }, {"photons_pass_cuts"});
        
        auto reco_photons_matched = truth_photons_from_axions.Define("reco_photons_matched",
        [&](RVec<Photon>& chosen_two, RVec<TruthParticle>& truth_photons_from_axions)
        {
            RVec<Photon> matchedPhotons;
            PtEtaPhiEVector tp1 = truth_photons_from_axions[0].Vector();
            PtEtaPhiEVector tp2 = truth_photons_from_axions[1].Vector();

            for (auto& rp: chosen_two)
            {
                if ((DeltaR(rp.Vector(), tp1) < 0.2 || DeltaR(rp.Vector(), tp2) < 0.2))
                {
                    matchedPhotons.push_back(rp);
                }
            }

            std::sort(matchedPhotons.begin(),matchedPhotons.end(),
            [](Photon& x, Photon& y)
            {
                return x.photon_pt > y.photon_pt;
            });

            return matchedPhotons;

        }, {"photons_pass_cuts", "truth_photons_from_axions"});
        
        auto two_reco_photons_matched = reco_photons_matched.Filter(
        [&](RVec<Photon>& reco_photons_matched)
        {
            return reco_photons_matched.size()==2;
            
        }, {"reco_photons_matched"});
        
        auto one_reco_photons_matched = reco_photons_matched.Filter(
        [&](RVec<Photon>& reco_photons_matched)
        {
            return reco_photons_matched.size()==1;
            
        }, {"reco_photons_matched"});
        
        Nodes.push_back(reco_photons_matched.Count());
        Nodes.push_back(two_reco_photons_matched.Count());
        Nodes.push_back(one_reco_photons_matched.Count());
    }
    
    ROOT::RDF::RunGraphs(Nodes); // running all computation nodes concurrently
    
    int count = 0;
    double total, both, one;
    
    for (int i = 0; i <= 3; i += 3)
    {
        total = static_cast<double>(*Nodes[i].GetResultPtr<ULong64_t>());
        both = static_cast<double>(*Nodes[i+1].GetResultPtr<ULong64_t>());
        one = static_cast<double>(*Nodes[i+2].GetResultPtr<ULong64_t>());
        
        std::cout << prefixes[count++] << " & " << (both / total)*100
        << " & " << (one / total)*100 << " & " << ((one + both) / total)*100
        << " & " << ((total - (one + both)) / total)*100 << R"--(\\ \hline )--" << '\n';
    }
    
    std::cout << R"--(\end{tabular}})--" << '\n';
    
    std::cout << "\n\n\n";
}

void Table5()
{
    std::vector<std::string> input_filenames = {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root", "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root",
    };
    
    constexpr std::array<float, 2> prefixes = {5.0, 1.0};
    
    auto findParentInChain = [](int targetBarcode, RVec<TruthParticle>& startParticles, RVec<TruthParticle>& truthChain)
    {
        RVec<TruthParticle> truthSelected;
        bool foundParent;
        if (truthChain.size() >= 1)
        {
            TruthParticle tp;
            for (auto& tpe: startParticles)
            {
                tp = tpe;
                while (true)
                {
                    if (tp.mc_parent_barcode == targetBarcode)
                    {
                        truthSelected.push_back(tp);
                        break;
                    }
                    else
                    {
                        foundParent = false;
                        for (auto& tmp: truthChain)
                        {
                            if (tp.mc_parent_barcode == tmp.mc_barcode)
                            {
                                tp = tmp;
                                foundParent = true;
                                break;
                            }
                        }
                        if (foundParent == false)
                        {
                            break;
                        }
                    }
                }
            }
        }
        return truthSelected;
    };
    std::vector<ROOT::RDF::RResultHandle> Nodes;
    std::cout << R"--(\section*{Table 5})--" << '\n';
    std::cout << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
    std::cout << R"--(\begin{tabular}{|c|c|c|c|c|})--" << '\n';
    std::cout << R"--(\hline)--" << '\n';
    std::cout << R"--($m_A$ \, (GeV) & Both photons matched (\%) & 1 photon matched (\%) & At least 1 photon matched (\%) & 0 photons matched (\%) \\ \hline )--" << '\n';
    for (auto& i: input_filenames)
    {
        SchottDataFrame df(MakeRDF({i}, 8));
        
        auto preselection = df.Filter(
        [&](const RVec<std::string>& trigger_passed_triggers, RVec<TruthParticle> truth_particles)
        {
            bool trigger_found = (std::find_first_of(trigger_passed_triggers.begin(), trigger_passed_triggers.end(), triggers.begin(), triggers.end()) != trigger_passed_triggers.end());

            if (!trigger_found)
            {
                return false;
            }
            
            truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
            [](TruthParticle& x)
            {
                return (abs(x.mc_pdg_id) != 11 && abs(x.mc_pdg_id) != 12 && abs(x.mc_pdg_id) != 13 &&
                        abs(x.mc_pdg_id) != 14 && abs(x.mc_pdg_id) != 15 && abs(x.mc_pdg_id) != 16 &&
                        abs(x.mc_pdg_id) != 17 && abs(x.mc_pdg_id) != 18);
                
            }), truth_particles.end());
            
            if (truth_particles.size() != 2)
            {
                return false;
            }
            
            if (DeltaR(truth_particles[0].Vector(), truth_particles[1].Vector()) <= 0.01)
            {
                return false;
            }
            
            if (truth_particles[0].mc_charge*truth_particles[1].mc_charge >= 0)
            {
                return false;
            }
            
            if (!((truth_particles[0].mc_pt > 27e3 && truth_particles[1].mc_pt > 20e3)
                                                   ||
                  (truth_particles[1].mc_pt > 27e3 && truth_particles[0].mc_pt > 20e3)))
            {
                return false;
            }
            
            PtEtaPhiEVector dilepton = truth_particles[0].Vector() + truth_particles[1].Vector();
            
            if ((dilepton.M() < 81e3) || (dilepton.M() > 101e3))
            {
                return false;
            }
            
            if ((truth_particles[0].Vector() + truth_particles[1].Vector()).Pt() <= 10e3)
            {
                return false;
            }
            
            return true;
            
        }, {"trigger_passed_triggers", "truth_particles"});
        
    //    std::cout << *preselection.Count() << '\n';
        
        auto truth_photons_from_axions = preselection.Define("truth_photons",[&](RVec<TruthParticle> truth_particles)
        {
            truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
            [](TruthParticle& x)
            {
                return (abs(x.mc_pdg_id) != 22);
                
            }), truth_particles.end());
            
            return truth_particles;
            
        }, {"truth_particles"})
        .Define("truth_axions", [&](RVec<TruthParticle> truth_particles)
        {
            truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
            [](TruthParticle& x)
            {
                return (abs(x.mc_pdg_id) != 35);
                
            }), truth_particles.end());
            
            return truth_particles;
            
        }, {"truth_particles"}).Define("truth_photons_from_axions",
        [&](RVec<TruthParticle>& truth_photons, RVec<TruthParticle>& truth_axions)
        {
            return findParentInChain(truth_axions[0].mc_barcode, truth_photons, truth_axions);
            
        }, {"truth_photons", "truth_axions"})
        .Filter([&] (RVec<TruthParticle>& truth_photons_from_axions)
        {
            return truth_photons_from_axions.size()==2;
        }, {"truth_photons_from_axions"})
        .Define("photons_pass_cuts",
        [&](RVec<Photon> photons)
        {
            photons.erase(std::remove_if(photons.begin(),photons.end(),
            [](Photon& x)
            {
                return ((abs(x.photon_eta) >= 2.37) || (abs(x.photon_eta) > 1.37 && abs(x.photon_eta) < 1.52));
            }), photons.end());
            
            return photons;
        }, {"photons"}).Filter([](RVec<Photon>& reco_photons_test)
        {
            RVec<Photon> reco_photons_matched = reco_photons_test;
            if (reco_photons_matched.size() == 1)
            {
                return reco_photons_matched[0].photon_pt > 20e3;
            }
            else if (reco_photons_matched.empty())
            {
                return false;
            }
            
            auto combs = Combinations(reco_photons_matched, 2);
            size_t length = combs[0].size();
            double delta_r, m, pt, X, best_X, pt1, pt2, chosen_delta_r;

            for (size_t i=0; i<length; i++)
            {
                delta_r = DeltaR(reco_photons_matched[combs[0][i]].Vector(), reco_photons_matched[combs[1][i]].Vector());
                m = (reco_photons_matched[combs[0][i]].Vector() + reco_photons_matched[combs[1][i]].Vector()).M();
                pt = (reco_photons_matched[combs[0][i]].Vector() + reco_photons_matched[combs[1][i]].Vector()).Pt();
                X = delta_r*(pt/(2.0*m));
                if (i==0 || abs(1-X) < abs(1-best_X))
                {
                    best_X = X;
                    pt1 = reco_photons_matched[combs[0][i]].photon_pt;
                    pt2 = reco_photons_matched[combs[1][i]].photon_pt;
                    chosen_delta_r = delta_r;
                }
            }
            if (pt1 > 10e3 && pt2 > 10e3 && best_X > 0.96 && best_X < 1.2 && chosen_delta_r < 1.5)
            {
                return false;
            }
            
            for (auto& p: reco_photons_matched)
            {
                if (p.photon_pt > 20e3)
                {
                    return true;
                }
            }
            return false;
            
        }, {"photons_pass_cuts"});
        
        auto reco_photons_matched = truth_photons_from_axions.Define("reco_photons_matched",
        [&](RVec<Photon>& merged, RVec<TruthParticle>& truth_photons_from_axions)
        {
            PtEtaPhiEVector tp1 = truth_photons_from_axions[0].Vector();
            PtEtaPhiEVector tp2 = truth_photons_from_axions[1].Vector();
            
            for (auto& rp: merged)
            {
                if ((DeltaR(rp.Vector(), tp1) < 0.2 && DeltaR(rp.Vector(), tp2) < 0.2))
                {
                    return RVec<TruthParticle>{truth_photons_from_axions[0], truth_photons_from_axions[1]};
                }
                else if (DeltaR(rp.Vector(), tp1) < 0.2)
                {
                    return RVec<TruthParticle>{truth_photons_from_axions[0]};
                }
                else if (DeltaR(rp.Vector(), tp2) < 0.2)
                {
                    return RVec<TruthParticle>{truth_photons_from_axions[1]};
                }
                else
                {
                    return RVec<TruthParticle>();
                }
            }
            return RVec<TruthParticle>();
            
        }, {"photons_pass_cuts", "truth_photons_from_axions"});
        
        auto two_reco_photons_matched = reco_photons_matched.Filter(
        [&](RVec<TruthParticle>& reco_photons_matched)
        {
            return reco_photons_matched.size()==2;
            
        }, {"reco_photons_matched"});
        
        auto one_reco_photons_matched = reco_photons_matched.Filter(
        [&](RVec<TruthParticle>& reco_photons_matched)
        {
            return reco_photons_matched.size()==1;
            
        }, {"reco_photons_matched"});
        
        Nodes.push_back(reco_photons_matched.Count());
        Nodes.push_back(two_reco_photons_matched.Count());
        Nodes.push_back(one_reco_photons_matched.Count());

    }
    
    ROOT::RDF::RunGraphs(Nodes); // running all computation nodes concurrently
    
    int count = 0;
    double total, both, one;
    
    for (int i = 0; i <= 3; i += 3)
    {
        total = static_cast<double>(*Nodes[i].GetResultPtr<ULong64_t>());
        both = static_cast<double>(*Nodes[i+1].GetResultPtr<ULong64_t>());
        one = static_cast<double>(*Nodes[i+2].GetResultPtr<ULong64_t>());
        
        std::cout << prefixes[count++] << " & " << (both / total)*100
        << " & " << (one / total)*100 << " & " << ((one + both) / total)*100
        << " & " << ((total - (one + both)) / total)*100 << R"--(\\ \hline )--" << '\n';
    }
    
    std::cout << R"--(\end{tabular}})--" << '\n';
    
    std::cout << "\n\n\n";
}*/

void Table14()
{
    std::vector<std::vector<std::string>> input_filenames =
    {
        //Z gamma background
        {
            "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617070._000001.LGNTuple.root",
            "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617064._000001.LGNTuple.root",
            "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617074._000001.LGNTuple.root"
        },
        //Jets
        {
            "/Users/edwardfinkelstein/ATLAS_axion/Jets/Zee_lightJet_0-70.root",
            "/Users/edwardfinkelstein/ATLAS_axion/Jets/Zee_lightJet_70-140.root",
            "/Users/edwardfinkelstein/ATLAS_axion/Jets/Zee_lightJet_140-280.root",
            "/Users/edwardfinkelstein/ATLAS_axion/Jets/Zee_cJet_0-70.root",
            "/Users/edwardfinkelstein/ATLAS_axion/Jets/Zee_cJet_70-140.root",
            "/Users/edwardfinkelstein/ATLAS_axion/Jets/Zee_cJet_140-280.root",
            "/Users/edwardfinkelstein/ATLAS_axion/Jets/Zee_bJet_0-70.root",
            "/Users/edwardfinkelstein/ATLAS_axion/Jets/Zee_bJet_70-140.root",
            "/Users/edwardfinkelstein/ATLAS_axion/Jets/Zee_bJet_140-280.root",
        },
    };
    
    std::vector<ROOT::RDF::RResultHandle> Nodes;
    
    for (auto& i: input_filenames)
    {
        SchottDataFrame df(MakeRDF(i, 8));
        
        auto two_leptons = df.Filter(
        [](RVec<Muon>& muons, RVec<Electron> electrons)
        {
            electrons.erase(std::remove_if(electrons.begin(),electrons.end(),
            [](Electron& ep)
            {
                return (!((ep.electron_pt/1e3 > 20) && (abs(ep.electron_eta) < 2.37) &&
                          (!((1.37 < abs(ep.electron_eta)) && (abs(ep.electron_eta) < 1.52)))
                          && (ep.electron_id_medium == 1)));
                
            }), electrons.end());
            
            return (electrons.size()==2 && muons.empty());
            
        }, {"muons", "electrons"});
        
        auto opp_charge = two_leptons.Define("di_electrons",
        [](RVec<Electron> electrons)
        {
            electrons.erase(std::remove_if(electrons.begin(),electrons.end(),
            [](Electron& ep)
            {
                return (!((ep.electron_pt/1e3 > 20) && (abs(ep.electron_eta) < 2.37) &&
                (!((1.37 < abs(ep.electron_eta)) && (abs(ep.electron_eta) < 1.52)))
                && (ep.electron_id_medium == 1)));

            }), electrons.end());
            
            return electrons;
            
        },{"electrons"})
        .Filter([](RVec<Electron> electrons)
        {
            return (electrons[0].electron_charge*electrons[1].electron_charge < 0);
            
        }, {"di_electrons"});
        
        auto leadingPt = opp_charge.Filter([](RVec<Electron>& electrons)
        {
            return ((electrons[0].electron_pt > 20e3 && electrons[1].electron_pt > 27e3) || (electrons[1].electron_pt > 20e3 && electrons[0].electron_pt > 27e3));
        }, {"di_electrons"});
        
        auto deltaR = leadingPt.Filter([] (RVec<Electron>& electrons)
        {
            return (DeltaR(electrons[0].Vector(), electrons[1].Vector()) > 0.01);
        }, {"di_electrons"});
        
        auto mass = deltaR.Filter([] (RVec<Electron>& electrons)
        {
            auto mass = (electrons[0].Vector() + electrons[1].Vector()).M()/1e3;
            return ((mass >= 81) && (mass <= 101));
        }, {"di_electrons"});
        
        auto ptCut = mass.Filter([] (RVec<Electron>& electrons)
        {
            auto pT = (electrons[0].Vector() + electrons[1].Vector()).Pt()/1e3;
            return pT > 10;
        }, {"di_electrons"});
        
        auto photon_passes_cuts = ptCut.Define("photons_pass_cuts",
        [&](RVec<TruthParticle> truth_particles)
        {
            truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
            [](TruthParticle& x)
            {
                return ((abs(x.mc_pdg_id) != 22) || (abs(x.mc_eta) >= 2.37));

            }), truth_particles.end());
            
            return truth_particles;
            
        }, {"truth_particles"});
        
        auto merged_reco_photons_matched = photon_passes_cuts.Filter(
        [&](RVec<TruthParticle>& reco_photons_test)
        {
            RVec<TruthParticle> reco_photons_matched = reco_photons_test;
            if (reco_photons_matched.size() == 1)
            {
                return reco_photons_matched[0].mc_pt > 20e3;
            }
            else if (reco_photons_matched.empty())
            {
                return false;
            }
            
            auto combs = Combinations(reco_photons_matched, 2);
            size_t length = combs[0].size();
            double delta_r, m, pt, X, best_X, pt1, pt2, chosen_delta_r;

            for (size_t i=0; i<length; i++)
            {
                delta_r = DeltaR(reco_photons_matched[combs[0][i]].Vector(), reco_photons_matched[combs[1][i]].Vector());
                m = (reco_photons_matched[combs[0][i]].Vector() + reco_photons_matched[combs[1][i]].Vector()).M();
                pt = (reco_photons_matched[combs[0][i]].Vector() + reco_photons_matched[combs[1][i]].Vector()).Pt();
                X = delta_r*(pt/(2.0*m));
                if (i==0 || abs(1-X) < abs(1-best_X))
                {
                    best_X = X;
                    pt1 = reco_photons_matched[combs[0][i]].mc_pt;
                    pt2 = reco_photons_matched[combs[1][i]].mc_pt;
                    chosen_delta_r = delta_r;
                }
            }
            if (pt1 > 10e3 && pt2 > 10e3 && best_X > 0.96 && best_X < 1.2 && chosen_delta_r < 1.5)
            {
                return false;
            }
            
            for (auto& p: reco_photons_matched)
            {
                if (p.mc_pt > 20e3)
                {
                    return true;
                }
            }
            return false;
            
        }, {"photons_pass_cuts"})
        .Define("merged_photon",
        [&](RVec<TruthParticle>& reco_photons_matched)
        {
            for (auto& p: reco_photons_matched)
            {
                if (p.mc_pt > 20e3)
                {
                    return p;
                }
            }
            
            return reco_photons_matched[0]; //jic the compiler complains
            
        }, {"photons_pass_cuts"}).Define("merged_pdg_id_origin",
        [](TruthParticle& merged_photon, RVec<TruthParticle>& truth_particles)
        {
            TruthParticle origin = merged_photon;
            int origin_id = 0;
            bool found;

            do
            {
            found = false;
            for (auto& i: truth_particles)
            {
                if (origin.mc_parent_barcode == i.mc_barcode)
                {
                    origin = i;
                    origin_id = i.mc_pdg_id;
                    found = true;
                }
            }
            } while (found);

            return origin_id;
        }, {"merged_photon", "truth_particles"});
        
        Nodes.push_back(merged_reco_photons_matched.Take<int, RVec<int>>("merged_pdg_id_origin"));
    }
    
    ROOT::RDF::RunGraphs(Nodes); // running all computation nodes concurrently
    
    std::cout << R"--(\section*{Table 14}\flushleft)--" << '\n';
    std::cout << R"--(\hspace{-3cm}\scalebox{0.65}{\setlength\extrarowheight{2pt}\renewcommand{\arraystretch}{1.5})--" << '\n';
    std::cout << R"--(\begin{tabular}{|c|c|c|c|})--" << '\n';
    std::cout << R"--(\hline)--" << '\n';
    std::cout << R"--(\multicolumn{4}{|c|}{\parbox{\linewidth}{\centering Merged Category \\ Background photon truth origin fractions}}\\[5 pt] \hline)--" << '\n';
    std::cout << R"--( {Pdg Id} & \multicolumn{1}{|c|}{\parbox{7.5cm}{\centering $Z\gamma$ \\ (\%)}} & \multicolumn{1}{|c|}{\parbox{7.5cm}{\centering $Z$ jets \\ (\%)}} & \multicolumn{1}{|c|}{\parbox{7.5cm}{\centering Total bkg \\ (\%)}} \\[5 pt]
        \hline)--" << '\n';
    
    std::unordered_map<int,int> merged_id_freqs_Z_gamma, merged_id_freqs_Z_jets, merged_id_freqs;
    int total_Z_gamma = 0, total_Z_jets = 0, total = 0;
    
    for (auto& i: *Nodes[0].GetResultPtr<RVec<int>>())
    {
        merged_id_freqs_Z_gamma[i]++;
        merged_id_freqs[i]++;
        total_Z_gamma++;
    }
    
    for (auto& i: *Nodes[1].GetResultPtr<RVec<int>>())
    {
        merged_id_freqs_Z_jets[i]++;
        merged_id_freqs[i]++;
        total_Z_jets++;
    }
    
    total = total_Z_gamma + total_Z_jets;
    
    for (auto& i: merged_id_freqs)
    {
        std::cout << i.first << " & " << std::setprecision(2) << std::fixed
        << 100*(merged_id_freqs_Z_gamma[i.first] / static_cast<double>(total_Z_gamma))
        << " & " << std::setprecision(2) << std::fixed
        << 100*(merged_id_freqs_Z_jets[i.first] / static_cast<double>(total_Z_jets))
        << " & " << std::setprecision(2) << std::fixed
        << 100*(merged_id_freqs[i.first] / static_cast<double>(total))
        
        << R"--( \\ \hline)--" << '\n';
    }
    
    std::cout << R"--(\end{tabular}})--" << "\n\n\n";
}
/*
void Table15()
{
    std::vector<std::string> input_filenames =
    {
        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617070._000001.LGNTuple.root", "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617064._000001.LGNTuple.root",
        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617074._000001.LGNTuple.root"
    };
    
    SchottDataFrame df(MakeRDF(input_filenames, 8));
    
    auto two_leptons = df.Filter(
    [](RVec<Muon>& muons, RVec<Electron> electrons)
    {
        electrons.erase(std::remove_if(electrons.begin(),electrons.end(),
        [](Electron& ep)
        {
            return (!((ep.electron_pt/1e3 > 20) && (abs(ep.electron_eta) < 2.37) &&
                      (!((1.37 < abs(ep.electron_eta)) && (abs(ep.electron_eta) < 1.52)))
                      && (ep.electron_id_medium == 1)));
            
        }), electrons.end());
        
        return (electrons.size()==2 && muons.empty());
        
    }, {"muons", "electrons"});
    
    auto opp_charge = two_leptons.Define("di_electrons",
    [](RVec<Electron> electrons)
    {
        electrons.erase(std::remove_if(electrons.begin(),electrons.end(),
        [](Electron& ep)
        {
            return (!((ep.electron_pt/1e3 > 20) && (abs(ep.electron_eta) < 2.37) &&
            (!((1.37 < abs(ep.electron_eta)) && (abs(ep.electron_eta) < 1.52)))
            && (ep.electron_id_medium == 1)));

        }), electrons.end());
        
        return electrons;
        
    },{"electrons"})
    .Filter([](RVec<Electron> electrons)
    {
        return (electrons[0].electron_charge*electrons[1].electron_charge < 0);
        
    }, {"di_electrons"});
    
    auto leadingPt = opp_charge.Filter([](RVec<Electron>& electrons)
    {
        return ((electrons[0].electron_pt > 20e3 && electrons[1].electron_pt > 27e3) || (electrons[1].electron_pt > 20e3 && electrons[0].electron_pt > 27e3));
    }, {"di_electrons"});
    
    auto deltaR = leadingPt.Filter([] (RVec<Electron>& electrons)
    {
        return (DeltaR(electrons[0].Vector(), electrons[1].Vector()) > 0.01);
    }, {"di_electrons"});
    
    auto mass = deltaR.Filter([] (RVec<Electron>& electrons)
    {
        auto mass = (electrons[0].Vector() + electrons[1].Vector()).M()/1e3;
        return ((mass >= 81) && (mass <= 101));
    }, {"di_electrons"});
    
    auto ptCut = mass.Filter([] (RVec<Electron>& electrons)
    {
        auto pT = (electrons[0].Vector() + electrons[1].Vector()).Pt()/1e3;
        return pT > 10;
    }, {"di_electrons"});
    
    auto photon_passes_cuts = ptCut.Define("photons_pass_cuts",
    [&](RVec<Photon> photons)
    {
        photons.erase(std::remove_if(photons.begin(),photons.end(),
        [](Photon& x)
        {
          return ((abs(x.photon_eta) >= 2.37) || (x.photon_pt <= 10e3) || (abs(x.photon_eta) > 1.37 && abs(x.photon_eta) < 1.52) || (!x.photon_id_loose));

        }), photons.end());

        return photons;
    }, {"photons"});
    
    auto merged_reco_photons_matched = photon_passes_cuts.Filter(
    [&](RVec<Photon>& reco_photons_test)
    {
        RVec<Photon> reco_photons_matched = reco_photons_test;
        if (reco_photons_matched.size() == 1)
        {
            return reco_photons_matched[0].photon_pt > 20e3;
        }
        else if (reco_photons_matched.empty())
        {
            return false;
        }
        
        auto combs = Combinations(reco_photons_matched, 2);
        size_t length = combs[0].size();
        double delta_r, m, pt, X, best_X, pt1, pt2, chosen_delta_r;

        for (size_t i=0; i<length; i++)
        {
            delta_r = DeltaR(reco_photons_matched[combs[0][i]].Vector(), reco_photons_matched[combs[1][i]].Vector());
            m = (reco_photons_matched[combs[0][i]].Vector() + reco_photons_matched[combs[1][i]].Vector()).M();
            pt = (reco_photons_matched[combs[0][i]].Vector() + reco_photons_matched[combs[1][i]].Vector()).Pt();
            X = delta_r*(pt/(2.0*m));
            if (i==0 || abs(1-X) < abs(1-best_X))
            {
                best_X = X;
                pt1 = reco_photons_matched[combs[0][i]].photon_pt;
                pt2 = reco_photons_matched[combs[1][i]].photon_pt;
                chosen_delta_r = delta_r;
            }
        }
        if (pt1 > 10e3 && pt2 > 10e3 && best_X > 0.96 && best_X < 1.2 && chosen_delta_r < 1.5)
        {
            return false;
        }
        
        for (auto& p: reco_photons_matched)
        {
            if (p.photon_pt > 20e3)
            {
                return true;
            }
        }
        return false;
        
    }, {"photons_pass_cuts"}).Define("merged_photon",
    [&](RVec<Photon>& reco_photons_matched)
    {
        for (auto& p: reco_photons_matched)
        {
            if (p.photon_pt > 20e3)
            {
              return p;
            }
        }

        return reco_photons_matched[0]; //jic the compiler complains

    }, {"photons_pass_cuts"}).Define("reconstructed_mass",
    [&](RVec<Electron>& di_electrons, Photon& merged_photon)
    {
        auto four_momentum = di_electrons[0].Vector() + di_electrons[1].Vector();
        
        return (four_momentum + merged_photon.Vector()).M()/1e3;
        
    }, {"di_electrons", "merged_photon"});
    
    auto pSB = merged_reco_photons_matched.Filter(
    [](double reconstructed_mass)
    {
        return (reconstructed_mass < 110) || (reconstructed_mass > 130);
    }, {"reconstructed_mass"});
    
    auto pSR = merged_reco_photons_matched.Filter(
    [](double reconstructed_mass)
    {
        return (reconstructed_mass >= 110) && (reconstructed_mass <= 130);
    }, {"reconstructed_mass"});
    
    auto SB = pSB.Filter(
    [](RVec<float>& Eratio)
    {
        return (!Any(Eratio <= 0.8));
    }, {"photon_shower_shape_e_ratio"});
    
    auto SR = pSR.Filter(
    [](RVec<float>& Eratio)
    {
        return (!Any(Eratio <= 0.8));
    }, {"photon_shower_shape_e_ratio"});
    
    auto pSB_count = pSB.Count();
    auto pSR_count = pSR.Count();
    auto SB_count = SB.Count();
    auto SR_count = SR.Count();
    
    std::cout << R"--(\section*{Table 15})--" << '\n';
    std::cout << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
    std::cout << R"--(\begin{tabular}{|c|c|})--" << '\n';
    std::cout << R"--(\hline)--" << '\n';
    std::cout << R"--({} & $Z\gamma$ \\ \hline )--" << '\n';
    std::cout << "pSB & " << *pSB_count << R"--(\\ \hline)--" << '\n';
    std::cout << "pSR & " << *pSR_count << R"--(\\ \hline)--" << '\n';
    std::cout << "SB & " << *SB_count << R"--(\\ \hline)--" << '\n';
    std::cout << "SR & " << *SR_count << R"--(\\ \hline)--" << '\n';
    std::cout << R"--(\end{tabular}})--" << "\n\n\n";
}

void Fig19()
{
    std::vector<std::string> input_filenames = {
        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root",
        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root"
    };
    
    std::vector<float> massPoints = {5, 1};
    std::ofstream out("Fig19.txt");
    std::vector<ROOT::RDF::RResultHandle> Nodes;
    
    auto findParentInChain = [](int targetBarcode, RVec<TruthParticle>& startParticles, RVec<TruthParticle>& truthChain)
    {
        RVec<TruthParticle> truthSelected;
        bool foundParent;
        if (truthChain.size() >= 1)
        {
            TruthParticle tp;
            for (auto& tpe: startParticles)
            {
                tp = tpe;
                while (true)
                {
                    if (tp.mc_parent_barcode == targetBarcode)
                    {
                        truthSelected.push_back(tp);
                        break;
                    }
                    else
                    {
                        foundParent = false;
                        for (auto& tmp: truthChain)
                        {
                            if (tp.mc_parent_barcode == tmp.mc_barcode)
                            {
                                tp = tmp;
                                foundParent = true;
                                break;
                            }
                        }
                        if (foundParent == false)
                        {
                            break;
                        }
                    }
                }
            }
        }
        return truthSelected;
    };
    
    std::vector<const char*> prefixes = {"no cat", "merged", "resolved"};
    std::vector<EColor> colors = {kCyan, kOrange, kBlue};
    
    for (auto& file: input_filenames)
    {
        SchottDataFrame df(MakeRDF({file}, 8));
        auto two_leptons = df.Filter(
        [](RVec<Muon>& muons, RVec<Electron> electrons)
        {
            electrons.erase(std::remove_if(electrons.begin(),electrons.end(),
            [](Electron& ep)
            {
                return (!((ep.electron_pt/1e3 > 20) && (abs(ep.electron_eta) < 2.37) &&
                          (!((1.37 < abs(ep.electron_eta)) && (abs(ep.electron_eta) < 1.52)))
                          && (ep.electron_id_medium == 1)));
                
            }), electrons.end());
            
            return (electrons.size()==2 && muons.empty());
            
        }, {"muons", "electrons"});
        
        auto opp_charge = two_leptons.Define("di_electrons",
        [](RVec<Electron> electrons)
        {
            electrons.erase(std::remove_if(electrons.begin(),electrons.end(),
            [](Electron& ep)
            {
                return (!((ep.electron_pt/1e3 > 20) && (abs(ep.electron_eta) < 2.37) &&
                (!((1.37 < abs(ep.electron_eta)) && (abs(ep.electron_eta) < 1.52)))
                && (ep.electron_id_medium == 1)));

            }), electrons.end());
            
            return electrons;
            
        },{"electrons"})
        .Filter([](RVec<Electron> electrons)
        {
            return (electrons[0].electron_charge*electrons[1].electron_charge < 0);
            
        }, {"di_electrons"});
        
        auto leadingPt = opp_charge.Filter([](RVec<Electron>& electrons)
        {
            return ((electrons[0].electron_pt > 20e3 && electrons[1].electron_pt > 27e3) || (electrons[1].electron_pt > 20e3 && electrons[0].electron_pt > 27e3));
        }, {"di_electrons"});
        
        auto deltaR = leadingPt.Filter([] (RVec<Electron>& electrons)
        {
            return (DeltaR(electrons[0].Vector(), electrons[1].Vector()) > 0.01);
        }, {"di_electrons"});
        
        auto mass = deltaR.Filter([] (RVec<Electron>& electrons)
        {
            auto mass = (electrons[0].Vector() + electrons[1].Vector()).M()/1e3;
            return ((mass >= 81) && (mass <= 101));
        }, {"di_electrons"});
        
        auto ptCut = mass.Filter([] (RVec<Electron>& electrons)
        {
            auto pT = (electrons[0].Vector() + electrons[1].Vector()).Pt()/1e3;
            return pT > 10;
        }, {"di_electrons"});
        
        auto photon_passes_cuts = ptCut.Define("photons_pass_cuts",
        [&](RVec<Photon> photons)
        {
            photons.erase(std::remove_if(photons.begin(),photons.end(),
            [](Photon& x)
            {
                return ((abs(x.photon_eta) >= 2.37) || (abs(x.photon_eta) > 1.37 && abs(x.photon_eta) < 1.52));

            }), photons.end());

            return photons;
        }, {"photons"});
        
        auto resolved = photon_passes_cuts.Define("chosen_two",
        [](RVec<Photon>& photons_pass_cuts)
        {
            RVec<Photon> x;
            if (photons_pass_cuts.size() < 2)
            {
                return x;
            }
            auto combs = Combinations(photons_pass_cuts, 2);
            size_t length = combs[0].size();
            double delta_r, m, pt, X, best_X, pt1, pt2, chosen_delta_r;

            for (size_t i=0; i<length; i++)
            {
                delta_r = DeltaR(photons_pass_cuts[combs[0][i]].Vector(), photons_pass_cuts[combs[1][i]].Vector());
                m = (photons_pass_cuts[combs[0][i]].Vector() + photons_pass_cuts[combs[1][i]].Vector()).M();
                pt = (photons_pass_cuts[combs[0][i]].Vector() + photons_pass_cuts[combs[1][i]].Vector()).Pt();
                X = delta_r*(pt/(2.0*m));
                if (i==0 || abs(1-X) < abs(1-best_X))
                {
                    best_X = X;
                    pt1 = photons_pass_cuts[combs[0][i]].photon_pt;
                    pt2 = photons_pass_cuts[combs[1][i]].photon_pt;
                    chosen_delta_r = delta_r;
                    x = {photons_pass_cuts[combs[0][i]], photons_pass_cuts[combs[1][i]]};
                }
            }
            if (pt1 > 10e3 && pt2 > 10e3 && best_X > 0.96 && best_X < 1.2 && chosen_delta_r < 1.5)
            {
                return x;
            }
            x.clear();
            return x;
        }, {"photons_pass_cuts"}).Filter(
        [&](RVec<Photon>& reco_photons_matched)
        {
            return (reco_photons_matched.size()==2);
        }, {"chosen_two"});
        
        auto merged = photon_passes_cuts.Filter(
        [&](RVec<Photon>& reco_photons_test)
        {
            RVec<Photon> reco_photons_matched = reco_photons_test;
            if (reco_photons_matched.size() == 1)
            {
                return reco_photons_matched[0].photon_pt > 20e3;
            }
            else if (reco_photons_matched.empty())
            {
                return false;
            }
            
            auto combs = Combinations(reco_photons_matched, 2);
            size_t length = combs[0].size();
            double delta_r, m, pt, X, best_X, pt1, pt2, chosen_delta_r;

            for (size_t i=0; i<length; i++)
            {
                delta_r = DeltaR(reco_photons_matched[combs[0][i]].Vector(), reco_photons_matched[combs[1][i]].Vector());
                m = (reco_photons_matched[combs[0][i]].Vector() + reco_photons_matched[combs[1][i]].Vector()).M();
                pt = (reco_photons_matched[combs[0][i]].Vector() + reco_photons_matched[combs[1][i]].Vector()).Pt();
                X = delta_r*(pt/(2.0*m));
                if (i==0 || abs(1-X) < abs(1-best_X))
                {
                    best_X = X;
                    pt1 = reco_photons_matched[combs[0][i]].photon_pt;
                    pt2 = reco_photons_matched[combs[1][i]].photon_pt;
                    chosen_delta_r = delta_r;
                }
            }
            if (pt1 > 10e3 && pt2 > 10e3 && best_X > 0.96 && best_X < 1.2 && chosen_delta_r < 1.5)
            {
                return false;
            }
            
            for (auto& p: reco_photons_matched)
            {
                if (p.photon_pt > 20e3)
                {
                    return true;
                }
            }
            return false;
            
        }, {"photons_pass_cuts"});
        
        Nodes.push_back(photon_passes_cuts.Count());
        Nodes.push_back(merged.Count());
        Nodes.push_back(resolved.Count());
    }
    
    ROOT::RDF::RunGraphs(Nodes); // running all computation nodes concurrently
    
    double total, numMerged, numResolved, noCat;
    
    for (int i = 0; i <= 3; i += 3)
    {
        total = static_cast<double>(*Nodes[i].GetResultPtr<ULong64_t>());
        numMerged = static_cast<double>(*Nodes[i+1].GetResultPtr<ULong64_t>());
        numResolved = static_cast<double>(*Nodes[i+2].GetResultPtr<ULong64_t>());
        noCat = total - (numMerged+numResolved);
        
        out << noCat / total << '\n' << numMerged / total  << '\n'
            << numResolved / total << '\n';
    }
    
    
    out.close();
    
    std::ifstream in("Fig19.txt");
    RVec<TH1D*> hists = {new TH1D("no cat", "no cat", 10, 0, 10), new TH1D("merged", "merged", 10, 0, 10), new TH1D("resolved", "resolved", 10, 0, 10)};
    auto hs = new THStack("hs1","");
    TCanvas* c1 = new TCanvas();
    TLegend* legend = new TLegend(0.65, 0.4, 0.85, 0.6);
    double val;
    for (int i = 0; in >> val; i++)
    {
        hists[i%3]->SetBinContent(hists[i%3]->FindFixBin(massPoints[i/3]), val);
        hists[i%3]->SetFillColor(colors[i%3]);
    }
    for (auto& hist: hists)
    {
        legend->AddEntry(hist, hist->GetTitle(), "f");
        hs->Add(hist);
    }
    hs->Draw("nostackb");
    hs->SetTitle(";m_{A} [GeV];Fraction of Events");
    hs->GetYaxis()->CenterTitle(true);
    hs->GetXaxis()->SetTitleOffset(0.85);
    hs->GetYaxis()->SetTitleFont(42);
    hs->GetXaxis()->SetTitleFont(42);
    hs->GetXaxis()->SetTitleSize(0.04);
    hs->GetYaxis()->SetTitleSize(0.04);
    hs->GetYaxis()->SetTitleOffset(1.1);
    hs->GetXaxis()->SetRangeUser(0,10);
    hs->SetMaximum(0.7);
    hs->GetXaxis()->SetTickLength(0);
    hs->GetXaxis()->SetLabelOffset(999);
    
    gStyle->SetOptStat(0);
    TLatex Tl;
    Tl.SetTextSize(0.04);
    Tl.DrawLatexNDC(0.5, 0.83, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.525, 0.05, "5.0");
    Tl.DrawLatexNDC(0.2, 0.05, "1.0");
    Tl.DrawLatexNDC(0.5, 0.73,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("Fig19.png");
}
*/
void Categorization()
{
    auto start_time = Clock::now();
//    Table4();
//    Table5();
    Table14();
//    Table15();
//    Fig19();
    auto end_time = Clock::now();
    std::cout << "Time difference: "
       << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count()/1e9 << " seconds" << std::endl;
}

int main()
{
    Categorization();
}
