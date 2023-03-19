#include <iostream>
#include <vector>
#include <chrono>
#include <algorithm>
#include <string>
#include <sstream>
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
#include "TStyle.h"
#include "TPaveText.h"
#include "TLatex.h"
#include "TLegend.h"
#include "Rtypes.h"
#include "ROOT/RDFHelpers.hxx"

#include "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/RDFObjects.h"
#include "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/MakeRDF.h"
#include "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/RDFevent.h"

using namespace ROOT::VecOps; // RVec, Combinations
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

void Table11MuonElectron()
{
    std::vector<std::vector<std::string>> input_filenames = {
        //{"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/muon17.3.root", "muon17.3part2.root"}
        {"Electrontest.root"}
    };
    
    constexpr std::array<const char*,1> Samples = {R"--(Signal $m_{\text{A}}$ = 5 GeV)--"};
    
    std::vector<ROOT::RDF::RResultHandle> Nodes;

    std::ostringstream os, ss;
    os << R"--(\section*{Table 11})--" << '\n';
    os << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
    os << R"--(\begin{tabular}{|c|c|c|c|c|c|c|})--" << '\n';
    os << R"--(\hline)--" << '\n';
    os << R"--(Sample & Total Events & pass preselection (PS) & photon $p_T$ + $\Delta R_{\gamma\gamma}$ cut & $X$ window & SR & SR-ID
           \\ \hline )--" << '\n';
    
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
    
    bool comb = true;
    std::string electrons_or_muons;

    for (auto& file: input_filenames)
    {
        SchottDataFrame df(MakeRDF(file, 8));
        
        auto trigger_selection = df.Filter(
        [](const RVec<std::string>& trigger_passed_triggers)
        {
            bool trigger_found = (std::find_first_of(trigger_passed_triggers.begin(), trigger_passed_triggers.end(), triggers.begin(), triggers.end()) != trigger_passed_triggers.end());

            if (!trigger_found)
            {
             return false;
            }
            return true;

        }, {"trigger_passed_triggers"});
        
        auto two_leptons = trigger_selection.Define("di_muons",
        [](RVec<Muon>& muons, RVec<int> muon_id_medium)
        {
            RVec<Muon> new_muons;
            new_muons.reserve(muons.size());
            
            for (auto i = 0; i < muons.size(); i++)
            {
                if  ((muons[i].muon_pt/1e3 <= 10) ||
                     (muons[i].muon_pt/1e3 <= 15 && muon_id_medium[i] == 1) ||
                     (abs(muons[i].muon_eta) >= 2.5) ||
                     (muon_id_medium[i] != 1))
                {
                    continue;
                }
                new_muons.push_back(muons[i]);
            }
            
            return new_muons;

        },{"muons", "muon_id_medium"})
        .Define("di_electrons",
        [](RVec<Electron> electrons)
        {
            electrons.erase(std::remove_if(electrons.begin(),electrons.end(),
            [](Electron& ep)
            {
                return ((ep.electron_pt/1e3 <= 20) || (abs(ep.electron_eta) >= 2.37)
                        || ((1.37 < abs(ep.electron_eta)) && (abs(ep.electron_eta) < 1.52))
                        || (ep.electron_id_medium != 1));
                

            }), electrons.end());

            return electrons;

        },{"electrons"})
        .Filter([&](RVec<Electron>& di_electrons, RVec<Muon>& di_muons)
        {
            if (!comb)
            {
                return
                ((di_electrons.size()==2 && di_muons.empty() && DeltaR(di_electrons[0].Vector(), di_electrons[1].Vector()) > 0.01)
                ||
                 (di_muons.size()==2 && di_electrons.empty() && DeltaR(di_muons[0].Vector(), di_muons[1].Vector()) > 0.01));
            }
            
            if (di_muons.size() < 2)
            {
                if (di_electrons.size() < 2) //need at least 2 electrons
                {
                    return false;
                }
                auto combs = Combinations(di_electrons, 2);
                auto length = combs[0].size();
                for (size_t i=0; i<length; i++)
                {
                    if (DeltaR(di_electrons[combs[0][i]].Vector(), di_electrons[combs[1][i]].Vector()) > 0.01)
                    {
                        electrons_or_muons = "use_electrons";
                        return true;
                    }
                }
                return false;
            }
            
            if (di_electrons.size() < 2)
            {
                if (di_muons.size() < 2) //need at least 2 muons
                {
                    return false;
                }
                auto combs = Combinations(di_muons, 2);
                auto length = combs[0].size();
                for (size_t i=0; i<length; i++)
                {
                    if (DeltaR(di_muons[combs[0][i]].Vector(), di_muons[combs[1][i]].Vector()) > 0.01)
                    {
                        electrons_or_muons = "use_muons";
                        return true;
                    }
                }
                return false;
            }
            
            return false;
            
            
        }, {"di_electrons", "di_muons"});
        
        auto opp_charge = two_leptons.Filter([&](RVec<Muon>& di_muons, RVec<Electron>& di_electrons)
        {
            if (!comb)
            {
                return ((di_electrons.empty() && di_muons[0].muon_charge*di_muons[1].muon_charge < 0)
                        ||
                       (di_muons.empty() && di_electrons[0].electron_charge*di_electrons[1].electron_charge < 0));
            }
            
            if (di_muons.size() < 2 && electrons_or_muons == "use_electrons")
            {
                if (di_electrons.size() < 2) //need at least 2 electrons
                {
                    return false;
                }
                auto combs = Combinations(di_electrons, 2);
                auto length = combs[0].size();
                for (size_t i=0; i<length; i++)
                {
                    if ((DeltaR(di_electrons[combs[0][i]].Vector(), di_electrons[combs[1][i]].Vector()) > 0.01)
                        && (di_electrons[combs[0][i]].electron_charge * di_electrons[combs[1][i]].electron_charge < 0))
                    {
                        return true;
                    }
                }
                return false;
            }
            
            if (di_electrons.size() < 2 && electrons_or_muons == "use_muons")
            {
                if (di_muons.size() < 2) //need at least 2 muons
                {
                    return false;
                }
                auto combs = Combinations(di_muons, 2);
                auto length = combs[0].size();
                for (size_t i=0; i<length; i++)
                {
                    if ((DeltaR(di_muons[combs[0][i]].Vector(), di_muons[combs[1][i]].Vector()) > 0.01)
                        && (di_muons[combs[0][i]].muon_charge * di_muons[combs[1][i]].muon_charge < 0))
                    {
                        return true;
                    }
                }
                return false;
            }
            
            return false;
            
        }, {"di_muons", "di_electrons"});
        
        auto leading_pt = opp_charge.Filter([&](RVec<Muon>& di_muons, RVec<Electron>& di_electrons)
        {
            if (!comb)
            {
                return (
                        (di_electrons.empty() && ((di_muons[0].muon_pt >= 20e3 && di_muons[1].muon_pt >= 27e3) || (di_muons[1].muon_pt >= 20e3 && di_muons[0].muon_pt >= 27e3)))
                ||
                        (di_muons.empty() && ((di_electrons[0].electron_pt >= 20e3 && di_electrons[1].electron_pt >= 27e3) || (di_electrons[1].electron_pt >= 20e3 && di_electrons[0].electron_pt >= 27e3)))
                        );
            }
            
            if (di_muons.size() < 2 && electrons_or_muons == "use_electrons")
            {
                if (di_electrons.size() < 2) //need at least 2 electrons
                {
                    return false;
                }
                auto combs = Combinations(di_electrons, 2);
                auto length = combs[0].size();
                for (size_t i=0; i<length; i++)
                {
                    if ((DeltaR(di_electrons[combs[0][i]].Vector(), di_electrons[combs[1][i]].Vector()) > 0.01)
                        && (di_electrons[combs[0][i]].electron_charge * di_electrons[combs[1][i]].electron_charge < 0)
                        && ((di_electrons[combs[0][i]].electron_pt >= 20e3 && di_electrons[combs[1][i]].electron_pt >= 27e3) || (di_electrons[combs[1][i]].electron_pt >= 20e3 && di_electrons[combs[0][i]].electron_pt >= 27e3)))
                    {
                        return true;
                    }
                }
                return false;
            }
            
            if (di_electrons.size() < 2 && electrons_or_muons == "use_muons")
            {
                if (di_muons.size() < 2) //need at least 2 muons
                {
                    return false;
                }
                auto combs = Combinations(di_muons, 2);
                auto length = combs[0].size();
                for (size_t i=0; i<length; i++)
                {
                    if ((DeltaR(di_muons[combs[0][i]].Vector(), di_muons[combs[1][i]].Vector()) > 0.01)
                        && (di_muons[combs[0][i]].muon_charge * di_muons[combs[1][i]].muon_charge < 0)
                        && ((di_muons[combs[0][i]].muon_pt >= 20e3 && di_muons[combs[1][i]].muon_pt >= 27e3) || (di_muons[combs[1][i]].muon_pt >= 20e3 && di_muons[combs[0][i]].muon_pt >= 27e3)))
                    {
                        return true;
                    }
                }
                return false;
            }
            
            return false;
                    
        }, {"di_muons", "di_electrons"});
        
//        auto delta_R = leading_pt.Filter([] (RVec<Electron>& electrons)
//        {
//            return (DeltaR(electrons[0].Vector(), electrons[1].Vector()) > 0.01);
//        }, {"di_electrons"});
        
        auto same_flavour = leading_pt.Filter([&] (RVec<Muon>& di_muons, RVec<Electron>& di_electrons)
        {
            if (!comb)
            {
                return (di_muons.size()==2 && di_electrons.empty() || di_electrons.size()==2 && di_muons.empty());
            }
            return (di_muons.size()>=2 && di_electrons.empty() || di_electrons.size()>=2 && di_muons.empty());
            
        }, {"di_muons", "di_electrons"});
        
        auto dilep = same_flavour.Define("dilep",[&] (RVec<Muon>& di_muons, RVec<Electron>& di_electrons)
        {
            if (!comb)
            {
                if (di_muons.empty())
                {
                    return (di_electrons[0].Vector() + di_electrons[1].Vector());
                }
                return (di_muons[0].Vector() + di_muons[1].Vector());
            }
            
            if (di_muons.size() < 2 && electrons_or_muons == "use_electrons")
            {
                auto combs = Combinations(di_electrons, 2);
                auto length = combs[0].size();
                for (size_t i=0; i<length; i++)
                {
                    if ((DeltaR(di_electrons[combs[0][i]].Vector(), di_electrons[combs[1][i]].Vector()) > 0.01)
                        && (di_electrons[combs[0][i]].electron_charge * di_electrons[combs[1][i]].electron_charge < 0)
                        && ((di_electrons[combs[0][i]].electron_pt >= 20e3 && di_electrons[combs[1][i]].electron_pt >= 27e3) || (di_electrons[combs[1][i]].electron_pt >= 20e3 && di_electrons[combs[0][i]].electron_pt >= 27e3)))
                    {
                        return di_electrons[combs[0][i]].Vector() + di_electrons[combs[1][i]].Vector();
                    }
                }
            }
            
            if (di_electrons.size() < 2 && electrons_or_muons == "use_muons")
            {
                auto combs = Combinations(di_muons, 2);
                auto length = combs[0].size();
                for (size_t i=0; i<length; i++)
                {
                    if ((DeltaR(di_muons[combs[0][i]].Vector(), di_muons[combs[1][i]].Vector()) > 0.01)
                        && (di_muons[combs[0][i]].muon_charge * di_muons[combs[1][i]].muon_charge < 0)
                        && ((di_muons[combs[0][i]].muon_pt >= 20e3 && di_muons[combs[1][i]].muon_pt >= 27e3) || (di_muons[combs[1][i]].muon_pt >= 20e3 && di_muons[combs[0][i]].muon_pt >= 27e3)))
                    {
                        return di_muons[combs[0][i]].Vector() + di_muons[combs[1][i]].Vector();
                    }
                }
            }
            
            return (di_muons[0].Vector() + di_muons[1].Vector());
            
        }, {"di_muons", "di_electrons" });
        
        auto mass = dilep.Filter([] (PtEtaPhiEVector& dilep)
        {
            auto mass = dilep.M()/1e3;
            return ((mass >= 81) && (mass <= 101));
        }, {"dilep"});
        
        auto ptCut = mass.Filter([] (PtEtaPhiEVector& dilep)
        {
            auto pT = dilep.Pt()/1e3;
            return pT > 10;
        }, {"dilep"});
        
        auto photonPtDeltaR = ptCut.Define("photonPtDeltaR",
        [&](RVec<Photon> photons)
        {
            photons.erase(std::remove_if(photons.begin(),photons.end(),
            [](Photon& x)
            {
                return ((abs(x.photon_eta) >= 2.37) || (x.photon_pt <= 10e3) || (abs(x.photon_eta) > 1.37 && abs(x.photon_eta) < 1.52));
                
            }), photons.end());
            return photons;
        }, {"photons"}).Filter(
       [&](RVec<Photon>& reco_photons_matched)
       {
            if (reco_photons_matched.size() < 2)
            {
              return false;
            }
            RVec<Photon> x;
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
                    x = {reco_photons_matched[combs[0][i]], reco_photons_matched[combs[1][i]]};
                }
            }
            return (chosen_delta_r < 1.5 && pt1 > 10e3 && pt2 > 10e3);

       }, {"photonPtDeltaR"});
                
        auto X_window = photonPtDeltaR.Define("chosen_two",
        [](RVec<Photon>& reco_photons_matched)
        {
            RVec<Photon> x;
            if (reco_photons_matched.size() < 2)
            {
                return x;
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
                    x = {reco_photons_matched[combs[0][i]], reco_photons_matched[combs[1][i]]};
                }
            }
            if (pt1 > 10e3 && pt2 > 10e3 && best_X > 0.96 && best_X < 1.2 && chosen_delta_r < 1.5)
            {
                return x;
            }
            x.clear();
            return x;
        }, {"photonPtDeltaR"}).Filter(
        [&](RVec<Photon>& reco_photons_matched)
        {
            return (reco_photons_matched.size()==2);
        }, {"chosen_two"});
                
        auto SR = X_window.Filter(
        [](RVec<Photon>& photons, RVec<Muon>& muons, RVec<Electron>& electrons)
        {
            PtEtaPhiEVector photonVec = photons[0].Vector() + photons[1].Vector();
            PtEtaPhiEVector leptonVec;
            if (muons.empty())
            {
                leptonVec = electrons[0].Vector() + electrons[1].Vector();
            }
            else
            {
                leptonVec = muons[0].Vector() + muons[1].Vector();
            }
            auto mass = (photonVec+leptonVec).M()/1e3;
            return ((mass >= 110) && (mass <= 140));
        },{"chosen_two","di_muons","di_electrons"});
                
        auto SR_ID = SR.Filter(
        [](RVec<Photon>& photons)
        {
            return (photons[0].photon_id_loose && photons[1].photon_id_loose);
        },{"chosen_two"});
        
        Nodes.push_back(df.Count());
        Nodes.push_back(ptCut.Count());
        Nodes.push_back(photonPtDeltaR.Count());
        Nodes.push_back(X_window.Count());
        Nodes.push_back(SR.Count());
        Nodes.push_back(SR_ID.Count());
    }
    
    ROOT::RDF::RunGraphs(Nodes); // running all computation nodes concurrently
    
//    0    1    2    3    4    5      Z-gamma
//
//    6    7    8    9    10   11     Z-gamma
//
//    12   13   14   15   16   17     Z-gamma
//
//    18   19   20   21   22   23     ma1
//
//    24   25   26   27   28   29     ma5
//
//    30   31   32   33   34   35     data
//
//    36   37   38   39   40   41     Z-jets
//
//    42   43   44   45   46   47     Z-jets
//
//    48   49   50   51   52   53     Z-jets
//
//    54   55   56   57   58   59     Z-jets
//
//    60   61   62   63   64   65     Z-jets
//
//    66   67   68   69   70   71     Z-jets
//
//    72   73   74   75   76   77     Z-jets
//
//    78   79   80   81   82   83     Z-jets
//
//    84   85   86   87   88   89     Z-jets

    std::cout << R"--(\section*{Table 11 Signal Ratios})--" << '\n';
    std::cout << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
    std::cout << R"--(\begin{tabular}{|c|c|c|c|c|c|c|})--" << '\n';
    std::cout << R"--(\hline)--" << '\n';
    std::cout << R"--(Sample & $\frac{\text{pass preselection (PS)}}{\text{pass preselection (PS)}}$ & $\frac{\text{photon } p_T \text{ + } \Delta R_{\gamma\gamma} \text{ cut }}{\text{pass preselection (PS)}}$ & $\frac{X \text{ window}}{\text{pass preselection (PS)}}$ & $\frac{\text{SR}}{\text{pass preselection (PS)}}$ & $\frac{\text{SR-ID}}{\text{pass preselection (PS)}}$
           \\ \hline )--" << '\n';
    
    ss << R"--(\section*{Table 11 Signal Ratios})--" << '\n';
    ss << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
    ss << R"--(\begin{tabular}{|c|c|c|c|c|c|c|})--" << '\n';
    ss << R"--(\hline)--" << '\n';
    ss << R"--(Sample & $\frac{\text{pass preselection (PS)}}{\text{pass preselection (PS)}}$ & $\frac{\text{photon } p_T \text{ + } \Delta R_{\gamma\gamma} \text{ cut }}{\text{pass preselection (PS)}}$ & $\frac{X \text{ window}}{\text{photon } p_T \text{ + } \Delta R_{\gamma\gamma} \text{ cut }}$ & $\frac{\text{SR}}{X \text{ window}}$ & $\frac{\text{SR-ID}}{\text{SR}}$
           \\ \hline )--" << '\n';
    
    for (int i=0, j=0; (i<1 && j <= 0); i++, j+=6)
    {
        os << Samples[i];
        
        os  << " & " << *Nodes[j].GetResultPtr<ULong64_t>()
            << " & " << *Nodes[j+1].GetResultPtr<ULong64_t>()
            << " & " << *Nodes[j+2].GetResultPtr<ULong64_t>()
            << " & " << *Nodes[j+3].GetResultPtr<ULong64_t>()
            << " & " << *Nodes[j+4].GetResultPtr<ULong64_t>()
            << " & " << *Nodes[j+5].GetResultPtr<ULong64_t>()
            << R"--( \\ \hline )--" << '\n';
            
            
            
            if (i==0) //5 GeV
            {
                std::cout << Samples[i] << " (Me) & " << 1 << " & " <<
                static_cast<double>(*Nodes[j+2].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
                << " & " <<
                static_cast<double>(*Nodes[j+3].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
                << " & " <<
                static_cast<double>(*Nodes[j+4].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
                << " & " <<
                static_cast<double>(*Nodes[j+5].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
                << R"--( \\ \hline )--" << '\n';
                
                std::cout << Samples[i] << " (Paper) & " << 1 << " & " <<
                3245.86 / 19536.68
                << " & " <<
                3052.81 / 19536.68
                << " & " <<
                2941.28 / 19536.68
                << " & " <<
                1959.68 / 19536.68
                << R"--( \\ \hline )--" << '\n';
                
                ss << Samples[i] << " (Me) & " << 1 << " & " <<
                static_cast<double>(*Nodes[j+2].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
                << " & " <<
                static_cast<double>(*Nodes[j+3].GetResultPtr<ULong64_t>()) / *Nodes[j+2].GetResultPtr<ULong64_t>()
                << " & " <<
                static_cast<double>(*Nodes[j+4].GetResultPtr<ULong64_t>()) / *Nodes[j+3].GetResultPtr<ULong64_t>()
                << " & " <<
                static_cast<double>(*Nodes[j+5].GetResultPtr<ULong64_t>()) / *Nodes[j+4].GetResultPtr<ULong64_t>()
                << R"--( \\ \hline )--" << '\n';
                
                ss << Samples[i] << " (Paper) & " << 1 << " & " <<
                3245.86 / 19536.68
                << " & " <<
                3052.81 / 3245.86
                << " & " <<
                2941.28 / 3052.81
                << " & " <<
                1959.68 / 2941.28
                << R"--( \\ \hline )--" << '\n';
            }
    }
    
    std::cout << R"--(\end{tabular}})--" << "\n\n\n";
    
    ss << R"--(\end{tabular}})--" << "\n\n\n";
    std::cout << ss.str() << '\n';
    
    
    
    os << R"--(\end{tabular}})--" << '\n';
    std::cout << os.str() << '\n';
}
/*
void Table11()
{
    std::vector<std::vector<std::string>> input_filenames = {
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/muon17.3.root", "muon17.3part2.root"}
    };
    
    constexpr std::array<const char*,1> Samples = {R"--(Signal $m_{\text{A}}$ = 5 GeV)--"};
    
    std::vector<ROOT::RDF::RResultHandle> Nodes;

    std::ostringstream os, ss;
    os << R"--(\section*{Table 11})--" << '\n';
    os << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
    os << R"--(\begin{tabular}{|c|c|c|c|c|c|c|})--" << '\n';
    os << R"--(\hline)--" << '\n';
    os << R"--(Sample & Total Events & pass preselection (PS) & photon $p_T$ + $\Delta R_{\gamma\gamma}$ cut & $X$ window & SR & SR-ID
           \\ \hline )--" << '\n';
    
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
    
 

    for (auto& file: input_filenames)
    {
        SchottDataFrame df(MakeRDF(file, 8));
        
        auto trigger_selection = df.Filter(
        [](const RVec<std::string>& trigger_passed_triggers)
        {
            bool trigger_found = (std::find_first_of(trigger_passed_triggers.begin(), trigger_passed_triggers.end(), triggers.begin(), triggers.end()) != trigger_passed_triggers.end());

            if (!trigger_found)
            {
             return false;
            }
            return true;

        }, {"trigger_passed_triggers"});
        
        auto two_leptons = df.Define("di_muons",
        [](RVec<Muon>& muons, RVec<int> muon_id_medium)
        {
            RVec<Muon> new_muons;
            new_muons.reserve(muons.size());
            
            for (auto i = 0; i < muons.size(); i++)
            {
                if  ((muons[i].muon_pt/1e3 <= 10) ||
                     (muons[i].muon_pt/1e3 <= 15 && muon_id_medium[i] == 1) ||
                     (abs(muons[i].muon_eta) >= 2.5) ||
                     (muon_id_medium[i] != 1))
                {
                    continue;
                }
                new_muons.push_back(muons[i]);
            }
            
            return new_muons;

        },{"muons", "muon_id_medium"}).Filter([](RVec<Electron>& electrons, RVec<Muon>& di_muons)
        {            
            return (di_muons.size()==2 && electrons.empty() && DeltaR(di_muons[0].Vector(), di_muons[1].Vector()) > 0.01);
            
        }, {"electrons", "di_muons"});
        
        auto opp_charge = two_leptons.Filter([](RVec<Muon>& di_muons)
        {
            return (di_muons[0].muon_charge*di_muons[1].muon_charge < 0);
            
        }, {"di_muons"});
        
        auto leading_pt = opp_charge.Filter([](RVec<Muon>& muons)
        {
            return ((muons[0].muon_pt >= 20e3 && muons[1].muon_pt >= 27e3) || (muons[1].muon_pt >= 20e3 && muons[0].muon_pt >= 27e3));
        }, {"di_muons"});
        
//        auto delta_R = leading_pt.Filter([] (RVec<Electron>& electrons)
//        {
//            return (DeltaR(electrons[0].Vector(), electrons[1].Vector()) > 0.01);
//        }, {"di_electrons"});
        
        auto same_flavour = leading_pt.Filter([] (RVec<Muon>& muons)
        {
            return true; //abs(electrons[0].electron_pdg_id) == abs(electrons[1].electron_pdg_id) == 11;
        }, {"di_muons"});
        
        auto dilep = same_flavour.Define("dilep",[] (RVec<Muon>& muons)
        {
            return (muons[0].Vector() + muons[1].Vector());
        }, {"di_muons"});
        
        auto mass = dilep.Filter([] (PtEtaPhiEVector& dilep)
        {
            auto mass = dilep.M()/1e3;
            return ((mass >= 81) && (mass <= 101));
        }, {"dilep"});
        
        auto ptCut = mass.Filter([] (PtEtaPhiEVector& dilep)
        {
            auto pT = dilep.Pt()/1e3;
            return pT > 10;
        }, {"dilep"});
        
        auto photonPtDeltaR = ptCut.Define("photonPtDeltaR",
        [&](RVec<Photon> photons)
        {
            photons.erase(std::remove_if(photons.begin(),photons.end(),
            [](Photon& x)
            {
                return ((abs(x.photon_eta) >= 2.37) || (x.photon_pt <= 10e3) || (abs(x.photon_eta) > 1.37 && abs(x.photon_eta) < 1.52));
                
            }), photons.end());
            return photons;
        }, {"photons"}).Filter(
       [&](RVec<Photon>& reco_photons_matched)
       {
            if (reco_photons_matched.size() < 2)
            {
              return false;
            }
            RVec<Photon> x;
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
                    x = {reco_photons_matched[combs[0][i]], reco_photons_matched[combs[1][i]]};
                }
            }
            return (chosen_delta_r < 1.5 && pt1 > 10e3 && pt2 > 10e3);

       }, {"photonPtDeltaR"});
                
        auto X_window = photonPtDeltaR.Define("chosen_two",
        [](RVec<Photon>& reco_photons_matched)
        {
            RVec<Photon> x;
            if (reco_photons_matched.size() < 2)
            {
                return x;
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
                    x = {reco_photons_matched[combs[0][i]], reco_photons_matched[combs[1][i]]};
                }
            }
            if (pt1 > 10e3 && pt2 > 10e3 && best_X > 0.96 && best_X < 1.2 && chosen_delta_r < 1.5)
            {
                return x;
            }
            x.clear();
            return x;
        }, {"photonPtDeltaR"}).Filter(
        [&](RVec<Photon>& reco_photons_matched)
        {
            return (reco_photons_matched.size()==2);
        }, {"chosen_two"});
                
        auto SR = X_window.Filter(
        [](RVec<Photon>& photons, RVec<Muon>& muons)
        {
            PtEtaPhiEVector photonVec = photons[0].Vector() + photons[1].Vector();
            PtEtaPhiEVector muonVec = muons[0].Vector() + muons[1].Vector();
            auto mass = (photonVec+muonVec).M()/1e3;
            return ((mass >= 110) && (mass <= 140));
        },{"chosen_two","di_muons"});
                
        auto SR_ID = SR.Filter(
        [](RVec<Photon>& photons)
        {
            return (photons[0].photon_id_loose && photons[1].photon_id_loose);
        },{"chosen_two"});
        
        Nodes.push_back(df.Count());
        Nodes.push_back(ptCut.Count());
        Nodes.push_back(photonPtDeltaR.Count());
        Nodes.push_back(X_window.Count());
        Nodes.push_back(SR.Count());
        Nodes.push_back(SR_ID.Count());
    }
    
    ROOT::RDF::RunGraphs(Nodes); // running all computation nodes concurrently
    
//    0    1    2    3    4    5      Z-gamma
//
//    6    7    8    9    10   11     Z-gamma
//
//    12   13   14   15   16   17     Z-gamma
//
//    18   19   20   21   22   23     ma1
//
//    24   25   26   27   28   29     ma5
//
//    30   31   32   33   34   35     data
//
//    36   37   38   39   40   41     Z-jets
//
//    42   43   44   45   46   47     Z-jets
//
//    48   49   50   51   52   53     Z-jets
//
//    54   55   56   57   58   59     Z-jets
//
//    60   61   62   63   64   65     Z-jets
//
//    66   67   68   69   70   71     Z-jets
//
//    72   73   74   75   76   77     Z-jets
//
//    78   79   80   81   82   83     Z-jets
//
//    84   85   86   87   88   89     Z-jets

    std::cout << R"--(\section*{Table 11 Signal Ratios})--" << '\n';
    std::cout << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
    std::cout << R"--(\begin{tabular}{|c|c|c|c|c|c|c|})--" << '\n';
    std::cout << R"--(\hline)--" << '\n';
    std::cout << R"--(Sample & $\frac{\text{pass preselection (PS)}}{\text{pass preselection (PS)}}$ & $\frac{\text{photon } p_T \text{ + } \Delta R_{\gamma\gamma} \text{ cut }}{\text{pass preselection (PS)}}$ & $\frac{X \text{ window}}{\text{pass preselection (PS)}}$ & $\frac{\text{SR}}{\text{pass preselection (PS)}}$ & $\frac{\text{SR-ID}}{\text{pass preselection (PS)}}$
           \\ \hline )--" << '\n';
    
    ss << R"--(\section*{Table 11 Signal Ratios})--" << '\n';
    ss << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
    ss << R"--(\begin{tabular}{|c|c|c|c|c|c|c|})--" << '\n';
    ss << R"--(\hline)--" << '\n';
    ss << R"--(Sample & $\frac{\text{pass preselection (PS)}}{\text{pass preselection (PS)}}$ & $\frac{\text{photon } p_T \text{ + } \Delta R_{\gamma\gamma} \text{ cut }}{\text{pass preselection (PS)}}$ & $\frac{X \text{ window}}{\text{photon } p_T \text{ + } \Delta R_{\gamma\gamma} \text{ cut }}$ & $\frac{\text{SR}}{X \text{ window}}$ & $\frac{\text{SR-ID}}{\text{SR}}$
           \\ \hline )--" << '\n';
    
    for (int i=0, j=0; (i<1 && j <= 0); i++, j+=6)
    {
        os << Samples[i];
        
        os  << " & " << *Nodes[j].GetResultPtr<ULong64_t>()
            << " & " << *Nodes[j+1].GetResultPtr<ULong64_t>()
            << " & " << *Nodes[j+2].GetResultPtr<ULong64_t>()
            << " & " << *Nodes[j+3].GetResultPtr<ULong64_t>()
            << " & " << *Nodes[j+4].GetResultPtr<ULong64_t>()
            << " & " << *Nodes[j+5].GetResultPtr<ULong64_t>()
            << R"--( \\ \hline )--" << '\n';
            
            
            
            if (i==0) //5 GeV
            {
                std::cout << Samples[i] << " (Me) & " << 1 << " & " <<
                static_cast<double>(*Nodes[j+2].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
                << " & " <<
                static_cast<double>(*Nodes[j+3].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
                << " & " <<
                static_cast<double>(*Nodes[j+4].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
                << " & " <<
                static_cast<double>(*Nodes[j+5].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
                << R"--( \\ \hline )--" << '\n';
                
                std::cout << Samples[i] << " (Paper) & " << 1 << " & " <<
                3245.86 / 19536.68
                << " & " <<
                3052.81 / 19536.68
                << " & " <<
                2941.28 / 19536.68
                << " & " <<
                1959.68 / 19536.68
                << R"--( \\ \hline )--" << '\n';
                
                ss << Samples[i] << " (Me) & " << 1 << " & " <<
                static_cast<double>(*Nodes[j+2].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
                << " & " <<
                static_cast<double>(*Nodes[j+3].GetResultPtr<ULong64_t>()) / *Nodes[j+2].GetResultPtr<ULong64_t>()
                << " & " <<
                static_cast<double>(*Nodes[j+4].GetResultPtr<ULong64_t>()) / *Nodes[j+3].GetResultPtr<ULong64_t>()
                << " & " <<
                static_cast<double>(*Nodes[j+5].GetResultPtr<ULong64_t>()) / *Nodes[j+4].GetResultPtr<ULong64_t>()
                << R"--( \\ \hline )--" << '\n';
                
                ss << Samples[i] << " (Paper) & " << 1 << " & " <<
                3245.86 / 19536.68
                << " & " <<
                3052.81 / 3245.86
                << " & " <<
                2941.28 / 3052.81
                << " & " <<
                1959.68 / 2941.28
                << R"--( \\ \hline )--" << '\n';
            }
    }
    
    std::cout << R"--(\end{tabular}})--" << "\n\n\n";
    
    ss << R"--(\end{tabular}})--" << "\n\n\n";
    std::cout << ss.str() << '\n';
    
    
    
    os << R"--(\end{tabular}})--" << '\n';
    std::cout << os.str() << '\n';
}
*/
void CutFlow()
{
    auto start_time = Clock::now();
    
//    Table3();
//    Table8();
//    Table11();
    Table11MuonElectron();
    
    auto end_time = Clock::now();
    std::cout << "Time difference: "
       << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count()/1e9 << " seconds" << std::endl;
    
}


int main()
{
    CutFlow();
}
