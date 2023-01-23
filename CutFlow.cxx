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

void Table3()
{
    std::vector<std::string> input_filenames = {
        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617070._000001.LGNTuple.root", "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617064._000001.LGNTuple.root",
        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617074._000001.LGNTuple.root",
        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root",
        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_data_test.root",
        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root"
    };
    
    std::vector<std::string> cutFlows;
    cutFlows.reserve(input_filenames.size());
    constexpr std::array<const char*,6> Samples = {R"--(pty2\_9\_17)--", R"--(pty\_17\_myy\_0\_80)--", R"--(pty\_17\_myy\_80)--", R"--(Signal $m_{\text{A}}$ = 1 GeV)--", "Data", R"--(Signal $m_{\text{A}}$ = 5 GeV)--"};
    
    std::ostringstream os;
    os << R"--(\section*{Table 3})--" << '\n';
    os << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
    os << R"--(\begin{tabular}{|c|c|c|c|c|c|c|c|})--" << '\n';
    os << R"--(\hline)--" << '\n';
    os << R"--(Sample & Before Preselection & 2 leptons
          & Opposite Charge & $p_{T}^{\text{leading}} > 27$ GeV, \; $p_{T}^{\text{sub-leading}} > 20$ GeV & $\Delta R > 0.01$ & dilep mass cut & dilep $p_{T}$ cut \\ \hline )--" << '\n';
    
    double beforePreselecZGamma = 0, twoLeptonsZGamma = 0, oppChargeZGamma = 0, leadingPtZGamma = 0, deltaRZGamma = 0, MassZGamma = 0, ptCutZGamma = 0;

    std::vector<ROOT::RDF::RResultHandle> Nodes;
    
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
        
        auto leading_pt = opp_charge.Filter([](RVec<Electron>& electrons)
        {
            return ((electrons[0].electron_pt > 20e3 && electrons[1].electron_pt > 27e3) || (electrons[1].electron_pt > 20e3 && electrons[0].electron_pt > 27e3));
        }, {"di_electrons"});
        
        auto delta_R = leading_pt.Filter([] (RVec<Electron>& electrons)
        {
            return (DeltaR(electrons[0].Vector(), electrons[1].Vector()) > 0.01);
        }, {"di_electrons"});
        
        auto mass = delta_R.Filter([] (RVec<Electron>& electrons)
        {
            auto mass = (electrons[0].Vector() + electrons[1].Vector()).M()/1e3;
            return ((mass >= 81) && (mass <= 101));
        }, {"di_electrons"});
        
        auto pt_cut = mass.Filter([] (RVec<Electron>& electrons)
        {
            auto pT = (electrons[0].Vector() + electrons[1].Vector()).Pt()/1e3;
            return pT > 10;
        }, {"di_electrons"});
        
        Nodes.push_back(df.Count());
        Nodes.push_back(two_leptons.Count());
        Nodes.push_back(opp_charge.Count());
        Nodes.push_back(leading_pt.Count());
        Nodes.push_back(delta_R.Count());
        Nodes.push_back(mass.Count());
        Nodes.push_back(pt_cut.Count());
    }
    
    ROOT::RDF::RunGraphs(Nodes); // running all computation nodes concurrently
    
    for (int i=0, j=0; i<6 && j <= 35; i++, j+=7)
    {
        os << Samples[i]
        << " & " << *Nodes[j].GetResultPtr<ULong64_t>()
        << " & " << *Nodes[j+1].GetResultPtr<ULong64_t>()
        << " & " << *Nodes[j+2].GetResultPtr<ULong64_t>()
        << " & " << *Nodes[j+3].GetResultPtr<ULong64_t>()
        << " & " << *Nodes[j+4].GetResultPtr<ULong64_t>()
        << " & " << *Nodes[j+5].GetResultPtr<ULong64_t>()
        << " & " << *Nodes[j+6].GetResultPtr<ULong64_t>()
        << R"--( \\ \hline )--" << '\n';
    }
    
    for (int i = 0; i <= 14; i += 7)
    {
        beforePreselecZGamma += *Nodes[i].GetResultPtr<ULong64_t>();
        twoLeptonsZGamma += *Nodes[i+1].GetResultPtr<ULong64_t>();
        oppChargeZGamma += *Nodes[i+2].GetResultPtr<ULong64_t>();
        leadingPtZGamma += *Nodes[i+3].GetResultPtr<ULong64_t>();
        deltaRZGamma += *Nodes[i+4].GetResultPtr<ULong64_t>();
        MassZGamma += *Nodes[i+5].GetResultPtr<ULong64_t>();
        ptCutZGamma += *Nodes[i+6].GetResultPtr<ULong64_t>();
    }

    os << R"--(Total $Z\gamma$ & )--" << beforePreselecZGamma << " & " << twoLeptonsZGamma
    << " & " << oppChargeZGamma << " & " << leadingPtZGamma << " & "
    << deltaRZGamma << " & " << MassZGamma << " & " << ptCutZGamma
    << R"--( \\ \hline )--" << '\n';
    
    os << R"--(\end{tabular}})--" << '\n';
    cutFlows.push_back(os.str());
    for (auto& i: cutFlows)
    {
        std::cout << i << "\n\n";
    }
}

void Table8()
{
    std::vector<std::string> input_filenames = {
        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617070._000001.LGNTuple.root", "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617064._000001.LGNTuple.root",
        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617074._000001.LGNTuple.root",
        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root",
        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_data_test.root",
        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root"
    };
    
    constexpr std::array<const char*,6> Samples = {R"--(pty2\_9\_17)--", R"--(pty\_17\_myy\_0\_80)--", R"--(pty\_17\_myy\_80)--", R"--(Signal $m_{\text{A}}$ = 1 GeV)--", "Data", R"--(Signal $m_{\text{A}}$ = 5 GeV)--"};
    
    double totalEvents = 0, resolvedEvents = 0, SREvents = 0, SBEvents = 0;
    std::vector<std::string> cutFlows;
    cutFlows.reserve(input_filenames.size());
    
    std::vector<ROOT::RDF::RResultHandle> Nodes;
    
    std::ostringstream os;
    os << R"--(\section*{Table 8})--" << '\n';
    os << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
    os << R"--(\begin{tabular}{|c|c|c|c|c|})--" << '\n';
    os << R"--(\hline)--" << '\n';
    os << R"--(Sample & Total Events & PS: Resolved & SB & SR
           \\ \hline )--" << '\n';
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
        
        auto resolved = ptCut.Define("photons_pass_cuts",
        [&](RVec<Photon> photons)
        {
            photons.erase(std::remove_if(photons.begin(),photons.end(),
            [](Photon& x)
            {
                return ((abs(x.photon_eta) >= 2.37) || (abs(x.photon_eta) > 1.37 && abs(x.photon_eta) < 1.52) || (!x.photon_id_loose));

            }), photons.end());
            return photons;
        }, {"photons"}).Filter(
        [&](RVec<Photon>& reco_photons_matched)
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

        auto SR = resolved.Filter(
        [&](RVec<Electron>& electrons, RVec<Photon>& reco_photons_matched)
        {
            auto mass = ((electrons[0].Vector()+electrons[1].Vector()).M() + (reco_photons_matched[0].Vector()+reco_photons_matched[1].Vector()).M())/1e3;
            return ((mass > 110) && (mass < 140));
        }, {"di_electrons", "photons_pass_cuts"});
        
        auto SB = resolved.Filter(
        [&](RVec<Electron>& electrons, RVec<Photon>& reco_photons_matched)
        {
            auto mass = ((electrons[0].Vector()+electrons[1].Vector()).M() + (reco_photons_matched[0].Vector()+reco_photons_matched[1].Vector()).M())/1e3;
            return (!((mass > 110) && (mass < 140)));
        }, {"di_electrons", "photons_pass_cuts"});
        
        Nodes.push_back(df.Count());
        Nodes.push_back(resolved.Count());
        Nodes.push_back(SB.Count());
        Nodes.push_back(SR.Count());
    }
    
    ROOT::RDF::RunGraphs(Nodes); // running all computation nodes concurrently
    
    for (int i=0, j=0; i<6 && j <= 20; i++, j+=4)
    {
        os << Samples[i]
        << " & " << *Nodes[j].GetResultPtr<ULong64_t>()
        << " & " << *Nodes[j+1].GetResultPtr<ULong64_t>()
        << " & " << *Nodes[j+2].GetResultPtr<ULong64_t>()
        << " & " << *Nodes[j+3].GetResultPtr<ULong64_t>()
        << R"--( \\ \hline )--" << '\n';
    }
    
    for (int i = 0; i <= 8; i += 4)
    {
        totalEvents += *Nodes[i].GetResultPtr<ULong64_t>();
        resolvedEvents += *Nodes[i+1].GetResultPtr<ULong64_t>();
        SBEvents += *Nodes[i+2].GetResultPtr<ULong64_t>();
        SREvents += *Nodes[i+3].GetResultPtr<ULong64_t>();
    }
    
    os << R"--(Total $Z\gamma$ & )--" << totalEvents << " & " << resolvedEvents
    << " & " << SBEvents << " & " << SREvents
    << R"--( \\ \hline )--" << '\n';
    
    os << R"--(\end{tabular}})--" << '\n';
    cutFlows.push_back(os.str());
    for (auto& i: cutFlows)
    {
        std::cout << i << "\n\n";
    }
}

void Table11()
{
    std::vector<std::string> input_filenames = {
        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617070._000001.LGNTuple.root", "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617064._000001.LGNTuple.root",
        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617074._000001.LGNTuple.root",
        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root",
        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_data_test.root",
        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root"
    };
    
    constexpr std::array<const char*,6> Samples = {R"--(pty2\_9\_17)--", R"--(pty\_17\_myy\_0\_80)--", R"--(pty\_17\_myy\_80)--", R"--(Signal $m_{\text{A}}$ = 1 GeV)--", "Data", R"--(Signal $m_{\text{A}}$ = 5 GeV)--"};
    
    std::vector<std::string> cutFlows;
    cutFlows.reserve(input_filenames.size());
    double totalEvents = 0, passPreselection = 0, photonPtDeltaRCount = 0, xWindow = 0,
    srCount = 0, srIDCount=0;
    std::vector<ROOT::RDF::RResultHandle> Nodes;

    std::ostringstream os;
    os << R"--(\section*{Table 11})--" << '\n';
    os << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
    os << R"--(\begin{tabular}{|c|c|c|c|c|c|c|})--" << '\n';
    os << R"--(\hline)--" << '\n';
    os << R"--(Sample & Total Events & pass preselection (PS) & photon $p_T$ + $\Delta R_{\gamma\gamma}$ cut & $X$ window & SR & SR-ID
           \\ \hline )--" << '\n';
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
           auto combs = Combinations(reco_photons_matched, 2);
           size_t length = combs[0].size();
           double delta_r;

           for (size_t i=0; i<length; i++)
           {
               delta_r = DeltaR(reco_photons_matched[combs[0][i]].Vector(), reco_photons_matched[combs[1][i]].Vector());
               if ((delta_r < 1.5))
               {
                   return true;
               }
           }
           return false;

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
            double delta_r, m, pt, X;

            for (size_t i=0; i<length; i++)
            {
                delta_r = DeltaR(reco_photons_matched[combs[0][i]].Vector(), reco_photons_matched[combs[1][i]].Vector());
                m = (reco_photons_matched[combs[0][i]].Vector() + reco_photons_matched[combs[1][i]].Vector()).M();
                pt = (reco_photons_matched[combs[0][i]].Vector() + reco_photons_matched[combs[1][i]].Vector()).Pt();
                X = delta_r*(pt/(2.0*m));
                if ((X > 0.96) && (X < 1.2))
                {
                    x = {reco_photons_matched[combs[0][i]], reco_photons_matched[combs[1][i]]};
                    return x;
                }
            }
            return x;
        }, {"photonPtDeltaR"}).Filter(
        [&](RVec<Photon>& reco_photons_matched)
        {
            return (reco_photons_matched.size()==2);
        }, {"chosen_two"});
                
        auto SR = X_window.Filter(
        [](RVec<Photon>& photons, RVec<Electron>& electrons)
        {
            PtEtaPhiEVector photonVec = photons[0].Vector() + photons[1].Vector();
            PtEtaPhiEVector electronVec = electrons[0].Vector() + electrons[1].Vector();
            auto mass = (photonVec+electronVec).M()/1e3;
            return ((mass >= 110) && (mass <= 140));
        },{"chosen_two","di_electrons"});
                
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
    
    for (int i=0, j=0; i<6 && j <= 36; i++, j+=6)
    {
        os << Samples[i]
        << " & " << *Nodes[j].GetResultPtr<ULong64_t>()
        << " & " << *Nodes[j+1].GetResultPtr<ULong64_t>()
        << " & " << *Nodes[j+2].GetResultPtr<ULong64_t>()
        << " & " << *Nodes[j+3].GetResultPtr<ULong64_t>()
        << " & " << *Nodes[j+4].GetResultPtr<ULong64_t>()
        << " & " << *Nodes[j+5].GetResultPtr<ULong64_t>()
        << R"--( \\ \hline )--" << '\n';
    }
    
    for (int i = 0; i <= 12; i += 6)
    {
        totalEvents += *Nodes[i].GetResultPtr<ULong64_t>();
        passPreselection += *Nodes[i+1].GetResultPtr<ULong64_t>();
        photonPtDeltaRCount += *Nodes[i+2].GetResultPtr<ULong64_t>();
        xWindow += *Nodes[i+3].GetResultPtr<ULong64_t>();
        srCount += *Nodes[i+4].GetResultPtr<ULong64_t>();
        srIDCount += *Nodes[i+5].GetResultPtr<ULong64_t>();
    }
    
    os << R"--(Total $Z\gamma$ & )--" << totalEvents << " & " << passPreselection
    << " & " << photonPtDeltaRCount << " & " << xWindow << " & " << srCount
    << " & " << srIDCount
    << R"--( \\ \hline )--" << '\n';
    
    os << R"--(\end{tabular}})--" << '\n';
    cutFlows.push_back(os.str());
    for (auto& i: cutFlows)
    {
        std::cout << i << "\n\n";
    }
}

void CutFlow()
{
    auto start_time = Clock::now();
    
    Table3();
    Table8();
    Table11();
    
    auto end_time = Clock::now();
    std::cout << "Time difference: "
       << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count()/1e9 << " nanoseconds" << std::endl;
    
}


int main()
{
    CutFlow();
}
