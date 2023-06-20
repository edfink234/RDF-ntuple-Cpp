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

void Fig34()
{
    std::vector<std::string> input_filenames = {
        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_data_test.root"
    };
    
    SchottDataFrame df(MakeRDF(input_filenames, 8));
    
    auto two_leptons = df
    .Filter([](const RVec<std::string>& trigger_passed_triggers)
    {
        bool trigger_found = (std::find_first_of(trigger_passed_triggers.begin(), trigger_passed_triggers.end(), triggers.begin(), triggers.end()) != trigger_passed_triggers.end());

        if (!trigger_found)
        {
            return false; //this event is filtered out
        }
        return true; //this event is kept because the trigger was found in its `trigger_passed_triggers` branch entry

    }, {"trigger_passed_triggers"})
    .Define("di_electrons", //the events that pass will have exactly 2 electrons that pass the following
    [](RVec<Electron> electrons)
    {
        //keep the electrons in each event that have pt > 20 GeV, |η| < 2.37,
        //|η| not between 1.37 and 1.52, and that satisfy a medium id criteria `electron_id_medium`
        electrons.erase(std::remove_if(electrons.begin(),electrons.end(),
        [](Electron& ep)
        {
             return (!((ep.electron_pt/1e3 > 20) && (std::abs(ep.electron_eta) < 2.37) &&
             (!((1.37 < std::abs(ep.electron_eta)) && (std::abs(ep.electron_eta) < 1.52)))
             && (ep.electron_id_medium == 1)));

        }), electrons.end());

        return electrons;

    },{"electrons"}).Filter(
    [](RVec<Muon>& muons, RVec<Electron> electrons)
    {
        return (electrons.size()==2 && muons.empty()); //keep events which have exactly 2 electrons for di_electrons and no muons

    }, {"muons", "di_electrons"});
    
    //new dataframe node: contains only the events from `two_leptons` whose electrons in the `di_electrons` branch have opposite charge
    auto opp_charge = two_leptons
    .Filter([](RVec<Electron>& electrons)
    {
        return (electrons[0].electron_charge*electrons[1].electron_charge < 0);

    }, {"di_electrons"});

    //new dataframe node: contains only the events from `opp_charge` that have 1 electron with pt > 20 GeV and the other with pt > 27 GeV
    auto leadingPt = opp_charge.Filter([](RVec<Electron>& electrons)
    {
        return ((electrons[0].electron_pt > 20e3 && electrons[1].electron_pt > 27e3) || (electrons[1].electron_pt > 20e3 && electrons[0].electron_pt > 27e3));
    }, {"di_electrons"});

    auto deltaR = leadingPt.Filter([] (RVec<Electron>& electrons)
    {
        return (DeltaR(electrons[0].Vector(), electrons[1].Vector()) > 0.01);
    }, {"di_electrons"});

    //new dataframe node: Contains a new column `dilep` in addition to the ones in `same_flavour` that stores the di-electron four-vector
    auto dilep = deltaR.Define("dilep",[] (RVec<Electron>& electrons)
    {
        return (electrons[0].Vector() + electrons[1].Vector());
    }, {"di_electrons"});

    //new dataframe node: contains only the events from `dilep` that have di-electron invariant mass between 81 and 101 GeV
    auto mass = dilep.Filter([] (PtEtaPhiEVector& dilep)
    {
        auto mass = dilep.M()/1e3;
        return ((mass >= 81) && (mass <= 101));
    }, {"dilep"});

    //new dataframe node: contains only the events from `mass` that have dilepton pT > 10 GeV
    auto ptCut = mass.Filter([] (PtEtaPhiEVector& dilep)
    {
        auto pT = dilep.Pt()/1e3;
        return pT > 10;
    }, {"dilep"});
    
    //end pre-selection -----------
    
    auto photon_passes_cuts = ptCut.Define("photons_pass_cuts",
    [&](RVec<Photon> photons)
    {
        photons.erase(std::remove_if(photons.begin(),photons.end(),
        [](Photon& x)
        {
            return ((std::abs(x.photon_eta) >= 2.37) || (std::abs(x.photon_eta) > 1.37 && std::abs(x.photon_eta) < 1.52));// || (!x.photon_id_loose));

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
            if (i==0 || ((std::abs(1-X) < std::abs(1-best_X)) and (delta_r < 1.5)))
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
    }, {"chosen_two"})
    .Define("higgs_mass_resolved",
    [&](RVec<Photon>& chosen_two, RVec<Electron>& electrons)
    {
        PtEtaPhiEVector gg = chosen_two[0].Vector()+chosen_two[1].Vector();
        PtEtaPhiEVector ll = electrons[0].Vector()+electrons[1].Vector();
        return (ll.M() + gg.M())/1e3;
    }, {"chosen_two", "di_electrons"})
    .Define("m_gg",
    [&](RVec<Photon>& chosen_two)
    {
        PtEtaPhiEVector gg = chosen_two[0].Vector()+chosen_two[1].Vector();
        return gg.M()/1e3;
    }, {"chosen_two"});
    
    auto SB = resolved.Filter(
    [&](double m_llgg)
    {
        return !((110 <= m_llgg) && (m_llgg <= 140));
        
    }, {"higgs_mass_resolved"});
    
    auto SR = resolved.Filter(
    [&](double m_llgg)
    {
        return ((110 <= m_llgg) && (m_llgg <= 140));
        
    }, {"higgs_mass_resolved"});
    
    auto di_ph_mm_SB = SB.Histo1D<double>({"data", "data", 10u, 0.6, 7.2}, "m_gg");
    auto di_ph_mm_SR = SR.Histo1D<double>({"data", "data", 10u, 0, 0}, "m_gg");
    
    TCanvas* c1 = new TCanvas();
    TLegend* legend = new TLegend(0.65, 0.4, 0.75, 0.6);
    di_ph_mm_SB->SetLineColor(kBlue);
    legend->AddEntry(&(*di_ph_mm_SB), di_ph_mm_SB->GetTitle(), "l");
    di_ph_mm_SB->SetTitle(";m_{#gamma#gamma} [GeV];Events");
    di_ph_mm_SB->GetYaxis()->CenterTitle(true);
    di_ph_mm_SB->GetXaxis()->SetTitleOffset(1.2);
    di_ph_mm_SB->Draw("HISTsame");
    gStyle->SetOptStat(0);
    TLatex Tl;
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.6, 0.8, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.6, 0.7,"#sqrt{s} = 13 TeV");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("Fig34A.pdf");
    
//    std::cout << SR.Count().GetValue() << '\n';
    c1 = new TCanvas();
    legend = new TLegend(0.65, 0.4, 0.75, 0.6);
    di_ph_mm_SR->SetLineColor(kBlue);
    legend->AddEntry(&(*di_ph_mm_SR), di_ph_mm_SR->GetTitle(), "l");
    di_ph_mm_SR->SetTitle(";m_{#gamma#gamma} [GeV];Events");
    di_ph_mm_SR->GetYaxis()->CenterTitle(true);
    di_ph_mm_SR->GetXaxis()->SetTitleOffset(1.2);
    di_ph_mm_SR->Draw("HISTsame");
    gStyle->SetOptStat(0);
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.6, 0.8, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.6, 0.7,"#sqrt{s} = 13 TeV");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("Fig34B.pdf");
    
}

//void Fig52()
//{
//    std::vector<std::vector<std::string>> input_filenames = {
//        //Z gamma background
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617070._000001.LGNTuple.root"},
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617064._000001.LGNTuple.root"},
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617074._000001.LGNTuple.root"},
//        //Data
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_data_test.root"},
//        //Jets
//        {"/Users/edwardfinkelstein/ATLAS_axion/Jets/Zee_lightJet_0-70.root"},
//        {"/Users/edwardfinkelstein/ATLAS_axion/Jets/Zee_lightJet_70-140.root"},
//        {"/Users/edwardfinkelstein/ATLAS_axion/Jets/Zee_lightJet_140-280.root"},
//        {"/Users/edwardfinkelstein/ATLAS_axion/Jets/Zee_cJet_0-70.root"},
//        {"/Users/edwardfinkelstein/ATLAS_axion/Jets/Zee_cJet_70-140.root"},
//        {"/Users/edwardfinkelstein/ATLAS_axion/Jets/Zee_cJet_140-280.root"},
//        {"/Users/edwardfinkelstein/ATLAS_axion/Jets/Zee_bJet_0-70.root"},
//        {"/Users/edwardfinkelstein/ATLAS_axion/Jets/Zee_bJet_70-140.root"},
//        {"/Users/edwardfinkelstein/ATLAS_axion/Jets/Zee_bJet_140-280.root"},
//    };
//
//    std::array<double,3> SFs = {((139e15)*(.871e-12)),((139e15)*(.199e-12)), ((139e15)*(.0345e-15))}; //numerators for Z-gamma bkg
//
//    std::array<double,9> JetNumeratorSFs = {((139e15)*(1.9828e-9)*(0.821204)),((139e15)*(110.64e-12)*(0.69275)),((139e15)*(40.645e-12)*(0.615906)),((139e15)*(1.9817e-9)*(0.1136684)),((139e15)*(110.47e-12)*(0.1912956)),((139e15)*(40.674e-12)*(0.2326772)),((139e15)*(1.9819e-9)*(0.0656969)),((139e15)*(110.53e-12)*(0.1158741)),((139e15)*(40.68e-12)*(0.1535215))}; //numerators for jet bkg
//
//    std::vector<const char*> prefixes = {"pty2_9_17", "pty_17_myy_0_80", "pty_17_myy_80", "data", "Zee_lightJet_0-70", "Zee_lightJet_70-140", "Zee_lightJet_140-280", "Zee_cJet_0-70", "Zee_cJet_70-140", "Zee_cJet_140-280", "Zee_bJet_0-70", "Zee_bJet_70-140", "Zee_bJet_140-280"};
//    std::vector<EColor> colors = {kBlue, kRed, kViolet, kGreen};
//    std::vector<EColor> Jetscolors = {kCyan, kOrange, kGreen, kYellow, kPink, kGray, kBlack, kSpring, kAzure};
//
//    std::vector<ROOT::RDF::RResultHandle> Nodes;
//    int count = 0;
//
//    for (auto& i: input_filenames)
//    {
//        SchottDataFrame df(MakeRDF(i, 8));
//
//        auto EventWeight = df.Define("EventWeight",
//        [](const RVec<float>& ei_event_weights_generator)
//        {
//            return  ((ei_event_weights_generator[0]) ? 1 / ei_event_weights_generator[0] : 1);
//
//        }, {"ei_event_weights_generator"});
//
//        auto two_leptons = df.Filter(
//        [](RVec<Muon>& muons, RVec<Electron> electrons)
//        {
//            electrons.erase(std::remove_if(electrons.begin(),electrons.end(),
//            [](Electron& ep)
//            {
//                return (!((ep.electron_pt/1e3 > 20) && (std::abs(ep.electron_eta) < 2.37) &&
//                          (!((1.37 < std::abs(ep.electron_eta)) && (std::abs(ep.electron_eta) < 1.52)))
//                          && (ep.electron_id_medium == 1)));
//
//            }), electrons.end());
//
//            return (electrons.size()==2 && muons.empty());
//
//        }, {"muons", "electrons"});
//
//        auto opp_charge = two_leptons.Define("di_electrons",
//        [](RVec<Electron> electrons)
//        {
//            electrons.erase(std::remove_if(electrons.begin(),electrons.end(),
//            [](Electron& ep)
//            {
//                return (!((ep.electron_pt/1e3 > 20) && (std::abs(ep.electron_eta) < 2.37) &&
//                (!((1.37 < std::abs(ep.electron_eta)) && (std::abs(ep.electron_eta) < 1.52)))
//                && (ep.electron_id_medium == 1)));
//
//            }), electrons.end());
//
//            return electrons;
//
//        },{"electrons"})
//        .Filter([](RVec<Electron> electrons)
//        {
//            return (electrons[0].electron_charge*electrons[1].electron_charge < 0);
//
//        }, {"di_electrons"});
//
//        auto leadingPt = opp_charge.Filter([](RVec<Electron>& electrons)
//        {
//            return ((electrons[0].electron_pt > 20e3 && electrons[1].electron_pt > 27e3) || (electrons[1].electron_pt > 20e3 && electrons[0].electron_pt > 27e3));
//        }, {"di_electrons"});
//
//        auto deltaR = leadingPt.Filter([] (RVec<Electron>& electrons)
//        {
//            return (DeltaR(electrons[0].Vector(), electrons[1].Vector()) > 0.01);
//        }, {"di_electrons"});
//
//        auto mass = deltaR.Filter([] (RVec<Electron>& electrons)
//        {
//            auto mass = (electrons[0].Vector() + electrons[1].Vector()).M()/1e3;
//            return ((mass >= 81) && (mass <= 101));
//        }, {"di_electrons"});
//
//        auto ptCut = mass.Filter([] (RVec<Electron>& electrons)
//        {
//            auto pT = (electrons[0].Vector() + electrons[1].Vector()).Pt()/1e3;
//            return pT > 10;
//        }, {"di_electrons"});
//
//        auto photon_passes_cuts = ptCut.Define("photons_pass_cuts",
//        [&](RVec<Photon> photons)
//        {
//            photons.erase(std::remove_if(photons.begin(),photons.end(),
//            [](Photon& x)
//            {
//              return ((std::abs(x.photon_eta) >= 2.37) || (std::abs(x.photon_eta) > 1.37 && std::abs(x.photon_eta) < 1.52) || (!x.photon_id_loose));
//
//            }), photons.end());
//
//            return photons;
//        }, {"photons"});
//
//        auto merged_reco_photons_matched = photon_passes_cuts.Filter(
//        [&](RVec<Photon>& reco_photons_test)
//        {
//            RVec<Photon> reco_photons_matched = reco_photons_test;
//            if (reco_photons_matched.size() == 1)
//            {
//                return reco_photons_matched[0].photon_pt > 20e3;
//            }
//            else if (reco_photons_matched.empty())
//            {
//                return false;
//            }
//
//            auto combs = Combinations(reco_photons_matched, 2);
//            size_t length = combs[0].size();
//            double delta_r, m, pt, X, best_X, pt1, pt2, chosen_delta_r;
//
//            for (size_t i=0; i<length; i++)
//            {
//                delta_r = DeltaR(reco_photons_matched[combs[0][i]].Vector(), reco_photons_matched[combs[1][i]].Vector());
//                m = (reco_photons_matched[combs[0][i]].Vector() + reco_photons_matched[combs[1][i]].Vector()).M();
//                pt = (reco_photons_matched[combs[0][i]].Vector() + reco_photons_matched[combs[1][i]].Vector()).Pt();
//                X = delta_r*(pt/(2.0*m));
//                if (i==0 || std::abs(1-X) < std::abs(1-best_X))
//                {
//                    best_X = X;
//                    pt1 = reco_photons_matched[combs[0][i]].photon_pt;
//                    pt2 = reco_photons_matched[combs[1][i]].photon_pt;
//                    chosen_delta_r = delta_r;
//                }
//            }
//            if (pt1 > 10e3 && pt2 > 10e3 && best_X > 0.96 && best_X < 1.2 && chosen_delta_r < 1.5)
//            {
//                return false;
//            }
//
//            for (auto& p: reco_photons_matched)
//            {
//                if (p.photon_pt > 20e3)
//                {
//                    return true;
//                }
//            }
//            return false;
//
//        }, {"photons_pass_cuts"})
//        .Define("merged_photon",
//        [](RVec<Photon>& reco_photons_matched)
//        {
//            for (auto& p: reco_photons_matched)
//            {
//                if (p.photon_pt > 20e3)
//                {
//                    return p;
//                }
//            }
//
//            return reco_photons_matched[0]; //jic the compiler complains
//
//        }, {"photons_pass_cuts"}).Define("merged_photon_pt",
//        [](Photon& p)
//        {
//            return p.photon_pt/1e3;
//        }, {"merged_photon"}).Define("reconstructed_mass",
//        [&](RVec<Electron>& di_electrons, Photon& merged_photon)
//        {
//            auto four_momentum = di_electrons[0].Vector() + di_electrons[1].Vector();
//
//            return (four_momentum + merged_photon.Vector()).M()/1e3;
//
//        }, {"di_electrons", "merged_photon"}).Filter(
//        [](double reconstructed_mass, RVec<float>& Eratio)
//        {
//            return ((!Any(Eratio <= 0.8)) && ((reconstructed_mass < 110) || (reconstructed_mass > 130)));
//        }, {"reconstructed_mass", "photon_shower_shape_e_ratio"}).Define("totEventWeight", [](RVec<float> photon_id_eff, RVec<float> photon_iso_eff, RVec<float> photon_trg_eff/*, RVec<float> ei_event_weights_generator*/)
//        {
//            auto ResizeVal = std::max({photon_id_eff.size(), photon_iso_eff.size(), photon_trg_eff.size()});
//            photon_id_eff.resize(ResizeVal,1);
//            photon_iso_eff.resize(ResizeVal,1);
//            photon_trg_eff.resize(ResizeVal,1);
//
//            return photon_id_eff*photon_iso_eff*photon_trg_eff;//*ei_event_weights_generator[0];
//
//        }, {"photon_id_eff", "photon_iso_eff", "photon_trg_eff",/* "ei_event_weights_generator"*/});
//
//        if (count <= 2 || count >= 4)
//        {
//            Nodes.push_back(merged_reco_photons_matched.Sum<RVec<float>>("totEventWeight"));
//            Nodes.push_back(EventWeight.Sum<float>("EventWeight"));
//        }
//
//        Nodes.push_back(merged_reco_photons_matched.Histo1D<double>({prefixes[count], prefixes[count++], 60u, 0, 165}, "merged_photon_pt", "totEventWeight"));
//    }
//
////    0   1   2   Z-gamma
////    3   4   5   Z-gamma
////    6   7   8   Z-gamma
////            9   Data
////    10  11  12  Z-jets
////    13  14  15  Z-jets
////    16  17  18  Z-jets
////    19  20  21  Z-jets
////    22  23  24  Z-jets
////    25  26  27  Z-jets
////    28  29  30  Z-jets
////    31  32  33  Z-jets
////    34  35  36  Z-jets
//
//    ROOT::RDF::RunGraphs(Nodes); // running all computation nodes concurrently
//
//    double factor;
//    count = 0;
//    for (auto& i: {0,3,6}) //Z-gamma
//    {
//        factor += (*Nodes[i].GetResultPtr<float>())*(SFs[count++] / *Nodes[i+1].GetResultPtr<float>());
//    }
//
//    for (int i = 10, j = 0; (i <= 34 && j <= 8); i += 3, j++) //Z-jets
//    {
//        factor += (*Nodes[i].GetResultPtr<float>())*(JetNumeratorSFs[j] / *Nodes[i+1].GetResultPtr<float>());
//    }
//
//    auto hs = new THStack("hs3","");
//    TCanvas* c1 = new TCanvas();
//    TLegend* legend = new TLegend(0.55, 0.2, 0.85, 0.6);
//    count=0;
//    for (auto& i: {2,5,8,9}) //Z-gamma
//    {
//        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0 && i != 9)
//        {
//            Nodes[i].GetResultPtr<TH1D>()->Scale(SFs[count] / *Nodes[i-1].GetResultPtr<float>());
//        }
//        Nodes[i].GetResultPtr<TH1D>()->SetFillColor(colors[count++]);
//        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "f");
//
//        if (i != 9)
//        {
//            hs->Add(&*Nodes[i].GetResultPtr<TH1D>());
//        }
//    }
//
//    for (int i = 12, j = 0; (i <= 36 && j <= 8); i += 3, j++) //Z-jets
//    {
//        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0)
//        {
//            Nodes[i].GetResultPtr<TH1D>()->Scale(JetNumeratorSFs[j] / *Nodes[i-1].GetResultPtr<float>());
//        }
//        Nodes[i].GetResultPtr<TH1D>()->SetFillColor(Jetscolors[j]);
//        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "f");
//        hs->Add(&*Nodes[i].GetResultPtr<TH1D>());
//    }
//
//    hs->Draw("HIST");
//    Nodes[9].GetResultPtr<TH1D>()->Draw("HISTsame");
//    hs->SetTitle(";photon p_{T} [GeV];Events");
//    hs->GetYaxis()->CenterTitle(true);
//    hs->GetXaxis()->SetTitleOffset(1.2);
//    hs->GetYaxis()->SetTitleOffset(1.35);
//    gStyle->SetOptStat(0);
//    TLatex Tl;
//    Tl.SetTextSize(0.03);
//    Tl.DrawLatexNDC(0.6, 0.8, "#it{ATLAS} Internal");
//    Tl.DrawLatexNDC(0.6, 0.7,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
//    legend->SetBorderSize(0);
//    legend->Draw("same");
//    c1->SaveAs("Fig52A.pdf");
//}

void Figs_34_52()
{
    auto start_time = Clock::now();
    Fig34();
//    Fig52();
    auto end_time = Clock::now();
    std::cout << "Time difference: "
       << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count()/1e9 << " nanoseconds" << std::endl;
}

int main()
{
    Figs_34_52();
}
