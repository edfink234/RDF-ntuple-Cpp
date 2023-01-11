#include <iostream>
#include <vector>
#include <chrono>
#include <algorithm>
#include <unordered_map>
#include <iomanip>
#include <array>
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

void fig27()
{
    auto hs = new THStack("hs","");
    std::vector<std::string> input_filenames = {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617070._000001.LGNTuple.root", "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617064._000001.LGNTuple.root",
        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_data_test.root"
    }; //"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617074._000001.LGNTuple.root"};
    
    std::vector<const char*> prefixes = {"pty2_9_17", "pty_17_myy_0_80", "data"};
//    "pty_17_myy_80"};
    std::vector<EColor> colors = {kBlue, kRed, kGreen};
    std::vector<ROOT::RDF::RResultPtr<TH1D>> histos;
    histos.reserve(3);
    int count = 0;
    for (auto& i: input_filenames)
    {
        SchottDataFrame df(MakeRDF({i}, 8));
        
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
        
        auto diphotons = photon_passes_cuts.Define("chosen_two",
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
                if ((delta_r < 1.5) && (X > 0.96) && (X < 1.2))
                {
                    x = {reco_photons_matched[combs[0][i]], reco_photons_matched[combs[1][i]]};
                    return x;
                }
            }
            return x;
        }, {"photons_pass_cuts"}).Filter(
        [&](RVec<Photon>& reco_photons_matched)
        {
            return (reco_photons_matched.size()==2);
        }, {"chosen_two"}).Define("mass",
        [&](RVec<Photon>& diph)
        {
          return (diph[0].Vector()+diph[1].Vector()).M()/1e3;
        }, {"chosen_two"});
        
//        std::cout << *(diphotons.Min<double>("mass")) << '\n';
//        std::cout << *(diphotons.Count()) << '\n';
        
        histos.push_back(diphotons.Histo1D<double>({prefixes[count], prefixes[count++], 100u, 0, 95}, "mass"));
        
        auto passed = diphotons.Count();
        std::cout << *passed << '\n';
    }
    
    
    double factor;
    TCanvas* c1 = new TCanvas();
    TLegend* legend = new TLegend(0.65, 0.4, 0.85, 0.6);
    count = 0;
    for (auto& h: histos)
    {
        h->SetFillColor(colors[count++]);
        legend->AddEntry(&(*h), h->GetTitle(), "f");
    
//            factor = h->Integral();
//            h->Scale(factor/h->Integral());
//            h->SetAxisRange(0., 82,"Y");
//            h->SetAxisRange(0., 500,"Y");
        if (&h != &histos.back())
            hs->Add(&*h);
//            h->Draw("HIST");
            
    }

    hs->Draw("HIST");
    hs->SetTitle(";m_{#gamma#gamma} [GeV];Events");
    hs->GetYaxis()->CenterTitle(true);
    hs->GetXaxis()->SetTitleOffset(1.2);
    histos[2]->Draw("same");
    gStyle->SetOptStat(0);
    TLatex Tl;
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.6, 0.8, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.6, 0.7,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("Fig27.png");

}

void fig28()
{
    std::vector<std::string> input_filenames = {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617070._000001.LGNTuple.root", "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617064._000001.LGNTuple.root",
        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_data_test.root"
    }; //"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617074._000001.LGNTuple.root"};
    
    std::vector<const char*> prefixes = {"pty2_9_17", "pty_17_myy_0_80", "data"};
//    "pty_17_myy_80"};
    std::vector<EColor> colors = {kBlue, kRed, kGreen};
    std::vector<ROOT::RDF::RResultPtr<TH1D>> histos_mass;
    histos_mass.reserve(3);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> histos_deltaR;
    histos_deltaR.reserve(3);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> histos_deltaPhi;
    histos_deltaPhi.reserve(3);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> histos_deltaEta;
    histos_deltaEta.reserve(3);
    int count = 0;
    for (auto& i: input_filenames)
    {
        SchottDataFrame df(MakeRDF({i}, 8));
        
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
        
        auto diphotons = photon_passes_cuts.Define("chosen_two",
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
                if ((delta_r < 1.5) && (X > 0.96) && (X < 1.2))
                {
                    x = {reco_photons_matched[combs[0][i]], reco_photons_matched[combs[1][i]]};
                    return x;
                }
            }
            return x;
        }, {"photons_pass_cuts"}).Filter(
        [&](RVec<Photon>& reco_photons_matched)
        {
            if (reco_photons_matched.size()!=2)
            {
                return false;
            }
            
            PtEtaPhiEVector rpm = reco_photons_matched[0].Vector()+reco_photons_matched[1].Vector();
            
            return rpm.M()/1e3 < 20;

        }, {"chosen_two"}).Define("mass",
        [&](RVec<Photon>& diph)
        {
            return (diph[0].Vector()+diph[1].Vector()).M()/1e3;
        }, {"chosen_two"}).Define("deltaR",
        [&](RVec<Photon>& diph)
        {
            return DeltaR(diph[0].Vector(),diph[1].Vector());
        }, {"chosen_two"}).Define("deltaPhi",
        [&](RVec<Photon>& diph)
        {
            return ROOT::Math::VectorUtil::DeltaPhi(diph[0].Vector(),diph[1].Vector());
        }, {"chosen_two"}).Define("deltaEta",
        [&](RVec<Photon>& diph)
        {
            return abs((diph[0].Vector() - diph[1].Vector()).Eta());
        }, {"chosen_two"});
        
//        std::cout << *(diphotons.Min<double>("mass")) << '\n';
//        std::cout << *(diphotons.Count()) << '\n';
        histos_mass.push_back(diphotons.Histo1D<double>({prefixes[count], prefixes[count], 100u, 0, 95}, "mass"));
        histos_deltaR.push_back(diphotons.Histo1D<double>({prefixes[count], prefixes[count], 100u, 0, 2.5}, "deltaR"));
        histos_deltaPhi.push_back(diphotons.Histo1D<double>({prefixes[count], prefixes[count], 100u, 0, 2}, "deltaPhi"));
        histos_deltaEta.push_back(diphotons.Histo1D<double>({prefixes[count], prefixes[count++], 100u, 0, 4}, "deltaEta"));
        
        auto passed = diphotons.Count();
        std::cout << *passed << '\n';
    }
    
    double factor;
    auto hs = new THStack("hs1","");
    TCanvas* c1 = new TCanvas();
    TLegend* legend = new TLegend(0.65, 0.4, 0.85, 0.6);
    count = 0;
    for (auto& h: histos_mass)
    {
        h->SetFillColor(colors[count++]);
        legend->AddEntry(&(*h), h->GetTitle(), "f");
        if (&h != &histos_mass.back())
            hs->Add(&*h);
    }

    hs->Draw("HIST");
    hs->SetTitle(";m_{#gamma#gamma} [GeV];Events");
    hs->GetYaxis()->CenterTitle(true);
    hs->GetXaxis()->SetTitleOffset(1.2);
    histos_mass[2]->Draw("same");
    gStyle->SetOptStat(0);
    TLatex Tl;
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.6, 0.8, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.6, 0.7,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("Fig28A.png");
    
    hs = new THStack("hs2","");
    c1 = new TCanvas();
    legend = new TLegend(0.65, 0.4, 0.85, 0.6);
    count = 0;
    for (auto& h: histos_deltaR)
    {
        h->SetFillColor(colors[count++]);
        legend->AddEntry(&(*h), h->GetTitle(), "f");
        if (&h != &histos_deltaR.back())
            hs->Add(&*h);
    }

    hs->Draw("HIST");
    hs->SetTitle(";#DeltaR_{#gamma#gamma};Events");
    hs->GetYaxis()->CenterTitle(true);
    hs->GetXaxis()->SetTitleOffset(1.2);
    histos_deltaR[2]->Draw("same");
    gStyle->SetOptStat(0);
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.6, 0.8, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.6, 0.7,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("Fig28B.png");
    
    hs = new THStack("hs3","");
    c1 = new TCanvas();
    legend = new TLegend(0.65, 0.4, 0.85, 0.6);
    count = 0;
    for (auto& h: histos_deltaPhi)
    {
        h->SetFillColor(colors[count++]);
        legend->AddEntry(&(*h), h->GetTitle(), "f");
        if (&h != &histos_deltaPhi.back())
            hs->Add(&*h);
    }

    hs->Draw("HIST");
    hs->SetTitle(";#Delta#phi_{#gamma#gamma};Events");
    hs->GetYaxis()->CenterTitle(true);
    hs->GetXaxis()->SetTitleOffset(1.2);
    histos_deltaPhi[2]->Draw("same");
    gStyle->SetOptStat(0);
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.6, 0.8, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.6, 0.7,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("Fig28C.png");
    
    hs = new THStack("hs4","");
    c1 = new TCanvas();
    legend = new TLegend(0.65, 0.4, 0.85, 0.6);
    count = 0;
    for (auto& h: histos_deltaEta)
    {
        h->SetFillColor(colors[count++]);
        legend->AddEntry(&(*h), h->GetTitle(), "f");
        if (&h != &histos_deltaEta.back())
            hs->Add(&*h);
    }

    hs->Draw("HIST");
    hs->SetTitle(";#Delta#eta_{#gamma#gamma};Events");
    hs->GetYaxis()->CenterTitle(true);
    hs->GetXaxis()->SetTitleOffset(1.2);
    histos_deltaEta[2]->Draw("same");
    gStyle->SetOptStat(0);
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.6, 0.8, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.6, 0.7,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("Fig28D.png");
}

void fig41()
{
    std::vector<std::string> input_filenames =
    {
        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root",
        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root", "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617070._000001.LGNTuple.root", "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617064._000001.LGNTuple.root",
        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_data_test.root",
      //  "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617074._000001.LGNTuple.root"
    };
    
    std::vector<const char*> prefixes = {"Sig m_{A} = 5 GeV", "Sig m_{A} = 1 GeV", "pty2_9_17", "pty_17_myy_0_80", "data"};//, "pty_17_myy_80"};
    std::vector<EColor> colors = {kBlack, kMagenta, kBlue, kRed, kGreen};//, kViolet};
    
    std::vector<ROOT::RDF::RResultPtr<TH1D>> histos_sig;
    histos_sig.reserve(2);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> histos_back_plus_data;
    histos_back_plus_data.reserve(3);//4);
    
    int count = 0;
    for (auto& i: input_filenames)
    {
        SchottDataFrame df(MakeRDF({i}, 8));
        
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
            
            reco_photons_matched.erase(std::remove_if(reco_photons_matched.begin(),reco_photons_matched.end(),
            [](Photon& x)
            {
                return x.photon_pt <= 10e3;

            }), reco_photons_matched.end());
            
            if (reco_photons_matched.size() == 1)
            {
                return true ? reco_photons_matched[0].photon_pt > 20e3 : false;
            }
            
            else if (reco_photons_matched.empty())
            {
                return false;
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
                if ((delta_r < 1.5) && (X > 0.96) && (X < 1.2))
                {
                    return false;
                }
            }
            
            for (auto& p: reco_photons_matched)
            {
                if (p.photon_pt > 20e3)
                {
                    return true;
                }
            }
            return false;
            
        }, {"photons_pass_cuts"})
        .Define("merged_photon",
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
            
        }, {"photons_pass_cuts"});
        
        auto dilepton_and_photon = merged_reco_photons_matched
        .Define("reconstructed_mass",[&](RVec<Electron>& di_electrons, Photon& merged_photon)
        {
            auto four_momentum = di_electrons[0].Vector() + di_electrons[1].Vector();
            
            return (four_momentum + merged_photon.Vector()).M()/1e3;
            
        }, {"di_electrons", "merged_photon"});
        
//        std::cout << *(diphotons.Min<double>("mass")) << '\n';
//        std::cout << *(diphotons.Count()) << '\n';
        
        if (count < 2) //signal only
        {
            histos_sig.push_back(dilepton_and_photon.Histo1D<double>({prefixes[count], prefixes[count++], 100u, 80, 200}, "reconstructed_mass"));
        }
        else
        {
            histos_back_plus_data.push_back(dilepton_and_photon.Histo1D<double>({prefixes[count], prefixes[count++], 100u, 80, 410}, "reconstructed_mass"));
        }

        auto passed = dilepton_and_photon.Count();
        std::cout << *passed << '\n';
    }
    count = 0;
    TCanvas* c1 = new TCanvas();
    TLegend* legend = new TLegend(0.65, 0.4, 0.85, 0.6);
    double factor;
    for (auto& h: histos_sig)
    {
        h->SetLineColor(colors[count++]);
        legend->AddEntry(&(*h), h->GetTitle(), "l");
        
        if (&h == &histos_sig.front())
        {
            factor = h->Integral();
            h->Scale(factor/h->Integral());
            h->SetTitle(";m_{ll#gamma} [GeV];Events");
            h->GetYaxis()->CenterTitle(true);
            h->GetXaxis()->SetTitleOffset(1.2);
            h->SetAxisRange(0., 200, "Y");
//            h->SetAxisRange(0., 1200,"Y");
            h->Draw("HIST");
        }
        else
        {
            h->Scale(factor/h->Integral());
            h->Draw("HISTsame");
            gPad->Modified(); gPad->Update();
        }
    }
    
    gStyle->SetOptStat(0);
    TLatex Tl;
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.6, 0.8, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.6, 0.7,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("Fig41B.png");
    
//    count = 0;
    auto hs = new THStack("hs3","");
    c1 = new TCanvas();
    legend = new TLegend(0.65, 0.4, 0.85, 0.6);
    for (auto& h: histos_back_plus_data)
    {
        h->SetFillColor(colors[count++]);
        legend->AddEntry(&(*h), h->GetTitle(), "f");
//        if (&h != &histos_back_plus_data.back())
        if (strcmp(h->GetTitle(),"data")!=0)
            hs->Add(&*h);
    }

    hs->Draw("HIST");
    histos_back_plus_data[2]->Draw("HISTsame");
    count=0;
    for (auto &h: histos_sig)
    {
        h->SetLineWidth(2);
        h->SetLineColor(colors[count]);
        std::cout << prefixes[count] << '\n';
        legend->AddEntry(&(*h), prefixes[count++], "l");
        
        factor = h->Integral();
        h->Scale(factor/h->Integral());
        h->SetTitle(";m_{ll#gamma} [GeV];Events");
        h->GetYaxis()->CenterTitle(true);
        h->GetXaxis()->SetTitleOffset(1.2);
//            h->SetAxisRange(0., 1200,"Y");
        h->DrawClone("HISTsame");
        gPad->Modified();
        gPad->Update();
        c1->Modified();
        c1->Update();
    }
    hs->SetMaximum(200);
    hs->SetTitle(";m_{ll#gamma} [GeV];Events");
    hs->GetYaxis()->CenterTitle(true);
    hs->GetXaxis()->SetTitleOffset(1.2);

//    histos_back_plus_data[2]->Draw("same");
    gStyle->SetOptStat(0);
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.6, 0.8, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.6, 0.7,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
    legend->SetBorderSize(0);
    legend->Draw("same");
    c1->SaveAs("Fig41A.png");
}

void fig48()
{
    std::vector<std::string> input_filenames =
//    {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617070._000001.LGNTuple.root", "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617064._000001.LGNTuple.root","/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_data_test.root","/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617074._000001.LGNTuple.root"};
    { "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617070._000001.LGNTuple.root", "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617064._000001.LGNTuple.root",
        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_data_test.root",
//        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617074._000001.LGNTuple.root"
    };
    
    std::vector<const char*> prefixes = {"pty2_9_17", "pty_17_myy_0_80", "data"};//, "pty_17_myy_80"};
    std::vector<EColor> colors = {kBlue, kRed, kGreen};//, kViolet};
    
    std::vector<ROOT::RDF::RResultPtr<TH1D>> histos_back_plus_data_mass;
    histos_back_plus_data_mass.reserve(3);//4);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> histos_back_plus_data_pt;
    histos_back_plus_data_pt.reserve(3);//4);
    
    int count = 0;
    for (auto& i: input_filenames)
    {
        SchottDataFrame df(MakeRDF({i}, 8));
        
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
            
            reco_photons_matched.erase(std::remove_if(reco_photons_matched.begin(),reco_photons_matched.end(),
            [](Photon& x)
            {
                return x.photon_pt <= 10e3;

            }), reco_photons_matched.end());
            
            if (reco_photons_matched.size() == 1)
            {
                return true ? reco_photons_matched[0].photon_pt > 20e3 : false;
            }
            
            else if (reco_photons_matched.empty())
            {
                return false;
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
                if ((delta_r < 1.5) && (X > 0.96) && (X < 1.2))
                {
                    return false;
                }
            }
            
            for (auto& p: reco_photons_matched)
            {
                if (p.photon_pt > 20e3)
                {
                    return true;
                }
            }
            return false;
            
        }, {"photons_pass_cuts"})
        .Define("merged_photon",
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
            
        }, {"photons_pass_cuts"});
        
        auto dilepton_and_photon = merged_reco_photons_matched
        .Define("dilepton_mass",[&](RVec<Electron>& di_electrons)
        {
            auto four_momentum = di_electrons[0].Vector() + di_electrons[1].Vector();
            
            return four_momentum.M()/1e3;
            
        }, {"di_electrons"})
        .Define("merged_photon_pt",[&](Photon& merged_photon)
        {
            return merged_photon.Vector().Pt()/1e3;
            
        }, {"merged_photon"});
        
//        std::cout << *(diphotons.Min<double>("mass")) << '\n';
//        std::cout << *(diphotons.Count()) << '\n';
        
        histos_back_plus_data_mass.push_back(dilepton_and_photon.Histo1D<double>({prefixes[count], prefixes[count], 100u, 60, 120}, "dilepton_mass"));
        histos_back_plus_data_pt.push_back(dilepton_and_photon.Histo1D<double>({prefixes[count], prefixes[count++], 100u, 0, 165}, "merged_photon_pt"));

        auto passed = dilepton_and_photon.Count();
        std::cout << *passed << '\n';
    }
    count = 0;
    TCanvas* c1 = new TCanvas();
    TLegend* legend = new TLegend(0.65, 0.4, 0.85, 0.6);
    auto hs = new THStack("hs1","");
    double factor;
    
    for (auto& h: histos_back_plus_data_mass)
    {
        h->SetFillColor(colors[count++]);
        legend->AddEntry(&(*h), h->GetTitle(), "f");
        if (strcmp(h->GetTitle(),"data")!=0)
            hs->Add(&*h);
    }
    hs->Draw("HIST");
    histos_back_plus_data_mass[2]->Draw("HISTsame");
    hs->SetTitle(";m_{ll} [GeV];Events");
    hs->GetYaxis()->CenterTitle(true);
    hs->GetXaxis()->SetTitleOffset(1.2);
    
    gStyle->SetOptStat(0);
    TLatex Tl;
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.6, 0.8, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.6, 0.7,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("Fig48A.png");
    
    count = 0;
    hs = new THStack("hs2","");
    c1 = new TCanvas();
    legend = new TLegend(0.65, 0.4, 0.85, 0.6);
    for (auto& h: histos_back_plus_data_pt)
    {
        h->SetFillColor(colors[count++]);
        legend->AddEntry(&(*h), h->GetTitle(), "f");
        if (count!=3)
            hs->Add(&*h);
    }
    hs->Draw("HIST");
    histos_back_plus_data_pt[2]->Draw("HISTsame");
    hs->SetTitle(";p_{t} [GeV];Events");
    hs->GetYaxis()->CenterTitle(true);
    hs->GetXaxis()->SetTitleOffset(1.2);
    
    gStyle->SetOptStat(0);
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.6, 0.8, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.6, 0.7,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("Fig48B.png");
}

void fig59()
{
    std::vector<std::string> input_filenames =
//    {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617070._000001.LGNTuple.root", "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617064._000001.LGNTuple.root","/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_data_test.root","/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617074._000001.LGNTuple.root"};
    {
        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root",
        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root", "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617070._000001.LGNTuple.root", "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617064._000001.LGNTuple.root",
        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_data_test.root",
      //  "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617074._000001.LGNTuple.root"
    };
    
    std::vector<const char*> prefixes = {"Sig m_{A} = 5 GeV", "Sig m_{A} = 1 GeV", "pty2_9_17", "pty_17_myy_0_80", "data"};//, "pty_17_myy_80"};
    std::vector<EColor> colors = {kBlack, kMagenta, kBlue, kRed, kGreen};//, kViolet};
    
    std::vector<ROOT::RDF::RResultPtr<TH1D>> histos_sig;
    histos_sig.reserve(2);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> histos_back_plus_data;
    histos_back_plus_data.reserve(3);//4);
    
    int count = 0;
    for (auto& i: input_filenames)
    {
        SchottDataFrame df(MakeRDF({i}, 8));
        
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
            
            reco_photons_matched.erase(std::remove_if(reco_photons_matched.begin(),reco_photons_matched.end(),
            [](Photon& x)
            {
                return x.photon_pt <= 10e3;

            }), reco_photons_matched.end());
            
            if (reco_photons_matched.size() == 1)
            {
                return true ? reco_photons_matched[0].photon_pt > 20e3 : false;
            }
            
            else if (reco_photons_matched.empty())
            {
                return false;
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
                if ((delta_r < 1.5) && (X > 0.96) && (X < 1.2))
                {
                    return false;
                }
            }
            
            for (auto& p: reco_photons_matched)
            {
                if (p.photon_pt > 20e3)
                {
                    return true;
                }
            }
            return false;
            
        }, {"photons_pass_cuts"})
        .Define("merged_photon",
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
            
        }, {"photons_pass_cuts"});
        
        auto dilepton_and_photon = merged_reco_photons_matched
        .Define("dilepton_mergedPhoton_deltaR",[&](RVec<Electron>& di_electrons, Photon& merged_photon)
        {
            auto dilepton_four_momentum = di_electrons[0].Vector() + di_electrons[1].Vector();
            return DeltaR(dilepton_four_momentum, merged_photon.Vector());
            
        }, {"di_electrons", "merged_photon"});
        
//        std::cout << *(diphotons.Min<double>("mass")) << '\n';
//        std::cout << *(diphotons.Count()) << '\n';
        
        if (count < 2) //signal only
        {
            histos_sig.push_back(dilepton_and_photon.Histo1D<double>({prefixes[count], prefixes[count++], 100u, 0, 4.7}, "dilepton_mergedPhoton_deltaR"));
        }
        else
        {
            histos_back_plus_data.push_back(dilepton_and_photon.Histo1D<double>({prefixes[count], prefixes[count++], 100u, 0, 4.7}, "dilepton_mergedPhoton_deltaR"));
        }
        auto passed = dilepton_and_photon.Count();
        std::cout << *passed << '\n';
    }
    count = 0;
    TCanvas* c1 = new TCanvas();
    TLegend* legend = new TLegend(0.2, 0.4, 0.4, 0.6);
    double factor;
    for (auto& h: histos_sig)
    {
        h->SetLineColor(colors[count++]);
        legend->AddEntry(&(*h), h->GetTitle(), "l");
        
        if (&h == &histos_sig.front())
        {
            factor = h->Integral();
            h->Scale(factor/h->Integral());
            h->SetTitle(";#DeltaR (ll#gamma);Events");
            h->GetYaxis()->CenterTitle(true);
            h->GetXaxis()->SetTitleOffset(1.2);
            h->SetAxisRange(0., 68, "Y");
//            h->SetAxisRange(0., 1200,"Y");
            h->Draw("HIST");
        }
        else
        {
            h->Scale(factor/h->Integral());
            h->Draw("HISTsame");
            gPad->Modified(); gPad->Update();
        }
    }
    
    gStyle->SetOptStat(0);
    TLatex Tl;
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.2, 0.8, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.2, 0.7,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("Fig59D.png");
    
    auto hs = new THStack("hs1","");
    c1 = new TCanvas();
    legend = new TLegend(0.2, 0.4, 0.4, 0.6);
    for (auto& h: histos_back_plus_data)
    {
        h->SetFillColor(colors[count++]);
        legend->AddEntry(&(*h), h->GetTitle(), "f");
        if (strcmp(h->GetTitle(),"data")!=0)
            hs->Add(&*h);
    }
    hs->Draw("HIST");
    histos_back_plus_data[2]->Draw("HISTsame");
    hs->SetTitle(";#DeltaR (ll#gamma);Events");
    hs->GetYaxis()->CenterTitle(true);
    hs->GetXaxis()->SetTitleOffset(1.2);
    
    gStyle->SetOptStat(0);
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.2, 0.8, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.2, 0.7,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("Fig59A.png");
    count=2;
    c1 = new TCanvas();
    legend = new TLegend(0.2, 0.4, 0.4, 0.6);
    hs = new THStack("hs2","");
    for (auto& h: histos_back_plus_data)
    {
        if (strcmp(h->GetTitle(),"data")!=0)
        {
            h->SetFillColor(colors[count++]);
            legend->AddEntry(&(*h), h->GetTitle(), "f");
            hs->Add(&*h);
        }
    }
    hs->Draw("HIST");
    count=0;
    for (auto& h: histos_sig)
    {
        h->SetLineColor(colors[count++]);
        h->SetLineWidth(2);
        legend->AddEntry(&(*h), prefixes[count-1], "l");
        
        h->SetTitle(";#DeltaR (ll#gamma);Events");
        h->GetYaxis()->CenterTitle(true);
        h->GetXaxis()->SetTitleOffset(1.2);
        h->SetAxisRange(0., 68, "Y");
//            h->SetAxisRange(0., 1200,"Y");
        h->Draw("HISTsame");
        gPad->Modified(); gPad->Update();
    }
    hs->SetTitle(";#DeltaR (ll#gamma);Events");
    hs->GetYaxis()->CenterTitle(true);
    hs->GetXaxis()->SetTitleOffset(1.2);
    
    gStyle->SetOptStat(0);
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.2, 0.8, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.2, 0.7,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("Fig59B.png");
    
    count=2;
    c1 = new TCanvas();
    legend = new TLegend(0.2, 0.4, 0.4, 0.6);
    count=2;
    factor = histos_sig.front()->Integral();
    hs = new THStack("hs3","");
    for (auto& h: histos_back_plus_data)
    {
        if (strcmp(h->GetTitle(),"data")!=0)
        {
            h->Scale(factor/(2*h->Integral()));
            h->SetFillColor(colors[count++]);
            legend->AddEntry(&(*h), h->GetTitle(), "f");
            hs->Add(&*h);
        }
    }
    hs->Draw("HIST");
    count=0;
    for (auto& h: histos_sig)
    {
        h->SetLineColor(colors[count++]);
        h->SetLineWidth(2);
        legend->AddEntry(&(*h), prefixes[count-1], "l");
        
        h->Scale(factor/h->Integral());
        h->SetTitle(";#DeltaR (ll#gamma);Events");
        h->GetYaxis()->CenterTitle(true);
        h->GetXaxis()->SetTitleOffset(1.2);
        h->SetAxisRange(0., 68, "Y");
//            h->SetAxisRange(0., 1200,"Y");
        h->Draw("HISTsame");
        gPad->Modified(); gPad->Update();
    }
   
    
    hs->SetTitle(";#DeltaR (ll#gamma);Events");
    hs->GetYaxis()->CenterTitle(true);
    hs->GetXaxis()->SetTitleOffset(1.2);
    
    gStyle->SetOptStat(0);
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.2, 0.8, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.2, 0.7,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("Fig59C.png");
}

void Table9()
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
    [&](RVec<TruthParticle> truth_particles)
    {
        truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
        [](TruthParticle& x)
        {
            return ((abs(x.mc_pdg_id) != 22) || (abs(x.mc_eta) >= 2.37) || (x.mc_pt <= 10e3));

        }), truth_particles.end());
        
        return truth_particles;
        
    }, {"truth_particles"});
    
    auto diphotons = photon_passes_cuts.Define("chosen_two",
    [](RVec<TruthParticle>& reco_photons_matched)
    {
        RVec<TruthParticle> x;
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
            if ((delta_r < 1.5) && (X > 0.96) && (X < 1.2))
            {
                x = {reco_photons_matched[combs[0][i]], reco_photons_matched[combs[1][i]]};
                return x;
            }
        }
        return x;
    }, {"photons_pass_cuts"}).Filter(
    [&](RVec<TruthParticle>& reco_photons_matched)
    {
        if (reco_photons_matched.size()!=2)
        {
            return false;
        }
        return true;
    }, {"chosen_two"}).Define("leading_photon",
    [&](RVec<TruthParticle>& reco_photons_matched)
    {
        return (reco_photons_matched[0].mc_pt > reco_photons_matched[1].mc_pt)
                ? reco_photons_matched[0]
                : reco_photons_matched[1];
    },{"chosen_two"}).Define("subleading_photon",
    [&](RVec<TruthParticle>& reco_photons_matched)
    {
        return (reco_photons_matched[0].mc_pt < reco_photons_matched[1].mc_pt)
                ? reco_photons_matched[0]
                : reco_photons_matched[1];
    },{"chosen_two"}).Define("leading_photon_pdg_id_origin",
    [](TruthParticle& leading_photon, RVec<TruthParticle>& truth_particles)
    {
        TruthParticle origin = leading_photon;
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
    }, {"leading_photon", "truth_particles"}).Define("subleading_photon_pdg_id_origin",
    [](TruthParticle& subleading_photon, RVec<TruthParticle>& truth_particles)
    {
        TruthParticle origin = subleading_photon;
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
    }, {"subleading_photon", "truth_particles"});
    
    auto leading_ids = diphotons.Take<int, RVec<int>>("leading_photon_pdg_id_origin");
    auto subleading_ids = diphotons.Take<int, RVec<int>>("subleading_photon_pdg_id_origin");
    
    std::unordered_map<int,int> leading_id_freqs, subleading_id_freqs;
    double total = *diphotons.Count();
    
    for (auto& i: *leading_ids)
    {
        leading_id_freqs[i]++;
    }
    
    for (auto& i: *subleading_ids)
    {
        subleading_id_freqs[i]++;
    }
    
    std::cout << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
    std::cout << R"--(\begin{tabular}{|c|c|})--" << '\n';
    std::cout << R"--(\hline)--" << '\n';
    std::cout << R"--(\multicolumn{2}{|c|}{$ee\gamma\gamma$ Leading Photon}\\ \hline)--" << '\n';
    std::cout << R"--(Pdg id & \% \\ \hline )--" << '\n';
    for (auto& i: leading_id_freqs)
    {
        std::cout << i.first << " & " << std::setprecision(2) << std::fixed
        << 100*(i.second/total) << R"--( \\ \hline)--" << '\n';
    }
    std::cout << R"--(\end{tabular}})--" << '\n';
    
    std::cout << "\n\n\n";
    
    std::cout << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
    std::cout << R"--(\begin{tabular}{|c|c|})--" << '\n';
    std::cout << R"--(\hline)--" << '\n';
    std::cout << R"--(\multicolumn{2}{|c|}{$ee\gamma\gamma$ Sub-Leading Photon}\\ \hline)--" << '\n';
    std::cout << R"--(Pdg id & \% \\ \hline )--" << '\n';
    for (auto& i: subleading_id_freqs)
    {
        std::cout << i.first << " & " << std::setprecision(2) << std::fixed
        << 100*(i.second/total) << R"--( \\ \hline)--" << '\n';
    }
    
    std::cout << R"--(\end{tabular}})--" << '\n';
    
    std::cout << "\n\n\n";
}

void Table10()
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
    [&](RVec<TruthParticle> truth_particles)
    {
        truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
        [](TruthParticle& x)
        {
            return ((abs(x.mc_pdg_id) != 22) || (abs(x.mc_eta) >= 2.37) || (x.mc_pt <= 10e3));

        }), truth_particles.end());
        
        return truth_particles;
        
    }, {"truth_particles"});
    
    auto diphotons = photon_passes_cuts.Define("chosen_two",
    [](RVec<TruthParticle>& reco_photons_matched)
    {
        RVec<TruthParticle> x;
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
            if ((delta_r < 1.5) && (X > 0.96) && (X < 1.2))
            {
                x = {reco_photons_matched[combs[0][i]], reco_photons_matched[combs[1][i]]};
                return x;
            }
        }
        return x;
    }, {"photons_pass_cuts"}).Filter(
    [&](RVec<TruthParticle>& reco_photons_matched)
    {
        if (reco_photons_matched.size()!=2)
        {
            return false;
        }
        return true;
    }, {"chosen_two"}).Define("leading_photon",
    [&](RVec<TruthParticle>& reco_photons_matched)
    {
        return (reco_photons_matched[0].mc_pt > reco_photons_matched[1].mc_pt)
                ? reco_photons_matched[0]
                : reco_photons_matched[1];
    },{"chosen_two"}).Define("subleading_photon",
    [&](RVec<TruthParticle>& reco_photons_matched)
    {
        return (reco_photons_matched[0].mc_pt < reco_photons_matched[1].mc_pt)
                ? reco_photons_matched[0]
                : reco_photons_matched[1];
    },{"chosen_two"}).Define("photon_pdg_id_origin",
    [](TruthParticle& leading_photon, TruthParticle& subleading_photon, RVec<TruthParticle>& truth_particles)
    {
        TruthParticle leading_origin = leading_photon;
        TruthParticle subleading_origin = subleading_photon;
        int leading_origin_id = 0, subleading_origin_id = 0;
        bool leading_found, subleading_found;
        
        do
        {
            leading_found = false;
            subleading_found = false;
            for (auto& i: truth_particles)
            {
                if (leading_origin.mc_parent_barcode == i.mc_barcode)
                {
                    leading_origin = i;
                    leading_origin_id = i.mc_pdg_id;
                    leading_found = true;
                }
                
                if (subleading_origin.mc_parent_barcode == i.mc_barcode)
                {
                    subleading_origin = i;
                    subleading_origin_id = i.mc_pdg_id;
                    subleading_found = true;
                }
            }
        } while (leading_found || subleading_found);
        
        return std::to_string(leading_origin_id)+"/"+std::to_string(subleading_origin_id);
    }, {"leading_photon", "subleading_photon", "truth_particles"});
    
    auto ids = diphotons.Take<std::string, RVec<std::string>>("photon_pdg_id_origin");
    
    std::unordered_map<std::string,int> id_freqs;
    double total = *diphotons.Count();
    
    for (auto& i: *ids)
    {
        id_freqs[i]++;
    }
    
    std::cout << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
    std::cout << R"--(\begin{tabular}{|c|c|})--" << '\n';
    std::cout << R"--(\hline)--" << '\n';
    std::cout << R"--(\multicolumn{2}{|c|}{$ee\gamma\gamma$ Leading/Subleading Photon Pairs}\\ \hline)--" << '\n';
    std::cout << R"--(Pdg id & \% \\ \hline )--" << '\n';
    for (auto& i: id_freqs)
    {
        std::cout << i.first << " & " << std::setprecision(2) << std::fixed
        << 100*(i.second/total) << R"--( \\ \hline)--" << '\n';
    }
    std::cout << R"--(\end{tabular}})--" << '\n';
    
    std::cout << "\n\n\n";
}

void Table16()
{
    std::vector<std::vector<std::string>> input_filenames =
    {
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617070._000001.LGNTuple.root"}, {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617064._000001.LGNTuple.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617074._000001.LGNTuple.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_data_test.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root"}
    };
    
    std::vector<std::string> prefixes = { R"--(pty2\_9\_17)--", R"--(pty\_17\_myy\_0\_80)--", R"--(pty\_17\_myy\_80)--", "data", R"--($\text{Sig } m_{A}$ = 5 GeV)--", R"--($\text{Sig } m_{A}$ = 1 GeV)--", R"--(Total Bkg)--"};
    
    std::ofstream out("Table16.txt");
    
    std::array<double,3> SFs = {((139e15)*(.871e-12))/150000.,((139e15)*(.199e-12))/150000., ((139e15)*(.0345e-15))/110465.};
       
    std::vector<ROOT::RDF::RResultHandle> Totals;
    
    for (auto& i: input_filenames)
    {
        SchottDataFrame df(MakeRDF(i,8));
//        std::cout << *df.Count() << '\n';
        
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
            
            reco_photons_matched.erase(std::remove_if(reco_photons_matched.begin(),reco_photons_matched.end(),
            [](Photon& x)
            {
                return x.photon_pt <= 10e3;

            }), reco_photons_matched.end());
            
            if (reco_photons_matched.size() == 1)
            {
                return true ? reco_photons_matched[0].photon_pt > 20e3 : false;
            }
            
            else if (reco_photons_matched.empty())
            {
                return false;
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
                if ((delta_r < 1.5) && (X > 0.96) && (X < 1.2))
                {
                    return false;
                }
            }
            
            for (auto& p: reco_photons_matched)
            {
                if (p.photon_pt > 20e3)
                {
                    return true;
                }
            }
            return false;
            
        }, {"photons_pass_cuts"})
        .Define("merged_photon",
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
            
        }, {"photons_pass_cuts"});
        
        auto dilepton_and_photon = merged_reco_photons_matched
        .Define("reconstructed_mass",[&](RVec<Electron>& di_electrons, Photon& merged_photon)
        {
            auto four_momentum = di_electrons[0].Vector() + di_electrons[1].Vector();
            
            return (four_momentum + merged_photon.Vector()).M()/1e3;
            
        }, {"di_electrons", "merged_photon"});
        
        auto pSB = dilepton_and_photon.Filter(
        [](double reconstructed_mass)
        {
            return (reconstructed_mass < 110) || (reconstructed_mass > 130);
        }, {"reconstructed_mass"});
        
        auto pSR = dilepton_and_photon.Filter(
        [](double reconstructed_mass)
        {
            return (reconstructed_mass >= 110) && (reconstructed_mass <= 130);
        }, {"reconstructed_mass"});
        
        
        Totals.push_back(merged_reco_photons_matched.Count());
        Totals.push_back(pSB.Count());
        Totals.push_back(pSR.Count());
        
//        if (count <= 2)
//        {
//            out << *merged_reco_photons_matched.Count()*SFs[count++] << '\n';
//        }
//        else
//        {
//            out << *merged_reco_photons_matched.Count() << '\n';
//            count++;
//        }
    }
    
    ROOT::RDF::RunGraphs(Totals);
    int count = 1, globCount = 0, otherCount = 0;
    std::vector<std::vector<double>> Vals(3);
    
    for (auto& i: Totals)
    {
        if (globCount < 2)
        {
//            std::cout << (*i.GetResultPtr<ULong64_t>())*(SFs[globCount]) << '\n';
            Vals[otherCount++].push_back((*i.GetResultPtr<ULong64_t>())*(SFs[globCount]));
        }
        else
        {
//            std::cout << *i.GetResultPtr<ULong64_t>() << '\n';
            Vals[otherCount++].push_back(*i.GetResultPtr<ULong64_t>());
        }
        if (count++%3==0)
        {
            otherCount = 0;
            globCount++;
        }
    }
    
    Vals[0].push_back(Vals[0][0]+Vals[0][1]+Vals[0][2]);
    Vals[1].push_back(Vals[1][0]+Vals[1][1]+Vals[1][2]);
    Vals[2].push_back(Vals[2][0]+Vals[2][1]+Vals[2][2]);
    
    std::cout << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
    std::cout << R"--(\begin{tabular}{|c|c|c|c|c|c|c|c|})--" << '\n';
    std::cout << R"--(\hline)--" << '\n';
//    std::vector<std::string> prefixes = { "pty2_9_17", "pty_17_myy_0_80", "pty_17_myy_80", "data", "Sig m_{A} = 5 GeV", "Sig m_{A} = 1 GeV"};
    std::cout << "{} & ";
    for (auto& i: prefixes)
    {
        if (i==prefixes.back())
        {
            std::cout << i << R"--(\\ \hline)--" << '\n';
        }
        else
        {
            std::cout << i << " & ";
        }
    }
    std::cout << '\n';
    count = 0;
    std::vector<std::string> rows = {"Full Reg", "pSB", "pSR"};
    for (auto& i: Vals)
    {
        std::cout << rows[count++] << " & ";
        for (auto& j: i)
        {
            if (j==i.back())
            {
                std::cout << j << R"--( \\ \hline)--" << '\n';
            }
            else
            {
                std::cout << j << " & ";
            }
            
//            std::cout << j << '\n';
        }
//        std::cout << R"--(\\ \hline)--" << '\n';
    }
    
    std::cout << R"--(\end{tabular}})--" << '\n';
    
    std::cout << "\n\n\n";
}

void Table19()
{
    std::vector<std::vector<std::string>> input_filenames =
    {
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617070._000001.LGNTuple.root"}, {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617064._000001.LGNTuple.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617074._000001.LGNTuple.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_data_test.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root"}
    };
    
    std::vector<std::string> prefixes = { R"--(pty2\_9\_17)--", R"--(pty\_17\_myy\_0\_80)--", R"--(pty\_17\_myy\_80)--", "data", R"--($\text{Sig } m_{A}$ = 5 GeV)--", R"--($\text{Sig } m_{A}$ = 1 GeV)--", R"--(Total Bkg)--"};
    
    std::ofstream out("Table16.txt");
    
    std::array<double,3> SFs = {((139e15)*(.871e-12))/150000.,((139e15)*(.199e-12))/150000., ((139e15)*(.0345e-15))/110465.};
       
    std::vector<ROOT::RDF::RResultHandle> Totals;
    
    for (auto& i: input_filenames)
    {
        SchottDataFrame df(MakeRDF(i,8));
//        std::cout << *df.Count() << '\n';
        
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
                return ((abs(x.photon_eta) >= 2.37) || (x.photon_pt <= 10e3) || (abs(x.photon_eta) > 1.37 && abs(x.photon_eta) < 1.52));

            }), photons.end());

            return photons;
        }, {"photons"});
        
        auto failed_resolved = photon_passes_cuts.Filter(
        [&](RVec<Photon>& reco_photons_test)
        {
            RVec<Photon> reco_photons_matched = reco_photons_test;
            
            reco_photons_matched.erase(std::remove_if(reco_photons_matched.begin(),reco_photons_matched.end(),
            [](Photon& x)
            {
                return x.photon_pt <= 10e3;

            }), reco_photons_matched.end());
            
            if (reco_photons_matched.size() == 1 || reco_photons_matched.empty())
            {
                return true;
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
                if ((delta_r < 1.5) && (X > 0.96) && (X < 1.2))
                {
                    return false;
                }
            }
            
            return true;
            
        }, {"photons_pass_cuts"});
        
        auto photon_pt_cut = failed_resolved.Filter(
        [&](RVec<Photon>& photon_passes_cuts)
        {
            for (auto& p: photon_passes_cuts)
            {
                if (p.photon_pt > 20e3)
                {
                    return true;
                }
            }
            
            return false;
            
        }, {"photons_pass_cuts"}).Define("merged_photon",
        [&](RVec<Photon>& photon_passes_cuts)
        {
            return photon_passes_cuts[0];
        }, {"photons_pass_cuts"});
        
        auto dilepton_and_photon = photon_pt_cut
        .Define("reconstructed_mass",[&](RVec<Electron>& di_electrons, Photon& merged_photon)
        {
            auto four_momentum = di_electrons[0].Vector() + di_electrons[1].Vector();
            
            return (four_momentum + merged_photon.Vector()).M()/1e3;
            
        }, {"di_electrons", "merged_photon"});
        
        auto pSR = dilepton_and_photon.Filter(
        [](double reconstructed_mass)
        {
            return (reconstructed_mass >= 110) && (reconstructed_mass <= 130);
        }, {"reconstructed_mass"});
        
        auto pSR_ID = pSR.Filter(
        [&](Photon& merged_photon)
        {
            return merged_photon.photon_id;
        }, {"merged_photon"});
        
        
        
        Totals.push_back(ptCut.Count()); //preselection
        Totals.push_back(failed_resolved.Count()); //failed_resolved
        Totals.push_back(photon_pt_cut.Count()); //photon_pt_cut
        Totals.push_back(pSR.Count()); //pSR
        Totals.push_back(pSR_ID.Count()); //pSR_ID
        
//        if (count <= 2)
//        {
//            out << *merged_reco_photons_matched.Count()*SFs[count++] << '\n';
//        }
//        else
//        {
//            out << *merged_reco_photons_matched.Count() << '\n';
//            count++;
//        }
    }
    
    ROOT::RDF::RunGraphs(Totals);
    int count = 1, globCount = 0, otherCount = 0;
    std::vector<std::vector<double>> Vals(5);
    
    for (auto& i: Totals)
    {
        if (globCount < 2)
        {
//            std::cout << (*i.GetResultPtr<ULong64_t>())*(SFs[globCount]) << '\n';
            Vals[otherCount++].push_back((*i.GetResultPtr<ULong64_t>())*(SFs[globCount]));
        }
        else
        {
//            std::cout << *i.GetResultPtr<ULong64_t>() << '\n';
            Vals[otherCount++].push_back(*i.GetResultPtr<ULong64_t>());
        }
        if (count++%5==0)
        {
            otherCount = 0;
            globCount++;
        }
    }
    
    Vals[0].push_back(Vals[0][0]+Vals[0][1]+Vals[0][2]);
    Vals[1].push_back(Vals[1][0]+Vals[1][1]+Vals[1][2]);
    Vals[2].push_back(Vals[2][0]+Vals[2][1]+Vals[2][2]);
    Vals[3].push_back(Vals[3][0]+Vals[3][1]+Vals[3][2]);
    Vals[4].push_back(Vals[4][0]+Vals[4][1]+Vals[4][2]);
    
    std::cout << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
    std::cout << R"--(\begin{tabular}{|c|c|c|c|c|c|c|c|})--" << '\n';
    std::cout << R"--(\hline)--" << '\n';
//    std::vector<std::string> prefixes = { "pty2_9_17", "pty_17_myy_0_80", "pty_17_myy_80", "data", "Sig m_{A} = 5 GeV", "Sig m_{A} = 1 GeV"};
    std::cout << "{} & ";
    for (auto& i: prefixes)
    {
        if (i==prefixes.back())
        {
            std::cout << i << R"--(\\ \hline)--" << '\n';
        }
        else
        {
            std::cout << i << " & ";
        }
    }
    std::cout << '\n';
    count = 0;
    std::vector<std::string> rows = {"pass preselection", "failed resolved category", R"--(photon $p_T$ cut)--", "pSR", "pSR-ID"};
    for (auto& i: Vals)
    {
        std::cout << rows[count++] << " & ";
        for (auto& j: i)
        {
            if (j==i.back())
            {
                std::cout << j << R"--( \\ \hline)--" << '\n';
            }
            else
            {
                std::cout << j << " & ";
            }
            
//            std::cout << j << '\n';
        }
//        std::cout << R"--(\\ \hline)--" << '\n';
    }
    
    std::cout << R"--(\end{tabular}})--" << '\n';
    
    std::cout << "\n\n\n";
}

void DataBackgroundComparison()
{
    auto start_time = Clock::now();
//    fig27();
//    fig28();
//    fig41();
//    fig48();
//    fig59();
//    Table9();
//    Table10();
//    Table16();
    Table19();
    auto end_time = Clock::now();
    std::cout << "Time difference: "
       << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count()/1e9 << " nanoseconds" << std::endl;
    
}

int main()
{
    DataBackgroundComparison();
}

ROOT.gInterpreter.Declare("""
using namespace ROOT::VecOps;
RVec<double> someFunc(RVec<double> l_tlv_pt, RVec<int>& l_pdg_id)
{
    RVec<double> new_pt;
    
    for (size_t i = 0; i < l_tlv_pt.size(); i++)
    {
        if (abs(l_pdg_id[i])==11)
        {
            new_pt.push_back(l_tlv_pt[i]);
        }
    }
    return new_pt;
}
""")
