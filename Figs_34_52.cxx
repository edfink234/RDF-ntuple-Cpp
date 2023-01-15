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
    
    auto resolved = photon_passes_cuts.Define("chosen_two",
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
    }, {"chosen_two"}).Filter(
    [&](RVec<Photon>& chosen_two, RVec<Electron>& electrons)
    {
        PtEtaPhiEVector gg = chosen_two[0].Vector()+chosen_two[1].Vector();
        PtEtaPhiEVector ll = electrons[0].Vector()+electrons[1].Vector();
        double m_llgg = (ll.M() + gg.M())/1e3;
        
        return !((110 <= m_llgg) && (m_llgg <= 140));
        
    }, {"chosen_two", "di_electrons"}).Define("m_gg",
    [&](RVec<Photon>& chosen_two)
    {
        PtEtaPhiEVector gg = chosen_two[0].Vector()+chosen_two[1].Vector();
        return gg.M()/1e3;
    }, {"chosen_two"});
    
    auto di_ph_mm = resolved.Histo1D<double>({"data", "data", 10u, 0.6, 7.2}, "m_gg");
    TCanvas* c1 = new TCanvas();
    TLegend* legend = new TLegend(0.65, 0.4, 0.75, 0.6);
    di_ph_mm->SetLineColor(kBlue);
    legend->AddEntry(&(*di_ph_mm), di_ph_mm->GetTitle(), "l");
    di_ph_mm->SetTitle(";m_{#gamma#gamma} [GeV];Events");
    di_ph_mm->GetYaxis()->CenterTitle(true);
    di_ph_mm->GetXaxis()->SetTitleOffset(1.2);
    di_ph_mm->Draw("same");
    gStyle->SetOptStat(0);
    TLatex Tl;
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.6, 0.8, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.6, 0.7,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("Fig34.png");
    
}

void Fig52()
{
    std::vector<std::string> input_filenames =
    {
        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617070._000001.LGNTuple.root", "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617064._000001.LGNTuple.root",
        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617074._000001.LGNTuple.root",
        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_data_test.root",
    };
    
    std::vector<const char*> prefixes = {"pty2_9_17", "pty_17_myy_0_80", "pty_17_myy_80", "data"};
    std::vector<EColor> colors = {kBlue, kRed, kViolet, kGreen};
    std::array<double,3> SFs = {((139e15)*(.871e-12))/150000.,((139e15)*(.199e-12))/150000., ((139e15)*(.0345e-15))/110465.};
    std::vector<ROOT::RDF::RResultPtr<TH1D>> histos_back_plus_data;
    histos_back_plus_data.reserve(4);
    std::vector<ROOT::RDF::RResultPtr<ULong64_t>> backCounts;
    backCounts.reserve(3);
    
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
        [](RVec<Photon>& reco_photons_matched)
        {
            for (auto& p: reco_photons_matched)
            {
                if (p.photon_pt > 20e3)
                {
                    return p;
                }
            }
            
            return reco_photons_matched[0]; //jic the compiler complains
            
        }, {"photons_pass_cuts"}).Define("merged_photon_pt",
        [](Photon& p)
        {
            return p.photon_pt/1e3;
        }, {"merged_photon"}).Define("reconstructed_mass",
        [&](RVec<Electron>& di_electrons, Photon& merged_photon)
        {
            auto four_momentum = di_electrons[0].Vector() + di_electrons[1].Vector();

            return (four_momentum + merged_photon.Vector()).M()/1e3;

        }, {"di_electrons", "merged_photon"}).Filter(
        [](double reconstructed_mass, RVec<float>& Eratio)
        {
            return ((!Any(Eratio <= 0.8)) && ((reconstructed_mass < 110) || (reconstructed_mass > 130)));
        }, {"reconstructed_mass", "photon_shower_shape_e_ratio"});

        if (count <= 2)
        {
            backCounts.push_back(merged_reco_photons_matched.Count());
        }
        
        histos_back_plus_data.push_back(merged_reco_photons_matched.Histo1D<double>({prefixes[count], prefixes[count++], 100u, 0, 165}, "merged_photon_pt"));
    }
    
    double factor;
    for (auto& i: backCounts)
    {
        factor += *i;
    }
    
    auto hs = new THStack("hs3","");
    TCanvas* c1 = new TCanvas();
    TLegend* legend = new TLegend(0.65, 0.4, 0.85, 0.6);
    count=0;
    for (auto& h: histos_back_plus_data)
    {
        if (h->Integral() != 0 && &h != &histos_back_plus_data.back())
        {
            h->Scale((factor/h->Integral())*SFs[count]);
        }
        h->SetFillColor(colors[count++]);
        legend->AddEntry(&(*h), h->GetTitle(), "f");
    
        if (&h != &histos_back_plus_data.back())
        {
            hs->Add(&*h);
        }
    }

    hs->Draw("HIST");
    histos_back_plus_data[3]->Draw("HISTsame");
    hs->SetTitle(";photon p_{T} [GeV];Events");
    hs->GetYaxis()->CenterTitle(true);
    hs->GetXaxis()->SetTitleOffset(1.2);
    gStyle->SetOptStat(0);
    TLatex Tl;
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.6, 0.8, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.6, 0.7,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
    legend->SetBorderSize(0);
    legend->Draw("same");
    c1->SaveAs("Fig52A.png");
}

void Figs_34_52()
{
    auto start_time = Clock::now();
    Fig34();
    Fig52();
    auto end_time = Clock::now();
    std::cout << "Time difference: "
       << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count()/1e9 << " nanoseconds" << std::endl;
}

int main()
{
    Figs_34_52();
}
