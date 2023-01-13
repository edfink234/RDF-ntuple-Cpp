#include <iostream>
#include <vector>
#include <chrono>
#include <algorithm>
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
#include "TStyle.h"
#include "TPaveText.h"
#include "TLatex.h"
#include "TLegend.h"
#include "Rtypes.h"
#include "THStack.h"

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
//Other mA5 file: /Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600909.PhPy8EG_AZNLO_ggH125_mA5p0_Cyy0p01_Czh1p0.merge.AOD.e8324_e7400_s3126_r10724_r10726_v2.root
void fig1A()
{
    std::vector<std::string> input_filenames =
    {
        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617070._000001.LGNTuple.root", "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617064._000001.LGNTuple.root",
        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617074._000001.LGNTuple.root",
        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root", "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root",
    };
    
    std::vector<ROOT::RDF::RResultPtr<ULong64_t>> backCounts;
    std::vector<ROOT::RDF::RResultPtr<TH1D>> histos, back_histos;
    histos.reserve(input_filenames.size());
    std::vector<const char*> prefixes = {"pty2_9_17", "pty_17_myy_0_80", "pty_17_myy_80", "sig m_{A} = 5 GeV", "sig m_{A} = 1 GeV"};
    std::array<double,3> SFs = {((139e15)*(.871e-12))/150000.,((139e15)*(.199e-12))/150000., ((139e15)*(.0345e-15))/110465.};
    std::vector<EColor> colors = {kBlue, kRed, kViolet, kBlack, kMagenta};
    TCanvas* c1 = new TCanvas();
    TLegend* legend = new TLegend(0.6, 0.4, 0.8, 0.6);
    int count = 0;
    
    for (auto& file: input_filenames)
    {
        SchottDataFrame df(MakeRDF({file}, 8));
//        df.Describe().Print();
//        exit(1);
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
            
        auto two_leptons = trigger_selection.Filter(
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
            
        auto dilep_mass = ptCut.Define("dilep_mass",[&](RVec<Electron>& electrons)
        {
            
            PtEtaPhiEVector four_momentum =
            electrons[0].Vector() + electrons[1].Vector();
            
            return four_momentum.M()/1e3;
            
        }, {"electrons"});
           
        if (count <= 2)
        {
//            break;
            back_histos.push_back(dilep_mass.Histo1D<double>({prefixes[count], prefixes[count++], 100u, 60, 120}, "dilep_mass"));
            backCounts.push_back(dilep_mass.Count());
//            auto passed = dilep_mass.Count();
//            std::cout << *passed << '\n';
        }
        else
        {
            histos.push_back(dilep_mass.Histo1D<double>({prefixes[count], prefixes[count++], 100u, 60, 120}, "dilep_mass"));
//            auto passed = dilep_mass.Count();
//            std::cout << *passed << '\n';
        }
    }
    count = 0;
    double factor;
    double total_back = 0;
    for (auto& i: backCounts)
    {
        total_back += *i;
    }
    int back_count = 0;
    factor = total_back;
    auto hs = new THStack("hs3","");
    for (auto& h: back_histos)
    {
        h->SetFillColor(colors[count]);
        legend->AddEntry(&(*h), h->GetTitle(), "f");

        h->Scale((factor/h->Integral())*SFs[count]);
        h->SetTitle(";m_{ll}  [GeV];Events");
        h->GetYaxis()->CenterTitle(true);
        h->SetAxisRange(0., 1000,"Y");
//        gPad->Modified(); gPad->Update();
        
        hs->Add(&*h);
        h->SetAxisRange(0., 1000,"Y");
        count++;
    }
    hs->Draw("HIST");
    
    for (auto& h: histos)
    {
        h->SetLineColor(colors[count++]);
        legend->AddEntry(&(*h), h->GetTitle(), "l");

        h->Scale(factor/h->Integral());
        h->SetTitle(";m_{ll}  [GeV];Events");
        h->GetYaxis()->CenterTitle(true);
//            h->SetAxisRange(0., 3300,"Y");
        h->Scale(factor/h->Integral());
        h->Draw("HISTsame");
        gPad->Modified(); gPad->Update();
    }
    
    hs->SetTitle(";m_{ll}  [GeV];Events");
    hs->GetYaxis()->CenterTitle(true);
    hs->GetXaxis()->SetTitleOffset(1.2);
    hs->SetMinimum(0);
    hs->SetMaximum(600.);
    
    gStyle->SetOptStat(0);
    TLatex Tl;
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.6, 0.83, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.6, 0.73,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("Fig1A.png");
}

void fig5()
{
    std::vector<std::string> input_filenames = { "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root",
    };

    std::vector<EColor> colors = {kRed, kBlue, kBlack};
    TCanvas* c1 = new TCanvas();
    TLegend* legend = new TLegend(0.6, 0.4, 0.8, 0.525);
    
    SchottDataFrame df(MakeRDF(input_filenames, 8));
    
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
    
    auto leading_truth_photons = preselection.Define("leading_truth_photons",[&](RVec<TruthParticle> truth_particles)
    {
        truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
        [](TruthParticle& x)
        {
            return (abs(x.mc_pdg_id) != 22);
            
        }), truth_particles.end());
        
        std::sort(truth_particles.begin(),truth_particles.end(),
        [](TruthParticle& x, TruthParticle& y)
        {
            return x.mc_pt > y.mc_pt;
        });
        
        RVec<TruthParticle> chosen = {truth_particles[0], truth_particles[1]};
        
        return chosen;
        
    }, {"truth_particles"})
    .Define("leading_truth_photon_pt", [&](RVec<TruthParticle> leading_truth_photons)
    {
        return leading_truth_photons[0].mc_pt/1e3;
    }, {"leading_truth_photons"})
    .Define("subleading_truth_photon_pt", [&](RVec<TruthParticle> leading_truth_photons)
    {
        return leading_truth_photons[1].mc_pt/1e3;
    }, {"leading_truth_photons"});
    
//    auto leading_truth_photons_passed = leading_truth_photons.Count();
//    std::cout << *leading_truth_photons_passed << '\n';
    
    auto leading_truth_photons_pt_10 = leading_truth_photons.Filter(
    [&](RVec<TruthParticle>& leading_truth_photons)
    {
        return (leading_truth_photons[0].mc_pt > 10e3 && leading_truth_photons[1].mc_pt > 10e3);
    }, {"leading_truth_photons"})
    .Define("leading_truth_photons_pt_10_dR",
    [&](RVec<TruthParticle>& leading_truth_photons)
    {
        return DeltaR(leading_truth_photons[0].Vector(), leading_truth_photons[1].Vector());
    }, {"leading_truth_photons"});
    
//    auto leading_truth_photons_pt_10_passed = leading_truth_photons_pt_10.Count();
//    std::cout << *leading_truth_photons_pt_10_passed << '\n';
    
    auto leading_truth_photons_pt_5 = leading_truth_photons.Filter(
    [&](RVec<TruthParticle>& leading_truth_photons)
    {
        return (leading_truth_photons[0].mc_pt > 5e3 && leading_truth_photons[1].mc_pt > 5e3);
    }, {"leading_truth_photons"})
    .Define("leading_truth_photons_pt_5_dR",
    [&](RVec<TruthParticle>& leading_truth_photons)
    {
        return DeltaR(leading_truth_photons[0].Vector(), leading_truth_photons[1].Vector());
    }, {"leading_truth_photons"});
    
//    auto leading_truth_photons_pt_5_passed = leading_truth_photons_pt_5.Count();
//    std::cout << *leading_truth_photons_pt_5_passed << '\n';
    
    auto leading_truth_photons_pt_0 = leading_truth_photons.Filter(
    [&](RVec<TruthParticle>& leading_truth_photons)
    {
        return (leading_truth_photons[0].mc_pt > 0 && leading_truth_photons[1].mc_pt > 0);
    }, {"leading_truth_photons"})
    .Define("leading_truth_photons_pt_0_dR",
    [&](RVec<TruthParticle>& leading_truth_photons)
    {
        return DeltaR(leading_truth_photons[0].Vector(), leading_truth_photons[1].Vector());
    }, {"leading_truth_photons"});
    
//    auto leading_truth_photons_pt_0_passed = leading_truth_photons_pt_0.Count();
//    std::cout << *leading_truth_photons_pt_0_passed << '\n';
    
    std::vector<ROOT::RDF::RResultPtr<TH1D>> histos =
    {
        leading_truth_photons.Histo1D<double>({"Leading Photon", "Leading Photon", 20u, 0, 50}, "leading_truth_photon_pt"),
        leading_truth_photons.Histo1D<double>({"Sub-Leading Photon", "Sub-Leading Photon", 20u, 0, 50}, "subleading_truth_photon_pt"),
        leading_truth_photons_pt_0.Histo1D<double>({"Photon p_{T} > 0 GeV", "Photon p_{T} > 0 GeV", 40u, 0, 0.32}, "leading_truth_photons_pt_0_dR"),
        leading_truth_photons_pt_5.Histo1D<double>({"Photon p_{T} > 5 GeV", "Photon p_{T} > 5 GeV", 40u, 0, 0.32}, "leading_truth_photons_pt_5_dR"),
        leading_truth_photons_pt_10.Histo1D<double>({"Photon p_{T} > 10 GeV", "Photon p_{T} > 10 GeV", 40u, 0, 0.32}, "leading_truth_photons_pt_10_dR"),
    };
    
    double factor = histos[0]->Integral();
//    std::cout << factor << "  " << histos[0]-> GetEntries();
    histos[0]->Scale(factor/factor);
    histos[0]->SetLineColor(colors[0]);
    histos[0]->Draw("HIST");
    legend->AddEntry(&(*histos[0]), histos[0]->GetTitle(), "l");
    histos[0]->SetTitle(";Truth photon p_{T}  [GeV];Events");
    histos[0]->SetTitleOffset(1.2);
    histos[0]->GetYaxis()->CenterTitle(true);
    histos[0]->SetAxisRange(0., 1800,"Y");
    
    histos[1]->Scale(factor/histos[1]->Integral());
    histos[1]->SetLineColor(colors[1]);
    histos[1]->Draw("HISTsame");
    legend->AddEntry(&(*histos[1]), histos[1]->GetTitle(), "l");
    
    gStyle->SetOptStat(0);
    TLatex Tl;
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.6, 0.83, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.6, 0.77,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
    Tl.DrawLatexNDC(0.6, 0.7,"ggF m_{A} = 1 GeV");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("Fig5A.png");

    c1 = new TCanvas();
    legend = new TLegend(0.6, 0.4, 0.8, 0.6);
    factor = histos[2]->Integral();
    
    histos[2]->Scale(factor/factor);
    histos[2]->SetLineColor(colors[0]);
    histos[2]->Draw("HIST");
    
    histos[2]->GetYaxis()->CenterTitle(true);
    histos[2]->SetAxisRange(0., 1200,"Y");
    
    legend->AddEntry(&(*histos[2]), histos[2]->GetTitle(), "l");
    histos[2]->SetTitle(";Truth #DeltaR_{#gamma#gamma}; Events");

    histos[3]->Scale(factor/histos[3]->Integral());
    histos[3]->SetLineColor(colors[1]);
    histos[3]->Draw("HISTsame");
    legend->AddEntry(&(*histos[3]), histos[3]->GetTitle(), "l");
    
    histos[4]->Scale(factor/histos[4]->Integral());
    histos[4]->SetLineColor(colors[2]);
    histos[4]->Draw("HISTsame");
    legend->AddEntry(&(*histos[4]), histos[4]->GetTitle(), "l");
    
    gStyle->SetOptStat(0);
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.6, 0.83, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.6, 0.77,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
    Tl.DrawLatexNDC(0.6, 0.7,"ggF m_{A} = 1 GeV");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("Fig5B.png");
}

void fig6()
{
    std::vector<std::string> input_filenames = { "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root",
    };

    std::vector<EColor> colors = {kRed, kBlue};
    TCanvas* c1 = new TCanvas();
    TLegend* legend = new TLegend(0.6, 0.4, 0.8, 0.525);
    
    SchottDataFrame df(MakeRDF(input_filenames, 8));
    
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
            return ((abs(x.photon_eta) >= 2.37) || (x.photon_pt <= 10e3) || (abs(x.photon_eta) > 1.37 && abs(x.photon_eta) < 1.52) || (!x.photon_id_loose));

        }), photons.end());
        
        return photons;
    }, {"photons"});
    
    auto reco_photons_matched = truth_photons_from_axions.Define("reco_photons_matched",
    [&](RVec<Photon>& photons_pass_cuts, RVec<TruthParticle>& truth_photons_from_axions)
    {
        RVec<Photon> matchedPhotons;
        PtEtaPhiEVector tp1 = truth_photons_from_axions[0].Vector();
        PtEtaPhiEVector tp2 = truth_photons_from_axions[1].Vector();
        
        for (auto& rp: photons_pass_cuts)
        {
            if (!(DeltaR(rp.Vector(), tp1) < 0.1 || DeltaR(rp.Vector(), tp2) < 0.1))
            {
                continue;
            }
            matchedPhotons.push_back(rp);
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
        
    }, {"reco_photons_matched"})
    .Define("leading_pt",[&] (RVec<Photon>& reco_photons_matched)
    {
        return reco_photons_matched[0].photon_pt/1e3;
    }, {"reco_photons_matched"})
    .Define("sub_leading_pt",[&] (RVec<Photon>& reco_photons_matched)
    {
        return reco_photons_matched[1].photon_pt/1e3;
    }, {"reco_photons_matched"});
    
    auto one_reco_photons_matched = reco_photons_matched.Filter(
    [&](RVec<Photon>& reco_photons_matched)
    {
        return reco_photons_matched.size()==1;
        
    }, {"reco_photons_matched"})
    .Define("leading_pt",[&] (RVec<Photon>& reco_photons_matched)
    {
        return reco_photons_matched[0].photon_pt/1e3;
    }, {"reco_photons_matched"});
    
    std::vector<ROOT::RDF::RResultPtr<TH1D>> histos =
    {
        two_reco_photons_matched.Histo1D<double>({"Leading Photon", "Leading Photon", 20u, 0, 50}, "leading_pt"),
        two_reco_photons_matched.Histo1D<double>({"Sub-Leading Photon", "Sub-Leading Photon", 20u, 0, 50}, "sub_leading_pt"),
        one_reco_photons_matched.Histo1D<double>({"Leading Photon", "Leading Photon", 40u, 0, 100}, "leading_pt"),
    };
    
    double factor = histos[0]->Integral();
//    std::cout << factor << "  " << histos[0]-> GetEntries();
    histos[0]->Scale(factor/factor);
    histos[0]->SetLineColor(colors[0]);
    histos[0]->Draw("HIST");
    legend->AddEntry(&(*histos[0]), histos[0]->GetTitle(), "l");
    histos[0]->SetTitle(";photon p_{T}  [GeV];Events");
    histos[0]->SetTitleOffset(1.2);
    histos[0]->GetYaxis()->CenterTitle(true);
    histos[0]->SetAxisRange(0., 80,"Y");
    
    histos[1]->Scale(factor/histos[1]->Integral());
    histos[1]->SetLineColor(colors[1]);
    histos[1]->Draw("HISTsame");
    legend->AddEntry(&(*histos[1]), histos[1]->GetTitle(), "l");
    
    gStyle->SetOptStat(0);
    TLatex Tl;
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.6, 0.83, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.6, 0.77,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
    Tl.DrawLatexNDC(0.6, 0.7,"ggF m_{A} = 1 GeV");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("Fig6A.png");
    
    c1 = new TCanvas();
    legend = new TLegend(0.6, 0.4, 0.8, 0.6);
    factor = histos[2]->Integral();
    histos[2]->Scale(factor/factor);
    histos[2]->SetLineColor(colors[0]);
    histos[2]->Draw("HIST");
    histos[2]->GetYaxis()->CenterTitle(true);
//    histos[2]->SetAxisRange(0., 1500,"Y");
    
    legend->AddEntry(&(*histos[2]), histos[2]->GetTitle(), "l");
    histos[2]->SetTitle(";photon p_{T}  [GeV];Events");
    
    gStyle->SetOptStat(0);
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.6, 0.83, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.6, 0.77,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
    Tl.DrawLatexNDC(0.6, 0.7,"ggF m_{A} = 1 GeV");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("Fig6B.png");
    
//    std::cout << *truth_photons_from_axions.Count() << '\n';
   
}

void fig8()
{
    std::vector<std::string> input_filenames = {
        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root",
        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root"
    };
    
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
    
    std::vector<ROOT::RDF::RResultPtr<TH1D>> histos_Z;
    histos_Z.reserve(input_filenames.size());
    std::vector<ROOT::RDF::RResultPtr<TH1D>> histos_a;
    histos_a.reserve(input_filenames.size());
    std::vector<ROOT::RDF::RResultPtr<TH1D>> histos_h;
    histos_h.reserve(input_filenames.size());
    std::vector<ROOT::RDF::RResultPtr<TH1D>> histos_deltaR;
    histos_deltaR.reserve(input_filenames.size());
    std::vector<ROOT::RDF::RResultPtr<TH1D>> histos_deltaPhi;
    histos_deltaPhi.reserve(input_filenames.size());
    std::vector<ROOT::RDF::RResultPtr<TH1D>> histos_deltaEta;
    histos_deltaEta.reserve(input_filenames.size());

    std::vector<const char*> prefixes = {"sig ggF m_{A} = 5 GeV", "sig ggF m_{A} = 1 GeV"};
    std::vector<EColor> colors = {kBlack, kMagenta};
    int count = 0;
    
    for (auto& file: input_filenames)
    {
        SchottDataFrame df(MakeRDF({file}, 8));
        
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
                return ((abs(x.photon_eta) >= 2.37) || (x.photon_pt <= 10e3) || (abs(x.photon_eta) > 1.37 && abs(x.photon_eta) < 1.52) || (!x.photon_id_loose));

            }), photons.end());
            
            return photons;
        }, {"photons"});
        
        auto reco_photons_matched = truth_photons_from_axions.Define("reco_photons_matched",
        [&](RVec<Photon>& photons_pass_cuts, RVec<TruthParticle>& truth_photons_from_axions)
        {
            RVec<Photon> matchedPhotons;
            PtEtaPhiEVector tp1 = truth_photons_from_axions[0].Vector();
            PtEtaPhiEVector tp2 = truth_photons_from_axions[1].Vector();
            
            for (auto& rp: photons_pass_cuts)
            {
                if (!(DeltaR(rp.Vector(), tp1) < 0.1 || DeltaR(rp.Vector(), tp2) < 0.1))
                {
                    continue;
                }
                matchedPhotons.push_back(rp);
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
        
//        std::cout << *two_reco_photons_matched.Count() << '\n';
        
        auto stable_truth_dileptons_and_diphotons = two_reco_photons_matched
        .Define("stable_truth_leptons",[&](RVec<TruthParticle> truth_particles)
        {
            truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
            [](TruthParticle& x)
            {
                return (((abs(x.mc_pdg_id) != 11 && abs(x.mc_pdg_id) != 12 && abs(x.mc_pdg_id) != 13 &&
                        abs(x.mc_pdg_id) != 14 && abs(x.mc_pdg_id) != 15 && abs(x.mc_pdg_id) != 16 &&
                         abs(x.mc_pdg_id) != 17 && abs(x.mc_pdg_id) != 18)) || (x.mc_status != 1));
                
            }), truth_particles.end());

            return truth_particles;
            
        }, {"truth_particles"})
        .Filter([&](RVec<TruthParticle>& stable_truth_leptons, RVec<TruthParticle>& truth_photons_from_axions)
        {
            return (stable_truth_leptons.size()==2 && truth_photons_from_axions.size()==2);
            
        }, {"stable_truth_leptons", "truth_photons_from_axions"})
        .Define("reconstructed_Z",[&](RVec<TruthParticle>& stable_truth_leptons)
        {
            auto four_momentum = stable_truth_leptons[0].Vector() + stable_truth_leptons[1].Vector();
            
            return four_momentum.Pt()/1e3;
            
        }, {"stable_truth_leptons"})
        .Define("reconstructed_a",[&](RVec<TruthParticle>& truth_photons_from_axions)
        {
            PtEtaPhiEVector four_momentum = truth_photons_from_axions[0].Vector() + truth_photons_from_axions[1].Vector();
            return four_momentum.Pt()/1e3;

        }, {"truth_photons_from_axions"})
        .Define("reconstructed_h",[&](RVec<TruthParticle>& truth_photons_from_axions, RVec<TruthParticle>& stable_truth_leptons)
        {
            PtEtaPhiEVector four_momentum_photons = truth_photons_from_axions[0].Vector() + truth_photons_from_axions[1].Vector();
            PtEtaPhiEVector four_momentum_leptons = stable_truth_leptons[0].Vector() + stable_truth_leptons[1].Vector();
            return (four_momentum_photons + four_momentum_leptons).Pt()/1e3;
        }, {"truth_photons_from_axions", "stable_truth_leptons"})
        .Define("delta_r_Z_a",[&](RVec<TruthParticle>& truth_photons_from_axions, RVec<TruthParticle>& stable_truth_leptons)
        {
            PtEtaPhiEVector four_momentum_photons = truth_photons_from_axions[0].Vector() + truth_photons_from_axions[1].Vector();
            PtEtaPhiEVector four_momentum_leptons = stable_truth_leptons[0].Vector() + stable_truth_leptons[1].Vector();
            return DeltaR(four_momentum_photons, four_momentum_leptons);
        }, {"truth_photons_from_axions", "stable_truth_leptons"})
        .Define("delta_phi_Z_a",[&](RVec<TruthParticle>& truth_photons_from_axions, RVec<TruthParticle>& stable_truth_leptons)
        {
            PtEtaPhiEVector four_momentum_photons = truth_photons_from_axions[0].Vector() + truth_photons_from_axions[1].Vector();
            PtEtaPhiEVector four_momentum_leptons = stable_truth_leptons[0].Vector() + stable_truth_leptons[1].Vector();

            return ROOT::Math::VectorUtil::DeltaPhi(four_momentum_photons, four_momentum_leptons);
        }, {"truth_photons_from_axions", "stable_truth_leptons"})
        .Define("delta_eta_Z_a",[&](RVec<TruthParticle>& truth_photons_from_axions, RVec<TruthParticle>& stable_truth_leptons)
        {
            PtEtaPhiEVector four_momentum_photons = truth_photons_from_axions[0].Vector() + truth_photons_from_axions[1].Vector();
            PtEtaPhiEVector four_momentum_leptons = stable_truth_leptons[0].Vector() + stable_truth_leptons[1].Vector();
            
            return abs((four_momentum_photons - four_momentum_leptons).Eta());
        }, {"truth_photons_from_axions", "stable_truth_leptons"});
        
//        std::cout << *stable_truth_dileptons_and_diphotons.Count() << '\n';
        
        histos_Z.push_back(stable_truth_dileptons_and_diphotons.Histo1D<double>({prefixes[count], prefixes[count], 50u, 0, 200}, "reconstructed_Z"));
        histos_a.push_back(stable_truth_dileptons_and_diphotons.Histo1D<double>({prefixes[count], prefixes[count], 50u, 0, 62}, "reconstructed_a"));
        histos_h.push_back(stable_truth_dileptons_and_diphotons.Histo1D<double>({prefixes[count], prefixes[count], 100u, 0, 217}, "reconstructed_h"));
        histos_deltaR.push_back(stable_truth_dileptons_and_diphotons.Histo1D<double>({prefixes[count], prefixes[count], 100u, 0, 6.4}, "delta_r_Z_a"));
        histos_deltaPhi.push_back(stable_truth_dileptons_and_diphotons.Histo1D<double>({prefixes[count], prefixes[count], 100u, 0, 6.25}, "delta_phi_Z_a"));
        histos_deltaEta.push_back(stable_truth_dileptons_and_diphotons.Histo1D<double>({prefixes[count], prefixes[count++], 100u, 0, 6.25}, "delta_eta_Z_a"));
        
        auto passed = stable_truth_dileptons_and_diphotons.Count();
        std::cout << *passed << '\n';
    }
    
    TCanvas* c1 = new TCanvas();
    TLegend* legend = new TLegend(0.6, 0.4, 0.8, 0.6);
    double factor;
    
    count = 0;
    for (auto& h: histos_Z)
    {
        h->SetLineColor(colors[count++]);
        legend->AddEntry(&(*h), h->GetTitle(), "l");
        
        if (&h == &histos_Z.front())
        {
            factor = h->Integral();
            h->Scale(factor/h->Integral());
            h->SetTitle(";ll p_{T} [GeV];Events");
            h->GetYaxis()->CenterTitle(true);
            h->SetAxisRange(0., 48,"Y");
//            h->SetAxisRange(0., 300,"Y");
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
    Tl.DrawLatexNDC(0.6, 0.83, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.6, 0.73,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("Fig8A.png");
    
    c1 = new TCanvas();
    legend = new TLegend(0.65, 0.45, 0.85, 0.65);
    count = 0;
    for (auto& h: histos_a)
    {
        h->SetLineColor(colors[count++]);
        legend->AddEntry(&(*h), h->GetTitle(), "l");
        
        if (&h == &histos_Z.front())
        {
            factor = h->Integral();
            h->Scale(factor/h->Integral());
            h->SetTitle(";#gamma#gamma p_{T} [GeV];Events");
            h->GetYaxis()->CenterTitle(true);
            h->SetAxisRange(0., 38,"Y");
//            h->SetAxisRange(0., 250,"Y");
            h->Draw("HIST");
        }
        else
        {
            h->Scale(factor/h->Integral());
            h->SetTitle(";#gamma#gamma p_{T} [GeV];Events");
            h->GetYaxis()->CenterTitle(true);
            h->SetAxisRange(0., 38,"Y");
//            h->SetAxisRange(0., 250,"Y");
            h->Draw("HISTsame");
            gPad->Modified(); gPad->Update();
        }
    }
    gStyle->SetOptStat(0);
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.6, 0.83, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.6, 0.73,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("Fig8B.png");
    
    c1 = new TCanvas();
    legend = new TLegend(0.65, 0.45, 0.85, 0.65);
    count = 0;
    for (auto& h: histos_h)
    {
        h->SetLineColor(colors[count++]);
        legend->AddEntry(&(*h), h->GetTitle(), "l");
        
        if (&h == &histos_Z.front())
        {
            factor = h->Integral();
            h->Scale(factor/h->Integral());
            h->SetTitle(";ll#gamma#gamma p_{T} [GeV];Events");
            h->GetYaxis()->CenterTitle(true);
            h->SetAxisRange(0., 25,"Y");
//            h->SetAxisRange(0., 155,"Y");
            h->Draw("HIST");
        }
        else
        {
            h->Scale(factor/h->Integral());
            h->SetTitle(";ll#gamma#gamma p_{T} [GeV];Events");
            h->GetYaxis()->CenterTitle(true);
            h->SetAxisRange(0., 25,"Y");
//            h->SetAxisRange(0., 155,"Y");
            h->Draw("HISTsame");
            gPad->Modified(); gPad->Update();
        }
    }
    gStyle->SetOptStat(0);
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.6, 0.83, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.6, 0.73,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("Fig8C.png");
    
    c1 = new TCanvas();
    legend = new TLegend(0.65, 0.45, 0.85, 0.65);
    count = 0;
    for (auto& h: histos_deltaR)
    {
        h->SetLineColor(colors[count++]);
        legend->AddEntry(&(*h), h->GetTitle(), "l");
        
        if (&h == &histos_Z.front())
        {
            factor = h->Integral();
            h->Scale(factor/h->Integral());
            h->SetTitle(";#DeltaR ll#gamma#gamma;Events");
            h->GetYaxis()->CenterTitle(true);
            h->SetAxisRange(0., 20,"Y");
//            h->SetAxisRange(0., 130,"Y");
            h->Draw("HIST");
        }
        else
        {
            h->Scale(factor/h->Integral());
            h->SetTitle(";#DeltaR ll#gamma#gamma;Events");
            h->GetYaxis()->CenterTitle(true);
            h->SetAxisRange(0., 20,"Y");
//            h->SetAxisRange(0., 130,"Y");
            h->Draw("HISTsame");
            gPad->Modified(); gPad->Update();
        }
    }
    gStyle->SetOptStat(0);
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.6, 0.83, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.6, 0.73,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("Fig8D.png");
    
    c1 = new TCanvas();
    legend = new TLegend(0.65, 0.45, 0.85, 0.65);
    count = 0;
    for (auto& h: histos_deltaPhi)
    {
        h->SetLineColor(colors[count++]);
        legend->AddEntry(&(*h), h->GetTitle(), "l");
        
        if (&h == &histos_Z.front())
        {
            factor = h->Integral();
            h->Scale(factor/h->Integral());
            h->SetTitle(";#Delta#phi ll#gamma#gamma;Events");
            h->GetYaxis()->CenterTitle(true);
            h->SetAxisRange(0., 35,"Y");
//            h->SetAxisRange(0., 200,"Y");
            h->Draw("HIST");
        }
        else
        {
            h->Scale(factor/h->Integral());
            h->SetTitle(";#Delta#phi ll#gamma#gamma;Events");
            h->GetYaxis()->CenterTitle(true);
            h->SetAxisRange(0., 35,"Y");
//            h->SetAxisRange(0., 200,"Y");
            h->Draw("HISTsame");
            gPad->Modified(); gPad->Update();
        }
    }
    gStyle->SetOptStat(0);
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.6, 0.83, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.6, 0.73,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("Fig8E.png");
    
    c1 = new TCanvas();
    legend = new TLegend(0.65, 0.45, 0.85, 0.65);
    count = 0;
    for (auto& h: histos_deltaEta)
    {
        h->SetLineColor(colors[count++]);
        legend->AddEntry(&(*h), h->GetTitle(), "l");
        
        if (&h == &histos_Z.front())
        {
            factor = h->Integral();
            h->Scale(factor/h->Integral());
            h->SetTitle(";#Delta#eta ll#gamma#gamma;Events");
            h->GetYaxis()->CenterTitle(true);
            h->SetAxisRange(0., 20,"Y");
//            h->SetAxisRange(0., 130,"Y");
            h->Draw("HIST");
        }
        else
        {
            h->Scale(factor/h->Integral());
            h->SetTitle(";#Delta#eta ll#gamma#gamma;Events");
            h->GetYaxis()->CenterTitle(true);
            h->SetAxisRange(0., 20,"Y");
//            h->SetAxisRange(0., 130,"Y");
            h->Draw("HISTsame");
            gPad->Modified(); gPad->Update();
        }
    }
    gStyle->SetOptStat(0);
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.6, 0.83, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.6, 0.73,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("Fig8F.png");
}

void fig10()
{
    std::vector<std::string> input_filenames = {
        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root",
        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root"
    };
    
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
    
    std::vector<ROOT::RDF::RResultPtr<TH1D>> histos_deltaR_low;
    histos_deltaR_low.reserve(input_filenames.size());
    std::vector<ROOT::RDF::RResultPtr<TH1D>> histos_deltaR_high;
    histos_deltaR_high.reserve(input_filenames.size());

    std::vector<const char*> prefixes = {"sig ggF m_{A} = 5 GeV", "sig ggF m_{A} = 1 GeV"};
    std::vector<EColor> colors = {kBlack, kMagenta};
    int count = 0;
    
    for (auto& file: input_filenames)
    {
        SchottDataFrame df(MakeRDF({file}, 8));
        
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
                return ((abs(x.photon_eta) >= 2.37) || (x.photon_pt <= 10e3) || (abs(x.photon_eta) > 1.37 && abs(x.photon_eta) < 1.52) || (!x.photon_id_loose));

            }), photons.end());
            
            return photons;
        }, {"photons"});
        
        auto reco_photons_matched = truth_photons_from_axions.Define("reco_photons_matched",
        [&](RVec<Photon>& photons_pass_cuts, RVec<TruthParticle>& truth_photons_from_axions)
        {
            RVec<Photon> matchedPhotons;
            PtEtaPhiEVector tp1 = truth_photons_from_axions[0].Vector();
            PtEtaPhiEVector tp2 = truth_photons_from_axions[1].Vector();
            
            for (auto& rp: photons_pass_cuts)
            {
                if (!(DeltaR(rp.Vector(), tp1) < 0.1 || DeltaR(rp.Vector(), tp2) < 0.1))
                {
                    continue;
                }
                matchedPhotons.push_back(rp);
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
        
//        std::cout << *two_reco_photons_matched.Count() << '\n';
        
        auto stable_truth_dileptons_and_diphotons = two_reco_photons_matched
        .Define("stable_truth_leptons",[&](RVec<TruthParticle> truth_particles)
        {
            truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
            [](TruthParticle& x)
            {
                return (((abs(x.mc_pdg_id) != 11 && abs(x.mc_pdg_id) != 12 && abs(x.mc_pdg_id) != 13 &&
                        abs(x.mc_pdg_id) != 14 && abs(x.mc_pdg_id) != 15 && abs(x.mc_pdg_id) != 16 &&
                         abs(x.mc_pdg_id) != 17 && abs(x.mc_pdg_id) != 18)) || (x.mc_status != 1));
                
            }), truth_particles.end());

            return truth_particles;
            
        }, {"truth_particles"})
        .Filter([&](RVec<TruthParticle>& stable_truth_leptons, RVec<TruthParticle>& truth_photons_from_axions)
        {
            return (stable_truth_leptons.size()==2 && truth_photons_from_axions.size()==2);
            
        }, {"stable_truth_leptons", "truth_photons_from_axions"})
        .Define("delta_r_gamma_gamma",[&](RVec<TruthParticle>& truth_photons_from_axions)
        {
            return DeltaR(truth_photons_from_axions[0].Vector(),truth_photons_from_axions[1].Vector());
            
        }, {"truth_photons_from_axions"});
        
//        std::cout << *stable_truth_dileptons_and_diphotons.Count() << '\n';
        if (count==0) //add more later
        {
            histos_deltaR_high.push_back(stable_truth_dileptons_and_diphotons.Histo1D<double>({prefixes[count], prefixes[count], 100u, 0, 1.05}, "delta_r_gamma_gamma"));
        }
        else if (count==1) //add more later
        {
            histos_deltaR_low.push_back(stable_truth_dileptons_and_diphotons.Histo1D<double>({prefixes[count], prefixes[count], 100u, 0, 1.05}, "delta_r_gamma_gamma"));
        }
        count++;
        
        auto passed = stable_truth_dileptons_and_diphotons.Count();
        std::cout << *passed << '\n';
    }
    
    TCanvas* c1 = new TCanvas();
    TLegend* legend = new TLegend(0.6, 0.4, 0.8, 0.6);
    double factor;
    
    count = 0;
    for (auto& h: histos_deltaR_high)
    {
        h->SetLineColor(colors[count++]);
        legend->AddEntry(&(*h), h->GetTitle(), "l");
        
        if (&h == &histos_deltaR_high.front())
        {
            factor = h->Integral();
            h->Scale(factor/h->Integral());
            h->SetTitle(";#DeltaR_{#gamma#gamma};Events");
            h->GetYaxis()->CenterTitle(true);
//            h->SetAxisRange(0., 300,"Y");
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
    Tl.DrawLatexNDC(0.6, 0.83, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.6, 0.73,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("Fig10A.png");
    
    c1 = new TCanvas();
    legend = new TLegend(0.65, 0.45, 0.85, 0.65);
    count = 0;
    for (auto& h: histos_deltaR_low)
    {
        h->SetLineColor(colors[count++]);
        legend->AddEntry(&(*h), h->GetTitle(), "l");
        
        if (&h == &histos_deltaR_low.front())
        {
            factor = h->Integral();
            h->Scale(factor/h->Integral());
            h->SetTitle(";#DeltaR_{#gamma#gamma};Events");
            h->GetYaxis()->CenterTitle(true);
//            h->SetAxisRange(0., 300,"Y");
            h->Draw("HIST");
        }
        else
        {
            h->Scale(factor/h->Integral());
            h->SetTitle(";#gamma#gamma p_{T} [GeV];Events");
            h->GetYaxis()->CenterTitle(true);
//            h->SetAxisRange(0., 300,"Y");
            h->Draw("HISTsame");
            gPad->Modified(); gPad->Update();
        }
    }
    gStyle->SetOptStat(0);
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.6, 0.83, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.6, 0.73,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("Fig10B.png");
}

void fig18()
{
    std::vector<std::string> input_filenames = {
        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root",
        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root"
    };
    
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
    
    std::vector<ROOT::RDF::RResultPtr<TH1D>> histos_mA5;
    histos_mA5.reserve(input_filenames.size());
    std::vector<ROOT::RDF::RResultPtr<TH1D>> histos_mA1;
    histos_mA1.reserve(input_filenames.size());

    std::vector<const char*> prefixes = {"Both photons match", "one photon match", "None match"};
    std::vector<EColor> colors = {kBlack, kRed, kBlue};
    int count = 0;
    
    for (auto& file: input_filenames)
    {
        SchottDataFrame df(MakeRDF({file}, 8));
        
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
                return ((abs(x.photon_eta) >= 2.37) || (x.photon_pt <= 10e3) || (abs(x.photon_eta) > 1.37 && abs(x.photon_eta) < 1.52) || (!x.photon_id_loose));

            }), photons.end());
            
            return photons;
        }, {"photons"});
        
        auto reco_photons_matched = truth_photons_from_axions.Define("reco_photons_matched",
        [&](RVec<Photon>& photons_pass_cuts, RVec<TruthParticle>& truth_photons_from_axions)
        {
            RVec<Photon> matchedPhotons;
            PtEtaPhiEVector tp1 = truth_photons_from_axions[0].Vector();
            PtEtaPhiEVector tp2 = truth_photons_from_axions[1].Vector();
            
            for (auto& rp: photons_pass_cuts)
            {
                if (!(DeltaR(rp.Vector(), tp1) < 0.1 || DeltaR(rp.Vector(), tp2) < 0.1))
                {
                    continue;
                }
                matchedPhotons.push_back(rp);
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
            
        }, {"reco_photons_matched"})
        .Define("truth_photons_from_axions_pt",
        [&](RVec<TruthParticle>& truth_photons_from_axions)
        {
            PtEtaPhiEVector four_momentum = truth_photons_from_axions[0].Vector() + truth_photons_from_axions[1].Vector();
            return four_momentum.Pt()/1e3;
        }, {"truth_photons_from_axions"});
        
        auto one_reco_photons_matched = reco_photons_matched.Filter(
        [&](RVec<Photon>& reco_photons_matched)
        {
            return reco_photons_matched.size()==1;
            
        }, {"reco_photons_matched"})
        .Define("truth_photons_from_axions_pt",
        [&](RVec<TruthParticle>& truth_photons_from_axions)
        {
            PtEtaPhiEVector four_momentum = truth_photons_from_axions[0].Vector() + truth_photons_from_axions[1].Vector();
            return four_momentum.Pt()/1e3;
        }, {"truth_photons_from_axions"});
        
        auto zero_reco_photons_matched = reco_photons_matched.Filter(
        [&](RVec<Photon>& reco_photons_matched)
        {
            return reco_photons_matched.empty();
            
        }, {"reco_photons_matched"})
        .Define("truth_photons_from_axions_pt",
        [&](RVec<TruthParticle>& truth_photons_from_axions)
        {
            PtEtaPhiEVector four_momentum = truth_photons_from_axions[0].Vector() + truth_photons_from_axions[1].Vector();
            return four_momentum.Pt()/1e3;
        }, {"truth_photons_from_axions"});
        
//        std::cout << *stable_truth_dileptons_and_diphotons.Count() << '\n';
        if (count==0) //add more later
        {
            histos_mA5.push_back(two_reco_photons_matched.Histo1D<double>({prefixes[0], prefixes[0], 100u, 0, 202}, "truth_photons_from_axions_pt"));
            histos_mA5.push_back(one_reco_photons_matched.Histo1D<double>({prefixes[1], prefixes[1], 100u, 0, 202}, "truth_photons_from_axions_pt"));
            histos_mA5.push_back(zero_reco_photons_matched.Histo1D<double>({prefixes[2], prefixes[2], 100u, 0, 202}, "truth_photons_from_axions_pt"));
        }
        else if (count==1) //add more later
        {
            histos_mA1.push_back(two_reco_photons_matched.Histo1D<double>({prefixes[0], prefixes[0], 100u, 0, 202}, "truth_photons_from_axions_pt"));
            histos_mA1.push_back(one_reco_photons_matched.Histo1D<double>({prefixes[1], prefixes[1], 100u, 0, 202}, "truth_photons_from_axions_pt"));
            histos_mA1.push_back(zero_reco_photons_matched.Histo1D<double>({prefixes[2], prefixes[2], 100u, 0, 202}, "truth_photons_from_axions_pt"));
        }
        count++;
        
        auto passed = preselection.Count();
        std::cout << *passed << '\n';
    }
    
    TCanvas* c1 = new TCanvas();
    TLegend* legend = new TLegend(0.55, 0.3, 0.75, 0.5);
    double factor;
    
    count = 0;
    for (auto& h: histos_mA5)
    {
        h->SetLineColor(colors[count++]);
        legend->AddEntry(&(*h), h->GetTitle(), "l");
        
        if (&h == &histos_mA5.front())
        {
            factor = h->Integral();
//            h->Scale(factor/h->Integral());
            h->SetTitle(";p_{T_{#gamma#gamma}} [GeV];Events");
            h->GetYaxis()->CenterTitle(true);
            h->GetXaxis()->SetTitleOffset(1.2);
            h->SetAxisRange(0., 170,"Y");
//            h->SetAxisRange(0., 1250,"Y");
            h->Draw("HIST");
        }
        else
        {
//            h->Scale(factor/h->Integral());
            h->Draw("HISTsame");
            gPad->Modified(); gPad->Update();
        }
    }
    
    gStyle->SetOptStat(0);
    TLatex Tl;
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.6, 0.8, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.6, 0.7,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
    Tl.DrawLatexNDC(0.6, 0.6,"m_{A} = 5 GeV, #DeltaR = 0.1");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("Fig18A.png");
    
    c1 = new TCanvas();
    legend = new TLegend(0.55, 0.3, 0.75, 0.5);
    count = 0;
    for (auto& h: histos_mA1)
    {
        h->SetLineColor(colors[count++]);
        legend->AddEntry(&(*h), h->GetTitle(), "l");
        
        if (&h == &histos_mA1.front())
        {
            factor = h->Integral();
//            h->Scale(factor/h->Integral());
            h->SetTitle(";p_{T_{#gamma#gamma}} [GeV];Events");
            h->GetYaxis()->CenterTitle(true);
            h->GetXaxis()->SetTitleOffset(1.2);
            h->SetAxisRange(0., 430,"Y");
            h->Draw("HIST");
        }
        else
        {
//            h->Scale(factor/h->Integral());
            h->SetTitle(";p_{T_{#gamma#gamma}} [GeV];Events");
            h->GetYaxis()->CenterTitle(true);
            h->GetXaxis()->SetTitleOffset(1.2);
            h->SetAxisRange(0., 430,"Y");
            h->Draw("HISTsame");
            gPad->Modified(); gPad->Update();
        }
    }
    gStyle->SetOptStat(0);
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.6, 0.8, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.6, 0.7,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
    Tl.DrawLatexNDC(0.6, 0.6,"m_{A} = 1 GeV, #DeltaR = 0.1");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("Fig18B.png");
}

void fig24()
{
    std::vector<std::string> input_filenames = {
        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root",
        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root",
    };
    
    std::vector<std::string> back_input_filenames = {
        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617070._000001.LGNTuple.root", "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617064._000001.LGNTuple.root",
        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617074._000001.LGNTuple.root",
    };
    
    TLegend* legend = new TLegend(0.6, 0.4, 0.8, 0.6);
    
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
    
    std::vector<ROOT::RDF::RResultPtr<TH1D>> histos_unweighted;
    histos_unweighted.reserve(input_filenames.size());
    std::vector<ROOT::RDF::RResultPtr<TH1D>> histos_weighted;
    histos_weighted.reserve(input_filenames.size());
    
    std::vector<ROOT::RDF::RResultPtr<ULong64_t>> backCounts;
    std::vector<ROOT::RDF::RResultPtr<TH1D>> back_histos_unweighted;
    back_histos_unweighted.reserve(back_input_filenames.size());
    std::vector<ROOT::RDF::RResultPtr<TH1D>> back_histos_weighted;
    back_histos_weighted.reserve(back_input_filenames.size());
    
    std::vector<const char*> prefixes = {"sig m_{A} = 5 GeV", "sig m_{A} = 1 GeV"};
    std::vector<const char*> back_prefixes = {"pty2_9_17", "pty_17_myy_0_80", "pty_17_myy_80"};
    
    std::array<double,3> SFs = {((139e15)*(.871e-12))/150000.,((139e15)*(.199e-12))/150000., ((139e15)*(.0345e-15))/110465.};
    
    std::vector<EColor> sig_colors = {kBlack, kMagenta};
    int count = 0;
    
    TCanvas* c1 = new TCanvas();
    legend = new TLegend(0.65, 0.4, 0.85, 0.6);
    
    for (auto& file: input_filenames)
    {
        SchottDataFrame df(MakeRDF({file}, 8));
        
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
                return ((abs(x.photon_eta) >= 2.37) || (x.photon_pt <= 10e3) || (abs(x.photon_eta) > 1.37 && abs(x.photon_eta) < 1.52) || (!x.photon_id_loose));

            }), photons.end());
            
            return photons;
        }, {"photons"});
        
        auto reco_photons_matched = truth_photons_from_axions.Define("reco_photons_matched",
        [&](RVec<Photon>& photons_pass_cuts, RVec<TruthParticle>& truth_photons_from_axions)
        {
            RVec<Photon> matchedPhotons;
            PtEtaPhiEVector tp1 = truth_photons_from_axions[0].Vector();
            PtEtaPhiEVector tp2 = truth_photons_from_axions[1].Vector();
            
            for (auto& rp: photons_pass_cuts)
            {
                if (!(DeltaR(rp.Vector(), tp1) < 0.1 || DeltaR(rp.Vector(), tp2) < 0.1))
                {
                    continue;
                }
                matchedPhotons.push_back(rp);
            }
            
            std::sort(matchedPhotons.begin(),matchedPhotons.end(),
            [](Photon& x, Photon& y)
            {
                return x.photon_pt > y.photon_pt;
            });
            
            return matchedPhotons;
            
        }, {"photons_pass_cuts", "truth_photons_from_axions"});
        
        auto two_reco_photons_matched = reco_photons_matched.Define("chosen_two",
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
        }, {"reco_photons_matched"}).Filter(
        [&](RVec<Photon>& reco_photons_matched)
        {
            return (reco_photons_matched.size()==2);
        }, {"chosen_two"})
        .Define("reco_photons_from_axions_deltaR",
        [&](RVec<Photon>& reco_photons_matched)
        {
            return DeltaR(reco_photons_matched[0].Vector(), reco_photons_matched[1].Vector());
            
        }, {"chosen_two"});

        histos_unweighted.push_back(two_reco_photons_matched.Histo1D<double>({prefixes[count], prefixes[count], 100u, 0, 0.8}, "reco_photons_from_axions_deltaR"));
        histos_weighted.push_back(two_reco_photons_matched.Histo1D<double>({prefixes[count], prefixes[count++], 100u, 0, 0.8}, "reco_photons_from_axions_deltaR"));
        
//        auto passed = preselection.Count();
//        std::cout << *passed << '\n';
    }
    
    count = 0;
    std::vector<EColor> back_colors = {kBlue, kRed, kViolet,};
    
    for (auto& file: back_input_filenames)
    {
        SchottDataFrame df(MakeRDF({file}, 8));
        
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
            
        auto two_leptons = trigger_selection.Filter(
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
        
        auto photons_pass_cuts = ptCut.Define("photons_pass_cuts",
        [&](RVec<Photon> photons)
        {
            photons.erase(std::remove_if(photons.begin(),photons.end(),
            [](Photon& x)
            {
                return ((abs(x.photon_eta) >= 2.37) || (x.photon_pt <= 10e3) || (abs(x.photon_eta) > 1.37 && abs(x.photon_eta) < 1.52) || (!x.photon_id_loose));

            }), photons.end());
            
            return photons;
        }, {"photons"});
        
        
        auto two_reco_photons_matched = photons_pass_cuts.Define("chosen_two",
        [](RVec<Photon>& photons_pass_cuts)
        {
            RVec<Photon> x;
            if (photons_pass_cuts.size() < 2)
            {
                return x;
            }
            auto combs = Combinations(photons_pass_cuts, 2);
            size_t length = combs[0].size();
            double delta_r, m, pt, X;

            for (size_t i=0; i<length; i++)
            {
                delta_r = DeltaR(photons_pass_cuts[combs[0][i]].Vector(), photons_pass_cuts[combs[1][i]].Vector());
                m = (photons_pass_cuts[combs[0][i]].Vector() + photons_pass_cuts[combs[1][i]].Vector()).M();
                pt = (photons_pass_cuts[combs[0][i]].Vector() + photons_pass_cuts[combs[1][i]].Vector()).Pt();
                X = delta_r*(pt/(2.0*m));
                if ((delta_r < 1.5) && (X > 0.96) && (X < 1.2))
                {
                    x = {photons_pass_cuts[combs[0][i]], photons_pass_cuts[combs[1][i]]};
                    return x;
                }
            }
            return x;
        }, {"photons_pass_cuts"}).Filter(
        [&](RVec<Photon>& two_reco_photons_matched)
        {
            return (two_reco_photons_matched.size()==2);
        }, {"chosen_two"})
        .Define("reco_photons_from_axions_deltaR",
        [&](RVec<Photon>& reco_photons_matched)
        {
            return DeltaR(reco_photons_matched[0].Vector(), reco_photons_matched[1].Vector());
            
        }, {"chosen_two"});

        back_histos_unweighted.push_back(two_reco_photons_matched.Histo1D<double>({back_prefixes[count], back_prefixes[count], 100u, 0, 0.8}, "reco_photons_from_axions_deltaR"));
        back_histos_weighted.push_back(two_reco_photons_matched.Histo1D<double>({back_prefixes[count], back_prefixes[count++], 100u, 0, 0.8}, "reco_photons_from_axions_deltaR"));
        
        backCounts.push_back(two_reco_photons_matched.Count());
    }
    
    double factor;
    double total_back = 0;
    for (size_t i = 0; i < backCounts.size(); i++)
    {
        total_back += *backCounts[i]*SFs[i];
    }
    int back_count = 0;
    factor = total_back;
    auto hs = new THStack("hs3","");
    for (auto& h: back_histos_weighted)
    {
        h->SetFillColor(back_colors[back_count]);
        legend->AddEntry(&(*h), h->GetTitle(), "f");
        if (h->Integral() != 0)
            h->Scale((factor/h->Integral())*SFs[back_count++]);
        h->SetTitle(";m_{ll}  [GeV];Events");
        h->GetYaxis()->CenterTitle(true);
        h->SetAxisRange(0., 1000,"Y");
//        gPad->Modified(); gPad->Update();
        
        hs->Add(&*h);
        h->SetAxisRange(0., 1000,"Y");
        count++;
    }
    hs->Draw("HIST");
    
    count = 0;
    for (auto& h: histos_weighted)
    {
        h->SetLineColor(sig_colors[count++]);
        legend->AddEntry(&(*h), h->GetTitle(), "l");
        
        if (&h == &histos_weighted.front())
        {
            if (h->Integral() != 0)
            {h->Scale(factor/h->Integral());}
            h->SetTitle(";#DeltaR_{#gamma#gamma} [GeV];Events");
            h->GetYaxis()->CenterTitle(true);
            h->GetXaxis()->SetTitleOffset(1.2);
            h->SetAxisRange(0., 82,"Y");
//            h->SetAxisRange(0., 500,"Y");

            h->Draw("HISTsame");
        }
        else
        {
            if (h->Integral() != 0)
            {h->Scale(factor/h->Integral());}
            h->Draw("HISTsame");
            gPad->Modified(); gPad->Update();
        }
    }
    
    hs->SetTitle(";#DeltaR_{#gamma#gamma} [GeV];Events");
    hs->GetYaxis()->CenterTitle(true);
    hs->GetXaxis()->SetTitleOffset(1.2);
    hs->SetMinimum(0.);
    hs->SetMaximum(29.);
    
    gStyle->SetOptStat(0);
    TLatex Tl;
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.6, 0.8, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.6, 0.7,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("Fig24B.png");
    
    c1 = new TCanvas();
    legend = new TLegend(0.65, 0.4, 0.85, 0.6);
    
    back_count=0;
    hs = new THStack("hs3","");
    for (auto& h: back_histos_unweighted)
    {
        h->SetFillColor(back_colors[back_count]);
        legend->AddEntry(&(*h), h->GetTitle(), "f");
        if (h->Integral() != 0)
            h->Scale((factor/h->Integral())*SFs[back_count++]);
        h->SetTitle(";m_{ll}  [GeV];Events");
        h->GetYaxis()->CenterTitle(true);
        h->SetAxisRange(0., 1000,"Y");
//        gPad->Modified(); gPad->Update();
        
        hs->Add(&*h);
        h->SetAxisRange(0., 1000,"Y");
        count++;
    }
    hs->Draw("HIST");
    count = 0;
    for (auto& h: histos_unweighted)
    {
        h->SetLineColor(sig_colors[count++]);
        legend->AddEntry(&(*h), h->GetTitle(), "l");
        
        if (&h == &histos_unweighted.front())
        {
            h->SetTitle(";#DeltaR_{#gamma#gamma} [GeV];Events");
            h->GetYaxis()->CenterTitle(true);
            h->GetXaxis()->SetTitleOffset(1.2);
            h->SetAxisRange(0., 36,"Y");
            h->Draw("HISTsame");
        }
        else
        {
            h->SetTitle(";#DeltaR_{#gamma#gamma} [GeV];Events");
            h->GetYaxis()->CenterTitle(true);
            h->GetXaxis()->SetTitleOffset(1.2);
            h->SetAxisRange(0., 36,"Y");
            h->Draw("HISTsame");
            gPad->Modified(); gPad->Update();
        }
    }
    hs->SetTitle(";#DeltaR_{#gamma#gamma} [GeV];Events");
    hs->GetYaxis()->CenterTitle(true);
    hs->GetXaxis()->SetTitleOffset(1.2);
    hs->SetMinimum(0.);
    hs->SetMaximum(36.);
    
    gStyle->SetOptStat(0);
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.6, 0.8, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.6, 0.7,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("Fig24A.png");
}

void fig54()
{
    std::vector<std::string> input_filenames = {
        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root",
        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root"
    };
    
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
    
    std::vector<ROOT::RDF::RResultPtr<TH1D>> histos_mass;
    histos_mass.reserve(input_filenames.size());
    std::vector<ROOT::RDF::RResultPtr<TH1D>> histos_deltaR;
    histos_deltaR.reserve(input_filenames.size());
    std::vector<ROOT::RDF::RResultPtr<TH1D>> histos_pt;
    histos_pt.reserve(input_filenames.size());
    std::vector<ROOT::RDF::RResultPtr<TH1D>> histos_E_Ratio;
    histos_E_Ratio.reserve(input_filenames.size());

    std::vector<const char*> prefixes = {"Sig m_{A} = 5 GeV", "Sig m_{A} = 1 GeV"};
    std::vector<EColor> colors = {kBlue, kRed};
    int count = 0;
    
    for (auto& file: input_filenames)
    {
        SchottDataFrame df(MakeRDF({file}, 8));
        
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
                return ((abs(x.photon_eta) >= 2.37) || (abs(x.photon_eta) > 1.37 && abs(x.photon_eta) < 1.52) || (!x.photon_id_loose));

            }), photons.end());
            
            return photons;
        }, {"photons"});
        
        auto reco_photons_matched = truth_photons_from_axions.Define("reco_photons_matched",
        [&](RVec<Photon>& photons_pass_cuts, RVec<TruthParticle>& truth_photons_from_axions)
        {
            RVec<Photon> matchedPhotons;
            PtEtaPhiEVector tp1 = truth_photons_from_axions[0].Vector();
            PtEtaPhiEVector tp2 = truth_photons_from_axions[1].Vector();
            
            for (auto& rp: photons_pass_cuts)
            {
                if (!(DeltaR(rp.Vector(), tp1) < 0.1 || DeltaR(rp.Vector(), tp2) < 0.1))
                {
                    continue;
                }
                matchedPhotons.push_back(rp);
            }
            
            std::sort(matchedPhotons.begin(),matchedPhotons.end(),
            [](Photon& x, Photon& y)
            {
                return x.photon_pt > y.photon_pt;
            });
            
            return matchedPhotons;
            
        }, {"photons_pass_cuts", "truth_photons_from_axions"});
        
        auto merged_reco_photons_matched = reco_photons_matched.Filter(
        [&](RVec<Photon>& reco_photons_test)
        {
            RVec<Photon> reco_photons_matched = reco_photons_test;
            
            reco_photons_matched.erase(std::remove_if(reco_photons_matched.begin(),reco_photons_matched.end(),
            [](Photon& x)
            {
                return x.photon_pt <= 10e3;

            }), reco_photons_matched.end());
            
            double delta_r = DeltaR(reco_photons_matched[0].Vector(), reco_photons_matched[1].Vector());
            double m = (reco_photons_matched[0].Vector() + reco_photons_matched[1].Vector()).M();
            double pt = (reco_photons_matched[0].Vector() + reco_photons_matched[1].Vector()).Pt();
            double X = delta_r*(pt/(2.0*m));
            
            if ((reco_photons_matched.size()==2)
            && (delta_r < 1.5)
            && (X > 0.96)
            && (X < 1.2))
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
            
        }, {"reco_photons_matched"})
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
            
        }, {"reco_photons_matched"});
        
        auto stable_truth_dilepton_and_photon = merged_reco_photons_matched
        .Define("stable_truth_leptons",[&](RVec<TruthParticle> truth_particles)
        {
            truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
            [](TruthParticle& x)
            {
                return (((abs(x.mc_pdg_id) != 11 && abs(x.mc_pdg_id) != 12 && abs(x.mc_pdg_id) != 13 &&
                        abs(x.mc_pdg_id) != 14 && abs(x.mc_pdg_id) != 15 && abs(x.mc_pdg_id) != 16 &&
                         abs(x.mc_pdg_id) != 17 && abs(x.mc_pdg_id) != 18)) || (x.mc_status != 1));
                
            }), truth_particles.end());

            return truth_particles;
            
        }, {"truth_particles"})
        .Filter([&](RVec<TruthParticle>& stable_truth_leptons, RVec<TruthParticle>& truth_photons_from_axions)
        {
            return (stable_truth_leptons.size()==2 && truth_photons_from_axions.size()==2);
            
        }, {"stable_truth_leptons", "truth_photons_from_axions"})
        .Define("reconstructed_mass",[&](RVec<TruthParticle>& stable_truth_leptons, Photon& merged_photon)
        {
            auto four_momentum = stable_truth_leptons[0].Vector() + stable_truth_leptons[1].Vector();
            
            return (four_momentum + merged_photon.Vector()).M()/1e3;
            
        }, {"stable_truth_leptons", "merged_photon"})
        .Define("reconstructed_deltaR",[&](RVec<TruthParticle>& stable_truth_leptons, Photon& merged_photon)
        {
            auto four_momentum = stable_truth_leptons[0].Vector() + stable_truth_leptons[1].Vector();
            
            return DeltaR(four_momentum,merged_photon.Vector());
            
        }, {"stable_truth_leptons", "merged_photon"})
        .Define("reconstructed_pt",[&](RVec<TruthParticle>& stable_truth_leptons, Photon& merged_photon)
        {
            auto four_momentum = stable_truth_leptons[0].Vector() + stable_truth_leptons[1].Vector();
            
            return (four_momentum + merged_photon.Vector()).Pt()/1e3;
            
        }, {"stable_truth_leptons", "merged_photon"});
                                            
//        std::cout << *stable_truth_dileptons_and_diphotons.Count() << '\n';

        histos_mass.push_back(stable_truth_dilepton_and_photon.Histo1D<double>({prefixes[count], prefixes[count], 100u, 80, 200}, "reconstructed_mass"));
        histos_deltaR.push_back(stable_truth_dilepton_and_photon.Histo1D<double>({prefixes[count], prefixes[count], 100u, 0, 4.7}, "reconstructed_deltaR"));
        histos_pt.push_back(stable_truth_dilepton_and_photon.Histo1D<double>({prefixes[count], prefixes[count], 100u, 0, 220}, "reconstructed_pt"));
        histos_E_Ratio.push_back(stable_truth_dilepton_and_photon.Histo1D<RVec<float>>({prefixes[count], prefixes[count++], 100u, 0, 1}, "photon_shower_shape_e_ratio"));
        
        auto passed = preselection.Count();
        std::cout << *passed << '\n';
    }
    
    TCanvas* c1 = new TCanvas();
    TLegend* legend = new TLegend(0.65, 0.4, 0.85, 0.6);
    double factor;
    
    count = 0;
    for (auto& h: histos_mass)
    {
        h->SetLineColor(colors[count++]);
        legend->AddEntry(&(*h), h->GetTitle(), "l");
        
        if (&h == &histos_mass.front())
        {
            factor = h->Integral();
            h->Scale(factor/h->Integral());
            h->SetTitle(";m_{ll#gamma} [GeV];Events");
            h->GetYaxis()->CenterTitle(true);
            h->GetXaxis()->SetTitleOffset(1.2);
            h->SetAxisRange(0., 155,"Y");
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
    c1->SaveAs("Fig54A.png");
    
    c1 = new TCanvas();
    legend = new TLegend(0.15, 0.5, 0.35, 0.7);
    count = 0;
    for (auto& h: histos_deltaR)
    {
        h->SetLineColor(colors[count++]);
        legend->AddEntry(&(*h), h->GetTitle(), "l");
        
        if (&h == &histos_deltaR.front())
        {
            factor = h->Integral();
            h->Scale(factor/h->Integral());
            h->SetTitle(";#DeltaR (ll,#gamma) ;Events");
            h->GetYaxis()->CenterTitle(true);
            h->GetXaxis()->SetTitleOffset(1.2);
            h->SetAxisRange(0., 29,"Y");
//            h->SetAxisRange(0., 210,"Y");
            h->Draw("HIST");
        }
        else
        {
            h->Scale(factor/h->Integral());
            h->SetTitle(";#DeltaR (ll,#gamma) ;Events");
            h->GetYaxis()->CenterTitle(true);
            h->GetXaxis()->SetTitleOffset(1.2);
            h->SetAxisRange(0., 29,"Y");
//            h->SetAxisRange(0., 210,"Y");
            h->Draw("HISTsame");
            gPad->Modified(); gPad->Update();
        }
    }
    gStyle->SetOptStat(0);
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.15, 0.85, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.15, 0.75,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("Fig54C.png");
    
    c1 = new TCanvas();
    legend = new TLegend(0.65, 0.4, 0.85, 0.6);
    count = 0;
    for (auto& h: histos_pt)
    {
        h->SetLineColor(colors[count++]);
        legend->AddEntry(&(*h), h->GetTitle(), "l");
        
        if (&h == &histos_pt.front())
        {
            factor = h->Integral();
            h->Scale(factor/h->Integral());
            h->SetTitle(";p_{T} (ll,#gamma) ;Events");
            h->GetYaxis()->CenterTitle(true);
            h->GetXaxis()->SetTitleOffset(1.2);
//            h->SetAxisRange(0., 450,"Y");
            h->Draw("HIST");
        }
        else
        {
            h->Scale(factor/h->Integral());
            h->SetTitle(";p_{T} (ll,#gamma) ;Events");
            h->GetYaxis()->CenterTitle(true);
            h->GetXaxis()->SetTitleOffset(1.2);
//            h->SetAxisRange(0., 450,"Y");
            h->Draw("HISTsame");
            gPad->Modified(); gPad->Update();
        }
    }
    gStyle->SetOptStat(0);
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.6, 0.8, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.6, 0.7,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("Fig54D.png");
    
    c1 = new TCanvas();
    legend = new TLegend(0.45, 0.4, 0.65, 0.6);
    count = 0;
    for (auto& h: histos_E_Ratio)
    {
        h->SetLineColor(colors[count++]);
        legend->AddEntry(&(*h), h->GetTitle(), "l");
        
        if (&h == &histos_E_Ratio.front())
        {
            factor = h->Integral();
            h->Scale(factor/h->Integral());
            h->SetTitle(";E_{ratio};Events");
            h->GetYaxis()->CenterTitle(true);
            h->GetXaxis()->SetTitleOffset(1.2);
//            h->SetAxisRange(0., 450,"Y");
            h->Draw("HIST");
        }
        else
        {
            h->Scale(factor/h->Integral());
            h->SetTitle(";E_{ratio};Events");
            h->GetYaxis()->CenterTitle(true);
            h->GetXaxis()->SetTitleOffset(1.2);
//            h->SetAxisRange(0., 450,"Y");
            h->Draw("HISTsame");
            gPad->Modified(); gPad->Update();
        }
    }
    gStyle->SetOptStat(0);
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.4, 0.8, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.4, 0.7,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("Fig54B.png");
}


void ControlPlotsSignalShapes()
{
    auto start_time = Clock::now();
//    fig1A();
//    fig5();
//    fig6();
//    fig8();
//    fig10();
//    fig18();
    fig24();
//    fig54();
    
    auto end_time = Clock::now();
    std::cout << "Time difference: "
       << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count()/1e9 << " nanoseconds" << std::endl;
    
}


int main()
{
    ControlPlotsSignalShapes();
}


