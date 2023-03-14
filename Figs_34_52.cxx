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

#include "RDFObjects.h"
#include "MakeRDF.h"
#include "RDFevent.h"

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
void Fig34()
{
    std::vector<std::vector<std::string>> input_filenames = {
        {"/home/common/Za/NTuples/Ntuple_data_test.root"}
    };
    
    SchottDataFrame df(MakeRDF(input_filenames[0], 8));
    
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
            return ((abs(x.photon_eta) >= 2.37) || (abs(x.photon_eta) > 1.37 && abs(x.photon_eta) < 1.52) || (!x.photon_id_loose));

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
    
}*/

void Fig52()
{
    std::vector<std::vector<std::string>> input_filenames = {
        //Z gamma background
        {"/home/common/Za/NTuples/Background/user.kschmied.364860.eegammagamma_pty2_9_17.deriv.DAOD_STDM3.e7057_s3126_r10724_p4092_LGNTuple.root/user.kschmied.31617070._000001.LGNTuple.root"},
        {"/home/common/Za/NTuples/Background/user.kschmied.364861.eegammagamma_pty_17_myy_0_80.deriv.DAOD_HIGG1D2.e7057_s3126_r10724_p4062_LGNTuple.root/user.kschmied.31617064._000001.LGNTuple.root"},
        {"/home/common/Za/NTuples/Background/user.kschmied.364862.eegammagamma_pty_17_myy_80.deriv.DAOD_HIGG1D2.e7057_s3126_r10724_p4204_LGNTuple.root/user.kschmied.31660711._000001.LGNTuple.root", "/home/common/Za/NTuples/Background/user.kschmied.364862.eegammagamma_pty_17_myy_80.deriv.DAOD_HIGG1D2.e7057_s3126_r10724_p4204_LGNTuple.root/user.kschmied.31660711._000002.LGNTuple.root"},
        //Data
        {"/home/common/Za/NTuples/Ntuple_data_test.root"},
        //Jets
        {
            "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000001.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000002.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000003.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000004.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000005.LGNTuple.root",
        },
        {
            "/home/common/Za/NTuples/Background/user.kschmied.364117.v11.Zee_70_140_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835715._000001.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364117.v11.Zee_70_140_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835715._000002.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364117.v11.Zee_70_140_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835715._000003.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364117.v11.Zee_70_140_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835715._000004.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364117.v11.Zee_70_140_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835715._000006.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364117.v11.Zee_70_140_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835715._000007.LGNTuple.root",
        },
        {
            "/home/common/Za/NTuples/Background/user.kschmied.364120.v11.Zee_140_280_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835749._000001.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364120.v11.Zee_140_280_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835749._000002.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364120.v11.Zee_140_280_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835749._000003.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364120.v11.Zee_140_280_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835749._000004.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364120.v11.Zee_140_280_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835749._000005.LGNTuple.root",
        },
        {
            "/home/common/Za/NTuples/Background/user.kschmied.364115.v11.Zee_0_70_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835769._000001.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364115.v11.Zee_0_70_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835769._000002.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364115.v11.Zee_0_70_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835769._000003.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364115.v11.Zee_0_70_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835769._000004.LGNTuple.root",
        },
        {
            "/home/common/Za/NTuples/Background/user.kschmied.364118.v11.Zee_70_140_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835773._000001.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364118.v11.Zee_70_140_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835773._000002.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364118.v11.Zee_70_140_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835773._000003.LGNTuple.root",
        },
        {
            "/home/common/Za/NTuples/Background/user.kschmied.364121.v11.Zee_140_280_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835802._000001.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364121.v11.Zee_140_280_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835802._000002.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364121.v11.Zee_140_280_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835802._000003.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364121.v11.Zee_140_280_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835802._000004.LGNTuple.root",
        },
        {
            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000001.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000002.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000003.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000004.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000005.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000006.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000007.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000008.LGNTuple.root",
        },
        {
            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000001.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000002.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000003.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000004.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000005.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000006.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000007.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000008.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000009.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000010.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000011.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000012.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000013.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000014.LGNTuple.root",
        },
        {
            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000001.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000002.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000003.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000004.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000005.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000006.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000007.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000008.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000009.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000010.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000011.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000012.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000013.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000014.LGNTuple.root",
        },
    };
    
    std::array<double,3> SFs = {((139e15)*(.871e-12)),((139e15)*(.199e-12)), ((139e15)*(.0345e-15))}; //numerators for Z-gamma bkg
    
    std::array<double,9> JetNumeratorSFs = {((139e15)*(1.9828e-9)*(0.821204)),((139e15)*(110.64e-12)*(0.69275)),((139e15)*(40.645e-12)*(0.615906)),((139e15)*(1.9817e-9)*(0.1136684)),((139e15)*(110.47e-12)*(0.1912956)),((139e15)*(40.674e-12)*(0.2326772)),((139e15)*(1.9819e-9)*(0.0656969)),((139e15)*(110.53e-12)*(0.1158741)),((139e15)*(40.68e-12)*(0.1535215))}; //numerators for jet bkg
    
    std::vector<const char*> prefixes = {"pty2_9_17", "pty_17_myy_0_80", "pty_17_myy_80", "data", "Zee_lightJet_0-70", "Zee_lightJet_70-140", "Zee_lightJet_140-280", "Zee_cJet_0-70", "Zee_cJet_70-140", "Zee_cJet_140-280", "Zee_bJet_0-70", "Zee_bJet_70-140", "Zee_bJet_140-280"};
    std::vector<EColor> colors = {kBlue, kRed, kViolet, kGreen};
    std::vector<EColor> Jetscolors = {kCyan, kOrange, kGreen, kYellow, kPink, kGray, kBlack, kSpring, kAzure};

    std::vector<ROOT::RDF::RResultHandle> Nodes;
    int count = 0;
    
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
        [&](RVec<Photon> photons)
        {
            photons.erase(std::remove_if(photons.begin(),photons.end(),
            [](Photon& x)
            {
              return ((abs(x.photon_eta) >= 2.37) || (abs(x.photon_eta) > 1.37 && abs(x.photon_eta) < 1.52) || (!x.photon_id_loose));

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

        if (count <= 2 || count >= 4)
        {
            Nodes.push_back(merged_reco_photons_matched.Count());
            Nodes.push_back(df.Count());
        }
        
        Nodes.push_back(merged_reco_photons_matched.Histo1D<double>({prefixes[count], prefixes[count++], 60u, 0, 165}, "merged_photon_pt"));
    }
    
//    0   1   2   Z-gamma
//    3   4   5   Z-gamma
//    6   7   8   Z-gamma
//            9   Data
//    10  11  12  Z-jets
//    13  14  15  Z-jets
//    16  17  18  Z-jets
//    19  20  21  Z-jets
//    22  23  24  Z-jets
//    25  26  27  Z-jets
//    28  29  30  Z-jets
//    31  32  33  Z-jets
//    34  35  36  Z-jets
    
    ROOT::RDF::RunGraphs(Nodes); // running all computation nodes concurrently
    
    double factor;
    count = 0;
    for (auto& i: {0,3,6}) //Z-gamma
    {
        factor += (*Nodes[i].GetResultPtr<ULong64_t>())*(SFs[count++] / *Nodes[i+1].GetResultPtr<ULong64_t>());
    }
    
    for (int i = 10, j = 0; (i <= 34 && j <= 8); i += 3, j++) //Z-jets
    {
        factor += (*Nodes[i].GetResultPtr<ULong64_t>())*(JetNumeratorSFs[j] / *Nodes[i+1].GetResultPtr<ULong64_t>());
    }
    
    auto hs = new THStack("hs3","");
    TCanvas* c1 = new TCanvas();
    TLegend* legend = new TLegend(0.55, 0.2, 0.85, 0.6);
    count=0;
    for (auto& i: {2,5,8,9}) //Z-gamma
    {
        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0 && i != 9)
        {
            Nodes[i].GetResultPtr<TH1D>()->Scale(SFs[count] / *Nodes[i-1].GetResultPtr<ULong64_t>());
        }
        Nodes[i].GetResultPtr<TH1D>()->SetFillColor(colors[count++]);
        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "f");
    
        if (i != 9)
        {
            hs->Add(&*Nodes[i].GetResultPtr<TH1D>());
        }
    }
    
    for (int i = 12, j = 0; (i <= 36 && j <= 8); i += 3, j++) //Z-jets
    {
        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0)
        {
            Nodes[i].GetResultPtr<TH1D>()->Scale(JetNumeratorSFs[j] / *Nodes[i-1].GetResultPtr<ULong64_t>());
        }
        Nodes[i].GetResultPtr<TH1D>()->SetFillColor(Jetscolors[j]);
        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "f");
        hs->Add(&*Nodes[i].GetResultPtr<TH1D>());
    }

    hs->Draw("HIST");
    Nodes[9].GetResultPtr<TH1D>()->Draw("HISTsame");
    hs->SetTitle(";photon p_{T} [GeV];Events");
    hs->GetYaxis()->CenterTitle(true);
    hs->GetXaxis()->SetTitleOffset(1.2);
    hs->GetYaxis()->SetTitleOffset(1.35);
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
//    Fig34();
    Fig52();
    auto end_time = Clock::now();
    std::cout << "Time difference: "
       << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count()/1e9 << " nanoseconds" << std::endl;
}

int main()
{
    Figs_34_52();
}
