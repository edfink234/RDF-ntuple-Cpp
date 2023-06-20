#include <iostream>
#include <vector>
#include <chrono>
#include <algorithm>
#include <string>
#include <array>
#include <cstdlib>
#include <iomanip>
#include <map>
#include <set>

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
#include "ROOT/RDFHelpers.hxx"

#include "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/RDFObjects.h"
#include "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/MakeRDF.h"
#include "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/RDFevent.h"

using namespace ROOT::VecOps; // RVec
using namespace ROOT::Math::VectorUtil; // DeltaR
using namespace ROOT::Math; // PtEtaPhiEVector
using ROOT::RDF::Experimental::RResultMap;
using ROOT::RDF::Experimental::VariationsFor;

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

float roundToOneDecimalPlace(float num)
{
    std::stringstream stream;
    stream << std::fixed << std::setprecision(1) << num;
    float rounded_num = std::stof(stream.str());
    return rounded_num;
}

//Other mA5 file: /Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600909.PhPy8EG_AZNLO_ggH125_mA5p0_Cyy0p01_Czh1p0.merge.AOD.e8324_e7400_s3126_r10724_r10726_v2.root


void fig1A()
{
    std::vector<std::vector<std::string>> input_filenames =
    {
        //Z gamma background
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617070._000001.LGNTuple.root"}, //pty2_9_17
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617064._000001.LGNTuple.root"}, //pty_17_myy_0_80
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617074._000001.LGNTuple.root"}, //pty_17_myy_80
        //Signal
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root"}, //5 GeV
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root"}, // 1 GeV
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600751.PhPy8EG_AZNLO_ggH125_mA2p0_v1.root"}, // 2 GeV
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600752.PhPy8EG_AZNLO_ggH125_mA3p0_v1.root"}, // 3 GeV
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600756.PhPy8EG_AZNLO_ggH125_mA9p0_v1.root"}, // 9 GeV
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Test_oneMassPoint.root"},
        //Jets
        {"/Users/edwardfinkelstein/ATLAS_axion/Jets/Zee_lightJet_0-70.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/Jets/Zee_lightJet_70-140.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/Jets/Zee_lightJet_140-280.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/Jets/Zee_cJet_0-70.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/Jets/Zee_cJet_70-140.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/Jets/Zee_cJet_140-280.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/Jets/Zee_bJet_0-70.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/Jets/Zee_bJet_70-140.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/Jets/Zee_bJet_140-280.root"},
    };

    std::array<double,9> JetNumeratorSFs = {((139e15)*(1.9828e-9)*(0.821204)),((139e15)*(110.64e-12)*(0.69275)),((139e15)*(40.645e-12)*(0.615906)),((139e15)*(1.9817e-9)*(0.1136684)),((139e15)*(110.47e-12)*(0.1912956)),((139e15)*(40.674e-12)*(0.2326772)),((139e15)*(1.9819e-9)*(0.0656969)),((139e15)*(110.53e-12)*(0.1158741)),((139e15)*(40.68e-12)*(0.1535215))}; //numerators for jet bkg

    std::vector<const char*> prefixes = {"pty2_9_17", "pty_17_myy_0_80", "pty_17_myy_80", "sig m_{A} = 5 GeV", "sig m_{A} = 1 GeV", "sig m_{A} = 2 GeV", "sig m_{A} = 3 GeV", "sig m_{A} = 9 GeV", "Zee_lightJet_0-70", "Zee_lightJet_70-140", "Zee_lightJet_140-280", "Zee_cJet_0-70", "Zee_cJet_70-140", "Zee_cJet_140-280", "Zee_bJet_0-70", "Zee_bJet_70-140", "Zee_bJet_140-280"};
    std::array<double,3> SFs = {((139e15)*(.871e-12)),((139e15)*(.199e-12)), ((139e15)*(.0345e-15))}; //numerators for Z-gamma bkg
    std::vector<EColor> colors = {static_cast<EColor>(kRed-7), static_cast<EColor>(kRed-9), static_cast<EColor>(kRed-10), kBlack, kMagenta, kTeal, static_cast<EColor>(kOrange+10), static_cast<EColor>(kSpring + 5)};
    std::vector<EColor> Jetscolors = {static_cast<EColor>(kGreen+0), static_cast<EColor>(kGreen-4), static_cast<EColor>(kGreen-7), static_cast<EColor>(kGreen-9), static_cast<EColor>(kGreen-10), static_cast<EColor>(kGreen+1), static_cast<EColor>(kGreen-3), static_cast<EColor>(kGreen-6), static_cast<EColor>(kGreen-8)};

    TCanvas* c1 = new TCanvas();
    TLegend* legend = new TLegend(0.625, 0.25, 0.88, 0.65);
    int count = 0;
    std::vector<ROOT::RDF::RResultHandle> Nodes;
    std::vector<ROOT::RDF::RResultHandle> VaryNodes;

    std::set<std::pair<std::string,std::string>> columns;

    for (auto& file: input_filenames)
    {
        SchottDataFrame df(MakeRDF(file, 8));
//        df.Describe().Print();
//        exit(1);
//        for (auto& i: df.GetColumnNames())
//        {
//            columns.insert(std::make_pair(i, df.GetColumnType(i)));
//        }

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
                return (!((ep.electron_pt/1e3 > 20) && (std::abs(ep.electron_eta) < 2.37) &&
                          (!((1.37 < std::abs(ep.electron_eta)) && (std::abs(ep.electron_eta) < 1.52)))
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
                return (!((ep.electron_pt/1e3 > 20) && (std::abs(ep.electron_eta) < 2.37) &&
                (!((1.37 < std::abs(ep.electron_eta)) && (std::abs(ep.electron_eta) < 1.52)))
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
            return pT >= 10;
        }, {"di_electrons"});

        auto dilep_mass = ptCut.Define("dilep_mass",[&](RVec<Electron>& electrons)
        {
            PtEtaPhiEVector four_momentum =
            electrons[0].Vector() + electrons[1].Vector();

            return four_momentum.M()/1e3;

        }, {"electrons"});

        if (count <= 2) //Z gamma
        {
            Nodes.push_back(dilep_mass.Histo1D<double>({prefixes[count], prefixes[count++], 60u, 60, 120}, "dilep_mass"));
            Nodes.push_back(dilep_mass.Count()); //Events that pass selection
            Nodes.push_back(df.Count()); //Total number of events in the sample
        }

        else if (count >= 3 && count <= 7) //signal
        {
            Nodes.push_back(dilep_mass.Histo1D<double>({prefixes[count], prefixes[count++], 60u, 60, 120}, "dilep_mass"));
        }

        else //Z+jets
        {
            Nodes.push_back(dilep_mass.Histo1D<double>({prefixes[count], prefixes[count++], 60u, 60, 120}, "dilep_mass"));
            Nodes.push_back(dilep_mass.Count()); //Events that pass selection
            Nodes.push_back(df.Count()); //Total number of events in the sample
        }
    }

//    for (auto& i: columns)
//    {
//        std::cout << std::setw(60) << i.second << std::setw(60) << i.first;
//        std::cout << '\n';
//    }

    ROOT::RDF::RunGraphs(Nodes); // running all computation nodes concurrently

//     Now I should be able to access my results as many times as desired without
//     penalty.

//    0  1  2
//    3  4  5
//    6  7  8
//    9
//    10
//    11
//    12
//    13
//    14 15 16
//    17 18 19
//    20 21 22
//    23 24 25
//    26 27 28
//    29 30 31
//    32 33 34
//    35 36 37
//    38 39 40

    int back_count = 0;
    double factor;
    double total_back = 0;

    std::cout << "Z+gamma\n=======\n";
    std::cout << std::setw(25) << "Sample"
    << std::setw(10) << "Events"
    << std::setw(20) << "Scale Factors"
    << std::setw(20) << "Weighted Events"
    << '\n';
    //calculating total Zgamma background by taking the weighted sum
    for (int j = 1, i=0; (j <= 7 && i<=2) ; j += 3, i++)
    {
        total_back += (*Nodes[j].GetResultPtr<ULong64_t>())*(SFs[back_count++] / *Nodes[j+1].GetResultPtr<ULong64_t>());

        std::cout << std::setw(25) << prefixes[i]
        << std::setw(10) << std::setprecision(2) << std::fixed << (*Nodes[j].GetResultPtr<ULong64_t>())
        << std::setw(20) << std::setprecision(2) << std::fixed << (SFs[back_count-1] / *Nodes[j+1].GetResultPtr<ULong64_t>())
        << std::setw(20) << std::setprecision(2) << std::fixed << (*Nodes[j].GetResultPtr<ULong64_t>())*(SFs[back_count-1] / *Nodes[j+1].GetResultPtr<ULong64_t>())
        << '\n';
    }
    std::cout << "\n\nWeighted Z-gamma sum = " << std::setprecision(2) << std::fixed
    << total_back << "\n\n\n";

    std::cout << "Z+jets\n=======\n";
    std::cout << std::setw(25) << "Sample"
    << std::setw(10) << "Events"
    << std::setw(20) << "Scale Factors"
    << std::setw(20) << "Weighted Events"
    << '\n';

    double jetCount = total_back;
    //then add Zjets to total background
    //Background counts from Z jets
    for (int i = 15, j = 0, k = 5; (i <= 39 && j <= 8); i += 3, j++, k++)
    {
        total_back += *Nodes[i].GetResultPtr<ULong64_t>()*(JetNumeratorSFs[j] / *Nodes[i+1].GetResultPtr<ULong64_t>());

        std::cout << std::setw(25) << prefixes[k]
        << std::setw(10) << std::setprecision(2) << std::fixed << *Nodes[i].GetResultPtr<ULong64_t>()
        << std::setw(20) << std::setprecision(2) << std::fixed << (JetNumeratorSFs[j] / *Nodes[i+1].GetResultPtr<ULong64_t>())
        << std::setw(20) << std::setprecision(2) << std::fixed << *Nodes[i].GetResultPtr<ULong64_t>()*(JetNumeratorSFs[j] / *Nodes[i+1].GetResultPtr<ULong64_t>())
        << '\n';

    }

    jetCount = total_back - jetCount;

    std::cout << "\n\nWeighted Z-Jets sum = " << std::setprecision(2) << std::fixed
    << jetCount << "\n\n\n";

    std::cout << "\n\nWeighted Total Bkg sum = " << std::setprecision(2) << std::fixed
    << total_back << "\n\n\n";


//     First, we calculate the total Zgamma background by taking the weighted sum, i.e.,
//     adding all of the counts of the individual Zgamma samples weighted by their
//     respective scaling factors.
//
//     so total_back will look something like this:
//
//     total_back = bc1*s1 + bc2*s2 + bc3*s3
//
//     Below, we loop over the background Zgamma samples. For each sample, we
//     now just have to scale it by its scale factor, and that way the total THStack
//     will add to total_back as we had intended.

    count = 0;
    factor = total_back;
    auto hs = new THStack("hs3","");

    //Z gamma background
    for (int i = 0; i <= 6; i += 3)
    {
        Nodes[i].GetResultPtr<TH1D>()->SetFillColor(colors[count]);
        Nodes[i].GetResultPtr<TH1D>()->SetLineColor(colors[count]);
        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "f");
        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0)
        {
            Nodes[i].GetResultPtr<TH1D>()->Scale(SFs[count]/(*Nodes[i+2].GetResultPtr<ULong64_t>()));
        }
        Nodes[i].GetResultPtr<TH1D>()->SetTitle(";m_{ll}  [GeV];Events");
        Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);

        hs->Add(&*Nodes[i].GetResultPtr<TH1D>());
        count++;
    }

    //Now adding on top Z+jets background
    //weighted ZJets
    for (int i = 14, j = 0; (i <= 38 && j <= 8); i += 3, j++)
    {
        Nodes[i].GetResultPtr<TH1D>()->SetFillColor(Jetscolors[j]);
        Nodes[i].GetResultPtr<TH1D>()->SetLineColor(Jetscolors[j]);
        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "f");
        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0)
        {
            Nodes[i].GetResultPtr<TH1D>()->Scale(JetNumeratorSFs[j] / (*Nodes[i+2].GetResultPtr<ULong64_t>()));
        }
        Nodes[i].GetResultPtr<TH1D>()->SetTitle(";m_{ll}  [GeV];Events");
        Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);

        hs->Add(&*Nodes[i].GetResultPtr<TH1D>());
    }

    hs->Draw("HIST");

    //Now processing signal sample results
    for (int i = 9; i <= 13; i++)
    {
        Nodes[i].GetResultPtr<TH1D>()->SetLineColor(colors[count++]);
        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "l");

        Nodes[i].GetResultPtr<TH1D>()->Scale(factor/Nodes[i].GetResultPtr<TH1D>()->Integral());
        Nodes[i].GetResultPtr<TH1D>()->SetTitle(";m_{ll}  [GeV];Events");
        Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);
        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0)
        {
//            std::cout << Nodes[i].GetResultPtr<TH1D>()->Integral()
//            << '\n';
            Nodes[i].GetResultPtr<TH1D>()->Scale((factor/Nodes[i].GetResultPtr<TH1D>()->Integral()));


        }
        Nodes[i].GetResultPtr<TH1D>()->Draw("HISTsame");
        gPad->Modified(); gPad->Update();
    }

    hs->SetTitle(";m_{ll}  [GeV];Events");
    hs->GetYaxis()->CenterTitle(true);
    hs->GetYaxis()->SetTitleOffset(1.3);
    hs->GetXaxis()->SetTitleOffset(1.2);
    hs->SetMinimum(0);
    hs->SetMaximum(6.1e6);

    gStyle->SetOptStat(0);
    TLatex Tl;
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.6, 0.83, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.6, 0.73,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("Fig1A.pdf");
}

void fig3()
{
    std::vector<std::vector<std::string>> input_filenames =
    {
        //Z gamma background
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617070._000001.LGNTuple.root"}, //pty2_9_17
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617064._000001.LGNTuple.root"}, //pty_17_myy_0_80
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617074._000001.LGNTuple.root"}, //pty_17_myy_80
        //Signal
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root"}, //5 GeV
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root"}, // 1 GeV
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600751.PhPy8EG_AZNLO_ggH125_mA2p0_v1.root"}, // 2 GeV
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600752.PhPy8EG_AZNLO_ggH125_mA3p0_v1.root"}, // 3 GeV
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600756.PhPy8EG_AZNLO_ggH125_mA9p0_v1.root"}, // 9 GeV
//        //Jets
        {"/Users/edwardfinkelstein/ATLAS_axion/Jets/Zee_lightJet_0-70.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/Jets/Zee_lightJet_70-140.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/Jets/Zee_lightJet_140-280.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/Jets/Zee_cJet_0-70.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/Jets/Zee_cJet_70-140.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/Jets/Zee_cJet_140-280.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/Jets/Zee_bJet_0-70.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/Jets/Zee_bJet_70-140.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/Jets/Zee_bJet_140-280.root"},
    };

    std::array<double,9> JetNumeratorSFs = {((139e15)*(1.9828e-9)*(0.821204)),((139e15)*(110.64e-12)*(0.69275)),((139e15)*(40.645e-12)*(0.615906)),((139e15)*(1.9817e-9)*(0.1136684)),((139e15)*(110.47e-12)*(0.1912956)),((139e15)*(40.674e-12)*(0.2326772)),((139e15)*(1.9819e-9)*(0.0656969)),((139e15)*(110.53e-12)*(0.1158741)),((139e15)*(40.68e-12)*(0.1535215))}; //numerators for jet bkg

    std::vector<const char*> prefixes = {"pty2_9_17", "pty_17_myy_0_80", "pty_17_myy_80", "sig m_{A} = 5 GeV", "sig m_{A} = 1 GeV", "sig m_{A} = 2 GeV", "sig m_{A} = 3 GeV", "sig m_{A} = 9 GeV", "Zee_lightJet_0-70", "Zee_lightJet_70-140", "Zee_lightJet_140-280", "Zee_cJet_0-70", "Zee_cJet_70-140", "Zee_cJet_140-280", "Zee_bJet_0-70", "Zee_bJet_70-140", "Zee_bJet_140-280"};
    std::array<double,3> SFs = {((139e15)*(.871e-12)),((139e15)*(.199e-12)), ((139e15)*(.0345e-15))}; //numerators for Z-gamma bkg
    std::vector<EColor> colors = {static_cast<EColor>(kRed-7), static_cast<EColor>(kRed-9), static_cast<EColor>(kRed-10), kBlack, kMagenta, kTeal, static_cast<EColor>(kOrange+10), static_cast<EColor>(kSpring + 5)};
    std::vector<EColor> Jetscolors = {static_cast<EColor>(kGreen+0), static_cast<EColor>(kGreen-4), static_cast<EColor>(kGreen-7), static_cast<EColor>(kGreen-9), static_cast<EColor>(kGreen-10), static_cast<EColor>(kGreen+1), static_cast<EColor>(kGreen-3), static_cast<EColor>(kGreen-6), static_cast<EColor>(kGreen-8)};

    TCanvas* c1 = new TCanvas();
    TLegend* legend = new TLegend(0.625, 0.25, 0.88, 0.65);
    int count = 0;
    std::vector<ROOT::RDF::RResultHandle> Nodes;
    std::vector<ROOT::RDF::RResultHandle> VaryNodes;

    std::set<std::pair<std::string,std::string>> columns;

    for (auto& file: input_filenames)
    {
        SchottDataFrame df(MakeRDF(file, 8));

        auto dilep_pt = df.Filter(
        [&](const RVec<Electron>& electrons)
        {
            return electrons.size() >= 2;
        }, {"electrons"})
        .Filter(
        [&](const RVec<float>& electron_pt)
        {
//            R__ASSERT(std::is_sorted(electron_pt.begin(), electron_pt.end(), [](const float x, const float y){return x >= y;}));

            return true;
        }, {"electron_pt"})
        .Define("dilep_pt",[&]( RVec<Electron> electrons)
        {
            std::sort(electrons.begin(), electrons.end(), //sorting electrons
            [](const Electron& x, const Electron& y)
            {
                return x.electron_pt >= y.electron_pt;
            });

            PtEtaPhiEVector four_momentum =
            electrons[0].Vector() + electrons[1].Vector();

            return four_momentum.Pt()/1e3;

        }, {"electrons"});

        if (count <= 2) //Z gamma
        {
            Nodes.push_back(dilep_pt.Histo1D<double>({prefixes[count], prefixes[count++], 100u, 0, 200}, "dilep_pt"));
            Nodes.push_back(dilep_pt.Count());
            Nodes.push_back(df.Count()); //Total number of events in the sample
        }

        else if (count >= 3 && count <= 7) //signal
        {
            Nodes.push_back(dilep_pt.Histo1D<double>({prefixes[count], prefixes[count++], 100u, 0, 200}, "dilep_pt"));
        }

        else //Z+jets
        {
            Nodes.push_back(dilep_pt.Histo1D<double>({prefixes[count], prefixes[count++], 100u, 0, 200}, "dilep_pt"));
            Nodes.push_back(dilep_pt.Count());
            Nodes.push_back(df.Count()); //Total number of events in the sample
        }
    }

    ROOT::RDF::RunGraphs(Nodes); // running all computation nodes concurrently

//     Now I should be able to access my results as many times as desired without
//     penalty.

//    0  1  2
//    3  4  5
//    6  7  8
//    9
//    10
//    11
//    12
//    13
//    14 15 16
//    17 18 19
//    20 21 22
//    23 24 25
//    26 27 28
//    29 30 31
//    32 33 34
//    35 36 37
//    38 39 40

    int back_count = 0;
    double factor;
    double total_back = 0;

    std::cout << "Z+gamma\n=======\n";
    std::cout << std::setw(25) << "Sample"
    << std::setw(10) << "Events"
    << std::setw(20) << "Scale Factors"
    << std::setw(20) << "Weighted Events"
    << '\n';
    //calculating total Zgamma background by taking the weighted sum
    for (int j = 1, i=0; (j <= 7 && i<=2) ; j += 3, i++)
    {
        total_back += (*Nodes[j].GetResultPtr<ULong64_t>())*(SFs[back_count++] / *Nodes[j+1].GetResultPtr<ULong64_t>());

        std::cout << std::setw(25) << prefixes[i]
        << std::setw(10) << std::setprecision(2) << std::fixed << (*Nodes[j].GetResultPtr<ULong64_t>())
        << std::setw(20) << std::setprecision(2) << std::fixed << (SFs[back_count-1] / *Nodes[j+1].GetResultPtr<ULong64_t>())
        << std::setw(20) << std::setprecision(2) << std::fixed << (*Nodes[j].GetResultPtr<ULong64_t>())*(SFs[back_count-1] / *Nodes[j+1].GetResultPtr<ULong64_t>())
        << '\n';
    }
    std::cout << "\n\nWeighted Z-gamma sum = " << std::setprecision(2) << std::fixed
    << total_back << "\n\n\n";

    std::cout << "Z+jets\n=======\n";
    std::cout << std::setw(25) << "Sample"
    << std::setw(10) << "Events"
    << std::setw(20) << "Scale Factors"
    << std::setw(20) << "Weighted Events"
    << '\n';

    double jetCount = total_back;
    //then add Zjets to total background
    //Background counts from Z jets
    for (int i = 15, j = 0, k = 5; (i <= 39 && j <= 8); i += 3, j++, k++)
    {
        total_back += *Nodes[i].GetResultPtr<ULong64_t>()*(JetNumeratorSFs[j] / *Nodes[i+1].GetResultPtr<ULong64_t>());

        std::cout << std::setw(25) << prefixes[k]
        << std::setw(10) << std::setprecision(2) << std::fixed << *Nodes[i].GetResultPtr<ULong64_t>()
        << std::setw(20) << std::setprecision(2) << std::fixed << (JetNumeratorSFs[j] / *Nodes[i+1].GetResultPtr<ULong64_t>())
        << std::setw(20) << std::setprecision(2) << std::fixed << *Nodes[i].GetResultPtr<ULong64_t>()*(JetNumeratorSFs[j] / *Nodes[i+1].GetResultPtr<ULong64_t>())
        << '\n';

    }

    jetCount = total_back - jetCount;

    std::cout << "\n\nWeighted Z-Jets sum = " << std::setprecision(2) << std::fixed
    << jetCount << "\n\n\n";

    std::cout << "\n\nWeighted Total Bkg sum = " << std::setprecision(2) << std::fixed
    << total_back << "\n\n\n";


//     First, we calculate the total Zgamma background by taking the weighted sum, i.e.,
//     adding all of the counts of the individual Zgamma samples weighted by their
//     respective scaling factors (s1, s2, s3).
//
//     so total_back will look something like this:
//
//     total_back = bc1*s1 + bc2*s2 + bc3*s3
//
//     Below, we loop over the background Zgamma samples. For each sample, we
//     now just have to scale it by its scale factor, and that way the total THStack
//     will add to total_back as we had intended.

    count = 0;
    factor = total_back;
    auto hs = new THStack("hs3","");

    //Z gamma background
    for (int i = 0; i <= 6; i += 3)
    {
        Nodes[i].GetResultPtr<TH1D>()->SetFillColor(colors[count]);
        Nodes[i].GetResultPtr<TH1D>()->SetLineColor(colors[count]);
        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "f");
        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0)
        {
            Nodes[i].GetResultPtr<TH1D>()->Scale(SFs[count]/(*Nodes[i+2].GetResultPtr<ULong64_t>()));
        }
        Nodes[i].GetResultPtr<TH1D>()->SetTitle(";p_{T_{ee}}  [GeV];Events");
        Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);

        hs->Add(&*Nodes[i].GetResultPtr<TH1D>());
        count++;
    }

    //Now adding on top Z+jets background
    //weighted ZJets
    for (int i = 14, j = 0; (i <= 38 && j <= 8); i += 3, j++)
    {
        Nodes[i].GetResultPtr<TH1D>()->SetFillColor(Jetscolors[j]);
        Nodes[i].GetResultPtr<TH1D>()->SetLineColor(Jetscolors[j]);
        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "f");
        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0)
        {
            Nodes[i].GetResultPtr<TH1D>()->Scale(JetNumeratorSFs[j] / (*Nodes[i+2].GetResultPtr<ULong64_t>()));
        }
        Nodes[i].GetResultPtr<TH1D>()->SetTitle(";p_{T_{ee}}  [GeV];Events");
        Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);

        hs->Add(&*Nodes[i].GetResultPtr<TH1D>());
    }

    hs->Draw("HIST");

    //Now processing signal sample results
    for (int i = 9; i <= 13; i++)
    {
        Nodes[i].GetResultPtr<TH1D>()->SetLineColor(colors[count++]);
        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "l");

        Nodes[i].GetResultPtr<TH1D>()->Scale(factor/Nodes[i].GetResultPtr<TH1D>()->Integral());
        Nodes[i].GetResultPtr<TH1D>()->SetTitle(";p_{T_{ee}}  [GeV];Events");
        Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);
        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0)
        {
//            std::cout << Nodes[i].GetResultPtr<TH1D>()->Integral()
//            << '\n';
            Nodes[i].GetResultPtr<TH1D>()->Scale((factor/Nodes[i].GetResultPtr<TH1D>()->Integral()));


        }
        Nodes[i].GetResultPtr<TH1D>()->Draw("HISTsame");
        gPad->Modified(); gPad->Update();
    }

    hs->SetTitle(";p_{T_{ee}}  [GeV];Events");
    hs->GetYaxis()->CenterTitle(true);
    hs->GetYaxis()->SetTitleOffset(1.3);
    hs->GetXaxis()->SetTitleOffset(1.2);
    hs->SetMinimum(0);
    hs->SetMaximum(8e6);

    gStyle->SetOptStat(0);
    TLatex Tl;
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.6, 0.83, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.6, 0.73,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("Fig3.pdf");
}


//void fig5()
//{
//    //No photon selection, only lepton preselection!
//    std::vector<std::string> input_filenames = {
//        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root",
////        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Test_oneMassPoint.root",
//    };
//    std::vector<EColor> colors = {kRed, kBlue, kBlack};
//    TCanvas* c1 = new TCanvas();
//    TLegend* legend = new TLegend(0.6, 0.4, 0.8, 0.525);
//
//    SchottDataFrame df(MakeRDF(input_filenames, 8));
//
//    auto preselection = df.Filter(
//    [&](const RVec<std::string>& trigger_passed_triggers, RVec<TruthParticle> truth_particles)
//    {
//        //Selecting the leptons and anti-leptons in the event
//        truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
//        [](TruthParticle& x)
//        {
//            //so, remove elements that are not leptons and/or don't satisfy pt > 20 GeV, |η| < 2.37, and 1.37 ≮ |η| ≮ 1.52 and/or don't have mc_status=1
//            return (std::abs(x.mc_pdg_id) != 11 && std::abs(x.mc_pdg_id) != 12 && std::abs(x.mc_pdg_id) != 13 &&
//                    std::abs(x.mc_pdg_id) != 14 && std::abs(x.mc_pdg_id) != 15 && std::abs(x.mc_pdg_id) != 16 &&
//                    std::abs(x.mc_pdg_id) != 17 && std::abs(x.mc_pdg_id) != 18)
//            ||
//            //Lepton preselection
//            (!((x.mc_pt/1e3 > 20) && (std::abs(x.mc_eta) < 2.37) &&
//                   (!((1.37 < std::abs(x.mc_eta)) && (std::abs(x.mc_eta) < 1.52)))))
//            ||
//            (x.mc_status != 1);
//
//        }), truth_particles.end());
//
//        //must be exactly two leptons/anti-leptons
//        if (truth_particles.size() != 2)
//        {
//            return false;
//        }
//
//        //must be same flavor
//        if (std::abs(truth_particles[0].mc_pdg_id) != std::abs(truth_particles[1].mc_pdg_id))
//        {
//            return false;
//        }
//
//        if (DeltaR(truth_particles[0].Vector(), truth_particles[1].Vector()) <= 0.01)
//        {
//            return false;
//        }
//
//        if (truth_particles[0].mc_charge*truth_particles[1].mc_charge >= 0)
//        {
//            return false;
//        }
//
//        if (!((truth_particles[0].mc_pt > 27e3 && truth_particles[1].mc_pt > 20e3)
//                                               ||
//              (truth_particles[1].mc_pt > 27e3 && truth_particles[0].mc_pt > 20e3)))
//        {
//            return false;
//        }
//
//        PtEtaPhiEVector dilepton = truth_particles[0].Vector() + truth_particles[1].Vector();
//
//        if ((dilepton.M() < 81e3) || (dilepton.M() > 101e3))
//        {
//            return false;
//        }
//
//        if ((truth_particles[0].Vector() + truth_particles[1].Vector()).Pt() <= 10e3)
//        {
//            return false;
//        }
//
//        return true;
//
//    }, {"trigger_passed_triggers", "truth_particles"});
//
//    auto leading_truth_photons = preselection.Define("leading_truth_photons",[&](RVec<TruthParticle> truth_particles)
//    {
//        truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
//        [](TruthParticle& x)
//        {
//            return (std::abs(x.mc_pdg_id) != 22);
//
//        }), truth_particles.end());
//
//        std::sort(truth_particles.begin(),truth_particles.end(),
//        [](TruthParticle& x, TruthParticle& y)
//        {
//            return x.mc_pt > y.mc_pt;
//        });
//
//        RVec<TruthParticle> chosen = {truth_particles[0], truth_particles[1]};
//
//        return chosen;
//
//    }, {"truth_particles"})
//    .Define("leading_truth_photon_pt", [&](RVec<TruthParticle> leading_truth_photons)
//    {
//        return leading_truth_photons[0].mc_pt/1e3;
//    }, {"leading_truth_photons"})
//    .Define("subleading_truth_photon_pt", [&](RVec<TruthParticle> leading_truth_photons)
//    {
//        return leading_truth_photons[1].mc_pt/1e3;
//    }, {"leading_truth_photons"});
//
//    auto leading_truth_photons_pt_10 = leading_truth_photons.Filter(
//    [&](RVec<TruthParticle>& leading_truth_photons)
//    {
//        return (leading_truth_photons[0].mc_pt > 10e3 && leading_truth_photons[1].mc_pt > 10e3);
//    }, {"leading_truth_photons"})
//    .Define("leading_truth_photons_pt_10_dR",
//    [&](RVec<TruthParticle>& leading_truth_photons)
//    {
//        return DeltaR(leading_truth_photons[0].Vector(), leading_truth_photons[1].Vector());
//    }, {"leading_truth_photons"});
//
//    auto leading_truth_photons_pt_5 = leading_truth_photons.Filter(
//    [&](RVec<TruthParticle>& leading_truth_photons)
//    {
//        return (leading_truth_photons[0].mc_pt > 5e3 && leading_truth_photons[1].mc_pt > 5e3);
//    }, {"leading_truth_photons"})
//    .Define("leading_truth_photons_pt_5_dR",
//    [&](RVec<TruthParticle>& leading_truth_photons)
//    {
//        return DeltaR(leading_truth_photons[0].Vector(), leading_truth_photons[1].Vector());
//    }, {"leading_truth_photons"});
//
//    auto leading_truth_photons_pt_0 = leading_truth_photons.Filter(
//    [&](RVec<TruthParticle>& leading_truth_photons)
//    {
//        return (leading_truth_photons[0].mc_pt > 0 && leading_truth_photons[1].mc_pt > 0);
//    }, {"leading_truth_photons"})
//    .Define("leading_truth_photons_pt_0_dR",
//    [&](RVec<TruthParticle>& leading_truth_photons)
//    {
//        return DeltaR(leading_truth_photons[0].Vector(), leading_truth_photons[1].Vector());
//    }, {"leading_truth_photons"});
//
//    std::vector<ROOT::RDF::RResultPtr<TH1D>> histos =
//    {
//        leading_truth_photons.Histo1D<double>({"Leading Photon", "Leading Photon", 20u, 0, 50}, "leading_truth_photon_pt"),
//        leading_truth_photons.Histo1D<double>({"Sub-Leading Photon", "Sub-Leading Photon", 20u, 0, 50}, "subleading_truth_photon_pt"),
//        leading_truth_photons_pt_0.Histo1D<double>({"Photon p_{T} > 0 GeV", "Photon p_{T} > 0 GeV", 40u, 0, 0.32}, "leading_truth_photons_pt_0_dR"),
//        leading_truth_photons_pt_5.Histo1D<double>({"Photon p_{T} > 5 GeV", "Photon p_{T} > 5 GeV", 40u, 0, 0.32}, "leading_truth_photons_pt_5_dR"),
//        leading_truth_photons_pt_10.Histo1D<double>({"Photon p_{T} > 10 GeV", "Photon p_{T} > 10 GeV", 40u, 0, 0.32}, "leading_truth_photons_pt_10_dR"),
//    };
//
//    double factor = histos[0]->Integral();
//    histos[0]->Scale(factor/factor);
//    histos[0]->SetLineColor(colors[0]);
//    histos[0]->Draw("HIST");
//    legend->AddEntry(&(*histos[0]), histos[0]->GetTitle(), "l");
//    histos[0]->SetTitle(";Truth photon p_{T}  [GeV];Events");
//    histos[0]->SetTitleOffset(1.2);
//    histos[0]->GetYaxis()->CenterTitle(true);
//    histos[0]->SetAxisRange(0., 12000,"Y");
//
//    histos[1]->Scale(factor/histos[1]->Integral());
//    histos[1]->SetLineColor(colors[1]);
//    histos[1]->Draw("HISTsame");
//    legend->AddEntry(&(*histos[1]), histos[1]->GetTitle(), "l");
//
//    gStyle->SetOptStat(0);
//    TLatex Tl;
//    Tl.SetTextSize(0.03);
//    Tl.DrawLatexNDC(0.6, 0.83, "#it{ATLAS} Internal");
//    Tl.DrawLatexNDC(0.6, 0.77,"#sqrt{s} = 13 TeV");
//    Tl.DrawLatexNDC(0.6, 0.7,"ggF m_{A} = 1 GeV");
//    legend->SetBorderSize(0);
//    legend->Draw();
//    c1->SaveAs("Fig5A.pdf");
//    c1 = new TCanvas();
//    legend = new TLegend(0.6, 0.4, 0.8, 0.6);
//    factor = histos[2]->Integral();
//
//    histos[2]->Scale(factor/factor);
//    histos[2]->SetLineColor(colors[0]);
//    histos[2]->Draw("HIST");
//
//    histos[2]->GetYaxis()->CenterTitle(true);
//    histos[2]->SetAxisRange(0., 7000,"Y");
//
//    legend->AddEntry(&(*histos[2]), histos[2]->GetTitle(), "l");
//    histos[2]->SetTitle(";Truth #DeltaR_{#gamma#gamma}; Events");
//    histos[3]->Scale(factor/histos[3]->Integral());
//    histos[3]->SetLineColor(colors[1]);
//    histos[3]->Draw("HISTsame");
//    legend->AddEntry(&(*histos[3]), histos[3]->GetTitle(), "l");
//
//    histos[4]->Scale(factor/histos[4]->Integral());
//    histos[4]->SetLineColor(colors[2]);
//    histos[4]->Draw("HISTsame");
//    legend->AddEntry(&(*histos[4]), histos[4]->GetTitle(), "l");
//
//    gStyle->SetOptStat(0);
//    Tl.SetTextSize(0.03);
//    Tl.DrawLatexNDC(0.6, 0.83, "#it{ATLAS} Internal");
//    Tl.DrawLatexNDC(0.6, 0.77,"#sqrt{s} = 13 TeV");
//    Tl.DrawLatexNDC(0.6, 0.7,"ggF m_{A} = 1 GeV");
//    legend->SetBorderSize(0);
//    legend->Draw();
//    c1->SaveAs("Fig5B.pdf");
//}


//void fig6()
//{
//    std::vector<std::string> input_filenames = {
//        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root",
////        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Test_oneMassPoint.root"
//    };
//
//    std::vector<EColor> colors = {kRed, kBlue, static_cast<EColor>(kOrange+2)};
//    TCanvas* c1 = new TCanvas();
//    TLegend* legend = new TLegend(0.6, 0.4, 0.8, 0.525);
//
//    SchottDataFrame df(MakeRDF(input_filenames, 8));
//
//    auto findParentInChain = [](int targetBarcode, RVec<TruthParticle>& startParticles, RVec<TruthParticle>& truthChain)
//    {
//        RVec<TruthParticle> truthSelected;
//        bool foundParent;
//        if (truthChain.size() >= 1)
//        {
//            TruthParticle tp;
//            for (auto& tpe: startParticles)
//            {
//                tp = tpe;
//                while (true)
//                {
//                    if (tp.mc_parent_barcode == targetBarcode)
//                    {
//                        truthSelected.push_back(tp);
//                        break;
//                    }
//                    else
//                    {
//                        foundParent = false;
//                        for (auto& tmp: truthChain)
//                        {
//                            if (tp.mc_parent_barcode == tmp.mc_barcode)
//                            {
//                                tp = tmp;
//                                foundParent = true;
//                                break;
//                            }
//                        }
//                        if (foundParent == false)
//                        {
//                            break;
//                        }
//                    }
//                }
//            }
//        }
//        return truthSelected;
//    };
//
//    auto preselection = df.Filter(
//    [&](const RVec<std::string>& trigger_passed_triggers, RVec<TruthParticle> truth_particles)
//    {
//
//        truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
//        [](TruthParticle& x)
//        {
//            //so, remove elements that are not leptons and/or don't satisfy pt > 20 GeV, |η| < 2.37, and 1.37 ≮ |η| ≮ 1.52 and/or don't have mc_status=1
//            return (std::abs(x.mc_pdg_id) != 11 && std::abs(x.mc_pdg_id) != 12 && std::abs(x.mc_pdg_id) != 13 &&
//                    std::abs(x.mc_pdg_id) != 14 && std::abs(x.mc_pdg_id) != 15 && std::abs(x.mc_pdg_id) != 16 &&
//                    std::abs(x.mc_pdg_id) != 17 && std::abs(x.mc_pdg_id) != 18)
//            ||
//            (!((x.mc_pt/1e3 > 20) && (std::abs(x.mc_eta) < 2.37) &&
//                   (!((1.37 < std::abs(x.mc_eta)) && (std::abs(x.mc_eta) < 1.52)))))
//            ||
//            (x.mc_status != 1);
//
//        }), truth_particles.end());
//
//        //must be exactly two leptons/anti-leptons
//        if (truth_particles.size() != 2)
//        {
//            return false;
//        }
//
//        //must be same flavor
//        if (std::abs(truth_particles[0].mc_pdg_id) != std::abs(truth_particles[1].mc_pdg_id))
//        {
//            return false;
//        }
//
//        if (DeltaR(truth_particles[0].Vector(), truth_particles[1].Vector()) <= 0.01)
//        {
//            return false;
//        }
//
//        if (truth_particles[0].mc_charge*truth_particles[1].mc_charge >= 0)
//        {
//            return false;
//        }
//
//        if (!((truth_particles[0].mc_pt > 27e3 && truth_particles[1].mc_pt > 20e3)
//                                               ||
//              (truth_particles[1].mc_pt > 27e3 && truth_particles[0].mc_pt > 20e3)))
//        {
//            return false;
//        }
//
//        PtEtaPhiEVector dilepton = truth_particles[0].Vector() + truth_particles[1].Vector();
//
//        if ((dilepton.M() < 81e3) || (dilepton.M() > 101e3))
//        {
//            return false;
//        }
//
//        if ((truth_particles[0].Vector() + truth_particles[1].Vector()).Pt() <= 10e3)
//        {
//            return false;
//        }
//
//        return true;
//
//    }, {"trigger_passed_triggers", "truth_particles"});
//
////    std::cout << *preselection.Count() << '\n';
//
//    auto truth_photons_from_axions = preselection.Define("truth_photons",[&](RVec<TruthParticle> truth_particles)
//    {
//        truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
//        [](TruthParticle& x)
//        {
//            return (std::abs(x.mc_pdg_id) != 22);
//
//        }), truth_particles.end());
//
//        return truth_particles;
//
//    }, {"truth_particles"})
//    .Define("truth_axions", [&](RVec<TruthParticle> truth_particles)
//    {
//        truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
//        [&](TruthParticle& x)
//        {
//            return (std::abs(x.mc_pdg_id) != 35 && std::abs(x.mc_pdg_id) != 36);
//
//        }), truth_particles.end());
//
//        return truth_particles;
//
//    }, {"truth_particles"}).Define("truth_photons_from_axions",
//    [&](RVec<TruthParticle>& truth_photons, RVec<TruthParticle>& truth_axions)
//    {
//        return findParentInChain(truth_axions[0].mc_barcode, truth_photons, truth_axions);
//
//    }, {"truth_photons", "truth_axions"})
//    .Filter([&] (RVec<TruthParticle>& truth_photons_from_axions)
//    {
//        return truth_photons_from_axions.size()==2;
//    }, {"truth_photons_from_axions"})
//    .Define("photons_pass_cut_indices",
//    [&](const RVec<Photon>& photons)
//    {
//        RVec<int> photon_indices;
//        photon_indices.reserve(photons.size());
//
//        for (int i = 0; i < photons.size(); i++)
//        {
//            if ((std::abs(photons[i].photon_eta) >= 2.37) or (photons[i].photon_pt <= 10e3) or (std::abs(photons[i].photon_eta) > 1.37 and std::abs(photons[i].photon_eta) < 1.52) or (not photons[i].photon_id_loose))
//            {
//                continue;
//            }
//            photon_indices.push_back(i);
//        }
//
//        return photon_indices;
//    }, {"photons"})
//    .Define("photons_pass_cuts",
//    [&](RVec<Photon>& photons, RVec<int>& photon_indices)
//    {
//        return Take(photons, photon_indices);
//    }, {"photons", "photons_pass_cut_indices"});
//
//    auto reco_photons_matched = truth_photons_from_axions
//    .Define("reco_photons_matched_indices",
//    [&](RVec<Photon>& photons_pass_cuts, RVec<TruthParticle>& truth_photons_from_axions)
//    {
//        RVec<int> reco_photons_matched_indices;
//        PtEtaPhiEVector tp1 = truth_photons_from_axions[0].Vector();
//        PtEtaPhiEVector tp2 = truth_photons_from_axions[1].Vector();
//
//        for (int i = 0; i < photons_pass_cuts.size(); i++)
//        {
//            if (not (DeltaR(photons_pass_cuts[i].Vector(), tp1) < 0.1 or DeltaR(photons_pass_cuts[i].Vector(), tp2) < 0.1))
//            {
//                continue;
//            }
//            reco_photons_matched_indices.push_back(i);
//        }
//
//
//
//        return reco_photons_matched_indices;
//
//    }, {"photons_pass_cuts", "truth_photons_from_axions"})
//    .Define("reco_photons_matched",
//    [&](RVec<Photon>& photons_pass_cuts, RVec<int>& reco_photons_matched_indices)
//    {
//        return Take(photons_pass_cuts, reco_photons_matched_indices);
//
//    }, {"photons_pass_cuts", "reco_photons_matched_indices"})
//    .Define("totEventWeight",
//    [](RVec<float> photon_id_eff, RVec<float> photon_iso_eff, RVec<float> photon_trg_eff/*, RVec<float> ei_event_weights_generator*/, RVec<int>& photons_pass_cut_indices, RVec<int>& reco_photons_matched_indices)
//    {
//        //Taking the values that correspond to the
//        auto ResizeVal = std::max({photon_id_eff.size(), photon_iso_eff.size(), photon_trg_eff.size()});
//        photon_id_eff.resize(ResizeVal,1);
//        photon_iso_eff.resize(ResizeVal,1);
//        photon_trg_eff.resize(ResizeVal,1);
////        R__ASSERT(ResizeVal);
//        return
//        Take(Take(photon_id_eff, photons_pass_cut_indices), reco_photons_matched_indices);
////        *
////        Take(Take(photon_iso_eff, photons_pass_cut_indices), reco_photons_matched_indices)*
////        Take(Take(photon_trg_eff, photons_pass_cut_indices), reco_photons_matched_indices);//*ei_event_weights_generator[0];
//
//    }, {"photon_id_eff", "photon_iso_eff", "photon_trg_eff"/*, "ei_event_weights_generator"*/, "photons_pass_cut_indices", "reco_photons_matched_indices"});
//
//    auto two_reco_photons_matched = reco_photons_matched.Filter(
//    [&](RVec<Photon>& reco_photons_matched)
//    {
//        return reco_photons_matched.size()==2;
//
//    }, {"reco_photons_matched"})
//    .Define("leading_pt",[&] (RVec<Photon>& reco_photons_matched)
//    {
//        if (reco_photons_matched[0].photon_pt/1e3 > reco_photons_matched[1].photon_pt/1e3)
//        {
//            return reco_photons_matched[0].photon_pt/1e3;
//        }
//        return reco_photons_matched[1].photon_pt/1e3;
//    }, {"reco_photons_matched"})
//    .Define("sub_leading_pt",[&] (RVec<Photon>& reco_photons_matched)
//    {
//        if (reco_photons_matched[0].photon_pt/1e3 < reco_photons_matched[1].photon_pt/1e3)
//        {
//            return reco_photons_matched[0].photon_pt/1e3;
//        }
//        return reco_photons_matched[1].photon_pt/1e3;
//    }, {"reco_photons_matched"});
//
//    auto one_reco_photons_matched = reco_photons_matched.Filter(
//    [&](RVec<Photon>& reco_photons_matched)
//    {
//        return reco_photons_matched.size()==1;
//
//    }, {"reco_photons_matched"})
//    .Define("leading_pt",[&] (RVec<Photon>& reco_photons_matched)
//    {
//        return reco_photons_matched[0].photon_pt/1e3;
//    }, {"reco_photons_matched"});
//
//
//    std::vector<ROOT::RDF::RResultPtr<TH1D>> histos =
//    {
//        two_reco_photons_matched.Histo1D<double>({"Leading Photon Two Matched", "Leading Photon Two Matched", 30u, 0, 50}, "leading_pt", "totEventWeight"),
//        two_reco_photons_matched.Histo1D<double>({"Sub-Leading Photon Two Matched", "Sub-Leading Photon Two Matched", 30u, 0, 50}, "sub_leading_pt", "totEventWeight"),
//        one_reco_photons_matched.Histo1D<double>({"Leading Photon One Matched", "Leading Photon One Matched", 40u, 0, 100}, "leading_pt", "totEventWeight"),
//    };
//
////    auto truth_photons_from_axions_Count = truth_photons_from_axions.Count();
////    std::cout << *truth_photons_from_axions_Count << '\n';
////    truth_photons_from_axions_Count = two_reco_photons_matched.Count();
////    std::cout << *truth_photons_from_axions_Count << '\n';
////    truth_photons_from_axions_Count = one_reco_photons_matched.Count();
////    std::cout << *truth_photons_from_axions_Count << '\n';
//
//    double factor = histos[0]->Integral();
//
////    histos[0]->Scale(factor/factor);
//    histos[0]->SetLineColor(colors[0]);
//    histos[0]->Draw("HIST");
//    legend->AddEntry(&(*histos[0]), histos[0]->GetTitle(), "l");
//    histos[0]->SetTitle(";photon p_{T}  [GeV];Events");
//    histos[0]->SetTitleOffset(1.2);
//    histos[0]->GetYaxis()->CenterTitle(true);
//    histos[0]->SetAxisRange(0., 800,"Y");
//
////    histos[1]->Scale(factor/histos[1]->Integral());
//    histos[1]->SetLineColor(colors[1]);
//    histos[1]->Draw("HISTsame");
//    legend->AddEntry(&(*histos[1]), histos[1]->GetTitle(), "l");
//
//    gStyle->SetOptStat(0);
//    TLatex Tl;
//    Tl.SetTextSize(0.03);
//    Tl.DrawLatexNDC(0.6, 0.83, "#it{ATLAS} Internal");
//    Tl.DrawLatexNDC(0.6, 0.77,"#sqrt{s} = 13 TeV");
//    Tl.DrawLatexNDC(0.6, 0.7,"ggF m_{A} = 1 GeV");
//    legend->SetBorderSize(0);
//    legend->Draw();
//    legend->SetTextSize(0.022);
//    c1->SaveAs("Fig6A.pdf");
//
//    c1 = new TCanvas();
//    legend = new TLegend(0.6, 0.4, 0.8, 0.6);
//    factor = histos[2]->Integral();
//    histos[2]->Scale(factor/factor);
//    histos[2]->SetLineColor(colors[2]);
//    histos[2]->Draw("HIST");
//    histos[2]->GetYaxis()->CenterTitle(true);
//
//    legend->AddEntry(&(*histos[2]), histos[2]->GetTitle(), "l");
//    histos[2]->SetTitle(";photon p_{T}  [GeV];Events");
//
//    gStyle->SetOptStat(0);
//    Tl.SetTextSize(0.03);
//    Tl.DrawLatexNDC(0.6, 0.83, "#it{ATLAS} Internal");
//    Tl.DrawLatexNDC(0.6, 0.77,"#sqrt{s} = 13 TeV");
//    Tl.DrawLatexNDC(0.6, 0.7,"ggF m_{A} = 1 GeV");
//    legend->SetBorderSize(0);
//    legend->SetTextSize(0.022);
//    legend->Draw();
//    c1->SaveAs("Fig6B.pdf");
//
//    c1 = new TCanvas();
//    legend = new TLegend(0.6, 0.4, 0.8, 0.6);
//
//    factor = histos[2]->Integral();
//
//    histos[0]->Draw("HIST");
//    histos[0]->Scale(factor/histos[0]->Integral());
//    legend->AddEntry(&(*histos[0]), "Leading Photon Two Matched", "l");
//    histos[0]->SetTitle(";photon p_{T}  [GeV];Events");
//    histos[0]->SetTitleOffset(1.2);
//    histos[0]->GetYaxis()->CenterTitle(true);
//    histos[0]->SetAxisRange(0., 8400,"Y");
//
//    histos[1]->SetLineColor(colors[1]);
//    histos[1]->Draw("HISTsame");
//    histos[1]->Scale(factor/histos[1]->Integral());
//    legend->AddEntry(&(*histos[1]), "Sub-Leading Photon Two Matched", "l");
//
//    histos[2]->SetLineColor(colors[2]);
//    histos[2]->Draw("HISTsame");
//    legend->AddEntry(&(*histos[2]), "Leading Photon One Matched", "l");
//
//    gStyle->SetOptStat(0);
//    Tl.SetTextSize(0.03);
//    Tl.DrawLatexNDC(0.6, 0.83, "#it{ATLAS} Internal");
//    Tl.DrawLatexNDC(0.6, 0.77,"#sqrt{s} = 13 TeV");
//    Tl.DrawLatexNDC(0.6, 0.7,"ggF m_{A} = 1 GeV");
//    legend->SetBorderSize(0);
//    legend->Draw();
//    legend->SetTextSize(0.022);
//    gPad->Modified(); gPad->Update();
//
//    c1->SaveAs("Fig6_A_and_B.pdf");
//}


//void fig8()
//{
//    std::vector<std::string> input_filenames = {
//        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600756.PhPy8EG_AZNLO_ggH125_mA9p0_v1.root",
//        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root"
//    };
//
//    auto findParentInChain = [](int targetBarcode, RVec<TruthParticle>& startParticles, RVec<TruthParticle>& truthChain)
//    {
//        RVec<TruthParticle> truthSelected;
//        bool foundParent;
//        if (truthChain.size() >= 1)
//        {
//            TruthParticle tp;
//            for (auto& tpe: startParticles)
//            {
//                tp = tpe;
//                while (true)
//                {
//                    if (tp.mc_parent_barcode == targetBarcode)
//                    {
//                        truthSelected.push_back(tp);
//                        break;
//                    }
//                    else
//                    {
//                        foundParent = false;
//                        for (auto& tmp: truthChain)
//                        {
//                            if (tp.mc_parent_barcode == tmp.mc_barcode)
//                            {
//                                tp = tmp;
//                                foundParent = true;
//                                break;
//                            }
//                        }
//                        if (foundParent == false)
//                        {
//                            break;
//                        }
//                    }
//                }
//            }
//        }
//        return truthSelected;
//    };
//
//    std::vector<ROOT::RDF::RResultHandle> Nodes;
//
//    std::vector<const char*> prefixes = {"sig ggF m_{A} = 9 GeV", "sig ggF m_{A} = 1 GeV"};
//    std::vector<EColor> colors = {kMagenta, kBlack};
//    int count = 0;
//
//    for (auto& file: input_filenames)
//    {
//        SchottDataFrame df(MakeRDF({file}, 8));
//
//        auto preselection = df.Filter(
//        [&](const RVec<std::string>& trigger_passed_triggers, RVec<TruthParticle> truth_particles)
//        {
//
//            truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
//            [](TruthParticle& x)
//            {
//                //so, remove elements that are not leptons and/or don't satisfy pt > 20 GeV, |η| < 2.37, and 1.37 ≮ |η| ≮ 1.52 and/or don't have mc_status=1
//                return (std::abs(x.mc_pdg_id) != 11 && std::abs(x.mc_pdg_id) != 12 && std::abs(x.mc_pdg_id) != 13 &&
//                        std::abs(x.mc_pdg_id) != 14 && std::abs(x.mc_pdg_id) != 15 && std::abs(x.mc_pdg_id) != 16 &&
//                        std::abs(x.mc_pdg_id) != 17 && std::abs(x.mc_pdg_id) != 18)
//                ||
//                (!((x.mc_pt/1e3 > 20) && (std::abs(x.mc_eta) < 2.37) &&
//                       (!((1.37 < std::abs(x.mc_eta)) && (std::abs(x.mc_eta) < 1.52)))))
//                ||
//                (x.mc_status != 1);
//
//            }), truth_particles.end());
//
//            //must be exactly two leptons/anti-leptons
//            if (truth_particles.size() != 2)
//            {
//                return false;
//            }
//
//            //must be same flavor
//            if (std::abs(truth_particles[0].mc_pdg_id) != std::abs(truth_particles[1].mc_pdg_id))
//            {
//                return false;
//            }
//
//            if (DeltaR(truth_particles[0].Vector(), truth_particles[1].Vector()) <= 0.01)
//            {
//                return false;
//            }
//
//            if (truth_particles[0].mc_charge*truth_particles[1].mc_charge >= 0)
//            {
//                return false;
//            }
//
//            if (!((truth_particles[0].mc_pt > 27e3 && truth_particles[1].mc_pt > 20e3)
//                                                   ||
//                  (truth_particles[1].mc_pt > 27e3 && truth_particles[0].mc_pt > 20e3)))
//            {
//                return false;
//            }
//
//            PtEtaPhiEVector dilepton = truth_particles[0].Vector() + truth_particles[1].Vector();
//
//            if ((dilepton.M() < 81e3) || (dilepton.M() > 101e3))
//            {
//                return false;
//            }
//
//            if ((truth_particles[0].Vector() + truth_particles[1].Vector()).Pt() <= 10e3)
//            {
//                return false;
//            }
//
//            return true;
//
//        }, {"trigger_passed_triggers", "truth_particles"});
//
//        auto truth_photons_from_axions = preselection.Define("truth_photons",[&](RVec<TruthParticle> truth_particles)
//        {
//            truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
//            [](TruthParticle& x)
//            {
//                return (std::abs(x.mc_pdg_id) != 22);
//
//            }), truth_particles.end());
//
//            return truth_particles;
//
//        }, {"truth_particles"})
//        .Define("truth_axions", [&](RVec<TruthParticle> truth_particles)
//        {
//            truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
//            [](TruthParticle& x)
//            {
//                return (std::abs(x.mc_pdg_id) != 35);
//
//            }), truth_particles.end());
//
//            return truth_particles;
//
//        }, {"truth_particles"}).Define("truth_photons_from_axions",
//        [&](RVec<TruthParticle>& truth_photons, RVec<TruthParticle>& truth_axions)
//        {
//            return findParentInChain(truth_axions[0].mc_barcode, truth_photons, truth_axions);
//
//        }, {"truth_photons", "truth_axions"})
//        .Filter([&] (RVec<TruthParticle>& truth_photons_from_axions)
//        {
//            return truth_photons_from_axions.size()==2;
//        }, {"truth_photons_from_axions"})
//        .Define("photons_pass_cut_indices",
//        [&](const RVec<Photon>& photons)
//        {
//            RVec<int> photon_indices;
//            photon_indices.reserve(photons.size());
//
//            for (int i = 0; i < photons.size(); i++)
//            {
//                if ((std::abs(photons[i].photon_eta) >= 2.37) or (photons[i].photon_pt <= 10e3) or (std::abs(photons[i].photon_eta) > 1.37 and std::abs(photons[i].photon_eta) < 1.52) or (not photons[i].photon_id_loose))
//                {
//                    continue;
//                }
//                photon_indices.push_back(i);
//            }
//
//            return photon_indices;
//        }, {"photons"})
//        .Define("photons_pass_cuts",
//        [&](RVec<Photon>& photons, RVec<int>& photon_indices)
//        {
//            return Take(photons, photon_indices);
//        }, {"photons", "photons_pass_cut_indices"});
//
//        auto reco_photons_matched = truth_photons_from_axions
//        .Define("reco_photons_matched_indices",
//        [&](RVec<Photon>& photons_pass_cuts, RVec<TruthParticle>& truth_photons_from_axions)
//        {
//            RVec<int> reco_photons_matched_indices;
//            PtEtaPhiEVector tp1 = truth_photons_from_axions[0].Vector();
//            PtEtaPhiEVector tp2 = truth_photons_from_axions[1].Vector();
//
//            for (int i = 0; i < photons_pass_cuts.size(); i++)
//            {
//                if (not (DeltaR(photons_pass_cuts[i].Vector(), tp1) < 0.1 or DeltaR(photons_pass_cuts[i].Vector(), tp2) < 0.1))
//                {
//                    continue;
//                }
//                reco_photons_matched_indices.push_back(i);
//            }
//
//            return reco_photons_matched_indices;
//
//        }, {"photons_pass_cuts", "truth_photons_from_axions"})
//        .Define("reco_photons_matched",
//        [&](RVec<Photon>& photons_pass_cuts, RVec<int>& reco_photons_matched_indices)
//        {
//            return Take(photons_pass_cuts, reco_photons_matched_indices);
//
//        }, {"photons_pass_cuts", "reco_photons_matched_indices"})
//        .Define("totEventWeight",
//        [](RVec<float> photon_id_eff, RVec<float> photon_iso_eff, RVec<float> photon_trg_eff/*, RVec<float> ei_event_weights_generator*/, RVec<int>& photons_pass_cut_indices, RVec<int>& reco_photons_matched_indices)
//        {
//            //Taking the values that correspond to the
//            auto ResizeVal = std::max({photon_id_eff.size(), photon_iso_eff.size(), photon_trg_eff.size()});
//            photon_id_eff.resize(ResizeVal,1);
//            photon_iso_eff.resize(ResizeVal,1);
//            photon_trg_eff.resize(ResizeVal,1);
//            return
//            Take(Take(photon_id_eff, photons_pass_cut_indices), reco_photons_matched_indices);
////            *
////            Take(Take(photon_iso_eff, photons_pass_cut_indices), reco_photons_matched_indices)*
////            Take(Take(photon_trg_eff, photons_pass_cut_indices), reco_photons_matched_indices);//*ei_event_weights_generator[0];
//
//        }, {"photon_id_eff", "photon_iso_eff", "photon_trg_eff"/*, "ei_event_weights_generator"*/, "photons_pass_cut_indices", "reco_photons_matched_indices"});
//
//        auto two_reco_photons_matched = reco_photons_matched.Filter(
//        [&](RVec<Photon>& reco_photons_matched)
//        {
//            return reco_photons_matched.size()==2;
//
//        }, {"reco_photons_matched"});
//
//        auto stable_truth_dileptons_and_diphotons = two_reco_photons_matched
//        .Define("stable_truth_leptons",[&](RVec<TruthParticle> truth_particles)
//        {
//            truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
//            [](TruthParticle& x)
//            {
//                //so, remove elements that are not leptons and/or don't satisfy pt > 20 GeV, |η| < 2.37, and 1.37 ≮ |η| ≮ 1.52 and/or don't have mc_status=1
//                return (std::abs(x.mc_pdg_id) != 11 && std::abs(x.mc_pdg_id) != 12 && std::abs(x.mc_pdg_id) != 13 &&
//                        std::abs(x.mc_pdg_id) != 14 && std::abs(x.mc_pdg_id) != 15 && std::abs(x.mc_pdg_id) != 16 &&
//                        std::abs(x.mc_pdg_id) != 17 && std::abs(x.mc_pdg_id) != 18)
//                ||
//                (!((x.mc_pt/1e3 > 20) && (std::abs(x.mc_eta) < 2.37) &&
//                       (!((1.37 < std::abs(x.mc_eta)) && (std::abs(x.mc_eta) < 1.52)))))
//                ||
//                (x.mc_status != 1);
//
//            }), truth_particles.end());
//
//            return truth_particles;
//
//        }, {"truth_particles"})
//        .Filter([&](RVec<TruthParticle>& stable_truth_leptons, RVec<TruthParticle>& truth_photons_from_axions)
//        {
//            return (stable_truth_leptons.size()==2 && truth_photons_from_axions.size()==2);
//
//        }, {"stable_truth_leptons", "truth_photons_from_axions"})
//        .Define("reconstructed_Z",[&](RVec<TruthParticle>& stable_truth_leptons)
//        {
//            auto four_momentum = stable_truth_leptons[0].Vector() + stable_truth_leptons[1].Vector();
//
//            return four_momentum.Pt()/1e3;
//
//        }, {"stable_truth_leptons"})
//        .Define("reconstructed_a",[&](RVec<TruthParticle>& truth_photons_from_axions)
//        {
//            PtEtaPhiEVector four_momentum = truth_photons_from_axions[0].Vector() + truth_photons_from_axions[1].Vector();
//            return four_momentum.Pt()/1e3;
//
//        }, {"truth_photons_from_axions"})
//        .Define("reconstructed_h",[&](RVec<TruthParticle>& truth_photons_from_axions, RVec<TruthParticle>& stable_truth_leptons)
//        {
//            PtEtaPhiEVector four_momentum_photons = truth_photons_from_axions[0].Vector() + truth_photons_from_axions[1].Vector();
//            PtEtaPhiEVector four_momentum_leptons = stable_truth_leptons[0].Vector() + stable_truth_leptons[1].Vector();
//            return (four_momentum_photons + four_momentum_leptons).Pt()/1e3;
//        }, {"truth_photons_from_axions", "stable_truth_leptons"})
//        .Define("delta_r_Z_a",[&](RVec<TruthParticle>& truth_photons_from_axions, RVec<TruthParticle>& stable_truth_leptons)
//        {
//            PtEtaPhiEVector four_momentum_photons = truth_photons_from_axions[0].Vector() + truth_photons_from_axions[1].Vector();
//            PtEtaPhiEVector four_momentum_leptons = stable_truth_leptons[0].Vector() + stable_truth_leptons[1].Vector();
//            return DeltaR(four_momentum_photons, four_momentum_leptons);
//        }, {"truth_photons_from_axions", "stable_truth_leptons"})
//        .Define("delta_phi_Z_a",[&](RVec<TruthParticle>& truth_photons_from_axions, RVec<TruthParticle>& stable_truth_leptons)
//        {
//            PtEtaPhiEVector four_momentum_photons = truth_photons_from_axions[0].Vector() + truth_photons_from_axions[1].Vector();
//            PtEtaPhiEVector four_momentum_leptons = stable_truth_leptons[0].Vector() + stable_truth_leptons[1].Vector();
//
//            return ROOT::Math::VectorUtil::DeltaPhi(four_momentum_photons, four_momentum_leptons);
//        }, {"truth_photons_from_axions", "stable_truth_leptons"})
//        .Define("delta_eta_Z_a",[&](RVec<TruthParticle>& truth_photons_from_axions, RVec<TruthParticle>& stable_truth_leptons)
//        {
//            PtEtaPhiEVector four_momentum_photons = truth_photons_from_axions[0].Vector() + truth_photons_from_axions[1].Vector();
//            PtEtaPhiEVector four_momentum_leptons = stable_truth_leptons[0].Vector() + stable_truth_leptons[1].Vector();
//
//            return std::abs((four_momentum_photons - four_momentum_leptons).Eta());
//        }, {"truth_photons_from_axions", "stable_truth_leptons"});
//
//        Nodes.push_back(stable_truth_dileptons_and_diphotons.Histo1D<double>({prefixes[count], prefixes[count], 50u, 0, 200}, "reconstructed_Z", "totEventWeight"));
//        Nodes.push_back(stable_truth_dileptons_and_diphotons.Histo1D<double>({prefixes[count], prefixes[count], 50u, 0, 62}, "reconstructed_a", "totEventWeight"));
//        Nodes.push_back(stable_truth_dileptons_and_diphotons.Histo1D<double>({prefixes[count], prefixes[count], 100u, 0, 217}, "reconstructed_h", "totEventWeight"));
//        Nodes.push_back(stable_truth_dileptons_and_diphotons.Histo1D<double>({prefixes[count], prefixes[count], 100u, 0, 6.4}, "delta_r_Z_a", "totEventWeight"));
//        Nodes.push_back(stable_truth_dileptons_and_diphotons.Histo1D<double>({prefixes[count], prefixes[count], 100u, 0, 6.25}, "delta_phi_Z_a", "totEventWeight"));
//        Nodes.push_back(stable_truth_dileptons_and_diphotons.Histo1D<double>({prefixes[count], prefixes[count++], 100u, 0, 6.25}, "delta_eta_Z_a", "totEventWeight"));
//    }
//
//    ROOT::RDF::RunGraphs(Nodes); // running all computation nodes concurrently
//
//    TCanvas* c1 = new TCanvas();
//    TLegend* legend = new TLegend(0.5, 0.3, 0.8, 0.6);
//    double factor;
//
//    count = 0;
//    for (auto& i: {0,6})
//    {
//        Nodes[i].GetResultPtr<TH1D>()->SetLineColor(colors[count++]);
//        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "l");
//
//        if (i==0)
//        {
//            factor = Nodes[i].GetResultPtr<TH1D>()->Integral();
//            Nodes[i].GetResultPtr<TH1D>()->Scale(factor/Nodes[i].GetResultPtr<TH1D>()->Integral());
//            Nodes[i].GetResultPtr<TH1D>()->SetTitle(";ll p_{T} [GeV];Events");
//            Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->SetTitleSize(1.25 * Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->GetTitleSize());
//            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->SetTitleSize(1.25 * Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->GetTitleSize());
//            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);
////            Nodes[i].GetResultPtr<TH1D>()->SetAxisRange(0., 100,"Y");
//            Nodes[i].GetResultPtr<TH1D>()->Draw("HIST");
//        }
//        else
//        {
//            Nodes[i].GetResultPtr<TH1D>()->Scale(factor/Nodes[i].GetResultPtr<TH1D>()->Integral());
//            Nodes[i].GetResultPtr<TH1D>()->Draw("HISTsame");
//            gPad->Modified(); gPad->Update();
//        }
//    }
//
//    gStyle->SetOptStat(0);
//    TLatex Tl;
//    Tl.SetTextSize(0.03);
//    Tl.DrawLatexNDC(0.6, 0.83, "#it{ATLAS} Internal");
//    Tl.DrawLatexNDC(0.6, 0.73,"#sqrt{s} = 13 TeV");
//    legend->SetBorderSize(0);
//    legend->Draw();
//    c1->SaveAs("Fig8A.pdf");
//
//    c1 = new TCanvas();
//    legend = new TLegend(0.55, 0.5, 0.85, 0.8);
//    count = 0;
//    for (auto& i: {1, 7})
//    {
//        Nodes[i].GetResultPtr<TH1D>()->SetLineColor(colors[count++]);
//        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "l");
//
//        if (i==1)
//        {
//            factor = Nodes[i].GetResultPtr<TH1D>()->Integral();
//            Nodes[i].GetResultPtr<TH1D>()->Scale(factor/Nodes[i].GetResultPtr<TH1D>()->Integral());
//            Nodes[i].GetResultPtr<TH1D>()->SetTitle(";#gamma#gamma p_{T} [GeV];Events");
//            Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->SetTitleSize(1.25 * Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->GetTitleSize());
//            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->SetTitleSize(1.25 * Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->GetTitleSize());
//            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);
//            Nodes[i].GetResultPtr<TH1D>()->SetAxisRange(0., 450,"Y");
//            Nodes[i].GetResultPtr<TH1D>()->Draw("HIST");
//        }
//        else
//        {
//            Nodes[i].GetResultPtr<TH1D>()->Scale(factor/Nodes[i].GetResultPtr<TH1D>()->Integral());
//            Nodes[i].GetResultPtr<TH1D>()->SetTitle(";#gamma#gamma p_{T} [GeV];Events");
//            Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->SetTitleSize(1.25 * Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->GetTitleSize());
//            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->SetTitleSize(1.25 * Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->GetTitleSize());
//            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);
//            Nodes[i].GetResultPtr<TH1D>()->SetAxisRange(0., 450,"Y");
//            Nodes[i].GetResultPtr<TH1D>()->Draw("HISTsame");
//            gPad->Modified(); gPad->Update();
//        }
//    }
//    gStyle->SetOptStat(0);
//    Tl.SetTextSize(0.03);
//    Tl.DrawLatexNDC(0.13, 0.83, "#it{ATLAS} Internal");
//    Tl.DrawLatexNDC(0.13, 0.73,"#sqrt{s} = 13 TeV");
//    legend->SetBorderSize(0);
//    legend->Draw();
//    c1->SaveAs("Fig8B.pdf");
//
//    c1 = new TCanvas();
//    legend = new TLegend(0.55, 0.35, 0.85, 0.65);
//    count = 0;
//    for (auto& i: {2, 8})
//    {
//        Nodes[i].GetResultPtr<TH1D>()->SetLineColor(colors[count++]);
//        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "l");
//
//        if (i == 2)
//        {
//            factor = Nodes[i].GetResultPtr<TH1D>()->Integral();
//            Nodes[i].GetResultPtr<TH1D>()->Scale(factor/Nodes[i].GetResultPtr<TH1D>()->Integral());
//            Nodes[i].GetResultPtr<TH1D>()->SetTitle(";ll#gamma#gamma p_{T} [GeV];Events");
//            Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->SetTitleSize(1.25 * Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->GetTitleSize());
//            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->SetTitleSize(1.25 * Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->GetTitleSize());
//            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);
//            Nodes[i].GetResultPtr<TH1D>()->SetAxisRange(0., 240,"Y");
//            Nodes[i].GetResultPtr<TH1D>()->Draw("HIST");
//        }
//        else
//        {
//            Nodes[i].GetResultPtr<TH1D>()->Scale(factor/Nodes[i].GetResultPtr<TH1D>()->Integral());
//            Nodes[i].GetResultPtr<TH1D>()->SetTitle(";ll#gamma#gamma p_{T} [GeV];Events");
//            Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->SetTitleSize(1.25 * Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->GetTitleSize());
//            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->SetTitleSize(1.25 * Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->GetTitleSize());
//            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);
//            Nodes[i].GetResultPtr<TH1D>()->SetAxisRange(0., 240,"Y");
//            Nodes[i].GetResultPtr<TH1D>()->Draw("HISTsame");
//            gPad->Modified(); gPad->Update();
//        }
//    }
//    gStyle->SetOptStat(0);
//    Tl.SetTextSize(0.03);
//    Tl.DrawLatexNDC(0.6, 0.83, "#it{ATLAS} Internal");
//    Tl.DrawLatexNDC(0.6, 0.73,"#sqrt{s} = 13 TeV");
//    legend->SetBorderSize(0);
//    legend->Draw();
//    c1->SaveAs("Fig8C.pdf");
//
//    c1 = new TCanvas();
//    legend = new TLegend(0.55, 0.35, 0.85, 0.65);
//    count = 0;
//    for (auto& i: {3, 9})
//    {
//        Nodes[i].GetResultPtr<TH1D>()->SetLineColor(colors[count++]);
//        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "l");
//
//        if (i == 3)
//        {
//            factor = Nodes[i].GetResultPtr<TH1D>()->Integral();
//            Nodes[i].GetResultPtr<TH1D>()->Scale(factor/Nodes[i].GetResultPtr<TH1D>()->Integral());
//            Nodes[i].GetResultPtr<TH1D>()->SetTitle(";#DeltaR ll#gamma#gamma;Events");
//            Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->SetTitleSize(1.25 * Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->GetTitleSize());
//            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->SetTitleSize(1.25 * Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->GetTitleSize());
//            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);
//            Nodes[i].GetResultPtr<TH1D>()->SetAxisRange(0., 199,"Y");
//            Nodes[i].GetResultPtr<TH1D>()->Draw("HIST");
//        }
//        else
//        {
//            Nodes[i].GetResultPtr<TH1D>()->Scale(factor/Nodes[i].GetResultPtr<TH1D>()->Integral());
//            Nodes[i].GetResultPtr<TH1D>()->SetTitle(";#DeltaR ll#gamma#gamma;Events");
//            Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->SetTitleSize(1.25 * Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->GetTitleSize());
//            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->SetTitleSize(1.25 * Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->GetTitleSize());
//            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);
//            Nodes[i].GetResultPtr<TH1D>()->SetAxisRange(0., 199,"Y");
//            Nodes[i].GetResultPtr<TH1D>()->Draw("HISTsame");
//            gPad->Modified(); gPad->Update();
//        }
//    }
//    gStyle->SetOptStat(0);
//    Tl.SetTextSize(0.03);
//    Tl.DrawLatexNDC(0.6, 0.83, "#it{ATLAS} Internal");
//    Tl.DrawLatexNDC(0.6, 0.73,"#sqrt{s} = 13 TeV");
//    legend->SetBorderSize(0);
//    legend->Draw();
//    c1->SaveAs("Fig8D.pdf");
//
//    c1 = new TCanvas();
//    legend = new TLegend(0.55, 0.35, 0.85, 0.65);
//    count = 0;
//    for (auto& i: {4,10})
//    {
//        Nodes[i].GetResultPtr<TH1D>()->SetLineColor(colors[count++]);
//        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "l");
//
//        if (i == 4)
//        {
//            factor = Nodes[i].GetResultPtr<TH1D>()->Integral();
//            Nodes[i].GetResultPtr<TH1D>()->Scale(factor/Nodes[i].GetResultPtr<TH1D>()->Integral());
//            Nodes[i].GetResultPtr<TH1D>()->SetTitle(";#Delta#phi ll#gamma#gamma;Events");
//            Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->SetTitleSize(1.25 * Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->GetTitleSize());
//            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->SetTitleSize(1.25 * Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->GetTitleSize());
//            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);
//            Nodes[i].GetResultPtr<TH1D>()->SetAxisRange(0., 163,"Y");
//            Nodes[i].GetResultPtr<TH1D>()->Draw("HIST");
//        }
//        else
//        {
//            Nodes[i].GetResultPtr<TH1D>()->Scale(factor/Nodes[i].GetResultPtr<TH1D>()->Integral());
//            Nodes[i].GetResultPtr<TH1D>()->SetTitle(";#Delta#phi ll#gamma#gamma;Events");
//            Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->SetTitleSize(1.25 * Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->GetTitleSize());
//            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->SetTitleSize(1.25 * Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->GetTitleSize());
//            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);
//            Nodes[i].GetResultPtr<TH1D>()->SetAxisRange(0., 163,"Y");
//            Nodes[i].GetResultPtr<TH1D>()->Draw("HISTsame");
//            gPad->Modified(); gPad->Update();
//        }
//    }
//    gStyle->SetOptStat(0);
//    Tl.SetTextSize(0.03);
//    Tl.DrawLatexNDC(0.6, 0.83, "#it{ATLAS} Internal");
//    Tl.DrawLatexNDC(0.6, 0.73,"#sqrt{s} = 13 TeV");
//    legend->SetBorderSize(0);
//    legend->Draw();
//    c1->SaveAs("Fig8E.pdf");
//
//    c1 = new TCanvas();
//    legend = new TLegend(0.55, 0.35, 0.85, 0.65);
//    count = 0;
//    for (auto& i: {5, 11})
//    {
//        Nodes[i].GetResultPtr<TH1D>()->SetLineColor(colors[count++]);
//        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "l");
//
//        if (i == 5)
//        {
//            factor = Nodes[i].GetResultPtr<TH1D>()->Integral();
//            Nodes[i].GetResultPtr<TH1D>()->Scale(factor/Nodes[i].GetResultPtr<TH1D>()->Integral());
//            Nodes[i].GetResultPtr<TH1D>()->SetTitle(";#Delta#eta ll#gamma#gamma;Events");
//            Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->SetTitleSize(1.25 * Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->GetTitleSize());
//            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->SetTitleSize(1.25 * Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->GetTitleSize());
//            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);
////            Nodes[i].GetResultPtr<TH1D>()->SetAxisRange(0., 45,"Y");
//            Nodes[i].GetResultPtr<TH1D>()->Draw("HIST");
//        }
//        else
//        {
//            Nodes[i].GetResultPtr<TH1D>()->Scale(factor/Nodes[i].GetResultPtr<TH1D>()->Integral());
//            Nodes[i].GetResultPtr<TH1D>()->SetTitle(";#Delta#eta ll#gamma#gamma;Events");
//            Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->SetTitleSize(1.25 * Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->GetTitleSize());
//            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->SetTitleSize(1.25 * Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->GetTitleSize());
//            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);
////            Nodes[i].GetResultPtr<TH1D>()->SetAxisRange(0., 45,"Y");
//            Nodes[i].GetResultPtr<TH1D>()->Draw("HISTsame");
//            gPad->Modified(); gPad->Update();
//        }
//    }
//    gStyle->SetOptStat(0);
//    Tl.SetTextSize(0.03);
//    Tl.DrawLatexNDC(0.6, 0.83, "#it{ATLAS} Internal");
//    Tl.DrawLatexNDC(0.6, 0.73,"#sqrt{s} = 13 TeV");
//    legend->SetBorderSize(0);
//    legend->Draw();
//    c1->SaveAs("Fig8F.pdf");
//}
//
//void fig10()
//{
//    std::vector<std::string> input_filenames =
//    {
//        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root", //1 GeV
//        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600751.PhPy8EG_AZNLO_ggH125_mA2p0_v1.root", // 2 GeV
//        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600752.PhPy8EG_AZNLO_ggH125_mA3p0_v1.root", // 3 GeV
//        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root", // 5 GeV
//        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600756.PhPy8EG_AZNLO_ggH125_mA9p0_v1.root", // 9 GeV
//    };
//
//    auto findParentInChain = [](int targetBarcode, RVec<TruthParticle>& startParticles, RVec<TruthParticle>& truthChain)
//    {
//        RVec<TruthParticle> truthSelected;
//        bool foundParent;
//        if (truthChain.size() >= 1)
//        {
//            TruthParticle tp;
//            for (auto& tpe: startParticles)
//            {
//                tp = tpe;
//                while (true)
//                {
//                    if (tp.mc_parent_barcode == targetBarcode)
//                    {
//                        truthSelected.push_back(tp);
//                        break;
//                    }
//                    else
//                    {
//                        foundParent = false;
//                        for (auto& tmp: truthChain)
//                        {
//                            if (tp.mc_parent_barcode == tmp.mc_barcode)
//                            {
//                                tp = tmp;
//                                foundParent = true;
//                                break;
//                            }
//                        }
//                        if (foundParent == false)
//                        {
//                            break;
//                        }
//                    }
//                }
//            }
//        }
//        return truthSelected;
//    };
//
//    std::vector<ROOT::RDF::RResultHandle> Nodes;
//
//    constexpr std::array<const char*, 5> prefixes = {"sig ggF m_{A} = 1 GeV", "sig ggF m_{A} = 2 GeV", "sig ggF m_{A} = 3 GeV", "sig ggF m_{A} = 5 GeV", "sig ggF m_{A} = 9 GeV"};
//    std::vector<EColor> colors = {kBlue, kBlack, static_cast<EColor>(kGreen + 2), kRed, kViolet, kMagenta};
//    int count = 0;
//
//    for (auto& file: input_filenames)
//    {
//        SchottDataFrame df(MakeRDF({file}, 8));
//
//        auto preselection = df.Filter(
//        [&](const RVec<std::string>& trigger_passed_triggers, RVec<TruthParticle> truth_particles)
//        {
//            //Selecting the leptons and anti-leptons in the event
//            truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
//            [](TruthParticle& x)
//            {
//                //so, remove elements that are not leptons and/or don't satisfy pt > 20 GeV, |η| < 2.37, and 1.37 ≮ |η| ≮ 1.52 and/or don't have mc_status=1
//                return (std::abs(x.mc_pdg_id) != 11 && std::abs(x.mc_pdg_id) != 12 && std::abs(x.mc_pdg_id) != 13 &&
//                        std::abs(x.mc_pdg_id) != 14 && std::abs(x.mc_pdg_id) != 15 && std::abs(x.mc_pdg_id) != 16 &&
//                        std::abs(x.mc_pdg_id) != 17 && std::abs(x.mc_pdg_id) != 18)
//                ||
//                //Lepton preselection
//                (!((x.mc_pt/1e3 > 20) && (std::abs(x.mc_eta) < 2.37) &&
//                       (!((1.37 < std::abs(x.mc_eta)) && (std::abs(x.mc_eta) < 1.52)))))
//                ||
//                (x.mc_status != 1);
//
//            }), truth_particles.end());
//
//            //must be exactly two leptons/anti-leptons
//            if (truth_particles.size() != 2)
//            {
//                return false;
//            }
//
//            //must be same flavor
//            if (std::abs(truth_particles[0].mc_pdg_id) != std::abs(truth_particles[1].mc_pdg_id))
//            {
//                return false;
//            }
//
//            if (DeltaR(truth_particles[0].Vector(), truth_particles[1].Vector()) <= 0.01)
//            {
//                return false;
//            }
//
//            if (truth_particles[0].mc_charge*truth_particles[1].mc_charge >= 0)
//            {
//                return false;
//            }
//
//            if (!((truth_particles[0].mc_pt > 27e3 && truth_particles[1].mc_pt > 20e3)
//                                                   ||
//                  (truth_particles[1].mc_pt > 27e3 && truth_particles[0].mc_pt > 20e3)))
//            {
//                return false;
//            }
//
//            PtEtaPhiEVector dilepton = truth_particles[0].Vector() + truth_particles[1].Vector();
//
//            if ((dilepton.M() < 81e3) || (dilepton.M() > 101e3))
//            {
//                return false;
//            }
//
//            if ((truth_particles[0].Vector() + truth_particles[1].Vector()).Pt() <= 10e3)
//            {
//                return false;
//            }
//
//            return true;
//
//        }, {"trigger_passed_triggers", "truth_particles"});
//
//        auto truth_photons_from_axions = preselection.Define("truth_photons",[&](RVec<TruthParticle> truth_particles)
//        {
//            truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
//            [](TruthParticle& x)
//            {
//                return (std::abs(x.mc_pdg_id) != 22);
//
//            }), truth_particles.end());
//
//            return truth_particles;
//
//        }, {"truth_particles"})
//        .Define("truth_axions", [&](RVec<TruthParticle> truth_particles)
//        {
//            truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
//            [](TruthParticle& x)
//            {
//                return (std::abs(x.mc_pdg_id) != 35);
//
//            }), truth_particles.end());
//
//            return truth_particles;
//
//        }, {"truth_particles"}).Define("truth_photons_from_axions",
//        [&](RVec<TruthParticle>& truth_photons, RVec<TruthParticle>& truth_axions)
//        {
//            return findParentInChain(truth_axions[0].mc_barcode, truth_photons, truth_axions);
//
//        }, {"truth_photons", "truth_axions"})
//        .Filter([&] (RVec<TruthParticle>& truth_photons_from_axions)
//        {
//            return truth_photons_from_axions.size()==2;
//        }, {"truth_photons_from_axions"})
//        .Define("photons_pass_cut_indices",
//        [&](const RVec<Photon>& photons)
//        {
//            RVec<int> photon_indices;
//            photon_indices.reserve(photons.size());
//
//            for (int i = 0; i < photons.size(); i++)
//            {
//                if ((std::abs(photons[i].photon_eta) >= 2.37) or (photons[i].photon_pt <= 10e3) or (std::abs(photons[i].photon_eta) > 1.37 and std::abs(photons[i].photon_eta) < 1.52) or (not photons[i].photon_id_loose))
//                {
//                    continue;
//                }
//                photon_indices.push_back(i);
//            }
//
//            return photon_indices;
//        }, {"photons"})
//        .Define("photons_pass_cuts",
//        [&](RVec<Photon>& photons, RVec<int>& photon_indices)
//        {
//            return Take(photons, photon_indices);
//        }, {"photons", "photons_pass_cut_indices"});
//
//        auto reco_photons_matched = truth_photons_from_axions
//        .Define("reco_photons_matched_indices",
//        [&](RVec<Photon>& photons_pass_cuts, RVec<TruthParticle>& truth_photons_from_axions)
//        {
//            RVec<int> reco_photons_matched_indices;
//            PtEtaPhiEVector tp1 = truth_photons_from_axions[0].Vector();
//            PtEtaPhiEVector tp2 = truth_photons_from_axions[1].Vector();
//
//            for (int i = 0; i < photons_pass_cuts.size(); i++)
//            {
//                if (not (DeltaR(photons_pass_cuts[i].Vector(), tp1) < 0.1 or DeltaR(photons_pass_cuts[i].Vector(), tp2) < 0.1))
//                {
//                    continue;
//                }
//                reco_photons_matched_indices.push_back(i);
//            }
//
//            return reco_photons_matched_indices;
//
//        }, {"photons_pass_cuts", "truth_photons_from_axions"})
//        .Define("reco_photons_matched",
//        [&](RVec<Photon>& photons_pass_cuts, RVec<int>& reco_photons_matched_indices)
//        {
//            return Take(photons_pass_cuts, reco_photons_matched_indices);
//
//        }, {"photons_pass_cuts", "reco_photons_matched_indices"})
//        .Define("totEventWeight",
//        [](RVec<float> photon_id_eff, RVec<float> photon_iso_eff, RVec<float> photon_trg_eff/*, RVec<float> ei_event_weights_generator*/, RVec<int>& photons_pass_cut_indices, RVec<int>& reco_photons_matched_indices)
//        {
//            //Taking the values that correspond to the
//            auto ResizeVal = std::max({photon_id_eff.size(), photon_iso_eff.size(), photon_trg_eff.size()});
//            photon_id_eff.resize(ResizeVal,1);
//            photon_iso_eff.resize(ResizeVal,1);
//            photon_trg_eff.resize(ResizeVal,1);
//    //        R__ASSERT(ResizeVal);
//            return
//            Take(Take(photon_id_eff, photons_pass_cut_indices), reco_photons_matched_indices);
////            *
////            Take(Take(photon_iso_eff, photons_pass_cut_indices), reco_photons_matched_indices)*
////            Take(Take(photon_trg_eff, photons_pass_cut_indices), reco_photons_matched_indices);//*ei_event_weights_generator[0];
//
//        }, {"photon_id_eff", "photon_iso_eff", "photon_trg_eff"/*, "ei_event_weights_generator"*/, "photons_pass_cut_indices", "reco_photons_matched_indices"});
//
//        auto two_reco_photons_matched = reco_photons_matched.Filter(
//        [&](RVec<Photon>& reco_photons_matched)
//        {
//            return reco_photons_matched.size()==2;
//
//        }, {"reco_photons_matched"});
//
//        auto stable_truth_dileptons_and_diphotons = two_reco_photons_matched
//        .Define("stable_truth_leptons",[&](RVec<TruthParticle> truth_particles)
//        {
//            truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
//            [](TruthParticle& x)
//            {
//                //so, remove elements that are not leptons and/or don't satisfy pt > 20 GeV, |η| < 2.37, and 1.37 ≮ |η| ≮ 1.52 and/or don't have mc_status=1
//                return (std::abs(x.mc_pdg_id) != 11 && std::abs(x.mc_pdg_id) != 12 && std::abs(x.mc_pdg_id) != 13 &&
//                        std::abs(x.mc_pdg_id) != 14 && std::abs(x.mc_pdg_id) != 15 && std::abs(x.mc_pdg_id) != 16 &&
//                        std::abs(x.mc_pdg_id) != 17 && std::abs(x.mc_pdg_id) != 18)
//                ||
//                (!((x.mc_pt/1e3 > 20) && (std::abs(x.mc_eta) < 2.37) &&
//                       (!((1.37 < std::abs(x.mc_eta)) && (std::abs(x.mc_eta) < 1.52)))))
//                ||
//                (x.mc_status != 1);
//
//            }), truth_particles.end());
//
//            return truth_particles;
//
//        }, {"truth_particles"})
//        .Filter([&](RVec<TruthParticle>& stable_truth_leptons, RVec<TruthParticle>& truth_photons_from_axions)
//        {
//            return (stable_truth_leptons.size()==2 && truth_photons_from_axions.size()==2);
//
//        }, {"stable_truth_leptons", "truth_photons_from_axions"})
//        .Define("delta_r_gamma_gamma",[&](RVec<TruthParticle>& truth_photons_from_axions)
//        {
//            return DeltaR(truth_photons_from_axions[0].Vector(),truth_photons_from_axions[1].Vector());
//
//        }, {"truth_photons_from_axions"});
//
//        Nodes.push_back(stable_truth_dileptons_and_diphotons.Histo1D<double>({prefixes[count], prefixes[count], 100u, 0, 1.05}, "delta_r_gamma_gamma", "totEventWeight"));
//
//        count++;
//
//    }
//
//    ROOT::RDF::RunGraphs(Nodes); // running all computation nodes concurrently
//
//    TCanvas* c1 = new TCanvas();
//    TLegend* legend = new TLegend(0.6, 0.4, 0.8, 0.6);
//    double factor;
//
//    count = 0;
//    for (auto& i: {0, 1, 2}) // 1 GeV, 2 GeV, 3 GeV
//    {
//        Nodes[i].GetResultPtr<TH1D>()->SetLineColor(colors[count++]);
//        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "l");
//
//        if (i == 0)
//        {
//            factor = Nodes[i].GetResultPtr<TH1D>()->Integral();
//            Nodes[i].GetResultPtr<TH1D>()->Scale(factor/Nodes[i].GetResultPtr<TH1D>()->Integral());
//            Nodes[i].GetResultPtr<TH1D>()->SetTitle(";#DeltaR_{#gamma#gamma};Events");
//            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);
//            Nodes[i].GetResultPtr<TH1D>()->Draw("HIST");
//        }
//        else
//        {
//            Nodes[i].GetResultPtr<TH1D>()->Scale(factor/Nodes[i].GetResultPtr<TH1D>()->Integral());
//            Nodes[i].GetResultPtr<TH1D>()->Draw("HISTsame");
//            gPad->Modified(); gPad->Update();
//        }
//    }
//
//    gStyle->SetOptStat(0);
//    TLatex Tl;
//    Tl.SetTextSize(0.03);
//    Tl.DrawLatexNDC(0.6, 0.83, "#it{ATLAS} Internal");
//    Tl.DrawLatexNDC(0.6, 0.73,"#sqrt{s} = 13 TeV");
//    legend->SetBorderSize(0);
//    legend->Draw();
//    c1->SaveAs("Fig10A.pdf");
//
//    c1 = new TCanvas();
//    legend = new TLegend(0.65, 0.45, 0.85, 0.65);
//
//    for (auto& i: {2, 3, 4}) // 3 GeV, 5 GeV, 9 GeV
//    {
//        Nodes[i].GetResultPtr<TH1D>()->SetLineColor(colors[count++]);
//        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "l");
//
//        if (i == 2)
//        {
//            factor = Nodes[i].GetResultPtr<TH1D>()->Integral();
//            Nodes[i].GetResultPtr<TH1D>()->Scale(factor/Nodes[i].GetResultPtr<TH1D>()->Integral());
//            Nodes[i].GetResultPtr<TH1D>()->SetTitle(";#DeltaR_{#gamma#gamma};Events");
//            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);
//            Nodes[i].GetResultPtr<TH1D>()->Draw("HIST");
//        }
//        else
//        {
//            Nodes[i].GetResultPtr<TH1D>()->Scale(factor/Nodes[i].GetResultPtr<TH1D>()->Integral());
//            Nodes[i].GetResultPtr<TH1D>()->SetTitle(";#gamma#gamma p_{T} [GeV];Events");
//            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);
//            Nodes[i].GetResultPtr<TH1D>()->Draw("HISTsame");
//            gPad->Modified(); gPad->Update();
//        }
//    }
//    gStyle->SetOptStat(0);
//    Tl.SetTextSize(0.03);
//    Tl.DrawLatexNDC(0.6, 0.83, "#it{ATLAS} Internal");
//    Tl.DrawLatexNDC(0.6, 0.73,"#sqrt{s} = 13 TeV");
//    legend->SetBorderSize(0);
//    legend->Draw();
//    c1->SaveAs("Fig10B.pdf");
//}
//
//
//void fig18()
//{
//    std::vector<std::string> input_filenames =
//    {
//        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root", //1 GeV
//        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600751.PhPy8EG_AZNLO_ggH125_mA2p0_v1.root", // 2 GeV
//        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600752.PhPy8EG_AZNLO_ggH125_mA3p0_v1.root", // 3 GeV
//        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root", // 5 GeV
//    };
//
//    auto findParentInChain = [](int targetBarcode, RVec<TruthParticle>& startParticles, RVec<TruthParticle>& truthChain)
//    {
//        RVec<TruthParticle> truthSelected;
//        bool foundParent;
//        if (truthChain.size() >= 1)
//        {
//            TruthParticle tp;
//            for (auto& tpe: startParticles)
//            {
//                tp = tpe;
//                while (true)
//                {
//                    if (tp.mc_parent_barcode == targetBarcode)
//                    {
//                        truthSelected.push_back(tp);
//                        break;
//                    }
//                    else
//                    {
//                        foundParent = false;
//                        for (auto& tmp: truthChain)
//                        {
//                            if (tp.mc_parent_barcode == tmp.mc_barcode)
//                            {
//                                tp = tmp;
//                                foundParent = true;
//                                break;
//                            }
//                        }
//                        if (foundParent == false)
//                        {
//                            break;
//                        }
//                    }
//                }
//            }
//        }
//        return truthSelected;
//    };
//
//    std::vector<ROOT::RDF::RResultHandle> Nodes;
//
//    std::vector<const char*> prefixes = {"Both photons match", "one photon match", "None match"};
//    std::vector<EColor> colors = {kBlack, kRed, kBlue};
//    int count = 0;
//
//    for (auto& file: input_filenames)
//    {
//        SchottDataFrame df(MakeRDF({file}, 8));
//
//        auto preselection = df.Filter(
//        [&](const RVec<std::string>& trigger_passed_triggers, RVec<TruthParticle> truth_particles)
//        {
//
//            //Selecting the leptons and anti-leptons in the event
//            truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
//            [](TruthParticle& x)
//            {
//                //so, remove elements that are not leptons and/or don't satisfy pt > 20 GeV, |η| < 2.37, and 1.37 ≮ |η| ≮ 1.52 and/or don't have mc_status=1
//                return (std::abs(x.mc_pdg_id) != 11 && std::abs(x.mc_pdg_id) != 12 && std::abs(x.mc_pdg_id) != 13 &&
//                        std::abs(x.mc_pdg_id) != 14 && std::abs(x.mc_pdg_id) != 15 && std::abs(x.mc_pdg_id) != 16 &&
//                        std::abs(x.mc_pdg_id) != 17 && std::abs(x.mc_pdg_id) != 18)
//                ||
//                //Lepton preselection
//                (!((x.mc_pt/1e3 > 20) && (std::abs(x.mc_eta) < 2.37) &&
//                       (!((1.37 < std::abs(x.mc_eta)) && (std::abs(x.mc_eta) < 1.52)))))
//                ||
//                (x.mc_status != 1);
//
//            }), truth_particles.end());
//
//            //must be exactly two leptons/anti-leptons
//            if (truth_particles.size() != 2)
//            {
//                return false;
//            }
//
//            //must be same flavor
//            if (std::abs(truth_particles[0].mc_pdg_id) != std::abs(truth_particles[1].mc_pdg_id))
//            {
//                return false;
//            }
//
//            if (DeltaR(truth_particles[0].Vector(), truth_particles[1].Vector()) <= 0.01)
//            {
//                return false;
//            }
//
//            if (truth_particles[0].mc_charge*truth_particles[1].mc_charge >= 0)
//            {
//                return false;
//            }
//
//            if (!((truth_particles[0].mc_pt > 27e3 && truth_particles[1].mc_pt > 20e3)
//                                                   ||
//                  (truth_particles[1].mc_pt > 27e3 && truth_particles[0].mc_pt > 20e3)))
//            {
//                return false;
//            }
//
//            PtEtaPhiEVector dilepton = truth_particles[0].Vector() + truth_particles[1].Vector();
//
//            if ((dilepton.M() < 81e3) || (dilepton.M() > 101e3))
//            {
//                return false;
//            }
//
//            if ((truth_particles[0].Vector() + truth_particles[1].Vector()).Pt() <= 10e3)
//            {
//                return false;
//            }
//
//            return true;
//
//        }, {"trigger_passed_triggers", "truth_particles"});
//
//        auto truth_photons_from_axions = preselection.Define("truth_photons",[&](RVec<TruthParticle> truth_particles)
//        {
//            truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
//            [](TruthParticle& x)
//            {
//                return (std::abs(x.mc_pdg_id) != 22);
//
//            }), truth_particles.end());
//
//            return truth_particles;
//
//        }, {"truth_particles"})
//        .Define("truth_axions", [&](RVec<TruthParticle> truth_particles)
//        {
//            truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
//            [](TruthParticle& x)
//            {
//                return (std::abs(x.mc_pdg_id) != 35);
//
//            }), truth_particles.end());
//
//            return truth_particles;
//
//        }, {"truth_particles"}).Define("truth_photons_from_axions",
//        [&](RVec<TruthParticle>& truth_photons, RVec<TruthParticle>& truth_axions)
//        {
//            return findParentInChain(truth_axions[0].mc_barcode, truth_photons, truth_axions);
//
//        }, {"truth_photons", "truth_axions"})
//        .Filter([&] (RVec<TruthParticle>& truth_photons_from_axions)
//        {
//            return truth_photons_from_axions.size()==2;
//        }, {"truth_photons_from_axions"})
//        .Define("photons_pass_cut_indices",
//        [&](const RVec<Photon>& photons)
//        {
//            RVec<int> photon_indices;
//            photon_indices.reserve(photons.size());
//
//            for (int i = 0; i < photons.size(); i++)
//            {
//                if (
//                    (std::abs(photons[i].photon_eta) >= 2.37) or
//                    (photons[i].photon_pt <= 10e3) or
//                    (std::abs(photons[i].photon_eta) > 1.37 and std::abs(photons[i].photon_eta) < 1.52) or
//                    (not photons[i].photon_id_loose)
//                    )
//                {
//                    continue;
//                }
//                photon_indices.push_back(i);
//            }
//
//            return photon_indices;
//        }, {"photons"})
//        .Define("photons_pass_cuts",
//        [&](RVec<Photon>& photons, RVec<int>& photon_indices)
//        {
//            return Take(photons, photon_indices);
//        }, {"photons", "photons_pass_cut_indices"});
//
//        auto reco_photons_matched = truth_photons_from_axions
//        .Define("reco_photons_matched_indices",
//        [&](RVec<Photon>& photons_pass_cuts, RVec<TruthParticle>& truth_photons_from_axions)
//        {
//            RVec<int> reco_photons_matched_indices;
//            PtEtaPhiEVector tp1 = truth_photons_from_axions[0].Vector();
//            PtEtaPhiEVector tp2 = truth_photons_from_axions[1].Vector();
//
//            for (int i = 0; i < photons_pass_cuts.size(); i++)
//            {
//                if (not (DeltaR(photons_pass_cuts[i].Vector(), tp1) < 0.1 or DeltaR(photons_pass_cuts[i].Vector(), tp2) < 0.1))
//                {
//                    continue;
//                }
//                reco_photons_matched_indices.push_back(i);
//            }
//
//            return reco_photons_matched_indices;
//
//        }, {"photons_pass_cuts", "truth_photons_from_axions"})
//        .Define("reco_photons_matched",
//        [&](RVec<Photon>& photons_pass_cuts, RVec<int>& reco_photons_matched_indices)
//        {
//            return Take(photons_pass_cuts, reco_photons_matched_indices);
//
//        }, {"photons_pass_cuts", "reco_photons_matched_indices"})
//        .Define("totEventWeight",
//        [](RVec<float> photon_id_eff, RVec<float> photon_iso_eff, RVec<float> photon_trg_eff/*, RVec<float> ei_event_weights_generator*/, RVec<int>& photons_pass_cut_indices, RVec<int>& reco_photons_matched_indices)
//        {
//            //Taking the values that correspond to the
//            auto ResizeVal = std::max({photon_id_eff.size(), photon_iso_eff.size(), photon_trg_eff.size()});
//            photon_id_eff.resize(ResizeVal,1);
//            photon_iso_eff.resize(ResizeVal,1);
//            photon_trg_eff.resize(ResizeVal,1);
//
//            return
//            Take(Take(photon_id_eff, photons_pass_cut_indices), reco_photons_matched_indices)*
//            Take(Take(photon_iso_eff, photons_pass_cut_indices), reco_photons_matched_indices)*
//            Take(Take(photon_trg_eff, photons_pass_cut_indices), reco_photons_matched_indices);//*ei_event_weights_generator[0];
//
//        }, {"photon_id_eff", "photon_iso_eff", "photon_trg_eff"/*, "ei_event_weights_generator"*/, "photons_pass_cut_indices", "reco_photons_matched_indices"});
//
//        auto two_reco_photons_matched = reco_photons_matched.Filter(
//        [&](RVec<Photon>& reco_photons_matched)
//        {
//            return reco_photons_matched.size()==2;
//
//        }, {"reco_photons_matched"})
//        .Define("truth_photons_from_axions_pt",
//        [&](RVec<TruthParticle>& truth_photons_from_axions)
//        {
//            PtEtaPhiEVector four_momentum = truth_photons_from_axions[0].Vector() + truth_photons_from_axions[1].Vector();
//            return four_momentum.Pt()/1e3;
//        }, {"truth_photons_from_axions"});
//
//        auto one_reco_photons_matched = reco_photons_matched.Filter(
//        [&](RVec<Photon>& reco_photons_matched)
//        {
//            return reco_photons_matched.size()==1;
//
//        }, {"reco_photons_matched"})
//        .Define("truth_photons_from_axions_pt",
//        [&](RVec<TruthParticle>& truth_photons_from_axions)
//        {
//            PtEtaPhiEVector four_momentum = truth_photons_from_axions[0].Vector() + truth_photons_from_axions[1].Vector();
//            return four_momentum.Pt()/1e3;
//        }, {"truth_photons_from_axions"});
//
//        auto zero_reco_photons_matched = reco_photons_matched.Filter(
//        [&](RVec<Photon>& reco_photons_matched)
//        {
//            return reco_photons_matched.empty();
//
//        }, {"reco_photons_matched"})
//        .Define("truth_photons_from_axions_pt",
//        [&](RVec<TruthParticle>& truth_photons_from_axions)
//        {
//            PtEtaPhiEVector four_momentum = truth_photons_from_axions[0].Vector() + truth_photons_from_axions[1].Vector();
//            return four_momentum.Pt()/1e3;
//        }, {"truth_photons_from_axions"});
//
//        Nodes.push_back(two_reco_photons_matched.Histo1D<double>({prefixes[0], prefixes[0], 100u, 0, 202}, "truth_photons_from_axions_pt"));
//        Nodes.push_back(one_reco_photons_matched.Histo1D<double>({prefixes[1], prefixes[1], 100u, 0, 202}, "truth_photons_from_axions_pt"));
//        Nodes.push_back(zero_reco_photons_matched.Histo1D<double>({prefixes[2], prefixes[2], 100u, 0, 202}, "truth_photons_from_axions_pt"));
//
//        count++;
//    }
//
//    ROOT::RDF::RunGraphs(Nodes); // running all computation nodes concurrently
//
//    TCanvas* c1 = new TCanvas();
//    TLegend* legend = new TLegend(0.55, 0.25, 0.85, 0.55);
//
//    count = 0;
//    for (int i = 0; i <= 2; i++)
//    {
//        Nodes[i].GetResultPtr<TH1D>()->SetLineColor(colors[count++]);
//        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "l");
//
//        if (i == 0)
//        {
//            Nodes[i].GetResultPtr<TH1D>()->SetTitle(";p_{T_{#gamma#gamma}} [GeV];Events");
//            Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->SetTitleSize(1.1 * Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->GetTitleSize());
//            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->SetTitleSize(1.25 * Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->GetTitleSize());
//            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);
//            Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->SetTitleOffset(1.1);
//            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->SetTitleOffset(1.1);
//            Nodes[i].GetResultPtr<TH1D>()->SetAxisRange(0., 2870,"Y");
//            Nodes[i].GetResultPtr<TH1D>()->Draw("HIST");
//        }
//        else
//        {
//            Nodes[i].GetResultPtr<TH1D>()->Draw("HISTsame");
//            gPad->Modified(); gPad->Update();
//        }
//        std::cout << i << ' ' << Nodes[i].GetResultPtr<TH1D>()->GetMean() << '\n';
//    }
//
//    gStyle->SetOptStat(0);
//    TLatex Tl;
//    Tl.SetTextSize(0.03);
//    Tl.DrawLatexNDC(0.6, 0.8, "#it{ATLAS} Internal");
//    Tl.DrawLatexNDC(0.6, 0.7,"#sqrt{s} = 13 TeV");
//    Tl.DrawLatexNDC(0.6, 0.6,"m_{A} = 1 GeV, #DeltaR = 0.1");
//    legend->SetBorderSize(0);
//    legend->Draw();
//    c1->SaveAs("Fig18A.pdf");
//
//    c1 = new TCanvas();
//    legend = new TLegend(0.55, 0.25, 0.85, 0.55);
//    count = 0;
//    for (int i = 3; i <= 5; i++)
//    {
//        Nodes[i].GetResultPtr<TH1D>()->SetLineColor(colors[count++]);
//        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "l");
//
//        if (i == 3)
//        {
//            Nodes[i].GetResultPtr<TH1D>()->SetTitle(";p_{T_{#gamma#gamma}} [GeV];Events");
//            Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->SetTitleSize(1.1 * Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->GetTitleSize());
//            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->SetTitleSize(1.25 * Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->GetTitleSize());
//            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);
//            Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->SetTitleOffset(1.1);
//            Nodes[i].GetResultPtr<TH1D>()->SetAxisRange(0., 850,"Y");
//            Nodes[i].GetResultPtr<TH1D>()->Draw("HIST");
//        }
//        else
//        {
//            Nodes[i].GetResultPtr<TH1D>()->SetTitle(";p_{T_{#gamma#gamma}} [GeV];Events");
//            Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->SetTitleSize(1.1 * Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->GetTitleSize());
//            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->SetTitleSize(1.25 * Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->GetTitleSize());
//            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);
//            Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->SetTitleOffset(1.1);
//            Nodes[i].GetResultPtr<TH1D>()->SetAxisRange(0., 850,"Y");
//            Nodes[i].GetResultPtr<TH1D>()->Draw("HISTsame");
//            gPad->Modified(); gPad->Update();
//        }
//        std::cout << i << ' ' << Nodes[i].GetResultPtr<TH1D>()->GetMean() << '\n';
//    }
//    gStyle->SetOptStat(0);
//    Tl.SetTextSize(0.03);
//    Tl.DrawLatexNDC(0.6, 0.8, "#it{ATLAS} Internal");
//    Tl.DrawLatexNDC(0.6, 0.7,"#sqrt{s} = 13 TeV");
//    Tl.DrawLatexNDC(0.6, 0.6,"m_{A} = 2 GeV, #DeltaR = 0.1");
//    legend->SetBorderSize(0);
//    legend->Draw();
//    c1->SaveAs("Fig18B.pdf");
//
//    c1 = new TCanvas();
//    legend = new TLegend(0.55, 0.25, 0.85, 0.55);
//    count = 0;
//    for (int i = 6; i <= 8; i++)
//    {
//        Nodes[i].GetResultPtr<TH1D>()->SetLineColor(colors[count++]);
//        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "l");
//
//        if (i == 6)
//        {
//            Nodes[i].GetResultPtr<TH1D>()->SetTitle(";p_{T_{#gamma#gamma}} [GeV];Events");
//            Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->SetTitleSize(1.1 * Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->GetTitleSize());
//            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->SetTitleSize(1.25 * Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->GetTitleSize());
//            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);
//            Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->SetTitleOffset(1.1);
//            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->SetTitleOffset(1.1);
//            Nodes[i].GetResultPtr<TH1D>()->SetAxisRange(0., 1850,"Y");
//            Nodes[i].GetResultPtr<TH1D>()->Draw("HIST");
//        }
//        else
//        {
//            Nodes[i].GetResultPtr<TH1D>()->SetTitle(";p_{T_{#gamma#gamma}} [GeV];Events");
//            Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->SetTitleSize(1.1 * Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->GetTitleSize());
//            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->SetTitleSize(1.25 * Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->GetTitleSize());
//            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);
//            Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->SetTitleOffset(1.1);
//            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->SetTitleOffset(1.1);
//            Nodes[i].GetResultPtr<TH1D>()->SetAxisRange(0., 1850,"Y");
//            Nodes[i].GetResultPtr<TH1D>()->Draw("HISTsame");
//            gPad->Modified(); gPad->Update();
//        }
//        std::cout << i << ' ' << Nodes[i].GetResultPtr<TH1D>()->GetMean() << '\n';
//    }
//    gStyle->SetOptStat(0);
//    Tl.SetTextSize(0.03);
//    Tl.DrawLatexNDC(0.6, 0.8, "#it{ATLAS} Internal");
//    Tl.DrawLatexNDC(0.6, 0.7,"#sqrt{s} = 13 TeV");
//    Tl.DrawLatexNDC(0.6, 0.6,"m_{A} = 3 GeV, #DeltaR = 0.1");
//    legend->SetBorderSize(0);
//    legend->Draw();
//    c1->SaveAs("Fig18C.pdf");
//
//    c1 = new TCanvas();
//    legend = new TLegend(0.55, 0.25, 0.85, 0.55);
//    count = 0;
//    for (int i = 9; i <= 11; i++)
//    {
//        Nodes[i].GetResultPtr<TH1D>()->SetLineColor(colors[count++]);
//        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "l");
//
//        if (i == 9)
//        {
//            Nodes[i].GetResultPtr<TH1D>()->SetTitle(";p_{T_{#gamma#gamma}} [GeV];Events");
//            Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->SetTitleSize(1.1 * Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->GetTitleSize());
//            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->SetTitleSize(1.25 * Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->GetTitleSize());
//            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);
//            Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->SetTitleOffset(1.1);
//            Nodes[i].GetResultPtr<TH1D>()->SetAxisRange(0., 950,"Y");
//            Nodes[i].GetResultPtr<TH1D>()->Draw("HIST");
//        }
//        else
//        {
//            Nodes[i].GetResultPtr<TH1D>()->SetTitle(";p_{T_{#gamma#gamma}} [GeV];Events");
//            Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->SetTitleSize(1.1 * Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->GetTitleSize());
//            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->SetTitleSize(1.25 * Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->GetTitleSize());
//            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);
//            Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->SetTitleOffset(1.1);
//            Nodes[i].GetResultPtr<TH1D>()->SetAxisRange(0., 950,"Y");
//            Nodes[i].GetResultPtr<TH1D>()->Draw("HISTsame");
//            gPad->Modified(); gPad->Update();
//        }
//        std::cout << i << ' ' << Nodes[i].GetResultPtr<TH1D>()->GetMean() << '\n';
//    }
//    gStyle->SetOptStat(0);
//    Tl.SetTextSize(0.03);
//    Tl.DrawLatexNDC(0.6, 0.8, "#it{ATLAS} Internal");
//    Tl.DrawLatexNDC(0.6, 0.7,"#sqrt{s} = 13 TeV");
//    Tl.DrawLatexNDC(0.6, 0.6,"m_{A} = 5 GeV, #DeltaR = 0.1");
//    legend->SetBorderSize(0);
//    legend->Draw();
//    c1->SaveAs("Fig18D.pdf");
//}

void fig24()
{
    std::vector<std::vector<std::string>> input_filenames = {
        //Signal samples
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600752.PhPy8EG_AZNLO_ggH125_mA3p0_v1.root"}, // 3 GeV
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root"}, //5 GeV
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600756.PhPy8EG_AZNLO_ggH125_mA9p0_v1.root"}, // 9 GeV
        //Z gamma samples
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617070._000001.LGNTuple.root"}, {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617064._000001.LGNTuple.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617074._000001.LGNTuple.root"},
        //Jet samples
        {"/Users/edwardfinkelstein/ATLAS_axion/Jets/Zee_lightJet_0-70.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/Jets/Zee_lightJet_70-140.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/Jets/Zee_lightJet_140-280.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/Jets/Zee_cJet_0-70.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/Jets/Zee_cJet_70-140.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/Jets/Zee_cJet_140-280.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/Jets/Zee_bJet_0-70.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/Jets/Zee_bJet_70-140.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/Jets/Zee_bJet_140-280.root"},
    };

    std::array<double,9> JetNumeratorSFs = {((139e15)*(1.9828e-9)*(0.821204)),((139e15)*(110.64e-12)*(0.69275)),((139e15)*(40.645e-12)*(0.615906)),((139e15)*(1.9817e-9)*(0.1136684)),((139e15)*(110.47e-12)*(0.1912956)),((139e15)*(40.674e-12)*(0.2326772)),((139e15)*(1.9819e-9)*(0.0656969)),((139e15)*(110.53e-12)*(0.1158741)),((139e15)*(40.68e-12)*(0.1535215))};

    std::vector<const char*> prefixes = {"sig m_{A} = 3 GeV", "sig m_{A} = 5 GeV", "sig m_{A} = 9 GeV", "pty2_9_17", "pty_17_myy_0_80", "pty_17_myy_80", "Zee_lightJet_0-70", "Zee_lightJet_70-140", "Zee_lightJet_140-280", "Zee_cJet_0-70", "Zee_cJet_70-140", "Zee_cJet_140-280", "Zee_bJet_0-70", "Zee_bJet_70-140", "Zee_bJet_140-280"};

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

    std::array<double,3> SFs = {((139e15)*(.871e-12)),((139e15)*(.199e-12)), ((139e15)*(.0345e-15))};

    std::vector<EColor> colors = {kBlack, kMagenta, static_cast<EColor>(kRed-7), static_cast<EColor>(kRed-9), static_cast<EColor>(kRed-10)};
    std::vector<EColor> Jetscolors = {static_cast<EColor>(kGreen+0), static_cast<EColor>(kGreen-4), static_cast<EColor>(kGreen-7), static_cast<EColor>(kGreen-9), static_cast<EColor>(kGreen-10), static_cast<EColor>(kGreen+1), static_cast<EColor>(kGreen-3), static_cast<EColor>(kGreen-6), static_cast<EColor>(kGreen-8)};
    int count = 0;

    for (auto& sample: input_filenames)
    {
        SchottDataFrame df(MakeRDF(sample, 8));

        auto EventWeight = df.Define("EventWeight", //mc generator weight
        [](const RVec<float>& ei_event_weights_generator)
        {
            return ((ei_event_weights_generator[0]) ? 1 / ei_event_weights_generator[0] : 1);

        }, {"ei_event_weights_generator"})
        .Define("totEventWeightVec", //vector of efficiencies, product of photon_id_eff, photon_iso_eff, and photon_trg_eff
        [](RVec<float> photon_id_eff, RVec<float> photon_iso_eff, RVec<float> photon_trg_eff/*, RVec<float> ei_event_weights_generator*/)
        {
            auto ResizeVal = std::max({photon_id_eff.size(), photon_iso_eff.size(), photon_trg_eff.size()});
            photon_id_eff.resize(ResizeVal,1);
            photon_iso_eff.resize(ResizeVal,1);
            photon_trg_eff.resize(ResizeVal,1);

            return photon_id_eff;//*photon_iso_eff*photon_trg_eff; //Returns a vector of efficiencies

        }, {"photon_id_eff", "photon_iso_eff", "photon_trg_eff",/* "ei_event_weights_generator"*/});

        //new dataframe node: contains only the events of newDf that pass the trigger cut
        auto trigger_selection = EventWeight
        .Filter([](const RVec<std::string>& trigger_passed_triggers)
        {
            bool trigger_found = (std::find_first_of(trigger_passed_triggers.begin(), trigger_passed_triggers.end(), triggers.begin(), triggers.end()) != trigger_passed_triggers.end());

            if (!trigger_found)
            {
                return false; //this event is filtered out
            }
            return true; //this event is kept because the trigger was found in its `trigger_passed_triggers` branch entry

        }, {"trigger_passed_triggers"});

        //new dataframe node: contains an additional column `di_electrons` and only the events that have 2 electrons and no muons
        auto two_leptons = trigger_selection
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
        [](const RVec<Muon>& muons, RVec<Electron> electrons)
        {
            return (electrons.size()==2 && muons.empty()); //keep events which have exactly 2 electrons for di_electrons and no muons

        }, {"muons", "di_electrons"});

        //new dataframe node: contains only the events from `two_leptons` whose electrons in the `di_electrons` branch have opposite charge
        auto opp_charge = two_leptons.Filter([](const RVec<Electron>& electrons)
        {
            return (electrons[0].electron_charge*electrons[1].electron_charge < 0);

        }, {"di_electrons"});

        //new dataframe node: contains only the events from `opp_charge` that have 1 electron with pt > 20 GeV and the other with pt > 27 GeV
        auto leading_pt = opp_charge.Filter([](const RVec<Electron>& electrons)
        {
            return ((electrons[0].electron_pt > 20e3 && electrons[1].electron_pt > 27e3) || (electrons[1].electron_pt > 20e3 && electrons[0].electron_pt > 27e3));
        }, {"di_electrons"});

        //new dataframe node: nothing changes here from `leading_pt`...
        auto same_flavour = leading_pt.Filter([] (const RVec<Electron>& electrons)
        {
            return true; //because true is always returned in this function
        }, {"di_electrons"});

        //new dataframe node: Contains a new column `dilep` in addition to the ones in `same_flavour` that stores the di-electron four-vector
        auto dilep = same_flavour.Define("dilep",[] (RVec<Electron>& electrons)
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

        //photon acceptance and id_loose cuts
        auto photons_pass_cuts = ptCut
        .Define("photons_pass_cut_indices", //new column that contains the good photon indices
        [&](RVec<Photon>& p) //p = photon
        {
            RVec<int> x; //indices of photons that passed the cuts
            for (auto i = 0; i < p.size(); i++)
            {
                //keep reco-photons that have |η| < 2.37, p_T > 10 GeV, and |η| not between 1.37 and 1.52
                if (not ((std::abs(p[i].photon_eta) >= 2.37) or (p[i].photon_pt <= 10e3) or (std::abs(p[i].photon_eta) > 1.37 and std::abs(p[i].photon_eta) < 1.52)))
                {
                    x.push_back(i);
                }
            }
            return x;

        }, {"photons"})
        .Define("photons_pass_cuts", //new column that contains the good photons corresponding to the good photon indices from above
        [&](RVec<Photon>& photons, RVec<int>& x)
        {
            return Take(photons, x); //Taking only the photons that passed the cuts in each event

        }, {"photons", "photons_pass_cut_indices"});

        auto resolved = photons_pass_cuts.Define("chosen_two_indices",
        [](RVec<Photon>& photons_pass_cuts)
        {
            RVec<unsigned long> x; //vector of indices
            if (photons_pass_cuts.size() < 2)
            {
                return x;
            }

            auto combs = Combinations(photons_pass_cuts, 2); //all combinations of 2 reco-photons
            size_t length = combs[0].size(); //number of combinations
            double delta_r, m, pt, X, best_X, pt1, pt2, chosen_delta_r;

            for (size_t i=0; i<length; i++) //looping through all of the possible combinations of photons in each event
            {
                delta_r = DeltaR(photons_pass_cuts[combs[0][i]].Vector(), photons_pass_cuts[combs[1][i]].Vector());
                m = (photons_pass_cuts[combs[0][i]].Vector() + photons_pass_cuts[combs[1][i]].Vector()).M();
                pt = (photons_pass_cuts[combs[0][i]].Vector() + photons_pass_cuts[combs[1][i]].Vector()).Pt();
                X = delta_r*(pt/(2.0*m));
                //if it's the first combination or if new X is closer to 1
                //than current best_X and ΔR
                //between the two reco-photons < 1.5, then update best_X, pt1, pt2,
                //and the corresponding reco-photon indices x
                if (i==0 || ((std::abs(1-X) < std::abs(1-best_X)) and (delta_r < 1.5)))
                {
                    best_X = X;
                    pt1 = photons_pass_cuts[combs[0][i]].photon_pt;
                    pt2 = photons_pass_cuts[combs[1][i]].photon_pt;
                    chosen_delta_r = delta_r;
                    x = {combs[0][i], combs[1][i]};
                }
            }
            if (chosen_delta_r < 1.5 && pt1 > 10e3 && pt2 > 10e3 && best_X > 0.96 && best_X < 1.2) //two photons corresponding to best_X must both have p_T > 10 GeV, ΔR < 1.5, and 0.96 < best_X < 1.2
            {
                return x;
            }

            x.clear();
            return x;
        }, {"photons_pass_cuts"}).Filter(//keep only events that have passed the first part of the resolved category
        [&](RVec<unsigned long>& indices)
        {
            return (indices.size()==2);

        }, {"chosen_two_indices"})
        .Define("chosen_two", //New column: consists of the good photons corresponding to the `chosen_two_indices` defined above
        [&](RVec<Photon>& reco_photons_matched, RVec<unsigned long>& indices)
        {
            return Take(reco_photons_matched, indices);

        }, {"photons_pass_cuts", "chosen_two_indices"})
        .Define("reco_photons_from_axions_deltaR",
        [&](RVec<Photon>& reco_photons_matched)
        {
            return DeltaR(reco_photons_matched[0].Vector(), reco_photons_matched[1].Vector());

        }, {"chosen_two"})
        .Define("totEventWeight", //New column: weight factor for events in RDF `resolved`
        [&](RVec<unsigned long> chosen_two_indices, RVec<int>& photons_pass_cut_indices, RVec<float>& photon_efficiencies)
        {
            //First, we take the elements from photon_efficiencies that correspond to the `photons_pass_cut_indices` defined
            //earlier. Then, from that resulting vector, we take the elements corresponding to the
            //`chosen_two_indices` defined above
            RVec<float> resolved_photon_efficiencies = Take(Take(photon_efficiencies, photons_pass_cut_indices), chosen_two_indices);
            //Now, multiply all of the elements of the vector we defined above
            float total = 1.0f;
            for (auto i: resolved_photon_efficiencies)
            {
                total *= i;
            }
            //and return the result
            return total;
        }, {"chosen_two_indices", "photons_pass_cut_indices", "totEventWeightVec"});

        Nodes.push_back(resolved.Histo1D<double>({prefixes[count], prefixes[count], 100u, 0, 0.8}, "reco_photons_from_axions_deltaR", "totEventWeight"));
        Nodes.push_back(resolved.Histo1D<double>({prefixes[count], prefixes[count++], 100u, 0, 0.8}, "reco_photons_from_axions_deltaR", "totEventWeight"));
        if (count > 3)
        {
            Nodes.push_back(resolved.Sum<float>("totEventWeight"));
            Nodes.push_back(EventWeight.Sum<float>("EventWeight")); //need total events for scale factor.
        }
    }

//    0       1
//    2       3
//    4       5
//    6       7       8       9
//    10      11      12      13
//    14      15      16      17
//    18      19      20      21
//    22      23      24      25
//    26      27      28      29
//    30      31      32      33
//    34      35      36      37
//    38      39      40      41
//    42      43      44      45
//    46      47      48      49
//    50      51      52      53

    ROOT::RDF::RunGraphs(Nodes); // running all computation nodes concurrently

    double factor;
    double total_back = 0;

    TCanvas* c1 = new TCanvas();
    TLegend* legend = new TLegend(0.4775, 0.43, 0.8, 0.75);

    //Background counts from Z gamma
    for (int i = 8, j = 0; (i <= 16 && j <= 2); i+=4, j++)
    {
        total_back += *Nodes[i].GetResultPtr<float>()*(SFs[j]/ *Nodes[i+1].GetResultPtr<float>());
    }

    //Background counts from Z jets
    for (int i = 20, j = 0; (i <= 52 && j <= 8); i += 4, j++)
    {
        total_back += *Nodes[i].GetResultPtr<float>()*(JetNumeratorSFs[j] / *Nodes[i+1].GetResultPtr<float>());
    }

    int back_count = 0;
    factor = total_back;
    auto hs = new THStack("hs3","");
    //weighted Zgamma
    for (auto& i: {6,10,14})
    {
        Nodes[i].GetResultPtr<TH1D>()->SetFillColor(colors[2+back_count]);
        Nodes[i].GetResultPtr<TH1D>()->SetLineColor(colors[2+back_count]);
        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "f");
        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0)
        {
            Nodes[i].GetResultPtr<TH1D>()->Scale(SFs[back_count++]/ *Nodes[i+3].GetResultPtr<float>());
        }
        Nodes[i].GetResultPtr<TH1D>()->SetTitle(";#DeltaR_{#gamma#gamma} [GeV];Events");
        Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);


        hs->Add(&*Nodes[i].GetResultPtr<TH1D>());
    }

    //weighted ZJets
    for (int i = 18, j = 0; (i <= 50 && j <= 8); i += 4, j++)
    {
        Nodes[i].GetResultPtr<TH1D>()->SetFillColor(Jetscolors[j]);
        Nodes[i].GetResultPtr<TH1D>()->SetLineColor(Jetscolors[j]);
        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "f");
        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0)
        {
            Nodes[i].GetResultPtr<TH1D>()->Scale(JetNumeratorSFs[j] / (*Nodes[i+3].GetResultPtr<float>()));
        }
        Nodes[i].GetResultPtr<TH1D>()->SetTitle(";#DeltaR_{#gamma#gamma} [GeV];Events");
        Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);


        hs->Add(&*Nodes[i].GetResultPtr<TH1D>());
    }

    hs->Draw("HIST");

    count = 0;
    for (auto& i: {0, 2, 4})
    {
        Nodes[i].GetResultPtr<TH1D>()->SetLineColor(colors[count++]);
        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "l");

        if (i == 0)
        {
            if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0)
            {
                Nodes[i].GetResultPtr<TH1D>()->Scale(factor/Nodes[i].GetResultPtr<TH1D>()->Integral());
            }
            Nodes[i].GetResultPtr<TH1D>()->SetTitle(";#DeltaR_{#gamma#gamma} [GeV];Events");
            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);

            Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->SetTitleOffset(1.2);

            Nodes[i].GetResultPtr<TH1D>()->Draw("HISTsame");
        }
        else
        {
            if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0)
            {
                Nodes[i].GetResultPtr<TH1D>()->Scale(factor/Nodes[i].GetResultPtr<TH1D>()->Integral());
            }
            Nodes[i].GetResultPtr<TH1D>()->Draw("HISTsame");
        }
    }

    hs->SetTitle(";#DeltaR_{#gamma#gamma} [GeV];Events");
    hs->GetYaxis()->CenterTitle(true);
    hs->GetXaxis()->SetTitleOffset(1.2);
    hs->GetYaxis()->SetTitleOffset(1.4);
    hs->SetMinimum(0.);
    hs->SetMaximum(24000);

    gStyle->SetOptStat(0);
    TLatex Tl;
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.6, 0.86, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.6, 0.78,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("Fig24B.pdf");

    c1 = new TCanvas();
    legend = new TLegend(0.525, 0.2125, 0.825, 0.6125);

    back_count=0;
    hs = new THStack("hs3","");
    //background Zgamma unweighted
    for (auto& i: {7, 11, 15})
    {
        Nodes[i].GetResultPtr<TH1D>()->SetFillColor(colors[2+back_count]);
        Nodes[i].GetResultPtr<TH1D>()->SetLineColor(colors[2+back_count]);
        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "f");
        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0)
        {
            Nodes[i].GetResultPtr<TH1D>()->Scale(SFs[back_count++]/(*Nodes[i+2].GetResultPtr<float>()));
            gPad->Modified(); gPad->Update();
        }
        Nodes[i].GetResultPtr<TH1D>()->SetTitle(";m_{ll}  [GeV];Events");
        Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);


        hs->Add(&*Nodes[i].GetResultPtr<TH1D>());
        count++;
    }

    for (int i = 19, j = 0; (i <= 51 && j <= 8); i += 4, j++)
    {
        Nodes[i].GetResultPtr<TH1D>()->SetFillColor(Jetscolors[j]);
        Nodes[i].GetResultPtr<TH1D>()->SetLineColor(Jetscolors[j]);
        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "f");
        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0)
        {
            Nodes[i].GetResultPtr<TH1D>()->Scale(JetNumeratorSFs[j] / (*Nodes[i+2].GetResultPtr<float>()));
            gPad->Modified(); gPad->Update();
        }
        Nodes[i].GetResultPtr<TH1D>()->SetTitle(";m_{ll}  [GeV];Events");
        Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);


        hs->Add(&*Nodes[i].GetResultPtr<TH1D>());
    }

    hs->Draw("HIST");
    count = 0;
    //signal
    for (auto& i: {1, 3, 5})
    {
        Nodes[i].GetResultPtr<TH1D>()->SetLineColor(colors[count++]);
        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "l");

        if (i == 1)
        {
            Nodes[i].GetResultPtr<TH1D>()->SetTitle(";#DeltaR_{#gamma#gamma} [GeV];Events");
            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);

            Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->SetTitleOffset(1.2);
            Nodes[i].GetResultPtr<TH1D>()->Draw("HISTsame");
            gPad->Modified(); gPad->Update();
        }
        else
        {
            Nodes[i].GetResultPtr<TH1D>()->SetTitle(";#DeltaR_{#gamma#gamma} [GeV];Events");
            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);

            Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->SetTitleOffset(1.2);
            Nodes[i].GetResultPtr<TH1D>()->Draw("HISTsame");
            gPad->Modified(); gPad->Update();
        }
    }
    hs->SetTitle(";#DeltaR_{#gamma#gamma} [GeV];Events");
    hs->GetYaxis()->CenterTitle(true);
    hs->GetXaxis()->SetTitleOffset(1.2);
    hs->GetYaxis()->SetTitleOffset(1.4);
    hs->SetMinimum(0.);
    hs->SetMaximum(24000);

    gStyle->SetOptStat(0);
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.6, 0.8, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.6, 0.7,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("Fig24A.pdf");
}

//void fig54()
//{
//    std::vector<std::string> input_filenames = {
//        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root", //1 GeV
//        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600752.PhPy8EG_AZNLO_ggH125_mA3p0_v1.root", //2 GeV
//        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600756.PhPy8EG_AZNLO_ggH125_mA9p0_v1.root" //3 GeV
//    };
//
//    auto findParentInChain = [](int targetBarcode, RVec<TruthParticle>& startParticles, RVec<TruthParticle>& truthChain)
//    {
//        RVec<TruthParticle> truthSelected;
//        bool foundParent;
//        if (truthChain.size() >= 1)
//        {
//            TruthParticle tp;
//            for (auto& tpe: startParticles)
//            {
//                tp = tpe;
//                while (true)
//                {
//                    if (tp.mc_parent_barcode == targetBarcode)
//                    {
//                        truthSelected.push_back(tp);
//                        break;
//                    }
//                    else
//                    {
//                        foundParent = false;
//                        for (auto& tmp: truthChain)
//                        {
//                            if (tp.mc_parent_barcode == tmp.mc_barcode)
//                            {
//                                tp = tmp;
//                                foundParent = true;
//                                break;
//                            }
//                        }
//                        if (foundParent == false)
//                        {
//                            break;
//                        }
//                    }
//                }
//            }
//        }
//        return truthSelected;
//    };
//
//    std::vector<ROOT::RDF::RResultHandle> Nodes;
//
//    std::vector<const char*> prefixes = {"Sig m_{A} = 1 GeV", "Sig m_{A} = 2 GeV", "Sig m_{A} = 3 GeV"};
//    std::vector<EColor> colors = {kRed, kBlue, kMagenta};
//    int count = 0;
//
//    for (auto& file: input_filenames)
//    {
//        SchottDataFrame df(MakeRDF({file}, 8));
//
//        auto preselection = df.Filter(
//        [&](const RVec<std::string>& trigger_passed_triggers, RVec<TruthParticle> truth_particles)
//        {
//
//            //Selecting the leptons and anti-leptons in the event
//            truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
//            [](TruthParticle& x)
//            {
//                //so, remove elements that are not leptons and/or don't satisfy pt > 20 GeV, |η| < 2.37, and 1.37 ≮ |η| ≮ 1.52 and/or don't have mc_status=1
//                return (std::abs(x.mc_pdg_id) != 11 && std::abs(x.mc_pdg_id) != 12 && std::abs(x.mc_pdg_id) != 13 &&
//                      std::abs(x.mc_pdg_id) != 14 && std::abs(x.mc_pdg_id) != 15 && std::abs(x.mc_pdg_id) != 16 &&
//                      std::abs(x.mc_pdg_id) != 17 && std::abs(x.mc_pdg_id) != 18)
//                ||
//                //Lepton preselection
//                (!((x.mc_pt/1e3 > 20) && (std::abs(x.mc_eta) < 2.37) &&
//                 (!((1.37 < std::abs(x.mc_eta)) && (std::abs(x.mc_eta) < 1.52)))))
//                ||
//                (x.mc_status != 1);
//
//            }), truth_particles.end());
//
//            //must be exactly two leptons/anti-leptons
//            if (truth_particles.size() != 2)
//            {
//                return false;
//            }
//
//            //must be same flavor
//            if (std::abs(truth_particles[0].mc_pdg_id) != std::abs(truth_particles[1].mc_pdg_id))
//            {
//                return false;
//            }
//
//            if (DeltaR(truth_particles[0].Vector(), truth_particles[1].Vector()) <= 0.01)
//            {
//                return false;
//            }
//
//            if (truth_particles[0].mc_charge*truth_particles[1].mc_charge >= 0)
//            {
//                return false;
//            }
//
//            if (!((truth_particles[0].mc_pt > 27e3 && truth_particles[1].mc_pt > 20e3)
//            ||
//            (truth_particles[1].mc_pt > 27e3 && truth_particles[0].mc_pt > 20e3)))
//            {
//                return false;
//            }
//
//            PtEtaPhiEVector dilepton = truth_particles[0].Vector() + truth_particles[1].Vector();
//
//            if ((dilepton.M() < 81e3) || (dilepton.M() > 101e3))
//            {
//                return false;
//            }
//
//            if ((truth_particles[0].Vector() + truth_particles[1].Vector()).Pt() <= 10e3)
//            {
//                return false;
//            }
//
//            return true;
//
//        }, {"trigger_passed_triggers", "truth_particles"});
//
//        auto truth_photons_from_axions = preselection.Define("truth_photons",[&](RVec<TruthParticle> truth_particles)
//        {
//            truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
//            [](TruthParticle& x)
//            {
//                return (std::abs(x.mc_pdg_id) != 22);
//
//            }), truth_particles.end());
//
//            return truth_particles;
//
//        }, {"truth_particles"})
//        .Define("truth_axions", [&](RVec<TruthParticle> truth_particles)
//        {
//            truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
//            [](TruthParticle& x)
//            {
//                return (std::abs(x.mc_pdg_id) != 35);
//
//            }), truth_particles.end());
//
//            return truth_particles;
//
//        }, {"truth_particles"}).Define("truth_photons_from_axions",
//        [&](RVec<TruthParticle>& truth_photons, RVec<TruthParticle>& truth_axions)
//        {
//            return findParentInChain(truth_axions[0].mc_barcode, truth_photons, truth_axions);
//
//        }, {"truth_photons", "truth_axions"})
//        .Filter([&] (RVec<TruthParticle>& truth_photons_from_axions)
//        {
//            return truth_photons_from_axions.size()==2;
//        }, {"truth_photons_from_axions"})
//        .Define("photons_pass_cut_indices",
//        [&](const RVec<Photon>& photons)
//        {
//            RVec<int> photon_indices;
//            photon_indices.reserve(photons.size());
//
//            for (int i = 0; i < photons.size(); i++)
//            {
//                if (
//                (std::abs(photons[i].photon_eta) >= 2.37) or
//                (photons[i].photon_pt <= 10e3) or
//                (std::abs(photons[i].photon_eta) > 1.37 and std::abs(photons[i].photon_eta) < 1.52) or
//                (not photons[i].photon_id_loose)
//                )
//                {
//                    continue;
//                }
//                photon_indices.push_back(i);
//            }
//
//            return photon_indices;
//        }, {"photons"})
//        .Define("photons_pass_cuts",
//        [&](RVec<Photon>& photons, RVec<int>& photon_indices)
//        {
//            return Take(photons, photon_indices);
//        }, {"photons", "photons_pass_cut_indices"});
//
//        auto reco_photons_matched = truth_photons_from_axions
//        .Define("reco_photons_matched_indices",
//        [&](RVec<Photon>& photons_pass_cuts, RVec<TruthParticle>& truth_photons_from_axions)
//        {
//            RVec<int> reco_photons_matched_indices;
//            PtEtaPhiEVector tp1 = truth_photons_from_axions[0].Vector();
//            PtEtaPhiEVector tp2 = truth_photons_from_axions[1].Vector();
//
//            for (int i = 0; i < photons_pass_cuts.size(); i++)
//            {
//                if (not (DeltaR(photons_pass_cuts[i].Vector(), tp1) < 0.1 or DeltaR(photons_pass_cuts[i].Vector(), tp2) < 0.1))
//                {
//                    continue;
//                }
//                reco_photons_matched_indices.push_back(i);
//            }
//
//            return reco_photons_matched_indices;
//
//        }, {"photons_pass_cuts", "truth_photons_from_axions"})
//        .Define("reco_photons_matched",
//        [&](RVec<Photon>& photons_pass_cuts, RVec<int>& reco_photons_matched_indices)
//        {
//            return Take(photons_pass_cuts, reco_photons_matched_indices);
//
//        }, {"photons_pass_cuts", "reco_photons_matched_indices"});
//
//        //New dataframe node: contains only events from `photon_passes_cuts` that fail the resolved category and pass the merged category
//        auto merged_reco_photons_matched = reco_photons_matched.Filter(
//        [&](RVec<Photon>& reco_photons_matched)
//        {
//           if (reco_photons_matched.size() == 1) // 1 photon in the event
//           {
//               return reco_photons_matched[0].photon_pt > 20e3; //event passes if photon pt > 20 GeV
//           }
//           else if (reco_photons_matched.empty())
//           {
//               return false; //fails if no photons in event
//           }
//
//           auto combs = Combinations(reco_photons_matched, 2); //all combinations of 2 reco-photons
//           size_t length = combs[0].size(); //number of combinations
//           double delta_r, m, pt, X, best_X, pt1, pt2, chosen_delta_r;
//
//           for (size_t i=0; i<length; i++) //looping through all of the possible combinations of photons in each event
//           {
//               delta_r = DeltaR(reco_photons_matched[combs[0][i]].Vector(), reco_photons_matched[combs[1][i]].Vector());
//               m = (reco_photons_matched[combs[0][i]].Vector() + reco_photons_matched[combs[1][i]].Vector()).M();
//               pt = (reco_photons_matched[combs[0][i]].Vector() + reco_photons_matched[combs[1][i]].Vector()).Pt();
//               X = delta_r*(pt/(2.0*m));
//               //if it's the first combination or if new X is closer to 1
//               //than current best_X and ΔR
//               //between the two reco-photons < 1.5, then update best_X, pt1, pt2,
//               //and the corresponding reco-photon indices x
//               if (i==0 || ((std::abs(1-X) < std::abs(1-best_X)) and (delta_r < 1.5)))
//               {
//                   best_X = X;
//                   pt1 = reco_photons_matched[combs[0][i]].photon_pt;
//                   pt2 = reco_photons_matched[combs[1][i]].photon_pt;
//                   chosen_delta_r = delta_r;
//               }
//           }
//           if (pt1 > 10e3 && pt2 > 10e3 && best_X > 0.96 && best_X < 1.2 && chosen_delta_r < 1.5) //two photons corresponding to best_X must both have p_T > 10 GeV, ΔR < 1.5, and 0.96 < best_X < 1.2 for resolved... but this is merged, so if it passes resolved it fails merged
//           {
//               return false;
//           }
//           //if we get to this point, it means we've failed resolved
//           for (auto& p: reco_photons_matched) //merged means at least 1 reco photon must have pt > 20 GeV
//           {
//               if (p.photon_pt > 20e3)
//               {
//                   return true; //passed merged if there's a reco-photon with pt > 20 GeV
//               }
//           }
//           return false; //failed merged
//
//        }, {"photons_pass_cuts"})
//        .Define("merged_photon_index", //new column: consists of the index corresponding to the photon that made the event be classified as merged
//        [&](const RVec<Photon>& rpm) //rpm = reco photons matched
//        {
//            for (auto i = 0; i < rpm.size(); i++)
//            {
//                if (rpm[i].photon_pt > 20e3)
//                {
//                    return i; //returning the index of the first photon that has photon_pt > 20 GeV
//                }
//            }
//            return 0; //jic the compiler complains, should not come to this
//
//        }, {"photons_pass_cuts"})
//        .Define("merged_photon", //new column: The reco-photon corresponding to `merged_photon_index`
//        [&](const RVec<Photon>& reco_photons_matched, int merged_photon_index)
//        {
//            return reco_photons_matched[merged_photon_index];
//
//        }, {"photons_pass_cuts", "merged_photon_index"})
//        //mpi = merged_photon_index
//        .Define("totEventWeight", //New column: weight factor for events in RDF `merged_reco_photons_matched`
//        [](RVec<int>& photons_pass_cut_indices, /*RVec<int>& reco_photons_matched_indices,*/ int mpi, RVec<float> photon_id_eff, RVec<float> photon_iso_eff, RVec<float> photon_trg_eff/*, RVec<float> ei_event_weights_generator*/)
//        {
//            //   ||  resize 3 vectors just in case they don't already have the same size  ||
//            //   \/                                                                       \/
//            auto ResizeVal = std::max({photon_id_eff.size(), photon_iso_eff.size(), photon_trg_eff.size()});
//            photon_id_eff.resize(ResizeVal,1);
//            photon_iso_eff.resize(ResizeVal,1);
//            photon_trg_eff.resize(ResizeVal,1);
//
//            //First, take the photons_pass_cut_indices indices
//            //from the respective photon_*_eff vectors,
//            //this corresponds to the photons from the set of all reco photons in
//            //the event that passed the "photon_passes_cuts" cuts. Then, from
//            //those photons, select the index mpi element that corresponds to
//            //the merged reco-photon
//
//            return Take(photon_id_eff, photons_pass_cut_indices)[mpi];
////            * Take(photon_iso_eff, photons_pass_cut_indices)[mpi]
////            * Take(photon_trg_eff, photons_pass_cut_indices)[mpi]; // a single number
//
//        }, {"photons_pass_cut_indices", /*"reco_photons_matched_indices",*/ "merged_photon_index", "photon_id_eff", "photon_iso_eff", "photon_trg_eff",/* "ei_event_weights_generator"*/});
//
//        auto stable_truth_dilepton_and_photon = merged_reco_photons_matched
//        .Define("stable_truth_lepton_indices",[&](RVec<TruthParticle>& truth_particles)
//        {
//            RVec<int> stable_truth_lepton_indices;
//
//            for (int i = 0; i < truth_particles.size(); i++)
//            {
//                if ((std::abs(truth_particles[i].mc_pdg_id) != 11 && std::abs(truth_particles[i].mc_pdg_id) != 12 && std::abs(truth_particles[i].mc_pdg_id) != 13 &&
//                     std::abs(truth_particles[i].mc_pdg_id) != 14 && std::abs(truth_particles[i].mc_pdg_id) != 15 && std::abs(truth_particles[i].mc_pdg_id) != 16 &&
//                     std::abs(truth_particles[i].mc_pdg_id) != 17 && std::abs(truth_particles[i].mc_pdg_id) != 18)
//                    ||
//                    (!((truth_particles[i].mc_pt/1e3 > 20) && (std::abs(truth_particles[i].mc_eta) < 2.37) &&
//                       (!((1.37 < std::abs(truth_particles[i].mc_eta)) && (std::abs(truth_particles[i].mc_eta) < 1.52)))))
//                    ||
//                    (truth_particles[i].mc_status != 1))
//                {
//                    continue;
//                }
//                stable_truth_lepton_indices.push_back(i);
//            }
//            return stable_truth_lepton_indices;
//
//        }, {"truth_particles"})
//        .Define("stable_truth_leptons",
//        [&](RVec<int>& stable_truth_lepton_indices, RVec<TruthParticle>& truth_particles)
//        {
//            return Take(truth_particles, stable_truth_lepton_indices);
//        }, {"stable_truth_lepton_indices", "truth_particles"})
//        .Filter([&](RVec<TruthParticle>& stable_truth_leptons, RVec<TruthParticle>& truth_photons_from_axions)
//        {
//            return (stable_truth_leptons.size()==2 && truth_photons_from_axions.size()==2);
//
//        }, {"stable_truth_leptons", "truth_photons_from_axions"})
//        .Define("reconstructed_mass",[&](RVec<TruthParticle>& stable_truth_leptons, Photon& merged_photon)
//        {
//            auto four_momentum = stable_truth_leptons[0].Vector() + stable_truth_leptons[1].Vector();
//
//            return (four_momentum + merged_photon.Vector()).M()/1e3;
//
//        }, {"stable_truth_leptons", "merged_photon"});
//
//        auto preSB = stable_truth_dilepton_and_photon.Filter(
//        [](double reconstructed_mass)
//        {
//            return ((reconstructed_mass <= 110) || (reconstructed_mass >= 130));
//        }, {"reconstructed_mass"});
//
//        auto SR = stable_truth_dilepton_and_photon.Filter(
//        [](double reconstructed_mass, RVec<float>& Eratio)
//        {
//            return ((reconstructed_mass > 110) && (reconstructed_mass < 130) && (!Any(Eratio <= 0.8)));
//        }, {"reconstructed_mass", "photon_shower_shape_e_ratio"})
//        .Define("reconstructed_pt",[&](RVec<TruthParticle>& stable_truth_leptons, Photon& merged_photon)
//        {
//            auto four_momentum = stable_truth_leptons[0].Vector() + stable_truth_leptons[1].Vector();
//
//            return (four_momentum + merged_photon.Vector()).Pt()/1e3;
//
//        }, {"stable_truth_leptons", "merged_photon"})
//        .Define("reconstructed_deltaR",[&](RVec<TruthParticle>& stable_truth_leptons, Photon& merged_photon)
//        {
//            auto four_momentum = stable_truth_leptons[0].Vector() + stable_truth_leptons[1].Vector();
//
//            return DeltaR(four_momentum,merged_photon.Vector());
//
//        }, {"stable_truth_leptons", "merged_photon"});
//
//        Nodes.push_back(stable_truth_dilepton_and_photon.Histo1D<double>({prefixes[count], prefixes[count], 100u, 80, 200}, "reconstructed_mass", "totEventWeight"));
//        Nodes.push_back(SR.Histo1D<double>({prefixes[count], prefixes[count], 100u, 0, 4.7}, "reconstructed_deltaR", "totEventWeight"));
//        Nodes.push_back(SR.Histo1D<double>({prefixes[count], prefixes[count], 100u, 0, 220}, "reconstructed_pt", "totEventWeight"));
//        Nodes.push_back(preSB.Histo1D<RVec<float>>({prefixes[count], prefixes[count++], 100u, 0, 1}, "photon_shower_shape_e_ratio", "totEventWeight")); //stable_truth_dilepton_and_photon -> preSB ...
//
//    }
//
//    ROOT::RDF::RunGraphs(Nodes); // running all computation nodes concurrently
//
//    TCanvas* c1 = new TCanvas();
//    TLegend* legend = new TLegend(0.65, 0.4, 0.85, 0.6);
//    double factor;
//
//    count = 0;
//    for (auto& i: {0,4,8})
//    {
//        Nodes[i].GetResultPtr<TH1D>()->SetLineColor(colors[count++]);
//        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "l");
//
//        if (i == 0)
//        {
//            factor = Nodes[i].GetResultPtr<TH1D>()->Integral();
//            Nodes[i].GetResultPtr<TH1D>()->Scale(factor/Nodes[i].GetResultPtr<TH1D>()->Integral());
//            Nodes[i].GetResultPtr<TH1D>()->SetTitle(";m_{ll#gamma} [GeV];Events");
//            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);
//            Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->SetTitleOffset(1.2);
//            Nodes[i].GetResultPtr<TH1D>()->SetAxisRange(0., 2750,"Y");
//            Nodes[i].GetResultPtr<TH1D>()->Draw("HIST");
//        }
//        else
//        {
//            Nodes[i].GetResultPtr<TH1D>()->Scale(factor/Nodes[i].GetResultPtr<TH1D>()->Integral());
//            Nodes[i].GetResultPtr<TH1D>()->Draw("HISTsame");
//            gPad->Modified(); gPad->Update();
//        }
//    }
//
//    gStyle->SetOptStat(0);
//    TLatex Tl;
//    Tl.SetTextSize(0.03);
//    Tl.DrawLatexNDC(0.6, 0.8, "#it{ATLAS} Internal");
//    Tl.DrawLatexNDC(0.6, 0.7,"#sqrt{s} = 13 TeV");
//    legend->SetBorderSize(0);
//    legend->Draw();
//    c1->SaveAs("Fig54A.pdf");
//
//    c1 = new TCanvas();
//    legend = new TLegend(0.15, 0.5, 0.35, 0.7);
//    count = 0;
//    for (auto& i: {1,5,9})
//    {
//        Nodes[i].GetResultPtr<TH1D>()->SetLineColor(colors[count++]);
//        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "l");
//
//        if (i == 1)
//        {
//            factor = Nodes[i].GetResultPtr<TH1D>()->Integral();
//            Nodes[i].GetResultPtr<TH1D>()->Scale(factor/Nodes[i].GetResultPtr<TH1D>()->Integral());
//            Nodes[i].GetResultPtr<TH1D>()->SetTitle(";#DeltaR (ll,#gamma) ;Events");
//            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);
//            Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->SetTitleOffset(1.2);
//            Nodes[i].GetResultPtr<TH1D>()->SetAxisRange(0., 350,"Y");
//            Nodes[i].GetResultPtr<TH1D>()->Draw("HIST");
//        }
//        else
//        {
//            Nodes[i].GetResultPtr<TH1D>()->Scale(factor/Nodes[i].GetResultPtr<TH1D>()->Integral());
//            Nodes[i].GetResultPtr<TH1D>()->SetTitle(";#DeltaR (ll,#gamma) ;Events");
//            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);
//            Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->SetTitleOffset(1.2);
//            Nodes[i].GetResultPtr<TH1D>()->SetAxisRange(0., 350,"Y");
//            Nodes[i].GetResultPtr<TH1D>()->Draw("HISTsame");
//            gPad->Modified(); gPad->Update();
//        }
//    }
//    gStyle->SetOptStat(0);
//    Tl.SetTextSize(0.03);
//    Tl.DrawLatexNDC(0.15, 0.85, "#it{ATLAS} Internal");
//    Tl.DrawLatexNDC(0.15, 0.75,"#sqrt{s} = 13 TeV");
//    legend->SetBorderSize(0);
//    legend->Draw();
//    c1->SaveAs("Fig54C.pdf");
//
//    c1 = new TCanvas();
//    legend = new TLegend(0.65, 0.4, 0.85, 0.6);
//    count = 0;
//    for (auto& i: {2,6,10})
//    {
//        Nodes[i].GetResultPtr<TH1D>()->SetLineColor(colors[count++]);
//        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "l");
//
//        if (i == 2)
//        {
//            factor = Nodes[i].GetResultPtr<TH1D>()->Integral();
//            Nodes[i].GetResultPtr<TH1D>()->Scale(factor/Nodes[i].GetResultPtr<TH1D>()->Integral());
//            Nodes[i].GetResultPtr<TH1D>()->SetTitle(";p_{T} (ll,#gamma) ;Events");
//            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);
//            Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->SetTitleOffset(1.2);
//            Nodes[i].GetResultPtr<TH1D>()->SetAxisRange(0., 560,"Y");
//            Nodes[i].GetResultPtr<TH1D>()->Draw("HIST");
//        }
//        else
//        {
//            Nodes[i].GetResultPtr<TH1D>()->Scale(factor/Nodes[i].GetResultPtr<TH1D>()->Integral());
//            Nodes[i].GetResultPtr<TH1D>()->SetTitle(";p_{T} (ll,#gamma) ;Events");
//            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);
//            Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->SetTitleOffset(1.2);
//            Nodes[i].GetResultPtr<TH1D>()->SetAxisRange(0., 560, "Y");
//            Nodes[i].GetResultPtr<TH1D>()->Draw("HISTsame");
//            gPad->Modified(); gPad->Update();
//        }
//    }
//    gStyle->SetOptStat(0);
//    Tl.SetTextSize(0.03);
//    Tl.DrawLatexNDC(0.6, 0.8, "#it{ATLAS} Internal");
//    Tl.DrawLatexNDC(0.6, 0.7,"#sqrt{s} = 13 TeV");
//    legend->SetBorderSize(0);
//    legend->Draw();
//    c1->SaveAs("Fig54D.pdf");
//
//    c1 = new TCanvas();
//    legend = new TLegend(0.45, 0.4, 0.65, 0.6);
//    count = 0;
//    for (auto& i: {3,7,11})
//    {
//
//        if (i == 3)
//        {
//            Nodes[i].GetResultPtr<TH1D>()->SetLineColor(colors[count++]);
//            legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "l");
//            factor = Nodes[i].GetResultPtr<TH1D>()->Integral();
//            Nodes[i].GetResultPtr<TH1D>()->Scale(factor/Nodes[i].GetResultPtr<TH1D>()->Integral());
//            Nodes[i].GetResultPtr<TH1D>()->SetTitle(";E_{ratio};Events");
//            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);
//            Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->SetTitleOffset(1.2);
//            Nodes[i].GetResultPtr<TH1D>()->Draw("HISTsame");
//        }
//        else
//        {
//            if (Nodes[i].GetResultPtr<TH1D>()->Integral() == 0)
//            {
//                continue;
//            }
//            Nodes[i].GetResultPtr<TH1D>()->SetLineColor(colors[count++]);
//            legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "l");
//            Nodes[i].GetResultPtr<TH1D>()->Scale(factor/Nodes[i].GetResultPtr<TH1D>()->Integral());
//            Nodes[i].GetResultPtr<TH1D>()->SetTitle(";E_{ratio};Events");
//            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);
//            Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->SetTitleOffset(1.2);
//            Nodes[i].GetResultPtr<TH1D>()->Draw("HISTsame");
//            gPad->Modified(); gPad->Update();
//        }
//    }
//    gStyle->SetOptStat(0);
//    Tl.SetTextSize(0.03);
//    Tl.DrawLatexNDC(0.4, 0.8, "#it{ATLAS} Internal");
//    Tl.DrawLatexNDC(0.4, 0.7,"#sqrt{s} = 13 TeV");
//    legend->SetBorderSize(0);
//    legend->Draw();
//    c1->SaveAs("Fig54B.pdf");
//
//}
//
//void fig29()
//{
//    auto findParentInChain = [](int targetBarcode, RVec<TruthParticle>& startParticles, RVec<TruthParticle>& truthChain)
//    {
//        RVec<TruthParticle> truthSelected;
//        bool foundParent;
//        if (truthChain.size() >= 1)
//        {
//            TruthParticle tp;
//            for (auto& tpe: startParticles)
//            {
//                tp = tpe;
//                while (true)
//                {
//                    if (tp.mc_parent_barcode == targetBarcode)
//                    {
//                        truthSelected.push_back(tp);
//                        break;
//                    }
//                    else
//                    {
//                        foundParent = false;
//                        for (auto& tmp: truthChain)
//                        {
//                            if (tp.mc_parent_barcode == tmp.mc_barcode)
//                            {
//                                tp = tmp;
//                                foundParent = true;
//                                break;
//                            }
//                        }
//                        if (foundParent == false)
//                        {
//                            break;
//                        }
//                    }
//                }
//            }
//        }
//        return truthSelected;
//    };
//
//    std::vector<std::string> input_filenames = {
////        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root/"
//        //"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Test_oneMassPoint.root",
////        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_FewMassPoints.root"
////        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_FewMassPoints_PhotonSFs.root"
////        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_MultiMass_Cayy0p001.root"
//        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_MultiMass_Cayy0p0001.root"
////        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_FewMassPoints_PhotonSFs.root"
////        "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_MultiMass_Cayy0p01.root",
////        {
////            "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600751.PhPy8EG_AZNLO_ggH125_mA2p0_v1.root",
////            "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600752.PhPy8EG_AZNLO_ggH125_mA3p0_v1.root",
////            "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600756.PhPy8EG_AZNLO_ggH125_mA9p0_v1.root"
////        }
//    };
//
//    TCanvas* c1 = new TCanvas();
//    TLegend* legend = new TLegend(0.6, 0.2, 0.85, 0.525);
//
//    SchottDataFrame df(MakeRDF(input_filenames, 8));
//
////    df.Describe().Print();
//
//    std::vector<float> massPoints = {0.2,0.5,1,3,5,10,10.1,20,20.1,29.5};
//    std::vector<ROOT::RDF::RResultHandle> Nodes;
//
//    std::cout << '\n';
//
//    auto preselection = df.Filter(
//    [&](const RVec<std::string>& trigger_passed_triggers, RVec<TruthParticle> truth_particles)
//    {
//        bool trigger_found = (std::find_first_of(trigger_passed_triggers.begin(), trigger_passed_triggers.end(), triggers.begin(), triggers.end()) != trigger_passed_triggers.end());
//        if (!trigger_found)
//        {
//            return false;
//        }
//
//        truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
//        [](TruthParticle& x)
//        {
//            return (std::abs(x.mc_pdg_id) != 11 && std::abs(x.mc_pdg_id) != 12 && std::abs(x.mc_pdg_id) != 13 &&
//                    std::abs(x.mc_pdg_id) != 14 && std::abs(x.mc_pdg_id) != 15 && std::abs(x.mc_pdg_id) != 16 &&
//                    std::abs(x.mc_pdg_id) != 17 && std::abs(x.mc_pdg_id) != 18);
//
//        }), truth_particles.end());
//
//        if (truth_particles.size() != 2)
//        {
//            return false;
//        }
//
//        if (DeltaR(truth_particles[0].Vector(), truth_particles[1].Vector()) <= 0.01)
//        {
//            return false;
//        }
//
//        if (truth_particles[0].mc_charge*truth_particles[1].mc_charge >= 0)
//        {
//            return false;
//        }
//
//        if (!((truth_particles[0].mc_pt > 27e3 && truth_particles[1].mc_pt > 20e3)
//                                               ||
//              (truth_particles[1].mc_pt > 27e3 && truth_particles[0].mc_pt > 20e3)))
//        {
//            return false;
//        }
//
//        PtEtaPhiEVector dilepton = truth_particles[0].Vector() + truth_particles[1].Vector();
//
//        if ((dilepton.M() < 81e3) || (dilepton.M() > 101e3))
//        {
//            return false;
//        }
//
//        if ((truth_particles[0].Vector() + truth_particles[1].Vector()).Pt() <= 10e3)
//        {
//            return false;
//        }
//
//        return true;
//
//    }, {"trigger_passed_triggers", "truth_particles"});
//
//    auto truth_photons_from_axions = df.Define("truth_photons",[&](RVec<TruthParticle> truth_particles)
//    {
//        truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
//        [](TruthParticle& x)
//        {
//            return (std::abs(x.mc_pdg_id) != 22);
//
//        }), truth_particles.end());
//
//        return truth_particles;
//
//    }, {"truth_particles"})
//    .Define("truth_axions", [&](RVec<TruthParticle> truth_particles)
//    {
//        truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
//        [](TruthParticle& x)
//        {
//            return (std::abs(x.mc_pdg_id) != 36);
//
//        }), truth_particles.end());
//
//        return truth_particles;
//
//    }, {"truth_particles"}).Define("truth_photons_from_axions",
//    [&](RVec<TruthParticle>& truth_photons, RVec<TruthParticle> truth_axions)
//    {
//        std::sort(truth_axions.begin(),truth_axions.end(),
//        [](TruthParticle& x, TruthParticle& y)
//        {
//            return x.mc_barcode < y.mc_barcode;
//        });
//
//        return findParentInChain(truth_axions[0].mc_barcode, truth_photons, truth_axions);
//
//    }, {"truth_photons", "truth_axions"})
//    .Define("truth_photons_from_axions_size",
//    [](RVec<TruthParticle>& truth_photons_from_axions)
//    {
//        return truth_photons_from_axions.size();
//    }, {"truth_photons_from_axions"})
////    .Filter([&] (RVec<TruthParticle>& truth_photons_from_axions)
////    {
////        return truth_photons_from_axions.size();// == 2;
//////        return true;
////    }, {"truth_photons_from_axions"})
//    .Define("diphoton_inv_mass", [&](RVec<TruthParticle> truth_photons_from_axions)
//    {
//        auto four_momentum = truth_photons_from_axions[0].Vector() + truth_photons_from_axions[1].Vector();
//        return four_momentum.M()/1e3f;
//    }, {"truth_photons"})
//    .Define("axion_masses", [&](RVec<TruthParticle>& truth_axions)
//    {
//        return truth_axions[0].mc_mass/1e3f;
//
//    }, {"truth_axions"})
//    .Define("axion_decay_times", [&](RVec<TruthParticle>& truth_axions)
//    {
//        return truth_axions[0].mc_decay_time;
//
//    }, {"truth_axions"});
//
//    auto numAxions = truth_photons_from_axions.Filter(
//    [](RVec<TruthParticle>& truth_particles)
//    {
//        return !truth_particles.empty();
//    }, {"truth_axions"});
//
//    for (auto& mass_point: massPoints)
//    {
//        auto mass_point_truth_photons_from_axions = truth_photons_from_axions.Filter([&]
//        (float axion_mass)
//        {
//            return (roundToOneDecimalPlace(axion_mass) == mass_point);
//
//        }, {"axion_masses"});
//
//        Nodes.push_back(mass_point_truth_photons_from_axions.Histo1D<std::vector<float>>({"MC decay lengths", "MC decay lengths", 60u, 0, 15000}, "mc_decay_length"));
//    }
//
//    ROOT::RDF::RunGraphs(Nodes); // running all computation nodes concurrently
//
//    auto vals = truth_photons_from_axions.Take<float, RVec<float>>("axion_masses");
//
//    auto number_of_truth_photons_in_each_event = truth_photons_from_axions.Take<size_t, RVec<size_t>>("truth_photons_from_axions_size");
//
//    ROOT::RDF::RResultPtr<ULong64_t> Count = truth_photons_from_axions.Count()
//    , AxionCount = numAxions.Count();
//
//    ROOT::RDF::RResultPtr<TH1D> histo = truth_photons_from_axions.Histo1D<double>({"Diphoton Invariant Mass", "Diphoton Invariant Mass", 60u, 0, 40}, "diphoton_inv_mass");
//
//    ROOT::RDF::RResultPtr<TH1D> mc_decay_time_histo = truth_photons_from_axions.Histo1D<float>({"MC decay times", "MC decay times", 60u, 0, 0}, "axion_decay_times");
//
//    ROOT::RDF::RResultPtr<TH1D> histo_axion_mass = truth_photons_from_axions.Histo1D<float>({"Axion Mass", "Axion Mass", 60u, 0, 0}, "axion_masses");
//
//    std::cout << *Count << '\n' << *AxionCount << '\n';
//
//    std::map<int, int> truth_photons_from_axion_counts;
//
//    for (auto& i: *number_of_truth_photons_in_each_event)
//    {
//        truth_photons_from_axion_counts[i]++;
//    }
//
//    for (auto& i: truth_photons_from_axion_counts)
//    {
//        std::cout <<
//        "number of events with " <<
//        i.first
//        << " truth photons from the axion = "
//        << i.second << '\n';
//    }
//
//    histo->SetLineColor(kGreen);
//    legend->AddEntry(&(*histo), histo->GetTitle(), "l");
//    histo->SetTitleOffset(1.2);
//    histo->GetYaxis()->CenterTitle(true);
//    histo->SetTitle(";m_{#gamma#gamma}  [GeV]; Events");
//    histo->Draw("HIST");
//
//    std::map<float, int> unique_vals;
//
//    for (auto& i: *vals)
//    {
////        unique_vals.emplace(roundToOneDecimalPlace(i));
//        unique_vals[roundToOneDecimalPlace(i)]++;
//    }
//    std::cout << (*vals).size() << '\n';
//    std::ofstream out("massPoints2.txt");
//    for (auto& i: unique_vals)
//    {
//        out << i.second << R"--( & \textcolor{green}{\checkmark} & \textcolor{green}{\checkmark} & \textcolor{green}{\checkmark} \\ \hline %)--" << i.first << '\n';
//    }
//    out.close();
//
////    gStyle->SetOptStat(0);
//    TLatex Tl;
//    Tl.SetTextSize(0.03);
//    Tl.DrawLatexNDC(0.6, 0.83, "#it{ATLAS} Internal");
//    Tl.DrawLatexNDC(0.6, 0.77,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
//    Tl.DrawLatexNDC(0.6, 0.7,"ggF m_{A} = {...} GeV");
//    legend->SetBorderSize(0);
//    legend->Draw();
//    c1->SaveAs("Fig29.pdf");
//
//    c1 = new TCanvas();
//    legend = new TLegend(0.5, 0.8, 0.25, 0.6);
//
//    mc_decay_time_histo->SetLineColor(kOrange);
//    legend->AddEntry(&(*mc_decay_time_histo), mc_decay_time_histo->GetTitle(), "l");
//    mc_decay_time_histo->SetTitle(";Decay times (wonder what the units are?);Events");
//    mc_decay_time_histo->SetTitleOffset(1.2);
//    mc_decay_time_histo->GetYaxis()->CenterTitle(true);
//    mc_decay_time_histo->SetAxisRange(0., 35e3, "X");
//
//    mc_decay_time_histo->Draw("HIST");
//
//    gStyle->SetOptStat(0);
//    Tl.SetTextSize(0.03);
//    Tl.DrawLatexNDC(0.6, 0.83, "#it{ATLAS} Internal");
//    Tl.DrawLatexNDC(0.6, 0.77,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
//    Tl.DrawLatexNDC(0.6, 0.7,"ggF m_{A} = {...} GeV");
//    legend->SetBorderSize(0);
//    legend->Draw();
//    c1->SaveAs("Decay_Times.pdf");
//
//    c1 = new TCanvas();
//    legend = new TLegend(0.2, 0.5, 0.45, 0.8);
//
//    histo_axion_mass->SetLineColor(kRed);
//    legend->AddEntry(&(*histo_axion_mass), histo_axion_mass->GetTitle(), "l");
//    histo_axion_mass->SetTitleOffset(1.2);
//    histo_axion_mass->GetYaxis()->CenterTitle(true);
//    histo_axion_mass->SetTitle(";m_{a}  [GeV]; Events");
//    histo_axion_mass->Draw("HIST");
//
////    gStyle->SetOptStat(0);
//    Tl.SetTextSize(0.03);
//    Tl.DrawLatexNDC(0.6, 0.83, "#it{ATLAS} Internal");
//    Tl.DrawLatexNDC(0.6, 0.77,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
//    Tl.DrawLatexNDC(0.6, 0.7,"ggF m_{A} = {...} GeV");
//    legend->SetBorderSize(0);
//    legend->Draw();
//    c1->SaveAs("Axion_masses.pdf");
//
//    for (int i = 0; i < Nodes.size(); i++)
//    {
//        if (Nodes[i].GetResultPtr<TH1D>()->GetEntries() == 0.0)
//        {
//            std::cout << "works\n";
//            continue;
//        }
//
//        c1 = new TCanvas();
//        legend = new TLegend(0.2, 0.5, 0.45, 0.8);
//
//        std::string title_string = std::to_string(massPoints[i]) + " GeV;decay length (Wonder what the units are);Events";
//
//        Nodes[i].GetResultPtr<TH1D>()->SetTitle(title_string.c_str());
//        Nodes[i].GetResultPtr<TH1D>()->SetTitleOffset(1.2);
//
//        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "l");
//
//        Nodes[i].GetResultPtr<TH1D>()->Draw("HIST");
//
//        gStyle->SetOptStat(0);
//        Tl.SetTextSize(0.03);
//        Tl.DrawLatexNDC(0.6, 0.83, "#it{ATLAS} Internal");
//        Tl.DrawLatexNDC(0.6, 0.77,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
//
//        std::string MassString = std::string("ggF m_{A} = {") + std::to_string(massPoints[i]) + "} GeV";
//        Tl.DrawLatexNDC(0.6, 0.7,MassString.c_str());
//        legend->SetBorderSize(0);
//        legend->Draw();
//
//        c1->SaveAs((title_string+".pdf").c_str());
//    }
//
//}

//void AxionMassEventCounter()
//{
//    std::vector<std::vector<std::string>> input_filenames =
//    {
//        {
//            "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_FewMassPoints_PhotonSFs.root"
//        },      //C_{ayy} = 1
//        {
//            "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_MultiMass_Cayy0p01.root",
//            "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_mA1p0_Cayy0p01.root"
//        },      //C_{ayy} = 0.01
//        {
//            "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_MultiMass_Cayy0p001.root",
//            "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_mA1p0_Cayy0p001.root"
//        },      //C_{ayy} = 0.001
//        {
//            "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_MultiMass_Cayy0p0001.root",
//                //C_{ayy} = 0.0001
//        }
//    };
//
//    std::vector<ROOT::RDF::RResultHandle> Nodes;
//
//    for (auto& sample: input_filenames)
//    {
//        SchottDataFrame df(MakeRDF(sample, 8));
//
//        auto truth_photons_from_axions = df.Define("truth_photons",
//        [&](RVec<TruthParticle> truth_particles)
//        {
//            truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
//            [](TruthParticle& x)
//            {
//                return (std::abs(x.mc_pdg_id) != 22);
//
//            }), truth_particles.end());
//
//            return truth_particles;
//
//        }, {"truth_particles"})
//        .Define("truth_axions",
//        [&](RVec<TruthParticle> truth_particles)
//        {
//            truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
//            [](TruthParticle& x)
//            {
//                return (std::abs(x.mc_pdg_id) != 36);
//
//            }), truth_particles.end());
//
//            return truth_particles;
//
//        }, {"truth_particles"})
//        .Define("axion_masses", [&](RVec<TruthParticle>& truth_axions)
//        {
//            return truth_axions[0].mc_mass/1e3f;
//
//        }, {"truth_axions"});
//
//        Nodes.push_back(truth_photons_from_axions.Take<float, RVec<float>>("axion_masses"));
//    }
//
//    ROOT::RDF::RunGraphs(Nodes); // running all computation nodes concurrently
//
//    std::map<float, int> unique_vals_1, unique_vals_0p01, unique_vals_0p001, unique_vals_0p0001;
//    std::set<float> allVals;
//
//    for (auto& i: *Nodes[0].GetResultPtr<RVec<float>>()) //C_ayy = 1
//    {
//        unique_vals_1[roundToOneDecimalPlace(i)]++;
//        allVals.emplace(roundToOneDecimalPlace(i));
//    }
//
//    for (auto& i: *Nodes[1].GetResultPtr<RVec<float>>()) //C_ayy = 0.01
//    {
//        unique_vals_0p01[roundToOneDecimalPlace(i)]++;
//        allVals.emplace(roundToOneDecimalPlace(i));
//    }
//
//    for (auto& i: *Nodes[2].GetResultPtr<RVec<float>>()) //C_ayy = 0.001
//    {
//        unique_vals_0p001[roundToOneDecimalPlace(i)]++;
//        allVals.emplace(roundToOneDecimalPlace(i));
//    }
//
//    for (auto& i: *Nodes[3].GetResultPtr<RVec<float>>()) //C_ayy = 0.0001
//    {
//        unique_vals_0p0001[roundToOneDecimalPlace(i)]++;
//        allVals.emplace(roundToOneDecimalPlace(i));
//    }
//
//    std::cout << allVals.size() << '\n';
//    std::ofstream out("massPoints2.txt");
//    std::vector<float> massPoints = {0.2, 0.5, 1, 2, 3, 4, 5, 7, 10, 15, 20, 25, 30};
//    for (auto& val: massPoints)
//    {
//        out << val << " & "
//        << unique_vals_1[val] << " & "
//        << unique_vals_0p01[val] << " & "
//        << unique_vals_0p001[val] << " & "
//        << unique_vals_0p0001[val] << R"--( \\ \hline %)--"
//        << val << '\n';
//    }
//    out.close();
//
//}
//
//void SignalShapes()
//{
//    auto findParentInChain = [](int targetBarcode, RVec<TruthParticle>& startParticles, RVec<TruthParticle>& truthChain)
//    {
//        RVec<TruthParticle> truthSelected;
//        bool foundParent;
//        if (truthChain.size() >= 1)
//        {
//            TruthParticle tp;
//            for (auto& tpe: startParticles)
//            {
//                tp = tpe;
//                while (true)
//                {
//                    if (tp.mc_parent_barcode == targetBarcode)
//                    {
//                        truthSelected.push_back(tp);
//                        break;
//                    }
//                    else
//                    {
//                        foundParent = false;
//                        for (auto& tmp: truthChain)
//                        {
//                            if (tp.mc_parent_barcode == tmp.mc_barcode)
//                            {
//                                tp = tmp;
//                                foundParent = true;
//                                break;
//                            }
//                        }
//                        if (foundParent == false)
//                        {
//                            break;
//                        }
//                    }
//                }
//            }
//        }
//        return truthSelected;
//    };
//
//    std::vector<std::vector<std::string>> input_filenames =
//    {
//        {
//            "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root"
//        },      // 5 GeV
//        {
//            "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600752.PhPy8EG_AZNLO_ggH125_mA3p0_v1.root"
//        },       // 3 GeV
//        {
//            "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root"
//        },      // 1 GeV
//        {
//            "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_FewMassPoints_PhotonSFs.root"
//        },      //C_{ayy} = 1
//        {
//            "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_MultiMass_Cayy0p01.root",
//            "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_mA1p0_Cayy0p01.root"
//        },      //C_{ayy} = 0.01
//        {
//            "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_MultiMass_Cayy0p001.root",
//            "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_mA1p0_Cayy0p001.root"
//        },     //C_{ayy} = 0.001
//        {
//            "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_MultiMass_Cayy0p0001.root"
//        }     //C_{ayy} = 0.0001
//    };
//
//    TCanvas* c1;
//    TLegend* legend;
//    TLatex Tl;
//    std::vector<float> massPoints = {5, 3, 1};
//    std::vector<ROOT::RDF::RResultHandle> Nodes;
//    int count = 0;
//
//    for (auto& input_file: input_filenames)
//    {
//        SchottDataFrame df(MakeRDF(input_file, 8));
//
//        auto truth_photons_from_axions = df.Define("truth_photons",[&](RVec<TruthParticle> truth_particles)
//        {
//            truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
//            [](TruthParticle& x)
//            {
//                return (std::abs(x.mc_pdg_id) != 22);
//
//            }), truth_particles.end());
//
//            return truth_particles;
//
//        }, {"truth_particles"})
//        .Define("truth_axions", [&](RVec<TruthParticle> truth_particles)
//        {
//            truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
//            [&](TruthParticle& x)
//            {
//                return (count <= 2 ? (std::abs(x.mc_pdg_id) != 35) : (std::abs(x.mc_pdg_id) != 36));
//
//            }), truth_particles.end());
//
//            return truth_particles;
//
//        }, {"truth_particles"}).Define("truth_photons_from_axions",
//        [&](RVec<TruthParticle>& truth_photons, RVec<TruthParticle> truth_axions)
//        {
//            std::sort(truth_axions.begin(),truth_axions.end(),
//            [](TruthParticle& x, TruthParticle& y)
//            {
//                return x.mc_barcode < y.mc_barcode;
//            });
//
//            return findParentInChain(truth_axions[0].mc_barcode, truth_photons, truth_axions);
//
//        }, {"truth_photons", "truth_axions"})
//        .Define("diphoton_inv_mass", [&](RVec<TruthParticle> truth_photons_from_axions)
//        {
//            auto four_momentum = truth_photons_from_axions[0].Vector() + truth_photons_from_axions[1].Vector();
//            return four_momentum.M()/1e3f;
//        }, {"truth_photons"})
//        .Define("axion_masses", [&](RVec<TruthParticle>& truth_axions)
//        {
//            return truth_axions[0].mc_mass/1e3f;
//
//        }, {"truth_axions"});
//
//        if (count <= 2)
//        {
//            Nodes.push_back(truth_photons_from_axions.Histo1D<double>({"", "", 120u, 0.2, 5.8}, "diphoton_inv_mass"));
//        }
//        else
//        {
//            for (auto& mass_point: massPoints)
//            {
//                auto mass_point_truth_photons_from_axions = truth_photons_from_axions.Filter([&]
//                (float axion_mass)
//                {
//                    return (roundToOneDecimalPlace(axion_mass) == mass_point);
//
//                }, {"axion_masses"});
//
//                Nodes.push_back(mass_point_truth_photons_from_axions.Histo1D<double>({"", "", 120u, 0.2, 5.8}, "diphoton_inv_mass"));
//            }
//        }
//        count++;
//    }
//
////    0           //1 GeV prompt
////    1           //3 GeV prompt
////    2           //5 GeV prompt
////    3  4   5    //1, 3, 5 GeV for Cayy=1
////    6  7   8    //1, 3, 5 GeV for Cayy=0.01
////    9  10  11   //1, 3, 5 GeV for Cayy=0.001
////    12 13  14   //1, 3, 5 GeV for Cayy=0.0001
//
//    ROOT::RDF::RunGraphs(Nodes); // running all computation nodes concurrently
//    gStyle->SetOptStat(0);
//    c1 = new TCanvas();
//    legend = new TLegend(0.25, 0.55, 0.5, 0.8);
//    constexpr std::array<EColor,3> colors = {kBlue, static_cast<EColor>(kGreen + 2), kViolet};
//    constexpr std::array<int,3> massVals = {5, 3, 1};
//    std::string legendTitle;
//    double normFactor;
//
//    //prompt samples
//
//    for (int i = 0; i <= 2; i++)
//    {
//        Nodes[i].GetResultPtr<TH1D>()->SetLineColor(colors[i]);
//        legendTitle = std::string("m_{a} = ") + std::to_string(massVals[i]) + "  GeV (prompt)";
//        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), legendTitle.c_str(), "l");
//        if (i == 0)
//        {
//            Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
//            Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->SetTitleSize(1.25 * Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->GetTitleSize());
//            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->SetTitleOffset(1.4);
//            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->SetTitle("Events");
//            Nodes[i].GetResultPtr<TH1D>()->SetMaximum(10200);
//            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);
//            Nodes[i].GetResultPtr<TH1D>()->Draw("HIST");
//            normFactor = Nodes[i].GetResultPtr<TH1D>()->Integral();
//        }
//        else
//        {
//            Nodes[i].GetResultPtr<TH1D>()->Scale(normFactor/Nodes[i].GetResultPtr<TH1D>()->Integral());
//            Nodes[i].GetResultPtr<TH1D>()->Draw("HISTsame");
//        }
//    }
//
//    Tl.SetTextSize(0.03);
//    Tl.DrawLatexNDC(0.6, 0.83, "#it{ATLAS} Internal");
//    Tl.DrawLatexNDC(0.6, 0.73,"#sqrt{s} = 13 TeV");
//    legend->SetBorderSize(0);
//    legend->Draw();
//    c1->SaveAs("Axion_masses_prompt.pdf");
//
//    //Displaced Samples
//    constexpr std::array<const char*, 4> legendLabels = {" (C_{a#gamma#gamma} = 1)", " (C_{a#gamma#gamma} = 0.01)", " (C_{a#gamma#gamma} = 0.001)", " (C_{a#gamma#gamma} = 0.0001)"};
//    constexpr std::array<const char*, 4> fileNames = {"Axion_masses_Cayy_1.pdf", "Axion_masses_Cayy_0p01.pdf", "Axion_masses_Cayy_0p001.pdf", "Axion_masses_Cayy_0p0001.pdf"};
//    constexpr std::array<float, 4> maxima = {103, 858, 700.5, 705};
//
//    for (int i = 3, k = 0; i <= 12; i += 3, k++)
//    {
//        c1 = new TCanvas();
//        legend = new TLegend(0.25, 0.55, 0.5, 0.85);
//        for (int j = 0; j <= 2; j++)
//        {
//            Nodes[i+j].GetResultPtr<TH1D>()->SetLineColor(colors[j]);
//            legendTitle = std::string("m_{a} = ") + std::to_string(massVals[j]) + "  GeV" + legendLabels[k];
//            legend->AddEntry(&(*Nodes[i+j].GetResultPtr<TH1D>()), legendTitle.c_str(), "l");
//            if (j == 0) //1 Gev
//            {
//                normFactor = Nodes[i+j].GetResultPtr<TH1D>()->Integral();
//                Nodes[i+j].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);
//                Nodes[i+j].GetResultPtr<TH1D>()->GetYaxis()->SetTitleOffset(1.2);
//                Nodes[i+j].GetResultPtr<TH1D>()->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
//                Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->SetTitleSize(1.25 * Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->GetTitleSize());
//                Nodes[i+j].GetResultPtr<TH1D>()->SetMaximum(maxima[k]);
//                Nodes[i+j].GetResultPtr<TH1D>()->GetYaxis()->SetTitle("Events");
//                Nodes[i+j].GetResultPtr<TH1D>()->Draw("HIST");
//            }
//            else
//            {
//                Nodes[i+j].GetResultPtr<TH1D>()->Scale(normFactor/Nodes[i+j].GetResultPtr<TH1D>()->Integral());
//                Nodes[i+j].GetResultPtr<TH1D>()->Draw("HISTsame");
//            }
//        }
//        Tl.SetTextSize(0.03);
//        Tl.DrawLatexNDC(0.6, 0.83, "#it{ATLAS} Internal");
//        Tl.DrawLatexNDC(0.6, 0.73,"#sqrt{s} = 13 TeV");
//        legend->SetBorderSize(0);
//        legend->Draw();
//        c1->SaveAs(fileNames[k]);
//    }
//
//}
//
//void SignalShapesDeltaR()
//{
//    auto findParentInChain = [](int targetBarcode, RVec<TruthParticle>& startParticles, RVec<TruthParticle>& truthChain)
//    {
//        RVec<TruthParticle> truthSelected;
//        bool foundParent;
//        if (truthChain.size() >= 1)
//        {
//            TruthParticle tp;
//            for (auto& tpe: startParticles)
//            {
//                tp = tpe;
//                while (true)
//                {
//                    if (tp.mc_parent_barcode == targetBarcode)
//                    {
//                        truthSelected.push_back(tp);
//                        break;
//                    }
//                    else
//                    {
//                        foundParent = false;
//                        for (auto& tmp: truthChain)
//                        {
//                            if (tp.mc_parent_barcode == tmp.mc_barcode)
//                            {
//                                tp = tmp;
//                                foundParent = true;
//                                break;
//                            }
//                        }
//                        if (foundParent == false)
//                        {
//                            break;
//                        }
//                    }
//                }
//            }
//        }
//        return truthSelected;
//    };
//
//    std::vector<std::vector<std::string>> input_filenames =
//    {
//        {
//            "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root"
//        },      // 5 GeV
//        {
//            "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600752.PhPy8EG_AZNLO_ggH125_mA3p0_v1.root"
//        },       // 3 GeV
//        {
//            "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root"
//        },      // 1 GeV
//        {
//            "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_FewMassPoints_PhotonSFs.root"
//        },      //C_{ayy} = 1
//        {
//            "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_MultiMass_Cayy0p01.root",
//            "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_mA1p0_Cayy0p01.root"
//        },      //C_{ayy} = 0.01
//        {
//            "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_MultiMass_Cayy0p001.root",
//            "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_mA1p0_Cayy0p001.root"
//        },     //C_{ayy} = 0.001
//        {
//            "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_MultiMass_Cayy0p0001.root"
//        }     //C_{ayy} = 0.0001
//    };
//
//    TCanvas* c1;
//    TLegend* legend;
//    TLatex Tl;
//    std::vector<float> massPoints = {5, 3, 1};
//    std::vector<ROOT::RDF::RResultHandle> Nodes;
//    int count = 0;
//
//    for (auto& input_file: input_filenames)
//    {
//        SchottDataFrame df(MakeRDF(input_file, 8));
//
//        auto truth_photons_from_axions = df.Define("truth_photons",[&](RVec<TruthParticle> truth_particles)
//        {
//            truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
//            [](TruthParticle& x)
//            {
//                return (std::abs(x.mc_pdg_id) != 22);
//
//            }), truth_particles.end());
//
//            return truth_particles;
//
//        }, {"truth_particles"})
//        .Define("truth_axions", [&](RVec<TruthParticle> truth_particles)
//        {
//            truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
//            [&](TruthParticle& x)
//            {
//                return (count <= 2 ? (std::abs(x.mc_pdg_id) != 35) : (std::abs(x.mc_pdg_id) != 36));
//
//            }), truth_particles.end());
//
//            return truth_particles;
//
//        }, {"truth_particles"}).Define("truth_photons_from_axions",
//        [&](RVec<TruthParticle>& truth_photons, RVec<TruthParticle> truth_axions)
//        {
//            std::sort(truth_axions.begin(),truth_axions.end(),
//            [](TruthParticle& x, TruthParticle& y)
//            {
//                return x.mc_barcode < y.mc_barcode;
//            });
//
//            return findParentInChain(truth_axions[0].mc_barcode, truth_photons, truth_axions);
//
//        }, {"truth_photons", "truth_axions"})
//        .Define("DeltaR_truth_photons_from_axion", [&](RVec<TruthParticle> truth_photons_from_axions)
//        {
//            return DeltaR(truth_photons_from_axions[0].Vector(), truth_photons_from_axions[1].Vector());
//
//        }, {"truth_photons"})
//        .Define("axion_masses", [&](RVec<TruthParticle>& truth_axions)
//        {
//            return truth_axions[0].mc_mass/1e3f;
//
//        }, {"truth_axions"});
//
//        if (count <= 2)
//        {
//            Nodes.push_back(truth_photons_from_axions.Histo1D<double>({"", "", 120u, 0, 1.2}, "DeltaR_truth_photons_from_axion"));
//        }
//        else
//        {
//            for (auto& mass_point: massPoints)
//            {
//                auto mass_point_truth_photons_from_axions = truth_photons_from_axions.Filter([&]
//                (float axion_mass)
//                {
//                    return (roundToOneDecimalPlace(axion_mass) == mass_point);
//
//                }, {"axion_masses"});
//
//                Nodes.push_back(mass_point_truth_photons_from_axions.Histo1D<double>({"", "", 120u, 0, 1.2}, "DeltaR_truth_photons_from_axion"));
//            }
//        }
//        count++;
//    }
//
////    0           //1 GeV prompt
////    1           //3 GeV prompt
////    2           //5 GeV prompt
////    3  4   5    //1, 3, 5 GeV for Cayy=1
////    6  7   8    //1, 3, 5 GeV for Cayy=0.01
////    9  10  11   //1, 3, 5 GeV for Cayy=0.001
////    12 13  14   //1, 3, 5 GeV for Cayy=0.0001
//
//    ROOT::RDF::RunGraphs(Nodes); // running all computation nodes concurrently
//    gStyle->SetOptStat(0);
//    c1 = new TCanvas();
//    legend = new TLegend(0.25, 0.55, 0.5, 0.8);
//    constexpr std::array<EColor,3> colors = {kBlue, static_cast<EColor>(kGreen + 2), kViolet};
//    constexpr std::array<int,3> massVals = {5, 3, 1};
//    std::string legendTitle;
//    double normFactor;
//
//    //prompt samples
//
//    for (int i = 0; i <= 2; i++)
//    {
//        Nodes[i].GetResultPtr<TH1D>()->SetLineColor(colors[i]);
//        legendTitle = std::string("m_{a} = ") + std::to_string(massVals[i]) + "  GeV (prompt)";
//        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), legendTitle.c_str(), "l");
//        if (i == 0)
//        {
//            Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->SetTitle("#Delta R_{#gamma#gamma}");
//            Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->SetTitleSize(1.25 * Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->GetTitleSize());
//            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->SetTitleOffset(1.4);
//            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->SetTitle("Events");
//            Nodes[i].GetResultPtr<TH1D>()->SetMaximum(770);
//            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);
//            Nodes[i].GetResultPtr<TH1D>()->Draw("HIST");
//            normFactor = Nodes[i].GetResultPtr<TH1D>()->Integral();
//        }
//        else
//        {
//            Nodes[i].GetResultPtr<TH1D>()->Scale(normFactor/Nodes[i].GetResultPtr<TH1D>()->Integral());
//            Nodes[i].GetResultPtr<TH1D>()->Draw("HISTsame");
//        }
//    }
//
//    Tl.SetTextSize(0.03);
//    Tl.DrawLatexNDC(0.6, 0.83, "#it{ATLAS} Internal");
//    Tl.DrawLatexNDC(0.6, 0.73,"#sqrt{s} = 13 TeV");
//    legend->SetBorderSize(0);
//    legend->Draw();
//    c1->SaveAs("DeltaR_prompt.pdf");
//
//    //Displaced Samples
//    constexpr std::array<const char*, 4> legendLabels = {" (C_{a#gamma#gamma} = 1)", " (C_{a#gamma#gamma} = 0.01)", " (C_{a#gamma#gamma} = 0.001)", " (C_{a#gamma#gamma} = 0.0001)"};
//    constexpr std::array<const char*, 4> fileNames = {"DeltaR_Cayy_1.pdf", "DeltaR_Cayy_0p01.pdf", "DeltaR_Cayy_0p001.pdf", "DeltaR_Cayy_0p0001.pdf"};
//    constexpr std::array<float, 4> maxima = {10.5, 66, 51, 198};
//
//    for (int i = 3, k = 0; i <= 12; i += 3, k++)
//    {
//        c1 = new TCanvas();
//        if (i == 12)
//        {
//            legend = new TLegend(0.35, 0.35, 0.6, 0.65);
//        }
//        else
//        {
//            legend = new TLegend(0.25, 0.55, 0.5, 0.85);
//        }
//
//        for (int j = 0; j <= 2; j++)
//        {
//            Nodes[i+j].GetResultPtr<TH1D>()->SetLineColor(colors[j]);
//            legendTitle = std::string("m_{a} = ") + std::to_string(massVals[j]) + "  GeV" + legendLabels[k];
//            legend->AddEntry(&(*Nodes[i+j].GetResultPtr<TH1D>()), legendTitle.c_str(), "l");
//            if (j == 0) //1 Gev
//            {
//                normFactor = Nodes[i+j].GetResultPtr<TH1D>()->Integral();
//                Nodes[i+j].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);
//                Nodes[i+j].GetResultPtr<TH1D>()->SetMaximum(maxima[k]);
//                Nodes[i+j].GetResultPtr<TH1D>()->GetYaxis()->SetTitleOffset(1.2);
//                Nodes[i+j].GetResultPtr<TH1D>()->GetXaxis()->SetTitle("#Delta R_{#gamma#gamma}");
//                Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->SetTitleSize(1.25 * Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->GetTitleSize());
//                Nodes[i+j].GetResultPtr<TH1D>()->GetYaxis()->SetTitle("Events");
//                Nodes[i+j].GetResultPtr<TH1D>()->Draw("HIST");
//            }
//            else
//            {
//                Nodes[i+j].GetResultPtr<TH1D>()->Scale(normFactor/Nodes[i+j].GetResultPtr<TH1D>()->Integral());
//                Nodes[i+j].GetResultPtr<TH1D>()->Draw("HISTsame");
//            }
//        }
//        Tl.SetTextSize(0.03);
//        Tl.DrawLatexNDC(0.6, 0.83, "#it{ATLAS} Internal");
//        Tl.DrawLatexNDC(0.6, 0.73,"#sqrt{s} = 13 TeV");
//        legend->SetBorderSize(0);
//        legend->Draw();
//        c1->SaveAs(fileNames[k]);
//    }
//
//}

//void ALP_Decay_Lengths()
//{
//    auto findParentInChain = [](int targetBarcode, RVec<TruthParticle>& startParticles, RVec<TruthParticle>& truthChain)
//    {
//        RVec<TruthParticle> truthSelected;
//        bool foundParent;
//        if (truthChain.size() >= 1)
//        {
//            TruthParticle tp;
//            for (auto& tpe: startParticles)
//            {
//                tp = tpe;
//                while (true)
//                {
//                    if (tp.mc_parent_barcode == targetBarcode)
//                    {
//                        truthSelected.push_back(tp);
//                        break;
//                    }
//                    else
//                    {
//                        foundParent = false;
//                        for (auto& tmp: truthChain)
//                        {
//                            if (tp.mc_parent_barcode == tmp.mc_barcode)
//                            {
//                                tp = tmp;
//                                foundParent = true;
//                                break;
//                            }
//                        }
//                        if (foundParent == false)
//                        {
//                            break;
//                        }
//                    }
//                }
//            }
//        }
//        return truthSelected;
//    };
//
//    std::vector<std::vector<std::string>> input_filenames =
//    {
//        {
//            "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_FewMassPoints_PhotonSFs.root"
//        },      //C_{ayy} = 1
//        {
//            "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_MultiMass_Cayy0p01.root",
//            "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_mA1p0_Cayy0p01.root"
//        },      //C_{ayy} = 0.01
//        {
//            "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_MultiMass_Cayy0p001.root",
//            "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_mA1p0_Cayy0p001.root"
//        },     //C_{ayy} = 0.001
//        {
//            "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_MultiMass_Cayy0p0001.root"
//            //C_{ayy} = 0.0001
//        }
//    };
//
//    std::vector<ROOT::RDF::RResultHandle> Nodes;
//    int count = 0;
//
//    for (auto& input_file: input_filenames)
//    {
//        SchottDataFrame df(MakeRDF(input_file, 8));
//
//        auto truth_photons_from_axions = df.Define("truth_photons",[&](RVec<TruthParticle> truth_particles)
//        {
//            truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
//            [](TruthParticle& x)
//            {
//                return (std::abs(x.mc_pdg_id) != 22);
//
//            }), truth_particles.end());
//
//            return truth_particles;
//
//        }, {"truth_particles"})
//        .Define("truth_axions", [&](RVec<TruthParticle> truth_particles)
//        {
//            truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
//            [&](TruthParticle& x)
//            {
//                return (std::abs(x.mc_pdg_id) != 36);
//
//            }), truth_particles.end());
//
//            return truth_particles;
//
//        }, {"truth_particles"}).Define("truth_photons_from_axions",
//        [&](RVec<TruthParticle>& truth_photons, RVec<TruthParticle> truth_axions)
//        {
//            std::sort(truth_axions.begin(),truth_axions.end(),
//            [](TruthParticle& x, TruthParticle& y)
//            {
//                return x.mc_barcode < y.mc_barcode;
//            });
//
//            return findParentInChain(truth_axions[0].mc_barcode, truth_photons, truth_axions);
//
//        }, {"truth_photons", "truth_axions"})
//        .Define("truth_higgs", [&](RVec<TruthParticle> truth_particles)
//        {
//            truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
//            [&](TruthParticle& x)
//            {
//                return (std::abs(x.mc_pdg_id) != 35);
//
//            }), truth_particles.end());
//
//            return truth_particles;
//
//        }, {"truth_particles"})
//        .Define("ALP_decay_length",
//        [&](RVec<TruthParticle>& truth_higgs, RVec<TruthParticle>& truth_axions, ROOT::VecOps::RVec<std::vector<int>>& mc_vertex_incoming_barcodes, RVec<float>& mc_vertex_x, RVec<float>& mc_vertex_y)
//        {
//            int startX = 0, startY = 0;
//            int truth_higgs_barcode = truth_higgs[truth_higgs.size() - 1].mc_barcode;
//
//            for (int i = 0; i < mc_vertex_incoming_barcodes.size(); i++)
//            {
//                if (mc_vertex_incoming_barcodes[i].size() > 0 && mc_vertex_incoming_barcodes[i][0] == truth_higgs_barcode)
//                {
//                    startX = mc_vertex_x[i];
//                    startY = mc_vertex_y[i];
//                    break;
//                }
//            }
//
//            double decayLength_a1 = 0;
//            for (int i = 0; i < mc_vertex_incoming_barcodes.size(); i++)
//            {
//                if (mc_vertex_incoming_barcodes[i].size() < 1)
//                {
//                    continue;
//                }
//                if (mc_vertex_incoming_barcodes[i][0] == truth_axions[0].mc_barcode)
//                {
//                    decayLength_a1 = sqrt(pow((mc_vertex_x[i] - startX),2) + pow((mc_vertex_y[i] - startY),2));
//                    break;
//                }
//            }
//            return decayLength_a1;
//        }, {"truth_higgs", "truth_axions", "mc_vertex_incoming_barcodes", "mc_vertex_x", "mc_vertex_y"})
//        .Define("axion_masses", [&](RVec<TruthParticle>& truth_axions)
//        {
//            return truth_axions[0].mc_mass/1e3f;
//
//        }, {"truth_axions"});
//
//        Nodes.push_back(truth_photons_from_axions.Histo1D<double>({"", "", 120u, 0, 0}, "ALP_decay_length"));
//
////        if (count <= 2)
////        {
////            Nodes.push_back(truth_photons_from_axions.Histo1D<double>({"", "", 120u, 0.2, 5.8}, "diphoton_inv_mass"));
////        }
////        else
////        {
////            for (auto& mass_point: massPoints)
////            {
////                auto mass_point_truth_photons_from_axions = truth_photons_from_axions.Filter([&]
////                (float axion_mass)
////                {
////                    return (roundToOneDecimalPlace(axion_mass) == mass_point);
////
////                }, {"axion_masses"});
////
////                Nodes.push_back(mass_point_truth_photons_from_axions.Histo1D<double>({"", "", 120u, 0.2, 5.8}, "diphoton_inv_mass"));
////            }
////        }
//        count++;
//    }
//
//    std::cout << Nodes.size() << '\n';
//    ROOT::RDF::RunGraphs(Nodes); // running all computation nodes concurrently
//
//    constexpr std::array<const char*, 4> titles = {"C_{a#gamma#gamma} = 1", "C_{a#gamma#gamma} = 0.01", "C_{a#gamma#gamma} = 0.001", "C_{a#gamma#gamma} = 0.0001"};
//    constexpr std::array<const char*, 4> filenames = {"Cayy_1_ALP_decay_lengths.pdf", "Cayy_0p01_ALP_decay_lengths.pdf", "Cayy_0p001_ALP_decay_lengths.pdf", "Cayy_0p0001_ALP_decay_lengths.pdf"};
//    constexpr std::array<EColor, 4> colors = {kBlue, static_cast<EColor>(kGreen + 2), kViolet, kRed};
//    TCanvas* c1;
//    TLegend* legend;
//    TLatex Tl;
//
//    for (int i = 0; i < Nodes.size(); i++)
//    {
//        c1 = new TCanvas();
//        legend = new TLegend(0.275, 0.55, 0.525, 0.8);
//        gStyle->SetOptStat(0);
//        Nodes[i].GetResultPtr<TH1D>()->SetLineColor(colors[i]);
//        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), titles[i], "l");
//
//        Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->SetTitle("ALP Decay Lengths  (mm)");
//        Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->SetTitleOffset(1.4);
//        Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->SetTitle("Events");
//        Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);
//        Nodes[i].GetResultPtr<TH1D>()->Draw("HISTsame");
//
//        Tl.SetTextSize(0.03);
//        Tl.DrawLatexNDC(0.6, 0.83, "#it{ATLAS} Internal");
//        Tl.DrawLatexNDC(0.6, 0.73,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
//        legend->SetBorderSize(0);
//        legend->Draw();
//        c1->SaveAs(filenames[i]);
//    }
//}


void ControlPlotsSignalShapes()
{
    auto start_time = Clock::now();
//    fig1A();
//    fig3();
//    fig5();
//    fig6();
//    fig8();
//    fig10();
//    fig18();
//    fig24();
//    fig54();
//    fig29();
//    AxionMassEventCounter();
    SignalShapes();
    SignalShapesDeltaR();
//    ALP_Decay_Lengths();
    
    auto end_time = Clock::now();
    std::cout << "Time difference: "
       << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count()/1e9 << " seconds" << std::endl;
    
}


int main()
{
    ControlPlotsSignalShapes();
}

