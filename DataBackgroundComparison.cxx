#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <iomanip>
#include <array>
#include <fstream>
#include <sstream>

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

float roundToOneDecimalPlace(float num) {
    std::stringstream stream;
    stream << std::fixed << std::setprecision(1) << num;
    float rounded_num = std::stof(stream.str());
    return rounded_num;
}

void fig26()
{
    std::vector<std::vector<std::string>> input_filenames = {
        //Signal
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root"}, //1 GeV
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600752.PhPy8EG_AZNLO_ggH125_mA3p0_v1.root"}, //2 GeV
        //Z gamma background
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617070._000001.LGNTuple.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617064._000001.LGNTuple.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617074._000001.LGNTuple.root"},
        //Data
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_data_test.root"},
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

    std::array<double,3> SFs = {((139e15)*(.871e-12)),((139e15)*(.199e-12)), ((139e15)*(.0345e-15))}; //numerators for Z-gamma bkg

    std::array<double,9> JetNumeratorSFs = {((139e15)*(1.9828e-9)*(0.821204)),((139e15)*(110.64e-12)*(0.69275)),((139e15)*(40.645e-12)*(0.615906)),((139e15)*(1.9817e-9)*(0.1136684)),((139e15)*(110.47e-12)*(0.1912956)),((139e15)*(40.674e-12)*(0.2326772)),((139e15)*(1.9819e-9)*(0.0656969)),((139e15)*(110.53e-12)*(0.1158741)),((139e15)*(40.68e-12)*(0.1535215))}; //numerators for jet bkg

    std::vector<const char*> prefixes = {"Sig m_{A} = 1 GeV", "Sig m_{A} = 2 GeV", "pty2_9_17", "pty_17_myy_0_80", "pty_17_myy_80", "data", "Zee_lightJet_0-70", "Zee_lightJet_70-140", "Zee_lightJet_140-280", "Zee_cJet_0-70", "Zee_cJet_70-140", "Zee_cJet_140-280", "Zee_bJet_0-70", "Zee_bJet_70-140", "Zee_bJet_140-280"};
    std::vector<EColor> colors = {kMagenta, static_cast<EColor>(kOrange+1), static_cast<EColor>(kRed-7), static_cast<EColor>(kRed-9), static_cast<EColor>(kRed-10), kBlack};
    std::vector<EColor> Jetscolors = {static_cast<EColor>(kGreen+0), static_cast<EColor>(kGreen-4), static_cast<EColor>(kGreen-7), static_cast<EColor>(kGreen-9), static_cast<EColor>(kGreen-10), static_cast<EColor>(kGreen+1), static_cast<EColor>(kGreen-3), static_cast<EColor>(kGreen-6), static_cast<EColor>(kGreen-8)};

    std::vector<ROOT::RDF::RResultHandle> Nodes;
    ROOT::RDF::RResultPtr<ULong64_t> SB_test, SR_test;

    int count = 0;
    for (auto& i: input_filenames)
    {
        SchottDataFrame df(MakeRDF(i, 8));

        auto EventWeight = df.Define("EventWeight", //mc generator weight
        [](const RVec<float>& ei_event_weights_generator)
        {
            return ((ei_event_weights_generator[0]) ? 1 / ei_event_weights_generator[0] : 1);

        }, {"ei_event_weights_generator"})
        .Define("totEventWeightVec", //vector of efficiencies, product of photon_id_eff, photon_iso_eff, and photon_trg_eff
        [](RVec<float> photon_id_eff/*, RVec<float> photon_iso_eff, RVec<float> photon_trg_eff, RVec<float> ei_event_weights_generator*/)
        {
//            auto ResizeVal = std::max({photon_id_eff.size(), photon_iso_eff.size(), photon_trg_eff.size()});
//            photon_id_eff.resize(ResizeVal,1);
//            photon_iso_eff.resize(ResizeVal,1);
//            photon_trg_eff.resize(ResizeVal,1);

            return photon_id_eff; // *photon_iso_eff*photon_trg_eff; //Returns a vector of efficiencies

        }, {"photon_id_eff"/*, "photon_iso_eff", "photon_trg_eff", "ei_event_weights_generator"*/});

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

        auto photon_passes_cuts = ptCut
        .Define("photons_pass_cut_indices",
        [&](const RVec<Photon>& photons)
        {
            RVec<int> photon_indices;
            photon_indices.reserve(photons.size());

            for (int i = 0; i < photons.size(); i++)
            {
                if (
                (std::abs(photons[i].photon_eta) >= 2.37) or
                (photons[i].photon_pt <= 10e3) or
                (std::abs(photons[i].photon_eta) > 1.37 and std::abs(photons[i].photon_eta) < 1.52) or
                (not photons[i].photon_id_loose)
                )
                {
                    continue;
                }
                photon_indices.push_back(i);
            }

            return photon_indices;
        }, {"photons"})
        .Define("photons_pass_cuts",
        [&](RVec<Photon>& photons, RVec<int>& photon_indices)
        {
            return Take(photons, photon_indices);
        }, {"photons", "photons_pass_cut_indices"});

        auto resolved = photon_passes_cuts.Define("chosen_two_indices",
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
        .Define("diphoton",
        [&](RVec<Photon>& reco_photons_matched)
        {
            return reco_photons_matched[0].Vector() + reco_photons_matched[1].Vector();
        }, {"chosen_two"})
        .Define("reco_higgs_mass",
        [&](PtEtaPhiEVector& diphoton, PtEtaPhiEVector& dilep)
        {
            return (dilep + diphoton).M() / 1e3;
        }, {"diphoton", "dilep"})
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

        if (count < 2) //signal only
        {
            Nodes.push_back(resolved.Histo1D<double>({prefixes[count], prefixes[count++], 100u, 80, 200}, "reco_higgs_mass", "totEventWeight"));
            if (count == 1)
            {
                auto SB_test_Filter = resolved.Filter(
                [](double reco_higgs_mass)
                {
                    return (!((reco_higgs_mass > 110) && (reco_higgs_mass < 140)));
                }, {"reco_higgs_mass"});
                SB_test = SB_test_Filter.Count();

                auto SR_test_Filter = resolved.Filter(
                [](double reco_higgs_mass)
                {
                    return ((reco_higgs_mass > 110) && (reco_higgs_mass < 140));
                }, {"reco_higgs_mass"});
                SR_test = SR_test_Filter.Count();
            }
        }
        else
        {
            if ((count >= 2 && count <= 4) || (count >= 6))
            {
                Nodes.push_back(resolved.Sum<float>("totEventWeight"));
                Nodes.push_back(EventWeight.Sum<float>("EventWeight"));
            }

            Nodes.push_back(resolved.Histo1D<double>({prefixes[count], prefixes[count++], 100u, 80, 200}, "reco_higgs_mass", "totEventWeight"));
        }
    }

//    0               ma1
//    1               ma2
//    2   3   4       Z-gamma
//    5   6   7       Z-gamma
//    8   9   10      Z-gamma
//    11              data
//    12  13  14      Z-jets
//    15  16  17      Z-jets
//    18  19  20      Z-jets
//    21  22  23      Z-jets
//    24  25  26      Z-jets
//    27  28  29      Z-jets
//    30  31  32      Z-jets
//    33  34  35      Z-jets
//    36  37  38      Z-jets

    ROOT::RDF::RunGraphs(Nodes); // running all computation nodes concurrently

    std::cout << "Sideband = " << *SB_test << '\n';
    std::cout << "Signal region = " << *SR_test << '\n';

    count = 0;
    TCanvas* c1 = new TCanvas();
    TLegend* legend = new TLegend(0.65, 0.4, 0.85, 0.6);
    double factor;
    //signal
    for (auto& i: {0,1})
    {
        Nodes[i].GetResultPtr<TH1D>()->SetLineColor(colors[count++]);
        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "l");

        if (i == 0)
        {
            factor = Nodes[i].GetResultPtr<TH1D>()->Integral();
            Nodes[i].GetResultPtr<TH1D>()->Scale(factor/Nodes[i].GetResultPtr<TH1D>()->Integral());
            Nodes[i].GetResultPtr<TH1D>()->SetTitle(";m_{ll#gamma#gamma} [GeV];Events");
            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);
            Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->SetTitleOffset(1.2);
            Nodes[i].GetResultPtr<TH1D>()->SetAxisRange(0., 71, "Y");
            Nodes[i].GetResultPtr<TH1D>()->Draw("HIST");
        }
        else
        {
            Nodes[i].GetResultPtr<TH1D>()->Scale(factor/Nodes[i].GetResultPtr<TH1D>()->Integral());
            Nodes[i].GetResultPtr<TH1D>()->Draw("HISTsame");
            gPad->Modified(); gPad->Update();
        }
    }

    gStyle->SetOptStat(0);
    TLatex Tl;
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.6, 0.8, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.6, 0.7,"#sqrt{s} = 13 TeV");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("Fig26B.pdf");

//    count = 0;
    factor = 0;
    int back_count = 0;
    for (auto i: {2,5,8})
    {
        factor += (*Nodes[i].GetResultPtr<float>())*(SFs[back_count++] / *Nodes[i+1].GetResultPtr<float>());
    }
    back_count = 0;
    for (int i = 12; i <= 36; i+=3)
    {
        factor += (*Nodes[i].GetResultPtr<float>())*(JetNumeratorSFs[back_count++] / *Nodes[i+1].GetResultPtr<float>());
    }

    auto hs = new THStack("hs3","");
    c1 = new TCanvas();
    legend = new TLegend(0.14, 0.45, 0.35, 0.8);

    //Z-gamma (4, 7, 10) and then data (11)

    Double_t data_integral = Nodes[11].GetResultPtr<TH1D>()->Integral();
    for (auto i: {4,7,10,11})
    {
        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0 && i != 11)
        {
            Nodes[i].GetResultPtr<TH1D>()->Scale(SFs[count-2] / *Nodes[i-1].GetResultPtr<float>());
            Nodes[i].GetResultPtr<TH1D>()->Scale(data_integral / factor);
        }

        if (i != 11)
        {
            Nodes[i].GetResultPtr<TH1D>()->SetLineColor(colors[count]);
            Nodes[i].GetResultPtr<TH1D>()->SetFillColor(colors[count++]);
            legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "f");
            hs->Add(&*Nodes[i].GetResultPtr<TH1D>());
        }
        else
        {
            Nodes[i].GetResultPtr<TH1D>()->SetLineColor(colors[count++]);
            legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "l");
        }
    }

    for (int i = 14, j = 0; (i <= 38 && j <= 8); i += 3, j++)
    {
        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0)
        {
            Nodes[i].GetResultPtr<TH1D>()->Scale(JetNumeratorSFs[j] / *Nodes[i-1].GetResultPtr<float>());
            Nodes[i].GetResultPtr<TH1D>()->Scale(data_integral / factor);
        }
        Nodes[i].GetResultPtr<TH1D>()->SetLineColor(Jetscolors[j]);
        Nodes[i].GetResultPtr<TH1D>()->SetFillColor(Jetscolors[j]);
        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "f");
        hs->Add(&*Nodes[i].GetResultPtr<TH1D>());
    }

    hs->Draw("HIST");
    Nodes[11].GetResultPtr<TH1D>()->Draw("E1same"); //draw data
    count=0;

    //signal: ma1, ma2
    for (auto i: {0,1})
    {
//        Nodes[i].GetResultPtr<TH1D>()->SetLineWidth(2);
        Nodes[i].GetResultPtr<TH1D>()->SetLineColor(colors[count]);
        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), prefixes[count++], "l");
        Nodes[i].GetResultPtr<TH1D>()->SetTitle(";m_{ll#gamma#gamma} [GeV];Events");
        Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);
        Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->SetTitleOffset(1.2);
        Nodes[i].GetResultPtr<TH1D>()->Scale(data_integral/Nodes[i].GetResultPtr<TH1D>()->Integral());
        std::cout << Nodes[i].GetResultPtr<TH1D>()->Integral() << '\n'
        << factor << '\n';
        Nodes[i].GetResultPtr<TH1D>()->DrawClone("HISTsame");
        gPad->Modified();
        gPad->Update();
        c1->Modified();
        c1->Update();
    }
    hs->SetMaximum(2);
    hs->SetTitle(";m_{ll#gamma#gamma} [GeV];Events");
    hs->GetYaxis()->CenterTitle(true);
    hs->GetXaxis()->SetTitleOffset(1.2);

    gStyle->SetOptStat(0);
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.6, 0.8, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.6, 0.7,"#sqrt{s} = 13 TeV");
    legend->SetBorderSize(0);
    legend->Draw("same");
    c1->SaveAs("Fig26A.pdf");
}

//void fig27()
//{
//    auto hs = new THStack("hs","");
//
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
//    std::vector<EColor> colors = {kBlue, kRed, kViolet, static_cast<EColor>(kGreen+2)};
//    std::vector<EColor> Jetscolors = {kCyan, kOrange, kGreen, kYellow, kPink, kGray, kBlack, kSpring, kAzure};
//
//    std::vector<ROOT::RDF::RResultHandle> Nodes;
//    int count = 0;
//    for (auto& i: input_filenames)
//    {
//        SchottDataFrame df(MakeRDF(i, 8));
//
//        auto EventWeight = df.Define("EventWeight", //mc generator weight
//        [](const RVec<float>& ei_event_weights_generator)
//        {
//            return ((ei_event_weights_generator[0]) ? 1 / ei_event_weights_generator[0] : 1);
//
//        }, {"ei_event_weights_generator"})
//        .Define("totEventWeightVec", //vector of efficiencies, product of photon_id_eff, photon_iso_eff, and photon_trg_eff
//        [](RVec<float> photon_id_eff/*, RVec<float> photon_iso_eff, RVec<float> photon_trg_eff, RVec<float> ei_event_weights_generator*/)
//        {
////            auto ResizeVal = std::max({photon_id_eff.size(), photon_iso_eff.size(), photon_trg_eff.size()});
////            photon_id_eff.resize(ResizeVal,1);
////            photon_iso_eff.resize(ResizeVal,1);
////            photon_trg_eff.resize(ResizeVal,1);
//
//            return photon_id_eff; // *photon_iso_eff*photon_trg_eff; //Returns a vector of efficiencies
//
//        }, {"photon_id_eff"/*, "photon_iso_eff", "photon_trg_eff", "ei_event_weights_generator"*/});
//
//        //new dataframe node: contains only the events of newDf that pass the trigger cut
//        auto trigger_selection = EventWeight
//        .Filter([](const RVec<std::string>& trigger_passed_triggers)
//        {
//            bool trigger_found = (std::find_first_of(trigger_passed_triggers.begin(), trigger_passed_triggers.end(), triggers.begin(), triggers.end()) != trigger_passed_triggers.end());
//
//            if (!trigger_found)
//            {
//                return false; //this event is filtered out
//            }
//            return true; //this event is kept because the trigger was found in its `trigger_passed_triggers` branch entry
//
//        }, {"trigger_passed_triggers"});
//
//        auto two_leptons = trigger_selection
//        .Define("di_electrons", //the events that pass will have exactly 2 electrons that pass the following
//        [](RVec<Electron> electrons)
//        {
//            //keep the electrons in each event that have pt > 20 GeV, |η| < 2.37,
//            //|η| not between 1.37 and 1.52, and that satisfy a medium id criteria `electron_id_medium`
//            electrons.erase(std::remove_if(electrons.begin(),electrons.end(),
//            [](Electron& ep)
//            {
//                 return (!((ep.electron_pt/1e3 > 20) && (std::abs(ep.electron_eta) < 2.37) &&
//                 (!((1.37 < std::abs(ep.electron_eta)) && (std::abs(ep.electron_eta) < 1.52)))
//                 && (ep.electron_id_medium == 1)));
//
//            }), electrons.end());
//
//            return electrons;
//
//        },{"electrons"}).Filter(
//        [](RVec<Muon>& muons, RVec<Electron> electrons)
//        {
//            return (electrons.size()==2 && muons.empty()); //keep events which have exactly 2 electrons for di_electrons and no muons
//
//        }, {"muons", "di_electrons"});
//
//        //new dataframe node: contains only the events from `two_leptons` whose electrons in the `di_electrons` branch have opposite charge
//        auto opp_charge = two_leptons
//        .Filter([](RVec<Electron>& electrons)
//        {
//            return (electrons[0].electron_charge*electrons[1].electron_charge < 0);
//
//        }, {"di_electrons"});
//
//        //new dataframe node: contains only the events from `opp_charge` that have 1 electron with pt > 20 GeV and the other with pt > 27 GeV
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
//        //new dataframe node: Contains a new column `dilep` in addition to the ones in `same_flavour` that stores the di-electron four-vector
//        auto dilep = deltaR.Define("dilep",[] (RVec<Electron>& electrons)
//        {
//            return (electrons[0].Vector() + electrons[1].Vector());
//        }, {"di_electrons"});
//
//        //new dataframe node: contains only the events from `dilep` that have di-electron invariant mass between 81 and 101 GeV
//        auto mass = dilep.Filter([] (PtEtaPhiEVector& dilep)
//        {
//            auto mass = dilep.M()/1e3;
//            return ((mass >= 81) && (mass <= 101));
//        }, {"dilep"});
//
//        //new dataframe node: contains only the events from `mass` that have dilepton pT > 10 GeV
//        auto ptCut = mass.Filter([] (PtEtaPhiEVector& dilep)
//        {
//            auto pT = dilep.Pt()/1e3;
//            return pT > 10;
//        }, {"dilep"});
//
//        //end pre-selection -----------
//
//        //photon acceptance and id_loose cuts
//        auto photons_pass_cuts = ptCut
//        .Define("photons_pass_cut_indices", //new column that contains the good photon indices
//        [&](RVec<Photon>& p) //p = photon
//        {
//            RVec<int> x; //indices of photons that passed the cuts
//            for (auto i = 0; i < p.size(); i++)
//            {
//                //keep reco-photons that have |η| < 2.37, p_T > 10 GeV, and |η| not between 1.37 and 1.52
//                if (not ((std::abs(p[i].photon_eta) >= 2.37) or (p[i].photon_pt <= 10e3) or (std::abs(p[i].photon_eta) > 1.37 and std::abs(p[i].photon_eta) < 1.52)))
//                {
//                    x.push_back(i);
//                }
//            }
//            return x;
//
//        }, {"photons"})
//        .Define("photons_pass_cuts", //new column that contains the good photons corresponding to the good photon indices from above
//        [&](RVec<Photon>& photons, RVec<int>& x)
//        {
//            return Take(photons, x); //Taking only the photons that passed the cuts in each event
//
//        }, {"photons", "photons_pass_cut_indices"});
//
//        auto resolved = photons_pass_cuts.Define("chosen_two_indices",
//        [](RVec<Photon>& photons_pass_cuts)
//        {
//            RVec<unsigned long> x; //vector of indices
//            if (photons_pass_cuts.size() < 2)
//            {
//                return x;
//            }
//
//            auto combs = Combinations(photons_pass_cuts, 2); //all combinations of 2 reco-photons
//            size_t length = combs[0].size(); //number of combinations
//            double delta_r, m, pt, X, best_X, pt1, pt2, chosen_delta_r;
//
//            for (size_t i=0; i<length; i++) //looping through all of the possible combinations of photons in each event
//            {
//                delta_r = DeltaR(photons_pass_cuts[combs[0][i]].Vector(), photons_pass_cuts[combs[1][i]].Vector());
//                m = (photons_pass_cuts[combs[0][i]].Vector() + photons_pass_cuts[combs[1][i]].Vector()).M();
//                pt = (photons_pass_cuts[combs[0][i]].Vector() + photons_pass_cuts[combs[1][i]].Vector()).Pt();
//                X = delta_r*(pt/(2.0*m));
//                //if it's the first combination or if new X is closer to 1
//                //than current best_X and ΔR
//                //between the two reco-photons < 1.5, then update best_X, pt1, pt2,
//                //and the corresponding reco-photon indices x
//                if (i==0 || ((std::abs(1-X) < std::abs(1-best_X)) and (delta_r < 1.5)))
//                {
//                    best_X = X;
//                    pt1 = photons_pass_cuts[combs[0][i]].photon_pt;
//                    pt2 = photons_pass_cuts[combs[1][i]].photon_pt;
//                    chosen_delta_r = delta_r;
//                    x = {combs[0][i], combs[1][i]};
//                }
//            }
//            if (chosen_delta_r < 1.5 && pt1 > 10e3 && pt2 > 10e3 && best_X > 0.96 && best_X < 1.2) //two photons corresponding to best_X must both have p_T > 10 GeV, ΔR < 1.5, and 0.96 < best_X < 1.2
//            {
//                return x;
//            }
//
//            x.clear();
//            return x;
//        }, {"photons_pass_cuts"}).Filter(//keep only events that have passed the first part of the resolved category
//        [&](RVec<unsigned long>& indices)
//        {
//            return (indices.size()==2);
//
//        }, {"chosen_two_indices"})
//        .Define("chosen_two", //New column: consists of the good photons corresponding to the `chosen_two_indices` defined above
//        [&](RVec<Photon>& reco_photons_matched, RVec<unsigned long>& indices)
//        {
//            return Take(reco_photons_matched, indices);
//
//        }, {"photons_pass_cuts", "chosen_two_indices"})
//        .Define("diphoton",
//        [&](RVec<Photon>& reco_photons_matched)
//        {
//            return reco_photons_matched[0].Vector() + reco_photons_matched[1].Vector();
//        }, {"chosen_two"})
//        .Define("reco_higgs_mass",
//        [&](PtEtaPhiEVector& diphoton, PtEtaPhiEVector& dilep)
//        {
//            return (dilep + diphoton).M() / 1e3;
//        }, {"diphoton", "dilep"})
//        .Filter(
//        [](double massVal)
//        {
//            return (!((massVal > 110) && (massVal < 140))); //sideband region
//        }, {"reco_higgs_mass"})
//        .Define("diphoton_mass",
//        [&](PtEtaPhiEVector& diphoton)
//        {
//            return diphoton.M() / 1e3;
//        }, {"diphoton"})
//        .Define("totEventWeight", //New column: weight factor for events in RDF `resolved`
//        [&](RVec<unsigned long> chosen_two_indices, RVec<int>& photons_pass_cut_indices, RVec<float>& photon_efficiencies)
//        {
//            //First, we take the elements from photon_efficiencies that correspond to the `photons_pass_cut_indices` defined
//            //earlier. Then, from that resulting vector, we take the elements corresponding to the
//            //`chosen_two_indices` defined above
//            RVec<float> resolved_photon_efficiencies = Take(Take(photon_efficiencies, photons_pass_cut_indices), chosen_two_indices);
//            //Now, multiply all of the elements of the vector we defined above
//            float total = 1.0f;
//            for (auto i: resolved_photon_efficiencies)
//            {
//                total *= i;
//            }
//            //and return the result
//            return total;
//        }, {"chosen_two_indices", "photons_pass_cut_indices", "totEventWeightVec"});
//
//        if (count <= 2 || count >= 4)
//        {
//            Nodes.push_back(resolved.Sum<float>("totEventWeight"));
//            Nodes.push_back(EventWeight.Sum<float>("EventWeight"));
//        }
//
//        Nodes.push_back(resolved.Histo1D<double>({prefixes[count], prefixes[count++], 60u, 0, 12}, "diphoton_mass", "totEventWeight"));
//
//    }
//
////    0   1   2   Z-gamma
////    3   4   5   Z-gamma
////    6   7   8   Z-gamma
////    9           data
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
//    int back_count = 0;
//    for (auto& i: {0,3,6})
//    {
//        factor += (*Nodes[i].GetResultPtr<float>())*(SFs[back_count++] / *Nodes[i+1].GetResultPtr<float>());
//    }
//
//    for (int i = 0, j = 10; (i <= 8 && j <= 34); i++, j+=3)
//    {
//        factor += (*Nodes[j].GetResultPtr<float>())*(JetNumeratorSFs[i] / *Nodes[j+1].GetResultPtr<float>());
//    }
//
//    TCanvas* c1 = new TCanvas();
//    TLegend* legend = new TLegend(0.325, 0.4, 0.6, 0.8);
//    count = 0;
//    for (auto& i: {2,5,8,9})
//    {
//        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0 && i != 9)
//        {
//            Nodes[i].GetResultPtr<TH1D>()->Scale(SFs[count] / *Nodes[i-1].GetResultPtr<float>() );
//        }
//        Nodes[i].GetResultPtr<TH1D>()->SetFillColor(colors[count++]);
//        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "f");
//
//        if (i != 9)
//        {
//            hs->Add(&*Nodes[i].GetResultPtr<TH1D>());
//        }
//    }
//    count = 0;
//    //Z-jets
//    for (int i = 12; i <= 36; i+=3)
//    {
//        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0)
//        {
//            Nodes[i].GetResultPtr<TH1D>()->Scale(JetNumeratorSFs[count] / *Nodes[i-1].GetResultPtr<float>() );
//        }
//        Nodes[i].GetResultPtr<TH1D>()->SetFillColor(Jetscolors[count++]);
//        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "f");
//
//        hs->Add(&*Nodes[i].GetResultPtr<TH1D>());
//    }
//
//    hs->Draw("HIST");
//    Nodes[9].GetResultPtr<TH1D>()->Draw("HISTsame");
//    hs->SetMinimum(0);
//    hs->SetMaximum(16300);
//    hs->SetTitle(";m_{#gamma#gamma} [GeV];Events");
//    hs->GetYaxis()->CenterTitle(true);
//    hs->GetYaxis()->SetTitleOffset(1.4);
//    hs->GetXaxis()->SetTitleOffset(1.2);
//    gStyle->SetOptStat(0);
//    TLatex Tl;
//    Tl.SetTextSize(0.03);
//    Tl.DrawLatexNDC(0.6, 0.8, "#it{ATLAS} Internal");
//    Tl.DrawLatexNDC(0.6, 0.7,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
//    legend->SetBorderSize(0);
//    legend->Draw();
//    c1->SaveAs("Fig27.pdf");
//}

void fig28()
{
    std::vector<std::vector<std::string>> input_filenames = {
        //Z gamma background
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617070._000001.LGNTuple.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617064._000001.LGNTuple.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617074._000001.LGNTuple.root"},
        //Data
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_data_test.root"},
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

    std::array<double,3> SFs = {((139e15)*(.871e-12)),((139e15)*(.199e-12)), ((139e15)*(.0345e-15))}; //numerators for Z-gamma bkg

    std::array<double,9> JetNumeratorSFs = {((139e15)*(1.9828e-9)*(0.821204)),((139e15)*(110.64e-12)*(0.69275)),((139e15)*(40.645e-12)*(0.615906)),((139e15)*(1.9817e-9)*(0.1136684)),((139e15)*(110.47e-12)*(0.1912956)),((139e15)*(40.674e-12)*(0.2326772)),((139e15)*(1.9819e-9)*(0.0656969)),((139e15)*(110.53e-12)*(0.1158741)),((139e15)*(40.68e-12)*(0.1535215))}; //numerators for jet bkg

    std::vector<const char*> prefixes = {"pty2_9_17", "pty_17_myy_0_80", "pty_17_myy_80", "data", "Zee_lightJet_0-70", "Zee_lightJet_70-140", "Zee_lightJet_140-280", "Zee_cJet_0-70", "Zee_cJet_70-140", "Zee_cJet_140-280", "Zee_bJet_0-70", "Zee_bJet_70-140", "Zee_bJet_140-280"};
    std::vector<EColor> colors = {static_cast<EColor>(kRed-7), static_cast<EColor>(kRed-9), static_cast<EColor>(kRed-10), kBlack};
    std::vector<EColor> Jetscolors = {static_cast<EColor>(kGreen+0), static_cast<EColor>(kGreen-4), static_cast<EColor>(kGreen-7), static_cast<EColor>(kGreen-9), static_cast<EColor>(kGreen-10), static_cast<EColor>(kGreen+1), static_cast<EColor>(kGreen-3), static_cast<EColor>(kGreen-6), static_cast<EColor>(kGreen-8)};

    std::vector<ROOT::RDF::RResultHandle> Nodes;

    int count = 0;
    for (auto& i: input_filenames)
    {
        SchottDataFrame df(MakeRDF(i, 8));

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

            return photon_id_eff;
            //*photon_iso_eff*photon_trg_eff; //Returns a vector of efficiencies

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
        .Define("diphoton",
        [&](RVec<Photon>& reco_photons_matched)
        {
            return reco_photons_matched[0].Vector() + reco_photons_matched[1].Vector();
        }, {"chosen_two"})
        .Define("reco_higgs_mass",
        [&](PtEtaPhiEVector& diphoton, PtEtaPhiEVector& dilep)
        {
            return (dilep + diphoton).M() / 1e3;
        }, {"diphoton", "dilep"})
        .Filter(
        [](double massVal)
        {
            return (!((massVal > 110) && (massVal < 140))); //sideband region
        }, {"reco_higgs_mass"})
        .Define("mass",
        [&](PtEtaPhiEVector& diphoton)
        {
            return diphoton.M() / 1e3;
        }, {"diphoton"})
        .Define("deltaR",
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
//            return std::abs((diph[0].Vector() - diph[1].Vector()).Eta());
            return diph[0].Vector().Eta() - diph[1].Vector().Eta();
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

        if (count <= 2 || count >= 4)
        {
            Nodes.push_back(resolved.Sum<float>("totEventWeight"));
            Nodes.push_back(EventWeight.Sum<float>("EventWeight"));
        }

        Nodes.push_back(resolved.Histo1D<double>({prefixes[count], prefixes[count], 60u, 0, 12}, "mass", "totEventWeight"));
        Nodes.push_back(resolved.Histo1D<double>({prefixes[count], prefixes[count], 60u, 0, 0.25}, "deltaR", "totEventWeight"));
        Nodes.push_back(resolved.Histo1D<double>({prefixes[count], prefixes[count], 60u, 0, 0.2}, "deltaPhi", "totEventWeight"));
        Nodes.push_back(resolved.Histo1D<double>({prefixes[count], prefixes[count++], 60u, 0, 0.04}, "deltaEta", "totEventWeight"));
    }

//    0    1    2    3    4    5      Z-gamma
//    6    7    8    9    10   11     Z-gamma
//    12   13   14   15   16   17     Z-gamma
//              18   19   20   21     data
//    22   23   24   25   26   27     Z-jets
//    28   29   30   31   32   33     Z-jets
//    34   35   36   37   38   39     Z-jets
//    40   41   42   43   44   45     Z-jets
//    46   47   48   49   50   51     Z-jets
//    52   53   54   55   56   57     Z-jets
//    58   59   60   61   62   63     Z-jets
//    64   65   66   67   68   69     Z-jets
//    70   71   72   73   74   75     Z-jets

    ROOT::RDF::RunGraphs(Nodes); // running all computation nodes concurrently

    double factor;
    int back_count = 0;
    for (auto& i: {0, 6, 12})
    {
        factor += (*Nodes[i].GetResultPtr<float>())*(SFs[back_count++] / *Nodes[i+1].GetResultPtr<float>());
    }

    for (int i = 22, j = 0; (i <= 70 && j<= 8); i += 6, j++)
    {
        factor += (*Nodes[i].GetResultPtr<float>())*(JetNumeratorSFs[j] / *Nodes[i+1].GetResultPtr<float>());
    }

    auto hs = new THStack("hs1","");
    TCanvas* c1 = new TCanvas();
    TLegend* legend = new TLegend(0.325, 0.4, 0.6, 0.8);
    count = 0;
    //Z-gamma
    Double_t data_integral = Nodes[18].GetResultPtr<TH1D>()->Integral();

    for (auto& i: {2,8,14,18})
    {
        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0 && i != 18)
        {
            Nodes[i].GetResultPtr<TH1D>()->Scale( (SFs[count]/ *Nodes[i-1].GetResultPtr<float>())  );
            Nodes[i].GetResultPtr<TH1D>()->Scale(data_integral / factor);
        }
        Nodes[i].GetResultPtr<TH1D>()->SetLineColor(colors[count]);
        Nodes[i].GetResultPtr<TH1D>()->SetFillColor(colors[count++]);

        if (i != 18)
        {
            legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "f");
            hs->Add(&*Nodes[i].GetResultPtr<TH1D>());
        }
        else
        {
            legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "l");
        }
    }
    //Z-jets
    count = 0;
    for (int i = 24; i <= 72; i += 6)
    {
        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0)
        {
            Nodes[i].GetResultPtr<TH1D>()->Scale(JetNumeratorSFs[count]/ *Nodes[i-1].GetResultPtr<float>());
            Nodes[i].GetResultPtr<TH1D>()->Scale(data_integral / factor);
        }
        Nodes[i].GetResultPtr<TH1D>()->SetLineColor(Jetscolors[count]);
        Nodes[i].GetResultPtr<TH1D>()->SetFillColor(Jetscolors[count++]);
        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "f");

        hs->Add(&*Nodes[i].GetResultPtr<TH1D>());
    }

    hs->Draw("HIST");
    Nodes[18].GetResultPtr<TH1D>()->Draw("E1Same"); //data
    hs->SetMinimum(0);
//    hs->SetMaximum(16300);
    hs->SetTitle(";m_{#gamma#gamma} [GeV];Events");
    hs->GetYaxis()->CenterTitle(true);
    hs->GetYaxis()->SetTitleOffset(1.4);
    hs->GetXaxis()->SetTitleOffset(1.2);
    gStyle->SetOptStat(0);
    TLatex Tl;
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.6, 0.8, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.6, 0.7,"#sqrt{s} = 13 TeV");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("Fig28A.pdf");

    hs = new THStack("hs2","");
    c1 = new TCanvas();
    legend = new TLegend(0.575, 0.3, 0.875, 0.675);
    count = 0;
    data_integral = Nodes[19].GetResultPtr<TH1D>()->Integral();
    //Z-gamma
    for (auto& i: {3,9,15,19})
    {
        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0 && i != 19)
        {
            Nodes[i].GetResultPtr<TH1D>()->Scale(SFs[count]/ *Nodes[i-2].GetResultPtr<float>());
            Nodes[i].GetResultPtr<TH1D>()->Scale(data_integral / factor);
        }
        Nodes[i].GetResultPtr<TH1D>()->SetLineColor(colors[count]);
        Nodes[i].GetResultPtr<TH1D>()->SetFillColor(colors[count++]);

        if (i != 19)
        {
            legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "f");
            hs->Add(&*Nodes[i].GetResultPtr<TH1D>());
        }
        else
        {
            legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "l");
        }
    }

    //Z-jets
    count = 0;
    for (int i = 25; i <= 73; i += 6)
    {
        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0)
        {
            Nodes[i].GetResultPtr<TH1D>()->Scale(JetNumeratorSFs[count]/ *Nodes[i-2].GetResultPtr<float>());
            Nodes[i].GetResultPtr<TH1D>()->Scale(data_integral / factor);
        }
        Nodes[i].GetResultPtr<TH1D>()->SetLineColor(Jetscolors[count]);
        Nodes[i].GetResultPtr<TH1D>()->SetFillColor(Jetscolors[count++]);
        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "f");

        hs->Add(&*Nodes[i].GetResultPtr<TH1D>());
    }

    hs->Draw("HIST");
    Nodes[19].GetResultPtr<TH1D>()->Draw("E1Same"); //data
    hs->SetMinimum(0);
//    hs->SetMaximum(7000);
    hs->SetTitle(";#DeltaR_{#gamma#gamma};Events");
    hs->GetYaxis()->CenterTitle(true);
    hs->GetXaxis()->SetTitleOffset(1.2);
    hs->GetYaxis()->SetTitleOffset(1.4);
    gStyle->SetOptStat(0);
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.6, 0.85, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.6, 0.775,"#sqrt{s} = 13 TeV");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("Fig28B.pdf");

    hs = new THStack("hs3","");
    c1 = new TCanvas();
    legend = new TLegend(0.625, 0.55, 0.88, 0.85);
    count = 0;
    data_integral = Nodes[20].GetResultPtr<TH1D>()->Integral();
    //Z-gamma
    for (auto& i: {4,10,16,20})
    {
        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0 && i != 20)
        {
            Nodes[i].GetResultPtr<TH1D>()->Scale(SFs[count]/ *Nodes[i-3].GetResultPtr<float>());
            Nodes[i].GetResultPtr<TH1D>()->Scale(data_integral / factor);
        }
        Nodes[i].GetResultPtr<TH1D>()->SetLineColor(colors[count]);
        Nodes[i].GetResultPtr<TH1D>()->SetFillColor(colors[count++]);

        if (i != 20)
        {
            legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "f");
            hs->Add(&*Nodes[i].GetResultPtr<TH1D>());
        }
        else
        {
            legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "l");
        }
    }

    //Z-jets
    count = 0;
    for (int i = 26; i <= 74; i += 6)
    {
        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0)
        {
            Nodes[i].GetResultPtr<TH1D>()->Scale(JetNumeratorSFs[count]/ *Nodes[i-3].GetResultPtr<float>());
            Nodes[i].GetResultPtr<TH1D>()->Scale(data_integral / factor);
        }
        Nodes[i].GetResultPtr<TH1D>()->SetLineColor(Jetscolors[count]);
        Nodes[i].GetResultPtr<TH1D>()->SetFillColor(Jetscolors[count++]);
        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "f");

        hs->Add(&*Nodes[i].GetResultPtr<TH1D>());
    }

    hs->Draw("HIST");
    Nodes[20].GetResultPtr<TH1D>()->Draw("E1Same"); //data
    hs->SetMinimum(0);
//    hs->SetMaximum(3400);
    hs->SetTitle(";#Delta#phi_{#gamma#gamma};Events");
    hs->GetYaxis()->CenterTitle(true);
    hs->GetYaxis()->SetTitleOffset(1.4);
    hs->GetXaxis()->SetTitleOffset(1.2);
    gStyle->SetOptStat(0);
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.28, 0.85, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.28, 0.775,"#sqrt{s} = 13 TeV");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("Fig28C.pdf");

    hs = new THStack("hs4","");
    c1 = new TCanvas();
    legend = new TLegend(0.375, 0.4, 0.625, 0.75);
    data_integral = Nodes[21].GetResultPtr<TH1D>()->Integral();
    std::cout<<"data_integral = " << data_integral << '\n';
    count = 0;
    for (auto& i: {5,11,17,21})
    {
        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0 && i != 21)
        {
            Nodes[i].GetResultPtr<TH1D>()->Scale(SFs[count]/ *Nodes[i-4].GetResultPtr<float>());
            Nodes[i].GetResultPtr<TH1D>()->Scale(data_integral / factor);
        }
        Nodes[i].GetResultPtr<TH1D>()->SetLineColor(colors[count]);
        Nodes[i].GetResultPtr<TH1D>()->SetFillColor(colors[count++]);

        if (i != 21)
        {
            legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "f");
            hs->Add(&*Nodes[i].GetResultPtr<TH1D>());
        }
        else
        {
            legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "l");
        }
    }

    //Z-jets
    count = 0;
    for (int i = 27; i <= 75; i += 6)
    {
        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0)
        {
            Nodes[i].GetResultPtr<TH1D>()->Scale(JetNumeratorSFs[count]/ *Nodes[i-4].GetResultPtr<float>());
            Nodes[i].GetResultPtr<TH1D>()->Scale(data_integral / factor);
        }
        Nodes[i].GetResultPtr<TH1D>()->SetLineColor(Jetscolors[count]);
        Nodes[i].GetResultPtr<TH1D>()->SetFillColor(Jetscolors[count++]);
        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "f");

        hs->Add(&*Nodes[i].GetResultPtr<TH1D>());
    }

    hs->Draw("HIST");
    Nodes[21].GetResultPtr<TH1D>()->Draw("E1Same"); //data
    hs->SetMinimum(0);
    hs->SetMaximum(1);
    hs->SetTitle(";#Delta#eta_{#gamma#gamma};Events");
    hs->GetYaxis()->CenterTitle(true);
    hs->GetYaxis()->SetTitleOffset(1.4);
    hs->GetXaxis()->SetTitleOffset(1.2);
    gStyle->SetOptStat(0);
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.575, 0.85, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.575, 0.775,"#sqrt{s} = 13 TeV");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("Fig28D.pdf");
}

void fig41()
{
    std::vector<std::vector<std::string>> input_filenames = {
        //Signal
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root"}, //1 GeV
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600752.PhPy8EG_AZNLO_ggH125_mA3p0_v1.root"}, //2 GeV
        //Z gamma background
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617070._000001.LGNTuple.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617064._000001.LGNTuple.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617074._000001.LGNTuple.root"},
        //Data
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_data_test.root"},
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

    std::array<double,3> SFs = {((139e15)*(.871e-12)),((139e15)*(.199e-12)), ((139e15)*(.0345e-15))}; //numerators for Z-gamma bkg

    std::array<double,9> JetNumeratorSFs = {((139e15)*(1.9828e-9)*(0.821204)),((139e15)*(110.64e-12)*(0.69275)),((139e15)*(40.645e-12)*(0.615906)),((139e15)*(1.9817e-9)*(0.1136684)),((139e15)*(110.47e-12)*(0.1912956)),((139e15)*(40.674e-12)*(0.2326772)),((139e15)*(1.9819e-9)*(0.0656969)),((139e15)*(110.53e-12)*(0.1158741)),((139e15)*(40.68e-12)*(0.1535215))}; //numerators for jet bkg

    std::vector<const char*> prefixes = {"Sig m_{A} = 1 GeV", "Sig m_{A} = 2 GeV", "pty2_9_17", "pty_17_myy_0_80", "pty_17_myy_80", "data", "Zee_lightJet_0-70", "Zee_lightJet_70-140", "Zee_lightJet_140-280", "Zee_cJet_0-70", "Zee_cJet_70-140", "Zee_cJet_140-280", "Zee_bJet_0-70", "Zee_bJet_70-140", "Zee_bJet_140-280"};
    std::vector<EColor> colors = {kMagenta, static_cast<EColor>(kOrange+1), static_cast<EColor>(kRed-7), static_cast<EColor>(kRed-9), static_cast<EColor>(kRed-10), kBlack};
    std::vector<EColor> Jetscolors = {static_cast<EColor>(kGreen+0), static_cast<EColor>(kGreen-4), static_cast<EColor>(kGreen-7), static_cast<EColor>(kGreen-9), static_cast<EColor>(kGreen-10), static_cast<EColor>(kGreen+1), static_cast<EColor>(kGreen-3), static_cast<EColor>(kGreen-6), static_cast<EColor>(kGreen-8)};

    std::vector<ROOT::RDF::RResultHandle> Nodes;

    int count = 0;
    for (auto& i: input_filenames)
    {
        SchottDataFrame df(MakeRDF(i, 8));

        auto EventWeight = df.Define("EventWeight", //mc generator weight
        [](const RVec<float>& ei_event_weights_generator)
        {
            return ((ei_event_weights_generator[0]) ? 1 / ei_event_weights_generator[0] : 1);

        }, {"ei_event_weights_generator"})
        .Define("totEventWeightVec", //vector of efficiencies, product of photon_id_eff, photon_iso_eff, and photon_trg_eff
        [](RVec<float> photon_id_eff/*, RVec<float> photon_iso_eff, RVec<float> photon_trg_eff, RVec<float> ei_event_weights_generator*/)
        {
//            auto ResizeVal = std::max({photon_id_eff.size(), photon_iso_eff.size(), photon_trg_eff.size()});
//            photon_id_eff.resize(ResizeVal,1);
//            photon_iso_eff.resize(ResizeVal,1);
//            photon_trg_eff.resize(ResizeVal,1);

            return photon_id_eff; // *photon_iso_eff*photon_trg_eff; //Returns a vector of efficiencies

        }, {"photon_id_eff"/*, "photon_iso_eff", "photon_trg_eff", "ei_event_weights_generator"*/});

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

        auto photon_passes_cuts = ptCut
        .Define("photons_pass_cut_indices",
        [&](const RVec<Photon>& photons)
        {
            RVec<int> photon_indices;
            photon_indices.reserve(photons.size());

            for (int i = 0; i < photons.size(); i++)
            {
                if (
                (std::abs(photons[i].photon_eta) >= 2.37) or
                (photons[i].photon_pt <= 10e3) or
                (std::abs(photons[i].photon_eta) > 1.37 and std::abs(photons[i].photon_eta) < 1.52) or
                (not photons[i].photon_id_loose)
                )
                {
                    continue;
                }
                photon_indices.push_back(i);
            }

            return photon_indices;
        }, {"photons"})
        .Define("photons_pass_cuts",
        [&](RVec<Photon>& photons, RVec<int>& photon_indices)
        {
            return Take(photons, photon_indices);
        }, {"photons", "photons_pass_cut_indices"});

        //New dataframe node: contains only events from `photon_passes_cuts` that fail the resolved category and pass the merged category
        auto merged_reco_photons_matched = photon_passes_cuts.Filter(
        [&](RVec<Photon>& reco_photons_matched)
        {
           if (reco_photons_matched.size() == 1) // 1 photon in the event
           {
               return reco_photons_matched[0].photon_pt > 20e3; //event passes if photon pt > 20 GeV
           }
           else if (reco_photons_matched.empty())
           {
               return false; //fails if no photons in event
           }

           auto combs = Combinations(reco_photons_matched, 2); //all combinations of 2 reco-photons
           size_t length = combs[0].size(); //number of combinations
           double delta_r, m, pt, X, best_X, pt1, pt2, chosen_delta_r;

           for (size_t i=0; i<length; i++) //looping through all of the possible combinations of photons in each event
           {
               delta_r = DeltaR(reco_photons_matched[combs[0][i]].Vector(), reco_photons_matched[combs[1][i]].Vector());
               m = (reco_photons_matched[combs[0][i]].Vector() + reco_photons_matched[combs[1][i]].Vector()).M();
               pt = (reco_photons_matched[combs[0][i]].Vector() + reco_photons_matched[combs[1][i]].Vector()).Pt();
               X = delta_r*(pt/(2.0*m));
               //if it's the first combination or if new X is closer to 1
               //than current best_X and ΔR
               //between the two reco-photons < 1.5, then update best_X, pt1, pt2,
               //and the corresponding reco-photon indices x
               if (i==0 || ((std::abs(1-X) < std::abs(1-best_X)) and (delta_r < 1.5)))
               {
                   best_X = X;
                   pt1 = reco_photons_matched[combs[0][i]].photon_pt;
                   pt2 = reco_photons_matched[combs[1][i]].photon_pt;
                   chosen_delta_r = delta_r;
               }
           }
           if (pt1 > 10e3 && pt2 > 10e3 && best_X > 0.96 && best_X < 1.2 && chosen_delta_r < 1.5) //two photons corresponding to best_X must both have p_T > 10 GeV, ΔR < 1.5, and 0.96 < best_X < 1.2 for resolved... but this is merged, so if it passes resolved it fails merged
           {
               return false;
           }
           //if we get to this point, it means we've failed resolved
           for (auto& p: reco_photons_matched) //merged means at least 1 reco photon must have pt > 20 GeV
           {
               if (p.photon_pt > 20e3)
               {
                   return true; //passed merged if there's a reco-photon with pt > 20 GeV
               }
           }
           return false; //failed merged

        }, {"photons_pass_cuts"})
        .Define("merged_photon_index", //new column: consists of the index corresponding to the photon that made the event be classified as merged
        [&](const RVec<Photon>& rpm) //rpm = reco photons matched
        {
            for (auto i = 0; i < rpm.size(); i++)
            {
                if (rpm[i].photon_pt > 20e3)
                {
                    return i; //returning the index of the first photon that has photon_pt > 20 GeV
                }
            }
            return 0; //jic the compiler complains, should not come to this

        }, {"photons_pass_cuts"})
        .Define("merged_photon", //new column: The reco-photon corresponding to `merged_photon_index`
        [&](const RVec<Photon>& reco_photons_matched, int merged_photon_index)
        {
            return reco_photons_matched[merged_photon_index];

        }, {"photons_pass_cuts", "merged_photon_index"});

        auto dilepton_and_photon = merged_reco_photons_matched
        .Define("reconstructed_mass",[&](RVec<Electron>& di_electrons, Photon& merged_photon)
        {
            auto four_momentum = di_electrons[0].Vector() + di_electrons[1].Vector();

            return (four_momentum + merged_photon.Vector()).M()/1e3;

        }, {"di_electrons", "merged_photon"})
        .Define("totEventWeight", [](RVec<float> totEventWeightVec, RVec<int>& photons_pass_cut_indices, int mpi)
        {
            return Take(totEventWeightVec, photons_pass_cut_indices)[mpi];

        }, {"totEventWeightVec", "photons_pass_cut_indices", "merged_photon_index"});

        if (count < 2) //signal only
        {
            Nodes.push_back(dilepton_and_photon.Histo1D<double>({prefixes[count], prefixes[count++], 100u, 80, 200}, "reconstructed_mass", "totEventWeight"));
        }
        else
        {
            if ((count >= 2 && count <= 4) || (count >= 6))
            {
                Nodes.push_back(dilepton_and_photon.Sum<float>("totEventWeight"));
                Nodes.push_back(EventWeight.Sum<float>("EventWeight"));
            }

            Nodes.push_back(dilepton_and_photon.Histo1D<double>({prefixes[count], prefixes[count++], 100u, 80, 200}, "reconstructed_mass", "totEventWeight"));
        }
    }

//    0               ma1
//    1               ma2
//    2   3   4       Z-gamma
//    5   6   7       Z-gamma
//    8   9   10      Z-gamma
//    11              data
//    12  13  14      Z-jets
//    15  16  17      Z-jets
//    18  19  20      Z-jets
//    21  22  23      Z-jets
//    24  25  26      Z-jets
//    27  28  29      Z-jets
//    30  31  32      Z-jets
//    33  34  35      Z-jets
//    36  37  38      Z-jets

    ROOT::RDF::RunGraphs(Nodes); // running all computation nodes concurrently

    count = 0;
    TCanvas* c1 = new TCanvas();
    TLegend* legend = new TLegend(0.65, 0.4, 0.85, 0.6);
    double factor;
    //signal
    for (auto& i: {0,1})
    {
        Nodes[i].GetResultPtr<TH1D>()->SetLineColor(colors[count++]);
        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "l");

        if (i == 0)
        {
            factor = Nodes[i].GetResultPtr<TH1D>()->Integral();
            Nodes[i].GetResultPtr<TH1D>()->Scale(factor/Nodes[i].GetResultPtr<TH1D>()->Integral());
            Nodes[i].GetResultPtr<TH1D>()->SetTitle(";m_{ll#gamma} [GeV];Events");
            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);
            Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->SetTitleOffset(1.2);
            Nodes[i].GetResultPtr<TH1D>()->SetAxisRange(0., 620, "Y");
            Nodes[i].GetResultPtr<TH1D>()->Draw("HIST");
        }
        else
        {
            Nodes[i].GetResultPtr<TH1D>()->Scale(factor/Nodes[i].GetResultPtr<TH1D>()->Integral());
            Nodes[i].GetResultPtr<TH1D>()->Draw("HISTsame");
            gPad->Modified(); gPad->Update();
        }
    }

    gStyle->SetOptStat(0);
    TLatex Tl;
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.6, 0.8, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.6, 0.7,"#sqrt{s} = 13 TeV");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("Fig41B.pdf");

//    count = 0;
    factor = 0;
    int back_count = 0;
    for (auto i: {2,5,8})
    {
        factor += (*Nodes[i].GetResultPtr<float>())*(SFs[back_count++] / *Nodes[i+1].GetResultPtr<float>());
    }
    back_count = 0;
    for (int i = 12; i <= 36; i+=3)
    {
        factor += (*Nodes[i].GetResultPtr<float>())*(JetNumeratorSFs[back_count++] / *Nodes[i+1].GetResultPtr<float>());
    }

    auto hs = new THStack("hs3","");
    c1 = new TCanvas();
    legend = new TLegend(0.5, 0.2, 0.85, 0.6);

    //Z-gamma (4, 7, 10) and then data (11)

    Double_t data_integral = Nodes[11].GetResultPtr<TH1D>()->Integral();
    for (auto i: {4,7,10,11})
    {
        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0 && i != 11)
        {
            Nodes[i].GetResultPtr<TH1D>()->Scale(SFs[count-2] / *Nodes[i-1].GetResultPtr<float>());
            Nodes[i].GetResultPtr<TH1D>()->Scale(data_integral / factor);
        }

        if (i != 11)
        {
            Nodes[i].GetResultPtr<TH1D>()->SetLineColor(colors[count]);
            Nodes[i].GetResultPtr<TH1D>()->SetFillColor(colors[count++]);
            legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "f");
            hs->Add(&*Nodes[i].GetResultPtr<TH1D>());
        }
        else
        {
            Nodes[i].GetResultPtr<TH1D>()->SetLineColor(colors[count++]);
            legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "l");
        }
    }

    for (int i = 14, j = 0; (i <= 38 && j <= 8); i += 3, j++)
    {
        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0)
        {
            Nodes[i].GetResultPtr<TH1D>()->Scale(JetNumeratorSFs[j] / *Nodes[i-1].GetResultPtr<float>());
            Nodes[i].GetResultPtr<TH1D>()->Scale(data_integral / factor);
        }
        Nodes[i].GetResultPtr<TH1D>()->SetLineColor(Jetscolors[j]);
        Nodes[i].GetResultPtr<TH1D>()->SetFillColor(Jetscolors[j]);
        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "f");
        hs->Add(&*Nodes[i].GetResultPtr<TH1D>());
    }

    hs->Draw("HIST");
    Nodes[11].GetResultPtr<TH1D>()->Draw("E1same"); //draw data
    count=0;

    //signal: ma1, ma2
    for (auto i: {0,1})
    {
//        Nodes[i].GetResultPtr<TH1D>()->SetLineWidth(2);
        Nodes[i].GetResultPtr<TH1D>()->SetLineColor(colors[count]);
        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), prefixes[count++], "l");
        Nodes[i].GetResultPtr<TH1D>()->SetTitle(";m_{ll#gamma} [GeV];Events");
        Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);
        Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->SetTitleOffset(1.2);
        Nodes[i].GetResultPtr<TH1D>()->Scale(data_integral/Nodes[i].GetResultPtr<TH1D>()->Integral());
        std::cout << Nodes[i].GetResultPtr<TH1D>()->Integral() << '\n'
        << factor << '\n';
        Nodes[i].GetResultPtr<TH1D>()->DrawClone("HISTsame");
        gPad->Modified();
        gPad->Update();
        c1->Modified();
        c1->Update();
    }
    hs->SetMaximum(0.9e2);
    hs->SetTitle(";m_{ll#gamma} [GeV];Events");
    hs->GetYaxis()->CenterTitle(true);
    hs->GetXaxis()->SetTitleOffset(1.2);

    gStyle->SetOptStat(0);
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.6, 0.8, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.6, 0.7,"#sqrt{s} = 13 TeV");
    legend->SetBorderSize(0);
    legend->Draw("same");
    c1->SaveAs("Fig41A.pdf");
}

void fig48()
{
    std::vector<std::vector<std::string>> input_filenames = {
        //Z gamma background
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617070._000001.LGNTuple.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617064._000001.LGNTuple.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617074._000001.LGNTuple.root"},
        //Data
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_data_test.root"},
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
    std::array<double,3> SFs = {((139e15)*(.871e-12)),((139e15)*(.199e-12)), ((139e15)*(.0345e-15))}; //numerators for Z-gamma bkg

    std::array<double,9> JetNumeratorSFs = {((139e15)*(1.9828e-9)*(0.821204)),((139e15)*(110.64e-12)*(0.69275)),((139e15)*(40.645e-12)*(0.615906)),((139e15)*(1.9817e-9)*(0.1136684)),((139e15)*(110.47e-12)*(0.1912956)),((139e15)*(40.674e-12)*(0.2326772)),((139e15)*(1.9819e-9)*(0.0656969)),((139e15)*(110.53e-12)*(0.1158741)),((139e15)*(40.68e-12)*(0.1535215))}; //numerators for jet bkg

    std::vector<const char*> prefixes = {"pty2_9_17", "pty_17_myy_0_80", "pty_17_myy_80", "data", "Zee_lightJet_0-70", "Zee_lightJet_70-140", "Zee_lightJet_140-280", "Zee_cJet_0-70", "Zee_cJet_70-140", "Zee_cJet_140-280", "Zee_bJet_0-70", "Zee_bJet_70-140", "Zee_bJet_140-280"};
    std::vector<EColor> colors = {static_cast<EColor>(kRed-7), static_cast<EColor>(kRed-9), static_cast<EColor>(kRed-10), kBlack};
    std::vector<EColor> Jetscolors = {static_cast<EColor>(kGreen+0), static_cast<EColor>(kGreen-4), static_cast<EColor>(kGreen-7), static_cast<EColor>(kGreen-9), static_cast<EColor>(kGreen-10), static_cast<EColor>(kGreen+1), static_cast<EColor>(kGreen-3), static_cast<EColor>(kGreen-6), static_cast<EColor>(kGreen-8)};

    std::vector<ROOT::RDF::RResultHandle> Nodes;

    int count = 0;
    for (auto& i: input_filenames)
    {
        SchottDataFrame df(MakeRDF(i, 8));

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

            return photon_id_eff;
            // *photon_iso_eff*photon_trg_eff; //Returns a vector of efficiencies

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

        auto photon_passes_cuts = ptCut
        .Define("photons_pass_cut_indices",
        [&](const RVec<Photon>& photons)
        {
            RVec<int> photon_indices;
            photon_indices.reserve(photons.size());

            for (int i = 0; i < photons.size(); i++)
            {
                if (
                (std::abs(photons[i].photon_eta) >= 2.37) or
                (photons[i].photon_pt <= 10e3) or
                (std::abs(photons[i].photon_eta) > 1.37 and std::abs(photons[i].photon_eta) < 1.52) or
                (not photons[i].photon_id_loose)
                )
                {
                    continue;
                }
                photon_indices.push_back(i);
            }

            return photon_indices;
        }, {"photons"})
        .Define("photons_pass_cuts",
        [&](RVec<Photon>& photons, RVec<int>& photon_indices)
        {
            return Take(photons, photon_indices);
        }, {"photons", "photons_pass_cut_indices"});

        //New dataframe node: contains only events from `photon_passes_cuts` that fail the resolved category and pass the merged category
        auto merged_reco_photons_matched = photon_passes_cuts.Filter(
        [&](RVec<Photon>& reco_photons_matched)
        {
           if (reco_photons_matched.size() == 1) // 1 photon in the event
           {
               return reco_photons_matched[0].photon_pt > 20e3; //event passes if photon pt > 20 GeV
           }
           else if (reco_photons_matched.empty())
           {
               return false; //fails if no photons in event
           }

           auto combs = Combinations(reco_photons_matched, 2); //all combinations of 2 reco-photons
           size_t length = combs[0].size(); //number of combinations
           double delta_r, m, pt, X, best_X, pt1, pt2, chosen_delta_r;

           for (size_t i=0; i<length; i++) //looping through all of the possible combinations of photons in each event
           {
               delta_r = DeltaR(reco_photons_matched[combs[0][i]].Vector(), reco_photons_matched[combs[1][i]].Vector());
               m = (reco_photons_matched[combs[0][i]].Vector() + reco_photons_matched[combs[1][i]].Vector()).M();
               pt = (reco_photons_matched[combs[0][i]].Vector() + reco_photons_matched[combs[1][i]].Vector()).Pt();
               X = delta_r*(pt/(2.0*m));
               //if it's the first combination or if new X is closer to 1
               //than current best_X and ΔR
               //between the two reco-photons < 1.5, then update best_X, pt1, pt2,
               //and the corresponding reco-photon indices x
               if (i==0 || ((std::abs(1-X) < std::abs(1-best_X)) and (delta_r < 1.5)))
               {
                   best_X = X;
                   pt1 = reco_photons_matched[combs[0][i]].photon_pt;
                   pt2 = reco_photons_matched[combs[1][i]].photon_pt;
                   chosen_delta_r = delta_r;
               }
           }
           if (pt1 > 10e3 && pt2 > 10e3 && best_X > 0.96 && best_X < 1.2 && chosen_delta_r < 1.5) //two photons corresponding to best_X must both have p_T > 10 GeV, ΔR < 1.5, and 0.96 < best_X < 1.2 for resolved... but this is merged, so if it passes resolved it fails merged
           {
               return false;
           }
           //if we get to this point, it means we've failed resolved
           for (auto& p: reco_photons_matched) //merged means at least 1 reco photon must have pt > 20 GeV
           {
               if (p.photon_pt > 20e3)
               {
                   return true; //passed merged if there's a reco-photon with pt > 20 GeV
               }
           }
           return false; //failed merged

        }, {"photons_pass_cuts"})
        .Define("merged_photon_index", //new column: consists of the index corresponding to the photon that made the event be classified as merged
        [&](const RVec<Photon>& rpm) //rpm = reco photons matched
        {
            for (auto i = 0; i < rpm.size(); i++)
            {
                if (rpm[i].photon_pt > 20e3)
                {
                    return i; //returning the index of the first photon that has photon_pt > 20 GeV
                }
            }
            return 0; //jic the compiler complains, should not come to this

        }, {"photons_pass_cuts"})
        .Define("merged_photon", //new column: The reco-photon corresponding to `merged_photon_index`
        [&](const RVec<Photon>& reco_photons_matched, int merged_photon_index)
        {
            return reco_photons_matched[merged_photon_index];

        }, {"photons_pass_cuts", "merged_photon_index"});

        auto dilepton_and_photon = merged_reco_photons_matched
        .Define("dilepton_mass",[&](RVec<Electron>& di_electrons)
        {
            auto four_momentum = di_electrons[0].Vector() + di_electrons[1].Vector();

            return four_momentum.M()/1e3;

        }, {"di_electrons"})
        .Define("merged_photon_pt",[&](Photon& merged_photon)
        {
            return merged_photon.Vector().Pt()/1e3;

        }, {"merged_photon"})
        .Define("reconstructed_mass",[&](RVec<Electron>& di_electrons, Photon& merged_photon)
        {
             auto four_momentum = di_electrons[0].Vector() + di_electrons[1].Vector();

             return (four_momentum + merged_photon.Vector()).M()/1e3;

        }, {"di_electrons", "merged_photon"})
        .Filter([](RVec<float>& Eratio, double reconstructed_mass)
        {
            return ((!Any(Eratio <= 0.8)) && ((reconstructed_mass <= 110) || (reconstructed_mass >= 130)));
        }, {"photon_shower_shape_e_ratio", "reconstructed_mass"})
        .Define("totEventWeight", [](RVec<float> totEventWeightVec, RVec<int>& photons_pass_cut_indices, int mpi)
        {
            return Take(totEventWeightVec, photons_pass_cut_indices)[mpi];

        }, {"totEventWeightVec", "photons_pass_cut_indices", "merged_photon_index"});

        if (count <= 2 || count >= 4) //background
        {
            Nodes.push_back(dilepton_and_photon.Sum<float>("totEventWeight"));
            Nodes.push_back(EventWeight.Sum<float>("EventWeight"));
        }

        //background and data
        Nodes.push_back(dilepton_and_photon.Histo1D<double>({prefixes[count], prefixes[count], 60u, 60, 120}, "dilepton_mass", "totEventWeight"));
        Nodes.push_back(dilepton_and_photon.Histo1D<double>({prefixes[count], prefixes[count++], 60u, 0, 165}, "merged_photon_pt", "totEventWeight"));
    }

//    0   1   2   3    //Z-gamma
//    4   5   6   7    //Z-gamma
//    8   9   10  11   //Z-gamma
//            12  13   //data
//    14  15  16  17   //Z-jets
//    18  19  20  21   //Z-jets
//    22  23  24  25   //Z-jets
//    26  27  28  29   //Z-jets
//    30  31  32  33   //Z-jets
//    34  35  36  37   //Z-jets
//    38  39  40  41   //Z-jets
//    42  43  44  45   //Z-jets
//    46  47  48  49   //Z-jets

    ROOT::RDF::RunGraphs(Nodes); // running all computation nodes concurrently

    count = 0;
    TCanvas* c1 = new TCanvas();
    TLegend* legend = new TLegend(0.15, 0.275, 0.475, 0.675);
    auto hs = new THStack("hs1","");
    double factor = 0;
    int back_count = 0;
    //Z-gamma
    for (auto& i: {0,4,8})
    {
        factor += (*Nodes[i].GetResultPtr<float>())*(SFs[back_count++]/ *Nodes[i+1].GetResultPtr<float>());
    }
    //Z-jets
    for (int i = 14, j = 0; (i <= 46 && j <= 8); i+=4, j++)
    {
        factor += (*Nodes[i].GetResultPtr<float>())*(JetNumeratorSFs[j]/ *Nodes[i+1].GetResultPtr<float>());
    }
    //Z-gamma
    
    Double_t data_integral = Nodes[12].GetResultPtr<TH1D>()->Integral();
    for (auto& i: {2,6,10,12})
    {
        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0 && i != 12)
        {
            Nodes[i].GetResultPtr<TH1D>()->Scale(SFs[count]/ *Nodes[i-1].GetResultPtr<float>());
            Nodes[i].GetResultPtr<TH1D>()->Scale(data_integral / factor);
        }

        if (i != 12)
        {
            Nodes[i].GetResultPtr<TH1D>()->SetLineColor(colors[count]);
            Nodes[i].GetResultPtr<TH1D>()->SetFillColor(colors[count++]);
            legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "f");
            hs->Add(&*Nodes[i].GetResultPtr<TH1D>());
        }
        else
        {
            Nodes[i].GetResultPtr<TH1D>()->SetLineColor(colors[count++]);
            legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "l");
        }
    }

    //Z-jets
    for (int i = 16, j = 0; (i <= 48 && j <= 8); i+=4, j++)
    {
        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0)
        {
            Nodes[i].GetResultPtr<TH1D>()->Scale(JetNumeratorSFs[j]/ *Nodes[i-1].GetResultPtr<float>());
            Nodes[i].GetResultPtr<TH1D>()->Scale(data_integral / factor);
        }
        Nodes[i].GetResultPtr<TH1D>()->SetLineColor(Jetscolors[j]);
        Nodes[i].GetResultPtr<TH1D>()->SetFillColor(Jetscolors[j]);
        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "f");
        hs->Add(&*Nodes[i].GetResultPtr<TH1D>());
    }

    hs->Draw("HIST");
    Nodes[12].GetResultPtr<TH1D>()->Draw("E1same");
    hs->SetTitle(";m_{ll} [GeV];Events");
    hs->GetYaxis()->CenterTitle(true);
    hs->GetXaxis()->SetTitleOffset(1.2);
    hs->GetYaxis()->SetTitleOffset(1.4);
    gStyle->SetOptStat(0);
    TLatex Tl;
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.6, 0.8, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.6, 0.7,"#sqrt{s} = 13 TeV");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("Fig48A.pdf");
    count = 0;
    hs = new THStack("hs2","");
    c1 = new TCanvas();
    legend = new TLegend(0.55, 0.2, 0.85, 0.6);
    //Z-gamma
    data_integral = Nodes[13].GetResultPtr<TH1D>()->Integral();
    for (auto& i: {3,7,11,13})
    {
        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0 && i != 13)
        {
            Nodes[i].GetResultPtr<TH1D>()->Scale(SFs[count]/ *Nodes[i-2].GetResultPtr<float>());
            Nodes[i].GetResultPtr<TH1D>()->Scale(data_integral / factor);
        }

        if (i != 13)
        {
            Nodes[i].GetResultPtr<TH1D>()->SetLineColor(colors[count]);
            Nodes[i].GetResultPtr<TH1D>()->SetFillColor(colors[count++]);
            legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "f");
            hs->Add(&*Nodes[i].GetResultPtr<TH1D>());
        }
        else
        {
            Nodes[i].GetResultPtr<TH1D>()->SetLineColor(colors[count++]);
            legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "l");
        }
    }
    //Z-Jets
    for (int i = 17, j = 0; (i <= 49 && j <= 8); i+=4, j++)
    {
        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0)
        {
            Nodes[i].GetResultPtr<TH1D>()->Scale(JetNumeratorSFs[j]/ *Nodes[i-2].GetResultPtr<float>());
            Nodes[i].GetResultPtr<TH1D>()->Scale(data_integral / factor);
        }
        Nodes[i].GetResultPtr<TH1D>()->SetLineColor(Jetscolors[j]);
        Nodes[i].GetResultPtr<TH1D>()->SetFillColor(Jetscolors[j]);
        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "f");
        hs->Add(&*Nodes[i].GetResultPtr<TH1D>());
    }

    hs->Draw("HIST");
    Nodes[13].GetResultPtr<TH1D>()->Draw("E1same");
    hs->SetTitle(";photon p_{T} [GeV];Events");
    hs->GetYaxis()->CenterTitle(true);
    hs->GetXaxis()->SetTitleOffset(1.2);
    hs->GetYaxis()->SetTitleOffset(1.4);
    gStyle->SetOptStat(0);
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.6, 0.8, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.6, 0.7,"#sqrt{s} = 13 TeV");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("Fig48B.pdf");
}

void fig59()
{
    std::vector<std::vector<std::string>> input_filenames = {
        //Signal
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root"}, //1 GeV
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600752.PhPy8EG_AZNLO_ggH125_mA3p0_v1.root"}, //2 GeV
        //Z gamma background
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617070._000001.LGNTuple.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617064._000001.LGNTuple.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617074._000001.LGNTuple.root"},
        //Data
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_data_test.root"},
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

    std::array<double,3> SFs = {((139e15)*(.871e-12)),((139e15)*(.199e-12)), ((139e15)*(.0345e-15))}; //numerators for Z-gamma bkg

    std::array<double,9> JetNumeratorSFs = {((139e15)*(1.9828e-9)*(0.821204)),((139e15)*(110.64e-12)*(0.69275)),((139e15)*(40.645e-12)*(0.615906)),((139e15)*(1.9817e-9)*(0.1136684)),((139e15)*(110.47e-12)*(0.1912956)),((139e15)*(40.674e-12)*(0.2326772)),((139e15)*(1.9819e-9)*(0.0656969)),((139e15)*(110.53e-12)*(0.1158741)),((139e15)*(40.68e-12)*(0.1535215))}; //numerators for jet bkg

    std::vector<const char*> prefixes = {"Sig m_{A} = 1 GeV", "Sig m_{A} = 2 GeV", "pty2_9_17", "pty_17_myy_0_80", "pty_17_myy_80", "data", "Zee_lightJet_0-70", "Zee_lightJet_70-140", "Zee_lightJet_140-280", "Zee_cJet_0-70", "Zee_cJet_70-140", "Zee_cJet_140-280", "Zee_bJet_0-70", "Zee_bJet_70-140", "Zee_bJet_140-280"};
    std::vector<EColor> colors = {kMagenta, kCyan, static_cast<EColor>(kRed-7), static_cast<EColor>(kRed-9), static_cast<EColor>(kRed-10), kBlack};
    std::vector<EColor> Jetscolors = {static_cast<EColor>(kGreen+0), static_cast<EColor>(kGreen-4), static_cast<EColor>(kGreen-7), static_cast<EColor>(kGreen-9), static_cast<EColor>(kGreen-10), static_cast<EColor>(kGreen+1), static_cast<EColor>(kGreen-3), static_cast<EColor>(kGreen-6), static_cast<EColor>(kGreen-8)};

    std::vector<ROOT::RDF::RResultHandle> Nodes;
    int count = 0;
    for (auto& i: input_filenames)
    {
        SchottDataFrame df(MakeRDF(i, 8));

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

            return photon_id_eff;
            // *photon_iso_eff*photon_trg_eff; //Returns a vector of efficiencies

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

        auto photon_passes_cuts = ptCut
        .Define("photons_pass_cut_indices",
        [&](const RVec<Photon>& photons)
        {
            RVec<int> photon_indices;
            photon_indices.reserve(photons.size());

            for (int i = 0; i < photons.size(); i++)
            {
                if (
                (std::abs(photons[i].photon_eta) >= 2.37) or
                (photons[i].photon_pt <= 10e3) or
                (std::abs(photons[i].photon_eta) > 1.37 and std::abs(photons[i].photon_eta) < 1.52) or
                (not photons[i].photon_id_loose)
                )
                {
                    continue;
                }
                photon_indices.push_back(i);
            }

            return photon_indices;
        }, {"photons"})
        .Define("photons_pass_cuts",
        [&](RVec<Photon>& photons, RVec<int>& photon_indices)
        {
            return Take(photons, photon_indices);
        }, {"photons", "photons_pass_cut_indices"});

        //New dataframe node: contains only events from `photon_passes_cuts` that fail the resolved category and pass the merged category
        auto merged_reco_photons_matched = photon_passes_cuts.Filter(
        [&](RVec<Photon>& reco_photons_matched)
        {
           if (reco_photons_matched.size() == 1) // 1 photon in the event
           {
               return reco_photons_matched[0].photon_pt > 20e3; //event passes if photon pt > 20 GeV
           }
           else if (reco_photons_matched.empty())
           {
               return false; //fails if no photons in event
           }

           auto combs = Combinations(reco_photons_matched, 2); //all combinations of 2 reco-photons
           size_t length = combs[0].size(); //number of combinations
           double delta_r, m, pt, X, best_X, pt1, pt2, chosen_delta_r;

           for (size_t i=0; i<length; i++) //looping through all of the possible combinations of photons in each event
           {
               delta_r = DeltaR(reco_photons_matched[combs[0][i]].Vector(), reco_photons_matched[combs[1][i]].Vector());
               m = (reco_photons_matched[combs[0][i]].Vector() + reco_photons_matched[combs[1][i]].Vector()).M();
               pt = (reco_photons_matched[combs[0][i]].Vector() + reco_photons_matched[combs[1][i]].Vector()).Pt();
               X = delta_r*(pt/(2.0*m));
               //if it's the first combination or if new X is closer to 1
               //than current best_X and ΔR
               //between the two reco-photons < 1.5, then update best_X, pt1, pt2,
               //and the corresponding reco-photon indices x
               if (i==0 || ((std::abs(1-X) < std::abs(1-best_X)) and (delta_r < 1.5)))
               {
                   best_X = X;
                   pt1 = reco_photons_matched[combs[0][i]].photon_pt;
                   pt2 = reco_photons_matched[combs[1][i]].photon_pt;
                   chosen_delta_r = delta_r;
               }
           }
           if (pt1 > 10e3 && pt2 > 10e3 && best_X > 0.96 && best_X < 1.2 && chosen_delta_r < 1.5) //two photons corresponding to best_X must both have p_T > 10 GeV, ΔR < 1.5, and 0.96 < best_X < 1.2 for resolved... but this is merged, so if it passes resolved it fails merged
           {
               return false;
           }
           //if we get to this point, it means we've failed resolved
           for (auto& p: reco_photons_matched) //merged means at least 1 reco photon must have pt > 20 GeV
           {
               if (p.photon_pt > 20e3)
               {
                   return true; //passed merged if there's a reco-photon with pt > 20 GeV
               }
           }
           return false; //failed merged

        }, {"photons_pass_cuts"})
        .Define("merged_photon_index", //new column: consists of the index corresponding to the photon that made the event be classified as merged
        [&](const RVec<Photon>& rpm) //rpm = reco photons matched
        {
            for (auto i = 0; i < rpm.size(); i++)
            {
                if (rpm[i].photon_pt > 20e3)
                {
                    return i; //returning the index of the first photon that has photon_pt > 20 GeV
                }
            }
            return 0; //jic the compiler complains, should not come to this

        }, {"photons_pass_cuts"})
        .Define("merged_photon", //new column: The reco-photon corresponding to `merged_photon_index`
        [&](const RVec<Photon>& reco_photons_matched, int merged_photon_index)
        {
            return reco_photons_matched[merged_photon_index];

        }, {"photons_pass_cuts", "merged_photon_index"});

        auto dilepton_and_photon = merged_reco_photons_matched
        .Define("dilepton_mergedPhoton_deltaR",[&](PtEtaPhiEVector& dilep, Photon& merged_photon)
        {
            return DeltaR(dilep, merged_photon.Vector());

        }, {"dilep", "merged_photon"})
        .Define("reconstructed_mass",[&](PtEtaPhiEVector& dilep, Photon& merged_photon)
        {
            return (dilep + merged_photon.Vector()).M()/1e3;

        }, {"dilep", "merged_photon"})
        .Define("totEventWeight", [](RVec<float> totEventWeightVec, RVec<int>& photons_pass_cut_indices, int mpi)
        {
            return Take(totEventWeightVec, photons_pass_cut_indices)[mpi];

        }, {"totEventWeightVec", "photons_pass_cut_indices", "merged_photon_index"});

        auto SB = dilepton_and_photon.Filter(
        [](RVec<float>& Eratio, double reconstructed_mass)
        {
            return ((!Any(Eratio <= 0.8)) && ((reconstructed_mass <= 110) || (reconstructed_mass >= 130)));
        }, {"photon_shower_shape_e_ratio", "reconstructed_mass"});

        auto SR = dilepton_and_photon.Filter(
        [](RVec<float>& Eratio, double reconstructed_mass)
        {
            return ((!Any(Eratio <= 0.8)) && ((reconstructed_mass > 110) && (reconstructed_mass < 130)));
        }, {"photon_shower_shape_e_ratio", "reconstructed_mass"});

        if ((count >= 2 && count <= 4) || (count >= 6)) //background: Z-gamma and Z-jets
        {
            Nodes.push_back(SB.Sum<float>("totEventWeight"));
            Nodes.push_back(SR.Sum<float>("totEventWeight"));
            Nodes.push_back(EventWeight.Sum<float>("EventWeight"));
        }

        Nodes.push_back(SB.Histo1D<double>({prefixes[count], prefixes[count], 100u, 0, 4.7}, "dilepton_mergedPhoton_deltaR", "totEventWeight"));
        Nodes.push_back(SR.Histo1D<double>({prefixes[count], prefixes[count++], 100u, 0, 4.7}, "dilepton_mergedPhoton_deltaR", "totEventWeight"));
    }

//                0   1       //ma5
//                2   3       //ma1
//    4   5   6   7   8       //Z-gamma
//    9   10  11  12  13      //Z-gamma
//    14  15  16  17  18      //Z-gamma
//                19  20      //data
//    21  22  23  24  25      //Z-jets
//    26  27  28  29  30      //Z-jets
//    31  32  33  34  35      //Z-jets
//    36  37  38  39  40      //Z-jets
//    41  42  43  44  45      //Z-jets
//    46  47  48  49  50      //Z-jets
//    51  52  53  54  55      //Z-jets
//    56  57  58  59  60      //Z-jets
//    61  62  63  64  65      //Z-jets

    ROOT::RDF::RunGraphs(Nodes); // running all computation nodes concurrently

    count = 0;
    TCanvas* c1 = new TCanvas();
    TLegend* legend = new TLegend(0.725, 0.5, 0.875, 0.7);
    double factor;
    for (auto& i: {0,2}) //sideband signal
    {
        Nodes[i].GetResultPtr<TH1D>()->SetLineColor(colors[count++]);
        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "l");

        if (i == 0)
        {
            Nodes[i].GetResultPtr<TH1D>()->SetTitle(";#DeltaR (ll#gamma);Events");
            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);
            Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->SetTitleOffset(1.2);
            Nodes[i].GetResultPtr<TH1D>()->SetAxisRange(0., 6, "Y");
            Nodes[i].GetResultPtr<TH1D>()->Draw("HIST");
        }
        else
        {
            Nodes[i].GetResultPtr<TH1D>()->Draw("HISTsame");
            gPad->Modified(); gPad->Update();
        }
    }

    gStyle->SetOptStat(0);
    TLatex Tl;
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.6, 0.85, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.6, 0.76,"#sqrt{s} = 13 TeV");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("Fig59D.pdf");

    factor = 0;
    count = 0;
    for (auto& i: {4,9,14}) //Z-gamma sideband region
    {
        factor += (*Nodes[i].GetResultPtr<float>())*(SFs[count++]/ *Nodes[i+2].GetResultPtr<float>());
    }

    for (int i = 21, j = 0; ( i <= 61 && j <= 8); i+=5, j++) //Z-jets sideband region
    {
        factor += (*Nodes[i].GetResultPtr<float>())*(JetNumeratorSFs[j]/ *Nodes[i+2].GetResultPtr<float>());
    }

    count = 0;
    auto hs = new THStack("hs1","");
    c1 = new TCanvas();
    legend = new TLegend(0.71, 0.435, 0.88, 0.785);
    //Z-gamma
    
    Double_t data_integral = Nodes[19].GetResultPtr<TH1D>()->Integral();
    for (auto& i: {7,12,17,19}) //Z-gamma & data sideband region
    {
        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0 && i != 19)
        {
            Nodes[i].GetResultPtr<TH1D>()->Scale(SFs[count]/ *Nodes[i-1].GetResultPtr<float>());
            Nodes[i].GetResultPtr<TH1D>()->Scale(data_integral / factor);
        }

        if (i != 19)
        {
            Nodes[i].GetResultPtr<TH1D>()->SetLineColor(colors[2+count]);
            Nodes[i].GetResultPtr<TH1D>()->SetFillColor(colors[2+count++]);
            legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "f");
            hs->Add(&*Nodes[i].GetResultPtr<TH1D>());
        }
        else
        {
            Nodes[i].GetResultPtr<TH1D>()->SetLineColor(colors[2+count]);
            legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "l");
        }
    }

    for (int i = 24, j = 0; (i <= 64 && j <= 8); i+=5, j++) //Z-jets sideband region
    {
        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0)
        {
            Nodes[i].GetResultPtr<TH1D>()->Scale(JetNumeratorSFs[j]/ *Nodes[i-1].GetResultPtr<float>());
            Nodes[i].GetResultPtr<TH1D>()->Scale(data_integral / factor);
        }
        Nodes[i].GetResultPtr<TH1D>()->SetLineColor(Jetscolors[j]);
        Nodes[i].GetResultPtr<TH1D>()->SetFillColor(Jetscolors[j]);
        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "f");
        hs->Add(&*Nodes[i].GetResultPtr<TH1D>());
    }

    hs->Draw("HIST");
    Nodes[19].GetResultPtr<TH1D>()->Draw("E1same");
    hs->SetTitle(";#DeltaR (ll#gamma);Events");
//    hs->SetMaximum(11500);
    hs->GetYaxis()->CenterTitle(true);
    hs->GetXaxis()->SetTitleOffset(1.2);
    hs->GetYaxis()->SetTitleOffset(1.35);

    gStyle->SetOptStat(0);
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.2, 0.82, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.2, 0.72,"#sqrt{s} = 13 TeV");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("Fig59A.pdf");

    factor = 0;
    count = 0;
    double sig_factor = Nodes[1].GetResultPtr<TH1D>()->Integral();
    std::vector<double> backScalings;

    for (auto& i: {5,10,15}) //Z-gamma SR
    {
        backScalings.push_back((*Nodes[i].GetResultPtr<float>())*(SFs[count]/ *Nodes[i+1].GetResultPtr<float>()));

        factor += (*Nodes[i].GetResultPtr<float>())*(SFs[count++]/ *Nodes[i+1].GetResultPtr<float>());
    }

    for (int i = 22, j = 0; ( i <= 62 && j <= 8); i+=5, j++) //Z-jets SR
    {
        backScalings.push_back((*Nodes[i].GetResultPtr<float>())*(JetNumeratorSFs[j]/ *Nodes[i+1].GetResultPtr<float>()));

        factor += (*Nodes[i].GetResultPtr<float>())*(JetNumeratorSFs[j]/ *Nodes[i+1].GetResultPtr<float>());
    }

    for (double& i: backScalings)
    {
        i = (i/factor)*sig_factor;
    }

    count=2;
    c1 = new TCanvas();
    legend = new TLegend(0.67, 0.5, 0.86, 0.8);
    hs = new THStack("hs2","");
    //Z-gamma
    
    data_integral = Nodes[20].GetResultPtr<TH1D>()->Integral();
    for (auto& i: {8,13,18,20}) //Z-gamma & data signal region
    {
        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0 && i != 20)
        {
            Nodes[i].GetResultPtr<TH1D>()->Scale(SFs[count-2]/ *Nodes[i-2].GetResultPtr<float>());
//            Nodes[i].GetResultPtr<TH1D>()->Scale(data_integral / factor);
        }
        if (i != 20)
        {
            hs->Add(static_cast<TH1D*>(Nodes[i].GetResultPtr<TH1D>()->Clone()));
            Nodes[i].GetResultPtr<TH1D>()->SetLineColor(colors[count]);
            Nodes[i].GetResultPtr<TH1D>()->SetFillColor(colors[count++]);
            legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "f");
        }
//        else
//        {
//            Nodes[i].GetResultPtr<TH1D>()->SetLineColor(colors[count++]);
//            legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "l");
//        }
    }

    for (int i = 25, j = 0; (i <= 65 && j <= 8); i+=5, j++) //Z-jets SR
    {
        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0)
        {
            Nodes[i].GetResultPtr<TH1D>()->Scale(JetNumeratorSFs[j]/ *Nodes[i-2].GetResultPtr<float>());
//            Nodes[i].GetResultPtr<TH1D>()->Scale(data_integral / factor);
        }
        Nodes[i].GetResultPtr<TH1D>()->SetLineColor(Jetscolors[j]);
        Nodes[i].GetResultPtr<TH1D>()->SetFillColor(Jetscolors[j]);
        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "f");
        hs->Add(static_cast<TH1D*>(Nodes[i].GetResultPtr<TH1D>()->Clone()));
    }

    hs->Draw("HIST");
    count=0;
    for (auto& i: {1,3}) // signal SR
    {
        Nodes[i].GetResultPtr<TH1D>()->SetLineColor(colors[count]);
        Nodes[i].GetResultPtr<TH1D>()->Scale(data_integral / Nodes[i].GetResultPtr<TH1D>()->Integral());
        Nodes[i].GetResultPtr<TH1D>()->SetLineWidth(2);
        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), prefixes[count++], "l");

        Nodes[i].GetResultPtr<TH1D>()->SetTitle(";#DeltaR (ll#gamma);Events");
        Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);
        Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->SetTitleOffset(1.2);
        Nodes[i].GetResultPtr<TH1D>()->SetAxisRange(0., 68, "Y");
        Nodes[i].GetResultPtr<TH1D>()->Draw("HISTsame");
        gPad->Modified(); gPad->Update();
    }
    hs->SetMaximum(2.05e4);
    hs->SetTitle(";#DeltaR (ll#gamma);Events");
    hs->GetYaxis()->CenterTitle(true);
    hs->GetXaxis()->SetTitleOffset(1.2);
    hs->GetYaxis()->SetTitleOffset(1.35);

    gStyle->SetOptStat(0);
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.37, 0.85, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.37, 0.75,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("Fig59B.pdf");

    count=2;
    c1 = new TCanvas();
    legend = new TLegend(0.7, 0.35, 0.875, 0.75);
    count=2;
    hs = new THStack("hs3","");
    int temp = 0;
    for (auto& i: {8,13,18,20}) //background and data SR
    {
        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0 && i != 20)
        {
            Nodes[i].GetResultPtr<TH1D>()->Scale(backScalings[temp++]/Nodes[i].GetResultPtr<TH1D>()->Integral());
        }

        if (i != 20)
        {
            hs->Add(&*Nodes[i].GetResultPtr<TH1D>());
            Nodes[i].GetResultPtr<TH1D>()->SetFillColor(colors[count++]);
            legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "f");
        }
    }
    temp = 3;
    for (int i = 25, j = 0; (i <= 65 && j <= 8); i+=5, j++) //Z-jets SR
    {
        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0)
        {
            Nodes[i].GetResultPtr<TH1D>()->Scale(backScalings[temp++]/Nodes[i].GetResultPtr<TH1D>()->Integral());
        }
        Nodes[i].GetResultPtr<TH1D>()->SetFillColor(Jetscolors[j]);
        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "f");

        hs->Add(&*Nodes[i].GetResultPtr<TH1D>());
    }

    hs->Draw("HIST");
    count=0;
    for (auto& i: {1,3}) //signal SR
    {
        Nodes[i].GetResultPtr<TH1D>()->SetLineColor(colors[count++]);
        Nodes[i].GetResultPtr<TH1D>()->SetLineWidth(2);
        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), prefixes[count-1], "l");

        Nodes[i].GetResultPtr<TH1D>()->Scale(sig_factor/Nodes[i].GetResultPtr<TH1D>()->Integral());
        Nodes[i].GetResultPtr<TH1D>()->SetTitle(";#DeltaR (ll#gamma);Events");
        Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);
        Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->SetTitleOffset(1.2);
        Nodes[i].GetResultPtr<TH1D>()->Draw("HISTsame");
        gPad->Modified(); gPad->Update();
    }
    hs->SetMaximum(120);
    hs->SetTitle(";#DeltaR (ll#gamma);Events");
    hs->GetYaxis()->CenterTitle(true);
    hs->GetXaxis()->SetTitleOffset(1.2);

    gStyle->SetOptStat(0);
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.2, 0.825, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.2, 0.725,"#sqrt{s} = 13 TeV");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("Fig59C.pdf");
}
//
//void Table9()
//{
//    std::vector<std::vector<std::string>> input_filenames =
//    {
//        //Z gamma background
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617070._000001.LGNTuple.root"}, //pty2_9_17
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617064._000001.LGNTuple.root"}, //pty_17_myy_0_80
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617074._000001.LGNTuple.root"}, //pty_17_myy_80
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
//    std::array<double,9> JetNumeratorSFs = {((139e15)*(1.9828e-9)*(0.821204)),((139e15)*(110.64e-12)*(0.69275)),((139e15)*(40.645e-12)*(0.615906)),((139e15)*(1.9817e-9)*(0.1136684)),((139e15)*(110.47e-12)*(0.1912956)),((139e15)*(40.674e-12)*(0.2326772)),((139e15)*(1.9819e-9)*(0.0656969)),((139e15)*(110.53e-12)*(0.1158741)),((139e15)*(40.68e-12)*(0.1535215))}; //numerators for jet bkg
//    std::array<double,3> SFs = {((139e15)*(.871e-12)),((139e15)*(.199e-12)), ((139e15)*(.0345e-15))}; //numerators for Z-gamma bkg
//
//    std::vector<ROOT::RDF::RResultHandle> Nodes;
//
//    for (auto& i: input_filenames)
//    {
//        SchottDataFrame df(MakeRDF(i, 8));
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
//        [&](RVec<TruthParticle> truth_particles)
//        {
//            truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
//            [](TruthParticle& x)
//            {
//                return ((std::abs(x.mc_pdg_id) != 22) || (std::abs(x.mc_eta) >= 2.37));
//
//            }), truth_particles.end());
//
//            return truth_particles;
//
//        }, {"truth_particles"});
//
//        auto diphotons = photon_passes_cuts.Define("chosen_two",
//        [](RVec<TruthParticle>& reco_photons_matched)
//        {
//            RVec<TruthParticle> x;
//            if (reco_photons_matched.size() < 2)
//            {
//                return x;
//            }
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
//                    pt1 = reco_photons_matched[combs[0][i]].mc_pt;
//                    pt2 = reco_photons_matched[combs[1][i]].mc_pt;
//                    chosen_delta_r = delta_r;
//                    x = {reco_photons_matched[combs[0][i]], reco_photons_matched[combs[1][i]]};
//                }
//            }
//            if (pt1 > 10e3 && pt2 > 10e3 && best_X > 0.96 && best_X < 1.2 && chosen_delta_r < 1.5)
//            {
//                return x;
//            }
//            x.clear(); //it doesn't pass, so return an empty vector
//            return x;
//        }, {"photons_pass_cuts"}).Filter(
//        [&](RVec<TruthParticle>& reco_photons_matched)
//        {
//            return reco_photons_matched.size() == 2;
//        }, {"chosen_two"}).Define("leading_photon",
//        [&](RVec<TruthParticle>& reco_photons_matched)
//        {
//            return (reco_photons_matched[0].mc_pt > reco_photons_matched[1].mc_pt)
//                    ? reco_photons_matched[0]
//                    : reco_photons_matched[1];
//        },{"chosen_two"}).Define("subleading_photon",
//        [&](RVec<TruthParticle>& reco_photons_matched)
//        {
//            return (reco_photons_matched[0].mc_pt < reco_photons_matched[1].mc_pt)
//                    ? reco_photons_matched[0]
//                    : reco_photons_matched[1];
//        },{"chosen_two"}).Define("leading_photon_pdg_id_origin",
//        [](TruthParticle& leading_photon, RVec<TruthParticle>& truth_particles)
//        {
//            TruthParticle origin = leading_photon;
//            int origin_id = 0;
//            bool found;
//
//            do
//            {
//                found = false;
//                for (auto& i: truth_particles)
//                {
//                    if (origin.mc_parent_barcode == i.mc_barcode)
//                    {
//                        origin = i;
//                        origin_id = i.mc_pdg_id;
//                        found = true;
//                    }
//                }
//            } while (found);
//
//            return origin_id;
//        }, {"leading_photon", "truth_particles"}).Define("subleading_photon_pdg_id_origin",
//        [](TruthParticle& subleading_photon, RVec<TruthParticle>& truth_particles)
//        {
//            TruthParticle origin = subleading_photon;
//            int origin_id = 0;
//            bool found;
//
//            do
//            {
//             found = false;
//             for (auto& i: truth_particles)
//             {
//                 if (origin.mc_parent_barcode == i.mc_barcode)
//                 {
//                     origin = i;
//                     origin_id = i.mc_pdg_id;
//                     found = true;
//                 }
//             }
//            } while (found);
//
//            return origin_id;
//        }, {"subleading_photon", "truth_particles"});
//
//        Nodes.push_back(diphotons.Take<int, RVec<int>>("leading_photon_pdg_id_origin"));
//        Nodes.push_back(diphotons.Take<int, RVec<int>>("subleading_photon_pdg_id_origin"));
//        Nodes.push_back(diphotons.Count());
//    }
//
//    ROOT::RDF::RunGraphs(Nodes); // running all computation nodes concurrently
//
//    std::unordered_map<int,int>
//    leading_id_freqs_Z_gamma, subleading_id_freqs_Z_gamma,
//    leading_id_freqs_Z_jets, subleading_id_freqs_Z_jets,
//    leading_id_freqs_total, subleading_id_freqs_total;
//
//    double total_Z_gamma = *Nodes[2].GetResultPtr<ULong64_t>();
//    double total_Z_jets = *Nodes[5].GetResultPtr<ULong64_t>();
//    double total = total_Z_gamma + total_Z_jets;
//
//    for (auto& i: *Nodes[0].GetResultPtr<RVec<int>>())
//    {
//        leading_id_freqs_Z_gamma[i]++;
//        leading_id_freqs_total[i]++;
//    }
//    for (auto& i: *Nodes[1].GetResultPtr<RVec<int>>())
//    {
//        subleading_id_freqs_Z_gamma[i]++;
//        subleading_id_freqs_total[i]++;
//    }
//    for (auto& i: *Nodes[3].GetResultPtr<RVec<int>>())
//    {
//        leading_id_freqs_Z_jets[i]++;
//        leading_id_freqs_total[i]++;
//    }
//    for (auto& i: *Nodes[4].GetResultPtr<RVec<int>>())
//    {
//        subleading_id_freqs_Z_jets[i]++;
//        subleading_id_freqs_total[i]++;
//    }
//
//    std::cout << R"--(\section*{Table 9}\flushleft)--" << '\n';
//        std::cout << R"--(\hspace{-3cm}\scalebox{0.65}{
//                          \setlength\extrarowheight{2pt}
//                          \renewcommand{\arraystretch}{1.5})--" << '\n';
//        std::cout << R"--(\begin{tabular}{|c|c|c|c|c|c|c|})--" << '\n';
//        std::cout << R"--(\hline)--" << '\n';
//        std::cout << R"--(\multicolumn{7}{|c|}{\parbox{\linewidth}{\centering Resolved Category \\ Background photon truth origin fractions}}\\[5 pt] \hline)--" << '\n';
//        std::cout << R"--( {} & \multicolumn{2}{|c|}{\parbox{7.5cm}{\centering $Z\gamma$ \\ (\%)}} & \multicolumn{2}{|c|}{\parbox{7.5cm}{\centering $Z$+jets \\ (\%)}} & \multicolumn{2}{|c|}{\parbox{7.5cm}{\centering Total bkg \\ (\%)}} \\[5 pt]
//            \hline)--" << '\n';
//        std::cout << R"--( {Pdg Id} & leading photon & subleading photon & leading photon & subleading photon & leading photon & subleading photon \\[5 pt]
//            \hline)--" << '\n';
//
//
//    for (auto& i: leading_id_freqs_total)
//    {
//        std::cout << i.first << " & " << std::setprecision(2) << std::fixed
//        << 100*(leading_id_freqs_Z_gamma[i.first]/total_Z_gamma)
//        << " & " << std::setprecision(2) << std::fixed
//        << 100*(subleading_id_freqs_Z_gamma[i.first]/total_Z_gamma)
//        << " & " << std::setprecision(2) << std::fixed
//        << 100*(leading_id_freqs_Z_jets[i.first]/total_Z_jets)
//        << " & " << std::setprecision(2) << std::fixed
//        << 100*(subleading_id_freqs_Z_jets[i.first]/total_Z_jets)
//        << " & " << std::setprecision(2) << std::fixed
//        << 100*(leading_id_freqs_total[i.first]/total)
//        << " & " << std::setprecision(2) << std::fixed
//        << 100*(subleading_id_freqs_total[i.first]/total)
//        << R"--( \\ \hline)--" << '\n';
//    }
//
//    std::cout << R"--(\end{tabular}})--" << '\n';
//
//    std::cout << "\n\n\n";
//}
//
//void Table10()
//{
//    std::vector<std::vector<std::string>> input_filenames =
//    {
//        //Z gamma background
//        {
//            "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617070._000001.LGNTuple.root",
//            "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617064._000001.LGNTuple.root",
//            "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617074._000001.LGNTuple.root"
//        },
//        //Jets
//        {
//            "/Users/edwardfinkelstein/ATLAS_axion/Jets/Zee_lightJet_0-70.root",
//            "/Users/edwardfinkelstein/ATLAS_axion/Jets/Zee_lightJet_70-140.root",
//            "/Users/edwardfinkelstein/ATLAS_axion/Jets/Zee_lightJet_140-280.root",
//            "/Users/edwardfinkelstein/ATLAS_axion/Jets/Zee_cJet_0-70.root",
//            "/Users/edwardfinkelstein/ATLAS_axion/Jets/Zee_cJet_70-140.root",
//            "/Users/edwardfinkelstein/ATLAS_axion/Jets/Zee_cJet_140-280.root",
//            "/Users/edwardfinkelstein/ATLAS_axion/Jets/Zee_bJet_0-70.root",
//            "/Users/edwardfinkelstein/ATLAS_axion/Jets/Zee_bJet_70-140.root",
//            "/Users/edwardfinkelstein/ATLAS_axion/Jets/Zee_bJet_140-280.root",
//        },
//    };
//
//    std::vector<ROOT::RDF::RResultHandle> Nodes;
//
//    for (auto& i: input_filenames)
//    {
//        SchottDataFrame df(MakeRDF(i, 8));
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
//        [&](RVec<TruthParticle> truth_particles)
//        {
//            truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
//            [](TruthParticle& x)
//            {
//                return ((std::abs(x.mc_pdg_id) != 22) || (std::abs(x.mc_eta) >= 2.37));
//
//            }), truth_particles.end());
//
//            return truth_particles;
//
//        }, {"truth_particles"});
//
//        auto diphotons = photon_passes_cuts.Define("chosen_two",
//        [](RVec<TruthParticle>& reco_photons_matched)
//        {
//            RVec<TruthParticle> x;
//            if (reco_photons_matched.size() < 2)
//            {
//                return x;
//            }
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
//                    pt1 = reco_photons_matched[combs[0][i]].mc_pt;
//                    pt2 = reco_photons_matched[combs[1][i]].mc_pt;
//                    chosen_delta_r = delta_r;
//                    x = {reco_photons_matched[combs[0][i]], reco_photons_matched[combs[1][i]]};
//                }
//            }
//            if (pt1 > 10e3 && pt2 > 10e3 && best_X > 0.96 && best_X < 1.2 && chosen_delta_r < 1.5)
//            {
//                return x;
//            }
//            x.clear();
//            return x; //it doesn't pass, so return an empty vector
//        }, {"photons_pass_cuts"}).Filter(
//        [&](RVec<TruthParticle>& reco_photons_matched)
//        {
//            return reco_photons_matched.size() == 2;
//        }, {"chosen_two"}).Define("leading_photon",
//        [&](RVec<TruthParticle>& reco_photons_matched)
//        {
//            return (reco_photons_matched[0].mc_pt > reco_photons_matched[1].mc_pt)
//                    ? reco_photons_matched[0]
//                    : reco_photons_matched[1];
//        },{"chosen_two"}).Define("subleading_photon",
//        [&](RVec<TruthParticle>& reco_photons_matched)
//        {
//            return (reco_photons_matched[0].mc_pt < reco_photons_matched[1].mc_pt)
//                    ? reco_photons_matched[0]
//                    : reco_photons_matched[1];
//        },{"chosen_two"}).Define("photon_pdg_id_origin",
//        [](TruthParticle& leading_photon, TruthParticle& subleading_photon, RVec<TruthParticle>& truth_particles)
//        {
//            TruthParticle leading_origin = leading_photon;
//            TruthParticle subleading_origin = subleading_photon;
//            int leading_origin_id = 0, subleading_origin_id = 0;
//            bool leading_found, subleading_found;
//
//            do
//            {
//                leading_found = false;
//                subleading_found = false;
//                for (auto& i: truth_particles)
//                {
//                    if (leading_origin.mc_parent_barcode == i.mc_barcode)
//                    {
//                        leading_origin = i;
//                        leading_origin_id = i.mc_pdg_id;
//                        leading_found = true;
//                    }
//
//                    if (subleading_origin.mc_parent_barcode == i.mc_barcode)
//                    {
//                        subleading_origin = i;
//                        subleading_origin_id = i.mc_pdg_id;
//                        subleading_found = true;
//                    }
//                }
//            } while (leading_found || subleading_found);
//
//            return std::to_string(leading_origin_id)+"/"+std::to_string(subleading_origin_id);
//        }, {"leading_photon", "subleading_photon", "truth_particles"});
//
//        Nodes.push_back(diphotons.Take<std::string, RVec<std::string>>("photon_pdg_id_origin"));
//        Nodes.push_back(diphotons.Count());
//
//    }
//
//    std::unordered_map<std::string,int>
//    id_freqs_Z_gamma, id_freqs_Z_jets, id_freqs_Z_total;
//
//    ROOT::RDF::RunGraphs(Nodes); // running all computation nodes concurrently
//
//    double total_Z_gamma = *Nodes[1].GetResultPtr<ULong64_t>();
//    double total_Z_jets = *Nodes[3].GetResultPtr<ULong64_t>();
//    double total = total_Z_gamma + total_Z_jets;
//
//    for (auto& i: *Nodes[0].GetResultPtr<RVec<std::string>>())
//    {
//        id_freqs_Z_gamma[i]++;
//        id_freqs_Z_total[i]++;
//    }
//
//    for (auto& i: *Nodes[2].GetResultPtr<RVec<std::string>>())
//    {
//        id_freqs_Z_jets[i]++;
//        id_freqs_Z_total[i]++;
//    }
//
//    std::cout << R"--(\section*{Table 10}\flushleft)--" << '\n';
//        std::cout << R"--(\hspace{-3cm}\scalebox{0.65}{
//                          \setlength\extrarowheight{2pt}
//                          \renewcommand{\arraystretch}{1.5})--" << '\n';
//        std::cout << R"--(\begin{tabular}{|c|c|c|c|})--" << '\n';
//        std::cout << R"--(\hline)--" << '\n';
//        std::cout << R"--(\multicolumn{4}{|c|}{\parbox{\linewidth}{\centering Resolved Category \\ Background photon truth origin fractions}}\\[5 pt] \hline)--" << '\n';
//        std::cout << R"--( \multicolumn{1}{|c|}{\parbox{7.5cm}{\centering
//            lead/sublead \; (Pdg Id) \\ photon pairs}} &
//            \multicolumn{1}{|c|}{\parbox{4.5cm}{\centering $Z\gamma$ \\ (\%)}} &
//            \multicolumn{1}{|c|}{\parbox{4.5cm}{\centering $Z$+jets \\ (\%)}} & \multicolumn{1}{|c|}{\parbox{4.5cm}{\centering Total bkg \\ (\%)}} \\[5 pt]
//            \hline)--" << '\n';
//
//    for (auto& i: id_freqs_Z_total)
//    {
//        std::cout << i.first << " & " << std::setprecision(2) << std::fixed
//        << 100*(id_freqs_Z_gamma[i.first]/total_Z_gamma)
//        << " & " << std::setprecision(2) << std::fixed
//        << 100*(id_freqs_Z_jets[i.first]/total_Z_jets)
//        << " & " << std::setprecision(2) << std::fixed
//        << 100*(id_freqs_Z_total[i.first]/total)
//        << R"--( \\ \hline)--" << '\n';
//    }
//    std::cout << R"--(\end{tabular}})--" << '\n';
//
//    std::cout << "\n\n\n";
//}
//
//void Table16()
//{
//    std::vector<std::vector<std::string>> input_filenames =
//    {
//        //Z-gamma
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617070._000001.LGNTuple.root"}, {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617064._000001.LGNTuple.root"},
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617074._000001.LGNTuple.root"},
//        //data
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_data_test.root"},
//        //signal
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root"}, //1 GeV
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600751.PhPy8EG_AZNLO_ggH125_mA2p0_v1.root"}, // 2 GeV
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600752.PhPy8EG_AZNLO_ggH125_mA3p0_v1.root"}, // 3 GeV
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root"}, //5 GeV
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600756.PhPy8EG_AZNLO_ggH125_mA9p0_v1.root"}, // 9 GeV
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
//    std::vector<std::string> prefixes = { R"--(pty2\_9\_17)--", R"--(pty\_17\_myy\_0\_80)--", R"--(pty\_17\_myy\_80)--", "data", R"--($\text{Sig } m_{A}$ = 1 GeV)--", R"--($\text{Sig } m_{A}$ = 2 GeV)--", R"--($\text{Sig } m_{A}$ = 3 GeV)--", R"--($\text{Sig } m_{A}$ = 5 GeV)--", R"--($\text{Sig } m_{A}$ = 9 GeV)--", R"--(Zee\_lightJet\_0-70)--", R"--(Zee\_lightJet\_70-140)--", R"--(Zee\_lightJet\_140-280)--", R"--(Zee\_cJet\_0-70)--", R"--(Zee\_cJet\_70-140)--", R"--(Zee\_cJet\_140-280)--", R"--(Zee\_bJet\_0-70)--", R"--(Zee\_bJet\_70-140)--", R"--(Zee\_bJet\_140-280)--"};
//
//    std::array<double,9> JetNumeratorSFs = {((139e15)*(1.9828e-9)*(0.821204)),((139e15)*(110.64e-12)*(0.69275)),((139e15)*(40.645e-12)*(0.615906)),((139e15)*(1.9817e-9)*(0.1136684)),((139e15)*(110.47e-12)*(0.1912956)),((139e15)*(40.674e-12)*(0.2326772)),((139e15)*(1.9819e-9)*(0.0656969)),((139e15)*(110.53e-12)*(0.1158741)),((139e15)*(40.68e-12)*(0.1535215))};
//    std::array<double,3> SFs = {((139e15)*(.871e-12)),((139e15)*(.199e-12)), ((139e15)*(.0345e-15))}; //numerators for Z-gamma bkg
//
//    std::vector<ROOT::RDF::RResultHandle> Totals;
//
//    for (auto& i: input_filenames)
//    {
//        SchottDataFrame df(MakeRDF(i,8));
////        std::cout << *df.Count() << '\n';
//
//        auto two_leptons = df
//        .Filter([](const RVec<std::string>& trigger_passed_triggers)
//        {
//            bool trigger_found = (std::find_first_of(trigger_passed_triggers.begin(), trigger_passed_triggers.end(), triggers.begin(), triggers.end()) != trigger_passed_triggers.end());
//
//            if (!trigger_found)
//            {
//                return false; //this event is filtered out
//            }
//            return true; //this event is kept because the trigger was found in its `trigger_passed_triggers` branch entry
//
//        }, {"trigger_passed_triggers"})
//        .Define("di_electrons", //the events that pass will have exactly 2 electrons that pass the following
//        [](RVec<Electron> electrons)
//        {
//            //keep the electrons in each event that have pt > 20 GeV, |η| < 2.37,
//            //|η| not between 1.37 and 1.52, and that satisfy a medium id criteria `electron_id_medium`
//            electrons.erase(std::remove_if(electrons.begin(),electrons.end(),
//            [](Electron& ep)
//            {
//                 return (!((ep.electron_pt/1e3 > 20) && (std::abs(ep.electron_eta) < 2.37) &&
//                 (!((1.37 < std::abs(ep.electron_eta)) && (std::abs(ep.electron_eta) < 1.52)))
//                 && (ep.electron_id_medium == 1)));
//
//            }), electrons.end());
//
//            return electrons;
//
//        },{"electrons"}).Filter(
//        [](RVec<Muon>& muons, RVec<Electron> electrons)
//        {
//            return (electrons.size()==2 && muons.empty()); //keep events which have exactly 2 electrons for di_electrons and no muons
//
//        }, {"muons", "di_electrons"});
//
//        //new dataframe node: contains only the events from `two_leptons` whose electrons in the `di_electrons` branch have opposite charge
//        auto opp_charge = two_leptons
//        .Filter([](RVec<Electron>& electrons)
//        {
//            return (electrons[0].electron_charge*electrons[1].electron_charge < 0);
//
//        }, {"di_electrons"});
//
//        //new dataframe node: contains only the events from `opp_charge` that have 1 electron with pt > 20 GeV and the other with pt > 27 GeV
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
//        //new dataframe node: Contains a new column `dilep` in addition to the ones in `same_flavour` that stores the di-electron four-vector
//        auto dilep = deltaR.Define("dilep",[] (RVec<Electron>& electrons)
//        {
//            return (electrons[0].Vector() + electrons[1].Vector());
//        }, {"di_electrons"});
//
//        //new dataframe node: contains only the events from `dilep` that have di-electron invariant mass between 81 and 101 GeV
//        auto mass = dilep.Filter([] (PtEtaPhiEVector& dilep)
//        {
//            auto mass = dilep.M()/1e3;
//            return ((mass >= 81) && (mass <= 101));
//        }, {"dilep"});
//
//        //new dataframe node: contains only the events from `mass` that have dilepton pT > 10 GeV
//        auto ptCut = mass.Filter([] (PtEtaPhiEVector& dilep)
//        {
//            auto pT = dilep.Pt()/1e3;
//            return pT > 10;
//        }, {"dilep"});
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
//                if (i==0 || ((std::abs(1-X) < std::abs(1-best_X)) and (delta_r < 1.5)))
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
//        [&](RVec<Photon>& reco_photons_matched)
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
//        }, {"photons_pass_cuts"});
//
//        auto dilepton_and_photon = merged_reco_photons_matched
//        .Define("reconstructed_mass",[&](RVec<Electron>& di_electrons, Photon& merged_photon)
//        {
//            auto four_momentum = di_electrons[0].Vector() + di_electrons[1].Vector();
//
//            return (four_momentum + merged_photon.Vector()).M()/1e3;
//
//        }, {"di_electrons", "merged_photon"});
//
//        auto pSB = dilepton_and_photon.Filter(
//        [](double reconstructed_mass)
//        {
//            return (reconstructed_mass < 110) || (reconstructed_mass > 130);
//        }, {"reconstructed_mass"});
//
//        auto pSR = dilepton_and_photon.Filter(
//        [](double reconstructed_mass)
//        {
//            return (reconstructed_mass >= 110) && (reconstructed_mass <= 130);
//        }, {"reconstructed_mass"});
//
//        auto SB = pSB.Filter(
//        [](RVec<float>& Eratio)
//        {
//            return (!Any(Eratio <= 0.8));
//        }, {"photon_shower_shape_e_ratio"});
//
//        auto SR = pSR.Filter(
//        [](RVec<float>& Eratio)
//        {
//            return (!Any(Eratio <= 0.8));
//        }, {"photon_shower_shape_e_ratio"});
//
//        Totals.push_back(df.Count());
//        Totals.push_back(merged_reco_photons_matched.Count());
//        Totals.push_back(pSB.Count());
//        Totals.push_back(pSR.Count());
//        Totals.push_back(SB.Count());
//        Totals.push_back(SR.Count());
//
//    }
//
////    0       1       2       3       4       5     //Z-gamma
////    6       7       8       9       10      11    //Z-gamma
////    12      13      14      15      16      17    //Z-gamma
////    18      19      20      21      22      23    //data
////    24      25      26      27      28      29    //ma1
////    30      31      32      33      34      35    //ma2
////    36      37      38      39      40      41    //ma3
////    42      43      44      45      46      47    //ma5
////    48      49      50      51      52      53    //ma9
////    54      55      56      57      58      59    //Z+jets
////    60      61      62      63      64      65    //Z+jets
////    66      67      68      69      70      71    //Z+jets
////    72      73      74      75      76      77    //Z+jets
////    78      79      80      81      82      83    //Z+jets
////    84      85      86      87      88      89    //Z+jets
////    90      91      92      93      94      95    //Z+jets
////    96      97      98      99      100     101   //Z+jets
////    102     103     104     105     106     107   //Z+jets
//
//    ROOT::RDF::RunGraphs(Totals); // running all computation nodes concurrently
//    std::ofstream out("Table16.txt");
//
//    out << R"--(\section*{Table 16})--" << '\n';
//    out << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
//    out << R"--(\setlength\extrarowheight{2pt})--" << '\n';
//    out << R"--(\renewcommand{\arraystretch}{1.5})--" << '\n';
//    out << R"--(\begin{tabular}{|c|c|c|c|c|c|})--" << '\n';
//    out << R"--(\hline)--" << '\n';
//    out << R"--({} & \textbf{Full Reg} & \textbf{pSB} & \textbf{pSR} & \textbf{SB} & \textbf{SR} \\ \hline)--" << '\n';
//
//    double total_Z_gamma[5] = {0}, total_Z_jets[5] = {0}, total[5] = {0};
//    double total_Z_gammaStatUnc[5] = {0}, total_Z_jetsStatUnc[5] = {0}, totalStatUnc[5] = {0};
//
//    for (int i = 0, j = 0; (i <= 102 && j <= 17); i += 6, j++)
//    {
//        out << prefixes[j];
//        if (j >= 0 && j <= 2)
//        {
//            out
//            << " & " << std::setprecision(2) << std::fixed << *Totals[i+1].GetResultPtr<ULong64_t>() * (SFs[j] / *Totals[i].GetResultPtr<ULong64_t>())
//            << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(*Totals[i+1].GetResultPtr<ULong64_t>()) * (SFs[j] / *Totals[i].GetResultPtr<ULong64_t>())
//            << " & " << std::setprecision(2) << std::fixed << *Totals[i+2].GetResultPtr<ULong64_t>() * (SFs[j] / *Totals[i].GetResultPtr<ULong64_t>())
//            << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(*Totals[i+2].GetResultPtr<ULong64_t>()) * (SFs[j] / *Totals[i].GetResultPtr<ULong64_t>())
//            << " & " << std::setprecision(2) << std::fixed << *Totals[i+3].GetResultPtr<ULong64_t>() * (SFs[j] / *Totals[i].GetResultPtr<ULong64_t>())
//            << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(*Totals[i+3].GetResultPtr<ULong64_t>()) * (SFs[j] / *Totals[i].GetResultPtr<ULong64_t>())
//            << " & " << std::setprecision(2) << std::fixed << *Totals[i+4].GetResultPtr<ULong64_t>() * (SFs[j] / *Totals[i].GetResultPtr<ULong64_t>())
//            << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(*Totals[i+4].GetResultPtr<ULong64_t>()) * (SFs[j] / *Totals[i].GetResultPtr<ULong64_t>())
//            << " & " << std::setprecision(2) << std::fixed << *Totals[i+5].GetResultPtr<ULong64_t>() * (SFs[j] / *Totals[i].GetResultPtr<ULong64_t>())
//            << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(*Totals[i+5].GetResultPtr<ULong64_t>()) * (SFs[j] / *Totals[i].GetResultPtr<ULong64_t>())
//            << R"--( \\ \hline)--" << '\n';
//        }
//
//        else if (j >= 9)
//        {
//            out
//            << " & " << std::setprecision(2) << std::fixed << *Totals[i+1].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j-9] / *Totals[i].GetResultPtr<ULong64_t>())
//            << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(*Totals[i+1].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[j-9] / *Totals[i].GetResultPtr<ULong64_t>())
//            << " & " << std::setprecision(2) << std::fixed << *Totals[i+2].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j-9] / *Totals[i].GetResultPtr<ULong64_t>())
//            << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(*Totals[i+2].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[j-9] / *Totals[i].GetResultPtr<ULong64_t>())
//            << " & " << std::setprecision(2) << std::fixed << *Totals[i+3].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j-9] / *Totals[i].GetResultPtr<ULong64_t>())
//            << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(*Totals[i+3].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[j-9] / *Totals[i].GetResultPtr<ULong64_t>())
//            << " & " << std::setprecision(2) << std::fixed << *Totals[i+4].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j-9] / *Totals[i].GetResultPtr<ULong64_t>())
//            << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(*Totals[i+4].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[j-9] / *Totals[i].GetResultPtr<ULong64_t>())
//            << " & " << std::setprecision(2) << std::fixed << *Totals[i+5].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j-9] / *Totals[i].GetResultPtr<ULong64_t>())
//            << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(*Totals[i+5].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[j-9] / *Totals[i].GetResultPtr<ULong64_t>())
//            << R"--( \\ \hline)--" << '\n';
//        }
//
//        else
//        {
//            out
//            << " & " << std::setprecision(2) << std::fixed << *Totals[i+1].GetResultPtr<ULong64_t>()
//            << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(*Totals[i+1].GetResultPtr<ULong64_t>())
//            << " & " << std::setprecision(2) << std::fixed << *Totals[i+2].GetResultPtr<ULong64_t>()
//            << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(*Totals[i+2].GetResultPtr<ULong64_t>())
//            << " & " << std::setprecision(2) << std::fixed << *Totals[i+3].GetResultPtr<ULong64_t>()
//            << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(*Totals[i+3].GetResultPtr<ULong64_t>())
//            << " & " << std::setprecision(2) << std::fixed << *Totals[i+4].GetResultPtr<ULong64_t>()
//            << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(*Totals[i+4].GetResultPtr<ULong64_t>())
//            << " & " << std::setprecision(2) << std::fixed << *Totals[i+5].GetResultPtr<ULong64_t>()
//            << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(*Totals[i+5].GetResultPtr<ULong64_t>())
//            << R"--( \\ \hline)--" << '\n';
//        }
//
//        if (i >= 0 && i <= 12) //Z-gamma
//        {
//            for (int k = 0; k < 5; k++)
//            {
//                total_Z_gamma[k] += *Totals[i+k+1].GetResultPtr<ULong64_t>() * (SFs[i/6] / *Totals[i].GetResultPtr<ULong64_t>());
//                total[k] += *Totals[i+k+1].GetResultPtr<ULong64_t>() * (SFs[i/6] / *Totals[i].GetResultPtr<ULong64_t>());
//
//                total_Z_gammaStatUnc[k] += pow(sqrt(*Totals[i+k+1].GetResultPtr<ULong64_t>()) * (SFs[i/6] / *Totals[i].GetResultPtr<ULong64_t>()),2);
//                totalStatUnc[k] += pow(sqrt(*Totals[i+k+1].GetResultPtr<ULong64_t>()) * (SFs[i/6] / *Totals[i].GetResultPtr<ULong64_t>()),2);
//            }
//        }
//
//        else if (i >= 54 && i <= 102) //Z-jets
//        {
//            for (int k = 0; k < 5; k++)
//            {
//                total_Z_jets[k] += *Totals[i+k+1].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i/6 - 9] / *Totals[i].GetResultPtr<ULong64_t>());
//                total[k] += *Totals[i+k+1].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i/6 - 9] / *Totals[i].GetResultPtr<ULong64_t>());
//
//                total_Z_jetsStatUnc[k] += pow(sqrt(*Totals[i+k+1].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i/6 - 9] / *Totals[i].GetResultPtr<ULong64_t>()),2);
//                totalStatUnc[k] += pow(sqrt(*Totals[i+k+1].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i/6 - 9] / *Totals[i].GetResultPtr<ULong64_t>()),2);
//            }
//        }
//    }
//
//    std::vector<std::string> totalPrefixes = {R"--(Total $Z\gamma\gamma$)--", R"--(Total $Z$+jets)--", R"--(Total Bkg)--"};
//
//    out << R"--(Total $Z\gamma\gamma$)--";
//    out << " & " << std::setprecision(2) << std::fixed << total_Z_gamma[0]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(total_Z_gammaStatUnc[0]);
//    out << " & " << std::setprecision(2) << std::fixed << total_Z_gamma[1]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(total_Z_gammaStatUnc[1]);
//    out << " & " << std::setprecision(2) << std::fixed << total_Z_gamma[2]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(total_Z_gammaStatUnc[2]);
//    out << " & " << std::setprecision(2) << std::fixed << total_Z_gamma[3]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(total_Z_gammaStatUnc[3]);
//    out << " & " << std::setprecision(2) << std::fixed << total_Z_gamma[4]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(total_Z_gammaStatUnc[4]);
//    out << R"--( \\ \hline)--" << '\n';
//
//    out << R"--(Total $Z$+jets)--";
//    out << " & " << std::setprecision(2) << std::fixed << total_Z_jets[0]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(total_Z_jetsStatUnc[0]);
//    out << " & " << std::setprecision(2) << std::fixed << total_Z_jets[1]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(total_Z_jetsStatUnc[1]);
//    out << " & " << std::setprecision(2) << std::fixed << total_Z_jets[2]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(total_Z_jetsStatUnc[2]);
//    out << " & " << std::setprecision(2) << std::fixed << total_Z_jets[3]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(total_Z_jetsStatUnc[3]);
//    out << " & " << std::setprecision(2) << std::fixed << total_Z_jets[4]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(total_Z_jetsStatUnc[4]);
//    out << R"--( \\ \hline)--" << '\n';
//
//    out << R"--(Total Bkg)--";
//    out << " & " << std::setprecision(2) << std::fixed << total[0]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(totalStatUnc[0]);
//    out << " & " << std::setprecision(2) << std::fixed << total[1]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(totalStatUnc[1]);
//    out << " & " << std::setprecision(2) << std::fixed << total[2]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(totalStatUnc[2]);
//    out << " & " << std::setprecision(2) << std::fixed << total[3]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(totalStatUnc[3]);
//    out << " & " << std::setprecision(2) << std::fixed << total[4]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(totalStatUnc[4]);
//    out << R"--( \\ \hline)--" << '\n';
//
//    out << R"--(\end{tabular}})--" << '\n';
//
//    out << "\n\n\n";
//
//    out.close();
//}
//
//void Table19()
//{
//    std::vector<std::vector<std::string>> input_filenames =
//    {
//        //Z-gamma
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617070._000001.LGNTuple.root"}, {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617064._000001.LGNTuple.root"},
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617074._000001.LGNTuple.root"},
//        //data
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_data_test.root"},
//        //signal
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root"}, //1 GeV
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600751.PhPy8EG_AZNLO_ggH125_mA2p0_v1.root"}, // 2 GeV
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600752.PhPy8EG_AZNLO_ggH125_mA3p0_v1.root"}, // 3 GeV
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root"}, //5 GeV
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600756.PhPy8EG_AZNLO_ggH125_mA9p0_v1.root"}, // 9 GeV
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
//    std::vector<std::string> prefixes = { R"--(pty2\_9\_17)--", R"--(pty\_17\_myy\_0\_80)--", R"--(pty\_17\_myy\_80)--", "data", R"--($\text{Sig } m_{A}$ = 1 GeV)--", R"--($\text{Sig } m_{A}$ = 2 GeV)--", R"--($\text{Sig } m_{A}$ = 3 GeV)--", R"--($\text{Sig } m_{A}$ = 5 GeV)--", R"--($\text{Sig } m_{A}$ = 9 GeV)--", R"--(Zee\_lightJet\_0-70)--", R"--(Zee\_lightJet\_70-140)--", R"--(Zee\_lightJet\_140-280)--", R"--(Zee\_cJet\_0-70)--", R"--(Zee\_cJet\_70-140)--", R"--(Zee\_cJet\_140-280)--", R"--(Zee\_bJet\_0-70)--", R"--(Zee\_bJet\_70-140)--", R"--(Zee\_bJet\_140-280)--"};
//
//    std::array<double,9> JetNumeratorSFs = {((139e15)*(1.9828e-9)*(0.821204)),((139e15)*(110.64e-12)*(0.69275)),((139e15)*(40.645e-12)*(0.615906)),((139e15)*(1.9817e-9)*(0.1136684)),((139e15)*(110.47e-12)*(0.1912956)),((139e15)*(40.674e-12)*(0.2326772)),((139e15)*(1.9819e-9)*(0.0656969)),((139e15)*(110.53e-12)*(0.1158741)),((139e15)*(40.68e-12)*(0.1535215))};
//    std::array<double,3> SFs = {((139e15)*(.871e-12)),((139e15)*(.199e-12)), ((139e15)*(.0345e-15))}; //numerators for Z-gamma bkg
//
//    std::vector<ROOT::RDF::RResultHandle> Totals;
//
//    for (auto& i: input_filenames)
//    {
//        SchottDataFrame df(MakeRDF(i,8));
////        std::cout << *df.Count() << '\n';
//
//        auto two_leptons = df
//        .Filter([](const RVec<std::string>& trigger_passed_triggers)
//        {
//            bool trigger_found = (std::find_first_of(trigger_passed_triggers.begin(), trigger_passed_triggers.end(), triggers.begin(), triggers.end()) != trigger_passed_triggers.end());
//
//            if (!trigger_found)
//            {
//                return false; //this event is filtered out
//            }
//            return true; //this event is kept because the trigger was found in its `trigger_passed_triggers` branch entry
//
//        }, {"trigger_passed_triggers"})
//        .Define("di_electrons", //the events that pass will have exactly 2 electrons that pass the following
//        [](RVec<Electron> electrons)
//        {
//            //keep the electrons in each event that have pt > 20 GeV, |η| < 2.37,
//            //|η| not between 1.37 and 1.52, and that satisfy a medium id criteria `electron_id_medium`
//            electrons.erase(std::remove_if(electrons.begin(),electrons.end(),
//            [](Electron& ep)
//            {
//                 return (!((ep.electron_pt/1e3 > 20) && (std::abs(ep.electron_eta) < 2.37) &&
//                 (!((1.37 < std::abs(ep.electron_eta)) && (std::abs(ep.electron_eta) < 1.52)))
//                 && (ep.electron_id_medium == 1)));
//
//            }), electrons.end());
//
//            return electrons;
//
//        },{"electrons"}).Filter(
//        [](RVec<Muon>& muons, RVec<Electron> electrons)
//        {
//            return (electrons.size()==2 && muons.empty()); //keep events which have exactly 2 electrons for di_electrons and no muons
//
//        }, {"muons", "di_electrons"});
//
//        //new dataframe node: contains only the events from `two_leptons` whose electrons in the `di_electrons` branch have opposite charge
//        auto opp_charge = two_leptons
//        .Filter([](RVec<Electron>& electrons)
//        {
//            return (electrons[0].electron_charge*electrons[1].electron_charge < 0);
//
//        }, {"di_electrons"});
//
//        //new dataframe node: contains only the events from `opp_charge` that have 1 electron with pt > 20 GeV and the other with pt > 27 GeV
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
//        //new dataframe node: Contains a new column `dilep` in addition to the ones in `same_flavour` that stores the di-electron four-vector
//        auto dilep = deltaR.Define("dilep",[] (RVec<Electron>& electrons)
//        {
//            return (electrons[0].Vector() + electrons[1].Vector());
//        }, {"di_electrons"});
//
//        //new dataframe node: contains only the events from `dilep` that have di-electron invariant mass between 81 and 101 GeV
//        auto mass = dilep.Filter([] (PtEtaPhiEVector& dilep)
//        {
//            auto mass = dilep.M()/1e3;
//            return ((mass >= 81) && (mass <= 101));
//        }, {"dilep"});
//
//        //new dataframe node: contains only the events from `mass` that have dilepton pT > 10 GeV
//        auto ptCut = mass.Filter([] (PtEtaPhiEVector& dilep)
//        {
//            auto pT = dilep.Pt()/1e3;
//            return pT > 10;
//        }, {"dilep"});
//
//        auto photon_passes_cuts = ptCut.Define("photons_pass_cuts",
//        [&](RVec<Photon> photons)
//        {
//            photons.erase(std::remove_if(photons.begin(),photons.end(),
//            [](Photon& x)
//            {
//                return ((std::abs(x.photon_eta) >= 2.37) || (std::abs(x.photon_eta) > 1.37 && std::abs(x.photon_eta) < 1.52));
//
//            }), photons.end());
//
//            return photons;
//        }, {"photons"});
//
//        auto failed_resolved = photon_passes_cuts.Filter(
//        [&](RVec<Photon>& reco_photons_test)
//        {
//            RVec<Photon> reco_photons_matched = reco_photons_test;
//            if (reco_photons_matched.size() == 1 || reco_photons_matched.empty())
//            {
//                return true;
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
//                if (i==0 || ((std::abs(1-X) < std::abs(1-best_X)) and (delta_r < 1.5)))
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
//            return true;
//
//        }, {"photons_pass_cuts"});
//
//        auto photon_pt_cut = failed_resolved.Filter(
//        [&](RVec<Photon>& photon_passes_cuts)
//        {
//            for (auto& p: photon_passes_cuts)
//            {
//                if (p.photon_pt > 20e3)
//                {
//                    return true;
//                }
//            }
//
//            return false;
//
//        }, {"photons_pass_cuts"}).Define("merged_photon",
//        [&](RVec<Photon>& photon_passes_cuts)
//        {
//            return photon_passes_cuts[0];
//        }, {"photons_pass_cuts"});
//
//        auto dilepton_and_photon = photon_pt_cut
//        .Define("reconstructed_mass",[&](RVec<Electron>& di_electrons, Photon& merged_photon)
//        {
//            auto four_momentum = di_electrons[0].Vector() + di_electrons[1].Vector();
//
//            return (four_momentum + merged_photon.Vector()).M()/1e3;
//
//        }, {"di_electrons", "merged_photon"});
//
//        auto pSR = dilepton_and_photon.Filter(
//        [](double reconstructed_mass)
//        {
//            return (reconstructed_mass >= 110) && (reconstructed_mass <= 130);
//        }, {"reconstructed_mass"});
//
//        auto SR = pSR.Filter(
//        [](RVec<float>& Eratio)
//        {
//            return (!Any(Eratio <= 0.8));
//        }, {"photon_shower_shape_e_ratio"});
//
//        auto SR_ID = SR.Filter(
//        [&](Photon& merged_photon)
//        {
//            return merged_photon.photon_id_loose;
//        }, {"merged_photon"});
//
//        Totals.push_back(df.Count());
//        Totals.push_back(ptCut.Count()); //preselection
//        Totals.push_back(failed_resolved.Count()); //failed_resolved
//        Totals.push_back(photon_pt_cut.Count()); //photon_pt_cut
//        Totals.push_back(pSR.Count()); //pSR
//        Totals.push_back(SR.Count()); //SR
//        Totals.push_back(SR_ID.Count()); //SR_ID
//    }
//
//    ROOT::RDF::RunGraphs(Totals);
//
////    0       1       2       3       4       5       6       //Z-gamma
////    7       8       9       10      11      12      13      //Z-gamma
////    14      15      16      17      18      19      20      //Z-gamma
////    21      22      23      24      25      26      27      //data
////    28      29      30      31      32      33      34      //ma1
////    35      36      37      38      39      40      41      //ma2
////    42      43      44      45      46      47      48      //ma3
////    49      50      51      52      53      54      55      //ma5
////    56      57      58      59      60      61      62      //ma9
////    63      64      65      66      67      68      69      //Z+jets
////    70      71      72      73      74      75      76      //Z+jets
////    77      78      79      80      81      82      83      //Z+jets
////    84      85      86      87      88      89      90      //Z+jets
////    91      92      93      94      95      96      97      //Z+jets
////    98      99      100     101     102     103     104     //Z+jets
////    105     106     107     108     109     110     111     //Z+jets
////    112     113     114     115     116     117     118     //Z+jets
////    119     120     121     122     123     124     125     //Z+jets
//
//    std::ofstream out("Table19.txt");
//    out << R"--(\section*{Table 19})--" << '\n';
//    out << R"--(\hspace{-3cm}\scalebox{0.55}{)--" << '\n';
//    out << R"--(\setlength\extrarowheight{2pt})--" << '\n';
//    out << R"--(\renewcommand{\arraystretch}{1.5})--" << '\n';
//    out << R"--(\begin{tabular}{|c|c|c|c|c|c|c|})--" << '\n';
//    out << R"--(\hline)--" << '\n';
//    out << R"--({} & \textbf{pass preselection} & \textbf{failed resolved category} & \textbf{photon} $\pmb{p_T}$ \textbf{cut} & \textbf{pSR} & \textbf{SR} & \textbf{SR-ID} \\ \hline)--" << '\n';
//
//    double total_Z_gamma[6] = {0}, total_Z_jets[6] = {0}, total[6] = {0};
//    double total_Z_gammaStatUnc[6] = {0}, total_Z_jetsStatUnc[6] = {0}, totalStatUnc[6] = {0};
//
//    for (int i = 0, j = 0; (i <= 119 && j <= 17); i += 7, j++)
//    {
//        out << prefixes[j];
//        if (j >= 0 && j <= 2)
//        {
//            out
//            << " & " << std::setprecision(2) << std::fixed << *Totals[i+1].GetResultPtr<ULong64_t>() * (SFs[j] / *Totals[i].GetResultPtr<ULong64_t>())
//            << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(*Totals[i+1].GetResultPtr<ULong64_t>()) * (SFs[j] / *Totals[i].GetResultPtr<ULong64_t>())
//            << " & " << std::setprecision(2) << std::fixed << *Totals[i+2].GetResultPtr<ULong64_t>() * (SFs[j] / *Totals[i].GetResultPtr<ULong64_t>())
//            << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(*Totals[i+2].GetResultPtr<ULong64_t>()) * (SFs[j] / *Totals[i].GetResultPtr<ULong64_t>())
//            << " & " << std::setprecision(2) << std::fixed << *Totals[i+3].GetResultPtr<ULong64_t>() * (SFs[j] / *Totals[i].GetResultPtr<ULong64_t>())
//            << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(*Totals[i+3].GetResultPtr<ULong64_t>()) * (SFs[j] / *Totals[i].GetResultPtr<ULong64_t>())
//            << " & " << std::setprecision(2) << std::fixed << *Totals[i+4].GetResultPtr<ULong64_t>() * (SFs[j] / *Totals[i].GetResultPtr<ULong64_t>())
//            << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(*Totals[i+4].GetResultPtr<ULong64_t>()) * (SFs[j] / *Totals[i].GetResultPtr<ULong64_t>())
//            << " & " << std::setprecision(2) << std::fixed << *Totals[i+5].GetResultPtr<ULong64_t>() * (SFs[j] / *Totals[i].GetResultPtr<ULong64_t>())
//            << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(*Totals[i+5].GetResultPtr<ULong64_t>()) * (SFs[j] / *Totals[i].GetResultPtr<ULong64_t>())
//            << " & " << std::setprecision(2) << std::fixed << *Totals[i+6].GetResultPtr<ULong64_t>() * (SFs[j] / *Totals[i].GetResultPtr<ULong64_t>())
//            << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(*Totals[i+6].GetResultPtr<ULong64_t>()) * (SFs[j] / *Totals[i].GetResultPtr<ULong64_t>())
//            << R"--( \\ \hline)--" << '\n';
//        }
//
//        else if (j >= 9)
//        {
//            out
//            << " & " << std::setprecision(2) << std::fixed << *Totals[i+1].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j-9] / *Totals[i].GetResultPtr<ULong64_t>())
//            << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(*Totals[i+1].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[j-9] / *Totals[i].GetResultPtr<ULong64_t>())
//            << " & " << std::setprecision(2) << std::fixed << *Totals[i+2].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j-9] / *Totals[i].GetResultPtr<ULong64_t>())
//            << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(*Totals[i+2].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[j-9] / *Totals[i].GetResultPtr<ULong64_t>())
//            << " & " << std::setprecision(2) << std::fixed << *Totals[i+3].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j-9] / *Totals[i].GetResultPtr<ULong64_t>())
//            << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(*Totals[i+1].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[j-9] / *Totals[i].GetResultPtr<ULong64_t>())
//            << " & " << std::setprecision(2) << std::fixed << *Totals[i+4].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j-9] / *Totals[i].GetResultPtr<ULong64_t>())
//            << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(*Totals[i+1].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[j-9] / *Totals[i].GetResultPtr<ULong64_t>())
//            << " & " << std::setprecision(2) << std::fixed << *Totals[i+5].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j-9] / *Totals[i].GetResultPtr<ULong64_t>())
//            << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(*Totals[i+1].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[j-9] / *Totals[i].GetResultPtr<ULong64_t>())
//            << " & " << std::setprecision(2) << std::fixed << *Totals[i+6].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j-9] / *Totals[i].GetResultPtr<ULong64_t>())
//            << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(*Totals[i+1].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[j-9] / *Totals[i].GetResultPtr<ULong64_t>())
//            << R"--( \\ \hline)--" << '\n';
//        }
//
//        else
//        {
//            out
//            << " & " << std::setprecision(2) << std::fixed << *Totals[i+1].GetResultPtr<ULong64_t>()
//            << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(*Totals[i+1].GetResultPtr<ULong64_t>())
//            << " & " << std::setprecision(2) << std::fixed << *Totals[i+2].GetResultPtr<ULong64_t>()
//            << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(*Totals[i+2].GetResultPtr<ULong64_t>())
//            << " & " << std::setprecision(2) << std::fixed << *Totals[i+3].GetResultPtr<ULong64_t>()
//            << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(*Totals[i+3].GetResultPtr<ULong64_t>())
//            << " & " << std::setprecision(2) << std::fixed << *Totals[i+4].GetResultPtr<ULong64_t>()
//            << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(*Totals[i+4].GetResultPtr<ULong64_t>())
//            << " & " << std::setprecision(2) << std::fixed << *Totals[i+5].GetResultPtr<ULong64_t>()
//            << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(*Totals[i+5].GetResultPtr<ULong64_t>())
//            << " & " << std::setprecision(2) << std::fixed << *Totals[i+6].GetResultPtr<ULong64_t>()
//            << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(*Totals[i+6].GetResultPtr<ULong64_t>())
//            << R"--( \\ \hline)--" << '\n';
//        }
//
//        if (i >= 0 && i <= 14) //Z-gamma
//        {
//            for (int k = 0; k < 6; k++)
//            {
//                total_Z_gamma[k] += *Totals[i+k+1].GetResultPtr<ULong64_t>() * (SFs[i/7] / *Totals[i].GetResultPtr<ULong64_t>());
//                total[k] += *Totals[i+k+1].GetResultPtr<ULong64_t>() * (SFs[i/7] / *Totals[i].GetResultPtr<ULong64_t>());
//
//                total_Z_gammaStatUnc[k] += pow(sqrt(*Totals[i+k+1].GetResultPtr<ULong64_t>()) * (SFs[i/7] / *Totals[i].GetResultPtr<ULong64_t>()),2);
//                totalStatUnc[k] += pow(sqrt(*Totals[i+k+1].GetResultPtr<ULong64_t>()) * (SFs[i/7] / *Totals[i].GetResultPtr<ULong64_t>()),2);
//            }
//        }
//
//        else if (i >= 63 && i <= 119) //Z-jets
//        {
//            for (int k = 0; k < 6; k++)
//            {
//                total_Z_jets[k] += *Totals[i+k+1].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i/7 - 9] / *Totals[i].GetResultPtr<ULong64_t>());
//                total[k] += *Totals[i+k+1].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i/7 - 9] / *Totals[i].GetResultPtr<ULong64_t>());
//
//                total_Z_jetsStatUnc[k] += pow(sqrt(*Totals[i+k+1].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i/7 - 9] / *Totals[i].GetResultPtr<ULong64_t>()),2);
//                totalStatUnc[k] += pow(sqrt(*Totals[i+k+1].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i/7 - 9] / *Totals[i].GetResultPtr<ULong64_t>()),2);
//            }
//        }
//    }
//
//    std::vector<std::string> totalPrefixes = {R"--(Total $Z\gamma\gamma$)--", R"--(Total $Z$+jets)--", R"--(Total Bkg)--"};
//
//    out << R"--(Total $Z\gamma\gamma$)--";
//    out << " & " << std::setprecision(2) << std::fixed << total_Z_gamma[0]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(total_Z_gammaStatUnc[0]);
//    out << " & " << std::setprecision(2) << std::fixed << total_Z_gamma[1]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(total_Z_gammaStatUnc[1]);
//    out << " & " << std::setprecision(2) << std::fixed << total_Z_gamma[2]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(total_Z_gammaStatUnc[2]);
//    out << " & " << std::setprecision(2) << std::fixed << total_Z_gamma[3]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(total_Z_gammaStatUnc[3]);
//    out << " & " << std::setprecision(2) << std::fixed << total_Z_gamma[4]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(total_Z_gammaStatUnc[4]);
//    out << " & " << std::setprecision(2) << std::fixed << total_Z_gamma[5]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(total_Z_gammaStatUnc[5]);
//    out << R"--( \\ \hline)--" << '\n';
//
//    out << R"--(Total $Z$+jets)--";
//    out << " & " << std::setprecision(2) << std::fixed << total_Z_jets[0]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(total_Z_jetsStatUnc[0]);
//    out << " & " << std::setprecision(2) << std::fixed << total_Z_jets[1]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(total_Z_jetsStatUnc[1]);
//    out << " & " << std::setprecision(2) << std::fixed << total_Z_jets[2]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(total_Z_jetsStatUnc[2]);
//    out << " & " << std::setprecision(2) << std::fixed << total_Z_jets[3]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(total_Z_jetsStatUnc[3]);
//    out << " & " << std::setprecision(2) << std::fixed << total_Z_jets[4]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(total_Z_jetsStatUnc[4]);
//    out << " & " << std::setprecision(2) << std::fixed << total_Z_jets[5]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(total_Z_jetsStatUnc[5]);
//
//    out << R"--( \\ \hline)--" << '\n';
//
//    out << R"--(Total Bkg)--";
//    out << " & " << std::setprecision(2) << std::fixed << total[0]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(totalStatUnc[0]);
//    out << " & " << std::setprecision(2) << std::fixed << total[1]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(totalStatUnc[1]);
//    out << " & " << std::setprecision(2) << std::fixed << total[2]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(totalStatUnc[2]);
//    out << " & " << std::setprecision(2) << std::fixed << total[3]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(totalStatUnc[3]);
//    out << " & " << std::setprecision(2) << std::fixed << total[4]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(totalStatUnc[4]);
//    out << " & " << std::setprecision(2) << std::fixed << total[5]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(totalStatUnc[5]);
//    out << R"--( \\ \hline)--" << '\n';
//
//    out << R"--(\end{tabular}})--" << '\n';
//
//    out << "\n\n\n";
//
//    out.close();
//}

//void Table16_Displaced_Axions()
//{
//    std::vector<std::vector<std::string>> input_filenames = {
//        //Signal
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root"},
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root"},
//        //Displaced Signal
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_FewMassPoints.root"},
//    };
//
//    constexpr std::array<const char*,15> prefixes = {R"--(Prompt Signal $m_{\text{A}}$ = 1 GeV)--", R"--(Prompt Signal $m_{\text{A}}$ = 5 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 0.2 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 0.5 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 1 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 3 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 5 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 10 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 10.1 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 20 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 20.1 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 29.5 GeV)--", };
//
//    std::vector<double> massPoints = {0.2,0.5,1,3,5,10,10.1,20,20.1,29.5};
//
//    std::vector<ROOT::RDF::RResultHandle> Totals;
//
//    int counter = 0;
//
//    for (auto& i: input_filenames)
//    {
//        SchottDataFrame df(MakeRDF(i,8));
////        std::cout << *df.Count() << '\n';
//
//        auto trigger_selection = df.Define("truth_axions", [&](RVec<TruthParticle> truth_particles)
//        {
//            truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
//            [](TruthParticle& x)
//            {
//               return (std::abs(x.mc_pdg_id) != 36);
//
//            }), truth_particles.end());
//
//            return truth_particles;
//
//        }, {"truth_particles"}).Define("axion_masses", [&](RVec<TruthParticle>& truth_axions)
//        {
//            return truth_axions[0].mc_mass/1e3f;
//
//        }, {"truth_axions"}).Filter(
//        [](const RVec<std::string>& trigger_passed_triggers)
//        {
//            bool trigger_found = (std::find_first_of(trigger_passed_triggers.begin(), trigger_passed_triggers.end(), triggers.begin(), triggers.end()) != trigger_passed_triggers.end());
//
//            if (!trigger_found)
//            {
//             return false;
//            }
//            return true;
//
//        }, {"trigger_passed_triggers"});
//
//        auto two_leptons = trigger_selection.Filter(
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
//        [&](RVec<Photon>& reco_photons_matched)
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
//        }, {"photons_pass_cuts"});
//
//        auto dilepton_and_photon = merged_reco_photons_matched
//        .Define("reconstructed_mass",[&](RVec<Electron>& di_electrons, Photon& merged_photon)
//        {
//            auto four_momentum = di_electrons[0].Vector() + di_electrons[1].Vector();
//
//            return (four_momentum + merged_photon.Vector()).M()/1e3;
//
//        }, {"di_electrons", "merged_photon"});
//
//        auto pSB = dilepton_and_photon.Filter(
//        [](double reconstructed_mass)
//        {
//            return (reconstructed_mass < 110) || (reconstructed_mass > 130);
//        }, {"reconstructed_mass"});
//
//        auto pSR = dilepton_and_photon.Filter(
//        [](double reconstructed_mass)
//        {
//            return (reconstructed_mass >= 110) && (reconstructed_mass <= 130);
//        }, {"reconstructed_mass"});
//
//        auto SB = pSB.Filter(
//        [](RVec<float>& Eratio)
//        {
//            return (!Any(Eratio <= 0.8));
//        }, {"photon_shower_shape_e_ratio"});
//
//        auto SR = pSR.Filter(
//        [](RVec<float>& Eratio)
//        {
//            return (!Any(Eratio <= 0.8));
//        }, {"photon_shower_shape_e_ratio"});
//
//        if (counter < 2)
//        {
//            Totals.push_back(trigger_selection.Count());
//            Totals.push_back(merged_reco_photons_matched.Count());
//            Totals.push_back(pSB.Count());
//            Totals.push_back(pSR.Count());
//            Totals.push_back(SB.Count());
//            Totals.push_back(SR.Count());
//        }
//        else
//        {
//            for (auto& mass_point: massPoints)
//            {
//                auto mass_point_trigger_selection = trigger_selection.Filter([&]
//                (float axion_mass)
//                {
//                    return (roundToOneDecimalPlace(axion_mass) == mass_point);
//
//                }, {"axion_masses"});
//
//                auto mass_point_merged_reco_photons_matched = merged_reco_photons_matched.Filter([&]
//                (float axion_mass)
//                {
//                    return (roundToOneDecimalPlace(axion_mass) == mass_point);
//
//                }, {"axion_masses"});
//
//                auto mass_point_pSB = pSB.Filter([&]
//                (float axion_mass)
//                {
//                    return (roundToOneDecimalPlace(axion_mass) == mass_point);
//
//                }, {"axion_masses"});
//
//                auto mass_point_pSR = pSR.Filter([&]
//                (float axion_mass)
//                {
//                    return (roundToOneDecimalPlace(axion_mass) == mass_point);
//
//                }, {"axion_masses"});
//
//                auto mass_point_SB = SB.Filter([&]
//                (float axion_mass)
//                {
//                    return (roundToOneDecimalPlace(axion_mass) == mass_point);
//
//                }, {"axion_masses"});
//
//                auto mass_point_SR = SR.Filter([&]
//                (float axion_mass)
//                {
//                    return (roundToOneDecimalPlace(axion_mass) == mass_point);
//
//                }, {"axion_masses"});
//
//                Totals.push_back(mass_point_trigger_selection.Count());
//                Totals.push_back(mass_point_merged_reco_photons_matched.Count());
//                Totals.push_back(mass_point_pSB.Count());
//                Totals.push_back(mass_point_pSR.Count());
//                Totals.push_back(mass_point_SB.Count());
//                Totals.push_back(mass_point_SR.Count());
//            }
//        }
//
//        counter++;
//    }
//
////    0    1    2    3    4    5      ma1
////    6    7    8    9    10   11     ma
////    12   13   14   15   16   17     displaced_axion_1
////    18   19   20   21   22   23     displaced_axion_2
////    24   25   26   27   28   29     displaced_axion_3
////    30   31   32   33   34   35     displaced_axion_4
////    36   37   38   39   40   41     displaced_axion_5
////    42   43   44   45   46   47     displaced_axion_6
////    48   49   50   51   52   53     displaced_axion_7
////    54   55   56   57   58   59     displaced_axion_8
////    60   61   62   63   64   65     displaced_axion_9
////    66   67   68   69   70   71     displaced_axion_10
//
//    ROOT::RDF::RunGraphs(Totals); // running all computation nodes concurrently
//
//    std::cout << R"--(\section*{Table 16 Prompt and Displaced Signal Samples})--" << '\n';
//    std::cout << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
//    std::cout << R"--(\begin{tabular}{|c|c|c|c|c|c|})--" << '\n';
//    std::cout << R"--(\hline)--" << '\n';
//    std::cout << R"--({} & Full Reg & pSB & pSR & SB & SR \\ \hline)--" << '\n';
//
//    for (int i = 0, j = 0; (i <= 66 && j < 12); i += 6, j++)
//    {
//        std::cout << prefixes[j];
//
//        std::cout
//        << " & " << std::setprecision(2) << std::fixed << *Totals[i+1].GetResultPtr<ULong64_t>()
//        << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(*Totals[i+1].GetResultPtr<ULong64_t>())
//        << " & " << std::setprecision(2) << std::fixed << *Totals[i+2].GetResultPtr<ULong64_t>()
//        << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(*Totals[i+2].GetResultPtr<ULong64_t>())
//        << " & " << std::setprecision(2) << std::fixed << *Totals[i+3].GetResultPtr<ULong64_t>()
//        << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(*Totals[i+3].GetResultPtr<ULong64_t>())
//        << " & " << std::setprecision(2) << std::fixed << *Totals[i+4].GetResultPtr<ULong64_t>()
//        << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(*Totals[i+4].GetResultPtr<ULong64_t>())
//        << " & " << std::setprecision(2) << std::fixed << *Totals[i+5].GetResultPtr<ULong64_t>()
//        << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(*Totals[i+5].GetResultPtr<ULong64_t>())
//        << R"--( \\ \hline)--" << '\n';
//    }
//
//    std::cout << R"--(\end{tabular}})--" << '\n';
//
//    std::cout << "\n\n\n";
//}
//
//void Table19_Displaced_Axions()
//{
//    std::vector<std::vector<std::string>> input_filenames = {
//        //Signal
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root"},
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root"},
//        //Displaced Signal
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_FewMassPoints.root"},
//    };
//
//    constexpr std::array<const char*,15> prefixes = {R"--(Prompt Signal $m_{\text{A}}$ = 1 GeV)--", R"--(Prompt Signal $m_{\text{A}}$ = 5 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 0.2 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 0.5 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 1 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 3 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 5 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 10 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 10.1 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 20 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 20.1 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 29.5 GeV)--", };
//
//    std::vector<double> massPoints = {0.2,0.5,1,3,5,10,10.1,20,20.1,29.5};
//
//    std::vector<ROOT::RDF::RResultHandle> Totals;
//
//    int counter = 0;
//
//    for (auto& i: input_filenames)
//    {
//        SchottDataFrame df(MakeRDF(i,8));
////        std::cout << *df.Count() << '\n';
//
//        auto trigger_selection = df.Define("truth_axions", [&](RVec<TruthParticle> truth_particles)
//        {
//            truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
//            [](TruthParticle& x)
//            {
//               return (std::abs(x.mc_pdg_id) != 36);
//
//            }), truth_particles.end());
//
//            return truth_particles;
//
//        }, {"truth_particles"}).Define("axion_masses", [&](RVec<TruthParticle>& truth_axions)
//        {
//            return truth_axions[0].mc_mass/1e3f;
//
//        }, {"truth_axions"}).Filter(
//        [](const RVec<std::string>& trigger_passed_triggers)
//        {
//            bool trigger_found = (std::find_first_of(trigger_passed_triggers.begin(), trigger_passed_triggers.end(), triggers.begin(), triggers.end()) != trigger_passed_triggers.end());
//
//            if (!trigger_found)
//            {
//             return false;
//            }
//            return true;
//
//        }, {"trigger_passed_triggers"});
//
//        auto two_leptons = trigger_selection.Filter(
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
//                return ((std::abs(x.photon_eta) >= 2.37) || (std::abs(x.photon_eta) > 1.37 && std::abs(x.photon_eta) < 1.52));
//
//            }), photons.end());
//
//            return photons;
//        }, {"photons"});
//
//        auto failed_resolved = photon_passes_cuts.Filter(
//        [&](RVec<Photon>& reco_photons_test)
//        {
//            RVec<Photon> reco_photons_matched = reco_photons_test;
//            if (reco_photons_matched.size() == 1 || reco_photons_matched.empty())
//            {
//                return true;
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
//            return true;
//
//        }, {"photons_pass_cuts"});
//
//        auto photon_pt_cut = failed_resolved.Filter(
//        [&](RVec<Photon>& photon_passes_cuts)
//        {
//            for (auto& p: photon_passes_cuts)
//            {
//                if (p.photon_pt > 20e3)
//                {
//                    return true;
//                }
//            }
//
//            return false;
//
//        }, {"photons_pass_cuts"}).Define("merged_photon",
//        [&](RVec<Photon>& photon_passes_cuts)
//        {
//            return photon_passes_cuts[0];
//        }, {"photons_pass_cuts"});
//
//        auto dilepton_and_photon = photon_pt_cut
//        .Define("reconstructed_mass",[&](RVec<Electron>& di_electrons, Photon& merged_photon)
//        {
//            auto four_momentum = di_electrons[0].Vector() + di_electrons[1].Vector();
//
//            return (four_momentum + merged_photon.Vector()).M()/1e3;
//
//        }, {"di_electrons", "merged_photon"});
//
//        auto pSR = dilepton_and_photon.Filter(
//        [](double reconstructed_mass)
//        {
//            return (reconstructed_mass >= 110) && (reconstructed_mass <= 130);
//        }, {"reconstructed_mass"});
//
//        auto SR = pSR.Filter(
//        [](RVec<float>& Eratio)
//        {
//            return (!Any(Eratio <= 0.8));
//        }, {"photon_shower_shape_e_ratio"});
//
//        auto SR_ID = SR.Filter(
//        [&](Photon& merged_photon)
//        {
//            return merged_photon.photon_id;
//        }, {"merged_photon"});
//
//        if (counter < 2)
//        {
//            Totals.push_back(trigger_selection.Count());
//            Totals.push_back(ptCut.Count()); //preselection
//            Totals.push_back(failed_resolved.Count()); //failed_resolved
//            Totals.push_back(photon_pt_cut.Count()); //photon_pt_cut
//            Totals.push_back(pSR.Count()); //pSR
//            Totals.push_back(SR.Count()); //SR
//            Totals.push_back(SR_ID.Count()); //SR_ID
//        }
//        else
//        {
//            for (auto& mass_point: massPoints)
//            {
//                auto mass_point_trigger_selection = trigger_selection.Filter([&]
//                (float axion_mass)
//                {
//                    return (roundToOneDecimalPlace(axion_mass) == mass_point);
//
//                }, {"axion_masses"});
//
//                auto mass_point_ptCut = ptCut.Filter([&]
//                (float axion_mass)
//                {
//                    return (roundToOneDecimalPlace(axion_mass) == mass_point);
//
//                }, {"axion_masses"});
//
//                auto mass_point_failed_resolved = failed_resolved.Filter([&]
//                (float axion_mass)
//                {
//                    return (roundToOneDecimalPlace(axion_mass) == mass_point);
//
//                }, {"axion_masses"});
//
//                auto mass_point_photon_pt_cut = photon_pt_cut.Filter([&]
//                (float axion_mass)
//                {
//                    return (roundToOneDecimalPlace(axion_mass) == mass_point);
//
//                }, {"axion_masses"});
//
//                auto mass_point_pSR = pSR.Filter([&]
//                (float axion_mass)
//                {
//                    return (roundToOneDecimalPlace(axion_mass) == mass_point);
//
//                }, {"axion_masses"});
//
//                auto mass_point_SR = SR.Filter([&]
//                (float axion_mass)
//                {
//                    return (roundToOneDecimalPlace(axion_mass) == mass_point);
//
//                }, {"axion_masses"});
//
//                auto mass_point_SR_ID = SR_ID.Filter([&]
//                (float axion_mass)
//                {
//                    return (roundToOneDecimalPlace(axion_mass) == mass_point);
//
//                }, {"axion_masses"});
//
//                Totals.push_back(mass_point_trigger_selection.Count());
//                Totals.push_back(mass_point_ptCut.Count()); //preselection
//                Totals.push_back(mass_point_failed_resolved.Count()); //failed_resolved
//                Totals.push_back(mass_point_photon_pt_cut.Count()); //photon_pt_cut
//                Totals.push_back(mass_point_pSR.Count()); //pSR
//                Totals.push_back(mass_point_SR.Count()); //SR
//                Totals.push_back(mass_point_SR_ID.Count()); //SR_ID
//            }
//        }
//
//        counter++;
//    }
//
//    ROOT::RDF::RunGraphs(Totals);
//
////    0    1    2    3    4    5    6          ma1
////    7    8    9    10   11   12   13         ma
////    14   15   16   17   18   19   20         displaced_axion_1
////    21   22   23   24   25   26   27         displaced_axion_2
////    28   29   30   31   32   33   34         displaced_axion_3
////    35   36   37   38   39   40   41         displaced_axion_4
////    42   43   44   45   46   47   48         displaced_axion_5
////    49   50   51   52   53   54   55         displaced_axion_6
////    56   57   58   59   60   61   62         displaced_axion_7
////    63   64   65   66   67   68   69         displaced_axion_8
////    70   71   72   73   74   75   76         displaced_axion_9
////    77   78   79   80   81   82   83         displaced_axion_10
//
//    std::cout << R"--(\section*{Table 19 Prompt and Displaced Signal Samples})--" << '\n';
//    std::cout << R"--(\hspace{-3cm}\scalebox{0.55}{)--" << '\n';
//    std::cout << R"--(\begin{tabular}{|c|c|c|c|c|c|c|})--" << '\n';
//    std::cout << R"--(\hline)--" << '\n';
//    std::cout << R"--({} & pass preselection & failed resolved category & photon $p_T$ cut & pSR & SR & SR-ID \\ \hline)--" << '\n';
//
//
//    for (int i = 0, j = 0; (i <= 77 && j < 12); i += 7, j++)
//    {
//        std::cout << prefixes[j];
//        std::cout
//        << " & " << std::setprecision(2) << std::fixed << *Totals[i+1].GetResultPtr<ULong64_t>()
//        << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(*Totals[i+1].GetResultPtr<ULong64_t>())
//        << " & " << std::setprecision(2) << std::fixed << *Totals[i+2].GetResultPtr<ULong64_t>()
//        << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(*Totals[i+2].GetResultPtr<ULong64_t>())
//        << " & " << std::setprecision(2) << std::fixed << *Totals[i+3].GetResultPtr<ULong64_t>()
//        << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(*Totals[i+3].GetResultPtr<ULong64_t>())
//        << " & " << std::setprecision(2) << std::fixed << *Totals[i+4].GetResultPtr<ULong64_t>()
//        << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(*Totals[i+4].GetResultPtr<ULong64_t>())
//        << " & " << std::setprecision(2) << std::fixed << *Totals[i+5].GetResultPtr<ULong64_t>()
//        << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(*Totals[i+5].GetResultPtr<ULong64_t>())
//        << " & " << std::setprecision(2) << std::fixed << *Totals[i+6].GetResultPtr<ULong64_t>()
//        << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(*Totals[i+6].GetResultPtr<ULong64_t>())
//        << R"--( \\ \hline)--" << '\n';
//
//    }
//
//    std::cout << R"--(\end{tabular}})--" << '\n';
//
//    std::cout << "\n\n\n";
//
//}

void DataBackgroundComparison()
{
    auto start_time = Clock::now();
    fig26();
//    fig27();
    fig28();
    fig41();
    fig48();
    fig59();
//    Table9();
//    Table10();
//    Table16();
//    Table19();
//    Table16_Displaced_Axions();
//    Table19_Displaced_Axions();
    
    auto end_time = Clock::now();
    std::cout << "Time difference: "
       << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count()/1e9 << " seconds" << std::endl;
}

int main()
{
    DataBackgroundComparison();
}

