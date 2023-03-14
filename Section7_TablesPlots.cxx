#include <iostream>
#include <vector>
#include <chrono>
#include <algorithm>
#include <string>
#include <sstream>
#include <iomanip>
#include <array>
#include <cstdlib>
#include <unordered_set>

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

void Table21()
{
    Event::systematics =
    {
//        "PH_EFF_ISO_Uncertainty",
//        "PH_EFF_ISO_Uncertainty",
//        "EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR",
//        "PRW_DATASF",
//        "MUON_EFF_RECO_SYS",
//        "MUON_EFF_ISO_SYS",
//        "MUON_EFF_TrigSystUncertainty",
//        "EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR",
//        "MUON_EFF_TrigStatUncertainty",
//        "MUON_EFF_RECO_STAT",
//        "MUON_EFF_TTVA_STAT",
//        "EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR",
//        "EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR",
//        "MUON_EFF_TTVA_SYS",
//        "MUON_EFF_ISO_STAT",
//        "MUON_SAGITTA_RHO",
//        "EG_RESOLUTION_ALL",
//        "EG_SCALE_ALL",
//        "MUON_MS",
//        "MUON_ID",
//        "EL_EFF_TriggerEff_TOTAL_1NPCOR_PLUS_UNCOR",
//        "MUON_SAGITTA_RESBIAS",
//        "MUON_SCALE",

//        "PH_EFF_ISO_Uncertainty__1down",
//        "EG_SCALE_ALL__1down",
//        "PH_EFF_ID_Uncertainty__1down",
//        "PH_EFF_TRIGGER_Uncertainty__1down",
        "EG_RESOLUTION_ALL__1up",
        "EG_SCALE_ALL__1up",
        "PH_EFF_ISO_Uncertainty__1up",
        "PH_EFF_ID_Uncertainty__1up",
        "PH_EFF_TRIGGER_Uncertainty__1up",
//        "EG_RESOLUTION_ALL__1down",
    };

    std::vector<std::vector<std::string>> input_filenames =
    {
        //Z-gamma
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617070._000001.LGNTuple.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617064._000001.LGNTuple.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617074._000001.LGNTuple.root"},
        //data
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_data_test.root"},
        //signal
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root"},
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

    std::array<double,9> JetNumeratorSFs = {((139e15)*(1.9828e-9)*(0.821204)),((139e15)*(110.64e-12)*(0.69275)),((139e15)*(40.645e-12)*(0.615906)),((139e15)*(1.9817e-9)*(0.1136684)),((139e15)*(110.47e-12)*(0.1912956)),((139e15)*(40.674e-12)*(0.2326772)),((139e15)*(1.9819e-9)*(0.0656969)),((139e15)*(110.53e-12)*(0.1158741)),((139e15)*(40.68e-12)*(0.1535215))};
    std::array<double,3> SFs = {((139e15)*(.871e-12)),((139e15)*(.199e-12)), ((139e15)*(.0345e-15))}; //numerators for Z-gamma bkg

    std::vector<std::string> prefixes = { R"--(pty2\_9\_17)--", R"--(pty\_17\_myy\_0\_80)--", R"--(pty\_17\_myy\_80)--", "data", R"--($\text{Sig } m_{A}$ = 5 GeV)--", R"--($\text{Sig } m_{A}$ = 1 GeV)--", R"--(Zee\_lightJet\_0-70)--", R"--(Zee\_lightJet\_70-140)--", R"--(Zee\_lightJet\_140-280)--", R"--(Zee\_cJet\_0-70)--", R"--(Zee\_cJet\_70-140)--", R"--(Zee\_cJet\_140-280)--", R"--(Zee\_bJet\_0-70)--", R"--(Zee\_bJet\_70-140)--", R"--(Zee\_bJet\_140-280)--"};

    std::vector<ROOT::RDF::RResultHandle> Totals;
//    std::vector<RResultMap<ULong64_t>> resultmaps;
    std::vector<RResultMap<float>> resultmaps;
    std::vector<RResultMap<float>> PH_EFF_resultmaps;
//    std::vector<ROOT::RDF::RResultPtr<float>> GeneratorWeightCounts;

    std::stringstream ss;

    for (auto& i: input_filenames)
    {
        SchottDataFrame df(MakeRDF(i,8));
//        std::cout << *df.Count() << '\n';
//        df.Describe().Print();
//        exit(1);

        auto EventWeight = df.Define("EventWeight",
        [](const RVec<float>& ei_event_weights_generator)
        {
            return  ((ei_event_weights_generator[0]) ? 1 / ei_event_weights_generator[0] : 1);

        }, {"ei_event_weights_generator"});

        auto two_leptons = df.Filter(
        [](const RVec<Muon>& muons, RVec<AbstractParticle> electrons)
        {
            electrons.erase(std::remove_if(electrons.begin(),electrons.end(),
            [](AbstractParticle& ep)
            {
                return (!((ep.electron_pt/1e3 > 20) && (abs(ep.electron_eta) < 2.37) &&
                          (!((1.37 < abs(ep.electron_eta)) && (abs(ep.electron_eta) < 1.52)))
                          && (ep.electron_id_medium == 1)));

            }), electrons.end());

            return (electrons.size()==2 && muons.empty());

        }, {"muons", "abstract_electrons"});

        auto opp_charge = two_leptons.Define("di_electrons",
        [](RVec<AbstractParticle> electrons)
        {
            electrons.erase(std::remove_if(electrons.begin(),electrons.end(),
            [](AbstractParticle& ep)
            {
                return (!((ep.electron_pt/1e3 > 20) && (abs(ep.electron_eta) < 2.37) &&
                (!((1.37 < abs(ep.electron_eta)) && (abs(ep.electron_eta) < 1.52)))
                && (ep.electron_id_medium == 1)));

            }), electrons.end());

            return electrons;

        },{"abstract_electrons"})
        .Filter([](const RVec<AbstractParticle>& electrons)
        {
            return (electrons[0].electron_charge*electrons[1].electron_charge < 0);

        }, {"di_electrons"});

        auto leadingPt = opp_charge.Filter([](const RVec<AbstractParticle>& electrons)
        {
            return ((electrons[0].electron_pt > 20e3 && electrons[1].electron_pt > 27e3) || (electrons[1].electron_pt > 20e3 && electrons[0].electron_pt > 27e3));
        }, {"di_electrons"});

        auto deltaR = leadingPt.Filter([] (const RVec<AbstractParticle>& electrons)
        {
            return (DeltaR(electrons[0].ElectronVector(), electrons[1].ElectronVector()) > 0.01);
        }, {"di_electrons"});

        auto mass = deltaR.Filter([] (const RVec<AbstractParticle>& electrons)
        {
            auto mass = (electrons[0].ElectronVector() + electrons[1].ElectronVector()).M()/1e3;
            return ((mass >= 81) && (mass <= 101));
        }, {"di_electrons"});

        auto ptCut = mass.Filter([] (const RVec<AbstractParticle>& electrons)
        {
            auto pT = (electrons[0].ElectronVector() + electrons[1].ElectronVector()).Pt()/1e3;
            return pT > 10;
        }, {"di_electrons"});

        auto photon_passes_cuts = ptCut.Define("photons_pass_cuts",
        [&](RVec<AbstractParticle> photons)
        {
            photons.erase(std::remove_if(photons.begin(),photons.end(),
            [](AbstractParticle& x)
            {
              return ((abs(x.photon_eta) >= 2.37) || (abs(x.photon_eta) > 1.37 && abs(x.photon_eta) < 1.52) || (!x.photon_id_loose));

            }), photons.end());

            return photons;
        }, {"abstract_photons"});

        auto merged_reco_photons_matched = photon_passes_cuts.Filter(
        [&](const RVec<AbstractParticle>& reco_photons_matched)
        {
//            RVec<AbstractParticle> reco_photons_matched = reco_photons_test;
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
                delta_r = DeltaR(reco_photons_matched[combs[0][i]].PhotonVector(), reco_photons_matched[combs[1][i]].PhotonVector());
                m = (reco_photons_matched[combs[0][i]].PhotonVector() + reco_photons_matched[combs[1][i]].PhotonVector()).M();
                pt = (reco_photons_matched[combs[0][i]].PhotonVector() + reco_photons_matched[combs[1][i]].PhotonVector()).Pt();
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
        [&](const RVec<AbstractParticle>& reco_photons_matched)
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
        .Define("reconstructed_mass",[&](const RVec<AbstractParticle>& di_electrons, const AbstractParticle& merged_photon)
        {
            auto four_momentum = di_electrons[0].ElectronVector() + di_electrons[1].ElectronVector();

            return (four_momentum + merged_photon.PhotonVector()).M()/1e3;

        }, {"di_electrons", "merged_photon"});

        auto pSB = dilepton_and_photon.Filter(
        [](const double reconstructed_mass)
        {
            return (reconstructed_mass < 110) || (reconstructed_mass > 130);
        }, {"reconstructed_mass"});

        auto pSR = dilepton_and_photon.Filter(
        [](const double reconstructed_mass)
        {
            return (reconstructed_mass >= 110) && (reconstructed_mass <= 130);
        }, {"reconstructed_mass"});

        auto SB = pSB.Filter(
        [](const RVec<float>& Eratio)
        {
            return (!Any(Eratio <= 0.8));
        }, {"photon_shower_shape_e_ratio"});

        auto SR = pSR.Filter(
        [](const RVec<float>& Eratio)
        {
            return (!Any(Eratio <= 0.8));
        }, {"photon_shower_shape_e_ratio"});

        auto totEventWeight = merged_reco_photons_matched
        .Define("totEventWeight", [](RVec<float> photon_id_eff, RVec<float> photon_iso_eff, RVec<float> photon_trg_eff/*, RVec<float> ei_event_weights_generator*/)
        {
            auto ResizeVal = std::max({photon_id_eff.size(), photon_iso_eff.size(), photon_trg_eff.size()});
            photon_id_eff.resize(ResizeVal,1);
            photon_iso_eff.resize(ResizeVal,1);
            photon_trg_eff.resize(ResizeVal,1);

            return photon_id_eff*photon_iso_eff*photon_trg_eff;//*ei_event_weights_generator[0];

        }, {"photon_id_eff", "photon_iso_eff", "photon_trg_eff",/* "ei_event_weights_generator"*/});

//        Totals.push_back(df.Count());
//        Totals.push_back(EventWeight.Sum<RVec<float>>("EventWeight"));
        Totals.push_back(EventWeight.Sum<float>("EventWeight"));
//        resultmaps.push_back(VariationsFor(merged_reco_photons_matched.Count()));
        resultmaps.push_back(VariationsFor(totEventWeight.Sum<RVec<float>>("totEventWeight")));
        PH_EFF_resultmaps.push_back(VariationsFor(totEventWeight.Sum<RVec<float>>("totEventWeight")));
    }

//EG_RESOLUTION_ALL__1down:           42864.7 ID 1.00524  ISO 1.02283  TRIG 1.04157
//EG_RESOLUTION_ALL__1up:             42751.5 ID 1.00524  ISO 1.02283  TRIG 1.04157
//EG_SCALE_ALL__1down:                42521.1 ID 1.00524  ISO 1.02283  TRIG 1.04157
//EG_SCALE_ALL__1up:                  43059.8 ID 1.00524  ISO 1.02283  TRIG 1.04157
//PH_EFF_ID_Uncertainty__1down:       42790.4 ID 0.992328 ISO 1.02283  TRIG 1.04157
//PH_EFF_ID_Uncertainty__1up:         42790.4 ID 1.01816  ISO 1.02283  TRIG 1.04157
//PH_EFF_ISO_Uncertainty__1down:      42790.4 ID 1.00524  ISO 1.00971  TRIG 1.04157
//PH_EFF_ISO_Uncertainty__1up:        42790.4 ID 1.00524  ISO 1.03595  TRIG 1.04157
//PH_EFF_TRIGGER_Uncertainty__1down:  42790.4 ID 1.00524  ISO 1.02283  TRIG 1.01646
//PH_EFF_TRIGGER_Uncertainty__1up:    42790.4 ID 1.00524  ISO 1.02283  TRIG 1.06668

    ROOT::RDF::RunGraphs(Totals); // running all computation nodes concurrently

//    for (auto& i: GeneratorWeightCounts)
//    {
//        std::cout << *i << '\n';
//    }
//    int count = 0;
//    for (auto& i: resultmaps)
//    {
//        std::cout << count++ << "\n==\n";
//        for (auto& j: i.GetKeys())
//        {
//            std::cout << j << '\n';
//        }
//        std::cout << '\n';
//    }
//          resultmaps
//          ----------
//    0    1    Z-gamma
//    2    3    Z-gamma
//    4    5    Z-gamma
//    6    7    data
//    8    9    signal
//    10   11   signal
//    12   13   Z-jets
//    14   15   Z-jets
//    16   17   Z-jets
//    18   19   Z-jets
//    20   21   Z-jets
//    22   23   Z-jets
//    24   25   Z-jets
//    26   27   Z-jets
//    28   29   Z-jets

    ss << R"--(\section*{Table 21})--" << '\n';
    ss << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
    ss << R"--(\setlength\extrarowheight{2pt}\renewcommand{\arraystretch}{1.5})--" << '\n';
    ss << R"--(\begin{tabular}{|c|c|c|c|c|c|})--" << '\n';
    ss << R"--(\hline)--" << '\n';

        ss << R"--(\multicolumn{6}{|c|}{\parbox{\linewidth}{\centering Merged Photon Category: Up Variations \\ (\% difference from nominal)}}\\[5 pt])--" << '\n';
    ss << R"--(\hline)--" << '\n';

    ss << R"--({Sample} & EG\_RESOLUTION\_ALL & EG\_SCALE\_ALL & PH\_EFF\_ISO\_Uncertainty & PH\_EFF\_ID\_Uncertainty & PH\_EFF\_TRIGGER\_Uncertainty \\ \hline)--" << '\n';

    double finalScaleVal;

    double ZgammaNominal = 0, ZgammaEG_RESOLUTION_ALL = 0,
    ZgammaEG_SCALE_ALL = 0, ZgammaPH_EFF_ISO_Uncertainty = 0,
    ZgammaPH_EFF_ID_Uncertainty = 0, ZgammaPH_EFF_TRIGGER_Uncertainty = 0,
    Zgamma_PH_EFF_Nominal = 0;

    double ZjetsNominal = 0, ZjetsEG_RESOLUTION_ALL = 0,
    ZjetsEG_SCALE_ALL = 0, ZjetsPH_EFF_ISO_Uncertainty = 0,
    ZjetsPH_EFF_ID_Uncertainty = 0, ZjetsPH_EFF_TRIGGER_Uncertainty = 0,
    Zjets_PH_EFF_Nominal = 0;

    double totbkgNominal = 0, totbkgEG_RESOLUTION_ALL = 0,
    totbkgEG_SCALE_ALL = 0, totbkgPH_EFF_ISO_Uncertainty = 0,
    totbkgPH_EFF_ID_Uncertainty = 0, totbkgPH_EFF_TRIGGER_Uncertainty = 0,
    totbkg_PH_EFF_Nominal = 0;

    for (auto i = 0; (i < PH_EFF_resultmaps.size()); i++)
    {
        double nominalVal = static_cast<double>(resultmaps[i]["nominal"]);
        double PH_EFF_nominalVal = static_cast<double>(PH_EFF_resultmaps[i]["nominal"]);

//        auto denominator = *Totals[i].GetResultPtr<ULong64_t>();
        auto denominator = *Totals[i].GetResultPtr<float>();

        if (i >= 0 && i <= 2) //Zgamma
        {
            finalScaleVal = SFs[i]/denominator;
            ZgammaNominal += finalScaleVal*nominalVal;
            Zgamma_PH_EFF_Nominal += finalScaleVal*PH_EFF_nominalVal;
            ZgammaEG_RESOLUTION_ALL += finalScaleVal*resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1up"];
            ZgammaEG_SCALE_ALL += finalScaleVal*resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1up"];
            ZgammaPH_EFF_ISO_Uncertainty += finalScaleVal*PH_EFF_resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1up"];
            ZgammaPH_EFF_ID_Uncertainty += finalScaleVal*PH_EFF_resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1up"];
            ZgammaPH_EFF_TRIGGER_Uncertainty += finalScaleVal*PH_EFF_resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up"];
        }

        else if (i >= 6)
        {
            finalScaleVal = JetNumeratorSFs[i-6]/denominator;
            ZjetsNominal += finalScaleVal*nominalVal;
            Zjets_PH_EFF_Nominal += finalScaleVal*PH_EFF_nominalVal;
            ZjetsEG_RESOLUTION_ALL += finalScaleVal*resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1up"];
            ZjetsEG_SCALE_ALL += finalScaleVal*resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1up"];
            ZjetsPH_EFF_ISO_Uncertainty += finalScaleVal*PH_EFF_resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1up"];
            ZjetsPH_EFF_ID_Uncertainty += finalScaleVal*PH_EFF_resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1up"];
            ZjetsPH_EFF_TRIGGER_Uncertainty += finalScaleVal*PH_EFF_resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up"];
        }

        ss << prefixes[i] << " & ";
        ss << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1up"]) ? ((resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1up"]-nominalVal)/nominalVal)*100.0 : 0.0)
        << " & " << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1up"]) ? ((resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1up"]-nominalVal)/nominalVal)*100.0 : 0.0)
        << " & " << std::setprecision(4) << std::fixed
        << ( (PH_EFF_nominalVal && PH_EFF_resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1up"]) ? ((PH_EFF_resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1up"]-PH_EFF_nominalVal)/PH_EFF_nominalVal)*100.0 : 0.0)
        << " & " << std::setprecision(4) << std::fixed
        << ( (PH_EFF_nominalVal && PH_EFF_resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1up"]) ? ((PH_EFF_resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1up"]-PH_EFF_nominalVal)/PH_EFF_nominalVal)*100.0 : 0.0)
        << " & " << std::setprecision(4) << std::fixed
        << ( (PH_EFF_nominalVal && PH_EFF_resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up"]) ? ((PH_EFF_resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up"]-PH_EFF_nominalVal)/PH_EFF_nominalVal)*100.0 : 0.0)
        << R"--( \\ \hline)--" << '\n';

        if (!PH_EFF_resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1up"])
        {
            std::cout << i << ": photon_iso_eff:PH_EFF_ISO_Uncertainty__1up not found\n";
        }
        else
        {
            std::cout << i << ':' << ' ' << PH_EFF_resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1up"] << ' ' << PH_EFF_nominalVal << '\n';
        }
        if (!PH_EFF_resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1up"])
        {
            std::cout << i << ": photon_id_eff:PH_EFF_ID_Uncertainty__1up not found\n";
        }
        else
        {
            std::cout << i << ':' << ' ' << PH_EFF_resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1up"] << ' ' << PH_EFF_nominalVal << '\n';
        }
        if (!PH_EFF_resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up"])
        {
            std::cout << i << ": photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up not found\n";
        }
        else
        {
            std::cout << i << ':' << ' ' << PH_EFF_resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up"] << ' ' << PH_EFF_nominalVal << '\n';
        }
    }

    ss << R"--(Total $Z\gamma$ & )--";
    ss << std::setprecision(4) << std::fixed
    << ( (ZgammaNominal && ZgammaEG_RESOLUTION_ALL) ? ((ZgammaEG_RESOLUTION_ALL-ZgammaNominal)/ZgammaNominal)*100.0 : 0.0)
    << " & " << std::setprecision(4) << std::fixed
    << ( (ZgammaNominal && ZgammaEG_SCALE_ALL) ? ((ZgammaEG_SCALE_ALL-ZgammaNominal)/ZgammaNominal)*100.0 : 0.0)
    << " & " << std::setprecision(4) << std::fixed
    << ( (Zgamma_PH_EFF_Nominal && ZgammaPH_EFF_ISO_Uncertainty) ? ((ZgammaPH_EFF_ISO_Uncertainty-Zgamma_PH_EFF_Nominal)/Zgamma_PH_EFF_Nominal)*100.0 : 0.0)
    << " & " << std::setprecision(4) << std::fixed
    << ( (Zgamma_PH_EFF_Nominal && ZgammaPH_EFF_ID_Uncertainty) ? ((ZgammaPH_EFF_ID_Uncertainty-Zgamma_PH_EFF_Nominal)/Zgamma_PH_EFF_Nominal)*100.0 : 0.0)
    << " & " << std::setprecision(4) << std::fixed
    << ( (Zgamma_PH_EFF_Nominal && ZgammaPH_EFF_TRIGGER_Uncertainty) ? ((ZgammaPH_EFF_TRIGGER_Uncertainty-Zgamma_PH_EFF_Nominal)/Zgamma_PH_EFF_Nominal)*100.0 : 0.0)
    << R"--( \\ \hline)--" << '\n';

    ss << R"--(Total $Z$ jets & )--";
    ss << std::setprecision(4) << std::fixed
    << ( (ZjetsNominal && ZjetsEG_RESOLUTION_ALL) ? ((ZjetsEG_RESOLUTION_ALL-ZjetsNominal)/ZjetsNominal)*100.0 : 0.0)
    << " & " << std::setprecision(4) << std::fixed
    << ( (ZjetsNominal && ZjetsEG_SCALE_ALL) ? ((ZjetsEG_SCALE_ALL-ZjetsNominal)/ZjetsNominal)*100.0 : 0.0)
    << " & " << std::setprecision(4) << std::fixed
    << ( (Zjets_PH_EFF_Nominal && ZjetsPH_EFF_ISO_Uncertainty) ? ((ZjetsPH_EFF_ISO_Uncertainty-Zjets_PH_EFF_Nominal)/Zjets_PH_EFF_Nominal)*100.0 : 0.0)
    << " & " << std::setprecision(4) << std::fixed
    << ( (Zjets_PH_EFF_Nominal && ZjetsPH_EFF_ID_Uncertainty) ? ((ZjetsPH_EFF_ID_Uncertainty-Zjets_PH_EFF_Nominal)/Zjets_PH_EFF_Nominal)*100.0 : 0.0)
    << " & " << std::setprecision(4) << std::fixed
    << ( (Zjets_PH_EFF_Nominal && ZjetsPH_EFF_TRIGGER_Uncertainty) ? ((ZjetsPH_EFF_TRIGGER_Uncertainty-Zjets_PH_EFF_Nominal)/Zjets_PH_EFF_Nominal)*100.0 : 0.0)
    << R"--( \\ \hline)--" << '\n';

    totbkgNominal = ZgammaNominal + ZjetsNominal;
    totbkg_PH_EFF_Nominal = Zgamma_PH_EFF_Nominal + Zjets_PH_EFF_Nominal;
    totbkgEG_RESOLUTION_ALL = ZgammaEG_RESOLUTION_ALL + ZjetsEG_RESOLUTION_ALL;
    totbkgEG_SCALE_ALL = ZgammaEG_SCALE_ALL + ZjetsEG_SCALE_ALL;
    totbkgPH_EFF_ISO_Uncertainty = ZgammaPH_EFF_ISO_Uncertainty + ZjetsPH_EFF_ISO_Uncertainty;
    totbkgPH_EFF_ID_Uncertainty = ZgammaPH_EFF_ID_Uncertainty + ZjetsPH_EFF_ID_Uncertainty;
    totbkgPH_EFF_TRIGGER_Uncertainty = ZgammaPH_EFF_TRIGGER_Uncertainty + ZjetsPH_EFF_TRIGGER_Uncertainty;

    ss << R"--(Total Bkg & )--";
    ss << std::setprecision(4) << std::fixed
    << ( (totbkgNominal && totbkgEG_RESOLUTION_ALL) ? ((totbkgEG_RESOLUTION_ALL-totbkgNominal)/totbkgNominal)*100.0 : 0.0)
    << " & " << std::setprecision(4) << std::fixed
    << ( (totbkgNominal && totbkgEG_SCALE_ALL) ? ((totbkgEG_SCALE_ALL-totbkgNominal)/totbkgNominal)*100.0 : 0.0)
    << " & " << std::setprecision(4) << std::fixed
    << ( (totbkg_PH_EFF_Nominal && totbkgPH_EFF_ISO_Uncertainty) ? ((totbkgPH_EFF_ISO_Uncertainty-totbkg_PH_EFF_Nominal)/totbkg_PH_EFF_Nominal)*100.0 : 0.0)
    << " & " << std::setprecision(4) << std::fixed
    << ( (totbkg_PH_EFF_Nominal && totbkgPH_EFF_ID_Uncertainty) ? ((totbkgPH_EFF_ID_Uncertainty-totbkg_PH_EFF_Nominal)/totbkg_PH_EFF_Nominal)*100.0 : 0.0)
    << " & " << std::setprecision(4) << std::fixed
    << ( (totbkg_PH_EFF_Nominal && totbkgPH_EFF_TRIGGER_Uncertainty) ? ((totbkgPH_EFF_TRIGGER_Uncertainty-totbkg_PH_EFF_Nominal)/totbkg_PH_EFF_Nominal)*100.0 : 0.0)
    << R"--( \\ \hline)--" << '\n';

    ss << R"--(\end{tabular}})--" << '\n';

    ss << "\n\n\n";

    std::cout << ss.str();
}

void Section7_TablesPlots()
{
    
    auto start_time = Clock::now();
    auto verbosity = ROOT::Experimental::RLogScopedVerbosity(ROOT::Detail::RDF::RDFLogChannel(), ROOT::Experimental::ELogLevel::kInfo);
    Table21();

    auto end_time = Clock::now();
    std::cout << "Time difference: "
       << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count()/1e9 << " seconds" << std::endl;
    
}


int main()
{
    Section7_TablesPlots();
}

//using Clock = std::chrono::high_resolution_clock;
//std::vector<double> times;
//using ROOT::RDF::Experimental::VariationsFor;
//std::vector<std::string> tags;
//for (int i = 0; i <= 10; i++){auto start_time = Clock::now(); ROOT::RDataFrame d(1e6); auto d1 = d.Define("x", [](){return ROOT::VecOps::RVec<double>(10, gRandom->Rndm());}); if (i==0){auto v1 = d1.Sum<ROOT::VecOps::RVec<double>>("x");std::cout << *v1 << '\n';auto end_time = Clock::now();times.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count()/1e9);continue;} tags.push_back(std::string(i,static_cast<char>(i+64))); auto v1 = d1.Vary("x", [&]{ROOT::VecOps::RVec<ROOT::VecOps::RVec<double>> x;for (int j = 1; j <= i; j++) {x.push_back(ROOT::VecOps::RVec<double>(10, gRandom->Rndm()));} return x;}, {}, tags).Sum<ROOT::VecOps::RVec<double>>("x"); auto sums = VariationsFor(v1);for (auto& i: sums.GetKeys()){std::cout << sums[i] << '\t';} std::cout << '\n'; auto end_time = Clock::now(); times.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count()/1e9);}
//for (int i = 0; i <= 10; i++){if (i==0){std::cout << std::setw(25) << std::fixed << "# of Variations" << std::setw(25) << std::fixed << "Time (s)" << '\n';}std::cout << std::setw(25) << std::fixed << i << std::setw(25) << std::fixed << times[i] << '\n';}



