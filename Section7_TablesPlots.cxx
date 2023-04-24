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
#include <functional>

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

float roundToOneDecimalPlace(float num) {
    std::stringstream stream;
    stream << std::fixed << std::setprecision(1) << num;
    float rounded_num = std::stof(stream.str());
    return rounded_num;
}

//TODO: Fix this file after finishing Inner Detector section
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
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root"},
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
    ZgammaPH_EFF_ID_Uncertainty = 0, ZgammaPH_EFF_TRIGGER_Uncertainty = 0;

    double ZjetsNominal = 0, ZjetsEG_RESOLUTION_ALL = 0,
    ZjetsEG_SCALE_ALL = 0, ZjetsPH_EFF_ISO_Uncertainty = 0,
    ZjetsPH_EFF_ID_Uncertainty = 0, ZjetsPH_EFF_TRIGGER_Uncertainty = 0;

    double totbkgNominal = 0, totbkgEG_RESOLUTION_ALL = 0,
    totbkgEG_SCALE_ALL = 0, totbkgPH_EFF_ISO_Uncertainty = 0,
    totbkgPH_EFF_ID_Uncertainty = 0, totbkgPH_EFF_TRIGGER_Uncertainty = 0;

    for (auto i = 0; (i < resultmaps.size()); i++)
    {
        double nominalVal = static_cast<double>(resultmaps[i]["nominal"]);

//        auto denominator = *Totals[i].GetResultPtr<ULong64_t>();
        auto denominator = *Totals[i].GetResultPtr<float>();

        if (i >= 0 && i <= 2) //Zgamma
        {
            finalScaleVal = SFs[i]/denominator;
            ZgammaNominal += finalScaleVal*nominalVal;

            ZgammaEG_RESOLUTION_ALL += finalScaleVal*resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1up"];
            ZgammaEG_SCALE_ALL += finalScaleVal*resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1up"];
            ZgammaPH_EFF_ISO_Uncertainty += finalScaleVal*resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1up"];
            ZgammaPH_EFF_ID_Uncertainty += finalScaleVal*resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1up"];
            ZgammaPH_EFF_TRIGGER_Uncertainty += finalScaleVal*resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up"];
        }

        else if (i >= 6)
        {
            finalScaleVal = JetNumeratorSFs[i-6]/denominator;
            ZjetsNominal += finalScaleVal*nominalVal;

            ZjetsEG_RESOLUTION_ALL += finalScaleVal*resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1up"];
            ZjetsEG_SCALE_ALL += finalScaleVal*resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1up"];
            ZjetsPH_EFF_ISO_Uncertainty += finalScaleVal*resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1up"];
            ZjetsPH_EFF_ID_Uncertainty += finalScaleVal*resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1up"];
            ZjetsPH_EFF_TRIGGER_Uncertainty += finalScaleVal*resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up"];
        }

        ss << prefixes[i] << " & ";
        ss << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1up"]) ? ((resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1up"]-nominalVal)/nominalVal)*100.0 : 0.0)
        << " & " << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1up"]) ? ((resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1up"]-nominalVal)/nominalVal)*100.0 : 0.0)
        << " & " << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1up"]) ? ((resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1up"]-nominalVal)/nominalVal)*100.0 : 0.0)
        << " & " << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1up"]) ? ((resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1up"]-nominalVal)/nominalVal)*100.0 : 0.0)
        << " & " << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up"]) ? ((resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up"]-nominalVal)/nominalVal)*100.0 : 0.0)
        << R"--( \\ \hline)--" << '\n';

        if (!resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1up"])
        {
            std::cout << i << ": photon_iso_eff:PH_EFF_ISO_Uncertainty__1up not found\n";
        }
        else
        {
            std::cout << i << ':' << ' ' << resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1up"] << ' ' << nominalVal << '\n';
        }
        if (!resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1up"])
        {
            std::cout << i << ": photon_id_eff:PH_EFF_ID_Uncertainty__1up not found\n";
        }
        else
        {
            std::cout << i << ':' << ' ' << resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1up"] << ' ' << nominalVal << '\n';
        }
        if (!resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up"])
        {
            std::cout << i << ": photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up not found\n";
        }
        else
        {
            std::cout << i << ':' << ' ' << resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up"] << ' ' << nominalVal << '\n';
        }
    }

    ss << R"--(Total $Z\gamma$ & )--";
    ss << std::setprecision(4) << std::fixed
    << ( (ZgammaNominal && ZgammaEG_RESOLUTION_ALL) ? ((ZgammaEG_RESOLUTION_ALL-ZgammaNominal)/ZgammaNominal)*100.0 : 0.0)
    << " & " << std::setprecision(4) << std::fixed
    << ( (ZgammaNominal && ZgammaEG_SCALE_ALL) ? ((ZgammaEG_SCALE_ALL-ZgammaNominal)/ZgammaNominal)*100.0 : 0.0)
    << " & " << std::setprecision(4) << std::fixed
    << ( (ZgammaNominal && ZgammaPH_EFF_ISO_Uncertainty) ? ((ZgammaPH_EFF_ISO_Uncertainty-ZgammaNominal)/ZgammaNominal)*100.0 : 0.0)
    << " & " << std::setprecision(4) << std::fixed
    << ( (ZgammaNominal && ZgammaPH_EFF_ID_Uncertainty) ? ((ZgammaPH_EFF_ID_Uncertainty-ZgammaNominal)/ZgammaNominal)*100.0 : 0.0)
    << " & " << std::setprecision(4) << std::fixed
    << ( (ZgammaNominal && ZgammaPH_EFF_TRIGGER_Uncertainty) ? ((ZgammaPH_EFF_TRIGGER_Uncertainty-ZgammaNominal)/ZgammaNominal)*100.0 : 0.0)
    << R"--( \\ \hline)--" << '\n';

    ss << R"--(Total $Z$ jets & )--";
    ss << std::setprecision(4) << std::fixed
    << ( (ZjetsNominal && ZjetsEG_RESOLUTION_ALL) ? ((ZjetsEG_RESOLUTION_ALL-ZjetsNominal)/ZjetsNominal)*100.0 : 0.0)
    << " & " << std::setprecision(4) << std::fixed
    << ( (ZjetsNominal && ZjetsEG_SCALE_ALL) ? ((ZjetsEG_SCALE_ALL-ZjetsNominal)/ZjetsNominal)*100.0 : 0.0)
    << " & " << std::setprecision(4) << std::fixed
    << ( (ZjetsNominal && ZjetsPH_EFF_ISO_Uncertainty) ? ((ZjetsPH_EFF_ISO_Uncertainty-ZjetsNominal)/ZjetsNominal)*100.0 : 0.0)
    << " & " << std::setprecision(4) << std::fixed
    << ( (ZjetsNominal && ZjetsPH_EFF_ID_Uncertainty) ? ((ZjetsPH_EFF_ID_Uncertainty-ZjetsNominal)/ZjetsNominal)*100.0 : 0.0)
    << " & " << std::setprecision(4) << std::fixed
    << ( (ZjetsNominal && ZjetsPH_EFF_TRIGGER_Uncertainty) ? ((ZjetsPH_EFF_TRIGGER_Uncertainty-ZjetsNominal)/ZjetsNominal)*100.0 : 0.0)
    << R"--( \\ \hline)--" << '\n';

    totbkgNominal = ZgammaNominal + ZjetsNominal;

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
    << ( (totbkgNominal && totbkgPH_EFF_ISO_Uncertainty) ? ((totbkgPH_EFF_ISO_Uncertainty-totbkgNominal)/totbkgNominal)*100.0 : 0.0)
    << " & " << std::setprecision(4) << std::fixed
    << ( (totbkgNominal && totbkgPH_EFF_ID_Uncertainty) ? ((totbkgPH_EFF_ID_Uncertainty-totbkgNominal)/totbkgNominal)*100.0 : 0.0)
    << " & " << std::setprecision(4) << std::fixed
    << ( (totbkgNominal && totbkgPH_EFF_TRIGGER_Uncertainty) ? ((totbkgPH_EFF_TRIGGER_Uncertainty-totbkgNominal)/totbkgNominal)*100.0 : 0.0)
    << R"--( \\ \hline)--" << '\n';

    ss << R"--(\end{tabular}})--" << '\n';

    ss << "\n\n\n";

    std::cout << ss.str();
}

void Table22()
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
//        "EG_RESOLUTION_ALL__1up",
//        "EG_SCALE_ALL__1up",
//        "PH_EFF_ISO_Uncertainty__1up",
//        "PH_EFF_ID_Uncertainty__1up",
//        "PH_EFF_TRIGGER_Uncertainty__1up",
        "EG_RESOLUTION_ALL__1down",
        "EG_SCALE_ALL__1down",
        "PH_EFF_ISO_Uncertainty__1down",
        "PH_EFF_ID_Uncertainty__1down",
        "PH_EFF_TRIGGER_Uncertainty__1down",
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
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root"},
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

    ss << R"--(\section*{Table 22})--" << '\n';
    ss << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
    ss << R"--(\setlength\extrarowheight{2pt}\renewcommand{\arraystretch}{1.5})--" << '\n';
    ss << R"--(\begin{tabular}{|c|c|c|c|c|c|})--" << '\n';
    ss << R"--(\hline)--" << '\n';

        ss << R"--(\multicolumn{6}{|c|}{\parbox{\linewidth}{\centering Merged Photon Category: Down Variations \\ (\% difference from nominal)}}\\[5 pt])--" << '\n';
    ss << R"--(\hline)--" << '\n';

    ss << R"--({Sample} & EG\_RESOLUTION\_ALL & EG\_SCALE\_ALL & PH\_EFF\_ISO\_Uncertainty & PH\_EFF\_ID\_Uncertainty & PH\_EFF\_TRIGGER\_Uncertainty \\ \hline)--" << '\n';

    double finalScaleVal;

    double ZgammaNominal = 0, ZgammaEG_RESOLUTION_ALL = 0,
    ZgammaEG_SCALE_ALL = 0, ZgammaPH_EFF_ISO_Uncertainty = 0,
    ZgammaPH_EFF_ID_Uncertainty = 0, ZgammaPH_EFF_TRIGGER_Uncertainty = 0;

    double ZjetsNominal = 0, ZjetsEG_RESOLUTION_ALL = 0,
    ZjetsEG_SCALE_ALL = 0, ZjetsPH_EFF_ISO_Uncertainty = 0,
    ZjetsPH_EFF_ID_Uncertainty = 0, ZjetsPH_EFF_TRIGGER_Uncertainty = 0;

    double totbkgNominal = 0, totbkgEG_RESOLUTION_ALL = 0,
    totbkgEG_SCALE_ALL = 0, totbkgPH_EFF_ISO_Uncertainty = 0,
    totbkgPH_EFF_ID_Uncertainty = 0, totbkgPH_EFF_TRIGGER_Uncertainty = 0;

    for (auto i = 0; (i < resultmaps.size()); i++)
    {
        double nominalVal = static_cast<double>(resultmaps[i]["nominal"]);

//        auto denominator = *Totals[i].GetResultPtr<ULong64_t>();
        auto denominator = *Totals[i].GetResultPtr<float>();

        if (i >= 0 && i <= 2) //Zgamma
        {
            finalScaleVal = SFs[i]/denominator;
            ZgammaNominal += finalScaleVal*nominalVal;

            ZgammaEG_RESOLUTION_ALL += finalScaleVal*resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1down"];
            ZgammaEG_SCALE_ALL += finalScaleVal*resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1down"];
            ZgammaPH_EFF_ISO_Uncertainty += finalScaleVal*resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1down"];
            ZgammaPH_EFF_ID_Uncertainty += finalScaleVal*resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1down"];
            ZgammaPH_EFF_TRIGGER_Uncertainty += finalScaleVal*resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down"];
        }

        else if (i >= 6)
        {
            finalScaleVal = JetNumeratorSFs[i-6]/denominator;
            ZjetsNominal += finalScaleVal*nominalVal;

            ZjetsEG_RESOLUTION_ALL += finalScaleVal*resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1down"];
            ZjetsEG_SCALE_ALL += finalScaleVal*resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1down"];
            ZjetsPH_EFF_ISO_Uncertainty += finalScaleVal*resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1down"];
            ZjetsPH_EFF_ID_Uncertainty += finalScaleVal*resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1down"];
            ZjetsPH_EFF_TRIGGER_Uncertainty += finalScaleVal*resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down"];
        }

        ss << prefixes[i] << " & ";
        ss << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1down"]) ? ((resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1down"]-nominalVal)/nominalVal)*100.0 : 0.0)
        << " & " << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1down"]) ? ((resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1down"]-nominalVal)/nominalVal)*100.0 : 0.0)
        << " & " << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1down"]) ? ((resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1down"]-nominalVal)/nominalVal)*100.0 : 0.0)
        << " & " << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1down"]) ? ((resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1down"]-nominalVal)/nominalVal)*100.0 : 0.0)
        << " & " << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down"]) ? ((resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down"]-nominalVal)/nominalVal)*100.0 : 0.0)
        << R"--( \\ \hline)--" << '\n';

        if (!resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1down"])
        {
            std::cout << i << ": photon_iso_eff:PH_EFF_ISO_Uncertainty__1down not found\n";
        }
        else
        {
            std::cout << i << ':' << ' ' << resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1down"] << ' ' << nominalVal << '\n';
        }
        if (!resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1down"])
        {
            std::cout << i << ": photon_id_eff:PH_EFF_ID_Uncertainty__1down not found\n";
        }
        else
        {
            std::cout << i << ':' << ' ' << resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1down"] << ' ' << nominalVal << '\n';
        }
        if (!resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down"])
        {
            std::cout << i << ": photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down not found\n";
        }
        else
        {
            std::cout << i << ':' << ' ' << resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down"] << ' ' << nominalVal << '\n';
        }
    }

    ss << R"--(Total $Z\gamma$ & )--";
    ss << std::setprecision(4) << std::fixed
    << ( (ZgammaNominal && ZgammaEG_RESOLUTION_ALL) ? ((ZgammaEG_RESOLUTION_ALL-ZgammaNominal)/ZgammaNominal)*100.0 : 0.0)
    << " & " << std::setprecision(4) << std::fixed
    << ( (ZgammaNominal && ZgammaEG_SCALE_ALL) ? ((ZgammaEG_SCALE_ALL-ZgammaNominal)/ZgammaNominal)*100.0 : 0.0)
    << " & " << std::setprecision(4) << std::fixed
    << ( (ZgammaNominal && ZgammaPH_EFF_ISO_Uncertainty) ? ((ZgammaPH_EFF_ISO_Uncertainty-ZgammaNominal)/ZgammaNominal)*100.0 : 0.0)
    << " & " << std::setprecision(4) << std::fixed
    << ( (ZgammaNominal && ZgammaPH_EFF_ID_Uncertainty) ? ((ZgammaPH_EFF_ID_Uncertainty-ZgammaNominal)/ZgammaNominal)*100.0 : 0.0)
    << " & " << std::setprecision(4) << std::fixed
    << ( (ZgammaNominal && ZgammaPH_EFF_TRIGGER_Uncertainty) ? ((ZgammaPH_EFF_TRIGGER_Uncertainty-ZgammaNominal)/ZgammaNominal)*100.0 : 0.0)
    << R"--( \\ \hline)--" << '\n';

    ss << R"--(Total $Z$ jets & )--";
    ss << std::setprecision(4) << std::fixed
    << ( (ZjetsNominal && ZjetsEG_RESOLUTION_ALL) ? ((ZjetsEG_RESOLUTION_ALL-ZjetsNominal)/ZjetsNominal)*100.0 : 0.0)
    << " & " << std::setprecision(4) << std::fixed
    << ( (ZjetsNominal && ZjetsEG_SCALE_ALL) ? ((ZjetsEG_SCALE_ALL-ZjetsNominal)/ZjetsNominal)*100.0 : 0.0)
    << " & " << std::setprecision(4) << std::fixed
    << ( (ZjetsNominal && ZjetsPH_EFF_ISO_Uncertainty) ? ((ZjetsPH_EFF_ISO_Uncertainty-ZjetsNominal)/ZjetsNominal)*100.0 : 0.0)
    << " & " << std::setprecision(4) << std::fixed
    << ( (ZjetsNominal && ZjetsPH_EFF_ID_Uncertainty) ? ((ZjetsPH_EFF_ID_Uncertainty-ZjetsNominal)/ZjetsNominal)*100.0 : 0.0)
    << " & " << std::setprecision(4) << std::fixed
    << ( (ZjetsNominal && ZjetsPH_EFF_TRIGGER_Uncertainty) ? ((ZjetsPH_EFF_TRIGGER_Uncertainty-ZjetsNominal)/ZjetsNominal)*100.0 : 0.0)
    << R"--( \\ \hline)--" << '\n';

    totbkgNominal = ZgammaNominal + ZjetsNominal;

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
    << ( (totbkgNominal && totbkgPH_EFF_ISO_Uncertainty) ? ((totbkgPH_EFF_ISO_Uncertainty-totbkgNominal)/totbkgNominal)*100.0 : 0.0)
    << " & " << std::setprecision(4) << std::fixed
    << ( (totbkgNominal && totbkgPH_EFF_ID_Uncertainty) ? ((totbkgPH_EFF_ID_Uncertainty-totbkgNominal)/totbkgNominal)*100.0 : 0.0)
    << " & " << std::setprecision(4) << std::fixed
    << ( (totbkgNominal && totbkgPH_EFF_TRIGGER_Uncertainty) ? ((totbkgPH_EFF_TRIGGER_Uncertainty-totbkgNominal)/totbkgNominal)*100.0 : 0.0)
    << R"--( \\ \hline)--" << '\n';

    ss << R"--(\end{tabular}})--" << '\n';

    ss << "\n\n\n";

    std::cout << ss.str();
}

void Table23()
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
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root"},
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

        auto resolved_reco_photons_matched = photon_passes_cuts.Define("chosen_two",
        [](const RVec<AbstractParticle>& reco_photons_matched)
        {
            RVec<AbstractParticle> x;
            if (reco_photons_matched.size() < 2)
            {
                return x;
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
                    x = {reco_photons_matched[combs[0][i]], reco_photons_matched[combs[1][i]]};
                }
            }
            if (pt1 > 10e3 && pt2 > 10e3 && best_X > 0.96 && best_X < 1.2 && chosen_delta_r < 1.5)
            {
                return x;
            }
            x.clear();
            return x;
        }, {"photons_pass_cuts"}).Filter(
        [&](RVec<AbstractParticle>& reco_photons_matched)
        {
            return (reco_photons_matched.size()==2);
        }, {"chosen_two"}).Define("reconstructed_mass",
        [&](RVec<AbstractParticle>& diph)
        {
            return (diph[0].PhotonVector()+diph[1].PhotonVector()).M()/1e3;
        }, {"chosen_two"});

        auto SB = resolved_reco_photons_matched.Filter(
        [](const double reconstructed_mass)
        {
            return (reconstructed_mass < 110) || (reconstructed_mass > 140);
        }, {"reconstructed_mass"});

        auto SR = resolved_reco_photons_matched.Filter(
        [](const double reconstructed_mass)
        {
            return (reconstructed_mass >= 110) && (reconstructed_mass <= 140);
        }, {"reconstructed_mass"});

        auto totEventWeight = resolved_reco_photons_matched
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

    ss << R"--(\section*{Table 23})--" << '\n';
    ss << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
    ss << R"--(\setlength\extrarowheight{2pt}\renewcommand{\arraystretch}{1.5})--" << '\n';
    ss << R"--(\begin{tabular}{|c|c|c|c|c|c|})--" << '\n';
    ss << R"--(\hline)--" << '\n';

        ss << R"--(\multicolumn{6}{|c|}{\parbox{\linewidth}{\centering Resolved Photon Category: Up Variations \\ (\% difference from nominal)}}\\[5 pt])--" << '\n';
    ss << R"--(\hline)--" << '\n';

    ss << R"--({Sample} & EG\_RESOLUTION\_ALL & EG\_SCALE\_ALL & PH\_EFF\_ISO\_Uncertainty & PH\_EFF\_ID\_Uncertainty & PH\_EFF\_TRIGGER\_Uncertainty \\ \hline)--" << '\n';

    double finalScaleVal;

    double ZgammaNominal = 0, ZgammaEG_RESOLUTION_ALL = 0,
    ZgammaEG_SCALE_ALL = 0, ZgammaPH_EFF_ISO_Uncertainty = 0,
    ZgammaPH_EFF_ID_Uncertainty = 0, ZgammaPH_EFF_TRIGGER_Uncertainty = 0;

    double ZjetsNominal = 0, ZjetsEG_RESOLUTION_ALL = 0,
    ZjetsEG_SCALE_ALL = 0, ZjetsPH_EFF_ISO_Uncertainty = 0,
    ZjetsPH_EFF_ID_Uncertainty = 0, ZjetsPH_EFF_TRIGGER_Uncertainty = 0;

    double totbkgNominal = 0, totbkgEG_RESOLUTION_ALL = 0,
    totbkgEG_SCALE_ALL = 0, totbkgPH_EFF_ISO_Uncertainty = 0,
    totbkgPH_EFF_ID_Uncertainty = 0, totbkgPH_EFF_TRIGGER_Uncertainty = 0;

    for (auto i = 0; (i < resultmaps.size()); i++)
    {
        double nominalVal = static_cast<double>(resultmaps[i]["nominal"]);
//        auto denominator = *Totals[i].GetResultPtr<ULong64_t>();
        auto denominator = *Totals[i].GetResultPtr<float>();

        if (i >= 0 && i <= 2) //Zgamma
        {
            finalScaleVal = SFs[i]/denominator;
            ZgammaNominal += finalScaleVal*nominalVal;

            ZgammaEG_RESOLUTION_ALL += finalScaleVal*resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1up"];
            ZgammaEG_SCALE_ALL += finalScaleVal*resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1up"];
            ZgammaPH_EFF_ISO_Uncertainty += finalScaleVal*resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1up"];
            ZgammaPH_EFF_ID_Uncertainty += finalScaleVal*resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1up"];
            ZgammaPH_EFF_TRIGGER_Uncertainty += finalScaleVal*resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up"];
        }

        else if (i >= 6)
        {
            finalScaleVal = JetNumeratorSFs[i-6]/denominator;
            ZjetsNominal += finalScaleVal*nominalVal;

            ZjetsEG_RESOLUTION_ALL += finalScaleVal*resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1up"];
            ZjetsEG_SCALE_ALL += finalScaleVal*resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1up"];
            ZjetsPH_EFF_ISO_Uncertainty += finalScaleVal*resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1up"];
            ZjetsPH_EFF_ID_Uncertainty += finalScaleVal*resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1up"];
            ZjetsPH_EFF_TRIGGER_Uncertainty += finalScaleVal*resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up"];
        }

        ss << prefixes[i] << " & ";
        ss << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1up"]) ? ((resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1up"]-nominalVal)/nominalVal)*100.0 : 0.0)
        << " & " << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1up"]) ? ((resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1up"]-nominalVal)/nominalVal)*100.0 : 0.0)
        << " & " << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1up"]) ? ((resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1up"]-nominalVal)/nominalVal)*100.0 : 0.0)
        << " & " << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1up"]) ? ((resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1up"]-nominalVal)/nominalVal)*100.0 : 0.0)
        << " & " << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up"]) ? ((resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up"]-nominalVal)/nominalVal)*100.0 : 0.0)
        << R"--( \\ \hline)--" << '\n';

        if (!resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1up"])
        {
            std::cout << i << ": photon_iso_eff:PH_EFF_ISO_Uncertainty__1up not found\n";
        }
        else
        {
            std::cout << i << ':' << ' ' << resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1up"] << ' ' << nominalVal << '\n';
        }
        if (!resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1up"])
        {
            std::cout << i << ": photon_id_eff:PH_EFF_ID_Uncertainty__1up not found\n";
        }
        else
        {
            std::cout << i << ':' << ' ' << resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1up"] << ' ' << nominalVal << '\n';
        }
        if (!resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up"])
        {
            std::cout << i << ": photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up not found\n";
        }
        else
        {
            std::cout << i << ':' << ' ' << resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up"] << ' ' << nominalVal << '\n';
        }
    }

    ss << R"--(Total $Z\gamma$ & )--";
    ss << std::setprecision(4) << std::fixed
    << ( (ZgammaNominal && ZgammaEG_RESOLUTION_ALL) ? ((ZgammaEG_RESOLUTION_ALL-ZgammaNominal)/ZgammaNominal)*100.0 : 0.0)
    << " & " << std::setprecision(4) << std::fixed
    << ( (ZgammaNominal && ZgammaEG_SCALE_ALL) ? ((ZgammaEG_SCALE_ALL-ZgammaNominal)/ZgammaNominal)*100.0 : 0.0)
    << " & " << std::setprecision(4) << std::fixed
    << ( (ZgammaNominal && ZgammaPH_EFF_ISO_Uncertainty) ? ((ZgammaPH_EFF_ISO_Uncertainty-ZgammaNominal)/ZgammaNominal)*100.0 : 0.0)
    << " & " << std::setprecision(4) << std::fixed
    << ( (ZgammaNominal && ZgammaPH_EFF_ID_Uncertainty) ? ((ZgammaPH_EFF_ID_Uncertainty-ZgammaNominal)/ZgammaNominal)*100.0 : 0.0)
    << " & " << std::setprecision(4) << std::fixed
    << ( (ZgammaNominal && ZgammaPH_EFF_TRIGGER_Uncertainty) ? ((ZgammaPH_EFF_TRIGGER_Uncertainty-ZgammaNominal)/ZgammaNominal)*100.0 : 0.0)
    << R"--( \\ \hline)--" << '\n';

    ss << R"--(Total $Z$ jets & )--";
    ss << std::setprecision(4) << std::fixed
    << ( (ZjetsNominal && ZjetsEG_RESOLUTION_ALL) ? ((ZjetsEG_RESOLUTION_ALL-ZjetsNominal)/ZjetsNominal)*100.0 : 0.0)
    << " & " << std::setprecision(4) << std::fixed
    << ( (ZjetsNominal && ZjetsEG_SCALE_ALL) ? ((ZjetsEG_SCALE_ALL-ZjetsNominal)/ZjetsNominal)*100.0 : 0.0)
    << " & " << std::setprecision(4) << std::fixed
    << ( (ZjetsNominal && ZjetsPH_EFF_ISO_Uncertainty) ? ((ZjetsPH_EFF_ISO_Uncertainty-ZjetsNominal)/ZjetsNominal)*100.0 : 0.0)
    << " & " << std::setprecision(4) << std::fixed
    << ( (ZjetsNominal && ZjetsPH_EFF_ID_Uncertainty) ? ((ZjetsPH_EFF_ID_Uncertainty-ZjetsNominal)/ZjetsNominal)*100.0 : 0.0)
    << " & " << std::setprecision(4) << std::fixed
    << ( (ZjetsNominal && ZjetsPH_EFF_TRIGGER_Uncertainty) ? ((ZjetsPH_EFF_TRIGGER_Uncertainty-ZjetsNominal)/ZjetsNominal)*100.0 : 0.0)
    << R"--( \\ \hline)--" << '\n';

    totbkgNominal = ZgammaNominal + ZjetsNominal;

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
    << ( (totbkgNominal && totbkgPH_EFF_ISO_Uncertainty) ? ((totbkgPH_EFF_ISO_Uncertainty-totbkgNominal)/totbkgNominal)*100.0 : 0.0)
    << " & " << std::setprecision(4) << std::fixed
    << ( (totbkgNominal && totbkgPH_EFF_ID_Uncertainty) ? ((totbkgPH_EFF_ID_Uncertainty-totbkgNominal)/totbkgNominal)*100.0 : 0.0)
    << " & " << std::setprecision(4) << std::fixed
    << ( (totbkgNominal && totbkgPH_EFF_TRIGGER_Uncertainty) ? ((totbkgPH_EFF_TRIGGER_Uncertainty-totbkgNominal)/totbkgNominal)*100.0 : 0.0)
    << R"--( \\ \hline)--" << '\n';

    ss << R"--(\end{tabular}})--" << '\n';

    ss << "\n\n\n";

    std::cout << ss.str();
}

void Table24()
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
//        "EG_RESOLUTION_ALL__1up",
//        "EG_SCALE_ALL__1up",
//        "PH_EFF_ISO_Uncertainty__1up",
//        "PH_EFF_ID_Uncertainty__1up",
//        "PH_EFF_TRIGGER_Uncertainty__1up",
        "EG_RESOLUTION_ALL__1down",
        "EG_SCALE_ALL__1down",
        "PH_EFF_ISO_Uncertainty__1down",
        "PH_EFF_ID_Uncertainty__1down",
        "PH_EFF_TRIGGER_Uncertainty__1down",
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
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root"},
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

        auto resolved_reco_photons_matched = photon_passes_cuts.Define("chosen_two",
        [](const RVec<AbstractParticle>& reco_photons_matched)
        {
            RVec<AbstractParticle> x;
            if (reco_photons_matched.size() < 2)
            {
                return x;
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
                    x = {reco_photons_matched[combs[0][i]], reco_photons_matched[combs[1][i]]};
                }
            }
            if (pt1 > 10e3 && pt2 > 10e3 && best_X > 0.96 && best_X < 1.2 && chosen_delta_r < 1.5)
            {
                return x;
            }
            x.clear();
            return x;
        }, {"photons_pass_cuts"}).Filter(
        [&](RVec<AbstractParticle>& reco_photons_matched)
        {
            return (reco_photons_matched.size()==2);
        }, {"chosen_two"}).Define("reconstructed_mass",
        [&](RVec<AbstractParticle>& diph)
        {
            return (diph[0].PhotonVector()+diph[1].PhotonVector()).M()/1e3;
        }, {"chosen_two"});

        auto SB = resolved_reco_photons_matched.Filter(
        [](const double reconstructed_mass)
        {
            return (reconstructed_mass < 110) || (reconstructed_mass > 140);
        }, {"reconstructed_mass"});

        auto SR = resolved_reco_photons_matched.Filter(
        [](const double reconstructed_mass)
        {
            return (reconstructed_mass >= 110) && (reconstructed_mass <= 140);
        }, {"reconstructed_mass"});

        auto totEventWeight = resolved_reco_photons_matched
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

    ss << R"--(\section*{Table 24})--" << '\n';
    ss << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
    ss << R"--(\setlength\extrarowheight{2pt}\renewcommand{\arraystretch}{1.5})--" << '\n';
    ss << R"--(\begin{tabular}{|c|c|c|c|c|c|})--" << '\n';
    ss << R"--(\hline)--" << '\n';

        ss << R"--(\multicolumn{6}{|c|}{\parbox{\linewidth}{\centering Resolved Photon Category: Down Variations \\ (\% difference from nominal)}}\\[5 pt])--" << '\n';
    ss << R"--(\hline)--" << '\n';

    ss << R"--({Sample} & EG\_RESOLUTION\_ALL & EG\_SCALE\_ALL & PH\_EFF\_ISO\_Uncertainty & PH\_EFF\_ID\_Uncertainty & PH\_EFF\_TRIGGER\_Uncertainty \\ \hline)--" << '\n';

    double finalScaleVal;

    double ZgammaNominal = 0, ZgammaEG_RESOLUTION_ALL = 0,
    ZgammaEG_SCALE_ALL = 0, ZgammaPH_EFF_ISO_Uncertainty = 0,
    ZgammaPH_EFF_ID_Uncertainty = 0, ZgammaPH_EFF_TRIGGER_Uncertainty = 0;

    double ZjetsNominal = 0, ZjetsEG_RESOLUTION_ALL = 0,
    ZjetsEG_SCALE_ALL = 0, ZjetsPH_EFF_ISO_Uncertainty = 0,
    ZjetsPH_EFF_ID_Uncertainty = 0, ZjetsPH_EFF_TRIGGER_Uncertainty = 0;

    double totbkgNominal = 0, totbkgEG_RESOLUTION_ALL = 0,
    totbkgEG_SCALE_ALL = 0, totbkgPH_EFF_ISO_Uncertainty = 0,
    totbkgPH_EFF_ID_Uncertainty = 0, totbkgPH_EFF_TRIGGER_Uncertainty = 0;

    for (auto i = 0; (i < resultmaps.size()); i++)
    {
        double nominalVal = static_cast<double>(resultmaps[i]["nominal"]);

//        auto denominator = *Totals[i].GetResultPtr<ULong64_t>();
        auto denominator = *Totals[i].GetResultPtr<float>();

        if (i >= 0 && i <= 2) //Zgamma
        {
            finalScaleVal = SFs[i]/denominator;
            ZgammaNominal += finalScaleVal*nominalVal;

            ZgammaEG_RESOLUTION_ALL += finalScaleVal*resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1down"];
            ZgammaEG_SCALE_ALL += finalScaleVal*resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1down"];
            ZgammaPH_EFF_ISO_Uncertainty += finalScaleVal*resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1down"];
            ZgammaPH_EFF_ID_Uncertainty += finalScaleVal*resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1down"];
            ZgammaPH_EFF_TRIGGER_Uncertainty += finalScaleVal*resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down"];
        }

        else if (i >= 6)
        {
            finalScaleVal = JetNumeratorSFs[i-6]/denominator;
            ZjetsNominal += finalScaleVal*nominalVal;

            ZjetsEG_RESOLUTION_ALL += finalScaleVal*resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1down"];
            ZjetsEG_SCALE_ALL += finalScaleVal*resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1down"];
            ZjetsPH_EFF_ISO_Uncertainty += finalScaleVal*resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1down"];
            ZjetsPH_EFF_ID_Uncertainty += finalScaleVal*resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1down"];
            ZjetsPH_EFF_TRIGGER_Uncertainty += finalScaleVal*resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down"];
        }

        ss << prefixes[i] << " & ";
        ss << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1down"]) ? ((resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1down"]-nominalVal)/nominalVal)*100.0 : 0.0)
        << " & " << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1down"]) ? ((resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1down"]-nominalVal)/nominalVal)*100.0 : 0.0)
        << " & " << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1down"]) ? ((resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1down"]-nominalVal)/nominalVal)*100.0 : 0.0)
        << " & " << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1down"]) ? ((resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1down"]-nominalVal)/nominalVal)*100.0 : 0.0)
        << " & " << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down"]) ? ((resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down"]-nominalVal)/nominalVal)*100.0 : 0.0)
        << R"--( \\ \hline)--" << '\n';

        if (!resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1down"])
        {
            std::cout << i << ": photon_iso_eff:PH_EFF_ISO_Uncertainty__1down not found\n";
        }
        else
        {
            std::cout << i << ':' << ' ' << resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1down"] << ' ' << nominalVal << '\n';
        }
        if (!resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1down"])
        {
            std::cout << i << ": photon_id_eff:PH_EFF_ID_Uncertainty__1down not found\n";
        }
        else
        {
            std::cout << i << ':' << ' ' << resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1down"] << ' ' << nominalVal << '\n';
        }
        if (!resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down"])
        {
            std::cout << i << ": photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down not found\n";
        }
        else
        {
            std::cout << i << ':' << ' ' << resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down"] << ' ' << nominalVal << '\n';
        }
    }

    ss << R"--(Total $Z\gamma$ & )--";
    ss << std::setprecision(4) << std::fixed
    << ( (ZgammaNominal && ZgammaEG_RESOLUTION_ALL) ? ((ZgammaEG_RESOLUTION_ALL-ZgammaNominal)/ZgammaNominal)*100.0 : 0.0)
    << " & " << std::setprecision(4) << std::fixed
    << ( (ZgammaNominal && ZgammaEG_SCALE_ALL) ? ((ZgammaEG_SCALE_ALL-ZgammaNominal)/ZgammaNominal)*100.0 : 0.0)
    << " & " << std::setprecision(4) << std::fixed
    << ( (ZgammaNominal && ZgammaPH_EFF_ISO_Uncertainty) ? ((ZgammaPH_EFF_ISO_Uncertainty-ZgammaNominal)/ZgammaNominal)*100.0 : 0.0)
    << " & " << std::setprecision(4) << std::fixed
    << ( (ZgammaNominal && ZgammaPH_EFF_ID_Uncertainty) ? ((ZgammaPH_EFF_ID_Uncertainty-ZgammaNominal)/ZgammaNominal)*100.0 : 0.0)
    << " & " << std::setprecision(4) << std::fixed
    << ( (ZgammaNominal && ZgammaPH_EFF_TRIGGER_Uncertainty) ? ((ZgammaPH_EFF_TRIGGER_Uncertainty-ZgammaNominal)/ZgammaNominal)*100.0 : 0.0)
    << R"--( \\ \hline)--" << '\n';

    ss << R"--(Total $Z$ jets & )--";
    ss << std::setprecision(4) << std::fixed
    << ( (ZjetsNominal && ZjetsEG_RESOLUTION_ALL) ? ((ZjetsEG_RESOLUTION_ALL-ZjetsNominal)/ZjetsNominal)*100.0 : 0.0)
    << " & " << std::setprecision(4) << std::fixed
    << ( (ZjetsNominal && ZjetsEG_SCALE_ALL) ? ((ZjetsEG_SCALE_ALL-ZjetsNominal)/ZjetsNominal)*100.0 : 0.0)
    << " & " << std::setprecision(4) << std::fixed
    << ( (ZjetsNominal && ZjetsPH_EFF_ISO_Uncertainty) ? ((ZjetsPH_EFF_ISO_Uncertainty-ZjetsNominal)/ZjetsNominal)*100.0 : 0.0)
    << " & " << std::setprecision(4) << std::fixed
    << ( (ZjetsNominal && ZjetsPH_EFF_ID_Uncertainty) ? ((ZjetsPH_EFF_ID_Uncertainty-ZjetsNominal)/ZjetsNominal)*100.0 : 0.0)
    << " & " << std::setprecision(4) << std::fixed
    << ( (ZjetsNominal && ZjetsPH_EFF_TRIGGER_Uncertainty) ? ((ZjetsPH_EFF_TRIGGER_Uncertainty-ZjetsNominal)/ZjetsNominal)*100.0 : 0.0)
    << R"--( \\ \hline)--" << '\n';

    totbkgNominal = ZgammaNominal + ZjetsNominal;

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
    << ( (totbkgNominal && totbkgPH_EFF_ISO_Uncertainty) ? ((totbkgPH_EFF_ISO_Uncertainty-totbkgNominal)/totbkgNominal)*100.0 : 0.0)
    << " & " << std::setprecision(4) << std::fixed
    << ( (totbkgNominal && totbkgPH_EFF_ID_Uncertainty) ? ((totbkgPH_EFF_ID_Uncertainty-totbkgNominal)/totbkgNominal)*100.0 : 0.0)
    << " & " << std::setprecision(4) << std::fixed
    << ( (totbkgNominal && totbkgPH_EFF_TRIGGER_Uncertainty) ? ((totbkgPH_EFF_TRIGGER_Uncertainty-totbkgNominal)/totbkgNominal)*100.0 : 0.0)
    << R"--( \\ \hline)--" << '\n';

    ss << R"--(\end{tabular}})--" << '\n';

    ss << "\n\n\n";

    std::cout << ss.str();
}

void Table21_Displaced_Axions()
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

    std::vector<std::vector<std::string>> input_filenames = {
        //Signal
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root"},
        //Displaced Signal
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_FewMassPoints.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_FewMassPoints_PhotonSFs.root"},
    };

    constexpr std::array<const char*,15> prefixes = {R"--(Prompt Signal $m_{\text{A}}$ = 1 GeV)--", R"--(Prompt Signal $m_{\text{A}}$ = 5 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 0.2 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 0.5 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 1 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 3 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 5 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 10 GeV)--", /*R"--(Displaced Signal $m_{\text{A}}$ = 10.1 GeV)--",*/ R"--(Displaced Signal $m_{\text{A}}$ = 20 GeV)--",/* R"--(Displaced Signal $m_{\text{A}}$ = 20.1 GeV)--",*/ R"--(Displaced Signal $m_{\text{A}}$ = 29.5 GeV)--", };

    std::vector<float> massPoints = {0.2, 0.5, 1, 3, 5, 10, 20, 29.5};//{0.2,0.5,1,3,5,10,10.1,20,20.1,29.5};
    std::vector<ROOT::RDF::RResultHandle> Totals, Nodes, evtNodes;
    std::vector<RResultMap<float>> resultmaps;
    std::stringstream ss, cs;
    
    int counter = 0;

    for (auto& i: input_filenames)
    {
        SchottDataFrame df(MakeRDF(i,8));

        auto EventWeight = df.Define("EventWeight",
        [](const RVec<float>& ei_event_weights_generator)
        {
            return  ((ei_event_weights_generator[0]) ? 1 / ei_event_weights_generator[0] : 1);

        }, {"ei_event_weights_generator"}).Define("truth_axions", [&](RVec<TruthParticle> truth_particles)
        {
            truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
            [](TruthParticle& x)
            {
                return (abs(x.mc_pdg_id) != 36); //axions are pdg_id = 36 in displaced samples, but we don't use them in the prompt samples, so it doesn't matter that in the prompt samples, the pdg_id of the axions is 35

            }), truth_particles.end());

            return truth_particles;

        }, {"truth_particles"}).Define("axion_masses", [&](RVec<TruthParticle>& truth_axions)
        {
            return truth_axions[0].mc_mass/1e3f;

        }, {"truth_axions"});

        auto trigger_selection = df.Define("truth_axions", [&](RVec<TruthParticle> truth_particles)
        {
            truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
            [](TruthParticle& x)
            {
               return (abs(x.mc_pdg_id) != 36); //axions are pdg_id = 36 in displaced samples, but we don't use them in the prompt samples, so it doesn't matter that in the prompt samples, the pdg_id of the axions is 35

            }), truth_particles.end());

            return truth_particles;

        }, {"truth_particles"}).Define("axion_masses", [&](RVec<TruthParticle>& truth_axions)
        {
            return truth_axions[0].mc_mass/1e3f;

        }, {"truth_axions"}).Filter(
        [](const RVec<std::string>& trigger_passed_triggers)
        {
            bool trigger_found = (std::find_first_of(trigger_passed_triggers.begin(), trigger_passed_triggers.end(), triggers.begin(), triggers.end()) != trigger_passed_triggers.end());

            if (!trigger_found)
            {
                return false;
            }
            return true;

        }, {"trigger_passed_triggers"});

        auto two_leptons = trigger_selection.Define("di_electrons",
        [](RVec<AbstractParticle> electrons)
        {
            electrons.erase(std::remove_if(electrons.begin(),electrons.end(),
            [](AbstractParticle& ep)
            {
                //keep electron if it has pt > 20 GeV, |eta| < 2.37, not (1.37 < |eta| < 1.52), and id_medium = 1

                return (!((ep.electron_pt/1e3 > 20) && (abs(ep.electron_eta) < 2.37) &&
                (!((1.37 < abs(ep.electron_eta)) && (abs(ep.electron_eta) < 1.52)))
                && (ep.electron_id_medium == 1)));

            }), electrons.end());

            return electrons;

        },{"abstract_electrons"}).Filter([](const RVec<Muon>& muons, RVec<AbstractParticle> electrons)
        {
            return (electrons.size()==2 && muons.empty());

        }, {"muons", "di_electrons"});

        auto opp_charge = two_leptons.Filter([](const RVec<AbstractParticle>& electrons)
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

        auto photon_passes_cuts = ptCut.Define("abstract_photons_pass_cut_indices",[&](RVec<AbstractParticle>& p) //p = photon
        {
            RVec<int> x; //indices of photons that pass cuts

            for (auto i = 0; i < p.size(); i++)
            {
                if (not ((abs(p[i].photon_eta) >= 2.37) or (abs(p[i].photon_eta) > 1.37 and abs(p[i].photon_eta) < 1.52) or (not p[i].photon_id_loose)))
                {
                    x.push_back(i);
                }
            }

            return x;
        }, {"abstract_photons"}).Define("photons_pass_cuts",
        [&](RVec<AbstractParticle>& photons, RVec<int>& x)
        {
            return Take(photons, x); //Taking only the photons that passed the cuts in each event

        }, {"abstract_photons", "abstract_photons_pass_cut_indices"}).Define("Event_number_and_photon",
        [&](const RVec<AbstractParticle>& reco_photons_matched, UInt_t ei_event_number)
        {
            return std::make_pair(ei_event_number, reco_photons_matched);
            
        }, {"photons_pass_cuts", "ei_event_number"});

        auto merged_reco_photons_matched = photon_passes_cuts.Filter(
        [&](const RVec<AbstractParticle>& reco_photons_matched)
        {
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

            for (size_t i=0; i<length; i++) //looping through all of the possible combinations of photons in each event
            {
                delta_r = DeltaR(reco_photons_matched[combs[0][i]].PhotonVector(), reco_photons_matched[combs[1][i]].PhotonVector());
                m = (reco_photons_matched[combs[0][i]].PhotonVector() + reco_photons_matched[combs[1][i]].PhotonVector()).M();
                pt = (reco_photons_matched[combs[0][i]].PhotonVector() + reco_photons_matched[combs[1][i]].PhotonVector()).Pt();
                X = delta_r*(pt/(2.0*m));
                if (i==0 || ((abs(1-X) < abs(1-best_X)) and (delta_r < 1.5)))
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
        .Define("merged_photon_index",
        [&](const RVec<AbstractParticle>& rpm) //rpm = reco photons matched
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
        .Define("merged_photon",
        [&](const RVec<AbstractParticle>& reco_photons_matched, int merged_photon_index)
        {
            return reco_photons_matched[merged_photon_index];

        }, {"photons_pass_cuts", "merged_photon_index"});

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

        //mpi = merged_photon_index
        auto totEventWeight = merged_reco_photons_matched.Define("totEventWeight", [](RVec<int>& x, int mpi, RVec<float> photon_id_eff, RVec<float> photon_iso_eff, RVec<float> photon_trg_eff/*, RVec<float> ei_event_weights_generator*/)
        {
            //   ||  jic they don't already have the same size  ||
            //   \/                                             \/
            auto ResizeVal = std::max({photon_id_eff.size(), photon_iso_eff.size(), photon_trg_eff.size()});
            photon_id_eff.resize(ResizeVal,1);
            photon_iso_eff.resize(ResizeVal,1);
            photon_trg_eff.resize(ResizeVal,1);

            //First, take the x indices from the respective photon_*_eff vectors,
            //this corresponds to the photons from the set of all reco photons in
            //the event that passed the "photon_passes_cuts" cuts. Then, from
            //those photons, select the index mpi element that corresponds to
            //the merged reco-photon

            return Take(photon_id_eff, x)[mpi] * Take(photon_iso_eff, x)[mpi] * Take(photon_trg_eff, x)[mpi];

        }, {"abstract_photons_pass_cut_indices", "merged_photon_index", "photon_id_eff", "photon_iso_eff", "photon_trg_eff",/* "ei_event_weights_generator"*/});

        if (counter < 2)
        {
            Totals.push_back(EventWeight.Sum<float>("EventWeight"));
            resultmaps.push_back(VariationsFor(totEventWeight.Sum<float>("totEventWeight")));
        }
        else
        {
            evtNodes.push_back(photon_passes_cuts.Take< std::pair<UInt_t, RVec<AbstractParticle>> , RVec<std::pair<UInt_t, RVec<AbstractParticle>>> >("Event_number_and_photon"));
            
            for (auto& mass_point: massPoints)
            {
                auto mass_point_EventWeight = EventWeight.Filter([&]
                (float axion_mass)
                {
                    return (roundToOneDecimalPlace(axion_mass) == mass_point);
//                    return ((axion_mass > mass_point - 0.1f) and (axion_mass < mass_point + 0.1f));
                }, {"axion_masses"});

                auto mass_point_totEventWeight = totEventWeight.Filter([&]
                (float axion_mass)
                {
                    return (roundToOneDecimalPlace(axion_mass) == mass_point);
//                    return ((axion_mass > mass_point - 0.1f) and (axion_mass < mass_point + 0.1f));
                }, {"axion_masses"});

                Totals.push_back(mass_point_EventWeight.Sum<float>("EventWeight"));
                resultmaps.push_back(VariationsFor(mass_point_totEventWeight.Sum<float>("totEventWeight")));
                Nodes.push_back(mass_point_totEventWeight.Take<UInt_t, RVec<UInt_t>>("ei_event_number"));
                
            }
        }

        counter++;
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
    
    std::cout << evtNodes.size() << '\n';
    
    for (int i = 0; i < evtNodes.size(); i++)
    {
        for (auto&j: *evtNodes[i].GetResultPtr<RVec<std::pair<UInt_t, RVec<AbstractParticle>>>>())
        {
            if (j.first == 370 and j.second.size() >= 2)
            {
                std::cout << j.first << " {";
//                for (auto& k: j.second)
//                {
//                    std::cout << "{photon pt = " << k.photon_pt <<
//                    ", photon eta = " << k.photon_eta <<
//                    ", photon phi = " << k.photon_phi <<
//                    "}, ";
//                }
//                std::cout << "} \n";
                
                auto combs = Combinations(j.second, 2); //combinations of photons
                
                size_t length = combs[0].size();

                double delta_r, m, pt, X, best_X, pt1, pt2, chosen_delta_r;

                for (size_t i=0; i<length; i++) //looping through all of the possible combinations of photons in each event
                {
                    std::cout << "{photon pt_1 = " << j.second[combs[0][i]].photon_pt
                    << ", photon eta_1 = " << j.second[combs[0][i]].photon_eta
                    << ", photon phi_1 = " << j.second[combs[0][i]].photon_phi
                    << "; photon pt_2 = " << j.second[combs[1][i]].photon_pt
                    << ", photon eta_2 = " << j.second[combs[1][i]].photon_eta
                    << ", photon phi_2 = " << j.second[combs[1][i]].photon_phi;
                    
                    delta_r = DeltaR(j.second[combs[0][i]].PhotonVector(), j.second[combs[1][i]].PhotonVector());
                    m = (j.second[combs[0][i]].PhotonVector() + j.second[combs[1][i]].PhotonVector()).M();
                    pt = (j.second[combs[0][i]].PhotonVector() + j.second[combs[1][i]].PhotonVector()).Pt();
                    X = delta_r*(pt/(2.0*m));
                    
                    std::cout << "; deltaR = " << delta_r <<
                    "; X = " << X << "}, ";
                    
                    if (i==0 || ((abs(1-X) < abs(1-best_X)) and (delta_r < 1.5)))
                    {
                        best_X = X;
                        pt1 = j.second[combs[0][i]].photon_pt;
                        pt2 = j.second[combs[1][i]].photon_pt;
                        chosen_delta_r = delta_r;
                    }
                }
                
                std::cout << " Best deltaR = " << chosen_delta_r
                << ", Best X = " << best_X << "}\n";
         
            }
        }
    }

    std::cout << "Merged Category Up\n==================\n";
    for (int i = 0; i < Nodes.size(); i++)
    {
        std::cout << prefixes[i+2] << '\n';
        std::cout << "==========\n";
        for (auto& j: *Nodes[i].GetResultPtr<RVec<UInt_t>>())
        {
            std::cout << j << ", ";
        }
        std::cout << "\n\n";
    }

//          resultmaps
//          ----------
//    0    1    ma1
//    2    3    ma5
//    4    5    displaced_axion_1
//    6    7    displaced_axion_2
//    8    9    displaced_axion_3
//    10   11   displaced_axion_4
//    12   13   displaced_axion_5
//    14   15   displaced_axion_6
//    16   17   displaced_axion_7
//    18   19   displaced_axion_8
//    20   21   displaced_axion_9
//    22   23   displaced_axion_10

    ss << R"--(\section*{Table 21 Prompt and Displaced Signal Samples})--" << '\n';
    ss << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
    ss << R"--(\setlength\extrarowheight{2pt}\renewcommand{\arraystretch}{1.5})--" << '\n';
    ss << R"--(\begin{tabular}{|c|c|c|c|c|c|})--" << '\n';
    ss << R"--(\hline)--" << '\n';

        ss << R"--(\multicolumn{6}{|c|}{\parbox{\linewidth}{\centering Merged Photon Category: Up Variations \\ (\% difference from nominal)}}\\[5 pt])--" << '\n';
    ss << R"--(\hline)--" << '\n';

    ss << R"--({Sample} & EG\_RESOLUTION\_ALL & EG\_SCALE\_ALL & PH\_EFF\_ISO\_Uncertainty & PH\_EFF\_ID\_Uncertainty & PH\_EFF\_TRIGGER\_Uncertainty \\ \hline)--" << '\n';
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    cs << R"--(\section*{Table 21 Prompt and Displaced Signal Samples})--" << '\n';
    cs << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
    cs << R"--(\setlength\extrarowheight{2pt}\renewcommand{\arraystretch}{1.5})--" << '\n';
    cs << R"--(\begin{tabular}{|c|c|c|c|c|c|})--" << '\n';
    cs << R"--(\hline)--" << '\n';

    cs << R"--(\multicolumn{6}{|c|}{\parbox{\linewidth}{\centering Merged Photon Category: Up Variations \\ (raw number of events)}}\\[5 pt])--" << '\n';
    cs << R"--(\hline)--" << '\n';

    cs << R"--({Sample} & EG\_RESOLUTION\_ALL & EG\_SCALE\_ALL & PH\_EFF\_ISO\_Uncertainty & PH\_EFF\_ID\_Uncertainty & PH\_EFF\_TRIGGER\_Uncertainty \\ \hline)--" << '\n';

    for (auto i = 0; (i < resultmaps.size()); i++)
    {
        float nominalVal = static_cast<float>(resultmaps[i]["nominal"]);

        ss << prefixes[i] << " & ";
        ss << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1up"]) ? ((resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1up"]-nominalVal)/nominalVal)*100.0 : 0.0)
        << " & " << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1up"]) ? ((resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1up"]-nominalVal)/nominalVal)*100.0 : 0.0)
        << " & " << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1up"]) ? ((resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1up"]-nominalVal)/nominalVal)*100.0 : 0.0)
        << " & " << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1up"]) ? ((resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1up"]-nominalVal)/nominalVal)*100.0 : 0.0)
        << " & " << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up"]) ? ((resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up"]-nominalVal)/nominalVal)*100.0 : 0.0)
        << R"--( \\ \hline)--" << '\n';
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        cs << prefixes[i] << " & ";
        cs << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1up"]) ? resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1up"] : 0.0)
        << " & " << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1up"]) ? resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1up"] : 0.0)
        << " & " << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1up"]) ? resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1up"] : 0.0)
        << " & " << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1up"]) ? resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1up"] : 0.0)
        << " & " << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up"]) ? resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up"] : 0.0)
        << R"--( \\ \hline)--" << '\n';

        if (!resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1up"])
        {
            std::cout << i << ": photon_iso_eff:PH_EFF_ISO_Uncertainty__1up not found\n";
        }
        else
        {
            std::cout << i << ':' << ' ' << resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1up"] << ' ' << nominalVal << '\n';
        }
        if (!resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1up"])
        {
            std::cout << i << ": photon_id_eff:PH_EFF_ID_Uncertainty__1up not found\n";
        }
        else
        {
            std::cout << i << ':' << ' ' << resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1up"] << ' ' << nominalVal << '\n';
        }
        if (!resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up"])
        {
            std::cout << i << ": photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up not found\n";
        }
        else
        {
            std::cout << i << ':' << ' ' << resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up"] << ' ' << nominalVal << '\n';
        }
    }

    ss << R"--(\end{tabular}})--" << '\n';
    cs << R"--(\end{tabular}})--" << '\n';

    ss << "\n\n\n";
    cs << "\n\n\n";

    std::cout << "\n\n\n" << ss.str();
    std::cout << "\n\n\n" << cs.str();

}

void Table22_Displaced_Axions()
{
    Event::systematics =
    {

        "EG_RESOLUTION_ALL__1down",
        "EG_SCALE_ALL__1down",
        "PH_EFF_ISO_Uncertainty__1down",
        "PH_EFF_ID_Uncertainty__1down",
        "PH_EFF_TRIGGER_Uncertainty__1down",
//        "EG_RESOLUTION_ALL__1down",
    };

    std::vector<std::vector<std::string>> input_filenames = {
        //Signal
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root"},
        //Displaced Signal
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_FewMassPoints.root"},
    };

    constexpr std::array<const char*,15> prefixes = {R"--(Prompt Signal $m_{\text{A}}$ = 1 GeV)--", R"--(Prompt Signal $m_{\text{A}}$ = 5 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 0.2 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 0.5 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 1 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 3 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 5 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 10 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 10.1 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 20 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 20.1 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 29.5 GeV)--", };

    std::vector<float> massPoints = {0.2,0.5,1,3,5,10,10.1,20,20.1,29.5};
    std::vector<ROOT::RDF::RResultHandle> Totals, Nodes;
    std::vector<RResultMap<float>> resultmaps;
    std::stringstream ss, cs;

    int counter = 0;

    for (auto& i: input_filenames)
    {
        SchottDataFrame df(MakeRDF(i,8));

        auto EventWeight = df.Define("EventWeight",
        [](const RVec<float>& ei_event_weights_generator)
        {
            return  ((ei_event_weights_generator[0]) ? 1 / ei_event_weights_generator[0] : 1);

        }, {"ei_event_weights_generator"}).Define("truth_axions", [&](RVec<TruthParticle> truth_particles)
        {
            truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
            [](TruthParticle& x)
            {
             return (abs(x.mc_pdg_id) != 36);

            }), truth_particles.end());

            return truth_particles;

        }, {"truth_particles"}).Define("axion_masses", [&](RVec<TruthParticle>& truth_axions)
        {
            return truth_axions[0].mc_mass/1e3f;

        }, {"truth_axions"});

        auto trigger_selection = df.Define("truth_axions", [&](RVec<TruthParticle> truth_particles)
        {
            truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
            [](TruthParticle& x)
            {
               return (abs(x.mc_pdg_id) != 36);

            }), truth_particles.end());

            return truth_particles;

        }, {"truth_particles"}).Define("axion_masses", [&](RVec<TruthParticle>& truth_axions)
        {
            return truth_axions[0].mc_mass/1e3f;

        }, {"truth_axions"}).Filter(
        [](const RVec<std::string>& trigger_passed_triggers)
        {
            bool trigger_found = (std::find_first_of(trigger_passed_triggers.begin(), trigger_passed_triggers.end(), triggers.begin(), triggers.end()) != trigger_passed_triggers.end());

            if (!trigger_found)
            {
                return false;
            }
            return true;

        }, {"trigger_passed_triggers"});

        auto two_leptons = trigger_selection.Define("di_electrons",
        [](RVec<AbstractParticle> electrons)
        {
            electrons.erase(std::remove_if(electrons.begin(),electrons.end(),
            [](AbstractParticle& ep)
            {
                //keep electron if it has pt > 20 GeV, |eta| < 2.37, not (1.37 < |eta| < 1.52), and id_medium = 1

                return (!((ep.electron_pt/1e3 > 20) && (abs(ep.electron_eta) < 2.37) &&
                (!((1.37 < abs(ep.electron_eta)) && (abs(ep.electron_eta) < 1.52)))
                && (ep.electron_id_medium == 1)));

            }), electrons.end());

            return electrons;

        },{"abstract_electrons"}).Filter([](const RVec<Muon>& muons, RVec<AbstractParticle> electrons)
        {
            return (electrons.size()==2 && muons.empty());

        }, {"muons", "di_electrons"});

        auto opp_charge = two_leptons.Filter([](const RVec<AbstractParticle>& electrons)
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

        auto photon_passes_cuts = ptCut.Define("abstract_photons_pass_cut_indices",[&](RVec<AbstractParticle>& p) //p = photon
        {
            RVec<int> x; //indices of photons that pass cuts

            for (auto i = 0; i < p.size(); i++)
            {
                if (not ((abs(p[i].photon_eta) >= 2.37) or (abs(p[i].photon_eta) > 1.37 and abs(p[i].photon_eta) < 1.52) or (not p[i].photon_id_loose)))
                {
                    x.push_back(i);
                }
            }

            return x;
        }, {"abstract_photons"}).Define("photons_pass_cuts",
        [&](RVec<AbstractParticle>& photons, RVec<int>& x)
        {
            return Take(photons, x); //Taking only the photons that passed the cuts in each event

        }, {"abstract_photons", "abstract_photons_pass_cut_indices"});

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

            for (size_t i=0; i<length; i++) //looping through all of the possible combinations of photons in each event
            {
                delta_r = DeltaR(reco_photons_matched[combs[0][i]].PhotonVector(), reco_photons_matched[combs[1][i]].PhotonVector());
                m = (reco_photons_matched[combs[0][i]].PhotonVector() + reco_photons_matched[combs[1][i]].PhotonVector()).M();
                pt = (reco_photons_matched[combs[0][i]].PhotonVector() + reco_photons_matched[combs[1][i]].PhotonVector()).Pt();
                X = delta_r*(pt/(2.0*m));
                if (i==0 || ((abs(1-X) < abs(1-best_X)) and (delta_r < 1.5)))
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
        .Define("merged_photon_index",
        [&](const RVec<AbstractParticle>& rpm) //rpm = reco photons matched
        {
            for (auto i = 0; i < rpm.size(); i++)
            {
                if (rpm[i].photon_pt > 20e3)
                {
                    return i; //returning the index of the first photon that has photon_pt > 20 GeV
                }
            }
            return 0; //jic the compiler complains

        }, {"photons_pass_cuts"})
        .Define("merged_photon",
        [&](const RVec<AbstractParticle>& reco_photons_matched, int merged_photon_index)
        {
            return reco_photons_matched[merged_photon_index];

        }, {"photons_pass_cuts", "merged_photon_index"});

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

        //mpi = merged_photon_index
        auto totEventWeight = merged_reco_photons_matched.Define("totEventWeight", [](RVec<int>& x, int mpi, RVec<float> photon_id_eff, RVec<float> photon_iso_eff, RVec<float> photon_trg_eff/*, RVec<float> ei_event_weights_generator*/)
        {
            //   ||  jic they don't already have the same size  ||
            //   \/                                             \/
            auto ResizeVal = std::max({photon_id_eff.size(), photon_iso_eff.size(), photon_trg_eff.size()});
            photon_id_eff.resize(ResizeVal,1);
            photon_iso_eff.resize(ResizeVal,1);
            photon_trg_eff.resize(ResizeVal,1);

            //First, take the x indices from the respective photon_*_eff vectors,
            //this corresponds to the photons from the set of all reco photons in
            //the event that passed the "photon_passes_cuts" cuts. Then, from
            //those photons, select the index mpi element that corresponds to
            //the merged reco-photon

            return Take(photon_id_eff, x)[mpi] * Take(photon_iso_eff, x)[mpi] * Take(photon_trg_eff, x)[mpi];

        }, {"abstract_photons_pass_cut_indices", "merged_photon_index", "photon_id_eff", "photon_iso_eff", "photon_trg_eff",/* "ei_event_weights_generator"*/});

        if (counter < 2)
        {
            Totals.push_back(EventWeight.Sum<float>("EventWeight"));
            resultmaps.push_back(VariationsFor(totEventWeight.Sum<float>("totEventWeight")));
        }
        else
        {
            for (auto& mass_point: massPoints)
            {
                auto mass_point_EventWeight = EventWeight.Filter([&]
                (float axion_mass)
                {
                    return (roundToOneDecimalPlace(axion_mass) == mass_point);
//                    return ((axion_mass > mass_point - 0.1f) and (axion_mass < mass_point + 0.1f));
                }, {"axion_masses"});

                auto mass_point_totEventWeight = totEventWeight.Filter([&]
                (float axion_mass)
                {
                    return (roundToOneDecimalPlace(axion_mass) == mass_point);
//                    return ((axion_mass > mass_point - 0.1f) and (axion_mass < mass_point + 0.1f));
                }, {"axion_masses"});

                Totals.push_back(mass_point_EventWeight.Sum<float>("EventWeight"));
                resultmaps.push_back(VariationsFor(mass_point_totEventWeight.Sum<float>("totEventWeight")));
                Nodes.push_back(mass_point_totEventWeight.Take<UInt_t, RVec<UInt_t>>("ei_event_number"));
            }
        }

        counter++;
    }

    ROOT::RDF::RunGraphs(Totals); // running all computation nodes concurrently

    std::cout << "Merged Category Down\n====================\n";
    for (int i = 0; i < Nodes.size(); i++)
    {
        std::cout << prefixes[i+2] << '\n';
        std::cout << "==========\n";
        for (auto& j: *Nodes[i].GetResultPtr<RVec<UInt_t>>())
        {
            std::cout << j << ", ";
        }
        std::cout << "\n\n";
    }

//          resultmaps
//          ----------
//    0    1    ma1
//    2    3    ma5
//    4    5    displaced_axion_1
//    6    7    displaced_axion_2
//    8    9    displaced_axion_3
//    10   11   displaced_axion_4
//    12   13   displaced_axion_5
//    14   15   displaced_axion_6
//    16   17   displaced_axion_7
//    18   19   displaced_axion_8
//    20   21   displaced_axion_9
//    22   23   displaced_axion_10

    ss << R"--(\section*{Table 22 Prompt and Displaced Signal Samples})--" << '\n';
    ss << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
    ss << R"--(\setlength\extrarowheight{2pt}\renewcommand{\arraystretch}{1.5})--" << '\n';
    ss << R"--(\begin{tabular}{|c|c|c|c|c|c|})--" << '\n';
    ss << R"--(\hline)--" << '\n';

        ss << R"--(\multicolumn{6}{|c|}{\parbox{\linewidth}{\centering Merged Photon Category: Down Variations \\ (\% difference from nominal)}}\\[5 pt])--" << '\n';
    ss << R"--(\hline)--" << '\n';

    ss << R"--({Sample} & EG\_RESOLUTION\_ALL & EG\_SCALE\_ALL & PH\_EFF\_ISO\_Uncertainty & PH\_EFF\_ID\_Uncertainty & PH\_EFF\_TRIGGER\_Uncertainty \\ \hline)--" << '\n';
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    cs << R"--(\section*{Table 22 Prompt and Displaced Signal Samples})--" << '\n';
    cs << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
    cs << R"--(\setlength\extrarowheight{2pt}\renewcommand{\arraystretch}{1.5})--" << '\n';
    cs << R"--(\begin{tabular}{|c|c|c|c|c|c|})--" << '\n';
    cs << R"--(\hline)--" << '\n';

    cs << R"--(\multicolumn{6}{|c|}{\parbox{\linewidth}{\centering Merged Photon Category: Down Variations \\ (raw number of events)}}\\[5 pt])--" << '\n';
    cs << R"--(\hline)--" << '\n';

    cs << R"--({Sample} & EG\_RESOLUTION\_ALL & EG\_SCALE\_ALL & PH\_EFF\_ISO\_Uncertainty & PH\_EFF\_ID\_Uncertainty & PH\_EFF\_TRIGGER\_Uncertainty \\ \hline)--" << '\n';

    for (auto i = 0; (i < resultmaps.size()); i++)
    {
        float nominalVal = static_cast<float>(resultmaps[i]["nominal"]);

        ss << prefixes[i] << " & ";
        ss << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1down"]) ? ((resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1down"]-nominalVal)/nominalVal)*100.0 : 0.0)
        << " & " << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1down"]) ? ((resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1down"]-nominalVal)/nominalVal)*100.0 : 0.0)
        << " & " << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1down"]) ? ((resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1down"]-nominalVal)/nominalVal)*100.0 : 0.0)
        << " & " << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1down"]) ? ((resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1down"]-nominalVal)/nominalVal)*100.0 : 0.0)
        << " & " << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down"]) ? ((resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down"]-nominalVal)/nominalVal)*100.0 : 0.0)
        << R"--( \\ \hline)--" << '\n';
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        cs << prefixes[i] << " & ";
        cs << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1down"]) ? resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1down"] : 0.0)
        << " & " << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1down"]) ? resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1down"] : 0.0)
        << " & " << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1down"]) ? resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1down"] : 0.0)
        << " & " << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1down"]) ? resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1down"] : 0.0)
        << " & " << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down"]) ? resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down"] : 0.0)
        << R"--( \\ \hline)--" << '\n';

        if (!resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1down"])
        {
            std::cout << i << ": photon_iso_eff:PH_EFF_ISO_Uncertainty__1down not found\n";
        }
        else
        {
            std::cout << i << ':' << ' ' << resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1down"] << ' ' << nominalVal << '\n';
        }
        if (!resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1down"])
        {
            std::cout << i << ": photon_id_eff:PH_EFF_ID_Uncertainty__1down not found\n";
        }
        else
        {
            std::cout << i << ':' << ' ' << resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1down"] << ' ' << nominalVal << '\n';
        }
        if (!resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down"])
        {
            std::cout << i << ": photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down not found\n";
        }
        else
        {
            std::cout << i << ':' << ' ' << resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down"] << ' ' << nominalVal << '\n';
        }
    }

    ss << R"--(\end{tabular}})--" << '\n';
    cs << R"--(\end{tabular}})--" << '\n';

    ss << "\n\n\n";
    cs << "\n\n\n";

    std::cout << "\n\n\n" << ss.str();
    std::cout << "\n\n\n" << cs.str();

}

void Table23_Displaced_Axions()
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

    std::vector<std::vector<std::string>> input_filenames = {
        //Signal
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root"},
        //Displaced Signal
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_FewMassPoints.root"},
    };

    constexpr std::array<const char*,15> prefixes = {R"--(Prompt Signal $m_{\text{A}}$ = 1 GeV)--", R"--(Prompt Signal $m_{\text{A}}$ = 5 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 0.2 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 0.5 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 1 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 3 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 5 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 10 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 10.1 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 20 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 20.1 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 29.5 GeV)--", };

    std::vector<float> massPoints = {0.2,0.5,1,3,5,10,10.1,20,20.1,29.5};
    std::vector<ROOT::RDF::RResultHandle> Totals, Nodes;
    std::vector<RResultMap<float>> resultmaps;
    std::stringstream ss, cs;

    int counter = 0;

    for (auto& i: input_filenames)
    {
        SchottDataFrame df(MakeRDF(i,8));

        auto EventWeight = df.Define("EventWeight",
        [](const RVec<float>& ei_event_weights_generator)
        {
            return  ((ei_event_weights_generator[0]) ? 1 / ei_event_weights_generator[0] : 1);

        }, {"ei_event_weights_generator"}).Define("truth_axions", [&](RVec<TruthParticle> truth_particles)
        {
            truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
            [](TruthParticle& x)
            {
             return (abs(x.mc_pdg_id) != 36);

            }), truth_particles.end());

            return truth_particles;

        }, {"truth_particles"}).Define("axion_masses", [&](RVec<TruthParticle>& truth_axions)
        {
            return truth_axions[0].mc_mass/1e3f;

        }, {"truth_axions"});

        auto trigger_selection = df.Define("truth_axions", [&](RVec<TruthParticle> truth_particles)
        {
            truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
            [](TruthParticle& x)
            {
               return (abs(x.mc_pdg_id) != 36);

            }), truth_particles.end());

            return truth_particles;

        }, {"truth_particles"}).Define("axion_masses", [&](RVec<TruthParticle>& truth_axions)
        {
            return truth_axions[0].mc_mass/1e3f;

        }, {"truth_axions"}).Filter(
        [](const RVec<std::string>& trigger_passed_triggers)
        {
            bool trigger_found = (std::find_first_of(trigger_passed_triggers.begin(), trigger_passed_triggers.end(), triggers.begin(), triggers.end()) != trigger_passed_triggers.end());

            if (!trigger_found)
            {
                return false;
            }
            return true;

        }, {"trigger_passed_triggers"});

        auto two_leptons = trigger_selection.Define("di_electrons",
        [](RVec<AbstractParticle> electrons)
        {
            electrons.erase(std::remove_if(electrons.begin(),electrons.end(),
            [](AbstractParticle& ep)
            {
                //keep electron if it has pt > 20 GeV, |eta| < 2.37, not (1.37 < |eta| < 1.52), and id_medium = 1

                return (!((ep.electron_pt/1e3 > 20) && (abs(ep.electron_eta) < 2.37) &&
                (!((1.37 < abs(ep.electron_eta)) && (abs(ep.electron_eta) < 1.52)))
                && (ep.electron_id_medium == 1)));

            }), electrons.end());

            return electrons;

        },{"abstract_electrons"}).Filter([](const RVec<Muon>& muons, RVec<AbstractParticle> electrons)
        {
            return (electrons.size()==2 && muons.empty());

        }, {"muons", "di_electrons"});

        auto opp_charge = two_leptons.Filter([](const RVec<AbstractParticle>& electrons)
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

        auto photon_passes_cuts = ptCut.Define("abstract_photons_pass_cut_indices",[&](RVec<AbstractParticle>& p) //p = photon
        {
            RVec<int> x; //indices of photons that pass cuts

            for (auto i = 0; i < p.size(); i++)
            {
                if (not ((abs(p[i].photon_eta) >= 2.37) or (abs(p[i].photon_eta) > 1.37 and abs(p[i].photon_eta) < 1.52) or (not p[i].photon_id_loose)))
                {
                    x.push_back(i);
                }
            }

            return x;
        }, {"abstract_photons"}).Define("photons_pass_cuts",
        [&](RVec<AbstractParticle>& photons, RVec<int>& x)
        {
            return Take(photons, x); //Taking only the photons that passed the cuts in each event

        }, {"abstract_photons", "abstract_photons_pass_cut_indices"});

        auto resolved_reco_photons_matched = photon_passes_cuts.Define("reco_photons_matched_indices",
        [](const RVec<AbstractParticle>& reco_photons_matched)
        {
            RVec<unsigned long> x; //vector of indices
            if (reco_photons_matched.size() < 2)
            {
                return x;
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
                if (i==0 || ((abs(1-X) < abs(1-best_X)) and (delta_r < 1.5)))
                {
                    best_X = X;
                    pt1 = reco_photons_matched[combs[0][i]].photon_pt;
                    pt2 = reco_photons_matched[combs[1][i]].photon_pt;
                    chosen_delta_r = delta_r;
                    x = {combs[0][i], combs[1][i]};
                }
            }
            if (pt1 > 10e3 && pt2 > 10e3 && best_X > 0.96 && best_X < 1.2 && chosen_delta_r < 1.5)
            {
                return x;
            }

            x.clear();
            return x;

        }, {"photons_pass_cuts"}).Filter(
        [&](RVec<unsigned long>& reco_photons_matched_indices)
        {
            return (reco_photons_matched_indices.size()==2);

        }, {"reco_photons_matched_indices"}).Define("chosen_two",
        [](const RVec<AbstractParticle>& reco_photons_matched, RVec<unsigned long>& reco_photons_matched_indices)
        {
            return Take(reco_photons_matched, reco_photons_matched_indices);

        }, {"photons_pass_cuts", "reco_photons_matched_indices"}).Define("reconstructed_mass",
        [&](RVec<AbstractParticle>& diph)
        {
            return (diph[0].PhotonVector()+diph[1].PhotonVector()).M()/1e3;

        }, {"chosen_two"});

        auto SB = resolved_reco_photons_matched.Filter(
        [](const double reconstructed_mass)
        {
            return (reconstructed_mass < 110) || (reconstructed_mass > 140);
        }, {"reconstructed_mass"});

        auto SR = resolved_reco_photons_matched.Filter(
        [](const double reconstructed_mass)
        {
            return (reconstructed_mass >= 110) && (reconstructed_mass <= 140);
        }, {"reconstructed_mass"});

        auto totEventWeight = resolved_reco_photons_matched //rpmi = reco_photons_matched_indices
        .Define("totEventWeight", [](RVec<int>& x, RVec<unsigned long>& rpmi, RVec<float> photon_id_eff, RVec<float> photon_iso_eff, RVec<float> photon_trg_eff/*, RVec<float> ei_event_weights_generator*/)
        {
            //   ||  jic they don't already have the same size  ||
            //   \/                                             \/
            auto ResizeVal = std::max({photon_id_eff.size(), photon_iso_eff.size(), photon_trg_eff.size()});
            photon_id_eff.resize(ResizeVal,1);
            photon_iso_eff.resize(ResizeVal,1);
            photon_trg_eff.resize(ResizeVal,1);

            //First, take the x indices from the respective photon_*_eff vectors,
            //this corresponds to the photons from the set of all reco photons in
            //the events that passed the "photon_passes_cuts" cuts. Then, from
            //those photons, select the indices rpmi elements that corresponds to
            //the two resolved reco-photons

            RVec<float> two_efficiencies = Take(Take(photon_id_eff,x),rpmi) * Take(Take(photon_iso_eff,x),rpmi) * Take(Take(photon_trg_eff,x),rpmi);//*ei_event_weights_generator[0];

            float total = 1.0f;

            for (auto i: two_efficiencies)
            {
                total *= i;
            }

            return total;

        }, {"abstract_photons_pass_cut_indices", "reco_photons_matched_indices", "photon_id_eff", "photon_iso_eff", "photon_trg_eff",/* "ei_event_weights_generator"*/});

        if (counter < 2)
        {
            Totals.push_back(EventWeight.Sum<float>("EventWeight"));
            resultmaps.push_back(VariationsFor(totEventWeight.Sum<float>("totEventWeight")));
        }
        else
        {
            for (auto& mass_point: massPoints)
            {
                auto mass_point_EventWeight = EventWeight.Filter([&]
                (float axion_mass)
                {
                    return (roundToOneDecimalPlace(axion_mass) == mass_point);
//                    return ((axion_mass > mass_point - 0.1f) and (axion_mass < mass_point + 0.1f));
                }, {"axion_masses"});

                auto mass_point_totEventWeight = totEventWeight.Filter([&]
                (float axion_mass)
                {
                    return (roundToOneDecimalPlace(axion_mass) == mass_point);
//                    return ((axion_mass > mass_point - 0.1f) and (axion_mass < mass_point + 0.1f));
                }, {"axion_masses"});

                Totals.push_back(mass_point_EventWeight.Sum<float>("EventWeight"));
                resultmaps.push_back(VariationsFor(mass_point_totEventWeight.Sum<float>("totEventWeight")));
                Nodes.push_back(mass_point_totEventWeight.Take<UInt_t, RVec<UInt_t>>("ei_event_number"));
            }
        }

        counter++;

    }

    ROOT::RDF::RunGraphs(Totals); // running all computation nodes concurrently

    std::cout << "Resolved Category Up\n====================\n";
    for (int i = 0; i < Nodes.size(); i++)
    {
        std::cout << prefixes[i+2] << '\n';
        std::cout << "==========\n";
        for (auto& j: *Nodes[i].GetResultPtr<RVec<UInt_t>>())
        {
            std::cout << j << ", ";
        }
        std::cout << "\n\n";
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

    ss << R"--(\section*{Table 23 Prompt and Displaced Signal Samples})--" << '\n';
    ss << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
    ss << R"--(\setlength\extrarowheight{2pt}\renewcommand{\arraystretch}{1.5})--" << '\n';
    ss << R"--(\begin{tabular}{|c|c|c|c|c|c|})--" << '\n';
    ss << R"--(\hline)--" << '\n';
    ss << R"--(\multicolumn{6}{|c|}{\parbox{\linewidth}{\centering Resolved Photon Category: Up Variations \\ (\% difference from nominal)}}\\[5 pt])--" << '\n';
    ss << R"--(\hline)--" << '\n';

    ss << R"--({Sample} & EG\_RESOLUTION\_ALL & EG\_SCALE\_ALL & PH\_EFF\_ISO\_Uncertainty & PH\_EFF\_ID\_Uncertainty & PH\_EFF\_TRIGGER\_Uncertainty \\ \hline)--" << '\n';
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    cs << R"--(\section*{Table 23 Prompt and Displaced Signal Samples})--" << '\n';
    cs << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
    cs << R"--(\setlength\extrarowheight{2pt}\renewcommand{\arraystretch}{1.5})--" << '\n';
    cs << R"--(\begin{tabular}{|c|c|c|c|c|c|})--" << '\n';
    cs << R"--(\hline)--" << '\n';
    cs << R"--(\multicolumn{6}{|c|}{\parbox{\linewidth}{\centering Resolved Photon Category: Up Variations \\ (raw number of events)}}\\[5 pt])--" << '\n';
    cs << R"--(\hline)--" << '\n';

    cs << R"--({Sample} & EG\_RESOLUTION\_ALL & EG\_SCALE\_ALL & PH\_EFF\_ISO\_Uncertainty & PH\_EFF\_ID\_Uncertainty & PH\_EFF\_TRIGGER\_Uncertainty \\ \hline)--" << '\n';

    for (auto i = 0; (i < resultmaps.size()); i++)
    {
        float nominalVal = static_cast<float>(resultmaps[i]["nominal"]);

        ss << prefixes[i] << " & ";
        ss << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1up"]) ? ((resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1up"]-nominalVal)/nominalVal)*100.0 : 0.0)
        << " & " << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1up"]) ? ((resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1up"]-nominalVal)/nominalVal)*100.0 : 0.0)
        << " & " << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1up"]) ? ((resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1up"]-nominalVal)/nominalVal)*100.0 : 0.0)
        << " & " << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1up"]) ? ((resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1up"]-nominalVal)/nominalVal)*100.0 : 0.0)
        << " & " << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up"]) ? ((resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up"]-nominalVal)/nominalVal)*100.0 : 0.0)
        << R"--( \\ \hline)--" << '\n';
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        cs << prefixes[i] << " & ";
        cs << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1up"]) ? resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1up"] : 0.0)
        << " & " << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1up"]) ? resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1up"] : 0.0)
        << " & " << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1up"]) ? resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1up"] : 0.0)
        << " & " << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1up"]) ? resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1up"] : 0.0)
        << " & " << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up"]) ? resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up"] : 0.0)
        << R"--( \\ \hline)--" << '\n';


        if (!resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1up"])
        {
            std::cout << i << ": photon_iso_eff:PH_EFF_ISO_Uncertainty__1up not found\n";
        }
        else
        {
            std::cout << i << ':' << ' ' << resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1up"] << ' ' << nominalVal << '\n';
        }
        if (!resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1up"])
        {
            std::cout << i << ": photon_id_eff:PH_EFF_ID_Uncertainty__1up not found\n";
        }
        else
        {
            std::cout << i << ':' << ' ' << resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1up"] << ' ' << nominalVal << '\n';
        }
        if (!resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up"])
        {
            std::cout << i << ": photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up not found\n";
        }
        else
        {
            std::cout << i << ':' << ' ' << resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up"] << ' ' << nominalVal << '\n';
        }
    }

    ss << R"--(\end{tabular}})--" << '\n';
    cs << R"--(\end{tabular}})--" << '\n';

    ss << "\n\n\n";
    cs << "\n\n\n";

    std::cout << "\n\n\n" << ss.str();
    std::cout << "\n\n\n" << cs.str();
}

void Table24_Displaced_Axions()
{
    Event::systematics =
    {
        "EG_RESOLUTION_ALL__1down",
        "EG_SCALE_ALL__1down",
        "PH_EFF_ISO_Uncertainty__1down",
        "PH_EFF_ID_Uncertainty__1down",
        "PH_EFF_TRIGGER_Uncertainty__1down",
//        "EG_RESOLUTION_ALL__1down",
    };

    std::vector<std::vector<std::string>> input_filenames = {
        //Signal
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root"},
        //Displaced Signal
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_FewMassPoints.root"},
    };

    constexpr std::array<const char*,15> prefixes = {R"--(Prompt Signal $m_{\text{A}}$ = 1 GeV)--", R"--(Prompt Signal $m_{\text{A}}$ = 5 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 0.2 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 0.5 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 1 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 3 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 5 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 10 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 10.1 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 20 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 20.1 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 29.5 GeV)--", };

    std::vector<float> massPoints = {0.2,0.5,1,3,5,10,10.1,20,20.1,29.5};
    std::vector<ROOT::RDF::RResultHandle> Totals, Nodes;
    std::vector<RResultMap<float>> resultmaps;
    std::stringstream ss, cs;

    int counter = 0;

    for (auto& i: input_filenames)
    {
        SchottDataFrame df(MakeRDF(i,8));

        auto EventWeight = df.Define("EventWeight",
        [](const RVec<float>& ei_event_weights_generator)
        {
            return  ((ei_event_weights_generator[0]) ? 1 / ei_event_weights_generator[0] : 1);

        }, {"ei_event_weights_generator"}).Define("truth_axions", [&](RVec<TruthParticle> truth_particles)
        {
            truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
            [](TruthParticle& x)
            {
             return (abs(x.mc_pdg_id) != 36);

            }), truth_particles.end());

            return truth_particles;

        }, {"truth_particles"}).Define("axion_masses", [&](RVec<TruthParticle>& truth_axions)
        {
            return truth_axions[0].mc_mass/1e3f;

        }, {"truth_axions"});

        auto trigger_selection = df.Define("truth_axions", [&](RVec<TruthParticle> truth_particles)
        {
            truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
            [](TruthParticle& x)
            {
               return (abs(x.mc_pdg_id) != 36);

            }), truth_particles.end());

            return truth_particles;

        }, {"truth_particles"}).Define("axion_masses", [&](RVec<TruthParticle>& truth_axions)
        {
            return truth_axions[0].mc_mass/1e3f;

        }, {"truth_axions"}).Filter(
        [](const RVec<std::string>& trigger_passed_triggers)
        {
            bool trigger_found = (std::find_first_of(trigger_passed_triggers.begin(), trigger_passed_triggers.end(), triggers.begin(), triggers.end()) != trigger_passed_triggers.end());

            if (!trigger_found)
            {
                return false;
            }
            return true;

        }, {"trigger_passed_triggers"});

        auto two_leptons = trigger_selection.Define("di_electrons",
        [](RVec<AbstractParticle> electrons)
        {
            electrons.erase(std::remove_if(electrons.begin(),electrons.end(),
            [](AbstractParticle& ep)
            {
                //keep electron if it has pt > 20 GeV, |eta| < 2.37, not (1.37 < |eta| < 1.52), and id_medium = 1

                return (!((ep.electron_pt/1e3 > 20) && (abs(ep.electron_eta) < 2.37) &&
                (!((1.37 < abs(ep.electron_eta)) && (abs(ep.electron_eta) < 1.52)))
                && (ep.electron_id_medium == 1)));

            }), electrons.end());

            return electrons;

        },{"abstract_electrons"}).Filter([](const RVec<Muon>& muons, RVec<AbstractParticle> electrons)
        {
            return (electrons.size()==2 && muons.empty());

        }, {"muons", "di_electrons"});

        auto opp_charge = two_leptons.Filter([](const RVec<AbstractParticle>& electrons)
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

        auto photon_passes_cuts = ptCut.Define("abstract_photons_pass_cut_indices",[&](RVec<AbstractParticle>& p) //p = photon
        {
            RVec<int> x; //indices of photons that pass cuts

            for (auto i = 0; i < p.size(); i++)
            {
                if (not ((abs(p[i].photon_eta) >= 2.37) or (abs(p[i].photon_eta) > 1.37 and abs(p[i].photon_eta) < 1.52) or (not p[i].photon_id_loose)))
                {
                    x.push_back(i);
                }
            }

            return x;
        }, {"abstract_photons"}).Define("photons_pass_cuts",
        [&](RVec<AbstractParticle>& photons, RVec<int>& x)
        {
            return Take(photons, x); //Taking only the photons that passed the cuts in each event

        }, {"abstract_photons", "abstract_photons_pass_cut_indices"});

        auto resolved_reco_photons_matched = photon_passes_cuts.Define("reco_photons_matched_indices",
        [](const RVec<AbstractParticle>& reco_photons_matched)
        {
            RVec<unsigned long> x; //vector of indices
            if (reco_photons_matched.size() < 2)
            {
                return x;
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
                if (i==0 || ((abs(1-X) < abs(1-best_X)) and (delta_r < 1.5)))
                {
                    best_X = X;
                    pt1 = reco_photons_matched[combs[0][i]].photon_pt;
                    pt2 = reco_photons_matched[combs[1][i]].photon_pt;
                    chosen_delta_r = delta_r;
                    x = {combs[0][i], combs[1][i]};
                }
            }
            if (pt1 > 10e3 && pt2 > 10e3 && best_X > 0.96 && best_X < 1.2 && chosen_delta_r < 1.5)
            {
                return x;
            }

            x.clear();
            return x;

        }, {"photons_pass_cuts"}).Filter(
        [&](RVec<unsigned long>& reco_photons_matched_indices)
        {
            return (reco_photons_matched_indices.size()==2);

        }, {"reco_photons_matched_indices"}).Define("chosen_two",
        [](const RVec<AbstractParticle>& reco_photons_matched, RVec<unsigned long>& reco_photons_matched_indices)
        {
            return Take(reco_photons_matched, reco_photons_matched_indices);

        }, {"photons_pass_cuts", "reco_photons_matched_indices"}).Define("reconstructed_mass",
        [&](RVec<AbstractParticle>& diph)
        {
            return (diph[0].PhotonVector()+diph[1].PhotonVector()).M()/1e3;

        }, {"chosen_two"});

        auto SB = resolved_reco_photons_matched.Filter(
        [](const double reconstructed_mass)
        {
            return (reconstructed_mass < 110) || (reconstructed_mass > 140);
        }, {"reconstructed_mass"});

        auto SR = resolved_reco_photons_matched.Filter(
        [](const double reconstructed_mass)
        {
            return (reconstructed_mass >= 110) && (reconstructed_mass <= 140);
        }, {"reconstructed_mass"});

        auto totEventWeight = resolved_reco_photons_matched //rpmi = reco_photons_matched_indices
        .Define("totEventWeight", [](RVec<int>& x, RVec<unsigned long>& rpmi, RVec<float> photon_id_eff, RVec<float> photon_iso_eff, RVec<float> photon_trg_eff/*, RVec<float> ei_event_weights_generator*/)
        {
            //   ||  jic they don't already have the same size  ||
            //   \/                                             \/
            auto ResizeVal = std::max({photon_id_eff.size(), photon_iso_eff.size(), photon_trg_eff.size()});
            photon_id_eff.resize(ResizeVal,1);
            photon_iso_eff.resize(ResizeVal,1);
            photon_trg_eff.resize(ResizeVal,1);

            //First, take the x indices from the respective photon_*_eff vectors,
            //this corresponds to the photons from the set of all reco photons in
            //the events that passed the "photon_passes_cuts" cuts. Then, from
            //those photons, select the indices rpmi elements that corresponds to
            //the two resolved reco-photons

            RVec<float> two_efficiencies = Take(Take(photon_id_eff,x),rpmi) * Take(Take(photon_iso_eff,x),rpmi) * Take(Take(photon_trg_eff,x),rpmi);//*ei_event_weights_generator[0];

            float total = 1.0f;

            for (auto i: two_efficiencies)
            {
                total *= i;
            }

            return total;

        }, {"abstract_photons_pass_cut_indices", "reco_photons_matched_indices", "photon_id_eff", "photon_iso_eff", "photon_trg_eff",/* "ei_event_weights_generator"*/});

        if (counter < 2)
        {
            Totals.push_back(EventWeight.Sum<float>("EventWeight"));
            resultmaps.push_back(VariationsFor(totEventWeight.Sum<float>("totEventWeight")));
        }
        else
        {
            for (auto& mass_point: massPoints)
            {
                auto mass_point_EventWeight = EventWeight.Filter([&]
                (float axion_mass)
                {
                    return (roundToOneDecimalPlace(axion_mass) == mass_point);
//                    return ((axion_mass > mass_point - 0.1f) and (axion_mass < mass_point + 0.1f));
                }, {"axion_masses"});

                auto mass_point_totEventWeight = totEventWeight.Filter([&]
                (float axion_mass)
                {
                    return (roundToOneDecimalPlace(axion_mass) == mass_point);
//                    return ((axion_mass > mass_point - 0.1f) and (axion_mass < mass_point + 0.1f));
                }, {"axion_masses"});

                Totals.push_back(mass_point_EventWeight.Sum<float>("EventWeight"));
                resultmaps.push_back(VariationsFor(mass_point_totEventWeight.Sum<float>("totEventWeight")));
                Nodes.push_back(mass_point_totEventWeight.Take<UInt_t, RVec<UInt_t>>("ei_event_number"));
            }
        }

        counter++;

    }

    ROOT::RDF::RunGraphs(Totals); // running all computation nodes concurrently

    std::cout << "Resolved Category Down\n======================\n";
    for (int i = 0; i < Nodes.size(); i++)
    {
        std::cout << prefixes[i+2] << '\n';
        std::cout << "==========\n";
        for (auto& j: *Nodes[i].GetResultPtr<RVec<UInt_t>>())
        {
            std::cout << j << ", ";
        }
        std::cout << "\n\n";
    }


    ss << R"--(\section*{Table 24 Prompt and Displaced Signal Samples})--" << '\n';
    ss << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
    ss << R"--(\setlength\extrarowheight{2pt}\renewcommand{\arraystretch}{1.5})--" << '\n';
    ss << R"--(\begin{tabular}{|c|c|c|c|c|c|})--" << '\n';
    ss << R"--(\hline)--" << '\n';
    ss << R"--(\multicolumn{6}{|c|}{\parbox{\linewidth}{\centering Resolved Photon Category: Down Variations \\ (\% difference from nominal)}}\\[5 pt])--" << '\n';
    ss << R"--(\hline)--" << '\n';

    ss << R"--({Sample} & EG\_RESOLUTION\_ALL & EG\_SCALE\_ALL & PH\_EFF\_ISO\_Uncertainty & PH\_EFF\_ID\_Uncertainty & PH\_EFF\_TRIGGER\_Uncertainty \\ \hline)--" << '\n';
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    cs << R"--(\section*{Table 24 Prompt and Displaced Signal Samples})--" << '\n';
    cs << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
    cs << R"--(\setlength\extrarowheight{2pt}\renewcommand{\arraystretch}{1.5})--" << '\n';
    cs << R"--(\begin{tabular}{|c|c|c|c|c|c|})--" << '\n';
    cs << R"--(\hline)--" << '\n';
    cs << R"--(\multicolumn{6}{|c|}{\parbox{\linewidth}{\centering Resolved Photon Category: Down Variations \\ (raw number of events)}}\\[5 pt])--" << '\n';
    cs << R"--(\hline)--" << '\n';

    cs << R"--({Sample} & EG\_RESOLUTION\_ALL & EG\_SCALE\_ALL & PH\_EFF\_ISO\_Uncertainty & PH\_EFF\_ID\_Uncertainty & PH\_EFF\_TRIGGER\_Uncertainty \\ \hline)--" << '\n';

    for (auto i = 0; (i < resultmaps.size()); i++)
    {
        float nominalVal = static_cast<float>(resultmaps[i]["nominal"]);

        ss << prefixes[i] << " & ";
        ss << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1down"]) ? ((resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1down"]-nominalVal)/nominalVal)*100.0 : 0.0)
        << " & " << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1down"]) ? ((resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1down"]-nominalVal)/nominalVal)*100.0 : 0.0)
        << " & " << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1down"]) ? ((resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1down"]-nominalVal)/nominalVal)*100.0 : 0.0)
        << " & " << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1down"]) ? ((resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1down"]-nominalVal)/nominalVal)*100.0 : 0.0)
        << " & " << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down"]) ? ((resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down"]-nominalVal)/nominalVal)*100.0 : 0.0)
        << R"--( \\ \hline)--" << '\n';
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        cs << prefixes[i] << " & ";
        cs << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1down"]) ? resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1down"] : 0.0)
        << " & " << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1down"]) ? resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1down"] : 0.0)
        << " & " << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1down"]) ? resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1down"] : 0.0)
        << " & " << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1down"]) ? resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1down"] : 0.0)
        << " & " << std::setprecision(4) << std::fixed
        << ( (nominalVal && resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down"]) ? resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down"] : 0.0)
        << R"--( \\ \hline)--" << '\n';


        if (!resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1down"])
        {
            std::cout << i << ": photon_iso_eff:PH_EFF_ISO_Uncertainty__1down not found\n";
        }
        else
        {
            std::cout << i << ':' << ' ' << resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1down"] << ' ' << nominalVal << '\n';
        }
        if (!resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1down"])
        {
            std::cout << i << ": photon_id_eff:PH_EFF_ID_Uncertainty__1down not found\n";
        }
        else
        {
            std::cout << i << ':' << ' ' << resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1down"] << ' ' << nominalVal << '\n';
        }
        if (!resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down"])
        {
            std::cout << i << ": photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down not found\n";
        }
        else
        {
            std::cout << i << ':' << ' ' << resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down"] << ' ' << nominalVal << '\n';
        }
    }

    ss << R"--(\end{tabular}})--" << '\n';
    cs << R"--(\end{tabular}})--" << '\n';

    ss << "\n\n\n";
    cs << "\n\n\n";

    std::cout << "\n\n\n" << ss.str();
    std::cout << "\n\n\n" << cs.str();
}

void Section7_TablesPlots()
{
    
    auto start_time = Clock::now();
//    auto verbosity = ROOT::Experimental::RLogScopedVerbosity(ROOT::Detail::RDF::RDFLogChannel(), ROOT::Experimental::ELogLevel::kInfo);
//    Table21();
//    Table22();
//    Table23();
//    Table24();
    Table21_Displaced_Axions();
//    Table22_Displaced_Axions();
//    Table23_Displaced_Axions();
//    Table24_Displaced_Axions();
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



