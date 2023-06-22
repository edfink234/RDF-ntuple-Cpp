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

//void Table21()
//{
//    Event::systematics =
//    {
////        "PH_EFF_ISO_Uncertainty",
////        "PH_EFF_ISO_Uncertainty",
////        "EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR",
////        "PRW_DATASF",
////        "MUON_EFF_RECO_SYS",
////        "MUON_EFF_ISO_SYS",
////        "MUON_EFF_TrigSystUncertainty",
////        "EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR",
////        "MUON_EFF_TrigStatUncertainty",
////        "MUON_EFF_RECO_STAT",
////        "MUON_EFF_TTVA_STAT",
////        "EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR",
////        "EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR",
////        "MUON_EFF_TTVA_SYS",
////        "MUON_EFF_ISO_STAT",
////        "MUON_SAGITTA_RHO",
////        "EG_RESOLUTION_ALL",
////        "EG_SCALE_ALL",
////        "MUON_MS",
////        "MUON_ID",
////        "EL_EFF_TriggerEff_TOTAL_1NPCOR_PLUS_UNCOR",
////        "MUON_SAGITTA_RESBIAS",
////        "MUON_SCALE",
//
////        "PH_EFF_ISO_Uncertainty__1down",
////        "EG_SCALE_ALL__1down",
////        "PH_EFF_ID_Uncertainty__1down",
////        "PH_EFF_TRIGGER_Uncertainty__1down",
//        "EG_RESOLUTION_ALL__1up",
//        "EG_SCALE_ALL__1up",
////        "PH_EFF_ISO_Uncertainty__1up",
//        "PH_EFF_ID_Uncertainty__1up",
////        "PH_EFF_TRIGGER_Uncertainty__1up",
//    };
//
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
//    std::array<double,9> JetNumeratorSFs = {((139e15)*(1.9828e-9)*(0.821204)),((139e15)*(110.64e-12)*(0.69275)),((139e15)*(40.645e-12)*(0.615906)),((139e15)*(1.9817e-9)*(0.1136684)),((139e15)*(110.47e-12)*(0.1912956)),((139e15)*(40.674e-12)*(0.2326772)),((139e15)*(1.9819e-9)*(0.0656969)),((139e15)*(110.53e-12)*(0.1158741)),((139e15)*(40.68e-12)*(0.1535215))};
//    std::array<double,3> SFs = {((139e15)*(.871e-12)),((139e15)*(.199e-12)), ((139e15)*(.0345e-15))}; //numerators for Z-gamma bkg
//
//    std::vector<std::string> prefixes = { R"--(pty2\_9\_17)--", R"--(pty\_17\_myy\_0\_80)--", R"--(pty\_17\_myy\_80)--", "data", R"--($\text{Sig } m_{A}$ = 1 GeV)--", R"--($\text{Sig } m_{A}$ = 2 GeV)--", R"--($\text{Sig } m_{A}$ = 3 GeV)--", R"--($\text{Sig } m_{A}$ = 5 GeV)--", R"--($\text{Sig } m_{A}$ = 9 GeV)--", R"--(Zee\_lightJet\_0-70)--", R"--(Zee\_lightJet\_70-140)--", R"--(Zee\_lightJet\_140-280)--", R"--(Zee\_cJet\_0-70)--", R"--(Zee\_cJet\_70-140)--", R"--(Zee\_cJet\_140-280)--", R"--(Zee\_bJet\_0-70)--", R"--(Zee\_bJet\_70-140)--", R"--(Zee\_bJet\_140-280)--"};
//
//    std::vector<ROOT::RDF::RResultHandle> Totals;
//    std::vector<RResultMap<float>> resultmaps;
//
//    std::stringstream ss;
//
//    for (auto& i: input_filenames)
//    {
//        SchottDataFrame df(MakeRDF(i,8));
//
//        auto EventWeight = df.Define("EventWeight", //mc generator weight
//        [](const RVec<float>& ei_event_weights_generator)
//        {
//            return ((ei_event_weights_generator[0]) ? 1 / ei_event_weights_generator[0] : 1);
//
//        }, {"ei_event_weights_generator"})
//        .Define("totEventWeightVec", //vector of efficiencies, product of photon_id_eff, photon_iso_eff, and photon_trg_eff
//        [](RVec<float> photon_id_eff, RVec<float> photon_iso_eff, RVec<float> photon_trg_eff/*, RVec<float> ei_event_weights_generator*/)
//        {
//            auto ResizeVal = std::max({photon_id_eff.size(), photon_iso_eff.size(), photon_trg_eff.size()});
//            photon_id_eff.resize(ResizeVal,1);
//            photon_iso_eff.resize(ResizeVal,1);
//            photon_trg_eff.resize(ResizeVal,1);
//
//            return photon_id_eff;
//            //*photon_iso_eff*photon_trg_eff; //Returns a vector of efficiencies
//
//        }, {"photon_id_eff", "photon_iso_eff", "photon_trg_eff",/* "ei_event_weights_generator"*/});
//
//        auto two_leptons = EventWeight
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
//        [](RVec<AbstractParticle> electrons)
//        {
//            //keep the electrons in each event that have pt > 20 GeV, |η| < 2.37,
//            //|η| not between 1.37 and 1.52, and that satisfy a medium id criteria `electron_id_medium`
//            electrons.erase(std::remove_if(electrons.begin(),electrons.end(),
//            [](AbstractParticle& ep)
//            {
//                 return (!((ep.electron_pt/1e3 > 20) && (std::abs(ep.electron_eta) < 2.37) &&
//                 (!((1.37 < std::abs(ep.electron_eta)) && (std::abs(ep.electron_eta) < 1.52)))
//                 && (ep.electron_id_medium == 1)));
//
//            }), electrons.end());
//
//            return electrons;
//
//        },{"abstract_electrons"}).Filter(
//        [](RVec<Muon>& muons, RVec<AbstractParticle>& electrons)
//        {
//            return (electrons.size()==2 && muons.empty()); //keep events which have exactly 2 electrons for di_electrons and no muons
//
//        }, {"muons", "di_electrons"});
//
//        //new dataframe node: contains only the events from `two_leptons` whose electrons in the `di_electrons` branch have opposite charge
//        auto opp_charge = two_leptons
//        .Filter([](const RVec<AbstractParticle>& electrons)
//        {
//            return (electrons[0].electron_charge*electrons[1].electron_charge < 0);
//
//        }, {"di_electrons"});
//
//        //new dataframe node: contains only the events from `opp_charge` that have 1 electron with pt > 20 GeV and the other with pt > 27 GeV
//        auto leadingPt = opp_charge.Filter([](const RVec<AbstractParticle>& electrons)
//        {
//            return ((electrons[0].electron_pt > 20e3 && electrons[1].electron_pt > 27e3) || (electrons[1].electron_pt > 20e3 && electrons[0].electron_pt > 27e3));
//        }, {"di_electrons"});
//
//        auto deltaR = leadingPt.Filter([] (const RVec<AbstractParticle>& electrons)
//        {
//            return (DeltaR(electrons[0].ElectronVector(), electrons[1].ElectronVector()) > 0.01);
//        }, {"di_electrons"});
//
//        //new dataframe node: Contains a new column `dilep` in addition to the ones in `same_flavour` that stores the di-electron four-vector
//        auto dilep = deltaR.Define("dilep",[] (RVec<AbstractParticle>& electrons)
//        {
//            return (electrons[0].ElectronVector() + electrons[1].ElectronVector());
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
//        auto photon_passes_cuts = ptCut
//        .Define("photons_pass_cut_indices",
//        [&](const RVec<AbstractParticle>& photons)
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
//        }, {"abstract_photons"})
//        .Define("photons_pass_cuts",
//        [&](RVec<AbstractParticle>& photons, RVec<int>& photon_indices)
//        {
//            return Take(photons, photon_indices);
//        }, {"abstract_photons", "photons_pass_cut_indices"});
//
//        auto merged_reco_photons_matched = photon_passes_cuts.Filter(
//        [&](const RVec<AbstractParticle>& reco_photons_matched)
//        {
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
//                delta_r = DeltaR(reco_photons_matched[combs[0][i]].PhotonVector(), reco_photons_matched[combs[1][i]].PhotonVector());
//                m = (reco_photons_matched[combs[0][i]].PhotonVector() + reco_photons_matched[combs[1][i]].PhotonVector()).M();
//                pt = (reco_photons_matched[combs[0][i]].PhotonVector() + reco_photons_matched[combs[1][i]].PhotonVector()).Pt();
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
//        .Define("merged_photon_index", //new column: consists of the index corresponding to the photon that made the event be classified as merged
//        [&](const RVec<AbstractParticle>& rpm) //rpm = reco photons matched
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
//        [&](const RVec<AbstractParticle>& reco_photons_matched, int merged_photon_index)
//        {
//            return reco_photons_matched[merged_photon_index];
//
//        }, {"photons_pass_cuts", "merged_photon_index"});
//
//        auto dilepton_and_photon = merged_reco_photons_matched
//        .Define("reconstructed_mass",[&](const RVec<AbstractParticle>& di_electrons, const AbstractParticle& merged_photon)
//        {
//            auto four_momentum = di_electrons[0].ElectronVector() + di_electrons[1].ElectronVector();
//
//            return (four_momentum + merged_photon.PhotonVector()).M()/1e3;
//
//        }, {"di_electrons", "merged_photon"});
//
//        auto pSB = dilepton_and_photon.Filter(
//        [](const double reconstructed_mass)
//        {
//            return (reconstructed_mass < 110) || (reconstructed_mass > 130);
//        }, {"reconstructed_mass"});
//
//        auto pSR = dilepton_and_photon.Filter(
//        [](const double reconstructed_mass)
//        {
//            return (reconstructed_mass >= 110) && (reconstructed_mass <= 130);
//        }, {"reconstructed_mass"});
//
//        auto SB = pSB.Filter(
//        [](const RVec<float>& Eratio)
//        {
//            return (!Any(Eratio <= 0.8));
//        }, {"photon_shower_shape_e_ratio"});
//
//        auto SR = pSR.Filter(
//        [](const RVec<float>& Eratio)
//        {
//            return (!Any(Eratio <= 0.8));
//        }, {"photon_shower_shape_e_ratio"});
//
//        auto totEventWeight = merged_reco_photons_matched
//        .Define("totEventWeight", [](RVec<float> totEventWeightVec, RVec<int>& photons_pass_cut_indices, int mpi)
//        {
//            return Take(totEventWeightVec, photons_pass_cut_indices)[mpi];
//
//        }, {"totEventWeightVec", "photons_pass_cut_indices", "merged_photon_index"});
//
//        Totals.push_back(EventWeight.Sum<float>("EventWeight"));
//
//        resultmaps.push_back(VariationsFor(totEventWeight.Sum<float>("totEventWeight")));
//    }
//
////EG_RESOLUTION_ALL__1down:           42864.7 ID 1.00524  ISO 1.02283  TRIG 1.04157
////EG_RESOLUTION_ALL__1up:             42751.5 ID 1.00524  ISO 1.02283  TRIG 1.04157
////EG_SCALE_ALL__1down:                42521.1 ID 1.00524  ISO 1.02283  TRIG 1.04157
////EG_SCALE_ALL__1up:                  43059.8 ID 1.00524  ISO 1.02283  TRIG 1.04157
////PH_EFF_ID_Uncertainty__1down:       42790.4 ID 0.992328 ISO 1.02283  TRIG 1.04157
////PH_EFF_ID_Uncertainty__1up:         42790.4 ID 1.01816  ISO 1.02283  TRIG 1.04157
////PH_EFF_ISO_Uncertainty__1down:      42790.4 ID 1.00524  ISO 1.00971  TRIG 1.04157
////PH_EFF_ISO_Uncertainty__1up:        42790.4 ID 1.00524  ISO 1.03595  TRIG 1.04157
////PH_EFF_TRIGGER_Uncertainty__1down:  42790.4 ID 1.00524  ISO 1.02283  TRIG 1.01646
////PH_EFF_TRIGGER_Uncertainty__1up:    42790.4 ID 1.00524  ISO 1.02283  TRIG 1.06668
//
//    ROOT::RDF::RunGraphs(Totals); // running all computation nodes concurrently
//
////    for (auto& i: GeneratorWeightCounts)
////    {
////        std::cout << *i << '\n';
////    }
////    int count = 0;
////    for (auto& i: resultmaps)
////    {
////        std::cout << count++ << "\n==\n";
////        for (auto& j: i.GetKeys())
////        {
////            std::cout << j << '\n';
////        }
////        std::cout << '\n';
////    }
////          resultmaps
////          ----------
////    0       1       //Z-gamma
////    2       3       //Z-gamma
////    4       5       //Z-gamma
////    6       7       //data
////    8       9       //ma1
////    10      11      //ma2
////    12      13      //ma3
////    14      15      //ma5
////    16      17      //ma9
////    18      19      //Z-jets
////    20      21      //Z-jets
////    22      23      //Z-jets
////    24      25      //Z-jets
////    26      27      //Z-jets
////    28      29      //Z-jets
////    30      31      //Z-jets
////    32      33      //Z-jets
////    34      35      //Z-jets
//
//    ss << R"--(\section*{Table 21})--" << '\n';
//    ss << R"--(\hspace{-3cm}\scalebox{0.8}{)--" << '\n';
//    ss << R"--(\setlength\extrarowheight{2pt}\renewcommand{\arraystretch}{1.5})--" << '\n';
////    ss << R"--(\begin{tabular}{|c|c|c|c|c|c|})--" << '\n';
//    ss << R"--(\begin{tabular}{|c|c|c|c|})--" << '\n';
//
//    ss << R"--(\hline)--" << '\n';
//
//        ss << R"--(\multicolumn{4}{|c|}{\parbox{\linewidth}{\centering Merged Photon Category: Up Variations \\ (\% difference from nominal)}}\\[5 pt])--" << '\n';
//    ss << R"--(\hline)--" << '\n';
//
////    ss << R"--({Sample} & EG\_RESOLUTION\_ALL & EG\_SCALE\_ALL & PH\_EFF\_ISO\_Uncertainty & PH\_EFF\_ID\_Uncertainty & PH\_EFF\_TRIGGER\_Uncertainty \\ \hline)--" << '\n';
//    ss << R"--({Sample} & EG\_RESOLUTION\_ALL & EG\_SCALE\_ALL & PH\_EFF\_ID\_Uncertainty \\ \hline)--" << '\n';
//
//
//    double finalScaleVal;
//
//    double ZgammaNominal = 0, ZgammaEG_RESOLUTION_ALL = 0,
//    ZgammaEG_SCALE_ALL = 0, ZgammaPH_EFF_ISO_Uncertainty = 0,
//    ZgammaPH_EFF_ID_Uncertainty = 0, ZgammaPH_EFF_TRIGGER_Uncertainty = 0;
//
//    double ZjetsNominal = 0, ZjetsEG_RESOLUTION_ALL = 0,
//    ZjetsEG_SCALE_ALL = 0, ZjetsPH_EFF_ISO_Uncertainty = 0,
//    ZjetsPH_EFF_ID_Uncertainty = 0, ZjetsPH_EFF_TRIGGER_Uncertainty = 0;
//
//    double totbkgNominal = 0, totbkgEG_RESOLUTION_ALL = 0,
//    totbkgEG_SCALE_ALL = 0, totbkgPH_EFF_ISO_Uncertainty = 0,
//    totbkgPH_EFF_ID_Uncertainty = 0, totbkgPH_EFF_TRIGGER_Uncertainty = 0;
//
//    for (auto i = 0; (i < resultmaps.size()); i++)
//    {
//        double nominalVal = static_cast<double>(resultmaps[i]["nominal"]);
//
////        auto denominator = *Totals[i].GetResultPtr<ULong64_t>();
//        auto denominator = *Totals[i].GetResultPtr<float>();
//
//        if (i >= 0 && i <= 2) //Zgamma
//        {
//            finalScaleVal = SFs[i]/denominator;
//            ZgammaNominal += finalScaleVal*nominalVal;
//
//            ZgammaEG_RESOLUTION_ALL += finalScaleVal*resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1up"];
//            ZgammaEG_SCALE_ALL += finalScaleVal*resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1up"];
////            ZgammaPH_EFF_ISO_Uncertainty += finalScaleVal*resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1up"];
//            ZgammaPH_EFF_ID_Uncertainty += finalScaleVal*resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1up"];
////            ZgammaPH_EFF_TRIGGER_Uncertainty += finalScaleVal*resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up"];
//        }
//
//        else if (i >= 9)
//        {
//            finalScaleVal = JetNumeratorSFs[i-9]/denominator;
//            ZjetsNominal += finalScaleVal*nominalVal;
//
//            ZjetsEG_RESOLUTION_ALL += finalScaleVal*resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1up"];
//            ZjetsEG_SCALE_ALL += finalScaleVal*resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1up"];
////            ZjetsPH_EFF_ISO_Uncertainty += finalScaleVal*resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1up"];
//            ZjetsPH_EFF_ID_Uncertainty += finalScaleVal*resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1up"];
////            ZjetsPH_EFF_TRIGGER_Uncertainty += finalScaleVal*resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up"];
//
//        }
//
//        ss << prefixes[i] << " & ";
//        ss << std::setprecision(4) << std::fixed
//        << ( (nominalVal && resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1up"]) ? ((resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1up"]-nominalVal)/nominalVal)*100.0 : 0.0)
//        << " & " << std::setprecision(4) << std::fixed
//        << ( (nominalVal && resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1up"]) ? ((resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1up"]-nominalVal)/nominalVal)*100.0 : 0.0)
////        << " & " << std::setprecision(4) << std::fixed
////        << ( (nominalVal && resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1up"]) ? ((resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1up"]-nominalVal)/nominalVal)*100.0 : 0.0)
//        << " & " << std::setprecision(4) << std::fixed
//        << ( (nominalVal && resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1up"]) ? ((resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1up"]-nominalVal)/nominalVal)*100.0 : 0.0)
////        << " & " << std::setprecision(4) << std::fixed
////        << ( (nominalVal && resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up"]) ? ((resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up"]-nominalVal)/nominalVal)*100.0 : 0.0)
//        << R"--( \\ \hline)--" << '\n';
//
//        std::cout << "Keys for index " << i << " = ";
//        for (auto& j: resultmaps[i].GetKeys())
//        {
//            std::cout << j << ", ";
//        }
//        std::cout << "\n\n";
//    }
//
//    ss << R"--(Total $Z\gamma\gamma$ & )--";
//    ss << std::setprecision(4) << std::fixed
//    << ( (ZgammaNominal && ZgammaEG_RESOLUTION_ALL) ? ((ZgammaEG_RESOLUTION_ALL-ZgammaNominal)/ZgammaNominal)*100.0 : 0.0)
//    << " & " << std::setprecision(4) << std::fixed
//    << ( (ZgammaNominal && ZgammaEG_SCALE_ALL) ? ((ZgammaEG_SCALE_ALL-ZgammaNominal)/ZgammaNominal)*100.0 : 0.0)
////    << " & " << std::setprecision(4) << std::fixed
////    << ( (ZgammaNominal && ZgammaPH_EFF_ISO_Uncertainty) ? ((ZgammaPH_EFF_ISO_Uncertainty-ZgammaNominal)/ZgammaNominal)*100.0 : 0.0)
//    << " & " << std::setprecision(4) << std::fixed
//    << ( (ZgammaNominal && ZgammaPH_EFF_ID_Uncertainty) ? ((ZgammaPH_EFF_ID_Uncertainty-ZgammaNominal)/ZgammaNominal)*100.0 : 0.0)
////    << " & " << std::setprecision(4) << std::fixed
////    << ( (ZgammaNominal && ZgammaPH_EFF_TRIGGER_Uncertainty) ? ((ZgammaPH_EFF_TRIGGER_Uncertainty-ZgammaNominal)/ZgammaNominal)*100.0 : 0.0)
//    << R"--( \\ \hline)--" << '\n';
//
//    ss << R"--(Total $Z$+jets & )--";
//    ss << std::setprecision(4) << std::fixed
//    << ( (ZjetsNominal && ZjetsEG_RESOLUTION_ALL) ? ((ZjetsEG_RESOLUTION_ALL-ZjetsNominal)/ZjetsNominal)*100.0 : 0.0)
//    << " & " << std::setprecision(4) << std::fixed
//    << ( (ZjetsNominal && ZjetsEG_SCALE_ALL) ? ((ZjetsEG_SCALE_ALL-ZjetsNominal)/ZjetsNominal)*100.0 : 0.0)
////    << " & " << std::setprecision(4) << std::fixed
////    << ( (ZjetsNominal && ZjetsPH_EFF_ISO_Uncertainty) ? ((ZjetsPH_EFF_ISO_Uncertainty-ZjetsNominal)/ZjetsNominal)*100.0 : 0.0)
//    << " & " << std::setprecision(4) << std::fixed
//    << ( (ZjetsNominal && ZjetsPH_EFF_ID_Uncertainty) ? ((ZjetsPH_EFF_ID_Uncertainty-ZjetsNominal)/ZjetsNominal)*100.0 : 0.0)
////    << " & " << std::setprecision(4) << std::fixed
////    << ( (ZjetsNominal && ZjetsPH_EFF_TRIGGER_Uncertainty) ? ((ZjetsPH_EFF_TRIGGER_Uncertainty-ZjetsNominal)/ZjetsNominal)*100.0 : 0.0)
//    << R"--( \\ \hline)--" << '\n';
//
//    totbkgNominal = ZgammaNominal + ZjetsNominal;
//
//    totbkgEG_RESOLUTION_ALL = ZgammaEG_RESOLUTION_ALL + ZjetsEG_RESOLUTION_ALL;
//    totbkgEG_SCALE_ALL = ZgammaEG_SCALE_ALL + ZjetsEG_SCALE_ALL;
//    totbkgPH_EFF_ISO_Uncertainty = ZgammaPH_EFF_ISO_Uncertainty + ZjetsPH_EFF_ISO_Uncertainty;
//    totbkgPH_EFF_ID_Uncertainty = ZgammaPH_EFF_ID_Uncertainty + ZjetsPH_EFF_ID_Uncertainty;
//    totbkgPH_EFF_TRIGGER_Uncertainty = ZgammaPH_EFF_TRIGGER_Uncertainty + ZjetsPH_EFF_TRIGGER_Uncertainty;
//
//    ss << R"--(Total Bkg & )--";
//    ss << std::setprecision(4) << std::fixed
//    << ( (totbkgNominal && totbkgEG_RESOLUTION_ALL) ? ((totbkgEG_RESOLUTION_ALL-totbkgNominal)/totbkgNominal)*100.0 : 0.0)
//    << " & " << std::setprecision(4) << std::fixed
//    << ( (totbkgNominal && totbkgEG_SCALE_ALL) ? ((totbkgEG_SCALE_ALL-totbkgNominal)/totbkgNominal)*100.0 : 0.0)
////    << " & " << std::setprecision(4) << std::fixed
////    << ( (totbkgNominal && totbkgPH_EFF_ISO_Uncertainty) ? ((totbkgPH_EFF_ISO_Uncertainty-totbkgNominal)/totbkgNominal)*100.0 : 0.0)
//    << " & " << std::setprecision(4) << std::fixed
//    << ( (totbkgNominal && totbkgPH_EFF_ID_Uncertainty) ? ((totbkgPH_EFF_ID_Uncertainty-totbkgNominal)/totbkgNominal)*100.0 : 0.0)
////    << " & " << std::setprecision(4) << std::fixed
////    << ( (totbkgNominal && totbkgPH_EFF_TRIGGER_Uncertainty) ? ((totbkgPH_EFF_TRIGGER_Uncertainty-totbkgNominal)/totbkgNominal)*100.0 : 0.0)
//    << R"--( \\ \hline)--" << '\n';
//
//    ss << R"--(\end{tabular}})--" << '\n';
//
//    ss << "\n\n\n";
//
//    std::ofstream out("Table21.txt");
//    out << ss.str() << '\n';
//    out.close();
////    std::cout << ss.str() << '\n';
//
//}

//void Table22()
//{
//    Event::systematics =
//    {
////        "PH_EFF_ISO_Uncertainty",
////        "PH_EFF_ISO_Uncertainty",
////        "EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR",
////        "PRW_DATASF",
////        "MUON_EFF_RECO_SYS",
////        "MUON_EFF_ISO_SYS",
////        "MUON_EFF_TrigSystUncertainty",
////        "EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR",
////        "MUON_EFF_TrigStatUncertainty",
////        "MUON_EFF_RECO_STAT",
////        "MUON_EFF_TTVA_STAT",
////        "EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR",
////        "EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR",
////        "MUON_EFF_TTVA_SYS",
////        "MUON_EFF_ISO_STAT",
////        "MUON_SAGITTA_RHO",
////        "EG_RESOLUTION_ALL",
////        "EG_SCALE_ALL",
////        "MUON_MS",
////        "MUON_ID",
////        "EL_EFF_TriggerEff_TOTAL_1NPCOR_PLUS_UNCOR",
////        "MUON_SAGITTA_RESBIAS",
////        "MUON_SCALE",
////        "EG_RESOLUTION_ALL__1up",
////        "EG_SCALE_ALL__1up",
////        "PH_EFF_ISO_Uncertainty__1up",
////        "PH_EFF_ID_Uncertainty__1up",
////        "PH_EFF_TRIGGER_Uncertainty__1up",
//        "EG_RESOLUTION_ALL__1down",
//        "EG_SCALE_ALL__1down",
////        "PH_EFF_ISO_Uncertainty__1down",
//        "PH_EFF_ID_Uncertainty__1down",
////        "PH_EFF_TRIGGER_Uncertainty__1down",
//    };
//
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
//    std::array<double,9> JetNumeratorSFs = {((139e15)*(1.9828e-9)*(0.821204)),((139e15)*(110.64e-12)*(0.69275)),((139e15)*(40.645e-12)*(0.615906)),((139e15)*(1.9817e-9)*(0.1136684)),((139e15)*(110.47e-12)*(0.1912956)),((139e15)*(40.674e-12)*(0.2326772)),((139e15)*(1.9819e-9)*(0.0656969)),((139e15)*(110.53e-12)*(0.1158741)),((139e15)*(40.68e-12)*(0.1535215))};
//    std::array<double,3> SFs = {((139e15)*(.871e-12)),((139e15)*(.199e-12)), ((139e15)*(.0345e-15))}; //numerators for Z-gamma bkg
//
//    std::vector<std::string> prefixes = { R"--(pty2\_9\_17)--", R"--(pty\_17\_myy\_0\_80)--", R"--(pty\_17\_myy\_80)--", "data", R"--($\text{Sig } m_{A}$ = 1 GeV)--", R"--($\text{Sig } m_{A}$ = 2 GeV)--", R"--($\text{Sig } m_{A}$ = 3 GeV)--", R"--($\text{Sig } m_{A}$ = 5 GeV)--", R"--($\text{Sig } m_{A}$ = 9 GeV)--", R"--(Zee\_lightJet\_0-70)--", R"--(Zee\_lightJet\_70-140)--", R"--(Zee\_lightJet\_140-280)--", R"--(Zee\_cJet\_0-70)--", R"--(Zee\_cJet\_70-140)--", R"--(Zee\_cJet\_140-280)--", R"--(Zee\_bJet\_0-70)--", R"--(Zee\_bJet\_70-140)--", R"--(Zee\_bJet\_140-280)--"};
//
//    std::vector<ROOT::RDF::RResultHandle> Totals;
//    std::vector<RResultMap<float>> resultmaps;
//
//    std::stringstream ss;
//
//    for (auto& i: input_filenames)
//    {
//        SchottDataFrame df(MakeRDF(i,8));
//
//        auto EventWeight = df.Define("EventWeight", //mc generator weight
//        [](const RVec<float>& ei_event_weights_generator)
//        {
//            return ((ei_event_weights_generator[0]) ? 1 / ei_event_weights_generator[0] : 1);
//
//        }, {"ei_event_weights_generator"})
//        .Define("totEventWeightVec", //vector of efficiencies, product of photon_id_eff, photon_iso_eff, and photon_trg_eff
//        [](RVec<float> photon_id_eff, RVec<float> photon_iso_eff, RVec<float> photon_trg_eff/*, RVec<float> ei_event_weights_generator*/)
//        {
//            auto ResizeVal = std::max({photon_id_eff.size(), photon_iso_eff.size(), photon_trg_eff.size()});
//            photon_id_eff.resize(ResizeVal,1);
//            photon_iso_eff.resize(ResizeVal,1);
//            photon_trg_eff.resize(ResizeVal,1);
//
//            return photon_id_eff;
//            //*photon_iso_eff*photon_trg_eff; //Returns a vector of efficiencies
//
//        }, {"photon_id_eff", "photon_iso_eff", "photon_trg_eff",/* "ei_event_weights_generator"*/});
//
//        auto two_leptons = EventWeight
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
//        [](RVec<AbstractParticle> electrons)
//        {
//            //keep the electrons in each event that have pt > 20 GeV, |η| < 2.37,
//            //|η| not between 1.37 and 1.52, and that satisfy a medium id criteria `electron_id_medium`
//            electrons.erase(std::remove_if(electrons.begin(),electrons.end(),
//            [](AbstractParticle& ep)
//            {
//                 return (!((ep.electron_pt/1e3 > 20) && (std::abs(ep.electron_eta) < 2.37) &&
//                 (!((1.37 < std::abs(ep.electron_eta)) && (std::abs(ep.electron_eta) < 1.52)))
//                 && (ep.electron_id_medium == 1)));
//
//            }), electrons.end());
//
//            return electrons;
//
//        },{"abstract_electrons"}).Filter(
//        [](RVec<Muon>& muons, RVec<AbstractParticle>& electrons)
//        {
//            return (electrons.size()==2 && muons.empty()); //keep events which have exactly 2 electrons for di_electrons and no muons
//
//        }, {"muons", "di_electrons"});
//
//        //new dataframe node: contains only the events from `two_leptons` whose electrons in the `di_electrons` branch have opposite charge
//        auto opp_charge = two_leptons
//        .Filter([](const RVec<AbstractParticle>& electrons)
//        {
//            return (electrons[0].electron_charge*electrons[1].electron_charge < 0);
//
//        }, {"di_electrons"});
//
//        //new dataframe node: contains only the events from `opp_charge` that have 1 electron with pt > 20 GeV and the other with pt > 27 GeV
//        auto leadingPt = opp_charge.Filter([](const RVec<AbstractParticle>& electrons)
//        {
//            return ((electrons[0].electron_pt > 20e3 && electrons[1].electron_pt > 27e3) || (electrons[1].electron_pt > 20e3 && electrons[0].electron_pt > 27e3));
//        }, {"di_electrons"});
//
//        auto deltaR = leadingPt.Filter([] (const RVec<AbstractParticle>& electrons)
//        {
//            return (DeltaR(electrons[0].ElectronVector(), electrons[1].ElectronVector()) > 0.01);
//        }, {"di_electrons"});
//
//        //new dataframe node: Contains a new column `dilep` in addition to the ones in `same_flavour` that stores the di-electron four-vector
//        auto dilep = deltaR.Define("dilep",[] (RVec<AbstractParticle>& electrons)
//        {
//            return (electrons[0].ElectronVector() + electrons[1].ElectronVector());
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
//        auto photon_passes_cuts = ptCut
//        .Define("photons_pass_cut_indices",
//        [&](const RVec<AbstractParticle>& photons)
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
//        }, {"abstract_photons"})
//        .Define("photons_pass_cuts",
//        [&](RVec<AbstractParticle>& photons, RVec<int>& photon_indices)
//        {
//            return Take(photons, photon_indices);
//        }, {"abstract_photons", "photons_pass_cut_indices"});
//
//        auto merged_reco_photons_matched = photon_passes_cuts.Filter(
//        [&](const RVec<AbstractParticle>& reco_photons_matched)
//        {
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
//                delta_r = DeltaR(reco_photons_matched[combs[0][i]].PhotonVector(), reco_photons_matched[combs[1][i]].PhotonVector());
//                m = (reco_photons_matched[combs[0][i]].PhotonVector() + reco_photons_matched[combs[1][i]].PhotonVector()).M();
//                pt = (reco_photons_matched[combs[0][i]].PhotonVector() + reco_photons_matched[combs[1][i]].PhotonVector()).Pt();
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
//        .Define("merged_photon_index", //new column: consists of the index corresponding to the photon that made the event be classified as merged
//        [&](const RVec<AbstractParticle>& rpm) //rpm = reco photons matched
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
//        [&](const RVec<AbstractParticle>& reco_photons_matched, int merged_photon_index)
//        {
//            return reco_photons_matched[merged_photon_index];
//
//        }, {"photons_pass_cuts", "merged_photon_index"});
//
//        auto dilepton_and_photon = merged_reco_photons_matched
//        .Define("reconstructed_mass",[&](const RVec<AbstractParticle>& di_electrons, const AbstractParticle& merged_photon)
//        {
//            auto four_momentum = di_electrons[0].ElectronVector() + di_electrons[1].ElectronVector();
//
//            return (four_momentum + merged_photon.PhotonVector()).M()/1e3;
//
//        }, {"di_electrons", "merged_photon"});
//
//        auto pSB = dilepton_and_photon.Filter(
//        [](const double reconstructed_mass)
//        {
//            return (reconstructed_mass < 110) || (reconstructed_mass > 130);
//        }, {"reconstructed_mass"});
//
//        auto pSR = dilepton_and_photon.Filter(
//        [](const double reconstructed_mass)
//        {
//            return (reconstructed_mass >= 110) && (reconstructed_mass <= 130);
//        }, {"reconstructed_mass"});
//
//        auto SB = pSB.Filter(
//        [](const RVec<float>& Eratio)
//        {
//            return (!Any(Eratio <= 0.8));
//        }, {"photon_shower_shape_e_ratio"});
//
//        auto SR = pSR.Filter(
//        [](const RVec<float>& Eratio)
//        {
//            return (!Any(Eratio <= 0.8));
//        }, {"photon_shower_shape_e_ratio"});
//
//        auto totEventWeight = merged_reco_photons_matched
//        .Define("totEventWeight", [](RVec<float> totEventWeightVec, RVec<int>& photons_pass_cut_indices, int mpi)
//        {
//            return Take(totEventWeightVec, photons_pass_cut_indices)[mpi];
//
//        }, {"totEventWeightVec", "photons_pass_cut_indices", "merged_photon_index"});
//
//        Totals.push_back(EventWeight.Sum<float>("EventWeight"));
//
//        resultmaps.push_back(VariationsFor(totEventWeight.Sum<float>("totEventWeight")));
//    }
//
////EG_RESOLUTION_ALL__1down:           42864.7 ID 1.00524  ISO 1.02283  TRIG 1.04157
////EG_RESOLUTION_ALL__1up:             42751.5 ID 1.00524  ISO 1.02283  TRIG 1.04157
////EG_SCALE_ALL__1down:                42521.1 ID 1.00524  ISO 1.02283  TRIG 1.04157
////EG_SCALE_ALL__1up:                  43059.8 ID 1.00524  ISO 1.02283  TRIG 1.04157
////PH_EFF_ID_Uncertainty__1down:       42790.4 ID 0.992328 ISO 1.02283  TRIG 1.04157
////PH_EFF_ID_Uncertainty__1up:         42790.4 ID 1.01816  ISO 1.02283  TRIG 1.04157
////PH_EFF_ISO_Uncertainty__1down:      42790.4 ID 1.00524  ISO 1.00971  TRIG 1.04157
////PH_EFF_ISO_Uncertainty__1up:        42790.4 ID 1.00524  ISO 1.03595  TRIG 1.04157
////PH_EFF_TRIGGER_Uncertainty__1down:  42790.4 ID 1.00524  ISO 1.02283  TRIG 1.01646
////PH_EFF_TRIGGER_Uncertainty__1up:    42790.4 ID 1.00524  ISO 1.02283  TRIG 1.06668
//
//    ROOT::RDF::RunGraphs(Totals); // running all computation nodes concurrently
//
////    for (auto& i: GeneratorWeightCounts)
////    {
////        std::cout << *i << '\n';
////    }
////    int count = 0;
////    for (auto& i: resultmaps)
////    {
////        std::cout << count++ << "\n==\n";
////        for (auto& j: i.GetKeys())
////        {
////            std::cout << j << '\n';
////        }
////        std::cout << '\n';
////    }
////          resultmaps
////          ----------
////    0       1       //Z-gamma
////    2       3       //Z-gamma
////    4       5       //Z-gamma
////    6       7       //data
////    8       9       //ma1
////    10      11      //ma2
////    12      13      //ma3
////    14      15      //ma5
////    16      17      //ma9
////    18      19      //Z-jets
////    20      21      //Z-jets
////    22      23      //Z-jets
////    24      25      //Z-jets
////    26      27      //Z-jets
////    28      29      //Z-jets
////    30      31      //Z-jets
////    32      33      //Z-jets
////    34      35      //Z-jets
//
//    ss << R"--(\section*{Table 22})--" << '\n';
//    ss << R"--(\hspace{-3cm}\scalebox{0.8}{)--" << '\n';
//    ss << R"--(\setlength\extrarowheight{2pt}\renewcommand{\arraystretch}{1.5})--" << '\n';
////    ss << R"--(\begin{tabular}{|c|c|c|c|c|c|})--" << '\n';
//    ss << R"--(\begin{tabular}{|c|c|c|c|})--" << '\n';
//    ss << R"--(\hline)--" << '\n';
//
//        ss << R"--(\multicolumn{4}{|c|}{\parbox{\linewidth}{\centering Merged Photon Category: Down Variations \\ (\% difference from nominal)}}\\[5 pt])--" << '\n';
//    ss << R"--(\hline)--" << '\n';
//
//    ss << R"--({Sample} & EG\_RESOLUTION\_ALL & EG\_SCALE\_ALL & PH\_EFF\_ID\_Uncertainty \\ \hline)--" << '\n';
//
//    double finalScaleVal;
//
//    double ZgammaNominal = 0, ZgammaEG_RESOLUTION_ALL = 0,
//    ZgammaEG_SCALE_ALL = 0, ZgammaPH_EFF_ISO_Uncertainty = 0,
//    ZgammaPH_EFF_ID_Uncertainty = 0, ZgammaPH_EFF_TRIGGER_Uncertainty = 0;
//
//    double ZjetsNominal = 0, ZjetsEG_RESOLUTION_ALL = 0,
//    ZjetsEG_SCALE_ALL = 0, ZjetsPH_EFF_ISO_Uncertainty = 0,
//    ZjetsPH_EFF_ID_Uncertainty = 0, ZjetsPH_EFF_TRIGGER_Uncertainty = 0;
//
//    double totbkgNominal = 0, totbkgEG_RESOLUTION_ALL = 0,
//    totbkgEG_SCALE_ALL = 0, totbkgPH_EFF_ISO_Uncertainty = 0,
//    totbkgPH_EFF_ID_Uncertainty = 0, totbkgPH_EFF_TRIGGER_Uncertainty = 0;
//
//    for (auto i = 0; (i < resultmaps.size()); i++)
//    {
//        double nominalVal = static_cast<double>(resultmaps[i]["nominal"]);
//
////        auto denominator = *Totals[i].GetResultPtr<ULong64_t>();
//        auto denominator = *Totals[i].GetResultPtr<float>();
//
//        if (i >= 0 && i <= 2) //Zgamma
//        {
//            finalScaleVal = SFs[i]/denominator;
//            ZgammaNominal += finalScaleVal*nominalVal;
//
//            ZgammaEG_RESOLUTION_ALL += finalScaleVal*resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1down"];
//            ZgammaEG_SCALE_ALL += finalScaleVal*resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1down"];
////            ZgammaPH_EFF_ISO_Uncertainty += finalScaleVal*resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1down"];
//            ZgammaPH_EFF_ID_Uncertainty += finalScaleVal*resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1down"];
////            ZgammaPH_EFF_TRIGGER_Uncertainty += finalScaleVal*resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down"];
//        }
//
//        else if (i >= 9)
//        {
//            finalScaleVal = JetNumeratorSFs[i-9]/denominator;
//            ZjetsNominal += finalScaleVal*nominalVal;
//
//            ZjetsEG_RESOLUTION_ALL += finalScaleVal*resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1down"];
//            ZjetsEG_SCALE_ALL += finalScaleVal*resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1down"];
////            ZjetsPH_EFF_ISO_Uncertainty += finalScaleVal*resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1down"];
//            ZjetsPH_EFF_ID_Uncertainty += finalScaleVal*resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1down"];
////            ZjetsPH_EFF_TRIGGER_Uncertainty += finalScaleVal*resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down"];
//        }
//
//        ss << prefixes[i] << " & ";
//        ss << std::setprecision(4) << std::fixed
//        << ( (nominalVal && resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1down"]) ? ((resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1down"]-nominalVal)/nominalVal)*100.0 : 0.0)
//        << " & " << std::setprecision(4) << std::fixed
//        << ( (nominalVal && resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1down"]) ? ((resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1down"]-nominalVal)/nominalVal)*100.0 : 0.0)
////        << " & " << std::setprecision(4) << std::fixed
////        << ( (nominalVal && resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1down"]) ? ((resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1down"]-nominalVal)/nominalVal)*100.0 : 0.0)
//        << " & " << std::setprecision(4) << std::fixed
//        << ( (nominalVal && resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1down"]) ? ((resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1down"]-nominalVal)/nominalVal)*100.0 : 0.0)
////        << " & " << std::setprecision(4) << std::fixed
////        << ( (nominalVal && resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down"]) ? ((resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down"]-nominalVal)/nominalVal)*100.0 : 0.0)
//        << R"--( \\ \hline)--" << '\n';
//
////        if (!resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1down"])
////        {
////            std::cout << i << ": photon_iso_eff:PH_EFF_ISO_Uncertainty__1down not found\n";
////        }
////        else
////        {
////            std::cout << i << ':' << ' ' << resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1down"] << ' ' << nominalVal << '\n';
////        }
//        if (!resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1down"])
//        {
//            std::cout << i << ": photon_id_eff:PH_EFF_ID_Uncertainty__1down not found\n";
//        }
//        else
//        {
//            std::cout << i << ':' << ' ' << resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1down"] << ' ' << nominalVal << '\n';
//        }
////        if (!resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down"])
////        {
////            std::cout << i << ": photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down not found\n";
////        }
////        else
////        {
////            std::cout << i << ':' << ' ' << resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down"] << ' ' << nominalVal << '\n';
////        }
//    }
//
//    ss << R"--(Total $Z\gamma\gamma$ & )--";
//    ss << std::setprecision(4) << std::fixed
//    << ( (ZgammaNominal && ZgammaEG_RESOLUTION_ALL) ? ((ZgammaEG_RESOLUTION_ALL-ZgammaNominal)/ZgammaNominal)*100.0 : 0.0)
//    << " & " << std::setprecision(4) << std::fixed
//    << ( (ZgammaNominal && ZgammaEG_SCALE_ALL) ? ((ZgammaEG_SCALE_ALL-ZgammaNominal)/ZgammaNominal)*100.0 : 0.0)
////    << " & " << std::setprecision(4) << std::fixed
////    << ( (ZgammaNominal && ZgammaPH_EFF_ISO_Uncertainty) ? ((ZgammaPH_EFF_ISO_Uncertainty-ZgammaNominal)/ZgammaNominal)*100.0 : 0.0)
//    << " & " << std::setprecision(4) << std::fixed
//    << ( (ZgammaNominal && ZgammaPH_EFF_ID_Uncertainty) ? ((ZgammaPH_EFF_ID_Uncertainty-ZgammaNominal)/ZgammaNominal)*100.0 : 0.0)
////    << " & " << std::setprecision(4) << std::fixed
////    << ( (ZgammaNominal && ZgammaPH_EFF_TRIGGER_Uncertainty) ? ((ZgammaPH_EFF_TRIGGER_Uncertainty-ZgammaNominal)/ZgammaNominal)*100.0 : 0.0)
//    << R"--( \\ \hline)--" << '\n';
//
//    ss << R"--(Total $Z$+jets & )--";
//    ss << std::setprecision(4) << std::fixed
//    << ( (ZjetsNominal && ZjetsEG_RESOLUTION_ALL) ? ((ZjetsEG_RESOLUTION_ALL-ZjetsNominal)/ZjetsNominal)*100.0 : 0.0)
//    << " & " << std::setprecision(4) << std::fixed
//    << ( (ZjetsNominal && ZjetsEG_SCALE_ALL) ? ((ZjetsEG_SCALE_ALL-ZjetsNominal)/ZjetsNominal)*100.0 : 0.0)
////    << " & " << std::setprecision(4) << std::fixed
////    << ( (ZjetsNominal && ZjetsPH_EFF_ISO_Uncertainty) ? ((ZjetsPH_EFF_ISO_Uncertainty-ZjetsNominal)/ZjetsNominal)*100.0 : 0.0)
//    << " & " << std::setprecision(4) << std::fixed
//    << ( (ZjetsNominal && ZjetsPH_EFF_ID_Uncertainty) ? ((ZjetsPH_EFF_ID_Uncertainty-ZjetsNominal)/ZjetsNominal)*100.0 : 0.0)
////    << " & " << std::setprecision(4) << std::fixed
////    << ( (ZjetsNominal && ZjetsPH_EFF_TRIGGER_Uncertainty) ? ((ZjetsPH_EFF_TRIGGER_Uncertainty-ZjetsNominal)/ZjetsNominal)*100.0 : 0.0)
//    << R"--( \\ \hline)--" << '\n';
//
//    totbkgNominal = ZgammaNominal + ZjetsNominal;
//
//    totbkgEG_RESOLUTION_ALL = ZgammaEG_RESOLUTION_ALL + ZjetsEG_RESOLUTION_ALL;
//    totbkgEG_SCALE_ALL = ZgammaEG_SCALE_ALL + ZjetsEG_SCALE_ALL;
//    totbkgPH_EFF_ISO_Uncertainty = ZgammaPH_EFF_ISO_Uncertainty + ZjetsPH_EFF_ISO_Uncertainty;
//    totbkgPH_EFF_ID_Uncertainty = ZgammaPH_EFF_ID_Uncertainty + ZjetsPH_EFF_ID_Uncertainty;
//    totbkgPH_EFF_TRIGGER_Uncertainty = ZgammaPH_EFF_TRIGGER_Uncertainty + ZjetsPH_EFF_TRIGGER_Uncertainty;
//
//    ss << R"--(Total Bkg & )--";
//    ss << std::setprecision(4) << std::fixed
//    << ( (totbkgNominal && totbkgEG_RESOLUTION_ALL) ? ((totbkgEG_RESOLUTION_ALL-totbkgNominal)/totbkgNominal)*100.0 : 0.0)
//    << " & " << std::setprecision(4) << std::fixed
//    << ( (totbkgNominal && totbkgEG_SCALE_ALL) ? ((totbkgEG_SCALE_ALL-totbkgNominal)/totbkgNominal)*100.0 : 0.0)
////    << " & " << std::setprecision(4) << std::fixed
////    << ( (totbkgNominal && totbkgPH_EFF_ISO_Uncertainty) ? ((totbkgPH_EFF_ISO_Uncertainty-totbkgNominal)/totbkgNominal)*100.0 : 0.0)
//    << " & " << std::setprecision(4) << std::fixed
//    << ( (totbkgNominal && totbkgPH_EFF_ID_Uncertainty) ? ((totbkgPH_EFF_ID_Uncertainty-totbkgNominal)/totbkgNominal)*100.0 : 0.0)
////    << " & " << std::setprecision(4) << std::fixed
////    << ( (totbkgNominal && totbkgPH_EFF_TRIGGER_Uncertainty) ? ((totbkgPH_EFF_TRIGGER_Uncertainty-totbkgNominal)/totbkgNominal)*100.0 : 0.0)
//    << R"--( \\ \hline)--" << '\n';
//
//    ss << R"--(\end{tabular}})--" << '\n';
//
//    ss << "\n\n\n";
//
//    std::ofstream out("Table22.txt");
//    out << ss.str() << '\n';
//    out.close();
//}
//
//void Table23()
//{
//    Event::systematics =
//    {
////        "PH_EFF_ISO_Uncertainty",
////        "PH_EFF_ISO_Uncertainty",
////        "EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR",
////        "PRW_DATASF",
////        "MUON_EFF_RECO_SYS",
////        "MUON_EFF_ISO_SYS",
////        "MUON_EFF_TrigSystUncertainty",
////        "EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR",
////        "MUON_EFF_TrigStatUncertainty",
////        "MUON_EFF_RECO_STAT",
////        "MUON_EFF_TTVA_STAT",
////        "EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR",
////        "EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR",
////        "MUON_EFF_TTVA_SYS",
////        "MUON_EFF_ISO_STAT",
////        "MUON_SAGITTA_RHO",
////        "EG_RESOLUTION_ALL",
////        "EG_SCALE_ALL",
////        "MUON_MS",
////        "MUON_ID",
////        "EL_EFF_TriggerEff_TOTAL_1NPCOR_PLUS_UNCOR",
////        "MUON_SAGITTA_RESBIAS",
////        "MUON_SCALE",
//
////        "PH_EFF_ISO_Uncertainty__1down",
////        "EG_SCALE_ALL__1down",
////        "PH_EFF_ID_Uncertainty__1down",
////        "PH_EFF_TRIGGER_Uncertainty__1down",
//        "EG_RESOLUTION_ALL__1up",
//        "EG_SCALE_ALL__1up",
////        "PH_EFF_ISO_Uncertainty__1up",
//        "PH_EFF_ID_Uncertainty__1up",
////        "PH_EFF_TRIGGER_Uncertainty__1up",
//    };
//
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
//
//    std::array<double,9> JetNumeratorSFs = {((139e15)*(1.9828e-9)*(0.821204)),((139e15)*(110.64e-12)*(0.69275)),((139e15)*(40.645e-12)*(0.615906)),((139e15)*(1.9817e-9)*(0.1136684)),((139e15)*(110.47e-12)*(0.1912956)),((139e15)*(40.674e-12)*(0.2326772)),((139e15)*(1.9819e-9)*(0.0656969)),((139e15)*(110.53e-12)*(0.1158741)),((139e15)*(40.68e-12)*(0.1535215))};
//    std::array<double,3> SFs = {((139e15)*(.871e-12)),((139e15)*(.199e-12)), ((139e15)*(.0345e-15))}; //numerators for Z-gamma bkg
//
//    std::vector<std::string> prefixes = { R"--(pty2\_9\_17)--", R"--(pty\_17\_myy\_0\_80)--", R"--(pty\_17\_myy\_80)--", "data", R"--($\text{Sig } m_{A}$ = 1 GeV)--", R"--($\text{Sig } m_{A}$ = 2 GeV)--", R"--($\text{Sig } m_{A}$ = 3 GeV)--", R"--($\text{Sig } m_{A}$ = 5 GeV)--", R"--($\text{Sig } m_{A}$ = 9 GeV)--", R"--(Zee\_lightJet\_0-70)--", R"--(Zee\_lightJet\_70-140)--", R"--(Zee\_lightJet\_140-280)--", R"--(Zee\_cJet\_0-70)--", R"--(Zee\_cJet\_70-140)--", R"--(Zee\_cJet\_140-280)--", R"--(Zee\_bJet\_0-70)--", R"--(Zee\_bJet\_70-140)--", R"--(Zee\_bJet\_140-280)--"};
//
//    std::vector<ROOT::RDF::RResultHandle> Totals;
//    std::vector<RResultMap<float>> resultmaps;
//
//    std::stringstream ss;
//
//    for (auto& i: input_filenames)
//    {
//        SchottDataFrame df(MakeRDF(i,8));
//
//        auto EventWeight = df.Define("EventWeight", //mc generator weight
//        [](const RVec<float>& ei_event_weights_generator)
//        {
//            return ((ei_event_weights_generator[0]) ? 1 / ei_event_weights_generator[0] : 1);
//
//        }, {"ei_event_weights_generator"})
//        .Define("totEventWeightVec", //vector of efficiencies, product of photon_id_eff, photon_iso_eff, and photon_trg_eff
//        [](RVec<float> photon_id_eff, RVec<float> photon_iso_eff, RVec<float> photon_trg_eff/*, RVec<float> ei_event_weights_generator*/)
//        {
//            auto ResizeVal = std::max({photon_id_eff.size(), photon_iso_eff.size(), photon_trg_eff.size()});
//            photon_id_eff.resize(ResizeVal,1);
//            photon_iso_eff.resize(ResizeVal,1);
//            photon_trg_eff.resize(ResizeVal,1);
//
//            return photon_id_eff;
//            //*photon_iso_eff*photon_trg_eff; //Returns a vector of efficiencies
//
//        }, {"photon_id_eff", "photon_iso_eff", "photon_trg_eff",/* "ei_event_weights_generator"*/});
//
//        auto two_leptons = EventWeight
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
//        [](RVec<AbstractParticle> electrons)
//        {
//            //keep the electrons in each event that have pt > 20 GeV, |η| < 2.37,
//            //|η| not between 1.37 and 1.52, and that satisfy a medium id criteria `electron_id_medium`
//            electrons.erase(std::remove_if(electrons.begin(),electrons.end(),
//            [](AbstractParticle& ep)
//            {
//                 return (!((ep.electron_pt/1e3 > 20) && (std::abs(ep.electron_eta) < 2.37) &&
//                 (!((1.37 < std::abs(ep.electron_eta)) && (std::abs(ep.electron_eta) < 1.52)))
//                 && (ep.electron_id_medium == 1)));
//
//            }), electrons.end());
//
//            return electrons;
//
//        },{"abstract_electrons"}).Filter(
//        [](RVec<Muon>& muons, RVec<AbstractParticle>& electrons)
//        {
//            return (electrons.size()==2 && muons.empty()); //keep events which have exactly 2 electrons for di_electrons and no muons
//
//        }, {"muons", "di_electrons"});
//
//        //new dataframe node: contains only the events from `two_leptons` whose electrons in the `di_electrons` branch have opposite charge
//        auto opp_charge = two_leptons
//        .Filter([](const RVec<AbstractParticle>& electrons)
//        {
//            return (electrons[0].electron_charge*electrons[1].electron_charge < 0);
//
//        }, {"di_electrons"});
//
//        //new dataframe node: contains only the events from `opp_charge` that have 1 electron with pt > 20 GeV and the other with pt > 27 GeV
//        auto leadingPt = opp_charge.Filter([](const RVec<AbstractParticle>& electrons)
//        {
//            return ((electrons[0].electron_pt > 20e3 && electrons[1].electron_pt > 27e3) || (electrons[1].electron_pt > 20e3 && electrons[0].electron_pt > 27e3));
//        }, {"di_electrons"});
//
//        auto deltaR = leadingPt.Filter([] (const RVec<AbstractParticle>& electrons)
//        {
//            return (DeltaR(electrons[0].ElectronVector(), electrons[1].ElectronVector()) > 0.01);
//        }, {"di_electrons"});
//
//        //new dataframe node: Contains a new column `dilep` in addition to the ones in `same_flavour` that stores the di-electron four-vector
//        auto dilep = deltaR.Define("dilep",[] (RVec<AbstractParticle>& electrons)
//        {
//            return (electrons[0].ElectronVector() + electrons[1].ElectronVector());
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
//        auto photon_passes_cuts = ptCut
//        .Define("photons_pass_cut_indices",
//        [&](const RVec<AbstractParticle>& photons)
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
//        }, {"abstract_photons"})
//        .Define("photons_pass_cuts",
//        [&](RVec<AbstractParticle>& photons, RVec<int>& photon_indices)
//        {
//            return Take(photons, photon_indices);
//        }, {"abstract_photons", "photons_pass_cut_indices"});
//
//        //New dataframe node: contains only events from `photon_passes_cuts` that pass the resolved category
//        auto resolved_reco_photons_matched = photon_passes_cuts.Define("chosen_two_indices",
//        [](RVec<AbstractParticle>& photons_pass_cuts)
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
//                delta_r = DeltaR(photons_pass_cuts[combs[0][i]].PhotonVector(), photons_pass_cuts[combs[1][i]].PhotonVector());
//                m = (photons_pass_cuts[combs[0][i]].PhotonVector() + photons_pass_cuts[combs[1][i]].PhotonVector()).M();
//                pt = (photons_pass_cuts[combs[0][i]].PhotonVector() + photons_pass_cuts[combs[1][i]].PhotonVector()).Pt();
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
//        [&](RVec<AbstractParticle>& reco_photons_matched, RVec<unsigned long>& indices)
//        {
//            return Take(reco_photons_matched, indices);
//
//        }, {"photons_pass_cuts", "chosen_two_indices"})
//        .Define("reconstructed_mass",
//        [&](RVec<AbstractParticle>& diph)
//        {
//            return (diph[0].PhotonVector()+diph[1].PhotonVector()).M()/1e3;
//        }, {"chosen_two"});
//
//        auto totEventWeight = resolved_reco_photons_matched.Define("totEventWeight", //New column: weight factor for events in RDF `resolved`
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
//        auto SB = resolved_reco_photons_matched.Filter(
//        [](const double reconstructed_mass)
//        {
//            return (reconstructed_mass < 110) || (reconstructed_mass > 140);
//        }, {"reconstructed_mass"});
//
//        auto SR = resolved_reco_photons_matched.Filter(
//        [](const double reconstructed_mass)
//        {
//            return (reconstructed_mass >= 110) && (reconstructed_mass <= 140);
//        }, {"reconstructed_mass"});
//
//
////        Totals.push_back(df.Count());
////        Totals.push_back(EventWeight.Sum<RVec<float>>("EventWeight"));
//        Totals.push_back(EventWeight.Sum<float>("EventWeight"));
////        resultmaps.push_back(VariationsFor(merged_reco_photons_matched.Count()));
//        resultmaps.push_back(VariationsFor(totEventWeight.Sum<float>("totEventWeight")));
//
//    }
//
////EG_RESOLUTION_ALL__1down:           42864.7 ID 1.00524  ISO 1.02283  TRIG 1.04157
////EG_RESOLUTION_ALL__1up:             42751.5 ID 1.00524  ISO 1.02283  TRIG 1.04157
////EG_SCALE_ALL__1down:                42521.1 ID 1.00524  ISO 1.02283  TRIG 1.04157
////EG_SCALE_ALL__1up:                  43059.8 ID 1.00524  ISO 1.02283  TRIG 1.04157
////PH_EFF_ID_Uncertainty__1down:       42790.4 ID 0.992328 ISO 1.02283  TRIG 1.04157
////PH_EFF_ID_Uncertainty__1up:         42790.4 ID 1.01816  ISO 1.02283  TRIG 1.04157
////PH_EFF_ISO_Uncertainty__1down:      42790.4 ID 1.00524  ISO 1.00971  TRIG 1.04157
////PH_EFF_ISO_Uncertainty__1up:        42790.4 ID 1.00524  ISO 1.03595  TRIG 1.04157
////PH_EFF_TRIGGER_Uncertainty__1down:  42790.4 ID 1.00524  ISO 1.02283  TRIG 1.01646
////PH_EFF_TRIGGER_Uncertainty__1up:    42790.4 ID 1.00524  ISO 1.02283  TRIG 1.06668
//
//    ROOT::RDF::RunGraphs(Totals); // running all computation nodes concurrently
//
////    for (auto& i: GeneratorWeightCounts)
////    {
////        std::cout << *i << '\n';
////    }
////    int count = 0;
////    for (auto& i: resultmaps)
////    {
////        std::cout << count++ << "\n==\n";
////        for (auto& j: i.GetKeys())
////        {
////            std::cout << j << '\n';
////        }
////        std::cout << '\n';
////    }
////          resultmaps
////          ----------
////    0       1       //Z-gamma
////    2       3       //Z-gamma
////    4       5       //Z-gamma
////    6       7       //data
////    8       9       //ma1
////    10      11      //ma2
////    12      13      //ma3
////    14      15      //ma5
////    16      17      //ma9
////    18      19      //Z-jets
////    20      21      //Z-jets
////    22      23      //Z-jets
////    24      25      //Z-jets
////    26      27      //Z-jets
////    28      29      //Z-jets
////    30      31      //Z-jets
////    32      33      //Z-jets
////    34      35      //Z-jets
//
//    ss << R"--(\section*{Table 23})--" << '\n';
//    ss << R"--(\hspace{-3cm}\scalebox{0.8}{)--" << '\n';
//    ss << R"--(\setlength\extrarowheight{2pt}\renewcommand{\arraystretch}{1.5})--" << '\n';
////    ss << R"--(\begin{tabular}{|c|c|c|c|c|c|})--" << '\n';
//    ss << R"--(\begin{tabular}{|c|c|c|c|})--" << '\n';
//    ss << R"--(\hline)--" << '\n';
//
//        ss << R"--(\multicolumn{4}{|c|}{\parbox{\linewidth}{\centering Resolved Photon Category: Up Variations \\ (\% difference from nominal)}}\\[5 pt])--" << '\n';
//    ss << R"--(\hline)--" << '\n';
//
////    ss << R"--({Sample} & EG\_RESOLUTION\_ALL & EG\_SCALE\_ALL & PH\_EFF\_ISO\_Uncertainty & PH\_EFF\_ID\_Uncertainty & PH\_EFF\_TRIGGER\_Uncertainty \\ \hline)--" << '\n';
//    ss << R"--({Sample} & EG\_RESOLUTION\_ALL & EG\_SCALE\_ALL & PH\_EFF\_ID\_Uncertainty \\ \hline)--" << '\n';
//
//    double finalScaleVal;
//
//    double ZgammaNominal = 0, ZgammaEG_RESOLUTION_ALL = 0,
//    ZgammaEG_SCALE_ALL = 0, ZgammaPH_EFF_ISO_Uncertainty = 0,
//    ZgammaPH_EFF_ID_Uncertainty = 0, ZgammaPH_EFF_TRIGGER_Uncertainty = 0;
//
//    double ZjetsNominal = 0, ZjetsEG_RESOLUTION_ALL = 0,
//    ZjetsEG_SCALE_ALL = 0, ZjetsPH_EFF_ISO_Uncertainty = 0,
//    ZjetsPH_EFF_ID_Uncertainty = 0, ZjetsPH_EFF_TRIGGER_Uncertainty = 0;
//
//    double totbkgNominal = 0, totbkgEG_RESOLUTION_ALL = 0,
//    totbkgEG_SCALE_ALL = 0, totbkgPH_EFF_ISO_Uncertainty = 0,
//    totbkgPH_EFF_ID_Uncertainty = 0, totbkgPH_EFF_TRIGGER_Uncertainty = 0;
//
//    for (auto i = 0; (i < resultmaps.size()); i++)
//    {
//        double nominalVal = static_cast<double>(resultmaps[i]["nominal"]);
////        auto denominator = *Totals[i].GetResultPtr<ULong64_t>();
//        auto denominator = *Totals[i].GetResultPtr<float>();
//
//        if (i >= 0 && i <= 2) //Zgamma
//        {
//            finalScaleVal = SFs[i]/denominator;
//            ZgammaNominal += finalScaleVal*nominalVal;
//
//            ZgammaEG_RESOLUTION_ALL += finalScaleVal*resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1up"];
//            ZgammaEG_SCALE_ALL += finalScaleVal*resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1up"];
////            ZgammaPH_EFF_ISO_Uncertainty += finalScaleVal*resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1up"];
//            ZgammaPH_EFF_ID_Uncertainty += finalScaleVal*resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1up"];
////            ZgammaPH_EFF_TRIGGER_Uncertainty += finalScaleVal*resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up"];
//        }
//
//        else if (i >= 9)
//        {
//            finalScaleVal = JetNumeratorSFs[i-9]/denominator;
//            ZjetsNominal += finalScaleVal*nominalVal;
//
//            ZjetsEG_RESOLUTION_ALL += finalScaleVal*resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1up"];
//            ZjetsEG_SCALE_ALL += finalScaleVal*resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1up"];
////            ZjetsPH_EFF_ISO_Uncertainty += finalScaleVal*resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1up"];
//            ZjetsPH_EFF_ID_Uncertainty += finalScaleVal*resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1up"];
////            ZjetsPH_EFF_TRIGGER_Uncertainty += finalScaleVal*resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up"];
//        }
//
//        ss << prefixes[i] << " & ";
//        ss << std::setprecision(4) << std::fixed
//        << ( (nominalVal && resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1up"]) ? ((resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1up"]-nominalVal)/nominalVal)*100.0 : 0.0)
//        << " & " << std::setprecision(4) << std::fixed
//        << ( (nominalVal && resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1up"]) ? ((resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1up"]-nominalVal)/nominalVal)*100.0 : 0.0)
////        << " & " << std::setprecision(4) << std::fixed
////        << ( (nominalVal && resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1up"]) ? ((resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1up"]-nominalVal)/nominalVal)*100.0 : 0.0)
//        << " & " << std::setprecision(4) << std::fixed
//        << ( (nominalVal && resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1up"]) ? ((resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1up"]-nominalVal)/nominalVal)*100.0 : 0.0)
////        << " & " << std::setprecision(4) << std::fixed
////        << ( (nominalVal && resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up"]) ? ((resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up"]-nominalVal)/nominalVal)*100.0 : 0.0)
//        << R"--( \\ \hline)--" << '\n';
//
////        if (!resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1up"])
////        {
////            std::cout << i << ": photon_iso_eff:PH_EFF_ISO_Uncertainty__1up not found\n";
////        }
////        else
////        {
////            std::cout << i << ':' << ' ' << resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1up"] << ' ' << nominalVal << '\n';
////        }
//        if (!resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1up"])
//        {
//            std::cout << i << ": photon_id_eff:PH_EFF_ID_Uncertainty__1up not found\n";
//        }
//        else
//        {
//            std::cout << i << ':' << ' ' << resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1up"] << ' ' << nominalVal << '\n';
//        }
////        if (!resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up"])
////        {
////            std::cout << i << ": photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up not found\n";
////        }
////        else
////        {
////            std::cout << i << ':' << ' ' << resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up"] << ' ' << nominalVal << '\n';
////        }
//    }
//
//    ss << R"--(Total $Z\gamma\gamma$ & )--";
//    ss << std::setprecision(4) << std::fixed
//    << ( (ZgammaNominal && ZgammaEG_RESOLUTION_ALL) ? ((ZgammaEG_RESOLUTION_ALL-ZgammaNominal)/ZgammaNominal)*100.0 : 0.0)
//    << " & " << std::setprecision(4) << std::fixed
//    << ( (ZgammaNominal && ZgammaEG_SCALE_ALL) ? ((ZgammaEG_SCALE_ALL-ZgammaNominal)/ZgammaNominal)*100.0 : 0.0)
////    << " & " << std::setprecision(4) << std::fixed
////    << ( (ZgammaNominal && ZgammaPH_EFF_ISO_Uncertainty) ? ((ZgammaPH_EFF_ISO_Uncertainty-ZgammaNominal)/ZgammaNominal)*100.0 : 0.0)
//    << " & " << std::setprecision(4) << std::fixed
//    << ( (ZgammaNominal && ZgammaPH_EFF_ID_Uncertainty) ? ((ZgammaPH_EFF_ID_Uncertainty-ZgammaNominal)/ZgammaNominal)*100.0 : 0.0)
////    << " & " << std::setprecision(4) << std::fixed
////    << ( (ZgammaNominal && ZgammaPH_EFF_TRIGGER_Uncertainty) ? ((ZgammaPH_EFF_TRIGGER_Uncertainty-ZgammaNominal)/ZgammaNominal)*100.0 : 0.0)
//    << R"--( \\ \hline)--" << '\n';
//
//    ss << R"--(Total $Z$+jets & )--";
//    ss << std::setprecision(4) << std::fixed
//    << ( (ZjetsNominal && ZjetsEG_RESOLUTION_ALL) ? ((ZjetsEG_RESOLUTION_ALL-ZjetsNominal)/ZjetsNominal)*100.0 : 0.0)
//    << " & " << std::setprecision(4) << std::fixed
//    << ( (ZjetsNominal && ZjetsEG_SCALE_ALL) ? ((ZjetsEG_SCALE_ALL-ZjetsNominal)/ZjetsNominal)*100.0 : 0.0)
////    << " & " << std::setprecision(4) << std::fixed
////    << ( (ZjetsNominal && ZjetsPH_EFF_ISO_Uncertainty) ? ((ZjetsPH_EFF_ISO_Uncertainty-ZjetsNominal)/ZjetsNominal)*100.0 : 0.0)
//    << " & " << std::setprecision(4) << std::fixed
//    << ( (ZjetsNominal && ZjetsPH_EFF_ID_Uncertainty) ? ((ZjetsPH_EFF_ID_Uncertainty-ZjetsNominal)/ZjetsNominal)*100.0 : 0.0)
////    << " & " << std::setprecision(4) << std::fixed
////    << ( (ZjetsNominal && ZjetsPH_EFF_TRIGGER_Uncertainty) ? ((ZjetsPH_EFF_TRIGGER_Uncertainty-ZjetsNominal)/ZjetsNominal)*100.0 : 0.0)
//    << R"--( \\ \hline)--" << '\n';
//
//    totbkgNominal = ZgammaNominal + ZjetsNominal;
//
//    totbkgEG_RESOLUTION_ALL = ZgammaEG_RESOLUTION_ALL + ZjetsEG_RESOLUTION_ALL;
//    totbkgEG_SCALE_ALL = ZgammaEG_SCALE_ALL + ZjetsEG_SCALE_ALL;
//    totbkgPH_EFF_ISO_Uncertainty = ZgammaPH_EFF_ISO_Uncertainty + ZjetsPH_EFF_ISO_Uncertainty;
//    totbkgPH_EFF_ID_Uncertainty = ZgammaPH_EFF_ID_Uncertainty + ZjetsPH_EFF_ID_Uncertainty;
//    totbkgPH_EFF_TRIGGER_Uncertainty = ZgammaPH_EFF_TRIGGER_Uncertainty + ZjetsPH_EFF_TRIGGER_Uncertainty;
//
//    ss << R"--(Total Bkg & )--";
//    ss << std::setprecision(4) << std::fixed
//    << ( (totbkgNominal && totbkgEG_RESOLUTION_ALL) ? ((totbkgEG_RESOLUTION_ALL-totbkgNominal)/totbkgNominal)*100.0 : 0.0)
//    << " & " << std::setprecision(4) << std::fixed
//    << ( (totbkgNominal && totbkgEG_SCALE_ALL) ? ((totbkgEG_SCALE_ALL-totbkgNominal)/totbkgNominal)*100.0 : 0.0)
////    << " & " << std::setprecision(4) << std::fixed
////    << ( (totbkgNominal && totbkgPH_EFF_ISO_Uncertainty) ? ((totbkgPH_EFF_ISO_Uncertainty-totbkgNominal)/totbkgNominal)*100.0 : 0.0)
//    << " & " << std::setprecision(4) << std::fixed
//    << ( (totbkgNominal && totbkgPH_EFF_ID_Uncertainty) ? ((totbkgPH_EFF_ID_Uncertainty-totbkgNominal)/totbkgNominal)*100.0 : 0.0)
////    << " & " << std::setprecision(4) << std::fixed
////    << ( (totbkgNominal && totbkgPH_EFF_TRIGGER_Uncertainty) ? ((totbkgPH_EFF_TRIGGER_Uncertainty-totbkgNominal)/totbkgNominal)*100.0 : 0.0)
//    << R"--( \\ \hline)--" << '\n';
//
//    ss << R"--(\end{tabular}})--" << '\n';
//
//    ss << "\n\n\n";
//
//    std::ofstream out("Table23.txt");
//    out << ss.str() << '\n';
//    out.close();
//}
//
//void Table24()
//{
//    Event::systematics =
//    {
////        "PH_EFF_ISO_Uncertainty",
////        "PH_EFF_ISO_Uncertainty",
////        "EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR",
////        "PRW_DATASF",
////        "MUON_EFF_RECO_SYS",
////        "MUON_EFF_ISO_SYS",
////        "MUON_EFF_TrigSystUncertainty",
////        "EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR",
////        "MUON_EFF_TrigStatUncertainty",
////        "MUON_EFF_RECO_STAT",
////        "MUON_EFF_TTVA_STAT",
////        "EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR",
////        "EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR",
////        "MUON_EFF_TTVA_SYS",
////        "MUON_EFF_ISO_STAT",
////        "MUON_SAGITTA_RHO",
////        "EG_RESOLUTION_ALL",
////        "EG_SCALE_ALL",
////        "MUON_MS",
////        "MUON_ID",
////        "EL_EFF_TriggerEff_TOTAL_1NPCOR_PLUS_UNCOR",
////        "MUON_SAGITTA_RESBIAS",
////        "MUON_SCALE",
////        "EG_RESOLUTION_ALL__1up",
////        "EG_SCALE_ALL__1up",
////        "PH_EFF_ISO_Uncertainty__1up",
////        "PH_EFF_ID_Uncertainty__1up",
////        "PH_EFF_TRIGGER_Uncertainty__1up",
//        "EG_RESOLUTION_ALL__1down",
//        "EG_SCALE_ALL__1down",
////        "PH_EFF_ISO_Uncertainty__1down",
//        "PH_EFF_ID_Uncertainty__1down",
////        "PH_EFF_TRIGGER_Uncertainty__1down",
//    };
//
//    std::vector<std::vector<std::string>> input_filenames =
//    {
//        //Z-gamma
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617070._000001.LGNTuple.root"},
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617064._000001.LGNTuple.root"},
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
//    std::array<double,9> JetNumeratorSFs = {((139e15)*(1.9828e-9)*(0.821204)),((139e15)*(110.64e-12)*(0.69275)),((139e15)*(40.645e-12)*(0.615906)),((139e15)*(1.9817e-9)*(0.1136684)),((139e15)*(110.47e-12)*(0.1912956)),((139e15)*(40.674e-12)*(0.2326772)),((139e15)*(1.9819e-9)*(0.0656969)),((139e15)*(110.53e-12)*(0.1158741)),((139e15)*(40.68e-12)*(0.1535215))};
//    std::array<double,3> SFs = {((139e15)*(.871e-12)),((139e15)*(.199e-12)), ((139e15)*(.0345e-15))}; //numerators for Z-gamma bkg
//
//    std::vector<std::string> prefixes = { R"--(pty2\_9\_17)--", R"--(pty\_17\_myy\_0\_80)--", R"--(pty\_17\_myy\_80)--", "data", R"--($\text{Sig } m_{A}$ = 1 GeV)--", R"--($\text{Sig } m_{A}$ = 2 GeV)--", R"--($\text{Sig } m_{A}$ = 3 GeV)--", R"--($\text{Sig } m_{A}$ = 5 GeV)--", R"--($\text{Sig } m_{A}$ = 9 GeV)--", R"--(Zee\_lightJet\_0-70)--", R"--(Zee\_lightJet\_70-140)--", R"--(Zee\_lightJet\_140-280)--", R"--(Zee\_cJet\_0-70)--", R"--(Zee\_cJet\_70-140)--", R"--(Zee\_cJet\_140-280)--", R"--(Zee\_bJet\_0-70)--", R"--(Zee\_bJet\_70-140)--", R"--(Zee\_bJet\_140-280)--"};
//
//    std::vector<ROOT::RDF::RResultHandle> Totals;
//    std::vector<RResultMap<float>> resultmaps;
//
//    std::stringstream ss;
//
//    for (auto& i: input_filenames)
//    {
//        SchottDataFrame df(MakeRDF(i,8));
//
//        auto EventWeight = df.Define("EventWeight", //mc generator weight
//        [](const RVec<float>& ei_event_weights_generator)
//        {
//            return ((ei_event_weights_generator[0]) ? 1 / ei_event_weights_generator[0] : 1);
//
//        }, {"ei_event_weights_generator"})
//        .Define("totEventWeightVec", //vector of efficiencies, product of photon_id_eff, photon_iso_eff, and photon_trg_eff
//        [](RVec<float> photon_id_eff, RVec<float> photon_iso_eff, RVec<float> photon_trg_eff/*, RVec<float> ei_event_weights_generator*/)
//        {
//            auto ResizeVal = std::max({photon_id_eff.size(), photon_iso_eff.size(), photon_trg_eff.size()});
//            photon_id_eff.resize(ResizeVal,1);
//            photon_iso_eff.resize(ResizeVal,1);
//            photon_trg_eff.resize(ResizeVal,1);
//
//            return photon_id_eff;
//            // *photon_iso_eff*photon_trg_eff; //Returns a vector of efficiencies
//
//        }, {"photon_id_eff", "photon_iso_eff", "photon_trg_eff",/* "ei_event_weights_generator"*/});
//
//        auto two_leptons = EventWeight
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
//        [](RVec<AbstractParticle> electrons)
//        {
//            //keep the electrons in each event that have pt > 20 GeV, |η| < 2.37,
//            //|η| not between 1.37 and 1.52, and that satisfy a medium id criteria `electron_id_medium`
//            electrons.erase(std::remove_if(electrons.begin(),electrons.end(),
//            [](AbstractParticle& ep)
//            {
//                 return (!((ep.electron_pt/1e3 > 20) && (std::abs(ep.electron_eta) < 2.37) &&
//                 (!((1.37 < std::abs(ep.electron_eta)) && (std::abs(ep.electron_eta) < 1.52)))
//                 && (ep.electron_id_medium == 1)));
//
//            }), electrons.end());
//
//            return electrons;
//
//        },{"abstract_electrons"}).Filter(
//        [](RVec<Muon>& muons, RVec<AbstractParticle>& electrons)
//        {
//            return (electrons.size()==2 && muons.empty()); //keep events which have exactly 2 electrons for di_electrons and no muons
//
//        }, {"muons", "di_electrons"});
//
//        //new dataframe node: contains only the events from `two_leptons` whose electrons in the `di_electrons` branch have opposite charge
//        auto opp_charge = two_leptons
//        .Filter([](const RVec<AbstractParticle>& electrons)
//        {
//            return (electrons[0].electron_charge*electrons[1].electron_charge < 0);
//
//        }, {"di_electrons"});
//
//        //new dataframe node: contains only the events from `opp_charge` that have 1 electron with pt > 20 GeV and the other with pt > 27 GeV
//        auto leadingPt = opp_charge.Filter([](const RVec<AbstractParticle>& electrons)
//        {
//            return ((electrons[0].electron_pt > 20e3 && electrons[1].electron_pt > 27e3) || (electrons[1].electron_pt > 20e3 && electrons[0].electron_pt > 27e3));
//        }, {"di_electrons"});
//
//        auto deltaR = leadingPt.Filter([] (const RVec<AbstractParticle>& electrons)
//        {
//            return (DeltaR(electrons[0].ElectronVector(), electrons[1].ElectronVector()) > 0.01);
//        }, {"di_electrons"});
//
//        //new dataframe node: Contains a new column `dilep` in addition to the ones in `same_flavour` that stores the di-electron four-vector
//        auto dilep = deltaR.Define("dilep",[] (RVec<AbstractParticle>& electrons)
//        {
//            return (electrons[0].ElectronVector() + electrons[1].ElectronVector());
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
//        auto photon_passes_cuts = ptCut
//        .Define("photons_pass_cut_indices",
//        [&](const RVec<AbstractParticle>& photons)
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
//        }, {"abstract_photons"})
//        .Define("photons_pass_cuts",
//        [&](RVec<AbstractParticle>& photons, RVec<int>& photon_indices)
//        {
//            return Take(photons, photon_indices);
//        }, {"abstract_photons", "photons_pass_cut_indices"});
//
//        //New dataframe node: contains only events from `photon_passes_cuts` that pass the resolved category
//        auto resolved_reco_photons_matched = photon_passes_cuts.Define("chosen_two_indices",
//        [](RVec<AbstractParticle>& photons_pass_cuts)
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
//                delta_r = DeltaR(photons_pass_cuts[combs[0][i]].PhotonVector(), photons_pass_cuts[combs[1][i]].PhotonVector());
//                m = (photons_pass_cuts[combs[0][i]].PhotonVector() + photons_pass_cuts[combs[1][i]].PhotonVector()).M();
//                pt = (photons_pass_cuts[combs[0][i]].PhotonVector() + photons_pass_cuts[combs[1][i]].PhotonVector()).Pt();
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
//        [&](RVec<AbstractParticle>& reco_photons_matched, RVec<unsigned long>& indices)
//        {
//            return Take(reco_photons_matched, indices);
//
//        }, {"photons_pass_cuts", "chosen_two_indices"})
//        .Define("reconstructed_mass",
//        [&](RVec<AbstractParticle>& diph)
//        {
//            return (diph[0].PhotonVector()+diph[1].PhotonVector()).M()/1e3;
//        }, {"chosen_two"});
//
//        auto totEventWeight = resolved_reco_photons_matched.Define("totEventWeight", //New column: weight factor for events in RDF `resolved`
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
//        auto SB = resolved_reco_photons_matched.Filter(
//        [](const double reconstructed_mass)
//        {
//            return (reconstructed_mass < 110) || (reconstructed_mass > 140);
//        }, {"reconstructed_mass"});
//
//        auto SR = resolved_reco_photons_matched.Filter(
//        [](const double reconstructed_mass)
//        {
//            return (reconstructed_mass >= 110) && (reconstructed_mass <= 140);
//        }, {"reconstructed_mass"});
//
//
////        Totals.push_back(df.Count());
////        Totals.push_back(EventWeight.Sum<RVec<float>>("EventWeight"));
//        Totals.push_back(EventWeight.Sum<float>("EventWeight"));
////        resultmaps.push_back(VariationsFor(merged_reco_photons_matched.Count()));
//        resultmaps.push_back(VariationsFor(totEventWeight.Sum<float>("totEventWeight")));
//
//    }
//
////EG_RESOLUTION_ALL__1down:           42864.7 ID 1.00524  ISO 1.02283  TRIG 1.04157
////EG_RESOLUTION_ALL__1up:             42751.5 ID 1.00524  ISO 1.02283  TRIG 1.04157
////EG_SCALE_ALL__1down:                42521.1 ID 1.00524  ISO 1.02283  TRIG 1.04157
////EG_SCALE_ALL__1up:                  43059.8 ID 1.00524  ISO 1.02283  TRIG 1.04157
////PH_EFF_ID_Uncertainty__1down:       42790.4 ID 0.992328 ISO 1.02283  TRIG 1.04157
////PH_EFF_ID_Uncertainty__1up:         42790.4 ID 1.01816  ISO 1.02283  TRIG 1.04157
////PH_EFF_ISO_Uncertainty__1down:      42790.4 ID 1.00524  ISO 1.00971  TRIG 1.04157
////PH_EFF_ISO_Uncertainty__1up:        42790.4 ID 1.00524  ISO 1.03595  TRIG 1.04157
////PH_EFF_TRIGGER_Uncertainty__1down:  42790.4 ID 1.00524  ISO 1.02283  TRIG 1.01646
////PH_EFF_TRIGGER_Uncertainty__1up:    42790.4 ID 1.00524  ISO 1.02283  TRIG 1.06668
//
//    ROOT::RDF::RunGraphs(Totals); // running all computation nodes concurrently
//
////    for (auto& i: GeneratorWeightCounts)
////    {
////        std::cout << *i << '\n';
////    }
////    int count = 0;
////    for (auto& i: resultmaps)
////    {
////        std::cout << count++ << "\n==\n";
////        for (auto& j: i.GetKeys())
////        {
////            std::cout << j << '\n';
////        }
////        std::cout << '\n';
////    }
////          resultmaps
////          ----------
////    0       1       //Z-gamma
////    2       3       //Z-gamma
////    4       5       //Z-gamma
////    6       7       //data
////    8       9       //ma1
////    10      11      //ma2
////    12      13      //ma3
////    14      15      //ma5
////    16      17      //ma9
////    18      19      //Z-jets
////    20      21      //Z-jets
////    22      23      //Z-jets
////    24      25      //Z-jets
////    26      27      //Z-jets
////    28      29      //Z-jets
////    30      31      //Z-jets
////    32      33      //Z-jets
////    34      35      //Z-jets
//
//    ss << R"--(\section*{Table 24})--" << '\n';
//    ss << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
//    ss << R"--(\setlength\extrarowheight{2pt}\renewcommand{\arraystretch}{1.5})--" << '\n';
////    ss << R"--(\begin{tabular}{|c|c|c|c|c|c|})--" << '\n';
//    ss << R"--(\begin{tabular}{|c|c|c|c|})--" << '\n';
//    ss << R"--(\hline)--" << '\n';
//
//        ss << R"--(\multicolumn{4}{|c|}{\parbox{\linewidth}{\centering Resolved Photon Category: Down Variations \\ (\% difference from nominal)}}\\[5 pt])--" << '\n';
//    ss << R"--(\hline)--" << '\n';
//
//    ss << R"--({Sample} & EG\_RESOLUTION\_ALL & EG\_SCALE\_ALL & PH\_EFF\_ID\_Uncertainty \\ \hline)--" << '\n';
//
//    double finalScaleVal;
//
//    double ZgammaNominal = 0, ZgammaEG_RESOLUTION_ALL = 0,
//    ZgammaEG_SCALE_ALL = 0, ZgammaPH_EFF_ISO_Uncertainty = 0,
//    ZgammaPH_EFF_ID_Uncertainty = 0, ZgammaPH_EFF_TRIGGER_Uncertainty = 0;
//
//    double ZjetsNominal = 0, ZjetsEG_RESOLUTION_ALL = 0,
//    ZjetsEG_SCALE_ALL = 0, ZjetsPH_EFF_ISO_Uncertainty = 0,
//    ZjetsPH_EFF_ID_Uncertainty = 0, ZjetsPH_EFF_TRIGGER_Uncertainty = 0;
//
//    double totbkgNominal = 0, totbkgEG_RESOLUTION_ALL = 0,
//    totbkgEG_SCALE_ALL = 0, totbkgPH_EFF_ISO_Uncertainty = 0,
//    totbkgPH_EFF_ID_Uncertainty = 0, totbkgPH_EFF_TRIGGER_Uncertainty = 0;
//
//    for (auto i = 0; (i < resultmaps.size()); i++)
//    {
//        double nominalVal = static_cast<double>(resultmaps[i]["nominal"]);
//
////        auto denominator = *Totals[i].GetResultPtr<ULong64_t>();
//        auto denominator = *Totals[i].GetResultPtr<float>();
//
//        if (i >= 0 && i <= 2) //Zgamma
//        {
//            finalScaleVal = SFs[i]/denominator;
//            ZgammaNominal += finalScaleVal*nominalVal;
//
//            ZgammaEG_RESOLUTION_ALL += finalScaleVal*resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1down"];
//            ZgammaEG_SCALE_ALL += finalScaleVal*resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1down"];
////            ZgammaPH_EFF_ISO_Uncertainty += finalScaleVal*resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1down"];
//            ZgammaPH_EFF_ID_Uncertainty += finalScaleVal*resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1down"];
////            ZgammaPH_EFF_TRIGGER_Uncertainty += finalScaleVal*resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down"];
//        }
//
//        else if (i >= 9)
//        {
//            finalScaleVal = JetNumeratorSFs[i-9]/denominator;
//            ZjetsNominal += finalScaleVal*nominalVal;
//
//            ZjetsEG_RESOLUTION_ALL += finalScaleVal*resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1down"];
//            ZjetsEG_SCALE_ALL += finalScaleVal*resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1down"];
////            ZjetsPH_EFF_ISO_Uncertainty += finalScaleVal*resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1down"];
//            ZjetsPH_EFF_ID_Uncertainty += finalScaleVal*resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1down"];
////            ZjetsPH_EFF_TRIGGER_Uncertainty += finalScaleVal*resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down"];
//        }
//
//        ss << prefixes[i] << " & ";
//        ss << std::setprecision(4) << std::fixed
//        << ( (nominalVal && resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1down"]) ? ((resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1down"]-nominalVal)/nominalVal)*100.0 : 0.0)
//        << " & " << std::setprecision(4) << std::fixed
//        << ( (nominalVal && resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1down"]) ? ((resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1down"]-nominalVal)/nominalVal)*100.0 : 0.0)
////        << " & " << std::setprecision(4) << std::fixed
////        << ( (nominalVal && resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1down"]) ? ((resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1down"]-nominalVal)/nominalVal)*100.0 : 0.0)
//        << " & " << std::setprecision(4) << std::fixed
//        << ( (nominalVal && resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1down"]) ? ((resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1down"]-nominalVal)/nominalVal)*100.0 : 0.0)
////        << " & " << std::setprecision(4) << std::fixed
////        << ( (nominalVal && resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down"]) ? ((resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down"]-nominalVal)/nominalVal)*100.0 : 0.0)
//        << R"--( \\ \hline)--" << '\n';
//
////        if (!resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1down"])
////        {
////            std::cout << i << ": photon_iso_eff:PH_EFF_ISO_Uncertainty__1down not found\n";
////        }
////        else
////        {
////            std::cout << i << ':' << ' ' << resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1down"] << ' ' << nominalVal << '\n';
////        }
//        if (!resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1down"])
//        {
//            std::cout << i << ": photon_id_eff:PH_EFF_ID_Uncertainty__1down not found\n";
//        }
//        else
//        {
//            std::cout << i << ':' << ' ' << resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1down"] << ' ' << nominalVal << '\n';
//        }
////        if (!resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down"])
////        {
////            std::cout << i << ": photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down not found\n";
////        }
////        else
////        {
////            std::cout << i << ':' << ' ' << resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down"] << ' ' << nominalVal << '\n';
////        }
//    }
//
//    ss << R"--(Total $Z\gamma\gamma$ & )--";
//    ss << std::setprecision(4) << std::fixed
//    << ( (ZgammaNominal && ZgammaEG_RESOLUTION_ALL) ? ((ZgammaEG_RESOLUTION_ALL-ZgammaNominal)/ZgammaNominal)*100.0 : 0.0)
//    << " & " << std::setprecision(4) << std::fixed
//    << ( (ZgammaNominal && ZgammaEG_SCALE_ALL) ? ((ZgammaEG_SCALE_ALL-ZgammaNominal)/ZgammaNominal)*100.0 : 0.0)
////    << " & " << std::setprecision(4) << std::fixed
////    << ( (ZgammaNominal && ZgammaPH_EFF_ISO_Uncertainty) ? ((ZgammaPH_EFF_ISO_Uncertainty-ZgammaNominal)/ZgammaNominal)*100.0 : 0.0)
//    << " & " << std::setprecision(4) << std::fixed
//    << ( (ZgammaNominal && ZgammaPH_EFF_ID_Uncertainty) ? ((ZgammaPH_EFF_ID_Uncertainty-ZgammaNominal)/ZgammaNominal)*100.0 : 0.0)
////    << " & " << std::setprecision(4) << std::fixed
////    << ( (ZgammaNominal && ZgammaPH_EFF_TRIGGER_Uncertainty) ? ((ZgammaPH_EFF_TRIGGER_Uncertainty-ZgammaNominal)/ZgammaNominal)*100.0 : 0.0)
//    << R"--( \\ \hline)--" << '\n';
//
//    ss << R"--(Total $Z$+jets & )--";
//    ss << std::setprecision(4) << std::fixed
//    << ( (ZjetsNominal && ZjetsEG_RESOLUTION_ALL) ? ((ZjetsEG_RESOLUTION_ALL-ZjetsNominal)/ZjetsNominal)*100.0 : 0.0)
//    << " & " << std::setprecision(4) << std::fixed
//    << ( (ZjetsNominal && ZjetsEG_SCALE_ALL) ? ((ZjetsEG_SCALE_ALL-ZjetsNominal)/ZjetsNominal)*100.0 : 0.0)
////    << " & " << std::setprecision(4) << std::fixed
////    << ( (ZjetsNominal && ZjetsPH_EFF_ISO_Uncertainty) ? ((ZjetsPH_EFF_ISO_Uncertainty-ZjetsNominal)/ZjetsNominal)*100.0 : 0.0)
//    << " & " << std::setprecision(4) << std::fixed
//    << ( (ZjetsNominal && ZjetsPH_EFF_ID_Uncertainty) ? ((ZjetsPH_EFF_ID_Uncertainty-ZjetsNominal)/ZjetsNominal)*100.0 : 0.0)
////    << " & " << std::setprecision(4) << std::fixed
////    << ( (ZjetsNominal && ZjetsPH_EFF_TRIGGER_Uncertainty) ? ((ZjetsPH_EFF_TRIGGER_Uncertainty-ZjetsNominal)/ZjetsNominal)*100.0 : 0.0)
//    << R"--( \\ \hline)--" << '\n';
//
//    totbkgNominal = ZgammaNominal + ZjetsNominal;
//
//    totbkgEG_RESOLUTION_ALL = ZgammaEG_RESOLUTION_ALL + ZjetsEG_RESOLUTION_ALL;
//    totbkgEG_SCALE_ALL = ZgammaEG_SCALE_ALL + ZjetsEG_SCALE_ALL;
//    totbkgPH_EFF_ISO_Uncertainty = ZgammaPH_EFF_ISO_Uncertainty + ZjetsPH_EFF_ISO_Uncertainty;
//    totbkgPH_EFF_ID_Uncertainty = ZgammaPH_EFF_ID_Uncertainty + ZjetsPH_EFF_ID_Uncertainty;
//    totbkgPH_EFF_TRIGGER_Uncertainty = ZgammaPH_EFF_TRIGGER_Uncertainty + ZjetsPH_EFF_TRIGGER_Uncertainty;
//
//    ss << R"--(Total Bkg & )--";
//    ss << std::setprecision(4) << std::fixed
//    << ( (totbkgNominal && totbkgEG_RESOLUTION_ALL) ? ((totbkgEG_RESOLUTION_ALL-totbkgNominal)/totbkgNominal)*100.0 : 0.0)
//    << " & " << std::setprecision(4) << std::fixed
//    << ( (totbkgNominal && totbkgEG_SCALE_ALL) ? ((totbkgEG_SCALE_ALL-totbkgNominal)/totbkgNominal)*100.0 : 0.0)
////    << " & " << std::setprecision(4) << std::fixed
////    << ( (totbkgNominal && totbkgPH_EFF_ISO_Uncertainty) ? ((totbkgPH_EFF_ISO_Uncertainty-totbkgNominal)/totbkgNominal)*100.0 : 0.0)
//    << " & " << std::setprecision(4) << std::fixed
//    << ( (totbkgNominal && totbkgPH_EFF_ID_Uncertainty) ? ((totbkgPH_EFF_ID_Uncertainty-totbkgNominal)/totbkgNominal)*100.0 : 0.0)
////    << " & " << std::setprecision(4) << std::fixed
////    << ( (totbkgNominal && totbkgPH_EFF_TRIGGER_Uncertainty) ? ((totbkgPH_EFF_TRIGGER_Uncertainty-totbkgNominal)/totbkgNominal)*100.0 : 0.0)
//    << R"--( \\ \hline)--" << '\n';
//
//    ss << R"--(\end{tabular}})--" << '\n';
//
//    ss << "\n\n\n";
//
//    std::ofstream out("Table24.txt");
//    out << ss.str() << '\n';
//    out.close();
//}

void Merged_Up_Higgs_Mass()
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
//        "PH_EFF_ISO_Uncertainty__1up",
        "PH_EFF_ID_Uncertainty__1up",
//        "PH_EFF_TRIGGER_Uncertainty__1up",
//        "EG_RESOLUTION_ALL__1down",
    };

    std::vector<std::vector<std::string>> input_filenames = {
        //Signal
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root"}, //1 GeV
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600751.PhPy8EG_AZNLO_ggH125_mA2p0_v1.root"}, // 2 GeV
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600752.PhPy8EG_AZNLO_ggH125_mA3p0_v1.root"}, // 3 GeV
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root"}, //5 GeV
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600756.PhPy8EG_AZNLO_ggH125_mA9p0_v1.root"}, // 9 GeV
        //Data
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_data_test.root"},
        //Z gamma background
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617070._000001.LGNTuple.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617064._000001.LGNTuple.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617074._000001.LGNTuple.root"},
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

    std::vector<const char*> prefixes = {"Sig m_{A} = 1 GeV", "Sig m_{A} = 2 GeV", "Sig m_{A} = 3 GeV", "Sig m_{A} = 5 GeV", "Sig m_{A} = 9 GeV", "data", "pty2_9_17", "pty_17_myy_0_80", "pty_17_myy_80", "Zee_lightJet_0-70", "Zee_lightJet_70-140", "Zee_lightJet_140-280", "Zee_cJet_0-70", "Zee_cJet_70-140", "Zee_cJet_140-280", "Zee_bJet_0-70", "Zee_bJet_70-140", "Zee_bJet_140-280"};

    std::vector<ROOT::RDF::RResultHandle> Nodes;

    int count = 0;

    std::vector<RResultMap<float>> resultmaps;
    std::vector<RResultMap<TH1D>> resultmapHists;

    auto verbosity = ROOT::Experimental::RLogScopedVerbosity(ROOT::Detail::RDF::RDFLogChannel(), ROOT::Experimental::ELogLevel::kInfo);

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
        [](RVec<AbstractParticle> electrons)
        {
            //keep the electrons in each event that have pt > 20 GeV, |η| < 2.37,
            //|η| not between 1.37 and 1.52, and that satisfy a medium id criteria `electron_id_medium`
            electrons.erase(std::remove_if(electrons.begin(),electrons.end(),
            [](AbstractParticle& ep)
            {
                 return (!((ep.electron_pt/1e3 > 20) && (std::abs(ep.electron_eta) < 2.37) &&
                 (!((1.37 < std::abs(ep.electron_eta)) && (std::abs(ep.electron_eta) < 1.52)))
                 && (ep.electron_id_medium == 1)));

            }), electrons.end());

            return electrons;

        },{"abstract_electrons"}).Filter(
        [](RVec<Muon>& muons, RVec<AbstractParticle>& electrons)
        {
            return (electrons.size()==2 && muons.empty()); //keep events which have exactly 2 electrons for di_electrons and no muons

        }, {"muons", "di_electrons"});

        //new dataframe node: contains only the events from `two_leptons` whose electrons in the `di_electrons` branch have opposite charge
        auto opp_charge = two_leptons
        .Filter([](RVec<AbstractParticle>& electrons)
        {
            return (electrons[0].electron_charge*electrons[1].electron_charge < 0);

        }, {"di_electrons"});

        //new dataframe node: contains only the events from `opp_charge` that have 1 electron with pt > 20 GeV and the other with pt > 27 GeV
        auto leadingPt = opp_charge.Filter([](RVec<AbstractParticle>& electrons)
        {
            return ((electrons[0].electron_pt > 20e3 && electrons[1].electron_pt > 27e3) || (electrons[1].electron_pt > 20e3 && electrons[0].electron_pt > 27e3));
        }, {"di_electrons"});

        auto deltaR = leadingPt.Filter([] (RVec<AbstractParticle>& electrons)
        {
            return (DeltaR(electrons[0].ElectronVector(), electrons[1].ElectronVector()) > 0.01);
        }, {"di_electrons"});

        //new dataframe node: Contains a new column `dilep` in addition to the ones in `same_flavour` that stores the di-electron four-vector
        auto dilep = deltaR.Define("dilep",[] (RVec<AbstractParticle>& electrons)
        {
            return (electrons[0].ElectronVector() + electrons[1].ElectronVector());
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
        [&](const RVec<AbstractParticle>& photons)
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
        }, {"abstract_photons"})
        .Define("photons_pass_cuts",
        [&](RVec<AbstractParticle>& photons, RVec<int>& photon_indices)
        {
            return Take(photons, photon_indices);
        }, {"abstract_photons", "photons_pass_cut_indices"});

        //New dataframe node: contains only events from `photon_passes_cuts` that fail the resolved category and pass the merged category
        auto merged_reco_photons_matched = photon_passes_cuts.Filter(
        [&](RVec<AbstractParticle>& reco_photons_matched)
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
               delta_r = DeltaR(reco_photons_matched[combs[0][i]].PhotonVector(), reco_photons_matched[combs[1][i]].PhotonVector());
               m = (reco_photons_matched[combs[0][i]].PhotonVector() + reco_photons_matched[combs[1][i]].PhotonVector()).M();
               pt = (reco_photons_matched[combs[0][i]].PhotonVector() + reco_photons_matched[combs[1][i]].PhotonVector()).Pt();
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
        .Define("merged_photon", //new column: The reco-photon corresponding to `merged_photon_index`
        [&](const RVec<AbstractParticle>& reco_photons_matched, int merged_photon_index)
        {
            return reco_photons_matched[merged_photon_index];

        }, {"photons_pass_cuts", "merged_photon_index"});

        auto dilepton_and_photon = merged_reco_photons_matched
        .Define("reconstructed_mass",[&](RVec<AbstractParticle>& di_electrons, AbstractParticle& merged_photon)
        {
            auto four_momentum = di_electrons[0].ElectronVector() + di_electrons[1].ElectronVector();

            return (four_momentum + merged_photon.PhotonVector()).M()/1e3;

        }, {"di_electrons", "merged_photon"})
        .Define("totEventWeight", [](RVec<float> totEventWeightVec, RVec<int>& photons_pass_cut_indices, int mpi)
        {
            return Take(totEventWeightVec, photons_pass_cut_indices)[mpi];

        }, {"totEventWeightVec", "photons_pass_cut_indices", "merged_photon_index"});

        if (count >= 6)
        {
            Nodes.push_back(EventWeight.Sum<float>("EventWeight"));
        }
        resultmapHists.push_back(VariationsFor(dilepton_and_photon.Histo1D<double>({prefixes[count], prefixes[count++], 100u, 80, 200}, "reconstructed_mass", "totEventWeight")));

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

    ROOT::RDF::RunGraphs(Nodes); // running all computation nodes concurrently

//    resultmapHists
//    ==============
//    0   //ma1
//    1   //ma2
//    2   //ma3
//    3   //ma5
//    4   //ma9
//    5   //data
//    6   //Z-gamma
//    7   //Z-gamma
//    8   //Z-gamma
//    9   //Z-jets
//    10  //Z-jets
//    11  //Z-jets
//    12  //Z-jets
//    13  //Z-jets
//    14  //Z-jets
//    15  //Z-jets
//    16  //Z-jets
//    17  //Z-jets
//
//    resultmaps
//    ==========
//
//    0       //Z-gamma
//    1       //Z-gamma
//    2       //Z-gamma
//    3       //Z-jets
//    4       //Z-jets
//    5       //Z-jets
//    6       //Z-jets
//    7       //Z-jets
//    8       //Z-jets
//    9       //Z-jets
//    10      //Z-jets
//    11      //Z-jets

    constexpr std::array<const char*, 3> systematics = {"photons_and_electrons:EG_RESOLUTION_ALL__1up", "photons_and_electrons:EG_SCALE_ALL__1up", "photon_id_eff:PH_EFF_ID_Uncertainty__1up"};
    constexpr std::array<const char*, 4> nominal_systematics =
    {
        "nominal",
        "photons_and_electrons:EG_RESOLUTION_ALL__1up",
        "photons_and_electrons:EG_SCALE_ALL__1up",
        "photon_id_eff:PH_EFF_ID_Uncertainty__1up",
    };
    constexpr std::array<const char*, 3> Systematics = {"EG_RESOLUTION_ALL__1up", "EG_SCALE_ALL__1up", "PH_EFF_ID_Uncertainty__1up"};
    constexpr std::array<const char*, 4> nominalSystematics = {"nominal", "EG_RESOLUTION_ALL__1up", "EG_SCALE_ALL__1up", "PH_EFF_ID_Uncertainty__1up"};
    constexpr std::array<const char*, 6> signalPlusDataSamples = {"Sig m_{A} = 1 GeV", "Sig m_{A} = 2 GeV", "Sig m_{A} = 3 GeV", "Sig m_{A} = 5 GeV", "Sig m_{A} = 9 GeV", "Data"};
    constexpr std::array<const char*, 6> signalPlusDataFileNames = {"prompt_ma1_merged_higgs_syst_up.pdf", "prompt_ma2_merged_higgs_syst_up.pdf", "prompt_ma3_merged_higgs_syst_up.pdf", "prompt_ma5_merged_higgs_syst_up.pdf", "prompt_ma9_merged_higgs_syst_up.pdf", "data_merged_higgs_syst_up.pdf"};
    constexpr std::array<EColor, 3> colors = {kBlue, kRed, static_cast<EColor>(kOrange+1)};
    constexpr std::array<EColor, 4> nominalColors = {kBlack, kBlue, kRed, static_cast<EColor>(kOrange+1)};
    constexpr std::array<float, 6> signalPlusDataResidualMaxima = {2.1, 1.15, 1.3, 0.89, 1.8, 1.05};

    TCanvas* c1;
    TLegend* legend;
    TLatex Tl;

    //signal+data
    for (int i = 0; i < 6; i++) //signal samples plus data
    {
        c1 = new TCanvas();
        legend = new TLegend(0.5, 0.2, 0.85, 0.6);
        c1->cd();
        TPad* pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 0.95);
        pad1->SetTopMargin(0.1);
        pad1->SetBottomMargin(0);
        pad1->Draw();
        pad1->cd();

        gStyle->SetOptStat(0);
        resultmapHists[i]["nominal"].SetLineColor(kBlack);
        resultmapHists[i]["nominal"].SetTitle(";m_{ll#gamma} [GeV];Events");
        resultmapHists[i]["nominal"].SetTitle(signalPlusDataSamples[i]);
        resultmapHists[i]["nominal"].GetXaxis()->SetTitle("m_{ll#gamma} [GeV]");
        resultmapHists[i]["nominal"].GetXaxis()->SetTitleOffset(1.2);
        resultmapHists[i]["nominal"].GetYaxis()->SetTitle("Events");
        resultmapHists[i]["nominal"].GetYaxis()->CenterTitle(true);
        legend->AddEntry(&resultmapHists[i]["nominal"], "nominal", "l");

        resultmapHists[i]["nominal"].Draw("HIST");

        for (int j = 0; j < 3; j++)
        {
            std::cout << "difference = "
            << resultmapHists[i][systematics[j]].Integral() - resultmapHists[i]["nominal"].Integral() << '\n';
            resultmapHists[i][systematics[j]].SetLineColor(colors[j]);
            legend->AddEntry(&resultmapHists[i][systematics[j]], Systematics[j], "l");
            resultmapHists[i][systematics[j]].Draw("HISTsame");
        }

        Tl.SetTextSize(0.03);
        Tl.DrawLatexNDC(0.6, 0.8, "#it{ATLAS} Internal");
        Tl.DrawLatexNDC(0.6, 0.7,"#sqrt{s} = 13 TeV");
        legend->SetBorderSize(0);
        legend->Draw("same");

        // Create and draw the error plot
        TH1F* errorHist;
        double error;

        c1->cd();
        TPad* pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.25);
        pad2->SetTopMargin(0);
        pad2->SetBottomMargin(0.25);
        pad2->Draw();
        pad2->cd();
        legend = new TLegend(0.675, 0.7, 0.875, 0.9);

        for (int j = 0; j < 3; j++)
        {
            errorHist = new TH1F(Form("errorHist_%d%d", i, j), "", resultmapHists[i]["nominal"].GetNbinsX(),
                                       resultmapHists[i]["nominal"].GetXaxis()->GetXmin(),
                                       resultmapHists[i]["nominal"].GetXaxis()->GetXmax());
            errorHist->SetLineColor(colors[j]);
            errorHist->SetLineWidth(1);
            errorHist->SetMarkerStyle(20);
            errorHist->SetMarkerSize(0.9);
            errorHist->GetXaxis()->SetLabelSize(0.1);
            errorHist->GetYaxis()->SetLabelSize(0.1);
            errorHist->GetXaxis()->SetTitleSize(0.1);
            errorHist->GetYaxis()->SetTitleSize(0.1);
            errorHist->SetMarkerColor(colors[j]);

            for (int bin = 1; bin <= resultmapHists[i]["nominal"].GetNbinsX(); bin++)
            {
                double nominal = resultmapHists[i]["nominal"].GetBinContent(bin);
                if (!nominal)
                {
                    error = 0;
                }
                else
                {
                    double systematic = resultmapHists[i][systematics[j]].GetBinContent(bin);
                    error = (systematic - nominal) / nominal;
                }
                errorHist->SetBinContent(bin, error);
            }

            if (j == 0)
            {
                errorHist->SetTitle(";m_{ll#gamma} [GeV];#frac{syst-nominal}{nominal}");
                errorHist->GetXaxis()->SetTitle("m_{ll#gamma} [GeV]");
                errorHist->GetXaxis()->SetTitleOffset(1.2);
                errorHist->GetYaxis()->SetTitle("#frac{syst-nominal}{nominal}");
                errorHist->GetYaxis()->SetTitleOffset(0.4);
                errorHist->GetYaxis()->CenterTitle(true);
                errorHist->SetMaximum(signalPlusDataResidualMaxima[i]);
                errorHist->Draw("histsame");

            }
            else
            {
                c1->cd();
                pad2->cd();
                errorHist->Draw("histsame");
            }

            legend->AddEntry(errorHist, Systematics[j], "l");

            c1->cd();
            c1->Update();

        }
        pad2->cd();
        legend->SetBorderSize(0);
        legend->Draw("same");
        c1->SaveAs(signalPlusDataFileNames[i]);
    }

    constexpr std::array<double,3> SFs = {((139e15)*(.871e-12)),((139e15)*(.199e-12)), ((139e15)*(.0345e-15))}; //numerators for Z-gamma bkg


    //Z-gamma
    c1 = new TCanvas();
    legend = new TLegend(0.65, 0.43, 0.875, 0.7);
    TH1D *hist, *nominalhist;
    gStyle->SetOptStat(0);
    c1->cd();
    TPad* pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 0.95);
    pad1->SetTopMargin(0.1);
    pad1->SetBottomMargin(0);
    pad1->Draw();
    pad1->cd();
    gStyle->SetOptStat(0);
    c1->cd();
    TPad* pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.25);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.25);
    pad2->Draw();
    pad2->cd();
    TLegend* errorLegend = new TLegend(0.475, 0.7, 0.675, 0.9);
    TH1F* errorHist;
    double error;

    for (int i = 0; i < 4; i++) //systematics: Z-gamma
    {
        pad1->cd();
        hist = new TH1D(Form("Z#gamma#gamma%d",i),Form("Z#gamma#gamma%d",i), 100u, 80, 200);
        for (int j = 0, k = 6; k <= 8; j++, k++)
        {
            hist->Add(&resultmapHists[k][nominal_systematics[i]], SFs[j] / *Nodes[j].GetResultPtr<float>());
        }
        hist->SetLineColor(nominalColors[i]);
        std::cout <<  hist->Integral() << '\n';

        if (i == 0)
        {
            nominalhist = static_cast<TH1D*>(hist->Clone());
            hist->SetTitle(";m_{ll#gamma} [GeV];Events");
            hist->SetTitle("Z#gamma#gamma");
            hist->GetXaxis()->SetTitle("m_{ll#gamma} [GeV]");
            hist->GetXaxis()->SetTitleOffset(1.2);
            hist->GetYaxis()->SetTitle("Events");
            hist->GetYaxis()->CenterTitle(true);
            hist->Draw("HIST");
        }
        else
        {
            
            hist->Draw("HISTsame");
            
            //Now do the error plot
            
            c1->cd();
            pad2->cd();
            
            
            
            errorHist = new TH1F(Form("errorHist_%d", i), "", nominalhist->GetNbinsX(),
                                 nominalhist->GetXaxis()->GetXmin(),
                                 nominalhist->GetXaxis()->GetXmax());
            errorHist->SetLineColor(colors[i-1]);
            errorHist->SetLineWidth(1);
            errorHist->SetMarkerStyle(20);
            errorHist->SetMarkerSize(0.9);
            errorHist->GetXaxis()->SetLabelSize(0.1);
            errorHist->GetYaxis()->SetLabelSize(0.1);
//            errorHist->GetYaxis()->SetLabelOffset(0.02);
//            errorHist->GetYaxis()->SetNoExponent();
            
            errorHist->GetXaxis()->SetTitleSize(0.1);
            errorHist->GetYaxis()->SetTitleSize(0.1);
            errorHist->SetMarkerColor(nominalColors[i]);

            for (int bin = 1; bin <= nominalhist->GetNbinsX(); bin++)
            {
                double nominal = nominalhist->GetBinContent(bin);
                if (!nominal)
                {
                    error = 0;
                }
                else
                {
                    double systematic = hist->GetBinContent(bin);
                    error = (systematic - nominal) / nominal;
                }
                errorHist->SetBinContent(bin, error);
            }

            if (i == 1)
            {
                errorHist->SetTitle(";m_{ll#gamma} [GeV];#frac{syst-nominal}{nominal}");
                errorHist->GetXaxis()->SetTitle("m_{ll#gamma} [GeV]");
                errorHist->GetXaxis()->SetTitleOffset(1.2);
                errorHist->GetYaxis()->SetTitle("#frac{syst-nominal}{nominal}");
                errorHist->GetYaxis()->SetTitleOffset(0.4);
                errorHist->GetYaxis()->CenterTitle(true);
                errorHist->SetMaximum(1.9);
                errorHist->SetMinimum(-0.85);
                errorHist->Draw("histsame");

            }
            else
            {
                c1->cd();
                pad2->cd();
                errorHist->Draw("histsame");
            }

            errorLegend->AddEntry(errorHist, nominalSystematics[i], "l");

            c1->cd();
            c1->Update();
            pad2->cd();
            errorLegend->SetBorderSize(0);
            errorLegend->Draw("same");

        }
        legend->AddEntry(hist, nominalSystematics[i], "l");
    }
    pad1->cd();
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.6, 0.85, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.6, 0.775,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
    legend->SetBorderSize(0);
    legend->Draw("same");

    c1->SaveAs("Zgamma_merged_higgs_syst_up.pdf");
    //Z-jets
    std::array<double,9> JetNumeratorSFs = {((139e15)*(1.9828e-9)*(0.821204)),((139e15)*(110.64e-12)*(0.69275)),((139e15)*(40.645e-12)*(0.615906)),((139e15)*(1.9817e-9)*(0.1136684)),((139e15)*(110.47e-12)*(0.1912956)),((139e15)*(40.674e-12)*(0.2326772)),((139e15)*(1.9819e-9)*(0.0656969)),((139e15)*(110.53e-12)*(0.1158741)),((139e15)*(40.68e-12)*(0.1535215))}; //numerators for jet bkg

    c1 = new TCanvas();
    legend = new TLegend(0.56, 0.25, 0.85, 0.65);
    gStyle->SetOptStat(0);
    c1->cd();
    pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 0.95);
    pad1->SetTopMargin(0.1);
    pad1->SetBottomMargin(0);
    pad1->Draw();
    pad1->cd();
    gStyle->SetOptStat(0);
    c1->cd();
    pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.25);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.25);
    pad2->Draw();
    pad2->cd();
    errorLegend = new TLegend(0.375, 0.7, 0.575, 0.9);

    for (int i = 0; i < 4; i++) //systematics: Z-jets
    {
        pad1->cd();
        hist = new TH1D(Form("Z+Jets%d",i),Form("Z+Jets%d",i), 100u, 80, 200);
        for (int j = 0, k = 9; k <= 17; j++, k++)
        {
            hist->Add(&resultmapHists[k][nominal_systematics[i]], JetNumeratorSFs[j] / *Nodes[j+3].GetResultPtr<float>());
        }
        hist->SetLineColor(nominalColors[i]);
        std::cout <<  hist->Integral() << '\n';

        if (i == 0)
        {
            nominalhist = static_cast<TH1D*>(hist->Clone());
            hist->SetTitle(";m_{ll#gamma} [GeV];Events");
            hist->SetTitle("Z+Jets");
            hist->GetXaxis()->SetTitle("m_{ll#gamma} [GeV]");
            hist->GetXaxis()->SetTitleOffset(1.2);
            hist->GetYaxis()->SetTitle("Events");
            hist->GetYaxis()->CenterTitle(true);
            hist->Draw("HIST");
        }
        else
        {
            hist->Draw("HISTsame");
            //Now do the error plot

            c1->cd();
            pad2->cd();



            errorHist = new TH1F(Form("errorHist_%d", i), "", nominalhist->GetNbinsX(),
                                 nominalhist->GetXaxis()->GetXmin(),
                                 nominalhist->GetXaxis()->GetXmax());
            errorHist->SetLineColor(colors[i-1]);
            errorHist->SetLineWidth(1);
            errorHist->SetMarkerStyle(20);
            errorHist->SetMarkerSize(0.9);
            errorHist->GetXaxis()->SetLabelSize(0.1);
            errorHist->GetYaxis()->SetLabelSize(0.1);
            errorHist->GetXaxis()->SetTitleSize(0.1);
            errorHist->GetYaxis()->SetTitleSize(0.1);
            errorHist->SetMarkerColor(nominalColors[i]);

            for (int bin = 1; bin <= nominalhist->GetNbinsX(); bin++)
            {
                double nominal = nominalhist->GetBinContent(bin);
                if (!nominal)
                {
                    error = 0;
                }
                else
                {
                    double systematic = hist->GetBinContent(bin);
                    error = (systematic - nominal) / nominal;
                }
                errorHist->SetBinContent(bin, error);
            }

            if (i == 1)
            {
                errorHist->SetTitle(";m_{ll#gamma} [GeV];#frac{syst-nominal}{nominal}");
                errorHist->GetXaxis()->SetTitle("m_{ll#gamma} [GeV]");
                errorHist->GetXaxis()->SetTitleOffset(1.2);
                errorHist->GetYaxis()->SetTitle("#frac{syst-nominal}{nominal}");
                errorHist->GetYaxis()->SetTitleOffset(0.4);
                errorHist->GetYaxis()->CenterTitle(true);
                errorHist->SetMaximum(1.4);
                errorHist->SetMinimum(-0.7);
                errorHist->Draw("histsame");

            }
            else
            {
                c1->cd();
                pad2->cd();
                errorHist->Draw("histsame");
            }

            errorLegend->AddEntry(errorHist, nominalSystematics[i], "l");

            c1->cd();
            c1->Update();
            pad2->cd();
            errorLegend->SetBorderSize(0);
            errorLegend->Draw("same");
        }
        legend->AddEntry(hist, nominalSystematics[i], "l");
    }
    pad1->cd();
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.6, 0.85, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.6, 0.775,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
    legend->SetBorderSize(0);
    legend->Draw("same");
    c1->SaveAs("Zjets_merged_higgs_syst_up.pdf");

}

//void Merged_Down_Higgs_Mass()
//{
//    Event::systematics =
//    {
////        "PH_EFF_ISO_Uncertainty",
////        "PH_EFF_ISO_Uncertainty",
////        "EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR",
////        "PRW_DATASF",
////        "MUON_EFF_RECO_SYS",
////        "MUON_EFF_ISO_SYS",
////        "MUON_EFF_TrigSystUncertainty",
////        "EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR",
////        "MUON_EFF_TrigStatUncertainty",
////        "MUON_EFF_RECO_STAT",
////        "MUON_EFF_TTVA_STAT",
////        "EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR",
////        "EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR",
////        "MUON_EFF_TTVA_SYS",
////        "MUON_EFF_ISO_STAT",
////        "MUON_SAGITTA_RHO",
////        "EG_RESOLUTION_ALL",
////        "EG_SCALE_ALL",
////        "MUON_MS",
////        "MUON_ID",
////        "EL_EFF_TriggerEff_TOTAL_1NPCOR_PLUS_UNCOR",
////        "MUON_SAGITTA_RESBIAS",
////        "MUON_SCALE",
////        "EG_RESOLUTION_ALL__1up",
////        "EG_SCALE_ALL__1up",
////        "PH_EFF_ISO_Uncertainty__1up",
////        "PH_EFF_ID_Uncertainty__1up",
////        "PH_EFF_TRIGGER_Uncertainty__1up",
//        "EG_RESOLUTION_ALL__1down",
//        "EG_SCALE_ALL__1down",
////        "PH_EFF_ISO_Uncertainty__1down",
//        "PH_EFF_ID_Uncertainty__1down",
////        "PH_EFF_TRIGGER_Uncertainty__1down",
//    };
//
//    std::vector<std::vector<std::string>> input_filenames = {
//        //Signal
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root"}, //1 GeV
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600751.PhPy8EG_AZNLO_ggH125_mA2p0_v1.root"}, // 2 GeV
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600752.PhPy8EG_AZNLO_ggH125_mA3p0_v1.root"}, // 3 GeV
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root"}, //5 GeV
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600756.PhPy8EG_AZNLO_ggH125_mA9p0_v1.root"}, // 9 GeV
//        //Data
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_data_test.root"},
//        //Z gamma background
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617070._000001.LGNTuple.root"},
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617064._000001.LGNTuple.root"},
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617074._000001.LGNTuple.root"},
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
//    std::vector<const char*> prefixes = {"Sig m_{A} = 1 GeV", "Sig m_{A} = 2 GeV", "Sig m_{A} = 3 GeV", "Sig m_{A} = 5 GeV", "Sig m_{A} = 9 GeV", "data", "pty2_9_17", "pty_17_myy_0_80", "pty_17_myy_80", "Zee_lightJet_0-70", "Zee_lightJet_70-140", "Zee_lightJet_140-280", "Zee_cJet_0-70", "Zee_cJet_70-140", "Zee_cJet_140-280", "Zee_bJet_0-70", "Zee_bJet_70-140", "Zee_bJet_140-280"};
//
//    std::vector<ROOT::RDF::RResultHandle> Nodes;
//
//    int count = 0;
//
//    std::vector<RResultMap<float>> resultmaps;
//    std::vector<RResultMap<TH1D>> resultmapHists;
//
//    auto verbosity = ROOT::Experimental::RLogScopedVerbosity(ROOT::Detail::RDF::RDFLogChannel(), ROOT::Experimental::ELogLevel::kInfo);
//
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
//        [](RVec<float> photon_id_eff, RVec<float> photon_iso_eff, RVec<float> photon_trg_eff/*, RVec<float> ei_event_weights_generator*/)
//        {
//            auto ResizeVal = std::max({photon_id_eff.size(), photon_iso_eff.size(), photon_trg_eff.size()});
//            photon_id_eff.resize(ResizeVal,1);
//            photon_iso_eff.resize(ResizeVal,1);
//            photon_trg_eff.resize(ResizeVal,1);
//
//            return photon_id_eff;
//            // *photon_iso_eff*photon_trg_eff; //Returns a vector of efficiencies
//
//        }, {"photon_id_eff", "photon_iso_eff", "photon_trg_eff",/* "ei_event_weights_generator"*/});
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
//        [](RVec<AbstractParticle> electrons)
//        {
//            //keep the electrons in each event that have pt > 20 GeV, |η| < 2.37,
//            //|η| not between 1.37 and 1.52, and that satisfy a medium id criteria `electron_id_medium`
//            electrons.erase(std::remove_if(electrons.begin(),electrons.end(),
//            [](AbstractParticle& ep)
//            {
//                 return (!((ep.electron_pt/1e3 > 20) && (std::abs(ep.electron_eta) < 2.37) &&
//                 (!((1.37 < std::abs(ep.electron_eta)) && (std::abs(ep.electron_eta) < 1.52)))
//                 && (ep.electron_id_medium == 1)));
//
//            }), electrons.end());
//
//            return electrons;
//
//        },{"abstract_electrons"}).Filter(
//        [](RVec<Muon>& muons, RVec<AbstractParticle>& electrons)
//        {
//            return (electrons.size()==2 && muons.empty()); //keep events which have exactly 2 electrons for di_electrons and no muons
//
//        }, {"muons", "di_electrons"});
//
//        //new dataframe node: contains only the events from `two_leptons` whose electrons in the `di_electrons` branch have opposite charge
//        auto opp_charge = two_leptons
//        .Filter([](RVec<AbstractParticle>& electrons)
//        {
//            return (electrons[0].electron_charge*electrons[1].electron_charge < 0);
//
//        }, {"di_electrons"});
//
//        //new dataframe node: contains only the events from `opp_charge` that have 1 electron with pt > 20 GeV and the other with pt > 27 GeV
//        auto leadingPt = opp_charge.Filter([](RVec<AbstractParticle>& electrons)
//        {
//            return ((electrons[0].electron_pt > 20e3 && electrons[1].electron_pt > 27e3) || (electrons[1].electron_pt > 20e3 && electrons[0].electron_pt > 27e3));
//        }, {"di_electrons"});
//
//        auto deltaR = leadingPt.Filter([] (RVec<AbstractParticle>& electrons)
//        {
//            return (DeltaR(electrons[0].ElectronVector(), electrons[1].ElectronVector()) > 0.01);
//        }, {"di_electrons"});
//
//        //new dataframe node: Contains a new column `dilep` in addition to the ones in `same_flavour` that stores the di-electron four-vector
//        auto dilep = deltaR.Define("dilep",[] (RVec<AbstractParticle>& electrons)
//        {
//            return (electrons[0].ElectronVector() + electrons[1].ElectronVector());
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
//        auto photon_passes_cuts = ptCut
//        .Define("photons_pass_cut_indices",
//        [&](const RVec<AbstractParticle>& photons)
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
//        }, {"abstract_photons"})
//        .Define("photons_pass_cuts",
//        [&](RVec<AbstractParticle>& photons, RVec<int>& photon_indices)
//        {
//            return Take(photons, photon_indices);
//        }, {"abstract_photons", "photons_pass_cut_indices"});
//
//        //New dataframe node: contains only events from `photon_passes_cuts` that fail the resolved category and pass the merged category
//        auto merged_reco_photons_matched = photon_passes_cuts.Filter(
//        [&](RVec<AbstractParticle>& reco_photons_matched)
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
//               delta_r = DeltaR(reco_photons_matched[combs[0][i]].PhotonVector(), reco_photons_matched[combs[1][i]].PhotonVector());
//               m = (reco_photons_matched[combs[0][i]].PhotonVector() + reco_photons_matched[combs[1][i]].PhotonVector()).M();
//               pt = (reco_photons_matched[combs[0][i]].PhotonVector() + reco_photons_matched[combs[1][i]].PhotonVector()).Pt();
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
//        [&](const RVec<AbstractParticle>& rpm) //rpm = reco photons matched
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
//        [&](const RVec<AbstractParticle>& reco_photons_matched, int merged_photon_index)
//        {
//            return reco_photons_matched[merged_photon_index];
//
//        }, {"photons_pass_cuts", "merged_photon_index"});
//
//        auto dilepton_and_photon = merged_reco_photons_matched
//        .Define("reconstructed_mass",[&](RVec<AbstractParticle>& di_electrons, AbstractParticle& merged_photon)
//        {
//            auto four_momentum = di_electrons[0].ElectronVector() + di_electrons[1].ElectronVector();
//
//            return (four_momentum + merged_photon.PhotonVector()).M()/1e3;
//
//        }, {"di_electrons", "merged_photon"})
//        .Define("totEventWeight", [](RVec<float> totEventWeightVec, RVec<int>& photons_pass_cut_indices, int mpi)
//        {
//            return Take(totEventWeightVec, photons_pass_cut_indices)[mpi];
//
//        }, {"totEventWeightVec", "photons_pass_cut_indices", "merged_photon_index"});
//
//        if (count >= 6)
//        {
//            Nodes.push_back(EventWeight.Sum<float>("EventWeight"));
//        }
//        resultmapHists.push_back(VariationsFor(dilepton_and_photon.Histo1D<double>({prefixes[count], prefixes[count++], 100u, 80, 200}, "reconstructed_mass", "totEventWeight")));
//
//    }
//
////EG_RESOLUTION_ALL__1down:           42864.7 ID 1.00524  ISO 1.02283  TRIG 1.04157
////EG_RESOLUTION_ALL__1up:             42751.5 ID 1.00524  ISO 1.02283  TRIG 1.04157
////EG_SCALE_ALL__1down:                42521.1 ID 1.00524  ISO 1.02283  TRIG 1.04157
////EG_SCALE_ALL__1up:                  43059.8 ID 1.00524  ISO 1.02283  TRIG 1.04157
////PH_EFF_ID_Uncertainty__1down:       42790.4 ID 0.992328 ISO 1.02283  TRIG 1.04157
////PH_EFF_ID_Uncertainty__1up:         42790.4 ID 1.01816  ISO 1.02283  TRIG 1.04157
////PH_EFF_ISO_Uncertainty__1down:      42790.4 ID 1.00524  ISO 1.00971  TRIG 1.04157
////PH_EFF_ISO_Uncertainty__1up:        42790.4 ID 1.00524  ISO 1.03595  TRIG 1.04157
////PH_EFF_TRIGGER_Uncertainty__1down:  42790.4 ID 1.00524  ISO 1.02283  TRIG 1.01646
////PH_EFF_TRIGGER_Uncertainty__1up:    42790.4 ID 1.00524  ISO 1.02283  TRIG 1.06668
//
//    ROOT::RDF::RunGraphs(Nodes); // running all computation nodes concurrently
//
////    resultmapHists
////    ==============
////    0   //ma1
////    1   //ma2
////    2   //ma3
////    3   //ma5
////    4   //ma9
////    5   //data
////    6   //Z-gamma
////    7   //Z-gamma
////    8   //Z-gamma
////    9   //Z-jets
////    10  //Z-jets
////    11  //Z-jets
////    12  //Z-jets
////    13  //Z-jets
////    14  //Z-jets
////    15  //Z-jets
////    16  //Z-jets
////    17  //Z-jets
////
////    resultmaps
////    ==========
////
////    0       //Z-gamma
////    1       //Z-gamma
////    2       //Z-gamma
////    3       //Z-jets
////    4       //Z-jets
////    5       //Z-jets
////    6       //Z-jets
////    7       //Z-jets
////    8       //Z-jets
////    9       //Z-jets
////    10      //Z-jets
////    11      //Z-jets
//
//    constexpr std::array<const char*, 3> systematics = {"photons_and_electrons:EG_RESOLUTION_ALL__1down", "photons_and_electrons:EG_SCALE_ALL__1down",
//        "photon_id_eff:PH_EFF_ID_Uncertainty__1down"};
//    constexpr std::array<const char*, 4> nominal_systematics =
//    {
//        "nominal",
//        "photons_and_electrons:EG_RESOLUTION_ALL__1down",
//        "photons_and_electrons:EG_SCALE_ALL__1down",
//        "photon_id_eff:PH_EFF_ID_Uncertainty__1down",
//    };
//    constexpr std::array<const char*, 3> Systematics = {"EG_RESOLUTION_ALL__1down", "EG_SCALE_ALL__1down", "PH_EFF_ID_Uncertainty__1down"};
//    constexpr std::array<const char*, 4> nominalSystematics = {"nominal", "EG_RESOLUTION_ALL__1down", "EG_SCALE_ALL__1down", "PH_EFF_ID_Uncertainty__1down"};
//    constexpr std::array<const char*, 6> signalPlusDataSamples = {"Sig m_{A} = 1 GeV", "Sig m_{A} = 2 GeV", "Sig m_{A} = 3 GeV", "Sig m_{A} = 5 GeV", "Sig m_{A} = 9 GeV", "Data"};
//    constexpr std::array<const char*, 6> signalPlusDataFileNames = {"prompt_ma1_merged_higgs_syst_down.pdf", "prompt_ma2_merged_higgs_syst_down.pdf", "prompt_ma3_merged_higgs_syst_down.pdf", "prompt_ma5_merged_higgs_syst_down.pdf", "prompt_ma9_merged_higgs_syst_down.pdf", "data_merged_higgs_syst_down.pdf"};
//    constexpr std::array<EColor, 3> colors = {kBlue, kRed, static_cast<EColor>(kOrange+1)};
//    constexpr std::array<EColor, 4> nominalColors = {kBlack, kBlue, kRed, static_cast<EColor>(kOrange+1)};
//    constexpr std::array<float, 6> signalPlusDataResidualMaxima = {2.1, 2.35, 1.05, 0.45, 1.2, 1.05};
//    constexpr std::array<float, 6> signalPlusDataResidualMinima = {-1.2, -1.2, -1.1, -1.1, -1.1, 0};
//
//    constexpr std::array<std::array<float, 4>, 6> signalPlusDataResidualLegend =
//    {{
//        {0.13, 0.76, 0.33, 0.96},
//        {0.5, 0.7, 0.7, 0.9},
//        {0.6, 0.7, 0.8, 0.9},
//        {0.6, 0.35, 0.8, 0.55},
//        {0.51, 0.7, 0.71, 0.9},
//        {0.5, 0.7, 0.7, 0.9},
//    }};
//
//    TCanvas* c1;
//    TLegend* legend;
//    TLatex Tl;
//
//    //signal+data
//    for (int i = 0; i < 6; i++) //signal samples plus data
//    {
//        c1 = new TCanvas();
//        legend = new TLegend(0.5, 0.2, 0.85, 0.6);
//        c1->cd();
//        TPad* pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 0.95);
//        pad1->SetTopMargin(0.1);
//        pad1->SetBottomMargin(0);
//        pad1->Draw();
//        pad1->cd();
//
//        gStyle->SetOptStat(0);
//        resultmapHists[i]["nominal"].SetLineColor(kBlack);
//        resultmapHists[i]["nominal"].SetTitle(";m_{ll#gamma} [GeV];Events");
//        resultmapHists[i]["nominal"].SetTitle(signalPlusDataSamples[i]);
//        resultmapHists[i]["nominal"].GetXaxis()->SetTitle("m_{ll#gamma} [GeV]");
//        resultmapHists[i]["nominal"].GetXaxis()->SetTitleOffset(1.2);
//        resultmapHists[i]["nominal"].GetYaxis()->SetTitle("Events");
//        resultmapHists[i]["nominal"].GetYaxis()->CenterTitle(true);
//        legend->AddEntry(&resultmapHists[i]["nominal"], "nominal", "l");
//
//        resultmapHists[i]["nominal"].Draw("HIST");
//
//        for (int j = 0; j < 3; j++)
//        {
//            std::cout << "difference = "
//            << resultmapHists[i][systematics[j]].Integral() - resultmapHists[i]["nominal"].Integral() << '\n';
//            resultmapHists[i][systematics[j]].SetLineColor(colors[j]);
//            legend->AddEntry(&resultmapHists[i][systematics[j]], Systematics[j], "l");
//            resultmapHists[i][systematics[j]].Draw("HISTsame");
//        }
//
//        Tl.SetTextSize(0.03);
//        Tl.DrawLatexNDC(0.6, 0.8, "#it{ATLAS} Internal");
//        Tl.DrawLatexNDC(0.6, 0.7,"#sqrt{s} = 13 TeV");
//        legend->SetBorderSize(0);
//        legend->Draw("same");
//
//        // Create and draw the error plot
//        TH1F* errorHist;
//        double error;
//
//        c1->cd();
//        TPad* pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.25);
//        pad2->SetTopMargin(0);
//        pad2->SetBottomMargin(0.25);
//        pad2->Draw();
//        pad2->cd();
//
//        legend = new TLegend(signalPlusDataResidualLegend[i][0], signalPlusDataResidualLegend[i][1], signalPlusDataResidualLegend[i][2], signalPlusDataResidualLegend[i][3]);
//
//        for (int j = 0; j < 3; j++)
//        {
//            errorHist = new TH1F(Form("errorHist_%d%d", i, j), "", resultmapHists[i]["nominal"].GetNbinsX(),
//                                       resultmapHists[i]["nominal"].GetXaxis()->GetXmin(),
//                                       resultmapHists[i]["nominal"].GetXaxis()->GetXmax());
//            errorHist->SetLineColor(colors[j]);
//            errorHist->SetLineWidth(1);
//            errorHist->SetMarkerStyle(20);
//            errorHist->SetMarkerSize(0.9);
//            errorHist->GetXaxis()->SetLabelSize(0.1);
//            errorHist->GetYaxis()->SetLabelSize(0.1);
//            errorHist->GetXaxis()->SetTitleSize(0.1);
//            errorHist->GetYaxis()->SetTitleSize(0.1);
//            errorHist->SetMarkerColor(colors[j]);
//
//            for (int bin = 1; bin <= resultmapHists[i]["nominal"].GetNbinsX(); bin++)
//            {
//                double nominal = resultmapHists[i]["nominal"].GetBinContent(bin);
//                if (!nominal)
//                {
//                    error = 0;
//                }
//                else
//                {
//                    double systematic = resultmapHists[i][systematics[j]].GetBinContent(bin);
//                    error = (systematic - nominal) / nominal;
//                }
//                errorHist->SetBinContent(bin, error);
//            }
//
//            if (j == 0)
//            {
//                errorHist->SetTitle(";m_{ll#gamma} [GeV];#frac{syst-nominal}{nominal}");
//                errorHist->GetXaxis()->SetTitle("m_{ll#gamma} [GeV]");
//                errorHist->GetXaxis()->SetTitleOffset(1.2);
//                errorHist->GetYaxis()->SetTitle("#frac{syst-nominal}{nominal}");
//                errorHist->GetYaxis()->SetTitleOffset(0.4);
//                errorHist->GetYaxis()->CenterTitle(true);
//                errorHist->SetMaximum(signalPlusDataResidualMaxima[i]);
//                errorHist->SetMinimum(signalPlusDataResidualMinima[i]);
//                errorHist->Draw("histsame");
//            }
//            else
//            {
//                c1->cd();
//                pad2->cd();
//                errorHist->Draw("histsame");
//            }
//
//            legend->AddEntry(errorHist, Systematics[j], "l");
//
//            c1->cd();
//            c1->Update();
//
//        }
//        pad2->cd();
//        legend->SetBorderSize(0);
//        legend->Draw("same");
//        c1->SaveAs(signalPlusDataFileNames[i]);
//    }
//
//    constexpr std::array<double,3> SFs = {((139e15)*(.871e-12)),((139e15)*(.199e-12)), ((139e15)*(.0345e-15))}; //numerators for Z-gamma bkg
//
//
//    //Z-gamma
//    c1 = new TCanvas();
//    legend = new TLegend(0.575, 0.4, 0.875, 0.7);
//    TH1D *hist, *nominalhist;
//    gStyle->SetOptStat(0);
//    c1->cd();
//    TPad* pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 0.95);
//    pad1->SetTopMargin(0.1);
//    pad1->SetBottomMargin(0);
//    pad1->Draw();
//    pad1->cd();
//    gStyle->SetOptStat(0);
//    c1->cd();
//    TPad* pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.25);
//    pad2->SetTopMargin(0);
//    pad2->SetBottomMargin(0.25);
//    pad2->Draw();
//    pad2->cd();
//    TLegend* errorLegend = new TLegend(0.325, 0.75, 0.525, 0.95);
//    TH1F* errorHist;
//    double error;
//
//    for (int i = 0; i < 4; i++) //systematics: Z-gamma
//    {
//        pad1->cd();
//        hist = new TH1D(Form("Z#gamma#gamma%d",i),Form("Z#gamma#gamma%d",i), 100u, 80, 200);
//        for (int j = 0, k = 6; k <= 8; j++, k++)
//        {
//            hist->Add(&resultmapHists[k][nominal_systematics[i]], SFs[j] / *Nodes[j].GetResultPtr<float>());
//        }
//        hist->SetLineColor(nominalColors[i]);
//        std::cout <<  hist->Integral() << '\n';
//
//        if (i == 0)
//        {
//            nominalhist = static_cast<TH1D*>(hist->Clone());
//            hist->SetTitle(";m_{ll#gamma} [GeV];Events");
//            hist->SetTitle("Z#gamma#gamma");
//            hist->GetXaxis()->SetTitle("m_{ll#gamma} [GeV]");
//            hist->GetXaxis()->SetTitleOffset(1.2);
//            hist->GetYaxis()->SetTitle("Events");
//            hist->GetYaxis()->CenterTitle(true);
//            hist->SetMaximum(33);
//            hist->Draw("HIST");
//        }
//        else
//        {
//            hist->Draw("HISTsame");
//
//            //Now do the error plot
//
//            c1->cd();
//            pad2->cd();
//
//
//
//            errorHist = new TH1F(Form("errorHist_%d", i), "", nominalhist->GetNbinsX(),
//                                 nominalhist->GetXaxis()->GetXmin(),
//                                 nominalhist->GetXaxis()->GetXmax());
//            errorHist->SetLineColor(colors[i-1]);
//            errorHist->SetLineWidth(1);
//            errorHist->SetMarkerStyle(20);
//            errorHist->SetMarkerSize(0.9);
//            errorHist->GetXaxis()->SetLabelSize(0.1);
//            errorHist->GetYaxis()->SetLabelSize(0.1);
//            errorHist->GetXaxis()->SetTitleSize(0.1);
//            errorHist->GetYaxis()->SetTitleSize(0.1);
//            errorHist->SetMarkerColor(nominalColors[i]);
//
//            for (int bin = 1; bin <= nominalhist->GetNbinsX(); bin++)
//            {
//                double nominal = nominalhist->GetBinContent(bin);
//                if (!nominal)
//                {
//                    error = 0;
//                }
//                else
//                {
//                    double systematic = hist->GetBinContent(bin);
//                    error = (systematic - nominal) / nominal;
//                }
//                errorHist->SetBinContent(bin, error);
//            }
//
//            if (i == 1)
//            {
//                errorHist->SetTitle(";m_{ll#gamma} [GeV];#frac{syst-nominal}{nominal}");
//                errorHist->GetXaxis()->SetTitle("m_{ll#gamma} [GeV]");
//                errorHist->GetXaxis()->SetTitleOffset(1.2);
//                errorHist->GetYaxis()->SetTitle("#frac{syst-nominal}{nominal}");
//                errorHist->GetYaxis()->SetTitleOffset(0.4);
//                errorHist->GetYaxis()->CenterTitle(true);
//                errorHist->SetMaximum(1.25);
//                errorHist->SetMinimum(-1.1);
//                errorHist->Draw("histsame");
//
//            }
//            else
//            {
//                c1->cd();
//                pad2->cd();
//                errorHist->Draw("histsame");
//            }
//
//            errorLegend->AddEntry(errorHist, nominalSystematics[i], "l");
//
//            c1->cd();
//            c1->Update();
//            pad2->cd();
//            errorLegend->SetBorderSize(0);
//            errorLegend->Draw("same");
//
//        }
//        legend->AddEntry(hist, nominalSystematics[i], "l");
//    }
//    pad1->cd();
//    Tl.SetTextSize(0.03);
//    Tl.DrawLatexNDC(0.6, 0.85, "#it{ATLAS} Internal");
//    Tl.DrawLatexNDC(0.6, 0.775,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
//    legend->SetBorderSize(0);
//    legend->Draw("same");
//
//    c1->SaveAs("Zgamma_merged_higgs_syst_down.pdf");
//
//    //Z-jets
//    std::array<double,9> JetNumeratorSFs = {((139e15)*(1.9828e-9)*(0.821204)),((139e15)*(110.64e-12)*(0.69275)),((139e15)*(40.645e-12)*(0.615906)),((139e15)*(1.9817e-9)*(0.1136684)),((139e15)*(110.47e-12)*(0.1912956)),((139e15)*(40.674e-12)*(0.2326772)),((139e15)*(1.9819e-9)*(0.0656969)),((139e15)*(110.53e-12)*(0.1158741)),((139e15)*(40.68e-12)*(0.1535215))}; //numerators for jet bkg
//
//    c1 = new TCanvas();
//    legend = new TLegend(0.51, 0.25, 0.85, 0.65);
//    gStyle->SetOptStat(0);
//    c1->cd();
//    pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 0.95);
//    pad1->SetTopMargin(0.1);
//    pad1->SetBottomMargin(0);
//    pad1->Draw();
//    pad1->cd();
//    gStyle->SetOptStat(0);
//    c1->cd();
//    pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.25);
//    pad2->SetTopMargin(0);
//    pad2->SetBottomMargin(0.25);
//    pad2->Draw();
//    pad2->cd();
//    errorLegend = new TLegend(0.3, 0.7, 0.5, 0.9);
//
//    for (int i = 0; i < 4; i++) //systematics: Z-jets
//    {
//        pad1->cd();
//        hist = new TH1D(Form("Z+Jets%d",i),Form("Z+Jets%d",i), 100u, 80, 200);
//        for (int j = 0, k = 9; k <= 17; j++, k++)
//        {
//            hist->Add(&resultmapHists[k][nominal_systematics[i]], JetNumeratorSFs[j] / *Nodes[j+3].GetResultPtr<float>());
//        }
//        hist->SetLineColor(nominalColors[i]);
//        std::cout <<  hist->Integral() << '\n';
//
//        if (i == 0)
//        {
//            nominalhist = static_cast<TH1D*>(hist->Clone());
//            hist->SetTitle(";m_{ll#gamma} [GeV];Events");
//            hist->SetTitle("Z+Jets");
//            hist->GetXaxis()->SetTitle("m_{ll#gamma} [GeV]");
//            hist->GetXaxis()->SetTitleOffset(1.2);
//            hist->GetYaxis()->SetTitle("Events");
//            hist->GetYaxis()->CenterTitle(true);
//            hist->Draw("HIST");
//        }
//        else
//        {
//
//            hist->Draw("HISTsame");
//            //Now do the error plot
//
//            c1->cd();
//            pad2->cd();
//            errorHist = new TH1F(Form("errorHist_%d", i), "", nominalhist->GetNbinsX(),
//                                 nominalhist->GetXaxis()->GetXmin(),
//                                 nominalhist->GetXaxis()->GetXmax());
//            errorHist->SetLineColor(colors[i-1]);
//            errorHist->SetLineWidth(1);
//            errorHist->SetMarkerStyle(20);
//            errorHist->SetMarkerSize(0.9);
//            errorHist->GetXaxis()->SetLabelSize(0.1);
//            errorHist->GetYaxis()->SetLabelSize(0.1);
//            errorHist->GetXaxis()->SetTitleSize(0.1);
//            errorHist->GetYaxis()->SetTitleSize(0.1);
//            errorHist->SetMarkerColor(nominalColors[i]);
//
//            for (int bin = 1; bin <= nominalhist->GetNbinsX(); bin++)
//            {
//                double nominal = nominalhist->GetBinContent(bin);
//                if (!nominal)
//                {
//                    error = 0;
//                }
//                else
//                {
//                    double systematic = hist->GetBinContent(bin);
//                    error = (systematic - nominal) / nominal;
//                }
//                errorHist->SetBinContent(bin, error);
//            }
//
//            if (i == 1)
//            {
//                errorHist->SetTitle(";m_{ll#gamma} [GeV];#frac{syst-nominal}{nominal}");
//                errorHist->GetXaxis()->SetTitle("m_{ll#gamma} [GeV]");
//                errorHist->GetXaxis()->SetTitleOffset(1.2);
//                errorHist->GetYaxis()->SetTitle("#frac{syst-nominal}{nominal}");
//                errorHist->GetYaxis()->SetTitleOffset(0.4);
//                errorHist->GetYaxis()->CenterTitle(true);
//                errorHist->SetMaximum(1.4);
//                errorHist->SetMinimum(-0.67);
//                errorHist->Draw("histsame");
//
//            }
//            else
//            {
//                c1->cd();
//                pad2->cd();
//                errorHist->Draw("histsame");
//            }
//
//            errorLegend->AddEntry(errorHist, nominalSystematics[i], "l");
//
//            c1->cd();
//            c1->Update();
//            pad2->cd();
//            errorLegend->SetBorderSize(0);
//            errorLegend->Draw("same");
//        }
//        legend->AddEntry(hist, nominalSystematics[i], "l");
//    }
//    pad1->cd();
//    Tl.SetTextSize(0.03);
//    Tl.DrawLatexNDC(0.6, 0.85, "#it{ATLAS} Internal");
//    Tl.DrawLatexNDC(0.6, 0.775,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
//    legend->SetBorderSize(0);
//    legend->Draw("same");
//    c1->SaveAs("Zjets_merged_higgs_syst_down.pdf");
//}
//
//void Resolved_Up_Higgs_Mass()
//{
//    Event::systematics =
//    {
////        "PH_EFF_ISO_Uncertainty",
////        "PH_EFF_ISO_Uncertainty",
////        "EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR",
////        "PRW_DATASF",
////        "MUON_EFF_RECO_SYS",
////        "MUON_EFF_ISO_SYS",
////        "MUON_EFF_TrigSystUncertainty",
////        "EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR",
////        "MUON_EFF_TrigStatUncertainty",
////        "MUON_EFF_RECO_STAT",
////        "MUON_EFF_TTVA_STAT",
////        "EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR",
////        "EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR",
////        "MUON_EFF_TTVA_SYS",
////        "MUON_EFF_ISO_STAT",
////        "MUON_SAGITTA_RHO",
////        "EG_RESOLUTION_ALL",
////        "EG_SCALE_ALL",
////        "MUON_MS",
////        "MUON_ID",
////        "EL_EFF_TriggerEff_TOTAL_1NPCOR_PLUS_UNCOR",
////        "MUON_SAGITTA_RESBIAS",
////        "MUON_SCALE",
//
////        "PH_EFF_ISO_Uncertainty__1down",
////        "EG_SCALE_ALL__1down",
////        "PH_EFF_ID_Uncertainty__1down",
////        "PH_EFF_TRIGGER_Uncertainty__1down",
//        "EG_RESOLUTION_ALL__1up",
//        "EG_SCALE_ALL__1up",
////        "PH_EFF_ISO_Uncertainty__1up",
//        "PH_EFF_ID_Uncertainty__1up",
////        "PH_EFF_TRIGGER_Uncertainty__1up",
////        "EG_RESOLUTION_ALL__1down",
//    };
//
//    std::vector<std::vector<std::string>> input_filenames = {
//        //Signal
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root"}, //1 GeV
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600751.PhPy8EG_AZNLO_ggH125_mA2p0_v1.root"}, // 2 GeV
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600752.PhPy8EG_AZNLO_ggH125_mA3p0_v1.root"}, // 3 GeV
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root"}, //5 GeV
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600756.PhPy8EG_AZNLO_ggH125_mA9p0_v1.root"}, // 9 GeV
//        //Data
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_data_test.root"},
//        //Z gamma background
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617070._000001.LGNTuple.root"},
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617064._000001.LGNTuple.root"},
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617074._000001.LGNTuple.root"},
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
//    std::vector<const char*> prefixes = {"Sig m_{A} = 1 GeV", "Sig m_{A} = 2 GeV", "Sig m_{A} = 3 GeV", "Sig m_{A} = 5 GeV", "Sig m_{A} = 9 GeV", "data", "pty2_9_17", "pty_17_myy_0_80", "pty_17_myy_80", "Zee_lightJet_0-70", "Zee_lightJet_70-140", "Zee_lightJet_140-280", "Zee_cJet_0-70", "Zee_cJet_70-140", "Zee_cJet_140-280", "Zee_bJet_0-70", "Zee_bJet_70-140", "Zee_bJet_140-280"};
//
//    std::vector<ROOT::RDF::RResultHandle> Nodes;
//
//    int count = 0;
//
//    std::vector<RResultMap<float>> resultmaps;
//    std::vector<RResultMap<TH1D>> resultmapHists;
//
//    auto verbosity = ROOT::Experimental::RLogScopedVerbosity(ROOT::Detail::RDF::RDFLogChannel(), ROOT::Experimental::ELogLevel::kInfo);
//
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
//        [](RVec<float> photon_id_eff, RVec<float> photon_iso_eff, RVec<float> photon_trg_eff/*, RVec<float> ei_event_weights_generator*/)
//        {
//            auto ResizeVal = std::max({photon_id_eff.size(), photon_iso_eff.size(), photon_trg_eff.size()});
//            photon_id_eff.resize(ResizeVal,1);
//            photon_iso_eff.resize(ResizeVal,1);
//            photon_trg_eff.resize(ResizeVal,1);
//
//            return photon_id_eff;
//            // *photon_iso_eff*photon_trg_eff; //Returns a vector of efficiencies
//
//        }, {"photon_id_eff", "photon_iso_eff", "photon_trg_eff",/* "ei_event_weights_generator"*/});
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
//        [](RVec<AbstractParticle> electrons)
//        {
//            //keep the electrons in each event that have pt > 20 GeV, |η| < 2.37,
//            //|η| not between 1.37 and 1.52, and that satisfy a medium id criteria `electron_id_medium`
//            electrons.erase(std::remove_if(electrons.begin(),electrons.end(),
//            [](AbstractParticle& ep)
//            {
//                 return (!((ep.electron_pt/1e3 > 20) && (std::abs(ep.electron_eta) < 2.37) &&
//                 (!((1.37 < std::abs(ep.electron_eta)) && (std::abs(ep.electron_eta) < 1.52)))
//                 && (ep.electron_id_medium == 1)));
//
//            }), electrons.end());
//
//            return electrons;
//
//        },{"abstract_electrons"}).Filter(
//        [](RVec<Muon>& muons, RVec<AbstractParticle>& electrons)
//        {
//            return (electrons.size()==2 && muons.empty()); //keep events which have exactly 2 electrons for di_electrons and no muons
//
//        }, {"muons", "di_electrons"});
//
//        //new dataframe node: contains only the events from `two_leptons` whose electrons in the `di_electrons` branch have opposite charge
//        auto opp_charge = two_leptons
//        .Filter([](RVec<AbstractParticle>& electrons)
//        {
//            return (electrons[0].electron_charge*electrons[1].electron_charge < 0);
//
//        }, {"di_electrons"});
//
//        //new dataframe node: contains only the events from `opp_charge` that have 1 electron with pt > 20 GeV and the other with pt > 27 GeV
//        auto leadingPt = opp_charge.Filter([](RVec<AbstractParticle>& electrons)
//        {
//            return ((electrons[0].electron_pt > 20e3 && electrons[1].electron_pt > 27e3) || (electrons[1].electron_pt > 20e3 && electrons[0].electron_pt > 27e3));
//        }, {"di_electrons"});
//
//        auto deltaR = leadingPt.Filter([] (RVec<AbstractParticle>& electrons)
//        {
//            return (DeltaR(electrons[0].ElectronVector(), electrons[1].ElectronVector()) > 0.01);
//        }, {"di_electrons"});
//
//        //new dataframe node: Contains a new column `dilep` in addition to the ones in `same_flavour` that stores the di-electron four-vector
//        auto dilep = deltaR.Define("dilep",[] (RVec<AbstractParticle>& electrons)
//        {
//            return (electrons[0].ElectronVector() + electrons[1].ElectronVector());
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
//        auto photon_passes_cuts = ptCut
//        .Define("photons_pass_cut_indices",
//        [&](const RVec<AbstractParticle>& photons)
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
//        }, {"abstract_photons"})
//        .Define("photons_pass_cuts",
//        [&](RVec<AbstractParticle>& photons, RVec<int>& photon_indices)
//        {
//            return Take(photons, photon_indices);
//        }, {"abstract_photons", "photons_pass_cut_indices"});
//
//        //New dataframe node: contains only events from `photon_passes_cuts` that pass the resolved category
//        auto resolved = photon_passes_cuts.Define("chosen_two_indices",
//        [](RVec<AbstractParticle>& photons_pass_cuts)
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
//                delta_r = DeltaR(photons_pass_cuts[combs[0][i]].PhotonVector(), photons_pass_cuts[combs[1][i]].PhotonVector());
//                m = (photons_pass_cuts[combs[0][i]].PhotonVector() + photons_pass_cuts[combs[1][i]].PhotonVector()).M();
//                pt = (photons_pass_cuts[combs[0][i]].PhotonVector() + photons_pass_cuts[combs[1][i]].PhotonVector()).Pt();
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
//        [&](RVec<AbstractParticle>& reco_photons_matched, RVec<unsigned long>& indices)
//        {
//            return Take(reco_photons_matched, indices);
//
//        }, {"photons_pass_cuts", "chosen_two_indices"})
//        .Define("diphoton",
//        [&](RVec<AbstractParticle>& reco_photons_matched)
//        {
//            return reco_photons_matched[0].PhotonVector() + reco_photons_matched[1].PhotonVector();
//        }, {"chosen_two"})
//        .Define("reco_higgs_mass",
//        [&](PtEtaPhiEVector& diphoton, PtEtaPhiEVector& dilep)
//        {
//            return (dilep + diphoton).M() / 1e3;
//        }, {"diphoton", "dilep"})
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
//        if (count >= 6)
//        {
//            Nodes.push_back(EventWeight.Sum<float>("EventWeight"));
//        }
//        resultmapHists.push_back(VariationsFor(resolved.Histo1D<double>({prefixes[count], prefixes[count++], 100u, 80, 200}, "reco_higgs_mass", "totEventWeight")));
//
//    }
//
////EG_RESOLUTION_ALL__1down:           42864.7 ID 1.00524  ISO 1.02283  TRIG 1.04157
////EG_RESOLUTION_ALL__1up:             42751.5 ID 1.00524  ISO 1.02283  TRIG 1.04157
////EG_SCALE_ALL__1down:                42521.1 ID 1.00524  ISO 1.02283  TRIG 1.04157
////EG_SCALE_ALL__1up:                  43059.8 ID 1.00524  ISO 1.02283  TRIG 1.04157
////PH_EFF_ID_Uncertainty__1down:       42790.4 ID 0.992328 ISO 1.02283  TRIG 1.04157
////PH_EFF_ID_Uncertainty__1up:         42790.4 ID 1.01816  ISO 1.02283  TRIG 1.04157
////PH_EFF_ISO_Uncertainty__1down:      42790.4 ID 1.00524  ISO 1.00971  TRIG 1.04157
////PH_EFF_ISO_Uncertainty__1up:        42790.4 ID 1.00524  ISO 1.03595  TRIG 1.04157
////PH_EFF_TRIGGER_Uncertainty__1down:  42790.4 ID 1.00524  ISO 1.02283  TRIG 1.01646
////PH_EFF_TRIGGER_Uncertainty__1up:    42790.4 ID 1.00524  ISO 1.02283  TRIG 1.06668
//
//    ROOT::RDF::RunGraphs(Nodes); // running all computation nodes concurrently
//
////    resultmapHists
////    ==============
////    0   //ma1
////    1   //ma2
////    2   //ma3
////    3   //ma5
////    4   //ma9
////    5   //data
////    6   //Z-gamma
////    7   //Z-gamma
////    8   //Z-gamma
////    9   //Z-jets
////    10  //Z-jets
////    11  //Z-jets
////    12  //Z-jets
////    13  //Z-jets
////    14  //Z-jets
////    15  //Z-jets
////    16  //Z-jets
////    17  //Z-jets
////
////    resultmaps
////    ==========
////
////    0       //Z-gamma
////    1       //Z-gamma
////    2       //Z-gamma
////    3       //Z-jets
////    4       //Z-jets
////    5       //Z-jets
////    6       //Z-jets
////    7       //Z-jets
////    8       //Z-jets
////    9       //Z-jets
////    10      //Z-jets
////    11      //Z-jets
//
//    constexpr std::array<const char*, 3> systematics =
//    {
//        "photons_and_electrons:EG_RESOLUTION_ALL__1up",
//        "photons_and_electrons:EG_SCALE_ALL__1up",
//        "photon_id_eff:PH_EFF_ID_Uncertainty__1up"
//    };
//    constexpr std::array<const char*, 4> nominal_systematics =
//    {
//        "nominal",
//        "photons_and_electrons:EG_RESOLUTION_ALL__1up",
//        "photons_and_electrons:EG_SCALE_ALL__1up",
//        "photon_id_eff:PH_EFF_ID_Uncertainty__1up",
//    };
//    constexpr std::array<const char*, 3> Systematics = {"EG_RESOLUTION_ALL__1up", "EG_SCALE_ALL__1up", "PH_EFF_ID_Uncertainty__1up"};
//    constexpr std::array<const char*, 4> nominalSystematics = {"nominal", "EG_RESOLUTION_ALL__1up", "EG_SCALE_ALL__1up", "PH_EFF_ID_Uncertainty__1up"};
//    constexpr std::array<const char*, 6> signalPlusDataSamples = {"Sig m_{A} = 1 GeV", "Sig m_{A} = 2 GeV", "Sig m_{A} = 3 GeV", "Sig m_{A} = 5 GeV", "Sig m_{A} = 9 GeV", "Data"};
//    constexpr std::array<const char*, 6> signalPlusDataFileNames = {"prompt_ma1_resolved_higgs_syst_up.pdf", "prompt_ma2_resolved_higgs_syst_up.pdf", "prompt_ma3_resolved_higgs_syst_up.pdf", "prompt_ma5_resolved_higgs_syst_up.pdf", "prompt_ma9_resolved_higgs_syst_up.pdf", "data_resolved_higgs_syst_up.pdf"};
//    constexpr std::array<EColor, 3> colors = {kBlue, kRed, static_cast<EColor>(kOrange+1)};
//    constexpr std::array<EColor, 4> nominalColors = {kBlack, kBlue, kRed, static_cast<EColor>(kOrange+1)};
//    constexpr std::array<int, 5> maxima = {62, 80, 270, 155, 155};
//    constexpr std::array<float, 6> signalPlusDataResidualMaxima = {1.2, 1.1, 2.15, 0.55, 1.2, 1.05};
//    constexpr std::array<float, 6> signalPlusDataResidualMinima = {-1.2, -1.2, -1.1, -1.15, -1.2, 0};
//    constexpr std::array<std::array<float, 4>, 6> signalPlusDataResidualLegend =
//    {{
//        {0.12, 0.74, 0.32, 0.94},
//        {0.12, 0.74, 0.32, 0.94},
//        {0.67, 0.74, 0.87, 0.94},
//        {0.12, 0.34, 0.32, 0.54},
//        {0.12, 0.74, 0.32, 0.94},
//        {0.15, 0.74, 0.35, 0.94},
//    }};
//
//    TCanvas* c1;
//    TLegend* legend;
//    TLatex Tl;
//
//    //signal+data
//    for (int i = 0; i < 6; i++) //signal samples plus data
//    {
//        c1 = new TCanvas();
//        if (i == 5)
//        {
//            legend = new TLegend(0.15, 0.35, 0.4, 0.7);
//        }
//        else
//        {
//            legend = new TLegend(0.5, 0.25, 0.85, 0.65);
//            resultmapHists[i]["nominal"].SetMaximum(maxima[i]);
//        }
//        c1->cd();
//        TPad* pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 0.95);
//        pad1->SetTopMargin(0.1);
//        pad1->SetBottomMargin(0);
//        pad1->Draw();
//        pad1->cd();
//
//        gStyle->SetOptStat(0);
//        resultmapHists[i]["nominal"].SetLineColor(kBlack);
//        resultmapHists[i]["nominal"].SetTitle(";m_{ll#gamma#gamma} [GeV];Events");
//        resultmapHists[i]["nominal"].SetTitle(signalPlusDataSamples[i]);
//        resultmapHists[i]["nominal"].GetXaxis()->SetTitle("m_{ll#gamma#gamma} [GeV]");
//        resultmapHists[i]["nominal"].GetXaxis()->SetTitleOffset(1.2);
//        resultmapHists[i]["nominal"].GetYaxis()->SetTitle("Events");
//        resultmapHists[i]["nominal"].GetYaxis()->CenterTitle(true);
//
//        legend->AddEntry(&resultmapHists[i]["nominal"], "nominal", "l");
//
//        resultmapHists[i]["nominal"].Draw("HIST");
//
//        for (int j = 0; j < 3; j++)
//        {
//            std::cout << "difference = "
//            << resultmapHists[i][systematics[j]].Integral() - resultmapHists[i]["nominal"].Integral() << '\n';
//            resultmapHists[i][systematics[j]].SetLineColor(colors[j]);
//            legend->AddEntry(&resultmapHists[i][systematics[j]], Systematics[j], "l");
//            resultmapHists[i][systematics[j]].Draw("HISTsame");
//        }
//
//        Tl.SetTextSize(0.03);
//        Tl.DrawLatexNDC(0.6, 0.8, "#it{ATLAS} Internal");
//        Tl.DrawLatexNDC(0.6, 0.7,"#sqrt{s} = 13 TeV");
//        legend->SetBorderSize(0);
//        legend->Draw("same");
//
//        // Create and draw the error plot
//        TH1F* errorHist;
//        double error;
//
//        c1->cd();
//        TPad* pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.25);
//        pad2->SetTopMargin(0);
//        pad2->SetBottomMargin(0.25);
//        pad2->Draw();
//        pad2->cd();
//        legend = new TLegend(signalPlusDataResidualLegend[i][0], signalPlusDataResidualLegend[i][1], signalPlusDataResidualLegend[i][2], signalPlusDataResidualLegend[i][3]);
//
//        for (int j = 0; j < 3; j++)
//        {
//            errorHist = new TH1F(Form("errorHist_%d%d", i, j), "", resultmapHists[i]["nominal"].GetNbinsX(),
//                                       resultmapHists[i]["nominal"].GetXaxis()->GetXmin(),
//                                       resultmapHists[i]["nominal"].GetXaxis()->GetXmax());
//            errorHist->SetLineColor(colors[j]);
//            errorHist->SetLineWidth(1);
//            errorHist->SetMarkerStyle(20);
//            errorHist->SetMarkerSize(0.9);
//            errorHist->GetXaxis()->SetLabelSize(0.1);
//            errorHist->GetYaxis()->SetLabelSize(0.1);
//            errorHist->GetXaxis()->SetTitleSize(0.1);
//            errorHist->GetYaxis()->SetTitleSize(0.1);
//            errorHist->SetMarkerColor(colors[j]);
//
//            for (int bin = 1; bin <= resultmapHists[i]["nominal"].GetNbinsX(); bin++)
//            {
//                double nominal = resultmapHists[i]["nominal"].GetBinContent(bin);
//                if (!nominal)
//                {
//                    error = 0;
//                }
//                else
//                {
//                    double systematic = resultmapHists[i][systematics[j]].GetBinContent(bin);
//                    error = (systematic - nominal) / nominal;
//                }
//                errorHist->SetBinContent(bin, error);
//            }
//
//            if (j == 0)
//            {
//                errorHist->SetTitle(";m_{ll#gamma#gamma} [GeV];#frac{syst-nominal}{nominal}");
//                errorHist->GetXaxis()->SetTitle("m_{ll#gamma#gamma} [GeV]");
//                errorHist->GetXaxis()->SetTitleOffset(1.2);
//                errorHist->GetYaxis()->SetTitle("#frac{syst-nominal}{nominal}");
//                errorHist->GetYaxis()->SetTitleOffset(0.4);
//                errorHist->GetYaxis()->CenterTitle(true);
//                errorHist->SetMaximum(signalPlusDataResidualMaxima[i]);
//                errorHist->SetMinimum(signalPlusDataResidualMinima[i]);
//                errorHist->Draw("histsame");
//            }
//            else
//            {
//                c1->cd();
//                pad2->cd();
//                errorHist->Draw("histsame");
//            }
//
//            legend->AddEntry(errorHist, Systematics[j], "l");
//
//            c1->cd();
//            c1->Update();
//
//        }
//        pad2->cd();
//        legend->SetBorderSize(0);
//        legend->Draw("same");
//        c1->SaveAs(signalPlusDataFileNames[i]);
//    }
//
//    constexpr std::array<double,3> SFs = {((139e15)*(.871e-12)),((139e15)*(.199e-12)), ((139e15)*(.0345e-15))}; //numerators for Z-gamma bkg
//
//
//    //Z-gamma
//    c1 = new TCanvas();
//    legend = new TLegend(0.15, 0.57, 0.42, 0.87);
//    TH1D *hist, *nominalhist;
//    gStyle->SetOptStat(0);
//    c1->cd();
//    TPad* pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 0.95);
//    pad1->SetTopMargin(0.1);
//    pad1->SetBottomMargin(0);
//    pad1->Draw();
//    pad1->cd();
//    gStyle->SetOptStat(0);
//    c1->cd();
//    TPad* pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.25);
//    pad2->SetTopMargin(0);
//    pad2->SetBottomMargin(0.25);
//    pad2->Draw();
//    pad2->cd();
//    TLegend* errorLegend = new TLegend(0.6, 0.74, 0.8, 0.94);
//    TH1F* errorHist;
//    double error;
//
//
//    for (int i = 0; i < 4; i++) //systematics: Z-gamma
//    {
//        pad1->cd();
//        hist = new TH1D(Form("Z#gamma#gamma%d",i),Form("Z#gamma#gamma%d",i), 100u, 80, 200);
//        for (int j = 0, k = 6; k <= 8; j++, k++)
//        {
//            hist->Add(&resultmapHists[k][nominal_systematics[i]], SFs[j] / *Nodes[j].GetResultPtr<float>());
//        }
//        hist->SetLineColor(nominalColors[i]);
//        std::cout <<  hist->Integral() << '\n';
//
//        if (i == 0)
//        {
//            nominalhist = static_cast<TH1D*>(hist->Clone());
//            hist->SetTitle(";m_{ll#gamma#gamma} [GeV];Events");
//            hist->SetTitle("Z#gamma#gamma");
//            hist->GetXaxis()->SetTitle("m_{ll#gamma#gamma} [GeV]");
//            hist->GetXaxis()->SetTitleOffset(1.2);
//            hist->GetYaxis()->SetTitle("Events");
//            hist->GetYaxis()->CenterTitle(true);
//            hist->SetMaximum(4);
//            hist->Draw("HIST");
//        }
//        else
//        {
//            hist->Draw("HISTsame");
//
//            //Now do the error plot
//
//            c1->cd();
//            pad2->cd();
//
//            errorHist = new TH1F(Form("errorHist_%d", i), "", nominalhist->GetNbinsX(),
//                                 nominalhist->GetXaxis()->GetXmin(),
//                                 nominalhist->GetXaxis()->GetXmax());
//            errorHist->SetLineColor(colors[i-1]);
//            errorHist->SetLineWidth(1);
//            errorHist->SetMarkerStyle(20);
//            errorHist->SetMarkerSize(0.9);
//            errorHist->GetXaxis()->SetLabelSize(0.1);
//            errorHist->GetYaxis()->SetLabelSize(0.1);
//            errorHist->GetXaxis()->SetTitleSize(0.1);
//            errorHist->GetYaxis()->SetTitleSize(0.1);
//            errorHist->SetMarkerColor(nominalColors[i]);
//
//            for (int bin = 1; bin <= nominalhist->GetNbinsX(); bin++)
//            {
//                double nominal = nominalhist->GetBinContent(bin);
//                if (!nominal)
//                {
//                    error = 0;
//                }
//                else
//                {
//                    double systematic = hist->GetBinContent(bin);
//                    error = (systematic - nominal) / nominal;
//                }
//                errorHist->SetBinContent(bin, error);
//            }
//
//            if (i == 1)
//            {
//                errorHist->SetTitle(";m_{ll#gamma#gamma} [GeV];#frac{syst-nominal}{nominal}");
//                errorHist->GetXaxis()->SetTitle("m_{ll#gamma#gamma} [GeV]");
//                errorHist->GetXaxis()->SetTitleOffset(1.2);
//                errorHist->GetYaxis()->SetTitle("#frac{syst-nominal}{nominal}");
//                errorHist->GetYaxis()->SetTitleOffset(0.4);
//                errorHist->GetYaxis()->CenterTitle(true);
//                errorHist->SetMaximum(9);
//                errorHist->SetMinimum(-1.2);
//                errorHist->Draw("histsame");
//
//            }
//            else
//            {
//                c1->cd();
//                pad2->cd();
//                errorHist->Draw("histsame");
//            }
//
//            errorLegend->AddEntry(errorHist, nominalSystematics[i], "l");
//
//            c1->cd();
//            c1->Update();
//            pad2->cd();
//            errorLegend->SetBorderSize(0);
//            errorLegend->Draw("same");
//        }
//        legend->AddEntry(hist, nominalSystematics[i], "l");
//    }
//    pad1->cd();
//    Tl.SetTextSize(0.03);
//    Tl.DrawLatexNDC(0.6, 0.85, "#it{ATLAS} Internal");
//    Tl.DrawLatexNDC(0.6, 0.775,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
//    legend->SetBorderSize(0);
//    legend->Draw("same");
//    c1->SaveAs("Zgamma_resolved_higgs_syst_up.pdf");
//
//    //Z-jets
//    std::array<double,9> JetNumeratorSFs = {((139e15)*(1.9828e-9)*(0.821204)),((139e15)*(110.64e-12)*(0.69275)),((139e15)*(40.645e-12)*(0.615906)),((139e15)*(1.9817e-9)*(0.1136684)),((139e15)*(110.47e-12)*(0.1912956)),((139e15)*(40.674e-12)*(0.2326772)),((139e15)*(1.9819e-9)*(0.0656969)),((139e15)*(110.53e-12)*(0.1158741)),((139e15)*(40.68e-12)*(0.1535215))}; //numerators for jet bkg
//
//    c1 = new TCanvas();
//    legend = new TLegend(0.51, 0.3, 0.85, 0.7);
//    gStyle->SetOptStat(0);
//    c1->cd();
//    pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 0.95);
//    pad1->SetTopMargin(0.1);
//    pad1->SetBottomMargin(0);
//    pad1->Draw();
//    pad1->cd();
//    gStyle->SetOptStat(0);
//    c1->cd();
//    pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.25);
//    pad2->SetTopMargin(0);
//    pad2->SetBottomMargin(0.25);
//    pad2->Draw();
//    pad2->cd();
//    errorLegend = new TLegend(0.12, 0.7, 0.32, 0.9);
//
//    for (int i = 0; i < 4; i++) //systematics: Z-jets
//    {
//        pad1->cd();
//        hist = new TH1D(Form("Z+Jets%d",i),Form("Z+Jets%d",i), 100u, 80, 200);
//        for (int j = 0, k = 9; k <= 17; j++, k++)
//        {
//            hist->Add(&resultmapHists[k][nominal_systematics[i]], JetNumeratorSFs[j] / *Nodes[j+3].GetResultPtr<float>());
//        }
//        hist->SetLineColor(nominalColors[i]);
//        std::cout <<  hist->Integral() << '\n';
//
//        if (i == 0)
//        {
//            nominalhist = static_cast<TH1D*>(hist->Clone());
//            hist->SetTitle(";m_{ll#gamma#gamma} [GeV];Events");
//            hist->SetTitle("Z+Jets");
//            hist->GetXaxis()->SetTitle("m_{ll#gamma#gamma} [GeV]");
//            hist->GetXaxis()->SetTitleOffset(1.2);
//            hist->GetYaxis()->SetTitle("Events");
//            hist->GetYaxis()->CenterTitle(true);
//            hist->SetMaximum(2000);
//            hist->Draw("HIST");
//        }
//        else
//        {
//            hist->Draw("HISTsame");
//
//            //Now do the error plot
//
//            c1->cd();
//            pad2->cd();
//            errorHist = new TH1F(Form("errorHist_%d", i), "", nominalhist->GetNbinsX(),
//                                 nominalhist->GetXaxis()->GetXmin(),
//                                 nominalhist->GetXaxis()->GetXmax());
//            errorHist->SetLineColor(colors[i-1]);
//            errorHist->SetLineWidth(1);
//            errorHist->SetMarkerStyle(20);
//            errorHist->SetMarkerSize(0.9);
//            errorHist->GetXaxis()->SetLabelSize(0.1);
//            errorHist->GetYaxis()->SetLabelSize(0.1);
//            errorHist->GetXaxis()->SetTitleSize(0.1);
//            errorHist->GetYaxis()->SetTitleSize(0.1);
//            errorHist->SetMarkerColor(nominalColors[i]);
//
//            for (int bin = 1; bin <= nominalhist->GetNbinsX(); bin++)
//            {
//                double nominal = nominalhist->GetBinContent(bin);
//                if (!nominal)
//                {
//                    error = 0;
//                }
//                else
//                {
//                    double systematic = hist->GetBinContent(bin);
//                    error = (systematic - nominal) / nominal;
//                }
//                errorHist->SetBinContent(bin, error);
//            }
//
//            if (i == 1)
//            {
//                errorHist->SetTitle(";m_{ll#gamma#gamma} [GeV];#frac{syst-nominal}{nominal}");
//                errorHist->GetXaxis()->SetTitle("m_{ll#gamma#gamma} [GeV]");
//                errorHist->GetXaxis()->SetTitleOffset(1.2);
//                errorHist->GetYaxis()->SetTitle("#frac{syst-nominal}{nominal}");
//                errorHist->GetYaxis()->SetTitleOffset(0.4);
//                errorHist->GetYaxis()->CenterTitle(true);
//                errorHist->SetMaximum(17.3);
//                errorHist->SetMinimum(-1.2);
//                errorHist->Draw("histsame");
//
//            }
//            else
//            {
//                c1->cd();
//                pad2->cd();
//                errorHist->Draw("histsame");
//            }
//
//            errorLegend->AddEntry(errorHist, nominalSystematics[i], "l");
//
//            c1->cd();
//            c1->Update();
//            pad2->cd();
//            errorLegend->SetBorderSize(0);
//            errorLegend->Draw("same");
//        }
//        legend->AddEntry(hist, nominalSystematics[i], "l");
//    }
//    pad1->cd();
//    Tl.SetTextSize(0.03);
//    Tl.DrawLatexNDC(0.6, 0.85, "#it{ATLAS} Internal");
//    Tl.DrawLatexNDC(0.6, 0.775,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
//    legend->SetBorderSize(0);
//    legend->Draw("same");
//    c1->SaveAs("Zjets_resolved_higgs_syst_up.pdf");
//
//}
//
//void Resolved_Down_Higgs_Mass()
//{
//    Event::systematics =
//    {
////        "PH_EFF_ISO_Uncertainty",
////        "PH_EFF_ISO_Uncertainty",
////        "EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR",
////        "PRW_DATASF",
////        "MUON_EFF_RECO_SYS",
////        "MUON_EFF_ISO_SYS",
////        "MUON_EFF_TrigSystUncertainty",
////        "EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR",
////        "MUON_EFF_TrigStatUncertainty",
////        "MUON_EFF_RECO_STAT",
////        "MUON_EFF_TTVA_STAT",
////        "EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR",
////        "EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR",
////        "MUON_EFF_TTVA_SYS",
////        "MUON_EFF_ISO_STAT",
////        "MUON_SAGITTA_RHO",
////        "EG_RESOLUTION_ALL",
////        "EG_SCALE_ALL",
////        "MUON_MS",
////        "MUON_ID",
////        "EL_EFF_TriggerEff_TOTAL_1NPCOR_PLUS_UNCOR",
////        "MUON_SAGITTA_RESBIAS",
////        "MUON_SCALE",
////        "EG_RESOLUTION_ALL__1up",
////        "EG_SCALE_ALL__1up",
////        "PH_EFF_ISO_Uncertainty__1up",
////        "PH_EFF_ID_Uncertainty__1up",
////        "PH_EFF_TRIGGER_Uncertainty__1up",
//        "EG_RESOLUTION_ALL__1down",
//        "EG_SCALE_ALL__1down",
////        "PH_EFF_ISO_Uncertainty__1down",
//        "PH_EFF_ID_Uncertainty__1down",
////        "PH_EFF_TRIGGER_Uncertainty__1down",
//    };
//
//    std::vector<std::vector<std::string>> input_filenames = {
//        //Signal
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root"}, //1 GeV
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600751.PhPy8EG_AZNLO_ggH125_mA2p0_v1.root"}, // 2 GeV
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600752.PhPy8EG_AZNLO_ggH125_mA3p0_v1.root"}, // 3 GeV
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root"}, //5 GeV
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600756.PhPy8EG_AZNLO_ggH125_mA9p0_v1.root"}, // 9 GeV
//        //Data
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_data_test.root"},
//        //Z gamma background
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617070._000001.LGNTuple.root"},
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617064._000001.LGNTuple.root"},
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617074._000001.LGNTuple.root"},
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
//    std::vector<const char*> prefixes = {"Sig m_{A} = 1 GeV", "Sig m_{A} = 2 GeV", "Sig m_{A} = 3 GeV", "Sig m_{A} = 5 GeV", "Sig m_{A} = 9 GeV", "data", "pty2_9_17", "pty_17_myy_0_80", "pty_17_myy_80", "Zee_lightJet_0-70", "Zee_lightJet_70-140", "Zee_lightJet_140-280", "Zee_cJet_0-70", "Zee_cJet_70-140", "Zee_cJet_140-280", "Zee_bJet_0-70", "Zee_bJet_70-140", "Zee_bJet_140-280"};
//
//    std::vector<ROOT::RDF::RResultHandle> Nodes;
//
//    int count = 0;
//
//    std::vector<RResultMap<float>> resultmaps;
//    std::vector<RResultMap<TH1D>> resultmapHists;
//
//    auto verbosity = ROOT::Experimental::RLogScopedVerbosity(ROOT::Detail::RDF::RDFLogChannel(), ROOT::Experimental::ELogLevel::kInfo);
//
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
//        [](RVec<float> photon_id_eff, RVec<float> photon_iso_eff, RVec<float> photon_trg_eff/*, RVec<float> ei_event_weights_generator*/)
//        {
//            auto ResizeVal = std::max({photon_id_eff.size(), photon_iso_eff.size(), photon_trg_eff.size()});
//            photon_id_eff.resize(ResizeVal,1);
//            photon_iso_eff.resize(ResizeVal,1);
//            photon_trg_eff.resize(ResizeVal,1);
//
//            return photon_id_eff;
//            // *photon_iso_eff*photon_trg_eff; //Returns a vector of efficiencies
//
//        }, {"photon_id_eff", "photon_iso_eff", "photon_trg_eff",/* "ei_event_weights_generator"*/});
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
//        [](RVec<AbstractParticle> electrons)
//        {
//            //keep the electrons in each event that have pt > 20 GeV, |η| < 2.37,
//            //|η| not between 1.37 and 1.52, and that satisfy a medium id criteria `electron_id_medium`
//            electrons.erase(std::remove_if(electrons.begin(),electrons.end(),
//            [](AbstractParticle& ep)
//            {
//                 return (!((ep.electron_pt/1e3 > 20) && (std::abs(ep.electron_eta) < 2.37) &&
//                 (!((1.37 < std::abs(ep.electron_eta)) && (std::abs(ep.electron_eta) < 1.52)))
//                 && (ep.electron_id_medium == 1)));
//
//            }), electrons.end());
//
//            return electrons;
//
//        },{"abstract_electrons"}).Filter(
//        [](RVec<Muon>& muons, RVec<AbstractParticle>& electrons)
//        {
//            return (electrons.size()==2 && muons.empty()); //keep events which have exactly 2 electrons for di_electrons and no muons
//
//        }, {"muons", "di_electrons"});
//
//        //new dataframe node: contains only the events from `two_leptons` whose electrons in the `di_electrons` branch have opposite charge
//        auto opp_charge = two_leptons
//        .Filter([](RVec<AbstractParticle>& electrons)
//        {
//            return (electrons[0].electron_charge*electrons[1].electron_charge < 0);
//
//        }, {"di_electrons"});
//
//        //new dataframe node: contains only the events from `opp_charge` that have 1 electron with pt > 20 GeV and the other with pt > 27 GeV
//        auto leadingPt = opp_charge.Filter([](RVec<AbstractParticle>& electrons)
//        {
//            return ((electrons[0].electron_pt > 20e3 && electrons[1].electron_pt > 27e3) || (electrons[1].electron_pt > 20e3 && electrons[0].electron_pt > 27e3));
//        }, {"di_electrons"});
//
//        auto deltaR = leadingPt.Filter([] (RVec<AbstractParticle>& electrons)
//        {
//            return (DeltaR(electrons[0].ElectronVector(), electrons[1].ElectronVector()) > 0.01);
//        }, {"di_electrons"});
//
//        //new dataframe node: Contains a new column `dilep` in addition to the ones in `same_flavour` that stores the di-electron four-vector
//        auto dilep = deltaR.Define("dilep",[] (RVec<AbstractParticle>& electrons)
//        {
//            return (electrons[0].ElectronVector() + electrons[1].ElectronVector());
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
//        auto photon_passes_cuts = ptCut
//        .Define("photons_pass_cut_indices",
//        [&](const RVec<AbstractParticle>& photons)
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
//        }, {"abstract_photons"})
//        .Define("photons_pass_cuts",
//        [&](RVec<AbstractParticle>& photons, RVec<int>& photon_indices)
//        {
//            return Take(photons, photon_indices);
//        }, {"abstract_photons", "photons_pass_cut_indices"});
//
//        //New dataframe node: contains only events from `photon_passes_cuts` that pass the resolved category
//        auto resolved = photon_passes_cuts.Define("chosen_two_indices",
//        [](RVec<AbstractParticle>& photons_pass_cuts)
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
//                delta_r = DeltaR(photons_pass_cuts[combs[0][i]].PhotonVector(), photons_pass_cuts[combs[1][i]].PhotonVector());
//                m = (photons_pass_cuts[combs[0][i]].PhotonVector() + photons_pass_cuts[combs[1][i]].PhotonVector()).M();
//                pt = (photons_pass_cuts[combs[0][i]].PhotonVector() + photons_pass_cuts[combs[1][i]].PhotonVector()).Pt();
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
//        [&](RVec<AbstractParticle>& reco_photons_matched, RVec<unsigned long>& indices)
//        {
//            return Take(reco_photons_matched, indices);
//
//        }, {"photons_pass_cuts", "chosen_two_indices"})
//        .Define("diphoton",
//        [&](RVec<AbstractParticle>& reco_photons_matched)
//        {
//            return reco_photons_matched[0].PhotonVector() + reco_photons_matched[1].PhotonVector();
//        }, {"chosen_two"})
//        .Define("reco_higgs_mass",
//        [&](PtEtaPhiEVector& diphoton, PtEtaPhiEVector& dilep)
//        {
//            return (dilep + diphoton).M() / 1e3;
//        }, {"diphoton", "dilep"})
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
//        if (count >= 6)
//        {
//            Nodes.push_back(EventWeight.Sum<float>("EventWeight"));
//        }
//        resultmapHists.push_back(VariationsFor(resolved.Histo1D<double>({prefixes[count], prefixes[count++], 100u, 80, 200}, "reco_higgs_mass", "totEventWeight")));
//
//    }
//
////EG_RESOLUTION_ALL__1down:           42864.7 ID 1.00524  ISO 1.02283  TRIG 1.04157
////EG_RESOLUTION_ALL__1up:             42751.5 ID 1.00524  ISO 1.02283  TRIG 1.04157
////EG_SCALE_ALL__1down:                42521.1 ID 1.00524  ISO 1.02283  TRIG 1.04157
////EG_SCALE_ALL__1up:                  43059.8 ID 1.00524  ISO 1.02283  TRIG 1.04157
////PH_EFF_ID_Uncertainty__1down:       42790.4 ID 0.992328 ISO 1.02283  TRIG 1.04157
////PH_EFF_ID_Uncertainty__1up:         42790.4 ID 1.01816  ISO 1.02283  TRIG 1.04157
////PH_EFF_ISO_Uncertainty__1down:      42790.4 ID 1.00524  ISO 1.00971  TRIG 1.04157
////PH_EFF_ISO_Uncertainty__1up:        42790.4 ID 1.00524  ISO 1.03595  TRIG 1.04157
////PH_EFF_TRIGGER_Uncertainty__1down:  42790.4 ID 1.00524  ISO 1.02283  TRIG 1.01646
////PH_EFF_TRIGGER_Uncertainty__1up:    42790.4 ID 1.00524  ISO 1.02283  TRIG 1.06668
//
//    ROOT::RDF::RunGraphs(Nodes); // running all computation nodes concurrently
//
////    resultmapHists
////    ==============
////    0   //ma1
////    1   //ma2
////    2   //ma3
////    3   //ma5
////    4   //ma9
////    5   //data
////    6   //Z-gamma
////    7   //Z-gamma
////    8   //Z-gamma
////    9   //Z-jets
////    10  //Z-jets
////    11  //Z-jets
////    12  //Z-jets
////    13  //Z-jets
////    14  //Z-jets
////    15  //Z-jets
////    16  //Z-jets
////    17  //Z-jets
////
////    resultmaps
////    ==========
////
////    0       //Z-gamma
////    1       //Z-gamma
////    2       //Z-gamma
////    3       //Z-jets
////    4       //Z-jets
////    5       //Z-jets
////    6       //Z-jets
////    7       //Z-jets
////    8       //Z-jets
////    9       //Z-jets
////    10      //Z-jets
////    11      //Z-jets
//
//    constexpr std::array<const char*, 3> systematics = {"photons_and_electrons:EG_RESOLUTION_ALL__1down", "photons_and_electrons:EG_SCALE_ALL__1down",
//        "photon_id_eff:PH_EFF_ID_Uncertainty__1down"};
//    constexpr std::array<const char*, 4> nominal_systematics =
//    {
//        "nominal",
//        "photons_and_electrons:EG_RESOLUTION_ALL__1down",
//        "photons_and_electrons:EG_SCALE_ALL__1down",
//        "photon_id_eff:PH_EFF_ID_Uncertainty__1down",
//    };
//    constexpr std::array<const char*, 3> Systematics = {"EG_RESOLUTION_ALL__1down", "EG_SCALE_ALL__1down", "PH_EFF_ID_Uncertainty__1down"};
//    constexpr std::array<const char*, 4> nominalSystematics = {"nominal", "EG_RESOLUTION_ALL__1down", "EG_SCALE_ALL__1down", "PH_EFF_ID_Uncertainty__1down"};
//    constexpr std::array<const char*, 6> signalPlusDataSamples = {"Sig m_{A} = 1 GeV", "Sig m_{A} = 2 GeV", "Sig m_{A} = 3 GeV", "Sig m_{A} = 5 GeV", "Sig m_{A} = 9 GeV", "Data"};
//    constexpr std::array<const char*, 6> signalPlusDataFileNames = {"prompt_ma1_resolved_higgs_syst_down.pdf", "prompt_ma2_resolved_higgs_syst_down.pdf", "prompt_ma3_resolved_higgs_syst_down.pdf", "prompt_ma5_resolved_higgs_syst_down.pdf", "prompt_ma9_resolved_higgs_syst_down.pdf", "data_resolved_higgs_syst_down.pdf"};
//    constexpr std::array<EColor, 3> colors = {kBlue, kRed, static_cast<EColor>(kOrange+1)};
//    constexpr std::array<EColor, 4> nominalColors = {kBlack, kBlue, kRed, static_cast<EColor>(kOrange+1)};
//    constexpr std::array<int, 5> maxima = {62, 80, 270, 155, 155};
//    constexpr std::array<float, 6> signalPlusDataResidualMaxima = {2.2, 2.6, 1.1, 1.05, 1.2, 1.05};
//    constexpr std::array<float, 6> signalPlusDataResidualMinima = {-1.1, -1.2, -1.1, -1.15, -1.2, 0};
//    constexpr std::array<std::array<float, 4>, 6> signalPlusDataResidualLegend =
//    {{
//        {0.12, 0.74, 0.32, 0.94},
//        {0.12, 0.74, 0.32, 0.94},
//        {0.67, 0.74, 0.87, 0.94},
//        {0.6, 0.74, 0.8, 0.94},
//        {0.12, 0.74, 0.32, 0.94},
//        {0.15, 0.74, 0.35, 0.94},
//    }};
//
//
//    TCanvas* c1;
//    TLegend* legend;
//    TLatex Tl;
//
//    //signal+data
//    for (int i = 0; i < 6; i++) //signal samples plus data
//    {
//        c1 = new TCanvas();
//        if (i == 5)
//        {
//            legend = new TLegend(0.15, 0.35, 0.35, 0.65);
//        }
//        else
//        {
//            legend = new TLegend(0.5, 0.25, 0.85, 0.65);
//        }
//        c1->cd();
//        TPad* pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 0.95);
//        pad1->SetTopMargin(0.1);
//        pad1->SetBottomMargin(0);
//        pad1->Draw();
//        pad1->cd();
//
//        gStyle->SetOptStat(0);
//        resultmapHists[i]["nominal"].SetLineColor(kBlack);
//        resultmapHists[i]["nominal"].SetTitle(";m_{ll#gamma#gamma} [GeV];Events");
//        resultmapHists[i]["nominal"].SetTitle(signalPlusDataSamples[i]);
//        resultmapHists[i]["nominal"].GetXaxis()->SetTitle("m_{ll#gamma#gamma} [GeV]");
//        resultmapHists[i]["nominal"].GetXaxis()->SetTitleOffset(1.2);
//        resultmapHists[i]["nominal"].GetYaxis()->SetTitle("Events");
//        resultmapHists[i]["nominal"].GetYaxis()->CenterTitle(true);
//
//        legend->AddEntry(&resultmapHists[i]["nominal"], "nominal", "l");
//
//        resultmapHists[i]["nominal"].Draw("HIST");
//
//        for (int j = 0; j < 3; j++)
//        {
//            std::cout << "difference = "
//            << resultmapHists[i][systematics[j]].Integral() - resultmapHists[i]["nominal"].Integral() << '\n';
//            resultmapHists[i][systematics[j]].SetLineColor(colors[j]);
//            legend->AddEntry(&resultmapHists[i][systematics[j]], Systematics[j], "l");
//            resultmapHists[i][systematics[j]].Draw("HISTsame");
//        }
//
//        Tl.SetTextSize(0.03);
//        Tl.DrawLatexNDC(0.6, 0.8, "#it{ATLAS} Internal");
//        Tl.DrawLatexNDC(0.6, 0.7,"#sqrt{s} = 13 TeV");
//        legend->SetBorderSize(0);
//        legend->Draw("same");
//
//        // Create and draw the error plot
//        TH1F* errorHist;
//        double error;
//
//        c1->cd();
//        TPad* pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.25);
//        pad2->SetTopMargin(0);
//        pad2->SetBottomMargin(0.25);
//        pad2->Draw();
//        pad2->cd();
//        legend = new TLegend(signalPlusDataResidualLegend[i][0], signalPlusDataResidualLegend[i][1], signalPlusDataResidualLegend[i][2], signalPlusDataResidualLegend[i][3]);
//
//        for (int j = 0; j < 3; j++)
//        {
//            errorHist = new TH1F(Form("errorHist_%d%d", i, j), "", resultmapHists[i]["nominal"].GetNbinsX(),
//                                       resultmapHists[i]["nominal"].GetXaxis()->GetXmin(),
//                                       resultmapHists[i]["nominal"].GetXaxis()->GetXmax());
//            errorHist->SetLineColor(colors[j]);
//            errorHist->SetLineWidth(1);
//            errorHist->SetMarkerStyle(20);
//            errorHist->SetMarkerSize(0.9);
//            errorHist->GetXaxis()->SetLabelSize(0.1);
//            errorHist->GetYaxis()->SetLabelSize(0.1);
//            errorHist->GetXaxis()->SetTitleSize(0.1);
//            errorHist->GetYaxis()->SetTitleSize(0.1);
//            errorHist->SetMarkerColor(colors[j]);
//
//            for (int bin = 1; bin <= resultmapHists[i]["nominal"].GetNbinsX(); bin++)
//            {
//                double nominal = resultmapHists[i]["nominal"].GetBinContent(bin);
//                if (!nominal)
//                {
//                    error = 0;
//                }
//                else
//                {
//                    double systematic = resultmapHists[i][systematics[j]].GetBinContent(bin);
//                    error = (systematic - nominal) / nominal;
//                }
//                errorHist->SetBinContent(bin, error);
//            }
//
//            if (j == 0)
//            {
//                errorHist->SetTitle(";m_{ll#gamma#gamma} [GeV];#frac{syst-nominal}{nominal}");
//                errorHist->GetXaxis()->SetTitle("m_{ll#gamma#gamma} [GeV]");
//                errorHist->GetXaxis()->SetTitleOffset(1.2);
//                errorHist->GetYaxis()->SetTitle("#frac{syst-nominal}{nominal}");
//                errorHist->GetYaxis()->SetTitleOffset(0.4);
//                errorHist->GetYaxis()->CenterTitle(true);
//                errorHist->SetMaximum(signalPlusDataResidualMaxima[i]);
//                errorHist->SetMinimum(signalPlusDataResidualMinima[i]);
//                errorHist->Draw("histsame");
//            }
//            else
//            {
//                c1->cd();
//                pad2->cd();
//                errorHist->Draw("histsame");
//            }
//
//            legend->AddEntry(errorHist, Systematics[j], "l");
//
//            c1->cd();
//            c1->Update();
//
//        }
//        pad2->cd();
//        legend->SetBorderSize(0);
//        legend->Draw("same");
//        c1->SaveAs(signalPlusDataFileNames[i]);
//    }
//
//    constexpr std::array<double,3> SFs = {((139e15)*(.871e-12)),((139e15)*(.199e-12)), ((139e15)*(.0345e-15))}; //numerators for Z-gamma bkg
//
//
//    //Z-gamma
//    c1 = new TCanvas();
//    legend = new TLegend(0.15, 0.57, 0.42, 0.87);
//    TH1D *hist, *nominalhist;
//    gStyle->SetOptStat(0);
//    c1->cd();
//    TPad* pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 0.95);
//    pad1->SetTopMargin(0.1);
//    pad1->SetBottomMargin(0);
//    pad1->Draw();
//    pad1->cd();
//    gStyle->SetOptStat(0);
//    c1->cd();
//    TPad* pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.25);
//    pad2->SetTopMargin(0);
//    pad2->SetBottomMargin(0.25);
//    pad2->Draw();
//    pad2->cd();
//    TLegend* errorLegend = new TLegend(0.6, 0.74, 0.8, 0.94);
//    TH1F* errorHist;
//    double error;
//
//    for (int i = 0; i < 4; i++) //systematics: Z-gamma
//    {
//        pad1->cd();
//        hist = new TH1D(Form("Z#gamma#gamma%d",i),Form("Z#gamma#gamma%d",i), 100u, 80, 200);
//        for (int j = 0, k = 6; k <= 8; j++, k++)
//        {
//            hist->Add(&resultmapHists[k][nominal_systematics[i]], SFs[j] / *Nodes[j].GetResultPtr<float>());
//        }
//        hist->SetLineColor(nominalColors[i]);
//        std::cout <<  hist->Integral() << '\n';
//
//        if (i == 0)
//        {
//            nominalhist = static_cast<TH1D*>(hist->Clone());
//            hist->SetTitle(";m_{ll#gamma#gamma} [GeV];Events");
//            hist->SetTitle("Z#gamma#gamma");
//            hist->GetXaxis()->SetTitle("m_{ll#gamma#gamma} [GeV]");
//            hist->GetXaxis()->SetTitleOffset(1.2);
//            hist->GetYaxis()->SetTitle("Events");
//            hist->GetYaxis()->CenterTitle(true);
//            hist->SetMaximum(5);
//            hist->Draw("HIST");
//        }
//        else
//        {
//            hist->Draw("HISTsame");
//
//            //Now do the error plot
//
//            c1->cd();
//            pad2->cd();
//
//            errorHist = new TH1F(Form("errorHist_%d", i), "", nominalhist->GetNbinsX(),
//                                 nominalhist->GetXaxis()->GetXmin(),
//                                 nominalhist->GetXaxis()->GetXmax());
//            errorHist->SetLineColor(colors[i-1]);
//            errorHist->SetLineWidth(1);
//            errorHist->SetMarkerStyle(20);
//            errorHist->SetMarkerSize(0.9);
//            errorHist->GetXaxis()->SetLabelSize(0.1);
//            errorHist->GetYaxis()->SetLabelSize(0.1);
//            errorHist->GetXaxis()->SetTitleSize(0.1);
//            errorHist->GetYaxis()->SetTitleSize(0.1);
//            errorHist->SetMarkerColor(nominalColors[i]);
//
//            for (int bin = 1; bin <= nominalhist->GetNbinsX(); bin++)
//            {
//                double nominal = nominalhist->GetBinContent(bin);
//                if (!nominal)
//                {
//                    error = 0;
//                }
//                else
//                {
//                    double systematic = hist->GetBinContent(bin);
//                    error = (systematic - nominal) / nominal;
//                }
//                errorHist->SetBinContent(bin, error);
//            }
//
//            if (i == 1)
//            {
//                errorHist->SetTitle(";m_{ll#gamma#gamma} [GeV];#frac{syst-nominal}{nominal}");
//                errorHist->GetXaxis()->SetTitle("m_{ll#gamma#gamma} [GeV]");
//                errorHist->GetXaxis()->SetTitleOffset(1.2);
//                errorHist->GetYaxis()->SetTitle("#frac{syst-nominal}{nominal}");
//                errorHist->GetYaxis()->SetTitleOffset(0.4);
//                errorHist->GetYaxis()->CenterTitle(true);
//                errorHist->SetMaximum(6);
//                errorHist->SetMinimum(-1.2);
//                errorHist->Draw("histsame");
//
//            }
//            else
//            {
//                c1->cd();
//                pad2->cd();
//                errorHist->Draw("histsame");
//            }
//
//            errorLegend->AddEntry(errorHist, nominalSystematics[i], "l");
//
//            c1->cd();
//            c1->Update();
//            pad2->cd();
//            errorLegend->SetBorderSize(0);
//            errorLegend->Draw("same");
//        }
//        legend->AddEntry(hist, nominalSystematics[i], "l");
//    }
//    pad1->cd();
//    Tl.SetTextSize(0.03);
//    Tl.DrawLatexNDC(0.6, 0.85, "#it{ATLAS} Internal");
//    Tl.DrawLatexNDC(0.6, 0.775,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
//    legend->SetBorderSize(0);
//    legend->Draw("same");
//    c1->SaveAs("Zgamma_resolved_higgs_syst_down.pdf");
//
//    //Z-jets
//    std::array<double,9> JetNumeratorSFs = {((139e15)*(1.9828e-9)*(0.821204)),((139e15)*(110.64e-12)*(0.69275)),((139e15)*(40.645e-12)*(0.615906)),((139e15)*(1.9817e-9)*(0.1136684)),((139e15)*(110.47e-12)*(0.1912956)),((139e15)*(40.674e-12)*(0.2326772)),((139e15)*(1.9819e-9)*(0.0656969)),((139e15)*(110.53e-12)*(0.1158741)),((139e15)*(40.68e-12)*(0.1535215))}; //numerators for jet bkg
//
//    c1 = new TCanvas();
//    legend = new TLegend(0.51, 0.37, 0.85, 0.75);
//    gStyle->SetOptStat(0);
//    c1->cd();
//    pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 0.95);
//    pad1->SetTopMargin(0.1);
//    pad1->SetBottomMargin(0);
//    pad1->Draw();
//    pad1->cd();
//    gStyle->SetOptStat(0);
//    c1->cd();
//    pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.25);
//    pad2->SetTopMargin(0);
//    pad2->SetBottomMargin(0.25);
//    pad2->Draw();
//    pad2->cd();
//    errorLegend = new TLegend(0.12, 0.7, 0.32, 0.9);
//
//    for (int i = 0; i < 4; i++) //systematics: Z-jets
//    {
//        pad1->cd();
//        hist = new TH1D(Form("Z+Jets%d",i),Form("Z+Jets%d",i), 100u, 80, 200);
//        for (int j = 0, k = 9; k <= 17; j++, k++)
//        {
//            hist->Add(&resultmapHists[k][nominal_systematics[i]], JetNumeratorSFs[j] / *Nodes[j+3].GetResultPtr<float>());
//        }
//        hist->SetLineColor(nominalColors[i]);
//        std::cout <<  hist->Integral() << '\n';
//
//        if (i == 0)
//        {
//            nominalhist = static_cast<TH1D*>(hist->Clone());
//            hist->SetTitle(";m_{ll#gamma#gamma} [GeV];Events");
//            hist->SetTitle("Z+Jets");
//            hist->GetXaxis()->SetTitle("m_{ll#gamma#gamma} [GeV]");
//            hist->GetXaxis()->SetTitleOffset(1.2);
//            hist->GetYaxis()->SetTitle("Events");
//            hist->GetYaxis()->CenterTitle(true);
//            hist->Draw("HIST");
//        }
//        else
//        {
//            hist->Draw("HISTsame");
//
//            //Now do the error plot
//
//            c1->cd();
//            pad2->cd();
//            errorHist = new TH1F(Form("errorHist_%d", i), "", nominalhist->GetNbinsX(),
//                                 nominalhist->GetXaxis()->GetXmin(),
//                                 nominalhist->GetXaxis()->GetXmax());
//            errorHist->SetLineColor(colors[i-1]);
//            errorHist->SetLineWidth(1);
//            errorHist->SetMarkerStyle(20);
//            errorHist->SetMarkerSize(0.9);
//            errorHist->GetXaxis()->SetLabelSize(0.1);
//            errorHist->GetYaxis()->SetLabelSize(0.1);
//            errorHist->GetXaxis()->SetTitleSize(0.1);
//            errorHist->GetYaxis()->SetTitleSize(0.1);
//            errorHist->SetMarkerColor(nominalColors[i]);
//
//            for (int bin = 1; bin <= nominalhist->GetNbinsX(); bin++)
//            {
//                double nominal = nominalhist->GetBinContent(bin);
//                if (!nominal)
//                {
//                    error = 0;
//                }
//                else
//                {
//                    double systematic = hist->GetBinContent(bin);
//                    error = (systematic - nominal) / nominal;
//                }
//                errorHist->SetBinContent(bin, error);
//            }
//
//            if (i == 1)
//            {
//                errorHist->SetTitle(";m_{ll#gamma#gamma} [GeV];#frac{syst-nominal}{nominal}");
//                errorHist->GetXaxis()->SetTitle("m_{ll#gamma#gamma} [GeV]");
//                errorHist->GetXaxis()->SetTitleOffset(1.2);
//                errorHist->GetYaxis()->SetTitle("#frac{syst-nominal}{nominal}");
//                errorHist->GetYaxis()->SetTitleOffset(0.4);
//                errorHist->GetYaxis()->CenterTitle(true);
//                errorHist->SetMaximum(47.3);
//                errorHist->SetMinimum(-1.2);
//                errorHist->Draw("histsame");
//
//            }
//            else
//            {
//                c1->cd();
//                pad2->cd();
//                errorHist->Draw("histsame");
//            }
//
//            errorLegend->AddEntry(errorHist, nominalSystematics[i], "l");
//
//            c1->cd();
//            c1->Update();
//            pad2->cd();
//            errorLegend->SetBorderSize(0);
//            errorLegend->Draw("same");
//        }
//        legend->AddEntry(hist, nominalSystematics[i], "l");
//    }
//
//    pad1->cd();
//    Tl.SetTextSize(0.03);
//    Tl.DrawLatexNDC(0.6, 0.85, "#it{ATLAS} Internal");
//    Tl.DrawLatexNDC(0.6, 0.775,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
//    legend->SetBorderSize(0);
//    legend->Draw("same");
//    c1->SaveAs("Zjets_resolved_higgs_syst_down.pdf");
//}


//void Table21_Displaced_Axions()
//{
//    Event::systematics =
//    {
////        "PH_EFF_ISO_Uncertainty",
////        "PH_EFF_ISO_Uncertainty",
////        "EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR",
////        "PRW_DATASF",
////        "MUON_EFF_RECO_SYS",
////        "MUON_EFF_ISO_SYS",
////        "MUON_EFF_TrigSystUncertainty",
////        "EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR",
////        "MUON_EFF_TrigStatUncertainty",
////        "MUON_EFF_RECO_STAT",
////        "MUON_EFF_TTVA_STAT",
////        "EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR",
////        "EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR",
////        "MUON_EFF_TTVA_SYS",
////        "MUON_EFF_ISO_STAT",
////        "MUON_SAGITTA_RHO",
////        "EG_RESOLUTION_ALL",
////        "EG_SCALE_ALL",
////        "MUON_MS",
////        "MUON_ID",
////        "EL_EFF_TriggerEff_TOTAL_1NPCOR_PLUS_UNCOR",
////        "MUON_SAGITTA_RESBIAS",
////        "MUON_SCALE",
//
////        "PH_EFF_ISO_Uncertainty__1down",
////        "EG_SCALE_ALL__1down",
////        "PH_EFF_ID_Uncertainty__1down",
////        "PH_EFF_TRIGGER_Uncertainty__1down",
//        "EG_RESOLUTION_ALL__1up",
//        "EG_SCALE_ALL__1up",
//        "PH_EFF_ISO_Uncertainty__1up",
//        "PH_EFF_ID_Uncertainty__1up",
//        "PH_EFF_TRIGGER_Uncertainty__1up",
////        "EG_RESOLUTION_ALL__1down",
//    };
//
//    std::vector<std::vector<std::string>> input_filenames = {
//        //Signal
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root"},
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root"},
//        //Displaced Signal
////        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_FewMassPoints.root"},
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_FewMassPoints_PhotonSFs.root"},
//    };
//
//    constexpr std::array<const char*,15> prefixes = {R"--(Prompt Signal $m_{\text{A}}$ = 1 GeV)--", R"--(Prompt Signal $m_{\text{A}}$ = 5 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 0.2 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 0.5 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 1 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 3 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 5 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 10 GeV)--", /*R"--(Displaced Signal $m_{\text{A}}$ = 10.1 GeV)--",*/ R"--(Displaced Signal $m_{\text{A}}$ = 20 GeV)--",/* R"--(Displaced Signal $m_{\text{A}}$ = 20.1 GeV)--",*/ R"--(Displaced Signal $m_{\text{A}}$ = 29.5 GeV)--", };
//
//    std::vector<float> massPoints = {0.2, 0.5, 1, 3, 5, 10, 20, 29.5};//{0.2,0.5,1,3,5,10,10.1,20,20.1,29.5};
//    std::vector<ROOT::RDF::RResultHandle> Totals, Nodes, evtNodes;
//    std::vector<RResultMap<float>> resultmaps;
//    std::stringstream ss, cs;
//
//    int counter = 0;
//
//    for (auto& i: input_filenames)
//    {
//        SchottDataFrame df(MakeRDF(i,8));
//
//        auto EventWeight = df.Define("EventWeight",
//        [](const RVec<float>& ei_event_weights_generator)
//        {
//            return  ((ei_event_weights_generator[0]) ? 1 / ei_event_weights_generator[0] : 1);
//
//        }, {"ei_event_weights_generator"}).Define("truth_axions", [&](RVec<TruthParticle> truth_particles)
//        {
//            truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
//            [](TruthParticle& x)
//            {
//                return (std::abs(x.mc_pdg_id) != 36); //axions are pdg_id = 36 in displaced samples, but we don't use them in the prompt samples, so it doesn't matter that in the prompt samples, the pdg_id of the axions is 35
//
//            }), truth_particles.end());
//
//            return truth_particles;
//
//        }, {"truth_particles"}).Define("axion_masses", [&](RVec<TruthParticle>& truth_axions)
//        {
//            return truth_axions[0].mc_mass/1e3f;
//
//        }, {"truth_axions"});
//
//        auto trigger_selection = df.Define("truth_axions", [&](RVec<TruthParticle> truth_particles)
//        {
//            truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
//            [](TruthParticle& x)
//            {
//               return (std::abs(x.mc_pdg_id) != 36); //axions are pdg_id = 36 in displaced samples, but we don't use them in the prompt samples, so it doesn't matter that in the prompt samples, the pdg_id of the axions is 35
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
//                return false;
//            }
//            return true;
//
//        }, {"trigger_passed_triggers"});
//
//        auto two_leptons = trigger_selection.Define("di_electrons",
//        [](RVec<AbstractParticle> electrons)
//        {
//            electrons.erase(std::remove_if(electrons.begin(),electrons.end(),
//            [](AbstractParticle& ep)
//            {
//                //keep electron if it has pt > 20 GeV, |eta| < 2.37, not (1.37 < |eta| < 1.52), and id_medium = 1
//
//                return (!((ep.electron_pt/1e3 > 20) && (std::abs(ep.electron_eta) < 2.37) &&
//                (!((1.37 < std::abs(ep.electron_eta)) && (std::abs(ep.electron_eta) < 1.52)))
//                && (ep.electron_id_medium == 1)));
//
//            }), electrons.end());
//
//            return electrons;
//
//        },{"abstract_electrons"}).Filter([](const RVec<Muon>& muons, RVec<AbstractParticle> electrons)
//        {
//            return (electrons.size()==2 && muons.empty());
//
//        }, {"muons", "di_electrons"});
//
//        auto opp_charge = two_leptons.Filter([](const RVec<AbstractParticle>& electrons)
//        {
//            return (electrons[0].electron_charge*electrons[1].electron_charge < 0);
//
//        }, {"di_electrons"});
//
//        auto leadingPt = opp_charge.Filter([](const RVec<AbstractParticle>& electrons)
//        {
//            return ((electrons[0].electron_pt > 20e3 && electrons[1].electron_pt > 27e3) || (electrons[1].electron_pt > 20e3 && electrons[0].electron_pt > 27e3));
//        }, {"di_electrons"});
//
//        auto deltaR = leadingPt.Filter([] (const RVec<AbstractParticle>& electrons)
//        {
//            return (DeltaR(electrons[0].ElectronVector(), electrons[1].ElectronVector()) > 0.01);
//        }, {"di_electrons"});
//
//        auto mass = deltaR.Filter([] (const RVec<AbstractParticle>& electrons)
//        {
//            auto mass = (electrons[0].ElectronVector() + electrons[1].ElectronVector()).M()/1e3;
//            return ((mass >= 81) && (mass <= 101));
//        }, {"di_electrons"});
//
//        auto ptCut = mass.Filter([] (const RVec<AbstractParticle>& electrons)
//        {
//            auto pT = (electrons[0].ElectronVector() + electrons[1].ElectronVector()).Pt()/1e3;
//            return pT > 10;
//        }, {"di_electrons"});
//
//        auto photon_passes_cuts = ptCut.Define("abstract_photons_pass_cut_indices",[&](RVec<AbstractParticle>& p) //p = photon
//        {
//            RVec<int> x; //indices of photons that pass cuts
//
//            for (auto i = 0; i < p.size(); i++)
//            {
//                if (not ((std::abs(p[i].photon_eta) >= 2.37) or (std::abs(p[i].photon_eta) > 1.37 and std::abs(p[i].photon_eta) < 1.52) or (not p[i].photon_id_loose)))
//                {
//                    x.push_back(i);
//                }
//            }
//
//            return x;
//        }, {"abstract_photons"}).Define("photons_pass_cuts",
//        [&](RVec<AbstractParticle>& photons, RVec<int>& x)
//        {
//            return Take(photons, x); //Taking only the photons that passed the cuts in each event
//
//        }, {"abstract_photons", "abstract_photons_pass_cut_indices"}).Define("Event_number_and_photon",
//        [&](const RVec<AbstractParticle>& reco_photons_matched, UInt_t ei_event_number)
//        {
//            return std::make_pair(ei_event_number, reco_photons_matched);
//
//        }, {"photons_pass_cuts", "ei_event_number"});
//
//        auto merged_reco_photons_matched = photon_passes_cuts.Filter(
//        [&](const RVec<AbstractParticle>& reco_photons_matched)
//        {
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
//            for (size_t i=0; i<length; i++) //looping through all of the possible combinations of photons in each event
//            {
//                delta_r = DeltaR(reco_photons_matched[combs[0][i]].PhotonVector(), reco_photons_matched[combs[1][i]].PhotonVector());
//                m = (reco_photons_matched[combs[0][i]].PhotonVector() + reco_photons_matched[combs[1][i]].PhotonVector()).M();
//                pt = (reco_photons_matched[combs[0][i]].PhotonVector() + reco_photons_matched[combs[1][i]].PhotonVector()).Pt();
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
//        .Define("merged_photon_index",
//        [&](const RVec<AbstractParticle>& rpm) //rpm = reco photons matched
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
//        .Define("merged_photon",
//        [&](const RVec<AbstractParticle>& reco_photons_matched, int merged_photon_index)
//        {
//            return reco_photons_matched[merged_photon_index];
//
//        }, {"photons_pass_cuts", "merged_photon_index"});
//
//        auto dilepton_and_photon = merged_reco_photons_matched
//        .Define("reconstructed_mass",[&](const RVec<AbstractParticle>& di_electrons, const AbstractParticle& merged_photon)
//        {
//            auto four_momentum = di_electrons[0].ElectronVector() + di_electrons[1].ElectronVector();
//
//            return (four_momentum + merged_photon.PhotonVector()).M()/1e3;
//
//        }, {"di_electrons", "merged_photon"});
//
//        auto pSB = dilepton_and_photon.Filter(
//        [](const double reconstructed_mass)
//        {
//            return (reconstructed_mass < 110) || (reconstructed_mass > 130);
//        }, {"reconstructed_mass"});
//
//        auto pSR = dilepton_and_photon.Filter(
//        [](const double reconstructed_mass)
//        {
//            return (reconstructed_mass >= 110) && (reconstructed_mass <= 130);
//        }, {"reconstructed_mass"});
//
//        auto SB = pSB.Filter(
//        [](const RVec<float>& Eratio)
//        {
//            return (!Any(Eratio <= 0.8));
//        }, {"photon_shower_shape_e_ratio"});
//
//        auto SR = pSR.Filter(
//        [](const RVec<float>& Eratio)
//        {
//            return (!Any(Eratio <= 0.8));
//        }, {"photon_shower_shape_e_ratio"});
//
//        //mpi = merged_photon_index
//        auto totEventWeight = merged_reco_photons_matched.Define("totEventWeight", [](RVec<int>& x, int mpi, RVec<float> photon_id_eff, RVec<float> photon_iso_eff, RVec<float> photon_trg_eff/*, RVec<float> ei_event_weights_generator*/)
//        {
//            //   ||  jic they don't already have the same size  ||
//            //   \/                                             \/
//            auto ResizeVal = std::max({photon_id_eff.size(), photon_iso_eff.size(), photon_trg_eff.size()});
//            photon_id_eff.resize(ResizeVal,1);
//            photon_iso_eff.resize(ResizeVal,1);
//            photon_trg_eff.resize(ResizeVal,1);
//
//            //First, take the x indices from the respective photon_*_eff vectors,
//            //this corresponds to the photons from the set of all reco photons in
//            //the event that passed the "photon_passes_cuts" cuts. Then, from
//            //those photons, select the index mpi element that corresponds to
//            //the merged reco-photon
//
//            return Take(photon_id_eff, x)[mpi] * Take(photon_iso_eff, x)[mpi] * Take(photon_trg_eff, x)[mpi];
//
//        }, {"abstract_photons_pass_cut_indices", "merged_photon_index", "photon_id_eff", "photon_iso_eff", "photon_trg_eff",/* "ei_event_weights_generator"*/});
//
//        if (counter < 2)
//        {
//            Totals.push_back(EventWeight.Sum<float>("EventWeight"));
//            resultmaps.push_back(VariationsFor(totEventWeight.Sum<float>("totEventWeight")));
//        }
//        else
//        {
//            evtNodes.push_back(photon_passes_cuts.Take< std::pair<UInt_t, RVec<AbstractParticle>> , RVec<std::pair<UInt_t, RVec<AbstractParticle>>> >("Event_number_and_photon"));
//
//            for (auto& mass_point: massPoints)
//            {
//                auto mass_point_EventWeight = EventWeight.Filter([&]
//                (float axion_mass)
//                {
//                    return (roundToOneDecimalPlace(axion_mass) == mass_point);
////                    return ((axion_mass > mass_point - 0.1f) and (axion_mass < mass_point + 0.1f));
//                }, {"axion_masses"});
//
//                auto mass_point_totEventWeight = totEventWeight.Filter([&]
//                (float axion_mass)
//                {
//                    return (roundToOneDecimalPlace(axion_mass) == mass_point);
////                    return ((axion_mass > mass_point - 0.1f) and (axion_mass < mass_point + 0.1f));
//                }, {"axion_masses"});
//
//                Totals.push_back(mass_point_EventWeight.Sum<float>("EventWeight"));
//                resultmaps.push_back(VariationsFor(mass_point_totEventWeight.Sum<float>("totEventWeight")));
//                Nodes.push_back(mass_point_totEventWeight.Take<UInt_t, RVec<UInt_t>>("ei_event_number"));
//
//            }
//        }
//
//        counter++;
//    }
//
////EG_RESOLUTION_ALL__1down:           42864.7 ID 1.00524  ISO 1.02283  TRIG 1.04157
////EG_RESOLUTION_ALL__1up:             42751.5 ID 1.00524  ISO 1.02283  TRIG 1.04157
////EG_SCALE_ALL__1down:                42521.1 ID 1.00524  ISO 1.02283  TRIG 1.04157
////EG_SCALE_ALL__1up:                  43059.8 ID 1.00524  ISO 1.02283  TRIG 1.04157
////PH_EFF_ID_Uncertainty__1down:       42790.4 ID 0.992328 ISO 1.02283  TRIG 1.04157
////PH_EFF_ID_Uncertainty__1up:         42790.4 ID 1.01816  ISO 1.02283  TRIG 1.04157
////PH_EFF_ISO_Uncertainty__1down:      42790.4 ID 1.00524  ISO 1.00971  TRIG 1.04157
////PH_EFF_ISO_Uncertainty__1up:        42790.4 ID 1.00524  ISO 1.03595  TRIG 1.04157
////PH_EFF_TRIGGER_Uncertainty__1down:  42790.4 ID 1.00524  ISO 1.02283  TRIG 1.01646
////PH_EFF_TRIGGER_Uncertainty__1up:    42790.4 ID 1.00524  ISO 1.02283  TRIG 1.06668
//
//    ROOT::RDF::RunGraphs(Totals); // running all computation nodes concurrently
//
//    std::cout << evtNodes.size() << '\n';
//
//    for (int i = 0; i < evtNodes.size(); i++)
//    {
//        for (auto&j: *evtNodes[i].GetResultPtr<RVec<std::pair<UInt_t, RVec<AbstractParticle>>>>())
//        {
//            if (j.first == 370 and j.second.size() >= 2)
//            {
//                std::cout << j.first << " {";
////                for (auto& k: j.second)
////                {
////                    std::cout << "{photon pt = " << k.photon_pt <<
////                    ", photon eta = " << k.photon_eta <<
////                    ", photon phi = " << k.photon_phi <<
////                    "}, ";
////                }
////                std::cout << "} \n";
//
//                auto combs = Combinations(j.second, 2); //combinations of photons
//
//                size_t length = combs[0].size();
//
//                double delta_r, m, pt, X, best_X, pt1, pt2, chosen_delta_r;
//
//                for (size_t i=0; i<length; i++) //looping through all of the possible combinations of photons in each event
//                {
//                    std::cout << "{photon pt_1 = " << j.second[combs[0][i]].photon_pt
//                    << ", photon eta_1 = " << j.second[combs[0][i]].photon_eta
//                    << ", photon phi_1 = " << j.second[combs[0][i]].photon_phi
//                    << "; photon pt_2 = " << j.second[combs[1][i]].photon_pt
//                    << ", photon eta_2 = " << j.second[combs[1][i]].photon_eta
//                    << ", photon phi_2 = " << j.second[combs[1][i]].photon_phi;
//
//                    delta_r = DeltaR(j.second[combs[0][i]].PhotonVector(), j.second[combs[1][i]].PhotonVector());
//                    m = (j.second[combs[0][i]].PhotonVector() + j.second[combs[1][i]].PhotonVector()).M();
//                    pt = (j.second[combs[0][i]].PhotonVector() + j.second[combs[1][i]].PhotonVector()).Pt();
//                    X = delta_r*(pt/(2.0*m));
//
//                    std::cout << "; deltaR = " << delta_r <<
//                    "; X = " << X << "}, ";
//
//                    if (i==0 || ((std::abs(1-X) < std::abs(1-best_X)) and (delta_r < 1.5)))
//                    {
//                        best_X = X;
//                        pt1 = j.second[combs[0][i]].photon_pt;
//                        pt2 = j.second[combs[1][i]].photon_pt;
//                        chosen_delta_r = delta_r;
//                    }
//                }
//
//                std::cout << " Best deltaR = " << chosen_delta_r
//                << ", Best X = " << best_X << "}\n";
//
//            }
//        }
//    }
//
//    std::cout << "Merged Category Up\n==================\n";
//    for (int i = 0; i < Nodes.size(); i++)
//    {
//        std::cout << prefixes[i+2] << '\n';
//        std::cout << "==========\n";
//        for (auto& j: *Nodes[i].GetResultPtr<RVec<UInt_t>>())
//        {
//            std::cout << j << ", ";
//        }
//        std::cout << "\n\n";
//    }
//
////          resultmaps
////          ----------
////    0    1    ma1
////    2    3    ma5
////    4    5    displaced_axion_1
////    6    7    displaced_axion_2
////    8    9    displaced_axion_3
////    10   11   displaced_axion_4
////    12   13   displaced_axion_5
////    14   15   displaced_axion_6
////    16   17   displaced_axion_7
////    18   19   displaced_axion_8
////    20   21   displaced_axion_9
////    22   23   displaced_axion_10
//
//    ss << R"--(\section*{Table 21 Prompt and Displaced Signal Samples})--" << '\n';
//    ss << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
//    ss << R"--(\setlength\extrarowheight{2pt}\renewcommand{\arraystretch}{1.5})--" << '\n';
//    ss << R"--(\begin{tabular}{|c|c|c|c|c|c|})--" << '\n';
//    ss << R"--(\hline)--" << '\n';
//
//        ss << R"--(\multicolumn{6}{|c|}{\parbox{\linewidth}{\centering Merged Photon Category: Up Variations \\ (\% difference from nominal)}}\\[5 pt])--" << '\n';
//    ss << R"--(\hline)--" << '\n';
//
//    ss << R"--({Sample} & EG\_RESOLUTION\_ALL & EG\_SCALE\_ALL & PH\_EFF\_ISO\_Uncertainty & PH\_EFF\_ID\_Uncertainty & PH\_EFF\_TRIGGER\_Uncertainty \\ \hline)--" << '\n';
//    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    cs << R"--(\section*{Table 21 Prompt and Displaced Signal Samples})--" << '\n';
//    cs << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
//    cs << R"--(\setlength\extrarowheight{2pt}\renewcommand{\arraystretch}{1.5})--" << '\n';
//    cs << R"--(\begin{tabular}{|c|c|c|c|c|c|})--" << '\n';
//    cs << R"--(\hline)--" << '\n';
//
//    cs << R"--(\multicolumn{6}{|c|}{\parbox{\linewidth}{\centering Merged Photon Category: Up Variations \\ (raw number of events)}}\\[5 pt])--" << '\n';
//    cs << R"--(\hline)--" << '\n';
//
//    cs << R"--({Sample} & EG\_RESOLUTION\_ALL & EG\_SCALE\_ALL & PH\_EFF\_ISO\_Uncertainty & PH\_EFF\_ID\_Uncertainty & PH\_EFF\_TRIGGER\_Uncertainty \\ \hline)--" << '\n';
//
//    for (auto i = 0; (i < resultmaps.size()); i++)
//    {
//        float nominalVal = static_cast<float>(resultmaps[i]["nominal"]);
//
//        ss << prefixes[i] << " & ";
//        ss << std::setprecision(4) << std::fixed
//        << ( (nominalVal && resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1up"]) ? ((resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1up"]-nominalVal)/nominalVal)*100.0 : 0.0)
//        << " & " << std::setprecision(4) << std::fixed
//        << ( (nominalVal && resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1up"]) ? ((resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1up"]-nominalVal)/nominalVal)*100.0 : 0.0)
//        << " & " << std::setprecision(4) << std::fixed
//        << ( (nominalVal && resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1up"]) ? ((resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1up"]-nominalVal)/nominalVal)*100.0 : 0.0)
//        << " & " << std::setprecision(4) << std::fixed
//        << ( (nominalVal && resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1up"]) ? ((resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1up"]-nominalVal)/nominalVal)*100.0 : 0.0)
//        << " & " << std::setprecision(4) << std::fixed
//        << ( (nominalVal && resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up"]) ? ((resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up"]-nominalVal)/nominalVal)*100.0 : 0.0)
//        << R"--( \\ \hline)--" << '\n';
//        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//        cs << prefixes[i] << " & ";
//        cs << std::setprecision(4) << std::fixed
//        << ( (nominalVal && resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1up"]) ? resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1up"] : 0.0)
//        << " & " << std::setprecision(4) << std::fixed
//        << ( (nominalVal && resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1up"]) ? resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1up"] : 0.0)
//        << " & " << std::setprecision(4) << std::fixed
//        << ( (nominalVal && resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1up"]) ? resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1up"] : 0.0)
//        << " & " << std::setprecision(4) << std::fixed
//        << ( (nominalVal && resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1up"]) ? resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1up"] : 0.0)
//        << " & " << std::setprecision(4) << std::fixed
//        << ( (nominalVal && resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up"]) ? resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up"] : 0.0)
//        << R"--( \\ \hline)--" << '\n';
//
//        if (!resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1up"])
//        {
//            std::cout << i << ": photon_iso_eff:PH_EFF_ISO_Uncertainty__1up not found\n";
//        }
//        else
//        {
//            std::cout << i << ':' << ' ' << resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1up"] << ' ' << nominalVal << '\n';
//        }
//        if (!resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1up"])
//        {
//            std::cout << i << ": photon_id_eff:PH_EFF_ID_Uncertainty__1up not found\n";
//        }
//        else
//        {
//            std::cout << i << ':' << ' ' << resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1up"] << ' ' << nominalVal << '\n';
//        }
//        if (!resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up"])
//        {
//            std::cout << i << ": photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up not found\n";
//        }
//        else
//        {
//            std::cout << i << ':' << ' ' << resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up"] << ' ' << nominalVal << '\n';
//        }
//    }
//
//    ss << R"--(\end{tabular}})--" << '\n';
//    cs << R"--(\end{tabular}})--" << '\n';
//
//    ss << "\n\n\n";
//    cs << "\n\n\n";
//
//    std::cout << "\n\n\n" << ss.str();
//    std::cout << "\n\n\n" << cs.str();
//
//}
//
//void Table22_Displaced_Axions()
//{
//    Event::systematics =
//    {
//
//        "EG_RESOLUTION_ALL__1down",
//        "EG_SCALE_ALL__1down",
//        "PH_EFF_ISO_Uncertainty__1down",
//        "PH_EFF_ID_Uncertainty__1down",
//        "PH_EFF_TRIGGER_Uncertainty__1down",
////        "EG_RESOLUTION_ALL__1down",
//    };
//
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
//    std::vector<float> massPoints = {0.2,0.5,1,3,5,10,10.1,20,20.1,29.5};
//    std::vector<ROOT::RDF::RResultHandle> Totals, Nodes;
//    std::vector<RResultMap<float>> resultmaps;
//    std::stringstream ss, cs;
//
//    int counter = 0;
//
//    for (auto& i: input_filenames)
//    {
//        SchottDataFrame df(MakeRDF(i,8));
//
//        auto EventWeight = df.Define("EventWeight",
//        [](const RVec<float>& ei_event_weights_generator)
//        {
//            return  ((ei_event_weights_generator[0]) ? 1 / ei_event_weights_generator[0] : 1);
//
//        }, {"ei_event_weights_generator"}).Define("truth_axions", [&](RVec<TruthParticle> truth_particles)
//        {
//            truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
//            [](TruthParticle& x)
//            {
//             return (std::abs(x.mc_pdg_id) != 36);
//
//            }), truth_particles.end());
//
//            return truth_particles;
//
//        }, {"truth_particles"}).Define("axion_masses", [&](RVec<TruthParticle>& truth_axions)
//        {
//            return truth_axions[0].mc_mass/1e3f;
//
//        }, {"truth_axions"});
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
//                return false;
//            }
//            return true;
//
//        }, {"trigger_passed_triggers"});
//
//        auto two_leptons = trigger_selection.Define("di_electrons",
//        [](RVec<AbstractParticle> electrons)
//        {
//            electrons.erase(std::remove_if(electrons.begin(),electrons.end(),
//            [](AbstractParticle& ep)
//            {
//                //keep electron if it has pt > 20 GeV, |eta| < 2.37, not (1.37 < |eta| < 1.52), and id_medium = 1
//
//                return (!((ep.electron_pt/1e3 > 20) && (std::abs(ep.electron_eta) < 2.37) &&
//                (!((1.37 < std::abs(ep.electron_eta)) && (std::abs(ep.electron_eta) < 1.52)))
//                && (ep.electron_id_medium == 1)));
//
//            }), electrons.end());
//
//            return electrons;
//
//        },{"abstract_electrons"}).Filter([](const RVec<Muon>& muons, RVec<AbstractParticle> electrons)
//        {
//            return (electrons.size()==2 && muons.empty());
//
//        }, {"muons", "di_electrons"});
//
//        auto opp_charge = two_leptons.Filter([](const RVec<AbstractParticle>& electrons)
//        {
//            return (electrons[0].electron_charge*electrons[1].electron_charge < 0);
//
//        }, {"di_electrons"});
//
//        auto leadingPt = opp_charge.Filter([](const RVec<AbstractParticle>& electrons)
//        {
//            return ((electrons[0].electron_pt > 20e3 && electrons[1].electron_pt > 27e3) || (electrons[1].electron_pt > 20e3 && electrons[0].electron_pt > 27e3));
//        }, {"di_electrons"});
//
//        auto deltaR = leadingPt.Filter([] (const RVec<AbstractParticle>& electrons)
//        {
//            return (DeltaR(electrons[0].ElectronVector(), electrons[1].ElectronVector()) > 0.01);
//        }, {"di_electrons"});
//
//        auto mass = deltaR.Filter([] (const RVec<AbstractParticle>& electrons)
//        {
//            auto mass = (electrons[0].ElectronVector() + electrons[1].ElectronVector()).M()/1e3;
//            return ((mass >= 81) && (mass <= 101));
//        }, {"di_electrons"});
//
//        auto ptCut = mass.Filter([] (const RVec<AbstractParticle>& electrons)
//        {
//            auto pT = (electrons[0].ElectronVector() + electrons[1].ElectronVector()).Pt()/1e3;
//            return pT > 10;
//        }, {"di_electrons"});
//
//        auto photon_passes_cuts = ptCut.Define("abstract_photons_pass_cut_indices",[&](RVec<AbstractParticle>& p) //p = photon
//        {
//            RVec<int> x; //indices of photons that pass cuts
//
//            for (auto i = 0; i < p.size(); i++)
//            {
//                if (not ((std::abs(p[i].photon_eta) >= 2.37) or (std::abs(p[i].photon_eta) > 1.37 and std::abs(p[i].photon_eta) < 1.52) or (not p[i].photon_id_loose)))
//                {
//                    x.push_back(i);
//                }
//            }
//
//            return x;
//        }, {"abstract_photons"}).Define("photons_pass_cuts",
//        [&](RVec<AbstractParticle>& photons, RVec<int>& x)
//        {
//            return Take(photons, x); //Taking only the photons that passed the cuts in each event
//
//        }, {"abstract_photons", "abstract_photons_pass_cut_indices"});
//
//        auto merged_reco_photons_matched = photon_passes_cuts.Filter(
//        [&](const RVec<AbstractParticle>& reco_photons_matched)
//        {
////            RVec<AbstractParticle> reco_photons_matched = reco_photons_test;
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
//            for (size_t i=0; i<length; i++) //looping through all of the possible combinations of photons in each event
//            {
//                delta_r = DeltaR(reco_photons_matched[combs[0][i]].PhotonVector(), reco_photons_matched[combs[1][i]].PhotonVector());
//                m = (reco_photons_matched[combs[0][i]].PhotonVector() + reco_photons_matched[combs[1][i]].PhotonVector()).M();
//                pt = (reco_photons_matched[combs[0][i]].PhotonVector() + reco_photons_matched[combs[1][i]].PhotonVector()).Pt();
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
//        .Define("merged_photon_index",
//        [&](const RVec<AbstractParticle>& rpm) //rpm = reco photons matched
//        {
//            for (auto i = 0; i < rpm.size(); i++)
//            {
//                if (rpm[i].photon_pt > 20e3)
//                {
//                    return i; //returning the index of the first photon that has photon_pt > 20 GeV
//                }
//            }
//            return 0; //jic the compiler complains
//
//        }, {"photons_pass_cuts"})
//        .Define("merged_photon",
//        [&](const RVec<AbstractParticle>& reco_photons_matched, int merged_photon_index)
//        {
//            return reco_photons_matched[merged_photon_index];
//
//        }, {"photons_pass_cuts", "merged_photon_index"});
//
//        auto dilepton_and_photon = merged_reco_photons_matched
//        .Define("reconstructed_mass",[&](const RVec<AbstractParticle>& di_electrons, const AbstractParticle& merged_photon)
//        {
//            auto four_momentum = di_electrons[0].ElectronVector() + di_electrons[1].ElectronVector();
//
//            return (four_momentum + merged_photon.PhotonVector()).M()/1e3;
//
//        }, {"di_electrons", "merged_photon"});
//
//        auto pSB = dilepton_and_photon.Filter(
//        [](const double reconstructed_mass)
//        {
//            return (reconstructed_mass < 110) || (reconstructed_mass > 130);
//        }, {"reconstructed_mass"});
//
//        auto pSR = dilepton_and_photon.Filter(
//        [](const double reconstructed_mass)
//        {
//            return (reconstructed_mass >= 110) && (reconstructed_mass <= 130);
//        }, {"reconstructed_mass"});
//
//        auto SB = pSB.Filter(
//        [](const RVec<float>& Eratio)
//        {
//            return (!Any(Eratio <= 0.8));
//        }, {"photon_shower_shape_e_ratio"});
//
//        auto SR = pSR.Filter(
//        [](const RVec<float>& Eratio)
//        {
//            return (!Any(Eratio <= 0.8));
//        }, {"photon_shower_shape_e_ratio"});
//
//        //mpi = merged_photon_index
//        auto totEventWeight = merged_reco_photons_matched.Define("totEventWeight", [](RVec<int>& x, int mpi, RVec<float> photon_id_eff, RVec<float> photon_iso_eff, RVec<float> photon_trg_eff/*, RVec<float> ei_event_weights_generator*/)
//        {
//            //   ||  jic they don't already have the same size  ||
//            //   \/                                             \/
//            auto ResizeVal = std::max({photon_id_eff.size(), photon_iso_eff.size(), photon_trg_eff.size()});
//            photon_id_eff.resize(ResizeVal,1);
//            photon_iso_eff.resize(ResizeVal,1);
//            photon_trg_eff.resize(ResizeVal,1);
//
//            //First, take the x indices from the respective photon_*_eff vectors,
//            //this corresponds to the photons from the set of all reco photons in
//            //the event that passed the "photon_passes_cuts" cuts. Then, from
//            //those photons, select the index mpi element that corresponds to
//            //the merged reco-photon
//
//            return Take(photon_id_eff, x)[mpi] * Take(photon_iso_eff, x)[mpi] * Take(photon_trg_eff, x)[mpi];
//
//        }, {"abstract_photons_pass_cut_indices", "merged_photon_index", "photon_id_eff", "photon_iso_eff", "photon_trg_eff",/* "ei_event_weights_generator"*/});
//
//        if (counter < 2)
//        {
//            Totals.push_back(EventWeight.Sum<float>("EventWeight"));
//            resultmaps.push_back(VariationsFor(totEventWeight.Sum<float>("totEventWeight")));
//        }
//        else
//        {
//            for (auto& mass_point: massPoints)
//            {
//                auto mass_point_EventWeight = EventWeight.Filter([&]
//                (float axion_mass)
//                {
//                    return (roundToOneDecimalPlace(axion_mass) == mass_point);
////                    return ((axion_mass > mass_point - 0.1f) and (axion_mass < mass_point + 0.1f));
//                }, {"axion_masses"});
//
//                auto mass_point_totEventWeight = totEventWeight.Filter([&]
//                (float axion_mass)
//                {
//                    return (roundToOneDecimalPlace(axion_mass) == mass_point);
////                    return ((axion_mass > mass_point - 0.1f) and (axion_mass < mass_point + 0.1f));
//                }, {"axion_masses"});
//
//                Totals.push_back(mass_point_EventWeight.Sum<float>("EventWeight"));
//                resultmaps.push_back(VariationsFor(mass_point_totEventWeight.Sum<float>("totEventWeight")));
//                Nodes.push_back(mass_point_totEventWeight.Take<UInt_t, RVec<UInt_t>>("ei_event_number"));
//            }
//        }
//
//        counter++;
//    }
//
//    ROOT::RDF::RunGraphs(Totals); // running all computation nodes concurrently
//
//    std::cout << "Merged Category Down\n====================\n";
//    for (int i = 0; i < Nodes.size(); i++)
//    {
//        std::cout << prefixes[i+2] << '\n';
//        std::cout << "==========\n";
//        for (auto& j: *Nodes[i].GetResultPtr<RVec<UInt_t>>())
//        {
//            std::cout << j << ", ";
//        }
//        std::cout << "\n\n";
//    }
//
////          resultmaps
////          ----------
////    0    1    ma1
////    2    3    ma5
////    4    5    displaced_axion_1
////    6    7    displaced_axion_2
////    8    9    displaced_axion_3
////    10   11   displaced_axion_4
////    12   13   displaced_axion_5
////    14   15   displaced_axion_6
////    16   17   displaced_axion_7
////    18   19   displaced_axion_8
////    20   21   displaced_axion_9
////    22   23   displaced_axion_10
//
//    ss << R"--(\section*{Table 22 Prompt and Displaced Signal Samples})--" << '\n';
//    ss << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
//    ss << R"--(\setlength\extrarowheight{2pt}\renewcommand{\arraystretch}{1.5})--" << '\n';
//    ss << R"--(\begin{tabular}{|c|c|c|c|c|c|})--" << '\n';
//    ss << R"--(\hline)--" << '\n';
//
//        ss << R"--(\multicolumn{6}{|c|}{\parbox{\linewidth}{\centering Merged Photon Category: Down Variations \\ (\% difference from nominal)}}\\[5 pt])--" << '\n';
//    ss << R"--(\hline)--" << '\n';
//
//    ss << R"--({Sample} & EG\_RESOLUTION\_ALL & EG\_SCALE\_ALL & PH\_EFF\_ISO\_Uncertainty & PH\_EFF\_ID\_Uncertainty & PH\_EFF\_TRIGGER\_Uncertainty \\ \hline)--" << '\n';
//    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    cs << R"--(\section*{Table 22 Prompt and Displaced Signal Samples})--" << '\n';
//    cs << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
//    cs << R"--(\setlength\extrarowheight{2pt}\renewcommand{\arraystretch}{1.5})--" << '\n';
//    cs << R"--(\begin{tabular}{|c|c|c|c|c|c|})--" << '\n';
//    cs << R"--(\hline)--" << '\n';
//
//    cs << R"--(\multicolumn{6}{|c|}{\parbox{\linewidth}{\centering Merged Photon Category: Down Variations \\ (raw number of events)}}\\[5 pt])--" << '\n';
//    cs << R"--(\hline)--" << '\n';
//
//    cs << R"--({Sample} & EG\_RESOLUTION\_ALL & EG\_SCALE\_ALL & PH\_EFF\_ISO\_Uncertainty & PH\_EFF\_ID\_Uncertainty & PH\_EFF\_TRIGGER\_Uncertainty \\ \hline)--" << '\n';
//
//    for (auto i = 0; (i < resultmaps.size()); i++)
//    {
//        float nominalVal = static_cast<float>(resultmaps[i]["nominal"]);
//
//        ss << prefixes[i] << " & ";
//        ss << std::setprecision(4) << std::fixed
//        << ( (nominalVal && resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1down"]) ? ((resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1down"]-nominalVal)/nominalVal)*100.0 : 0.0)
//        << " & " << std::setprecision(4) << std::fixed
//        << ( (nominalVal && resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1down"]) ? ((resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1down"]-nominalVal)/nominalVal)*100.0 : 0.0)
//        << " & " << std::setprecision(4) << std::fixed
//        << ( (nominalVal && resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1down"]) ? ((resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1down"]-nominalVal)/nominalVal)*100.0 : 0.0)
//        << " & " << std::setprecision(4) << std::fixed
//        << ( (nominalVal && resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1down"]) ? ((resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1down"]-nominalVal)/nominalVal)*100.0 : 0.0)
//        << " & " << std::setprecision(4) << std::fixed
//        << ( (nominalVal && resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down"]) ? ((resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down"]-nominalVal)/nominalVal)*100.0 : 0.0)
//        << R"--( \\ \hline)--" << '\n';
//        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//        cs << prefixes[i] << " & ";
//        cs << std::setprecision(4) << std::fixed
//        << ( (nominalVal && resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1down"]) ? resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1down"] : 0.0)
//        << " & " << std::setprecision(4) << std::fixed
//        << ( (nominalVal && resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1down"]) ? resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1down"] : 0.0)
//        << " & " << std::setprecision(4) << std::fixed
//        << ( (nominalVal && resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1down"]) ? resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1down"] : 0.0)
//        << " & " << std::setprecision(4) << std::fixed
//        << ( (nominalVal && resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1down"]) ? resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1down"] : 0.0)
//        << " & " << std::setprecision(4) << std::fixed
//        << ( (nominalVal && resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down"]) ? resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down"] : 0.0)
//        << R"--( \\ \hline)--" << '\n';
//
//        if (!resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1down"])
//        {
//            std::cout << i << ": photon_iso_eff:PH_EFF_ISO_Uncertainty__1down not found\n";
//        }
//        else
//        {
//            std::cout << i << ':' << ' ' << resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1down"] << ' ' << nominalVal << '\n';
//        }
//        if (!resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1down"])
//        {
//            std::cout << i << ": photon_id_eff:PH_EFF_ID_Uncertainty__1down not found\n";
//        }
//        else
//        {
//            std::cout << i << ':' << ' ' << resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1down"] << ' ' << nominalVal << '\n';
//        }
//        if (!resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down"])
//        {
//            std::cout << i << ": photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down not found\n";
//        }
//        else
//        {
//            std::cout << i << ':' << ' ' << resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down"] << ' ' << nominalVal << '\n';
//        }
//    }
//
//    ss << R"--(\end{tabular}})--" << '\n';
//    cs << R"--(\end{tabular}})--" << '\n';
//
//    ss << "\n\n\n";
//    cs << "\n\n\n";
//
//    std::cout << "\n\n\n" << ss.str();
//    std::cout << "\n\n\n" << cs.str();
//
//}
//
//void Table23_Displaced_Axions()
//{
//    Event::systematics =
//    {
////        "PH_EFF_ISO_Uncertainty",
////        "PH_EFF_ISO_Uncertainty",
////        "EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR",
////        "PRW_DATASF",
////        "MUON_EFF_RECO_SYS",
////        "MUON_EFF_ISO_SYS",
////        "MUON_EFF_TrigSystUncertainty",
////        "EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR",
////        "MUON_EFF_TrigStatUncertainty",
////        "MUON_EFF_RECO_STAT",
////        "MUON_EFF_TTVA_STAT",
////        "EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR",
////        "EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR",
////        "MUON_EFF_TTVA_SYS",
////        "MUON_EFF_ISO_STAT",
////        "MUON_SAGITTA_RHO",
////        "EG_RESOLUTION_ALL",
////        "EG_SCALE_ALL",
////        "MUON_MS",
////        "MUON_ID",
////        "EL_EFF_TriggerEff_TOTAL_1NPCOR_PLUS_UNCOR",
////        "MUON_SAGITTA_RESBIAS",
////        "MUON_SCALE",
//
////        "PH_EFF_ISO_Uncertainty__1down",
////        "EG_SCALE_ALL__1down",
////        "PH_EFF_ID_Uncertainty__1down",
////        "PH_EFF_TRIGGER_Uncertainty__1down",
//        "EG_RESOLUTION_ALL__1up",
//        "EG_SCALE_ALL__1up",
//        "PH_EFF_ISO_Uncertainty__1up",
//        "PH_EFF_ID_Uncertainty__1up",
//        "PH_EFF_TRIGGER_Uncertainty__1up",
////        "EG_RESOLUTION_ALL__1down",
//    };
//
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
//    std::vector<float> massPoints = {0.2,0.5,1,3,5,10,10.1,20,20.1,29.5};
//    std::vector<ROOT::RDF::RResultHandle> Totals, Nodes;
//    std::vector<RResultMap<float>> resultmaps;
//    std::stringstream ss, cs;
//
//    int counter = 0;
//
//    for (auto& i: input_filenames)
//    {
//        SchottDataFrame df(MakeRDF(i,8));
//
//        auto EventWeight = df.Define("EventWeight",
//        [](const RVec<float>& ei_event_weights_generator)
//        {
//            return  ((ei_event_weights_generator[0]) ? 1 / ei_event_weights_generator[0] : 1);
//
//        }, {"ei_event_weights_generator"}).Define("truth_axions", [&](RVec<TruthParticle> truth_particles)
//        {
//            truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
//            [](TruthParticle& x)
//            {
//             return (std::abs(x.mc_pdg_id) != 36);
//
//            }), truth_particles.end());
//
//            return truth_particles;
//
//        }, {"truth_particles"}).Define("axion_masses", [&](RVec<TruthParticle>& truth_axions)
//        {
//            return truth_axions[0].mc_mass/1e3f;
//
//        }, {"truth_axions"});
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
//                return false;
//            }
//            return true;
//
//        }, {"trigger_passed_triggers"});
//
//        auto two_leptons = trigger_selection.Define("di_electrons",
//        [](RVec<AbstractParticle> electrons)
//        {
//            electrons.erase(std::remove_if(electrons.begin(),electrons.end(),
//            [](AbstractParticle& ep)
//            {
//                //keep electron if it has pt > 20 GeV, |eta| < 2.37, not (1.37 < |eta| < 1.52), and id_medium = 1
//
//                return (!((ep.electron_pt/1e3 > 20) && (std::abs(ep.electron_eta) < 2.37) &&
//                (!((1.37 < std::abs(ep.electron_eta)) && (std::abs(ep.electron_eta) < 1.52)))
//                && (ep.electron_id_medium == 1)));
//
//            }), electrons.end());
//
//            return electrons;
//
//        },{"abstract_electrons"}).Filter([](const RVec<Muon>& muons, RVec<AbstractParticle> electrons)
//        {
//            return (electrons.size()==2 && muons.empty());
//
//        }, {"muons", "di_electrons"});
//
//        auto opp_charge = two_leptons.Filter([](const RVec<AbstractParticle>& electrons)
//        {
//            return (electrons[0].electron_charge*electrons[1].electron_charge < 0);
//
//        }, {"di_electrons"});
//
//        auto leadingPt = opp_charge.Filter([](const RVec<AbstractParticle>& electrons)
//        {
//            return ((electrons[0].electron_pt > 20e3 && electrons[1].electron_pt > 27e3) || (electrons[1].electron_pt > 20e3 && electrons[0].electron_pt > 27e3));
//        }, {"di_electrons"});
//
//        auto deltaR = leadingPt.Filter([] (const RVec<AbstractParticle>& electrons)
//        {
//            return (DeltaR(electrons[0].ElectronVector(), electrons[1].ElectronVector()) > 0.01);
//        }, {"di_electrons"});
//
//        auto mass = deltaR.Filter([] (const RVec<AbstractParticle>& electrons)
//        {
//            auto mass = (electrons[0].ElectronVector() + electrons[1].ElectronVector()).M()/1e3;
//            return ((mass >= 81) && (mass <= 101));
//        }, {"di_electrons"});
//
//        auto ptCut = mass.Filter([] (const RVec<AbstractParticle>& electrons)
//        {
//            auto pT = (electrons[0].ElectronVector() + electrons[1].ElectronVector()).Pt()/1e3;
//            return pT > 10;
//        }, {"di_electrons"});
//
//        auto photon_passes_cuts = ptCut.Define("abstract_photons_pass_cut_indices",[&](RVec<AbstractParticle>& p) //p = photon
//        {
//            RVec<int> x; //indices of photons that pass cuts
//
//            for (auto i = 0; i < p.size(); i++)
//            {
//                if (not ((std::abs(p[i].photon_eta) >= 2.37) or (std::abs(p[i].photon_eta) > 1.37 and std::abs(p[i].photon_eta) < 1.52) or (not p[i].photon_id_loose)))
//                {
//                    x.push_back(i);
//                }
//            }
//
//            return x;
//        }, {"abstract_photons"}).Define("photons_pass_cuts",
//        [&](RVec<AbstractParticle>& photons, RVec<int>& x)
//        {
//            return Take(photons, x); //Taking only the photons that passed the cuts in each event
//
//        }, {"abstract_photons", "abstract_photons_pass_cut_indices"});
//
//        auto resolved_reco_photons_matched = photon_passes_cuts.Define("reco_photons_matched_indices",
//        [](const RVec<AbstractParticle>& reco_photons_matched)
//        {
//            RVec<unsigned long> x; //vector of indices
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
//                delta_r = DeltaR(reco_photons_matched[combs[0][i]].PhotonVector(), reco_photons_matched[combs[1][i]].PhotonVector());
//                m = (reco_photons_matched[combs[0][i]].PhotonVector() + reco_photons_matched[combs[1][i]].PhotonVector()).M();
//                pt = (reco_photons_matched[combs[0][i]].PhotonVector() + reco_photons_matched[combs[1][i]].PhotonVector()).Pt();
//                X = delta_r*(pt/(2.0*m));
//                if (i==0 || ((std::abs(1-X) < std::abs(1-best_X)) and (delta_r < 1.5)))
//                {
//                    best_X = X;
//                    pt1 = reco_photons_matched[combs[0][i]].photon_pt;
//                    pt2 = reco_photons_matched[combs[1][i]].photon_pt;
//                    chosen_delta_r = delta_r;
//                    x = {combs[0][i], combs[1][i]};
//                }
//            }
//            if (pt1 > 10e3 && pt2 > 10e3 && best_X > 0.96 && best_X < 1.2 && chosen_delta_r < 1.5)
//            {
//                return x;
//            }
//
//            x.clear();
//            return x;
//
//        }, {"photons_pass_cuts"}).Filter(
//        [&](RVec<unsigned long>& reco_photons_matched_indices)
//        {
//            return (reco_photons_matched_indices.size()==2);
//
//        }, {"reco_photons_matched_indices"}).Define("chosen_two",
//        [](const RVec<AbstractParticle>& reco_photons_matched, RVec<unsigned long>& reco_photons_matched_indices)
//        {
//            return Take(reco_photons_matched, reco_photons_matched_indices);
//
//        }, {"photons_pass_cuts", "reco_photons_matched_indices"}).Define("reconstructed_mass",
//        [&](RVec<AbstractParticle>& diph)
//        {
//            return (diph[0].PhotonVector()+diph[1].PhotonVector()).M()/1e3;
//
//        }, {"chosen_two"});
//
//        auto SB = resolved_reco_photons_matched.Filter(
//        [](const double reconstructed_mass)
//        {
//            return (reconstructed_mass < 110) || (reconstructed_mass > 140);
//        }, {"reconstructed_mass"});
//
//        auto SR = resolved_reco_photons_matched.Filter(
//        [](const double reconstructed_mass)
//        {
//            return (reconstructed_mass >= 110) && (reconstructed_mass <= 140);
//        }, {"reconstructed_mass"});
//
//        auto totEventWeight = resolved_reco_photons_matched //rpmi = reco_photons_matched_indices
//        .Define("totEventWeight", [](RVec<int>& x, RVec<unsigned long>& rpmi, RVec<float> photon_id_eff, RVec<float> photon_iso_eff, RVec<float> photon_trg_eff/*, RVec<float> ei_event_weights_generator*/)
//        {
//            //   ||  jic they don't already have the same size  ||
//            //   \/                                             \/
//            auto ResizeVal = std::max({photon_id_eff.size(), photon_iso_eff.size(), photon_trg_eff.size()});
//            photon_id_eff.resize(ResizeVal,1);
//            photon_iso_eff.resize(ResizeVal,1);
//            photon_trg_eff.resize(ResizeVal,1);
//
//            //First, take the x indices from the respective photon_*_eff vectors,
//            //this corresponds to the photons from the set of all reco photons in
//            //the events that passed the "photon_passes_cuts" cuts. Then, from
//            //those photons, select the indices rpmi elements that corresponds to
//            //the two resolved reco-photons
//
//            RVec<float> two_efficiencies = Take(Take(photon_id_eff,x),rpmi) * Take(Take(photon_iso_eff,x),rpmi) * Take(Take(photon_trg_eff,x),rpmi);//*ei_event_weights_generator[0];
//
//            float total = 1.0f;
//
//            for (auto i: two_efficiencies)
//            {
//                total *= i;
//            }
//
//            return total;
//
//        }, {"abstract_photons_pass_cut_indices", "reco_photons_matched_indices", "photon_id_eff", "photon_iso_eff", "photon_trg_eff",/* "ei_event_weights_generator"*/});
//
//        if (counter < 2)
//        {
//            Totals.push_back(EventWeight.Sum<float>("EventWeight"));
//            resultmaps.push_back(VariationsFor(totEventWeight.Sum<float>("totEventWeight")));
//        }
//        else
//        {
//            for (auto& mass_point: massPoints)
//            {
//                auto mass_point_EventWeight = EventWeight.Filter([&]
//                (float axion_mass)
//                {
//                    return (roundToOneDecimalPlace(axion_mass) == mass_point);
////                    return ((axion_mass > mass_point - 0.1f) and (axion_mass < mass_point + 0.1f));
//                }, {"axion_masses"});
//
//                auto mass_point_totEventWeight = totEventWeight.Filter([&]
//                (float axion_mass)
//                {
//                    return (roundToOneDecimalPlace(axion_mass) == mass_point);
////                    return ((axion_mass > mass_point - 0.1f) and (axion_mass < mass_point + 0.1f));
//                }, {"axion_masses"});
//
//                Totals.push_back(mass_point_EventWeight.Sum<float>("EventWeight"));
//                resultmaps.push_back(VariationsFor(mass_point_totEventWeight.Sum<float>("totEventWeight")));
//                Nodes.push_back(mass_point_totEventWeight.Take<UInt_t, RVec<UInt_t>>("ei_event_number"));
//            }
//        }
//
//        counter++;
//
//    }
//
//    ROOT::RDF::RunGraphs(Totals); // running all computation nodes concurrently
//
//    std::cout << "Resolved Category Up\n====================\n";
//    for (int i = 0; i < Nodes.size(); i++)
//    {
//        std::cout << prefixes[i+2] << '\n';
//        std::cout << "==========\n";
//        for (auto& j: *Nodes[i].GetResultPtr<RVec<UInt_t>>())
//        {
//            std::cout << j << ", ";
//        }
//        std::cout << "\n\n";
//    }
//
////EG_RESOLUTION_ALL__1down:           42864.7 ID 1.00524  ISO 1.02283  TRIG 1.04157
////EG_RESOLUTION_ALL__1up:             42751.5 ID 1.00524  ISO 1.02283  TRIG 1.04157
////EG_SCALE_ALL__1down:                42521.1 ID 1.00524  ISO 1.02283  TRIG 1.04157
////EG_SCALE_ALL__1up:                  43059.8 ID 1.00524  ISO 1.02283  TRIG 1.04157
////PH_EFF_ID_Uncertainty__1down:       42790.4 ID 0.992328 ISO 1.02283  TRIG 1.04157
////PH_EFF_ID_Uncertainty__1up:         42790.4 ID 1.01816  ISO 1.02283  TRIG 1.04157
////PH_EFF_ISO_Uncertainty__1down:      42790.4 ID 1.00524  ISO 1.00971  TRIG 1.04157
////PH_EFF_ISO_Uncertainty__1up:        42790.4 ID 1.00524  ISO 1.03595  TRIG 1.04157
////PH_EFF_TRIGGER_Uncertainty__1down:  42790.4 ID 1.00524  ISO 1.02283  TRIG 1.01646
////PH_EFF_TRIGGER_Uncertainty__1up:    42790.4 ID 1.00524  ISO 1.02283  TRIG 1.06668
//
//    ss << R"--(\section*{Table 23 Prompt and Displaced Signal Samples})--" << '\n';
//    ss << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
//    ss << R"--(\setlength\extrarowheight{2pt}\renewcommand{\arraystretch}{1.5})--" << '\n';
//    ss << R"--(\begin{tabular}{|c|c|c|c|c|c|})--" << '\n';
//    ss << R"--(\hline)--" << '\n';
//    ss << R"--(\multicolumn{6}{|c|}{\parbox{\linewidth}{\centering Resolved Photon Category: Up Variations \\ (\% difference from nominal)}}\\[5 pt])--" << '\n';
//    ss << R"--(\hline)--" << '\n';
//
//    ss << R"--({Sample} & EG\_RESOLUTION\_ALL & EG\_SCALE\_ALL & PH\_EFF\_ISO\_Uncertainty & PH\_EFF\_ID\_Uncertainty & PH\_EFF\_TRIGGER\_Uncertainty \\ \hline)--" << '\n';
//    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    cs << R"--(\section*{Table 23 Prompt and Displaced Signal Samples})--" << '\n';
//    cs << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
//    cs << R"--(\setlength\extrarowheight{2pt}\renewcommand{\arraystretch}{1.5})--" << '\n';
//    cs << R"--(\begin{tabular}{|c|c|c|c|c|c|})--" << '\n';
//    cs << R"--(\hline)--" << '\n';
//    cs << R"--(\multicolumn{6}{|c|}{\parbox{\linewidth}{\centering Resolved Photon Category: Up Variations \\ (raw number of events)}}\\[5 pt])--" << '\n';
//    cs << R"--(\hline)--" << '\n';
//
//    cs << R"--({Sample} & EG\_RESOLUTION\_ALL & EG\_SCALE\_ALL & PH\_EFF\_ISO\_Uncertainty & PH\_EFF\_ID\_Uncertainty & PH\_EFF\_TRIGGER\_Uncertainty \\ \hline)--" << '\n';
//
//    for (auto i = 0; (i < resultmaps.size()); i++)
//    {
//        float nominalVal = static_cast<float>(resultmaps[i]["nominal"]);
//
//        ss << prefixes[i] << " & ";
//        ss << std::setprecision(4) << std::fixed
//        << ( (nominalVal && resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1up"]) ? ((resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1up"]-nominalVal)/nominalVal)*100.0 : 0.0)
//        << " & " << std::setprecision(4) << std::fixed
//        << ( (nominalVal && resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1up"]) ? ((resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1up"]-nominalVal)/nominalVal)*100.0 : 0.0)
//        << " & " << std::setprecision(4) << std::fixed
//        << ( (nominalVal && resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1up"]) ? ((resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1up"]-nominalVal)/nominalVal)*100.0 : 0.0)
//        << " & " << std::setprecision(4) << std::fixed
//        << ( (nominalVal && resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1up"]) ? ((resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1up"]-nominalVal)/nominalVal)*100.0 : 0.0)
//        << " & " << std::setprecision(4) << std::fixed
//        << ( (nominalVal && resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up"]) ? ((resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up"]-nominalVal)/nominalVal)*100.0 : 0.0)
//        << R"--( \\ \hline)--" << '\n';
//        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//        cs << prefixes[i] << " & ";
//        cs << std::setprecision(4) << std::fixed
//        << ( (nominalVal && resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1up"]) ? resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1up"] : 0.0)
//        << " & " << std::setprecision(4) << std::fixed
//        << ( (nominalVal && resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1up"]) ? resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1up"] : 0.0)
//        << " & " << std::setprecision(4) << std::fixed
//        << ( (nominalVal && resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1up"]) ? resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1up"] : 0.0)
//        << " & " << std::setprecision(4) << std::fixed
//        << ( (nominalVal && resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1up"]) ? resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1up"] : 0.0)
//        << " & " << std::setprecision(4) << std::fixed
//        << ( (nominalVal && resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up"]) ? resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up"] : 0.0)
//        << R"--( \\ \hline)--" << '\n';
//
//
//        if (!resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1up"])
//        {
//            std::cout << i << ": photon_iso_eff:PH_EFF_ISO_Uncertainty__1up not found\n";
//        }
//        else
//        {
//            std::cout << i << ':' << ' ' << resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1up"] << ' ' << nominalVal << '\n';
//        }
//        if (!resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1up"])
//        {
//            std::cout << i << ": photon_id_eff:PH_EFF_ID_Uncertainty__1up not found\n";
//        }
//        else
//        {
//            std::cout << i << ':' << ' ' << resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1up"] << ' ' << nominalVal << '\n';
//        }
//        if (!resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up"])
//        {
//            std::cout << i << ": photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up not found\n";
//        }
//        else
//        {
//            std::cout << i << ':' << ' ' << resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up"] << ' ' << nominalVal << '\n';
//        }
//    }
//
//    ss << R"--(\end{tabular}})--" << '\n';
//    cs << R"--(\end{tabular}})--" << '\n';
//
//    ss << "\n\n\n";
//    cs << "\n\n\n";
//
//    std::cout << "\n\n\n" << ss.str();
//    std::cout << "\n\n\n" << cs.str();
//}
//
//void Table24_Displaced_Axions()
//{
//    Event::systematics =
//    {
//        "EG_RESOLUTION_ALL__1down",
//        "EG_SCALE_ALL__1down",
//        "PH_EFF_ISO_Uncertainty__1down",
//        "PH_EFF_ID_Uncertainty__1down",
//        "PH_EFF_TRIGGER_Uncertainty__1down",
////        "EG_RESOLUTION_ALL__1down",
//    };
//
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
//    std::vector<float> massPoints = {0.2,0.5,1,3,5,10,10.1,20,20.1,29.5};
//    std::vector<ROOT::RDF::RResultHandle> Totals, Nodes, NoCatNodes;
//    std::vector<RResultMap<float>> resultmaps;
//    std::stringstream ss, cs;
//
//    int counter = 0;
//
//    for (auto& i: input_filenames)
//    {
//        SchottDataFrame df(MakeRDF(i,8));
//
//        auto EventWeight = df.Define("EventWeight",
//        [](const RVec<float>& ei_event_weights_generator)
//        {
//            return  ((ei_event_weights_generator[0]) ? 1 / ei_event_weights_generator[0] : 1);
//
//        }, {"ei_event_weights_generator"}).Define("truth_axions", [&](RVec<TruthParticle> truth_particles)
//        {
//            truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
//            [](TruthParticle& x)
//            {
//             return (std::abs(x.mc_pdg_id) != 36);
//
//            }), truth_particles.end());
//
//            return truth_particles;
//
//        }, {"truth_particles"}).Define("axion_masses", [&](RVec<TruthParticle>& truth_axions)
//        {
//            return truth_axions[0].mc_mass/1e3f;
//
//        }, {"truth_axions"});
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
//                return false;
//            }
//            return true;
//
//        }, {"trigger_passed_triggers"});
//
//        auto two_leptons = trigger_selection.Define("di_electrons",
//        [](RVec<AbstractParticle> electrons)
//        {
//            electrons.erase(std::remove_if(electrons.begin(),electrons.end(),
//            [](AbstractParticle& ep)
//            {
//                //keep electron if it has pt > 20 GeV, |eta| < 2.37, not (1.37 < |eta| < 1.52), and id_medium = 1
//
//                return (!((ep.electron_pt/1e3 > 20) && (std::abs(ep.electron_eta) < 2.37) &&
//                (!((1.37 < std::abs(ep.electron_eta)) && (std::abs(ep.electron_eta) < 1.52)))
//                && (ep.electron_id_medium == 1)));
//
//            }), electrons.end());
//
//            return electrons;
//
//        },{"abstract_electrons"}).Filter([](const RVec<Muon>& muons, RVec<AbstractParticle> electrons)
//        {
//            return (electrons.size()==2 && muons.empty());
//
//        }, {"muons", "di_electrons"});
//
//        auto opp_charge = two_leptons.Filter([](const RVec<AbstractParticle>& electrons)
//        {
//            return (electrons[0].electron_charge*electrons[1].electron_charge < 0);
//
//        }, {"di_electrons"});
//
//        auto leadingPt = opp_charge.Filter([](const RVec<AbstractParticle>& electrons)
//        {
//            return ((electrons[0].electron_pt > 20e3 && electrons[1].electron_pt > 27e3) || (electrons[1].electron_pt > 20e3 && electrons[0].electron_pt > 27e3));
//        }, {"di_electrons"});
//
//        auto deltaR = leadingPt.Filter([] (const RVec<AbstractParticle>& electrons)
//        {
//            return (DeltaR(electrons[0].ElectronVector(), electrons[1].ElectronVector()) > 0.01);
//        }, {"di_electrons"});
//
//        auto mass = deltaR.Filter([] (const RVec<AbstractParticle>& electrons)
//        {
//            auto mass = (electrons[0].ElectronVector() + electrons[1].ElectronVector()).M()/1e3;
//            return ((mass >= 81) && (mass <= 101));
//        }, {"di_electrons"});
//
//        auto ptCut = mass.Filter([] (const RVec<AbstractParticle>& electrons)
//        {
//            auto pT = (electrons[0].ElectronVector() + electrons[1].ElectronVector()).Pt()/1e3;
//            return pT > 10;
//        }, {"di_electrons"});
//
//        auto photon_passes_cuts = ptCut.Define("abstract_photons_pass_cut_indices",[&](RVec<AbstractParticle>& p) //p = photon
//        {
//            RVec<int> x; //indices of photons that pass cuts
//
//            for (auto i = 0; i < p.size(); i++)
//            {
//                if (not ((std::abs(p[i].photon_eta) >= 2.37) or (std::abs(p[i].photon_eta) > 1.37 and std::abs(p[i].photon_eta) < 1.52) or (not p[i].photon_id_loose)))
//                {
//                    x.push_back(i);
//                }
//            }
//
//            return x;
//        }, {"abstract_photons"}).Define("photons_pass_cuts",
//        [&](RVec<AbstractParticle>& photons, RVec<int>& x)
//        {
//            return Take(photons, x); //Taking only the photons that passed the cuts in each event
//
//        }, {"abstract_photons", "abstract_photons_pass_cut_indices"});
//
//        auto resolved_reco_photons_matched = photon_passes_cuts.Define("reco_photons_matched_indices",
//        [](const RVec<AbstractParticle>& reco_photons_matched)
//        {
//            RVec<unsigned long> x; //vector of indices
//            if (reco_photons_matched.size() < 2)
//            {
//                return x;
//            }
//            auto combs = Combinations(reco_photons_matched, 2);
//            size_t length = combs[0].size();
//            double delta_r, m, pt, X, best_X, pt1, pt2, chosen_delta_r;
//
//            for (size_t i = 0; i < length; i++)
//            {
//                delta_r = DeltaR(reco_photons_matched[combs[0][i]].PhotonVector(), reco_photons_matched[combs[1][i]].PhotonVector());
//                m = (reco_photons_matched[combs[0][i]].PhotonVector() + reco_photons_matched[combs[1][i]].PhotonVector()).M();
//                pt = (reco_photons_matched[combs[0][i]].PhotonVector() + reco_photons_matched[combs[1][i]].PhotonVector()).Pt();
//                X = delta_r*(pt/(2.0*m));
//                if (i==0 || ((std::abs(1-X) < std::abs(1-best_X)) and (delta_r < 1.5)))
//                {
//                    best_X = X;
//                    pt1 = reco_photons_matched[combs[0][i]].photon_pt;
//                    pt2 = reco_photons_matched[combs[1][i]].photon_pt;
//                    chosen_delta_r = delta_r;
//                    x = {combs[0][i], combs[1][i]};
//                }
//            }
//            if (pt1 > 10e3 && pt2 > 10e3 && best_X > 0.96 && best_X < 1.2 && chosen_delta_r < 1.5)
//            {
//                return x;
//            }
//
//            x.clear();
//            return x;
//
//        }, {"photons_pass_cuts"}).Filter(
//        [&](RVec<unsigned long>& reco_photons_matched_indices)
//        {
//            return (reco_photons_matched_indices.size()==2);
//
//        }, {"reco_photons_matched_indices"}).Define("chosen_two",
//        [](const RVec<AbstractParticle>& reco_photons_matched, RVec<unsigned long>& reco_photons_matched_indices)
//        {
//            return Take(reco_photons_matched, reco_photons_matched_indices);
//
//        }, {"photons_pass_cuts", "reco_photons_matched_indices"}).Define("reconstructed_mass",
//        [&](RVec<AbstractParticle>& diph)
//        {
//            return (diph[0].PhotonVector()+diph[1].PhotonVector()).M()/1e3;
//
//        }, {"chosen_two"});
//
//        auto SB = resolved_reco_photons_matched.Filter(
//        [](const double reconstructed_mass)
//        {
//            return (reconstructed_mass < 110) || (reconstructed_mass > 140);
//        }, {"reconstructed_mass"});
//
//        auto SR = resolved_reco_photons_matched.Filter(
//        [](const double reconstructed_mass)
//        {
//            return (reconstructed_mass >= 110) && (reconstructed_mass <= 140);
//        }, {"reconstructed_mass"});
//
//        auto totEventWeight = resolved_reco_photons_matched //rpmi = reco_photons_matched_indices
//        .Define("totEventWeight", [](RVec<int>& x, RVec<unsigned long>& rpmi, RVec<float> photon_id_eff, RVec<float> photon_iso_eff, RVec<float> photon_trg_eff/*, RVec<float> ei_event_weights_generator*/)
//        {
//            //   ||  jic they don't already have the same size  ||
//            //   \/                                             \/
//            auto ResizeVal = std::max({photon_id_eff.size(), photon_iso_eff.size(), photon_trg_eff.size()});
//            photon_id_eff.resize(ResizeVal,1);
//            photon_iso_eff.resize(ResizeVal,1);
//            photon_trg_eff.resize(ResizeVal,1);
//
//            //First, take the x indices from the respective photon_*_eff vectors,
//            //this corresponds to the photons from the set of all reco photons in
//            //the events that passed the "photon_passes_cuts" cuts. Then, from
//            //those photons, select the indices rpmi elements that corresponds to
//            //the two resolved reco-photons
//
//            RVec<float> two_efficiencies = Take(Take(photon_id_eff,x),rpmi) * Take(Take(photon_iso_eff,x),rpmi) * Take(Take(photon_trg_eff,x),rpmi);//*ei_event_weights_generator[0];
//
//            float total = 1.0f;
//
//            for (auto i: two_efficiencies)
//            {
//                total *= i;
//            }
//
//            return total;
//
//        }, {"abstract_photons_pass_cut_indices", "reco_photons_matched_indices", "photon_id_eff", "photon_iso_eff", "photon_trg_eff",/* "ei_event_weights_generator"*/});
//
//        if (counter < 2)
//        {
//            Totals.push_back(EventWeight.Sum<float>("EventWeight"));
//            resultmaps.push_back(VariationsFor(totEventWeight.Sum<float>("totEventWeight")));
//            NoCatNodes.push_back(photon_passes_cuts.Count());
//        }
//        else
//        {
//            for (auto& mass_point: massPoints)
//            {
//                auto mass_point_EventWeight = EventWeight.Filter([&]
//                (float axion_mass)
//                {
//                    return (roundToOneDecimalPlace(axion_mass) == mass_point);
////                    return ((axion_mass > mass_point - 0.1f) and (axion_mass < mass_point + 0.1f));
//                }, {"axion_masses"});
//
//                auto mass_point_totEventWeight = totEventWeight.Filter([&]
//                (float axion_mass)
//                {
//                    return (roundToOneDecimalPlace(axion_mass) == mass_point);
////                    return ((axion_mass > mass_point - 0.1f) and (axion_mass < mass_point + 0.1f));
//                }, {"axion_masses"});
//
//                Totals.push_back(mass_point_EventWeight.Sum<float>("EventWeight"));
//                resultmaps.push_back(VariationsFor(mass_point_totEventWeight.Sum<float>("totEventWeight")));
//                Nodes.push_back(mass_point_totEventWeight.Take<UInt_t, RVec<UInt_t>>("ei_event_number"));
//            }
//        }
//
//        counter++;
//
//    }
//
//    ROOT::RDF::RunGraphs(Totals); // running all computation nodes concurrently
//
//    for (int i = 0; i < NoCatNodes.size(); i++)
//    {
//        std::cout << *NoCatNodes[i].GetResultPtr<ULong64_t>() << '\n';
//    }
//
//    std::cout << "Resolved Category Down\n======================\n";
//    for (int i = 0; i < Nodes.size(); i++)
//    {
//        std::cout << prefixes[i+2] << '\n';
//        std::cout << "==========\n";
//        for (auto& j: *Nodes[i].GetResultPtr<RVec<UInt_t>>())
//        {
//            std::cout << j << ", ";
//        }
//        std::cout << "\n\n";
//    }
//
//
//    ss << R"--(\section*{Table 24 Prompt and Displaced Signal Samples})--" << '\n';
//    ss << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
//    ss << R"--(\setlength\extrarowheight{2pt}\renewcommand{\arraystretch}{1.5})--" << '\n';
//    ss << R"--(\begin{tabular}{|c|c|c|c|c|c|})--" << '\n';
//    ss << R"--(\hline)--" << '\n';
//    ss << R"--(\multicolumn{6}{|c|}{\parbox{\linewidth}{\centering Resolved Photon Category: Down Variations \\ (\% difference from nominal)}}\\[5 pt])--" << '\n';
//    ss << R"--(\hline)--" << '\n';
//
//    ss << R"--({Sample} & EG\_RESOLUTION\_ALL & EG\_SCALE\_ALL & PH\_EFF\_ISO\_Uncertainty & PH\_EFF\_ID\_Uncertainty & PH\_EFF\_TRIGGER\_Uncertainty \\ \hline)--" << '\n';
//    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    cs << R"--(\section*{Table 24 Prompt and Displaced Signal Samples})--" << '\n';
//    cs << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
//    cs << R"--(\setlength\extrarowheight{2pt}\renewcommand{\arraystretch}{1.5})--" << '\n';
//    cs << R"--(\begin{tabular}{|c|c|c|c|c|c|})--" << '\n';
//    cs << R"--(\hline)--" << '\n';
//    cs << R"--(\multicolumn{6}{|c|}{\parbox{\linewidth}{\centering Resolved Photon Category: Down Variations \\ (raw number of events)}}\\[5 pt])--" << '\n';
//    cs << R"--(\hline)--" << '\n';
//
//    cs << R"--({Sample} & EG\_RESOLUTION\_ALL & EG\_SCALE\_ALL & PH\_EFF\_ISO\_Uncertainty & PH\_EFF\_ID\_Uncertainty & PH\_EFF\_TRIGGER\_Uncertainty \\ \hline)--" << '\n';
//
//    for (auto i = 0; (i < resultmaps.size()); i++)
//    {
//        float nominalVal = static_cast<float>(resultmaps[i]["nominal"]);
//
//        ss << prefixes[i] << " & ";
//        ss << std::setprecision(4) << std::fixed
//        << ( (nominalVal && resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1down"]) ? ((resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1down"]-nominalVal)/nominalVal)*100.0 : 0.0)
//        << " & " << std::setprecision(4) << std::fixed
//        << ( (nominalVal && resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1down"]) ? ((resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1down"]-nominalVal)/nominalVal)*100.0 : 0.0)
//        << " & " << std::setprecision(4) << std::fixed
//        << ( (nominalVal && resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1down"]) ? ((resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1down"]-nominalVal)/nominalVal)*100.0 : 0.0)
//        << " & " << std::setprecision(4) << std::fixed
//        << ( (nominalVal && resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1down"]) ? ((resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1down"]-nominalVal)/nominalVal)*100.0 : 0.0)
//        << " & " << std::setprecision(4) << std::fixed
//        << ( (nominalVal && resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down"]) ? ((resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down"]-nominalVal)/nominalVal)*100.0 : 0.0)
//        << R"--( \\ \hline)--" << '\n';
//        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//        cs << prefixes[i] << " & ";
//        cs << std::setprecision(4) << std::fixed
//        << ( (nominalVal && resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1down"]) ? resultmaps[i]["photons_and_electrons:EG_RESOLUTION_ALL__1down"] : 0.0)
//        << " & " << std::setprecision(4) << std::fixed
//        << ( (nominalVal && resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1down"]) ? resultmaps[i]["photons_and_electrons:EG_SCALE_ALL__1down"] : 0.0)
//        << " & " << std::setprecision(4) << std::fixed
//        << ( (nominalVal && resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1down"]) ? resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1down"] : 0.0)
//        << " & " << std::setprecision(4) << std::fixed
//        << ( (nominalVal && resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1down"]) ? resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1down"] : 0.0)
//        << " & " << std::setprecision(4) << std::fixed
//        << ( (nominalVal && resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down"]) ? resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down"] : 0.0)
//        << R"--( \\ \hline)--" << '\n';
//
//
//        if (!resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1down"])
//        {
//            std::cout << i << ": photon_iso_eff:PH_EFF_ISO_Uncertainty__1down not found\n";
//        }
//        else
//        {
//            std::cout << i << ':' << ' ' << resultmaps[i]["photon_iso_eff:PH_EFF_ISO_Uncertainty__1down"] << ' ' << nominalVal << '\n';
//        }
//        if (!resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1down"])
//        {
//            std::cout << i << ": photon_id_eff:PH_EFF_ID_Uncertainty__1down not found\n";
//        }
//        else
//        {
//            std::cout << i << ':' << ' ' << resultmaps[i]["photon_id_eff:PH_EFF_ID_Uncertainty__1down"] << ' ' << nominalVal << '\n';
//        }
//        if (!resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down"])
//        {
//            std::cout << i << ": photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down not found\n";
//        }
//        else
//        {
//            std::cout << i << ':' << ' ' << resultmaps[i]["photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down"] << ' ' << nominalVal << '\n';
//        }
//    }
//
//    ss << R"--(\end{tabular}})--" << '\n';
//    cs << R"--(\end{tabular}})--" << '\n';
//
//    ss << "\n\n\n";
//    cs << "\n\n\n";
//
//    std::cout << "\n\n\n" << ss.str();
//    std::cout << "\n\n\n" << cs.str();
//}

void Section7_TablesPlots()
{
    
    auto start_time = Clock::now();
//    auto verbosity = ROOT::Experimental::RLogScopedVerbosity(ROOT::Detail::RDF::RDFLogChannel(), ROOT::Experimental::ELogLevel::kInfo);
//    Table21();
//    Table22();
//    Table23();
//    Table24();
    Merged_Up_Higgs_Mass();
//    Merged_Down_Higgs_Mass();
//    Resolved_Up_Higgs_Mass();
//    Resolved_Down_Higgs_Mass();
//    Table21_Displaced_Axions();
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



