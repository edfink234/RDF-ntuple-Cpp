#include <iostream>
#include <vector>
#include <chrono>
#include <algorithm>
#include <string>
#include <sstream>
#include <iomanip>
#include <array>
#include <cstdlib>
#include <cmath>
#include <unordered_set>
#include <unordered_map>
#include <map>

#include <ROOT/RLogger.hxx>
#include "Math/VectorUtil.h"
#include "ROOT/RDataFrame.hxx"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "Rtypes.h"
#include "ROOT/RDFHelpers.hxx"
#include "TAttAxis.h"

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

static std::unordered_map<float, float> a, c, e, g;
static std::unordered_map<float, std::unordered_map<float, float>> b, d, f, h;

float roundToOneDecimalPlace(float num) {
    std::stringstream stream;
    stream << std::fixed << std::setprecision(1) << num;
    float rounded_num = std::stof(stream.str());
    return rounded_num;
}

//Calculates the uncertainty in the ratio ε_p / ε_ll
double unc(double eff_p, double eff_ll, double N_p, double N_ll)
{
    return (eff_ll / eff_p) * sqrt(((1 - eff_p)/(N_p * eff_p)) + ((1 - eff_ll)/(N_ll * eff_ll)));
}

//void Table3()
//{
//    std::vector<std::vector<std::string>> input_filenames = {
//        //Z gamma background
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617070._000001.LGNTuple.root"},
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617064._000001.LGNTuple.root"},
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617074._000001.LGNTuple.root"},
//        //Signal
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root"}, // 1 GeV
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root"}, // 5 GeV
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600751.PhPy8EG_AZNLO_ggH125_mA2p0_v1.root"}, // 2 GeV
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600752.PhPy8EG_AZNLO_ggH125_mA3p0_v1.root"}, // 3 GeV
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600756.PhPy8EG_AZNLO_ggH125_mA9p0_v1.root"}, // 9 GeV
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
////    Event::systematics =
////    {
//////        "PH_EFF_ISO_Uncertainty",
//////        "PH_EFF_ISO_Uncertainty",
//////        "EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR",
//////        "PRW_DATASF",
//////        "MUON_EFF_RECO_SYS",
//////        "MUON_EFF_ISO_SYS",
//////        "MUON_EFF_TrigSystUncertainty",
//////        "EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR",
//////        "MUON_EFF_TrigStatUncertainty",
//////        "MUON_EFF_RECO_STAT",
//////        "MUON_EFF_TTVA_STAT",
//////        "EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR",
//////        "EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR",
//////        "MUON_EFF_TTVA_SYS",
//////        "MUON_EFF_ISO_STAT",
//////        "MUON_SAGITTA_RHO",
//////        "EG_RESOLUTION_ALL",
//////        "EG_SCALE_ALL",
//////        "MUON_MS",
//////        "MUON_ID",
//////        "EL_EFF_TriggerEff_TOTAL_1NPCOR_PLUS_UNCOR",
//////        "MUON_SAGITTA_RESBIAS",
//////        "MUON_SCALE",
////
//////        "PH_EFF_ISO_Uncertainty__1down",
//////        "PH_EFF_ID_Uncertainty__1up",
//////        "PH_EFF_TRIGGER_Uncertainty__1up",
////        "EG_SCALE_ALL__1down",
////        "EG_SCALE_ALL__1up",
//////        "PH_EFF_ID_Uncertainty__1down",
//////        "PH_EFF_TRIGGER_Uncertainty__1down",
//////        "PH_EFF_ISO_Uncertainty__1up",
////        "EG_RESOLUTION_ALL__1up",
////        "EG_RESOLUTION_ALL__1down",
////    };
//
//    std::array<double,9> JetNumeratorSFs = {((139e15)*(1.9828e-9)*(0.821204)),((139e15)*(110.64e-12)*(0.69275)),((139e15)*(40.645e-12)*(0.615906)),((139e15)*(1.9817e-9)*(0.1136684)),((139e15)*(110.47e-12)*(0.1912956)),((139e15)*(40.674e-12)*(0.2326772)),((139e15)*(1.9819e-9)*(0.0656969)),((139e15)*(110.53e-12)*(0.1158741)),((139e15)*(40.68e-12)*(0.1535215))}; //numerators for jet bkg
//    std::array<double,3> SFs = {((139e15)*(.871e-12)),((139e15)*(.199e-12)), ((139e15)*(.0345e-15))}; //numerators for Z-gamma bkg
//
//    constexpr std::array<const char*,18> Samples = {R"--(pty2\_9\_17)--", R"--(pty\_17\_myy\_0\_80)--", R"--(pty\_17\_myy\_80)--", R"--(Signal $m_{\text{A}}$ = 1 GeV)--", R"--(Signal $m_{\text{A}}$ = 5 GeV)--", R"--(Signal $m_{\text{A}}$ = 2 GeV)--", R"--(Signal $m_{\text{A}}$ = 3 GeV)--", R"--(Signal $m_{\text{A}}$ = 9 GeV)--", "Data", R"--(Zee\_lightJet\_0-70)--", R"--(Zee\_lightJet\_70-140)--", R"--(Zee\_lightJet\_140-280)--", R"--(Zee\_cJet\_0-70)--", R"--(Zee\_cJet\_70-140)--", R"--(Zee\_cJet\_140-280)--", R"--(Zee\_bJet\_0-70)--", R"--(Zee\_bJet\_70-140)--", R"--(Zee\_bJet\_140-280)--"};
//
//    std::ostringstream os;
//    os << R"--(\section*{Table 3})--" << '\n';
//    os << R"--(\hspace{-3cm}\scalebox{0.45}{)--" << '\n';
//    os << R"--(\setlength\extrarowheight{2pt})--" << '\n';
//    os << R"--(\renewcommand{\arraystretch}{1.5})--" << '\n';
//    os << R"--(\begin{tabular}{|c|c|c|c|c|c|c|c|})--" << '\n';
//    os << R"--(\hline)--" << '\n';
//    os << R"--(\textbf{Sample} & \textbf{Before Preselection} & \textbf{2 leptons}
//          & \textbf{Opposite Charge} & $\pmb{p_{T}}^{\textbf{leading}}\pmb{>}\text{27}$ \textbf{GeV}, \; $\pmb{p_{T}}^{\textbf{sub-leading}}\pmb{>}\textbf{20}$ \textbf{GeV} & \textbf{Same flavour} & \textbf{dilep mass cut} & \textbf{dilep} $\pmb{p_{T}}$ \textbf{cut} \\ \hline )--" << '\n';
//
//    double beforePreselecZGamma = 0, twoLeptonsZGamma = 0, oppChargeZGamma = 0, leadingPtZGamma = 0, deltaRZGamma = 0, MassZGamma = 0, ptCutZGamma = 0;
//    double beforePreselecZGammaStatUnc = 0, twoLeptonsZGammaStatUnc = 0, oppChargeZGammaStatUnc = 0, leadingPtZGammaStatUnc = 0, deltaRZGammaStatUnc = 0, MassZGammaStatUnc = 0, ptCutZGammaStatUnc = 0;
//
//    double beforePreselecZJets = 0, twoLeptonsZJets = 0, oppChargeZJets = 0, leadingPtZJets = 0, deltaRZJets = 0, MassZJets = 0, ptCutZJets = 0;
//    double beforePreselecZJetsStatUnc = 0, twoLeptonsZJetsStatUnc = 0, oppChargeZJetsStatUnc = 0, leadingPtZJetsStatUnc = 0, deltaRZJetsStatUnc = 0, MassZJetsStatUnc = 0, ptCutZJetsStatUnc = 0;
//
//    std::vector<ROOT::RDF::RResultHandle> Nodes;
////    std::vector<ROOT::RDF::RResultHandle> tempNodes;
////
////    std::vector<RResultMap<ULong64_t>> resultmaps;
//
//    for (auto& file: input_filenames)
//    {
//        SchottDataFrame df(MakeRDF(file, 8));
//
////        df.Describe().Print();
////        exit(1);
//
//        auto trigger_selection = df.Filter(
//        [](const RVec<std::string>& trigger_passed_triggers)
//        {
//            bool trigger_found = (std::find_first_of(trigger_passed_triggers.begin(), trigger_passed_triggers.end(), triggers.begin(), triggers.end()) != trigger_passed_triggers.end());
//
//            if (!trigger_found)
//            {
//              return false;
//            }
//            return true;
//
//        }, {"trigger_passed_triggers"});
//
//        auto two_leptons = df.Define("di_electrons",
//        [](RVec<Electron> electrons)
//        {
//            electrons.erase(std::remove_if(electrons.begin(),electrons.end(),
//            [](Electron& ep)
//            {
//                return (!((ep.electron_pt/1e3 > 20) && (std::abs(ep.electron_eta) < 2.37) &&
//                                          (!((1.37 < std::abs(ep.electron_eta)) && (std::abs(ep.electron_eta) < 1.52)))
//                                          && (ep.electron_id_medium == 1)));
//
//            }), electrons.end());
//
//            return electrons;
//
//        },{"electrons"}).Filter([](RVec<Muon>& muons, RVec<Electron> di_electrons)
//        {
//            return (di_electrons.size()==2 && muons.empty() && DeltaR(di_electrons[0].Vector(), di_electrons[1].Vector()) > 0.01);
//
//        }, {"muons", "di_electrons"});
//
//        auto opp_charge = two_leptons.Filter([](RVec<Electron> di_electrons)
//        {
//            return (di_electrons[0].electron_charge*di_electrons[1].electron_charge < 0);
//
//        }, {"di_electrons"});
//
//        auto leading_pt = opp_charge.Filter([](RVec<Electron>& electrons)
//        {
//            return ((electrons[0].electron_pt >= 20e3 && electrons[1].electron_pt >= 27e3) || (electrons[1].electron_pt >= 20e3 && electrons[0].electron_pt >= 27e3));
//        }, {"di_electrons"});
//
////        auto delta_R = leading_pt.Filter([] (RVec<Electron>& electrons)
////        {
////            return (DeltaR(electrons[0].Vector(), electrons[1].Vector()) > 0.01);
////        }, {"di_electrons"});
//
//        auto same_flavour = leading_pt.Filter([] (RVec<Electron>& electrons)
//        {
//            return true; //std::abs(electrons[0].electron_pdg_id) == std::abs(electrons[1].electron_pdg_id) == 11;
//        }, {"di_electrons"});
//
//        auto dilep = same_flavour.Define("dilep",[] (RVec<Electron>& electrons)
//        {
//            return (electrons[0].Vector() + electrons[1].Vector());
//        }, {"di_electrons"});
//
//        auto mass = dilep.Filter([] (PtEtaPhiEVector& dilep)
//        {
//            auto mass = dilep.M()/1e3;
//            return ((mass >= 81) && (mass <= 101));
//        }, {"dilep"});
//
//        auto pt_cut = mass.Filter([] (PtEtaPhiEVector& dilep)
//        {
//            auto pT = dilep.Pt()/1e3;
//            return pT > 10;
//        }, {"dilep"});
//
//        Nodes.push_back(df.Count());
//        Nodes.push_back(two_leptons.Count());
//        Nodes.push_back(opp_charge.Count());
//        Nodes.push_back(leading_pt.Count());
//        Nodes.push_back(same_flavour.Count());
//        Nodes.push_back(mass.Count());
//        Nodes.push_back(pt_cut.Count());
//
////        resultmaps.push_back(VariationsFor(df.Count()));
////        resultmaps.push_back(VariationsFor(two_leptons.Count()));
////        resultmaps.push_back(VariationsFor(opp_charge.Count()));
////        resultmaps.push_back(VariationsFor(leading_pt.Count()));
////        resultmaps.push_back(VariationsFor(same_flavour.Count()));
////        resultmaps.push_back(VariationsFor(mass.Count()));
////        resultmaps.push_back(VariationsFor(pt_cut.Count()));
////
////        tempNodes.push_back(df.Take<RVec<std::vector<std::string>>, RVec<RVec<std::vector<std::string>>>>("electron_syst_name"));
////        tempNodes.push_back(df.Take<RVec<std::vector<std::string>>, RVec<RVec<std::vector<std::string>>>>("photon_syst_name"));
//    }
//
//    constexpr std::array<const char*,7> Cuts = {"total", "two leptons", "opposite charge", "leading pt", "same flavour", "mass", "pt cut"};
//
//    ROOT::RDF::RunGraphs(Nodes); // running all computation nodes concurrently
//
////    0     1     2     3     4     5     6       Z-gamma
////    7     8     9     10    11    12    13      Z-gamma
////    14    15    16    17    18    19    20      Z-gamma
////    21    22    23    24    25    26    27      ma1
////    28    29    30    31    32    33    34      ma5
////    35    36    37    38    39    40    41      ma2
////    42    43    44    45    46    47    48      ma3
////    49    50    51    52    53    54    55      ma9
////    56    57    58    59    60    61    62      data
////    63    64    65    66    67    68    69      Z-jets
////    70    71    72    73    74    75    76      Z-jets
////    77    78    79    80    81    82    83      Z-jets
////    84    85    86    87    88    89    90      Z-jets
////    91    92    93    94    95    96    97      Z-jets
////    98    99    100   101   102   103   104     Z-jets
////    105   106   107   108   109   110   111     Z-jets
////    112   113   114   115   116   117   118     Z-jets
////    119   120   121   122   123   124   125     Z-jets
//
////    std::unordered_set<std::string> uniqueSystematics;
////    std::unordered_set<std::string> ZGammaSystematics;
////    std::unordered_set<std::string> SignalSystematics;
////    std::unordered_set<std::string> DataSystematics;
////    std::unordered_set<std::string> ZJetsSystematics;
////
////    for (auto& i: tempNodes)
////    {
////        for (auto& j: *i.GetResultPtr<RVec<RVec<std::vector<std::string>>>>())
////        {
////            for (auto& k: j)
////            {
////                for (auto& l: k)
////                {
////                    uniqueSystematics.insert(l);
////                }
////            }
////        }
////    }
////
////    int counter = 0;
////    for (auto& i: tempNodes)
////    {
////        for (auto& j: (*i.GetResultPtr<RVec<RVec<std::vector<std::string>>>>())[0])
////        {
////            for (auto& k: j)
////            {
////                if (counter >= 0 && counter <= 5) //Z-gamma
////                {
////                    ZGammaSystematics.insert(k);
////                }
////
////                else if (counter >= 6 && counter <= 9) //Signal
////                {
////                    SignalSystematics.insert(k);
////                }
////
////                else if (counter == 10) //Data
////                {
////                    DataSystematics.insert(k);
////                }
////
////                else
////                {
////                    ZJetsSystematics.insert(k);
////                }
////            }
////        }
////        counter++;
////    }
//////
////    std::cout << "ZGammaSystematics\n=================\n";
////    for (auto& i: ZGammaSystematics)
////    {
////        std::cout << i << '\n';
////    }
////    std::cout << "\n\n\n\n";
////
////    std::cout << "SignalSystematics\n=================\n";
////    for (auto& i: SignalSystematics)
////    {
////        std::cout << i << '\n';
////    }
////    std::cout << "\n\n\n\n";
////
////    std::cout << "DataSystematics\n===============\n";
////    for (auto& i: DataSystematics)
////    {
////        std::cout << i << '\n';
////    }
////    std::cout << "\n\n\n\n";
////
////    std::cout << "ZJetsSystematics\n================\n";
////    for (auto& i: ZJetsSystematics)
////    {
////        std::cout << i << '\n';
////    }
////    std::cout << "\n\n\n\n";
//
////
////    std::cout << "\n\n";
////    std::cout << resultmaps.size() << '\n';
////    for (auto i = 0; i < resultmaps.size(); i++)
////    {
////        if (i % 7 == 0)
////        {
////            std::cout << Samples[i/7] << "\n============================\n\n";
////        }
////        std::cout << Cuts[i%7] << "\n===============\n";
////        for (auto& var: resultmaps[i].GetKeys())
////        {
////            std::cout << std::setw(44) << var <<
////            std::setw(44) << resultmaps[i][var] << '\n';
////        }
////        std::cout << '\n';
////    }
//
//    std::cout << R"--(\section*{Table 3 Signal Ratios})--" << '\n';
//    std::cout << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
//    std::cout << R"--(\begin{tabular}{|c|c|c|c|c|c|c|})--" << '\n';
//    std::cout << R"--(\hline)--" << '\n';
//    std::cout << R"--(Sample & $\frac{\text{2 leptons}}{\text{2 leptons}}$
//          & $\frac{\text{Opposite Charge}}{\text{2 leptons}}$ & $\frac{p_{T}^{\text{leading}} > 27\text{ GeV, } \, p_{T}^{\text{sub-leading}} > 20 \text{ GeV }}{\text{2 leptons}}$ & $\frac{\text{Same flavour}}{\text{2 leptons}}$ & $\frac{\text{dilep mass cut}}{\text{2 leptons}}$ & $\frac{\text{dilep }p_{T} \text{ cut}}{\text{2 leptons}}$ \\ \hline )--" << '\n';
//
//    os.setf(std::ios::fixed);
//    os.precision(2);
//
//    for (int i=0, j=0; (i<18 && j <= 119); i++, j+=7)
//    {
//        os << Samples[i];
//        if (i >= 0 && i <= 2) //Z-gamma
//        {
//            os  << " & " << *Nodes[j].GetResultPtr<ULong64_t>() * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j].GetResultPtr<ULong64_t>()) * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << " & " << *Nodes[j+1].GetResultPtr<ULong64_t>() * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+1].GetResultPtr<ULong64_t>()) * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << " & " << *Nodes[j+2].GetResultPtr<ULong64_t>() * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+2].GetResultPtr<ULong64_t>()) * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << " & " << *Nodes[j+3].GetResultPtr<ULong64_t>() * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+3].GetResultPtr<ULong64_t>()) * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << " & " << *Nodes[j+4].GetResultPtr<ULong64_t>() * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+4].GetResultPtr<ULong64_t>()) * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << " & " << *Nodes[j+5].GetResultPtr<ULong64_t>() * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+5].GetResultPtr<ULong64_t>()) * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << " & " << *Nodes[j+6].GetResultPtr<ULong64_t>() * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+6].GetResultPtr<ULong64_t>()) * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << R"--( \\ \hline )--" << '\n';
//        }
//
//        else if (i >= 9) //Z-jets
//        {
//            os  << " & " << *Nodes[j].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << " & " << *Nodes[j+1].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+1].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << " & " << *Nodes[j+2].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+2].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << " & " << *Nodes[j+3].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+3].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << " & " << *Nodes[j+4].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+4].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << " & " << *Nodes[j+5].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+5].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << " & " << *Nodes[j+6].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+6].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << R"--( \\ \hline )--" << '\n';
//        }
//        else
//        {
//            os  << " & " << *Nodes[j].GetResultPtr<ULong64_t>()
//            << R"--($\, \pm \,$)--" << sqrt(*Nodes[j].GetResultPtr<ULong64_t>())
//            << " & " << *Nodes[j+1].GetResultPtr<ULong64_t>()
//            << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+1].GetResultPtr<ULong64_t>())
//            << " & " << *Nodes[j+2].GetResultPtr<ULong64_t>()
//            << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+2].GetResultPtr<ULong64_t>())
//            << " & " << *Nodes[j+3].GetResultPtr<ULong64_t>()
//            << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+3].GetResultPtr<ULong64_t>())
//            << " & " << *Nodes[j+4].GetResultPtr<ULong64_t>()
//            << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+4].GetResultPtr<ULong64_t>())
//            << " & " << *Nodes[j+5].GetResultPtr<ULong64_t>()
//            << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+5].GetResultPtr<ULong64_t>())
//            << " & " << *Nodes[j+6].GetResultPtr<ULong64_t>()
//            << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+6].GetResultPtr<ULong64_t>())
//            << R"--( \\ \hline )--" << '\n';
//
//            if (i==3) //1 GeV
//            {
//                std::cout << Samples[i] << " (Me) & " << 1 << " & " <<
//                static_cast<double>(*Nodes[j+2].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
//                << " & " <<
//                static_cast<double>(*Nodes[j+3].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
//                << " & " <<
//                static_cast<double>(*Nodes[j+4].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
//                << " & " <<
//                static_cast<double>(*Nodes[j+5].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
//                << " & " <<
//                static_cast<double>(*Nodes[j+6].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
//                << R"--( \\ \hline )--" << '\n';
//
//                std::cout << Samples[i] << " (Paper) & " << 1 << " & " <<
//                21505.03 / 21606.75
//                << " & " <<
//                21375.48 / 21606.75
//                << " & " <<
//                21375.33 / 21606.75
//                << " & " <<
//                20543.06 / 21606.75
//                << " & " <<
//                19516.87 / 21606.75
//                << R"--( \\ \hline )--" << '\n';
//            }
//
//            if (i==4) //5 GeV
//            {
//                std::cout << Samples[i] << " (Me) & " << 1 << " & " <<
//                static_cast<double>(*Nodes[j+2].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
//                << " & " <<
//                static_cast<double>(*Nodes[j+3].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
//                << " & " <<
//                static_cast<double>(*Nodes[j+4].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
//                << " & " <<
//                static_cast<double>(*Nodes[j+5].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
//                << " & " <<
//                static_cast<double>(*Nodes[j+6].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
//                << R"--( \\ \hline )--" << '\n';
//
//                std::cout << Samples[i] << " (Paper) & " << 1 << " & " <<
//                21556.03 / 21655.38
//                << " & " <<
//                21424.29 / 21655.38
//                << " & " <<
//                21424.28 / 21655.38
//                << " & " <<
//                20585.09 / 21655.38
//                << " & " <<
//                19536.68 / 21655.38
//                << R"--( \\ \hline )--" << '\n';
//            }
//        }
//    }
//
//    std::cout << R"--(\end{tabular}})--" << "\n\n\n";
//
//    for (int i = 0, j = 0; (i <= 14 && j <= 2); i += 7, j++)
//    {
//        beforePreselecZGamma += *Nodes[i].GetResultPtr<ULong64_t>() * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
//        twoLeptonsZGamma += *Nodes[i+1].GetResultPtr<ULong64_t>() * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
//        oppChargeZGamma += *Nodes[i+2].GetResultPtr<ULong64_t>() * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
//        leadingPtZGamma += *Nodes[i+3].GetResultPtr<ULong64_t>() * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
//        deltaRZGamma += *Nodes[i+4].GetResultPtr<ULong64_t>() * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
//        MassZGamma += *Nodes[i+5].GetResultPtr<ULong64_t>() * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
//        ptCutZGamma += *Nodes[i+6].GetResultPtr<ULong64_t>() * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
//
//        beforePreselecZGammaStatUnc += pow(sqrt(*Nodes[i].GetResultPtr<ULong64_t>()) * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
//        twoLeptonsZGammaStatUnc += pow(sqrt(*Nodes[i+1].GetResultPtr<ULong64_t>()) * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
//        oppChargeZGammaStatUnc += pow(sqrt(*Nodes[i+2].GetResultPtr<ULong64_t>()) * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
//        leadingPtZGammaStatUnc += pow(sqrt(*Nodes[i+3].GetResultPtr<ULong64_t>()) * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
//        deltaRZGammaStatUnc += pow(sqrt(*Nodes[i+4].GetResultPtr<ULong64_t>()) * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
//        MassZGammaStatUnc += pow(sqrt(*Nodes[i+5].GetResultPtr<ULong64_t>()) * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
//        ptCutZGammaStatUnc += pow(sqrt(*Nodes[i+6].GetResultPtr<ULong64_t>()) * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
//    }
//
//    for (int i = 63, j = 0; (i <= 119 && j <= 8); i += 7, j++)
//    {
//        beforePreselecZJets += *Nodes[i].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
//        twoLeptonsZJets += *Nodes[i+1].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
//        oppChargeZJets += *Nodes[i+2].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
//        leadingPtZJets += *Nodes[i+3].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
//        deltaRZJets += *Nodes[i+4].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
//        MassZJets += *Nodes[i+5].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
//        ptCutZJets += *Nodes[i+6].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
//
//        beforePreselecZJetsStatUnc += pow(sqrt(*Nodes[i].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
//        twoLeptonsZJetsStatUnc += pow(sqrt(*Nodes[i+1].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
//        oppChargeZJetsStatUnc += pow(sqrt(*Nodes[i+2].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
//        leadingPtZJetsStatUnc += pow(sqrt(*Nodes[i+3].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
//        deltaRZJetsStatUnc += pow(sqrt(*Nodes[i+4].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
//        MassZJetsStatUnc += pow(sqrt(*Nodes[i+5].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
//        ptCutZJetsStatUnc += pow(sqrt(*Nodes[i+6].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
//    }
//
//    os << R"--(Total $Z\gamma$ & )--"
//    << beforePreselecZGamma << R"--($\, \pm \,$)--" << sqrt(beforePreselecZGammaStatUnc) << " & "
//    << twoLeptonsZGamma << R"--($\, \pm \,$)--" << sqrt(twoLeptonsZGammaStatUnc) << " & "
//    << oppChargeZGamma << R"--($\, \pm \,$)--" << sqrt(oppChargeZGammaStatUnc) << " & "
//    << leadingPtZGamma << R"--($\, \pm \,$)--" << sqrt(leadingPtZGammaStatUnc) << " & "
//    << deltaRZGamma << R"--($\, \pm \,$)--" << sqrt(deltaRZGammaStatUnc) << " & "
//    << MassZGamma << R"--($\, \pm \,$)--" << sqrt(MassZGammaStatUnc) << " & "
//    << ptCutZGamma << R"--($\, \pm \,$)--" << sqrt(ptCutZGammaStatUnc)
//    << R"--( \\ \hline )--" << '\n';
//
//    os << R"--(Total $Z$+jets & )--"
//    << beforePreselecZJets << R"--($\, \pm \,$)--" << sqrt(beforePreselecZJetsStatUnc) << " & "
//    << twoLeptonsZJets << R"--($\, \pm \,$)--" << sqrt(twoLeptonsZJetsStatUnc) << " & "
//    << oppChargeZJets << R"--($\, \pm \,$)--" << sqrt(oppChargeZJetsStatUnc) << " & "
//    << leadingPtZJets << R"--($\, \pm \,$)--" << sqrt(leadingPtZJetsStatUnc) << " & "
//    << deltaRZJets << R"--($\, \pm \,$)--" << sqrt(deltaRZJetsStatUnc) << " & "
//    << MassZJets << R"--($\, \pm \,$)--" << sqrt(MassZJetsStatUnc) << " & "
//    << ptCutZJets << R"--($\, \pm \,$)--" << sqrt(ptCutZJetsStatUnc)
//    << R"--( \\ \hline )--" << '\n';
//
//    os << R"--(Total Bkg & )--"
//    << beforePreselecZJets+beforePreselecZGamma << R"--($\, \pm \,$)--" <<
//    sqrt(beforePreselecZJetsStatUnc+beforePreselecZGammaStatUnc) << " & "
//    << twoLeptonsZJets+twoLeptonsZGamma << R"--($\, \pm \,$)--" <<
//    sqrt(twoLeptonsZJetsStatUnc+twoLeptonsZGammaStatUnc) << " & "
//    << oppChargeZJets+oppChargeZGamma << R"--($\, \pm \,$)--" <<
//    sqrt(oppChargeZJetsStatUnc+oppChargeZGammaStatUnc) << " & "
//    << leadingPtZJets+leadingPtZGamma << R"--($\, \pm \,$)--" <<
//    sqrt(leadingPtZJetsStatUnc+leadingPtZGammaStatUnc) << " & "
//    << deltaRZJets+deltaRZGamma << R"--($\, \pm \,$)--" <<
//    sqrt(deltaRZJetsStatUnc+deltaRZGammaStatUnc) << " & "
//    << MassZJets+MassZGamma << R"--($\, \pm \,$)--" <<
//    sqrt(MassZJetsStatUnc+MassZGammaStatUnc) << " & "
//    << ptCutZJets+ptCutZGamma << R"--($\, \pm \,$)--" <<
//    sqrt(ptCutZJetsStatUnc+ptCutZGammaStatUnc) << R"--( \\ \hline )--" << '\n';
//
//    os << R"--(\end{tabular}})--" << '\n';
//    std::cout << os.str() << '\n';
//}
//
//void Table8()
//{
//    std::vector<std::vector<std::string>> input_filenames =
//    {
//        //Z gamma background
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617070._000001.LGNTuple.root"},
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617064._000001.LGNTuple.root"},
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617074._000001.LGNTuple.root"},
//        //Signal
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root"}, // 1 GeV
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root"}, // 5 GeV
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600751.PhPy8EG_AZNLO_ggH125_mA2p0_v1.root"}, // 2 GeV
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600752.PhPy8EG_AZNLO_ggH125_mA3p0_v1.root"}, // 3 GeV
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600756.PhPy8EG_AZNLO_ggH125_mA9p0_v1.root"}, // 9 GeV
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
//    std::array<double,9> JetNumeratorSFs = {((139e15)*(1.9828e-9)*(0.821204)),((139e15)*(110.64e-12)*(0.69275)),((139e15)*(40.645e-12)*(0.615906)),((139e15)*(1.9817e-9)*(0.1136684)),((139e15)*(110.47e-12)*(0.1912956)),((139e15)*(40.674e-12)*(0.2326772)),((139e15)*(1.9819e-9)*(0.0656969)),((139e15)*(110.53e-12)*(0.1158741)),((139e15)*(40.68e-12)*(0.1535215))}; //numerators for jet bkg
//    std::array<double,3> SFs = {((139e15)*(.871e-12)),((139e15)*(.199e-12)), ((139e15)*(.0345e-15))}; //numerators for Z-gamma bkg
//
//    constexpr std::array<const char*,18> Samples = {R"--(pty2\_9\_17)--", R"--(pty\_17\_myy\_0\_80)--", R"--(pty\_17\_myy\_80)--", R"--(Signal $m_{\text{A}}$ = 1 GeV)--", R"--(Signal $m_{\text{A}}$ = 5 GeV)--", R"--(Signal $m_{\text{A}}$ = 2 GeV)--", R"--(Signal $m_{\text{A}}$ = 3 GeV)--", R"--(Signal $m_{\text{A}}$ = 9 GeV)--", "Data", R"--(Zee\_lightJet\_0-70)--", R"--(Zee\_lightJet\_70-140)--", R"--(Zee\_lightJet\_140-280)--", R"--(Zee\_cJet\_0-70)--", R"--(Zee\_cJet\_70-140)--", R"--(Zee\_cJet\_140-280)--", R"--(Zee\_bJet\_0-70)--", R"--(Zee\_bJet\_70-140)--", R"--(Zee\_bJet\_140-280)--"};
//
//    double totalEventsZgamma = 0, resolvedEventsZgamma = 0, SREventsZgamma = 0, SBEventsZgamma = 0;
//    double totalEventsZgammaStatUnc = 0, resolvedEventsZgammaStatUnc = 0, SREventsZgammaStatUnc = 0, SBEventsZgammaStatUnc = 0;
//    double totalEventsZJets = 0, resolvedEventsZJets = 0, SREventsZJets = 0, SBEventsZJets = 0;
//    double totalEventsZJetsStatUnc = 0, resolvedEventsZJetsStatUnc = 0, SREventsZJetsStatUnc = 0, SBEventsZJetsStatUnc = 0;
//
//    std::vector<ROOT::RDF::RResultHandle> Nodes;
//
//    std::ostringstream os;
//    os << R"--(\section*{Table 8})--" << '\n';
//    os << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
//    os << R"--(\setlength\extrarowheight{2pt})--" << '\n';
//    os << R"--(\renewcommand{\arraystretch}{1.5})--" << '\n';
//    os << R"--(\begin{tabular}{|c|c|c|c|c|})--" << '\n';
//    os << R"--(\hline)--" << '\n';
//    os << R"--(\textbf{Sample} & \textbf{Total Events} & \textbf{PS: Resolved} & \textbf{SB} & \textbf{SR}
//           \\ \hline )--" << '\n';
//
//    for (auto& file: input_filenames)
//    {
//        SchottDataFrame df(MakeRDF(file, 8));
//
//        auto trigger_selection = df.Filter(
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
//        [](const RVec<Muon>& muons, RVec<Electron> electrons)
//        {
//            return (electrons.size()==2 && muons.empty()); //keep events which have exactly 2 electrons for di_electrons and no muons
//
//        }, {"muons", "di_electrons"});
//
//        //new dataframe node: contains only the events from `two_leptons` whose electrons in the `di_electrons` branch have opposite charge
//        auto opp_charge = two_leptons.Filter([](const RVec<Electron>& electrons)
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
//        //new dataframe node: contains only the events from `dilep` that have di-electron invariant mass between 81 and 101 GeV
//        auto mass = deltaR.Filter([] (RVec<Electron>& electrons)
//        {
//            auto mass = (electrons[0].Vector() + electrons[1].Vector()).M()/1e3;
//            return ((mass >= 81) && (mass <= 101));
//        }, {"di_electrons"});
//
//        //new dataframe node: contains only the events from `mass` that have dilepton pT > 10 GeV
//        auto ptCut = mass.Filter([] (RVec<Electron>& electrons)
//        {
//            auto pT = (electrons[0].Vector() + electrons[1].Vector()).Pt()/1e3;
//            return pT > 10;
//        }, {"di_electrons"});
//
//        //end pre-selection -----------
//
//        //photon cuts
//        auto resolved = ptCut.Define("photons_pass_cuts",
//        [&](RVec<Photon> photons)
//        {
//            photons.erase(std::remove_if(photons.begin(),photons.end(),
//            [](Photon& x)
//            {
//                return ((std::abs(x.photon_eta) >= 2.37) || (std::abs(x.photon_eta) > 1.37 && std::abs(x.photon_eta) < 1.52) || (!x.photon_id_loose));
//
//            }), photons.end());
//            return photons;
//        }, {"photons"}).Filter(
//        [&](RVec<Photon>& reco_photons_matched)
//        {
//            if (reco_photons_matched.size() < 2)
//            {
//                return false;
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
//                //if it's the first combination or if new X is closer to 1
//                //than current best_X and ΔR
//                //between the two reco-photons < 1.5, then update best_X, pt1, pt2,
//                //and the corresponding reco-photons x
//                if (i==0 || ((std::abs(1-X) < std::abs(1-best_X)) and (delta_r < 1.5)))
//                {
//                    best_X = X;
//                    pt1 = reco_photons_matched[combs[0][i]].photon_pt;
//                    pt2 = reco_photons_matched[combs[1][i]].photon_pt;
//                    chosen_delta_r = delta_r;
//                }
//            }
//            if (chosen_delta_r < 1.5 && pt1 > 10e3 && pt2 > 10e3 && best_X > 0.96 && best_X < 1.2) //two photons corresponding to best_X must both have p_T > 10 GeV, ΔR < 1.5, and 0.96 < best_X < 1.2
//            {
//                return true;
//            }
//            return false;
//
//        }, {"photons_pass_cuts"})
//        .Define("mass", [&](RVec<Electron>& electrons, RVec<Photon>& reco_photons_matched)
//        {
//            return ((electrons[0].Vector()+electrons[1].Vector()).M() + (reco_photons_matched[0].Vector()+reco_photons_matched[1].Vector()).M())/1e3;
//        }, {"di_electrons", "photons_pass_cuts"});
//
//        auto SB = resolved.Filter(
//        [&](double mass)
//        {
//            return (!((mass > 110) && (mass < 140)));
//        }, {"mass"});
//
//        auto SR = resolved.Filter(
//        [&](double mass)
//        {
//            return ((mass > 110) && (mass < 140));
//        }, {"mass"});
//
//        Nodes.push_back(df.Count());
//        Nodes.push_back(resolved.Count());
//        Nodes.push_back(SB.Count());
//        Nodes.push_back(SR.Count());
//    }
//
//    ROOT::RDF::RunGraphs(Nodes); // running all computation nodes concurrently
//
//    std::cout << Nodes.size() << '\n';
//
////    0       1       2       3   //Z-gamma
////    4       5       6       7   //Z-gamma
////    8       9       10      11  //Z-gamma
////    12      13      14      15  //ma1
////    16      17      18      19  //ma2
////    20      21      22      23  //ma3
////    24      25      26      27  //ma5
////    28      29      30      31  //ma9
////    32      33      34      35  //data
////    36      37      38      39  //Z+jets
////    40      41      42      43  //Z+jets
////    44      45      46      47  //Z+jets
////    48      49      50      51  //Z+jets
////    52      53      54      55  //Z+jets
////    56      57      58      59  //Z+jets
////    60      61      62      63  //Z+jets
////    64      65      66      67  //Z+jets
////    68      69      70      71  //Z+jets
//
//    os.setf(std::ios::fixed);
//    os.precision(2);
//
//    for (int i=0, j=0; (i<18 && j <= 68); i++, j+=4)
//    {
//        os << Samples[i];
//        if (i >= 0 && i <= 2)
//        {
//            os  << " & " << *Nodes[j].GetResultPtr<ULong64_t>() * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j].GetResultPtr<ULong64_t>()) * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << " & " << *Nodes[j+1].GetResultPtr<ULong64_t>() * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+1].GetResultPtr<ULong64_t>()) * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << " & " << *Nodes[j+2].GetResultPtr<ULong64_t>() * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+2].GetResultPtr<ULong64_t>()) * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << " & " << *Nodes[j+3].GetResultPtr<ULong64_t>() * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+3].GetResultPtr<ULong64_t>()) * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << R"--( \\ \hline )--" << '\n';
//        }
//        else if (i >= 9)
//        {
//            os  << " & " << *Nodes[j].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << " & " << *Nodes[j+1].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+1].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << " & " << *Nodes[j+2].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+2].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << " & " << *Nodes[j+3].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+3].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << R"--( \\ \hline )--" << '\n';
//        }
//
//        else
//        {
//            os  << " & " << *Nodes[j].GetResultPtr<ULong64_t>()
//                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j].GetResultPtr<ULong64_t>())
//                << " & " << *Nodes[j+1].GetResultPtr<ULong64_t>()
//                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+1].GetResultPtr<ULong64_t>())
//                << " & " << *Nodes[j+2].GetResultPtr<ULong64_t>()
//                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+2].GetResultPtr<ULong64_t>())
//                << " & " << *Nodes[j+3].GetResultPtr<ULong64_t>()
//                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+3].GetResultPtr<ULong64_t>())
//                << R"--( \\ \hline )--" << '\n';
//        }
//    }
//
//    for (int i = 0, j = 0; (i <= 8 && j <= 2); i += 4, j++)
//    {
//        totalEventsZgamma += *Nodes[i].GetResultPtr<ULong64_t>() * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
//        resolvedEventsZgamma += *Nodes[i+1].GetResultPtr<ULong64_t>() * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
//        SBEventsZgamma += *Nodes[i+2].GetResultPtr<ULong64_t>() * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
//        SREventsZgamma += *Nodes[i+3].GetResultPtr<ULong64_t>() * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
//
//        totalEventsZgammaStatUnc += pow(sqrt(*Nodes[i].GetResultPtr<ULong64_t>()) * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
//        resolvedEventsZgammaStatUnc += pow(sqrt(*Nodes[i+1].GetResultPtr<ULong64_t>()) * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
//        SBEventsZgammaStatUnc += pow(sqrt(*Nodes[i+2].GetResultPtr<ULong64_t>()) * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
//        SREventsZgammaStatUnc += pow(sqrt(*Nodes[i+3].GetResultPtr<ULong64_t>()) * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
//    }
//
//    for (int i = 36, j = 0; (i <= 68 && j <= 8); i += 4, j++)
//    {
//        totalEventsZJets += *Nodes[i].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
//        resolvedEventsZJets += *Nodes[i+1].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
//        SBEventsZJets += *Nodes[i+2].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
//        SREventsZJets += *Nodes[i+3].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
//
//        totalEventsZJetsStatUnc += pow(sqrt(*Nodes[i].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
//        resolvedEventsZJetsStatUnc += pow(sqrt(*Nodes[i+1].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
//        SBEventsZJetsStatUnc += pow(sqrt(*Nodes[i+2].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
//        SREventsZJetsStatUnc += pow(sqrt(*Nodes[i+3].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
//    }
//
//    os << R"--(Total $Z\gamma$ & )--" << totalEventsZgamma
//    << R"--($\, \pm \,$)--" << sqrt(totalEventsZgammaStatUnc)
//    << " & " << resolvedEventsZgamma
//    << R"--($\, \pm \,$)--" << sqrt(resolvedEventsZgammaStatUnc)
//    << " & " << SBEventsZgamma
//    << R"--($\, \pm \,$)--" << sqrt(SBEventsZgammaStatUnc)
//    << " & " << SREventsZgamma
//    << R"--($\, \pm \,$)--" << sqrt(SREventsZgammaStatUnc)
//    << R"--( \\ \hline )--" << '\n';
//
//    os << R"--(Total $Z$+jets & )--" << totalEventsZJets
//    << R"--($\, \pm \,$)--" << sqrt(totalEventsZJetsStatUnc)
//    << " & " << resolvedEventsZJets
//    << R"--($\, \pm \,$)--" << sqrt(resolvedEventsZJetsStatUnc)
//    << " & " << SBEventsZJets
//    << R"--($\, \pm \,$)--" << sqrt(SBEventsZJetsStatUnc)
//    << " & " << SREventsZJets
//    << R"--($\, \pm \,$)--" << sqrt(SREventsZJetsStatUnc)
//    << R"--( \\ \hline )--" << '\n';
//
//    os << R"--(Total Bkg & )--" << totalEventsZJets+totalEventsZgamma
//    << R"--($\, \pm \,$)--" << sqrt(totalEventsZJetsStatUnc+totalEventsZgammaStatUnc)
//    << " & " << resolvedEventsZJets+resolvedEventsZgamma
//    << R"--($\, \pm \,$)--" << sqrt(resolvedEventsZJetsStatUnc+resolvedEventsZgammaStatUnc)
//    << " & " << SBEventsZJets+SBEventsZgamma
//    << R"--($\, \pm \,$)--" << sqrt(SBEventsZJetsStatUnc+SBEventsZgammaStatUnc)
//    << " & " << SREventsZJets+SREventsZgamma
//    << R"--($\, \pm \,$)--" << sqrt(SREventsZJetsStatUnc+SREventsZgammaStatUnc)
//    << R"--( \\ \hline )--" << '\n';
//
//    os << R"--(\end{tabular}})--" << '\n';
//    std::cout << os.str() << '\n';
//}
//
//void Table11()
//{
//    std::vector<std::vector<std::string>> input_filenames =
//    {
//        //Z gamma background
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617070._000001.LGNTuple.root"},
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617064._000001.LGNTuple.root"},
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617074._000001.LGNTuple.root"},
//        //Signal
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root"}, //1 GeV
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600751.PhPy8EG_AZNLO_ggH125_mA2p0_v1.root"}, // 2 GeV
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600752.PhPy8EG_AZNLO_ggH125_mA3p0_v1.root"}, // 3 GeV
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root"}, //5 GeV
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600756.PhPy8EG_AZNLO_ggH125_mA9p0_v1.root"}, // 9 GeV
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
//    constexpr std::array<const char*,18> Samples = {R"--(pty2\_9\_17)--", R"--(pty\_17\_myy\_0\_80)--", R"--(pty\_17\_myy\_80)--", R"--(Signal $m_{\text{A}}$ = 1 GeV)--", R"--(Signal $m_{\text{A}}$ = 5 GeV)--", R"--(Signal $m_{\text{A}}$ = 2 GeV)--", R"--(Signal $m_{\text{A}}$ = 3 GeV)--", R"--(Signal $m_{\text{A}}$ = 9 GeV)--", "Data", R"--(Zee\_lightJet\_0-70)--", R"--(Zee\_lightJet\_70-140)--", R"--(Zee\_lightJet\_140-280)--", R"--(Zee\_cJet\_0-70)--", R"--(Zee\_cJet\_70-140)--", R"--(Zee\_cJet\_140-280)--", R"--(Zee\_bJet\_0-70)--", R"--(Zee\_bJet\_70-140)--", R"--(Zee\_bJet\_140-280)--"};
//
//    std::array<double,9> JetNumeratorSFs = {((139e15)*(1.9828e-9)*(0.821204)),((139e15)*(110.64e-12)*(0.69275)),((139e15)*(40.645e-12)*(0.615906)),((139e15)*(1.9817e-9)*(0.1136684)),((139e15)*(110.47e-12)*(0.1912956)),((139e15)*(40.674e-12)*(0.2326772)),((139e15)*(1.9819e-9)*(0.0656969)),((139e15)*(110.53e-12)*(0.1158741)),((139e15)*(40.68e-12)*(0.1535215))}; //numerators for jet bkg
//    std::array<double,3> SFs = {((139e15)*(.871e-12)),((139e15)*(.199e-12)), ((139e15)*(.0345e-15))}; //numerators for Z-gamma bkg
//
//    double totalEventsZgamma = 0, passPreselectionZgamma = 0, photonPtDeltaRCountZgamma = 0, xWindowZgamma = 0,
//    srCountZgamma = 0, srIDCountZgamma = 0;
//    double totalEventsZgammaStatUnc = 0, passPreselectionZgammaStatUnc = 0, photonPtDeltaRCountZgammaStatUnc = 0, xWindowZgammaStatUnc = 0, srCountZgammaStatUnc = 0, srIDCountZgammaStatUnc = 0;
//
//    double totalEventsZjets = 0, passPreselectionZjets = 0, photonPtDeltaRCountZjets = 0, xWindowZjets = 0,
//        srCountZjets = 0, srIDCountZjets = 0;
//    double totalEventsZjetsStatUnc = 0, passPreselectionZjetsStatUnc = 0, photonPtDeltaRCountZjetsStatUnc = 0, xWindowZjetsStatUnc = 0, srCountZjetsStatUnc = 0, srIDCountZjetsStatUnc = 0;
//
//    std::vector<ROOT::RDF::RResultHandle> Nodes;
//
//    std::ostringstream os, ss;
//
//    os.setf(std::ios::fixed);
//    os.precision(2);
//
//    os << R"--(\section*{Table 11})--" << '\n';
//    os << R"--(\hspace{-3cm}\scalebox{0.55}{)--" << '\n';
//    os << R"--(\setlength\extrarowheight{2pt})--" << '\n';
//    os << R"--(\renewcommand{\arraystretch}{1.5})--" << '\n';
//    os << R"--(\begin{tabular}{|c|c|c|c|c|c|c|})--" << '\n';
//    os << R"--(\hline)--" << '\n';
//    os << R"--(\textbf{Sample} & \textbf{Total Events} & \textbf{pass preselection (PS)} & \textbf{photon} $\pmb{p_T}$ \textbf{+} $\pmb{\Delta R_{\gamma\gamma}}$ \textbf{cut} & $\pmb{X}$ \textbf{window} & \textbf{SR} & \textbf{SR-ID}
//           \\ \hline )--" << '\n';
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
//    for (auto& file: input_filenames)
//    {
//        SchottDataFrame df(MakeRDF(file, 8));
//
//        auto trigger_selection = df.Filter(
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
//        [](const RVec<Muon>& muons, RVec<Electron> electrons)
//        {
//            return (electrons.size()==2 && muons.empty()); //keep events which have exactly 2 electrons for di_electrons and no muons
//
//        }, {"muons", "di_electrons"});
//
//        //new dataframe node: contains only the events from `two_leptons` whose electrons in the `di_electrons` branch have opposite charge
//        auto opp_charge = two_leptons.Filter([](const RVec<Electron>& electrons)
//        {
//            return (electrons[0].electron_charge*electrons[1].electron_charge < 0);
//
//        }, {"di_electrons"});
//
//        //new dataframe node: contains only the events from `opp_charge` that have 1 electron with pt > 20 GeV and the other with pt > 27 GeV
//        auto leading_pt = opp_charge.Filter([](RVec<Electron>& electrons)
//        {
//            return ((electrons[0].electron_pt > 20e3 && electrons[1].electron_pt > 27e3) || (electrons[1].electron_pt > 20e3 && electrons[0].electron_pt > 27e3));
//        }, {"di_electrons"});
//
////        auto delta_R = leading_pt.Filter([] (RVec<Electron>& electrons)
////        {
////            return (DeltaR(electrons[0].Vector(), electrons[1].Vector()) > 0.01);
////        }, {"di_electrons"});
//
//        auto same_flavour = leading_pt.Filter([] (RVec<Electron>& electrons)
//        {
//            return true;
//        }, {"di_electrons"});
//
//        auto dilep = same_flavour.Define("dilep",[] (RVec<Electron>& electrons)
//        {
//            return (electrons[0].Vector() + electrons[1].Vector());
//        }, {"di_electrons"});
//
//        //new dataframe node: contains only the events from `dilep` that have di-electron invariant mass between 81 and 101 GeV
//        auto mass = dilep.Filter([] (RVec<Electron>& electrons)
//        {
//            auto mass = (electrons[0].Vector() + electrons[1].Vector()).M()/1e3;
//            return ((mass >= 81) && (mass <= 101));
//        }, {"di_electrons"});
//
//        //new dataframe node: contains only the events from `mass` that have dilepton pT > 10 GeV
//        auto ptCut = mass.Filter([] (RVec<Electron>& electrons)
//        {
//            auto pT = (electrons[0].Vector() + electrons[1].Vector()).Pt()/1e3;
//            return pT > 10;
//        }, {"di_electrons"});
//
//        //end pre-selection -----------
//
//        //photon cuts
//        auto photonPtDeltaR = ptCut.Define("photonPtDeltaR",
//        [&](RVec<Photon> photons)
//        {
//            photons.erase(std::remove_if(photons.begin(),photons.end(),
//            [](Photon& x)
//            {
//                return ((std::abs(x.photon_eta) >= 2.37) || (x.photon_pt <= 10e3) || (std::abs(x.photon_eta) > 1.37 && std::abs(x.photon_eta) < 1.52));
//
//            }), photons.end());
//            return photons;
//        }, {"photons"}).Filter(
//       [&](RVec<Photon>& reco_photons_matched)
//       {
//            if (reco_photons_matched.size() < 2)
//            {
//              return false;
//            }
//            RVec<Photon> x;
//            auto combs = Combinations(reco_photons_matched, 2);
//            size_t length = combs[0].size();
//            double delta_r, m, pt, X, best_X, pt1, pt2, chosen_delta_r;
//            for (size_t i=0; i<length; i++)
//            {
//                delta_r = DeltaR(reco_photons_matched[combs[0][i]].Vector(), reco_photons_matched[combs[1][i]].Vector());
//                m = (reco_photons_matched[combs[0][i]].Vector() + reco_photons_matched[combs[1][i]].Vector()).M();
//                pt = (reco_photons_matched[combs[0][i]].Vector() + reco_photons_matched[combs[1][i]].Vector()).Pt();
//                X = delta_r*(pt/(2.0*m));
//                //if it's the first combination or if new X is closer to 1
//                //than current best_X and ΔR
//                //between the two reco-photons < 1.5, then update best_X, pt1, pt2,
//                //and the corresponding reco-photons x
//                if (i==0 || ((std::abs(1-X) < std::abs(1-best_X)) and (delta_r < 1.5)))
//                {
//                    best_X = X;
//                    pt1 = reco_photons_matched[combs[0][i]].photon_pt;
//                    pt2 = reco_photons_matched[combs[1][i]].photon_pt;
//                    chosen_delta_r = delta_r;
//                    x = {reco_photons_matched[combs[0][i]], reco_photons_matched[combs[1][i]]};
//                }
//            }
//            return (chosen_delta_r < 1.5 && pt1 > 10e3 && pt2 > 10e3);
//
//       }, {"photonPtDeltaR"});
//
//        auto X_window = photonPtDeltaR.Define("chosen_two",
//        [](RVec<Photon>& reco_photons_matched)
//        {
//            RVec<Photon> x;
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
//                if (i==0 || ((std::abs(1-X) < std::abs(1-best_X)) and (delta_r < 1.5)))
//                {
//                    best_X = X;
//                    pt1 = reco_photons_matched[combs[0][i]].photon_pt;
//                    pt2 = reco_photons_matched[combs[1][i]].photon_pt;
//                    chosen_delta_r = delta_r;
//                    x = {reco_photons_matched[combs[0][i]], reco_photons_matched[combs[1][i]]};
//                }
//            }
//            if (pt1 > 10e3 && pt2 > 10e3 && best_X > 0.96 && best_X < 1.2 && chosen_delta_r < 1.5)
//            {
//                return x;
//            }
//            x.clear();
//            return x;
//        }, {"photonPtDeltaR"}).Filter(
//        [&](RVec<Photon>& reco_photons_matched)
//        {
//            return (reco_photons_matched.size()==2);
//        }, {"chosen_two"});
//
//        auto SR = X_window.Filter(
//        [](RVec<Photon>& photons, RVec<Electron>& electrons)
//        {
//            PtEtaPhiEVector photonVec = photons[0].Vector() + photons[1].Vector();
//            PtEtaPhiEVector electronVec = electrons[0].Vector() + electrons[1].Vector();
//            auto mass = (photonVec+electronVec).M()/1e3;
//            return ((mass >= 110) && (mass <= 140));
//        },{"chosen_two","di_electrons"});
//
//        auto SR_ID = SR.Filter(
//        [](RVec<Photon>& photons)
//        {
//            return (photons[0].photon_id_loose && photons[1].photon_id_loose);
//        },{"chosen_two"});
//
//        Nodes.push_back(df.Count());
//        Nodes.push_back(ptCut.Count());
//        Nodes.push_back(photonPtDeltaR.Count());
//        Nodes.push_back(X_window.Count());
//        Nodes.push_back(SR.Count());
//        Nodes.push_back(SR_ID.Count());
//    }
//
//    ROOT::RDF::RunGraphs(Nodes); // running all computation nodes concurrently
//
////    0       1       2       3       4       5       //Z-gamma
////    6       7       8       9       10      11      //Z-gamma
////    12      13      14      15      16      17      //Z-gamma
////    18      19      20      21      22      23      //ma1
////    24      25      26      27      28      29      //ma2
////    30      31      32      33      34      35      //ma3
////    36      37      38      39      40      41      //ma5
////    42      43      44      45      46      47      //ma9
////    48      49      50      51      52      53      //data
////    54      55      56      57      58      59      //Z+jets
////    60      61      62      63      64      65      //Z+jets
////    66      67      68      69      70      71      //Z+jets
////    72      73      74      75      76      77      //Z+jets
////    78      79      80      81      82      83      //Z+jets
////    84      85      86      87      88      89      //Z+jets
////    90      91      92      93      94      95      //Z+jets
////    96      97      98      99      100     101     //Z+jets
////    102     103     104     105     106     107     //Z+jets
//
//    std::cout << R"--(\section*{Table 11 Signal Ratios})--" << '\n';
//    std::cout << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
//    std::cout << R"--(\begin{tabular}{|c|c|c|c|c|c|c|})--" << '\n';
//    std::cout << R"--(\hline)--" << '\n';
//    std::cout << R"--(Sample & $\frac{\text{pass preselection (PS)}}{\text{pass preselection (PS)}}$ & $\frac{\text{photon } p_T \text{ + } \Delta R_{\gamma\gamma} \text{ cut }}{\text{pass preselection (PS)}}$ & $\frac{X \text{ window}}{\text{pass preselection (PS)}}$ & $\frac{\text{SR}}{\text{pass preselection (PS)}}$ & $\frac{\text{SR-ID}}{\text{pass preselection (PS)}}$
//           \\ \hline )--" << '\n';
//
//    ss << R"--(\section*{Table 11 Signal Ratios})--" << '\n';
//    ss << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
//    ss << R"--(\begin{tabular}{|c|c|c|c|c|c|c|})--" << '\n';
//    ss << R"--(\hline)--" << '\n';
//    ss << R"--(Sample & $\frac{\text{pass preselection (PS)}}{\text{pass preselection (PS)}}$ & $\frac{\text{photon } p_T \text{ + } \Delta R_{\gamma\gamma} \text{ cut }}{\text{pass preselection (PS)}}$ & $\frac{X \text{ window}}{\text{photon } p_T \text{ + } \Delta R_{\gamma\gamma} \text{ cut }}$ & $\frac{\text{SR}}{X \text{ window}}$ & $\frac{\text{SR-ID}}{\text{SR}}$
//           \\ \hline )--" << '\n';
//
//    for (int i=0, j=0; (i<18 && j <= 102); i++, j+=6)
//    {
//        os << Samples[i];
//
//        if (i >= 0 && i <= 2)
//        {
//            os  << " & " << *Nodes[j].GetResultPtr<ULong64_t>() * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << R"--($\, \pm \,$)--" <<  sqrt(*Nodes[j].GetResultPtr<ULong64_t>()) * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << " & " << *Nodes[j+1].GetResultPtr<ULong64_t>() * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << R"--($\, \pm \,$)--" <<  sqrt(*Nodes[j+1].GetResultPtr<ULong64_t>()) * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << " & " << *Nodes[j+2].GetResultPtr<ULong64_t>() * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << R"--($\, \pm \,$)--" <<  sqrt(*Nodes[j+2].GetResultPtr<ULong64_t>()) * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << " & " << *Nodes[j+3].GetResultPtr<ULong64_t>() * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << R"--($\, \pm \,$)--" <<  sqrt(*Nodes[j+3].GetResultPtr<ULong64_t>()) * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << " & " << *Nodes[j+4].GetResultPtr<ULong64_t>() * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << R"--($\, \pm \,$)--" <<  sqrt(*Nodes[j+4].GetResultPtr<ULong64_t>()) * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << " & " << *Nodes[j+5].GetResultPtr<ULong64_t>() * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << R"--($\, \pm \,$)--" <<  sqrt(*Nodes[j+5].GetResultPtr<ULong64_t>()) * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << R"--( \\ \hline )--" << '\n';
//        }
//        else if (i >= 9)
//        {
//            os  << " & " << *Nodes[j].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << " & " << *Nodes[j+1].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+1].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << " & " << *Nodes[j+2].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+2].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << " & " << *Nodes[j+3].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+3].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << " & " << *Nodes[j+4].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+4].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << " & " << *Nodes[j+5].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+5].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
//                << R"--( \\ \hline )--" << '\n';
//        }
//        else
//        {
//            os  << " & " << *Nodes[j].GetResultPtr<ULong64_t>()
//                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j].GetResultPtr<ULong64_t>())
//                << " & " << *Nodes[j+1].GetResultPtr<ULong64_t>()
//                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+1].GetResultPtr<ULong64_t>())
//                << " & " << *Nodes[j+2].GetResultPtr<ULong64_t>()
//                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+2].GetResultPtr<ULong64_t>())
//                << " & " << *Nodes[j+3].GetResultPtr<ULong64_t>()
//                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+3].GetResultPtr<ULong64_t>())
//                << " & " << *Nodes[j+4].GetResultPtr<ULong64_t>()
//                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+4].GetResultPtr<ULong64_t>())
//                << " & " << *Nodes[j+5].GetResultPtr<ULong64_t>()
//                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+5].GetResultPtr<ULong64_t>())
//                << R"--( \\ \hline )--" << '\n';
//
//            if (i==3) //1 GeV
//            {
//                std::cout << Samples[i] << " (Me) & " << 1 << " & " <<
//                static_cast<double>(*Nodes[j+2].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
//                << " & " <<
//                static_cast<double>(*Nodes[j+3].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
//                << " & " <<
//                static_cast<double>(*Nodes[j+4].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
//                << " & " <<
//                static_cast<double>(*Nodes[j+5].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
//                << R"--( \\ \hline )--" << '\n';
//
//                std::cout << Samples[i] << " (Paper) & " << 1 << " & " <<
//                2343.65 / 19516.87
//                << " & " <<
//                2204.22 / 19516.87
//                << " & " <<
//                2076.73 / 19516.87
//                << " & " <<
//                333.21 / 19516.87
//                << R"--( \\ \hline )--" << '\n';
//
//                ss << Samples[i] << " (Me) & " << 1 << " & " <<
//                static_cast<double>(*Nodes[j+2].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
//                << " & " <<
//                static_cast<double>(*Nodes[j+3].GetResultPtr<ULong64_t>()) / *Nodes[j+2].GetResultPtr<ULong64_t>()
//                << " & " <<
//                static_cast<double>(*Nodes[j+4].GetResultPtr<ULong64_t>()) / *Nodes[j+3].GetResultPtr<ULong64_t>()
//                << " & " <<
//                static_cast<double>(*Nodes[j+5].GetResultPtr<ULong64_t>()) / *Nodes[j+4].GetResultPtr<ULong64_t>()
//                << R"--( \\ \hline )--" << '\n';
//
//                ss << Samples[i] << " (Paper) & " << 1 << " & " <<
//                2343.65 / 19516.87
//                << " & " <<
//                2204.22 / 2343.65
//                << " & " <<
//                2076.73 / 2204.22
//                << " & " <<
//                333.21 / 2076.73
//                << R"--( \\ \hline )--" << '\n';
//
//            }
//
//            if (i==6) //5 GeV
//            {
//                std::cout << Samples[i] << " (Me) & " << 1 << " & " <<
//                static_cast<double>(*Nodes[j+2].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
//                << " & " <<
//                static_cast<double>(*Nodes[j+3].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
//                << " & " <<
//                static_cast<double>(*Nodes[j+4].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
//                << " & " <<
//                static_cast<double>(*Nodes[j+5].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
//                << R"--( \\ \hline )--" << '\n';
//
//                std::cout << Samples[i] << " (Paper) & " << 1 << " & " <<
//                3245.86 / 19536.68
//                << " & " <<
//                3052.81 / 19536.68
//                << " & " <<
//                2941.28 / 19536.68
//                << " & " <<
//                1959.68 / 19536.68
//                << R"--( \\ \hline )--" << '\n';
//
//                ss << Samples[i] << " (Me) & " << 1 << " & " <<
//                static_cast<double>(*Nodes[j+2].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
//                << " & " <<
//                static_cast<double>(*Nodes[j+3].GetResultPtr<ULong64_t>()) / *Nodes[j+2].GetResultPtr<ULong64_t>()
//                << " & " <<
//                static_cast<double>(*Nodes[j+4].GetResultPtr<ULong64_t>()) / *Nodes[j+3].GetResultPtr<ULong64_t>()
//                << " & " <<
//                static_cast<double>(*Nodes[j+5].GetResultPtr<ULong64_t>()) / *Nodes[j+4].GetResultPtr<ULong64_t>()
//                << R"--( \\ \hline )--" << '\n';
//
//                ss << Samples[i] << " (Paper) & " << 1 << " & " <<
//                3245.86 / 19536.68
//                << " & " <<
//                3052.81 / 3245.86
//                << " & " <<
//                2941.28 / 3052.81
//                << " & " <<
//                1959.68 / 2941.28
//                << R"--( \\ \hline )--" << '\n';
//            }
//        }
//    }
//
//    std::cout << R"--(\end{tabular}})--" << "\n\n\n";
//
//    ss << R"--(\end{tabular}})--" << "\n\n\n";
//    std::cout << ss.str() << '\n';
//
//    for (int i = 0, j = 0; (i <= 12 && j <= 2); i += 6, j++)
//    {
//        totalEventsZgamma += *Nodes[i].GetResultPtr<ULong64_t>() * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
//        passPreselectionZgamma += *Nodes[i+1].GetResultPtr<ULong64_t>() * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
//        photonPtDeltaRCountZgamma += *Nodes[i+2].GetResultPtr<ULong64_t>() * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
//        xWindowZgamma += *Nodes[i+3].GetResultPtr<ULong64_t>() * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
//        srCountZgamma += *Nodes[i+4].GetResultPtr<ULong64_t>() * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
//        srIDCountZgamma += *Nodes[i+5].GetResultPtr<ULong64_t>() * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
//
//        totalEventsZgammaStatUnc += pow(sqrt(*Nodes[i].GetResultPtr<ULong64_t>()) * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
//        passPreselectionZgammaStatUnc += pow(sqrt(*Nodes[i+1].GetResultPtr<ULong64_t>()) * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
//        photonPtDeltaRCountZgammaStatUnc += pow(sqrt(*Nodes[i+2].GetResultPtr<ULong64_t>()) * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
//        xWindowZgammaStatUnc += pow(sqrt(*Nodes[i+3].GetResultPtr<ULong64_t>()) * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
//        srCountZgammaStatUnc += pow(sqrt(*Nodes[i+4].GetResultPtr<ULong64_t>()) * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
//        srIDCountZgammaStatUnc += pow(sqrt(*Nodes[i+5].GetResultPtr<ULong64_t>()) * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
//    }
//
//    for (int i = 54, j = 0; (i <= 102 && j <= 8); i += 6, j++)
//    {
//        totalEventsZjets += *Nodes[i].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
//        passPreselectionZjets += *Nodes[i+1].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
//        photonPtDeltaRCountZjets += *Nodes[i+2].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
//        xWindowZjets += *Nodes[i+3].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
//        srCountZjets += *Nodes[i+4].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
//        srIDCountZjets += *Nodes[i+5].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
//
//        totalEventsZjetsStatUnc += pow(sqrt(*Nodes[i].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
//        passPreselectionZjetsStatUnc += pow(sqrt(*Nodes[i+1].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
//        photonPtDeltaRCountZjetsStatUnc += pow(sqrt(*Nodes[i+2].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
//        xWindowZjetsStatUnc += pow(sqrt(*Nodes[i+3].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
//        srCountZjetsStatUnc += pow(sqrt(*Nodes[i+4].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
//        srIDCountZjetsStatUnc += pow(sqrt(*Nodes[i+5].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
//    }
//
//    os << R"--(Total $Z\gamma$ & )--" << totalEventsZgamma
//    << R"--($\, \pm \,$)--" << sqrt(totalEventsZgammaStatUnc)
//    << " & " << passPreselectionZgamma
//    << R"--($\, \pm \,$)--" << sqrt(passPreselectionZgammaStatUnc)
//    << " & " << photonPtDeltaRCountZgamma
//    << R"--($\, \pm \,$)--" << sqrt(photonPtDeltaRCountZgammaStatUnc)
//    << " & " << xWindowZgamma
//    << R"--($\, \pm \,$)--" << sqrt(xWindowZgammaStatUnc)
//    << " & " << srCountZgamma
//    << R"--($\, \pm \,$)--" << sqrt(srCountZgammaStatUnc)
//    << " & " << srIDCountZgamma
//    << R"--($\, \pm \,$)--" << sqrt(srIDCountZgammaStatUnc)
//    << R"--( \\ \hline )--" << '\n';
//
//    os << R"--(Total $Z$+jets & )--" << totalEventsZjets
//    << R"--($\, \pm \,$)--" << sqrt(totalEventsZjetsStatUnc)
//    << " & " << passPreselectionZjets
//    << R"--($\, \pm \,$)--" << sqrt(passPreselectionZjetsStatUnc)
//    << " & " << photonPtDeltaRCountZjets
//    << R"--($\, \pm \,$)--" << sqrt(photonPtDeltaRCountZjetsStatUnc)
//    << " & " << xWindowZjets
//    << R"--($\, \pm \,$)--" << sqrt(xWindowZjetsStatUnc)
//    << " & " << srCountZjets
//    << R"--($\, \pm \,$)--" << sqrt(srCountZjetsStatUnc)
//    << " & " << srIDCountZjets
//    << R"--($\, \pm \,$)--" << sqrt(srIDCountZjetsStatUnc)
//    << R"--( \\ \hline )--" << '\n';
//
//    os << R"--(Total Bkg & )--" << totalEventsZjets+totalEventsZgamma
//    << R"--($\, \pm \,$)--" << sqrt(totalEventsZjetsStatUnc+totalEventsZgammaStatUnc)
//    << " & " << passPreselectionZjets+passPreselectionZgamma
//    << R"--($\, \pm \,$)--" << sqrt(passPreselectionZjetsStatUnc+passPreselectionZgammaStatUnc)
//    << " & " << photonPtDeltaRCountZjets+photonPtDeltaRCountZgamma
//    << R"--($\, \pm \,$)--" << sqrt(photonPtDeltaRCountZjetsStatUnc+photonPtDeltaRCountZgammaStatUnc)
//    << " & " << xWindowZjets+xWindowZgamma
//    << R"--($\, \pm \,$)--" << sqrt(xWindowZjetsStatUnc+xWindowZgammaStatUnc)
//    << " & " << srCountZjets+srCountZgamma
//    << R"--($\, \pm \,$)--" << sqrt(srCountZjetsStatUnc+srCountZgammaStatUnc)
//    << " & " << srIDCountZjets+srIDCountZgamma
//    << R"--($\, \pm \,$)--" << sqrt(srIDCountZjetsStatUnc+srIDCountZgammaStatUnc)
//        << R"--( \\ \hline )--" << '\n';
//
//    os << R"--(\end{tabular}})--" << '\n';
//    std::cout << os.str() << '\n';
//}
//
//void Table3_Displaced_Axions()
//{
//    std::vector<std::vector<std::string>> input_filenames =
//    {
//        //Displaced Signal
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
//    };
//
//    constexpr std::array<const char*,15> Samples = {R"--(Prompt Signal $m_{\text{A}}$ = 1 GeV)--", R"--(Prompt Signal $m_{\text{A}}$ = 5 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 0.2 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 0.5 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 1 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 3 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 5 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 10 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 10.1 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 20 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 20.1 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 29.5 GeV)--", };
//
//    std::vector<double> massPoints = {0.2,0.5,1,3,5,10,10.1,20,20.1,29.5};
//
//    std::ostringstream os;
//    os << R"--(\section*{Table 3 Prompt and Displaced Signal Samples})--" << '\n';
//    os << R"--(\hspace{-3cm}\scalebox{0.45}{)--" << '\n';
//    os << R"--(\begin{tabular}{|c|c|c|c|c|c|c|c|})--" << '\n';
//    os << R"--(\hline)--" << '\n';
//    os << R"--(Sample & Before Preselection & 2 leptons
//          & Opposite Charge & $p_{T}^{\text{leading}} > 27$ GeV, \; $p_{T}^{\text{sub-leading}} > 20$ GeV & Same flavour & dilep mass cut & dilep $p_{T}$ cut \\ \hline )--" << '\n';
//
//    std::vector<ROOT::RDF::RResultHandle> Nodes;
//
//    int counter = 0;
//
//    for (auto& file: input_filenames)
//    {
//        SchottDataFrame df(MakeRDF(file, 8));
//
//        auto trigger_selection = df.Define("truth_axions", [&](RVec<TruthParticle> truth_particles)
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
//              return false;
//            }
//            return true;
//
//        }, {"trigger_passed_triggers"});
//
//        auto two_leptons = trigger_selection.Define("di_electrons",
//        [](RVec<Electron> electrons)
//        {
//            electrons.erase(std::remove_if(electrons.begin(),electrons.end(),
//            [](Electron& ep)
//            {
//                return (!((ep.electron_pt/1e3 > 20) && (std::abs(ep.electron_eta) < 2.37) &&
//                                          (!((1.37 < std::abs(ep.electron_eta)) && (std::abs(ep.electron_eta) < 1.52)))
//                                          && (ep.electron_id_medium == 1)));
//
//            }), electrons.end());
//
//            return electrons;
//
//        },{"electrons"}).Filter([](RVec<Muon>& muons, RVec<Electron> di_electrons)
//        {
//            return (di_electrons.size()==2 && muons.empty() && DeltaR(di_electrons[0].Vector(), di_electrons[1].Vector()) > 0.01);
//
//        }, {"muons", "di_electrons"});
//
//        auto opp_charge = two_leptons.Filter([](RVec<Electron> di_electrons)
//        {
//            return (di_electrons[0].electron_charge*di_electrons[1].electron_charge < 0);
//
//        }, {"di_electrons"});
//
//        auto leading_pt = opp_charge.Filter([](RVec<Electron>& electrons)
//        {
//            return ((electrons[0].electron_pt >= 20e3 && electrons[1].electron_pt >= 27e3) || (electrons[1].electron_pt >= 20e3 && electrons[0].electron_pt >= 27e3));
//        }, {"di_electrons"});
//
//        auto same_flavour = leading_pt.Filter([] (RVec<Electron>& electrons)
//        {
//            return true;
//        }, {"di_electrons"});
//
//        auto dilep = same_flavour.Define("dilep",[] (RVec<Electron>& electrons)
//        {
//            return (electrons[0].Vector() + electrons[1].Vector());
//        }, {"di_electrons"});
//
//        auto mass = dilep.Filter([] (PtEtaPhiEVector& dilep)
//        {
//            auto mass = dilep.M()/1e3;
//            return ((mass >= 81) && (mass <= 101));
//        }, {"dilep"});
//
//        auto pt_cut = mass.Filter([] (PtEtaPhiEVector& dilep)
//        {
//            auto pT = dilep.Pt()/1e3;
//            return pT > 10;
//        }, {"dilep"});
//
//        if (counter < 2)
//        {
//            Nodes.push_back(trigger_selection.Count());
//            Nodes.push_back(two_leptons.Count());
//            Nodes.push_back(opp_charge.Count());
//            Nodes.push_back(leading_pt.Count());
//            Nodes.push_back(same_flavour.Count());
//            Nodes.push_back(mass.Count());
//            Nodes.push_back(pt_cut.Count());
//        }
//
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
//                auto mass_point_two_leptons = two_leptons.Filter([&]
//                (float axion_mass)
//                {
//                    return (roundToOneDecimalPlace(axion_mass) == mass_point);
//
//                }, {"axion_masses"});
//
//                auto mass_point_opp_charge = opp_charge.Filter([&]
//                (float axion_mass)
//                {
//                    return (roundToOneDecimalPlace(axion_mass) == mass_point);
//
//                }, {"axion_masses"});
//
//                auto mass_point_leading_pt = leading_pt.Filter([&]
//                (float axion_mass)
//                {
//                    return (roundToOneDecimalPlace(axion_mass) == mass_point);
//
//                }, {"axion_masses"});
//
//                auto mass_point_same_flavour = same_flavour.Filter([&]
//                (float axion_mass)
//                {
//                    return (roundToOneDecimalPlace(axion_mass) == mass_point);
//
//                }, {"axion_masses"});
//
//                auto mass_point_mass = mass.Filter([&]
//                (float axion_mass)
//                {
//                    return (roundToOneDecimalPlace(axion_mass) == mass_point);
//
//                }, {"axion_masses"});
//
//                auto mass_point_pt_cut = pt_cut.Filter([&]
//                (float axion_mass)
//                {
//                    return (roundToOneDecimalPlace(axion_mass) == mass_point);
//
//                }, {"axion_masses"});
//
//                Nodes.push_back(mass_point_trigger_selection.Count());
//                Nodes.push_back(mass_point_two_leptons.Count());
//                Nodes.push_back(mass_point_opp_charge.Count());
//                Nodes.push_back(mass_point_leading_pt.Count());
//                Nodes.push_back(mass_point_same_flavour.Count());
//                Nodes.push_back(mass_point_mass.Count());
//                Nodes.push_back(mass_point_pt_cut.Count());
//
//            }
//        }
//
//        counter++;
//    }
//
//    constexpr std::array<const char*,7> Cuts = {"total", "two leptons", "opposite charge", "leading pt", "same flavour", "mass", "pt cut"};
//
//    std::cout << Nodes.size();
//
//    ROOT::RDF::RunGraphs(Nodes); // running all computation nodes concurrently
//
////   0    1   2   3   4   5   6                          ma1
////   7    8   9   10  11  12  13                         ma5
////   14   15  16  17  18  19  20                         displaced_axion_1
////   21   22  23  24  25  26  27                         displaced_axion_2
////   28   29  30  31  32  33  34                         displaced_axion_3
////   35   36  37  38  39  40  41                         displaced_axion_4
////   42   43  44  45  46  47  48                         displaced_axion_5
////   49   50  51  52  53  54  55                         displaced_axion_6
////   56   57  58  59  60  61  62                         displaced_axion_7
////   63   64  65  66  67  68  69                         displaced_axion_8
////   70   71  72  73  74  75  76                         displaced_axion_9
////   77   78  79  80  81  82  83                         displaced_axion_10
//
//    std::cout << R"--(\section*{Table 3 Signal Ratios})--" << '\n';
//    std::cout << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
//    std::cout << R"--(\begin{tabular}{|c|c|c|c|c|c|c|})--" << '\n';
//    std::cout << R"--(\hline)--" << '\n';
//    std::cout << R"--(Sample & $\frac{\text{2 leptons}}{\text{2 leptons}}$
//          & $\frac{\text{Opposite Charge}}{\text{2 leptons}}$ & $\frac{p_{T}^{\text{leading}} > 27\text{ GeV, } \, p_{T}^{\text{sub-leading}} > 20 \text{ GeV }}{\text{2 leptons}}$ & $\frac{\text{Same flavour}}{\text{2 leptons}}$ & $\frac{\text{dilep mass cut}}{\text{2 leptons}}$ & $\frac{\text{dilep }p_{T} \text{ cut}}{\text{2 leptons}}$ \\ \hline )--" << '\n';
//
//    os.setf(std::ios::fixed);
//    os.precision(2);
//
//    for (int i=0, j=0; (i<12 && j <= 77); i++, j+=7)
//    {
//        os << Samples[i];
//        os  << " & " << *Nodes[j].GetResultPtr<ULong64_t>()
//        << R"--($\, \pm \,$)--" << sqrt(*Nodes[j].GetResultPtr<ULong64_t>())
//        << " & " << *Nodes[j+1].GetResultPtr<ULong64_t>()
//        << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+1].GetResultPtr<ULong64_t>())
//        << " & " << *Nodes[j+2].GetResultPtr<ULong64_t>()
//        << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+2].GetResultPtr<ULong64_t>())
//        << " & " << *Nodes[j+3].GetResultPtr<ULong64_t>()
//        << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+3].GetResultPtr<ULong64_t>())
//        << " & " << *Nodes[j+4].GetResultPtr<ULong64_t>()
//        << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+4].GetResultPtr<ULong64_t>())
//        << " & " << *Nodes[j+5].GetResultPtr<ULong64_t>()
//        << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+5].GetResultPtr<ULong64_t>())
//        << " & " << *Nodes[j+6].GetResultPtr<ULong64_t>()
//        << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+6].GetResultPtr<ULong64_t>())
//        << R"--( \\ \hline )--" << '\n';
//
//        if (i==0) //1 GeV
//        {
//            std::cout << Samples[i] << " (Me) & " << 1 << " & " <<
//            static_cast<double>(*Nodes[j+2].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
//            << " & " <<
//            static_cast<double>(*Nodes[j+3].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
//            << " & " <<
//            static_cast<double>(*Nodes[j+4].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
//            << " & " <<
//            static_cast<double>(*Nodes[j+5].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
//            << " & " <<
//            static_cast<double>(*Nodes[j+6].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
//            << R"--( \\ \hline )--" << '\n';
//
//            std::cout << Samples[i] << " (Paper) & " << 1 << " & " <<
//            21505.03 / 21606.75
//            << " & " <<
//            21375.48 / 21606.75
//            << " & " <<
//            21375.33 / 21606.75
//            << " & " <<
//            20543.06 / 21606.75
//            << " & " <<
//            19516.87 / 21606.75
//            << R"--( \\ \hline )--" << '\n';
//        }
//
//        else if (i==1) //5 GeV
//        {
//            std::cout << Samples[i] << " (Me) & " << 1 << " & " <<
//            static_cast<double>(*Nodes[j+2].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
//            << " & " <<
//            static_cast<double>(*Nodes[j+3].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
//            << " & " <<
//            static_cast<double>(*Nodes[j+4].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
//            << " & " <<
//            static_cast<double>(*Nodes[j+5].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
//            << " & " <<
//            static_cast<double>(*Nodes[j+6].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
//            << R"--( \\ \hline )--" << '\n';
//
//            std::cout << Samples[i] << " (Paper) & " << 1 << " & " <<
//            21556.03 / 21655.38
//            << " & " <<
//            21424.29 / 21655.38
//            << " & " <<
//            21424.28 / 21655.38
//            << " & " <<
//            20585.09 / 21655.38
//            << " & " <<
//            19536.68 / 21655.38
//            << R"--( \\ \hline )--" << '\n';
//        }
//
//        else
//        {
//            std::cout << Samples[i] << " Displaced (Me) & " << 1 << " & " <<
//            static_cast<double>(*Nodes[j+2].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
//            << " & " <<
//            static_cast<double>(*Nodes[j+3].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
//            << " & " <<
//            static_cast<double>(*Nodes[j+4].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
//            << " & " <<
//            static_cast<double>(*Nodes[j+5].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
//            << " & " <<
//            static_cast<double>(*Nodes[j+6].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
//            << R"--( \\ \hline )--" << '\n';
//        }
//    }
//
//    std::cout << R"--(\end{tabular}})--" << "\n\n\n";
//
//    os << R"--(\end{tabular}})--" << '\n';
//    std::cout << os.str() << '\n';
//}
//
//void Table8_Displaced_Axions()
//{
//    std::vector<std::vector<std::string>> input_filenames = {
//        //Signal
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root"},
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root"},
//        //Displaced Signal
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_FewMassPoints.root"},
//    };
//
//    constexpr std::array<const char*,15> Samples = {R"--(Prompt Signal $m_{\text{A}}$ = 1 GeV)--", R"--(Prompt Signal $m_{\text{A}}$ = 5 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 0.2 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 0.5 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 1 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 3 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 5 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 10 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 10.1 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 20 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 20.1 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 29.5 GeV)--", };
//
//    std::vector<double> massPoints = {0.2,0.5,1,3,5,10,10.1,20,20.1,29.5};
//
//    std::vector<ROOT::RDF::RResultHandle> Nodes;
//
//    std::ostringstream os;
//    os << R"--(\section*{Table 8 Prompt and Displaced Signal Samples})--" << '\n';
//    os << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
//    os << R"--(\begin{tabular}{|c|c|c|c|c|})--" << '\n';
//    os << R"--(\hline)--" << '\n';
//    os << R"--(Sample & Total Events & PS: Resolved & SB & SR
//           \\ \hline )--" << '\n';
//
//    int counter = 0;
//
//    for (auto& file: input_filenames)
//    {
//        SchottDataFrame df(MakeRDF(file, 8));
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
//        auto resolved = ptCut.Define("photons_pass_cuts",
//        [&](RVec<Photon> photons)
//        {
//            photons.erase(std::remove_if(photons.begin(),photons.end(),
//            [](Photon& x)
//            {
//                return ((std::abs(x.photon_eta) >= 2.37) || (std::abs(x.photon_eta) > 1.37 && std::abs(x.photon_eta) < 1.52) || (!x.photon_id_loose));
//
//            }), photons.end());
//            return photons;
//        }, {"photons"}).Filter(
//        [&](RVec<Photon>& reco_photons_matched)
//        {
//            if (reco_photons_matched.size() < 2)
//            {
//                return false;
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
//                    pt1 = reco_photons_matched[combs[0][i]].photon_pt;
//                    pt2 = reco_photons_matched[combs[1][i]].photon_pt;
//                    chosen_delta_r = delta_r;
//                }
//            }
//            if (pt1 > 10e3 && pt2 > 10e3 && best_X > 0.96 && best_X < 1.2 && chosen_delta_r < 1.5)
//            {
//                return true;
//            }
//            return false;
//
//        }, {"photons_pass_cuts"});
//
//        auto SR = resolved.Filter(
//        [&](RVec<Electron>& electrons, RVec<Photon>& reco_photons_matched)
//        {
//            auto mass = ((electrons[0].Vector()+electrons[1].Vector()).M() + (reco_photons_matched[0].Vector()+reco_photons_matched[1].Vector()).M())/1e3;
//            return ((mass > 110) && (mass < 140));
//        }, {"di_electrons", "photons_pass_cuts"});
//
//        auto SB = resolved.Filter(
//        [&](RVec<Electron>& electrons, RVec<Photon>& reco_photons_matched)
//        {
//            auto mass = ((electrons[0].Vector()+electrons[1].Vector()).M() + (reco_photons_matched[0].Vector()+reco_photons_matched[1].Vector()).M())/1e3;
//            return (!((mass > 110) && (mass < 140)));
//        }, {"di_electrons", "photons_pass_cuts"});
//
//        if (counter < 2)
//        {
//            Nodes.push_back(trigger_selection.Count());
//            Nodes.push_back(resolved.Count());
//            Nodes.push_back(SB.Count());
//            Nodes.push_back(SR.Count());
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
//                auto mass_point_resolved = resolved.Filter([&]
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
//                Nodes.push_back(mass_point_trigger_selection.Count());
//                Nodes.push_back(mass_point_resolved.Count());
//                Nodes.push_back(mass_point_SB.Count());
//                Nodes.push_back(mass_point_SR.Count());
//            }
//        }
//
//        counter++;
//    }
//
//    ROOT::RDF::RunGraphs(Nodes); // running all computation nodes concurrently
//
////       0    1   2   3                           ma1
////       4    5   6   7                           ma5
////       8    9   10  11                          displaced_axion_1
////       12   13  14  15                          displaced_axion_2
////       16   17  18  19                          displaced_axion_3
////       20   21  22  23                          displaced_axion_4
////       24   25  26  27                          displaced_axion_5
////       28   29  30  31                          displaced_axion_6
////       32   33  34  35                          displaced_axion_7
////       36   37  38  39                          displaced_axion_8
////       40   41  42  43                          displaced_axion_9
////       44   45  46  47                          displaced_axion_10
//    os.setf(std::ios::fixed);
//    os.precision(2);
//
//    for (int i=0, j=0; (i<12 && j <= 44); i++, j+=4)
//    {
//        os << Samples[i];
//
//        os  << " & " << *Nodes[j].GetResultPtr<ULong64_t>()
//            << R"--($\, \pm \,$)--" << sqrt(*Nodes[j].GetResultPtr<ULong64_t>())
//            << " & " << *Nodes[j+1].GetResultPtr<ULong64_t>()
//            << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+1].GetResultPtr<ULong64_t>())
//            << " & " << *Nodes[j+2].GetResultPtr<ULong64_t>()
//            << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+2].GetResultPtr<ULong64_t>())
//            << " & " << *Nodes[j+3].GetResultPtr<ULong64_t>()
//            << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+3].GetResultPtr<ULong64_t>())
//            << R"--( \\ \hline )--" << '\n';
//
//    }
//
//    os << R"--(\end{tabular}})--" << '\n';
//    std::cout << os.str() << '\n';
//}
//
//
//void Table11_Displaced_Axions()
//{
//    std::vector<std::vector<std::string>> input_filenames = {
//        //Signal
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root"},
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root"},
//        //Displaced Signal
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_FewMassPoints.root"},
//    };
//
//    constexpr std::array<const char*,15> Samples = {R"--(Prompt Signal $m_{\text{A}}$ = 1 GeV)--", R"--(Prompt Signal $m_{\text{A}}$ = 5 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 0.2 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 0.5 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 1 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 3 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 5 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 10 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 10.1 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 20 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 20.1 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 29.5 GeV)--", };
//
//    std::vector<double> massPoints = {0.2,0.5,1,3,5,10,10.1,20,20.1,29.5};
//
//    std::vector<ROOT::RDF::RResultHandle> Nodes;
//
//    std::ostringstream os, ss;
//
//    os.setf(std::ios::fixed);
//    os.precision(2);
//
//    os << R"--(\section*{Table 11 Prompt and Displaced Signal Samples})--" << '\n';
//    os << R"--(\hspace{-3cm}\scalebox{0.55}{)--" << '\n';
//    os << R"--(\begin{tabular}{|c|c|c|c|c|c|c|})--" << '\n';
//    os << R"--(\hline)--" << '\n';
//    os << R"--(Sample & Total Events & pass preselection (PS) & photon $p_T$ + $\Delta R_{\gamma\gamma}$ cut & $X$ window & SR & SR-ID
//           \\ \hline )--" << '\n';
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
//    int counter = 0;
//
//    for (auto& file: input_filenames)
//    {
//        SchottDataFrame df(MakeRDF(file, 8));
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
//        auto two_leptons = trigger_selection.Define("di_electrons",
//        [](RVec<Electron> electrons)
//        {
//            electrons.erase(std::remove_if(electrons.begin(),electrons.end(),
//            [](Electron& ep)
//            {
//                return ((ep.electron_pt/1e3 <= 20) || (std::abs(ep.electron_eta) >= 2.37)
//                        || ((1.37 < std::abs(ep.electron_eta)) && (std::abs(ep.electron_eta) < 1.52))
//                        || (ep.electron_id_medium != 1));
//
//
//            }), electrons.end());
//
//            return electrons;
//
//        },{"electrons"}).Filter([](RVec<Muon>& muons, RVec<Electron> di_electrons)
//        {
//            return (di_electrons.size()==2 && muons.empty() && DeltaR(di_electrons[0].Vector(), di_electrons[1].Vector()) > 0.01);
//
//        }, {"muons", "di_electrons"});
//
//        auto opp_charge = two_leptons.Filter([](RVec<Electron> di_electrons)
//        {
//            return (di_electrons[0].electron_charge*di_electrons[1].electron_charge < 0);
//
//        }, {"di_electrons"});
//
//        auto leading_pt = opp_charge.Filter([](RVec<Electron>& electrons)
//        {
//            return ((electrons[0].electron_pt >= 20e3 && electrons[1].electron_pt >= 27e3) || (electrons[1].electron_pt >= 20e3 && electrons[0].electron_pt >= 27e3));
//        }, {"di_electrons"});
//
////        auto delta_R = leading_pt.Filter([] (RVec<Electron>& electrons)
////        {
////            return (DeltaR(electrons[0].Vector(), electrons[1].Vector()) > 0.01);
////        }, {"di_electrons"});
//
//        auto same_flavour = leading_pt.Filter([] (RVec<Electron>& electrons)
//        {
//            return true; //std::abs(electrons[0].electron_pdg_id) == std::abs(electrons[1].electron_pdg_id) == 11;
//        }, {"di_electrons"});
//
//        auto dilep = same_flavour.Define("dilep",[] (RVec<Electron>& electrons)
//        {
//            return (electrons[0].Vector() + electrons[1].Vector());
//        }, {"di_electrons"});
//
//        auto mass = dilep.Filter([] (PtEtaPhiEVector& dilep)
//        {
//            auto mass = dilep.M()/1e3;
//            return ((mass >= 81) && (mass <= 101));
//        }, {"dilep"});
//
//        auto ptCut = mass.Filter([] (PtEtaPhiEVector& dilep)
//        {
//            auto pT = dilep.Pt()/1e3;
//            return pT > 10;
//        }, {"dilep"});
//
//        auto photonPtDeltaR = ptCut.Define("photonPtDeltaR",
//        [&](RVec<Photon> photons)
//        {
//            photons.erase(std::remove_if(photons.begin(),photons.end(),
//            [](Photon& x)
//            {
//                return ((std::abs(x.photon_eta) >= 2.37) || (x.photon_pt <= 10e3) || (std::abs(x.photon_eta) > 1.37 && std::abs(x.photon_eta) < 1.52));
//
//            }), photons.end());
//            return photons;
//        }, {"photons"}).Filter(
//       [&](RVec<Photon>& reco_photons_matched)
//       {
//            if (reco_photons_matched.size() < 2)
//            {
//              return false;
//            }
//            RVec<Photon> x;
//            auto combs = Combinations(reco_photons_matched, 2);
//            size_t length = combs[0].size();
//            double delta_r, m, pt, X, best_X, pt1, pt2, chosen_delta_r;
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
//                    x = {reco_photons_matched[combs[0][i]], reco_photons_matched[combs[1][i]]};
//                }
//            }
//            return (chosen_delta_r < 1.5 && pt1 > 10e3 && pt2 > 10e3);
//
//       }, {"photonPtDeltaR"});
//
//        auto X_window = photonPtDeltaR.Define("chosen_two",
//        [](RVec<Photon>& reco_photons_matched)
//        {
//            RVec<Photon> x;
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
//                    pt1 = reco_photons_matched[combs[0][i]].photon_pt;
//                    pt2 = reco_photons_matched[combs[1][i]].photon_pt;
//                    chosen_delta_r = delta_r;
//                    x = {reco_photons_matched[combs[0][i]], reco_photons_matched[combs[1][i]]};
//                }
//            }
//            if (pt1 > 10e3 && pt2 > 10e3 && best_X > 0.96 && best_X < 1.2 && chosen_delta_r < 1.5)
//            {
//                return x;
//            }
//            x.clear();
//            return x;
//        }, {"photonPtDeltaR"}).Filter(
//        [&](RVec<Photon>& reco_photons_matched)
//        {
//            return (reco_photons_matched.size()==2);
//        }, {"chosen_two"});
//
//        auto SR = X_window.Filter(
//        [](RVec<Photon>& photons, RVec<Electron>& electrons)
//        {
//            PtEtaPhiEVector photonVec = photons[0].Vector() + photons[1].Vector();
//            PtEtaPhiEVector electronVec = electrons[0].Vector() + electrons[1].Vector();
//            auto mass = (photonVec+electronVec).M()/1e3;
//            return ((mass >= 110) && (mass <= 140));
//        },{"chosen_two","di_electrons"});
//
//        auto SR_ID = SR.Filter(
//        [](RVec<Photon>& photons)
//        {
//            return (photons[0].photon_id_loose && photons[1].photon_id_loose);
//        },{"chosen_two"});
//
//        if (counter < 2)
//        {
//            Nodes.push_back(trigger_selection.Count());
//            Nodes.push_back(ptCut.Count());
//            Nodes.push_back(photonPtDeltaR.Count());
//            Nodes.push_back(X_window.Count());
//            Nodes.push_back(SR.Count());
//            Nodes.push_back(SR_ID.Count());
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
//                auto mass_point_photonPtDeltaR = photonPtDeltaR.Filter([&]
//                (float axion_mass)
//                {
//                    return (roundToOneDecimalPlace(axion_mass) == mass_point);
//
//                }, {"axion_masses"});
//
//                auto mass_point_X_window = X_window.Filter([&]
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
//                Nodes.push_back(mass_point_trigger_selection.Count());
//                Nodes.push_back(mass_point_ptCut.Count());
//                Nodes.push_back(mass_point_photonPtDeltaR.Count());
//                Nodes.push_back(mass_point_X_window.Count());
//                Nodes.push_back(mass_point_SR.Count());
//                Nodes.push_back(mass_point_SR_ID.Count());
//            }
//        }
//
//        counter++;
//    }
//
//    ROOT::RDF::RunGraphs(Nodes); // running all computation nodes concurrently
//
////    0    1   2   3   4   5                            ma1
////    6    7   8   9   10  11                           ma5
////    12   13  14  15  16  17                           displaced_axion_1
////    18   19  20  21  22  23                           displaced_axion_2
////    24   25  26  27  28  29                           displaced_axion_3
////    30   31  32  33  34  35                           displaced_axion_4
////    36   37  38  39  40  41                           displaced_axion_5
////    42   43  44  45  46  47                           displaced_axion_6
////    48   49  50  51  52  53                           displaced_axion_7
////    54   55  56  57  58  59                           displaced_axion_8
////    60   61  62  63  64  65                           displaced_axion_9
////    66   67  68  69  70  71                           displaced_axion_10
//
//    std::cout << R"--(\section*{Table 11 Signal Ratios})--" << '\n';
//    std::cout << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
//    std::cout << R"--(\begin{tabular}{|c|c|c|c|c|c|c|})--" << '\n';
//    std::cout << R"--(\hline)--" << '\n';
//    std::cout << R"--(Sample & $\frac{\text{pass preselection (PS)}}{\text{pass preselection (PS)}}$ & $\frac{\text{photon } p_T \text{ + } \Delta R_{\gamma\gamma} \text{ cut }}{\text{pass preselection (PS)}}$ & $\frac{X \text{ window}}{\text{pass preselection (PS)}}$ & $\frac{\text{SR}}{\text{pass preselection (PS)}}$ & $\frac{\text{SR-ID}}{\text{pass preselection (PS)}}$
//           \\ \hline )--" << '\n';
//
//    ss << R"--(\section*{Table 11 Signal Ratios})--" << '\n';
//    ss << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
//    ss << R"--(\begin{tabular}{|c|c|c|c|c|c|c|})--" << '\n';
//    ss << R"--(\hline)--" << '\n';
//    ss << R"--(Sample & $\frac{\text{pass preselection (PS)}}{\text{pass preselection (PS)}}$ & $\frac{\text{photon } p_T \text{ + } \Delta R_{\gamma\gamma} \text{ cut }}{\text{pass preselection (PS)}}$ & $\frac{X \text{ window}}{\text{photon } p_T \text{ + } \Delta R_{\gamma\gamma} \text{ cut }}$ & $\frac{\text{SR}}{X \text{ window}}$ & $\frac{\text{SR-ID}}{\text{SR}}$
//           \\ \hline )--" << '\n';
//
//    for (int i=0, j=0; (i<12 && j <= 66); i++, j+=6)
//    {
//        os << Samples[i];
//
//        os  << " & " << *Nodes[j].GetResultPtr<ULong64_t>()
//            << R"--($\, \pm \,$)--" << sqrt(*Nodes[j].GetResultPtr<ULong64_t>())
//            << " & " << *Nodes[j+1].GetResultPtr<ULong64_t>()
//            << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+1].GetResultPtr<ULong64_t>())
//            << " & " << *Nodes[j+2].GetResultPtr<ULong64_t>()
//            << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+2].GetResultPtr<ULong64_t>())
//            << " & " << *Nodes[j+3].GetResultPtr<ULong64_t>()
//            << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+3].GetResultPtr<ULong64_t>())
//            << " & " << *Nodes[j+4].GetResultPtr<ULong64_t>()
//            << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+4].GetResultPtr<ULong64_t>())
//            << " & " << *Nodes[j+5].GetResultPtr<ULong64_t>()
//            << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+5].GetResultPtr<ULong64_t>())
//            << R"--( \\ \hline )--" << '\n';
//
//        if (i==0) //1 GeV
//        {
//            std::cout << Samples[i] << " (Me) & " << 1 << " & " <<
//            static_cast<double>(*Nodes[j+2].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
//            << " & " <<
//            static_cast<double>(*Nodes[j+3].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
//            << " & " <<
//            static_cast<double>(*Nodes[j+4].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
//            << " & " <<
//            static_cast<double>(*Nodes[j+5].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
//            << R"--( \\ \hline )--" << '\n';
//
//            std::cout << Samples[i] << " (Paper) & " << 1 << " & " <<
//            2343.65 / 19516.87
//            << " & " <<
//            2204.22 / 19516.87
//            << " & " <<
//            2076.73 / 19516.87
//            << " & " <<
//            333.21 / 19516.87
//            << R"--( \\ \hline )--" << '\n';
//
//            ss << Samples[i] << " (Me) & " << 1 << " & " <<
//            static_cast<double>(*Nodes[j+2].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
//            << " & " <<
//            static_cast<double>(*Nodes[j+3].GetResultPtr<ULong64_t>()) / *Nodes[j+2].GetResultPtr<ULong64_t>()
//            << " & " <<
//            static_cast<double>(*Nodes[j+4].GetResultPtr<ULong64_t>()) / *Nodes[j+3].GetResultPtr<ULong64_t>()
//            << " & " <<
//            static_cast<double>(*Nodes[j+5].GetResultPtr<ULong64_t>()) / *Nodes[j+4].GetResultPtr<ULong64_t>()
//            << R"--( \\ \hline )--" << '\n';
//
//            ss << Samples[i] << " (Paper) & " << 1 << " & " <<
//            2343.65 / 19516.87
//            << " & " <<
//            2204.22 / 2343.65
//            << " & " <<
//            2076.73 / 2204.22
//            << " & " <<
//            333.21 / 2076.73
//            << R"--( \\ \hline )--" << '\n';
//
//        }
//
//        else if (i==1) //5 GeV
//        {
//            std::cout << Samples[i] << " (Me) & " << 1 << " & " <<
//            static_cast<double>(*Nodes[j+2].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
//            << " & " <<
//            static_cast<double>(*Nodes[j+3].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
//            << " & " <<
//            static_cast<double>(*Nodes[j+4].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
//            << " & " <<
//            static_cast<double>(*Nodes[j+5].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
//            << R"--( \\ \hline )--" << '\n';
//
//            std::cout << Samples[i] << " (Paper) & " << 1 << " & " <<
//            3245.86 / 19536.68
//            << " & " <<
//            3052.81 / 19536.68
//            << " & " <<
//            2941.28 / 19536.68
//            << " & " <<
//            1959.68 / 19536.68
//            << R"--( \\ \hline )--" << '\n';
//
//            ss << Samples[i] << " (Me) & " << 1 << " & " <<
//            static_cast<double>(*Nodes[j+2].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
//            << " & " <<
//            static_cast<double>(*Nodes[j+3].GetResultPtr<ULong64_t>()) / *Nodes[j+2].GetResultPtr<ULong64_t>()
//            << " & " <<
//            static_cast<double>(*Nodes[j+4].GetResultPtr<ULong64_t>()) / *Nodes[j+3].GetResultPtr<ULong64_t>()
//            << " & " <<
//            static_cast<double>(*Nodes[j+5].GetResultPtr<ULong64_t>()) / *Nodes[j+4].GetResultPtr<ULong64_t>()
//            << R"--( \\ \hline )--" << '\n';
//
//            ss << Samples[i] << " (Paper) & " << 1 << " & " <<
//            3245.86 / 19536.68
//            << " & " <<
//            3052.81 / 3245.86
//            << " & " <<
//            2941.28 / 3052.81
//            << " & " <<
//            1959.68 / 2941.28
//            << R"--( \\ \hline )--" << '\n';
//        }
//
//        else
//        {
//            std::cout << Samples[i] << " Displaced (Me) & " << 1 << " & " <<
//            static_cast<double>(*Nodes[j+2].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
//            << " & " <<
//            static_cast<double>(*Nodes[j+3].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
//            << " & " <<
//            static_cast<double>(*Nodes[j+4].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
//            << " & " <<
//            static_cast<double>(*Nodes[j+5].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
//            << R"--( \\ \hline )--" << '\n';
//
//            ss << Samples[i] << " Displaced (Me) & " << 1 << " & " <<
//            static_cast<double>(*Nodes[j+2].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
//            << " & " <<
//            static_cast<double>(*Nodes[j+3].GetResultPtr<ULong64_t>()) / *Nodes[j+2].GetResultPtr<ULong64_t>()
//            << " & " <<
//            static_cast<double>(*Nodes[j+4].GetResultPtr<ULong64_t>()) / *Nodes[j+3].GetResultPtr<ULong64_t>()
//            << " & " <<
//            static_cast<double>(*Nodes[j+5].GetResultPtr<ULong64_t>()) / *Nodes[j+4].GetResultPtr<ULong64_t>()
//            << R"--( \\ \hline )--" << '\n';
//        }
//
//    }
//
//    std::cout << R"--(\end{tabular}})--" << "\n\n\n";
//
//    ss << R"--(\end{tabular}})--" << "\n\n\n";
//    std::cout << ss.str() << '\n';
//
//    os << R"--(\end{tabular}})--" << '\n';
//    std::cout << os.str() << '\n';
//}

void Coupling_and_Systematics_resolved(std::unordered_map<float, float>& resolved_prompt = a, std::unordered_map<float, std::unordered_map<float, float>>& resolved_long_lived = b, std::unordered_map<float, float>& resolved_prompt_N = e, std::unordered_map<float, std::unordered_map<float, float>>& resolved_long_lived_N = f)
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
        "EG_RESOLUTION_ALL__1up",
        "EG_SCALE_ALL__1up",
        "PH_EFF_ISO_Uncertainty__1up",
        "PH_EFF_ID_Uncertainty__1up",
        "PH_EFF_TRIGGER_Uncertainty__1up",

        "EG_RESOLUTION_ALL__1down",
        "EG_SCALE_ALL__1down",
        "PH_EFF_ISO_Uncertainty__1down",
        "PH_EFF_ID_Uncertainty__1down",
        "PH_EFF_TRIGGER_Uncertainty__1down",
    };

    std::vector<std::vector<std::string>> input_filenames =
    {
        //Prompt Signal
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root"}, //1 GeV
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root"}, //5 GeV
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600751.PhPy8EG_AZNLO_ggH125_mA2p0_v1.root"}, // 2 GeV
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600752.PhPy8EG_AZNLO_ggH125_mA3p0_v1.root"}, // 3 GeV
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600756.PhPy8EG_AZNLO_ggH125_mA9p0_v1.root"}, // 9 GeV
        //Displaced Signal
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_FewMassPoints.root"},
        {
            "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_FewMassPoints_PhotonSFs.root"
        },      //C_{ayy} = 1
        {
            "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_MultiMass_Cayy0p01.root",
            "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_mA1p0_Cayy0p01.root"
        },      //C_{ayy} = 0.01
        {
            "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_MultiMass_Cayy0p001.root",
            "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_mA1p0_Cayy0p001.root"
        },     //C_{ayy} = 0.001
    };

    std::vector<float> massPoints = {0.2, 0.5, 1, 2, 3, 4, 5, 7, 8.1, 10, 15, 20, 25, 25.2, 29.5, 30,}; //mass points for displaced samples
    std::vector<float> massPoints_prompt = {1, 5, 2, 3, 9};
    std::vector<float> ALP_photon_couplings = {1.0f, 0.01f, 0.001f};
    std::vector<RResultMap<float>> resultmaps, prompt_resultmaps; //long-lived and prompt resultmaps for sys variations
    std::vector<ROOT::RDF::RResultHandle> Nodes; //Hold nominal nodes so we can use RunGraphs to run the computational graph from multiple RDF objects (samples)

    //Function that returns the particles from `startParticles` that have a parent with `targetBarcode`, given a `truthChain` to search through
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

    int counter = 0; //keep track of each file

    for (auto& file: input_filenames)
    {
        SchottDataFrame df(MakeRDF(file, 8)); //Starting dataframe, I'm running with 8 threads

        //New dataframe node that has additional columns
        auto newDf = df.Define("EventWeight", //mc generator weight
        [](const RVec<float>& ei_event_weights_generator)
        {
            return  ((ei_event_weights_generator[0]) ? 1 / ei_event_weights_generator[0] : 1);

        }, {"ei_event_weights_generator"})
        .Define("totEventWeight", //vector of efficiencies, product of photon_id_eff, photon_iso_eff, and photon_trg_eff
        [](RVec<float> photon_id_eff, RVec<float> photon_iso_eff, RVec<float> photon_trg_eff/*, RVec<float> ei_event_weights_generator*/)
        {
            auto ResizeVal = std::max({photon_id_eff.size(), photon_iso_eff.size(), photon_trg_eff.size()});
            photon_id_eff.resize(ResizeVal,1);
            photon_iso_eff.resize(ResizeVal,1);
            photon_trg_eff.resize(ResizeVal,1);

            return photon_id_eff*photon_iso_eff*photon_trg_eff; //Returns a vector of efficiencies

        }, {"photon_id_eff", "photon_iso_eff", "photon_trg_eff",/* "ei_event_weights_generator"*/})
        .Define("totEventWeightFactor", //product of vector above
        [](RVec<float>& photon_efficiencies)
        {
            float total = 1.0f;

            for (auto i: photon_efficiencies)
            {
                total *= i;
            }

            return total; //one efficiency

        }, {"totEventWeight"})
        .Define("truth_axions", //truth ALPs in each event, should only be 1 ALP per event
        [&](RVec<TruthParticle> truth_particles)
        {
            truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
            [](TruthParticle& x)
            {
                return (std::abs(x.mc_pdg_id) != 36 && std::abs(x.mc_pdg_id) != 35); //in prompt samples, ALP pdg id = 35, in long-lived, ALP pdg id = 36

            }), truth_particles.end());

            return truth_particles;

        }, {"truth_particles"})
        .Define("axion_masses", //masses of the truth_axions, should only be 1 ALP per event
        [&](RVec<TruthParticle>& truth_axions)
        {
            for (auto& particle: truth_axions)
            {
                if (counter <= 4 && particle.mc_pdg_id == 35) //then it's a prompt ALP
                {
                    return particle.mc_mass/1e3f;
                }
                else if (counter > 4 && particle.mc_pdg_id == 36) //then it's a displaced ALP
                {
                    return particle.mc_mass/1e3f;
                }
            }

            return 0.0f;

        }, {"truth_axions"});

        //new dataframe node: contains only the events of newDf that pass the trigger cut
        auto trigger_selection = newDf
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
        [](const RVec<Muon>& muons, RVec<AbstractParticle> electrons)
        {
            return (electrons.size()==2 && muons.empty()); //keep events which have exactly 2 electrons for di_electrons and no muons

        }, {"muons", "di_electrons"});

        //new dataframe node: contains only the events from `two_leptons` whose electrons in the `di_electrons` branch have opposite charge
        auto opp_charge = two_leptons.Filter([](const RVec<AbstractParticle>& electrons)
        {
            return (electrons[0].electron_charge*electrons[1].electron_charge < 0);

        }, {"di_electrons"});

        //new dataframe node: contains only the events from `opp_charge` that have 1 electron with pt > 20 GeV and the other with pt > 27 GeV
        auto leading_pt = opp_charge.Filter([](const RVec<AbstractParticle>& electrons)
        {
            return ((electrons[0].electron_pt > 20e3 && electrons[1].electron_pt > 27e3) || (electrons[1].electron_pt > 20e3 && electrons[0].electron_pt > 27e3));
        }, {"di_electrons"});

        //new dataframe node: nothing changes here from `leading_pt`...
        auto same_flavour = leading_pt.Filter([] (const RVec<AbstractParticle>& electrons)
        {
            return true; //because true is always returned in this function
        }, {"di_electrons"});

        //new dataframe node: Contains a new column `dilep` in addition to the ones in `same_flavour` that stores the di-electron four-vector
        auto dilep = same_flavour.Define("dilep",[] (RVec<AbstractParticle>& electrons)
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

        //end pre-selection -----------

        //photon acceptance and id_loose cuts
        auto photonPtDeltaR = ptCut
        .Define("photons_pass_cut_indices", //new column that contains the good photon indices
        [&](RVec<AbstractParticle>& p) //p = photon
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

        }, {"abstract_photons"})
        .Define("photonPtDeltaR", //new column that contains the good photons corresponding to the good photon indices from above
        [&](RVec<AbstractParticle>& photons, RVec<int>& x)
        {
            return Take(photons, x); //Taking only the photons that passed the cuts in each event

        }, {"abstract_photons", "photons_pass_cut_indices"})
        .Define("photonPtDeltaR_photons_pass_cut_indices", //new column for indices of the `photonPtDeltaR` photons that pass the photonPtDeltaR part of the resolved category
        [&](RVec<AbstractParticle>& reco_photons_matched)
        {
            RVec<unsigned long> x; //vector of indices
            if (reco_photons_matched.size() < 2)
            {
                return x;
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
                    x = {combs[0][i], combs[1][i]};
                }
            }
            if (chosen_delta_r < 1.5 && pt1 > 10e3 && pt2 > 10e3) //two photons corresponding to best_X must both have p_T > 10 GeV and ΔR < 1.5
            {
                return x;
            }

            x.clear();
            return x;

        }, {"photonPtDeltaR"}).Filter(//keep only events that have passed the first part of the resolved category
        [&](RVec<unsigned long>& indices)
        {
            return (indices.size()==2);

        }, {"photonPtDeltaR_photons_pass_cut_indices"})
        .Define("photonPtDeltaR_TotalEventWeightFactor", //New column: weight factor for events in RDF `photonPtDeltaR`
        [](RVec<unsigned long>& photonPtDeltaR_photons_pass_cut_indices, RVec<int>& photons_pass_cut_indices, RVec<float>& photon_efficiencies)
        {
            //First, we take the elements from photon_efficiencies that correspond to the `photons_pass_cut_indices` defined
            //earlier. Then, from that resulting vector, we take the elements corresponding to the
            //`photonPtDeltaR_photons_pass_cut_indices` defined above
            RVec<float> photonPtDeltaR_photon_efficiencies = Take(Take(photon_efficiencies, photons_pass_cut_indices), photonPtDeltaR_photons_pass_cut_indices);
            //Now, multiply all of the elements of the vector we defined above
            float total = 1.0f;
            for (auto i: photonPtDeltaR_photon_efficiencies)
            {
                total *= i;
            }
            //and return the result
            return total;
        }, {"photonPtDeltaR_photons_pass_cut_indices", "photons_pass_cut_indices", "totEventWeight"});

        //last part of the resolved category
        auto X_window = photonPtDeltaR
        .Define("X_window_photons_pass_cut_indices", //new column for indices of the `photonPtDeltaR` photons that pass the photonPtDeltaR and X_window part of the resolved category
        [](RVec<AbstractParticle>& reco_photons_matched)
        {
            RVec<unsigned long> x; //vector of indices
            if (reco_photons_matched.size() < 2)
            {
                return x;
            }
            auto combs = Combinations(reco_photons_matched, 2); //all combinations of 2 reco-photons
            size_t length = combs[0].size(); //all combinations
            double delta_r, m, pt, X, best_X, pt1, pt2, chosen_delta_r;

            for (size_t i=0; i<length; i++)
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
                    x = {combs[0][i], combs[1][i]};
                }
            }
            if (pt1 > 10e3 && pt2 > 10e3 && best_X > 0.96 && best_X < 1.2 && chosen_delta_r < 1.5) //two photons corresponding to best_X must both have p_T > 10 GeV, ΔR < 1.5, and 0.96 < best_X < 1.2
            {
                return x;
            }
            x.clear();
            return x;
        }, {"photonPtDeltaR"}).Filter(//keep only events that have passed the resolved category
        [&](RVec<unsigned long>& X_window_photons_pass_cut_indices)
        {
            return (X_window_photons_pass_cut_indices.size()==2);

        }, {"X_window_photons_pass_cut_indices"})
        .Define("chosen_two", //New column: consists of the good photons corresponding to the `X_window_photons_pass_cut_indices` defined above
        [&](RVec<AbstractParticle>& reco_photons_matched, RVec<unsigned long>& X_window_photons_pass_cut_indices)
        {
            return Take(reco_photons_matched, X_window_photons_pass_cut_indices);

        }, {"photonPtDeltaR", "X_window_photons_pass_cut_indices"})
        .Define("X_window_TotalEventWeightFactor", //New column: weight factor for events in RDF `X_window`
        [&](RVec<unsigned long>& X_window_photons_pass_cut_indices, RVec<int>& photons_pass_cut_indices, RVec<float>& photon_efficiencies)
        {
            //First, we take the elements from photon_efficiencies that correspond to the `photons_pass_cut_indices` defined
            //earlier. Then, from that resulting vector, we take the elements corresponding to the
            //`X_window_photons_pass_cut_indices` defined above
            RVec<float> X_window_photon_efficiencies = Take(Take(photon_efficiencies, photons_pass_cut_indices), X_window_photons_pass_cut_indices);
            //Now, multiply all of the elements of the vector we defined above
            float total = 1.0f;
            for (auto i: X_window_photon_efficiencies)
            {
                total *= i;
            }
            //and return the result
            return total;
        }, {"X_window_photons_pass_cut_indices", "photons_pass_cut_indices", "totEventWeight"});

        //New dataframe node: Take only the events in `X_window` that have di-electron+di-photon invariant mass between 110 and 140 GeV
        auto SR = X_window.Filter(
        [](RVec<AbstractParticle>& photons, RVec<AbstractParticle>& electrons)
        {
            PtEtaPhiEVector photonVec = photons[0].PhotonVector() + photons[1].PhotonVector();
            PtEtaPhiEVector electronVec = electrons[0].ElectronVector() + electrons[1].ElectronVector();
            auto mass = (photonVec+electronVec).M()/1e3;
            return ((mass >= 110) && (mass <= 140));
        },{"chosen_two","di_electrons"});

        //New dataframe node: Take only events in `SR` that have both of the `chosen_two` photons satisfying a loose id criteria
        auto SR_ID = SR.Filter(
        [](RVec<AbstractParticle>& photons)
        {
            return (photons[0].photon_id_loose && photons[1].photon_id_loose);
        },{"chosen_two"});

        if (counter < 5) //or however many prompt signal samples there are, each file has 1 ALP mass point
        {
            //These RResultHandles aren't used, per-se, but they allow us to call RunGraphs on the `Nodes` vector so that all of the RDataFrame graphs are run concurrently
            Nodes.push_back(newDf.Count());
            Nodes.push_back(trigger_selection.Count());
            Nodes.push_back(ptCut.Count());
            Nodes.push_back(photonPtDeltaR.Count());
            Nodes.push_back(X_window.Count());
            Nodes.push_back(SR.Count());
            Nodes.push_back(SR_ID.Count());

            //prompt_resultmaps is the vector of RResultMaps that holds the
            //nominal and variations for the following actions

            //VariationsFor registers systematic variations, `Sum` calculates the sum of the
            //event weights for each of the RDF nodes that we are interested in below:
            prompt_resultmaps.push_back(VariationsFor(newDf.Sum<float>("totEventWeightFactor")));
            prompt_resultmaps.push_back(VariationsFor(trigger_selection.Sum<float>("totEventWeightFactor")));
            prompt_resultmaps.push_back(VariationsFor(ptCut.Sum<float>("totEventWeightFactor")));
            prompt_resultmaps.push_back(VariationsFor(photonPtDeltaR.Sum<float>("photonPtDeltaR_TotalEventWeightFactor")));
            prompt_resultmaps.push_back(VariationsFor(X_window.Sum<float>("X_window_TotalEventWeightFactor")));
            prompt_resultmaps.push_back(VariationsFor(SR.Sum<float>("X_window_TotalEventWeightFactor")));
            prompt_resultmaps.push_back(VariationsFor(SR_ID.Sum<float>("X_window_TotalEventWeightFactor")));
        }
        else //long-lived samples, each file has multiple mass points
        {
            for (auto& mass_point: massPoints) //loop over mass points
            {
                //create versions of the RDF nodes that contain only the events
                //with ALP mass `mass_point`
                auto mass_point_newDf = newDf.Filter([&]
                (float axion_mass)
                {
                    return (roundToOneDecimalPlace(axion_mass) == mass_point);

                }, {"axion_masses"});

                auto mass_point_trigger_selection = trigger_selection.Filter([&]
                (float axion_mass)
                {
                    return (roundToOneDecimalPlace(axion_mass) == mass_point);

                }, {"axion_masses"});

                auto mass_point_ptCut = ptCut.Filter([&]
                (float axion_mass)
                {
                    return (roundToOneDecimalPlace(axion_mass) == mass_point);

                }, {"axion_masses"});

                auto mass_point_photonPtDeltaR = photonPtDeltaR.Filter([&]
                (float axion_mass)
                {
                    return (roundToOneDecimalPlace(axion_mass) == mass_point);

                }, {"axion_masses"});

                auto mass_point_X_window = X_window.Filter([&]
                (float axion_mass)
                {
                    return (roundToOneDecimalPlace(axion_mass) == mass_point);

                }, {"axion_masses"});

                auto mass_point_SR = SR.Filter([&]
                (float axion_mass)
                {
                    return (roundToOneDecimalPlace(axion_mass) == mass_point);

                }, {"axion_masses"});

                auto mass_point_SR_ID = SR_ID.Filter([&]
                (float axion_mass)
                {
                    return (roundToOneDecimalPlace(axion_mass) == mass_point);

                }, {"axion_masses"});

                //These RResultHandles aren't used, per-se, but they allow us to call RunGraphs on the `Nodes` vector so that all of the RDataFrame graphs are run concurrently
                Nodes.push_back(mass_point_newDf.Count());
                Nodes.push_back(mass_point_trigger_selection.Count());
                Nodes.push_back(mass_point_ptCut.Count());
                Nodes.push_back(mass_point_photonPtDeltaR.Count());
                Nodes.push_back(mass_point_X_window.Count());
                Nodes.push_back(mass_point_SR.Count());
                Nodes.push_back(mass_point_SR_ID.Count());

                //prompt_resultmaps is the vector of RResultMaps that holds the
                //nominal and variations for the following actions

                //VariationsFor registers systematic variations, `Sum` calculates the sum of the
                //event weights for each of the RDF nodes that we are interested in below:
                resultmaps.push_back(VariationsFor(mass_point_newDf.Sum<float>("totEventWeightFactor")));
                resultmaps.push_back(VariationsFor(mass_point_trigger_selection.Sum<float>("totEventWeightFactor")));
                resultmaps.push_back(VariationsFor(mass_point_ptCut.Sum<float>("totEventWeightFactor")));
                resultmaps.push_back(VariationsFor(mass_point_photonPtDeltaR.Sum<float>("photonPtDeltaR_TotalEventWeightFactor")));
                resultmaps.push_back(VariationsFor(mass_point_X_window.Sum<float>("X_window_TotalEventWeightFactor")));
                resultmaps.push_back(VariationsFor(mass_point_SR.Sum<float>("X_window_TotalEventWeightFactor")));
                resultmaps.push_back(VariationsFor(mass_point_SR_ID.Sum<float>("X_window_TotalEventWeightFactor")));
            }
        }

        counter++; //keep track of which file we're on
    }

    ROOT::RDF::RunGraphs(Nodes); // running all computation nodes concurrently, this will trigger the resultmaps as well.

    //systematics' names
    std::vector<std::string> syst_indices = {"nominal", "photons_and_electrons:EG_RESOLUTION_ALL__1down", "photons_and_electrons:EG_SCALE_ALL__1down", "photon_iso_eff:PH_EFF_ISO_Uncertainty__1down", "photon_id_eff:PH_EFF_ID_Uncertainty__1down", "photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down", "photons_and_electrons:EG_RESOLUTION_ALL__1up", "photons_and_electrons:EG_SCALE_ALL__1up", "photon_iso_eff:PH_EFF_ISO_Uncertainty__1up", "photon_id_eff:PH_EFF_ID_Uncertainty__1up", "photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up"};

    TLatex Tl;
    TCanvas* c1;
    TH2D* histo2d;
    std::string name, title;
    std::unordered_map<float,float> efficiencies = {{1.0f, 1.0f}, {0.01f, 1.0f}, {0.001f, 0.284f}};

    //for each systematic
    for (auto& syst_index: syst_indices)
    {
        ///displaced case
        int start_index = 0;

        //create a plot
        c1 = new TCanvas();
        c1->SetLogy(); //Set log scale for y-axis
        title = syst_index + std::string(" (displaced resolved)") + ";m_{a} [GeV];C_{a#gamma#gamma}" + ";Efficiency";
        name = std::string("C_{a#gamma#gamma}") + syst_index;
        double y_bins[] = {0.001, 0.008, 0.0099, 0.08, 0.1, 0.8, 1, 8, 10}; //y bins

        histo2d = new TH2D(name.c_str(), title.c_str(), 299, 0.1, 31, 8, y_bins); //numbers are: x_bins, x_min, x_max, y_bins, y_min, y_max
        histo2d->GetYaxis()->SetMoreLogLabels();
        histo2d->GetYaxis()->SetNoExponent();
        histo2d->SetOption("LOGZ");
        histo2d->GetXaxis()->SetTitleSize(0.04);
        histo2d->GetYaxis()->SetTitleSize(0.04);
        histo2d->GetZaxis()->SetTitleSize(0.04);
        histo2d->GetYaxis()->SetTitleOffset(0.7);
        histo2d->GetZaxis()->SetTitleOffset(0.6);
        histo2d->GetYaxis()->CenterTitle(true);
        histo2d->GetZaxis()->SetTitle("Efficiency");
        //for each coupling
        for (auto& coupling: ALP_photon_couplings)
        {
            //i is indexing each mass point, and j is indexing the
            //group of results corresponding to mass point i
            for (int i=0, j=start_index; (i < massPoints.size()); i++, j+=7)
            {
                int bin_number = histo2d->FindFixBin(massPoints[i], coupling);
                float efficiency;

                auto keys = resultmaps[j].GetKeys(); //gets systematic variation names for the node `mass_point_newDf`
                auto keys5 = resultmaps[j+5].GetKeys(); //gets systematic variation names for the node `mass_point_SR`

                //if syst_index is a systematic variation of the total # of events and
                //the ones that passed the resolved category, then assign the efficiency
                //the corresponding value, else assign it 0
                if ((std::find(keys.begin(), keys.end(), syst_index) != keys.end())
                    &&
                    (std::find(keys5.begin(), keys5.end(), syst_index) != keys5.end()))
                {
                    efficiency = resultmaps[j][syst_index] ? resultmaps[j+5][syst_index] / resultmaps[j][syst_index] : 0.0;
                }
                else if (std::find(keys5.begin(), keys5.end(), syst_index) != keys5.end())
                {
                    efficiency = resultmaps[j]["nominal"] ? resultmaps[j+5][syst_index] / resultmaps[j]["nominal"] : 0.0;
                }
                else //then something's wrong :(
                {
                    throw std::runtime_error(std::string("Something's wrong at ") + syst_index + ": " + std::to_string(massPoints[i]) + " GeV, coupling = " + std::to_string(coupling));
                }
                
                efficiency *= efficiencies[coupling]; //mutiply efficiency by Pythia Filter Efficiency

                if (syst_index == "nominal") //nominal case: Store the efficiency in the unordered_map
                {
                    resolved_long_lived[coupling][massPoints[i]] = efficiency;
                    resolved_long_lived_N[coupling][massPoints[i]] = resultmaps[j]["nominal"];

                    if (massPoints[i] == 1)
                    {
                        std::cout << resolved_long_lived_N[coupling][massPoints[i]] << '\n';
                    }
                }
                histo2d->SetBinContent(bin_number, efficiency);
            }
            start_index += 112; //to the next displaced file (with a different coupling)
        }

        gStyle->SetPalette(1); //heat-map color style
        histo2d->Draw("COLZ"); //draw with color-bar
        gStyle->SetOptStat(0);
        Tl.SetTextSize(0.03);
        Tl.DrawLatexNDC(0.7, 0.83, "#it{ATLAS} Internal");
        Tl.DrawLatexNDC(0.7, 0.73,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
        c1->SetCanvasSize(2.5*c1->GetWw(), c1->GetWh());
        gPad->SetLeftMargin(0.06);
        title = std::string("Cayy") + "DisplacedResolved" + syst_index + ".pdf";
        c1->SaveAs(title.c_str());
        //open -a Safari C_{a#gamma#gamma}_displaced_resolved_* C_{a#gamma#gamma}_prompt_resolved_*
        delete c1;
        delete histo2d;

        ///Prompt case

        //create a new plot
        c1 = new TCanvas();
        c1->SetLogy(); //Set log scale for y-axis
        title = syst_index + std::string(" (prompt resolved)") + ";m_{a} [GeV];C_{a#gamma#gamma};Efficiency";
        name = std::string("C_{a#gamma#gamma}") + syst_index;
        double prompt_y_bins[] = {0.008, 0.0099, 0.08, 0.1}; //y bins
        histo2d = new TH2D(name.c_str(), title.c_str(), 25, 0.1, 10, 3, prompt_y_bins); //numbers are: x_bins, x_min, x_max, y_bins, y_min, y_max
        histo2d->GetXaxis()->SetTitleSize(0.04);
        histo2d->GetYaxis()->SetTitleSize(0.04);
        histo2d->GetZaxis()->SetTitleSize(0.04);
        histo2d->GetYaxis()->SetTitleOffset(0.7);
        histo2d->GetZaxis()->SetTitleOffset(0.6);
        histo2d->GetYaxis()->CenterTitle(true);
        histo2d->GetZaxis()->SetTitle("Efficiency");

        //i is indexing each mass point, and j is indexing the
        //group of results corresponding to mass point i
        for (int i=0, j=0; (i < massPoints_prompt.size()); i++, j+=7)
        {
            int bin_number = histo2d->FindFixBin(massPoints_prompt[i], ALP_photon_couplings[1]);
            float efficiency;

            auto keys = prompt_resultmaps[j].GetKeys();
            auto keys5 = prompt_resultmaps[j+5].GetKeys();

            //if syst_index is a systematic variation of the total # of events and
            //the ones that passed the resolved category, then assign the efficiency
            //the corresponding value, else assign it 0
            if ((std::find(keys.begin(), keys.end(), syst_index) != keys.end())
                &&
                (std::find(keys5.begin(), keys5.end(), syst_index) != keys5.end()))
            {
                efficiency = prompt_resultmaps[j][syst_index] ? prompt_resultmaps[j+5][syst_index] / prompt_resultmaps[j][syst_index] : 0.0;
            }
            else if (std::find(keys5.begin(), keys5.end(), syst_index) != keys5.end())
            {
                efficiency = prompt_resultmaps[j]["nominal"] ? prompt_resultmaps[j+5][syst_index] / prompt_resultmaps[j]["nominal"] : 0.0;
            }
            else //then something's wrong :(
            {
                throw std::runtime_error(std::string("Something's wrong at ") + syst_index + ": " + std::to_string(massPoints[i]) + " GeV");
            }

            if (syst_index == "nominal") //nominal case: Store the efficiency in the unordered_map
            {
                resolved_prompt[massPoints_prompt[i]] = efficiency; //efficiency
                resolved_prompt_N[massPoints_prompt[i]] = prompt_resultmaps[j]["nominal"]; //total number of events
            }

            histo2d->SetBinContent(bin_number, efficiency);
        }

        gStyle->SetPalette(1); //heat-map color style
        histo2d->Draw("COLZ"); //draw with color-bar

        gStyle->SetOptStat(0);
        Tl.SetTextSize(0.03);
        Tl.DrawLatexNDC(0.6, 0.83, "#it{ATLAS} Internal");
        Tl.DrawLatexNDC(0.6, 0.73,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
        c1->SetCanvasSize(2.5*c1->GetWw(), c1->GetWh());
        gPad->SetLeftMargin(0.06);
        title = std::string("CayyPromptResolved") + syst_index + ".pdf";
        c1->SaveAs(title.c_str());

    }

    ///loop over the results for the long-lived nominal case for different couplings to see the results
    for (auto& coupling: resolved_long_lived)
    {
        std::cout << coupling.first << ": ";
        for (auto& masspoint: coupling.second)
        {
            std::cout << '{' << masspoint.first << ','
            << masspoint.second << "}, ";
        }
        std::cout << '\n';
    }

}

void Coupling_and_Systematics_merged(std::unordered_map<float, float>& merged_prompt = c, std::unordered_map<float, std::unordered_map<float, float>>& merged_long_lived = d, std::unordered_map<float, float>& merged_prompt_N = g, std::unordered_map<float, std::unordered_map<float, float>>& merged_long_lived_N = h)
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
        "EG_RESOLUTION_ALL__1up",
        "EG_SCALE_ALL__1up",
        "PH_EFF_ISO_Uncertainty__1up",
        "PH_EFF_ID_Uncertainty__1up",
        "PH_EFF_TRIGGER_Uncertainty__1up",

        "EG_RESOLUTION_ALL__1down",
        "EG_SCALE_ALL__1down",
        "PH_EFF_ISO_Uncertainty__1down",
        "PH_EFF_ID_Uncertainty__1down",
        "PH_EFF_TRIGGER_Uncertainty__1down",
    };

    std::vector<std::vector<std::string>> input_filenames =
    {
        //Prompt Signal
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root"}, //1 GeV
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root"}, //5 GeV
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600751.PhPy8EG_AZNLO_ggH125_mA2p0_v1.root"}, // 2 GeV
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600752.PhPy8EG_AZNLO_ggH125_mA3p0_v1.root"}, // 3 GeV
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600756.PhPy8EG_AZNLO_ggH125_mA9p0_v1.root"}, // 9 GeV
        //Displaced Signal
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_FewMassPoints.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_FewMassPoints_PhotonSFs.root"}, //C_{ayy} = 1
        {
            "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_MultiMass_Cayy0p01.root",
            "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_mA1p0_Cayy0p01.root"
        },      //C_{ayy} = 0.01
        {
            "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_MultiMass_Cayy0p001.root",
            "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_mA1p0_Cayy0p001.root"
        },     //C_{ayy} = 0.001
    };

    std::vector<float> massPoints = {0.2, 0.5, 1, 2, 3, 4, 5, 7, 8.1, 10, 15, 20, 25, 25.2, 29.5, 30,}; //mass points for displaced samples
    std::vector<float> massPoints_prompt = {1, 5, 2, 3, 9};
    std::vector<float> ALP_photon_couplings = {1.0f, 0.01f, 0.001f};
    std::vector<RResultMap<float>> resultmaps, prompt_resultmaps; //long-lived and prompt resultmaps for sys variations
    std::vector<ROOT::RDF::RResultHandle> Nodes; //Hold nominal nodes so we can use RunGraphs to run the computational graph from multiple RDF objects (samples)

    //Function that returns the particles from `startParticles` that have a parent with `targetBarcode`, given a `truthChain` to search through
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

    int counter = 0; //keep track of each file

    for (auto& file: input_filenames)
    {
        SchottDataFrame df(MakeRDF(file, 8)); //Starting dataframe, I'm running with 8 threads

        //New dataframe node that has additional columns
        auto newDf = df.Define("EventWeight", //mc generator weight
        [](const RVec<float>& ei_event_weights_generator)
        {
            return  ((ei_event_weights_generator[0]) ? 1 / ei_event_weights_generator[0] : 1);

        }, {"ei_event_weights_generator"})
        .Define("totEventWeightVec", //vector of efficiencies, product of photon_id_eff, photon_iso_eff, and photon_trg_eff
        [](RVec<float> photon_id_eff, RVec<float> photon_iso_eff, RVec<float> photon_trg_eff/*, RVec<float> ei_event_weights_generator*/)
        {
            auto ResizeVal = std::max({photon_id_eff.size(), photon_iso_eff.size(), photon_trg_eff.size()});
            photon_id_eff.resize(ResizeVal,1);
            photon_iso_eff.resize(ResizeVal,1);
            photon_trg_eff.resize(ResizeVal,1);

            return photon_id_eff*photon_iso_eff*photon_trg_eff; //Returns a vector of efficiencies

        }, {"photon_id_eff", "photon_iso_eff", "photon_trg_eff",/* "ei_event_weights_generator"*/})
        .Define("totEventWeightFactor", //product of vector above
        [](RVec<float>& photon_efficiencies)
        {
            float total = 1.0f;

            for (auto i: photon_efficiencies)
            {
                total *= i;
            }

            return total; //one efficiency

        }, {"totEventWeightVec"})
        .Define("truth_axions", //truth ALPs in each event, should only be 1 ALP per event
        [&](RVec<TruthParticle> truth_particles)
        {
            truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
            [](TruthParticle& x)
            {
                return (std::abs(x.mc_pdg_id) != 36 && std::abs(x.mc_pdg_id) != 35); //in prompt samples, ALP pdg id = 35, in long-lived, ALP pdg id = 36

            }), truth_particles.end());

            return truth_particles;

        }, {"truth_particles"})
        .Define("axion_masses", //masses of the truth_axions, should only be 1 ALP per event
        [&](RVec<TruthParticle>& truth_axions)
        {
            for (auto& particle: truth_axions)
            {
                if (counter <= 4 && particle.mc_pdg_id == 35) //then it's a prompt ALP
                {
                    return particle.mc_mass/1e3f;
                }
                else if (counter > 4 && particle.mc_pdg_id == 36) //then it's a displaced ALP
                {
                    return particle.mc_mass/1e3f;
                }
            }
            return 0.0f;

        }, {"truth_axions"});

        //new dataframe node: contains only the events of newDf that pass the trigger cut
        auto trigger_selection = newDf
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
        [](const RVec<Muon>& muons, RVec<AbstractParticle> electrons)
        {
            return (electrons.size()==2 && muons.empty()); //keep events which have exactly 2 electrons for di_electrons and no muons

        }, {"muons", "di_electrons"});

        //new dataframe node: contains only the events from `two_leptons` whose electrons in the `di_electrons` branch have opposite charge
        auto opp_charge = two_leptons.Filter([](const RVec<AbstractParticle>& electrons)
        {
            return (electrons[0].electron_charge*electrons[1].electron_charge < 0);

        }, {"di_electrons"});

        //new dataframe node: contains only the events from `opp_charge` that have 1 electron with pt > 20 GeV and the other with pt > 27 GeV
        auto leading_pt = opp_charge.Filter([](const RVec<AbstractParticle>& electrons)
        {
            return ((electrons[0].electron_pt > 20e3 && electrons[1].electron_pt > 27e3) || (electrons[1].electron_pt > 20e3 && electrons[0].electron_pt > 27e3));
        }, {"di_electrons"});

        //new dataframe node: nothing changes here from `leading_pt`...
        auto same_flavour = leading_pt.Filter([] (const RVec<AbstractParticle>& electrons)
        {
            return true; //because true is always returned in this function
        }, {"di_electrons"});

        //new dataframe node: Contains a new column `dilep` in addition to the ones in `same_flavour` that stores the di-electron four-vector
        auto dilep = same_flavour.Define("dilep",[] (RVec<AbstractParticle>& electrons)
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

        //end pre-selection -----------

        //photon acceptance and id_loose cuts
        auto photon_passes_cuts = ptCut
        .Define("abstract_photons_pass_cut_indices", //new column that contains the good photon indices
        [&](RVec<AbstractParticle>& p) //p = photon
        {
            RVec<int> x; //indices of photons that pass cuts

            for (auto i = 0; i < p.size(); i++)
            {
                //keep reco-photons that have |η| < 2.37, p_T > 10 GeV, and |η| not between 1.37 and 1.52
                if (not ((std::abs(p[i].photon_eta) >= 2.37) or (p[i].photon_pt <= 10e3) or (std::abs(p[i].photon_eta) > 1.37 and std::abs(p[i].photon_eta) < 1.52)))
                {
                    x.push_back(i);
                }
            }

            return x;
        }, {"abstract_photons"})
        .Define("photons_pass_cuts", //new column that contains the good photons corresponding to the good photon indices from above
        [&](RVec<AbstractParticle>& photons, RVec<int>& x)
        {
            return Take(photons, x); //Taking only the photons that passed the cuts in each event

        }, {"abstract_photons", "abstract_photons_pass_cut_indices"});

        //New dataframe node: contains only events from `photon_passes_cuts` that fail the resolved category and pass the merged category
        auto merged_reco_photons_matched = photon_passes_cuts.Filter(
        [&](const RVec<AbstractParticle>& reco_photons_matched)
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

        }, {"photons_pass_cuts", "merged_photon_index"})
        //mpi = merged_photon_index
        .Define("totEventWeight", //New column: weight factor for events in RDF `merged_reco_photons_matched`
        [](RVec<int>& x, int mpi, RVec<float> photon_id_eff, RVec<float> photon_iso_eff, RVec<float> photon_trg_eff/*, RVec<float> ei_event_weights_generator*/)
        {
            //   ||  resize 3 vectors just in case they don't already have the same size  ||
            //   \/                                                                       \/
            auto ResizeVal = std::max({photon_id_eff.size(), photon_iso_eff.size(), photon_trg_eff.size()});
            photon_id_eff.resize(ResizeVal,1);
            photon_iso_eff.resize(ResizeVal,1);
            photon_trg_eff.resize(ResizeVal,1);

            //First, take the x indices from the respective photon_*_eff vectors,
            //this corresponds to the photons from the set of all reco photons in
            //the event that passed the "photon_passes_cuts" cuts. Then, from
            //those photons, select the index mpi element that corresponds to
            //the merged reco-photon

            return Take(photon_id_eff, x)[mpi] * Take(photon_iso_eff, x)[mpi] * Take(photon_trg_eff, x)[mpi]; // a single number

        }, {"abstract_photons_pass_cut_indices", "merged_photon_index", "photon_id_eff", "photon_iso_eff", "photon_trg_eff",/* "ei_event_weights_generator"*/})
        .Define("reconstructed_mass", //new column: dilep+merged invariant mass
        [&](const RVec<AbstractParticle>& di_electrons, const AbstractParticle& merged_photon)
        {
            auto four_momentum = di_electrons[0].ElectronVector() + di_electrons[1].ElectronVector();

            return (four_momentum + merged_photon.PhotonVector()).M()/1e3;

        }, {"di_electrons", "merged_photon"});

        //new dataframe node: keep only events from `merged_reco_photons_matched` where
        //the dilep+merged invariant mass is not between 110 and 130 GeV inclusive
        auto pSB = merged_reco_photons_matched.Filter(
        [](const double reconstructed_mass)
        {
            return (reconstructed_mass < 110) || (reconstructed_mass > 130);
        }, {"reconstructed_mass"});

        //new dataframe node: keep only events from `merged_reco_photons_matched` where
        //the dilep+merged invariant mass is between 110 and 130 GeV inclusive
        auto pSR = merged_reco_photons_matched.Filter(
        [](const double reconstructed_mass)
        {
            return (reconstructed_mass >= 110) && (reconstructed_mass <= 130);
        }, {"reconstructed_mass"});

        //new dataframe node: keep only events from `pSB` that have all E-ratios > 0.8
        auto SB = pSB.Filter(
        [](const RVec<float>& Eratio)
        {
            return (!Any(Eratio <= 0.8));
        }, {"photon_shower_shape_e_ratio"});

        //new dataframe node: keep only events from `pSR` that have all E-ratios > 0.8
        auto SR = pSR.Filter(
        [](const RVec<float>& Eratio)
        {
            return (!Any(Eratio <= 0.8));
        }, {"photon_shower_shape_e_ratio"});

        if (counter < 5) //or however many prompt signal samples there are, each file has 1 ALP mass point
        {
            //These RResultHandles aren't used, per-se, but they allow us to call RunGraphs on the `Nodes` vector so that all of the RDataFrame graphs are run concurrently
            Nodes.push_back(newDf.Count());
            Nodes.push_back(merged_reco_photons_matched.Count());
            Nodes.push_back(pSB.Count());
            Nodes.push_back(pSR.Count());
            Nodes.push_back(SB.Count());
            Nodes.push_back(SR.Count());

            prompt_resultmaps.push_back(VariationsFor(newDf.Sum<float>("totEventWeightFactor")));
            prompt_resultmaps.push_back(VariationsFor(merged_reco_photons_matched.Sum<float>("totEventWeight")));
            prompt_resultmaps.push_back(VariationsFor(pSB.Sum<float>("totEventWeight")));
            prompt_resultmaps.push_back(VariationsFor(pSR.Sum<float>("totEventWeight")));
            prompt_resultmaps.push_back(VariationsFor(SB.Sum<float>("totEventWeight")));
            prompt_resultmaps.push_back(VariationsFor(SR.Sum<float>("totEventWeight")));
        }
        else
        {
            for (auto& mass_point: massPoints)
            {
                auto mass_point_newDf = newDf.Filter([&]
                (float axion_mass)
                {
                    return (roundToOneDecimalPlace(axion_mass) == mass_point);

                }, {"axion_masses"});

                auto mass_point_merged_reco_photons_matched = merged_reco_photons_matched.Filter([&]
                (float axion_mass)
                {
                    return (roundToOneDecimalPlace(axion_mass) == mass_point);

                }, {"axion_masses"});

                auto mass_point_pSB = pSB.Filter([&]
                (float axion_mass)
                {
                    return (roundToOneDecimalPlace(axion_mass) == mass_point);

                }, {"axion_masses"});

                auto mass_point_pSR = pSR.Filter([&]
                (float axion_mass)
                {
                    return (roundToOneDecimalPlace(axion_mass) == mass_point);

                }, {"axion_masses"});

                auto mass_point_SB = SB.Filter([&]
                (float axion_mass)
                {
                    return (roundToOneDecimalPlace(axion_mass) == mass_point);

                }, {"axion_masses"});

                auto mass_point_SR = SR.Filter([&]
                (float axion_mass)
                {
                    return (roundToOneDecimalPlace(axion_mass) == mass_point);

                }, {"axion_masses"});

                Nodes.push_back(mass_point_newDf.Count());
                Nodes.push_back(mass_point_merged_reco_photons_matched.Count());
                Nodes.push_back(mass_point_pSB.Count());
                Nodes.push_back(mass_point_pSR.Count());
                Nodes.push_back(mass_point_SB.Count());
                Nodes.push_back(mass_point_SR.Count());

                resultmaps.push_back(VariationsFor(mass_point_newDf.Sum<float>("totEventWeightFactor")));
                resultmaps.push_back(VariationsFor(mass_point_merged_reco_photons_matched.Sum<float>("totEventWeight")));
                resultmaps.push_back(VariationsFor(mass_point_pSB.Sum<float>("totEventWeight")));
                resultmaps.push_back(VariationsFor(mass_point_pSR.Sum<float>("totEventWeight")));
                resultmaps.push_back(VariationsFor(mass_point_SB.Sum<float>("totEventWeight")));
                resultmaps.push_back(VariationsFor(mass_point_SR.Sum<float>("totEventWeight")));

            }
        }
        counter++;
    }

    ROOT::RDF::RunGraphs(Nodes); // running all computation nodes concurrently, this will trigger resultmaps as well.

    std::vector<std::string> syst_indices = {"nominal", "photons_and_electrons:EG_RESOLUTION_ALL__1down", "photons_and_electrons:EG_SCALE_ALL__1down", "photon_iso_eff:PH_EFF_ISO_Uncertainty__1down", "photon_id_eff:PH_EFF_ID_Uncertainty__1down", "photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down", "photons_and_electrons:EG_RESOLUTION_ALL__1up", "photons_and_electrons:EG_SCALE_ALL__1up", "photon_iso_eff:PH_EFF_ISO_Uncertainty__1up", "photon_id_eff:PH_EFF_ID_Uncertainty__1up", "photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up"};

    TLatex Tl;
    TCanvas* c1;
    TH2D* histo2d;
    std::string name, title;
    std::unordered_map<float,float> efficiencies = {{1.0f, 1.0f}, {0.01f, 1.0f}, {0.001f, 0.284f}};

    //for each systematic
    for (auto& syst_index: syst_indices)
    {
        ///displaced case
        int start_index = 0;

        //create a plot
        c1 = new TCanvas();
        c1->SetLogy(); //Set log scale for y-axis
        title = syst_index + std::string(" (displaced merged)") + ";m_{a} [GeV];C_{a#gamma#gamma};Efficiency";
        name = std::string("C_{a#gamma#gamma}") + syst_index;
        double y_bins[] = {0.001, 0.008, 0.0099, 0.08, 0.1, 0.8, 1, 8, 10}; //y bins
        histo2d = new TH2D(name.c_str(), title.c_str(), 299, 0.1, 31, 8, y_bins);
        histo2d->GetYaxis()->SetMoreLogLabels();
        histo2d->GetYaxis()->SetNoExponent();
        histo2d->SetOption("LOGZ");
        histo2d->GetXaxis()->SetTitleSize(0.04);
        histo2d->GetYaxis()->SetTitleSize(0.04);
        histo2d->GetZaxis()->SetTitleSize(0.04);

        histo2d->GetYaxis()->SetTitleOffset(0.7);
        histo2d->GetZaxis()->SetTitleOffset(0.6);

        histo2d->GetYaxis()->CenterTitle(true);
        histo2d->GetZaxis()->SetTitle("Efficiency");

        //for each coupling
        for (auto& coupling: ALP_photon_couplings)
        {
            //i is indexing each mass point, and j is indexing the
            //group of results corresponding to mass point i
            for (int i=0, j=start_index; (i < massPoints.size()); i++, j+=6)
            {
                int bin_number = histo2d->FindFixBin(massPoints[i], coupling);
                float efficiency;
                
                auto keys = resultmaps[j].GetKeys();
                auto keys5 = resultmaps[j+5].GetKeys();

    //            if syst_index is a systematic variation of the total # of events and
    //            the ones that passed the resolved category, then assign the efficiency
    //            the corresponding value
                if ((std::find(keys.begin(), keys.end(), syst_index) != keys.end())
                    &&
                    (std::find(keys5.begin(), keys5.end(), syst_index) != keys5.end()))
                {
                    efficiency = resultmaps[j][syst_index] ? resultmaps[j+5][syst_index] / resultmaps[j][syst_index] : 0.0;
                }
                else if (std::find(keys5.begin(), keys5.end(), syst_index) != keys5.end())
                {
                    efficiency = resultmaps[j]["nominal"] ? resultmaps[j+5][syst_index] / resultmaps[j]["nominal"] : 0.0;
                }
                else //then something's wrong :(
                {
                    efficiency = -1;
                }
                
                efficiency *= efficiencies[coupling]; //mutiply efficiency by Pythia Filter Efficiency

                if (syst_index == "nominal")
                {
                    merged_long_lived_N[coupling][massPoints[i]] = resultmaps[j]["nominal"];
                }
                histo2d->SetBinContent(bin_number, efficiency);
            }
            start_index += 96; //to the next displaced file (with a different coupling)
        }

        gStyle->SetPalette(1); //heat-map color style
        histo2d->Draw("COLZ"); //draw with color-bar

        gStyle->SetOptStat(0);
        Tl.SetTextSize(0.03);
        Tl.DrawLatexNDC(0.7, 0.83, "#it{ATLAS} Internal");
        Tl.DrawLatexNDC(0.7, 0.73,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
        c1->SetCanvasSize(2.5*c1->GetWw(), c1->GetWh());
        gPad->SetLeftMargin(0.06);
        title = std::string("CayyDisplacedMerged") + syst_index + ".pdf";
        c1->SaveAs(title.c_str());
        //open -a Safari C_{a#gamma#gamma}_displaced_merged_* C_{a#gamma#gamma}_prompt_merged_*
        delete c1;
        delete histo2d;

        ///prompt case

        //create a plot
        c1 = new TCanvas();
        c1->SetLogy(); //Set log scale for y-axis
        title = syst_index + std::string(" (prompt merged)") + ";m_{a} [GeV];C_{a#gamma#gamma};Efficiency";
        name = std::string("C_{a#gamma#gamma}") + syst_index;
        double prompt_y_bins[] = {0.008, 0.0099, 0.08, 0.1}; //y bins
        histo2d = new TH2D(name.c_str(), title.c_str(), 25, 0.1, 10, 3, prompt_y_bins);
        histo2d->GetXaxis()->SetTitleSize(0.04);
        histo2d->GetYaxis()->SetTitleSize(0.04);
        histo2d->GetZaxis()->SetTitleSize(0.04);
        histo2d->SetOption("LOGZ");
        histo2d->GetYaxis()->SetTitleOffset(0.7);
        histo2d->GetZaxis()->SetTitleOffset(0.6);

        histo2d->GetYaxis()->CenterTitle(true);
        histo2d->GetZaxis()->SetTitle("Efficiency");

        //i is indexing each mass point, and j is indexing the
        //group of results corresponding to mass point i
        for (int i=0, j=0; (i < massPoints_prompt.size()); i++, j+=6)
        {
            int bin_number = histo2d->FindFixBin(massPoints_prompt[i], ALP_photon_couplings[1]);
            float efficiency;

            auto keys = prompt_resultmaps[j].GetKeys();
            auto keys5 = prompt_resultmaps[j+5].GetKeys();

//            if syst_index is a systematic variation of the total # of events and
//            the ones that passed the resolved category, then assign the efficiency
//            the corresponding value
            if ((std::find(keys.begin(), keys.end(), syst_index) != keys.end())
                &&
                (std::find(keys5.begin(), keys5.end(), syst_index) != keys5.end()))
            {
                efficiency = prompt_resultmaps[j][syst_index] ? prompt_resultmaps[j+5][syst_index] / prompt_resultmaps[j][syst_index] : 0.0;
            }
            else if (std::find(keys5.begin(), keys5.end(), syst_index) != keys5.end())
            {
                efficiency = prompt_resultmaps[j]["nominal"] ? prompt_resultmaps[j+5][syst_index] / prompt_resultmaps[j]["nominal"] : 0.0;
            }
            else //then something's wrong :(
            {
                efficiency = -1;
            }

            if (syst_index == "nominal")
            {
                merged_prompt[massPoints_prompt[i]] = efficiency;
                merged_prompt_N[massPoints[i]] = resultmaps[j]["nominal"];
            }

            histo2d->SetBinContent(bin_number, efficiency);
        }

        gStyle->SetPalette(1); //heat-map color style
        histo2d->Draw("COLZ"); //draw with color-bar

        gStyle->SetOptStat(0);
        Tl.SetTextSize(0.03);
        Tl.DrawLatexNDC(0.6, 0.83, "#it{ATLAS} Internal");
        Tl.DrawLatexNDC(0.6, 0.73,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
        c1->SetCanvasSize(2.5*c1->GetWw(), c1->GetWh());
        gPad->SetLeftMargin(0.06);
        title = std::string("CayyPromptMerged") + syst_index + ".pdf";
        c1->SaveAs(title.c_str());
    }
}

void LimitPlot()
{
    // >= 1 GeV -> resolved branching ratios, prompt effciencies, long-lived efficiencies
    // <= 1 GeV -> merged branching ratios, prompt effciencies, long-lived efficiencies

    //For Coupling_and_Systematics_merged(), return 1 GeV prompt efficiency and 1 GeV long-lived efficiency
    //For Coupling_and_Systematics_resolved(), return 1, 5 GeV prompt efficiency and 1, 5 GeV long-lived efficiency

    //key is mass, value is efficiency
    std::unordered_map<float, float> resolved_prompt_eff, merged_prompt_eff, resolved_prompt_N, merged_prompt_N;
    std::unordered_map<float, std::unordered_map<float, float>> resolved_long_lived_eff, merged_long_lived_eff, resolved_long_lived_N, merged_long_lived_N;

    Coupling_and_Systematics_resolved(resolved_prompt_eff, resolved_long_lived_eff, resolved_prompt_N, resolved_long_lived_N);
    Coupling_and_Systematics_merged(merged_prompt_eff, merged_long_lived_eff, merged_prompt_N, merged_long_lived_N);

    std::cout << "\n\n\n";
    std::cout << "Resolved Prompt\n===============\n";
    for (auto& i: resolved_prompt_eff)
    {
        std::cout << i.first << " GeV: " << i.second << '\n';
    }
    std::cout << "\n\n\n";
    std::cout << "Resolved Long Lived\n===================\n";
    for (auto& i: resolved_long_lived_eff)
    {
        std::cout << "\nCoupling = " << i.first << "\n======\n";
        for (auto& j: i.second)
        {
            std::cout << j.first << " GeV: " << j.second << '\n';
        }
    }
    std::cout << "\n\n\n";
    std::cout << "Merged Prompt\n=============\n";
    for (auto& i: merged_prompt_eff)
    {
        std::cout << i.first << " GeV: " << i.second << '\n';
    }
    std::cout << "\n\n\n";
    std::cout << "Merged Long Lived\n=================\n";
    for (auto& i: merged_long_lived_eff)
    {
        std::cout << "\nCoupling = " << i.first << "\n======\n";
        for (auto& j: i.second)
        {
            std::cout << j.first << " GeV: " << j.second << '\n';
        }
    }
    std::cout << "\n\n\n";

    std::map<float, float> BR_prompt_merged = {{1.0f, 0.0324852f}, {2.0f, 0.02782953f}};
    std::map<float, float> BR_prompt_resolved = {{2.0f, 0.0198092f}, {3.0f, 0.00362767f}, {5.0f, 0.001756f}};

    std::map<float, std::map<float, float>> BR_ll_merged, BR_ll_resolved, BR_ll_merged_unc, BR_ll_resolved_unc;
    std::vector<float> ALP_photon_couplings = {1.0f, 0.01f, 0.001f};

    ///Now, time to populate the branching ratios for the long-lived case
    for (auto& coupling: ALP_photon_couplings)
    {
        //looping over the prompot branching ratios for the merged category (<= 2 GeV);
        for (auto& i: BR_prompt_merged)
        {
            // check if efficiencies for prompt and displaced merged for the mass point i.first were calculated
            if (merged_prompt_eff.count(i.first) && merged_long_lived_eff[coupling].count(i.first))
            {
                // BR_ll = BR_prompt * (eff_prompt / eff_long_lived)
                BR_ll_merged[coupling][i.first] = i.second * (merged_prompt_eff[i.first] / merged_long_lived_eff[coupling][i.first]);
                BR_ll_merged_unc[coupling][i.first] = unc(merged_prompt_eff[i.first], merged_long_lived_eff[coupling][i.first], merged_prompt_N[i.first], merged_long_lived_N[coupling][i.first]);

            }
        }

        //looping over the prompt branching ratios for the resolved category (>= 2 GeV)
        for (auto& i: BR_prompt_resolved)
        {
            // check if efficiencies for prompt and displaced resolved for the mass point i.first were calculated
            if (resolved_prompt_eff.count(i.first) && resolved_long_lived_eff[coupling].count(i.first))
            {
                // BR_ll = BR_prompt * (eff_prompt / eff_long_lived)
                BR_ll_resolved[coupling][i.first] = i.second * (resolved_prompt_eff[i.first] / resolved_long_lived_eff[coupling][i.first]);
                BR_ll_resolved_unc[coupling][i.first] = unc(resolved_prompt_eff[i.first], resolved_long_lived_eff[coupling][i.first], resolved_prompt_N[i.first], resolved_long_lived_N[coupling][i.first]);
            }
        }
    }

    std::cout << "Branching Ratios Long Lived\n=================\n";

    std::cout << "Merged\n======\n";
    for (auto& i: BR_ll_merged)
    {
        std::cout << "\nCoupling = " << i.first << "\n=======\n";
        for (auto& j: i.second)
        {
            std::cout << j.first << " GeV: " << j.second << '\n';
        }
    }
    std::cout << "\n\n\n";
    std::cout << "Resolved\n========\n";
    for (auto& i: BR_ll_resolved)
    {
        std::cout << "\nCoupling = " << i.first << "\n=======\n";
        for (auto& j: i.second)
        {
            std::cout << j.first << " GeV: " << j.second << '\n';
        }
    }
    std::cout << "\n\n\n";

    std::cout << "Branching Ratios Prompt\n=================\n";
    for (auto& i: BR_prompt_merged)
    {
        std::cout << i.first << " GeV: " << i.second << '\n';
    }
    for (auto& i: BR_prompt_resolved)
    {
        std::cout << i.first << " GeV: " << i.second << '\n';
    }

    std::cout << "\n\n\n";

    // Create a TCanvas to hold the plot
    TCanvas* canvas = new TCanvas("canvas", "Lines Plot", 800, 600);
    canvas->SetLogy();
    // Create an empty multi-graph to hold all the lines
    TMultiGraph* mg = new TMultiGraph();
    // Create a legend for the limit plot
    TLegend *legend = new TLegend(0.6, 0.6, 0.9, 0.9);
    // no fill color
    legend->SetFillColor(0);
    legend->SetTextSize(0.03); // Adjust the font size here
    //colors for the ALP_photon_couplings
    std::vector<EColor> colors = {kRed, kBlue, kGreen};
    //index for colors and ALP_photon_couplings
    int index = 0;
    // Iterate over the outer unordered_map BR_ll_merged to get the couplings and inner unordered_maps (masses and branching ratios)
    for (const auto& Coupling : BR_ll_merged)
    {
        const float coupling = Coupling.first;
        const std::map<float, float>& innerMap = Coupling.second;

        // Create arrays to hold the x and y values for the current line
        const size_t nPoints = innerMap.size();
        float* xValues = new float[nPoints];
        float* yValues = new float[nPoints];
        float* yErrValues = new float[nPoints];

        // Iterate over the inner unordered_map to get the x and y values
        size_t i = 0;
        for (const auto& innerPair : innerMap)
        {
            if (not std::isinf(innerPair.second))
            {
                xValues[i] = innerPair.first; //mass
                yValues[i] = innerPair.second; //branching ratio
                yErrValues[i] = BR_ll_merged_unc[coupling][innerPair.first]; //uncertainty on branching ratio
                if (std::isinf(BR_ll_merged_unc[coupling][innerPair.first]))
                {
                    yErrValues[i] = 0;
                }
                std::cout << coupling << ' ' << xValues[i] << ' ' << yValues[i]  << ' ' << yErrValues[i] << '\n';
                ++i;
            }
            else
            {
                std::cout << "inf!\n";
            }
        }

        // Create a TGraph for the current line
        TGraphErrors* graph = new TGraphErrors(i, xValues, yValues, 0, yErrValues);

        // Set line color, style, and width (adjust as desired)
        graph->SetLineColor(colors[index]); //red, blue, green
        graph->SetMarkerColor(colors[index]); //red, blue, green
        graph->SetMarkerStyle(20+index); // filled circle, square, diamond
        index++;
        graph->SetLineWidth(2);

        // Add the graph to the multi-graph
        mg->Add(graph);

        // Clean up dynamic memory
        delete[] xValues;
        delete[] yValues;
    }
    index = 0;
    // Iterate over the outer unordered_map BR_ll_resolved to get the couplings and inner unordered_maps (masses and branching ratios, and uncertainties)
    for (const auto& Coupling : BR_ll_resolved)
    {
        const float coupling = Coupling.first;
        const std::map<float, float>& innerMap = Coupling.second;

        // Create arrays to hold the x and y values for the current line
        const size_t nPoints = innerMap.size();
        float* xValues = new float[nPoints];
        float* yValues = new float[nPoints];
        float* yErrValues = new float[nPoints];

        // Iterate over the inner unordered_map to get the x and y values
        size_t i = 0;
        for (const auto& innerPair : innerMap)
        {
            if (not std::isinf(innerPair.second))
            {
                xValues[i] = innerPair.first;
                yValues[i] = innerPair.second;
                yErrValues[i] = BR_ll_resolved_unc[coupling][innerPair.first];
                std::cout << coupling << ' ' << xValues[i] << ' ' << yValues[i]  << ' ' << yErrValues[i] << '\n';
                ++i;
            }
            else
            {
                std::cout << "inf!\n";
            }
        }
        std::cout << '\n';

        // Create a TGraph for the current line
        TGraphErrors* graph = new TGraphErrors(i, xValues, yValues, 0, yErrValues);

        // Set line color, style, and width (adjust as desired)
        graph->SetLineColor(colors[index]); //red, blue, green
        graph->SetMarkerColor(colors[index]); //red, blue, green
        graph->SetMarkerStyle(20+index); // filled circle, square, diamond
        index++;
        graph->SetLineWidth(2);

        // Add the graph to the multi-graph
        mg->Add(graph);
        legend->AddEntry(graph, Form("Limit for C_{a#gamma#gamma} = %.3f", coupling), "p");

        // Clean up dynamic memory
        delete[] xValues;
        delete[] yValues;
    }

    //now plot the prompt case from the paper draft
    float massValues[] = {1.0f, 2.0f, 2.0f, 3.0f, 5.0f};
    float promptBR[] = {0.0324852f, 0.02782953f, 0.0198092f, 0.00362767f, 0.001756f};

    TGraph* graph = new TGraph(5, massValues, promptBR);
    // Set line color, style, and width (adjust as desired)
    graph->SetLineColor(kBlack); //red, blue, green
    graph->SetMarkerColor(kBlack); //red, blue, green
    graph->SetMarkerStyle(47); // filled cross
    graph->SetLineWidth(2);

    // Add the graph to the multi-graph
    mg->Add(graph);
    legend->AddEntry(graph, "Paper Draft Limit", "p");

    mg->SetTitle(";m_{a}  [GeV];Br(H#rightarrow Za) #times Br(a#rightarrow#gamma#gamma)");
    mg->GetYaxis()->SetTitleOffset(1.4); //By default, the title offset is 1.0
    mg->GetXaxis()->SetTitleOffset(1.1); //By default, the title offset is 1.0
    // Draw the multi-graph
    mg->Draw("ALP");  // "ALP" means draw lines with points

    legend->Draw();
    // Update the canvas
    canvas->Update();
    // Save the plot
    canvas->SaveAs("Limits.pdf");


}

void CutFlow()
{
    auto start_time = Clock::now();

//    Table3();
//    Table8();
//    Table11();
//    Table3_Displaced_Axions();
//    Table8_Displaced_Axions();
//    Table11_Displaced_Axions();
//    Coupling();
//    Coupling_and_Systematics_resolved();
//    Coupling_and_Systematics_merged();
    LimitPlot();

    auto end_time = Clock::now();
    std::cout << "Time difference: "
       << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count()/1e9 << " seconds" << std::endl;
}


int main()
{
    CutFlow();
}


//{photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up, photon_id_eff:PH_EFF_ID_Uncertainty__1up, photon_iso_eff:PH_EFF_ISO_Uncertainty__1up, photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down, photon_id_eff:PH_EFF_ID_Uncertainty__1down, photon_iso_eff:PH_EFF_ISO_Uncertainty__1down}
//
//{photons_and_electrons:EG_SCALE_ALL__1up}
//{photons_and_electrons:EG_RESOLUTION_ALL__1up}
//{photons_and_electrons:EG_RESOLUTION_ALL__1down}
//{photons_and_electrons:EG_SCALE_ALL__1down}


