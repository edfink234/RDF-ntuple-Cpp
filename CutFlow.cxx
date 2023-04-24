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

float roundToOneDecimalPlace(float num) {
    std::stringstream stream;
    stream << std::fixed << std::setprecision(1) << num;
    float rounded_num = std::stof(stream.str());
    return rounded_num;
}


void Table3()
{
    std::vector<std::vector<std::string>> input_filenames = {
        //Z gamma background
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617070._000001.LGNTuple.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617064._000001.LGNTuple.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617074._000001.LGNTuple.root"},
        //Signal
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root"},
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
////        "PH_EFF_ID_Uncertainty__1up",
////        "PH_EFF_TRIGGER_Uncertainty__1up",
//        "EG_SCALE_ALL__1down",
//        "EG_SCALE_ALL__1up",
////        "PH_EFF_ID_Uncertainty__1down",
////        "PH_EFF_TRIGGER_Uncertainty__1down",
////        "PH_EFF_ISO_Uncertainty__1up",
//        "EG_RESOLUTION_ALL__1up",
//        "EG_RESOLUTION_ALL__1down",
//    };
    
    std::array<double,9> JetNumeratorSFs = {((139e15)*(1.9828e-9)*(0.821204)),((139e15)*(110.64e-12)*(0.69275)),((139e15)*(40.645e-12)*(0.615906)),((139e15)*(1.9817e-9)*(0.1136684)),((139e15)*(110.47e-12)*(0.1912956)),((139e15)*(40.674e-12)*(0.2326772)),((139e15)*(1.9819e-9)*(0.0656969)),((139e15)*(110.53e-12)*(0.1158741)),((139e15)*(40.68e-12)*(0.1535215))}; //numerators for jet bkg
    std::array<double,3> SFs = {((139e15)*(.871e-12)),((139e15)*(.199e-12)), ((139e15)*(.0345e-15))}; //numerators for Z-gamma bkg
    
    constexpr std::array<const char*,15> Samples = {R"--(pty2\_9\_17)--", R"--(pty\_17\_myy\_0\_80)--", R"--(pty\_17\_myy\_80)--", R"--(Signal $m_{\text{A}}$ = 1 GeV)--", R"--(Signal $m_{\text{A}}$ = 5 GeV)--", "Data", R"--(Zee\_lightJet\_0-70)--", R"--(Zee\_lightJet\_70-140)--", R"--(Zee\_lightJet\_140-280)--", R"--(Zee\_cJet\_0-70)--", R"--(Zee\_cJet\_70-140)--", R"--(Zee\_cJet\_140-280)--", R"--(Zee\_bJet\_0-70)--", R"--(Zee\_bJet\_70-140)--", R"--(Zee\_bJet\_140-280)--"};
    
    std::ostringstream os;
    os << R"--(\section*{Table 3})--" << '\n';
    os << R"--(\hspace{-3cm}\scalebox{0.45}{)--" << '\n';
    os << R"--(\begin{tabular}{|c|c|c|c|c|c|c|c|})--" << '\n';
    os << R"--(\hline)--" << '\n';
    os << R"--(Sample & Before Preselection & 2 leptons
          & Opposite Charge & $p_{T}^{\text{leading}} > 27$ GeV, \; $p_{T}^{\text{sub-leading}} > 20$ GeV & Same flavour & dilep mass cut & dilep $p_{T}$ cut \\ \hline )--" << '\n';
    
    double beforePreselecZGamma = 0, twoLeptonsZGamma = 0, oppChargeZGamma = 0, leadingPtZGamma = 0, deltaRZGamma = 0, MassZGamma = 0, ptCutZGamma = 0;
    double beforePreselecZGammaStatUnc = 0, twoLeptonsZGammaStatUnc = 0, oppChargeZGammaStatUnc = 0, leadingPtZGammaStatUnc = 0, deltaRZGammaStatUnc = 0, MassZGammaStatUnc = 0, ptCutZGammaStatUnc = 0;

    double beforePreselecZJets = 0, twoLeptonsZJets = 0, oppChargeZJets = 0, leadingPtZJets = 0, deltaRZJets = 0, MassZJets = 0, ptCutZJets = 0;
    double beforePreselecZJetsStatUnc = 0, twoLeptonsZJetsStatUnc = 0, oppChargeZJetsStatUnc = 0, leadingPtZJetsStatUnc = 0, deltaRZJetsStatUnc = 0, MassZJetsStatUnc = 0, ptCutZJetsStatUnc = 0;
    
    std::vector<ROOT::RDF::RResultHandle> Nodes;
//    std::vector<ROOT::RDF::RResultHandle> tempNodes;
//
//    std::vector<RResultMap<ULong64_t>> resultmaps;
    
    for (auto& file: input_filenames)
    {
        SchottDataFrame df(MakeRDF(file, 8));
        
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
        
        auto two_leptons = df.Define("di_electrons",
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

        },{"electrons"}).Filter([](RVec<Muon>& muons, RVec<Electron> di_electrons)
        {
            return (di_electrons.size()==2 && muons.empty() && DeltaR(di_electrons[0].Vector(), di_electrons[1].Vector()) > 0.01);
            
        }, {"muons", "di_electrons"});
        
        auto opp_charge = two_leptons.Filter([](RVec<Electron> di_electrons)
        {
            return (di_electrons[0].electron_charge*di_electrons[1].electron_charge < 0);
            
        }, {"di_electrons"});
        
        auto leading_pt = opp_charge.Filter([](RVec<Electron>& electrons)
        {
            return ((electrons[0].electron_pt >= 20e3 && electrons[1].electron_pt >= 27e3) || (electrons[1].electron_pt >= 20e3 && electrons[0].electron_pt >= 27e3));
        }, {"di_electrons"});
        
//        auto delta_R = leading_pt.Filter([] (RVec<Electron>& electrons)
//        {
//            return (DeltaR(electrons[0].Vector(), electrons[1].Vector()) > 0.01);
//        }, {"di_electrons"});
        
        auto same_flavour = leading_pt.Filter([] (RVec<Electron>& electrons)
        {
            return true; //abs(electrons[0].electron_pdg_id) == abs(electrons[1].electron_pdg_id) == 11;
        }, {"di_electrons"});
        
        auto dilep = same_flavour.Define("dilep",[] (RVec<Electron>& electrons)
        {
            return (electrons[0].Vector() + electrons[1].Vector());
        }, {"di_electrons"});
        
        auto mass = dilep.Filter([] (PtEtaPhiEVector& dilep)
        {
            auto mass = dilep.M()/1e3;
            return ((mass >= 81) && (mass <= 101));
        }, {"dilep"});
        
        auto pt_cut = mass.Filter([] (PtEtaPhiEVector& dilep)
        {
            auto pT = dilep.Pt()/1e3;
            return pT > 10;
        }, {"dilep"});
        
        Nodes.push_back(df.Count());
        Nodes.push_back(two_leptons.Count());
        Nodes.push_back(opp_charge.Count());
        Nodes.push_back(leading_pt.Count());
        Nodes.push_back(same_flavour.Count());
        Nodes.push_back(mass.Count());
        Nodes.push_back(pt_cut.Count());
        
//        resultmaps.push_back(VariationsFor(df.Count()));
//        resultmaps.push_back(VariationsFor(two_leptons.Count()));
//        resultmaps.push_back(VariationsFor(opp_charge.Count()));
//        resultmaps.push_back(VariationsFor(leading_pt.Count()));
//        resultmaps.push_back(VariationsFor(same_flavour.Count()));
//        resultmaps.push_back(VariationsFor(mass.Count()));
//        resultmaps.push_back(VariationsFor(pt_cut.Count()));
//
//        tempNodes.push_back(df.Take<RVec<std::vector<std::string>>, RVec<RVec<std::vector<std::string>>>>("electron_syst_name"));
//        tempNodes.push_back(df.Take<RVec<std::vector<std::string>>, RVec<RVec<std::vector<std::string>>>>("photon_syst_name"));
    }

    constexpr std::array<const char*,7> Cuts = {"total", "two leptons", "opposite charge", "leading pt", "same flavour", "mass", "pt cut"};
    
    ROOT::RDF::RunGraphs(Nodes); // running all computation nodes concurrently
    
//     0      1    2     3     4     5     6      Z-gamma
//     7      8    9     10    11    12    13     Z-gamma
//     14    15    16    17    18    19    20     Z-gamma
//     21    22    23    24    25    26    27     ma1
//     28    29    30    31    32    33    34     ma5
//     35    36    37    38    39    40    41     data
//     42    43    44    45    46    47    48     Z-jets
//     49    50    51    52    53    54    55     Z-jets
//     56    57    58    59    60    61    62     Z-jets
//     63    64    65    66    67    68    69     Z-jets
//     70    71    72    73    74    75    76     Z-jets
//     77    78    79    80    81    82    83     Z-jets
//     84    85    86    87    88    89    90     Z-jets
//     91    92    93    94    95    96    97     Z-jets
//     98    99    100   101   102   103   104    Z-jets
    
//    std::unordered_set<std::string> uniqueSystematics;
//    std::unordered_set<std::string> ZGammaSystematics;
//    std::unordered_set<std::string> SignalSystematics;
//    std::unordered_set<std::string> DataSystematics;
//    std::unordered_set<std::string> ZJetsSystematics;
//
//    for (auto& i: tempNodes)
//    {
//        for (auto& j: *i.GetResultPtr<RVec<RVec<std::vector<std::string>>>>())
//        {
//            for (auto& k: j)
//            {
//                for (auto& l: k)
//                {
//                    uniqueSystematics.insert(l);
//                }
//            }
//        }
//    }
//
//    int counter = 0;
//    for (auto& i: tempNodes)
//    {
//        for (auto& j: (*i.GetResultPtr<RVec<RVec<std::vector<std::string>>>>())[0])
//        {
//            for (auto& k: j)
//            {
//                if (counter >= 0 && counter <= 5) //Z-gamma
//                {
//                    ZGammaSystematics.insert(k);
//                }
//
//                else if (counter >= 6 && counter <= 9) //Signal
//                {
//                    SignalSystematics.insert(k);
//                }
//
//                else if (counter == 10) //Data
//                {
//                    DataSystematics.insert(k);
//                }
//
//                else
//                {
//                    ZJetsSystematics.insert(k);
//                }
//            }
//        }
//        counter++;
//    }
////
//    std::cout << "ZGammaSystematics\n=================\n";
//    for (auto& i: ZGammaSystematics)
//    {
//        std::cout << i << '\n';
//    }
//    std::cout << "\n\n\n\n";
//
//    std::cout << "SignalSystematics\n=================\n";
//    for (auto& i: SignalSystematics)
//    {
//        std::cout << i << '\n';
//    }
//    std::cout << "\n\n\n\n";
//
//    std::cout << "DataSystematics\n===============\n";
//    for (auto& i: DataSystematics)
//    {
//        std::cout << i << '\n';
//    }
//    std::cout << "\n\n\n\n";
//
//    std::cout << "ZJetsSystematics\n================\n";
//    for (auto& i: ZJetsSystematics)
//    {
//        std::cout << i << '\n';
//    }
//    std::cout << "\n\n\n\n";
    
//
//    std::cout << "\n\n";
//    std::cout << resultmaps.size() << '\n';
//    for (auto i = 0; i < resultmaps.size(); i++)
//    {
//        if (i % 7 == 0)
//        {
//            std::cout << Samples[i/7] << "\n============================\n\n";
//        }
//        std::cout << Cuts[i%7] << "\n===============\n";
//        for (auto& var: resultmaps[i].GetKeys())
//        {
//            std::cout << std::setw(44) << var <<
//            std::setw(44) << resultmaps[i][var] << '\n';
//        }
//        std::cout << '\n';
//    }
    
    std::cout << R"--(\section*{Table 3 Signal Ratios})--" << '\n';
    std::cout << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
    std::cout << R"--(\begin{tabular}{|c|c|c|c|c|c|c|})--" << '\n';
    std::cout << R"--(\hline)--" << '\n';
    std::cout << R"--(Sample & $\frac{\text{2 leptons}}{\text{2 leptons}}$
          & $\frac{\text{Opposite Charge}}{\text{2 leptons}}$ & $\frac{p_{T}^{\text{leading}} > 27\text{ GeV, } \, p_{T}^{\text{sub-leading}} > 20 \text{ GeV }}{\text{2 leptons}}$ & $\frac{\text{Same flavour}}{\text{2 leptons}}$ & $\frac{\text{dilep mass cut}}{\text{2 leptons}}$ & $\frac{\text{dilep }p_{T} \text{ cut}}{\text{2 leptons}}$ \\ \hline )--" << '\n';
    
    os.setf(std::ios::fixed);
    os.precision(2);
    
    for (int i=0, j=0; (i<15 && j <= 98); i++, j+=7)
    {
        os << Samples[i];
        if (i >= 0 && i <= 2) //Z-gamma
        {
            os  << " & " << *Nodes[j].GetResultPtr<ULong64_t>() * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j].GetResultPtr<ULong64_t>()) * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+1].GetResultPtr<ULong64_t>() * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+1].GetResultPtr<ULong64_t>()) * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+2].GetResultPtr<ULong64_t>() * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+2].GetResultPtr<ULong64_t>()) * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+3].GetResultPtr<ULong64_t>() * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+3].GetResultPtr<ULong64_t>()) * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+4].GetResultPtr<ULong64_t>() * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+4].GetResultPtr<ULong64_t>()) * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+5].GetResultPtr<ULong64_t>() * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+5].GetResultPtr<ULong64_t>()) * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+6].GetResultPtr<ULong64_t>() * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+6].GetResultPtr<ULong64_t>()) * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--( \\ \hline )--" << '\n';
        }
        
        else if (i >= 6) //Z-jets
        {
            os  << " & " << *Nodes[j].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-6] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i-6] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+1].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-6] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+1].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i-6] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+2].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-6] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+2].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i-6] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+3].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-6] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+3].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i-6] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+4].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-6] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+4].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i-6] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+5].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-6] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+5].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i-6] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+6].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-6] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+6].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i-6] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--( \\ \hline )--" << '\n';
        }
        else
        {
            os  << " & " << *Nodes[j].GetResultPtr<ULong64_t>()
            << R"--($\, \pm \,$)--" << sqrt(*Nodes[j].GetResultPtr<ULong64_t>())
            << " & " << *Nodes[j+1].GetResultPtr<ULong64_t>()
            << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+1].GetResultPtr<ULong64_t>())
            << " & " << *Nodes[j+2].GetResultPtr<ULong64_t>()
            << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+2].GetResultPtr<ULong64_t>())
            << " & " << *Nodes[j+3].GetResultPtr<ULong64_t>()
            << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+3].GetResultPtr<ULong64_t>())
            << " & " << *Nodes[j+4].GetResultPtr<ULong64_t>()
            << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+4].GetResultPtr<ULong64_t>())
            << " & " << *Nodes[j+5].GetResultPtr<ULong64_t>()
            << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+5].GetResultPtr<ULong64_t>())
            << " & " << *Nodes[j+6].GetResultPtr<ULong64_t>()
            << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+6].GetResultPtr<ULong64_t>())
            << R"--( \\ \hline )--" << '\n';
            
            if (i==3) //1 GeV
            {
                std::cout << Samples[i] << " (Me) & " << 1 << " & " <<
                static_cast<double>(*Nodes[j+2].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
                << " & " <<
                static_cast<double>(*Nodes[j+3].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
                << " & " <<
                static_cast<double>(*Nodes[j+4].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
                << " & " <<
                static_cast<double>(*Nodes[j+5].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
                << " & " <<
                static_cast<double>(*Nodes[j+6].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
                << R"--( \\ \hline )--" << '\n';
                
                std::cout << Samples[i] << " (Paper) & " << 1 << " & " <<
                21505.03 / 21606.75
                << " & " <<
                21375.48 / 21606.75
                << " & " <<
                21375.33 / 21606.75
                << " & " <<
                20543.06 / 21606.75
                << " & " <<
                19516.87 / 21606.75
                << R"--( \\ \hline )--" << '\n';
            }
            
            if (i==4) //5 GeV
            {
                std::cout << Samples[i] << " (Me) & " << 1 << " & " <<
                static_cast<double>(*Nodes[j+2].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
                << " & " <<
                static_cast<double>(*Nodes[j+3].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
                << " & " <<
                static_cast<double>(*Nodes[j+4].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
                << " & " <<
                static_cast<double>(*Nodes[j+5].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
                << " & " <<
                static_cast<double>(*Nodes[j+6].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
                << R"--( \\ \hline )--" << '\n';
                
                std::cout << Samples[i] << " (Paper) & " << 1 << " & " <<
                21556.03 / 21655.38
                << " & " <<
                21424.29 / 21655.38
                << " & " <<
                21424.28 / 21655.38
                << " & " <<
                20585.09 / 21655.38
                << " & " <<
                19536.68 / 21655.38
                << R"--( \\ \hline )--" << '\n';
            }
        }
    }
    
    std::cout << R"--(\end{tabular}})--" << "\n\n\n";
    
    for (int i = 0, j = 0; (i <= 14 && j <= 2); i += 7, j++)
    {
        beforePreselecZGamma += *Nodes[i].GetResultPtr<ULong64_t>() * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
        twoLeptonsZGamma += *Nodes[i+1].GetResultPtr<ULong64_t>() * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
        oppChargeZGamma += *Nodes[i+2].GetResultPtr<ULong64_t>() * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
        leadingPtZGamma += *Nodes[i+3].GetResultPtr<ULong64_t>() * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
        deltaRZGamma += *Nodes[i+4].GetResultPtr<ULong64_t>() * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
        MassZGamma += *Nodes[i+5].GetResultPtr<ULong64_t>() * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
        ptCutZGamma += *Nodes[i+6].GetResultPtr<ULong64_t>() * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
        
        beforePreselecZGammaStatUnc += pow(sqrt(*Nodes[i].GetResultPtr<ULong64_t>()) * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
        twoLeptonsZGammaStatUnc += pow(sqrt(*Nodes[i+1].GetResultPtr<ULong64_t>()) * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
        oppChargeZGammaStatUnc += pow(sqrt(*Nodes[i+2].GetResultPtr<ULong64_t>()) * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
        leadingPtZGammaStatUnc += pow(sqrt(*Nodes[i+3].GetResultPtr<ULong64_t>()) * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
        deltaRZGammaStatUnc += pow(sqrt(*Nodes[i+4].GetResultPtr<ULong64_t>()) * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
        MassZGammaStatUnc += pow(sqrt(*Nodes[i+5].GetResultPtr<ULong64_t>()) * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
        ptCutZGammaStatUnc += pow(sqrt(*Nodes[i+6].GetResultPtr<ULong64_t>()) * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
    }
    
    for (int i = 42, j = 0; (i <= 98 && j <= 8); i += 7, j++)
    {
        beforePreselecZJets += *Nodes[i].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
        twoLeptonsZJets += *Nodes[i+1].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
        oppChargeZJets += *Nodes[i+2].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
        leadingPtZJets += *Nodes[i+3].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
        deltaRZJets += *Nodes[i+4].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
        MassZJets += *Nodes[i+5].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
        ptCutZJets += *Nodes[i+6].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
        
        beforePreselecZJetsStatUnc += pow(sqrt(*Nodes[i].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
        twoLeptonsZJetsStatUnc += pow(sqrt(*Nodes[i+1].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
        oppChargeZJetsStatUnc += pow(sqrt(*Nodes[i+2].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
        leadingPtZJetsStatUnc += pow(sqrt(*Nodes[i+3].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
        deltaRZJetsStatUnc += pow(sqrt(*Nodes[i+4].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
        MassZJetsStatUnc += pow(sqrt(*Nodes[i+5].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
        ptCutZJetsStatUnc += pow(sqrt(*Nodes[i+6].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
    }

    os << R"--(Total $Z\gamma$ & )--"
    << beforePreselecZGamma << R"--($\, \pm \,$)--" << sqrt(beforePreselecZGammaStatUnc) << " & "
    << twoLeptonsZGamma << R"--($\, \pm \,$)--" << sqrt(twoLeptonsZGammaStatUnc) << " & "
    << oppChargeZGamma << R"--($\, \pm \,$)--" << sqrt(oppChargeZGammaStatUnc) << " & "
    << leadingPtZGamma << R"--($\, \pm \,$)--" << sqrt(leadingPtZGammaStatUnc) << " & "
    << deltaRZGamma << R"--($\, \pm \,$)--" << sqrt(deltaRZGammaStatUnc) << " & "
    << MassZGamma << R"--($\, \pm \,$)--" << sqrt(MassZGammaStatUnc) << " & "
    << ptCutZGamma << R"--($\, \pm \,$)--" << sqrt(ptCutZGammaStatUnc)
    << R"--( \\ \hline )--" << '\n';
    
    os << R"--(Total $Z$ jets & )--"
    << beforePreselecZJets << R"--($\, \pm \,$)--" << sqrt(beforePreselecZJetsStatUnc) << " & "
    << twoLeptonsZJets << R"--($\, \pm \,$)--" << sqrt(twoLeptonsZJetsStatUnc) << " & "
    << oppChargeZJets << R"--($\, \pm \,$)--" << sqrt(oppChargeZJetsStatUnc) << " & "
    << leadingPtZJets << R"--($\, \pm \,$)--" << sqrt(leadingPtZJetsStatUnc) << " & "
    << deltaRZJets << R"--($\, \pm \,$)--" << sqrt(deltaRZJetsStatUnc) << " & "
    << MassZJets << R"--($\, \pm \,$)--" << sqrt(MassZJetsStatUnc) << " & "
    << ptCutZJets << R"--($\, \pm \,$)--" << sqrt(ptCutZJetsStatUnc)
    << R"--( \\ \hline )--" << '\n';
    
    os << R"--(Total Bkg & )--"
    << beforePreselecZJets+beforePreselecZGamma << R"--($\, \pm \,$)--" <<
    sqrt(beforePreselecZJetsStatUnc+beforePreselecZGammaStatUnc) << " & "
    << twoLeptonsZJets+twoLeptonsZGamma << R"--($\, \pm \,$)--" <<
    sqrt(twoLeptonsZJetsStatUnc+twoLeptonsZGammaStatUnc) << " & "
    << oppChargeZJets+oppChargeZGamma << R"--($\, \pm \,$)--" <<
    sqrt(oppChargeZJetsStatUnc+oppChargeZGammaStatUnc) << " & "
    << leadingPtZJets+leadingPtZGamma << R"--($\, \pm \,$)--" <<
    sqrt(leadingPtZJetsStatUnc+leadingPtZGammaStatUnc) << " & "
    << deltaRZJets+deltaRZGamma << R"--($\, \pm \,$)--" <<
    sqrt(deltaRZJetsStatUnc+deltaRZGammaStatUnc) << " & "
    << MassZJets+MassZGamma << R"--($\, \pm \,$)--" <<
    sqrt(MassZJetsStatUnc+MassZGammaStatUnc) << " & "
    << ptCutZJets+ptCutZGamma << R"--($\, \pm \,$)--" <<
    sqrt(ptCutZJetsStatUnc+ptCutZGammaStatUnc) << R"--( \\ \hline )--" << '\n';
    
    os << R"--(\end{tabular}})--" << '\n';
    std::cout << os.str() << '\n';
}

void Table8()
{
    std::vector<std::vector<std::string>> input_filenames = {
        //Z gamma background
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617070._000001.LGNTuple.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617064._000001.LGNTuple.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617074._000001.LGNTuple.root"},
        //Signal
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root"},
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
    
    std::array<double,9> JetNumeratorSFs = {((139e15)*(1.9828e-9)*(0.821204)),((139e15)*(110.64e-12)*(0.69275)),((139e15)*(40.645e-12)*(0.615906)),((139e15)*(1.9817e-9)*(0.1136684)),((139e15)*(110.47e-12)*(0.1912956)),((139e15)*(40.674e-12)*(0.2326772)),((139e15)*(1.9819e-9)*(0.0656969)),((139e15)*(110.53e-12)*(0.1158741)),((139e15)*(40.68e-12)*(0.1535215))}; //numerators for jet bkg
    std::array<double,3> SFs = {((139e15)*(.871e-12)),((139e15)*(.199e-12)), ((139e15)*(.0345e-15))}; //numerators for Z-gamma bkg
    
    constexpr std::array<const char*,15> Samples = {R"--(pty2\_9\_17)--", R"--(pty\_17\_myy\_0\_80)--", R"--(pty\_17\_myy\_80)--", R"--(Signal $m_{\text{A}}$ = 1 GeV)--", R"--(Signal $m_{\text{A}}$ = 5 GeV)--", "Data", R"--(Zee\_lightJet\_0-70)--", R"--(Zee\_lightJet\_70-140)--", R"--(Zee\_lightJet\_140-280)--", R"--(Zee\_cJet\_0-70)--", R"--(Zee\_cJet\_70-140)--", R"--(Zee\_cJet\_140-280)--", R"--(Zee\_bJet\_0-70)--", R"--(Zee\_bJet\_70-140)--", R"--(Zee\_bJet\_140-280)--"};
    
    double totalEventsZgamma = 0, resolvedEventsZgamma = 0, SREventsZgamma = 0, SBEventsZgamma = 0;
    double totalEventsZgammaStatUnc = 0, resolvedEventsZgammaStatUnc = 0, SREventsZgammaStatUnc = 0, SBEventsZgammaStatUnc = 0;
    double totalEventsZJets = 0, resolvedEventsZJets = 0, SREventsZJets = 0, SBEventsZJets = 0;
    double totalEventsZJetsStatUnc = 0, resolvedEventsZJetsStatUnc = 0, SREventsZJetsStatUnc = 0, SBEventsZJetsStatUnc = 0;
    
    std::vector<ROOT::RDF::RResultHandle> Nodes;
    
    std::ostringstream os;
    os << R"--(\section*{Table 8})--" << '\n';
    os << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
    os << R"--(\begin{tabular}{|c|c|c|c|c|})--" << '\n';
    os << R"--(\hline)--" << '\n';
    os << R"--(Sample & Total Events & PS: Resolved & SB & SR
           \\ \hline )--" << '\n';
    for (auto& file: input_filenames)
    {
        SchottDataFrame df(MakeRDF(file, 8));
        
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
        
        auto resolved = ptCut.Define("photons_pass_cuts",
        [&](RVec<Photon> photons)
        {
            photons.erase(std::remove_if(photons.begin(),photons.end(),
            [](Photon& x)
            {
                return ((abs(x.photon_eta) >= 2.37) || (abs(x.photon_eta) > 1.37 && abs(x.photon_eta) < 1.52) || (!x.photon_id_loose));

            }), photons.end());
            return photons;
        }, {"photons"}).Filter(
        [&](RVec<Photon>& reco_photons_matched)
        {
            if (reco_photons_matched.size() < 2)
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
                return true;
            }
            return false;

        }, {"photons_pass_cuts"});

        auto SR = resolved.Filter(
        [&](RVec<Electron>& electrons, RVec<Photon>& reco_photons_matched)
        {
            auto mass = ((electrons[0].Vector()+electrons[1].Vector()).M() + (reco_photons_matched[0].Vector()+reco_photons_matched[1].Vector()).M())/1e3;
            return ((mass > 110) && (mass < 140));
        }, {"di_electrons", "photons_pass_cuts"});
        
        auto SB = resolved.Filter(
        [&](RVec<Electron>& electrons, RVec<Photon>& reco_photons_matched)
        {
            auto mass = ((electrons[0].Vector()+electrons[1].Vector()).M() + (reco_photons_matched[0].Vector()+reco_photons_matched[1].Vector()).M())/1e3;
            return (!((mass > 110) && (mass < 140)));
        }, {"di_electrons", "photons_pass_cuts"});
        
        Nodes.push_back(df.Count());
        Nodes.push_back(resolved.Count());
        Nodes.push_back(SB.Count());
        Nodes.push_back(SR.Count());
    }
    
    ROOT::RDF::RunGraphs(Nodes); // running all computation nodes concurrently
    
//     0    1    2    3   Z-gamma
//
//     4    5    6    7   Z-gamma
//
//     8    9    10   11  Z-gamma
//
//     12   13   14   15  ma1
//
//     16   17   18   19  ma5
//
//     20   21   22   23  data
//
//     24   25   26   27  Z-jets
//
//     28   29   30   31  Z-jets
//
//     32   33   34   35  Z-jets
//
//     36   37   38   39  Z-jets
//
//     40   41   42   43  Z-jets
//
//     44   45   46   47  Z-jets
//
//     48   49   50   51  Z-jets
//
//     52   53   54   55  Z-jets
//
//     56   57   58   59  Z-jets
    os.setf(std::ios::fixed);
    os.precision(2);
 
    for (int i=0, j=0; (i<15 && j <= 56); i++, j+=4)
    {
        os << Samples[i];
        if (i >= 0 && i <= 2)
        {
            os  << " & " << *Nodes[j].GetResultPtr<ULong64_t>() * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j].GetResultPtr<ULong64_t>()) * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+1].GetResultPtr<ULong64_t>() * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+1].GetResultPtr<ULong64_t>()) * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+2].GetResultPtr<ULong64_t>() * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+2].GetResultPtr<ULong64_t>()) * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+3].GetResultPtr<ULong64_t>() * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+3].GetResultPtr<ULong64_t>()) * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--( \\ \hline )--" << '\n';
        }
        
        else if (i >= 6)
        {
            os  << " & " << *Nodes[j].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-6] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i-6] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+1].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-6] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+1].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i-6] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+2].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-6] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+2].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i-6] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+3].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-6] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+3].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i-6] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--( \\ \hline )--" << '\n';
        }
        
        else
        {
            os  << " & " << *Nodes[j].GetResultPtr<ULong64_t>()
                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+1].GetResultPtr<ULong64_t>()
                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+1].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+2].GetResultPtr<ULong64_t>()
                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+2].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+3].GetResultPtr<ULong64_t>()
                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+3].GetResultPtr<ULong64_t>())
                << R"--( \\ \hline )--" << '\n';
        }
    }
    
    for (int i = 0, j = 0; (i <= 8 && j <= 2); i += 4, j++)
    {
        totalEventsZgamma += *Nodes[i].GetResultPtr<ULong64_t>() * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
        resolvedEventsZgamma += *Nodes[i+1].GetResultPtr<ULong64_t>() * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
        SBEventsZgamma += *Nodes[i+2].GetResultPtr<ULong64_t>() * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
        SREventsZgamma += *Nodes[i+3].GetResultPtr<ULong64_t>() * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
        
        totalEventsZgammaStatUnc += pow(sqrt(*Nodes[i].GetResultPtr<ULong64_t>()) * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
        resolvedEventsZgammaStatUnc += pow(sqrt(*Nodes[i+1].GetResultPtr<ULong64_t>()) * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
        SBEventsZgammaStatUnc += pow(sqrt(*Nodes[i+2].GetResultPtr<ULong64_t>()) * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
        SREventsZgammaStatUnc += pow(sqrt(*Nodes[i+3].GetResultPtr<ULong64_t>()) * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
    }
    
    for (int i = 24, j = 0; (i <= 56 && j <= 8); i += 4, j++)
    {
        totalEventsZJets += *Nodes[i].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
        resolvedEventsZJets += *Nodes[i+1].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
        SBEventsZJets += *Nodes[i+2].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
        SREventsZJets += *Nodes[i+3].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
        
        totalEventsZJetsStatUnc += pow(sqrt(*Nodes[i].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
        resolvedEventsZJetsStatUnc += pow(sqrt(*Nodes[i+1].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
        SBEventsZJetsStatUnc += pow(sqrt(*Nodes[i+2].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
        SREventsZJetsStatUnc += pow(sqrt(*Nodes[i+3].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
    }
    
    os << R"--(Total $Z\gamma$ & )--" << totalEventsZgamma
    << R"--($\, \pm \,$)--" << sqrt(totalEventsZgammaStatUnc)
    << " & " << resolvedEventsZgamma
    << R"--($\, \pm \,$)--" << sqrt(resolvedEventsZgammaStatUnc)
    << " & " << SBEventsZgamma
    << R"--($\, \pm \,$)--" << sqrt(SBEventsZgammaStatUnc)
    << " & " << SREventsZgamma
    << R"--($\, \pm \,$)--" << sqrt(SREventsZgammaStatUnc)
    << R"--( \\ \hline )--" << '\n';
    
    os << R"--(Total $Z$ jets & )--" << totalEventsZJets
    << R"--($\, \pm \,$)--" << sqrt(totalEventsZJetsStatUnc)
    << " & " << resolvedEventsZJets
    << R"--($\, \pm \,$)--" << sqrt(resolvedEventsZJetsStatUnc)
    << " & " << SBEventsZJets
    << R"--($\, \pm \,$)--" << sqrt(SBEventsZJetsStatUnc)
    << " & " << SREventsZJets
    << R"--($\, \pm \,$)--" << sqrt(SREventsZJetsStatUnc)
    << R"--( \\ \hline )--" << '\n';
    
    os << R"--(Total Bkg & )--" << totalEventsZJets+totalEventsZgamma
    << R"--($\, \pm \,$)--" << sqrt(totalEventsZJetsStatUnc+totalEventsZgammaStatUnc)
    << " & " << resolvedEventsZJets+resolvedEventsZgamma
    << R"--($\, \pm \,$)--" << sqrt(resolvedEventsZJetsStatUnc+resolvedEventsZgammaStatUnc)
    << " & " << SBEventsZJets+SBEventsZgamma
    << R"--($\, \pm \,$)--" << sqrt(SBEventsZJetsStatUnc+SBEventsZgammaStatUnc)
    << " & " << SREventsZJets+SREventsZgamma
    << R"--($\, \pm \,$)--" << sqrt(SREventsZJetsStatUnc+SREventsZgammaStatUnc)
    << R"--( \\ \hline )--" << '\n';
    
    os << R"--(\end{tabular}})--" << '\n';
    std::cout << os.str() << '\n';
}

void Table11()
{
    std::vector<std::vector<std::string>> input_filenames = {
        //Z gamma background
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617070._000001.LGNTuple.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617064._000001.LGNTuple.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/user.kschmied.31617074._000001.LGNTuple.root"},
        //Signal
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root"},
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/sample.root", "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/sample2.root"},
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
    
    constexpr std::array<const char*,15> Samples = {R"--(pty2\_9\_17)--", R"--(pty\_17\_myy\_0\_80)--", R"--(pty\_17\_myy\_80)--", R"--(Signal $m_{\text{A}}$ = 1 GeV)--", R"--(Signal $m_{\text{A}}$ = 5 GeV)--", "Data", R"--(Zee\_lightJet\_0-70)--", R"--(Zee\_lightJet\_70-140)--", R"--(Zee\_lightJet\_140-280)--", R"--(Zee\_cJet\_0-70)--", R"--(Zee\_cJet\_70-140)--", R"--(Zee\_cJet\_140-280)--", R"--(Zee\_bJet\_0-70)--", R"--(Zee\_bJet\_70-140)--", R"--(Zee\_bJet\_140-280)--"};
    
    std::array<double,9> JetNumeratorSFs = {((139e15)*(1.9828e-9)*(0.821204)),((139e15)*(110.64e-12)*(0.69275)),((139e15)*(40.645e-12)*(0.615906)),((139e15)*(1.9817e-9)*(0.1136684)),((139e15)*(110.47e-12)*(0.1912956)),((139e15)*(40.674e-12)*(0.2326772)),((139e15)*(1.9819e-9)*(0.0656969)),((139e15)*(110.53e-12)*(0.1158741)),((139e15)*(40.68e-12)*(0.1535215))}; //numerators for jet bkg
    std::array<double,3> SFs = {((139e15)*(.871e-12)),((139e15)*(.199e-12)), ((139e15)*(.0345e-15))}; //numerators for Z-gamma bkg
    
    double totalEventsZgamma = 0, passPreselectionZgamma = 0, photonPtDeltaRCountZgamma = 0, xWindowZgamma = 0,
    srCountZgamma = 0, srIDCountZgamma = 0;
    double totalEventsZgammaStatUnc = 0, passPreselectionZgammaStatUnc = 0, photonPtDeltaRCountZgammaStatUnc = 0, xWindowZgammaStatUnc = 0, srCountZgammaStatUnc = 0, srIDCountZgammaStatUnc = 0;
    
    double totalEventsZjets = 0, passPreselectionZjets = 0, photonPtDeltaRCountZjets = 0, xWindowZjets = 0,
        srCountZjets = 0, srIDCountZjets = 0;
    double totalEventsZjetsStatUnc = 0, passPreselectionZjetsStatUnc = 0, photonPtDeltaRCountZjetsStatUnc = 0, xWindowZjetsStatUnc = 0, srCountZjetsStatUnc = 0, srIDCountZjetsStatUnc = 0;
    
    std::vector<ROOT::RDF::RResultHandle> Nodes;

    std::ostringstream os, ss;
    
    os.setf(std::ios::fixed);
    os.precision(2);
    
    os << R"--(\section*{Table 11})--" << '\n';
    os << R"--(\hspace{-3cm}\scalebox{0.55}{)--" << '\n';
    os << R"--(\begin{tabular}{|c|c|c|c|c|c|c|})--" << '\n';
    os << R"--(\hline)--" << '\n';
    os << R"--(Sample & Total Events & pass preselection (PS) & photon $p_T$ + $\Delta R_{\gamma\gamma}$ cut & $X$ window & SR & SR-ID
           \\ \hline )--" << '\n';
    
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

    for (auto& file: input_filenames)
    {
        SchottDataFrame df(MakeRDF(file, 8));
        
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
        
        auto two_leptons = trigger_selection.Define("di_electrons",
        [](RVec<Electron> electrons)
        {
            electrons.erase(std::remove_if(electrons.begin(),electrons.end(),
            [](Electron& ep)
            {
                return ((ep.electron_pt/1e3 <= 20) || (abs(ep.electron_eta) >= 2.37)
                        || ((1.37 < abs(ep.electron_eta)) && (abs(ep.electron_eta) < 1.52))
                        || (ep.electron_id_medium != 1));
                

            }), electrons.end());

            return electrons;

        },{"electrons"}).Filter([](RVec<Muon>& muons, RVec<Electron> di_electrons)
        {
            return (di_electrons.size()==2 && muons.empty() && DeltaR(di_electrons[0].Vector(), di_electrons[1].Vector()) > 0.01);
            
        }, {"muons", "di_electrons"});
        
        auto opp_charge = two_leptons.Filter([](RVec<Electron> di_electrons)
        {
            return (di_electrons[0].electron_charge*di_electrons[1].electron_charge < 0);
            
        }, {"di_electrons"});
        
        auto leading_pt = opp_charge.Filter([](RVec<Electron>& electrons)
        {
            return ((electrons[0].electron_pt >= 20e3 && electrons[1].electron_pt >= 27e3) || (electrons[1].electron_pt >= 20e3 && electrons[0].electron_pt >= 27e3));
        }, {"di_electrons"});
        
//        auto delta_R = leading_pt.Filter([] (RVec<Electron>& electrons)
//        {
//            return (DeltaR(electrons[0].Vector(), electrons[1].Vector()) > 0.01);
//        }, {"di_electrons"});
        
        auto same_flavour = leading_pt.Filter([] (RVec<Electron>& electrons)
        {
            return true; //abs(electrons[0].electron_pdg_id) == abs(electrons[1].electron_pdg_id) == 11;
        }, {"di_electrons"});
        
        auto dilep = same_flavour.Define("dilep",[] (RVec<Electron>& electrons)
        {
            return (electrons[0].Vector() + electrons[1].Vector());
        }, {"di_electrons"});
        
        auto mass = dilep.Filter([] (PtEtaPhiEVector& dilep)
        {
            auto mass = dilep.M()/1e3;
            return ((mass >= 81) && (mass <= 101));
        }, {"dilep"});
        
        auto ptCut = mass.Filter([] (PtEtaPhiEVector& dilep)
        {
            auto pT = dilep.Pt()/1e3;
            return pT > 10;
        }, {"dilep"});
        
        auto photonPtDeltaR = ptCut.Define("photonPtDeltaR",
        [&](RVec<Photon> photons)
        {
            photons.erase(std::remove_if(photons.begin(),photons.end(),
            [](Photon& x)
            {
                return ((abs(x.photon_eta) >= 2.37) || (x.photon_pt <= 10e3) || (abs(x.photon_eta) > 1.37 && abs(x.photon_eta) < 1.52));
                
            }), photons.end());
            return photons;
        }, {"photons"}).Filter(
       [&](RVec<Photon>& reco_photons_matched)
       {
            if (reco_photons_matched.size() < 2)
            {
              return false;
            }
            RVec<Photon> x;
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
                    x = {reco_photons_matched[combs[0][i]], reco_photons_matched[combs[1][i]]};
                }
            }
            return (chosen_delta_r < 1.5 && pt1 > 10e3 && pt2 > 10e3);

       }, {"photonPtDeltaR"});
                
        auto X_window = photonPtDeltaR.Define("chosen_two",
        [](RVec<Photon>& reco_photons_matched)
        {
            RVec<Photon> x;
            if (reco_photons_matched.size() < 2)
            {
                return x;
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
                    x = {reco_photons_matched[combs[0][i]], reco_photons_matched[combs[1][i]]};
                }
            }
            if (pt1 > 10e3 && pt2 > 10e3 && best_X > 0.96 && best_X < 1.2 && chosen_delta_r < 1.5)
            {
                return x;
            }
            x.clear();
            return x;
        }, {"photonPtDeltaR"}).Filter(
        [&](RVec<Photon>& reco_photons_matched)
        {
            return (reco_photons_matched.size()==2);
        }, {"chosen_two"});
                
        auto SR = X_window.Filter(
        [](RVec<Photon>& photons, RVec<Electron>& electrons)
        {
            PtEtaPhiEVector photonVec = photons[0].Vector() + photons[1].Vector();
            PtEtaPhiEVector electronVec = electrons[0].Vector() + electrons[1].Vector();
            auto mass = (photonVec+electronVec).M()/1e3;
            return ((mass >= 110) && (mass <= 140));
        },{"chosen_two","di_electrons"});
                
        auto SR_ID = SR.Filter(
        [](RVec<Photon>& photons)
        {
            return (photons[0].photon_id_loose && photons[1].photon_id_loose);
        },{"chosen_two"});
        
        Nodes.push_back(df.Count());
        Nodes.push_back(ptCut.Count());
        Nodes.push_back(photonPtDeltaR.Count());
        Nodes.push_back(X_window.Count());
        Nodes.push_back(SR.Count());
        Nodes.push_back(SR_ID.Count());
    }
    
    ROOT::RDF::RunGraphs(Nodes); // running all computation nodes concurrently
    
//    0    1    2    3    4    5      Z-gamma
//
//    6    7    8    9    10   11     Z-gamma
//
//    12   13   14   15   16   17     Z-gamma
//
//    18   19   20   21   22   23     ma1
//
//    24   25   26   27   28   29     ma5
//
//    30   31   32   33   34   35     data
//
//    36   37   38   39   40   41     Z-jets
//
//    42   43   44   45   46   47     Z-jets
//
//    48   49   50   51   52   53     Z-jets
//
//    54   55   56   57   58   59     Z-jets
//
//    60   61   62   63   64   65     Z-jets
//
//    66   67   68   69   70   71     Z-jets
//
//    72   73   74   75   76   77     Z-jets
//
//    78   79   80   81   82   83     Z-jets
//
//    84   85   86   87   88   89     Z-jets

    std::cout << R"--(\section*{Table 11 Signal Ratios})--" << '\n';
    std::cout << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
    std::cout << R"--(\begin{tabular}{|c|c|c|c|c|c|c|})--" << '\n';
    std::cout << R"--(\hline)--" << '\n';
    std::cout << R"--(Sample & $\frac{\text{pass preselection (PS)}}{\text{pass preselection (PS)}}$ & $\frac{\text{photon } p_T \text{ + } \Delta R_{\gamma\gamma} \text{ cut }}{\text{pass preselection (PS)}}$ & $\frac{X \text{ window}}{\text{pass preselection (PS)}}$ & $\frac{\text{SR}}{\text{pass preselection (PS)}}$ & $\frac{\text{SR-ID}}{\text{pass preselection (PS)}}$
           \\ \hline )--" << '\n';
    
    ss << R"--(\section*{Table 11 Signal Ratios})--" << '\n';
    ss << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
    ss << R"--(\begin{tabular}{|c|c|c|c|c|c|c|})--" << '\n';
    ss << R"--(\hline)--" << '\n';
    ss << R"--(Sample & $\frac{\text{pass preselection (PS)}}{\text{pass preselection (PS)}}$ & $\frac{\text{photon } p_T \text{ + } \Delta R_{\gamma\gamma} \text{ cut }}{\text{pass preselection (PS)}}$ & $\frac{X \text{ window}}{\text{photon } p_T \text{ + } \Delta R_{\gamma\gamma} \text{ cut }}$ & $\frac{\text{SR}}{X \text{ window}}$ & $\frac{\text{SR-ID}}{\text{SR}}$
           \\ \hline )--" << '\n';
    
    for (int i=0, j=0; (i<15 && j <= 84); i++, j+=6)
    {
        os << Samples[i];
        
        if (i >= 0 && i <= 2)
        {
            os  << " & " << *Nodes[j].GetResultPtr<ULong64_t>() * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--($\, \pm \,$)--" <<  sqrt(*Nodes[j].GetResultPtr<ULong64_t>()) * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+1].GetResultPtr<ULong64_t>() * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--($\, \pm \,$)--" <<  sqrt(*Nodes[j+1].GetResultPtr<ULong64_t>()) * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+2].GetResultPtr<ULong64_t>() * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--($\, \pm \,$)--" <<  sqrt(*Nodes[j+2].GetResultPtr<ULong64_t>()) * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+3].GetResultPtr<ULong64_t>() * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--($\, \pm \,$)--" <<  sqrt(*Nodes[j+3].GetResultPtr<ULong64_t>()) * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+4].GetResultPtr<ULong64_t>() * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--($\, \pm \,$)--" <<  sqrt(*Nodes[j+4].GetResultPtr<ULong64_t>()) * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+5].GetResultPtr<ULong64_t>() * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--($\, \pm \,$)--" <<  sqrt(*Nodes[j+5].GetResultPtr<ULong64_t>()) * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--( \\ \hline )--" << '\n';
        }
        else if (i >= 6)
        {
            os  << " & " << *Nodes[j].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-6] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i-6] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+1].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-6] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+1].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i-6] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+2].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-6] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+2].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i-6] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+3].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-6] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+3].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i-6] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+4].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-6] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+4].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i-6] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+5].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-6] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+5].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i-6] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--( \\ \hline )--" << '\n';
        }
        else
        {
            os  << " & " << *Nodes[j].GetResultPtr<ULong64_t>()
                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+1].GetResultPtr<ULong64_t>()
                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+1].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+2].GetResultPtr<ULong64_t>()
                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+2].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+3].GetResultPtr<ULong64_t>()
                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+3].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+4].GetResultPtr<ULong64_t>()
                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+4].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+5].GetResultPtr<ULong64_t>()
                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+5].GetResultPtr<ULong64_t>())
                << R"--( \\ \hline )--" << '\n';
            
            if (i==3) //1 GeV
            {
                std::cout << Samples[i] << " (Me) & " << 1 << " & " <<
                static_cast<double>(*Nodes[j+2].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
                << " & " <<
                static_cast<double>(*Nodes[j+3].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
                << " & " <<
                static_cast<double>(*Nodes[j+4].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
                << " & " <<
                static_cast<double>(*Nodes[j+5].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
                << R"--( \\ \hline )--" << '\n';
                
                std::cout << Samples[i] << " (Paper) & " << 1 << " & " <<
                2343.65 / 19516.87
                << " & " <<
                2204.22 / 19516.87
                << " & " <<
                2076.73 / 19516.87
                << " & " <<
                333.21 / 19516.87
                << R"--( \\ \hline )--" << '\n';
                
                ss << Samples[i] << " (Me) & " << 1 << " & " <<
                static_cast<double>(*Nodes[j+2].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
                << " & " <<
                static_cast<double>(*Nodes[j+3].GetResultPtr<ULong64_t>()) / *Nodes[j+2].GetResultPtr<ULong64_t>()
                << " & " <<
                static_cast<double>(*Nodes[j+4].GetResultPtr<ULong64_t>()) / *Nodes[j+3].GetResultPtr<ULong64_t>()
                << " & " <<
                static_cast<double>(*Nodes[j+5].GetResultPtr<ULong64_t>()) / *Nodes[j+4].GetResultPtr<ULong64_t>()
                << R"--( \\ \hline )--" << '\n';
                
                ss << Samples[i] << " (Paper) & " << 1 << " & " <<
                2343.65 / 19516.87
                << " & " <<
                2204.22 / 2343.65
                << " & " <<
                2076.73 / 2204.22
                << " & " <<
                333.21 / 2076.73
                << R"--( \\ \hline )--" << '\n';
                
            }
            
            if (i==4) //5 GeV
            {
                std::cout << Samples[i] << " (Me) & " << 1 << " & " <<
                static_cast<double>(*Nodes[j+2].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
                << " & " <<
                static_cast<double>(*Nodes[j+3].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
                << " & " <<
                static_cast<double>(*Nodes[j+4].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
                << " & " <<
                static_cast<double>(*Nodes[j+5].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
                << R"--( \\ \hline )--" << '\n';
                
                std::cout << Samples[i] << " (Paper) & " << 1 << " & " <<
                3245.86 / 19536.68
                << " & " <<
                3052.81 / 19536.68
                << " & " <<
                2941.28 / 19536.68
                << " & " <<
                1959.68 / 19536.68
                << R"--( \\ \hline )--" << '\n';
                
                ss << Samples[i] << " (Me) & " << 1 << " & " <<
                static_cast<double>(*Nodes[j+2].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
                << " & " <<
                static_cast<double>(*Nodes[j+3].GetResultPtr<ULong64_t>()) / *Nodes[j+2].GetResultPtr<ULong64_t>()
                << " & " <<
                static_cast<double>(*Nodes[j+4].GetResultPtr<ULong64_t>()) / *Nodes[j+3].GetResultPtr<ULong64_t>()
                << " & " <<
                static_cast<double>(*Nodes[j+5].GetResultPtr<ULong64_t>()) / *Nodes[j+4].GetResultPtr<ULong64_t>()
                << R"--( \\ \hline )--" << '\n';
                
                ss << Samples[i] << " (Paper) & " << 1 << " & " <<
                3245.86 / 19536.68
                << " & " <<
                3052.81 / 3245.86
                << " & " <<
                2941.28 / 3052.81
                << " & " <<
                1959.68 / 2941.28
                << R"--( \\ \hline )--" << '\n';
            }
        }
    }
    
    std::cout << R"--(\end{tabular}})--" << "\n\n\n";
    
    ss << R"--(\end{tabular}})--" << "\n\n\n";
    std::cout << ss.str() << '\n';
    
    for (int i = 0, j = 0; (i <= 12 && j <= 2); i += 6, j++)
    {
        totalEventsZgamma += *Nodes[i].GetResultPtr<ULong64_t>() * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
        passPreselectionZgamma += *Nodes[i+1].GetResultPtr<ULong64_t>() * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
        photonPtDeltaRCountZgamma += *Nodes[i+2].GetResultPtr<ULong64_t>() * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
        xWindowZgamma += *Nodes[i+3].GetResultPtr<ULong64_t>() * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
        srCountZgamma += *Nodes[i+4].GetResultPtr<ULong64_t>() * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
        srIDCountZgamma += *Nodes[i+5].GetResultPtr<ULong64_t>() * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
        
        totalEventsZgammaStatUnc += pow(sqrt(*Nodes[i].GetResultPtr<ULong64_t>()) * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
        passPreselectionZgammaStatUnc += pow(sqrt(*Nodes[i+1].GetResultPtr<ULong64_t>()) * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
        photonPtDeltaRCountZgammaStatUnc += pow(sqrt(*Nodes[i+2].GetResultPtr<ULong64_t>()) * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
        xWindowZgammaStatUnc += pow(sqrt(*Nodes[i+3].GetResultPtr<ULong64_t>()) * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
        srCountZgammaStatUnc += pow(sqrt(*Nodes[i+4].GetResultPtr<ULong64_t>()) * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
        srIDCountZgammaStatUnc += pow(sqrt(*Nodes[i+5].GetResultPtr<ULong64_t>()) * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
    }
    
    for (int i = 36, j = 0; (i <= 84 && j <= 8); i += 6, j++)
    {
        totalEventsZjets += *Nodes[i].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
        passPreselectionZjets += *Nodes[i+1].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
        photonPtDeltaRCountZjets += *Nodes[i+2].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
        xWindowZjets += *Nodes[i+3].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
        srCountZjets += *Nodes[i+4].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
        srIDCountZjets += *Nodes[i+5].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
        
        totalEventsZjetsStatUnc += pow(sqrt(*Nodes[i].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
        passPreselectionZjetsStatUnc += pow(sqrt(*Nodes[i+1].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
        photonPtDeltaRCountZjetsStatUnc += pow(sqrt(*Nodes[i+2].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
        xWindowZjetsStatUnc += pow(sqrt(*Nodes[i+3].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
        srCountZjetsStatUnc += pow(sqrt(*Nodes[i+4].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
        srIDCountZjetsStatUnc += pow(sqrt(*Nodes[i+5].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>()),2);
    }
    
    os << R"--(Total $Z\gamma$ & )--" << totalEventsZgamma
    << R"--($\, \pm \,$)--" << sqrt(totalEventsZgammaStatUnc)
    << " & " << passPreselectionZgamma
    << R"--($\, \pm \,$)--" << sqrt(passPreselectionZgammaStatUnc)
    << " & " << photonPtDeltaRCountZgamma
    << R"--($\, \pm \,$)--" << sqrt(photonPtDeltaRCountZgammaStatUnc)
    << " & " << xWindowZgamma
    << R"--($\, \pm \,$)--" << sqrt(xWindowZgammaStatUnc)
    << " & " << srCountZgamma
    << R"--($\, \pm \,$)--" << sqrt(srCountZgammaStatUnc)
    << " & " << srIDCountZgamma
    << R"--($\, \pm \,$)--" << sqrt(srIDCountZgammaStatUnc)
    << R"--( \\ \hline )--" << '\n';
    
    os << R"--(Total $Z$ jets & )--" << totalEventsZjets
    << R"--($\, \pm \,$)--" << sqrt(totalEventsZjetsStatUnc)
    << " & " << passPreselectionZjets
    << R"--($\, \pm \,$)--" << sqrt(passPreselectionZjetsStatUnc)
    << " & " << photonPtDeltaRCountZjets
    << R"--($\, \pm \,$)--" << sqrt(photonPtDeltaRCountZjetsStatUnc)
    << " & " << xWindowZjets
    << R"--($\, \pm \,$)--" << sqrt(xWindowZjetsStatUnc)
    << " & " << srCountZjets
    << R"--($\, \pm \,$)--" << sqrt(srCountZjetsStatUnc)
    << " & " << srIDCountZjets
    << R"--($\, \pm \,$)--" << sqrt(srIDCountZjetsStatUnc)
    << R"--( \\ \hline )--" << '\n';
    
    os << R"--(Total Bkg & )--" << totalEventsZjets+totalEventsZgamma
    << R"--($\, \pm \,$)--" << sqrt(totalEventsZjetsStatUnc+totalEventsZgammaStatUnc)
    << " & " << passPreselectionZjets+passPreselectionZgamma
    << R"--($\, \pm \,$)--" << sqrt(passPreselectionZjetsStatUnc+passPreselectionZgammaStatUnc)
    << " & " << photonPtDeltaRCountZjets+photonPtDeltaRCountZgamma
    << R"--($\, \pm \,$)--" << sqrt(photonPtDeltaRCountZjetsStatUnc+photonPtDeltaRCountZgammaStatUnc)
    << " & " << xWindowZjets+xWindowZgamma
    << R"--($\, \pm \,$)--" << sqrt(xWindowZjetsStatUnc+xWindowZgammaStatUnc)
    << " & " << srCountZjets+srCountZgamma
    << R"--($\, \pm \,$)--" << sqrt(srCountZjetsStatUnc+srCountZgammaStatUnc)
    << " & " << srIDCountZjets+srIDCountZgamma
    << R"--($\, \pm \,$)--" << sqrt(srIDCountZjetsStatUnc+srIDCountZgammaStatUnc)
        << R"--( \\ \hline )--" << '\n';
    
    os << R"--(\end{tabular}})--" << '\n';
    std::cout << os.str() << '\n';
}

void Table3_Displaced_Axions()
{
    std::vector<std::vector<std::string>> input_filenames = {
        //Signal
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root"},
        //Displaced Signal
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_FewMassPoints.root"},
    };
    
    constexpr std::array<const char*,15> Samples = {R"--(Prompt Signal $m_{\text{A}}$ = 1 GeV)--", R"--(Prompt Signal $m_{\text{A}}$ = 5 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 0.2 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 0.5 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 1 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 3 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 5 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 10 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 10.1 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 20 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 20.1 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 29.5 GeV)--", };
    
    std::vector<double> massPoints = {0.2,0.5,1,3,5,10,10.1,20,20.1,29.5};
    
    std::ostringstream os;
    os << R"--(\section*{Table 3 Prompt and Displaced Signal Samples})--" << '\n';
    os << R"--(\hspace{-3cm}\scalebox{0.45}{)--" << '\n';
    os << R"--(\begin{tabular}{|c|c|c|c|c|c|c|c|})--" << '\n';
    os << R"--(\hline)--" << '\n';
    os << R"--(Sample & Before Preselection & 2 leptons
          & Opposite Charge & $p_{T}^{\text{leading}} > 27$ GeV, \; $p_{T}^{\text{sub-leading}} > 20$ GeV & Same flavour & dilep mass cut & dilep $p_{T}$ cut \\ \hline )--" << '\n';
    
    std::vector<ROOT::RDF::RResultHandle> Nodes;
    
    int counter = 0;
    
    for (auto& file: input_filenames)
    {
        SchottDataFrame df(MakeRDF(file, 8));
        
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

        },{"electrons"}).Filter([](RVec<Muon>& muons, RVec<Electron> di_electrons)
        {
            return (di_electrons.size()==2 && muons.empty() && DeltaR(di_electrons[0].Vector(), di_electrons[1].Vector()) > 0.01);
            
        }, {"muons", "di_electrons"});
        
        auto opp_charge = two_leptons.Filter([](RVec<Electron> di_electrons)
        {
            return (di_electrons[0].electron_charge*di_electrons[1].electron_charge < 0);
            
        }, {"di_electrons"});
        
        auto leading_pt = opp_charge.Filter([](RVec<Electron>& electrons)
        {
            return ((electrons[0].electron_pt >= 20e3 && electrons[1].electron_pt >= 27e3) || (electrons[1].electron_pt >= 20e3 && electrons[0].electron_pt >= 27e3));
        }, {"di_electrons"});
        
        auto same_flavour = leading_pt.Filter([] (RVec<Electron>& electrons)
        {
            return true;
        }, {"di_electrons"});
        
        auto dilep = same_flavour.Define("dilep",[] (RVec<Electron>& electrons)
        {
            return (electrons[0].Vector() + electrons[1].Vector());
        }, {"di_electrons"});
        
        auto mass = dilep.Filter([] (PtEtaPhiEVector& dilep)
        {
            auto mass = dilep.M()/1e3;
            return ((mass >= 81) && (mass <= 101));
        }, {"dilep"});
        
        auto pt_cut = mass.Filter([] (PtEtaPhiEVector& dilep)
        {
            auto pT = dilep.Pt()/1e3;
            return pT > 10;
        }, {"dilep"});
        
        if (counter < 2)
        {
            Nodes.push_back(trigger_selection.Count());
            Nodes.push_back(two_leptons.Count());
            Nodes.push_back(opp_charge.Count());
            Nodes.push_back(leading_pt.Count());
            Nodes.push_back(same_flavour.Count());
            Nodes.push_back(mass.Count());
            Nodes.push_back(pt_cut.Count());
        }
        
        else
        {
            for (auto& mass_point: massPoints)
            {
                auto mass_point_trigger_selection = trigger_selection.Filter([&]
                (float axion_mass)
                {
                    return (roundToOneDecimalPlace(axion_mass) == mass_point);
                    
                }, {"axion_masses"});
                
                auto mass_point_two_leptons = two_leptons.Filter([&]
                (float axion_mass)
                {
                    return (roundToOneDecimalPlace(axion_mass) == mass_point);
                    
                }, {"axion_masses"});
                
                auto mass_point_opp_charge = opp_charge.Filter([&]
                (float axion_mass)
                {
                    return (roundToOneDecimalPlace(axion_mass) == mass_point);
                    
                }, {"axion_masses"});
                
                auto mass_point_leading_pt = leading_pt.Filter([&]
                (float axion_mass)
                {
                    return (roundToOneDecimalPlace(axion_mass) == mass_point);
                    
                }, {"axion_masses"});
                
                auto mass_point_same_flavour = same_flavour.Filter([&]
                (float axion_mass)
                {
                    return (roundToOneDecimalPlace(axion_mass) == mass_point);
                    
                }, {"axion_masses"});
                
                auto mass_point_mass = mass.Filter([&]
                (float axion_mass)
                {
                    return (roundToOneDecimalPlace(axion_mass) == mass_point);
                    
                }, {"axion_masses"});
                
                auto mass_point_pt_cut = pt_cut.Filter([&]
                (float axion_mass)
                {
                    return (roundToOneDecimalPlace(axion_mass) == mass_point);
                    
                }, {"axion_masses"});
                
                Nodes.push_back(mass_point_trigger_selection.Count());
                Nodes.push_back(mass_point_two_leptons.Count());
                Nodes.push_back(mass_point_opp_charge.Count());
                Nodes.push_back(mass_point_leading_pt.Count());
                Nodes.push_back(mass_point_same_flavour.Count());
                Nodes.push_back(mass_point_mass.Count());
                Nodes.push_back(mass_point_pt_cut.Count());
                
            }
        }
        
        counter++;
    }

    constexpr std::array<const char*,7> Cuts = {"total", "two leptons", "opposite charge", "leading pt", "same flavour", "mass", "pt cut"};
    
    std::cout << Nodes.size();
    
    ROOT::RDF::RunGraphs(Nodes); // running all computation nodes concurrently
    
//   0    1   2   3   4   5   6                          ma1
//   7    8   9   10  11  12  13                         ma5
//   14   15  16  17  18  19  20                         displaced_axion_1
//   21   22  23  24  25  26  27                         displaced_axion_2
//   28   29  30  31  32  33  34                         displaced_axion_3
//   35   36  37  38  39  40  41                         displaced_axion_4
//   42   43  44  45  46  47  48                         displaced_axion_5
//   49   50  51  52  53  54  55                         displaced_axion_6
//   56   57  58  59  60  61  62                         displaced_axion_7
//   63   64  65  66  67  68  69                         displaced_axion_8
//   70   71  72  73  74  75  76                         displaced_axion_9
//   77   78  79  80  81  82  83                         displaced_axion_10

    std::cout << R"--(\section*{Table 3 Signal Ratios})--" << '\n';
    std::cout << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
    std::cout << R"--(\begin{tabular}{|c|c|c|c|c|c|c|})--" << '\n';
    std::cout << R"--(\hline)--" << '\n';
    std::cout << R"--(Sample & $\frac{\text{2 leptons}}{\text{2 leptons}}$
          & $\frac{\text{Opposite Charge}}{\text{2 leptons}}$ & $\frac{p_{T}^{\text{leading}} > 27\text{ GeV, } \, p_{T}^{\text{sub-leading}} > 20 \text{ GeV }}{\text{2 leptons}}$ & $\frac{\text{Same flavour}}{\text{2 leptons}}$ & $\frac{\text{dilep mass cut}}{\text{2 leptons}}$ & $\frac{\text{dilep }p_{T} \text{ cut}}{\text{2 leptons}}$ \\ \hline )--" << '\n';
    
    os.setf(std::ios::fixed);
    os.precision(2);
        
    for (int i=0, j=0; (i<12 && j <= 77); i++, j+=7)
    {
        os << Samples[i];
        os  << " & " << *Nodes[j].GetResultPtr<ULong64_t>()
        << R"--($\, \pm \,$)--" << sqrt(*Nodes[j].GetResultPtr<ULong64_t>())
        << " & " << *Nodes[j+1].GetResultPtr<ULong64_t>()
        << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+1].GetResultPtr<ULong64_t>())
        << " & " << *Nodes[j+2].GetResultPtr<ULong64_t>()
        << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+2].GetResultPtr<ULong64_t>())
        << " & " << *Nodes[j+3].GetResultPtr<ULong64_t>()
        << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+3].GetResultPtr<ULong64_t>())
        << " & " << *Nodes[j+4].GetResultPtr<ULong64_t>()
        << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+4].GetResultPtr<ULong64_t>())
        << " & " << *Nodes[j+5].GetResultPtr<ULong64_t>()
        << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+5].GetResultPtr<ULong64_t>())
        << " & " << *Nodes[j+6].GetResultPtr<ULong64_t>()
        << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+6].GetResultPtr<ULong64_t>())
        << R"--( \\ \hline )--" << '\n';
        
        if (i==0) //1 GeV
        {
            std::cout << Samples[i] << " (Me) & " << 1 << " & " <<
            static_cast<double>(*Nodes[j+2].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
            << " & " <<
            static_cast<double>(*Nodes[j+3].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
            << " & " <<
            static_cast<double>(*Nodes[j+4].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
            << " & " <<
            static_cast<double>(*Nodes[j+5].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
            << " & " <<
            static_cast<double>(*Nodes[j+6].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
            << R"--( \\ \hline )--" << '\n';
            
            std::cout << Samples[i] << " (Paper) & " << 1 << " & " <<
            21505.03 / 21606.75
            << " & " <<
            21375.48 / 21606.75
            << " & " <<
            21375.33 / 21606.75
            << " & " <<
            20543.06 / 21606.75
            << " & " <<
            19516.87 / 21606.75
            << R"--( \\ \hline )--" << '\n';
        }
        
        else if (i==1) //5 GeV
        {
            std::cout << Samples[i] << " (Me) & " << 1 << " & " <<
            static_cast<double>(*Nodes[j+2].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
            << " & " <<
            static_cast<double>(*Nodes[j+3].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
            << " & " <<
            static_cast<double>(*Nodes[j+4].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
            << " & " <<
            static_cast<double>(*Nodes[j+5].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
            << " & " <<
            static_cast<double>(*Nodes[j+6].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
            << R"--( \\ \hline )--" << '\n';
            
            std::cout << Samples[i] << " (Paper) & " << 1 << " & " <<
            21556.03 / 21655.38
            << " & " <<
            21424.29 / 21655.38
            << " & " <<
            21424.28 / 21655.38
            << " & " <<
            20585.09 / 21655.38
            << " & " <<
            19536.68 / 21655.38
            << R"--( \\ \hline )--" << '\n';
        }
        
        else
        {
            std::cout << Samples[i] << " Displaced (Me) & " << 1 << " & " <<
            static_cast<double>(*Nodes[j+2].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
            << " & " <<
            static_cast<double>(*Nodes[j+3].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
            << " & " <<
            static_cast<double>(*Nodes[j+4].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
            << " & " <<
            static_cast<double>(*Nodes[j+5].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
            << " & " <<
            static_cast<double>(*Nodes[j+6].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
            << R"--( \\ \hline )--" << '\n';
        }
    }

    std::cout << R"--(\end{tabular}})--" << "\n\n\n";

    os << R"--(\end{tabular}})--" << '\n';
    std::cout << os.str() << '\n';
}

void Table8_Displaced_Axions()
{
    std::vector<std::vector<std::string>> input_filenames = {
        //Signal
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root"},
        //Displaced Signal
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_FewMassPoints.root"},
    };
    
    constexpr std::array<const char*,15> Samples = {R"--(Prompt Signal $m_{\text{A}}$ = 1 GeV)--", R"--(Prompt Signal $m_{\text{A}}$ = 5 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 0.2 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 0.5 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 1 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 3 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 5 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 10 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 10.1 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 20 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 20.1 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 29.5 GeV)--", };
    
    std::vector<double> massPoints = {0.2,0.5,1,3,5,10,10.1,20,20.1,29.5};
    
    std::vector<ROOT::RDF::RResultHandle> Nodes;
    
    std::ostringstream os;
    os << R"--(\section*{Table 8 Prompt and Displaced Signal Samples})--" << '\n';
    os << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
    os << R"--(\begin{tabular}{|c|c|c|c|c|})--" << '\n';
    os << R"--(\hline)--" << '\n';
    os << R"--(Sample & Total Events & PS: Resolved & SB & SR
           \\ \hline )--" << '\n';
    
    int counter = 0;
    
    for (auto& file: input_filenames)
    {
        SchottDataFrame df(MakeRDF(file, 8));
        
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
        
        auto resolved = ptCut.Define("photons_pass_cuts",
        [&](RVec<Photon> photons)
        {
            photons.erase(std::remove_if(photons.begin(),photons.end(),
            [](Photon& x)
            {
                return ((abs(x.photon_eta) >= 2.37) || (abs(x.photon_eta) > 1.37 && abs(x.photon_eta) < 1.52) || (!x.photon_id_loose));

            }), photons.end());
            return photons;
        }, {"photons"}).Filter(
        [&](RVec<Photon>& reco_photons_matched)
        {
            if (reco_photons_matched.size() < 2)
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
                return true;
            }
            return false;

        }, {"photons_pass_cuts"});

        auto SR = resolved.Filter(
        [&](RVec<Electron>& electrons, RVec<Photon>& reco_photons_matched)
        {
            auto mass = ((electrons[0].Vector()+electrons[1].Vector()).M() + (reco_photons_matched[0].Vector()+reco_photons_matched[1].Vector()).M())/1e3;
            return ((mass > 110) && (mass < 140));
        }, {"di_electrons", "photons_pass_cuts"});
        
        auto SB = resolved.Filter(
        [&](RVec<Electron>& electrons, RVec<Photon>& reco_photons_matched)
        {
            auto mass = ((electrons[0].Vector()+electrons[1].Vector()).M() + (reco_photons_matched[0].Vector()+reco_photons_matched[1].Vector()).M())/1e3;
            return (!((mass > 110) && (mass < 140)));
        }, {"di_electrons", "photons_pass_cuts"});
        
        if (counter < 2)
        {
            Nodes.push_back(trigger_selection.Count());
            Nodes.push_back(resolved.Count());
            Nodes.push_back(SB.Count());
            Nodes.push_back(SR.Count());
        }
        else
        {
            for (auto& mass_point: massPoints)
            {
                auto mass_point_trigger_selection = trigger_selection.Filter([&]
                (float axion_mass)
                {
                    return (roundToOneDecimalPlace(axion_mass) == mass_point);
                    
                }, {"axion_masses"});
                
                auto mass_point_resolved = resolved.Filter([&]
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
                
                Nodes.push_back(mass_point_trigger_selection.Count());
                Nodes.push_back(mass_point_resolved.Count());
                Nodes.push_back(mass_point_SB.Count());
                Nodes.push_back(mass_point_SR.Count());
            }
        }
        
        counter++;
    }
    
    ROOT::RDF::RunGraphs(Nodes); // running all computation nodes concurrently
    
//       0    1   2   3                           ma1
//       4    5   6   7                           ma5
//       8    9   10  11                          displaced_axion_1
//       12   13  14  15                          displaced_axion_2
//       16   17  18  19                          displaced_axion_3
//       20   21  22  23                          displaced_axion_4
//       24   25  26  27                          displaced_axion_5
//       28   29  30  31                          displaced_axion_6
//       32   33  34  35                          displaced_axion_7
//       36   37  38  39                          displaced_axion_8
//       40   41  42  43                          displaced_axion_9
//       44   45  46  47                          displaced_axion_10
    os.setf(std::ios::fixed);
    os.precision(2);
 
    for (int i=0, j=0; (i<12 && j <= 44); i++, j+=4)
    {
        os << Samples[i];
        
        os  << " & " << *Nodes[j].GetResultPtr<ULong64_t>()
            << R"--($\, \pm \,$)--" << sqrt(*Nodes[j].GetResultPtr<ULong64_t>())
            << " & " << *Nodes[j+1].GetResultPtr<ULong64_t>()
            << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+1].GetResultPtr<ULong64_t>())
            << " & " << *Nodes[j+2].GetResultPtr<ULong64_t>()
            << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+2].GetResultPtr<ULong64_t>())
            << " & " << *Nodes[j+3].GetResultPtr<ULong64_t>()
            << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+3].GetResultPtr<ULong64_t>())
            << R"--( \\ \hline )--" << '\n';

    }
    
    os << R"--(\end{tabular}})--" << '\n';
    std::cout << os.str() << '\n';
}


void Table11_Displaced_Axions()
{
    std::vector<std::vector<std::string>> input_filenames = {
        //Signal
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root"},
        //Displaced Signal
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_FewMassPoints.root"},
    };
    
    constexpr std::array<const char*,15> Samples = {R"--(Prompt Signal $m_{\text{A}}$ = 1 GeV)--", R"--(Prompt Signal $m_{\text{A}}$ = 5 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 0.2 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 0.5 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 1 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 3 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 5 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 10 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 10.1 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 20 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 20.1 GeV)--", R"--(Displaced Signal $m_{\text{A}}$ = 29.5 GeV)--", };
    
    std::vector<double> massPoints = {0.2,0.5,1,3,5,10,10.1,20,20.1,29.5};
    
    std::vector<ROOT::RDF::RResultHandle> Nodes;

    std::ostringstream os, ss;
    
    os.setf(std::ios::fixed);
    os.precision(2);
    
    os << R"--(\section*{Table 11 Prompt and Displaced Signal Samples})--" << '\n';
    os << R"--(\hspace{-3cm}\scalebox{0.55}{)--" << '\n';
    os << R"--(\begin{tabular}{|c|c|c|c|c|c|c|})--" << '\n';
    os << R"--(\hline)--" << '\n';
    os << R"--(Sample & Total Events & pass preselection (PS) & photon $p_T$ + $\Delta R_{\gamma\gamma}$ cut & $X$ window & SR & SR-ID
           \\ \hline )--" << '\n';
    
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
    
    int counter = 0;

    for (auto& file: input_filenames)
    {
        SchottDataFrame df(MakeRDF(file, 8));
        
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
        [](RVec<Electron> electrons)
        {
            electrons.erase(std::remove_if(electrons.begin(),electrons.end(),
            [](Electron& ep)
            {
                return ((ep.electron_pt/1e3 <= 20) || (abs(ep.electron_eta) >= 2.37)
                        || ((1.37 < abs(ep.electron_eta)) && (abs(ep.electron_eta) < 1.52))
                        || (ep.electron_id_medium != 1));
                

            }), electrons.end());

            return electrons;

        },{"electrons"}).Filter([](RVec<Muon>& muons, RVec<Electron> di_electrons)
        {
            return (di_electrons.size()==2 && muons.empty() && DeltaR(di_electrons[0].Vector(), di_electrons[1].Vector()) > 0.01);
            
        }, {"muons", "di_electrons"});
        
        auto opp_charge = two_leptons.Filter([](RVec<Electron> di_electrons)
        {
            return (di_electrons[0].electron_charge*di_electrons[1].electron_charge < 0);
            
        }, {"di_electrons"});
        
        auto leading_pt = opp_charge.Filter([](RVec<Electron>& electrons)
        {
            return ((electrons[0].electron_pt >= 20e3 && electrons[1].electron_pt >= 27e3) || (electrons[1].electron_pt >= 20e3 && electrons[0].electron_pt >= 27e3));
        }, {"di_electrons"});
        
//        auto delta_R = leading_pt.Filter([] (RVec<Electron>& electrons)
//        {
//            return (DeltaR(electrons[0].Vector(), electrons[1].Vector()) > 0.01);
//        }, {"di_electrons"});
        
        auto same_flavour = leading_pt.Filter([] (RVec<Electron>& electrons)
        {
            return true; //abs(electrons[0].electron_pdg_id) == abs(electrons[1].electron_pdg_id) == 11;
        }, {"di_electrons"});
        
        auto dilep = same_flavour.Define("dilep",[] (RVec<Electron>& electrons)
        {
            return (electrons[0].Vector() + electrons[1].Vector());
        }, {"di_electrons"});
        
        auto mass = dilep.Filter([] (PtEtaPhiEVector& dilep)
        {
            auto mass = dilep.M()/1e3;
            return ((mass >= 81) && (mass <= 101));
        }, {"dilep"});
        
        auto ptCut = mass.Filter([] (PtEtaPhiEVector& dilep)
        {
            auto pT = dilep.Pt()/1e3;
            return pT > 10;
        }, {"dilep"});
        
        auto photonPtDeltaR = ptCut.Define("photonPtDeltaR",
        [&](RVec<Photon> photons)
        {
            photons.erase(std::remove_if(photons.begin(),photons.end(),
            [](Photon& x)
            {
                return ((abs(x.photon_eta) >= 2.37) || (x.photon_pt <= 10e3) || (abs(x.photon_eta) > 1.37 && abs(x.photon_eta) < 1.52));
                
            }), photons.end());
            return photons;
        }, {"photons"}).Filter(
       [&](RVec<Photon>& reco_photons_matched)
       {
            if (reco_photons_matched.size() < 2)
            {
              return false;
            }
            RVec<Photon> x;
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
                    x = {reco_photons_matched[combs[0][i]], reco_photons_matched[combs[1][i]]};
                }
            }
            return (chosen_delta_r < 1.5 && pt1 > 10e3 && pt2 > 10e3);

       }, {"photonPtDeltaR"});
                
        auto X_window = photonPtDeltaR.Define("chosen_two",
        [](RVec<Photon>& reco_photons_matched)
        {
            RVec<Photon> x;
            if (reco_photons_matched.size() < 2)
            {
                return x;
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
                    x = {reco_photons_matched[combs[0][i]], reco_photons_matched[combs[1][i]]};
                }
            }
            if (pt1 > 10e3 && pt2 > 10e3 && best_X > 0.96 && best_X < 1.2 && chosen_delta_r < 1.5)
            {
                return x;
            }
            x.clear();
            return x;
        }, {"photonPtDeltaR"}).Filter(
        [&](RVec<Photon>& reco_photons_matched)
        {
            return (reco_photons_matched.size()==2);
        }, {"chosen_two"});
                
        auto SR = X_window.Filter(
        [](RVec<Photon>& photons, RVec<Electron>& electrons)
        {
            PtEtaPhiEVector photonVec = photons[0].Vector() + photons[1].Vector();
            PtEtaPhiEVector electronVec = electrons[0].Vector() + electrons[1].Vector();
            auto mass = (photonVec+electronVec).M()/1e3;
            return ((mass >= 110) && (mass <= 140));
        },{"chosen_two","di_electrons"});
                
        auto SR_ID = SR.Filter(
        [](RVec<Photon>& photons)
        {
            return (photons[0].photon_id_loose && photons[1].photon_id_loose);
        },{"chosen_two"});
        
        if (counter < 2)
        {
            Nodes.push_back(trigger_selection.Count());
            Nodes.push_back(ptCut.Count());
            Nodes.push_back(photonPtDeltaR.Count());
            Nodes.push_back(X_window.Count());
            Nodes.push_back(SR.Count());
            Nodes.push_back(SR_ID.Count());
        }
        else
        {
            for (auto& mass_point: massPoints)
            {
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
                
                Nodes.push_back(mass_point_trigger_selection.Count());
                Nodes.push_back(mass_point_ptCut.Count());
                Nodes.push_back(mass_point_photonPtDeltaR.Count());
                Nodes.push_back(mass_point_X_window.Count());
                Nodes.push_back(mass_point_SR.Count());
                Nodes.push_back(mass_point_SR_ID.Count());
            }
        }
        
        counter++;
    }
    
    ROOT::RDF::RunGraphs(Nodes); // running all computation nodes concurrently
    
//    0    1   2   3   4   5                            ma1
//    6    7   8   9   10  11                           ma5
//    12   13  14  15  16  17                           displaced_axion_1
//    18   19  20  21  22  23                           displaced_axion_2
//    24   25  26  27  28  29                           displaced_axion_3
//    30   31  32  33  34  35                           displaced_axion_4
//    36   37  38  39  40  41                           displaced_axion_5
//    42   43  44  45  46  47                           displaced_axion_6
//    48   49  50  51  52  53                           displaced_axion_7
//    54   55  56  57  58  59                           displaced_axion_8
//    60   61  62  63  64  65                           displaced_axion_9
//    66   67  68  69  70  71                           displaced_axion_10

    std::cout << R"--(\section*{Table 11 Signal Ratios})--" << '\n';
    std::cout << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
    std::cout << R"--(\begin{tabular}{|c|c|c|c|c|c|c|})--" << '\n';
    std::cout << R"--(\hline)--" << '\n';
    std::cout << R"--(Sample & $\frac{\text{pass preselection (PS)}}{\text{pass preselection (PS)}}$ & $\frac{\text{photon } p_T \text{ + } \Delta R_{\gamma\gamma} \text{ cut }}{\text{pass preselection (PS)}}$ & $\frac{X \text{ window}}{\text{pass preselection (PS)}}$ & $\frac{\text{SR}}{\text{pass preselection (PS)}}$ & $\frac{\text{SR-ID}}{\text{pass preselection (PS)}}$
           \\ \hline )--" << '\n';
    
    ss << R"--(\section*{Table 11 Signal Ratios})--" << '\n';
    ss << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
    ss << R"--(\begin{tabular}{|c|c|c|c|c|c|c|})--" << '\n';
    ss << R"--(\hline)--" << '\n';
    ss << R"--(Sample & $\frac{\text{pass preselection (PS)}}{\text{pass preselection (PS)}}$ & $\frac{\text{photon } p_T \text{ + } \Delta R_{\gamma\gamma} \text{ cut }}{\text{pass preselection (PS)}}$ & $\frac{X \text{ window}}{\text{photon } p_T \text{ + } \Delta R_{\gamma\gamma} \text{ cut }}$ & $\frac{\text{SR}}{X \text{ window}}$ & $\frac{\text{SR-ID}}{\text{SR}}$
           \\ \hline )--" << '\n';
    
    for (int i=0, j=0; (i<12 && j <= 66); i++, j+=6)
    {
        os << Samples[i];
        
        os  << " & " << *Nodes[j].GetResultPtr<ULong64_t>()
            << R"--($\, \pm \,$)--" << sqrt(*Nodes[j].GetResultPtr<ULong64_t>())
            << " & " << *Nodes[j+1].GetResultPtr<ULong64_t>()
            << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+1].GetResultPtr<ULong64_t>())
            << " & " << *Nodes[j+2].GetResultPtr<ULong64_t>()
            << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+2].GetResultPtr<ULong64_t>())
            << " & " << *Nodes[j+3].GetResultPtr<ULong64_t>()
            << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+3].GetResultPtr<ULong64_t>())
            << " & " << *Nodes[j+4].GetResultPtr<ULong64_t>()
            << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+4].GetResultPtr<ULong64_t>())
            << " & " << *Nodes[j+5].GetResultPtr<ULong64_t>()
            << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+5].GetResultPtr<ULong64_t>())
            << R"--( \\ \hline )--" << '\n';
            
        if (i==0) //1 GeV
        {
            std::cout << Samples[i] << " (Me) & " << 1 << " & " <<
            static_cast<double>(*Nodes[j+2].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
            << " & " <<
            static_cast<double>(*Nodes[j+3].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
            << " & " <<
            static_cast<double>(*Nodes[j+4].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
            << " & " <<
            static_cast<double>(*Nodes[j+5].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
            << R"--( \\ \hline )--" << '\n';
            
            std::cout << Samples[i] << " (Paper) & " << 1 << " & " <<
            2343.65 / 19516.87
            << " & " <<
            2204.22 / 19516.87
            << " & " <<
            2076.73 / 19516.87
            << " & " <<
            333.21 / 19516.87
            << R"--( \\ \hline )--" << '\n';
            
            ss << Samples[i] << " (Me) & " << 1 << " & " <<
            static_cast<double>(*Nodes[j+2].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
            << " & " <<
            static_cast<double>(*Nodes[j+3].GetResultPtr<ULong64_t>()) / *Nodes[j+2].GetResultPtr<ULong64_t>()
            << " & " <<
            static_cast<double>(*Nodes[j+4].GetResultPtr<ULong64_t>()) / *Nodes[j+3].GetResultPtr<ULong64_t>()
            << " & " <<
            static_cast<double>(*Nodes[j+5].GetResultPtr<ULong64_t>()) / *Nodes[j+4].GetResultPtr<ULong64_t>()
            << R"--( \\ \hline )--" << '\n';
            
            ss << Samples[i] << " (Paper) & " << 1 << " & " <<
            2343.65 / 19516.87
            << " & " <<
            2204.22 / 2343.65
            << " & " <<
            2076.73 / 2204.22
            << " & " <<
            333.21 / 2076.73
            << R"--( \\ \hline )--" << '\n';
            
        }
        
        else if (i==1) //5 GeV
        {
            std::cout << Samples[i] << " (Me) & " << 1 << " & " <<
            static_cast<double>(*Nodes[j+2].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
            << " & " <<
            static_cast<double>(*Nodes[j+3].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
            << " & " <<
            static_cast<double>(*Nodes[j+4].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
            << " & " <<
            static_cast<double>(*Nodes[j+5].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
            << R"--( \\ \hline )--" << '\n';
            
            std::cout << Samples[i] << " (Paper) & " << 1 << " & " <<
            3245.86 / 19536.68
            << " & " <<
            3052.81 / 19536.68
            << " & " <<
            2941.28 / 19536.68
            << " & " <<
            1959.68 / 19536.68
            << R"--( \\ \hline )--" << '\n';
            
            ss << Samples[i] << " (Me) & " << 1 << " & " <<
            static_cast<double>(*Nodes[j+2].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
            << " & " <<
            static_cast<double>(*Nodes[j+3].GetResultPtr<ULong64_t>()) / *Nodes[j+2].GetResultPtr<ULong64_t>()
            << " & " <<
            static_cast<double>(*Nodes[j+4].GetResultPtr<ULong64_t>()) / *Nodes[j+3].GetResultPtr<ULong64_t>()
            << " & " <<
            static_cast<double>(*Nodes[j+5].GetResultPtr<ULong64_t>()) / *Nodes[j+4].GetResultPtr<ULong64_t>()
            << R"--( \\ \hline )--" << '\n';
            
            ss << Samples[i] << " (Paper) & " << 1 << " & " <<
            3245.86 / 19536.68
            << " & " <<
            3052.81 / 3245.86
            << " & " <<
            2941.28 / 3052.81
            << " & " <<
            1959.68 / 2941.28
            << R"--( \\ \hline )--" << '\n';
        }
        
        else
        {
            std::cout << Samples[i] << " Displaced (Me) & " << 1 << " & " <<
            static_cast<double>(*Nodes[j+2].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
            << " & " <<
            static_cast<double>(*Nodes[j+3].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
            << " & " <<
            static_cast<double>(*Nodes[j+4].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
            << " & " <<
            static_cast<double>(*Nodes[j+5].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
            << R"--( \\ \hline )--" << '\n';
            
            ss << Samples[i] << " Displaced (Me) & " << 1 << " & " <<
            static_cast<double>(*Nodes[j+2].GetResultPtr<ULong64_t>()) / *Nodes[j+1].GetResultPtr<ULong64_t>()
            << " & " <<
            static_cast<double>(*Nodes[j+3].GetResultPtr<ULong64_t>()) / *Nodes[j+2].GetResultPtr<ULong64_t>()
            << " & " <<
            static_cast<double>(*Nodes[j+4].GetResultPtr<ULong64_t>()) / *Nodes[j+3].GetResultPtr<ULong64_t>()
            << " & " <<
            static_cast<double>(*Nodes[j+5].GetResultPtr<ULong64_t>()) / *Nodes[j+4].GetResultPtr<ULong64_t>()
            << R"--( \\ \hline )--" << '\n';
        }
        
    }
    
    std::cout << R"--(\end{tabular}})--" << "\n\n\n";
    
    ss << R"--(\end{tabular}})--" << "\n\n\n";
    std::cout << ss.str() << '\n';
    
    os << R"--(\end{tabular}})--" << '\n';
    std::cout << os.str() << '\n';
}


//void Coupling()
//{
//    std::vector<std::vector<std::string>> input_filenames = {
//        //Signal
////        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root"},
////        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root"},
//        //Displaced Signal
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_FewMassPoints.root"},
//    };
//
//    std::vector<double> massPoints = {0.2,0.5,1,3,5,10,10.1,20,20.1,29.5};
//    std::vector<double> ALP_photon_couplings = {1.0};
//    std::vector<ROOT::RDF::RResultHandle> Nodes;
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
//        auto newDf = df.Define("EventWeight",
//        [](const RVec<float>& ei_event_weights_generator)
//        {
//            return  ((ei_event_weights_generator[0]) ? 1 / ei_event_weights_generator[0] : 1);
//
//        }, {"ei_event_weights_generator"}).Define("totEventWeight", [](RVec<float> photon_id_eff, RVec<float> photon_iso_eff, RVec<float> photon_trg_eff/*, RVec<float> ei_event_weights_generator*/)
//        {
//            auto ResizeVal = std::max({photon_id_eff.size(), photon_iso_eff.size(), photon_trg_eff.size()});
//            photon_id_eff.resize(ResizeVal,1);
//            photon_iso_eff.resize(ResizeVal,1);
//            photon_trg_eff.resize(ResizeVal,1);
//
//            return photon_id_eff*photon_iso_eff*photon_trg_eff;//*ei_event_weights_generator[0];
//
//        }, {"photon_id_eff", "photon_iso_eff", "photon_trg_eff",/* "ei_event_weights_generator"*/}).Define("truth_axions", [&](RVec<TruthParticle> truth_particles)
//        {
//           truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
//           [](TruthParticle& x)
//           {
//              return (abs(x.mc_pdg_id) != 36);
//
//           }), truth_particles.end());
//
//           return truth_particles;
//
//        }, {"truth_particles"}).Define("axion_masses", [&](RVec<TruthParticle>& truth_axions)
//        {
//           return truth_axions[0].mc_mass/1e3f;
//
//        }, {"truth_axions"});
//
//        auto trigger_selection = newDf
//        .Filter([](const RVec<std::string>& trigger_passed_triggers)
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
//        [](RVec<Electron> electrons)
//        {
//            electrons.erase(std::remove_if(electrons.begin(),electrons.end(),
//            [](Electron& ep)
//            {
//                return ((ep.electron_pt/1e3 <= 20) || (abs(ep.electron_eta) >= 2.37)
//                        || ((1.37 < abs(ep.electron_eta)) && (abs(ep.electron_eta) < 1.52))
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
//            return true; //abs(electrons[0].electron_pdg_id) == abs(electrons[1].electron_pdg_id) == 11;
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
//                return ((abs(x.photon_eta) >= 2.37) || (x.photon_pt <= 10e3) || (abs(x.photon_eta) > 1.37 && abs(x.photon_eta) < 1.52));
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
//                if (i==0 || abs(1-X) < abs(1-best_X))
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
//                if (i==0 || abs(1-X) < abs(1-best_X))
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
//        for (auto& mass_point: massPoints)
//        {
//            auto mass_point_newDf = newDf.Filter([&]
//            (float axion_mass)
//            {
//                return (roundToOneDecimalPlace(axion_mass) == mass_point);
//
//            }, {"axion_masses"});
//
//            auto mass_point_trigger_selection = trigger_selection.Filter([&]
//            (float axion_mass)
//            {
//                return (roundToOneDecimalPlace(axion_mass) == mass_point);
//
//            }, {"axion_masses"});
//
//            auto mass_point_ptCut = ptCut.Filter([&]
//            (float axion_mass)
//            {
//                return (roundToOneDecimalPlace(axion_mass) == mass_point);
//
//            }, {"axion_masses"});
//
//            auto mass_point_photonPtDeltaR = photonPtDeltaR.Filter([&]
//            (float axion_mass)
//            {
//                return (roundToOneDecimalPlace(axion_mass) == mass_point);
//
//            }, {"axion_masses"});
//
//            auto mass_point_X_window = X_window.Filter([&]
//            (float axion_mass)
//            {
//                return (roundToOneDecimalPlace(axion_mass) == mass_point);
//
//            }, {"axion_masses"});
//
//            auto mass_point_SR = SR.Filter([&]
//            (float axion_mass)
//            {
//                return (roundToOneDecimalPlace(axion_mass) == mass_point);
//
//            }, {"axion_masses"});
//
//            auto mass_point_SR_ID = SR_ID.Filter([&]
//            (float axion_mass)
//            {
//                return (roundToOneDecimalPlace(axion_mass) == mass_point);
//
//            }, {"axion_masses"});
//
////            Nodes.push_back(mass_point_newDf.Count());
////            Nodes.push_back(mass_point_trigger_selection.Count());
////            Nodes.push_back(mass_point_ptCut.Count());
////            Nodes.push_back(mass_point_photonPtDeltaR.Count());
////            Nodes.push_back(mass_point_X_window.Count());
////            Nodes.push_back(mass_point_SR.Count());
////            Nodes.push_back(mass_point_SR_ID.Count());
//
//
//
//
//            Nodes.push_back(mass_point_newDf.Sum<RVec<float>>("totEventWeight"));
//            Nodes.push_back(mass_point_trigger_selection.Sum<RVec<float>>("totEventWeight"));
//            Nodes.push_back(mass_point_ptCut.Sum<RVec<float>>("totEventWeight"));
//            Nodes.push_back(mass_point_photonPtDeltaR.Sum<RVec<float>>("totEventWeight"));
//            Nodes.push_back(mass_point_X_window.Sum<RVec<float>>("totEventWeight"));
//            Nodes.push_back(mass_point_SR.Sum<RVec<float>>("totEventWeight"));
//            Nodes.push_back(mass_point_SR_ID.Sum<RVec<float>>("totEventWeight"));
//
//        }
//
//        counter++;
//    }
//
//    ROOT::RDF::RunGraphs(Nodes); // running all computation nodes concurrently
//    std::cout << Nodes.size() << '\n';
//
//    gStyle->SetPalette(1);
//    TCanvas* c1 = new TCanvas();
////    gPad->SetFrameLineWidth(3);
////    TLegend* legend = new TLegend(0.6, 0.2, 0.85, 0.525);
//    TH2D* histo2d = new TH2D("C_{a#gamma#gamma}", ";m_{a} [GeV];C_{a#gamma#gamma};Efficiency", 295, 0.2, 29.5, 5, 0, 2);
//    histo2d->GetXaxis()->SetTitleSize(0.04);
//    histo2d->GetYaxis()->SetTitleSize(0.04);
//    histo2d->GetZaxis()->SetTitleSize(0.04);
//
//    histo2d->GetYaxis()->SetTitleOffset(0.7);
//    histo2d->GetZaxis()->SetTitleOffset(0.6);
//
//    histo2d->GetYaxis()->CenterTitle(true);
//    histo2d->GetZaxis()->SetTitle("Efficiency");
//
//    std::cout << "here?\n\n";
//    for (int i=0, j=0; (i < massPoints.size()); i++, j+=7)
//    {
////        std::cout << i << ' ' << j << ' ' << Nodes.size() << '\n';
//
//        int bin_number = histo2d->FindFixBin(massPoints[i], ALP_photon_couplings[0]);
//
//        //efficiency is the number of events that passed divided by the total for a given mass point
////        double efficiency = *Nodes[j].GetResultPtr<ULong64_t>()
////        ? static_cast<double>(*Nodes[j+5].GetResultPtr<ULong64_t>()) / static_cast<double>(*Nodes[j].GetResultPtr<ULong64_t>())
////        : 0.0;
//        float efficiency = *Nodes[j].GetResultPtr<float>()
//        ? (*Nodes[j+5].GetResultPtr<float>()) / (*Nodes[j].GetResultPtr<float>())
//        : 0.0;
//
//        {
//            std::cout << massPoints[i] << ' '
//            << efficiency << '\n';
//        }
//
//        histo2d->SetBinContent(bin_number, efficiency);
//    }
////    std::cout << "here?\n\n";
//    histo2d->Draw("COLZ");
//
////    histo2d->GetZaxis()->SetTitle("Efficiency");
//
//    TLatex Tl;
//    gStyle->SetOptStat(0);
//    Tl.SetTextSize(0.03);
//    Tl.DrawLatexNDC(0.6, 0.83, "#it{ATLAS} Internal");
//    Tl.DrawLatexNDC(0.6, 0.73,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
//
////    c1->Resize();
////    gPad->SetLeftMargin(0.07);
//    c1->SetCanvasSize(2.5*c1->GetWw(), c1->GetWh());
//    gPad->SetLeftMargin(0.06);
////    gPad->SetRightMargin(0.2);
//    c1->SaveAs("C_{a#gamma#gamma}.pdf");
//
//}

void Coupling_and_Systematics()
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

    
    
    std::vector<std::vector<std::string>> input_filenames = {
        //Signal
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root"},
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root"},
        //Displaced Signal
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_FewMassPoints.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_FewMassPoints_PhotonSFs.root"},
    };
    
    std::vector<float> massPoints = {0.2, 0.5, 1, 3, 5, 10, 20, 29.5};//{0.2,0.5,1,3,5,10,10.1,20,20.1,29.5};
    std::vector<double> ALP_photon_couplings = {1.0};
    std::vector<RResultMap<float>> resultmaps;
    std::vector<ROOT::RDF::RResultHandle> Nodes;
    
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
    
    int counter = 0;
    
    for (auto& file: input_filenames)
    {
        SchottDataFrame df(MakeRDF(file, 8));
    
        auto newDf = df.Define("EventWeight",
        [](const RVec<float>& ei_event_weights_generator)
        {
            return  ((ei_event_weights_generator[0]) ? 1 / ei_event_weights_generator[0] : 1);

        }, {"ei_event_weights_generator"}).Define("totEventWeight", [](RVec<float> photon_id_eff, RVec<float> photon_iso_eff, RVec<float> photon_trg_eff/*, RVec<float> ei_event_weights_generator*/)
        {
            auto ResizeVal = std::max({photon_id_eff.size(), photon_iso_eff.size(), photon_trg_eff.size()});
            photon_id_eff.resize(ResizeVal,1);
            photon_iso_eff.resize(ResizeVal,1);
            photon_trg_eff.resize(ResizeVal,1);
            
            return photon_id_eff*photon_iso_eff*photon_trg_eff; //Returns a vector of efficiencies

        }, {"photon_id_eff", "photon_iso_eff", "photon_trg_eff",/* "ei_event_weights_generator"*/}).Define("totEventWeightFactor",[](RVec<float>& photon_efficiencies)
        {
            float total = 1.0f;
            
            for (auto i: photon_efficiencies)
            {
                total *= i;
            }
            
            return total; //1 efficiency
            
        }, {"totEventWeight"}).Define("truth_axions", [&](RVec<TruthParticle> truth_particles)
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
        
        auto trigger_selection = newDf
        .Filter([](const RVec<std::string>& trigger_passed_triggers)
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
                 return (!((ep.electron_pt/1e3 > 20) && (abs(ep.electron_eta) < 2.37) &&
                 (!((1.37 < abs(ep.electron_eta)) && (abs(ep.electron_eta) < 1.52)))
                 && (ep.electron_id_medium == 1)));

            }), electrons.end());

            return electrons;

        },{"abstract_electrons"}).Filter(
        [](const RVec<Muon>& muons, RVec<AbstractParticle> electrons)
        {
            return (electrons.size()==2 && muons.empty());

        }, {"muons", "di_electrons"});
        
        auto opp_charge = two_leptons.Filter([](const RVec<AbstractParticle>& electrons)
        {
            return (electrons[0].electron_charge*electrons[1].electron_charge < 0);

        }, {"di_electrons"});
        
        auto leading_pt = opp_charge.Filter([](const RVec<AbstractParticle>& electrons)
        {
            return ((electrons[0].electron_pt > 20e3 && electrons[1].electron_pt > 27e3) || (electrons[1].electron_pt > 20e3 && electrons[0].electron_pt > 27e3));
        }, {"di_electrons"});
        
//        auto delta_R = leading_pt.Filter([] (RVec<Electron>& electrons)
//        {
//            return (DeltaR(electrons[0].Vector(), electrons[1].Vector()) > 0.01);
//        }, {"di_electrons"});
        
        auto same_flavour = leading_pt.Filter([] (const RVec<AbstractParticle>& electrons)
        {
            return true; //abs(electrons[0].electron_pdg_id) == abs(electrons[1].electron_pdg_id) == 11;
        }, {"di_electrons"});
        
        auto dilep = same_flavour.Define("dilep",[] (RVec<AbstractParticle>& electrons)
        {
            return (electrons[0].ElectronVector() + electrons[1].ElectronVector());
        }, {"di_electrons"});
        
        auto mass = dilep.Filter([] (PtEtaPhiEVector& dilep)
        {
            auto mass = dilep.M()/1e3;
            return ((mass >= 81) && (mass <= 101));
        }, {"dilep"});
        
        auto ptCut = mass.Filter([] (PtEtaPhiEVector& dilep)
        {
            auto pT = dilep.Pt()/1e3;
            return pT > 10;
        }, {"dilep"});
        
        auto photonPtDeltaR = ptCut.Define("photons_pass_cut_indices",
        [&](RVec<AbstractParticle>& p) //p = photon
        {
            RVec<int> x; //indices of photons that passed the cuts
            for (auto i = 0; i < p.size(); i++)
            {
                if (not ((abs(p[i].photon_eta) >= 2.37) || (p[i].photon_pt <= 10e3) || (abs(p[i].photon_eta) > 1.37 && abs(p[i].photon_eta) < 1.52)))
                {
                    x.push_back(i);
                }
            }
            return x;
            
        }, {"abstract_photons"}).Define("photonPtDeltaR",
        [&](RVec<AbstractParticle>& photons, RVec<int>& x)
        {
            return Take(photons, x); //Taking only the photons that passed the cuts in each event
            
        }, {"abstract_photons", "photons_pass_cut_indices"}).Define("photonPtDeltaR_photons_pass_cut_indices",
        [&](RVec<AbstractParticle>& reco_photons_matched)
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
            if (chosen_delta_r < 1.5 && pt1 > 10e3 && pt2 > 10e3)
            {
                return x;
            }
            
            x.clear();
            return x;

        }, {"photonPtDeltaR"}).Filter(
        [&](RVec<unsigned long>& indices)
        {
            return (indices.size()==2);
          
        }, {"photonPtDeltaR_photons_pass_cut_indices"}).Define("photonPtDeltaR_TotalEventWeightFactor",
        [](RVec<unsigned long>& photonPtDeltaR_photons_pass_cut_indices, RVec<int>& photons_pass_cut_indices, RVec<float>& photon_efficiencies)
        {
            RVec<float> photonPtDeltaR_photon_efficiencies = Take(Take(photon_efficiencies, photons_pass_cut_indices), photonPtDeltaR_photons_pass_cut_indices);
            float total = 1.0f;
            
            for (auto i: photonPtDeltaR_photon_efficiencies)
            {
                total *= i;
            }
            
            return total;
        }, {"photonPtDeltaR_photons_pass_cut_indices", "photons_pass_cut_indices", "totEventWeight"});

        auto X_window = photonPtDeltaR.Define("X_window_photons_pass_cut_indices",
        [](RVec<AbstractParticle>& reco_photons_matched)
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
        }, {"photonPtDeltaR"}).Filter(
        [&](RVec<unsigned long>& X_window_photons_pass_cut_indices)
        {
            return (X_window_photons_pass_cut_indices.size()==2);
            
        }, {"X_window_photons_pass_cut_indices"})
        .Define("chosen_two",[&](RVec<AbstractParticle>& reco_photons_matched, RVec<unsigned long>& X_window_photons_pass_cut_indices)
        {
            return Take(reco_photons_matched, X_window_photons_pass_cut_indices);
            
        }, {"photonPtDeltaR", "X_window_photons_pass_cut_indices"})
        .Define("X_window_TotalEventWeightFactor",
        [&](RVec<unsigned long>& X_window_photons_pass_cut_indices, RVec<int>& photons_pass_cut_indices, RVec<float>& photon_efficiencies)
        {
            RVec<float> X_window_photon_efficiencies = Take(Take(photon_efficiencies, photons_pass_cut_indices), X_window_photons_pass_cut_indices);
            float total = 1.0f;
            
            for (auto i: X_window_photon_efficiencies)
            {
                total *= i;
            }
            
            return total;
        }, {"X_window_photons_pass_cut_indices", "photons_pass_cut_indices", "totEventWeight"});
                
        auto SR = X_window.Filter(
        [](RVec<AbstractParticle>& photons, RVec<AbstractParticle>& electrons)
        {
            PtEtaPhiEVector photonVec = photons[0].PhotonVector() + photons[1].PhotonVector();
            PtEtaPhiEVector electronVec = electrons[0].ElectronVector() + electrons[1].ElectronVector();
            auto mass = (photonVec+electronVec).M()/1e3;
            return ((mass >= 110) && (mass <= 140));
        },{"chosen_two","di_electrons"});
                
        auto SR_ID = SR.Filter(
        [](RVec<AbstractParticle>& photons)
        {
            return (photons[0].photon_id_loose && photons[1].photon_id_loose);
        },{"chosen_two"});
        
        for (auto& mass_point: massPoints)
        {
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
            
            Nodes.push_back(mass_point_newDf.Count());
            Nodes.push_back(mass_point_trigger_selection.Count());
            Nodes.push_back(mass_point_ptCut.Count());
            Nodes.push_back(mass_point_photonPtDeltaR.Count());
            Nodes.push_back(mass_point_X_window.Count());
            Nodes.push_back(mass_point_SR.Count());
            Nodes.push_back(mass_point_SR_ID.Count());
            
            resultmaps.push_back(VariationsFor(mass_point_newDf.Sum<float>("totEventWeightFactor")));
            resultmaps.push_back(VariationsFor(mass_point_trigger_selection.Sum<float>("totEventWeightFactor")));
            resultmaps.push_back(VariationsFor(mass_point_ptCut.Sum<float>("totEventWeightFactor")));
            resultmaps.push_back(VariationsFor(mass_point_photonPtDeltaR.Sum<float>("photonPtDeltaR_TotalEventWeightFactor")));
            resultmaps.push_back(VariationsFor(mass_point_X_window.Sum<float>("X_window_TotalEventWeightFactor")));
            resultmaps.push_back(VariationsFor(mass_point_SR.Sum<float>("X_window_TotalEventWeightFactor")));
            resultmaps.push_back(VariationsFor(mass_point_SR_ID.Sum<float>("X_window_TotalEventWeightFactor")));

        }
        counter++;
    }
    
    ROOT::RDF::RunGraphs(Nodes); // running all computation nodes concurrently, this will trigger resultmaps as well.
    
//    for (auto& j: resultmaps[0].GetKeys())
//    {
//        std::cout << j << '\n';
//    }
//
//    std::cout << '\n';
//
//    for (auto& j: resultmaps[2].GetKeys())
//    {
//        std::cout << j << '\n';
//    }
//
//    std::cout << '\n';
    

    std::vector<std::string> syst_indices = {"nominal", "photons_and_electrons:EG_RESOLUTION_ALL__1down", "photons_and_electrons:EG_SCALE_ALL__1down", "photon_iso_eff:PH_EFF_ISO_Uncertainty__1down", "photon_id_eff:PH_EFF_ID_Uncertainty__1down", "photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down", "photons_and_electrons:EG_RESOLUTION_ALL__1up", "photons_and_electrons:EG_SCALE_ALL__1up", "photon_iso_eff:PH_EFF_ISO_Uncertainty__1up", "photon_id_eff:PH_EFF_ID_Uncertainty__1up", "photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up"};

    TLatex Tl;
    TCanvas* c1;
    TH2D* histo2d;

    //for each systematic
    for (auto& syst_index: syst_indices)
    {
        //create a plot
        c1 = new TCanvas();
        std::string title = syst_index + ";m_{a} [GeV];C_{a#gamma#gamma};Efficiency";
        std::string name = std::string("C_{a#gamma#gamma}") + syst_index;
        histo2d = new TH2D(name.c_str(), title.c_str(), 295, 0.2, 29.5, 5, 0, 2);
        histo2d->GetXaxis()->SetTitleSize(0.04);
        histo2d->GetYaxis()->SetTitleSize(0.04);
        histo2d->GetZaxis()->SetTitleSize(0.04);
        
        histo2d->GetYaxis()->SetTitleOffset(0.7);
        histo2d->GetZaxis()->SetTitleOffset(0.6);
        
        histo2d->GetYaxis()->CenterTitle(true);
        histo2d->GetZaxis()->SetTitle("Efficiency");

        //i is indexing each mass point, and j is indexing the
        //group of results corresponding to mass point i
        for (int i=0, j=0; (i < massPoints.size()); i++, j+=7)
        {
            int bin_number = histo2d->FindFixBin(massPoints[i], ALP_photon_couplings[0]);
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
            
            {
                std::cout << syst_index << ": " << massPoints[i] << " GeV, "
                << efficiency << '\n';
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
        title = std::string("C_{a#gamma#gamma}") + std::string(syst_index) + ".png";
        c1->SaveAs(title.c_str());
    }
}

void Coupling_and_Systematics_merged()
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

    
    
    std::vector<std::vector<std::string>> input_filenames = {
        //Signal
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root"},
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4.root"},
        //Displaced Signal
//        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_FewMassPoints.root"},
        {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/ZaSignal_FewMassPoints_PhotonSFs.root"},
    };
    
    std::vector<float> massPoints = {0.2, 0.5, 1, 3, 5, 10, 20, 29.5};//{0.2,0.5,1,3,5,10,10.1,20,20.1,29.5};
    std::vector<double> ALP_photon_couplings = {1.0};
    std::vector<RResultMap<float>> resultmaps;
    std::vector<ROOT::RDF::RResultHandle> Nodes;
    
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
    
    int counter = 0;
    
    for (auto& file: input_filenames)
    {
        SchottDataFrame df(MakeRDF(file, 8));
    
        auto newDf = df.Define("EventWeight",
        [](const RVec<float>& ei_event_weights_generator)
        {
            return  ((ei_event_weights_generator[0]) ? 1 / ei_event_weights_generator[0] : 1);

        }, {"ei_event_weights_generator"}).Define("totEventWeightVec", [](RVec<float> photon_id_eff, RVec<float> photon_iso_eff, RVec<float> photon_trg_eff/*, RVec<float> ei_event_weights_generator*/)
        {
            auto ResizeVal = std::max({photon_id_eff.size(), photon_iso_eff.size(), photon_trg_eff.size()});
            photon_id_eff.resize(ResizeVal,1);
            photon_iso_eff.resize(ResizeVal,1);
            photon_trg_eff.resize(ResizeVal,1);
            
            return photon_id_eff*photon_iso_eff*photon_trg_eff; //Returns a vector of efficiencies

        }, {"photon_id_eff", "photon_iso_eff", "photon_trg_eff",/* "ei_event_weights_generator"*/}).Define("totEventWeightFactor",[](RVec<float>& photon_efficiencies)
        {
            float total = 1.0f;
            
            for (auto i: photon_efficiencies)
            {
                total *= i;
            }
            
            return total; //1 efficiency
            
        }, {"totEventWeightVec"}).Define("truth_axions", [&](RVec<TruthParticle> truth_particles)
        {
            truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
            [](TruthParticle& x)
            {
                return (abs(x.mc_pdg_id) != 36); //axions are pdg id 36 in these samples

            }), truth_particles.end());

            return truth_particles;

        }, {"truth_particles"}).Define("axion_masses", [&](RVec<TruthParticle>& truth_axions)
        {
            return truth_axions[0].mc_mass/1e3f;

        }, {"truth_axions"});
        
        auto trigger_selection = newDf
        .Filter([](const RVec<std::string>& trigger_passed_triggers)
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
                 return (!((ep.electron_pt/1e3 > 20) && (abs(ep.electron_eta) < 2.37) &&
                 (!((1.37 < abs(ep.electron_eta)) && (abs(ep.electron_eta) < 1.52)))
                 && (ep.electron_id_medium == 1)));

            }), electrons.end());

            return electrons;

        },{"abstract_electrons"}).Filter(
        [](const RVec<Muon>& muons, RVec<AbstractParticle> electrons)
        {
            return (electrons.size()==2 && muons.empty());

        }, {"muons", "di_electrons"});
        
        auto opp_charge = two_leptons.Filter([](const RVec<AbstractParticle>& electrons)
        {
            return (electrons[0].electron_charge*electrons[1].electron_charge < 0);

        }, {"di_electrons"});
        
        auto leading_pt = opp_charge.Filter([](const RVec<AbstractParticle>& electrons)
        {
            return ((electrons[0].electron_pt > 20e3 && electrons[1].electron_pt > 27e3) || (electrons[1].electron_pt > 20e3 && electrons[0].electron_pt > 27e3));
        }, {"di_electrons"});
        
//        auto delta_R = leading_pt.Filter([] (RVec<Electron>& electrons)
//        {
//            return (DeltaR(electrons[0].Vector(), electrons[1].Vector()) > 0.01);
//        }, {"di_electrons"});
        
        auto same_flavour = leading_pt.Filter([] (const RVec<AbstractParticle>& electrons)
        {
            return true; //abs(electrons[0].electron_pdg_id) == abs(electrons[1].electron_pdg_id) == 11;
        }, {"di_electrons"});
        
        auto dilep = same_flavour.Define("dilep",[] (RVec<AbstractParticle>& electrons)
        {
            return (electrons[0].ElectronVector() + electrons[1].ElectronVector());
        }, {"di_electrons"});
        
        auto mass = dilep.Filter([] (PtEtaPhiEVector& dilep)
        {
            auto mass = dilep.M()/1e3;
            return ((mass >= 81) && (mass <= 101));
        }, {"dilep"});
        
        auto ptCut = mass.Filter([] (PtEtaPhiEVector& dilep)
        {
            auto pT = dilep.Pt()/1e3;
            return pT > 10;
        }, {"dilep"});
        
        //end pre-selection -----------
        
        //photon acceptance and id_loose cuts
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

        }, {"photons_pass_cuts", "merged_photon_index"})
        //mpi = merged_photon_index
        .Define("totEventWeight", [](RVec<int>& x, int mpi, RVec<float> photon_id_eff, RVec<float> photon_iso_eff, RVec<float> photon_trg_eff/*, RVec<float> ei_event_weights_generator*/)
        {
            //   ||  resize 3 vectors jic they don't already have the same size  ||
            //   \/                                                              \/
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
        .Define("reconstructed_mass",[&](const RVec<AbstractParticle>& di_electrons, const AbstractParticle& merged_photon)
        {
            auto four_momentum = di_electrons[0].ElectronVector() + di_electrons[1].ElectronVector();

            return (four_momentum + merged_photon.PhotonVector()).M()/1e3;

        }, {"di_electrons", "merged_photon"});
        
        auto pSB = merged_reco_photons_matched.Filter(
        [](const double reconstructed_mass)
        {
            return (reconstructed_mass < 110) || (reconstructed_mass > 130);
        }, {"reconstructed_mass"});

        auto pSR = merged_reco_photons_matched.Filter(
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
        counter++;
    }
    
    ROOT::RDF::RunGraphs(Nodes); // running all computation nodes concurrently, this will trigger resultmaps as well.
    
    std::vector<std::string> syst_indices = {"nominal", "photons_and_electrons:EG_RESOLUTION_ALL__1down", "photons_and_electrons:EG_SCALE_ALL__1down", "photon_iso_eff:PH_EFF_ISO_Uncertainty__1down", "photon_id_eff:PH_EFF_ID_Uncertainty__1down", "photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1down", "photons_and_electrons:EG_RESOLUTION_ALL__1up", "photons_and_electrons:EG_SCALE_ALL__1up", "photon_iso_eff:PH_EFF_ISO_Uncertainty__1up", "photon_id_eff:PH_EFF_ID_Uncertainty__1up", "photon_trg_eff:PH_EFF_TRIGGER_Uncertainty__1up"};

    TLatex Tl;
    TCanvas* c1;
    TH2D* histo2d;

    //for each systematic
    for (auto& syst_index: syst_indices)
    {
        //create a plot
        c1 = new TCanvas();
        std::string title = syst_index + ";m_{a} [GeV];C_{a#gamma#gamma};Efficiency";
        std::string name = std::string("C_{a#gamma#gamma}") + syst_index;
        histo2d = new TH2D(name.c_str(), title.c_str(), 295, 0.2, 29.5, 5, 0, 2);
        histo2d->GetXaxis()->SetTitleSize(0.04);
        histo2d->GetYaxis()->SetTitleSize(0.04);
        histo2d->GetZaxis()->SetTitleSize(0.04);
        
        histo2d->GetYaxis()->SetTitleOffset(0.7);
        histo2d->GetZaxis()->SetTitleOffset(0.6);
        
        histo2d->GetYaxis()->CenterTitle(true);
        histo2d->GetZaxis()->SetTitle("Efficiency");

        //i is indexing each mass point, and j is indexing the
        //group of results corresponding to mass point i
        for (int i=0, j=0; (i < massPoints.size()); i++, j+=6)
        {
            int bin_number = histo2d->FindFixBin(massPoints[i], ALP_photon_couplings[0]);
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
            
            {
                std::cout << syst_index << ": " << massPoints[i] << " GeV, "
                << efficiency << '\n';
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
        title = std::string("C_{a#gamma#gamma}") + std::string(syst_index) + "_merged.pdf";
        c1->SaveAs(title.c_str());
    }
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
//    Coupling_and_Systematics();
    Coupling_and_Systematics_merged();
    
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


