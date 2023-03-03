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
    
    std::vector<std::string> cutFlows;
    cutFlows.reserve(input_filenames.size());
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
//        "PH_EFF_ID_Uncertainty__1up",
//        "PH_EFF_TRIGGER_Uncertainty__1up",
        "EG_SCALE_ALL__1down",
        "EG_SCALE_ALL__1up",
//        "PH_EFF_ID_Uncertainty__1down",
//        "PH_EFF_TRIGGER_Uncertainty__1down",
//        "PH_EFF_ISO_Uncertainty__1up",
        "EG_RESOLUTION_ALL__1up",
        "EG_RESOLUTION_ALL__1down",
    };
    
    std::array<double,9> JetNumeratorSFs = {((139e15)*(1.9828e-9)*(0.821204)),((139e15)*(110.64e-12)*(0.69275)),((139e15)*(40.645e-12)*(0.615906)),((139e15)*(1.9817e-9)*(0.1136684)),((139e15)*(110.47e-12)*(0.1912956)),((139e15)*(40.674e-12)*(0.2326772)),((139e15)*(1.9819e-9)*(0.0656969)),((139e15)*(110.53e-12)*(0.1158741)),((139e15)*(40.68e-12)*(0.1535215))}; //numerators for jet bkg
    std::array<double,3> SFs = {((139e15)*(.871e-12)),((139e15)*(.199e-12)), ((139e15)*(.0345e-15))}; //numerators for Z-gamma bkg
    
    constexpr std::array<const char*,15> Samples = {R"--(pty2\_9\_17)--", R"--(pty\_17\_myy\_0\_80)--", R"--(pty\_17\_myy\_80)--", R"--(Signal $m_{\text{A}}$ = 1 GeV)--", R"--(Signal $m_{\text{A}}$ = 5 GeV)--", "Data", R"--(Zee\_lightJet\_0-70)--", R"--(Zee\_lightJet\_70-140)--", R"--(Zee\_lightJet\_140-280)--", R"--(Zee\_cJet\_0-70)--", R"--(Zee\_cJet\_70-140)--", R"--(Zee\_cJet\_140-280)--", R"--(Zee\_bJet\_0-70)--", R"--(Zee\_bJet\_70-140)--", R"--(Zee\_bJet\_140-280)--"};
    
    std::ostringstream os;
    os << R"--(\section*{Table 3})--" << '\n';
    os << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
    os << R"--(\begin{tabular}{|c|c|c|c|c|c|c|c|})--" << '\n';
    os << R"--(\hline)--" << '\n';
    os << R"--(Sample & Before Preselection & 2 leptons
          & Opposite Charge & $p_{T}^{\text{leading}} > 27$ GeV, \; $p_{T}^{\text{sub-leading}} > 20$ GeV & Same flavour & dilep mass cut & dilep $p_{T}$ cut \\ \hline )--" << '\n';
    
    double beforePreselecZGamma = 0, twoLeptonsZGamma = 0, oppChargeZGamma = 0, leadingPtZGamma = 0, deltaRZGamma = 0, MassZGamma = 0, ptCutZGamma = 0;

    double beforePreselecZJets = 0, twoLeptonsZJets = 0, oppChargeZJets = 0, leadingPtZJets = 0, deltaRZJets = 0, MassZJets = 0, ptCutZJets = 0;
    
    std::vector<ROOT::RDF::RResultHandle> Nodes;
    std::vector<ROOT::RDF::RResultHandle> tempNodes;
    
    std::vector<RResultMap<ULong64_t>> resultmaps;
    
    for (auto& file: input_filenames)
    {
        SchottDataFrame df(MakeRDF(file, 8));
        
//        df.Describe().Print();
//        exit(1);
        
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
        
        resultmaps.push_back(VariationsFor(df.Count()));
        resultmaps.push_back(VariationsFor(two_leptons.Count()));
        resultmaps.push_back(VariationsFor(opp_charge.Count()));
        resultmaps.push_back(VariationsFor(leading_pt.Count()));
        resultmaps.push_back(VariationsFor(same_flavour.Count()));
        resultmaps.push_back(VariationsFor(mass.Count()));
        resultmaps.push_back(VariationsFor(pt_cut.Count()));
        
        tempNodes.push_back(df.Take<RVec<std::vector<std::string>>, RVec<RVec<std::vector<std::string>>>>("electron_syst_name"));
        tempNodes.push_back(df.Take<RVec<std::vector<std::string>>, RVec<RVec<std::vector<std::string>>>>("photon_syst_name"));
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
    
    std::unordered_set<std::string> uniqueSystematics;
    std::unordered_set<std::string> ZGammaSystematics;
    std::unordered_set<std::string> SignalSystematics;
    std::unordered_set<std::string> DataSystematics;
    std::unordered_set<std::string> ZJetsSystematics;
//
    for (auto& i: tempNodes)
    {
        for (auto& j: *i.GetResultPtr<RVec<RVec<std::vector<std::string>>>>())
        {
            for (auto& k: j)
            {
                for (auto& l: k)
                {
                    uniqueSystematics.insert(l);
                }
            }
        }
    }
    
    int counter = 0;
    for (auto& i: tempNodes)
    {
        for (auto& j: (*i.GetResultPtr<RVec<RVec<std::vector<std::string>>>>())[0])
        {
            for (auto& k: j)
            {
                if (counter >= 0 && counter <= 5) //Z-gamma
                {
                    ZGammaSystematics.insert(k);
                }
            
                else if (counter >= 6 && counter <= 9) //Signal
                {
                    SignalSystematics.insert(k);
                }
                
                else if (counter == 10) //Data
                {
                    DataSystematics.insert(k);
                }
                
                else
                {
                    ZJetsSystematics.insert(k);
                }
            }
        }
        counter++;
    }
//
    std::cout << "ZGammaSystematics\n=================\n";
    for (auto& i: ZGammaSystematics)
    {
        std::cout << i << '\n';
    }
    std::cout << "\n\n\n\n";
    
    std::cout << "SignalSystematics\n=================\n";
    for (auto& i: SignalSystematics)
    {
        std::cout << i << '\n';
    }
    std::cout << "\n\n\n\n";
    
    std::cout << "DataSystematics\n===============\n";
    for (auto& i: DataSystematics)
    {
        std::cout << i << '\n';
    }
    std::cout << "\n\n\n\n";
    
    std::cout << "ZJetsSystematics\n================\n";
    for (auto& i: ZJetsSystematics)
    {
        std::cout << i << '\n';
    }
    std::cout << "\n\n\n\n";
    
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
    
    for (int i=0, j=0; (i<15 && j <= 98); i++, j+=7)
    {
        os << Samples[i];
        if (i >= 0 && i <= 2) //Z-gamma
        {
            os  << " & " << *Nodes[j].GetResultPtr<ULong64_t>() * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+1].GetResultPtr<ULong64_t>() * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+2].GetResultPtr<ULong64_t>() * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+3].GetResultPtr<ULong64_t>() * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+4].GetResultPtr<ULong64_t>() * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+5].GetResultPtr<ULong64_t>() * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+6].GetResultPtr<ULong64_t>() * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--( \\ \hline )--" << '\n';
        }
        
        else if (i >= 6) //Z-jets
        {
            os  << " & " << *Nodes[j].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-6] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+1].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-6] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+2].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-6] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+3].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-6] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+4].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-6] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+5].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-6] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+6].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-6] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--( \\ \hline )--" << '\n';
        }
        else
        {
            os  << " & " << *Nodes[j].GetResultPtr<ULong64_t>()
                << " & " << *Nodes[j+1].GetResultPtr<ULong64_t>()
                << " & " << *Nodes[j+2].GetResultPtr<ULong64_t>()
                << " & " << *Nodes[j+3].GetResultPtr<ULong64_t>()
                << " & " << *Nodes[j+4].GetResultPtr<ULong64_t>()
                << " & " << *Nodes[j+5].GetResultPtr<ULong64_t>()
                << " & " << *Nodes[j+6].GetResultPtr<ULong64_t>()
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
    }

    os << R"--(Total $Z\gamma$ & )--" << beforePreselecZGamma << " & " << twoLeptonsZGamma
    << " & " << oppChargeZGamma << " & " << leadingPtZGamma << " & "
    << deltaRZGamma << " & " << MassZGamma << " & " << ptCutZGamma
    << R"--( \\ \hline )--" << '\n';
    
    os << R"--(Total $Z$ jets & )--" << beforePreselecZJets << " & " << twoLeptonsZJets
    << " & " << oppChargeZJets << " & " << leadingPtZJets << " & "
    << deltaRZJets << " & " << MassZJets << " & " << ptCutZJets
    << R"--( \\ \hline )--" << '\n';
    
    os << R"--(Total Bkg & )--" << beforePreselecZJets+beforePreselecZGamma << " & " << twoLeptonsZJets+twoLeptonsZGamma
    << " & " << oppChargeZJets+oppChargeZGamma << " & " << leadingPtZJets+leadingPtZGamma << " & "
    << deltaRZJets+deltaRZGamma << " & " << MassZJets+MassZGamma << " & " << ptCutZJets+ptCutZGamma
    << R"--( \\ \hline )--" << '\n';
    
    os << R"--(\end{tabular}})--" << '\n';
    cutFlows.push_back(os.str());
    for (auto& i: cutFlows)
    {
        std::cout << i << "\n\n";
    }
}
/*
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
    double totalEventsZJets = 0, resolvedEventsZJets = 0, SREventsZJets = 0, SBEventsZJets = 0;
    std::vector<std::string> cutFlows;
    cutFlows.reserve(input_filenames.size());
    
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
 
    for (int i=0, j=0; (i<15 && j <= 56); i++, j+=4)
    {
        os << Samples[i];
        if (i >= 0 && i <= 2)
        {
            os  << " & " << *Nodes[j].GetResultPtr<ULong64_t>() * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+1].GetResultPtr<ULong64_t>() * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+2].GetResultPtr<ULong64_t>() * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+3].GetResultPtr<ULong64_t>() * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--( \\ \hline )--" << '\n';
        }
        
        else if (i >= 6)
        {
            os  << " & " << *Nodes[j].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-6] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+1].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-6] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+2].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-6] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+3].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-6] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--( \\ \hline )--" << '\n';
        }
        
        else
        {
            os  << " & " << *Nodes[j].GetResultPtr<ULong64_t>()
                << " & " << *Nodes[j+1].GetResultPtr<ULong64_t>()
                << " & " << *Nodes[j+2].GetResultPtr<ULong64_t>()
                << " & " << *Nodes[j+3].GetResultPtr<ULong64_t>()
                << R"--( \\ \hline )--" << '\n';
        }
    }
    
    for (int i = 0, j = 0; (i <= 8 && j <= 2); i += 4, j++)
    {
        totalEventsZgamma += *Nodes[i].GetResultPtr<ULong64_t>() * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
        resolvedEventsZgamma += *Nodes[i+1].GetResultPtr<ULong64_t>() * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
        SBEventsZgamma += *Nodes[i+2].GetResultPtr<ULong64_t>() * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
        SREventsZgamma += *Nodes[i+3].GetResultPtr<ULong64_t>() * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
    }
    
    for (int i = 24, j = 0; (i <= 56 && j <= 8); i += 4, j++)
    {
        totalEventsZJets += *Nodes[i].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
        resolvedEventsZJets += *Nodes[i+1].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
        SBEventsZJets += *Nodes[i+2].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
        SREventsZJets += *Nodes[i+3].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
    }
    
    os << R"--(Total $Z\gamma$ & )--" << totalEventsZgamma << " & " << resolvedEventsZgamma
    << " & " << SBEventsZgamma << " & " << SREventsZgamma
    << R"--( \\ \hline )--" << '\n';
    
    os << R"--(Total $Z$ jets & )--" << totalEventsZJets << " & " << resolvedEventsZJets
        << " & " << SBEventsZJets << " & " << SREventsZJets
        << R"--( \\ \hline )--" << '\n';
    
    os << R"--(Total Bkg & )--" << totalEventsZJets+totalEventsZgamma << " & " << resolvedEventsZJets+resolvedEventsZgamma
        << " & " << SBEventsZJets+SBEventsZgamma << " & " << SREventsZJets+SREventsZgamma
        << R"--( \\ \hline )--" << '\n';
    
    os << R"--(\end{tabular}})--" << '\n';
    cutFlows.push_back(os.str());
    for (auto& i: cutFlows)
    {
        std::cout << i << "\n\n";
    }
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
    
    std::vector<std::string> cutFlows;
    cutFlows.reserve(input_filenames.size());
    double totalEventsZgamma = 0, passPreselectionZgamma = 0, photonPtDeltaRCountZgamma = 0, xWindowZgamma = 0,
    srCountZgamma = 0, srIDCountZgamma = 0;
    double totalEventsZjets = 0, passPreselectionZjets = 0, photonPtDeltaRCountZjets = 0, xWindowZjets = 0,
        srCountZjets = 0, srIDCountZjets = 0;
    std::vector<ROOT::RDF::RResultHandle> Nodes;

    std::ostringstream os;
    os << R"--(\section*{Table 11})--" << '\n';
    os << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
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
    
    auto TruthMatching = [](SchottDataFrame df)
    {
        
    };
    
    for (auto& file: input_filenames)
    {
        SchottDataFrame df(MakeRDF(file, 8));
        
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
    
    for (int i=0, j=0; (i<15 && j <= 84); i++, j+=6)
    {
        os << Samples[i];
        
        if (i >= 0 && i <= 2)
        {
            os  << " & " << *Nodes[j].GetResultPtr<ULong64_t>() * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+1].GetResultPtr<ULong64_t>() * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+2].GetResultPtr<ULong64_t>() * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+3].GetResultPtr<ULong64_t>() * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+4].GetResultPtr<ULong64_t>() * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+5].GetResultPtr<ULong64_t>() * (SFs[i] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--( \\ \hline )--" << '\n';
        }
        else if (i >= 6)
        {
            os  << " & " << *Nodes[j].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-6] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+1].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-6] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+2].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-6] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+3].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-6] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+4].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-6] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+5].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-6] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--( \\ \hline )--" << '\n';
        }
        else
        {
            os  << " & " << *Nodes[j].GetResultPtr<ULong64_t>()
                << " & " << *Nodes[j+1].GetResultPtr<ULong64_t>()
                << " & " << *Nodes[j+2].GetResultPtr<ULong64_t>()
                << " & " << *Nodes[j+3].GetResultPtr<ULong64_t>()
                << " & " << *Nodes[j+4].GetResultPtr<ULong64_t>()
                << " & " << *Nodes[j+5].GetResultPtr<ULong64_t>()
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
            }
        }
    }
    
    std::cout << R"--(\end{tabular}})--" << "\n\n\n";
    
    for (int i = 0, j = 0; (i <= 12 && j <= 2); i += 6, j++)
    {
        totalEventsZgamma += *Nodes[i].GetResultPtr<ULong64_t>() * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
        passPreselectionZgamma += *Nodes[i+1].GetResultPtr<ULong64_t>() * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
        photonPtDeltaRCountZgamma += *Nodes[i+2].GetResultPtr<ULong64_t>() * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
        xWindowZgamma += *Nodes[i+3].GetResultPtr<ULong64_t>() * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
        srCountZgamma += *Nodes[i+4].GetResultPtr<ULong64_t>() * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
        srIDCountZgamma += *Nodes[i+5].GetResultPtr<ULong64_t>() * (SFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
    }
    
    for (int i = 36, j = 0; (i <= 84 && j <= 8); i += 6, j++)
    {
        totalEventsZjets += *Nodes[i].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
        passPreselectionZjets += *Nodes[i+1].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
        photonPtDeltaRCountZjets += *Nodes[i+2].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
        xWindowZjets += *Nodes[i+3].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
        srCountZjets += *Nodes[i+4].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
        srIDCountZjets += *Nodes[i+5].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j] / *Nodes[i].GetResultPtr<ULong64_t>());
    }
    
    os << R"--(Total $Z\gamma$ & )--" << totalEventsZgamma
        << " & " << passPreselectionZgamma
        << " & " << photonPtDeltaRCountZgamma
        << " & " << xWindowZgamma
        << " & " << srCountZgamma
        << " & " << srIDCountZgamma
        << R"--( \\ \hline )--" << '\n';
    
    os << R"--(Total $Z$ jets & )--" << totalEventsZjets
        << " & " << passPreselectionZjets
        << " & " << photonPtDeltaRCountZjets
        << " & " << xWindowZjets
        << " & " << srCountZjets
        << " & " << srIDCountZjets
        << R"--( \\ \hline )--" << '\n';
    
    os << R"--(Total Bkg & )--" << totalEventsZjets+totalEventsZgamma
        << " & " << passPreselectionZjets+passPreselectionZgamma
        << " & " << photonPtDeltaRCountZjets+photonPtDeltaRCountZgamma
        << " & " << xWindowZjets+xWindowZgamma
        << " & " << srCountZjets+srCountZgamma
        << " & " << srIDCountZjets+srIDCountZgamma
        << R"--( \\ \hline )--" << '\n';
    
    os << R"--(\end{tabular}})--" << '\n';
    cutFlows.push_back(os.str());
    for (auto& i: cutFlows)
    {
        std::cout << i << "\n\n";
    }
}
*/
void CutFlow()
{
    auto start_time = Clock::now();
    
    Table3();
//    Table8();
//    Table11();
    
    auto end_time = Clock::now();
    std::cout << "Time difference: "
       << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count()/1e9 << " seconds" << std::endl;
    
}


int main()
{
    CutFlow();
}
