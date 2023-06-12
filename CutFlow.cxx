#include <iostream>
#include <vector>
#include <chrono>
#include <algorithm>
#include <string>
#include <sstream>
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
#include "ROOT/RDFHelpers.hxx"

#include "RDFObjects.h"
#include "MakeRDF.h"
#include "RDFevent.h"

using namespace ROOT::VecOps; // RVec, Combinations
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

void Table3()
{
    std::vector<std::vector<std::string>> input_filenames = {
        //Z gamma background
        {"/home/common/Za/NTuples/Background/user.kschmied.364860.eegammagamma_pty2_9_17.deriv.DAOD_STDM3.e7057_s3126_r10724_p4092_LGNTuple.root/user.kschmied.31617070._000001.LGNTuple.root"},
        {"/home/common/Za/NTuples/Background/user.kschmied.364861.eegammagamma_pty_17_myy_0_80.deriv.DAOD_HIGG1D2.e7057_s3126_r10724_p4062_LGNTuple.root/user.kschmied.31617064._000001.LGNTuple.root"},
        {"/home/common/Za/NTuples/Background/user.kschmied.364862.eegammagamma_pty_17_myy_80.deriv.DAOD_HIGG1D2.e7057_s3126_r10724_p4204_LGNTuple.root/user.kschmied.31660711._000001.LGNTuple.root", "/home/common/Za/NTuples/Background/user.kschmied.364862.eegammagamma_pty_17_myy_80.deriv.DAOD_HIGG1D2.e7057_s3126_r10724_p4204_LGNTuple.root/user.kschmied.31660711._000002.LGNTuple.root"},
        //Signal
        {"/home/common/Za/NTuples/Signal/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root"}, //1 GeV
        {"/home/common/Za/NTuples/Ntuple_MC_Za_m5p0_v4.root"}, //5 GeV
        {"/home/common/Za/NTuples/Signal/mc16_13TeV.600751.PhPy8EG_AZNLO_ggH125_mA2p0_v1.root"}, // 2 GeV
        {"/home/common/Za/NTuples/Signal/mc16_13TeV.600752.PhPy8EG_AZNLO_ggH125_mA3p0_v1.root"}, // 3 GeV
        {"/home/common/Za/NTuples/Signal/mc16_13TeV.600756.PhPy8EG_AZNLO_ggH125_mA9p0_v1.root"}, // 9 GeV
        //Data
        {"/home/common/Za/NTuples/Ntuple_data_test.root"},
        //Jets
        {
            "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000001.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000002.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000003.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000004.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000005.LGNTuple.root"
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
            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000008.LGNTuple.root"
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
    
    std::array<double,9> JetNumeratorSFs = {((139e15)*(1.9828e-9)*(0.821204)),((139e15)*(110.64e-12)*(0.69275)),((139e15)*(40.645e-12)*(0.615906)),((139e15)*(1.9817e-9)*(0.1136684)),((139e15)*(110.47e-12)*(0.1912956)),((139e15)*(40.674e-12)*(0.2326772)),((139e15)*(1.9819e-9)*(0.0656969)),((139e15)*(110.53e-12)*(0.1158741)),((139e15)*(40.68e-12)*(0.1535215))}; //numerators for jet bkg
    std::array<double,3> SFs = {((139e15)*(.871e-12)),((139e15)*(.199e-12)), ((139e15)*(.0345e-15))}; //numerators for Z-gamma bkg

    constexpr std::array<const char*,18> Samples = {R"--(pty2\_9\_17)--", R"--(pty\_17\_myy\_0\_80)--", R"--(pty\_17\_myy\_80)--", R"--(Signal $m_{\text{A}}$ = 1 GeV)--", R"--(Signal $m_{\text{A}}$ = 5 GeV)--", R"--(Signal $m_{\text{A}}$ = 2 GeV)--", R"--(Signal $m_{\text{A}}$ = 3 GeV)--", R"--(Signal $m_{\text{A}}$ = 9 GeV)--", "Data", R"--(Zee\_lightJet\_0-70)--", R"--(Zee\_lightJet\_70-140)--", R"--(Zee\_lightJet\_140-280)--", R"--(Zee\_cJet\_0-70)--", R"--(Zee\_cJet\_70-140)--", R"--(Zee\_cJet\_140-280)--", R"--(Zee\_bJet\_0-70)--", R"--(Zee\_bJet\_70-140)--", R"--(Zee\_bJet\_140-280)--"};

    std::ostringstream os;
    os << R"--(\section*{Table 3})--" << '\n';
    os << R"--(\hspace{-3cm}\scalebox{0.45}{)--" << '\n';
    os << R"--(\setlength\extrarowheight{2pt})--" << '\n';
    os << R"--(\renewcommand{\arraystretch}{1.5})--" << '\n';
    os << R"--(\begin{tabular}{|c|c|c|c|c|c|c|c|})--" << '\n';
    os << R"--(\hline)--" << '\n';
    os << R"--(\textbf{Sample} & \textbf{Before Preselection} & \textbf{2 leptons}
          & \textbf{Opposite Charge} & $\pmb{p_{T}}^{\textbf{leading}}\pmb{>}\text{27}$ \textbf{GeV}, \; $\pmb{p_{T}}^{\textbf{sub-leading}}\pmb{>}\textbf{20}$ \textbf{GeV} & \textbf{Same flavour} & \textbf{dilep mass cut} & \textbf{dilep} $\pmb{p_{T}}$ \textbf{cut} \\ \hline )--" << '\n';

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
                return (!((ep.electron_pt/1e3 > 20) && (std::abs(ep.electron_eta) < 2.37) &&
                                          (!((1.37 < std::abs(ep.electron_eta)) && (std::abs(ep.electron_eta) < 1.52)))
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
            return true; //std::abs(electrons[0].electron_pdg_id) == std::abs(electrons[1].electron_pdg_id) == 11;
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

//    0     1     2     3     4     5     6       Z-gamma
//    7     8     9     10    11    12    13      Z-gamma
//    14    15    16    17    18    19    20      Z-gamma
//    21    22    23    24    25    26    27      ma1
//    28    29    30    31    32    33    34      ma5
//    35    36    37    38    39    40    41      ma2
//    42    43    44    45    46    47    48      ma3
//    49    50    51    52    53    54    55      ma9
//    56    57    58    59    60    61    62      data
//    63    64    65    66    67    68    69      Z-jets
//    70    71    72    73    74    75    76      Z-jets
//    77    78    79    80    81    82    83      Z-jets
//    84    85    86    87    88    89    90      Z-jets
//    91    92    93    94    95    96    97      Z-jets
//    98    99    100   101   102   103   104     Z-jets
//    105   106   107   108   109   110   111     Z-jets
//    112   113   114   115   116   117   118     Z-jets
//    119   120   121   122   123   124   125     Z-jets

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

    for (int i=0, j=0; (i<18 && j <= 119); i++, j+=7)
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

        else if (i >= 9) //Z-jets
        {
            os  << " & " << *Nodes[j].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+1].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+1].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+2].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+2].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+3].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+3].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+4].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+4].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+5].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+5].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+6].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+6].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
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

    for (int i = 63, j = 0; (i <= 119 && j <= 8); i += 7, j++)
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

    os << R"--(Total $Z\gamma\gamma$ & )--"
    << beforePreselecZGamma << R"--($\, \pm \,$)--" << sqrt(beforePreselecZGammaStatUnc) << " & "
    << twoLeptonsZGamma << R"--($\, \pm \,$)--" << sqrt(twoLeptonsZGammaStatUnc) << " & "
    << oppChargeZGamma << R"--($\, \pm \,$)--" << sqrt(oppChargeZGammaStatUnc) << " & "
    << leadingPtZGamma << R"--($\, \pm \,$)--" << sqrt(leadingPtZGammaStatUnc) << " & "
    << deltaRZGamma << R"--($\, \pm \,$)--" << sqrt(deltaRZGammaStatUnc) << " & "
    << MassZGamma << R"--($\, \pm \,$)--" << sqrt(MassZGammaStatUnc) << " & "
    << ptCutZGamma << R"--($\, \pm \,$)--" << sqrt(ptCutZGammaStatUnc)
    << R"--( \\ \hline )--" << '\n';

    os << R"--(Total $Z$+jets & )--"
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
    std::ofstream out("Table3.txt");
    out << os.str() << '\n';
    out.close();
}

void Table8()
{
    std::vector<std::vector<std::string>> input_filenames = {
        //Z gamma background
        {"/home/common/Za/NTuples/Background/user.kschmied.364860.eegammagamma_pty2_9_17.deriv.DAOD_STDM3.e7057_s3126_r10724_p4092_LGNTuple.root/user.kschmied.31617070._000001.LGNTuple.root"},
        {"/home/common/Za/NTuples/Background/user.kschmied.364861.eegammagamma_pty_17_myy_0_80.deriv.DAOD_HIGG1D2.e7057_s3126_r10724_p4062_LGNTuple.root/user.kschmied.31617064._000001.LGNTuple.root"},
        {"/home/common/Za/NTuples/Background/user.kschmied.364862.eegammagamma_pty_17_myy_80.deriv.DAOD_HIGG1D2.e7057_s3126_r10724_p4204_LGNTuple.root/user.kschmied.31660711._000001.LGNTuple.root", "/home/common/Za/NTuples/Background/user.kschmied.364862.eegammagamma_pty_17_myy_80.deriv.DAOD_HIGG1D2.e7057_s3126_r10724_p4204_LGNTuple.root/user.kschmied.31660711._000002.LGNTuple.root"},
        //Signal
        {"/home/common/Za/NTuples/Signal/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root"}, //1 GeV
        {"/home/common/Za/NTuples/Ntuple_MC_Za_m5p0_v4.root"}, //5 GeV
        {"/home/common/Za/NTuples/Signal/mc16_13TeV.600751.PhPy8EG_AZNLO_ggH125_mA2p0_v1.root"}, // 2 GeV
        {"/home/common/Za/NTuples/Signal/mc16_13TeV.600752.PhPy8EG_AZNLO_ggH125_mA3p0_v1.root"}, // 3 GeV
        {"/home/common/Za/NTuples/Signal/mc16_13TeV.600756.PhPy8EG_AZNLO_ggH125_mA9p0_v1.root"}, // 9 GeV
        //Data
        {"/home/common/Za/NTuples/Ntuple_data_test.root"},
        //Jets
        {
            "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000001.LGNTuple.root", "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000002.LGNTuple.root", "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000003.LGNTuple.root", "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000004.LGNTuple.root", "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000005.LGNTuple.root"
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
            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000008.LGNTuple.root"
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
    
    std::array<double,9> JetNumeratorSFs = {((139e15)*(1.9828e-9)*(0.821204)),((139e15)*(110.64e-12)*(0.69275)),((139e15)*(40.645e-12)*(0.615906)),((139e15)*(1.9817e-9)*(0.1136684)),((139e15)*(110.47e-12)*(0.1912956)),((139e15)*(40.674e-12)*(0.2326772)),((139e15)*(1.9819e-9)*(0.0656969)),((139e15)*(110.53e-12)*(0.1158741)),((139e15)*(40.68e-12)*(0.1535215))}; //numerators for jet bkg
    std::array<double,3> SFs = {((139e15)*(.871e-12)),((139e15)*(.199e-12)), ((139e15)*(.0345e-15))}; //numerators for Z-gamma bkg

    constexpr std::array<const char*,18> Samples = {R"--(pty2\_9\_17)--", R"--(pty\_17\_myy\_0\_80)--", R"--(pty\_17\_myy\_80)--", R"--(Signal $m_{\text{A}}$ = 1 GeV)--", R"--(Signal $m_{\text{A}}$ = 5 GeV)--", R"--(Signal $m_{\text{A}}$ = 2 GeV)--", R"--(Signal $m_{\text{A}}$ = 3 GeV)--", R"--(Signal $m_{\text{A}}$ = 9 GeV)--", "Data", R"--(Zee\_lightJet\_0-70)--", R"--(Zee\_lightJet\_70-140)--", R"--(Zee\_lightJet\_140-280)--", R"--(Zee\_cJet\_0-70)--", R"--(Zee\_cJet\_70-140)--", R"--(Zee\_cJet\_140-280)--", R"--(Zee\_bJet\_0-70)--", R"--(Zee\_bJet\_70-140)--", R"--(Zee\_bJet\_140-280)--"};

    double totalEventsZgamma = 0, resolvedEventsZgamma = 0, SREventsZgamma = 0, SBEventsZgamma = 0;
    double totalEventsZgammaStatUnc = 0, resolvedEventsZgammaStatUnc = 0, SREventsZgammaStatUnc = 0, SBEventsZgammaStatUnc = 0;
    double totalEventsZJets = 0, resolvedEventsZJets = 0, SREventsZJets = 0, SBEventsZJets = 0;
    double totalEventsZJetsStatUnc = 0, resolvedEventsZJetsStatUnc = 0, SREventsZJetsStatUnc = 0, SBEventsZJetsStatUnc = 0;

    std::vector<ROOT::RDF::RResultHandle> Nodes;

    std::ostringstream os;
    os << R"--(\section*{Table 8})--" << '\n';
    os << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
    os << R"--(\setlength\extrarowheight{2pt})--" << '\n';
    os << R"--(\renewcommand{\arraystretch}{1.5})--" << '\n';
    os << R"--(\begin{tabular}{|c|c|c|c|c|})--" << '\n';
    os << R"--(\hline)--" << '\n';
    os << R"--(\textbf{Sample} & \textbf{Total Events} & \textbf{PS: Resolved} & \textbf{SB} & \textbf{SR}
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
        auto leadingPt = opp_charge.Filter([](RVec<Electron>& electrons)
        {
            return ((electrons[0].electron_pt > 20e3 && electrons[1].electron_pt > 27e3) || (electrons[1].electron_pt > 20e3 && electrons[0].electron_pt > 27e3));
        }, {"di_electrons"});

        auto deltaR = leadingPt.Filter([] (RVec<Electron>& electrons)
        {
            return (DeltaR(electrons[0].Vector(), electrons[1].Vector()) > 0.01);
        }, {"di_electrons"});

        //new dataframe node: contains only the events from `dilep` that have di-electron invariant mass between 81 and 101 GeV
        auto mass = deltaR.Filter([] (RVec<Electron>& electrons)
        {
            auto mass = (electrons[0].Vector() + electrons[1].Vector()).M()/1e3;
            return ((mass >= 81) && (mass <= 101));
        }, {"di_electrons"});

        //new dataframe node: contains only the events from `mass` that have dilepton pT > 10 GeV
        auto ptCut = mass.Filter([] (RVec<Electron>& electrons)
        {
            auto pT = (electrons[0].Vector() + electrons[1].Vector()).Pt()/1e3;
            return pT > 10;
        }, {"di_electrons"});

        //end pre-selection -----------

        //photon cuts
        auto resolved = ptCut.Define("photons_pass_cuts",
        [&](RVec<Photon> photons)
        {
            photons.erase(std::remove_if(photons.begin(),photons.end(),
            [](Photon& x)
            {
                return ((std::abs(x.photon_eta) >= 2.37) || (std::abs(x.photon_eta) > 1.37 && std::abs(x.photon_eta) < 1.52) || (!x.photon_id_loose));

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
                //if it's the first combination or if new X is closer to 1
                //than current best_X and ΔR
                //between the two reco-photons < 1.5, then update best_X, pt1, pt2,
                //and the corresponding reco-photons x
                if (i==0 || ((std::abs(1-X) < std::abs(1-best_X)) and (delta_r < 1.5)))
                {
                    best_X = X;
                    pt1 = reco_photons_matched[combs[0][i]].photon_pt;
                    pt2 = reco_photons_matched[combs[1][i]].photon_pt;
                    chosen_delta_r = delta_r;
                }
            }
            if (chosen_delta_r < 1.5 && pt1 > 10e3 && pt2 > 10e3 && best_X > 0.96 && best_X < 1.2) //two photons corresponding to best_X must both have p_T > 10 GeV, ΔR < 1.5, and 0.96 < best_X < 1.2
            {
                return true;
            }
            return false;

        }, {"photons_pass_cuts"})
        .Define("mass", [&](RVec<Electron>& electrons, RVec<Photon>& reco_photons_matched)
        {
            return ((electrons[0].Vector()+electrons[1].Vector()) + (reco_photons_matched[0].Vector()+reco_photons_matched[1].Vector())).M()/1e3;
        }, {"di_electrons", "photons_pass_cuts"});

        auto SB = resolved.Filter(
        [&](double mass)
        {
            return (!((mass > 110) && (mass < 140)));
        }, {"mass"});

        auto SR = resolved.Filter(
        [&](double mass)
        {
            return ((mass > 110) && (mass < 140));
        }, {"mass"});

        Nodes.push_back(df.Count());
        Nodes.push_back(resolved.Count());
        Nodes.push_back(SB.Count());
        Nodes.push_back(SR.Count());
    }

    ROOT::RDF::RunGraphs(Nodes); // running all computation nodes concurrently

    std::cout << Nodes.size() << '\n';

//    0       1       2       3   //Z-gamma
//    4       5       6       7   //Z-gamma
//    8       9       10      11  //Z-gamma
//    12      13      14      15  //ma1
//    16      17      18      19  //ma2
//    20      21      22      23  //ma3
//    24      25      26      27  //ma5
//    28      29      30      31  //ma9
//    32      33      34      35  //data
//    36      37      38      39  //Z+jets
//    40      41      42      43  //Z+jets
//    44      45      46      47  //Z+jets
//    48      49      50      51  //Z+jets
//    52      53      54      55  //Z+jets
//    56      57      58      59  //Z+jets
//    60      61      62      63  //Z+jets
//    64      65      66      67  //Z+jets
//    68      69      70      71  //Z+jets

    os.setf(std::ios::fixed);
    os.precision(2);

    for (int i=0, j=0; (i<18 && j <= 68); i++, j+=4)
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
        else if (i >= 9)
        {
            os  << " & " << *Nodes[j].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+1].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+1].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+2].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+2].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+3].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+3].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
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

    for (int i = 36, j = 0; (i <= 68 && j <= 8); i += 4, j++)
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

    os << R"--(Total $Z\gamma\gamma$ & )--" << totalEventsZgamma
    << R"--($\, \pm \,$)--" << sqrt(totalEventsZgammaStatUnc)
    << " & " << resolvedEventsZgamma
    << R"--($\, \pm \,$)--" << sqrt(resolvedEventsZgammaStatUnc)
    << " & " << SBEventsZgamma
    << R"--($\, \pm \,$)--" << sqrt(SBEventsZgammaStatUnc)
    << " & " << SREventsZgamma
    << R"--($\, \pm \,$)--" << sqrt(SREventsZgammaStatUnc)
    << R"--( \\ \hline )--" << '\n';

    os << R"--(Total $Z$+jets & )--" << totalEventsZJets
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
    std::ofstream out("Table8.txt");
    out << os.str() << '\n';
    out.close();
}

void Table11()
{
    std::vector<std::vector<std::string>> input_filenames = {
        //Z gamma background
        {"/home/common/Za/NTuples/Background/user.kschmied.364860.eegammagamma_pty2_9_17.deriv.DAOD_STDM3.e7057_s3126_r10724_p4092_LGNTuple.root/user.kschmied.31617070._000001.LGNTuple.root"},
        {"/home/common/Za/NTuples/Background/user.kschmied.364861.eegammagamma_pty_17_myy_0_80.deriv.DAOD_HIGG1D2.e7057_s3126_r10724_p4062_LGNTuple.root/user.kschmied.31617064._000001.LGNTuple.root"},
        {"/home/common/Za/NTuples/Background/user.kschmied.364862.eegammagamma_pty_17_myy_80.deriv.DAOD_HIGG1D2.e7057_s3126_r10724_p4204_LGNTuple.root/user.kschmied.31660711._000001.LGNTuple.root", "/home/common/Za/NTuples/Background/user.kschmied.364862.eegammagamma_pty_17_myy_80.deriv.DAOD_HIGG1D2.e7057_s3126_r10724_p4204_LGNTuple.root/user.kschmied.31660711._000002.LGNTuple.root"},
        //Signal
        {"/home/common/Za/NTuples/Signal/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root"}, //1 GeV
        {"/home/common/Za/NTuples/Signal/mc16_13TeV.600751.PhPy8EG_AZNLO_ggH125_mA2p0_v1.root"}, // 2 GeV
        {"/home/common/Za/NTuples/Signal/mc16_13TeV.600752.PhPy8EG_AZNLO_ggH125_mA3p0_v1.root"}, // 3 GeV
        {"/home/common/Za/NTuples/Ntuple_MC_Za_m5p0_v4.root"}, //5 GeV
        {"/home/common/Za/NTuples/Signal/mc16_13TeV.600756.PhPy8EG_AZNLO_ggH125_mA9p0_v1.root"}, // 9 GeV
        //Data
        {"/home/common/Za/NTuples/Ntuple_data_test.root"},
        //Jets
        {
            "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000001.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000002.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000003.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000004.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000005.LGNTuple.root"
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
            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000008.LGNTuple.root"
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
    
    constexpr std::array<const char*,18> Samples = {R"--(pty2\_9\_17)--", R"--(pty\_17\_myy\_0\_80)--", R"--(pty\_17\_myy\_80)--", R"--(Signal $m_{\text{A}}$ = 1 GeV)--", R"--(Signal $m_{\text{A}}$ = 5 GeV)--", R"--(Signal $m_{\text{A}}$ = 2 GeV)--", R"--(Signal $m_{\text{A}}$ = 3 GeV)--", R"--(Signal $m_{\text{A}}$ = 9 GeV)--", "Data", R"--(Zee\_lightJet\_0-70)--", R"--(Zee\_lightJet\_70-140)--", R"--(Zee\_lightJet\_140-280)--", R"--(Zee\_cJet\_0-70)--", R"--(Zee\_cJet\_70-140)--", R"--(Zee\_cJet\_140-280)--", R"--(Zee\_bJet\_0-70)--", R"--(Zee\_bJet\_70-140)--", R"--(Zee\_bJet\_140-280)--"};

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
    os << R"--(\setlength\extrarowheight{2pt})--" << '\n';
    os << R"--(\renewcommand{\arraystretch}{1.5})--" << '\n';
    os << R"--(\begin{tabular}{|c|c|c|c|c|c|c|})--" << '\n';
    os << R"--(\hline)--" << '\n';
    os << R"--(\textbf{Sample} & \textbf{Total Events} & \textbf{pass preselection (PS)} & \textbf{photon} $\pmb{p_T}$ \textbf{+} $\pmb{\Delta R_{\gamma\gamma}}$ \textbf{cut} & $\pmb{X}$ \textbf{window} & \textbf{SR} & \textbf{SR-ID}
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
        auto leading_pt = opp_charge.Filter([](RVec<Electron>& electrons)
        {
            return ((electrons[0].electron_pt > 20e3 && electrons[1].electron_pt > 27e3) || (electrons[1].electron_pt > 20e3 && electrons[0].electron_pt > 27e3));
        }, {"di_electrons"});

//        auto delta_R = leading_pt.Filter([] (RVec<Electron>& electrons)
//        {
//            return (DeltaR(electrons[0].Vector(), electrons[1].Vector()) > 0.01);
//        }, {"di_electrons"});

        auto same_flavour = leading_pt.Filter([] (RVec<Electron>& electrons)
        {
            return true;
        }, {"di_electrons"});

        auto dilep = same_flavour.Define("dilep",[] (RVec<Electron>& electrons)
        {
            return (electrons[0].Vector() + electrons[1].Vector());
        }, {"di_electrons"});

        //new dataframe node: contains only the events from `dilep` that have di-electron invariant mass between 81 and 101 GeV
        auto mass = dilep.Filter([] (RVec<Electron>& electrons)
        {
            auto mass = (electrons[0].Vector() + electrons[1].Vector()).M()/1e3;
            return ((mass >= 81) && (mass <= 101));
        }, {"di_electrons"});

        //new dataframe node: contains only the events from `mass` that have dilepton pT > 10 GeV
        auto ptCut = mass.Filter([] (RVec<Electron>& electrons)
        {
            auto pT = (electrons[0].Vector() + electrons[1].Vector()).Pt()/1e3;
            return pT > 10;
        }, {"di_electrons"});

        //end pre-selection -----------

        //photon cuts
        auto photonPtDeltaR = ptCut.Define("photonPtDeltaR",
        [&](RVec<Photon> photons)
        {
            photons.erase(std::remove_if(photons.begin(),photons.end(),
            [](Photon& x)
            {
                return ((std::abs(x.photon_eta) >= 2.37) || (x.photon_pt <= 10e3) || (std::abs(x.photon_eta) > 1.37 && std::abs(x.photon_eta) < 1.52));

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
                //if it's the first combination or if new X is closer to 1
                //than current best_X and ΔR
                //between the two reco-photons < 1.5, then update best_X, pt1, pt2,
                //and the corresponding reco-photons x
                if (i==0 || ((std::abs(1-X) < std::abs(1-best_X)) and (delta_r < 1.5)))
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
                if (i==0 || ((std::abs(1-X) < std::abs(1-best_X)) and (delta_r < 1.5)))
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

//    0       1       2       3       4       5       //Z-gamma
//    6       7       8       9       10      11      //Z-gamma
//    12      13      14      15      16      17      //Z-gamma
//    18      19      20      21      22      23      //ma1
//    24      25      26      27      28      29      //ma2
//    30      31      32      33      34      35      //ma3
//    36      37      38      39      40      41      //ma5
//    42      43      44      45      46      47      //ma9
//    48      49      50      51      52      53      //data
//    54      55      56      57      58      59      //Z+jets
//    60      61      62      63      64      65      //Z+jets
//    66      67      68      69      70      71      //Z+jets
//    72      73      74      75      76      77      //Z+jets
//    78      79      80      81      82      83      //Z+jets
//    84      85      86      87      88      89      //Z+jets
//    90      91      92      93      94      95      //Z+jets
//    96      97      98      99      100     101     //Z+jets
//    102     103     104     105     106     107     //Z+jets

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

    for (int i=0, j=0; (i<18 && j <= 102); i++, j+=6)
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
        else if (i >= 9)
        {
            os  << " & " << *Nodes[j].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+1].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+1].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+2].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+2].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+3].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+3].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+4].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+4].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
                << " & " << *Nodes[j+5].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
                << R"--($\, \pm \,$)--" << sqrt(*Nodes[j+5].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i-9] / *Nodes[j].GetResultPtr<ULong64_t>())
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

            if (i==6) //5 GeV
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

    for (int i = 54, j = 0; (i <= 102 && j <= 8); i += 6, j++)
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

    os << R"--(Total $Z\gamma\gamma$ & )--" << totalEventsZgamma
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

    os << R"--(Total $Z$+jets & )--" << totalEventsZjets
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
    
    std::ofstream out("Table11.txt");
    out << os.str() << '\n';
    out.close();
}

void CutFlow()
{
    auto start_time = Clock::now();
    
    Table3();
    Table8();
    Table11();
    
    auto end_time = Clock::now();
    std::cout << "Time difference: "
       << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count()/1e9 << " nanoseconds" << std::endl;
    
}


int main()
{
    CutFlow();
}
