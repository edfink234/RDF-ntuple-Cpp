#include <iostream>
#include <vector>
#include <chrono>
#include <algorithm>
#include <unordered_map>
#include <iomanip>
#include <array>
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

//void fig27()
//{
//    auto hs = new THStack("hs","");
//
//    std::vector<std::vector<std::string>> input_filenames = {
//        //Z gamma background
//        {"/home/common/Za/NTuples/Background/user.kschmied.364860.eegammagamma_pty2_9_17.deriv.DAOD_STDM3.e7057_s3126_r10724_p4092_LGNTuple.root/user.kschmied.31617070._000001.LGNTuple.root"},
//        {"/home/common/Za/NTuples/Background/user.kschmied.364861.eegammagamma_pty_17_myy_0_80.deriv.DAOD_HIGG1D2.e7057_s3126_r10724_p4062_LGNTuple.root/user.kschmied.31617064._000001.LGNTuple.root"},
//        {"/home/common/Za/NTuples/Background/user.kschmied.364862.eegammagamma_pty_17_myy_80.deriv.DAOD_HIGG1D2.e7057_s3126_r10724_p4204_LGNTuple.root/user.kschmied.31660711._000001.LGNTuple.root", "/home/common/Za/NTuples/Background/user.kschmied.364862.eegammagamma_pty_17_myy_80.deriv.DAOD_HIGG1D2.e7057_s3126_r10724_p4204_LGNTuple.root/user.kschmied.31660711._000002.LGNTuple.root"},
//        //Data
//        {"/home/common/Za/NTuples/Ntuple_data_test.root"},
//        //Jets
//        {
//            "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000001.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000002.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000003.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000004.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000005.LGNTuple.root"
//        },
//        {
//            "/home/common/Za/NTuples/Background/user.kschmied.364117.v11.Zee_70_140_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835715._000001.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364117.v11.Zee_70_140_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835715._000002.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364117.v11.Zee_70_140_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835715._000003.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364117.v11.Zee_70_140_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835715._000004.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364117.v11.Zee_70_140_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835715._000006.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364117.v11.Zee_70_140_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835715._000007.LGNTuple.root",
//        },
//        {
//            "/home/common/Za/NTuples/Background/user.kschmied.364120.v11.Zee_140_280_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835749._000001.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364120.v11.Zee_140_280_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835749._000002.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364120.v11.Zee_140_280_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835749._000003.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364120.v11.Zee_140_280_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835749._000004.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364120.v11.Zee_140_280_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835749._000005.LGNTuple.root",
//        },
//        {
//            "/home/common/Za/NTuples/Background/user.kschmied.364115.v11.Zee_0_70_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835769._000001.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364115.v11.Zee_0_70_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835769._000002.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364115.v11.Zee_0_70_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835769._000003.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364115.v11.Zee_0_70_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835769._000004.LGNTuple.root",
//        },
//        {
//            "/home/common/Za/NTuples/Background/user.kschmied.364118.v11.Zee_70_140_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835773._000001.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364118.v11.Zee_70_140_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835773._000002.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364118.v11.Zee_70_140_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835773._000003.LGNTuple.root",
//        },
//        {
//            "/home/common/Za/NTuples/Background/user.kschmied.364121.v11.Zee_140_280_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835802._000001.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364121.v11.Zee_140_280_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835802._000002.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364121.v11.Zee_140_280_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835802._000003.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364121.v11.Zee_140_280_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835802._000004.LGNTuple.root",
//        },
//        {
//            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000001.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000002.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000003.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000004.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000005.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000006.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000007.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000008.LGNTuple.root",
//        },
//        {
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000001.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000002.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000003.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000004.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000005.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000006.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000007.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000008.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000009.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000010.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000011.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000012.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000013.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000014.LGNTuple.root",
//        },
//        {
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000001.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000002.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000003.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000004.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000005.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000006.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000007.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000008.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000009.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000010.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000011.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000012.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000013.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000014.LGNTuple.root",
//        },
//    };
//
//    std::array<double,3> SFs = {((139e15)*(.871e-12)),((139e15)*(.199e-12)), ((139e15)*(.0345e-15))}; //numerators for Z-gamma bkg
//
//    std::array<double,9> JetNumeratorSFs = {((139e15)*(1.9828e-9)*(0.821204)),((139e15)*(110.64e-12)*(0.69275)),((139e15)*(40.645e-12)*(0.615906)),((139e15)*(1.9817e-9)*(0.1136684)),((139e15)*(110.47e-12)*(0.1912956)),((139e15)*(40.674e-12)*(0.2326772)),((139e15)*(1.9819e-9)*(0.0656969)),((139e15)*(110.53e-12)*(0.1158741)),((139e15)*(40.68e-12)*(0.1535215))}; //numerators for jet bkg
//
//    std::vector<const char*> prefixes = {"pty2_9_17", "pty_17_myy_0_80", "pty_17_myy_80", "data", "Zee_lightJet_0-70", "Zee_lightJet_70-140", "Zee_lightJet_140-280", "Zee_cJet_0-70", "Zee_cJet_70-140", "Zee_cJet_140-280", "Zee_bJet_0-70", "Zee_bJet_70-140", "Zee_bJet_140-280"};
//    std::vector<EColor> colors = {kBlue, kRed, kViolet, kGreen};
//    std::vector<EColor> Jetscolors = {kCyan, kOrange, kGreen, kYellow, kPink, kGray, kBlack, kSpring, kAzure};
//
//    std::vector<ROOT::RDF::RResultHandle> Nodes;
//    int count = 0;
//    for (auto& i: input_filenames)
//    {
//        SchottDataFrame df(MakeRDF(i, 24));
//
//        auto EventWeight = df.Define("EventWeight",
//        [](const RVec<float>& ei_event_weights_generator)
//        {
//            return  ((ei_event_weights_generator[0]) ? 1 / ei_event_weights_generator[0] : 1);
//
//        }, {"ei_event_weights_generator"});
//
//
//        auto two_leptons = EventWeight.Filter(
//        [](RVec<Muon>& muons, RVec<Electron> electrons)
//        {
//            electrons.erase(std::remove_if(electrons.begin(),electrons.end(),
//            [](Electron& ep)
//            {
//                return (!((ep.electron_pt/1e3 > 20) && (abs(ep.electron_eta) < 2.37) &&
//                          (!((1.37 < abs(ep.electron_eta)) && (abs(ep.electron_eta) < 1.52)))
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
//                return (!((ep.electron_pt/1e3 > 20) && (abs(ep.electron_eta) < 2.37) &&
//                (!((1.37 < abs(ep.electron_eta)) && (abs(ep.electron_eta) < 1.52)))
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
//          photons.erase(std::remove_if(photons.begin(),photons.end(),
//          [](Photon& x)
//          {
//              return ((abs(x.photon_eta) >= 2.37) || (abs(x.photon_eta) > 1.37 && abs(x.photon_eta) < 1.52) || (!x.photon_id_loose));
//
//          }), photons.end());
//
//          return photons;
//        }, {"photons"});
//
//        auto diphotons = photon_passes_cuts.Define("chosen_two",
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
//        }, {"photons_pass_cuts"}).Filter(
//        [&](RVec<Photon>& reco_photons_matched)
//        {
//            return (reco_photons_matched.size()==2);
//        }, {"chosen_two"}).Define("mass",
//        [&](RVec<Photon>& diph)
//        {
//          return (diph[0].Vector()+diph[1].Vector()).M()/1e3;
//        }, {"chosen_two"}).Filter([](double massVal) //Sideband
//        {
//            return (!((massVal > 110) && (massVal < 140)));
//        }, {"mass"}).Define("totEventWeight", [](RVec<float> photon_id_eff, RVec<float> photon_iso_eff, RVec<float> photon_trg_eff/*, RVec<float> ei_event_weights_generator*/)
//        {
//            auto ResizeVal = std::max({photon_id_eff.size(), photon_iso_eff.size(), photon_trg_eff.size()});
//            photon_id_eff.resize(ResizeVal,1);
//            photon_iso_eff.resize(ResizeVal,1);
//            photon_trg_eff.resize(ResizeVal,1);
//
//            return photon_id_eff*photon_iso_eff*photon_trg_eff;//*ei_event_weights_generator[0];
//
//        }, {"photon_id_eff", "photon_iso_eff", "photon_trg_eff",/* "ei_event_weights_generator"*/});
//
//        if (count <= 2 || count >= 4)
//        {
//            Nodes.push_back(diphotons.Sum<RVec<float>>("totEventWeight"));
//            Nodes.push_back(EventWeight.Sum<float>("EventWeight"));
//        }
//
//        Nodes.push_back(diphotons.Histo1D<double>({prefixes[count], prefixes[count++], 60u, 0, 12}, "mass", "totEventWeight"));
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
//    hs->SetMaximum(4000);
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
        SchottDataFrame df(MakeRDF(i, 24));

        auto EventWeight = df.Define("EventWeight",
        [](const RVec<float>& ei_event_weights_generator)
        {
            return  ((ei_event_weights_generator[0]) ? 1 / ei_event_weights_generator[0] : 1);

        }, {"ei_event_weights_generator"});

        auto two_leptons = EventWeight.Filter(
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

        auto diphotons = photon_passes_cuts.Define("chosen_two",
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
        }, {"photons_pass_cuts"}).Filter(
        [&](RVec<Photon>& reco_photons_matched)
        {
            if (reco_photons_matched.size()!=2)
            {
                return false;
            }

            PtEtaPhiEVector rpm = reco_photons_matched[0].Vector()+reco_photons_matched[1].Vector();

            return rpm.M()/1e3 < 2; //looks like their applying this mass cut, but I don't think they say it...

        }, {"chosen_two"}).Define("mass",
        [&](RVec<Photon>& diph)
        {
            return (diph[0].Vector()+diph[1].Vector()).M()/1e3;
        }, {"chosen_two"}).Filter([](double massVal)
        {
           return (!((massVal > 110) && (massVal < 140)));
        }, {"mass"}).Define("deltaR",
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
            return abs((diph[0].Vector() - diph[1].Vector()).Eta());
        }, {"chosen_two"})
        .Define("totEventWeight", [](RVec<float> photon_id_eff, RVec<float> photon_iso_eff, RVec<float> photon_trg_eff/*, RVec<float> ei_event_weights_generator*/)
        {
            auto ResizeVal = std::max({photon_id_eff.size(), photon_iso_eff.size(), photon_trg_eff.size()});
            photon_id_eff.resize(ResizeVal,1);
            photon_iso_eff.resize(ResizeVal,1);
            photon_trg_eff.resize(ResizeVal,1);

            return photon_id_eff*photon_iso_eff*photon_trg_eff;//*ei_event_weights_generator[0];

        }, {"photon_id_eff", "photon_iso_eff", "photon_trg_eff",/* "ei_event_weights_generator"*/});

        if (count <= 2 || count >= 4)
        {
            Nodes.push_back(diphotons.Sum<RVec<float>>("totEventWeight"));
            Nodes.push_back(EventWeight.Sum<float>("EventWeight"));
        }

        Nodes.push_back(diphotons.Histo1D<double>({prefixes[count], prefixes[count], 100u, 0, 12}, "mass", "totEventWeight"));
        Nodes.push_back(diphotons.Histo1D<double>({prefixes[count], prefixes[count], 100u, 0, 0.25}, "deltaR", "totEventWeight"));
        Nodes.push_back(diphotons.Histo1D<double>({prefixes[count], prefixes[count], 100u, 0, 0.2}, "deltaPhi", "totEventWeight"));
        Nodes.push_back(diphotons.Histo1D<double>({prefixes[count], prefixes[count++], 100u, 0, 0.04}, "deltaEta", "totEventWeight"));
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
    TLegend* legend = new TLegend(0.475, 0.2, 0.85, 0.6);
    count = 0;
    //Z-gamma
    for (auto& i: {2,8,14,18})
    {
        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0 && i != 18)
        {
            Nodes[i].GetResultPtr<TH1D>()->Scale(SFs[count]/ *Nodes[i-1].GetResultPtr<float>());
        }
        Nodes[i].GetResultPtr<TH1D>()->SetFillColor(colors[count++]);
        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "f");

        if (i != 18)
        {
            hs->Add(&*Nodes[i].GetResultPtr<TH1D>());
        }
    }
    //Z-jets
    count = 0;
    for (int i = 24; i <= 72; i += 6)
    {
        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0)
        {
            Nodes[i].GetResultPtr<TH1D>()->Scale(JetNumeratorSFs[count]/ *Nodes[i-1].GetResultPtr<float>());
        }
        Nodes[i].GetResultPtr<TH1D>()->SetFillColor(Jetscolors[count++]);
        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "f");

        hs->Add(&*Nodes[i].GetResultPtr<TH1D>());
    }

    hs->Draw("HIST");
    Nodes[18].GetResultPtr<TH1D>()->Draw("HISTsame");
    hs->SetMinimum(0);
    hs->SetMaximum(2300);
    hs->SetTitle(";m_{#gamma#gamma} [GeV];Events");
    hs->GetYaxis()->CenterTitle(true);
    hs->GetXaxis()->SetTitleOffset(1.2);
    gStyle->SetOptStat(0);
    TLatex Tl;
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.6, 0.8, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.6, 0.7,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("Fig28A.pdf");

    hs = new THStack("hs2","");
    c1 = new TCanvas();
    legend = new TLegend(0.5, 0.2, 0.875, 0.6);
    count = 0;
    //Z-gamma
    for (auto& i: {3,9,15,19})
    {
        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0 && i != 19)
        {
            Nodes[i].GetResultPtr<TH1D>()->Scale(SFs[count]/ *Nodes[i-2].GetResultPtr<float>());
        }
        Nodes[i].GetResultPtr<TH1D>()->SetFillColor(colors[count++]);
        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "f");

        if (i != 19)
        {
            hs->Add(&*Nodes[i].GetResultPtr<TH1D>());
        }
    }

    //Z-jets
    count = 0;
    for (int i = 25; i <= 73; i += 6)
    {
        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0)
        {
            Nodes[i].GetResultPtr<TH1D>()->Scale(JetNumeratorSFs[count]/ *Nodes[i-2].GetResultPtr<float>());
        }
        Nodes[i].GetResultPtr<TH1D>()->SetFillColor(Jetscolors[count++]);
        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "f");

        hs->Add(&*Nodes[i].GetResultPtr<TH1D>());
    }

    hs->Draw("HIST");
    Nodes[19].GetResultPtr<TH1D>()->Draw("HISTsame");
    hs->SetMinimum(0);
    hs->SetMaximum(1000);
    hs->SetTitle(";#DeltaR_{#gamma#gamma};Events");
    hs->GetYaxis()->CenterTitle(true);
    hs->GetXaxis()->SetTitleOffset(1.2);
    gStyle->SetOptStat(0);
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.6, 0.8, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.6, 0.7,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("Fig28B.pdf");

    hs = new THStack("hs3","");
    c1 = new TCanvas();
    legend = new TLegend(0.575, 0.2125, 0.9, 0.6125);
    count = 0;
    for (auto& i: {4,10,16,20})
    {
        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0 && i != 20)
        {
            Nodes[i].GetResultPtr<TH1D>()->Scale(SFs[count]/ *Nodes[i-3].GetResultPtr<float>());
        }
        Nodes[i].GetResultPtr<TH1D>()->SetFillColor(colors[count++]);
        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "f");

        if (i != 20)
        {
            hs->Add(&*Nodes[i].GetResultPtr<TH1D>());
        }
    }

    //Z-jets
    count = 0;
    for (int i = 26; i <= 74; i += 6)
    {
        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0)
        {
            Nodes[i].GetResultPtr<TH1D>()->Scale(JetNumeratorSFs[count]/ *Nodes[i-3].GetResultPtr<float>());
        }
        Nodes[i].GetResultPtr<TH1D>()->SetFillColor(Jetscolors[count++]);
        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "f");

        hs->Add(&*Nodes[i].GetResultPtr<TH1D>());
    }

    hs->Draw("HIST");
    Nodes[20].GetResultPtr<TH1D>()->Draw("HISTsame");
    hs->SetMinimum(0);
    hs->SetMaximum(500);
    hs->SetTitle(";#Delta#phi_{#gamma#gamma};Events");
    hs->GetYaxis()->CenterTitle(true);
    hs->GetXaxis()->SetTitleOffset(1.2);
    gStyle->SetOptStat(0);
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.6, 0.8, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.6, 0.7,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("Fig28C.pdf");

    hs = new THStack("hs4","");
    c1 = new TCanvas();
    legend = new TLegend(0.125, 0.4, 0.5, 0.8);
    count = 0;
    for (auto& i: {5,11,17,21})
    {
        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0 && i != 21)
        {
            Nodes[i].GetResultPtr<TH1D>()->Scale(SFs[count]/ *Nodes[i-4].GetResultPtr<float>());
        }
        Nodes[i].GetResultPtr<TH1D>()->SetFillColor(colors[count++]);
        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "f");

        if (i != 21)
        {
            hs->Add(&*Nodes[i].GetResultPtr<TH1D>());
        }
    }

    //Z-jets
    count = 0;
    for (int i = 27; i <= 75; i += 6)
    {
        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0)
        {
            Nodes[i].GetResultPtr<TH1D>()->Scale(JetNumeratorSFs[count]/ *Nodes[i-4].GetResultPtr<float>());
        }
        Nodes[i].GetResultPtr<TH1D>()->SetFillColor(Jetscolors[count++]);
        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "f");

        hs->Add(&*Nodes[i].GetResultPtr<TH1D>());
    }

    hs->Draw("HIST");
    Nodes[21].GetResultPtr<TH1D>()->Draw("HISTsame");
    hs->SetMinimum(0);
    hs->SetMaximum(218);
    hs->SetTitle(";#Delta#eta_{#gamma#gamma};Events");
    hs->GetYaxis()->CenterTitle(true);
    hs->GetXaxis()->SetTitleOffset(1.2);
    gStyle->SetOptStat(0);
    Tl.SetTextSize(0.03);
    Tl.DrawLatexNDC(0.6, 0.8, "#it{ATLAS} Internal");
    Tl.DrawLatexNDC(0.6, 0.7,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("Fig28D.pdf");
}

//void fig41()
//{
//    std::vector<std::vector<std::string>> input_filenames =
//    {
//        //Signal
//        {"/home/common/Za/NTuples/Ntuple_MC_Za_m5p0_v4.root"},
//        {"/home/common/Haa/ntuples/Za/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root"},
//        //Z gamma background
//        {"/home/common/Za/NTuples/Background/user.kschmied.364860.eegammagamma_pty2_9_17.deriv.DAOD_STDM3.e7057_s3126_r10724_p4092_LGNTuple.root/user.kschmied.31617070._000001.LGNTuple.root"},
//        {"/home/common/Za/NTuples/Background/user.kschmied.364861.eegammagamma_pty_17_myy_0_80.deriv.DAOD_HIGG1D2.e7057_s3126_r10724_p4062_LGNTuple.root/user.kschmied.31617064._000001.LGNTuple.root"},
//        {"/home/common/Za/NTuples/Background/user.kschmied.364862.eegammagamma_pty_17_myy_80.deriv.DAOD_HIGG1D2.e7057_s3126_r10724_p4204_LGNTuple.root/user.kschmied.31660711._000001.LGNTuple.root", "/home/common/Za/NTuples/Background/user.kschmied.364862.eegammagamma_pty_17_myy_80.deriv.DAOD_HIGG1D2.e7057_s3126_r10724_p4204_LGNTuple.root/user.kschmied.31660711._000002.LGNTuple.root"},
//        //Data
//        {"/home/common/Za/NTuples/Ntuple_data_test.root"},
//        //Jets
//        {
//            "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000001.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000002.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000003.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000004.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000005.LGNTuple.root",
//        },
//        {
//            "/home/common/Za/NTuples/Background/user.kschmied.364117.v11.Zee_70_140_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835715._000001.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364117.v11.Zee_70_140_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835715._000002.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364117.v11.Zee_70_140_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835715._000003.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364117.v11.Zee_70_140_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835715._000004.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364117.v11.Zee_70_140_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835715._000006.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364117.v11.Zee_70_140_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835715._000007.LGNTuple.root",
//        },
//        {
//            "/home/common/Za/NTuples/Background/user.kschmied.364120.v11.Zee_140_280_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835749._000001.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364120.v11.Zee_140_280_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835749._000002.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364120.v11.Zee_140_280_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835749._000003.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364120.v11.Zee_140_280_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835749._000004.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364120.v11.Zee_140_280_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835749._000005.LGNTuple.root",
//        },
//        {
//            "/home/common/Za/NTuples/Background/user.kschmied.364115.v11.Zee_0_70_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835769._000001.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364115.v11.Zee_0_70_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835769._000002.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364115.v11.Zee_0_70_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835769._000003.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364115.v11.Zee_0_70_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835769._000004.LGNTuple.root",
//        },
//        {
//            "/home/common/Za/NTuples/Background/user.kschmied.364118.v11.Zee_70_140_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835773._000001.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364118.v11.Zee_70_140_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835773._000002.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364118.v11.Zee_70_140_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835773._000003.LGNTuple.root",
//        },
//        {
//            "/home/common/Za/NTuples/Background/user.kschmied.364121.v11.Zee_140_280_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835802._000001.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364121.v11.Zee_140_280_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835802._000002.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364121.v11.Zee_140_280_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835802._000003.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364121.v11.Zee_140_280_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835802._000004.LGNTuple.root",
//        },
//        {
//            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000001.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000002.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000003.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000004.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000005.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000006.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000007.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000008.LGNTuple.root",
//        },
//        {
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000001.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000002.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000003.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000004.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000005.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000006.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000007.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000008.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000009.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000010.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000011.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000012.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000013.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000014.LGNTuple.root",
//        },
//        {
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000001.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000002.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000003.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000004.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000005.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000006.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000007.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000008.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000009.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000010.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000011.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000012.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000013.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000014.LGNTuple.root",
//        },
//    };
//
//    std::array<double,3> SFs = {((139e15)*(.871e-12)),((139e15)*(.199e-12)), ((139e15)*(.0345e-15))}; //numerators for Z-gamma bkg
//
//    std::array<double,9> JetNumeratorSFs = {((139e15)*(1.9828e-9)*(0.821204)),((139e15)*(110.64e-12)*(0.69275)),((139e15)*(40.645e-12)*(0.615906)),((139e15)*(1.9817e-9)*(0.1136684)),((139e15)*(110.47e-12)*(0.1912956)),((139e15)*(40.674e-12)*(0.2326772)),((139e15)*(1.9819e-9)*(0.0656969)),((139e15)*(110.53e-12)*(0.1158741)),((139e15)*(40.68e-12)*(0.1535215))}; //numerators for jet bkg
//
//    std::vector<const char*> prefixes = {"Sig m_{A} = 5 GeV", "Sig m_{A} = 1 GeV", "pty2_9_17", "pty_17_myy_0_80", "pty_17_myy_80", "data", "Zee_lightJet_0-70", "Zee_lightJet_70-140", "Zee_lightJet_140-280", "Zee_cJet_0-70", "Zee_cJet_70-140", "Zee_cJet_140-280", "Zee_bJet_0-70", "Zee_bJet_70-140", "Zee_bJet_140-280"};
//    std::vector<EColor> colors = {kBlack, kMagenta, kBlue, kRed, kViolet, kGreen};
//    std::vector<EColor> Jetscolors = {kCyan, kOrange, kGreen, kYellow, kPink, kGray, kBlack, kSpring, kAzure};
//
//    std::vector<ROOT::RDF::RResultHandle> Nodes;
//
//    int count = 0;
//    for (auto& i: input_filenames)
//    {
//        SchottDataFrame df(MakeRDF(i, 24));
//
//        auto EventWeight = df.Define("EventWeight",
//        [](const RVec<float>& ei_event_weights_generator)
//        {
//            return  ((ei_event_weights_generator[0]) ? 1 / ei_event_weights_generator[0] : 1);
//
//        }, {"ei_event_weights_generator"});
//
//        auto two_leptons = EventWeight.Filter(
//        [](RVec<Muon>& muons, RVec<Electron> electrons)
//        {
//            electrons.erase(std::remove_if(electrons.begin(),electrons.end(),
//            [](Electron& ep)
//            {
//                return (!((ep.electron_pt/1e3 > 20) && (abs(ep.electron_eta) < 2.37) &&
//                          (!((1.37 < abs(ep.electron_eta)) && (abs(ep.electron_eta) < 1.52)))
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
//                return (!((ep.electron_pt/1e3 > 20) && (abs(ep.electron_eta) < 2.37) &&
//                (!((1.37 < abs(ep.electron_eta)) && (abs(ep.electron_eta) < 1.52)))
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
//              return ((abs(x.photon_eta) >= 2.37) || (abs(x.photon_eta) > 1.37 && abs(x.photon_eta) < 1.52) || (!x.photon_id_loose));
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
//                if (i==0 || abs(1-X) < abs(1-best_X))
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
//        }, {"di_electrons", "merged_photon"})
//        .Define("totEventWeight", [](RVec<float> photon_id_eff, RVec<float> photon_iso_eff, RVec<float> photon_trg_eff/*, RVec<float> ei_event_weights_generator*/)
//        {
//            auto ResizeVal = std::max({photon_id_eff.size(), photon_iso_eff.size(), photon_trg_eff.size()});
//            photon_id_eff.resize(ResizeVal,1);
//            photon_iso_eff.resize(ResizeVal,1);
//            photon_trg_eff.resize(ResizeVal,1);
//
//            return photon_id_eff*photon_iso_eff*photon_trg_eff;//*ei_event_weights_generator[0];
//
//        }, {"photon_id_eff", "photon_iso_eff", "photon_trg_eff",/* "ei_event_weights_generator"*/});
//
//        if (count < 2) //signal only
//        {
//            Nodes.push_back(dilepton_and_photon.Histo1D<double>({prefixes[count], prefixes[count++], 100u, 80, 200}, "reconstructed_mass", "totEventWeight"));
//        }
//        else
//        {
//            if ((count >= 2 && count <= 4) || (count >= 6))
//            {
//                Nodes.push_back(dilepton_and_photon.Sum<RVec<float>>("totEventWeight"));
//                Nodes.push_back(EventWeight.Sum<float>("EventWeight"));
//            }
//
//            Nodes.push_back(dilepton_and_photon.Histo1D<double>({prefixes[count], prefixes[count++], 100u, 80, 410}, "reconstructed_mass", "totEventWeight"));
//        }
//    }
//
////    0               ma5
////    1               ma1
////    2   3   4       Z-gamma
////    5   6   7       Z-gamma
////    8   9   10      Z-gamma
////    11              data
////    12  13  14      Z-jets
////    15  16  17      Z-jets
////    18  19  20      Z-jets
////    21  22  23      Z-jets
////    24  25  26      Z-jets
////    27  28  29      Z-jets
////    30  31  32      Z-jets
////    33  34  35      Z-jets
////    36  37  38      Z-jets
//
//    ROOT::RDF::RunGraphs(Nodes); // running all computation nodes concurrently
//
//    count = 0;
//    TCanvas* c1 = new TCanvas();
//    TLegend* legend = new TLegend(0.65, 0.4, 0.85, 0.6);
//    double factor;
//    //signal
//    for (auto& i: {0,1})
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
//            Nodes[i].GetResultPtr<TH1D>()->SetAxisRange(0., 250, "Y");
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
//    Tl.DrawLatexNDC(0.6, 0.7,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
//    legend->SetBorderSize(0);
//    legend->Draw();
//    c1->SaveAs("Fig41B.pdf");
//
////    count = 0;
//    factor = 0;
//    int back_count = 0;
//    for (auto& i: {2,5,8})
//    {
//        factor += (*Nodes[i].GetResultPtr<float>())*(SFs[back_count++] / *Nodes[i+1].GetResultPtr<float>());
//    }
//
//    auto hs = new THStack("hs3","");
//    c1 = new TCanvas();
//    legend = new TLegend(0.5, 0.2, 0.85, 0.6);
//
//    //Z-gamma
//    for (auto& i: {4,7,10,11})
//    {
//        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0 && i != 11)
//        {
//            Nodes[i].GetResultPtr<TH1D>()->Scale(SFs[count-2] / *Nodes[i-1].GetResultPtr<float>());
//        }
//        Nodes[i].GetResultPtr<TH1D>()->SetFillColor(colors[count++]);
//        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "f");
//
//        if (i != 11)
//        {
//            hs->Add(&*Nodes[i].GetResultPtr<TH1D>());
//        }
//    }
//
//    for (int i = 14, j = 0; (i <= 38 && j <= 8); i += 3, j++)
//    {
//        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0)
//        {
//            Nodes[i].GetResultPtr<TH1D>()->Scale(JetNumeratorSFs[j] / *Nodes[i-1].GetResultPtr<float>());
//        }
//        Nodes[i].GetResultPtr<TH1D>()->SetFillColor(Jetscolors[j]);
//        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "f");
//        hs->Add(&*Nodes[i].GetResultPtr<TH1D>());
//    }
//
//    hs->Draw("HIST");
//    Nodes[11].GetResultPtr<TH1D>()->Draw("HISTsame");
//    count=0;
//
//    for (auto &i: {0,1})
//    {
//        Nodes[i].GetResultPtr<TH1D>()->SetLineWidth(2);
//        Nodes[i].GetResultPtr<TH1D>()->SetLineColor(colors[count]);
//        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), prefixes[count++], "l");
//        Nodes[i].GetResultPtr<TH1D>()->SetTitle(";m_{ll#gamma} [GeV];Events");
//        Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);
//        Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->SetTitleOffset(1.2);
//        Nodes[i].GetResultPtr<TH1D>()->DrawClone("HISTsame");
//        gPad->Modified();
//        gPad->Update();
//        c1->Modified();
//        c1->Update();
//    }
//    hs->SetMaximum(6.8e5);
//    hs->SetTitle(";m_{ll#gamma} [GeV];Events");
//    hs->GetYaxis()->CenterTitle(true);
//    hs->GetXaxis()->SetTitleOffset(1.2);
//
//    gStyle->SetOptStat(0);
//    Tl.SetTextSize(0.03);
//    Tl.DrawLatexNDC(0.6, 0.8, "#it{ATLAS} Internal");
//    Tl.DrawLatexNDC(0.6, 0.7,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
//    legend->SetBorderSize(0);
//    legend->Draw("same");
//    c1->SaveAs("Fig41A.pdf");
//}
//
//void fig48()
//{
//    std::vector<std::vector<std::string>> input_filenames = {
//        //Z gamma background
//        {"/home/common/Za/NTuples/Background/user.kschmied.364860.eegammagamma_pty2_9_17.deriv.DAOD_STDM3.e7057_s3126_r10724_p4092_LGNTuple.root/user.kschmied.31617070._000001.LGNTuple.root"},
//        {"/home/common/Za/NTuples/Background/user.kschmied.364861.eegammagamma_pty_17_myy_0_80.deriv.DAOD_HIGG1D2.e7057_s3126_r10724_p4062_LGNTuple.root/user.kschmied.31617064._000001.LGNTuple.root"},
//        {"/home/common/Za/NTuples/Background/user.kschmied.364862.eegammagamma_pty_17_myy_80.deriv.DAOD_HIGG1D2.e7057_s3126_r10724_p4204_LGNTuple.root/user.kschmied.31660711._000001.LGNTuple.root", "/home/common/Za/NTuples/Background/user.kschmied.364862.eegammagamma_pty_17_myy_80.deriv.DAOD_HIGG1D2.e7057_s3126_r10724_p4204_LGNTuple.root/user.kschmied.31660711._000002.LGNTuple.root"},
//        //Data
//        {"/home/common/Za/NTuples/Ntuple_data_test.root"},
//        //Jets
//        {
//            "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000001.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000002.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000003.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000004.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000005.LGNTuple.root",
//        },
//        {
//            "/home/common/Za/NTuples/Background/user.kschmied.364117.v11.Zee_70_140_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835715._000001.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364117.v11.Zee_70_140_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835715._000002.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364117.v11.Zee_70_140_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835715._000003.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364117.v11.Zee_70_140_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835715._000004.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364117.v11.Zee_70_140_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835715._000006.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364117.v11.Zee_70_140_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835715._000007.LGNTuple.root",
//        },
//        {
//            "/home/common/Za/NTuples/Background/user.kschmied.364120.v11.Zee_140_280_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835749._000001.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364120.v11.Zee_140_280_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835749._000002.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364120.v11.Zee_140_280_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835749._000003.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364120.v11.Zee_140_280_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835749._000004.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364120.v11.Zee_140_280_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835749._000005.LGNTuple.root",
//        },
//        {
//            "/home/common/Za/NTuples/Background/user.kschmied.364115.v11.Zee_0_70_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835769._000001.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364115.v11.Zee_0_70_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835769._000002.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364115.v11.Zee_0_70_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835769._000003.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364115.v11.Zee_0_70_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835769._000004.LGNTuple.root",
//        },
//        {
//            "/home/common/Za/NTuples/Background/user.kschmied.364118.v11.Zee_70_140_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835773._000001.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364118.v11.Zee_70_140_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835773._000002.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364118.v11.Zee_70_140_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835773._000003.LGNTuple.root",
//        },
//        {
//            "/home/common/Za/NTuples/Background/user.kschmied.364121.v11.Zee_140_280_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835802._000001.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364121.v11.Zee_140_280_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835802._000002.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364121.v11.Zee_140_280_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835802._000003.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364121.v11.Zee_140_280_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835802._000004.LGNTuple.root",
//        },
//        {
//            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000001.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000002.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000003.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000004.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000005.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000006.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000007.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000008.LGNTuple.root",
//        },
//        {
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000001.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000002.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000003.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000004.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000005.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000006.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000007.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000008.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000009.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000010.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000011.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000012.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000013.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000014.LGNTuple.root",
//        },
//        {
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000001.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000002.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000003.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000004.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000005.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000006.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000007.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000008.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000009.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000010.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000011.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000012.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000013.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000014.LGNTuple.root",
//        },
//    };
//    std::array<double,3> SFs = {((139e15)*(.871e-12)),((139e15)*(.199e-12)), ((139e15)*(.0345e-15))}; //numerators for Z-gamma bkg
//
//    std::array<double,9> JetNumeratorSFs = {((139e15)*(1.9828e-9)*(0.821204)),((139e15)*(110.64e-12)*(0.69275)),((139e15)*(40.645e-12)*(0.615906)),((139e15)*(1.9817e-9)*(0.1136684)),((139e15)*(110.47e-12)*(0.1912956)),((139e15)*(40.674e-12)*(0.2326772)),((139e15)*(1.9819e-9)*(0.0656969)),((139e15)*(110.53e-12)*(0.1158741)),((139e15)*(40.68e-12)*(0.1535215))}; //numerators for jet bkg
//
//    std::vector<const char*> prefixes = {"pty2_9_17", "pty_17_myy_0_80", "pty_17_myy_80", "data", "Zee_lightJet_0-70", "Zee_lightJet_70-140", "Zee_lightJet_140-280", "Zee_cJet_0-70", "Zee_cJet_70-140", "Zee_cJet_140-280", "Zee_bJet_0-70", "Zee_bJet_70-140", "Zee_bJet_140-280"};
//    std::vector<EColor> colors = {kBlue, kRed, kViolet, kGreen};
//    std::vector<EColor> Jetscolors = {kCyan, kOrange, kGreen, kYellow, kPink, kGray, kBlack, kSpring, kAzure};
//
//    std::vector<ROOT::RDF::RResultHandle> Nodes;
//
//    int count = 0;
//    for (auto& i: input_filenames)
//    {
//        SchottDataFrame df(MakeRDF(i, 24));
//
//        auto EventWeight = df.Define("EventWeight",
//        [](const RVec<float>& ei_event_weights_generator)
//        {
//            return  ((ei_event_weights_generator[0]) ? 1 / ei_event_weights_generator[0] : 1);
//
//        }, {"ei_event_weights_generator"});
//
//        auto two_leptons = EventWeight.Filter(
//        [](RVec<Muon>& muons, RVec<Electron> electrons)
//        {
//            electrons.erase(std::remove_if(electrons.begin(),electrons.end(),
//            [](Electron& ep)
//            {
//                return (!((ep.electron_pt/1e3 > 20) && (abs(ep.electron_eta) < 2.37) &&
//                          (!((1.37 < abs(ep.electron_eta)) && (abs(ep.electron_eta) < 1.52)))
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
//                return (!((ep.electron_pt/1e3 > 20) && (abs(ep.electron_eta) < 2.37) &&
//                (!((1.37 < abs(ep.electron_eta)) && (abs(ep.electron_eta) < 1.52)))
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
//              return ((abs(x.photon_eta) >= 2.37) || (abs(x.photon_eta) > 1.37 && abs(x.photon_eta) < 1.52) || (!x.photon_id_loose));
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
//                if (i==0 || abs(1-X) < abs(1-best_X))
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
//        .Define("dilepton_mass",[&](RVec<Electron>& di_electrons)
//        {
//            auto four_momentum = di_electrons[0].Vector() + di_electrons[1].Vector();
//
//            return four_momentum.M()/1e3;
//
//        }, {"di_electrons"})
//        .Define("merged_photon_pt",[&](Photon& merged_photon)
//        {
//            return merged_photon.Vector().Pt()/1e3;
//
//        }, {"merged_photon"})
//        .Define("reconstructed_mass",[&](RVec<Electron>& di_electrons, Photon& merged_photon)
//        {
//         auto four_momentum = di_electrons[0].Vector() + di_electrons[1].Vector();
//
//         return (four_momentum + merged_photon.Vector()).M()/1e3;
//
//        }, {"di_electrons", "merged_photon"})
//        .Filter([](RVec<float>& Eratio, double reconstructed_mass)
//        {
//            return ((!Any(Eratio <= 0.8)) && ((reconstructed_mass <= 110) || (reconstructed_mass >= 130)));
//        }, {"photon_shower_shape_e_ratio", "reconstructed_mass"})
//        .Define("totEventWeight", [](RVec<float> photon_id_eff, RVec<float> photon_iso_eff, RVec<float> photon_trg_eff/*, RVec<float> ei_event_weights_generator*/)
//        {
//            auto ResizeVal = std::max({photon_id_eff.size(), photon_iso_eff.size(), photon_trg_eff.size()});
//            photon_id_eff.resize(ResizeVal,1);
//            photon_iso_eff.resize(ResizeVal,1);
//            photon_trg_eff.resize(ResizeVal,1);
//
//            return photon_id_eff*photon_iso_eff*photon_trg_eff;//*ei_event_weights_generator[0];
//
//        }, {"photon_id_eff", "photon_iso_eff", "photon_trg_eff",/* "ei_event_weights_generator"*/});
//
//        if (count <= 2 || count >= 4)
//        {
//            Nodes.push_back(dilepton_and_photon.Sum<RVec<float>>("totEventWeight"));
//            Nodes.push_back(EventWeight.Sum<float>("EventWeight"));
//        }
//
//        Nodes.push_back(dilepton_and_photon.Histo1D<double>({prefixes[count], prefixes[count], 60u, 60, 120}, "dilepton_mass", "totEventWeight"));
//        Nodes.push_back(dilepton_and_photon.Histo1D<double>({prefixes[count], prefixes[count++], 60u, 0, 165}, "merged_photon_pt", "totEventWeight"));
//    }
//
////    0   1   2   3    //Z-gamma
////    4   5   6   7    //Z-gamma
////    8   9   10  11   //Z-gamma
////            12  13   //data
////    14  15  16  17   //Z-jets
////    18  19  20  21   //Z-jets
////    22  23  24  25   //Z-jets
////    26  27  28  29   //Z-jets
////    30  31  32  33   //Z-jets
////    34  35  36  37   //Z-jets
////    38  39  40  41   //Z-jets
////    42  43  44  45   //Z-jets
////    46  47  48  49   //Z-jets
//
//    ROOT::RDF::RunGraphs(Nodes); // running all computation nodes concurrently
//
//    count = 0;
//    TCanvas* c1 = new TCanvas();
//    TLegend* legend = new TLegend(0.15, 0.275, 0.475, 0.675);
//    auto hs = new THStack("hs1","");
//    double factor = 0;
//    int back_count = 0;
//    //Z-gamma
//    for (auto& i: {0,4,8})
//    {
//        factor += (*Nodes[i].GetResultPtr<float>())*(SFs[back_count++]/ *Nodes[i+1].GetResultPtr<float>());
//    }
//    //Z-jets
//    for (int i = 14, j = 0; (i <= 46 && j <= 8); i+=4, j++)
//    {
//        factor += (*Nodes[i].GetResultPtr<float>())*(JetNumeratorSFs[j]/ *Nodes[i+1].GetResultPtr<float>());
//    }
//    //Z-gamma
//    for (auto& i: {2,6,10,12})
//    {
//        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0 && i != 12)
//        {
//            Nodes[i].GetResultPtr<TH1D>()->Scale(SFs[count]/ *Nodes[i-1].GetResultPtr<float>());
//        }
//        Nodes[i].GetResultPtr<TH1D>()->SetFillColor(colors[count++]);
//        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "f");
//
//        if (i != 12)
//        {
//            hs->Add(&*Nodes[i].GetResultPtr<TH1D>());
//        }
//    }
//
//    //Z-jets
//    for (int i = 16, j = 0; (i <= 48 && j <= 8); i+=4, j++)
//    {
//        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0)
//        {
//            Nodes[i].GetResultPtr<TH1D>()->Scale(JetNumeratorSFs[j]/ *Nodes[i-1].GetResultPtr<float>());
//        }
//        Nodes[i].GetResultPtr<TH1D>()->SetFillColor(Jetscolors[j]);
//        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "f");
//
//        hs->Add(&*Nodes[i].GetResultPtr<TH1D>());
//    }
//
//    hs->Draw("HIST");
//    Nodes[12].GetResultPtr<TH1D>()->Draw("HISTsame");
//    hs->SetTitle(";m_{ll} [GeV];Events");
//    hs->GetYaxis()->CenterTitle(true);
//    hs->GetXaxis()->SetTitleOffset(1.2);
//    hs->GetYaxis()->SetTitleOffset(1.4);
//    gStyle->SetOptStat(0);
//    TLatex Tl;
//    Tl.SetTextSize(0.03);
//    Tl.DrawLatexNDC(0.6, 0.8, "#it{ATLAS} Internal");
//    Tl.DrawLatexNDC(0.6, 0.7,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
//    legend->SetBorderSize(0);
//    legend->Draw();
//    c1->SaveAs("Fig48A.pdf");
//    count = 0;
//    hs = new THStack("hs2","");
//    c1 = new TCanvas();
//    legend = new TLegend(0.55, 0.2, 0.85, 0.6);
//    //Z-gamma
//    for (auto& i: {3,7,11,13})
//    {
//        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0 && i != 13)
//        {
//            Nodes[i].GetResultPtr<TH1D>()->Scale(SFs[count]/ *Nodes[i-2].GetResultPtr<float>());
//        }
//        Nodes[i].GetResultPtr<TH1D>()->SetFillColor(colors[count++]);
//        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "f");
//
//        if (i != 13)
//        {
//            hs->Add(&*Nodes[i].GetResultPtr<TH1D>());
//        }
//    }
//    //Z-Jets
//    for (int i = 17, j = 0; (i <= 49 && j <= 8); i+=4, j++)
//    {
//        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0)
//        {
//            Nodes[i].GetResultPtr<TH1D>()->Scale(JetNumeratorSFs[j]/ *Nodes[i-2].GetResultPtr<float>());
//        }
//        Nodes[i].GetResultPtr<TH1D>()->SetFillColor(Jetscolors[j]);
//        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "f");
//
//        hs->Add(&*Nodes[i].GetResultPtr<TH1D>());
//    }
//
//    hs->Draw("HIST");
//    Nodes[13].GetResultPtr<TH1D>()->Draw("HISTsame");
//    hs->SetTitle(";photon p_{T} [GeV];Events");
//    hs->GetYaxis()->CenterTitle(true);
//    hs->GetXaxis()->SetTitleOffset(1.2);
//    hs->GetYaxis()->SetTitleOffset(1.4);
//    gStyle->SetOptStat(0);
//    Tl.SetTextSize(0.03);
//    Tl.DrawLatexNDC(0.6, 0.8, "#it{ATLAS} Internal");
//    Tl.DrawLatexNDC(0.6, 0.7,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
//    legend->SetBorderSize(0);
//    legend->Draw();
//    c1->SaveAs("Fig48B.pdf");
//}
//
//void fig59()
//{
//    std::vector<std::vector<std::string>> input_filenames = {
//        //Signal
//        {"/home/common/Za/NTuples/Ntuple_MC_Za_m5p0_v4.root"},
//        {"/home/common/Haa/ntuples/Za/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root"},
//        //Z gamma background
//        {"/home/common/Za/NTuples/Background/user.kschmied.364860.eegammagamma_pty2_9_17.deriv.DAOD_STDM3.e7057_s3126_r10724_p4092_LGNTuple.root/user.kschmied.31617070._000001.LGNTuple.root"},
//        {"/home/common/Za/NTuples/Background/user.kschmied.364861.eegammagamma_pty_17_myy_0_80.deriv.DAOD_HIGG1D2.e7057_s3126_r10724_p4062_LGNTuple.root/user.kschmied.31617064._000001.LGNTuple.root"},
//        {"/home/common/Za/NTuples/Background/user.kschmied.364862.eegammagamma_pty_17_myy_80.deriv.DAOD_HIGG1D2.e7057_s3126_r10724_p4204_LGNTuple.root/user.kschmied.31660711._000001.LGNTuple.root", "/home/common/Za/NTuples/Background/user.kschmied.364862.eegammagamma_pty_17_myy_80.deriv.DAOD_HIGG1D2.e7057_s3126_r10724_p4204_LGNTuple.root/user.kschmied.31660711._000002.LGNTuple.root"},
//        //Data
//        {"/home/common/Za/NTuples/Ntuple_data_test.root"},
//        //Jets
//        {
//            "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000001.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000002.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000003.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000004.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000005.LGNTuple.root"
//        },
//        {
//            "/home/common/Za/NTuples/Background/user.kschmied.364117.v11.Zee_70_140_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835715._000001.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364117.v11.Zee_70_140_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835715._000002.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364117.v11.Zee_70_140_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835715._000003.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364117.v11.Zee_70_140_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835715._000004.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364117.v11.Zee_70_140_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835715._000006.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364117.v11.Zee_70_140_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835715._000007.LGNTuple.root",
//        },
//        {
//            "/home/common/Za/NTuples/Background/user.kschmied.364120.v11.Zee_140_280_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835749._000001.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364120.v11.Zee_140_280_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835749._000002.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364120.v11.Zee_140_280_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835749._000003.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364120.v11.Zee_140_280_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835749._000004.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364120.v11.Zee_140_280_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835749._000005.LGNTuple.root",
//        },
//        {
//            "/home/common/Za/NTuples/Background/user.kschmied.364115.v11.Zee_0_70_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835769._000001.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364115.v11.Zee_0_70_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835769._000002.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364115.v11.Zee_0_70_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835769._000003.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364115.v11.Zee_0_70_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835769._000004.LGNTuple.root",
//        },
//        {
//            "/home/common/Za/NTuples/Background/user.kschmied.364118.v11.Zee_70_140_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835773._000001.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364118.v11.Zee_70_140_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835773._000002.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364118.v11.Zee_70_140_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835773._000003.LGNTuple.root",
//        },
//        {
//            "/home/common/Za/NTuples/Background/user.kschmied.364121.v11.Zee_140_280_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835802._000001.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364121.v11.Zee_140_280_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835802._000002.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364121.v11.Zee_140_280_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835802._000003.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364121.v11.Zee_140_280_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835802._000004.LGNTuple.root",
//        },
//        {
//            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000001.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000002.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000003.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000004.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000005.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000006.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000007.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000008.LGNTuple.root",
//        },
//        {
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000001.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000002.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000003.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000004.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000005.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000006.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000007.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000008.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000009.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000010.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000011.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000012.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000013.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000014.LGNTuple.root",
//        },
//        {
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000001.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000002.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000003.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000004.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000005.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000006.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000007.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000008.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000009.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000010.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000011.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000012.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000013.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000014.LGNTuple.root",
//        },
//    };
//
//    std::array<double,3> SFs = {((139e15)*(.871e-12)),((139e15)*(.199e-12)), ((139e15)*(.0345e-15))}; //numerators for Z-gamma bkg
//
//    std::array<double,9> JetNumeratorSFs = {((139e15)*(1.9828e-9)*(0.821204)),((139e15)*(110.64e-12)*(0.69275)),((139e15)*(40.645e-12)*(0.615906)),((139e15)*(1.9817e-9)*(0.1136684)),((139e15)*(110.47e-12)*(0.1912956)),((139e15)*(40.674e-12)*(0.2326772)),((139e15)*(1.9819e-9)*(0.0656969)),((139e15)*(110.53e-12)*(0.1158741)),((139e15)*(40.68e-12)*(0.1535215))}; //numerators for jet bkg
//
//    std::vector<const char*> prefixes = {"Sig m_{A} = 5 GeV", "Sig m_{A} = 1 GeV", "pty2_9_17", "pty_17_myy_0_80", "pty_17_myy_80", "data", "Zee_lightJet_0-70", "Zee_lightJet_70-140", "Zee_lightJet_140-280", "Zee_cJet_0-70", "Zee_cJet_70-140", "Zee_cJet_140-280", "Zee_bJet_0-70", "Zee_bJet_70-140", "Zee_bJet_140-280"};
//    std::vector<EColor> colors = {kBlack, kMagenta, kBlue, kRed, kViolet, kGreen};
//    std::vector<EColor> Jetscolors = {kCyan, kOrange, kGreen, kYellow, kPink, kGray, kBlack, kSpring, kAzure};
//
//    std::vector<ROOT::RDF::RResultHandle> Nodes;
//    int count = 0;
//    for (auto& i: input_filenames)
//    {
//        SchottDataFrame df(MakeRDF(i, 24));
//
//        auto EventWeight = df.Define("EventWeight",
//        [](const RVec<float>& ei_event_weights_generator)
//        {
//            return  ((ei_event_weights_generator[0]) ? 1 / ei_event_weights_generator[0] : 1);
//
//        }, {"ei_event_weights_generator"});
//
//
//        auto two_leptons = EventWeight.Filter(
//        [](RVec<Muon>& muons, RVec<Electron> electrons)
//        {
//            electrons.erase(std::remove_if(electrons.begin(),electrons.end(),
//            [](Electron& ep)
//            {
//                return (!((ep.electron_pt/1e3 > 20) && (abs(ep.electron_eta) < 2.37) &&
//                          (!((1.37 < abs(ep.electron_eta)) && (abs(ep.electron_eta) < 1.52)))
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
//                return (!((ep.electron_pt/1e3 > 20) && (abs(ep.electron_eta) < 2.37) &&
//                (!((1.37 < abs(ep.electron_eta)) && (abs(ep.electron_eta) < 1.52)))
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
//              return ((abs(x.photon_eta) >= 2.37) || (abs(x.photon_eta) > 1.37 && abs(x.photon_eta) < 1.52) || (!x.photon_id_loose));
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
//                if (i==0 || abs(1-X) < abs(1-best_X))
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
//        .Define("dilepton_mergedPhoton_deltaR",[&](RVec<Electron>& di_electrons, Photon& merged_photon)
//        {
//            auto dilepton_four_momentum = di_electrons[0].Vector() + di_electrons[1].Vector();
//            return DeltaR(dilepton_four_momentum, merged_photon.Vector());
//
//        }, {"di_electrons", "merged_photon"})
//        .Define("reconstructed_mass",[&](RVec<Electron>& di_electrons, Photon& merged_photon)
//        {
//            auto four_momentum = di_electrons[0].Vector() + di_electrons[1].Vector();
//
//            return (four_momentum + merged_photon.Vector()).M()/1e3;
//
//        }, {"di_electrons", "merged_photon"})
//        .Define("totEventWeight", [](RVec<float> photon_id_eff, RVec<float> photon_iso_eff, RVec<float> photon_trg_eff/*, RVec<float> ei_event_weights_generator*/)
//        {
//            auto ResizeVal = std::max({photon_id_eff.size(), photon_iso_eff.size(), photon_trg_eff.size()});
//            photon_id_eff.resize(ResizeVal,1);
//            photon_iso_eff.resize(ResizeVal,1);
//            photon_trg_eff.resize(ResizeVal,1);
//
//            return photon_id_eff*photon_iso_eff*photon_trg_eff;//*ei_event_weights_generator[0];
//
//        }, {"photon_id_eff", "photon_iso_eff", "photon_trg_eff",/* "ei_event_weights_generator"*/});
//
//        auto SB = dilepton_and_photon.Filter(
//        [](RVec<float>& Eratio, double reconstructed_mass)
//        {
//            return ((!Any(Eratio <= 0.8)) && ((reconstructed_mass <= 110) || (reconstructed_mass >= 130)));
//        }, {"photon_shower_shape_e_ratio", "reconstructed_mass"});
//
//        auto SR = dilepton_and_photon.Filter(
//        [](RVec<float>& Eratio, double reconstructed_mass)
//        {
//            return ((!Any(Eratio <= 0.8)) && ((reconstructed_mass > 110) && (reconstructed_mass < 130)));
//        }, {"photon_shower_shape_e_ratio", "reconstructed_mass"});
//
//        if ((count >= 2 && count <= 4) || (count >= 6)) //background: Z-gamma and Z-jets
//        {
//            Nodes.push_back(SB.Sum<RVec<float>>("totEventWeight"));
//            Nodes.push_back(SR.Sum<RVec<float>>("totEventWeight"));
//            Nodes.push_back(EventWeight.Sum<float>("EventWeight"));
//        }
//
//        Nodes.push_back(SB.Histo1D<double>({prefixes[count], prefixes[count], 100u, 0, 4.7}, "dilepton_mergedPhoton_deltaR", "totEventWeight"));
//        Nodes.push_back(SR.Histo1D<double>({prefixes[count], prefixes[count++], 100u, 0, 4.7}, "dilepton_mergedPhoton_deltaR", "totEventWeight"));
//    }
//
////                0   1       //ma5
////                2   3       //ma1
////    4   5   6   7   8       //Z-gamma
////    9   10  11  12  13      //Z-gamma
////    14  15  16  17  18      //Z-gamma
////                19  20      //data
////    21  22  23  24  25      //Z-jets
////    26  27  28  29  30      //Z-jets
////    31  32  33  34  35      //Z-jets
////    36  37  38  39  40      //Z-jets
////    41  42  43  44  45      //Z-jets
////    46  47  48  49  50      //Z-jets
////    51  52  53  54  55      //Z-jets
////    56  57  58  59  60      //Z-jets
////    61  62  63  64  65      //Z-jets
//
//    ROOT::RDF::RunGraphs(Nodes); // running all computation nodes concurrently
//
//    count = 0;
//    TCanvas* c1 = new TCanvas();
//    TLegend* legend = new TLegend(0.675, 0.5, 0.875, 0.7);
//    double factor;
//    for (auto& i: {0,2}) //sideband signal
//    {
//        Nodes[i].GetResultPtr<TH1D>()->SetLineColor(colors[count++]);
//        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "l");
//
//        if (i == 0)
//        {
//            Nodes[i].GetResultPtr<TH1D>()->SetTitle(";#DeltaR (ll#gamma);Events");
//            Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);
//            Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->SetTitleOffset(1.2);
//            Nodes[i].GetResultPtr<TH1D>()->SetAxisRange(0., 10, "Y");
//            Nodes[i].GetResultPtr<TH1D>()->Draw("HIST");
//        }
//        else
//        {
//            Nodes[i].GetResultPtr<TH1D>()->Draw("HISTsame");
//            gPad->Modified(); gPad->Update();
//        }
//    }
//
//    gStyle->SetOptStat(0);
//    TLatex Tl;
//    Tl.SetTextSize(0.03);
//    Tl.DrawLatexNDC(0.6, 0.85, "#it{ATLAS} Internal");
//    Tl.DrawLatexNDC(0.6, 0.76,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
//    legend->SetBorderSize(0);
//    legend->Draw();
//    c1->SaveAs("Fig59D.pdf");
//
//    factor = 0;
//    count = 0;
//    for (auto& i: {4,9,14}) //Z-gamma sideband region
//    {
//        factor += (*Nodes[i].GetResultPtr<float>())*(SFs[count++]/ *Nodes[i+2].GetResultPtr<float>());
//    }
//
//    for (int i = 21, j = 0; ( i <= 61 && j <= 8); i+=5, j++) //Z-jets sideband region
//    {
//        factor += (*Nodes[i].GetResultPtr<float>())*(JetNumeratorSFs[j]/ *Nodes[i+2].GetResultPtr<float>());
//    }
//
//    count = 0;
//    auto hs = new THStack("hs1","");
//    c1 = new TCanvas();
//    legend = new TLegend(0.71, 0.435, 0.88, 0.785);
//    for (auto& i: {7,12,17,19}) //Z-gamma & data sideband region
//    {
//        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0 && i != 19)
//        {
//            Nodes[i].GetResultPtr<TH1D>()->Scale(SFs[count]/ *Nodes[i-1].GetResultPtr<float>());
//        }
//        Nodes[i].GetResultPtr<TH1D>()->SetFillColor(colors[2+count++]);
//        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "f");
//
//        if (i != 19)
//        {
//            hs->Add(&*Nodes[i].GetResultPtr<TH1D>());
//        }
//    }
//
//    for (int i = 24, j = 0; (i <= 64 && j <= 8); i+=5, j++) //Z-jets sideband region
//    {
//        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0)
//        {
//            Nodes[i].GetResultPtr<TH1D>()->Scale(JetNumeratorSFs[j]/ *Nodes[i-1].GetResultPtr<float>());
//        }
//        Nodes[i].GetResultPtr<TH1D>()->SetFillColor(Jetscolors[j]);
//        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "f");
//
//        hs->Add(&*Nodes[i].GetResultPtr<TH1D>());
//    }
//
//    hs->Draw("HIST");
//    Nodes[19].GetResultPtr<TH1D>()->Draw("HISTsame");
//    hs->SetTitle(";#DeltaR (ll#gamma);Events");
//    hs->SetMaximum(2.5e4);
//    hs->GetYaxis()->CenterTitle(true);
//    hs->GetXaxis()->SetTitleOffset(1.2);
//    hs->GetYaxis()->SetTitleOffset(1.35);
//
//    gStyle->SetOptStat(0);
//    Tl.SetTextSize(0.03);
//    Tl.DrawLatexNDC(0.2, 0.82, "#it{ATLAS} Internal");
//    Tl.DrawLatexNDC(0.2, 0.72,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
//    legend->SetBorderSize(0);
//    legend->Draw();
//    c1->SaveAs("Fig59A.pdf");
//
//    factor = 0;
//    count = 0;
//    double sig_factor = Nodes[1].GetResultPtr<TH1D>()->Integral();
//    std::vector<double> backScalings;
//
//    for (auto& i: {5,10,15}) //Z-gamma SR
//    {
//        backScalings.push_back((*Nodes[i].GetResultPtr<float>())*(SFs[count]/ *Nodes[i+1].GetResultPtr<float>()));
//
//        factor += (*Nodes[i].GetResultPtr<float>())*(SFs[count++]/ *Nodes[i+1].GetResultPtr<float>());
//    }
//
//    for (int i = 22, j = 0; ( i <= 62 && j <= 8); i+=5, j++) //Z-jets SR
//    {
//        backScalings.push_back((*Nodes[i].GetResultPtr<float>())*(JetNumeratorSFs[j]/ *Nodes[i+1].GetResultPtr<float>()));
//
//        factor += (*Nodes[i].GetResultPtr<float>())*(JetNumeratorSFs[j]/ *Nodes[i+1].GetResultPtr<float>());
//    }
//
//    for (double& i: backScalings)
//    {
//        i = (i/factor)*sig_factor;
//    }
//
//    count=2;
//    c1 = new TCanvas();
//    legend = new TLegend(0.66, 0.4, 0.86, 0.8);
//    hs = new THStack("hs2","");
//    for (auto& i: {8,13,18,20}) //Z-gamma & data signal region
//    {
//        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0 && i != 20)
//        {
//            Nodes[i].GetResultPtr<TH1D>()->Scale(SFs[count-2]/ *Nodes[i-2].GetResultPtr<float>());
//        }
//        if (i != 20)
//        {
//            hs->Add(static_cast<TH1D*>(Nodes[i].GetResultPtr<TH1D>()->Clone()));
//            Nodes[i].GetResultPtr<TH1D>()->SetFillColor(colors[count++]);
//            legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "f");
//        }
//    }
//
//    for (int i = 25, j = 0; (i <= 65 && j <= 8); i+=5, j++) //Z-jets SR
//    {
//        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0)
//        {
//            Nodes[i].GetResultPtr<TH1D>()->Scale(JetNumeratorSFs[j]/ *Nodes[i-2].GetResultPtr<float>());
//        }
//        Nodes[i].GetResultPtr<TH1D>()->SetFillColor(Jetscolors[j]);
//        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "f");
//
//        hs->Add(static_cast<TH1D*>(Nodes[i].GetResultPtr<TH1D>()->Clone()));
//    }
//
//    hs->Draw("HIST");
//    count=0;
//    for (auto& i: {1,3}) // signal SR
//    {
//        Nodes[i].GetResultPtr<TH1D>()->SetLineColor(colors[count]);
//        Nodes[i].GetResultPtr<TH1D>()->SetLineWidth(2);
//        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), prefixes[count++], "l");
//
//        Nodes[i].GetResultPtr<TH1D>()->SetTitle(";#DeltaR (ll#gamma);Events");
//        Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);
//        Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->SetTitleOffset(1.2);
//        Nodes[i].GetResultPtr<TH1D>()->SetAxisRange(0., 68, "Y");
//        Nodes[i].GetResultPtr<TH1D>()->Draw("HISTsame");
//        gPad->Modified(); gPad->Update();
//    }
//    hs->SetMaximum(5e4);
//    hs->SetTitle(";#DeltaR (ll#gamma);Events");
//    hs->GetYaxis()->CenterTitle(true);
//    hs->GetXaxis()->SetTitleOffset(1.2);
//    hs->GetYaxis()->SetTitleOffset(1.35);
//
//    gStyle->SetOptStat(0);
//    Tl.SetTextSize(0.03);
//    Tl.DrawLatexNDC(0.35, 0.85, "#it{ATLAS} Internal");
//    Tl.DrawLatexNDC(0.35, 0.75,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
//    legend->SetBorderSize(0);
//    legend->Draw();
//    c1->SaveAs("Fig59B.pdf");
//
//    count=2;
//    c1 = new TCanvas();
//    legend = new TLegend(0.7, 0.35, 0.875, 0.75);
//    count=2;
//    hs = new THStack("hs3","");
//    int temp = 0;
//    for (auto& i: {8,13,18,20}) //background and data SR
//    {
//        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0 && i != 20)
//        {
//            Nodes[i].GetResultPtr<TH1D>()->Scale(backScalings[temp++]/Nodes[i].GetResultPtr<TH1D>()->Integral());
//        }
//
//        if (i != 20)
//        {
//            hs->Add(&*Nodes[i].GetResultPtr<TH1D>());
//            Nodes[i].GetResultPtr<TH1D>()->SetFillColor(colors[count++]);
//            legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "f");
//        }
//    }
//    temp = 3;
//    for (int i = 25, j = 0; (i <= 65 && j <= 8); i+=5, j++) //Z-jets SR
//    {
//        if (Nodes[i].GetResultPtr<TH1D>()->Integral() != 0)
//        {
//            Nodes[i].GetResultPtr<TH1D>()->Scale(backScalings[temp++]/Nodes[i].GetResultPtr<TH1D>()->Integral());
//        }
//        Nodes[i].GetResultPtr<TH1D>()->SetFillColor(Jetscolors[j]);
//        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), Nodes[i].GetResultPtr<TH1D>()->GetTitle(), "f");
//
//        hs->Add(&*Nodes[i].GetResultPtr<TH1D>());
//    }
//
//    hs->Draw("HIST");
//    count=0;
//    for (auto& i: {1,3}) //signal SR
//    {
//        Nodes[i].GetResultPtr<TH1D>()->SetLineColor(colors[count++]);
//        Nodes[i].GetResultPtr<TH1D>()->SetLineWidth(2);
//        legend->AddEntry(&(*Nodes[i].GetResultPtr<TH1D>()), prefixes[count-1], "l");
//
//        Nodes[i].GetResultPtr<TH1D>()->Scale(sig_factor/Nodes[i].GetResultPtr<TH1D>()->Integral());
//        Nodes[i].GetResultPtr<TH1D>()->SetTitle(";#DeltaR (ll#gamma);Events");
//        Nodes[i].GetResultPtr<TH1D>()->GetYaxis()->CenterTitle(true);
//        Nodes[i].GetResultPtr<TH1D>()->GetXaxis()->SetTitleOffset(1.2);
//        Nodes[i].GetResultPtr<TH1D>()->Draw("HISTsame");
//        gPad->Modified(); gPad->Update();
//    }
//    hs->SetMaximum(68);
//    hs->SetTitle(";#DeltaR (ll#gamma);Events");
//    hs->GetYaxis()->CenterTitle(true);
//    hs->GetXaxis()->SetTitleOffset(1.2);
//
//    gStyle->SetOptStat(0);
//    Tl.SetTextSize(0.03);
//    Tl.DrawLatexNDC(0.2, 0.825, "#it{ATLAS} Internal");
//    Tl.DrawLatexNDC(0.2, 0.725,"#sqrt{s} = 13 TeV  #int L #bullet dt = 139 fb^{-1}");
//    legend->SetBorderSize(0);
//    legend->Draw();
//    c1->SaveAs("Fig59C.pdf");
//}
/*
void Table9()
{
    std::vector<std::vector<std::string>> input_filenames =
    {
        //Z gamma background
        {
            "/home/common/Za/NTuples/Background/user.kschmied.364860.eegammagamma_pty2_9_17.deriv.DAOD_STDM3.e7057_s3126_r10724_p4092_LGNTuple.root/user.kschmied.31617070._000001.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364861.eegammagamma_pty_17_myy_0_80.deriv.DAOD_HIGG1D2.e7057_s3126_r10724_p4062_LGNTuple.root/user.kschmied.31617064._000001.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364862.eegammagamma_pty_17_myy_80.deriv.DAOD_HIGG1D2.e7057_s3126_r10724_p4204_LGNTuple.root/user.kschmied.31660711._000001.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364862.eegammagamma_pty_17_myy_80.deriv.DAOD_HIGG1D2.e7057_s3126_r10724_p4204_LGNTuple.root/user.kschmied.31660711._000002.LGNTuple.root",
        },
        //Jets
        {
            "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000001.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000002.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000003.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000004.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000005.LGNTuple.root",

            "/home/common/Za/NTuples/Background/user.kschmied.364117.v11.Zee_70_140_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835715._000001.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364117.v11.Zee_70_140_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835715._000002.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364117.v11.Zee_70_140_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835715._000003.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364117.v11.Zee_70_140_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835715._000004.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364117.v11.Zee_70_140_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835715._000006.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364117.v11.Zee_70_140_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835715._000007.LGNTuple.root",

            "/home/common/Za/NTuples/Background/user.kschmied.364120.v11.Zee_140_280_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835749._000001.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364120.v11.Zee_140_280_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835749._000002.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364120.v11.Zee_140_280_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835749._000003.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364120.v11.Zee_140_280_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835749._000004.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364120.v11.Zee_140_280_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835749._000005.LGNTuple.root",


            "/home/common/Za/NTuples/Background/user.kschmied.364115.v11.Zee_0_70_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835769._000001.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364115.v11.Zee_0_70_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835769._000002.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364115.v11.Zee_0_70_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835769._000003.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364115.v11.Zee_0_70_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835769._000004.LGNTuple.root",


            "/home/common/Za/NTuples/Background/user.kschmied.364118.v11.Zee_70_140_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835773._000001.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364118.v11.Zee_70_140_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835773._000002.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364118.v11.Zee_70_140_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835773._000003.LGNTuple.root",


            "/home/common/Za/NTuples/Background/user.kschmied.364121.v11.Zee_140_280_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835802._000001.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364121.v11.Zee_140_280_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835802._000002.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364121.v11.Zee_140_280_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835802._000003.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364121.v11.Zee_140_280_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835802._000004.LGNTuple.root",


            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000001.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000002.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000003.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000004.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000005.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000006.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000007.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000008.LGNTuple.root",


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
    
    std::vector<ROOT::RDF::RResultHandle> Nodes;
    
    for (auto& i: input_filenames)
    {
        SchottDataFrame df(MakeRDF(i, 24));
        
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
        [&](RVec<TruthParticle> truth_particles)
        {
            truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
            [](TruthParticle& x)
            {
                return ((abs(x.mc_pdg_id) != 22) || (abs(x.mc_eta) >= 2.37));

            }), truth_particles.end());
            
            return truth_particles;
            
        }, {"truth_particles"});
        
        auto diphotons = photon_passes_cuts.Define("chosen_two",
        [](RVec<TruthParticle>& reco_photons_matched)
        {
            RVec<TruthParticle> x;
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
                    pt1 = reco_photons_matched[combs[0][i]].mc_pt;
                    pt2 = reco_photons_matched[combs[1][i]].mc_pt;
                    chosen_delta_r = delta_r;
                    x = {reco_photons_matched[combs[0][i]], reco_photons_matched[combs[1][i]]};
                }
            }
            if (pt1 > 10e3 && pt2 > 10e3 && best_X > 0.96 && best_X < 1.2 && chosen_delta_r < 1.5)
            {
                return x;
            }
            x.clear(); //it doesn't pass, so return an empty vector
            return x;
        }, {"photons_pass_cuts"}).Filter(
        [&](RVec<TruthParticle>& reco_photons_matched)
        {
            return reco_photons_matched.size() == 2;
        }, {"chosen_two"}).Define("leading_photon",
        [&](RVec<TruthParticle>& reco_photons_matched)
        {
            return (reco_photons_matched[0].mc_pt > reco_photons_matched[1].mc_pt)
                    ? reco_photons_matched[0]
                    : reco_photons_matched[1];
        },{"chosen_two"}).Define("subleading_photon",
        [&](RVec<TruthParticle>& reco_photons_matched)
        {
            return (reco_photons_matched[0].mc_pt < reco_photons_matched[1].mc_pt)
                    ? reco_photons_matched[0]
                    : reco_photons_matched[1];
        },{"chosen_two"}).Define("leading_photon_pdg_id_origin",
        [](TruthParticle& leading_photon, RVec<TruthParticle>& truth_particles)
        {
            TruthParticle origin = leading_photon;
            int origin_id = 0;
            bool found;
            
            do
            {
                found = false;
                for (auto& i: truth_particles)
                {
                    if (origin.mc_parent_barcode == i.mc_barcode)
                    {
                        origin = i;
                        origin_id = i.mc_pdg_id;
                        found = true;
                    }
                }
            } while (found);
            
            return origin_id;
        }, {"leading_photon", "truth_particles"}).Define("subleading_photon_pdg_id_origin",
        [](TruthParticle& subleading_photon, RVec<TruthParticle>& truth_particles)
        {
            TruthParticle origin = subleading_photon;
            int origin_id = 0;
            bool found;

            do
            {
             found = false;
             for (auto& i: truth_particles)
             {
                 if (origin.mc_parent_barcode == i.mc_barcode)
                 {
                     origin = i;
                     origin_id = i.mc_pdg_id;
                     found = true;
                 }
             }
            } while (found);

            return origin_id;
        }, {"subleading_photon", "truth_particles"});
        
        Nodes.push_back(diphotons.Take<int, RVec<int>>("leading_photon_pdg_id_origin"));
        Nodes.push_back(diphotons.Take<int, RVec<int>>("subleading_photon_pdg_id_origin"));
        Nodes.push_back(diphotons.Count());
    }
    
    ROOT::RDF::RunGraphs(Nodes); // running all computation nodes concurrently
    
    std::unordered_map<int,int>
    leading_id_freqs_Z_gamma, subleading_id_freqs_Z_gamma,
    leading_id_freqs_Z_jets, subleading_id_freqs_Z_jets,
    leading_id_freqs_total, subleading_id_freqs_total;
    
    double total_Z_gamma = *Nodes[2].GetResultPtr<ULong64_t>();
    double total_Z_jets = *Nodes[5].GetResultPtr<ULong64_t>();
    double total = total_Z_gamma + total_Z_jets;
    
    for (auto& i: *Nodes[0].GetResultPtr<RVec<int>>())
    {
        leading_id_freqs_Z_gamma[i]++;
        leading_id_freqs_total[i]++;
    }
    for (auto& i: *Nodes[1].GetResultPtr<RVec<int>>())
    {
        subleading_id_freqs_Z_gamma[i]++;
        subleading_id_freqs_total[i]++;
    }
    for (auto& i: *Nodes[3].GetResultPtr<RVec<int>>())
    {
        leading_id_freqs_Z_jets[i]++;
        leading_id_freqs_total[i]++;
    }
    for (auto& i: *Nodes[4].GetResultPtr<RVec<int>>())
    {
        subleading_id_freqs_Z_jets[i]++;
        subleading_id_freqs_total[i]++;
    }
    
    std::cout << R"--(\section*{Table 9}\flushleft)--" << '\n';
        std::cout << R"--(\hspace{-3cm}\scalebox{0.65}{
                          \setlength\extrarowheight{2pt}
                          \renewcommand{\arraystretch}{1.5})--" << '\n';
        std::cout << R"--(\begin{tabular}{|c|c|c|c|c|c|c|})--" << '\n';
        std::cout << R"--(\hline)--" << '\n';
        std::cout << R"--(\multicolumn{7}{|c|}{\parbox{\linewidth}{\centering Resolved Category \\ Background photon truth origin fractions}}\\[5 pt] \hline)--" << '\n';
        std::cout << R"--( {} & \multicolumn{2}{|c|}{\parbox{7.5cm}{\centering $Z\gamma$ \\ (\%)}} & \multicolumn{2}{|c|}{\parbox{7.5cm}{\centering $Z$ jets \\ (\%)}} & \multicolumn{2}{|c|}{\parbox{7.5cm}{\centering Total bkg \\ (\%)}} \\[5 pt]
            \hline)--" << '\n';
        std::cout << R"--( {Pdg Id} & leading photon & subleading photon & leading photon & subleading photon & leading photon & subleading photon \\[5 pt]
            \hline)--" << '\n';
    
    
    for (auto& i: leading_id_freqs_total)
    {
        std::cout << i.first << " & " << std::setprecision(2) << std::fixed
        << 100*(leading_id_freqs_Z_gamma[i.first]/total_Z_gamma)
        << " & " << std::setprecision(2) << std::fixed
        << 100*(subleading_id_freqs_Z_gamma[i.first]/total_Z_gamma)
        << " & " << std::setprecision(2) << std::fixed
        << 100*(leading_id_freqs_Z_jets[i.first]/total_Z_jets)
        << " & " << std::setprecision(2) << std::fixed
        << 100*(subleading_id_freqs_Z_jets[i.first]/total_Z_jets)
        << " & " << std::setprecision(2) << std::fixed
        << 100*(leading_id_freqs_total[i.first]/total)
        << " & " << std::setprecision(2) << std::fixed
        << 100*(subleading_id_freqs_total[i.first]/total)
        << R"--( \\ \hline)--" << '\n';
    }
    
    std::cout << R"--(\end{tabular}})--" << '\n';
    
    std::cout << "\n\n\n";
}
 

void Table10()
{
    std::vector<std::vector<std::string>> input_filenames =
    {
        //Z gamma background
        {
            "/home/common/Za/NTuples/Background/user.kschmied.364860.eegammagamma_pty2_9_17.deriv.DAOD_STDM3.e7057_s3126_r10724_p4092_LGNTuple.root/user.kschmied.31617070._000001.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364861.eegammagamma_pty_17_myy_0_80.deriv.DAOD_HIGG1D2.e7057_s3126_r10724_p4062_LGNTuple.root/user.kschmied.31617064._000001.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364862.eegammagamma_pty_17_myy_80.deriv.DAOD_HIGG1D2.e7057_s3126_r10724_p4204_LGNTuple.root/user.kschmied.31660711._000001.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364862.eegammagamma_pty_17_myy_80.deriv.DAOD_HIGG1D2.e7057_s3126_r10724_p4204_LGNTuple.root/user.kschmied.31660711._000002.LGNTuple.root",
        },
        //Jets
        {
            "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000001.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000002.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000003.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000004.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000005.LGNTuple.root",

            "/home/common/Za/NTuples/Background/user.kschmied.364117.v11.Zee_70_140_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835715._000001.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364117.v11.Zee_70_140_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835715._000002.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364117.v11.Zee_70_140_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835715._000003.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364117.v11.Zee_70_140_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835715._000004.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364117.v11.Zee_70_140_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835715._000006.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364117.v11.Zee_70_140_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835715._000007.LGNTuple.root",

            "/home/common/Za/NTuples/Background/user.kschmied.364120.v11.Zee_140_280_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835749._000001.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364120.v11.Zee_140_280_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835749._000002.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364120.v11.Zee_140_280_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835749._000003.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364120.v11.Zee_140_280_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835749._000004.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364120.v11.Zee_140_280_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835749._000005.LGNTuple.root",


            "/home/common/Za/NTuples/Background/user.kschmied.364115.v11.Zee_0_70_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835769._000001.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364115.v11.Zee_0_70_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835769._000002.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364115.v11.Zee_0_70_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835769._000003.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364115.v11.Zee_0_70_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835769._000004.LGNTuple.root",


            "/home/common/Za/NTuples/Background/user.kschmied.364118.v11.Zee_70_140_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835773._000001.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364118.v11.Zee_70_140_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835773._000002.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364118.v11.Zee_70_140_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835773._000003.LGNTuple.root",


            "/home/common/Za/NTuples/Background/user.kschmied.364121.v11.Zee_140_280_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835802._000001.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364121.v11.Zee_140_280_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835802._000002.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364121.v11.Zee_140_280_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835802._000003.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364121.v11.Zee_140_280_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835802._000004.LGNTuple.root",


            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000001.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000002.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000003.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000004.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000005.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000006.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000007.LGNTuple.root",
            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000008.LGNTuple.root",


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
    
    std::vector<ROOT::RDF::RResultHandle> Nodes;
    
    for (auto& i: input_filenames)
    {
        SchottDataFrame df(MakeRDF(i, 24));
        
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
        [&](RVec<TruthParticle> truth_particles)
        {
            truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
            [](TruthParticle& x)
            {
                return ((abs(x.mc_pdg_id) != 22) || (abs(x.mc_eta) >= 2.37));

            }), truth_particles.end());
            
            return truth_particles;
            
        }, {"truth_particles"});
        
        auto diphotons = photon_passes_cuts.Define("chosen_two",
        [](RVec<TruthParticle>& reco_photons_matched)
        {
            RVec<TruthParticle> x;
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
                    pt1 = reco_photons_matched[combs[0][i]].mc_pt;
                    pt2 = reco_photons_matched[combs[1][i]].mc_pt;
                    chosen_delta_r = delta_r;
                    x = {reco_photons_matched[combs[0][i]], reco_photons_matched[combs[1][i]]};
                }
            }
            if (pt1 > 10e3 && pt2 > 10e3 && best_X > 0.96 && best_X < 1.2 && chosen_delta_r < 1.5)
            {
                return x;
            }
            x.clear();
            return x; //it doesn't pass, so return an empty vector
        }, {"photons_pass_cuts"}).Filter(
        [&](RVec<TruthParticle>& reco_photons_matched)
        {
            return reco_photons_matched.size() == 2;
        }, {"chosen_two"}).Define("leading_photon",
        [&](RVec<TruthParticle>& reco_photons_matched)
        {
            return (reco_photons_matched[0].mc_pt > reco_photons_matched[1].mc_pt)
                    ? reco_photons_matched[0]
                    : reco_photons_matched[1];
        },{"chosen_two"}).Define("subleading_photon",
        [&](RVec<TruthParticle>& reco_photons_matched)
        {
            return (reco_photons_matched[0].mc_pt < reco_photons_matched[1].mc_pt)
                    ? reco_photons_matched[0]
                    : reco_photons_matched[1];
        },{"chosen_two"}).Define("photon_pdg_id_origin",
        [](TruthParticle& leading_photon, TruthParticle& subleading_photon, RVec<TruthParticle>& truth_particles)
        {
            TruthParticle leading_origin = leading_photon;
            TruthParticle subleading_origin = subleading_photon;
            int leading_origin_id = 0, subleading_origin_id = 0;
            bool leading_found, subleading_found;
            
            do
            {
                leading_found = false;
                subleading_found = false;
                for (auto& i: truth_particles)
                {
                    if (leading_origin.mc_parent_barcode == i.mc_barcode)
                    {
                        leading_origin = i;
                        leading_origin_id = i.mc_pdg_id;
                        leading_found = true;
                    }
                    
                    if (subleading_origin.mc_parent_barcode == i.mc_barcode)
                    {
                        subleading_origin = i;
                        subleading_origin_id = i.mc_pdg_id;
                        subleading_found = true;
                    }
                }
            } while (leading_found || subleading_found);
            
            return std::to_string(leading_origin_id)+"/"+std::to_string(subleading_origin_id);
        }, {"leading_photon", "subleading_photon", "truth_particles"});
        
        Nodes.push_back(diphotons.Take<std::string, RVec<std::string>>("photon_pdg_id_origin"));
        Nodes.push_back(diphotons.Count());
        
    }
    
    std::unordered_map<std::string,int>
    id_freqs_Z_gamma, id_freqs_Z_jets, id_freqs_Z_total;
        
    ROOT::RDF::RunGraphs(Nodes); // running all computation nodes concurrently
    
    double total_Z_gamma = *Nodes[1].GetResultPtr<ULong64_t>();
    double total_Z_jets = *Nodes[3].GetResultPtr<ULong64_t>();
    double total = total_Z_gamma + total_Z_jets;
    
    for (auto& i: *Nodes[0].GetResultPtr<RVec<std::string>>())
    {
        id_freqs_Z_gamma[i]++;
        id_freqs_Z_total[i]++;
    }
    
    for (auto& i: *Nodes[2].GetResultPtr<RVec<std::string>>())
    {
        id_freqs_Z_jets[i]++;
        id_freqs_Z_total[i]++;
    }
    
    std::cout << R"--(\section*{Table 10}\flushleft)--" << '\n';
        std::cout << R"--(\hspace{-3cm}\scalebox{0.65}{
                          \setlength\extrarowheight{2pt}
                          \renewcommand{\arraystretch}{1.5})--" << '\n';
        std::cout << R"--(\begin{tabular}{|c|c|c|c|})--" << '\n';
        std::cout << R"--(\hline)--" << '\n';
        std::cout << R"--(\multicolumn{4}{|c|}{\parbox{\linewidth}{\centering Resolved Category \\ Background photon truth origin fractions}}\\[5 pt] \hline)--" << '\n';
        std::cout << R"--( \multicolumn{1}{|c|}{\parbox{7.5cm}{\centering
            lead/sublead \; (Pdg Id) \\ photon pairs}} &
            \multicolumn{1}{|c|}{\parbox{4.5cm}{\centering $Z\gamma$ \\ (\%)}} &
            \multicolumn{1}{|c|}{\parbox{4.5cm}{\centering $Z$ jets \\ (\%)}} & \multicolumn{1}{|c|}{\parbox{4.5cm}{\centering Total bkg \\ (\%)}} \\[5 pt]
            \hline)--" << '\n';
    
    for (auto& i: id_freqs_Z_total)
    {
        std::cout << i.first << " & " << std::setprecision(2) << std::fixed
        << 100*(id_freqs_Z_gamma[i.first]/total_Z_gamma)
        << " & " << std::setprecision(2) << std::fixed
        << 100*(id_freqs_Z_jets[i.first]/total_Z_jets)
        << " & " << std::setprecision(2) << std::fixed
        << 100*(id_freqs_Z_total[i.first]/total)
        << R"--( \\ \hline)--" << '\n';
    }
    std::cout << R"--(\end{tabular}})--" << '\n';
    
    std::cout << "\n\n\n";
}
*/
//void Table16()
//{
//    std::vector<std::vector<std::string>> input_filenames =
//    {
//        //Z gamma background
//        {"/home/common/Za/NTuples/Background/user.kschmied.364860.eegammagamma_pty2_9_17.deriv.DAOD_STDM3.e7057_s3126_r10724_p4092_LGNTuple.root/user.kschmied.31617070._000001.LGNTuple.root"},
//        {"/home/common/Za/NTuples/Background/user.kschmied.364861.eegammagamma_pty_17_myy_0_80.deriv.DAOD_HIGG1D2.e7057_s3126_r10724_p4062_LGNTuple.root/user.kschmied.31617064._000001.LGNTuple.root"},
//        {"/home/common/Za/NTuples/Background/user.kschmied.364862.eegammagamma_pty_17_myy_80.deriv.DAOD_HIGG1D2.e7057_s3126_r10724_p4204_LGNTuple.root/user.kschmied.31660711._000001.LGNTuple.root", "/home/common/Za/NTuples/Background/user.kschmied.364862.eegammagamma_pty_17_myy_80.deriv.DAOD_HIGG1D2.e7057_s3126_r10724_p4204_LGNTuple.root/user.kschmied.31660711._000002.LGNTuple.root"},
//        //Data
//        {"/home/common/Za/NTuples/Ntuple_data_test.root"},
//        //Signal
//        {"/home/common/Za/NTuples/Ntuple_MC_Za_m5p0_v4.root"},
//        {"/home/common/Haa/ntuples/Za/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root"},
//        //Jets
//        {
//            "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000001.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000002.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000003.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000004.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000005.LGNTuple.root",
//        },
//        {
//            "/home/common/Za/NTuples/Background/user.kschmied.364117.v11.Zee_70_140_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835715._000001.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364117.v11.Zee_70_140_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835715._000002.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364117.v11.Zee_70_140_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835715._000003.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364117.v11.Zee_70_140_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835715._000004.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364117.v11.Zee_70_140_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835715._000006.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364117.v11.Zee_70_140_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835715._000007.LGNTuple.root",
//        },
//        {
//            "/home/common/Za/NTuples/Background/user.kschmied.364120.v11.Zee_140_280_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835749._000001.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364120.v11.Zee_140_280_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835749._000002.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364120.v11.Zee_140_280_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835749._000003.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364120.v11.Zee_140_280_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835749._000004.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364120.v11.Zee_140_280_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835749._000005.LGNTuple.root",
//        },
//        {
//            "/home/common/Za/NTuples/Background/user.kschmied.364115.v11.Zee_0_70_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835769._000001.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364115.v11.Zee_0_70_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835769._000002.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364115.v11.Zee_0_70_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835769._000003.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364115.v11.Zee_0_70_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835769._000004.LGNTuple.root",
//        },
//        {
//            "/home/common/Za/NTuples/Background/user.kschmied.364118.v11.Zee_70_140_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835773._000001.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364118.v11.Zee_70_140_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835773._000002.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364118.v11.Zee_70_140_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835773._000003.LGNTuple.root",
//        },
//        {
//            "/home/common/Za/NTuples/Background/user.kschmied.364121.v11.Zee_140_280_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835802._000001.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364121.v11.Zee_140_280_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835802._000002.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364121.v11.Zee_140_280_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835802._000003.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364121.v11.Zee_140_280_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835802._000004.LGNTuple.root",
//        },
//        {
//            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000001.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000002.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000003.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000004.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000005.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000006.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000007.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000008.LGNTuple.root",
//        },
//        {
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000001.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000002.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000003.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000004.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000005.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000006.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000007.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000008.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000009.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000010.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000011.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000012.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000013.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000014.LGNTuple.root",
//        },
//        {
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000001.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000002.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000003.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000004.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000005.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000006.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000007.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000008.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000009.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000010.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000011.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000012.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000013.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000014.LGNTuple.root",
//        },
//    };
//
//    std::array<double,9> JetNumeratorSFs = {((139e15)*(1.9828e-9)*(0.821204)),((139e15)*(110.64e-12)*(0.69275)),((139e15)*(40.645e-12)*(0.615906)),((139e15)*(1.9817e-9)*(0.1136684)),((139e15)*(110.47e-12)*(0.1912956)),((139e15)*(40.674e-12)*(0.2326772)),((139e15)*(1.9819e-9)*(0.0656969)),((139e15)*(110.53e-12)*(0.1158741)),((139e15)*(40.68e-12)*(0.1535215))};
//    std::array<double,3> SFs = {((139e15)*(.871e-12)),((139e15)*(.199e-12)), ((139e15)*(.0345e-15))}; //numerators for Z-gamma bkg
//
//    std::vector<std::string> prefixes = { R"--(pty2\_9\_17)--", R"--(pty\_17\_myy\_0\_80)--", R"--(pty\_17\_myy\_80)--", "data", R"--($\text{Sig } m_{A}$ = 5 GeV)--", R"--($\text{Sig } m_{A}$ = 1 GeV)--", R"--(Zee\_lightJet\_0-70)--", R"--(Zee\_lightJet\_70-140)--", R"--(Zee\_lightJet\_140-280)--", R"--(Zee\_cJet\_0-70)--", R"--(Zee\_cJet\_70-140)--", R"--(Zee\_cJet\_140-280)--", R"--(Zee\_bJet\_0-70)--", R"--(Zee\_bJet\_70-140)--", R"--(Zee\_bJet\_140-280)--"};
//
//    std::vector<ROOT::RDF::RResultHandle> Totals;
//
//    for (auto& i: input_filenames)
//    {
//        SchottDataFrame df(MakeRDF(i,24));
////        std::cout << *df.Count() << '\n';
//
//        auto two_leptons = df.Filter(
//        [](RVec<Muon>& muons, RVec<Electron> electrons)
//        {
//            electrons.erase(std::remove_if(electrons.begin(),electrons.end(),
//            [](Electron& ep)
//            {
//                return (!((ep.electron_pt/1e3 > 20) && (abs(ep.electron_eta) < 2.37) &&
//                          (!((1.37 < abs(ep.electron_eta)) && (abs(ep.electron_eta) < 1.52)))
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
//                return (!((ep.electron_pt/1e3 > 20) && (abs(ep.electron_eta) < 2.37) &&
//                (!((1.37 < abs(ep.electron_eta)) && (abs(ep.electron_eta) < 1.52)))
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
//              return ((abs(x.photon_eta) >= 2.37) || (abs(x.photon_eta) > 1.37 && abs(x.photon_eta) < 1.52) || (!x.photon_id_loose));
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
//                if (i==0 || abs(1-X) < abs(1-best_X))
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
////    0    1    2    3    4    5      Z-gamma
////    6    7    8    9    10   11     Z-gamma
////    12   13   14   15   16   17     Z-gamma
////    18   19   20   21   22   23     data
////    24   25   26   27   28   29     signal
////    30   31   32   33   34   35     signal
////    36   37   38   39   40   41     Z-jets
////    42   43   44   45   46   47     Z-jets
////    48   49   50   51   52   53     Z-jets
////    54   55   56   57   58   59     Z-jets
////    60   61   62   63   64   65     Z-jets
////    66   67   68   69   70   71     Z-jets
////    72   73   74   75   76   77     Z-jets
////    78   79   80   81   82   83     Z-jets
////    84   85   86   87   88   89     Z-jets
//
//
//    ROOT::RDF::RunGraphs(Totals); // running all computation nodes concurrently
//
//    std::cout << R"--(\section*{Table 16})--" << '\n';
//    std::cout << R"--(\hspace{-3cm}\scalebox{0.65}{)--" << '\n';
//    std::cout << R"--(\begin{tabular}{|c|c|c|c|c|c|})--" << '\n';
//    std::cout << R"--(\hline)--" << '\n';
//    std::cout << R"--({} & Full Reg & pSB & pSR & SB & SR \\ \hline)--" << '\n';
//
//    double total_Z_gamma[5] = {0}, total_Z_jets[5] = {0}, total[5] = {0};
//    double total_Z_gammaStatUnc[5] = {0}, total_Z_jetsStatUnc[5] = {0}, totalStatUnc[5] = {0};
//
//    for (int i = 0, j = 0; (i <= 84 && j <= 14); i += 6, j++)
//    {
//        std::cout << prefixes[j];
//        if (j >= 0 && j <= 2)
//        {
//            std::cout
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
//        else if (j >= 6)
//        {
//            std::cout
//            << " & " << std::setprecision(2) << std::fixed << *Totals[i+1].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j-6] / *Totals[i].GetResultPtr<ULong64_t>())
//            << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(*Totals[i+1].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[j-6] / *Totals[i].GetResultPtr<ULong64_t>())
//            << " & " << std::setprecision(2) << std::fixed << *Totals[i+2].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j-6] / *Totals[i].GetResultPtr<ULong64_t>())
//            << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(*Totals[i+2].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[j-6] / *Totals[i].GetResultPtr<ULong64_t>())
//            << " & " << std::setprecision(2) << std::fixed << *Totals[i+3].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j-6] / *Totals[i].GetResultPtr<ULong64_t>())
//            << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(*Totals[i+3].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[j-6] / *Totals[i].GetResultPtr<ULong64_t>())
//            << " & " << std::setprecision(2) << std::fixed << *Totals[i+4].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j-6] / *Totals[i].GetResultPtr<ULong64_t>())
//            << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(*Totals[i+4].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[j-6] / *Totals[i].GetResultPtr<ULong64_t>())
//            << " & " << std::setprecision(2) << std::fixed << *Totals[i+5].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j-6] / *Totals[i].GetResultPtr<ULong64_t>())
//            << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(*Totals[i+5].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[j-6] / *Totals[i].GetResultPtr<ULong64_t>())
//            << R"--( \\ \hline)--" << '\n';
//        }
//
//        else
//        {
//            std::cout
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
//        else if (i >= 36 && i <= 84) //Z-jets
//        {
//            for (int k = 0; k < 5; k++)
//            {
//                total_Z_jets[k] += *Totals[i+k+1].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i/6 - 6] / *Totals[i].GetResultPtr<ULong64_t>());
//                total[k] += *Totals[i+k+1].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i/6 - 6] / *Totals[i].GetResultPtr<ULong64_t>());
//
//                total_Z_jetsStatUnc[k] += pow(sqrt(*Totals[i+k+1].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i/6 - 6] / *Totals[i].GetResultPtr<ULong64_t>()),2);
//                totalStatUnc[k] += pow(sqrt(*Totals[i+k+1].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i/6 - 6] / *Totals[i].GetResultPtr<ULong64_t>()),2);
//            }
//        }
//    }
//
//    std::vector<std::string> totalPrefixes = {R"--(Total $Z\gamma$)--", R"--(Total $Z$ jets)--", R"--(Total Bkg)--"};
//
//    std::cout << R"--(Total $Z\gamma$)--";
//    std::cout << " & " << std::setprecision(2) << std::fixed << total_Z_gamma[0]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(total_Z_gammaStatUnc[0]);
//    std::cout << " & " << std::setprecision(2) << std::fixed << total_Z_gamma[1]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(total_Z_gammaStatUnc[1]);
//    std::cout << " & " << std::setprecision(2) << std::fixed << total_Z_gamma[2]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(total_Z_gammaStatUnc[2]);
//    std::cout << " & " << std::setprecision(2) << std::fixed << total_Z_gamma[3]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(total_Z_gammaStatUnc[3]);
//    std::cout << " & " << std::setprecision(2) << std::fixed << total_Z_gamma[4]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(total_Z_gammaStatUnc[4]);
//    std::cout << R"--( \\ \hline)--" << '\n';
//
//    std::cout << R"--(Total $Z$ jets)--";
//    std::cout << " & " << std::setprecision(2) << std::fixed << total_Z_jets[0]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(total_Z_jetsStatUnc[0]);
//    std::cout << " & " << std::setprecision(2) << std::fixed << total_Z_jets[1]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(total_Z_jetsStatUnc[1]);
//    std::cout << " & " << std::setprecision(2) << std::fixed << total_Z_jets[2]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(total_Z_jetsStatUnc[2]);
//    std::cout << " & " << std::setprecision(2) << std::fixed << total_Z_jets[3]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(total_Z_jetsStatUnc[3]);
//    std::cout << " & " << std::setprecision(2) << std::fixed << total_Z_jets[4]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(total_Z_jetsStatUnc[4]);
//    std::cout << R"--( \\ \hline)--" << '\n';
//
//    std::cout << R"--(Total Bkg)--";
//    std::cout << " & " << std::setprecision(2) << std::fixed << total[0]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(totalStatUnc[0]);
//    std::cout << " & " << std::setprecision(2) << std::fixed << total[1]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(totalStatUnc[1]);
//    std::cout << " & " << std::setprecision(2) << std::fixed << total[2]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(totalStatUnc[2]);
//    std::cout << " & " << std::setprecision(2) << std::fixed << total[3]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(totalStatUnc[3]);
//    std::cout << " & " << std::setprecision(2) << std::fixed << total[4]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(totalStatUnc[4]);
//    std::cout << R"--( \\ \hline)--" << '\n';
//
//    std::cout << R"--(\end{tabular}})--" << '\n';
//
//    std::cout << "\n\n\n";
//}
//
//void Table19()
//{
//    std::vector<std::vector<std::string>> input_filenames =
//    {
//        //Z gamma background
//        {"/home/common/Za/NTuples/Background/user.kschmied.364860.eegammagamma_pty2_9_17.deriv.DAOD_STDM3.e7057_s3126_r10724_p4092_LGNTuple.root/user.kschmied.31617070._000001.LGNTuple.root"},
//        {"/home/common/Za/NTuples/Background/user.kschmied.364861.eegammagamma_pty_17_myy_0_80.deriv.DAOD_HIGG1D2.e7057_s3126_r10724_p4062_LGNTuple.root/user.kschmied.31617064._000001.LGNTuple.root"},
//        {"/home/common/Za/NTuples/Background/user.kschmied.364862.eegammagamma_pty_17_myy_80.deriv.DAOD_HIGG1D2.e7057_s3126_r10724_p4204_LGNTuple.root/user.kschmied.31660711._000001.LGNTuple.root", "/home/common/Za/NTuples/Background/user.kschmied.364862.eegammagamma_pty_17_myy_80.deriv.DAOD_HIGG1D2.e7057_s3126_r10724_p4204_LGNTuple.root/user.kschmied.31660711._000002.LGNTuple.root"},
//        //Data
//        {"/home/common/Za/NTuples/Ntuple_data_test.root"},
//        //Signal
//        {"/home/common/Za/NTuples/Ntuple_MC_Za_m5p0_v4.root"},
//        {"/home/common/Haa/ntuples/Za/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v3.root"},
//        //Jets
//        {
//            "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000001.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000002.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000003.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000004.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364114.v11.Zee_0_70_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835698._000005.LGNTuple.root",
//        },
//        {
//            "/home/common/Za/NTuples/Background/user.kschmied.364117.v11.Zee_70_140_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835715._000001.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364117.v11.Zee_70_140_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835715._000002.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364117.v11.Zee_70_140_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835715._000003.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364117.v11.Zee_70_140_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835715._000004.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364117.v11.Zee_70_140_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835715._000006.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364117.v11.Zee_70_140_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835715._000007.LGNTuple.root",
//        },
//        {
//            "/home/common/Za/NTuples/Background/user.kschmied.364120.v11.Zee_140_280_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835749._000001.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364120.v11.Zee_140_280_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835749._000002.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364120.v11.Zee_140_280_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835749._000003.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364120.v11.Zee_140_280_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835749._000004.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364120.v11.Zee_140_280_CVetoBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835749._000005.LGNTuple.root",
//        },
//        {
//            "/home/common/Za/NTuples/Background/user.kschmied.364115.v11.Zee_0_70_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835769._000001.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364115.v11.Zee_0_70_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835769._000002.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364115.v11.Zee_0_70_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835769._000003.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364115.v11.Zee_0_70_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835769._000004.LGNTuple.root",
//        },
//        {
//            "/home/common/Za/NTuples/Background/user.kschmied.364118.v11.Zee_70_140_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835773._000001.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364118.v11.Zee_70_140_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835773._000002.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364118.v11.Zee_70_140_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835773._000003.LGNTuple.root",
//        },
//        {
//            "/home/common/Za/NTuples/Background/user.kschmied.364121.v11.Zee_140_280_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835802._000001.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364121.v11.Zee_140_280_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835802._000002.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364121.v11.Zee_140_280_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835802._000003.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364121.v11.Zee_140_280_CFilterBVeto.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835802._000004.LGNTuple.root",
//        },
//        {
//            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000001.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000002.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000003.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000004.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000005.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000006.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000007.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364116.v11.Zee_0_70_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835825._000008.LGNTuple.root",
//        },
//        {
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000001.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000002.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000003.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000004.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000005.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000006.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000007.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000008.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000009.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000010.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000011.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000012.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000013.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364119.v11.Zee_70_140_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835838._000014.LGNTuple.root",
//        },
//        {
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000001.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000002.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000003.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000004.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000005.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000006.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000007.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000008.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000009.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000010.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000011.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000012.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000013.LGNTuple.root",
//            "/home/common/Za/NTuples/Background/user.kschmied.364122.v11.Zee_140_280_BFilter.DAOD_STDM3.e5299_s3126_r10724_p4252_LGNTuple.root/user.kschmied.31835843._000014.LGNTuple.root",
//        },
//    };
//
//    std::vector<std::string> prefixes = { R"--(pty2\_9\_17)--", R"--(pty\_17\_myy\_0\_80)--", R"--(pty\_17\_myy\_80)--", "data", R"--($\text{Sig } m_{A}$ = 5 GeV)--", R"--($\text{Sig } m_{A}$ = 1 GeV)--", R"--(Zee\_lightJet\_0-70)--", R"--(Zee\_lightJet\_70-140)--", R"--(Zee\_lightJet\_140-280)--", R"--(Zee\_cJet\_0-70)--", R"--(Zee\_cJet\_70-140)--", R"--(Zee\_cJet\_140-280)--", R"--(Zee\_bJet\_0-70)--", R"--(Zee\_bJet\_70-140)--", R"--(Zee\_bJet\_140-280)--"};
//
//    std::array<double,9> JetNumeratorSFs = {((139e15)*(1.9828e-9)*(0.821204)),((139e15)*(110.64e-12)*(0.69275)),((139e15)*(40.645e-12)*(0.615906)),((139e15)*(1.9817e-9)*(0.1136684)),((139e15)*(110.47e-12)*(0.1912956)),((139e15)*(40.674e-12)*(0.2326772)),((139e15)*(1.9819e-9)*(0.0656969)),((139e15)*(110.53e-12)*(0.1158741)),((139e15)*(40.68e-12)*(0.1535215))};
//    std::array<double,3> SFs = {((139e15)*(.871e-12)),((139e15)*(.199e-12)), ((139e15)*(.0345e-15))}; //numerators for Z-gamma bkg
//
//    std::vector<ROOT::RDF::RResultHandle> Totals;
//
//    for (auto& i: input_filenames)
//    {
//        SchottDataFrame df(MakeRDF(i,24));
////        std::cout << *df.Count() << '\n';
//
//        auto two_leptons = df.Filter(
//        [](RVec<Muon>& muons, RVec<Electron> electrons)
//        {
//            electrons.erase(std::remove_if(electrons.begin(),electrons.end(),
//            [](Electron& ep)
//            {
//                return (!((ep.electron_pt/1e3 > 20) && (abs(ep.electron_eta) < 2.37) &&
//                          (!((1.37 < abs(ep.electron_eta)) && (abs(ep.electron_eta) < 1.52)))
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
//                return (!((ep.electron_pt/1e3 > 20) && (abs(ep.electron_eta) < 2.37) &&
//                (!((1.37 < abs(ep.electron_eta)) && (abs(ep.electron_eta) < 1.52)))
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
//                return ((abs(x.photon_eta) >= 2.37) || (abs(x.photon_eta) > 1.37 && abs(x.photon_eta) < 1.52));
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
//                if (i==0 || abs(1-X) < abs(1-best_X))
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
////    0    1    2    3    4    5    6     Z-gamma
////    7    8    9    10   11   12   13    Z-gamma
////    14   15   16   17   18   19   20    Z-gamma
////    21   22   23   24   25   26   27    data
////    28   29   30   31   32   33   34    signal
////    35   36   37   38   39   40   41    signal
////    42   43   44   45   46   47   48    Z-jets
////    49   50   51   52   53   54   55    Z-jets
////    56   57   58   59   60   61   62    Z-jets
////    63   64   65   66   67   68   69    Z-jets
////    70   71   72   73   74   75   76    Z-jets
////    77   78   79   80   81   82   83    Z-jets
////    84   85   86   87   88   89   90    Z-jets
////    91   92   93   94   95   96   97    Z-jets
////    98   99   100  101  102  103  104   Z-jets
//
//    std::cout << R"--(\section*{Table 19})--" << '\n';
//    std::cout << R"--(\hspace{-3cm}\scalebox{0.55}{)--" << '\n';
//    std::cout << R"--(\begin{tabular}{|c|c|c|c|c|c|c|})--" << '\n';
//    std::cout << R"--(\hline)--" << '\n';
//    std::cout << R"--({} & pass preselection & failed resolved category & photon $p_T$ cut & pSR & SR & SR-ID \\ \hline)--" << '\n';
//
//    double total_Z_gamma[6] = {0}, total_Z_jets[6] = {0}, total[6] = {0};
//    double total_Z_gammaStatUnc[6] = {0}, total_Z_jetsStatUnc[6] = {0}, totalStatUnc[6] = {0};
//
//    for (int i = 0, j = 0; (i <= 98 && j <= 14); i += 7, j++)
//    {
//        std::cout << prefixes[j];
//        if (j >= 0 && j <= 2)
//        {
//            std::cout
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
//        else if (j >= 6)
//        {
//            std::cout
//            << " & " << std::setprecision(2) << std::fixed << *Totals[i+1].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j-6] / *Totals[i].GetResultPtr<ULong64_t>())
//            << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(*Totals[i+1].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[j-6] / *Totals[i].GetResultPtr<ULong64_t>())
//            << " & " << std::setprecision(2) << std::fixed << *Totals[i+2].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j-6] / *Totals[i].GetResultPtr<ULong64_t>())
//            << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(*Totals[i+2].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[j-6] / *Totals[i].GetResultPtr<ULong64_t>())
//            << " & " << std::setprecision(2) << std::fixed << *Totals[i+3].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j-6] / *Totals[i].GetResultPtr<ULong64_t>())
//            << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(*Totals[i+1].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[j-6] / *Totals[i].GetResultPtr<ULong64_t>())
//            << " & " << std::setprecision(2) << std::fixed << *Totals[i+4].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j-6] / *Totals[i].GetResultPtr<ULong64_t>())
//            << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(*Totals[i+1].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[j-6] / *Totals[i].GetResultPtr<ULong64_t>())
//            << " & " << std::setprecision(2) << std::fixed << *Totals[i+5].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j-6] / *Totals[i].GetResultPtr<ULong64_t>())
//            << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(*Totals[i+1].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[j-6] / *Totals[i].GetResultPtr<ULong64_t>())
//            << " & " << std::setprecision(2) << std::fixed << *Totals[i+6].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[j-6] / *Totals[i].GetResultPtr<ULong64_t>())
//            << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(*Totals[i+1].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[j-6] / *Totals[i].GetResultPtr<ULong64_t>())
//            << R"--( \\ \hline)--" << '\n';
//        }
//
//        else
//        {
//            std::cout
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
//        else if (i >= 42 && i <= 98) //Z-jets
//        {
//            for (int k = 0; k < 6; k++)
//            {
//                total_Z_jets[k] += *Totals[i+k+1].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i/7 - 6] / *Totals[i].GetResultPtr<ULong64_t>());
//                total[k] += *Totals[i+k+1].GetResultPtr<ULong64_t>() * (JetNumeratorSFs[i/7 - 6] / *Totals[i].GetResultPtr<ULong64_t>());
//
//                total_Z_jetsStatUnc[k] += pow(sqrt(*Totals[i+k+1].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i/7 - 6] / *Totals[i].GetResultPtr<ULong64_t>()),2);
//                totalStatUnc[k] += pow(sqrt(*Totals[i+k+1].GetResultPtr<ULong64_t>()) * (JetNumeratorSFs[i/7 - 6] / *Totals[i].GetResultPtr<ULong64_t>()),2);
//            }
//        }
//    }
//
//    std::vector<std::string> totalPrefixes = {R"--(Total $Z\gamma$)--", R"--(Total $Z$ jets)--", R"--(Total Bkg)--"};
//
//    std::cout << R"--(Total $Z\gamma$)--";
//    std::cout << " & " << std::setprecision(2) << std::fixed << total_Z_gamma[0]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(total_Z_gammaStatUnc[0]);
//    std::cout << " & " << std::setprecision(2) << std::fixed << total_Z_gamma[1]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(total_Z_gammaStatUnc[1]);
//    std::cout << " & " << std::setprecision(2) << std::fixed << total_Z_gamma[2]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(total_Z_gammaStatUnc[2]);
//    std::cout << " & " << std::setprecision(2) << std::fixed << total_Z_gamma[3]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(total_Z_gammaStatUnc[3]);
//    std::cout << " & " << std::setprecision(2) << std::fixed << total_Z_gamma[4]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(total_Z_gammaStatUnc[4]);
//    std::cout << " & " << std::setprecision(2) << std::fixed << total_Z_gamma[5]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(total_Z_gammaStatUnc[5]);
//    std::cout << R"--( \\ \hline)--" << '\n';
//
//    std::cout << R"--(Total $Z$ jets)--";
//    std::cout << " & " << std::setprecision(2) << std::fixed << total_Z_jets[0]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(total_Z_jetsStatUnc[0]);
//    std::cout << " & " << std::setprecision(2) << std::fixed << total_Z_jets[1]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(total_Z_jetsStatUnc[1]);
//    std::cout << " & " << std::setprecision(2) << std::fixed << total_Z_jets[2]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(total_Z_jetsStatUnc[2]);
//    std::cout << " & " << std::setprecision(2) << std::fixed << total_Z_jets[3]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(total_Z_jetsStatUnc[3]);
//    std::cout << " & " << std::setprecision(2) << std::fixed << total_Z_jets[4]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(total_Z_jetsStatUnc[4]);
//    std::cout << " & " << std::setprecision(2) << std::fixed << total_Z_jets[5]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(total_Z_jetsStatUnc[5]);
//
//    std::cout << R"--( \\ \hline)--" << '\n';
//
//    std::cout << R"--(Total Bkg)--";
//    std::cout << " & " << std::setprecision(2) << std::fixed << total[0]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(totalStatUnc[0]);
//    std::cout << " & " << std::setprecision(2) << std::fixed << total[1]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(totalStatUnc[1]);
//    std::cout << " & " << std::setprecision(2) << std::fixed << total[2]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(totalStatUnc[2]);
//    std::cout << " & " << std::setprecision(2) << std::fixed << total[3]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(totalStatUnc[3]);
//    std::cout << " & " << std::setprecision(2) << std::fixed << total[4]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(totalStatUnc[4]);
//    std::cout << " & " << std::setprecision(2) << std::fixed << total[5]
//    << R"--($\, \pm \,$)--" << std::setprecision(2) << std::fixed << sqrt(totalStatUnc[5]);
//    std::cout << R"--( \\ \hline)--" << '\n';
//
//    std::cout << R"--(\end{tabular}})--" << '\n';
//
//    std::cout << "\n\n\n";
//}

void DataBackgroundComparison()
{
    auto start_time = Clock::now();
//    fig27();
    fig28();
//    fig41();
//    fig48();
//    fig59();
//    Table9();
//    Table10();
//    Table16();
//    Table19();
    auto end_time = Clock::now();
    std::cout << "Time difference: "
       << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count()/1e9 << " seconds" << std::endl;
    
}

int main()
{
    DataBackgroundComparison();
}

