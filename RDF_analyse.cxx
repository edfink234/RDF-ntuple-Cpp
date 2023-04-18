#include <iostream>
#include <vector>
#include <chrono>
#include <algorithm>
#include <string>
#include <cstdlib>
#include <iomanip>
#include <memory>

#include <ROOT/RLogger.hxx>
#include "Math/VectorUtil.h"
//#include "ROOT/RDF/RResultMap.hxx"
#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RDataFrame.hxx"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH1D.h"

#include "RDFObjects.h"
#include "MakeRDF.h"
#include "RDFevent.h"

using namespace ROOT::VecOps; // RVec
using namespace ROOT::Math::VectorUtil; // DeltaR
using namespace ROOT::Math; // PtEtaPhiEVector

using ROOT::RDF::Experimental::RResultMap;
using ROOT::RDF::Experimental::VariationsFor;
//using ct = typename ROOT::RDF::RDisplay;
//using ct = typename ROOT::RDF::Experimental::RResultMap<RDisplay>;

using Clock = std::chrono::high_resolution_clock;

//should use full-screen on this
template<typename T>
void printNominalAndVariedObjects(SchottDataFrame& df, const char* column, int rows)
{
    auto v = df.Take<RVec<T>>(column);
    auto vs = VariationsFor(v);
    int diffCount = 0;
    const std::vector<std::string>& keys = vs.GetKeys();
    std::string temp;
    if (vs.GetKeys().empty())
    {
        std::cout << "Error, no keys!\n";
        return;
    }
    for (int i = 0; i < rows; i++)
    {
        std::cout << "\n\n";
        for (int j = 0; j < vs[keys[0]][i].size(); ++j)
        {
            temp = cling::printValue(&vs[keys[0]][i][j]);
            for (const auto& k: keys)
            {
                std::cout <<
                cling::printValue(&vs[k][i][j]) << "     ";
                if (temp != cling::printValue(&vs[k][i][j]))
                {
                    ++diffCount;
                }
            }
            std::cout << '\n';
        }
        std::cout << '\n';
    }
    
    std::cout << "\nNumber of different objects = " << diffCount << '\n';
}

bool photon_selection(Photon& photon)
{
    if (!photon.photon_id)
    {
        return false;
    }
    
    if (photon.photon_pt < 2500)
    {
        return false;
    }
    
    if (abs(photon.photon_eta) > 2.37)
    {
        return false;
    }
    
    if ((1.37 < abs(photon.photon_eta)) && (abs(photon.photon_eta) < 1.52))
    {
        return false;
    }
    return true;
}

bool truth_photon_selection(TruthParticle& truth_photon)
{
    if (truth_photon.mc_pt < 2500)
    {
        return false;
    }
    
    if (abs(truth_photon.mc_eta) > 2.37)
    {
        return false;
    }
    
    if ((1.37 < abs(truth_photon.mc_eta)) && (abs(truth_photon.mc_eta) < 1.52))
    {
        return false;
    }
    return true;
}

bool track_selection(Track& track)
{
    if (track.track_pt < 100)
    {
        return false;
    }
    
    if (abs(track.track_eta) > 2.5)
    {
        return false;
    }
    
    return true;
}
//https://root-forum.cern.ch/t/rdataframe-count-and-report-re-looping-over-whole-dataframe/46592
void RDF_analyse()
{
    auto start_time = Clock::now();
    auto verbosity = ROOT::Experimental::RLogScopedVerbosity(ROOT::Detail::RDF::RDFLogChannel(), ROOT::Experimental::ELogLevel::kInfo);
    
    Event::systematics = {"EG_RESOLUTION_ALL__1down"};
    std::vector<std::vector<std::string>> input_filenames = {{"/home/common/Za/NTuples/Ntuple_MC_Za_m5p0_v4.root"}};

    SchottDataFrame df(MakeRDF(input_filenames[0]));
    
//    df.Describe().Print();
    exit(1);
    
//    printNominalAndVariedObjects<Electron>(df, "electrons", 10);
//    printNominalAndVariedObjects<Photon>(df, "photons", 10);
    bool mc = true;
    
    auto firstTriggerCut = df.Filter(
    [&](const RVec<std::string>& trigger_passed_triggers)
    {
//        return !(!mc && std::find(trigger_passed_triggers.begin(), trigger_passed_triggers.end(), "HLT_hi_upc_FgapAC3_hi_gg_upc_L1TAU1_TE4_VTE200") == trigger_passed_triggers.end());
        
//         The commented out return statement above is
//         equivalent to the one below
        
        return (mc || std::find(trigger_passed_triggers.begin(), trigger_passed_triggers.end(), "HLT_hi_upc_FgapAC3_hi_gg_upc_L1TAU1_TE4_VTE200") != trigger_passed_triggers.end());
        
//        Proof: This Python code doesn't give an AssertionError
//        ======================================================
//        >>> choices = (True, False)
//        >>> for x in choices:
//        ...     for y in choices:
//        ...             assert( not (not x and not y) == (x or y))
        
    }, {"trigger_passed_triggers"}, "firstTriggerCut");
                                     
//    auto passed = firstTriggerCut.Count();
//    std::cout << *passed << '\n';
    
    auto selected_photons_and_tracks = firstTriggerCut.Define("selected_photons", [&](RVec<Photon> photons)
    {
        photons.erase(std::remove_if(photons.begin(),photons.end(),
        [](Photon& x)
        {
            return !photon_selection(x);
            
        }), photons.end());
        
        return photons;
    }, {"photons"})
    .Define("selected_tracks", [&](RVec<Track> tracks)
    {
       tracks.erase(std::remove_if(tracks.begin(),tracks.end(),
       [](Track& x)
       {
           return !track_selection(x);
           
       }), tracks.end());
       
       return tracks;
    }, {"tracks"});

    auto truth_photons = firstTriggerCut.Define("truth_photons", [&](RVec<TruthParticle> truth_particles)
    {
        truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
        [](TruthParticle& x)
        {
            return (x.mc_pdg_id!=22 || !truth_photon_selection(x));
            
        }), truth_particles.end());
        
        return truth_particles;
    }, {"truth_particles"});
    
    auto truth_candidates = truth_photons.Filter(
    [&](RVec<TruthParticle>& truth_photons)
    {
        if (truth_photons.size() == 2)
        {
            PtEtaPhiEVector four_momentum =
            truth_photons[0].Vector() + truth_photons[1].Vector();
            
            if (four_momentum.M() > 5000)
            {
                return true;
            }
        }
        return false;
    }, {"truth_photons"}, "truth_candidates")
    .Define("truth_candidates_pt", [&](RVec<TruthParticle>& truth_photons)
    {
        return RVec<float>({static_cast<float>(truth_photons[0].mc_pt/1e3), static_cast<float>(truth_photons[1].mc_pt/1e3)});
        
    }, {"truth_photons"})
    .Define("truth_candidates_mass", [&](RVec<TruthParticle>& truth_photons)
    {
        PtEtaPhiEVector four_momentum =
        truth_photons[0].Vector() + truth_photons[1].Vector();
        
        return four_momentum.M()/1e3;
    }, {"truth_photons"});
    
    auto diphotons_and_tracks = selected_photons_and_tracks.Filter(
    [&](RVec<Photon>& selected_photons)
    {
        return (selected_photons.size() == 2);
       
    }, {"selected_photons"}, "diphotons_and_tracks")
    .Define("diphotons_pt", [&](RVec<Photon>& diphotons)
    {
        return RVec<float>({static_cast<float>(diphotons[0].photon_pt/1e3), static_cast<float>(diphotons[1].photon_pt/1e3)});
        
    }, {"selected_photons"})
    .Define("diphotons_mass", [&](RVec<Photon>& diphotons)
    {
        PtEtaPhiEVector four_momentum =
        diphotons[0].Vector() + diphotons[1].Vector();
        
        return four_momentum.M()/1e3;
    }, {"selected_photons"})
    .Define("num_tracks", [&](RVec<Track>& tracks)
    {
        return tracks.size();
        
    }, {"selected_tracks"})
    .Define("selected_tracks_pt", [&](RVec<Track>& tracks)
    {
        RVec<float> selected_tracks_pt;
        selected_tracks_pt.reserve(tracks.size());
        for (auto& track: tracks)
        {
            selected_tracks_pt.push_back(track.track_pt/1e3);
        }
        return selected_tracks_pt;
        
    }, {"selected_tracks"})
    .Define("selected_tracks_eta", [&](RVec<Track>& tracks)
    {
        RVec<float> selected_tracks_eta;
        selected_tracks_eta.reserve(tracks.size());
        for (auto& track: tracks)
        {
            selected_tracks_eta.push_back(track.track_eta);
        }
        return selected_tracks_eta;
        
    }, {"selected_tracks"});
    
    auto no_tracks_inv_mass = diphotons_and_tracks.Filter(
    [&](RVec<Track>& selected_tracks, RVec<Photon>& diphotons)
    {
        PtEtaPhiEVector four_momentum =
        diphotons[0].Vector() + diphotons[1].Vector();
        
        return (selected_tracks.empty() && four_momentum.M() > 5000);
     
    }, {"selected_tracks", "selected_photons"}, "no_tracks_inv_mass")
    .Define("diphotons_pt_no_tracks_inv_mass", [&](RVec<Photon>& diphotons)
    {
        return RVec<float>({static_cast<float>(diphotons[0].photon_pt/1e3), static_cast<float>(diphotons[1].photon_pt/1e3)});
        
    }, {"selected_photons"})
    .Define("diphotons_mass_no_tracks_inv_mass", [&](RVec<Photon>& diphotons)
    {
        PtEtaPhiEVector four_momentum =
        diphotons[0].Vector() + diphotons[1].Vector();
        
        return four_momentum.M()/1e3;
    }, {"selected_photons"})
    .Define("num_tracks_no_tracks_inv_mass", [&](RVec<Track>& tracks)
    {
        return tracks.size();
        
    }, {"selected_tracks"})
    .Define("selected_tracks_pt_no_tracks_inv_mass", [&](RVec<Track>& tracks)
    {
        RVec<float> selected_tracks_pt;
        selected_tracks_pt.reserve(tracks.size());
        for (auto& track: tracks)
        {
            selected_tracks_pt.push_back(track.track_pt/1e3);
        }
        return selected_tracks_pt;
        
    }, {"selected_tracks"})
    .Define("selected_tracks_eta_no_tracks_inv_mass", [&](RVec<Track>& tracks)
    {
        RVec<float> selected_tracks_eta;
        selected_tracks_eta.reserve(tracks.size());
        for (auto& track: tracks)
        {
            selected_tracks_eta.push_back(track.track_eta);
        }
        return selected_tracks_eta;
        
    }, {"selected_tracks"});
    
    
    auto no_tracks = diphotons_and_tracks.Filter(
    [&](RVec<Track>& selected_tracks)
    {
        return (selected_tracks.empty());
        
    }, {"selected_tracks"}, "no_tracks")
    .Define("diphotons_pt_no_tracks", [&](RVec<Photon>& diphotons)
    {
        return RVec<float>({static_cast<float>(diphotons[0].photon_pt/1e3), static_cast<float>(diphotons[1].photon_pt/1e3)});
        
    }, {"selected_photons"})
    .Define("diphotons_mass_no_tracks", [&](RVec<Photon>& diphotons)
    {
        PtEtaPhiEVector four_momentum =
        diphotons[0].Vector() + diphotons[1].Vector();
        
        return four_momentum.M()/1e3;
    }, {"selected_photons"})
    .Define("num_tracks_no_tracks", [&](RVec<Track>& tracks)
    {
        return tracks.size();
        
    }, {"selected_tracks"})
    .Define("selected_tracks_pt_no_tracks", [&](RVec<Track>& tracks)
    {
        RVec<float> selected_tracks_pt;
        selected_tracks_pt.reserve(tracks.size());
        for (auto& track: tracks)
        {
            selected_tracks_pt.push_back(track.track_pt/1e3);
        }
        return selected_tracks_pt;
        
    }, {"selected_tracks"})
    .Define("selected_tracks_eta_no_tracks", [&](RVec<Track>& tracks)
    {
        RVec<float> selected_tracks_eta;
        selected_tracks_eta.reserve(tracks.size());
        for (auto& track: tracks)
        {
            selected_tracks_eta.push_back(track.track_eta);
        }
        return selected_tracks_eta;
        
    }, {"selected_tracks"});

//    auto passed = truth_candidates.Count();
//    std::cout << *passed << '\n';
    
//    std::cout << '\n' << selected_photons.Display<RVec<Photon>,RVec<Photon>> ({"photons", "selected_photons"},100)->AsString() << '\n';
    
    std::vector<ROOT::RDF::RResultPtr<TH1D>> histos =
    {
        truth_candidates.Histo1D<RVec<float>>({"TruthRecoPhotonPt", "TruthRecoPhotonPt", 20u, 0, 25}, "truth_candidates_pt"),
        truth_candidates.Histo1D<RVec<float>>({"TruthRecoPhotonPtTight", "TruthRecoPhotonPtTight", 100u, 0, 10}, "truth_candidates_pt"),
        truth_candidates.Histo1D<double>({"TruthCandidateMass", "TruthCandidateMass", 250u, 0, 50}, "truth_candidates_mass"),
        truth_candidates.Histo1D<double>({"TruthCandidateMassLarge", "TruthCandidateMassLarge", 200u, 0, 200}, "truth_candidates_mass"),
        truth_candidates.Histo1D<double>({"TruthCandidateMassLargeFine", "TruthCandidateMassLargeFine", 600u, 0, 200}, "truth_candidates_mass"),
        
        diphotons_and_tracks.Histo1D<RVec<float>>({"00NoCutsRecoPhotonPt", "00NoCutsRecoPhotonPt", 20u, 0, 25}, "diphotons_pt"),
        diphotons_and_tracks.Histo1D<RVec<float>>({"00NoCutsRecoPhotonPtTight", "00NoCutsRecoPhotonPtTight", 100u, 0, 10}, "diphotons_pt"),
        diphotons_and_tracks.Histo1D<double>({"00NoCutsCandidateMass", "00NoCutsCandidateMass", 250u, 0, 50}, "diphotons_mass"),
        diphotons_and_tracks.Histo1D<double>({"00NoCutsCandidateMassLarge", "00NoCutsCandidateMassLarge", 200u, 0, 200}, "diphotons_mass"),
        diphotons_and_tracks.Histo1D<double>({"00NoCutsCandidateMassLargeFine", "00NoCutsCandidateMassLargeFine", 600u, 0, 200}, "diphotons_mass"),
        diphotons_and_tracks.Histo1D<unsigned long>({"00NoCutsTrackingNumTracks", "00NoCutsTrackingNumTracks", 20u, 0, 20}, "num_tracks"),
        diphotons_and_tracks.Histo1D<RVec<float>>({"00NoCutsTrackingTrackPt", "00NoCutsTrackingTrackPt", 140u, 0, 7}, "selected_tracks_pt"),
        diphotons_and_tracks.Histo1D<RVec<float>>({"00NoCutsTrackingTrackEta", "00NoCutsTrackingTrackEta", 50u, -2.5, 2.5}, "selected_tracks_eta"),
        
        no_tracks_inv_mass.Histo1D<RVec<float>>({"02MassCutRecoPhotonPt", "02MassCutRecoPhotonPt", 20u, 0, 25}, "diphotons_pt_no_tracks_inv_mass"),
        no_tracks_inv_mass.Histo1D<RVec<float>>({"02MassCutRecoPhotonPtTight", "02MassCutRecoPhotonPtTight", 100u, 0, 10}, "diphotons_pt_no_tracks_inv_mass"),
        no_tracks_inv_mass.Histo1D<double>({"02MassCutCandidateMass", "02MassCutCandidateMass", 250u, 0, 50}, "diphotons_mass_no_tracks_inv_mass"),
        no_tracks_inv_mass.Histo1D<double>({"02MassCutCandidateMassLarge", "02MassCutCandidateMassLarge", 200u, 0, 200}, "diphotons_mass_no_tracks_inv_mass"),
        no_tracks_inv_mass.Histo1D<double>({"02MassCutCandidateMassLargeFine", "02MassCutCandidateMassLargeFine", 600u, 0, 200}, "diphotons_mass_no_tracks_inv_mass"),
        no_tracks_inv_mass.Histo1D<unsigned long>({"02MassCutTrackingNumTracks", "02MassCutTrackingNumTracks", 20u, 0, 20}, "num_tracks_no_tracks_inv_mass"),
        no_tracks_inv_mass.Histo1D<RVec<float>>({"02MassCutTrackingTrackPt", "02MassCutTrackingTrackPt", 140u, 0, 7}, "selected_tracks_pt_no_tracks_inv_mass"),
        no_tracks_inv_mass.Histo1D<RVec<float>>({"02MassCutTrackingTrackEta", "02MassCutTrackingTrackEta", 50u, -2.5, 2.5}, "selected_tracks_eta_no_tracks_inv_mass"),
        
        no_tracks.Histo1D<RVec<float>>({"01NoTracksRecoPhotonPt", "01NoTracksRecoPhotonPt", 20u, 0, 25}, "diphotons_pt_no_tracks"),
        no_tracks.Histo1D<RVec<float>>({"01NoTracksRecoPhotonPtTight", "01NoTracksRecoPhotonPtTight", 100u, 0, 10}, "diphotons_pt_no_tracks"),
        no_tracks.Histo1D<double>({"01NoTracksCandidateMass", "01NoTracksCandidateMass", 250u, 0, 50}, "diphotons_mass_no_tracks"),
        no_tracks.Histo1D<double>({"01NoTracksCandidateMassLarge", "01NoTracksCandidateMassLarge", 200u, 0, 200}, "diphotons_mass_no_tracks"),
        no_tracks.Histo1D<double>({"01NoTracksCandidateMassLargeFine", "01NoTracksCandidateMassLargeFine", 600u, 0, 200}, "diphotons_mass_no_tracks"),
        no_tracks.Histo1D<unsigned long>({"01NoTracksTrackingNumTracks", "01NoTracksTrackingNumTracks", 20u, 0, 20}, "num_tracks_no_tracks"),
        no_tracks.Histo1D<RVec<float>>({"01NoTracksTrackingTrackPt", "01NoTracksTrackingTrackPt", 140u, 0, 7}, "selected_tracks_pt_no_tracks"),
        no_tracks.Histo1D<RVec<float>>({"01NoTracksTrackingTrackEta", "01NoTracksTrackingTrackEta", 50u, -2.5, 2.5}, "selected_tracks_eta_no_tracks"),
        
    };
    
    std::vector<RResultMap<TH1D>> resultmaps;
    std::vector<std::string> histNames;
    resultmaps.reserve(53);
    histNames.reserve(53);

    for (auto& h: histos)
    {
        resultmaps.push_back(VariationsFor(h));
        histNames.push_back(h->GetName());
    }

    TCanvas *c1;
    std::string str;
    for (auto i = 0; i < resultmaps.size(); i++)
    {
        for (auto& var: resultmaps[i].GetKeys())
        {
            c1 = new TCanvas("","",800, 700);
            resultmaps[i][var].Draw("same");
            str = var;
            str.erase(std::remove(str.begin(), str.end(), ':'), str.end());
            c1->SaveAs((str+histNames[i]+".png").c_str());
            if (histNames[i] == "02MassCutTrackingNumTracks")
            {
//                std::cout << *no_tracks_inv_mass.Count() << '\n';
                std::cout << "# events for " << var << " = "
                << resultmaps[i][var].GetEntries() << '\n';
            }
        }
    }
    
//    df.Describe().Print();
    auto end_time = Clock::now();
    std::cout << "\nTime difference:"
       << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count() << " nanoseconds" << std::endl;
    
}


int main()
{
    RDF_analyse();
}



