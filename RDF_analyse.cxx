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

//should use full-screen on this.
//just a helper function to print out variations for a
//column if they exist. Not used in the event-loop, more like
//a sanity check.
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

//Select reco-photons that pass a photon id criterion, pt >= 2.5 GeV,
//|eta| <= 2.37, |eta| <= 1.37 or |eta| >= 1.52 (interlayer of end-cap
//and barrel calorimeters in the ATLAS detector)
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
    
    if (std::abs(photon.photon_eta) > 2.37)
    {
        return false;
    }
    
    if ((1.37 < std::abs(photon.photon_eta)) && (std::abs(photon.photon_eta) < 1.52))
    {
        return false;
    }
    return true;
}

//Select truth-photons that pass a photon id criterion, pt >= 2.5 GeV,
//|eta| <= 2.37, |eta| <= 1.37 or |eta| >= 1.52 (interlayer of end-cap
//and barrel calorimeters in the ATLAS detector)
bool truth_photon_selection(TruthParticle& truth_photon)
{
    if (truth_photon.mc_pt < 2500)
    {
        return false;
    }
    
    if (std::abs(truth_photon.mc_eta) > 2.37)
    {
        return false;
    }
    
    if ((1.37 < std::abs(truth_photon.mc_eta)) && (std::abs(truth_photon.mc_eta) < 1.52))
    {
        return false;
    }
    return true;
}

//Select tracks that have pt >= 0.1 GeV and |eta| <= 2.5
bool track_selection(Track& track)
{
    if (track.track_pt < 100)
    {
        return false;
    }
    
    if (std::abs(track.track_eta) > 2.5)
    {
        return false;
    }
    
    return true;
}
//https://root-forum.cern.ch/t/rdataframe-count-and-report-re-looping-over-whole-dataframe/46592
void RDF_analyse()
{
    auto start_time = Clock::now();
//    Prints out the info of the event loop(s) running (in this case 1 event
//    loop because we do it correctly :) ) in a verbose manner
    auto verbosity = ROOT::Experimental::RLogScopedVerbosity(ROOT::Detail::RDF::RDFLogChannel(), ROOT::Experimental::ELogLevel::kInfo);
    
//    List of systematics for this analysis
    Event::systematics = {"EG_RESOLUTION_ALL__1down"};
//    Input file names
    std::vector<std::string> input_filenames = {"/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/Ntuple_MC_Za_mA5p0_v4"};
//    Construct the RDataFrame
    SchottDataFrame df(MakeRDF(input_filenames));
    
//    df.Describe().Print();
//    exit(1);
//    
//    printNominalAndVariedObjects<Electron>(df, "electrons", 10);
//    printNominalAndVariedObjects<Photon>(df, "photons", 10);
    bool mc = true;
    
//    Book a filtration of events for those that don't have a trigger
//    and if mc is false. So this filter does nothing!
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
    
//    Define a new column for reco-photons that pass photon_selection (see above)
    auto selected_photons_and_tracks = firstTriggerCut.Define("selected_photons", [&](RVec<Photon> photons)
    {
        photons.erase(std::remove_if(photons.begin(),photons.end(),
        [](Photon& x)
        {
            return !photon_selection(x);
            
        }), photons.end());
        
        return photons;
    }, {"photons"})
//    Define a new column for tracks that pass track_selection (see above)
    .Define("selected_tracks", [&](RVec<Track> tracks)
    {
       tracks.erase(std::remove_if(tracks.begin(),tracks.end(),
       [](Track& x)
       {
           return !track_selection(x);
           
       }), tracks.end());
       
       return tracks;
    }, {"tracks"});

//    Define a new column for truth-photons that pass truth_photon_selection
//    (see above)
    auto truth_photons = firstTriggerCut.Define("truth_photons", [&](RVec<TruthParticle> truth_particles)
    {
        truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
        [](TruthParticle& x)
        {
            return (x.mc_pdg_id!=22 || !truth_photon_selection(x));
            
        }), truth_particles.end());
        
        return truth_particles;
    }, {"truth_particles"});
    
//    Book a filtration of events: keep events that have exactly
//    2 selected truth photons and whose diphoton invariant mass > 5 GeV
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
//    For the events that passed the above Filter, define a new column
//    that holds the pt of the selected truth photons
    .Define("truth_candidates_pt", [&](RVec<TruthParticle>& truth_photons)
    {
        return RVec<float>({static_cast<float>(truth_photons[0].mc_pt/1e3), static_cast<float>(truth_photons[1].mc_pt/1e3)});
        
    }, {"truth_photons"})
//    Do the same for the diphoton invariant mass
    .Define("truth_candidates_mass", [&](RVec<TruthParticle>& truth_photons)
    {
        PtEtaPhiEVector four_momentum =
        truth_photons[0].Vector() + truth_photons[1].Vector();
        
        return four_momentum.M()/1e3;
    }, {"truth_photons"});
    
//    Book a new filtration of events: remove events that don't have
//    exactly two selected photons
    auto diphotons_and_tracks = selected_photons_and_tracks.Filter(
    [&](RVec<Photon>& selected_photons)
    {
        return (selected_photons.size() == 2);
       
    }, {"selected_photons"}, "diphotons_and_tracks")
//    For the events that passed the above Filter, define a new column to hold
//    the pt of the two photons per event
    .Define("diphotons_pt", [&](RVec<Photon>& diphotons)
    {
        return RVec<float>({static_cast<float>(diphotons[0].photon_pt/1e3), static_cast<float>(diphotons[1].photon_pt/1e3)});
        
    }, {"selected_photons"})
//    Do the same for the diphoton invariant mass
    .Define("diphotons_mass", [&](RVec<Photon>& diphotons)
    {
        PtEtaPhiEVector four_momentum =
        diphotons[0].Vector() + diphotons[1].Vector();
        
        return four_momentum.M()/1e3;
    }, {"selected_photons"})
//    Do the same for the number of tracks in each of the remaining events
    .Define("num_tracks", [&](RVec<Track>& tracks)
    {
        return tracks.size();
        
    }, {"selected_tracks"})
//    Do the same for the pt of each track in each event
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
//    Do the same for the eta of each track in each event
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

//    For the diphotons_and_tracks node that we defined above, book a
//    filter for events of this node that will discard events with
//    selected tracks and discard events with diphoton invariant mass
//    > 5 GeV
    auto no_tracks_inv_mass = diphotons_and_tracks.Filter(
    [&](RVec<Track>& selected_tracks, RVec<Photon>& diphotons)
    {
        PtEtaPhiEVector four_momentum =
        diphotons[0].Vector() + diphotons[1].Vector();
        
        return (selected_tracks.empty() && four_momentum.M() > 5000);
     
    }, {"selected_tracks", "selected_photons"}, "no_tracks_inv_mass")
//    For events that pass the above Filter, define a new column for the
//    pt of the two photons
    .Define("diphotons_pt_no_tracks_inv_mass", [&](RVec<Photon>& diphotons)
    {
        return RVec<float>({static_cast<float>(diphotons[0].photon_pt/1e3), static_cast<float>(diphotons[1].photon_pt/1e3)});
        
    }, {"selected_photons"})
//    Do the same for the diphoton invariant mass
    .Define("diphotons_mass_no_tracks_inv_mass", [&](RVec<Photon>& diphotons)
    {
        PtEtaPhiEVector four_momentum =
        diphotons[0].Vector() + diphotons[1].Vector();
        
        return four_momentum.M()/1e3;
    }, {"selected_photons"})
//    Do the same for the number of tracks in each of the remaining events
    .Define("num_tracks_no_tracks_inv_mass", [&](RVec<Track>& tracks)
    {
        return tracks.size();
        
    }, {"selected_tracks"})
//    Do the same for the pt of each track in each event
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
//    Do the same for the pt of each track in each event
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
    
//    For the diphotons_and_tracks node defined above, create a new node that will
//    hold a lazy removal of events that have selected_tracks not empty
    auto no_tracks = diphotons_and_tracks.Filter(
    [&](RVec<Track>& selected_tracks)
    {
        return (selected_tracks.empty());
        
    }, {"selected_tracks"}, "no_tracks")
//    For the events that passed the above Filter, define a new column for the
//    pt of the two selected photons in each of the remaining events
    .Define("diphotons_pt_no_tracks", [&](RVec<Photon>& diphotons)
    {
        return RVec<float>({static_cast<float>(diphotons[0].photon_pt/1e3), static_cast<float>(diphotons[1].photon_pt/1e3)});
        
    }, {"selected_photons"})
//    Do the same for the diphoton invariant mass
    .Define("diphotons_mass_no_tracks", [&](RVec<Photon>& diphotons)
    {
        PtEtaPhiEVector four_momentum =
        diphotons[0].Vector() + diphotons[1].Vector();
        
        return four_momentum.M()/1e3;
    }, {"selected_photons"})
//    Do the same for the number of tracks in each of the remaining events
    .Define("num_tracks_no_tracks", [&](RVec<Track>& tracks)
    {
        return tracks.size();
        
    }, {"selected_tracks"})
//    Do the same for the pt of each track in each event
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
//    Do the same for the eta of each track in each event
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
    
//    Histograms that we want to eventually draw. Note, the event-loop
//    has NOT been triggered yet!
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
    
//    vector of RResultMaps that will hold the nominal and varied results for
//    each histogram (technically, any 'action', but we're just doing histograms here)
    std::vector<RResultMap<TH1D>> resultmaps;
//    Vector to store the name of each histogram
    std::vector<std::string> histNames;
    resultmaps.reserve(53); //good practice to reserve capacity for vectors
    histNames.reserve(53);  //good practice to reserve capacity for vectors

//    Looping over all histos, storing maps of nominal/varied
//    results and histogram names
    for (auto& h: histos)
    {
        resultmaps.push_back(VariationsFor(h));
        histNames.push_back(h->GetName());
    }

    TCanvas *c1;
    std::string str;
    for (auto i = 0; i < resultmaps.size(); i++) //For each histogram
    {
        for (auto& var: resultmaps[i].GetKeys()) //For each variation of said histogram
        {
            c1 = new TCanvas("","",800, 700);
            resultmaps[i][var].Draw("same"); //Draw the histogram
            str = var;
            str.erase(std::remove(str.begin(), str.end(), ':'), str.end());
            c1->SaveAs((str+histNames[i]+".png").c_str()); //save it
            if (histNames[i] == "02MassCutTrackingNumTracks")
            {
//                std::cout << *no_tracks_inv_mass.Count() << '\n';
                std::cout << "# events for " << var << " = "
                << resultmaps[i][var].GetEntries() << '\n'; //print out some interesting info while we're at it :)
            }
        }
    }
    
//    df.Describe().Print();
    auto end_time = Clock::now(); //done
    std::cout << "\nTime difference: "
       << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count()/1e9 << " seconds" << std::endl;
    
    
}


int main()
{
    RDF_analyse();
}



