#include <vector>
#include <string>
#include <algorithm>

#include "TChain.h"
#include "TH1F.h"
#include "TBranch.h"
#include "TInterpreter.h"

//#include "RDFObjects.h"
//#include "MakeRDF.h"
//#include "RDFevent.h"

#include "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/RDFObjects.h"
#include "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/MakeRDF.h"
#include "/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/RDFevent.h"

TChain* RDFTree::__chain;
TChain* RDFTree::__event_info_chain;

using namespace ROOT::VecOps;

/*
 Creates an RDataFrame out of the physics and full_event_info TTrees,
 then adds columns for all of the objects defined in RDFObjects.h.
 The resulting dataframe is returned.
*/
SchottDataFrame MakeRDF(const std::vector<std::string>& files, short numThreads)
{
    if (numThreads > 0) //Enabling Multi-threading 
    {
        ROOT::EnableImplicitMT(numThreads);
    }
    static bool loaded = false;
    
    if (!loaded)
    {
        //load objects
//        gInterpreter->LoadMacro("RDFObjects.h");
        gInterpreter->LoadMacro("/Users/edwardfinkelstein/ATLAS_axion/ntupleC++_v2/RDFObjects.h");
        //then load printValue overloads so objects can be printed
        gInterpreter->Declare("std::string cling::printValue(TruthParticle *);");
        gInterpreter->Declare("std::string cling::printValue(Electron *);");
        gInterpreter->Declare("std::string cling::printValue(Muon *);");
        gInterpreter->Declare("std::string cling::printValue(Photon *);");
        gInterpreter->Declare("std::string cling::printValue(Cluster *);");
        gInterpreter->Declare("std::string cling::printValue(Track *);");
        loaded = true;
    }
    RDFTree::__chain = new TChain("physics");
    RDFTree::__event_info_chain = new TChain("full_event_info");
    
    RDFTree::__chain->Reset();
    RDFTree::__event_info_chain->Reset();
    for (const auto& f: files)
    {
        RDFTree::__chain->Add(f.c_str());
        RDFTree::__event_info_chain->Add(f.c_str());
    }
    RDFTree::__chain->AddFriend(RDFTree::__event_info_chain);

    ROOT::RDataFrame df(*RDFTree::__chain);
    
    auto NewDf = df.Define("truth_particles",[&](RVec<int>& mc_pdg_id, RVec<int>& mc_barcode, RVec<int>& mc_parent_barcode, RVec<int>& mc_status, RVec<float>& mc_pt, RVec<float>& mc_charge, RVec<float>& mc_eta, RVec<float>& mc_phi, RVec<float>& mc_e, RVec<float>& mc_mass)
    {
        RVec<TruthParticle> x;
        x.reserve(mc_pt.size());
        TruthParticle temp;
        for (size_t i = 0; i < mc_pt.size(); i++)
        {
            temp.mc_pdg_id =  mc_pdg_id[i];
            temp.mc_barcode =  mc_barcode[i];
            temp.mc_parent_barcode =  mc_parent_barcode[i];
            temp.mc_status =  mc_status[i];
            temp.mc_pt =  mc_pt[i];
            temp.mc_charge =  mc_charge[i];
            temp.mc_eta =  mc_eta[i];
            temp.mc_phi =  mc_phi[i];
            temp.mc_e =  mc_e[i];
            temp.mc_mass =  mc_mass[i];
            x.push_back(temp);
        }
        return x;
    }, {"mc_pdg_id", "mc_barcode", "mc_parent_barcode", "mc_status", "mc_pt", "mc_charge", "mc_eta", "mc_phi", "mc_e", "mc_mass"})
    .Define("electrons",[&](RVec<float>& electron_charge, RVec<float>& electron_pt, RVec<float>& electron_e, RVec<float>& electron_eta, RVec<float>& electron_phi, /*RVec<float>& electron_id,*/ RVec<float>& electron_isolation, RVec<float>& electron_d0, RVec<float>& electron_z0 , RVec<int>& electron_id_medium)
    {
        RVec<Electron> x;
        x.reserve(electron_pt.size());
        Electron temp;
        for (size_t i = 0; i < electron_pt.size(); i++)
        {
            temp.electron_charge =  electron_charge[i];
            temp.electron_pt =  electron_pt[i];
            temp.electron_e =  electron_e[i];
            temp.electron_eta =  electron_eta[i];
            temp.electron_phi =  electron_phi[i];
//            temp.electron_id =  electron_id[i];
            temp.electron_isolation =  electron_isolation[i];
            temp.electron_d0 =  electron_d0[i];
            temp.electron_z0 =  electron_z0[i];
            temp.electron_id_medium =  electron_id_medium[i];
            x.push_back(temp);
        }
        return x;
    }, {"electron_charge", "electron_pt", "electron_e", "electron_eta", "electron_phi", /*"electron_id",*/ "electron_isolation", "electron_d0", "electron_z0", "electron_id_medium",})
    .Define("muons",[&](RVec<int>& muon_charge, RVec<float>& muon_pt, RVec<float>& muon_e, RVec<float>& muon_eta, RVec<float>& muon_phi)
    {
        RVec<Muon> x;
        x.reserve(muon_pt.size());
        Muon temp;
        for (size_t i = 0; i < muon_pt.size(); i++)
        {
            temp.muon_charge =  muon_charge[i];
            temp.muon_pt =  muon_pt[i];
            temp.muon_e =  muon_e[i];
            temp.muon_eta =  muon_eta[i];
            temp.muon_phi =  muon_phi[i];
            x.push_back(temp);
        }
        return x;
    }, {"muon_charge", "muon_pt", "muon_e", "muon_eta", "muon_phi"})
    .Define("photons",[&](RVec<float>& photon_pt, RVec<float>& photon_e, RVec<float>& photon_eta, RVec<float>& photon_phi,  RVec<float>& photon_etcone40, RVec<int>& photon_id, RVec<int>& photon_id_loose, RVec<int>& photon_id_tight, RVec<float>& photon_cluster_eta_be_2 /*, RVec<int>& photon_id_nn */)
    {
        RVec<Photon> x;
        x.reserve(photon_pt.size());
        Photon temp;
        for (size_t i = 0; i < photon_pt.size(); i++)
        {
            temp.photon_pt =  photon_pt[i];
            temp.photon_e =  photon_e[i];
            temp.photon_eta =  photon_eta[i];
            temp.photon_phi =  photon_phi[i];
            temp.photon_etcone40 =  photon_etcone40[i];
            temp.photon_id =  photon_id[i];
            temp.photon_id_loose =  photon_id_loose[i];
            temp.photon_id_tight =  photon_id_tight[i];
            temp.photon_cluster_eta_be_2 =  photon_cluster_eta_be_2[i];
//            temp.photon_id_nn = photon_id_nn[i];
            x.push_back(temp);
        }
        return x;
    }, {"photon_pt", "photon_e", "photon_eta", "photon_phi", "photon_etcone40", "photon_id", "photon_id_loose", "photon_id_tight", "photon_cluster_eta_be_2", /*"photon_id_nn"*/})
    .Define("clusters",[&](RVec<float>& cluster_pt, RVec<float>& cluster_phi, RVec<float>& cluster_e, RVec<float>& cluster_eta)
    {
        RVec<Cluster> x;
        x.reserve(cluster_pt.size());
        Cluster temp;
        for (size_t i = 0; i < cluster_pt.size(); i++)
        {
            temp.cluster_pt =  cluster_pt[i];
            temp.cluster_e =  cluster_e[i];
            temp.cluster_eta =  cluster_eta[i];
            temp.cluster_phi =  cluster_phi[i];
            x.push_back(temp);
        }
        return x;
    }, {"cluster_pt", "cluster_phi", "cluster_e", "cluster_eta"})
    .Define("tracks",[&](RVec<int>& track_type, RVec<float>& track_pt, RVec<float>& track_eta, RVec<float>& track_phi, /*RVec<float>& track_e,*/ RVec<float>& track_charge, RVec<int>& track_num_pixel_hits, RVec<int>& track_num_sct_hits)
    {
        RVec<Track> x;
        x.reserve(track_pt.size());
        Track temp;
        int temp_track_type;
        for (size_t i = 0; i < track_pt.size(); i++)
        {
            if (track_type[i]==0)
            {
                temp.track_pt =  track_pt[i];
                temp.track_eta =  track_eta[i];
                temp.track_phi =  track_phi[i];
    //            temp.track_e =  track_e[i];
                temp.track_charge =  track_charge[i];
                temp.track_num_pixel_hits =  track_num_pixel_hits[i];
                temp.track_num_sct_hits =  track_num_sct_hits[i];
                x.push_back(temp);
            }
        }
        return x;
    }, {"track_type", "track_pt", "track_eta", "track_phi", /*"track_e",*/ "track_charge", "track_num_pixel_hits", "track_num_sct_hits"})
    .Vary("photons",[&](RVec<Photon>& photons, RVec<float>& photon_pt, RVec<std::vector<std::string>>& photon_syst_name, RVec<std::vector<float>>& photon_syst_pt, RVec<std::vector<float>>& photon_syst_e)
    {
        int index;
        RVec<RVec<Photon>> variedPhotons;
        variedPhotons.reserve(Event::systematics.size());
        RVec<Photon> x;
        auto length = photon_pt.size();
        x.reserve(length);
        Photon temp;
        
        for (auto& systematic: Event::systematics)
        {
            for (size_t i = 0; i < photon_pt.size(); i++)
            {
                auto it = std::find(photon_syst_name[i].begin(), photon_syst_name[i].end(), systematic);
                if (it == photon_syst_name[i].end())
                {
                    continue;
                }
                
                index = it - photon_syst_name[i].begin();
                temp.photon_pt = photon_syst_pt[i][index];
                temp.photon_e = photon_syst_e[i][index];
                temp.photon_eta =  photons[i].photon_eta;
                temp.photon_phi =  photons[i].photon_phi;
                temp.photon_etcone40 =  photons[i].photon_etcone40;
                temp.photon_id =  photons[i].photon_id;
                temp.photon_id_loose =  photons[i].photon_id_loose;
                temp.photon_id_tight =  photons[i].photon_id_tight;
                temp.photon_cluster_eta_be_2 =  photons[i].photon_cluster_eta_be_2;
    //            temp.photon_id_nn = photons[i].photon_id_nn;
                
                x.push_back(temp);
            }
            variedPhotons.push_back(x);
            x.clear();
        }
        
        return variedPhotons;
        
    }, {"photons", "photon_pt", "photon_syst_name", "photon_syst_pt", "photon_syst_e"}, Event::systematics)
    .Vary("electrons",[&](RVec<Electron>& electrons, RVec<float>& electron_pt, RVec<std::vector<std::string>>& electron_syst_name, RVec<std::vector<float>>& electron_syst_pt, RVec<std::vector<float>>& electron_syst_e)
    {
        int index;
        RVec<RVec<Electron>> variedElectrons;
        variedElectrons.reserve(Event::systematics.size());
        RVec<Electron> x;
        auto length = electron_pt.size();
        x.reserve(length);
        Electron temp;
        
        for (auto& systematic: Event::systematics)
        {
            for (size_t i = 0; i < electron_pt.size(); i++)
            {
                auto it = std::find(electron_syst_name[i].begin(), electron_syst_name[i].end(), systematic);
                if (it == electron_syst_name[i].end())
                {
                    continue;
                }
                
                index = it - electron_syst_name[i].begin();
                temp.electron_pt = electron_syst_pt[i][index];
                temp.electron_e = electron_syst_e[i][index];
                temp.electron_charge =  electrons[i].electron_charge;
                temp.electron_pt =  electrons[i].electron_pt;
                temp.electron_e =  electrons[i].electron_e;
                temp.electron_eta =  electrons[i].electron_eta;
                temp.electron_phi =  electrons[i].electron_phi;
    //            temp.electron_id =  electrons[i].electron_id;
                temp.electron_isolation =  electrons[i].electron_isolation;
                temp.electron_d0 =  electrons[i].electron_d0;
                temp.electron_z0 =  electrons[i].electron_z0;
    //            temp.electron_id_medium =  electrons[i].electron_id_medium;
                
                x.push_back(temp);
            }
            variedElectrons.push_back(x);
            x.clear();
        }
        
        return variedElectrons;
        
    }, {"electrons", "electron_pt", "electron_syst_name", "electron_syst_pt", "electron_syst_e"}, Event::systematics);
    
    return NewDf;
}

