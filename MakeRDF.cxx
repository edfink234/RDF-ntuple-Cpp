#include <vector>
#include <string>
#include <algorithm>

#include "TChain.h"
#include "TH1F.h"
#include "TBranch.h"
#include "TInterpreter.h"
#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RDF/RDatasetSpec.hxx"

#include "RDFObjects.h"
#include "MakeRDF.h"
#include "RDFevent.h"

using namespace ROOT::VecOps;

//Function to make an RVec of photons. We define it outside the call
//to Define because we want to use it more than once
template <typename T>
RVec<T> MakePhotons (const RVec<float>& photon_pt, const RVec<float>& photon_e, const RVec<float>& photon_eta, const RVec<float>& photon_phi, const RVec<float>& photon_etcone40, const RVec<int>& photon_id, const RVec<int>& photon_id_loose, const RVec<int>& photon_id_tight, const RVec<float>& photon_cluster_eta_be_2 /*, const RVec<int>& photon_id_nn */)
{
    RVec<T> x;
    x.reserve(photon_pt.size());
    T temp;
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
}

//Function to make an RVec of electrons. We define it outside the call
//to Define because we want to use it more than once
template <typename T>
RVec<T> MakeElectrons (const RVec<float>& electron_charge, const RVec<float>& electron_pt, const RVec<float>& electron_e, const RVec<float>& electron_eta, const RVec<float>& electron_phi, /*const RVec<float>& electron_id,*/ const RVec<float>& electron_isolation, const RVec<float>& electron_d0, const RVec<float>& electron_z0 , const RVec<int>& electron_id_medium)
{
    RVec<T> x;
    x.reserve(electron_pt.size());
    T temp;
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
}

//Function to register systematic variations for photons. We define it outside the call
//to Vary because we want to use it more than once
template <typename T>
RVec<RVec<T>> VaryPhotons(const RVec<T>& photons, const RVec<float>& photon_pt, const RVec<std::vector<std::string>>& photon_syst_name, const RVec<std::vector<float>>& photon_syst_pt, const RVec<std::vector<float>>& photon_syst_e)
{
    int index;
    RVec<RVec<T>> variedPhotons;
    variedPhotons.reserve(Event::systematics.size());
    RVec<T> x;
    auto length = photon_pt.size();
    x.reserve(length);
    T temp;
    
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
}

//Function to register systematic variations for electrons. We define it outside the call
//to Vary because we want to use it more than once
template <typename T>
RVec<RVec<T>> VaryElectrons(const RVec<T>& electrons, const RVec<float>& electron_pt, const RVec<std::vector<std::string>>& electron_syst_name, const RVec<std::vector<float>>& electron_syst_pt, const RVec<std::vector<float>>& electron_syst_e)
{
    int index;
    RVec<RVec<T>> variedElectrons;
    variedElectrons.reserve(Event::systematics.size());
    RVec<T> x;
    auto length = electron_pt.size();
    x.reserve(length);
    T temp;
    
    for (auto& systematic: Event::systematics)
    {
        for (size_t i = 0; i < electron_pt.size(); i++)
        {
            auto it = std::find(electron_syst_name[i].begin(), electron_syst_name[i].end(), systematic);
            if (it == electron_syst_name[i].end()) //not found
            {
                continue;
            }
            
            index = it - electron_syst_name[i].begin();
            temp.electron_pt = electron_syst_pt[i][index];
            temp.electron_e = electron_syst_e[i][index];
            temp.electron_charge =  electrons[i].electron_charge;
            temp.electron_eta =  electrons[i].electron_eta;
            temp.electron_phi =  electrons[i].electron_phi;
//            temp.electron_id =  electrons[i].electron_id;
            temp.electron_isolation =  electrons[i].electron_isolation;
            temp.electron_d0 =  electrons[i].electron_d0;
            temp.electron_z0 =  electrons[i].electron_z0;
            temp.electron_id_medium =  electrons[i].electron_id_medium;
            
            x.push_back(temp);
        }
        variedElectrons.push_back(x);
        x.clear();
    }
    
    return variedElectrons;
}

/*
 Creates an RDataFrame out of the physics and full_event_info TTrees,
 then adds columns for all of the objects defined in RDFObjects.h.
 The resulting dataframe is returned.
*/
SchottDataFrame MakeRDF(const std::vector<std::string>& files, short numThreads)
{
    if (numThreads > 0) //Enabling Multi-threading if specified
    {
        ROOT::EnableImplicitMT(numThreads);
    }
    static bool loaded = false;
    if (!loaded) //Only load objects once, or else an error occurs
    {
//        Uncomment the following line!
//        gInterpreter->LoadMacro("RDFObjects.h");
        
//        Delete the following line!
        gInterpreter->LoadMacro("RDFObjects.h");

//        Then load printValue overloads so objects can be printed
        gInterpreter->Declare("std::string cling::printValue(TruthParticle *);");
        gInterpreter->Declare("std::string cling::printValue(Electron *);");
        gInterpreter->Declare("std::string cling::printValue(Muon *);");
        gInterpreter->Declare("std::string cling::printValue(Photon *);");
        gInterpreter->Declare("std::string cling::printValue(Cluster *);");
        gInterpreter->Declare("std::string cling::printValue(Track *);");
        loaded = true;
    }
    
    ROOT::RDF::Experimental::RDatasetSpec x;
    ROOT::RDF::Experimental::RSample r("myPhysicsDf", "physics", files);
    x.AddSample(r);
    x.WithGlobalFriends("full_event_info", files);
    ROOT::RDataFrame df(x);

//    Now, create a new RDF representation that also holds the columns
//    for the objects, as well as Variations (if we choose to use them),
//    that are created using Define and Vary respectively
    auto NewDf = df.Define("truth_particles",[&](const RVec<int>& mc_pdg_id, const RVec<int>& mc_barcode, const RVec<int>& mc_parent_barcode, const RVec<int>& mc_status, const RVec<float>& mc_pt, const RVec<float>& mc_charge, const RVec<float>& mc_eta, const RVec<float>& mc_phi, const RVec<float>& mc_e, const RVec<float>& mc_mass)
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
    .Define("electrons",MakeElectrons<Electron>, {"electron_charge", "electron_pt", "electron_e", "electron_eta", "electron_phi", /*"electron_id",*/ "electron_isolation", "electron_d0", "electron_z0", "electron_id_medium",})
    .Define("muons",[&](const RVec<int>& muon_charge, const RVec<float>& muon_pt, const RVec<float>& muon_e, const RVec<float>& muon_eta, const RVec<float>& muon_phi)
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
    .Define("photons",MakePhotons<Photon>, {"photon_pt", "photon_e", "photon_eta", "photon_phi", "photon_etcone40", "photon_id", "photon_id_loose", "photon_id_tight", "photon_cluster_eta_be_2", /*"photon_id_nn"*/})
    .Define("clusters",[&](const RVec<float>& cluster_pt, const RVec<float>& cluster_phi, const RVec<float>& cluster_e, const RVec<float>& cluster_eta)
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
    .Define("tracks",[&](const RVec<int>& track_type, const RVec<float>& track_pt, const RVec<float>& track_eta, const RVec<float>& track_phi, /*const RVec<float>& track_e,*/ const RVec<float>& track_charge, const RVec<int>& track_num_pixel_hits, const RVec<int>& track_num_sct_hits)
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
    .Define("abstract_photons",MakePhotons<AbstractParticle>, {"photon_pt", "photon_e", "photon_eta", "photon_phi", "photon_etcone40", "photon_id", "photon_id_loose", "photon_id_tight", "photon_cluster_eta_be_2", /*"photon_id_nn"*/})
    .Define("abstract_electrons",MakeElectrons<AbstractParticle>, {"electron_charge", "electron_pt", "electron_e", "electron_eta", "electron_phi", /*"electron_id",*/ "electron_isolation", "electron_d0", "electron_z0", "electron_id_medium",})
    .Vary("photons", VaryPhotons<Photon>, {"photons", "photon_pt", "photon_syst_name", "photon_syst_pt", "photon_syst_e"}, Event::systematics)
    .Vary("electrons", VaryElectrons<Electron>, {"electrons", "electron_pt", "electron_syst_name", "electron_syst_pt", "electron_syst_e"}, Event::systematics)
    .Vary({"abstract_photons", "abstract_electrons"},
    [](const RVec<AbstractParticle>& photons, const RVec<float>& photon_pt, const RVec<std::vector<std::string>>& photon_syst_name, const RVec<std::vector<float>>& photon_syst_pt, const RVec<std::vector<float>>& photon_syst_e, const RVec<AbstractParticle>& electrons, const RVec<float>& electron_pt, const RVec<std::vector<std::string>>& electron_syst_name, const RVec<std::vector<float>>& electron_syst_pt, const RVec<std::vector<float>>& electron_syst_e)
    {
        return RVec<RVec<RVec<AbstractParticle>>>
        {
            VaryPhotons<AbstractParticle>(photons, photon_pt, photon_syst_name, photon_syst_pt, photon_syst_e),
            VaryElectrons<AbstractParticle>(electrons, electron_pt, electron_syst_name, electron_syst_pt, electron_syst_e)
        };
    }, {"abstract_photons", "photon_pt", "photon_syst_name", "photon_syst_pt", "photon_syst_e", "abstract_electrons", "electron_pt", "electron_syst_name", "electron_syst_pt", "electron_syst_e"}, Event::systematics, "photons_and_electrons")
    .Vary("photon_id_eff",[&](const RVec<float>& photon_id_eff, const RVec<std::vector<std::string>>& photon_syst_name, const RVec<std::vector<float>>& photon_syst_id_eff)
    {
        RVec<RVec<float>> varied_photon_id_eff;
        int index;
        varied_photon_id_eff.reserve(Event::systematics.size());
        RVec<float> x;
        auto length = photon_id_eff.size();
        x.reserve(length);
        float temp;
        
        for (auto& systematic: Event::systematics)
        {
            for (size_t i = 0; i < photon_id_eff.size(); i++)
            {
                auto it = std::find(photon_syst_name[i].begin(), photon_syst_name[i].end(), systematic);
                if (it == photon_syst_name[i].end())
                {
                    continue;
                }
                
                index = it - photon_syst_name[i].begin();
                temp = photon_syst_id_eff[i][index];
                x.push_back(temp);
            }
            varied_photon_id_eff.push_back(x);
            x.clear();
        }
        return varied_photon_id_eff;
    }, {"photon_id_eff", "photon_syst_name", "photon_syst_id_eff"}, Event::systematics)
    .Vary("photon_iso_eff",[&](const RVec<float>& photon_iso_eff, const RVec<std::vector<std::string>>& photon_syst_name, const RVec<std::vector<float>>& photon_syst_iso_eff)
    {
        RVec<RVec<float>> varied_photon_iso_eff;
        int index;
        varied_photon_iso_eff.reserve(Event::systematics.size());
        RVec<float> x;
        auto length = photon_iso_eff.size();
        x.reserve(length);
        float temp;
        
        for (auto& systematic: Event::systematics)
        {
            for (size_t i = 0; i < photon_iso_eff.size(); i++)
            {
                auto it = std::find(photon_syst_name[i].begin(), photon_syst_name[i].end(), systematic);
                if (it == photon_syst_name[i].end())
                {
                    continue;
                }
                
                index = it - photon_syst_name[i].begin();
                temp = photon_syst_iso_eff[i][index];
                x.push_back(temp);
            }
            varied_photon_iso_eff.push_back(x);
            x.clear();
        }
        return varied_photon_iso_eff;
    }, {"photon_iso_eff", "photon_syst_name", "photon_syst_iso_eff"}, Event::systematics)
    .Vary("photon_trg_eff",[&](const RVec<float>& photon_trg_eff, const RVec<std::vector<std::string>>& photon_syst_name, const RVec<std::vector<float>>& photon_syst_trg_eff)
    {
        RVec<RVec<float>> varied_photon_trg_eff;
        int index;
        varied_photon_trg_eff.reserve(Event::systematics.size());
        RVec<float> x;
        auto length = photon_trg_eff.size();
        x.reserve(length);
        float temp;
        
        for (auto& systematic: Event::systematics)
        {
            for (size_t i = 0; i < photon_trg_eff.size(); i++)
            {
                auto it = std::find(photon_syst_name[i].begin(), photon_syst_name[i].end(), systematic);
                if (it == photon_syst_name[i].end())
                {
                    continue;
                }
                
                index = it - photon_syst_name[i].begin();
                temp = photon_syst_trg_eff[i][index];
                x.push_back(temp);
            }
            varied_photon_trg_eff.push_back(x);
            x.clear();
        }
        return varied_photon_trg_eff;
    }, {"photon_trg_eff", "photon_syst_name", "photon_syst_trg_eff"}, Event::systematics);
    
    return NewDf; //Return the RDF with all the info we want!
}
