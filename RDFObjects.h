#ifndef RDFOBJECTSH
#define RDFOBJECTSH

#include <string>

#include "Math/Vector4D.h"
#include "ROOT/RDataFrame.hxx"

using namespace ROOT::Math;

struct TruthParticle final
{
    static const std::string PREFIX;
    
    int   mc_pdg_id;
    int   mc_barcode;
    int   mc_parent_barcode;
    int   mc_status;
    float mc_pt;
    float mc_charge;
    float mc_eta;
    float mc_phi;
    float mc_e;
    float mc_mass;
    
    PtEtaPhiEVector Vector()
    {
        return PtEtaPhiEVector(mc_pt, mc_eta, mc_phi, mc_e);
    }
};

struct Electron final
{
    static const std::string PREFIX;
    static const int PDG_ID;
    
    float electron_charge;
    float electron_pt;
    float electron_e;
    float electron_eta;
    float electron_phi;
//    int   electron_id;
    float electron_isolation;
    float electron_d0;
    float electron_z0;
//    int   electron_id_medium;
    
    PtEtaPhiEVector Vector()
    {
        return PtEtaPhiEVector(electron_pt, electron_eta, electron_phi, electron_e);
    }
};

struct Muon final
{
    static const std::string PREFIX;
    static const int PDG_ID;
    
    int muon_charge;
    float muon_pt;
    float muon_e;
    float muon_eta;
    float muon_phi;
    
    PtEtaPhiEVector Vector()
    {
        return PtEtaPhiEVector(muon_pt, muon_eta, muon_phi, muon_e);
    }
};

struct Photon final
{
    static const std::string PREFIX;
    static const int PDG_ID;

    float photon_pt;
    float photon_e;
    float photon_eta;
    float photon_phi;
    float photon_etcone40;
    int   photon_id;
    int   photon_id_loose;
    int   photon_id_tight;
    float photon_cluster_eta_be_2;
    int   photon_id_nn;
    
    PtEtaPhiEVector Vector()
    {
        return PtEtaPhiEVector(photon_pt, photon_eta, photon_phi, photon_e);
    }
};

struct Cluster final
{
    static const std::string PREFIX;
    
    float cluster_pt;
    float cluster_eta;
    float cluster_phi;
    float cluster_e;
    
    PtEtaPhiEVector Vector()
    {
        return PtEtaPhiEVector(cluster_pt, cluster_eta, cluster_phi, cluster_e);
    }
};

struct Track final
{
    static const std::string PREFIX;
    
    float track_pt;
    float track_charge;
    float track_eta;
    float track_phi;
//    float track_e;
    int   track_num_pixel_hits;
    int   track_num_sct_hits;
    
    PtEtaPhiEVector Vector()
    {
        return PtEtaPhiEVector(track_pt, track_eta, track_phi, 0 /*track_e*/);
    }
};

namespace cling
{
    std::string printValue(TruthParticle *p);
    std::string printValue(Electron *p);
    std::string printValue(Muon *p);
    std::string printValue(Photon *p);
    std::string printValue(Cluster *p);
    std::string printValue(Track *p);
}

#endif

