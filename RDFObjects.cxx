#include <string>
#include <sstream>

#include "RDFObjects.h"

const std::string TruthParticle::PREFIX = "mc";
const std::string Electron::PREFIX = "electron";
const std::string Photon::PREFIX = "photon";
const std::string Muon::PREFIX = "muon";
const std::string Track::PREFIX = "track";
const std::string Cluster::PREFIX = "cluster";

const int Electron::PDG_ID = 11;
const int Muon::PDG_ID = 13;
const int Photon::PDG_ID = 22;

std::string cling::printValue(TruthParticle *p)
{
    std::ostringstream os;
    os << "TruthParticle( pdg id = " << p->mc_pdg_id <<
    ", pT = " << p->mc_pt  << ", charge = " << p->mc_charge <<
    ", eta = " << p->mc_eta << ", phi = "    << p->mc_phi << " )";
   
    return os.str();
}

std::string cling::printValue(Electron *p)
{
    std::ostringstream os;
    os << "Electron( pdg id = " << p->electron_charge*-1*Electron::PDG_ID <<
    ", pT = " << p->electron_pt  << ", charge = " << p->electron_charge <<
    ", eta = " << p->electron_eta << ", phi = "   << p->electron_phi << " )";
   
    return os.str();
}

std::string cling::printValue(Muon *p)
{
    std::ostringstream os;
    os << "Muon( pdg id = " << p->muon_charge*-1*Muon::PDG_ID <<
    ", pT = " << p->muon_pt  << ", charge = " << p->muon_charge <<
    ", eta = " << p->muon_eta << ", phi = "   << p->muon_phi << " )";
   
    return os.str();
}

std::string cling::printValue(Photon *p)
{
    std::ostringstream os;
    os << "Photon( pdg id = " << Photon::PDG_ID <<
    ", pT = " << p->photon_pt  << ", eta = " << p->photon_eta
    << ", phi = "   << p->photon_phi << " )";
   
    return os.str();
}

std::string cling::printValue(Cluster *p)
{
    std::ostringstream os;
    os << "Cluster( pT = " << p->cluster_pt <<
    ", eta = " << p->cluster_eta  << ", phi = " << p->cluster_phi << " )";
   
    return os.str();
}

std::string cling::printValue(Track *p)
{
    std::ostringstream os;
    os << "Track( pT = " << p->track_pt <<
    ", eta = " << p->track_eta  << ", phi = " << p->track_phi
    << ", charge = " << p->track_charge << " )";
   
    return os.str();
}
