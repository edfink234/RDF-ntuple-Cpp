#ifndef RDFEVENTH
#define RDFEVENTH

#include <string>

namespace Event
{
    extern bool cache_truth;
    extern bool load_reco;
    extern bool load_photons;
    extern bool load_electrons;
    extern bool load_muons;
    extern bool load_clusters;
    extern bool load_tracks;
    extern bool load_triggers;
    extern std::string systematic;
};

#endif

