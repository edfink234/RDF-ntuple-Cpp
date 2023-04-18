#ifndef MAKERDFH
#define MAKERDFH

#include <vector>
#include <string>
#include <memory>
#include "ROOT/RDataFrame.hxx"

//Alias for head RNode; return type of MakeRDF
using SchottDataFrame = ROOT::RDF::RInterface<ROOT::Detail::RDF::RLoopManager>;

SchottDataFrame MakeRDF(const std::vector<std::string>&, short numThreads = -1);

#endif
