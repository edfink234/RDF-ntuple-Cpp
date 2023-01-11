#ifndef MAKERDFH
#define MAKERDFH

#include <vector>
#include <string>
#include <memory>
#include "ROOT/RDataFrame.hxx"

using SchottDataFrame = ROOT::RDF::RInterface<ROOT::Detail::RDF::RLoopManager>;

class RDFTree;
SchottDataFrame MakeRDF(const std::vector<std::string>&, short numThreads = -1);
//TODO: Also explore RDatasetSpec: https://root.cern/doc/master/classROOT_1_1RDataFrame.html#a1fc7cc9a6eae595ccdc56b92cba21346
class RDFTree
{
private:
    static TChain* __chain;
    static TChain* __event_info_chain;
    friend SchottDataFrame MakeRDF(const std::vector<std::string>&, short numThreads);
};

#endif
