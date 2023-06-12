/*
 Run on Schott-lab-bsm! Here's how to set it up!
 1. git clone -b Schott-lab-branch https://github.com/edfink234/RDF-ntuple-Cpp.git
 2. cd RDF-ntuple-Cpp/
 3. source /cvmfs/sft.cern.ch/lcg/views/LCG_102/x86_64-centos7-gcc11-opt/setup.sh
 4. root --version
    a. Should say something like the following:
        ROOT Version: 6.26/04
        Built for linuxx8664gcc on Jun 07 2022, 16:01:16
        From tags/v6-26-04@v6-26-04
 5. In the directory RDF-ntuple-Cpp/, run the following two commands:
    a. rootcling -v4 -f mydict.cxx  -rmf libmydict.rootmap -rml libmydict.so  LinkDef.h
    b. g++ -shared -o libmydict.so mydict.cxx `root-config --cflags --libs` -fPIC
 6. Now compile with the following:
    a. g++ -g -o exampleAnalysis exampleAnalysis.cxx MakeRDF.cxx RDFObjects.cxx RDFevent.cxx -O2 $(root-config --libs --cflags) -std=c++17
 7. And run:
    a. ./exampleAnalysis
 */

/*
 Analysis: Reconstruct Higgs boson invariant mass from channel H->Za->llyy
 */

#include "MakeRDF.h"
#include "RDFObjects.h"
#include "RDFevent.h"
#include "ROOT/RDFHelpers.hxx" //VariationsFor
#include "TH1.h"   // If you are using histograms
#include "TH1D.h"
#include "TFile.h"
#include <iostream>

using ROOT::VecOps::RVec;
using ROOT::RDF::Experimental::VariationsFor;

void exampleAnalysis()
{
    //Create RDataFrame
//    std::string file = "Ntuple_MC_Za_mA5p0_v4.root";
    std::string file = "/home/common/Za/NTuples/Ntuple_MC_Za_m5p0_v4.root";
    Event::systematics = {"EG_RESOLUTION_ALL__1up", "EG_RESOLUTION_ALL__1down",}; //systematic variations
    SchottDataFrame df(MakeRDF({file}, 24)); //Make df with 24 threads: Schott-lab-bsm has 24 processors!
    
    //Transformations
    //---------------
    
    //electron selection
    auto selected_electrons = df.Define("selected_electrons", //create a column for selected reco-electrons
    [](RVec<AbstractParticle> electrons) //pass by value because we want to return a modified version
    {
        electrons.erase(std::remove_if(electrons.begin(),electrons.end(),
        [](AbstractParticle& electron)
        {
            //remove electron if it has |η| > 2.5, p_T <= 20 GeV, or doesn't pass a medium id criteria
            return std::abs(electron.electron_eta) > 2.5 or electron.electron_pt <= 20e3 or not electron.electron_id_medium;
            
        }), electrons.end());
        
        return electrons; //return the selected electrons in each event
        
    }, {"abstract_electrons"});
    
    auto electron_cut = selected_electrons.Filter( //Get rid of events that don't have exactly two opp-charged reco-electrons
    [](const RVec<AbstractParticle>& electrons)
    {
        return electrons.size() == 2 and electrons[0].electron_charge*electrons[1].electron_charge < 0;
    }, {"selected_electrons"});

    auto electron_mass_cut = electron_cut
    .Define("dielectron",
    [](RVec<AbstractParticle>& electrons)
    {
        return (electrons[0].ElectronVector() + electrons[1].ElectronVector());
    }, {"selected_electrons"})
    .Define("di_electron_inv_mass", //from the remaining events, create a new column di_electron_inv_mass
    [](PtEtaPhiEVector& dilep)
    {
        return dilep.M() / 1e3;
    }, {"dielectron"})
    .Filter( //Then get rid of events that don't have 66 GeV <= m_ee <= 116 GeV
    [](double mass)
    {
        return (66 <= mass and mass <= 116);
    }, {"di_electron_inv_mass"});
    
    //photon selectrion
    auto selected_photons = electron_mass_cut.Define("selected_photons", //create a column for selected reco-photons
    [](RVec<AbstractParticle> photons) //pass by value because we want to return a modified version
    {
        photons.erase(std::remove_if(photons.begin(),photons.end(),
        [](AbstractParticle& photon)
        {
            //remove photon if it has |η| > 2.5, p_T <= 10 GeV, or doesn't pass a loose id criteria
            return std::abs(photon.photon_eta) > 2.5 or photon.photon_pt <= 10e3 or not photon.photon_id_loose;
            
        }), photons.end());
        
        return photons; //return the selected photons in each event
        
    }, {"abstract_photons"});
    
    auto photon_cut = selected_photons.Filter( //Get rid of events that don't have exactly two opp-charged reco-photons
    [](const RVec<AbstractParticle>& photons)
    {
        return photons.size() == 2;
    }, {"selected_photons"});

    auto photon_mass = photon_cut.Define("diphoton", //create a new column for the diphoton 4-vector
    [](RVec<AbstractParticle>& photons)
    {
        return (photons[0].PhotonVector() + photons[1].PhotonVector());
    }, {"selected_photons"})
    .Define("di_photon_inv_mass", //using the diphoton 4-vector, create a new column di_photon_inv_mass
    [](PtEtaPhiEVector& diphoton)
    {
        return diphoton.M() / 1e3;
    }, {"diphoton"});
    
    auto di_photon_plus_lep = photon_mass.Define("higgs_mass", //create a new column for the higgs boson mass
    [](PtEtaPhiEVector& dielectron, PtEtaPhiEVector& diphoton)
    {
        return (dielectron + diphoton).M() / 1e3;
    }, {"dielectron", "diphoton"})
    .Filter( //and keep events where the invariant higgs mass is between 110 and 140 GeV.
    [](double higgs_mass)
    {
        return higgs_mass >= 110 and higgs_mass <= 140;
    }, {"higgs_mass"});

    //Book actions
    //------------
    
    auto Z_boson_mass = electron_mass_cut.Histo1D<double>({"Reconstructed_Z_boson_mass", "Reconstructed_Z_boson_mass", 60, 66, 116}, {"di_electron_inv_mass"});
    auto ALP_mass = photon_mass.Histo1D<double>({"Reconstructed_ALP_mass", "Reconstructed_ALP_mass", 60, 0, 10}, {"di_photon_inv_mass"});
    auto higgs_mass = di_photon_plus_lep.Histo1D<double>({"Reconstructed_Higgs_mass", "Reconstructed_Higgs_mass", 60, 110, 140}, {"higgs_mass"});
    
    auto varied_higgs_masses = VariationsFor(higgs_mass);
    
    //Access results
    //--------------
    Z_boson_mass->Draw();
    TFile* output_file = TFile::Open("exampleAnalysis_output.root", "RECREATE");
    Z_boson_mass->Write();
    ALP_mass->Write();
    
    for (auto& varied_higgs_mass: varied_higgs_masses.GetKeys())
    {
        varied_higgs_masses[varied_higgs_mass].SetTitle(varied_higgs_mass.c_str());
        varied_higgs_masses[varied_higgs_mass].Write();
    }
    output_file->Close();
}

int main()
{
    exampleAnalysis();
    return 0;
}
