#include "fastjet-3.5.1/include/fastjet/ClusterSequenceArea.hh"
#include <iostream>
#include "./StJetTrackEvent.h"
#include "./StJetTrackEventLinkDef.h"
#include "./functions.h"
#include <algorithm>
#include "fastjet-3.5.1/include/fastjet/config.h"
#include "fastjet-3.5.1/include/fastjet/PseudoJet.hh"
#include "fastjet-3.5.1/include/fastjet/Selector.hh"
#include "fastjet-3.5.1/include/fastjet/tools/Subtractor.hh"
//#include "StPhysicalHelixD.hh"
#include "fastjet-3.5.1/include/fastjet/tools/JetMedianBackgroundEstimator.hh"
using namespace fastjet;
using namespace std;

void JetFinder()
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /*
    /
    /                                                           JET FINDER CODE 
    /
    /
    /
    /
    *//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // To do: define jet radius, background radius etc
    #define jet_radius 0.4

    // define variables 
    PseudoJet jets_fiducial;//contains the selected jets in the end

    Double_t ghost_maxrap = 2.5; // Fiducial cut for background estimation
    Double_t eta_acceptance = 2.5;
    // choose a jet definition
    JetDefinition jet_def(antikt_algorithm, jet_radius);
    // jet area definition
    GhostedAreaSpec area_spec(ghost_maxrap);
    AreaDefinition area_def(active_area_explicit_ghosts,GhostedAreaSpec(ghost_maxrap,1,0.01));

    ClusterSequenceArea clust_seq_hard(vec_particles, jet_def, area_def);
    // run the clustering, extract the jets
    double ptmin = 0.15;
    vector<PseudoJet> jets_all = sorted_by_pt(clust_seq_hard.inclusive_jets(ptmin));
    Selector Fiducial_cut_selector = SelectorAbsEtaMax(eta_acceptance - jet_radius); // Fiducial cut for jets
    jets_fiducial = Fiducial_cut_selector(jets_all);

    // background estimation
    // Standard rho calculation with kT algorithm and N hardest jets removed

    Int_t Rem_n_hardest = 2; 
    JetDefinition jet_def_bkgd(kt_algorithm, eBkg_R); // <--
    AreaDefinition area_def_bkgd(active_area_explicit_ghosts,GhostedAreaSpec(ghost_maxrap,1,0.01));
    Int_t Rem_n_hardest_use = eRem_n_hardest;
    Selector selector = SelectorAbsEtaMax(0.9 - eBkg_R) * (!SelectorNHardest(Rem_n_hardest)); // <--

    JetMedianBackgroundEstimator bkgd_estimator(selector, jet_def_bkgd, area_def_bkgd); // <--
    Subtractor subtractor(&bkgd_estimator);
    bkgd_estimator.set_particles(vec_PJ_particles);
    Double_t jet_rho   = bkgd_estimator.rho();


    // loop over jets  
    Double_t Jet_area_cut = 0.65*TMath::Pi()*TMath::Power(jet_radius,2);
    for(Int_t i_jet = 0; i_jet < (Int_t)jets_fiducial.size(); i_jet++)
    {
        Float_t jet_pt         = jets_fiducial[i_jet].perp();
        if(jet_pt <= 0.0) continue;
        Float_t jet_area       = jets_fiducial[i_jet].area();
        Float_t jet_pt_sub     = jets_fiducial[i_jet].perp() - jet_rho*jet_area;
        Float_t jet_eta        = jets_fiducial[i_jet].eta();
        Float_t jet_phi        = jets_fiducial[i_jet].phi();
        Float_t jet_phi_std    = jets_fiducial[i_jet].phi_std();
        if(jet_area < Jet_area_cut) continue; 


        // constituent track loop
        vector<PseudoJet> jet_constituents = jets_fiducial[i_jet].constituents();

        for(Int_t i_constituent = 0; i_constituent < (Int_t)jet_constituents.size(); i_constituent++)
        {
            Float_t jet_const_pt  = jet_constituents[i_constituent].perp();
            if(jet_const_pt <= 0.0) continue;
            Float_t jet_const_phi = jet_constituents[i_constituent].phi();
            Float_t jet_const_rap = jet_constituents[i_constituent].rap();
            Float_t jet_const_Phi = jet_constituents[i_constituent].phi_std();
            Float_t jet_const_eta = jet_constituents[i_constituent].eta();
            Int_t   user_index    = jet_constituents[i_constituent].user_index();

        }
    }

}


int main () {
  vector<PseudoJet> particles;
  // an event with three particles:   px    py  pz      E
  particles.push_back( PseudoJet(   99.0,  0.1,  0, 100.0) ); 
  particles.push_back( PseudoJet(    4.0, -0.1,  0,   5.0) ); 
  particles.push_back( PseudoJet(  -99.0,    0,  0,  99.0) );

  // choose a jet definition
  double R = 0.7;
  JetDefinition jet_def(antikt_algorithm, R);

  // run the clustering, extract the jets
  ClusterSequence cs(particles, jet_def);
  vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());

  // print out some infos
  cout << "Clustering with " << jet_def.description() << endl;

  // print the jets
  cout <<   "        pt y phi" << endl;
  for (unsigned i = 0; i < jets.size(); i++) {
    cout << "jet " << i << ": "<< jets[i].pt() << " " 
                   << jets[i].rap() << " " << jets[i].phi() << endl;
    vector<PseudoJet> constituents = jets[i].constituents();
    for (unsigned j = 0; j < constituents.size(); j++) {
      cout << "    constituent " << j << "'s pt: " << constituents[j].pt()
           << endl;
    }
  }
} 