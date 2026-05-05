#include "./StJetTrackEvent.h"
#include "./StJetTrackEventLinkDef.h"
#include "./functions.h"
#include "TClassRef.h"
#include "Math/Vector4D.h"
#include <algorithm>
#include "./PythiaEvent.h"
#include "./PythiaEvent_LinkDef.h"
//#include "Fastjet/fastjet-3.5.1/include/fastjet/ClusterSequence.hh"
#include "Fastjet/fastjet-3.5.1/include/fastjet/ClusterSequenceArea.hh"
#include "Fastjet/fastjet-3.5.1/include/fastjet/config.h"
#include "Fastjet/fastjet-3.5.1/include/fastjet/PseudoJet.hh"
#include "Fastjet/fastjet-3.5.1/include/fastjet/Selector.hh"
#include "Fastjet/fastjet-3.5.1/include/fastjet/tools/Subtractor.hh"
//#include "StPhysicalHelixD.hh"
#include "Fastjet/fastjet-3.5.1/include/fastjet/tools/JetMedianBackgroundEstimator.hh"
using namespace fastjet;
using namespace std;

const char PYTHIA_TREE[]   = "PythiaEvent";
const char PYTHIA_EVENT_BRANCH[] = "Events";

vector<PseudoJet> FindJets(vector<PseudoJet> vec_particles);
double GetPhiByProb(void);
double GetEtaByProb(void);
double GetPtByProb(void);
double GetNByProb(void);

/*Parameters:
should be clear

Additional Info:

Event Plane is always at \phi = 0
*/
Int_t FlowMC(Int_t i_NoOfEvents = -1, double i_PtLowerLimit = 40)
{

    gStyle->SetOptStat(0);
    SetRootGraphicStyle();

    #define DEF_BinningPerUnit 100

    //other variables
    TH1D* h1D_JetMultiplicityInPsi = new TH1D("h1D_JetMultiplicityInPsi", "h1D_JetMultiplicityInPsi", DEF_BinningPerUnit*2*Pi, -Pi, Pi);
    TF1 *PhiDist=new TF1("PhiDist","(1/(2*pi))*(1+0.05*cos(2*x))",-Pi,Pi);
    TRandom TR_Eta;
    TR_Eta.SetSeed(0);

    TChain              *inputPYTHIA;
    PythiaEvent         *PYTHIAEvent;
    PythiaParticle      *PYTHIAParticle;
    Long64_t            file_entries;

    //open the root file
    TFile *file = TFile::Open("merge_mult_track_pT.root");

    if (!file || file->IsZombie()) {
        std::cout << "Error while opening the file!" << std::endl;
        return 0;
    }

    //get statistics histograms from root file
    TH1D* h1D_NDist = (TH1D*)file->Get("h_mult_particles_used_SE_ME_0");
    if(!h1D_NDist){cout << "Multiplicity histogram not found!" << endl; return 0;}
    TH1D* h1D_PtDist = (TH1D*)file->Get("h_particle_pT");
    if(!h1D_PtDist){cout << "Pt histogram not found!" << endl; return 0;}
    
    //------------------------------------
    // Read input data pythia tree
    TString pinputdirPYTHIA = "/media/nicolas-wirth/INTENSO/PYTHIA_pp/";
    TString SEListPYTHIA = "/filelist_PYTHIA_pp.txt";
    //TString SEListPYTHIA = "/filelist_PYTHIA_pp_OnlyHigh.txt";
    TString Pythia_List =  "/media/nicolas-wirth/INTENSO/PYTHIA_pp";
    Pythia_List+=SEListPYTHIA;
    Int_t file_loop = 0;
    if (!Pythia_List.IsNull())   // if input file is ok
    {
        cout << "Open file list " << Pythia_List << endl;
        ifstream in(Pythia_List);  // input stream
        if(in)
        {
            cout << "file list is ok" << endl;
            inputPYTHIA = new TChain(PYTHIA_TREE , PYTHIA_TREE );
            char str[255];       // char array for each file name
            Long64_t entries_save = 0;
            while(in)
            {
                in.getline(str,255);  // take the lines of the file list
                if(str[0] != 0)
                {
                    TString addfile = str;
                    addfile = pinputdirPYTHIA+addfile;
                    //cout << "addfile: " << addfile.Data() << endl;
                    Long64_t file_entries;////////////////////////////////////////////////////////////////////////////////////////////////////////
                    //if(file_loop == file_select) inputPYTHIA ->AddFile(addfile.Data(),-1, PYTHIA_TREE );
                    inputPYTHIA->AddFile(addfile.Data(), -1, PYTHIA_TREE);
                    file_entries = inputPYTHIA->GetEntries();
                    cout << "File added to data chain: " << addfile.Data() << " with " << file_entries-entries_save << " entries" << endl;
                    entries_save = file_entries;
                    file_loop++;
                }
            }
            file_entries   = inputPYTHIA->GetEntries();
            cout << "Number of entries in chain: " << file_entries << endl;

            PYTHIAEvent = new PythiaEvent();
            inputPYTHIA  ->SetBranchAddress( PYTHIA_EVENT_BRANCH, &PYTHIAEvent);
        }
    }
    //------------------------------------
    printf("StJetAnalysis::InitPythiaTree() finished \n");

    TString str_out = "./FlowMC_Output.root";
    TFile* OutputFile = new TFile(str_out.Data(),"RECREATE");
    OutputFile->cd();
    cout << "Output file is " << str_out << endl;
    TDirectory *EventHistos = OutputFile->mkdir("Event Histos");
    TDirectory *Results = OutputFile->mkdir("Results");
    
    if(i_NoOfEvents != -1) file_entries = i_NoOfEvents;

    //loop over events
    for(Int_t iEvent = 0; iEvent < file_entries; iEvent++)
    {

        vector<PseudoJet> ParticleVector;//holds all the particles in the event
        TH2D* h2D_phi_vs_eta = new TH2D("h2D_phi_vs_eta","h2D_phi_vs_eta", DEF_BinningPerUnit * 2 * Pi, -Pi, Pi, DEF_BinningPerUnit * 2 * 0.9, -0.9, 0.9);
        Int_t N = h1D_NDist->GetRandom();

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /*
        /
        /                                                     
        /                                            GENERATE PBPB BACKGROUND BASED ON PROBABILITY DIST WEIGHTED MC PARTICLES
        /                                                     
        /                                                     
        /                                                     
        /
        /                                                    
        /
        /
        /
        *//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
        for(Int_t iBGPart = 0; iBGPart < N; iBGPart++)
        {
            double Phi = PhiDist->GetRandom();
            double Eta = TR_Eta.Uniform(-0.9, 0.9);
            double Pt = h1D_PtDist->GetRandom();
            //fill histogram
            h2D_phi_vs_eta->Fill(Phi, Eta, Pt);

            //feed jetfinder
            double p_x = Pt/(sqrt(1+pow(tan(Phi), 2)));//ONLY FOR E CALCULATION - SIGNS NOT NECESSARILY CORRECT
            double p_y = sqrt(pow(Pt,2) - pow(p_x,2));//ONLY FOR E CALCULATION - SIGNS NOT NECESSARILY CORRECT
            double p_z = (pow(p_x,2) + pow(p_y,2))/((1/pow(tanh(Eta),2))-1);//ONLY FOR E CALCULATION - SIGNS NOT NECESSARILY CORRECT
            double E = sqrt(pow(p_x,2) + pow(p_y,2) + pow(p_z,2) + pow(0.1349766,2));
            ROOT::Math::PtEtaPhiEVector v(Pt, Eta, Phi, E);
            p_x = v.Px();
            p_y = v.Py();
            p_z = v.Pz();
            ParticleVector.push_back(PseudoJet(p_x, p_y, p_z, E));

        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /*
        /
        /                                                     
        /                                            OVERLAY WITH PYTHIA EVENT TO HAVE SIMULATED PBPB EVENT
        /                                                     
        /                                                     
        /                                                     
        /
        /                                                    
        /
        /
        /
        *//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
        /* 
        printf("StJetAnalysis::LoopEventPythia() started \n, ");
        cout << "------------------------iEvent PYTHIA: " << iEvent << endl;
        */
        if (!inputPYTHIA->GetEntry( iEvent )) // take the event -> information is stored in event
            cout << "Error" << endl;

        if (iEvent != 0  &&  iEvent % 100 == 0)
            cout << "." << flush;
        if (iEvent != 0  &&  iEvent % 1000 == 0)
        {
            if((file_entries-0) > 0)
            {
                Double_t event_percent = 100.0*((Double_t)(iEvent-0))/((Double_t)(file_entries-0));
                cout << " " << iEvent << " (" << event_percent << "%) " << "\n" << "==> Processing data " << flush;
            }
        }

        //--------------------
        // Event information
        Int_t   N_Particles     = PYTHIAEvent->getNumParticles();
        if(N_Particles <= 0)
        {
            printf("WARNING: iEvent: %d has no entries! \n",iEvent);
        }

        for(Int_t i_Particle = 0; i_Particle < N_Particles; i_Particle++)
        {    
            // Particle information
            // Particle Level
            PYTHIAParticle = PYTHIAEvent  ->getParticle(i_Particle);
            TLorentzVector TLV_pythia_particle = PYTHIAParticle -> get_TLV_part();

            Double_t pt_track  = TLV_pythia_particle.Pt();
            //if((fabs(pt_track) < 0.15)) continue;
            Double_t phi_track = TLV_pythia_particle.Phi();
            Double_t eta_track = TLV_pythia_particle.Eta();
            Double_t px_track  = TLV_pythia_particle.Px();
            Double_t py_track  = TLV_pythia_particle.Py();
            Double_t pz_track  = TLV_pythia_particle.Pz();
            Double_t E_track   = TLV_pythia_particle.E();

            // track cuts 
            if((fabs(pt_track) < 0.15)) continue;
            if((fabs(eta_track)) > 0.9) continue;

            //write particle into histogram
            h2D_phi_vs_eta->Fill(phi_track, eta_track, pt_track);

            ParticleVector.push_back(PseudoJet(px_track, py_track, pz_track, E_track));

        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /*
        /
        /                                                     
        /                                           FIND JETS WITH THE ANTI-KT ALGORITHM   
        /                                                     
        /                                                     
        /                                                     
        /
        /                                                    
        /
        /
        /
        *//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        vector<PseudoJet> JetVector = FindJets(ParticleVector);

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /*
        /
        /                                                     
        /                                           ONLY WRITE JETS THAT HAVE P_T > XXX INTO THE HISTOGRAM  
        /                                                     
        /                                                     
        /                                                     
        /
        /                                                    
        /
        /
        /
        *//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        for(Int_t i=0; i<JetVector.size(); i++)
        {
            if(JetVector[i].pt() >= i_PtLowerLimit) h1D_JetMultiplicityInPsi->Fill(JetVector[i].phi_std());
        }

        //write event results
        EventHistos->cd();
        h2D_phi_vs_eta->Write();
        delete h2D_phi_vs_eta;
    }

    //write global results
    Results->cd();
    h1D_JetMultiplicityInPsi->Write();
    OutputFile ->Close();
    
    return 1;
}


vector<PseudoJet> FindJets(vector<PseudoJet> vec_particles)
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
    vector<PseudoJet> jets_fiducial;//contains the selected jets in the end

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

    /*
    Int_t Rem_n_hardest = 2; 
    Int_t eRem_n_hardest = 2;////???????
    double eBkg_R = 0.4;
    JetDefinition jet_def_bkgd(kt_algorithm, eBkg_R); // <--
    AreaDefinition area_def_bkgd(active_area_explicit_ghosts,GhostedAreaSpec(ghost_maxrap,1,0.01));
    Int_t Rem_n_hardest_use = eRem_n_hardest;
    Selector selector = SelectorAbsEtaMax(0.9 - eBkg_R) * (!SelectorNHardest(Rem_n_hardest)); // <--

    JetMedianBackgroundEstimator bkgd_estimator(selector, jet_def_bkgd, area_def_bkgd); // <--
    Subtractor subtractor(&bkgd_estimator);
    bkgd_estimator.set_particles(vec_particles);
    Double_t jet_rho   = bkgd_estimator.rho();
    */

    // loop over jets  
    Double_t Jet_area_cut = 0.65*TMath::Pi()*TMath::Power(jet_radius,2);
    for(Int_t i_jet = 0; i_jet < (Int_t)jets_fiducial.size(); i_jet++)
    {
        Float_t jet_pt         = jets_fiducial[i_jet].perp();
        if(jet_pt <= 0.0) continue;
        Float_t jet_area       = jets_fiducial[i_jet].area();
        //Float_t jet_pt_sub     = jets_fiducial[i_jet].perp() - jet_rho*jet_area;
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

    return jets_fiducial;

}
