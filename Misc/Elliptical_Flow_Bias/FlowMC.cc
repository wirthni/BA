#include "./StJetTrackEvent.h"
#include "./StJetTrackEventLinkDef.h"
#include "./functions.h"
#include "TClassRef.h"
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

//returns phi of the particle based on the elliptical flow probablity distribution
double GetPhiByProb(void)
{
    TRandom TR_Phi;
    TR_Phi.SetSeed(0);
    return TR_Phi.Uniform(-Pi, Pi);
}
double GetEtaByProb(void)
{
    TRandom TR_Eta;
    TR_Eta.SetSeed(0);
    return TR_Eta.Uniform(-0.9, 0.9);
}
double GetNByProb(void)
{
    TRandom TR_N;
    TR_N.SetSeed(0);
    return TR_N.Uniform(3000, 4000);
}
double GetPtByProb(void)
{
    TRandom TR_Pt;
    TR_Pt.SetSeed(0);
    return TR_Pt.Uniform(20, 200);
}

vector<PseudoJet> FindJets(vector<PseudoJet> vec_particles);

/*Parameters:

*/
void FlowMC(void)
{

    #define DEF_BinningPerUnit 200

    TChain              *inputPYTHIA;
    PythiaEvent         *PYTHIAEvent;
    PythiaParticle      *PYTHIAParticle;
    Long64_t            file_entries;

    
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

    TString str_out = "./CMS_PYTHIA_pp_Output_Cuts_30_60_MinSep1_5.root";
    TFile* OutputFile = new TFile(str_out.Data(),"RECREATE");
    OutputFile->cd();
    cout << "Output file is " << str_out << endl;
    
    //loop over events
    for(Int_t iEvent = 0; iEvent < file_entries; iEvent++)
    {

        vector<PseudoJet> ParticleVector;//holds all the particles in the event

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
        
        TH2D* h2D_phi_vs_eta = new TH2D("h2D_phi_vs_eta","h2D_phi_vs_eta", DEF_BinningPerUnit * 2 * 0.9, -0.9, 0.9, DEF_BinningPerUnit * 2 * Pi, -Pi, Pi);
        Int_t N = GetNByProb();

        for(Int_t iBGPart = 0; iBGPart < N; iBGPart++)
        {
            h2D_phi_vs_eta->Fill(GetPhiByProb(), GetEtaByProb(), GetPtByProb());
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

            // track cuts 
            if((fabs(pt_track) < 0.15)) continue;
            if((fabs(eta_track)) > 0.9) continue;

            //write particle into histogram
            h2D_phi_vs_eta->Fill(phi_track, eta_track, pt_track);

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



        //DELETE HISTOGRAMS
    }
    
}

