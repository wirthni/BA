#include "./ReadTree.h"
#include <Math/Vector4D.h>
#include <unordered_map>
#include <vector>
#include "functions.h"
#include "StJetTrackEvent.h"
#include "StJetTrackEventLinkDef.h"

Int_t ReadTree(TString InputFile, TString i_OutputFile)
{
    #define DEF_BinningPerUnit 100

    cout << "ReadTree started" << endl;

    // Make histograms for a global event analysis
    TH1D* h1d_mult_event   = new TH1D("h1d_mult_event","h1d_mult_event", 500,0,5000);
    TH1D* h1d_occ_event    = new TH1D("h1d_occ_event","h1d_occ_event", 10000,0,10000);
    TH1D* h1d_occft0_event = new TH1D("h1d_occft0_event","h1d_occft0_event", 100000,0,100000);
    TH1D* h1d_pTsum        = new TH1D("h1d_pTsum","h1d_pTsum", 2000,0,2000);
    TH1D* h1d_TOF          = new TH1D("h1d_TOF","h1d_TOF", 60000,-3000000,3000000);
    TH1D* h1d_cent         = new TH1D("h1d_cent","h1d_cent", 80,0,40);
    TH1D* h1d_psi2         = new TH1D("h1d_psi2","h1d_psi2", 800,-4.0,4.0);
    TH1D* h1d_pT           = new TH1D("h1d_pT","h1d_pT", 1000,0,200);
    TH1D* h1d_eta          = new TH1D("h1d_eta","h1d_eta", 200,-1.0,1.0);
    TH1D* h1d_phi          = new TH1D("h1d_phi","h1d_phi", 800,-4.0,4.0);

    TH2D* h1d_mult_vs_cent = new TH2D("h1d_mult_vs_cent","h1d_mult_vs_cent", 500,0,5000,80,0,40);

    //Generate output file
    TString str_out =  "./" + i_OutputFile;
    str_out += ".root";
    TFile* OutputFile = new TFile(str_out.Data(),"RECREATE");
    OutputFile->cd();
    cout << "Output file is " << str_out << endl;
    TDirectory *EventHistos = OutputFile->mkdir("Event Histos");
    TDirectory *Results = OutputFile->mkdir("Results");

    //Get data from tree file
    TFile* file = TFile::Open(InputFile);

    TTree * collisions = nullptr;
    TTree * tracks = nullptr;

    // loop over all directories and print name
    cout << "Process file" << endl;
    TIter next(file->GetListOfKeys());
    TKey *key;
    while ((key = (TKey *)next())) {
        TClass *cl = gROOT->GetClass(key->GetClassName());
        if (!cl->InheritsFrom("TDirectory"))
            continue;
        TDirectory *dir = (TDirectory *)key->ReadObj();
        TIter next2(dir->GetListOfKeys());
        TKey *key2;
        while ((key2 = (TKey *)next2())) {
            TClass *cl2 = gROOT->GetClass(key2->GetClassName());
            if (!cl2->InheritsFrom("TTree"))
                continue;
            TTree *tree = (TTree *)key2->ReadObj();
            if (strcmp(tree->GetName(), "O2tablecol") == 0) {
                collisions = tree;
            } else if (strcmp(tree->GetName(), "O2tabletrack") == 0) {
                tracks = tree;
            }
        }
        cout << "Processing directory " << dir->GetName() << endl;

        // Get basic information from the dataframe
        Int_t entries_col   =  collisions->GetEntries();
        Int_t entries_track =  tracks->GetEntries();
        cout << "------------------------ COLLISION TABLE -----------------------------------" << endl;
        cout <<  entries_col << " entries " << endl;
        cout << "------------------------ TRACK TABLE -----------------------------------" << endl;
        cout <<  entries_track << " entries " << endl;
        cout << "-----------------------------------------------------------" << endl;

        std::unordered_map<int, std::vector<int>> CollisionMap; //holds the Collision ID and underneath the track indices inside of the whole trackChain (0,1,2,...,2644) for example
        Int_t T_CollisionID;
        int E_RunNumber, E_Multiplicity, E_Occupancy;
        float E_VertX, E_VertY, E_VertZ, E_Centrality, E_OccupancyFT0;
        short E_psi2_Substitute, E_psi3_Substitute, T_Charge, T_DCAxy_Substitute, T_DCAz_Substitute;
        unsigned short T_dEdX_Substitute;
        ULong64_t T_P_Substitute;
        
        collisions->SetBranchAddress("fRn",&E_RunNumber);
        collisions->SetBranchAddress("fCent",&E_Centrality);
        collisions->SetBranchAddress("fMult",&E_Multiplicity);
        collisions->SetBranchAddress("fOccu",&E_Occupancy);
        collisions->SetBranchAddress("fOccuft0",&E_OccupancyFT0);
        collisions->SetBranchAddress("fVertexX",&E_VertX);
        collisions->SetBranchAddress("fVertexY",&E_VertY);
        collisions->SetBranchAddress("fVertexZ",&E_VertZ);
        collisions->SetBranchAddress("fPsi2",&E_psi2_Substitute);
        collisions->SetBranchAddress("fPsi3",&E_psi3_Substitute);

        tracks->SetBranchAddress("fIndexTableCols", &T_CollisionID);
        tracks->SetBranchAddress("fCharge", &T_Charge);
        tracks->SetBranchAddress("fP", &T_P_Substitute);
        tracks->SetBranchAddress("fDedx", &T_dEdX_Substitute);
        tracks->SetBranchAddress("fDcaxy", &T_DCAxy_Substitute);
        tracks->SetBranchAddress("fDcaz", &T_DCAz_Substitute);

        // Make an unordered map particle->collision
        bool scrawl = false;
        for(uint64_t i_Track = 0; i_Track < entries_track; i_Track ++)
        {
            tracks->GetEntry(i_Track);
            CollisionMap[T_CollisionID].push_back(i_Track);
        }
        cout << "Number of different collisions in the tracks: " << CollisionMap.size() << endl;

        for (int j = 0; j < entries_col; j++) {

            // Get event properties
            collisions->GetEntry(j);
            int E_CollisionID = j;
            E_RunNumber;
            E_Centrality;
            E_Multiplicity;
            E_Occupancy;
            E_OccupancyFT0;
            E_VertX;
            E_VertY;
            E_VertZ;
            float E_Psi2 = static_cast<float>(E_psi2_Substitute)/1000.0f;
            float E_Psi3 = static_cast<float>(E_psi3_Substitute)/1000.0f;
        
            
            TH2D* h2D_phi_vs_eta = new TH2D("h2D_phi_vs_eta","h2D_phi_vs_eta", DEF_BinningPerUnit * 2 * Pi, -Pi, Pi, DEF_BinningPerUnit * 2 * 0.9, -0.9, 0.9);

            // loop over jets for this event
            for (int k = 0; k < CollisionMap[E_CollisionID].size(); k++)
            {
                // Get track properties
                tracks->GetEntry(CollisionMap[E_CollisionID].at(k));
                T_Charge;
                float T_px = get_px_Particle(T_P_Substitute);
                float T_py = get_py_Particle(T_P_Substitute);
                float T_pz = get_pz_Particle(T_P_Substitute);
                float T_dEdX = static_cast<float>(T_dEdX_Substitute)/10.0f;
                float T_DCAxy = static_cast<float>(T_DCAxy_Substitute)/100.0f;
                float T_DCAz = static_cast<float>(T_DCAz_Substitute)/100.0f;
                float T_E = sqrt(pow(T_px, 2) + pow(T_py, 2) + pow(T_pz, 2) + pow(0.1349766, 2));
                ROOT::Math::PxPyPzEVector ParticleLV(T_px, T_py, T_pz, T_E);
                float T_Eta = ParticleLV.Eta();
                float T_Phi = ParticleLV.Phi();
                float T_Pt = sqrt(pow(T_px, 2) + pow(T_py, 2));

                h2D_phi_vs_eta->Fill(T_Phi, T_Eta, T_Pt);
                
            }

            EventHistos->cd();
            h2D_phi_vs_eta->Write();
            delete(h2D_phi_vs_eta);

        }

        OutputFile->Close();
        
        return 1;

    }
    return 1;
}