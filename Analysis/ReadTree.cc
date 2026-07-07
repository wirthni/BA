#include "./ReadTree.h"
#include <Math/Vector4D.h>
#include <unordered_map>
#include <vector>
#include "functions.h"
#include "StJetTrackEvent.h"
#include "StJetTrackEventLinkDef.h"

Int_t ReadTree(TString InputFile, TString OutputFile)
{
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

        //NOW WE HAVE BOTH TREES
        Int_t entries_col   =  collisions->GetEntries();
        Int_t entries_track =  tracks->GetEntries();


        cout << "------------------------ COLLISION TABLE -----------------------------------" << endl;
        cout <<  entries_col << " entries " << endl;
        cout << "------------------------ TRACK TABLE -----------------------------------" << endl;
        cout <<  entries_track << " entries " << endl;
        cout << "-----------------------------------------------------------" << endl;


        cout << "Read and fill Tree" << endl;

        std::unordered_map<int, std::vector<int>> CollisionMap; //holds the Collision ID and underneath the track indices inside of the whole trackChain (0,1,2,...,2644) for example

        Long64_t gi;
        Int_t track_CollisionID;
        int32_t trackid;
        int rn, mult, occu;
        float vertx, verty, vertz, cent, occuft0, tof;
        short psi2, psi3, charge, dcaxy, dcaz;
        unsigned short dedx,length;
        ULong64_t p;

        //collisions->SetBranchAddress("fGI",&gi);
        //collisions->SetBranchAddress("fCi", &colID_event);
        
        collisions->SetBranchAddress("fRn",&rn);
        collisions->SetBranchAddress("fCent",&cent);
        collisions->SetBranchAddress("fMult",&mult);
        collisions->SetBranchAddress("fOccu",&occu);
        collisions->SetBranchAddress("fOccuft0",&occuft0);
        collisions->SetBranchAddress("fVertexX",&vertx);
        collisions->SetBranchAddress("fVertexY",&verty);
        collisions->SetBranchAddress("fVertexZ",&vertz);
        collisions->SetBranchAddress("fPsi2",&psi2);
        collisions->SetBranchAddress("fPsi3",&psi3);

        tracks->SetBranchAddress("fIndexTableCols", &track_CollisionID);
        tracks->SetBranchAddress("fCharge", &charge);
        tracks->SetBranchAddress("fP", &p);
        tracks->SetBranchAddress("fDedx", &dedx);
        tracks->SetBranchAddress("fDcaxy", &dcaxy);
        tracks->SetBranchAddress("fDcaz", &dcaz); 

        // Make an unordered map particle->collision
        cout << "Mapping" << endl;
        bool scrawl = false;
        for(uint64_t i_Track = 0; i_Track < entries_track; i_Track ++)
        {
            tracks->GetEntry(i_Track);
            int collisionID = track_CollisionID;
            CollisionMap[collisionID].push_back(i_Track);

            if(collisionID == 0 && !scrawl) 
            {
                cout << "A particle with collision ID = 0 occurs in this DF!" << endl;
                scrawl = true;
            }

        }

        cout << "Number of different collisions in the tracks: " << CollisionMap.size() << endl;

        //LOOP
        
        cout << "collision loop" << endl;
        for (int j = 0; j < entries_col; j++) {
            collisions->GetEntry(j);
            int collisionID = j;

            if(j < 50)cout << "Collision " << j << " - Multiplicity / actual Number:" << mult << "/" << CollisionMap[collisionID].size() << endl;

            // loop over jets for this event
            for (int k = 0; k < CollisionMap[collisionID].size(); k++)
            {
                tracks->GetEntry(CollisionMap[collisionID].at(k));
                //if(j < 5 && k < 10)cout << "Track " << k << " has P " << p << endl;
            }

        }

        

        cout << endl << endl;

    }
    return 1;
}