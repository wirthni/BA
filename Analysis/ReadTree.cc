#include "./ReadTree.h"

Int_t ReadTree(TString InputFile, TString i_OutputFile)
{
    cout << "ReadTree started" << endl;

    /*USER VARIABLES*/

    #define DEF_BinningPerUnit 100

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

    /*END USER VARIABLES*/

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
        cout << "------------------------ COLLISION TABLE -------------------------------" << endl;
        cout <<  entries_col << " entries " << endl;
        cout << "------------------------ TRACK TABLE -----------------------------------" << endl;
        cout <<  entries_track << " entries " << endl;
        cout << "------------------------------------------------------------------------" << endl;

        std::unordered_map<int, std::vector<int>> CollisionMap; //holds the Collision ID and underneath the track indices inside of the whole trackChain (0,1,2,...,2644) for example

        // Make an unordered map particle->collision
        for(uint64_t i_Track = 0; i_Track < entries_track; i_Track ++)
        {
            Particle track = Particle::Read(tracks, i_Track);
            CollisionMap[track.ColID].push_back(i_Track);
        }
        cout << "Have mapped the particles to the collisions.\nNumber of different collisions in the tracks: " << CollisionMap.size() << endl;

        for (int i_Collision = 0; i_Collision < entries_col; i_Collision++) {

            Collision collision = Collision::Read(collisions, i_Collision);
            if(CollisionMap[collision.ColID] == 0)
            {
                cout << "No tracks in collision " << collision.ColID << "!. Skipping." << endl;
                continue;
            }
            
            /*USER EVENT CODE*/
            h1d_mult_event->Fill(collision.Multiplicity);
            TH2D* h2D_phi_vs_eta = new TH2D("h2D_phi_vs_eta","h2D_phi_vs_eta", DEF_BinningPerUnit * 2 * Pi, -Pi, Pi, DEF_BinningPerUnit * 2 * 0.9, -0.9, 0.9);
            /*END USER EVENT CODE*/

            // loop over tracks for this event
            for (int i_Track = 0; i_Track < CollisionMap[collision.ColID].size(); i_Track++)
            {
                Particle track = Particle::Read(tracks, CollisionMap[collision.ColID].at(i_Track));

                /*USER TRACK CODE*/
                h2D_phi_vs_eta->Fill(track.Phi, track.Eta, track.Pt);
                /*END USER TRACK CODE*/
            }

            /*USER EVENT CODE*/
            EventHistos->cd();
            h2D_phi_vs_eta->Write();
            delete(h2D_phi_vs_eta);
            /*END USER EVENT CODE*/

        }

        /*USER FINALIZE DF CODE*/
        
        /*END USER FINALIZE DF CODE*/

    }

    /*USER FINISH CODE*/
    Results->cd();
    h1d_mult_event->Write();
    OutputFile->Close();
    /*END USER FINISH CODE*/

    return 1;
}