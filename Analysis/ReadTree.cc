#include "./ReadTree.h"

Int_t ReadTree(TString InputFile, TString i_OutputFile)
{
    cout << "ReadTree started" << endl;

    /*USER VARIABLES*/

    #define DEF_BinningPerUnit 100

    double pT_thresh = 20.0;   // threshold that was used to produce the trees
    int total_events = 0;
    int total_skip   = 0;

    // Make histograms for a global event analysis
    TH1D* h1d_mult_event   = new TH1D("h1d_mult_event","h1d_mult_event", 500,0,30000);
    TH1D* h1d_Ntrack_event = new TH1D("h1d_Ntrack_event","h1d_Ntrack_event", 2500,0,5000);
    TH1D* h1d_occ_event    = new TH1D("h1d_occ_event","h1d_occ_event", 1000,0,10000);
    TH1D* h1d_occft0_event = new TH1D("h1d_occft0_event","h1d_occft0_event", 100000,0,100000);
    TH1D* h1d_pTsum        = new TH1D("h1d_pTsum","h1d_pTsum", 2000,0,8000);
    TH1D* h1d_cent         = new TH1D("h1d_cent","h1d_cent", 15,0,15);
    TH1D* h1d_psi2         = new TH1D("h1d_psi2","h1d_psi2", 800,-4.0,4.0);
    TH1D* h1d_psi3         = new TH1D("h1d_psi3","h1d_psi3", 800,-4.0,4.0);
    TH2D* h1d_mult_vs_cent = new TH2D("h1d_mult_vs_cent","h1d_mult_vs_cent", 15,0,15,500,0,30000);

    // Track histos
    TH1D* h1d_pT           = new TH1D("h1d_pT","h1d_pT", 1000,0,200);
    TH1D* h1d_eta          = new TH1D("h1d_eta","h1d_eta", 200,-1.0,1.0);
    TH1D* h1d_phi          = new TH1D("h1d_phi","h1d_phi", 800,-4.0,4.0);


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
        cout << "+++++++++++++++++++++++++++++ Processing directory " << dir->GetName() << " +++++++++++++++++++++++++" << endl;

        // Get basic information from the dataframe
        Int_t entries_col   =  collisions->GetEntries();
        Int_t entries_track =  tracks->GetEntries();
        cout << "------------------- COLLISION TABLE --------------------------" << endl;
        cout <<  entries_col << " entries " << endl;
        cout << "------------------- TRACK TABLE ------------------------------" << endl;
        cout <<  entries_track << " entries " << endl;
        cout << "--------------------------------------------------------------" << endl;

        total_events += entries_col;

        // Definitions
        std::unordered_map<int, std::vector<int>> CollisionMap; //holds the Collision ID and underneath the track indices inside of the whole trackChain (0,1,2,...,2644) for example
        std::unordered_map<int, bool> HighpTMap;    // to check whether the collision has a high pT track
        std::unordered_set<int> finishedIDs;
        int currentID = -1;
        bool grouped  = true;
        int max_ID    = 0;


        // Make an unordered map particle -> collision
        for(uint64_t i_Track = 0; i_Track < entries_track; i_Track ++)
        {
            Particle track = Particle::Read(tracks, i_Track);
            //cout << "track " << i_Track << ", col ID " << track.ColID << endl;
            CollisionMap[track.ColID].push_back(i_Track);
            if(track.Pt > pT_thresh)HighpTMap[track.ColID] = true;

            if(track.ColID > max_ID) max_ID = track.ColID;

            //------ this checks whether tracks are grouped -------
            if(track.ColID != currentID)
            {
                if(currentID != -1) finishedIDs.insert(currentID);

                // Have we already seen and finished this ID before?
                if(finishedIDs.count(track.ColID))
                {
                    cout << "Collision ID " << track.ColID << " appears again at track "  << i_Track << endl;

                    grouped = false;
                    break;
                }

                currentID = track.ColID;
            }
            //---------------------------------------------------
        }
        cout << "Have mapped the particles to the collisions.\nNumber of different collisions in the tracks: " << CollisionMap.size() << endl;
        cout << "Grouped = " << (grouped ? "YES" : "NO") << endl;

        //------- More sanity checks -------------------------------------------------------
        if(max_ID+1 > entries_col){
            cout << "[ERROR] max ID is bigger than N_collision!" << endl;
            break;
        }
        if(CollisionMap.size() != entries_col){
            cout << "[ERROR] Number of indices doesnt match number of collisions! There should be " << entries_col << " but there are " << CollisionMap.size() << endl;
            break;
        }
        //----------------------- Checks done ---------------------------------------------

        cout << "Start event loop" << endl;
        int N_events_skipped = 0;

        for (int i_Collision = 0; i_Collision < entries_col; i_Collision++) {
        //for (int i_Collision = 46; i_Collision < 50; i_Collision++) {

            Collision collision = Collision::Read(collisions, i_Collision);
            if(CollisionMap[collision.ColID].size() == 0){
                cout << "No tracks in collision " << collision.ColID << "! Skipping." << endl;
                N_events_skipped++;
                continue;
            }
            if(!HighpTMap[collision.ColID]){
                //cout << "No high pT tracks in collision " << collision.ColID << "! Skipping." << endl;
                N_events_skipped++;
                continue;
            }
            // Fill event histos
            h1d_mult_event   ->Fill(collision.Multiplicity);
            h1d_occ_event    ->Fill(collision.Occupancy);
            h1d_occft0_event ->Fill(collision.OccupancyFT0);
            h1d_cent         ->Fill(collision.Centrality);
            h1d_psi2         ->Fill(collision.Psi2);
            h1d_psi3         ->Fill(collision.Psi3);
            h1d_mult_vs_cent ->Fill(collision.Centrality,collision.Multiplicity);
            h1d_Ntrack_event ->Fill(CollisionMap[collision.ColID].size());
            
            /*USER EVENT CODE*/
            //TH2D* h2D_phi_vs_eta = new TH2D("h2D_phi_vs_eta","h2D_phi_vs_eta", DEF_BinningPerUnit * 2 * Pi, -Pi, Pi, DEF_BinningPerUnit * 2 * 0.9, -0.9, 0.9);
            /*END USER EVENT CODE*/

            // loop over tracks for this event

            float pT_sum = 0.0;
            for (int i_Track = 0; i_Track < CollisionMap[collision.ColID].size(); i_Track++)
            {
                Particle track = Particle::Read(tracks, CollisionMap[collision.ColID].at(i_Track));
                h1d_pT -> Fill(track.Pt);
                h1d_eta-> Fill(track.Eta);
                h1d_phi-> Fill(track.Phi);
                pT_sum += track.Pt;
                /*USER TRACK CODE*/
                //h2D_phi_vs_eta->Fill(track.Phi, track.Eta, track.Pt);
                /*END USER TRACK CODE*/
            }

            h1d_pTsum ->Fill(pT_sum);


            /*USER EVENT CODE*/
            //EventHistos->cd();
            //h2D_phi_vs_eta->Write();
            //delete(h2D_phi_vs_eta);
            /*END USER EVENT CODE*/

        }

        cout << "Skipped " << N_events_skipped << "/" << entries_col << " events" << endl;
        total_skip += N_events_skipped;
        /*USER FINALIZE DF CODE*/

        /*END USER FINALIZE DF CODE*/

    }

    /*USER FINISH CODE*/
    cout << "Total number of events: " << total_events << ", of which skipped " << total_skip << endl;
    cout << "Write ..." << endl;
    EventHistos->cd();

    h1d_mult_event   ->Write();
    h1d_Ntrack_event ->Write();
    h1d_occ_event    ->Write();
    h1d_occft0_event ->Write();
    h1d_cent         ->Write();
    h1d_psi2         ->Write();
    h1d_psi3         ->Write();
    h1d_mult_vs_cent ->Write();
    h1d_pTsum        ->Write();

    h1d_pT           ->Write();
    h1d_eta          ->Write();
    h1d_phi          ->Write();

    OutputFile->Close();
    /*END USER FINISH CODE*/

    return 1;
}