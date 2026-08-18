#include "./ReadTree.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
FUNCTION NAME:
ReadTree

INFO:
Reads in a RUN 3 AO2D.root file and writes general dataset information into an output file. 
If executed once, the output .root file contains the Q vector correction histograms

PARAMETERS:
i_InputFile: Path and filename to the RUN 3 AO2D file
i_OutputFile: Free-to-choose filename for the generated output file (.root)

OUTPUT:
[Run Number]_Calibration.root
*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
Int_t ReadTree(TString i_InputFile)
{
    ReadTreeInit();

    PrintInfo("ReadTree started");

    #define DEF_BinningPerUnit 100

    double pT_thresh = 20.0;   // threshold that was used to produce the trees
    int total_events = 0;
    int total_skip   = 0;
    int RunNumber = 0;

    // Make histograms for a global event analysis
    TH1D* h1d_mult_event        = new TH1D("h1d_mult_event","h1d_mult_event", 500,0,30000);
    TH1D* h1d_Ntrack_event      = new TH1D("h1d_Ntrack_event","h1d_Ntrack_event", 2500,0,5000);
    TH1D* h1d_occ_event         = new TH1D("h1d_occ_event","h1d_occ_event", 1000,0,10000);
    TH1D* h1d_occft0_event      = new TH1D("h1d_occft0_event","h1d_occft0_event", 100000,0,100000);
    TH1D* h1d_pTsum             = new TH1D("h1d_pTsum","h1d_pTsum", 2000,0,8000);
    TH1D* h1d_cent              = new TH1D("h1d_cent","h1d_cent", 15,0,15);
    TH1D* h1d_psi2              = new TH1D("h1d_psi2","h1d_psi2", 800,-4.0,4.0);
    TH1D* h1d_psi2_corrected    = new TH1D("h1d_psi2_corrected","h1d_psi2_corrected", 800,-4.0,4.0);
    TH1D* h1d_psi3              = new TH1D("h1d_psi3","h1d_psi3", 800,-4.0,4.0);
    TH2D* h1d_mult_vs_cent      = new TH2D("h1d_mult_vs_cent","h1d_mult_vs_cent", 15,0,15,500,0,30000);
    
    // Track histos
    TH1D* h1d_pT                = new TH1D("h1d_pT","h1d_pT", 1000,0,200);
    TH1D* h1d_eta               = new TH1D("h1d_eta","h1d_eta", 200,-1.0,1.0);
    TH1D* h1d_phi               = new TH1D("h1d_phi","h1d_phi", 800,-4.0,4.0);

    //Get data from tree file
    TFile* file = TFile::Open(i_InputFile);

    TTree * collisions = nullptr;
    TTree * tracks = nullptr;

    PrintInfo("Write Q Vector correction");
    TIter next_corr(file->GetListOfKeys());
    TKey *key;
    while ((key = (TKey *)next_corr())) {
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
        FillQVectorCorrectionHistos_QVectorCorrectionRetrieved = false;//for this DF, the information has not been retrieved yet
        GetEventPlaneInfo_QVectorCorrectionRetrieved = false;//for this DF, the information has not been retrieved yet


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


        for (int i_Collision = 0; i_Collision < entries_col; i_Collision++) {

            PrintProgress(i_Collision, entries_col);

            Collision collision = Collision::Read(collisions, i_Collision);
            if(CollisionMap[collision.ColID].size() == 0){
                cout << "No tracks in collision " << collision.ColID << "! Skipping." << endl;
                continue;
            }
            if(!HighpTMap[collision.ColID]){
                continue;
            }

            RunNumber = collision.RunNumber;
            FillQVectorCorrectionHistos(collision.ColID, tracks , CollisionMap);

        }

    }
    
    PrintInfo("Generate Correction File");
    TString RunNumberAsString;
    RunNumberAsString.Form("%d", RunNumber);

    TString QVectorCorrectionFileName =  "./" + (TString)(RunNumberAsString.Data()) + "_Calibration.root";
    //look if the file already exists from another AO2D file in the same run
    TFile* QVectorCorrectionFile = TFile::Open(QVectorCorrectionFileName);
    if(!QVectorCorrectionFile)
    {
        QVectorCorrectionFile = new TFile(QVectorCorrectionFileName.Data(),"RECREATE");
        QVectorCorrectionFile->cd();
        cout << "Q vector correction file is " << QVectorCorrectionFileName << endl;
        h2D_qx_qy_for_EP_neg_eta->Write();
        h2D_qx_qy_for_EP_pos_eta->Write();
        QVectorCorrectionFile->Close();
    }
    else
    {
        cout << "ADD TO EXISTING Q VECTOR CORRECTION FILE" << endl;
        QVectorCorrectionFile->cd();
        TH2D* Summed_h2D_qx_qy_for_EP_neg_eta = (TH2D*)QVectorCorrectionFile->Get("h2D_qx_qy_for_EP_neg_eta");
        Summed_h2D_qx_qy_for_EP_neg_eta->Add(h2D_qx_qy_for_EP_neg_eta);
        
        TH2D* Summed_h2D_qx_qy_for_EP_pos_eta = (TH2D*)QVectorCorrectionFile->Get("h2D_qx_qy_for_EP_pos_eta");
        Summed_h2D_qx_qy_for_EP_pos_eta->Add(h2D_qx_qy_for_EP_pos_eta);

        QVectorCorrectionFile = new TFile(QVectorCorrectionFileName.Data(),"RECREATE");
        Summed_h2D_qx_qy_for_EP_pos_eta->Write();
        Summed_h2D_qx_qy_for_EP_neg_eta->Write();
        QVectorCorrectionFile->Close();

    }

    PrintInfo("Loop Q vector calibrated collisions");
    TIter next_loop(file->GetListOfKeys());
    while ((key = (TKey *)next_loop())) {
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
        FillQVectorCorrectionHistos_QVectorCorrectionRetrieved = false;//for this DF, the information has not been retrieved yet
        GetEventPlaneInfo_QVectorCorrectionRetrieved = false;//for this DF, the information has not been retrieved yet


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

        int N_events_skipped = 0;

        for (int i_Collision = 0; i_Collision < entries_col; i_Collision++) {

            PrintProgress(i_Collision, entries_col);

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
            h1d_psi2_corrected->Fill(GetEventPlaneInfo(collision.ColID, tracks , CollisionMap, QVectorCorrectionFileName));

            // loop over tracks for this event
            float pT_sum = 0.0;
            for (int i_Track = 0; i_Track < CollisionMap[collision.ColID].size(); i_Track++)
            {
                Particle track = Particle::Read(tracks, CollisionMap[collision.ColID].at(i_Track));
                h1d_pT -> Fill(track.Pt);
                h1d_eta-> Fill(track.Eta);
                h1d_phi-> Fill(track.Phi);
                pT_sum += track.Pt;
            }

            h1d_pTsum ->Fill(pT_sum);

        }

        /*USER FINALIZE DF CODE*/
        cout << "Skipped " << N_events_skipped << "/" << entries_col << " events" << endl;
        total_skip += N_events_skipped;
        /*END USER FINALIZE DF CODE*/

    }

    /*USER FINISH CODE*/
    cout << "Total number of events: " << total_events << ", of which skipped " << total_skip << endl;

    PrintInfo("Generate Output File");

    //look if the file already exists from another AO2D file in the same run
    TString OutputFileName =  "./" + (TString)(RunNumberAsString.Data()) + "_Overview.root";
    TFile* OutputFile = TFile::Open(OutputFileName);
    if(!OutputFile)
    {
        OutputFile = new TFile(OutputFileName.Data(),"RECREATE");
        OutputFile->cd();
        cout << "Output file is " << OutputFile << endl;

        h1d_mult_event          ->Write();
        h1d_Ntrack_event        ->Write();
        h1d_occ_event           ->Write();
        h1d_occft0_event        ->Write();
        h1d_cent                ->Write();
        h1d_psi2                ->Write();
        h1d_psi2_corrected      ->Write();
        h1d_psi3                ->Write();
        h1d_mult_vs_cent        ->Write();
        h1d_pTsum               ->Write();
        h1d_pT                  ->Write();
        h1d_eta                 ->Write();
        h1d_phi                 ->Write();
        //global histograms
        vec_h_psi_test[0]->Write();
        vec_h_psi_test[1]->Write();
        h_EP_peak[0]->Write();
        h_EP_peak[1]->Write();
        h2D_QA_vs_QB_peak[0]->Write();
        h2D_QA_vs_QB_tail[0]->Write();
        h2D_QA_vs_QB_peak[1]->Write();
        h2D_QA_vs_QB_tail[1]->Write();
        h_Psi_SE_EP_resolution->Write();
        h2D_qx_qy_for_EP_pos_eta->Write();
        h2D_qx_qy_for_EP_neg_eta->Write();
        h2D_Psi_pos_vs_Psi_neg->Write();
        h2D_Q_vec_peak->Write();
        h2D_Q_vec_tail->Write();

        OutputFile->Close();
    }
    else
    {
        cout << "ADD TO EXISTING OUTPUT FILE" << endl;
        OutputFile->cd();
        
        TH1D* Summed_h1d_mult_event = (TH1D*)OutputFile->Get("h1d_mult_event");
        cout << "Before it had entries: " << Summed_h1d_mult_event->GetEntries() << endl;
        Summed_h1d_mult_event->Add(h1d_mult_event);
        cout << "Then it had entries: " << Summed_h1d_mult_event->GetEntries() << endl;

        TH1D* Summed_h1d_Ntrack_event = (TH1D*)OutputFile->Get("h1d_Ntrack_event");
        Summed_h1d_Ntrack_event->Add(h1d_Ntrack_event);
        

        TH1D* Summed_h1d_occ_event = (TH1D*)OutputFile->Get("h1d_occ_event");
        Summed_h1d_occ_event->Add(h1d_occ_event);
        

        TH1D* Summed_h1d_occft0_event = (TH1D*)OutputFile->Get("h1d_occft0_event");
        Summed_h1d_occft0_event->Add(h1d_occft0_event);
        

        TH1D* Summed_h1d_cent = (TH1D*)OutputFile->Get("h1d_cent");
        Summed_h1d_cent->Add(h1d_cent);
        

        TH1D* Summed_h1d_psi2 = (TH1D*)OutputFile->Get("h1d_psi2");
        Summed_h1d_psi2->Add(h1d_psi2);
        

        TH1D* Summed_h1d_psi2_corrected = (TH1D*)OutputFile->Get("h1d_psi2_corrected");
        Summed_h1d_psi2_corrected->Add(h1d_psi2_corrected);
        

        TH1D* Summed_h1d_psi3 = (TH1D*)OutputFile->Get("h1d_psi3");
        Summed_h1d_psi3->Add(h1d_psi3);
        

        TH2D* Summed_h1d_mult_vs_cent = (TH2D*)OutputFile->Get("h1d_mult_vs_cent");
        Summed_h1d_mult_vs_cent->Add(h1d_mult_vs_cent);
        

        TH1D* Summed_h1d_pTsum = (TH1D*)OutputFile->Get("h1d_pTsum");
        Summed_h1d_pTsum->Add(h1d_pTsum);
        

        TH1D* Summed_h1d_pT = (TH1D*)OutputFile->Get("h1d_pT");
        Summed_h1d_pT->Add(h1d_pT);
        

        TH1D* Summed_h1d_eta = (TH1D*)OutputFile->Get("h1d_eta");
        Summed_h1d_eta->Add(h1d_eta);
        

        TH1D* Summed_h1d_phi = (TH1D*)OutputFile->Get("h1d_phi");
        Summed_h1d_phi->Add(h1d_phi);
        

        //global histograms
        TH1D* Summed_h_psi_test_pos = (TH1D*)OutputFile->Get("h_psi_test_pos");
        Summed_h_psi_test_pos->Add(vec_h_psi_test[0]);
        

        TH1D* Summed_h_psi_test_neg = (TH1D*)OutputFile->Get("h_psi_test_neg");
        Summed_h_psi_test_neg->Add(vec_h_psi_test[1]);
        

        TH1D* Summed_h_EP_peak_eta_pos = (TH1D*)OutputFile->Get("h_EP_peak_eta_pos");
        Summed_h_EP_peak_eta_pos->Add(h_EP_peak[0]);
        

        TH1D* Summed_h_EP_peak_eta_neg = (TH1D*)OutputFile->Get("h_EP_peak_eta_neg");
        Summed_h_EP_peak_eta_neg->Add(h_EP_peak[1]);
        

        TH2D* Summed_h2D_QA_vs_QB_peak_eta_pos = (TH2D*)OutputFile->Get("h2D_QA_vs_QB_peak_eta_pos");
        Summed_h2D_QA_vs_QB_peak_eta_pos->Add(h2D_QA_vs_QB_peak[0]);
        

        TH2D* Summed_h2D_QA_vs_QB_tail_eta_pos = (TH2D*)OutputFile->Get("h2D_QA_vs_QB_tail_eta_pos");
        Summed_h2D_QA_vs_QB_tail_eta_pos->Add(h2D_QA_vs_QB_tail[0]);
        

        TH2D* Summed_h2D_QA_vs_QB_peak_eta_neg = (TH2D*)OutputFile->Get("h2D_QA_vs_QB_peak_eta_neg");
        Summed_h2D_QA_vs_QB_peak_eta_neg->Add(h2D_QA_vs_QB_peak[1]);
        

        TH2D* Summed_h2D_QA_vs_QB_tail_eta_neg = (TH2D*)OutputFile->Get("h2D_QA_vs_QB_tail_eta_neg");
        Summed_h2D_QA_vs_QB_tail_eta_neg->Add(h2D_QA_vs_QB_tail[1]);
        

        TH1D* Summed_h_Psi_SE_EP_resolution = (TH1D*)OutputFile->Get("h_Psi_SE_EP_resolution");
        Summed_h_Psi_SE_EP_resolution->Add(h_Psi_SE_EP_resolution);
        

        TH2D* Summed_h2D_qx_qy_for_EP_pos_eta = (TH2D*)OutputFile->Get("h2D_qx_qy_for_EP_pos_eta");
        Summed_h2D_qx_qy_for_EP_pos_eta->Add(h2D_qx_qy_for_EP_pos_eta);
        

        TH2D* Summed_h2D_qx_qy_for_EP_neg_eta = (TH2D*)OutputFile->Get("h2D_qx_qy_for_EP_neg_eta");
        Summed_h2D_qx_qy_for_EP_neg_eta->Add(h2D_qx_qy_for_EP_neg_eta);
        

        TH2D* Summed_h2D_Psi_pos_vs_Psi_neg = (TH2D*)OutputFile->Get("h2D_Psi_pos_vs_Psi_neg");
        Summed_h2D_Psi_pos_vs_Psi_neg->Add(h2D_Psi_pos_vs_Psi_neg);
        

        TH2D* Summed_h2D_Q_vec_peak = (TH2D*)OutputFile->Get("h2D_Q_vec_peak");
        Summed_h2D_Q_vec_peak->Add(h2D_Q_vec_peak);
        

        TH2D* Summed_h2D_Q_vec_tail = (TH2D*)OutputFile->Get("h2D_Q_vec_tail");
        Summed_h2D_Q_vec_tail->Add(h2D_Q_vec_tail);
        

        OutputFile = new TFile(OutputFileName.Data(),"RECREATE");

        Summed_h1d_mult_event->Write();
        Summed_h1d_Ntrack_event->Write();
        Summed_h1d_occ_event->Write();
        Summed_h1d_occft0_event->Write();
        Summed_h1d_cent->Write();
        Summed_h1d_psi2->Write();
        Summed_h1d_psi2_corrected->Write();
        Summed_h1d_psi3->Write();
        Summed_h1d_mult_vs_cent->Write();
        Summed_h1d_pTsum->Write();
        Summed_h1d_pT->Write();
        Summed_h1d_eta->Write();
        Summed_h1d_phi->Write();
        Summed_h_psi_test_pos->Write();
        Summed_h_psi_test_neg->Write();
        Summed_h_EP_peak_eta_pos->Write();
        Summed_h_EP_peak_eta_neg->Write();
        Summed_h2D_QA_vs_QB_peak_eta_pos->Write();
        Summed_h2D_QA_vs_QB_tail_eta_pos->Write();
        Summed_h2D_QA_vs_QB_peak_eta_neg->Write();
        Summed_h2D_QA_vs_QB_tail_eta_neg->Write();
        Summed_h_Psi_SE_EP_resolution->Write();
        Summed_h2D_qx_qy_for_EP_pos_eta->Write();
        Summed_h2D_qx_qy_for_EP_neg_eta->Write();
        Summed_h2D_Psi_pos_vs_Psi_neg->Write();
        Summed_h2D_Q_vec_peak->Write();
        Summed_h2D_Q_vec_tail->Write();

        OutputFile->Close();

    }
    

    return 1;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
FUNCTION NAME:
ReadTreeInit

INFO:
Is called from ReadTree() and necessary for variable initialization

PARAMETERS:
none
*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

Int_t ReadTreeInit()
{
    // Define histograms
    vec_h_psi_test[0]  = new TH1D("h_psi_test_pos","h_psi_test_pos",100,-180.0,180.0);
    vec_h_psi_test[1]  = new TH1D("h_psi_test_neg","h_psi_test_neg",100,-180.0,180.0);
    h_EP_peak[0] = new TH1D("h_EP_peak_eta_pos","h_EP_peak_eta_pos",180,-90,90);
    h_EP_peak[1] = new TH1D("h_EP_peak_eta_neg","h_EP_peak_eta_neg",180,-90,90);
    h2D_QA_vs_QB_peak[0] = new TH2D("h2D_QA_vs_QB_peak_eta_pos","h2D_QA_vs_QB_peak_eta_pos",1000,-0.5,0.5,1000,-0.5,0.5);
    h2D_QA_vs_QB_tail[0] = new TH2D("h2D_QA_vs_QB_tail_eta_pos","h2D_QA_vs_QB_tail_eta_pos",1000,-0.5,0.5,1000,-0.5,0.5);
    h2D_QA_vs_QB_peak[1] = new TH2D("h2D_QA_vs_QB_peak_eta_neg","h2D_QA_vs_QB_peak_eta_neg",1000,-0.5,0.5,1000,-0.5,0.5);
    h2D_QA_vs_QB_tail[1] = new TH2D("h2D_QA_vs_QB_tail_eta_neg","h2D_QA_vs_QB_tail_eta_neg",1000,-0.5,0.5,1000,-0.5,0.5);
    h_Psi_SE_EP_resolution = new TH1D("h_Psi_SE_EP_resolution","h_Psi_SE_EP_resolution",2000,-1,1);
    h2D_qx_qy_for_EP_pos_eta = new TH2D("h2D_qx_qy_for_EP_pos_eta","h2D_qx_qy_for_EP_pos_eta",1000,-0.5,0.5,1000,-0.5,0.5);
    h2D_qx_qy_for_EP_neg_eta = new TH2D("h2D_qx_qy_for_EP_neg_eta","h2D_qx_qy_for_EP_neg_eta",1000,-0.5,0.5,1000,-0.5,0.5);
    h2D_Psi_pos_vs_Psi_neg = new TH2D("h2D_Psi_pos_vs_Psi_neg","h2D_Psi_pos_vs_Psi_neg",180,-90,90,180,-90,90);
    h2D_Q_vec_peak = new TH2D("h2D_Q_vec_peak","h2D_Q_vec_peak",1000,-0.5,0.5,1000,-0.5,0.5);
    h2D_Q_vec_tail = new TH2D("h2D_Q_vec_tail","h2D_Q_vec_tail",1000,-0.5,0.5,1000,-0.5,0.5);

    //initialize variables
    
    GetEventPlaneInfo_QVectorCorrectionRetrieved = false;
    FillQVectorCorrectionHistos_QVectorCorrectionRetrieved = false;
    return 1;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
FUNCTION NAME:
GetEventPlaneInfo

INFO:
Is called from the application. Returns the Q-vector corrected event plane for an event

PARAMETERS:
i_CollisionID: Collision of the to-investigate event
i_TrackTree: Pointer to the tree of all tracks in the dataset
i_CollisionMap: Hands over the collision map from ReadTrees()
i_QVectorCorrectionFile: File and path to the Q vector correction file (generated by MakeQVectorCorrectionFile)
*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
float GetEventPlaneInfo(int i_CollisionID, TTree* i_TrackTree, std::unordered_map<int, std::vector<int>> i_CollisionMap, TString i_QVectorCorrectionFile)
{
    #define DEF_DCA_Cutoff 10.0

    TVector3 Qvec_eta_pos, Qvec_eta_neg;
    vector<TVector3>vec_TV3_Qvec_eta[2];
    vec_TV3_Qvec_eta[0].clear();
    vec_TV3_Qvec_eta[1].clear();

    for(int i_Track = 0; i_Track < i_CollisionMap[i_CollisionID].size(); i_Track++)
    {
        Particle track = Particle::Read(i_TrackTree, i_CollisionMap[i_CollisionID].at(i_Track));
        Float_t dca_xy = track.DCAxy;
        Float_t dca_z = track.DCAz;
        int charge = track.Charge;
        Float_t px_track = track.px;
        Float_t py_track = track.py;
        Float_t pz_track = track.pz;
        Float_t pt_track = track.Pt;
        Float_t E_track = track.E;
        Float_t eta_track = track.Eta;

        TLorentzVector TLV_Particle_prim;
        TLV_Particle_prim.SetPxPyPzE(px_track,py_track,pz_track,E_track);

        if(fabs(dca_xy) > DEF_DCA_Cutoff || fabs(dca_z) > DEF_DCA_Cutoff) continue;
        if((fabs(eta_track)) > 0.9) continue;

        if(pt_track < 3.0 && fabs(dca_xy) < 1.0 && fabs(dca_z) < 1.0)
        {
            if(eta_track < 0.9 && eta_track > 0.1)
            {
                Qvec_eta_pos.SetXYZ(px_track,py_track,0.0);
                vec_TV3_Qvec_eta[0].push_back(Qvec_eta_pos);
            }
            if(eta_track > -0.9 && eta_track < -0.1)
            {
                Qvec_eta_neg.SetXYZ(px_track,py_track,0.0);
                vec_TV3_Qvec_eta[1].push_back(Qvec_eta_neg);
            }
        }

    }

    float Psi_full_r  = 0.0;
    TVector3 TV3_sum_Qvec_eta[2];
    Double_t harmonic = 2.0;
    Double_t Psi_pos_neg[2];
    Double_t Psi_pos_neg_eta[2];
    Psi_pos_neg[0] = 0.0;
    Psi_pos_neg[1] = 0.0;
    Psi_pos_neg_eta[0] = 0.0;
    Psi_pos_neg_eta[1] = 0.0;
    Double_t QA_vec[2]= {0.0};
    Double_t QB_vec[2]= {0.0};
    Double_t Q_vec[2]= {0.0};

    if(!GetEventPlaneInfo_QVectorCorrectionRetrieved)
    {
        //open Q vector correction file
        TFile *QVectorCorrectionFile = TFile::Open(i_QVectorCorrectionFile);

        if (!QVectorCorrectionFile || QVectorCorrectionFile->IsZombie()) {
            PrintInfo("Error while opening the Q vector correction file! Program will be halted!");
            while(1);
        }

        h2D_qx_qy_for_EP_pos_eta = (TH2D*)QVectorCorrectionFile->Get("h2D_qx_qy_for_EP_pos_eta");
        if(!h2D_qx_qy_for_EP_pos_eta)
        {
            PrintInfo("Error while opening the Q vector correction histogram! Program will be halted!");
            while(1);
        }
        h2D_qx_qy_for_EP_neg_eta = (TH2D*)QVectorCorrectionFile->Get("h2D_qx_qy_for_EP_neg_eta");
        if(!h2D_qx_qy_for_EP_neg_eta)
        {
            PrintInfo("Error while opening the Q vector correction histogram! Program will be halted!");
            while(1);
        }

        Qvec_correction_qx[0] = h2D_qx_qy_for_EP_pos_eta -> GetMean(1);
        Qvec_correction_qy[0] = h2D_qx_qy_for_EP_pos_eta -> GetMean(2);
        Qvec_correction_qx[1] = h2D_qx_qy_for_EP_neg_eta -> GetMean(1);
        Qvec_correction_qy[1] = h2D_qx_qy_for_EP_neg_eta -> GetMean(2);

        GetEventPlaneInfo_QVectorCorrectionRetrieved = true;

    }

    for(Int_t i_eta_pos_neg = 0; i_eta_pos_neg <= 1; i_eta_pos_neg++)
    {
        TV3_sum_Qvec_eta[i_eta_pos_neg].SetXYZ(0.0,0.0,0.0);
        Double_t Psi_num = 0.0;
        Double_t Psi_den = 0.0;
        Double_t weight  = 1.0;
        Double_t track_phi = 0.0;
        Double_t weight_inv_phi = 1.0;
        Double_t bin_content  = 0.0;
        Int_t    N_tracks = 0;


        for(Int_t i_Qvec = 0; i_Qvec < (Int_t)vec_TV3_Qvec_eta[i_eta_pos_neg].size(); i_Qvec++)
        {
            TV3_sum_Qvec_eta[i_eta_pos_neg] += vec_TV3_Qvec_eta[i_eta_pos_neg][i_Qvec];
            weight    = vec_TV3_Qvec_eta[i_eta_pos_neg][i_Qvec].Perp(); // pt weight
            track_phi = vec_TV3_Qvec_eta[i_eta_pos_neg][i_Qvec].Phi(); // track azimuthal angle

            Psi_den += ((weight_inv_phi*weight*TMath::Sin(harmonic*track_phi)) - Qvec_correction_qy[i_eta_pos_neg]);
            Psi_num += ((weight_inv_phi*weight*TMath::Cos(harmonic*track_phi)) - Qvec_correction_qx[i_eta_pos_neg]);
            N_tracks++;
        }

        Psi_pos_neg[i_eta_pos_neg] = TMath::RadToDeg()*TMath::ATan2(Psi_num,Psi_den)/harmonic;
    }

    if(Psi_pos_neg[0] - Psi_pos_neg[1] > 90.0) Psi_pos_neg[1]  -= 180.0;
    else
    {
        if(Psi_pos_neg[1] - Psi_pos_neg[0] > 90.0) Psi_pos_neg[0]  -= 180.0;
    }

    Psi_full_r = (Psi_pos_neg[0] + Psi_pos_neg[1])/2.0;

    if(Psi_full_r > 90.0)  Psi_full_r -= 180.0;
    else if(Psi_full_r < -90.0) Psi_full_r += 180.0;
    
    //return average event plane
    return 2*Pi*Psi_full_r/360.0;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
FUNCTION NAME:
FillQVectorCorrectionHistos

INFO:
Is called from ReadTree(), calculates the q(Q) vector to an event and writes the Q vector corrected event data
into overview histograms

PARAMETERS:
i_CollisionID: Collision of the to-investigate event
i_TrackTree: Pointer to the tree of all tracks in the dataset
i_CollisionMap: Hands over the collision map from ReadTrees()
i_MakeQVectorCorrection: if true, the Q vector correction is made. If false, the resulting EP info is Q vector corrected
OUTPUT:
The function returns the Q vector corrected event plane angle for the current event (currently only psi_2)
*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
float FillQVectorCorrectionHistos(int i_CollisionID, TTree* i_TrackTree, std::unordered_map<int, std::vector<int>> i_CollisionMap)
{
    #define DEF_DCA_Cutoff 10.0

    TVector3 Qvec_eta_pos, Qvec_eta_neg;
    vector<TVector3>vec_TV3_Qvec_eta[2];
    vec_TV3_Qvec_eta[0].clear();
    vec_TV3_Qvec_eta[1].clear();

    for(int i_Track = 0; i_Track < i_CollisionMap[i_CollisionID].size(); i_Track++)
    {
        Particle track = Particle::Read(i_TrackTree, i_CollisionMap[i_CollisionID].at(i_Track));
        Float_t dca_xy = track.DCAxy;
        Float_t dca_z = track.DCAz;
        int charge = track.Charge;
        Float_t px_track = track.px;
        Float_t py_track = track.py;
        Float_t pz_track = track.pz;
        Float_t pt_track = track.Pt;
        Float_t E_track = track.E;
        Float_t eta_track = track.Eta;

        TLorentzVector TLV_Particle_prim;
        TLV_Particle_prim.SetPxPyPzE(px_track,py_track,pz_track,E_track);

        if(fabs(dca_xy) > DEF_DCA_Cutoff || fabs(dca_z) > DEF_DCA_Cutoff) continue;
        if((fabs(eta_track)) > 0.9) continue;

        if(pt_track < 3.0 && fabs(dca_xy) < 1.0 && fabs(dca_z) < 1.0)
        {
            if(eta_track < 0.9 && eta_track > 0.1)
            {
                Qvec_eta_pos.SetXYZ(px_track,py_track,0.0);
                vec_TV3_Qvec_eta[0].push_back(Qvec_eta_pos);
            }
            if(eta_track > -0.9 && eta_track < -0.1)
            {
                Qvec_eta_neg.SetXYZ(px_track,py_track,0.0);
                vec_TV3_Qvec_eta[1].push_back(Qvec_eta_neg);
            }
        }

    }

    float Psi_full_r  = 0.0;
    TVector3 TV3_sum_Qvec_eta[2];
    Double_t harmonic = 2.0;
    Double_t Psi_pos_neg[2];
    Double_t Psi_pos_neg_eta[2];
    Psi_pos_neg[0] = 0.0;
    Psi_pos_neg[1] = 0.0;
    Psi_pos_neg_eta[0] = 0.0;
    Psi_pos_neg_eta[1] = 0.0;
    Double_t QA_vec[2]= {0.0};
    Double_t QB_vec[2]= {0.0};
    Double_t Q_vec[2]= {0.0};
    Qvec_correction_qx[0] = 0;
    Qvec_correction_qx[1] = 0;
    Qvec_correction_qy[0] = 0;
    Qvec_correction_qy[1] = 0;

    for(Int_t i_eta_pos_neg = 0; i_eta_pos_neg <= 1; i_eta_pos_neg++)
    {
        TV3_sum_Qvec_eta[i_eta_pos_neg].SetXYZ(0.0,0.0,0.0);
        Double_t Psi_num = 0.0;
        Double_t Psi_den = 0.0;
        Double_t weight  = 1.0;
        Double_t track_phi = 0.0;
        Double_t weight_inv_phi = 1.0;
        Double_t bin_content  = 0.0;
        Int_t    N_tracks = 0;


        for(Int_t i_Qvec = 0; i_Qvec < (Int_t)vec_TV3_Qvec_eta[i_eta_pos_neg].size(); i_Qvec++)
        {
            TV3_sum_Qvec_eta[i_eta_pos_neg] += vec_TV3_Qvec_eta[i_eta_pos_neg][i_Qvec];
            weight    = vec_TV3_Qvec_eta[i_eta_pos_neg][i_Qvec].Perp(); // pt weight
            track_phi = vec_TV3_Qvec_eta[i_eta_pos_neg][i_Qvec].Phi(); // track azimuthal angle

            Psi_den += ((weight_inv_phi*weight*TMath::Sin(harmonic*track_phi)) - Qvec_correction_qy[i_eta_pos_neg]);
            Psi_num += ((weight_inv_phi*weight*TMath::Cos(harmonic*track_phi)) - Qvec_correction_qx[i_eta_pos_neg]);
            N_tracks++;
        }

        if(i_eta_pos_neg == 0)
        {
            QA_vec[0] = Psi_num/N_tracks;
            QB_vec[0] = Psi_den/N_tracks;
            Q_vec[0]  = TMath::Sqrt((QA_vec[0]*QA_vec[0])+(QB_vec[0]*QB_vec[0]));
            h2D_qx_qy_for_EP_pos_eta   ->Fill(Psi_num/N_tracks,Psi_den/N_tracks);
        }
        else if(i_eta_pos_neg == 1)
        {
            QA_vec[1] = Psi_num/N_tracks;
            QB_vec[1] = Psi_den/N_tracks;
            Q_vec[1]  = TMath::Sqrt((QA_vec[1]*QA_vec[1])+(QB_vec[1]*QB_vec[1]));
            h2D_qx_qy_for_EP_neg_eta   ->Fill(Psi_num/N_tracks,Psi_den/N_tracks);
        }
        
        Psi_pos_neg[i_eta_pos_neg] = TMath::RadToDeg()*TMath::ATan2(Psi_num,Psi_den)/harmonic;
    }


    h2D_Psi_pos_vs_Psi_neg ->Fill(Psi_pos_neg[1],Psi_pos_neg[0]);
    vec_h_psi_test[0]->Fill(Psi_pos_neg[0]);
    vec_h_psi_test[1]->Fill(Psi_pos_neg[1]);
    h_Psi_SE_EP_resolution ->Fill(TMath::Cos(2*(Psi_pos_neg[1]-Psi_pos_neg[0])));


    Psi_pos_neg_eta[0] = Psi_pos_neg[0];
    Psi_pos_neg_eta[1] = Psi_pos_neg[1];
    Double_t Psi_diff =  Psi_pos_neg[0] - Psi_pos_neg[1];
    if(Psi_diff < -90.0) Psi_diff += 180.0;
    if(Psi_diff > +90.0) Psi_diff -= 180.0;
    if(fabs(Psi_diff < 20))
    {
        h_EP_peak[0]->Fill(Psi_pos_neg[0]);
        h_EP_peak[1]->Fill(Psi_pos_neg[1]);
    }


    for(Int_t i_eta = 0 ; i_eta < 2; i_eta++)
    {
        if(fabs(Psi_diff < 20))  h2D_QA_vs_QB_peak[i_eta]->Fill(QA_vec[i_eta], QB_vec[i_eta]);
        if(fabs(Psi_diff > 60))  h2D_QA_vs_QB_tail[i_eta]->Fill(QA_vec[i_eta], QB_vec[i_eta]);
    }
    
    if(fabs(Psi_diff < 20)) h2D_Q_vec_peak->Fill(Q_vec[0], Q_vec[1]);
    if(fabs(Psi_diff > 60)) h2D_Q_vec_tail->Fill(Q_vec[0], Q_vec[1]);


    if(Psi_pos_neg[0] - Psi_pos_neg[1] > 90.0) Psi_pos_neg[1]  -= 180.0;
    else
    {
        if(Psi_pos_neg[1] - Psi_pos_neg[0] > 90.0) Psi_pos_neg[0]  -= 180.0;
    }

    Psi_full_r = (Psi_pos_neg[0] + Psi_pos_neg[1])/2.0;

    if(Psi_full_r > 90.0)  Psi_full_r -= 180.0;
    else if(Psi_full_r < -90.0) Psi_full_r += 180.0;

    //return average event plane
    return 2*Pi*Psi_full_r/360.0;
}

void PrintProgress(int i_Iteration, int i_Total)
{

    int barWidth = 100;
    int ItemsPerLine = i_Total/barWidth;

    if( ItemsPerLine != 0 && i_Iteration % ItemsPerLine == 0)
    {
        std::cout << "[";
        int pos = i_Iteration/ItemsPerLine;
        for (int i = 0; i < barWidth; ++i) {
            if (i < pos) std::cout << "=";
            else if (i == pos) std::cout << ">";
            else std::cout << " ";
        }
        cout << "] " << int((float)i_Iteration/(float)i_Total * 100.0) << " %\r";
        std::cout.flush();

        std::cout << std::endl;
    }

}

void PrintInfo(TString i_String)
{
    cout << endl;
    cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
    cout << i_String << endl;
    cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
}