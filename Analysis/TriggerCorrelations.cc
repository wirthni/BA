#include "./ReadTree.cc"

using namespace fastjet;
using namespace std;

vector<PseudoJet> FindJets(const vector<PseudoJet> vec_particles);

/* Parameters:
i_InputFile: Format xxx.root
*/
Int_t TriggerCorrelations(TString i_InputFile = "in.root")
{
    //CONFIGURE
    bool UseJetTrigger = false;

    gStyle->SetOptStat(0);
    SetRootGraphicStyle();

    #define DEF_BinningPerUnit 25
    #define DEF_JetRadius 0.2
    #define DEF_JetLeadingPt 20.0
    #define DEF_JetSubleadingPt 10.0
    #define DEF_HadLeadingPt 15.0
    #define DEF_HadSubleadingPt 8.0
    #define DEF_BackgroundLimit 7.0 //GeV
    #define DEF_MaxPartilclesPerJet 200
    #define DEF_MaxJetPt 150.0
    #define DEF_AxisLabelSize 0.05
    #define DEF_HistoTitleSize 0.1
    #define DEF_DijetMinPhiSeparation Pi/2.0
    #define DEF_1DCorrelationPhiTolerance Pi/2.0

    int TotalEvents = 0;
    int ProcessedEvents = 0;
    int EventsSkipped   = 0;
    int RunNumber = -1;
    int LargeGapEventCounter = 0;
    int SmallGapEventCounter = 0;

    /*END USER VARIABLES*/

    //Get data from tree file
    TFile* file = TFile::Open(i_InputFile);
    if(!file || file->IsZombie())
    {
        cout << "Input File could not be opened. Wrong path?" << endl;
        return 0;
    }

    TTree * collisions = nullptr;
    TTree * tracks = nullptr;

    //APPLICATION VARIABLES
    Int_t NBins[6] = {72,72,36,72,100,100};
    Double_t XMin[6] = {-Pi, -Pi, -0.9, -1.8, 10, 0};
    Double_t XMax[6] = {Pi, Pi, 0.9, 1.8, 100, 100};

    //6D phase space array containing the occurence of kinematic variables on indices 0...5
    /*
    0: phi_leading - psi_2              in range    -pi...pi        with NBins      72      with resolution     5 deg
    1: phi_leading - phi_subleading                 -pi...pi                        72                          5 deg
    2: eta_leading                                  -0.9...0.9                      36                          0.05
    3: eta_leading - eta_subleading                 -1.8...1.8                      72                          0.05
    4: p_T_leading                                  10...100                        100                         1
    5: p_T_leading - p_T_subleading                 0...100                         100                         1
    */
    THnSparseI* TriggerPhaseSpace = new THnSparseI("TriggerPhaseSpace", "TriggerPhaseSpace", 6, NBins, XMin, XMax);
    TH1F* h1F_Psi2 = new TH1F("h1F_Psi2", "h1F_Psi2", 100, -Pi/2, Pi/2);

    // loop over all directories and print name
    cout << "Process file" << endl;
    TIter next(file->GetListOfKeys());
    TKey *key;
    while ((key = (TKey *)next())) {
        TString DFName = key->GetName();
        //if(DFName != i_DFName) continue;
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

        TotalEvents += entries_col;

        // Definitions
        std::unordered_map<int, std::vector<int>> CollisionMap; //holds the Collision ID and underneath the track indices inside of the whole trackChain (0,1,2,...,2644) for example
        std::unordered_set<int> finishedIDs;
        int currentID = -1;
        bool grouped  = true;
        int max_ID    = 0;
        int N_events_skipped = 0;

        //Dijet ana variables
        vector<PseudoJet> ParticleVector[entries_col];
        vector<PseudoJet> HighPtParticles[entries_col];
        vector<PseudoJet> Jets[entries_col];

        // Make an unordered map particle -> collision
        for(uint64_t i_Track = 0; i_Track < entries_track; i_Track ++)
        {
            Particle track = Particle::Read(tracks, i_Track);
            CollisionMap[track.ColID].push_back(i_Track);
            //prepare particles for jetfinder
            ParticleVector[track.ColID].push_back(PseudoJet(track.px, track.py, track.pz, track.E));
            if(track.Pt > DEF_BackgroundLimit)
            {
                HighPtParticles[track.ColID].push_back(PseudoJet(track.px, track.py, track.pz, track.E));
            }
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

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /*
        /
        /                                                     
        /                                                           ANALYSIS LOOP 
        /                                                     
        /                                                     
        /                                                     
        /
        /                                                    
        /
        /
        /
        *//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        PrintInfo("Start Trigger Analysis");

        for(int iEvent = 0; iEvent < entries_col; iEvent++)
        {
            
            PrintProgress(iEvent, entries_col);

            Collision collision = Collision::Read(collisions, iEvent);
            if(CollisionMap[collision.ColID].size() == 0 || ParticleVector[collision.ColID].size() == 0){
                cout << "No tracks in collision " << collision.ColID << "! Skipping." << endl;
                N_events_skipped++;
                continue;
            }

            RunNumber = collision.RunNumber;

            //Q vector correction to get correct event plane
            TString RunNumberAsString;
            RunNumberAsString.Form("%d", RunNumber);
            //differ between the LHC passes
            TString QVectorCorrectionFileName;
            if(i_InputFile.Contains("LHC23"))
            {
                QVectorCorrectionFileName =  "./RUN3_Data/QVecCalibRoots/LHC23_PbPb_pass5/" + (TString)RunNumberAsString.Data() + "_Calibration.root";
            }
            else if(i_InputFile.Contains("LHC24"))
            {
                QVectorCorrectionFileName =  "./RUN3_Data/QVecCalibRoots/LHC24ar_pass3/" + (TString)RunNumberAsString.Data() + "_Calibration.root";
            }
            else if(i_InputFile.Contains("LHC25"))
            {
                QVectorCorrectionFileName =  "./RUN3_Data/QVecCalibRoots/LHC25_PbPb_pass1/" + (TString)RunNumberAsString.Data() + "_Calibration.root";
            }
            
            //check if the file exists
            if(gSystem->AccessPathName(QVectorCorrectionFileName))
            {
                cout << "Q vector correction file does not exist!" << endl;
                return 0;
            }
            float Psi2_Corrected = GetEventPlaneInfo(collision.ColID, tracks, CollisionMap, QVectorCorrectionFileName);

            //event specific variables
            bool m_MirrorEvent = false;
            bool m_LargeGapEvent = false;

            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            /*
            /
            /                                                     
            /                                             FIND HIGH PT HADRON TRIGGERS   
            /                                                     
            /                                                     
            /                                                     
            /
            /                                                    
            /
            /
            /
            *//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            if(UseJetTrigger)
            {
                Jets[collision.ColID] = FindJets(ParticleVector[collision.ColID]);
            }
            else
            {
                if(HighPtParticles[collision.ColID].size() == 0) continue;
                do{

                    Int_t HighestPtIndex = 0;

                    //search for highest momentum and choose corresponding particle as a primary vertex
                    for(Int_t i=0; i<HighPtParticles[collision.ColID].size(); i++)
                    {
                        if(HighPtParticles[collision.ColID][i].pt() > HighPtParticles[collision.ColID][HighestPtIndex].pt())
                        {
                            HighestPtIndex = i;
                        }
                    }

                    //push particle into jet vector 
                    Jets[collision.ColID].push_back(HighPtParticles[collision.ColID][HighestPtIndex]);

                    //pop jet particle
                    HighPtParticles[collision.ColID].erase(HighPtParticles[collision.ColID].begin() + HighestPtIndex);

                    //look for other high pt particles in a radius of DEF_JetRadius because they belong to the same jet
                    for(Int_t i=0; i<HighPtParticles[collision.ColID].size(); i++)
                    {

                        double DeltaEta = Jets[collision.ColID].back().eta() - HighPtParticles[collision.ColID][i].eta();
                        double DeltaPhi = Jets[collision.ColID].back().phi() - HighPtParticles[collision.ColID][i].phi();

                        //normalize
                        if(DeltaPhi > Pi) DeltaPhi -= 2*Pi;
                        else if(DeltaPhi < -Pi) DeltaPhi += 2*Pi;

                        double DeltaR = TMath::Sqrt(DeltaEta*DeltaEta + DeltaPhi*DeltaPhi);

                        if(DeltaR <= DEF_JetRadius)
                        {
                            //pop candidate
                            HighPtParticles[collision.ColID].erase(HighPtParticles[collision.ColID].begin()+i);
                            i = -1;//start again for safety
                        }
                    }

                }while(HighPtParticles[collision.ColID].size() > 0);
            }

            if(Jets[collision.ColID].size() < 1) continue; 

            //sort high pt hadrons by size to have same structure as FindJets return type
            Jets[collision.ColID] = sorted_by_pt(Jets[collision.ColID]);

            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            /*
            /
            /                                                     
            /                                                      SORT OUT SINGLE JET EVENTS
            /                                                     
            /                                                     
            /                                                     
            /
            /                                                    
            /
            /
            /
            *//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            if(Jets[collision.ColID].size() < 2) continue;

            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            /*
            /
            /                                                     
            /                                                      COLLECT DATA 
            /                                                     
            /                                                     
            /                                                     
            /
            /                                                    
            /
            /
            /
            *//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            float Kinematics[6];

            Kinematics[0] = Jets[collision.ColID][0].phi() - Psi2_Corrected;
            while(Kinematics[0] > Pi) Kinematics[0] -= 2*Pi;
            while(Kinematics[0] < -Pi) Kinematics[0] += 2*Pi;

            Kinematics[1] = Jets[collision.ColID][0].phi() - Jets[collision.ColID][1].phi();
            while(Kinematics[1] > Pi) Kinematics[1] -= 2*Pi;
            while(Kinematics[1] < -Pi) Kinematics[1] += 2*Pi;

            Kinematics[2] = Jets[collision.ColID][0].eta();

            Kinematics[3] = Jets[collision.ColID][0].eta() - Jets[collision.ColID][1].eta();

            Kinematics[4] = Jets[collision.ColID][0].pt();

            Kinematics[5] = Jets[collision.ColID][0].pt() - Jets[collision.ColID][1].pt();

            TriggerPhaseSpace->Fill(Kinematics[0], Kinematics[1], Kinematics[2], Kinematics[3], Kinematics[4], Kinematics[5]);

            h1F_Psi2->Fill(Psi2_Corrected);
        }
    
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /*
    /
    /                                                     
    /                                                               EDIT HISTOGRAMS
    /                                                     
    /                                                     
    /                                                     
    /
    /                                                    
    /
    /
    /
    *//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    TriggerPhaseSpace->GetAxis(0)->SetTitle("#phi_{Leading} - #psi_{2}");
    TriggerPhaseSpace->GetAxis(1)->SetTitle("#phi_{Leading} - #phi_{Subleading}");
    TriggerPhaseSpace->GetAxis(2)->SetTitle("#eta_{Leading}");
    TriggerPhaseSpace->GetAxis(3)->SetTitle("#eta_{Leading} - #eta_{Subleading}");
    TriggerPhaseSpace->GetAxis(4)->SetTitle("p_{T}^{Leading}");
    TriggerPhaseSpace->GetAxis(5)->SetTitle("p_{T}^{Leading} - p_{T}^{Subleading}");

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /*
    /
    /                                                     
    /                                                         WRITE OUTPUT FILE
    /                                                     
    /                                                     
    /                                                     
    /
    /                                                    
    /
    /
    /
    *//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //Generate output file
    TString RunNumberAsString;
    RunNumberAsString.Form("%d", RunNumber);
    TString str_out =  "./" + (TString)RunNumberAsString.Data() + "_TriggerCorrelations.root";

    TFile* OutputFile;

    if(gSystem->AccessPathName(str_out))
    {
        OutputFile = new TFile(str_out,"RECREATE");
        OutputFile->cd();
        cout << "Output file is " << str_out << endl;
        
        TriggerPhaseSpace->Write();
        h1F_Psi2->Write();

        OutputFile->Close();
    }
    else
    {
        OutputFile = TFile::Open(str_out.Data());
        cout << "ADD TO EXISTING OUTPUT FILE" << endl;
        OutputFile->cd();

        THnSparseI* TriggerPhaseSpaceNew = (THnSparseI*) OutputFile->Get("TriggerPhaseSpace");
        TriggerPhaseSpaceNew->Add(TriggerPhaseSpace);
        TriggerPhaseSpaceNew->Write();

        TH1F* h1F_Psi2New = (TH1F*) OutputFile->Get("h1F_Psi2");
        h1F_Psi2New->Add(h1F_Psi2);
        h1F_Psi2New->Write();

        OutputFile->Close();

    }

    return 1;
}

vector<PseudoJet> FindJets(const vector<PseudoJet> vec_particles)
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

    //DEFINITIONS
    #define jet_radius 0.2
    Double_t ghost_maxrap = 1.1;
    Double_t eta_acceptance = 0.9;
    double ptmin = 0.15;
    JetDefinition jet_def(antikt_algorithm, jet_radius);
    GhostedAreaSpec area_spec(ghost_maxrap);
    AreaDefinition area_def(active_area_explicit_ghosts,GhostedAreaSpec(ghost_maxrap,1,0.01));

    //CLUSTER JETS
    ClusterSequenceArea clust_seq_hard(vec_particles, jet_def, area_def);
    vector<PseudoJet> jets_all = sorted_by_pt(clust_seq_hard.inclusive_jets(ptmin));
    Selector Fiducial_cut_selector = SelectorAbsEtaMax(eta_acceptance - jet_radius);
    vector<PseudoJet> FiducialJets = Fiducial_cut_selector(jets_all);

    //ESTIMATE AND SUBTRACT BACKGROUND
    Int_t Rem_n_hardest = 2;
    double eBkg_R = 0.2;
    JetDefinition jet_def_bkgd(kt_algorithm, eBkg_R);
    AreaDefinition area_def_bkgd(active_area_explicit_ghosts,GhostedAreaSpec(ghost_maxrap,1,0.01));
    ClusterSequenceArea ClustSeqBg(vec_particles, jet_def_bkgd, area_def_bkgd);
    Selector selector = SelectorAbsEtaMax(0.9 - eBkg_R) * (!SelectorNHardest(Rem_n_hardest));
    JetMedianBackgroundEstimator bkgd_estimator(selector, jet_def_bkgd, area_def_bkgd);
    bkgd_estimator.set_cluster_sequence(ClustSeqBg);
    Subtractor subtractor(&bkgd_estimator);
    subtractor.set_use_rho_m(true);
    vector<PseudoJet> SubtractedJets = subtractor(FiducialJets);

    //APPEND INFO AND KICK OUT JETS WITH TOO SMALL AREA
    Double_t Jet_area_cut = 0.56*TMath::Pi()*TMath::Power(jet_radius,2);
    for(int iJet = SubtractedJets.size() - 1; iJet >= 0; iJet--)
    {
        if(SubtractedJets[iJet].area() < Jet_area_cut){SubtractedJets.erase(SubtractedJets.begin() + iJet);}
        //SubtractedJets[iJet].set_user_info(new MyUserInfo(SubtractedJets[iJet].constituents().size()));

    }

    return sorted_by_pt(SubtractedJets);

}