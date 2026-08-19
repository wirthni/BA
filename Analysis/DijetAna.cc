#include "./ReadTree.cc"

using namespace fastjet;
using namespace std;

/* Parameters:
i_InputFile: Format xxx.root
*/
Int_t DijetAna(TString i_InputFile = "in.root")
{
    //CONFIGURE
    bool i_LocalDebugRun = true;
    bool i_UseHadronInstead = false;
    int i_NoOfEvents = -1;

    gStyle->SetOptStat(0);
    SetRootGraphicStyle();

    #define DEF_BinningPerUnit 100
    #define DEF_JetRadius 0.2
    #define DEF_OutputEventOverviews false
    #define DEF_HadLeadingPt 15.0
    #define DEF_HadSubleadingPt 10.0
    #define DEF_BackgroundLimit 7.0 //GeV
    #define DEF_MaxPartilclesPerJet 200
    #define DEF_MaxJetPt 150.0
    #define DEF_AxisLabelSize 0.05
    #define DEF_HistoTitleSize 0.1
    #define DEF_DijetMinPhiSeparation Pi/2.0

    int TotalEvents = 0;
    int ProcessedEvents = 0;
    int EventsSkipped   = 0;
    int RunNumber = -1;    

    Long64_t LargeGapEventCounter = 0;
    Long64_t SmallGapEventCounter = 0;

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
    TH2F* h2F_2DCorrelation_eta_vs_dphi_pT_0_1[2];//[0]: SmallGap, [1]: LargeGap
    TH2F* h2F_2DCorrelation_eta_vs_dphi_pT_1_2[2];
    TH2F* h2F_2DCorrelation_eta_vs_dphi_pT_2_4[2];
    TH2F* h2F_2DCorrelation_eta_vs_dphi_pT_4_6[2];
    TH1F* h1F_1DCorrelation_eta_pT_0_1[2];
    TH1F* h1F_1DCorrelation_eta_pT_1_2[2];
    TH1F* h1F_1DCorrelation_eta_pT_2_4[2];
    TH1F* h1F_1DCorrelation_eta_pT_4_6[2];
    
    h2F_2DCorrelation_eta_vs_dphi_pT_0_1[0] = new TH2F("h2F_2DCorrelation_eta_vs_dphi_pT_0_1_SmallGap", "h2F_2DCorrelation_eta_vs_dphi_pT_0_1_SmallGap", DEF_BinningPerUnit*1.8, -0.9, 0.9, DEF_BinningPerUnit*2*Pi, -Pi/2, 3*Pi/2);
    h2F_2DCorrelation_eta_vs_dphi_pT_1_2[0] = new TH2F("h2F_2DCorrelation_eta_vs_dphi_pT_1_2_SmallGap", "h2F_2DCorrelation_eta_vs_dphi_pT_1_2_SmallGap", DEF_BinningPerUnit*1.8, -0.9, 0.9, DEF_BinningPerUnit*2*Pi, -Pi/2, 3*Pi/2);
    h2F_2DCorrelation_eta_vs_dphi_pT_2_4[0] = new TH2F("h2F_2DCorrelation_eta_vs_dphi_pT_2_4_SmallGap", "h2F_2DCorrelation_eta_vs_dphi_pT_2_4_SmallGap", DEF_BinningPerUnit*1.8, -0.9, 0.9, DEF_BinningPerUnit*2*Pi, -Pi/2, 3*Pi/2);
    h2F_2DCorrelation_eta_vs_dphi_pT_4_6[0] = new TH2F("h2F_2DCorrelation_eta_vs_dphi_pT_4_6_SmallGap", "h2F_2DCorrelation_eta_vs_dphi_pT_4_6_SmallGap", DEF_BinningPerUnit*1.8, -0.9, 0.9, DEF_BinningPerUnit*2*Pi, -Pi/2, 3*Pi/2);
    h1F_1DCorrelation_eta_pT_0_1[0] = new TH1F("h1F_1DCorrelation_eta_pT_0_1_SmallGap", "h1F_1DCorrelation_eta_pT_0_1_SmallGap", DEF_BinningPerUnit*1.8, -0.9, 0.9);
    h1F_1DCorrelation_eta_pT_1_2[0] = new TH1F("h1F_1DCorrelation_eta_pT_1_2_SmallGap", "h1F_1DCorrelation_eta_pT_1_2_SmallGap", DEF_BinningPerUnit*1.8, -0.9, 0.9);
    h1F_1DCorrelation_eta_pT_2_4[0] = new TH1F("h1F_1DCorrelation_eta_pT_2_4_SmallGap", "h1F_1DCorrelation_eta_pT_2_4_SmallGap", DEF_BinningPerUnit*1.8, -0.9, 0.9);
    h1F_1DCorrelation_eta_pT_4_6[0] = new TH1F("h1F_1DCorrelation_eta_pT_4_6_SmallGap", "h1F_1DCorrelation_eta_pT_4_6_SmallGap", DEF_BinningPerUnit*1.8, -0.9, 0.9);
    
    h2F_2DCorrelation_eta_vs_dphi_pT_0_1[1] = new TH2F("h2F_2DCorrelation_eta_vs_dphi_pT_0_1_LargeGap", "h2F_2DCorrelation_eta_vs_dphi_pT_0_1_LargeGap", DEF_BinningPerUnit*1.8, -0.9, 0.9, DEF_BinningPerUnit*2*Pi, -Pi/2, 3*Pi/2);
    h2F_2DCorrelation_eta_vs_dphi_pT_1_2[1] = new TH2F("h2F_2DCorrelation_eta_vs_dphi_pT_1_2_LargeGap", "h2F_2DCorrelation_eta_vs_dphi_pT_1_2_LargeGap", DEF_BinningPerUnit*1.8, -0.9, 0.9, DEF_BinningPerUnit*2*Pi, -Pi/2, 3*Pi/2);
    h2F_2DCorrelation_eta_vs_dphi_pT_2_4[1] = new TH2F("h2F_2DCorrelation_eta_vs_dphi_pT_2_4_LargeGap", "h2F_2DCorrelation_eta_vs_dphi_pT_2_4_LargeGap", DEF_BinningPerUnit*1.8, -0.9, 0.9, DEF_BinningPerUnit*2*Pi, -Pi/2, 3*Pi/2);
    h2F_2DCorrelation_eta_vs_dphi_pT_4_6[1] = new TH2F("h2F_2DCorrelation_eta_vs_dphi_pT_4_6_LargeGap", "h2F_2DCorrelation_eta_vs_dphi_pT_4_6_LargeGap", DEF_BinningPerUnit*1.8, -0.9, 0.9, DEF_BinningPerUnit*2*Pi, -Pi/2, 3*Pi/2);
    h1F_1DCorrelation_eta_pT_0_1[1] = new TH1F("h1F_1DCorrelation_eta_pT_0_1_LargeGap", "h1F_1DCorrelation_eta_pT_0_1_LargeGap", DEF_BinningPerUnit*1.8, -0.9, 0.9);
    h1F_1DCorrelation_eta_pT_1_2[1] = new TH1F("h1F_1DCorrelation_eta_pT_1_2_LargeGap", "h1F_1DCorrelation_eta_pT_1_2_LargeGap", DEF_BinningPerUnit*1.8, -0.9, 0.9);
    h1F_1DCorrelation_eta_pT_2_4[1] = new TH1F("h1F_1DCorrelation_eta_pT_2_4_LargeGap", "h1F_1DCorrelation_eta_pT_2_4_LargeGap", DEF_BinningPerUnit*1.8, -0.9, 0.9);
    h1F_1DCorrelation_eta_pT_4_6[1] = new TH1F("h1F_1DCorrelation_eta_pT_4_6_LargeGap", "h1F_1DCorrelation_eta_pT_4_6_LargeGap", DEF_BinningPerUnit*1.8, -0.9, 0.9);
    

    int AnalyzedEventsCounter = 0;

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
        vector<PseudoJet> HighPtParticles[entries_col];
        vector<PseudoJet> Jets[entries_col];
        int LargeGapEventCounter = 0;
        int SmallGapEventCounter = 0;

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

        PrintInfo("Start Diffusion Wake Analysis");

        for(int iEvent = 0; iEvent < entries_col; iEvent++)
        {
            
            PrintProgress(iEvent, entries_col);

            Collision collision = Collision::Read(collisions, iEvent);
            if(CollisionMap[collision.ColID].size() == 0 || ParticleVector[collision.ColID].size() == 0){
                cout << "No tracks in collision " << collision.ColID << "! Skipping." << endl;
                N_events_skipped++;
                continue;
            }

            //event specific variables
            bool m_MirrorEvent = false;
            bool m_LargeGapEvent = false;

            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            /*
            /
            /                                                     
            /                                             CLUSTER JETS BY FINDING HIGH PT PARTICLES 
            /                                                     
            /                                                     
            /                                                     
            /
            /                                                    
            /
            /
            /
            *//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            if(HighPtParticles[collision.ColID].size() == 0) continue;

            do{

                Int_t HighestPtIndex = 0;

                //search for highest momentum and choose corresponding particle as a primary vertex
                for(Int_t i=0; i<HighPtParticles[collision.ColID].size(); i++)
                {
                    if(HighPtParticles[i][2] > HighPtParticles[HighestPtIndex][2])
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
                    double DeltaPhi = Jets[collision.ColID].back().phi_std() - HighPtParticles[collision.ColID].phi_std();

                    //normalize
                    if(DeltaPhi > Pi) DeltaPhi -= 2*Pi;
                    else if(DeltaPhi < -Pi) DeltaPhi += 2*Pi;

                    double DeltaR = TMath::Sqrt(DeltaEta*DeltaEta + DeltaPhi*DeltaPhi);

                    if(DeltaR <= DEF_JetRadius)
                    {
                        //pop candidate
                        HighPtParticles[collision.ColID].erase(HighPtParticles.begin()+i);
                        i = -1;//start again for safety
                    }
                }

            }while(HighPtParticles.size() > 0);

            //sort high pt hadrons by size to have same structure as FindJets return type
            Jets[collision.ColID] = sorted_by_pt(Jets[collision.ColID]);
            cout << "Found jets: " << Jets[collision.ColID].size() << endl;
            //make eta cut
            //Selector Fiducial_cut_selector = SelectorAbsEtaMax(0.9 - DEF_JetRadius); // Fiducial cut for jets
            //JetVector = Fiducial_cut_selector(JetVector);


            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            /*
            /
            /                                                     
            /                                                      SORT EVENTS INTO THE TWO CATEGORIES
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

            if(Jets[collision.ColID][0].pt() < DEF_HadLeadingPt) continue;

            if(Jets[collision.ColId][1].pt() < DEF_HadSubleadingPt) continue;

            if(abs(Jets[collision.ColID][0].eta()) > 0.9 - DEF_JetRadius) continue;

            if(abs(Jets[collision.ColID][1].eta()) > 0.9 - DEF_JetRadius) continue;

            float PhiSeparation = Jets[collision.ColID][0].phi_std() - Jets[collision.ColID][1].phi_std();
            if(PhiSeparation > Pi) PhiSeparation -= 2*Pi;
            else if(PhiSeparation < -Pi) PhiSeparation += 2*Pi;

            if(abs(PhiSeparation) < DEF_DijetMinPhiSeparation) continue;

            if(Jets[collision.ColID][0].eta() < 0) m_MirrorEvent = true;

            if(Jets[collision.ColID][0].eta() * Jets[collision.ColID][1].eta() < 0) m_LargeGapEvent = true;



            for(const auto& particle : ParticleVector[collision.ColID])
            {

                //float EtaDistance = particle.eta() - RecoilSiteEta;
                float PhiDistance = particle.phi_std() - RecoilSitePhi;
                if(PhiDistance > Pi) PhiDistance -= 2*Pi;
                else if(PhiDistance < -Pi) PhiDistance += 2*Pi;

                if(abs(PhiDistance) <= DEF_MaxCorrelationDeltaPhi)
                {
                    NoOfParticlesInRecoilAreaCounter++;
                    h2F_CorrelationRecoil_eta_vs_pT->Fill(particle.eta(), particle.pt());
                }
            }

            //Fill reference correlation histogram
            TH2F* h2F_CorrelationReference_eta_vs_pT = new TH2F("h2F_CorrelationReference_eta_vs_pT", "h2F_CorrelationReference_eta_vs_pT", 50, -0.9, 0.9, 30, 0, 3);

            for(const auto& particle : ParticleVector[EventProperties[iEvent].SuitableReferenceEvent])
            {
                //float EtaDistance = particle.eta() - RecoilSiteEta;
                float PhiDistance = particle.phi_std() - RecoilSitePhi;
                if(PhiDistance > Pi) PhiDistance -= 2*Pi;
                else if(PhiDistance < -Pi) PhiDistance += 2*Pi;

                if(abs(PhiDistance) <= DEF_MaxCorrelationDeltaPhi)
                {
                    NoOfParticlesInReferenceAreaCounter++;
                    h2F_CorrelationReference_eta_vs_pT->Fill(particle.eta(), particle.pt());
                }
            }
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /*
        /
        /                                                     
        /                                 print overall information on the past analysis  
        /                                                     
        /                                                     
        /                                                     
        /
        /                                                    
        /
        /
        /
        *//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


        cout << "OVERALL ANALYSIS STATS" << endl;
        cout << "Total number of events: " << TotalEvents << ", of which skipped " << EventsSkipped << endl;
        
        EventsSkipped += N_events_skipped;

        if(i_LocalDebugRun) break;
    }

    h2F_CorrelationDifference_eta_vs_pT->Scale(1./(float)AnalyzedEventsCounter);
    h2F_CorrelationDifference_eta_vs_pT->Scale(1./(h2F_CorrelationDifference_eta_vs_pT->GetXaxis()->GetBinWidth(1) * h2F_CorrelationDifference_eta_vs_pT->GetYaxis()->GetBinWidth(1)));

    //Generate output file
    TString RunNumberAsString;
    RunNumberAsString.Form("%d", RunNumber);
    TString str_out =  "./" + (TString)RunNumberAsString.Data() + "_SingleJetAna.root";

    TFile* OutputFile;
    if(gSystem->AccessPathName(str_out) || i_LocalDebugRun)
    {
        OutputFile = new TFile(str_out,"RECREATE");
        OutputFile->cd();
        cout << "Output file is " << str_out << endl;
        h1F_NoOfParticlesInRecoilArea->Write();
        h1F_NoOfParticlesInReferenceArea->Write();
        h1F_N_ev_minus_N_ref->Write();
        h1F_NoOfHighPtParticles->Write();
        h2F_CorrelationDifference_eta_vs_pT->Write();

        OutputFile->Close();
    }
    else
    {
        OutputFile = TFile::Open(str_out.Data());
        cout << "ADD TO EXISTING OUTPUT FILE" << endl;
        OutputFile->cd();

        TH1F* Summed_h1F_NoOfParticlesInRecoilArea = (TH1F*) OutputFile->Get("h1F_NoOfParticlesInRecoilArea");
        Summed_h1F_NoOfParticlesInRecoilArea->Add(h1F_NoOfParticlesInRecoilArea);

        TH1F* Summed_h1F_NoOfParticlesInReferenceArea = (TH1F*) OutputFile->Get("h1F_NoOfParticlesInReferenceArea");
        Summed_h1F_NoOfParticlesInReferenceArea->Add(h1F_NoOfParticlesInReferenceArea);

        TH1F* Summed_h1F_N_ev_minus_N_ref = (TH1F*) OutputFile->Get("h1F_N_ev_minus_N_ref");
        Summed_h1F_N_ev_minus_N_ref->Add(h1F_N_ev_minus_N_ref);

        TH1F* Summed_h1F_NoOfHighPtParticles = (TH1F*) OutputFile->Get("h1F_NoOfHighPtParticles");
        Summed_h1F_NoOfHighPtParticles->Add(h1F_NoOfHighPtParticles);

        TH1F* Summed_h2F_CorrelationDifference_eta_vs_pT = (TH1F*) OutputFile->Get("h2F_CorrelationDifference_eta_vs_pT");
        Summed_h2F_CorrelationDifference_eta_vs_pT->Add(h2F_CorrelationDifference_eta_vs_pT);


        OutputFile = new TFile(str_out.Data(),"RECREATE");
        Summed_h1F_NoOfParticlesInRecoilArea->Write();
        Summed_h1F_NoOfParticlesInReferenceArea->Write();
        Summed_h1F_N_ev_minus_N_ref->Write();
        Summed_h1F_NoOfHighPtParticles->Write();
        Summed_h2F_CorrelationDifference_eta_vs_pT->Write();

        OutputFile->Close();

    }

    return 1;
}