#include "./ReadTree.cc"

using namespace fastjet;
using namespace std;

vector<PseudoJet> FindJets(const vector<PseudoJet> vec_particles);

/* Parameters:
i_InputFile: Format xxx.root
*/
Int_t DijetAna(TString i_InputFile = "in.root")
{
    //CONFIGURE
    bool UseJetTrigger = true;

    gStyle->SetOptStat(0);
    SetRootGraphicStyle();

    #define DEF_BinningPerUnit 100
    #define DEF_JetRadius 0.3
    #define DEF_JetLeadingPt 30.0
    #define DEF_JetSubleadingPt 15.0
    #define DEF_HadLeadingPt 25.0
    #define DEF_HadSubleadingPt 15.0
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
    TH1I* CounterHisto = new TH1I("CounterHisto", "CounterHisto", 1, 0, 1);

    TH2F* h2F_2DCorrelation_eta_vs_dphi_pT_0_1[2];//[0]: SmallGap, [1]: LargeGap
    TH2F* h2F_2DCorrelation_eta_vs_dphi_pT_1_2[2];
    TH2F* h2F_2DCorrelation_eta_vs_dphi_pT_2_4[2];
    TH2F* h2F_2DCorrelation_eta_vs_dphi_pT_4_6[2];
    TH1F* h1F_1DCorrelation_eta_pT_0_1[2];
    TH1F* h1F_1DCorrelation_eta_pT_1_2[2];
    TH1F* h1F_1DCorrelation_eta_pT_2_4[2];
    TH1F* h1F_1DCorrelation_eta_pT_4_6[2];
    TH1F* h1F_1DCorrelationDifference_eta_pT_0_1;
    TH1F* h1F_1DCorrelationDifference_eta_pT_1_2;
    TH1F* h1F_1DCorrelationDifference_eta_pT_2_4;
    TH1F* h1F_1DCorrelationDifference_eta_pT_4_6;
    
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
    
    h1F_1DCorrelationDifference_eta_pT_0_1 = new TH1F("h1F_1DCorrelationDifference_eta_pT_0_1", "h1F_1DCorrelationDifference_eta_pT_0_1", DEF_BinningPerUnit*1.8, -0.9, 0.9);
    h1F_1DCorrelationDifference_eta_pT_1_2 = new TH1F("h1F_1DCorrelationDifference_eta_pT_1_2", "h1F_1DCorrelationDifference_eta_pT_1_2", DEF_BinningPerUnit*1.8, -0.9, 0.9);
    h1F_1DCorrelationDifference_eta_pT_2_4 = new TH1F("h1F_1DCorrelationDifference_eta_pT_2_4", "h1F_1DCorrelationDifference_eta_pT_2_4", DEF_BinningPerUnit*1.8, -0.9, 0.9);
    h1F_1DCorrelationDifference_eta_pT_4_6 = new TH1F("h1F_1DCorrelationDifference_eta_pT_4_6", "h1F_1DCorrelationDifference_eta_pT_4_6", DEF_BinningPerUnit*1.8, -0.9, 0.9);

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

            RunNumber = collision.RunNumber;

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

            if(!UseJetTrigger && Jets[collision.ColID][0].pt() < DEF_HadLeadingPt) continue;

            if(UseJetTrigger && Jets[collision.ColID][0].pt() < DEF_JetLeadingPt) continue;

            if(!UseJetTrigger && Jets[collision.ColID][1].pt() < DEF_HadSubleadingPt) continue;

            if(UseJetTrigger && Jets[collision.ColID][1].pt() < DEF_JetSubleadingPt) continue;

            if(abs(Jets[collision.ColID][0].eta()) > 0.9 - DEF_JetRadius) continue;

            if(abs(Jets[collision.ColID][1].eta()) > 0.9 - DEF_JetRadius) continue;

            float PhiSeparation = Jets[collision.ColID][0].phi() - Jets[collision.ColID][1].phi();
            if(PhiSeparation > Pi) PhiSeparation -= 2*Pi;
            else if(PhiSeparation < -Pi) PhiSeparation += 2*Pi;

            if(abs(PhiSeparation) < DEF_DijetMinPhiSeparation) continue;

            if(Jets[collision.ColID][0].eta() < 0) m_MirrorEvent = true;

            if(Jets[collision.ColID][0].eta() * Jets[collision.ColID][1].eta() < 0) m_LargeGapEvent = true;

            (m_LargeGapEvent) ? (LargeGapEventCounter++) : (SmallGapEventCounter++);

            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            /*
            /
            /                                                     
            /                                                      FILL CORRELATION HISTOGRAMS
            /                                                     
            /                                                     
            /                                                     
            /
            /                                                    
            /
            /
            /
            *//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            for(const auto& particle : ParticleVector[collision.ColID])
            {
                float Sgn = 1;
                (m_MirrorEvent) ? (Sgn = -1) : (Sgn = 1);

                float PhiDistance = particle.phi() - Jets[collision.ColID][0].phi();
                if(PhiDistance > 3*Pi/2) PhiDistance -= 2*Pi;
                else if(PhiDistance < -Pi/2) PhiDistance += 2*Pi;

                if(particle.pt() < 1.0)
                {
                    h2F_2DCorrelation_eta_vs_dphi_pT_0_1[m_LargeGapEvent]->Fill(Sgn*particle.eta(), PhiDistance);
                    if(abs(PhiDistance) <= DEF_1DCorrelationPhiTolerance) h1F_1DCorrelation_eta_pT_0_1[m_LargeGapEvent]->Fill(Sgn*particle.eta());
                }
                else if(particle.pt() < 2.0 && particle.pt() >= 1.0)
                {
                    h2F_2DCorrelation_eta_vs_dphi_pT_1_2[m_LargeGapEvent]->Fill(Sgn*particle.eta(), PhiDistance);
                    if(abs(PhiDistance) <= DEF_1DCorrelationPhiTolerance) h1F_1DCorrelation_eta_pT_1_2[m_LargeGapEvent]->Fill(Sgn*particle.eta());
                }
                else if(particle.pt() < 4.0 && particle.pt() >= 2.0)
                {
                    h2F_2DCorrelation_eta_vs_dphi_pT_2_4[m_LargeGapEvent]->Fill(Sgn*particle.eta(), PhiDistance);
                    if(abs(PhiDistance) <= DEF_1DCorrelationPhiTolerance) h1F_1DCorrelation_eta_pT_2_4[m_LargeGapEvent]->Fill(Sgn*particle.eta());
                }
                else if(particle.pt() < 60 && particle.pt() >= 4.0)
                {
                    h2F_2DCorrelation_eta_vs_dphi_pT_4_6[m_LargeGapEvent]->Fill(Sgn*particle.eta(), PhiDistance);
                    if(abs(PhiDistance) <= DEF_1DCorrelationPhiTolerance) h1F_1DCorrelation_eta_pT_4_6[m_LargeGapEvent]->Fill(Sgn*particle.eta());
                }
                
            }

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

    //divide by event numbers
    if(LargeGapEventCounter > 0 || SmallGapEventCounter > 0)
    {
        for(int iGap = 0; iGap <=1; iGap++)
        {
            float ScaleFactor;
            (iGap) ? (ScaleFactor = 1./(float)LargeGapEventCounter) : (ScaleFactor = 1./(float)SmallGapEventCounter);
            h2F_2DCorrelation_eta_vs_dphi_pT_0_1[iGap]->Scale(ScaleFactor);
            h2F_2DCorrelation_eta_vs_dphi_pT_1_2[iGap]->Scale(ScaleFactor);
            h2F_2DCorrelation_eta_vs_dphi_pT_2_4[iGap]->Scale(ScaleFactor);
            h2F_2DCorrelation_eta_vs_dphi_pT_4_6[iGap]->Scale(ScaleFactor);
        
            h1F_1DCorrelation_eta_pT_0_1[iGap]->Scale(ScaleFactor);
            h1F_1DCorrelation_eta_pT_1_2[iGap]->Scale(ScaleFactor);
            h1F_1DCorrelation_eta_pT_2_4[iGap]->Scale(ScaleFactor);
            h1F_1DCorrelation_eta_pT_4_6[iGap]->Scale(ScaleFactor);
        
        }
    }

    //produce differences
    h1F_1DCorrelationDifference_eta_pT_0_1->Add(h1F_1DCorrelation_eta_pT_0_1[1], +1);
    h1F_1DCorrelationDifference_eta_pT_0_1->Add(h1F_1DCorrelation_eta_pT_0_1[0], -1);

    h1F_1DCorrelationDifference_eta_pT_1_2->Add(h1F_1DCorrelation_eta_pT_1_2[1], +1);
    h1F_1DCorrelationDifference_eta_pT_1_2->Add(h1F_1DCorrelation_eta_pT_1_2[0], -1);

    h1F_1DCorrelationDifference_eta_pT_2_4->Add(h1F_1DCorrelation_eta_pT_2_4[1], +1);
    h1F_1DCorrelationDifference_eta_pT_2_4->Add(h1F_1DCorrelation_eta_pT_2_4[0], -1);

    h1F_1DCorrelationDifference_eta_pT_4_6->Add(h1F_1DCorrelation_eta_pT_4_6[1], +1);
    h1F_1DCorrelationDifference_eta_pT_4_6->Add(h1F_1DCorrelation_eta_pT_4_6[0], -1);

    //scale by bin width
    for(int iGap = 0; iGap <=1; iGap++)
    {
        h2F_2DCorrelation_eta_vs_dphi_pT_0_1[iGap]->Scale(1./((float)h2F_2DCorrelation_eta_vs_dphi_pT_0_1[iGap]->GetXaxis()->GetBinWidth(1) * (float)h2F_2DCorrelation_eta_vs_dphi_pT_0_1[iGap]->GetYaxis()->GetBinWidth(1)));
        h2F_2DCorrelation_eta_vs_dphi_pT_1_2[iGap]->Scale(1./((float)h2F_2DCorrelation_eta_vs_dphi_pT_1_2[iGap]->GetXaxis()->GetBinWidth(1) * (float)h2F_2DCorrelation_eta_vs_dphi_pT_1_2[iGap]->GetYaxis()->GetBinWidth(1)));
        h2F_2DCorrelation_eta_vs_dphi_pT_2_4[iGap]->Scale(1./((float)h2F_2DCorrelation_eta_vs_dphi_pT_2_4[iGap]->GetXaxis()->GetBinWidth(1) * (float)h2F_2DCorrelation_eta_vs_dphi_pT_2_4[iGap]->GetYaxis()->GetBinWidth(1)));
        h2F_2DCorrelation_eta_vs_dphi_pT_4_6[iGap]->Scale(1./((float)h2F_2DCorrelation_eta_vs_dphi_pT_4_6[iGap]->GetXaxis()->GetBinWidth(1) * (float)h2F_2DCorrelation_eta_vs_dphi_pT_4_6[iGap]->GetYaxis()->GetBinWidth(1)));
    
        h1F_1DCorrelation_eta_pT_0_1[iGap]->Scale(1./(float)h1F_1DCorrelation_eta_pT_0_1[iGap]->GetXaxis()->GetBinWidth(1));
        h1F_1DCorrelation_eta_pT_1_2[iGap]->Scale(1./(float)h1F_1DCorrelation_eta_pT_1_2[iGap]->GetXaxis()->GetBinWidth(1));
        h1F_1DCorrelation_eta_pT_2_4[iGap]->Scale(1./(float)h1F_1DCorrelation_eta_pT_2_4[iGap]->GetXaxis()->GetBinWidth(1));
        h1F_1DCorrelation_eta_pT_4_6[iGap]->Scale(1./(float)h1F_1DCorrelation_eta_pT_4_6[iGap]->GetXaxis()->GetBinWidth(1));
    
    }
    h1F_1DCorrelationDifference_eta_pT_0_1->Scale(1./(float)h1F_1DCorrelationDifference_eta_pT_0_1->GetXaxis()->GetBinWidth(1));
    h1F_1DCorrelationDifference_eta_pT_1_2->Scale(1./(float)h1F_1DCorrelationDifference_eta_pT_1_2->GetXaxis()->GetBinWidth(1));
    h1F_1DCorrelationDifference_eta_pT_2_4->Scale(1./(float)h1F_1DCorrelationDifference_eta_pT_2_4->GetXaxis()->GetBinWidth(1));
    h1F_1DCorrelationDifference_eta_pT_4_6->Scale(1./(float)h1F_1DCorrelationDifference_eta_pT_4_6->GetXaxis()->GetBinWidth(1));

    //set axis titles
    for(int iGap = 0; iGap <=1; iGap++)
    {
        h2F_2DCorrelation_eta_vs_dphi_pT_0_1[iGap]->SetTitle("");
        h2F_2DCorrelation_eta_vs_dphi_pT_0_1[iGap]->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
        h2F_2DCorrelation_eta_vs_dphi_pT_0_1[iGap]->GetXaxis()->SetTitle("#eta");
        h2F_2DCorrelation_eta_vs_dphi_pT_0_1[iGap]->GetYaxis()->SetTitleSize(DEF_AxisLabelSize);
        h2F_2DCorrelation_eta_vs_dphi_pT_0_1[iGap]->GetYaxis()->SetTitle("#Delta #phi [rad]");
        h2F_2DCorrelation_eta_vs_dphi_pT_0_1[iGap]->GetZaxis()->SetTitleSize(DEF_AxisLabelSize);
        h2F_2DCorrelation_eta_vs_dphi_pT_0_1[iGap]->GetZaxis()->SetTitle("#frac{1}{N_{Events}} #frac{d N}{d #Delta #phi d #eta}");

        h2F_2DCorrelation_eta_vs_dphi_pT_1_2[iGap]->SetTitle("");
        h2F_2DCorrelation_eta_vs_dphi_pT_1_2[iGap]->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
        h2F_2DCorrelation_eta_vs_dphi_pT_1_2[iGap]->GetXaxis()->SetTitle("#eta");
        h2F_2DCorrelation_eta_vs_dphi_pT_1_2[iGap]->GetYaxis()->SetTitleSize(DEF_AxisLabelSize);
        h2F_2DCorrelation_eta_vs_dphi_pT_1_2[iGap]->GetYaxis()->SetTitle("#Delta #phi [rad]");
        h2F_2DCorrelation_eta_vs_dphi_pT_1_2[iGap]->GetZaxis()->SetTitleSize(DEF_AxisLabelSize);
        h2F_2DCorrelation_eta_vs_dphi_pT_1_2[iGap]->GetZaxis()->SetTitle("#frac{1}{N_{Events}} #frac{d N}{d #Delta #phi d #eta}");

        h2F_2DCorrelation_eta_vs_dphi_pT_2_4[iGap]->SetTitle("");
        h2F_2DCorrelation_eta_vs_dphi_pT_2_4[iGap]->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
        h2F_2DCorrelation_eta_vs_dphi_pT_2_4[iGap]->GetXaxis()->SetTitle("#eta");
        h2F_2DCorrelation_eta_vs_dphi_pT_2_4[iGap]->GetYaxis()->SetTitleSize(DEF_AxisLabelSize);
        h2F_2DCorrelation_eta_vs_dphi_pT_2_4[iGap]->GetYaxis()->SetTitle("#Delta #phi [rad]");
        h2F_2DCorrelation_eta_vs_dphi_pT_2_4[iGap]->GetZaxis()->SetTitleSize(DEF_AxisLabelSize);
        h2F_2DCorrelation_eta_vs_dphi_pT_2_4[iGap]->GetZaxis()->SetTitle("#frac{1}{N_{Events}} #frac{d N}{d #Delta #phi d #eta}");

        h2F_2DCorrelation_eta_vs_dphi_pT_4_6[iGap]->SetTitle("");
        h2F_2DCorrelation_eta_vs_dphi_pT_4_6[iGap]->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
        h2F_2DCorrelation_eta_vs_dphi_pT_4_6[iGap]->GetXaxis()->SetTitle("#eta");
        h2F_2DCorrelation_eta_vs_dphi_pT_4_6[iGap]->GetYaxis()->SetTitleSize(DEF_AxisLabelSize);
        h2F_2DCorrelation_eta_vs_dphi_pT_4_6[iGap]->GetYaxis()->SetTitle("#Delta #phi [rad]");
        h2F_2DCorrelation_eta_vs_dphi_pT_4_6[iGap]->GetZaxis()->SetTitleSize(DEF_AxisLabelSize);
        h2F_2DCorrelation_eta_vs_dphi_pT_4_6[iGap]->GetZaxis()->SetTitle("#frac{1}{N_{Events}} #frac{d N}{d #Delta #phi d #eta}");

        h1F_1DCorrelation_eta_pT_0_1[iGap]->SetTitle("");
        h1F_1DCorrelation_eta_pT_0_1[iGap]->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
        h1F_1DCorrelation_eta_pT_0_1[iGap]->GetXaxis()->SetTitle("#eta");
        h1F_1DCorrelation_eta_pT_0_1[iGap]->GetYaxis()->SetTitleSize(DEF_AxisLabelSize);
        h1F_1DCorrelation_eta_pT_0_1[iGap]->GetYaxis()->SetTitle("#frac{1}{N_{Events}} #frac{d N}{d #eta}");

        h1F_1DCorrelation_eta_pT_1_2[iGap]->SetTitle("");
        h1F_1DCorrelation_eta_pT_1_2[iGap]->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
        h1F_1DCorrelation_eta_pT_1_2[iGap]->GetXaxis()->SetTitle("#eta");
        h1F_1DCorrelation_eta_pT_1_2[iGap]->GetYaxis()->SetTitleSize(DEF_AxisLabelSize);
        h1F_1DCorrelation_eta_pT_1_2[iGap]->GetYaxis()->SetTitle("#frac{1}{N_{Events}} #frac{d N}{d #eta}");

        h1F_1DCorrelation_eta_pT_2_4[iGap]->SetTitle("");
        h1F_1DCorrelation_eta_pT_2_4[iGap]->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
        h1F_1DCorrelation_eta_pT_2_4[iGap]->GetXaxis()->SetTitle("#eta");
        h1F_1DCorrelation_eta_pT_2_4[iGap]->GetYaxis()->SetTitleSize(DEF_AxisLabelSize);
        h1F_1DCorrelation_eta_pT_2_4[iGap]->GetYaxis()->SetTitle("#frac{1}{N_{Events}} #frac{d N}{d #eta}");

        h1F_1DCorrelation_eta_pT_4_6[iGap]->SetTitle("");
        h1F_1DCorrelation_eta_pT_4_6[iGap]->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
        h1F_1DCorrelation_eta_pT_4_6[iGap]->GetXaxis()->SetTitle("#eta");
        h1F_1DCorrelation_eta_pT_4_6[iGap]->GetYaxis()->SetTitleSize(DEF_AxisLabelSize);
        h1F_1DCorrelation_eta_pT_4_6[iGap]->GetYaxis()->SetTitle("#frac{1}{N_{Events}} #frac{d N}{d #eta}");
    
    }

    h1F_1DCorrelationDifference_eta_pT_0_1->SetTitle("");
    h1F_1DCorrelationDifference_eta_pT_0_1->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
    h1F_1DCorrelationDifference_eta_pT_0_1->GetXaxis()->SetTitle("#eta");
    h1F_1DCorrelationDifference_eta_pT_0_1->GetYaxis()->SetTitleSize(DEF_AxisLabelSize);
    h1F_1DCorrelationDifference_eta_pT_0_1->GetYaxis()->SetTitle("#frac{1}{N_{Large Gap}} #frac{d N_{Large Gap}}{d #eta} - #frac{1}{N_{Small Gap}} #frac{d N_{SmallGap}}{d #eta}");

    h1F_1DCorrelationDifference_eta_pT_1_2->SetTitle("");
    h1F_1DCorrelationDifference_eta_pT_1_2->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
    h1F_1DCorrelationDifference_eta_pT_1_2->GetXaxis()->SetTitle("#eta");
    h1F_1DCorrelationDifference_eta_pT_1_2->GetYaxis()->SetTitleSize(DEF_AxisLabelSize);
    h1F_1DCorrelationDifference_eta_pT_1_2->GetYaxis()->SetTitle("#frac{1}{N_{Large Gap}} #frac{d N_{Large Gap}}{d #eta} - #frac{1}{N_{Small Gap}} #frac{d N_{SmallGap}}{d #eta}");

    h1F_1DCorrelationDifference_eta_pT_2_4->SetTitle("");
    h1F_1DCorrelationDifference_eta_pT_2_4->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
    h1F_1DCorrelationDifference_eta_pT_2_4->GetXaxis()->SetTitle("#eta");
    h1F_1DCorrelationDifference_eta_pT_2_4->GetYaxis()->SetTitleSize(DEF_AxisLabelSize);
    h1F_1DCorrelationDifference_eta_pT_2_4->GetYaxis()->SetTitle("#frac{1}{N_{Large Gap}} #frac{d N_{Large Gap}}{d #eta} - #frac{1}{N_{Small Gap}} #frac{d N_{SmallGap}}{d #eta}");

    h1F_1DCorrelationDifference_eta_pT_4_6->SetTitle("");
    h1F_1DCorrelationDifference_eta_pT_4_6->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
    h1F_1DCorrelationDifference_eta_pT_4_6->GetXaxis()->SetTitle("#eta");
    h1F_1DCorrelationDifference_eta_pT_4_6->GetYaxis()->SetTitleSize(DEF_AxisLabelSize);
    h1F_1DCorrelationDifference_eta_pT_4_6->GetYaxis()->SetTitle("#frac{1}{N_{Large Gap}} #frac{d N_{Large Gap}}{d #eta} - #frac{1}{N_{Small Gap}} #frac{d N_{SmallGap}}{d #eta}");


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /*
    /
    /                                                     
    /                                                               PRINT ANALYSIS STATS
    /                                                     
    /                                                     
    /                                                     
    /
    /                                                    
    /
    /
    /
    *//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    PrintInfo("Overall Analysis Stats");

    cout << "Accepted Events / Total Events: " << 100.0 * ((float)(LargeGapEventCounter + SmallGapEventCounter) / (float)TotalEvents) << " percent" << endl;

    //save to counter histo
    CounterHisto->Fill("Input Events", TotalEvents);
    CounterHisto->Fill("Large Gap Events", LargeGapEventCounter);
    CounterHisto->Fill("Small Gap Events", SmallGapEventCounter);

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
    TString str_out =  "./" + (TString)RunNumberAsString.Data() + "_DijetAna.root";

    TFile* OutputFile;
    if(gSystem->AccessPathName(str_out))
    {
        OutputFile = new TFile(str_out,"RECREATE");
        OutputFile->cd();
        cout << "Output file is " << str_out << endl;

        for(int iGap = 0; iGap <=1; iGap++)
        {
            h2F_2DCorrelation_eta_vs_dphi_pT_0_1[iGap]->Write();
            h2F_2DCorrelation_eta_vs_dphi_pT_1_2[iGap]->Write();
            h2F_2DCorrelation_eta_vs_dphi_pT_2_4[iGap]->Write();
            h2F_2DCorrelation_eta_vs_dphi_pT_4_6[iGap]->Write();

            h1F_1DCorrelation_eta_pT_0_1[iGap]->Write();
            h1F_1DCorrelation_eta_pT_1_2[iGap]->Write();
            h1F_1DCorrelation_eta_pT_2_4[iGap]->Write();
            h1F_1DCorrelation_eta_pT_4_6[iGap]->Write();
        }
        h1F_1DCorrelationDifference_eta_pT_0_1->Write();
        h1F_1DCorrelationDifference_eta_pT_1_2->Write();
        h1F_1DCorrelationDifference_eta_pT_2_4->Write();
        h1F_1DCorrelationDifference_eta_pT_4_6->Write();
        CounterHisto->Write();

        OutputFile->Close();
    }
    else
    {
        OutputFile = TFile::Open(str_out.Data());
        cout << "ADD TO EXISTING OUTPUT FILE" << endl;
        OutputFile->cd();

        //get
        TH2F* h2F_2DCorrelation_eta_vs_dphi_pT_0_1_SmallGap_Sum = (TH2F*) OutputFile->Get("h2F_2DCorrelation_eta_vs_dphi_pT_0_1_SmallGap");
        h2F_2DCorrelation_eta_vs_dphi_pT_0_1_SmallGap_Sum->Add(h2F_2DCorrelation_eta_vs_dphi_pT_0_1[0]);
        h2F_2DCorrelation_eta_vs_dphi_pT_0_1_SmallGap_Sum->Scale(0.5);

        TH2F* h2F_2DCorrelation_eta_vs_dphi_pT_1_2_SmallGap_Sum = (TH2F*) OutputFile->Get("h2F_2DCorrelation_eta_vs_dphi_pT_1_2_SmallGap");
        h2F_2DCorrelation_eta_vs_dphi_pT_1_2_SmallGap_Sum->Add(h2F_2DCorrelation_eta_vs_dphi_pT_1_2[0]);
        h2F_2DCorrelation_eta_vs_dphi_pT_1_2_SmallGap_Sum->Scale(0.5);

        TH2F* h2F_2DCorrelation_eta_vs_dphi_pT_2_4_SmallGap_Sum = (TH2F*) OutputFile->Get("h2F_2DCorrelation_eta_vs_dphi_pT_2_4_SmallGap");
        h2F_2DCorrelation_eta_vs_dphi_pT_2_4_SmallGap_Sum->Add(h2F_2DCorrelation_eta_vs_dphi_pT_2_4[0]);
        h2F_2DCorrelation_eta_vs_dphi_pT_2_4_SmallGap_Sum->Scale(0.5);

        TH2F* h2F_2DCorrelation_eta_vs_dphi_pT_4_6_SmallGap_Sum = (TH2F*) OutputFile->Get("h2F_2DCorrelation_eta_vs_dphi_pT_4_6_SmallGap");
        h2F_2DCorrelation_eta_vs_dphi_pT_4_6_SmallGap_Sum->Add(h2F_2DCorrelation_eta_vs_dphi_pT_4_6[0]);
        h2F_2DCorrelation_eta_vs_dphi_pT_4_6_SmallGap_Sum->Scale(0.5);

        TH2F* h2F_2DCorrelation_eta_vs_dphi_pT_0_1_LargeGap_Sum = (TH2F*) OutputFile->Get("h2F_2DCorrelation_eta_vs_dphi_pT_0_1_LargeGap");
        h2F_2DCorrelation_eta_vs_dphi_pT_0_1_LargeGap_Sum->Add(h2F_2DCorrelation_eta_vs_dphi_pT_0_1[1]);
        h2F_2DCorrelation_eta_vs_dphi_pT_0_1_LargeGap_Sum->Scale(0.5);

        TH2F* h2F_2DCorrelation_eta_vs_dphi_pT_1_2_LargeGap_Sum = (TH2F*) OutputFile->Get("h2F_2DCorrelation_eta_vs_dphi_pT_1_2_LargeGap");
        h2F_2DCorrelation_eta_vs_dphi_pT_1_2_LargeGap_Sum->Add(h2F_2DCorrelation_eta_vs_dphi_pT_1_2[1]);
        h2F_2DCorrelation_eta_vs_dphi_pT_1_2_LargeGap_Sum->Scale(0.5);

        TH2F* h2F_2DCorrelation_eta_vs_dphi_pT_2_4_LargeGap_Sum = (TH2F*) OutputFile->Get("h2F_2DCorrelation_eta_vs_dphi_pT_2_4_LargeGap");
        h2F_2DCorrelation_eta_vs_dphi_pT_2_4_LargeGap_Sum->Add(h2F_2DCorrelation_eta_vs_dphi_pT_2_4[1]);
        h2F_2DCorrelation_eta_vs_dphi_pT_2_4_LargeGap_Sum->Scale(0.5);

        TH2F* h2F_2DCorrelation_eta_vs_dphi_pT_4_6_LargeGap_Sum = (TH2F*) OutputFile->Get("h2F_2DCorrelation_eta_vs_dphi_pT_4_6_LargeGap");
        h2F_2DCorrelation_eta_vs_dphi_pT_4_6_LargeGap_Sum->Add(h2F_2DCorrelation_eta_vs_dphi_pT_4_6[1]);
        h2F_2DCorrelation_eta_vs_dphi_pT_4_6_LargeGap_Sum->Scale(0.5);

        TH1F* h1F_1DCorrelation_eta_pT_0_1_SmallGap_Sum = (TH1F*) OutputFile->Get("h1F_1DCorrelation_eta_pT_0_1_SmallGap");
        h1F_1DCorrelation_eta_pT_0_1_SmallGap_Sum->Add(h1F_1DCorrelation_eta_pT_0_1[0]);
        h1F_1DCorrelation_eta_pT_0_1_SmallGap_Sum->Scale(0.5);

        TH1F* h1F_1DCorrelation_eta_pT_1_2_SmallGap_Sum = (TH1F*) OutputFile->Get("h1F_1DCorrelation_eta_pT_1_2_SmallGap");
        h1F_1DCorrelation_eta_pT_1_2_SmallGap_Sum->Add(h1F_1DCorrelation_eta_pT_1_2[0]);
        h1F_1DCorrelation_eta_pT_1_2_SmallGap_Sum->Scale(0.5);

        TH1F* h1F_1DCorrelation_eta_pT_2_4_SmallGap_Sum = (TH1F*) OutputFile->Get("h1F_1DCorrelation_eta_pT_2_4_SmallGap");
        h1F_1DCorrelation_eta_pT_2_4_SmallGap_Sum->Add(h1F_1DCorrelation_eta_pT_2_4[0]);
        h1F_1DCorrelation_eta_pT_2_4_SmallGap_Sum->Scale(0.5);

        TH1F* h1F_1DCorrelation_eta_pT_4_6_SmallGap_Sum = (TH1F*) OutputFile->Get("h1F_1DCorrelation_eta_pT_4_6_SmallGap");
        h1F_1DCorrelation_eta_pT_4_6_SmallGap_Sum->Add(h1F_1DCorrelation_eta_pT_4_6[0]);
        h1F_1DCorrelation_eta_pT_4_6_SmallGap_Sum->Scale(0.5);

        TH1F* h1F_1DCorrelation_eta_pT_0_1_LargeGap_Sum = (TH1F*) OutputFile->Get("h1F_1DCorrelation_eta_pT_0_1_LargeGap");
        h1F_1DCorrelation_eta_pT_0_1_LargeGap_Sum->Add(h1F_1DCorrelation_eta_pT_0_1[1]);
        h1F_1DCorrelation_eta_pT_0_1_LargeGap_Sum->Scale(0.5);

        TH1F* h1F_1DCorrelation_eta_pT_1_2_LargeGap_Sum = (TH1F*) OutputFile->Get("h1F_1DCorrelation_eta_pT_1_2_LargeGap");
        h1F_1DCorrelation_eta_pT_1_2_LargeGap_Sum->Add(h1F_1DCorrelation_eta_pT_1_2[1]);
        h1F_1DCorrelation_eta_pT_1_2_LargeGap_Sum->Scale(0.5);

        TH1F* h1F_1DCorrelation_eta_pT_2_4_LargeGap_Sum = (TH1F*) OutputFile->Get("h1F_1DCorrelation_eta_pT_2_4_LargeGap");
        h1F_1DCorrelation_eta_pT_2_4_LargeGap_Sum->Add(h1F_1DCorrelation_eta_pT_2_4[1]);
        h1F_1DCorrelation_eta_pT_2_4_LargeGap_Sum->Scale(0.5);

        TH1F* h1F_1DCorrelation_eta_pT_4_6_LargeGap_Sum = (TH1F*) OutputFile->Get("h1F_1DCorrelation_eta_pT_4_6_LargeGap");
        h1F_1DCorrelation_eta_pT_4_6_LargeGap_Sum->Add(h1F_1DCorrelation_eta_pT_4_6[1]);
        h1F_1DCorrelation_eta_pT_4_6_LargeGap_Sum->Scale(0.5);
        
        TH1F* h1F_1DCorrelationDifference_eta_pT_0_1_Sum = (TH1F*) OutputFile->Get("h1F_1DCorrelationDifference_eta_pT_0_1");
        h1F_1DCorrelationDifference_eta_pT_0_1_Sum->Add(h1F_1DCorrelationDifference_eta_pT_0_1);
        h1F_1DCorrelationDifference_eta_pT_0_1_Sum->Scale(0.5);

        TH1F* h1F_1DCorrelationDifference_eta_pT_1_2_Sum = (TH1F*) OutputFile->Get("h1F_1DCorrelationDifference_eta_pT_1_2");
        h1F_1DCorrelationDifference_eta_pT_1_2_Sum->Add(h1F_1DCorrelationDifference_eta_pT_1_2);
        h1F_1DCorrelationDifference_eta_pT_1_2_Sum->Scale(0.5);

        TH1F* h1F_1DCorrelationDifference_eta_pT_2_4_Sum = (TH1F*) OutputFile->Get("h1F_1DCorrelationDifference_eta_pT_2_4");
        h1F_1DCorrelationDifference_eta_pT_2_4_Sum->Add(h1F_1DCorrelationDifference_eta_pT_2_4);
        h1F_1DCorrelationDifference_eta_pT_2_4_Sum->Scale(0.5);

        TH1F* h1F_1DCorrelationDifference_eta_pT_4_6_Sum = (TH1F*) OutputFile->Get("h1F_1DCorrelationDifference_eta_pT_4_6");
        h1F_1DCorrelationDifference_eta_pT_4_6_Sum->Add(h1F_1DCorrelationDifference_eta_pT_4_6);
        h1F_1DCorrelationDifference_eta_pT_4_6_Sum->Scale(0.5);

        TH1I* CounterHisto_Sum = (TH1I*) OutputFile->Get("CounterHisto");
        CounterHisto_Sum->Add(CounterHisto);

        OutputFile = new TFile(str_out.Data(),"RECREATE");
        OutputFile->cd();
        h2F_2DCorrelation_eta_vs_dphi_pT_0_1_SmallGap_Sum->Write();
        h2F_2DCorrelation_eta_vs_dphi_pT_1_2_SmallGap_Sum->Write();
        h2F_2DCorrelation_eta_vs_dphi_pT_2_4_SmallGap_Sum->Write();
        h2F_2DCorrelation_eta_vs_dphi_pT_4_6_SmallGap_Sum->Write();
        h2F_2DCorrelation_eta_vs_dphi_pT_0_1_LargeGap_Sum->Write();
        h2F_2DCorrelation_eta_vs_dphi_pT_1_2_LargeGap_Sum->Write();
        h2F_2DCorrelation_eta_vs_dphi_pT_2_4_LargeGap_Sum->Write();
        h2F_2DCorrelation_eta_vs_dphi_pT_4_6_LargeGap_Sum->Write();
        h1F_1DCorrelation_eta_pT_0_1_SmallGap_Sum->Write();
        h1F_1DCorrelation_eta_pT_1_2_SmallGap_Sum->Write();
        h1F_1DCorrelation_eta_pT_2_4_SmallGap_Sum->Write();
        h1F_1DCorrelation_eta_pT_4_6_SmallGap_Sum->Write();
        h1F_1DCorrelation_eta_pT_0_1_LargeGap_Sum->Write();
        h1F_1DCorrelation_eta_pT_1_2_LargeGap_Sum->Write();
        h1F_1DCorrelation_eta_pT_2_4_LargeGap_Sum->Write();
        h1F_1DCorrelation_eta_pT_4_6_LargeGap_Sum->Write();
        h1F_1DCorrelationDifference_eta_pT_0_1_Sum->Write();
        h1F_1DCorrelationDifference_eta_pT_1_2_Sum->Write();
        h1F_1DCorrelationDifference_eta_pT_2_4_Sum->Write();
        h1F_1DCorrelationDifference_eta_pT_4_6_Sum->Write();
        CounterHisto_Sum->Write();

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