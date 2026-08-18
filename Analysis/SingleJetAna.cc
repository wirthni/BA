#include "./ReadTree.cc"

using namespace fastjet;
using namespace std;

vector<PseudoJet> FindJets(const vector<PseudoJet> vec_particles);

/* Parameters:
i_InputFile: Format xxx.root
*/
Int_t SingleJetAna(TString i_InputFile = "in.root")
{
    //CONFIGURE
    bool i_LocalDebugRun = true;
    bool i_MakeCrosscheckAnalysis = false;
    bool i_UseHadronInstead = false;
    int i_NoOfEvents = -1;

    gStyle->SetOptStat(0);
    SetRootGraphicStyle();

    #define DEF_BinningPerUnit 100
    #define DEF_JetRadius 0.2
    #define DEF_OutputEventOverviews false
    #define DEF_HadLeadingPt 20.0
    #define DEF_BackgroundLimit 7.0 //GeV
    #define DEF_CorrelationMinPt 0.0
    #define DEF_CorrelationMaxPt 2.0
    #define DEF_MaxPartilclesPerJet 200
    #define DEF_MaxJetPt 150.0
    #define DEF_AxisLabelSize 0.05
    #define DEF_HistoTitleSize 0.1
    #define DEF_EventTolerance_Multiplicity 10
    #define DEF_EventTolerance_Psi2 Pi/8
    #define DEF_RecoilCorrelationFrameSize 0.5
    #define DEF_MaxParticlePtInCorrelation 5.0
    #define DEF_Epsilon 0.001

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
    TH3F* h3F_dphi_vs_deta_vs_pt_CorrelationDifference = new TH3F("h3F_dphi_vs_deta_vs_pt_CorrelationDifference", "h3F_dphi_vs_deta_vs_pt_CorrelationDifference",
                DEF_BinningPerUnit * DEF_RecoilCorrelationFrameSize / 2, -DEF_RecoilCorrelationFrameSize, DEF_RecoilCorrelationFrameSize,
                DEF_BinningPerUnit * DEF_RecoilCorrelationFrameSize / 2, -DEF_RecoilCorrelationFrameSize, DEF_RecoilCorrelationFrameSize,
                 0.1 * DEF_BinningPerUnit * DEF_MaxParticlePtInCorrelation, 0, DEF_MaxParticlePtInCorrelation);//temporarily filled
    TH2F* h2F_RadialCorrelationDifference = new TH2F("h2F_RadialCorrelationDifference", "h2F_RadialCorrelationDifference", 25, 0, DEF_RecoilCorrelationFrameSize, 30, 0, 3);
    TH1F* h1F_NoOfParticlesInRecoilArea = new TH1F("h1F_NoOfParticlesInRecoilArea", "h1F_NoOfParticlesInRecoilArea", 100, 0, 300);
    TH1F* h1F_NoOfParticlesInReferenceArea = new TH1F("h1F_NoOfParticlesInReferenceArea", "h1F_NoOfParticlesInReferenceArea", 100, 0, 300);
    TH1F* h1F_N_ev_minus_N_ref = new TH1F("h1F_N_ev_minus_N_ref", "h1F_N_ev_minus_N_ref", 20, -10, 10);
    TH1F* h1F_NoOfHighPtParticles = new TH1F("h1F_NoOfHighPtParticles", "h1F_NoOfHighPtParticles", 50, 0, 50);
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
        vector<PseudoJet> ParticleVector[entries_col];//holds particles j = 0,1,2,...,N for event i at position [i][j]
        struct struct_EventProperties
        {
            int Multiplicity = 0;
            float Psi2 = -100;
            vector<PseudoJet> HighPtParticles;
            bool IsHighPtJetEvent = false;
            int SuitableReferenceEvent = -1;
            
        };
        struct_EventProperties EventProperties[entries_col];

        // Make an unordered map particle -> collision
        for(uint64_t i_Track = 0; i_Track < entries_track; i_Track ++)
        {
            Particle track = Particle::Read(tracks, i_Track);
            CollisionMap[track.ColID].push_back(i_Track);
            //prepare particles for jetfinder
            ParticleVector[track.ColID].push_back(PseudoJet(track.px, track.py, track.pz, track.E));
            if(track.Pt > DEF_BackgroundLimit)
            {
                EventProperties[track.ColID].HighPtParticles.push_back(PseudoJet(track.px, track.py, track.pz, track.E));
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

        //APPLICATION VARIABLES
        int HighPtEventCounter = 0;
        int PairedEventsCounter = 0;

        //get run number to retrieve event plane correction file
        Collision collision = Collision::Read(collisions, 0);
        RunNumber = collision.RunNumber;
        TString RunNumberAsString;
        RunNumberAsString.Form("%d", RunNumber);

        //differ between the LHC passes
        TString QVectorCorrectionFileName;
        if(i_InputFile.Contains("LHC23"))
        {
            QVectorCorrectionFileName =  "./RUN3_Data/QVecCalibRoots/LHC23_PbPb_pass5/" + (TString)RunNumberAsString.Data() + "_Calibration.root";
        }
        else if(i_InputFile.Contains("LHC25"))
        {
            QVectorCorrectionFileName =  "./RUN3_Data/QVecCalibRoots/LHC25_PbPb_pass1/" + (TString)RunNumberAsString.Data() + "_Calibration.root";
        }

        cout << "Q vector correction file is " << QVectorCorrectionFileName << endl;
        

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /*
        /
        /                                                     
        /                                           Get EventProperties for each event to check which events
        /                                           are comparable for background subtraction  
        /                                                     
        /                                                     
        /                                                     
        /
        /                                                    
        /
        /
        /
        *//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        PrintInfo("Get event properties");

        if(i_NoOfEvents != -1) entries_col = i_NoOfEvents;
        for (int iEvent = 0; iEvent < entries_col; iEvent++) {

            PrintProgress(iEvent, entries_col);

            Collision collision = Collision::Read(collisions, iEvent);
            if(CollisionMap[collision.ColID].size() == 0 || ParticleVector[collision.ColID].size() == 0){
                cout << "No tracks in collision " << collision.ColID << "! Skipping." << endl;
                N_events_skipped++;
                continue;
            }

            //sort high pt particles by pt
            if(EventProperties[iEvent].HighPtParticles.size() >= 1)EventProperties[iEvent].HighPtParticles = sorted_by_pt(EventProperties[iEvent].HighPtParticles);

            EventProperties[iEvent].Multiplicity = CollisionMap[collision.ColID].size();

            //event plane correction
            EventProperties[iEvent].Psi2 = GetEventPlaneInfo(collision.ColID, tracks, CollisionMap, QVectorCorrectionFileName);

            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            /*
            /
            /                                                     
            /                                           Check if the current event is suitable for diffusion wake ana  
            /                                                     
            /                                                     
            /                                                     
            /
            /                                                    
            /
            /
            /
            *//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            if(EventProperties[iEvent].HighPtParticles.size() < 1) continue;

            if(EventProperties[iEvent].HighPtParticles[0].pt() < DEF_HadLeadingPt) continue;

            if(abs(EventProperties[iEvent].HighPtParticles[0].eta()) > 0.9 - DEF_RecoilCorrelationFrameSize) continue;

            //check that there is no jet close to the direct recoil site
            float RecoilSiteEta = EventProperties[iEvent].HighPtParticles[0].eta();
            if(i_MakeCrosscheckAnalysis) RecoilSiteEta*= -1;

            float RecoilSitePhi;
            if(EventProperties[iEvent].HighPtParticles[0].phi_std() >= 0)
            {
                RecoilSitePhi = EventProperties[iEvent].HighPtParticles[0].phi_std() - Pi;
            }
            else
            {
                RecoilSitePhi = EventProperties[iEvent].HighPtParticles[0].phi_std() + Pi;
            }
            for(int iJet = 0; iJet < EventProperties[iEvent].HighPtParticles.size(); iJet++)
            {
                float EtaDistance = EventProperties[iEvent].HighPtParticles[iJet].eta() - RecoilSiteEta;
                float PhiDistance = EventProperties[iEvent].HighPtParticles[iJet].phi_std() - RecoilSitePhi;
                if(PhiDistance > Pi) PhiDistance -= 2*Pi;
                else if(PhiDistance < -Pi) PhiDistance += 2*Pi;
                float Distance = sqrt(pow(EtaDistance, 2) + pow(PhiDistance, 2));
                if(Distance < 2*DEF_JetRadius) continue;
            }

            EventProperties[iEvent].IsHighPtJetEvent = true;
            HighPtEventCounter++;
            h1F_NoOfHighPtParticles->Fill(EventProperties[iEvent].HighPtParticles.size());

        }

        cout << "High Pt Events/Total Events: " << HighPtEventCounter << "/" << entries_col << ". That is " << ((float)HighPtEventCounter/(float)entries_col)*100.0 << " percent." << endl;

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /*
        /
        /                                                     
        /                                     Check which event is suitable for background subtraction for each
        /                                                           high pt jet event  
        /                                                     
        /                                                     
        /                                                     
        /
        /                                                    
        /
        /
        /
        *//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        PrintInfo("Check possible reference events");

        for(int iEvent = 0; iEvent < entries_col; iEvent++)
        {
            PrintProgress(iEvent, entries_col);

            if(!EventProperties[iEvent].IsHighPtJetEvent) continue;

            for(int iReferenceEvent = 0; iReferenceEvent < entries_col; iReferenceEvent++)
            {
                bool BadReferenceEvent = false;

                if(iReferenceEvent == iEvent) continue;

                if(abs(EventProperties[iReferenceEvent].Multiplicity - EventProperties[iEvent].Multiplicity) > DEF_EventTolerance_Multiplicity) continue;

                float PsiDiff = EventProperties[iReferenceEvent].Psi2 - EventProperties[iEvent].Psi2;
                if(PsiDiff > Pi) PsiDiff -= 2*Pi;
                else if(PsiDiff < -Pi) PsiDiff += 2*Pi;

                if(abs(PsiDiff) > DEF_EventTolerance_Psi2) continue;

                //check that there is no jet in the reference event where the reference data is taken
                float RecoilSiteEta = EventProperties[iEvent].HighPtParticles[0].eta();
                if(i_MakeCrosscheckAnalysis) RecoilSiteEta*= -1;

                float RecoilSitePhi;
                if(EventProperties[iEvent].HighPtParticles[0].phi_std() >= 0)
                {
                    RecoilSitePhi = EventProperties[iEvent].HighPtParticles[0].phi_std() - Pi;
                }
                else
                {
                    RecoilSitePhi = EventProperties[iEvent].HighPtParticles[0].phi_std() + Pi;
                }
                for(int iJet = 0; iJet < EventProperties[iReferenceEvent].HighPtParticles.size(); iJet++)
                {
                    //check recoil site is free
                    float EtaDistance = EventProperties[iReferenceEvent].HighPtParticles[iJet].eta() - RecoilSiteEta;
                    float PhiDistance = EventProperties[iReferenceEvent].HighPtParticles[iJet].phi_std() - RecoilSitePhi;
                    if(PhiDistance > Pi) PhiDistance -= 2*Pi;
                    else if(PhiDistance < -Pi) PhiDistance += 2*Pi;
                    float Distance = sqrt(pow(EtaDistance, 2) + pow(PhiDistance, 2));
                    if(Distance < 2*DEF_JetRadius)
                    {
                        BadReferenceEvent = true;
                        break;
                    }

                    //check jet site is free
                    EtaDistance = EventProperties[iReferenceEvent].HighPtParticles[iJet].eta() - EventProperties[iEvent].HighPtParticles[0].eta();
                    PhiDistance = EventProperties[iReferenceEvent].HighPtParticles[iJet].phi_std() - EventProperties[iEvent].HighPtParticles[0].phi_std();
                    if(PhiDistance > Pi) PhiDistance -= 2*Pi;
                    else if(PhiDistance < -Pi) PhiDistance += 2*Pi;
                    Distance = sqrt(pow(EtaDistance, 2) + pow(PhiDistance, 2));
                    if(Distance < 2*DEF_JetRadius)
                    {
                        BadReferenceEvent = true;
                        break;
                    }

                }
                if(BadReferenceEvent) continue;

                EventProperties[iEvent].SuitableReferenceEvent = iReferenceEvent;
                PairedEventsCounter++;
                break;
            }

        }

        cout << "Paired Events/High Pt Events: " << PairedEventsCounter << "/" << HighPtEventCounter << ". That is " << ((float)PairedEventsCounter/(float)HighPtEventCounter)*100.0 << " percent." << endl;

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /*
        /
        /                                                     
        /                                                           ANALYSIS EVENT LOOP  
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
            //DEBUG
            
            PrintProgress(iEvent, entries_col);

            Collision collision = Collision::Read(collisions, iEvent);
            if(CollisionMap[collision.ColID].size() == 0 || ParticleVector[collision.ColID].size() == 0){
                cout << "No tracks in collision " << collision.ColID << "! Skipping." << endl;
                N_events_skipped++;
                continue;
            }

            if(!EventProperties[iEvent].IsHighPtJetEvent) continue;
            if(EventProperties[iEvent].SuitableReferenceEvent == -1) continue;

            //Fill correlation histogram
            TH3F* h3F_dphi_vs_deta_vs_pt_RecoilCorrelation = new TH3F("h3F_dphi_vs_deta_vs_pt_RecoilCorrelation", "h3F_dphi_vs_deta_vs_pt_RecoilCorrelation",
                DEF_BinningPerUnit * DEF_RecoilCorrelationFrameSize / 2, -DEF_RecoilCorrelationFrameSize, DEF_RecoilCorrelationFrameSize,
                DEF_BinningPerUnit * DEF_RecoilCorrelationFrameSize / 2, -DEF_RecoilCorrelationFrameSize, DEF_RecoilCorrelationFrameSize,
                 0.1 * DEF_BinningPerUnit * DEF_MaxParticlePtInCorrelation, 0, DEF_MaxParticlePtInCorrelation);
            
            TH2F* h2F_RadialCorrelationRecoil = new TH2F("h2F_RadialCorrelationRecoil", "h2F_RadialCorrelationRecoil", 25, 0, DEF_RecoilCorrelationFrameSize, 30, 0, 3);
    
            float RecoilSiteEta = EventProperties[iEvent].HighPtParticles[0].eta();
            if(i_MakeCrosscheckAnalysis) RecoilSiteEta*= -1;

            float RecoilSitePhi;
            if(EventProperties[iEvent].HighPtParticles[0].phi_std() >= 0)
            {
                RecoilSitePhi = EventProperties[iEvent].HighPtParticles[0].phi_std() - Pi;
            }
            else
            {
                RecoilSitePhi = EventProperties[iEvent].HighPtParticles[0].phi_std() + Pi;
            }

            int NoOfParticlesInRecoilAreaCounter = 0;
            int NoOfParticlesInReferenceAreaCounter = 0;

            h1F_N_ev_minus_N_ref->Fill(EventProperties[iEvent].Multiplicity - EventProperties[EventProperties[iEvent].SuitableReferenceEvent].Multiplicity);

            for(const auto& particle : ParticleVector[collision.ColID])
            {

                float EtaDistance = particle.eta() - RecoilSiteEta;
                float PhiDistance = particle.phi_std() - RecoilSitePhi;
                if(PhiDistance > Pi) PhiDistance -= 2*Pi;
                else if(PhiDistance < -Pi) PhiDistance += 2*Pi;
                float Distance = sqrt(powf(EtaDistance, 2) + powf(PhiDistance, 2));
                if(abs(EtaDistance) <= DEF_RecoilCorrelationFrameSize && abs(PhiDistance) <= DEF_RecoilCorrelationFrameSize)
                {
                    h3F_dphi_vs_deta_vs_pt_RecoilCorrelation->Fill(PhiDistance, EtaDistance, particle.pt());
                    if(particle.pt() <= 3.0 && Distance <= DEF_RecoilCorrelationFrameSize && Distance > DEF_Epsilon)
                    {
                        NoOfParticlesInRecoilAreaCounter++;
                        h2F_RadialCorrelationRecoil->Fill(Distance, particle.pt(), 1.0/(2*Pi*Distance));
                    }
                }
            }

            //Fill reference correlation histogram
            TH3F* h3F_dphi_vs_deta_vs_pt_ReferenceCorrelation = new TH3F("h3F_dphi_vs_deta_vs_pt_ReferenceCorrelation", "h3F_dphi_vs_deta_vs_pt_ReferenceCorrelation",
                DEF_BinningPerUnit * DEF_RecoilCorrelationFrameSize / 2, -DEF_RecoilCorrelationFrameSize, DEF_RecoilCorrelationFrameSize,
                DEF_BinningPerUnit * DEF_RecoilCorrelationFrameSize / 2, -DEF_RecoilCorrelationFrameSize, DEF_RecoilCorrelationFrameSize,
                 0.1 * DEF_BinningPerUnit * DEF_MaxParticlePtInCorrelation, 0, DEF_MaxParticlePtInCorrelation);
            TH2F* h2F_RadialCorrelationReference = new TH2F("h2F_RadialCorrelationReference", "h2F_RadialCorrelationReference", 25, 0, DEF_RecoilCorrelationFrameSize, 30, 0, 3);

            for(const auto& particle : ParticleVector[EventProperties[iEvent].SuitableReferenceEvent])
            {
                float EtaDistance = particle.eta() - RecoilSiteEta;
                float PhiDistance = particle.phi_std() - RecoilSitePhi;
                if(PhiDistance > Pi) PhiDistance -= 2*Pi;
                else if(PhiDistance < -Pi) PhiDistance += 2*Pi;
                float Distance = sqrt(powf(EtaDistance, 2) + powf(PhiDistance, 2));
                if(abs(EtaDistance) <= DEF_RecoilCorrelationFrameSize && abs(PhiDistance) <= DEF_RecoilCorrelationFrameSize)
                {
                    h3F_dphi_vs_deta_vs_pt_ReferenceCorrelation->Fill(PhiDistance, EtaDistance, particle.pt());
                    if(particle.pt() <= 3.0 && Distance <= DEF_RecoilCorrelationFrameSize && Distance > DEF_Epsilon)
                    {
                        NoOfParticlesInReferenceAreaCounter++;
                        h2F_RadialCorrelationReference->Fill(Distance, particle.pt(), 1.0/(2*Pi*Distance));
                    } 
                }
            }

            h1F_NoOfParticlesInRecoilArea->Fill(NoOfParticlesInRecoilAreaCounter);
            h1F_NoOfParticlesInReferenceArea->Fill(NoOfParticlesInReferenceAreaCounter);

            //subtract the reference from the recoil correlation
            h3F_dphi_vs_deta_vs_pt_RecoilCorrelation->Add(h3F_dphi_vs_deta_vs_pt_ReferenceCorrelation, -1);
            h3F_dphi_vs_deta_vs_pt_CorrelationDifference->Add(h3F_dphi_vs_deta_vs_pt_RecoilCorrelation);

            h2F_RadialCorrelationDifference->Add(h2F_RadialCorrelationRecoil, +1);
            h2F_RadialCorrelationDifference->Add(h2F_RadialCorrelationReference, -1);

            delete(h3F_dphi_vs_deta_vs_pt_RecoilCorrelation);
            delete(h3F_dphi_vs_deta_vs_pt_ReferenceCorrelation);
            delete(h2F_RadialCorrelationRecoil);
            delete(h2F_RadialCorrelationReference);

            AnalyzedEventsCounter++;
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

    h3F_dphi_vs_deta_vs_pt_CorrelationDifference->Scale(1./(float)AnalyzedEventsCounter);
    h3F_dphi_vs_deta_vs_pt_CorrelationDifference->Scale(1./(h3F_dphi_vs_deta_vs_pt_CorrelationDifference->GetXaxis()->GetBinWidth(1) * h3F_dphi_vs_deta_vs_pt_CorrelationDifference->GetYaxis()->GetBinWidth(1)));
    h2F_RadialCorrelationDifference->Scale(1./(float)AnalyzedEventsCounter);
    h2F_RadialCorrelationDifference->Scale(1./(h2F_RadialCorrelationDifference->GetXaxis()->GetBinWidth(1) * h2F_RadialCorrelationDifference->GetYaxis()->GetBinWidth(1)));

    h3F_dphi_vs_deta_vs_pt_CorrelationDifference->SetTitle("Correlation Difference at the Recoil Site");
    h3F_dphi_vs_deta_vs_pt_CorrelationDifference->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
    h3F_dphi_vs_deta_vs_pt_CorrelationDifference->GetXaxis()->SetTitle("#phi^{Track} - #phi^{Leading Jet} [rad]");
    h3F_dphi_vs_deta_vs_pt_CorrelationDifference->GetYaxis()->SetTitleSize(DEF_AxisLabelSize);
    h3F_dphi_vs_deta_vs_pt_CorrelationDifference->GetYaxis()->SetTitle("#eta^{Track} - #eta^{Leading Jet}");
    h3F_dphi_vs_deta_vs_pt_CorrelationDifference->GetZaxis()->SetTitleSize(DEF_AxisLabelSize);
    h3F_dphi_vs_deta_vs_pt_CorrelationDifference->GetZaxis()->SetTitle("p_T^{Track}[GeV]");

    //Generate output file
    TString RunNumberAsString;
    RunNumberAsString.Form("%d", RunNumber);
    TString str_out =  "./" + (TString)RunNumberAsString.Data() + "_SingleJetAna.root";

    TFile* OutputFile = TFile::Open(str_out);
    if(!OutputFile)
    {
        OutputFile = new TFile(str_out.Data(),"RECREATE");
        OutputFile->cd();
        cout << "Output file is " << str_out << endl;
        h3F_dphi_vs_deta_vs_pt_CorrelationDifference->Write();
        h1F_NoOfParticlesInRecoilArea->Write(),
        h1F_NoOfParticlesInReferenceArea->Write();
        h1F_N_ev_minus_N_ref->Write();
        h1F_NoOfHighPtParticles->Write();
        h2F_RadialCorrelationDifference->Write();

        OutputFile->Close();
    }
    else
    {
        cout << "ADD TO EXISTING OUTPUT FILE" << endl;
        OutputFile->cd();

        TH3F* Summed_h3F_dphi_vs_deta_vs_pt_CorrelationDifference = (TH3F*) OutputFile->Get("h3F_dphi_vs_deta_vs_pt_CorrelationDifference");
        Summed_h3F_dphi_vs_deta_vs_pt_CorrelationDifference->Add(h3F_dphi_vs_deta_vs_pt_CorrelationDifference);

        TH1F* Summed_h1F_NoOfParticlesInRecoilArea = (TH1F*) OutputFile->Get("h1F_NoOfParticlesInRecoilArea");
        Summed_h1F_NoOfParticlesInRecoilArea->Add(h1F_NoOfParticlesInRecoilArea);

        TH1F* Summed_h1F_NoOfParticlesInReferenceArea = (TH1F*) OutputFile->Get("h1F_NoOfParticlesInReferenceArea");
        Summed_h1F_NoOfParticlesInReferenceArea->Add(h1F_NoOfParticlesInReferenceArea);

        TH1F* Summed_h1F_N_ev_minus_N_ref = (TH1F*) OutputFile->Get("h1F_N_ev_minus_N_ref");
        Summed_h1F_N_ev_minus_N_ref->Add(h1F_N_ev_minus_N_ref);

        TH1F* Summed_h1F_NoOfHighPtParticles = (TH1F*) OutputFile->Get("h1F_NoOfHighPtParticles");
        Summed_h1F_NoOfHighPtParticles->Add(h1F_NoOfHighPtParticles);

        TH1F* Summed_h2F_RadialCorrelationDifference = (TH1F*) OutputFile->Get("h2F_RadialCorrelationDifference");
        Summed_h2F_RadialCorrelationDifference->Add(h2F_RadialCorrelationDifference);


        OutputFile = new TFile(str_out.Data(),"RECREATE");
        Summed_h3F_dphi_vs_deta_vs_pt_CorrelationDifference->Write();
        Summed_h1F_NoOfParticlesInRecoilArea->Write();
        Summed_h1F_NoOfParticlesInReferenceArea->Write();
        Summed_h1F_N_ev_minus_N_ref->Write();
        Summed_h1F_NoOfHighPtParticles->Write();
        Summed_h2F_RadialCorrelationDifference->Write();

        OutputFile->Close();

    }

    return 1;
}