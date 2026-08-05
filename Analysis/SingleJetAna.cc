#include "./SingleJetAna.h"

class MyUserInfo : public fastjet::PseudoJet::UserInfoBase {
public:
    MyUserInfo(int NoOfParticles) 
        : m_NoOfParticles(NoOfParticles){}

    int getNoOfParticles() const { return m_NoOfParticles; }

private:
    int m_NoOfParticles;
};

/* Parameters:
i_InputFile: Format xxx.root
i_OutputFile: Format xxx
*/
Int_t SingleJetAna(TString i_InputFile = "in.root", TString i_DFName = "DF_xxx", int i_NoOfEvents = -1)
{
    //CONFIGURE
    bool i_UseHadronInstead = false;

    gStyle->SetOptStat(0);
    SetRootGraphicStyle();

    #define DEF_BinningPerUnit 100
    #define DEF_JetRadius 0.2
    #define DEF_OutputEventOverviews true 
    #define DEF_JetLeadingPt 30
    #define DEF_HadLeadingPt 10.0
    #define DEF_HadSubleadingPt 7.0
    #define DEF_CorrelationMinPt 1.0
    #define DEF_CorrelationMaxPt 2.0
    #define DEF_MaxPartilclesPerJet 200
    #define DEF_MaxJetPt 150.0
    #define DEF_AxisLabelSize 0.05
    #define DEF_HistoTitleSize 0.1
    #define DEF_EventTolerance_Multiplicity 50
    #define DEF_EventTolerance_Psi2 Pi/16

    int TotalEvents = 0;
    int ProcessedEvents = 0;
    int EventsSkipped   = 0;

    //set jet/had pt limit
    double LeadingPtLimit = DEF_JetLeadingPt;
    if(i_UseHadronInstead)
    {
        LeadingPtLimit = DEF_HadLeadingPt;
    }

    Long64_t LargeGapEventCounter = 0;
    Long64_t SmallGapEventCounter = 0;

    /*END USER VARIABLES*/

    //Generate output file
    TString str_out = i_DFName;
    str_out += ".root";
    TFile* OutputFile = new TFile(str_out.Data(),"RECREATE");
    OutputFile->cd();
    cout << "Output file is " << str_out << endl;
    TDirectory *EventHistos = OutputFile->mkdir("Event Histos");
    TDirectory *Results = OutputFile->mkdir("Results");

    //Get data from tree file
    TFile* file = TFile::Open(i_InputFile);
    if(!file || file->IsZombie())
    {
        cout << "Input File could not be opened. Wrong path?" << endl;
        return 0;
    }

    TTree * collisions = nullptr;
    TTree * tracks = nullptr;

    // loop over all directories and print name
    cout << "Process file" << endl;
    TIter next(file->GetListOfKeys());
    TKey *key;
    while ((key = (TKey *)next())) {
        TString DFName = key->GetName();
        if(DFName != i_DFName) continue;
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

        // Make an unordered map particle -> collision
        for(uint64_t i_Track = 0; i_Track < entries_track; i_Track ++)
        {
            Particle track = Particle::Read(tracks, i_Track);
            //cout << "track " << i_Track << ", col ID " << track.ColID << endl;
            CollisionMap[track.ColID].push_back(i_Track);
            //prepare particles for jetfinder
            ParticleVector[track.ColID].push_back(PseudoJet(track.px, track.py, track.pz, track.E));

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

        //APPLICATION VARIABLES
        struct struct_EventProperties
        {
            int Multiplicity = 0;
            float Psi2 = -100;
            vector<PseudoJet> JetCoordinates;
            bool IsHighPtJetEvent = false;
            int SuitableReferenceEvent = -1;
            
        };
        vector<struct_EventProperties> EventProperties;
        int HighPtEventCounter = 0;
        int PairedEventsCounter = 0;

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

        if(i_NoOfEvents != -1) entries_col = i_NoOfEvents;
        for (int i_Collision = 0; i_Collision < entries_col; i_Collision++) {

            Collision collision = Collision::Read(collisions, i_Collision);
            if(CollisionMap[collision.ColID].size() == 0 || ParticleVector[collision.ColID].size() == 0){
                cout << "No tracks in collision " << collision.ColID << "! Skipping." << endl;
                N_events_skipped++;
                continue;
            }

            struct_EventProperties CurrentEventProperties;

            CurrentEventProperties.Multiplicity = CollisionMap[collision.ColID].size();
            CurrentEventProperties.Psi2 = collision.Psi2;

            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            /*
            /
            /                                                     
            /                                           Find Jets with the anti-k_t algorithm  
            /                                                     
            /                                                     
            /                                                     
            /
            /                                                    
            /
            /
            /
            *//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            vector<PseudoJet> Jets = FindJets(ParticleVector[collision.ColID], DEF_JetLeadingPt);

            CurrentEventProperties.JetCoordinates = Jets;

            //output event overviews for debug reasons
            if(DEF_OutputEventOverviews)
            {
                TH2D* h2D_eta_vs_phi = new TH2D("h2D_eta_vs_phi","h2D_eta_vs_phi",DEF_BinningPerUnit * 1.8,-0.9,0.9,DEF_BinningPerUnit * 2 * Pi,-Pi,Pi);//temporarily filled

                for(int iParticle = 0; iParticle < ParticleVector[collision.ColID].size(); iParticle++)
                {
                    h2D_eta_vs_phi->Fill(ParticleVector[collision.ColID][iParticle].eta(), ParticleVector[collision.ColID][iParticle].phi_std(), ParticleVector[collision.ColID][iParticle].pt());
                }

                EventHistos->cd();
                h2D_eta_vs_phi->Write();
                delete(h2D_eta_vs_phi);

            }

            cout << "Jets in Event: " << Jets.size() << endl;

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

            //push back, because if event is not suitable, the event is in the list nevertheless for potential background subtraction
            EventProperties.push_back(CurrentEventProperties);

            if(Jets.size() < 1) continue;

            if(Jets[0].pt() < DEF_JetLeadingPt) continue;

            

            //check that there is no jet close to the direct recoil site
            float RecoilSiteEta = Jets[0].eta();
            float RecoilSitePhi;
            if(Jets[0].phi_std() >= 0)
            {
                RecoilSitePhi = Jets[0].phi_std() - Pi;
            }
            else
            {
                RecoilSitePhi = Jets[0].phi_std() + Pi;
            }
            for(int iJet = 0; iJet < Jets.size(); iJet++)
            {
                float EtaDistance = Jets[iJet].eta() - RecoilSiteEta;
                float PhiDistance = Jets[iJet].phi_std() - RecoilSitePhi;
                if(PhiDistance > Pi) PhiDistance -= 2*Pi;
                else if(PhiDistance < -Pi) PhiDistance += 2*Pi;
                float Distance = sqrt(pow(EtaDistance, 2) + pow(PhiDistance, 2));
                if(Distance < 2*DEF_JetRadius) continue;
            }

            CurrentEventProperties.IsHighPtJetEvent = true;
            HighPtEventCounter++;
            EventProperties.pop_back();
            EventProperties.push_back(CurrentEventProperties);

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

        for(int iEvent = 0; iEvent < EventProperties.size(); iEvent++)
        {
            if(!EventProperties[iEvent].IsHighPtJetEvent) continue;

            for(int iReferenceEvent = 0; iReferenceEvent < EventProperties.size(); iReferenceEvent++)
            {
                bool BadReferenceEvent = false;

                if(iReferenceEvent == iEvent) continue;

                if(abs(EventProperties[iReferenceEvent].Multiplicity - EventProperties[iEvent].Multiplicity) > DEF_EventTolerance_Multiplicity) continue;

                if(abs(EventProperties[iReferenceEvent].Psi2 - EventProperties[iEvent].Psi2) > DEF_EventTolerance_Psi2) continue;

                //check that there is no jet in the reference event where the reference data is taken
                float RecoilSiteEta = EventProperties[iEvent].JetCoordinates[0].eta();
                float RecoilSitePhi;
                if(EventProperties[iEvent].JetCoordinates[0].phi_std() >= 0)
                {
                    RecoilSitePhi = EventProperties[iEvent].JetCoordinates[0].phi_std() - Pi;
                }
                else
                {
                    RecoilSitePhi = EventProperties[iEvent].JetCoordinates[0].phi_std() + Pi;
                }
                for(int iJet = 0; iJet < EventProperties[iReferenceEvent].JetCoordinates.size(); iJet++)
                {
                    float EtaDistance = EventProperties[iReferenceEvent].JetCoordinates[iJet].eta() - RecoilSiteEta;
                    float PhiDistance = EventProperties[iReferenceEvent].JetCoordinates[iJet].phi_std() - RecoilSitePhi;
                    if(PhiDistance > Pi) PhiDistance -= 2*Pi;
                    else if(PhiDistance < -Pi) PhiDistance += 2*Pi;
                    float Distance = sqrt(pow(EtaDistance, 2) + pow(PhiDistance, 2));
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
        /                                     Count how many high pt events have got a suitable reference event  
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
        cout << "Skipped because of unusual event behaviour: " << N_events_skipped << "/" << entries_col << " events" << endl;
        cout << "In the end " << ((float)(SmallGapEventCounter + LargeGapEventCounter))/((float)entries_col) * 100 << " percent of events have been accepted!" << endl;
        cout << "Total number of events: " << TotalEvents << ", of which skipped " << EventsSkipped << endl;
        
        EventsSkipped += N_events_skipped;

    }

    OutputFile->Close();

    return 1;
}



vector<PseudoJet> FindJets(const vector<PseudoJet> vec_particles, float i_LowPtLimit)
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

    // To do: define jet radius, background radius etc
    #define jet_radius 0.2
    // define variables 
    vector<PseudoJet> jets_fiducial;//contains the selected jets in the end

    Double_t ghost_maxrap = 1.1; // Fiducial cut for background estimation //1.1
    Double_t eta_acceptance = 0.9;
    // choose a jet definition
    JetDefinition jet_def(antikt_algorithm, jet_radius);
    // jet area definition
    GhostedAreaSpec area_spec(ghost_maxrap);
    AreaDefinition area_def(active_area_explicit_ghosts,GhostedAreaSpec(ghost_maxrap,1,0.01));

    ClusterSequenceArea clust_seq_hard(vec_particles, jet_def, area_def);
    // run the clustering, extract the jets
    double ptmin = 0.15;
    vector<PseudoJet> jets_all = sorted_by_pt(clust_seq_hard.inclusive_jets(ptmin));
    Selector Fiducial_cut_selector = SelectorAbsEtaMax(eta_acceptance - jet_radius); // Fiducial cut for jets
    jets_fiducial = Fiducial_cut_selector(jets_all);

    // background estimation
    // Standard rho calculation with kT algorithm and N hardest jets removed

    Int_t Rem_n_hardest = 2; 
    double eBkg_R = 0.2;
    JetDefinition jet_def_bkgd(kt_algorithm, eBkg_R); // <--
    AreaDefinition area_def_bkgd(active_area_explicit_ghosts,GhostedAreaSpec(ghost_maxrap,1,0.01));
    Selector selector = SelectorAbsEtaMax(0.9 - eBkg_R) * (!SelectorNHardest(Rem_n_hardest)); // <--

    JetMedianBackgroundEstimator bkgd_estimator(selector, jet_def_bkgd, area_def_bkgd); // <--
    Subtractor subtractor(&bkgd_estimator);
    bkgd_estimator.set_particles(vec_particles);
    Double_t jet_rho   = bkgd_estimator.rho();

    vector<PseudoJet> AcceptedJets;

    // loop over jets  
    Double_t Jet_area_cut = 0.56*TMath::Pi()*TMath::Power(jet_radius,2);
    for(Int_t i_jet = 0; i_jet < (Int_t)jets_fiducial.size(); i_jet++)
    {
        Float_t jet_pt         = jets_fiducial[i_jet].perp();
        if(jet_pt <= 0.0) continue;
        Float_t jet_area       = jets_fiducial[i_jet].area();
        Float_t jet_pt_sub     = jets_fiducial[i_jet].perp() - jet_rho*jet_area;
        if(jet_pt_sub < i_LowPtLimit) continue;
        Float_t jet_eta        = jets_fiducial[i_jet].eta();
        Float_t jet_phi        = jets_fiducial[i_jet].phi();
        Float_t jet_phi_std    = jets_fiducial[i_jet].phi_std();
        //if(jet_area < Jet_area_cut) continue;  this is a problem apparently

        // Store jet size in user information
        jets_fiducial[i_jet].set_user_info(new MyUserInfo(jets_fiducial[i_jet].constituents().size()));

        AcceptedJets.push_back(jets_fiducial[i_jet]);

        /*
        // constituent track loop
        vector<PseudoJet> jet_constituents = AcceptedJets[i_jet].constituents();

        for(Int_t i_constituent = 0; i_constituent < (Int_t)jet_constituents.size(); i_constituent++)
        {
            Float_t jet_const_pt  = jet_constituents[i_constituent].perp();
            if(jet_const_pt <= 0.0) continue;
            Float_t jet_const_phi = jet_constituents[i_constituent].phi();
            Float_t jet_const_rap = jet_constituents[i_constituent].rap();
            Float_t jet_const_Phi = jet_constituents[i_constituent].phi_std();
            Float_t jet_const_eta = jet_constituents[i_constituent].eta();
            Int_t   user_index    = jet_constituents[i_constituent].user_index();

        }
        */
    }

    return sorted_by_pt(AcceptedJets);
}