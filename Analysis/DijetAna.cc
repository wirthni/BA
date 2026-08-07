#include "./DijetAna.h"

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
Int_t DijetAna(TString i_InputFile = "in.root", TString i_DFName = "DF_xxx")
{
    //CONFIGURE
    bool i_UseHadronInstead = false;

    gStyle->SetOptStat(0);
    SetRootGraphicStyle();

    #define DEF_BinningPerUnit 100
    #define DEF_JetRadius 0.2
    #define DEF_OutputEventOverviews false 
    #define DEF_JetLeadingPt 40
    #define DEF_JetSubleadingPt 20
    #define DEF_HadLeadingPt 10.0
    #define DEF_HadSubleadingPt 7.0
    #define DEF_CorrelationMinPt 1.0
    #define DEF_CorrelationMaxPt 2.0
    #define DEF_MaxPartilclesPerJet 200
    #define DEF_MaxJetPt 150.0
    #define DEF_AxisLabelSize 0.05
    #define DEF_HistoTitleSize 0.1

    TH2D* h2D_pt_vs_eta_LargeGap = new TH2D("h2D_pt_vs_eta_LargeGap","h2D_pt_vs_eta_LargeGap",100, 0, 10, DEF_BinningPerUnit * 2 * 0.9,-0.9,0.9);
    TH2D* h2D_pt_vs_eta_SmallGap = new TH2D("h2D_pt_vs_eta_SmallGap","h2D_pt_vs_eta_SmallGap",100, 0, 10, DEF_BinningPerUnit * 2 * 0.9,-0.9,0.9);
    TH2D* h2D_eta_vs_dphi_LargeGap = new TH2D("h2D_eta_vs_dphi_LargeGap","h2D_eta_vs_dphi_LargeGap",DEF_BinningPerUnit * 2 * 0.9, -0.9, 0.9, DEF_BinningPerUnit * 2 * Pi, -Pi/2, 3*Pi/2);
    TH2D* h2D_eta_vs_dphi_SmallGap = new TH2D("h2D_eta_vs_dphi_SmallGap","h2D_eta_vs_dphi_SmallGap",DEF_BinningPerUnit * 2 * 0.9, -0.9, 0.9, DEF_BinningPerUnit * 2 * Pi, -Pi/2, 3*Pi/2);
    TH1D* h1D_JetPopulation_Leading = new TH1D("h1D_JetPopulation_Leading", "h1D_JetPopulation_Leading", DEF_MaxPartilclesPerJet, 0, DEF_MaxPartilclesPerJet);
    TH1D* h1D_JetPopulation_Subleading = new TH1D("h1D_JetPopulation_Subleading", "h1D_JetPopulation_Subleading", DEF_MaxPartilclesPerJet, 0, DEF_MaxPartilclesPerJet);
    TH1D* h1D_JetPt_Leading = new TH1D("h1D_JetPt_Leading", "h1D_JetPt_Leading", DEF_MaxJetPt, 0, DEF_MaxJetPt);
    TH1D* h1D_JetPt_Subleading = new TH1D("h1D_JetPt_Subleading", "h1D_JetPt_Subleading", DEF_MaxJetPt, 0, DEF_MaxJetPt);
    TH2D* h2D_eta_vs_phi_JetCoordinates_Leading_NoCut = new TH2D("h2D_eta_vs_phi_JetCoordinates_Leading_NoCut","h2D_eta_vs_phi_JetCoordinates_Leading_NoCut",DEF_BinningPerUnit * 2 * 0.9, -0.9, 0.9, DEF_BinningPerUnit * 2 * Pi, -Pi, Pi); // contains only the jet coordinates before applying the cut
    TH2D* h2D_ParticleCorrelation = new TH2D("h2D_ParticleCorrelation","h2D_ParticleCorrelation",DEF_BinningPerUnit * 2 * 0.9, -0.9, 0.9, DEF_BinningPerUnit * 2 * Pi, -Pi, Pi); // Contains all particles

    int TotalEvents = 0;
    int ProcessedEvents = 0;
    int EventsSkipped   = 0;

    //set jet/had pt limit
    double LeadingPtLimit = DEF_JetLeadingPt;
    double SubleadingPtLimit = DEF_JetSubleadingPt;
    if(i_UseHadronInstead)
    {
        LeadingPtLimit = DEF_HadLeadingPt;
        SubleadingPtLimit = DEF_HadSubleadingPt;
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

        for (int i_Collision = 0; i_Collision < entries_col; i_Collision++) {

            ProcessedEvents++;

            Collision collision = Collision::Read(collisions, i_Collision);
            if(CollisionMap[collision.ColID].size() == 0 || ParticleVector[collision.ColID].size() == 0){
                cout << "No tracks in collision " << collision.ColID << "! Skipping." << endl;
                N_events_skipped++;
                continue;
            }

            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            /*
            /
            /                                                     
            /                                           FIND JETS WITH THE ANTI-KT ALGORITHM OR FIND HIGH PT HADRONS   
            /                                                     
            /                                                     
            /                                                     
            /
            /                                                    
            /
            /
            /
            *//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            vector<PseudoJet> JetVector;

            if(i_UseHadronInstead)
            {
                vector<PseudoJet> ParticleVectorCopy = ParticleVector[collision.ColID];
                int NoOfFoundJets = 0;
                do{

                    Int_t HighestPtIndex = 0;

                    //search for highest momentum and choose corresponding particle as a primary vertex
                    for(Int_t i=0; i<ParticleVectorCopy.size(); i++)
                    {
                        if(ParticleVectorCopy[i].pt() > ParticleVectorCopy[HighestPtIndex].pt())
                        {
                            HighestPtIndex = i;
                        }
                    }

                    //push particle into jet vector
                    JetVector.push_back(ParticleVectorCopy[HighestPtIndex]);

                    //pop jet particle
                    ParticleVectorCopy.erase(ParticleVectorCopy.begin() + HighestPtIndex);

                    //look for other high pt particles in a radius of DEF_JetRadius because they belong to the same jet
                    for(Int_t i=0; i<ParticleVectorCopy.size(); i++)
                    {

                        double DeltaEta = abs(JetVector.back().eta() - ParticleVectorCopy[i].eta());
                        double DeltaPhi = abs(JetVector.back().phi_std() - ParticleVectorCopy[i].phi_std());

                        //since a jet at pi and another one at -pi are equivalent, normalize
                        if(DeltaPhi > Pi) DeltaPhi -= 2* Pi;

                        double DeltaR = TMath::Sqrt(DeltaEta*DeltaEta + DeltaPhi*DeltaPhi);

                        if(DeltaR <= DEF_JetRadius)
                        {
                            //pop candidate
                            ParticleVectorCopy.erase(ParticleVectorCopy.begin()+i);
                            i = -1;//start again for safe
                        }
                    }

                }while(ParticleVectorCopy.size() > 0 && NoOfFoundJets < 2);

                //sort high pt hadrons by size to have same structure as FindJets return type
                JetVector = sorted_by_pt(JetVector);
                //make eta cut
                Selector Fiducial_cut_selector = SelectorAbsEtaMax(0.9 - DEF_JetRadius); // Fiducial cut for jets
                JetVector = Fiducial_cut_selector(JetVector);

            }
            else
            {
                JetVector = FindJets(ParticleVector[collision.ColID]);
            }

            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            /*
            /
            /                                                     
            /                                           FILTER EVENTS BASED ON THE CONSTRAINTS AND DIFFER BETW. SMALL/LARGE GAP
            /                       -DIJET EVENT
            /                       -LEADING JET PT > DEF_LeadingPt, SUBLEADING JET PT > DEF_SubleadingPt
            /                       -Eta_Jet1 > 0
            /                       -Delta Phi(Jet1, Jet2) > pi/2
            /                       -(LARGE GAP / SMALL GAP) Eta_Jet1 * Eta_Jet2 (</>) 0        
            /                                                     
            /                                                     
            /
            /                                                    
            /
            /
            /
            *//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            if(JetVector.size() >= 1) h2D_eta_vs_phi_JetCoordinates_Leading_NoCut->Fill(JetVector[0].eta(), JetVector[0].phi_std());

            if(JetVector.size() < 2 )continue;
            
            if(JetVector[0].eta() < 0) continue;
            
            if(JetVector[0].pt() <= LeadingPtLimit || JetVector[1].pt() <= SubleadingPtLimit) continue;
            
            double PhiSeparation = abs(JetVector[0].phi_std() - JetVector[1].phi_std());
            if(PhiSeparation > Pi) PhiSeparation -= 2*Pi;
            if(abs(PhiSeparation) < Pi/2) continue;

            bool LargeGap = false;
            if(JetVector[0].eta() * JetVector[1].eta() < 0) LargeGap = true;

            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            /*
            /
            /                                                     
            /                                           GET PARTICLE MULTIPLICITES IN THE REGIONS OF INTEREST
            /                       
            /                       
            /                       
            /                       
            /                             
            /                                                     
            /                                                     
            /
            /                                                    
            /
            /
            /
            *//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            
            #if DEF_OutputEventOverviews
            TH2D* h2D_phi_vs_eta = new TH2D("h2D_phi_vs_eta","h2D_phi_vs_eta", DEF_BinningPerUnit * 2 * Pi, -Pi, Pi, DEF_BinningPerUnit * 2 * 0.9, -0.9, 0.9);
            #endif

            // Fill leading and subleading jet pt
            h1D_JetPt_Leading->Fill(JetVector[0].pt());
            h1D_JetPt_Subleading->Fill(JetVector[1].pt());

            // Prepare / Fill histogram on how many particles contribute to each jet
            int LeadingJetParticleCounter = 0;
            int SubleadingJetParticleCounter = 0;

            if(!i_UseHadronInstead)
            {
                h1D_JetPopulation_Leading->Fill(JetVector[0].user_info<MyUserInfo>().getNoOfParticles());
                h1D_JetPopulation_Subleading->Fill(JetVector[1].user_info<MyUserInfo>().getNoOfParticles());
            }

            if(LargeGap)
            {
                
                for(const auto& particle : ParticleVector[collision.ColID])
                {
                    //fill particle correlation for this event
                    h2D_ParticleCorrelation->Fill(particle.eta(), particle.phi_std());

                    // Fill analysis plot
                    double DeltaPhi = particle.phi_std() - JetVector[0].phi_std();
                    if(DeltaPhi > Pi) DeltaPhi -= 2*Pi;
                    else if(DeltaPhi < -Pi) DeltaPhi += 2*Pi;
                    if(abs(DeltaPhi) <= Pi/2) h2D_pt_vs_eta_LargeGap->Fill(particle.pt(), particle.eta());

                    // Fill event overview
                    #if DEF_OutputEventOverviews
                        h2D_phi_vs_eta->Fill(particle.phi_std(), particle.eta(), particle.pt());
                    #endif

                    // Fill correlation relative to leading jet in phi
                    if(DeltaPhi < -Pi/2) DeltaPhi += 2*Pi;
                    else if(DeltaPhi > 3*Pi/2) DeltaPhi -= 2*Pi;
                    if(particle.pt() >= DEF_CorrelationMinPt && particle.pt() <= DEF_CorrelationMaxPt)h2D_eta_vs_dphi_LargeGap->Fill(particle.eta(), DeltaPhi);

                    //in case of hadron analysis, look how many particles contribute to the leading jet
                    if(i_UseHadronInstead)
                    {
                        double DeltaEta;
                        //leading jet
                        DeltaPhi = particle.phi_std() - JetVector[0].phi_std();
                        if(DeltaPhi > Pi) DeltaPhi -= 2*Pi;
                        else if(DeltaPhi < -Pi) DeltaPhi += 2*Pi;
                        DeltaEta = particle.eta() - JetVector[0].eta();
                        if(sqrt(pow(DeltaPhi, 2) + pow(DeltaEta, 2)) <= DEF_JetRadius) LeadingJetParticleCounter++;

                        //subleading jet
                        DeltaPhi = particle.phi_std() - JetVector[1].phi_std();
                        if(DeltaPhi > Pi) DeltaPhi -= 2*Pi;
                        else if(DeltaPhi < -Pi) DeltaPhi += 2*Pi;
                        DeltaEta = particle.eta() - JetVector[1].eta();
                        if(sqrt(pow(DeltaPhi, 2) + pow(DeltaEta, 2)) <= DEF_JetRadius) SubleadingJetParticleCounter++;
                    }

                }
                LargeGapEventCounter++;
            }
            else
            {
                for(const auto& particle : ParticleVector[collision.ColID])
                {
                    //fill particle correlation for this event
                    h2D_ParticleCorrelation->Fill(particle.eta(), particle.phi_std());
                    
                    // Fill analysis plot
                    double DeltaPhi = particle.phi_std() - JetVector[0].phi_std();
                    if(DeltaPhi > Pi) DeltaPhi -= 2*Pi;
                    else if(DeltaPhi < -Pi) DeltaPhi += 2*Pi;
                    if(abs(DeltaPhi) <= Pi/2) h2D_pt_vs_eta_SmallGap->Fill(particle.pt(), particle.eta());

                    // Fill event overview
                    #if DEF_OutputEventOverviews
                        h2D_phi_vs_eta->Fill(particle.phi_std(), particle.eta(), particle.pt());
                    #endif

                    // Fill correlation relative to leading jet in phi
                    if(DeltaPhi < -Pi/2) DeltaPhi += 2*Pi;
                    else if(DeltaPhi > 3*Pi/2) DeltaPhi -= 2*Pi;
                    if(particle.pt() >= DEF_CorrelationMinPt && particle.pt() <= DEF_CorrelationMaxPt)h2D_eta_vs_dphi_SmallGap->Fill(particle.eta(), DeltaPhi);

                    //in case of hadron analysis, look how many particles contribute to the leading jet
                    if(i_UseHadronInstead)
                    {
                        double DeltaEta;
                        //leading jet
                        DeltaPhi = particle.phi_std() - JetVector[0].phi_std();
                        if(DeltaPhi > Pi) DeltaPhi -= 2*Pi;
                        else if(DeltaPhi < -Pi) DeltaPhi += 2*Pi;
                        DeltaEta = particle.eta() - JetVector[0].eta();
                        if(sqrt(pow(DeltaPhi, 2) + pow(DeltaEta, 2)) <= DEF_JetRadius) LeadingJetParticleCounter++;

                        //subleading jet
                        DeltaPhi = particle.phi_std() - JetVector[1].phi_std();
                        if(DeltaPhi > Pi) DeltaPhi -= 2*Pi;
                        else if(DeltaPhi < -Pi) DeltaPhi += 2*Pi;
                        DeltaEta = particle.eta() - JetVector[1].eta();
                        if(sqrt(pow(DeltaPhi, 2) + pow(DeltaEta, 2)) <= DEF_JetRadius) SubleadingJetParticleCounter++;
                    }

                }
                SmallGapEventCounter++;
            }

            if(i_UseHadronInstead)
            {
                h1D_JetPopulation_Leading->Fill(LeadingJetParticleCounter);
                h1D_JetPopulation_Subleading->Fill(SubleadingJetParticleCounter);
            }

            //write event results
            #if DEF_OutputEventOverviews
                EventHistos->cd();
                h2D_phi_vs_eta->Write();
                delete h2D_phi_vs_eta;
            #endif


            ////////////////////////////////////////////////////////////////////////////

        }

        cout << "Skipped because of unusual event behaviour: " << N_events_skipped << "/" << entries_col << " events" << endl;
        cout << "In the end " << ((float)(SmallGapEventCounter + LargeGapEventCounter))/((float)entries_col) * 100 << " percent of events have been accepted!" << endl;
        cout << "Total number of events: " << TotalEvents << ", of which skipped " << EventsSkipped << endl;
        
        EventsSkipped += N_events_skipped;

    }

    // Normalize histograms by number of events
    h2D_pt_vs_eta_LargeGap->Scale(1./(double)LargeGapEventCounter);
    h2D_pt_vs_eta_SmallGap->Scale(1./(double)SmallGapEventCounter);

    h2D_eta_vs_dphi_LargeGap->Scale(1./(double)LargeGapEventCounter);
    h2D_eta_vs_dphi_SmallGap->Scale(1./(double)SmallGapEventCounter);

    h1D_JetPopulation_Leading->Scale(1./(double)(LargeGapEventCounter + SmallGapEventCounter));
    h1D_JetPopulation_Subleading->Scale(1./(double)(LargeGapEventCounter + SmallGapEventCounter));

    h1D_JetPt_Leading->Scale(1./(double)(LargeGapEventCounter + SmallGapEventCounter));
    h1D_JetPt_Subleading->Scale(1./(double)(LargeGapEventCounter + SmallGapEventCounter));

    // Rebin Histograms if necessary
    h2D_eta_vs_dphi_LargeGap->Rebin2D(8,8);
    h2D_eta_vs_dphi_SmallGap->Rebin2D(8,8);
    h2D_eta_vs_dphi_LargeGap->Scale(1.0/64.0);
    h2D_eta_vs_dphi_SmallGap->Scale(1.0/64.0);

    // Divide by bin widths
    h2D_eta_vs_dphi_LargeGap->Scale(1./(h2D_eta_vs_dphi_LargeGap->GetXaxis()->GetBinWidth(0) * h2D_eta_vs_dphi_LargeGap->GetYaxis()->GetBinWidth(0)));
    h2D_eta_vs_dphi_SmallGap->Scale(1./(h2D_eta_vs_dphi_SmallGap->GetXaxis()->GetBinWidth(0) * h2D_eta_vs_dphi_SmallGap->GetYaxis()->GetBinWidth(0)));

    h2D_pt_vs_eta_LargeGap->Scale(1./(h2D_pt_vs_eta_LargeGap->GetXaxis()->GetBinWidth(0) * h2D_pt_vs_eta_LargeGap->GetYaxis()->GetBinWidth(0)));
    h2D_pt_vs_eta_SmallGap->Scale(1./(h2D_pt_vs_eta_SmallGap->GetXaxis()->GetBinWidth(0) * h2D_pt_vs_eta_SmallGap->GetYaxis()->GetBinWidth(0)));

    h1D_JetPopulation_Leading->Scale(1./h1D_JetPopulation_Leading->GetXaxis()->GetBinWidth(0));
    h1D_JetPopulation_Subleading->Scale(1./h1D_JetPopulation_Subleading->GetXaxis()->GetBinWidth(0));
    
    h1D_JetPt_Leading->Scale(1./h1D_JetPt_Leading->GetXaxis()->GetBinWidth(0));
    h1D_JetPt_Subleading->Scale(1./h1D_JetPt_Subleading->GetXaxis()->GetBinWidth(0));

    h2D_eta_vs_phi_JetCoordinates_Leading_NoCut->Scale(1./(h2D_eta_vs_phi_JetCoordinates_Leading_NoCut->GetXaxis()->GetBinWidth(0) * h2D_eta_vs_phi_JetCoordinates_Leading_NoCut->GetYaxis()->GetBinWidth(0)));

    //normalize histograms if necessary

    //put axis labels on all histograms
    h2D_pt_vs_eta_LargeGap->SetTitle("Near-side correlation of particles with p_{t}, large gaps");
    h2D_pt_vs_eta_LargeGap->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
    h2D_pt_vs_eta_LargeGap->GetXaxis()->SetTitle("p_{t}[GeV]");
    h2D_pt_vs_eta_LargeGap->GetYaxis()->SetTitleSize(DEF_AxisLabelSize);
    h2D_pt_vs_eta_LargeGap->GetYaxis()->SetTitle("#eta");
    h2D_pt_vs_eta_LargeGap->GetZaxis()->SetTitleSize(DEF_AxisLabelSize);
    h2D_pt_vs_eta_LargeGap->GetZaxis()->SetTitle("#frac{1}{N_{Events}} #frac{d N}{d p_{t} d #eta}");

    h2D_pt_vs_eta_SmallGap->SetTitle("Near-side correlation of particles with p_{t}, small gaps");
    h2D_pt_vs_eta_SmallGap->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
    h2D_pt_vs_eta_SmallGap->GetXaxis()->SetTitle("p_{t}[GeV]");
    h2D_pt_vs_eta_SmallGap->GetYaxis()->SetTitleSize(DEF_AxisLabelSize);
    h2D_pt_vs_eta_SmallGap->GetYaxis()->SetTitle("#eta");
    h2D_pt_vs_eta_SmallGap->GetZaxis()->SetTitleSize(DEF_AxisLabelSize);
    h2D_pt_vs_eta_SmallGap->GetZaxis()->SetTitle("#frac{1}{N_{Events}} #frac{d N}{d p_{t} d #eta}");

    h2D_eta_vs_dphi_LargeGap->SetTitle("Particle correlation abs. in #eta, rel. in #phi, large gaps. p_{t} #in [1,2] GeV");
    h2D_eta_vs_dphi_LargeGap->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
    h2D_eta_vs_dphi_LargeGap->GetXaxis()->SetTitle("#eta");
    h2D_eta_vs_dphi_LargeGap->GetYaxis()->SetTitleSize(DEF_AxisLabelSize);
    h2D_eta_vs_dphi_LargeGap->GetYaxis()->SetTitle("#Delta #phi[rad]");
    h2D_eta_vs_dphi_LargeGap->GetZaxis()->SetTitleSize(DEF_AxisLabelSize);
    h2D_eta_vs_dphi_LargeGap->GetZaxis()->SetTitle("#frac{1}{N_{Events}} #frac{d N}{d #phi d #eta}");

    h2D_eta_vs_dphi_SmallGap->SetTitle("Particle correlation abs. in #eta, rel. in #phi, small gaps. p_{t} #in [1,2] GeV");
    h2D_eta_vs_dphi_SmallGap->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
    h2D_eta_vs_dphi_SmallGap->GetXaxis()->SetTitle("#eta");
    h2D_eta_vs_dphi_SmallGap->GetYaxis()->SetTitleSize(DEF_AxisLabelSize);
    h2D_eta_vs_dphi_SmallGap->GetYaxis()->SetTitle("#Delta #phi[rad]");
    h2D_eta_vs_dphi_SmallGap->GetZaxis()->SetTitleSize(DEF_AxisLabelSize);
    h2D_eta_vs_dphi_SmallGap->GetZaxis()->SetTitle("#frac{1}{N_{Events}} #frac{d N}{d #phi d #eta}");

    h1D_JetPopulation_Leading->SetTitle("Population of the leading jet");
    h1D_JetPopulation_Leading->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
    h1D_JetPopulation_Leading->GetXaxis()->SetTitle("Leading Jet Population");
    h1D_JetPopulation_Leading->GetYaxis()->SetTitleSize(DEF_AxisLabelSize);
    h1D_JetPopulation_Leading->GetYaxis()->SetTitle("a. u.");

    h1D_JetPopulation_Subleading->SetTitle("Population of the subleading jet");
    h1D_JetPopulation_Subleading->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
    h1D_JetPopulation_Subleading->GetXaxis()->SetTitle("Subleading Jet Population");
    h1D_JetPopulation_Subleading->GetYaxis()->SetTitleSize(DEF_AxisLabelSize);
    h1D_JetPopulation_Subleading->GetYaxis()->SetTitle("a. u.");

    h1D_JetPt_Leading->SetTitle("Pt of the leading jet");
    h1D_JetPt_Leading->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
    h1D_JetPt_Leading->GetXaxis()->SetTitle("p_{t}^{leading jet} [GeV]");
    h1D_JetPt_Leading->GetYaxis()->SetTitleSize(DEF_AxisLabelSize);
    h1D_JetPt_Leading->GetYaxis()->SetTitle("a. u.");

    h1D_JetPt_Subleading->SetTitle("Pt of the subleading jet");
    h1D_JetPt_Subleading->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
    h1D_JetPt_Subleading->GetXaxis()->SetTitle("p_{t}^{subleading jet} [GeV]");
    h1D_JetPt_Subleading->GetYaxis()->SetTitleSize(DEF_AxisLabelSize);
    h1D_JetPt_Subleading->GetYaxis()->SetTitle("a. u.");

    h2D_eta_vs_phi_JetCoordinates_Leading_NoCut->SetTitle("Fastjet-returned leading jet coordinates in Hybrid Events w/o cut");
    h2D_eta_vs_phi_JetCoordinates_Leading_NoCut->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
    h2D_eta_vs_phi_JetCoordinates_Leading_NoCut->GetXaxis()->SetTitle("#eta");
    h2D_eta_vs_phi_JetCoordinates_Leading_NoCut->GetYaxis()->SetTitleSize(DEF_AxisLabelSize);
    h2D_eta_vs_phi_JetCoordinates_Leading_NoCut->GetYaxis()->SetTitle("#phi[rad]");
    h2D_eta_vs_phi_JetCoordinates_Leading_NoCut->GetZaxis()->SetTitleSize(DEF_AxisLabelSize);
    h2D_eta_vs_phi_JetCoordinates_Leading_NoCut->GetZaxis()->SetTitle("#frac{1}{N_{Events}} #frac{d N}{d #phi d #eta}");

    //write global results
    OutputFile->cd();
    Results->cd();
    h2D_pt_vs_eta_LargeGap->Write();
    h2D_pt_vs_eta_SmallGap->Write();
    h2D_eta_vs_dphi_LargeGap->Write();
    h2D_eta_vs_dphi_SmallGap->Write();
    h1D_JetPopulation_Leading->Write();
    h1D_JetPopulation_Subleading->Write();
    h1D_JetPt_Leading->Write();
    h1D_JetPt_Subleading->Write();


    cout << "Write ..." << endl;

    OutputFile->Close();

    return 1;
}

Int_t DijetAna_DiffWakeAna(const TString DataFile, double i_PtRange, double i_LowPtCut) {

    #define DEF_AxisLabelSize 0.05
    #define DEF_HistoTitleSize 0.1
    #define DEF_Rebin 16

    gStyle->SetOptStat(0);
    SetRootGraphicStyle();

    //open the root file
    TFile *file = TFile::Open(DataFile);

    if (!file || file->IsZombie()) {
        std::cout << "Error while opening the file!" << std::endl;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /*
    /
    /                                   ANALYZE THE DIRECT JET RECOIL SITES AND THE CORRESPONDING BACKGROUND 
    /
    /
    /
    /
    /
    *//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //load histograms from root file
    TH2D* h2D_pt_vs_eta_SmallGap = (TH2D*)file->Get("h2D_pt_vs_eta_SmallGap");
    if(!h2D_pt_vs_eta_SmallGap){ cout << "Small Gap histogram not found!" << endl; return 0;}
    TH2D* h2D_pt_vs_eta_LargeGap = (TH2D*)file->Get("h2D_pt_vs_eta_LargeGap");
    if(!h2D_pt_vs_eta_LargeGap){ cout << "Large Gap histogram not found!" << endl; return 0;}
    

    //get X(pt) binning
    double BinsPerMomentum = h2D_pt_vs_eta_SmallGap->GetNbinsX() / h2D_pt_vs_eta_SmallGap->GetXaxis()->GetXmax();

    //configure and draw canvas
    TCanvas* can_ParticleMultiplicities = new TCanvas("ParticleMultiplicities","ParticleMultiplicities",1000,1000);
    can_ParticleMultiplicities->Divide(3,1, 0, 0);

    //SMALL GAP
    can_ParticleMultiplicities->cd(1);

    TH1D* h1D_PartMult_SmallGap = h2D_pt_vs_eta_SmallGap->ProjectionY("h1D_PartMult_SmallGap", i_LowPtCut * BinsPerMomentum, (i_LowPtCut + i_PtRange) * BinsPerMomentum);
    cout << "Lower bin: " << i_LowPtCut*BinsPerMomentum << endl;
    cout << "Upper bin: " << (i_LowPtCut + i_PtRange) * BinsPerMomentum << endl;
    h1D_PartMult_SmallGap->Scale(h2D_pt_vs_eta_SmallGap->GetXaxis()->GetBinWidth(1));
    

    //markings
    h1D_PartMult_SmallGap->SetTitle("Small Gap");
    h1D_PartMult_SmallGap->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
    h1D_PartMult_SmallGap->GetXaxis()->SetTitle("#eta");
    h1D_PartMult_SmallGap->GetYaxis()->SetTitleSize(DEF_AxisLabelSize);
    h1D_PartMult_SmallGap->GetYaxis()->SetTitle("#frac{1}{N_{Small Gap Events}} #frac{d N}{d #eta}");
    h1D_PartMult_SmallGap->Rebin(DEF_Rebin);
    h1D_PartMult_SmallGap->Scale(1./DEF_Rebin);
    h1D_PartMult_SmallGap->DrawCopy();

    can_ParticleMultiplicities->cd(2);

    TH1D* h1D_PartMult_LargeGap = h2D_pt_vs_eta_LargeGap->ProjectionY("h1D_PartMult_LargeGap", i_LowPtCut * BinsPerMomentum, (i_LowPtCut + i_PtRange) * BinsPerMomentum);
    h1D_PartMult_LargeGap->Scale(h2D_pt_vs_eta_LargeGap->GetXaxis()->GetBinWidth(1));
   
    
    h1D_PartMult_LargeGap->SetTitle("Large Gap");
    h1D_PartMult_LargeGap->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
    h1D_PartMult_LargeGap->GetXaxis()->SetTitle("#eta");
    h1D_PartMult_LargeGap->GetYaxis()->SetTitleSize(DEF_AxisLabelSize);
    h1D_PartMult_LargeGap->GetYaxis()->SetTitle("#frac{1}{N_{Large Gap Events}} #frac{d N}{d #eta}");
    h1D_PartMult_LargeGap->Rebin(DEF_Rebin);
    h1D_PartMult_LargeGap->Scale(1./DEF_Rebin);
    h1D_PartMult_LargeGap->DrawCopy();

    //draw difference
    can_ParticleMultiplicities->cd(3);

    h1D_PartMult_LargeGap->Add(h1D_PartMult_SmallGap, -1);
    h1D_PartMult_LargeGap->SetTitle("Difference");
    h1D_PartMult_LargeGap->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
    h1D_PartMult_LargeGap->GetXaxis()->SetTitle("#eta");
    h1D_PartMult_LargeGap->GetYaxis()->SetTitleSize(DEF_AxisLabelSize);
    h1D_PartMult_LargeGap->GetYaxis()->SetTitle("#frac{1}{N_{Large Gap Events}} #frac{d N}{d #eta} - #frac{1}{N_{Small Gap Events}} #frac{d N}{d #eta}");
    h1D_PartMult_LargeGap->DrawCopy();

    

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
    for(int iJet = 0; iJet < SubtractedJets.size(); iJet++)
    {
        if(SubtractedJets[iJet].area() < Jet_area_cut){SubtractedJets.erase(SubtractedJets.begin() + iJet);}
        SubtractedJets[iJet].set_user_info(new MyUserInfo(SubtractedJets[iJet].constituents().size()));

    }

    return sorted_by_pt(SubtractedJets);

}