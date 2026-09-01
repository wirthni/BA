#include "./StJetTrackEvent.h"
#include "./StJetTrackEventLinkDef.h"
#include "./functions.h"
#include "TClassRef.h"
#include "Math/Vector4D.h"
#include <algorithm>
#include <cstdlib>
#include "./PythiaEvent.h"
#include "./PythiaEvent_LinkDef.h"
#include "Fastjet/fastjet-3.5.1/include/fastjet/config.h"
#include "Fastjet/fastjet-3.5.1/include/fastjet/ClusterSequenceArea.hh"
#include "Fastjet/fastjet-3.5.1/include/fastjet/PseudoJet.hh"
#include "Fastjet/fastjet-3.5.1/include/fastjet/Selector.hh"
#include "Fastjet/fastjet-3.5.1/include/fastjet/tools/Subtractor.hh"
#include "Fastjet/fastjet-3.5.1/include/fastjet/tools/JetMedianBackgroundEstimator.hh"

using namespace fastjet;
using namespace std;
const char PYTHIA_TREE[]   = "PythiaEvent";
const char PYTHIA_EVENT_BRANCH[] = "Events";

vector<PseudoJet> FindJets(const vector<PseudoJet> vec_particles);
double GetPhiByProb(void);
double GetEtaByProb(void);
double GetPtByProb(void);
double GetNByProb(void);

double GetRandomF(double (*i_FuncPtr)(double, double[]), double i_LowerLim, double i_UpperLim, double i_FuncParams[]);
double Function_NormGauss(double x, double params[]);
double Function_PhiByFlow(double x, double params[]);
double Function_v2v3ByPtAndEta(double x, double params[]);
template <typename T> int sgn(T val);
void PrintProgress(int i_Iteration, int i_Total);
void PrintInfo(TString i_String);

/* Parameters:
i_InputFile: Format xxx.root
*/
Int_t FlowMC(TString i_PathToPYTHIAFiles = "default", const Int_t i_NoOfEvents = 5000,bool i_UsePYTHIAData = true, bool i_GenerateBackground = true, bool i_AllowV2 = false, bool i_AllowV3 = false, TString i_OutputFileName = "Default", bool i_UseHadronInstadOfJet = false, bool i_UsePYTHIAMultipleTimes = false, uint8_t i_InputDataDivisions = 1, uint16_t i_InputDataDivisionNumber = 1)
{
    //CONFIGURE
    bool i_MakeCrosscheckAnalysis = false;

    gStyle->SetOptStat(0);
    SetRootGraphicStyle();

    #define DEF_BinningPerUnit 100
    #define DEF_JetRadius 0.2
    #define DEF_OutputEventOverviews false
    #define DEF_HadLeadingPt 100.0
    #define DEF_BackgroundLimit 4.0 //GeV
    #define DEF_CorrelationMinPt 0.0
    #define DEF_CorrelationMaxPt 2.0
    #define DEF_MaxPartilclesPerJet 200
    #define DEF_MaxJetPt 150.0
    #define DEF_AxisLabelSize 0.05
    #define DEF_AxisNumbersSize 0.04
    #define DEF_HistoTitleSize 0.1
    #define DEF_EventTolerance_Multiplicity 10
    #define DEF_EventTolerance_Psi2 Pi/10
    #define DEF_RecoilCorrelationFrameSize 0.5
    #define DEF_MaxParticlePtInCorrelation 5.0
    #define DEF_BroadCorrelationDeltaPhi Pi/4
    #define DEF_NarrowCorrelationDeltaPhi Pi/8

    int TotalEvents = 0;
    int ProcessedEvents = 0;
    int EventsSkipped   = 0;
    int RunNumber = -1;    

    Long64_t LargeGapEventCounter = 0;
    Long64_t SmallGapEventCounter = 0;

    /*END USER VARIABLES*/

    //APPLICATION VARIABLES
    TH2F* h2F_BroadCorrelationDifference_eta_vs_pT = new TH2F("h2F_BroadCorrelationDifference_eta_vs_pT", "h2F_BroadCorrelationDifference_eta_vs_pT", 50, -0.9, 0.9, 30, 0, 3);
    TH2F* h2F_NarrowCorrelationDifference_eta_vs_pT = new TH2F("h2F_NarrowCorrelationDifference_eta_vs_pT", "h2F_NarrowCorrelationDifference_eta_vs_pT", 50, -0.9, 0.9, 30, 0, 3);
    TH1F* h1F_TrackPtSpectrumRecoilArea = new TH1F("h1F_TrackPtSpectrumRecoilArea", "h1F_TrackPtSpectrumRecoilArea", 2000, 0, 15);
    TH1F* h1F_TrackPtSpectrumReferenceArea = new TH1F("h1F_TrackPtSpectrumReferenceArea", "h1F_TrackPtSpectrumReferenceArea", 2000, 0, 15);
    TH1F* h1F_TrackPtSpectrumDifference = new TH1F("h1F_TrackPtSpectrumDifference", "h1F_TrackPtSpectrumDifference", 2000, 0, 15);
    TH1F* h1F_N_ev_minus_N_ref = new TH1F("h1F_N_ev_minus_N_ref", "h1F_N_ev_minus_N_ref", 100, -50, 50);
    TH1F* h1F_NoOfHighPtParticles = new TH1F("h1F_NoOfHighPtParticles", "h1F_NoOfHighPtParticles", 50, 0, 50);

    TH1F* h1F_LowPtRelToHadron = new TH1F("h1F_LowPtRelToHadron", "h1F_LowPtRelToHadron", 50, -Pi, Pi);
    int AnalyzedEventsCounter = 0;

    // Definitions
    std::unordered_map<int, std::vector<int>> CollisionMap; //holds the Collision ID and underneath the track indices inside of the whole trackChain (0,1,2,...,2644) for example
    std::unordered_set<int> finishedIDs;
    int currentID = -1;
    bool grouped  = true;
    int max_ID    = 0;
    int N_events_skipped = 0;

    //Dijet ana variables
    

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /*
    /
    /                                                     
    /                                                   Prepare everything for event generation 
    /                                                     
    /                                                     
    /                                                     
    /
    /                                                    
    /
    /
    /
    *//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    TRandom TR_Eta;
    TR_Eta.SetSeed(0);
    TRandom TR_V2Angle;
    TR_V2Angle.SetSeed(0);
    TRandom TR_V3Angle;
    TR_V3Angle.SetSeed(0);

    //Get N and p_t statistics
    TFile *file = TFile::Open("merge_mult_track_pT.root");
    if (!file || file->IsZombie()) {
        std::cout << "Error while opening the file!" << std::endl;
        return 0;
    }
    TH1D* h1D_NDist = (TH1D*)file->Get("h_mult_particles_used_SE_ME_0");
    if(!h1D_NDist){cout << "Multiplicity histogram not found!" << endl; return 0;}
    TH1D* h1D_PtDist = (TH1D*)file->Get("h_particle_pT");
    if(!h1D_PtDist){cout << "Pt histogram not found!" << endl; return 0;}

    TChain              *inputPYTHIA;
    PythiaEvent         *PYTHIAEvent;
    PythiaParticle      *PYTHIAParticle;
    Long64_t            file_entries;

    TString pinputdirPYTHIA = i_PathToPYTHIAFiles;
    //TString SEListPYTHIA = "filelist_PYTHIA_pp.txt";
    TString SEListPYTHIA = "filelist_PYTHIA_pp_OnlyHigh.txt";
    TString Pythia_List =  i_PathToPYTHIAFiles;
    Pythia_List+=SEListPYTHIA;
    Int_t file_loop = 0;
    if (!Pythia_List.IsNull())   // if input file is ok
    {
        cout << "Open file list " << Pythia_List << endl;
        ifstream in(Pythia_List);  // input stream
        if(in)
        {
            cout << "file list is ok" << endl;
            inputPYTHIA = new TChain(PYTHIA_TREE , PYTHIA_TREE );
            char str[255];       // char array for each file name
            Long64_t entries_save = 0;
            while(in)
            {
                in.getline(str,255);  // take the lines of the file list
                if(str[0] != 0)
                {
                    TString addfile = str;
                    addfile = pinputdirPYTHIA+addfile;
                    //cout << "addfile: " << addfile.Data() << endl;
                    Long64_t file_entries;
                    //if(file_loop == file_select) inputPYTHIA ->AddFile(addfile.Data(),-1, PYTHIA_TREE );
                    inputPYTHIA->AddFile(addfile.Data(), -1, PYTHIA_TREE);
                    file_entries = inputPYTHIA->GetEntries();
                    cout << "File added to data chain: " << addfile.Data() << " with " << file_entries-entries_save << " entries" << endl;
                    entries_save = file_entries;
                    file_loop++;
                }
            }
            file_entries   = inputPYTHIA->GetEntries();
            cout << "Number of entries in chain: " << file_entries << endl;

            PYTHIAEvent = new PythiaEvent();
            inputPYTHIA  ->SetBranchAddress( PYTHIA_EVENT_BRANCH, &PYTHIAEvent);
        }
        else
        {
            cout << "PYTHIA file has not been opened correclty. Maybe something wrong with the path?" << endl;
            return 0;
        }
    }

    TRandom3 RandomFirstEvent(0);
    Int_t FirstEventNo = RandomFirstEvent.Integer(file_entries - i_NoOfEvents);

    struct struct_EventProperties
    {
        int Multiplicity = 0;
        float Psi2 = -100;
        vector<PseudoJet> HighPtParticles;
        bool IsHighPtJetEvent = false;
        int SuitableReferenceEvent = -1;
        
    };
    struct_EventProperties EventProperties[i_NoOfEvents];
    
    vector<PseudoJet> ParticleVector[i_NoOfEvents];//remembers randomly generated particle coordinates (phi, eta, pt) for relevant cuts later in the analysis


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /*
    /
    /                                                     
    /                                                   Get PYTHIA data  
    /                                                     
    /                                                     
    /                                                     
    /
    /                                                    
    /
    /
    /
    *//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    PrintInfo("Get PYTHIA Data");

    for(int iEvent = 0; iEvent < i_NoOfEvents; iEvent++)
    {
        if (!inputPYTHIA->GetEntry( FirstEventNo + iEvent ))
        {
            cout << "Error" << endl;
            continue;
        }

        //--------------------
        // Event information
        Int_t   N_Particles     = PYTHIAEvent->getNumParticles();
        if(N_Particles <= 0)
        {
            printf("WARNING: iEvent: %d has no entries! \n",iEvent);
        }

        for(Int_t i_Particle = 0; i_Particle < N_Particles; i_Particle++)
        {    
            // Particle information
            // Particle Level
            PYTHIAParticle = PYTHIAEvent  ->getParticle(i_Particle);
            TLorentzVector TLV_pythia_particle = PYTHIAParticle -> get_TLV_part();

            Double_t pt_track  = TLV_pythia_particle.Pt();
            if((fabs(pt_track) < 0.15)) continue;
            Double_t phi_track = TLV_pythia_particle.Phi();
            Double_t eta_track = TLV_pythia_particle.Eta();
            Double_t px_track  = TLV_pythia_particle.Px();
            Double_t py_track  = TLV_pythia_particle.Py();
            Double_t pz_track  = TLV_pythia_particle.Pz();
            Double_t m_track = 0.134977;
            Double_t E_track = sqrt(pow(TLV_pythia_particle.P(),2) + pow(m_track, 2));

            // track cuts 
            if((fabs(pt_track) < 0.15)) continue;
            if((fabs(eta_track)) > 0.9) continue;

            ParticleVector[iEvent].push_back(PseudoJet(px_track, py_track, pz_track, E_track));

            if(pt_track >= DEF_HadLeadingPt)
            {
                EventProperties[iEvent].HighPtParticles.push_back(PseudoJet(px_track, py_track, pz_track, E_track));
            }
            
        }
    
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /*
    /
    /                                                     
    /                              GENERATE PBPB BACKGROUND BASED ON PROBABILITY DIST WEIGHTED MC PARTICLES
    /                                                     
    /                                                     
    /                                                     
    /
    /                                                    
    /
    /
    /
    *//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    PrintInfo("Generate Background");

    for(int iEvent = 0; iEvent < i_NoOfEvents; iEvent++)
    {

        PrintProgress(iEvent, i_NoOfEvents);

        Int_t N = h1D_NDist->GetRandom();
        double V2_Angle = TR_V2Angle.Uniform(-Pi/2, Pi/2);
        double V3_Angle = TR_V3Angle.Uniform(-Pi/2, Pi/2);

        EventProperties[iEvent].Psi2 = V2_Angle;

        for(Int_t iBGPart = 0; iBGPart < N; iBGPart++)
        {
            //get particle pt
            double Pt;
            do{
                Pt = h1D_PtDist->GetRandom();
            }while(fabs(Pt) < 0.15);
            //get particle eta
            double Eta = TR_Eta.Uniform(-0.9, 0.9);
            //get elliptic flow strength from pt
            double pseudoparams[2] = {Pt,Eta};

            //flow parameter v_2 is generated based on RUN 2 statistics
            double FlowParams = Function_v2v3ByPtAndEta(Pt, pseudoparams);

            double FlowParamsArray[4] = {FlowParams * i_AllowV2, FlowParams * i_AllowV3, V2_Angle, V3_Angle};

            double Phi = GetRandomF(Function_PhiByFlow, -Pi, +Pi, FlowParamsArray);
            if(Phi > Pi) Phi -= 2*Pi;
            else if(Phi < -Pi) Phi += 2*Pi;

            //log particles
            ROOT::Math::PtEtaPhiMVector p4(Pt, Eta, Phi, 0.1349766); //0.1349766
            ParticleVector[iEvent].push_back(PseudoJet(p4.Px(), p4.Py(), p4.Pz(), p4.E()));

            if(Pt >= DEF_HadLeadingPt)
            {
                EventProperties[iEvent].HighPtParticles.push_back(PseudoJet(p4.Px(), p4.Py(), p4.Pz(), p4.E()));
            }

        }
    }

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

    int HighPtEventCounter = 0;

    for (int iEvent = 0; iEvent < i_NoOfEvents; iEvent++) {

        bool m_BadEvent = false;

        if(EventProperties[iEvent].HighPtParticles.size() < 1) continue;

        EventProperties[iEvent].HighPtParticles = sorted_by_pt(EventProperties[iEvent].HighPtParticles);

        EventProperties[iEvent].Multiplicity = ParticleVector[iEvent].size();

        if(abs(EventProperties[iEvent].HighPtParticles[0].eta()) > 0.9) continue;

        //check that there is no hard hadron close to the recoil site
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
            float PhiDistance = EventProperties[iEvent].HighPtParticles[iJet].phi_std() - RecoilSitePhi;
            if(PhiDistance > Pi) PhiDistance -= 2*Pi;
            else if(PhiDistance < -Pi) PhiDistance += 2*Pi;
            if(PhiDistance < DEF_BroadCorrelationDeltaPhi)
            {
                m_BadEvent = true;
                break;
            }
        }

        if(m_BadEvent) continue;
        
        EventProperties[iEvent].IsHighPtJetEvent = true;
        HighPtEventCounter++;
        h1F_NoOfHighPtParticles->Fill(EventProperties[iEvent].HighPtParticles.size());

    }

    cout << "High Pt Events/Total Events: " << HighPtEventCounter << "/" << i_NoOfEvents << ". That is " << ((float)HighPtEventCounter/(float)i_NoOfEvents)*100.0 << " percent." << endl;

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
    #if false
    
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

                //check that there is no jet close to the recoil site
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
                    float PhiDistance = EventProperties[iReferenceEvent].HighPtParticles[iJet].phi_std() - RecoilSitePhi;
                    if(PhiDistance > Pi) PhiDistance -= 2*Pi;
                    else if(PhiDistance < -Pi) PhiDistance += 2*Pi;
                    if(abs(PhiDistance) < DEF_BroadCorrelationDeltaPhi)
                    {
                        BadReferenceEvent = true;
                        break;
                    }

                    //check jet site is free
                    PhiDistance = EventProperties[iReferenceEvent].HighPtParticles[iJet].phi_std() - EventProperties[iEvent].HighPtParticles[0].phi_std();
                    if(PhiDistance > Pi) PhiDistance -= 2*Pi;
                    else if(PhiDistance < -Pi) PhiDistance += 2*Pi;
                    if(abs(PhiDistance) < DEF_BroadCorrelationDeltaPhi)
                    {
                        BadReferenceEvent = true;
                        break;
                    }

                }
                if(BadReferenceEvent) continue;

                EventProperties[iEvent].SuitableReferenceEvent = iReferenceEvent;
                break;
            }

            h1F_N_ev_minus_N_ref->Fill(EventProperties[iEvent].Multiplicity - EventProperties[EventProperties[iEvent].SuitableReferenceEvent].Multiplicity);


        }
        
    cout << "Paired Events/High Pt Events: " << PairedEventsCounter << "/" << HighPtEventCounter << ". That is " << ((float)PairedEventsCounter/(float)HighPtEventCounter)*100.0 << " percent." << endl;

    #endif

    
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

    for(int iEvent = 0; iEvent < i_NoOfEvents; iEvent++)
    {

        PrintProgress(iEvent, i_NoOfEvents);

        if(!EventProperties[iEvent].IsHighPtJetEvent) continue;
        //if(EventProperties[iEvent].SuitableReferenceEvent == -1) continue;

        //Fill correlation histogram
        TH2F* h2F_BroadCorrelationRecoil_eta_vs_pT = new TH2F("h2F_BroadCorrelationRecoil_eta_vs_pT", "h2F_BroadCorrelationRecoil_eta_vs_pT", 50, -0.9, 0.9, 30, 0, 3);
        TH2F* h2F_NarrowCorrelationRecoil_eta_vs_pT = new TH2F("h2F_NarrowCorrelationRecoil_eta_vs_pT", "h2F_NarrowCorrelationRecoil_eta_vs_pT", 50, -0.9, 0.9, 30, 0, 3);

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

        int DEBUG_EventCounter = 0;
        int DEBUG_RefCounter = 0;

        for(const auto& particle : ParticleVector[iEvent])
        {

            //float EtaDistance = particle.eta() - RecoilSiteEta;
            float PhiDistance = particle.phi_std() - RecoilSitePhi;
            if(PhiDistance > Pi) PhiDistance -= 2*Pi;
            else if(PhiDistance < -Pi) PhiDistance += 2*Pi;

            if(particle.pt() < 2.0) h1F_LowPtRelToHadron->Fill(PhiDistance);

            if(abs(PhiDistance) <= DEF_BroadCorrelationDeltaPhi)
            {
                h2F_BroadCorrelationRecoil_eta_vs_pT->Fill(particle.eta(), particle.pt());
                if(abs(PhiDistance) <= DEF_NarrowCorrelationDeltaPhi)
                {
                    h2F_NarrowCorrelationRecoil_eta_vs_pT->Fill(particle.eta(), particle.pt());
                    h1F_TrackPtSpectrumRecoilArea->Fill(particle.pt());//, 1./(float) EventProperties[collision.ColID].Multiplicity);
                }
            }
        }

        //Fill reference correlation histogram
        TH2F* h2F_BroadCorrelationReference_eta_vs_pT = new TH2F("h2F_BroadCorrelationReference_eta_vs_pT", "h2F_BroadCorrelationReference_eta_vs_pT", 50, -0.9, 0.9, 30, 0, 3);
        TH2F* h2F_NarrowCorrelationReference_eta_vs_pT = new TH2F("h2F_NarrowCorrelationReference_eta_vs_pT", "h2F_NarrowCorrelationReference_eta_vs_pT", 50, -0.9, 0.9, 30, 0, 3);

        TRandom3 RandomEvent(0);
        TRandom3 RandomParticle(0);

        //mixed event method
        for(int iMEParticle = 0; iMEParticle < ParticleVector[iEvent].size()/8; iMEParticle++)
        {
            int ME_Event_No = RandomEvent.Integer(i_NoOfEvents);
            bool ParticleFound = false;

            do
            {
                int Particle_No = RandomParticle.Integer(ParticleVector[ME_Event_No].size());
                
                PseudoJet track = ParticleVector[ME_Event_No][Particle_No];

                float PhiDistance = track.phi_std() - RecoilSitePhi;
                if(PhiDistance > Pi) PhiDistance -= 2*Pi;
                else if(PhiDistance < -Pi) PhiDistance += 2*Pi;

                if(abs(PhiDistance) < DEF_NarrowCorrelationDeltaPhi)
                {
                    h2F_NarrowCorrelationReference_eta_vs_pT->Fill(track.eta(), track.pt());
                    h1F_TrackPtSpectrumReferenceArea->Fill(track.pt());
                    ParticleFound = true;
                }

            } while (ParticleFound == false);

        }

        /*
        for(const auto& particle : ParticleVector[EventProperties[iEvent].SuitableReferenceEvent])
        {
            //float EtaDistance = particle.eta() - RecoilSiteEta;
            float PhiDistance = particle.phi_std() - RecoilSitePhi;
            if(PhiDistance > Pi) PhiDistance -= 2*Pi;
            else if(PhiDistance < -Pi) PhiDistance += 2*Pi;

            if(particle.pt() < 2.0) DEBUG_RefCounter++;

            if(abs(PhiDistance) <= DEF_BroadCorrelationDeltaPhi)
            {
                h2F_BroadCorrelationReference_eta_vs_pT->Fill(particle.eta(), particle.pt());
                if(abs(PhiDistance) <= DEF_NarrowCorrelationDeltaPhi)
                {
                    h2F_NarrowCorrelationReference_eta_vs_pT->Fill(particle.eta(), particle.pt());
                    h1F_TrackPtSpectrumReferenceArea->Fill(particle.pt());//, 1./(float) EventProperties[EventProperties[iEvent].SuitableReferenceEvent].Multiplicity);
                }
            }
        }
        */

        h2F_BroadCorrelationDifference_eta_vs_pT->Add(h2F_BroadCorrelationRecoil_eta_vs_pT, +1);
        h2F_BroadCorrelationDifference_eta_vs_pT->Add(h2F_BroadCorrelationReference_eta_vs_pT, -1);

        h2F_NarrowCorrelationDifference_eta_vs_pT->Add(h2F_NarrowCorrelationRecoil_eta_vs_pT, +1);
        h2F_NarrowCorrelationDifference_eta_vs_pT->Add(h2F_NarrowCorrelationReference_eta_vs_pT, -1);

        delete(h2F_BroadCorrelationRecoil_eta_vs_pT);
        delete(h2F_BroadCorrelationReference_eta_vs_pT);
        delete(h2F_NarrowCorrelationRecoil_eta_vs_pT);
        delete(h2F_NarrowCorrelationReference_eta_vs_pT);

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


    h2F_BroadCorrelationDifference_eta_vs_pT->Scale(1./(float)AnalyzedEventsCounter);
    h2F_BroadCorrelationDifference_eta_vs_pT->Scale(1./(h2F_BroadCorrelationDifference_eta_vs_pT->GetXaxis()->GetBinWidth(1) * h2F_BroadCorrelationDifference_eta_vs_pT->GetYaxis()->GetBinWidth(1)));
    h2F_NarrowCorrelationDifference_eta_vs_pT->Scale(1./(float)AnalyzedEventsCounter);
    h2F_NarrowCorrelationDifference_eta_vs_pT->Scale(1./(h2F_NarrowCorrelationDifference_eta_vs_pT->GetXaxis()->GetBinWidth(1) * h2F_NarrowCorrelationDifference_eta_vs_pT->GetYaxis()->GetBinWidth(1)));

    h1F_TrackPtSpectrumDifference->Add(h1F_TrackPtSpectrumRecoilArea, 1);
    h1F_TrackPtSpectrumDifference->Add(h1F_TrackPtSpectrumReferenceArea, -1);

    //Generate output file
    TString InputNumberAsString;
    InputNumberAsString.Form("%d", i_InputDataDivisionNumber);
    TString str_out =  "./Flow_SingleJetAna_" + InputNumberAsString +".root";

    TFile* OutputFile;
    if(gSystem->AccessPathName(str_out))
    {
        OutputFile = new TFile(str_out,"RECREATE");
        OutputFile->cd();
        cout << "Output file is " << str_out << endl;
        h1F_TrackPtSpectrumRecoilArea->Write();
        h1F_TrackPtSpectrumReferenceArea->Write();
        h1F_N_ev_minus_N_ref->Write();
        h1F_NoOfHighPtParticles->Write();
        h2F_BroadCorrelationDifference_eta_vs_pT->Write();
        h2F_NarrowCorrelationDifference_eta_vs_pT->Write();
        h1F_TrackPtSpectrumDifference->Write();

        h1F_LowPtRelToHadron->Write();

        OutputFile->Close();
    }
    else
    {
        OutputFile = TFile::Open(str_out.Data());
        cout << "ADD TO EXISTING OUTPUT FILE" << endl;
        OutputFile->cd();

        TH1F* Summed_h1F_LowPtRelToHadron = (TH1F*) OutputFile->Get("h1F_LowPtRelToHadron");
        Summed_h1F_LowPtRelToHadron->Add(h1F_LowPtRelToHadron);
        Summed_h1F_LowPtRelToHadron->Scale(0.5);

        TH1F* Summed_h1F_TrackPtSpectrumRecoilArea = (TH1F*) OutputFile->Get("h1F_TrackPtSpectrumRecoilArea");
        Summed_h1F_TrackPtSpectrumRecoilArea->Add(h1F_TrackPtSpectrumRecoilArea);
        Summed_h1F_TrackPtSpectrumRecoilArea->Scale(0.5);

        TH1F* Summed_h1F_TrackPtSpectrumReferenceArea = (TH1F*) OutputFile->Get("h1F_TrackPtSpectrumReferenceArea");
        Summed_h1F_TrackPtSpectrumReferenceArea->Add(h1F_TrackPtSpectrumReferenceArea);
        Summed_h1F_TrackPtSpectrumReferenceArea->Scale(0.5);

        TH1F* Summed_h1F_TrackPtSpectrumDifference = (TH1F*) OutputFile->Get("h1F_TrackPtSpectrumDifference");
        Summed_h1F_TrackPtSpectrumDifference->Add(h1F_TrackPtSpectrumDifference);
        Summed_h1F_TrackPtSpectrumDifference->Scale(0.5);

        TH1F* Summed_h1F_N_ev_minus_N_ref = (TH1F*) OutputFile->Get("h1F_N_ev_minus_N_ref");
        Summed_h1F_N_ev_minus_N_ref->Add(h1F_N_ev_minus_N_ref);
        Summed_h1F_N_ev_minus_N_ref->Scale(0.5);

        TH1F* Summed_h1F_NoOfHighPtParticles = (TH1F*) OutputFile->Get("h1F_NoOfHighPtParticles");
        Summed_h1F_NoOfHighPtParticles->Add(h1F_NoOfHighPtParticles);
        Summed_h1F_NoOfHighPtParticles->Scale(0.5);

        TH1F* Summed_h2F_BroadCorrelationDifference_eta_vs_pT = (TH1F*) OutputFile->Get("h2F_BroadCorrelationDifference_eta_vs_pT");
        Summed_h2F_BroadCorrelationDifference_eta_vs_pT->Add(h2F_BroadCorrelationDifference_eta_vs_pT);
        Summed_h2F_BroadCorrelationDifference_eta_vs_pT->Scale(0.5);

        TH1F* Summed_h2F_NarrowCorrelationDifference_eta_vs_pT = (TH1F*) OutputFile->Get("h2F_NarrowCorrelationDifference_eta_vs_pT");
        Summed_h2F_NarrowCorrelationDifference_eta_vs_pT->Add(h2F_NarrowCorrelationDifference_eta_vs_pT);
        Summed_h2F_NarrowCorrelationDifference_eta_vs_pT->Scale(0.5);

        OutputFile = new TFile(str_out.Data(),"RECREATE");
        Summed_h1F_TrackPtSpectrumRecoilArea->Write();
        Summed_h1F_TrackPtSpectrumReferenceArea->Write();
        Summed_h1F_TrackPtSpectrumDifference->Write();
        Summed_h1F_N_ev_minus_N_ref->Write();
        Summed_h1F_NoOfHighPtParticles->Write();
        Summed_h2F_BroadCorrelationDifference_eta_vs_pT->Write();
        Summed_h2F_NarrowCorrelationDifference_eta_vs_pT->Write();
        Summed_h1F_LowPtRelToHadron->Write();

        OutputFile->Close();

    }

    return 1;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
/
/                                                   FUNCTIONS FOR BACKGROUND GENERATION
/
/
/
/
*//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


float RandomFloat(float a, float b) {
    float random = ((float) rand()) / (float) RAND_MAX;
    float diff = b - a;
    float r = random * diff;
    return a + r;
}

#include <random>
double GetRandomF(double (*i_FuncPtr)(double, double[]), double i_LowerLim, double i_UpperLim, double i_FuncParams[])
{
    /*
    ------------------------------------------------------------------------------------------------------------------------------------------------------------------
    RETURNS RANDOM NUMBERS BASED ON A WEIGHT FUNCTION WITH AN ACCURACY OF #define ACCURACY
    ------------------------------------------------------------------------------------------------------------------------------------------------------------------
    INPUT PARAMETERS:
    i_FuncPtr is a function pointer to a (non)normalized function that is the probability density
    i_LowerLim is the lowest number possible that can be returned
    i_UpperLim is the highest number possible that can be returned
    ACCURACY is the number of times the function is sampled within (i_LowerLim, i_UpperLim). Higher values give more accurate results, but are computationally costly
    ------------------------------------------------------------------------------------------------------------------------------------------------------------------
    EXAMPLE CODE:
    int main()
    {
        double FlowParamsArray[2] = {0.25, 0};
        double Val = GetRandomF(Function_NormGauss, -Pi, +Pi, FlowParamsArray); 
        cout << Val << endl;
        return 1;
    }

    double Function_NormGauss(double x, double params[])
    {    
        double Sigma = params[0];
        double Mean = params[1];
        return exp(-pow(x-Mean,2)/(2*pow(Sigma,2)))* 1/(sqrt(2*Pi*pow(Sigma, 2)));
    }

    double GetRandomF(...)
    {
        (...)
    }
    -------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Author: Nicolas Wirth, Nicolas.Wirth@cern.ch
    Utilize and modify freely
    -------------------------------------------------------------------------------------------------------------------------------------------------------------------
    */

    #define ACCURACY 100

    double TotalProb = 0;
    double IntervalLength = (i_UpperLim - i_LowerLim)/ACCURACY;
    for(int8_t i=0; i<=ACCURACY-1; i++)
    {
        TotalProb += i_FuncPtr(i_LowerLim + (0.5+i)*IntervalLength, i_FuncParams) * (i_UpperLim - i_LowerLim)/ACCURACY;
    }

    double Lim = TotalProb * (double)rand() / RAND_MAX;
    for(int8_t i=0; i<=ACCURACY+10; i++)
    {
        if(Lim <=0) return i_LowerLim + (0.5+i)*IntervalLength + RandomFloat(-IntervalLength/2, +IntervalLength/2);
        Lim -= i_FuncPtr(i_LowerLim + (0.5+i)*IntervalLength, i_FuncParams) * (i_UpperLim - i_LowerLim)/ACCURACY;
    }

    return i_LowerLim + (0.5+(ACCURACY-1))*IntervalLength + RandomFloat(-IntervalLength/2, +IntervalLength/2);
}

double Function_NormGauss(double x, double params[])
{
    double Sigma = params[0];
    double Mean = params[1];
    return exp(-pow(x-Mean,2)/(2*pow(Sigma,2)))* 1/(sqrt(2*Pi*pow(Sigma, 2)));
}

double Function_PhiByFlow(double x, double params[])
{
    double v_2 = params[0];
    double v_3 = params[1];
    double v2_Angle = params[2];
    double v3_Angle = params[3];
    return (1/(2*pi))*(1+2*(v_2*cos(2*(x-v2_Angle)) + v_3*cos(3*(x-v3_Angle))));
}

double Function_v2v3ByPtAndEta(double x, double params[])
{

    double Pt = params[0];
    double Eta = params[1];
    double FlowParams;
    //an approximation for pions(+) @Nadine Gruenwald
    FlowParams = 0.063 * pow(x, 1.6) * exp(-0.47*x);

    //now reduce based on eta, a fit to https://www.hepdata.net/record/ins1456145 from -1.25 <= x <= 1.25
    FlowParams *= (1.0 - 0.105*pow(Eta, 2));

    return FlowParams;
}

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

void PrintProgress(int i_Iteration, int i_Total)
{

    int barWidth = 10;
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