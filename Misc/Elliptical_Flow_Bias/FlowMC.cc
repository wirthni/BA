#include "./StJetTrackEvent.h"
#include "./StJetTrackEventLinkDef.h"
#include "./functions.h"
#include "TClassRef.h"
#include "Math/Vector4D.h"
#include <algorithm>
#include "./PythiaEvent.h"
#include "./PythiaEvent_LinkDef.h"
//#include "Fastjet/fastjet-3.5.1/include/fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/config.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/Subtractor.hh"
//#include "StPhysicalHelixD.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
using namespace fastjet;
using namespace std;

const char PYTHIA_TREE[]   = "PythiaEvent";
const char PYTHIA_EVENT_BRANCH[] = "Events";

vector<PseudoJet> FindJets(vector<PseudoJet> vec_particles);
double GetPhiByProb(void);
double GetEtaByProb(void);
double GetPtByProb(void);
double GetNByProb(void);

double GetRandomF(double (*i_FuncPtr)(double, double[]), double i_LowerLim, double i_UpperLim, double i_FuncParams[]);
double Function_NormGauss(double x, double params[]);
double Function_PhiByFlow(double x, double params[]);
double Function_FlowByPtAndEta(double x, double params[]);

class MyUserInfo : public fastjet::PseudoJet::UserInfoBase {
public:
    MyUserInfo(int NoOfParticles) 
        : m_NoOfParticles(NoOfParticles){}

    int getNoOfParticles() const { return m_NoOfParticles; }

private:
    int m_NoOfParticles;
};

/*Parameters:
i_NoOfEvents optional to limit the total number of events looped over. If you want to separate the input data, use the last two input parameters, not this one. Leave -1
i_PtLowerLimit is the jet threshold at which it is accepted as a jet
i_AllowV2 if false V2=0 always, else it is generated based on a pion v2 depending on pt
i_AllowV3 if false V3=0 always else ...
i_OutputFileName should be of form "NAME". File type is appended automatically
i_InputDataDivisions number of times the input files should be separated (e.g. 10) ->only 1/10 of the input data is processed
i_InputDataDivisionNumber 1->work off the first 1/10, 2->work off the second 1/10, ... . Goes from 1...i_InputDataDivisions

Additional Info:

Event Plane is always at \phi = 0
*/
Int_t FlowMC(TString i_PathToPYTHIAFiles = "default", Int_t i_NoOfEvents = -1, bool i_GenerateBackground = true, bool i_AllowV2 = false, bool i_AllowV3 = false, TString i_OutputFileName = "Default", bool i_UseHadronInstadOfJet = false, bool i_UsePYTHIAMultipleTimes = false, uint8_t i_InputDataDivisions = 1, uint16_t i_InputDataDivisionNumber = 1)
{

    gStyle->SetOptStat(0);
    SetRootGraphicStyle();

    #define DEF_PYTHIAOversampling 20
    #define DEF_BinningPerUnit 100 //100 is okay
    #define DEF_JetRadius 0.2
    #define DEF_UsePYTHIAJetData true
    #define DEF_OutputEventOverviews false 
    #define DEF_JetLeadingPt 40.0
    #define DEF_JetSubleadingPt 20.0
    #define DEF_HadLeadingPt 25.0
    #define DEF_HadSubleadingPt 15.0
    #define DEF_CorrelationMinPt 1.0
    #define DEF_CorrelationMaxPt 2.0
    #define DEF_MaxPartilclesPerJet 200
    #define DEF_MaxJetPt 150.0

    //other variables
    TH1D* h1D_JetMultiplicityInPsi = new TH1D("h1D_JetMultiplicityInPsi", "h1D_JetMultiplicityInPsi", DEF_BinningPerUnit*2*Pi, -Pi, Pi);

    TH2D* h2D_pt_vs_eta_LargeGap = new TH2D("h2D_pt_vs_eta_LargeGap","h2D_pt_vs_eta_LargeGap",100, 0, 10, DEF_BinningPerUnit * 2 * 0.9,-0.9,0.9);
    TH2D* h2D_pt_vs_eta_SmallGap = new TH2D("h2D_pt_vs_eta_SmallGap","h2D_pt_vs_eta_SmallGap",100, 0, 10, DEF_BinningPerUnit * 2 * 0.9,-0.9,0.9);

    TH2D* h2D_eta_vs_dphi_LargeGap = new TH2D("h2D_eta_vs_dphi_LargeGap","h2D_eta_vs_dphi_LargeGap",DEF_BinningPerUnit * 2 * 0.9, -0.9, 0.9, DEF_BinningPerUnit * 2 * Pi, -Pi/2, 3*Pi/2);
    TH2D* h2D_eta_vs_dphi_SmallGap = new TH2D("h2D_eta_vs_dphi_SmallGap","h2D_eta_vs_dphi_SmallGap",DEF_BinningPerUnit * 2 * 0.9, -0.9, 0.9, DEF_BinningPerUnit * 2 * Pi, -Pi/2, 3*Pi/2);

    //DEBUG REASONS
    TH2D* h2D_deta_vs_dphi_LargeGap = new TH2D("h2D_deta_vs_dphi_LargeGap","h2D_deta_vs_dphi_LargeGap",DEF_BinningPerUnit * 4 * 0.9, -1.8, 1.8, DEF_BinningPerUnit * 2 * Pi, -Pi/2, 3*Pi/2); // contains particle correlations relative to the leading jet
    TH2D* h2D_deta_vs_dphi_SmallGap = new TH2D("h2D_deta_vs_dphi_SmallGap","h2D_deta_vs_dphi_SmallGap",DEF_BinningPerUnit * 4 * 0.9, -1.8, 1.8, DEF_BinningPerUnit * 2 * Pi, -Pi/2, 3*Pi/2); // contains particle correlations relative to the leading jet
    TH1D* h1D_JetPopulation_Leading = new TH1D("h1D_JetPopulation_Leading", "h1D_JetPopulation_Leading", DEF_MaxPartilclesPerJet, 0, DEF_MaxPartilclesPerJet);
    TH1D* h1D_JetPopulation_Subleading = new TH1D("h1D_JetPopulation_Subleading", "h1D_JetPopulation_Subleading", DEF_MaxPartilclesPerJet, 0, DEF_MaxPartilclesPerJet);
    TH1D* h1D_JetPt_Leading = new TH1D("h1D_JetPt_Leading", "h1D_JetPt_Leading", DEF_MaxJetPt, 0, DEF_MaxJetPt);
    TH1D* h1D_JetPt_Subleading = new TH1D("h1D_JetPt_Subleading", "h1D_JetPt_Subleading", DEF_MaxJetPt, 0, DEF_MaxJetPt);  
    TH2D* h2D_eta_vs_phi_JetCoordinates_Leading_NoCut = new TH2D("h2D_eta_vs_phi_JetCoordinates_Leading_NoCut","h2D_eta_vs_phi_JetCoordinates_Leading_NoCut",DEF_BinningPerUnit * 2 * 0.9, -0.9, 0.9, DEF_BinningPerUnit * 2 * Pi, -Pi, Pi); // contains only the jet coordinates before applying the cut
    TH2D* h2D_eta_vs_phi_JetCoordinates_Leading_OnlyEtaCut = new TH2D("h2D_eta_vs_phi_JetCoordinates_Leading_OnlyEtaCut","h2D_eta_vs_phi_JetCoordinates_Leading_OnlyEtaCut",DEF_BinningPerUnit * 2 * 0.9, -0.9, 0.9, DEF_BinningPerUnit * 2 * Pi, -Pi, Pi); // contains only the jet coordinates before applying the cut
    TH2D* h2D_eta_vs_phi_JetCoordinates_Leading_LargeGapCut = new TH2D("h2D_eta_vs_phi_JetCoordinates_Leading_LargeGapCut","h2D_eta_vs_phi_JetCoordinates_Leading_LargeGapCut",DEF_BinningPerUnit * 2 * 0.9, -0.9, 0.9, DEF_BinningPerUnit * 2 * Pi, -Pi, Pi); // contains only the jet coordinates after applying the large gap cut
    TH2D* h2D_eta_vs_phi_JetCoordinates_Leading_SmallGapCut = new TH2D("h2D_eta_vs_phi_JetCoordinates_Leading_SmallGapCut","h2D_eta_vs_phi_JetCoordinates_Leading_SmallGapCut",DEF_BinningPerUnit * 2 * 0.9, -0.9, 0.9, DEF_BinningPerUnit * 2 * Pi, -Pi, Pi); // contains only the jet coordinates after applying the small gap cut
    
    //eta x phi x pt
    TH2D* h2D_GeneratedParticles = new TH2D("h2D_GeneratedParticles","h2D_GeneratedParticles",DEF_BinningPerUnit * 2 * 0.9, -0.9, 0.9, DEF_BinningPerUnit * 2 * Pi, -Pi, Pi); // contains only the jet coordinates before applying the cut

    TRandom TR_Eta;
    TR_Eta.SetSeed(0);
    TRandom TR_EventPlaneCoeffCorrelation;
    TR_EventPlaneCoeffCorrelation.SetSeed(0);
    TRandom TR_EventPlaneAngle;
    TR_EventPlaneAngle.SetSeed(0);
    TChain              *inputPYTHIA;
    PythiaEvent         *PYTHIAEvent;
    PythiaParticle      *PYTHIAParticle;
    Long64_t            file_entries;

    Long64_t            LargeGapEventCounter = 0;
    Long64_t            SmallGapEventCounter = 0;

    //set jet/had pt limit
    double LeadingPtLimit = DEF_JetLeadingPt;
    double SubleadingPtLimit = DEF_JetSubleadingPt;
    if(i_UseHadronInstadOfJet)
    {
        LeadingPtLimit = DEF_HadLeadingPt;
        SubleadingPtLimit = DEF_HadSubleadingPt;
    }

    //open the root file
    TFile *file = TFile::Open("merge_mult_track_pT.root");

    if (!file || file->IsZombie()) {
        std::cout << "Error while opening the file!" << std::endl;
        return 0;
    }

    //get statistics histograms from root file
    TH1D* h1D_NDist = (TH1D*)file->Get("h_mult_particles_used_SE_ME_0");
    if(!h1D_NDist){cout << "Multiplicity histogram not found!" << endl; return 0;}
    TH1D* h1D_PtDist = (TH1D*)file->Get("h_particle_pT");
    if(!h1D_PtDist){cout << "Pt histogram not found!" << endl; return 0;}

    if(DEF_UsePYTHIAJetData)
    {
        /*
        GSI:
        TString pinputdirPYTHIA = "./PYTHIA_Data/";
        //TString SEListPYTHIA = "/filelist_PYTHIA_pp.txt";
        TString SEListPYTHIA = "/filelist_PYTHIA_pp_OnlyHigh.txt";
        TString Pythia_List =  "./PYTHIA_Data";
        */

        //------------------------------------
        // Read input data pythia tree
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
        //------------------------------------
        printf("StJetAnalysis::InitPythiaTree() finished \n");

    }
    TString str_out =  "./" + i_OutputFileName;
    if(i_InputDataDivisions != 1)
    {
        str_out += Form("_%d", i_InputDataDivisionNumber);
    }
    str_out += ".root";
    TFile* OutputFile = new TFile(str_out.Data(),"RECREATE");
    OutputFile->cd();
    cout << "Output file is " << str_out << endl;
    TDirectory *EventHistos = OutputFile->mkdir("Event Histos");
    TDirectory *Results = OutputFile->mkdir("Results");
    
    //determine the events to loop over
    if(i_NoOfEvents != -1) file_entries = i_NoOfEvents;
    Int_t FirstEventNo = 0;
    Int_t LastEventNo = file_entries;
    if(i_InputDataDivisions != 1)
    {
        FirstEventNo = (i_InputDataDivisionNumber - 1) * ceil((float)file_entries / (float)i_InputDataDivisions);
        LastEventNo = i_InputDataDivisionNumber * floor((float)file_entries / (float)i_InputDataDivisions);
    }

    cout << "First event number is " << FirstEventNo << endl;
    cout << "Last event number is " << LastEventNo << endl;

    short Oversampling;
    if(i_UsePYTHIAMultipleTimes) Oversampling = DEF_PYTHIAOversampling;
    else Oversampling = 1;

    cout << "Oversampling in total: " << Oversampling << endl;

    //loop over events
    for(short iRepeat = 1; iRepeat <= Oversampling; iRepeat++)
    {
        cout << "Oversampling turn " << iRepeat << endl;

        for(Int_t iEvent = FirstEventNo; iEvent <= LastEventNo; iEvent++)
        {
            //event output
            if (iEvent != 0  &&  iEvent % 100 == 0)
                cout << "." << flush;
            if (iEvent != 0  &&  iEvent % 1000 == 0)
            {
                if((file_entries-0) > 0)
                {
                    Double_t event_percent = 100.0*((Double_t)(iEvent-FirstEventNo + (iRepeat-1) * (LastEventNo - FirstEventNo)))/((Double_t)((LastEventNo - FirstEventNo)*Oversampling));
                    cout << " " << iEvent << " (" << event_percent << "%) " << "\n" << "==> Processing data " << flush;
                }
            }

            vector<PseudoJet> ParticleVector;//holds all the particles in the event
                    
            Int_t N = h1D_NDist->GetRandom();
            double V2_V3_Phase = TR_EventPlaneCoeffCorrelation.Uniform(-Pi, Pi);
            vector<array<double, 4>> ParticleMemory;//remembers randomly generated particle coordinates (phi, eta, pt) for relevant cuts later in the analysis
            vector<array<double, 4>> HighPtParticles;//save for later in case of hadron analysis. All particles that have pt > DEF_SubleadingPt

            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            /*
            /
            /                                                     
            /                                            GENERATE PBPB BACKGROUND BASED ON PROBABILITY DIST WEIGHTED MC PARTICLES
            /                                                     
            /                                                     
            /                                                     
            /
            /                                                    
            /
            /
            /
            *//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            if(i_GenerateBackground)
            {
                for(Int_t iBGPart = 0; iBGPart < N; iBGPart++)
                {
                    //get particle pt
                    double Pt = h1D_PtDist->GetRandom();
                    //get particle eta
                    double Eta = TR_Eta.Uniform(-0.9, 0.9);
                    //get elliptic flow strength from pt
                    double pseudoparams[2] = {Pt,Eta};
                    double FlowParams = Function_FlowByPtAndEta(Pt, pseudoparams);
                    double FlowParamsArray[3] = {FlowParams * i_AllowV2, FlowParams * i_AllowV3, V2_V3_Phase};

                    double Phi = GetRandomF(Function_PhiByFlow, -Pi, +Pi, FlowParamsArray);
                    if(Phi > Pi) Phi -= 2*Pi;
                    else if(Phi < -Pi) Phi += 2*Pi;
                    //remember particle
                    ParticleMemory.push_back({Phi, Eta, Pt});

                    h2D_GeneratedParticles->Fill(Eta, Phi);

                    //feed jetfinder
                    double p_x = Pt/(sqrt(1+pow(tan(Phi), 2)));//ONLY FOR E CALCULATION - SIGNS NOT NECESSARILY CORRECT
                    double p_y = sqrt(pow(Pt,2) - pow(p_x,2));//ONLY FOR E CALCULATION - SIGNS NOT NECESSARILY CORRECT
                    double p_z = (pow(p_x,2) + pow(p_y,2))/((1/pow(tanh(Eta),2))-1);//ONLY FOR E CALCULATION - SIGNS NOT NECESSARILY CORRECT
                    double E = sqrt(pow(p_x,2) + pow(p_y,2) + pow(p_z,2) + pow(0.1349766,2));
                    ROOT::Math::PtEtaPhiEVector v(Pt, Eta, Phi, E);
                    p_x = v.Px();
                    p_y = v.Py();
                    p_z = v.Pz();

                    if(i_UseHadronInstadOfJet && Pt >= SubleadingPtLimit) HighPtParticles.push_back({Phi, Eta, Pt});
                    else if(!i_UseHadronInstadOfJet)ParticleVector.push_back(PseudoJet(p_x, p_y, p_z, E));

                }
            }

            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            /*
            /
            /                                                     
            /                                            OVERLAY WITH PYTHIA EVENT TO HAVE SIMULATED PBPB EVENT
            /                                                     
            /                                                     
            /                                                     
            /
            /                                                    
            /
            /
            /
            *//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            if(DEF_UsePYTHIAJetData)
            {

                if (!inputPYTHIA->GetEntry( iEvent ))
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
                    //if((fabs(pt_track) < 0.15)) continue;
                    Double_t phi_track = TLV_pythia_particle.Phi();
                    Double_t eta_track = TLV_pythia_particle.Eta();
                    Double_t px_track  = TLV_pythia_particle.Px();
                    Double_t py_track  = TLV_pythia_particle.Py();
                    Double_t pz_track  = TLV_pythia_particle.Pz();
                    Double_t E_track   = TLV_pythia_particle.E();

                    // track cuts 
                    if((fabs(pt_track) < 0.15)) continue;
                    if((fabs(eta_track)) > 0.9) continue;

                    ParticleMemory.push_back({phi_track, eta_track, pt_track});

                    h2D_GeneratedParticles->Fill(eta_track, phi_track);

                    if(i_UseHadronInstadOfJet && pt_track >= SubleadingPtLimit) HighPtParticles.push_back({phi_track, eta_track, pt_track});
                    else if(!i_UseHadronInstadOfJet) ParticleVector.push_back(PseudoJet(px_track, py_track, pz_track, E_track));

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

                if(i_UseHadronInstadOfJet)
                {
                    if(HighPtParticles.size() == 0) continue;

                    do{

                        Int_t HighestPtIndex = 0;

                        //search for highest momentum and choose corresponding particle as a primary vertex
                        for(Int_t i=0; i<HighPtParticles.size(); i++)
                        {
                            if(HighPtParticles[i][2] > HighPtParticles[HighestPtIndex][2])
                            {
                                HighestPtIndex = i;
                            }
                        }

                        //push particle into jet vector
                        ROOT::Math::PtEtaPhiMVector p4(HighPtParticles[HighestPtIndex][2], HighPtParticles[HighestPtIndex][1], HighPtParticles[HighestPtIndex][0], 0.1349766); 
                        JetVector.push_back(PseudoJet(p4.Px(), p4.Py(), p4.Pz(), p4.E()));

                        //pop jet particle
                        HighPtParticles.erase(HighPtParticles.begin() + HighestPtIndex);

                        //look for other high pt particles in a radius of DEF_JetRadius because they belong to the same jet
                        for(Int_t i=0; i<HighPtParticles.size(); i++)
                        {

                            double DeltaEta = abs(JetVector.back().eta() - HighPtParticles[i][1]);
                            double DeltaPhi = abs(JetVector.back().phi_std() - HighPtParticles[i][0]);

                            //since a jet at pi and another one at -pi are equivalent, normalize
                            if(DeltaPhi > Pi) DeltaPhi -= 2* Pi;

                            double DeltaR = TMath::Sqrt(DeltaEta*DeltaEta + DeltaPhi*DeltaPhi);

                            if(DeltaR <= DEF_JetRadius)
                            {
                                //pop candidate
                                HighPtParticles.erase(HighPtParticles.begin()+i);
                                i = -1;//start again for safe
                            }
                        }

                    }while(HighPtParticles.size() > 0);

                    //sort high pt hadrons by size to have same structure as FindJets return type
                    JetVector = sorted_by_pt(JetVector);
                    //make eta cut
                    Selector Fiducial_cut_selector = SelectorAbsEtaMax(0.9 - DEF_JetRadius); // Fiducial cut for jets
                    JetVector = Fiducial_cut_selector(JetVector);

                }
                else
                {
                    JetVector = FindJets(ParticleVector);
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

                if(JetVector.size() >= 1)h2D_eta_vs_phi_JetCoordinates_Leading_NoCut->Fill(JetVector[0].rap(), JetVector[0].phi_std());

                if(JetVector.size() < 2 )continue;

                if(JetVector[0].eta() < 0) continue;

                h2D_eta_vs_phi_JetCoordinates_Leading_OnlyEtaCut->Fill(JetVector[0].eta(), JetVector[0].phi_std());

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

                // Fill leading and subleading jet population
                if(!i_UseHadronInstadOfJet)
                {
                    h1D_JetPopulation_Leading->Fill(JetVector[0].user_info<MyUserInfo>().getNoOfParticles());
                    h1D_JetPopulation_Subleading->Fill(JetVector[1].user_info<MyUserInfo>().getNoOfParticles());
                }

                if(LargeGap)
                {
                    for(const auto& particle : ParticleMemory)
                    {
                        // Fill analysis plot
                        double DeltaPhi = particle[0] - JetVector[0].phi_std();
                        if(DeltaPhi > Pi) DeltaPhi -= 2*Pi;
                        else if(DeltaPhi < -Pi) DeltaPhi += 2*Pi;
                        if(abs(DeltaPhi) <= Pi/2) h2D_pt_vs_eta_LargeGap->Fill(particle[2], particle[1]);

                        // Fill event overview
                        #if DEF_OutputEventOverviews
                            h2D_phi_vs_eta->Fill(particle[0], particle[1], particle[2]);
                        #endif

                        // Fill correlation relative to leading jet in phi
                        if(DeltaPhi < -Pi/2) DeltaPhi += 2*Pi;
                        else if(DeltaPhi > 3*Pi/2) DeltaPhi -= 2*Pi;
                        if(particle[2] >= DEF_CorrelationMinPt && particle[2] <= DEF_CorrelationMaxPt)h2D_eta_vs_dphi_LargeGap->Fill(particle[1], DeltaPhi);

                        // Fill correlation relative to leading jet
                        if(DeltaPhi < -Pi/2) DeltaPhi += 2*Pi;
                        else if(DeltaPhi > 3*Pi/2) DeltaPhi -= 2*Pi;
                        h2D_deta_vs_dphi_LargeGap->Fill(particle[1] - JetVector[0].eta(), DeltaPhi);

                    }

                    // Fill large gap cut jet coordinates
                    h2D_eta_vs_phi_JetCoordinates_Leading_LargeGapCut->Fill(JetVector[0].eta(), JetVector[0].phi_std());

                    LargeGapEventCounter++;
                }
                else
                {
                    for(const auto& particle : ParticleMemory)
                    {
                        // Fill analysis plot
                        double DeltaPhi = particle[0] - JetVector[0].phi_std();
                        if(DeltaPhi > Pi) DeltaPhi -= 2*Pi;
                        else if(DeltaPhi < -Pi) DeltaPhi += 2*Pi;
                        if(abs(DeltaPhi) <= Pi/2) h2D_pt_vs_eta_SmallGap->Fill(particle[2], particle[1]);

                        // Fill event overview
                        #if DEF_OutputEventOverviews
                            h2D_phi_vs_eta->Fill(particle[0], particle[1], particle[2]);
                        #endif

                        // Fill correlation relative to leading jet in phi
                        if(DeltaPhi < -Pi/2) DeltaPhi += 2*Pi;
                        else if(DeltaPhi > 3*Pi/2) DeltaPhi -= 2*Pi;
                        if(particle[2] >= DEF_CorrelationMinPt && particle[2] <= DEF_CorrelationMaxPt)h2D_eta_vs_dphi_SmallGap->Fill(particle[1], DeltaPhi);

                        // Fill correlation relative to leading jet
                        if(DeltaPhi < -Pi/2) DeltaPhi += 2*Pi;
                        else if(DeltaPhi > 3*Pi/2) DeltaPhi -= 2*Pi;
                        h2D_deta_vs_dphi_SmallGap->Fill(particle[1] - JetVector[0].eta(), DeltaPhi);

                    }

                    // Fill small gap cut jet coordinates
                    h2D_eta_vs_phi_JetCoordinates_Leading_SmallGapCut->Fill(JetVector[0].eta(), JetVector[0].phi_std());

                    SmallGapEventCounter++;
                }

                //write event results
                #if DEF_OutputEventOverviews
                    EventHistos->cd();
                    h2D_phi_vs_eta->Write();
                    delete h2D_phi_vs_eta;
                #endif

            }

        }

    }

    cout << "In the end " << ((double)(SmallGapEventCounter + LargeGapEventCounter)/(double)((LastEventNo - FirstEventNo) * Oversampling)) * 100 << " percent of events were taken into account." << endl;

    // Normalize histograms by number of events
    h2D_pt_vs_eta_LargeGap->Scale(1./(double)LargeGapEventCounter);
    h2D_pt_vs_eta_SmallGap->Scale(1./(double)SmallGapEventCounter);

    h2D_eta_vs_dphi_LargeGap->Scale(1./(double)LargeGapEventCounter);
    h2D_eta_vs_dphi_SmallGap->Scale(1./(double)SmallGapEventCounter);

    h1D_JetPopulation_Leading->Scale(1./(double)(LargeGapEventCounter + SmallGapEventCounter));
    h1D_JetPopulation_Subleading->Scale(1./(double)(LargeGapEventCounter + SmallGapEventCounter));

    h1D_JetPt_Leading->Scale(1./(double)(LargeGapEventCounter + SmallGapEventCounter));
    h1D_JetPt_Subleading->Scale(1./(double)(LargeGapEventCounter + SmallGapEventCounter));

    h2D_GeneratedParticles->Scale(1./(double)(LargeGapEventCounter + SmallGapEventCounter));

    // Rebin Histograms if necessary
    h2D_eta_vs_dphi_LargeGap->Rebin2D(8,8);
    h2D_eta_vs_dphi_SmallGap->Rebin2D(8,8);
    h2D_eta_vs_dphi_LargeGap->Scale(1.0/8.0);
    h2D_eta_vs_dphi_SmallGap->Scale(1.0/8.0);

    // Divide by bin widths
    h2D_eta_vs_dphi_LargeGap->Scale(1./(h2D_eta_vs_dphi_LargeGap->GetXaxis()->GetBinWidth(0) * h2D_eta_vs_dphi_LargeGap->GetYaxis()->GetBinWidth(0)));
    h2D_eta_vs_dphi_SmallGap->Scale(1./(h2D_eta_vs_dphi_SmallGap->GetXaxis()->GetBinWidth(0) * h2D_eta_vs_dphi_SmallGap->GetYaxis()->GetBinWidth(0)));

    h2D_pt_vs_eta_LargeGap->Scale(1./(h2D_pt_vs_eta_LargeGap->GetXaxis()->GetBinWidth(0) * h2D_pt_vs_eta_LargeGap->GetYaxis()->GetBinWidth(0)));
    h2D_pt_vs_eta_SmallGap->Scale(1./(h2D_pt_vs_eta_SmallGap->GetXaxis()->GetBinWidth(0) * h2D_pt_vs_eta_SmallGap->GetYaxis()->GetBinWidth(0)));

    h2D_deta_vs_dphi_LargeGap->Scale(1./(h2D_deta_vs_dphi_LargeGap->GetXaxis()->GetBinWidth(0) * h2D_deta_vs_dphi_LargeGap->GetYaxis()->GetBinWidth(0)));
    h2D_deta_vs_dphi_SmallGap->Scale(1./(h2D_deta_vs_dphi_SmallGap->GetXaxis()->GetBinWidth(0) * h2D_deta_vs_dphi_SmallGap->GetYaxis()->GetBinWidth(0)));

    h1D_JetPopulation_Leading->Scale(1./h1D_JetPopulation_Leading->GetXaxis()->GetBinWidth(0));
    h1D_JetPopulation_Subleading->Scale(1./h1D_JetPopulation_Subleading->GetXaxis()->GetBinWidth(0));
    
    h1D_JetPt_Leading->Scale(1./h1D_JetPt_Leading->GetXaxis()->GetBinWidth(0));
    h1D_JetPt_Subleading->Scale(1./h1D_JetPt_Subleading->GetXaxis()->GetBinWidth(0));

    h2D_eta_vs_phi_JetCoordinates_Leading_NoCut->Scale(1./(h2D_eta_vs_phi_JetCoordinates_Leading_NoCut->GetXaxis()->GetBinWidth(0) * h2D_eta_vs_phi_JetCoordinates_Leading_NoCut->GetYaxis()->GetBinWidth(0)));
    h2D_eta_vs_phi_JetCoordinates_Leading_OnlyEtaCut->Scale(1./(h2D_eta_vs_phi_JetCoordinates_Leading_OnlyEtaCut->GetXaxis()->GetBinWidth(0) * h2D_eta_vs_phi_JetCoordinates_Leading_OnlyEtaCut->GetYaxis()->GetBinWidth(0)));
    h2D_eta_vs_phi_JetCoordinates_Leading_LargeGapCut->Scale(1./(h2D_eta_vs_phi_JetCoordinates_Leading_LargeGapCut->GetXaxis()->GetBinWidth(0) * h2D_eta_vs_phi_JetCoordinates_Leading_LargeGapCut->GetYaxis()->GetBinWidth(0)));
    h2D_eta_vs_phi_JetCoordinates_Leading_SmallGapCut->Scale(1./(h2D_eta_vs_phi_JetCoordinates_Leading_SmallGapCut->GetXaxis()->GetBinWidth(0) * h2D_eta_vs_phi_JetCoordinates_Leading_SmallGapCut->GetYaxis()->GetBinWidth(0)));

    h2D_GeneratedParticles->Scale(1./(h2D_GeneratedParticles->GetXaxis()->GetBinWidth(0)));

    //write global results
    Results->cd();
    h2D_pt_vs_eta_LargeGap->Write();
    h2D_pt_vs_eta_SmallGap->Write();
    h2D_eta_vs_dphi_LargeGap->Write();
    h2D_eta_vs_dphi_SmallGap->Write();
    h2D_deta_vs_dphi_LargeGap->Write();
    h2D_deta_vs_dphi_SmallGap->Write();
    h1D_JetPopulation_Leading->Write();
    h1D_JetPopulation_Subleading->Write();
    h1D_JetPt_Leading->Write();
    h1D_JetPt_Subleading->Write();
    h2D_eta_vs_phi_JetCoordinates_Leading_NoCut->Write();
    h2D_eta_vs_phi_JetCoordinates_Leading_OnlyEtaCut->Write();
    h2D_eta_vs_phi_JetCoordinates_Leading_LargeGapCut->Write();
    h2D_eta_vs_phi_JetCoordinates_Leading_SmallGapCut->Write();
    h2D_GeneratedParticles->Write();
    OutputFile ->Close();
    
    return 1;
}

/*PARAMETERS
DataFile is data file that was generated by FlowMC. Form: "NAME.root"
i_PtRange: e.g. particles within (1.0, 2.0)GeV/c should be analyzed -> i_PtRange = 1
i_LowPtCut: e.g. particles within (1.0, 2.0)GeV/c should be analyzed ->i_LowPtCut = 1
*/

Int_t FlowMC_Ana(const TString DataFile_SmallGap, const TString DataFile_LargeGap, double i_PtRange, double i_LowPtCut) {

    #define DEF_AxisLabelSize 0.05
    #define DEF_HistoTitleSize 0.1
    #define DEF_Rebin 16

    gStyle->SetOptStat(0);
    SetRootGraphicStyle();

    //open the root file
    TFile *file_SmallGap = TFile::Open(DataFile_SmallGap);

    if (!file_SmallGap || file_SmallGap->IsZombie()) {
        std::cout << "Error while opening the small gap file!" << std::endl;
    }

    TFile *file_LargeGap = TFile::Open(DataFile_LargeGap);

    if (!file_LargeGap || file_LargeGap->IsZombie()) {
        std::cout << "Error while opening the large gap file!" << std::endl;
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
    TH2D* h2D_pt_vs_eta_SmallGap = (TH2D*)file_SmallGap->Get("h2D_pt_vs_eta_SmallGap");
    if(!h2D_pt_vs_eta_SmallGap){ cout << "Small Gap histogram not found!" << endl; return 0;}
    TH2D* h2D_pt_vs_eta_LargeGap = (TH2D*)file_LargeGap->Get("h2D_pt_vs_eta_LargeGap");
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


    cout << "DONE!" << endl;
    file_SmallGap->Close();
    file_LargeGap->Close();

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
    double Phase = params[2];
    return (1/(2*pi))*(1+2*(v_2*cos(2*x) + v_3*cos(3*(x+Phase))));
}

double Function_FlowByPtAndEta(double x, double params[])
{

    double Pt = params[0];
    double Eta = params[1];
    double FlowParams;
    //an approximation for pions(+), out of paper http://arxiv.org/abs/1301.2348v1
    if(x<3.0) FlowParams =  -0.015*pow((x-3),2) + 0.14;
    else FlowParams = 0.14;

    //now reduce based on eta, just an approximation
    FlowParams *= exp(-0.005*pow(Eta,2));

    return FlowParams;
}