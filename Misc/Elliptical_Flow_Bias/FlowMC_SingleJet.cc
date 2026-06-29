#include "./StJetTrackEvent.h"
#include "./StJetTrackEventLinkDef.h"
#include "./functions.h"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/config.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/Subtractor.hh"
//#include "StPhysicalHelixD.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include <algorithm>
#include "TClassRef.h"
#include "Math/Vector4D.h"
#include "./PythiaEvent.h"
#include "./PythiaEvent_LinkDef.h"
//#include "Fastjet/fastjet-3.5.1/include/fastjet/ClusterSequence.hh"
using namespace std;
using namespace fastjet;
#if false
#include "fastjet/config.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/Subtractor.hh"
#include "StPhysicalHelixD.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#endif

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


//JetEventNo gives back the ith event in which a jet with p_t > DEF_HighPtLimit was found (0 is the first one)

int8_t ReadTrees(Int_t i_BinningPerUnit, TString i_PathToPYTHIAFiles = "default", Int_t i_NoOfEvents = -1, bool i_GenerateBackground = true, bool i_AllowV2 = false, bool i_AllowV3 = false, TString i_OutputFileName = "Default", bool i_UseHadronInstadOfJet = false, bool i_UsePYTHIAMultipleTimes = false, uint8_t i_InputDataDivisions = 1, uint16_t i_InputDataDivisionNumber = 1)
{
    //definitions to control data output. Uncomment if necessary

    #define DEF_Output_h2D_JetsDijets
    #define DEF_Output_h2D_JetsAvgSpatialDistribution
    //#define DEF_Output_h2D_CorrelatedDijets
    //#define DEF_Output_ShowJetFinderResults


    #define DEF_HighMomentumJetLimit 80                     //to classify jets based on their primarys pt momentum (GeV)
    #define DEF_MiddleMomentumJetLimit 50
    #define DEF_LowMomentumJetLimit 30
    #define DEF_SecondaryJetLimit 7
    #define DEF_HighestExpectableJetMomentum 230
    #define DEF_HighestExpectableRecoilMomentum 20
    #define DEF_MomentumBinning 1000
    #define DEF_JetRadius 0.4
    #define DEF_DijetCorrelationPhiTolerance 0.1            //allow this tolerance in phi to tell if two jets are correlated
    #define DEF_JetProfileLowPtCut 0                        //only particles above this limit are taken into account while analyzing the jet profile/shape
    //#define DEF_GetCrosscheckInsteadOfRecoil                //if this is active, not the recoil site but a site perpendicular to the jet axis is investigated(eta is not inverted)
    #define DEF_PYTHIAOversampling 5
    #define DEF_BinningPerUnit 100
    #define DEF_JetRadius 0.2
    #define DEF_LeadingPt 40
    #define DEF_SubleadingPt 30
    #define DEF_UsePYTHIAJetData true
    #define DEF_OutputEventOverviews true  

    //OTHER USEFUL VARIABLES
    vector<PseudoJet> v_JetWaitingForBackgroundCapture;    //holds (phi, eta, pt) of the last direct jet events to remember where to take the background in the next suitable event
    int32_t DEBUG_NoOfLoggedJets = 0;
    int32_t DEBUG_NoOfLoggedBg = 0;
    int32_t DEBUG_NoOfJetEvents = 0;

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

    printf("StJetAnalysis::ReadTrees() started \n");

    SetRootGraphicStyle();

    //open the background statistics histograms
    TFile *file = TFile::Open("merge_mult_track_pT.root");

    if (!file || file->IsZombie()) {
        std::cout << "Error while opening the background histograms file!" << std::endl;
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

    //create jet directories
    TDirectory *JetsDir = OutputFile->mkdir("Jets");
    TDirectory *DijetsDir = OutputFile->mkdir("Dijets");
    TDirectory *SuitableDir = OutputFile->mkdir("Suitable Events");
    TDirectory *PtVsEtaDir = OutputFile->mkdir("Pt_vs_Eta");
    TDirectory *PtVsPhiDir = OutputFile->mkdir("Pt_vs_Phi");

    // define histograms etc.
    vector<TH2D*> h2D_pt_vs_eta_JET;
    vector<TH2D*> h2D_pt_vs_phi_JET;
    //jet cross section plots
    TH2D* h2D_pt_vs_eta_Summed_LowPt = new TH2D("h2D_pt_vs_eta_Summed_Lowpt","h2D_pt_vs_eta_Summed_Lowpt",i_BinningPerUnit * 2,-1,1,DEF_MomentumBinning,0,DEF_HighestExpectableJetMomentum);//[0]...[2] contains low, middle, high energetic jets
    TH2D* h2D_pt_vs_phi_Summed_LowPt = new TH2D("h2D_pt_vs_phi_Summed_Lowpt","h2D_pt_vs_phi_Summed_Lowpt",i_BinningPerUnit * 2 * Pi,-Pi,Pi,DEF_MomentumBinning,0,DEF_HighestExpectableJetMomentum);//[0]...[2] contains low, middle, high energetic jets
    TH2D* h2D_pt_vs_eta_Summed_MiddlePt = new TH2D("h2D_pt_vs_eta_Summed_Middlept","h2D_pt_vs_eta_Summed_Middlept",i_BinningPerUnit * 2,-1,1,DEF_MomentumBinning,0,DEF_HighestExpectableJetMomentum);//[0]...[2] contains low, middle, high energetic jets
    TH2D* h2D_pt_vs_phi_Summed_MiddlePt = new TH2D("h2D_pt_vs_phi_Summed_Middlept","h2D_pt_vs_phi_Summed_Middlept",i_BinningPerUnit * 2 * Pi,-Pi,Pi,DEF_MomentumBinning,0,DEF_HighestExpectableJetMomentum);//[0]...[2] contains low, middle, high energetic jets
    TH2D* h2D_pt_vs_eta_Summed_HighPt = new TH2D("h2D_pt_vs_eta_Summed_Highpt","h2D_pt_vs_eta_Summed_Highpt",i_BinningPerUnit * 2,-1,1,DEF_MomentumBinning,0,DEF_HighestExpectableJetMomentum);//[0]...[2] contains low, middle, high energetic jets
    TH2D* h2D_pt_vs_phi_Summed_HighPt = new TH2D("h2D_pt_vs_phi_Summed_Highpt","h2D_pt_vs_phi_Summed_Highpt",i_BinningPerUnit * 2 * Pi,-Pi,Pi,DEF_MomentumBinning,0,DEF_HighestExpectableJetMomentum);//[0]...[2] contains low, middle, high energetic jets
    //jet recoil section plots
    TH2D* h2D_pt_vs_eta_recoil_Summed_LowPt = new TH2D("h2D_pt_vs_eta_recoil_Summed_LowPt","h2D_pt_vs_eta_recoil_Summed_LowPt",i_BinningPerUnit * 2,-1,1,DEF_MomentumBinning,0,DEF_HighestExpectableRecoilMomentum);//[0]...[2] contains low, middle, high energetic jets
    TH2D* h2D_pt_vs_phi_recoil_Summed_LowPt = new TH2D("h2D_pt_vs_phi_recoil_Summed_LowPt","h2D_pt_vs_phi_recoil_Summed_LowPt",i_BinningPerUnit * 2 * Pi,-Pi,Pi,DEF_MomentumBinning,0,DEF_HighestExpectableRecoilMomentum);//[0]...[2] contains low, middle, high energetic jets
    TH2D* h2D_pt_vs_eta_recoil_Summed_MiddlePt = new TH2D("h2D_pt_vs_eta_recoil_Summed_MiddlePt","h2D_pt_vs_eta_recoil_Summed_MiddlePt",i_BinningPerUnit * 2,-1,1,DEF_MomentumBinning,0,DEF_HighestExpectableRecoilMomentum);//[0]...[2] contains low, middle, high energetic jets
    TH2D* h2D_pt_vs_phi_recoil_Summed_MiddlePt = new TH2D("h2D_pt_vs_phi_recoil_Summed_MiddlePt","h2D_pt_vs_phi_recoil_Summed_MiddlePt",i_BinningPerUnit * 2 * Pi,-Pi,Pi,DEF_MomentumBinning,0,DEF_HighestExpectableRecoilMomentum);//[0]...[2] contains low, middle, high energetic jets
    TH2D* h2D_pt_vs_eta_recoil_Summed_HighPt = new TH2D("h2D_pt_vs_eta_recoil_Summed_HighPt","h2D_pt_vs_eta_recoil_Summed_HighPt",i_BinningPerUnit * 2,-1,1,DEF_MomentumBinning,0,DEF_HighestExpectableRecoilMomentum);//[0]...[2] contains low, middle, high energetic jets
    TH2D* h2D_pt_vs_phi_recoil_Summed_HighPt = new TH2D("h2D_pt_vs_phi_recoil_Summed_HighPt","h2D_pt_vs_phi_recoil_Summed_HighPt",i_BinningPerUnit * 2 * Pi,-Pi,Pi,DEF_MomentumBinning,0,DEF_HighestExpectableRecoilMomentum);//[0]...[2] contains low, middle, high energetic jets
    //reference background plots
    TH2D* h2D_pt_vs_eta_Background_Summed = new TH2D("h2D_pt_vs_eta_Background_Summed","h2D_pt_vs_eta_Background_Summed",i_BinningPerUnit * 2,-1,1,DEF_MomentumBinning,0,DEF_HighestExpectableRecoilMomentum);//contains background of non-jet events
    TH2D* h2D_pt_vs_phi_Background_Summed = new TH2D("h2D_pt_vs_phi_Background_Summed","h2D_pt_vs_phi_Background_Summed",i_BinningPerUnit * 2 * Pi,-Pi,Pi,DEF_MomentumBinning,0,DEF_HighestExpectableRecoilMomentum);//contains background of non-jet events
    


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

    for(short iRepeat = 1; iRepeat <= Oversampling; iRepeat++)
    {
        cout << "Oversampling turn " << iRepeat << endl;

        // loop over events 
        for(Int_t iEvent = FirstEventNo; iEvent <= LastEventNo; iEvent++)
        {
            //event output
            if (iEvent != 0  &&  iEvent % 100 == 0)
                cout << "." << flush;
            if (iEvent != 0  &&  iEvent % 1000 == 0)
            {
                if((file_entries-0) > 0)
                {
                    Double_t event_percent = 100.0*((Double_t)(iEvent-0))/((Double_t)(file_entries-0));
                    cout << " " << iEvent << " (" << event_percent << "%) " << "\n" << "==> Processing data " << flush;
                }
            }

            //Event variables
            vector<PseudoJet> ParticleMemory;//remembers randomly generated particle coordinates (phi, eta, pt) for relevant cuts later in the analysis
            vector<PseudoJet> HighPtParticles;//save for later in case of hadron analysis. All particles that have pt > DEF_SubleadingPt

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
                Int_t N = h1D_NDist->GetRandom();
                double V2_V3_Phase = TR_EventPlaneCoeffCorrelation.Uniform(-Pi, Pi);

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

                    //feed jetfinder
                    double p_x = Pt/(sqrt(1+pow(tan(Phi), 2)));//ONLY FOR E CALCULATION - SIGNS NOT NECESSARILY CORRECT
                    double p_y = sqrt(pow(Pt,2) - pow(p_x,2));//ONLY FOR E CALCULATION - SIGNS NOT NECESSARILY CORRECT
                    double p_z = (pow(p_x,2) + pow(p_y,2))/((1/pow(tanh(Eta),2))-1);//ONLY FOR E CALCULATION - SIGNS NOT NECESSARILY CORRECT
                    double E = sqrt(pow(p_x,2) + pow(p_y,2) + pow(p_z,2) + pow(0.1349766,2));
                    ROOT::Math::PtEtaPhiEVector v(Pt, Eta, Phi, E);
                    p_x = v.Px();
                    p_y = v.Py();
                    p_z = v.Pz();

                    if(i_UseHadronInstadOfJet && Pt >= DEF_SubleadingPt) HighPtParticles.push_back({Phi, Eta, Pt});
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

            if (!inputPYTHIA->GetEntry( iEvent )) // take the event -> information is stored in event
                cout << "Error" << endl;

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

                if(i_UseHadronInstadOfJet && pt_track >= DEF_SubleadingPt) HighPtParticles.push_back({phi_track, eta_track, pt_track});
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
            //ParticleMemory.push_back({phi_track, eta_track, pt_track});

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


            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            bool m_SuitableEvent = true;

            //if there are more than two jets, the event is not suitable for recoil ana
            if(JetVector.size() > 1) m_SuitableEvent = false;
            //check if inside of acceptance, else not suitable
            if(abs(JetVector[0].eta()) > 0.9 - DEF_JetRadius) m_SuitableEvent = false;


            #ifdef DEF_Output_h2D_JetsAvgSpatialDistribution

                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                /*
                /
                /               IF A PREVIOUS JET EVENTS WAITS FOR BACKGROUND, CHECK IF THE CURRENT EVENT IS SUITABLE FOR THAT
                /               AND IF YES, TAKE BACKGROUND DATA. SUITABLE MEANS THAT IN THE ...
                /
                /
                /
                /
                *//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                if(v_JetWaitingForBackgroundCapture.size() > 0)
                {
                    double DeltaEta, DeltaPhi, DeltaR;
                    double EtaRecoilSite, PhiRecoilSite;
                    bool m_NotSuitable;

                    for(Int_t i=0; i<v_JetWaitingForBackgroundCapture.size(); i++)
                    {
                        m_NotSuitable = false;

                        //get coordinates of the direct recoil site
                        #ifdef DEF_GetCrosscheckInsteadOfRecoil
                            EtaRecoilSite = (v_JetWaitingForBackgroundCapture[i].eta() + (0.9 - DEF_JetRadius));
                            if(EtaRecoilSite > (0.9 - DEF_JetRadius)) EtaRecoilSite -= 2*(0.9 - DEF_JetRadius);//normalize
                            (v_JetWaitingForBackgroundCapture[i].phi_std() > 0) ? (PhiRecoilSite = v_JetWaitingForBackgroundCapture[i].phi_std() - Pi/2) : (PhiRecoilSite = v_JetWaitingForBackgroundCapture[i].phi_std() + Pi/2);
                        #else    
                        EtaRecoilSite = v_JetWaitingForBackgroundCapture[i].eta();
                        (v_JetWaitingForBackgroundCapture[i].phi_std() > 0) ? (PhiRecoilSite = v_JetWaitingForBackgroundCapture[i].phi_std() - Pi) : (PhiRecoilSite = v_JetWaitingForBackgroundCapture[i].phi_std() + Pi);
                        #endif
                        


                        //loop all the jets of the current event
                        for(Int_t j=0; j<JetVector.size(); j++)
                        {
                            //check that there is no jet in the current event close to the jet waiting for background
                            DeltaEta = abs(JetVector[j].eta() - v_JetWaitingForBackgroundCapture[i].eta());
                            DeltaPhi = abs(JetVector[j].phi_std() - v_JetWaitingForBackgroundCapture[i].phi_std());

                            if(DeltaPhi > Pi) DeltaPhi -= 2* Pi;//normalize

                            DeltaR = TMath::Sqrt(DeltaEta*DeltaEta + DeltaPhi*DeltaPhi);

                            if(DeltaR <= 2*DEF_JetRadius)
                            {
                                m_NotSuitable = true;
                                break;
                            }

                            //check that there is no jet in the current event close to the direct recoil site of the jet waiting for background
                            DeltaEta = abs(JetVector[j].eta() - EtaRecoilSite);
                            DeltaPhi = abs(JetVector[j].phi_std() - PhiRecoilSite);

                            if(DeltaPhi > Pi) DeltaPhi -= 2* Pi;//normalize

                            DeltaR = TMath::Sqrt(DeltaEta*DeltaEta + DeltaPhi*DeltaPhi);

                            if(DeltaR <= 2*DEF_JetRadius)
                            {
                                m_NotSuitable = true;
                                break;
                            }
                        }
                        
                        //dataset is suitable for jet i waiting for background capture
                        if(!m_NotSuitable)
                        {
                            TH2D* h2D_pt_vs_eta_Background = new TH2D("h2D_pt_vs_eta_Background","h2D_pt_vs_eta_Background",i_BinningPerUnit * 2,-1,1,DEF_MomentumBinning,0,DEF_HighestExpectableRecoilMomentum);//contains background of non-jet events
                            TH2D* h2D_pt_vs_phi_Background = new TH2D("h2D_pt_vs_phi_Background","h2D_pt_vs_phi_Background",i_BinningPerUnit * 2 * Pi,-Pi,Pi,DEF_MomentumBinning,0,DEF_HighestExpectableRecoilMomentum);//contains background of non-jet events
        

                            //take background data at the direct recoil site
                            for(Int_t i_Particle = 0; i_Particle < ParticleMemory.size(); i_Particle++)
                            {
                                
                                //only take particles into account that are in DEF_JetRadius to the direct recoil site
                                DeltaEta = abs(ParticleMemory[i_Particle][1] - EtaRecoilSite);
                                DeltaPhi = abs(ParticleMemory[i_Particle][0] - PhiRecoilSite);

                                if(DeltaPhi > Pi) DeltaPhi -= 2* Pi;//normalize

                                DeltaR = TMath::Sqrt(DeltaEta*DeltaEta + DeltaPhi*DeltaPhi);

                                if(DeltaR > DEF_JetRadius)continue;

                                h2D_pt_vs_eta_Background->Fill(ParticleMemory[i_Particle][1] - EtaRecoilSite, ParticleMemory[i_Particle][2]);
                                h2D_pt_vs_phi_Background->Fill(ParticleMemory[i_Particle][0] - PhiRecoilSite, ParticleMemory[i_Particle][2]);
                                
                            }       
                            //delete event from waiting list
                            v_JetWaitingForBackgroundCapture.erase(v_JetWaitingForBackgroundCapture.begin() + i);
                            v_JetWaitingForBackgroundCapture.shrink_to_fit();
                            

                            //append histogram to jet histogram vector
                            h2D_pt_vs_eta_Background_Summed->Add(h2D_pt_vs_eta_Background);
                            h2D_pt_vs_phi_Background_Summed->Add(h2D_pt_vs_phi_Background);

                            DEBUG_NoOfLoggedBg++;

                            delete h2D_pt_vs_eta_Background;
                            delete h2D_pt_vs_phi_Background;

                            break;

                        }
                        
                    }

                }

                if(m_SuitableEvent)
                {

                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    /*
                    /
                    /                               IF THE CURRENT EVENT IS A SUITABLE JET EVENT, GET THE SPATIAL JET
                    /                                           DISTRIBUTION IN ETA&PHI DIRECTIONS 
                    /
                    /
                    /
                    /
                    *//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    double EtaDistance;
                    double PhiDistance;
                    double ParticleDistance;
                    Int_t MomentumClassOfCurrentJet;//flag for the current jet to differ between low-middle-high energetic jets
                    double EtaRecoilSite1;
                    double PhiRecoilSite1;

                    TH2D* h2D_pt_vs_eta = new TH2D("h2D_pt_vs_eta","h2D_pt_vs_eta",i_BinningPerUnit * 2,-1,1,500,0,DEF_HighestExpectableJetMomentum);//temporarily filled
                    TH2D* h2D_pt_vs_phi = new TH2D("h2D_pt_vs_phi","h2D_pt_vs_phi",i_BinningPerUnit * 2 * Pi,-Pi,Pi,500,0,DEF_HighestExpectableJetMomentum);//temporarily filled
                    TH2D* h2D_pt_vs_eta_recoil = new TH2D("h2D_pt_vs_eta_recoil","h2D_pt_vs_eta_recoil",i_BinningPerUnit * 2,-1,1,DEF_MomentumBinning,0,DEF_HighestExpectableRecoilMomentum);//temporarily filled
                    TH2D* h2D_pt_vs_phi_recoil = new TH2D("h2D_pt_vs_phi_recoil","h2D_pt_vs_phi_recoil",i_BinningPerUnit * 2 * Pi,-Pi,Pi,DEF_MomentumBinning,0,DEF_HighestExpectableRecoilMomentum);//temporarily filled

                    
                    //get the momentum class of the trigger particle
                    if(DEF_LowMomentumJetLimit < JetVector[0].pt() && JetVector[0].pt() <= DEF_MiddleMomentumJetLimit) MomentumClassOfCurrentJet = DEF_LowMomentumJetLimit;
                    else if(DEF_MiddleMomentumJetLimit < JetVector[0].pt() && JetVector[0].pt() <= DEF_HighMomentumJetLimit) MomentumClassOfCurrentJet = DEF_MiddleMomentumJetLimit;
                    else if(DEF_HighMomentumJetLimit < JetVector[0].pt()) MomentumClassOfCurrentJet = DEF_HighMomentumJetLimit;

                    //get coordinates of the direct recoil site
                    #ifdef DEF_GetCrosscheckInsteadOfRecoil
                        EtaRecoilSite1 = (JetVector[0].eta() + (0.9 - DEF_JetRadius));
                        if(EtaRecoilSite1 > (0.9 - DEF_JetRadius)) EtaRecoilSite1 -= 2*(0.9 - DEF_JetRadius);//normalize
                        (v_HighPtCoordinates[0].phi_std() > 0) ? (PhiRecoilSite1 = v_HighPtCoordinates[0].phi_std() - Pi/2) : (PhiRecoilSite1 = v_HighPtCoordinates[0].phi_std() + Pi/2);
                    #else    
                    EtaRecoilSite1 = JetVector[0].eta();
                    (JetVector[0].phi_std() > 0) ? (PhiRecoilSite1 = JetVector[0].phi_std() - Pi) : (PhiRecoilSite1 = JetVector[0].phi_std() + Pi);
                    #endif

                    //queue vector for background capture
                    v_JetWaitingForBackgroundCapture.push_back(JetVector[0]);

                    for(Int_t i_Particle = 0; i_Particle < ParticleMemory.size(); i_Particle++)
                    {

                        //skip if particle is the jet particle itself
                        if(JetVector[0].eta() == ParticleMemory[i_Particle][1] && JetVector[0].phi_std() == ParticleMemory[i_Particle][0]) continue;

                        //JET SITE
                        EtaDistance = abs(JetVector[0].eta() - ParticleMemory[i_Particle][1]);
                        PhiDistance = abs(JetVector[0].phi_std() - ParticleMemory[i_Particle][0]);
                        if(PhiDistance > Pi) PhiDistance -= 2* Pi;
                        ParticleDistance = TMath::Sqrt(EtaDistance*EtaDistance + PhiDistance*PhiDistance);

                        if(ParticleDistance < DEF_JetRadius)
                        {
                            //get data in eta direction
                            h2D_pt_vs_eta->Fill(ParticleMemory[i_Particle][1]-JetVector[0].eta(), ParticleMemory[i_Particle][2]);
                            //get data in phi direction
                            h2D_pt_vs_phi->Fill(ParticleMemory[i_Particle][0]-JetVector[0].phi_std(), ParticleMemory[i_Particle][2]);

                            continue;
                        }
                            
                        //RECOIL SITE
                        EtaDistance = abs(EtaRecoilSite1 - ParticleMemory[i_Particle][1]);
                        PhiDistance = abs(PhiRecoilSite1 - ParticleMemory[i_Particle][0]);
                        if(PhiDistance > Pi) PhiDistance -= 2* Pi;
                        ParticleDistance = TMath::Sqrt(EtaDistance*EtaDistance + PhiDistance*PhiDistance);

                        if(ParticleDistance < DEF_JetRadius)
                        {
                            //get data in eta direction
                            h2D_pt_vs_eta_recoil->Fill(ParticleMemory[i_Particle][1]-EtaRecoilSite1, ParticleMemory[i_Particle][2]);
                            //get data in phi direction
                            h2D_pt_vs_phi_recoil->Fill(ParticleMemory[i_Particle][0]-PhiRecoilSite1, ParticleMemory[i_Particle][2]);
                        }
                    }

                    //write into the correct histogram
                    if(MomentumClassOfCurrentJet == DEF_LowMomentumJetLimit)
                    {
                        //JET
                        h2D_pt_vs_eta_Summed_LowPt->Add(h2D_pt_vs_eta);
                        h2D_pt_vs_phi_Summed_LowPt->Add(h2D_pt_vs_phi);
                        //RECOIL
                        h2D_pt_vs_eta_recoil_Summed_LowPt->Add(h2D_pt_vs_eta_recoil);
                        h2D_pt_vs_phi_recoil_Summed_LowPt->Add(h2D_pt_vs_phi_recoil);

                    }
                    else if(MomentumClassOfCurrentJet == DEF_MiddleMomentumJetLimit)
                    {

                        h2D_pt_vs_eta_Summed_MiddlePt->Add(h2D_pt_vs_eta);
                        h2D_pt_vs_phi_Summed_MiddlePt->Add(h2D_pt_vs_phi);

                        h2D_pt_vs_eta_recoil_Summed_MiddlePt->Add(h2D_pt_vs_eta_recoil);
                        h2D_pt_vs_phi_recoil_Summed_MiddlePt->Add(h2D_pt_vs_phi_recoil);
                    }
                    else if(MomentumClassOfCurrentJet == DEF_HighMomentumJetLimit)
                    {

                        h2D_pt_vs_eta_Summed_HighPt->Add(h2D_pt_vs_eta);
                        h2D_pt_vs_phi_Summed_HighPt->Add(h2D_pt_vs_phi);

                        h2D_pt_vs_eta_recoil_Summed_HighPt->Add(h2D_pt_vs_eta_recoil);
                        h2D_pt_vs_phi_recoil_Summed_HighPt->Add(h2D_pt_vs_phi_recoil);

                    }

                    //delte histograms
                    delete h2D_pt_vs_eta;
                    delete h2D_pt_vs_phi;
                    delete h2D_pt_vs_eta_recoil;
                    delete h2D_pt_vs_phi_recoil;

                }

            #endif

            #ifdef DEF_Output_h2D_JetsDijets

                //take data of jet event
                if(m_SuitableEvent)
                {

                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    /*
                    /
                    /               GENERATE HISTOGRAMS OF EACH SINGLE/DIJET EVENT TO GET A VISUAL IDEA OF WHAT IS HAPPENING
                    /
                    /
                    /
                    /
                    *//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                    //generate new histogram
                    TH2D* h2D_phi_vs_eta = new TH2D("h2D_phi_vs_eta","h2D_phi_vs_eta",i_BinningPerUnit * 2,-1,1,i_BinningPerUnit * 2 * Pi,-Pi,Pi);//temporarily filled

                    // loop over tracks within the event 
                    for(Int_t i_Particle = 0; i_Particle < ParticleMemory.size(); i_Particle++)
                    {
                        h2D_phi_vs_eta->Fill(ParticleMemory[i_Particle][1], ParticleMemory[i_Particle][0],ParticleMemory[i_Particle][2]);
                        
                    }

                    //axes
                    h2D_phi_vs_eta->GetXaxis()->SetTitle("eta");
                    h2D_phi_vs_eta->GetYaxis()->SetTitle("phi[rad]"); 

                    //write output file
                    JetsDir->cd();
                    h2D_phi_vs_eta->Write();

                    //append histogram to dijet histogram vector if it is a dijet
                    if(JetVector.size() > 1)
                    {
                        DijetsDir->cd();
                        h2D_phi_vs_eta->Write();
                    }

                    if(m_SuitableEvent)
                    {
                        SuitableDir->cd();
                        h2D_phi_vs_eta->Write();
                    }

                    delete h2D_phi_vs_eta;

                }

            #endif

            #ifdef DEF_Output_ShowJetFinderResults

                
                

            #endif

        }   // end of event loop 

    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /*
    /
    /                                   OUTPUT THE RESULTS OF THE JET SPATIAL DISTRIBUTION PLOTS 
    /
    /
    /
    /
    /
    *//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    #ifdef DEF_Output_h2D_JetsAvgSpatialDistribution

        //mark axes
        h2D_pt_vs_eta_Summed_LowPt->GetXaxis()->SetTitle("Delta eta");
        h2D_pt_vs_eta_Summed_LowPt->GetYaxis()->SetTitle("pt[GeV/c]");
        h2D_pt_vs_phi_Summed_LowPt->GetXaxis()->SetTitle("Delta phi[rad]");
        h2D_pt_vs_phi_Summed_LowPt->GetYaxis()->SetTitle("pt[GeV/c]");

        h2D_pt_vs_eta_Summed_MiddlePt->GetXaxis()->SetTitle("Delta eta");
        h2D_pt_vs_eta_Summed_MiddlePt->GetYaxis()->SetTitle("pt[GeV/c]");
        h2D_pt_vs_phi_Summed_MiddlePt->GetXaxis()->SetTitle("Delta phi[rad]");
        h2D_pt_vs_phi_Summed_MiddlePt->GetYaxis()->SetTitle("pt[GeV/c]");

        h2D_pt_vs_eta_Summed_HighPt->GetXaxis()->SetTitle("Delta eta");
        h2D_pt_vs_eta_Summed_HighPt->GetYaxis()->SetTitle("pt[GeV/c]");
        h2D_pt_vs_phi_Summed_HighPt->GetXaxis()->SetTitle("Delta phi[rad]");
        h2D_pt_vs_phi_Summed_HighPt->GetYaxis()->SetTitle("pt[GeV/c]");

        h2D_pt_vs_eta_recoil_Summed_LowPt->GetXaxis()->SetTitle("Delta eta");
        h2D_pt_vs_eta_recoil_Summed_LowPt->GetYaxis()->SetTitle("pt[GeV/c]");
        h2D_pt_vs_phi_recoil_Summed_LowPt->GetXaxis()->SetTitle("Delta phi[rad]");
        h2D_pt_vs_phi_recoil_Summed_LowPt->GetYaxis()->SetTitle("pt[GeV/c]");

        h2D_pt_vs_eta_recoil_Summed_MiddlePt->GetXaxis()->SetTitle("Delta eta");
        h2D_pt_vs_eta_recoil_Summed_MiddlePt->GetYaxis()->SetTitle("pt[GeV/c]");
        h2D_pt_vs_phi_recoil_Summed_MiddlePt->GetXaxis()->SetTitle("Delta phi[rad]");
        h2D_pt_vs_phi_recoil_Summed_MiddlePt->GetYaxis()->SetTitle("pt[GeV/c]");

        h2D_pt_vs_eta_recoil_Summed_HighPt->GetXaxis()->SetTitle("Delta eta");
        h2D_pt_vs_eta_recoil_Summed_HighPt->GetYaxis()->SetTitle("pt[GeV/c]");
        h2D_pt_vs_phi_recoil_Summed_HighPt->GetXaxis()->SetTitle("Delta phi[rad]");
        h2D_pt_vs_phi_recoil_Summed_HighPt->GetYaxis()->SetTitle("pt[GeV/c]");

        h2D_pt_vs_eta_Background_Summed->GetXaxis()->SetTitle("eta");
        h2D_pt_vs_eta_Background_Summed->GetYaxis()->SetTitle("pt[GeV/c]");
        h2D_pt_vs_phi_Background_Summed->GetXaxis()->SetTitle("phi[rad]");
        h2D_pt_vs_phi_Background_Summed->GetYaxis()->SetTitle("pt[GeV/c]");

        //normalize plots

        //write plots
        PtVsEtaDir->cd();
        h2D_pt_vs_eta_Summed_LowPt->Write();
        h2D_pt_vs_eta_Summed_MiddlePt->Write();
        h2D_pt_vs_eta_Summed_HighPt->Write();
        h2D_pt_vs_eta_recoil_Summed_LowPt->Write();
        h2D_pt_vs_eta_recoil_Summed_MiddlePt->Write();
        h2D_pt_vs_eta_recoil_Summed_HighPt->Write();
        h2D_pt_vs_eta_Background_Summed->Write();

        PtVsPhiDir->cd();
        h2D_pt_vs_phi_Summed_LowPt->Write();
        h2D_pt_vs_phi_Summed_MiddlePt->Write();
        h2D_pt_vs_phi_Summed_HighPt->Write();
        h2D_pt_vs_phi_recoil_Summed_LowPt->Write();
        h2D_pt_vs_phi_recoil_Summed_MiddlePt->Write();
        h2D_pt_vs_phi_recoil_Summed_HighPt->Write();
        h2D_pt_vs_phi_Background_Summed->Write();

        //make sure this is small to ensure reliable results. Depends on how many jet-free areas could be found
        cout << "In the end " << v_JetWaitingForBackgroundCapture.size() << " recoil distributions were not compensated for background" << endl;


    #endif

    OutputFile ->Close();

    cout << "Logged Jets: " << DEBUG_NoOfLoggedJets << endl << "Logged Bg: " << DEBUG_NoOfLoggedBg << endl;

    cout << "#Jet events: " << DEBUG_NoOfJetEvents << endl;


    printf("StJetAnalysis::ReadTrees() finished \n");

}


void AnaRecoil(const TString DataFile, Int_t i_ProjectionBins, Int_t i_BinLow, Int_t i_BinMid, Int_t i_BinHigh, const TString i_PrintDiff) {

    #define i_StartIntegral -0.5 
    #define i_StopIntegral -0.1

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
    TH2D* h2D_pt_vs_eta_Background = (TH2D*)file->Get("Pt_vs_Eta/h2D_pt_vs_eta_Background_Summed");
    if(!h2D_pt_vs_eta_Background) cout << "Background histogram not found!" << endl;
    TH2D* h2D_pt_vs_eta_recoil_LowPt = (TH2D*)file->Get("Pt_vs_Eta/h2D_pt_vs_eta_recoil_Summed_LowPt");
    if(!h2D_pt_vs_eta_recoil_LowPt) cout << "LowPt histogram not found!" << endl;
    TH2D* h2D_pt_vs_eta_recoil_MiddlePt = (TH2D*)file->Get("Pt_vs_Eta/h2D_pt_vs_eta_recoil_Summed_MiddlePt");
    if(!h2D_pt_vs_eta_recoil_MiddlePt) cout << "MiddlePt histogram not found!" << endl;
    TH2D* h2D_pt_vs_eta_recoil_HighPt = (TH2D*)file->Get("Pt_vs_Eta/h2D_pt_vs_eta_recoil_Summed_HighPt");
    if(!h2D_pt_vs_eta_recoil_HighPt) cout << "HighPt histogram not found!" << endl;

    //get Y binning
    double MomentumPerBin = (h2D_pt_vs_eta_Background->GetYaxis()->GetXmax() /h2D_pt_vs_eta_Background->GetNbinsY());

    //configure and draw canvas
    TCanvas* can_JetRecoilSitesAndBg = new TCanvas("JetRecoilSitesAndBg","JetRecoilSitesAndBg",1000,1000);
    can_JetRecoilSitesAndBg->Divide(3,2);

    //loop through three horizontal cuts with different binning
    for(int i=1; i<=3; i++)
    {
        Int_t FirstBin;
        Int_t LastBin;

        if(i==1) FirstBin = i_BinLow;
        else if(i==2) FirstBin = i_BinMid;
        else if(i==3) FirstBin = i_BinHigh;

        LastBin = FirstBin + i_ProjectionBins;

        //generate projections
        TH1D* h1D_pt_vs_eta_Background_Projection = h2D_pt_vs_eta_Background->ProjectionX("h1D_pt_vs_eta_Background_Projection", FirstBin, LastBin);
        TH1D* h1D_pt_vs_eta_recoil_LowPt_Projection = h2D_pt_vs_eta_recoil_LowPt->ProjectionX("h1D_pt_vs_eta_recoil_LowPt_Projection", FirstBin, LastBin);
        TH1D* h1D_pt_vs_eta_recoil_MiddlePt_Projection = h2D_pt_vs_eta_recoil_MiddlePt->ProjectionX("h1D_pt_vs_eta_recoil_MiddlePt_Projection", FirstBin, LastBin);
        TH1D* h1D_pt_vs_eta_recoil_HighPt_Projection = h2D_pt_vs_eta_recoil_HighPt->ProjectionX("h1D_pt_vs_eta_recoil_HighPt_Projection", FirstBin, LastBin);

        //get bins
        Int_t NoOfXBins = h1D_pt_vs_eta_Background_Projection->GetNbinsX();
        double BinsPerEta = (double)(NoOfXBins)/(1.8);
        Int_t Int1_StartBin = (Int_t)(BinsPerEta * (i_StartIntegral + 0.9));
        Int_t Int1_EndBin = (Int_t)(BinsPerEta * (i_StopIntegral + 0.9));
        Int_t Int2_StartBin = (Int_t)(NoOfXBins - Int1_EndBin);
        Int_t Int2_EndBin = (Int_t)(NoOfXBins - Int1_StartBin);

        //integrate 
        double BackgroundIntegral = h1D_pt_vs_eta_Background_Projection->Integral(Int1_StartBin, Int1_EndBin) + h1D_pt_vs_eta_Background_Projection->Integral(Int2_StartBin, Int2_EndBin);

        double LowIntegral = h1D_pt_vs_eta_recoil_LowPt_Projection->Integral(Int1_StartBin, Int1_EndBin) + h1D_pt_vs_eta_recoil_LowPt_Projection->Integral(Int2_StartBin, Int2_EndBin);
        double MiddleIntegral = h1D_pt_vs_eta_recoil_MiddlePt_Projection->Integral(Int1_StartBin, Int1_EndBin) + h1D_pt_vs_eta_recoil_MiddlePt_Projection->Integral(Int2_StartBin, Int2_EndBin);
        double HighIntegral = h1D_pt_vs_eta_recoil_HighPt_Projection->Integral(Int1_StartBin, Int1_EndBin) + h1D_pt_vs_eta_recoil_HighPt_Projection->Integral(Int2_StartBin, Int2_EndBin);

        h1D_pt_vs_eta_recoil_LowPt_Projection->Scale(BackgroundIntegral/LowIntegral);
        h1D_pt_vs_eta_recoil_MiddlePt_Projection->Scale(BackgroundIntegral/MiddleIntegral);
        h1D_pt_vs_eta_recoil_HighPt_Projection->Scale(BackgroundIntegral/HighIntegral);

        can_JetRecoilSitesAndBg->cd(i);

        //axes
        if(i==1)
        {
            h1D_pt_vs_eta_Background_Projection->GetYaxis()->SetTitle("Multiplicity");
        }

        else
        {
            h1D_pt_vs_eta_Background_Projection->GetYaxis()->SetTitle("");
            h1D_pt_vs_eta_Background_Projection->GetYaxis()->SetTitle("");
        }
        h1D_pt_vs_eta_Background_Projection->GetXaxis()->SetTitle("#Delta #eta");

        double LowestMomentum = MomentumPerBin * FirstBin;
        double HighestMomentum = MomentumPerBin * (FirstBin + i_ProjectionBins);

        //histograms
        h1D_pt_vs_eta_Background_Projection->SetTitle(Form("Multiplicity for particles with p #in [%.3f,%.3f] on off-side", LowestMomentum, HighestMomentum));
        h1D_pt_vs_eta_Background_Projection->SetMarkerStyle(kFullSquare);
        h1D_pt_vs_eta_Background_Projection->SetNdivisions(5);
        TH1* hbg = (TH1*)h1D_pt_vs_eta_Background_Projection->DrawCopy("HIST P");
        h1D_pt_vs_eta_recoil_LowPt_Projection->SetTitle("LowPt(30...50GeV/c) Distribution");
        h1D_pt_vs_eta_recoil_LowPt_Projection->SetMarkerStyle(kFullCircle);
        h1D_pt_vs_eta_recoil_LowPt_Projection->SetMarkerColor(kRed);
        TH1* hlow = (TH1*)h1D_pt_vs_eta_recoil_LowPt_Projection->DrawCopy("HIST P SAME");
        h1D_pt_vs_eta_recoil_MiddlePt_Projection->SetTitle("MiddlePt(50...80GeV/c) Distribution");
        h1D_pt_vs_eta_recoil_MiddlePt_Projection->SetMarkerStyle(kPlus);
        h1D_pt_vs_eta_recoil_MiddlePt_Projection->SetMarkerColor(kBlue);
        TH1* hmid = (TH1*)h1D_pt_vs_eta_recoil_MiddlePt_Projection->DrawCopy("HIST P SAME");
        h1D_pt_vs_eta_recoil_HighPt_Projection->SetTitle("HighPt(>80GeV/c) Distribution");
        h1D_pt_vs_eta_recoil_HighPt_Projection->SetMarkerStyle(kStar);
        h1D_pt_vs_eta_recoil_HighPt_Projection->SetMarkerColor(kGreen);
        TH1* hhigh = (TH1*)h1D_pt_vs_eta_recoil_HighPt_Projection->DrawCopy("HIST P SAME");
        //draw legend
        TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
        leg->AddEntry(hbg,"Background","p");
        leg->AddEntry(hlow,"LowPt(30...50GeV/c) jets","p");
        leg->AddEntry(hmid,"MiddlePt(50...80GeV/c) jets","p");
        leg->AddEntry(hhigh, "HighPt(>80GeV/c) jets","p");
        leg->Draw();

        can_JetRecoilSitesAndBg->cd(i+3);

        //get the differences
        TH1D* h1D_MultDif_vs_eta_LowPt = (TH1D*)h1D_pt_vs_eta_recoil_LowPt_Projection->Clone("h1D_MultDif_vs_eta_LowPt");
        h1D_MultDif_vs_eta_LowPt->Add(h1D_pt_vs_eta_Background_Projection, -1);
        TH1D* h1D_MultDif_vs_eta_MiddlePt = (TH1D*)h1D_pt_vs_eta_recoil_MiddlePt_Projection->Clone("h1D_MultDif_vs_eta_MiddlePt");
        h1D_MultDif_vs_eta_MiddlePt->Add(h1D_pt_vs_eta_Background_Projection, -1);
        TH1D* h1D_MultDif_vs_eta_HighPt = (TH1D*)h1D_pt_vs_eta_recoil_HighPt_Projection->Clone("h1D_MultDif_vs_eta_HighPt");
        h1D_MultDif_vs_eta_HighPt->Add(h1D_pt_vs_eta_Background_Projection, -1);

        if(i_PrintDiff == "HIGH"){
        h1D_MultDif_vs_eta_HighPt->SetNdivisions(5);
        h1D_MultDif_vs_eta_HighPt->GetXaxis()->SetTitle("#Delta #eta");
        if(i==1)
        {
            h1D_MultDif_vs_eta_HighPt->GetYaxis()->SetTitle("#Delta Multiplicity");
        }
        else
        {
            h1D_MultDif_vs_eta_HighPt->GetYaxis()->SetTitle("");
        }
        h1D_MultDif_vs_eta_HighPt->SetTitle("DIfference in Multiplicity behind HighPt jets relative to BG");
        h1D_MultDif_vs_eta_HighPt->Rebin(8);
        h1D_MultDif_vs_eta_HighPt->DrawCopy("E P");
        }
        else if(i_PrintDiff == "MIDDLE"){
        h1D_MultDif_vs_eta_MiddlePt->SetNdivisions(5);
        h1D_MultDif_vs_eta_MiddlePt->GetXaxis()->SetTitle("#Delta #eta");
        if(i==1)
        {
            h1D_MultDif_vs_eta_MiddlePt->GetYaxis()->SetTitle("#Delta Multiplicity");
        }
        else
        {
            h1D_MultDif_vs_eta_MiddlePt->GetYaxis()->SetTitle("");
        }
        h1D_MultDif_vs_eta_MiddlePt->SetTitle("DIfference in Multiplicity behind MiddlePt jets relative to BG");
        h1D_MultDif_vs_eta_MiddlePt->Rebin(8);
        h1D_MultDif_vs_eta_MiddlePt->DrawCopy("E P");
        }
        else if(i_PrintDiff == "LOW"){
        h1D_MultDif_vs_eta_LowPt->SetNdivisions(5);
        h1D_MultDif_vs_eta_LowPt->GetXaxis()->SetTitle("#Delta #eta");
        if(i==1)
        {
            h1D_MultDif_vs_eta_LowPt->GetYaxis()->SetTitle("#Delta Multiplicity");
        }
        else
        {
            h1D_MultDif_vs_eta_LowPt->GetYaxis()->SetTitle("");
        }
        h1D_MultDif_vs_eta_LowPt->SetTitle("DIfference in Multiplicity behind LowPt jets relative to BG");
        h1D_MultDif_vs_eta_LowPt->Rebin(8);
        h1D_MultDif_vs_eta_LowPt->DrawCopy("E P");
        }

    }

    cout << "DONE!" << endl;
    file->Close();
}


void AnaRecoilALL(const TString DataFile, Int_t i_ProjectionBins, Int_t i_BinLow, Int_t i_BinMid, Int_t i_BinHigh) {

    #define DEF_AxisLabelSize 0.05
    #define DEF_HistoTitleSize 0.1
    #define DEF_InfoBoxX 0.95
    #define DEF_InfoBoxY 0.9
    #define DEF_LastColumnInfoXOffset -0.15
    #define DEF_InfoBoxLineDist 0.1
    #define DEF_LatexTextSize 0.04

    #define i_StartIntegral -0.5 
    #define i_StopIntegral -0.2

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
    TH2D* h2D_pt_vs_eta_Background = (TH2D*)file->Get("Pt_vs_Eta/h2D_pt_vs_eta_Background_Summed");
    if(!h2D_pt_vs_eta_Background) cout << "Background histogram not found!" << endl;
    TH2D* h2D_pt_vs_eta_recoil_LowPt = (TH2D*)file->Get("Pt_vs_Eta/h2D_pt_vs_eta_recoil_Summed_LowPt");
    if(!h2D_pt_vs_eta_recoil_LowPt) cout << "LowPt histogram not found!" << endl;
    TH2D* h2D_pt_vs_eta_recoil_MiddlePt = (TH2D*)file->Get("Pt_vs_Eta/h2D_pt_vs_eta_recoil_Summed_MiddlePt");
    if(!h2D_pt_vs_eta_recoil_MiddlePt) cout << "MiddlePt histogram not found!" << endl;
    TH2D* h2D_pt_vs_eta_recoil_HighPt = (TH2D*)file->Get("Pt_vs_Eta/h2D_pt_vs_eta_recoil_Summed_HighPt");
    if(!h2D_pt_vs_eta_recoil_HighPt) cout << "HighPt histogram not found!" << endl;

    //get Y binning
    double MomentumPerBin = (h2D_pt_vs_eta_Background->GetYaxis()->GetXmax() /h2D_pt_vs_eta_Background->GetNbinsY());

    //configure and draw canvas
    TCanvas* can_JetRecoilSitesAndBg = new TCanvas("JetRecoilSitesAndBg","JetRecoilSitesAndBg",1000,1000);
    can_JetRecoilSitesAndBg->Divide(3,2);

    //add all distributions to a new histogram
    TH2D* h2D_pt_vs_eta_recoil_All = (TH2D*)h2D_pt_vs_eta_recoil_LowPt->Clone();
    h2D_pt_vs_eta_recoil_All->Add(h2D_pt_vs_eta_recoil_MiddlePt);
    h2D_pt_vs_eta_recoil_All->Add(h2D_pt_vs_eta_recoil_HighPt);

    //scale histograms
    h2D_pt_vs_eta_recoil_All->Scale(1./h2D_pt_vs_eta_recoil_All->GetEntries());
    h2D_pt_vs_eta_Background->Scale(1./h2D_pt_vs_eta_Background->GetEntries());

    //loop through three horizontal cuts with different binning
    for(int i=1; i<=3; i++)
    {
        

        Int_t FirstBin;
        Int_t LastBin;

        if(i==1) FirstBin = i_BinLow;
        else if(i==2) FirstBin = i_BinMid;
        else if(i==3) FirstBin = i_BinHigh;

        LastBin = FirstBin + i_ProjectionBins;

        //generate projections
        TH1D* h1D_pt_vs_eta_Background_Projection = h2D_pt_vs_eta_Background->ProjectionX("h1D_pt_vs_eta_Background_Projection", FirstBin, LastBin);
        TH1D* h1D_pt_vs_eta_recoil_All_Projection = h2D_pt_vs_eta_recoil_All->ProjectionX("h1D_pt_vs_eta_recoil_All_Projection", FirstBin, LastBin);

        //get bins
        Int_t NoOfXBins = h1D_pt_vs_eta_Background_Projection->GetNbinsX();
        double BinsPerEta = (double)(NoOfXBins)/(1.8);
        Int_t Int1_StartBin = (Int_t)(BinsPerEta * (i_StartIntegral + 0.9));
        Int_t Int1_EndBin = (Int_t)(BinsPerEta * (i_StopIntegral + 0.9));
        Int_t Int2_StartBin = (Int_t)(NoOfXBins - Int1_EndBin);
        Int_t Int2_EndBin = (Int_t)(NoOfXBins - Int1_StartBin);

        //integrate 
        double BackgroundIntegral = h1D_pt_vs_eta_Background_Projection->Integral(Int1_StartBin, Int1_EndBin) + h1D_pt_vs_eta_Background_Projection->Integral(Int2_StartBin, Int2_EndBin);

        double Integral = h1D_pt_vs_eta_recoil_All_Projection->Integral(Int1_StartBin, Int1_EndBin) + h1D_pt_vs_eta_recoil_All_Projection->Integral(Int2_StartBin, Int2_EndBin);

        h1D_pt_vs_eta_recoil_All_Projection->Scale(BackgroundIntegral/Integral);

        can_JetRecoilSitesAndBg->cd(i);

        gPad->SetMargin(0.16*(i==1), 0.15*(i==3), 0.15, 0.02);

        //markings
        h1D_pt_vs_eta_Background_Projection->SetTitle("");
        h1D_pt_vs_eta_Background_Projection->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
        if(i==3) h1D_pt_vs_eta_Background_Projection->GetXaxis()->SetTitle("#Delta #eta");
        else h1D_pt_vs_eta_Background_Projection->GetXaxis()->SetTitle("");
        if(i==1)h1D_pt_vs_eta_Background_Projection->GetYaxis()->SetTitle("Multiplicity");
        else h1D_pt_vs_eta_Background_Projection->GetYaxis()->SetTitle("");
        

        double LowestMomentum = MomentumPerBin * FirstBin;
        double HighestMomentum = MomentumPerBin * (FirstBin + i_ProjectionBins);

        //histograms
        h1D_pt_vs_eta_Background_Projection->SetMarkerStyle(kFullSquare);
        h1D_pt_vs_eta_Background_Projection->SetNdivisions(5);
        TH1* hbg = (TH1*)h1D_pt_vs_eta_Background_Projection->DrawCopy("HIST P");
        h1D_pt_vs_eta_recoil_All_Projection->SetTitle("All jets Distribution");
        h1D_pt_vs_eta_recoil_All_Projection->SetMarkerStyle(kFullCircle);
        h1D_pt_vs_eta_recoil_All_Projection->SetMarkerColor(kRed);
        TH1* hall = (TH1*)h1D_pt_vs_eta_recoil_All_Projection->DrawCopy("HIST P SAME");
        //draw legend
        TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
        leg->AddEntry(hbg,"Background","p");
        leg->AddEntry(hall,"All jets","p");
        leg->Draw();

        can_JetRecoilSitesAndBg->cd(i+3);

        //get the differences
        TH1D* h1D_MultDif_vs_eta_All = (TH1D*)h1D_pt_vs_eta_recoil_All_Projection->Clone("h1D_MultDif_vs_eta_All");
        h1D_MultDif_vs_eta_All->Add(h1D_pt_vs_eta_Background_Projection, -1);

        h1D_MultDif_vs_eta_All->SetNdivisions(5);
        h1D_MultDif_vs_eta_All->GetXaxis()->SetTitle("#Delta #eta");
        if(i==1)
        {
            h1D_MultDif_vs_eta_All->GetYaxis()->SetTitle("#Delta Multiplicity");
        }
        else
        {
            h1D_MultDif_vs_eta_All->GetYaxis()->SetTitle("");
        }
        h1D_MultDif_vs_eta_All->SetTitle("DIfference in Multiplicity behind all jets relative to BG");
        h1D_MultDif_vs_eta_All->DrawCopy("E P");

    }

    cout << "DONE!" << endl;
    file->Close();
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
    #define DEF_MinJetAcceptancePt 20

    // define variables 
    vector<PseudoJet> jets_fiducial;//contains the selected jets in the end

    Double_t ghost_maxrap = 1.1; // Fiducial cut for background estimation
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

    Int_t Rem_n_hardest = 1; 
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
        if(jet_area < Jet_area_cut) continue; 
        if(jet_pt_sub < DEF_MinJetAcceptancePt) continue;

        AcceptedJets.push_back(jets_fiducial[i_jet]);

        // constituent track loop
        vector<PseudoJet> jet_constituents = jets_fiducial[i_jet].constituents();

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
    }

    return AcceptedJets;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
/
/                                                FUNCTIONS FOR BACKGROUND GENERATION 
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