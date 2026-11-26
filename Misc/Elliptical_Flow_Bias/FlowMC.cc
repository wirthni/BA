/*
Author: Nicolas Wirth, Nicolas.Wirth@cern.ch
*/

#include "./StJetTrackEvent.h"
#include "./StJetTrackEventLinkDef.h"
#include "./functions.h"
#include "TClassRef.h"
#include "Math/Vector4D.h"
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include "Fastjet/fastjet-3.5.1/include/fastjet/config.h"
#include "Fastjet/fastjet-3.5.1/include/fastjet/ClusterSequenceArea.hh"
#include "Fastjet/fastjet-3.5.1/include/fastjet/PseudoJet.hh"
#include "Fastjet/fastjet-3.5.1/include/fastjet/Selector.hh"
#include "Fastjet/fastjet-3.5.1/include/fastjet/tools/Subtractor.hh"
#include "Fastjet/fastjet-3.5.1/include/fastjet/tools/JetMedianBackgroundEstimator.hh"
using namespace fastjet;
using namespace std;

vector<PseudoJet> FindJets(vector<PseudoJet> vec_particles);
double GetPhiByProb(void);
double GetEtaByProb(void);
double GetPtByProb(void);
double GetNByProb(void);

double GetRandomF(double (*i_FuncPtr)(double, double[]), double i_LowerLim, double i_UpperLim, double i_FuncParams[]);
double Function_NormGauss(double x, double params[]);
double Function_PhiByFlow(double x, double params[]);
double Function_v2v3ByPtAndEta(double x, double params[]);
double Function_DirectedFlowByEta(double Eta);
template <typename T> int sgn(T val);

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
*/
Int_t FlowMC(int i_NoOfEvents = -1, bool i_AllowV2 = false, bool i_AllowV3 = false, TString i_OutputFileName = "Default", int i_InstanceNo = 1)
{

    gStyle->SetOptStat(0);
    SetRootGraphicStyle();

    long NoOfPosPulls = 0;
    long NoOfNegPulls = 0;

    #define DEF_BinningPerUnit 100
    #define DEF_JetRadius 0.2
    #define DEF_OutputEventOverviews false
    #define DEF_HadLeadingPt 15.0
    #define DEF_HadSubleadingPt 8.0
    #define DEF_CorrelationMinPt 0.0
    #define DEF_CorrelationMaxPt 1.0
    #define DEF_MaxPartilclesPerJet 200
    #define DEF_MaxJetPt 150.0
    #define DEF_AxisLabelSize 0.05
    #define DEF_HistoTitleSize 0.1
    #define DEF_BackToBackSmearingSigma 50

    //other variables
    TH2D* h2D_pt_vs_eta_LargeGap = new TH2D("h2D_pt_vs_eta_LargeGap","h2D_pt_vs_eta_LargeGap",100, 0, 10, DEF_BinningPerUnit * 2 * 0.9,-0.9,0.9);
    TH2D* h2D_pt_vs_eta_SmallGap = new TH2D("h2D_pt_vs_eta_SmallGap","h2D_pt_vs_eta_SmallGap",100, 0, 10, DEF_BinningPerUnit * 2 * 0.9,-0.9,0.9);
    TH2D* h2D_eta_vs_dphi_LargeGap = new TH2D("h2D_eta_vs_dphi_LargeGap","h2D_eta_vs_dphi_LargeGap",DEF_BinningPerUnit * 2 * 0.9, -0.9, 0.9, DEF_BinningPerUnit * 2 * Pi, -Pi/2, 3*Pi/2);
    TH2D* h2D_eta_vs_dphi_SmallGap = new TH2D("h2D_eta_vs_dphi_SmallGap","h2D_eta_vs_dphi_SmallGap",DEF_BinningPerUnit * 2 * 0.9, -0.9, 0.9, DEF_BinningPerUnit * 2 * Pi, -Pi/2, 3*Pi/2);
    TH2D* h2D_deta_vs_dphi_LargeGap = new TH2D("h2D_deta_vs_dphi_LargeGap","h2D_deta_vs_dphi_LargeGap",DEF_BinningPerUnit * 4 * 0.9, -1.8, 1.8, DEF_BinningPerUnit * 2 * Pi, -Pi/2, 3*Pi/2); // contains particle correlations relative to the leading jet
    TH2D* h2D_deta_vs_dphi_SmallGap = new TH2D("h2D_deta_vs_dphi_SmallGap","h2D_deta_vs_dphi_SmallGap",DEF_BinningPerUnit * 4 * 0.9, -1.8, 1.8, DEF_BinningPerUnit * 2 * Pi, -Pi/2, 3*Pi/2); // contains particle correlations relative to the leading jet
    TH1D* h1D_JetPt_Leading = new TH1D("h1D_JetPt_Leading", "h1D_JetPt_Leading", DEF_MaxJetPt, 0, DEF_MaxJetPt);
    TH1D* h1D_JetPt_Subleading = new TH1D("h1D_JetPt_Subleading", "h1D_JetPt_Subleading", DEF_MaxJetPt, 0, DEF_MaxJetPt);
    TH2D* h2D_eta_vs_phi_JetCoordinates_Leading_NoCut = new TH2D("h2D_eta_vs_phi_JetCoordinates_Leading_NoCut","h2D_eta_vs_phi_JetCoordinates_Leading_NoCut",DEF_BinningPerUnit * 2 * 0.9, -0.9, 0.9, DEF_BinningPerUnit * 2 * Pi, -Pi, Pi); // contains only the jet coordinates before applying the cut
    TH2D* h2D_eta_vs_phi_JetCoordinates_Leading_OnlyEtaCut = new TH2D("h2D_eta_vs_phi_JetCoordinates_Leading_OnlyEtaCut","h2D_eta_vs_phi_JetCoordinates_Leading_OnlyEtaCut",DEF_BinningPerUnit * 2 * 0.9, -0.9, 0.9, DEF_BinningPerUnit * 2 * Pi, -Pi, Pi); // contains only the jet coordinates before applying the cut
    TH2D* h2D_eta_vs_phi_JetCoordinates_Leading_LargeGapCut = new TH2D("h2D_eta_vs_phi_JetCoordinates_Leading_LargeGapCut","h2D_eta_vs_phi_JetCoordinates_Leading_LargeGapCut",DEF_BinningPerUnit * 2 * 0.9, -0.9, 0.9, DEF_BinningPerUnit * 2 * Pi, -Pi, Pi); // contains only the jet coordinates after applying the large gap cut
    TH2D* h2D_eta_vs_phi_JetCoordinates_Leading_SmallGapCut = new TH2D("h2D_eta_vs_phi_JetCoordinates_Leading_SmallGapCut","h2D_eta_vs_phi_JetCoordinates_Leading_SmallGapCut",DEF_BinningPerUnit * 2 * 0.9, -0.9, 0.9, DEF_BinningPerUnit * 2 * Pi, -Pi, Pi); // contains only the jet coordinates after applying the small gap cut
    TH2D* h2D_GeneratedParticles = new TH2D("h2D_GeneratedParticles","h2D_GeneratedParticles",DEF_BinningPerUnit * 2 * 0.9, -0.9, 0.9, DEF_BinningPerUnit * 2 * Pi, -Pi, Pi); // contains only the jet coordinates before applying the cut
    TH2D* h2D_Eta_vs_v1 = new TH2D("h2D_Eta_vs_v1", "h2D_Eta_vs_v1", 20, -0.9, 0.9, 100, -0.2, 0.2);
    TH1D* h1D_AllJetsEta_PYTHIA = new TH1D("h1D_AllJetsEta_PYTHIA", "h1D_AllJetsEta_PYTHIA", 20, -0.9, 0.9);
    TH1D* h1D_AllJetsEta_HYBRID = new TH1D("h1D_AllJetsEta_HYBRID", "h1D_AllJetsEta_HYBRID", 20, -0.9, 0.9);
    TH2F* h2F_SoftHadronRelEP = new TH2F("h2F_SoftHadronRelEP", "h2F_SoftHadronRelEP", DEF_BinningPerUnit * 2 * Pi, -Pi/2, 3*Pi/2, DEF_BinningPerUnit * 2 * 0.9, -0.9, 0.9);
    TH1F* h1F_PhiLeadMinusPsi_LG = new TH1F("h1F_PhiLeadMinusPsi_LG", "h1F_PhiLeadMinusPsi_LG", DEF_BinningPerUnit * 2 * Pi, -Pi, +Pi);
    TH1F* h1F_PhiLeadMinusPsi_SG = new TH1F("h1F_PhiLeadMinusPsi_SG", "h1F_PhiLeadMinusPsi_SG", DEF_BinningPerUnit * 2 * Pi, -Pi, +Pi);
    

    Long64_t            LargeGapEventCounter = 0;
    Long64_t            SmallGapEventCounter = 0;

    //Get N and p_t statistics
    TFile *file = TFile::Open("merge_mult_track_pT.root");
    if (!file || file->IsZombie()) {
        std::cout << "Error while opening the file!" << std::endl;
        return 0;
    }
    std::cout << "Found N and pT file!" << std::endl;

    //get multiplicity and pt statistics for background
    TH1D* h1D_NDist = (TH1D*)file->Get("h_mult_particles_used_SE_ME_0");
    if(!h1D_NDist){cout << "Multiplicity histogram not found!" << endl; return 0;}
    TH1D* h1D_PtDist = (TH1D*)file->Get("h_particle_pT");
    if(!h1D_PtDist){cout << "Pt histogram not found!" << endl; return 0;}
    std::cout << "Found N and pT statistics!" << std::endl;

    //trigger phase space information
    TFile *PhaseSpaceFile = TFile::Open("M__TriggerCorrelations.root");
    if (!PhaseSpaceFile || PhaseSpaceFile->IsZombie()) {
        std::cout << "Error while opening the phase space file!" << std::endl;
        return 0;
    }
    std::cout << "Found trigger phase space file!" << std::endl;

    THnSparseI* PhaseSpaceHisto = nullptr;
    PhaseSpaceFile->GetObject("TriggerPhaseSpace", PhaseSpaceHisto);

    if(!PhaseSpaceHisto){cout << "Phase space histogram not found!" << endl; return 0;}
    Double_t *TriggerPhaseSpaceTuple = new Double_t[6];
    std::cout << "Found trigger phase space histogram!" << std::endl;

    //ALICE directed flow measurement for background
    TFile *DirectedFlowFile = TFile::Open("v1Oddeta58.root");
    if (!DirectedFlowFile || DirectedFlowFile->IsZombie()) {
        std::cout << "Error while opening the directed flow file!" << std::endl;
        return 0;
    }
    std::cout << "Found directed flow file!" << std::endl;
    TH1D* h1D_v1_pT_Slice[7];
    h1D_v1_pT_Slice[5] = (TH1D*) DirectedFlowFile->Get("hOdd_00-pos08");
    h1D_v1_pT_Slice[1] = (TH1D*) DirectedFlowFile->Get("hOdd_neg08-00");
    h1D_v1_pT_Slice[6] = (TH1D*) DirectedFlowFile->Get("hOdd_pos04-pos08");
    h1D_v1_pT_Slice[0] = (TH1D*) DirectedFlowFile->Get("hOdd_neg08-neg04");
    h1D_v1_pT_Slice[4] = (TH1D*) DirectedFlowFile->Get("hOdd_00-pos04");
    h1D_v1_pT_Slice[2] = (TH1D*) DirectedFlowFile->Get("hOdd_neg04-00");
    //h1D_v1_pT_Slice[3] stays zero
    for(int i=0; i<=6; i++)
    {
        if (i==3) continue;
        if(!h1D_v1_pT_Slice[i]){cout << "Flow histograms not found!" << endl; return 0;}
    }
    int pTBins = h1D_v1_pT_Slice[0]->GetXaxis()->GetNbins();
    std::cout << "Found directed flow statistics!" << std::endl;
    
    //get pT binning
    Double_t pTBinEdges[pTBins + 1];
    for(int iBin=1; iBin <= pTBins; iBin++)
    {
        pTBinEdges[iBin-1] = h1D_v1_pT_Slice[0]->GetXaxis()->GetBinLowEdge(iBin);
    }
    pTBinEdges[pTBins] = h1D_v1_pT_Slice[0]->GetXaxis()->GetBinUpEdge(pTBins);

    TH2D* h2D_v1_pT_Eta = new TH2D("h2D_v1_pT_Eta", "h2D_v1_pT_Eta",
                                    pTBins, pTBinEdges,
                                    7, -.7, .7);

    for (int iEtaBin = 1; iEtaBin <= 7; iEtaBin++)
    {
        if (iEtaBin == 4) continue;
        for(int ipTBin = 1; ipTBin <= pTBins; ipTBin++)
        {
            float Val = h1D_v1_pT_Slice[iEtaBin-1]->GetBinContent(ipTBin);
            if(iEtaBin < 4) Val *= -1.0f;
            h2D_v1_pT_Eta->SetBinContent(ipTBin, iEtaBin, Val);
        }
    }

    TString str_out =  "./" + i_OutputFileName + Form("_%d", i_InstanceNo) + ".root";
    TFile* OutputFile = new TFile(str_out.Data(),"RECREATE");
    OutputFile->cd();
    cout << "Output file is " << str_out << endl;
    TDirectory *EventHistos = OutputFile->mkdir("Event Histos");
    TDirectory *Results = OutputFile->mkdir("Results");

    TRandom3 TR_Eta(0);
    TRandom3 TR_V2Angle(0);
    TRandom3 TR_V3Angle(0);
    TRandom3 TR_N(0);
    TRandom3 TR_pT(0);

    for(Int_t iEvent = 0; iEvent < i_NoOfEvents; iEvent++)
    {
        //event output
        if (iEvent != 0  &&  iEvent % 100 == 0)
            cout << "." << flush;
        if (iEvent != 0  &&  iEvent % 1000 == 0)
        {
            Double_t event_percent = 100.0*((Double_t)(iEvent))/((Double_t)(i_NoOfEvents));
            cout << " " << iEvent << " (" << event_percent << "%) " << "\n" << "==> Processing data " << flush;
        }
                
        float Q_x = 0;
        float Q_y = 0;
        float Psi = 0;
        float DistanceToClosestV2MaximumBeforeBG;
        float DistanceToClosestV2MaximumAfterBG;
        float DistanceToClosestV3MaximumBeforeBG;
        float DistanceToClosestV3MaximumAfterBG;

        
        Int_t N = h1D_NDist->GetRandom(&TR_N);
        double V2_Angle = TR_V2Angle.Uniform(-Pi, Pi);
        double V3_Angle = TR_V3Angle.Uniform(-Pi, Pi);
        vector<array<double, 4>> ParticleMemory;//remembers randomly generated particle coordinates (phi, eta, pt) for relevant cuts later in the analysis
        vector<array<double, 4>> HighPtParticles;//save for later in case of hadron analysis. All particles that have pt > DEF_SubleadingPt
        vector<PseudoJet> ParticleVector;

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /*
        /
        /                                                     
        /                                                       GENERATE TWO JET TRIGGERS
        /                                                     
        /                                                     
        /                                                     
        /
        /                                                    
        /
        /
        /
        *//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        //get trigger hadron kinematics
        //PhaseSpaceHisto->GetRandom(TriggerPhaseSpaceTuple, kTRUE);

        //generate leading hadron
        double Pt;
        do{
            Pt = h1D_PtDist->GetRandom(&TR_pT);
        }while(fabs(Pt) <= DEF_HadLeadingPt);
        //get particle eta
        double Eta = TR_Eta.Uniform(-0.9, 0.9);
        //get elliptic flow strength from pt
        double pseudoparams[2] = {Pt,Eta};
        double v_2_and_v_3 = Function_v2v3ByPtAndEta(Pt, pseudoparams);
    
        /* Decision Matrix: 6 possibilities
        Eta: > 1
        Pt: < 1 >
        */
        double v_1;
        if(abs(Eta) > 0.6 && Pt < 0.2)
        {
            v_1 = h2D_v1_pT_Eta->Interpolate(0.2, 0.6*sgn(Eta));
            v_1 *= abs(Eta)/0.6;
        }
        else if(abs(Eta) > 0.6 && Pt >= 0.2 && Pt <= 9.5)
        {
            v_1 = h2D_v1_pT_Eta->Interpolate(Pt, 0.6*sgn(Eta));
            v_1 *= abs(Eta)/0.6;
        }
        else if(abs(Eta) > 0.6 && Pt > 9.5)
        {
            v_1 = h2D_v1_pT_Eta->Interpolate(9.5, 0.6*sgn(Eta));
            v_1 *= abs(Eta)/0.6;
        }
        else if(abs(Eta) <= 0.6 && Pt < 0.2)
        {
            v_1 = h2D_v1_pT_Eta->Interpolate(0.2, Eta);
        }
        else if(abs(Eta) <= 0.6 && Pt >= 0.2 && Pt <= 9.5)
        {
            v_1 = h2D_v1_pT_Eta->Interpolate(Pt, Eta);
        }
        else //eta is ok, Pt too big
        {
            v_1 = h2D_v1_pT_Eta->Interpolate(9.5, Eta);
        }
    
        h2D_Eta_vs_v1->Fill(Eta, v_1);

        double FlowParamsArray[5] = {v_1, v_2_and_v_3 * i_AllowV2, v_2_and_v_3 * i_AllowV3, V2_Angle, V3_Angle};

        double Phi = GetRandomF(Function_PhiByFlow, -Pi, +Pi, FlowParamsArray);
        if(Phi > Pi) Phi -= 2*Pi;
        else if(Phi < -Pi) Phi += 2*Pi;
        //remember particle
        ParticleMemory.push_back({Phi, Eta, Pt});

        h2D_GeneratedParticles->Fill(Eta, Phi);

        HighPtParticles.push_back({Phi, Eta, Pt});

        
        //generate subleading hadron
        do{
            Pt = h1D_PtDist->GetRandom(&TR_pT);
        }while(fabs(Pt) <= DEF_HadSubleadingPt);
        //get particle eta
        Eta = TR_Eta.Uniform(-0.9, 0.9);

        TRandom3 QCDFluctRNG(0);
        Phi += (Pi + QCDFluctRNG.Gaus(0.0, DEF_BackToBackSmearingSigma));
        while(Phi > Pi) Phi -= 2*Pi;
        while(Phi < -Pi) Phi += 2*Pi;
        //remember particle
        ParticleMemory.push_back({Phi, Eta, Pt});

        h2D_GeneratedParticles->Fill(Eta, Phi);

        HighPtParticles.push_back({Phi, Eta, Pt});


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

        for(Int_t iBGPart = 0; iBGPart < N-2; iBGPart++)
        {
            //get particle pt
            double Pt;
            do{
                Pt = h1D_PtDist->GetRandom(&TR_pT);
            }while(fabs(Pt) <= 0.15);
            if (Pt >= DEF_HadSubleadingPt) continue;//dont draw them, we dont need them
            //get particle eta
            double Eta = TR_Eta.Uniform(-0.9, 0.9);
            //get elliptic flow strength from pt
            double pseudoparams[2] = {Pt,Eta};
            double v_2_and_v_3 = Function_v2v3ByPtAndEta(Pt, pseudoparams);
        
            //get directed flow
            double v_1;
            if(abs(Eta) <= 0.6 && Pt >= 0.2)//both are in range
            {
                v_1 = h2D_v1_pT_Eta->Interpolate(Pt, Eta);
            }
            else if(Pt >= 0.2)//eta is too high
            {
                v_1 = h2D_v1_pT_Eta->Interpolate(Pt, 0.6*sgn(Eta));
                v_1 *= abs(Eta)/0.6;
            }
            else//both are out of range
            {
                v_1 = h2D_v1_pT_Eta->Interpolate(0.2, 0.6*sgn(Eta));
                v_1 *= abs(Eta)/0.6;
            }
        
            h2D_Eta_vs_v1->Fill(Eta, v_1);

            double FlowParamsArray[5] = {v_1, v_2_and_v_3 * i_AllowV2, v_2_and_v_3 * i_AllowV3, V2_Angle, V3_Angle};

            double Phi = GetRandomF(Function_PhiByFlow, -Pi, +Pi, FlowParamsArray);
            if(Phi > Pi) Phi -= 2*Pi;
            else if(Phi < -Pi) Phi += 2*Pi;
            //remember particle
            ParticleMemory.push_back({Phi, Eta, Pt});

            h2D_GeneratedParticles->Fill(Eta, Phi);

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

        if(HighPtParticles.size() == 0) continue;

        vector<PseudoJet> JetVector;

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
                double DeltaPhi = JetVector.back().phi_std() - HighPtParticles[i][0];

                //since a jet at pi and another one at -pi are equivalent, normalize
                if(DeltaPhi > Pi) DeltaPhi -= 2* Pi;
                else if(DeltaPhi < -Pi) DeltaPhi += 2*Pi;

                double DeltaR = TMath::Sqrt(DeltaEta*DeltaEta + DeltaPhi*DeltaPhi);

                if(DeltaR <= DEF_JetRadius)
                {
                    //pop candidate
                    HighPtParticles.erase(HighPtParticles.begin()+i);
                    i = -1;//start again for safety
                }
            }

        }while(HighPtParticles.size() > 0);

        //sort high pt hadrons by size to have same structure as FindJets return type
        JetVector = sorted_by_pt(JetVector);
        //make eta cut
        Selector Fiducial_cut_selector = SelectorAbsEtaMax(0.9 - DEF_JetRadius); // Fiducial cut for jets
        JetVector = Fiducial_cut_selector(JetVector);

        if(JetVector.size() > 0)
        {
            //log where jets have been found in eta
            for(int iJet = 0; iJet < JetVector.size(); iJet++)
            {
                h1D_AllJetsEta_HYBRID->Fill(JetVector[iJet].eta());
            }
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

        bool m_MirrorEvent = false;

        if(JetVector.size() < 1) continue;

        h2D_eta_vs_phi_JetCoordinates_Leading_NoCut->Fill(JetVector[0].eta(), JetVector[0].phi_std());

        if(JetVector.size() < 2 )continue;

        if(JetVector[0].eta() < 0) m_MirrorEvent = true;

        if(abs(JetVector[0].eta()) > 0.9 - DEF_JetRadius) continue;

        if(abs(JetVector[1].eta()) > 0.9 - DEF_JetRadius) continue;

        h2D_eta_vs_phi_JetCoordinates_Leading_OnlyEtaCut->Fill(JetVector[0].eta(), JetVector[0].phi_std());

        if(JetVector[0].pt() <= DEF_HadLeadingPt || JetVector[1].pt() <= DEF_HadSubleadingPt) continue;

        double PhiSeparation = JetVector[0].phi_std() - JetVector[1].phi_std();

        if(PhiSeparation > Pi) PhiSeparation -= 2*Pi;
        else if(PhiSeparation < -Pi) PhiSeparation += 2*Pi;
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
        TH2D* h2D_eta_vs_phi = new TH2D("h2D_eta_vs_phi","h2D_eta_vs_phi", DEF_BinningPerUnit * 2 * 0.9, -0.9, 0.9, DEF_BinningPerUnit * 2 * Pi, -Pi, Pi);
        #endif

        // Fill leading and subleading jet pt
        h1D_JetPt_Leading->Fill(JetVector[0].pt());
        h1D_JetPt_Subleading->Fill(JetVector[1].pt());

        // Prepare / Fill histogram on how many particles contribute to each jet
        int LeadingJetParticleCounter = 0;
        int SubleadingJetParticleCounter = 0;

        float Sgn = 1;
        if(m_MirrorEvent) Sgn = -1;

        if(LargeGap)
        {
            
            for(const auto& particle : ParticleMemory)
            {
                // Fill analysis plot
                double DeltaPhi = particle[0] - JetVector[0].phi_std();
                if(DeltaPhi > Pi) DeltaPhi -= 2*Pi;
                else if(DeltaPhi < -Pi) DeltaPhi += 2*Pi;
                if(abs(DeltaPhi) <= Pi/2) h2D_pt_vs_eta_LargeGap->Fill(particle[2], particle[1] * Sgn);

                // Fill event overview
                #if DEF_OutputEventOverviews
                    h2D_eta_vs_phi->Fill(particle[1], particle[0], particle[2]);
                #endif

                // Fill correlation relative to leading jet in phi
                if(DeltaPhi < -Pi/2) DeltaPhi += 2*Pi;
                else if(DeltaPhi > 3*Pi/2) DeltaPhi -= 2*Pi;
                if(particle[2] >= DEF_CorrelationMinPt && particle[2] <= DEF_CorrelationMaxPt)h2D_eta_vs_dphi_LargeGap->Fill(Sgn * particle[1], DeltaPhi);

                // Fill correlation relative to leading jet
                if(DeltaPhi < -Pi/2) DeltaPhi += 2*Pi;
                else if(DeltaPhi > 3*Pi/2) DeltaPhi -= 2*Pi;
                if(particle[2] >= DEF_CorrelationMinPt && particle[2] <= DEF_CorrelationMaxPt) h2D_deta_vs_dphi_LargeGap->Fill(Sgn * particle[1] - JetVector[0].eta(), DeltaPhi);

                // Fill correlation relative to event plane
                DeltaPhi = particle[0] - V2_Angle;
                if(DeltaPhi < -Pi/2) DeltaPhi += 2*Pi;
                else if(DeltaPhi > 3*Pi/2) DeltaPhi -= 2*Pi;
                if(particle[2] >= DEF_CorrelationMinPt && particle[2] <= DEF_CorrelationMaxPt)h2F_SoftHadronRelEP->Fill(DeltaPhi, particle[1]);


            }

            // Fill large gap cut jet coordinates
            h2D_eta_vs_phi_JetCoordinates_Leading_LargeGapCut->Fill(JetVector[0].eta(), JetVector[0].phi_std());

            double DeltaPhi = JetVector[0].phi_std() - V2_Angle;
            if(DeltaPhi > Pi) DeltaPhi -= 2*Pi;
            else if(DeltaPhi < -Pi) DeltaPhi += 2*Pi;

            //in case of negative eta, invert in phi because of direct flow asymmetry
            if (JetVector[0].eta() < 0)
            {
                DeltaPhi += Pi;
                if(DeltaPhi > Pi) DeltaPhi -= 2*Pi;
                else if(DeltaPhi < -Pi) DeltaPhi += 2*Pi;
            }

            h1F_PhiLeadMinusPsi_LG->Fill(DeltaPhi);

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
                if(abs(DeltaPhi) <= Pi/2) h2D_pt_vs_eta_SmallGap->Fill(particle[2], particle[1] * Sgn);

                // Fill event overview
                #if DEF_OutputEventOverviews
                    h2D_eta_vs_phi->Fill(particle[1], particle[0], particle[2]);
                #endif

                // Fill correlation relative to leading jet in phi
                if(DeltaPhi < -Pi/2) DeltaPhi += 2*Pi;
                else if(DeltaPhi > 3*Pi/2) DeltaPhi -= 2*Pi;
                if(particle[2] >= DEF_CorrelationMinPt && particle[2] <= DEF_CorrelationMaxPt)h2D_eta_vs_dphi_SmallGap->Fill(Sgn * particle[1], DeltaPhi);

                // Fill correlation relative to leading jet
                if(DeltaPhi < -Pi/2) DeltaPhi += 2*Pi;
                else if(DeltaPhi > 3*Pi/2) DeltaPhi -= 2*Pi;
                if(particle[2] >= DEF_CorrelationMinPt && particle[2] <= DEF_CorrelationMaxPt) h2D_deta_vs_dphi_SmallGap->Fill(Sgn * particle[1] - JetVector[0].eta(), DeltaPhi);

                // Fill correlation relative to event plane
                DeltaPhi = particle[0] - V2_Angle;
                if(DeltaPhi < -Pi/2) DeltaPhi += 2*Pi;
                else if(DeltaPhi > 3*Pi/2) DeltaPhi -= 2*Pi;
                if(particle[2] >= DEF_CorrelationMinPt && particle[2] <= DEF_CorrelationMaxPt)h2F_SoftHadronRelEP->Fill(DeltaPhi, particle[1]);

            }

            // Fill small gap cut jet coordinates
            h2D_eta_vs_phi_JetCoordinates_Leading_SmallGapCut->Fill(JetVector[0].eta(), JetVector[0].phi_std());

            double DeltaPhi = JetVector[0].phi_std() - V2_Angle;
            if(DeltaPhi > Pi) DeltaPhi -= 2*Pi;
            else if(DeltaPhi < -Pi) DeltaPhi += 2*Pi;

            //in case of negative eta, invert in phi because of direct flow asymmetry
            if (JetVector[0].eta() < 0)
            {
                DeltaPhi += Pi;
                if(DeltaPhi > Pi) DeltaPhi -= 2*Pi;
                else if(DeltaPhi < -Pi) DeltaPhi += 2*Pi;
            }

            h1F_PhiLeadMinusPsi_SG->Fill(DeltaPhi);

            SmallGapEventCounter++;
        }

        //write event results
        #if DEF_OutputEventOverviews
            EventHistos->cd();
            h2D_eta_vs_phi->SetTitle("");
            h2D_eta_vs_phi->GetYaxis()->SetTitleSize(DEF_AxisLabelSize);
            h2D_eta_vs_phi->GetYaxis()->SetLabelSize(DEF_AxisLabelSize);
            h2D_eta_vs_phi->GetYaxis()->SetTitle("#phi [rad]");
            h2D_eta_vs_phi->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
            h2D_eta_vs_phi->GetXaxis()->SetLabelSize(DEF_AxisLabelSize);
            h2D_eta_vs_phi->GetXaxis()->SetTitle("#eta");
            h2D_eta_vs_phi->GetZaxis()->SetTitleSize(DEF_AxisLabelSize);
            h2D_eta_vs_phi->GetZaxis()->SetLabelSize(DEF_AxisLabelSize);
            h2D_eta_vs_phi->GetZaxis()->SetTitle("#frac{d p_{T}}{d #phi d #eta}");
            h2D_eta_vs_phi->Write();
            delete h2D_eta_vs_phi;
        #endif

    }

    cout << "In the end " << ((double)(SmallGapEventCounter + LargeGapEventCounter)/(double)(i_NoOfEvents)) * 100 << " percent of events were taken into account." << endl;

    // Normalize histograms by number of events
    h2D_pt_vs_eta_LargeGap->Scale(1./(double)LargeGapEventCounter);
    h2D_pt_vs_eta_SmallGap->Scale(1./(double)SmallGapEventCounter);

    h2D_eta_vs_dphi_LargeGap->Scale(1./(double)LargeGapEventCounter);
    h2D_eta_vs_dphi_SmallGap->Scale(1./(double)SmallGapEventCounter);

    h1D_JetPt_Leading->Scale(1./(double)(LargeGapEventCounter + SmallGapEventCounter));
    h1D_JetPt_Subleading->Scale(1./(double)(LargeGapEventCounter + SmallGapEventCounter));

    h2D_GeneratedParticles->Scale(1./(double)(LargeGapEventCounter + SmallGapEventCounter));

    // Rebin Histograms if necessary
    h2D_eta_vs_dphi_LargeGap->Rebin2D(8,8);
    h2D_eta_vs_dphi_SmallGap->Rebin2D(8,8);
    h2D_eta_vs_dphi_LargeGap->Scale(1.0/8.0);
    h2D_eta_vs_dphi_SmallGap->Scale(1.0/8.0);

    h1F_PhiLeadMinusPsi_LG->Rebin(32);
    h1F_PhiLeadMinusPsi_SG->Rebin(32);

    // Divide by bin widths
    h2D_eta_vs_dphi_LargeGap->Scale(1./(h2D_eta_vs_dphi_LargeGap->GetXaxis()->GetBinWidth(0) * h2D_eta_vs_dphi_LargeGap->GetYaxis()->GetBinWidth(0)));
    h2D_eta_vs_dphi_SmallGap->Scale(1./(h2D_eta_vs_dphi_SmallGap->GetXaxis()->GetBinWidth(0) * h2D_eta_vs_dphi_SmallGap->GetYaxis()->GetBinWidth(0)));

    h2D_pt_vs_eta_LargeGap->Scale(1./(h2D_pt_vs_eta_LargeGap->GetXaxis()->GetBinWidth(0) * h2D_pt_vs_eta_LargeGap->GetYaxis()->GetBinWidth(0)));
    h2D_pt_vs_eta_SmallGap->Scale(1./(h2D_pt_vs_eta_SmallGap->GetXaxis()->GetBinWidth(0) * h2D_pt_vs_eta_SmallGap->GetYaxis()->GetBinWidth(0)));

    h2D_deta_vs_dphi_LargeGap->Scale(1./(h2D_deta_vs_dphi_LargeGap->GetXaxis()->GetBinWidth(0) * h2D_deta_vs_dphi_LargeGap->GetYaxis()->GetBinWidth(0)));
    h2D_deta_vs_dphi_SmallGap->Scale(1./(h2D_deta_vs_dphi_SmallGap->GetXaxis()->GetBinWidth(0) * h2D_deta_vs_dphi_SmallGap->GetYaxis()->GetBinWidth(0)));

    h1D_JetPt_Leading->Scale(1./h1D_JetPt_Leading->GetXaxis()->GetBinWidth(0));
    h1D_JetPt_Subleading->Scale(1./h1D_JetPt_Subleading->GetXaxis()->GetBinWidth(0));

    h2D_eta_vs_phi_JetCoordinates_Leading_NoCut->Scale(1./(h2D_eta_vs_phi_JetCoordinates_Leading_NoCut->GetXaxis()->GetBinWidth(0) * h2D_eta_vs_phi_JetCoordinates_Leading_NoCut->GetYaxis()->GetBinWidth(0)));
    h2D_eta_vs_phi_JetCoordinates_Leading_OnlyEtaCut->Scale(1./(h2D_eta_vs_phi_JetCoordinates_Leading_OnlyEtaCut->GetXaxis()->GetBinWidth(0) * h2D_eta_vs_phi_JetCoordinates_Leading_OnlyEtaCut->GetYaxis()->GetBinWidth(0)));
    h2D_eta_vs_phi_JetCoordinates_Leading_LargeGapCut->Scale(1./(h2D_eta_vs_phi_JetCoordinates_Leading_LargeGapCut->GetXaxis()->GetBinWidth(0) * h2D_eta_vs_phi_JetCoordinates_Leading_LargeGapCut->GetYaxis()->GetBinWidth(0)));
    h2D_eta_vs_phi_JetCoordinates_Leading_SmallGapCut->Scale(1./(h2D_eta_vs_phi_JetCoordinates_Leading_SmallGapCut->GetXaxis()->GetBinWidth(0) * h2D_eta_vs_phi_JetCoordinates_Leading_SmallGapCut->GetYaxis()->GetBinWidth(0)));

    h2D_GeneratedParticles->Scale(1./(h2D_GeneratedParticles->GetXaxis()->GetBinWidth(0)));

    h1D_AllJetsEta_PYTHIA->Scale(1./h1D_AllJetsEta_PYTHIA->GetXaxis()->GetBinWidth(1));
    h1D_AllJetsEta_HYBRID->Scale(1./h1D_AllJetsEta_HYBRID->GetXaxis()->GetBinWidth(1));

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

    h2D_eta_vs_dphi_LargeGap->SetTitle("Particle correlation abs. in #eta, rel. in #phi, large gaps. p_{T} #in [0,1] GeV");
    h2D_eta_vs_dphi_LargeGap->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
    h2D_eta_vs_dphi_LargeGap->GetXaxis()->SetTitle("#eta");
    h2D_eta_vs_dphi_LargeGap->GetYaxis()->SetTitleSize(DEF_AxisLabelSize);
    h2D_eta_vs_dphi_LargeGap->GetYaxis()->SetTitle("#Delta #phi[rad]");
    h2D_eta_vs_dphi_LargeGap->GetZaxis()->SetTitleSize(DEF_AxisLabelSize);
    h2D_eta_vs_dphi_LargeGap->GetZaxis()->SetTitle("#frac{1}{N_{Events}} #frac{d N}{d #phi d #eta}");

    h2D_eta_vs_dphi_SmallGap->SetTitle("Particle correlation abs. in #eta, rel. in #phi, small gaps. p_{T} #in [0,1] GeV");
    h2D_eta_vs_dphi_SmallGap->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
    h2D_eta_vs_dphi_SmallGap->GetXaxis()->SetTitle("#eta");
    h2D_eta_vs_dphi_SmallGap->GetYaxis()->SetTitleSize(DEF_AxisLabelSize);
    h2D_eta_vs_dphi_SmallGap->GetYaxis()->SetTitle("#Delta #phi[rad]");
    h2D_eta_vs_dphi_SmallGap->GetZaxis()->SetTitleSize(DEF_AxisLabelSize);
    h2D_eta_vs_dphi_SmallGap->GetZaxis()->SetTitle("#frac{1}{N_{Events}} #frac{d N}{d #phi d #eta}");

    h2D_deta_vs_dphi_LargeGap->SetTitle("Particle correlation rel. in #eta, rel. in #phi, large gaps. p_{t} #in [1,2] GeV");
    h2D_deta_vs_dphi_LargeGap->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
    h2D_deta_vs_dphi_LargeGap->GetXaxis()->SetTitle("#Delta #eta");
    h2D_deta_vs_dphi_LargeGap->GetYaxis()->SetTitleSize(DEF_AxisLabelSize);
    h2D_deta_vs_dphi_LargeGap->GetYaxis()->SetTitle("#Delta #phi[rad]");
    h2D_deta_vs_dphi_LargeGap->GetZaxis()->SetTitleSize(DEF_AxisLabelSize);
    h2D_deta_vs_dphi_LargeGap->GetZaxis()->SetTitle("#frac{1}{N_{Events}} #frac{d N}{d #phi d #eta}");

    h2D_deta_vs_dphi_SmallGap->SetTitle("Particle correlation rel. in #eta, rel. in #phi, small gaps. p_{t} #in [1,2] GeV");
    h2D_deta_vs_dphi_SmallGap->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
    h2D_deta_vs_dphi_SmallGap->GetXaxis()->SetTitle("#Delta #eta");
    h2D_deta_vs_dphi_SmallGap->GetYaxis()->SetTitleSize(DEF_AxisLabelSize);
    h2D_deta_vs_dphi_SmallGap->GetYaxis()->SetTitle("#Delta #phi[rad]");
    h2D_deta_vs_dphi_SmallGap->GetZaxis()->SetTitleSize(DEF_AxisLabelSize);
    h2D_deta_vs_dphi_SmallGap->GetZaxis()->SetTitle("#frac{1}{N_{Events}} #frac{d N}{d #phi d #eta}");

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

    h2D_eta_vs_phi_JetCoordinates_Leading_OnlyEtaCut->SetTitle("Fastjet-returned leading jet coordinates in Hybrid Events with #eta_{Leading} >= 0");
    h2D_eta_vs_phi_JetCoordinates_Leading_OnlyEtaCut->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
    h2D_eta_vs_phi_JetCoordinates_Leading_OnlyEtaCut->GetXaxis()->SetTitle("#eta");
    h2D_eta_vs_phi_JetCoordinates_Leading_OnlyEtaCut->GetYaxis()->SetTitleSize(DEF_AxisLabelSize);
    h2D_eta_vs_phi_JetCoordinates_Leading_OnlyEtaCut->GetYaxis()->SetTitle("#phi[rad]");
    h2D_eta_vs_phi_JetCoordinates_Leading_OnlyEtaCut->GetZaxis()->SetTitleSize(DEF_AxisLabelSize);
    h2D_eta_vs_phi_JetCoordinates_Leading_OnlyEtaCut->GetZaxis()->SetTitle("#frac{1}{N_{Events}} #frac{d N}{d #phi d #eta}");

    h2D_eta_vs_phi_JetCoordinates_Leading_LargeGapCut->SetTitle("Fastjet-returned leading jet coordinates in Hybrid Events. Large Gaps");
    h2D_eta_vs_phi_JetCoordinates_Leading_LargeGapCut->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
    h2D_eta_vs_phi_JetCoordinates_Leading_LargeGapCut->GetXaxis()->SetTitle("#eta");
    h2D_eta_vs_phi_JetCoordinates_Leading_LargeGapCut->GetYaxis()->SetTitleSize(DEF_AxisLabelSize);
    h2D_eta_vs_phi_JetCoordinates_Leading_LargeGapCut->GetYaxis()->SetTitle("#phi[rad]");
    h2D_eta_vs_phi_JetCoordinates_Leading_LargeGapCut->GetZaxis()->SetTitleSize(DEF_AxisLabelSize);
    h2D_eta_vs_phi_JetCoordinates_Leading_LargeGapCut->GetZaxis()->SetTitle("#frac{1}{N_{Events}} #frac{d N}{d #phi d #eta}");

    h2D_eta_vs_phi_JetCoordinates_Leading_SmallGapCut->SetTitle("Fastjet-returned leading jet coordinates in Hybrid Events. Small Gaps");
    h2D_eta_vs_phi_JetCoordinates_Leading_SmallGapCut->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
    h2D_eta_vs_phi_JetCoordinates_Leading_SmallGapCut->GetXaxis()->SetTitle("#eta");
    h2D_eta_vs_phi_JetCoordinates_Leading_SmallGapCut->GetYaxis()->SetTitleSize(DEF_AxisLabelSize);
    h2D_eta_vs_phi_JetCoordinates_Leading_SmallGapCut->GetYaxis()->SetTitle("#phi[rad]");
    h2D_eta_vs_phi_JetCoordinates_Leading_SmallGapCut->GetZaxis()->SetTitleSize(DEF_AxisLabelSize);
    h2D_eta_vs_phi_JetCoordinates_Leading_SmallGapCut->GetZaxis()->SetTitle("#frac{1}{N_{Events}} #frac{d N}{d #phi d #eta}");

    h2D_GeneratedParticles->SetTitle("Particle Correlation in Hybrid Events. #psi_{2} = 0, #psi_{3} random.");
    h2D_GeneratedParticles->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
    h2D_GeneratedParticles->GetXaxis()->SetTitle("#eta");
    h2D_GeneratedParticles->GetYaxis()->SetTitleSize(DEF_AxisLabelSize);
    h2D_GeneratedParticles->GetYaxis()->SetTitle("#phi[rad]");
    h2D_GeneratedParticles->GetZaxis()->SetTitleSize(DEF_AxisLabelSize);
    h2D_GeneratedParticles->GetZaxis()->SetTitle("#frac{1}{N_{Events}} #frac{d N}{d #phi d #eta}");

    h1D_AllJetsEta_PYTHIA->SetTitle("Pseudorapidity of all PYTHIA jets");
    h1D_AllJetsEta_PYTHIA->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
    h1D_AllJetsEta_PYTHIA->GetXaxis()->SetTitle("#eta");
    h1D_AllJetsEta_PYTHIA->GetYaxis()->SetTitleSize(DEF_AxisLabelSize);
    h1D_AllJetsEta_PYTHIA->GetYaxis()->SetTitle("counts");

    h1D_AllJetsEta_HYBRID->SetTitle("Pseudorapidity of all HYBRID jets");
    h1D_AllJetsEta_HYBRID->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
    h1D_AllJetsEta_HYBRID->GetXaxis()->SetTitle("#eta");
    h1D_AllJetsEta_HYBRID->GetYaxis()->SetTitleSize(DEF_AxisLabelSize);
    h1D_AllJetsEta_HYBRID->GetYaxis()->SetTitle("counts");

    h2D_v1_pT_Eta->SetTitle("Directed flow of hadrons");
    h2D_v1_pT_Eta->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
    h2D_v1_pT_Eta->GetXaxis()->SetTitle("p_{T}[GeV]");
    h2D_v1_pT_Eta->GetYaxis()->SetTitleSize(DEF_AxisLabelSize);
    h2D_v1_pT_Eta->GetYaxis()->SetTitle("#eta");

    h1F_PhiLeadMinusPsi_LG->SetTitle("Correlation of the leading hadron with the event plane, Large Gap");
    h1F_PhiLeadMinusPsi_LG->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
    h1F_PhiLeadMinusPsi_LG->GetXaxis()->SetTitle("#phi_{Lead} - #psi_{2}");
    h1F_PhiLeadMinusPsi_LG->GetYaxis()->SetTitleSize(DEF_AxisLabelSize);
    h1F_PhiLeadMinusPsi_LG->GetYaxis()->SetTitle("#");

    h1F_PhiLeadMinusPsi_SG->SetTitle("Correlation of the leading hadron with the event plane, Small Gap");
    h1F_PhiLeadMinusPsi_SG->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
    h1F_PhiLeadMinusPsi_SG->GetXaxis()->SetTitle("#phi_{Lead} - #psi_{2}");
    h1F_PhiLeadMinusPsi_SG->GetYaxis()->SetTitleSize(DEF_AxisLabelSize);
    h1F_PhiLeadMinusPsi_SG->GetYaxis()->SetTitle("#");

    //write global results
    Results->cd();
    h2D_pt_vs_eta_LargeGap->Write();
    h2D_pt_vs_eta_SmallGap->Write();
    h2D_eta_vs_dphi_LargeGap->Write();
    h2D_eta_vs_dphi_SmallGap->Write();
    h2D_deta_vs_dphi_LargeGap->Write();
    h2D_deta_vs_dphi_SmallGap->Write();
    h1D_JetPt_Leading->Write();
    h1D_JetPt_Subleading->Write();
    h2D_eta_vs_phi_JetCoordinates_Leading_NoCut->Write();
    h2D_eta_vs_phi_JetCoordinates_Leading_OnlyEtaCut->Write();
    h2D_eta_vs_phi_JetCoordinates_Leading_LargeGapCut->Write();
    h2D_eta_vs_phi_JetCoordinates_Leading_SmallGapCut->Write();
    h2D_GeneratedParticles->Write();

    //DEBUG
    h2D_Eta_vs_v1->Write();
    h1D_AllJetsEta_PYTHIA->Write();
    h1D_AllJetsEta_HYBRID->Write();
    h2D_v1_pT_Eta->Write();
    h2F_SoftHadronRelEP->Write();
    h1F_PhiLeadMinusPsi_LG->Write();
    h1F_PhiLeadMinusPsi_SG->Write();

    OutputFile ->Close();
    
    return 1;
}


/*PARAMETERS
DataFile is data file that was generated by FlowMC. Form: "NAME.root"
i_PtRange: e.g. particles within (1.0, 2.0)GeV/c should be analyzed -> i_PtRange = 1
i_LowPtCut: e.g. particles within (1.0, 2.0)GeV/c should be analyzed ->i_LowPtCut = 1
*/

Int_t FlowMC_Ana(const TString DataFile, float DiffRange, bool IsJetMeasurement, int LeadPt, int SublPt, int V1Mult = -1) {

    #define DEF_AxisLabelSize 0.07
    #define DEF_AxisTitleSize 0.08
    #define DEF_HistoTitleSize 0.1
    #define DEF_Rebin 20
    #define DEF_LegendFontSize 0.055
    #define DEF_Margin_Top 0.05
    #define DEF_Margin_Bottom 0.23
    #define DEF_Margin_Left 0.20
    #define DEF_Margin_Left_ForMiddlePanel 0.1
    #define DEF_Margin_Right_ForMiddlePanel 0.01

    float LowPtCuts[4] = {0.0f, 1.0, 2.0f, 4.0f};
    float HighPtCuts[4] = {1.0f, 2.0f, 4.0f, 6.0f};

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
    if(!h2D_pt_vs_eta_SmallGap)
    {
        h2D_pt_vs_eta_SmallGap = (TH2D*)file->Get("Results/h2D_pt_vs_eta_SmallGap");

        if(!h2D_pt_vs_eta_SmallGap){ cout << "Small Gap histogram not found!" << endl; return 0;}
    }
    TH2D* h2D_pt_vs_eta_LargeGap = (TH2D*)file->Get("h2D_pt_vs_eta_LargeGap");
    if(!h2D_pt_vs_eta_LargeGap)
    {
        h2D_pt_vs_eta_LargeGap = (TH2D*)file->Get("Results/h2D_pt_vs_eta_LargeGap");

        if(!h2D_pt_vs_eta_LargeGap){ cout << "Large Gap histogram not found!" << endl; return 0;}
    }
    

    //get X(pt) binning
    double BinsPerMomentum = h2D_pt_vs_eta_SmallGap->GetNbinsX() / h2D_pt_vs_eta_SmallGap->GetXaxis()->GetXmax();

    //configure and draw canvas
    TCanvas* can_ParticleMultiplicities = new TCanvas("ParticleMultiplicities","ParticleMultiplicities",1000,1000);
    can_ParticleMultiplicities->Divide(2,4, 0, 0);

    for(int iRow=0; iRow<=3; iRow++)
    {
        //ABSOLUTE DISTROS
        can_ParticleMultiplicities->cd(1 + (2*iRow));

        TH1D* h1D_PartMult_SmallGap = h2D_pt_vs_eta_SmallGap->ProjectionY("h1D_PartMult_SmallGap", LowPtCuts[iRow] * BinsPerMomentum, HighPtCuts[iRow] * BinsPerMomentum);
        h1D_PartMult_SmallGap->Scale(1./(float)((HighPtCuts[iRow] - LowPtCuts[iRow]) * BinsPerMomentum));
        
        gPad->SetMargin(DEF_Margin_Left, DEF_Margin_Right_ForMiddlePanel, DEF_Margin_Bottom*(iRow==3), DEF_Margin_Top*(iRow==0));
        //markings
        h1D_PartMult_SmallGap->SetTitle("");
        h1D_PartMult_SmallGap->SetTitleSize(DEF_HistoTitleSize);
        h1D_PartMult_SmallGap->GetXaxis()->SetTitleSize(DEF_AxisTitleSize);
        h1D_PartMult_SmallGap->GetYaxis()->SetTitleSize(DEF_AxisTitleSize);
        h1D_PartMult_SmallGap->GetXaxis()->SetLabelSize(DEF_AxisLabelSize);
        h1D_PartMult_SmallGap->GetYaxis()->SetLabelSize(DEF_AxisLabelSize);
        (iRow==0) ? (h1D_PartMult_SmallGap->GetYaxis()->SetTitle("C_{Large/Small Gap}")) : (h1D_PartMult_SmallGap->GetYaxis()->SetTitle(""));
        h1D_PartMult_SmallGap->Rebin(DEF_Rebin);
        h1D_PartMult_SmallGap->Scale(1./DEF_Rebin);
        h1D_PartMult_SmallGap->SetMarkerStyle(kCircle);
        h1D_PartMult_SmallGap->SetMarkerColor(kBlue);
        TH1* SmallCopy = (TH1*)h1D_PartMult_SmallGap->DrawCopy();

        TH1D* h1D_PartMult_LargeGap = h2D_pt_vs_eta_LargeGap->ProjectionY("h1D_PartMult_LargeGap", LowPtCuts[iRow] * BinsPerMomentum, HighPtCuts[iRow] * BinsPerMomentum);
        h1D_PartMult_LargeGap->Scale(1./(float)((HighPtCuts[iRow] - LowPtCuts[iRow]) * BinsPerMomentum));
    
        h1D_PartMult_LargeGap->SetTitle("");
        h1D_PartMult_LargeGap->SetTitleSize(DEF_HistoTitleSize);
        h1D_PartMult_LargeGap->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
        h1D_PartMult_LargeGap->GetYaxis()->SetTitleSize(DEF_AxisLabelSize);
        h1D_PartMult_LargeGap->GetXaxis()->SetLabelSize(DEF_AxisLabelSize);
        h1D_PartMult_LargeGap->GetYaxis()->SetLabelSize(DEF_AxisLabelSize);
        h1D_PartMult_LargeGap->Rebin(DEF_Rebin);
        h1D_PartMult_LargeGap->Scale(1./DEF_Rebin);
        h1D_PartMult_LargeGap->SetMarkerStyle(kCircle);
        h1D_PartMult_LargeGap->SetMarkerColor(kRed);
        TH1* LargeCopy = (TH1*)h1D_PartMult_LargeGap->DrawCopy("SAME");


        TLegend *leg = new TLegend(0.25,0.77,0.26,0.93);
        if(iRow == 0)
        {
            leg->AddEntry((TObject*)nullptr,"This thesis","");
            leg->AddEntry((TObject*)nullptr,"Pb--Pb Simulation @ #sqrt{s_{NN}} = 5.02 TeV","");
            if(IsJetMeasurement) leg->AddEntry((TObject*)nullptr,Form("Dijet, R = 0.2, |#eta_{jet}| #leq 0.7, p_{T} #geq (%.1d, %.1d) GeV/c", LeadPt, SublPt), "");
            else leg->AddEntry((TObject*)nullptr, Form("Dihadron, |#eta_{Lead./Subl. Hadron}| #leq 0.7, p_{T} #geq (%.1d, %.1d) GeV/c", LeadPt, SublPt), "");
            leg->AddEntry((TObject*)nullptr,Form("|#phi - #phi_{Leading}| #leq #pi/2"),"");
            leg->AddEntry((TObject*)nullptr,"statistical errors only", "");
            if(V1Mult != -1) leg->AddEntry((TObject*)nullptr,Form("%.1d times bigger v_{1}(#eta)", V1Mult), "");
        }
        leg->AddEntry((TObject*)nullptr,Form("%.1f #leq p_{T} #leq %.1f GeV/c", LowPtCuts[iRow], HighPtCuts[iRow]),"");
        leg->SetLineColor(10);
        leg->SetTextSize(DEF_LegendFontSize); 
        leg->Draw();

        if(iRow==0)
        {
            TLegend *leg2 = new TLegend(0.6,0.6,0.8,0.9);
            leg2->AddEntry(LargeCopy,"Large gap events","p");
            leg2->AddEntry(SmallCopy,"Small gap events","p");
            leg2->SetLineColor(10);
            leg2->SetTextSize(DEF_LegendFontSize); 
            leg2->Draw();
        }

        //draw difference
        can_ParticleMultiplicities->cd(2 + (2*iRow));

        gPad->SetMargin(DEF_Margin_Left, DEF_Margin_Right_ForMiddlePanel, DEF_Margin_Bottom*(iRow==3), DEF_Margin_Top*(iRow==0));
        
        h1D_PartMult_LargeGap->Add(h1D_PartMult_SmallGap, -1);

        h1D_PartMult_LargeGap->SetTitle("");
        (iRow==0) ? (h1D_PartMult_LargeGap->GetYaxis()->SetTitle("C_{Large Gap} - C_{Small Gap}")) : (h1D_PartMult_LargeGap->GetYaxis()->SetTitle(""));
        h1D_PartMult_LargeGap->SetTitleSize(DEF_HistoTitleSize);
        (iRow == 3) ? (h1D_PartMult_LargeGap->GetXaxis()->SetTitle("#eta")) : (h1D_PartMult_LargeGap->GetXaxis()->SetTitle(""));
        h1D_PartMult_LargeGap->GetXaxis()->SetTitleSize(DEF_AxisTitleSize);
        h1D_PartMult_LargeGap->GetYaxis()->SetTitleSize(DEF_AxisTitleSize);
        h1D_PartMult_LargeGap->GetXaxis()->SetLabelSize(DEF_AxisLabelSize);
        h1D_PartMult_LargeGap->GetYaxis()->SetLabelSize(DEF_AxisLabelSize);
        h1D_PartMult_LargeGap->GetYaxis()->SetMaxDigits(6);
        //h1D_PartMult_LargeGap->GetYaxis()->SetRangeUser(-DiffRange, DiffRange);
        h1D_PartMult_LargeGap->SetMarkerStyle(kCircle);
        h1D_PartMult_LargeGap->SetMarkerColor(kBlack);
        h1D_PartMult_LargeGap->DrawCopy();

        TLine *line = new TLine(-.9, 0.0, .9, 0.0);
        line->SetLineColor(kGray);
        line->Draw("SAME");

        h1D_PartMult_LargeGap->DrawCopy("SAME");

        //get chi^2 of null hypothesis test
        int NBins = 5;//h1F_LargeGap[iRow]->GetXaxis()->GetNbins();
        float Sum = 0;

        for(int iBin = 1; iBin <= NBins; iBin++)
        {
            float Num = pow(h1D_PartMult_LargeGap->GetBinContent(iBin), 2);
            float Den = pow(h1D_PartMult_LargeGap->GetBinError(iBin),2);
            cout << Num/Den << endl;
            Sum += Num/Den;
        }

        Sum /= NBins-1;

        TLegend *ChiLeg = new TLegend(0.7,0.2,0.8,0.35);
        ChiLeg->AddEntry((TObject*)nullptr, "Null hypothesis","");
        ChiLeg->AddEntry((TObject*)nullptr, "for #eta #leq 0","");
        ChiLeg->AddEntry((TObject*)nullptr,Form("#Chi^{2} = %.2f", Sum),"");
        ChiLeg->SetLineColor(10);
        ChiLeg->SetTextSize(DEF_LegendFontSize); 
        ChiLeg->Draw();

    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /*
    /
    /                                   ANALYZE THE RECONSTRUCTED JET SHIFT DUE TO BACKGROUND FLOW 
    /
    /
    /
    /
    /
    *//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    #define LowPtJetFrom 20.0f //GeV
    #define MiddlePtJetFrom 40.0f //GeV
    #define HighPtJetFrom 70.0f //GeV
    #define DEF_ShiftBins 60
    #define DEF_Rebin 3
    #define DEF_AxisLabelSize 0.05
    #define DEF_AxisNumbersSize 0.04
    #define DEF_LegendFontSize 0.025

    //load histograms from root file
    TH3F* h3F_ShiftV2 = (TH3F*)file->Get("h3D_LeadDist_vs_PhiShift_vs_Pt_V2");
    if(!h3F_ShiftV2){ cout << "V2 shift histogram not found!" << endl; return 0;}
    TH3F* h3F_ShiftV3 = (TH3F*)file->Get("h3D_LeadDist_vs_PhiShift_vs_Pt_V3");
    if(!h3F_ShiftV3){ cout << "V3 shift histogram not found!" << endl; return 0;}
    
    //V2 SHIFT
    float ZBinWidth = h3F_ShiftV2->GetZaxis()->GetBinWidth(1);
    
    //project
    TAxis *z = h3F_ShiftV2->GetZaxis();
    h3F_ShiftV2->GetZaxis()->SetRange(z->FindBin((float)LowPtJetFrom), z->FindBin((float)MiddlePtJetFrom) - 1);
    TH2F* h2F_LowPtJetsProjection = (TH2F*) h3F_ShiftV2->Project3D("yx");
    h2F_LowPtJetsProjection->SetName("h2F_LowPtJetsProjection");
    //h2F_LowPtJetsProjection->GetXaxis()->SetRangeUser(-Pi/2, Pi/2);
    //h2F_LowPtJetsProjection->GetYaxis()->SetRangeUser(-Pi/2, Pi/2);
    h3F_ShiftV2->GetZaxis()->SetRange(z->FindBin((float)MiddlePtJetFrom), z->FindBin((float)HighPtJetFrom) - 1);
    TH2F* h2F_MiddlePtJetsProjection = (TH2F*) h3F_ShiftV2->Project3D("yx");
    h2F_MiddlePtJetsProjection->SetName("h2F_MiddlePtJetsProjection");
    h3F_ShiftV2->GetZaxis()->SetRange(z->FindBin((float)HighPtJetFrom), z->GetNbins());
    TH2F* h2F_HighPtJetsProjection = (TH2F*) h3F_ShiftV2->Project3D("yx");
    h2F_HighPtJetsProjection->SetName("h2F_HighPtJetsProjection");

    TCanvas* can_V2Shift = new TCanvas("V2 Shift","V2 Shift",1000,1000);

    TH1D* h1D_HighPtJetsShift = new TH1D("h1D_HighPtJetsShift", "h1D_HighPtJetsShift", DEF_ShiftBins, -Pi/2, +Pi/2);
    TH1D* h1D_MiddlePtJetsShift = new TH1D("h1D_MiddlePtJetsShift", "h1D_MiddlePtJetsShift", DEF_ShiftBins, -Pi/2, +Pi/2);
    TH1D* h1D_LowPtJetsShift = new TH1D("h1D_LowPtJetsShift", "h1D_LowPtJetsShift", DEF_ShiftBins, -Pi/2, +Pi/2);
    
    float XBinWidth = h2F_LowPtJetsProjection->GetXaxis()->GetBinWidth(1);
    float BinsPerProjection = (Pi/(float)DEF_ShiftBins)/XBinWidth;


    for(int iBin = 0; iBin < DEF_ShiftBins; iBin++)
    {
        //Low pt
        TH1D* TempProjection = h2F_LowPtJetsProjection->ProjectionY("TempProjection", (iBin)*BinsPerProjection, (iBin+1)*BinsPerProjection);
        h1D_LowPtJetsShift->SetBinContent(iBin+1, TempProjection->GetMean());
        h1D_LowPtJetsShift->SetBinError(iBin+1, TempProjection->GetMeanError());
        
        //Middle pt
        TempProjection = h2F_MiddlePtJetsProjection->ProjectionY("TempProjection", iBin*BinsPerProjection, (iBin+1)*BinsPerProjection);
        h1D_MiddlePtJetsShift->SetBinContent(iBin+1, TempProjection->GetMean());
        h1D_MiddlePtJetsShift->SetBinError(iBin+1, TempProjection->GetMeanError());

        //High pt
        TempProjection = h2F_HighPtJetsProjection->ProjectionY("TempProjection", iBin*BinsPerProjection, (iBin+1)*BinsPerProjection);
        h1D_HighPtJetsShift->SetBinContent(iBin+1, TempProjection->GetMean());
        h1D_HighPtJetsShift->SetBinError(iBin+1, TempProjection->GetMeanError());

        delete(TempProjection);
    }

    TF1  *Fctn_Sinewave_Low = new TF1("sine_low","[0]*sin(2*x)",-Pi,Pi);
    TF1  *Fctn_Sinewave_Mid = new TF1("sine_mid","[0]*sin(2*x)",-Pi,Pi);
    TF1  *Fctn_Sinewave_Hig = new TF1("sine_hig","[0]*sin(2*x)",-Pi,Pi);
    float FitParams[3];//holds the sine fit param of the three jet classes [0]: Low, [1]: Middle, [2]: High

    //markings
    h1D_LowPtJetsShift->SetTitle("");
    h1D_LowPtJetsShift->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
    h1D_LowPtJetsShift->GetXaxis()->SetLabelSize(DEF_AxisNumbersSize);
    h1D_LowPtJetsShift->GetXaxis()->SetTitle("d [rad]");
    h1D_LowPtJetsShift->GetXaxis()->CenterTitle();
    h1D_LowPtJetsShift->GetYaxis()->SetTitleSize(DEF_AxisLabelSize);
    h1D_LowPtJetsShift->GetYaxis()->SetLabelSize(DEF_AxisNumbersSize);
    h1D_LowPtJetsShift->GetYaxis()->SetTitle("<#Delta #phi> [rad]");
    h1D_LowPtJetsShift->GetYaxis()->CenterTitle();
    h1D_LowPtJetsShift->SetMarkerStyle(kFullCircle);
    h1D_LowPtJetsShift->SetMarkerColor(kRed);
    h1D_LowPtJetsShift->SetLineColor(kRed);
    h1D_LowPtJetsShift->Rebin(DEF_Rebin);
    h1D_LowPtJetsShift->Scale(1./DEF_Rebin);
    TFitResultPtr FitResultPointer = h1D_LowPtJetsShift->Fit("sine_low", "S");
    //TH1D *hint = new TH1D("hint", "Fit with .95 conf.band", DEF_ShiftBins, -Pi/2, +Pi/2);
    //TVirtualFitter::GetFitter()->GetConfidenceIntervals(hint, 0.95);
    //hint->SetFillColor(kBlue-9);
    TF1 *Temp = (TF1*)h1D_LowPtJetsShift->GetListOfFunctions()->FindObject("sine_low");
    Temp->SetLineColor(kRed);
    FitParams[0] = Temp->GetParameter(0);

    //h1D_MiddlePtJetsShift->SetTitle(Form("Jets with p_{t} #in [%.3f,%.3f]", MiddlePtJetFrom, HighPtJetFrom));
    h1D_MiddlePtJetsShift->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
    h1D_MiddlePtJetsShift->GetXaxis()->SetLabelSize(DEF_AxisNumbersSize);
    h1D_MiddlePtJetsShift->GetXaxis()->SetTitle("Jet distance in #phi to closest V_{2} maximum [rad]");
    h1D_MiddlePtJetsShift->GetYaxis()->SetTitleSize(DEF_AxisLabelSize);
    h1D_MiddlePtJetsShift->GetYaxis()->SetLabelSize(DEF_AxisNumbersSize);
    h1D_MiddlePtJetsShift->GetYaxis()->SetTitle("Jet shift in #phi [rad]");
    h1D_MiddlePtJetsShift->SetMarkerStyle(kFullCircle);
    h1D_MiddlePtJetsShift->SetMarkerColor(kGreen);
    h1D_MiddlePtJetsShift->SetLineColor(kGreen);
    h1D_MiddlePtJetsShift->Rebin(DEF_Rebin);
    h1D_MiddlePtJetsShift->Scale(1./DEF_Rebin);
    h1D_MiddlePtJetsShift->Fit("sine_mid");
    TF1 *Temp2 = (TF1*)h1D_MiddlePtJetsShift->GetListOfFunctions()->FindObject("sine_mid");
    Temp2->SetLineColor(kGreen);
    FitParams[1] = Temp2->GetParameter(0);

    
    //h1D_HighPtJetsShift->SetTitle(Form("Jets with p_{t} #in [%.3f,#infty]", HighPtJetFrom));
    h1D_HighPtJetsShift->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
    h1D_HighPtJetsShift->GetXaxis()->SetLabelSize(DEF_AxisNumbersSize);
    h1D_HighPtJetsShift->GetXaxis()->SetTitle("Jet distance to closest V_{2} maximum in #phi [rad]");
    h1D_HighPtJetsShift->GetYaxis()->SetTitleSize(DEF_AxisLabelSize);
    h1D_HighPtJetsShift->GetYaxis()->SetLabelSize(DEF_AxisNumbersSize);
    h1D_HighPtJetsShift->GetYaxis()->SetTitle("Jet shift in #phi [rad]");
    h1D_HighPtJetsShift->SetMarkerStyle(kFullCircle);
    h1D_HighPtJetsShift->SetMarkerColor(kBlue);
    h1D_HighPtJetsShift->SetLineColor(kBlue);
    h1D_HighPtJetsShift->Rebin(DEF_Rebin);
    h1D_HighPtJetsShift->Scale(1./DEF_Rebin);
    h1D_HighPtJetsShift->Fit("sine_hig");
    TF1 *Temp3 = (TF1*)h1D_HighPtJetsShift->GetListOfFunctions()->FindObject("sine_hig");
    Temp3->SetLineColor(kBlue);
    FitParams[2] = Temp3->GetParameter(0);

    TH1* Low_Lgd = (TH1*)h1D_LowPtJetsShift->DrawCopy("P E1");
    //hint->DrawCopy("e3 same");
    TH1* Mid_Lgd = (TH1*)h1D_MiddlePtJetsShift->DrawCopy("SAME P E1");
    TH1* High_Lgd = (TH1*)h1D_HighPtJetsShift->DrawCopy("SAME P E1");

    TLegend *leg = new TLegend(0.7,0.7,0.85,0.85);
    leg->AddEntry((TObject*)nullptr,Form("Jet p_{T} of reconstructed jet before adding background"),"");
    leg->AddEntry(Low_Lgd,Form("Jets with %.1f #leq p_{T,Leading} #leq %.1f GeV/c", LowPtJetFrom, MiddlePtJetFrom),"p");
    leg->AddEntry(Mid_Lgd,Form("Jets with %.1f #leq p_{T,Leading} #leq %.1f GeV/c", MiddlePtJetFrom, HighPtJetFrom),"p");
    leg->AddEntry(High_Lgd,Form("Jets with %.1f #leq p_{T,Leading} < #infty GeV/c", HighPtJetFrom),"p");
    leg->SetLineColor(10);
    leg->SetTextSize(DEF_LegendFontSize); 
    leg->Draw();

    TLegend *info = new TLegend(0.7,0.7,0.85,0.85);
    info->AddEntry((TObject*)nullptr,Form("This thesis"),"");
    info->AddEntry((TObject*)nullptr,Form("200M events from PYTHIA p-p + artificial background (BG)"),"");
    info->AddEntry((TObject*)nullptr,Form("BG: MC generated, based on ALICE 2018 Pb-Pb, 0-10%, #sqrt{s_{NN}} = 5.02TeV/2.76TeV"),"");
    info->AddEntry((TObject*)nullptr,Form("|#eta_{Leading}| #leq 0.7 = 0.9 - R"),"");
    info->AddEntry((TObject*)nullptr,Form("Jet clustering: R=0.2, Anti-k_{t}"),"");
    info->AddEntry((TObject*)nullptr,Form("Only jets that were not shifted by more than #pi/2"),"");
    info->AddEntry((TObject*)nullptr,Form("Fitted function: A*sin(2d) with parameter A. Fit results:"),"");
    info->AddEntry((TObject*)nullptr,Form("A(Low p_{T,Leading}, Middle p_{T,L}, High p_{T,L}) = (%.3e, %.3e, %.3e)", FitParams[0], FitParams[1], FitParams[2]),"");
    info->SetLineColor(10);
    info->SetTextSize(DEF_LegendFontSize); 
    info->Draw();

    
    cout << "DONE!" << endl;
    file->Close();

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
        SubtractedJets[iJet].set_user_info(new MyUserInfo(SubtractedJets[iJet].constituents().size()));

    }

    return sorted_by_pt(SubtractedJets);

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
    //srand(static_cast<unsigned int>(std::time(0)));
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
    
    -------------------------------------------------------------------------------------------------------------------------------------------------------------------
    */

    #define ACCURACY 100

    double TotalProb = 0;
    double IntervalLength = (i_UpperLim - i_LowerLim)/ACCURACY;
    for(int8_t i=0; i<=ACCURACY-1; i++)
    {
        TotalProb += i_FuncPtr(i_LowerLim + (0.5+i)*IntervalLength, i_FuncParams) * (i_UpperLim - i_LowerLim)/ACCURACY;
    }

    //srand(static_cast<unsigned int>(std::time(0)));
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
    double v_1 = params[0];
    double v_2 = params[1];
    double v_3 = params[2];
    double v2_Angle = params[3];
    double v3_Angle = params[4];
    return (1/(2*pi))*(1+2*(v_1 * cos(x-v2_Angle) + v_2*cos(2*(x-v2_Angle)) + v_3*cos(3*(x-v3_Angle))));
}

double Function_DirectedFlowByEta(double Eta)
{
    return (pow(10, -3) * -0.7043 * Eta);
}

double Function_v2v3ByPtAndEta(double x, double params[])
{

    double Pt = params[0];
    double Eta = params[1];
    double FlowParams;
    //an approximation for pions(+) @Nadine Gruenwald
    FlowParams = 0.063 * pow(x, 1.6) * exp(-0.47*x);

    //now reduce based on eta, a fit to https://www.hepdata.net/record/ins1456145 from -1.25 <= x <= 1.25
    FlowParams *= 1 - 0.15*abs(Eta);

    return FlowParams;
}

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}