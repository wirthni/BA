#include "./StJetTrackEvent.h"
#include "./StJetTrackEventLinkDef.h"
#include "./functions.h"
#include "TClassRef.h"
#include "Math/Vector4D.h"
#include <algorithm>
#include "./PythiaEvent.h"
#include "./PythiaEvent_LinkDef.h"
//#include "Fastjet/fastjet-3.5.1/include/fastjet/ClusterSequence.hh"
#include "Fastjet/fastjet-3.5.1/include/fastjet/ClusterSequenceArea.hh"
#include "Fastjet/fastjet-3.5.1/include/fastjet/config.h"
#include "Fastjet/fastjet-3.5.1/include/fastjet/PseudoJet.hh"
#include "Fastjet/fastjet-3.5.1/include/fastjet/Selector.hh"
#include "Fastjet/fastjet-3.5.1/include/fastjet/tools/Subtractor.hh"
//#include "StPhysicalHelixD.hh"
#include "Fastjet/fastjet-3.5.1/include/fastjet/tools/JetMedianBackgroundEstimator.hh"
using namespace fastjet;
using namespace std;

const char PYTHIA_TREE[]   = "PythiaEvent";
const char PYTHIA_EVENT_BRANCH[] = "Events";

vector<PseudoJet> FindJets(vector<PseudoJet> vec_particles);
double GetPhiByProb(void);
double GetEtaByProb(void);
double GetPtByProb(void);
double GetNByProb(void);

/*Parameters:
should be clear

Additional Info:

Event Plane is always at \phi = 0
*/
Int_t FlowMC(Int_t i_NoOfEvents = -1, double i_PtLowerLimit = 40)
{

    gStyle->SetOptStat(0);
    SetRootGraphicStyle();

    #define DEF_BinningPerUnit 100
    #define DEF_UsePYTHIAJetData true

    //other variables
    TH1D* h1D_JetMultiplicityInPsi = new TH1D("h1D_JetMultiplicityInPsi", "h1D_JetMultiplicityInPsi", DEF_BinningPerUnit*2*Pi, -Pi, Pi);
    TF1 *PhiDistLow=new TF1("PhiDist","(1/(2*pi))*(1+0.05*cos(2*x))",-Pi,Pi);
    TF1 *PhiDistMiddle=new TF1("PhiDist","(1/(2*pi))*(1+0.1*cos(2*x))",-Pi,Pi);
    TF1 *PhiDistHigh=new TF1("PhiDist","(1/(2*pi))*(1+0.14*cos(2*x))",-Pi,Pi);

    TH2D* h2D_pt_vs_eta_LargeGap = new TH2D("h2D_pt_vs_eta_LargeGap","h2D_pt_vs_eta_LargeGap",100, 0, 10, DEF_BinningPerUnit * 2 * 0.9,-0.9,0.9);
    TH2D* h2D_pt_vs_eta_SmallGap = new TH2D("h2D_pt_vs_eta_SmallGap","h2D_pt_vs_eta_SmallGap",100, 0, 10, DEF_BinningPerUnit * 2 * 0.9,-0.9,0.9);

    TH2D* h2D_eta_vs_deltaphi_LargeGap = new TH2D("h2D_eta_vs_deltaphi_LargeGap","h2D_eta_vs_deltaphi_LargeGap",DEF_BinningPerUnit * 2 * 0.9, -0.9, 0.9, DEF_BinningPerUnit * 2 * Pi, -Pi/2, 3*Pi/2);
    TH2D* h2D_eta_vs_deltaphi_SmallGap = new TH2D("h2D_eta_vs_deltaphi_SmallGap","h2D_eta_vs_deltaphi_SmallGap",DEF_BinningPerUnit * 2 * 0.9, -0.9, 0.9, DEF_BinningPerUnit * 2 * Pi, -Pi/2, 3*Pi/2);

    TRandom TR_Eta;
    TR_Eta.SetSeed(0);

    TChain              *inputPYTHIA;
    PythiaEvent         *PYTHIAEvent;
    PythiaParticle      *PYTHIAParticle;
    Long64_t            file_entries;

    Long64_t            LargeGapEventCounter = 0;
    Long64_t            SmallGapEventCounter = 0;

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
        //------------------------------------
        // Read input data pythia tree
        TString pinputdirPYTHIA = "/media/nicolas-wirth/INTENSO/PYTHIA_pp/";
        //TString SEListPYTHIA = "/filelist_PYTHIA_pp.txt";
        TString SEListPYTHIA = "/filelist_PYTHIA_pp_OnlyHigh.txt";
        TString Pythia_List =  "/media/nicolas-wirth/INTENSO/PYTHIA_pp";
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
                        Long64_t file_entries;////////////////////////////////////////////////////////////////////////////////////////////////////////
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
        }
        //------------------------------------
        printf("StJetAnalysis::InitPythiaTree() finished \n");

    }

    TString str_out = "./FlowMC_Output.root";
    TFile* OutputFile = new TFile(str_out.Data(),"RECREATE");
    OutputFile->cd();
    cout << "Output file is " << str_out << endl;
    TDirectory *EventHistos = OutputFile->mkdir("Event Histos");
    TDirectory *Results = OutputFile->mkdir("Results");
    
    if(i_NoOfEvents != -1) file_entries = i_NoOfEvents;

    //loop over events
    for(Int_t iEvent = 0; iEvent < file_entries; iEvent++)
    {

        vector<PseudoJet> ParticleVector;//holds all the particles in the event
               
                
        Int_t N = h1D_NDist->GetRandom() / 5;//approx max 1000 particles
        vector<array<double, 4>> ParticleMemory;//remembers randomly generated particle coordinates (phi, eta, pt) for relevant cuts later in the analysis

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
        
        for(Int_t iBGPart = 0; iBGPart < N; iBGPart++)
        {

            double Pt = h1D_PtDist->GetRandom();
            //get elliptic flow from pt
            double Phi;
            if(Pt < 1) Phi = PhiDistLow->GetRandom();
            else if(Pt < 2) Phi = PhiDistMiddle->GetRandom();
            else Phi = PhiDistHigh->GetRandom();
            double Eta = TR_Eta.Uniform(-0.9, 0.9);

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
            ParticleVector.push_back(PseudoJet(p_x, p_y, p_z, E));

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
            if (!inputPYTHIA->GetEntry( iEvent )) // take the event -> information is stored in event
                cout << "Error" << endl;

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

                ParticleVector.push_back(PseudoJet(px_track, py_track, pz_track, E_track));

            }

        }
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /*
        /
        /                                                     
        /                                           FIND JETS WITH THE ANTI-KT ALGORITHM   
        /                                                     
        /                                                     
        /                                                     
        /
        /                                                    
        /
        /
        /
        *//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        //R = 0.2
        vector<PseudoJet> JetVector = FindJets(ParticleVector);

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /*
        /
        /                                                     
        /                                           FILTER EVENTS BASED ON THE CONSTRAINTS AND DIFFER BETW. SMALL/LARGE GAP
        /                       -DIJET EVENT
        /                       -LEADING JET PT > 40GeV, SUBLEADING JET PT > 20GeV
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

        if(true)
        {
            if(JetVector.size() != 2) continue;

            if(JetVector[0].pt() <= 40 || JetVector[1].pt() <= 20) continue;

            if(JetVector[0].eta() <= 0) continue;

            if(abs(JetVector[0].phi_std() - JetVector[1].phi_std()) <= Pi/2) continue;

        }

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

        TH2D* h2D_phi_vs_eta = new TH2D("h2D_phi_vs_eta","h2D_phi_vs_eta", DEF_BinningPerUnit * 2 * Pi, -Pi, Pi, DEF_BinningPerUnit * 2 * 0.9, -0.9, 0.9);

        if(LargeGap)
        {
            for(const auto& particle : ParticleMemory)
            {
                //only take particles that are around the leading jet
                double DeltaPhi = particle[0] - JetVector[0].phi_std();
                if(DeltaPhi > Pi) DeltaPhi -= 2*Pi;
                else if(DeltaPhi < -Pi) DeltaPhi += 2*Pi;
                if(abs(DeltaPhi) <= Pi/2) h2D_pt_vs_eta_LargeGap->Fill(particle[2], particle[1]);

                //write particle into histogram
                h2D_phi_vs_eta->Fill(particle[0], particle[1], particle[2]);

                //normalize delta phi
                if(DeltaPhi < -Pi/2) DeltaPhi += 2*Pi;
                else if(DeltaPhi > 3*Pi/2) DeltaPhi -= 2*Pi;
                h2D_eta_vs_deltaphi_LargeGap->Fill(particle[1], DeltaPhi);
            }
            LargeGapEventCounter++;
        }
        else
        {
            for(const auto& particle : ParticleMemory)
            {
                //only take particles that are around the leading jet
                double DeltaPhi = particle[0] - JetVector[0].phi_std();
                if(DeltaPhi > Pi) DeltaPhi -= 2*Pi;
                else if(DeltaPhi < -Pi) DeltaPhi += 2*Pi;
                if(abs(DeltaPhi) <= Pi/2) h2D_pt_vs_eta_SmallGap->Fill(particle[2], particle[1]);

                //write particle into histogram
                h2D_phi_vs_eta->Fill(particle[0], particle[1], particle[2]);

                //normalize delta phi
                if(DeltaPhi < -Pi/2) DeltaPhi += 2*Pi;
                else if(DeltaPhi > 3*Pi/2) DeltaPhi -= 2*Pi;
                h2D_eta_vs_deltaphi_SmallGap->Fill(particle[1], DeltaPhi);
            }
            SmallGapEventCounter++;
        }

        //write event results
        EventHistos->cd();
        h2D_phi_vs_eta->Write();
        delete h2D_phi_vs_eta;
    }

    cout << "In the end " << ((double)(SmallGapEventCounter + LargeGapEventCounter)/(double)file_entries) * 100 << "percent of events were taken into account." << endl;

    //normalize histograms
    h2D_pt_vs_eta_LargeGap->Scale(1./(double)LargeGapEventCounter);
    h2D_pt_vs_eta_SmallGap->Scale(1./(double)SmallGapEventCounter);

    //rebin histograms
    h2D_eta_vs_deltaphi_LargeGap->Rebin2D(8,8);
    h2D_eta_vs_deltaphi_SmallGap->Rebin2D(8,8);

    //write global results
    Results->cd();
    h2D_pt_vs_eta_LargeGap->Write();
    h2D_pt_vs_eta_SmallGap->Write();
    h2D_eta_vs_deltaphi_LargeGap->Write();
    h2D_eta_vs_deltaphi_SmallGap->Write();
    OutputFile ->Close();
    
    return 1;
}

Int_t FlowMC_Ana(const TString DataFile, Int_t i_PtRange, Int_t i_LowPtCut) {

    #define DEF_AxisLabelSize 0.05
    #define DEF_HistoTitleSize 0.1

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
    TH2D* h2D_pt_vs_eta_SmallGap = (TH2D*)file->Get("Results/h2D_pt_vs_eta_SmallGap");
    if(!h2D_pt_vs_eta_SmallGap){ cout << "Small Gap histogram not found!" << endl; return 0;}
    TH2D* h2D_pt_vs_eta_LargeGap = (TH2D*)file->Get("Results/h2D_pt_vs_eta_LargeGap");
    if(!h2D_pt_vs_eta_LargeGap){ cout << "Large Gap histogram not found!" << endl; return 0;}
    

    //get X(pt) binning
    double BinsPerMomentum = h2D_pt_vs_eta_SmallGap->GetNbinsX() / h2D_pt_vs_eta_SmallGap->GetXaxis()->GetXmax();

    //configure and draw canvas
    TCanvas* can_ParticleMultiplicities = new TCanvas("ParticleMultiplicities","ParticleMultiplicities",1000,1000);
    can_ParticleMultiplicities->Divide(3,1, 0, 0);

    //SMALL GAP
    can_ParticleMultiplicities->cd(1);

    TH1D* h1D_PartMult_SmallGap = h2D_pt_vs_eta_SmallGap->ProjectionY("h1D_PartMult_SmallGap", i_LowPtCut * BinsPerMomentum, (i_LowPtCut + i_PtRange) * BinsPerMomentum);

    //markings
    h1D_PartMult_SmallGap->SetTitle("Small Gap");
    h1D_PartMult_SmallGap->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
    h1D_PartMult_SmallGap->GetXaxis()->SetTitle("#eta");
    h1D_PartMult_SmallGap->GetYaxis()->SetTitleSize(DEF_AxisLabelSize);
    h1D_PartMult_SmallGap->GetYaxis()->SetTitle("#frac{1}{N_{Small Gap Events}} #frac{d N}{d #eta}");
    h1D_PartMult_SmallGap->Rebin(8);
    h1D_PartMult_SmallGap->DrawCopy();

    can_ParticleMultiplicities->cd(2);

    TH1D* h1D_PartMult_LargeGap = h2D_pt_vs_eta_LargeGap->ProjectionY("h1D_PartMult_LargeGap", i_LowPtCut * BinsPerMomentum, (i_LowPtCut + i_PtRange) * BinsPerMomentum);
    h1D_PartMult_LargeGap->SetTitle("Large Gap");
    h1D_PartMult_LargeGap->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
    h1D_PartMult_LargeGap->GetXaxis()->SetTitle("#eta");
    h1D_PartMult_LargeGap->GetYaxis()->SetTitleSize(DEF_AxisLabelSize);
    h1D_PartMult_LargeGap->GetYaxis()->SetTitle("#frac{1}{N_{Large Gap Events}} #frac{d N}{d #eta}");
    h1D_PartMult_LargeGap->Rebin(8);
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
    file->Close();

    return 1;
}

vector<PseudoJet> FindJets(vector<PseudoJet> vec_particles)
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

    Double_t ghost_maxrap = 0.9; // Fiducial cut for background estimation
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
    Selector Pt_selector = SelectorPtMin(20.0);
    Selector N_Hardest = SelectorNHardest(2);
    Selector Selector_All = Pt_selector && Fiducial_cut_selector;// && N_Hardest;

    jets_fiducial = Selector_All(jets_all);

    // background estimation
    // Standard rho calculation with kT algorithm and N hardest jets removed

    
    Int_t Rem_n_hardest = 2; 
    Int_t eRem_n_hardest = 2;////???????
    double eBkg_R = 0.4;
    JetDefinition jet_def_bkgd(kt_algorithm, eBkg_R); // <--
    AreaDefinition area_def_bkgd(active_area_explicit_ghosts,GhostedAreaSpec(ghost_maxrap,1,0.005));
    Int_t Rem_n_hardest_use = eRem_n_hardest;
    Selector Selector_Bkg = SelectorAbsEtaMax(0.9 - eBkg_R) * (!SelectorNHardest(Rem_n_hardest)); // <--

    JetMedianBackgroundEstimator bkgd_estimator(Selector_Bkg, jet_def_bkgd, area_def_bkgd); // <--
    Subtractor subtractor(&bkgd_estimator);
    bkgd_estimator.set_particles(vec_particles);
    Double_t jet_rho   = bkgd_estimator.rho();
    

    // loop over jets  
    Double_t Jet_area_cut = 0.56*TMath::Pi()*TMath::Power(jet_radius,2);
    for(Int_t i_jet = 0; i_jet < (Int_t)jets_fiducial.size(); i_jet++)
    {
        Float_t jet_pt         = jets_fiducial[i_jet].perp();
        if(jet_pt <= 0.0) continue;
        Float_t jet_area       = jets_fiducial[i_jet].area();
        //Float_t jet_pt_sub     = jets_fiducial[i_jet].perp() - jet_rho*jet_area;
        Float_t jet_eta        = jets_fiducial[i_jet].eta();
        Float_t jet_phi        = jets_fiducial[i_jet].phi();
        Float_t jet_phi_std    = jets_fiducial[i_jet].phi_std();
        if(jet_area < Jet_area_cut) continue; 


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

    return jets_fiducial;

}
