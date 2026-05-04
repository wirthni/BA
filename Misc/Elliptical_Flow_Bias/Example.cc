#include "./StJetTrackEvent.h"
#include "./StJetTrackEventLinkDef.h"
#include "./functions.h"
#include <algorithm>

#if false
#include "fastjet/config.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/Subtractor.hh"
#include "StPhysicalHelixD.hh"
using namepsace fastjet;
using namespace std;
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#endif

const char JETTRACK_EVENT_TREE[]   = "JetTrackEvent";
const char JETTRACK_EVENT_BRANCH[] = "Events";

/*Parameters:
i_Acceptance: e.g. 2.5 like CMS
i_Range: Eta range (plusminus) in which jets are generated (e.g. 5)    
i_Sigma: Sigma of Gauss, e.g. 0.25
i_MXXParticles_XX: Number of particles in one jet
i_MinAsymmDistance: Lower eta limit for asymmetric jet, e.g. 1-1.5 in CMS paper
i_Scale:scale result, probably unnecessary
*/
void pp_vs_pbpb(double i_Acceptance, double i_Range, Int_t i_Events, double i_Sigma, Int_t i_MinParticles_pp, Int_t i_MaxParticles_pp, Int_t i_MinParticles_pbpb, Int_t i_MaxParticles_pbpb, double i_MinAsymmDistance, Int_t i_Scale)
{

    gStyle->SetOptStat(0);
    SetRootGraphicStyle();

    #define ParticleArraySize 150
    #define Bins 300
    #define MinSecondaryParticlePt 1
    #define MaxSecondaryParticlePt 7
    #define MinPrimaryParticlePt_pp 30
    #define MaxPrimaryParticlePt_pp 80
    #define MinPrimaryParticlePt_pbpb 15
    #define MaxPrimaryParticlePt_pbpb 50
    #define ENUM_pp 0
    #define ENUM_pbpb 1

    //histograms design
    #define DEF_AxisLabelSize 0.04
    #define DEF_HistoTitleSize 0.1
    #define DEF_InfoBoxX 0.95
    #define DEF_InfoBoxY 0.9
    #define DEF_LastColumnInfoXOffset -0.15
    #define DEF_InfoBoxLineDist 0.03
    #define DEF_LatexTextSize 0.03

    TRandom TR_JetPos1;
    TRandom TR_JetPos2;
    TR_JetPos1.SetSeed(0);
    TR_JetPos2.SetSeed(0);

    TH1D* h1D_SymmDist_pp;
    TH1D* h1D_SymmDist_pbpb;
    TH1D* h1D_AsymmDist_pp;
    TH1D* h1D_AsymmDist_pbpb;
    TH1D* h1D_Difference_pp;
    TH1D* h1D_Difference_pbpb;

    //get DeltaEta Distribution
    TH1D* h1D_DeltaEtaDist = new TH1D("DeltaEtaDist", "DeltaEtaDist", Bins, -i_Range, i_Range);
    for(Int_t i=0; i<i_Events; i++)
    {
        double Jet1Pos = TR_JetPos1.Uniform(-i_Range, i_Range);
        double Jet2Pos = TR_JetPos2.Uniform(-i_Range, i_Range);

        if(abs(Jet1Pos) > i_Acceptance || abs(Jet2Pos) > i_Acceptance) continue;

        h1D_DeltaEtaDist->Fill(Jet2Pos - (-Jet1Pos));
    }

    //plot h1D_DeltaEtaDist
    TCanvas* can_SeparationDist = new TCanvas("can_SeparationDist","can_SeparationDist",1000,1000);
    gPad->SetMargin(0.16, 0.15, 0.15, 0.02);
    h1D_DeltaEtaDist->Scale(1/h1D_DeltaEtaDist->Integral());
    h1D_DeltaEtaDist->Rebin(4);
    h1D_DeltaEtaDist->SetMarkerStyle(kFullCircle);
    h1D_DeltaEtaDist->SetMarkerColor(kBlack);
    h1D_DeltaEtaDist->SetTitle("");
    h1D_DeltaEtaDist->GetXaxis()->SetTitle("#Delta #eta^{jet_{1}, jet_{2}}");
    h1D_DeltaEtaDist->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
    h1D_DeltaEtaDist->GetXaxis()->SetLabelSize(DEF_AxisLabelSize);
    h1D_DeltaEtaDist->GetYaxis()->SetTitle("#frac{dp(#Delta #eta^{jet_{1}, jet_{2}})}{d #Delta #eta^{jet_{1}, jet_{2}}}");
    h1D_DeltaEtaDist->GetYaxis()->SetTitleSize(DEF_AxisLabelSize);
    h1D_DeltaEtaDist->GetYaxis()->SetLabelSize(DEF_AxisLabelSize);
    h1D_DeltaEtaDist->Draw();

    //do p-p in the first cycle and pb-pb in the second
    for(Int_t k=0; k<=1; k++)
    {
        if(k==ENUM_pp)
        {
            h1D_SymmDist_pp = new TH1D("SymmDist_pp", "SymmDist_pp", Bins, -i_Acceptance, i_Acceptance);
            h1D_AsymmDist_pp = new TH1D("AsymmDist_pp", "AsymmDist_pp", Bins, -i_Acceptance, i_Acceptance);
        }
        else
        {
            h1D_SymmDist_pbpb = new TH1D("SymmDist_pbpb", "SymmDist_pbpb", Bins, -i_Acceptance, i_Acceptance);
            h1D_AsymmDist_pbpb = new TH1D("AsymmDist_pbpb", "AsymmDist_pbpb", Bins, -i_Acceptance, i_Acceptance);
        }
    
        Int_t SymmetricEventCounter = 0;
        Int_t AsymmetricEventCounter = 0;

        TRandom TR_JetSeparation;
        TRandom TR_ParticlesInJet;
        TRandom TR_X_Particle1;
        TRandom TR_Pt_Particle1;
        TRandom TR_Pt_Particle2;
        TRandom TR_X_Particle2;
        TRandom TR_Pt_Leading1;
        TRandom TR_Pt_Leading2;
        TR_JetSeparation.SetSeed(0);
        TR_ParticlesInJet.SetSeed(0);
        TR_X_Particle1.SetSeed(0);
        TR_X_Particle2.SetSeed(0);
        TR_Pt_Particle1.SetSeed(0);
        TR_Pt_Particle2.SetSeed(0);
        TR_Pt_Leading1.SetSeed(0);
        TR_Pt_Leading2.SetSeed(0);


        for(Int_t i=0; i<i_Events; i++)
        {
            //SEPARATION IS X_SUBLEADING - X_LEADING

            //get separation of the two jets by probability density
            double Separation = h1D_DeltaEtaDist->GetRandom();
            if(Separation > 0) Separation *= -1;

            Int_t ParticlesInJet1, ParticlesInJet2;

            if(k==ENUM_pp)
            {
                ParticlesInJet1 = TR_ParticlesInJet.Integer(i_MaxParticles_pp - i_MinParticles_pp) + i_MinParticles_pp;
                ParticlesInJet2 = TR_ParticlesInJet.Integer(i_MaxParticles_pp - i_MinParticles_pp) + i_MinParticles_pp;
            }
            else
            {
                ParticlesInJet1 = TR_ParticlesInJet.Integer(i_MaxParticles_pbpb - i_MinParticles_pbpb) + i_MinParticles_pbpb;
                ParticlesInJet2 = TR_ParticlesInJet.Integer(i_MaxParticles_pbpb - i_MinParticles_pbpb) + i_MinParticles_pbpb;
            }

            double cg1, cg2;
            double ParticleInfo1[ParticleArraySize][2];//0: x position, 1: pt
            double ParticleInfo2[ParticleArraySize][2];//0: x position, 1: pt
            //produce two realistic discrete almost-gaussian distributions and the leading particle
            if(k == ENUM_pp)//pp
            {
                //JET1
                for(Int_t i=0; i<ParticlesInJet1; i++)
                {
                    ParticleInfo1[i][0] = TR_X_Particle1.Gaus(0, i_Sigma);
                    ParticleInfo1[i][1] = TR_Pt_Particle1.Uniform(MinSecondaryParticlePt, MaxSecondaryParticlePt);
                }
                ParticleInfo1[ParticlesInJet1][0] = 0;//leading particle
                ParticleInfo1[ParticlesInJet1][1] = TR_Pt_Leading1.Uniform(MinPrimaryParticlePt_pp, MaxPrimaryParticlePt_pp);

                cg1 = 0;
                double SumOfPts1 = 0;
                for(Int_t i=0; i<ParticlesInJet1+1; i++)//loop also over the leading particle
                {
                    cg1 += ParticleInfo1[i][0] * ParticleInfo1[i][1];//weight according to pt of the particle
                    SumOfPts1 += ParticleInfo1[i][1];
                }
                cg1 /= SumOfPts1;

                //JET2
                for(Int_t i=0; i<ParticlesInJet2; i++)
                {
                    ParticleInfo2[i][0] = TR_X_Particle2.Gaus(Separation, i_Sigma);
                    ParticleInfo2[i][1] = TR_Pt_Particle2.Uniform(MinSecondaryParticlePt, MaxSecondaryParticlePt);
                }
                ParticleInfo2[ParticlesInJet2][0] = Separation;//leading particle
                ParticleInfo2[ParticlesInJet2][1] = TR_Pt_Leading2.Uniform(MinPrimaryParticlePt_pp, MaxPrimaryParticlePt_pp);

                cg2 = 0;
                double SumOfPts2 = 0;
                for(Int_t i=0; i<ParticlesInJet2+1; i++)//loop also over the leading particle
                {
                    cg2 += ParticleInfo2[i][0] * ParticleInfo2[i][1];//weight according to pt of the particle
                    SumOfPts2 += ParticleInfo2[i][1];
                }
                cg2 /= SumOfPts2;
            }
            else//pbpb
            {
                //JET1
                for(Int_t i=0; i<ParticlesInJet1; i++)
                {
                    ParticleInfo1[i][0] = TR_X_Particle1.Gaus(0, i_Sigma);
                    ParticleInfo1[i][1] = TR_Pt_Particle1.Uniform(MinSecondaryParticlePt, MaxSecondaryParticlePt);
                }
                ParticleInfo1[ParticlesInJet1][0] = 0;//leading particle
                ParticleInfo1[ParticlesInJet1][1] = TR_Pt_Leading1.Uniform(MinPrimaryParticlePt_pbpb, MaxPrimaryParticlePt_pbpb);

                cg1 = 0;
                double SumOfPts1 = 0;
                for(Int_t i=0; i<ParticlesInJet1+1; i++)//loop also over the leading particle
                {
                    cg1 += ParticleInfo1[i][0] * ParticleInfo1[i][1];//weight according to pt of the particle
                    SumOfPts1 += ParticleInfo1[i][1];
                }
                cg1 /= SumOfPts1;

                //JET2
                for(Int_t i=0; i<ParticlesInJet2; i++)
                {
                    ParticleInfo2[i][0] = TR_X_Particle2.Gaus(Separation, i_Sigma);
                    ParticleInfo2[i][1] = TR_Pt_Particle2.Uniform(MinSecondaryParticlePt, MaxSecondaryParticlePt);
                }
                ParticleInfo2[ParticlesInJet2][0] = Separation;//leading particle
                ParticleInfo2[ParticlesInJet2][1] = TR_Pt_Leading2.Uniform(MinPrimaryParticlePt_pbpb, MaxPrimaryParticlePt_pbpb);

                cg2 = 0;
                double SumOfPts2 = 0;
                for(Int_t i=0; i<ParticlesInJet2+1; i++)//loop also over the leading particle
                {
                    cg2 += ParticleInfo2[i][0] * ParticleInfo2[i][1];//weight according to pt of the particle
                    SumOfPts2 += ParticleInfo2[i][1];
                }
                cg2 /= SumOfPts2;
            }

            //categorize into the two cases and write
            double RealSeparation = cg2 - cg1;
            //double RealSeparation = Separation;
            if(RealSeparation > -0.5)
            {
                SymmetricEventCounter += 1;
                //fill symmetrichistogram
                for(Int_t i=0; i<ParticlesInJet1+1; i++)
                {
                    if(k==ENUM_pp) h1D_SymmDist_pp->Fill(ParticleInfo1[i][0]);
                    else h1D_SymmDist_pbpb->Fill(ParticleInfo1[i][0]);
                }
            }
            else if(RealSeparation < i_MinAsymmDistance && RealSeparation > i_MinAsymmDistance-1.0)
            {
                AsymmetricEventCounter += 1;
                //fill asymmetric histogram
                for(Int_t i=0; i<ParticlesInJet1+1; i++)
                {
                    if(k==ENUM_pp) h1D_AsymmDist_pp->Fill(ParticleInfo1[i][0]);
                    else h1D_AsymmDist_pbpb->Fill(ParticleInfo1[i][0]);
                }
            }

        }//EVENTLOOP END

        //subtract distributions from each other

        if(SymmetricEventCounter > 0 && k==ENUM_pp)h1D_SymmDist_pp->Scale((double)1/(double)SymmetricEventCounter);
        else if(SymmetricEventCounter > 0 && k==ENUM_pbpb)h1D_SymmDist_pbpb->Scale((double)1/(double)SymmetricEventCounter);
        if(AsymmetricEventCounter > 0 && k==ENUM_pp)h1D_AsymmDist_pp->Scale((double)1/(double)AsymmetricEventCounter);
        else if(AsymmetricEventCounter > 0 && k==ENUM_pbpb)h1D_AsymmDist_pbpb->Scale((double)1/(double)AsymmetricEventCounter);

        if(k==ENUM_pp)
        {
            h1D_Difference_pp = (TH1D*)h1D_AsymmDist_pp->Clone("Difference in pp");
            h1D_Difference_pp->Add(h1D_SymmDist_pp, -1);
        }
        else
        {
            h1D_Difference_pbpb = (TH1D*)h1D_AsymmDist_pbpb->Clone("Difference in pbpb");
            h1D_Difference_pbpb->Add(h1D_SymmDist_pbpb, -1);
        }

    }

    //TCanvas* can_Difference = new TCanvas("Difference","Difference",2500,1000);

    /* //REFERENCE PLOT FROM CMS PAPER
    const int n = 14;
    double x[n] = {
    -2.37, -2.07, -1.81, -1.59, -1.19, -0.796, -0.418, -0.209,
    4.44e-16, 0.398, 0.786, 1.3, 1.65, 2.0
    };
    double y[n] = {
    0.183, -0.00469, -0.296, -0.596, -0.737, -0.596, -0.183,
    0.164, 0.502, 0.559, 0.202, 0.0235, 0.00469, -0.0704
    };

    TGraph *ReferenceGraph = new TGraph(n, x, y); */
    TCanvas* can_Comparison = new TCanvas("Comparison","Comparison",2500,1000);

    h1D_Difference_pbpb->SetTitle("");
    h1D_Difference_pbpb->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
    h1D_Difference_pbpb->GetXaxis()->SetLabelSize(DEF_AxisLabelSize);
    h1D_Difference_pbpb->GetXaxis()->SetTitle("#Delta #eta^{ch, jet_{1}}");
    h1D_Difference_pbpb->GetXaxis()->SetRangeUser(-i_Acceptance, i_Acceptance);
    h1D_Difference_pbpb->GetYaxis()->SetTitleSize(DEF_AxisLabelSize);
    h1D_Difference_pbpb->GetYaxis()->SetLabelSize(DEF_AxisLabelSize);
    h1D_Difference_pbpb->GetYaxis()->SetTitle("R^{Asym} - R^{Symm}");
    h1D_Difference_pbpb->SetMarkerColor(kRed);
    h1D_Difference_pbpb->SetMarkerStyle(kFullSquare);
    TH1* pbpb = h1D_Difference_pbpb->DrawCopy("P");
    h1D_Difference_pp->SetTitle("p-p");
    h1D_Difference_pp->SetMarkerColor(kBlue);
    h1D_Difference_pp->SetMarkerStyle(kFullCircle);
    TH1* pp = (TH1*)h1D_Difference_pp->DrawCopy("P SAME");

    TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
    leg->SetTextSize(DEF_LatexTextSize);
    leg->AddEntry(pp,"pp","p");
    leg->AddEntry(pbpb,"pbpb","p");
    leg->Draw();

    TLatex latex_diff;
    latex_diff.SetNDC();
    latex_diff.SetTextSize(DEF_LatexTextSize);
    latex_diff.SetTextAlign(31);
    latex_diff.DrawLatex(DEF_InfoBoxX, DEF_InfoBoxY, Form("Acceptance = %.1f", i_Acceptance));
    latex_diff.DrawLatex(DEF_InfoBoxX, DEF_InfoBoxY - DEF_InfoBoxLineDist, "Event selection like mentioned");
    latex_diff.DrawLatex(DEF_InfoBoxX, DEF_InfoBoxY - 2*DEF_InfoBoxLineDist, "1e7 events");
    latex_diff.DrawLatex(DEF_InfoBoxX, DEF_InfoBoxY - 3*DEF_InfoBoxLineDist, Form("#sigma = %.2f", i_Sigma));
    latex_diff.DrawLatex(DEF_InfoBoxX, DEF_InfoBoxY - 4*DEF_InfoBoxLineDist, Form("N^{pp} #in {%u, ..., %u }, N^{pbpb} #in {%u, ..., %u }", i_MinParticles_pp, i_MaxParticles_pp, i_MinParticles_pbpb, i_MaxParticles_pbpb));



}

/* void ToyModel(double i_Acceptance, double i_Range, Int_t i_Events, double i_Sigma, Int_t i_MinParticles_pp, Int_t i_MaxParticles_pp, Int_t i_MinParticles_pbpb, Int_t i_MaxParticles_pbpb, double MinAsymmDistance, Int_t i_Scale)
{

    #define Bins 300

    TRandom TR_JetPos1;
    TRandom TR_JetPos2;
    TRandom TR_PNum1;
    TRandom TR_PNum2;
    TR_JetPos1.SetSeed(0);
    TR_JetPos2.SetSeed(0);
    TR_PNum1.SetSeed(0);
    TR_PNum2.SetSeed(0);

    //get DeltaEta Distribution
    TH1D* h1D_DeltaEtaDist = new TH1D("DeltaEtaDist", "DeltaEtaDist", Bins, -i_Range, i_Range);
    for(Int_t i=0; i<i_Events; i++)
    {
        double Jet1Pos = TR_JetPos1.Uniform(-i_Range, i_Range);
        double Jet2Pos = TR_JetPos2.Uniform(-i_Range, i_Range);

        if(abs(Jet1Pos) > i_Acceptance || abs(Jet2Pos) > i_Acceptance) continue;

        h1D_DeltaEtaDist->Fill(Jet1Pos - Jet2Pos);
    }
    



    #define AtLeastParticles 2
    #define MaxParticles 60

    Int_t SymmetricEventCounter = 0;
    Int_t AsymmetricEventCounter = 0;

    TRandom TR_JetSeparation;
    TRandom TR_ParticlesInJet;
    TRandom TR_X_Particle1;
    TRandom TR_X_Particle2;
    TR_JetSeparation.SetSeed(0);
    TR_ParticlesInJet.SetSeed(0);
    TR_X_Particle1.SetSeed(0);
    TR_X_Particle2.SetSeed(0);

    //distribution on near-side
    TH1D* h1D_SymmDist = new TH1D("SymmDist", "SymmDist", Bins, -i_Acceptance, i_Acceptance);
    TH1D* h1D_AsymmDist = new TH1D("AsymmDist", "AsymmDist", Bins, -i_Acceptance, i_Acceptance);
    TH1D* h1D_Difference = new TH1D("Difference", "Difference", Bins, -i_Acceptance, i_Acceptance);

    for(Int_t i=0; i<i_Events; i++)
    {
        //SEPARATION IS X_SUBLEADING - X_LEADING

        //get separation of the two jets by probability density
        double Separation = h1D_DeltaEtaDist->GetRandom();
        if(Separation > 0) Separation *= -1;
        if(k==ENUM_pp)
        {
            Int_t ParticlesInJet1 = TR_ParticlesInJet.Integer(i_MaxParticles - i_MinParticles) + i_MinParticles;
            Int_t ParticlesInJet2 = TR_ParticlesInJet.Integer(i_MaxParticles - i_MinParticles) + i_MinParticles;
        }

        //produce two realistic discrete almost-gaussian distributions
        //JET1
        double x_Particle1[MaxParticles];
        for(Int_t i=0; i<ParticlesInJet1; i++)
        {
            x_Particle1[i] = TR_X_Particle1.Gaus(0, i_Sigma);
        }
        double cg1 = 0;
        for(Int_t i=0; i<ParticlesInJet1; i++)
        {
            cg1 += x_Particle1[i];
        }
        cg1 /= ParticlesInJet1;
        //JET2
        double x_Particle2[MaxParticles];
        for(Int_t i=0; i<ParticlesInJet2; i++)
        {
            x_Particle2[i] = TR_X_Particle2.Gaus(Separation, i_Sigma);
        }
        double cg2 = 0;
        for(Int_t i=0; i<ParticlesInJet2; i++)
        {
            cg2 += x_Particle2[i];
        }
        cg2 /= ParticlesInJet2;

        //categorize into the two cases and write
        double RealSeparation = cg2 - cg1;
        if(RealSeparation > -0.5)
        {
            SymmetricEventCounter += 1;
            //fill symmetrichistogram
            for(Int_t i=0; i<ParticlesInJet1; i++)
            {
                h1D_SymmDist->Fill(x_Particle1[i]);
            }
        }
        else if(RealSeparation < MinAsymmDistance && RealSeparation > MinAsymmDistance-1.0)
        {
            AsymmetricEventCounter += 1;
            //fill asymmetric histogram
            for(Int_t i=0; i<ParticlesInJet1; i++)
            {
                h1D_AsymmDist->Fill(x_Particle1[i]);
            }
        }

    }//EVENTLOOP END

    //subtract distributions from each other

    h1D_SymmDist->Scale((double)1/(double)SymmetricEventCounter);
    h1D_AsymmDist->Scale((double)1/(double)AsymmetricEventCounter);

    h1D_Difference->Add(h1D_AsymmDist);
    h1D_Difference->Add(h1D_SymmDist, -1);


    TCanvas* can_Difference = new TCanvas("Difference","Difference",2500,1000);

    //REFERENCE PLOT
    const int n = 14;
    double x[n] = {
    -2.37, -2.07, -1.81, -1.59, -1.19, -0.796, -0.418, -0.209,
    4.44e-16, 0.398, 0.786, 1.3, 1.65, 2.0
    };
    double y[n] = {
    0.183, -0.00469, -0.296, -0.596, -0.737, -0.596, -0.183,
    0.164, 0.502, 0.559, 0.202, 0.0235, 0.00469, -0.0704
    };

    TGraph *ReferenceGraph = new TGraph(n, x, y);

    h1D_Difference->Scale(i_Scale);
    h1D_Difference->Draw();
    ReferenceGraph->Draw("SAME");

} */


//JetEventNo gives back the ith event in which a jet with p_t > DEF_HighPtLimit was found (0 is the first one)
/* void Disprove1(long NoEvents, double Cutoff, Int_t AtLeastParticles, double Sigma)
{

    #define Bins 500
    #define HistoRange 5

    TRandom TR_ParticlesInJet;
    TRandom TR_X_Particle;
    TRandom TR_Mean;
    TR_ParticlesInJet.SetSeed(0);

    //histograms
    TH1D* h1D_WithCutoff = new TH1D("h1D_WithCutoff", "h1D_WithCutoff", Bins, -HistoRange, HistoRange);
    TH1D* h1D_WithoutCutoff = new TH1D("h1D_WithoutCutoff", "h1D_WithoutCutoff", Bins, -HistoRange, HistoRange);

    for(long i_Event=0; i_Event<NoEvents; i_Event++)
    {
        //generate random jet position
        double Mean = TR_Mean.Uniform(-10, 0);

        //generate random number of particles
        Int_t ParticlesInJet = TR_ParticlesInJet.Integer(16) + AtLeastParticles;

        //generate jet particles
        double x_Particle[ParticlesInJet];
        for(Int_t i=0; i<ParticlesInJet; i++)
        {
            x_Particle[i] = TR_X_Particle.Gaus(Mean, Sigma);
        }

        //get jet center of gravity
        double cg = 0;
        for(Int_t i=0; i<ParticlesInJet; i++)
        {
            cg += x_Particle[i];
        }
        cg /= ParticlesInJet;

        if(cg <=Cutoff)
        {
            for(Int_t i=0; i<ParticlesInJet; i++)
            {
                h1D_WithCutoff->Fill(x_Particle[i]);
            }
        }

        for(Int_t i=0; i<ParticlesInJet; i++)
        {
            h1D_WithoutCutoff->Fill(x_Particle[i]);
        }
        
    }

    TCanvas* can_DisproveCMS = new TCanvas("DisproveCMS","DisproveCMS",2500,1000);
    can_DisproveCMS->Divide(3);

    can_DisproveCMS->cd(1);
    h1D_WithCutoff->Draw();
    can_DisproveCMS->cd(2);
    h1D_WithoutCutoff->Draw();
    
    //subtract distributions
    TH1D* h1D_Difference = new TH1D("Difference", "Difference", Bins, -HistoRange, HistoRange);
    h1D_Difference->Add(h1D_WithCutoff, 1);
    h1D_Difference->Add(h1D_WithoutCutoff, -1);
    can_DisproveCMS->cd(3);
    double Difference = h1D_Difference->Integral();
    double NoOfParticles = h1D_WithoutCutoff->Integral();
    double RelDif = Difference / NoOfParticles;
    h1D_Difference->SetTitle(Form("Difference (Integrated relative diff: %.3f)",RelDif));
    h1D_Difference->Draw();
} */