#include "./StJetTrackEvent.h"
#include "./StJetTrackEventLinkDef.h"
#include "./functions.h"
#include "TPaveText.h"
#include <algorithm>
using namespace std;
#if false
#include "fastjet/config.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/Subtractor.hh"
#include "StPhysicalHelixD.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#endif

Int_t Analysis_SimulationStats()
{
    #define DEF_AxisLabelSize 0.06
    #define DEF_AxisTitleSize 0.06
    #define DEF_HistoTitleSize 0.1
    #define DEF_LegendTextSize 0.06
    #define DEF_Rebin 16

    gStyle->SetOptStat(0);
    SetRootGraphicStyle();

    //open the root file
    TFile *RUN2File = TFile::Open("/home/nicolas-wirth/Dokumente/BA/Misc/Elliptical_Flow_Bias/merge_mult_track_pT.root");

    if (!RUN2File || RUN2File->IsZombie()) {
        std::cout << "Error while opening the file!" << std::endl;
        return 0;
    }

    TH1D* h1D_NDist = (TH1D*)RUN2File->Get("h_mult_particles_used_SE_ME_0");
    if(!h1D_NDist){cout << "Multiplicity histogram not found!" << endl; return 0;}
    TH1D* h1D_PtDist = (TH1D*)RUN2File->Get("h_particle_pT");
    if(!h1D_PtDist){cout << "Pt histogram not found!" << endl; return 0;}
    TProfile* prof_Pt_v2 = (TProfile*)RUN2File->Get("h_profile_v2_track_SE");
    if(!prof_Pt_v2){cout << "v_2(p_T) histogram not found!" << endl; return 0;}
    
    cout << "Retrieved Histograms!" << endl;

    //configure and draw canvas
    TCanvas* Canvas = new TCanvas("Canvas","Canvas",1000,1000);
    Canvas->Divide(1,4, 0, 0);

    //Multiplicity
    Canvas->cd(1);
    h1D_NDist->SetTitle("");
    h1D_NDist->GetXaxis()->SetTitleSize(DEF_AxisTitleSize);
    h1D_NDist->GetXaxis()->SetLabelSize(DEF_AxisLabelSize);
    h1D_NDist->GetXaxis()->SetTitle("Multiplicity");
    h1D_NDist->GetXaxis()->SetRangeUser(1600, 3500);
    h1D_NDist->GetYaxis()->SetTitleSize(DEF_AxisTitleSize);
    h1D_NDist->GetYaxis()->SetLabelSize(DEF_AxisLabelSize);
    h1D_NDist->GetYaxis()->SetTitle("a. u.");
    h1D_NDist->Rebin(32);
    h1D_NDist->Scale(1./h1D_NDist->Integral());
    h1D_NDist->SetMarkerStyle(kStar);
    h1D_NDist->SetMarkerColor(kBlue);
    h1D_NDist->DrawCopy("P E");

    TLegend *leg1 = new TLegend(0.0,0.8,0.3,1.0);
    leg1->SetTextSize(DEF_LegendTextSize);
    leg1->AddEntry((TObject*)nullptr,Form("0-10%% ALICE Pb-Pb"),"");
    leg1->AddEntry((TObject*)nullptr,Form("#sqrt{s_{NN}} = 5.02 TeV, |v_{z}| < 8 cm"),"");
    leg1->AddEntry(h1D_NDist,"Data","p");
    leg1->SetLineColor(10);
    leg1->Draw();

    //p_T
    Canvas->cd(2);
    h1D_PtDist->SetTitle("");
    h1D_PtDist->GetXaxis()->SetTitleSize(DEF_AxisTitleSize);
    h1D_PtDist->GetXaxis()->SetLabelSize(DEF_AxisLabelSize);
    h1D_PtDist->GetXaxis()->SetTitle("p_{T} [GeV]");
    h1D_PtDist->GetYaxis()->SetTitleSize(DEF_AxisTitleSize);
    h1D_PtDist->GetYaxis()->SetLabelSize(DEF_AxisLabelSize);
    h1D_PtDist->GetYaxis()->SetTitle("a. u.");
    h1D_PtDist->Rebin(32);
    h1D_PtDist->Scale(1./h1D_PtDist->Integral());
    h1D_PtDist->SetMarkerStyle(kStar);
    h1D_PtDist->SetMarkerColor(kBlue);
    h1D_PtDist->DrawCopy("P E");

    TLegend *leg2 = new TLegend(0.0,0.8,0.3,1.0);
    leg2->SetTextSize(DEF_LegendTextSize);
    leg2->AddEntry((TObject*)nullptr,Form("0-10%% ALICE Pb-Pb"),"");
    leg2->AddEntry((TObject*)nullptr,Form("#sqrt{s_{NN}} = 5.02 TeV, |v_{z}| < 8 cm"),"");
    leg2->AddEntry(h1D_PtDist,"Data","p");
    leg2->SetLineColor(10);
    leg2->Draw();

    //v_2(p_T)
    Canvas->cd(3);

    Int_t n = 40;
    Double_t x[n], y[n];
    for (Int_t i=0;i<n;i++) {
        x[i] = 0.5*i;
        y[i] = 0.063*pow(x[i],1.6)*exp(-0.47*x[i]);
    }
    // create graph
    TGraph *Graph_v2_pT  = new TGraph(n,x,y);
    // draw the graph with axis, continuous line, and put

    Graph_v2_pT->SetTitle("");
    Graph_v2_pT->GetXaxis()->SetTitleSize(DEF_AxisTitleSize);
    Graph_v2_pT->GetXaxis()->SetLabelSize(DEF_AxisLabelSize);
    Graph_v2_pT->GetXaxis()->SetTitle("p_{T} [GeV]");
    Graph_v2_pT->GetXaxis()->SetRangeUser(0,20);
    Graph_v2_pT->GetYaxis()->SetTitleSize(DEF_AxisTitleSize);
    Graph_v2_pT->GetYaxis()->SetLabelSize(DEF_AxisLabelSize);
    Graph_v2_pT->GetYaxis()->SetTitle("v_{2}");
    Graph_v2_pT->Draw("AC");
    prof_Pt_v2->SetMarkerStyle(kStar);
    prof_Pt_v2->SetMarkerColor(kBlue);
    prof_Pt_v2->DrawCopy("SAME P E");

    TLegend *leg3 = new TLegend(0.0,0.8,0.3,1.0);
    leg3->SetTextSize(DEF_LegendTextSize);
    leg3->AddEntry((TObject*)nullptr,Form("0-10%% ALICE Pb-Pb"),"");
    leg3->AddEntry((TObject*)nullptr,Form("#sqrt{s_{NN}} = 5.02 TeV, |v_{z}| < 8 cm"),"");
    leg3->AddEntry(prof_Pt_v2,"Data","p");
    leg3->AddEntry(Graph_v2_pT,"Fit","l");
    leg3->SetLineColor(10);
    leg3->Draw();
    

    //v_2(eta)
    Canvas->cd(4);
    
    Double_t OrigX[6] = {-1.25, -0.75, -0.25, 0.25, 0.75, 1.25};
    Double_t OrigY[6] = {0.0213, 0.0235, 0.0251, 0.0253, 0.0231, 0.0213};
    Double_t eX[6] = {0};
    Double_t eY[6] = {2.01e-3, 2.38e-3, 2.642e-3, 2.67e-3, 2.281e-3, 2.06e-3};
    TGraphErrors* Graph_Orig = new TGraphErrors(6, OrigX, OrigY, eX, eY);
    Graph_Orig->SetMarkerColor(kBlue);
    Graph_Orig->SetTitle("");
    Graph_Orig->GetXaxis()->SetTitleSize(DEF_AxisTitleSize);
    Graph_Orig->GetXaxis()->SetLabelSize(DEF_AxisLabelSize);
    Graph_Orig->GetXaxis()->SetTitle("#eta");
    Graph_Orig->GetXaxis()->SetRangeUser(-1.5, 1.5);
    Graph_Orig->GetYaxis()->SetTitleSize(DEF_AxisTitleSize);
    Graph_Orig->GetYaxis()->SetLabelSize(DEF_AxisLabelSize);
    Graph_Orig->GetYaxis()->SetTitle("v_{2}");
    Graph_Orig->SetLineColor(kBlue);
    Graph_Orig->Draw("A*");

    auto ParabFit = new TF1("parab", "[0] - [1]*sqrt(x^2)", -1.5, 1.5);
    // create graph
    Graph_Orig->Fit("parab");
    gPad->Update();
    TF1 *fit = (TF1*)Graph_Orig->GetListOfFunctions()->FindObject("parab");
    fit->SetLineColor(kBlack);

    TLegend *leg4 = new TLegend(0.0,0.8,0.3,1.0);
    leg4->SetTextSize(DEF_LegendTextSize);
    leg4->AddEntry((TObject*)nullptr,Form("0-5%% ALICE Pb-Pb"),"");
    leg4->AddEntry((TObject*)nullptr,Form("#sqrt{s_{NN}} = 2.76 TeV, |v_{z}| < 10 cm"),"");
    leg4->AddEntry((TObject*)nullptr,Form("0.2 #leq p_{T}^{track} #leq 5.0 GeV"),"");
    leg4->AddEntry(Graph_Orig,"Data","p");
    leg4->AddEntry(fit,"Fit","l");
    leg4->SetLineColor(10);
    leg4->Draw();

    return 1;
}

Int_t DijetAna(const TString DataFile, float DiffRange) {

    #define DEF_AxisLabelSize 0.07
    #define DEF_AxisTitleSize 0.08
    #define DEF_HistoTitleSize 0.1
    #define DEF_Rebin 25
    #define DEF_LegendFontSize 0.055
    #define DEF_Margin_Top 0.05
    #define DEF_Margin_Bottom 0.23
    #define DEF_Margin_Left 0.18
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
    TH1F* h1F_SmallGap[4];
    h1F_SmallGap[0] = (TH1F*)file->Get("h1F_1DCorrelation_eta_pT_0_1_SmallGap");
    h1F_SmallGap[1] = (TH1F*)file->Get("h1F_1DCorrelation_eta_pT_1_2_SmallGap");
    h1F_SmallGap[2] = (TH1F*)file->Get("h1F_1DCorrelation_eta_pT_2_4_SmallGap");
    h1F_SmallGap[3] = (TH1F*)file->Get("h1F_1DCorrelation_eta_pT_4_6_SmallGap");

    TH1F* h1F_LargeGap[4];
    h1F_LargeGap[0] = (TH1F*)file->Get("h1F_1DCorrelation_eta_pT_0_1_LargeGap");
    h1F_LargeGap[1] = (TH1F*)file->Get("h1F_1DCorrelation_eta_pT_1_2_LargeGap");
    h1F_LargeGap[2] = (TH1F*)file->Get("h1F_1DCorrelation_eta_pT_2_4_LargeGap");
    h1F_LargeGap[3] = (TH1F*)file->Get("h1F_1DCorrelation_eta_pT_4_6_LargeGap");

    for(int iRow=0; iRow<=3; iRow++)
    {
        if(!h1F_LargeGap[iRow] || !h1F_SmallGap[iRow])
        {
            cout << "Histograms not found" << endl;
            return 0;
        }
    }

    //configure and draw canvas
    TCanvas* can_ParticleMultiplicities = new TCanvas("ParticleMultiplicities","ParticleMultiplicities",1000,1000);
    can_ParticleMultiplicities->Divide(2,4, 0, 0);

    for(int iRow=0; iRow<=3; iRow++)
    {
        //ABSOLUTE DISTROS
        can_ParticleMultiplicities->cd(1 + (2*iRow));

        gPad->SetMargin(DEF_Margin_Left, DEF_Margin_Right_ForMiddlePanel, DEF_Margin_Bottom*(iRow==3), DEF_Margin_Top*(iRow==0));
        //markings
        h1F_SmallGap[iRow]->SetTitle("");
        h1F_SmallGap[iRow]->SetTitleSize(DEF_HistoTitleSize);
        h1F_SmallGap[iRow]->GetXaxis()->SetTitleSize(DEF_AxisTitleSize);
        h1F_SmallGap[iRow]->GetYaxis()->SetTitleSize(DEF_AxisTitleSize);
        h1F_SmallGap[iRow]->GetXaxis()->SetLabelSize(DEF_AxisLabelSize);
        h1F_SmallGap[iRow]->GetYaxis()->SetLabelSize(DEF_AxisLabelSize);
        (iRow==0) ? (h1F_SmallGap[iRow]->GetYaxis()->SetTitle("C_{Large/Small Gap}")) : (h1F_SmallGap[iRow]->GetYaxis()->SetTitle(""));
        h1F_SmallGap[iRow]->Rebin(DEF_Rebin);
        h1F_SmallGap[iRow]->Scale(1./DEF_Rebin);
        h1F_SmallGap[iRow]->SetMarkerStyle(kCircle);
        h1F_SmallGap[iRow]->SetMarkerColor(kBlue);
        TH1* SmallCopy = (TH1*)h1F_SmallGap[iRow]->DrawCopy();

        h1F_LargeGap[iRow]->SetTitle("");
        h1F_LargeGap[iRow]->SetTitleSize(DEF_HistoTitleSize);
        h1F_LargeGap[iRow]->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
        h1F_LargeGap[iRow]->GetYaxis()->SetTitleSize(DEF_AxisLabelSize);
        h1F_LargeGap[iRow]->GetXaxis()->SetLabelSize(DEF_AxisLabelSize);
        h1F_LargeGap[iRow]->GetYaxis()->SetLabelSize(DEF_AxisLabelSize);
        h1F_LargeGap[iRow]->Rebin(DEF_Rebin);
        h1F_LargeGap[iRow]->Scale(1./DEF_Rebin);
        h1F_LargeGap[iRow]->SetMarkerStyle(kCircle);
        h1F_LargeGap[iRow]->SetMarkerColor(kRed);
        TH1* LargeCopy = (TH1*)h1F_LargeGap[iRow]->DrawCopy("SAME");


        TLegend *leg = new TLegend(0.25,0.77,0.26,0.93);
        leg->AddEntry((TObject*)nullptr,"ALICE Pb-Pb 0-10%  #sqrt{s_{NN}}=5.02 TeV","");
        leg->AddEntry((TObject*)nullptr,Form("|#phi - #phi_{Leading}| #leq #pi/2, %.1f #leq p_{T} #leq %.1f GeV", LowPtCuts[iRow], HighPtCuts[iRow]),"");
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
        
        h1F_LargeGap[iRow]->Add(h1F_SmallGap[iRow], -1);

        h1F_LargeGap[iRow]->SetTitle("");
        (iRow==0) ? (h1F_LargeGap[iRow]->GetYaxis()->SetTitle("C_{Large Gap} - C_{Small Gap}")) : (h1F_LargeGap[iRow]->GetYaxis()->SetTitle(""));
        h1F_LargeGap[iRow]->SetTitleSize(DEF_HistoTitleSize);
        (iRow == 3) ? (h1F_LargeGap[iRow]->GetXaxis()->SetTitle("#eta")) : (h1F_LargeGap[iRow]->GetXaxis()->SetTitle(""));
        h1F_LargeGap[iRow]->GetXaxis()->SetTitleSize(DEF_AxisTitleSize);
        h1F_LargeGap[iRow]->GetYaxis()->SetTitleSize(DEF_AxisTitleSize);
        h1F_LargeGap[iRow]->GetXaxis()->SetLabelSize(DEF_AxisLabelSize);
        h1F_LargeGap[iRow]->GetYaxis()->SetLabelSize(DEF_AxisLabelSize);
        h1F_LargeGap[iRow]->SetMarkerStyle(kCircle);
        h1F_LargeGap[iRow]->SetMarkerColor(kBlack);
        h1F_LargeGap[iRow]->DrawCopy();
        h1F_LargeGap[iRow]->GetYaxis()->SetRangeUser(-DiffRange, DiffRange);
        h1F_LargeGap[iRow]->DrawCopy();
    }

    return 1;

}