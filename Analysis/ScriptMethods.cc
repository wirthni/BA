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
    #define DEF_AxisLabelSize 0.04
    #define DEF_AxisTitleSize 0.05
    #define DEF_HistoTitleSize 0.1
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
    Canvas->Divide(2,2, 0, 0);

    //Multiplicity
    Canvas->cd(1);
    h1D_NDist->SetTitle("");
    h1D_NDist->GetXaxis()->SetTitleSize(DEF_AxisTitleSize);
    h1D_NDist->GetXaxis()->SetLabelSize(DEF_AxisLabelSize);
    h1D_NDist->GetXaxis()->SetTitle("Multiplicity");
    h1D_NDist->GetYaxis()->SetTitleSize(DEF_AxisTitleSize);
    h1D_NDist->GetYaxis()->SetLabelSize(DEF_AxisLabelSize);
    h1D_NDist->GetYaxis()->SetTitle("a. u.");
    h1D_NDist->DrawCopy();

    //p_T
    Canvas->cd(2);
    h1D_PtDist->SetTitle("");
    h1D_PtDist->GetXaxis()->SetTitleSize(DEF_AxisTitleSize);
    h1D_PtDist->GetXaxis()->SetLabelSize(DEF_AxisLabelSize);
    h1D_PtDist->GetXaxis()->SetTitle("p_{T} [GeV]");
    h1D_PtDist->GetYaxis()->SetTitleSize(DEF_AxisTitleSize);
    h1D_PtDist->GetYaxis()->SetLabelSize(DEF_AxisLabelSize);
    h1D_PtDist->GetYaxis()->SetTitle("a. u.");
    h1D_PtDist->DrawCopy();

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
    prof_Pt_v2->DrawCopy("SAME HIST");
    

    //v_2(eta)
    Canvas->cd(4);
    
    Double_t OrigX[6] = {-1.25, -0.75, -0.25, 0.25, 0.75, 1.25};
    Double_t OrigY[6] = {0.0213, 0.0235, 0.0251, 0.0253, 0.0231, 0.0213};
    TGraph* Graph_Orig = new TGraph(6, OrigX, OrigY);
    Graph_Orig->SetMarkerColor(kBlue);
    Graph_Orig->SetTitle("");
    Graph_Orig->GetXaxis()->SetTitleSize(DEF_AxisTitleSize);
    Graph_Orig->GetXaxis()->SetLabelSize(DEF_AxisLabelSize);
    Graph_Orig->GetXaxis()->SetTitle("#eta");
    Graph_Orig->GetXaxis()->SetRangeUser(-1.5, 1.5);
    Graph_Orig->GetYaxis()->SetTitleSize(DEF_AxisTitleSize);
    Graph_Orig->GetYaxis()->SetTitle("v_{2}(p_{T},#eta)/v_{2}(p_{T})");
    Graph_Orig->SetLineColor(kBlue);
    Graph_Orig->Draw("A*");

    auto ParabFit = new TF1("parab", "[0] - [1]*sqrt(x^2)", -1.5, 1.5);
    // create graph
    Graph_Orig->Fit("parab");

    return 1;
}

void Analysis_ALICEDepletion(const TString DataFile, Int_t i_ProjectionBins, Int_t i_BinLow, Int_t i_BinMid, Int_t i_BinHigh) {

    #define DEF_AxisLabelSize 0.05
    #define DEF_HistoTitleSize 0.1
    #define DEF_InfoBoxX 0.95
    #define DEF_InfoBoxY 0.8
    #define DEF_LastColumnInfoXOffset -0.15
    #define DEF_InfoBoxLineDist 0.08
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
    can_JetRecoilSitesAndBg->Divide(3,2,0,0);

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

        gPad->SetMargin(0.12*(i==1) + 0.07, 0.04, 0.12, 0.1);

        //markings
        h1D_pt_vs_eta_Background_Projection->SetTitle("");
        h1D_pt_vs_eta_Background_Projection->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
        h1D_pt_vs_eta_Background_Projection->GetXaxis()->SetLabelSize(DEF_AxisLabelSize);
        if(i==3) h1D_pt_vs_eta_Background_Projection->GetXaxis()->SetTitle("#Delta #eta");
        else h1D_pt_vs_eta_Background_Projection->GetXaxis()->SetTitle("");
        if(i==1)h1D_pt_vs_eta_Background_Projection->GetYaxis()->SetTitle("Multiplicity");
        else h1D_pt_vs_eta_Background_Projection->GetYaxis()->SetTitle("");
        h1D_pt_vs_eta_Background_Projection->GetYaxis()->SetLabelSize(DEF_AxisLabelSize);
        h1D_pt_vs_eta_Background_Projection->GetYaxis()->SetTitleSize(DEF_AxisLabelSize);
        

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
        if(i==1)
        {
            TLegend *leg = new TLegend(0.0,0.8,0.3,1.0);
            leg->SetTextSize(DEF_AxisLabelSize);
            leg->AddEntry(hbg,"Background","p");
            leg->AddEntry(hall,"All jets","p");
            leg->Draw();
        }
        //draw additional info
        TLatex latex_dist;
        latex_dist.SetNDC();
        latex_dist.SetTextSize(DEF_LatexTextSize);
        latex_dist.SetTextAlign(31);
        latex_dist.DrawLatex(DEF_InfoBoxX, DEF_InfoBoxY, Form("p_{t}^{recoil} #in [%.1f,%.1f] #frac{GeV}{c}", LowestMomentum, HighestMomentum));
        latex_dist.DrawLatex(DEF_InfoBoxX, DEF_InfoBoxY - DEF_InfoBoxLineDist, "p_{t}^{jet} > 7 #frac{GeV}{c}");
        latex_dist.DrawLatex(DEF_InfoBoxX, DEF_InfoBoxY - 2*DEF_InfoBoxLineDist, "Pb-Pb, R = 0.4");

        can_JetRecoilSitesAndBg->cd(i+3);

        gPad->SetMargin(0.12*(i==1) + 0.07, 0.04, 0.12, 0.1);

        //get the differences
        TH1D* h1D_MultDif_vs_eta_All = (TH1D*)h1D_pt_vs_eta_recoil_All_Projection->Clone("h1D_MultDif_vs_eta_All");
        h1D_MultDif_vs_eta_All->Add(h1D_pt_vs_eta_Background_Projection, -1);

        h1D_MultDif_vs_eta_All->SetNdivisions(5);
        h1D_MultDif_vs_eta_All->GetXaxis()->SetLabelSize(DEF_AxisLabelSize);
        h1D_MultDif_vs_eta_All->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
        h1D_MultDif_vs_eta_All->GetYaxis()->SetLabelSize(DEF_AxisLabelSize);
        h1D_MultDif_vs_eta_All->GetYaxis()->SetTitleSize(DEF_AxisLabelSize);
        if(i==3) h1D_MultDif_vs_eta_All->GetXaxis()->SetTitle("#Delta #eta");
        else h1D_MultDif_vs_eta_All->GetXaxis()->SetTitle("");
        h1D_MultDif_vs_eta_All->SetTitle("");
        if(i==1)
        {
            h1D_MultDif_vs_eta_All->GetYaxis()->SetTitle("#Delta Multiplicity");
        }
        else
        {
            h1D_MultDif_vs_eta_All->GetYaxis()->SetTitle("");
        }
        h1D_MultDif_vs_eta_All->Rebin(8);
        h1D_MultDif_vs_eta_All->SetMarkerStyle(kFullCircle);
        h1D_MultDif_vs_eta_All->SetMarkerColor(kBlue);
        h1D_MultDif_vs_eta_All->DrawCopy("E P");

        //draw additional info
        TLatex latex_diff;
        latex_diff.SetNDC();
        latex_diff.SetTextSize(DEF_LatexTextSize);
        latex_diff.SetTextAlign(31);
        latex_diff.DrawLatex(DEF_InfoBoxX, DEF_InfoBoxY, Form("p_{t}^{recoil} #in [%.1f,%.1f] #frac{GeV}{c}", LowestMomentum, HighestMomentum));
        latex_diff.DrawLatex(DEF_InfoBoxX, DEF_InfoBoxY - DEF_InfoBoxLineDist, "p_{t}^{jet} > 7 #frac{GeV}{c}");
        latex_diff.DrawLatex(DEF_InfoBoxX, DEF_InfoBoxY - 2*DEF_InfoBoxLineDist, "Pb-Pb, R = 0.4");

    }

    cout << "DONE!" << endl;
    file->Close();
}

void Analysis_RecoilMultiplicitySingle(const TString DataFile)
{
    #define DEF_AxisLabelSize 0.05
    #define DEF_HistoTitleSize 0.1
    #define DEF_InfoBoxX 0.95
    #define DEF_InfoBoxY 0.9
    #define DEF_LastColumnInfoXOffset -0.15
    #define DEF_InfoBoxLineDist 0.1
    #define DEF_LatexTextSize 0.04

    TFile *file = TFile::Open(DataFile);
    if (!file || file->IsZombie()) {
        std::cout << "Error while opening the file!" << std::endl;
    }
    gStyle->SetTextSize(0.05);
    gStyle->SetOptStat(0);
    SetRootGraphicStyle();

    TH2D* h2D_ptvseta[3];

    h2D_ptvseta[0] = (TH2D*)file->Get("Pt_vs_Eta/h2D_pt_vs_eta_recoil_Summed_LowPt");
    h2D_ptvseta[1] = (TH2D*)file->Get("Pt_vs_Eta/h2D_pt_vs_eta_recoil_Summed_MiddlePt");
    h2D_ptvseta[2] = (TH2D*)file->Get("Pt_vs_Eta/h2D_pt_vs_eta_recoil_Summed_HighPt");

    h2D_ptvseta[0]->Add(h2D_ptvseta[1], 1);
    h2D_ptvseta[0]->Add(h2D_ptvseta[2], 1);

    TCanvas* can_Distributions = new TCanvas("Jet and Recoil Site Multiplicities","Jet and Recoil Site Multiplicities",1000,1000);

    h2D_ptvseta[0]->SetTitle("");
    h2D_ptvseta[0]->GetYaxis()->SetRangeUser(0, 60);
    h2D_ptvseta[0]->GetXaxis()->SetRangeUser(-0.5, 0.5);

    h2D_ptvseta[0]->Scale(1./h2D_ptvseta[0]->GetEntries());

    h2D_ptvseta[0]->GetYaxis()->SetTitle("p_{t} [#frac{GeV}{c}]");
    h2D_ptvseta[0]->GetYaxis()->SetTitleSize(DEF_AxisLabelSize);
    h2D_ptvseta[0]->GetYaxis()->SetLabelSize(DEF_AxisLabelSize);

    h2D_ptvseta[0]->GetXaxis()->SetTitle("#Delta #eta");
    h2D_ptvseta[0]->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
    h2D_ptvseta[0]->GetXaxis()->SetLabelSize(DEF_AxisLabelSize);

    h2D_ptvseta[0]->Draw("COLZ");

    TLatex latex_eta;
    latex_eta.SetNDC();
    latex_eta.SetTextSize(DEF_LatexTextSize);
    latex_eta.SetTextAlign(31);

    latex_eta.DrawLatex(DEF_InfoBoxX, DEF_InfoBoxY, "#frac{1}{N} #frac{d^{2}N}{d#eta dp_{t}}");

}


void Analysis_RecoilMultiplicityDistribution(const TString DataFile)
{
    #define DEF_AxisLabelSize 0.05
    #define DEF_HistoTitleSize 0.1
    #define DEF_InfoBoxX 0.95
    #define DEF_InfoBoxY 0.9
    #define DEF_LastColumnInfoXOffset -0.15
    #define DEF_InfoBoxLineDist 0.1
    #define DEF_LatexTextSize 0.04

    TFile *file = TFile::Open(DataFile);
    if (!file || file->IsZombie()) {
        std::cout << "Error while opening the file!" << std::endl;
    }
    gStyle->SetTextSize(0.05);
    gStyle->SetOptStat(0);
    SetRootGraphicStyle();

    TH2D* h2D_ptvseta[3];
    TH2D* h2D_ptvsphi[3];

    h2D_ptvseta[0] = (TH2D*)file->Get("Pt_vs_Eta/h2D_pt_vs_eta_recoil_Summed_LowPt");
    h2D_ptvseta[1] = (TH2D*)file->Get("Pt_vs_Eta/h2D_pt_vs_eta_recoil_Summed_MiddlePt");
    h2D_ptvseta[2] = (TH2D*)file->Get("Pt_vs_Eta/h2D_pt_vs_eta_recoil_Summed_HighPt");

    h2D_ptvsphi[0] = (TH2D*)file->Get("Pt_vs_Phi/h2D_pt_vs_phi_recoil_Summed_LowPt");
    h2D_ptvsphi[1] = (TH2D*)file->Get("Pt_vs_Phi/h2D_pt_vs_phi_recoil_Summed_MiddlePt");
    h2D_ptvsphi[2] = (TH2D*)file->Get("Pt_vs_Phi/h2D_pt_vs_phi_recoil_Summed_HighPt");

    TCanvas* can_Distributions = new TCanvas("Jet and Recoil Site Multiplicities","Jet and Recoil Site Multiplicities",1000,1000);
    can_Distributions->Divide(3,2,0,0);

    for(Int_t i=1; i<=3; i++)
    {
        // -------- ETA ROW --------
        can_Distributions->cd(i);
        gPad->SetMargin(0.16*(i==1), 0.15*(i==3), 0.15, 0.02);

        h2D_ptvseta[i-1]->SetTitle("");
        h2D_ptvseta[i-1]->GetYaxis()->SetRangeUser(0, 60);
        h2D_ptvseta[i-1]->GetXaxis()->SetRangeUser(-0.5, 0.5);

        h2D_ptvseta[i-1]->Scale(1./h2D_ptvseta[i-1]->GetEntries());

        // Axis handling
        if(i==1){
            h2D_ptvseta[i-1]->GetYaxis()->SetTitle("p_{t} [#frac{GeV}{c}]");
            h2D_ptvseta[i-1]->GetYaxis()->SetTitleSize(DEF_AxisLabelSize);
        } else {
            h2D_ptvseta[i-1]->GetYaxis()->SetTitle("");
        }
        h2D_ptvseta[i-1]->GetYaxis()->SetLabelSize(DEF_AxisLabelSize);

        if(i<3){
            h2D_ptvseta[i-1]->GetXaxis()->SetTitle("");
        } else {
            h2D_ptvseta[i-1]->GetXaxis()->SetTitle("#Delta #eta");
            h2D_ptvseta[i-1]->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
        }
        h2D_ptvseta[i-1]->GetXaxis()->SetLabelSize(DEF_AxisLabelSize);


        if(i==3)h2D_ptvseta[i-1]->Draw("COLZ");
        else h2D_ptvseta[i-1]->Draw();

        TLatex latex_eta;
        latex_eta.SetNDC();
        latex_eta.SetTextSize(DEF_LatexTextSize);
        latex_eta.SetTextAlign(31);

        if(i==1){
            latex_eta.DrawLatex(DEF_InfoBoxX, DEF_InfoBoxY, "log#left(#frac{1}{N_{low}} #frac{d^{2}N}{d#eta dp_{t}}#right)");
            latex_eta.DrawLatex(DEF_InfoBoxX, DEF_InfoBoxY - DEF_InfoBoxLineDist, "Leading p_{t} #in [30, 50) #frac{GeV}{c}");
        } else if(i==2){
            latex_eta.DrawLatex(DEF_InfoBoxX, DEF_InfoBoxY, "log#left(#frac{1}{N_{middle}} #frac{d^{2}N}{d#eta dp_{t}}#right)");
            latex_eta.DrawLatex(DEF_InfoBoxX, DEF_InfoBoxY - DEF_InfoBoxLineDist, "Leading p_{t} #in [50, 80) #frac{GeV}{c}");
        } else {
            latex_eta.DrawLatex(DEF_InfoBoxX + DEF_LastColumnInfoXOffset, DEF_InfoBoxY, "log#left(#frac{1}{N_{high}} #frac{d^{2}N}{d#eta dp_{t}}#right)");
            latex_eta.DrawLatex(DEF_InfoBoxX + DEF_LastColumnInfoXOffset, DEF_InfoBoxY - DEF_InfoBoxLineDist, "Leading p_{t} #in [80, #infty) #frac{GeV}{c}");
        }

        // -------- PHI ROW --------
        can_Distributions->cd(i+3);
        gPad->SetMargin(0.16*(i==1), 0.15*(i==3), 0.15, 0.02);

        h2D_ptvsphi[i-1]->SetTitle("");
        h2D_ptvsphi[i-1]->GetYaxis()->SetRangeUser(0, 60);
        h2D_ptvsphi[i-1]->GetXaxis()->SetRangeUser(-0.5, 0.5);

        h2D_ptvsphi[i-1]->Scale(1./h2D_ptvsphi[i-1]->GetEntries());

        if(i==1){
            h2D_ptvsphi[i-1]->GetYaxis()->SetTitle("p_{t} [#frac{GeV}{c}]");
            h2D_ptvsphi[i-1]->GetYaxis()->SetTitleSize(DEF_AxisLabelSize);
        } else {
            h2D_ptvsphi[i-1]->GetYaxis()->SetLabelSize(0);
            h2D_ptvsphi[i-1]->GetYaxis()->SetTitle("");
        }

        h2D_ptvsphi[i-1]->GetYaxis()->SetLabelSize(DEF_AxisLabelSize);

        if(i<3){
            //h2D_ptvsphi[i-1]->GetXaxis()->SetLabelSize(0);
            h2D_ptvsphi[i-1]->GetXaxis()->SetTitle("");
        } else {
            h2D_ptvsphi[i-1]->GetXaxis()->SetTitle("#Delta #phi [rad]");
            h2D_ptvsphi[i-1]->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
        }

        h2D_ptvsphi[i-1]->GetXaxis()->SetLabelSize(DEF_AxisLabelSize);

        //h2D_ptvsphi[i-1]->SetMaximum(0.013);

        if(i==3)h2D_ptvsphi[i-1]->Draw("COLZ");
        else h2D_ptvsphi[i-1]->Draw();

        TLatex latex_phi;
        latex_phi.SetNDC();
        latex_phi.SetTextSize(DEF_LatexTextSize);
        latex_phi.SetTextAlign(31);

        if(i==1){
            latex_phi.DrawLatex(DEF_InfoBoxX, DEF_InfoBoxY, "log#left(#frac{1}{N_{low}} #frac{d^{2}N}{d#phi dp_{t}}#right)");
            latex_phi.DrawLatex(DEF_InfoBoxX, DEF_InfoBoxY - DEF_InfoBoxLineDist, "Leading p_{t} #in [30, 50) #frac{GeV}{c}");
        } else if(i==2){
            latex_phi.DrawLatex(DEF_InfoBoxX, DEF_InfoBoxY, "log#left(#frac{1}{N_{middle}} #frac{d^{2}N}{d#phi dp_{t}}#right)");
            latex_phi.DrawLatex(DEF_InfoBoxX, DEF_InfoBoxY - DEF_InfoBoxLineDist, "Leading p_{t} #in [50, 80) #frac{GeV}{c}");
        } else {
            latex_phi.DrawLatex(DEF_InfoBoxX + DEF_LastColumnInfoXOffset, DEF_InfoBoxY, "log#left(#frac{1}{N_{high}} #frac{d^{2}N}{d#phi dp_{t}}#right)");
            latex_phi.DrawLatex(DEF_InfoBoxX + DEF_LastColumnInfoXOffset, DEF_InfoBoxY - DEF_InfoBoxLineDist, "Leading p_{t} #in [80, #infty) #frac{GeV}{c}");
        }
    }
}

void Analysis_JetMultiplicityDistribution(const TString DataFile)
{
    #define DEF_AxisLabelSize 0.05
    #define DEF_HistoTitleSize 0.1
    #define DEF_InfoBoxX 0.99
    #define DEF_InfoBoxY 0.9
    #define DEF_LastColumnInfoXOffset -0.15
    #define DEF_InfoBoxLineDist 0.1
    #define DEF_LatexTextSize 0.04

    TFile *file = TFile::Open(DataFile);
    if (!file || file->IsZombie()) {
        std::cout << "Error while opening the file!" << std::endl;
    }

    gStyle->SetOptStat(0);
    SetRootGraphicStyle();

    TH2D* h2D_ptvseta[3];
    TH2D* h2D_ptvsphi[3];

    h2D_ptvseta[0] = (TH2D*)file->Get("Pt_vs_Eta/h2D_pt_vs_eta_Summed_Lowpt");
    if(!h2D_ptvseta[0]) cout << "FILE NOT FOUND" << endl;
    h2D_ptvseta[1] = (TH2D*)file->Get("Pt_vs_Eta/h2D_pt_vs_eta_Summed_Middlept");
    if(!h2D_ptvseta[1]) cout << "FILE NOT FOUND" << endl;
    h2D_ptvseta[2] = (TH2D*)file->Get("Pt_vs_Eta/h2D_pt_vs_eta_Summed_Highpt");
    if(!h2D_ptvseta[2]) cout << "FILE NOT FOUND" << endl;

    h2D_ptvsphi[0] = (TH2D*)file->Get("Pt_vs_Phi/h2D_pt_vs_phi_Summed_Lowpt");
    if(!h2D_ptvsphi[0]) cout << "FILE NOT FOUND" << endl;
    h2D_ptvsphi[1] = (TH2D*)file->Get("Pt_vs_Phi/h2D_pt_vs_phi_Summed_Middlept");
    if(!h2D_ptvsphi[1]) cout << "FILE NOT FOUND" << endl;
    h2D_ptvsphi[2] = (TH2D*)file->Get("Pt_vs_Phi/h2D_pt_vs_phi_Summed_Highpt");
    if(!h2D_ptvsphi[2]) cout << "FILE NOT FOUND" << endl;

    TCanvas* can_Distributions = new TCanvas("Jet and Recoil Site Multiplicities","Jet and Recoil Site Multiplicities",1000,1000);
    can_Distributions->Divide(3,2,0,0);

    for(Int_t i=1; i<=3; i++)
    {
        // -------- ETA ROW --------
        can_Distributions->cd(i);
        gPad->SetMargin(0.16*(i==1), 0.15*(i==3), 0.15, 0.02);

        h2D_ptvseta[i-1]->SetTitle("");
        h2D_ptvseta[i-1]->GetYaxis()->SetRangeUser(0, 60);
        h2D_ptvseta[i-1]->GetXaxis()->SetRangeUser(-0.5, 0.5);

        h2D_ptvseta[i-1]->Scale(1./h2D_ptvseta[i-1]->GetEntries());

        // Axis handling
        if(i==1){
            h2D_ptvseta[i-1]->GetYaxis()->SetTitle("p_{t} [#frac{GeV}{c}]");
            h2D_ptvseta[i-1]->GetYaxis()->SetTitleSize(DEF_AxisLabelSize);
        } else {
            h2D_ptvseta[i-1]->GetYaxis()->SetLabelSize(0);
            h2D_ptvseta[i-1]->GetYaxis()->SetTitle("");
        }
        h2D_ptvseta[i-1]->GetYaxis()->SetLabelSize(DEF_AxisLabelSize);

        if(i<3){
            //h2D_ptvseta[i-1]->GetXaxis()->SetLabelSize(0);
            h2D_ptvseta[i-1]->GetXaxis()->SetTitle("");
        } else {
            h2D_ptvseta[i-1]->GetXaxis()->SetTitle("#Delta #eta");
            h2D_ptvseta[i-1]->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
        }
        h2D_ptvseta[i-1]->GetXaxis()->SetLabelSize(DEF_AxisLabelSize);

        if(i==3)h2D_ptvseta[i-1]->Draw("COLZ");
        else h2D_ptvseta[i-1]->Draw();

        TLatex latex_eta;
        latex_eta.SetNDC();
        latex_eta.SetTextSize(DEF_LatexTextSize);
        latex_eta.SetTextAlign(31);

        if(i==1){
            latex_eta.DrawLatex(DEF_InfoBoxX, DEF_InfoBoxY, "log#left(#frac{1}{N_{low}} #frac{d^{2}N}{d#eta dp_{t}}#right)");
            latex_eta.DrawLatex(DEF_InfoBoxX, DEF_InfoBoxY - DEF_InfoBoxLineDist, "Leading p_{t} #in [30, 50) #frac{GeV}{c}");
        } else if(i==2){
            latex_eta.DrawLatex(DEF_InfoBoxX, DEF_InfoBoxY, "log#left(#frac{1}{N_{middle}} #frac{d^{2}N}{d#eta dp_{t}}#right)");
            latex_eta.DrawLatex(DEF_InfoBoxX, DEF_InfoBoxY - DEF_InfoBoxLineDist, "Leading p_{t} #in [50, 80) #frac{GeV}{c}");
        } else {
            latex_eta.DrawLatex(DEF_InfoBoxX + DEF_LastColumnInfoXOffset, DEF_InfoBoxY, "log#left(#frac{1}{N_{high}} #frac{d^{2}N}{d#eta dp_{t}}#right)");
            latex_eta.DrawLatex(DEF_InfoBoxX + DEF_LastColumnInfoXOffset, DEF_InfoBoxY - DEF_InfoBoxLineDist, "Leading p_{t} #in [80, #infty) #frac{GeV}{c}");
        }

        // -------- PHI ROW --------
        can_Distributions->cd(i+3);
        gPad->SetMargin(0.16*(i==1), 0.15*(i==3), 0.15, 0.02);

        h2D_ptvsphi[i-1]->SetTitle("");
        h2D_ptvsphi[i-1]->GetYaxis()->SetRangeUser(0, 60);
        h2D_ptvsphi[i-1]->GetXaxis()->SetRangeUser(-0.5, 0.5);

        h2D_ptvsphi[i-1]->Scale(1./h2D_ptvsphi[i-1]->GetEntries());

        if(i==1){
            h2D_ptvsphi[i-1]->GetYaxis()->SetTitle("p_{t} [#frac{GeV}{c}]");
            h2D_ptvsphi[i-1]->GetYaxis()->SetTitleSize(DEF_AxisLabelSize);
        } else {
            h2D_ptvsphi[i-1]->GetYaxis()->SetLabelSize(0);
            h2D_ptvsphi[i-1]->GetYaxis()->SetTitle("");
        }
        h2D_ptvsphi[i-1]->GetYaxis()->SetLabelSize(DEF_AxisLabelSize);

        if(i<3){
            //h2D_ptvsphi[i-1]->GetXaxis()->SetLabelSize(0);
            h2D_ptvsphi[i-1]->GetXaxis()->SetTitle("");
        } else {
            h2D_ptvsphi[i-1]->GetXaxis()->SetTitle("#Delta #phi [rad]");
            h2D_ptvsphi[i-1]->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
        }
        h2D_ptvsphi[i-1]->GetXaxis()->SetLabelSize(DEF_AxisLabelSize);

        if(i==3)h2D_ptvsphi[i-1]->Draw("COLZ");
        else h2D_ptvsphi[i-1]->Draw();

        TLatex latex_phi;
        latex_phi.SetNDC();
        latex_phi.SetTextSize(DEF_LatexTextSize);
        latex_phi.SetTextAlign(31);

        if(i==1){
            latex_phi.DrawLatex(DEF_InfoBoxX, DEF_InfoBoxY, "log#left(#frac{1}{N_{low}} #frac{d^{2}N}{d#phi dp_{t}}#right)");
            latex_phi.DrawLatex(DEF_InfoBoxX, DEF_InfoBoxY - DEF_InfoBoxLineDist, "Leading p_{t} #in [30, 50) #frac{GeV}{c}");
        } else if(i==2){
            latex_phi.DrawLatex(DEF_InfoBoxX, DEF_InfoBoxY, "log#left(#frac{1}{N_{middle}} #frac{d^{2}N}{d#phi dp_{t}}#right)");
            latex_phi.DrawLatex(DEF_InfoBoxX, DEF_InfoBoxY - DEF_InfoBoxLineDist, "Leading p_{t} #in [50, 80) #frac{GeV}{c}");
        } else {
            latex_phi.DrawLatex(DEF_InfoBoxX + DEF_LastColumnInfoXOffset, DEF_InfoBoxY, "log#left(#frac{1}{N_{high}} #frac{d^{2}N}{d#phi dp_{t}}#right)");
            latex_phi.DrawLatex(DEF_InfoBoxX + DEF_LastColumnInfoXOffset, DEF_InfoBoxY - DEF_InfoBoxLineDist, "Leading p_{t} #in [80, #infty) #frac{GeV}{c}");
        }
    }
}


void Analysis_EtaVsPhi_Dijet(const TString DataFile) {


    //open the root file
    TFile *file = TFile::Open(DataFile);

    if (!file || file->IsZombie()) {
        std::cout << "Error while opening the file!" << std::endl;
    }

    //load histograms from root file
    TH2D* h2D_eta_vs_phi = (TH2D*)file->Get("Dijets/h2D_phi_vs_eta;1390");
    if(!h2D_eta_vs_phi) cout << "Background histogram not found!" << endl;

    //configure and draw canvas
    TCanvas* can = new TCanvas("can","can",1000,1000);

    gStyle->SetOptStat(0);

    h2D_eta_vs_phi->SetTitle("");
    h2D_eta_vs_phi->SetMaximum(50);

    h2D_eta_vs_phi->GetXaxis()->SetLabelSize(0.03);
    h2D_eta_vs_phi->GetXaxis()->SetNdivisions(8);
    h2D_eta_vs_phi->GetXaxis()->SetTitleSize(0.03);
    h2D_eta_vs_phi->GetXaxis()->SetTitle("#eta");
    h2D_eta_vs_phi->GetXaxis()->SetRangeUser(-0.9, 0.9);

    h2D_eta_vs_phi->GetYaxis()->SetLabelSize(0.03);
    h2D_eta_vs_phi->GetYaxis()->SetNdivisions(8);
    h2D_eta_vs_phi->GetYaxis()->SetTitleSize(0.03);
    h2D_eta_vs_phi->GetYaxis()->SetTitle("#phi [rad]");
    h2D_eta_vs_phi->GetYaxis()->SetRangeUser(-Pi, Pi);

    h2D_eta_vs_phi->GetZaxis()->SetLabelSize(0.03);
    h2D_eta_vs_phi->GetZaxis()->SetNdivisions(5);
    h2D_eta_vs_phi->GetZaxis()->SetTitleSize(0.03);
    h2D_eta_vs_phi->GetZaxis()->SetTitle("p_{t} [#frac{GeV}{c}]");

    h2D_eta_vs_phi->Draw("lego2, fb");

}
