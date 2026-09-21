#include "./ReadTree.cc"

using namespace fastjet;
using namespace std;

static const int arr_color_jet[40] = {kOrange,kGreen+2,kMagenta,kBlue,kCyan,kOrange,kGreen+2,kMagenta,kBlue,kCyan,kOrange,kGreen+2,kMagenta,kBlue,kCyan,kOrange,kGreen+2,kMagenta,kBlue,kCyan,kOrange,kGreen+2,kMagenta,kBlue,kCyan,kOrange,kGreen+2,kMagenta,kBlue,kCyan,kOrange,kGreen+2,kMagenta,kBlue,kCyan,kOrange,kGreen+2,kMagenta,kBlue,kCyan};

#include "fastjet/tools/JetMedianBackgroundEstimator.hh"

void Plot_DiJets(int event = 1313, int rebinX = 2, int rebinY = 2, float theta = 45.0, float phi = 70.0)
{
    printf("Plot_DiJets started for event: %d \n",event);

    double min_jet_pt_sub_draw_label = 40.0;

    SetRootGraphicStyle();
    TFile* inputfile = TFile::Open("./ResultRoots/OF_HighPtRecoilDistr_BGCompensated_wrong.root");
    HistName = "/Dijets/h2D_phi_vs_eta;";
    HistName += event;
    TH2D* h2D_dijet = (TH2D*)inputfile->Get(HistName.Data());
    h2D_dijet ->Rebin2D(rebinX,rebinY);



    //------------------------------------------------
    // jet finder
    double eJet_R = 0.3;
    double eBkg_R = 0.3;
    int eRem_n_hardest = 2;
    double min_jet_pt_select = 0.0;
    vector<PseudoJet> vec_PJ_particles;
    vector<PseudoJet> jets_fiducial[2]; // [not smeared, smeared] only for PYTHIA
    TLorentzVector TLV_part;
    printf("Fill pseudo jet \n");
    int particle_counter = 0;
    for(int ibinx = 1; ibinx <= h2D_dijet->GetNbinsX(); ibinx++)
    {
        for(int ibiny = 1; ibiny <= h2D_dijet->GetNbinsY(); ibiny++)
        {
            double eta = h2D_dijet ->GetXaxis()->GetBinCenter(ibinx);
            double phi = h2D_dijet ->GetYaxis()->GetBinCenter(ibiny);
            double pt  = h2D_dijet ->GetBinContent(ibinx,ibiny);

            if(pt < 0.15) continue;


            TLV_part.SetPtEtaPhiM(pt,eta,phi,0.139);
            PseudoJet Fill_PseudoJet(TLV_part.Px(),TLV_part.Py(),TLV_part.Pz(),TLV_part.E());
            vec_PJ_particles.push_back(Fill_PseudoJet);
            //printf("particle: %d, eta: %4.3f, phi: %4.3f, pt: %4.3f \n",particle_counter,eta,phi,pt);
            particle_counter++;
        }
    }




    //--------------------
    printf("Pseudo jet filled \n");
    JetDefinition jet_def(antikt_algorithm, eJet_R);
    // jet area definition
    printf("Define ghost maxrap \n");
    Double_t ghost_maxrap = 1.0; // Fiducial cut for background estimation
    GhostedAreaSpec area_spec(ghost_maxrap);
    printf("Define  area \n");
    AreaDefinition area_def(active_area_explicit_ghosts,GhostedAreaSpec(ghost_maxrap,1,0.01));
    printf("Define ClusterSequenceArea \n");
    ClusterSequenceArea clust_seq_hard(vec_PJ_particles, jet_def, area_def);

    // run the clustering, extract the jets
    printf("Start of clustering \n");
    double ptmin = 0.15;
    vector<PseudoJet> jets_all = sorted_by_pt(clust_seq_hard.inclusive_jets(ptmin));
    //Selector Fiducial_cut_selector = SelectorAbsEtaMax(0.9 - eJet_R); // Fiducial cut for jets    10.0 -> 1.0
    Selector Fiducial_cut_selector = SelectorAbsEtaMax(5.0); // Fiducial cut for jets    10.0 -> 1.0
    jets_fiducial[0] = Fiducial_cut_selector(jets_all);
    //--------------------




    //--------------------
    // background estimation
    Double_t jet_rho_array[2];
    //cout << "Define JetDefinition" << endl;
    JetDefinition jet_def_bkgd(kt_algorithm, eBkg_R); // <--
    //JetDefinition jet_def_bkgd(antikt_algorithm, jet_R); // test
    //cout << "Define AreaDefinition" << endl;
    AreaDefinition area_def_bkgd(active_area_explicit_ghosts,GhostedAreaSpec(ghost_maxrap,1,0.01));
    //AreaDefinition area_def_bkgd(active_area,GhostedAreaSpec(ghost_maxrap,1,0.005));
    //cout << "Define selector" << endl;
    Int_t Rem_n_hardest_use = eRem_n_hardest;
    Selector selector = SelectorAbsEtaMax(1.0) * (!SelectorNHardest(Rem_n_hardest_use)); // <--
    //Selector selector = SelectorAbsEtaMax(1.0 - jet_R); // test

    //cout << "Define JetMedianBackgroundEstimator" << endl;
    JetMedianBackgroundEstimator bkgd_estimator(selector, jet_def_bkgd, area_def_bkgd); // <--
    //JetMedianBackgroundEstimator bkgd_estimator(selector, jet_def, area_def); // test
    //cout << "Define Subtractor" << endl;

    Subtractor subtractor(&bkgd_estimator);
    //cout << "Define bkgd_estimator" << endl;
    bkgd_estimator.set_particles(vec_PJ_particles);

    //cout << "Calculate jet_rho and jet_sigma" << endl;
    //Double_t jet_rho   = bkgd_estimator.rho() - 0.212;   // TO BE CHANGED, equivallent to a 60 MeV shift for R = 0.3 (rho*pi*R*2 = 60 MeV)
    Double_t jet_rho   = bkgd_estimator.rho();
    jet_rho_array[0] = jet_rho;
    //cout << "jet_sigma" << endl;
    Double_t jet_sigma = bkgd_estimator.sigma();
    //printf("rho: %4.3f, sigma: %4.3f \n",jet_rho,jet_sigma);
    //--------------------



    vector<TH2D*> vec_h2D_dijet_select;
    vector<double> vec_jet_pt_sub;
    vector<double> vec_jet_pt;

    for(Int_t i_jet = 0; i_jet < (Int_t)jets_fiducial[0].size(); i_jet++)
    {
        Float_t jet_pt     = jets_fiducial[0][i_jet].perp();
        if(jet_pt <= 0.0) continue;
        Float_t jet_area   = jets_fiducial[0][i_jet].area();
        Float_t jet_pt_sub = jets_fiducial[0][i_jet].perp() - jet_rho*jet_area;
        Float_t jet_eta    = jets_fiducial[0][i_jet].eta();
        Float_t jet_phi    = jets_fiducial[0][i_jet].phi();

        //if(jet_pt_sub > min_jet_pt_select)
        {
            HistName = "vec_h2D_dijet_select_";
            HistName += i_jet;
            vec_h2D_dijet_select.push_back((TH2D*)h2D_dijet->Clone(HistName.Data()));
            vec_h2D_dijet_select[vec_h2D_dijet_select.size()-1] ->Reset();
            vec_jet_pt_sub.push_back(jet_pt_sub);
            vec_jet_pt.push_back(jet_pt);

            vector<PseudoJet> jet_constituents = jets_fiducial[0][i_jet].constituents();
            int N_constituents = (Int_t)jet_constituents.size();
            printf("Selected jet: %d, jet_pt: %4.3f, jet_pt-rho*area: %4.3f, N_constituents: %d \n ",i_jet,jet_pt,jet_pt_sub,N_constituents);

            for(Int_t i_constituent = 0; i_constituent < N_constituents; i_constituent++)
            {
                Float_t jet_const_pt  = jet_constituents[i_constituent].perp();
                if(jet_const_pt < 0.15) continue;
                Float_t jet_const_phi = jet_constituents[i_constituent].phi();
                if(jet_const_phi > TMath::Pi()) jet_const_phi -= 2.0*TMath::Pi();
                Float_t jet_const_eta = jet_constituents[i_constituent].eta();
                Int_t   user_index    = jet_constituents[i_constituent].user_index();

                int eta_bin = vec_h2D_dijet_select[vec_h2D_dijet_select.size()-1] ->GetXaxis()->FindBin(jet_const_eta);
                int phi_bin = vec_h2D_dijet_select[vec_h2D_dijet_select.size()-1] ->GetYaxis()->FindBin(jet_const_phi);
                vec_h2D_dijet_select[vec_h2D_dijet_select.size()-1] ->SetBinContent(eta_bin,phi_bin,jet_const_pt);
                h2D_dijet                                           ->SetBinContent(eta_bin,phi_bin,0);

                //printf("i_constituent: %d, eta: %4.3f, phi: %4.3f, pt: %4.3f \n",i_constituent,jet_const_eta,jet_const_phi,jet_const_pt);
            }
        }
    }
    //------------------------------------------------




    //------------------------------------------------
    // Create a THStack for 2D histograms
    THStack *hs = new THStack("hs", "Superimposed LEGO Plots");

    h2D_dijet ->GetXaxis()->SetTitle("#eta");
    h2D_dijet ->GetYaxis()->SetTitle("#phi (rad)");
    h2D_dijet ->GetZaxis()->SetTitle("#it{p}^{track}_{T} (GeV/#it{c})");
    h2D_dijet ->GetZaxis()->SetLabelSize(0.055);
    h2D_dijet ->GetZaxis()->SetTitleSize(0.055);
    h2D_dijet ->GetZaxis()->SetTitleOffset(0.7);
    h2D_dijet ->GetZaxis()->SetNdivisions(505,'N');
    //h2D_dijet ->GetZaxis()->SetRangeUser(0.0,15.0);
    //TCanvas* can_h2D_dijet = Draw_2D_histo_and_canvas(h2D_dijet,"can_h2D_dijet",1200,500,0.15,60.0,"LEGO2 FBBB"); // TH2D* hist, TString name, Int_t x_size, Int_t y_size,Double_t min_val, Double_t max_val, TString option
    //TCanvas* can_h2D_dijet = Draw_2D_histo_and_canvas(h2D_dijet,"can_h2D_dijet",1200,500,0.15,60.0,"LEGO 0 FBBB"); // TH2D* hist, TString name, Int_t x_size, Int_t y_size,Double_t min_val, Double_t max_val, TString option
    //can_h2D_dijet ->cd()->SetLogz(0);
    //can_h2D_dijet ->SetTheta(theta); // 45
    //can_h2D_dijet ->SetPhi(phi); // 60

    //can_h2D_dijet ->cd();

    h2D_dijet->SetFillColor(kGray);
    hs->Add(h2D_dijet);

    for(int i_jet = 0; i_jet < vec_h2D_dijet_select.size(); i_jet++)
    {
        double jet_pt_sub = vec_jet_pt_sub[i_jet];
        double jet_pt     = vec_jet_pt[i_jet];

        vec_h2D_dijet_select[i_jet] ->SetFillColor(arr_color_jet[i_jet]);
        if(jet_pt_sub >= min_jet_pt_sub_draw_label) vec_h2D_dijet_select[i_jet] ->SetFillColor(kRed);
        //vec_h2D_dijet_select[i_jet] ->Draw("same LEGO1 0");
        hs->Add(vec_h2D_dijet_select[i_jet]);
    }

    // Draw the stack in lego style
    TCanvas* can_h2D_dijet = new TCanvas("can_h2D_dijet","can_h2D_dijet",10,10,1200,500);
    can_h2D_dijet->SetFillColor(10);
    can_h2D_dijet->SetTopMargin(0.05);
    can_h2D_dijet->SetBottomMargin(0.2);
    can_h2D_dijet->SetRightMargin(0.22);
    can_h2D_dijet->SetLeftMargin(0.2);
    can_h2D_dijet->SetTicks(1,1);
    can_h2D_dijet->SetGrid(0,0);
    can_h2D_dijet->cd()->SetLogz(0);
    can_h2D_dijet->SetTheta(theta); // 45
    can_h2D_dijet->SetPhi(phi); // 60

    printf("Draw stack \n");

    //hs->SetStats(0);

    hs->SetTitle("");


    hs->Draw("lego1 0 FBBB");
    hs->GetXaxis()->SetNdivisions(505,'N');
    hs->GetYaxis()->SetNdivisions(505,'N');
    hs->GetXaxis()->CenterTitle();
    hs->GetYaxis()->CenterTitle();
    hs->GetXaxis()->SetTitle("#eta");
    hs->GetYaxis()->SetTitle("#varphi (rad)");
    hs->GetZaxis()->SetTitle("#it{p}^{track}_{T} (GeV/#it{c})");
    hs->GetZaxis()->SetLabelSize(0.055);
    hs->GetZaxis()->SetTitleSize(0.055);
    hs->GetZaxis()->SetTitleOffset(0.7);
    hs->GetZaxis()->SetNdivisions(505,'N');


    can_h2D_dijet ->Update();


    // 1. Ensure the canvas has been drawn and updated
    // This forces ROOT to generate the 3D geometry and projection matrix

    for(Int_t i_jet = 0; i_jet < (Int_t)vec_h2D_dijet_select.size(); i_jet++)
    {
        double jet_pt_sub = vec_jet_pt_sub[i_jet];
        double jet_pt     = vec_jet_pt[i_jet];

        if(jet_pt_sub < min_jet_pt_sub_draw_label) continue;

        double max_pt = 0.0;
        int ibinx_max, ibiny_max;
        for(int ibinx = 1; ibinx <= vec_h2D_dijet_select[i_jet]->GetNbinsX(); ibinx++)
        {
            for(int ibiny = 1; ibiny <= vec_h2D_dijet_select[i_jet]->GetNbinsY(); ibiny++)
            {
                double eta =  vec_h2D_dijet_select[i_jet]->GetXaxis()->GetBinCenter(ibinx);
                double phi =  vec_h2D_dijet_select[i_jet]->GetYaxis()->GetBinCenter(ibiny);
                double pt  =  vec_h2D_dijet_select[i_jet]->GetBinContent(ibinx,ibiny);

                if(pt > max_pt)
                {
                    max_pt = pt;
                    ibinx_max = ibinx;
                    ibiny_max = ibiny;
                }
            }
        }


        // 2. Get the physical center coordinates (X, Y) of the bin
        double x_val = vec_h2D_dijet_select[i_jet]->GetXaxis()->GetBinCenter(ibinx_max);
        double y_val = vec_h2D_dijet_select[i_jet]->GetYaxis()->GetBinCenter(ibiny_max);

        // For a LEGO plot, the Z value determines the height in 3D space
        double z_val = vec_h2D_dijet_select[i_jet]->GetBinContent(ibinx_max, ibiny_max);

        printf("pos: {%4.3f, %4.3f, %4.3f} \n",x_val,y_val,z_val);

        // 3. Create an array holding the 3D world coordinates
        double world_coords[3] = {x_val, y_val, z_val};
        int pixel_coords[2] = {0, 0};

        // 4. Extract the 3D view object from the current pad
        TView *view = can_h2D_dijet->GetView();

        if (view) {
            // Convert 3D world coordinates to 2D absolute pixel coordinates
            view->WCtoNDC(world_coords, world_coords); // Optional step depending on ROOT version, converts to Normalized Device Coords

            // Direct pixel conversion helper method
            pixel_coords[0] = can_h2D_dijet->XtoAbsPixel(world_coords[0]);
            pixel_coords[1] = can_h2D_dijet->YtoAbsPixel(world_coords[1]);

            std::cout << "Pixel X: " << pixel_coords[0] << ", Pixel Y: " << pixel_coords[1] << std::endl;
        } else {
            std::cout << "Error: Pad view not initialized. Make sure to draw with \"LEGO\" and call can_h2D_dijet->Update() first." << std::endl;
        }

        // 2. Convert pixels to user X and Y coordinates
        Double_t userX, userY;
        can_h2D_dijet->AbsPixeltoXY(pixel_coords[0]-15, pixel_coords[1]-15, userX, userY);

        //jet_pt_sub
        sprintf(NoP,"%4.1f",jet_pt);
        HistName = NoP;
        HistName += " GeV/#it{c}";
        plotTopLegend((char*)HistName.Data(),userX,userY,0.04,1,0.0,42,0,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1

        sprintf(NoP,"%4.1f",jet_pt_sub);
        HistName = NoP;
        HistName += " GeV/#it{c}";
        //plotTopLegend((char*)HistName.Data(),userX,userY,0.04,1,0.0,42,0,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
    }
    //------------------------------------------------

}

Int_t CompareDijetResults()
{

    #define DEF_AxisLabelSize 0.07
    #define DEF_AxisTitleSize 0.08
    #define DEF_HistoTitleSize 0.1
    #define DEF_LegendFontSize 0.055
    #define DEF_Margin_Top 0.05
    #define DEF_Margin_Bottom 0.23
    #define DEF_Margin_Left 0.18
    #define DEF_Margin_Left_ForMiddlePanel 0.1
    #define DEF_Margin_Right_ForMiddlePanel 0.01

    float LowPtCuts[4] = {0.0f, 1.0, 2.0f, 4.0f};
    float HighPtCuts[4] = {1.0f, 2.0f, 4.0f, 6.0f};

    double xbins[11] = {-0.9, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 0.9};
    Int_t nbins = 10;

    gStyle->SetOptStat(0);
    SetRootGraphicStyle();

    //open the root file
    TFile *file = TFile::Open("ResultRoots/DijetAna/DijetAna_Id_Jet_Pt_40_20.root");

    if (!file || file->IsZombie()) {
        std::cout << "Error while opening the file!" << std::endl;
    }

    //open the root file
    TFile *File_Yong = TFile::Open("ResultRoots/DijetAna_Yongzhen/OutputHadronEta_dj40_sub20.root");

    if (!File_Yong || File_Yong->IsZombie()) {
        std::cout << "Error while opening the file!" << std::endl;
        return 0;
    }

    TFile *File_Hydro = TFile::Open("ResultRoots/DijetAna_HydroSimulation/OutputHadronEta_ld40_sub20.root");

    if (!File_Hydro || File_Hydro->IsZombie()) {
        std::cout << "Error while opening the file!" << std::endl;
        return 0;
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
    TH1F* h1D_Hydro[4];
    h1D_Hydro[1] = (TH1F*)File_Hydro->Get("trkpt_100_200/hetadiffS_100_200");
    h1D_Hydro[2] = (TH1F*)File_Hydro->Get("trkpt_200_400/hetadiffS_200_400");


    TH1D* h1F_Yongzhen[4];
    h1F_Yongzhen[0] = (TH1D*)File_Yong->Get("trkpt_015_100/hetadiffS_015_100");
    h1F_Yongzhen[1] = (TH1D*)File_Yong->Get("trkpt_100_200/hetadiffS_100_200");
    h1F_Yongzhen[2] = (TH1D*)File_Yong->Get("trkpt_200_400/hetadiffS_200_400");
    h1F_Yongzhen[3] = (TH1D*)File_Yong->Get("trkpt_400_600/hetadiffS_400_600");
    if(!(h1F_Yongzhen[0] && h1F_Yongzhen[1] && h1F_Yongzhen[2] && h1F_Yongzhen[3])){ cout << "histograms not found!" << endl; return 0;}

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

    //get bin width before
    float OriginalBinWidth = h1F_SmallGap[0]->GetXaxis()->GetBinWidth(1);

    //configure and draw canvas
    TCanvas* can_ParticleMultiplicities = new TCanvas("ParticleMultiplicities","ParticleMultiplicities",1000,1000);
    can_ParticleMultiplicities->Divide(1,4, 0, 0);

    for(int iRow=0; iRow<=3; iRow++)
    {

        can_ParticleMultiplicities->cd(1 + iRow);

        gPad->SetMargin(DEF_Margin_Left, DEF_Margin_Right_ForMiddlePanel, DEF_Margin_Bottom*(iRow==3), DEF_Margin_Top*(iRow==0));

        h1F_SmallGap[iRow]->SetTitle("");
        h1F_SmallGap[iRow]->SetTitleSize(DEF_HistoTitleSize);
        h1F_SmallGap[iRow]->GetXaxis()->SetTitleSize(DEF_AxisTitleSize);
        h1F_SmallGap[iRow]->GetYaxis()->SetTitleSize(DEF_AxisTitleSize);
        h1F_SmallGap[iRow]->GetXaxis()->SetLabelSize(DEF_AxisLabelSize);
        h1F_SmallGap[iRow]->GetYaxis()->SetLabelSize(DEF_AxisLabelSize);
        (iRow==0) ? (h1F_SmallGap[iRow]->GetYaxis()->SetTitle("C_{Large/Small Gap}")) : (h1F_SmallGap[iRow]->GetYaxis()->SetTitle(""));
        h1F_SmallGap[iRow] = (TH1F*) h1F_SmallGap[iRow]->Rebin(nbins, Form("h1F_SmallGap_%d_rebinned", iRow), xbins);
        for (Int_t i = 1; i <= h1F_SmallGap[iRow]->GetNbinsX(); i++) {
            float width   = h1F_SmallGap[iRow]->GetBinWidth(i);
            float content = h1F_SmallGap[iRow]->GetBinContent(i);
            float error   = h1F_SmallGap[iRow]->GetBinError(i);

            h1F_SmallGap[iRow]->SetBinContent(i, content / (width/OriginalBinWidth));
            h1F_SmallGap[iRow]->SetBinError(i, error / (width/OriginalBinWidth));
        }
        h1F_SmallGap[iRow]->SetMarkerStyle(kCircle);
        h1F_SmallGap[iRow]->SetMarkerColor(kBlue);

        h1F_LargeGap[iRow]->SetTitle("");
        h1F_LargeGap[iRow]->SetTitleSize(DEF_HistoTitleSize);
        h1F_LargeGap[iRow]->GetXaxis()->SetTitleSize(DEF_AxisLabelSize);
        h1F_LargeGap[iRow]->GetYaxis()->SetTitleSize(DEF_AxisLabelSize);
        h1F_LargeGap[iRow]->GetXaxis()->SetLabelSize(DEF_AxisLabelSize);
        h1F_LargeGap[iRow]->GetYaxis()->SetLabelSize(DEF_AxisLabelSize);
        h1F_LargeGap[iRow] = (TH1F*) h1F_LargeGap[iRow]->Rebin(nbins, Form("h1F_LargeGap_%d_rebinned", iRow), xbins);
        for (Int_t i = 1; i <= h1F_LargeGap[iRow]->GetNbinsX(); i++) {
            float width   = h1F_LargeGap[iRow]->GetBinWidth(i);
            cout << "width " << i << " = " << width << endl;
            float content = h1F_LargeGap[iRow]->GetBinContent(i);
            float error   = h1F_LargeGap[iRow]->GetBinError(i);

            cout << content << "," << width << "," << content / width << endl;

            h1F_LargeGap[iRow]->SetBinContent(i, content / (width/OriginalBinWidth));
            h1F_LargeGap[iRow]->SetBinError(i, error / (width/OriginalBinWidth));
        }
        h1F_LargeGap[iRow]->SetMarkerStyle(kCircle);
        h1F_LargeGap[iRow]->SetMarkerColor(kRed);
        
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
        h1F_LargeGap[iRow]->SetLineColor(kBlack);
        h1F_LargeGap[iRow]->DrawCopy("");

        //ALICE 0-30% reference
        h1F_Yongzhen[iRow]->SetMarkerStyle(kCircle);
        h1F_Yongzhen[iRow]->SetMarkerColor(kBlue);
        h1F_Yongzhen[iRow]->SetLineColor(kBlue);
        h1F_Yongzhen[iRow]->DrawCopy("SAME");


        //hydro prediction
        if(iRow==1 || iRow==2)
        {
            TGraphAsymmErrors *band = new TGraphAsymmErrors(h1D_Hydro[iRow]);

            for (int i = 0; i < h1D_Hydro[iRow]->GetNbinsX(); ++i) {
                double x  = h1D_Hydro[iRow]->GetBinCenter(i + 1);
                double ex = h1D_Hydro[iRow]->GetBinWidth(i + 1) / 2.0;

                double y  = h1D_Hydro[iRow]->GetBinContent(i + 1);
                double ey = h1D_Hydro[iRow]->GetBinError(i + 1);

                band->SetPoint(i, x, y);
                band->SetPointError(i, ex, ex, ey, ey);
            }

            band->SetFillColorAlpha(kGray, 0.5);
            band->SetLineColor(kGray);
            band->SetLineWidth(2);
            
            band->Draw("3, SAME");

        }

        gPad->RedrawAxis();


        TLine *line = new TLine(-.9, 0.0, .9, 0.0);
        line->SetLineColor(kBlack);
        line->Draw("SAME");

        TH1F* Pointer1 = (TH1F*) h1F_LargeGap[iRow]->DrawCopy("same");
        TH1F* Pointer2 = (TH1F*) h1F_Yongzhen[iRow]->DrawCopy("SAME");

        TGraphAsymmErrors *Pointer3 = new TGraphAsymmErrors(h1F_LargeGap[iRow]);
        Pointer3->SetFillColorAlpha(kGray, 0.05);
        Pointer3->SetLineColor(kGray);

        //draw legend

        TLegend *leg = new TLegend(0.25,0.77,0.26,0.93);
        if(iRow==0)
        {
            leg->AddEntry((TObject*)nullptr,"This thesis", "");
            leg->AddEntry(Pointer1,"5.36 TeV ALICE Pb--Pb 0-10% data","p");
            leg->AddEntry(Pointer2,"5.36 TeV ALICE Pb--Pb 0-30% data (external)","p");
            leg->AddEntry(Pointer3,"CoLBT-hydro ALICE prediction","f");
            leg->AddEntry((TObject*)nullptr,Form("Dijet, R = 0.2, |#eta_{jet}| #leq 0.7, p_{T} #geq (%.1d, %.1d) GeV", 40, 20), "");
            leg->AddEntry((TObject*)nullptr,"|#phi_{Hadron} - #phi_{Leading}| #leq #pi/2", "");
            leg->AddEntry((TObject*)nullptr,"statistical errors only", "");
        }
        leg->AddEntry((TObject*)nullptr,Form("%.1f #leq p_{T, Hadron} #leq %.1f GeV", LowPtCuts[iRow], HighPtCuts[iRow]),"");
        leg->SetLineColor(10);
        leg->SetTextSize(DEF_LegendFontSize); 
        leg->Draw();

    }

    return 1;
}

Int_t Print3DPlots(TString i_InputFile = "in.root")
{
    //CONFIGURE

    gStyle->SetOptStat(0);
    SetRootGraphicStyle();

    #define DEF_BinningPerUnit 100
    #define DEF_AxisLabelSize 0.04
    #define DEF_AxisTitleSize 0.04
    #define DEF_HistoTitleSize 0.1

    int TotalEvents = 0;
    int ProcessedEvents = 0;
    int EventsSkipped   = 0;
    int RunNumber = -1;
    int LargeGapEventCounter = 0;
    int SmallGapEventCounter = 0;

    /*END USER VARIABLES*/

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
        //if(DFName != i_DFName) continue;
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
        vector<PseudoJet> ParticleVector[entries_col];

        // Make an unordered map particle -> collision
        for(uint64_t i_Track = 0; i_Track < entries_track; i_Track ++)
        {
            Particle track = Particle::Read(tracks, i_Track);
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

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /*
        /
        /                                                     
        /                                                           ANALYSIS LOOP 
        /                                                     
        /                                                     
        /                                                     
        /
        /                                                    
        /
        /
        /
        *//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        PrintInfo("Start generating 3D events");

        TFile* OutputFile = new TFile("3DEvents.root","RECREATE");
        OutputFile->cd();

        for(int iEvent = 74; iEvent < entries_col; iEvent++)
        {

            Collision collision = Collision::Read(collisions, iEvent);
            if(CollisionMap[collision.ColID].size() == 0 || ParticleVector[collision.ColID].size() == 0){
                cout << "No tracks in collision " << collision.ColID << "! Skipping." << endl;
                N_events_skipped++;
                continue;
            }

            TH2D* h2D_Phi_vs_Eta = new TH2D("h2D_Phi_vs_Eta,", "h2D_Phi_vs_Eta", 2*Pi*DEF_BinningPerUnit, -Pi, +Pi, 2*0.9*DEF_BinningPerUnit, -0.9, 0.9);

            for(int iTrack = 0; iTrack < ParticleVector[collision.ColID].size(); iTrack++)
            {
                h2D_Phi_vs_Eta->Fill(ParticleVector[collision.ColID][iTrack].phi_std(), ParticleVector[collision.ColID][iTrack].eta(), ParticleVector[collision.ColID][iTrack].pt());
            }

            h2D_Phi_vs_Eta->SetTitle("");
            h2D_Phi_vs_Eta->GetXaxis()->SetTitle("#phi [rad]");
            h2D_Phi_vs_Eta->GetXaxis()->SetTitleSize(DEF_AxisTitleSize);
            h2D_Phi_vs_Eta->GetXaxis()->SetLabelSize(DEF_AxisLabelSize);
            h2D_Phi_vs_Eta->GetXaxis()->SetTitleOffset(1.4);
            h2D_Phi_vs_Eta->GetYaxis()->SetTitle("#eta");
            h2D_Phi_vs_Eta->GetYaxis()->SetTitleSize(DEF_AxisTitleSize);
            h2D_Phi_vs_Eta->GetYaxis()->SetLabelSize(DEF_AxisLabelSize);
            h2D_Phi_vs_Eta->GetYaxis()->SetTitleOffset(1.4);
            h2D_Phi_vs_Eta->GetZaxis()->SetTitle("p_{T} [GeV/c]");
            h2D_Phi_vs_Eta->GetZaxis()->SetTitleSize(DEF_AxisTitleSize);
            h2D_Phi_vs_Eta->GetZaxis()->SetLabelSize(DEF_AxisLabelSize);

            h2D_Phi_vs_Eta->Write();
            h2D_Phi_vs_Eta->DrawCopy("lego2, fb");

            delete(h2D_Phi_vs_Eta);

            break;

        }

        //only do the first DF, more is not necessary
        break;
    
    }

    return 1;
}

Int_t Analysis_SimulationStats()
{
    #define DEF_AxisLabelSize 0.09
    #define DEF_AxisTitleSize 0.09
    #define DEF_HistoTitleSize 0.1
    #define DEF_LegendTextSize 0.07
    #define DEF_Rebin 16
    #define DEF_Margin_Top 0.05
    #define DEF_Margin_Bottom 0.20
    #define DEF_Margin_Left 0.1
    #define DEF_Margin_Right 0.1
    #define DEF_TitleOffset 0.5

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
    Canvas->Divide(1,5, 0, 0);

    //Multiplicity
    Canvas->cd(1);
    gPad->SetMargin(DEF_Margin_Left, DEF_Margin_Right, DEF_Margin_Bottom, DEF_Margin_Top);
    h1D_NDist->SetTitle("");
    h1D_NDist->GetXaxis()->SetTitleSize(DEF_AxisTitleSize);
    h1D_NDist->GetXaxis()->SetLabelSize(DEF_AxisLabelSize);
    h1D_NDist->GetXaxis()->SetTitle("Multiplicity");
    h1D_NDist->GetXaxis()->SetRangeUser(1600, 3500);
    h1D_NDist->GetYaxis()->SetTitleSize(DEF_AxisTitleSize);
    h1D_NDist->GetYaxis()->SetLabelSize(DEF_AxisLabelSize);
    h1D_NDist->GetYaxis()->SetTitleOffset(DEF_TitleOffset);
    h1D_NDist->GetYaxis()->SetTitle("yield");
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
    gPad->SetMargin(DEF_Margin_Left, DEF_Margin_Right, DEF_Margin_Bottom, DEF_Margin_Top);
    h1D_PtDist->SetTitle("");
    h1D_PtDist->GetXaxis()->SetTitleSize(DEF_AxisTitleSize);
    h1D_PtDist->GetXaxis()->SetLabelSize(DEF_AxisLabelSize);
    h1D_PtDist->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    h1D_PtDist->GetYaxis()->SetTitleSize(DEF_AxisTitleSize);
    h1D_PtDist->GetYaxis()->SetLabelSize(DEF_AxisLabelSize);
    h1D_PtDist->GetYaxis()->SetTitle("yield");
    h1D_PtDist->GetYaxis()->SetTitleOffset(DEF_TitleOffset);
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
    gPad->SetMargin(DEF_Margin_Left, DEF_Margin_Right, DEF_Margin_Bottom, DEF_Margin_Top);
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
    Graph_v2_pT->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    Graph_v2_pT->GetXaxis()->SetRangeUser(0,20);
    Graph_v2_pT->GetYaxis()->SetTitleSize(DEF_AxisTitleSize);
    Graph_v2_pT->GetYaxis()->SetLabelSize(DEF_AxisLabelSize);
    Graph_v2_pT->GetYaxis()->SetTitle("v_{2}");
    Graph_v2_pT->GetYaxis()->SetTitleOffset(DEF_TitleOffset);
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
    gPad->SetMargin(DEF_Margin_Left, DEF_Margin_Right, DEF_Margin_Bottom, DEF_Margin_Top);
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
    Graph_Orig->GetYaxis()->SetTitleOffset(DEF_TitleOffset);
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
    leg4->AddEntry((TObject*)nullptr,Form("0.2 #leq p_{T}^{track} #leq 5.0 GeV/c"),"");
    leg4->AddEntry(Graph_Orig,"Data","p");
    leg4->AddEntry(fit,"Fit","l");
    leg4->SetLineColor(10);
    leg4->Draw();

    Canvas->cd(5);
    gPad->SetMargin(DEF_Margin_Left, DEF_Margin_Right, DEF_Margin_Bottom, DEF_Margin_Top);
    #define XErr 0.16
    Double_t OrigX2[5] = {-0.64, -0.32, 0.0, 0.32, 0.64};
    Double_t OrigY2[5] = {5.843e-4, 2.641e-4, 4.290e-5, -1.414e-4, -5.401e-4};
    Double_t eX2[5] = {XErr, XErr, XErr, XErr, XErr};
    Double_t eY2[5] = {1.626e-4, 1.282e-4, 1.333e-4, 1.280e-4, 1.221e-4};
    TGraphErrors* Graph_Orig2 = new TGraphErrors(5, OrigX2, OrigY2, eX2, eY2);
    Graph_Orig2->SetMarkerColor(kBlue);
    Graph_Orig2->SetTitle("");
    Graph_Orig2->GetXaxis()->SetTitleSize(DEF_AxisTitleSize);
    Graph_Orig2->GetXaxis()->SetLabelSize(DEF_AxisLabelSize);
    Graph_Orig2->GetXaxis()->SetTitle("#eta");
    //Graph_Orig2->GetXaxis()->SetRangeUser(-1.5, 1.5);
    Graph_Orig2->GetYaxis()->SetTitleSize(DEF_AxisTitleSize);
    Graph_Orig2->GetYaxis()->SetLabelSize(DEF_AxisLabelSize);
    Graph_Orig2->GetYaxis()->SetTitle("v_{1}");
    Graph_Orig2->GetYaxis()->SetTitleOffset(DEF_TitleOffset);
    Graph_Orig2->SetLineColor(kBlue);
    Graph_Orig2->Draw("A*");

    auto LinFit = new TF1("lin", "[0]*x", -0.64, 0.64);
    LinFit->SetParameters(0. -8.4e-4);
    // create graph
    Graph_Orig2->Fit("lin");
    gPad->Update();

    TF1 *fit2 = (TF1*)Graph_Orig2->GetListOfFunctions()->FindObject("lin");
    fit2->SetLineColor(kBlack);

    TLegend *leg5 = new TLegend(0.0,0.8,0.3,1.0);
    leg5->SetTextSize(DEF_LegendTextSize);
    leg5->AddEntry((TObject*)nullptr,Form("10-20%% ALICE Pb-Pb"),"");
    leg5->AddEntry((TObject*)nullptr,Form("#sqrt{s_{NN}} = 2.76 TeV"),"");
    leg5->AddEntry((TObject*)nullptr,Form("0.15 #leq p_{T}^{track} #leq 20.0 GeV/c"),"");
    leg5->AddEntry(Graph_Orig2,"Data","p");
    leg5->AddEntry(fit2,"Fit","l");
    leg5->SetLineColor(10);
    leg5->Draw();


    return 1;
}

Int_t DijetAna(const TString DataFile, float DiffRange, bool IsJetMeasurement, int LeadPt, int SublPt) {

    #define DEF_AxisLabelSize 0.07
    #define DEF_AxisTitleSize 0.08
    #define DEF_HistoTitleSize 0.1
    #define DEF_Rebin 20
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
        if(iRow==0)
        {
            leg->AddEntry((TObject*)nullptr,"This thesis","");
            leg->AddEntry((TObject*)nullptr,"ALICE Pb-Pb 0-10%  #sqrt{s_{NN}}=5.36 TeV","");
            if(IsJetMeasurement) leg->AddEntry((TObject*)nullptr,Form("Dijet, R = 0.2, |#eta_{jet}| #leq 0.7, p_{T} #geq (%.1d, %.1d) GeV/c", LeadPt, SublPt), "");
            else leg->AddEntry((TObject*)nullptr, Form("Dihadron, |#eta_{Lead./Subl. Hadron}| #leq 0.7, p_{T} #geq (%.1d, %.1d) GeV/c", LeadPt, SublPt), "");
            leg->AddEntry((TObject*)nullptr,"|#phi_{Hadron} - #phi_{Leading}| #leq #pi/2", "");
            leg->AddEntry((TObject*)nullptr,"statistical errors only", "");
        }
        leg->AddEntry((TObject*)nullptr,Form("%.1f #leq p_{T, Hadron} #leq %.1f GeV/c", LowPtCuts[iRow], HighPtCuts[iRow]),"");
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
        h1F_LargeGap[iRow]->DrawCopy("");

        TLine *line = new TLine(-.9, 0.0, .9, 0.0);
        line->SetLineColor(kGray);
        line->Draw("SAME");

        h1F_LargeGap[iRow]->DrawCopy("same");

        //get chi^2 of null hypothesis test
        int NBins = 5;//h1F_LargeGap[iRow]->GetXaxis()->GetNbins();
        float Sum = 0;

        for(int iBin = 1; iBin <= NBins; iBin++)
        {
            float Num = pow(h1F_LargeGap[iRow]->GetBinContent(iBin), 2);
            float Den = pow(h1F_LargeGap[iRow]->GetBinError(iBin),2);
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

    return 1;

}