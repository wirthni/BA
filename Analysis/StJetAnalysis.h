#ifndef StJetAnalysis_hh
#define StJetAnalysis_hh

#define USEEVE

//#include "StarClassLibrary/SystemOfUnits.h" // ALEX
#include "TString.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCutG.h"
#include "TCanvas.h"
#include "TTree.h"
//#include "StMessMgr.h" // ALEX
//#include <iostream.h>
#include "TChain.h"
#include "TNtuple.h"
#include "TROOT.h"
#include "TFile.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TStyle.h"
#include "TGaxis.h"
#include <fstream>
#include "TMath.h"
#include "TColor.h"
#include "TF1.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TExec.h"
#include "TPolyMarker.h"
#include "TVirtualPad.h"
#include "TPolyLine.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TH3D.h"
#include "TPolyLine3D.h"
#include "TPolyMarker3D.h"
#include "TGLViewer.h"
#include "TGLSAViewer.h"
#include "TGLCamera.h"
#include "TGLPerspectiveCamera.h"
#include "TGFrame.h"
#include "TGLUtil.h"
#include "TGLLightSet.h"
#include "TGLCameraOverlay.h"
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoVolume.h"
#include "TGeoMatrix.h"

//#include "StJetAnalysisLinkDef.h"



using namespace std;

//#include "/home/ceres/schmah/ALICE/Jet_Analysis/StJetTrackEventLinkDef.h"
//#include "/home/ceres/schmah/ALICE/Jet_Analysis/StJetTrackEvent.h"
#include "../StJetTrackEvent.h"
//#include "StJetTrackEvent.h"
#include "../StJetTrackEventLinkDef.h"
#include "../PythiaEvent.h"
#include "../PythiaEvent_LinkDef.h"

ClassImp(StEMCal)
ClassImp(StJetTrackParticle)
ClassImp(StJetTrackEvent)
ClassImp(PythiaEvent)
ClassImp(PythiaParticle)
// #include "StarClassLibrary/StThreeVectorF.hh" // ALEX
//#include "StRoot/StRefMultCorr/StRefMultCorr.h" // ALEX
// #include "StRoot/StRefMultCorr/CentralityMaker.h" // ALEX

#include "fastjet/config.h"             // will allow a test for FJ3
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/Subtractor.hh"


//#include "StPhysicalHelixD.hh" // ALEX

using namespace fastjet;
using namespace std;


#include "fastjet/tools/JetMedianBackgroundEstimator.hh"

#if defined(USEEVE)
#include "TEveBox.h"
#include "TEveArrow.h"
#include <TEveManager.h>
#include "TEveLine.h"
#include "TEvePointSet.h"
#endif

typedef std::vector< std::vector<Int_t> >  Two_dim_Int_vector;
typedef std::vector< std::vector< std::vector<Int_t> > >  Three_dim_Int_vector;
typedef std::vector< std::vector< std::vector<StJetTrackEvent> > >  Three_dim_StJetTrackEvent_vector;

static const Double_t Pi = TMath::Pi();


#define KNRM  "\x1B[0m"
#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"
#define KYEL  "\x1B[33m"
#define KBLU  "\x1B[34m"
#define KMAG  "\x1B[35m"
#define KCYN  "\x1B[36m"
#define KWHT  "\x1B[37m"



class StJetAnalysis {
public:
    StJetAnalysis();
    ~StJetAnalysis();



    vector<Double_t> Get_Helix_params_from_kine(TLorentzVector TLV_particle, TVector3 TV3_vertex, Double_t charge);
    void  setGeomDir(TString iGeomDir) {GeomDir = iGeomDir;};
    void  setInListDir(TString iInListDir) {InListDir = iInListDir;};
    void  setSEList(TString iSEList) {SEList = iSEList;};
    void  setSEListPYTHIA(TString iSEListPYTHIA) {SEListPYTHIA = iSEListPYTHIA;};
    void  setInputDir(TString ipinputdir) {pinputdir = ipinputdir;};
    void  setInputDirPYTHIA(TString ipinputdirPYTHIA) {pinputdirPYTHIA = ipinputdirPYTHIA;};
    void  setStartEvent(Long64_t istart_event) {start_event_use = istart_event;};
    void  setStopEvent(Long64_t istop_event) {stop_event_use = istop_event;};
    void  setOutputfile(TString ioutputfile_name) {outputfile_name = ioutputfile_name;};
    void  setEnergy(Int_t iBeamTimeNum) {eBeamTimeNum = iBeamTimeNum;};
    void  setRandom(Int_t iRandom) {eRandom = iRandom;};
    void  setMode(Int_t iMode) {eMode = iMode;};
    void  setIn_Mode(Int_t iIn_Mode) {eIn_Mode = iIn_Mode;};
    void  setOutdir(TString iOutdir) {eOutdir = iOutdir;};
    void  setOutFileName(TString iOutFileName) {eOutFileName = iOutFileName;};
    void  setSampleHisto(TString iSampleHist) {eSampleHist = iSampleHist;};
    void  setSE_Et_Histo(TString iSE_Et_Hist) {eSE_Et_Hist = iSE_Et_Hist;};
    void  setME_Et_Histo(TString iME_Et_Hist) {eME_Et_Hist = iME_Et_Hist;};
    void  setSuffix(TString iSuffix) {eSuffix = iSuffix;};
    void  setReCenteringFile(TString ire_centering_name) {re_centering_name = ire_centering_name;};
    void  setPYTHIA_eff_factor(Double_t iPYTHIA_eff_factor) {ePYTHIA_eff_factor = iPYTHIA_eff_factor;};
    void  setBins(Int_t iz_bin, Int_t imult_bin, Int_t iPsi_bin) {ez_bin = iz_bin; emult_bin = imult_bin; ePsi_bin = iPsi_bin;};
    void  setCentrality(Int_t iCentrality) {eCentrality = iCentrality;};
    void  setPrim_Glob(Int_t iflab_prim_glob) {eflab_prim_glob = iflab_prim_glob;};
    void  setJet_R(Double_t iJet_R) {eJet_R  = iJet_R;};
    void  setDoJetR(Int_t idoJetR)  {doRjet = idoJetR;};
    void  setBkg_R(Double_t iBkg_R) {eBkg_R  = iBkg_R;};
    void  setRem_n_hardest(Int_t iRem_n_hardest) {eRem_n_hardest = iRem_n_hardest;};
    void  set_Trigger_pt(Double_t iTrigger_pt) {eTrigger_pt = iTrigger_pt;};
    void  set_N_vertex_mult_Psi_bins(Int_t N_z_vertex_bins_set, Int_t N_mult_bins_set, Int_t N_Psi_bins_set) {N_z_vertex_bins = N_z_vertex_bins_set; N_mult_bins = N_mult_bins_set; N_Psi_bins = N_Psi_bins_set;};
    void  set_N_sub_Psi_bins(Int_t set_EP_sub_split) {EP_sub_split = set_EP_sub_split;};
    void  set_N_sub_z_vertex_bins(Int_t set_z_vertex_sub_split) {z_vertex_sub_split = set_z_vertex_sub_split;};
    void  set_ME_p_splitting(Double_t set_p_threshold, Int_t set_N_split, Int_t set_flag_equiv_split) {p_threshold = set_p_threshold; N_split = set_N_split; flag_equiv_split = set_flag_equiv_split;};
    void  set_DoPythia(Int_t set_do_pythia) {doPYTHIA = set_do_pythia;};
    void  set_DoFullSim(Int_t set_do_full_sim) {doFullSim = set_do_full_sim;};
    void  set_DoRdet(Int_t set_doRdet) {doRdet = set_doRdet;};
    void  set_DoEmbeddingPythiaJets(Int_t set_doEmbeddingPythiaJets, Int_t set_doParticleOrDetectorLevel, Int_t set_doEmbeddingInSE_ME) {doEmbeddingPythiaJets = set_doEmbeddingPythiaJets; doParticleOrDetectorLevel = set_doParticleOrDetectorLevel; doEmbeddingInSE_ME = set_doEmbeddingInSE_ME;};
    void  set_DoEmbeddingParticleOrDetectorLevelJets(Int_t set_doEmbeddingParticleOrDetectorLevelJets) {doEmbeddingParticleOrDetectorLevelJets = set_doEmbeddingParticleOrDetectorLevelJets;};
    void  set_PYTHIA_tree(Int_t set_PYTHIA_tree) { file_select_use = set_PYTHIA_tree;};
    void  set_syst_unc_reco_eff(Int_t set_syst_unc_recon_eff) { syst_unc_recon_eff = set_syst_unc_recon_eff;};
    void  set_DoMC(Int_t set_do_MC) {doMC = set_do_MC;};
    void  setDCA_cut(Double_t iDCA_cut) {eDCA_cut = iDCA_cut;};
    void  setDoCout(Int_t set_doCout) {doCout = set_doCout;};
    void  setThreshMatch(Double_t set_threshold_match) {threshold_match = set_threshold_match;};
    void  setDoWidthScaling(Int_t set_do_width_scaling) {doWidthScalingME = set_do_width_scaling;};
    //void  setMEScalingFactor(Double_t set_ME_scaling_factor) {MEScalingFactor = set_ME_scaling_factor;};
    void  setMEScalingFactor(vector<vector<Double_t>> set_ME_scaling_factor) {MEScalingFactor = set_ME_scaling_factor;};
    void  setMEScalingShift(vector<vector<Double_t>> set_ME_scaling_shift) {MEScalingShift = set_ME_scaling_shift;};
    void  setHighpTTrackCut(Double_t set_HighpTTrackCut) {HighpTTrackCut = set_HighpTTrackCut;};
    void  setDoV0EP(Int_t set_DoV0EP) {doEPV0 = set_DoV0EP;};


    void  SplitTrees(Int_t s_z_vertex_bin, Int_t s_mult_bin, Int_t s_Psi_bin);
    void  MakeJets(Int_t graphics, Int_t ME_flag, Int_t EP_bin, Int_t i_jetR, Int_t Psi3_bin, Int_t pT_sum_bin);
    void  MakeME(Int_t i_ME_z_vertex_bin,Int_t i_ME_mult_bin, Int_t i_ME_EP_bin, Int_t i_ME_EP_sub_bin, Int_t i_ME_z_vertex_sub_bin);
    //Double_t MakeEventPlane(Int_t do_loop_event);
    void  MakeEventPlane(Int_t do_loop_event, Double_t &Psi_full_r,Double_t &Psi_full_r_v3, Int_t &N_particles_TPC_refit_r, Int_t doEP_v2, Int_t doEP_v3);
    Int_t LoopEvent(Int_t LoopEvent, Int_t graphics, Int_t do_MakeJets);
    Int_t LoopEventMC(Int_t LoopEvent, Int_t graphics, Int_t do_MakeJets);
    Int_t LoopEventMCpp(Int_t LoopEvent,Int_t do_MakeJets);
    void  setMultHist(TString iMultHistName) {eMultHistName = iMultHistName;}
    Int_t LoopEventPYTHIA(Int_t LoopEventPythia);

    Int_t Get_Mult_bin(Int_t N_Particles, Int_t i_z_vertex);
    Int_t Get_EP_bin();
    Int_t Get_Psi3_bin();
    Int_t Get_v2_bin(Double_t v2_value, Int_t EP_bin, Int_t mult_bin, Int_t zvertex_bin);
    Int_t Get_z_vertex_bin(Float_t prim_z_vertex);
    Int_t Get_sub_EP_bin(Int_t EP_split, Int_t EP_bin);
    Int_t Get_sub_z_vertex_bin(Float_t prim_z_vertex, Int_t z_vertex_split, Int_t z_vertex_bin);
    Int_t Get_pTsum_bin(Int_t mult_bin, Int_t zvertex_bin);
    //----------------------------------------------------------------------
    // Function for track reconstruction efficiency
    void      Read_ABC_params_eff_track_rec_function();             //read root file  for Global_ABC
    Double_t* Parameter_eff_track_rec_function(Int_t cent, Int_t runID,  Int_t PID, Int_t use_all);  //return A,B,C,runID for different centraltiy
    void      setEff_track_rec_parameter_file(TString iEff_file) {eEff_file = iEff_file;};
    Int_t     Get_runID_range_Index(Int_t runid);
    //----------------------------------------------------------------------


    void Init3DGraphics();
    void InitSplit(Int_t s_z_vertex_bin, Int_t s_mult_bin, Int_t s_Psi_bin);
    void InitJet();
    Long64_t ReadData();
    Long64_t InitPythiaTree(Int_t file_select);
    Long64_t InitMCTree();
    void WriteJet();
    Int_t Make();
    Int_t Finish();
    void Plot_test();
    void NormME();
    void DoHistograms();
    void DoEmbeddingPythiaJets();
    Int_t DoParticleLevel();
    void DoDetectorLevel();


    Double_t CalculateEventPlaneChi(int iDet, int iCentr){
    std::vector<std::vector<double>> res(2,std::vector<double>(8));

    //iDet=0 -> V0A
    res[0][0] = 0.393689;
    res[0][1] = 0.559982;
    res[0][2] = 0.67293;
    res[0][3] = 0.693582;
    res[0][4] = 0.65589;
    res[0][5] = 0.572881;
    res[0][6] = 0.368911;
    res[0][7] = 0.133216;

    //iDet=1 -> V0C
    res[1][0] = 0.489577;
    res[1][1] = 0.659353;
    res[1][2] = 0.76015;
    res[1][3] = 0.784079;
    res[1][4] = 0.750566;
    res[1][5] = 0.673186;
    res[1][6] = 0.459845;
    res[1][7] = 0.17636;

    Double_t chi (2.) , delta (1.) ;
    Double_t con (( TMath :: Sqrt ( TMath :: Pi () ) ) /(2.* TMath :: Sqrt (2) ) ) ;
    for ( Int_t i(0) ; i < 15; i++) {
      chi = (( con *chi * TMath :: Exp ( -chi *chi /4.) *( TMath :: BesselI0 (chi *chi /4.) + TMath :: BesselI1 (chi*chi /4.) ) ) < res[iDet][iCentr]) ? chi + delta : chi - delta ;
      delta = delta / 2.;
    }
    return chi;
  };

    //Double_t Split_function(Double_t x_val, Double_t par);


private:

    TRandom ran;
    TRandom ran_mult;
    TRandom3 r3;
    TString HistName;
    char NoP[50];
    TRandom ran_gen;
    Int_t N_z_vertex_bins, N_mult_bins, N_Psi_bins;
    Double_t p_threshold;
    Int_t N_split;
    Int_t flag_equiv_split;
    Int_t EP_sub_split;
    Int_t z_vertex_sub_split;
    Int_t doPYTHIA;
    Int_t doMC;
    Int_t doFullSim;
    Int_t doRdet;
    Int_t doEmbeddingInSE_ME;
    Int_t doEmbeddingPythiaJets;
    Int_t doParticleOrDetectorLevel;
    Int_t i_event_pythia_embedding = 0;
    Int_t doEmbeddingParticleOrDetectorLevelJets;
    Int_t N_events_used_for_jet_analysis;
    Double_t parameter_event[6];
    vector<Double_t> vec_dca_xy;
    vector<Double_t> vec_dca_z;
    Int_t flag_SEorME = 0;
    Int_t doEPV0;

    Int_t doWidthScalingME;
    //Double_t MEScalingFactor;
    vector<vector<Double_t>> MEScalingFactor;
    vector<vector<Double_t>> MEScalingShift;

    Int_t doCout;

    Double_t HighpTTrackCut;

    vector< vector< vector<StJetTrackEvent*> > > JetTrackEvent_Fill;
    vector< vector< vector<TTree*> > >           Tree_JetTrackEvent_Fill;
    vector< vector< vector<TFile*> > >           vec_Outputfiles;
    TFile* outputfile_jet;
    TFile* Outputfile_Histo;

    TVector3 Qvec_eta_pos, Qvec_eta_neg;
    vector< vector<TVector3> > vec_TV3_Qvec_eta;
    Double_t Psi_full_global;
    Double_t Psi_full_global_v3;
    Double_t pT_sum_SE_global;
    Double_t Psi_pos_neg[2];
    Double_t Psi_pos_neg_eta[2];
    Int_t N_EP_bins;
    Int_t N_particles_TPC_refit_used;
    Int_t mult_used_for_jets;
    Int_t N_radius_parameters_used;

    vector<PseudoJet> vec_PJ_particles;
    vector<PseudoJet> vec_PJ_particles_recon;
    vector<PseudoJet> jets_pythia;
    vector<PseudoJet> vec_PJ_partiles_truth;
    vector<PseudoJet> vec_PJ_partiles_detector;

    Int_t eflab_prim_glob, eRem_n_hardest;
    Double_t eJet_R, eBkg_R, eTrigger_pt, Jet_area_cut, doRjet;
    Double_t jet_rho;
    Double_t eDCA_cut;
    TString eSuffix, InListDir, GeomDir, eMultHistName;
    Double_t ePYTHIA_eff_factor;
    Int_t ez_bin, emult_bin, ePsi_bin, eCentrality;
    TString eOutdir, eOutFileName, eSampleHist, eSE_Et_Hist, eME_Et_Hist;
    Int_t eMode, eIn_Mode;
    Int_t eBeamTimeNum;
    Int_t eRandom;
    Long64_t file_entries_SE_ME[2];
    TString SEList,SEListPYTHIA, pinputdir, pinputdirPYTHIA, outputfile_name, re_centering_name, OutputfileName;
    TChain* input_SE_ME[2];

    //-----------------------------------
    // PYTHIA
    // PYTHIA response matrix calc.
    TChain           *inputPYTHIA;
    PythiaEvent      *PYTHIAEvent;
    PythiaParticle   *PYTHIAParticle;
    Long64_t          file_entries_PYTHIA;
    Long64_t          file_entries_MCpp;
    
    TH1D* h_jet_pT_inclusive_pythia;
    TH1D* h_jet_pT_biased_pythia;
    TH1D* h_jet_pT_pythia_max_pT;
    TH1D* h_delta_pT_embed_pythia_jets;

    TH1D* h_delta_pT_embed_pythia_jets_7GeV;
    TH1D* h_delta_pT_embed_pythia_jets_6GeV;
    TH1D* h_delta_pT_embed_pythia_jets_5GeV;
    TH1D* h_delta_pT_embed_pythia_jets_4GeV;
    TH1D* h_delta_pT_embed_pythia_jets_3GeV;
    TH1D* h_delta_pT_embed_pythia_jets_2_5GeV;
    TH1D* h_delta_pT_embed_pythia_jets_2GeV;

    TH2D* h2D_pT_particle_pythia_jets_vs_pT_recon;

    TH2D* h2D_pT_particle_pythia_jets_vs_pT_recon_7GeV;
    TH2D* h2D_pT_particle_pythia_jets_vs_pT_recon_6GeV;
    TH2D* h2D_pT_particle_pythia_jets_vs_pT_recon_5GeV;
    TH2D* h2D_pT_particle_pythia_jets_vs_pT_recon_4GeV;
    TH2D* h2D_pT_particle_pythia_jets_vs_pT_recon_3GeV;
    TH2D* h2D_pT_particle_pythia_jets_vs_pT_recon_2_5GeV;
    TH2D* h2D_pT_particle_pythia_jets_vs_pT_recon_2GeV;

    TH1D* h_jet_pT_7GeV_bias_pythia;
    TH1D* h_jet_pT_6GeV_bias_pythia;
    TH1D* h_jet_pT_5GeV_bias_pythia;
    TH1D* h_jet_pT_4GeV_bias_pythia;
    TH1D* h_jet_pT_3GeV_bias_pythia;
    TH1D* h_jet_pT_2_5GeV_bias_pythia;
    TH1D* h_jet_pT_2GeV_bias_pythia;

    TH1D* h_jet_pT_pythia_max_pT_bias[7];

    TH2D* h2D_particle_level_vs_detector_level;
    TH2D* h2D_pythia_particle_level_jet_pT_vs_area;
    TH2D* h2D_pythia_particle_level_jet_pT_vs_mult;

    vector<TLorentzVector> vec_tlv_particle_level_jet;
    vector<TLorentzVector> vec_tlv_pythia_particles;

    TGraphAsymmErrors* tgae_efficiency;
    TGraphAsymmErrors* tgae_pT_smearing;

    TGraph* tg_pileup_TOF;
    TGraph* tg_pileup_V0A;
    TGraph* tg_pileup_V0C;

    vector<Int_t> vec_user_index_particle_level;

    TProfile* h_profile_jet_reco_eff;
    vector<TProfile*> vec_h_profile_jet_reco_eff;
    TProfile* profile_efficiency; // input

    TH1D* h_eff_uncertainty; // input

    TRandom ran_pythia;

    Double_t pT_jet_pythia_embedded_global;
    Double_t pT_jet_pythia_embedded_global_detector_level;
    Double_t pT_jet_pythia_reconstructed;
    Double_t eta_jet_pythia_reconstructed;
    Double_t phi_jet_pythia_reconstructed;
    Double_t px_jet_pythia_reconstructed;
    Double_t py_jet_pythia_reconstructed;
    Double_t pz_jet_pythia_reconstructed;
    Int_t syst_unc_recon_eff;
    Int_t doDetectorOrParticleLevel;
    Int_t threshold_embedding;
    Int_t file_select_use;
    Int_t real_event;
    Double_t threshold_match;

    TH1D* h_N_events_PYTHIA;

    Int_t particle_level_done;

    //-----------------------------------
    // Closure test

    
    //-----------------------------------
    // MC data
    TProfile* p_eff_vs_pt;
    TProfile2D* p_eff_eta_vs_phi[5];
    TH2D* h2D_MC_pT_vs_recon_pT;
    TH2D* h2D_delta_pt_vs_pt;
    TH1D* h_eta_phi_diff_best;
    TH1D* h_eta_phi_diff_all;
    TH1D* h_eta_diff_best;
    TH1D* h_phi_diff_best;
    TH1D* h_pT_diff_best;
    TH1D* h_phi_MC_all;
    vector<TLorentzVector> tlv_reco_particle;
    vector<TLorentzVector> tlv_MC_particle;
    vector<TH1D*> vec_h_MC_particle_vertex;
    TH1D* h_pT_gen;
    TH2D* h2D_MC_pT_vs_recon_pT_TEST;
    TH1D* h_eta_phi_diff_best_TEST;
    TH2D* h2D_eta_phi_truth;
    TH2D* h2D_eta_phi_detector;
    TH1D* h_pt_truth;
    TH1D* h_pt_detector;
    //-----------------------------------


    TChain* input_PYTHIA[11]; // PYTHIA hard bins, only for mode 312 -> embedding
    StJetTrackEvent     *JetTrackEvent_PYTHIA[11]; // PYTHIA hard bins, only for mode 312 -> embedding
    StJetTrackParticle  *JetTrackParticle_PYTHIA[11]; // PYTHIA hard bins, only for mode 312 -> embedding
    Long64_t start_event_use, stop_event_use;
    TFile* Outputfile;
    TFile* Inputfile;
    TFile* Inputfile_Et[2];
    StJetTrackEvent     *JetTrackEvent;
    StJetTrackParticle  *JetTrackParticle;
    StMCParticle        *JetTrackMCParticle;
    StEMCal             *JetEMCalParticle;
    StJetTrackParticle  *JetTrackParticle_Fill;
    StEMCal             *JetEMCalParticle_Fill;
    // mixing
    vector<Int_t> vec_ME_mult_bins; // multiplicity
    vector<Int_t> vec_bin_limits;   // bin limits for mult
    vector<vector<Int_t>> vec_ME_mult_bins_new; // multiplicity
    vector<vector<Int_t>> vec_bin_limits_new;   // bin limits for mult
    vector<Int_t> vec_ME_mult_used_bins; // actually number of particles after cuts used multiplicity
    vector<Int_t> vec_ME_EP_bins; // event plane
    vector<Int_t> vec_ME_EP_Psi3_bins; // event plane Psi3 bins
    vector<Int_t> vec_ME_z_vertex_bins; // z-veretx
    vector<vector<Int_t>> vec_ME_pT_sum_bins; // pT sum bins
    vector<vector<vector<vector<Double_t>>>> vec_ME_v2_bins; // v2 bins
    vector<vector<vector<TH1D*>>> vec_h_v2_bins; // histo v2 bins
    vector<vector<vector<Double_t>>> vec_ME_pT_sum_bins_new;  // new pT sum bins
    vector< vector< vector< vector<Int_t> > > > vec_ME_event_index_mult;
    vector< vector< vector< vector< vector< vector<Int_t> > > > > > vec_ME_event_index_mult_EP_sub_bins;
    vector<Double_t> vec_trigger_pT;
    vector<vector<vector<TH1D*>>>  vec_h_mult_SE_input; // input mult SE for mixing
    vector<vector<vector<vector<vector<TH1D*>>>>>  vec_h_mult_SE_input_new; // input mult SE for mixing new in pT sum and Psi3
    vector<TH1D*>  vec_h_mult_SE_input_tree_splitting;
    vector<TH2D*>  vec_h2D_mult_dist_vs_pT_sum;

    TProfile* h_parameters;
    TProfile* h_profile_v2_track_SE;
    TProfile* h_profile_v2_track_ME;
    TProfile* h_profile_v2_track_SE_shift;
    TProfile* h_profile_v2_track_ME_shift;
    TH2D* h2D_dEdx_vs_p;
    //-----------------------------------
    // Jet histos
    vector<vector<vector<TH1D*>>> vec_h_jet_sub_SE;
    vector<vector<vector<TH1D*>>> vec_h_jet_sub_ME;
    vector<vector<vector<TH1D*>>> vec_h_h_jet_SE_project;
    vector<vector<vector<TH1D*>>> vec_h_h_jet_ME_project;
    vector<vector<vector<TH1D*>>> vec_h_jet_SE;
    vector<vector<vector<TH1D*>>> vec_h_jet_ME;

    vector<vector<vector<TH1D*>>> vec_h_jet_pT_sub_SE_7GeV;
    vector<vector<vector<TH1D*>>> vec_h_jet_pT_sub_ME_7GeV;
    vector<vector<vector<TH1D*>>> vec_h_jet_pT_sub_SE_6GeV;
    vector<vector<vector<TH1D*>>> vec_h_jet_pT_sub_ME_6GeV;
    vector<vector<vector<TH1D*>>> vec_h_jet_pT_sub_SE_5GeV;
    vector<vector<vector<TH1D*>>> vec_h_jet_pT_sub_ME_5GeV;
    vector<vector<vector<TH1D*>>> vec_h_jet_pT_sub_SE_4GeV;
    vector<vector<vector<TH1D*>>> vec_h_jet_SE_4GeV;
    vector<vector<vector<TH1D*>>> vec_h_jet_pT_sub_ME_4GeV;
    vector<vector<vector<TH1D*>>> vec_h_jet_pT_sub_SE_3GeV;
    vector<vector<vector<TH1D*>>> vec_h_jet_pT_sub_ME_3GeV;
    vector<vector<vector<TH1D*>>> vec_h_jet_pT_sub_SE_2GeV;
    vector<vector<vector<TH1D*>>> vec_h_jet_pT_sub_ME_2GeV;

    vector<vector<vector<TH1D*>>> vec_h_jet_pT_sub_SE_1GeV;
    vector<vector<vector<TH1D*>>> vec_h_jet_pT_sub_ME_1GeV;
    vector<vector<vector<TH1D*>>> vec_h_jet_pT_sub_SE_2_5GeV;
    vector<vector<vector<TH1D*>>> vec_h_jet_pT_sub_ME_2_5GeV;

    vector<vector<vector<TH1D*>>> vec_h_jet_pT_sub_ME_between_1_10GeV;
    vector<vector<vector<TH1D*>>> vec_h_jet_pT_sub_ME_between_2_10GeV;
    vector<vector<vector<TH1D*>>> vec_h_jet_pT_sub_ME_between_2_5_10GeV;
    vector<vector<vector<TH1D*>>> vec_h_jet_pT_sub_ME_between_3_10GeV;
    vector<vector<vector<TH1D*>>> vec_h_jet_pT_sub_ME_between_4_10GeV;
    vector<vector<vector<TH1D*>>> vec_h_jet_ME_between_4_10GeV;
    vector<vector<vector<TH1D*>>> vec_h_jet_pT_sub_ME_between_5_10GeV;
    vector<vector<vector<TH1D*>>> vec_h_jet_pT_sub_ME_between_6_10GeV;
    vector<vector<vector<TH1D*>>> vec_h_jet_pT_sub_ME_between_7_10GeV;

    vector<vector<TH1D*>> vec_h_rho_sub_SE;
    vector<vector<TH1D*>> vec_h_rho_sub_ME;

    vector<TH1D*> h_area_sub_SE;
    vector<TH1D*> h_area_sub_ME;
    vector<TH1D*> h_area_SE;
    vector<TH1D*> h_area_ME;

    vector<vector<vector<TH1D*>>> vec_h_jet_sub_SE_wo_area_cut;
    vector<vector<vector<TH1D*>>> vec_h_jet_sub_SE_wo_area_cut_raw;
    vector<vector<vector<TH1D*>>> vec_h_jet_sub_SE_wo_high_pt_cut;
    vector<vector<vector<TH1D*>>> vec_h_jet_SE_wo_high_pt_cut;
    vector<vector<vector<TH1D*>>> vec_h_jet_sub_SE_w_high_pt_cut_200;
    vector<vector<vector<TH1D*>>> vec_h_jet_sub_SE_w_high_pt_cut_100;

    // incl. tests
    vector<vector<vector<TH1D*>>> vec_h_jet_pT_ME_without_10GeV;
    vector<vector<vector<TH1D*>>> vec_h_jet_pT_ME_without_5GeV;
    vector<vector<vector<TH1D*>>> vec_h_jet_pT_ME_without_4_5GeV;
    vector<vector<vector<TH1D*>>> vec_h_jet_pT_ME_without_4GeV;
    vector<vector<vector<TH1D*>>> vec_h_jet_pT_ME_without_3_5GeV;
    vector<vector<vector<TH1D*>>> vec_h_jet_pT_ME_without_3GeV;
    vector<vector<vector<TH1D*>>> vec_h_jet_pT_ME_without_2_5GeV;
    vector<vector<vector<TH1D*>>> vec_h_jet_pT_ME_without_2GeV;

    vector<vector<vector<TH1D*>>> vec_h_jet_pT_ME_10GeV;
    vector<vector<vector<TH1D*>>> vec_h_jet_pT_ME_5GeV;
    vector<vector<vector<TH1D*>>> vec_h_jet_pT_ME_4_5GeV;
    vector<vector<vector<TH1D*>>> vec_h_jet_pT_ME_4GeV;
    vector<vector<vector<TH1D*>>> vec_h_jet_pT_ME_3_5GeV;
    vector<vector<vector<TH1D*>>> vec_h_jet_pT_ME_3GeV;
    vector<vector<vector<TH1D*>>> vec_h_jet_pT_ME_2_5GeV;
    vector<vector<vector<TH1D*>>> vec_h_jet_pT_ME_2GeV;

    vector<vector<vector<TH1D*>>> vec_h_jet_pT_sub_ME_without_10GeV;
    vector<vector<vector<TH1D*>>> vec_h_jet_pT_sub_ME_without_5GeV;
    vector<vector<vector<TH1D*>>> vec_h_jet_pT_sub_ME_without_4_5GeV;
    vector<vector<vector<TH1D*>>> vec_h_jet_pT_sub_ME_without_4GeV;
    vector<vector<vector<TH1D*>>> vec_h_jet_pT_sub_ME_without_3_5GeV;
    vector<vector<vector<TH1D*>>> vec_h_jet_pT_sub_ME_without_3GeV;
    vector<vector<vector<TH1D*>>> vec_h_jet_pT_sub_ME_without_2_5GeV;
    vector<vector<vector<TH1D*>>> vec_h_jet_pT_sub_ME_without_2GeV;

    vector<vector<vector<TH1D*>>> vec_h_jet_pT_sub_ME_10GeV;
    vector<vector<vector<TH1D*>>> vec_h_jet_pT_sub_ME_4_5GeV;
    vector<vector<vector<TH1D*>>> vec_h_jet_pT_sub_ME_3_5GeV;

    vector<TH1D*> h_jet_track_pT_SE;
    vector<TH1D*> h_jet_track_pT_ME;

    vector<TH1D*> h_N_jet_SE;
    vector<TH1D*> h_N_jet_ME;

    vector<TH2D*> h2D_eta_phi_jet_SE;
    vector<TH2D*> h2D_eta_phi_jet_ME;

    // h-jet tests
    vector<vector<vector<TH1D*>>> vec_h1D_h_jet_sub_pT_SE;
    vector<vector<vector<TH1D*>>> vec_h1D_h_jet_sub_pT_ME;
    vector<vector<vector<vector<TH1D*>>>> vec_h1D_h_jet_sub_pT_SE_trigger_int;
    vector<vector<vector<vector<TH1D*>>>> vec_h1D_h_jet_sub_pT_ME_trigger_int;

    //-----------------------------
    vector<TH2D*> h2D_pT_sum_vs_mult_SE;
    vector<TH2D*> h2D_pT_sum_vs_mult_ME;
    vector<TH2D*> h2D_rho_vs_mult_SE;
    vector<TH2D*> h2D_rho_vs_mult_ME;
    vector<TH2D*> h2D_pT_sum_vs_rho_SE;
    vector<TH2D*> h2D_pT_sum_vs_rho_ME;
    vector<TH1D*> vec_h_pt_sum_SE;
    vector<TH1D*> vec_h_pt_sum_ME;
   
    TH1D* h_jet_pT_sub_profile_SE;
    TH1D* h_jet_pT_sub_profile_ME;
    TH1D* h_minus_rho_area_SE;
    TH1D* h_minus_rho_area_ME;
    
    TH2D* h2D_dNdphi_vs_phi_psi;

    vector<TH1D*> h_jet_const_pT_high_pT_jets;
    

    TH1D* h_mult_particles;
    TH1D* h_mult_particles_ITS;
    TH1D* h_mult_particles_ITS_1_2;
    TH1D* h_mult_particles_ITS_3_4;
    TH1D* h_mult_particles_ITS_5_6;
    TH1D* h_mult_particles_ITS_1_or_2;
    TH1D* h_mult_particles_ITS_3_or_4;
    TH1D* h_mult_particles_ITS_5_or_6;
    vector<TH2D*> h2D_jet_pT_vs_max_const_pT_ME;
    vector<TH2D*> h2D_jet_pT_vs_max_const_pT_SE;
    TH2D* h2D_N_const_vs_jet_pT_sub_ME;
    TH2D* h2D_N_const_vs_jet_pT_sub_SE;

    vector<vector<vector<TH1D*>>> vec_h_EP_bins_SE_input;

    TH1D* h_Psi_full_v3;
    TH2D* h2D_Psi_pos_vs_Psi_neg_v3;

    TH2D* h2D_qx_qy_for_EP_V0A;
    TH2D* h2D_qx_qy_for_EP_V0C;
    TH2D* h2D_qx_qy_for_EP_V0A_input;
    TH2D* h2D_qx_qy_for_EP_V0C_input;
    TH2D* h2D_qx_qy_for_EP_v3_V0A;
    TH2D* h2D_qx_qy_for_EP_v3_V0C;


    TH1D* h_phi_minus_Psi2_SE;
    TH1D* h_phi_minus_Psi2_all_SE;
    TH1D* h_cos_phi_minus_Psi2_SE;
    TH1D* h_phi_minus_Psi2_ME;
    TH1D* h_phi_minus_Psi2_all_ME;
    TH1D* h_cos_phi_minus_Psi2_ME;

    TH2D* h2D_v2_vs_pT_SE;
    TH2D* h2D_phi_minus_Psi2_vs_pT_SE;
    TH1D* h_v2_dist_8_25GeV_SE;
    TH1D* h_EbyE_v2_8_25_GeV_SE;
    TH1D* h_EbyE_v2_8_25_GeV_SE_V0;
    TH1D* h_EbyE_v2_8_25_GeV_SE_acc_pos;
    TH1D* h_EbyE_v2_8_25_GeV_SE_acc_neg;
    TH1D* h_EbyE_v2_8_25_GeV_SE_acc_pos_ss;
    TH1D* h_EbyE_v2_8_25_GeV_SE_acc_neg_ss;
    TH1D* h_EbyE_v2_st_3_GeV_SE;
    TH1D* h_EbyE_v2_3_5_GeV_SE;
    vector<TH1D*> h_EbyE_v2_full_pT_SE;
    TH2D* h2D_phi_vs_v2_ebye;
    TH1D* h_phi_dist_v2_lt_0_98;
    TH2D* h2D_phi_eta_v2_lt_0_98;

    TH2D* h2D_v2_vs_pT_ME;
    TH2D* h2D_phi_minus_Psi2_vs_pT_ME;
    TH1D* h_v2_dist_8_25GeV_ME;
    TH1D* h_EbyE_v2_8_25_GeV_ME;
    TH1D* h_EbyE_v2_3_5_GeV_ME;
    TH1D* h_EbyE_v2_st_3_GeV_ME;
    vector<TH1D*> h_EbyE_v2_full_pT_ME;
    Int_t flagEPPosNeg = 0;

    TF1* gauss_func;

    TH1D* h_V0A_amplitude;
    TH1D* h_V0C_amplitude;

    TH1D* h_Psi2_V0A;
    TH1D* h_Psi2_V0C;
    TH1D* h_Psi2_V0;
    TH2D* h2D_corr_EP_TPC_V0A;
    TH2D* h2D_corr_EP_TPC_V0C;
    TH2D* h2D_corr_EP_TPC_V0;

    vector<vector<TH1D*>> h_shiftCoeff;
    vector<vector<TH1D*>> h_shiftCoeff_V0A;
    vector<vector<TH1D*>> h_shiftCoeff_V0C;

    double dSin2psi;
    double dCos2psi;
    double dSin4psi;
    double dCos4psi;

    double dSin2psiA;
    double dCos2psiA;
    double dSin4psiA;
    double dCos4psiA;

    double dSin2psiC;
    double dCos2psiC;
    double dSin4psiC;
    double dCos4psiC;


    TH2D* h2D_particle_mult_n_tof_hits[3];
    TH2D* h2D_particle_mult_n_trd_tracklets[3];
    TH2D* h2D_particle_mult_fNumEMCal[3];
    TH2D* h2D_particle_mult_n_V0A[3];
    TH2D* h2D_particle_mult_n_V0C[3];
    TH2D* h2D_particle_mult_128_vs_tof;
    TH2D* h2D_TRD_tracklets_vs_n_tof_hits;

    TH1D* h_mult_particles_hybrid;
    TH1D* h_mult_particles_used_SE_ME[2];
    TH2D* h2D_mult_particles_hybrid_vs_mult;
    TH1D* h_input_mult_hist;
    TH2D* h2D_input_mult_particles_used_vs_mult;
    vector<TH1D*> vec_h_input_mult_hist;
    TH2D* h2D_QA_vs_QB_peak[2];
    TH2D* h2D_QA_vs_QB_tail[2];
    TH2D* h2D_mult_vs_Psi_diff;
    TH1D* h_dca_xy_peak;
    TH1D* h_dca_z_peak;
    TH1D* h_dca_xy_tail;
    TH1D* h_dca_z_tail;
    TH2D* h2D_Q_vec_peak;
    TH2D* h2D_Q_vec_tail;
    TH1D* h_EP_peak[2];

    TH2D* h2D_Psi_pos_vs_Psi_neg;
    TH2D* h2D_Psi_pos_vs_Psi_neg_ME;
    TH2D* h2D_Psi_V0A_vs_V0C;
    TH2D* h2D_jet_SE_vs_dPsi_phi;
    TH2D* h2D_jet_sub_SE_vs_dPsi_phi;
    TH2D* h2D_jet_sub_EP_SE_vs_dPsi_phi;
    TH2D* h2D_jet_ME_vs_dPsi_phi;
    TH2D* h2D_jet_sub_ME_vs_dPsi_phi;
    TH2D* h2D_jet_sub_EP_ME_vs_dPsi_phi;
    vector<TH2D*> vec_h2D_jet_SE_vs_dPsi_phi;
    TH1D* h_Psi_full;
    TH1D* h_Psi2_full;
    TH1D* h_Psi_full_ME;
    TH1D* h_Psi_full_ME_rejected;

    TH1D* h_Psi_V0A_V0C_EP_resolution;
    TH1D* h_Psi_SE_EP_resolution;
    TH1D* h_Psi_ME_EP_resolution;

    vector<TH2D*> vec_h2D_jet_sub_vs_mult_SE;
    vector<TH2D*> vec_h2D_jet_sub_vs_mult_ME;
    vector<TH2D*> vec_h2D_h_jet_sub_vs_mult_SE;
    vector<TH2D*> vec_h2D_h_jet_sub_vs_mult_ME;
    vector<TH2D*> vec_h2D_jet_rho_vs_mult_SE;
    vector<TH2D*> vec_h2D_jet_rho_vs_mult_ME;
    vector<vector<TH2D*>> vec_h2D_h_jet_sub_vs_mult_SE_diff_trigger_pT;
    vector<vector<TH2D*>> vec_h2D_h_jet_sub_vs_mult_ME_diff_trigger_pT;
    vector<vector<TH1D*>> vec_h_h_jet_sub_vs_mult_SE_diff_trigger_pT_project;
    vector<vector<TH1D*>> vec_h_h_jet_sub_vs_mult_ME_diff_trigger_pT_project;

    vector<vector<TH1D*>> vec_h_mult_SE;
    vector<vector<TH1D*>> vec_h_mult_ME;



    TH2D* h2D_dca_xy_vs_pT;
    TH2D* h2D_dca_z_vs_pT;
    TH2D* h2D_dca_xyz_vs_pT;
    TH2D* h2D_chi2_vs_pT;
    TH1D* h_particle_pT;
    TH1D* h_particle_pT_ME;
    TH1D* h_trigger_pT;
    vector<TH1D*> vec_h_trigger_pT;
    TH1D* h_input_track_mult;
    TH1D* h_input_phi_track;
    TH1D* h_phi_track_for_EP_pos_eta;
    TH1D* h_phi_track_for_EP_neg_eta;
    TH1D* h_input_phi_track_for_EP_pos_eta;
    TH1D* h_input_phi_track_for_EP_neg_eta;
    TH2D* h2D_qx_qy_for_EP_pos_eta;
    TH2D* h2D_qx_qy_for_EP_neg_eta;
    TH2D* h2D_qx_qy_for_EP_pos_eta_ME;
    TH2D* h2D_qx_qy_for_EP_neg_eta_ME;
    TH2D* h2D_qx_qy_for_EP_pos_eta_v3;
    TH2D* h2D_qx_qy_for_EP_neg_eta_v3;
    TH2D* h2D_input_qx_qy_for_EP_pos_eta;
    TH2D* h2D_input_qx_qy_for_EP_neg_eta;
    TH2D* h2D_input_qx_qy_for_EP_pos_eta_ME;
    TH2D* h2D_input_qx_qy_for_EP_neg_eta_ME;
    TH2D* h2D_input_qx_qy_for_EP_pos_eta_v3;
    TH2D* h2D_input_qx_qy_for_EP_neg_eta_v3;
    TH2D* h2D_eta_phi_ITS_3_hits_N_TOF_hits_left;
    TH2D* h2D_eta_phi_ITS_3_hits_N_TOF_hits_right;
    TH2D* h2D_particle_mult_n_tof_hits_left;
    TH2D* h2D_particle_mult_n_tof_hits_right;
    TH2D* h2D_z_vertex_vs_xy_vertex;
    TH2D* h2D_z_vertex_vs_xy_vertex_left;
    TH2D* h2D_z_vertex_vs_xy_vertex_right;
    TH2D* h2D_z_vertex_vs_particle_mult_3_ITS_hits;
    TH2D* h2D_z_vertex_vs_particle_mult_TPC_refit;
    TH1D* h_z_vertex_dist;
    TH2D* h2D_z_vertex_3ITS_hits_not_1[3];
    TH2D* h2D_z_vertex_3ITS_hits_not_2[3];
    TH2D* h2D_z_vertex_3ITS_hits_not_3[3];
    TH2D* h2D_z_vertex_3ITS_hits_not_4[3];
    TH2D* h2D_z_vertex_3ITS_hits_not_5[3];
    TH2D* h2D_z_vertex_3ITS_hits_not_6[3];

    vector<TH1D*> vec_h_phi_track;
    vector<TH1D*> vec_h_psi_test;
    vector<TH1D*> vec_h_psi3_test;
    vector<TH2D*> vec_h2D_phi_eta_track;
    vector<TH2D*> vec_h2D_phi_eta_track_pos;
    vector<TH2D*> vec_h2D_phi_eta_track_neg;
    vector<TH2D*> vec_h2D_hybrid_phi_eta_track;
    vector<TH1D*> vec_h_hybrid_phi_track;
    

    // centrality
    TH1D* h_centrality_V0M;
    TH1D* h_centrality_V0M_cuts;

    Int_t N_events_excluded = 0;

    TH2D* h2D_phi_track_bits;
    TH2D* h2D_z_vertex_vs_particle_mult_hybrid;

    TH1D* h_z_dca;
    TH1D* h_xy_dca;
    TH1D* h_track_pT;
    TH1D* h_TPC_chi2;
    TH1D* h_z_dca_eta_0_84;
    TH1D* h_xy_dca_eta_0_84;
    TH1D* h_TPC_chi2_eta_0_84;

    TH1D* vec_h_phi_track_hybrid;
    TH1D* vec_h_phi_track_512;
    TH1D* vec_h_phi_track_256;           
    TH2D* vec_h2D_phi_eta_track_hybrid;
    TH2D* vec_h2D_phi_eta_track_512;
    TH2D* vec_h2D_phi_eta_track_256;

    TRandom random_e;
    TProfile* PX_jet_SE_vs_dPsi_phi[10];
    TProfile* PX_jet_ME_vs_dPsi_phi[10];

    TH1D* h_parameter_event_mult[2];
    TH1D* h_parameter_event_prim_z_vertex[2];   
    TH1D* h_parameter_event_TOF_hits[2];        
    TH1D* h_parameter_event_pile_up[2];         
    TH1D* h_parameter_event_n_prim[2];          
    TH1D* h_parameter_event_N_TRD_tracklets[2];
    vector<Double_t> vec_phi_track;
    Int_t run_ID;
    TH1D* h_run_ID_jet_pT_minus_35GeV_SE;
    TH1D* h_run_ID;
    TH1D* h_phi_track_jet_pT_minus_35GeV_SE;
    TH1D* h_phi_const_jet_pT_minus_35GeV_SE;
    TH1D* h_phi_const_jet_pT_minus_35GeV_SE_eta_pos;
    TH1D* h_phi_const_jet_pT_minus_35GeV_SE_eta_neg;
    TH1D* h_phi_track_SE;
    TH1D* h_phi_track_ME;
    TH2D* h2D_phi_track_eta_ME;
    vector<TH1D*> h_jet_sub_SE_eta_neg;
    vector<TH1D*> h_jet_sub_SE_eta_pos;
    vector<TH1D*> h_jet_sub_ME_eta_neg;
    vector<TH1D*> h_jet_sub_ME_eta_pos;

    TH1D* h_run_ID_low_rho;
    TH1D* h_mult_low_rho;
    TH1D* h_TOF_hits_low_rho;
    TH1D* h_TRD_tracklets_low_rho;
    TH1D* h_run_ID_norm_rho;
    TH1D* h_mult_norm_rho;
    TH1D* h_TOF_hits_norm_rho;
    TH1D* h_TRD_tracklets_norm_rho;

    TGraph* graph_cut_pileupTOF;
    TGraph* graph_cut_pileupTRD;
    TGraph* graph_cut_pileupV0A;
    TGraph* graph_cut_pileupV0C;

    //------------------
    // Different radius parameters 29.12
    vector<Double_t> vec_jet_radii;

    //----------------------------------------------------------------------
    // Parameters for track reconstruction efficiency
    // Centrality: 0-5,5-10,10-20,20-30,30-40,40-50,50-60,60-70,70-80
    // runID: Index 0 is for all run ranges, then starting from 1: 12138025,12145021,12152017,12154022,12165032,12171017
    // PID: p, anti-p", pi+, pi-, K+, K-
    Double_t Global_ABC_Parameters_A[9][7][6];  //centrality(9),runID(7),PID(6)
    Double_t Global_ABC_Parameters_B[9][7][6];  //centrality(9),runID(7),PID(6)
    Double_t Global_ABC_Parameters_C[9][7][6];  //centrality(9),runID(7),PID(6)
    TString eEff_file;
    //----------------------------------------------------------------------


    //----------------------------------------------------------------------
    vector< vector<TH1D*> > vec_TH1D_TRD_geometry; // store for all 540 chambers the 8 corner vertices per detector
    TEveLine* TEveLine_beam_axis = NULL;
    TEvePointSet* TEveP_primary_vertex;
    TEveLine* TPL3D_helix = NULL;
    vector<TEveLine*> vec_TPL3D_helix;
    vector<TEveLine*> vec_TPL3D_helix_inner;
    vector<TEveLine*> vec_TPL3D_helix_hull;
    vector<TEveBox*> vec_eve_TRD_detector_box;
    vector<TEveBox*> vec_eve_EMCal_cluster;
    vector<TEveBox*> vec_eve_Jets;
    vector<TEveBox*> vec_eve_tracks_jets;
    //----------------------------------------------------------------------


    ClassDef(StJetAnalysis,1)

};

#endif
