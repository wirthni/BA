#include <Math/Vector4D.h>
#include <unordered_map>
#include <vector>
#include "functions.h"
#include "StJetTrackEvent.h"
#include "StJetTrackEventLinkDef.h"
#include "TVector3.h"

// Histos that are accessed by different methods and have to be available globally
TH2D* h2D_qx_qy_for_EP_pos_eta;
TH2D* h2D_qx_qy_for_EP_neg_eta;
TH2D* h2D_Psi_pos_vs_Psi_neg;
TH1D* vec_h_psi_test[2];
TH1D* h_Psi_SE_EP_resolution;
TH1D* h_EP_peak[2];
TH2D* h2D_QA_vs_QB_peak[2];
TH2D* h2D_QA_vs_QB_tail[2];
TH2D* h2D_Q_vec_peak;
TH2D* h2D_Q_vec_tail;

//Variables that are accessed by different methods and have to be available globally
Double_t Qvec_correction_qx[2] = {0.0};
Double_t Qvec_correction_qy[2] = {0.0};
bool FillQVectorCorrectionHistos_QVectorCorrectionRetrieved;
bool GetEventPlaneInfo_QVectorCorrectionRetrieved;

Int_t ReadTreeInit();
void PrintProgress(int iEvent, int i_Total);
void PrintInfo(TString i_String);
float FillQVectorCorrectionHistos(int i_CollisionID, TTree* i_TrackTree, std::unordered_map<int, std::vector<int>> i_CollisionMap);
float GetEventPlaneInfo(int i_CollisionID, TTree* i_TrackTree, std::unordered_map<int, std::vector<int>> i_CollisionMap, TString i_QVectorCorrectionFile);
Float_t get_px_Particle(Long64_t Particle_p);
Float_t get_py_Particle(Long64_t Particle_p);
Float_t get_pz_Particle(Long64_t Particle_p);

class Collision {
public:
    int RunNumber, ColID, Multiplicity, Occupancy;
    float VertX, VertY, VertZ, Centrality, Psi2, Psi3, OccupancyFT0;

    static Collision Read(TTree* tree, Long64_t entry) {

        //fetched variables
        static int RunNumber, ColID, Multiplicity, Occupancy;
        static float VertX, VertY, VertZ, Centrality, OccupancyFT0;
        static Short_t Psi2_Substitute, Psi3_Substitute;
        //derived variables
        static float Psi2, Psi3;
        //other variables
        static TTree* currentTree = nullptr;

        if (tree != currentTree) {
            tree->SetBranchAddress("fRn",&RunNumber);
            tree->SetBranchAddress("fCent",&Centrality);
            tree->SetBranchAddress("fMult",&Multiplicity);
            tree->SetBranchAddress("fOccu",&Occupancy);
            tree->SetBranchAddress("fOccuft0",&OccupancyFT0);
            tree->SetBranchAddress("fVertexX",&VertX);
            tree->SetBranchAddress("fVertexY",&VertY);
            tree->SetBranchAddress("fVertexZ",&VertZ);
            tree->SetBranchAddress("fPsi2",&Psi2_Substitute);
            tree->SetBranchAddress("fPsi3",&Psi3_Substitute);
            currentTree = tree;
        }
        tree->GetEntry(entry);

        ColID = entry;
        Psi2 = static_cast<float>(Psi2_Substitute)/1000.0f;
        Psi3 = static_cast<float>(Psi3_Substitute)/1000.0f;


        return Collision{ColID, RunNumber, Centrality, Multiplicity, Occupancy, OccupancyFT0, VertX, VertY, VertZ, Psi2, Psi3};
    }

private:
    Collision(int ColID_, int RunNumber_, float Centrality_, int Multiplicity_, int Occupancy_, float OccupancyFT0_, float VertX_, float VertY_, float VertZ_, float Psi2_, float Psi3_)
    : ColID(ColID_), RunNumber(RunNumber_), Centrality(Centrality_), Multiplicity(Multiplicity_), Occupancy(Occupancy_), VertX(VertX_), VertY(VertY_), VertZ(VertZ_), Psi2(Psi2_), Psi3(Psi3_) {}
};

class Particle {
public:
    int ColID, Charge;
    float px, py, pz, Pt, E, dEdX, DCAxy, DCAz, Eta, Phi;

    static Particle Read(TTree* tree, Long64_t entry) {

        //fetched variables
        static int ColID;
        static Short_t Charge;
        static ULong64_t P_Substitute;
        static unsigned short dEdX_Substitute;
        static short DCAxy_Substitute, DCAz_Substitute;
        //derived variables
        static float px, py, pz, E, dEdX, DCAxy, DCAz, Eta, Phi, Pt;
        //other variables
        static TTree* currentTree = nullptr;

        if (tree != currentTree) {
            tree->SetBranchAddress("fIndexTableCols", &ColID);
            tree->SetBranchAddress("fCharge", &Charge);
            tree->SetBranchAddress("fP", &P_Substitute);
            tree->SetBranchAddress("fDedx", &dEdX_Substitute);
            tree->SetBranchAddress("fDcaxy", &DCAxy_Substitute);
            tree->SetBranchAddress("fDcaz", &DCAz_Substitute);
            currentTree = tree;
        }
        tree->GetEntry(entry);

        px = get_px_Particle(P_Substitute);
        py = get_py_Particle(P_Substitute);
        pz = get_pz_Particle(P_Substitute);

        dEdX = static_cast<float>(dEdX_Substitute)/10.0f;
        DCAxy = static_cast<float>(DCAxy_Substitute)/100.0f;
        DCAz = static_cast<float>(DCAz_Substitute)/100.0f;
        E = sqrt(pow(px, 2) + pow(py, 2) + pow(pz, 2) + pow(0.1349766, 2));
        ROOT::Math::PxPyPzEVector ParticleLV(px, py, pz, E);
        Eta = ParticleLV.Eta();
        Phi = ParticleLV.Phi();
        Pt = sqrt(pow(px, 2) + pow(py, 2));


        return Particle{px, py, pz, Pt, E, dEdX, DCAxy, DCAz, Eta, Phi, ColID, Charge};
    }

private:
    Particle(float px_, float py_, float pz_, float Pt_, float E_, float dEdX_, float DCAxy_, float DCAz_, float Eta_, float Phi_, int ColID_, Short_t Charge_)
        : px(px_), py(py_), pz(pz_), Pt(Pt_), E(E_), dEdX(dEdX_), DCAxy(DCAxy_), DCAz(DCAz_), Eta(Eta_), Phi(Phi_), ColID(ColID_), Charge(Charge_) {}
};


Float_t get_px_Particle(Long64_t Particle_p)
{
    Long64_t Particle_px = 0;
    for(Int_t i_bit = 0; i_bit < 20 ;i_bit++)
    {
        if((Particle_p & ((ULong64_t) 1 <<  i_bit))) Particle_px |= (Long64_t)1 << i_bit;
    }
    if((Particle_p & ((ULong64_t) 1 <<  20)))  Particle_px = (-1)*Particle_px;
    return ((Float_t)Particle_px)/6000.0;
}

Float_t get_py_Particle(Long64_t Particle_p)
{
    Long64_t Particle_py = 0;
    for(Int_t i_bit = 0; i_bit < 20 ;i_bit++)
    {
        if((Particle_p & ((ULong64_t) 1 <<  (i_bit+21)))) Particle_py |= (Long64_t)1 << i_bit;
    }
    if((Particle_p & ((ULong64_t) 1 <<  41)))  Particle_py = (-1)*Particle_py;
    return ((Float_t)Particle_py)/6000.0;
}

Float_t get_pz_Particle(Long64_t Particle_p)
{
    Long64_t Particle_pz = 0;
    for(Int_t i_bit = 0; i_bit < 20 ;i_bit++)
    {
        if((Particle_p & ((ULong64_t) 1 <<  (i_bit+42)))) Particle_pz |= (Long64_t)1 << i_bit;
    }
    if((Particle_p & ((ULong64_t) 1 <<  62)))  Particle_pz = (-1)*Particle_pz;
    return ((Float_t)Particle_pz)/6000.0;
}