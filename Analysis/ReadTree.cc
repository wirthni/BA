#include "./ReadTree.h"
#include <Math/Vector4D.h>
#include <unordered_map>
#include <vector>
#include "functions.h"
#include "StJetTrackEvent.h"
#include "StJetTrackEventLinkDef.h"

Int_t ReadTree(TString InputFile, TString OutputFile)
{
    cout << "ReadTree started" << endl;

    // Make histograms for a global event analysis
    TH1D* h1d_mult_event   = new TH1D("h1d_mult_event","h1d_mult_event", 500,0,5000);
    TH1D* h1d_occ_event    = new TH1D("h1d_occ_event","h1d_occ_event", 10000,0,10000);
    TH1D* h1d_occft0_event = new TH1D("h1d_occft0_event","h1d_occft0_event", 100000,0,100000);
    TH1D* h1d_pTsum        = new TH1D("h1d_pTsum","h1d_pTsum", 2000,0,2000);
    TH1D* h1d_TOF          = new TH1D("h1d_TOF","h1d_TOF", 60000,-3000000,3000000);
    TH1D* h1d_cent         = new TH1D("h1d_cent","h1d_cent", 80,0,40);
    TH1D* h1d_psi2         = new TH1D("h1d_psi2","h1d_psi2", 800,-4.0,4.0);
    TH1D* h1d_pT           = new TH1D("h1d_pT","h1d_pT", 1000,0,200);
    TH1D* h1d_eta          = new TH1D("h1d_eta","h1d_eta", 200,-1.0,1.0);
    TH1D* h1d_phi          = new TH1D("h1d_phi","h1d_phi", 800,-4.0,4.0);

    TH2D* h1d_mult_vs_cent = new TH2D("h1d_mult_vs_cent","h1d_mult_vs_cent", 500,0,5000,80,0,40);

    TFile* file = TFile::Open(InputFile);

    TChain* colChain   = new TChain("O2tablecol");
    TChain* trackChain = new TChain("O2tabletrack");

    TIter next(file->GetListOfKeys());
    TKey *key;

    while ((key = (TKey*)next()))
    {

        if (TString(key->GetClassName()) == "TDirectoryFile")
        {
            TString dirName = key->GetName();

            // Optional: Nur Ordner die mit "DF_" beginnen
            if (dirName.BeginsWith("DF_"))
            {
                //TString treePathCol   = Form("Trees.root/%s/O2tablecol", dirName.Data());
                //TString treePathTrack = Form("Trees.root/%s/O2tabletrack", dirName.Data());
                TString treePathCol   = Form(InputFile + "/%s/O2tablecol", dirName.Data());
                TString treePathTrack = Form(InputFile + "/%s/O2tabletrack", dirName.Data());
                colChain  ->Add(treePathCol);
                trackChain->Add(treePathTrack);
                cout << "Added " << treePathCol << endl;
            }
        }
    }

    Int_t entries_col   =    colChain->GetEntries();
    Int_t entries_track =  trackChain->GetEntries();


    cout << "------------------------ COLLISION TABLE -----------------------------------" << endl;
    cout <<  entries_col << " entries " << endl;
    cout << "------------------------ TRACK TABLE -----------------------------------" << endl;
    cout <<  entries_track << " entries " << endl;
    cout << "-----------------------------------------------------------" << endl;

    std::unordered_map<Long64_t, std::vector<Long64_t>> CollisionMap; //holds the Collision ID and underneath the track indices inside of the whole trackChain (0,1,2,...,2644) for example

    Long64_t gi;
    Int_t track_CollisionID;
    int32_t trackid;
    int rn, mult, occu;
    float vertx, verty, vertz, cent, occuft0, tof;
    short psi2, psi3, charge, dcaxy, dcaz;
    unsigned short dedx,length;
    ULong64_t p;

    //colChain->SetBranchAddress("fGI",&gi);
    //colChain->SetBranchAddress("fCi", &colID_event);
    colChain->SetBranchAddress("fRn",&rn);
    colChain->SetBranchAddress("fCent",&cent);
    colChain->SetBranchAddress("fMult",&mult);
    colChain->SetBranchAddress("fOccu",&occu);
    colChain->SetBranchAddress("fOccuft0",&occuft0);
    colChain->SetBranchAddress("fVertexX",&vertx);
    colChain->SetBranchAddress("fVertexY",&verty);
    colChain->SetBranchAddress("fVertexZ",&vertz);
    colChain->SetBranchAddress("fPsi2",&psi2);
    colChain->SetBranchAddress("fPsi3",&psi3);

    trackChain->SetBranchAddress("fIndexTableCols", &track_CollisionID);
    trackChain->SetBranchAddress("fCharge", &charge);
    trackChain->SetBranchAddress("fP", &p);
    trackChain->SetBranchAddress("fDedx", &dedx);
    trackChain->SetBranchAddress("fDcaxy", &dcaxy);
    trackChain->SetBranchAddress("fDcaz", &dcaz);

    Long64_t Counter = 0;
    // Make an unordered map particle->collision
    for(uint64_t i_Track = 0; i_Track < entries_track; i_Track ++) {

        trackChain->GetEntry(i_Track);
        CollisionMap[track_CollisionID].push_back(i_Track);
        if(track_CollisionID == 40) Counter++;
    }

    cout << Counter << " tracks have ID 0!" << endl;

    //LOOP
    for (int j = 0; j < 120; j++) {

        colChain->GetEntry(j);
        int collisionID = j;

        cout << "centrality of event " << j << " is " << cent << " and number of particles is " << CollisionMap[j].size() << " and mult is " << mult << endl;
    

    }

    return 1;
}