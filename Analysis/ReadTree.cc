#include "./ReadTree.h"
#include <Math/Vector4D.h>
#include "functions.h"
#include "StJetTrackEvent.h"
#include "StJetTrackEventLinkDef.h"

void ReadTree(TString InputFile, TString OutputFile)
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

    //Get data from tree file
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
                TString treePathCol   = Form("Trees.root/%s/O2tablecol", dirName.Data());
                TString treePathTrack = Form("Trees.root/%s/O2tabletrack", dirName.Data());
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


    // Don't open all trees but only open individual DFs and read trees from there
    Long64_t gi;
    Long64_t track_CollisionID;
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
    //trackChain->SetBranchAddress("fLength", &length);

    //get the first particle indices of each event and save it into FirstParticleIndex
    vector<uint32_t> FirstParticleIndex;//holds the first particle index of the first event, second event, ... in positions (0, 1, 2, ...)
    FirstParticleIndex.push_back(0);
    trackChain->GetEntry(0);
    uint64_t LastColID = colID_track;

    for(uint64_t i_Track = 1; i_Track < entries_track; i_Track ++)
    {
        trackChain->GetEntry(i_Track);
        if(colID_track != LastColID)
        {
            LastColID = colID_track;
            FirstParticleIndex.push_back(i_Track);
        }
    }
    FirstParticleIndex.push_back(entries_track);

    cout << "Have identified First Particles" << endl;

    //loop over events
    for(uint16_t i_Event = 0; i_Event < entries_col; i_Event ++)
    {
        //process global event information
        colChain -> GetEntry(i_Event);
        //cout << "-----------------------------------------------------------" << endl;
        cout << "------ Collision number " << i_Event << " of run " << rn << " with multiplicity " << mult << ", centrality " << cent << ", z vertex " << vertz << endl;

        cout << "Collision ID event " << colID_event << endl;

        h1d_mult_event  ->Fill(mult);
        h1d_cent        ->Fill(cent);
        h1d_psi2        ->Fill((float)psi2/1000);
        h1d_occft0_event->Fill(occuft0);
        h1d_occ_event   ->Fill(occu);
        h1d_mult_vs_cent->Fill(mult,cent);

        cout << "number of tracks in this event is " << (FirstParticleIndex[i_Event+1] - FirstParticleIndex[i_Event]) << endl;

        // INSERT EVENT ANALYIS HERE


        //loop over tracks in event i_Event
        for(uint16_t i_Track = FirstParticleIndex[i_Event-1]; i_Track < FirstParticleIndex[i_Event]; i_Track ++)
        {
            //obtain particle track quantities
            trackChain->GetEntry(i_Track);
            Float_t px = get_px_Particle(p);
            Float_t py = get_py_Particle(p);
            Float_t pz = get_pz_Particle(p);
            Float_t m = 0.1349766;
            ROOT::Math::PxPyPzMVector vec_track(px,py,pz,m);

            //INSERT PARTICLE ANALYSIS HERE

            if(vec_track.Pt() > 20.0) cout << "High track " << vec_track.Pt() << " GeV/c" << endl;

            h1d_pT ->Fill(vec_track.Pt());
            h1d_eta->Fill(vec_track.Eta());
            h1d_phi->Fill(vec_track.Phi());

        }

    }


    //--------- Write -----------------------------------
    cout << "Write" << endl;
    TFile* output_file = new TFile(OutputFile,"RECREATE");    // Save output without scaling
    output_file->cd();

    // Track properties
    h1d_pT ->Write();
    h1d_eta->Write();
    h1d_phi->Write();
    //h1d_TOF->Write();

    // Event properties
    h1d_mult_event  ->Write();
    h1d_pTsum       ->Write();
    h1d_cent        ->Write();
    h1d_psi2        ->Write();
    h1d_occft0_event->Write();
    h1d_occ_event   ->Write();
    h1d_mult_vs_cent->Write();

    output_file->Close();

    cout << "ReadTree end" << endl;
}