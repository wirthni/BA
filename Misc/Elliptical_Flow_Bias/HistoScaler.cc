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

int Scale(TString Filename, TString HistoName, double Factor)
{
    TFile *file = TFile::Open(Filename);

    if (!file || file->IsZombie()) {
        std::cout << "Error while opening the file!" << std::endl;
        return 0;
    }

    //get statistics histograms from root file
    TH2D* Histo = (TH2D*)file->Get(HistoName);
    file->Close();
    if(!Histo){cout << "Histogram not found!" << endl; return 0;}
    

    Histo->Scale(Factor);

    TString str_out = Filename + ".root";
    TFile* OutputFile = new TFile(str_out.Data(),"RECREATE");
    OutputFile->cd();
    Histo->Write();

    return 1;
}