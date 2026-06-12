#include "./StJetTrackEvent.h"
#include "./StJetTrackEventLinkDef.h"
#include "./functions.h"
#include "TClassRef.h"
#include "Math/Vector4D.h"
#include <algorithm>

using namespace std;

/*Parameters:
i_RootFileName: "ROOTFILE_x.root" -> ROOTFILE
i_PathToHistosInRootFile: -> subfolder/
i_Histoname: ->h2D_...
i_NoOfHistos: clear
*/
Int_t HistoMerger(TString i_RootFileName, TString i_PathToHistosInRootFile, TString i_HistoName, int i_NoOfHistos, bool i_ScaleByNumber)
{

    gStyle->SetOptStat(0);
    SetRootGraphicStyle();

    TString str_out = "./Merger_" + i_RootFileName + "_" + i_HistoName + ".root";
    TFile* OutputFile = new TFile(str_out.Data(),"RECREATE");
    OutputFile->cd();
    cout << "Output file is " << str_out << endl;

    //get first one to have correct type
    TString DataFile = "./DataToMerge/" + i_RootFileName + Form("_%d", 1) + ".root";
    TFile *file1 = TFile::Open(DataFile);
    if (!file1 || file1->IsZombie()) {
        std::cout << "Error while opening file " << 1 << "!" << endl;
        return 0;
    }

    TH2D* h2D_MergerResult = (TH2D*)file1->Get(i_PathToHistosInRootFile + i_HistoName);

    for(int i=2; i<=i_NoOfHistos; i++)
    {
        DataFile = "./DataToMerge/" + i_RootFileName + Form("_%d", i) + ".root";

        TFile *file = TFile::Open(DataFile);

        if (!file || file->IsZombie()) {
            std::cout << "Error while opening file " << i << "!" << endl;
            return 0;
        }

        TH2D* h2D_MergerScrawl = (TH2D*)file->Get(i_PathToHistosInRootFile + i_HistoName);

        h2D_MergerResult->Add(h2D_MergerScrawl);

    }

    if(i_ScaleByNumber) h2D_MergerResult->Scale(1./double(i_NoOfHistos));

    OutputFile->cd();
    h2D_MergerResult->Write();

    OutputFile->Close();

    return 1;
}