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
Int_t HistoMerger(TString i_RootFileName, TString i_PathToHistosInRootFile, int i_NoOfHistos, bool i_ScaleByNumber)
{

    gStyle->SetOptStat(0);
    SetRootGraphicStyle();

    //other analysis variables
    vector<TString> HistoVector;

    //name output file
    TString str_out = "./M_" + i_RootFileName + ".root";
    TFile* OutputFile = new TFile(str_out.Data(),"RECREATE");
    cout << "Output file is " << str_out << endl;

    //open first root file
    TString DataFile = "./DataToMerge/" + i_RootFileName + Form("_%d", 1) + ".root";
    TFile *file1 = TFile::Open(DataFile);
    if (!file1 || file1->IsZombie()) {
        std::cout << "Error while opening file " << 1 << "!" << endl;
        return 0;
    }
    cout << "Getting list of histograms from first datafile " << DataFile << endl;

    TDirectory *Directory = file1->GetDirectory(i_PathToHistosInRootFile);
    TIter next(Directory->GetListOfKeys());
    TKey *key;
    while ((key = (TKey *)next())) {
        HistoVector.push_back(key->GetName());
    }

    cout << "Found histograms" << endl;

    for(int i=0; i<HistoVector.size(); i++)
    {
        cout << HistoVector[i] << endl;
    }
    cout << endl;

    for(int i=0; i<HistoVector.size(); i++)
    {
        cout << "Now merging " << HistoVector[i] << endl;

        //take general histo info from first file
        TH2D* h2D_MergerResult = (TH2D*)Directory->Get(HistoVector[i]);

        //append the other files 2,3,...
        for(int j=2; j<=i_NoOfHistos; j++)
        {
            DataFile = "./DataToMerge/" + i_RootFileName + Form("_%d", j) + ".root";

            TFile *file = TFile::Open(DataFile);

            if (!file || file->IsZombie()) {
                std::cout << "Error while opening file " << j << "!" << endl;
                return 0;
            }

            TH2D* h2D_MergerScrawl = (TH2D*)file->Get(i_PathToHistosInRootFile + "/" + HistoVector[i]);

            h2D_MergerResult->Add(h2D_MergerScrawl);

            file->Close();

        }

        if(i_ScaleByNumber) h2D_MergerResult->Scale(1./double(i_NoOfHistos));

        OutputFile->cd();
        h2D_MergerResult->Write();

    }

    OutputFile->Close();

    return 1;

}