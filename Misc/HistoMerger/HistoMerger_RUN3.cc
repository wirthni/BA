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
Int_t HistoMerger(int i_FromRunNumber, int i_ToRunNumber, TString i_Suffix)
{

    gStyle->SetOptStat(0);
    SetRootGraphicStyle();

    //other analysis variables
    vector<TString> HistoVector;
    int NoOfCombinedFiles = 1;

    //name output file
    TString str_out = "./M_" + i_Suffix + ".root";
    TFile* OutputFile = new TFile(str_out.Data(),"RECREATE");
    cout << "Output file is " << str_out << endl;

    //open first root file
    TString FirstRunNumberAsString;
    FirstRunNumberAsString.Form("%d", i_FromRunNumber);
    TString DataFile = "DataToMerge/" + (TString)FirstRunNumberAsString.Data() + i_Suffix + ".root";
    TFile *file1 = TFile::Open(DataFile);
    if (!file1 || file1->IsZombie()) {
        std::cout << "Error while opening file " << 1 << "!" << endl;
        return 0;
    }
    cout << "Getting list of histograms from first datafile " << DataFile << endl;

    TIter next(file1->GetListOfKeys());
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
        int CorruptFilesCounter = 0;
        cout << "Now merging " << HistoVector[i] << endl;

        //take general histo info from first file
        TH3F* h2D_MergerResult = (TH3F*)file1->Get(HistoVector[i]);
        if(!h2D_MergerResult || h2D_MergerResult->GetEntries() == 0 || h2D_MergerResult->IsZombie() || h2D_MergerResult->GetEntries() != h2D_MergerResult->GetEntries())
        {
             cout << "First file has corrupt data, cannot retrieve histogram!" << endl;
             continue;
        }

        //append the other files 2,3,...
        for(int j = i_FromRunNumber+1; j <= i_ToRunNumber; j++)
        {

            TString RunNumberAsString;
            RunNumberAsString.Form("%d", j);
            DataFile = "DataToMerge/" + (TString)RunNumberAsString.Data() + i_Suffix + ".root";

            if(gSystem->AccessPathName(DataFile))
            {
                continue;
            }
            
            TFile *file = TFile::Open(DataFile);
            if (!file || file->IsZombie())
            {
                cout << "Skipping corrupt file!" << endl;
                delete(file);
                continue;
            }
            

            TH3F* h2D_MergerScrawl = (TH3F*)file->Get( HistoVector[i]);

            if(!h2D_MergerScrawl || h2D_MergerScrawl->GetEntries() == 0 || h2D_MergerScrawl->IsZombie() || h2D_MergerScrawl->GetEntries() != h2D_MergerScrawl->GetEntries() || h2D_MergerScrawl->GetMean(1) != h2D_MergerScrawl->GetMean(1) )
            {
                CorruptFilesCounter++;
                cout << "Skipping corrupt or empty histo in file " << j << endl;
                file->Close();
                continue;
            }

            h2D_MergerResult->Add(h2D_MergerScrawl);
            NoOfCombinedFiles++;

            file->Close();
            delete(file);

        }

        h2D_MergerResult->Scale(1./(double)(NoOfCombinedFiles));

        OutputFile->cd();
        h2D_MergerResult->Write();

    }

    OutputFile->Close();

    return 1;

}