#include <TString.h>
#include <TFile.h>
#include <TKey.h>
#include <TList.h>
#include <TCollection.h>
#include <TSystem.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <array>

// -----------------------------------------------------------------------------
// DijetAna_ClusterJob
//
// i_InputFileList ist eine .txt Datei, in der pro Zeile der Pfad zu einer
// .root Datei steht (z.B. "filelist.txt"). Jede dieser .root Dateien enthaelt
// (top-level) Directories, die nach dem Schema "DF_xxx" benannt sind.
// Fuer jedes gefundene DF_xxx Directory wird ein Eintrag in einen
// std::vector<std::array<TString, 2>> geschrieben:
//   entry[0] = Dateiname (voller Pfad) der .root Datei
//   entry[1] = Name des DF_xxx Directories
// (Hinweis: rohe Arrays wie TString[2] sind nicht kopier-/zuweisbar und
//  daher als Vector-Elementtyp ungueltig - std::array loest das sauber.)
// -----------------------------------------------------------------------------

// Sammelt alle DF_xxx Directorynamen aus einer einzelnen root Datei
// und haengt (Dateiname, DF-Name) Paare an io_DFList an.
void CollectDFsFromFile(const TString &i_FileName, std::vector<std::array<TString, 2>> &io_DFList)
{
    TFile *file = TFile::Open(i_FileName, "READ");
    if (!file || file->IsZombie())
    {
        std::cerr << "Warnung: Konnte Datei nicht oeffnen: " << i_FileName << std::endl;
        if (file) delete file;
        return;
    }

    TList *keyList = file->GetListOfKeys();
    if (!keyList)
    {
        std::cerr << "Warnung: Keine Keys in Datei: " << i_FileName << std::endl;
        file->Close();
        delete file;
        return;
    }

    TIter next(keyList);
    TKey *key = nullptr;
    while ((key = (TKey*)next()))
    {
        TString className = key->GetClassName();
        TString objName    = key->GetName();

        // Nur Directories (TDirectory / TDirectoryFile) beruecksichtigen
        if (className == "TDirectoryFile" || className == "TDirectory")
        {
            if (objName.BeginsWith("DF_"))
            {
                std::array<TString, 2> entry;
                entry[0] = i_FileName;
                entry[1] = objName;
                io_DFList.push_back(entry);
            }
        }
    }

    file->Close();
    delete file;
}

// Liest die Filelist (eine .txt Datei, ein .root Dateiname pro Zeile) ein
// und gibt den Vektor aller (Dateiname, DF-Name) Paare zurueck.
std::vector<std::array<TString, 2>> BuildDFList(const TString &i_InputFileListDirectory,
                                                  const TString &i_InputFileList,
                                                  Int_t i_NoOfDFs)
{
    std::vector<std::array<TString, 2>> dfList;

    TString listPath = i_InputFileListDirectory + "/" + i_InputFileList;

    if (!listPath.EndsWith(".txt"))
    {
        std::cerr << "Hinweis: '" << listPath
                   << "' hat keine .txt Endung - wird trotzdem als Textdatei gelesen."
                   << std::endl;
    }

    std::ifstream inFile(listPath.Data());

    if (!inFile.is_open())
    {
        std::cerr << "Fehler: Konnte Filelist nicht oeffnen: " << listPath << std::endl;
        return dfList;
    }

    std::string line;
    while (std::getline(inFile, line))
    {
        if (line.empty()) continue;

        TString fileName(line.c_str());
        fileName = fileName.Strip(TString::kBoth);
        if (fileName.IsNull()) continue;
        if (fileName.BeginsWith("#")) continue; // Kommentarzeilen ueberspringen

        CollectDFsFromFile(i_InputFileListDirectory + fileName, dfList);

        // Falls eine maximale Anzahl an DFs vorgegeben ist, hier abbrechen
        if (i_NoOfDFs > 0 && (Int_t)dfList.size() >= i_NoOfDFs)
        {
            dfList.resize(i_NoOfDFs);
            break;
        }
    }

    inFile.close();
    return dfList;
}

Int_t DijetAna_ClusterJob(TString i_InputFileListDirectory = "./",
                           TString i_InputFileList = "filelist.txt",
                           TString i_OutputFile = "DijetAna",
                           Int_t i_NoOfDFs = -1,
                           bool i_UseHadronInstead = false,
                           int i_NoOfNodes = 1)
{
    std::cout << "=== DijetAna_ClusterJob gestartet ===" << std::endl;
    std::cout << "InputFileListDirectory : " << i_InputFileListDirectory << std::endl;
    std::cout << "InputFileList          : " << i_InputFileList          << std::endl;
    std::cout << "OutputFile             : " << i_OutputFile             << std::endl;
    std::cout << "NoOfDFs                : " << i_NoOfDFs                << std::endl;
    std::cout << "UseHadronInstead       : " << i_UseHadronInstead       << std::endl;
    std::cout << "NoOfNodes              : " << i_NoOfNodes              << std::endl;

    std::vector<std::array<TString, 2>> dfList = BuildDFList(i_InputFileListDirectory,
                                                               i_InputFileList,
                                                               i_NoOfDFs);

    std::cout << "Gefundene DF-Eintraege: " << dfList.size() << std::endl;
    for (size_t i = 0; i < dfList.size(); ++i)
    {
        std::cout << "  [" << i << "] Datei: " << dfList[i][0]
                   << "  DF: " << dfList[i][1] << std::endl;
    }

    if (dfList.empty())
    {
        std::cerr << "Fehler: Keine DF_xxx Directories gefunden." << std::endl;
        return 1;
    }

    // Ab hier: dfList[i][0] (Dateiname) / dfList[i][1] (DF-Name)
    // fuer die weitere Analyse verwenden, z.B. Aufteilung auf i_NoOfNodes
    // Knoten oder Verarbeitung von Hadron- statt Parton-Level (i_UseHadronInstead).

    return 0;
}