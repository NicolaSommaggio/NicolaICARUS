#include "slice_struct.h"
void openFile()
{
    TFile * infile = TFile::Open("varMC_STANDARD.root","READ");
    infile -> cd();
        
    TTree * intree = (TTree*)infile->Get("tree");

    //intree->Print();
    //intree->GetBranch("slice")->Print();
    //intree->GetBranch("slice")->GetClassName();

    _slice *slc = nullptr;

    intree -> SetBranchAddress("slice",&slc);

    TH1D * protons_dedx = new TH1D("protons_dedx","",100,0,30);

    for(int i = 0; i < intree->GetEntries(); i++)
    {
        intree -> GetEntry(i);

        for(const auto &p : slc->_protons)
        {
            for(int hit = 0; hit < p._dedx.size(); hit++)
            {
                protons_dedx -> Fill(p._dedx[hit]);
            }
        }
    }

    protons_dedx -> Draw();
}