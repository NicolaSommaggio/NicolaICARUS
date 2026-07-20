#include "slice_struct.h"

void CHECK_SELECTION()
{

    TFile * infile = TFile::Open("tree_outfile_OFFBEAM.root","READ");

    TTree * intree = (TTree*)infile->Get("tree");

    _slice *slc = nullptr;

    intree -> SetBranchAddress("slice",&slc);

    TH1D *h_mu_pro_angle = new TH1D("h_mu_pro_angle","",100,-1,1);
    TH1D *h_mu_pro_angle_cut = new TH1D("h_mu_pro_angle_cut","",100,-1,1);
    TH1D *h_p_as_pro_ris = new TH1D("h_p_as_pro_ris","",100,0,1);

    for(int i = 0; i < intree -> GetEntries(); i ++)
    {
        intree -> GetEntry(i);

        h_mu_pro_angle -> Fill(slc->_mu_pro_angle);

        if(slc->_mu_pro_angle > -0.85){h_mu_pro_angle_cut -> Fill(slc->_mu_pro_angle);}

        for(const auto &pro : slc -> _protons)
        {
            h_p_as_pro_ris -> Fill();
        }
        
    }

    TCanvas * c1 = new TCanvas("c1");
    h_mu_pro_angle -> Draw();
    TCanvas * c2 = new TCanvas("c2");
    h_mu_pro_angle_cut -> Draw();

}