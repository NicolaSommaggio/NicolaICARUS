#include "slice_struct.h"

void CHECK_SELECTION()
{

    TFile * infile = TFile::Open("tree_outfile_OFFBEAM.root","READ");

    TTree * intree = (TTree*)infile->Get("tree");

    _slice *slc;

    intree -> SetBranchAddress("slice",&slc);

    TH1D *h_mu_pro_angle = new TH1D("h_mu_pro_angle","",100,-1,1);
    TH1D *h_mu_pro_angle_cut = new TH1D("h_mu_pro_angle_cut","",100,-1,1);

    for(int i = 0; i < intree -> GetEntries(); i ++)
    {
        intree -> GetEntry(i);

        h_mu_pro_angle -> Fill(slc->_mu_pro_angle);

        if(slc->_mu_pro_angle > -0.85){h_mu_pro_angle_cut -> Fill(slc->_mu_pro_angle);}
        
    }

    h_mu_pro_angle -> Draw();
    h_mu_pro_angle_cut -> Draw();

}