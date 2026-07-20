#include "slice_struct.h"

bool is1muNp_true(std::string true_class)
{
  if(
    true_class == "bad_slice" || 
    true_class == "not_in_FV" ||
    true_class == "not_in_Active" ||
    true_class == "not_contained" ||
    true_class == "cosmic" ||
    true_class == "reco_vertex_isNAN" ||
    true_class == "true_vertex_isNAN" ||
    true_class == "not_numu" || 
    true_class == "N/D"
  )return false;

  std::vector<int> true_c;
  for(char c : true_class){true_c.push_back(c - '0');}

  if(true_c[0]==1 && true_c[1]>0 && true_c[2]==0 and true_c[3]==0 and true_c[4]==0)return true;

  return false;
}


void CHECK_SELECTION()
{
    std::array<std::string,2> samples = {"OFFBEAM","MC"};

    TFile * outfile = new TFile("outfile.root","RECREATE");

    for(const auto &s : samples)
    {
        TFile * infile = TFile::Open(Form("tree_outfile_%s_withp.root",s.c_str()),"READ");
        infile -> cd();
        
        TTree * intree = (TTree*)infile->Get("tree");

        //intree->Print();
        //intree->GetBranch("slice")->Print();
        //intree->GetBranch("slice")->GetClassName();

        _slice *slc = nullptr;

        intree -> SetBranchAddress("slice",&slc);

        TH1D *h_mu_pro_angle = new TH1D(Form("h_mu_pro_angle_%s",s.c_str()),"",100,-1,1);
        TH1D *h_mu_pro_angle_cut = new TH1D(Form("h_mu_pro_angle_cut_%s",s.c_str()),"",100,-1,1);
        TH1D *h_p_as_pro = new TH1D(Form("h_p_as_pro_%s",s.c_str()),"",100,0,1);
        TH1D *h_p_as_pro_true1muNp = new TH1D(Form("h_p_as_pro_true1muNp_%s",s.c_str()),"",100,0,1);

        for(int i = 0; i < intree -> GetEntries(); i ++)
        {
            intree -> GetEntry(i);

            h_mu_pro_angle -> Fill(slc->_mu_pro_angle);

            if(slc->_mu_pro_angle > -0.85){h_mu_pro_angle_cut -> Fill(slc->_mu_pro_angle);}

            for (const auto &pro : slc->_protons)
            {
                //for(const auto &p : pro._proba){cout << p << " ";} cout << endl;
                
                if(pro._proba[2] > pro._proba[3]){h_p_as_pro->Fill(pro._proba[2]);}
                else{h_p_as_pro->Fill(pro._proba[3]);}
                
            }

            if(is1muNp_true(slc -> _true_class))
            {
                for (const auto &pro : slc->_protons)
                {
                    if(pro._proba[2] > pro._proba[3]){h_p_as_pro_true1muNp->Fill(pro._proba[2]);}
                    else{h_p_as_pro_true1muNp->Fill(pro._proba[3]);}
                }
            }
        }

        outfile -> cd();

        h_mu_pro_angle->Scale(1./h_mu_pro_angle->Integral());
        h_mu_pro_angle_cut->Scale(1./h_mu_pro_angle_cut->Integral());
        h_p_as_pro->Scale(1./h_p_as_pro->Integral());
        h_p_as_pro_true1muNp->Scale(1./h_p_as_pro_true1muNp->Integral());

        h_mu_pro_angle -> Write();
        h_mu_pro_angle_cut -> Write();
        h_p_as_pro -> Write();
        h_p_as_pro_true1muNp -> Write();


    }

}