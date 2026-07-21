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

    //double DEP_E_CUT = 50.;
    //double PROBA_CUT = 0.5;
    //double ANGLE_CUT = -0.85;

    double POT_MC = 2.43007e+20;
    double POT_OFFBEAM = 2.11e+20;

    std::array<std::string,2> samples = {"OFFBEAM","MC"};

    TFile * outfile = new TFile("outfile.root","RECREATE");

    TGraph *gr_eff = new TGraph();
    gr_eff->SetName("gr_eff");
    TGraph *gr_pur = new TGraph();
    gr_pur->SetName("gr_pur");
    TGraph *gr_eff_times_pur = new TGraph();
    gr_eff_times_pur->SetName("gr_eff_times_pur");


    //std::array<double,9> VEC_DEP_E_CUT = {0., 25, 30, 35, 40, 45, 50, 55, 60};
    //std::array<double,13> VEC_PROBA_CUT = {0., 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.56, 0.57, 0.58, 0.59, 0.6, 0.65};
    //std::array<double,10> VEC_ANGLE_CUT = {-1., -0.95, -0.93, -0.9, -0.91, -0.89, -0.87, -0.85, -0.83, -0.81};
    //std::array<double,1> VEC_ANGLE_CUT = {-0.85};

    std::array<double,1> VEC_DEP_E_CUT = {0.};
    std::array<double,1> VEC_PROBA_CUT = {0.};
    std::array<double,1> VEC_ANGLE_CUT = {-1.};

    int iteration = 0;

    for(const auto &PROBA_CUT : VEC_PROBA_CUT)
    {
        for(const auto &DEP_E_CUT : VEC_DEP_E_CUT)
        {
            for(const auto &ANGLE_CUT : VEC_ANGLE_CUT)
            {
                iteration++;

                double tot_selected_MC = 0;
                double tot_selected_true_MC = 0;
                double tot_selected_OFFBEAM = 0;
                
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
                    TH1D *h_depE = new TH1D(Form("h_depE_%s",s.c_str()),"",100,0,200);
                    TH1D *h_depE_true1muNp = new TH1D(Form("h_depE_true1muNp_%s",s.c_str()),"",100,0,200);
                    TH2D *h_depE_vs_p = new TH2D(Form("h_depE_vs_p_%s",s.c_str()),"",100,0,1,100,0,200);
                    TH1D *h_crlongtrkdiry = new TH1D(Form("h_crlongtrkdiry_%s",s.c_str()),"",100,-5.,5.);
                    TH1D *h_crlongtrkdiry_true1muNp = new TH1D(Form("h_crlongtrkdiry_true1muNp_%s",s.c_str()),"",100,-5.,5.);
                    TH1D *h_theta_xw = new TH1D(Form("h_theta_xw_%s",s.c_str()),"",100,0,90);
                    TH1D *h_theta_xw_true1muNp = new TH1D(Form("h_theta_xw_true1muNp_%s",s.c_str()),"",100,0,90);
                    TH2D *h_theta_xw_vs_length = new TH2D(Form("h_theta_xw_vs_length_%s",s.c_str()),"",100,0,90,100,0,100);
                    TH2D *h_theta_xw_vs_length_true1muNp = new TH2D(Form("h_theta_xw_vs_length_true1muNp_%s",s.c_str()),"",100,0,90,100,0,100);

                    int tot_selected = 0;
                    int tot_selected_true = 0;
                    int tot_selected_true_precut = 0;

                    for(int i = 0; i < intree -> GetEntries(); i ++)
                    {
                        intree -> GetEntry(i);

                        h_mu_pro_angle -> Fill(slc->_mu_pro_angle);

                        h_crlongtrkdiry -> Fill(slc->_crlongtrkdiry);

                        if(slc->_mu_pro_angle > -0.85){h_mu_pro_angle_cut -> Fill(slc->_mu_pro_angle);}

                        int num_protons = 0;
                        for (const auto &pro : slc->_protons)
                        {
                            //for(const auto &p : pro._proba){cout << p << " ";} cout << endl;
                
                            if(pro._proba[2] > pro._proba[3]){h_p_as_pro->Fill(pro._proba[2]);}
                            else{h_p_as_pro->Fill(pro._proba[3]);}

                            h_depE -> Fill(pro._depE);

                            h_theta_xw -> Fill(std::abs(pro._theta_xw)*180/M_PI);

                            h_theta_xw_vs_length -> Fill(std::abs(pro._theta_xw)*180/M_PI,pro._length);

                            double best_p_proba = (pro._proba[2] > pro._proba[3]) ? pro._proba[2] : pro._proba[3];

                            h_depE_vs_p->Fill(best_p_proba,pro._depE);

                            //if(best_p_proba > 0.5)num_protons++;
                            if(best_p_proba > PROBA_CUT)
                            {
                                if(pro._proba[2] > pro._proba[3])
                                {
                                    if(pro._depE > DEP_E_CUT)num_protons++;
                                }
                                else num_protons++;
                            }
                        }
                        if(num_protons > 0 && slc->_mu_pro_angle > ANGLE_CUT)tot_selected++;

                        int num_protons_true_slice = 0;
                        if(is1muNp_true(slc -> _true_class))
                        {

                            h_crlongtrkdiry_true1muNp -> Fill(slc->_crlongtrkdiry);

                            for (const auto &pro : slc->_protons)
                            {
                                if(pro._proba[2] > pro._proba[3]){h_p_as_pro_true1muNp->Fill(pro._proba[2]);}
                                else{h_p_as_pro_true1muNp->Fill(pro._proba[3]);}

                                h_depE_true1muNp -> Fill(pro._depE);

                                h_theta_xw_true1muNp -> Fill(std::abs(pro._theta_xw)*180/M_PI);

                                h_theta_xw_vs_length_true1muNp -> Fill(std::abs(pro._theta_xw)*180/M_PI,pro._length);

                                double best_p_proba = (pro._proba[2] > pro._proba[3]) ? pro._proba[2] : pro._proba[3];
                    
                                //if(best_p_proba > 0.5)num_protons_true_slice++;
                                if(best_p_proba > PROBA_CUT)
                                {
                                    if(pro._proba[2] > pro._proba[3])
                                    {
                                        if(pro._depE > DEP_E_CUT)num_protons_true_slice++;
                                    }
                                    else num_protons_true_slice++;
                                }
                            }

                            tot_selected_true_precut++;
                        }
                        if(num_protons_true_slice > 0 && slc->_mu_pro_angle > ANGLE_CUT)tot_selected_true++;
                    }

                    //cout << "*** " << s << " ***" << endl;
                    //cout << "NO CUT : " << intree -> GetEntries() << " " << tot_selected_true_precut << " | EFF : " << tot_selected_true_precut / 29503. << endl;
                    //cout << "WITH CUTS : " << tot_selected << " " << tot_selected_true << " | EFF : " << tot_selected_true / 29503. << endl;
                    if(s == "MC")
                    {
                        tot_selected_MC = tot_selected;
                        tot_selected_true_MC = tot_selected_true;
                    }
                    else if(s == "OFFBEAM")
                    {
                        tot_selected_OFFBEAM = tot_selected;
                    }


                    outfile -> cd();

                    h_mu_pro_angle->Scale(1./h_mu_pro_angle->Integral());
                    h_mu_pro_angle_cut->Scale(1./h_mu_pro_angle_cut->Integral());
                    h_p_as_pro->Scale(1./h_p_as_pro->Integral());
                    h_p_as_pro_true1muNp->Scale(1./h_p_as_pro_true1muNp->Integral());
                    h_depE -> Scale(1./h_depE->Integral());
                    h_depE_true1muNp -> Scale(1./h_depE_true1muNp->Integral());
                    h_crlongtrkdiry_true1muNp -> Scale(1./h_crlongtrkdiry_true1muNp->Integral());
                    h_crlongtrkdiry -> Scale(1./h_crlongtrkdiry->Integral());
                    h_theta_xw -> Scale(1./h_theta_xw->Integral());
                    h_theta_xw_true1muNp -> Scale(1./h_theta_xw_true1muNp->Integral());

                    h_mu_pro_angle -> Write();
                    h_mu_pro_angle_cut -> Write();
                    h_p_as_pro -> Write();
                    h_p_as_pro_true1muNp -> Write();
                    h_depE -> Write();
                    h_depE_true1muNp -> Write();
                    h_depE_vs_p -> Write();
                    h_crlongtrkdiry_true1muNp -> Write();
                    h_crlongtrkdiry -> Write();
                    h_theta_xw -> Write();
                    h_theta_xw_true1muNp -> Write();
                    h_theta_xw_vs_length_true1muNp -> Write();
                    h_theta_xw_vs_length -> Write();

                }

                double eff = tot_selected_true_MC / 29503.;
                double offbeam_rescaled = POT_MC / POT_OFFBEAM * tot_selected_OFFBEAM;
                double pur = tot_selected_true_MC / (offbeam_rescaled + tot_selected_MC);
                //cout << "DEP_E CUT : " << DEP_E_CUT << " PROBA_CUT : " << PROBA_CUT << " ANGLE_CUT : " << ANGLE_CUT << " | eff : " << eff << " pur : " <<  pur << endl;
                //cout << "PROBA_CUT : " << PROBA_CUT << " DEP_E CUT : " << DEP_E_CUT << " | eff : " << eff << " pur : " <<  pur << endl;
                //cout << "PROBA_CUT : " << PROBA_CUT << " | eff : " << eff << " pur : " <<  pur << endl;
                //cout << "PROBA_CUT : " << PROBA_CUT << " ANGLE_CUT : " << ANGLE_CUT  << " | eff : " << eff << " pur : " <<  pur << endl;

                gr_eff -> SetPoint(iteration, PROBA_CUT, eff);
                gr_pur -> SetPoint(iteration, PROBA_CUT, pur);
                gr_eff_times_pur -> SetPoint(iteration, PROBA_CUT, eff*pur);
            }
        }
    }

    outfile->cd();
    gr_eff -> Write();
    gr_pur -> Write();
    gr_eff_times_pur -> Write();
    

}