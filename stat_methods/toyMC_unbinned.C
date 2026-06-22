std::vector<double> x_var;

Double_t MU1_FIT_D_PLUS_PHI_PI = 1.86969e+00;
Double_t MU2_FIT_D_PLUS_PHI_PI = 1.87035e+00;
Double_t SIGMA1_FIT_D_PLUS_PHI_PI = 5.79484e-03;
Double_t SIGMA2_FIT_D_PLUS_PHI_PI = 8.15639e-03;
Double_t FRACTION_FIT_D_PLUS_PHI_PI = 8.83747e-01;

Double_t M0_FIT_DS_PLUS_PHI_MU_NU = 1.97801e+00;
Double_t C_FIT_DS_PLUS_PHI_MU_NU = -3.82391e+00;
Double_t P_FIT_DS_PLUS_PHI_MU_NU = 1.41325e+00;

Double_t MU1_FIT_DS_PLUS_PHI_PI = 1.96843e+00;
Double_t MU2_FIT_DS_PLUS_PHI_PI = 1.96807e+00;
Double_t SIGMA1_FIT_DS_PLUS_PHI_PI = 6.08795e-03;
Double_t SIGMA2_FIT_DS_PLUS_PHI_PI = 7.87780e-03;
Double_t FRACTION_FIT_DS_PLUS_PHI_PI = 8.65940e-01;

Double_t MU1_FIT_MASS_DS_PLUS_TAU_NU = 1.77687e+00;
Double_t MU2_FIT_MASS_DS_PLUS_TAU_NU = 1.77748e+00;
Double_t SIGMA1_FIT_MASS_DS_PLUS_TAU_NU = 5.42707e-03;
Double_t SIGMA2_FIT_MASS_DS_PLUS_TAU_NU = 7.16697e-03;
Double_t FRACTION_FIT_MASS_DS_PLUS_TAU_NU = 7.88093e-01;

Double_t TAU_FIT_COMBINATORIAL = 1.03105;

Double_t TOT_MASS_LOW = 1.6;
Double_t TOT_MASS_HIGH = 2.1;

//Double_t f_D_PhiPi = 0.0214727;
//Double_t f_Ds_PhiNuMu = 0.652866;
//Double_t f_Ds_PhiPi = 0.04305;

Double_t f_D_PhiPi = 0.0215021;
Double_t f_Ds_PhiMuNu = 0.653769;
Double_t f_Ds_PhiPi = 0.0431105;


Double_t EFF_D_PhiPi = 53699./50e6;
Double_t EFF_Ds_PhiMuNu = 95217./5e6;
Double_t EFF_Ds_PhiPi = 32106./50e6;
Double_t EFF_Sig = 86069./1e6;

Double_t SIGMA_EFF_D_PhiPi = std::sqrt((EFF_D_PhiPi*(1-EFF_D_PhiPi))/50e6);
Double_t SIGMA_EFF_Ds_PhiMuNu = std::sqrt((EFF_Ds_PhiMuNu*(1-EFF_Ds_PhiMuNu))/5e6);
Double_t SIGMA_EFF_Ds_PhiPi = std::sqrt((EFF_Ds_PhiPi*(1-EFF_Ds_PhiPi))/50e6);
Double_t SIGMA_EFF_Sig = std::sqrt((EFF_Sig*(1-EFF_Sig))/1e6);

Double_t BR_D_PhiPi = 2.69e-3;
Double_t BR_Ds_PhiMuNu = 2.24e-2;
Double_t BR_Ds_PhiPi = 2.25e-2;
Double_t BR_Ds_TauNu = 5.39e-2;
Double_t BR_Phi_KK = 49.9e-2;

Double_t SIGMA_BR_D_PhiPi = 0.08e-3;
Double_t SIGMA_BR_Ds_PhiMuNu = 0.11e-2;
Double_t SIGMA_BR_Ds_PhiPi = 0.05e-2;
Double_t SIGMA_BR_Ds_TauNu = 0.09e-2;
Double_t SIGMA_BR_Phi_KK = 0.5e-2;

Double_t CS_pp_D = 834;
Double_t CS_pp_Ds = 353;

Double_t SIGMA_CS_pp_D = 2;
Double_t SIGMA_CS_pp_Ds = 9;

Double_t R_STEP = 0.01;
//Double_t INTEGRATION_STEP = 0.01;
Double_t f_Sig_STEP = 0.00005;

Double_t Gauss_pdf(Double_t *x, Double_t *par)
{
  return TMath::Gaus(x[0],par[0],par[1],1);
}

Double_t Gauss_pdf_bin(Double_t *x, Double_t *par)
{
  return par[2] * par[3] * TMath::Gaus(x[0],par[0],par[1],1);
}

Double_t sample_val(Double_t mu, Double_t sigma)
{
  Double_t par[2] = {mu,sigma};

  TF1 *f = new TF1("f",Gauss_pdf,mu - 5 * sigma, mu + 5 * sigma, 2);
  f->SetParameters(mu,sigma);
  Double_t extracted_val = f->GetRandom(mu - 5 * sigma,mu + 5 * sigma);

  delete f;

  return extracted_val;
}

Double_t Likelihood_ratio(Double_t *x, Double_t *par)
{
  Double_t val;
  if(x[0] >= 0) val = TMath::Exp(-1./2 * std::pow((x[0] - par[0]),2) / std::pow(par[1],2) );
  else val = TMath::Exp((x[0] * par[0] - par[0] * par[0] / 2) / std::pow(par[1],2));

  return val;
}

std::pair<Double_t,Double_t> get_x1_x2_at_R(Double_t R, Double_t mu, Double_t sigma)
{
  Double_t x_right = mu + std::sqrt(-2 * sigma * sigma * std::log(R));
  
  Double_t x_left_1 = mu - std::sqrt(-2 * sigma * sigma * std::log(R));

  Double_t x_left_2 = sigma * sigma * std::log(R)/mu + mu/2.;

  Double_t x_left;
  if(x_left_1 > 0){x_left = x_left_1;}
  else{x_left = x_left_2;}

  return {x_left,x_right};
}

/*
Double_t numeric_gaus_integral(Double_t x1, Double_t x2, Double_t mu, Double_t sigma)
{
  Double_t sum = 0;
  Double_t pars[2] = {mu,sigma};

  for(Double_t x = x1; x <= x2; x += INTEGRATION_STEP)
  {
    sum = sum + Gauss_pdf(&x,pars);
  }

  return sum * INTEGRATION_STEP;
}
*/

std::pair<Double_t, Double_t> get_acceptance_interval(Double_t mu, Double_t sigma)
{
  for(Double_t R = 0.99; R > 0; R -= R_STEP )
  {
    std::pair<Double_t,Double_t> x_interval = get_x1_x2_at_R(R,mu,sigma);

    //Double_t integral = numeric_gaus_integral(x_interval.first, x_interval.second, mu, sigma);

    Double_t integral = ROOT::Math::normal_cdf(x_interval.second,sigma,mu) - ROOT::Math::normal_cdf(x_interval.first,sigma,mu);

    if(integral >= 0.9){return x_interval;} 
  }

  return {-1,-1};
}

#define analysis_cxx
#include "fit_unbinned_TotalSpectrum_pseudo.h"
#include "fit_unbinned_Gauss.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <stdio.h>
#include <stdlib.h> 
#include <iostream>

using namespace std;


void toy() 
{

  TFile * outfile = TFile::Open("outfile.root","RECREATE");
  TFile * dummy = TFile::Open("dummy.root","RECREATE");
  // TFile * data = TFile::Open("data.root","RECREATE");

  // TTree *tree = new TTree("sim_tree", "Simulation Tree");


  gRandom ->SetSeed (121001);
  const double mMin = 1.60;
  const double mMax = 2.10;
  const int nBins = 100;
  const double binWidth = (mMax - mMin) / nBins;
  const int nToys = 500;

  const double nTotTrue = 74439; // 1000.0;

  std::vector<double> lower_belt;
  std::vector<double> upper_belt;
  std::vector<double> f_sig_belt;
  std::vector<double> lower_belt_BR;
  std::vector<double> upper_belt_BR;
  std::vector<double> f_sig_belt_BR;

  std::vector<double> lower_belt_RTau;
  std::vector<double> upper_belt_RTau;
  std::vector<double> f_sig_belt_RTau;

  double m;
  // tree->Branch("m", &m);

int point = -1;
// Double_t fSigTrue = 0.01;
// for(Double_t fSigTrue = 0.001; fSigTrue <= 0.1; fSigTrue += f_Sig_STEP)
for(Double_t fSigTrue = 0.0001; fSigTrue <= 0.005; fSigTrue += f_Sig_STEP)
{

  std::vector<double> fsig;
  std::vector<double> fdiffsig;
  std::vector<double> fpullsig;
  std::vector<double> fsig_over_fphiMuNu;

  std::vector<double> BR_Sig;
  std::vector<double> R_tau;
  std::vector<double> sigma_f_sig;

  f_sig_belt.push_back(fSigTrue);

  cout << "PROCESSING f_Sig = " << fSigTrue << endl;
  point ++;

  // Known model used to generate the pseudo - experiments .
  // TF1 :: GetRandom uses only the shape of the function.
  TF1 generatorModel(Form("generatorModel_%d",point), total_mass_spectrum_pdf , mMin , mMax, 4);
  generatorModel.SetParameters(f_D_PhiPi, f_Ds_PhiMuNu, f_Ds_PhiPi, fSigTrue);
  TH1D hfSigFit(Form("hfSigFit_%d",point), "Fitted signal yield;#hat{f}_{sig};Pseudo experiments", 50, 0., 0.);
  TH1D hdifffSigFit(Form("hdifffSigFit_%d",point), "Fitted signal yield - Simulated signal yield ; #hat{f_{sig}}-f_{sig} ; Pseudo experiments", 50, 0., 0.);
  TH1D hpullSigFit(Form("hpullSigFit_%d",point), "(Fitted signal yield - Simulated signal yield)/Sigma fitted signal yield ; (#hat{f_{sig}}-f_{sig})/#sigma_{#hat{f_{sig}}} ; Pseudo experiments", 200, 0., 0.);
  TH1D hfSig_over_fPhiMuNu(Form("hfSig_over_fPhiMuNu_%d",point),"Fitted signal yield / Fitted #var_phi#rightarrow#mu#nu yield ; hat{f_{sig}}/#hat{f_{#var_phi#rightarrow#mu#nu}} ; Pseudo experiments", 100, 0., 0.);
  TH1D hfBRSig(Form("hfBRSig_%d",point),"BR(#tau^+#rightarrow#var_phi#mu^+) ; BR(#tau^+#rightarrow#var_phi#mu^+) ; Pseudo experiments", 100, 0., 0.);
  TH1D hfR_tau(Form("hfR_tau_%d",point),"R_{#tau} ; R_{#tau} ; Pseudo experiments", 100, 0., 0.);
  TH1D *h_sigma_f_Sig = new TH1D(Form("h_sigma_f_Sig_%d",point),"",100,0.,0.); 

  for (int iToy = 0; iToy < nToys; ++iToy) 
  {
    if (iToy % 100 == 0){
      cout << "PROCESSING " << iToy << " TOY" << endl;

    }

    TH1D *hToy = new TH1D(Form("hToy_%d_%d",point,iToy), "", nBins, mMin, mMax);
    // const int nObs = gRandom ->Poisson(nTotTrue); // Extended toy generation .
    const int nObs = nTotTrue;

    std::vector<double> vec_generated;
    for (int i = 0; i < nObs; ++i) 
    {
      m = generatorModel.GetRandom(mMin , mMax);
      hToy -> Fill(m);
      vec_generated.push_back(m);
      // tree -> Fill();
    }


    std::tuple<Double_t,Double_t,Double_t,Double_t> tuple_fSigfromFIT = fit_unbinned_TotalSpectrum(
                            vec_generated, 
                            {0.021, 0.654, 0.043, 0.01}, 
                            {0.0001, 0.0002, 0.0001, 0.0001}, 
                            //{0., 0., 0., 0.0001},
			                      {0., 0., 0., 0.}, 
                            {0., 0., 0., 0.}, 
                            {"f_D_PhiPi","f_DS_PhiMuNu","f_DS_PhiPi","f_DS_TauNu"}, 
                            hToy, 
                            dummy, 
                            Form("Generated_Invariant_Mass_Spectrum_%d_%d",point,iToy), 
                            mMin, 
                            mMax
                          );


    Double_t fSigfromFIT     = std::get<0>(tuple_fSigfromFIT);
    Double_t err_fSigfromFIT = std::get<1>(tuple_fSigfromFIT);
    Double_t fPhiMuNufromFIT = std::get<2>(tuple_fSigfromFIT);
    Double_t err_fPhiMuNufromFIT = std::get<3>(tuple_fSigfromFIT);

    hfSigFit.Fill(fSigfromFIT);
    fsig.push_back(fSigfromFIT);
    hdifffSigFit.Fill(fSigfromFIT-fSigTrue);
    fdiffsig.push_back(fSigfromFIT-fSigTrue);
    hpullSigFit.Fill((fSigfromFIT-fSigTrue)/err_fSigfromFIT);
    fpullsig.push_back((fSigfromFIT-fSigTrue)/err_fSigfromFIT);
    hfSig_over_fPhiMuNu.Fill(fSigfromFIT/fPhiMuNufromFIT);
    fsig_over_fphiMuNu.push_back(fSigfromFIT/fPhiMuNufromFIT);

    double BR = ( fSigfromFIT / ((1 - fSigfromFIT) * fPhiMuNufromFIT) ) * ( ( sample_val(BR_Ds_PhiMuNu,SIGMA_BR_Ds_PhiMuNu) * sample_val(EFF_Ds_PhiMuNu,SIGMA_EFF_Ds_PhiMuNu) ) / ( sample_val(BR_Ds_TauNu,SIGMA_BR_Ds_TauNu) * sample_val(EFF_Sig,SIGMA_EFF_Sig) ) );

    hfBRSig.Fill(BR);
    BR_Sig.push_back(BR);
    hfR_tau.Fill( (fSigfromFIT * sample_val(EFF_Ds_PhiMuNu,SIGMA_EFF_Ds_PhiMuNu)) / (fPhiMuNufromFIT * sample_val(EFF_Sig,SIGMA_EFF_Sig)) );
    R_tau.push_back( (fSigfromFIT * sample_val(EFF_Ds_PhiMuNu,SIGMA_EFF_Ds_PhiMuNu)) / (fPhiMuNufromFIT * sample_val(EFF_Sig,SIGMA_EFF_Sig)) );

    sigma_f_sig.push_back(err_fSigfromFIT);
    h_sigma_f_Sig -> Fill(err_fSigfromFIT);

    delete hToy;

  }

    outfile -> cd();


  std::tuple <Double_t,Double_t> results_fit_fsig = fit_unbinned_Gauss( fsig, 
                                                                  {hfSigFit.GetMean() , hfSigFit.GetStdDev()}, 
                                                                  {0.00001, 0.00001}, 
                                                                  {0., 0.}, 
                                                                  {0., 0.}, 
                                                                  {"mu", "sigma"}, 
                                                                  &hfSigFit, 
                                                                  outfile, 
                                                                  Form("fsig_%d",point), 
                                                                  hfSigFit.GetXaxis()->GetXmin(), 
                                                                  hfSigFit.GetXaxis()->GetXmax()
                                                                  );

  cout << "\n\n **** FIT f_Sig : MU = " << std::get<0>(results_fit_fsig) << " SIGMA : " <<  std::get<1>(results_fit_fsig) << " ****\n\n" << endl;
  
  Double_t mu_f_sig = std::get<0>(results_fit_fsig);
  Double_t sigma_mu_f_sig = std::get<1>(results_fit_fsig);
  
  std::pair<Double_t,Double_t> acceptance_interval = get_acceptance_interval(mu_f_sig,sigma_mu_f_sig);

  lower_belt.push_back(acceptance_interval.first);
  upper_belt.push_back(acceptance_interval.second);

  std::tuple <Double_t,Double_t> results_fit_residuals = fit_unbinned_Gauss(  fdiffsig, 
                                                                        {hdifffSigFit.GetMean(), hdifffSigFit.GetStdDev()}, 
                                                                        {0.00001, 0.00001}, 
                                                                        {0., 0.}, 
                                                                        {0., 0.}, 
                                                                        {"mu", "sigma"}, 
                                                                        &hdifffSigFit, 
                                                                        outfile, 
                                                                        Form("residuals_%d",point), 
                                                                        hdifffSigFit.GetXaxis()->GetXmin(), 
                                                                        hdifffSigFit.GetXaxis()->GetXmax()
                                                                      );

  cout << "\n\n **** FIT RESIDUALS : MU = " << std::get<0>(results_fit_residuals) << " SIGMA : " <<  std::get<1>(results_fit_residuals) << " ****\n\n" << endl;

  std::tuple <Double_t,Double_t> results_fit_pulls = fit_unbinned_Gauss(  fpullsig, 
                                                                    {hpullSigFit.GetMean(), hpullSigFit.GetStdDev()}, 
                                                                    {0.00001, 0.00001}, 
                                                                    {0., 0.}, 
                                                                    {0., 0.}, 
                                                                    {"mu", "sigma"},
                                                                    &hpullSigFit, 
                                                                    outfile, 
                                                                    Form("pulls_%d",point), 
                                                                    hpullSigFit.GetXaxis()->GetXmin(), 
                                                                    hpullSigFit.GetXaxis()->GetXmax()
                                                                  );  

  cout << "\n\n **** FIT PULLS : MU = " << std::get<0>(results_fit_pulls) << " SIGMA : " <<  std::get<1>(results_fit_pulls) << " ****\n\n" << endl;                                                             

  
  std::tuple <Double_t,Double_t> results_fit_fSig_over_fPhiMuNu = fit_unbinned_Gauss(  fsig_over_fphiMuNu, 
                                                                    {hfSig_over_fPhiMuNu.GetMean(), hfSig_over_fPhiMuNu.GetStdDev()}, 
                                                                    {0.00001, 0.00001}, 
                                                                    {0., 0.}, 
                                                                    {0., 0.}, 
                                                                    {"mu", "sigma"},
                                                                    &hfSig_over_fPhiMuNu, 
                                                                    outfile, 
                                                                    Form("fSig_over_fPhiMuNu_%d",point), 
                                                                    hfSig_over_fPhiMuNu.GetXaxis()->GetXmin(), 
                                                                    hfSig_over_fPhiMuNu.GetXaxis()->GetXmax()
                                                                  ); 
  cout << "\n\n **** FIT fSig over fPhiMuNu : MU = " << std::get<0>(results_fit_fSig_over_fPhiMuNu) << " SIGMA : " <<  std::get<1>(results_fit_fSig_over_fPhiMuNu) << " ****\n\n" << endl;                                                                
  

  
  std::tuple <Double_t,Double_t> results_fit_BRSig = fit_unbinned_Gauss(  BR_Sig, 
                                                                    {hfBRSig.GetMean(), hfBRSig.GetStdDev()}, 
                                                                    {0.00001, 0.00001}, 
                                                                    {0., 0.}, 
                                                                    {0., 0.}, 
                                                                    {"mu", "sigma"},
                                                                    &hfBRSig, 
                                                                    outfile, 
                                                                    Form("BRSig_%d",point), 
                                                                    hfBRSig.GetXaxis()->GetXmin(), 
                                                                    hfBRSig.GetXaxis()->GetXmax()
                                                                  );  

  cout << "\n\n **** FIT BR : MU = " << std::get<0>(results_fit_BRSig) << " SIGMA : " <<  std::get<1>(results_fit_BRSig) << " ****\n\n" << endl;                                                                  
                  
  Double_t mu_BRSig = std::get<0>(results_fit_BRSig);
  Double_t sigma_mu_BRSig = std::get<1>(results_fit_BRSig);
  
  std::pair<Double_t,Double_t> acceptance_interval_BR = get_acceptance_interval(mu_BRSig,sigma_mu_BRSig);

  lower_belt_BR.push_back(acceptance_interval_BR.first);
  upper_belt_BR.push_back(acceptance_interval_BR.second);
  f_sig_belt_BR.push_back(mu_BRSig);

  
  
  std::tuple <Double_t,Double_t> results_fit_R_tau = fit_unbinned_Gauss(  R_tau, 
                                                                    {hfR_tau.GetMean(), hfR_tau.GetStdDev()}, 
                                                                    {0.00001, 0.00001}, 
                                                                    {0., 0.}, 
                                                                    {0., 0.}, 
                                                                    {"mu", "sigma"},
                                                                    &hfR_tau, 
                                                                    outfile, 
                                                                    Form("R_tau_%d",point), 
                                                                    hfR_tau.GetXaxis()->GetXmin(), 
                                                                    hfR_tau.GetXaxis()->GetXmax()
                                                                  );  
                                                                  
  cout << "\n\n **** FIT R_TAU : MU = " << std::get<0>(results_fit_R_tau) << " SIGMA : " <<  std::get<1>(results_fit_R_tau) << " ****\n\n" << endl; 
  
  Double_t mu_RTauSig = std::get<0>(results_fit_R_tau);
  Double_t sigma_mu_RTauSig = std::get<1>(results_fit_R_tau);
  
  std::pair<Double_t,Double_t> acceptance_interval_RTau = get_acceptance_interval(mu_RTauSig,sigma_mu_RTauSig);

  lower_belt_RTau.push_back(acceptance_interval_RTau.first);
  upper_belt_RTau.push_back(acceptance_interval_RTau.second);
  f_sig_belt_RTau.push_back(mu_RTauSig);

  std::tuple <Double_t,Double_t> results_fit_sigma_f_sigma = fit_unbinned_Gauss(  sigma_f_sig, 
                                                                    {h_sigma_f_Sig->GetMean(), h_sigma_f_Sig->GetStdDev()}, 
                                                                    {0.00001, 0.00001}, 
                                                                    {0., 0.}, 
                                                                    {0., 0.}, 
                                                                    {"mu", "sigma"},
                                                                    h_sigma_f_Sig, 
                                                                    outfile, 
                                                                    Form("sigma_fsig_%d",point), 
                                                                    h_sigma_f_Sig->GetXaxis()->GetXmin(), 
                                                                    h_sigma_f_Sig->GetXaxis()->GetXmax()
                                                                  );  
                                                                  
  cout << "\n\n **** FIT SIGMA f_Sig : MU = " << std::get<0>(results_fit_sigma_f_sigma) << " SIGMA : " <<  std::get<1>(results_fit_sigma_f_sigma) << " ****\n\n" << endl; 

  /*
  TF1 *fit_gaus_binned = new TF1(Form("fit_gaus_binned_%d",point),Gauss_pdf_bin,-10,10,4);
  fit_gaus_binned -> SetParameters(0.05, 1., hpullSigFit.Integral(), hpullSigFit.GetBinWidth(1));
  fit_gaus_binned -> FixParameter(3,hpullSigFit.GetBinWidth(1));
  hpullSigFit.Fit(fit_gaus_binned,"L");

  hpullSigFit.Write();
  */

}
// data->cd();
// tree->Write();
// data->Close();
TGraph * g_lower_belt = new TGraph(lower_belt.size(),f_sig_belt.data(),lower_belt.data());
TGraph * g_upper_belt = new TGraph(upper_belt.size(),f_sig_belt.data(),upper_belt.data());

TSpline3 *s_lower_belt = new TSpline3("spline_lower_belt",g_lower_belt);
TSpline3 *s_upper_belt = new TSpline3("spline_upper_belt",g_upper_belt);

TCanvas *c_belt = new TCanvas("c_belt","");
c_belt->cd();
g_lower_belt->SetMarkerStyle(7);
g_lower_belt->SetMarkerSize(2);
g_upper_belt->SetMarkerStyle(7);
g_upper_belt->SetMarkerSize(2);
s_lower_belt->SetLineColor(kBlue);
s_upper_belt->SetLineColor(kBlue);
s_lower_belt->SetLineWidth(2);
s_upper_belt->SetLineWidth(2);
g_upper_belt->Draw("AP");
g_lower_belt->Draw("SAME");
s_lower_belt->Draw("SAME");
s_upper_belt->Draw("SAME");
g_upper_belt->SetMinimum(lower_belt.data()[0]-0.1*lower_belt.data()[0]);

g_lower_belt->Write("f_sig_g_lower_belt");
g_upper_belt->Write("f_sig_g_upper_belt");
s_lower_belt->Write("f_sig_s_lower_belt");
s_upper_belt->Write("f_sig_s_upper_belt");
c_belt->Write();

TGraph * g_lower_belt_BR = new TGraph(lower_belt_BR.size(),f_sig_belt_BR.data(),lower_belt_BR.data());
TGraph * g_upper_belt_BR = new TGraph(upper_belt_BR.size(),f_sig_belt_BR.data(),upper_belt_BR.data());

TSpline3 *s_lower_belt_BR = new TSpline3("spline_lower_belt_BR",g_lower_belt_BR);
TSpline3 *s_upper_belt_BR = new TSpline3("spline_upper_belt_BR",g_upper_belt_BR);

TCanvas *c_belt_BR = new TCanvas("c_belt_BR","");
c_belt_BR->cd();
g_lower_belt_BR->SetMarkerStyle(7);
g_lower_belt_BR->SetMarkerSize(2);
g_upper_belt_BR->SetMarkerStyle(7);
g_upper_belt_BR->SetMarkerSize(2);
s_lower_belt_BR->SetLineColor(kBlue);
s_upper_belt_BR->SetLineColor(kBlue);
s_lower_belt_BR->SetLineWidth(2);
s_upper_belt_BR->SetLineWidth(2);
g_upper_belt_BR->Draw("AP");
g_lower_belt_BR->Draw("SAME");
s_lower_belt_BR->Draw("SAME");
s_upper_belt_BR->Draw("SAME");
g_upper_belt_BR->SetMinimum(lower_belt_BR.data()[0]*1.1);

g_lower_belt_BR->Write("g_lower_belt_BR");
g_upper_belt_BR->Write("g_upper_belt_BR");
s_lower_belt_BR->Write("s_lower_belt_BR");
s_upper_belt_BR->Write("s_upper_belt_BR");
c_belt_BR->Write();

TGraph * g_lower_belt_RTau = new TGraph(lower_belt_RTau.size(),f_sig_belt_RTau.data(),lower_belt_RTau.data());
TGraph * g_upper_belt_RTau = new TGraph(upper_belt_RTau.size(),f_sig_belt_RTau.data(),upper_belt_RTau.data());

TSpline3 *s_lower_belt_RTau = new TSpline3("spline_lower_belt_RTau",g_lower_belt_RTau);
TSpline3 *s_upper_belt_RTau = new TSpline3("spline_upper_belt_RTau",g_upper_belt_RTau);

TCanvas *c_belt_RTau = new TCanvas("c_belt_RTau","");
c_belt_RTau->cd();
g_lower_belt_RTau->SetMarkerStyle(7);
g_lower_belt_RTau->SetMarkerSize(2);
g_upper_belt_RTau->SetMarkerStyle(7);
g_upper_belt_RTau->SetMarkerSize(2);
s_lower_belt_RTau->SetLineColor(kBlue);
s_upper_belt_RTau->SetLineColor(kBlue);
s_lower_belt_RTau->SetLineWidth(2);
s_upper_belt_RTau->SetLineWidth(2);
g_upper_belt_RTau->Draw("AP");
g_lower_belt_RTau->Draw("SAME");
s_lower_belt_RTau->Draw("SAME");
s_upper_belt_RTau->Draw("SAME");
g_upper_belt_RTau->SetMinimum(lower_belt_RTau.data()[0]*1.1);

c_belt_RTau->Write();


g_lower_belt_RTau->Write("g_lower_belt_RTau");
g_upper_belt_RTau->Write("g_upper_belt_RTau");
s_lower_belt_RTau->Write("s_lower_belt_RTau");
s_upper_belt_RTau->Write("s_upper_belt_RTau");

outfile->Close();
dummy->Close();

}
