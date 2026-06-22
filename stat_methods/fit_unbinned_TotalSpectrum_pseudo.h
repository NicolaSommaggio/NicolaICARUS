Double_t _invariant_mass_pdf_no_norm(Double_t x, Double_t *par)
{
    Double_t frac   = par[0];
    Double_t mu1    = par[1];
    Double_t mu2    = par[2];
    Double_t sigma1 = par[3];
    Double_t sigma2 = par[4];

    return frac * TMath::Gaus(x, mu1, sigma1, kTRUE)
         + (1.0 - frac) * TMath::Gaus(x, mu2, sigma2, kTRUE);
}

Double_t _invariant_mass_pdf_norm(Double_t *par,
                                  Double_t xmin,
                                  Double_t xmax)
{
    Double_t frac   = par[0];
    Double_t mu1    = par[1];
    Double_t mu2    = par[2];
    Double_t sigma1 = par[3];
    Double_t sigma2 = par[4];

    const Double_t sqrt2 = TMath::Sqrt(2.0);

    auto gaussIntegral = [&](Double_t mu, Double_t sigma)
    {
        return 0.5 * (
            TMath::Erf((xmax - mu)/(sqrt2*sigma))
          - TMath::Erf((xmin - mu)/(sqrt2*sigma))
        );
    };

    return frac * gaussIntegral(mu1, sigma1)
         + (1.0 - frac) * gaussIntegral(mu2, sigma2);
}

Double_t _invariant_mass_pdf(Double_t x,
                             Double_t *par)
{
    Double_t norm = _invariant_mass_pdf_norm(par, TOT_MASS_LOW, TOT_MASS_HIGH);

    if (norm <= 0.0)
        return 0.0;

    return _invariant_mass_pdf_no_norm(x, par) / norm;
}


double _pdf_Argus ( double *x , double * par )
{
  const double m = x [0];
  //const double m0 = par [0];
  //const double c = par [1];
  //const double p = par [2];

  const double m0 = M0_FIT_DS_PLUS_PHI_MU_NU;
  const double c = par [0];
  const double p = par [1];

  if ( m >= m0 ) {
    return 0.0;
  }

  double L_taumass_range = TOT_MASS_LOW;
  double H_taumass_range = TOT_MASS_HIGH;

  const double xL = 1.0 - ( L_taumass_range / m0 ) *( L_taumass_range / m0 ) ;
  const double xH = 1.0 - ( H_taumass_range / m0 ) *( H_taumass_range / m0 ) ;

  const double gammaA = ROOT :: Math :: tgamma (1.0 + p ) ;
  const double dL = gammaA * ROOT :: Math :: inc_gamma_c (1.0 + p , -c * xL ) ;
  const double dH = gammaA * ROOT :: Math :: inc_gamma_c (1.0 + p , -c * xH ) ;
  const double norm = ( m0 * m0 ) /(2.0* c * std :: pow ( -c , p ) ) * ( dL - dH ) ;

  const double u = 1.0 - ( m / m0 ) *( m / m0 ) ;
  return m * std :: pow (u , p ) * std :: exp ( c * u ) / norm ;

}

Double_t _combinatorial_pdf(Double_t x, Double_t *par)
{
  Double_t tau = par[0];

  Double_t val = 1./tau * TMath::Exp(-1 * x * (1./tau)); 
  
  Double_t norm = - 1 * TMath::Exp(-1 * TOT_MASS_HIGH * (1./tau)) + TMath::Exp(-1 * TOT_MASS_LOW * (1./tau));

  return val / norm;
}


Double_t total_mass_spectrum_pdf(Double_t *xx, Double_t *par)
{
  Double_t x = xx[0];
  Double_t f_D_PhiPi = par[0]; 
  Double_t f_DS_PhiMuNu = par[1];
  Double_t f_DS_PhiPi = par[2];
  Double_t f_DS_TauNu = par[3];
  //Double_t tau = par[4];

  Double_t par_D_PhiPi[5] = {
                             FRACTION_FIT_D_PLUS_PHI_PI,
                             MU1_FIT_D_PLUS_PHI_PI, 
                             MU2_FIT_D_PLUS_PHI_PI,
                             SIGMA1_FIT_D_PLUS_PHI_PI,
                             SIGMA2_FIT_D_PLUS_PHI_PI
                            };

  Double_t par_DS_PhiMuNu[2] = {
                                // M0_FIT_DS_PLUS_PHI_MU_NU,
                                C_FIT_DS_PLUS_PHI_MU_NU,
                                P_FIT_DS_PLUS_PHI_MU_NU
                              };

  Double_t par_DS_PhiPi[5] = {
                              FRACTION_FIT_DS_PLUS_PHI_PI,
                              MU1_FIT_DS_PLUS_PHI_PI,
                              MU2_FIT_DS_PLUS_PHI_PI,
                              SIGMA1_FIT_DS_PLUS_PHI_PI,
                              SIGMA2_FIT_DS_PLUS_PHI_PI
                            };

  Double_t par_DS_TauNu[5] = {
                              FRACTION_FIT_MASS_DS_PLUS_TAU_NU,
                              MU1_FIT_MASS_DS_PLUS_TAU_NU,
                              MU2_FIT_MASS_DS_PLUS_TAU_NU,
                              SIGMA1_FIT_MASS_DS_PLUS_TAU_NU,
                              SIGMA2_FIT_MASS_DS_PLUS_TAU_NU
                            };

  Double_t val = f_DS_TauNu * _invariant_mass_pdf(x,par_DS_TauNu) + (1-f_DS_TauNu) *
                 (
                    f_D_PhiPi * _invariant_mass_pdf(x,par_D_PhiPi) +
                    f_DS_PhiMuNu * _pdf_Argus(&x,par_DS_PhiMuNu) +
                    f_DS_PhiPi * _invariant_mass_pdf(x,par_DS_PhiPi) +
                    (1 - f_D_PhiPi - f_DS_PhiMuNu - f_DS_PhiPi) * _combinatorial_pdf(x,&TAU_FIT_COMBINATORIAL)
                 );
  return val;

}

Double_t pdf_proj_total_mass_spectrum(Double_t *x, Double_t *par, Int_t max, Double_t bin_width)
{
  return (max*bin_width)*total_mass_spectrum_pdf(x,par);
}

void fcn_total_mass_spectrum(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t )
{
  Double_t Like=0;

  for(size_t i=0; i<x_var.size(); ++i)
  {
      Double_t p = total_mass_spectrum_pdf(&x_var[i],par);
      Like += TMath::Log(p);
  };

  f= - 2. * Like;


}

std::tuple<Double_t,Double_t,Double_t,Double_t> fit_unbinned_TotalSpectrum(std::vector<double> input_xvar, std::vector<double> vstart, std::vector<double> step, std::vector<double> low, std::vector<double> high,  std::vector<std::string> par_names, TH1D *h_time, TFile *file, std::string name, double xmin, double xmax)
{
     x_var = input_xvar;

     const int nparam = 4;


     const int n_fit = x_var.size();

      TMinuit *my_gMinuit = new TMinuit(nparam);  //initialize TMinuit with a maximum of 5 params

      my_gMinuit->SetPrintLevel(-1);

      my_gMinuit->SetFCN(fcn_total_mass_spectrum);      // set the FCN
        
      Double_t arglist[2];
      Int_t ierflg = 0;  // Error return code: 0 if the command was correctly executed, >0 otherwise. 
     
      for(int par = 0; par < nparam; par ++)
      {
        my_gMinuit->mnparm(par, par_names[par].c_str(), vstart[par], step[par], low[par], high[par], ierflg);
      }

      arglist[0] = 5000.;//500;
      arglist[1] = 0.1;

      my_gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
      my_gMinuit->mnexcm("HESSE", arglist ,0, ierflg);

      // Print results
      Double_t amin,edm,errdef;
      Int_t nvpar,nparx,icstat;
      my_gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);     

      if(!(icstat == 3 && edm < 1e-3)){ cout << "FIT NOT CONVERGED" << endl; }

      for (int par = 0; par < nparam; par++)
      {
          TString pname;
          Double_t val, err, lowb, upb;
          Int_t ivar;

          my_gMinuit->mnpout(par, pname, val, err, lowb, upb, ivar);
      }

      
      Double_t pars[nparam];
      Double_t pars_errors[nparam];

      Int_t ivar=0;
      TString chnam;
      
      for(int par = 0; par < nparam; par++)
      {
          my_gMinuit->mnpout(par, chnam, pars[par], pars_errors[par], low[par], high[par], ivar);
      }

      /*
      
      std::vector<double> x_points;
      std::vector<double> y_points;

      for(int xi=0; xi<10000; xi++)
      { 
            double x = xmin + xi * (xmax-xmin)/10000;

            //if(TMath::Abs(x - 1.77699) <= 3 * 0.00581298)continue; //BLINDING

            x_points.push_back(x);
            y_points.push_back(pdf_proj_total_mass_spectrum(&x,pars,n_fit,h_time->GetBinWidth(1)));
      }

      TGraph * fit_function_plot = new TGraph(x_points.size(),x_points.data(),y_points.data());

      TCanvas *c_fit = new TCanvas(Form("fit_%s",name.c_str()), "");
      c_fit->cd();

      const float split = 0.3;

      // -------------------- TOP PAD --------------------
      TPad *pad_fit = new TPad("pad_fit", "", 0., split, 1., 1.);
      pad_fit->SetBottomMargin(0.0);
      pad_fit->Draw();
      pad_fit->cd();

      gStyle->SetOptStat(0);
      gPad->SetTicks(1,1);

      // //h_time->GetXaxis()->SetRangeUser(0., 10.0);
      // h_time->SetLineColor(kBlack);
      // h_time->SetMarkerStyle(8);
      // h_time->SetMarkerSize(0.5);
      // h_time ->GetYaxis() -> SetTitleOffset(0.85);
      // //h_time->GetYaxis()->SetTitle("counts / 0.1 [fs / 410.3]");
      // h_time->GetYaxis()->SetTitleSize(0.06);
      // h_time->GetYaxis()->SetLabelSize(0.05);
      // h_time->Draw("PE");

      TH1D *h_plot = (TH1D*)h_time->Clone("h_plot");

      for (int i = 1; i <= h_plot->GetNbinsX(); i++)
      {
          double x = h_plot->GetBinCenter(i);

          //if (fabs(x - 1.77699) <= 3 * 0.00581298)
          //{
              //h_plot->SetBinContent(i, 0);
              //h_plot->SetBinError(i, 0);
          //}
      }
      h_plot->SetLineColor(kBlack);
      h_plot->SetMarkerStyle(8);
      h_plot->SetMarkerSize(0.5);
      h_plot ->GetYaxis() -> SetTitleOffset(0.85);
      h_plot->GetYaxis()->SetTitleSize(0.06);
      h_plot->GetYaxis()->SetLabelSize(0.05);
      h_plot->Draw("PE");
      
      // overlay fit

      fit_function_plot->SetLineColor(kRed);
      fit_function_plot->SetMarkerColor(kRed);
      fit_function_plot->SetLineWidth(2);
      fit_function_plot->Draw("L SAME");

      // IMPORTANT: go back to canvas
      c_fit->cd();

      // -------------------- BOTTOM PAD (PULLS) --------------------
      TPad *pad_pull = new TPad("pad_pull", "", 0., 0., 1., split);
      pad_pull->SetTopMargin(0.0);
      pad_pull->SetBottomMargin(0.3);
      pad_pull->Draw();
      pad_pull->cd();
        
      int n = h_time->GetNbinsX();

      TH1D *h_pull = new TH1D(Form("h_pull_%s",name.c_str()), "", 20, -5, 5);      
        
      TGraphErrors *pulls = new TGraphErrors(n);
      pulls->SetMarkerStyle(7);
        
      //pulls->GetXaxis()->SetTitle("t / 410.3 fs");
      pulls->GetYaxis()->SetTitle("Pull");
        
      int ip = 0;

      for (int i = 0; i < n; i++)
      {
          double x  = h_time->GetBinCenter(i+1);

          //if (fabs(x - 1.77699) <= 3 * 0.00581298) continue;

          double y  = h_time->GetBinContent(i+1);
          double ey = h_time->GetBinError(i+1);
      
          double y_fit = pdf_proj_total_mass_spectrum(&x, pars, n_fit, h_time->GetBinWidth(1));
      
          // avoid division by zero
          double pull = 0;
          if (ey > 0) pull = (y - y_fit) / ey;
      
          pulls->SetPoint(ip, x, pull);
          pulls->SetPointError(ip, 0., 1.0);
          ip++;
          h_pull->Fill(pull);

      }
      pulls->Set(ip);

      double xmin_pull = h_time->GetXaxis()->GetXmin();
      double xmax_pull = h_time->GetXaxis()->GetXmax();

      
      pulls -> SetMarkerStyle(8);
      pulls -> SetTitle("");
      pulls -> SetMarkerSize(0.4);
      pulls -> GetYaxis() -> SetTitleOffset(0.3);
      //pulls -> GetXaxis()-> SetTitle("t [fs / 410.3]");
      pulls -> GetYaxis()-> SetTitle("pulls");
      pulls -> GetYaxis()-> SetTitleSize(0.13);
      pulls -> GetXaxis()-> SetTitleSize(0.13);
      pulls -> GetYaxis()-> SetLabelSize(0.11);
      pulls -> GetXaxis()-> SetLabelSize(0.11);
      pulls->Draw("APE");
            
      // optional: axis range for visibility
      pulls->SetMinimum(-5);
      pulls->SetMaximum(5);
      pulls->GetXaxis()->SetLimits(xmin_pull, xmax_pull);



      // baseline at 0
      TLine *line_zero = new TLine(xmin_pull, 0., xmax_pull, 0.);
      
      line_zero->SetLineStyle(2);
      line_zero->Draw("same");
      // final update

      c_fit -> Update();


      file -> cd();
      c_fit -> Write();
      //c_fit -> SaveAs("PLOT_REPORT/FIT_BKG_TIME.pdf");
      */

      std::tuple<Double_t,Double_t,Double_t,Double_t> returned_vals = {pars[3],pars_errors[3],pars[1],pars_errors[1]};
      return returned_vals;
}