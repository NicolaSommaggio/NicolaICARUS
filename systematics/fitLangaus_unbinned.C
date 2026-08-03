double X_LOW = 0;
double X_HIGH = 30;
int N_POINTS = 1000;
double BIN_WIDTH = 0.25;

std::vector<double> xvar;

double langaufun(double *x, double *par) {
 
   //Fit parameters:
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Width (sigma) of convoluted Gaussian function
   //
   //In the Landau distribution (represented by the CERNLIB approximation),
   //the maximum is located at x=-0.22278298 with the location parameter=0.
   //This shift is corrected within this function, so that the actual
   //maximum is identical to the MP parameter.
 
      // Numeric constants
      double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
      double mpshift  = -0.22278298;       // Landau maximum location
 
      // Control constants
      double np = 200.0;      // number of convolution steps
      double sc =   5.0;      // convolution extends to +-sc Gaussian sigmas
 
      // Variables
      double xx;
      double mpc;
      double fland;
      double sum = 0.0;
      double xlow,xupp;
      double step;
      double i;
 
 
      // MP shift correction
      mpc = par[1] - mpshift * par[0];
 
      // Range of convolution integral
      xlow = x[0] - sc * par[2];
      xupp = x[0] + sc * par[2];
 
      step = (xupp-xlow) / np;
 
      // Convolution integral of Landau and Gaussian by sum
      for(i=1.0; i<=np/2; i++) {
         xx = xlow + (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0],1) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[2],1);
 
         xx = xupp - (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0],1) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[2],1);
      }
 
      return (step * sum * invsq2pi / par[2]);
}

/*
double get_langaus_norm(double *par, double xmin, double xmax, int npoints)
{
    double step = (xmax - xmin) / npoints;
    double sum = 0;
    for (int i = 0; i <= npoints; ++i) 
    {
        double xx = xmin + i * step;
        sum += langaufun(&xx, par);
    }
    return step * sum;
}
*/

double get_langaus_norm(double *par, double xmin, double xmax, int npoints)
{
    if (npoints % 2 == 1) npoints++; // Simpson richiede n pari
    double step = (xmax - xmin) / npoints;
    double sum = langaufun(&xmin, par);
    double xend = xmax;
    sum += langaufun(&xend, par);
    for (int i = 1; i < npoints; ++i) {
        double xx = xmin + i * step;
        double val = langaufun(&xx, par);
        sum += (i % 2 == 0 ? 2.0 : 4.0) * val;
    }
    return step / 3.0 * sum;
}


double gauss_pdf(double x, double *par)
{
    const Double_t sqrt2 = TMath::Sqrt(2.0);

    // *** 1 GAUSSIAN ***
    double gnorm = 0.5 * (
            TMath::Erf((X_HIGH - par[0])/(sqrt2*par[1]))
          - TMath::Erf((X_LOW - par[0])/(sqrt2*par[1]))
        );
    return TMath::Gaus(x, par[0], par[1], 1) / gnorm;

    //return TMath::Gaus(x, par[0], par[1], 1);

    // *** 2 GAUSSIANS ***
    //double gnorm1 = 0.5 * (
    //        TMath::Erf((X_HIGH - par[0])/(sqrt2*par[1]))
    //      - TMath::Erf((X_LOW - par[0])/(sqrt2*par[1]))
    //);

    //double gnorm2 = 0.5 * (
    //        TMath::Erf((X_HIGH - par[2])/(sqrt2*par[3]))
    //      - TMath::Erf((X_LOW - par[2])/(sqrt2*par[3]))
    //);

    //return par[5] * TMath::Gaus(x, par[0], par[1], 1) / gnorm1 + (1-par[5]) * TMath::Gaus(x, par[2], par[3], 1) / gnorm2;
    

}

/*
double fit_pdf(double *x, double *par, double langau_norm)
{
    // *** LANDAU + 1 GAUSSIAN ***
    //double landau_par[3] = {par[0], par[1], par[2]};
    //double gauss_par[2]  = {par[3], par[4]};

    //double langau_val = langaufun(x, landau_par) / langau_norm;
    //double gauss_val  = gauss_pdf(x[0], gauss_par);

    //return par[5] * langau_val + (1-par[5]) * gauss_val;

    // *** LANDAU + 2 GAUSSIANS ***
    //double landau_par[3] = {par[0], par[1], par[2]};
    //double gauss_par[5]  = {par[3], par[4], par[5], par[6], par[7]};

    //double langau_val = langaufun(x, landau_par) / langau_norm;
    //double gauss_val  = gauss_pdf(x[0], gauss_par);

    //return par[8] * langau_val + (1-par[8]) * gauss_val;

}
*/

double fit_pdf(double *x, double *par, double langau_norm1, double langau_norm2)
{
    // *** LANDAU + 1 GAUSSIAN ***
    //double landau_par[3] = {par[0], par[1], par[2]};
    //double gauss_par[2]  = {par[3], par[4]};

    //double langau_val = langaufun(x, landau_par) / langau_norm;
    //double gauss_val  = gauss_pdf(x[0], gauss_par);

    //return par[5] * langau_val + (1-par[5]) * gauss_val;

    // *** LANDAU + 2 GAUSSIANS ***
    //double landau_par[3] = {par[0], par[1], par[2]};
    //double gauss_par[5]  = {par[3], par[4], par[5], par[6], par[7]};

    //double langau_val = langaufun(x, landau_par) / langau_norm;
    //double gauss_val  = gauss_pdf(x[0], gauss_par);

    //return par[8] * langau_val + (1-par[8]) * gauss_val;

    // *** 2 LANDAU
    double landau_par1[3] = {par[0], par[1], par[2]};
    double landau_par2[3]  = {par[3], par[4], par[5]};

    double langau_val1 = langaufun(x, landau_par1) / langau_norm1;
    double langau_val2 = langaufun(x, landau_par2) / langau_norm2;

    return par[6] * langau_val1 + (1-par[6]) * langau_val2;

}


double pdf_proj(double *x, double *par, double langau_norm1, double langau_norm2, int nmax, double bin_width)
{
    //return (nmax * bin_width) * fit_pdf(x, par, langau_norm);
    return fit_pdf(x, par, langau_norm1, langau_norm2);
    
}

double landau_proj(double *x, double *par, double langau_norm, int nmax, double bin_width)
{
    //return (nmax * bin_width) * langaufun(x, par) / langau_norm;
    return langaufun(x, par) / langau_norm;
}

double gauss_proj(double *x, double *par, int nmax, double bin_width)
{
    //return (nmax * bin_width) * gauss_pdf(x[0], par);
    return gauss_pdf(x[0], par);
}


double pdf_proj_light(double *x, double *par, double langau_norm1, double langau_norm2)
{
    //return (nmax * bin_width) * fit_pdf(x, par, langau_norm);
    return fit_pdf(x, par, langau_norm1, langau_norm2);
    
}

double landau_proj_light(double *x, double *par, double langau_norm)
{
    //return (nmax * bin_width) * langaufun(x, par) / langau_norm;
    return langaufun(x, par) / langau_norm;
}

double gauss_proj_light(double *x, double *par)
{
    //return (nmax * bin_width) * gauss_pdf(x[0], par);
    return gauss_pdf(x[0], par);
}

/*
void fcn(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t)
{
    double landau_par[3] = {par[0], par[1], par[2]};
    double langau_norm = get_langaus_norm(landau_par, X_LOW, X_HIGH, N_POINTS);

    double Like = 0;
    for (size_t i = 0; i < xvar.size(); ++i) {
        double p = fit_pdf(&xvar[i], par, langau_norm);
        if (p <= 0) { f = 1e30; return; }  // guard against log(0)/log(negative)
        Like += log(p);
    }
    f = -2.0 * Like;
}
*/

/*
double fitFunction(double *x, double *par)
{
    // par[0..2] = Landau (width, MPV, sigma conv.)
    // par[3..4] = Gauss (mean, sigma)
    // par[5]    = frazione Landau vs Gauss
    // par[6]    = scala (N_eventi * bin_width)

    double landau_par[3] = {par[0], par[1], par[2]};
    double langau_norm = get_langaus_norm(landau_par, X_LOW, X_HIGH, N_POINTS);

    return par[6] * BIN_WIDTH * fit_pdf(x, par, langau_norm);
}
*/

double fitFunction(double *x, double *par)
{
    static double last_par1[3] = {-1, -1, -1};
    static double cached_norm1 = 1.0;

    bool changed1 = (par[0] != last_par1[0] || par[1] != last_par1[1] || par[2] != last_par1[2]);
    if (changed1) {
        double landau_par1[3] = {par[0], par[1], par[2]};
        cached_norm1 = get_langaus_norm(landau_par1, X_LOW, X_HIGH, N_POINTS);
        last_par1[0] = par[0]; last_par1[1] = par[1]; last_par1[2] = par[2];
    }

    static double last_par2[3] = {-1, -1, -1};
    static double cached_norm2 = 1.0;

    bool changed2 = (par[3] != last_par2[0] || par[4] != last_par2[1] || par[5] != last_par2[2]);
    if (changed2) {
        double landau_par2[3] = {par[3], par[4], par[5]};
        cached_norm2 = get_langaus_norm(landau_par2, X_LOW, X_HIGH, N_POINTS);
        last_par2[0] = par[3]; last_par2[1] = par[4]; last_par2[2] = par[5];
    }

    //return par[6] * BIN_WIDTH * fit_pdf(x, par, cached_norm);
    return fit_pdf(x, par, cached_norm1, cached_norm2);
}


/*
        const int nparam = 6;
        TMinuit *my_gMinuit = new TMinuit(nparam);  //initialize TMinuit with a maximum of 5 params
        //gMinuit->SetPrintLevel(-1);
        my_gMinuit->SetFCN(fcn);      // set the FCN
        
        Double_t arglist[2];
        Int_t ierflg = 0;  // Error return code: 0 if the command was correctly executed, >0 otherwise. 
        
        // Set starting values and step sizes for parameters
        Double_t vstart[nparam] = {0.1, dedx_expected, 0.1, 2.5, 1., 0.5};
        Double_t step[nparam]   = {0.001, 0.001, 0.001, 0.001, 0.001, 0.001};   //step 0 li rende costanti
     
        my_gMinuit->mnparm(0, "landau_width", vstart[0], step[0], 0., 0., ierflg);
        my_gMinuit->mnparm(1, "landau_mpv", vstart[1], step[1], 0., 0., ierflg);
        my_gMinuit->mnparm(2, "gauss_width", vstart[2], step[2], 0., 0., ierflg);
        my_gMinuit->mnparm(3, "interacting_gauss_mean", vstart[3], step[3], 0., 0., ierflg);
        my_gMinuit->mnparm(4, "interacting_gauss_width", vstart[4], step[4], 0., 0., ierflg);
        my_gMinuit->mnparm(5, "fraction", vstart[5], step[5], 0., 1., ierflg);

        arglist[0] = 500.;//500;
        arglist[1] = 0.1;

        my_gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
      
        // Print results
        Double_t amin,edm,errdef;
        Int_t nvpar,nparx,icstat;
        my_gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);     
        my_gMinuit->mnprin(3,amin);
   
        Double_t landau_width, landau_mpv, gauss_width, interacting_gauss_mean, interacting_gauss_width, fraction;
        Double_t e_landau_width, e_landau_mpv, e_gauss_width, e_interacting_gauss_mean, e_interacting_gauss_width, e_fraction;
        Int_t ivar=0;
        Double_t bnd1, bnd2;
        TString chnam;
        
        my_gMinuit->mnpout(0, chnam, landau_width, e_landau_width, bnd1, bnd2, ivar);
        my_gMinuit->mnpout(1, chnam, landau_mpv, e_landau_mpv, bnd1, bnd2, ivar);
        my_gMinuit->mnpout(2, chnam, gauss_width, e_gauss_width, bnd1, bnd2, ivar);
        my_gMinuit->mnpout(3, chnam, interacting_gauss_mean, e_interacting_gauss_mean, bnd1, bnd2, ivar);
        my_gMinuit->mnpout(4, chnam, interacting_gauss_width, e_interacting_gauss_width, bnd1, bnd2, ivar);
        my_gMinuit->mnpout(5, chnam, fraction, e_fraction, bnd1, bnd2, ivar);

        std::vector<double> x_points;
        std::vector<double> y_points;

        Double_t pars[nparam] = {landau_width, landau_mpv, gauss_width, interacting_gauss_mean, interacting_gauss_width, fraction};

        Double_t landau_pars[3] = {landau_width,landau_mpv,gauss_width};
        double landau_norm = get_langaus_norm(landau_pars,0,30,10000);

        for(int xi=0; xi<10000; xi++)
        {
            double xmin = 0;
            double xmax = 10;
        
            double x = xmin + xi * (xmax-xmin)/10000;

            x_points.push_back(x);
            y_points.push_back(pdf_proj(&x,pars,landau_norm,xvar.size(),h_dedx_PRO_MC->GetBinWidth(1)));
        }

        TGraph * fit_function_plot = new TGraph(10000,x_points.data(),y_points.data());

        TCanvas *c_fit = new TCanvas("c_fit", "fit");
        c_fit->cd();
        h_dedx_PRO_MC -> Draw();
        fit_function_plot -> Draw("SAME");

        outfile -> cd();

        c_fit ->Write();

        */

        /*
        TFile *dedx_template = TFile::Open("RefCurvesChi2.root","READ");
        TH1D *dedx_range_pro = (TH1D*)dedx_template->Get("dedx_range_pro");
        double dedx_expected = dedx_range_pro -> GetBinContent(dedx_range_pro -> FindBin(rr));

        BIN_WIDTH = h_dedx_PRO_MC->GetBinWidth(1);

        cout << "START FITTING WITH: RANGE = [" << X_LOW << " " << X_HIGH << "]" << " BIN WIDTH " << endl;

        TF1 *f1 = new TF1("f1", fitFunction, X_LOW, X_HIGH, 7);

        f1->SetParNames("LandauWidth", "LandauMPV", "ConvSigma",
                "GaussMean", "GaussSigma", "FracLandau", "Scale");

        // valori iniziali
        f1->SetParameters(0.5, dedx_expected, 0.3, 2.5, 1.0, 0.5, h_dedx_PRO_MC->Integral());

        // limiti sensati
        //f1->SetParLimits(0, 0., 0.);
        //f1->SetParLimits(1, 0., 0.);
        //f1->SetParLimits(2, 0., 0.);
        //f1->SetParLimits(3, 0., 0.);
        //f1->SetParLimits(4, 0., 0.);
        //f1->SetParLimits(5, 0., 0.);
        //f1->SetParLimits(6, 0., 0.);
        //f1->SetParLimits(7, 0., 0.);

        h_dedx_PRO_MC->Fit(f1,"L"); // L = binned log-likelihood, R = usa il range della TF1

        outfile -> cd();
        h_dedx_PRO_MC -> Write();
        */
        
