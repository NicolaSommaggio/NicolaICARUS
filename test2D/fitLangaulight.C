
double FindXAtY(TH1D* h, double yTarget, double mpv, bool searchRight) {
    int binMPV = h->FindBin(mpv);
    int step = searchRight ? 1 : -1;
    int bin = binMPV;

    while (bin > 1 && bin < h->GetNbinsX()) {
        double yVal = h->GetBinContent(bin);
        if (yVal <= yTarget) return h->GetBinCenter(bin);
        bin += step;
    }
    return h->GetBinCenter(bin);
}

std::vector<double> FWHM(TF1 *f, double mpv) {

double y_half = f->Eval(mpv) / 2.0;
double x1 = f->GetX(y_half, f->GetXmin(), mpv);
double x2 = f->GetX(y_half, mpv, f->GetXmax());
double width_half_max = x2 - x1;

std::vector<double> vector_active;
vector_active.push_back(x1);
vector_active.push_back(x2);
vector_active.push_back(width_half_max);
vector_active.push_back(y_half);

return vector_active;
//assumendo la distribuzione come gaussiana, sigma = width_half_max/2.354

}

/*
class RooDoubleSidedCB : public RooAbsPdf {
public:
  RooDoubleSidedCB(const char* name, const char* title,
                   RooAbsReal& _x,
                   RooAbsReal& _mu, RooAbsReal& _width,
                   RooAbsReal& _a1, RooAbsReal& _p1,
                   RooAbsReal& _a2, RooAbsReal& _p2)
    : RooAbsPdf(name, title),
      x("x", "x", this, _x),
      mu("mu", "mu", this, _mu),
      width("width", "width", this, _width),
      a1("a1", "a1", this, _a1),
      p1("p1", "p1", this, _p1),
      a2("a2", "a2", this, _a2),
      p2("p2", "p2", this, _p2) {}

  RooDoubleSidedCB() {}

  RooDoubleSidedCB(const RooDoubleSidedCB& other, const char* name=0)
    : RooAbsPdf(other, name),
      x("x", this, other.x),
      mu("mu", this, other.mu),
      width("width", this, other.width),
      a1("a1", this, other.a1),
      p1("p1", this, other.p1),
      a2("a2", this, other.a2),
      p2("p2", this, other.p2) {}

  virtual TObject* clone(const char* newname) const override {
    return new RooDoubleSidedCB(*this, newname);
  }

  inline virtual ~RooDoubleSidedCB() {}

protected:
  RooRealProxy x, mu, width, a1, p1, a2, p2;

  Double_t evaluate() const override {
    double u   = (x - mu)/width;
    double A1  = pow(p1/fabs(a1), p1) * exp(-0.5 * a1*a1);
    double B1  = p1/fabs(a1) - fabs(a1);
    double A2  = pow(p2/fabs(a2), p2) * exp(-0.5 * a2*a2);
    double B2  = p2/fabs(a2) - fabs(a2);

    if (u < -a1) return A1 * pow(B1 - u, -p1);
    else if (u < a2) return exp(-0.5 * u * u);
    else return A2 * pow(B2 + u, -p2);
  }
};
*/

/*
Double_t DoubleSidedCB2(double x, double mu, double width, double a1, double p1, double a2, double p2)
{
  double u   = (x-mu)/width;
  double A1  = TMath::Power(p1/TMath::Abs(a1),p1)*TMath::Exp(-a1*a1/2);
  double A2  = TMath::Power(p2/TMath::Abs(a2),p2)*TMath::Exp(-a2*a2/2);
  double B1  = p1/TMath::Abs(a1) - TMath::Abs(a1);
  double B2  = p2/TMath::Abs(a2) - TMath::Abs(a2);

  double result(1);
  if      (u<-a1) result *= A1*TMath::Power(B1-u,-p1);
  else if (u<a2)  result *= TMath::Exp(-u*u/2);
  else            result *= A2*TMath::Power(B2+u,-p2);
  return result;
}


double DoubleSidedCB(double* x, double *par)
{
  return(par[0] * DoubleSidedCB2(x[0], par[1],par[2],par[3],par[4],par[5],par[6]));
}
*/


std::vector<double> fitLangau(TH1D* h, TFile *fitPrintFile, const char * fitname, std::string fitPrintDirectory ){

    
        fitPrintFile->cd();
        TDirectory *d = (TDirectory*)fitPrintFile->Get(fitPrintDirectory.c_str()); 
        d->cd();

        std::vector<double> fitParameters;
        RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

        //TFile *FitOutfile = TFile::Open(filename, "UPDATE");

        int binMax=0;
        double maxVal=0;
        for(int j=1; j<h->GetNbinsX(); j++)
        {
            double y=h->GetBinContent(j);
            if(y>maxVal)
            {
                maxVal=y;
                binMax=j;
            }
        }

        Double_t startMpv=h->GetBinCenter(binMax), startSigmaL = 0.1 , sigmaGaus = h->GetStdDev();

        RooRealVar x("x", "x", h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
        RooDataHist dataHist("dataHist", "Dataset from TH1", x, RooFit::Import(*h));

        //LANDAU
        RooRealVar mpv("mpv", "Most Probable Value", startMpv, startMpv-5., startMpv+5.);
        RooRealVar sigmaL("sigmaL", "Landau width", startSigmaL, 0.02, 5.5);
        RooLandau landau("landau", "Landau PDF", x, mpv, sigmaL);

        //MEAN 0 GAUSSIAN (SMEARING)
        RooRealVar meanG("meanG", "Mean of Gaussian", 0);
        RooRealVar sigmaG("sigmaG", "Gaussian sigma", sigmaGaus, 0.1, 4.5);
        RooGaussian gauss("gauss", "Gaussian PDF", x, meanG, sigmaG);

        //LANDAU X GAUSSIAN CONVOLUTION
        RooFFTConvPdf langau("langau", "Landau convoluted with Gaussian", x, landau, gauss);

        x.setRange("fitRange", 1.5, 6.); 

        RooFitResult* fitRes = langau.fitTo(dataHist, RooFit::Save(), RooFit::Range("fitRange"), RooFit::PrintLevel(-1), RooFit::PrintEvalErrors(-1), RooFit::Warnings(false), RooFit::Verbose(false));
        

        // Save parameters
        fitParameters.push_back(mpv.getVal());
        fitParameters.push_back(mpv.getError());
        fitParameters.push_back(sigmaG.getVal());
        fitParameters.push_back(sigmaG.getError());
        fitParameters.push_back(sigmaL.getVal());
        fitParameters.push_back(sigmaL.getError());
        
        RooPlot* frame = x.frame();

        TCanvas* cFit = new TCanvas(Form("fit_%s",fitname));
        cFit->Divide(1,2);
        cFit->cd(1); 

        dataHist.plotOn(frame, RooFit::Name("data"));
 
        langau.plotOn(frame, RooFit::Name("fit"), RooFit::Range("fitRange"), RooFit::NormRange("fitRange"), RooFit::LineColor(kBlue));
        langau.paramOn(frame, RooFit::Layout(0.4, 0.9, 0.9));

        frame->Draw();

        cFit->cd(2);
        frame->GetXaxis()->SetTitle(""); 

        RooHist *hresid = frame->residHist("data", "fit");
        hresid->GetYaxis()->SetTitle("Data - Fit");
        hresid->GetXaxis()->SetTitle("dE/dx");
        hresid->Draw("AP"); 

        
        RooChi2Var chi2("chi2", "chi2", langau, dataHist);
        int ndf = dataHist.numEntries() - langau.getParameters(dataHist)->getSize();
        cout << "chi2/ndof (ndof=Nbins) " << chi2.getVal()/ndf << endl;
        

        fitParameters.push_back(chi2.getVal()/ndf);
        
        cFit->Write(fitname, TObject::kOverwrite);
            
        delete cFit;
        return fitParameters;

}


struct FitOutput
{
    std::vector<double> fitResults;
    TF1* fitFunction= nullptr;
    std::map<std::string, TF1*> componentsFunctions;
    double fraction;
};


FitOutput fitLangaulight(TH1D* h, double x_low, double x_high, TDirectory *fitPrintFile, std::string modelChoice, std::string name){

    
        fitPrintFile->cd();
        FitOutput FitOut;
        std::vector<double> fitParameters;
        RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

        //TFile *FitOutfile = TFile::Open(filename, "UPDATE");

        int binMax=0;
        double maxVal=0;
        for(int j=1; j<h->GetNbinsX(); j++)
        {
            double y=h->GetBinContent(j);
            if(y>maxVal)
            {
                maxVal=y;
                binMax=j;
            }
        }

        Double_t startMpv=h->GetBinCenter(binMax), startSigmaL = 0.1 , sigmaGaus = h->GetStdDev();

        RooRealVar x("x", "x", h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
        RooDataHist dataHist("dataHist", "Dataset from TH1", x, RooFit::Import(*h));

        //LANDAU
        RooRealVar mpv("mpv", "Most Probable Value", startMpv, startMpv-5., startMpv+5.);
        RooRealVar sigmaL("sigmaL", "Landau width", startSigmaL, 0.02, 5.5);
        RooLandau landau("landau", "Landau PDF", x, mpv, sigmaL);

        //MEAN 0 GAUSSIAN (SMEARING)
        RooRealVar meanG("meanG", "Mean of Gaussian", 0);
        RooRealVar sigmaG("sigmaG", "Gaussian sigma", sigmaGaus, 0.1, 4.5);
        RooGaussian gauss("gauss", "Gaussian PDF", x, meanG, sigmaG);

        //LANDAU X GAUSSIAN CONVOLUTION
        RooFFTConvPdf langau("langau", "Landau convoluted with Gaussian", x, landau, gauss);

        Double_t startSigma=0.2, starta1=1., startp1=1., starta2=1., startp2=1.;
        
        //DOUBLE SIDE CRYSTAL BALL
        RooRealVar mpv_dscb("mpv_dscb", "Most Probable Value", startMpv, startMpv-5., startMpv+5);
        RooRealVar sigmaLeft("sigmaLeft", "Gaussian Lsigma dscb", startSigma, 0.25, 4.);
        RooRealVar sigmaRight("sigmaRight", "Gaussain Rsigma dscb", startSigma, 0.25, 4. );
        RooRealVar a1("a1", "a1", starta1, -5., 5.);
        RooRealVar p1("p1", "p1", startp1, 0.,5.);
        RooRealVar a2("a2", "a2", starta2, -5., 5.);
        RooRealVar p2("p2", "p2", startp2, 0.,5.);

        RooRealVar p("p","p", 1, 0.01 , 10.);
        RooRealVar alpha("alpha", "alpha", -1., -4., -0.);

        RooCrystalBall dcbPdf("dcbPdf","DoubleSidedCB",x,mpv_dscb,sigmaLeft, sigmaRight,a1,p1,a2,p2);
        //RooDoubleSidedCB dcbPdf("dcbPdf", "Custom DSCB", x, mpv_dscb,sigma,a1,p1,a2,p2);
        RooCrystalBall CBPdf("CBPdf", "Single Side Crystal Ball", x,mpv_dscb,sigmaLeft,alpha,p,false);
        RooFFTConvPdf CBPdfConv("CBPdfConv", "CBPdf x gauss", x, CBPdf, gauss);

        //PASSING MUON CONTRIBUTION
        RooRealVar meanMIP("mean_mip", "Most Probable Value for passing", 1.8, 0.5, 5.);
        RooRealVar sigmaMIP("sigmaMIP", " width for passing", 1., 0., 10.);
        RooGaussian curveMIP("curveMIP", "curve PDF for passing", x, meanMIP, sigmaMIP);
        
        RooRealVar fraction("fraction", "fraction of stopping", 0.5, 0. ,1.);

        //DOUBLE SIDE CRYSTAL BALL X GAUSSIAN CONVOLUTION
        RooFFTConvPdf model("model", "dcbpdf x gauss", x, dcbPdf, gauss);

        //DOUBLE SIDE CRYSTAL BALL X GAUSSIAN CONVOLUTION PLUS PASSING CONTRIBUTION
        RooAddPdf modelplusPassing("model_plus_passing", "dcbpdf x gauss + Passing", RooArgList(model, curveMIP), RooArgList(fraction) );

        x.setRange("fitRange", x_low, x_high); 

        // Select model
        RooAbsPdf* chosenModel = nullptr;
        if (modelChoice == "dscb") 
        {
            chosenModel = &dcbPdf;
        } 
        else if (modelChoice=="CB")
        {
            chosenModel = &CBPdfConv;
        }
        else if (modelChoice == "langau") 
        {
            chosenModel = &langau;
        }
        else if (modelChoice== "model_plus_passing") 
        {
            chosenModel = &modelplusPassing;
        }

        x.setRange("fitRange", x_low, x_high);

        RooFitResult* fitRes = chosenModel->fitTo(dataHist, RooFit::Save(), RooFit::Range("fitRange"), RooFit::PrintLevel(-1), RooFit::PrintEvalErrors(-1), RooFit::Warnings(false), RooFit::Verbose(false));
        
        //TF1* tf1_langau = chosenModel->asTF(RooArgList(x), RooArgList());
        //std::vector<double> vec_larghezza = FWHM(tf1_langau,tf1_langau->GetMaximumX(x_low, x_high));
        //fitParameters.push_back(vec_larghezza[2]);

        // Save parameters depending on model
        if (modelChoice == "dscb" || modelChoice == "model_plus_passing" || modelChoice=="CB") 
        {
            fitParameters.push_back(mpv_dscb.getVal());
            fitParameters.push_back(mpv_dscb.getError());
            fitParameters.push_back(sigmaG.getVal());
            fitParameters.push_back(sigmaG.getError());
            fitParameters.push_back(sigmaLeft.getVal());
            fitParameters.push_back(sigmaLeft.getError());
            fitParameters.push_back(startMpv);
        } 
        else if (modelChoice == "langau") 
        {
            fitParameters.push_back(mpv.getVal());
            fitParameters.push_back(mpv.getError());
            fitParameters.push_back(sigmaG.getVal());
            fitParameters.push_back(sigmaG.getError());
            fitParameters.push_back(sigmaL.getVal());
            fitParameters.push_back(sigmaL.getError());
            fitParameters.push_back(startMpv);
        }


        RooPlot* frame = x.frame();

        TCanvas* cFit = new TCanvas(Form("fit_%s",name.c_str()));
        cFit->Divide(1,2);
        cFit->cd(1); 

        dataHist.plotOn(frame, RooFit::Name("data"));
 
        chosenModel->plotOn(frame, RooFit::Name("fit"), RooFit::Range("fitRange"), RooFit::NormRange("fitRange"), RooFit::LineColor(kBlue));
        chosenModel->paramOn(frame, RooFit::Layout(0.4, 0.9, 0.9));

        frame->Draw();

        cFit->cd(2);
        frame->GetXaxis()->SetTitle(""); 

        RooHist *hresid = frame->residHist("data", "fit");
        hresid->GetYaxis()->SetTitle("Data - Fit");
        hresid->GetXaxis()->SetTitle("dE/dx");
        hresid->Draw("AP"); 

        
        RooChi2Var chi2("chi2", "chi2", *chosenModel, dataHist);
        int ndf = dataHist.numEntries() - chosenModel->getParameters(dataHist)->getSize();
        //cout << "chi2/ndof (ndof=Nbins)" << chi2.getVal()/ndf << endl;
        

        fitParameters.push_back(chi2.getVal()/ndf);
        
        cFit->Write(0, TObject::kOverwrite);
        //cFit->Write();
    

        delete cFit;

        FitOut.fitResults = fitParameters;
        FitOut.fitFunction = nullptr; // fitLangaulight non crea tf1

        return FitOut;

}




FitOutput fitLangaulight_unbinned(TH1D* h, std::vector<double> v, double x_low, double x_high, TDirectory *fitPrintFile, std::string modelChoice , int nbin){

    
        fitPrintFile->cd();
        FitOutput FitOut;
        std::vector<double> fitParameters;
        RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

        //TFile *FitOutfile = TFile::Open(filename, "UPDATE");

        int binMax=0;
        double maxVal=0;
        for(int j=1; j<h->GetNbinsX(); j++)
        {
            double y=h->GetBinContent(j);
            if(y>maxVal)
            {
                maxVal=y;
                binMax=j;
            }
        }

        Double_t startMpv=h->GetBinCenter(binMax), startSigmaL = 0.4 , sigmaGaus = h->GetStdDev();

        RooRealVar x("x", "x", h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
        //RooDataHist dataHist("dataHist", "Dataset from TH1", x, RooFit::Import(*h));
        RooDataSet data("data", "Imported Data", RooArgSet(x));

        for(int i =0; i<int(v.size()); i++)
        {

            x.setVal(v[i]);
            data.add(RooArgSet(x));
        }
 
        //LANDAU
        RooRealVar mpv("mpv", "Most Probable Value", startMpv, startMpv-5., startMpv+5.);
        RooRealVar sigmaL("sigmaL", "Landau width", startSigmaL, 0.01, 1.5);
        RooLandau landau("landau", "Landau PDF", x, mpv, sigmaL);

        //MEAN 0 GAUSSIAN (SMEARING)
        RooRealVar meanG("meanG", "Mean of Gaussian", 0);
        RooRealVar sigmaG("sigmaG", "Gaussian sigma", sigmaGaus, 0.1, 4.5);
        RooGaussian gauss("gauss", "Gaussian PDF", x, meanG, sigmaG);

        //LANDAU X GAUSSIAN CONVOLUTION
        RooFFTConvPdf langau("langau", "Landau convoluted with Gaussian", x, landau, gauss);

        //ADDITIONAL LANDAU
        RooRealVar mpv_addL("mpv_addL", "Most Probable Value", 2.8, 2.65, 3.5);
        //RooFormulaVar mpv_addL("mpv_addL", "2*@0", RooArgList(mpv));
        RooRealVar sigmaL_addL("sigmaL_addL", "Landau width", 0.1, 0.01, 2.);
        RooLandau landau_add("landau_add", "Landau PDF", x, mpv_addL, sigmaL_addL);

        //ADDITIONAL MEAN 0 GAUSSIAN (SMEARING)
        RooRealVar meanG_addG("meanG_addG", "Mean of Gaussian", 0);
        RooRealVar sigmaG_addG("sigmaG_addG", "Gaussian sigma", 0.2, 0.05, 2.);
        RooGaussian gauss_add("gauss_add", "Gaussian PDF", x, meanG_addG, sigmaG_addG);

        //ADDITIONAL LANDAU X GAUSSIAN COVOLUTION
        RooFFTConvPdf langau_add("langau_add", "Additional Landau convoluted with Gaussian", x, landau_add, gauss_add);
        
        RooRealVar fraction("fraction", "fraction of main Landau", 0.9, 0.85 ,0.96);

        RooAddPdf double_langau("double_langau", "double_langau", RooArgList(langau, langau_add), RooArgList(fraction) );


        Double_t startSigma=0.3, starta1=1., startp1=1., starta2=1., startp2=1.;
        //DOUBLE SIDE CRYSTAL BALL
        RooRealVar mpv_dscb("mpv_dscb", "Most Probable Value", startMpv, startMpv-5., startMpv+5);
        RooRealVar sigma("sigma", "Gaussian sigma dscb", startSigma, 0.05, 5.);
        RooRealVar a1("a1", "a1", starta1, -10.,10.);
        RooRealVar p1("p1", "p1", startp1, -10, 10.);
        RooRealVar a2("a2", "a2", starta2, -10.,10.);
        RooRealVar p2("p2", "p2", startp2, -10.,10.);
        RooCrystalBall dcbPdf("dcbPdf","DoubleSidedCB",x,mpv_dscb,sigma,a1,p1,a2,p2);

        //DOUBLE SIDE CRYSTAL BALL X GAUSSIAN CONVOLUTION
        RooFFTConvPdf dscbxgauss("dscbxgauss", "dcbpdf x gauss", x, dcbPdf, gauss);

        x.setRange("fitRange", x_low, x_high); 

        // Select model
        RooAbsPdf* chosenModel = nullptr;
        if (modelChoice == "dscb") {chosenModel = &dcbPdf;} 
        else if (modelChoice == "langau") {chosenModel = &langau;}
        else if (modelChoice== "dscbxgauss") {chosenModel = &dscbxgauss;}
        else if (modelChoice == "double_langau") {chosenModel = &double_langau;}


        x.setRange("fitRange", x_low, x_high);

        RooFitResult* fitRes = chosenModel->fitTo(data, RooFit::Save(), RooFit::Range("fitRange"), RooFit::PrintLevel(-1), RooFit::PrintEvalErrors(-1), RooFit::Warnings(false), RooFit::Verbose(false));

        std::vector<double> xvals;
        std::vector<double> yvals;
        //ofstream csmoothness(Form("csmoothness_%d.txt", nbin));

        TF1* tf1;
        std::map<std::string, TF1*> components;
        double frac;
        if(modelChoice == "dscb")tf1 = chosenModel->asTF(x, RooArgList(mpv_dscb, sigma, a1, p1, a2, p2));
        else if(modelChoice == "langau")tf1 = chosenModel->asTF(x, RooArgList(mpv, sigmaL, meanG, sigmaG));
        else if(modelChoice == "double_langau")
        {
            tf1 = chosenModel->asTF(x, RooArgList(mpv, sigmaL, meanG, sigmaG, mpv_addL, sigmaL_addL, meanG_addG, sigmaG_addG)); 
            components["main_landau"] = langau.asTF(x, RooArgList(mpv, sigmaL, meanG, sigmaG));
            components["background_landau"] = langau_add.asTF(x, RooArgList(mpv_addL, sigmaL_addL, meanG_addG, sigmaG_addG));
            frac = fraction.getVal(); 
        }
        else if(modelChoice == "dscbxgauss")tf1 = chosenModel->asTF(x, RooArgList(mpv_dscb, sigma, a1, p1, a2, p2, meanG, sigmaG, gauss));
        tf1->SetNpx(10000);


        

        /*
        for (double val = h->GetXaxis()->GetXmin(); val <= h->GetXaxis()->GetXmax(); val += 0.0001) 
        {
            x.setVal(val);
            xvals.push_back(val);
            yvals.push_back(chosenModel->getVal());
            csmoothness << val << " " << chosenModel->getVal() << endl;
        }
        TSpline3* spline = new TSpline3(Form("spline%d.txt", nbin), xvals.data(), yvals.data(), 300000);
        */

        // Save parameters depending on model
        if (modelChoice == "dscb" || modelChoice == "dscbxgauss" ) 
        {
            fitParameters.push_back(mpv_dscb.getVal());
            fitParameters.push_back(mpv_dscb.getError());
            fitParameters.push_back(sigmaG.getVal());
            fitParameters.push_back(sigmaG.getError());
            fitParameters.push_back(sigma.getVal());
            fitParameters.push_back(sigma.getError());
            fitParameters.push_back(startMpv);
        } 
        else if (modelChoice == "langau") 
        {
            fitParameters.push_back(mpv.getVal()); //0
            fitParameters.push_back(mpv.getError()); //1
            fitParameters.push_back(sigmaG.getVal()); //2
            fitParameters.push_back(sigmaG.getError()); //3
            fitParameters.push_back(sigmaL.getVal()); //4
            fitParameters.push_back(sigmaL.getError()); //5 
            fitParameters.push_back(startMpv); //6 
        }
        else if (modelChoice == "double_langau") 
        {
            fitParameters.push_back(mpv.getVal()); //0
            fitParameters.push_back(mpv.getError()); //1
            fitParameters.push_back(sigmaG.getVal()); //2
            fitParameters.push_back(sigmaG.getError()); //3
            fitParameters.push_back(sigmaL.getVal()); //4
            fitParameters.push_back(sigmaL.getError()); //5 
            fitParameters.push_back(sigmaG_addG.getVal()); //6
            fitParameters.push_back(sigmaG_addG.getError()); //7
            fitParameters.push_back(sigmaL_addL.getVal()); //8
            fitParameters.push_back(sigmaL_addL.getError()); //9 
            fitParameters.push_back(startMpv); //10
        }


        RooPlot* frame = x.frame();

        TCanvas* cFit = new TCanvas(Form("fit_%d",nbin ));
        cFit->Divide(1,2);
        cFit->cd(1); 

        //data.plotOn(frame, RooFit::Name("data"), RooFit::Binning(300));
        RooDataHist dataHist_plot("dataHist_plot", "Binned Data", x, RooFit::Import(*h));
        dataHist_plot.plotOn(frame, RooFit::Name("data"));
 
        chosenModel->plotOn(frame, RooFit::Name("fit"), RooFit::Range("fitRange"), RooFit::NormRange("fitRange"), RooFit::LineColor(kBlue));

        chosenModel->plotOn(frame, RooFit::Components("langau"), RooFit::Name("langau"), RooFit::LineColor(kRed));
        if(modelChoice=="double_langau")
        {
            chosenModel->plotOn(frame, RooFit::Components("langau_add"), RooFit::Name("langau_add"), RooFit::LineColor(kOrange));
        }

        chosenModel->paramOn(frame, RooFit::Layout(0.4, 0.9, 0.9));

        std::vector<double> vec_larghezza = FWHM(tf1,tf1->GetMaximumX(x_low, x_high));


        frame->Draw();

        cFit->cd(2);
        frame->GetXaxis()->SetTitle(""); 

        RooHist *hresid = frame->residHist("data", "fit");
        hresid->GetYaxis()->SetTitle("Data - Fit");
        hresid->GetXaxis()->SetTitle("dE/dx");
        hresid->Draw("AP"); 
        

        fitParameters.push_back(fitRes->minNll()); //7 -- 11
        fitParameters.push_back(fitRes->status()); //8 -- 12
        fitParameters.push_back(fitRes->covQual()); //9 -- 13 

        //fitParameters.push_back(vec_larghezza[2]);  
        //fitParameters.push_back(meanG.getVal()); 
        
        cFit->Write(0, TObject::kOverwrite);
        //cFit->Write();
    

        delete cFit;

        FitOut.fitResults = fitParameters;
        FitOut.fitFunction = tf1 ? (TF1*)tf1->Clone() : nullptr;
        std::map<std::string, TF1*> final_components;
        final_components["main_landau"] = components["main_landau"] ? (TF1*)components["main_landau"]->Clone() : nullptr;
        final_components["background_landau"] = components["background_landau"] ? (TF1*)components["background_landau"]->Clone() : nullptr;
        FitOut.componentsFunctions = final_components;
        FitOut.fraction = frac;

        return FitOut;

}
    
