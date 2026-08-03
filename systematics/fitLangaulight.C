
struct FitResults
{
    double landau_width = -999;
    double landau_mpv = -999;
    double gauss_width = -999;
    double int_landau_width = -999;
    double int_landau_mpv = -999;
    double int_gauss_width = -999;
    double fraction = -999;
    double e_landau_width = -999;
    double e_landau_mpv = -999;
    double e_gauss_width = -999;
    double e_int_landau_width = -999;
    double e_int_landau_mpv = -999;
    double e_int_gauss_width = -999;
    double e_fraction = -999;
    double ndf = -999;
    double chi2 = -999; 
    double status = -999;
    double covqual = -999;
    TF1* fitFunction;
};

FitResults fitDoubleLangau_unbinned(TH1D* h, std::vector<double> v, std::vector<double> weights, double x_low, double x_high, TFile *fitPrintFile , int nbin, TH1D *dedx_range_expected){

        std::vector<double> interacting_start_value = {6., 5.25, 4.5, 4.24, 3.75, 3.75, 3.5, 3.25, 3.25, 3.25, 3.25, 3., 3., 3., 2.75, 2.75, 2.75, 2.75, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5};
    
        fitPrintFile->cd();
        FitResults FitOut;
        RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

        double startMpv=dedx_range_expected -> GetBinContent(dedx_range_expected->FindBin(nbin));
        double startSigmaL = 1.; 
        double sigmaGaus = 1.;
        double int_startSigmaL = 1.5;
        double int_sigmaGaus = 1.5;

        RooRealVar x("x", "x", h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
        RooRealVar w("w", "weight", 0, 1e6);
        RooDataHist dataHist("dataHist", "Dataset from TH1", x, RooFit::Import(*h));
        RooDataSet data("data", "Imported Data", RooArgSet(x, w), RooFit::WeightVar(w));

        for(int i = 0; i < int(v.size()); i++)
        {
            x.setVal(v[i]);
            w.setVal(weights[i]);
            data.add(RooArgSet(x, w), weights[i]);
        }
 
        //LANDAU
        RooRealVar mpv("mpv", "Most Probable Value", startMpv, startMpv-2.5, startMpv+2.5);
        RooRealVar sigmaL("sigmaL", "Landau width", startSigmaL, 0.001, 10.);
        RooLandau landau("landau", "Landau PDF", x, mpv, sigmaL);

        //MEAN 0 GAUSSIAN (SMEARING)
        RooRealVar meanG("meanG", "Mean of Gaussian", 0);
        RooRealVar sigmaG("sigmaG", "Gaussian sigma", sigmaGaus, 0.001, 10.);
        RooGaussian gauss("gauss", "Gaussian PDF", x, meanG, sigmaG);

        //LANDAU X GAUSSIAN CONVOLUTION
        RooFFTConvPdf langau("langau", "Interacting Landau convoluted with Gaussian", x, landau, gauss);

        //LANDAU INTERACTING
        RooRealVar int_mpv("int_mpv", "Interacting Most Probable Value", interacting_start_value[nbin-1], 1., startMpv);
        RooRealVar int_sigmaL("int_sigmaL", "Interacting Landau width", int_startSigmaL, 0.001, 10.);
        RooLandau int_landau("int_landau", "Interacting Landau PDF", x, int_mpv, int_sigmaL);

        //MEAN 0 GAUSSIAN (SMEARING) INTERACTING
        RooRealVar int_meanG("int_meanG", "Interacting Mean of Gaussian", 0);
        RooRealVar int_sigmaG("int_sigmaG", "Interacting Gaussian sigma", int_sigmaGaus, 0.001, 10.);
        RooGaussian int_gauss("int_gauss", "Interacting Gaussian PDF", x, int_meanG, int_sigmaG);

        //LANDAU X GAUSSIAN CONVOLUTION INTERACTING
        RooFFTConvPdf int_langau("int_langau", "Interacting Landau convoluted with Gaussian", x, int_landau, int_gauss);
        
        RooRealVar frac("frac", "fraction of signal", 0.5, 0, 1);

        RooAddPdf model("model", "double langau", RooArgList(langau, int_langau), RooArgList(frac));
    
        x.setRange("fitRange", x_low, x_high); 

        // Select model
        RooAbsPdf* chosenModel = &model;

        x.setRange("fitRange", x_low, x_high);

        RooFitResult* fitRes = chosenModel->fitTo(data, RooFit::Save(), RooFit::Range("fitRange"), RooFit::SumW2Error(true), RooFit::PrintLevel(0), RooFit::PrintEvalErrors(-1), RooFit::Warnings(false), RooFit::Verbose(false));

        FitOut.landau_width = sigmaL.getVal();
        FitOut.landau_mpv = mpv.getVal();
        FitOut.gauss_width = sigmaG.getVal();
        FitOut.int_landau_width = int_sigmaL.getVal();
        FitOut.int_landau_mpv = int_mpv.getVal();
        FitOut.int_gauss_width = int_sigmaG.getVal();
        FitOut.fraction = frac.getVal();
        FitOut.e_landau_width = sigmaL.getError();
        FitOut.e_landau_mpv = mpv.getError();
        FitOut.e_gauss_width = sigmaG.getError();
        FitOut.e_int_landau_width = int_sigmaL.getError();
        FitOut.e_int_landau_mpv = int_mpv.getError();
        FitOut.e_int_gauss_width = int_sigmaG.getError();
        FitOut.e_fraction = frac.getError();

        TF1* tf1 = chosenModel->asTF(RooArgList(x), RooArgList(mpv, sigmaL, sigmaG, int_mpv, int_sigmaL, int_sigmaG, frac));
        if (tf1) tf1->SetNpx(10000);

        
        RooPlot* frame = x.frame();

        TCanvas* cFit = new TCanvas(Form("fit_%d",nbin ));
        cFit->Divide(1,2);
        cFit->cd(1); 

        data.plotOn(frame, RooFit::Name("data"), RooFit::Binning(h->GetNbinsX()));
 
        chosenModel->plotOn(frame, RooFit::Name("fit"), RooFit::Range("fitRange"), RooFit::NormRange("fitRange"), RooFit::LineColor(kBlue));

        // Componente "non interagente" (langau) 
        chosenModel->plotOn(frame, RooFit::Name("langau_comp"), RooFit::Components("langau"), RooFit::Range("fitRange"), RooFit::NormRange("fitRange"), RooFit::LineColor(kGreen+2), RooFit::LineStyle(kDashed));

        // Componente "interagente" (int_langau)
        chosenModel->plotOn(frame, RooFit::Name("int_langau_comp"), RooFit::Components("int_langau"), RooFit::Range("fitRange"), RooFit::NormRange("fitRange"), RooFit::LineColor(kRed), RooFit::LineStyle(kDashed));


        chosenModel->paramOn(frame, RooFit::Layout(0.4, 0.9, 0.9));

        frame->Draw();

        cFit->cd(2);
        frame->GetXaxis()->SetTitle(""); 

        RooHist *hresid = frame->residHist("data", "fit");
        hresid->GetYaxis()->SetTitle("Data - Fit");
        hresid->GetXaxis()->SetTitle("dE/dx");
        hresid->Draw("AP"); 
        
        RooAbsReal* chi2 = chosenModel->createChi2(dataHist, RooFit::DataError(RooAbsData::Poisson));
        int ndf = dataHist.numEntries() - model.getParameters(dataHist)->getSize();
         
        FitOut.ndf = ndf;
        FitOut.chi2 = chi2->getVal();
        FitOut.status = fitRes->status();
        FitOut.covqual = fitRes->covQual(); 
        
        cFit->Write(0, TObject::kOverwrite);

        delete cFit;

        FitOut.fitFunction = tf1 ? (TF1*)tf1->Clone() : nullptr;

        return FitOut;

}
    

FitResults fitSingleLangau_unbinned(TH1D* h, std::vector<double> v, std::vector<double> weights, double x_low, double x_high, TFile *fitPrintFile , int nbin, TH1D *dedx_range_expected){

    
        fitPrintFile->cd();
        FitResults FitOut;
        RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

        double startMpv=dedx_range_expected -> GetBinContent(dedx_range_expected->FindBin(nbin));
        double startSigmaL = 0.05; 
        double sigmaGaus = 0.5;

        RooRealVar x("x", "x", h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
        RooRealVar w("w", "weight", 0, 1e6);
        RooDataHist dataHist("dataHist", "Dataset from TH1", x, RooFit::Import(*h));
        RooDataSet data("data", "Imported Data", RooArgSet(x, w), RooFit::WeightVar(w));

        for(int i = 0; i < int(v.size()); i++)
        {
            x.setVal(v[i]);
            w.setVal(weights[i]);
            data.add(RooArgSet(x, w), weights[i]);
        }
 
        //LANDAU
        RooRealVar mpv("mpv", "Most Probable Value", startMpv, startMpv-2.5, startMpv+2.5);
        RooRealVar sigmaL("sigmaL", "Landau width", startSigmaL, 0.001, 5.);
        RooLandau landau("landau", "Landau PDF", x, mpv, sigmaL);

        //MEAN 0 GAUSSIAN (SMEARING)
        RooRealVar meanG("meanG", "Mean of Gaussian", 0);
        RooRealVar sigmaG("sigmaG", "Gaussian sigma", sigmaGaus, 0.001, 5.);
        RooGaussian gauss("gauss", "Gaussian PDF", x, meanG, sigmaG);

        //LANDAU X GAUSSIAN CONVOLUTION
        RooFFTConvPdf model("langau", "Interacting Landau convoluted with Gaussian", x, landau, gauss);
    
        x.setRange("fitRange", x_low, x_high); 

        // Select model
        RooAbsPdf* chosenModel = &model;

        x.setRange("fitRange", x_low, x_high);

        RooFitResult* fitRes = chosenModel->fitTo(data, RooFit::Save(), RooFit::Range("fitRange"), RooFit::SumW2Error(true), RooFit::PrintLevel(0), RooFit::PrintEvalErrors(-1), RooFit::Warnings(false), RooFit::Verbose(false));

        FitOut.landau_width = sigmaL.getVal();
        FitOut.landau_mpv = mpv.getVal();
        FitOut.gauss_width = sigmaG.getVal();
        FitOut.e_landau_width = sigmaL.getError();
        FitOut.e_landau_mpv = mpv.getError();
        FitOut.e_gauss_width = sigmaG.getError();

        TF1* tf1 = chosenModel->asTF(RooArgList(x), RooArgList(mpv, sigmaL, sigmaG));
        if (tf1) tf1->SetNpx(10000);

        
        RooPlot* frame = x.frame();

        TCanvas* cFit = new TCanvas(Form("fit_%d",nbin ));
        cFit->Divide(1,2);
        cFit->cd(1); 

        data.plotOn(frame, RooFit::Name("data"), RooFit::Binning(h->GetNbinsX()));
 
        chosenModel->plotOn(frame, RooFit::Name("fit"), RooFit::Range("fitRange"), RooFit::NormRange("fitRange"), RooFit::LineColor(kBlue));

        chosenModel->paramOn(frame, RooFit::Layout(0.4, 0.9, 0.9));

        frame->Draw();

        cFit->cd(2);
        frame->GetXaxis()->SetTitle(""); 

        RooHist *hresid = frame->residHist("data", "fit");
        hresid->GetYaxis()->SetTitle("Data - Fit");
        hresid->GetXaxis()->SetTitle("dE/dx");
        hresid->Draw("AP"); 
        
        RooAbsReal* chi2 = chosenModel->createChi2(dataHist, RooFit::DataError(RooAbsData::Poisson));
        int ndf = dataHist.numEntries() - model.getParameters(dataHist)->getSize();
         
        FitOut.ndf = ndf;
        FitOut.chi2 = chi2->getVal();
        FitOut.status = fitRes->status();
        FitOut.covqual = fitRes->covQual(); 
        
        cFit->Write(0, TObject::kOverwrite);

        delete cFit;

        FitOut.fitFunction = tf1 ? (TF1*)tf1->Clone() : nullptr;

        return FitOut;

}