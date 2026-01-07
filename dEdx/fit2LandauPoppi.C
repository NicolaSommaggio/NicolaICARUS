#include <iostream>
#include <vector>
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"

struct fitPar{
    double mpv;
    double mpvMax;
    double mpvMin;
    double witdth;
    double witdthMax;
    double witdthMin;
    double rangeMin;
    double rangeMax;
};

struct HistStats {
    double mean;
    double mpv;
    double fwhm;
};

struct fitValues {
    double mpv;
    double mpv_err;
    double sigmaL;
    double sigmaL_err;
    double sigmaG;
    double sigmaG_err;
};

std::map<std::pair<int,int>,fitPar> dataFitPar;
std::map<std::pair<int,int>,fitPar> cvFitPar;
std::map<std::pair<int,int>,fitPar> varFitPar;

int nbins=2000;

TCanvas* makedQdXAnalysis(TH1D* histo_dqdx, int plane, int TPC, fitValues& theseFittedValues){
    
    std::cout<<"Histogram Entries "<<histo_dqdx->GetEntries()<<std::endl;
    RooRealVar dQdX("dQdX", "dQdX", 0., 6000.);
    dQdX.setBins(nbins);

    double meanL=700., minL=500., maxL=850.;
    double lowRange = 400., upRange = 3000.;
    if(plane==2){
        minL=650.;
        maxL=800.;
        lowRange=500.;
        upRange=1500.;
    }
    // Landau + Gauss
    RooRealVar mpv("mpv", "MPV of Landau", meanL, minL, maxL);
    RooRealVar sigmaL("sigmaL", "Width of Landau", 40., 0., 100.);
    RooLandau landau("landau", "Landau distribution", dQdX, mpv, sigmaL);
    RooRealVar mean("mean", "Mean of Gaussian", 0.);
    RooRealVar sigmaG("sigmaG", "Width of Gaussian", 0., 400.);
    RooGaussian gauss("gauss", "Gaussian distribution", dQdX, mean, sigmaG);
    RooFFTConvPdf landauGauss("landauGauss", "Landau (x) Gaussian", dQdX, landau, gauss);
    RooDataHist dataHist("dataHist", "dataHist", dQdX, RooFit::Import(*histo_dqdx));
    dQdX.setRange("fitRange", lowRange, upRange); // Fit solo tra 10 e 50

    landauGauss.fitTo(dataHist, RooFit::Range("fitRange"), RooFit::PrintLevel(-1), RooFit::Verbose(kFALSE), RooFit::Warnings(kFALSE), RooFit::PrintEvalErrors(0));
    RooPlot* frame = dQdX.frame();
    dataHist.plotOn(frame);
    landauGauss.plotOn(frame, RooFit::Range("fitRange"), RooFit::NormRange("fitRange"),RooFit::LineStyle(kDashed), RooFit::LineColor(kOrange));
    landauGauss.paramOn(frame, RooFit::Layout(0.6, 0.9, 0.9));
    TCanvas* c = new TCanvas("c", "Fit Landau + Gaussian", 800, 600);
    frame->Draw();
    theseFittedValues = {mpv.getVal(), mpv.getError(), sigmaL.getVal(), sigmaL.getError(), sigmaG.getVal(), sigmaG.getError()};
    return c;
}

TCanvas* makedEdXAnalysis_2L(TH1D* histo_dqdx, int plane, int TPC, fitValues& theseFittedValues){
    TCanvas* c = new TCanvas("c", "Fit Landau + Gaussian", 800, 600);
    std::cout<<"Histogram Entries "<<histo_dqdx->GetEntries()<<std::endl;
    RooRealVar dQdX("dEdX", "dEdX", 0., 15.);
    dQdX.setBins(nbins);

    double meanL=1.8, minL=1.5, maxL=2.2;
    double lowRange = 1., upRange = 6.;
    if(plane==1){
        minL=1.2;
        lowRange = 0.8;
//        upRange = 6.;
    }
    if(plane==2){
        minL=1.6;
        maxL=2.0;
        lowRange = 1.2;
        upRange = 10.;
    }
    if(plane==0){

        lowRange = 0.8;
        upRange = 8.;
    }
    /*
    if(plane==2){
        minL=650.;
        maxL=800.;
        lowRange=500.;
        upRange=1500.;
    } */

    // Landau + Gauss
    RooRealVar mpv("mpv", "MPV of Landau", meanL, minL, maxL);
    RooRealVar sigmaL("sigmaL", "Width of Landau", 0.1, 0., 1.);
    RooLandau landau("landau", "Landau distribution", dQdX, mpv, sigmaL);
    RooRealVar mean("mean", "Mean of Gaussian", 0.);
    RooRealVar sigmaG("sigmaG", "Width of Gaussian", 0., 1.);
    RooGaussian gauss("gauss", "Gaussian distribution", dQdX, mean, sigmaG);
    RooFFTConvPdf landauGauss("landauGauss", "Landau (x) Gaussian", dQdX, landau, gauss);

    //RooRealVar mpv_bkg("mpv_bkg", "mpv_bkg", 2*meanL, 3., 4.);
    RooFormulaVar mpv_bkg("mpv_bkg", "mpv_bkg", "2*@0", RooArgList(mpv));
    RooRealVar sigmaL_bkg("sigmaL_bkg", "Width of Landau", 0.3, 0., 0.7);
    RooLandau landau_bkg("landau_bkg", "Landau distribution", dQdX, mpv_bkg, sigmaL_bkg);
    RooRealVar mean_bkg("mean_bkg", "Mean of Gaussian bkg", 0.);
    RooRealVar sigmaG_bkg("sigmaG_bkg", "Width of Gaussian", 0.8, 0.3, 2.);
    RooGaussian gauss_bkg("gauss_bkg", "Gaussian distribution", dQdX, mean_bkg, sigmaG_bkg);
    RooFFTConvPdf landauGauss_bkg("landauGauss_bkg", "Landau (x) Gaussian", dQdX, landau_bkg, gauss_bkg);

    RooRealVar mean_l("mean_l", "Mean of Gaussian bkg", meanL, 1.2, maxL);
    RooRealVar sigmaG_l("sigmaG_l", "Width of Gaussian", 0., 2.);
    RooLandau landau_l("landau_l", "Landau distribution", dQdX, mean_l, sigmaG_l);
    
    RooRealVar mean_ll("mean_ll", "Mean of Gaussian bkg", meanL, 1.2, maxL);
    RooRealVar sigmaG_ll("sigmaG_ll", "Width of Gaussian", 0., 2.);
    RooGaussian gauss_ll("gauss_ll", "Gaussian distribution", dQdX, mean_ll, sigmaG_ll);

    RooRealVar frac("frac", "Fraction of Landau+Gauss in total model", 0.85, 0.8, 1.);

    RooRealVar frac2("frac2", "Fraction of Gauss_l small in total model", 0.01, 0., 0.2);

    RooRealVar frac3("frac3", "Fraction of Gauss_ll small in total model", 0.01, 0., 0.2);

    RooAddPdf model("model", "Landau+Gauss + Gauss + schifo", RooArgList(landauGauss, landauGauss_bkg), RooArgList(frac));

    RooDataHist dataHist("dataHist", "dataHist", dQdX, RooFit::Import(*histo_dqdx));
    dQdX.setRange("fitRange", lowRange, upRange);
    model.fitTo(dataHist, RooFit::Range("fitRange"), RooFit::PrintLevel(-1), RooFit::Verbose(kFALSE), RooFit::Warnings(kFALSE), RooFit::PrintEvalErrors(0));
    RooPlot* frame = dQdX.frame();
    dataHist.plotOn(frame);
    model.plotOn(frame, RooFit::Range("fitRange"), RooFit::NormRange("fitRange"));
    model.plotOn(frame, RooFit::Components("landauGauss"), RooFit::LineStyle(kDashed), RooFit::LineColor(kRed));
    model.plotOn(frame, RooFit::Components("landauGauss_bkg"), RooFit::LineStyle(kDashed), RooFit::LineColor(kGreen));
    model.plotOn(frame, RooFit::Components("landau_l"), RooFit::LineStyle(kDashed), RooFit::LineColor(kOrange));
    model.plotOn(frame, RooFit::Components("gauss_ll"), RooFit::LineStyle(kDashed), RooFit::LineColor(kOrange));
    
    model.paramOn(frame, RooFit::Layout(0.4, 0.9, 0.9));
    frame->Draw();
    theseFittedValues = {mpv.getVal(), mpv.getError(), sigmaL.getVal(), sigmaL.getError(), sigmaG.getVal(), sigmaG.getError()};
    return c;
}

TCanvas* makedQdXAnalysis_2L(TH1D* histo_dqdx, int plane, int TPC, fitValues& theseFittedValues){
    
    std::cout<<"Histogram Entries "<<histo_dqdx->GetEntries()<<std::endl;
    RooRealVar dQdX("dQdX", "dQdX", 0., 6000.);
    dQdX.setBins(nbins);

    double meanL=700., minL=500., maxL=850.;
    double lowRange = 300., upRange = 2500.;
    if(plane==2){
        minL=650.;
        maxL=800.;
        lowRange=300.;
        upRange=2500.;
    }
    // Landau + Gauss
    RooRealVar mpv("mpv", "MPV of Landau", meanL, minL, maxL);
    RooRealVar sigmaL("sigmaL", "Width of Landau", 40., 0., 100.);
    RooLandau landau("landau", "Landau distribution", dQdX, mpv, sigmaL);
    RooRealVar mean("mean", "Mean of Gaussian", 0.);
    RooRealVar sigmaG("sigmaG", "Width of Gaussian", 0., 200.);
    RooGaussian gauss("gauss", "Gaussian distribution", dQdX, mean, sigmaG);
    RooFFTConvPdf landauGauss("landauGauss", "Landau (x) Gaussian", dQdX, landau, gauss);

    //RooRealVar mpv_bkg("mpv_bkg", "mpv_bkg", 2*meanL, 1000., 1600);
    RooFormulaVar mpv_bkg("mpv_bkg", "mpv_bkg", "2*@0", RooArgList(mpv));
    RooRealVar sigmaL_bkg("sigmaL_bkg", "Width of Landau",40., 0., 200.);
    RooLandau landau_bkg("landau_bkg", "Landau distribution", dQdX, mpv_bkg, sigmaL);
    RooRealVar mean_bkg("mean_bkg", "Mean of Gaussian bkg", 0.);
    RooRealVar sigmaG_bkg("sigmaG_bkg", "Width of Gaussian", 0., 300.);
    RooGaussian gauss_bkg("gauss_bkg", "Gaussian distribution", dQdX, mean_bkg, sigmaG_bkg);
    RooFFTConvPdf landauGauss_bkg("landauGauss_bkg", "Landau (x) Gaussian", dQdX, landau_bkg, gauss_bkg);
    RooRealVar frac("frac", "Fraction of Landau+Gauss in total model", 0.85, 0.8, 1.);
    RooAddPdf model("model", "Landau+Gauss + Gauss + Gauss", RooArgList(landauGauss, landauGauss_bkg), RooArgList(frac));

    RooDataHist dataHist("dataHist", "dataHist", dQdX, RooFit::Import(*histo_dqdx));
    dQdX.setRange("fitRange", lowRange, upRange);
    model.fitTo(dataHist, RooFit::Range("fitRange"), RooFit::PrintLevel(-1), RooFit::Verbose(kFALSE), RooFit::Warnings(kFALSE), RooFit::PrintEvalErrors(0));
    RooPlot* frame = dQdX.frame();
    dataHist.plotOn(frame);
    model.plotOn(frame, RooFit::Range("fitRange"), RooFit::NormRange("fitRange"));
    model.plotOn(frame, RooFit::Components("landauGauss"), RooFit::LineStyle(kDashed), RooFit::LineColor(kRed));
    model.plotOn(frame, RooFit::Components("landauGauss_bkg"), RooFit::LineStyle(kDashed), RooFit::LineColor(kGreen));
    model.paramOn(frame, RooFit::Layout(0.4, 0.9, 0.9));
    TCanvas* c = new TCanvas("c", "Fit Landau + Gaussian", 800, 600);
    frame->Draw();
    theseFittedValues = {mpv.getVal(), mpv.getError(), sigmaL.getVal(), sigmaL.getError(), sigmaG.getVal(), sigmaG.getError()};
    return c;
}

TCanvas* makedEdXAnalysis(TH1D* histo_dqdx, int plane, int TPC, fitValues& theseFittedValues){
    
    std::cout<<"Histogram Entries "<<histo_dqdx->GetEntries()<<std::endl;
    RooRealVar dQdX("dEdX", "dEdX", 0., 15.);
    dQdX.setBins(nbins);

    double meanL=1.8, minL=1.5, maxL=2.2;
    double lowRange = 1., upRange = 4.;
    /*
    if(plane==2){
        minL=650.;
        maxL=800.;
        lowRange=500.;
        upRange=1500.;
    } */
    // Landau + Gauss
    RooRealVar mpv("mpv", "MPV of Landau", meanL, minL, maxL);
    RooRealVar sigmaL("sigmaL", "Width of Landau", 0.1, 0., 100.);
    RooLandau landau("landau", "Landau distribution", dQdX, mpv, sigmaL);
    RooRealVar mean("mean", "Mean of Gaussian", 0.);
    RooRealVar sigmaG("sigmaG", "Width of Gaussian", 0., 10.);
    RooGaussian gauss("gauss", "Gaussian distribution", dQdX, mean, sigmaG);
    RooFFTConvPdf landauGauss("landauGauss", "Landau (x) Gaussian", dQdX, landau, gauss);
    RooDataHist dataHist("dataHist", "dataHist", dQdX, RooFit::Import(*histo_dqdx));
    dQdX.setRange("fitRange", lowRange, upRange); // Fit solo tra 10 e 50

    landauGauss.fitTo(dataHist, RooFit::Range("fitRange"), RooFit::PrintLevel(-1), RooFit::Verbose(kFALSE), RooFit::Warnings(kFALSE), RooFit::PrintEvalErrors(0));
    RooPlot* frame = dQdX.frame();
    dataHist.plotOn(frame);
    landauGauss.plotOn(frame, RooFit::Range("fitRange"), RooFit::NormRange("fitRange"),RooFit::LineStyle(kDashed), RooFit::LineColor(kOrange));
    landauGauss.paramOn(frame, RooFit::Layout(0.6, 0.9, 0.9));
    TCanvas* c = new TCanvas("c", "Fit Landau + Gaussian", 800, 600);
    frame->Draw();
    theseFittedValues = {mpv.getVal(), mpv.getError(), sigmaL.getVal(), sigmaL.getError(), sigmaG.getVal(), sigmaG.getError()};
    return c;
}

fitValues returnValues(std::string etype, int variation, int tpc, int plane, int type, int selection){

    fitValues theseFittedValues = {0,0,0,0,0,0};
    // type: 0 data, 1 CV, 2 Var
    if (etype != "dQdX" && etype != "dEdX") {
        std::cerr << "Errore: etype 'dQdX' o 'dEdX'!" << std::endl;
        return theseFittedValues;
    }
    if (plane < 0 || plane > 2) {
        std::cerr << "Errore: plane 0 (I), 1 (II) o 2 (C)!" << std::endl;
        return theseFittedValues;
    }
    if (tpc < 0 || tpc > 3) {
        std::cerr << "Errore: TPC 0 (EE), 1 (EW), 2 (WE) o 3 (WW)!" << std::endl;
        return theseFittedValues;
    }
    std::vector<std::string> vartype = {"CathodeBending_histograms.root", "GainHi_histograms.root", "GainVar_histograms.root", "HighLT_histograms.root", "HighTPCNoise_histograms.root", "HighTPCCohNoise_histograms.root", "IndGap1WireFil_histograms.root", "LowTPCNoise_histograms.root", "LowTPCCohNoise_histograms.root", "NullVar_histograms.root","RecoMod_histograms.root","ScintQE_histograms.root","LowLT_histograms.root","TPCYZ_histograms.root"};
    std::vector<std::string> varlabels = {"_CathodeBending.root", "_GainHi.root", "_GainVar.root", "_HighLT.root", "_HighTPCNoise.root", "_HighTPCCohNoise.root", "_IndGap1WireFil.root", "_LowTPCNoise.root", "_LowTPCCohNoise.root", "_NullVar.root","_RecoMod.root","_ScintQE.root","_LowLT.root","_TPCYZ.root"};

    TFile* fileData = TFile::Open("Data_histograms.root", "READ");
    //TFile* fileMC = TFile::Open("CV_histograms.root", "READ");
    //TFile* fileVar = TFile::Open(vartype[variation].c_str(), "READ");

    std::vector<std::string> planeNames = {"I", "II", "C"};
    std::vector<std::string> TPCNames = {"EE", "EW", "WE", "WW"};
    std::vector<std::string> SelectionNames = {"0","1","2"};
    std::vector<std::string> dataType = {"data","CV","variation"};
    
    std::string histName = etype + "_" + TPCNames[tpc] + "_" + planeNames[plane] + "_" + SelectionNames[selection];
    
    //if (!fileData || !fileMC || !fileVar) {
    //    std::cerr << "Error file" << std::endl;
    //    return;
    //}

    TH1D* histD = (TH1D*)fileData->Get(histName.c_str());
    //TH1D* histMC = (TH1D*)fileMC->Get(histName.c_str());
    //TH1D* histVar = (TH1D*)fileVar->Get(histName.c_str());

    //if (!histD || !histMC || !histVar) {
    //    std::cerr << "Error histo" << std::endl;
    //    return;
    //}
    TCanvas* fittedHisto;
    if(plane==2){
        fittedHisto = (etype == "dQdX") ? makedQdXAnalysis_2L(histD, plane, tpc, theseFittedValues) : makedEdXAnalysis_2L(histD, plane, tpc, theseFittedValues);
    /*
    if(type==0){ // data
        fittedHisto = (etype == "dQdX") ? makedQdXAnalysis_2L(histD, plane, tpc) : makedEdXAnalysis_2L(histD, plane, tpc);
    } else if(type==1){ // CV
        fittedHisto = (etype == "dQdX") ? makedQdXAnalysis_2L(histMC, plane, tpc) : makedEdXAnalysis_2L(histMC, plane, tpc);
    } else if(type==2){ // Var
        fittedHisto = (etype == "dQdX") ? makedQdXAnalysis_2L(histVar, plane, tpc) : makedEdXAnalysis_2L(histVar, plane, tpc);
    }*/
    } else if(plane==1){
        fittedHisto = (etype == "dQdX") ? makedQdXAnalysis_2L(histD, plane, tpc, theseFittedValues) : makedEdXAnalysis_2L(histD, plane, tpc, theseFittedValues);
        /*
        if(type==0){ // data
            fittedHisto = (etype == "dQdX") ? makedQdXAnalysis(histD, plane, tpc) : makedEdXAnalysis(histD, plane, tpc);
        } else if(type==1){ // CV
            fittedHisto = (etype == "dQdX") ? makedQdXAnalysis(histMC, plane, tpc) : makedEdXAnalysis(histMC, plane, tpc);
        } else if(type==2){ // Var
            fittedHisto = (etype == "dQdX") ? makedQdXAnalysis(histVar, plane, tpc) : makedEdXAnalysis(histVar, plane, tpc);
        }*/ 
    } else if(plane==0){
        fittedHisto = (etype == "dQdX") ? makedQdXAnalysis_2L(histD, plane, tpc, theseFittedValues) : makedEdXAnalysis_2L(histD, plane, tpc, theseFittedValues);
        /*
        if(type==0){ // data
            fittedHisto = (etype == "dQdX") ? makedQdXAnalysis(histD, plane, tpc) : makedEdXAnalysis(histD, plane, tpc);
        } else if(type==1){ // CV
            fittedHisto = (etype == "dQdX") ? makedQdXAnalysis(histMC, plane, tpc) : makedEdXAnalysis(histMC, plane, tpc);
        } else if(type==2){ // Var
            fittedHisto = (etype == "dQdX") ? makedQdXAnalysis(histVar, plane, tpc) : makedEdXAnalysis(histVar, plane, tpc);
        }*/ 
    }
    std::string imgNameVar = "Fit/"+histName + "_"+ dataType[type];
    if(type==2) imgNameVar = imgNameVar + "_" + varlabels[variation];
    imgNameVar = imgNameVar + ".root";
    fittedHisto->SaveAs(imgNameVar.c_str());
    return theseFittedValues;
}

void getValues(){
    std::string etype = "dEdX";
    int type=0; // data
    int variation=0; // variation
    for(int sel=0; sel<3; sel++){
        for(int plane=0; plane<3; plane++){
            for(int tpc=0; tpc<4; tpc++){
                fitValues thisfit = returnValues(etype, variation, tpc, plane, type, sel);
                std::cout<<etype<<" , TPC: "<<tpc<<" , plane "<<plane<<" , selection "<<sel<<":\n";
                std::cout<<"MPV "<<thisfit.mpv<<" p/m "<<thisfit.mpv_err<<" SigmaL "<<thisfit.sigmaL<<" p/m "<<thisfit.sigmaL_err<<" SigmaG "<<thisfit.sigmaG<<" p/m "<<thisfit.sigmaG_err<<std::endl;
                
            }
        }
    }
    etype = "dQdX";
    for(int sel=0; sel<3; sel++){
        for(int plane=0; plane<3; plane++){
            for(int tpc=0; tpc<4; tpc++){
                fitValues thisfit = returnValues(etype, variation, tpc, plane, type, sel);
                std::cout<<etype<<" , TPC: "<<tpc<<" , plane "<<plane<<" , selection "<<sel<<":\n";
                std::cout<<"MPV "<<thisfit.mpv<<" p/m "<<thisfit.mpv_err<<" SigmaL "<<thisfit.sigmaL<<" p/m "<<thisfit.sigmaL_err<<" SigmaG "<<thisfit.sigmaG<<" p/m "<<thisfit.sigmaG_err<<std::endl;
                
            }
        }
    }
}


HistStats AnalyzeHistogram(TH1D* h) {
    if (!h) {
        std::cerr << "Errore: l'istogramma è nullo!" << std::endl;
        return {0, 0, 0};
    }
    double mean = h->GetMean();
    int binMax = h->GetMaximumBin();
    double mpv = h->GetXaxis()->GetBinCenter(binMax);
    double halfMax = h->GetMaximum() / 2.0;
    int binLeft = binMax, binRight = binMax;
    while (binLeft > 1 && h->GetBinContent(binLeft) > halfMax) {
        binLeft--;
    }
    double xLeft = h->GetXaxis()->GetBinCenter(binLeft);
    while (binRight < h->GetNbinsX() && h->GetBinContent(binRight) > halfMax) {
        binRight++;
    }
    double xRight = h->GetXaxis()->GetBinCenter(binRight);
    double fwhm = xRight - xLeft;
    return {mean, mpv, fwhm};
}

std::string formatNumber(double value) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(6) << value;
    std::string str = oss.str();
    std::replace(str.begin(), str.end(), '.', ','); // Sostituiamo il punto con la virgola
    return str;
}


void doHistoStats() {
    std::vector<std::string> vartype = {
        "CathodeBending_histograms.root", "GainHi_histograms.root", "GainVar_histograms.root",
        "HighLT_histograms.root", "HighTPCNoise_histograms.root", "HighTPCCohNoise_histograms.root",
        "IndGap1WireFil_histograms.root", "LowTPCNoise_histograms.root", "LowTPCCohNoise_histograms.root",
        "NullVar_histograms.root", "RecoMod_histograms.root", "ScintQE_histograms.root",
        "LowLT_histograms.root", "TPCYZ_histograms.root"
    };

    TFile* fileData = TFile::Open("Data_histograms.root", "READ");
    TFile* fileMC = TFile::Open("CV_histograms.root", "READ");

    std::vector<std::string> planeNames = {"I", "II", "C"};
    std::vector<std::string> TPCNames = {"EE", "EW", "WE", "WW"};

    for (std::string etype : {"dEdX", "dQdX"}) { // Ciclo su dEdX e dQdX
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 3; j++) {
                std::string histName = etype + "_" + TPCNames[i] + "_" + planeNames[j];
                std::string filename = etype + "_" + std::to_string(i) + "_" + std::to_string(j) + ".txt";

                // Apriamo il file di output
                std::ofstream outFile(filename);
                if (!outFile.is_open()) {
                    std::cerr << "Errore nell'apertura del file: " << filename << std::endl;
                    continue;
                }

                // Scriviamo l'intestazione
                outFile << "Source\tMean\tMPV\tFWHM\n";

                // Analizziamo gli istogrammi Data e MC
                TH1D* histD = (TH1D*)fileData->Get(histName.c_str());
                TH1D* histMC = (TH1D*)fileMC->Get(histName.c_str());

                HistStats dataStat = AnalyzeHistogram(histD);
                HistStats mcStat = AnalyzeHistogram(histMC);

                // Scriviamo i dati di Data e MC nel file, con numeri formattati con la virgola
                outFile << "Data\t" << formatNumber(dataStat.mean) << "\t"
                        << formatNumber(dataStat.mpv) << "\t"
                        << formatNumber(dataStat.fwhm) << "\n";

                outFile << "MC\t" << formatNumber(mcStat.mean) << "\t"
                        << formatNumber(mcStat.mpv) << "\t"
                        << formatNumber(mcStat.fwhm) << "\n";

                // Analizziamo le variazioni dei dati
                for (int k = 0; k < 14; k++) {
                    TFile* fileVar = TFile::Open(vartype[k].c_str(), "READ");
                    TH1D* histVar = (TH1D*)fileVar->Get(histName.c_str());

                    if (!histVar) {
                        std::cerr << "Errore: istogramma " << histName << " non trovato in " << vartype[k] << std::endl;
                        continue;
                    }

                    HistStats varStat = AnalyzeHistogram(histVar);
                    outFile << vartype[k] << "\t" 
                            << formatNumber(varStat.mean) << "\t"
                            << formatNumber(varStat.mpv) << "\t"
                            << formatNumber(varStat.fwhm) << "\n";

                    fileVar->Close(); // Chiudiamo il file ROOT
                }

                outFile.close(); // Chiudiamo il file di output
                std::cout << "File salvato: " << filename << std::endl;
            }
        }
    }

    fileData->Close();
    fileMC->Close();
}