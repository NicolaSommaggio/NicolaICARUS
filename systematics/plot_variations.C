#include "slice_struct.h"
#include "fitLangaus_unbinned.C"

void plotEk()
{
    TFile * infile_MC = TFile::Open("varMC_STANDARD_KE.root","READ"); 

    TTree * intree = (TTree*)infile_MC->Get("tree");

    _slice *slc = nullptr;

    intree -> SetBranchAddress("slice",&slc);

    TH1D *ek_calo_CAF = new TH1D("ek_calo_CAF","",100,0.,1000.); 
    TH1D *ek_calo = new TH1D("ek_calo","",100,0.,1000.);
    TH1D *diff = new TH1D("diff","",100,-0.02,0.02);

    for(int i = 0; i < intree -> GetEntries(); i ++)
    {
        intree -> GetEntry(i);

        for (const auto &pro : slc->_protons)
        {
            if(pro._length > 25.)continue;

            ek_calo_CAF -> Fill(pro._KE[0]);

            //cout << "CAF: " << pro._KE[0] << endl;

            double ek = 0;
            for(int hit = 0; hit < pro._dedx.size(); hit ++)
            {
                ek = ek + pro._dedx[hit] * pro._pitch[hit];
            }

            //cout << "CALO: " << ek << endl;

            ek_calo -> Fill(ek);

            diff -> Fill(ek-pro._KE[0]);
        }
    }
    TFile *outfile = new TFile("outfile_ke.root","RECREATE");
    outfile -> cd();
    ek_calo_CAF -> Write();
    ek_calo -> Write();
    diff -> Write();

}

double p_gauss_pdf(double *x, double *par)
{
    const Double_t sqrt2 = TMath::Sqrt(2.0);

    // *** 1 GAUSSIAN ***
    double gnorm = 0.5 * (
            TMath::Erf((X_HIGH - par[0])/(sqrt2*par[1]))
          - TMath::Erf((X_LOW - par[0])/(sqrt2*par[1]))
        );
    return TMath::Gaus(x[0], par[0], par[1], 1) / gnorm;
}

void PLOT_VARIATIONS()
{
    gROOT->SetBatch(kTRUE);

    std::vector<std::string> samples = {"MC","OFFBEAM"};
    std::vector<double> POTS = {2.43007e+20, 2.11e+20};
    TFile * outfile = new TFile("OUT_VARIATIONS_STANDARD_RUN2_DATA.test.root","RECREATE"); //-->HERE
            
    //TFile * infile_MC = TFile::Open("varMC_STANDARD_RUN2.root","READ");  
    //TFile * infile_OFFBEAM = TFile::Open("varMC_STANDARD_RUN2_OFFBEAM.root","READ");        
    //TFile * infile_DATA = TFile::Open("varMC_STANDARD_RUN2_DATA.root","READ");

    TFile * infile_MC = TFile::Open("ROOT_TREES_DEDX/RUN2/varMC_STANDARD_RUN2_CUT.root","READ");  
    TFile * infile_OFFBEAM = TFile::Open("ROOT_TREES_DEDX/RUN2/varMC_STANDARD_RUN2_OFFBEAM_CUT.root","READ");        
    TFile * infile_DATA = TFile::Open("ROOT_TREES_DEDX/RUN2/varMC_STANDARD_RUN2_DATA_CUT.root","READ");

    std::vector<TFile*> infile = {infile_DATA, infile_MC, infile_OFFBEAM};

    std::vector<int> xhigh_pro = {30, 30, 30, 25, 20, 20, 20, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 10, 10, 10, 10, 10, 10};
    std::vector<int> nbins_pro = {90, 120, 120, 100, 80, 80, 80, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 30, 30, 30, 30, 30, 30}; //-->HERE
    //std::vector<int> nbins_pro = {150, 150, 150, 125, 100, 100, 100, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 50, 50, 50, 50, 50, 50}; //-->HERE

    std::vector<double> interacting_start_value = {6., 5.25, 4.5, 4.24, 3.75, 3.75, 3.5, 3.25, 3.25, 3.25, 3.25, 3., 3., 3., 2.75, 2.75, 2.75, 2.75, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5};

    ofstream fit("fitDATA.test.txt"); //-->HERE

    fit << "rr ";
    fit << "landau_width e_landau_width ";
    fit << "landau_mpv e_landau_mpv ";
    fit << "gauss_width e_gauss_width ";
    fit << "int_landau_width e_int_landau_width ";
    fit << "int_landau_mpv e_int_landau_mpv ";
    fit << "int_gauss_width e_int_gauss_width ";
    fit << "fraction e_fraction ";
    fit << "chi2 ndof expected_dedx" << endl; 

    for(int rr = 1; rr < 25; rr++)
    {

        cout << "PROCESSING RR " << rr << endl;

        cout << "   DATA" << endl;

        TH1D * h_dedx_MU_DATA = new TH1D(Form("h_dedx_rr_%d_MU_DATA",rr),"",nbins_pro[rr-1],0,xhigh_pro[rr-1]);
        h_dedx_MU_DATA->Sumw2();
        TH1D * h_dedx_PRO_DATA = new TH1D(Form("h_dedx_rr_%d_PRO_DATA",rr),"",nbins_pro[rr-1],0,xhigh_pro[rr-1]);
        h_dedx_PRO_DATA->Sumw2();

        
        infile[0] -> cd();
        
        TTree * intree_DATA = (TTree*)infile[0]->Get("tree");

        _slice *slc_DATA = nullptr;

        intree_DATA -> SetBranchAddress("slice",&slc_DATA);

        for(int i = 0; i < intree_DATA -> GetEntries(); i ++)
        {
            intree_DATA -> GetEntry(i);

            for(int hit = 0; hit < slc_DATA->_mu._dedx.size(); hit++)
            {
                if(slc_DATA->_mu._rr[hit] >= rr && slc_DATA->_mu._rr[hit] < rr + 1)
                {
                    h_dedx_MU_DATA -> Fill(slc_DATA->_mu._dedx[hit]);
                }
            }

            for (const auto &pro : slc_DATA->_protons)
            {
                for(int hit = 0; hit < pro._dedx.size(); hit ++)
                {
                    if(pro._rr[hit] >= rr && pro._rr[hit] < rr + 1)
                    {
                        h_dedx_PRO_DATA-> Fill(pro._dedx[hit]);
                    }
                }
            }
        }


        TH1D * h_dedx_MU_MC = new TH1D(Form("h_dedx_rr_%d_MU_MC",rr),"",nbins_pro[rr-1],0,xhigh_pro[rr-1]);
        h_dedx_MU_MC->Sumw2();
        TH1D * h_dedx_PRO_MC = new TH1D(Form("h_dedx_rr_%d_PRO_MC",rr),"",nbins_pro[rr-1],0,xhigh_pro[rr-1]);
        h_dedx_PRO_MC->Sumw2();

        int sample_idx = -1;
        for(const auto &s : samples)
        {
            cout << Form("   %s",s.c_str()) << endl;

            sample_idx++;

            infile[sample_idx] -> cd();
        
            TTree * intree = (TTree*)infile[sample_idx]->Get("tree");

            //intree->Print();
            //intree->GetBranch("slice")->Print();
            //intree->GetBranch("slice")->GetClassName();

            _slice *slc = nullptr;

            intree -> SetBranchAddress("slice",&slc);

            for(int i = 0; i < intree -> GetEntries(); i ++)
            {
                intree -> GetEntry(i);

                for(int hit = 0; hit < slc->_mu._dedx.size(); hit ++)
                {
                    if(slc->_mu._rr[hit] >= rr && slc->_mu._rr[hit] < rr + 1)
                    {
                        h_dedx_MU_MC -> Fill(slc->_mu._dedx[hit],POTS[sample_idx]);
                    }
                }

                for (const auto &pro : slc->_protons)
                {
                    for(int hit = 0; hit < pro._dedx.size(); hit ++)
                    {
                        if(pro._rr[hit] >= rr && pro._rr[hit] < rr + 1)
                        {
                            h_dedx_PRO_MC-> Fill(pro._dedx[hit],POTS[sample_idx]);
                            //xvar.push_back(pro._dedx[hit]);
                        }
                    }
                }
            }
        }

        h_dedx_MU_DATA -> Scale(1./h_dedx_MU_DATA->Integral("width"));
        h_dedx_PRO_DATA -> Scale(1./h_dedx_PRO_DATA->Integral("width"));
        h_dedx_MU_MC -> Scale(1./h_dedx_MU_MC->Integral("width")); 
        h_dedx_PRO_MC -> Scale(1./h_dedx_PRO_MC->Integral("width"));

        //h_dedx_MU_DATA->Write(0,TObject::kOverwrite);
        //h_dedx_PRO_DATA->Write(0,TObject::kOverwrite);
        //h_dedx_MU_MC->Write(0,TObject::kOverwrite);
        //h_dedx_PRO_MC->Write(0,TObject::kOverwrite);

        X_HIGH = xhigh_pro[rr -1];
        X_LOW = 0.;

        TFile *dedx_template = TFile::Open("RefCurvesChi2.root","READ");
        TH1D *dedx_range_pro = (TH1D*)dedx_template->Get("dedx_range_pro");
        double dedx_expected = dedx_range_pro -> GetBinContent(dedx_range_pro -> FindBin(rr));

        if(rr == 1){dedx_expected = dedx_expected - 2;}

        BIN_WIDTH = h_dedx_PRO_MC->GetBinWidth(1);

        cout << "START FITTING WITH: RANGE = [" << X_LOW << " " << X_HIGH << "]" << " BIN WIDTH " << BIN_WIDTH << endl;

        // *** LANDAU + 2 GAUSSIANS ***
        /* 
            TF1 *f1 = new TF1(Form("f1_rr%d",rr), fitFunction, X_LOW, X_HIGH, 9);

            f1->SetParNames("LandauWidth", 
                            "LandauMPV", 
                            "ConvSigma",
                            "GaussMean1", 
                            "GaussSigma1", 
                            "GaussMean2", 
                            "GaussWidth2", 
                            "TwoGaussFraction", 
                            "FracLandau" );

            // valori iniziali
            double initial_landau_width = 1.5;
            double initial_gauss_width = 1.;
            double initial_interacting_gauss_sigma = 3.;
            double initial_fraction = 0.5;
            double initial_interacting_gauss_fraction = 0.5;
            f1->SetParameters(initial_landau_width, 
                              dedx_expected, 
                              initial_gauss_width, 
                              interacting_start_value[rr-1], 
                              initial_interacting_gauss_sigma, 
                              interacting_start_value[rr-1], 
                              initial_interacting_gauss_sigma, 
                              initial_interacting_gauss_fraction,
                              initial_fraction);

            cout << "FIT INITIAL PARAMETERS: " << initial_landau_width << " " << dedx_expected << " " << initial_gauss_width << " " << interacting_start_value[rr-1] << " " << initial_interacting_gauss_sigma << " " << interacting_start_value[rr-1] << " " << initial_interacting_gauss_sigma << " " << initial_interacting_gauss_fraction << " " << initial_fraction << endl;

            h_dedx_PRO_MC->Fit(f1,"0"); // L = binned log-likelihood, R = usa il range della TF1

            double landau_width = f1 -> GetParameter(0);
            double landau_mpv = f1 -> GetParameter(1);
            double gauss_width = f1 -> GetParameter(2);
            double interacting_gauss_mean1 = f1 -> GetParameter(3);
            double interacting_gauss_width1 = f1 -> GetParameter(4);
            double interacting_gauss_mean2 = f1 -> GetParameter(5);
            double interacting_gauss_width2 = f1 -> GetParameter(6);
            double interacting_gauss_fraction = f1 -> GetParameter(7);
            double fraction = f1 -> GetParameter(8);
            double norm = 1.;

            //Double_t pars[7] = {landau_width, landau_mpv, gauss_width, interacting_gauss_mean, interacting_gauss_width, fraction, norm};
            Double_t pars[9] = {landau_width, landau_mpv, gauss_width, interacting_gauss_mean1, interacting_gauss_width1, interacting_gauss_mean2, interacting_gauss_width2, interacting_gauss_fraction, fraction};

            Double_t landau_pars[3] = {landau_width,landau_mpv,gauss_width};

            Double_t gauss_pars[5] = {interacting_gauss_mean1, interacting_gauss_width1, interacting_gauss_mean2, interacting_gauss_width2, interacting_gauss_fraction};

        */

        // *** LANDAU + 1 GAUSSIAN ***
        /* 
            TF1 *f1 = new TF1(Form("f1_rr%d",rr), fitFunction, X_LOW, X_HIGH, 6);

            f1->SetParNames("LandauWidth", 
                            "LandauMPV", 
                            "ConvSigma",
                            "GaussMean", 
                            "GaussSigma", 
                            "FracLandau" );

            // valori iniziali
            double initial_landau_width = 1.5;
            double initial_gauss_width = 1.;
            double initial_interacting_gauss_sigma = 3.;
            double initial_fraction = 0.5;
            f1->SetParameters(initial_landau_width, 
                              dedx_expected, 
                              initial_gauss_width, 
                              interacting_start_value[rr-1], 
                              initial_interacting_gauss_sigma, 
                              initial_fraction);

            cout << "FIT INITIAL PARAMETERS: " << initial_landau_width << " " << dedx_expected << " " << initial_gauss_width << " " << interacting_start_value[rr-1] << " " << initial_interacting_gauss_sigma << " " << initial_fraction << endl;

            h_dedx_PRO_MC->Fit(f1,"0"); // L = binned log-likelihood, R = usa il range della TF1

            double landau_width = f1 -> GetParameter(0);
            double landau_mpv = f1 -> GetParameter(1);
            double gauss_width = f1 -> GetParameter(2);
            double interacting_gauss_mean1 = f1 -> GetParameter(3);
            double interacting_gauss_width1 = f1 -> GetParameter(4);
            double fraction = f1 -> GetParameter(5);
            double norm = 1.;

            //Double_t pars[6] = {landau_width, landau_mpv, gauss_width, interacting_gauss_mean, interacting_gauss_width, fraction};

            Double_t landau_pars[3] = {landau_width,landau_mpv,gauss_width};

            Double_t gauss_pars[2] = {interacting_gauss_mean1, interacting_gauss_width1};

        */

        // *** 2 LANDAU ***

            TF1 *f1 = new TF1(Form("f1_rr%d",rr), fitFunction, X_LOW, X_HIGH, 7);

            f1->SetParNames("LandauWidth1", 
                            "LandauMPV1", 
                            "ConvSigma1",
                            "LandauWidth2", 
                            "LandauMPV2", 
                            "ConvSigma2", 
                            "FracLandau" );

            // valori iniziali
            double initial_landau_width = 1.;
            double initial_gauss_width = 1.;
            double initial_interacting_landau_width = 1.5;
            double initial_interacting_gauss_width = 3.;
            double initial_fraction = 0.5;
            f1->SetParameters(initial_landau_width, 
                              dedx_expected, 
                              initial_gauss_width, 
                              initial_interacting_landau_width,
                              interacting_start_value[rr-1], 
                              initial_interacting_gauss_width, 
                              initial_fraction);

            cout << "FIT INITIAL PARAMETERS: " << initial_landau_width << " " << dedx_expected << " " << initial_gauss_width << " " << initial_interacting_landau_width << " " << interacting_start_value[rr-1] << " " << initial_interacting_gauss_width << " " << initial_fraction << endl;

            f1->SetParLimits(0, 0.05, 2.);
            f1->SetParLimits(1, dedx_expected-2.5, dedx_expected+2.5);
            f1->SetParLimits(2, 0.05, 5.);
            f1->SetParLimits(3, 0.05, 5.);
            f1->SetParLimits(4, 0., dedx_expected);
            f1->SetParLimits(5, 0.05, 5.);
            f1->SetParLimits(6, 0., 1.);


            h_dedx_PRO_DATA->Fit(f1,"0R"); // L = binned log-likelihood, R = usa il range della TF1 //-->HERE

            double landau_width1 = f1 -> GetParameter(0);
            double landau_mpv1 = f1 -> GetParameter(1);
            double gauss_width1 = f1 -> GetParameter(2);
            double landau_width2 = f1 -> GetParameter(3);
            double landau_mpv2 = f1 -> GetParameter(4);
            double gauss_width2 = f1 -> GetParameter(5);
            double fraction = f1 -> GetParameter(6);
            double norm = 1.;

            fit << rr << " ";
            fit << landau_width1 << " " << f1 -> GetParError(0) << " ";
            fit << landau_mpv1 << " " << f1 -> GetParError(1) << " "; 
            fit << gauss_width1 << " " << f1 -> GetParError(2) << " "; 
            fit << landau_width2 << " " << f1 -> GetParError(3) << " ";
            fit << landau_mpv2 << " " << f1 -> GetParError(4) << " "; 
            fit << gauss_width2 << " " << f1 -> GetParError(5) << " ";
            fit << fraction << " " << f1 -> GetParError(6) << " ";
            fit << f1->GetChisquare() << " " << f1->GetNDF() << " " << dedx_expected << endl;

            Double_t pars[7] = {landau_width1, landau_mpv1, gauss_width1, landau_width2, landau_mpv2, gauss_width2, fraction};

            Double_t landau_pars1[3] = {landau_width1,landau_mpv1,gauss_width1};

            Double_t landau_pars2[3] = {landau_width2,landau_mpv2,gauss_width2};


        std::vector<double> x_points;
        std::vector<double> y_points;
        std::vector<double> y_points_landau;
        std::vector<double> y_points_gauss;

        double landau_norm1 = get_langaus_norm(landau_pars1,X_LOW,X_HIGH,1000);
        double landau_norm2 = get_langaus_norm(landau_pars2,X_LOW,X_HIGH,1000);

        for(int xi=0; xi<1000; xi++)
        {
            double xmin = X_LOW;
            double xmax = X_HIGH;
        
            double x = xmin + xi * (xmax-xmin)/1000;

            x_points.push_back(x);
            y_points.push_back( pdf_proj_light(&x, pars, landau_norm1, landau_norm2) );

            y_points_landau.push_back( fraction * landau_proj_light(&x, landau_pars1, landau_norm1) );
            
            //y_points_gauss.push_back(
            //    (1-fraction) * gauss_proj(&x,
            //        gauss_pars,
            //        (1-fraction) * norm,
            //        h_dedx_PRO_MC->GetBinWidth(1)));
            
            y_points_gauss.push_back( (1-fraction) * landau_proj_light(&x, landau_pars2, landau_norm2) );
        }

        TGraph * fit_function_plot = new TGraph(1000,x_points.data(),y_points.data());
        TGraph * fit_landau_plot = new TGraph(1000,x_points.data(),y_points_landau.data());
        TGraph * fit_gauss_plot = new TGraph(1000,x_points.data(),y_points_gauss.data());

        TCanvas *c_fit = new TCanvas(Form("c_fit_%d",rr), "fit");
        c_fit->cd();
        h_dedx_PRO_DATA -> Draw("P"); //-->HERE
        fit_function_plot -> SetLineColor(kRed);
        fit_function_plot -> SetMarkerColor(kRed);
        fit_function_plot -> SetLineWidth(2);
        fit_function_plot -> Draw("SAME");
        fit_landau_plot -> SetLineColor(kOrange);
        fit_landau_plot -> SetMarkerColor(kOrange);
        fit_landau_plot -> SetLineStyle(kDashed);
        fit_landau_plot -> Draw("SAME");
        fit_gauss_plot -> SetLineColor(kViolet);
        fit_gauss_plot -> SetMarkerColor(kViolet);
        fit_gauss_plot -> SetLineStyle(kDashed);
        fit_gauss_plot -> Draw("SAME");

        outfile -> cd();
        c_fit -> Write();

        h_dedx_PRO_DATA -> Write(Form("hist_rr_%d",rr));
        fit_function_plot -> Write(Form("tf1_rr_%d",rr));

        //h_dedx_PRO_MC -> Write();
    }
}