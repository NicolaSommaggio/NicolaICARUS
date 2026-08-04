#include "slice_struct.h"
#include "fitLangaus_unbinned.C"
#include "fitLangaulight.C"

double p_gauss_pdf(double *x, double *par)
{
    const Double_t sqrt2 = TMath::Sqrt(2.0);

    // *** 1 GAUSSIAN ***
    //double gnorm = 0.5 * (
    //        TMath::Erf((X_HIGH - par[0])/(sqrt2*par[1]))
    //      - TMath::Erf((X_LOW - par[0])/(sqrt2*par[1]))
    //    );
    return TMath::Gaus(x[0], par[0], par[1], 1);
}

double medianDedxRRcut(const std::vector<double>& dedx, const std::vector<double>& rr, double rr_cut = 5.0)
{
    std::vector<double> selected;
    selected.reserve(dedx.size());

    for (size_t i = 0; i < dedx.size(); ++i)
    { 
        if (rr[i] < rr_cut)selected.push_back(dedx[i]);
    }

    if (selected.empty())return -1;

    std::sort(selected.begin(), selected.end());

    size_t n = selected.size();
    if(n % 2 == 0){return 0.5 * (selected[n/2 - 1] + selected[n/2]);}
    else{return selected[n/2];}
}

std::string DATASET = "MC";
std::string PARTICLE = "pro";

void PLOT_VARIATIONS()
{
    gROOT->SetBatch(kTRUE);

    std::vector<std::string> samples = {"MC","OFFBEAM"};
    //std::vector<double> POTS = {2.43007e+20, 2.11e+20};
    std::vector<double> POTS = {2.43007, 2.11};
    TFile * outfile = new TFile(Form("OUT_VARIATIONS_STANDARD_RUN2_%s_%s.root",PARTICLE.c_str(),DATASET.c_str()),"RECREATE"); 
            
    TFile * infile_MC = TFile::Open("varMC_STANDARD_CUT.root","READ");  
    TFile * infile_OFFBEAM = TFile::Open("varMC_OFFBEAM_CUT.root","READ");        
    TFile * infile_DATA = TFile::Open("varMC_DATA_CUT.root","READ");

    std::vector<TFile*> infile = {infile_DATA, infile_MC, infile_OFFBEAM};

    ofstream fit(Form("fit%s_%s.txt",DATASET.c_str(), PARTICLE.c_str())); 

    TH1D *mediana_mc = new TH1D("mediana_mc","",100,0,30);
    TH1D *mediana_data = new TH1D("mediana_data","",100,0,30);
    TH1D *mediana_mc_pro = new TH1D("mediana_mc_pro","",100,0,30);
    TH1D *mediana_data_pro = new TH1D("mediana_data_pro","",100,0,30);

    TFile *files_cm[2] = {infile_MC,infile_OFFBEAM};
    for(TFile* &in : files_cm)
    {
        TTree * intree = (TTree*)in->Get("tree");

        _slice *slc = nullptr;

        intree -> SetBranchAddress("slice",&slc);

        //double POT_temp = (in == infile_MC) ? 2.43007e+20 : 2.11e+20;
        double POT_temp = (in == infile_MC) ? 2.43007 : 2.11;

        for(int i = 0; i < intree -> GetEntries(); i ++)
        {
            intree -> GetEntry(i);

            mediana_mc->Fill(medianDedxRRcut(slc->_mu._dedx,slc->_mu._rr,5.), POT_temp);

            for (const auto &pro : slc->_protons)
            {
                mediana_mc_pro->Fill(medianDedxRRcut(pro._dedx,pro._rr,5.), POT_temp);
            }
        }

        
    } 

    
    TTree * intree = (TTree*)infile_DATA->Get("tree");

    _slice *slc = nullptr;

    intree -> SetBranchAddress("slice",&slc);

    for(int i = 0; i < intree -> GetEntries(); i ++)
    {
        intree -> GetEntry(i);

        mediana_data->Fill(medianDedxRRcut(slc->_mu._dedx,slc->_mu._rr,5.));

        for (const auto &pro : slc->_protons)
        {
            mediana_data_pro->Fill(medianDedxRRcut(pro._dedx,pro._rr,5.));
        }
    }
 

    fit << "rr ";
    fit << "landau_width e_landau_width ";
    fit << "landau_mpv e_landau_mpv ";
    fit << "gauss_width e_gauss_width ";
    fit << "int_landau_width e_int_landau_width ";
    fit << "int_landau_mpv e_int_landau_mpv ";
    fit << "int_gauss_width e_int_gauss_width ";
    fit << "fraction e_fraction ";
    fit << "chi2 ndof expected_dedx status cov_matrix_quality" << endl; 

    int rrs[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 21, 25};
    //for(int rr = 1; rr < 25; rr++)
    for(int &rr : rrs)
    {
        std::vector<double> dedx_thisbin_data_mu;
        std::vector<double> dedx_thisbin_mc_mu;
        std::vector<double> dedx_thisbin_data_pro;
        std::vector<double> dedx_thisbin_mc_pro;

        std::vector<double> weights_thisbin_data_mu;
        std::vector<double> weights_thisbin_mc_mu;
        std::vector<double> weights_thisbin_data_pro;
        std::vector<double> weights_thisbin_mc_pro;

        cout << "PROCESSING RR " << rr << endl;

        cout << "   DATA" << endl;

        TH1D * h_dedx_MU_DATA = new TH1D(Form("h_dedx_rr_%d_MU_DATA",rr),"",100,0,30.);
        h_dedx_MU_DATA->Sumw2();
        TH1D * h_dedx_PRO_DATA = new TH1D(Form("h_dedx_rr_%d_PRO_DATA",rr),"",100,0,30.);
        h_dedx_PRO_DATA->Sumw2();
        
        infile[0] -> cd();
        
        TTree * intree_DATA = (TTree*)infile[0]->Get("tree");

        _slice *slc_DATA = nullptr;

        intree_DATA -> SetBranchAddress("slice",&slc_DATA);

        for(int i = 0; i < intree_DATA -> GetEntries(); i ++)
        {
            intree_DATA -> GetEntry(i);

            if(medianDedxRRcut(slc_DATA->_mu._dedx,slc_DATA->_mu._rr,5.)>3.3)
            {
                for(int hit = 0; hit < slc_DATA->_mu._dedx.size(); hit++)
                {
                    if(slc_DATA->_mu._rr[hit] >= rr && slc_DATA->_mu._rr[hit] < rr + 1)
                    {
                        h_dedx_MU_DATA -> Fill(slc_DATA->_mu._dedx[hit]);
                        dedx_thisbin_data_mu.push_back(slc_DATA->_mu._dedx[hit]);
                        weights_thisbin_data_mu.push_back(1.);
                    }
                }
            }
            
            for (const auto &pro : slc_DATA->_protons)
            {
                for(int hit = 0; hit < pro._dedx.size(); hit ++)
                {
                    if(pro._rr[hit] >= rr && pro._rr[hit] < rr + 1)
                    {
                        h_dedx_PRO_DATA-> Fill(pro._dedx[hit]);
                        dedx_thisbin_data_pro.push_back(pro._dedx[hit]);
                        weights_thisbin_data_pro.push_back(1.);
                    }
                }
            }
        }


        TH1D * h_dedx_MU_MC = new TH1D(Form("h_dedx_rr_%d_MU_MC",rr),"",100,0,30.);
        h_dedx_MU_MC->Sumw2();
        TH1D * h_dedx_PRO_MC = new TH1D(Form("h_dedx_rr_%d_PRO_MC",rr),"",100,0,30.);
        h_dedx_PRO_MC->Sumw2();

        int sample_idx = 0;
        for(const auto &s : samples)
        {
            cout << Form("   %s",s.c_str()) << endl;

            sample_idx++;

            infile[sample_idx] -> cd();
        
            TTree * intree = (TTree*)infile[sample_idx]->Get("tree");

            _slice *slc = nullptr;

            intree -> SetBranchAddress("slice",&slc);

            for(int i = 0; i < intree -> GetEntries(); i ++)
            {
                intree -> GetEntry(i);

                std::vector<double> temp_dedx;
                std::vector<double> temp_rr;

                if(medianDedxRRcut(slc->_mu._dedx,slc->_mu._rr,5.)>3.3)
                {

                    for(int hit = 0; hit < slc->_mu._dedx.size(); hit ++)
                    {
                        if(slc->_mu._rr[hit] >= rr && slc->_mu._rr[hit] < rr + 1)
                        {
                            h_dedx_MU_MC -> Fill(slc->_mu._dedx[hit],POTS[sample_idx]);
                            dedx_thisbin_mc_mu.push_back(slc->_mu._dedx[hit]);
                            weights_thisbin_mc_mu.push_back(POTS[sample_idx]);
                        }
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
                            dedx_thisbin_mc_pro.push_back(pro._dedx[hit]);
                            weights_thisbin_mc_pro.push_back(POTS[sample_idx]);
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

        TFile *dedx_template = TFile::Open("RefCurvesChi2.root","READ");
        TH1D *dedx_range = (TH1D*)dedx_template->Get(Form("dedx_range_%s",PARTICLE.c_str()));
        double dedx_expected = dedx_range -> GetBinContent(dedx_range -> FindBin(rr));

        //TF1 *gauss_tf1 = new TF1("gauss_tf1", p_gauss_pdf, X_LOW, X_HIGH, 2);
        //gauss_tf1 -> SetParameters(dedx_expected, 1.5);
        //h_dedx_MU_MC -> Fit(gauss_tf1, "R");

        TH1D *h_fit;
        std::vector<double> x_fit;
        std::vector<double> w_fit;
        if(PARTICLE == "mu")
        {
            h_fit = (DATASET == "DATA") ? h_dedx_MU_DATA : h_dedx_MU_MC;
            x_fit = (DATASET == "DATA") ? dedx_thisbin_data_mu : dedx_thisbin_mc_mu;
            w_fit = (DATASET == "DATA") ? weights_thisbin_data_mu : weights_thisbin_mc_mu;
        }
        else
        {
            h_fit = (DATASET == "DATA") ? h_dedx_PRO_DATA : h_dedx_PRO_MC;
            x_fit = (DATASET == "DATA") ? dedx_thisbin_data_pro : dedx_thisbin_mc_pro;
            w_fit = (DATASET == "DATA") ? weights_thisbin_data_pro : weights_thisbin_mc_pro;
        }
        

        if(PARTICLE == "mu")
        {
            X_LOW = (rr == 1) ? 4.5 : h_fit -> GetMean() - 1.3 * h_fit -> GetStdDev();
            //X_LOW = h_fit -> GetMean() - 1.3 * h_fit -> GetStdDev();
            X_HIGH = (rr == 1) ? 30. : 20.;
        }
        else{X_LOW = 0.; X_HIGH = 20;}
        
        BIN_WIDTH = h_fit->GetBinWidth(1);

        cout << "START FITTING WITH: RANGE = [" << X_LOW << " " << X_HIGH << "]" << " BIN WIDTH " << BIN_WIDTH << endl;

        FitResults fitres;
        if(PARTICLE == "mu")fitres = fitSingleLangau_unbinned(h_fit, x_fit, w_fit, X_LOW, X_HIGH, outfile, rr, dedx_range);
        else 
        {
            double best_chi2 = 10000;
            double best_xlow = 0; 
            for(double xlow = 0; xlow <= 1.5; xlow = xlow + 0.25)
            {
                fitres = fitDoubleLangau_unbinned(h_fit, x_fit, w_fit, X_LOW, X_HIGH, outfile, rr, dedx_range);
                if(fitres.chi2 < best_chi2 && fitres.status == 0 && fitres.covqual >= 2)
                {
                    best_xlow = xlow;
                }

                cout << "xlow: " << xlow << " chi2: " << fitres.chi2 << " status: " << fitres.status << " covqual: " << fitres.covqual << endl;
            }

            fitres = fitDoubleLangau_unbinned(h_fit, x_fit, w_fit, best_xlow, X_HIGH, outfile, rr, dedx_range);

        }

            fit << rr << " ";
            fit << fitres.landau_width << " " << fitres.e_landau_width << " ";
            fit << fitres.landau_mpv << " " << fitres.e_landau_mpv << " ";
            fit << fitres.gauss_width << " " << fitres.e_gauss_width << " ";
            fit << fitres.int_landau_width << " " << fitres.e_int_landau_width << " ";
            fit << fitres.int_landau_mpv << " " << fitres.e_int_landau_mpv << " ";
            fit << fitres.int_gauss_width << " " << fitres.e_int_gauss_width << " ";
            fit << fitres.fraction << " " << fitres.e_fraction << " ";
            fit << fitres.chi2 << " " << fitres.ndf << " " << dedx_range->GetBinContent(dedx_range->FindBin(rr)) << " " << fitres.status << " " << fitres.covqual << endl;

    }

    outfile -> cd();
    mediana_data -> Scale(1./mediana_data->Integral("width"));
    mediana_mc -> Scale(1./mediana_mc->Integral("width"));
    mediana_data -> Write();
    mediana_mc -> Write();

    mediana_data_pro -> Scale(1./mediana_data_pro->Integral("width"));
    mediana_mc_pro -> Scale(1./mediana_mc_pro->Integral("width"));
    mediana_data_pro -> Write();
    mediana_mc_pro -> Write();
}


void create_spline()
{
    ifstream spline_table("spline_table_new.txt");

    TGraph *spline_points = new TGraph();

    std::string line;
    while (std::getline(spline_table, line))
    {
        std::stringstream ss(line);
        double x, y;
        ss >> x >> y;
        spline_points->SetPoint(spline_points->GetN(), x, y);
    }

    TSpline3 *spline = new TSpline3("dedx_corr_spline", spline_points);

    TCanvas *c = new TCanvas();
    c->cd();
    spline_points->Draw("AP");
    spline->SetLineColor(kRed);
    spline->Draw("same");
    c->Draw();

    TFile *outfile = new TFile("spline_file_new.root","RECREATE");
    outfile -> cd();
    spline -> Write();
}