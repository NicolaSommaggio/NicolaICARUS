#include "ReadTree.C"
#include "fitLangaulight.C"
//#include "daughtersInfo.h"


//#include "MacroPoppiLandau.C"

// /storage/gpfs_data/icarus/plain/data/prod/run2-v09_84_00_01-202403-cnaf/concat-caf/control-sample-100

int findBinMpv(TH1D* h)
{
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
    return binMax;
}

void macroTPC(std::string dati_o_mc, double mu_len_cond, double rr_low, double rr_high){


 //gROOT->ProcessLine(".L libdaughtersInfo.so");
    EventsData dat = load_data(dati_o_mc , "muon", "Np", "full");//QUI

    cout << dat.tree->GetEntries() << " total tracks" << endl;

    TH1D *muon_dEdx_EE = new TH1D("dEdx_EE", "dEdx_EE", 100, 0., 10.);
    TH1D *muon_dEdx_EW = new TH1D("dEdx_EW", "dEdx_EW", 100, 0., 10.);
    TH1D *muon_dEdx_WE = new TH1D("dEdx_WE", "dEdx_WE", 100, 0., 10.);
    TH1D *muon_dEdx_WW = new TH1D("dEdx_WW", "dEdx_WW", 100, 0., 10.);

    std::vector<double> v_dEdx_EE;
    std::vector<double> v_dEdx_EW;
    std::vector<double> v_dEdx_WE;
    std::vector<double> v_dEdx_WW;

    //TH1D *pitch = new TH1D("pitch", "",  200, 0, 10);
    TH1D *pitch = new TH1D("pitch", "",  80, 0, 6);

    for(int track=0; track<dat.tree->GetEntries(); track++)
    {
        dat.tree->GetEntry(track);

        if(dat.track.len_reco > mu_len_cond)
        {
            for(int hit=0; hit<int(dat.track.rr->size()); hit++)
            {                
                if(dat.track.rr->at(hit)> rr_low && dat.track.rr->at(hit)< rr_high)
                {
                    if(dat.track.hitx->at(hit)<-210.)
                    {
                        muon_dEdx_EE->Fill(dat.track.dE->at(hit));
                        v_dEdx_EE.push_back(dat.track.dE->at(hit));
                    }
                    else if(dat.track.hitx->at(hit)>-210.&&dat.track.hitx->at(hit)<0.)
                    {
                        muon_dEdx_EW->Fill(dat.track.dE->at(hit));
                        v_dEdx_EW.push_back(dat.track.dE->at(hit));
                    }
                    else if(dat.track.hitx->at(hit)>0.&&dat.track.hitx->at(hit)<210.)
                    {
                        muon_dEdx_WE->Fill(dat.track.dE->at(hit));
                        v_dEdx_WE.push_back(dat.track.dE->at(hit));
                    }
                    else if(dat.track.hitx->at(hit)>210.)
                    {
                        muon_dEdx_WW->Fill(dat.track.dE->at(hit));
                        v_dEdx_WW.push_back(dat.track.dE->at(hit));
                    }
                }    
            }
        }//condition on L_mu
    }


    for(int track=0; track<dat.tree->GetEntries(); track++)
    {
        dat.tree->GetEntry(track);

        for(int hit=0; hit<int(dat.track.rr->size()); hit++)
        { 
            pitch->Fill(dat.track.pitch->at(hit));
        }

    }

    
    
    TH2D *dedx_vs_range = new TH2D("dedx_vs_range", "", 500, 0, 25, 300, 0, 30);

    for(int track=0; track<dat.tree->GetEntries(); track++)
    {
        dat.tree->GetEntry(track);

            for(int hit=0; hit<int(dat.track.rr->size()); hit++)
            { 
                if(dat.track.rr->at(hit)<= 25)
                {
                    dedx_vs_range -> Fill(dat.track.rr->at(hit), dat.track.dE->at(hit));
            
                }
            }
    }
    
   
    std::string dataType;
    if (dati_o_mc == "dati"){dataType="DATI";}
    if (dati_o_mc == "mc"){dataType = "MC";}
    if (dati_o_mc == "dati2d"){dataType = "DATI_2D";}
    if (dati_o_mc == "dati2d10"){dataType= "DATI_2D_10";}
    if (dati_o_mc == "dati10"){dataType= "DATI_10";}

    

    TH1D *muon_dEdx_EE_noNorm = (TH1D*)muon_dEdx_EE->Clone("muon_dEdx_EE_noNorm");
    TH1D *muon_dEdx_EW_noNorm = (TH1D*)muon_dEdx_EW->Clone("muon_dEdx_EW_noNorm");
    TH1D *muon_dEdx_WE_noNorm = (TH1D*)muon_dEdx_WE->Clone("muon_dEdx_WE_noNorm");
    TH1D *muon_dEdx_WW_noNorm = (TH1D*)muon_dEdx_WW->Clone("muon_dEdx_WW_noNorm");


    
    //muon_dEdx_EE->Scale(1./muon_dEdx_EE->GetBinContent(findBinMpv(muon_dEdx_EE)));
    //muon_dEdx_EW->Scale(1./muon_dEdx_EW->GetBinContent(findBinMpv(muon_dEdx_EW)));
    //muon_dEdx_WE->Scale(1./muon_dEdx_WE->GetBinContent(findBinMpv(muon_dEdx_WE)));
    //muon_dEdx_WW->Scale(1./muon_dEdx_WW->GetBinContent(findBinMpv(muon_dEdx_WW)));

    muon_dEdx_EE->Scale(1./muon_dEdx_EE->Integral());
    muon_dEdx_EW->Scale(1./muon_dEdx_EW->Integral());
    muon_dEdx_WE->Scale(1./muon_dEdx_WE->Integral());
    muon_dEdx_WW->Scale(1./muon_dEdx_WW->Integral());

    TFile *outTPC = TFile::Open("outTPC_dati_1d_vs_2d.root", "UPDATE");
    TDirectory *d = (TDirectory*)outTPC->Get(dataType.c_str()); 
    d->cd();

    dedx_vs_range -> Write(0, TObject::kOverwrite);

    muon_dEdx_EE->Write(0, TObject::kOverwrite);
    muon_dEdx_EW->Write(0, TObject::kOverwrite);
    muon_dEdx_WE->Write(0, TObject::kOverwrite);
    muon_dEdx_WW->Write(0, TObject::kOverwrite);

    pitch->Scale(1./pitch->Integral());
    pitch->Write(0, TObject::kOverwrite);


    //std::vector<double> resEE = fitLangau(muon_dEdx_EE_noNorm, outTPC, "muon_dEdx_EE", dataType);
    //std::vector<double> resEW = fitLangau(muon_dEdx_EW_noNorm, outTPC, "muon_dEdx_EW", dataType);
    //std::vector<double> resWE = fitLangau(muon_dEdx_WE_noNorm, outTPC, "muon_dEdx_WE", dataType);
    //std::vector<double> resWW = fitLangau(muon_dEdx_WW_noNorm, outTPC, "muon_dEdx_WW", dataType);

    FitOutput fitOut_EE = fitLangaulight_unbinned(muon_dEdx_EE_noNorm, v_dEdx_EE, 1.25, 6. , d, "langau", 0);
    FitOutput fitOut_EW = fitLangaulight_unbinned(muon_dEdx_EW_noNorm, v_dEdx_EW, 1.25, 6. , d, "langau", 1);
    FitOutput fitOut_WE = fitLangaulight_unbinned(muon_dEdx_WE_noNorm, v_dEdx_WE, 1.25, 6. , d, "langau", 2);
    FitOutput fitOut_WW = fitLangaulight_unbinned(muon_dEdx_WW_noNorm, v_dEdx_WW, 1.4, 6. , d, "langau", 3);

    std::vector<double> resEE;
    std::vector<double> resEW;
    std::vector<double> resWE;
    std::vector<double> resWW;

    for(int i=0; i<7; i++)
    {
        if(i==6)
        {
            resEE.push_back(fitOut_EE.fitResults[9]);
            resEW.push_back(fitOut_EW.fitResults[9]);
            resWE.push_back(fitOut_WE.fitResults[9]);
            resWW.push_back(fitOut_WW.fitResults[9]);
        }
        else
        {
            resEE.push_back(fitOut_EE.fitResults[i]);
            resEW.push_back(fitOut_EW.fitResults[i]);
            resWE.push_back(fitOut_WE.fitResults[i]);
            resWW.push_back(fitOut_WW.fitResults[i]);
        }
    }

    cout <<  mu_len_cond << " " << rr_low << " " << rr_high << endl;
    cout << "#EE\tEW\tWE\tWW" << endl;
    cout << "#mpv\tmpv_err\tsigmaG\tsigmaG_err\tsigmaL\tsigmaL_err\tchi2" << endl; 

    cout << resEE[0] << "\t" << resEE[1] << "\t" << resEE[2] << "\t" << resEE[3] << "\t" << resEE[4] << "\t" << resEE[5] << "\t" << resEE[6] << endl;
    cout << resEW[0] << "\t" << resEW[1] << "\t" << resEW[2] << "\t" << resEW[3] << "\t" << resEW[4] << "\t" << resEW[5] << "\t" << resEW[6] << endl;
    cout << resWE[0] << "\t" << resWE[1] << "\t" << resWE[2] << "\t" << resWE[3] << "\t" << resWE[4] << "\t" << resWE[5] << "\t" << resWE[6] << endl;
    cout << resWW[0] << "\t" << resWW[1] << "\t" << resWW[2] << "\t" << resWW[3] << "\t" << resWW[4] << "\t" << resWW[5] << "\t" << resWW[6] << endl;
   
    //ofstream mean("meanMC.txt"); // QUI

    
    cout << mu_len_cond << " "<<  rr_low << " " << rr_high << endl;
    cout << "#mean\tsigma_mean" << endl;
    cout << muon_dEdx_EE_noNorm->GetMean() << " " <<  muon_dEdx_EE_noNorm->GetMeanError() << " " << muon_dEdx_EE_noNorm->GetEntries() << endl;
    cout << muon_dEdx_EW_noNorm->GetMean() << " " <<  muon_dEdx_EW_noNorm->GetMeanError() << " " << muon_dEdx_EW_noNorm->GetEntries() << endl;
    cout << muon_dEdx_WE_noNorm->GetMean() << " " <<  muon_dEdx_WE_noNorm->GetMeanError() << " " << muon_dEdx_WE_noNorm->GetEntries() << endl;
    cout << muon_dEdx_WW_noNorm->GetMean() << " " <<  muon_dEdx_WW_noNorm->GetMeanError() << " " << muon_dEdx_WW_noNorm->GetEntries() << endl;  


    TCanvas *c = new TCanvas("TPCs_dEdx_muon", "TPC_dEdx_muon");
    gStyle->SetOptStat(0);
    muon_dEdx_EW->Draw("hist");
    muon_dEdx_EW->GetXaxis()->SetTitle("dE/dx [MeV/cm]");
    muon_dEdx_EW->GetYaxis()->SetTitle("counts (area normalized)");
    muon_dEdx_EW->SetLineWidth(2);
    muon_dEdx_EE->Draw("same hist");
    muon_dEdx_EE->SetLineColor(kRed);
    muon_dEdx_EE->SetLineWidth(2);
    muon_dEdx_WE->Draw("same hist");
    muon_dEdx_WE->SetLineColor(kMagenta);
    muon_dEdx_WE->SetLineWidth(2);
    muon_dEdx_WW->Draw("same hist");
    muon_dEdx_WW->SetLineColor(kGreen);
    muon_dEdx_WW->SetLineWidth(2);
    TLegend *l = new TLegend();
    l->AddEntry(muon_dEdx_EE, "EE");
    l->AddEntry(muon_dEdx_EW, "EW");
    l->AddEntry(muon_dEdx_WE, "WE");
    l->AddEntry(muon_dEdx_WW, "WW");

    TLatex* latex = new TLatex();
    latex->SetNDC();              
    latex->SetTextSize(0.03);     
    TString text1 = Form("EE: mpv = %.4f #pm %.4f", resEE[0], resEE[1]);
    TString text2 = Form("EE: sigmaG = %.4f #pm %.4f", resEE[2], resEE[3]);
    TString text3 = Form("EE: sigmaL = %.4f #pm %.4f", resEE[4], resEE[5]);

    TString text4 = Form("EW: mpv = %.4f #pm %.4f", resEW[0], resEW[1]);
    TString text5 = Form("EW: sigmaG = %.4f #pm %.4f", resEW[2], resEW[3]);
    TString text6 = Form("EW: sigmaL = %.4f #pm %.4f", resEW[4], resEW[5]);

    TString text7 = Form("WE: mpv = %.4f #pm %.4f", resWE[0], resWE[1]);
    TString text8 = Form("WE: sigmaG = %.4f #pm %.4f", resWE[2], resWE[3]);
    TString text9 = Form("WE: sigmaL = %.4f #pm %.4f", resWE[4], resWE[5]);

    TString text10 = Form("WW: mpv = %.4f #pm %.4f", resWW[0], resWW[1]);
    TString text11 = Form("WW: sigmaG = %.4f #pm %.4f", resWW[2], resWW[3]);
    TString text12 = Form("WW: sigmaL = %.4f #pm %.4f", resWW[4], resWW[5]);

    latex->DrawLatex(0.55, 0.8, text1);
    latex->DrawLatex(0.55, 0.77, text2);
    latex->DrawLatex(0.55, 0.74, text3);
    latex->DrawLatex(0.55, 0.71, text4);
    latex->DrawLatex(0.55, 0.68, text5);
    latex->DrawLatex(0.55, 0.65, text6);
    latex->DrawLatex(0.55, 0.62, text7);
    latex->DrawLatex(0.55, 0.59, text8);
    latex->DrawLatex(0.55, 0.56, text9);
    latex->DrawLatex(0.55, 0.53, text10);
    latex->DrawLatex(0.55, 0.50, text11);
    latex->DrawLatex(0.55, 0.47, text12);

    l->Draw();
    c->Write(0, TObject::kOverwrite);

    outTPC->Close();

    delete c;


}