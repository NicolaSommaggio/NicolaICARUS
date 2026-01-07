#include "ReadTree.C"
#include "fitLangaulight.C"
#include <cmath>

// root -l -q 'macroTPC2.C("mc", 1.4, 6.)'
void macroTPC2(std::string dati_o_mc){

    EventsData dat = load_data(dati_o_mc , "muon", "Np", "full");

    cout << dat.tree->GetEntries() << " total tracks" << endl;

    TH2D *dedx_vs_range = new TH2D("dedx_vs_range", "", 500, 0, 25, 300, 0, 30);

    TH1D *dedx_rr_100_150 = new TH1D("dedx_rr_100_150", "", 200, 0, 10);
    TH1D *dedx_rr_20_25 = new TH1D("dedx_rr_20_25", "", 200, 0, 10);
    TH1D *dedx_rr_100_150_corretto = new TH1D("dedx_rr_100_150_corretto", "", 200, 0, 10);

    for(int track=0; track<dat.tree->GetEntries(); track++)
    {
        dat.tree->GetEntry(track);

            for(int hit=0; hit<int(dat.track.rr->size()); hit++)
            { 
                if(dat.track.rr->at(hit)> 16. && dat.track.rr->at(hit)< 25.)
                {
                    dedx_vs_range -> Fill(dat.track.rr->at(hit), dat.track.dE->at(hit));
                }

                if(dat.track.len_reco>150.)
                {
                    if(dat.track.rr->at(hit)> 100. && dat.track.rr->at(hit)<150.)
                    {
                        dedx_rr_100_150 -> Fill(dat.track.dE->at(hit));
                        dedx_rr_100_150_corretto->Fill(dat.track.dE->at(hit)*1.045027006);
                    }
                }

                if(dat.track.rr->at(hit)> 20. && dat.track.rr->at(hit)<25.)
                {
                    dedx_rr_20_25 -> Fill(dat.track.dE->at(hit));
                }
            }
    }


    cout << dedx_rr_100_150->GetMean() << endl;
    

    TFile * outTPC = TFile::Open("outTPC.root", "UPDATE");
    std::string folder;
    if(dati_o_mc=="mc"){folder="MC";}
    else if(dati_o_mc=="dati"){folder="DATI";}

    TDirectory *d = (TDirectory*)outTPC->Get(folder.c_str());
    d->cd();

    dedx_rr_100_150->Scale(1./dedx_rr_100_150->Integral());
    dedx_rr_20_25->Scale(1./dedx_rr_20_25->Integral());
    dedx_rr_100_150_corretto->Scale(1./dedx_rr_100_150_corretto->Integral());

    dedx_rr_100_150->Write(0,TObject::kOverwrite);
    dedx_rr_20_25->Write(0,TObject::kOverwrite);

    if(folder=="MC"){dedx_rr_100_150_corretto->Write(0,TObject::kOverwrite);}

    outTPC->Close();
    delete outTPC;
    
}

ofstream med9435;
void look_at_median()
{

    TFile *f = new TFile("median_new.root","RECREATE");
    //std::array<std::string,2> dataTypes = {"Full9435_1d","Full9435_2d"};
    std::array<std::string,2> dataTypes = {"data_1d_run2_full_v09_89","dati2d"};


    for(auto dati_o_mc : dataTypes)
    {
        //if(dati_o_mc=="Full9435_1d"){med9435.open("med9435_1D.txt");}
        //else if(dati_o_mc=="Full9435_2d"){med9435.open("med9435_2D.txt");}

        if(dati_o_mc=="data_1d_run2_full_v09_89"){med9435.open("Run2_1D_v09_89_full.txt");}
        else if(dati_o_mc=="dati2d"){med9435.open("Run2_2D_full.txt");}

        TH1D *h_median = new TH1D(Form("hmedian_%s",dati_o_mc.c_str()), "", 100, 0, 20 );

        EventsData dat = load_data(dati_o_mc , "muon");
    
        cout << dat.tree->GetEntries() << " total tracks" << endl;

        std::vector<double> trackdEdx; 

        for(int track=0; track<dat.tree->GetEntries(); track++)
        {
            dat.tree->GetEntry(track);

            for(int hit=0; hit<int(dat.track.rr->size()); hit++)
            {
                if(dat.track.rr->at(hit)<5.)
                {
                    trackdEdx.push_back(dat.track.dE->at(hit));
                }
            }

            if (trackdEdx.size() >= 2) {trackdEdx.erase(trackdEdx.end() - 2, trackdEdx.end());} //removing the last 2 hits to compute the median
    
            double med =0;
            if(trackdEdx.size()>0){
                med= mediana(trackdEdx); 
                h_median->Fill(med);
            }
            else med=-1;

            trackdEdx.clear();

            med9435 << dat.EvtInfo.run << " " << dat.EvtInfo.evt << " " <<  dat.track.len_reco << " " << dat.track.end_reco->at(0) << " " << dat.track.end_reco->at(1) << " " << dat.track.end_reco->at(2) << " " << med << endl;
        
        }
    med9435.close();

    f->cd();
    h_median->Scale(1./h_median->Integral());
    h_median->Write(0,TObject::kOverwrite);
    }
    

    
}

bool is_present(std::vector<int> v, int valore)
{
    int N=int(v.size());
    for(int i=0; i<N; i++)
    {
        if(v[i]==valore)return true;
    }
    return false;
}

ofstream dump_mpv_1d2d;

void confrontoDatiMClight(std::string particle, std::string dati_o_mc, const char *filename)
{

    //if(dati_o_mc=="dati"){dump_mpv_1d2d.open("dump_mpv_1d_mult1.txt");}
    //else if(dati_o_mc=="dati2d10"){dump_mpv_1d2d.open("dump_mpv_2d_mult1.txt");}

    EventsData dat = load_data(dati_o_mc , particle, "Np", "full");
    
    cout << dat.tree->GetEntries() << " total tracks" << endl;

    //std::string data_dir = dati_o_mc + "_mult1";

    TFile *f = TFile::Open(filename , "UPDATE");
    TDirectory *d= (TDirectory*)f->Get(particle.c_str());
    TDirectory *d1= (TDirectory*)d->Get(dati_o_mc.c_str()); 
    //TDirectory *d1= (TDirectory*)d->Get(data_dir.c_str()); 

    d1->cd();

    //TFile * f = new TFile("hmedian_new.root", "RECREATE");

    std::vector<std::vector<double>> dE;
    std::vector<std::vector<double>> rr;
    std::vector<double> track_length;
    

    int Nbins_rr=100;
    double low_rr=0;
    double high_rr=30;
    

    TH2D *dedx_range= new TH2D("dedx_range", " ", Nbins_rr, low_rr, high_rr, 300, 0, 30);

    TH2D *dEdx_vs_median = new TH2D("dEdx_vs_median", "", 300,0,30, 200,0,20);

    TH1D *h_median = new TH1D("h_median", "", 200, 0, 20 );

    //TH1D *dEdx_range_media_reference = new TH1D("dEdx_range_media_reference", "", Nbins_rr, low_rr, high_rr);

    TH1D *hdirx = new TH1D("dirx", "", 200, -1, 1);
    TH1D *hdiry = new TH1D("diry", "", 200, -1, 1);
    TH1D *hdirz = new TH1D("dirz", "", 200, -1, 1);

    TH1D *h_theta_drift = new TH1D("htheta_drif", "", 180, 0, 90 );

    //TH1D *h_theta_drift_cut = new TH1D("htheta_drif_cut", "", 180, 0, 90 );

    //dEdx TRA 25 E 30 CM
    //TH1D *hdEdx = new TH1D("hdEdx", "", 200, 0., 10.);
    //TH1D *hdEdx = new TH1D("hdEdx", "", 80, 0., 6.);
    TH1D *hdEdx = new TH1D("hdEdx", "", 200, 0., 15.);
    TH1D *hdEdx_onlymult1 = new TH1D("hdEdx_onlymult1", "", 200, 0, 15);
    TH1D *hdEdx_multmag1 = new TH1D("hdEdx_multmag1", "", 200, 0, 15);
    TH1D *hdQdx = new TH1D("hdQdx", "", 500, 0., 2000.);
    TH1D *hdQdx_onlymult1 = new TH1D("hdQdx_onlymult1", "", 500, 0., 2000.);
    TH1D *hdQdx_multmag1 = new TH1D("hdQdx_multmag1", "", 500, 0., 2000.);
    TH1D *hdEdx_ind1 = new TH1D("hdEdx_ind1", "", 200, 0., 15.);
    TH1D *hdEdx_ind2 = new TH1D("hdEdx_ind2", "", 200, 0., 15.);
    TH2D *dEdx_range_ind1 = new TH2D("dEdx_range_ind1","", 100,0,30,300,0,30);
    TH2D *dEdx_range_ind2 = new TH2D("dEdx_range_ind2", "", 100,0,30,300,0,30);

    std::vector<double> trackdEdx; //the one in the last 5 cm that will be used to compute the median
    std::vector<double> rr_temp;
    std::vector<double> dEdx_temp;

    for(int track=0; track<dat.tree->GetEntries(); track++)
    {
        dat.tree->GetEntry(track);

        hdirx->Fill(dat.track.dirx);
        hdiry->Fill(dat.track.diry);
        hdirz->Fill(dat.track.dirz);

        double theta_drift = 180./M_PI*acos(abs(dat.track.dirx));
        h_theta_drift->Fill(theta_drift);

        //std::vector<double> dQdx_splitted_hit;
        //std::vector<double> pitch_splitted_hit;
        //std::vector<double> mult_splitted_hit;
        //std::vector<double> sumadc_splitted_hit;

        if(dat.track.len_reco>=30)
        {
            for(int hit=0; hit<int(dat.track.rr->size()); hit++)
            {
                if(dat.track.rr->at(hit)>=25 && dat.track.rr->at(hit)<=30)
                {
                    hdEdx->Fill(dat.track.dE->at(hit));
                    hdQdx->Fill(dat.track.dQdx->at(hit));

                    if(dat.track.mult->at(hit)==1)
                    {  
                        hdEdx_onlymult1->Fill(dat.track.dE->at(hit));
                        hdQdx_onlymult1->Fill(dat.track.dQdx->at(hit));
                    }
                    else if(dat.track.mult->at(hit)!=1)
                    {
                        hdEdx_multmag1->Fill(dat.track.dE->at(hit));
                        hdQdx_multmag1->Fill(dat.track.dQdx->at(hit));
                    }

                    //if(dat.track.mult->at(hit)>1)
                    //{  
                        //cout << dat.track.dQdx->at(hit) << " " << dat.track.pitch->at(hit) << " " << dat.track.mult->at(hit) << " " << dat.track.sumadc->at(hit) << endl;
                        //dQdx_splitted_hit.push_back(dat.track.dQdx->at(hit));
                        //pitch_splitted_hit.push_back(dat.track.pitch->at(hit));
                        //mult_splitted_hit.push_back(dat.track.mult->at(hit));
                        //sumadc_splitted_hit.push_back(dat.track.sumadc->at(hit));
                    //}
                }
            }
        }

        /*
        std::vector<int> processed_indeces;
        for(int j=0; j<int(dQdx_splitted_hit.size()); j++)
        {
            //cout << dQdx_splitted_hit[j] << " " << pitch_splitted_hit[j] << " " << sumadc_splitted_hit[j] << " " << mult_splitted_hit[j] << endl;

            processed_indeces.push_back(j);
            for(int i=0; i<int(dQdx_splitted_hit.size()); i++)
            {
                if(!is_present(processed_indeces,i) && sumadc_splitted_hit[i]==sumadc_splitted_hit[j])
                {
                    //cout << "amica trovata" << endl;
                    processed_indeces.push_back(i);
                    //cout << dQdx_splitted_hit[i] << " " << pitch_splitted_hit[i] << " " << sumadc_splitted_hit[i] << " " << mult_splitted_hit[i] << endl;
                }
            }
            //cout << endl; 
        }
        */
        

        //looking also to ind1 and 1nd2
        if(dati_o_mc!="dati" && dati_o_mc!="dati10")
        {
            if(dat.track.len_reco>=30)
            {
                //cout << "nuova traccia ind1" << endl;
                for(int hit=0; hit<int(dat.track.rr_ind1->size()); hit++)
                {
                    if(dat.track.rr_ind1->at(hit)>=25. && dat.track.rr_ind1->at(hit)<30.)
                    {
                        hdEdx_ind1->Fill(dat.track.dEdx_ind1->at(hit));
                    }
                }

                //cout << "nuova traccia ind2" << endl;
                for(int hit=0; hit<int(dat.track.rr_ind2->size()); hit++)
                {
                    if(dat.track.rr_ind2->at(hit)>=25. && dat.track.rr_ind2->at(hit)<30.)
                    {
                        hdEdx_ind2->Fill(dat.track.dEdx_ind2->at(hit));
                    }
                }
            }

            //plot banana
            for(int hit=0; hit<int(dat.track.rr_ind1->size()); hit++){dEdx_range_ind1->Fill(dat.track.rr_ind1->at(hit), dat.track.dEdx_ind1->at(hit));}
            for(int hit=0; hit<int(dat.track.rr_ind2->size()); hit++){dEdx_range_ind2->Fill(dat.track.rr_ind2->at(hit), dat.track.dEdx_ind2->at(hit));}
        }


        // DA QUI IN POI SI FANNO I FIT
        for(int hit=0; hit<int(dat.track.rr->size()); hit++)
        {
            //if(dat.track.mult->at(hit)>1)continue;

            dedx_range->Fill(dat.track.rr->at(hit), dat.track.dE->at(hit));

            rr_temp.push_back(dat.track.rr->at(hit));
            dEdx_temp.push_back(dat.track.dE->at(hit));

            if(dat.track.rr->at(hit)<5.)
            {
                trackdEdx.push_back(dat.track.dE->at(hit));
            }
        }

        if (trackdEdx.size() >= 2) {trackdEdx.erase(trackdEdx.end() - 2, trackdEdx.end());} //removing the last 2 hits to compute the median
    
        double med =0;
        if(trackdEdx.size()>0) 
        {
            med= mediana(trackdEdx);
            h_median->Fill(med);
            for(int i=0; i<int(dEdx_temp.size()); i++)
            {
                dEdx_vs_median->Fill(dEdx_temp[i], med);
            }
        }

        double threshold;
        if(particle=="muon" && (dati_o_mc=="dati" || dati_o_mc=="dati10" || dati_o_mc=="mc")){threshold=3.2;}
        else if(particle=="muon" && (dati_o_mc=="dati2d" || dati_o_mc=="dati2d10")){threshold=3.5;}
        else if(particle=="proton"){threshold=8.;}

        if(med>=threshold)
        {
            dE.push_back(dEdx_temp);
            rr.push_back(rr_temp);
            track_length.push_back(dat.track.len_reco);
        }


        trackdEdx.clear();
        rr_temp.clear();
        dEdx_temp.clear();
    }

    hdEdx->Scale(1./hdEdx->Integral());
    hdEdx->Write(0,TObject::kOverwrite);
    h_median->Scale(1./h_median->Integral());
    h_median->Write(0,TObject::kOverwrite);
    dedx_range->Write(0,TObject::kOverwrite);
    dEdx_vs_median->Write(0,TObject::kOverwrite);
    hdirx->Scale(1./hdirx->Integral());
    hdiry->Scale(1./hdiry->Integral());
    hdirz->Scale(1./hdirz->Integral());
    h_theta_drift->Scale(1./h_theta_drift->Integral());
    hdirx->Write(0,TObject::kOverwrite);
    hdiry->Write(0,TObject::kOverwrite);
    hdirz->Write(0,TObject::kOverwrite);
    h_theta_drift->Write(0,TObject::kOverwrite);
    hdQdx -> Scale(1./hdQdx->Integral());
    hdEdx_ind1 -> Scale(1./hdEdx_ind1->Integral());
    hdEdx_ind2 -> Scale(1./hdEdx_ind2->Integral());
    hdQdx->Write(0,TObject::kOverwrite);
    hdEdx_ind1->Write(0,TObject::kOverwrite);
    hdEdx_ind2->Write(0,TObject::kOverwrite);
    dEdx_range_ind1->Write(0,TObject::kOverwrite);
    dEdx_range_ind2->Write(0,TObject::kOverwrite);

    hdEdx_onlymult1->Scale(1./hdEdx_onlymult1->Integral());
    hdQdx_onlymult1->Scale(1./hdQdx_onlymult1->Integral());
    hdEdx_multmag1->Scale(1./hdEdx_multmag1->Integral());
    hdQdx_multmag1->Scale(1./hdQdx_multmag1->Integral());

    hdEdx_onlymult1->Write(0,TObject::kOverwrite);
    hdQdx_onlymult1->Write(0,TObject::kOverwrite);
    hdEdx_multmag1->Write(0,TObject::kOverwrite);
    hdQdx_multmag1->Write(0,TObject::kOverwrite);

    cout << particle << " " << dati_o_mc << " " << rr.size() << " tracks passed" << endl;

    /*
    TGraphErrors * graph = new TGraphErrors();

    //loop su tutti i bin

    std::array<double,9> low_extremes={15., 14., 12.5, 11.5, 10.5, 10., 9.5, 9.1, 8.5 };

    for(int bin=3; bin<=dedx_range->GetNbinsX(); bin++) //cutting all hits with a RR<0.6 cm, starting from bin center RR=0.75cm 
    {

        int NBIN=600;

        double x_low;
        double x_high;

        if(particle=="muon")
        {
            if(bin<30)
            {
                NBIN=150;
            }
        }

        if(particle=="proton")
        {
            if(bin<30)
            {
                NBIN=200;
            }
            if(bin==3){NBIN=100;}
        }


        double dEdx_high_limit;
        if(particle=="muon"){dEdx_high_limit=30.;}
        else if(particle=="proton"){dEdx_high_limit=40.;}
        TH1D* dEdx_i = new TH1D(Form("dEdx_i_%d", bin), "", NBIN, 0., dEdx_high_limit);


        std::vector<double> dEdx_fit;

        for(int trk=0; trk<int(rr.size()); trk++ )
        {
            for(int hit =0; hit<int(rr[trk].size()); hit++)
            {
                if(rr[trk][hit]<bin*0.3 && rr[trk][hit]>=(bin-1)*0.3)
                {
                    dEdx_i->Fill(dE[trk][hit]);
                    dEdx_fit.push_back(dE[trk][hit]);
                }
            }
        
        }

        int binMax=0;
        double maxVal=0;
        for(int j=1; j<dEdx_i->GetNbinsX(); j++)
        {
            double y=dEdx_i->GetBinContent(j);
            if(y>maxVal)
            {
                maxVal=y;
                binMax=j;
            }
        }

        if(particle=="proton")
        {
            x_low=dEdx_i->GetBinCenter(binMax)-0.60*dEdx_i->GetStdDev();
            x_high=dEdx_i->GetBinCenter(binMax)+3*dEdx_i->GetStdDev();
            if(bin<12)
            {
                x_low=low_extremes[bin-3];
            }
        }
        if(particle=="muon")
        {
            x_low=dEdx_i->GetBinCenter(binMax)-0.50*dEdx_i->GetStdDev();
            x_high=dEdx_i->GetBinCenter(binMax)+3*dEdx_i->GetStdDev();
            if(bin<100){x_low=dEdx_i->GetBinCenter(binMax)-0.75*dEdx_i->GetStdDev();}
            if(bin==3){x_low=dEdx_i->GetBinCenter(binMax)-1.2*dEdx_i->GetStdDev();}
        }

        std::string model = "langau";

        FitOutput fitOut = fitLangaulight_unbinned(dEdx_i, dEdx_fit, x_low, x_high , d1, model, bin);
        //FitOutput fitOut = fitLangaulight_unbinned(dEdx_i, dEdx_fit, 0.5, 20., d1, model, bin);

        cout << bin << " " << "mpv: " << fitOut.fitResults[0] << " ± " << fitOut.fitResults[1] << " Gauss sigma: " << fitOut.fitResults[2] << " ± " << fitOut.fitResults[3];
        cout << " Landau sigma: " << fitOut.fitResults[4] << " ± " << fitOut.fitResults[5] << " minNll: " << fitOut.fitResults[7] << " status: " << fitOut.fitResults[8];
        cout << " covQual: " << fitOut.fitResults[9] << endl;

        dump_mpv_1d2d << dedx_range->GetXaxis()->GetBinCenter(bin) << " " << fitOut.fitResults[0] << " " << fitOut.fitResults[1] << " " << fitOut.fitResults[2] << " " << fitOut.fitResults[3] << " " << fitOut.fitResults[4] << " " << fitOut.fitResults[5] << " " << fitOut.fitResults[7] << " " << fitOut.fitResults[8] << " " << fitOut.fitResults[9] << endl;

        //TF1 *LandauGauss = fitOut.fitFunction;

        TF1 * fitfun;
        fitfun = fitOut.fitFunction;
        double norm=fitfun->Integral(fitfun->GetXmin(), fitfun->GetXmax(),1e-3);  
        TF1* langau_pdf = new TF1(
            Form("scaled_density_rr%d", bin),
            [fitfun, norm](double* x, double*) { return 1./norm * fitfun->Eval(x[0]); },
            fitfun->GetXmin(),
            fitfun->GetXmax(),
            0 // nessun parametro
        );
        langau_pdf->SetNpx(1000);
        langau_pdf->Write(Form("density_rr%.d", bin), TObject::kOverwrite);

        //cout << bin << " " << fitOut.fitResults[0] << " " << LandauGauss->Integral(dEdx_i->GetXaxis()->GetXmin(), dEdx_i->GetXaxis()->GetXmax() ) << endl;
        dEdx_i->Scale(1./dEdx_i->Integral("width"));
        dEdx_i->Write(0,TObject::kOverwrite);

        int n = graph->GetN();
        graph->SetPoint(n, dedx_range->GetXaxis()->GetBinCenter(bin), fitOut.fitResults[0]);
        graph->SetPointError(n, (dedx_range->GetXaxis()->GetBinWidth(bin)) / 2, fitOut.fitResults[1]);

    }
    */

    
    //graph->Write(0,TObject::kOverwrite);

    

    f->Close();
    delete f;
    
}




void confrontoEintDatiMC(std::string dati_o_mc, std::string particle)
{
    EventsData dat = load_data(dati_o_mc , particle, "Np", "25only");
    
    cout << dat.tree->GetEntries() << " total tracks" << endl;

    TH1D *Eint = new TH1D(Form("Eint_%s_%s", particle.c_str(), dati_o_mc.c_str() ), "", 40, 0., 120. );
    TH1D *Mediana_Eint = new TH1D(Form("Mediana_Eint_%s_%s", particle.c_str(), dati_o_mc.c_str() ), "", 40, 0., 120. );

    for(int track=0; track<dat.tree->GetEntries(); track++)
    {
        dat.tree->GetEntry(track);
        std::vector<double> trackdEdx;

        //if(distanza(dat.track.end_reco, dat.track.end_true)>=1.){continue;}
        for(int hit=0; hit<int(dat.track.rr->size()); hit++)
        { 
            if(dat.track.rr->at(hit)<5.)
                {
                    trackdEdx.push_back(dat.track.dE->at(hit));
                }
                
            if(dat.track.rr->at(hit)<12.7 && dat.track.rr->at(hit)>12.3)
            {
                Eint->Fill(dat.track.Eint->at(hit));

            }
        }

        double med =0;
        if(trackdEdx.size()>0) 
        {
            cout << "ok" << endl;
            med= mediana(trackdEdx);
        }

        if(med>3.2)
        {
            cout << "ok" << endl;
            for(int hit=0; hit<int(dat.track.rr->size()); hit++)
            {
                if(dat.track.rr->at(hit)<12.7 && dat.track.rr->at(hit)>12.3)
                {
                    Mediana_Eint->Fill(dat.track.Eint->at(hit));
            
                }    
            }
        }
    }

    TFile *f = TFile::Open("confrontoEint.root", "UPDATE");
    Eint->Scale(1./Eint->Integral());
    Eint->Write(0,TObject::kOverwrite);
    Mediana_Eint->Scale(1./Mediana_Eint->Integral());
    Mediana_Eint->Write(0, TObject::kOverwrite);

    f->Close();
    delete f;
    
}

Double_t correction_function(Double_t *x, Double_t *par)
{
    if(x[0]<=1.64000) return 0.0354575;
    else if(x[0]>=15.0000) return -0.0164076;
    else return exp(par[0]+par[1]*x[0])+par[2]+par[3]*exp(-0.5*pow((x[0]-par[4])/par[5],2));
}

Double_t correction_function_25only(Double_t *x, Double_t *par)
{
    if(x[0]<=2) return 0.0127159;
    else if (x[0]>=15) return -0.0235248;
    else return exp(par[0]+par[1]*x[0])+par[2]+par[3]*exp(-0.5*pow((x[0]-par[4])/par[5],2));
}

Double_t correction_muon(Double_t *x, Double_t *par)
{
    if(x[0]<=2.) return 0.0120827;
    else if(x[0]>=5.) return 0.0373432;
    else return par[0]*(x[0]-par[1])*(x[0]-par[1])+par[2]*x[0]+par[3];
}

Double_t correction_proton(Double_t *x, Double_t *par)
{
    if(x[0]<=4.2) return 0.0125211;
    else if(x[0]>=15.) return  -0.0288104;
    else return par[0]+par[1]*TMath::Log(1+par[2]*(x[0]-par[3])*(x[0]-par[3]));
}


void confrontoDatiMC(std::string particle, std::string dati_o_mc, std::string corrected, const char *filename, std::string range)
{

      gROOT->ProcessLine(".L libdaughtersInfo.so");
    std::ofstream check_factor_mu;
    std::ofstream check_factor_pro;
    if(particle=="muon" && dati_o_mc=="mc" && corrected=="corrected")check_factor_mu.open("check_factor_mu.txt");
    if(particle=="proton" && dati_o_mc=="mc" && corrected=="corrected")check_factor_pro.open("check_factor_pro.txt");

    TF1 *fitMCcorr;
    if(range==""){
        fitMCcorr= new TF1("fitMCcorr", correction_function, 0., 16., 6);
        fitMCcorr->SetParameters(3.13069e-01, -2.08971e+00, -1.81601e-02, 4.87814e-02, 7.17860e+00, -3.03200e+00);
    } //full
    else if(range=="25only"){
        fitMCcorr= new TF1("fitMCcorr", correction_function_25only, 0., 16., 6);
        fitMCcorr->SetParameters(-1.08107e+00, -1.38302e+00, -2.84635e-02, 5.87297e-02, 7.17814e+00, -3.51461e+00);
    } //25only
    else if(range=="specific")
    {
        if(particle=="muon")
        {
            fitMCcorr= new TF1("fitMCcorr", correction_muon, 0., 16. , 4);
            fitMCcorr->SetParameters(7.65506e-03, 1.59985e+00, -2.06713e-02, 5.21996e-02);
        }//only muon and 25 cm rr
        if(particle=="proton")
        {
            fitMCcorr= new TF1("fitMCcorr", correction_proton, 0., 16., 4);
            fitMCcorr->SetParameters(3.04412e-02, -3.33271e-02, 8.03482e-02, 7.17696e+00);
        }//only proton and 25 cm rr
    }

    EventsData dat = load_data(dati_o_mc , particle, "Np", "full");
    
    cout << dat.tree->GetEntries() << " total tracks" << endl;

    TFile *f = TFile::Open(filename , "UPDATE");
    TDirectory *d= (TDirectory*)f->Get(particle.c_str());
    TDirectory *d1= (TDirectory*)d->Get(dati_o_mc.c_str()); 
    d1->cd();

    TGraph *graph_func = new TGraph();
    for(int i=0; i<16000; i++)
    {
        graph_func->SetPoint(i, i/1000., fitMCcorr->Eval(i/1000.));
    }
    graph_func->SetMarkerStyle(7);
    graph_func->Write("correction_function_graph", TObject::kOverwrite);

    std::vector<std::vector<double>> dE;
    std::vector<std::vector<double>> rr;
    std::vector<double> track_length;

    
    std::vector<std::vector<double>> dE_tot;
    std::vector<std::vector<double>> rr_tot;
    std::vector<double> track_length_tot;
    

    int Nbins_rr;
    double low_rr, high_rr;
    if(particle=="muon" && range==""){Nbins_rr=500; low_rr=0.; high_rr=150; }
    else if(particle=="proton" && range==""){Nbins_rr=200; low_rr=0; high_rr=60;}
    else if(range=="25only" || range=="specific"){Nbins_rr=85; low_rr=0.; high_rr=25.5; }
    

    TH2D *dedx_range= new TH2D("dedx_range", " ", Nbins_rr, low_rr, high_rr, 300, 0, 30);

    TH2D *dedx_range_25only = new TH2D("dedx_range_25only", " ", 500, 0, 25, 300, 0, 30);

    TH1D *hit_number_last5cm = new TH1D("hit_number_last5cm", " ", 100, 0, 100);

    TH2D *dEdx_vs_median = new TH2D("dEdx_vs_median", "", 300,0,30, 200,0,20);

    TH1D *h_median = new TH1D("h_median", "", 200, 0, 20 );

    TH1D *scarto = new TH1D("scarto", "", 200, -10, 10);

    TH1D *dEdx_range_media_reference = new TH1D("dEdx_range_media_reference", "", Nbins_rr, low_rr, high_rr);

    TH1D *hdirx = new TH1D("dirx", "", 200, -1, 1);
    TH1D *hdiry = new TH1D("diry", "", 200, -1, 1);
    TH1D *hdirz = new TH1D("dirz", "", 200, -1, 1);

    TH1D *h_theta_drift = new TH1D("htheta_drif", "", 180, 0, 90 );

    TH1D *h_theta_drift_cut = new TH1D("htheta_drif_cut", "", 180, 0, 90 );

    for(int track=0; track<dat.tree->GetEntries(); track++)
    {
        dat.tree->GetEntry(track);

        hdirx->Fill(dat.track.dirx);
        hdiry->Fill(dat.track.diry);
        hdirz->Fill(dat.track.dirz);

        double theta_drift = 180./M_PI*acos(abs(dat.track.dirx));
        h_theta_drift->Fill(theta_drift);
        //if(theta_drift>30)continue;

        h_theta_drift_cut->Fill(theta_drift);

        std::vector<double> trackdEdx;
        std::vector<double> rr_temp;
        std::vector<double> dEdx_temp;

        //cout << dat.track.rr->size() << endl;

        for(int hit=0; hit<int(dat.track.rr->size())-2; hit++)//removing the first two hits
        { 
            double factor=1;
            if(dati_o_mc=="mc" && corrected=="corrected"){factor = (2+fitMCcorr->Eval(dat.track.dE->at(hit)))/(2-fitMCcorr->Eval(dat.track.dE->at(hit)));}
            dedx_range->Fill(dat.track.rr->at(hit), dat.track.dE->at(hit)*factor);
            if(particle=="proton" && dati_o_mc=="mc" && corrected=="corrected") {check_factor_pro << dat.track.dE->at(hit) << " " << factor << endl;}
            if(particle=="muon" && dati_o_mc=="mc" && corrected=="corrected"){check_factor_mu << dat.track.dE->at(hit) << " " << factor << endl;}
            dedx_range_25only->Fill(dat.track.rr->at(hit), dat.track.dE->at(hit)*factor);
            rr_temp.push_back(dat.track.rr->at(hit));
            dEdx_temp.push_back(dat.track.dE->at(hit)*factor);

            //cout << dat.track.rr->at(hit) << " " ;

            if(dat.track.rr->at(hit)<5.)
            {
                trackdEdx.push_back(dat.track.dE->at(hit)*factor);
            }
        }

        //senza la condizione sulla mediana
        
        dE_tot.push_back(dEdx_temp);
        rr_tot.push_back(rr_temp);
        track_length_tot.push_back(dat.track.len_reco);
        

        hit_number_last5cm->Fill(trackdEdx.size());


        double med =0;
        if(trackdEdx.size()>0) 
        {
            med= mediana(trackdEdx);
            h_median->Fill(med);
            for(int i=0; i<int(dEdx_temp.size()); i++)
            {
                dEdx_vs_median->Fill(dEdx_temp[i], med);
            }
        }

        double threshold;
        if(particle=="muon"){threshold=3.2;}
        else if(particle=="proton"){threshold=8.;}

        if(med>=threshold)
        {
            dE.push_back(dEdx_temp);
            rr.push_back(rr_temp);
            track_length.push_back(dat.track.len_reco);
        }


        trackdEdx.clear();
        rr_temp.clear();
        dEdx_temp.clear();
    }

    h_median->Scale(1./h_median->Integral());
    h_median->Write(0,TObject::kOverwrite);
    hit_number_last5cm->Write(0,TObject::kOverwrite);
    dedx_range->Write(0,TObject::kOverwrite);
    dedx_range_25only->Write(0,TObject::kOverwrite);
    dEdx_vs_median->Write(0,TObject::kOverwrite);
    hdirx->Scale(1./hdirx->Integral());
    hdiry->Scale(1./hdiry->Integral());
    hdirz->Scale(1./hdirz->Integral());
    h_theta_drift->Scale(1./h_theta_drift->Integral());
    h_theta_drift_cut->Scale(1./h_theta_drift_cut->Integral());
    

    hdirx->Write(0,TObject::kOverwrite);
    hdiry->Write(0,TObject::kOverwrite);
    hdirz->Write(0,TObject::kOverwrite);
    h_theta_drift->Write(0,TObject::kOverwrite);
    h_theta_drift_cut->Write(0,TObject::kOverwrite);

    cout << particle << " " << dati_o_mc << " " << rr.size() << " tracks passed" << endl;

    TGraphErrors * graph = new TGraphErrors();

    std::string range_scritto;
    if(range=="")range_scritto=range;
    else range_scritto=(range+"_");
    //ofstream confrontoDatiMCdump(Form("/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/MPV/%sconfronto_%s_%s%s.txt", range_scritto.c_str(), particle.c_str(), dati_o_mc.c_str(), ("_"+corrected).c_str() ));
    //ofstream confrontoDatiMCdump("text.txt");
    TH2D *dEdx_range_median_mag4 = new TH2D("dEdx_range_median_mag4", "", 500, 0, 25, 300, 0, 30 );

    TH1D *dEdx_100_150 = new TH1D("dEdx_100_150", "", 200, 0., 10.);
    TH1D *dEdx_100_150_tot = new TH1D("dEdx_100_150_tot", "", 200, 0., 10.);

    //loop su tutti i bin
    //if(dati_o_mc=="dati"){MPVsfitDATI.clear();}
    //if(dati_o_mc=="mc"){MPVsfitMC.clear();}

    for(int trk=0; trk<int(rr.size()); trk++ )
        {
            for(int hit =0; hit<int(rr[trk].size()); hit++)
            {
                dEdx_range_median_mag4->Fill(rr[trk][hit],dE[trk][hit]);
            }
        
        }

    std::array<double,9> low_extremes={15., 14., 12.5, 11.5, 10.5, 10., 9.5, 9.1, 8.5 };

    for(int bin=3; bin<=dedx_range->GetNbinsX(); bin++) //cutting all hits with a RR<0.6 cm, starting from bin center RR=0.75cm 
    {

        int NBIN=600;

        double x_low;
        double x_high;

        if(dati_o_mc=="mc" && particle=="muon")
        {
            if(bin<30)
            {
                NBIN=150;
            }
        }
        else if(dati_o_mc=="dati" && particle=="muon")
        {
            if(bin<30)
            {
                NBIN=150;
            }
        }

        if(dati_o_mc=="mc" && particle=="proton")
        {
            if(bin<30)
            {
                NBIN=200;
            }
            if(bin==3){NBIN=100;}
        }
        else if(dati_o_mc=="dati" && particle=="proton")
        {
            if(bin<30)
            {
                NBIN=200;
            }
            if(bin==3){NBIN=100;}
        }


        double dEdx_high_limit;
        if(particle=="muon"){dEdx_high_limit=30.;}
        else if(particle=="proton"){dEdx_high_limit=40.;}
        TH1D* dEdx_i = new TH1D(Form("dEdx_i_%d", bin), "", NBIN, 0., dEdx_high_limit);


        std::vector<double> dEdx_fit;

        for(int trk=0; trk<int(rr.size()); trk++ )
        {
            for(int hit =0; hit<int(rr[trk].size()); hit++)
            {
                if(rr[trk][hit]<bin*0.3 && rr[trk][hit]>=(bin-1)*0.3)
                {
                    dEdx_i->Fill(dE[trk][hit]);
                    //dEdx_range_median_mag4->Fill(rr[trk][hit],dE[trk][hit]);
                    dEdx_fit.push_back(dE[trk][hit]);
                }
            }
        
        }

        dEdx_range_media_reference->SetBinContent(bin, dEdx_i->GetMean());
        dEdx_range_media_reference->SetBinError(bin, dEdx_i->GetStdDev());

        if(particle=="proton")
        {
            int binMax=0;
            double maxVal=0;
            for(int j=1; j<dEdx_i->GetNbinsX(); j++)
            {
                double y=dEdx_i->GetBinContent(j);
                if(y>maxVal)
                {
                    maxVal=y;
                    binMax=j;
                }
            }

            x_low=dEdx_i->GetBinCenter(binMax)-0.60*dEdx_i->GetStdDev();
            x_high=dEdx_i->GetBinCenter(binMax)+3*dEdx_i->GetStdDev();
            if(bin<12)
            {
                x_low=low_extremes[bin-3];
            }
        }
        if(particle=="muon")
        {
            int binMax=0;
            double maxVal=0;
            for(int j=1; j<dEdx_i->GetNbinsX(); j++)
            {
                double y=dEdx_i->GetBinContent(j);
                if(y>maxVal)
                {
                    maxVal=y;
                    binMax=j;
                }
            }
            x_low=dEdx_i->GetBinCenter(binMax)-0.50*dEdx_i->GetStdDev();
            x_high=dEdx_i->GetBinCenter(binMax)+3*dEdx_i->GetStdDev();
            if(bin<100){x_low=dEdx_i->GetBinCenter(binMax)-0.75*dEdx_i->GetStdDev();}
            if(bin==3){x_low=dEdx_i->GetBinCenter(binMax)-1.2*dEdx_i->GetStdDev();}
        }

        std::string model = "langau";


        FitOutput fitOut = fitLangaulight_unbinned(dEdx_i, dEdx_fit, x_low, x_high , d1, model, bin);

        //TF1 *LandauGauss = fitOut.fitFunction;

        TF1 * fitfun;
        fitfun = fitOut.fitFunction;
        double norm=fitfun->Integral(fitfun->GetXmin(), fitfun->GetXmax(),1e-3);  
        TF1* langau_pdf = new TF1(
            Form("scaled_density_rr%d", bin),
            [fitfun, norm](double* x, double*) { return 1./norm * fitfun->Eval(x[0]); },
            fitfun->GetXmin(),
            fitfun->GetXmax(),
            0 // nessun parametro
        );
        langau_pdf->SetNpx(1000);
        langau_pdf->Write(Form("density_rr%.d", bin), TObject::kOverwrite);

        //cout << bin << " " << fitOut.fitResults[0] << " " << LandauGauss->Integral(dEdx_i->GetXaxis()->GetXmin(), dEdx_i->GetXaxis()->GetXmax() ) << endl;
        dEdx_i->Scale(1./dEdx_i->Integral("width"));
        dEdx_i->Write(0,TObject::kOverwrite);

        /*
        TCanvas *canvas = new TCanvas(Form("dEdx %s %s rr [%f , %f] cm", particle.c_str(), dati_o_mc.c_str(), (bin-1)*0.3 , bin*0.3));
        dEdx_i->Draw("p");
        dEdx_i->SetMarkerStyle(7);
        dEdx_i->SetLineWidth(2);
        LandauGauss->SetLineColor(kRed);
        LandauGauss->Draw("same");
        std::vector<double> vec_larghezza = FWHM(LandauGauss,LandauGauss->GetMaximumX(x_low, x_high));
        TLine* l1 = new TLine(vec_larghezza[0], 0, vec_larghezza[0], vec_larghezza[3]);
        TLine* l2 = new TLine(vec_larghezza[1], 0, vec_larghezza[1], vec_larghezza[3]);
        TLine* l3 = new TLine(vec_larghezza[0], vec_larghezza[3], vec_larghezza[1], vec_larghezza[3]);
        l1->SetLineStyle(kDashed);
        l2->SetLineStyle(kDashed);
        l3->SetLineStyle(kDashed);
        l1->SetLineColor(kOrange);
        l2->SetLineColor(kOrange);
        l3->SetLineColor(kOrange);
        l1->Draw("same");
        l2->Draw("same");
        l3->Draw("same");
        //int bin_mpv = dEdx_i->FindBin(fitOut.fitResults[0]);
        TArrow *a = new TArrow(fitOut.fitResults[0], 0, fitOut.fitResults[0], LandauGauss->Eval(fitOut.fitResults[0]));
        a->SetLineStyle(kDashed);
        a->SetLineColor(kGreen);
        a->Draw("same");
        //int bin_mpv_start = dEdx_i->FindBin(fitOut.fitResults[6]);
        //TArrow *a1 = new TArrow(fitOut.fitResults[6], 0, fitOut.fitResults[6], dEdx_i->GetBinContent(bin_mpv_start));
        //a1->SetLineStyle(kDashed);
        //a1->SetLineColor(kGreen);
        //a1->Draw("same");
        canvas->Write(0,TObject::kOverwrite);
        */

        //confrontoDatiMCdump << dedx_range->GetXaxis()->GetBinCenter(bin) << " " <<  fitResults[0] << " " << fitResults[1] << " " << fitResults[7] << " " << fitResults[8] << " " << fitResults[9] << endl; 
        //confrontoDatiMCdump << dedx_range->GetXaxis()->GetBinCenter(bin) << " " <<  fitResults[6] << " " << dEdx_i->GetStdDev()  << endl; 
        //confrontoDatiMCdump << dedx_range->GetXaxis()->GetBinCenter(bin) << " " <<  fitOut.fitResults[0] << " " << fitOut.fitResults[1] << " " << fitOut.fitResults[2] << " " << fitOut.fitResults[3] << " " << fitOut.fitResults[4] << " " << fitOut.fitResults[5] << " " << dEdx_i->GetMean() << " "<< dEdx_i->GetStdDev() << " " <<  fitOut.fitResults[10] << endl;

        //if(dati_o_mc=="dati"){MPVsfitDATI.push_back(fitResults[0]);}
        //if(dati_o_mc=="mc"){MPVsfitMC.push_back(fitResults[0]);}

        scarto->Fill(fitOut.fitResults[0]-fitOut.fitResults[6]);

    int n = graph->GetN();
    graph->SetPoint(n, dedx_range->GetXaxis()->GetBinCenter(bin), fitOut.fitResults[0]);
    graph->SetPointError(n, (dedx_range->GetXaxis()->GetBinWidth(bin)) / 2, fitOut.fitResults[1]);

    //delete dEdx_i;
    //delete fitMCcorr;
    }

    
    //controllo - dEdx per rr tra 100 e 150 per tutte le tracce
    for(int trk=0; trk<int(rr_tot.size()); trk++)
            {
                if(track_length_tot[trk]>=150.)
                {
                    for(int hit=0; hit<int(rr_tot[trk].size()); hit++)
                    {
                        if(rr_tot[trk][hit]<150. && rr_tot[trk][hit]>100.)
                        {
                            dEdx_100_150_tot->Fill(dE_tot[trk][hit]);
                        }
                    }
                }
            } 

    //controllo - dEdx per rr tra 100 e 150 per le tracce filtrate con la mediana
    for(int trk=0; trk<int(rr.size()); trk++)
        {
            if(track_length[trk]>=150.)
            {
                for(int hit=0; hit<int(rr[trk].size()); hit++)
                {
                    if(rr[trk][hit]<150. && rr[trk][hit]>100.)
                    {
                        dEdx_100_150->Fill(dE[trk][hit]);
                    }
                }
            }
        }    
    
    scarto->Write(0,TObject::kOverwrite);
    dEdx_100_150->Scale(1./dEdx_100_150->Integral());
    dEdx_100_150->Write(0,TObject::kOverwrite);
    dEdx_100_150_tot->Scale(1./dEdx_100_150_tot->Integral());
    dEdx_100_150_tot->Write(0,TObject::kOverwrite);
    dEdx_range_median_mag4->Write(0,TObject::kOverwrite);
    dEdx_range_media_reference->Write(0,TObject::kOverwrite);
    graph->Write(0,TObject::kOverwrite);

    f->Close();
    delete f;
    
}


/*
void visual(std::string particle)
{
    
    TFile *f = TFile::Open("ConfrontoDatiMC.root", "READ");
    TDirectory *d= (TDirectory*)f->Get(particle.c_str());
    TDirectory *d1mc= (TDirectory*)d->Get("mc"); 
    TH2D *h2 = (TH2D*)d1mc->Get("dedx_range");
    
    TFile *file = new TFile("ConfrontoHisto.root", "RECREATE");


    for(int i=3; i<h2->GetNbinsX(); i++)
    {

        cout << h2->GetNbinsX() << " " << MPVsfitMC.size() << " " << MPVsfitDATI.size() <<  endl;
        std::string nome = Form("dEdx_i_%d", i);
        double rr=(i-1)*0.3+0.15;

            d1mc->cd();
            TH1D *hmc = (TH1D*)d1mc->Get(nome.c_str()); 
            hmc->Scale(1./hmc->Integral());
            double mpvmc = MPVsfitMC[i-3];

            TDirectory *d1dati= (TDirectory*)d->Get("dati"); 
            d1dati->cd();
            TH1D *hdati = (TH1D*)d1dati->Get(nome.c_str()); 
            hdati->Scale(1./hdati->Integral());
            double mpvdati = MPVsfitDATI[i-3];

            TCanvas *canvas = new TCanvas(Form("Confronto dEdx %s rr [%f , %f] cm", particle.c_str(), (i-1)*0.3 , i*0.3));
            hmc->Draw("hist");
            hmc->SetLineWidth(2);
            hmc->SetLineColor(kRed);
            hdati->Draw("hist same");
            hdati->SetLineWidth(2);
            hdati->SetLineColor(kBlue);
            //hdati->SetLineStyle(kDashed);
            int bin_mpv_mc = hmc->FindBin(mpvmc);
            TArrow *amc = new TArrow(mpvmc, 0, mpvmc, hmc->GetBinContent(bin_mpv_mc));
            amc->SetLineStyle(kDashed);
            amc->SetLineColor(kRed);
            amc->Draw("same");
            int bin_mpv_dati = hmc->FindBin(mpvdati);
            TArrow *adati = new TArrow(mpvdati, 0, mpvdati, hdati->GetBinContent(bin_mpv_dati));
            adati->SetLineStyle(kDashed);
            adati->SetLineColor(kBlue);
            adati->Draw("same");
            TLegend *l = new TLegend();
            l->AddEntry(hmc, "MC");
            l->AddEntry(hdati, "DATI");
            l->AddEntry(amc, "MPV MC");
            l->AddEntry(adati, "MPV DATI");
            l->Draw("same");
            

            file->cd();

            canvas->Write(0,TObject::kOverwrite);

    }

    
}
*/

/*
TGraph* derivata()
{

    TH2D* h2d = histo2d("muon");

    double dx = h2d->GetXaxis()->GetBinWidth(1); //assuming equal spaces bins
    TGraph *g = new TGraph();

    std::vector<double> derivate;
    std::vector<double> rr;
    std::vector<double> sigma;

    for(int nbin=1+2; nbin<=h2d->GetNbinsX()-2; nbin++)
    {
        TH1D* dEdxim2 = h2d->ProjectionY(Form("dEdxim2_%d",nbin-2), nbin-2, nbin-2 );
        TH1D* dEdxim1 = h2d->ProjectionY(Form("dEdxim1_%d",nbin-1), nbin-1, nbin-1);
        TH1D* dEdxip1 = h2d->ProjectionY(Form("dEdxip1_%d", nbin+1), nbin+1, nbin+1);
        TH1D* dEdxip2 = h2d->ProjectionY(Form("dEdxip2_%d", nbin+2), nbin+2, nbin+2);

        TH1D* dEdx_i = h2d->ProjectionY(Form("dEdx_i_%d", nbin), nbin, nbin);

        double der= (-dEdxip2->GetMean()+8*dEdxip1->GetMean()-8*dEdxim1->GetMean()+dEdxim2->GetMean())/(12*dx);
        g->SetPoint(g->GetN(), h2d->GetXaxis()->GetBinCenter(nbin), der);
        rr.push_back(h2d->GetXaxis()->GetBinCenter(nbin));
        derivate.push_back(der);
        sigma.push_back(dEdx_i->GetStdDev());
    }

    double x_low=h2d->GetXaxis()->GetXmax()+2*h2d->GetNbinsX();
    double x_high= h2d->GetXaxis()->GetXmax()-2*h2d->GetNbinsX();
    TSpline3 *spline3 = new TSpline3("spline_della_derivata",g);

    g->SetMarkerStyle(7);

    TCanvas *c = new TCanvas();
    g->Draw("A p");
    spline3->Draw("same");
    spline3->SetLineWidth(2);
    spline3->SetLineColor(kRed);
    


    return g;
}



struct InfoVariation
{
    std::vector<double> rr;
    std::vector<double> sigma;
    std::vector<double> derivate;
    std::vector<double> dx;
};

InfoVariation derivata_light()
{
    TH2D* h2d = histo2d("muon");

    std::vector<double> derivate;
    std::vector<double> rr;
    std::vector<double> sigma;
    std::vector<double> dx_vec;

    TGraph *g = new TGraph();

    TFile *file = new TFile("fitDTCBbinning.root", "RECREATE");

    for(int nbin=1+1; nbin<=h2d->GetNbinsX()-1; nbin++)
    {
        double dx = h2d->GetXaxis()->GetBinWidth(nbin);
        dx_vec.push_back(dx);

        

        TH1D* dEdxim1 = h2d->ProjectionY(Form("dEdxim1_%d",nbin-1), nbin-1, nbin-1);
        TH1D* dEdxip1 = h2d->ProjectionY(Form("dEdxip1_%d", nbin+1), nbin+1, nbin+1);
        TH1D* dEdx_i = h2d->ProjectionY(Form("dEdx_i_%d", nbin), nbin, nbin);

        double x1 = h2d->GetXaxis()->GetBinCenter(nbin);
        double x0 = h2d->GetXaxis()->GetBinCenter(nbin-1);
        double x2 = h2d->GetXaxis()->GetBinCenter(nbin+1);


        double L0 = ((2 * x1 - x2 - x1) / ((x0 - x1) * (x0 - x2)));
        double L1 = ((2 * x1 - x0 - x2) / ((x1 - x0) * (x1 - x2)));
        double L2 = ((2 * x1 - x0 - x1) / ((x2 - x0) * (x2 - x1)));

        double der= dEdxim1->GetMean() * L0 + dEdx_i->GetMean() * L1 + dEdxip1->GetMean() * L2;
        g->SetPoint(g->GetN(), h2d->GetXaxis()->GetBinCenter(nbin), der);
        rr.push_back(h2d->GetXaxis()->GetBinCenter(nbin));
        derivate.push_back(der);

        
        std::vector<double> fitResults = fitLangaulight(dEdx_i, 0.8, 15., file, "dscb", nbin);
        sigma.push_back(fitResults[4]);

        cout << "bin= " << nbin << " rr= " << h2d->GetXaxis()->GetBinCenter(nbin) << " derivata= " << der << " sigma fit= " << fitResults[4] << " sigma composed= " << std::sqrt(pow(fitResults[4],2)+pow(fitResults[2],2)) << " rms= " << dEdx_i->GetStdDev() << " step= " << fitResults[4]/der<<  endl;
    }

    g->Write();

    InfoVariation InfoV;
    InfoV.rr=rr;
    InfoV.sigma=sigma;
    InfoV.derivate=derivate;
    InfoV.dx=dx_vec;

    file->Close();
    delete file;

    return InfoV;

}


int findIndex(const std::vector<double>& vec, double val) {
    int closestIndex = -1;
    double minDiff = std::numeric_limits<double>::max();

    for (size_t i = 0; i < vec.size(); ++i) {
        double diff = std::abs(vec[i] - val);
        if (diff < minDiff) {
            minDiff = diff;
            closestIndex = static_cast<int>(i);
        }
    }

    return closestIndex;
}

//std::vector<std::vector<double>> FindBinEdges()
void FindBinEdges()
{

    InfoVariation InfoV = derivata_light();

    std::vector <double> coppia_startRR_interval;
    std::vector<std::vector <double>> BinEdges;
    
    double interval = 0;
    interval = InfoV.sigma[0]/(InfoV.derivate[0]);
    coppia_startRR_interval.push_back(InfoV.rr[0]-InfoV.dx[0]/2);
    coppia_startRR_interval.push_back(interval);
    BinEdges.push_back(coppia_startRR_interval);
    coppia_startRR_interval.clear();

    int i=0;
    while( i <int(InfoV.rr.size()))
    {
        interval = InfoV.sigma[i]/(InfoV.derivate[i]);
        double next_point=InfoV.rr[i]+interval;
        int j=findIndex(InfoV.rr, next_point);
        
        cout << i << " " << InfoV.rr[i] << " " << interval << " " << j << " " ;

        if(j<=i)
        {
            cout << "j<i " << endl;

            coppia_startRR_interval.push_back(InfoV.rr[i]);//dimezza la larghezza del bin
            coppia_startRR_interval.push_back(interval);  
            BinEdges.push_back(coppia_startRR_interval);
            coppia_startRR_interval.clear();
            i++;     
        }
        else if(j>i)
        {
            i=j;
            interval = InfoV.sigma[i]/TMath::Abs(InfoV.derivate[i]);

            cout << interval << " " << endl;

            coppia_startRR_interval.push_back(InfoV.rr[i]-InfoV.dx[i]/2);
            coppia_startRR_interval.push_back(interval);
            BinEdges.push_back(coppia_startRR_interval);
            coppia_startRR_interval.clear();
        }
        if(j>=(InfoV.rr.size()-2)){break;}
        
    }

    coppia_startRR_interval.push_back(25);
    coppia_startRR_interval.push_back(0);
    BinEdges.push_back(coppia_startRR_interval);
    coppia_startRR_interval.clear();

    cout << endl;

    for(int k=0; k<int(BinEdges.size()); k++)
    {
        cout << BinEdges[k][0] << " " << BinEdges[k][1] << endl;
    }

    //return BinEdges;

}
*/

void dEdxAngle()
{
    std::array<double,4> angles ={0., 50., 70., 90.};
    std::array<std::string,2> datatype = {"mc", "dati"};
    std::array<std::string,2> particle = {"muon", "proton"};
    double threshold_med;

    TFile *f = new TFile("dEdxAngle.root", "RECREATE");

    for(double ang : angles)
    {
        TDirectory *d = (TDirectory*)f->mkdir(Form("%f",ang));
        d->cd();
        for(std::string data : datatype)
        {
            TDirectory *d1 = (TDirectory*)d->mkdir(data.c_str());
            d1->cd();
            for(std::string par : particle)
            {
                cout << ang << " " << data << " " << par << endl;
                TDirectory *d2 = (TDirectory*)d1->mkdir(par.c_str());
            
                threshold_med = (par=="muon") ? 3.2 : 8.;

                EventsData dat = load_data(data, par, "Np", "full");
                std::string dE_branch_name = (par=="muon") ? "dE_mu" : "dE_pro";
                std::string rr_branch_name = (par=="muon") ? "rr_mu" : "rr_pro";
                dat.tree->SetBranchStatus("*", 0);
                dat.tree->SetBranchStatus(dE_branch_name.c_str(), 1);
                dat.tree->SetBranchStatus(rr_branch_name.c_str(), 1);
                dat.tree->SetBranchStatus("dirx", 1);
                d2->cd();
                TTree *tree = new TTree("tree_fitres","tree_fitres");
                double mpv ,mpv_err,resrange;
                tree->Branch("rr", &resrange );
                tree->Branch("mpv", &mpv);
                tree->Branch("mpv_err", &mpv_err);
                TH1D *htheta_drift = new TH1D("htheta_drift", "", 90, 0, 90);
                TH1D *htheta_drift_cut = new TH1D("htheta_drift_cut", "", 90, 0, 90);
                TH1D *hmedian = new TH1D("hmedian", "", 300, 0, 30);
                TH1D *hmedian_cut = new TH1D("hmedian_cut", "", 300, 0, 30);

                TH1D *MPV_range = new TH1D(Form("MPV_rr_%d_%s_%s", int(ang), data.c_str(), par.c_str() ), Form("MPV_rr_%d_%s_%s", int(ang), data.c_str(), par.c_str() ), 50, 0, 50 );

                for(int rr=2; rr<=50; rr++)
                {
                    TH1D *dEdx = new TH1D(Form("dEdx_rr%d", rr), Form("dEdx_rr%d", rr), 250, 0, 30 );
                    std::vector<double> dEdx_vec;

                    for(int track=0; track<dat.tree->GetEntries(); track++)
                    {
                        dat.tree->GetEntry(track);

                        std::vector<double> dEdx_temp;
                        for(int hit=0; hit<int(dat.track.dE->size())-2; hit++)
                        {
                            if(dat.track.rr->at(hit)<5.)
                            dEdx_temp.push_back(dat.track.dE->at(hit));
                        }

                        if(dEdx_temp.size()==0)continue;
                        double med = mediana(dEdx_temp);
                        hmedian->Fill(med);
                        if(med<threshold_med)continue;
                        hmedian_cut->Fill(med);

                        for(int hit=0; hit<int(dat.track.dE->size()); hit++)
                        {
                            if(dat.track.rr->at(hit)<rr && dat.track.rr->at(hit)>=rr-1)
                            {
                                double theta_drift = 180./M_PI*acos(abs(dat.track.dirx));
                                htheta_drift->Fill(theta_drift);
                                if(ang!=angles[0])
                                {
                                    if(theta_drift<ang && theta_drift >=ang-30)
                                    {   
                                        htheta_drift_cut->Fill(theta_drift);
                                        dEdx->Fill(dat.track.dE->at(hit));
                                        dEdx_vec.push_back(dat.track.dE->at(hit));
                                    }
                                }
                                else {dEdx->Fill(dat.track.dE->at(hit)); dEdx_vec.push_back(dat.track.dE->at(hit));}    
                            }
                            
                        }
                    }


                    int binMax=0;
                    double maxVal=0;
                    for(int j=1; j<dEdx->GetNbinsX(); j++)
                    {
                        double y=dEdx->GetBinContent(j);
                        if(y>maxVal)
                        {
                        maxVal=y;
                        binMax=j;
                        }
                    }
                    cout << dEdx->GetEntries() << " ";
                    double x_low=dEdx->GetBinCenter(binMax)-0.80*dEdx->GetStdDev();
                    double x_high=dEdx->GetBinCenter(binMax)+4*dEdx->GetStdDev(); 
                    FitOutput fitout = fitLangaulight_unbinned(dEdx,dEdx_vec, x_low, x_high, d2, "langau", rr); 
                    mpv = fitout.fitResults[0];
                    mpv_err = fitout.fitResults[1];
                    resrange=rr-0.5;
                    MPV_range->SetBinContent(rr,mpv );
                    MPV_range->SetBinError(rr,mpv_err);
                    cout << resrange << " " << mpv << " " << mpv_err << endl;
                    tree->Fill();
                    //dEdx->Scale(1./dEdx->Integral());
                    //dEdx->Write(Form("dEdx_rr%d", rr),TObject::kOverwrite);

                    delete dEdx;

                }
                MPV_range->Write(0,TObject::kOverwrite);
                tree->Write(0,TObject::kOverwrite);
                htheta_drift->Scale(1./htheta_drift->Integral());
                htheta_drift_cut->Scale(1./htheta_drift_cut->Integral());
                htheta_drift->Write(0,TObject::kOverwrite);
                htheta_drift_cut->Write(0,TObject::kOverwrite);
                hmedian->Scale(1./hmedian->Integral());
                hmedian_cut->Scale(1./hmedian_cut->Integral());
                hmedian->Write(0,TObject::kOverwrite);
                hmedian_cut->Write(0,TObject::kOverwrite);  
            }
        }
    }
}


Double_t inverted_modified_box(Double_t *x, Double_t *par)
{
    Double_t electric_field = par[0];
    Double_t density = par[1];
    Double_t Wion = par[2];
    Double_t gain = par[3];

    Double_t alpha = par[4];
    Double_t beta = par[5];
    Double_t R = beta/(electric_field*density);
    Double_t exponent = gain*R*Wion;
    return (TMath::Exp(x[0]*exponent)-alpha)/R; 
}

Double_t inverted_ellipsoidal_modified_box(Double_t *x, Double_t *par)
{
    Double_t electric_field = par[0];
    Double_t density = par[1];
    Double_t Wion = par[2];
    Double_t gain = par[3];

    Double_t alpha = par[4];
    Double_t beta = par[5];
    Double_t r = par[6];
    Double_t R = beta/((electric_field*density)*std::sqrt(std::pow(TMath::Sin(x[1]),2)+pow(TMath::Cos(x[1]),2)/r));
    Double_t exponent = gain*R*Wion;
    return (TMath::Exp(x[0]*exponent)-alpha)/R;
}

void whichRecombination()
{
    TF1 * f_inverted_modified_box = new TF1("inverted_modified_box",inverted_modified_box, 0., 10000, 6);
    f_inverted_modified_box->SetParameters(0.5, 1.383, 26.3e-6, 75., 0.93, 0.212);
    TF2  *f_inverted_ellipsoidal_box_model = new TF2("inverted_ellipsoidal_box_model",inverted_ellipsoidal_modified_box, 0, 10000, 0, TMath::Pi()/2, 7 );  
    f_inverted_ellipsoidal_box_model->SetParameters(0.5, 1.383, 26.3e-6, 75., 0.904, 0.204, 1.25);

    EventsData dat = load_data("mc", "muon", "Np", "full");
    
    TH2D *dQdx_recomb_mu_mc = new TH2D("dQdx_recomb_mu_mc", "", 10000, 0, 10000, 600, 0, 60);
    TH2D *g_modified_box = new TH2D("g_modified_box", "", 10000, 0, 10000, 600, 0, 60);
    TH2D *g_elipsoidal = new TH2D("g_elipsoidal", "", 10000, 0, 10000, 600, 0, 60);
    for(int track=0; track<dat.tree->GetEntries(); track++)
    {

        dat.tree->GetEntry(track);
        double phi = acos(abs(dat.track.dirx));

        for(int hit=0; hit<int(dat.track.dQdx->size()); hit++)
        {
            dQdx_recomb_mu_mc->Fill(dat.track.dQdx->at(hit), dat.track.dE->at(hit));
            g_modified_box->Fill(dat.track.dQdx->at(hit),f_inverted_modified_box->Eval(dat.track.dQdx->at(hit))  );
            g_elipsoidal->Fill(dat.track.dQdx->at(hit), f_inverted_ellipsoidal_box_model->Eval(dat.track.dQdx->at(hit), phi) );
        }
    }

    TCanvas *c = new TCanvas();
    dQdx_recomb_mu_mc->Draw("colz");
    dQdx_recomb_mu_mc->GetXaxis()->SetTitle("dQ/dx [ADC]");
    dQdx_recomb_mu_mc->GetYaxis()->SetTitle("dE/dx [MeV/cm]");
    g_modified_box->Draw("scatt same");
    g_modified_box->SetMarkerStyle(7);
    g_modified_box->SetMarkerColor(kRed);
    g_modified_box->SetLineColor(kRed);
    g_elipsoidal->Draw("scatt same");
    g_elipsoidal->SetMarkerStyle(7);
    g_elipsoidal->SetMarkerColor(kGreen);
    g_elipsoidal->SetLineColor(kGreen);
    TLegend *l = new TLegend();
    l->AddEntry(g_modified_box, "modified box model");
    l->AddEntry(g_elipsoidal, "ellipsoidal box model");
    l->Draw("same");
    TFile *f = new TFile("recombinantion_model.root", "RECREATE");
    c->Write(0,TObject::kOverwrite);
}



/*
void plotdQdx()
{
    EventsData dat = load_data("mc", "proton", "Np", "full");

    TH1D *h1 = new TH1D("", "", 10000, 0, 10000);
    for(int track=0; track<dat.tree->GetEntries(); track++)
    {
        dat.tree->GetEntry(track);
        for(int hit=0; hit<int(dat.track.dQdx->size()); hit++)
        {
            h1->Fill(dat.track.dQdx->at(hit));
        }
    }
    h1->Draw();
}
*/

int FindBinMax(TH1D *h)
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

void medianAngle()
{
    std::array<double,4> angles ={0., 50., 70., 90.};
    std::array<std::string,2> datatype = {"mc", "dati"};
    std::array<std::string,2> particle = {"muon", "proton"};
    double threshold_med;

    TFile *f = TFile::Open("dEdxAngle.root", "UPDATE");

    for(double ang : angles)
    {
        for(std::string data : datatype)
        {
            for(std::string par : particle)
            {
                cout << ang << " " << data << " " << par << endl;
            
                threshold_med = (par=="muon") ? 3.2 : 8.;

                EventsData dat = load_data(data, par, "Np", "full");
                std::string dE_branch_name = (par=="muon") ? "dE_mu" : "dE_pro";
                std::string rr_branch_name = (par=="muon") ? "rr_mu" : "rr_pro";
                dat.tree->SetBranchStatus("*", 0);
                dat.tree->SetBranchStatus(dE_branch_name.c_str(), 1);
                dat.tree->SetBranchStatus(rr_branch_name.c_str(), 1);
                dat.tree->SetBranchStatus("dirx", 1);

                TH1D *hmedian = new TH1D("hmedian", "", 300, 0, 30);
                TH1D *hmedian_cut = new TH1D("hmedian_cut", "", 300, 0, 30);

                for(int rr=2; rr<=50; rr++)
                {

                    for(int track=0; track<dat.tree->GetEntries(); track++)
                    {
                        dat.tree->GetEntry(track);

                        double theta_drift = 180./M_PI*acos(abs(dat.track.dirx));
                        if(ang!=angles[0])
                        {
                            if(theta_drift<ang && theta_drift >=ang-30)
                            {
                                std::vector<double> dEdx_temp;
                                for(int hit=0; hit<int(dat.track.dE->size())-2; hit++)
                                {
                                    if(dat.track.rr->at(hit)<5.)
                                    dEdx_temp.push_back(dat.track.dE->at(hit));
                                }

                                if(dEdx_temp.size()==0)continue;
                                double med = mediana(dEdx_temp);
                                hmedian->Fill(med);
                            }
                        }
                        else
                        {
                            std::vector<double> dEdx_temp;
                            for(int hit=0; hit<int(dat.track.dE->size())-2; hit++)
                            {
                                if(dat.track.rr->at(hit)<5.)
                                dEdx_temp.push_back(dat.track.dE->at(hit));
                            }

                            if(dEdx_temp.size()==0)continue;
                            double med = mediana(dEdx_temp);
                            hmedian->Fill(med);
                        }
                    }

                }
                hmedian->Scale(1./hmedian->Integral()); 
                f->cd();
                hmedian->Write(Form("mediana_%s_%s_%d", par.c_str(), data.c_str(), int(ang)),TObject::kOverwrite);
            }
        }
    }
}

// kBlue = 600
// kBlack = 1
// kAzure = 860
// kRed = 632
// kOrange = 800
// kMagenta = 616
// kGreen = 416 
// kPink = 900

////////// plot MPV in funzione dell'angolo /////////////////////////////////////////////////////////////////////////////////

void plotMedian()
{
    TFile *f = TFile::Open("dEdxAngle.root", "UPDATE");

    std::array<std::string,2> particles ={"muon", "proton"};
    std::array<std::string,2> datas ={"dati", "mc"};
    
    for(auto particle : particles)
    {
        for(auto data : datas)
        {
            auto h0= (TH1D*)f->Get(Form("mediana_%s_%s_0", particle.c_str(), data.c_str()));
            auto h50= (TH1D*)f->Get(Form("mediana_%s_%s_50", particle.c_str(), data.c_str()));
            auto h70= (TH1D*)f->Get(Form("mediana_%s_%s_70", particle.c_str(), data.c_str()));
            auto h90= (TH1D*)f->Get(Form("mediana_%s_%s_90", particle.c_str(), data.c_str()));

            h0->Rebin(4);
            h50->Rebin(4);
            h70->Rebin(4);
            h90->Rebin(4);

            TCanvas *c = new TCanvas();

            h0->Draw("hist");
            h0->GetYaxis()->SetRangeUser(0.,0.5);
            h0->GetYaxis()->SetTitle("counts (area normalized)");
            h0->GetXaxis()->SetTitle("median [MeV/cm]");
            h0->SetLineColor(600);
            h0->SetLineWidth(2);
            h50->Draw("hist same");
            h50->SetLineColor(632);
            h50->SetLineWidth(2);
            h70->Draw("hist same");
            h70->SetLineColor(800);
            h70->SetLineWidth(2);
            h90->Draw("hist same");
            h90->SetLineColor(416);
            h90->SetLineWidth(2);
            TLegend *l = new TLegend();
            l->AddEntry(h0, "no conditions");
            l->AddEntry(h50, "0^{#circ} #leq #theta < 50^{#circ}");
            l->AddEntry(h70, "50^{#circ} #leq #theta < 70^{#circ}");
            l->AddEntry(h90, "70^{#circ} #leq #theta < 90^{#circ}");
            l->Draw("same");

            cout << Form("mediane_%s_%s", particle.c_str(), data.c_str()) << " " << h0->GetBinCenter(FindBinMax(h0)) << " " << h50->GetBinCenter(FindBinMax(h50)) << endl;  
            c->Write(Form("mediane_%s_%s", particle.c_str(), data.c_str()), TObject::kOverwrite);
        }
    }
}
    

std::pair<double, double> diff(double x, double y, double err_x, double err_y) {
    // Calcolo della differenza normalizzata
    double R = 2.0 * (x - y) / (x + y);

    // Derivate parziali
    double dR_dx = 2.0 * y / std::pow(x + y, 2);
    double dR_dy = -2.0 * x / std::pow(x + y, 2);

    // Propagazione dell'errore
    double sigma_R = std::sqrt(std::pow(dR_dx * err_x, 2) + std::pow(dR_dy * err_y, 2));

    return std::make_pair(R, sigma_R);
}

std::pair<double, double> diffdEdx(double x, double y, double err_x, double err_y) {
    // Calcolo della differenza normalizzata
    double R = 2.0 * (x - y) / (x + y);

    // Derivate parziali
    double dR_dx = 2.0 * y / std::pow(x + y, 2);
    double dR_dy = -2.0 * x / std::pow(x + y, 2);

    // Propagazione dell'errore
    double sigma_R = 2*std::sqrt(std::pow(dR_dx * err_x, 2) + std::pow(dR_dy * err_y, 2));

    return std::make_pair(R, sigma_R);
}

TH1D* GetHisto(TTree *tree_dati, TTree *tree_mc, int color, std::string angle, std::string particle)
{
    double rr_dati;
    double mpv_dati;
    double mpv_err_dati;

    double rr_mc;
    double mpv_mc;
    double mpv_err_mc;

    tree_dati->SetBranchAddress("rr", &rr_dati);
    tree_dati->SetBranchAddress("mpv", &mpv_dati);
    tree_dati->SetBranchAddress("mpv_err", &mpv_err_dati);

    tree_mc->SetBranchAddress("rr", &rr_mc);
    tree_mc->SetBranchAddress("mpv", &mpv_mc);
    tree_mc->SetBranchAddress("mpv_err", &mpv_err_mc);

    TH1D * h = new TH1D(Form("h_%s_%s", angle.c_str(), particle.c_str() ), "", 50, 0, 50);
    for(int i=2; i<=h->GetNbinsX(); i++)
    {
        tree_dati->GetEntry(i-2);
        tree_mc->GetEntry(i-2);
        if(particle=="muon" && angle=="0" && rr_mc==40.5 ){mpv_mc=1.82483; }
        if(particle=="muon" && angle=="0" && rr_mc==46.5 ){mpv_mc=1.78585;}
        if(particle=="muon" && angle=="50" && rr_mc==29.5 ){mpv_mc=1.911225;}
        if(particle=="muon" && angle=="50" && rr_mc==48.5 ){mpv_mc=1.742055;}
        if(particle=="muon" && angle=="90" && rr_dati==34.5 ){mpv_dati=1.91752;}
        if(particle=="proton" && i==2)continue;
        std::pair<double, double> differenza = diff(mpv_dati, mpv_mc, mpv_err_dati, mpv_err_mc);
        h->SetBinContent(i, differenza.first);
        h->SetBinError(i, differenza.second);
    }

    h->SetLineWidth(2);
    h->SetMarkerColor(color);
    h->SetLineColor(color);
    h->SetMarkerStyle(7);

    return h;

}

TGraphErrors* GetGraph(TTree *tree_dati, TTree *tree_mc, int color, std::string angle, std::string particle)
{   
    double rr_dati;
    double mpv_dati;
    double mpv_err_dati;

    double rr_mc;
    double mpv_mc;
    double mpv_err_mc;

    tree_dati->SetBranchAddress("rr", &rr_dati);
    tree_dati->SetBranchAddress("mpv", &mpv_dati);
    tree_dati->SetBranchAddress("mpv_err", &mpv_err_dati);

    tree_mc->SetBranchAddress("rr", &rr_mc);
    tree_mc->SetBranchAddress("mpv", &mpv_mc);
    tree_mc->SetBranchAddress("mpv_err", &mpv_err_mc);

    TGraphErrors *g = new TGraphErrors();
    for(int i=0; i<tree_dati->GetEntries(); i++)
    {
        tree_dati->GetEntry(i);
        tree_mc->GetEntry(i);
        if(particle=="muon" && angle=="0" && rr_mc==40.5 ){mpv_mc=1.82483;}
        if(particle=="muon" && angle=="0" && rr_mc==46.5 ){mpv_mc=1.78585;}
        if(particle=="muon" && angle=="50" && rr_mc==29.5 ){mpv_mc=1.911225;}
        if(particle=="muon" && angle=="50" && rr_mc==48.5 ){mpv_mc=1.742055;}
        if(particle=="muon" && angle=="90" && rr_dati==34.5 ){mpv_dati=1.91752;}
        if(particle=="proton" && i==0)continue;
        std::pair<double, double> differenza = diff(mpv_dati, mpv_mc, mpv_err_dati, mpv_err_mc);
        g->SetPoint(i, (mpv_dati+mpv_mc)/2, differenza.first);
        g->SetPointError(i, std::sqrt(pow(mpv_err_dati,2)+pow(mpv_err_mc,2)), differenza.second);
    }

    g->SetLineWidth(2);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    g->SetMarkerStyle(7);

    return g;

}

void plotAngle()
{
    TFile *f = TFile::Open("dEdxAngle.root", "UPDATE");
    TTree *tree_dati_0 = (TTree*)f->Get("0.000000/dati/muon/tree_fitres");
    TTree *tree_mc_0 = (TTree*)f->Get("0.000000/mc/muon/tree_fitres");
    TTree *tree_dati_50 = (TTree*)f->Get("50.000000/dati/muon/tree_fitres");
    TTree *tree_mc_50 = (TTree*)f->Get("50.000000/mc/muon/tree_fitres");
    TTree *tree_dati_70 = (TTree*)f->Get("70.000000/dati/muon/tree_fitres");
    TTree *tree_mc_70 = (TTree*)f->Get("70.000000/mc/muon/tree_fitres");
    TTree *tree_dati_90 = (TTree*)f->Get("90.000000/dati/muon/tree_fitres");
    TTree *tree_mc_90 = (TTree*)f->Get("90.000000/mc/muon/tree_fitres");
    TTree *tree_pro_dati_0 = (TTree*)f->Get("0.000000/dati/proton/tree_fitres");
    TTree *tree_pro_mc_0 = (TTree*)f->Get("0.000000/mc/proton/tree_fitres");
    TTree *tree_pro_dati_50 = (TTree*)f->Get("50.000000/dati/proton/tree_fitres");
    TTree *tree_pro_mc_50 = (TTree*)f->Get("50.000000/mc/proton/tree_fitres");
    TTree *tree_pro_dati_70 = (TTree*)f->Get("70.000000/dati/proton/tree_fitres");
    TTree *tree_pro_mc_70 = (TTree*)f->Get("70.000000/mc/proton/tree_fitres");
    TTree *tree_pro_dati_90 = (TTree*)f->Get("90.000000/dati/proton/tree_fitres");
    TTree *tree_pro_mc_90 = (TTree*)f->Get("90.000000/mc/proton/tree_fitres");

    TH1D* h0 = GetHisto(tree_dati_0,tree_mc_0,860,"0","muon" );
    TH1D* h50 = GetHisto(tree_dati_50,tree_mc_50,632,"50","muon");
    TH1D* h70 = GetHisto(tree_dati_70,tree_mc_70,800,"70","muon");
    TH1D* h90 = GetHisto(tree_dati_90,tree_mc_90, 416 ,"90","muon");
    TH1D* h0_pro = GetHisto(tree_pro_dati_0,tree_pro_mc_0,860,"0","proton" );
    TH1D* h50_pro = GetHisto(tree_pro_dati_50,tree_pro_mc_50,632,"50","proton");
    TH1D* h70_pro = GetHisto(tree_pro_dati_70,tree_pro_mc_70,800,"70","proton");
    TH1D* h90_pro = GetHisto(tree_pro_dati_90,tree_pro_mc_90, 416,"90","proton" );

    TGraphErrors* g0 = GetGraph(tree_dati_0,tree_mc_0,860,"0","muon" );
    TGraphErrors* g50 = GetGraph(tree_dati_50,tree_mc_50,632,"0","muon");
    TGraphErrors* g70 = GetGraph(tree_dati_70,tree_mc_70,800,"0","muon");
    TGraphErrors* g90 = GetGraph(tree_dati_90,tree_mc_90, 416,"0","muon" );
    TGraphErrors* g0_pro = GetGraph(tree_pro_dati_0,tree_pro_mc_0,860,"0","proton" );
    TGraphErrors* g50_pro = GetGraph(tree_pro_dati_50,tree_pro_mc_50,632,"0","proton");
    TGraphErrors* g70_pro = GetGraph(tree_pro_dati_70,tree_pro_mc_70,800,"0","proton");
    TGraphErrors* g90_pro = GetGraph(tree_pro_dati_90,tree_pro_mc_90, 416 ,"0","proton");

    g0->Write("g0");
    g50->Write("g50");
    g70->Write("g70");
    g90->Write("g90");
    g0_pro->Write("g0_pro");
    g50_pro->Write("g50_pro");
    g70_pro->Write("g70_pro");
    g90_pro->Write("g90_pro");
    h0->Write("h0");
    h50->Write("h50");
    h70->Write("h70");
    h90->Write("h90");
    h0_pro->Write("h0_pro");
    h50_pro->Write("h50_pro");
    h70_pro->Write("h70_pro");
    h90_pro->Write("h90_pro");


    TCanvas *c = new TCanvas();
    h0->Draw("P");
    h0->GetXaxis()->SetTitle("Residual Range [cm]");
    h0->GetYaxis()->SetTitle("(Data-MC)/average");
    h50->Draw("P same");
    h70->Draw("P same");
    h90->Draw("P same");
    TLegend *l = new TLegend();
    l->AddEntry(h0, "no conditions");
    l->AddEntry(h50, "0^{#circ} #leq #theta < 50^{#circ}");
    l->AddEntry(h70, "50^{#circ} #leq #theta < 70^{#circ}");
    l->AddEntry(h90, "70^{#circ} #leq #theta < 90^{#circ}");
    l->Draw("same");
    c->Write("canvas_mu",TObject::kOverwrite);

    TCanvas *c_pro = new TCanvas();
    h0_pro->Draw("P");
    h0_pro->GetXaxis()->SetTitle("Residual Range [cm]");
    h0_pro->GetYaxis()->SetTitle("(Data-MC)/average");
    h50_pro->Draw("P same");
    h70_pro->Draw("P same");
    h90_pro->Draw("P same");
    TLegend *l_pro = new TLegend();
    l_pro->AddEntry(h0_pro, "no conditions");
    l_pro->AddEntry(h50_pro, "0^{#circ} #leq #theta < 50^{#circ}");
    l_pro->AddEntry(h70_pro, "50^{#circ} #leq #theta < 70^{#circ}");
    l_pro->AddEntry(h90_pro, "70^{#circ} #leq #theta < 90^{#circ}");
    l_pro->Draw("same");
    c_pro->Write("canvas_pro",TObject::kOverwrite);

    TCanvas *cg = new TCanvas();
    g0->Draw("AP");
    g0->GetXaxis()->SetTitle("average [MeV/cm]");
    g0->GetYaxis()->SetTitle("(Data-MC)/average");
    g50->Draw("P same");
    g70->Draw("P same");
    g90->Draw("P same");
    TLegend *lg = new TLegend();
    lg->AddEntry(g0, "no conditions");
    lg->AddEntry(g50, "0^{#circ} #leq #theta < 50^{#circ}");
    lg->AddEntry(g70, "50^{#circ} #leq #theta < 70^{#circ}");
    lg->AddEntry(g90, "70^{#circ} #leq #theta < 90^{#circ}");
    lg->Draw("same");
    cg->Write("gcanvas_mu",TObject::kOverwrite);

    TCanvas *cg_pro = new TCanvas();
    g0_pro->Draw("AP");
    g0_pro->GetXaxis()->SetTitle("average [MeV/cm]");
    g0_pro->GetYaxis()->SetTitle("(Data-MC)/average");
    g50_pro->Draw("P same");
    g70_pro->Draw("P same");
    g90_pro->Draw("P same");
    TLegend *lg_pro = new TLegend();
    lg_pro->AddEntry(g0_pro, "no conditions");
    lg_pro->AddEntry(g50_pro, "0^{#circ} #leq #theta < 50^{#circ}");
    lg_pro->AddEntry(g70_pro, "50^{#circ} #leq #theta < 70^{#circ}");
    lg_pro->AddEntry(g90_pro, "70^{#circ} #leq #theta < 90^{#circ}");
    lg_pro->Draw("same");
    cg_pro->Write("gcanvas_pro",TObject::kOverwrite);

    TCanvas * gdiff_dEdx_tot = new TCanvas();
    g0->Draw("AP");
    g0->GetXaxis()->SetLimits(0., 16.);
    g0->GetXaxis()->SetTitle("average [MeV/cm]");
    g0->GetYaxis()->SetTitle("(Data-MC)/average");
    g50->Draw("P same");
    g70->Draw("P same");
    g90->Draw("P same");
    g0_pro->Draw("P same");
    g50_pro->Draw("P same");
    g70_pro->Draw("P same");
    g90_pro->Draw("P same");
    TLegend *ltot = new TLegend();
    ltot->AddEntry(g0, "no conditions");
    ltot->AddEntry(g50, "0^{#circ} #leq #theta < 50^{#circ}");
    ltot->AddEntry(g70, "50^{#circ} #leq #theta < 70^{#circ}");
    ltot->AddEntry(g90, "70^{#circ} #leq #theta < 90^{#circ}");
    ltot->Draw("same");
    gdiff_dEdx_tot->Write("tot_diff_dEdx", TObject::kOverwrite);


}


Double_t oscillationProb(Double_t *x, Double_t *par)
{
    //Double_t first = pow( TMath::Sin(2*par[0]) ,2);
    Double_t first = par[0];
    Double_t second = pow(TMath::Sin((1.27*par[1]*par[2])/x[0]),2); //x is neutrino energy, unit in km/GeV= m/MeV
    return second;
}

void plotOscillationProb()
{
    TF1 *f = new TF1("oscProb", oscillationProb, 0., 50, 3);
    f->SetParameters(0.005,0.5, 730);

    TF1 *f1 = new TF1("oscProb1", oscillationProb, 0., 50, 3);
    f1->SetParameters(0.005,0.1, 730);

    TF1 *f2 = new TF1("oscProb2", oscillationProb, 0., 50, 3);
    f2->SetParameters(0.005,0.01, 730);


    TCanvas *c = new TCanvas();
    f->SetNpx(100000);
    f1->SetNpx(100000);
    f2->SetNpx(100000); 
    f->GetXaxis()->SetTitle("neutrino energy [GeV]");
    f->GetYaxis()->SetTitle("sin^{2}(#frac{1.27*#Delta m^{2}*730}/{E_{#nu}})");
    f->Draw();
    f1->SetLineColor(kAzure);
    f1->SetLineStyle(kDashed);
    f1->SetLineWidth(2);
    f2->SetLineWidth(2);
    f1->Draw("same");
    f2->SetLineColor(kGreen);
    f2->SetLineStyle(kDashed);
    f2->Draw("same");

    TLegend *l = new TLegend();
    l->AddEntry(f, "#Delta m^{2} = 0.5");
    l->AddEntry(f1, "#Delta m^{2} = 0.1");
    l->AddEntry(f2, "#Delta m^{2} = 0.01");
    l->Draw("same");
    
    TFile *file = new TFile("oscProb.root", "RECREATE");
    c->Write();

    //cout << pow(TMath::Sin(2*2.027403614),2) << endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
struct HistInfo
{
    double mean;
    double sigma;
    double mpv;
    double sigmaG;
    double sigmaL;
    double start_rr;
    double rr_width;
    int rrNbins;
};

struct Hist
{
    HistInfo histInfo;
    TH1D *histo_pdf = nullptr;
};



Hist dEdxFixedRR(TH2D *dedx_range_endpoint_min1, TDirectory *d, int rrBin){

    d->cd();
    TH1D* hProj = dedx_range_endpoint_min1->ProjectionY("hProj", rrBin, rrBin);

    Hist hist;

    hist.histInfo.rrNbins = dedx_range_endpoint_min1->GetNbinsX();


    hist.histInfo.mean = hProj->GetMean();
    hist.histInfo.sigma = hProj->GetStdDev();
    hist.histInfo.start_rr = dedx_range_endpoint_min1->GetXaxis()->GetBinLowEdge(rrBin);
    hist.histInfo.rr_width = dedx_range_endpoint_min1->GetXaxis()->GetBinWidth(rrBin);

    TH1D *dEdx_prob_density = new TH1D(Form("dEdx_probability_density%d",rrBin), "", hProj->GetNbinsX(), hProj->GetXaxis()->GetXmin(), hProj->GetXaxis()->GetXmax());
    for(int bin=1; bin<=dEdx_prob_density->GetNbinsX(); bin++)
    {
        double prob=0;
        prob = hProj->GetBinContent(bin)/(hProj->GetEntries()*hProj->GetBinWidth(bin));
        dEdx_prob_density->SetBinContent(bin, prob );
    }

    hist.histo_pdf= dEdx_prob_density;

    TCanvas *canvas_dEdx_prob_density = new TCanvas();
    dEdx_prob_density->Draw("p");
    dEdx_prob_density->SetMarkerSize(7);
    canvas_dEdx_prob_density->Write(Form("dEdx_probability_density%d",rrBin ), TObject::kOverwrite);

    delete hProj;
    delete canvas_dEdx_prob_density;

    return hist;
}


void likelihood(std::string particle_hypothesis, std::string particle)
{
    TFile *file_di_output = TFile::Open("OUTPUT_likelihood.root", "UPDATE");
    
    std::string directory_name;
    if(particle_hypothesis=="muon"){directory_name="muon_hypothesis";}
    else if(particle_hypothesis=="proton"){directory_name="proton_hypothesis";}
    TDirectory *d_hyp = (TDirectory*)file_di_output->Get(directory_name.c_str());

    TH2D * h2 = histo2d(particle_hypothesis, d_hyp);

    TGraph * gderivata = derivata(); //ricordardi che ci va il th2d qui dentro
    d_hyp->cd();
    TCanvas *canvas_derivata = new TCanvas("canvas_derivata");
    gderivata->Draw("p");
    canvas_derivata->Write(0,TObject::kOverwrite);

    int h2Nbins = h2->GetNbinsX();

    EventsData dat = load_data("dati", particle, "Np", "25only");
    cout << dat.tree->GetEntries() << " total tracks" << endl;

    TH1D *single_probability = new TH1D("single_probability", "", 50, 0, 2.);

    TH1D *hLogL = new TH1D("hLogL", "", 200, 0., 0.);

    std::vector<Hist> vec_dEdxFixedRR;
    vec_dEdxFixedRR.reserve(h2Nbins);

    for(int rrbin = 1; rrbin <= h2Nbins; rrbin++) 
    {
        vec_dEdxFixedRR.push_back(dEdxFixedRR(h2, d_hyp, rrbin));
    }


    for(int track=0; track<dat.tree->GetEntries(); track++)
    {
        double like=1;

        for(int rrbin =1; rrbin<=h2Nbins; rrbin++)
        {

            
            cout << "mean " << hist.histInfo.mean << endl;
            cout << "sigma " << hist.histInfo.sigma << endl;
            cout << "start rr (included) " << hist.histInfo.start_rr << endl;
            cout << "rr interval width (upper edge excluded) " << hist.histInfo.rr_width << endl;
            

            Hist &hist = vec_dEdxFixedRR[rrbin - 1];
       
            dat.tree->GetEntry(track);

            for(int hit=0; hit<int(dat.track.rr->size()); hit++)
            { 
                if(dat.track.rr->at(hit)>=hist.histInfo.start_rr && dat.track.rr->at(hit)<hist.histInfo.start_rr + hist.histInfo.rr_width && dat.track.rr->at(hit)<10.2)
                {
                    int bin_lik = hist.histo_pdf->FindBin(dat.track.dE->at(hit));
                    single_probability->Fill(hist.histo_pdf->GetBinContent(bin_lik));

                    like=like*hist.histo_pdf->GetBinContent(bin_lik);
                }
            }   

        }

        //hLogL->Fill(2*TMath::Log(like));
        hLogL->Fill(like);

    }
    
    TDirectory *d = (TDirectory*)file_di_output->Get(particle.c_str());
    d->cd();

    single_probability->Scale(1./single_probability->Integral());
    single_probability->Write(0,TObject::kOverwrite);    
    
    hLogL->Scale(1./hLogL->Integral());
    hLogL->Write(0,TObject::kOverwrite);

    file_di_output->Close();
    delete file_di_output;
}
*/

double propagateError(double x, double dx, double y, double dy) {
    double z = (2 * (x - y)) / (x + y);

    // Derivate parziali
    double dz_dx = (4 * y) / std::pow(x + y, 2);
    double dz_dy = (-4 * x) / std::pow(x + y, 2);

    // Propagazione dell'errore
    double dz = std::sqrt(std::pow(dz_dx * dx, 2) + std::pow(dz_dy * dy, 2));

    return dz;
}

TGraphErrors* histoRatio(TH1D* dEdx_dati, TH1D* dEdx_mc)
{
    TGraphErrors *ratio = new TGraphErrors();
    for(int i=1; i<dEdx_dati->GetNbinsX(); i++)
        {
            double x = dEdx_dati->GetBinCenter(i); 
            double ydat = dEdx_dati->GetBinContent(i);
            double ydat_err = dEdx_dati->GetBinError(i);
            double ymc = dEdx_mc->GetBinContent(i);
            double ymc_err = dEdx_mc->GetBinError(i);

            if(ymc==0)continue;
            double val = ydat/ymc;
            //cout << val << endl;
            double val_err = std::sqrt(std::pow(1./ymc*ydat_err,2)+std::pow(ydat/ymc/ymc*ymc_err,2));

            ratio->SetPoint(i-1,x,val);
            ratio->SetPointError(i-1,0,val_err);
        }

    return ratio;
}

Double_t Correction_muon(Double_t *x, Double_t *par)
{
    if(x[0]<=2.) return 0.0120827;
    else if(x[0]>=5.) return 0.0373432;
    else return par[0]*(x[0]-par[1])*(x[0]-par[1])+par[2]*x[0]+par[3];
}


TCanvas* canvasHisto( TH1D *h1, TH1D *h2, const char* title1, const char* title2, const char* doption1, const char* doption2, Color_t color1, Color_t color2, int width1, int width2, const char* xtitle, const char* ytitle ){
    TCanvas *c = new TCanvas();
    h1->Draw(doption1);
    h1->SetLineWidth(width1);
    h1->SetLineColor(color1);
    h1->GetXaxis()->SetTitle(xtitle);
    h1->GetYaxis()->SetTitle(ytitle);
    h2->Draw(doption2);
    h2->SetLineWidth(width2);
    h2->SetLineColor(color2);

    TLegend *l = new TLegend();
    l->AddEntry(h1,title1);
    l->AddEntry(h2,title2);
    l->Draw("same");

    return c;
}

TCanvas* canvasGraph( TGraphErrors *g1, TGraphErrors *g2, const char* title1, const char* title2, const char* doption1, const char* doption2, Color_t color1, Color_t color2, int style1, int style2, const char* xtitle, const char* ytitle, TF1* func ){

    TCanvas *c = new TCanvas();
    g1->Draw(doption1);
    g1->SetMarkerStyle(style1);
    g1->SetLineColor(color1);
    g1->SetMarkerColor(color1);
    g1->GetXaxis()->SetTitle(xtitle);
    g1->GetYaxis()->SetTitle(ytitle);
    g2->Draw(doption2);
    g2->SetMarkerStyle(style2);
    g2->SetMarkerColor(color2);
    g2->SetLineColor(color2);

    func->SetLineColor(kGray);
    func->SetLineStyle(kDashed);
    func->Draw("same");

    TLegend *l = new TLegend();
    l->AddEntry(g1,title1);
    l->AddEntry(g2,title2);
    l->Draw("same");

    return c;
}

std::array<double,2> scartofunc(double x, double y, double sigmax, double sigmay)
{
    std::array<double,2> returned;
    returned[0]=(x-y)/x*100.;
    returned[1]=100.*std::sqrt(std::pow((x-(x-y))*sigmax/x/x,2)+std::pow(sigmay/x,2));
    return returned;
}

void ratioDatiMC()
{
    //ofstream scarto("scarto.txt");
    TH1D *dEdx_dati = new TH1D("dEdx_dati", "", 500, 0, 2500);
    TH2D *checkfactor = new TH2D("check_factor", "", 1000, 0, 15, 1000, -2,2);

    TF1* fscale = new TF1("fscale", Correction_muon, 4, 0,16);
    fscale->SetParameters(7.65506e-03, 1.59985e+00, -2.06713e-02, 5.21996e-02);

    EventsData dati = load_data("dati", "muon", "Np", "full");
    for(int track=0; track<dati.tree->GetEntries(); track++)
    {
        dati.tree->GetEntry(track);
        if(dati.track.len_reco<150.)continue;
        for(int hit=0; hit<dati.track.dQdx->size(); hit++)
        {
            if(dati.track.rr->at(hit)>=100. && dati.track.rr->at(hit)<=150.)
            {
                dEdx_dati->Fill(dati.track.dQdx->at(hit));
            }
        }
    }

    dEdx_dati->Scale(1./dEdx_dati->Integral());
    double x_low_dati, x_high_dati, x_low_dati_err, x_high_dati_err;
    double peak_data = dEdx_dati->GetBinCenter(FindBinMax(dEdx_dati));
    double peak_data_error = dEdx_dati->GetBinWidth(FindBinMax(dEdx_dati))/2.;
    x_low_dati = FindXAtY(dEdx_dati, dEdx_dati->GetBinContent(FindBinMax(dEdx_dati))/2., dEdx_dati->GetBinCenter(FindBinMax(dEdx_dati)), false );
    x_high_dati = FindXAtY(dEdx_dati, dEdx_dati->GetBinContent(FindBinMax(dEdx_dati))/2., dEdx_dati->GetBinCenter(FindBinMax(dEdx_dati)), true );
    x_low_dati_err = dEdx_dati->GetBinWidth(2)/2.;
    x_high_dati_err = dEdx_dati->GetBinWidth(2)/2.;
    cout << dEdx_dati->GetBinCenter(FindBinMax(dEdx_dati)) << " " << x_low_dati << " " << x_high_dati << endl;


    TF1 *line1 = new TF1("line1", "1", 0.,5000);
    TCanvas *dati_vs_mc = new TCanvas();
    dati_vs_mc->Divide(0,2);
    dati_vs_mc->cd(1);
    dEdx_dati->Draw("hist");
    dEdx_dati->SetLineWidth(2);
    dEdx_dati->SetLineColor(kBlack);
    TLegend *l = new TLegend();
    l->AddEntry(dEdx_dati, "dati");

    dati_vs_mc->cd(2);
    line1->Draw();
    line1->SetLineColor(kGray);
    line1->SetLineStyle(kDashed);

    TLegend *l_canvas_diff = new TLegend();
    TF1 *line0 = new TF1("line0", "0", 0.99, 1.07);
    TCanvas *canvas_difference = new TCanvas();
    canvas_difference->Divide(0,3);
    canvas_difference->cd(1);
    /*line0->Draw();
    line0->SetLineColor(kBlack);
    line0->SetLineStyle(kDashed);
    line0->SetLineWidth(2);
    line0->GetYaxis()->SetRangeUser(-50,50);
    line0->GetXaxis()->SetTitle("scaling factor");
    line0->GetYaxis()->SetTitle("peak position difference (%)");
    canvas_difference->cd(2);
    line0->Draw();
    line0->SetLineColor(kBlack);
    line0->SetLineStyle(kDashed);
    line0->SetLineWidth(2);
    line0->GetYaxis()->SetRangeUser(-50,50);
    line0->GetXaxis()->SetTitle("scaling factor");
    line0->GetYaxis()->SetTitle("x left difference (%)");
    canvas_difference->cd(3);
    line0->Draw();
    line0->SetLineColor(kBlack);
    line0->SetLineStyle(kDashed);
    line0->SetLineWidth(2);
    line0->GetYaxis()->SetRangeUser(-50,50);
    line0->GetXaxis()->SetTitle("scaling factor");
    line0->GetYaxis()->SetTitle("x right difference (%)");
    */

    int color=0;
    std::array<Color_t,11> colors = {kCyan, kRed, kAzure, kOrange, kGreen, kMagenta, kYellow, kGreen+3, kGray, kRed-9, kViolet+1};
    for(double scaling=1.00; scaling<1.33; scaling+=0.03)
    {
        TH1D *dEdx_mc = new TH1D(Form("dEdx_mc_scaled_%.2f",scaling), "", 500, 0, 2500);
        TGraphErrors *peak_diff = new TGraphErrors();
        TGraphErrors *low_diff = new TGraphErrors();
        TGraphErrors *high_diff = new TGraphErrors();
        EventsData mc = load_data("mc", "muon", "Np", "full");
        for(int track=0; track<mc.tree->GetEntries(); track++)
        {
            mc.tree->GetEntry(track);
            if(mc.track.len_reco<150.)continue;
            for(int hit=0; hit<mc.track.dQdx->size(); hit++)
            {
                if(mc.track.rr->at(hit)>=100. && mc.track.rr->at(hit)<=150.)
                {
                    //double factor = (2+fscale->Eval(mc.track.dE->at(hit)))/(2-fscale->Eval(mc.track.dE->at(hit)));
                    dEdx_mc->Fill(mc.track.dQdx->at(hit)*scaling);
        
                    //checkfactor->Fill(mc.track.dE->at(hit),factor);
                }
            }
        }
    


    
    dEdx_mc->Scale(1./dEdx_mc->Integral());

    //TGraphErrors *ratio = new TGraphErrors();
    //ratio = histoRatio(dEdx_dati,dEdx_mc);

    double x_low_mc, x_high_mc, x_low_mc_err, x_high_mc_err;
    x_low_mc = FindXAtY(dEdx_mc, dEdx_mc->GetBinContent(FindBinMax(dEdx_mc))/2., dEdx_mc->GetBinCenter(FindBinMax(dEdx_mc)), false );
    x_high_mc = FindXAtY(dEdx_mc, dEdx_mc->GetBinContent(FindBinMax(dEdx_mc))/2., dEdx_mc->GetBinCenter(FindBinMax(dEdx_mc)), true );
    x_low_mc_err = dEdx_mc->GetBinWidth(1)/2.;
    x_high_mc_err = dEdx_mc->GetBinWidth(1)/2.;
    
    double peak_mc = dEdx_mc->GetBinCenter(FindBinMax(dEdx_mc));
    double peak_mc_error = dEdx_mc->GetBinWidth(FindBinMax(dEdx_mc))/2.;
    cout << dEdx_mc->GetBinCenter(FindBinMax(dEdx_mc)) << " " << x_low_mc << " " << x_high_mc <<  " " << scartofunc(peak_data,peak_mc,peak_data_error,peak_mc_error)[0] << " " << (x_low_dati-x_low_mc)/x_low_dati*100 << " " << (x_high_dati-x_high_mc)/x_high_dati*100 << endl;

    //scarto << scartofunc(peak_data,peak_mc,peak_data_error,peak_mc_error)[0] << " " << scartofunc(peak_data,peak_mc,peak_data_error,peak_mc_error)[1] << " ";
    //scarto << scartofunc(x_low_dati,x_low_mc,x_low_dati_err,x_low_mc_err)[0] << " " << scartofunc(x_low_dati,x_low_mc,x_low_dati_err,x_low_mc_err)[1] << " ";
    //scarto << scartofunc(x_high_dati,x_high_mc,x_high_dati_err,x_high_mc_err)[0] << " " << scartofunc(x_high_dati,x_high_mc,x_high_dati_err,x_high_mc_err)[1] << endl;

    peak_diff->SetPoint(0,scaling,scartofunc(peak_data,peak_mc,peak_data_error,peak_mc_error)[0]);
    peak_diff->SetPointError(0,0,scartofunc(peak_data,peak_mc,peak_data_error,peak_mc_error)[1]);
    low_diff->SetPoint(0,scaling,scartofunc(x_low_dati,x_low_mc,x_low_dati_err,x_low_mc_err)[0]);
    low_diff->SetPointError(0,0,scartofunc(x_low_dati,x_low_mc,x_low_dati_err,x_low_mc_err)[1]);
    high_diff->SetPoint(0,scaling,scartofunc(x_high_dati,x_high_mc,x_high_dati_err,x_high_mc_err)[0]);
    high_diff->SetPointError(0,0,scartofunc(x_high_dati,x_high_mc,x_high_dati_err,x_high_mc_err)[1]);

    peak_diff->SetMarkerStyle(7);
    peak_diff->SetMarkerColor(colors[color]);
    peak_diff->SetLineColor(colors[color]);
    low_diff->SetMarkerStyle(7);
    low_diff->SetMarkerColor(colors[color]);
    low_diff->SetLineColor(colors[color]);
    high_diff->SetMarkerStyle(7);
    high_diff->SetMarkerColor(colors[color]);
    high_diff->SetLineColor(colors[color]);

    canvas_difference->cd(1);
    peak_diff->Draw("P same");
    l_canvas_diff->AddEntry(peak_diff,Form("mc x %.2f", scaling));
    canvas_difference->cd(2);
    low_diff->Draw("P same");
    canvas_difference->cd(3);
    high_diff->Draw("P same");  
      

    dati_vs_mc->cd(1);
    dEdx_mc->Draw("hist same");
    dEdx_mc->SetLineWidth(2);
    dEdx_mc->SetLineColor(colors[color]);
    
    l->AddEntry(dEdx_mc,Form("mc x %.2f", scaling));

    
    TGraphErrors *ratio = new TGraphErrors();
    ratio = histoRatio(dEdx_dati,dEdx_mc);
    dati_vs_mc->cd(2);
    ratio->Draw("P same");
    ratio->SetMarkerStyle(7);
    ratio->SetMarkerColor(colors[color]);
    ratio->SetLineColor(colors[color]);
    color++;

    }
    dati_vs_mc->cd(1);
    l->Draw("same");
    canvas_difference->cd(1);
    l_canvas_diff->Draw("same");
    TFile *f = new TFile("plot_output.root", "RECREATE");
    dati_vs_mc->Write("dati_vs_mc");
    canvas_difference->Write("canvas_difference");

}
