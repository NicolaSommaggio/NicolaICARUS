#include "ReadTree.C"

void confrontoDatiMClight(std::string particle, std::string configration, const char *filename)
{
    EventsData dat = load_data(configration , particle);
    
    cout << dat.tree->GetEntries() << " " << particle << " tracks" << endl;

    TFile *f = TFile::Open(filename , "UPDATE");
    TDirectory *d= (TDirectory*)f->Get(particle.c_str());
    TDirectory *d1= (TDirectory*)d->Get(dati_o_mc.c_str()); 

    d1->cd();

    std::vector<std::vector<double>> dE;
    std::vector<std::vector<double>> rr;  

    int Nbins_rr=100;
    double low_rr=0;
    double high_rr=30;
    

    TH2D *dedx_range= new TH2D("dedx_range", " ", Nbins_rr, low_rr, high_rr, 300, 0, 30);
    TH1D *h_median = new TH1D("h_median", "", 200, 0, 20 );

    TH1D *hdirx = new TH1D("dirx", "", 200, -1, 1);
    TH1D *hdiry = new TH1D("diry", "", 200, -1, 1);
    TH1D *hdirz = new TH1D("dirz", "", 200, -1, 1);

    TH1D *h_theta_drift = new TH1D("htheta_drif", "", 180, 0, 90 );

    //looking at rr between 25 and 30 cm
    TH1D *hdEdx = new TH1D("hdEdx", "", 300, 0., 30.);
    TH1D *hdEdx_onlymult1 = new TH1D("hdEdx_onlymult1", "", 300, 0, 30);
    TH1D *hdEdx_multmag1 = new TH1D("hdEdx_multmag1", "", 300, 0, 30);
    TH1D *hdQdx = new TH1D("hdQdx", "", 500, 0., 2000.);
    TH1D *hdQdx_onlymult1 = new TH1D("hdQdx_onlymult1", "", 500, 0., 2000.);
    TH1D *hdQdx_multmag1 = new TH1D("hdQdx_multmag1", "", 500, 0., 2000.);
    TH1D *hdEdx_ind1 = new TH1D("hdEdx_ind1", "", 300, 0., 30.);
    TH1D *hdEdx_ind2 = new TH1D("hdEdx_ind2", "", 300, 0., 30.);
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
                }
            }
        }
        
        //looking also to ind1 and 1nd2
        if(dat.track.len_reco>=30)
        {
            for(int hit=0; hit<int(dat.track.rr_ind1->size()); hit++)
            {
                if(dat.track.rr_ind1->at(hit)>=25. && dat.track.rr_ind1->at(hit)<30.){hdEdx_ind1->Fill(dat.track.dEdx_ind1->at(hit));}
            }

            for(int hit=0; hit<int(dat.track.rr_ind2->size()); hit++)
            {
                if(dat.track.rr_ind2->at(hit)>=25. && dat.track.rr_ind2->at(hit)<30.){hdEdx_ind2->Fill(dat.track.dEdx_ind2->at(hit));}
            }
        }

        //plot banana
        for(int hit=0; hit<int(dat.track.rr_ind1->size()); hit++){dEdx_range_ind1->Fill(dat.track.rr_ind1->at(hit), dat.track.dEdx_ind1->at(hit));}
        for(int hit=0; hit<int(dat.track.rr_ind2->size()); hit++){dEdx_range_ind2->Fill(dat.track.rr_ind2->at(hit), dat.track.dEdx_ind2->at(hit));}


        // median and filter on the median
        for(int hit=0; hit<int(dat.track.rr->size()); hit++)
        {
            dedx_range->Fill(dat.track.rr->at(hit), dat.track.dE->at(hit));

            rr_temp.push_back(dat.track.rr->at(hit));
            dEdx_temp.push_back(dat.track.dE->at(hit));

            if(dat.track.rr->at(hit)<5.){trackdEdx.push_back(dat.track.dE->at(hit));}
        }

        if (trackdEdx.size() >= 2) {trackdEdx.erase(trackdEdx.end() - 2, trackdEdx.end());} //removing the last 2 hits to compute the median
    
        double med =0;
        if(trackdEdx.size()>0) { med= mediana(trackdEdx); h_median->Fill(med);}

        double threshold=0;
        if(particle=="muon" && configuration=="mc1D"){threshold=3.2;}
        else if(particle=="muon" && (configuration=="mc2D" || configuration=="mc2D_DNN_PT1" || configuration=="mc2D_DNN_PT10")){threshold=3.5;}
        else if(particle=="proton" && configuration=="mc1D"){threshold=8.;}
        else if(particle=="muon" && (configuration=="mc2D" || configuration=="mc2D_DNN_PT1" || configuration=="mc2D_DNN_PT10")){threshold=8.;}
        //still to look for pions

        if(med>=threshold)
        {
            dE.push_back(dEdx_temp);
            rr.push_back(rr_temp);
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

    cout << particle << " " << dati_o_mc << " " << rr.size() << " tracks with a NON MIP-like median" << endl;
    cout << "****************************************************" << endl;
    //dEdx distributions in each rr bin
    cout << 
    for(double rr=1.5; rr<=25; rr+=0.5)
    {
        cout << "processing rr " << rr << endl;
        TDirectory *rrdir = (TDirectory*)d1->mkdir(Form("%.1f", rr));
        rrdir->cd();

        int Nbins_coll = 300;
        int Nbins_ind1 = 300;
        int Nbins_ind2 = 300;

        std::vector<double> binlowedges_coll;
        std::vector<double> binlowedges_ind1;
        std::vector<double> binlowedges_ind2;
        for(int i=1; i<=Nbins_coll+1; i++)
        {
            binlowedges_coll.push_back((i-1)*30./Nbins_coll);
            binlowedges_ind1.push_back((i-1)*30./Nbins_ind1);
            binlowedges_ind2.push_back((i-1)*30./Nbins_ind2);
        }

        TH1D *dEdx_coll;
        TH1D *dEdx_ind1;
        TH1D *dEdx_ind2;

        dEdx_coll = new TH1D(Form("dEdx_coll_rr_%.1f", rr),"", Nbins_coll, binlowedges_coll.data());
        dEdx_ind1 = new TH1D(Form("dEdx_ind1_rr_%.1f", rr),"", Nbins_ind1, binlowedges_ind1.data());
        dEdx_ind2 = new TH1D(Form("dEdx_ind2_rr_%.1f", rr),"", Nbins_ind2, binlowedges_ind2.data());


        for(int track=0; track<dat.tree->GetEntries(); track++)
        {
            dat.tree->GetEntry(track);

            //collection plane
            for(int hit=0; hit<int(dat.track.dE->size()); hit++)
            {
                if(dat.track.rr->at(hit)<rr && dat.track.rr->at(hit)>=rr-0.5){dEdx_coll->Fill(dat.track.dE->at(hit));}   
            }   

            //induction 1 plane
            for(int hit=0; hit<int(dat.track.dEdx_ind1->size()); hit++)
            {
                if(dat.track.rr_ind1->at(hit)<rr && dat.track.rr_ind1->at(hit)>=rr-0.5){dEdx_ind1->Fill(dat.track.dEdx_ind1->at(hit));}   
            } 

            //induction 2 plane
            for(int hit=0; hit<int(dat.track.dEdx_ind2->size()); hit++)
            {
                if(dat.track.rr_ind2->at(hit)<rr && dat.track.rr_ind2->at(hit)>=rr-0.5){dEdx_ind2->Fill(dat.track.dEdx_ind2->at(hit));}   
            }
        }   

        dEdx_coll->Scale(1./dEdx_coll->Integral());
        dEdx_coll->Write(0,TObject::kOverwrite);
        dEdx_ind1->Scale(1./dEdx_ind1->Integral());
        dEdx_ind1->Write(0,TObject::kOverwrite);
        dEdx_ind2->Scale(1./dEdx_ind2->Integral());
        dEdx_ind2->Write(0,TObject::kOverwrite);  

    }

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


void dEdxdistro()
{
    //gROOT->ProcessLine(".L libdaughtersInfo.so");
    std::array<std::string,3> particles={"muon", "proton", "pion"};
    TFile *f = TFile::Open("dEdx_distro_general.root", "UPDATE");
    for(const auto& par : particles)
    {
        EventsData dat = load_data("mc2d_general", par);

        cout << dat.tree->GetEntries() << " " << par << "tracks" << endl;
        TDirectory *pardir = (TDirectory*)f->Get(Form("%s",par.c_str()));
        
        for(double rr=1.5; rr<=25; rr+=0.5)
        {
            cout << "processing rr " << rr << endl;
            TDirectory *rrdir = (TDirectory*)pardir->mkdir(Form("%.1f", rr));
            rrdir->cd();

            int Nbins_coll = 300;
            int Nbins_ind1 = 300;
            int Nbins_ind2 = 300;

            std::vector<double> binlowedges_coll;
            std::vector<double> binlowedges_ind1;
            std::vector<double> binlowedges_ind2;
            for(int i=1; i<=Nbins_coll+1; i++)
            {
                binlowedges_coll.push_back((i-1)*30./Nbins_coll);
                binlowedges_ind1.push_back((i-1)*30./Nbins_ind1);
                binlowedges_ind2.push_back((i-1)*30./Nbins_ind2);
            }

            TH1D *dEdx_coll;
            TH1D *dEdx_ind1;
            TH1D *dEdx_ind2;

            dEdx_coll = new TH1D(Form("dEdx_coll_rr_%.1f", rr),"", Nbins_coll, binlowedges_coll.data());
            dEdx_ind1 = new TH1D(Form("dEdx_ind1_rr_%.1f", rr),"", Nbins_ind1, binlowedges_ind1.data());
            dEdx_ind2 = new TH1D(Form("dEdx_ind2_rr_%.1f", rr),"", Nbins_ind2, binlowedges_ind2.data());


            for(int track=0; track<dat.tree->GetEntries(); track++)
            {
                dat.tree->GetEntry(track);

                //collection plane
                for(int hit=0; hit<int(dat.track.dE->size()); hit++)
                {
                    if(dat.track.rr->at(hit)<rr && dat.track.rr->at(hit)>=rr-0.5)
                    {
                        dEdx_coll->Fill(dat.track.dE->at(hit));
                    }   
                }   

                //induction 1 plane
                for(int hit=0; hit<int(dat.track.dEdx_ind1->size()); hit++)
                {
                    if(dat.track.rr_ind1->at(hit)<rr && dat.track.rr_ind1->at(hit)>=rr-0.5)
                    {
                        dEdx_ind1->Fill(dat.track.dEdx_ind1->at(hit));
                    }   
                } 

                //induction 2 plane
                for(int hit=0; hit<int(dat.track.dEdx_ind2->size()); hit++)
                {
                    if(dat.track.rr_ind2->at(hit)<rr && dat.track.rr_ind2->at(hit)>=rr-0.5)
                    {
                        dEdx_ind2->Fill(dat.track.dEdx_ind2->at(hit));
                    }   
                }
            }   

            dEdx_coll->Scale(1./dEdx_coll->Integral());
            dEdx_coll->Write(0,TObject::kOverwrite);

            dEdx_ind1->Scale(1./dEdx_ind1->Integral());
            dEdx_ind1->Write(0,TObject::kOverwrite);

            dEdx_ind2->Scale(1./dEdx_ind2->Integral());
            dEdx_ind2->Write(0,TObject::kOverwrite);  

        }
    }
}