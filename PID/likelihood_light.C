#include "../dEdx/fitLangaulight.C"
#include "../dEdx/ReadTree.C"


std::vector<double> chi2_ALG(std::vector<double> &dEdx,std::vector<double> &RR, double rr_min, double rr_max)
{
    TFile* file = TFile::Open("/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/dEdx/THdedx.root");
    auto dedx_range_pro = (TProfile*)file->Get("dedx_range_pro");
    auto dedx_range_ka  = (TProfile*)file->Get("dedx_range_ka");
    auto dedx_range_pi  = (TProfile*)file->Get("dedx_range_pi");
    auto dedx_range_mu  = (TProfile*)file->Get("dedx_range_mu");

    double threshold = 0.5;
    double max_rr = rr_max;
    double min_rr = rr_min;

    std::vector<float> trkdedx;
    std::vector<float> trkres;
    std::vector<double> vpida;
    

    for(std::size_t i(0); i<dEdx.size(); ++i){
      if(i==0 || i==dEdx.size()-1)continue;
        if(RR[i]<max_rr && RR[i]>rr_min ){trkdedx.push_back(dEdx[i]);trkres.push_back(RR[i]);}       
    }


    int npt = 0;
    double chi2pro = 0;
    double chi2ka = 0;
    double chi2pi = 0;
    double chi2mu = 0;
    double avgdedx = 0;
    double PIDA = 0;

    int used_trkres = 0;
    for (unsigned i = 0; i < trkdedx.size(); ++i) { 
      avgdedx += trkdedx[i];
      if (trkres[i] < 26) {
        PIDA += trkdedx[i] * pow(trkres[i], 0.42);
        vpida.push_back(trkdedx[i] * pow(trkres[i], 0.42));
        used_trkres++;
      }
      if (trkdedx[i] > 100 || trkdedx[i]<threshold) continue; //protect against large pulse height
      int bin = dedx_range_pro->FindBin(trkres[i]);
      if (bin >= 1 && bin <= dedx_range_pro->GetNbinsX()) {
        double bincpro = dedx_range_pro->GetBinContent(bin);
        if (bincpro < 1e-6) { //for 0 bin content, using neighboring bins
          bincpro =
            (dedx_range_pro->GetBinContent(bin - 1) + dedx_range_pro->GetBinContent(bin + 1)) / 2;
        }
        double bincka = dedx_range_ka->GetBinContent(bin);
        if (bincka < 1e-6) {
          bincka =
            (dedx_range_ka->GetBinContent(bin - 1) + dedx_range_ka->GetBinContent(bin + 1)) / 2;
        }
        double bincpi = dedx_range_pi->GetBinContent(bin);
        if (bincpi < 1e-6) {
          bincpi =
            (dedx_range_pi->GetBinContent(bin - 1) + dedx_range_pi->GetBinContent(bin + 1)) / 2;
        }
        double bincmu = dedx_range_mu->GetBinContent(bin);
        if (bincmu < 1e-6) {
          bincmu =
            (dedx_range_mu->GetBinContent(bin - 1) + dedx_range_mu->GetBinContent(bin + 1)) / 2;
        }
        double binepro = dedx_range_pro->GetBinError(bin);
        if (binepro < 1e-6) {
          binepro =
            (dedx_range_pro->GetBinError(bin - 1) + dedx_range_pro->GetBinError(bin + 1)) / 2;
        }
        double bineka = dedx_range_ka->GetBinError(bin);
        if (bineka < 1e-6) {
          bineka = (dedx_range_ka->GetBinError(bin - 1) + dedx_range_ka->GetBinError(bin + 1)) / 2;
        }
        double binepi = dedx_range_pi->GetBinError(bin);
        if (binepi < 1e-6) {
          binepi = (dedx_range_pi->GetBinError(bin - 1) + dedx_range_pi->GetBinError(bin + 1)) / 2;
        }
        double binemu = dedx_range_mu->GetBinError(bin);
        if (binemu < 1e-6) {
          binemu = (dedx_range_mu->GetBinError(bin - 1) + dedx_range_mu->GetBinError(bin + 1)) / 2;
        }
        //double errke = 0.05*trkdedx[i];   //5% KE resolution
        double errdedx = 0.04231 + 0.0001783 * trkdedx[i] * trkdedx[i]; //resolution on dE/dx
        errdedx *= trkdedx[i];
        chi2pro += pow((trkdedx[i] - bincpro) / std::sqrt(pow(binepro, 2) + pow(errdedx, 2)), 2);
        chi2ka += pow((trkdedx[i] - bincka) / std::sqrt(pow(bineka, 2) + pow(errdedx, 2)), 2);
        chi2pi += pow((trkdedx[i] - bincpi) / std::sqrt(pow(binepi, 2) + pow(errdedx, 2)), 2);
        chi2mu += pow((trkdedx[i] - bincmu) / std::sqrt(pow(binemu, 2) + pow(errdedx, 2)), 2);
        //std::cout<<i<<" "<<trkdedx[i]<<" "<<trkres[i]<<" "<<bincpro<<std::endl;
        ++npt;
      }
    } //hits
        std::vector<double> chi2s {chi2mu/npt,chi2pro/npt,chi2ka/npt,chi2pi/npt};

    file->Close();
    delete file;
    return chi2s;
}

bool has_valid_hits(std::vector<double> rr_temp, std::vector<double> dEdx_temp)
{
    bool valid=false;
    if(rr_temp.size()<3)return false;
    for(int i=1; i<rr_temp.size()-1; i++)
    {
        if(dEdx_temp[i]<=100 && dEdx_temp[i]>=0.5) return true;
    }
    return false;
}


int find_idx(double rr)
{
    if(rr==25.)return 47;
    double a = std::round(rr);
    double max=-1;
    double min=-1;
    if(a>rr)
    {
        min=a-0.5;
        max=a;
    }
    else
    {
        min=a;
        max=a+0.5;
    }
    
    int n = (max-1.5)/0.5;

    return n;
}


TGraph* ROCbuilder(std::vector<std::pair<double,double>> pts, const char* name)
{
    TGraph * rocCurve = new TGraph();
    rocCurve->SetName(name);
    std::sort(pts.begin(), pts.end(),[](const std::pair<double,double>& a, const std::pair<double,double>& b){ return a.first < b.first; });

    std::vector<std::pair<double,double>> best;
    double current_x = pts[0].first;
    std::vector<double> temp_sig_eff;
    temp_sig_eff.push_back(pts[0].second);

    for(int i=1; i<pts.size(); i++)
    {
        if(pts[i].first==current_x) temp_sig_eff.push_back(pts[i].second);
        else
        {
            best.push_back({current_x , *std::max_element(temp_sig_eff.begin() , temp_sig_eff.end())});
            temp_sig_eff.clear();
            current_x=pts[i].first;
            temp_sig_eff.push_back(pts[i].second);
        }
    }
    best.push_back({current_x, *std::max_element(temp_sig_eff.begin(), temp_sig_eff.end())});
    
    for(int p=0; p<best.size(); p++) rocCurve->SetPoint(p, best[p].first, best[p].second);

    rocCurve->SetPoint(rocCurve->GetN(),1,0);
    
    return rocCurve;

}


//ofstream controlloROC("controlloROC.txt");
std::array<TGraph*,4> GetROCcurveChi2(std::string name, int N1, int N2,  double T_min, double T_max, const std::array<std::vector<double>,2>& v_signal, const std::array<std::vector<double>,2>& v_background1, const std::array<std::vector<double>,2>& v_background2)
{
    int NT1 = N1;
    int NT2 = N2;
    double Tmin = T_min;
    double Tmax = T_max;

    std::vector<std::pair<double,double>> pts_mu_pro;
    std::vector<std::pair<double,double>> pts_mu_pi;
    std::vector<std::pair<double,double>> pts_pro_pi;
    std::vector<std::pair<double,double>> pts_mu_other;

    TGraph *g_mu_pro = new TGraph();
    g_mu_pro->SetName(Form("new_mu_pro_%s",name.c_str()));
    TGraph *g_mu_pi = new TGraph();
    g_mu_pi->SetName(Form("new_mu_pi_%s",name.c_str()));
    TGraph *g_pro_pi = new TGraph();
    g_pro_pi->SetName(Form("new_pro_pi_%s",name.c_str()));
    TGraph *g_mu_other = new TGraph();
    g_mu_other->SetName(Form("new_mu_other_%s",name.c_str()));

    for (int i = 0; i <= NT1; i++)
    {
        for(int j=0; j<=NT2; j++)
        {   
            double thr = Tmin + i * (Tmax - Tmin) / (NT1);
            double thr1 = Tmin + j* (Tmax - Tmin) / (NT2);

            double count_s = 0;
            double count_b1 = 0;
            double count_b2 = 0;
            double tot_s=v_signal[0].size();
            double tot_b1=v_background1[0].size();
            double tot_b2=v_background2[0].size();

            for(int entry=0; entry<v_signal[0].size(); entry++) if (v_signal[0][entry] < thr && v_signal[1][entry] > thr1){count_s++;} 

            for(int entry=0; entry<v_background1[0].size(); entry++) if (v_background1[0][entry] < thr && v_background1[1][entry] > thr1){count_b1++;} 
            
            for(int entry=0; entry<v_background2[0].size(); entry++) if (v_background2[0][entry] < thr && v_background2[1][entry] > thr1){count_b2++;}
            
            pts_mu_pro.push_back({1-count_b1/tot_b1,count_s/tot_s});
            pts_mu_pi.push_back({1-count_b2/tot_b2,count_s/tot_s});
            pts_pro_pi.push_back({1-count_b2/tot_b2,count_b1/tot_b1});
            pts_mu_other.push_back({1-(count_b1+count_b2)/(tot_b1+tot_b2),count_s/tot_s});
            
            int point = g_mu_pro->GetN();
            g_mu_pro->SetPoint(point,1-count_b1/tot_b1,count_s/tot_s);
            g_mu_pi->SetPoint(point,1-count_b2/tot_b2,count_s/tot_s);
            g_pro_pi->SetPoint(point,1-count_b2/tot_b2,count_b1/tot_b1);
            g_mu_other->SetPoint(point,1-(count_b1+count_b2)/(tot_b1+tot_b2),count_s/tot_s);
        }
    }


    //per ogni valore di bkg rejection prendo la miglior signal eff
    //std::array<TGraph*,4> rocCurves = {ROCbuilder(pts_mu_pro,Form("new_mu_pro_%s",name.c_str())), ROCbuilder(pts_mu_pi,Form("new_mu_pi_%s",name.c_str())), ROCbuilder(pts_pro_pi,Form("new_pro_pi_%s",name.c_str())), ROCbuilder(pts_mu_other,Form("new_mu_other_%s",name.c_str()))};

    //return rocCurves;

    return {g_mu_pro,g_mu_pi,g_pro_pi,g_mu_other};

}



std::array<TGraph*,4> GetROCcurveChi2new(TFile* f, std::string name, const std::array<std::vector<double>,2>& v_signal, const std::array<std::vector<double>,2>& v_background1, const std::array<std::vector<double>,2>& v_background2)
{
    int nbinsx = 200;
    int nbinsy = 500;
    TH2D * hBKG1 = new TH2D(Form("hBKG1_%s",name.c_str()), "", nbinsx,0,200,nbinsy,0,500);
    TH2D * hBKG2 = new TH2D(Form("hBKG2_%s",name.c_str()), "", nbinsx,0,200,nbinsy,0,500);
    TH2D * hBKGtot = new TH2D(Form("hBKGtot_%s",name.c_str()), "", nbinsx,0,200,nbinsy,0,500);
    TH2D * hSIG = new TH2D(Form("hSIG_%s",name.c_str()), "", nbinsx,0,200,nbinsy,0,500);

    for(int i=0; i<v_signal[0].size(); i++){hSIG->Fill(v_signal[0][i], v_signal[1][i]);}
    for(int i=0; i<v_background1[0].size(); i++){hBKG1->Fill(v_background1[0][i], v_background1[1][i]); hBKGtot->Fill(v_background1[0][i], v_background1[1][i]);}
    for(int i=0; i<v_background2[0].size(); i++){hBKG2->Fill(v_background2[0][i], v_background2[1][i]); hBKGtot->Fill(v_background2[0][i], v_background2[1][i]);}

    f->cd();
    TDirectory *d = (TDirectory*)f->Get("chi2");
    d->cd();
    hBKG1->Write(0,TObject::kOverwrite);
    hBKG2->Write(0,TObject::kOverwrite);
    hBKGtot->Write(0,TObject::kOverwrite);
    hSIG->Write(0,TObject::kOverwrite);

    double tot_signal = hSIG->Integral();
    double tot_bkg1 = hBKG1->Integral();
    double tot_bkg2 = hBKG2->Integral();
    double tot_bkgtot = hBKGtot->Integral();

    //TGraph *g_mu_pro = new TGraph();
    //g_mu_pro->SetName(Form("mu_pro_%s",name.c_str()));
    //TGraph *g_mu_pi = new TGraph();
    //g_mu_pi->SetName(Form("mu_pi_%s",name.c_str()));
    //TGraph *g_pro_pi = new TGraph();
    //g_pro_pi->SetName(Form("pro_pi%s",name.c_str()));
    //TGraph *g_mu_other = new TGraph();
    //g_mu_other->SetName(Form("mu_other_%s",name.c_str()));

    std::vector<std::pair<double,double>> pts_mu_pro;
    std::vector<std::pair<double,double>> pts_mu_pi;
    std::vector<std::pair<double,double>> pts_pro_pi;
    std::vector<std::pair<double,double>> pts_mu_other;
    

    for(int ix=1; ix<nbinsx; ix++)
    {
        for(int iy=1; iy<nbinsy; iy++)
        {
            double counts_signal = hSIG->Integral(1,ix,iy,nbinsy);
            double counts_bkg1 = hBKG1->Integral(1,ix,iy,nbinsy);
            double counts_bkg2 = hBKG2->Integral(1,ix,iy,nbinsy);
            double counts_bkgtot = hBKGtot->Integral(1,ix,iy,nbinsy);

            //int point = g_mu_pro->GetN();
            //g_mu_pro->SetPoint(point,1-counts_bkg1/tot_bkg1,counts_signal/tot_signal);
            //g_mu_pi->SetPoint(point,1-counts_bkg2/tot_bkg2,counts_signal/tot_signal);
            //g_pro_pi->SetPoint(point,1-counts_bkg2/tot_bkg2,counts_bkg1/tot_bkg1);
            //g_mu_other->SetPoint(point,1-counts_bkgtot/tot_bkgtot,counts_signal/tot_signal);

            pts_mu_pro.push_back({1-counts_bkg1/tot_bkg1,counts_signal/tot_signal});
            pts_mu_pi.push_back({1-counts_bkg2/tot_bkg2,counts_signal/tot_signal});
            pts_pro_pi.push_back({1-counts_bkg2/tot_bkg2,counts_bkg1/tot_bkg1});
            pts_mu_other.push_back({1-counts_bkgtot/tot_bkgtot,counts_signal/tot_signal});
        }
    }

    std::array<TGraph*,4> rocCurves = {ROCbuilder(pts_mu_pro,Form("mu_pro_%s",name.c_str())), ROCbuilder(pts_mu_pi,Form("mu_pi_%s",name.c_str())), ROCbuilder(pts_pro_pi,Form("pro_pi_%s",name.c_str())), ROCbuilder(pts_mu_other,Form("mu_other_%s",name.c_str()))};

    return rocCurves;
    //return {g_mu_pro,g_mu_pi,g_pro_pi,g_mu_other};

}


std::array<TGraph*,4> getROCcurvesLR(int N, const std::vector<std::vector<double>> &LRmu, const std::vector<std::vector<double>> &LRpro, const std::vector<std::vector<double>> &LRpi, std::string plane)
{
    //each entry of LRmu, LRpro and LRpi is a vector of three element containing the likelihood ratios muon/proton, muon/pion, proton/pion
    
    double Tmin=-1;
    double Tmax=1;

    std::vector<std::pair<double,double>> pts_mu_pro;
    std::vector<std::pair<double,double>> pts_mu_pi;
    std::vector<std::pair<double,double>> pts_pro_pi;
    std::vector<std::pair<double,double>> pts_mu_other;

    for(int t=0; t<=N; t++)
    {
        double thr = Tmin + t * (Tmax - Tmin) / (N);
        int count_s=0; 
        int count_b=0;
        int count_b2=0;
        double eff_b=0;
        double eff_s=0;
        double eff_tot=0;

        for(const auto &v : LRmu) if(v[0]<thr) count_s++;
        eff_s = double(count_s)/LRmu.size();
        for(auto v : LRpro) if (v[0]<thr) count_b++;
        eff_b = double(count_b)/LRpro.size();
        pts_mu_pro.push_back({1 - eff_b,eff_s});
        //ROCC_mu_pro->SetPoint(t, 1 - eff_b, eff_s);
        for(auto v : LRpi) if (v[0]<thr) count_b2++;
        eff_tot = double(count_b + count_b2)/(LRpro.size()+LRpi.size());
        pts_mu_other.push_back({1 - eff_tot,eff_s});

        count_s=0;
        count_b=0; 
        eff_b=0; 
        eff_s=0;

        for(const auto &v : LRmu) if (v[1]<thr) count_s++;
        eff_s = double(count_s)/LRmu.size();
        for(auto v : LRpi) if (v[1]<thr) count_b++;
        eff_b = double(count_b)/LRpi.size();
        pts_mu_pi.push_back({1 - eff_b,eff_s});
        //ROCC_mu_pi->SetPoint(t, 1 - eff_b, eff_s);
        count_s=0;
        count_b=0; 
        eff_b=0; 
        eff_s=0;

        for(const auto &v : LRpro) if (v[2]<thr) count_s++;
        eff_s = double(count_s)/LRpro.size();
        for(auto v : LRpi) if (v[2]<thr) count_b++;
        eff_b = double(count_b)/LRpi.size();
        pts_pro_pi.push_back({1 - eff_b,eff_s});
        //ROCC_pro_pi->SetPoint(t, 1 - eff_b, eff_s);
    }
    std::array<TGraph*,4> rocks = {
        ROCbuilder(pts_mu_pro, Form("rocCurve_%s_muon_proton", plane.c_str())),
        ROCbuilder(pts_mu_pi, Form("rocCurve_%s_muon_pion", plane.c_str())),
        ROCbuilder(pts_pro_pi, Form("rocCurve_%s_proton_pion", plane.c_str())),
        ROCbuilder(pts_mu_other, Form("rocCurve_%s_muon_other", plane.c_str()))
    };

    return rocks;
}

void hprobDensities()
{
    //gROOT->ProcessLine(".L libdaughtersInfo.so");
    std::array<std::string,3> particles={"muon", "proton", "pion"};
    TFile *f = TFile::Open("newHISTOdistro2d_general.root", "UPDATE");
    for(const auto& par : particles)
    {
        EventsData dat = load_data("mc2d_general", par);
        /*
        std::string dE_branch_name = (par=="muon") ? "dE_mu" : "dE_pro";
        std::string rr_branch_name = (par=="muon") ? "rr_mu" : "rr_pro";
        std::string dE_ind1_branch_name = (par=="muon") ? "dE_mu_ind1" : "dE_pro_ind1";
        std::string dE_ind2_branch_name = (par=="muon") ? "dE_mu_ind2" : "dE_pro_ind2";
        std::string rr_ind1_branch_name = (par=="muon") ? "rr_mu_ind1" : "rr_pro_ind1";
        std::string rr_ind2_branch_name = (par=="muon") ? "rr_mu_ind2" : "rr_pro_ind2";
        dat.tree->SetBranchStatus("*", 0);
        dat.tree->SetBranchStatus(dE_branch_name.c_str(), 1);
        dat.tree->SetBranchStatus(rr_branch_name.c_str(), 1);
        dat.tree->SetBranchStatus(dE_ind1_branch_name.c_str(), 1);
        dat.tree->SetBranchStatus(dE_ind2_branch_name.c_str(), 1);
        dat.tree->SetBranchStatus(rr_ind1_branch_name.c_str(), 1);
        dat.tree->SetBranchStatus(rr_ind2_branch_name.c_str(), 1);
        */

        cout << dat.tree->GetEntries() << " " << par << "tracks" << endl;
        TDirectory *pardir = (TDirectory*)f->Get(Form("%s",par.c_str()));
        
        for(double rr=1.5; rr<=25; rr+=0.5)
        {
            cout << "processing rr " << rr << endl;
            TDirectory *rrdir = (TDirectory*)pardir->Get(Form("%.1f", rr));
            rrdir->cd();

            int Nbins_coll = 300;
            int Nbins_coll_ind1 = 300;
            int Nbins_coll_ind1_ind2 = 300;

            std::vector<double> binlowedges_coll;
            std::vector<double> binlowedges_coll_ind1;
            std::vector<double> binlowedges_coll_ind1_ind2;
            for(int i=1; i<=Nbins_coll+1; i++)
            {
                binlowedges_coll.push_back((i-1)*30./Nbins_coll);
                binlowedges_coll_ind1.push_back((i-1)*30./Nbins_coll_ind1);
                binlowedges_coll_ind1_ind2.push_back((i-1)*30./Nbins_coll_ind1_ind2);
            }
            bool keep_coll = true;
            bool keep_coll_ind1 = true;
            bool keep_coll_ind1_ind2 = true;

            TH1D *dEdx_coll;
            TH1D *dEdx_coll_ind1;
            TH1D *dEdx_coll_ind1_ind2;

            //int iterations=1;
            //collection
            while(keep_coll)
            {
                //cout << "iterazione numero " << iterations << endl;
                //for(int j=0; j<binlowedges_coll.size(); j++){cout << j << " " << binlowedges_coll[j] << " ";}

                keep_coll = false;
                dEdx_coll = new TH1D(Form("dEdx_coll_rr_%.1f", rr),"", Nbins_coll, binlowedges_coll.data());

                for(int track=0; track<dat.tree->GetEntries(); track++)
                {
                    dat.tree->GetEntry(track);

                    if(par == "pion" && dat.track.end_process!=3)continue;
                    //cout << dat.track.end_process << endl;
                    for(int hit=0; hit<int(dat.track.dE->size()); hit++)
                    {
                        if(dat.track.rr->at(hit)<rr && dat.track.rr->at(hit)>=rr-0.5)
                        {
                            dEdx_coll->Fill(dat.track.dE->at(hit));
                        }   
                    }   
                }   


                for(int bin=Nbins_coll-1; bin>0; bin--)
                {
                    //cout << "index= " << bin << " " << "index+1 bin entries= " << dEdx_coll->GetBinContent(bin+1) << endl;
                    if(dEdx_coll->GetBinContent(bin+1) < 1) // bin+1 perché TH1D è 1-indexed
                    {
                        keep_coll = true;
                        binlowedges_coll.erase(binlowedges_coll.begin()+bin);
                        //for(int j=0; j<binlowedges_coll.size(); j++){cout << j << " " << binlowedges_coll[j] << " ";}
                        Nbins_coll--;
                    }
                }
                if(dEdx_coll->GetBinContent(1)<1){binlowedges_coll.erase(binlowedges_coll.begin()+1);}
                if(keep_coll){delete dEdx_coll;}
                //if(keep_coll)iterations+=1;
            }

            
            //collection + induction 1
            while(keep_coll_ind1)
            {
                keep_coll_ind1 = false;
                dEdx_coll_ind1 = new TH1D(Form("dEdx_coll_ind1_rr_%.1f", rr),"", Nbins_coll_ind1, binlowedges_coll_ind1.data());

                for(int track=0; track<dat.tree->GetEntries(); track++)
                {
                    dat.tree->GetEntry(track);

                    if(par == "pion" && dat.track.end_process!=3)continue;
                    //collection plane
                    for(int hit=0; hit<int(dat.track.dE->size()); hit++)
                    {
                        if(dat.track.rr->at(hit)<rr && dat.track.rr->at(hit)>=rr-0.5)
                        {
                            dEdx_coll_ind1->Fill(dat.track.dE->at(hit));
                        }   
                    } 
                    
                    //induction 1 plane
                    for(int hit=0; hit<int(dat.track.dEdx_ind1->size()); hit++)
                    {
                        if(dat.track.rr_ind1->at(hit)<rr && dat.track.rr_ind1->at(hit)>=rr-0.5)
                        {
                            dEdx_coll_ind1->Fill(dat.track.dEdx_ind1->at(hit));
                        }   
                    } 
                }   

                for(int bin=Nbins_coll_ind1-1; bin>=0; bin--)
                {
                    if(dEdx_coll_ind1->GetBinContent(bin+1) < 1) // bin+1 perché TH1D è 1-indexed
                    {
                        keep_coll_ind1 = true;
                        binlowedges_coll_ind1.erase(binlowedges_coll_ind1.begin()+bin);
                        Nbins_coll_ind1--;
                    }
                }
                if(dEdx_coll_ind1->GetBinContent(1)<1){binlowedges_coll_ind1.erase(binlowedges_coll_ind1.begin()+1);}
                if(keep_coll_ind1){delete dEdx_coll_ind1;}
            }

            //collection + induction 1 + induction 2
            while(keep_coll_ind1_ind2)
            {
                keep_coll_ind1_ind2 = false;
                dEdx_coll_ind1_ind2 = new TH1D(Form("dEdx_coll_ind1_ind2_rr_%.1f", rr),"", Nbins_coll_ind1_ind2, binlowedges_coll_ind1_ind2.data());

                for(int track=0; track<dat.tree->GetEntries(); track++)
                {
                    dat.tree->GetEntry(track);

                    if(par == "pion" && dat.track.end_process!=3)continue;
                    //collection plane
                    for(int hit=0; hit<int(dat.track.dE->size()); hit++)
                    {
                        if(dat.track.rr->at(hit)<rr && dat.track.rr->at(hit)>=rr-0.5)
                        {
                            dEdx_coll_ind1_ind2->Fill(dat.track.dE->at(hit));
                        }   
                    }
                    
                    //induction 1 plane
                    for(int hit=0; hit<int(dat.track.dEdx_ind1->size()); hit++)
                    {
                        if(dat.track.rr_ind1->at(hit)<rr && dat.track.rr_ind1->at(hit)>=rr-0.5)
                        {
                            dEdx_coll_ind1_ind2->Fill(dat.track.dEdx_ind1->at(hit));
                        }   
                    }

                    //induction 2 plane
                    for(int hit=0; hit<int(dat.track.dEdx_ind2->size()); hit++)
                    {
                        if(dat.track.rr_ind2->at(hit)<rr && dat.track.rr_ind2->at(hit)>=rr-0.5)
                        {
                            dEdx_coll_ind1_ind2->Fill(dat.track.dEdx_ind2->at(hit));
                        }   
                    }
                }   

                for(int bin=Nbins_coll_ind1_ind2-1; bin>=0; bin--)
                {
                    if(dEdx_coll_ind1_ind2->GetBinContent(bin+1) < 1) // bin+1 perché TH1D è 1-indexed
                    {
                        keep_coll_ind1_ind2 = true;
                        binlowedges_coll_ind1_ind2.erase(binlowedges_coll_ind1_ind2.begin()+bin);
                        Nbins_coll_ind1_ind2--;
                    }
                }
                if(dEdx_coll_ind1_ind2->GetBinContent(1)<1){binlowedges_coll_ind1_ind2.erase(binlowedges_coll_ind1_ind2.begin()+1);}
                if(keep_coll_ind1_ind2){delete dEdx_coll_ind1_ind2;}
            }
            

            dEdx_coll->Scale(1./dEdx_coll->Integral());
            dEdx_coll->Write(0,TObject::kOverwrite);

            dEdx_coll_ind1->Scale(1./dEdx_coll_ind1->Integral());
            dEdx_coll_ind1->Write(0,TObject::kOverwrite);

            dEdx_coll_ind1_ind2->Scale(1./dEdx_coll_ind1_ind2->Integral());
            dEdx_coll_ind1_ind2->Write(0,TObject::kOverwrite);  

        }
    }
}


void dEdxdistro()
{
    //gROOT->ProcessLine(".L libdaughtersInfo.so");
    std::array<std::string,3> particles={"muon", "proton", "pion"};
    //TFile *f = TFile::Open("newHISTOdistro2d_general.root", "UPDATE");
    TFile *f = new TFile("dEdxdistro.root","RECREATE");
    for(const auto& par : particles)
    {
        EventsData dat = load_data("mc2d_general", par);

        cout << dat.tree->GetEntries() << " " << par << "tracks" << endl;
        TDirectory *pardir = (TDirectory*)f->mkdir(Form("%s",par.c_str()));
        
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

                if(par=="pion" && dat.track.end_process!=3)continue;

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




void likelihood()
{    
    TFile * f = TFile::Open("newHISTOdistro2d_general.root", "READ");

    std::vector<TH1D*> probDensitiesMU_coll;
    std::vector<TH1D*> probDensitiesPRO_coll;
    std::vector<TH1D*> probDensitiesPI_coll;
    std::vector<TH1D*> probDensitiesMU_coll_ind1;
    std::vector<TH1D*> probDensitiesPRO_coll_ind1;
    std::vector<TH1D*> probDensitiesPI_coll_ind1;
    std::vector<TH1D*> probDensitiesMU_coll_ind1_ind2;
    std::vector<TH1D*> probDensitiesPRO_coll_ind1_ind2;
    std::vector<TH1D*> probDensitiesPI_coll_ind1_ind2;
    std::array<std::string,3> particle = {"muon", "proton", "pion"};
    for(auto &par : particle)
    {
        for(double i=1.5; i<=25.0; i+=0.5)
        {
            TDirectory *dpar = (TDirectory*)f->Get(Form("%s", par.c_str()));
            TDirectory *d = (TDirectory*)dpar->Get(Form("%.1f",i));

            TH1D *dEdx_coll = (TH1D*)d->Get(Form("dEdx_coll_rr_%.1f", i));
            if(par=="muon")probDensitiesMU_coll.push_back(dEdx_coll);
            if(par=="proton") probDensitiesPRO_coll.push_back(dEdx_coll);
            if(par=="pion") probDensitiesPI_coll.push_back(dEdx_coll);

            TH1D *dEdx_coll_ind1 = (TH1D*)d->Get(Form("dEdx_coll_ind1_rr_%.1f", i));
            if(par=="muon")probDensitiesMU_coll_ind1.push_back(dEdx_coll_ind1);
            if(par=="proton") probDensitiesPRO_coll_ind1.push_back(dEdx_coll_ind1);
            if(par=="pion") probDensitiesPI_coll_ind1.push_back(dEdx_coll_ind1);

            TH1D *dEdx_coll_ind1_ind2 = (TH1D*)d->Get(Form("dEdx_coll_ind1_ind2_rr_%.1f", i));
            if(par=="muon")probDensitiesMU_coll_ind1_ind2.push_back(dEdx_coll_ind1_ind2);
            if(par=="proton") probDensitiesPRO_coll_ind1_ind2.push_back(dEdx_coll_ind1_ind2);
            if(par=="pion") probDensitiesPI_coll_ind1_ind2.push_back(dEdx_coll_ind1_ind2);
        }
    }

    std::array<std::vector<TH1D*>,3> probDensities_coll={probDensitiesMU_coll,probDensitiesPRO_coll,probDensitiesPI_coll};
    std::array<std::vector<TH1D*>,3> probDensities_coll_ind1={probDensitiesMU_coll_ind1,probDensitiesPRO_coll_ind1,probDensitiesPI_coll_ind1};
    std::array<std::vector<TH1D*>,3> probDensities_coll_ind1_ind2={probDensitiesMU_coll_ind1_ind2,probDensitiesPRO_coll_ind1_ind2,probDensitiesPI_coll_ind1_ind2};

    TFile *outfile = new TFile("newnew_LR.root","RECREATE");
    TDirectory *d = (TDirectory*)outfile->mkdir("likelihood");

    std::vector<std::vector<double>> v_likelihood_ratio_coll_mu;  //each entry contains 3 likelihood ratios for muon tracks: mu_hyp/pro_hyp, mu_hyp/pi_hyp, pro_hyp/pi_hyp
    std::vector<std::vector<double>> v_likelihood_ratio_coll_ind1_mu;
    std::vector<std::vector<double>> v_likelihood_ratio_coll_ind1_ind2_mu;
    std::vector<std::vector<double>> v_likelihood_ratio_coll_pro;
    std::vector<std::vector<double>> v_likelihood_ratio_coll_ind1_pro;
    std::vector<std::vector<double>> v_likelihood_ratio_coll_ind1_ind2_pro;
    std::vector<std::vector<double>> v_likelihood_ratio_coll_pi;
    std::vector<std::vector<double>> v_likelihood_ratio_coll_ind1_pi;
    std::vector<std::vector<double>> v_likelihood_ratio_coll_ind1_ind2_pi;

    for(const auto &par : particle)
    {
        std::vector<TH1D*> likelihood_ratio_coll;
        std::vector<TH1D*> likelihood_ratio_coll_ind1;
        std::vector<TH1D*> likelihood_ratio_coll_ind1_ind2;

        std::array<std::string,3> hyps={"MUhyp_meno_PROhyp", "MUhyp_meno_PIhyp", "PROhyp_meno_PIhyp"};

        for(auto hyp : hyps)
        {
            TH1D *h1 = new TH1D(Form("%s_likelihood_ratio_coll_%s",par.c_str(),hyp.c_str()),"",100,-1,1);
            likelihood_ratio_coll.push_back(h1);
            TH1D *h2 = new TH1D(Form("%s_likelihood_ratio_coll_ind1_%s",par.c_str(),hyp.c_str()),"",100,-1,1);
            likelihood_ratio_coll_ind1.push_back(h2);
            TH1D *h3 = new TH1D(Form("%s_likelihood_ratio_coll_ind1_ind2_%s",par.c_str(),hyp.c_str()),"",100,-1,1);
            likelihood_ratio_coll_ind1_ind2.push_back(h3);
        }
        
        EventsData dat = load_data("mc2d_general", par);
        cout << dat.tree->GetEntries() << " " <<  par << " tracks" << endl;

        for(int track=0; track<dat.tree->GetEntries(); track++)
        {

          dat.tree->GetEntry(track);

          if(par=="pion" && dat.track.end_process!=3)continue;
            
            std::array<double,3> dEdx_lkl_coll; //3 particle hypothesis
            std::array<double,3> dEdx_lkl_coll_ind1;
            std::array<double,3> dEdx_lkl_coll_ind1_ind2;

            bool hasValidHits_coll=false;
            bool hasValidHits_ind1=false;
            bool hasValidHits_ind2=false;

            //ciclo sulla particle hypothesis
            for(int j=0; j<3; j++)
            {
                //collection plane
                auto density_coll = probDensities_coll[j];
                double log_lkh=0;
        
                for(int hit=0; hit<dat.track.dE->size(); hit++)
                {
                    if(dat.track.rr->at(hit)<25 && dat.track.rr->at(hit)>1. && dat.track.dE->at(hit)<30. /*&& dat.track.dE->at(hit)>1.*/)
                    {
                        hasValidHits_coll=true;
                        int idx=find_idx(dat.track.rr->at(hit));
                        int bin = density_coll[idx]->FindBin(dat.track.dE->at(hit));
                        log_lkh = log_lkh + (-1*std::log(density_coll[idx]->GetBinContent(bin)/density_coll[idx]->GetBinWidth(bin)));
                    }
                }
                if(hasValidHits_coll) dEdx_lkl_coll[j]=log_lkh;

                //collection plane + induction 1 plane
                auto density_coll_ind1 = probDensities_coll_ind1[j];
        
                for(int hit=0; hit<dat.track.dEdx_ind1->size(); hit++)
                {
                    if(dat.track.rr_ind1->at(hit)<25 && dat.track.rr_ind1->at(hit)>1. && dat.track.dEdx_ind1->at(hit)<30. /*&& dat.track.dEdx_ind1->at(hit)>1.*/)
                    {
                        hasValidHits_ind1=true;
                        int idx=find_idx(dat.track.rr_ind1->at(hit));
                        int bin = density_coll_ind1[idx]->FindBin(dat.track.dEdx_ind1->at(hit));
                        log_lkh = log_lkh + (-1*std::log(density_coll_ind1[idx]->GetBinContent(bin)/density_coll_ind1[idx]->GetBinWidth(bin)));
 
                    }
                }
                if(hasValidHits_coll || hasValidHits_ind1) dEdx_lkl_coll_ind1[j]=log_lkh;

                //collection plane + induction 1 plane + induction 2 plane
                auto density_coll_ind1_ind2 = probDensities_coll_ind1_ind2[j];
        
                for(int hit=0; hit<dat.track.dEdx_ind2->size(); hit++)
                {
                    if(dat.track.rr_ind2->at(hit)<25 && dat.track.rr_ind2->at(hit)>1. && dat.track.dEdx_ind2->at(hit)<30. /*&& dat.track.dEdx_ind2->at(hit)>1.*/)
                    {
                        hasValidHits_ind2=true;
                        int idx=find_idx(dat.track.rr_ind2->at(hit));
                        int bin = density_coll_ind1_ind2[idx]->FindBin(dat.track.dEdx_ind2->at(hit));
                        log_lkh = log_lkh + (-1*std::log(density_coll_ind1_ind2[idx]->GetBinContent(bin)/density_coll_ind1_ind2[idx]->GetBinWidth(bin)));
 
                    }
                }
                if(hasValidHits_coll || hasValidHits_ind1 || hasValidHits_ind2) dEdx_lkl_coll_ind1_ind2[j]=log_lkh;
               
            }//cycle over particle hypothesis


            //GETTING LIKELIHOOD RATIOS FOR EACH PARTICLE
            std::vector<double> lr_coll_temp;
            std::vector<double> lr_coll_ind1_temp;
            std::vector<double> lr_coll_ind1_ind2_temp;

            int index_hist=0;
            for(int k=0; k<3; k++)
            {
                for(int t=k+1; t<3; t++)
                {
                    //filling the istograms
                    if(hasValidHits_coll) likelihood_ratio_coll[index_hist]->Fill(std::atan((dEdx_lkl_coll[k]-dEdx_lkl_coll[t])/3.)*180/M_PI/90);
                    if(hasValidHits_coll || hasValidHits_ind1) likelihood_ratio_coll_ind1[index_hist]->Fill(std::atan((dEdx_lkl_coll_ind1[k]-dEdx_lkl_coll_ind1[t])/3.)*180/M_PI/90);
                    if(hasValidHits_coll || hasValidHits_ind1 || hasValidHits_ind2) likelihood_ratio_coll_ind1_ind2[index_hist]->Fill(std::atan((dEdx_lkl_coll_ind1_ind2[k]-dEdx_lkl_coll_ind1_ind2[t])/3.)*180/M_PI/90);
                    index_hist++;

                    //the likelihood ratios are: muon - proton, muon - pion, proton - pion
                    if(hasValidHits_coll) lr_coll_temp.push_back(std::atan((dEdx_lkl_coll[k]-dEdx_lkl_coll[t])/3.)*180/M_PI/90);
                    if(hasValidHits_coll || hasValidHits_ind1) lr_coll_ind1_temp.push_back(std::atan((dEdx_lkl_coll_ind1[k]-dEdx_lkl_coll_ind1[t])/3.)*180/M_PI/90);
                    if(hasValidHits_coll || hasValidHits_ind1 || hasValidHits_ind2) lr_coll_ind1_ind2_temp.push_back(std::atan((dEdx_lkl_coll_ind1_ind2[k]-dEdx_lkl_coll_ind1_ind2[t])/3.)*180/M_PI/90);
                }
            }    

            if(par=="muon")
            {   if(hasValidHits_coll)v_likelihood_ratio_coll_mu.push_back(lr_coll_temp); 
                if(hasValidHits_coll || hasValidHits_ind1)v_likelihood_ratio_coll_ind1_mu.push_back(lr_coll_ind1_temp); 
                if(hasValidHits_coll || hasValidHits_ind1 || hasValidHits_ind2)v_likelihood_ratio_coll_ind1_ind2_mu.push_back(lr_coll_ind1_ind2_temp);
            }
            if(par=="proton")
            {
                if(hasValidHits_coll)v_likelihood_ratio_coll_pro.push_back(lr_coll_temp); 
                if(hasValidHits_coll || hasValidHits_ind1)v_likelihood_ratio_coll_ind1_pro.push_back(lr_coll_ind1_temp); 
                if(hasValidHits_coll || hasValidHits_ind1 || hasValidHits_ind2)v_likelihood_ratio_coll_ind1_ind2_pro.push_back(lr_coll_ind1_ind2_temp);
            }
            if(par=="pion")
            {   
                if(hasValidHits_coll)v_likelihood_ratio_coll_pi.push_back(lr_coll_temp); 
                if(hasValidHits_coll || hasValidHits_ind1)v_likelihood_ratio_coll_ind1_pi.push_back(lr_coll_ind1_temp); 
                if(hasValidHits_coll || hasValidHits_ind1 || hasValidHits_ind2)v_likelihood_ratio_coll_ind1_ind2_pi.push_back(lr_coll_ind1_ind2_temp);
            }
        }

        
        outfile->cd();
        d->cd();
        for(int histo=0; histo<3; histo++)
        {
            likelihood_ratio_coll[histo]->Scale(1./likelihood_ratio_coll[histo]->Integral());
            likelihood_ratio_coll_ind1[histo]->Scale(1./likelihood_ratio_coll_ind1[histo]->Integral());
            likelihood_ratio_coll_ind1_ind2[histo]->Scale(1./likelihood_ratio_coll_ind1_ind2[histo]->Integral());
            likelihood_ratio_coll[histo]->Write(0,TObject::kOverwrite);
            likelihood_ratio_coll_ind1[histo]->Write(0,TObject::kOverwrite);
            likelihood_ratio_coll_ind1_ind2[histo]->Write(0,TObject::kOverwrite);
        }

    }//cycle on particle type

    //making roc curves
    
    std::array<TGraph*,4> roc_coll = getROCcurvesLR(20000, v_likelihood_ratio_coll_mu, v_likelihood_ratio_coll_pro, v_likelihood_ratio_coll_pi, "coll");
    std::array<TGraph*,4> roc_coll_ind1 = getROCcurvesLR(20000,v_likelihood_ratio_coll_ind1_mu, v_likelihood_ratio_coll_ind1_pro,v_likelihood_ratio_coll_ind1_pi, "coll_ind1");
    std::array<TGraph*,4> roc_coll_ind1_ind2 = getROCcurvesLR(20000,v_likelihood_ratio_coll_ind1_ind2_mu,v_likelihood_ratio_coll_ind1_ind2_pro,v_likelihood_ratio_coll_ind1_ind2_pi,"coll_ind1_ind2");

    d->cd();
    std::array<std::array<TGraph*,4>,3> rocs ={roc_coll,roc_coll_ind1,roc_coll_ind1_ind2};
    for(int i=0; i<3; i++){ for(int j=0; j<4; j++) rocs[i][j]->Write(0,TObject::kOverwrite);}
    
}





void chi2()
{

    //TFile *f = TFile::Open("likelihood_ratios.root","UPDATE");
    TFile *f = new TFile("newnew_LR.root","RECREATE");
    //TDirectory *d = (TDirectory*)f->Get("chi2");

    std::array<std::string,3> pars = {"muon","proton","pion"};
    std::array<std::array<std::vector<double>,2>,3> chi2_coll;
    std::array<std::array<std::vector<double>,2>,3> chi2_coll_ind1;
    std::array<std::array<std::vector<double>,2>,3> chi2_coll_ind1_ind2;


    for(int par=0; par<3; par++)
    {
        TH1D *hchi2_coll_asmu = new TH1D(Form("%s_hchi2_coll_asmu",pars[par].c_str()), "", 200, 0, 200);
        TH1D *hchi2_coll_aspro = new TH1D(Form("%s_hchi2_coll_aspro",pars[par].c_str()), "", 500, 0, 500);
        TH1D *hchi2_coll_ind1_asmu = new TH1D(Form("%s_hchi2_coll_ind1_asmu",pars[par].c_str()), "", 200, 0, 200);
        TH1D *hchi2_coll_ind1_aspro = new TH1D(Form("%s_hchi2_coll_ind1_aspro",pars[par].c_str()), "", 500, 0, 500);
        TH1D *hchi2_coll_ind1_ind2_asmu = new TH1D(Form("%s_hchi2_coll_ind1_ind2_asmu",pars[par].c_str()), "", 200, 0, 200);
        TH1D *hchi2_coll_ind1_ind2_aspro = new TH1D(Form("%s_hchi2_coll_ind1_ind2_aspro",pars[par].c_str()), "", 500, 0, 500);
        
        //TH1D *h_diff_chi2_coll = new TH1D(Form("%s_h_diff_chi2_coll",pars[par].c_str()),"",400,-400,400);
        //TH1D *h_diff_chi2_coll_ind1 = new TH1D(Form("%s_h_diff_chi2_coll_ind1",pars[par].c_str()),"",400,-400,400);
        //TH1D *h_diff_chi2_coll_ind1_ind2 = new TH1D(Form("%s_h_diff_chi2_coll_ind1_ind2",pars[par].c_str()),"",400,-400,400);

        EventsData dat = load_data("mc2d_general",pars[par]);
        cout << dat.tree->GetEntries() << " " << pars[par] << " tracks" << endl;
      
        for(int track=0; track<dat.tree->GetEntries(); track++)
        {
            dat.tree->GetEntry(track);

            if(pars[par]=="pion" && dat.track.end_process!=3)continue;

            std::vector<double> dEdx;
            std::vector<double> rr;

            //collection
            for(int hit=0; hit<dat.track.dE->size(); hit++)
            {
                dEdx.push_back(dat.track.dE->at(hit));
                rr.push_back(dat.track.rr->at(hit));
            }   
            if(has_valid_hits(rr,dEdx))
            {
                for(int i=0; i<2; i++) chi2_coll[par][i].push_back(chi2_ALG(dEdx,rr,0,25)[i]);
                //cout << pars[par] << " " << chi2_ALG(dEdx,rr,0,25)[0] << " " << chi2_ALG(dEdx,rr,0,25)[1] << endl;
                hchi2_coll_asmu->Fill(chi2_ALG(dEdx,rr,0,25)[0]);
                hchi2_coll_aspro->Fill(chi2_ALG(dEdx,rr,0,25)[1]);
                //h_diff_chi2_coll->Fill(chi2_ALG(dEdx,rr,0,25)[0]-chi2_ALG(dEdx,rr,0,25)[1]);
            }

            //collection + induction 1
            for(int hit=0; hit<dat.track.dEdx_ind1->size(); hit++)
            {
                dEdx.push_back(dat.track.dEdx_ind1->at(hit));
                rr.push_back(dat.track.rr_ind1->at(hit));
            }
            if(has_valid_hits(rr,dEdx))
            {
                for(int i=0; i<2; i++) chi2_coll_ind1[par][i].push_back(chi2_ALG(dEdx,rr,0,25)[i]);
                hchi2_coll_ind1_asmu->Fill(chi2_ALG(dEdx,rr,0,25)[0]);
                hchi2_coll_ind1_aspro->Fill(chi2_ALG(dEdx,rr,0,25)[1]);
                //h_diff_chi2_coll_ind1->Fill(chi2_ALG(dEdx,rr,0,25)[0]-chi2_ALG(dEdx,rr,0,25)[1]);
            }

            //collection + induction 1 + induction 2
            for(int hit=0; hit<dat.track.dEdx_ind2->size(); hit++)
            {
                dEdx.push_back(dat.track.dEdx_ind2->at(hit));
                rr.push_back(dat.track.rr_ind2->at(hit));
            }
            if(has_valid_hits(rr,dEdx))
            {
                for(int i=0; i<2; i++) chi2_coll_ind1_ind2[par][i].push_back(chi2_ALG(dEdx,rr,0,25)[i]);
                hchi2_coll_ind1_ind2_asmu->Fill(chi2_ALG(dEdx,rr,0,25)[0]);
                hchi2_coll_ind1_ind2_aspro->Fill(chi2_ALG(dEdx,rr,0,25)[1]);
                //h_diff_chi2_coll_ind1_ind2->Fill(chi2_ALG(dEdx,rr,0,25)[0]-chi2_ALG(dEdx,rr,0,25)[1]);
            }
        }

        f->cd();
        hchi2_coll_asmu->Scale(1./hchi2_coll_asmu->Integral());
        hchi2_coll_aspro->Scale(1./hchi2_coll_aspro->Integral());
        hchi2_coll_ind1_asmu->Scale(1./hchi2_coll_ind1_asmu->Integral());
        hchi2_coll_ind1_aspro->Scale(1./hchi2_coll_ind1_aspro->Integral());
        hchi2_coll_ind1_ind2_asmu->Scale(1./hchi2_coll_ind1_ind2_asmu->Integral());
        hchi2_coll_ind1_ind2_aspro->Scale(1./hchi2_coll_ind1_ind2_aspro->Integral());
        hchi2_coll_asmu->Write(0,TObject::kOverwrite);
        hchi2_coll_aspro->Write(0,TObject::kOverwrite);
        hchi2_coll_ind1_asmu->Write(0,TObject::kOverwrite);
        hchi2_coll_ind1_aspro->Write(0,TObject::kOverwrite);
        hchi2_coll_ind1_ind2_asmu->Write(0,TObject::kOverwrite);
        hchi2_coll_ind1_ind2_aspro->Write(0,TObject::kOverwrite);

        //h_diff_chi2_coll->Scale(1./h_diff_chi2_coll->Integral());
        //h_diff_chi2_coll_ind1->Scale(1./h_diff_chi2_coll_ind1->Integral());
        //h_diff_chi2_coll_ind1_ind2->Scale(1./h_diff_chi2_coll_ind1_ind2->Integral());
        //h_diff_chi2_coll->Write(0,TObject::kOverwrite);
        //h_diff_chi2_coll_ind1->Write(0,TObject::kOverwrite);
        //h_diff_chi2_coll_ind1_ind2->Write(0,TObject::kOverwrite);
    }
    
    //std::array<TGraph*,4> rocCurve_chi2_coll = GetROCcurveChi2new(f,"rocCurve_chi2_coll", chi2_coll[0], chi2_coll[1], chi2_coll[2]);
    //std::array<TGraph*,4> rocCurve_chi2_coll_ind1 = GetROCcurveChi2new(f,"rocCurve_chi2_coll_ind1", chi2_coll_ind1[0], chi2_coll_ind1[1], chi2_coll_ind1[2]);
    //std::array<TGraph*,4> rocCurve_chi2_coll_ind1_ind2 = GetROCcurveChi2new(f,"rocCurve_chi2_coll_ind1_ind2", chi2_coll_ind1_ind2[0], chi2_coll_ind1_ind2[1], chi2_coll_ind1_ind2[2]);

    std::array<TGraph*,4> rocCurve_chi2_coll = GetROCcurveChi2("rocCurve_chi2_coll", 200, 500, 0, 500,chi2_coll[0], chi2_coll[1], chi2_coll[2]);
    //std::array<TGraph*,4> rocCurve_chi2_coll_ind1 = GetROCcurveChi2("rocCurve_chi2_coll_ind1", 200, 500, 0, 500, chi2_coll_ind1[0], chi2_coll_ind1[1], chi2_coll_ind1[2]);
    //std::array<TGraph*,4> rocCurve_chi2_coll_ind1_ind2 = GetROCcurveChi2("rocCurve_chi2_coll_ind1_ind2", 200, 500, 0, 500, chi2_coll_ind1_ind2[0], chi2_coll_ind1_ind2[1], chi2_coll_ind1_ind2[2]);

    f->cd();
    for(int i=0; i<4; i++){rocCurve_chi2_coll[i]->Write(0,TObject::kOverwrite); /*rocCurve_chi2_coll_ind1[i]->Write(0,TObject::kOverwrite); rocCurve_chi2_coll_ind1_ind2[i]->Write(0,TObject::kOverwrite);*/}

}