#include "make_prob_densities.C"


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

/*
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

std::array<TGraph*,4> getROCcurvesLR(
    int N, 
    const std::vector<std::array<double,15>> &LRmu_class0, 
    const std::vector<std::array<double,15>> &LRmu_class1, 
    const std::vector<std::array<double,15>> &LRpro_class2,
    const std::vector<std::array<double,15>> &LRpro_class3,
    const std::vector<std::array<double,15>> &LRpi_class4,
    const std::vector<std::array<double,15>> &LRpi_class5,
)
{
    
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

        for(const auto &v : LRmu_class0) if(v[0]<thr) count_s++;
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
*/

int reco_classification(std::array<double,15> likelihood_ratios)
{
    int reco_class=-1;
    //muon rising
    if( likelihood_ratios[0]<=0 && 
        likelihood_ratios[1]<=0 &&
        //likelihood_ratios[2]<=0 &&
        likelihood_ratios[3]<=0 &&
        likelihood_ratios[4]<=0 &&
        likelihood_ratios[5]<=0 &&
        //likelihood_ratios[6]<=0 &&
        likelihood_ratios[7]>0 &&
        likelihood_ratios[8]<=0 &&
        likelihood_ratios[9]>0 &&
        likelihood_ratios[10]>0 &&
        //likelihood_ratios[11]>0 &&
        //likelihood_ratios[12]>0 &&
        //likelihood_ratios[13]>0 &&
        likelihood_ratios[14]<=0 
    )reco_class = 0;

    //muon mip
    if( likelihood_ratios[0]>0 && 
        likelihood_ratios[1]<=0 &&
        //likelihood_ratios[2]<=0 &&
        likelihood_ratios[3]<=0 &&
        likelihood_ratios[4]>0 &&
        likelihood_ratios[5]<=0 &&
        //likelihood_ratios[6]<=0 &&
        likelihood_ratios[7]<=0 &&
        likelihood_ratios[8]<=0 &&
        likelihood_ratios[9]>0 &&
        likelihood_ratios[10]>0 &&
        //likelihood_ratios[11]>0 &&
        //likelihood_ratios[12]>0 &&
        //likelihood_ratios[13]>0 &&
        likelihood_ratios[14]>0 
    )reco_class = 1;

    //proton rising
    if( likelihood_ratios[0]>0 && 
        likelihood_ratios[1]>0 &&
        //likelihood_ratios[2]>0 &&
        likelihood_ratios[3]>0 &&
        likelihood_ratios[4]>0 &&
        likelihood_ratios[5]>0 &&
        //likelihood_ratios[6]>0 &&
        likelihood_ratios[7]<=0 &&
        likelihood_ratios[8]>0 &&
        likelihood_ratios[9]<=0 &&
        likelihood_ratios[10]<=0 &&
        //likelihood_ratios[11]<=0 &&
        //likelihood_ratios[12]<=0 &&
        //likelihood_ratios[13]<=0 &&
        likelihood_ratios[14]>0 
    )reco_class = 2;

    //proton interacting
    if( likelihood_ratios[0]>0 && 
        likelihood_ratios[1]>0 &&
        //likelihood_ratios[2]>0 &&
        likelihood_ratios[3]>0 &&
        likelihood_ratios[4]>0 &&
        likelihood_ratios[5]>0 &&
        //likelihood_ratios[6]>0 &&
        likelihood_ratios[7]<=0 &&
        likelihood_ratios[8]>0 &&
        likelihood_ratios[9]>0 &&
        likelihood_ratios[10]<=0 &&
        //likelihood_ratios[11]<=0 &&
        //likelihood_ratios[12]<=0 &&
        //likelihood_ratios[13]<=0 &&
        likelihood_ratios[14]>0 
    )reco_class = 3;

    //pion rising
    if( likelihood_ratios[0]<=0 && 
        likelihood_ratios[1]<=0 &&
        //likelihood_ratios[2]<=0 &&
        likelihood_ratios[3]>0 &&
        likelihood_ratios[4]>0 &&
        likelihood_ratios[5]<=0 &&
        //likelihood_ratios[6]<=0 &&
        likelihood_ratios[7]>0 &&
        likelihood_ratios[8]>0 &&
        likelihood_ratios[9]>0 &&
        likelihood_ratios[10]>0 &&
        //likelihood_ratios[11]>0 &&
        //likelihood_ratios[12]>0 &&
        //likelihood_ratios[13]>0 &&
        likelihood_ratios[14]<=0 
    )reco_class = 4;

    //pion interacting
    if( likelihood_ratios[0]>0 && 
        likelihood_ratios[1]<=0 &&
        //likelihood_ratios[2]<=0 &&
        likelihood_ratios[3]<=0 &&
        likelihood_ratios[4]>0 &&
        likelihood_ratios[5]<=0 &&
        //likelihood_ratios[6]<=0 &&
        likelihood_ratios[7]<=0 &&
        likelihood_ratios[8]>0 &&
        likelihood_ratios[9]>0 &&
        likelihood_ratios[10]>0 &&
        //likelihood_ratios[11]>0 &&
        //likelihood_ratios[12]>0 &&
        //likelihood_ratios[13]>0 &&
        likelihood_ratios[14]>0 
    )reco_class = 5;

    return reco_class;
}

void likelihood()
{    
    TH2D * confusion_matrix = new TH2D("confusion_mactrix","",6,-0.5,5.5,6,-0.5,5.5);
    std::ofstream dump_lkl_ratios("dump_lkl_ratios_withdepE.txt");
    dump_lkl_ratios << "class" << " ";
    for(int i=0; i<15; i++){dump_lkl_ratios << "lr" << i << " ";}
    dump_lkl_ratios << "depE" << endl;
    dump_lkl_ratios << endl;
    TFile * f = TFile::Open("HISTO_prob_densities_6_classes.root", "READ");

    std::vector<TH1D*> probDensitiesMU_class0;
    std::vector<TH1D*> probDensitiesMU_class1;
    std::vector<TH1D*> probDensitiesPRO_class2;
    std::vector<TH1D*> probDensitiesPRO_class3;
    std::vector<TH1D*> probDensitiesPI_class4;
    std::vector<TH1D*> probDensitiesPI_class5;

    std::array<std::vector<TH1D*>,6> probDensities={
        probDensitiesMU_class0,
        probDensitiesMU_class1,
        probDensitiesPRO_class2,
        probDensitiesPRO_class3,
        probDensitiesPI_class4,
        probDensitiesPI_class5
    };

    ofstream check_prob_distro("check_prob_distro.txt");
    check_prob_distro << "class rr bin bin_content dedx_centroid" << endl;
    std::vector<std::pair<std::string,int>> classes = {{"muon_class0",0}, {"muon_class1",1}, {"proton_class2",2}, {"proton_class3",3}, {"pion_class4",4}, {"pion_class5",5}};

    for(const auto &clas : classes)
    {
        for(double i=1.5; i<=25.0; i+=0.5)
        {
            TDirectory *dclass = (TDirectory*)f->Get(clas.first.c_str());
            TDirectory *d = (TDirectory*)dclass->Get(Form("%.1f",i));

            TH1D *dEdx_coll = (TH1D*)d->Get(Form("dEdx_coll_rr_%.1f", i));

            for(int bin=1; bin<=dEdx_coll->GetNbinsX(); bin++)
            {
                check_prob_distro << clas.first << " " << i << " " << bin << " " << dEdx_coll->GetBinContent(bin) << " " << dEdx_coll->GetBinCenter(bin) << endl;
            }
            

            probDensities[clas.second].push_back(dEdx_coll);
        }
    }
    
    std::array<std::string,3> particles = {"muon", "proton", "pion"};

    TFile *outfile = new TFile("outfile_LR_6_classes.root","RECREATE");

    std::vector<std::array<double,15>> v_likelihood_ratios_mu_class0;
    std::vector<std::array<double,15>> v_likelihood_ratios_mu_class1;
    std::vector<std::array<double,15>> v_likelihood_ratios_pro_class2;
    std::vector<std::array<double,15>> v_likelihood_ratios_pro_class3;
    std::vector<std::array<double,15>> v_likelihood_ratios_pi_class4;
    std::vector<std::array<double,15>> v_likelihood_ratios_pi_class5;
 
    int class_index=-1;
    for(const auto &par : particles)
    {
        for(int subclass=0; subclass<2; subclass++)
        {
            int tracks_with_no_valid_hits = 0;
            class_index++;

            TDirectory *d = (TDirectory*)outfile->mkdir(Form("%s",classes[class_index].first.c_str()));
            std::vector<TH1D*> likelihood_ratio;

            std::array<std::string,15> hyps={
                "MU_class0_meno_MU_class1", 
                "MU_class0_meno_PRO_class2",
                "MU_class0_meno_PRO_class3",
                "MU_class0_meno_PI_class4",
                "MU_class0_meno_PI_class5",
                "MU_class1_meno_PRO_class2",
                "MU_class1_meno_PRO_class3",
                "MU_class1_meno_PI_class4",
                "MU_class1_meno_PI_class5",
                "PRO_class2_meno_PRO_class3",
                "PRO_class2_meno_PI_class4",
                "PRO_class2_meno_PI_class5",
                "PRO_class3_meno_PI_class4",
                "PRO_class3_meno_PI_class5",
                "PI_class4_meno_PI_class5"
            };
            //since arctan is an odd function, the negative of the likelihood ration is equal to the inverse likelihood ratio

            for(auto hyp : hyps)
            {
                TH1D *h1 = new TH1D(Form("%s_class%d_likelihood_ratio_%s", par.c_str(), classes[class_index].second, hyp.c_str()),"",100,-1,1);
                likelihood_ratio.push_back(h1);
            }
        
            EventsData dat = load_data("mc2d_general_contained", par);
            //cout << dat.tree->GetEntries() << " " <<  par << " tracks" << endl;

            for(int track=0; track<dat.tree->GetEntries(); track++)
            {
                std::array<double,15> likelihood_ratios;

                //choosing only track with odd index
                //if(track % 2 ==0)continue;

                //true selection
            
                dat.tree->GetEntry(track);

                int true_class = true_selection(dat,par);
                if(true_class!=classes[class_index].second)continue;
            
                std::array<double,6> lkl; //it contains the likelihoods in the 6 hypotheses for the current track

                bool hasValidHits=false;
                //ciclo sulla particle hypothesis
                for(int j=0; j<6; j++)
                {
                    auto density = probDensities[j];
                    double log_lkh=0;
        
                    for(int hit=0; hit<dat.track.dE->size(); hit++)
                    {
                        if(dat.track.rr->at(hit)<25 && dat.track.rr->at(hit)>1. && dat.track.dE->at(hit)<30.)
                        {
                            hasValidHits=true;
                            int idx=find_idx(dat.track.rr->at(hit));
                            int bin = density[idx]->FindBin(dat.track.dE->at(hit));
                            log_lkh = log_lkh + (-1*std::log(density[idx]->GetBinContent(bin)/density[idx]->GetBinWidth(bin)));
                            if(std::isinf(log_lkh)){ cout << dat.track.rr->at(hit) << " " << dat.track.dE->at(hit) << " | " << density[idx]->GetBinContent(bin) << " " << std::log(density[idx]->GetBinContent(bin)) << " " << density[idx]->GetBinWidth(bin) << " | " << par << " " << subclass << " " << idx << "(" << dat.track.rr->at(hit) << ")" << " " << bin << endl; break;}
                        }
                    }
                    lkl[j]=log_lkh;
                }//cycle over particle hypothesis
                if(!hasValidHits){tracks_with_no_valid_hits++; continue;}

                //GETTING LIKELIHOOD RATIOS FOR EACH TRACK
                int index_hist=-1;
                //dump_lkl_ratios << classes[class_index].second << " ";
                dump_lkl_ratios << classes[class_index].first << " ";
                for(int k=0; k<6; k++)
                {
                    for(int t=k+1; t<6; t++)
                    {
                        index_hist++;

                        //if(isnan(std::atan((lkl[k]-lkl[t])/3.)*180/M_PI/90))cout << "nan value ecountered " << lkl[k] << " " << lkl[t] << endl;
                        dump_lkl_ratios << std::atan((lkl[k]-lkl[t])/3.)*180/M_PI/90 << " ";
                        //all the likelihood ratios
                        likelihood_ratios[index_hist]=(std::atan((lkl[k]-lkl[t])/3.)*180/M_PI/90);
                        
                        //filling the istogram
                        likelihood_ratio[index_hist]->Fill(std::atan((lkl[k]-lkl[t])/3.)*180/M_PI/90);
                    }
                }  

                double dep_energy=0;
                for(int hit=0; hit<dat.track.dE->size(); hit++)
                {
                    if(dat.track.rr->at(hit)>5.)continue;
                    dep_energy = dep_energy + dat.track.dE->at(hit)*dat.track.pitch->at(hit);
                }
                dump_lkl_ratios << dep_energy << endl;

                //SELECTION
                int reco_class=reco_classification(likelihood_ratios);
                confusion_matrix->Fill(true_class,reco_class);

                //storing likelihood ratios
                if(par=="muon")
                {
                    if(subclass==0)v_likelihood_ratios_mu_class0.push_back(likelihood_ratios);
                    else v_likelihood_ratios_mu_class1.push_back(likelihood_ratios);
                }
                if(par=="proton")
                {
                    if(subclass==0)v_likelihood_ratios_pro_class2.push_back(likelihood_ratios);
                    else v_likelihood_ratios_pro_class3.push_back(likelihood_ratios);
                }
                if(par=="pion")
                {
                    if(subclass==0)v_likelihood_ratios_pi_class4.push_back(likelihood_ratios);
                    else v_likelihood_ratios_pi_class5.push_back(likelihood_ratios);
                }  
            }//cycle over all tracks

            outfile->cd();
            d->cd();
            for(int histo=0; histo<15; histo++)
            {
                likelihood_ratio[histo]->Scale(1./likelihood_ratio[histo]->Integral());
                likelihood_ratio[histo]->Write(0,TObject::kOverwrite);
            }
        
            cout << tracks_with_no_valid_hits << " " << par << " of subclass " << subclass << " with no valid hits in COLL plane" << endl;

        }//cycle over subclasses
 
    }//cycle on particle type

    //confusion matrix
    for(int bin=1; bin<=6; bin++)
    {
        double total=0;
        for(int binbin=1; binbin<=6; binbin++)
        {
            total = total + confusion_matrix->GetBinContent(bin, binbin);
        }
        for(int binbin=1; binbin<=6; binbin++)
        {
            confusion_matrix->SetBinContent(bin,binbin,confusion_matrix->GetBinContent(bin, binbin)/total);
        }
    }
    outfile->cd();
    confusion_matrix->Write(0,TObject::kOverwrite);

    //making roc curves
    //std::array<TGraph*,4> roc_coll = getROCcurvesLR(20000, v_likelihood_ratio_coll_mu, v_likelihood_ratio_coll_pro, v_likelihood_ratio_coll_pi, "coll");
    //std::array<TGraph*,4> roc_coll_ind1 = getROCcurvesLR(20000,v_likelihood_ratio_coll_ind1_mu, v_likelihood_ratio_coll_ind1_pro,v_likelihood_ratio_coll_ind1_pi, "coll_ind1");
    //std::array<TGraph*,4> roc_coll_ind1_ind2 = getROCcurvesLR(20000,v_likelihood_ratio_coll_ind1_ind2_mu,v_likelihood_ratio_coll_ind1_ind2_pro,v_likelihood_ratio_coll_ind1_ind2_pi,"coll_ind1_ind2");

    //d->cd();
    //std::array<std::array<TGraph*,4>,3> rocs ={roc_coll,roc_coll_ind1,roc_coll_ind1_ind2};
    //sfor(int i=0; i<3; i++){ for(int j=0; j<4; j++) rocs[i][j]->Write(0,TObject::kOverwrite);}
    
}

void likelihoodRM()
{    
    TH2D * confusion_matrix = new TH2D("confusion_mactrix","",6,-0.5,5.5,6,-0.5,5.5);
    std::ofstream dump_lkl_ratios("dump_lkl_ratios_RM_gbdt_withdepE.txt");
    //dump_lkl_ratios << "class" << " ";
    //for(int i=0; i<15; i++){dump_lkl_ratios << "lr" << i << " ";}
    //dump_lkl_ratios << endl;
    TFile * f = TFile::Open("HISTO_prob_densities_RM_6_classes.root", "READ");

    std::vector<TH1D*> probDensitiesMU_class0;
    std::vector<TH1D*> probDensitiesMU_class1;
    std::vector<TH1D*> probDensitiesPRO_class2;
    std::vector<TH1D*> probDensitiesPRO_class3;
    std::vector<TH1D*> probDensitiesPI_class4;
    std::vector<TH1D*> probDensitiesPI_class5;

    std::array<std::vector<TH1D*>,6> probDensities={
        probDensitiesMU_class0,
        probDensitiesMU_class1,
        probDensitiesPRO_class2,
        probDensitiesPRO_class3,
        probDensitiesPI_class4,
        probDensitiesPI_class5
    };

    //ofstream check_prob_distro("check_prob_distro.txt");
    //check_prob_distro << "class rr bin bin_content dedx_centroid" << endl;
    std::vector<std::pair<std::string,int>> classes = {{"muon_class0",0}, {"muon_class1",1}, {"proton_class2",2}, {"proton_class3",3}, {"pion_class4",4}, {"pion_class5",5}};

    for(const auto &clas : classes)
    {
        for(double i=1.5; i<=25.0; i+=0.5)
        {
            TDirectory *dclass = (TDirectory*)f->Get(clas.first.c_str());
            TDirectory *d = (TDirectory*)dclass->Get(Form("%.1f",i));

            TH1D *dEdx_coll = (TH1D*)d->Get(Form("dEdx_coll_rr_%.1f", i));

            //for(int bin=1; bin<=dEdx_coll->GetNbinsX(); bin++)
            //{
                //check_prob_distro << clas.first << " " << i << " " << bin << " " << dEdx_coll->GetBinContent(bin) << " " << dEdx_coll->GetBinCenter(bin) << endl;
            //}
            

            probDensities[clas.second].push_back(dEdx_coll);
        }
    }
    
    std::array<std::string,3> particles = {"muon", "proton", "pion"};

    TFile *outfile = new TFile("outfile_LR_RM_6_classes.root","RECREATE");

    std::vector<std::array<double,15>> v_likelihood_ratios_mu_class0;
    std::vector<std::array<double,15>> v_likelihood_ratios_mu_class1;
    std::vector<std::array<double,15>> v_likelihood_ratios_pro_class2;
    std::vector<std::array<double,15>> v_likelihood_ratios_pro_class3;
    std::vector<std::array<double,15>> v_likelihood_ratios_pi_class4;
    std::vector<std::array<double,15>> v_likelihood_ratios_pi_class5;
 
    int class_index=-1;
    for(const auto &par : particles)
    {
        for(int subclass=0; subclass<2; subclass++)
        {
            int tracks_with_no_valid_hits = 0;
            class_index++;

            TDirectory *d = (TDirectory*)outfile->mkdir(Form("%s",classes[class_index].first.c_str()));
            std::vector<TH1D*> likelihood_ratio;

            std::array<std::string,15> hyps={
                "MU_class0_meno_MU_class1", 
                "MU_class0_meno_PRO_class2",
                "MU_class0_meno_PRO_class3",
                "MU_class0_meno_PI_class4",
                "MU_class0_meno_PI_class5",
                "MU_class1_meno_PRO_class2",
                "MU_class1_meno_PRO_class3",
                "MU_class1_meno_PI_class4",
                "MU_class1_meno_PI_class5",
                "PRO_class2_meno_PRO_class3",
                "PRO_class2_meno_PI_class4",
                "PRO_class2_meno_PI_class5",
                "PRO_class3_meno_PI_class4",
                "PRO_class3_meno_PI_class5",
                "PI_class4_meno_PI_class5"
            };
            //since arctan is an odd function, the negative of the likelihood ration is equal to the inverse likelihood ratio

            for(auto hyp : hyps)
            {
                TH1D *h1 = new TH1D(Form("%s_class%d_likelihood_ratio_%s", par.c_str(), classes[class_index].second, hyp.c_str()),"",100,-1,1);
                likelihood_ratio.push_back(h1);
            }
        
            EventsData dat = load_data("mc2d_general_contained", par);
            //cout << dat.tree->GetEntries() << " " <<  par << " tracks" << endl;

            for(int track=0; track<dat.tree->GetEntries(); track++)
            {
                std::array<double,15> likelihood_ratios;

                //choosing only track with odd index
                //if(track % 2 ==0)continue;

                //true selection
            
                dat.tree->GetEntry(track);

                int true_class = true_selection_RM(dat,par);
                if(true_class!=classes[class_index].second)continue;
            
                std::array<double,6> lkl; //it contains the likelihoods in the 6 hypotheses for the current track

                bool hasValidHits=false;
                //ciclo sulla particle hypothesis
                for(int j=0; j<6; j++)
                {
                    auto density = probDensities[j];
                    double log_lkh=0;
        
                    for(int hit=0; hit<dat.track.rm_coll->size(); hit++)
                    {
                        if(dat.track.rr_rm_coll->at(hit)<25 && dat.track.rr_rm_coll->at(hit)>1. && dat.track.rm_coll->at(hit)<30.)
                        {
                            hasValidHits=true;
                            int idx=find_idx(dat.track.rr_rm_coll->at(hit));
                            int bin = density[idx]->FindBin(dat.track.rm_coll->at(hit));
                            log_lkh = log_lkh + (-1*std::log(density[idx]->GetBinContent(bin)/density[idx]->GetBinWidth(bin)));
                            if(std::isinf(log_lkh)){ cout << dat.track.rr_rm_coll->at(hit) << " " << dat.track.rm_coll->at(hit) << " | " << density[idx]->GetBinContent(bin) << " " << std::log(density[idx]->GetBinContent(bin)) << " " << density[idx]->GetBinWidth(bin) << " | " << par << " " << subclass << " " << idx << "(" << dat.track.rr_rm_coll->at(hit) << ")" << " " << bin << endl; break;}
                        }
                    }
                    lkl[j]=log_lkh;
                }//cycle over particle hypothesis
                if(!hasValidHits){tracks_with_no_valid_hits++; continue;}

                //GETTING LIKELIHOOD RATIOS FOR EACH TRACK
                int index_hist=-1;
                dump_lkl_ratios << classes[class_index].second << " ";
                //dump_lkl_ratios << classes[class_index].first << " ";
                for(int k=0; k<6; k++)
                {
                    for(int t=k+1; t<6; t++)
                    {
                        index_hist++;

                        //if(isnan(std::atan((lkl[k]-lkl[t])/3.)*180/M_PI/90))cout << "nan value ecountered " << lkl[k] << " " << lkl[t] << endl;
                        dump_lkl_ratios << std::atan((lkl[k]-lkl[t])/3.)*180/M_PI/90 << " ";
                        //all the likelihood ratios
                        likelihood_ratios[index_hist]=(std::atan((lkl[k]-lkl[t])/3.)*180/M_PI/90);
                        
                        //filling the istogram
                        likelihood_ratio[index_hist]->Fill(std::atan((lkl[k]-lkl[t])/3.)*180/M_PI/90);
                    }
                }  
                double dep_energy=0;
                for(int hit=0; hit<dat.track.dE->size(); hit++)
                {
                    if(dat.track.rr->at(hit)>5.)continue;
                    dep_energy = dep_energy + dat.track.dE->at(hit)*dat.track.pitch->at(hit);
                }
                dump_lkl_ratios << dep_energy << endl;
                //dump_lkl_ratios << endl;

                //SELECTION
                int reco_class=reco_classification(likelihood_ratios);
                confusion_matrix->Fill(true_class,reco_class);

                //storing likelihood ratios
                if(par=="muon")
                {
                    if(subclass==0)v_likelihood_ratios_mu_class0.push_back(likelihood_ratios);
                    else v_likelihood_ratios_mu_class1.push_back(likelihood_ratios);
                }
                if(par=="proton")
                {
                    if(subclass==0)v_likelihood_ratios_pro_class2.push_back(likelihood_ratios);
                    else v_likelihood_ratios_pro_class3.push_back(likelihood_ratios);
                }
                if(par=="pion")
                {
                    if(subclass==0)v_likelihood_ratios_pi_class4.push_back(likelihood_ratios);
                    else v_likelihood_ratios_pi_class5.push_back(likelihood_ratios);
                }  
            }//cycle over all tracks

            outfile->cd();
            d->cd();
            for(int histo=0; histo<15; histo++)
            {
                likelihood_ratio[histo]->Scale(1./likelihood_ratio[histo]->Integral());
                likelihood_ratio[histo]->Write(0,TObject::kOverwrite);
            }
        
            cout << tracks_with_no_valid_hits << " " << par << " of subclass " << subclass << " with no valid hits in COLL plane" << endl;

        }//cycle over subclasses
 
    }//cycle on particle type

    //confusion matrix
    for(int bin=1; bin<=6; bin++)
    {
        double total=0;
        for(int binbin=1; binbin<=6; binbin++)
        {
            total = total + confusion_matrix->GetBinContent(bin, binbin);
        }
        for(int binbin=1; binbin<=6; binbin++)
        {
            confusion_matrix->SetBinContent(bin,binbin,confusion_matrix->GetBinContent(bin, binbin)/total);
        }
    }
    outfile->cd();
    confusion_matrix->Write(0,TObject::kOverwrite);

    //making roc curves
    //std::array<TGraph*,4> roc_coll = getROCcurvesLR(20000, v_likelihood_ratio_coll_mu, v_likelihood_ratio_coll_pro, v_likelihood_ratio_coll_pi, "coll");
    //std::array<TGraph*,4> roc_coll_ind1 = getROCcurvesLR(20000,v_likelihood_ratio_coll_ind1_mu, v_likelihood_ratio_coll_ind1_pro,v_likelihood_ratio_coll_ind1_pi, "coll_ind1");
    //std::array<TGraph*,4> roc_coll_ind1_ind2 = getROCcurvesLR(20000,v_likelihood_ratio_coll_ind1_ind2_mu,v_likelihood_ratio_coll_ind1_ind2_pro,v_likelihood_ratio_coll_ind1_ind2_pi,"coll_ind1_ind2");

    //d->cd();
    //std::array<std::array<TGraph*,4>,3> rocs ={roc_coll,roc_coll_ind1,roc_coll_ind1_ind2};
    //sfor(int i=0; i<3; i++){ for(int j=0; j<4; j++) rocs[i][j]->Write(0,TObject::kOverwrite);}
    
}