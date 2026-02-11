#include "/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/dEdx/ReadTree.C"

void make_prob_densities()
{
    gROOT->ProcessLine("make_file()");
    gROOT->ProcessLine("hprobDensities()");
    gROOT->ProcessLine(".q");
}

void make_prob_densitiesRM()
{
    gROOT->ProcessLine("make_file_RM()");
    gROOT->ProcessLine("hprobDensities_RM()");
    gROOT->ProcessLine(".q");
}

void make_file()
{
    TFile * f = new TFile("HISTO_prob_densities_6_classes.root","RECREATE");
    TDirectory * d_mu_class0 = (TDirectory*)f->mkdir("muon_class0");
    TDirectory * d_mu_class1 = (TDirectory*)f->mkdir("muon_class1");
    TDirectory * d_pro_class2 = (TDirectory*)f->mkdir("proton_class2");
    TDirectory * d_pro_class3 = (TDirectory*)f->mkdir("proton_class3");
    TDirectory * d_pi_class4 = (TDirectory*)f->mkdir("pion_class4");
    TDirectory * d_pi_class5 = (TDirectory*)f->mkdir("pion_class5");

    std::array<TDirectory*,6> dirs={d_mu_class0, d_mu_class1, d_pro_class2, d_pro_class3, d_pi_class4, d_pi_class5};
    for(const auto &dir : dirs)
    {
        for(double rr=1.5; rr<=25; rr+=0.5)
        {
            TDirectory *rrdir = (TDirectory*)dir->mkdir(Form("%.1f", rr));
        }
    }
    f->Close();

}

void make_file_RM()
{
    TFile * f = new TFile("HISTO_prob_densities_RM_6_classes.root","RECREATE");
    TDirectory * d_mu_class0 = (TDirectory*)f->mkdir("muon_class0");
    TDirectory * d_mu_class1 = (TDirectory*)f->mkdir("muon_class1");
    TDirectory * d_pro_class2 = (TDirectory*)f->mkdir("proton_class2");
    TDirectory * d_pro_class3 = (TDirectory*)f->mkdir("proton_class3");
    TDirectory * d_pi_class4 = (TDirectory*)f->mkdir("pion_class4");
    TDirectory * d_pi_class5 = (TDirectory*)f->mkdir("pion_class5");

    std::array<TDirectory*,6> dirs={d_mu_class0, d_mu_class1, d_pro_class2, d_pro_class3, d_pi_class4, d_pi_class5};
    for(const auto &dir : dirs)
    {
        for(double rr=1.5; rr<=25; rr+=0.5)
        {
            TDirectory *rrdir = (TDirectory*)dir->mkdir(Form("%.1f", rr));
        }
    }
    f->Close();

}

int true_selection(EventsData dat, std::string par)
{
    if(dat.track.rr->size()==0)return -1;

    //controllo che vertice reco e true siano abbastanza vicini per eliminare la componente di cosmici
    TVector3 vertex_true;
    vertex_true.SetXYZ(dat.vertex.true_vertex->at(0), dat.vertex.true_vertex->at(1), dat.vertex.true_vertex->at(2));
    TVector3 vertex_reco;
    vertex_reco.SetXYZ(dat.vertex.reco_vertex->at(0), dat.vertex.reco_vertex->at(1), dat.vertex.reco_vertex->at(2));
    if((vertex_true-vertex_reco).Mag()>100.)return -1;

    //compute distance between the endpoints
    TVector3 end_true;
    end_true.SetXYZ(dat.track.end_true->at(0), dat.track.end_true->at(1), dat.track.end_true->at(2));
    TVector3 end_hit;
    double endx = dat.track.hitx->at((int)dat.track.hitx->size()-1);
    double endy = dat.track.hity->at((int)dat.track.hity->size()-1);
    double endz = dat.track.hitz->at((int)dat.track.hitz->size()-1);
    end_hit.SetXYZ(endx,endy,endz);
    double end_distance = (end_hit-end_true).Mag();

    //looking at the energy match
    double daughter_electron_energy_match = -1;
    double daughter_proton_energy_match = -1;
    double daughter_pion_energy_match = -1;

    for(int match=0; match<dat.track.pdg_match->size(); match++)
    {       
        if(std::abs(dat.track.pdg_match->at(match))==2212 && dat.track.is_daughter->at(match)==true)
        {
            if(dat.track.energy_match->at(match)>daughter_proton_energy_match)daughter_proton_energy_match = dat.track.energy_match->at(match);
        }
        if(std::abs(dat.track.pdg_match->at(match))==211 && dat.track.is_daughter->at(match)==true)
        {
            if(dat.track.energy_match->at(match)>daughter_pion_energy_match)daughter_pion_energy_match = dat.track.energy_match->at(match);
        }
        if(std::abs(dat.track.pdg_match->at(match))==11 && dat.track.is_daughter->at(match)==true)
        {
            if(dat.track.energy_match->at(match)>daughter_electron_energy_match)daughter_electron_energy_match = dat.track.energy_match->at(match);
        }
    }   

    int true_class=-1;
    //TRUE CLASSIFICATION
    if(par=="muon")
    {
        if(dat.track.hit_completeness>=0.5 && dat.track.hit_purity>=0.5)
        {
            if(daughter_electron_energy_match<=0.01)
            {
                if(end_distance<=3.)true_class = 0;       
                else true_class = 1;
            }
            else true_class=6;
        } 
        else true_class=-1;
    }
    if(par=="proton")
    {
        if(dat.track.hit_completeness>=0.3 && dat.track.hit_purity>=0.3)
        {
            if
            ( 
                (dat.track.end_process==54 && end_distance <= 5.) ||
                (dat.track.end_process!=54 && (daughter_proton_energy_match>0.055 || daughter_pion_energy_match>0.055))
            ) true_class = 2;
            else true_class = 3;
        }
        else true_class=-1;
    }
    if(par=="pion")
    {
        if(dat.track.hit_completeness>=0.3 && dat.track.hit_purity>=0.3)
        {
            if(( dat.track.end_process==3 || dat.track.end_process==45 ) && end_distance <= 1.5) true_class = 4;
            else true_class = 5;
        }
        else true_class=-1;
    }

    return true_class;
}

void hprobDensities()
{
    std::vector<std::pair<std::string,int>> classes = {{"class0",0}, {"class1",1}, {"class2",2}, {"class3",3}, {"class4",4}, {"class5",5}};
    std::array<std::string,3> particles={"muon", "proton", "pion"};
    TFile *f = TFile::Open("HISTO_prob_densities_6_classes.root", "UPDATE");
    for(const auto& par : particles)
    {
        std::array<std::pair<std::string,int>,2> subclasses;
        if(par=="muon"){subclasses[0] = classes[0]; subclasses[1] = classes[1];}
        if(par=="proton"){subclasses[0] = classes[2]; subclasses[1] = classes[3];}
        if(par=="pion"){subclasses[0] = classes[4]; subclasses[1] = classes[5];}

        int used_tracks=0;
        for(const auto &subclass : subclasses)
        {
            cout << "processing " << par << " " << subclass.first << endl;
            EventsData dat = load_data("mc2d_general_contained", par);

            TH2D *dedx_vs_range = new TH2D("dedx_vs_range","",300,0,30,300,0,30);
    
            for(int track=0; track<dat.tree->GetEntries(); track++)
            {
                dat.tree->GetEntry(track);
                if(true_selection(dat,par)!=subclass.second)continue;

                for(int hit=0; hit<dat.track.rr->size(); hit++)
                {dedx_vs_range->Fill(dat.track.rr->at(hit),dat.track.dE->at(hit));}
                
            }

            cout << dat.tree->GetEntries() << " " << par << " tracks" << endl;
            TDirectory *pardir = (TDirectory*)f->Get(Form("%s_%s",par.c_str(),subclass.first.c_str()));
        
            for(double rr=1.5; rr<=25; rr+=0.5)
            {
                cout << "processing rr " << rr << endl;
                TDirectory *rrdir = (TDirectory*)pardir->Get(Form("%.1f", rr));
                rrdir->cd();

                int Nbins_coll = 100;

                std::vector<double> binlowedges_coll;
                for(int i=1; i<=Nbins_coll+1; i++)
                {
                    binlowedges_coll.push_back((i-1)*30./Nbins_coll);
                }
                bool keep_coll = true;

                TH1D *dEdx_coll;

                while(keep_coll)
                {
                    keep_coll = false;
                    dEdx_coll = new TH1D(Form("dEdx_coll_rr_%.1f", rr),"", Nbins_coll, binlowedges_coll.data());

                    for(int track=0; track<dat.tree->GetEntries(); track++)
                    {
                        //choosing only track with even index
                        //if(track % 2 !=0)continue;
                        dat.tree->GetEntry(track);
                        int true_class = true_selection(dat,par);

                        //CHOOSING THE CLASS
                        if(true_class!=subclass.second)continue;

                        //FILLING THE dE/dx HISTOGRAM
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
                        if(dEdx_coll->GetBinContent(bin+1) < 1) // bin+1 perché TH1D è 1-indexed
                        {
                            keep_coll = true;
                            binlowedges_coll.erase(binlowedges_coll.begin()+bin);
                            Nbins_coll--;
                        }
                    }
                    if(dEdx_coll->GetBinContent(1)<1){binlowedges_coll.erase(binlowedges_coll.begin()+1);}
                    if(keep_coll){delete dEdx_coll;}
                }
            
                dEdx_coll->Scale(1./dEdx_coll->Integral());
                dEdx_coll->Write(0,TObject::kOverwrite);
            }
            pardir->cd();
            dedx_vs_range->Write();
        }
    }
}

int true_selection_RM(EventsData dat, std::string par)
{
    if(dat.track.rr_rm_coll->size()==0 || dat.track.rr->size()==0)return -1;

    //controllo che vertice reco e true siano abbastanza vicini per eliminare la componente di cosmici
    TVector3 vertex_true;
    vertex_true.SetXYZ(dat.vertex.true_vertex->at(0), dat.vertex.true_vertex->at(1), dat.vertex.true_vertex->at(2));
    TVector3 vertex_reco;
    vertex_reco.SetXYZ(dat.vertex.reco_vertex->at(0), dat.vertex.reco_vertex->at(1), dat.vertex.reco_vertex->at(2));
    if((vertex_true-vertex_reco).Mag()>100.)return -1;

    //compute distance between the endpoints
    TVector3 end_true;
    end_true.SetXYZ(dat.track.end_true->at(0), dat.track.end_true->at(1), dat.track.end_true->at(2));
    TVector3 end_hit;
    double endx = dat.track.hitx->at((int)dat.track.hitx->size()-1);
    double endy = dat.track.hity->at((int)dat.track.hity->size()-1);
    double endz = dat.track.hitz->at((int)dat.track.hitz->size()-1);
    end_hit.SetXYZ(endx,endy,endz);
    double end_distance = (end_hit-end_true).Mag();

    //looking at the energy match
    double daughter_electron_energy_match = -1;
    double daughter_proton_energy_match = -1;
    double daughter_pion_energy_match = -1;

    for(int match=0; match<dat.track.pdg_match->size(); match++)
    {       
        if(std::abs(dat.track.pdg_match->at(match))==2212 && dat.track.is_daughter->at(match)==true)
        {
            if(dat.track.energy_match->at(match)>daughter_proton_energy_match)daughter_proton_energy_match = dat.track.energy_match->at(match);
        }
        if(std::abs(dat.track.pdg_match->at(match))==211 && dat.track.is_daughter->at(match)==true)
        {
            if(dat.track.energy_match->at(match)>daughter_pion_energy_match)daughter_pion_energy_match = dat.track.energy_match->at(match);
        }
        if(std::abs(dat.track.pdg_match->at(match))==11 && dat.track.is_daughter->at(match)==true)
        {
            if(dat.track.energy_match->at(match)>daughter_electron_energy_match)daughter_electron_energy_match = dat.track.energy_match->at(match);
        }
    }   

    int true_class=-1;
    //TRUE CLASSIFICATION
    if(par=="muon")
    {
        if(dat.track.hit_completeness>=0.5 && dat.track.hit_purity>=0.5)
        {
            if(daughter_electron_energy_match<=0.01)
            {
                if(end_distance<=3.)true_class = 0;       
                else true_class = 1;
            }
            else true_class=6;
        } 
        else true_class=-1;
    }
    if(par=="proton")
    {
        if(dat.track.hit_completeness>=0.3 && dat.track.hit_purity>=0.3)
        {
            if
            ( 
                (dat.track.end_process==54 && end_distance <= 5.) ||
                (dat.track.end_process!=54 && (daughter_proton_energy_match>0.055 || daughter_pion_energy_match>0.055))
            ) true_class = 2;
            else true_class = 3;
        }
        else true_class=-1;
    }
    if(par=="pion")
    {
        if(dat.track.hit_completeness>=0.3 && dat.track.hit_purity>=0.3)
        {
            if(( dat.track.end_process==3 || dat.track.end_process==45 ) && end_distance <= 1.5) true_class = 4;
            else true_class = 5;
        }
        else true_class=-1;
    }

    return true_class;
}

void hprobDensities_RM()
{
    std::vector<std::pair<std::string,int>> classes = {{"class0",0}, {"class1",1}, {"class2",2}, {"class3",3}, {"class4",4}, {"class5",5}};
    std::array<std::string,3> particles={"muon", "proton", "pion"};
    TFile *f = TFile::Open("HISTO_prob_densities_RM_6_classes.root", "UPDATE");
    for(const auto& par : particles)
    {
        std::array<std::pair<std::string,int>,2> subclasses;
        if(par=="muon"){subclasses[0] = classes[0]; subclasses[1] = classes[1];}
        if(par=="proton"){subclasses[0] = classes[2]; subclasses[1] = classes[3];}
        if(par=="pion"){subclasses[0] = classes[4]; subclasses[1] = classes[5];}

        for(const auto &subclass : subclasses)
        {
            cout << "processing " << par << " " << subclass.first << endl;
            EventsData dat = load_data("mc2d_general_contained", par);

            TH2D *rm_vs_range = new TH2D("rm_vs_range","",300,0,30,150,0,30);
    
            for(int track=0; track<dat.tree->GetEntries(); track++)
            {
                dat.tree->GetEntry(track);
                if(true_selection_RM(dat,par)!=subclass.second)continue;

                for(int hit=0; hit<dat.track.rm_coll->size(); hit++)
                {rm_vs_range->Fill(dat.track.rr_rm_coll->at(hit),dat.track.rm_coll->at(hit));}
                
            }

            cout << dat.tree->GetEntries() << " " << par << " tracks" << endl;
            TDirectory *pardir = (TDirectory*)f->Get(Form("%s_%s",par.c_str(),subclass.first.c_str()));
        
            for(double rr=1.5; rr<=25; rr+=0.5)
            {
                cout << "processing rr " << rr << endl;
                TDirectory *rrdir = (TDirectory*)pardir->Get(Form("%.1f", rr));
                rrdir->cd();

                int Nbins_coll = 100;

                std::vector<double> binlowedges_coll;
                for(int i=1; i<=Nbins_coll+1; i++)
                {
                    binlowedges_coll.push_back((i-1)*30./Nbins_coll);
                }
                bool keep_coll = true;

                TH1D *dEdx_coll;

                while(keep_coll)
                {
                    keep_coll = false;
                    dEdx_coll = new TH1D(Form("dEdx_coll_rr_%.1f", rr),"", Nbins_coll, binlowedges_coll.data());

                    for(int track=0; track<dat.tree->GetEntries(); track++)
                    {
                        //choosing only track with even index
                        //if(track % 2 !=0)continue;

                        dat.tree->GetEntry(track);
                        int true_class = true_selection_RM(dat, par);

                        //CHOOSING THE CLASS
                        if(true_class!=subclass.second)continue;

                        //FILLING THE dE/dx HISTOGRAM
                        for(int hit=0; hit<int(dat.track.rm_coll->size()); hit++)
                        {
                            if(dat.track.rr_rm_coll->at(hit)<rr && dat.track.rr_rm_coll->at(hit)>=rr-0.5)
                            {
                                dEdx_coll->Fill(dat.track.rm_coll->at(hit));
                            }   
                        }   
                    }   

                    for(int bin=Nbins_coll-1; bin>0; bin--)
                    {
                        if(dEdx_coll->GetBinContent(bin+1) < 1) // bin+1 perché TH1D è 1-indexed
                        {
                            keep_coll = true;
                            binlowedges_coll.erase(binlowedges_coll.begin()+bin);
                            Nbins_coll--;
                        }
                    }
                    if(dEdx_coll->GetBinContent(1)<1){binlowedges_coll.erase(binlowedges_coll.begin()+1);}
                    if(keep_coll){delete dEdx_coll;}
                }
            
                dEdx_coll->Scale(1./dEdx_coll->Integral());
                dEdx_coll->Write(0,TObject::kOverwrite);
            }
            pardir->cd();
            rm_vs_range->Write();
        }
    }
}