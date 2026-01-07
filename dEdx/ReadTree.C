
//#include <unordered_map>
//#include "daughtersInfo.h"
#pragma link C++ class std::unordered_map<int,int>+;


struct tracksInfo
{
    std::vector<double> *dE=nullptr;
    std::vector<double> *rr=nullptr;
    std::vector<double> *mult=nullptr;
    std::vector<double> *width=nullptr;
    std::vector<double> *pitch=nullptr;
    std::vector<double> *Eint=nullptr;
    std::vector<double> *rr_invertito=nullptr;
    std::vector<double> *start_true=nullptr;
    std::vector<double> *start_reco=nullptr;
    std::vector<double> *end_true=nullptr;
    std::vector<double> *end_reco=nullptr;
    std::vector<double> *genMomentum=nullptr;
    std::vector<double> *integral=nullptr;
    std::vector<double> *sumadc=nullptr;
    std::vector<double> *hitx=nullptr;
    std::vector<double> *hity=nullptr;
    std::vector<double> *hitz=nullptr;
    std::vector<double> *dQdx=nullptr;
    std::vector<double> *wire=nullptr;

    std::vector<double> *dEdx_ind1=nullptr;
    std::vector<double> *dEdx_ind2=nullptr;
    std::vector<double> *rr_ind1=nullptr;
    std::vector<double> *rr_ind2=nullptr;

    std::vector<double> *chi2_coll= nullptr;
    std::vector<double> *chi2_ind1= nullptr;
    std::vector<double> *chi2_ind2= nullptr;

    int bestplane;
    int end_process;
    double len_reco;
    double len_true; 
    double dirx;
    double diry;
    double dirz;  
    //daughters_info 
    unsigned int ndaughters_reco;
    //std::vector<daughtersInfo> *daughters=nullptr;
 
};

struct vertexInfo {
    std::vector<double>* true_vertex = nullptr;
    std::vector<double>* reco_vertex = nullptr;
};

struct RunEvtInfo {
    unsigned int run;
    unsigned int subrun;
    unsigned int evt;
    bool ismc;
    int which_slice;
};

struct EventsData {
    tracksInfo track;
    vertexInfo vertex;
    RunEvtInfo EvtInfo;
    TTree* tree= nullptr;
    TFile* file= nullptr;
};

EventsData load_data(const std::string& inputstringFile, const std::string& inputstringParticle, const std::string& Np="Np", const std::string& fullRR="full", bool corrected=false, bool true1muNp=true, const std::string& option="") {

    
    EventsData dat;

    std::string datafile;
    std::string treename;

    //if(inputstringFile=="mc" && Np=="Np" && fullRR=="25only" && corrected==false && true1muNp==true && option=="" ){datafile="/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/datiMC_1muNp.root";};
    //if(inputstringFile=="dati" && Np=="Np" && fullRR=="25only" ){datafile="/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/datiDAT_1muNp.root";};
    //if(inputstringFile=="mc" && Np=="1p" && fullRR=="25only" && corrected==false && true1muNp==true && option==""){datafile="/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/datiMC_1mu1p.root";};
    //if(inputstringFile=="dati" && Np=="1p" && fullRR=="25only"){datafile="/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/datiDAT_1mu1p.root";};
    //if(inputstringFile=="nc" && Np=="1p" && fullRR=="25only" && corrected==false && true1muNp==true && option==""){datafile="/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/datiNC_1mu1p.root";};
    if(inputstringFile=="mc" && fullRR=="full" && corrected==false && true1muNp==true && option==""){datafile="/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/datafiles/datiMC_1muNp_fullRR.root";};
    //if(inputstringFile=="mc" && fullRR=="full" && corrected==false && true1muNp==true && option=="puro"){datafile="/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/datiMC_1muNp_puro_NUOVO.root";};
    if(inputstringFile=="dati" && fullRR=="full"){datafile="/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/datafiles/datiDAT_1muNp_fullRR.root";};
    if(inputstringFile=="mc" && fullRR=="full" && corrected==false && true1muNp==true && option=="withCOS"){datafile="/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/datafiles/datiMC_withCOS_1muNp_full.root";}
    if(inputstringFile=="mc" && fullRR=="full" && corrected==false && true1muNp==true && option=="YZvar"){datafile="/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/datafiles/datiMC_YZvar_1muNp_full.root";}

    //if(inputstringFile=="mc" /*&& fullRR=="25only"*/ && corrected==true && true1muNp==true && option==""){datafile="/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/datiMC_1muNp_chi2MOD.root";} 
    //if(inputstringFile=="mc" /*&& fullRR=="25only"*/ && corrected==true && true1muNp==true && option=="withCOS"){datafile="/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/datiMC_withCOS_1muNp_chi2MOD.root";}
    //if(inputstringFile=="mc" /*&& fullRR=="25only"*/ && corrected==true && true1muNp==true && option=="YZvar"){datafile="/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/datiMC_YZvar_1muNp_chi2MOD.root";}
    //if(inputstringFile=="mc" /*&& fullRR=="25only"*/ && corrected==true && true1muNp==false && option==""){datafile="/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/datiMC_1muNp_chi2MOD_notTrue.root";}
    //if(inputstringFile=="mc" /*&& fullRR=="25only"*/ && corrected==true && true1muNp==false && option=="withCOS"){datafile="/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/datiMC_withCOS_chi2MOD_notTrue.root";}
    //if(inputstringFile=="mc" /*&& fullRR=="25only"*/ && corrected==true && true1muNp==false && option=="YZvar"){datafile="/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/datiMC_YZvar_1muNp_chi2MOD_notTrue.root";}

    //if(inputstringFile=="mc" && corrected==false && true1muNp==false && option =="YZvar"){datafile="/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/datiMC_YZvar_notTrue.root";}

    if(inputstringFile=="dati2d"){datafile="/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/datafiles/dati100per100_nuovi_caf_2d.root";}
    if(inputstringFile=="dati2d10"){datafile="/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/datafiles/dati10percento_nuovi_caf_2d.root";}
    if(inputstringFile=="Full9435_1d"){datafile="/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/datafiles/fullRun9435_1D_new.root";}
    if(inputstringFile=="Full9435_2d"){datafile="/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/datafiles/fullRun9435_2D.root";}
    if(inputstringFile=="data_1d_run2_full_v09_89"){datafile="/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/datafiles/new100per100Run2.root";}
    if(inputstringFile=="dati10"){datafile="/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/datafiles/10_percent_data1d.root";}
    if(inputstringFile=="Full9435_2d_bestplane"){datafile="/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/datafiles/BestPlane_fullRun9435_2D.root";}
    if(inputstringFile=="Full9435_2d_Coll"){datafile="/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/datafiles/Coll_fullRun9435_2D.root";}
    if(inputstringFile=="Full9435_2d_Ind1"){datafile="/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/datafiles/Ind1_fullRun9435_2D.root";}
    if(inputstringFile=="Full9435_2d_Ind2"){datafile="/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/datafiles/Ind2_fullRun9435_2D.root";}
    if(inputstringFile=="mc2d"){datafile="/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/datafiles/mc2d.root";}
    if(inputstringFile=="mc2d_general"){datafile="/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/datafiles/mc2d_general.root";}
    if(inputstringFile=="mc2dprova"){datafile="/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/datafiles/mc2dprova.root";}
    if(inputstringFile=="mc2d_general_contained"){datafile="/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/datafiles/mc2d_general_contained.root";}

    cout << "Using " << datafile << " data sample" << endl;

    TFile *DataFile= TFile::Open(datafile.c_str());    

    if(inputstringParticle=="muon"){treename="DATAtreeMU";}
    if(inputstringParticle=="proton"){treename="DATAtreePRO";}
    if(inputstringParticle=="pion"){treename="DATAtreePI";}
    TTree *tree = (TTree*)DataFile->Get(treename.c_str());

    if(inputstringParticle=="muon")
    {

        tree->SetBranchAddress("bestplane_mu", &dat.track.bestplane);
        tree->SetBranchAddress("chi2_mu_coll", &dat.track.chi2_coll);
        tree->SetBranchAddress("chi2_mu_ind1", &dat.track.chi2_ind1);
        tree->SetBranchAddress("chi2_mu_ind2", &dat.track.chi2_ind2);

        tree->SetBranchAddress("wire_mu", &dat.track.wire);

        tree->SetBranchAddress("dE_mu_ind1", &dat.track.dEdx_ind1);
        tree->SetBranchAddress("dE_mu_ind2", &dat.track.dEdx_ind2);
        tree->SetBranchAddress("rr_mu_ind1", &dat.track.rr_ind1);
        tree->SetBranchAddress("rr_mu_ind2", &dat.track.rr_ind2);

        tree->SetBranchAddress("dE_mu", &dat.track.dE);
        tree->SetBranchAddress("rr_mu", &dat.track.rr);
        tree->SetBranchAddress("mult_mu", &dat.track.mult);
        tree->SetBranchAddress("width_mu", &dat.track.width);
        tree->SetBranchAddress("pitch_mu", &dat.track.pitch);
        tree->SetBranchAddress("Eint_mu", &dat.track.Eint);
        tree->SetBranchAddress("rr_mu_invertito", &dat.track.rr_invertito);
        
        tree->SetBranchAddress("integral_mu", &dat.track.integral);
        tree->SetBranchAddress("sumadc_mu", &dat.track.sumadc);
        tree->SetBranchAddress("hitxMU", &dat.track.hitx);
        tree->SetBranchAddress("hityMU", &dat.track.hity);
        tree->SetBranchAddress("hitzMU", &dat.track.hitz);

        tree->SetBranchAddress("dQdx", &dat.track.dQdx);

        tree->SetBranchAddress("vertice_reco", &dat.vertex.reco_vertex);
        tree->SetBranchAddress("end_reco_mu", &dat.track.end_reco);
        tree->SetBranchAddress("start_reco_mu", &dat.track.start_reco);
        tree->SetBranchAddress("len_reco_mu", &dat.track.len_reco);

        tree->SetBranchAddress("dirx", &dat.track.dirx);
        tree->SetBranchAddress("diry", &dat.track.diry);
        tree->SetBranchAddress("dirz", &dat.track.dirz);

        //daughters
        tree->SetBranchAddress("ndaughters_mu",&dat.track.ndaughters_reco);
        //tree->SetBranchAddress("mu_daughters_info", &dat.track.daughters);

        //general informations

        tree->SetBranchAddress("run", &dat.EvtInfo.run );
        tree->SetBranchAddress("subrun", &dat.EvtInfo.subrun);
        tree->SetBranchAddress("evt", &dat.EvtInfo.evt);
        tree->SetBranchAddress("ismc", &dat.EvtInfo.ismc);
        tree->SetBranchAddress("which_slice", &dat.EvtInfo.which_slice);
    }

    if(inputstringParticle=="proton")
    {
        tree->SetBranchAddress("vertice_reco", &dat.vertex.reco_vertex);
        tree->SetBranchAddress("bestplane_pro", &dat.track.bestplane);
        tree->SetBranchAddress("chi2_pro_coll", &dat.track.chi2_coll);
        tree->SetBranchAddress("chi2_pro_ind1", &dat.track.chi2_ind1);
        tree->SetBranchAddress("chi2_pro_ind2", &dat.track.chi2_ind2);

        tree->SetBranchAddress("dE_pro_ind1", &dat.track.dEdx_ind1);
        tree->SetBranchAddress("dE_pro_ind2", &dat.track.dEdx_ind2);
        tree->SetBranchAddress("rr_pro_ind1", &dat.track.rr_ind1);
        tree->SetBranchAddress("rr_pro_ind2", &dat.track.rr_ind2);

        tree->SetBranchAddress("dE_pro", &dat.track.dE);
        tree->SetBranchAddress("rr_pro", &dat.track.rr);
        tree->SetBranchAddress("mult_pro", &dat.track.mult);
        tree->SetBranchAddress("width_pro", &dat.track.width);
        tree->SetBranchAddress("pitch_pro", &dat.track.pitch);
        tree->SetBranchAddress("Eint_pro", &dat.track.Eint);
        tree->SetBranchAddress("rr_pro_invertito", &dat.track.rr_invertito);

        tree->SetBranchAddress("integral_pro", &dat.track.integral);
        tree->SetBranchAddress("sumadc_pro", &dat.track.sumadc);

        tree->SetBranchAddress("hitxPRO", &dat.track.hitx);
        tree->SetBranchAddress("hityPRO", &dat.track.hity);
        tree->SetBranchAddress("hitzPRO", &dat.track.hitz);

        tree->SetBranchAddress("dQdx", &dat.track.dQdx);

        //tree->SetBranchAddress("vertice_reco", &dat.vertex.reco_vertex);
        tree->SetBranchAddress("end_reco_pro", &dat.track.end_reco);
        tree->SetBranchAddress("start_reco_pro", &dat.track.start_reco);
        tree->SetBranchAddress("len_reco_pro", &dat.track.len_reco);

        tree->SetBranchAddress("dirx", &dat.track.dirx);
        tree->SetBranchAddress("diry", &dat.track.diry);
        tree->SetBranchAddress("dirz", &dat.track.dirz);

        //tree->SetBranchAddress("ndaughters_pro", &dat.track.ndaughters_reco);

        //daughters
        tree->SetBranchAddress("ndaughters_pro",&dat.track.ndaughters_reco);
        //tree->SetBranchAddress("pro_daughters_info", &dat.track.daughters);

        //general informations

        tree->SetBranchAddress("run", &dat.EvtInfo.run );
        tree->SetBranchAddress("subrun", &dat.EvtInfo.subrun);
        tree->SetBranchAddress("evt", &dat.EvtInfo.evt);
        tree->SetBranchAddress("ismc", &dat.EvtInfo.ismc);
        tree->SetBranchAddress("which_slice", &dat.EvtInfo.which_slice);
    }

    if(inputstringParticle=="pion")
    {
        tree->SetBranchAddress("vertice_reco", &dat.vertex.reco_vertex);
        tree->SetBranchAddress("bestplane_pi", &dat.track.bestplane);
        tree->SetBranchAddress("chi2_pi_coll", &dat.track.chi2_coll);
        tree->SetBranchAddress("chi2_pi_ind1", &dat.track.chi2_ind1);
        tree->SetBranchAddress("chi2_pi_ind2", &dat.track.chi2_ind2);

        tree->SetBranchAddress("dE_pi_ind1", &dat.track.dEdx_ind1);
        tree->SetBranchAddress("dE_pi_ind2", &dat.track.dEdx_ind2);
        tree->SetBranchAddress("rr_pi_ind1", &dat.track.rr_ind1);
        tree->SetBranchAddress("rr_pi_ind2", &dat.track.rr_ind2);

        tree->SetBranchAddress("dE_pi", &dat.track.dE);
        tree->SetBranchAddress("rr_pi", &dat.track.rr);

        tree->SetBranchAddress("hitxPI", &dat.track.hitx);
        tree->SetBranchAddress("hityPI", &dat.track.hity);
        tree->SetBranchAddress("hitzPI", &dat.track.hitz);

        //tree->SetBranchAddress("vertice_reco", &dat.vertex.reco_vertex);
        tree->SetBranchAddress("end_reco_pi", &dat.track.end_reco);
        tree->SetBranchAddress("start_reco_pi", &dat.track.start_reco);
        tree->SetBranchAddress("len_reco_pi", &dat.track.len_reco);

        tree->SetBranchAddress("dirx", &dat.track.dirx);
        tree->SetBranchAddress("diry", &dat.track.diry);
        tree->SetBranchAddress("dirz", &dat.track.dirz);

        //general informations

        tree->SetBranchAddress("run", &dat.EvtInfo.run );
        tree->SetBranchAddress("subrun", &dat.EvtInfo.subrun);
        tree->SetBranchAddress("evt", &dat.EvtInfo.evt);
        tree->SetBranchAddress("ismc", &dat.EvtInfo.ismc);
        tree->SetBranchAddress("which_slice", &dat.EvtInfo.which_slice);
    }

    if ((inputstringFile == "mc" || inputstringFile == "mc2d" || inputstringFile == "mc2d_general" || inputstringFile == "mc2dprova" || inputstringFile == "mc2d_general_contained") && (inputstringParticle=="muon")) 
    {
        tree->SetBranchAddress("vertice_true", &dat.vertex.true_vertex);
        tree->SetBranchAddress("end_true_mu", &dat.track.end_true);
        tree->SetBranchAddress("start_true_mu", &dat.track.start_true);
        tree->SetBranchAddress("len_true_mu", &dat.track.len_true);
        tree->SetBranchAddress("gen_momentum_mu", &dat.track.genMomentum);
        tree->SetBranchAddress("end_process_mu",&dat.track.end_process);
    }

    if((inputstringFile == "mc" || inputstringFile == "mc2d" || inputstringFile == "mc2d_general" || inputstringFile == "mc2dprova" || inputstringFile == "mc2d_general_contained") && (inputstringParticle=="proton"))
    {
        tree->SetBranchAddress("vertice_true", &dat.vertex.true_vertex);
        tree->SetBranchAddress("end_true_pro", &dat.track.end_true);
        tree->SetBranchAddress("start_true_pro", &dat.track.start_true);
        tree->SetBranchAddress("len_true_pro", &dat.track.len_true);
        tree->SetBranchAddress("gen_momentum_pro", &dat.track.genMomentum);
        tree->SetBranchAddress("end_process_pro", &dat.track.end_process);
    }

    if((inputstringFile == "mc" || inputstringFile == "mc2d" || inputstringFile == "mc2d_general" || inputstringFile == "mc2dprova" || inputstringFile == "mc2d_general_contained") && (inputstringParticle=="pion"))
    {
        tree->SetBranchAddress("vertice_true", &dat.vertex.true_vertex);
        tree->SetBranchAddress("end_true_pi", &dat.track.end_true);
        tree->SetBranchAddress("start_true_pi", &dat.track.start_true);
        tree->SetBranchAddress("len_true_pi", &dat.track.len_true);
        tree->SetBranchAddress("gen_momentum_pi", &dat.track.genMomentum);
        tree->SetBranchAddress("end_process_pi", &dat.track.end_process);
    }

    dat.tree = tree;
    dat.file = DataFile;
    return dat;
    
}

double distanza(std::vector<double> *p1, std::vector<double> *p2)
{
    return sqrt(pow(p1->at(0)-p2->at(0),2)+pow(p1->at(1)-p2->at(1),2)+pow(p1->at(2)-p2->at(2),2)); 
}

double mediana(std::vector<double> dummy)
{
    std::sort(dummy.begin(), dummy.end());
    double mediana=0;
    int idx=0;
    if(int(dummy.size())%2==1)
    {
        idx=(dummy.size()-1)/2;
        mediana=dummy.at(idx);
    }
    if(dummy.size()%2==0)
    {
        idx =dummy.size()/2-1;
        mediana=( dummy.at(dummy.size()/2) + dummy.at(idx) )/2.;
    }
    return mediana;
}

void selecet1mu1p(EventsData &dat_pro)
{
    std::vector<int> IDvector;
    std::unordered_map<int,int> count;
    int N1muNp=0; 
    int Nptracks=0;

    for(int i=0; i<dat_pro.tree->GetEntries(); i++)
    {
        dat_pro.tree->GetEntry(i);
        IDvector.push_back(dat_pro.EvtInfo.which_slice);
    }

    for(int id : IDvector){count[id]++;}

    for(const auto& it : count)
    {
        int id = it.first;
        int counts= it.second;
        if(counts>1)
        {
            N1muNp+=1;
            Nptracks= Nptracks+counts-1;
        }
    } 
    std::cout << N1muNp << " 1muN(>1)p events" << endl;
    std::cout << Nptracks << endl;


}