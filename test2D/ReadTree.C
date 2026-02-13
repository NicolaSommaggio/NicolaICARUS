
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

    double trackscore;
 
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

    if(inputstringFile=="1d"){datafile="/exp/icarus/data/users/nsommagg/test1d.root";}
    if(inputstringFile=="2d"){datafile="/exp/icarus/data/users/nsommagg/test2d.root";}
    if(inputstringFile=="2d_DNN_NOPulseTrains"){datafile="/exp/icarus/data/users/nsommagg/test2d_DNN_NOPulseTrains.root";}
    if(inputstringFile=="2d_DNN_PulseTrains"){datafile="/exp/icarus/data/users/nsommagg/test2d_DNN_PulseTrains.root";}
    if(inputstringFile=="1d_YZ"){datafile="/exp/icarus/data/users/nsommagg/test1d_YZ.root";}
    if(inputstringFile=="2d_YZ"){datafile="/exp/icarus/data/users/nsommagg/test2d_YZ.root";}
    if(inputstringFile=="2d_DNN_NOPulseTrains_YZ"){datafile="/exp/icarus/data/users/nsommagg/test2d_DNN_NOPulseTrains_YZ.root";}
    if(inputstringFile=="2d_DNN_PulseTrains_YZ"){datafile="/exp/icarus/data/users/nsommagg/test2d_DNN_PulseTrains_YZ.root";}

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
        tree->SetBranchAddress("trackscore_mu", &dat.track.trackscore);

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
    }

    if(inputstringParticle=="proton")
    {
        tree->SetBranchAddress("vertice_reco", &dat.vertex.reco_vertex);
        tree->SetBranchAddress("bestplane_pro", &dat.track.bestplane);
        tree->SetBranchAddress("chi2_pro_coll", &dat.track.chi2_coll);
        tree->SetBranchAddress("chi2_pro_ind1", &dat.track.chi2_ind1);
        tree->SetBranchAddress("chi2_pro_ind2", &dat.track.chi2_ind2);
        tree->SetBranchAddress("trackscore_pro", &dat.track.trackscore);

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
    }

    if(inputstringParticle=="pion")
    {
        tree->SetBranchAddress("vertice_reco", &dat.vertex.reco_vertex);
        tree->SetBranchAddress("bestplane_pi", &dat.track.bestplane);
        tree->SetBranchAddress("chi2_pi_coll", &dat.track.chi2_coll);
        tree->SetBranchAddress("chi2_pi_ind1", &dat.track.chi2_ind1);
        tree->SetBranchAddress("chi2_pi_ind2", &dat.track.chi2_ind2);
        tree->SetBranchAddress("trackscore_pi", &dat.track.trackscore);

        tree->SetBranchAddress("dQdx_pi", &dat.track.dQdx);
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
