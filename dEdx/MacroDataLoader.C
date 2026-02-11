#pragma once 

#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"


//#include "helper_stitch_simulz0.h"
//#include "helper_1muNp_puro.h"
#include "helper_1muNp.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTreeReader.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "THStack.h"
#include <TLine.h>
#include <fstream>
#include <sstream>
#include <numeric>
#include <algorithm>


using namespace ana;

void MacroDataLoader(){

  
//10 % DEI DATI
//const std::string fdata= "/storage/gpfs_data/icarus/local/users/cfarnese/new_studies_Large_Prod_CNAF/Prescaled_DATA_bnbmaj/*flat.caf.root";

//1 file
//const std::string fdata= "/storage/gpfs_data/icarus/local/users/cfarnese/new_studies_Large_Prod_CNAF/Prescaled_DATA_bnbmaj/run2-v09_84_00_01-202403-cnaf-run933X-run930X.concatflat.caf.root" ;

//100% DEI DATI
//const std::string fdata= "/storage/gpfs_data/icarus/local/users/cfarnese/test_creation_cafs_controlsample/prod-at-cnaf-final-100/results/*flat.caf.root";

//path MC senza cosmici 
//const std::string fdata = "/storage/gpfs_data/icarus/plain/user/cfarnese/test_genie_v0984_largestat_nuonly/concat_singleneu_*.flat.caf.root";

//1 file MCsenza cosmici
//const std::string fdata= "/storage/gpfs_data/icarus/plain/user/cfarnese/test_genie_v0984_largestat_nuonly/concat_singleneu_2.flat.caf.root";

//YZ variations 
//const std::string fdata = "/storage/gpfs_data/icarus/local/users/cfarnese/New_CAFS_Variations_Nov2024/YZ/CV_yz_FNAL_10new_*.flat.caf.root";

//path MC + cosmici 
//const std::string fdata = "/storage/gpfs_data/icarus/local/users/cfarnese/new_studies_Large_Prod_CNAF/Complete_MC_final/*run*.concatflat.caf.root";


// path nuovi caf 2d deconv 100%
//const std::string fdata= "/storage/gpfs_data/icarus/local/users/cfarnese/Run2_production_Oct2025_flatcafs_unblind_concat/*.root";

//path nuovi caf 2d 10%
//const std::string fdata= "/storage/gpfs_data/icarus/local/users/cfarnese/Run2_production_Oct2025_2D_flatcafs_prescaled_concat/*.root";

// path nuovi caf 2d deconv 1 file
//const std::string fdata= "/storage/gpfs_data/icarus/local/users/cfarnese/Run2_production_Oct2025_flatcafs_unblind_concat/run100xx_flatcaf_unblind.root";


//100 per 100 RUN 2 - run 9435 2D deconv.
//const std::string fdata= "/storage/gpfs_data/icarus/local/users/cfarnese/9435_spring_production_2025/9435_2025_FNAL/cafs_calibrated/*.root";

//100 per 100 RUN 2 - run 9435 1D deconv.
//const std::string fdata= "/storage/gpfs_data/icarus/local/users/cfarnese/9435_1muNp_icarusonly/9435_2025_FNAL_2/*.root";

// tutto run 2 1d versione 09 89
//const std::string fdata = "/storage/gpfs_data/icarus/plain/data/prod/run2-v09_89_01_01p03-202412-fnal/flatcafs_allevents/Run2*_v0989p03_Run*.flat.caf.root";

//MC 2D
const std::string fdata = "/storage/gpfs_data/icarus/local/users/cfarnese/production_mc_2025A_ICARUS_Overlays_BNB_MC_RUN2_summer_2025_v10_06_00_04p04/MC_overlay_neutrino_stage1_flat_cafs_v10_06_00_04p04_concat.root";

//new productions
//1d

TFile* file = TFile::Open("/storage/gpfs_data/icarus/local/users/marterop/sbnana_v09_78_06/mc_test/dEdxrestemplates.root");
auto dedx_range_pro = (TProfile*)file->Get("dedx_range_pro");
auto dedx_range_mu  = (TProfile*)file->Get("dedx_range_mu");
auto dedx_range_ka  = (TProfile*)file->Get("dedx_range_ka");
auto dedx_range_pi  = (TProfile*)file->Get("dedx_range_pi");

SpectrumLoader loader(fdata);       //CAF that I produced with all dedx vs rr   

const Binning kBinz = Binning::Simple(300,0,30);

Spectrum s1("", kBinz, loader, DataLoader ,kCRTPMTNeutrino );

loader.Go();

double factor = s1.POT();  

TH1D* h1 = s1.ToTH1(factor);

std::string filename;


//filename = "../datafiles/datiMC_1muNp_fullRR.root";
//filename = "../datafiles/datiDAT_1muNp_fullRR.root";
//filename = "../datafiles/dati100per100_nuovi_caf_2d.root";
//filename = "../datafiles/dati10percento_nuovi_caf_2d.root";

//filename = "../datafiles/BestPlane_fullRun9435_2D.root"; //ovviamente solo 2d
//filename = "../datafiles/Coll_fullRun9435_2D.root";
//filename = "../datafiles/Ind1_fullRun9435_2D.root";
//filename = "../datafiles/Ind2_fullRun9435_2D.root";

//filename = "../datafiles/fullRun9435_1D_new.root";

//filename = "../datafiles/10_percent_data1d.root";

//filename = "../datafiles/fullRun9435_2D.root";

//filename = "../datafiles/new100per100Run2.root";

//filename = "../datafiles/mc2d.root";

//filename = "../datafiles/mc2d_general.root";
filename = "../datafiles/mc2d_general_contained.root";

TFile *f = new TFile(filename.c_str(), "RECREATE");

TTree *tree_mu= new TTree("DATAtreeMU", "");
TTree *tree_pro = new TTree("DATAtreePRO", "");
TTree *tree_pi = new TTree("DATAtreePI","");

std::vector<double> vec1mu;        
std::vector<double> vec2mu;
std::vector<double> vec3mu; 
std::vector<double> vec4mu;
std::vector<int> vec5mu;    
std::vector<double> vec6mu;
std::vector<double> vec7mu;

std::vector<double> vec1pro;        
std::vector<double> vec2pro; 
std::vector<double> vec3pro; 
std::vector<double> vec4pro;
std::vector<int> vec5pro; 
std::vector<double> vec6pro;
std::vector<double> vec7pro;

std::vector<double> vecdEdxPI;        
std::vector<double> vecrrPI; 
std::vector<double> vecpitchPI;

std::vector<double> vec_dEdx_ind1_MU;
std::vector<double> vec_dEdx_ind2_MU;
std::vector<double> vec_rr_ind1_MU;
std::vector<double> vec_rr_ind2_MU;

std::vector<double> vec_dEdx_ind1_PRO;
std::vector<double> vec_dEdx_ind2_PRO;
std::vector<double> vec_rr_ind1_PRO;
std::vector<double> vec_rr_ind2_PRO;

std::vector<double> vec_dEdx_ind1_PI;
std::vector<double> vec_dEdx_ind2_PI;
std::vector<double> vec_rr_ind1_PI;
std::vector<double> vec_rr_ind2_PI;

std::vector<double> vec_rr_bestplane_MU;
std::vector<double> vec_dEdx_bestplane_MU;
std::vector<double> vec_rr_bestplane_PRO;
std::vector<double> vec_dEdx_bestplane_PRO;
std::vector<double> vec_rr_bestplane_PI;
std::vector<double> vec_dEdx_bestplane_PI;

std::vector<double> vec_pitch_bestplane_MU;
std::vector<double> vec_pitch_bestplane_PRO;
std::vector<double> vec_pitch_bestplane_PI;

std::vector<double> vec_hitx_bestplane_MU;
std::vector<double> vec_hity_bestplane_MU;
std::vector<double> vec_hitz_bestplane_MU;
std::vector<double> vec_hitx_bestplane_PRO;
std::vector<double> vec_hity_bestplane_PRO;
std::vector<double> vec_hitz_bestplane_PRO;
std::vector<double> vec_hitx_bestplane_PI;
std::vector<double> vec_hity_bestplane_PI;
std::vector<double> vec_hitz_bestplane_PI;

std::vector<double> chi2_mu_coll;
std::vector<double> chi2_mu_ind1;
std::vector<double> chi2_mu_ind2;

std::vector<double> chi2_pro_coll;
std::vector<double> chi2_pro_ind1;
std::vector<double> chi2_pro_ind2;

std::vector<double> chi2_pi_coll;
std::vector<double> chi2_pi_ind1;
std::vector<double> chi2_pi_ind2;

int bestplane_mu;
int bestplane_pro;
int bestplane_pi;

unsigned int run_mu;
unsigned int subrun_mu;
unsigned int evt_mu;
bool ismc_mu;

unsigned int run_pro;
unsigned int subrun_pro;
unsigned int evt_pro;
bool ismc_pro;

unsigned int run_pi;
unsigned int subrun_pi;
unsigned int evt_pi;
bool ismc_pi;

double dirx_mu;
double diry_mu;
double dirz_mu;

double dirx_pro;
double diry_pro;
double dirz_pro;

std::vector<double> vertice_true_mu;
std::vector<double> vertice_reco_mu;
std::vector<double> vertice_true_pro;
std::vector<double> vertice_reco_pro;
std::vector<double> vertice_true_pi;
std::vector<double> vertice_reco_pi;

std::vector<double> end_reco_mu;
std::vector<double> end_true_mu;
std::vector<double> start_reco_mu;
std::vector<double> start_true_mu;

std::vector<double> end_reco_pro;
std::vector<double> end_true_pro;
std::vector<double> start_reco_pro;
std::vector<double> start_true_pro;

std::vector<double> end_reco_pi;
std::vector<double> end_true_pi;
std::vector<double> start_reco_pi;
std::vector<double> start_true_pi;

std::vector<double> genMomentumMU;
std::vector<double> genMomentumPRO;

std::vector<double> integral_pro;
std::vector<double> integral_mu;
std::vector<double> sumadc_mu;
std::vector<double> sumadc_pro;

std::vector<double> hitxMU;
std::vector<double> hityMU;
std::vector<double> hitzMU;

std::vector<double> hitxPRO;
std::vector<double> hityPRO;
std::vector<double> hitzPRO;

std::vector<double> hitxPI;
std::vector<double> hityPI;
std::vector<double> hitzPI;

std::vector<double> dQdxMU;
std::vector<double> dQdxPRO;

std::vector<double> wireMU;

unsigned int ndaughters_mu_reco;
unsigned int ndaughters_pro_reco;

//std::vector<daughtersInfo> muons_d;
//std::vector<daughtersInfo> protons_d;

double len_true_mu;
double len_reco_mu;
double len_true_pro;
double len_reco_pro;
double len_true_pi;
double len_reco_pi;

int end_proc_mu, end_proc_pro, end_proc_pi;

std::vector<double> rr_rolling_median_mu_coll;
std::vector<double> rolling_median_mu_coll;
std::vector<double> rr_rolling_median_pro_coll;
std::vector<double> rolling_median_pro_coll;
std::vector<double> rr_rolling_median_pi_coll;
std::vector<double> rolling_median_pi_coll;

std::vector<double> rr_rolling_median_mu_bestplane;
std::vector<double> rolling_median_mu_bestplane;
std::vector<double> rr_rolling_median_pro_bestplane;
std::vector<double> rolling_median_pro_bestplane;
std::vector<double> rr_rolling_median_pi_bestplane;
std::vector<double> rolling_median_pi_bestplane;

std::vector<double> rr_rolling_median_mu_ind1;
std::vector<double> rolling_median_mu_ind1;
std::vector<double> rr_rolling_median_pro_ind1;
std::vector<double> rolling_median_pro_ind1;
std::vector<double> rr_rolling_median_pi_ind1;
std::vector<double> rolling_median_pi_ind1;

std::vector<double> rr_rolling_median_mu_ind2;
std::vector<double> rolling_median_mu_ind2;
std::vector<double> rr_rolling_median_pro_ind2;
std::vector<double> rolling_median_pro_ind2;
std::vector<double> rr_rolling_median_pi_ind2;
std::vector<double> rolling_median_pi_ind2;

bool has_Tsecondaries_mu;
bool has_Tsecondaries_pro;
bool has_Tsecondaries_pi;

double visE_mu;
double visE_pro;
double visE_pi;

double hitP_mu;
double hitP_pro;
double hitP_pi;
double hitC_mu;
double hitC_pro;
double hitC_pi;

double energyP_mu;
double energyP_pro;
double energyP_pi;
double energyC_mu;
double energyC_pro;
double energyC_pi;

std::vector<int> vec_pdg_matches_mu;
std::vector<double> vec_energy_matches_mu;
std::vector<double> vec_end_points_distance_matches_mu;
std::vector<int> vec_pdg_matches_pro;
std::vector<double> vec_energy_matches_pro;
std::vector<double> vec_end_points_distance_matches_pro;
std::vector<int> vec_pdg_matches_pi;
std::vector<double> vec_energy_matches_pi;
std::vector<double> vec_end_points_distance_matches_pi;

std::vector<double> vec_pitch_ind1_MU;
std::vector<double> vec_pitch_ind2_MU;
std::vector<double> vec_pitch_ind1_PRO;
std::vector<double> vec_pitch_ind2_PRO;
std::vector<double> vec_pitch_ind1_PI;
std::vector<double> vec_pitch_ind2_PI;

std::vector<bool> vec_is_daughter_mu;
std::vector<bool> vec_is_daughter_pro;
std::vector<bool> vec_is_daughter_pi;

int temp_start_process_mu;
int temp_start_process_pro;
int temp_start_process_pi;

double temp_energy_at_first_hit_mu;
double temp_energy_at_first_hit_pro;
double temp_energy_at_first_hit_pi;

double temp_energy_at_last_hit_mu;
double temp_energy_at_last_hit_pro;
double temp_energy_at_last_hit_pi;

//DEFINING BRANCHES

//general run evt info
tree_mu->Branch("run", &run_mu);
tree_mu->Branch("subrun", &subrun_mu);
tree_mu->Branch("evt", &evt_mu);
tree_mu->Branch("ismc", &ismc_mu);

tree_pro->Branch("run", &run_pro);
tree_pro->Branch("subrun", &subrun_pro);
tree_pro->Branch("evt", &evt_pro);
tree_pro->Branch("ismc", &ismc_pro);

tree_pi->Branch("run", &run_pi);
tree_pi->Branch("subrun", &subrun_pi);
tree_pi->Branch("evt", &evt_pi);
tree_pi->Branch("ismc", &ismc_pi);

//muons
//calo info
tree_mu->Branch("bestplane_mu", &bestplane_mu);
tree_mu->Branch("chi2_mu_coll", &chi2_mu_coll);
tree_mu->Branch("chi2_mu_ind1", &chi2_mu_ind1);
tree_mu->Branch("chi2_mu_ind2", &chi2_mu_ind2);

tree_mu->Branch("rr_rm_mu_coll", &rr_rolling_median_mu_coll);
tree_mu->Branch("rm_mu_coll", &rolling_median_mu_coll);
tree_mu->Branch("rr_rm_mu_ind1", &rr_rolling_median_mu_ind1);
tree_mu->Branch("rm_mu_ind1", &rolling_median_mu_ind1);
tree_mu->Branch("rr_rm_mu_ind2", &rr_rolling_median_mu_ind2);
tree_mu->Branch("rm_mu_ind2", &rolling_median_mu_ind2);
tree_mu->Branch("rr_rm_mu_bestplane", &rr_rolling_median_mu_bestplane);
tree_mu->Branch("rm_mu_bestplane", &rolling_median_mu_bestplane);

tree_mu->Branch("wire_mu", &wireMU);
tree_mu->Branch("dE_mu", &vec1mu);
tree_mu->Branch("Eint_mu", &vec2mu);
tree_mu->Branch("rr_mu", &vec3mu);
tree_mu->Branch("rr_mu_invertito", &vec4mu );
tree_mu->Branch("mult_mu", &vec5mu);
tree_mu->Branch("width_mu", &vec6mu);
tree_mu->Branch("pitch_mu", &vec7mu);
tree_mu->Branch("integral_mu", &integral_mu);
tree_mu->Branch("sumadc_mu", &sumadc_mu);
tree_mu->Branch("hitxMU", &hitxMU);
tree_mu->Branch("hityMU", &hityMU);
tree_mu->Branch("hitzMU", &hitzMU);
tree_mu->Branch("dQdx", &dQdxMU);

tree_mu->Branch("dE_mu_ind1", &vec_dEdx_ind1_MU);
tree_mu->Branch("dE_mu_ind2", &vec_dEdx_ind2_MU);
tree_mu->Branch("rr_mu_ind1", &vec_rr_ind1_MU);
tree_mu->Branch("rr_mu_ind2", &vec_rr_ind2_MU);

tree_mu->Branch("dE_mu_bestplane", &vec_dEdx_bestplane_MU);
tree_mu->Branch("rr_mu_bestplane", &vec_rr_bestplane_MU);
tree_mu->Branch("pitch_mu_bestplane", &vec_pitch_bestplane_MU);
tree_mu->Branch("hitx_mu_bestplane", &vec_hitx_bestplane_MU);
tree_mu->Branch("hity_mu_bestplane", &vec_hity_bestplane_MU);
tree_mu->Branch("hitz_mu_bestplane", &vec_hitz_bestplane_MU);

tree_mu->Branch("pitch_mu_ind1", &vec_pitch_ind1_MU);
tree_mu->Branch("pitch_mu_ind2", &vec_pitch_ind2_MU);

//track info
tree_mu->Branch("vertice_reco", &vertice_reco_mu);
tree_mu->Branch("end_reco_mu", &end_reco_mu);
tree_mu->Branch("start_reco_mu", &start_reco_mu);
tree_mu->Branch("len_reco_mu", &len_reco_mu);
tree_mu->Branch("dirx", &dirx_mu);
tree_mu->Branch("diry", &diry_mu);
tree_mu->Branch("dirz", &dirz_mu);

//daughters info
tree_mu->Branch("ndaughters_mu", &ndaughters_mu_reco);
//tree_mu->Branch("mu_daughters_info", &muons_d);

//protons
//calo info
tree_pro->Branch("bestplane_pro", &bestplane_pro);
tree_pro->Branch("chi2_pro_coll", &chi2_pro_coll);
tree_pro->Branch("chi2_pro_ind1", &chi2_pro_ind1);
tree_pro->Branch("chi2_pro_ind2", &chi2_pro_ind2);

tree_pro->Branch("rr_rm_pro_coll", &rr_rolling_median_pro_coll);
tree_pro->Branch("rm_pro_coll", &rolling_median_pro_coll);
tree_pro->Branch("rr_rm_pro_ind1", &rr_rolling_median_pro_ind1);
tree_pro->Branch("rm_pro_ind1", &rolling_median_pro_ind1);
tree_pro->Branch("rr_rm_pro_ind2", &rr_rolling_median_pro_ind2);
tree_pro->Branch("rm_pro_ind2", &rolling_median_pro_ind2);
tree_pro->Branch("rr_rm_pro_bestplane", &rr_rolling_median_pro_bestplane);
tree_pro->Branch("rm_pro_bestplane", &rolling_median_pro_bestplane);

tree_pro->Branch("dE_pro", &vec1pro);
tree_pro->Branch("Eint_pro", &vec2pro);
tree_pro->Branch("rr_pro", &vec3pro);
tree_pro->Branch("rr_pro_invertito", &vec4pro );
tree_pro->Branch("mult_pro", &vec5pro);
tree_pro->Branch("width_pro", &vec6pro);
tree_pro->Branch("pitch_pro", &vec7pro);
tree_pro->Branch("integral_pro", &integral_pro);
tree_pro->Branch("sumadc_pro", &sumadc_pro);
tree_pro->Branch("hitxPRO", &hitxPRO);
tree_pro->Branch("hityPRO", &hityPRO);
tree_pro->Branch("hitzPRO", &hitzPRO);
tree_pro->Branch("dQdx", &dQdxPRO);

tree_pro->Branch("dE_pro_ind1", &vec_dEdx_ind1_PRO);
tree_pro->Branch("dE_pro_ind2", &vec_dEdx_ind2_PRO);
tree_pro->Branch("rr_pro_ind1", &vec_rr_ind1_PRO);
tree_pro->Branch("rr_pro_ind2", &vec_rr_ind2_PRO);

tree_pro->Branch("dE_pro_bestplane", &vec_dEdx_bestplane_PRO);
tree_pro->Branch("rr_pro_bestplane", &vec_rr_bestplane_PRO);
tree_pro->Branch("pitch_pro_bestplane", &vec_pitch_bestplane_PRO);
tree_pro->Branch("hitx_pro_bestplane", &vec_hitx_bestplane_PRO);
tree_pro->Branch("hity_pro_bestplane", &vec_hity_bestplane_PRO);
tree_pro->Branch("hitz_pro_bestplane", &vec_hitz_bestplane_PRO);

tree_pro->Branch("pitch_pro_ind1", &vec_pitch_ind1_PRO);
tree_pro->Branch("pitch_pro_ind2", &vec_pitch_ind2_PRO);

//track info
tree_pro->Branch("vertice_reco", &vertice_reco_pro);
tree_pro->Branch("end_reco_pro", &end_reco_pro);
tree_pro->Branch("start_reco_pro", &start_reco_pro);
tree_pro->Branch("len_reco_pro", &len_reco_pro);
tree_pro->Branch("dirx", &dirx_pro);
tree_pro->Branch("diry", &diry_pro);
tree_pro->Branch("dirz", &dirz_pro);

//daughters info
tree_pro->Branch("ndaughters_pro", &ndaughters_pro_reco);
//tree_pro->Branch("pro_daughters_info", &protons_d);



//pions
//calo info
tree_pi->Branch("bestplane_pi", &bestplane_pi);
tree_pi->Branch("chi2_pi_coll", &chi2_pi_coll);
tree_pi->Branch("chi2_pi_ind1", &chi2_pi_ind1);
tree_pi->Branch("chi2_pi_ind2", &chi2_pi_ind2);

tree_pi->Branch("rr_rm_pi_coll", &rr_rolling_median_pi_coll);
tree_pi->Branch("rm_pi_coll", &rolling_median_pi_coll);
tree_pi->Branch("rr_rm_pi_ind1", &rr_rolling_median_pi_ind1);
tree_pi->Branch("rm_pi_ind1", &rolling_median_pi_ind1);
tree_pi->Branch("rr_rm_pi_ind2", &rr_rolling_median_pi_ind2);
tree_pi->Branch("rm_pi_ind2", &rolling_median_pi_ind2);
tree_pi->Branch("rr_rm_pi_bestplane", &rr_rolling_median_pi_bestplane);
tree_pi->Branch("rm_pi_bestplane", &rolling_median_pi_bestplane);

tree_pi->Branch("dE_pi", &vecdEdxPI);
tree_pi->Branch("rr_pi", &vecrrPI);
tree_pi->Branch("pitch_pi",&vecpitchPI);
tree_pi->Branch("hitxPI", &hitxPI);
tree_pi->Branch("hityPI", &hityPI);
tree_pi->Branch("hitzPI", &hitzPI);

tree_pi->Branch("dE_pi_ind1", &vec_dEdx_ind1_PI);
tree_pi->Branch("dE_pi_ind2", &vec_dEdx_ind2_PI);
tree_pi->Branch("rr_pi_ind1", &vec_rr_ind1_PI);
tree_pi->Branch("rr_pi_ind2", &vec_rr_ind2_PI);

tree_pi->Branch("dE_pi_bestplane", &vec_dEdx_bestplane_PI);
tree_pi->Branch("rr_pi_bestplane", &vec_rr_bestplane_PI);
tree_pi->Branch("pitch_pi_bestplane", &vec_pitch_bestplane_PI);
tree_pi->Branch("hitx_pi_bestplane", &vec_hitx_bestplane_PI);
tree_pi->Branch("hity_pi_bestplane", &vec_hity_bestplane_PI);
tree_pi->Branch("hitz_pi_bestplane", &vec_hitz_bestplane_PI);

tree_pi->Branch("pitch_pi_ind1", &vec_pitch_ind1_PI);
tree_pi->Branch("pitch_pi_ind2", &vec_pitch_ind2_PI);

//track info
tree_pi->Branch("vertice_reco", &vertice_reco_pi);
tree_pi->Branch("end_reco_pi", &end_reco_pi);
tree_pi->Branch("start_reco_pi", &start_reco_pi);
tree_pi->Branch("len_reco_pi", &len_reco_pi);


if(is_mc)
{
//muons
tree_mu->Branch("vertice_true", &vertice_true_mu);
tree_mu->Branch("end_true_mu", &end_true_mu);
tree_mu->Branch("start_true_mu", &start_true_mu);
tree_mu->Branch("len_true_mu", &len_true_mu);
tree_mu->Branch("gen_momentum_mu", &genMomentumMU);
tree_mu->Branch("end_process_mu", &end_proc_mu);
tree_mu->Branch("visE_mu",&visE_mu);
tree_mu->Branch("hit_purity_mu",&hitP_mu);
tree_mu->Branch("hit_completeness_mu",&hitC_mu);
tree_mu->Branch("energy_purity_mu",&energyP_mu);
tree_mu->Branch("energy_completeness_mu",&energyC_mu);
tree_mu->Branch("pdg_matches_mu",&vec_pdg_matches_mu);
tree_mu->Branch("energy_matches_mu",&vec_energy_matches_mu);
tree_mu->Branch("is_daughter_mu",&vec_is_daughter_mu);
tree_mu->Branch("end_points_distance_matches_mu",&vec_end_points_distance_matches_mu);
tree_mu->Branch("start_process_mu",&temp_start_process_mu);
tree_mu->Branch("energy_at_first_hit_mu",&temp_energy_at_first_hit_mu);
tree_mu->Branch("energy_at_last_hit_mu",&temp_energy_at_last_hit_mu);
tree_mu->Branch("has_true_secondaries_mu",&has_Tsecondaries_mu);

//protons
tree_pro->Branch("vertice_true", &vertice_true_pro);
tree_pro->Branch("end_true_pro", &end_true_pro);
tree_pro->Branch("start_true_pro", &start_true_pro);
tree_pro->Branch("len_true_pro", &len_true_pro);
tree_pro->Branch("gen_momentum_pro", &genMomentumPRO);
tree_pro->Branch("end_process_pro", &end_proc_pro);
tree_pro->Branch("visE_pro",&visE_pro);
tree_pro->Branch("hit_purity_pro",&hitP_pro);
tree_pro->Branch("hit_completeness_pro",&hitC_pro);
tree_pro->Branch("energy_purity_pro",&energyP_pro);
tree_pro->Branch("energy_completeness_pro",&energyC_pro);
tree_pro->Branch("pdg_matches_pro",&vec_pdg_matches_pro);
tree_pro->Branch("energy_matches_pro",&vec_energy_matches_pro);
tree_pro->Branch("is_daughter_pro",&vec_is_daughter_pro);
tree_pro->Branch("end_points_distance_matches_pro",&vec_end_points_distance_matches_pro);
tree_pro->Branch("has_true_secondaries_pro",&has_Tsecondaries_pro);
tree_pro->Branch("start_process_pro",&temp_start_process_pro);
tree_pro->Branch("energy_at_first_hit_pro",&temp_energy_at_first_hit_pro);
tree_pro->Branch("energy_at_last_hit_pro",&temp_energy_at_last_hit_pro);

//pions
tree_pi->Branch("vertice_true", &vertice_true_pi);
tree_pi->Branch("end_true_pi", &end_true_pi);
tree_pi->Branch("start_true_pi", &start_true_pi);
tree_pi->Branch("len_true_pi", &len_true_pi);
tree_pi->Branch("end_process_pi", &end_proc_pi);
tree_pi->Branch("visE_pi",&visE_pi);
tree_pi->Branch("hit_purity_pi",&hitP_pi);
tree_pi->Branch("hit_completeness_pi",&hitC_pi);
tree_pi->Branch("energy_purity_pi",&energyP_pi);
tree_pi->Branch("energy_completeness_pi",&energyC_pi);
tree_pi->Branch("pdg_matches_pi",&vec_pdg_matches_pi);
tree_pi->Branch("energy_matches_pi",&vec_energy_matches_pi);
tree_pi->Branch("is_daughter_pi",&vec_is_daughter_pi);
tree_pi->Branch("end_points_distance_matches_pi", &vec_end_points_distance_matches_pi);
tree_pi->Branch("has_true_secondaries_pi",&has_Tsecondaries_pi);
tree_pi->Branch("start_process_pi",&temp_start_process_pi);
tree_pi->Branch("energy_at_first_hit_pi",&temp_energy_at_first_hit_pi);
tree_pi->Branch("energy_at_last_hit_pi",&temp_energy_at_last_hit_pi);

}

//FILLING BRANCHES VARIABLES MUONS
for(int i=0; i<int(CHIdeMU.size()); i++)
{
   bestplane_mu = CHIbest_plane_mu[i];
   chi2_mu_coll = chi_quadro_coll_MU[i];
   chi2_mu_ind1 = chi_quadro_ind1_MU[i];
   chi2_mu_ind2 = chi_quadro_ind2_MU[i];

   rr_rolling_median_mu_coll = rr_rm_mu_coll[i];
   rolling_median_mu_coll = rm_mu_coll[i];
   rr_rolling_median_mu_ind1 = rr_rm_mu_ind1[i];
   rolling_median_mu_ind1 = rm_mu_ind1[i];
   rr_rolling_median_mu_ind2 = rr_rm_mu_ind2[i];
   rolling_median_mu_ind2 = rm_mu_ind2[i];
   rr_rolling_median_mu_bestplane = rr_rm_mu_bestplane[i];
   rolling_median_mu_bestplane = rm_mu_bestplane[i];

   wireMU = CHIwireMU[i];

   vec_dEdx_ind1_MU = dEdx_ind1_MU[i];
   vec_dEdx_ind2_MU = dEdx_ind2_MU[i];
   vec_rr_ind1_MU = rr_ind1_MU[i];
   vec_rr_ind2_MU = rr_ind2_MU[i];

   vec_dEdx_bestplane_MU = dEdx_bestplane_MU[i];
   vec_rr_bestplane_MU = rr_bestplane_MU[i];
   vec_pitch_bestplane_MU = pitch_bestplane_MU[i];
   vec_hitx_bestplane_MU = hitx_bestplane_MU[i];
   vec_hity_bestplane_MU = hity_bestplane_MU[i];
   vec_hitz_bestplane_MU = hitz_bestplane_MU[i];

   vec_pitch_ind1_MU = pitch_ind1_MU[i];
   vec_pitch_ind2_MU = pitch_ind2_MU[i];

   vec1mu = CHIdeMU[i]; 

   vec2mu = CHIEintMUallt[i];
   vec3mu = CHIrrMUninv[i];
   vec4mu = CHIrrMU[i];
   vec5mu = CHImultMU[i];
   vec6mu = CHIwidthMU[i];
   vec7mu = CHIpitchMU[i];
   integral_mu = CHIintegralMU[i];
   sumadc_mu = CHIsumadcMU[i];
   hitxMU = CHIxMU[i];
   hityMU = CHIyMU[i];
   hitzMU = CHIzMU[i];
   end_reco_mu =       CHIend3DrecoMU[i];
   start_reco_mu =   CHIstart3DrecoMU[i]; 
   len_reco_mu =       CHItlRECOmu[i]; 
   dQdxMU = CHIdQdxMU[i];

   dirx_mu = vdirx_mu[i];
   diry_mu = vdiry_mu[i];
   dirz_mu = vdirz_mu[i];

   //daughters
   //muons_d=muons_daughter[i];
   ndaughters_mu_reco=MUndaughters_reco[i];

   //vertex
   vertice_reco_mu = CHIvtxRECOmu[i]; 

   //event information
   run_mu = vrun_mu[i];
   subrun_mu = vsubrun_mu[i];
   evt_mu = vevt_mu[i];
   ismc_mu = visMC_mu[i];

   if(is_mc)
   {
      vertice_true_mu = CHIvtxMCmu[i];    
      end_true_mu =       CHIend3DmcMU[i];
      start_true_mu =    CHIstart3DmcMU[i]; 
      len_true_mu =       CHItlMCmu[i];
      genMomentumMU = CHIgenMomentumMU[i];
      end_proc_mu = end_process_mu[i];
      has_Tsecondaries_mu = has_true_secondaries_MU[i];
      visE_mu = visible_energy_mu[i];
      hitP_mu = hit_purity_mu[i];
      hitC_mu = hit_completeness_mu[i];
      energyP_mu = energy_purity_mu[i];
      energyC_mu = energy_completeness_mu[i];
      vec_pdg_matches_mu = pdg_matches_mu[i];
      vec_energy_matches_mu = energy_matches_mu[i];
      vec_is_daughter_mu = is_daughter_mu[i];
      vec_end_points_distance_matches_mu = end_points_distance_matches_mu[i];
      temp_start_process_mu = start_process_mu[i];
      temp_energy_at_first_hit_mu = energy_at_first_hit_mu[i];
      temp_energy_at_last_hit_mu = energy_at_last_hit_mu[i];
   }
   tree_mu->Fill();
}

cout << "ok filling muons" << endl;

//FILLING BRANCHES VARIABLES PROTONS
for(int i=0; i<int(CHIdePRO.size()); i++)
{
   vertice_reco_pro = CHIvtxRECOpro[i]; 

   vec1pro = CHIdePRO[i];         
   vec2pro = CHIEintPROallt[i];
   vec3pro = CHIrrPROninv[i];
   vec4pro = CHIrrPRO[i];
   vec5pro = CHImultPRO[i];
   vec6pro = CHIwidthPRO[i];
   vec7pro = CHIpitchPRO[i];
   end_reco_pro =      CHIend3DrecoPRO[i];
   start_reco_pro =       CHIstart3DrecoPRO[i];
   len_reco_pro =      CHItlRECOpro[i];
   hitxPRO = CHIxPRO[i];
   hityPRO = CHIyPRO[i];
   hitzPRO = CHIzPRO[i];
   sumadc_pro = CHIsumadcPRO[i];
   integral_pro = CHIintegralPRO[i];
   dQdxPRO = CHIdQdxPRO[i];

   dirx_pro = vdirx_pro[i];
   diry_pro = vdiry_pro[i];
   dirz_pro = vdirz_pro[i];

   bestplane_pro = CHIbest_plane_pro[i];
   chi2_pro_coll = chi_quadro_coll_PRO[i];
   chi2_pro_ind1 = chi_quadro_ind1_PRO[i];
   chi2_pro_ind2 = chi_quadro_ind2_PRO[i];

   rr_rolling_median_pro_coll = rr_rm_pro_coll[i];
   rolling_median_pro_coll = rm_pro_coll[i];
   rr_rolling_median_pro_ind1 = rr_rm_pro_ind1[i];
   rolling_median_pro_ind1 = rm_pro_ind1[i];
   rr_rolling_median_pro_ind2 = rr_rm_pro_ind2[i];
   rolling_median_pro_ind2 = rm_pro_ind2[i];
   rr_rolling_median_pro_bestplane = rr_rm_pro_bestplane[i];
   rolling_median_pro_bestplane = rm_pro_bestplane[i];

   vec_pitch_ind1_PRO = pitch_ind1_PRO[i];
   vec_pitch_ind2_PRO = pitch_ind2_PRO[i];

   vec_dEdx_ind1_PRO = dEdx_ind1_PRO[i];
   vec_dEdx_ind2_PRO = dEdx_ind2_PRO[i];
   vec_rr_ind1_PRO = rr_ind1_PRO[i];
   vec_rr_ind2_PRO = rr_ind2_PRO[i];

   vec_dEdx_bestplane_PRO = dEdx_bestplane_PRO[i];
   vec_rr_bestplane_PRO = rr_bestplane_PRO[i];
   vec_pitch_bestplane_PRO = pitch_bestplane_PRO[i];
   vec_hitx_bestplane_PRO = hitx_bestplane_PRO[i];
   vec_hity_bestplane_PRO = hity_bestplane_PRO[i];
   vec_hitz_bestplane_PRO = hitz_bestplane_PRO[i];

   //daughters
   //protons_d=protons_daughter[i];
   ndaughters_pro_reco=PROndaughters_reco[i];

   //general infomations
   run_pro = vrun_pro[i];
   subrun_pro = vsubrun_pro[i];
   evt_pro = vevt_pro[i];
   ismc_pro = visMC_pro[i];

   if(is_mc)
   {
      vertice_true_pro = CHIvtxMCpro[i];
      end_true_pro =      CHIend3DmcPRO[i];
      start_true_pro =       CHIstart3DmcPRO[i];
      len_true_pro =      CHItlMCpro[i];
      genMomentumPRO = CHIgenMomentumPRO[i];
      end_proc_pro = end_process_pro[i];
      has_Tsecondaries_pro = has_true_secondaries_PRO[i];
      visE_pro = visible_energy_pro[i];
      hitP_pro = hit_purity_pro[i];
      hitC_pro = hit_completeness_pro[i];
      energyP_pro = energy_purity_pro[i];
      energyC_pro = energy_completeness_pro[i];
      vec_pdg_matches_pro = pdg_matches_pro[i];
      vec_energy_matches_pro = energy_matches_pro[i];
      vec_is_daughter_pro = is_daughter_pro[i];
      vec_end_points_distance_matches_pro = end_points_distance_matches_pro[i];
      temp_start_process_pro = start_process_pro[i];
      temp_energy_at_first_hit_pro = energy_at_first_hit_pro[i];
      temp_energy_at_last_hit_pro = energy_at_last_hit_pro[i];
   }
   tree_pro->Fill();
}

cout << "ok filling protons" << endl;

//FILLING BRANCHES VARIABLES PIONS
for(int i=0; i<int(CHIdePI.size()); i++)
{
   vertice_reco_pi = CHIvtxRECOpi[i]; 

   vecdEdxPI = CHIdePI[i];         
   vecrrPI = CHIrrPIninv[i];
   vecpitchPI = CHIpitchPI[i];
   end_reco_pi =CHIend3DrecoPI[i];
   start_reco_pi =       CHIstart3DrecoPI[i];
   len_reco_pi =      CHItlRECOpi[i];
   hitxPI = CHIxPI[i];
   hityPI = CHIyPI[i];
   hitzPI = CHIzPI[i];

   bestplane_pi = CHIbest_plane_pi[i];
   chi2_pi_coll = chi_quadro_coll_PI[i];
   chi2_pi_ind1 = chi_quadro_ind1_PI[i];
   chi2_pi_ind2 = chi_quadro_ind2_PI[i];

   rr_rolling_median_pi_coll = rr_rm_pi_coll[i];
   rolling_median_pi_coll = rm_pi_coll[i];
   rr_rolling_median_pi_ind1 = rr_rm_pi_ind1[i];
   rolling_median_pi_ind1 = rm_pi_ind1[i];
   rr_rolling_median_pi_ind2 = rr_rm_pi_ind2[i];
   rolling_median_pi_ind2 = rm_pi_ind2[i];
   rr_rolling_median_pi_bestplane = rr_rm_pi_bestplane[i];
   rolling_median_pi_bestplane = rm_pi_bestplane[i];

   vec_pitch_ind1_PI = pitch_ind1_PI[i];
   vec_pitch_ind2_PI = pitch_ind2_PI[i];

   vec_dEdx_ind1_PI = dEdx_ind1_PI[i];
   vec_dEdx_ind2_PI = dEdx_ind2_PI[i];
   vec_rr_ind1_PI = rr_ind1_PI[i];
   vec_rr_ind2_PI = rr_ind2_PI[i];

   vec_dEdx_bestplane_PI = dEdx_ind2_PI[i];
   vec_rr_bestplane_PI = rr_ind1_PI[i];
   vec_pitch_bestplane_PI = pitch_bestplane_PI[i];
   vec_hitx_bestplane_PI = hitx_bestplane_PI[i];
   vec_hity_bestplane_PI = hity_bestplane_PI[i];
   vec_hitz_bestplane_PI = hitz_bestplane_PI[i];

   //general infomations
   run_pi = vrun_pi[i];
   subrun_pi = vsubrun_pi[i];
   evt_pi = vevt_pi[i];
   ismc_pi = visMC_pi[i];

   if(is_mc)
   {
      vertice_true_pi = CHIvtxMCpi[i];
      end_true_pi =      CHIend3DmcPI[i];
      start_true_pi =       CHIstart3DmcPI[i];
      len_true_pi =      CHItlMCpi[i];
      end_proc_pi = end_process_pi[i];
      has_Tsecondaries_pi = has_true_secondaries_PI[i];
      visE_pi = visible_energy_pi[i];
      hitP_pi = hit_purity_pi[i];
      hitC_pi = hit_completeness_pi[i];
      energyP_pi = energy_purity_pi[i];
      energyC_pi = energy_completeness_pi[i];
      vec_pdg_matches_pi = pdg_matches_pi[i];
      vec_energy_matches_pi = energy_matches_pi[i];
      vec_is_daughter_pi = is_daughter_pi[i];
      vec_end_points_distance_matches_pi = end_points_distance_matches_pi[i];
      temp_start_process_pi = start_process_pi[i];
      temp_energy_at_first_hit_pi = energy_at_first_hit_pi[i];
      temp_energy_at_last_hit_pi = energy_at_last_hit_pi[i];
   }
   tree_pi->Fill();
}

cout << "ok filling pions" << endl;


tree_mu->Write(0, TObject::kOverwrite);
tree_pro->Write(0, TObject::kOverwrite);
tree_pi->Write(0, TObject::kOverwrite);

dedx_range_mu->Write(0,TObject::kOverwrite);
dedx_range_pro->Write(0,TObject::kOverwrite);
dedx_range_pi->Write(0,TObject::kOverwrite);
dedx_range_ka->Write(0,TObject::kOverwrite);

f->Close();

}