#pragma once 

#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"


#include "helper_1muNp_true_particles.h"
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

//1d
const std::string fdata = "msotgia_v10_06_00_07_BNB_yzsim_wcdnn_pulseTrains_numu_caf";

//ofstream checkMichel("checkMichel.txt", std::ios::app);
//checkMichel << "processing file " << fdata << endl;

//2d
//const std::string fdata = "msotgia_v10_06_00_07_BNB_yzsim_wcdnn_pulseTrains_numu_caf";

SpectrumLoader loader(fdata);    

const Binning kBinz     = Binning::Simple(300,0,30);

Spectrum s1("", kBinz, loader, DataLoader ,kCRTPMTNeutrino );

loader.Go();

double factor = s1.POT();  

TH1D* h1 = s1.ToTH1(factor);


TH1D *hflashGateTime = new TH1D("flashGateTime", "flashGateTime", 400, -2,2);
for(int i=0; i<int(flashGateTime.size()); i++){hflashGateTime->Fill(flashGateTime[i]);}

TH1D *htriggerWithinGate = new TH1D("triggerWithinGate","triggerWithinGate", 400,-2,-2);
for(int i=0; i<int(triggerWithinGate.size()); i++){htriggerWithinGate->Fill(triggerWithinGate[i]);}

std::string filename;

filename = "/exp/icarus/data/users/nsommagg/test2d_DNN_PulseTrains_YZ.root";

TFile *f = new TFile(filename.c_str(), "RECREATE");
hflashGateTime->Write(0,TObject::kOverwrite);
htriggerWithinGate->Write(0,TObject::kOverwrite);

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
std::vector<double> dQdxPI;

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

double trkscoreMU;
double trkscorePRO;
double trkscorePI;

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
tree_mu->Branch("trackscore_mu", &trkscoreMU);

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
tree_pro->Branch("trackscore_pro", &trkscorePRO);

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
tree_pi->Branch("trackscore_pi", &trkscorePI);

tree_pi->Branch("dQdx_pi", &dQdxPI);
tree_pi->Branch("dE_pi", &vecdEdxPI);
tree_pi->Branch("rr_pi", &vecrrPI);
tree_pi->Branch("hitxPI", &hitxPRO);
tree_pi->Branch("hityPI", &hityPRO);
tree_pi->Branch("hitzPI", &hitzPRO);

tree_pi->Branch("dE_pi_ind1", &vec_dEdx_ind1_PI);
tree_pi->Branch("dE_pi_ind2", &vec_dEdx_ind2_PI);
tree_pi->Branch("rr_pi_ind1", &vec_rr_ind1_PI);
tree_pi->Branch("rr_pi_ind2", &vec_rr_ind2_PI);

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
//protons
tree_pro->Branch("vertice_true", &vertice_true_pro);
tree_pro->Branch("end_true_pro", &end_true_pro);
tree_pro->Branch("start_true_pro", &start_true_pro);
tree_pro->Branch("len_true_pro", &len_true_pro);
tree_pro->Branch("gen_momentum_pro", &genMomentumPRO);
tree_pro->Branch("end_process_pro", &end_proc_pro);
//pions
tree_pi->Branch("vertice_true", &vertice_true_pi);
tree_pi->Branch("end_true_pi", &end_true_pi);
tree_pi->Branch("start_true_pi", &start_true_pi);
tree_pi->Branch("len_true_pi", &len_true_pi);
tree_pi->Branch("end_process_pi", &end_proc_pi);
}

for(int i=0; i<int(CHIdeMU.size()); i++)
{
      //muons
      bestplane_mu = CHIbest_plane_mu[i];
      chi2_mu_coll = chi_quadro_coll_MU[i];
      chi2_mu_ind1 = chi_quadro_ind1_MU[i];
      chi2_mu_ind2 = chi_quadro_ind2_MU[i];

      wireMU = CHIwireMU[i];

      vec_dEdx_ind1_MU = dEdx_ind1_MU[i];
      //cout << "ind1 " << dEdx_ind1_MU[i].size() << endl;
      vec_dEdx_ind2_MU = dEdx_ind2_MU[i];
      //cout << "ind2 " << dEdx_ind2_MU[i].size() << endl;

      vec_rr_ind1_MU = rr_ind1_MU[i];
      vec_rr_ind2_MU = rr_ind2_MU[i];

      vec1mu = CHIdeMU[i]; 
      //cout << "coll " << CHIdeMU[i].size() << endl;
    
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
      //vertex
      vertice_true_mu = CHIvtxMCmu[i];   

      //muons
      end_true_mu =       CHIend3DmcMU[i];
      start_true_mu =    CHIstart3DmcMU[i]; 
      len_true_mu =       CHItlMCmu[i];
      genMomentumMU = CHIgenMomentumMU[i];
      end_proc_mu = end_process_mu[i];

   }

   tree_mu->Fill();

}

cout << "ok filling muon" << endl;

for(int i=0; i<int(CHIdePRO.size()); i++)
{

   //protons
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


   vec_dEdx_ind1_PRO = dEdx_ind1_PRO[i];
   vec_dEdx_ind2_PRO = dEdx_ind2_PRO[i];
   vec_rr_ind1_PRO = rr_ind1_PRO[i];
   vec_rr_ind2_PRO = rr_ind2_PRO[i];


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
      //protons
      vertice_true_pro = CHIvtxMCpro[i];
      end_true_pro =      CHIend3DmcPRO[i];
      start_true_pro =       CHIstart3DmcPRO[i];
      len_true_pro =      CHItlMCpro[i];
      genMomentumPRO = CHIgenMomentumPRO[i];
      end_proc_pro = end_process_pro[i];
   }

   tree_pro->Fill();

}

cout << "ok filling protons" << endl;

//pions
for(int i=0; i<int(CHIdePI.size()); i++)
{
   vertice_reco_pi = CHIvtxRECOpi[i]; 

   dQdxPI = CHIdQdxPI[i];
   vecdEdxPI = CHIdePI[i];         
   vecrrPI = CHIrrPIninv[i];
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

   vec_dEdx_ind1_PI = dEdx_ind1_PI[i];
   vec_dEdx_ind2_PI = dEdx_ind2_PI[i];
   vec_rr_ind1_PI = rr_ind1_PI[i];
   vec_rr_ind2_PI = rr_ind2_PI[i];

   //general infomations
   run_pi = vrun_pi[i];
   subrun_pi = vsubrun_pi[i];
   evt_pi = vevt_pi[i];
   ismc_pi = visMC_pi[i];

   if(is_mc)
   {
      //pions
      vertice_true_pi = CHIvtxMCpi[i];
      end_true_pi =      CHIend3DmcPI[i];
      start_true_pi =       CHIstart3DmcPI[i];
      len_true_pi =      CHItlMCpi[i];
      end_proc_pi = end_process_pi[i];
      //cout << end_process_pi[i] << endl;
   }



   tree_pi->Fill();

}

cout << "ok filling pions" << endl;

tree_mu->Write(0, TObject::kOverwrite);
tree_pro->Write(0, TObject::kOverwrite);
tree_pi->Write(0, TObject::kOverwrite);

f->Close();

}