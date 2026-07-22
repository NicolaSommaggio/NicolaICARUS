#pragma once 

#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"

//#include "helper_selection_newPID.h"
#include "helper_selection_newPID_light.h"
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

void macro_selection_newPID(){

//MC 2D
//const std::string fdata = "production_mc_2025A_ICARUS_Overlays_BNB_MC_RUN2_summer_2025_v10_06_00_06p01_flatcaf";

//mc low stat
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlay_neutrino_stage1_flat_cafs_v10_06_00_04p04_concat.root";

//10% data run2
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/Data_prescaled_Run2_calibrated/*.root";

//large MC (30 % used to construct probability densities)
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026_[1]/*flat.caf.root";

//test file
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026_3/group_6664_1.flat.caf.root";

//offbeam
//const std::string fdata = "Icaruspro_2025_wcdnn_production_Reproc_Run2_SBN_v10_06_00_01p05_offbeambnbmajority_flatcaf_unblind";

//respun 
//const std::string fdata = "/pnfs/sbn/data/sbn_fd/poms_production/2024A_ICARUS_MC_Sys_NuCos/2024A_MC_Sys_NuCos_respunCV_2ndV/mc/reconstructed/icaruscode_v09_89_01_02p01/flatcaf/*/*/*.root";

//high LT
//const std::string fdata = "/pnfs/sbn/data/sbn_fd/poms_production/2024A_ICARUS_MC_Sys_NuCos/2024A_MC_Sys_NuCos_HighLT_2ndV_correct/mc/reconstructed/icaruscode_v09_89_01_01p03/flatcaf/*/*/*.root";

//null
//const std::string fdata = "/pnfs/sbn/data/sbn_fd/poms_production/2024A_ICARUS_MC_Sys_NuCos/2024A_MC_Sys_NuCos_NullVar_2ndV/mc/reconstructed/icaruscode_v09_89_01_01p03/flatcaf/*/*/*.root";

//path MC + cosmics 1D
//const std::string fdata = "/storage/gpfs_data/icarus/local/users/cfarnese/new_studies_Large_Prod_CNAF/Complete_MC_final/*run*.concatflat.caf.root";

//path 10 % data 1d
//const std::string fdata= "/storage/gpfs_data/icarus/local/users/cfarnese/new_studies_Large_Prod_CNAF/Prescaled_DATA_bnbmaj/*flat.caf.root";

//POT OFFBEAM = 2.11e+20
//POT MC CLASSIC = 1.5342e+20 
//POT MC GRAY = 4.88297e+20
//POT MC GRAY PDF = 1.22871e+20 
//POT MC GRAY TRAINING = 1.22419e+20 
//POT MC GRAY SELECTION = 2.43007e+20 

//19197 muon rising
//6510 muon mip
//61268 proton rising
//16810 proton interacting
//3245 pion rising
//8499 pion interacting

//POT MC LOW STAT = 1.64645e+19 --> TRAINING
//POT MC CLASSIC 1/3 = 5.01111e+19 --> PROB DISTRO
//POT DATA_PRESCALED = 1.97615e+19 / 1.99201e+19

//--> POT DATA NEW
//1.90381e+19 POT from hdr differs from 1.91916e+19 POT from the TotalPOT histogram!
//1.90381e+19 POT over 5111954 readouts


//SELECTION FOLDER 
//const std::string fdata = "/exp/icarus/data/users/nsommagg/SELECTION_FOLDER/*.root";

// 1 file SELECTION FOLDER
//const std::string fdata = "/exp/icarus/data/users/nsommagg/SELECTION_FOLDER/84005943_0_out1.flat.caf.root";

//ICARUS MC
//const std::string fdata = "/pnfs/icarus/scratch/users/gputnam/Ar23+_iterE/ICARUSSpringMC/*/*flat.caf.root";

//ICARUS Offbeam
//const std::string fdata = "Icaruspro_2025_wcdnn_production_Reproc_Run2_SBN_v10_06_00_01p05_offbeambnbmajority_flatcaf_unblind";

//DATI
const std::string fdata = "Icaruspro_2025_wcdnn_production_Reproc_Run2_SBN_2_v10_06_00_06p03_bnbmajority_flatcaf_prescaled";

//MC per SBND
//const std::string fdata = "/pnfs/sbn/scratch/users/sungbino/sbnd/v10_06_00_09/mc_1e20";

//CNAF
//TFile* file = TFile::Open("/storage/gpfs_data/icarus/local/users/marterop/sbnana_v09_78_06/mc_test/dEdxrestemplates.root");
//auto dedx_range_pro = (TProfile*)file->Get("dedx_range_pro");
//auto dedx_range_ka  = (TProfile*)file->Get("dedx_range_ka");
//auto dedx_range_pi  = (TProfile*)file->Get("dedx_range_pi");
//auto dedx_range_mu  = (TProfile*)file->Get("dedx_range_mu");

//FNAL
TFile* file = TFile::Open("/exp/icarus/data/users/nsommagg/RefCurvesChi2.root");
auto dedx_range_pro = (TProfile*)file->Get("dedx_range_pro");
auto dedx_range_ka  = (TProfile*)file->Get("dedx_range_ka");
auto dedx_range_pi  = (TProfile*)file->Get("dedx_range_pi");
auto dedx_range_mu  = (TProfile*)file->Get("dedx_range_mu");


SpectrumLoader loader(fdata);         

const Binning kBinz = Binning::Simple(300,0,30);

//Spectrum s1("", kBinz, loader, selection_newPID ,kNoSpillCut );
//Spectrum s1("", kBinz, loader, DATAlikelihood ,kCRTPMTNeutrino );
//Spectrum s1("", kBinz, loader, kStitch ,kNoSpillCut );
//Spectrum s1("", kBinz, loader, selection ,kNoSpillCut );
//Spectrum s1("", kBinz, loader, fchi2dump ,kNoSpillCut );
//Spectrum s1("", kBinz, loader, fdebug1mu0p0pi ,kNoSpillCut );
//Spectrum s1("", kBinz, loader, fdump_pred_proba ,kNoSpillCut );
//Spectrum s1("", kBinz, loader, files_check ,kNoSpillCut );
//Spectrum s1("", kBinz, loader, dump_BDT_vars ,kNoSpillCut );
Spectrum s1("", kBinz, loader, dedx_var ,kNoSpillCut );

//muons.clear();
//pions.clear();
//protons.clear();

loader.Go();

double factor = s1.POT();  

TH1D* h1 = s1.ToTH1(factor);

TFile *tree_outfile = new TFile("varMC_DATA.root","RECREATE");
TTree * tree = new TTree("tree","tree");

_slice thislice;
tree -> Branch("slice",&thislice);

for(const auto &slc : _slices)
{
    thislice = slc;
    tree -> Fill();
}

tree_outfile -> cd();
tree -> Write(0,TObject::kOverwrite);
tree_outfile -> Close();


/*
cout << nmuons << " " << npions << " " << nprotons << endl;

TFile *outfile = new TFile("pred_proba_tot.root","RECREATE");

TTree *muons_tree = new TTree("muons","");

TTree *pions_tree = new TTree("pions","");

TTree *protons_tree = new TTree("protons","");

int run_mu;
int evt_mu;
int slice_index_mu;
int bestplane_mu;
int nhit_mu;
int pdg_mu;
double length_mu;
double depE_mu;
int predicted_class_mu;
std::array<double,3> dir_mu;
std::vector<double> v_pred_proba_mu;
std::vector<double> rr_mu;
std::vector<double> dedx_mu;

int run_pro;
int evt_pro;
int slice_index_pro;
int bestplane_pro;
int nhit_pro;
int pdg_pro;
double length_pro;
double depE_pro;
int predicted_class_pro;
std::array<double,3> dir_pro;
std::vector<double> v_pred_proba_pro;
std::vector<double> rr_pro;
std::vector<double> dedx_pro;

int run_pi;
int evt_pi;
int slice_index_pi;
int bestplane_pi;
int nhit_pi;
int pdg_pi;
double length_pi;
double depE_pi;
int predicted_class_pi;
std::array<double,3> dir_pi;
std::vector<double> v_pred_proba_pi;
std::vector<double> rr_pi;
std::vector<double> dedx_pi;

muons_tree->Branch("run",&run_mu);
muons_tree->Branch("evt",&evt_mu);
muons_tree->Branch("slice_index",&slice_index_mu);
muons_tree->Branch("bestplane",&bestplane_mu);
muons_tree->Branch("nhit",&nhit_mu);
muons_tree->Branch("pdg",&pdg_mu);
muons_tree->Branch("length",&length_mu);
muons_tree->Branch("depE",&depE_mu);
muons_tree->Branch("dir",&dir_mu);
muons_tree->Branch("v_pred_proba",&v_pred_proba_mu);
muons_tree->Branch("rr",&rr_mu);
muons_tree->Branch("dedx",&dedx_mu);

pions_tree->Branch("run",&run_pi);
pions_tree->Branch("evt",&evt_pi);
pions_tree->Branch("slice_index",&slice_index_pi);
pions_tree->Branch("bestplane",&bestplane_pi);
pions_tree->Branch("nhit",&nhit_pi);
pions_tree->Branch("pdg",&pdg_pi);
pions_tree->Branch("length",&length_pi);
pions_tree->Branch("depE",&depE_pi);
pions_tree->Branch("dir",&dir_pi);
pions_tree->Branch("v_pred_proba",&v_pred_proba_pi);
pions_tree->Branch("rr",&rr_pi);
pions_tree->Branch("dedx",&dedx_pi);

protons_tree->Branch("run",&run_pro);
protons_tree->Branch("evt",&evt_pro);  
protons_tree->Branch("slice_index",&slice_index_pro);
protons_tree->Branch("bestplane",&bestplane_pro);
protons_tree->Branch("nhit",&nhit_pro);
protons_tree->Branch("pdg",&pdg_pro);
protons_tree->Branch("length",&length_pro);
protons_tree->Branch("depE",&depE_pro);
protons_tree->Branch("dir",&dir_pro);
protons_tree->Branch("v_pred_proba",&v_pred_proba_pro);
protons_tree->Branch("rr",&rr_pro);
protons_tree->Branch("dedx",&dedx_pro);

for(int i=0; i<(int)muons.size(); i++)
{
    run_mu = muons[i].run;
    evt_mu = muons[i].evt;
    slice_index_mu = muons[i].slice_index;
    bestplane_mu = muons[i].bestplane;
    nhit_mu = muons[i].nhit;
    pdg_mu = muons[i].pdg;
    length_mu = muons[i].length;
    depE_mu = muons[i].depE;
    dir_mu = muons[i].dir;
    v_pred_proba_mu = muons[i].v_pred_proba;
    rr_mu = muons[i].rr;
    dedx_mu = muons[i].dedx;

    muons_tree->Fill();
}

for(int i=0; i<(int)pions.size(); i++)
{
    run_pi = pions[i].run;
    evt_pi = pions[i].evt;
    slice_index_pi = pions[i].slice_index;
    bestplane_pi = pions[i].bestplane;
    nhit_pi = pions[i].nhit;
    pdg_pi = pions[i].pdg;
    length_pi = pions[i].length;
    depE_pi = pions[i].depE;
    dir_pi = pions[i].dir;
    v_pred_proba_pi = pions[i].v_pred_proba;
    rr_pi = pions[i].rr;
    dedx_pi = pions[i].dedx;

    pions_tree->Fill();
}

for(int i=0; i<(int)protons.size(); i++)
{
    run_pro = protons[i].run;
    evt_pro = protons[i].evt;
    slice_index_pro = protons[i].slice_index;
    bestplane_pro = protons[i].bestplane;
    nhit_pro = protons[i].nhit;
    pdg_pro = protons[i].pdg;
    length_pro = protons[i].length;
    depE_pro = protons[i].depE;
    dir_pro = protons[i].dir;
    v_pred_proba_pro = protons[i].v_pred_proba;
    rr_pro = protons[i].rr;
    dedx_pro = protons[i].dedx;

    protons_tree->Fill();
}

muons_tree->Write(0,TObject::kOverwrite);
pions_tree->Write(0,TObject::kOverwrite);
protons_tree->Write(0,TObject::kOverwrite);
outfile->Close();

*/

//counts how many 1muNp selected

/* //HERE 
double selected_1muNp = 0;
double selected_1muNp_correct = 0;

for(int i = 0; i<(int)slices_reco_class.size(); i++)
{
    if(is1muNp(slices_reco_class[i]))
    {
        selected_1muNp++;
        if(is1muNp_true(slices_true_class[i]))selected_1muNp_correct++;
    }
}

cout << tot_1muNp << " tot true 1muNp" << endl;
cout << selected_1muNp << " tot selected 1muNp" << endl;
if(ismc)
{
    cout << selected_1muNp_correct << " correct" << endl;
    cout << "efficiency " << selected_1muNp_correct/tot_1muNp << endl;
    cout << "purity " << selected_1muNp_correct/selected_1muNp << endl;

    //cout << endl << "********" << endl << endl;

    //cout << tot_1mu0p0pi_slices << " 1mu0p0pi slices " << tot_NON_protons_1mu0p0pi_slices << " with no protons" << endl;
    //cout << tot_NON_protons_1mu0p0pi_protons << " misclassified protons in 1mu0p0pi slices" << endl;
}
*/

//stitching only with daughters

/*
TH1D *h_end_distance_stitched = new TH1D("end_distance_stitched","",200,0,200);
TH1D *h_end_distance = new TH1D("end_distance","",200,0,200);
TH2D *h_len_reco_stitched_vs_len_true = new TH2D("len_reco_stitched_vs_len_true","",50,0,500,50,0,500);
TH2D *h_len_reco_vs_len_true = new TH2D("len_reco_vs_len_true","",50,0,500,50,0,500);
TH1D *h_len_diff_stitched = new TH1D("len_diff_stitched","",100,-5,5);
TH1D *h_len_diff = new TH1D("len_diff","",100,-5,5);

int tot_correctly_stitched = 0;
int tot_stiched = 0;
for(int i=0; i<(int)len_reco_stitched.size(); i++)
{
    double end_distance_stitched = std::sqrt(std::pow(end_reco_stitched[i][0]-end_true[i][0],2)+std::pow(end_reco_stitched[i][1]-end_true[i][1],2)+std::pow(end_reco_stitched[i][2]-end_true[i][2],2));
    double end_distance = std::sqrt(std::pow(end_reco_before_stitching[i][0]-end_true[i][0],2)+std::pow(end_reco_before_stitching[i][1]-end_true[i][1],2)+std::pow(end_reco_before_stitching[i][2]-end_true[i][2],2));

    h_end_distance_stitched->Fill(end_distance_stitched);
    h_end_distance->Fill(end_distance);
    h_len_reco_stitched_vs_len_true->Fill(len_reco_stitched[i],len_true[i]);
    h_len_reco_vs_len_true->Fill(len_reco_before_stitching[i],len_true[i]);
    h_len_diff_stitched->Fill((len_reco_stitched[i]-len_true[i])/len_true[i]);
    h_len_diff->Fill((len_reco_before_stitching[i]-len_true[i])/len_true[i]);

    if(correctly_stitched[i]==true){tot_correctly_stitched++;}
    tot_stiched++;
}

TFile * stitching_file = new TFile("stitching_file_new.root","RECREATE");
h_len_reco_stitched_vs_len_true->Write();
h_len_reco_vs_len_true->Write();
h_end_distance_stitched->Scale(1./h_end_distance_stitched->Integral());
h_end_distance->Scale(1./h_end_distance->Integral());
h_len_diff_stitched->Scale(1./h_len_diff_stitched->Integral());
h_len_diff->Scale(1./h_len_diff->Integral());
h_end_distance_stitched->Write();
h_end_distance->Write();
h_len_diff_stitched->Write();
h_len_diff->Write();

cout << tot_correctly_stitched << " correctly stitched out of " << tot_stiched << " stitched tracks" << endl;

*/
}
