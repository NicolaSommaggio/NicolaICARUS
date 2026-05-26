#pragma once 

#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"


//#include "helper_stitch_simulz0.h"
//#include "helper_1muNp_puro.h"
//#include "create_data_struct.h"
#include "create_data_struct_multiplane.h"
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

// carica la libreria

using namespace ana;

void MacroDataStruct(){
//MC 2D
//const std::string fdata = "production_mc_2025A_ICARUS_Overlays_BNB_MC_RUN2_summer_2025_v10_06_00_06p01_flatcaf";

//mc low stat
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlay_neutrino_stage1_flat_cafs_v10_06_00_04p04_concat.root";

//10 % statistics 1.47154e+19 POT
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026/group_6657.flat.caf.root";

//100% statistics concat
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026/*flat.caf.root";
//divided in 3
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026_1/*flat.caf.root";

//30 percent stat
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026_1/*flat.caf.root";

//1
//5.01111e+19 POT over 0 readouts
//tot 1muNp interactions: 6652

//2
//5.1204e+19 POT over 0 readouts
//tot 1muNp interactions: 6721

//3
//5.2105e+19 POT over 0 readouts
//tot 1muNp interactions: 6733


// 11 files concat
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026_3/group_6664.flat.caf.root";  // 13G 1.66188e+19 POT 
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026_2/group_6662.flat.caf.root";  // 12G 1.58954e+19 POT 
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026_1/group_6656_2.flat.caf.root";  // 11G 1.47245e+19 POT
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026_1/group_6657.flat.caf.root";  // 11G 1.47154e+19 POT
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026_2/group_6661.flat.caf.root";  // 11G 1.35995e+19 POT
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026_3/group_6663.flat.caf.root";  //9.2G 1.24961e+19 POT
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026_3/group_6660.flat.caf.root";  //9.1G 1.24326e+19 POT
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026_2/group_6658.flat.caf.root";  //8.2G 1.12881e+19 POT
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026_1/group_6656_1.flat.caf.root";  //7.8G 1.05475e+19 POT
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026_2/group_667X.flat.caf.root";  //7.7G 1.04211e+19 POT
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026_1/group_6656_3.flat.caf.root";  //7.5G 1.01238e+19 POT
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026_3/group_668X.flat.caf.root";  //4.2G 5.70136e+18 POT
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026_3/group_6659.flat.caf.root";  //3.6G 4.85613e+18 POT

//10% run2
//const std::string fdata = "Icaruspro_2025_wcdnn_production_Reproc_Run2_SBN_wCalib_v10_06_00_04p03_bnbmajority_flatcaf_prescaled";
//concat
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/Data_prescaled_Run2_calibrated/group_1X_data_prescaled.flat.caf.root";
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/Data_prescaled_Run2_calibrated/group_2X1234_data_prescaled.flat.caf.root";
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/Data_prescaled_Run2_calibrated/group_2X5678_data_prescaled.flat.caf.root"; //4.45726e+18 POT
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/Data_prescaled_Run2_calibrated/group_3X_data_prescaled.flat.caf.root";
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/Data_prescaled_Run2_calibrated/group_X2345X_data_prescaled.flat.caf.root";
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/Data_prescaled_Run2_calibrated/group_X6789X_data_prescaled.flat.caf.root";

//file che contengono RUN 9435
//const std::string fdata = "Icaruspro_2025_wcdnn_production_Reproc_Run2_SBN_wCalib_v10_06_00_04p03_bnbmajority_flatcaf_unblind";

//mpv from nu 2.36005e+19 POT
//const std::string fdata = "/exp/icarus/data/users/nsommagg/MPVMPR/nu_1/*.root";

//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026_[3]/*flat.caf.root";

const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026_3/group_6659.flat.caf.root";

//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026_4/group_6664_24.flat.caf.root";

SpectrumLoader loader(fdata);       //CAF that I produced with all dedx vs rr   

const Binning kBinz = Binning::Simple(300,0,30);

std::string filename;
filename = "/exp/icarus/data/users/nsommagg/DataStructDebug.root";
//filename = "/exp/icarus/data/users/nsommagg/data_struct_multiplane/run9435.root";
TFile *f = new TFile(filename.c_str(), "RECREATE");

tree = new TTree("tree","");

tree->Branch("slice",&thisslice);

cout << "using file: " << fdata << endl;
cout << "writing the data struct in: " << filename << endl;
//cout << "writing neutrino from spill in: " << txt_name << endl;

//TFile * f_prob_densities_coll_150 = TFile::Open("/exp/icarus/data/users/nsommagg/HISTO_prob_densities_6_classes_30percent_mediumbinning.root", "READ");
//TFile * f_prob_densities_ind1_150 = TFile::Open("/exp/icarus/data/users/nsommagg/HISTO_prob_densities_6_classes_30percent_mediumbinning_IND1.root", "READ");
//TFile * f_prob_densities_ind2_150 = TFile::Open("/exp/icarus/data/users/nsommagg/HISTO_prob_densities_6_classes_30percent_mediumbinning_IND2.root", "READ");

//TFile * f_prob_densities_coll = TFile::Open("/exp/icarus/data/users/nsommagg/HISTO_prob_densities_6_classes_30percent_stat.root", "READ");
//TFile * f_prob_densities_ind1 = TFile::Open("/exp/icarus/data/users/nsommagg/HISTO_prob_densities_6_classes_30percent_IND1.root", "READ");
//TFile * f_prob_densities_ind2 = TFile::Open("/exp/icarus/data/users/nsommagg/HISTO_prob_densities_6_classes_30percent_IND2.root", "READ");

//prob_d_coll_150 = load_prob_densities("coll",f_prob_densities_coll_150);
//prob_d_ind1_150 = load_prob_densities("ind1",f_prob_densities_ind1_150);
//prob_d_ind2_150 = load_prob_densities("ind2",f_prob_densities_ind2_150);

//prob_d_coll = load_prob_densities("coll",f_prob_densities_coll);
//prob_d_ind1 = load_prob_densities("ind1",f_prob_densities_ind1);
//prob_d_ind2 = load_prob_densities("ind2",f_prob_densities_ind2);


Spectrum s1("", kBinz, loader, DataLoader ,kCRTPMTNeutrino );

loader.Go();

double factor = s1.POT();  

TH1D* h1 = s1.ToTH1(factor);

f->cd();
tree->Write(0,TObject::kOverwrite);
f->Close();

//cout << "tot 1muNp interactions: " << tot_1muNp << endl;

}