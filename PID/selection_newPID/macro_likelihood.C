#pragma once 

#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"

//#include "likelihood.h"
#include "check_daughter.h"
//#include "likelihood_old.h"
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

void macro_likelihood(){

  
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

//MC 2D small stat
//const std::string fdata = "/storage/gpfs_data/icarus/local/users/cfarnese/production_mc_2025A_ICARUS_Overlays_BNB_MC_RUN2_summer_2025_v10_06_00_04p04/MC_overlay_neutrino_stage1_flat_cafs_v10_06_00_04p04_concat.root";

//MC 2D large stat
//const std::string fdata = "production_mc_2025A_ICARUS_Overlays_BNB_MC_RUN2_summer_2025_v10_06_00_06p01_flatcaf";
const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026/*flat.caf.root";
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026_1/*flat.caf.root";
// 1 file
//const std::string fdata = "/pnfs/sbn/data/sbn_fd/poms_production/mc/2025A_ICARUS_Overlays_BNB_MC_RUN2/summer_2025/v10_06_00_06p01/flatcaf/39/51/overlay_neutrino_stage1_66812186_2689.flat.caf-009e8064-a5da-49d8-8634-3acebbf4b82f.root";
// small stat MC
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlay_neutrino_stage1_flat_cafs_v10_06_00_04p04_concat.root";

//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026_[1,2,3]/*flat.caf.root";

//MPVMPR
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/cafs_MPV_MPR_Nicola/*.root";
//1 file
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/cafs_MPV_MPR_Nicola/mc-nicola-v10-01_51_20260309T110245-G4_20260309T122412-DetSim_20260309T153647-MCstage0_20260309T153815-MCstage1.flat.caf.root";
//large stat
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/cafs_MPV_MPR_Nicola/4/*.root";
//large stat low momentum
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/cafs_MPV_MPR_Nicola/5/*.root";
//from nu
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/cafs_MPV_MPR_Nicola/nu_1/*.root";

//const std::string fdata = "/exp/icarus/data/users/nsommagg/MPVMPR/nu_1_tot/*.root";

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


//TFile * f_prob_densities_coll_150 = TFile::Open("/exp/icarus/data/users/nsommagg/HISTO_prob_densities_6_classes_30percent_mediumbinning.root", "READ");
//TFile * f_prob_densities_ind1_150 = TFile::Open("/exp/icarus/data/users/nsommagg/HISTO_prob_densities_6_classes_30percent_mediumbinning_IND1.root", "READ");
//TFile * f_prob_densities_ind2_150 = TFile::Open("/exp/icarus/data/users/nsommagg/HISTO_prob_densities_6_classes_30percent_mediumbinning_IND2.root", "READ");

TFile * f_prob_densities_coll = TFile::Open("/exp/icarus/data/users/nsommagg/HISTO_prob_densities_6_classes_30percent_stat.root", "READ");
TFile * f_prob_densities_ind1 = TFile::Open("/exp/icarus/data/users/nsommagg/HISTO_prob_densities_6_classes_30percent_IND1.root", "READ");
TFile * f_prob_densities_ind2 = TFile::Open("/exp/icarus/data/users/nsommagg/HISTO_prob_densities_6_classes_30percent_IND2.root", "READ");

//prob_d_coll_150 = load_prob_densities("coll",f_prob_densities_coll_150);
//prob_d_ind1_150 = load_prob_densities("ind1",f_prob_densities_ind1_150);
//prob_d_ind2_150 = load_prob_densities("ind2",f_prob_densities_ind2_150);

prob_d_coll = load_prob_densities("coll",f_prob_densities_coll);
prob_d_ind1 = load_prob_densities("ind1",f_prob_densities_ind1);
prob_d_ind2 = load_prob_densities("ind2",f_prob_densities_ind2);

SpectrumLoader loader(fdata);         

const Binning kBinz = Binning::Simple(300,0,30);

Spectrum s1("", kBinz, loader, likelihood_dump ,kCRTPMTNeutrino );

loader.Go();

double factor = s1.POT();  

TH1D* h1 = s1.ToTH1(factor);

/*
TFile *file_genp_MPVMPR = new TFile("file_genp_MPVMPR.root","RECREATE");
TH1D *h_genp_mu = new TH1D("genp_module_MU","genp module MU [MeV]",200, 0, 1000);
TH1D *h_genp_pro = new TH1D("genp_module_PRO","genp module PRO [MeV]",200, 0, 1000);
TH1D *h_genp_pi = new TH1D("genp_module_PI","genp module PI [MeV]",200, 0, 1000);

for(int track=0; track<(int)genp_mu.size(); track++)
{
    h_genp_mu->Fill(genp_mu[track]*1000);
}
for(int track=0; track<(int)genp_pro.size(); track++)
{
    h_genp_pro->Fill(genp_pro[track]*1000);
}
for(int track=0; track<(int)genp_pi.size(); track++)
{
    h_genp_pi->Fill(genp_pi[track]*1000);
}

h_genp_mu->Write();
h_genp_pro->Write();
h_genp_pi->Write();
*/
cout << n_mu << " " << n_pro << " " << n_pi << endl << endl;

cout << n_c0 << " muon rising" << endl;
cout << n_c1 << " muon mip" << endl;
cout << n_c2 << " proton rising" << endl;
cout << n_c3 << " proton interacting" << endl;
cout << n_c4 << " pion rising" << endl;
cout << n_c5 << " pion interacting" << endl;

}