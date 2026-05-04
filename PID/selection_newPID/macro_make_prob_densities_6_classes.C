#pragma once 

#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"

#include "make_prob_densities_6_classes.h"
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

void make_file()
{
    TFile * f = new TFile("/exp/icarus/data/users/nsommagg/HISTO_prob_densities_6_classes_30percent_IND1.root","RECREATE");
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

void macro_make_prob_densities_6_classes(){

  
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
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026/*flat.caf.root";
const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026_1/*flat.caf.root";
// 1 file
//const std::string fdata = "/pnfs/sbn/data/sbn_fd/poms_production/mc/2025A_ICARUS_Overlays_BNB_MC_RUN2/summer_2025/v10_06_00_06p01/flatcaf/39/51/overlay_neutrino_stage1_66812186_2689.flat.caf-009e8064-a5da-49d8-8634-3acebbf4b82f.root";
//mc small stat
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlay_neutrino_stage1_flat_cafs_v10_06_00_04p04_concat.root";

//30 percent stat
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026_1/*flat.caf.root";

// 11 files concat
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026/group_6656.flat.caf.root";
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026/group_6657.flat.caf.root";
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026/group_6658.flat.caf.root";
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026/group_6659.flat.caf.root";
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026/group_6660.flat.caf.root";
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026/group_6661.flat.caf.root";
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026/group_6662.flat.caf.root";
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026/group_6663.flat.caf.root";
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026/group_6664.flat.caf.root";
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026/group_667X.flat.caf.root";
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026/group_668X.flat.caf.root";


//MPVMPR from NU
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

SpectrumLoader loader(fdata);         

const Binning kBinz = Binning::Simple(300,0,30);

Spectrum s1("", kBinz, loader, select_true_class ,kCRTPMTNeutrino );

loader.Go();

double factor = s1.POT();  

TH1D* h1 = s1.ToTH1(factor);

make_file();

std::vector<std::pair<std::string,int>> classes = {{"muon_class0",0}, {"muon_class1",1}, {"proton_class2",2}, {"proton_class3",3}, {"pion_class4",4}, {"pion_class5",5}};
std::vector<std::vector<std::vector<double>>> rr_classes ={rr_mu_class0, rr_mu_class1, rr_pro_class2, rr_pro_class3, rr_pi_class4, rr_pi_class5};
std::vector<std::vector<std::vector<double>>> dedx_classes ={dedx_mu_class0, dedx_mu_class1, dedx_pro_class2, dedx_pro_class3, dedx_pi_class4, dedx_pi_class5};

for(int i=0; i<6; i++){cout << rr_classes[i].size() << endl;}

TFile *f = TFile::Open("/exp/icarus/data/users/nsommagg/HISTO_prob_densities_6_classes_30percent_IND1.root", "UPDATE");

for(const auto &clas : classes)
{
    cout << "processing " << clas.first << endl;

    TH2D *dedx_vs_range = new TH2D("dedx_vs_range","",300,0,30,300,0,30);
    
    for(int track=0; track<(int)rr_classes[clas.second].size(); track++)
    {

        for(int hit=0; hit<(int)rr_classes[clas.second][track].size(); hit++)
        {dedx_vs_range->Fill(rr_classes[clas.second][track][hit],dedx_classes[clas.second][track][hit]);}
                
    }

    cout << rr_classes[clas.second].size() << " " << clas.first << " tracks" << endl;
    TDirectory *pardir = (TDirectory*)f->Get(clas.first.c_str());
        
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

        //for(int i=0; i<(int)binlowedges_coll.size(); i++)
        //{
            //cout << binlowedges_coll[i] << " ";
        //}
        //cout << endl;

        while(keep_coll)
        {
            keep_coll = false;
            dEdx_coll = new TH1D(Form("dEdx_ind1_rr_%.1f", rr),"", Nbins_coll, binlowedges_coll.data());

            for(int track=0; track<(int)rr_classes[clas.second].size(); track++)
            {
                //FILLING THE dE/dx HISTOGRAM
                for(int hit=0; hit<(int)rr_classes[clas.second][track].size(); hit++)
                {
                    if(rr_classes[clas.second][track][hit]<rr && rr_classes[clas.second][track][hit]>=rr-0.5)
                    {
                        dEdx_coll->Fill(dedx_classes[clas.second][track][hit]);
                    }   
                }   
            }   

            for(int bin=Nbins_coll-1; bin>0; bin--)
            {
                //cout << bin+1 << " " << dEdx_coll->GetBinContent(bin+1) << " ";
                if(dEdx_coll->GetBinContent(bin+1) < 1) // bin+1 perché TH1D è 1-indexed
                {
                    keep_coll = true;
                    binlowedges_coll.erase(binlowedges_coll.begin()+bin);
                    Nbins_coll--;
                }
                //for(int i=0; i<(int)binlowedges_coll.size(); i++)
                //{
                    //cout << binlowedges_coll[i] << " ";
                //}
                //cout << endl;
            }
            //cout << "contenuto primo bin " << dEdx_coll->GetBinContent(1) << endl;
            if(dEdx_coll->GetBinContent(1)<1){binlowedges_coll.erase(binlowedges_coll.begin()+1);keep_coll = true;}
            //for(int i=0; i<(int)binlowedges_coll.size(); i++)
            //{
                //cout << binlowedges_coll[i] << " ";
            //}
            //cout << endl;
            if(keep_coll){delete dEdx_coll;}
        }
        
        dEdx_coll->Scale(1./dEdx_coll->Integral());
        dEdx_coll->Write(0,TObject::kOverwrite);
    }
    pardir->cd();
    dedx_vs_range->Write();

}//loop over all true classes

}