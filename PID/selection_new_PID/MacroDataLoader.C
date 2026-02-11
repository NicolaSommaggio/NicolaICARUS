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

//MC 2D
const std::string fdata = "/storage/gpfs_data/icarus/local/users/cfarnese/production_mc_2025A_ICARUS_Overlays_BNB_MC_RUN2_summer_2025_v10_06_00_04p04/MC_overlay_neutrino_stage1_flat_cafs_v10_06_00_04p04_concat.root";

SpectrumLoader loader(fdata);       //CAF that I produced with all dedx vs rr   

const Binning kBinz = Binning::Simple(300,0,30);

Spectrum s1("", kBinz, loader, DataLoader ,kCRTPMTNeutrino );

loader.Go();

double factor = s1.POT();  

TH1D* h1 = s1.ToTH1(factor);

std::string filename;


filename = "data_struct_mc2d.root";

TFile *f = new TFile(filename.c_str(), "RECREATE");

TTree * tree = new TTree("tree","");

RecoSlice thisslice;

tree->Branch("slice",&thisslice);

for(const auto &temp_slice : slices)
{
   thisslice = temp_slice;
   tree->Fill();
}

tree->Write(0,TObject::kOverwrite);
f->Close();

}