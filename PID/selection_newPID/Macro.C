#pragma once 

#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"


//#include "helper_stitch_simulz0.h"
//#include "helper_1muNp_puro.h"
#include "helper_macro.h"
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

void Macro(){
//MC 2D
//const std::string fdata = "production_mc_2025A_ICARUS_Overlays_BNB_MC_RUN2_summer_2025_v10_06_00_06p01_flatcaf";

//mc low stat
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlay_neutrino_stage1_flat_cafs_v10_06_00_04p04_concat.root";

//10 % statistics 1.47154e+19 POT
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026/group_6657.flat.caf.root";

//100% statistics concat
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026_*/*.root";
//divided in 3
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026_3/*flat.caf.root";

//30 percent stat
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026_1/*flat.caf.root";

//1
//5.01111e+19 POT over 0 readouts
//tot 1muNp interactions: 6655

//2
//5.1204e+19 POT over 0 readouts
//tot 1muNp interactions: 6727

//3
//5.2105e+19 POT over 0 readouts
//tot 1muNp interactions: 6742


// 11 files concat
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026_3/group_6664.flat.caf.root";  // 13G 
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026_2/group_6662.flat.caf.root";  // 12G
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026_1/group_6656_2.flat.caf.root";  // 11G
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026_1/group_6657.flat.caf.root";  // 11G
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026_2/group_6661.flat.caf.root";  // 11G
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026_3/group_6663.flat.caf.root";  //9.2G
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026_3/group_6660.flat.caf.root";  //9.1G
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026_2/group_6658.flat.caf.root";  //8.2G
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026_1/group_6656_1.flat.caf.root";  //7.8G
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026_2/group_667X.flat.caf.root";  //7.7G
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026_1/group_6656_3.flat.caf.root";  //7.5G
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026_3/group_668X.flat.caf.root";  //4.2G
//const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026_3/group_6659.flat.caf.root";  //3.6G

//const std::string fdata ="/pnfs/sbn/data/sbn_fd/poms_production/mc/2025A_ICARUS_Overlays_BNB_MC_RUN2/summer_2025/v10_06_00_04p04/flatcaf/*/*/overlay_*flat.caf*.root";


SpectrumLoader loader(fdata);       //CAF that I produced with all dedx vs rr   

const Binning kBinz = Binning::Simple(300,0,30);

Spectrum s1("", kBinz, loader, DataLoader ,kCRTPMTNeutrino );

loader.Go();

double factor = s1.POT();  

TH1D* h1 = s1.ToTH1(factor);

cout << "tot 1muNp interactions: " << tot_1muNp << endl;

}